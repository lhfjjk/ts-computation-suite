#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""mlip_preopt.py -- DFT 之前先用 OrbMol 预优化几何（去碰撞 / 放松构象）

为什么需要它
  chain 按几何摆位构造出来的构象经常带严重原子重叠（实测 H···C 0.83 A）。
  ORCA 直接跑会：① LEANSCF 段错误（A1-face2 就是 0xc0000005）或
  ② 在极平 PES 上极限环振荡、撞 MaxIter（A2-face2 已空转 10 h）。
  先用 MLIP 在保持反应坐标约束的前提下把几何放松，再交给 DFT。

两种用法
  1) 独立 CLI（必须用 MLIP venv 的 python，因为要 torch/orb）：
       C:\\Tools\\mlip\\env\\Scripts\\python.exe mlip_preopt.py \\
           --in A1-face2.inp --out A1-face2-pre.xyz \\
           --charge 0 --spin 1 --constraint 58,60,2.15 --freeze-constrained --json rep.json
  2) 链内挂载（系统 python 也能 import；torch/ase 只在子进程里用）：
       import mlip_preopt
       mlip_preopt.enable(M, mode="clash", threshold=0.95)
     enable() 会包装 M.run：每个作业送 ORCA 前先检查最短接触，
     低于阈值才调 MLIP 预优化；失败/校验不过就原样退回，不影响链。

约定
  * 约束编号与 ORCA %geom Constraints 一致：0 基。
  * 预优化默认冻结所有约束原子（Ni 与两个配位端），只让其余部分松弛，
    这样新成键的朝向不会被换掉。
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CLI = os.path.join(HERE, "mlip_preopt.py")
MLIP_PY = os.environ.get("MLIP_PYTHON", sys.executable)
MODEL = "orb_v3_conservative_omol"


# --------------------------------------------------------------------------- #
# 轻量结构工具（不依赖 ase，链进程也能用）
# --------------------------------------------------------------------------- #
def parse_orca_constraints(text: str) -> list[tuple[int, int, float]]:
    """从 ORCA %geom Constraints 文本里抽 { B i j d C }。编号原样（0 基）。"""
    return [(int(a), int(b), float(d))
            for a, b, d in re.findall(r"\{\s*B\s+(\d+)\s+(\d+)\s+([0-9.]+)\s+C\s*\}", text or "")]


def min_contact(atoms: list) -> tuple[float, tuple[int, int]]:
    """全局最短原子间距离。正常键都 >=0.96 A（H-O），所以 <0.95 A 必为重叠。

    不用成键判据：严重重叠(0.83 A)比"正常键"还短，按共价半径分类会把它误判成键。
    """
    X = np.array([a[1:] for a in atoms], dtype=float)
    D = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2)
    np.fill_diagonal(D, 1e9)
    i, j = np.unravel_index(np.argmin(D), D.shape)
    return float(D[i, j]), (int(i), int(j))


def _write_xyz(path: str, atoms: list, comment: str = "") -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write("%d\n%s\n" % (len(atoms), comment))
        for a in atoms:
            f.write("%-3s %16.9f %16.9f %16.9f\n" % (a[0], a[1], a[2], a[3]))


def _read_xyz(path: str) -> list:
    with open(path, encoding="utf-8") as f:
        lines = f.read().splitlines()
    n = int(lines[0].strip())
    out = []
    for line in lines[2:2 + n]:
        p = line.split()
        out.append((p[0], float(p[1]), float(p[2]), float(p[3])))
    return out


# --------------------------------------------------------------------------- #
# 链内挂载
# --------------------------------------------------------------------------- #
def enable(M, mode: str = "clash", threshold: float = 0.95,
           always: bool = False, max_rmsd: float = 2.0,
           timeout: float = 900, verbose: bool = True):
    """包装 M.run(job, atoms, mult, cons, uks)：送 DFT 前按需做 MLIP 预优化。

    mode='clash'  : 只有最短接触 < threshold 才预优化（默认，最安全）
    mode='always' : 每个作业都预优化
    """
    if getattr(M, "_mlip_preopt_enabled", False):
        return
    orig_run = M.run
    work = getattr(M, "WORK", os.getcwd())

    def _log(msg):
        try:
            with open(os.path.join(work, "mlip-preopt.log"), "a", encoding="utf-8") as f:
                f.write(time.strftime("%Y-%m-%dT%H:%M:%S") + "|" + msg + "\n")
        except Exception:
            pass
        if verbose:
            print("[mlip-preopt] " + msg, flush=True)

    def run(job, atoms, mult, cons="", uks=False):
        new_atoms = atoms
        try:
            d, pair = min_contact(atoms)
            run_it = always or (mode == "clash" and d < threshold)
            if not run_it:
                _log("%s|SKIP|min-contact=%.2f" % (job, d))
            else:
                new_atoms, rep = _preopt_via_subprocess(
                    atoms, mult, cons, work, job,
                    max_rmsd=max_rmsd, timeout=timeout)
                if rep.get("ok"):
                    _log("%s|OK|contact %.2f->%.2f|rmsd=%.2f|steps=%s|%.1fs"
                         % (job, rep["contact_before"], rep["contact_after"],
                            rep["rmsd"], rep.get("steps"), rep.get("seconds", 0)))
                else:
                    new_atoms = atoms
                    _log("%s|REJECT|%s|contact=%.2f" % (job, rep.get("reason", "?"), d))
        except Exception as exc:
            new_atoms = atoms
            _log("%s|ERROR|%s" % (job, exc))
        return orig_run(job, new_atoms, mult, cons, uks)

    M.run = run
    M._mlip_preopt_enabled = True
    _log("enabled mode=%s threshold=%.2f always=%s" % (mode, threshold, always))


def _preopt_via_subprocess(atoms, mult, cons, work, tag,
                           max_rmsd=2.0, timeout=900):
    """写临时 xyz -> 调 MLIP venv 的 CLI -> 读回。任何失败都抛异常由上层退回。"""
    if not os.path.exists(MLIP_PY):
        raise RuntimeError("MLIP venv 不存在: " + MLIP_PY)
    cons_list = parse_orca_constraints(cons)
    tmpdir = tempfile.mkdtemp(prefix="mlip_preopt_", dir=work)
    tin = os.path.join(tmpdir, tag + ".xyz")
    tout = os.path.join(tmpdir, tag + "-pre.xyz")
    tjson = os.path.join(tmpdir, tag + ".json")
    _write_xyz(tin, atoms, "constructed by chain")
    cmd = [MLIP_PY, CLI, "--in", tin, "--out", tout, "--json", tjson,
           "--charge", "0", "--spin", str(mult), "--model", MODEL,
           "--max-rmsd", str(max_rmsd), "--freeze-shell", "1"]
    for (i, j, d) in cons_list:
        cmd += ["--constraint", "%d,%d,%s" % (i, j, d)]
    proc = None
    for attempt in range(2):          # 偶发 0xC0000005（torch 加载竞态），重试一次
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        if os.path.exists(tjson):
            break
        time.sleep(2)
    if not os.path.exists(tjson):
        raise RuntimeError("MLIP preopt 失败(退出码 %s): %s"
                           % (proc.returncode, (proc.stderr or proc.stdout or "")[-300:]))
    with open(tjson, encoding="utf-8") as f:
        rep = json.load(f)
    rep["constraints"] = cons_list
    if not rep.get("ok"):
        return atoms, rep
    return _read_xyz(tout), rep


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _read_structure(path):
    try:
        from ase.io import read
        return read(path)
    except Exception:
        return _read_orca_inp(path)


def _read_orca_inp(path):
    from ase import Atoms
    lines = open(path, encoding="utf-8", errors="ignore").read().splitlines()
    k = next(i for i, l in enumerate(lines) if l.strip().startswith("* xyz"))
    parts = lines[k].split()
    charge, mult = int(parts[2]), int(parts[3])
    syms, xyz = [], []
    for line in lines[k + 1:]:
        if not line.strip():
            break
        p = line.split()
        if len(p) >= 4 and p[0][0].isalpha():
            syms.append(p[0])
            xyz.append([float(p[1]), float(p[2]), float(p[3])])
    at = Atoms(symbols=syms, positions=xyz)
    at.info["charge"], at.info["spin"] = charge, mult
    return at


def _coord_sets(atoms, cutoff=2.6, exclude=frozenset()):
    """金属的配位原子索引集合（按索引比较，避免距离舍入造成的假变化）。"""
    from mlip_tools import METALS
    sym = atoms.get_chemical_symbols()
    out = {}
    for i, s in enumerate(sym):
        if s not in METALS:
            continue
        # 只比对重原子；2.6 A 处的 Ni...H 接触不是配位，不该触发"配位层变了"
        out["%s%d" % (s, i + 1)] = frozenset(
            j for j in range(len(atoms))
            if j != i and atoms.numbers[j] > 1 and j not in exclude
            and atoms.get_distance(i, j) < cutoff)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="MLIP 预优化（OrbMol）")
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--json")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--spin", type=int, default=1)
    ap.add_argument("--model", default=MODEL)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--constraint", action="append", default=[],
                    help="i,j,d（0 基，与 ORCA %%geom 一致），可重复")
    ap.add_argument("--freeze", default="", help="额外冻结原子，逗号分隔（0 基）")
    ap.add_argument("--freeze-shell", type=int, default=0,
                    help="冻结金属第 N 层配位壳（1 = 直接配位层）")
    ap.add_argument("--freeze-constrained", action="store_true",
                    help="冻结所有约束原子（保住新成键的朝向）")
    ap.add_argument("--fmax", type=float, default=0.15)
    ap.add_argument("--steps", type=int, default=500)
    ap.add_argument("--max-rmsd", type=float, default=2.0)
    ap.add_argument("--min-contact", type=float, default=0.95)
    args = ap.parse_args(argv)

    import time as _t
    from mlip_tools import (load_calculator, set_charge_spin, _attach,
                            make_constraints, run_opt, kabsch_rmsd)

    atoms = _read_structure(args.inp)
    sym = atoms.get_chemical_symbols()
    X0 = atoms.positions.copy()
    cons = [(int(a), int(b), float(d)) for a, b, d in
            (c.split(",") for c in args.constraint)]
    freeze = [int(x) for x in args.freeze.split(",") if x.strip()]
    if args.freeze_shell:
        from mlip_tools import metal_shell
        from mlip_tools import METALS
        # 约束原子必须能沿约束方向移动，不能冻；金属本身可以冻
        cons_atoms = {i for (i, _j, _d) in cons} | {j for (_i, j, _d) in cons}
        # ASE 的多约束是按顺序投影的，FixInternals 在后会把已冻结的原子挪走。
        # 所以约束原子（含金属）一律不冻，避免两套约束互相打架。
        sym_all = atoms.get_chemical_symbols()
        shell = [k for k in metal_shell(atoms, shells=args.freeze_shell)
                 if k not in cons_atoms or sym_all[k] in METALS]
        freeze += shell
    if args.freeze_constrained:
        for (i, j, _d) in cons:
            freeze += [i, j]
    freeze = sorted(set(freeze))

    def _contact(at):
        lst = [(at.get_chemical_symbols()[k], *at.positions[k]) for k in range(len(at))]
        return min_contact(lst)

    contact_before, (i0, j0) = _contact(atoms)

    calc = load_calculator("orb", args.model, args.device)
    set_charge_spin(atoms, args.charge, args.spin, "orb")
    cobj = make_constraints(atoms, 0, cons, freeze)
    if cobj:
        atoms.set_constraint(cobj)
    _attach(atoms, calc)
    t0 = _t.perf_counter()
    conv, secs, nsteps = run_opt(atoms, args.fmax, args.steps)
    X1 = atoms.positions.copy()
    contact_after, (i1, j1) = _contact(atoms)
    rmsd = kabsch_rmsd(X0, X1)

    reason = ""
    if contact_after < args.min_contact:
        reason = "仍有重叠 %.2f A" % contact_after
    elif rmsd > args.max_rmsd:
        reason = "RMSD %.2f A 超过上限" % rmsd
    else:
        # 只检查"原来的配体有没有掉"：自由基靠近让配位数变大是正常的，
        # 但原有的金属-配体键不许被拉长超过 0.4 A。
        from ase import Atoms as _A
        a0 = _A(symbols=sym, positions=X0)
        a1 = _A(symbols=sym, positions=X1)
        from mlip_tools import METALS
        bad = []
        for m in range(len(sym)):
            if sym[m] not in METALS:
                continue
            for j in range(len(sym)):
                if j == m or a0.numbers[j] <= 1:
                    continue
                d0 = a0.get_distance(m, j)
                if d0 < 2.6 and a1.get_distance(m, j) > d0 + 0.4:
                    bad.append("%s%d(%s%d %.2f->%.2f)"
                               % (sym[m], m + 1, sym[j], j + 1, d0,
                                  a1.get_distance(m, j)))
        if bad:
            reason = "配体掉了: " + ", ".join(bad[:3])
        # 约束键必须真的到位
        for (i, j, d) in cons:
            if abs(a1.get_distance(i, j) - d) > 0.15:
                reason = "约束未满足 %s%d-%s%d %.2f vs %.2f" % (
                    sym[i], i + 1, sym[j], j + 1, a1.get_distance(i, j), d)

    report = dict(ok=(reason == ""), reason=reason, n_atoms=len(atoms),
                  charge=args.charge, spin=args.spin,
                  contact_before=contact_before,
                  contact_pair=[sym[i0] + str(i0 + 1), sym[j0] + str(j0 + 1)],
                  contact_after=contact_after,
                  contact_pair_after=[sym[i1] + str(i1 + 1), sym[j1] + str(j1 + 1)],
                  rmsd=rmsd, converged=bool(conv), steps=int(nsteps),
                  seconds=round(secs, 1), freeze=freeze,
                  constraints=[list(c) for c in cons],
                  energy_eV=float(atoms.get_potential_energy()))
    if report["ok"]:
        _write_xyz(args.out, [(sym[k], *atoms.positions[k]) for k in range(len(atoms))],
                   "mlip preopt: contact %.2f -> %.2f A, rmsd %.2f A"
                   % (contact_before, contact_after, rmsd))
    if args.json:
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(report, f, indent=2, ensure_ascii=False)
    print(json.dumps({k: report[k] for k in
                      ("ok", "reason", "contact_before", "contact_after",
                       "rmsd", "converged", "steps", "seconds")},
                     ensure_ascii=False))
    return 0 if report["ok"] else 2


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""mlip_tools.py -- 基础机器学习势（MLIP）本地工具

后端
  orb  : OrbMol / orb-v3-conservative-omol（OMol25 训练，支持 charge/spin 条件）[推荐]
  mace : MACE-MP-0（材料数据训练，89 元素，无 charge/spin）

子命令
  sp    单点能 + 受力
  opt   几何优化（可冻结金属配位壳 / 固定键长）
  scan  沿一根键做约束扫描（相当于 ORCA 的 %geom Constraints 扫描）
  conf  构象搜索（随机转动可旋转键 -> MLIP 约束优化 -> 去重排序）

约定
  * 原子编号一律从 1 开始（与 ORCA 输入一致）；--zero-based 可切成 ASE 风格。
  * 文件读写用 ASE，支持 .xyz / .gjf / .com 等。
  * 能量内部 eV，打印时附 kcal/mol。

示例
  python mlip_tools.py sp complex.xyz --charge 0 --spin 2
  python mlip_tools.py opt complex.xyz -o complex.mlip.xyz --freeze-shell 1
  python mlip_tools.py scan complex.xyz --bond 58,74 --start 2.05 --stop 4.25 --n 8 \
         --freeze-shell 1 --fix 58,4,2.12 --fix 58,6,2.08
  python mlip_tools.py conf complex.xyz -o conf_run --n 64 --keep 8 --freeze-shell 1
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

try:  # Windows 控制台默认 GBK
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

from ase import Atoms
from ase.constraints import FixAtoms, FixInternals
from ase.io import read, write
from ase.optimize import BFGS, FIRE

EV2KCAL = 23.060548867
HARTREE2KCAL = 627.509474
METALS = {
    "Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Al", "Ga", "In", "Sn", "Pb", "Bi",
}


# --------------------------------------------------------------------------- #
# 模型
# --------------------------------------------------------------------------- #
def load_calculator(backend: str, model: str, device: str,
                    precision: str = "float32-high",
                    dispersion: bool | None = None):
    """返回 ASE calculator。charge/spin 由 set_charge_spin() 写到 atoms.info。"""
    if backend == "orb":
        from orb_models.forcefield import pretrained
        from orb_models.forcefield.calculator import ORBCalculator

        if not hasattr(pretrained, model):
            available = [n for n in dir(pretrained) if n.startswith("orb")]
            raise SystemExit(f"未知 orb 模型 {model!r}；可用: {available}")
        loader = getattr(pretrained, model)
        # Windows 上 torch.compile 通常没有可用编译器，显式关掉
        m = loader(device=device, precision=precision, compile=False)
        calc = ORBCalculator(m, device=device)
        return calc

    if backend == "mace":
        from mace.calculators import mace_mp

        kwargs = dict(model=model, device=device, default_dtype="float64")
        if dispersion:
            kwargs["dispersion"] = True
        try:
            return mace_mp(**kwargs)
        except Exception as exc:  # torch-dftd 未装等
            if dispersion:
                print(f"[warn] 带色散加载失败({exc})，退回无色散", file=sys.stderr)
                kwargs.pop("dispersion", None)
                return mace_mp(**kwargs)
            raise

    raise SystemExit(f"未知 backend {backend!r}")


def set_charge_spin(atoms: Atoms, charge: int, spin: int, backend: str) -> None:
    """OrbMol 需要 atoms.info['charge'] / ['spin']；MACE 忽略。"""
    if backend == "orb":
        atoms.info["charge"] = int(charge)
        atoms.info["spin"] = int(spin)


# --------------------------------------------------------------------------- #
# 几何工具
# --------------------------------------------------------------------------- #
def _bond_matrix(atoms: Atoms, mult: float = 1.15) -> np.ndarray:
    from ase.neighborlist import NeighborList, natural_cutoffs

    cutoffs = natural_cutoffs(atoms, mult=mult)
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True, skin=0.0)
    nl.update(atoms)
    return nl.get_connectivity_matrix(sparse=False)


def metal_shell(atoms: Atoms, shells: int = 1, mult: float = 1.25) -> list[int]:
    """金属原子 + 前 shells 层配位壳的原子索引（0 基）。"""
    m = _bond_matrix(atoms, mult=mult)
    n = len(atoms)
    metal = [i for i, s in enumerate(atoms.get_chemical_symbols()) if s in METALS]
    if not metal:
        return []
    selected = set(metal)
    frontier = set(metal)
    for _ in range(max(0, shells)):
        nxt: set[int] = set()
        for i in frontier:
            nxt.update(int(j) for j in np.nonzero(m[i])[0])
        nxt -= selected
        selected |= nxt
        frontier = nxt
    return sorted(selected)


def heavy_indices(atoms: Atoms) -> list[int]:
    return [i for i, z in enumerate(atoms.numbers) if z > 1]


def kabsch_rmsd(pos_a: np.ndarray, pos_b: np.ndarray, align: list[int] | None = None) -> float:
    """对齐（默认全体）后在指定原子集上算 RMSD。"""
    a = pos_a.copy()
    b = pos_b.copy()
    idx = np.array(align if align is not None else range(len(a)))
    ca, cb = a[idx].mean(0), b[idx].mean(0)
    h = (a[idx] - ca).T @ (b[idx] - cb)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    a = (a - ca) @ rot.T + cb
    return float(np.sqrt(((a - b) ** 2).sum(1).mean()))


def rotate_fragment(atoms: Atoms, pivot: int, other: int, angle_deg: float,
                    moving: list[int]) -> None:
    """绕 pivot->other 轴旋转 moving 原子（就地）。"""
    axis = atoms.positions[other] - atoms.positions[pivot]
    norm = np.linalg.norm(axis)
    if norm < 1e-8:
        return
    axis /= norm
    th = math.radians(angle_deg)
    k = axis
    rot = np.array([
        [math.cos(th) + k[0] ** 2 * (1 - math.cos(th)),
         k[0] * k[1] * (1 - math.cos(th)) - k[2] * math.sin(th),
         k[0] * k[2] * (1 - math.cos(th)) + k[1] * math.sin(th)],
        [k[1] * k[0] * (1 - math.cos(th)) + k[2] * math.sin(th),
         math.cos(th) + k[1] ** 2 * (1 - math.cos(th)),
         k[1] * k[2] * (1 - math.cos(th)) - k[0] * math.sin(th)],
        [k[2] * k[0] * (1 - math.cos(th)) - k[1] * math.sin(th),
         k[2] * k[1] * (1 - math.cos(th)) + k[0] * math.sin(th),
         math.cos(th) + k[2] ** 2 * (1 - math.cos(th))],
    ])
    origin = atoms.positions[pivot].copy()
    sel = np.array(moving, dtype=int)
    atoms.positions[sel] = (atoms.positions[sel] - origin) @ rot.T + origin


def _component(adj: np.ndarray, start: int, banned: tuple[int, int]) -> set[int]:
    seen = {start}
    stack = [start]
    a, b = banned
    while stack:
        i = stack.pop()
        for j in np.nonzero(adj[i])[0]:
            j = int(j)
            if (i == a and j == b) or (i == b and j == a):
                continue
            if j not in seen:
                seen.add(j)
                stack.append(j)
    return seen


def rotatable_bonds(atoms: Atoms, mult: float = 1.15) -> list[tuple[int, int]]:
    """可旋转单键：两端非氢非金属、各有 >=2 个重原子邻居、不在环内、非酰胺 C-N。"""
    adj = _bond_matrix(atoms, mult=mult)
    sym = atoms.get_chemical_symbols()
    out = []
    n = len(atoms)
    for i in range(n):
        for j in range(i + 1, n):
            if not adj[i, j]:
                continue
            if atoms.numbers[i] <= 1 or atoms.numbers[j] <= 1:
                continue
            if sym[i] in METALS or sym[j] in METALS:
                continue
            # 重原子邻居数（含彼此）
            deg_i = sum(1 for k in np.nonzero(adj[i])[0] if atoms.numbers[int(k)] > 1)
            deg_j = sum(1 for k in np.nonzero(adj[j])[0] if atoms.numbers[int(k)] > 1)
            if deg_i < 2 or deg_j < 2:
                continue
            # 环内（去掉这条边仍连通）
            if j in _component(adj, i, (i, j)):
                continue
            # 酰胺 C-N 近似刚性
            amide = False
            for c, nn in ((i, j), (j, i)):
                if sym[c] == "C" and sym[nn] == "N":
                    o_double = any(
                        sym[int(k)] == "O" and atoms.get_distance(c, int(k)) < 1.30
                        for k in np.nonzero(adj[c])[0]
                    )
                    if o_double:
                        amide = True
            if amide:
                continue
            out.append((i, j))
    return out


# --------------------------------------------------------------------------- #
# 优化
# --------------------------------------------------------------------------- #
def make_constraints(atoms: Atoms, freeze_shell: int = 0,
                     fix_bonds: list[tuple[int, int, float]] | None = None,
                     freeze_idx: list[int] | None = None):
    cons = []
    idx: list[int] = []
    if freeze_shell:
        idx += metal_shell(atoms, shells=freeze_shell)
    if freeze_idx:
        idx += freeze_idx
    # ASE 按列表顺序依次投影：FixInternals 必须在 FixAtoms 之前，
    # 否则它会把已冻结的原子挪走去满足键长。
    if fix_bonds:
        cons.append(FixInternals(bonds=[[d, [i, j]] for (i, j, d) in fix_bonds]))
    if idx:
        cons.append(FixAtoms(indices=sorted(set(idx))))
    return cons


def set_bond_distance(atoms: Atoms, i: int, j: int, target: float) -> None:
    """平移较小一侧的原子，使 i-j 距离等于 target（不改变内部几何）。"""
    d = atoms.get_distance(i, j)
    if d < 1e-8:
        return
    vec = atoms.positions[j] - atoms.positions[i]
    shift = vec * (target - d) / d
    adj = _bond_matrix(atoms)
    side_j = _component(adj, j, (i, j))
    side_i = _component(adj, i, (i, j))
    if len(side_j) <= len(side_i):
        sel = np.array(sorted(side_j), dtype=int)
        atoms.positions[sel] += shift
    else:
        sel = np.array(sorted(side_i), dtype=int)
        atoms.positions[sel] -= shift


def run_opt(atoms: Atoms, fmax: float = 0.05, steps: int = 400,
            algo: str = "bfgs", logfile: str | None = None):
    opt_cls = BFGS if algo == "bfgs" else FIRE
    opt = opt_cls(atoms, logfile=logfile, trajectory=None)
    t0 = time.perf_counter()
    conv = opt.run(fmax=fmax, steps=steps)
    return bool(conv), time.perf_counter() - t0, int(opt.nsteps)


def _attach(atoms: Atoms, calc) -> None:
    # ASE 的结果缓存只比较几何，不看 atoms.info；换自旋/电荷时必须显式清空
    calc.reset()
    atoms.calc = calc


def _forces_max(atoms: Atoms) -> float:
    f = atoms.get_forces()
    return float(np.sqrt((f ** 2).sum(1)).max()) if len(f) else 0.0


# --------------------------------------------------------------------------- #
# 子命令
# --------------------------------------------------------------------------- #
def cmd_sp(args) -> int:
    calc = load_calculator(args.backend, args.model, args.device,
                           args.precision, args.dispersion)
    rows = []
    for path in args.files:
        atoms = read(path)
        set_charge_spin(atoms, args.charge, args.spin, args.backend)
        _attach(atoms, calc)
        t0 = time.perf_counter()
        e = atoms.get_potential_energy()
        fmax = _forces_max(atoms)
        dt = time.perf_counter() - t0
        rows.append(dict(file=str(path), n_atoms=len(atoms), energy_eV=e,
                         fmax_eV_A=fmax, seconds=dt))
        print(f"{Path(path).name:40s} E={e:.6f} eV  fmax={fmax:.3f} eV/A  {dt*1000:.0f} ms")
    if args.json:
        Path(args.json).write_text(json.dumps(rows, indent=2, ensure_ascii=False), encoding="utf-8")
    if len(rows) > 1:
        e0 = min(r["energy_eV"] for r in rows)
        print("\n相对能量（最低点为 0）:")
        for r in rows:
            print(f"  {Path(r['file']).name:40s} {r['energy_eV']-e0:+.3f} eV "
                  f"= {(r['energy_eV']-e0)*EV2KCAL:+.2f} kcal/mol")
    return 0


def cmd_opt(args) -> int:
    calc = load_calculator(args.backend, args.model, args.device,
                           args.precision, args.dispersion)
    atoms = read(args.file)
    set_charge_spin(atoms, args.charge, args.spin, args.backend)
    pos0 = atoms.positions.copy()
    e0 = None
    cons = make_constraints(atoms, args.freeze_shell, _parse_fix(args.fix),
                            _parse_idx(args.freeze))
    if cons:
        atoms.set_constraint(cons)
    _attach(atoms, calc)
    e0 = atoms.get_potential_energy()
    conv, dt, nsteps = run_opt(atoms, args.fmax, args.steps, args.algo,
                               args.logfile if args.logfile != "-" else "-")
    e1 = atoms.get_potential_energy()
    sym = atoms.get_chemical_symbols()
    rmsd = kabsch_rmsd(pos0, atoms.positions)
    print(f"收敛: {conv}  步数: {nsteps}  用时: {dt:.1f} s")
    print(f"E: {e0:.6f} -> {e1:.6f} eV  (dE = {(e1-e0)*EV2KCAL:+.2f} kcal/mol)")
    print(f"全原子 RMSD(优化前后): {rmsd:.3f} A   fmax: {_forces_max(atoms):.4f} eV/A")
    for m in [i for i, s in enumerate(sym) if s in METALS]:
        partners = [(j, atoms.get_distance(m, j)) for j in range(len(atoms))
                    if j != m and atoms.get_distance(m, j) < args.shell_cutoff]
        if partners:
            print(f"  {sym[m]}{m+1} 配位: " + "  ".join(f"{sym[j]}{j+1}-{d:.2f}" for j, d in partners))
    out = Path(args.out) if args.out else Path(args.file).with_suffix(".mlip.xyz")
    out.parent.mkdir(parents=True, exist_ok=True)
    atoms.info["mlip_energy_eV"] = float(e1)
    write(out, atoms)
    print(f"写出: {out}")
    return 0


def cmd_scan(args) -> int:
    calc = load_calculator(args.backend, args.model, args.device,
                           args.precision, args.dispersion)
    atoms = read(args.file)
    set_charge_spin(atoms, args.charge, args.spin, args.backend)
    i, j = _idx(args.bond[0]), _idx(args.bond[1])
    targets = np.linspace(args.start, args.stop, args.n)
    base_fix = _parse_fix(args.fix)
    cons = make_constraints(atoms, args.freeze_shell, base_fix, _parse_idx(args.freeze))
    if cons:
        atoms.set_constraint(cons)
    _attach(atoms, calc)
    rows = []
    for k, d in enumerate(targets):
        set_bond_distance(atoms, i, j, float(d))
        cons = make_constraints(atoms, args.freeze_shell,
                                base_fix + [(i, j, float(d))], _parse_idx(args.freeze))
        atoms.set_constraint(cons)
        conv, dt, nsteps = run_opt(atoms, args.fmax, args.steps, args.algo, None)
        e = atoms.get_potential_energy()
        actual = atoms.get_distance(i, j)
        imax = _forces_max(atoms)
        rows.append(dict(point=k + 1, target=float(d), actual=actual, energy_eV=e,
                         converged=conv, steps=nsteps, fmax=imax))
        print(f"  p{k+1}  d={actual:.4f} A  E={e:.6f} eV  conv={conv} steps={nsteps} "
              f"fmax={imax:.3f}  ({dt:.1f}s)")
    e0 = min(r["energy_eV"] for r in rows)
    print("\n相对能量:")
    for r in rows:
        print(f"  d={r['actual']:.3f}  {(r['energy_eV']-e0)*EV2KCAL:+7.2f} kcal/mol")
    if args.csv:
        with open(args.csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"\nCSV: {args.csv}")
    if args.out:
        outdir = Path(args.out)
        outdir.mkdir(parents=True, exist_ok=True)
        for r in rows:
            pass
        # 逐点几何：重新跑一遍太贵，直接存最后几何不可靠 -> 用轨迹
        # 扫描过程中的几何在 atoms 上已覆盖，这里存最终点
        write(outdir / "final.xyz", atoms)
    return 0


def cmd_conf(args) -> int:
    calc = load_calculator(args.backend, args.model, args.device,
                           args.precision, args.dispersion)
    base = read(args.file)
    set_charge_spin(base, args.charge, args.spin, args.backend)
    rng = np.random.default_rng(args.seed)
    bonds = rotatable_bonds(base)
    core = metal_shell(base, shells=args.freeze_shell) if args.freeze_shell else []
    if not bonds:
        raise SystemExit("没找到可旋转键")
    print(f"可旋转键 {len(bonds)} 根，核心冻结原子 {len(core)} 个，采样 {args.n} 个构象")

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    results = []
    t_start = time.perf_counter()
    for k in range(args.n):
        at = base.copy()
        set_charge_spin(at, args.charge, args.spin, args.backend)
        adj = _bond_matrix(at)
        for (i, j) in bonds:
            # 旋转原子较少的一侧
            side_j = _component(adj, j, (i, j))
            side_i = _component(adj, i, (i, j))
            if len(side_j) <= len(side_i):
                pivot, other, moving = i, j, sorted(side_j)
            else:
                pivot, other, moving = j, i, sorted(side_i)
            if core and all(x in core for x in moving):
                continue  # 核心内部不转
            rotate_fragment(at, pivot, other, rng.uniform(-180, 180), moving)
        cons = [FixAtoms(indices=core)] if core else []
        if cons:
            at.set_constraint(cons)
        _attach(at, calc)
        conv, dt, nsteps = run_opt(at, args.fmax, args.steps, args.algo, None)
        e = at.get_potential_energy()
        results.append(dict(k=k + 1, energy=e, conv=conv, steps=nsteps, atoms=at,
                            seconds=dt))
        if (k + 1) % max(1, args.n // 10) == 0:
            print(f"  ... {k+1}/{args.n}  E_min={min(r['energy'] for r in results):.4f} eV")

    results.sort(key=lambda r: r["energy"])
    kept = []
    for r in results:
        dup = False
        for q in kept:
            if kabsch_rmsd(r["atoms"].positions, q["atoms"].positions, align=core or None) < args.rmsd_cut:
                dup = True
                break
        if not dup:
            kept.append(r)
        if len(kept) >= args.keep:
            break

    e_ref = kept[0]["energy"] if kept else 0.0
    print(f"\n保留 {len(kept)} 个去重构象（RMSD 阈值 {args.rmsd_cut} A），"
          f"总用时 {time.perf_counter()-t_start:.0f} s")
    for rank, r in enumerate(kept, 1):
        name = f"conf_{rank:02d}_E{r['energy']-e_ref:+.2f}eV.xyz"
        r["atoms"].info["mlip_energy_eV"] = float(r["energy"])
        write(outdir / name, r["atoms"])
        print(f"  {name}  dE={(r['energy']-e_ref)*EV2KCAL:+.2f} kcal/mol  conv={r['conv']}")
    with open(outdir / "summary.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["rank", "energy_eV", "dE_kcal", "converged", "steps", "file"])
        for rank, r in enumerate(kept, 1):
            w.writerow([rank, r["energy"], (r["energy"] - e_ref) * EV2KCAL,
                        r["conv"], r["steps"], f"conf_{rank:02d}.xyz"])
    print(f"输出目录: {outdir}")
    return 0


# --------------------------------------------------------------------------- #
# 参数
# --------------------------------------------------------------------------- #
ZERO_BASED = False


def _idx(v: int) -> int:
    """CLI 默认 1 基（坐标块风格）；--zero-based 时按 ASE/ORCA %geom 风格。"""
    return int(v) if ZERO_BASED else int(v) - 1


def _parse_fix(items) -> list[tuple[int, int, float]]:
    out = []
    for it in items or []:
        a, b, d = it.split(",")
        out.append((_idx(int(a)), _idx(int(b)), float(d)))
    return out


def _parse_idx(items) -> list[int]:
    out = []
    for it in items or []:
        out.append(_idx(int(it)))
    return out


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="基础 MLIP 本地工具（OrbMol / MACE-MP-0）")
    p.add_argument("--backend", default="orb", choices=["orb", "mace"])
    p.add_argument("--model", default="orb_v3_conservative_omol")
    p.add_argument("--device", default="cuda")
    p.add_argument("--precision", default="float32-high")
    p.add_argument("--charge", type=int, default=0)
    p.add_argument("--spin", type=int, default=1, help="自旋多重度（OrbMol 用）")
    p.add_argument("--dispersion", action="store_true", help="MACE 加 D3 色散")
    p.add_argument("--zero-based", action="store_true",
                   help="原子编号按 0 基（ORCA %%geom Constraints 就是 0 基）")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("sp", help="单点能")
    s.add_argument("files", nargs="+")
    s.add_argument("--json")
    s.set_defaults(func=cmd_sp)

    s = sub.add_parser("opt", help="几何优化")
    s.add_argument("file")
    s.add_argument("-o", "--out")
    s.add_argument("--fmax", type=float, default=0.05)
    s.add_argument("--steps", type=int, default=400)
    s.add_argument("--algo", default="bfgs", choices=["bfgs", "fire"])
    s.add_argument("--freeze-shell", type=int, default=0, help="冻结金属配位壳层数")
    s.add_argument("--freeze", nargs="*", help="额外冻结的原子编号(1 基)")
    s.add_argument("--fix", nargs="*", help="固定键长 i,j,d(Å)，1 基")
    s.add_argument("--shell-cutoff", type=float, default=2.8)
    s.add_argument("--logfile", default=None)
    s.set_defaults(func=cmd_opt)

    s = sub.add_parser("scan", help="沿键约束扫描")
    s.add_argument("file")
    s.add_argument("--bond", nargs=2, type=int, required=True)
    s.add_argument("--start", type=float, required=True)
    s.add_argument("--stop", type=float, required=True)
    s.add_argument("--n", type=int, default=8)
    s.add_argument("--fmax", type=float, default=0.08)
    s.add_argument("--steps", type=int, default=250)
    s.add_argument("--algo", default="bfgs", choices=["bfgs", "fire"])
    s.add_argument("--freeze-shell", type=int, default=0)
    s.add_argument("--freeze", nargs="*")
    s.add_argument("--fix", nargs="*")
    s.add_argument("--csv")
    s.add_argument("-o", "--out")
    s.set_defaults(func=cmd_scan)

    s = sub.add_parser("conf", help="构象搜索")
    s.add_argument("file")
    s.add_argument("-o", "--out", default="conf_run")
    s.add_argument("--n", type=int, default=32)
    s.add_argument("--keep", type=int, default=8)
    s.add_argument("--seed", type=int, default=20261028)
    s.add_argument("--fmax", type=float, default=0.10)
    s.add_argument("--steps", type=int, default=250)
    s.add_argument("--algo", default="bfgs", choices=["bfgs", "fire"])
    s.add_argument("--freeze-shell", type=int, default=1)
    s.add_argument("--rmsd-cut", type=float, default=0.40)
    s.set_defaults(func=cmd_conf)
    return p


def main(argv=None) -> int:
    global ZERO_BASED
    args = build_parser().parse_args(argv)
    ZERO_BASED = bool(args.zero_based)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())

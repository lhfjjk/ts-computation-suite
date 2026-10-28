#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ts_mlip_screen.py -- 在 DFT 精修之前，用 MLIP 把 TS 候选筛一遍。

插在 ts-generation-workflow 的 L2（OptTS）之前：

    L0 声明 ──► [手搭 / GOAT-TS / UniTS-Gen] ──► ★本脚本★ ──► L1 参考点 ──► L2 OptTS ──► L3 IRC

它做三件事（都是几何层面的，不碰能量结论）：
  1. 合理性：单点最大受力、最短接触；
  2. 是不是"垒"：沿声明的反应坐标把候选往两侧各推一段，再让 MLIP 自由弛豫；
     两次落到**同一个**极小点 ⇒ 候选在盆地里，不可能是 TS（直接淘汰）；
  3. 是不是**目标通道**：两个端点的键图差值应当恰好是「形成 forming 对、断开 breaking 对」；
     差值里没有这两对、或出现计划外的键 ⇒ 判为 WRONG-CHANNEL。

  判决：BARRIER-LIKE / BASIN / WRONG-CHANNEL / PES-MISMATCH。
  PES-MISMATCH = 候选比两侧极小点高 60+ kcal，说明这个模型的势能面不支持这条反应
  （典型：BS 破缺单重态、强开壳层金属）。此时**不要采信 BASIN**，直接去 DFT 判。

用法
  C:/Tools/mlip/env/Scripts/python.exe ts_mlip_screen.py \
      --cand guesses.xyz --forming 21,72 --breaking 1,21 \
      --charge 0 --spin 2 --out screen/

  # 有参考端点时，额外报两次弛豫各自落在谁身上
  ... --reactant enc.xyz --product prod.xyz

输出
  screen/<candidate>.json        每次判定的完整数据
  screen/<candidate>-toward-product.xyz / -toward-reactant.xyz   两个端点几何
  控制台汇总表：verdict / MLIP 垒估计 / 键图差

⚠️ 边界（与 RESULTS.md 一致）
  * MLIP 的**能量**只作粗排，不能写进论文；本脚本的判决只看**几何**（键图/落点）。
  * BS 破缺单重态无法用 MLIP 表达（它只有离散自旋标签）：这类体系请用 spin=1，
    且把结果当"几何提示"，不要采信它的能量。
  * 这里不做频率。判"几个虚频"仍然只能靠 DFT 的 Freq。
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)   # ts_common.py / mlip_tools.py 与本脚本同目录

import ts_common as TC  # noqa: E402
from mlip_tools import (EV2KCAL, load_calculator, set_charge_spin, _attach,  # noqa: E402
                        run_opt)


def _to_ase(frame):
    from ase import Atoms
    return Atoms(symbols=[a[0] for a in frame], positions=[a[1] for a in frame])


def _to_frame(atoms):
    return [(s, tuple(p)) for s, p in zip(atoms.get_chemical_symbols(), atoms.positions)]


def min_contact(frame):
    """ts_common 的 min_contact 返回 (i, j, d, bonded) 或 None。"""
    mc = TC.min_contact(frame)
    return (float(mc[2]), (mc[0], mc[1])) if mc else (float("nan"), (-1, -1))


def push_pairs(atoms, pairs, delta):
    """把每对距离改变 delta（Å）：两端各沿连线移动一半。配对可以共享原子。"""
    for (a, b) in pairs:
        i, j = a - 1, b - 1
        v = atoms.positions[j] - atoms.positions[i]
        d = float(np.linalg.norm(v))
        if d < 1e-6:
            continue
        u = v / d
        atoms.positions[i] -= 0.5 * delta * u
        atoms.positions[j] += 0.5 * delta * u


def relax(frame, calc, spin, charge, fmax, steps):
    at = _to_ase(frame)
    set_charge_spin(at, charge, spin, "orb")
    _attach(at, calc)
    e0 = float(at.get_potential_energy())
    conv, secs, n = run_opt(at, fmax, steps)
    e1 = float(at.get_potential_energy())
    return at, e0, e1, conv, n, secs


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="TS 候选的 MLIP 预筛")
    ap.add_argument("--cand", nargs="+", required=True,
                    help="候选几何（可多帧 xyz，可多个文件）")
    ap.add_argument("--forming", nargs="+", default=[], help="成键对 i,j（1 基，可重复）")
    ap.add_argument("--breaking", nargs="+", default=[], help="断键对 i,j（1 基，可重复）")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--spin", type=int, default=1, help="自旋多重度（BS 单重态填 1）")
    ap.add_argument("--model", default="orb_v3_conservative_omol")
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--push", type=float, default=0.35, help="往每侧推的键长变化（Å）")
    ap.add_argument("--fmax", type=float, default=0.15)
    ap.add_argument("--steps", type=int, default=400)
    ap.add_argument("--bs", action="store_true",
                    help="BS 破缺单重态：只报几何，不做连通性判决（模型无法表达该电子态）")
    ap.add_argument("--reactant")
    ap.add_argument("--product")
    ap.add_argument("--out", default="ts_screen")
    ap.add_argument("--json")
    args = ap.parse_args(argv)

    pairs_f = [tuple(int(x) for x in p.split(",")) for p in args.forming]
    pairs_b = [tuple(int(x) for x in p.split(",")) for p in args.breaking]
    if not pairs_f and not pairs_b:
        raise SystemExit("至少要给 --forming 或 --breaking（1 基）")

    frames = []
    for path in args.cand:
        for k, fr in enumerate(TC.read_xyz_frames(path)):
            frames.append((Path(path).stem if len(frames) == 0 else
                           "%s_%d" % (Path(path).stem, k), fr))
    ref_r = TC.read_xyz_frames(args.reactant, last_only=True)[0] if args.reactant else None
    ref_p = TC.read_xyz_frames(args.product, last_only=True)[0] if args.product else None

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    calc = load_calculator("orb", args.model, args.device)

    rows = []
    for name, fr in frames:
        rec = {"candidate": name, "n_atoms": len(fr)}
        at0 = _to_ase(fr)
        set_charge_spin(at0, args.charge, args.spin, "orb")
        _attach(at0, calc)
        e_ts = float(at0.get_potential_energy())
        fmax_ts = float(np.sqrt((at0.get_forces() ** 2).sum(1)).max())
        cmin, cpair = min_contact(fr)
        rec.update(e_ts=e_ts, fmax_ts=fmax_ts, contact=float(cmin))

        endpoints = {}
        for tag, sign in (("product", +1), ("reactant", -1)):
            at = _to_ase(fr)
            set_charge_spin(at, args.charge, args.spin, "orb")
            # 往产物侧：成键对缩短、断键对拉长；反应物侧反过来
            push_pairs(at, pairs_f, -sign * args.push)
            push_pairs(at, pairs_b, +sign * args.push)
            at2, _, e_end, conv, n, secs = relax(_to_frame(at), calc, args.spin,
                                                 args.charge, args.fmax, args.steps)
            endpoints[tag] = (at2, e_end, conv, n, secs)
            TC.write_xyz(str(outdir / ("%s-toward-%s.xyz" % (name, tag))),
                         _to_frame(at2), "MLIP descent toward %s" % tag)

        (ae, ee, ce, ne, _), (ar, er, cr, nr, _) = endpoints["product"], endpoints["reactant"]
        same = TC.kabsch_rmsd(_to_frame(ae), _to_frame(ar)) < 0.35
        gained, lost = TC.connectivity_delta(_to_frame(ar), _to_frame(ae))
        want_g = {tuple(sorted(p)) for p in pairs_f}
        want_l = {tuple(sorted(p)) for p in pairs_b}
        got_g = {tuple(sorted(p)) for p in gained}
        got_l = {tuple(sorted(p)) for p in lost}

        barrier = float((e_ts - min(ee, er)) * EV2KCAL)
        if args.bs:
            verdict = "SKIP-BS"
        elif barrier > 60.0:
            # 候选点在 MLIP 面上比两侧极小点高 60+ kcal：不是"盆地里"，而是
            # 这个模型的势能面根本没有这条反应（典型：BS 单重态、强开壳层金属体系）。
            verdict = "PES-MISMATCH"
        elif same:
            verdict = "BASIN"
        elif want_g <= got_g and want_l <= got_l:
            verdict = "BARRIER-LIKE"
        else:
            verdict = "WRONG-CHANNEL"

        rec.update(verdict=verdict,
                   same_endpoint=bool(same),
                   barrier_mlip_kcal=barrier,
                   gained=gained, lost=lost,
                   missing_forming=sorted(want_g - got_g),
                   missing_breaking=sorted(want_l - got_l),
                   conv=[bool(ce), bool(cr)], steps=[ne, nr],
                   e_end_product=ee, e_end_reactant=er)
        if ref_r is not None:
            rec["rmsd_to_reactant"] = float(TC.kabsch_rmsd(ref_r, _to_frame(ar), heavy_only=True))
        if ref_p is not None:
            rec["rmsd_to_product"] = float(TC.kabsch_rmsd(ref_p, _to_frame(ae), heavy_only=True))
        rows.append(rec)

        with open(outdir / ("%s.json" % name), "w", encoding="utf-8") as f:
            json.dump(rec, f, indent=2, ensure_ascii=False)

    print("%-22s %-14s %8s %8s %7s  %s"
          % ("candidate", "verdict", "fmax", "contact", "dE#", "bond-graph delta"))
    for r in rows:
        print("%-22s %-14s %8.2f %8.2f %7.1f  +%s -%s"
              % (r["candidate"], r["verdict"], r["fmax_ts"], r["contact"],
                 r["barrier_mlip_kcal"], r["gained"][:3], r["lost"][:3]))
    print("\n提示：BARRIER-LIKE 才值得进 L1/L2；BASIN 直接淘汰；"
          "WRONG-CHANNEL 先回 L0-a 核对反应对/结构身份。")
    if args.json:
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(rows, f, indent=2, ensure_ascii=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

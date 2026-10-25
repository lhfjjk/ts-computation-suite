#!/usr/bin/env python3
"""过渡态寻找工作流：分级生成 + 带门禁的串行执行。

实现本技能的四级阶梯（L1→L2→L3→L5），并把每道门禁
程序化。设计目标：**即使没找到真 TS，L1 的参考点也一定产出可用数据**；
**每种失败都被路由到特定补救**，而不是"重跑一次试试"。

生成四个作业（工序写进作业名，门禁据名判读）：
  {n}-L1-constOptFreq  L1 约束优化+频率 —— 产出 .hess 与模式谱
  {n}-L2-TSopt         L2 真鞍点：OptTS+Freq，InHess Read 读 L1 的 Hessian
  {n}-L3-relax         L3 从 TS 向下弛豫（IRC 的廉价替代），看连到哪个极小点
  sp-{n}-L2-TS         L5 高一级单点，**保留 L2 的电子态设置**

执行时逐级过门禁，任一不过就停在那一级并给出诊断与出口建议。

用法：
  python ts_ladder.py --name P1-bHE --xyz start.xyz \\
      --forming 56,57 --breaking 21,57 \\
      --hold 56,57=1.400 --hold 21,57=1.750 \\
      --charge 0 --mult 2 --solvent Toluene \\
      --method TPSS --basis def2-SVP --disp D3BJ \\
      --outdir ladder_P1bHE

  # BS 破缺对称面（注意 FlipSpin 用 1-based，与约束的 0-based 不同）
  python ts_ladder.py ... --uks --flipspin 72 --finalms 0 --moinp start.gbw
"""

import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ts_common import read_xyz_frames, xyz_repr, formula  # noqa: E402


def route(method, basis, disp, extra=()):
    parts = ["!"]
    for kw in (method, disp, basis, "RIJCOSX", "TightSCF"):
        if kw and kw.lower() != "none":
            parts.append(kw)
    parts.extend(extra)
    return " ".join(parts)


def blocks(args, moinp):
    out = ""
    if args.solvent:
        out += f'\n%cpcm smd true\n   smdsolvent "{args.solvent}"\nend\n'
    out += f"\n%pal nprocs {args.nprocs} end\n%maxcore {args.maxcore}\n"
    if args.mulliken:
        out += "\n%output\n   Print[P_Mulliken] 3\nend\n"
    if args.uks or args.flipspin:
        out += f'\n%scf\n   moinp "{moinp}"\n'
        if args.flipspin:
            out += f"   FlipSpin {args.flipspin}\n"     # 1-based
        if args.finalms is not None:
            out += f"   FinalMs {args.finalms}\n"
        out += "end\n"
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", required=True, help="作业名前缀(工序会追加)")
    ap.add_argument("--xyz", required=True, help="起点几何")
    ap.add_argument("--forming", action="append", default=[], metavar="I,J",
                    help="正在形成的键(1-based),可重复")
    ap.add_argument("--breaking", action="append", default=[], metavar="I,J",
                    help="正在断裂的键(1-based),可重复")
    ap.add_argument("--hold", action="append", default=[], metavar="I,J=D",
                    help="L1 约束:把 i-j 钉在 D 埃(1-based),可重复")
    ap.add_argument("--constraint-cutoff", type=float, default=None,
                    help="L1 约束用'当前几何距离'时的不动点(默认取 --hold)")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--mult", type=int, default=1)
    ap.add_argument("--method", default="TPSS")
    ap.add_argument("--basis", default="def2-SVP")
    ap.add_argument("--disp", default="D3BJ")
    ap.add_argument("--sp-method", default="TPSS0")
    ap.add_argument("--sp-basis", default="def2-TZVP")
    ap.add_argument("--sp-disp", default="D4")
    ap.add_argument("--solvent", default="Toluene")
    ap.add_argument("--maxiter", type=int, default=300)
    ap.add_argument("--nprocs", type=int, default=12)
    ap.add_argument("--maxcore", type=int, default=1100)
    ap.add_argument("--uks", action="store_true", help="开壳层/BS 面")
    ap.add_argument("--flipspin", type=int, default=None,
                    help="FlipSpin 原子号(**1-based**)")
    ap.add_argument("--finalms", type=int, default=None)
    ap.add_argument("--moinp", default=None, help="L1 的初始 gbw")
    ap.add_argument("--mulliken", action="store_true", default=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    frames = read_xyz_frames(args.xyz)
    if not frames:
        sys.exit(f"读不到几何: {args.xyz}")
    atoms = frames[0]
    n = len(atoms)
    print(f"起点: {args.xyz}  {n} 原子  {formula(atoms)}")

    def parse_pair(spec):
        a, b = [int(x) for x in spec.replace(" ", ",").split(",")]
        return a, b

    form = [parse_pair(s) for s in args.forming]
    brk = [parse_pair(s) for s in args.breaking]
    holds = []
    for spec in args.hold:
        ij, d = spec.split("=")
        a, b = parse_pair(ij)
        holds.append((a, b, float(d)))
    if not holds:
        sys.exit("必须给 --hold i,j=D 指定 L1 的约束（TS 合理值）")

    # ---- L0-a 结构身份回显 ----
    print("\n=== L0-a 结构身份核对（逐个看元素符号）===")
    for tag, pairs in (("形成", form), ("断裂", brk)):
        for a, b in pairs:
            if not (1 <= a <= n and 1 <= b <= n):
                sys.exit(f"--{tag} {a}-{b} 越界(共 {n} 原子)")
            print(f"  {tag} {a:>3}-{b:<3} -> {atoms[a-1][0]}{a}-{atoms[b-1][0]}{b}"
                  f"  -> ORCA 约束 {{ B {a-1} {b-1} <d> C }}")
    if args.flipspin:
        if not (1 <= args.flipspin <= n):
            sys.exit(f"--flipspin {args.flipspin} 越界(共 {n} 原子)")
        print(f"  FlipSpin {args.flipspin} (**1-based**) -> "
              f"{atoms[args.flipspin-1][0]}{args.flipspin}")
    print(f"  电荷/多重度 = {args.charge} {args.mult}"
          + ("   [UKS/BS 面]" if args.uks or args.flipspin else ""))

    os.makedirs(args.outdir, exist_ok=True)
    L1 = f"{args.name}-L1-constOptFreq"
    L2 = f"{args.name}-L2-TSopt"
    L3 = f"{args.name}-L3-relax"
    SP = f"sp-{args.name}-L2-TS"
    # 每个作业在自己的子目录里跑(cwd = 该目录),所以跨作业引用一律 ../<作业>/<文件>
    jobs = {}

    # ---- L1: 约束 Opt + Freq ----
    hold_lines = "\n".join(f"    {{ B {a-1} {b-1} {d:.3f} C }}" for a, b, d in holds)
    jobs[L1] = (route(args.method, args.basis, args.disp, ("Opt", "Freq"))
                + blocks(args, args.moinp or "scratch.gbw")
                + f"\n%geom\n  MaxIter {args.maxiter}\n  Constraints\n"
                  f"{hold_lines}\n  end\nend\n"
                + f"\n* xyz {args.charge} {args.mult}\n{xyz_repr(atoms)}\n*\n")

    # ---- L2: OptTS + Freq，读 L1 的 Hessian ----
    jobs[L2] = (route(args.method, args.basis, args.disp, ("OptTS", "Freq"))
                + blocks(args, f"../{L1}/{L1}.gbw")
                + f'\n%geom\n  MaxIter {args.maxiter}\n  InHess Read\n'
                  f'  InHessName "../{L1}/{L1}.hess"\nend\n'
                + f"\n* xyzfile {args.charge} {args.mult} ../{L1}/{L1}.xyz\n")

    # ---- L3: 向下弛豫（IRC 的廉价替代）----
    jobs[L3] = (route(args.method, args.basis, args.disp, ("Opt",))
                + blocks(args, f"../{L2}/{L2}.gbw")
                + f"\n%geom\n  MaxIter {args.maxiter}\nend\n"
                + f"\n* xyzfile {args.charge} {args.mult} ../{L2}/{L2}.xyz\n")

    # ---- L5: 高一级单点（保留电子态设置！）----
    jobs[SP] = (route(args.sp_method, args.sp_basis, args.sp_disp)
                + blocks(args, f"../{L2}/{L2}.gbw")
                + f"\n* xyzfile {args.charge} {args.mult} ../{L2}/{L2}.xyz\n")

    for job, text in jobs.items():
        d = os.path.join(args.outdir, job)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, job + ".inp"), "w", encoding="utf-8") as fh:
            fh.write(text)
        # 起点几何一并放进每个目录，便于人工查看/重跑
        with open(os.path.join(d, "start.xyz"), "w", encoding="utf-8") as fh:
            fh.write(f"{n}\nstart\n" + xyz_repr(atoms) + "\n")

    print(f"\n生成 {len(jobs)} 个作业于 {args.outdir}/")
    for job in jobs:
        print(f"  {job}/")
    print("\nL1 若用 BS/UKS，需把起点的 gbw 放进 L1 目录（--moinp 指定文件名）。")
    print("串行执行并逐级过门禁：")
    print(f"  python ts_ladder_run.py --dir {args.outdir} "
          f"--reactive " + " --reactive ".join(
              f"{a},{b}" for a, b in (form + brk)))
    return 0


if __name__ == "__main__":
    sys.exit(main())

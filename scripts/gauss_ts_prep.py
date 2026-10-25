#!/usr/bin/env python3
"""Generate a Gaussian 16 two-step TS ladder + IRC + single point.

The ladder mirrors the ORCA L1/L2/L3/L5 protocol:

  01_const_optfreq.gjf   constrained Opt + Freq at the guessed TS distance
                         (produces the Hessian in the .chk)
  02_ts_optfreq.gjf      Opt=(TS,ReadFC,NoEigen) + Freq, reading 01's chk
  03_irc.gjf             IRC=(ReadFC) from the TS chk
  04_sp.gjf              higher-level single point at the TS geometry

`ReadFC` is the Gaussian analogue of ORCA's `InHess Read`: it uses the
analytic Hessian from step 01 instead of the model Hessian. If your Gaussian
version or checkpoint does not contain a usable Hessian, use --ts-route calcfc.

Examples
--------
python gauss_ts_prep.py --xyz guess.xyz --name g6TS --charge 0 --mult 2 \
    --method UB3LYP --basis "6-31G(d)" --disp "em=gd3bj" --solvent toluene \
    --constraints "1,2,1.500;3,4,2.100" --nproc 24 --mem 48GB \
    --outdir gauss_g6TS --slurm --partition community

The generated job names are self-describing; run them in order.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ts_common import read_xyz_frames, xyz_repr  # noqa: E402


def parse_constraints(spec: str):
    out = []
    if not spec:
        return out
    for item in spec.split(";"):
        item = item.strip()
        if not item:
            continue
        parts = [x for x in item.replace(" ", ",").split(",") if x]
        if len(parts) != 3:
            raise SystemExit(f"bad --constraints item: {item!r}; use i,j,d")
        out.append((int(parts[0]), int(parts[1]), float(parts[2])))
    return out


def route(method, basis, disp, solvent, extra, tail):
    bits = [f"#p {method}/{basis}", tail]
    if disp and disp.lower() != "none":
        bits.append(disp)
    if solvent and solvent.lower() != "none":
        bits.append(f"scrf=(smd,solvent={solvent})")
    bits.append("int=ultrafine")
    if extra:
        bits.extend(extra.split())
    return " ".join(bits)


def header(name, nproc, mem, chk, oldchk=None):
    lines = [f"%mem={mem}", f"%nprocshared={nproc}"]
    if oldchk:
        lines.append(f"%oldchk={oldchk}")
    lines.append(f"%chk={chk}")
    return "\n".join(lines)


def geom_section(atoms, charge, mult):
    return f"\n{charge} {mult}\n{xyz_repr(atoms)}\n"


def modredundant(cons):
    return "\n".join(f"B {i} {j} {d:.3f} F" for i, j, d in cons)


def make_l1(atoms, args, cons):
    # modredundant input follows the geometry, separated by a blank line
    return (
        header("01_const_optfreq", args.nproc, args.mem, "01_const_optfreq.chk")
        + "\n"
        + route(args.method, args.basis, args.disp, args.solvent, args.extra_route,
                f"opt=(calcfc,modredundant,maxcycles={args.maxiter}) freq")
        + "\n\nTS step 1: constrained Opt+Freq (reference point, not a TS)\n"
        + geom_section(atoms, args.charge, args.mult)
        + "\n"
        + modredundant(cons)
        + "\n\n"
    )


def make_l2(args):
    ts_tail = (f"opt=(ts,readfc,noeigen,maxcycles={args.maxiter}) freq"
               if args.ts_route == "readfc" else
               f"opt=(ts,calcfc,noeigen,maxcycles={args.maxiter}) freq")
    return (
        header("02_ts_optfreq", args.nproc, args.mem, "02_ts_optfreq.chk",
               oldchk="01_const_optfreq.chk")
        + "\n"
        + route(args.method, args.basis, args.disp, args.solvent, args.extra_route,
                ts_tail + " geom=check guess=read")
        + "\n\nTS step 2: OptTS+Freq, Hessian from step 1\n"
        + f"\n{args.charge} {args.mult}\n\n"
    )


def make_irc(args):
    return (
        header("03_irc", args.nproc, args.mem, "03_irc.chk",
               oldchk="02_ts_optfreq.chk")
        + "\n"
        + route(args.method, args.basis, args.disp, args.solvent, args.extra_route,
                f"irc=(readfc,maxpoints={args.irc_points},stepsize={args.irc_step}) "
                "geom=check guess=read")
        + "\n\nIRC both directions from the TS geometry\n"
        + f"\n{args.charge} {args.mult}\n\n"
    )


def make_sp(args):
    return (
        header("04_sp", args.nproc, args.mem, "04_sp.chk",
               oldchk="02_ts_optfreq.chk")
        + "\n"
        + route(args.sp_method, args.sp_basis, args.sp_disp, args.solvent,
                args.extra_route, "geom=check guess=read")
        + "\n\nHigher-level single point at the TS geometry\n"
        + f"\n{args.charge} {args.mult}\n\n"
    )


def run_chain(args, has_irc, has_sp):
    lines = [
        "#!/bin/bash",
        "set -u",
        'cd "$(dirname "$0")"',
        'G16="${G16_EXE:-g16}"',
        'export GAUSS_SCRDIR="${GAUSS_SCRDIR:-/tmp}"',
        "stamp() { echo \"$(date +%Y-%m-%dT%H:%M:%S)|$1|$2|exit=$3\" >> chain.events.log; }",
        "",
    ]
    steps = [("01_const_optfreq", [])]
    steps.append(("02_ts_optfreq", ["01_const_optfreq.chk"]))
    if has_irc:
        steps.append(("03_irc", ["02_ts_optfreq.chk"]))
    if has_sp:
        steps.append(("04_sp", ["02_ts_optfreq.chk"]))
    for step, requires in steps:
        lines.append(f'if [ ! -f "{step}.gjf" ]; then '
                     f'stamp {step} SKIP-MISSING-INPUT 1; exit 1; fi')
        for req in requires:
            lines.append(f'if [ ! -f "{req}" ]; then '
                         f'stamp {step} MISSING-{req} 1; '
                         f'echo "required file {req} not found - stop"; exit 1; fi')
        lines += [
            f'echo "=== {step} ==="',
            f'"$G16" {step}.gjf > {step}.log 2>&1',
            'rc=$?',
            f'if ! grep -q "Normal termination" {step}.log; then',
            f'  stamp {step} FAIL $rc; echo "{step} FAILED"; exit 1;',
            "fi",
            f'stamp {step} DONE $rc',
            "",
        ]
    lines += ['echo CHAIN-ALL-DONE > chain-state.txt',
              'stamp CHAIN ALL-DONE 0',
              'echo "chain complete"']
    return "\n".join(lines) + "\n"


def slurm(args):
    if args.g16root:
        profile = (f"export g16root='{args.g16root}'\n"
                   'source "${g16root}/g16/bsd/g16.profile"')
    else:
        profile = ("# Set g16root/source here if Gaussian is not on PATH.\n"
                   "# export g16root='<your g16root>'\n"
                   '# source "${g16root}/g16/bsd/g16.profile"')
    return f"""#!/bin/bash
#SBATCH -J {args.name}
#SBATCH -N 1
#SBATCH -n {args.nproc}
#SBATCH -c 1
#SBATCH -p {args.partition}
#SBATCH -o slurm-%j.out
#SBATCH --time={args.time}

export PATH=/bin:/usr/bin:/usr/local/bin:/sbin:/usr/sbin:/usr/local/sbin
export LD_LIBRARY_PATH=/lib64:/usr/lib64:/usr/local/lib64:/usr/lib64/atlas
{profile}
export GAUSS_SCRDIR={args.scratch}
export G16_EXE="${{G16_EXE:-g16}}"

date
bash run_chain.sh
echo "EXIT=$?"
date
"""


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xyz", required=True, help="starting geometry (first frame)")
    ap.add_argument("--name", required=True, help="job name prefix")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--mult", type=int, default=1)
    ap.add_argument("--method", default="UB3LYP", help="e.g. UB3LYP, UM06-2X, B3LYP")
    ap.add_argument("--basis", default="6-31G(d)", help='e.g. "6-31G(d)", def2-SVP')
    ap.add_argument("--disp", default="em=gd3bj", help="em=gd3bj, em=gd3bj, or none")
    ap.add_argument("--solvent", default="toluene", help="SMD solvent, or none")
    ap.add_argument("--extra-route", default="", help="extra route keywords")
    ap.add_argument("--constraints", default="",
                    help='L1 fixed distances: "i,j,d;i,j,d" (1-based)')
    ap.add_argument("--maxiter", type=int, default=220)
    ap.add_argument("--ts-route", choices=["readfc", "calcfc"], default="readfc",
                    help="readfc uses step-1 Hessian; calcfc recomputes it")
    ap.add_argument("--irc-points", type=int, default=50)
    ap.add_argument("--irc-step", type=int, default=10)
    ap.add_argument("--no-irc", action="store_true")
    ap.add_argument("--no-sp", action="store_true")
    ap.add_argument("--sp-method", default="UM06-2X")
    ap.add_argument("--sp-basis", default="6-311+G(d,p)")
    ap.add_argument("--sp-disp", default="em=gd3bj")
    ap.add_argument("--nproc", type=int, default=24)
    ap.add_argument("--mem", default="48GB")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--slurm", action="store_true", help="write submit.sbatch")
    ap.add_argument("--partition", default="community")
    ap.add_argument("--time", default="48:00:00")
    ap.add_argument("--g16root", default="",
                    help="Gaussian installation root for the sbatch template; "
                         "leave empty if g16 is on PATH")
    ap.add_argument("--scratch", default="/tmp",
                    help="node-local GAUSS_SCRDIR in the sbatch template")
    args = ap.parse_args()

    frames = read_xyz_frames(args.xyz)
    if not frames:
        sys.exit(f"cannot read geometry: {args.xyz}")
    atoms = frames[0]
    cons = parse_constraints(args.constraints)
    if not cons:
        sys.exit("give at least one L1 constraint, e.g. --constraints '1,2,1.500'")
    n = len(atoms)
    for i, j, _ in cons:
        if not (1 <= i <= n and 1 <= j <= n):
            sys.exit(f"constraint {i},{j} out of range (n={n})")

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    files = {
        "01_const_optfreq.gjf": make_l1(atoms, args, cons),
        "02_ts_optfreq.gjf": make_l2(args),
    }
    if not args.no_irc:
        files["03_irc.gjf"] = make_irc(args)
    if not args.no_sp:
        files["04_sp.gjf"] = make_sp(args)
    for name, text in files.items():
        (out / name).write_text(text, encoding="utf-8")
    # remove stale steps from a previous invocation with different flags
    if args.no_irc:
        (out / "03_irc.gjf").unlink(missing_ok=True)
    if args.no_sp:
        (out / "04_sp.gjf").unlink(missing_ok=True)
    if not args.slurm:
        (out / "submit.sbatch").unlink(missing_ok=True)
    (out / "run_chain.sh").write_text(
        run_chain(args, not args.no_irc, not args.no_sp), encoding="utf-8")
    if args.slurm:
        (out / "submit.sbatch").write_text(slurm(args), encoding="utf-8")

    print(f"generated Gaussian ladder in {out}/")
    for name in files:
        print(f"  {name}")
    print(f"  run_chain.sh" + ("  submit.sbatch" if args.slurm else ""))
    print(f"\nL1 constraints (1-based):")
    for i, j, d in cons:
        print(f"  B {i} {j} {d:.3f} F")
    print("\nRun locally:   bash run_chain.sh")
    if args.slurm:
        print("Submit:        sbatch submit.sbatch")
        if not args.g16root:
            print("WARNING: no --g16root; edit submit.sbatch or ensure g16 is on PATH.")
    print("Verify:        python gauss_verify.py --out 02_ts_optfreq.log")
    return 0


if __name__ == "__main__":
    sys.exit(main())

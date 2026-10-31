#!/usr/bin/env python3
"""Turn candidate TS geometries (xyz) into ready-to-run ORCA input files.

This is the ORCA counterpart of the UniTS repository's `units-xyz2gjf`, which
emits Gaussian input. Two task layouts are offered, mirroring the paper's
`direct_ts` and `two_step_ts`:

  direct    one job:  OptTS + Freq from the guess as-is.
            Cheapest, and what a good ML guess is meant to support.

  two_step  the field-verified protocol for floppy or weakly-bound systems:
              01_cons_optfreq  constrained Opt + Freq (hard distance constraints
                               hold the reactive bonds) -> yields a .hess file
              02_optts         OptTS + Freq started from that point WITH
                               `InHess Read`, so the optimizer follows the real
                               transition mode instead of a model-Hessian
                               lowest mode
              03_irc           IRC from the TS Hessian, both directions
            Slower, but it is the layout that survives real ligand frameworks.

Index convention (a classic silent error):
  `--reactive` takes 1-BASED atom numbers, as everywhere else in this skill and
  as chemists count them. ORCA's own coordinate block and `{ B i j ... }`
  constraint entries are 0-BASED, so the script converts and prints the mapping
  together with the element symbols, e.g.
      declared 21-57  ->  Ni21-H57  ->  ORCA { B 20 56 ... }
  Always read that line back before submitting.

Usage:
  python ts_orca_prep.py --xyz pass.xyz --reactive 21,57 --charge 0 --mult 2 \
      --outdir ts_jobs --mode two_step
  python ts_orca_prep.py --xyz pass.xyz --reactive 21,57 --method TPSS \
      --basis def2-SVP --disp D3BJ --solvent Toluene --mode direct
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ts_common import (  # noqa: E402
    formula, read_xyz_frames, xyz_repr,
)

# Elements whose treatment the paper's default level handles poorly.
D_BLOCK = {"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
           "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
           "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"}


def route(method, basis, disp, extra=()):
    """Assemble the ORCA simple-input line."""
    parts = ["!"]
    for kw in (method, disp, basis, "RIJCOSX", "TightSCF"):
        if kw and kw.lower() != "none":
            parts.append(kw)
    parts.extend(extra)
    return " ".join(parts)


def solvent_block(solvent, smd=True):
    if not solvent:
        return ""
    if smd:
        return (f'\n%cpcm\n   smd true\n   smdsolvent "{solvent}"\nend\n')
    return f'\n%cpcm\n   solvent "{solvent}"\nend\n'


def pal_block(nprocs, maxcore):
    return f"\n%pal nprocs {nprocs} end\n%maxcore {maxcore}\n"


def coord_block(atoms, charge, mult):
    return (f"\n* xyz {charge} {mult}\n{xyz_repr(atoms)}\n*\n")


def geom_block(maxiter, constraints=None, inhess=None):
    """%geom block; `constraints` is a list of (i0, j0, distance, 0based)."""
    lines = ["\n%geom", f"  MaxIter {maxiter}"]
    if inhess:
        lines.append("  InHess Read")
        lines.append(f'  InHessName "{inhess}"')
    if constraints:
        lines.append("  Constraints")
        for i0, j0, d, _ in constraints:
            lines.append(f"    {{ B {i0} {j0} {d:.3f} C }}")
        lines.append("  end")
    lines.append("end\n")
    return "\n".join(lines)


def build_direct(atoms, args, pairs0):
    return (route(args.method, args.basis, args.disp, ("OptTS", "Freq"))
            + solvent_block(args.solvent)
            + pal_block(args.nprocs, args.maxcore)
            + geom_block(args.maxiter)
            + coord_block(atoms, args.charge, args.mult))


def build_cons_optfreq(atoms, args, cons):
    return (route(args.method, args.basis, args.disp, ("Opt", "Freq"))
            + solvent_block(args.solvent)
            + pal_block(args.nprocs, args.maxcore)
            + geom_block(args.maxiter, constraints=cons)
            + coord_block(atoms, args.charge, args.mult))


def build_optts(atoms, args, hess_name):
    return (route(args.method, args.basis, args.disp, ("OptTS", "Freq"))
            + solvent_block(args.solvent)
            + pal_block(args.nprocs, args.maxcore)
            + geom_block(args.maxiter, inhess=hess_name)
            + coord_block(atoms, args.charge, args.mult))


def build_irc(atoms, args, hess_name, geom_file):
    # Verified against the ORCA 6.1.1 binary: the [IRC] block reads
    # InitHess / Hess_Filename / Direction / MaxIter / Trajectory; IRC requires
    # a Hessian file and cannot be combined with internal constraints.
    #
    # The coordinates are read from the previous step's optimised geometry
    # (ORCA writes <job>.xyz when an optimisation finishes), NOT from the
    # starting guess: an IRC must begin at the converged TS, and the Hessian
    # passed in must be the one computed at that same geometry.
    return (route(args.method, args.basis, args.disp, ("IRC",))
            + solvent_block(args.solvent)
            + pal_block(args.nprocs, args.maxcore)
            + "\n%irc\n  MaxIter 20\n  InitHess Read\n"
              f'  Hess_Filename "{hess_name}"\n  Direction both\n'
              "  PrintLevel 2\n  Trajectory True\nend\n"
            + f"\n* xyzfile {args.charge} {args.mult} {geom_file}\n")


def run_chain_script(jobname, steps, preflight=None):
    """Sequential per-candidate chain; never run two ORCA jobs at once.

    `steps` is a list of (jobname, [required files]) so a missing prerequisite
    (the .hess / .xyz produced by the previous step) aborts the chain instead
    of letting ORCA run with stale or absent inputs.
    """
    lines = ["#!/bin/bash",
             "# generated by ts_orca_prep.py - run jobs strictly in sequence",
             "cd \"$(dirname \"$0\")\"",
             'export PATH="${ORCA_DIR:+$ORCA_DIR:}$PATH"',
             'if ! command -v orca >/dev/null 2>&1; then',
             '  echo "orca not found: set ORCA_DIR or put orca on PATH" >&2; exit 1',
             'fi',
             f"LOG={jobname}_chain.events.log",
             "stamp() { echo \"$(date +%Y-%m-%dT%H:%M:%S)|$1|$2|exit=$3\" >> $LOG; }",
             ""]
    first = steps[0][0]
    if preflight:
        lines += [
            f'python "{preflight}" --job {first} --dir . >/dev/null 2>&1',
            "if [ $? -ne 0 ]; then stamp PREFLIGHT FAIL 1; "
            "echo 'preflight FAIL - not launching'; exit 1; fi",
            "stamp PREFLIGHT OK 0",
            ""]
    for step, requires in steps:
        lines.append(f'if [ ! -f "{step}.inp" ]; then '
                     f'stamp {step} SKIP-MISSING-INPUT 1; exit 1; fi')
        for req in requires:
            lines.append(f'if [ ! -f "{req}" ]; then '
                         f'stamp {step} MISSING-{req} 1; '
                         f'echo "required file {req} not found - stop"; exit 1; fi')
        lines += [
            f'echo "=== {step} ==="',
            f'orca {step}.inp > {step}.out 2>&1',
            'rc=$?',
            f'if ! grep -q "ORCA TERMINATED NORMALLY" {step}.out; then',
            f'  stamp {step} FAIL $rc; echo "{step} FAILED"; exit 1;',
            "fi",
            f'stamp {step} DONE $rc',
            ""]
    lines += ['echo CHAIN-ALL-DONE > chain-state.txt',
              'stamp CHAIN ALL-DONE 0',
              'echo "chain complete"']
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xyz", required=True,
                    help="candidate geometry file (multi-frame allowed)")
    ap.add_argument("--reactive", action="append", required=True, metavar="I,J",
                    help="reactive atom pair(s) to constrain, 1-BASED, repeatable")
    ap.add_argument("--set-distance", action="append", default=None,
                    metavar="I,J,D", help="override the constraint distance for a "
                                          "pair, 1-based, repeatable")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--mult", type=int, default=1)
    ap.add_argument("--method", default="B3LYP")
    ap.add_argument("--basis", default="def2-SVP")
    ap.add_argument("--disp", default="D3BJ",
                    help="dispersion keyword (D3BJ, D4, ... or 'none')")
    ap.add_argument("--solvent", default=None,
                    help='SMD solvent name, e.g. "Toluene" (default: gas phase)')
    ap.add_argument("--mode", choices=["direct", "two_step"], default="two_step")
    ap.add_argument("--maxiter", type=int, default=300)
    ap.add_argument("--nprocs", type=int, default=12)
    ap.add_argument("--maxcore", type=int, default=1100)
    ap.add_argument("--outdir", default="ts_jobs")
    ap.add_argument("--only", default=None,
                    help="restrict to these candidate indices, e.g. 0,3,5")
    ap.add_argument("--preflight", default=None,
                    help="path to preflight_check.py (bundled in scripts/) to "
                         "wire into the generated runner")
    args = ap.parse_args()

    try:
        frames = read_xyz_frames(args.xyz)
    except ValueError as exc:
        sys.exit(f"cannot read --xyz: {exc}")
    if not frames:
        sys.exit(f"no xyz frames read from {args.xyz}")

    pairs = []
    for spec in args.reactive:
        parts = [int(x) for x in spec.replace(" ", ",").split(",") if x]
        if len(parts) != 2:
            sys.exit(f"--reactive needs an I,J pair, got {spec!r}")
        if min(parts) < 1:
            sys.exit(f"--reactive is 1-based; got {spec!r}")
        pairs.append(tuple(parts))

    overrides = {}
    for spec in (args.set_distance or []):
        parts = [float(x) for x in spec.replace(" ", ",").split(",") if x]
        if len(parts) != 3:
            sys.exit(f"--set-distance needs I,J,D; got {spec!r}")
        overrides[(int(parts[0]), int(parts[1]))] = parts[2]

    want = None
    if args.only:
        want = {int(x) for x in args.only.replace(" ", ",").split(",") if x}

    n_atoms = len(frames[0])
    for i, f in enumerate(frames):
        if len(f) != n_atoms:
            sys.exit(f"frame {i} has {len(f)} atoms, frame 0 has {n_atoms}; "
                     "all candidates must describe the same system")

    # --- report the index mapping loudly: this is where silent errors happen ---
    print(f"candidate file : {args.xyz}  ({len(frames)} frame(s), "
          f"{n_atoms} atoms, {formula(frames[0])})")
    print(f"level of theory: {route(args.method, args.basis, args.disp).replace('! ', '')}"
          + (f'  SMD: {args.solvent}' if args.solvent else "  gas phase"))
    print("reactive pairs (1-based -> ORCA 0-based constraint):")
    for a, b in pairs:
        if not (1 <= a <= n_atoms and 1 <= b <= n_atoms):
            sys.exit(f"reactive pair {a}-{b} out of range for {n_atoms} atoms")
        sym_a = frames[0][a - 1][0]
        sym_b = frames[0][b - 1][0]
        print(f"    declared {a:>3}-{b:<3}  ->  {sym_a}{a}-{sym_b}{b}"
              f"  ->  {{ B {a - 1} {b - 1} <d> C }}")
    # transition-metal caution
    tms = sorted({s for s, _ in frames[0] if s in D_BLOCK})
    if tms:
        print()
        print(f"  ! structure contains {', '.join(tms)}. The paper's dataset level "
              "(B3LYP-D3BJ/def2-SVP)\n    is a cheap level chosen for dataset "
              "consistency and is NOT the recommended\n    choice for transition-"
              "metal systems. Per the local ORCA functional\n    guide, "
              "use a pure/low-HF functional for optimisation (TPSS-D3BJ or\n    "
              "r2SCAN-3c) and a higher single-point level (TPSS0-D4 or PBE0-D4 "
              "with\n    def2-TZVP). Optimisation, frequency and IRC must all stay "
              "at one level.")

    os.makedirs(args.outdir, exist_ok=True)
    made = []
    for idx, atoms in enumerate(frames):
        if want is not None and idx not in want:
            continue
        tag = f"cand{idx:02d}"
        d = os.path.join(args.outdir, tag)
        os.makedirs(d, exist_ok=True)

        # constraint list with distances taken from THIS candidate's geometry
        cons = []
        for a, b in pairs:
            import math
            dist = math.dist(atoms[a - 1][1], atoms[b - 1][1])
            key = (a, b)
            if key not in overrides and (b, a) in overrides:
                key = (b, a)
            if key in overrides:
                dist = overrides[key]
            cons.append((a - 1, b - 1, dist, True))

        with open(os.path.join(d, f"{tag}.xyz"), "w", encoding="utf-8") as fh:
            fh.write(f"{len(atoms)}\n{tag}\n")
            for sym, (x, y, z) in atoms:
                fh.write(f"{sym:<2} {x:16.6f} {y:16.6f} {z:16.6f}\n")

        if args.mode == "direct":
            steps = [("01_optts", [])]
            with open(os.path.join(d, "01_optts.inp"), "w", encoding="utf-8") as fh:
                fh.write(build_direct(atoms, args, pairs))
        else:
            steps = [("01_cons_optfreq", []),
                     ("02_optts", ["01_cons_optfreq.hess"]),
                     ("03_irc", ["02_optts.hess", "02_optts.xyz"])]
            with open(os.path.join(d, "01_cons_optfreq.inp"), "w",
                      encoding="utf-8") as fh:
                fh.write(build_cons_optfreq(atoms, args, cons))
            with open(os.path.join(d, "02_optts.inp"), "w", encoding="utf-8") as fh:
                fh.write(build_optts(atoms, args, "01_cons_optfreq.hess"))
            with open(os.path.join(d, "03_irc.inp"), "w", encoding="utf-8") as fh:
                fh.write(build_irc(atoms, args, "02_optts.hess", "02_optts.xyz"))

        with open(os.path.join(d, "run_chain.sh"), "w", encoding="utf-8") as fh:
            fh.write(run_chain_script(tag, steps, args.preflight))
        made.append((tag, d, steps, cons))

    # top-level sequential runner over all candidates
    with open(os.path.join(args.outdir, "run_all.sh"), "w", encoding="utf-8") as fh:
        fh.write("#!/bin/bash\n"
                 "# runs candidates strictly one after another: never two ORCA\n"
                 "# jobs in parallel, and never interrupt a running job.\n"
                 "cd \"$(dirname \"$0\")\"\n")
        for tag, _, _, _ in made:
            fh.write(f'bash {tag}/run_chain.sh\n')
        fh.write('echo ALL-CANDIDATES-DONE > all-done.txt\n')
    print()
    print(f"wrote {len(made)} candidate job set(s) under {args.outdir}/")
    for tag, d, steps, cons in made:
        cs = ", ".join(f"{i + 1}-{j + 1}@{dd:.3f}A" for i, j, dd, _ in cons)
        order = " -> ".join(s for s, _ in steps)
        print(f"  {tag}: {order}   constraints: {cs}")
    print()
    print(f"run everything in sequence:  bash {args.outdir}/run_all.sh")
    return 0


if __name__ == "__main__":
    sys.exit(main())

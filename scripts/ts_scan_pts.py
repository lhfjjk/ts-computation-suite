#!/usr/bin/env python3
"""Per-point independent constrained scan for ORCA (no %geom Scan).

Why this exists
---------------
ORCA 6.1.1 can fail an entire relaxed scan if *any* scan point hits MaxIter:
the archive step then fails repeatedly and the remaining points never run.
Running each point as its own `! Opt` job, starting from the previous
converged geometry, makes a failed point local instead of fatal.

The output is a barrier profile, not a TS. Use the highest converged point
as the starting geometry for the L1/L2 TS ladder (ts_ladder.py).

Examples
--------
# generate only
python ts_scan_pts.py --xyz reactant.xyz --pair 59,75 --start 2.05 --end 4.25 \
    --points 8 --charge 0 --mult 2 --method TPSS --basis def2-SVP \
    --disp D3BJ --solvent Toluene \
    --constraints "59,5,2.12;59,7,2.08;59,60,2.42;59,81,3.40" \
    --outdir scan_r2

# generate and run sequentially (one ORCA at a time)
python ts_scan_pts.py ... --outdir scan_r2 --run

# re-parse an existing scan
python ts_scan_pts.py --summary scan_r2

Notes
-----
* `--pair`, `--constraints` are 1-based (converted to ORCA's 0-based).
* Only points whose output contains THE OPTIMIZATION HAS CONVERGED are used
  as the next starting geometry; failed points are recorded and skipped.
* The reaction-coordinate constraint is `{ B i j d C }`; add spectator
  constraints with `--constraints`.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ts_common import read_xyz_frames, write_xyz  # noqa: E402

HARTREE = 627.509474


def parse_pair(spec: str):
    a, b = [int(x) for x in spec.replace(" ", ",").split(",")]
    return a, b


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


def route(method: str, basis: str, disp: str, extra: str, rijcosx: bool = True):
    bits = ["! Opt", method]
    if disp and disp.lower() != "none":
        bits.append(disp)
    bits.append(basis)
    if rijcosx:
        bits.append("RIJCOSX")
    if extra:
        bits.extend(extra.split())
    return " ".join(bits)


def solvent_block(solvent: str | None):
    if not solvent:
        return ""
    return f"\n%cpcm\n  smd true\n  SMDsolvent \"{solvent}\"\nend\n"


def geom_block(maxiter: int, cons):
    lines = ["\n%geom", f"  MaxIter {maxiter}", "  Constraints"]
    for i, j, d in cons:
        lines.append(f"    {{ B {i - 1} {j - 1} {d:.3f} C }}")
    lines.extend(["  end", "end\n"])
    return "\n".join(lines)


def make_input(atoms, args, cons, geom_name: str):
    n = len(atoms)
    text = (
        route(args.method, args.basis, args.disp, args.extra, not args.no_rijcosx)
        + solvent_block(args.solvent)
        + f"\n%pal nprocs {args.nprocs} end\n%maxcore {args.maxcore}\n"
        + geom_block(args.maxiter, cons)
        + f"\n* xyzfile {args.charge} {args.mult} {geom_name}\n"
    )
    return text


def last_energy(text: str):
    vals = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", text)
    return float(vals[-1]) if vals else None


def last_geometry(text: str):
    """Return [(element, x, y, z)] from the last CARTESIAN COORDINATES block."""
    blocks = list(re.finditer(r"CARTESIAN COORDINATES \(ANGSTROEM\)", text))
    if not blocks:
        return None
    start = blocks[-1].end()
    lines = text[start:].splitlines()
    atoms = []
    for line in lines:
        if not line.strip():
            if atoms:
                break
            continue
        parts = line.split()
        if len(parts) >= 4 and re.match(r"^[A-Za-z]{1,3}$", parts[0]):
            try:
                atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
            except ValueError:
                break
        elif atoms:
            break
    return atoms or None


def distance(atoms, i, j):
    if not atoms or not (1 <= i <= len(atoms) and 1 <= j <= len(atoms)):
        return None
    a, b = atoms[i - 1], atoms[j - 1]
    return math.sqrt((a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2 + (a[3] - b[3]) ** 2)


def converged(text: str):
    return "THE OPTIMIZATION HAS CONVERGED" in text


def normal(text: str):
    return "ORCA TERMINATED NORMALLY" in text


def read_point(point_dir: Path):
    out = point_dir / "point.out"
    if not out.exists():
        return None
    text = out.read_text(errors="replace")
    return {
        "e": last_energy(text),
        "conv": converged(text),
        "normal": normal(text),
        "geom": last_geometry(text),
    }


def summarize(outdir: Path):
    rows = []
    for p in sorted(outdir.glob("p*")):
        if not p.is_dir():
            continue
        meta = {}
        mf = p / "point.json"
        if mf.exists():
            try:
                meta = json.loads(mf.read_text(encoding="utf-8"))
            except Exception:
                meta = {}
        r = read_point(p)
        if r is None:
            rows.append({"point": p.name, "d": meta.get("d"), "e": None,
                         "conv": False, "note": "missing"})
            continue
        rows.append({
            "point": p.name,
            "d": meta.get("d"),
            "e": r["e"],
            "conv": r["conv"],
            "normal": r["normal"],
            "d_actual": meta.get("pair"),
            "note": "" if r["conv"] else "not-converged",
        })
    # actual distance from the last geometry
    for row in rows:
        p = outdir / row["point"]
        r = read_point(p)
        meta = {}
        mf = p / "point.json"
        if mf.exists():
            try:
                meta = json.loads(mf.read_text(encoding="utf-8"))
            except Exception:
                pass
        pair = meta.get("pair")
        if r and r.get("geom") and pair:
            row["d_actual"] = distance(r["geom"], pair[0], pair[1])

    # relative energy: first converged point if available, else first point
    ref_row = next((r for r in rows if r["conv"] and r["e"] is not None), None)
    if ref_row is None:
        ref_row = next((r for r in rows if r["e"] is not None), None)
    ref = ref_row["e"] if ref_row else None
    ref_point = ref_row["point"] if ref_row else ""
    for row in rows:
        row["e_rel"] = (row["e"] - ref) * HARTREE if (row["e"] is not None and ref is not None) else None
        row["ref"] = ref_point

    tsv = outdir / "scan_summary.tsv"
    with tsv.open("w", encoding="utf-8") as fh:
        fh.write("point\td_requested\td_actual\tE_h\tE_rel_kcal\tconv\tref\tnote\n")
        for r in rows:
            def fmt(x, nd=6):
                return "" if x is None else (f"{x:.{nd}f}" if isinstance(x, float) else str(x))
            fh.write("\t".join([
                r["point"], fmt(r["d"], 3), fmt(r.get("d_actual"), 3),
                fmt(r["e"], 8), fmt(r.get("e_rel"), 2), "Y" if r["conv"] else "N",
                r.get("ref", ""), r.get("note", ""),
            ]) + "\n")
    (outdir / "scan_summary.json").write_text(
        json.dumps(rows, indent=2, ensure_ascii=False), encoding="utf-8")

    n_conv = sum(1 for r in rows if r["conv"])
    state = "SCAN-DONE" if n_conv == len(rows) and rows else (
        "SCAN-PARTIAL" if n_conv else "SCAN-FAILED")
    (outdir / "scan_state.txt").write_text(
        f"{state}\nconverged={n_conv}/{len(rows)}\n", encoding="utf-8")
    return rows, state


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xyz", help="starting geometry (first frame is used)")
    ap.add_argument("--pair", help="reaction-coordinate atom pair i,j (1-based)")
    ap.add_argument("--start", type=float, help="start distance (A)")
    ap.add_argument("--end", type=float, help="end distance (A)")
    ap.add_argument("--points", type=int, default=8, help="number of points")
    ap.add_argument("--constraints", default="",
                    help='extra fixed distances: "i,j,d;i,j,d" (1-based)')
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--mult", type=int, default=1)
    ap.add_argument("--method", default="TPSS")
    ap.add_argument("--basis", default="def2-SVP")
    ap.add_argument("--disp", default="D3BJ", help="D3BJ, D4, or none")
    ap.add_argument("--solvent", default=None, help="SMD solvent name, e.g. Toluene")
    ap.add_argument("--extra", default="", help="extra route keywords")
    ap.add_argument("--no-rijcosx", action="store_true",
                    help="omit RIJCOSX (for methods/bases where it is unavailable)")
    ap.add_argument("--maxiter", type=int, default=220)
    ap.add_argument("--nprocs", type=int, default=12)
    ap.add_argument("--maxcore", type=int, default=1100)
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--orca", default=os.environ.get("ORCA_EXE", ""),
                    help="ORCA executable (default: $ORCA_EXE, then orca on PATH)")
    ap.add_argument("--run", action="store_true", help="run the points sequentially")
    ap.add_argument("--resume", action="store_true", help="skip converged points")
    ap.add_argument("--summary", metavar="DIR", help="only re-parse an existing scan")
    ap.add_argument("--timeout", type=int, default=0, help="per-point timeout (s), 0=none")
    args = ap.parse_args()

    if args.summary:
        rows, state = summarize(Path(args.summary))
        print(f"{state}: {len(rows)} points -> {Path(args.summary) / 'scan_summary.tsv'}")
        return 0

    if not args.outdir:
        ap.error("--outdir is required unless --summary is used")
    if not (args.xyz and args.pair and args.start is not None and args.end is not None):
        ap.error("--xyz, --pair, --start, --end are required unless --summary is used")
    if args.points < 2:
        ap.error("--points must be >= 2")

    atoms = read_xyz_frames(args.xyz, last_only=False)
    if not atoms:
        sys.exit(f"cannot read geometry: {args.xyz}")
    atoms = atoms[0]
    n = len(atoms)
    i, j = parse_pair(args.pair)
    if not (1 <= i <= n and 1 <= j <= n):
        sys.exit(f"--pair {i},{j} out of range (n={n})")

    cons_extra = parse_constraints(args.constraints)
    for a, b, _ in cons_extra:
        if not (1 <= a <= n and 1 <= b <= n):
            sys.exit(f"--constraints {a},{b} out of range (n={n})")

    ds = [args.start + (args.end - args.start) * k / (args.points - 1)
          for k in range(args.points)]
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # warn if the starting geometry does not match the requested start distance
    from_start = distance([(a[0], a[1][0], a[1][1], a[1][2]) for a in atoms], i, j)
    if from_start is not None and abs(from_start - args.start) > 0.6:
        print(f"WARNING: start geometry has d({i},{j})={from_start:.3f} A, "
              f"but --start={args.start:.3f}. The first point will be constrained to "
              f"{args.start:.3f}; check that this is intentional.")

    for k, d in enumerate(ds, 1):
        pdir = outdir / f"p{k:02d}"
        pdir.mkdir(parents=True, exist_ok=True)
        cons = [(i, j, d)] + cons_extra
        geom_name = "start.xyz" if k == 1 else "prev.xyz"
        (pdir / "point.inp").write_text(
            make_input(atoms, args, cons, geom_name), encoding="utf-8")
        if k == 1:
            write_xyz(pdir / "start.xyz", atoms, "scan start")
        (pdir / "point.json").write_text(
            json.dumps({"d": d, "pair": [i, j], "charge": args.charge,
                        "mult": args.mult}, indent=2), encoding="utf-8")

    # convenience shell runner (the Python --run path is primary)
    sh = outdir / "run_scan.sh"
    this_script = Path(__file__).resolve()
    lines = ["#!/bin/bash", "set -u", 'cd "$(dirname "$0")"',
             'ORCA="${ORCA_EXE:-orca}"', "prev=''", ""]
    for k in range(1, len(ds) + 1):
        p = f"p{k:02d}"
        lines.append(f"# --- {p}: d={ds[k-1]:.3f} A ---")
        lines.append(f'if [ -z "$prev" ]; then cp p01/start.xyz {p}/prev.xyz; '
                     f'else cp "$prev/point.xyz" {p}/prev.xyz; fi')
        lines.append(f'(cd {p} && "$ORCA" point.inp > point.out 2>&1)')
        lines.append(f'if grep -q "ORCA TERMINATED NORMALLY" {p}/point.out && '
                     f'grep -q "THE OPTIMIZATION HAS CONVERGED" {p}/point.out; then')
        lines.append(f'  prev={p}; echo "{p} OK"')
        lines.append('else')
        lines.append(f'  echo "{p} FAILED (kept previous geometry)"')
        lines.append('fi')
        lines.append("")
    lines.append(f'python "{this_script}" --summary .')
    sh.write_text("\n".join(lines) + "\n", encoding="utf-8")
    try:
        sh.chmod(0o755)
    except OSError:
        pass

    print(f"generated {len(ds)} points in {outdir}/")
    print(f"  d = {ds[0]:.3f} ... {ds[-1]:.3f} A, pair {i}-{j}, {n} atoms")
    print(f"  run: python {Path(__file__).name} --summary {outdir}   (parse only)")
    print(f"       bash {sh}                                      (sequential ORCA)")

    if not args.run:
        return 0

    orca = args.orca or shutil.which("orca")
    if not orca:
        sys.exit("ORCA executable not found: pass --orca or set ORCA_EXE")
    print(f"using ORCA: {orca}")
    first_start = outdir / "p01" / "start.xyz"
    prev_dir = None
    for k, d in enumerate(ds, 1):
        pdir = outdir / f"p{k:02d}"
        out = pdir / "point.out"
        if args.resume and out.exists():
            r = read_point(pdir)
            if r and r["conv"]:
                print(f"{pdir.name}: already converged, skip")
                prev_dir = pdir
                continue
        if prev_dir is None:
            shutil.copyfile(first_start, pdir / "prev.xyz")
        else:
            shutil.copyfile(prev_dir / "point.xyz", pdir / "prev.xyz")
        print(f"{pdir.name}: d={d:.3f} ...", flush=True)
        with out.open("w", encoding="utf-8") as fh:
            try:
                subprocess.run([orca, "point.inp"], cwd=pdir, stdout=fh,
                               stderr=subprocess.STDOUT,
                               timeout=args.timeout if args.timeout else None)
            except subprocess.TimeoutExpired:
                print(f"{pdir.name}: TIMEOUT")
        r = read_point(pdir)
        if r and r["conv"] and r["normal"]:
            print(f"{pdir.name}: CONVERGED  E={r['e']}")
            prev_dir = pdir
        else:
            print(f"{pdir.name}: FAILED (next point uses {prev_dir.name if prev_dir else 'start'})")

    rows, state = summarize(outdir)
    print(f"{state}: {sum(1 for r in rows if r['conv'])}/{len(rows)} converged")
    print(f"summary: {outdir / 'scan_summary.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

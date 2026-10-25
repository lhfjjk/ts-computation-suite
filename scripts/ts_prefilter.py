#!/usr/bin/env python3
"""Screen candidate transition-state guesses before spending DFT time on them.

Designed for the output of UniTS-Gen (many independently sampled 3D TS guesses
from one 2D graph) and for hand-built guesses alike. The point is the same:
a single DFT OptTS/Freq run costs minutes to hours, so discard the physically
impossible candidates up front instead of discovering them afterwards.

Screening rules applied per candidate:
  * hard overlap        any interatomic distance < 1.00 A  -> FAIL
                        (the shortest real bond is C-H ~1.09 A)
  * structural clash    d < 0.75 * (rcov_i + rcov_j)        -> FAIL
                        (the clash criterion defined in the UniTS paper)
  * reactive pairs      distances of the declared forming/breaking bonds
                        -> reported, plus a WARN when a bond to be formed is
                        still far apart (> --far-cutoff, default 3.0 A)
  * RMSD vs reference   Kabsch-aligned RMSD to a known TS, if --ref is given
  * connectivity        bonds gained/lost vs the reference, if --ref is given
                        (guards against plausible-looking but rearranged guesses)
  * formula             identity check against --ref

Usage:
  python ts_prefilter.py --cand guesses.xyz --reactive 21,57
  python ts_prefilter.py --cand guesses.xyz --ref known_ts.xyz --reactive 21,57
  python ts_prefilter.py --cand guesses.xyz --ref ts.xyz --reactive 21,57 \
      --write-pass passing.xyz --json report.json
"""

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ts_common import (  # noqa: E402
    BOND_FACTOR, CLASH_FACTOR, COVALENT_RADII, HARD_OVERLAP,
    clash_report, connectivity_delta, formula, kabsch_rmsd, min_contact,
    reactive_distances, read_xyz_frames, write_xyz,
)


def screen(cand, ref=None, pairs=None, far_cutoff=3.0, forming_set=frozenset()):
    """Evaluate one candidate structure. Returns a result dict."""
    r = {"n_atoms": len(cand), "formula": formula(cand),
         "problems": [], "warnings": [], "info": {}}

    # 1. hard overlaps
    cl = clash_report(cand)
    r["info"]["n_hard_overlaps"] = cl["n_hard_overlaps"]
    r["info"]["n_clash_pairs"] = cl["n_clash_pairs"]
    if cl["n_hard_overlaps"]:
        ex = ", ".join(f"{i + 1}-{j + 1} {d:.3f}A" for i, j, d in cl["hard_overlaps"][:3])
        r["problems"].append(
            f"{cl['n_hard_overlaps']} interatomic distance(s) below {HARD_OVERLAP} A "
            f"({ex}) - physically impossible, discard")

    # 2. clashes per the paper's criterion
    if cl["n_clash_pairs"]:
        ex = ", ".join(f"{i + 1}-{j + 1} {d:.3f}A (limit {lim:.3f})"
                       for i, j, d, lim, _ in cl["clash_pairs"][:3])
        r["warnings"].append(
            f"{cl['n_clash_pairs']} structural clash pair(s) by the "
            f"d < {CLASH_FACTOR}*(r_i+r_j) criterion ({ex})")

    # 3. closest contact, to spot near-miss collapses
    mc = min_contact(cand)
    if mc:
        i, j, d, bonded = mc
        r["info"]["min_contact"] = {
            "atoms": [i + 1, j + 1], "distance_A": round(d, 3), "bonded": bonded,
            "elements": [cand[i][0], cand[j][0]],
        }
        if not bonded and d < 1.9:
            r["warnings"].append(
                f"closest non-bonded contact is {cand[i][0]}{i + 1}-{cand[j][0]}{j + 1} "
                f"at {d:.3f} A - inspect whether the optimizer will fuse it into a bond")

    # 4. reactive pair distances
    if pairs:
        rd = reactive_distances(cand, pairs)
        r["info"]["reactive_distances"] = [
            {"atoms": [a, b], "distance_A": round(d, 3),
             "role": "forming" if frozenset((a, b)) in forming_set else
                     ("breaking" if forming_set else "unspecified"),
             "elements": [cand[a - 1][0], cand[b - 1][0]]} for a, b, d in rd]
        for a, b, d in rd:
            # The upper bound only makes sense for a bond that is being FORMED:
            # a bond being broken is legitimately stretched at the TS (a leaving
            # group 3+ A away is normal). Warn only when the caller told us this
            # pair is forming, or when no role was given at all (then say so).
            is_forming = frozenset((a, b)) in forming_set
            if d > far_cutoff and (is_forming or not forming_set):
                role = "declared forming" if is_forming else "role not specified"
                r["warnings"].append(
                    f"reactive pair {a}-{b} is {d:.3f} A apart (> {far_cutoff} A, "
                    f"{role}) - a TS optimizer may not close this gap; consider "
                    "the two-step protocol with a distance constraint")
            elif d > far_cutoff:
                r["info"].setdefault("noted_far_pairs", []).append(
                    f"{a}-{b}:{d:.3f}A (declared breaking - expected at a TS)")

    # 5. reference comparison
    if ref is not None:
        if len(ref) != len(cand):
            r["warnings"].append(
                f"atom count differs from reference ({len(cand)} vs {len(ref)}) - "
                "RMSD and connectivity comparison skipped")
        else:
            if [s for s, _ in cand] != [s for s, _ in ref]:
                r["warnings"].append(
                    "atom ORDER or elements differ from the reference; RMSD is only "
                    "meaningful for identical orderings")
            else:
                rmsd_all = kabsch_rmsd(ref, cand)
                rmsd_heavy = kabsch_rmsd(ref, cand, heavy_only=True)
                r["info"]["rmsd_all_A"] = round(rmsd_all, 3)
                r["info"]["rmsd_heavy_A"] = round(rmsd_heavy, 3)
                gained, lost = connectivity_delta(ref, cand)
                r["info"]["bonds_gained"] = gained
                r["info"]["bonds_lost"] = lost
                if gained or lost:
                    g = f"gained {gained[:4]}" if gained else ""
                    l = f"lost {lost[:4]}" if lost else ""
                    r["warnings"].append(
                        "connectivity differs from the reference "
                        f"({g}{' ' if g and l else ''}{l}) - the guess may describe a "
                        "different reactive site")

    r["status"] = "FAIL" if r["problems"] else ("WARN" if r["warnings"] else "PASS")
    return r


def score(r):
    """Ranking key: clean structures first, then RMSD, then reactive distances.

    Deliberately lexicographic so a clean candidate always beats a flawed one.
    """
    info = r["info"]
    clean = 0 if r["status"] == "PASS" else (1 if r["status"] == "WARN" else 2)
    rmsd = info.get("rmsd_heavy_A", info.get("rmsd_all_A", None))
    # no reference: prefer the guess whose reactive bonds are closest to bonding
    if rmsd is None:
        rd = info.get("reactive_distances") or []
        spread = max((x["distance_A"] for x in rd), default=0.0)
    else:
        spread = rmsd
    return (clean, spread)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cand", required=True,
                    help="candidate xyz file (multi-frame allowed)")
    ap.add_argument("--ref", default=None,
                    help="reference TS xyz for RMSD / connectivity comparison")
    ap.add_argument("--reactive", action="append", default=None, metavar="I,J",
                    help="declared reactive atom pair, 1-based, repeatable "
                         "(the bond(s) forming/breaking)")
    ap.add_argument("--forming", action="append", default=None, metavar="I,J",
                    help="subset of --reactive pairs that are being FORMED, "
                         "1-based, repeatable. --far-cutoff is applied only to "
                         "these (a bond being broken is legitimately stretched "
                         "at a TS). Omit to apply the cutoff to every pair.")
    ap.add_argument("--far-cutoff", type=float, default=3.0,
                    help="warn when a bond to be FORMED is farther than this "
                         "(default 3.0 A)")
    ap.add_argument("--write-pass", default=None,
                    help="write the passing candidates to this xyz file")
    ap.add_argument("--json", default=None, help="write the full report as JSON")
    ap.add_argument("--top", type=int, default=10,
                    help="how many candidates to list in detail (default 10)")
    args = ap.parse_args()

    pairs = []
    for spec in (args.reactive or []):
        parts = [int(x) for x in spec.replace(" ", ",").split(",") if x]
        if len(parts) != 2:
            sys.exit(f"--reactive needs an I,J pair, got {spec!r}")
        pairs.append(tuple(parts))

    forming = set()
    for spec in (args.forming or []):
        parts = [int(x) for x in spec.replace(" ", ",").split(",") if x]
        if len(parts) != 2:
            sys.exit(f"--forming needs an I,J pair, got {spec!r}")
        if pairs and tuple(parts) not in pairs and tuple(reversed(parts)) not in pairs:
            sys.exit(f"--forming {spec!r} is not among the --reactive pairs")
        forming.add(frozenset(parts))
    forming = frozenset(forming)

    try:
        frames = read_xyz_frames(args.cand)
    except ValueError as exc:
        sys.exit(f"cannot read --cand: {exc}")
    if not frames:
        sys.exit(f"no xyz frames read from {args.cand}")
    ref = None
    if args.ref:
        try:
            rf = read_xyz_frames(args.ref, last_only=True)
        except ValueError as exc:
            sys.exit(f"cannot read --ref: {exc}")
        if not rf:
            sys.exit(f"no xyz frame read from {args.ref}")
        ref = rf[0]
        print(f"reference: {args.ref}  ({len(ref)} atoms, {formula(ref)})")
    print(f"candidates: {len(frames)} frame(s) from {args.cand}")
    if pairs:
        print("declared reactive pairs (1-based): "
              + ", ".join(f"{a}-{b}" for a, b in pairs)
              + ("   [forming: " + ", ".join(
                  f"{min(s)}-{max(s)}" for s in sorted(forming,
                                                       key=lambda x: min(x)))
                 + "]" if forming else "   [no --forming given: the far-cutoff "
                                       "is applied to every pair]"))
    print()

    results = []
    for idx, cand in enumerate(frames):
        r = screen(cand, ref, pairs, args.far_cutoff, forming)
        r["index"] = idx
        results.append(r)

    order = sorted(range(len(results)), key=lambda i: score(results[i]))
    print(f"{'rank':<5}{'cand':<6}{'status':<7}{'RMSD(heavy)':<13}"
          f"{'min d':<8}{'clash':<7}{'reactive dist (A)':<26}problems / warnings")
    for rank, i in enumerate(order[:args.top], 1):
        r = results[i]
        info = r["info"]
        rmsd = info.get("rmsd_heavy_A")
        rmsd_s = f"{rmsd:.3f} A" if rmsd is not None else "-"
        mc = info.get("min_contact") or {}
        mcs = f"{mc.get('distance_A', '-')}"
        rd = info.get("reactive_distances") or []
        rds = ", ".join(f"{x['atoms'][0]}-{x['atoms'][1]}:{x['distance_A']}"
                        for x in rd) or "-"
        if len(rds) > 25:
            rds = rds[:24] + "~"
        msgs = (r["problems"] + r["warnings"])
        msg = msgs[0] if msgs else ""
        print(f"{rank:<5}{r['index']:<6}{r['status']:<7}{rmsd_s:<13}"
              f"{mcs:<8}{info.get('n_clash_pairs', 0):<7}{rds:<26}{msg[:60]}")

    counts = {}
    for r in results:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    print()
    print("summary: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))

    n_fail = counts.get("FAIL", 0)
    if n_fail:
        print(f"  {n_fail} candidate(s) FAILed a physical-plausibility check and "
              "should not be sent to DFT.")
    n_pass = counts.get("PASS", 0)
    if n_pass == 0:
        print("  no clean candidate: consider regenerating with more samples "
              "(the ML route gets markedly better with ~10 samples per reaction).")

    if args.write_pass:
        keep = [frames[i] for i in order if results[i]["status"] == "PASS"]
        if keep:
            with open(args.write_pass, "w", encoding="utf-8") as fh:
                for k, atoms in enumerate(keep):
                    fh.write(f"{len(atoms)}\npass-{k}\n")
                    for sym, (x, y, z) in atoms:
                        fh.write(f"{sym:<2} {x:16.6f} {y:16.6f} {z:16.6f}\n")
            print(f"  wrote {len(keep)} passing candidate(s) to {args.write_pass}")

    if args.json:
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump({"candidates": results,
                       "ranking": [results[i]["index"] for i in order]},
                      fh, indent=2)
        print(f"  report -> {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

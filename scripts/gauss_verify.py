#!/usr/bin/env python3
"""Verify a Gaussian 16 frequency / IRC / single-point output.

This is the Gaussian counterpart of ts_verify.py. It checks the things that
can be checked from the text output:

  * normal vs error termination
  * number and size of imaginary frequencies
  * geometry optimization convergence
  * SCF energy, thermal free energy, <S**2>
  * IRC completion (when --mode irc)

It cannot decide whether the imaginary mode is the reaction coordinate you
intended. For that, run `freq=hpmodes` and inspect the mode in GaussView, or
compare the key bond lengths at the two IRC endpoints.

Examples
--------
python gauss_verify.py --out 02_ts_optfreq.log
python gauss_verify.py --out 03_irc.log --mode irc
python gauss_verify.py --out 04_sp.log --mode sp --json

Exit code is 0 for SADDLE-1 (freq mode), a normally terminated SP, or a
completed IRC; 1 otherwise. This makes it usable as a shell gate.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path


def last(pattern, text, cast=float, default=None):
    vals = re.findall(pattern, text)
    if not vals:
        return default
    try:
        return cast(vals[-1])
    except (TypeError, ValueError):
        return default


def parse(path: Path):
    text = path.read_text(errors="replace")
    normal = "Normal termination" in text
    error = "Error termination" in text

    freqs = []
    for line in text.splitlines():
        if "Frequencies --" in line:
            tail = line.split("Frequencies --", 1)[1]
            freqs.extend(float(v) for v in re.findall(r"-?\d+\.\d+", tail))
    negative = [v for v in freqs if v < 0]
    significant = [v for v in negative if v < -5.0]

    e_scf = last(r"SCF Done:\s+E\([^)]+\)\s+=\s+(-?\d+\.\d+)", text)
    g = last(r"Sum of electronic and thermal Free Energies=\s+(-?\d+\.\d+)", text)
    zpe = last(r"Zero-point correction=\s+(-?\d+\.\d+)", text)
    s2 = last(r"<S\*\*2>=\s*(-?\d+\.\d+)", text)
    stationary = ("Stationary point found" in text
                  or "Optimization completed" in text)
    irc_done = ("IRC calculation complete" in text
                or "IRC has converged" in text
                or "Calculation of IRC complete" in text)
    irc_max = "Maximum number of IRC steps" in text

    return {
        "file": str(path),
        "normal": normal,
        "error": error,
        "n_freq": len(freqs),
        "negative": negative,
        "significant_negative": significant,
        "e_scf": e_scf,
        "g": g,
        "zpe": zpe,
        "s2": s2,
        "stationary": stationary,
        "irc_done": irc_done,
        "irc_max_steps": irc_max,
    }


def classify(r, mode):
    if not r["normal"] or r["error"]:
        return "FAIL"
    if mode == "sp":
        return "SP-OK"
    if mode == "irc":
        if r["irc_done"]:
            return "IRC-OK"
        if r["irc_max_steps"]:
            return "IRC-MAX-STEPS"
        return "IRC-INCOMPLETE"
    # frequency mode
    if r["n_freq"] == 0:
        return "NO-FREQ"
    sig = r["significant_negative"]
    if len(sig) == 0:
        return "MINIMUM"
    if len(sig) == 1:
        return "SADDLE-1"
    if 2 <= len(sig) <= 3 and all(-25.0 <= v <= -5.0 for v in sig):
        return "FLAT-VALLEY"
    return f"ORDER-{len(sig)}-SADDLE"


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, help="Gaussian .log file")
    ap.add_argument("--mode", choices=["freq", "irc", "sp"], default="freq")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    path = Path(args.out)
    if not path.exists():
        sys.exit(f"file not found: {path}")
    r = parse(path)
    verdict = classify(r, args.mode)
    r["verdict"] = verdict
    r["mode"] = args.mode

    if args.json:
        print(json.dumps(r, indent=2, ensure_ascii=False))
    else:
        print(f"file: {r['file']}")
        print(f"termination: {'NORMAL' if r['normal'] else 'NOT NORMAL'}"
              + ("  [Error termination]" if r["error"] else ""))
        print(f"stationary point: {r['stationary']}")
        print(f"frequencies parsed: {r['n_freq']}")
        print(f"negative frequencies: {r['negative']}")
        print(f"E(SCF): {r['e_scf']}")
        print(f"G: {r['g']}   ZPE: {r['zpe']}   <S**2>: {r['s2']}")
        if args.mode == "irc":
            print(f"IRC complete: {r['irc_done']}   max steps: {r['irc_max_steps']}")
        print(f"VERDICT: {verdict}")
        if verdict in ("ORDER-2-SADDLE", "ORDER-3-SADDLE"):
            print("  -> multiple significant imaginary modes; do not report as a "
                  "first-order TS without a decision (see references/ts-protocol.md).")
        if verdict == "FLAT-VALLEY":
            print("  -> shallow imaginary modes; likely a flat transition valley. "
                  "Use the last geometry's Freq as a reference barrier and say so.")
        if verdict == "SADDLE-1":
            print("  -> still check the mode identity and the IRC endpoints; "
                  "a saddle is not necessarily the reaction you intended.")

    ok = verdict == "SADDLE-1" or (args.mode == "sp" and verdict == "SP-OK") \
        or (args.mode == "irc" and verdict == "IRC-OK")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

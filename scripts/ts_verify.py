#!/usr/bin/env python3
"""TS verification for ORCA 6.x output.

Turns an ORCA .out into a verdict on whether the stationary point is a genuine
first-order saddle for the intended reaction coordinate.

Checks implemented (rationale in SKILL.md):
  1. normal termination / optimization convergence flags
  2. count of *significant* imaginary frequencies (default: < -20 cm**-1)
  3. for each imaginary mode, the atom pairs whose interatomic distance changes
     most -> tells you which bonds are forming/breaking
  4. whether the intended reactive pair is among those pairs (this is the
     programmatic version of the "filter out mere group rotations" step)
  5. electronic energy / Gibbs energy extraction
  6. optional IRC output check

Only the standard library is used.

Usage:
  python ts_verify.py --out job.out
  python ts_verify.py --out job.out --reactive 21,56 --reactive 20,71
  python ts_verify.py --out irc.out --irc
  python ts_verify.py --out job.out --json
"""

import argparse
import json
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ts_common import (ATOMIC_MASSES, COVALENT_RADII, HARTREE_KCAL)  # noqa: E402

# --- constants -------------------------------------------------------------

# Significant imaginary frequency threshold (cm**-1). The user's field-tested
# rule: shallow -15..-25 values are broad flat-valley artefacts; -600..-700
# means a genuine saddle but possibly far from the unconstrained stationary
# point. Numeric noise never exceeds a few cm**-1.
SIG_IMAG = -20.0        # counted towards the saddle order
NOISE_FLOOR = -5.0      # above this an "imaginary" value is numerical noise
FLAT_VALLEY_LO = -25.0  # field-verified broad flat transition valley band
FLAT_VALLEY_HI = -15.0


# --- parsing helpers -------------------------------------------------------

def _last(pattern, text, flags=0, cast=float, default=None):
    """Return the LAST match (critical: ORCA repeats energy lines every cycle,
    and the first match belongs to the initial geometry)."""
    hits = re.findall(pattern, text, flags)
    if not hits:
        return default
    try:
        return cast(hits[-1])
    except (TypeError, ValueError):
        return default


def parse_geometry(text):
    """Last 'CARTESIAN COORDINATES (ANGSTROEM)' block -> [(sym, (x,y,z)), ...]."""
    blocks = re.findall(
        r"CARTESIAN COORDINATES \(ANGSTROEM\)\s*-*\s*\n(.*?)\n\s*\n",
        text, re.S)
    if not blocks:
        return []
    atoms = []
    for line in blocks[-1].splitlines():
        parts = line.split()
        if len(parts) >= 4:
            try:
                atoms.append((parts[0], (float(parts[1]), float(parts[2]),
                                         float(parts[3]))))
            except ValueError:
                continue
    return atoms


def parse_frequencies(text):
    """-> [(mode_idx, value_cm-1, is_marked_imaginary), ...] (harmonic block)."""
    out = []
    for m in re.finditer(
            r"^\s*(\d+):\s+(-?\d+\.\d+)\s*cm\*\*-1(\s+\*\*\*imaginary mode\*\*\*)?",
            text, re.M):
        idx = int(m.group(1))
        val = float(m.group(2))
        marked = bool(m.group(3))
        out.append((idx, val, marked))
    return out


def parse_normal_modes(text):
    """Parse the NORMAL MODES matrix.

    Layout (verified against ORCA 6.1.1 output):
      header row: a run of consecutive mode indices, e.g. ' 0 1 2 3 4 5'
      data rows : row_index followed by one value per column (mode)
    row_index r is the Cartesian component: atom = r // 3, comp = r % 3.

    The matrix is printed in blocks of (usually) 6 modes, and the LAST block is
    often short (e.g. 3 columns for the 36 remaining modes of a 39-mode job), so
    the header run must be accepted at any length >= 1 -- requiring a fixed
    width silently dropped the trailing block and lost 6 modes. The section is
    bounded by the following IR SPECTRUM heading rather than by a fixed
    character window, which used to truncate large molecules.

    Returns {mode_index: [(vx,vy,vz) per atom]}.
    """
    start = text.find("NORMAL MODES")
    if start < 0:
        return {}
    end = text.find("IR SPECTRUM", start)
    body = text[start:end] if end > 0 else text[start:]

    blocks = []          # [(col_modes, {row_index: [values]})]
    cur = None
    for line in body.splitlines():
        s = line.strip()
        if not s:
            continue
        toks = s.split()
        if all(re.fullmatch(r"\d+", t) for t in toks):
            idxs = [int(t) for t in toks]
            if idxs == list(range(idxs[0], idxs[0] + len(idxs))):
                cur = (idxs, {})
                blocks.append(cur)
                continue
        if cur is not None and re.fullmatch(r"\d+", toks[0]):
            try:
                vals = [float(v) for v in toks[1:]]
            except ValueError:
                continue
            if len(vals) == len(cur[0]):
                cur[1][int(toks[0])] = vals

    modes = {}
    for col_modes, data in blocks:
        if not data:
            continue
        n_rows = max(data) + 1
        if n_rows % 3:
            continue
        n_atoms = n_rows // 3
        for c, m_idx in enumerate(col_modes):
            vecs = []
            for a in range(n_atoms):
                try:
                    v = (data[3 * a][c], data[3 * a + 1][c], data[3 * a + 2][c])
                except (KeyError, IndexError):
                    vecs = []
                    break
                vecs.append(v)
            if vecs:
                modes[m_idx] = vecs
    return modes


def bond_stretch_analysis(mode_vecs, geometry, top_n=8, cutoff=3.2):
    """Which interatomic distances change most in this mode?

    Only pairs within `cutoff` (A) are considered - far-apart pairs project
    huge but meaningless "distance changes" onto the mode vector, so a
    chemically relevant window is required. 3.2 A covers a forming bond;
    anything longer is not a bond and not a reaction coordinate.

    Returns [(i, j, delta_d_signed, fraction, current_d, is_bonded), ...]
    sorted by |delta_d| desc. Atom indices are 0-based.

    ⚠️ THE SIGN IS ARBITRARY. A normal-mode eigenvector is defined only up to
    an overall sign, and ORCA may print either one, so `delta_d_signed` cannot
    tell you whether a bond is forming or breaking. Read the *magnitudes*: a
    pair with fraction 1.00 changes its distance an order of magnitude more
    than one at 0.10, and that ranking is what identifies the reaction
    coordinate. To assign forming vs breaking, compare against the relaxed
    scan, the IRC endpoints, or the reactant/product geometries - never
    against the mode's own sign. An earlier version of this script labelled
    the sign as forming/breaking, which was wrong on both counts (the
    threshold was inverted AND the quantity is gauge-dependent).

    The ORCA normal-mode vectors are mass-weighted (M_ii = 1/sqrt(m_i)) and
    normalised; they are un-weighted here with sqrt(m) so the projected
    displacement has a physical length interpretation.
    """
    n = len(mode_vecs)
    if n == 0 or len(geometry) != n:
        return []
    # un-weight
    cart = []
    for (sym, _), v in zip(geometry, mode_vecs):
        m = ATOMIC_MASSES.get(sym, 12.0)
        s = math.sqrt(m)
        cart.append((v[0] * s, v[1] * s, v[2] * s))

    out = []
    for i in range(n):
        for j in range(i + 1, n):
            si, ri = geometry[i]
            sj, rj = geometry[j]
            d = math.dist(ri, rj)
            if d < 1e-6 or d > cutoff:
                continue
            u = [(rj[k] - ri[k]) / d for k in range(3)]
            # change in bond length ~ u . (v_j - v_i)
            dd = sum(u[k] * (cart[j][k] - cart[i][k]) for k in range(3))
            rbond = COVALENT_RADII.get(si, 0.8) + COVALENT_RADII.get(sj, 0.8)
            out.append((i, j, dd, d, d < 1.3 * rbond))
    if not out:
        return []
    out.sort(key=lambda t: -abs(t[2]))
    scale = abs(out[0][2]) or 1.0
    return [(i, j, dd, abs(dd) / scale, d, bonded)
            for i, j, dd, d, bonded in out[:top_n]]


def atomic_displacement_ranking(mode_vecs, geometry, top_n=8):
    """Per-atom displacement magnitude as a fraction of the largest one."""
    if not mode_vecs:
        return []
    mags = []
    for (sym, _), v in zip(geometry, mode_vecs):
        m = ATOMIC_MASSES.get(sym, 12.0)
        s = math.sqrt(m)
        mags.append((math.sqrt(sum((c * s) ** 2 for c in v)), sym))
    if not mags:
        return []
    mx = max(m for m, _ in mags) or 1.0
    order = sorted(range(len(mags)), key=lambda i: -mags[i][0])
    return [(i, mags[i][1], mags[i][0] / mx) for i in order[:top_n]]


# --- main verification -----------------------------------------------------

def verify(path, reactive_pairs=None, irc=False):
    text = open(path, encoding="utf-8", errors="ignore").read()
    res = {"file": path, "checks": [], "warnings": [], "errors": []}

    def note(ok, label, detail=""):
        res["checks"].append({"ok": ok, "label": label, "detail": detail})

    # 1. completion
    terminated = "ORCA TERMINATED NORMALLY" in text
    note(terminated, "normal termination",
         "ORCA TERMINATED NORMALLY" if terminated else "termination marker NOT found")
    if not terminated:
        res["errors"].append("ORCA did not terminate normally - job was killed "
                             "or crashed; do not interpret this output")

    # 2. energies
    e_elec = _last(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", text)
    e_smd = _last(r"Total Energy after SMD CDS correction\s*=\s*(-?\d+\.\d+)", text)
    g = _last(r"Final Gibbs free energy[\s.:]*(-?\d+\.\d+)", text)
    res["energy"] = {
        "final_single_point_Eh": e_elec,
        "smd_corrected_Eh": e_smd,
        "gibbs_Eh": g,
    }

    # 3. convergence flags
    opt_converged = "THE OPTIMIZATION HAS CONVERGED" in text
    opt_requested = bool(re.search(r"^\s*!.*\bOpt\b", text, re.M)) or \
        bool(re.search(r"\bOptTS\b", text))
    if opt_requested:
        note(opt_converged, "optimization converged",
             "converged" if opt_converged
             else "NOT converged - the reported geometry is not a stationary point")
        if not opt_converged:
            res["warnings"].append(
                "Optimization did not converge. Per the project rule, an "
                "unconverged OptTS endpoint must pass a connectivity check "
                "before any interpretation, and negative results must be "
                "reported as such.")

    # 4. IRC mode
    if irc:
        conv = "THE IRC HAS CONVERGED" in text
        note(conv, "IRC converged",
             "converged" if conv else "IRC did NOT converge")
        res["irc"] = {
            "converged": conv,
            "max_iterations_reached": "MAXIMUM NUMBER OF ITERATIONS REACHED - STOPPING IRC RUN" in text,
        }

    # 5. frequencies
    freqs = parse_frequencies(text)
    if not freqs:
        res["frequencies"] = None
        note(False, "vibrational analysis present",
             "no VIBRATIONAL FREQUENCIES block found - run Freq to verify "
             "the stationary point")
        res["verdict"] = "NO-FREQ"
        return res

    # drop the 6 (or 5 for linear) zero modes; keep everything else
    zero_modes = [f for f in freqs if abs(f[1]) <= 1.0]
    n_zero = len(zero_modes)
    real = [f for f in freqs if abs(f[1]) > 1.0]
    positive = [f for f in real if f[1] > 0]
    # mutually exclusive bands over the imaginary values
    sig_imag = [f for f in real if f[1] <= SIG_IMAG]                 # counted
    shallow_imag = [f for f in real if SIG_IMAG < f[1] <= NOISE_FLOOR]  # not counted
    # the field-verified "broad flat valley" signature band (-25..-15 cm-1);
    # annotated separately because it spans the counting threshold
    flat_valley = [f for f in real
                   if FLAT_VALLEY_LO < f[1] < FLAT_VALLEY_HI]

    res["frequencies"] = {
        "n_modes_total": len(freqs),
        "n_zero_modes": n_zero,
        "imaginary_significant": [[i, v] for i, v, _ in sig_imag],
        "imaginary_shallow": [[i, v] for i, v, _ in shallow_imag],
        "flat_valley_band": [[i, v] for i, v, _ in flat_valley],
        "lowest_real": round(min((v for _, v, _ in positive), default=0.0), 2),
    }

    expected_zero = 5 if re.search(r"linear molecule", text, re.I) else 6
    if n_zero not in (expected_zero, expected_zero + 1):
        res["warnings"].append(
            f"{n_zero} zero modes found (expected {expected_zero}); check "
            "whether the geometry is a proper stationary point.")

    # 6. imaginary mode character
    geometry = parse_geometry(text)
    modes = parse_normal_modes(text)
    # Invariant: ORCA prints one Cartesian component row per 3N, so a complete
    # normal-mode matrix must yield exactly 3N mode vectors. Anything less means
    # the parser lost a block, and any mode looked up from it may be missing.
    if geometry:
        expected = 3 * len(geometry)
        res["normal_modes"] = {"parsed": len(modes), "expected_3N": expected}
        if len(modes) != expected:
            res["warnings"].append(
                f"normal-mode matrix incomplete: parsed {len(modes)} of "
                f"{expected} (3N) mode vectors - do NOT trust mode analysis "
                "for modes that are missing; the parser may have dropped a block")
    res["mode_analysis"] = []
    reactive_atoms = set()
    for pair in (reactive_pairs or []):
        reactive_atoms.update(pair)

    for idx, val, _ in sig_imag:
        entry = {"mode_index": idx, "frequency_cm-1": val}
        vecs = modes.get(idx)
        if vecs and geometry and len(vecs) == len(geometry):
            entry["bond_stretch"] = [
                {"atoms": [i + 1, j + 1], "delta_d_rel": round(frac, 3),
                 "current_distance_A": round(d, 3),
                 "bonded_at_start": bonded,
                 "signed_projection_arbitrary": round(dd, 4)}
                for i, j, dd, frac, d, bonded in bond_stretch_analysis(vecs, geometry)
            ]
            entry["bond_stretch_note"] = (
                "ranked by how much each interatomic distance changes in this "
                "mode; the SIGN is arbitrary (eigenvector gauge) and does NOT "
                "mean forming vs breaking")
            entry["top_atoms"] = [
                {"atom": i + 1, "element": sym, "displacement_rel": round(frac, 3)}
                for i, sym, frac in atomic_displacement_ranking(vecs, geometry, 6)
            ]
            # does the intended reactive pair participate?
            pairs_seen = {frozenset((b["atoms"][0], b["atoms"][1]))
                          for b in entry["bond_stretch"][:5]}
            if reactive_atoms:
                hit = any(frozenset(p) in pairs_seen for p in reactive_pairs)
                involved = bool(reactive_atoms & {a["atom"] for a in entry["top_atoms"]})
                entry["matches_reactive_pair"] = hit
                entry["involves_reactive_atoms"] = involved
                if not hit and not involved:
                    res["warnings"].append(
                        f"Imaginary mode {idx} ({val} cm-1) does NOT involve the "
                        "declared reactive atoms - suspect a mere group rotation "
                        "or an unintended reaction coordinate")
        else:
            entry["bond_stretch"] = None
            entry["note"] = ("normal-mode matrix unavailable or size mismatch; "
                             "inspect the NORMAL MODES block manually")
        res["mode_analysis"].append(entry)

    # 7. verdict
    n_sig = len(sig_imag)
    n_flat = len(flat_valley)
    if not terminated:
        verdict = "FAIL"
    elif n_sig == 1:
        verdict = "SADDLE-1"
        if irc and res.get("irc", {}).get("converged") is False:
            verdict = "SADDLE-1-IRC-UNCONVERGED"
    elif n_sig == 0:
        verdict = "MINIMUM"
        res["warnings"].append(
            "No significant imaginary frequency: this is a local MINIMUM, not a "
            "transition state. For a barrier claim you need a genuine TS, or an "
            "explicit statement that a constrained reference point is being used.")
    elif 2 <= n_sig <= 3 and n_flat >= n_sig:
        verdict = "FLAT-VALLEY"
        res["warnings"].append(
            "2-3 shallow imaginary modes, all inside the -15..-25 cm-1 band: the "
            "signature of a broad flat transition valley. Repeated OptTS restarts "
            "are usually wasted effort - take the endpoint geometry, use its Freq "
            "energy as a reference barrier, and say so explicitly in the writeup.")
    else:
        verdict = f"ORDER-{n_sig}-SADDLE"
        res["warnings"].append(
            f"{n_sig} significant imaginary modes: not a first-order saddle. "
            "Either a higher-order saddle, or the structure sits on a "
            "symmetry-degenerate/flat surface. Inspect each mode's character "
            "below before using this point for anything.")
    res["verdict"] = verdict
    return res


def fmt(res, show_modes=True):
    lines = []
    v = res["verdict"]
    lines.append(f"=== {res['file']} ===")
    lines.append(f"VERDICT: {v}")
    e = res.get("energy", {})
    if e.get("final_single_point_Eh") is not None:
        lines.append(f"  E(SP)   = {e['final_single_point_Eh']:.6f} Eh")
    if e.get("smd_corrected_Eh") is not None:
        lines.append(f"  E(SMD)  = {e['smd_corrected_Eh']:.6f} Eh")
    if e.get("gibbs_Eh") is not None:
        lines.append(f"  G       = {e['gibbs_Eh']:.6f} Eh")
    f = res.get("frequencies")
    if f:
        lines.append(f"  modes={f['n_modes_total']}  zero={f['n_zero_modes']}  "
                     f"lowest positive={f['lowest_real']} cm-1")
        if f["imaginary_significant"]:
            s = ", ".join(f"{i}:{v} " for i, v in f["imaginary_significant"])
            lines.append(f"  imaginary (significant): {s}cm-1")
        if f["imaginary_shallow"]:
            s = ", ".join(f"{i}:{v} " for i, v in f["imaginary_shallow"])
            lines.append(f"  imaginary (shallow, not counted): {s}cm-1")
        if f["flat_valley_band"]:
            s = ", ".join(f"{i}:{v} " for i, v in f["flat_valley_band"])
            lines.append(f"  in flat-valley band (-15..-25): {s}cm-1")
    if res.get("irc"):
        i = res["irc"]
        lines.append(f"  IRC: converged={i['converged']}"
                     + ("  (hit iteration limit)" if i["max_iterations_reached"] else ""))
    if show_modes and res.get("mode_analysis"):
        for m in res["mode_analysis"]:
            lines.append(f"  --- imaginary mode {m['mode_index']} "
                         f"({m['frequency_cm-1']} cm-1) ---")
            if m.get("bond_stretch"):
                lines.append("      pairs whose distance changes most "
                             "(rank by magnitude; sign is meaningless):")
                for b in m["bond_stretch"][:5]:
                    a1, a2 = b["atoms"]
                    tag = "bond" if b["bonded_at_start"] else "nonbond"
                    lines.append(f"      {a1:>4}-{a2:<4} rel={b['delta_d_rel']:.3f}  "
                                 f"d={b['current_distance_A']} A  [{tag}]")
            if m.get("top_atoms"):
                s = ", ".join(f"{t['element']}{t['atom']}({t['displacement_rel']:.2f})"
                              for t in m["top_atoms"])
                lines.append(f"      moving atoms: {s}")
            if "matches_reactive_pair" in m:
                lines.append(f"      matches declared reactive pair: "
                             f"{m['matches_reactive_pair']}")
    for c in res["checks"]:
        lines.append(f"  [{'ok ' if c['ok'] else 'FAIL'}] {c['label']}"
                     + (f" - {c['detail']}" if c["detail"] else ""))
    for w in res["warnings"]:
        lines.append(f"  ! {w}")
    for er in res["errors"]:
        lines.append(f"  X {er}")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", nargs="+", required=True, help="ORCA .out file(s)")
    ap.add_argument("--reactive", action="append", default=None,
                    metavar="I,J",
                    help="declared reactive atom pair, 1-based, repeatable "
                         "(the bond(s) forming/breaking)")
    ap.add_argument("--irc", action="store_true",
                    help="also expect/report IRC convergence")
    ap.add_argument("--json", action="store_true", help="emit JSON")
    args = ap.parse_args()

    pairs = []
    for spec in (args.reactive or []):
        parts = [int(x) for x in re.split(r"[,\s]+", spec.strip()) if x]
        if len(parts) != 2:
            sys.exit(f"--reactive needs an I,J pair, got {spec!r}")
        pairs.append(tuple(parts))

    results = [verify(p, pairs, args.irc) for p in args.out]
    if args.json:
        print(json.dumps(results, indent=2, default=str))
    else:
        for r in results:
            print(fmt(r))
            print()

    # batch summary: rank by Gibbs if several files
    ranked = [r for r in results if r.get("energy", {}).get("gibbs_Eh") is not None]
    if len(ranked) > 1:
        ranked.sort(key=lambda r: r["energy"]["gibbs_Eh"])
        ref = ranked[0]["energy"]["gibbs_Eh"]
        print("=== ranking by G (relative to the lowest) ===")
        for r in ranked:
            dg = (r["energy"]["gibbs_Eh"] - ref) * HARTREE_KCAL
            print(f"  dG={dg:>8.2f} kcal/mol  {r['verdict']:<12} {r['file']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

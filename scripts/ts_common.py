#!/usr/bin/env python3
"""Shared helpers for the TS computation suite.

Element tables, xyz I/O, geometry checks (clash ratio, bonding), and a
Kabsch-aligned RMSD. Standard library only; NumPy is used when present.
"""

import math
import re

HARTREE_KCAL = 627.509474

COVALENT_RADII = {
    "H": 0.31, "He": 0.28, "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76,
    "N": 0.71, "O": 0.66, "F": 0.57, "Ne": 0.58, "Na": 1.66, "Mg": 1.41,
    "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Ar": 1.06,
    "K": 2.03, "Ca": 1.76, "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39,
    "Mn": 1.39, "Fe": 1.32, "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22,
    "Ga": 1.22, "Ge": 1.20, "As": 1.19, "Se": 1.20, "Br": 1.20, "Kr": 1.16,
    "Rb": 2.20, "Sr": 1.95, "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54,
    "Tc": 1.47, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Cd": 1.44,
    "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39, "Xe": 1.40,
    "Cs": 2.44, "Ba": 2.15, "La": 2.07, "Ce": 2.04, "Lu": 1.87, "Hf": 1.75,
    "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44, "Ir": 1.41, "Pt": 1.36,
    "Au": 1.36, "Hg": 1.32, "Tl": 1.45, "Pb": 1.46, "Bi": 1.48, "Po": 1.40,
    "At": 1.50, "Rn": 1.50,
}

ATOMIC_MASSES = {
    "H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81,
    "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180,
    "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974,
    "S": 32.06, "Cl": 35.45, "Ar": 39.948, "K": 39.098, "Ca": 40.078,
    "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996, "Mn": 54.938,
    "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38,
    "Ga": 69.723, "Ge": 72.630, "As": 74.922, "Se": 78.971, "Br": 79.904,
    "Kr": 83.798, "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224,
    "Nb": 92.906, "Mo": 95.95, "Tc": 98.0, "Ru": 101.07, "Rh": 102.91,
    "Pd": 106.42, "Ag": 107.87, "Cd": 112.41, "In": 114.82, "Sn": 118.71,
    "Sb": 121.76, "Te": 127.60, "I": 126.90, "Xe": 131.29, "Cs": 132.91,
    "Ba": 137.33, "La": 138.91, "Ce": 140.12, "Lu": 174.97, "Hf": 178.49,
    "Ta": 180.95, "W": 183.84, "Re": 186.21, "Os": 190.23, "Ir": 192.22,
    "Pt": 195.08, "Au": 196.97, "Hg": 200.59, "Tl": 204.38, "Pb": 207.2,
    "Bi": 208.98, "Po": 209.0, "At": 210.0, "Rn": 222.0,
}

# Bonding / clash thresholds
BOND_FACTOR = 1.3        # bonded if d < BOND_FACTOR * (rcov_i + rcov_j)
CLASH_FACTOR = 0.75      # clash if d < CLASH_FACTOR * (rcov_i + rcov_j)
# (the clash factor is the criterion defined in the UniTS paper)
HARD_OVERLAP = 1.00      # any pair closer than this is physically impossible


# --- xyz I/O ---------------------------------------------------------------

def read_xyz_frames(path, last_only=False):
    """Read one or more xyz frames -> [ [(sym, (x,y,z)), ...], ... ].

    Frames are concatenated frames in a single file (as UniTS saves guesses).

    A frame whose coordinate block is shorter than its declared atom count is a
    hard error, not something to pad or truncate: silently returning a
    "molecule with a few atoms missing" would let a partially written file pass
    the physical-plausibility screen and be sent to DFT.

    Headerless coordinate dumps (element + 3 floats per line, no atom-count
    line) are accepted as a single frame - the project keeps several of these in
    `00-geom/`, exported by hand from the original Gaussian workflow. The
    discrimination is unambiguous because a standard first line is a bare
    integer while a coordinate line starts with an element symbol.
    """
    lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
    frames = []
    i = 0
    n = len(lines)
    while i < n:
        while i < n and not lines[i].strip():
            i += 1
        if i >= n:
            break
        try:
            natoms = int(lines[i].strip())
        except ValueError:
            if not frames and i == 0:
                return _read_headerless(path, lines)
            raise ValueError(
                f"{path}: line {i + 1} should be an atom count, "
                f"found {lines[i].strip()!r}")
        if natoms <= 0:
            raise ValueError(f"{path}: line {i + 1} declares {natoms} atoms")
        if i + 2 + natoms > n:
            raise ValueError(
                f"{path}: frame starting at line {i + 1} declares {natoms} atoms "
                f"but only {n - (i + 2)} coordinate line(s) follow - the file is "
                "truncated; refusing to analyse a partial structure")
        atoms = []
        for k in range(i + 2, i + 2 + natoms):
            p = lines[k].split()
            if len(p) < 4:
                raise ValueError(
                    f"{path}: line {k + 1} has {len(p)} field(s), expected "
                    "element + 3 coordinates")
            try:
                atoms.append((p[0], (float(p[1]), float(p[2]), float(p[3]))))
            except ValueError:
                raise ValueError(
                    f"{path}: line {k + 1} has non-numeric coordinates: "
                    f"{lines[k].strip()!r}")
        frames.append(atoms)
        i = i + 2 + natoms
    if last_only and frames:
        return [frames[-1]]
    return frames


def _read_headerless(path, lines):
    """Fallback for coordinate dumps with no atom-count/comment header."""
    atoms = []
    for k, line in enumerate(lines):
        s = line.strip()
        if not s:
            continue
        p = s.split()
        if len(p) < 4:
            raise ValueError(
                f"{path}: not a valid xyz file - line 1 is neither an atom count "
                f"nor a coordinate line, and line {k + 1} has only {len(p)} "
                f"field(s)")
        try:
            atoms.append((p[0], (float(p[1]), float(p[2]), float(p[3]))))
        except ValueError:
            raise ValueError(
                f"{path}: not a valid xyz file (no atom-count header) and line "
                f"{k + 1} is not element+3-coordinates: {s!r}")
    if not atoms:
        raise ValueError(f"{path}: no coordinates found")
    return [atoms]


def write_xyz(path, atoms, comment=""):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"{len(atoms)}\n{comment}\n")
        for sym, (x, y, z) in atoms:
            fh.write(f"{sym:<2} {x:16.6f} {y:16.6f} {z:16.6f}\n")


def xyz_repr(atoms):
    return "\n".join(f"{s:<2} {x:16.6f} {y:16.6f} {z:16.6f}"
                     for s, (x, y, z) in atoms)


# --- geometry checks -------------------------------------------------------

def pair_distances(atoms):
    """[(i, j, d, rsum), ...] for every unique pair, 0-based indices."""
    out = []
    for i in range(len(atoms)):
        si, ri = atoms[i]
        for j in range(i + 1, len(atoms)):
            sj, rj = atoms[j]
            d = math.dist(ri, rj)
            out.append((i, j, d,
                        COVALENT_RADII.get(si, 0.8) + COVALENT_RADII.get(sj, 0.8)))
    return out


def min_contact(atoms):
    """(i, j, d, bonded) of the closest non-trivial pair, or None."""
    best = None
    for i, j, d, rsum in pair_distances(atoms):
        if best is None or d < best[2]:
            best = (i, j, d, d < BOND_FACTOR * rsum)
    return best


def clash_report(atoms):
    """Structural-clash analysis using the clash criterion from the UniTS paper.

    Returns dict with the clash ratio for this single structure (0.0 or 1.0),
    the offending pairs, and the hard-overlap pairs.
    """
    clashes = []
    overlaps = []
    for i, j, d, rsum in pair_distances(atoms):
        if d < CLASH_FACTOR * rsum:
            clashes.append((i, j, d, CLASH_FACTOR * rsum,
                            d < BOND_FACTOR * rsum))
        if d < HARD_OVERLAP:
            overlaps.append((i, j, d))
    return {
        "has_clash": bool(clashes),
        "n_clash_pairs": len(clashes),
        "clash_pairs": clashes,
        "n_hard_overlaps": len(overlaps),
        "hard_overlaps": overlaps,
    }


def bond_graph(atoms):
    """Set of frozenset({i,j}) deemed bonded by the covalent-radius criterion."""
    return {frozenset((i, j)) for i, j, d, rsum in pair_distances(atoms)
            if d < BOND_FACTOR * rsum}


def reactive_distances(atoms, pairs):
    """[(i, j, d), ...] for the declared reactive pairs (1-based in, 0-based out)."""
    out = []
    for a, b in pairs:
        if 1 <= a <= len(atoms) and 1 <= b <= len(atoms):
            out.append((a, b, math.dist(atoms[a - 1][1], atoms[b - 1][1])))
    return out


# --- Kabsch RMSD -----------------------------------------------------------

def kabsch_rmsd(ref_atoms, cand_atoms, heavy_only=False):
    """RMSD after optimal rigid-body alignment (rotation + translation).

    Both structures must have the same atoms in the same order (this is the
    case for UniTS guesses: same graph, same ordering).
    """
    idx = list(range(len(ref_atoms)))
    if heavy_only:
        idx = [k for k in idx if ref_atoms[k][0] != "H"]
    if not idx:
        return None
    A = [ref_atoms[k][1] for k in idx]
    B = [cand_atoms[k][1] for k in idx]
    return _kabsch(A, B)


def _kabsch(A, B):
    n = len(A)
    ca = [sum(p[k] for p in A) / n for k in range(3)]
    cb = [sum(p[k] for p in B) / n for k in range(3)]
    A2 = [[p[k] - ca[k] for k in range(3)] for p in A]
    B2 = [[p[k] - cb[k] for k in range(3)] for p in B]

    try:
        import numpy as np
        P = np.array(A2)
        Q = np.array(B2)
        H = P.T @ Q
        U, _, Vt = np.linalg.svd(H)
        d = np.sign(np.linalg.det(Vt.T @ U.T))
        D = np.diag([1.0, 1.0, d])
        R = Vt.T @ D @ U.T
        diff = (R @ P.T).T - Q
        return float(np.sqrt((diff ** 2).sum() / n))
    except ImportError:
        pass

    # Pure-python fallback: Horn's quaternion method.
    S = [[sum(A2[i][a] * B2[i][b] for i in range(n)) for b in range(3)]
         for a in range(3)]
    K = [
        [S[0][0] + S[1][1] + S[2][2], S[1][2] - S[2][1], S[2][0] - S[0][2],
         S[0][1] - S[1][0]],
        [S[1][2] - S[2][1], S[0][0] - S[1][1] - S[2][2], S[0][1] + S[1][0],
         S[2][0] + S[0][2]],
        [S[2][0] - S[0][2], S[0][1] + S[1][0], -S[0][0] + S[1][1] - S[2][2],
         S[1][2] + S[2][1]],
        [S[0][1] - S[1][0], S[2][0] + S[0][2], S[1][2] + S[2][1],
         -S[0][0] - S[1][1] + S[2][2]],
    ]
    q = _jacobi_largest_eigvec(K)
    R = _quat_to_matrix(q)
    ss = 0.0
    for i in range(n):
        rp = [sum(R[a][b] * A2[i][b] for b in range(3)) for a in range(3)]
        ss += sum((rp[k] - B2[i][k]) ** 2 for k in range(3))
    return math.sqrt(ss / n)


def _jacobi_largest_eigvec(M, iters=100):
    """Largest-eigenvalue eigenvector of a real symmetric 4x4 matrix."""
    n = len(M)
    A = [row[:] for row in M]
    V = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(iters):
        p, q, mx = 0, 1, 0.0
        for i in range(n):
            for j in range(i + 1, n):
                if abs(A[i][j]) > mx:
                    mx, p, q = abs(A[i][j]), i, j
        if mx < 1e-12:
            break
        app, aqq, apq = A[p][p], A[q][q], A[p][q]
        theta = 0.5 * math.atan2(2 * apq, aqq - app) if (aqq - app) != 0 else math.pi / 4
        c, s = math.cos(theta), math.sin(theta)
        for k in range(n):
            akp, akq = A[k][p], A[k][q]
            A[k][p] = c * akp - s * akq
            A[k][q] = s * akp + c * akq
        for k in range(n):
            apk, aqk = A[p][k], A[q][k]
            A[p][k] = c * apk - s * aqk
            A[q][k] = s * apk + c * aqk
        for k in range(n):
            vkp, vkq = V[k][p], V[k][q]
            V[k][p] = c * vkp - s * vkq
            V[k][q] = s * vkp + c * vkq
    best = max(range(n), key=lambda i: A[i][i])
    return [V[i][best] for i in range(n)]


def _quat_to_matrix(q):
    w, x, y, z = q
    return [
        [w * w + x * x - y * y - z * z, 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), w * w - x * x + y * y - z * z, 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), w * w - x * x - y * y + z * z],
    ]


def formula(atoms):
    """Hill-notation formula of a structure, for identity checks."""
    from collections import Counter
    c = Counter(s for s, _ in atoms)
    out = ""
    if "C" in c:
        out += f"C{c.pop('C')}" + (f"H{c.pop('H')}" if "H" in c else "")
    for sym in sorted(c):
        out += f"{sym}{c[sym]}"
    return out


def connectivity_delta(ref_atoms, cand_atoms):
    """Bonds present in one structure but not the other, as 1-based pairs."""
    a = bond_graph(ref_atoms)
    b = bond_graph(cand_atoms)
    gained = sorted(tuple(sorted(p)) for p in (b - a))
    lost = sorted(tuple(sorted(p)) for p in (a - b))
    return ([ (i + 1, j + 1) for i, j in gained ],
            [ (i + 1, j + 1) for i, j in lost ])

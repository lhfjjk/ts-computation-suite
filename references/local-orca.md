# Local ORCA execution (with an optional MLIP helper)

This reference is for running the protocol on a workstation with ORCA 6.1,
plus an optional learned-potential environment for geometry repair, conformer
folding, and pre-screening. It assumes the five-stage protocol in
[ts-protocol.md](ts-protocol.md); only the engine-specific parts are here.

## 1. Environment

- ORCA: `orca` on `PATH`, or set `TS_ORCA_EXE`/`ORCA_EXE`; the scripts also
  accept `--orca`. Example: `C:\ORCA_6.1.1\orca.exe` or
  `/opt/orca/orca`. For generated **bash** chains, set `ORCA_DIR` to a
  shell-style path (Git Bash: `/c/ORCA_6.1.1`; Linux: `/opt/orca`).
- MLIP (optional): a separate Python environment with `torch` + `orb-models`
  (or `mace-torch`), pointed to by `MLIP_PYTHON`. Example:
  `C:\Tools\mlip\env\Scripts\python.exe`. The core scripts do not need it;
  only the MLIP commands below do. Check with
  `$MLIP_PYTHON -c "import torch, orb_models"`.
- Put the ORCA directory on `PATH`, or set `ORCA_DIR` for generated chains.
- Use a pure-ASCII working directory when possible. Some viewers and older
  scripts fail silently on Chinese paths; if you must use one, verify every
  output file was actually written.

## 2. Level selection

| Use | Recommendation |
|---|---|
| Geometry + frequency, 3d metal | **TPSS-D3BJ/def2-SVP** (or r²SCAN-3c) + `RIJCOSX` + SMD |
| Geometry + frequency, organic | B3LYP-D3BJ/def2-SVP is acceptable; TPSS is also fine |
| Single point, production | **TPSS0-D4 or PBE0-D4 / def2-TZVP** |
| Key barrier, highest confidence | DLPNO-CCSD(T) single point on the DFT geometry |
| 3d metal, avoid | B3LYP for energetics; M06-2X for transition metals |

Rules:

- Optimization, frequency, and IRC must use the **same** functional, basis,
  grid, dispersion, and solvent.
- Single-point level must be **≥** optimization level.
- Add D3/D4 everywhere or nowhere; do not mix.
- For open-shell metals, enumerate all accessible spin states, optimize each,
  and compare **G** (not just E).

## 3. Input templates

### 3.1 Geometry + frequency

```text
! Opt Freq TPSS D3BJ def2-SVP RIJCOSX
%pal nprocs 12 end
%maxcore 1100
%cpcm
  smd true
  SMDsolvent "Toluene"
end
%geom
  MaxIter 220
end
* xyzfile 0 2 start.xyz
```

### 3.2 Constrained reference point (L1)

```text
! Opt Freq TPSS D3BJ def2-SVP RIJCOSX
%pal nprocs 12 end
%maxcore 1100
%cpcm
  smd true
  SMDsolvent "Toluene"
end
%geom
  MaxIter 220
  Constraints
    { B 20 56 1.750 C }
    { B 18 71 2.100 C }
  end
end
* xyzfile 0 2 start.xyz
```

The constrained point is a **reference point**, not a TS. Run `Freq` here so
the `.hess` exists for L2. Check the lowest modes before using it.

### 3.3 OptTS with the L1 Hessian (L2)

```text
! OptTS Freq TPSS D3BJ def2-SVP RIJCOSX
%pal nprocs 12 end
%maxcore 1100
%cpcm
  smd true
  SMDsolvent "Toluene"
end
%geom
  MaxIter 160
  InHess Read
  InHessName "P1-L1-constOptFreq.hess"
end
* xyzfile 0 2 P1-L1-constOptFreq.xyz
```

Without `InHess Read`, ORCA uses an Almlöf model Hessian and may climb the
wrong mode for 200–300 cycles. This is the most common "OptTS never converges"
cause in flexible systems.

### 3.4 IRC

```text
! IRC TPSS D3BJ def2-SVP RIJCOSX
%pal nprocs 12 end
%maxcore 1100
%cpcm
  smd true
  SMDsolvent "Toluene"
end
%irc
  MaxIter 20
  InitHess Read
  Hess_Filename "P1-L2-TSopt.hess"
  Direction both
  PrintLevel 2
  Trajectory True
end
* xyzfile 0 2 P1-L2-TSopt.xyz
```

IRC must start from the **converged TS geometry** and read the Hessian from
that same geometry. Internal constraints cannot be combined with IRC/NEB.

### 3.5 Broken-symmetry open-shell

```text
! UKS Opt Freq TPSS D3BJ def2-SVP RIJCOSX
%pal nprocs 12 end
%maxcore 1100
%moinp "reference.gbw"
%scf
  FlipSpin 75
  FinalMs 0
end
```

- `FlipSpin` is **1-based**; constraints and Mulliken rows are **0-based**.
- A BS single point must repeat `UKS + MORead + moinp + FlipSpin + FinalMs`.
  Without them ORCA runs RKS and reports a closed-shell singlet at the BS
  geometry.
- Report ⟨S²⟩ with every open-shell energy. A BS singlet near 0.2–0.4 is a
  different solution from a closed-shell singlet at 0.0.

## 4. Per-point scan (the ORCA 6.1.1-safe alternative)

Do not use `%geom Scan` for a long profile if a single point can hit
`MaxIter`; one failed point can destroy the whole scan. Use
`scripts/ts_scan_pts.py`:

```bash
python scripts/ts_scan_pts.py --xyz start.xyz --pair 59,75 \
    --start 2.05 --end 4.25 --points 8 --charge 0 --mult 2 \
    --method TPSS --basis def2-SVP --disp D3BJ --solvent Toluene \
    --constraints "59,5,2.12;59,7,2.08;59,60,2.42;59,81,3.40" \
    --outdir scan_r2 --run
```

It generates one `! Opt` input per point, starts each point from the previous
**converged** geometry, skips failed points, and writes
`scan_summary.tsv` with requested/actual distance, energy, convergence flag,
and relative energy. Use the highest converged point as the L1 guess.

If you must use `%geom Scan` for a quick look, remember:

```text
%geom
  Scan
    B 0 20 = 4.0, 2.0, 11
  end
  MaxIter 220
end
```

The `Scan` block needs its own `end`; otherwise `MaxIter` is swallowed and
ORCA exits with `Unknown error in GEOM block`.

## 5. Energy extraction

ORCA writes several energies. Use the **last** match:

```python
E_final = float(re.findall(
    r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", text)[-1])
E_cds = float(re.findall(
    r"Total Energy after SMD CDS correction\s*=\s*(-?\d+\.\d+)", text)[-1])
G = float(re.findall(
    r"Final Gibbs free energy[\s.:]*(-?\d+\.\d+)", text)[-1])
```

SMD convention (pick one and keep it for the whole study):

- Use `E_CDS` with thermal correction `G − E_CDS` (**recommended**), or
- use `E_FINAL` with `G − E_FINAL`.

Never use `E_FINAL` with `G − E_CDS`; that counts dispersion twice. The
difference is not cosmetic: a real BDFE shifted from +21.48 to +32.92 kcal/mol
when the two sides used different conventions.

For minima-vs-saddle comparisons use **ΔE**, not ΔG. The imaginary frequency
is excluded from the partition function, so the thermal correction of a TS is
not comparable to a minimum's.

## 6. Frequency and mode checks

- ORCA prints **3N** frequencies, including 6 translations/rotations. Do not
  expect 3N−6.
- Imaginary frequencies are negative. A genuine first-order TS has exactly one
  significant imaginary frequency (typically below −20 cm⁻¹).
- Use `scripts/ts_verify.py` to count imaginary modes and project each mode
  onto atom-pair distances. It reports which pairs move most; the **sign** of
  the eigenvector is arbitrary and must not be read as forming/breaking.
- Localization: the transferring atom's share of the total displacement.
  A true TS is often ~40%; `<15%` is a delocalized soft mode and should not be
  trusted as the reaction coordinate.
- A `SADDLE-1` that matches the declared pair can still be the wrong channel.
  Read the key bond lengths and run IRC.

## 7. Process discipline

- **Never run two ORCA jobs at once** on one workstation.
- Long jobs: launch detached (`Start-Process` on Windows), write stdout to a
  file, do not pipe to `tail` (the pipe hides the real exit code).
- Waiters must check both `orca.exe` **and** any project `run-*.sh` runner.
  A single "no orca.exe" check can slip into the gap between two chained jobs.
- Do not kill processes by a loose command-line substring; it can match the
  waiter itself or an unrelated generator. Match the exact script name and
  print the PIDs before killing.
- If a job is running, do not interrupt it to "fix" an input. Wait, or run a
  separate check that does not touch the running job.

## 8. MLIP pre-optimization and conformer folding

The learned potential replaces the old semiempirical folding step. Both commands
need the `MLIP_PYTHON` environment.

```bash
# fold a flexible ligand/chain: random rotamers -> MLIP relaxation -> dedup
$MLIP_PYTHON scripts/mlip_tools.py --charge 0 --spin 1 conf start.xyz \
    -o conf_run --n 32 --keep 5 --freeze-shell 1

# repair a malformed start (overlaps, stretched bonds) before DFT
$MLIP_PYTHON scripts/mlip_tools.py --charge 0 --spin 1 opt start.xyz \
    -o start.mlip.xyz --freeze-shell 1
```

Global options (`--charge`, `--spin`, `--device`, `--backend`) go **before**
the subcommand; the subcommand options follow it.

Conformer ranking uses the MLIP energies only to decide which 3–5 structures
reach DFT — re-optimize those at the DFT level before comparing them. Do
**not** quote MLIP energies or bond lengths in a paper: on this project's
Ni/Fe systems, conformer/stereoisomer energy differences from the potential
were off by ~2.6 kcal mol⁻¹, broken-symmetry singlets are outside its domain,
and spin-state gaps are not reliable.

Without an MLIP environment, use ORCA GOAT or another conformer generator;
the protocol does not depend on any one folding engine.

## 9. Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `OptTS` runs 200+ cycles, energy drifts up | Missing `InHess Read`; model Hessian wrong mode | Re-run L1 `Opt+Freq`, then L2 with `InHess Read` |
| Energy decreases monotonically for >50 cycles | Start is in the reactant/product basin | Stop; return to P1/P2; do not increase MaxIter |
| Two displaced runs return the same geometry | Real second-order saddle or flat valley | Stop; accept with disclosure, rigidify the spectator, or IRC along the reactive mode |
| `Freq` block missing | `OptTS` did not converge; ORCA skips frequencies | Run `Freq` on the last geometry, or raise MaxIter only if the mode is right |
| Imaginary mode does not touch the reactive pair | Constraint point's lowest mode is a spectator | Change constraint value/conformer; do not run L2 |
| Single point energy differs by ~0.09 Eh | SMD/CDS convention mixed | Recompute both sides with one convention |
| BS single point gives a closed-shell singlet | Missing UKS/MORead/FlipSpin | Repeat the full BS setup in the single point |
| Scan dies after one point | ORCA 6.1.1 scan MaxIter defect | Use `ts_scan_pts.py` per-point jobs |

For the full failure catalogue see [pitfalls.md](pitfalls.md).

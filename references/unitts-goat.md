# Candidate generation: UniTS-Gen, GOAT-TS, or manual

The five-stage protocol starts with a TS guess. This reference covers the three
practical sources of that guess. The rest of the suite (prefilter, L1/L2,
validation) is independent of which source you use.

| | UniTS-Gen (ML) | ORCA GOAT-TS (built in) | Manual / template |
|---|---|---|---|
| Install | Python env + torch + weights (~3 GB) | none (ORCA 6.1+) | none |
| Cost | generation is seconds; every candidate needs DFT | exploration costs CPU; the uphill stage can use ORCA's built-in GFN2-xTB (no separate install) | cheapest if the target is known |
| Strength | learned prior of what a TS looks like; change the reactive atoms to enumerate channels | physically rigorous exploration; can specify reacting atoms/bonds | fastest when you already know the step |
| Weakness | single-shot ~42% correct; must multi-sample; out of distribution for large BS organometallic | can be slow for 70+ flexible atoms | may miss conformers/channels |

**Decision guide:**

- Target step known, just need a guess → **manual + two-step protocol** first.
- Want to check whether you missed a channel/conformer → **GOAT-TS**.
- Need to enumerate reaction sites/stereochemistry, or GOAT is too slow →
  **UniTS-Gen multi-sample**.
- Any route → **prefilter → L1/L2 → P4 validation**. Never trust the generator.

---

## 1. UniTS-Gen

### 1.1 What it is

UniTS-Lib is a dataset of TS geometries; UniTS-Gen is a diffusion model that
takes a 2D molecular graph plus reactive-atom indices and generates 3D TS
guesses. The published workflow reports ~68.7% of single samples converging to
a first-order saddle, and ~41.9% passing human filtering + IRC. Multi-sampling
is the intended use: 10 samples per case, then prefilter and DFT.

### 1.2 Environment

The reference environment is Python 3.11 + PyTorch 2.4 + PyG + e3nn + RDKit +
QCBot + MolOP. The upstream project is tested on Ubuntu 22.04. A GPU is useful
but not required for small systems.

Generic installation outline (replace `<env>` and `<repo>`):

```bash
python3.11 -m venv <env>
<env>/bin/python -m pip install -U pip wheel setuptools

# code (use the current upstream instructions; tarball avoids git URL rewrites)
curl -sSL -o units.tar.gz https://codeload.github.com/licheng-xu-echo/UniTS/tar.gz/refs/heads/main
tar xzf units.tar.gz
curl -sSL -o qcbot.tar.gz https://codeload.github.com/licheng-xu-echo/QCBot/tar.gz/refs/heads/main
tar xzf qcbot.tar.gz

# MolOP ships inside the UniTS repo
cd UniTS && unzip -o MolOP-main.zip
cd MolOP-main && <env>/bin/python -m pip install -e . && cd ..
cd ../QCBot && <env>/bin/python -m pip install . && cd ../UniTS
<env>/bin/python -m pip install . \
  -f https://data.pyg.org/whl/torch-2.4.0+cu124.html \
  --extra-index-url https://download.pytorch.org/whl/cu124

# weights (modelscope CLI; no `python -m` entry)
<env>/bin/pip install modelscope
<env>/bin/modelscope download \
  --model 'XuLiCheng2025/UniTS-Gen-v1' --include 'model_path/units_hiegnn/*' \
  --local_dir ./units
# the CLI resolves the model from the installed package directory; copy the
# weights there as well if your install expects that layout
```

**Windows torch note:** if `import torch` fails with
`WinError 126 ... fbgemm.dll`, the wheel expects `libomp140.x86_64.dll` but only
ships `libiomp5md.dll`. A pragmatic workaround is a same-name copy of
`libiomp5md.dll`; it is not the official fix. Prefer a current MSVC/LLVM
runtime or a Linux environment for production.

### 1.3 CLI

**Index basis:** UniTS uses **0-based** `reactive_atom_idx`; this suite's
`--reactive` options use **1-based**. Convert explicitly.

```bash
# from SMILES
units-infer-smiles \
  --smiles 'CC(=O)C' \
  --reactive_atom_idx 5,12 \
  --charge 0 --multi 2 \
  --model_type units_hiegnn \
  --num_samples 10 \
  --output_dir ./ts_guess

# from xyz (the xyz is only used to rebuild the 2D graph)
units-infer-xyz \
  --xyz reactant.xyz \
  --reactive_atom_idx 0,3,5 \
  --charge 0 --multi 2 \
  --model_type units_hiegnn \
  --num_samples 10 \
  --output_dir ./ts_guess
```

Useful flags: `--seed`, `--batch_size`, `--save_full_trajectory True`,
`--save_combined_xyz`. `--multi` is the spin multiplicity; a wrong value puts
you on the wrong spin surface.

Output layout is typically:

```text
ts_guess/
  reaction_000/
    gen_0.xyz ... gen_9.xyz
    gen_all.xyz
```

### 1.4 Prefilter and DFT

```bash
python scripts/ts_prefilter.py --cand ts_guess/reaction_000/gen_all.xyz \
    --reactive 21,57 --forming 21,57 --write-pass pass.xyz
```

Then feed only the passing candidates to `ts_orca_prep.py` or `ts_ladder.py`.
The prefilter checks hard overlaps, clash ratios, Kabsch RMSD, and connectivity
changes. `--far-cutoff` only makes sense for forming bonds; a breaking bond is
allowed to be long.

### 1.5 Known limitations

- Single-shot correctness is ~42%; always multi-sample.
- The model is out of distribution for large, flexible, broken-symmetry
  organometallic systems. In one 84-atom Ni BS test, 4 candidates passed the
  geometric prefilter but their heavy-atom RMSD to the reference was
  3.6–4.4 Å; none had both the Ni–C bond and the incipient C···C contact.
- It does not validate connectivity, spin state, or IRC. Those remain the
  user's job.

---

## 2. ORCA GOAT-TS

ORCA 6.1+ includes GOAT (Global Optimization of Atomistic sTructures) with
`GOAT-TS` and `GOAT-REACT` modes. It is a physically rigorous alternative to an
ML prior and requires no extra install.

### 2.1 Verified syntax

These were verified against the ORCA 6.1.1 binary with a "no coordinates → stop"
probe:

- `! GOAT-TS`, `! GOAT`, `! GOAT-REACT` are valid simple-input keywords.
- Accepted `%GOAT` parameters: `USEGOAT`, `GFNUPHILL`, `MAXCORESOPT`,
  `MAXITER`, `NWORKERS`, `REACTIVEATOMS`, `RANDOMSEED`, `MAXEN`.
- `REACTIVEATOMS` uses braces: `REACTIVEATOMS {0 4 5} end`; comma/space-only
  forms are rejected. Indices are **0-based**.
- `GFNUPHILL GFN2XTB` (or `GFN1XTB`) uses ORCA's built-in GFN2-xTB for the
  cheap uphill stage; it needs no separate xtb binary.
- Rejected in `%GOAT`: `NSTEPS`, `RUNMODE`, `NPROCS`, `MAXCORE`, `TEMP`,
  `DOSOLV`, `OPTIMIZER`, `PRINTLEVEL`, `FIXEDATOMS`, `ENSEMBLE`, `MAXCONF`,
  `ETOL`, `GFN`, `COARSE`, `UPHILL`, `SCAN_GOAT_*`.

**Warning:** strings found in the binary (`SCAN_GOAT_REACTIVE_ATOMS`,
`Scan_GOAT_TOPO_BREAK`, etc.) can be internal function names, not user
parameters. Probe the input parser before relying on them.

### 2.2 Template

```text
! GOAT-TS B3LYP D3BJ def2-SVP TightSCF
%pal nprocs 12 end
%maxcore 1100

%GOAT
  USEGOAT true
  MAXITER 30
  REACTIVEATOMS {0 4 5} end
end

* xyz 0 1
<coords>
*
```

Use ORCA's built-in GFN2-xTB uphill stage (`GFNUPHILL GFN2XTB`) for a cheaper
exploration, then re-optimize the candidate at the DFT level. GOAT is a search
tool, not a substitute for `OptTS` + `Freq` + IRC.

### 2.3 Cost and constraints

- GOAT can be slow for 70+ flexible atoms; it is best for small/medium systems
  or for checking whether a channel exists.
- IRC and NEB cannot be combined with internal constraints in ORCA. Whether
  GOAT and internal constraints are mutually exclusive is not established;
  test it before relying on it.
- After GOAT returns a candidate, run the normal P2/P3/P4 gates. GOAT does not
  tell you whether the candidate is the channel you intended.

---

## 3. Manual / template

When the target step is known, a hand-built guess plus the two-step protocol is
often the fastest and most reliable route, especially for large organometallic
systems where ML generation is out of distribution.

Build the guess by:

1. taking the reactant/encounter geometry;
2. moving the reactive atoms to plausible TS distances (e.g. forming bond
   ~2.0–2.4 Å, breaking bond ~1.8–2.2 Å);
3. checking connectivity and clashes with `ts_prefilter.py`;
4. running L1 (constrained Opt+Freq) → L2 (OptTS + Hessian) → P4.

The guess does not need to be perfect; it needs to be in the right basin and
have the right reactive atoms. That is what the L1 gate checks.

---

## 4. Reporting an AI/automated candidate

Report:

- tool and version (UniTS-Gen commit/weights, ORCA version for GOAT);
- number of candidates, prefilter pass rate;
- how many DFT jobs converged, how many were first-order saddles;
- how many passed IRC and the key bond-length check;
- human interventions (manual guesses, changed constraints/spin).

Do not write "AI found the TS" without this breakdown. The distinction between
"generated a candidate" and "validated a TS" is the whole point.

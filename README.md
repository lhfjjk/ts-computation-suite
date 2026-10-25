# TS Computation Suite

A portable, validation-first protocol for **AI-assisted transition-state (TS)
discovery**. It is built for chemists who can generate TS candidates with a
model (UniTS-Gen, ORCA GOAT-TS, or a hand-built guess) but need a reliable way
to decide whether a candidate is a real, connected, first-order saddle for the
reaction they intended.

The entry point is [SKILL.md](SKILL.md). Repository: https://github.com/lhfjjk/ts-computation-suite

The suite is self-contained: no
author-specific paths, accounts, or auxiliary skills are required.

## What it is

- A four-stage protocol: **P0 chemical preflight → P1 barrier scan → P2
  constrained reference + Hessian → P3 OptTS → P4 validation**.
- Local ORCA/xtb execution and Gaussian/Slurm execution, with templates.
- Scripts for candidate prefiltering, per-point independent scans, ORCA/Gaussian
  input generation, imaginary-frequency/IRC checks, and an environment doctor.
- A failure catalogue and printable checklists distilled from real Ni/Fe
  catalysis projects.

## What it is not

- It is not a generative model. Candidate generation can come from UniTS-Gen,
  ORCA GOAT-TS, or a manual guess; see
  [references/unitts-goat.md](references/unitts-goat.md).
- It does not replace electronic-structure judgment: spin states, connectivity,
  mode identity, and IRC endpoints still need to be checked.
- It does not guarantee that a candidate is the target channel. A `SADDLE-1`
  with the right atom pair can still be the wrong reaction.

## Quick start

```bash
# 1) environment check (no quantum-chemistry job is run)
python scripts/doctor.py

# 2) local ORCA: per-point independent scan
python scripts/ts_scan_pts.py --xyz reactant.xyz --pair 21,57 \
    --start 1.40 --end 2.60 --points 6 --charge 0 --mult 2 \
    --method TPSS --basis def2-SVP --disp D3BJ --solvent Toluene \
    --outdir scan_P1 --run

# 3) ORCA L1 -> L2 -> L3 -> SP ladder
python scripts/ts_ladder.py --name P1 --xyz scan_P1/p04/point.xyz \
    --forming 21,57 --hold 21,57=1.75 --charge 0 --mult 2 \
    --method TPSS --basis def2-SVP --solvent Toluene --outdir ladder_P1
python scripts/ts_ladder_run.py --dir ladder_P1 --reactive 21,57

# 4) Gaussian/Slurm route
python scripts/gauss_ts_prep.py --xyz guess.xyz --name g6TS \
    --charge 0 --mult 2 --method UB3LYP --basis "6-31G(d)" \
    --disp "em=gd3bj" --solvent toluene \
    --constraints "1,2,1.500" --nproc 24 --mem 48GB \
    --outdir gauss_g6TS --slurm --partition community
```

See [references/setup.md](references/setup.md) for installation, environment
variables, and a five-minute smoke test.

## Repository layout

```text
SKILL.md                  entry point and protocol summary
references/               setup, protocol, local ORCA, HPC, pitfalls,
                          checklists, UniTS/GOAT, background, AI role
scripts/                  doctor, prefilter, scan, ladders, verifiers
assets/templates/         ORCA/Gaussian/Slurm input templates
```

## Design principles

1. Move failure into cheap calculations.
2. Reuse Hessians, conformers, and spin states.
3. Stop early when the evidence says the target is not a first-order saddle.
4. Keep the output reproducible: report method, spin, frequencies, IRC, and
   uncertainty.

## Citation

A companion preprint describing the protocol and case evidence is in
preparation for ChemRxiv; the draft source and PDF are in
[`paper/`](paper/). Until it is available, cite this repository:

```text
Li, H.-F. TS Computation Suite: a validation-first protocol for AI-assisted
transition-state discovery. GitHub repository, 2026.
```

## License

MIT. See [LICENSE](LICENSE).

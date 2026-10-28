# Changelog

All notable changes to this skill package are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project adheres to
[Semantic Versioning](https://semver.org/).

The version of record is the `version:` field in [`CITATION.cff`](CITATION.cff).
Every release carries a `vX.Y.Z` tag and a GitHub Release whose notes are copied
from the section below.

## [0.2.2] - 2026-10-28

### Changed

- Repositioned the package as an **AI-assisted computational toolkit whose core
  function is transition-state discovery**. The README and the `SKILL.md`
  description now also cover reaction energies, bond dissociation free energies,
  conformer ensembles, spin-state checks, geometry repair, and job
  orchestration — the same gates and rules serve those jobs.
- `docs/INTRO-VIDEO-SCRIPT.md` rewritten around that scope: a new section on the
  non-TS uses (one energy convention, conformer ensembles first, spin-state
  enumeration), with an updated 30-second cut, generator prompt, shot list and
  fact-check table.

## [0.2.1] - 2026-10-28

### Added

- `docs/INTRO-VIDEO-SCRIPT.md` — a three-minute narration script, a 30-second
  short version, a ready-to-paste prompt for an AI video generator, and a
  fact-check table for the on-screen numbers.

### Fixed

- Documentation now says **five stages (P0–P4)** everywhere. The README,
  `SKILL.md`, the cover letter, and three reference files still described the
  protocol as four-stage while the manuscript and the P-labeling used five; the
  manuscript's revised text is the reference.

## [0.2.0] - 2026-10-28

### Added

- **P0.5 learned-potential pre-screen** (`scripts/ts_mlip_screen.py`). Pushes a
  candidate along the declared reaction coordinate in both directions on a
  foundation MLIP surface and classifies it as `BARRIER-LIKE`, `BASIN`,
  `WRONG-CHANNEL`, or `PES-MISMATCH`. `--bs` restricts the output to geometry
  columns for broken-symmetry singlets.
- `scripts/mlip_preopt.py` — geometry repair: removes hard overlaps and pulls
  the declared bonds to their target lengths while holding the coordination
  sphere; includes an optional chain-runner hook.
- `scripts/mlip_tools.py` — shared OrbMol/MACE backend, constraints, conformer
  search.
- Measured applicability domain in the manuscript: five DFT-validated
  broken-symmetry saddles all fell outside the learned surface (83–85
  kcal mol⁻¹, `PES-MISMATCH`); one DFT minimum was correctly classified as
  `BASIN`.

### Changed

- `SKILL.md` — the protocol is now "four mandatory stages + one optional
  accelerator"; new P0.5 section with verdict table, applicability domain, and
  the note that the stage does not save DFT cycles.
- `README.md` — documents the optional MLIP requirements (`MLIP_PYTHON`) and the
  versioning policy.
- `paper/manuscript.tex` — new P0.5 paragraph, case-evidence paragraph, rows in
  the stop-rule, division-of-labor and case tables, updated workflow figure, and
  three new references (Husistein & Reiher 2026; Orb; OMol25).

### Notes

- The pre-screen does not reduce DFT cycles: ORCA projects the constraints onto
  the input geometry before the first optimization cycle (measured, Ni–C\*
  3.000 → 2.150 Å), and the remaining cycles are ligand-conformer relaxation.
- No positive control exists yet for `BARRIER-LIKE`; verify the first case by
  hand before trusting the verdict.

## [0.1.0] - 2026-10-25

### Added

- Initial release: the P0–P4 validation-first protocol, ORCA and Gaussian ladder
  generators and runners, output verifiers, a structure preflight checker, an
  environment doctor, and the ChemRxiv draft with its supporting information.

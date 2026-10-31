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

- **The core function:** a five-stage protocol for transition-state discovery —
  **P0 chemical preflight → P1 barrier scan → P2 constrained reference +
  Hessian → P3 OptTS → P4 validation**.
- **The same gates for other jobs:** reaction energies and bond dissociation
  free energies (constrained reference points, frequency/thermochemistry
  checks, one fixed energy convention), conformer ensembles, spin-state
  enumeration, and geometry repair before DFT.
- **The plumbing:** local ORCA execution and Gaussian/Slurm execution with
  templates, an optional learned-potential layer (geometry repair, conformer
  folding, candidate pre-screen), candidate prefiltering, per-point independent
  scans, input generation, imaginary-frequency/IRC/⟨S²⟩ checks, and an
  environment doctor.
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

## Quick start (agent deployment)

1. Put this folder in your agent's skill directory:
   - ZCode: `~/.zcode/skills/ts-computation-suite/`
   - Codex: `~/.codex/skills/ts-computation-suite/`

   ```bash
   git clone https://github.com/lhfjjk/ts-computation-suite \
       ~/.zcode/skills/ts-computation-suite
   ```

2. Start a new session so the agent discovers the skill.

3. Invoke it by name (`$ts-computation-suite`) or with a request such as:

   > Use `$ts-computation-suite` to validate this TS candidate. Start with the
   > P0 declarations, run the five mandatory stages (and the optional MLIP
   > pre-screen where it applies), and report the imaginary frequency, mode
   > localization, IRC endpoints, and spin diagnostics.

4. On first use, run `python scripts/doctor.py`. If a tool is not found, set
   `TS_ORCA_EXE`, `ORCA_DIR`, or `TS_G16_EXE`; see
   [references/setup.md](references/setup.md). The optional MLIP layer (repair,
   conformer folding, P0.5 screen) additionally needs PyTorch and `orb-models`
   in a separate environment, pointed to by `MLIP_PYTHON`; the rest of the
   package does not.

5. The agent reads [SKILL.md](SKILL.md) for the protocol and
   [references/](references/) for details. No manual command sequence is
   required. The agent should not submit a DFT job until the P0 declarations
   are complete.

For manual execution, the exact ORCA/Gaussian commands are in
[references/local-orca.md](references/local-orca.md) and
[references/hpc-slurm.md](references/hpc-slurm.md).

## Repository layout

```text
SKILL.md                  entry point and protocol summary
references/               setup, protocol, local ORCA, HPC, pitfalls,
                          checklists, UniTS/GOAT, background, AI role
scripts/                  doctor, prefilter, scan, ladders, verifiers,
                          optional MLIP repair/conformer/screen tools
assets/templates/         ORCA/Gaussian/Slurm input templates
```

## Design principles

1. Move failure into cheap calculations — including the learned potential: use
   MLIP repair and screening before DFT wherever its applicability domain allows.
2. Reuse Hessians, conformers, and spin states.
3. Stop early when the evidence says the target is not a first-order saddle.
4. Keep the output reproducible: report method, spin, frequencies, IRC, and
   uncertainty.

## Versioning

Releases follow [Semantic Versioning](https://semver.org/):

- **MAJOR** — a gate's pass/fail semantics or a script's command-line interface
  changes incompatibly.
- **MINOR** — a new stage, script, or capability that is backward compatible.
- **PATCH** — fixes, documentation, and wording only.

The version of record is `version:` in [`CITATION.cff`](CITATION.cff). A release
is only complete when, in the same change:

1. `version:` and `date-released:` are bumped in `CITATION.cff`;
2. a dated section is added to [CHANGELOG.md](CHANGELOG.md);
3. the commit is tagged `vX.Y.Z`;
4. a GitHub Release is published, with the changelog section as its notes.

The protocol text, the figure, and the scripts that implement it are versioned
together: the paper in [`paper/`](paper/) ships as part of the same release.

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

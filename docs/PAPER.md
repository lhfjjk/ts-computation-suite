# The companion preprint

The manuscript is **not** shipped inside this skill package. It lives in its own
repository:

| | |
|---|---|
| Repository (local) | `C:\Users\Administrator\.zcode\repos\ts-computation-suite-paper` |
| Extracted from | this repository @ `a78b77b` (v0.4.2), subtree `paper/` |
| Canonical working copy | the author's LaTeX project (`chemrxiv-ai-ts/`); its `submission/` build is the QA-passed ChemRxiv version (13 pages, SI 6 pages) |

**Why it was split out.** The paper added ~1.2 MB of PDFs and LaTeX sources to a
package whose job is to run calculations; agents load `SKILL.md` and the
`references/` files, never `paper/`. The two also aged at different speeds: after
v0.4.0 replaced `xtb` with the learned potential in the skill, the shipped
manuscript text still read "ORCA/xtb" in six places.

**If you are an agent reading this:** do not quote numbers from the preprint
snapshot. Take the protocol from `SKILL.md` and `references/`, and cite the paper
only after its ChemRxiv record exists.

## Where each artifact lives now

| Artifact | Location |
|---|---|
| protocol summary | `SKILL.md` (this repository) |
| full protocol, stop rules, checklists | `references/ts-protocol.md`, `references/checklists.md` |
| ORCA/Gaussian templates | `assets/templates/`, `references/local-orca.md`, `references/hpc-slurm.md` |
| manuscript source + PDFs + figure generator | the separate paper repository |
| figure generator (copy that ships with the skill) | `assets/templates/` is input-only; the figure is built from the paper repository's `make_figure.py` |

---
name: ts-computation-suite
description: >-
  End-to-end transition-state (TS) computation for chemists: accelerated
  five-stage protocol (chemical preflight, barrier scan, constrained
  reference/Hessian, OptTS, validation), local ORCA/xtb runs, Gaussian/Slurm
  HPC submissions, candidate prefiltering, imaginary-frequency/IRC checks, and
  AI-assisted candidate generation. Use when planning or running a TS search,
  choosing between local and HPC execution, writing OptTS/opt=ts inputs,
  diagnosing a failed TS, or comparing barriers and selectivities. Not for
  reactions without a defined elementary step, or for excited-state/nonadiabatic
  dynamics.
metadata:
  short-description: Accelerated TS search with local ORCA and Slurm Gaussian
---

# TS Computation Suite

A practical, engine-aware kit for finding and validating transition states
without burning days on a wrong guess. It packages the protocol, the local
ORCA/xtb commands, the Gaussian/Slurm route, the failure-mode gates, and the
AI-assisted parts that actually worked in real projects.

**The one idea:** acceleration does not come from a magic keyword. It comes
from moving failure into cheap calculations, reusing Hessians and conformers,
and stopping early when the evidence says the target is not a first-order
saddle.

**What is not bundled:** a generative model or a cluster account. Candidate
generation can come from UniTS-Gen, ORCA GOAT-TS, or a hand-built guess; this
suite provides the prefilter, the ladder, and the validation layer. The
UniTS-Gen/GOAT interfaces are documented in
[references/unitts-goat.md](references/unitts-goat.md); the AI division of
labor is in [references/ai-role.md](references/ai-role.md).

**First run:** `python scripts/doctor.py` checks Python, compiles every script,
locates ORCA/Gaussian/xtb/Slurm/SSH, and runs four generator/parser self-tests.
Setup and environment variables are in
[references/setup.md](references/setup.md).

## 中文速查

- **先写清**：反应对、电荷/多重度、所有可及自旋态、连通性、构象、参考态。
- **本地 ORCA**：`ts_prefilter.py` → `ts_scan_pts.py --run` → `ts_ladder.py`
  → `ts_ladder_run.py` → `ts_verify.py`。
- **超算 Gaussian**：`gauss_ts_prep.py --slurm` → `sbatch submit.sbatch`
  → `gauss_verify.py`。
- **看到就停**：P2 最低模不牵动反应对；OptTS 能量单调下降 >50 周期；两次位移
  回到同一鞍点；2–3 个 −15…−25 cm⁻¹ 浅虚频；虚频定域度 <15%；`SADDLE-1`
  但关键键长不对。
- **详细文档**：[协议](references/ts-protocol.md)、
  [本地 ORCA](references/local-orca.md)、
  [超算 Gaussian](references/hpc-slurm.md)、
  [避坑](references/pitfalls.md)、
  [AI 作用](references/ai-role.md)。

## Scope and boundary

Use this skill when you already know the elementary step you want:

- you can name the bonds that form and break;
- you know the charge and multiplicity, and have enumerated the accessible
  spin surfaces;
- you have a reasonable reactant/encounter geometry and a candidate product.

If the mechanism itself is unknown (competing channels, unknown spin surface,
unknown protonation state), first enumerate the chemistry. A faster OptTS on
the wrong channel is still the wrong answer.

Out of scope: excited-state/nonadiabatic dynamics, strong multireference
systems where a single DFT surface is not meaningful, and transport-limited
reactors where the intrinsic barrier is not the observed rate.

## Before any calculation: five declarations

Write these down before generating an input. Most failed TS projects fail here.

1. **Reactive pairs** — which bonds form, which break, and their actual
   distances in the starting geometry. A "reactive pair" at 4.5 Å is usually a
   mistake.
2. **Charge, multiplicity, spin surface** — for a metal, list every accessible
   dⁿ state (e.g. d⁷ doublet/quartet, d⁸ singlet/triplet) plus any broken-symmetry
   solution. Find the TS on each surface; do not default to one.
3. **Connectivity** — build a bond graph with `1.3 × (r_i + r_j)` and check that
   the intended bonds are present and no unplanned bond has formed. Three
   real projects were lost to a "reactant" that had already rearranged.
4. **Conformers** — for flexible ligands/chains, run a cheap xtb/GFN-FF
   conformer search and keep 3–5 low-energy conformers. A single-conformer
   R/S energy difference can be conformer-dominated (observed RMSD 1.62 Å).
5. **Reference state** — the reactant/encounter complex and the product, at the
   same level as the TS. A barrier without a defined reference is not a barrier.

## Choose the execution route

| Situation | Route | Why |
|---|---|---|
| Small/medium system, fast iteration, ORCA available | **Local ORCA/xtb** | Cheapest feedback; best diagnostics (`ts_verify.py`, mode projection) |
| Production runs, Gaussian-only method, larger basis | **HPC Gaussian/Slurm** | Throughput and memory; use `gauss_ts_prep.py` |
| Cluster has ORCA and you want the same protocol | **HPC ORCA/Slurm** | Use `ts_orca_prep.py` or `ts_ladder.py`, wrap with the Slurm template |
| 100+ atoms, many conformers | **Local pre-screen + HPC production** | Do not burn queue time on candidates that fail cheap gates |
| Unknown channel / no TS guess | **GOAT-TS or a relaxed scan first** | Do not start with a blind OptTS |

Local ORCA details: [references/local-orca.md](references/local-orca.md).
Gaussian/Slurm details: [references/hpc-slurm.md](references/hpc-slurm.md).

## The protocol: five mandatory stages (P0–P4) + one optional accelerator

### P0 — Chemical preflight (0.5–2 h)

Declarations above, plus a candidate prefilter. Use
`scripts/ts_prefilter.py` to reject hard overlaps, clashes, wrong connectivity,
and candidates far from the intended reactive distance. If the candidate came
from a generative model, this is where most of them should die.

### P0.5 — Learned-potential pre-screen (optional; preferred where applicable) (seconds per candidate)

Before spending hours on P1, run the candidate through a foundation machine-learned
potential: `scripts/ts_mlip_screen.py`. It works at the geometry level only.

1. **Repair** — remove hard overlaps and pull the declared bonds to their target
   lengths (`scripts/mlip_preopt.py`). Measured: 0.83–0.93 Å overlaps → 2.07–2.97 Å,
   coordination sphere unchanged, radical face rotated by only ~5°.
2. **Screen** — push the candidate along the declared reaction coordinate in both
   directions and relax on the MLIP surface.

| verdict | meaning | action |
|---|---|---|
| `BARRIER-LIKE` | two distinct endpoints; bond-graph delta = the declared forming/breaking pairs | proceed to P1 |
| `BASIN` | both descents reach the same minimum | drop the candidate |
| `WRONG-CHANNEL` | endpoints differ, but not by the declared pairs | return to P0 and re-check the reactive pairs |
| `PES-MISMATCH` | candidate lies >60 kcal above both endpoints | the model does not support this reaction — go to P1, do not use this verdict |
| `SKIP-BS` | `--bs` declared (broken-symmetry singlet) | geometry columns only |

**Applicability domain — measured, not assumed.** Five DFT-validated saddles
(S_H2 coupling, −574…−581 cm⁻¹, 84–98 atoms, broken-symmetry singlets) were all
`PES-MISMATCH` (83–85 kcal): the MLIP surface relaxes the biradical straight to the
coupled product. One DFT minimum was correctly `BASIN`. So on broken-symmetry
organometallics this stage is a **geometry-repair tool only** — it never earns a "go".
On closed-shell or high-spin systems it may gate candidates, but no positive control
exists yet: verify the first case by hand.

**It does not save DFT cycles by satisfying constraints.** ORCA projects the constraints
onto the input geometry before cycle 1 (measured: Ni–C* 3.000 → 2.150 Å in the first
printed cycle), and the 142–265 cycles that follow are ligand-conformer relaxation.
The pre-screen repairs bad starts; it does not replace conformer search.

Requirements: `torch` + `orb-models` in a separate environment, passed via
`MLIP_PYTHON`. The core scripts stay standard-library only; this stage is strictly optional.

### P1 — Barrier-region scan (1–20 h)

Do **not** start with a blind `OptTS`. Scan the reaction coordinate with 4–8
points. For ORCA, use `scripts/ts_scan_pts.py`: each point is an independent
`! Opt` job, constrained at a fixed distance, starting from the previous
converged geometry. This avoids the ORCA 6.1.1 failure mode where one point
hitting `MaxIter` destroys the entire `%geom Scan`.

Use only `CONVERGED` points. The highest converged point is the TS guess; the
whole curve tells you whether the channel is monotonic, has a shallow bump
(possible conformer change), or is a flat valley.

For Gaussian, the same idea: a sequence of `opt=modredundant` constrained
optimizations, each reading the previous geometry. Do not use a single
relaxed scan if a point failure would lose the whole job.

### P2 — Constrained reference + Hessian (2–6 h)

At the highest point, run a constrained `Opt + Freq` at the chosen distance.
This produces two things:

- a **reference point** (not a TS) that is already usable for relative
  comparisons if P3 fails;
- a **`.hess` / checkpoint Hessian** for the next step.

**Gate:** the lowest modes at this point must involve the declared reactive
pair. If they do not, do not run OptTS; change the constraint, the conformer,
or the guess.

### P3 — OptTS (4–12 h)

ORCA: `! OptTS Freq` with `%geom InHess Read / InHessName "<P2>.hess"`.
Gaussian: `opt=(ts,readfc,noeigen) freq` with `geom=check guess=read` from the
L1 checkpoint.

The Hessian reuse is the single biggest reliability gain. Without it, ORCA
falls back to a model Hessian and can climb the wrong mode for hundreds of
cycles.

If a small displacement is needed, move along the imaginary mode but scale by
the **target bond-length change (about 0.05 Å)**, not by the maximum atomic
displacement. A 0.18 Å maximum displacement can change the reactive bond by
−0.31 Å and throw the geometry into the reactant basin.

### P4 — Validation (2–6 h)

Do not report a barrier until all of these pass:

- exactly one significant imaginary frequency (typically below −20 cm⁻¹);
- the imaginary mode is localized on the reactive pair (transfer atom
  contributes a substantial fraction of the displacement; <15% is suspicious);
- the key bond lengths match the intended step (a `SADDLE-1` that matches the
  declared pair can still be a different channel);
- IRC both directions (or a downhill relaxation) connects the intended
  reactant and product;
- ⟨S²⟩ and Mulliken spin are consistent with the geometry's electronic state;
- energies use the same level for optimization/frequency/IRC, and ΔE for
  minima-vs-saddle comparisons.

Full protocol and stop rules: [references/ts-protocol.md](references/ts-protocol.md).

## Quick start

### Local ORCA

```bash
# 1) prefilter candidate guesses
python scripts/ts_prefilter.py --cand guesses.xyz --reactive 21,57 \
    --forming 21,57 --write-pass passing.xyz

# 1b) optional MLIP pre-screen (needs torch + orb-models; MLIP_PYTHON points to that env)
$MLIP_PYTHON scripts/ts_mlip_screen.py --cand passing.xyz \
    --forming 21,57 --charge 0 --spin 2 --out screen/

# 2) per-point independent scan (generates and runs)
python scripts/ts_scan_pts.py --xyz passing.xyz --pair 21,57 \
    --start 1.40 --end 2.60 --points 6 --charge 0 --mult 2 \
    --method TPSS --basis def2-SVP --disp D3BJ --solvent Toluene \
    --constraints "1,2,1.95" --outdir scan_P1 --run

# 3) L1 -> L2 -> L3 -> SP ladder
python scripts/ts_ladder.py --name P1 --xyz scan_P1/p04/point.xyz \
    --forming 21,57 --hold 21,57=1.75 --charge 0 --mult 2 \
    --method TPSS --basis def2-SVP --solvent Toluene --outdir ladder_P1
python scripts/ts_ladder_run.py --dir ladder_P1 --reactive 21,57

# 4) verify
python scripts/ts_verify.py --out ladder_P1/P1-L2-TSopt/P1-L2-TSopt.out \
    --reactive 21,57
```

### HPC Gaussian

```bash
python scripts/gauss_ts_prep.py --xyz guess.xyz --name g6TS \
    --charge 0 --mult 2 --method UB3LYP --basis "6-31G(d)" \
    --disp "em=gd3bj" --solvent toluene \
    --constraints "1,2,1.500" --nproc 24 --mem 48GB \
    --outdir gauss_g6TS --slurm --partition community

# on the cluster
cd gauss_g6TS && sbatch submit.sbatch

# after 02 finishes
python scripts/gauss_verify.py --out 02_ts_optfreq.log
```

## Stop rules

These are the highest-value part of the kit. They prevent 10–20 h dead ends.

| # | Symptom | Action |
|---|---|---|
| 1 | P2 lowest mode does not involve the reactive pair | Do not run OptTS; change constraint/conformer |
| 2 | OptTS energy decreases monotonically for >50 cycles and the reactive bond is unchanged | Wrong basin; stop and return to P1/P2 |
| 3 | Two displaced OptTS runs return to the same geometry/energy (<0.01 kcal) and the soft mode persists | Stop; do not run a third. Accept with disclosure, rigidify the spectator, or run IRC along the reactive mode |
| 4 | 2–3 shallow imaginary modes around −15…−25 cm⁻¹ | Flat transition valley; use the last geometry's Freq as a reference barrier and say so |
| 5 | Imaginary-mode localization <15% | Delocalized soft mode, not the target coordinate |
| 6 | `SADDLE-1` and the declared pair matches, but key bond lengths are wrong | Wrong channel; do not report this barrier |
| 7 | xtb and DFT disagree | DFT wins; xtb is for folding only |
| 8 | IRC/relaxation endpoint is not the intended reactant/product | Return to P1/P2; no barrier claim |
| 9 | Single-point ⟨S²⟩ differs from the geometry job | Wrong electronic state; redo with the same UKS/BS settings |
| 10 | R/S difference without a conformer search | Not a stereoselectivity result; do the ensemble first |

## Energy and spin rules

- **Optimization / frequency / IRC must use the same functional, basis, grid,
  dispersion, and solvent.** Do not switch levels inside a path.
- **Single-point level ≥ optimization level.** For 3d metals, optimize with
  TPSS-D3BJ/def2-SVP (or r²SCAN-3c) and single-point with TPSS0-D4 or
  PBE0-D4/def2-TZVP; use DLPNO-CCSD(T) for key barriers if needed.
- **ΔE for minima-vs-saddle comparisons.** A saddle's imaginary frequency is
  excluded from the partition function, so its thermal correction is not
  comparable to a minimum's. Use ΔG only within the same type of stationary
  point.
- **ORCA SMD:** choose one convention and keep it. If you use `E_CDS`, the
  thermal correction must be `G − E_CDS`; if you use `E_FINAL`, use
  `G − E_FINAL`. Mixing them double-counts or omits dispersion and can shift a
  barrier by 5–11 kcal/mol.
- **Broken-symmetry single points must inherit `UKS + MORead + moinp +
  FlipSpin + FinalMs`.** Otherwise ORCA defaults to RKS and you compute a
  closed-shell singlet at an open-shell geometry.
- **Report ⟨S²⟩ with every open-shell number.** A BS singlet near 0.2–0.4 is
  not the same solution as a closed-shell singlet at 0.0.

## Scripts

| Script | Role |
|---|---|
| `scripts/doctor.py` | Environment check + four generator/parser self-tests |
| `scripts/ts_prefilter.py` | Candidate prefilter: clash, overlap, Kabsch RMSD, connectivity |
| `scripts/ts_scan_pts.py` | ORCA per-point independent constrained scan (no `%geom Scan`) |
| `scripts/ts_ladder.py` | Generate ORCA L1 constrained Opt+Freq, L2 OptTS, L3 relax, SP |
| `scripts/ts_ladder_run.py` | Run the ORCA ladder serially with gates and routing |
| `scripts/ts_verify.py` | ORCA verdict: imaginary count, mode atoms, IRC |
| `scripts/ts_orca_prep.py` | Turn candidate xyz files into ORCA inputs |
| `scripts/gauss_ts_prep.py` | Generate Gaussian L1/L2/IRC/SP ladder + sbatch |
| `scripts/gauss_verify.py` | Gaussian verdict: frequencies, termination, thermo, ⟨S²⟩ |
| `scripts/preflight_check.py` | Structure sanity gate used by the ladder runner |
| `scripts/ts_mlip_screen.py` | Optional MLIP pre-screen: basin / barrier / channel verdict at geometry level |
| `scripts/mlip_preopt.py` | Optional MLIP geometry repair: remove overlaps, satisfy declared bond constraints |
| `scripts/mlip_tools.py` | Shared MLIP backend (OrbMol / MACE calculators, constraints, conformer search) |

## Reporting checklist

A TS result is only as good as its evidence trail. Report:

- method, basis, dispersion, solvent model, grid, and ORCA/Gaussian version;
- charge, multiplicity, spin surface, and ⟨S²⟩;
- the imaginary frequency and which atoms move;
- the key bond lengths at the TS;
- the IRC/relaxation endpoints and their connectivity;
- the reference state used for the barrier;
- conformer search and Boltzmann averaging, if selectivity is claimed;
- the uncertainty introduced by the method (do not report 0.1 kcal
  differences as if they were meaningful when the method error is several
  kcal/mol).

## References

Read only what the current task needs:

- [references/setup.md](references/setup.md) — install the skill, set
  environment variables, run the doctor, five-minute smoke test.
- [references/ts-protocol.md](references/ts-protocol.md) — full five-stage
  protocol, stop rules, route decision tree, time budget.
- [references/local-orca.md](references/local-orca.md) — ORCA 6.1 input
  templates, per-point scans, energy/frequency parsing, process discipline.
- [references/hpc-slurm.md](references/hpc-slurm.md) — Gaussian/Slurm
  templates, a worked cluster example, checkpoint semantics, queue/cores.
- [references/pitfalls.md](references/pitfalls.md) — the accumulated failure
  catalogue (structure, spin, energy, TS diagnosis, scan, HPC, candidate
  generation).
- [references/checklists.md](references/checklists.md) — printable P0–P4,
  HPC submission, and reporting checklists.
- [references/unitts-goat.md](references/unitts-goat.md) — UniTS-Gen
  installation/CLI, GOAT-TS syntax, multi-sample/prefilter strategy, limits.
- [references/background.md](references/background.md) — thermodynamics vs
  kinetics, energy span, flux branching, and why the protocol uses ΔE/IRC.
- [references/ai-role.md](references/ai-role.md) — what AI/ML actually
  contributes to TS finding, with evidence and limits.

The suite is self-contained: it does not require the author's local paths, the
older `ts-generation-workflow` skill, or a specific cluster account. Replace
the example cluster settings with your own.

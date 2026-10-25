# Printable checklists

Copy these into a lab notebook or a job log. Every unchecked box is a reason
not to trust the number yet.

## P0 — declaration sheet

```text
Project / date:
Elementary step (one sentence):
Bonds forming (atom indices, element symbols, starting distance):
Bonds breaking (atom indices, element symbols, starting distance):
Charge / multiplicity:
Accessible spin states (d^n), all enumerated?   [ ] yes [ ] no
Broken-symmetry solution needed?                [ ] yes [ ] no
Reactant/encounter geometry source:
Product geometry source:
Conformer search done?                          [ ] yes [ ] no  (method: ____)
Reference state for the barrier:
Method / basis / dispersion / solvent / grid:
ORCA or Gaussian version:
```

## P1 — barrier scan

- [ ] Reaction coordinate is the chemically correct one (metal-centered
      substitution vs C–C distance vs H transfer).
- [ ] 4–8 points span the barrier region, not just the product side.
- [ ] Each point is an independent `! Opt` (ORCA) or `opt=modredundant`
      (Gaussian), starting from the previous **converged** geometry.
- [ ] Only `CONVERGED` points are used.
- [ ] The highest converged point is recorded as the L1 guess.
- [ ] The curve is checked for a shallow double bump (conformer/channel change).
- [ ] R/S or substrate comparisons use identical constraints, solvent, spin,
      and MaxIter; the inputs were diffed.

## P2 — constrained reference + Hessian

- [ ] The constrained point converged.
- [ ] `Freq` ran at the constrained point; `.hess` / `.chk` exists.
- [ ] The lowest modes involve the declared reactive pair.
- [ ] If the lowest mode is a spectator, the constraint/conformer was changed
      **before** L2.
- [ ] The constrained point is labeled "reference point", not TS.

## P3 — OptTS

- [ ] `InHess Read` (ORCA) or `readfc` (Gaussian) is present.
- [ ] The Hessian comes from the same geometry/level as the start.
- [ ] MaxIter is appropriate for the flexibility (120–160 soft, 220–300 rigid).
- [ ] If displaced: the displacement was scaled by the target bond-length
      change (~0.05 Å), not the maximum atomic displacement.
- [ ] The energy trace is not monotonically decreasing into a basin.
- [ ] If it did not converge, a separate `Freq` was run on the last geometry.

## P4 — validation

- [ ] Exactly one significant imaginary frequency.
- [ ] The imaginary mode moves the declared reactive pair.
- [ ] Mode localization is not `<15%` (delocalized soft mode).
- [ ] Key bond lengths match the intended step.
- [ ] IRC both directions (or downhill relaxation) reaches the intended
      reactant and product.
- [ ] ⟨S²⟩ and Mulliken spin are consistent with the geometry job.
- [ ] The single point inherits the BS settings (`UKS/MORead/FlipSpin/FinalMs`).
- [ ] Optimization/frequency/IRC are at the same level.
- [ ] ΔE used for minimum-vs-saddle; ΔG only within the same type.
- [ ] SMD/CDS convention is consistent across all species.

## HPC submission

- [ ] Submitted to an allowed partition (not a default that your account
      cannot use).
- [ ] `%nprocshared` = `#SBATCH -n`.
- [ ] `%mem` is realistic for the node (e.g. 48GB for 24 cores).
- [ ] `GAUSS_SCRDIR` is node-local (`/tmp`), not `/home`.
- [ ] `g16root`/`source` is set, or `g16` is on `PATH`.
- [ ] The job directory is the intended one; job names are self-describing.
- [ ] Upstream `.chk`/`.hess` exists before the dependent step starts.
- [ ] The job is not running on the login node.
- [ ] `Normal termination` **and** the expected deliverable (frequencies,
      convergence table) are present.

## Reporting

- [ ] Method, basis, dispersion, solvent, grid, program version.
- [ ] Charge, multiplicity, spin surface, ⟨S²⟩.
- [ ] Imaginary frequency and moving atoms.
- [ ] Key bond lengths at the TS.
- [ ] IRC/relaxation endpoints and connectivity.
- [ ] Reference state for the barrier.
- [ ] Conformer search + Boltzmann averaging if selectivity is claimed.
- [ ] Method uncertainty stated; no 0.1 kcal claims inside a 2–5 kcal error.
- [ ] Constrained reference points labeled as such.
- [ ] AI/automated candidates reported with tool, version, candidate counts,
      pass rates, and validation level.

## Red flags — stop and re-check

- The "reactant" has an unplanned new bond.
- The lowest mode at P2 is a ligand torsion.
- The OptTS energy decreases for >50 cycles.
- Two displaced runs return to the same geometry/energy.
- The imaginary mode is delocalized (`<15%` on the transferring atom).
- `SADDLE-1` but a key bond length is wrong.
- The IRC endpoint is not the intended product.
- The single-point ⟨S²⟩ differs from the geometry job.
- The R/S number comes from one conformer.
- The AI-generated candidate was never prefilted or IRC-checked.

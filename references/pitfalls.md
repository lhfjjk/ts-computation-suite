# Pitfalls catalogue

These are the failure modes that cost real time in the source projects. Each
row is written as **symptom → root cause → fix**, so it can be used as a
checklist rather than read linearly.

## A. Structure and mechanism

| Symptom | Root cause | Fix |
|---|---|---|
| The "reactant" is actually a rearrangement product | The guess was built by hand or by a model and never checked | Build a bond graph with `1.3 × (r_i + r_j)`; check every intended bond and every unplanned bond. Three real cases: 1,2-acyl migration, ketone addition + β-H elimination, spontaneous O–Cα closure (1.442 Å) |
| A script calibrated on one geometry gives nonsense on another | Different atom count/composition (91 vs 98 atoms; missing ethyl) | Compare element counts and atom indices before reusing index-based scripts |
| RDKit assigns the opposite R/S label at a metal-bound carbon | CIP priorities for metals are not handled the way the project expects | Assign the metal-carbon configuration manually or with the project convention; do not trust RDKit CIP for metal-bound carbon |
| The 2D figure and the 3D geometry disagree | Figure built from a different build file than the calculation | Draw from the actual optimized geometry; check atom count and connectivity |

## B. Spin and electronic state

| Symptom | Root cause | Fix |
|---|---|---|
| BDFE/barrier systematically off by ~15 kcal | Reference state was run as a singlet; the ground state is a triplet | Enumerate all accessible dⁿ states; read `HFTyp` in the output, not the input; optimize each state |
| A BS single point gives a closed-shell singlet | Missing `UKS + MORead + moinp + FlipSpin + FinalMs`; ORCA defaults to RKS | Repeat the full BS setup in every single point; check ⟨S²⟩ matches the geometry job |
| `FlipSpin` flips the wrong atom | `FlipSpin` is 1-based; constraints and Mulliken rows are 0-based | Convert explicitly; print the element symbol to confirm |
| Two calculations at the same geometry give different ⟨S²⟩ | Different SCF solutions (BS vs closed-shell vs triplet) | Report ⟨S²⟩; do not compare energies across different solution types |
| "No TS on this spin surface" becomes a rule | A failed OptTS was treated as evidence of absence | Re-run that surface with the correct coordinate and Hessian; report a negative result with settings (coordinate, Hessian source, MaxIter, endpoint connectivity) |

## C. Energy and thermodynamics

| Symptom | Root cause | Fix |
|---|---|---|
| Reaction energy shifts by 5–11 kcal between scripts | SMD/CDS convention mixed: `E_FINAL` vs `E_CDS` | Pick one: `E_CDS + (G − E_CDS)` or `E_FINAL + (G − E_FINAL)`; never mix |
| The first energy in the output is used | `re.search` instead of the last match | Always take the last `FINAL SINGLE POINT ENERGY` / `SCF Done` / `HF=` |
| ΔE and ΔG have opposite signs | Saddle's imaginary frequency is excluded from the partition function | Use ΔE for minimum-vs-saddle; use ΔG only within the same stationary-point type |
| A single point is lower-level than the optimization | Level hierarchy violated | Single-point level ≥ optimization level; same functional/basis for all compared species |
| ORCA G and E don't match the expected correction | ORCA's G uses the dispersion-containing `E_FINAL` | Verify with the `G-E(el)` line; derive the correction before using it |

## D. TS diagnosis

| Symptom | Root cause | Fix |
|---|---|---|
| `SADDLE-1` and the declared pair matches, but the chemistry is wrong | The mode is a different channel with a coincidental atom match | Read the key bond lengths; run IRC. Real case: −1363 cm⁻¹ with Ni–Cα 3.712 Å and Ni–O 1.971 Å = enolization, not β-H elimination |
| The eigenvector sign is read as forming/breaking | Eigenvector sign is arbitrary | Only read which pairs move most; determine direction from IRC/relaxation or reactant/product geometries |
| A soft mode looks like the reaction coordinate | Delocalized spectator mode | Compute the transferring atom's share of the displacement; `<15%` is suspicious (real TS ~40%) |
| 2–3 shallow imaginary modes around −15…−25 cm⁻¹ | Flat transition valley, real physics | Do not restart repeatedly; use the last geometry's Freq as a reference barrier and say so |
| Two displaced OptTS runs return to the same geometry | The point is a second-order saddle or a flat valley; the optimizer slides back | Stop after the second repeat. Accept with disclosure, rigidify the spectator (e.g. tBu → Me/Et), or IRC along the reactive mode |
| A second soft mode is dismissed as a grid artifact | It may be real | Recompute the frequency with a finer grid (`DefGrid3`); if it persists, it is real |
| `OptTS` runs 200–300 cycles and never converges | Missing `InHess Read`; the model Hessian points along the wrong mode | Run L1 constrained `Opt+Freq`, then L2 with `InHess Read / InHessName` |
| A displacement of 0.18 Å throws the geometry into the reactant basin | The displacement was scaled by maximum atomic displacement, not by the reactive bond change | Scale by the target bond-length change (~0.05 Å); check the reactive bond before running |
| `OptTS` finishes but there is no frequency block | ORCA skips frequencies when the optimization did not converge | Run `Freq` on the last geometry; do not just raise MaxIter |
| A generated candidate has the right formula but the wrong TS | The model was trained on a different chemistry/size | Prefilter with RMSD/connectivity; large BS organometallic systems are out of distribution |

## E. Scan and optimization

| Symptom | Root cause | Fix |
|---|---|---|
| A relaxed scan dies and the whole profile is lost | ORCA 6.1.1: any point hitting `MaxIter` breaks the archive/cleanup chain | Use per-point independent `! Opt` jobs (`ts_scan_pts.py`); one failure stays local |
| `%geom Scan` returns `Unknown error in GEOM block` | The `Scan` block is missing its own `end`, so `MaxIter` is swallowed | Give `Scan ... end` its own `end`; keep `MaxIter` outside |
| R and S profiles are not comparable | Constraint values differ (e.g. 2.13/2.39 vs 2.08/2.42) | Diff the two inputs; only stereochemistry/substituent may differ; constraints, solvent, spin, MaxIter must match |
| A failed scan point is used as a data point | The parser accepted non-converged steps | Only accept `CONVERGED` points; mark failed points and skip them |
| The profile is monotonic but the TS is missing | The scan never reached the barrier region, or the coordinate is wrong | Check the coordinate definition (metal-centered substitution vs C–C distance); extend the scan |
| Increasing `MaxIter` does not help | The starting geometry is in the wrong basin | Stop; return to P1/P2; a larger MaxIter only makes the wrong answer more expensive |
| MLIP energies or bond lengths are quoted in a paper | A learned potential is a screening tool, not evidence; on this project's Ni/Fe systems conformer/stereoisomer differences were off by ~2.6 kcal mol⁻¹ | Use the MLIP for repair, conformer folding, and screening only; all reported numbers from DFT |

## F. Process and monitoring

| Symptom | Root cause | Fix |
|---|---|---|
| Two ORCA jobs run at once | A waiter checked only `orca.exe` and slipped into the gap between chained jobs | Check for `orca.exe` **and** the project `run-*.sh`; require idle for 2 × 30 s |
| A kill command kills its own shell | The command line contained the script name it was matching | Match the exact script name and print PIDs first; do not put the literal name in the `echo` text |
| A Chinese-path grep returns nothing | Encoding/locale issue in the shell | Use ASCII paths or verify with Python; do not trust an empty grep |
| Orphan waiters start ORCA later | Killing the queue parent did not kill its children | Enumerate `python.exe` PIDs/parents; remove exact PIDs; check again after restart |
| A job looks stalled because CPU time went backwards | Module switch (e.g. `leanscf` → `scfgrad`) creates a new process counter | Judge health by `.out/.tmp` mtime and process names, not by a negative CPU delta |
| A long job's log is empty after a kill | Output was piped to `tail`; the pipe hid the exit code and buffered the log | Write stdout to a file; do not pipe long jobs |

## G. HPC / Gaussian

| Symptom | Root cause | Fix |
|---|---|---|
| `allocation failure: User's group not permitted` | Submitted to the `theory` partition | Use `community` (or the partition your account is allowed) |
| The job waits for days | `small` partition, or a 32-core request | Use `community` with 24 cores |
| The job is killed for memory | `%mem` over-request on a shared node; the cluster may not enforce it | Set `%mem` to 48GB for 24 cores / 32GB for 16; keep `GAUSS_SCRDIR` on `/tmp` |
| `Normal termination` but no frequencies | `freq` missing, or the job was a single point | Check the route line; add `freq` |
| `.chk` read error in the TS step | `%oldchk` points to an overwritten/missing checkpoint | Use step-specific checkpoint names; rerun the upstream step |
| `readfc` does not work | The checkpoint has no analytic Hessian | Use `opt=(ts,calcfc,noeigen) freq` (recompute the Hessian) |
| Gaussian runs on the login node | Direct `g16` invocation | Always submit through `sbatch` |
| A 32-core job waits ~21 h | Single-node 32-core requests need a whole idle node | Use 24 cores; the per-job speedup does not repay the queue |
| The node runs out of memory | `%mem` is over-requested and the cluster does not enforce it | 48GB for 24 cores / 32GB for 16; keep scratch on `/tmp` |
| The job "succeeds" but is a single point | `Opt` was omitted from the route line | Check the route; a `Normal termination` alone is not evidence of an optimization |
| `readfc` fails or reads no Hessian | The checkpoint has no analytic Hessian, or `%oldchk` is wrong | Copy the L1 chk to the L2 chk name, or use `opt=(ts,calcfc,noeigen)` |
| The TS step reads an overwritten checkpoint | A later job reused the same `%chk` name | Use step-specific checkpoint names; never overwrite a checkpoint still needed as `%oldchk` |

## H. Reporting and interpretation

| Symptom | Root cause | Fix |
|---|---|---|
| An R/S energy difference is reported as stereoselectivity | Single conformer; the difference is conformer-dominated (RMSD 1.62 Å) | Do a conformer search, optimize the low-energy conformers, Boltzmann-average; otherwise report it as a single-conformer result |
| A constrained point is reported as a TS | It is a reference point, not a stationary point | Label it "constrained reference point"; only a one-imaginary-mode stationary point is a TS |
| A 0.1 kcal difference is treated as meaningful | Method error is 2–5 kcal/mol | Report the method uncertainty; use differences only when they exceed the error |
| An open-shell energy is reported without ⟨S²⟩ | Spin contamination/solution identity is hidden | Report ⟨S²⟩ and the spin surface with every open-shell number |
| An AI-generated TS is reported without validation | The model's output was treated as a result | Report the tool/version, number of candidates, prefilter pass rate, and the DFT/IRC validation |

## I. Candidate generation (UniTS-Gen / GOAT)

| Symptom | Root cause | Fix |
|---|---|---|
| The generated candidate points at the wrong atom | UniTS uses 0-based `reactive_atom_idx`; this suite's `--reactive` is 1-based | Convert explicitly; print the element symbols to confirm |
| A large BS organometallic candidate looks plausible but is not the target TS | The model is out of distribution for 80+ atom, flexible, broken-symmetry systems | Prefilter with RMSD/connectivity; prefer a manual guess + L1/L2 for these systems |
| A `%GOAT` parameter is rejected as unknown | Binary strings can be internal function names, not user parameters | Use the verified list in `unitts-goat.md`; probe the parser before relying on a new keyword |
| GOAT is very slow or does not finish | 70+ flexible atoms; exploration cost | Use it for small/medium systems or channel discovery; use a manual guess for production |
| A generated candidate is treated as a validated TS | Generation and validation are different steps | Run prefilter → L1 → L2 → Freq/IRC; report the candidate counts and pass rate |

---

## The one-page "everything looks fine but the answer is wrong" check

1. Did you read the **last** energy, not the first?
2. Is the single-point level ≥ the optimization level, with the same
   functional/basis/dispersion/solvent on both sides?
3. Did you compare ΔE for a minimum-vs-saddle, not ΔG?
4. Does the BS single point carry the full `UKS/MORead/FlipSpin/FinalMs` setup?
5. Does the imaginary mode touch the atoms you think it does, and is it
   localized?
6. Do the key bond lengths match the reaction you intended?
7. Does IRC connect the reactant and product you intended?
8. Is the R/S difference conformationally averaged?
9. Is the TS on every accessible spin surface, or did you default to one?
10. Did you label a constrained reference point as a TS?

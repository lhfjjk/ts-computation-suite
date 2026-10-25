# Gaussian / Slurm HPC execution

This reference covers the production route: Gaussian 16 on a Slurm cluster.
The protocol is the same as the local ORCA route; only the engine and the
submission layer change.

## 1. When to move to HPC

Move a TS job to the cluster when one of these is true:

- the system is too large for the workstation memory or wall-clock budget;
- you need a Gaussian-only method (e.g. a functional/basis combination not
  available in ORCA);
- you need many conformers/spin surfaces and want them to run in parallel
  across separate jobs;
- a local job would block the workstation for more than a day.

Do **not** move a job to HPC before the cheap gates: connectivity, conformer,
reactive distance, spin surface, and a rough barrier region. Queue time is
more expensive than local pre-screening.

## 2. Portable Slurm template

```bash
#!/bin/bash
#SBATCH -J <name>
#SBATCH -N 1
#SBATCH -n 24                 # must match %nprocshared
#SBATCH -c 1
#SBATCH -p <partition>
#SBATCH -o slurm-%j.out
#SBATCH --time=48:00:00

export PATH=/bin:/usr/bin:/usr/local/bin:/sbin:/usr/sbin:/usr/local/sbin
export LD_LIBRARY_PATH=/lib64:/usr/lib64:/usr/local/lib64:/usr/lib64/atlas
export g16root='<g16root>'
source "${g16root}/g16/bsd/g16.profile"
export GAUSS_SCRDIR=/tmp      # node-local scratch, never /home

date
g16 01_const_optfreq.gjf
g16 02_ts_optfreq.gjf
date
```

- `%nprocshared` in the `.gjf` and `#SBATCH -n` in the `.sbatch` must match.
- Set `%mem` to a realistic value. Many clusters do not enforce memory, so an
  over-request can OOM a shared node without being rejected.
- `GAUSS_SCRDIR` must be node-local (`/tmp`); pointing it at `/home` causes
  slow I/O and quota problems.

`gauss_ts_prep.py --slurm` writes this file with the placeholders filled.

## 3. Concrete example: FJIRSM gibbs cluster

The table below is a worked example from one account. Replace the placeholders
with your cluster's values; none of these values are required by the scripts.

| Item | Example value |
|---|---|
| SSH | `ssh <alias>` (host `<host>`, login node `<login>`, user `<user>`) |
| Scheduler | Slurm |
| Preferred partition | **`community`** (a general partition; check your site) |
| Avoid | a partition your account cannot use; a small partition with a long queue; a `debug` partition with too few cores |
| Cores | **24** (16→24 is ~48% faster at +1.2% core-hours; 32 cores often queues much longer) |
| Memory | `%mem=48GB` for 24 cores, `%mem=32GB` for 16 cores |
| Submit root | `<cluster-workdir>/<name>-<MMDD>/` (example: `/home/<user>/...`) |
| Gaussian | `<g16root>` (the site's Gaussian 16 installation) |
| Scratch | `GAUSS_SCRDIR=/tmp` |

The account's `qjob.sh` helper pushes a local directory to the TS root,
submits, and optionally watches/retrieves:

```bash
./qjob.sh <local-dir> <sbatch-file> [--watch] [--name <system>]
```

If you write your own `sbatch`, check that `-p community` is present; the
cluster's default sample script points at `theory`, which this account cannot
use.

## 4. Gaussian TS ladder

`gauss_ts_prep.py` generates the four jobs. The important route lines are:

### L1 — constrained reference + Hessian

```text
%mem=48GB
%nprocshared=24
%chk=01_const_optfreq.chk
#p UB3LYP/6-31G(d) opt=(calcfc,modredundant,maxcycles=220) freq em=gd3bj scrf=(smd,solvent=toluene) int=ultrafine

TS step 1: constrained Opt+Freq (reference point, not a TS)

0 2
<coordinates>

B 1 2 1.500 F
B 3 4 2.100 F

```

`B i j d F` freezes the distance at `d`. The `.chk` now contains the
analytic Hessian.

### L2 — OptTS reading the L1 Hessian

```text
%mem=48GB
%nprocshared=24
%oldchk=01_const_optfreq.chk
%chk=02_ts_optfreq.chk
#p UB3LYP/6-31G(d) opt=(ts,readfc,noeigen,maxcycles=220) freq em=gd3bj scrf=(smd,solvent=toluene) int=ultrafine geom=check guess=read

TS step 2: OptTS+Freq, Hessian from step 1

0 2

```

`readfc` is the Gaussian analogue of ORCA's `InHess Read`. If the checkpoint
does not contain a usable Hessian, use `calcfc` instead (recomputes it; slower
but robust):

```text
#p UB3LYP/6-31G(d) opt=(ts,calcfc,noeigen,maxcycles=220) freq ...
```

### L3 — IRC

```text
%oldchk=02_ts_optfreq.chk
%chk=03_irc.chk
#p UB3LYP/6-31G(d) irc=(readfc,maxpoints=50,stepsize=10) em=gd3bj scrf=(smd,solvent=toluene) int=ultrafine geom=check guess=read

IRC both directions from the TS geometry

0 2

```

### L4 — higher-level single point

```text
%oldchk=02_ts_optfreq.chk
%chk=04_sp.chk
#p UM06-2X/6-311+G(d,p) em=gd3bj scrf=(smd,solvent=toluene) int=ultrafine geom=check guess=read

Higher-level single point at the TS geometry

0 2

```

## 5. Checkpoint semantics

- `%chk=<file>` is the checkpoint this job writes.
- `%oldchk=<file>` is the previous job's checkpoint to read.
- `geom=check` reads the geometry from the checkpoint; `guess=read` reads the
  wavefunction.
- `readfc` reads the force constants (Hessian) from the checkpoint.
- Keep the `.chk` files in the same directory as the `.gjf`; `run_chain.sh`
  checks that the required upstream file exists before starting the next job.
- If a job is restarted, do not overwrite a `.chk` that is still needed as
  `%oldchk`. Copy it to a step-specific name first.

## 6. Output parsing

Use the last match, not the first:

```python
import re
text = open("02_ts_optfreq.log", errors="replace").read()

freqs = []
for line in text.splitlines():
    if "Frequencies --" in line:
        freqs += [float(x) for x in re.findall(r"-?\d+\.\d+",
                                                line.split("Frequencies --", 1)[1])]
imag = [v for v in freqs if v < 0]

e_scf = float(re.findall(r"SCF Done:\s+E\([^)]+\)\s+=\s+(-?\d+\.\d+)", text)[-1])
g     = float(re.findall(
    r"Sum of electronic and thermal Free Energies=\s+(-?\d+\.\d+)", text)[-1])
s2    = float(re.findall(r"<S\*\*2>=\s+(-?\d+\.\d+)", text)[-1])

normal = "Normal termination" in text
stationary = ("Stationary point found" in text
              or "Optimization completed" in text)
```

`gauss_verify.py` implements these checks and returns a shell-friendly exit
code. It cannot assign the imaginary mode to a bond; for that use
`freq=hpmodes` and inspect the mode in GaussView, or compare the key bond
lengths at the two IRC endpoints.

## 7. Queue and resource rules

- **24 cores** is the measured sweet spot on `community`; 32 cores often
  queues for ~21 h, which outweighs the per-job speedup.
- `#SBATCH -n` = `%nprocshared`; `-N 1` for a single node.
- `%mem` is not enforced by every cluster; set it to what the node can
  actually provide (48GB for 24 cores is a safe default; 32GB for 16).
- Do not run Gaussian on the login node. Everything through `sbatch`.
- For a chain of jobs, submit one sbatch that runs `run_chain.sh`, or submit
  each step with a dependency (`--dependency=afterok:<jobid>`). Do not start
  step 2 before step 1's `.chk` exists.
- If the queue is long, split the work into independent conformer/spin jobs
  and submit them together; do not run one giant job that blocks everything.

## 8. ORCA on Slurm

If the cluster has ORCA, the same protocol applies; only the submission layer
changes:

1. Generate the inputs locally with `ts_orca_prep.py` or `ts_ladder.py`.
2. Wrap each step with `assets/templates/slurm_orca.sbatch.tmpl` (set
   `ORCA_EXE`, `INPUT`, `OUTPUT`).
3. Submit steps with `--dependency=afterok:<jobid>` so step 2 does not start
   before step 1's `.hess` exists.
4. Do **not** use `ts_ladder_run.py` as the HPC driver; it waits on a desktop
   ORCA process and is intended for a workstation.

On the FJIRSM cluster, `/usr/bin/orca` may be a screen reader rather than the
quantum-chemistry binary. Load the real module (`module avail orca`), then
check `orca --version` before submitting. If there is no ORCA module, use the
Gaussian route instead of guessing a binary path.

## 9. Common HPC failures

| Symptom | Cause | Fix |
|---|---|---|
| `allocation failure: User's group not permitted` | Submitted to `theory` | Use `community` |
| Job waits for days | `small` partition, or 32-core request | Use `community` + 24 cores |
| `Error termination` in `l301`/`l502` | Basis/ECP mismatch, memory, or SCF | Check the route and `%mem`; try `calcfc`/`int=ultrafine`; inspect the last lines |
| `Normal termination` but no frequencies | `freq` missing, or the job was a single point | Add `freq`; check the route line |
| Geometry optimization converged but TS has 2 imaginary modes | Second-order saddle / flat valley | Follow the stop rules; do not just resubmit |
| `.chk` read error | `%oldchk` points to a missing/overwritten file | Keep step-specific checkpoints; rerun the upstream step |
| `GAUSS_SCRDIR` quota/I/O error | Scratch on `/home` | Set `GAUSS_SCRDIR=/tmp` |

For the chemistry-level failures (wrong channel, spin state, energy
convention), see [pitfalls.md](pitfalls.md).

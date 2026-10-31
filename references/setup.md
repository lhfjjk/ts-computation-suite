# Setup for a new user

This skill is designed to be copied to another machine and used without the
author's local paths, accounts, or auxiliary skills.

## 1. Requirements

| Component | Required? | Notes |
|---|---|---|
| Python | yes | 3.9+; the scripts use only the standard library |
| ORCA 6.x | for the ORCA route | `orca` on `PATH`, or set `TS_ORCA_EXE` |
| MLIP environment | optional | geometry repair, conformer folding, pre-screen; `torch` + `orb-models` (or `mace-torch`), pointed to by `MLIP_PYTHON` |
| Gaussian 16 | for the Gaussian route | `g16` on `PATH`, or set `TS_G16_EXE` |
| Slurm | for the HPC route | `sbatch`/`squeue`; not needed for local runs |
| SSH | optional | only for a remote cluster; `TS_SSH_ALIAS` |
| NumPy | optional | only some helper paths use it; not required |

No `pip install` is needed for the core scripts. Do not install the skill into
a project directory; put it in the agent's skill directory.

## 2. Install the skill

Copy the whole `ts-computation-suite/` folder into the skill directory used by
your agent, for example:

```bash
cp -r ts-computation-suite ~/.codex/skills/
# or, for ZCode:
cp -r ts-computation-suite ~/.zcode/skills/
```

Then start a new session and invoke it by name (`$ts-computation-suite`) or by
asking for a TS search. If your agent does not auto-discover new folders, point
it at `SKILL.md` directly.

## 3. Configure the environment

The scripts read these variables; all are optional if the tool is on `PATH`.

```bash
export TS_ORCA_EXE=/path/to/orca          # or ORCA_EXE
export ORCA_DIR=/path/to/orca-directory   # used by generated run_chain.sh
export MLIP_PYTHON=/path/to/mlip-env/bin/python   # torch + orb-models
export TS_G16_EXE=g16
export TS_SBATCH=sbatch
export TS_SSH_ALIAS=mycluster             # an entry in ~/.ssh/config
```

`ORCA_DIR` is used inside generated **bash** runners. On Windows Git Bash use a
POSIX-style path (`/c/ORCA_6.1.1`), not `C:\ORCA_6.1.1`. `TS_ORCA_EXE` can be a
native Windows path because the Python runner passes it as an argument.

For a cluster, also decide:

- partition/queue name;
- 24 cores vs 32 (24 is usually the better wall-clock/core-hour trade-off);
- `%mem` per core;
- a node-local scratch directory (`GAUSS_SCRDIR=/tmp`, never `/home`);
- the submit directory and job-name convention.

Put these in your shell profile or in the command line; the skill does not
write to your shell configuration.

## 4. Run the doctor

```bash
python scripts/doctor.py
python scripts/doctor.py --json
python scripts/doctor.py --strict    # fail if neither ORCA nor Gaussian exists
```

It checks Python, compiles every script, locates ORCA/Gaussian/Slurm/SSH,
verifies that `MLIP_PYTHON` (when set) can import `torch` and `orb_models`, and
runs four self-tests:

1. generate a per-point ORCA scan;
2. generate a Gaussian L1/L2 ladder;
3. parse a synthetic one-imaginary-frequency Gaussian log (expect `SADDLE-1`);
4. parse a synthetic two-imaginary-frequency log (expect `ORDER-2-SADDLE`).

No quantum-chemistry job is run. A passing doctor means the package is intact;
it does not mean ORCA/Gaussian is configured correctly for your system.

## 5. Five-minute smoke test

Before a production run, do one small calculation with a known answer:

```bash
# 1) generate a 3-point scan (no ORCA run)
python scripts/ts_scan_pts.py --xyz water.xyz --pair 1,2 \
    --start 1.10 --end 0.90 --points 3 --charge 0 --mult 1 \
    --method HF --basis def2-SVP --disp none --no-rijcosx \
    --outdir smoke_scan

# 2) run it (requires ORCA)
python scripts/ts_scan_pts.py --summary smoke_scan
# or: bash smoke_scan/run_scan.sh

# 3) generate a Gaussian ladder (no Gaussian run)
python scripts/gauss_ts_prep.py --xyz water.xyz --name smoke \
    --charge 0 --mult 1 --method HF --basis STO-3G --disp none \
    --solvent none --constraints "1,2,1.000" --outdir smoke_gauss \
    --no-irc --no-sp

# 4) inspect the inputs, then run one step
#    bash smoke_gauss/run_chain.sh   (requires Gaussian)
```

For a real system, replace water with a small molecule whose TS you already
know, and check that the imaginary frequency and the key bond lengths match
the known result.

## 6. Directory map

```text
ts-computation-suite/
  SKILL.md                  entry point and protocol summary
  references/
    setup.md                this file
    ts-protocol.md          five-stage protocol + stop rules
    local-orca.md           ORCA execution + optional MLIP helper
    hpc-slurm.md            Gaussian/Slurm execution
    pitfalls.md             failure catalogue
    checklists.md           printable P0-P4/HPC/reporting checklists
    unitts-goat.md          candidate generation (UniTS-Gen / GOAT-TS)
    background.md           thermodynamics/kinetics and AI context
    ai-role.md              what AI contributes and what it does not
  scripts/
    doctor.py               environment check + self-tests
    ts_prefilter.py         candidate filter
    ts_scan_pts.py          per-point independent ORCA scan
    ts_ladder.py            ORCA L1/L2/L3/SP generator
    ts_ladder_run.py        ORCA ladder runner
    ts_verify.py            ORCA TS verdict
    ts_orca_prep.py         candidate xyz -> ORCA inputs
    gauss_ts_prep.py        Gaussian L1/L2/IRC/SP generator + sbatch
    gauss_verify.py         Gaussian TS verdict
    preflight_check.py      structure sanity gate
    ts_mlip_screen.py       optional MLIP pre-screen (P0.5)
    mlip_preopt.py          optional MLIP geometry repair
    mlip_tools.py           optional MLIP backend: sp/opt/scan/conf
  assets/templates/         ORCA/Gaussian/Slurm templates
```

## 7. Platform notes

- **Windows:** use Git Bash or PowerShell; avoid paths with spaces for
  `run_chain.sh`. The ladder runner detects ORCA with PowerShell on Windows and
  `pgrep` on Linux.
- **Linux/macOS:** the scripts are standard-library Python; make sure `orca` and
  `g16` are on `PATH` or set the variables above. The MLIP helper additionally
  needs an environment with `torch` + `orb-models`.
- **Chinese/non-ASCII paths:** ORCA, Gaussian, and some viewers can fail
  silently. Use an ASCII working directory for production; if you must use a
  non-ASCII path, verify every output file was written.
- **No cluster account is bundled.** The HPC reference uses one cluster as a
  worked example; replace partitions, paths, and modules with your own.

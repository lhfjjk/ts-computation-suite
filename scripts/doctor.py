#!/usr/bin/env python3
"""Check a TS-computation environment and run parser/generator self-tests.

This is the first command a new user should run. It does not run quantum
chemistry; it checks that the scripts are intact, that ORCA/Gaussian/xtb are
discoverable, and that the input generators and output parsers work.

Examples
--------
python scripts/doctor.py
python scripts/doctor.py --json
python scripts/doctor.py --skip-self-test
python scripts/doctor.py --strict          # fail if neither ORCA nor Gaussian

Exit code: 0 if all required checks and the self-test pass; 1 otherwise.
Missing HPC tools (sbatch/ssh) are warnings, not failures.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
SKILL = HERE.parent


def resolve(name, env_names):
    for env in env_names:
        val = os.environ.get(env)
        if val:
            return val, f"${env}"
    found = shutil.which(name)
    if found:
        return found, "PATH"
    return None, "not found"


def compile_all():
    bad = []
    for p in sorted(HERE.glob("*.py")):
        try:
            compile(p.read_text(encoding="utf-8"), str(p), "exec")
        except SyntaxError as exc:
            bad.append((p.name, str(exc)))
    return bad


def run(cmd, cwd=None):
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, env=env)


def self_test():
    checks = []
    with tempfile.TemporaryDirectory(prefix="ts_doctor_") as td:
        td = Path(td)
        water = td / "water.xyz"
        water.write_text(
            "3\nwater\nO 0.000000 0.000000 0.000000\n"
            "H 0.957200 0.000000 0.000000\n"
            "H -0.240000 0.927000 0.000000\n", encoding="utf-8")

        scan = td / "scan"
        r = run([sys.executable, str(HERE / "ts_scan_pts.py"),
                 "--xyz", str(water), "--pair", "1,2",
                 "--start", "1.10", "--end", "0.90", "--points", "3",
                 "--charge", "0", "--mult", "1", "--method", "HF",
                 "--basis", "def2-SVP", "--disp", "none", "--no-rijcosx",
                 "--outdir", str(scan)])
        checks.append(("ts_scan_pts.py generation", r.returncode == 0
                       and (scan / "p01" / "point.inp").exists(),
                       r.stderr.strip()[:200]))

        gauss = td / "gauss"
        r = run([sys.executable, str(HERE / "gauss_ts_prep.py"),
                 "--xyz", str(water), "--name", "selftest",
                 "--charge", "0", "--mult", "1", "--method", "HF",
                 "--basis", "STO-3G", "--disp", "none", "--solvent", "none",
                 "--constraints", "1,2,1.000", "--outdir", str(gauss),
                 "--no-irc", "--no-sp"])
        checks.append(("gauss_ts_prep.py generation", r.returncode == 0
                       and (gauss / "02_ts_optfreq.gjf").exists(),
                       r.stderr.strip()[:200]))

        one = td / "ts_one.log"
        one.write_text(
            " SCF Done:  E(RHF) =  -76.000000000     A.U. after    1 cycles\n"
            " Stationary point found.\n"
            " Frequencies --   -500.0000    10.0000    20.0000\n"
            " Sum of electronic and thermal Free Energies=  -75.900000\n"
            " <S**2>=  0.0000\n"
            " Normal termination of Gaussian 16\n", encoding="utf-8")
        r = run([sys.executable, str(HERE / "gauss_verify.py"),
                 "--out", str(one)])
        checks.append(("gauss_verify.py SADDLE-1", r.returncode == 0
                       and "SADDLE-1" in r.stdout, r.stderr.strip()[:200]))

        two = td / "ts_two.log"
        two.write_text(
            " SCF Done:  E(RHF) =  -76.000000000     A.U. after    1 cycles\n"
            " Stationary point found.\n"
            " Frequencies --   -500.0000   -100.0000    10.0000\n"
            " Sum of electronic and thermal Free Energies=  -75.900000\n"
            " <S**2>=  0.0000\n"
            " Normal termination of Gaussian 16\n", encoding="utf-8")
        r = run([sys.executable, str(HERE / "gauss_verify.py"),
                 "--out", str(two)])
        checks.append(("gauss_verify.py ORDER-2", r.returncode == 1
                       and "ORDER-2-SADDLE" in r.stdout, r.stderr.strip()[:200]))

    return checks


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--skip-self-test", action="store_true")
    ap.add_argument("--strict", action="store_true",
                    help="fail if neither ORCA nor Gaussian is available")
    args = ap.parse_args()

    tools = {
        "ORCA": resolve("orca", ["TS_ORCA_EXE", "ORCA_EXE"]),
        "Gaussian": resolve("g16", ["TS_G16_EXE"]),
        "xtb": resolve("xtb", ["TS_XTB_EXE"]),
        "sbatch": resolve("sbatch", ["TS_SBATCH"]),
        "ssh": resolve("ssh", ["TS_SSH"]),
    }
    py_ok = sys.version_info >= (3, 9)
    compile_bad = compile_all()
    tests = [] if args.skip_self_test else self_test()

    result = {
        "python": {
            "version": sys.version.split()[0],
            "ok": py_ok,
        },
        "scripts_compile": {
            "ok": not compile_bad,
            "failures": compile_bad,
        },
        "tools": {
            name: {"path": path, "source": source, "found": bool(path)}
            for name, (path, source) in tools.items()
        },
        "ssh_alias": os.environ.get("TS_SSH_ALIAS", ""),
        "self_test": [
            {"name": name, "ok": ok, "detail": detail}
            for name, ok, detail in tests
        ],
    }

    required_ok = py_ok and not compile_bad and all(ok for _, ok, _ in tests)
    if args.strict and not (tools["ORCA"][0] or tools["Gaussian"][0]):
        result["strict_failure"] = "neither ORCA nor Gaussian was found"
        required_ok = False

    if args.json:
        print(json.dumps(result, indent=2, ensure_ascii=False))
    else:
        print("TS computation suite: environment check")
        print(f"  python: {sys.version.split()[0]}  {'OK' if py_ok else 'NEED >=3.9'}")
        print(f"  scripts: {'OK' if not compile_bad else 'COMPILE FAILURES'}")
        for name, (path, source) in tools.items():
            print(f"  {name:9s}: {path or 'not found'}  ({source})")
        if result["ssh_alias"]:
            print(f"  ssh alias: {result['ssh_alias']}")
        if not args.skip_self_test:
            print("  self-test:")
            for name, ok, detail in tests:
                print(f"    [{'OK' if ok else 'FAIL'}] {name}"
                      + (f"  {detail}" if detail and not ok else ""))
        print("\nSet these if a tool is not found:")
        print("  TS_ORCA_EXE=/path/to/orca   ORCA_DIR=/path/to/orca-dir")
        print("  TS_G16_EXE=g16              (and gauss_ts_prep.py --g16root ...)")
        print("  TS_XTB_EXE=/path/to/xtb")
        print("  TS_SSH_ALIAS=<ssh-config-alias>")
        print("\nSee references/setup.md for the full setup guide.")

    return 0 if required_ok else 1


if __name__ == "__main__":
    sys.exit(main())

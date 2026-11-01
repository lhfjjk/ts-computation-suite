#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""check_docs.py -- documentation consistency gate for this skill package.

Catches the failure mode that v0.4.2 fixed: SKILL.md and references/ drifting
apart (a stage that exists in one file and not the other), script inventory
mismatches, broken relative links, and a version that disagrees with CHANGELOG.

Standard library only. Exit 0 = clean, 1 = problems.
Usage: python scripts/check_docs.py [--verbose]
"""
import argparse
import io
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SKILL = "SKILL.md"
PROTOCOL = os.path.join("references", "ts-protocol.md")
CHECKLISTS = os.path.join("references", "checklists.md")
CHANGELOG = "CHANGELOG.md"
CITATION = "CITATION.cff"

problems = []
notes = []


def read(rel):
    with io.open(os.path.join(ROOT, rel), encoding="utf-8", errors="replace") as fh:
        return fh.read()


def check_links():
    """Every relative markdown link must resolve."""
    bad = 0
    for base, dirs, files in os.walk(ROOT):
        dirs[:] = [d for d in dirs if d not in (".git", "__pycache__", ".github")]
        for name in files:
            if not name.endswith(".md"):
                continue
            path = os.path.join(base, name)
            text = io.open(path, encoding="utf-8", errors="replace").read()
            for m in re.finditer(r"\[[^\]]+\]\(([^)\s]+)\)", text):
                target = m.group(1)
                if target.startswith(("http://", "https://", "mailto:", "#")):
                    continue
                target = target.split("#", 1)[0]
                if not target:
                    continue
                resolved = os.path.join(base, target.replace("/", os.sep))
                if not os.path.exists(resolved):
                    problems.append(
                        f"broken link: {os.path.relpath(path, ROOT)} -> {m.group(1)}")
                    bad += 1
    notes.append(f"links: {bad} broken")


def check_script_inventory():
    """Every scripts/*.py must appear in the SKILL.md scripts table."""
    skill = read(SKILL)
    script_dir = os.path.join(ROOT, "scripts")
    scripts = sorted(f for f in os.listdir(script_dir) if f.endswith(".py"))
    missing = [s for s in scripts if f"scripts/{s}" not in skill]
    for s in missing:
        problems.append(f"script not listed in SKILL.md: scripts/{s}")
    notes.append(f"scripts: {len(scripts)} found, {len(missing)} unlisted")


def check_reference_inventory():
    """Every references/*.md must be linked from SKILL.md."""
    skill = read(SKILL)
    ref_dir = os.path.join(ROOT, "references")
    refs = sorted(f for f in os.listdir(ref_dir) if f.endswith(".md"))
    missing = [r for r in refs if f"references/{r}" not in skill]
    for r in missing:
        problems.append(f"reference not linked from SKILL.md: references/{r}")
    notes.append(f"references: {len(refs)} found, {len(missing)} unlinked")


def check_stage_coverage():
    """The protocol and the checklists must carry every stage SKILL.md names."""
    skill = read(SKILL)
    stages = sorted(set(re.findall(r"\bP[0-9](?:\.[0-9])?\b", skill)))
    for rel in (PROTOCOL, CHECKLISTS):
        text = read(rel)
        missing = [s for s in stages if not re.search(r"\b" + re.escape(s) + r"\b", text)]
        for s in missing:
            problems.append(f"{rel} does not mention stage {s} (SKILL.md does)")
    notes.append(f"stages named in SKILL.md: {', '.join(stages)}")


def check_version():
    """CITATION version == top CHANGELOG heading."""
    cit = read(CITATION)
    m = re.search(r"^version:\s*(\S+)", cit, re.M)
    if not m:
        problems.append("CITATION.cff has no version field")
        return
    version = m.group(1)
    changelog = read(CHANGELOG)
    top = re.search(r"^##\s*\[([0-9]+\.[0-9]+\.[0-9]+)\]", changelog, re.M)
    if not top:
        problems.append("CHANGELOG.md has no versioned section")
        return
    if top.group(1) != version:
        problems.append(
            f"version mismatch: CITATION.cff {version} vs CHANGELOG {top.group(1)}")
    if f"[{version}]" not in changelog:
        problems.append(f"CHANGELOG has no entry for {version}")
    notes.append(f"version: {version}")


def check_keyword_evidence():
    """GOAT-TS may only appear where the text flags it as parser-only."""
    skill = read(SKILL)
    if "GOAT-TS" in skill and "parser" not in skill:
        problems.append("SKILL.md mentions GOAT-TS without the parser-only caveat")
    ref = read(os.path.join("references", "unitts-goat.md"))
    if "GOAT-REACT" not in ref:
        problems.append("references/unitts-goat.md does not mention the documented GOAT-REACT")
    notes.append("keywords: GOAT evidence note present")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    for fn in (check_links, check_script_inventory, check_reference_inventory,
               check_stage_coverage, check_version, check_keyword_evidence):
        before = len(problems)
        fn()
        if args.verbose and len(problems) == before:
            print(f"  [OK] {fn.__name__}")

    if args.verbose:
        for n in notes:
            print(f"  {n}")
    if problems:
        print(f"documentation check: {len(problems)} problem(s)")
        for p in problems:
            print(f"  [FAIL] {p}")
        return 1
    print(f"documentation check: OK ({len(notes)} checks)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

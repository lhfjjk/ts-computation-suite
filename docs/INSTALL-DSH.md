# Installing this skill in DeepSeek Harness (DSH)

DSH discovers skills by scanning `~/.dsh/skills/` for a directory containing
`SKILL.md`. There is no registry file and no restart: a new directory appears in
the skill catalog on the next turn, for **every** Bot and session on the machine.

## 1. Install

Keep one copy in a git checkout and link it in, so `git pull` updates the skill:

```powershell
# Windows (junction - no admin rights needed)
$repo = "$env:USERPROFILE\.zcode\skills\ts-computation-suite"
$link = "$env:USERPROFILE\.dsh\skills\ts-computation-suite"
New-Item -ItemType Junction -Path $link -Target $repo
```

```bash
# Linux / macOS
ln -s ~/code/ts-computation-suite ~/.dsh/skills/ts-computation-suite
```

For a frozen copy (no live updates), copy the tree instead of linking:

```powershell
Copy-Item -Recurse $repo $link
```

## 2. Verify

```powershell
Test-Path "$env:USERPROFILE\.dsh\skills\ts-computation-suite\SKILL.md"   # True
python "$env:USERPROFILE\.dsh\skills\ts-computation-suite\scripts\doctor.py"
```

The skill then appears in the session catalog by name (`ts-computation-suite`).
Ask for it explicitly the first time ("use the ts-computation-suite skill") so
the agent reads `SKILL.md` instead of improvising a TS workflow.

## 3. Update

```powershell
cd $repo; git pull        # the junction follows automatically
```

## 4. Remove

```powershell
Remove-Item "$env:USERPROFILE\.dsh\skills\ts-computation-suite" -Force
```

Removing a junction does not touch the target, but deleting or moving the
repository breaks the link, so keep the checkout in a stable location.

## 5. Optional: the MLIP stage

P0.5 needs `torch` + `orb-models` in a separate environment; point `MLIP_PYTHON`
at that interpreter. On Windows that is `$env:MLIP_PYTHON`, not `$MLIP_PYTHON`.
The rest of the suite runs on the standard library alone, so a missing MLIP
environment disables one optional stage, not the protocol.

#!/usr/bin/env python3
"""Workflow figure for the ChemRxiv manuscript.

Design constraints (2026-10-28 revision, after the 10-28 regression to the old 5.62 in / 6.1 pt spec)
    * canvas width == the printed width, so a nominal 7.6 pt line is a printed 7.6 pt line
      (FIG_W = 6.3 in; height keeps the original 7.2 : 4.6 aspect);
    * minimum printed size >= 7.2 pt: card title 9.0, card body 7.6, bottom band 8.4,
      top bracket 7.8, bottom italic note 7.2;
    * every card keeps its stage name and the meaning of its gate; wording is shortened,
      nothing is dropped (the MLIP pre-screen gate of P0 is part of the content);
    * no "Figure 1" label and no panel letters inside the artwork;
    * saved WITHOUT bbox_inches, so the canvas really is FIG_W x FIG_H inches.

Run from the submission directory:  python make_figure.py --check
Writes figures/workflow.pdf, figures/workflow.png (600 dpi) and
figures/_preview_6p3in.png (150 dpi, the printed width), then mirrors the pdf/png
byte-identically into ../figures (the copy manuscript.tex includes); with --check it
also verifies every text bbox against its card box, lists the nominal font sizes,
asserts that no required phrase was lost and prints size + SHA256 of every output.
"""
import argparse
import hashlib
import os
import shutil
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

matplotlib.rcParams["pdf.fonttype"] = 42   # embed TrueType, not Type 3 (publisher-safe)
matplotlib.rcParams["ps.fonttype"] = 42

FIG_W, FIG_H = 6.3, 4.03          # 7.2 : 4.6 aspect at the printed width
FS_TITLE, FS_BODY = 9.0, 7.6
FS_BAND, FS_BRACKET, FS_NOTE = 8.4, 7.8, 7.2

PALETTE = {
    "ai":   ("#DCEBFA", "#4A78B0"),
    "phys": ("#FCE8D5", "#C77A2E"),
    "dft":  ("#DFF0E0", "#4E8B57"),
    "val":  ("#E8E0F3", "#6B5B95"),
    "hum":  ("#ECECEC", "#666666"),
    "red":  "#B03A2E",
    "ink":  "#333333",
}

fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=300)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(-0.8, 100.8)
ax.set_ylim(-0.8, 64.8)
ax.axis("off")

CARD_PATCHES = []


def rounded(x, y, w, h, fc, ec, lw=1.0, radius=1.2, z=2):
    p = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.45,rounding_size={radius}",
        facecolor=fc, edgecolor=ec, linewidth=lw, zorder=z,
    )
    ax.add_patch(p)
    return p


def card(x, y, w, h, title, lines, key):
    fc, ec = PALETTE[key]
    p = rounded(x, y, w, h, fc, ec)
    t1 = ax.text(x + w / 2, y + h - 1.9, title, ha="center", va="top",
                 fontsize=FS_TITLE, fontweight="bold", color=PALETTE["ink"],
                 linespacing=1.12, zorder=3)
    t2 = ax.text(x + w / 2, y + h - 8.7, "\n".join(lines), ha="center", va="top",
                 fontsize=FS_BODY, color=PALETTE["ink"], linespacing=1.62, zorder=3)
    CARD_PATCHES.append((p, [t1, t2]))
    return p


def arrow(x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=7.5,
        linewidth=1.1, color="#555555", shrinkA=0, shrinkB=0, zorder=1,
    ))


# ---- top cards -----------------------------------------------------------
y, h, w, gap = 29.5, 28.0, 17.2, 3.2
xs = [1.0, 1.0 + w + gap, 1.0 + 2 * (w + gap), 1.0 + 3 * (w + gap),
      1.0 + 4 * (w + gap)]

card(xs[0], y, w, h,
     "1. AI candidate\ngeneration",
     ["UniTS-Gen / GOAT", "manual guess", "multi-sample,", "not one shot",
      "candidate =", "hypothesis"],
     "ai")
card(xs[1], y, w, h,
     "2. P0\npreflight",
     ["reactive pairs", "connectivity, clash", "RMSD / spin sanity",
      "MLIP pre-screen", "(where valid)", "most candidates", "die here"],
     "phys")
card(xs[2], y, w, h,
     "3. P1–P2\nreference",
     ["per-point scan", "constrained", "Opt + Freq", "Hessian for P3",
      "lowest-mode", "gate"],
     "dft")
card(xs[3], y, w, h,
     "4. P3\nOptTS",
     ["OptTS + InHess", "Read", "displace", "along mode",
      "0.05 Å bond change", "stop after 2 repeats"],
     "dft")
card(xs[4], y, w, h,
     "5. P4\nvalidation",
     ["one imaginary", "mode", "mode localization", "key bond lengths",
      "IRC / ⟨S²⟩ /", "energy"],
     "val")

for i in range(4):
    arrow(xs[i] + w + 0.15, y + h / 2, xs[i + 1] - 0.15, y + h / 2)

# ---- validation-first bracket -------------------------------------------
x0, x1 = xs[1] - 0.2, xs[4] + w + 0.2
yb = 60.6
ax.plot([x0, x1], [yb, yb], color=PALETTE["red"], lw=1.1, zorder=1)
for xb in (x0, x1):
    ax.plot([xb, xb], [yb, yb - 1.2], color=PALETTE["red"], lw=1.1, zorder=1)
ax.text((x0 + x1) / 2, yb + 0.9,
        "validation-first: cheap gates before expensive DFT",
        ha="center", va="bottom", fontsize=FS_BRACKET, color=PALETTE["red"],
        fontweight="bold")

# ---- bottom band ---------------------------------------------------------
yb2, hb = 12.0, 9.2
segs = [
    (1.0, 24.0, "AI proposes", "ai"),
    (25.0, 24.0, "physics filters", "phys"),
    (49.0, 30.0, "DFT validates", "dft"),
    (79.0, 20.8, "human decides", "hum"),
]
for x, ww, label, key in segs:
    fc, ec = PALETTE[key]
    rounded(x, yb2, ww, hb, fc, ec, lw=0.9, radius=0.9, z=2)
    ax.text(x + ww / 2, yb2 + hb / 2, label, ha="center", va="center",
            fontsize=FS_BAND, fontweight="bold", color=PALETTE["ink"], zorder=3)

ax.text(50, 6.2,
        "The generator can be replaced; the validation layer stays the same.",
        ha="center", va="bottom", fontsize=FS_NOTE, style="italic",
        color="#555555")


# ---- verification --------------------------------------------------------
REQUIRED = [
    # card titles
    "1. AI candidate generation", "2. P0 preflight", "3. P1–P2 reference",
    "4. P3 OptTS", "5. P4 validation",
    # card 1
    "UniTS-Gen / GOAT", "manual guess", "multi-sample,", "not one shot",
    "candidate = hypothesis",
    # card 2
    "reactive pairs", "connectivity, clash", "RMSD / spin sanity",
    "MLIP pre-screen", "(where valid)", "most candidates die here",
    # card 3
    "per-point scan", "constrained", "Opt + Freq", "Hessian for P3",
    "lowest-mode", "gate",
    # card 4 (wording must stay complete, breaking lines is allowed)
    "OptTS + InHess Read", "displace along mode", "0.05 Å bond change",
    "stop after 2 repeats",
    # card 5
    "one imaginary mode", "mode localization", "key bond lengths",
    "IRC / ⟨S²⟩ /", "energy",
    # band, bracket, note
    "AI proposes", "physics filters", "DFT validates", "human decides",
    "validation-first: cheap gates before expensive DFT",
    "The generator can be replaced; the validation layer stays the same.",
]


def _norm(s):
    return " ".join(s.split())


def missing_phrases():
    corpus = _norm(" ".join(t.get_text() for t in ax.texts))
    return [p for p in REQUIRED if _norm(p) not in corpus]


def _bbox_in_data(text):
    fig.canvas.draw()
    bb = text.get_window_extent(renderer=fig.canvas.get_renderer())
    inv = ax.transData.inverted()
    (x0d, y0d), (x1d, y1d) = inv.transform([[bb.x0, bb.y0], [bb.x1, bb.y1]])
    return x0d, y0d, x1d, y1d


def verify():
    problems = []
    fig.canvas.draw()
    inv = ax.transData.inverted()
    for patch, texts in CARD_PATCHES:
        pr = patch.get_window_extent(renderer=fig.canvas.get_renderer())
        (px0, py0), (px1, py1) = inv.transform([[pr.x0, pr.y0], [pr.x1, pr.y1]])
        bx0, bx1 = min(px0, px1), max(px0, px1)
        by0, by1 = min(py0, py1), max(py0, py1)
        for t in texts:
            x0d, y0d, x1d, y1d = _bbox_in_data(t)
            if x0d < bx0 - 0.05 or x1d > bx1 + 0.05 or y0d < by0 - 0.05 or y1d > by1 + 0.05:
                problems.append((t.get_text().split("\n")[0][:28],
                                 round(x1d - x0d, 1), round(bx1 - bx0, 1)))
    for t in ax.texts:
        x0d, y0d, x1d, y1d = _bbox_in_data(t)
        if x0d < ax.get_xlim()[0] or x1d > ax.get_xlim()[1] \
           or y0d < ax.get_ylim()[0] or y1d > ax.get_ylim()[1]:
            problems.append(("OUTSIDE CANVAS: " + t.get_text()[:28], 0, 0))
    fonts = sorted({round(t.get_fontsize(), 2) for t in ax.texts})
    return problems, fonts


OUT_SUBDIR = "figures"        # the copies that ship with the submission package
MIRROR_DIR = "../figures"     # what manuscript.tex actually includes; kept byte-identical
PNG_DPI = 600
PREVIEW_DPI = 150


def _sha(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest().upper()


def main():
    ap = argparse.ArgumentParser(description="Render the workflow figure at its printed size.")
    ap.add_argument("--check", action="store_true",
                    help="print the in-box verification, the nominal font sizes and the hashes")
    args = ap.parse_args()

    pdf = os.path.join(OUT_SUBDIR, "workflow.pdf")
    png = os.path.join(OUT_SUBDIR, "workflow.png")
    prev = os.path.join(OUT_SUBDIR, "_preview_6p3in.png")
    fig.savefig(pdf, format="pdf")
    fig.savefig(png, format="png", dpi=PNG_DPI)
    # a preview at exactly the printed width, so the printed font sizes can be judged
    fig.savefig(prev, format="png", dpi=PREVIEW_DPI)
    print("wrote %s, %s (%.2f x %.2f in, png %d dpi) and %s (%d dpi)"
          % (pdf, png, FIG_W, FIG_H, PNG_DPI, prev, PREVIEW_DPI))

    # manuscript.tex includes ../figures/workflow.pdf: keep the two copies byte-identical
    if os.path.isdir(MIRROR_DIR):
        for name in ("workflow.pdf", "workflow.png"):
            src = os.path.join(OUT_SUBDIR, name)
            dst = os.path.join(MIRROR_DIR, name)
            shutil.copyfile(src, dst)
            if _sha(src) != _sha(dst):
                print("MIRROR MISMATCH: %s" % dst)
                return 5
            print("mirrored %-13s -> %s (byte-identical)" % (name, MIRROR_DIR))
    else:
        print("mirror dir %r not found, skipped" % MIRROR_DIR)

    if args.check:
        problems, fonts = verify()
        missing = missing_phrases()
        print("canvas             : %.3f x %.3f in" % (FIG_W, FIG_H))
        print("nominal font sizes : %s  (min %.2f pt)" % (fonts, min(fonts)))
        print("text outside its box: %d" % len(problems))
        for p in problems:
            print("   line %-30s width %.1f units vs box %.1f units" % p)
        print("required phrases missing: %s" % (missing if missing else "none"))
        from PIL import Image
        for path in (pdf, png, prev, os.path.join(MIRROR_DIR, "workflow.pdf")):
            if not os.path.exists(path):
                continue
            extra = ""
            if path.endswith(".png"):
                im = Image.open(path)
                extra = "  %dx%d px" % im.size
            print("%-32s %8d B  %s%s" % (path.replace("\\", "/"), os.path.getsize(path),
                                         _sha(path), extra))
        return 0 if (not problems and not missing) else 4
    return 0


if __name__ == "__main__":
    sys.exit(main())

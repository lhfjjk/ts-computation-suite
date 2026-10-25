#!/usr/bin/env python3
"""Workflow figure for the ChemRxiv manuscript."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

FIG_W, FIG_H = 7.2, 4.6
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H), dpi=300)
ax.set_xlim(0, 100)
ax.set_ylim(0, 64)
ax.axis("off")

PALETTE = {
    "ai":   ("#DCEBFA", "#4A78B0"),
    "phys": ("#FCE8D5", "#C77A2E"),
    "dft":  ("#DFF0E0", "#4E8B57"),
    "val":  ("#E8E0F3", "#6B5B95"),
    "hum":  ("#ECECEC", "#666666"),
    "red":  "#B03A2E",
    "ink":  "#333333",
}


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
    rounded(x, y, w, h, fc, ec)
    ax.text(x + w / 2, y + h - 2.2, title, ha="center", va="top",
            fontsize=7.9, fontweight="bold", color=PALETTE["ink"],
            linespacing=1.15, zorder=3)
    ax.text(x + w / 2, y + h - 8.6, "\n".join(lines), ha="center", va="top",
            fontsize=6.1, color=PALETTE["ink"], linespacing=1.4, zorder=3)


def arrow(x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=7.5,
        linewidth=1.1, color="#555555", shrinkA=0, shrinkB=0, zorder=1,
    ))


# ---- top cards -----------------------------------------------------------
y, h, w, gap = 22.0, 32.0, 17.2, 3.2
xs = [1.0, 1.0 + w + gap, 1.0 + 2 * (w + gap), 1.0 + 3 * (w + gap),
      1.0 + 4 * (w + gap)]

card(xs[0], y, w, h,
     "1. AI candidate\ngeneration",
     ["UniTS-Gen / GOAT", "manual guess", "multi-sample,\nnot one shot",
      "candidate =\nhypothesis"],
     "ai")
card(xs[1], y, w, h,
     "2. Physics\nprefilter",
     ["reactive pairs", "connectivity /\nclash", "RMSD / spin\nsanity",
      "most candidates\ndie here"],
     "phys")
card(xs[2], y, w, h,
     "3. P1–P2\nreference",
     ["per-point\nindependent scan", "constrained\nOpt + Freq",
      "Hessian for\nnext step", "lowest-mode gate"],
     "dft")
card(xs[3], y, w, h,
     "4. P3\nOptTS",
     ["OptTS +\nInHess Read", "displace along\nmode",
      "0.05 Å target\nbond change", "stop after two\nrepeats"],
     "dft")
card(xs[4], y, w, h,
     "5. P4\nvalidation",
     ["one imaginary\nmode", "mode localization", "key bond lengths",
      "IRC / S² /\nenergy"],
     "val")

for i in range(4):
    arrow(xs[i] + w + 0.15, y + h / 2, xs[i + 1] - 0.15, y + h / 2)

# ---- validation-first bracket -------------------------------------------
x0, x1 = xs[1] - 0.2, xs[4] + w + 0.2
yb = 57.0
ax.plot([x0, x1], [yb, yb], color=PALETTE["red"], lw=1.1, zorder=1)
for xb in (x0, x1):
    ax.plot([xb, xb], [yb, yb - 1.2], color=PALETTE["red"], lw=1.1, zorder=1)
ax.text((x0 + x1) / 2, yb + 1.1,
        "validation-first: cheap gates before expensive DFT",
        ha="center", va="bottom", fontsize=7.2, color=PALETTE["red"],
        fontweight="bold")

# ---- bottom band ---------------------------------------------------------
yb2, hb = 6.0, 8.6
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
            fontsize=7.6, fontweight="bold", color=PALETTE["ink"], zorder=3)

ax.text(50, 1.8,
        "The generator can be replaced; the validation layer stays the same.",
        ha="center", va="bottom", fontsize=6.6, style="italic",
        color="#555555")

# ---- save ---------------------------------------------------------------
fig.savefig("figures/workflow.pdf", bbox_inches="tight", pad_inches=0.02)
fig.savefig("figures/workflow.png", bbox_inches="tight", pad_inches=0.02,
            dpi=300)
print("wrote figures/workflow.pdf and figures/workflow.png")

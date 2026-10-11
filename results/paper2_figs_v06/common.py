"""Shared style for the v0.7 figures. One palette across every figure:
orange = early-passage target, blue = young-donor target, gray = not called /
no-effect controls, black = the four validated factors and reference lines."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ORANGE, BLUE, BLACK = "#D55E00", "#0072B2", "#000000"
LIGHT, MID, DARK = "#C8C8C8", "#8C8C8C", "#4D4D4D"
FS = 8
W = 7.5
FOUR = ["CRA_E2F3", "CRA_EZH2", "CRI_STAT3", "CRI_ZFX"]
NAME = {"CRA_E2F3": "E2F3", "CRA_EZH2": "EZH2", "CRI_STAT3": "STAT3", "CRI_ZFX": "ZFX"}

plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": FS, "axes.titlesize": FS,
    "axes.labelsize": FS, "xtick.labelsize": FS - 0.5, "ytick.labelsize": FS - 0.5,
    "legend.fontsize": FS - 0.5, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "savefig.dpi": 300, "pdf.fonttype": 42, "legend.frameon": False,
})


def rows(path):
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def f(v):
    return float(v)


def tag(ax, letter, title, y=1.04):
    """Panel letter and title, anchored to the axes' top-left corner in points,
    so every panel's heading sits at the same place whatever its width."""
    lift = title.count("\n") * (FS + 0.5) * 1.25 * 0.75   # keep the letter level with the title's first line
    ax.annotate(letter, (0, y), xycoords="axes fraction", xytext=(-8, lift), textcoords="offset points",
                ha="right", va="bottom", fontweight="bold", fontsize=FS + 3)
    ax.annotate(title, (0, y), xycoords="axes fraction", ha="left", va="bottom",
                fontweight="bold", fontsize=FS + 0.5, linespacing=1.25)


def save(fig, out, name):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{name}.png"); fig.savefig(out / f"{name}.pdf")
    plt.close(fig)


def jitter(i, n, width=0.5):
    """Spread n points evenly across a band of the given width around 0."""
    return 0.0 if n == 1 else -width / 2 + width * i / (n - 1)

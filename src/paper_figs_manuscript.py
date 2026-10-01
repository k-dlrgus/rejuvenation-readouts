"""Manuscript Figures 1 and 3, 65-donor fibroblast atlas only.

Reads stored scores and stored intervals. Does not refit the ruler and does not
recompute gated statistics.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "paper_figs" / "manuscript"

BLUE = "#0072B2"
ORANGE = "#D55E00"
PINK = "#CC79A7"
BLACK = "#000000"

# Stored true-donor ruler Spearman (results/paper_figs/cxg_true_donor_results.csv).
ATLAS_RHO = 0.43338041757870716


def _font():
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 9,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    })


def _panel(fig, xy, letter):
    fig.text(xy[0], xy[1], letter, fontsize=12, fontweight="bold", ha="left", va="top")


def _open_spines(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=3.0, width=0.8)


def _signed(x, digits=3):
    """Fixed-point with an explicit sign and a unicode minus."""
    s = f"{x:+.{digits}f}"
    return s.replace("-", "\u2212")


def _load():
    transfer = pd.read_csv(ROOT / "results" / "fibro" / "stage1_transfer.csv")
    raw = transfer[(transfer.regime == "raw") & (transfer.method == "ridge")].copy()
    b1 = raw[raw.direction == "B1→C1"].iloc[0]
    c1 = raw[raw.direction == "C1→B1"].iloc[0]

    inst = pd.read_csv(ROOT / "results" / "md3" / "t2_instruments.csv")
    gse_ruler = inst[inst.tag == "gse226189_ruler_check"].iloc[0]
    gse_md = inst[inst.tag == "gse226189_md_check"].iloc[0]
    gse_s3 = inst[inst.tag == "GSE226189_S3"].iloc[0]

    donor_res = pd.read_csv(ROOT / "results" / "paper_figs" / "cxg_true_donor_results.csv")
    atlas = donor_res[(donor_res.analysis == "true_donor")].copy()
    atlas_rho = atlas[atlas.contrast == "ruler_spearman"].iloc[0]
    atlas_md = atlas[atlas.contrast == "ruler_minus_MD_AddModuleScore"].iloc[0]
    atlas_s3 = atlas[atlas.contrast == "ruler_minus_age_up_minus_age_down"].iloc[0]
    if abs(float(atlas_rho.rho) - ATLAS_RHO) > 1e-12:
        raise SystemExit("STOP: stored 65-donor rho is not 0.43338041757870716")
    if int(atlas_rho.n_units) != 65:
        raise SystemExit("STOP: atlas ruler row is not 65 donors")

    delta = pd.read_csv(ROOT / "results" / "md3" / "t2_delta_rho.csv")
    gse_d_md = delta[delta.tag == "gse226189_ruler_minus_MD"].iloc[0]
    gse_d_s3 = delta[delta.tag == "gse226189_ruler_minus_S3"].iloc[0]

    gse_pts = pd.read_csv(ROOT / "results" / "fibro2" / "ta_GSE226189_near_miss_bulk_scores.csv")
    if len(gse_pts) != 82:
        raise SystemExit(f"STOP: GSE226189 score table has {len(gse_pts)} rows, expected 82")

    don_pts = pd.read_csv(ROOT / "results" / "paper_figs" / "cxg_true_donor_scores.csv")
    if len(don_pts) != 65:
        raise SystemExit(f"STOP: true-donor score table has {len(don_pts)} rows, expected 65")
    # Points behind the stored rho. Rank correlation of the stored means; not a new fit.
    rho_pts = pd.Series(don_pts.age_score).corr(don_pts.age_years, method="spearman")
    if abs(float(rho_pts) - float(atlas_rho.rho)) > 1e-12:
        raise SystemExit("STOP: 65-donor points do not reproduce the stored rho")

    n_train = int(json.loads((ROOT / "results" / "fibro" / "freeze.json").read_text(encoding="utf-8"))["n_donors"])
    if n_train != 652:
        raise SystemExit(f"STOP: frozen ruler n_donors is {n_train}, expected 652")

    return dict(
        b1=b1, c1=c1, gse_ruler=gse_ruler, gse_md=gse_md, gse_s3=gse_s3,
        atlas_rho=atlas_rho, atlas_md=atlas_md, atlas_s3=atlas_s3,
        gse_d_md=gse_d_md, gse_d_s3=gse_d_s3,
        gse_pts=gse_pts, don_pts=don_pts, n_train=n_train,
    )


def figure1(d):
    fig = plt.figure(figsize=(7.0, 5.6))
    _panel(fig, (0.012, 0.985), "A")
    _panel(fig, (0.012, 0.455), "B")
    _panel(fig, (0.545, 0.455), "C")

    ax = fig.add_axes([0.300, 0.593, 0.440, 0.300])
    rows = [
        ("GTEx fibroblasts: site B1 \u2192 C1", f"{int(d['b1'].n_train)} \u2192 {int(d['b1'].n_test)}",
         float(d["b1"].rho), float(d["b1"].rho_ci_lo), float(d["b1"].rho_ci_hi)),
        ("GTEx fibroblasts: site C1 \u2192 B1", f"{int(d['c1'].n_train)} \u2192 {int(d['c1'].n_test)}",
         float(d["c1"].rho), float(d["c1"].rho_ci_lo), float(d["c1"].rho_ci_hi)),
        ("Frozen ruler \u2192 GSE226189 (bulk)", f"{d['n_train']} \u2192 {int(d['gse_ruler'].n_donors)}",
         float(d["gse_ruler"].rho), float(d["gse_ruler"].rho_ci_lo), float(d["gse_ruler"].rho_ci_hi)),
        ("Frozen ruler \u2192 fibroblast atlas", f"{d['n_train']} \u2192 65",
         float(d["atlas_rho"].rho), float(d["atlas_rho"].rho_ci_lo), float(d["atlas_rho"].rho_ci_hi)),
    ]
    ys = np.arange(len(rows))[::-1]
    for y, (label, _n, rho, lo, hi) in zip(ys, rows):
        ax.errorbar(
            rho, y, xerr=[[rho - lo], [hi - rho]], fmt="o", color=BLUE, ecolor=BLUE,
            elinewidth=1.5, capsize=3.2, capthick=1.3, ms=5.5, zorder=3,
        )
    ax.axvline(0, color="black", lw=0.8, zorder=1)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows])
    ax.set_xlim(0.0, 0.80)
    ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8])
    ax.set_ylim(-0.55, len(rows) - 0.45)
    ax.set_xlabel("Spearman \u03c1, frozen score vs donor age (held-out donors)")
    _open_spines(ax)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)

    trans = ax.get_yaxis_transform()
    ax.text(1.06, len(rows) - 0.45, "train \u2192 test n", transform=trans,
            ha="left", va="bottom", fontsize=7.5, fontweight="bold", clip_on=False)
    ax.text(1.38, len(rows) - 0.45, "\u03c1", transform=trans,
            ha="left", va="bottom", fontsize=7.5, fontweight="bold", clip_on=False)
    for y, (_label, n, rho, _lo, _hi) in zip(ys, rows):
        ax.text(1.06, y, n, transform=trans, ha="left", va="center", fontsize=7.5, clip_on=False)
        ax.text(1.38, y, _signed(rho), transform=trans, ha="left", va="center", fontsize=7.5, clip_on=False)

    axb = fig.add_axes([0.0995, 0.0994, 0.3790, 0.3012])
    g = d["gse_pts"]
    axb.scatter(g.age_years, g.age_score, s=16, c=BLUE, linewidths=0, zorder=3)
    axb.set_xlim(18, 92)
    axb.set_ylim(-14.6, 1.6)
    axb.set_xticks([20, 30, 40, 50, 60, 70, 80, 90])
    axb.set_xlabel("Donor age (years)")
    axb.set_ylabel("Frozen fibroblast score (a.u.)")
    axb.set_title("GSE226189, bulk (82 donors)", pad=4)
    axb.text(0.03, 0.96, "\u03c1 = 0.549", transform=axb.transAxes, ha="left", va="top", fontsize=8)
    _open_spines(axb)

    axc = fig.add_axes([0.5914, 0.0994, 0.3790, 0.3012])
    p = d["don_pts"]
    ph = p.true_donor.str.startswith("Reynolds_2021")
    if int(ph.sum()) != 5 or not np.allclose(p.loc[ph, "age_years"], 42.0):
        raise SystemExit("STOP: expected 5 Reynolds donors recorded at age 42")
    other = p.loc[~ph]
    axc.scatter(other.age_years, other.age_score, s=18, c=BLUE, linewidths=0, zorder=3, label="other donors (60)")
    axc.scatter(
        p.loc[ph, "age_years"], p.loc[ph, "age_score"], s=28, facecolors="none", edgecolors=ORANGE,
        linewidths=1.15, zorder=4, label="5 donors, age recorded as 42\n(author range 25\u201360)",
    )
    axc.set_xlim(14, 84)
    axc.set_ylim(-11.2, 8.4)
    axc.set_xticks([20, 30, 40, 50, 60, 70, 80])
    axc.set_xlabel("Donor age (years)")
    axc.set_title("Fibroblast atlas, 10x (65 donors)", pad=4)
    axc.text(0.03, 0.96, "\u03c1 = 0.433", transform=axc.transAxes, ha="left", va="top", fontsize=8)
    leg = axc.legend(
        frameon=False, fontsize=6.2, loc="lower left", handletextpad=0.4, borderpad=0.15,
        labelspacing=0.35,
    )
    for t in leg.get_texts():
        t.set_linespacing(1.05)
    _open_spines(axc)

    path = OUT / "figure1.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return path


def figure3(d):
    fig = plt.figure(figsize=(7.0, 3.1))
    _panel(fig, (0.012, 0.985), "A")
    _panel(fig, (0.500, 0.985), "B")

    ax = fig.add_axes([0.0995, 0.1989, 0.2210, 0.6419])
    specs = [
        (0, "Frozen\nruler", float(d["gse_ruler"].rho), float(d["gse_ruler"].rho_ci_lo),
         float(d["gse_ruler"].rho_ci_hi), BLUE, "0.55"),
        (1, "MD", float(d["gse_md"].rho), float(d["gse_md"].rho_ci_lo),
         float(d["gse_md"].rho_ci_hi), ORANGE, "0.08"),
        (2, "Age-up \u2212\nage-down", float(d["gse_s3"].rho), float(d["gse_s3"].rho_ci_lo),
         float(d["gse_s3"].rho_ci_hi), PINK, "0.04"),
    ]
    for x, _lab, rho, lo, hi, color, lab in specs:
        ax.errorbar(
            x, rho, yerr=[[rho - lo], [hi - rho]], fmt="o", color=color, ecolor=color,
            elinewidth=1.6, capsize=3.5, capthick=1.4, ms=6.5, zorder=3,
        )
        ax.text(x, hi + 0.035, lab, ha="center", va="bottom", fontsize=7.5, color=color)
    ax.axhline(0, color="black", lw=0.8, zorder=1)
    ax.set_xlim(-0.55, 2.55)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["Frozen\nruler", "MD", "Age-up \u2212\nage-down"])
    ax.set_ylim(-0.32, 0.92)
    ax.set_yticks([-0.2, 0.0, 0.2, 0.4, 0.6, 0.8])
    ax.set_ylabel("Spearman \u03c1 with donor age (95% CI)")
    ax.set_title("GSE226189, 82 donors", pad=3)
    _open_spines(ax)

    axb = fig.add_axes([0.575, 0.1989, 0.230, 0.6957])
    # Top to bottom. Grey "as run" rows are omitted.
    items = [
        ("Ruler \u2013 MD", None),
        ("atlas, 65 donors", d["atlas_md"]),
        ("GSE226189, 82 donors", d["gse_d_md"]),
        ("Ruler \u2013 age-up/down", None),
        ("atlas, 65 donors", d["atlas_s3"]),
        ("GSE226189, 82 donors", d["gse_d_s3"]),
    ]
    ys = np.arange(len(items))[::-1]
    labels = []
    for y, (lab, rec) in zip(ys, items):
        labels.append(lab if rec is None else "  " + lab)
        if rec is None:
            continue
        est, lo, hi = float(rec.delta_rho), float(rec.delta_ci_lo), float(rec.delta_ci_hi)
        axb.errorbar(
            est, y, xerr=[[est - lo], [hi - est]], fmt="o", color=BLACK, ecolor=BLACK,
            elinewidth=1.35, capsize=3.0, capthick=1.15, ms=5.2, zorder=3,
        )
        axb.text(
            1.03, y, f"{_signed(est)} [{_signed(lo)}, {_signed(hi)}]",
            transform=axb.get_yaxis_transform(), ha="left", va="center", fontsize=6.3, clip_on=False,
        )
    axb.axvline(0, color="0.35", ls=(0, (3, 2)), lw=0.8, zorder=1)
    axb.set_yticks(ys)
    axb.set_yticklabels(labels)
    for tick, (_lab, rec) in zip(axb.get_yticklabels(), items):
        if rec is None:
            tick.set_fontweight("bold")
    axb.set_xlim(-0.25, 0.95)
    axb.set_xticks([-0.2, 0.0, 0.2, 0.4, 0.6, 0.8])
    axb.set_ylim(-0.55, len(items) - 0.35)
    axb.set_xlabel("Paired \u0394\u03c1 (95% percentile CI, B = 200)")
    _open_spines(axb)
    axb.tick_params(axis="y", length=0)
    axb.spines["left"].set_visible(False)

    path = OUT / "figure3.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return path


def _check_png(path, size):
    from PIL import Image
    im = Image.open(path)
    if im.size != size:
        raise SystemExit(f"STOP: {path.name} is {im.size}, expected {size}")
    dpi = im.info.get("dpi")
    if dpi is None or abs(dpi[0] - 300) > 0.1 or abs(dpi[1] - 300) > 0.1:
        raise SystemExit(f"STOP: {path.name} dpi is {dpi}, expected 300")


def main():
    _font()
    OUT.mkdir(parents=True, exist_ok=True)
    d = _load()
    p1 = figure1(d)
    p3 = figure3(d)
    _check_png(p1, (2100, 1680))
    _check_png(p3, (2100, 930))
    print(p1.relative_to(ROOT).as_posix())
    print(p3.relative_to(ROOT).as_posix())


if __name__ == "__main__":
    main()

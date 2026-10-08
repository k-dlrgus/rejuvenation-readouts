"""Final paper figures in one colour scheme: Fig1-Fig5 and FigS1.

Only colours differ from the source figures. Every figure is re-plotted from the stored tables the
source figure was drawn from; nothing is refit or recomputed.

  Fig1  = results/newstory/fig_N1_md_own_genes        (src/newstory.py figure_n1)
  Fig2  = results/newstory/fig_N2_any_genes           (src/newstory.py figure_n2)
  Fig3  = results/paper_figs/manuscript/figure1.png   (src/paper_figs_manuscript.py figure1, 65 donors)
  Fig4  = results/paper_figs/manuscript/figure3.png   (src/paper_figs_manuscript.py figure3, 65 donors)
  Fig5  = results/newstory/fig_N5_ruler_same_mistake  (src/newstory.py figure_n5)
  FigS1 = GTEx youth-axis figure (Supplementary Figure S1 of
          paper/rejuvenation_readouts_preprint_9.docx; Figure 4 in the superseded preprint_3). No script for
          it is in the repository; the layout is rebuilt from that image, the data are
          results/toward/stats.csv, results/toward/c3_folds.csv and results/toward/summary.json.

Colours: blue = frozen ruler and all genes (23,485); vermillion = MD; reddish purple = age-up/age-down
lists and age genes; black and greys = everything else. With both donors in a panel, 96-year-old
donor = black and 22-year-old donor = mid grey.

Preprint v1.2 layout: every figure is drawn at the text width of
paper/rejuvenation_readouts_preprint_v1.2_draft.docx (7.5 in: US Letter, 0.5 in margins) so it is
inserted at 100%, at most 3.0 in tall for one row of panels and 4.5 in for two. One style throughout:
Arial, every label 7 pt, panel letters 9 pt bold, one axis/reference line width, one data line width,
one marker size. Only size, layout and styling changed from the v12 figures in results/figures_final/
(src/figures_final.py at 77b74e0); src/paper_figs_v1_2_check.py shows the plotted data are identical.

Writes only results/paper_figs_v1_2/ (PNG 300 dpi, PDF and plotting-data CSV per figure).

Usage: python src/figures_final.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import RESULTS, ROOT  # noqa: E402

OUT = RESULTS / "paper_figs_v1_2"
NEWSTORY = RESULTS / "newstory"
TOWARD = RESULTS / "toward"
FIBRO3_EXTRAP = RESULTS / "fibro3" / "t2_extrap.csv"
FIBRO3_CORR = RESULTS / "fibro3" / "t2_correlation.csv"

AGED_LINE, YOUNG_LINE = "GM00731", "GM23815"
LINES = (AGED_LINE, YOUNG_LINE)

BLUE = "#0072B2"
VERMILLION = "#D55E00"
PURPLE = "#CC79A7"
BLACK = "#000000"
GREY = "#999999"
LIGHTGREY = "#C8C8C8"
DARKGREY = "#555555"
DONOR_COLOR = {AGED_LINE: BLACK, YOUNG_LINE: GREY}

SPACE_COLOR = dict(MD=VERMILLION, AGE=PURPLE, FULL=BLUE)
SPACE_WORDS = dict(MD="MD genes (205)", AGE="Age genes (3,089)", FULL="All genes (23,485)")
RAND_SPACES = ("MD", "AGE")
F_VALUES = (0.0, 0.10, 0.25, 0.50)
STATE_WORDS = dict(Fibroblast_d0="Day 0", PartialReprog="Partial\nreprogramming",
                   EarlyPluripotency="Early\npluripotency", Pluripotency="Pluripotency",
                   NonReprog="Non-\nreprogrammed")

# Fig3 / Fig4: stored true-donor ruler Spearman (results/paper_figs/cxg_true_donor_results.csv).
ATLAS_RHO = 0.43338041757870716

TEXT_W = 7.5            # in; figure width = docx text width
MAX_H = {1: 3.0, 2: 4.5}  # in; height limit by number of panel rows
FS = 7                  # every label, tick label, title, legend entry and annotation
FS_LETTER = 9           # panel letters, bold
LW_AXIS = 0.6           # spines, ticks, and reference lines (zero lines, dashed guides, circle)
LW_DATA = 1.0           # data lines and error bars
MS = 4.0                # every data marker (scatter: s = MS ** 2)
MS_BG = 2.0             # Fig2 background cloud of 400 random gene sets
MS_STAR = 6.0           # target star in the plane panels; a star at MS reads smaller than a circle
MEW = 0.8               # edge width of open markers
CAP = 2.0               # error-bar cap size

STYLE_RC = {
    "font.family": "Arial", "font.size": FS, "axes.labelsize": FS, "axes.titlesize": FS,
    "xtick.labelsize": FS, "ytick.labelsize": FS, "legend.fontsize": FS,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": LW_AXIS, "xtick.major.width": LW_AXIS, "ytick.major.width": LW_AXIS,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "lines.linewidth": LW_DATA,
    "lines.markersize": MS, "lines.markeredgewidth": MEW, "errorbar.capsize": CAP,
    "axes.titlepad": 4, "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 300,
    "savefig.facecolor": "white", "figure.facecolor": "white",
}


def use_style(rc=STYLE_RC):
    matplotlib.rc_file_defaults()
    plt.rcParams.update(rc)


def new_figure(height, rows):
    if height > MAX_H[rows]:
        raise SystemExit(f"STOP: {height} in is over the {MAX_H[rows]} in limit for {rows} row(s)")
    return plt.figure(figsize=(TEXT_W, height))


def axes_in(fig, left, bottom, width, height):
    """Axes placed in inches from the lower-left corner of the figure."""
    w, h = fig.get_size_inches()
    return fig.add_axes([left / w, bottom / h, width / w, height / h])


def panel_letters(fig, rows, pad=0.06):
    """Letters at the left edge of each panel's labelled area, one baseline per row of panels."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    _, h = fig.get_size_inches()
    for row in rows:
        boxes = [ax.get_tightbbox(r) for ax, _ in row]
        top = max(inv.transform((0, b.y1))[1] for b in boxes)
        for (ax, s), b in zip(row, boxes):
            x0 = inv.transform((b.x0, 0))[0]
            fig.text(x0, top + pad / h, s, fontsize=FS_LETTER, fontweight="bold", ha="left", va="bottom")


def save(fig, stem, data):
    w, h = fig.get_size_inches()
    if abs(w - TEXT_W) > 1e-9:
        raise SystemExit(f"STOP: {stem} is {w} in wide, not the text width {TEXT_W} in")
    bb = fig.get_tightbbox(fig.canvas.get_renderer())
    if bb.x0 < 0 or bb.y0 < 0 or bb.x1 > w or bb.y1 > h:
        raise SystemExit(f"STOP: {stem} draws outside the figure: {bb.bounds} vs {w} x {h} in")
    fig.savefig(OUT / f"{stem}.png", dpi=300)
    fig.savefig(OUT / f"{stem}.pdf")
    plt.close(fig)
    pd.DataFrame(data).to_csv(OUT / f"{stem}_data.csv", index=False)
    print(f"wrote {(OUT / stem).relative_to(ROOT).as_posix()}.{{png,pdf}} and _data.csv")


# --------------------------------------------------------------------------- newstory helpers
def plane_width(height, xlim, ylim):
    """Width in inches of an equal-aspect plane panel of the given height."""
    return height * (xlim[1] - xlim[0]) / (ylim[1] - ylim[0])


def draw_plane_frame(ax, xlim, ylim, circle_label_xy=(1.3, 0.5)):
    from matplotlib.patches import Circle
    ax.add_patch(Circle((1.0, 0.0), 1.0, fill=False, ls=(0, (4, 2.5)), lw=LW_AXIS, ec=GREY, zorder=1))
    if circle_label_xy:
        ax.text(*circle_label_xy, "closer than\nat start", ha="center", va="center", fontsize=FS, color=DARKGREY)
    ax.axhline(0, color=BLACK, lw=LW_AXIS, zorder=0)
    ax.plot([0], [0], marker="s", ms=MS, mfc="white", mec=BLACK, mew=MEW, ls="none", zorder=5, clip_on=False)
    ax.plot([1], [0], marker="*", ms=MS_STAR, mfc=BLACK, mec=BLACK, ls="none", zorder=5, clip_on=False)
    ax.annotate("start", (0, 0), xytext=(-4, 3), textcoords="offset points", fontsize=FS, ha="right", va="bottom")
    ax.annotate("target", (1, 0), xytext=(4.5, 3), textcoords="offset points", fontsize=FS, ha="left", va="bottom")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Progress along start-to-target axis\n(axis lengths)")
    ax.set_ylabel("Sideways distance (axis lengths)")


# --------------------------------------------------------------------------- Fig1
def fig1():
    use_style()
    st = pd.read_csv(NEWSTORY / "partA_md_by_state.csv")
    same = pd.read_csv(NEWSTORY / "partB_same_distances.csv")
    mix = pd.read_csv(NEWSTORY / "partB_mixture.csv")
    data = []

    fig = new_figure(2.8, rows=1)
    plane_xlim, plane_ylim = (-0.75, 2.05), (0, 4.25)
    # A
    ax = axes_in(fig, 0.55, 0.95, 3.1, 1.55)
    order = ("Fibroblast_d0", "PartialReprog", "EarlyPluripotency", "Pluripotency")
    aged = st[st.cell_line == AGED_LINE].set_index("state").loc[list(order)]
    young0 = st[(st.cell_line == YOUNG_LINE) & (st.state == "Fibroblast_d0")].iloc[0]
    x = np.arange(len(order))
    ax.axhline(young0.joint_mean, color=DARKGREY, ls=(0, (4, 2.5)), lw=LW_AXIS, zorder=1)
    ax.text(len(order) - 0.55, young0.joint_mean + 0.012, "22-year-old donor, day 0", ha="right", va="bottom",
            fontsize=FS, color=DARKGREY)
    ax.errorbar(x, aged.joint_mean, yerr=[aged.joint_mean - aged.joint_ci_lo, aged.joint_ci_hi - aged.joint_mean],
                fmt="o", ms=MS, color=VERMILLION, ecolor=VERMILLION, elinewidth=LW_DATA, capsize=CAP, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels([STATE_WORDS[s].replace("\n", " ") for s in order], rotation=30, ha="right",
                       rotation_mode="anchor")
    ax.set_xlim(-0.5, len(order) - 0.5)
    ax.set_ylim(-0.05, 0.6)
    ax.set_ylabel("MD score (mean, 95% CI)")
    ax.set_xlabel("Aged donor's cells")
    ax_a = ax
    for s, r in aged.iterrows():
        data.append(dict(panel="A", series="aged donor", label=s, x=float(order.index(s)), y=r.joint_mean,
                         y_lo=r.joint_ci_lo, y_hi=r.joint_ci_hi, n_cells=int(r.n_cells)))
    data.append(dict(panel="A", series="young donor day 0 (dashed line)", label="Fibroblast_d0", x=np.nan,
                     y=young0.joint_mean, y_lo=young0.joint_ci_lo, y_hi=young0.joint_ci_hi,
                     n_cells=int(young0.n_cells)))

    # B
    plane_h = 2.0
    ax = axes_in(fig, 4.9, 0.5, plane_width(plane_h, plane_xlim, plane_ylim), plane_h)
    draw_plane_frame(ax, plane_xlim, plane_ylim, circle_label_xy=(1.0, 0.57))
    md = same[same.space == "MD"].set_index("test")
    mm = mix[mix.space == "MD"].sort_values("f")
    for _, r in mm.iterrows():
        ax.plot(r.progress, r.perp, marker="o", ms=MS, ls="none", mec=LIGHTGREY, mew=MEW,
                mfc="white" if r.f == 0 else LIGHTGREY, zorder=3)
        data.append(dict(panel="B", series="mixture", label=f"f={r.f:.2f}", x=r.progress, y=r.perp))
    ax.annotate("young-cell\nmixtures\n(0–50%)", xy=(mm.progress.iloc[1], mm.perp.iloc[1]), xytext=(-0.22, 0.78),
                fontsize=FS, color=DARKGREY, ha="center", va="bottom",
                arrowprops=dict(arrowstyle="-", color=GREY, lw=LW_AXIS, shrinkA=1, shrinkB=3))
    f = md.loc["forward"]
    ax.plot(f.progress, f.perp, marker="o", ms=MS, color=VERMILLION, ls="none", zorder=4)
    ax.text(f.progress + 0.12, f.perp, "Partial\nreprogramming", fontsize=FS, color=BLACK, va="center")
    for test, word in (("Pluri", "Pluripotency"), ("NonReprog", "Non-reprogrammed")):
        r = md.loc[test]
        ax.plot(r.progress, r.perp, marker="o", ms=MS, mfc="white", mec=BLACK, mew=MEW, ls="none", zorder=4)
        ax.text(r.progress + (0.1 if test == "NonReprog" else -0.1), r.perp + (-0.02 if test == "NonReprog" else 0.0),
                word, fontsize=FS, ha="left" if test == "NonReprog" else "right", va="center")
    r = md.loc["reverse"]
    ax.plot(r.progress, r.perp, marker="o", ms=MS, color=VERMILLION, ls="none", zorder=4)
    ax.text(r.progress, r.perp + 0.2, "reverse", fontsize=FS, color=BLACK, ha="center", va="bottom")
    ax.set_xlabel("Progress toward young donor\n(axis lengths)")
    panel_letters(fig, [[(ax_a, "A"), (ax, "B")]])
    for test in ("forward", "reverse", "Pluri", "NonReprog"):
        r = md.loc[test]
        data.append(dict(panel="B", series=test, label=test, x=r.progress, y=r.perp, n_cells=int(r.n_cells),
                         note=f"final/start={r.final_over_start:.4f}"))
    save(fig, "Fig1", data)


# --------------------------------------------------------------------------- Fig2
def fig2():
    use_style()
    from matplotlib.lines import Line2D
    same = pd.read_csv(NEWSTORY / "partB_same_distances.csv")
    mix = pd.read_csv(NEWSTORY / "partB_mixture.csv")
    rnd = pd.read_csv(NEWSTORY / "partB_random_sets_per_set.csv")
    data = []
    rand_style = dict(MD=dict(marker="o", color=LIGHTGREY, ms=MS_BG, word="random 205-gene sets (200)"),
                      AGE=dict(marker="^", color=LIGHTGREY, ms=MS_BG, word="random 3,089-gene sets (200)"))

    fig = new_figure(3.0, rows=1)
    xlim, ylim = (-1.15, 2.05), (0, 3.05)
    plane_h, bottom = 1.8, 0.95
    plane_w = plane_width(plane_h, xlim, ylim)
    lefts = (0.5, 0.5 + plane_w + 0.5)
    letters = []
    for k, (test, pre, letter) in enumerate((("forward", "fwd", "A"), ("reverse", "rev", "B"))):
        ax = axes_in(fig, lefts[k], bottom, plane_w, plane_h)
        draw_plane_frame(ax, xlim, ylim)
        for sp in RAND_SPACES:
            rs = rnd[rnd.space_size == sp]
            sty = rand_style[sp]
            ax.plot(rs[f"{pre}_progress"], rs[f"{pre}_perp"], ls="none", marker=sty["marker"], ms=sty["ms"],
                    mfc=sty["color"], mec="none", alpha=0.9, zorder=2)
            for _, r in rs.iterrows():
                data.append(dict(panel=letter, series=f"random {sp}-size", label=f"set {int(r.set_index)}",
                                 x=r[f"{pre}_progress"], y=r[f"{pre}_perp"], note=f"delta={r[f'{pre}_delta']:.4f}"))
        for sp in ("FULL", "AGE", "MD"):
            r = same[(same.space == sp) & (same.test == test)].iloc[0]
            ax.plot(r.progress, r.perp, marker="o", ms=MS, color=SPACE_COLOR[sp], mec="white", mew=MEW,
                    ls="none", zorder=4)
            data.append(dict(panel=letter, series=f"real {sp}", label=SPACE_WORDS[sp], x=r.progress, y=r.perp,
                             n_cells=int(r.n_cells), note=f"final/start={r.final_over_start:.4f}"))
        if test == "forward":
            ax.set_xlabel("Progress toward young donor\n(axis lengths)")
        else:
            ax.set_xlabel("Progress toward aged donor\n(axis lengths)")
            ax.text(0.97, 0.97, "reverse", transform=ax.transAxes, fontsize=FS, color=BLACK,
                    ha="right", va="top")
        letters.append((ax, letter))

    c_left = lefts[1] + plane_w + 0.75
    ax = axes_in(fig, c_left, bottom, TEXT_W - 0.1 - c_left, plane_h)
    for sp in ("MD", "AGE", "FULL"):
        mm = mix[mix.space == sp].sort_values("f")
        ax.plot(mm.f, mm.rel_delta, color=SPACE_COLOR[sp], lw=LW_DATA, zorder=2)
        for _, r in mm.iterrows():
            ax.plot(r.f, r.rel_delta, marker="o", ms=MS, ls="none", mec=SPACE_COLOR[sp], mew=MEW,
                    mfc="white" if r.f == 0 else SPACE_COLOR[sp], zorder=3)
            data.append(dict(panel="C", series=f"mixture {sp}", label=f"f={r.f:.2f}", x=r.f, y=r.rel_delta,
                             note=f"delta={r.delta:.4f}"))
    ax.axhline(0, color=BLACK, lw=LW_AXIS)
    ax.set_xticks(F_VALUES)
    ax.set_xticklabels(["0", "0.10", "0.25", "0.50"])
    ax.set_xlim(-0.03, 0.53)
    ax.set_xlabel("Fraction of young cells in mixture")
    ax.set_ylabel("Change in distance ÷ starting distance")
    letters.append((ax, "C"))
    panel_letters(fig, [letters])

    handles = [Line2D([], [], marker="o", ls="none", ms=MS, color=SPACE_COLOR[sp], label=SPACE_WORDS[sp])
               for sp in ("MD", "AGE", "FULL")]
    # Legend symbols for the background cloud at MS: at MS_BG they are too small to read (v12 also enlarged them).
    handles += [Line2D([], [], marker=rand_style[sp]["marker"], ls="none", ms=MS, mfc=rand_style[sp]["color"],
                       mec="none", label=rand_style[sp]["word"]) for sp in RAND_SPACES]
    handles += [Line2D([], [], marker="s", ls="none", ms=MS, mfc="white", mec=BLACK, mew=MEW, label="start"),
                Line2D([], [], marker="*", ls="none", ms=MS_STAR, color=BLACK, label="target (young donor, day 0; "
                       "reverse: aged donor)"),
                Line2D([], [], ls=(0, (4, 2.5)), lw=LW_AXIS, color=GREY, label="closer than at start")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.0),
               columnspacing=1.2, handletextpad=0.4)
    save(fig, "Fig2", data)


# --------------------------------------------------------------------------- Fig5
def fig5():
    use_style()
    from matplotlib.lines import Line2D
    from scipy.stats import spearmanr
    tw = pd.read_csv(TOWARD / "stats.csv", float_precision="round_trip")
    ex = pd.read_csv(FIBRO3_EXTRAP, float_precision="round_trip")
    corr = pd.read_csv(FIBRO3_CORR, float_precision="round_trip")
    data = []

    fig = new_figure(3.0, rows=1)
    # A
    ax = axes_in(fig, 0.55, 0.85, 2.3, 1.82)
    order = ("Fibroblast", "PartialReprog", "EarlyPluripotency", "Pluripotency")
    aged = tw[tw.cell_line == AGED_LINE].set_index("state").loc[list(order)]
    young0 = tw[(tw.cell_line == YOUNG_LINE) & (tw.state == "Fibroblast")].iloc[0]
    x = np.arange(len(order))
    ax.axhline(young0.ruler_score_context, color=DARKGREY, ls=(0, (4, 2.5)), lw=LW_AXIS, zorder=1)
    ax.text(len(order) - 0.55, young0.ruler_score_context + 0.35, "22-year-old donor, day 0", ha="right",
            va="bottom", fontsize=FS, color=DARKGREY)
    ax.plot(x, aged.ruler_score_context, marker="o", ms=MS, ls="none", color=BLUE, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels([STATE_WORDS["Fibroblast_d0" if s == "Fibroblast" else s].replace("\n", " ") for s in order],
                       rotation=30, ha="right", rotation_mode="anchor")
    ax.set_xlim(-0.5, len(order) - 0.5)
    ax.set_ylim(-12, 5)
    ax.set_ylabel("Frozen ruler score (lower = scored younger)")
    ax.set_xlabel("Aged donor's cells")
    ax_a = ax
    for s, r in aged.iterrows():
        data.append(dict(panel="A", series="aged donor", label=s, x=float(order.index(s)), y=r.ruler_score_context,
                         n_cells=int(r.n_cells), note="Fibroblast = day-0 cells only" if s == "Fibroblast" else ""))
    data.append(dict(panel="A", series="young donor day 0 (dashed line)", label="Fibroblast", x=np.nan,
                     y=young0.ruler_score_context, n_cells=int(young0.n_cells), note="Fibroblast = day-0 cells only"))

    # B
    ax = axes_in(fig, 3.6, 0.6, TEXT_W - 0.1 - 3.6, 2.07)
    sc = ex[~ex.below_min_cells.astype(bool) & np.isfinite(ex.age_score.to_numpy(float))]
    ref = corr[(corr.distance == "nn_euclidean") & (corr.scope == "all_donor_cluster_timepoint_n_ge_20")].iloc[0]
    rho_chk = float(spearmanr(sc.nn_euclidean, sc.age_score).statistic)
    if len(sc) != int(ref.n) or abs(rho_chk - float(ref.rho)) > 1e-12:
        raise RuntimeError(f"t2_extrap rows (n={len(sc)}, rho={rho_chk}) do not match t2_correlation.csv")
    donor_word = {AGED_LINE: "96-y donor", YOUNG_LINE: "22-y donor"}
    day_marker = {0: "o", 3: "s", 7: "^", 10: "D"}
    day_ms = {d: MS for d in day_marker}
    for _, r in sc.iterrows():
        d = int(r.day)
        ax.plot(r.nn_euclidean, r.age_score, marker=day_marker[d], ms=day_ms[d], ls="none", mec="none",
                color=DONOR_COLOR[r.cell_line], zorder=3)
        data.append(dict(panel="B", series=donor_word[r.cell_line], label=f"{r.cell_line} day {d} cluster {int(r.cluster)}",
                         x=r.nn_euclidean, y=r.age_score, n_cells=int(r.n_cells), note=f"day={d}"))
    ax.text(0.98, 0.98, f"ρ = {ref.rho:.3f}\n95% CI [{ref.rho_ci_lo:.3f}, {ref.rho_ci_hi:.3f}]\n"
            f"n = {int(ref.n)} pseudobulks".replace("-", "−"), transform=ax.transAxes, ha="right", va="top",
            fontsize=FS)
    ax.set_ylim(-20, 10)
    ax.set_xlabel("Nearest-neighbor distance to GTEx training donors\n(PCA, k = 50)")
    ax.set_ylabel("Frozen-ruler score (a.u.)")
    donor_h = [Line2D([], [], marker="o", ls="none", ms=MS, mec="none", color=DONOR_COLOR[l], label=donor_word[l])
               for l in LINES]
    day_h = [Line2D([], [], marker=day_marker[d], ls="none", ms=day_ms[d], mec="none", color=LIGHTGREY,
                    label=f"day {d}") for d in day_marker]
    leg = ax.legend(handles=donor_h, loc="lower left", frameon=False, handletextpad=0.3, borderaxespad=0.3)
    ax.add_artist(leg)
    ax.legend(handles=day_h, loc="lower left", bbox_to_anchor=(0.36, 0.0), ncol=2, frameon=False,
              handletextpad=0.3, columnspacing=1.0, borderaxespad=0.3)
    data.append(dict(panel="B", series="annotation", label="spearman nn_euclidean", x=np.nan, y=float(ref.rho),
                     y_lo=float(ref.rho_ci_lo), y_hi=float(ref.rho_ci_hi), n_cells=np.nan,
                     note=f"n={int(ref.n)}; from {FIBRO3_CORR.relative_to(ROOT).as_posix()}"))
    panel_letters(fig, [[(ax_a, "A"), (ax, "B")]])
    save(fig, "Fig5", data)


# --------------------------------------------------------------------------- manuscript helpers (Fig3, Fig4)
def _open_spines(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=2.5, width=LW_AXIS)


def _signed(x, digits=3):
    """Fixed-point with an explicit sign and a unicode minus."""
    s = f"{x:+.{digits}f}"
    return s.replace("-", "\u2212")


def _load_manuscript():
    transfer = pd.read_csv(RESULTS / "fibro" / "stage1_transfer.csv")
    raw = transfer[(transfer.regime == "raw") & (transfer.method == "ridge")].copy()
    b1 = raw[raw.direction == "B1→C1"].iloc[0]
    c1 = raw[raw.direction == "C1→B1"].iloc[0]

    inst = pd.read_csv(RESULTS / "md3" / "t2_instruments.csv")
    gse_ruler = inst[inst.tag == "gse226189_ruler_check"].iloc[0]
    gse_md = inst[inst.tag == "gse226189_md_check"].iloc[0]
    gse_s3 = inst[inst.tag == "GSE226189_S3"].iloc[0]

    donor_res = pd.read_csv(RESULTS / "paper_figs" / "cxg_true_donor_results.csv")
    atlas = donor_res[(donor_res.analysis == "true_donor")].copy()
    atlas_rho = atlas[atlas.contrast == "ruler_spearman"].iloc[0]
    atlas_md = atlas[atlas.contrast == "ruler_minus_MD_AddModuleScore"].iloc[0]
    atlas_s3 = atlas[atlas.contrast == "ruler_minus_age_up_minus_age_down"].iloc[0]
    if abs(float(atlas_rho.rho) - ATLAS_RHO) > 1e-12:
        raise SystemExit("STOP: stored 65-donor rho is not 0.43338041757870716")
    if int(atlas_rho.n_units) != 65:
        raise SystemExit("STOP: atlas ruler row is not 65 donors")

    delta = pd.read_csv(RESULTS / "md3" / "t2_delta_rho.csv")
    gse_d_md = delta[delta.tag == "gse226189_ruler_minus_MD"].iloc[0]
    gse_d_s3 = delta[delta.tag == "gse226189_ruler_minus_S3"].iloc[0]

    gse_pts = pd.read_csv(RESULTS / "fibro2" / "ta_GSE226189_near_miss_bulk_scores.csv")
    if len(gse_pts) != 82:
        raise SystemExit(f"STOP: GSE226189 score table has {len(gse_pts)} rows, expected 82")

    don_pts = pd.read_csv(RESULTS / "paper_figs" / "cxg_true_donor_scores.csv")
    if len(don_pts) != 65:
        raise SystemExit(f"STOP: true-donor score table has {len(don_pts)} rows, expected 65")
    rho_pts = pd.Series(don_pts.age_score).corr(don_pts.age_years, method="spearman")
    if abs(float(rho_pts) - float(atlas_rho.rho)) > 1e-12:
        raise SystemExit("STOP: 65-donor points do not reproduce the stored rho")

    n_train = int(json.loads((RESULTS / "fibro" / "freeze.json").read_text(encoding="utf-8"))["n_donors"])
    if n_train != 652:
        raise SystemExit(f"STOP: frozen ruler n_donors is {n_train}, expected 652")

    return dict(
        b1=b1, c1=c1, gse_ruler=gse_ruler, gse_md=gse_md, gse_s3=gse_s3,
        atlas_rho=atlas_rho, atlas_md=atlas_md, atlas_s3=atlas_s3,
        gse_d_md=gse_d_md, gse_d_s3=gse_d_s3,
        gse_pts=gse_pts, don_pts=don_pts, n_train=n_train,
    )


# --------------------------------------------------------------------------- Fig3
def fig3(d):
    use_style()
    data = []
    # Two rows: A (transfer forest with its n and rho columns) on top, B and C side by side below.
    fig = new_figure(4.5, rows=2)

    a_left, a_width = 2.0, 3.2
    ax = axes_in(fig, a_left, 2.85, a_width, 1.3)
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
    for y, (label, n, rho, lo, hi) in zip(ys, rows):
        ax.errorbar(
            rho, y, xerr=[[rho - lo], [hi - rho]], fmt="o", color=BLUE, ecolor=BLUE,
            elinewidth=LW_DATA, capsize=CAP, capthick=LW_DATA, ms=MS, zorder=3,
        )
        data.append(dict(panel="A", series="ruler transfer", label=label, x=rho, x_lo=lo, x_hi=hi, y=float(y),
                         note=f"train \u2192 test n: {n}"))
    ax.axvline(0, color="black", lw=LW_AXIS, zorder=1)
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
    col_n, col_rho = 1 + 0.2 / a_width, 1 + 1.1 / a_width  # table columns 0.2 in and 1.1 in right of the axes
    ax.text(col_n, len(rows) - 0.45, "train \u2192 test n", transform=trans,
            ha="left", va="bottom", fontsize=FS, fontweight="bold", clip_on=False)
    ax.text(col_rho, len(rows) - 0.45, "\u03c1", transform=trans,
            ha="left", va="bottom", fontsize=FS, fontweight="bold", clip_on=False)
    for y, (_label, n, rho, _lo, _hi) in zip(ys, rows):
        ax.text(col_n, y, n, transform=trans, ha="left", va="center", fontsize=FS, clip_on=False)
        ax.text(col_rho, y, _signed(rho), transform=trans, ha="left", va="center", fontsize=FS, clip_on=False)

    axb = axes_in(fig, 0.55, 0.42, 2.95, 1.5)
    g = d["gse_pts"]
    axb.scatter(g.age_years, g.age_score, s=MS ** 2, c=BLUE, linewidths=0, zorder=3)
    for _, r in g.iterrows():
        data.append(dict(panel="B", series="GSE226189 donor", label=r["sample"], x=float(r.age_years),
                         y=float(r.age_score)))
    axb.set_xlim(18, 92)
    axb.set_ylim(-14.6, 1.6)
    axb.set_xticks([20, 30, 40, 50, 60, 70, 80, 90])
    axb.set_xlabel("Donor age (years)")
    axb.set_ylabel("Frozen fibroblast score (a.u.)")
    axb.set_title("GSE226189, bulk (82 donors)")
    axb.text(0.03, 0.96, "\u03c1 = 0.549", transform=axb.transAxes, ha="left", va="top", fontsize=FS)
    data.append(dict(panel="B", series="annotation", label="ruler spearman", y=float(d["gse_ruler"].rho),
                     y_lo=float(d["gse_ruler"].rho_ci_lo), y_hi=float(d["gse_ruler"].rho_ci_hi),
                     note="results/md3/t2_instruments.csv gse226189_ruler_check"))
    _open_spines(axb)

    axc = axes_in(fig, 4.3, 0.42, 2.95, 1.5)
    p = d["don_pts"]
    ph = p.true_donor.str.startswith("Reynolds_2021")
    if int(ph.sum()) != 5 or not np.allclose(p.loc[ph, "age_years"], 42.0):
        raise SystemExit("STOP: expected 5 Reynolds donors recorded at age 42")
    other = p.loc[~ph]
    axc.scatter(other.age_years, other.age_score, s=MS ** 2, c=BLUE, linewidths=0, zorder=3, label="other donors (60)")
    axc.scatter(
        p.loc[ph, "age_years"], p.loc[ph, "age_score"], s=MS ** 2, facecolors="none", edgecolors=BLACK,
        linewidths=MEW, zorder=4, label="5 donors, age recorded as 42\n(author range 25\u201360)",
    )
    for _, r in p.iterrows():
        placeholder = r.true_donor.startswith("Reynolds_2021")
        data.append(dict(panel="C", series="donor, age recorded as 42" if placeholder else "other donor",
                         label=r.true_donor, x=float(r.age_years), y=float(r.age_score),
                         n=int(r.n_pseudobulks), note=f"n_cells={int(r.n_cells)}"))
    axc.set_xlim(14, 84)
    axc.set_ylim(-11.2, 8.4)
    axc.set_xticks([20, 30, 40, 50, 60, 70, 80])
    # v12 tick step; at this panel height the automatic locator would switch to steps of 5.
    axc.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(2.5))
    axc.set_xlabel("Donor age (years)")
    axc.set_title("Fibroblast atlas, 10x (65 donors)")
    axc.text(0.03, 0.96, "\u03c1 = 0.433", transform=axc.transAxes, ha="left", va="top", fontsize=FS)
    data.append(dict(panel="C", series="annotation", label="ruler spearman", y=float(d["atlas_rho"].rho),
                     y_lo=float(d["atlas_rho"].rho_ci_lo), y_hi=float(d["atlas_rho"].rho_ci_hi),
                     n=int(d["atlas_rho"].n_units), note="results/paper_figs/cxg_true_donor_results.csv true_donor"))
    leg = axc.legend(
        frameon=False, fontsize=FS, loc="lower left", handletextpad=0.4, borderpad=0.15,
        labelspacing=0.35,
    )
    for t in leg.get_texts():
        t.set_linespacing(1.05)
    _open_spines(axc)
    panel_letters(fig, [[(ax, "A")], [(axb, "B"), (axc, "C")]])
    save(fig, "Fig3", data)


# --------------------------------------------------------------------------- Fig4
def fig4(d):
    use_style()
    data = []
    fig = new_figure(2.8, rows=1)

    ax = axes_in(fig, 0.6, 0.5, 2.0, 1.95)
    specs = [
        (0, "Frozen\nruler", float(d["gse_ruler"].rho), float(d["gse_ruler"].rho_ci_lo),
         float(d["gse_ruler"].rho_ci_hi), BLUE, "0.55"),
        (1, "MD", float(d["gse_md"].rho), float(d["gse_md"].rho_ci_lo),
         float(d["gse_md"].rho_ci_hi), VERMILLION, "0.08"),
        (2, "Age-up \u2212\nage-down", float(d["gse_s3"].rho), float(d["gse_s3"].rho_ci_lo),
         float(d["gse_s3"].rho_ci_hi), PURPLE, "0.04"),
    ]
    for x, name, rho, lo, hi, color, lab in specs:
        ax.errorbar(
            x, rho, yerr=[[rho - lo], [hi - rho]], fmt="o", color=color, ecolor=color,
            elinewidth=LW_DATA, capsize=CAP, capthick=LW_DATA, ms=MS, zorder=3,
        )
        ax.text(x, hi + 0.035, lab, ha="center", va="bottom", fontsize=FS, color=BLACK)
        data.append(dict(panel="A", series=name.replace("\n", " "), label=lab, x=float(x), y=rho, y_lo=lo, y_hi=hi,
                         n=82))
    ax.axhline(0, color="black", lw=LW_AXIS, zorder=1)
    ax.set_xlim(-0.55, 2.55)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["Frozen\nruler", "MD", "Age-up \u2212\nage-down"])
    ax.set_ylim(-0.32, 0.92)
    ax.set_yticks([-0.2, 0.0, 0.2, 0.4, 0.6, 0.8])
    ax.set_ylabel("Spearman \u03c1 with donor age (95% CI)")
    ax.set_title("GSE226189, 82 donors")
    _open_spines(ax)

    axb = axes_in(fig, 4.35, 0.5, 1.6, 2.1)
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
    group = None
    for y, (lab, rec) in zip(ys, items):
        labels.append(lab if rec is None else "  " + lab)
        if rec is None:
            group = lab
            continue
        est, lo, hi = float(rec.delta_rho), float(rec.delta_ci_lo), float(rec.delta_ci_hi)
        axb.errorbar(
            est, y, xerr=[[est - lo], [hi - est]], fmt="o", color=BLACK, ecolor=BLACK,
            elinewidth=LW_DATA, capsize=CAP, capthick=LW_DATA, ms=MS, zorder=3,
        )
        axb.text(
            1.03, y, f"{_signed(est)} [{_signed(lo)}, {_signed(hi)}]",
            transform=axb.get_yaxis_transform(), ha="left", va="center", fontsize=FS, clip_on=False,
        )
        data.append(dict(panel="B", series=group, label=lab, x=est, x_lo=lo, x_hi=hi, y=float(y)))
    axb.axvline(0, color="0.35", ls=(0, (3, 2)), lw=LW_AXIS, zorder=1)
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
    panel_letters(fig, [[(ax, "A"), (axb, "B")]])
    save(fig, "Fig4", data)


# --------------------------------------------------------------------------- FigS1
def figS1():
    use_style()
    tw = pd.read_csv(TOWARD / "stats.csv", float_precision="round_trip")
    folds = pd.read_csv(TOWARD / "c3_folds.csv", float_precision="round_trip").sort_values("fold")
    c3 = json.loads((TOWARD / "summary.json").read_text(encoding="utf-8"))["c3"]
    if len(folds) != c3["n_folds"] or abs(folds.cos.min() - c3["min_cos"]) > 1e-12:
        raise SystemExit("STOP: c3_folds.csv does not match results/toward/summary.json")
    bar = float(c3["bar"])
    data = []

    order = ["Fibroblast", "PartialReprog", "EarlyPluripotency", "Pluripotency"]
    xs = [0, 1, 2, 3]
    x_non = 4.3
    ticks = ["Fibro.\nday 0", "Partial\nreprog.", "Early\npluri.", "Pluri-\npotent", "Non-\nreprog."]
    donor_label = {AGED_LINE: "96-y donor (GM00731)", YOUNG_LINE: "22-y donor (GM23815)"}

    fig = new_figure(4.5, rows=2)
    pw, ph = 2.95, 1.4
    axs = np.array([[axes_in(fig, 0.6, 2.6, pw, ph), axes_in(fig, 4.4, 2.6, pw, ph)],
                    [axes_in(fig, 0.6, 0.45, pw, ph), axes_in(fig, 4.4, 0.45, pw, ph)]])
    specs = [
        (axs[0, 0], "A", "cos_S", "Alignment with youth axis", "cos(displacement, GTEx young \u2212 old)", 1),
        (axs[0, 1], "B", "delta_young", "Change in distance to young", "\u0394 distance to GTEx young centroid", 0),
        (axs[1, 0], "C", "ruler_score_context", "One-dimensional age score", "Frozen-ruler score (a.u.; context)", 0),
    ]
    for ax, letter, col, title, ylabel, first in specs:
        for line in LINES:
            dl = tw[tw.cell_line == line].set_index("state")
            y = dl.loc[order[first:], col].to_numpy(float)
            ax.plot(xs[first:], y, "-o", color=DONOR_COLOR[line], lw=LW_DATA, ms=MS, mew=MEW, label=donor_label[line])
            y_non = float(dl.loc["NonReprog", col])
            ax.plot([x_non], [y_non], "o", mfc="white", mec=DONOR_COLOR[line], ms=MS, mew=MEW)
            for x, st, v in zip(xs[first:] + [x_non], order[first:] + ["NonReprog"], list(y) + [y_non]):
                data.append(dict(panel=letter, series=donor_label[line], label=st, x=float(x), y=float(v),
                                 n_cells=int(dl.loc[st, "n_cells"]),
                                 note=f"results/toward/stats.csv {col}" + ("; open marker" if st == "NonReprog" else "")))
        if col != "ruler_score_context":
            ax.axhline(0, color=BLACK, lw=LW_AXIS)
        ax.set_xticks(xs + [x_non])
        ax.set_xticklabels(ticks)
        ax.set_xlim(-0.4, 4.7)
        ax.set_title(title)
        ax.set_ylabel(ylabel)
    axs[0, 0].set_ylim(-0.02, 0.16)
    axs[0, 0].set_yticks(np.arange(0, 0.1501, 0.025))
    axs[0, 0].legend(loc="upper left", frameon=False, fontsize=FS)

    ax = axs[1, 1]
    fx = 1 + (folds.fold.to_numpy(float) - 2) * 0.06
    ax.plot(fx, folds.cos, "o", color="#333333", ms=MS, mew=MEW, ls="none")
    # Same anchor as v12; folds 1, 3 and 4 lie within one text height of each other, so their labels
    # sit above or below the anchor instead of centred on it.
    label_va = {1: "bottom", 3: "top", 4: "bottom"}
    for x, (_, r) in zip(fx, folds.iterrows()):
        ax.text(x + 0.03, r.cos, f"{r.cos:.3f}", va=label_va.get(int(r.fold), "center"), ha="left", fontsize=FS)
        data.append(dict(panel="D", series="GTEx fold", label=f"fold {int(r.fold)}", x=float(x), y=float(r.cos),
                         n_cells=np.nan, note=f"n_train={int(r.n_train)}; n_test={int(r.n_test)}"))
    ax.axhline(bar, color=DARKGREY, ls="--", lw=LW_AXIS)
    ax.text(1.48, bar - 0.03, f"pre-registered bar {bar:g}", ha="right", va="top", fontsize=FS, color=DARKGREY)
    data.append(dict(panel="D", series="pre-registered bar (dashed line)", label="c3 bar", x=np.nan, y=bar,
                     n_cells=np.nan, note="results/toward/summary.json c3.bar"))
    ax.set_xlim(0.5, 1.5)
    ax.set_ylim(0, 1.0)
    ax.set_xticks([])
    ax.set_title("Positive control within GTEx (5 folds)")
    ax.set_ylabel("Fold cosine")
    panel_letters(fig, [[(axs[0, 0], "A"), (axs[0, 1], "B")], [(axs[1, 0], "C"), (axs[1, 1], "D")]])
    save(fig, "FigS1", data)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig1()
    fig2()
    d = _load_manuscript()
    fig3(d)
    fig4(d)
    fig5()
    figS1()


if __name__ == "__main__":
    main()

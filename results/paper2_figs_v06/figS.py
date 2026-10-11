"""Supplementary Figures S2-S4 (v0.9; S1 is drawn by fig12.py). s4() is kept for the record but not used. usage: python figS.py <run3 dir> <run6 dir> <out dir>"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from common import *  # noqa: E402,F403
from matplotlib.lines import Line2D  # noqa: E402

RUN3, RUN6, OUT = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
ROWS = [("A3", "Passage", ORANGE), ("A1", "Donor age", BLUE),
        ("A4", "Donor pair (22 vs 96 y)", DARK), ("A2", "Frozen age ruler", DARK)]


RULER_CATS = [("called", ORANGE), ("younger only", DARK), ("neither", LIGHT)]


def ruler_panel(ax, o, letter, title):
    R = rows(RUN3 / "FigS1_data.csv")
    P = [r for r in R if r["origin"] == o and r["perturbation"]]
    for name, col in RULER_CATS[::-1]:
        if name == "called":
            s = [r for r in P if r["called"] == "True"]
        elif name == "younger only":
            s = [r for r in P if r["called"] != "True" and r["ruler_scores_younger"] == "True"]
        else:
            s = [r for r in P if r["called"] != "True" and r["ruler_scores_younger"] != "True"]
        ax.scatter([f(r["delta"]) for r in s], [f(r["ruler_change"]) for r in s], s=8, color=col, lw=0)
    for r in P:
        if r["perturbation"] in FOUR:
            ax.plot(f(r["delta"]), f(r["ruler_change"]), "o", mfc="white", mec=BLACK, ms=5, mew=1.0, zorder=5)
    ax.axhline(0, color=BLACK, lw=0.7); ax.axvline(0, color=BLACK, lw=0.7)
    ax.set_xlabel("Change in distance to early-passage cells")
    ax.set_ylabel("Change in age-ruler score\n(negative = scores younger)")
    n = {r["item"]: int(f(r["value"])) for r in R if r["item"] and r["origin"] == o}
    tag(ax, letter, title.format(k=n["n_younger_not_called"], n=n["n_younger"]))


RULER_LEGEND = [Line2D([], [], ls="", marker="o", color=ORANGE, ms=5, label="called closer to early passage"),
                Line2D([], [], ls="", marker="o", color=DARK, ms=5, label="ruler says younger, not called"),
                Line2D([], [], ls="", marker="o", color=LIGHT, ms=5, label="neither"),
                Line2D([], [], ls="", marker="o", mfc="white", mec=BLACK, ms=5, label="the four validated factors")]


def s5():
    """Directional age score against distance, own-control start (was S1A in v0.7)."""
    fig = plt.figure(figsize=(W, 2.9))
    ax = fig.add_axes([0.11, 0.18, 0.36, 0.64])
    ruler_panel(ax, "i", "", "{k} of {n} perturbations that score younger\ndo not move closer to early-passage cells")
    ax.legend(handles=RULER_LEGEND, loc="upper left", bbox_to_anchor=(1.04, 1.0), handletextpad=0.3)
    save(fig, OUT, "FigS4")


def s2():
    R = rows(RUN6 / "FigS2_data.csv")
    fig = plt.figure(figsize=(W, 4.4))
    axA = fig.add_axes([0.20, 0.62, 0.28, 0.29]); axB = fig.add_axes([0.70, 0.62, 0.28, 0.29])
    axC = fig.add_axes([0.10, 0.10, 0.36, 0.33]); axD = fig.add_axes([0.62, 0.10, 0.36, 0.33])
    order = [("S2A_k20", "A"), ("S2A_k20", "I"), ("S1", "A"), ("S1", "I")]
    labels = ["Activation screen,\n20 components", "Interference screen,\n20 components",
              "Activation screen,\nfull gene space", "Interference screen,\nfull gene space"]
    for ax, letter, title, xl, ref in (
            (axA, "A", "Offset between screen controls\nand late-passage cells", "Offset, as a multiple of the passage gap", 1.0),
            (axB, "B", "Direction of that offset", "Cosine with the passage direction", 0.0)):
        for r in R:
            if r["panel"] != letter:
                continue
            y = order.index((r["space"], r["library"]))
            full = r["space"] == "S1"
            ax.plot([f(r["lo_2_5"]), f(r["hi_97_5"])], [y, y], color=MID if full else BLACK, lw=1.0)
            ax.plot(f(r["thinned_median"]), y, "o", ms=5, mew=1.0, mec=MID if full else BLACK,
                    mfc="white" if full else BLACK)
        ax.axvline(ref, color=MID, ls="--", lw=0.8)
        ax.set_yticks(range(4)); ax.set_yticklabels(labels, fontsize=FS - 1); ax.set_ylim(3.5, -0.5)
        ax.set_xlabel(xl)
        tag(ax, letter, title)
    v = [f(r["n_cells"]) for r in R if r["panel"] == "C" and r["modality"] in ("CRA", "CRI")]
    axC.hist(v, bins=range(30, 215, 10), color=LIGHT, ec="white", lw=0.6)
    top = axC.get_ylim()[1]
    for r in R:
        if r["panel"] == "C" and r["perturbation"] in FOUR:
            n = f(r["n_cells"])
            axC.plot(n, top * 1.02, "v", color=BLACK, ms=5, clip_on=False)
            axC.text(n, top * 1.08, NAME[r["perturbation"]], rotation=90, fontsize=FS - 1.5, ha="center", va="bottom")
    axC.set_ylim(0, top * 1.35)
    axC.set_xlim(0, 220)
    axC.set_xlabel("Cells per perturbation (minimum 30)"); axC.set_ylabel(f"Perturbations (of {len(v)})")
    tag(axC, "C", "The four factors are among the largest groups")
    D = [r for r in R if r["panel"] == "D"]
    for i, r in enumerate(D):
        p, a = f(r["present"]), f(r["absent"])
        axD.barh(i, p, color=LIGHT, height=0.6)
        axD.barh(i, a, left=p, color="white", ec=DARK, lw=0.7, height=0.6)
        axD.text(300, i, f"{int(p):,} measured; {f(r['pct_gap_sq_unmeasured']):.1f}% of the\ndonor-age difference unmeasured",
                 va="center", fontsize=FS - 1.5)
    axD.set_yticks(range(len(D))); axD.set_yticklabels([r["dataset"].replace(" (corrected)", "") for r in D])
    axD.set_ylim(len(D) - 0.5, -0.5); axD.set_xlim(0, 23485)
    axD.set_xticks([0, 10000, 20000]); axD.set_xticklabels(["0", "10,000", "20,000"])
    axD.set_xlabel("Reference genes (of 23,485) measured")
    tag(axD, "D", "Reference genes measured per dataset")
    save(fig, OUT, "FigS2")


def s4():
    R = rows(RUN3 / "FigS3_data.csv")
    fig = plt.figure(figsize=(W, 2.3))
    ax = fig.add_axes([0.20, 0.25, 0.42, 0.62])
    src = {"C2": ("Cell-cycle gene list (94 genes)", "^"), "C1_gtex": ("Bulk GTEx donor samples (652)", "v")}
    order = [k for k, _, _ in ROWS]
    keys = list(src)
    for r in R:
        i = order.index(r["aging_axis"]); col = ROWS[i][2]
        dy = (keys.index(r["prolif_axis"]) - 0.5) * 0.18
        ax.plot(f(r["angle_from_parallel"]), i + dy, src[r["prolif_axis"]][1], ms=5.5, mec=col, mfc=col)
    ax.set_yticks(range(4)); ax.set_yticklabels([lab for _, lab, _ in ROWS])
    for t, (_, _, col) in zip(ax.get_yticklabels(), ROWS):
        t.set_color(col)
    ax.set_ylim(3.5, -0.5); ax.set_xlim(0, 90); ax.set_xticks([0, 15, 30, 45, 60, 75, 90])
    ax.set_xticklabels(["0\nsame", "15", "30", "45", "60", "75", "90\nunrelated"])
    ax.set_xlabel("Angle to the cell-division direction (degrees)")
    h = [Line2D([], [], ls="", marker=m, mec=DARK, mfc=DARK, ms=5.5, label=lab) for lab, m in src.values()]
    ax.legend(handles=h, title="Cell-division direction from", title_fontsize=FS - 0.5, loc="upper left",
              bbox_to_anchor=(1.02, 1.0), alignment="left", handletextpad=0.3)
    tag(ax, "", "Directions not built from single cells")
    save(fig, OUT, "FigS4_unused")


def s3():
    """Every main result from the late-passage start (was S4 and S1B in v0.7)."""
    R2 = rows(RUN3 / "Fig2_data.csv"); R3 = rows(RUN6 / "Fig3_data.csv")
    fig = plt.figure(figsize=(W, 5.6))
    axA = fig.add_axes([0.10, 0.62, 0.22, 0.30]); axB = fig.add_axes([0.42, 0.62, 0.20, 0.30])
    axC = fig.add_axes([0.72, 0.62, 0.24, 0.30]); axC.axis("off")
    axD = fig.add_axes([0.10, 0.10, 0.28, 0.30]); axE = fig.add_axes([0.50, 0.10, 0.26, 0.30])
    sty = {"S2A_k20": (BLACK, BLACK), "S1": (MID, "white")}
    depths = [96, 107, 128, 157]
    off = {("S1", "CRA"): -0.27, ("S1", "CRI"): -0.12, ("S2A_k20", "CRA"): 0.12, ("S2A_k20", "CRI"): 0.27}
    for r in R2:
        if r["panel"] == "A" and r["origin"] == "ii":
            x = depths.index(int(f(r["n_cells"]))) + off[(r["space"], r["library"])]
            ec, fc = sty[r["space"]]
            axA.plot([x, x], [f(r["lo_2_5"]), f(r["hi_97_5"])], color=ec, lw=0.9)
            axA.plot(x, f(r["null_delta_median"]), "o", ms=4, mew=1.0, mec=ec, mfc=fc)
    axA.axhline(0, color=BLACK, lw=0.7); axA.set_xticks(range(4)); axA.set_xticklabels([str(d) for d in depths])
    axA.set_xlabel("Cells per control sample"); axA.set_ylabel("Change in distance to\nearly-passage cells")
    tag(axA, "A", "No-effect samples")
    for space in ("S1", "S2A_k20"):
        for i, dz in enumerate([0.25, 0.5, 1.0]):
            ys = sorted(f(r["power"]) for r in R2 if r["panel"] == "B" and r["drawn"] == "True" and r["origin"] == "ii"
                        and r["space"] == space and f(r["dose"]) == dz)
            ec, fc = sty[space]
            for j, yv in enumerate(ys):
                axB.plot(i + (-0.14 if space == "S1" else 0.14) + jitter(j, len(ys), 0.16), yv, "o", ms=3.5,
                         mew=1.0, mec=ec, mfc=fc)
    axB.axhline(0.8, color=MID, ls="--", lw=0.8); axB.set_ylim(-0.04, 1.06)
    axB.set_xticks(range(3)); axB.set_xticklabels(["25%", "50%", "100%"]); axB.set_xlim(-0.5, 2.5)
    axB.set_xlabel("Planted effect"); axB.set_ylabel("Fraction detected")
    tag(axB, "B", "Planted effects")
    calls = {(r["column"], r["factor"]): r for r in R2 if r["panel"] == "C"}
    cols = [("Full-gene (ii)", 0.44, "full\ngene"), ("S2A_k20 (ii)", 0.86, "20\ncomponents")]
    for c, x, lab in cols:
        axC.text(x, 0.82, lab, ha="center", va="bottom", fontsize=FS - 1, transform=axC.transAxes)
    for y, k in zip([0.68, 0.54, 0.40, 0.26], FOUR):
        axC.text(0.28, y, NAME[k], ha="right", va="center", transform=axC.transAxes)
        for c, x, _ in cols:
            on = calls[(c, k)]["called"] == "True"
            axC.plot(x, y, "s", ms=8, mew=1.0, mec=ORANGE if on else MID, mfc=ORANGE if on else "white",
                     transform=axC.transAxes, clip_on=False)
    axC.text(0.28, 0.10, "of 340", ha="right", va="center", fontsize=FS - 0.5, transform=axC.transAxes)
    for c, x, _ in cols:
        axC.text(x, 0.10, str(int(f(calls[(c, "all")]["n_called"]))), ha="center", va="center", transform=axC.transAxes)
    tag(axC, "C", "Factors called")
    cats = ["NONE", "P1"]
    for i, k in enumerate(cats):
        for lib, mk, dx in (("CRA", "s", -0.15), ("CRI", "^", 0.15)):
            ys = sorted(f(r["power_50"]) for r in R3 if r["panel"] == "B" and r["axis_removed"] == k
                        and r["origin"] == "ii" and r["library"] == lib)
            for j, yv in enumerate(ys):
                axD.plot(i + dx + jitter(j, len(ys), 0.2), yv, mk, ms=3.5, mec=BLACK, mfc=BLACK)
    axD.axhline(0.8, color=MID, ls="--", lw=0.8); axD.set_ylim(-0.04, 1.06); axD.set_xlim(-0.6, 1.6)
    axD.set_xticks(range(2)); axD.set_xticklabels(["nothing\nremoved", "cell division removed\n(own screen's control cells)"], fontsize=FS - 1)
    axD.set_ylabel("Fraction detected")
    tag(axD, "D", "After removing cell division")
    h = [Line2D([], [], ls="", marker="o", mec=BLACK, mfc=BLACK, ms=5, label="20-component space"),
         Line2D([], [], ls="", marker="o", mec=MID, mfc="white", mew=1.0, ms=5, label="full gene space"),
         Line2D([], [], ls="", marker="s", mec=ORANGE, mfc=ORANGE, ms=6, label="called")]
    fig.legend(handles=h, loc="lower left", bbox_to_anchor=(0.06, 0.475), ncol=3, handletextpad=0.3, columnspacing=1.6)
    ruler_panel(axE, "ii", "E", "Directional score: {k} of {n} scoring\nyounger do not move closer")
    axE.legend(handles=RULER_LEGEND, loc="upper left", bbox_to_anchor=(1.03, 1.0), handletextpad=0.3, fontsize=FS - 1)
    save(fig, OUT, "FigS3")


for fn in (s2, s3, s5):
    fn()
print("S2-S4 written")

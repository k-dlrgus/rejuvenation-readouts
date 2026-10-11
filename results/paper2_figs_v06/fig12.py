"""Figure 1 and Supplementary Figure S1 (v0.8; S1 was Figure 2 in v0.7). usage: python fig12.py <run3 dir> <out dir>"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from common import *  # noqa: E402,F403
from matplotlib.patches import FancyBboxPatch  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

RUN3, OUT = Path(sys.argv[1]), sys.argv[2]


def fig1():
    d = {r["item"]: int(f(r["value"])) for r in rows(RUN3 / "Fig1_data.csv")}
    fig = plt.figure(figsize=(W, 2.85))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 42); ax.axis("off")

    def box(x, y, w, h, title, lines, ec):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.2,rounding_size=0.8",
                                    ec=ec, fc="white", lw=1.0))
        ax.text(x + 1.2, y + h - 1.4, title, ha="left", va="top", fontweight="bold", fontsize=FS)
        ax.text(x + 1.2, y + h - 5.2, "\n".join(lines), ha="left", va="top", fontsize=FS - 0.5, linespacing=1.5)

    ax.text(1, 41.3, f"Sengstack et al. Perturb-seq data: one neonatal dermal fibroblast line, {d['all cells']:,} cells",
            ha="left", va="top", fontsize=FS + 0.5, fontweight="bold")
    box(1, 23.5, 30, 13, f"Wild-type cells ({d['WT group']:,})",
        [f"Early passage, PD14   {d['PD14']:,}", f"Middle passage, PD26   {d['PD26']:,}",
         f"Late passage, PD32   {d['PD32']:,}"], DARK)
    box(35, 23.5, 30, 13, f"CRISPR activation screen ({d['CRA group']:,})",
        ["200 transcription factors switched on", f"  {d['perturbations CRA']} with at least 30 cells, analyzed", f"{d['NT CRA']:,} non-targeting control cells"], DARK)
    box(69, 23.5, 30, 13, f"CRISPR interference screen ({d['CRI group']:,})",
        ["200 transcription factors switched off", f"  {d['perturbations CRI']} with at least 30 cells, analyzed", f"{d['NT CRI']:,} non-targeting control cells"], DARK)
    box(1, 1.2, 47, 9.5, "Target 1: early-passage cells",
        ["Same experiment. Young = PD14, old = PD32.", "Figures 2 and 3"], ORANGE)
    box(52, 1.2, 47, 9.5, "Target 2: young human donors",
        [f"GTEx fibroblasts: {d['GTEx young 20-39']} donors aged 20–39 (young),",
         f"{d['GTEx old 60-79']} aged 60–79 (old). Figure 2"], BLUE)

    arr = dict(arrowstyle="-|>", lw=1.0, mutation_scale=9)
    ax.annotate("", xy=(12, 11.4), xytext=(12, 23.1), arrowprops=dict(color=DARK, **arr))
    ax.text(13.5, 17.3, "defines young\nand old", fontsize=FS - 1, color=DARK, va="center")
    for x in (50, 84):
        ax.plot([x, x], [23.1, 19.5], color=DARK, lw=1.0)
    ax.plot([30, 84], [19.5, 19.5], color=DARK, lw=1.0)
    ax.annotate("", xy=(30, 11.4), xytext=(30, 19.5), arrowprops=dict(color=ORANGE, **arr))
    ax.annotate("", xy=(75, 11.4), xytext=(75, 19.5), arrowprops=dict(color=BLUE, **arr))
    ax.text(52.5, 18.6, "each perturbation's effect, relative to the control cells of its screen,\n"
            "is tested against both targets", fontsize=FS - 1, color=DARK, va="top", ha="center", linespacing=1.3)
    save(fig, OUT, "Fig1")


def fig2():
    R = rows(RUN3 / "Fig2_data.csv")
    fig = plt.figure(figsize=(W, 3.0))
    axA = fig.add_axes([0.085, 0.25, 0.25, 0.60])
    axB = fig.add_axes([0.43, 0.25, 0.22, 0.60])
    axC = fig.add_axes([0.73, 0.25, 0.26, 0.60]); axC.axis("off")
    sty = {"S2A_k20": dict(color=BLACK, mfc=BLACK), "S1": dict(color=MID, mfc="white")}

    depths = [96, 107, 128, 157]
    off = {("S1", "CRA"): -0.27, ("S1", "CRI"): -0.12, ("S2A_k20", "CRA"): 0.12, ("S2A_k20", "CRI"): 0.27}
    for r in R:
        if r["panel"] != "A" or r["origin"] != "i":
            continue
        x = depths.index(int(f(r["n_cells"]))) + off[(r["space"], r["library"])]
        s = sty[r["space"]]
        axA.plot([x, x], [f(r["lo_2_5"]), f(r["hi_97_5"])], color=s["color"], lw=0.9)
        axA.plot(x, f(r["null_delta_median"]), "o", ms=4.2, mew=1.0, mec=s["color"], mfc=s["mfc"])
    axA.axhline(0, color=BLACK, lw=0.7)
    axA.set_xticks(range(4)); axA.set_xticklabels([str(v) for v in depths]); axA.set_xlim(-0.6, 3.6)
    axA.set_xlabel("Cells per control sample")
    axA.set_ylabel("Change in distance to\nearly-passage cells")
    axA.text(3.55, 52, "full gene space:\nalways farther", color=DARK, fontsize=FS - 1, ha="right", va="top")
    axA.text(3.55, -12, "20-component space:\ncentred on zero", color=BLACK, fontsize=FS - 1, ha="right", va="top")
    axA.set_ylim(-24, 54)
    tag(axA, "A", "Control samples with no effect")

    doses = [0.25, 0.5, 1.0]
    for space in ("S1", "S2A_k20"):
        pts = [r for r in R if r["panel"] == "B" and r["drawn"] == "True" and r["origin"] == "i" and r["space"] == space]
        for i, dz in enumerate(doses):
            ys = sorted(f(r["power"]) for r in pts if f(r["dose"]) == dz)
            base = i + (-0.14 if space == "S1" else 0.14)
            for j, yv in enumerate(ys):
                axB.plot(base + jitter(j, len(ys), 0.16), yv, "o", ms=3.8, mew=1.0,
                         mec=sty[space]["color"], mfc=sty[space]["mfc"])
    axB.axhline(0.8, color=MID, ls="--", lw=0.8)
    axB.set_xticks(range(3)); axB.set_xticklabels(["25%", "50%", "100%"]); axB.set_xlim(-0.5, 2.5)
    axB.set_ylim(-0.04, 1.06)
    axB.set_xlabel("Size of planted effect\n(share of the passage gap)")
    axB.set_ylabel("Fraction detected")
    tag(axB, "B", "Planted effect of known size")

    calls = {(r["column"], r["factor"]): r for r in R if r["panel"] == "C"}
    cols = [("Full-gene (i)", 0.48, "full gene\nspace"), ("S2A_k20 (i)", 0.80, "20-component\nspace")]
    for c, x, lab in cols:
        axC.text(x, 0.80, lab, ha="center", va="bottom", fontsize=FS - 0.5, transform=axC.transAxes)
    ys = [0.66, 0.53, 0.40, 0.27]
    for y, fac in zip(ys, FOUR):
        axC.text(0.30, y, NAME[fac], ha="right", va="center", fontsize=FS, transform=axC.transAxes)
        for c, x, _ in cols:
            on = calls[(c, fac)]["called"] == "True"
            axC.plot(x, y, "s", ms=9, mew=1.0, mec=ORANGE if on else MID, mfc=ORANGE if on else "white",
                     transform=axC.transAxes, clip_on=False)
    axC.text(0.30, 0.12, "Called, of 340", ha="right", va="center", fontsize=FS - 0.5, transform=axC.transAxes)
    for c, x, _ in cols:
        axC.text(x, 0.12, str(int(f(calls[(c, "all")]["n_called"]))), ha="center", va="center",
                 fontsize=FS, transform=axC.transAxes)
    axC.text(0.30, 0.01, "Expected by chance", ha="right", va="center", fontsize=FS - 0.5, color=DARK,
             transform=axC.transAxes)
    axC.text(0.64, 0.01, f"{f(calls[(cols[0][0], 'all')]['chance_expected']):.0f}", ha="center", va="center",
             fontsize=FS, color=DARK, transform=axC.transAxes)
    tag(axC, "C", "Validated factors called")

    h = [Line2D([], [], ls="", marker="o", mec=BLACK, mfc=BLACK, ms=5, label="20-component space"),
         Line2D([], [], ls="", marker="o", mec=MID, mfc="white", mew=1.0, ms=5, label="full gene space"),
         Line2D([], [], ls="", marker="s", mec=ORANGE, mfc=ORANGE, ms=6, label="called"),
         Line2D([], [], ls="", marker="s", mec=MID, mfc="white", ms=6, label="not called")]
    fig.legend(handles=h, loc="lower left", bbox_to_anchor=(0.06, -0.01), ncol=4, handletextpad=0.3, columnspacing=1.6)
    save(fig, OUT, "FigS1")


fig1()
fig2()
print("Fig1, FigS1 written")

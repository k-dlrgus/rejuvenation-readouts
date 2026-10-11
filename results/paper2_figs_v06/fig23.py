"""Figures 2 and 3 (v0.8; Figure 2 was Figure 4 in v0.7).
usage: python fig34.py <run3 dir> <run6 dir> <seng3 d4_table.csv> <axang d3b_table.csv> <out dir>"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from common import *  # noqa: E402,F403
from matplotlib.lines import Line2D  # noqa: E402

RUN3, RUN6, S3TAB, AXTAB, OUT = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5]

ROWS = [("A3", "Passage\n(Sengstack, PD32 vs PD14)", ORANGE), ("A1", "Donor age\n(GTEx, 60–79 vs 20–39 y)", BLUE),
        ("A4", "Donor pair\n(Lu, 96 vs 22 y)", DARK), ("A2", "Age ruler\n(trained on GTEx)", DARK)]
CELLS = {"C1_south": ("Southard control cells", "o"),
         "C1_seng_CRA": ("Sengstack activation-screen controls", "s"),
         "C1_seng_CRI": ("Sengstack interference-screen controls", "^"),
         "C1_lu": ("Lu day-0 fibroblasts", "D"),
         "C2": ("Cell-cycle gene list (no cells)", "*"),
         "C1_gtex": ("GTEx bulk tissue samples", "v")}


def angle_panel(ax, data, sources):
    order = [k for k, _, _ in ROWS]
    keys = list(sources)
    for r in data:
        if r["prolif_axis"] not in sources:
            continue
        i = order.index(r["aging_axis"])
        col = ROWS[i][2]
        dy = (keys.index(r["prolif_axis"]) - (len(keys) - 1) / 2) * 0.13
        hollow = r["prolif_axis"] in ("C2", "C1_gtex")
        ax.plot(f(r["angle_from_parallel"]), i + dy, sources[r["prolif_axis"]][1], ms=5.5 if hollow else 5,
                mec=col, mfc="white" if hollow else col, mew=1.1)
    ax.set_yticks(range(len(ROWS)))
    ax.set_yticklabels([lab for _, lab, _ in ROWS])
    for t, (_, _, col) in zip(ax.get_yticklabels(), ROWS):
        t.set_color(col)
    ax.set_ylim(len(ROWS) - 0.5, -0.5)
    ax.set_xlim(0, 90); ax.set_xticks([0, 15, 30, 45, 60, 75, 90])
    ax.set_xticklabels(["0\nsame", "15", "30", "45", "60", "75", "90\nunrelated"])
    ax.set_xlabel("Angle to the cell-division direction (degrees)")
    h = [Line2D([], [], ls="", marker=m, mec=DARK, mfc="white" if k in ("C2", "C1_gtex") else DARK, mew=1.1, ms=5, label=lab)
         for k, (lab, m) in sources.items()]
    ax.legend(handles=h, title="Cell-division direction from", title_fontsize=FS - 0.5, loc="upper left",
              bbox_to_anchor=(1.02, 1.0), alignment="left", handletextpad=0.3)


def fig3():
    R = rows(RUN6 / "Fig3_data.csv")
    fig = plt.figure(figsize=(W, 5.55))
    axA = fig.add_axes([0.23, 0.67, 0.40, 0.27])
    axB = fig.add_axes([0.09, 0.12, 0.36, 0.36])
    axC = fig.add_axes([0.60, 0.12, 0.38, 0.36])

    gl = [r for r in rows(RUN3 / "FigS3_data.csv") if r["prolif_axis"] in ("C2", "C1_gtex")]
    angle_panel(axA, [r for r in R if r["panel"] == "A"] + gl, CELLS)
    tag(axA, "A", "Passage lies closest to cell division measured in single cells")

    cats = [("NONE", "Nothing removed"), ("P1", "Cell division removed\n(own screen's control-cell direction,\n19.7° / 20.6° from passage)")]
    for i, (key, _) in enumerate(cats):
        for lib, mk, dx in (("CRA", "s", -0.15), ("CRI", "^", 0.15)):
            ys = sorted(f(r["power_50"]) for r in R if r["panel"] == "B" and r["axis_removed"] == key
                        and r["origin"] == "i" and r["library"] == lib)
            for j, yv in enumerate(ys):
                axB.plot(i + dx + jitter(j, len(ys), 0.2), yv, mk, ms=3.8, mec=BLACK, mfc=BLACK)
    axB.axhline(0.8, color=MID, ls="--", lw=0.8)
    axB.set_xticks(range(len(cats))); axB.set_xticklabels([lab for _, lab in cats], fontsize=FS - 1)
    axB.set_xlim(-0.6, len(cats) - 0.4); axB.set_ylim(-0.04, 1.06)
    axB.set_ylabel("Fraction of planted effects detected")
    axB.legend(handles=[Line2D([], [], ls="", marker="s", color=BLACK, ms=4, label="activation screen"),
                        Line2D([], [], ls="", marker="^", color=BLACK, ms=4, label="interference screen")],
               loc="center", bbox_to_anchor=(0.5, 0.5), handletextpad=0.2, fontsize=FS - 1)
    tag(axB, "B", "Removing cell division removes\nthe passage signal")

    C = [r for r in R if r["panel"] == "C" and r["perturbation"]]
    for called, col in ((False, MID), (True, ORANGE)):
        g = [r for r in C if (r["called_delta_i"] == "True") == called]
        axC.scatter([f(r["prolif_shift"]) for r in g], [f(r["delta_i"]) for r in g], s=8, color=col, lw=0, zorder=2)
    nC = sum(r["called_delta_i"] == "True" for r in C)
    axC.legend(handles=[Line2D([], [], ls="", marker="o", color=ORANGE, ms=4, label=f"called closer ({nC})"),
                        Line2D([], [], ls="", marker="o", color=MID, ms=4, label=f"not called ({len(C) - nC})")],
               loc="upper right", bbox_to_anchor=(1.0, 0.88), handletextpad=0.2, fontsize=FS - 1)
    pos = {r["perturbation"]: (f(r["prolif_shift"]), f(r["delta_i"])) for r in C}
    lab = {"CRA_EZH2": (0.17, 4.0), "CRA_E2F3": (0.22, -3.6), "CRI_STAT3": (-0.14, -15.0), "CRI_ZFX": (-0.07, -21.0)}
    for k in FOUR:
        x, y = pos[k]
        axC.plot(x, y, "o", mfc="white", mec=BLACK, ms=5.5, mew=1.1, zorder=5)
        axC.annotate(NAME[k], (x, y), xytext=lab[k], fontsize=FS - 0.5, va="center",
                     arrowprops=dict(arrowstyle="-", lw=0.6, color=BLACK, shrinkA=1, shrinkB=3), zorder=6)
    rho = [r for r in R if r["panel"] == "C" and r["item"]][0]
    axC.axhline(0, color=BLACK, lw=0.7); axC.axvline(0, color=BLACK, lw=0.7)
    axC.text(0.97, 0.96, f"Spearman ρ = {f(rho['value']):.2f}", transform=axC.transAxes, ha="right", va="top")
    axC.set_xlabel("Cell-cycle gene score, relative to controls\n(higher = more dividing)")
    axC.set_ylabel("Change in distance to early-passage\ncells (negative = closer)")
    tag(axC, "C", "Perturbations that raise cell-cycle genes\nmove toward early passage")
    save(fig, OUT, "Fig3")


def fig4():
    s3 = {r["perturbation"]: r for r in rows(S3TAB)}
    ax_ = {r["perturbation"]: r for r in rows(AXTAB)}
    pts = []
    for p in s3:
        a, b = s3[p], ax_[p]
        pts.append((p, f(a["delta_i"]) - f(a["floor_delta_i"]), f(b["delta"]) - f(b["floor"]),
                    a["called_delta_i"] == "True", b["called"] == "True"))
    fig = plt.figure(figsize=(W, 3.3))
    ax = fig.add_axes([0.10, 0.15, 0.50, 0.78])
    groups = [("neither", lambda e, d: not e and not d, MID, "o", 9),
              ("early", lambda e, d: e and not d, ORANGE, "o", 11),
              ("donor", lambda e, d: d and not e, BLUE, "s", 15),
              ("both", lambda e, d: e and d, BLACK, "o", 11)]
    n = {}
    for key, test, col, mk, sz in groups:
        g = [q for q in pts if test(q[3], q[4])]
        n[key] = len(g)
        ax.scatter([q[1] for q in g], [q[2] for q in g], s=sz, color=col, marker=mk, lw=0, zorder=3)
    ax.axvline(0, color=BLACK, lw=0.8); ax.axhline(0, color=BLACK, lw=0.8)
    xl, yl = ax.get_xlim(), ax.get_ylim()
    ax.axvspan(xl[0], 0, color=ORANGE, alpha=0.06, lw=0); ax.axhspan(yl[0], 0, color=BLUE, alpha=0.06, lw=0)
    ax.set_xlim(xl); ax.set_ylim(yl)
    lab = {"CRI_ZFX": (-16.3, 20.8), "CRA_EZH2": (-5.0, 20.8), "CRA_E2F3": (-16.3, 8.0), "CRI_STAT3": (-14.0, 3.0)}
    for p, x, y, _, _ in pts:
        if p in FOUR:
            ax.plot(x, y, "o", mfc="white", mec=BLACK, ms=5.5, mew=1.1, zorder=5)
            ax.annotate(NAME[p], (x, y), xytext=lab[p], fontsize=FS - 0.5, va="center",
                        arrowprops=dict(arrowstyle="-", lw=0.6, color=BLACK, shrinkA=1, shrinkB=3), zorder=6)
    ax.set_xlabel("Change in distance to early-passage cells,\nrelative to the chance threshold")
    ax.set_ylabel("Change in distance to young donors,\nrelative to the chance threshold")
    ax.text(xl[0] + 0.5, yl[0] + 0.35, "closer to both: none", fontsize=FS - 1, color=DARK, va="bottom")
    h = [Line2D([], [], ls="", marker="o", color=ORANGE, ms=5, label=f"closer to early-passage cells only ({n['early']})"),
         Line2D([], [], ls="", marker="s", color=BLUE, ms=5, label=f"closer to young donors only ({n['donor']})"),
         Line2D([], [], ls="", marker="o", color=BLACK, ms=5, label=f"closer to both ({n['both']})"),
         Line2D([], [], ls="", marker="o", color=MID, ms=5, label=f"neither ({n['neither']})"),
         Line2D([], [], ls="", marker="o", mfc="white", mec=BLACK, ms=5.5, label="the four validated factors")]
    ax.legend(handles=h, loc="upper left", bbox_to_anchor=(1.03, 1.0), handletextpad=0.3,
              title="340 perturbations", title_fontsize=FS, alignment="left")
    ax.text(1.04, 0.38, "Zero on each axis is the chance threshold\nfor that perturbation's own group size.\n"
            "Left of the vertical line: closer to early-\npassage cells than chance allows. Below\n"
            "the horizontal line: closer to young donors.", transform=ax.transAxes, fontsize=FS - 1,
            color=DARK, va="top", linespacing=1.35)
    save(fig, OUT, "Fig2")


fig3()
fig4()
print("Fig2, Fig3 written")

"""PAPER2FIGS: figures for the second paper, drawn from stored SENG/AXANG results only.

A drawing task. No statistic is computed, nothing is refit. Every plotted value is read from a
stored file under results/senginv/, results/seng2/, results/seng3/, results/seng4/ or
results/axang/; if a value is not stored the script stops with `value_not_stored`.
The plan and the author's answers are in results/paper2_figs/PLAN.flag.

  Fig1  design schematic (numbers from results/senginv/, seng2, seng3, axang)
  Fig2  A null delta, full-gene vs S2A_k20         results/seng3/d1_null_validation.csv
        B planted-approach power, both spaces      results/seng3/d2_plant.csv
        C the four factors called or not           results/seng2/d2_table.csv, results/seng3/d4_table.csv
  Fig3  A angle from parallel, A1-A4 vs cell C1   results/axang/d1_angles.csv,
                                                   results/axang/d1_angles_corrected_C1lu.csv (C1_lu)
        B power at a 50% plant, residualized       results/seng4/d2_plant_residualized.csv
        C proliferation shift vs delta (i)         results/seng3/d4_table.csv
  Fig4  delta to early-passage vs delta to GTEx young
                                                   results/seng3/d4_table.csv x results/axang/d3b_table.csv
  FigS1 frozen-ruler change vs delta, (i) and (ii) results/seng3/d4_table.csv
  FigS2 library offsets, cells per perturbation, ruler-gene coverage
  FigS3 as Fig3A for C2 and C1_gtex

"Angle from parallel" is arccos|cos| (AXANG amendment F10). d1_angles.csv stores the signed
angle, so it is shown as min(angle, 180 - angle): a display conversion, at the author's
instruction. C1_lu is read from angle_from_parallel_corrected in the corrected CSV.

Style: the constants and helpers of src/figures_final.py (results/paper_figs_v1_2/), imported.

Also writes results/paper2_figs/paper2_numbers.md: each listed manuscript number, the stored file
and line it is read from, and MATCH or MISMATCH.

Usage: python src/paper2_figs.py
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figures_final as ff  # noqa: E402  style constants and layout helpers only
from config import RESULTS, ROOT  # noqa: E402

plt = ff.plt
from matplotlib.lines import Line2D  # noqa: E402

OUT = RESULTS / "paper2_figs"
SENGINV, SENG2, SENG3, SENG4, AXANG = (RESULTS / d for d in ("senginv", "seng2", "seng3", "seng4", "axang"))

FOUR = ("CRA_E2F3", "CRA_EZH2", "CRI_STAT3", "CRI_ZFX")
LIBS = ("CRA", "CRI")
DEPTHS = (96, 107, 128, 157)
LIB_MARKER = dict(CRA="o", CRI="s")
SPACE_WORDS = dict(S1="Full-gene", S2A_k20="S2A_k20")
SPACE_COLOR = dict(S1=ff.BLUE, S2A_k20=ff.VERMILLION)
CALLED, NOT_CALLED = ff.VERMILLION, ff.GREY
AGING_WORDS = dict(A1="A1 GTEx donor age", A2="A2 frozen ruler", A3="A3 Sengstack passage",
                   A4="A4 Lu donor pair")
PROLIF_WORDS = dict(C1_south="C1 Southard NT (14,013 cells)", C1_seng_CRA="C1 Sengstack CRA NT (333)",
                    C1_seng_CRI="C1 Sengstack CRI NT (739)", C1_lu="C1 Lu day 0 (12,803), corrected",
                    C2="C2 Tirosh gene list (94 genes)", C1_gtex="C1 GTEx donors (652, bulk)")
PROLIF_STYLE = dict(C1_south=(ff.BLACK, "o"), C1_seng_CRA=(ff.VERMILLION, "o"),
                    C1_seng_CRI=(ff.VERMILLION, "s"), C1_lu=(ff.BLUE, "D"),
                    C2=(ff.PURPLE, "^"), C1_gtex=(ff.DARKGREY, "v"))

# PAPER2FIGS-fix (src/paper2_figs_fix_run.py) sets FIX = True: Fig1 routes both libraries to both
# targets, and the Fig3C / Fig4 factor labels are placed by place_labels. False draws PAPER2FIGS as committed.
FIX = False
FIG4_MANUAL = None  # run 3: name -> (x, y) text anchor in data coordinates, fixed in its PLAN.flag
FIG3_MANUAL = None  # fix2: name -> (x, y) text anchor in data coordinates, fixed in its PLAN.flag
LABEL_FALLBACK = False  # run 3: an exhausted Fig3C search keeps the 7a1efb9 offset instead of stopping
LABEL_CHECK = {}  # filled by place_labels: panel -> list of per-label overlap results


def stop(key, msg):
    raise SystemExit(f"STOP {key}: {msg}")


def rel(p):
    return Path(p).relative_to(ROOT).as_posix()


# --------------------------------------------------------------------------- stored-value readers
def json_at(path, *keys):
    """Value at a key path of a stored JSON file, and the line it sits on."""
    obj = json.loads(Path(path).read_text(encoding="utf-8"))
    for k in keys:
        if isinstance(obj, list):
            obj = obj[k]
        elif k not in obj:
            stop("value_not_stored", f"{rel(path)} has no key path {keys}")
        else:
            obj = obj[k]
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    i = 0
    for k in keys:  # walk down the pretty-printed file one key at a time
        if isinstance(k, int):
            continue
        while i < len(lines) and f'"{k}":' not in lines[i]:
            i += 1
        if i == len(lines):
            stop("value_not_stored", f"{rel(path)}: key {k} not found on any line")
    return obj, i + 1


def csv_rows(path, **eq):
    df = pd.read_csv(path)
    df["_line"] = np.arange(len(df)) + 2  # header is line 1; no field in these files spans lines
    m = np.ones(len(df), bool)
    for k, v in eq.items():
        m &= (df[k] == v).to_numpy()
    if not m.any():
        stop("value_not_stored", f"{rel(path)} has no row {eq}")
    return df[m]


def csv_at(path, col, **eq):
    r = csv_rows(path, **eq)
    if len(r) != 1:
        stop("value_not_stored", f"{rel(path)} row {eq} is not unique ({len(r)})")
    return r[col].iloc[0], int(r["_line"].iloc[0])


def md_line(path, needle, start=1):
    for i, l in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), 1):
        if i >= start and needle in l:
            return i, l
    stop("value_not_stored", f"{rel(path)} does not contain {needle!r}")


def from_parallel(signed_deg):
    """AXANG F10: angle from parallel = arccos|cos| = min(theta, 180 - theta). Display conversion."""
    return float(min(signed_deg, 180.0 - signed_deg))


# --------------------------------------------------------------------------- shared drawing bits
def legend_handles(items):
    hs = []
    for label, color, marker, filled in items:
        hs.append(Line2D([], [], ls="none", marker=marker, ms=ff.MS, mec=color, mew=ff.MEW,
                         mfc=color if filled else "white", label=label))
    return hs


def dot(ax, x, y, color, marker, filled, zorder=3, ms=ff.MS):
    ax.plot([x], [y], ls="none", marker=marker, ms=ms, mec=color, mew=ff.MEW,
            mfc=color if filled else "white", zorder=zorder)


def tidy(ax):
    ax.tick_params(length=2.5, width=ff.LW_AXIS)


def _label_geometry(fig, ax, xy_all, r_pt, pad_pt):
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    P = ax.transData.transform(np.asarray(xy_all, float))
    rad = (r_pt + pad_pt) * fig.dpi / 72.0

    def hits_points(bb):
        cx = np.clip(P[:, 0], bb.x0, bb.x1)
        cy = np.clip(P[:, 1], bb.y0, bb.y1)
        return int(((P[:, 0] - cx) ** 2 + (P[:, 1] - cy) ** 2 < rad ** 2).sum())
    return rend, hits_points


def _text_box(t, rend):
    from matplotlib.text import Text
    t.update_positions(rend)
    return Text.get_window_extent(t, rend)  # the text box only, not the leader line


def check_labels(fig, ax, panel, annots, xy_all, avoid=(), r_pt=2.6, pad_pt=1.0, how=""):
    """Log, per label, how many point discs (radius r_pt + pad_pt), other labels and `avoid` artists its
    text box overlaps, and whether it lies inside the axes. Logged in LABEL_CHECK[panel]; never stops."""
    rend, hits_points = _label_geometry(fig, ax, xy_all, r_pt, pad_pt)
    ab = ax.get_window_extent(rend)
    boxes = [_text_box(t, rend) for t in annots]
    blocked = [a.get_window_extent(rend) for a in avoid]
    log = []
    for i, (t, bb) in enumerate(zip(annots, boxes)):
        log.append(dict(panel=panel, label=t.get_text(), placement=how[i] if isinstance(how, list) else how,
                        points_overlapped=hits_points(bb),
                        labels_overlapped=sum(bb.overlaps(o) for j, o in enumerate(boxes) if j != i),
                        avoid_overlapped=sum(bb.overlaps(o) for o in blocked),
                        inside_axes=bool(bb.x0 >= ab.x0 and bb.x1 <= ab.x1 and bb.y0 >= ab.y0 and bb.y1 <= ab.y1)))
    LABEL_CHECK[panel] = log


def place_labels(fig, ax, panel, labels, xy_all, avoid=(), r_pt=2.6, pad_pt=1.0, fallback=None, manual=None):
    """Place each label (name, (x, y)) where its text box overlaps no data point, no earlier label
    and no artist in `avoid`, at the shortest leader from a fixed candidate list (deterministic).
    Points are discs of radius r_pt (the 4 pt open highlight marker, 2 pt + edge, rounded up),
    padded by pad_pt. If no candidate clears a label: STOP, or, with `fallback` (name -> (dx, dy)
    in points), use that offset and log it. The result is checked and logged by check_labels."""
    rend, hits_points = _label_geometry(fig, ax, xy_all, r_pt, pad_pt)
    ab = ax.get_window_extent(rend)
    blocked = [a.get_window_extent(rend) for a in avoid]
    cands = [(d * np.cos(np.radians(a)), d * np.sin(np.radians(a)))
             for d in (8, 11, 14, 18, 22, 27, 33, 40, 48) for a in range(0, 360, 15)]
    lead = dict(arrowstyle="-", lw=ff.LW_AXIS * 0.8, color=ff.BLACK, shrinkA=0, shrinkB=1.5)
    placed, annots, how = [], [], []
    manual = manual or {}
    for name, (x, y) in labels:
        best = None
        if name in manual:  # fixed by hand in a PLAN.flag; checked, never searched
            t = ax.annotate(name, (x, y), xytext=manual[name], textcoords="data", ha="left", va="center",
                            zorder=5, arrowprops=lead)
            placed.append(_text_box(t, rend))
            annots.append(t)
            how.append(f"by hand, text at {manual[name]} (data)")
            continue
        for dx, dy in cands:
            t = ax.annotate(name, (x, y), xytext=(dx, dy), textcoords="offset points",
                            ha="left" if dx >= 0 else "right", va="center", zorder=5, arrowprops=lead)
            bb = _text_box(t, rend)
            if (bb.x0 >= ab.x0 and bb.x1 <= ab.x1 and bb.y0 >= ab.y0 and bb.y1 <= ab.y1
                    and hits_points(bb) == 0 and not any(bb.overlaps(o) for o in placed + blocked)):
                best = (t, bb, f"search offset ({dx:.2f}, {dy:.2f}) pt")
                break
            t.remove()
        if best is None:
            if fallback is None:
                stop("label_unplaceable", f"{panel}: no candidate offset clears {name}")
            dx, dy = fallback[name]
            t = ax.annotate(name, (x, y), xytext=(dx, dy), textcoords="offset points",
                            ha="left" if dx > 0 else "right", va="center", zorder=5, arrowprops=lead)
            best = (t, _text_box(t, rend), f"search exhausted; 7a1efb9 offset ({dx}, {dy}) pt")
        placed.append(best[1])
        annots.append(best[0])
        how.append(best[2])
    check_labels(fig, ax, panel, annots, xy_all, avoid, r_pt, pad_pt, how)


# --------------------------------------------------------------------------- Fig1
def fig1():
    ff.use_style()
    meta = SENGINV / "addendum_meta.json"
    n_all, _ = json_at(meta, "n_cells")
    n_wt, _ = json_at(meta, "n_match_barcode")
    add = SENGINV / "ADDENDUM_SENGINV.md"
    grp = {}
    for g in ("barcode `_A`", "barcode `_I`", "`PD14`", "`PD26`", "`PD32`"):
        line, txt = md_line(add, f"| {g} |")
        grp[g] = int(txt.split("|")[2])
    nt_cra, _ = json_at(SENG2 / "d1_summary.json", "sizes", "NT_CRA")
    nt_cri, _ = json_at(SENG2 / "d1_summary.json", "sizes", "NT_CRI")
    n_pert_cra, _ = json_at(SENG3 / "d4_summary.json", "summary", "delta_i", "CRA", "n")
    n_pert_cri, _ = json_at(SENG3 / "d4_summary.json", "summary", "delta_i", "CRI", "n")
    st0 = json.loads((AXANG / "stage0_summary.json").read_text(encoding="utf-8"))
    n_young = next(c["stored"] for c in st0["checks"] if c["item"] == "anchors n_young")
    n_old = next(c["stored"] for c in st0["checks"] if c["item"] == "anchors n_old")
    if grp["`PD14`"] + grp["`PD26`"] + grp["`PD32`"] != n_wt or n_wt + grp["barcode `_A`"] + grp["barcode `_I`"] != n_all:
        stop("value_not_stored", "Fig1 group sizes do not add up to the stored totals")

    fig = ff.new_figure(2.2, 1)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 7.5)
    ax.set_ylim(0, 2.2)
    ax.axis("off")
    from matplotlib.patches import FancyBboxPatch

    def box(x, y, w, h, title, lines, ec=ff.BLACK):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.06",
                                    fc="white", ec=ec, lw=ff.LW_AXIS))
        ax.text(x + 0.08, y + h - 0.08, title, ha="left", va="top", fontweight="bold")
        ax.text(x + 0.08, y + h - 0.30, "\n".join(lines), ha="left", va="top", linespacing=1.35)

    ax.text(0.05, 2.15, f"Sengstack Perturb-seq, fibroblast line, {n_all:,} cells in three barcode groups",
            ha="left", va="top")
    box(0.05, 1.22, 2.25, 0.8, f"WT, no guide ({n_wt:,} cells)",
        [f"PD14 {grp['`PD14`']:,}   early passage",
         f"PD26 {grp['`PD26`']:,}",
         f"PD32 {grp['`PD32`']:,}   late passage"])
    box(0.05, 0.08, 2.25, 0.95, f"CRA, activation ({grp['barcode `_A`']:,} cells)",
        [f"{n_pert_cra} perturbations", f"non-targeting (NT) {nt_cra:,} cells",
         "displacement = perturbation − own NT"], ec=ff.VERMILLION)
    box(2.45, 0.08, 2.25, 0.95, f"CRI, interference ({grp['barcode `_I`']:,} cells)",
        [f"{n_pert_cri} perturbations", f"non-targeting (NT) {nt_cri:,} cells",
         "displacement = perturbation − own NT"], ec=ff.VERMILLION)
    box(5.05, 1.22, 2.4, 0.8, "Target 1: early-passage cells",
        ["PD14 (half A), same line and experiment", "aging axis a = PD14 − PD32 (half A)",
         "Figs 2–4: change in distance to PD14"])
    box(5.05, 0.08, 2.4, 0.95, "Target 2: young human donors",
        [f"GTEx donors aged 20–39 (n = {n_young})", f"origin: GTEx donors 60–79 (n = {n_old})",
         "displacement applied to the old centroid", "Fig 4: change in distance to young"])
    kw = dict(arrowstyle="-|>", lw=ff.LW_AXIS, color=ff.BLACK, mutation_scale=7)
    ax.annotate("", (5.02, 1.75), (2.33, 1.75), arrowprops=kw)
    ax.text(3.65, 1.80, "PD14 and PD32 define the target and the aging axis", ha="center", va="bottom",
            color=ff.DARKGREY)
    if not FIX:
        ax.annotate("", (5.02, 1.38), (4.73, 0.85), arrowprops=kw)
        ax.annotate("", (5.02, 0.55), (4.73, 0.55), arrowprops=kw)
    else:
        # Both libraries feed one bus at x = 4.88 that leads to both targets.
        ln = dict(color=ff.BLACK, lw=ff.LW_AXIS, solid_capstyle="butt", solid_joinstyle="miter")
        bus_x, y_cra, y_cri, y_t1, y_t2 = 4.88, 1.125, 0.55, 1.45, 0.55
        ax.plot([1.175, 1.175, bus_x], [1.03, y_cra, y_cra], **ln)        # CRA: up from its box, along the gap
        ax.plot([4.73, bus_x], [y_cri, y_cri], **ln)                     # CRI: out of its right edge
        ax.plot([bus_x, bus_x], [y_t2, y_t1], **ln)                      # the bus
        ax.annotate("", (5.02, y_t1), (bus_x, y_t1), arrowprops=kw)      # to Target 1
        ax.annotate("", (5.02, y_t2), (bus_x, y_t2), arrowprops=kw)      # to Target 2
        cap = ax.text(3.60, y_cra + 0.02, "displacements go to both targets", ha="center",
                      va="bottom", color=ff.DARKGREY)
        # G3: boxes as drawn (pad 0.02), connector segments and the caption must not meet
        boxes = dict(WT=(0.05, 1.22, 2.25, 0.8), CRA=(0.05, 0.08, 2.25, 0.95), CRI=(2.45, 0.08, 2.25, 0.95),
                     T1=(5.05, 1.22, 2.4, 0.8), T2=(5.05, 0.08, 2.4, 0.95))
        boxes = {k: (x - 0.02, y - 0.02, x + w + 0.02, y + h + 0.02) for k, (x, y, w, h) in boxes.items()}
        segs = dict(CRA_up=((1.175, 1.03), (1.175, y_cra)), CRA_run=((1.175, y_cra), (bus_x, y_cra)),
                    CRI_out=((4.73, y_cri), (bus_x, y_cri)), bus=((bus_x, y_t2), (bus_x, y_t1)),
                    to_T1=((bus_x, y_t1), (5.02, y_t1)), to_T2=((bus_x, y_t2), (5.02, y_t2)))
        allowed = dict(CRA_up={"CRA"}, CRI_out={"CRI"}, to_T1={"T1"}, to_T2={"T2"})  # endpoints at own box edge

        def meets(seg, b, eps=1e-9):
            (x0, y0), (x1, y1) = seg
            return min(x0, x1) < b[2] - eps and max(x0, x1) > b[0] + eps and                 min(y0, y1) < b[3] - eps and max(y0, y1) > b[1] + eps

        fig.canvas.draw()
        cb = cap.get_window_extent().transformed(ax.transData.inverted())
        cbox = (cb.x0, cb.y0, cb.x1, cb.y1)
        log = []
        for k, sg in segs.items():
            hit = [b for b, e in boxes.items() if meets(sg, e) and b not in allowed.get(k, set())]
            log.append(dict(panel="Fig1", item=f"connector {k}", crosses=";".join(hit) or "none"))
        hit = [b for b, e in boxes.items() if cbox[0] < e[2] and cbox[2] > e[0] and cbox[1] < e[3] and cbox[3] > e[1]]
        hit += [k for k, ((x0, y0), (x1, y1)) in segs.items()
                if min(x0, x1) < cbox[2] and max(x0, x1) > cbox[0] and min(y0, y1) < cbox[3] and max(y0, y1) > cbox[1]]
        log.append(dict(panel="Fig1", item="caption 'displacements go to both targets'",
                        crosses=";".join(hit) or "none",
                        extent=f"x {cbox[0]:.3f}-{cbox[2]:.3f}, y {cbox[1]:.3f}-{cbox[3]:.3f} in"))
        LABEL_CHECK["Fig1"] = log
    data = [dict(item=k, value=v, source=s) for k, v, s in [
        ("all cells", n_all, "results/senginv/addendum_meta.json n_cells"),
        ("WT group", n_wt, "results/senginv/addendum_meta.json n_match_barcode"),
        ("CRA group", grp["barcode `_A`"], "results/senginv/ADDENDUM_SENGINV.md group table"),
        ("CRI group", grp["barcode `_I`"], "results/senginv/ADDENDUM_SENGINV.md group table"),
        ("PD14", grp["`PD14`"], "results/senginv/ADDENDUM_SENGINV.md group table"),
        ("PD26", grp["`PD26`"], "results/senginv/ADDENDUM_SENGINV.md group table"),
        ("PD32", grp["`PD32`"], "results/senginv/ADDENDUM_SENGINV.md group table"),
        ("NT CRA", nt_cra, "results/seng2/d1_summary.json sizes.NT_CRA"),
        ("NT CRI", nt_cri, "results/seng2/d1_summary.json sizes.NT_CRI"),
        ("perturbations CRA", n_pert_cra, "results/seng3/d4_summary.json summary.delta_i.CRA.n"),
        ("perturbations CRI", n_pert_cri, "results/seng3/d4_summary.json summary.delta_i.CRI.n"),
        ("GTEx young 20-39", n_young, "results/axang/stage0_summary.json anchors n_young"),
        ("GTEx old 60-79", n_old, "results/axang/stage0_summary.json anchors n_old")]]
    save(fig, "Fig1", data)


# --------------------------------------------------------------------------- Fig2
def fig2():
    ff.use_style()
    fig = ff.new_figure(2.9, 1)
    axA = ff.axes_in(fig, 0.62, 0.55, 1.95, 1.9)
    axB = ff.axes_in(fig, 3.05, 0.55, 1.85, 1.9)
    axC = ff.axes_in(fig, 5.75, 0.55, 1.55, 1.65)
    data = []

    # A: effect-free null delta
    src = SENG3 / "d1_null_validation.csv"
    nv = csv_rows(src)
    off = {}
    k = 0
    for sp in ("S2A_k20", "S1"):
        for lib in LIBS:
            for org in ("i", "ii"):
                off[(sp, lib, org)] = -0.35 + k * 0.1
                k += 1
    for _, r in nv.iterrows():
        x = DEPTHS.index(int(r.n_cells)) + off[(r.space, r.library, r.origin)]
        c = SPACE_COLOR[r.space]
        axA.plot([x, x], [r.eval_lo, r.eval_hi], color=c, lw=ff.LW_AXIS, zorder=2)
        dot(axA, x, r.eval_median, c, LIB_MARKER[r.library], r.origin == "i")
        data.append(dict(panel="A", space=r.space, library=r.library, origin=r.origin, n_cells=r.n_cells,
                         null_delta_median=r.eval_median, lo_2_5=r.eval_lo, hi_97_5=r.eval_hi,
                         firing_rate=r.firing_rate, n_eval_negative_of_200=r.n_eval_negative,
                         source=f"{rel(src)}:{r._line}"))
    axA.axhline(0, color=ff.BLACK, lw=ff.LW_AXIS, zorder=0)
    axA.set_xticks(range(4), [str(d) for d in DEPTHS])
    axA.set_xlim(-0.6, 3.6)
    axA.set_xlabel("Cells per draw")
    axA.set_ylabel("Change in distance to PD14,\neffect-free draws (median, 95% range)")
    axA.set_title("Effect-free null", loc="left")
    axA.legend(handles=legend_handles([("S2A_k20", SPACE_COLOR["S2A_k20"], "o", True),
                                       ("Full-gene", SPACE_COLOR["S1"], "o", True),
                                       ("CRA", ff.DARKGREY, "o", True), ("CRI", ff.DARKGREY, "s", True),
                                       ("origin (i)", ff.DARKGREY, "o", True),
                                       ("origin (ii)", ff.DARKGREY, "o", False)]),
               loc="upper right", frameon=False, ncol=2, handletextpad=0.2, columnspacing=0.6,
               borderaxespad=0.1)
    axA.set_ylim(-14, 72)

    # B: planted-approach power, both spaces from seng3 (author's answer Q1)
    src = SENG3 / "d2_plant.csv"
    pl = csv_rows(src)
    doses = (0.25, 0.5, 1.0)
    gx = {("S2A_k20", "i"): -0.27, ("S2A_k20", "ii"): -0.09, ("S1", "i"): 0.09, ("S1", "ii"): 0.27}
    lib_dx = dict(CRA=-0.035, CRI=0.035)
    for _, r in pl.iterrows():
        data.append(dict(panel="B", space=r.space, library=r.library, origin=r.origin, n_cells=r.n_cells,
                         dose=r.dose, power=r.power, drawn=r.dose in doses, source=f"{rel(src)}:{r._line}"))
        if r.dose not in doses:
            continue
        x = doses.index(r.dose) + gx[(r.space, r.origin)] + lib_dx[r.library]
        dot(axB, x, r.power, SPACE_COLOR[r.space], LIB_MARKER[r.library], r.origin == "i", ms=ff.MS * 0.85)
    axB.axhline(0.8, color=ff.GREY, lw=ff.LW_AXIS, ls=(0, (4, 2.5)), zorder=0)
    axB.text(-0.48, 0.81, "0.80", color=ff.DARKGREY, va="bottom", ha="left")
    axB.set_xticks(range(3), ["25%", "50%", "100%"])
    axB.set_xlim(-0.5, 2.5)
    axB.set_ylim(-0.03, 1.05)
    axB.set_xlabel("Planted approach (share of the aging gap)")
    axB.set_ylabel("Power (share of 200 draws called)")
    axB.set_title("Planted approach, 96–157 cells", loc="left")

    # C: the four factors, called or not
    cols = [("Full-gene\n(i)", SENG2 / "d2_table.csv", "called_delta_i", SENG2 / "d2_summary.json", "delta_i"),
            ("Full-gene\n(ii)", SENG2 / "d2_table.csv", "called_delta_ii", SENG2 / "d2_summary.json", "delta_ii"),
            ("S2A_k20\n(i)", SENG3 / "d4_table.csv", "called_delta_i", SENG3 / "d4_summary.json", "delta_i"),
            ("S2A_k20\n(ii)", SENG3 / "d4_table.csv", "called_delta_ii", SENG3 / "d4_summary.json", "delta_ii")]
    for j, (lab, tab, col, summ, key) in enumerate(cols):
        sp = "S1" if j < 2 else "S2A_k20"
        for i, f in enumerate(FOUR):
            v, line = csv_at(tab, col, perturbation=f)
            dot(axC, j, i, SPACE_COLOR[sp], "s", bool(v), ms=ff.MS * 1.4)
            data.append(dict(panel="C", column=lab.replace("\n", " "), factor=f, called=bool(v),
                             source=f"{rel(tab)}:{line} {col}"))
        n_called, line = json_at(summ, "summary", key, "called")
        n_tot, _ = json_at(summ, "summary", key, "n")
        chance, _ = json_at(summ, "summary", key, "chance_expected")
        axC.text(j, 3.95, f"{n_called}", ha="center", va="center")
        data.append(dict(panel="C", column=lab.replace("\n", " "), factor="all", n_called=n_called, n=n_tot,
                         chance_expected=chance, source=f"{rel(summ)}:{line}"))
    axC.text(-0.75, 3.95, f"called of {n_tot}", ha="right", va="center")
    axC.text(-0.75, 4.55, f"chance {chance:.1f}", ha="right", va="center", color=ff.DARKGREY)
    for j0, sp in ((0.5, "S1"), (2.5, "S2A_k20")):
        axC.text(j0, -1.05, SPACE_WORDS[sp], ha="center", va="bottom", color=SPACE_COLOR[sp])
    axC.set_xlim(-0.6, 3.6)
    axC.set_ylim(4.8, -0.5)
    axC.set_yticks(range(4), FOUR)
    axC.set_xticks(range(4), ["(i)", "(ii)", "(i)", "(ii)"])
    axC.xaxis.tick_top()
    axC.tick_params(axis="x", length=0, pad=2)
    axC.tick_params(axis="y", length=0)
    for s in axC.spines.values():
        s.set_visible(False)
    axC.set_xlabel("filled = called (below floor and closer)", color=ff.DARKGREY)
    ff.panel_letters(fig, [[(axA, "A"), (axB, "B"), (axC, "C")]])
    save(fig, "Fig2", data)


# --------------------------------------------------------------------------- angle helpers (Fig3A, S3)
def angles_s2a(aging, prolif):
    """Angle from parallel in S2A_k20, (value, source). C1_lu from the corrected CSV only."""
    out = {}
    corr = AXANG / "d1_angles_corrected_C1lu.csv"
    raw = AXANG / "d1_angles.csv"
    for a in aging:
        for p in prolif:
            if p == "C1_lu":
                v, line = csv_at(corr, "angle_from_parallel_corrected", space="S2A_k20", axis_a=a, axis_b=p)
                out[(a, p)] = (float(v), f"{rel(corr)}:{line} angle_from_parallel_corrected", None, None)
            else:
                th, line = csv_at(raw, "angle_deg", space="S2A_k20", axis_a=a, axis_b=p)
                c, _ = csv_at(raw, "cos", space="S2A_k20", axis_a=a, axis_b=p)
                out[(a, p)] = (from_parallel(th), f"{rel(raw)}:{line} min(angle_deg, 180-angle_deg)",
                               float(th), float(c))
    return out


def draw_angles(ax, vals, prolif, aging=("A1", "A2", "A3", "A4")):
    dy = np.linspace(-0.27, 0.27, len(prolif))
    for i, a in enumerate(aging):
        for k, p in enumerate(prolif):
            color, marker = PROLIF_STYLE[p]
            dot(ax, vals[(a, p)][0], i + dy[k], color, marker, True)
    for x in (45,):
        ax.axvline(x, color=ff.GREY, lw=ff.LW_AXIS, ls=(0, (4, 2.5)), zorder=0)
    ax.set_xlim(0, 90)
    ax.set_xticks([0, 15, 30, 45, 60, 75, 90])
    ax.set_ylim(len(aging) - 0.5, -0.5)
    ax.set_yticks(range(len(aging)), [AGING_WORDS[a] for a in aging])
    ax.set_xlabel("Angle from parallel, S2A_k20 (degrees)\n0 = same direction, 90 = orthogonal")


def angle_rows(fig, panel, vals):
    return [dict(panel=panel, aging_axis=a, prolif_axis=p, angle_from_parallel=v[0], signed_angle=v[2],
                 cos=v[3], source=v[1]) for (a, p), v in vals.items()]


# --------------------------------------------------------------------------- Fig3
def fig3():
    ff.use_style()
    fig = ff.new_figure(4.5, 2)
    axA = ff.axes_in(fig, 1.2, 2.75, 2.9, 1.35)
    axB = ff.axes_in(fig, 0.55, 0.5, 2.6, 1.6)
    axC = ff.axes_in(fig, 4.45, 0.5, 2.9, 1.6)
    data = []

    prolif = ("C1_south", "C1_seng_CRA", "C1_seng_CRI", "C1_lu")
    vals = angles_s2a(("A1", "A2", "A3", "A4"), prolif)
    draw_angles(axA, vals, prolif)
    axA.set_title("Aging axes vs cell-derived cell-division axes", loc="left")
    axA.legend(handles=legend_handles([(PROLIF_WORDS[p], *PROLIF_STYLE[p], True) for p in prolif]),
               loc="upper left", bbox_to_anchor=(1.03, 1.0), frameon=False, handletextpad=0.2)
    data += angle_rows(fig, "A", vals)

    # B: power at a 50% plant
    src = SENG4 / "d2_plant_residualized.csv"
    rp = csv_rows(src, dose=0.5)
    groups = ["SENG3", "NONE", "P1", "P2", "P3"]
    gwords = ["SENG3\nfull NT", "none", "P1", "P2", "P3"]
    lib_dx = dict(CRA=-0.12, CRI=0.12)
    org_dx = dict(i=-0.05, ii=0.05)
    for _, r in rp.iterrows():
        x = groups.index(r.axis) + lib_dx[r.library] + org_dx[r.origin]
        dot(axB, x, r.power, ff.BLACK if r.axis != "NONE" else ff.DARKGREY, LIB_MARKER[r.library],
            r.origin == "i", ms=ff.MS * 0.85)
        data.append(dict(panel="B", axis_removed=r.axis, library=r.library, n_cells=r.n_cells, origin=r.origin,
                         power_50=r.power, source=f"{rel(src)}:{r._line} power"))
        if r.axis == "NONE":
            x0 = lib_dx[r.library] + org_dx[r.origin]
            dot(axB, x0, r.seng3_power_unresidualized, ff.VERMILLION, LIB_MARKER[r.library], r.origin == "i",
                ms=ff.MS * 0.85)
            data.append(dict(panel="B", axis_removed="SENG3 (unresidualized, full NT origin)",
                             library=r.library, n_cells=r.n_cells, origin=r.origin,
                             power_50=r.seng3_power_unresidualized,
                             source=f"{rel(src)}:{r._line} seng3_power_unresidualized"))
    axB.axhline(0.8, color=ff.GREY, lw=ff.LW_AXIS, ls=(0, (4, 2.5)), zorder=0)
    axB.set_xticks(range(5), gwords)
    axB.set_xlim(-0.5, 4.5)
    axB.set_ylim(-0.03, 1.05)
    axB.set_xlabel("Cell-division axis projected out")
    axB.set_ylabel("Power at a 50% planted approach")
    axB.set_title("Residualizing removes the signal", loc="left")
    axB.legend(handles=legend_handles([("CRA", ff.DARKGREY, "o", True), ("CRI", ff.DARKGREY, "s", True),
                                       ("origin (i)", ff.DARKGREY, "o", True),
                                       ("origin (ii)", ff.DARKGREY, "o", False)]),
               loc="center right", frameon=False, handletextpad=0.2, borderaxespad=0.1)

    # C: proliferation shift vs delta (i)
    src = SENG3 / "d4_table.csv"
    t = csv_rows(src)
    for called in (False, True):
        s = t[t.called_delta_i == called]
        axC.scatter(s.prolif_shift, s.delta_i, s=(ff.MS * 0.75) ** 2, c=CALLED if called else NOT_CALLED,
                    lw=0, zorder=3 if called else 2)
    for f, (dx, dy) in zip(FOUR, [(10, -2), (8, 10), (-8, 9), (-10, -4)]):
        r = t[t.perturbation == f].iloc[0]
        if not FIX:
            axC.annotate(f, (r.prolif_shift, r.delta_i), xytext=(dx, dy), textcoords="offset points",
                         ha="left" if dx > 0 else "right", va="center",
                         arrowprops=dict(arrowstyle="-", lw=ff.LW_AXIS * 0.8, color=ff.BLACK, shrinkA=0,
                                         shrinkB=1.5))
        dot(axC, r.prolif_shift, r.delta_i, ff.BLACK, "o", False, zorder=4)
    for _, r in t.iterrows():
        data.append(dict(panel="C", perturbation=r.perturbation, prolif_shift=r.prolif_shift, delta_i=r.delta_i,
                         called_delta_i=bool(r.called_delta_i), source=f"{rel(src)}:{r._line}"))
    rho, line = json_at(SENG3 / "d5_confounds.json", "spearman", "prolif_shift_vs_delta_i", "rho")
    data.append(dict(panel="C", item="stored Spearman rho", value=rho,
                     source=f"{rel(SENG3 / 'd5_confounds.json')}:{line}"))
    n_called, _ = json_at(SENG3 / "d4_summary.json", "summary", "delta_i", "called")
    rho_text = axC.text(0.97, 0.97, f"Spearman ρ = {rho:.4f}\n(stored; n = {len(t)})", transform=axC.transAxes,
             ha="right", va="top")
    axC.axhline(0, color=ff.BLACK, lw=ff.LW_AXIS, zorder=0)
    axC.axvline(0, color=ff.BLACK, lw=ff.LW_AXIS, zorder=0)
    axC.set_xlabel("Cell-division shift (perturbation − NT)")
    axC.set_ylabel("Change in distance to PD14, origin (i)")
    axC.set_title(f"All {len(t)} perturbations, S2A_k20", loc="left")
    axC.legend(handles=legend_handles([(f"called ({n_called})", CALLED, "o", True),
                                       ("not called", NOT_CALLED, "o", True)]),
               loc="lower left", frameon=False, handletextpad=0.2, borderaxespad=0.1)
    if FIX:
        four = [(f, tuple(t.loc[t.perturbation == f, ["prolif_shift", "delta_i"]].iloc[0])) for f in FOUR]
        place_labels(fig, axC, "Fig3C", four, t[["prolif_shift", "delta_i"]].to_numpy(),
                     avoid=[axC.get_legend(), rho_text],
                     fallback=dict(zip(FOUR, [(10, -2), (8, 10), (-8, 9), (-10, -4)])) if LABEL_FALLBACK else None, manual=FIG3_MANUAL)
    ff.panel_letters(fig, [[(axA, "A")], [(axB, "B"), (axC, "C")]])
    save(fig, "Fig3", data)


# --------------------------------------------------------------------------- Fig4
def fig4():
    ff.use_style()
    s3 = csv_rows(SENG3 / "d4_table.csv")
    ax3 = csv_rows(AXANG / "d3b_table.csv")
    j = s3.merge(ax3, on="perturbation", how="inner", suffixes=("_s3", "_ax"))
    if len(j) != 340 or len(s3) != 340 or len(ax3) != 340 \
            or (j.modality_s3 != j.modality_ax).any() or (j.n_cells_s3 != j.n_cells_ax).any():
        stop("join_mismatch", f"join gives {len(j)} rows; modality/n_cells agree: "
             f"{(j.modality_s3 == j.modality_ax).all()}/{(j.n_cells_s3 == j.n_cells_ax).all()}")
    cats = [("both tests", True, True, ff.BLACK, "o", True),
            ("early-passage test only", True, False, CALLED, "o", True),
            ("GTEx test only", False, True, ff.BLUE, "s", True),
            ("neither", False, False, NOT_CALLED, "o", True)]
    fig = ff.new_figure(3.0, 1)
    ax = ff.axes_in(fig, 0.6, 0.5, 3.4, 2.35)
    data = []
    counts = {}
    for lab, c_s3, c_ax, color, marker, filled in cats[::-1]:
        s = j[(j.called_delta_i == c_s3) & (j.called == c_ax)]
        counts[lab] = len(s)
        ax.scatter(s.delta_i, s.delta, s=(ff.MS * 0.75) ** 2, c=color, marker=marker, lw=0,
                   zorder=2 if lab == "neither" else 3)
    for f, (dx, dy) in zip(FOUR, [(10, -7), (10, 6), (10, -6), (10, 4)]):
        r = j[j.perturbation == f].iloc[0]
        dot(ax, r.delta_i, r.delta, ff.BLACK, "o", False, zorder=4)
        if not FIX:
            ax.annotate(f, (r.delta_i, r.delta), xytext=(dx, dy), textcoords="offset points",
                        ha="left" if dx > 0 else "right", va="center",
                        arrowprops=dict(arrowstyle="-", lw=ff.LW_AXIS * 0.8, color=ff.BLACK, shrinkA=0,
                                        shrinkB=1.5))
    ax.axhline(0, color=ff.BLACK, lw=ff.LW_AXIS, zorder=0)
    ax.axvline(0, color=ff.BLACK, lw=ff.LW_AXIS, zorder=0)
    ax.set_xlabel("Change in distance to early-passage cells (PD14)\nS2A_k20, origin (i); negative = closer")
    ax.set_ylabel("Change in distance to GTEx young donors\nS2A_k20; negative = closer")
    ax.legend(handles=legend_handles([(f"called on {lab} ({counts[lab]})", col, m, fl)
                                      for lab, _, _, col, m, fl in cats]),
              loc="upper left", bbox_to_anchor=(1.03, 1.0), frameon=False, handletextpad=0.2,
              title="Marker: called (below floor and closer)", title_fontsize=ff.FS, alignment="left")
    s3_called, l1 = json_at(SENG3 / "d4_summary.json", "summary", "delta_i", "called")
    ax_called, l2 = json_at(AXANG / "d3_summary.json", "summary", "called")
    ax_closer, l3 = json_at(AXANG / "d3_summary.json", "summary", "decreased")
    chance, _ = json_at(AXANG / "d3_summary.json", "summary", "chance_expected")
    fig.text(4.43 / 7.5, 0.62 / 3.0,
             f"Stored counts of {len(j)}, chance expectation {chance:.1f} each:\n"
             f"early-passage test: {s3_called} called\n"
             f"GTEx young test: {ax_closer} closer, {ax_called} called\n"
             "The GTEx test transplants a single-cell displacement\nonto a bulk GTEx origin.",
             ha="left", va="bottom")
    if FIX:
        four = [(f, tuple(j.loc[j.perturbation == f, ["delta_i", "delta"]].iloc[0])) for f in FOUR]
        if FIG4_MANUAL is None:
            place_labels(fig, ax, "Fig4", four, j[["delta_i", "delta"]].to_numpy())
        else:
            annots = [ax.annotate(f, xy, xytext=FIG4_MANUAL[f], textcoords="data", ha="left", va="center",
                                  zorder=5, arrowprops=dict(arrowstyle="-", lw=ff.LW_AXIS * 0.8, color=ff.BLACK,
                                                            shrinkA=0, shrinkB=1.5))
                      for f, xy in four]
            check_labels(fig, ax, "Fig4", annots, j[["delta_i", "delta"]].to_numpy(),
                         how=[f"by hand, text at {FIG4_MANUAL[f]} (data)" for f, _ in four])
    for _, r in j.iterrows():
        data.append(dict(perturbation=r.perturbation, modality=r.modality_s3, n_cells=r.n_cells_s3,
                         delta_early_passage_i=r.delta_i, called_early_passage_i=bool(r.called_delta_i),
                         delta_gtex_young=r.delta, called_gtex_young=bool(r.called),
                         source=f"{rel(SENG3 / 'd4_table.csv')}:{r._line_s3}; {rel(AXANG / 'd3b_table.csv')}:{r._line_ax}"))
    data += [dict(perturbation="all", item=k, value=v, source=s) for k, v, s in [
        ("seng3 called (i)", s3_called, f"{rel(SENG3 / 'd4_summary.json')}:{l1}"),
        ("axang called", ax_called, f"{rel(AXANG / 'd3_summary.json')}:{l2}"),
        ("axang decreased", ax_closer, f"{rel(AXANG / 'd3_summary.json')}:{l3}")]]
    save(fig, "Fig4", data)


# --------------------------------------------------------------------------- S1
def figS1():
    ff.use_style()
    src = SENG3 / "d4_table.csv"
    t = csv_rows(src)
    fig = ff.new_figure(2.9, 1)
    axs = [ff.axes_in(fig, 0.6, 0.5, 2.3, 1.85), ff.axes_in(fig, 3.45, 0.5, 2.3, 1.85)]
    cats = [("younger, not called", True, False, ff.BLUE), ("younger and called", True, True, ff.BLACK),
            ("called, not younger", False, True, CALLED), ("neither", False, False, NOT_CALLED)]
    data = []
    for ax, org in zip(axs, ("i", "ii")):
        summ = SENG3 / "d4_summary.json"
        n_y, l1 = json_at(summ, "ruler", org, "n_younger")
        n_ync, l2 = json_at(summ, "ruler", org, "n_younger_not_called")
        for lab, yo, ca, color in cats[::-1]:
            s = t[(t.ruler_scores_younger == yo) & (t[f"called_delta_{org}"] == ca)]
            ax.scatter(s[f"delta_{org}"], s.ruler_change, s=(ff.MS * 0.75) ** 2, c=color, lw=0,
                       zorder=2 if lab == "neither" else 3)
        for f in FOUR:
            r = t[t.perturbation == f].iloc[0]
            dot(ax, r[f"delta_{org}"], r.ruler_change, ff.BLACK, "o", False, zorder=4)
        ax.axhline(0, color=ff.BLACK, lw=ff.LW_AXIS, zorder=0)
        ax.axvline(0, color=ff.BLACK, lw=ff.LW_AXIS, zorder=0)
        ax.set_xlabel(f"Change in distance to PD14, origin ({org})")
        ax.set_ylabel("Frozen-ruler score change\n(negative = scores younger)")
        ax.set_title(f"Origin ({org}): {n_y} of {len(t)} score younger,\n{n_ync} of them not called",
                     loc="left")
        data.append(dict(origin=org, item="n_younger", value=n_y, source=f"{rel(summ)}:{l1}"))
        data.append(dict(origin=org, item="n_younger_not_called", value=n_ync, source=f"{rel(summ)}:{l2}"))
        for _, r in t.iterrows():
            data.append(dict(origin=org, perturbation=r.perturbation, delta=r[f"delta_{org}"],
                             ruler_change=r.ruler_change, ruler_scores_younger=bool(r.ruler_scores_younger),
                             called=bool(r[f"called_delta_{org}"]), source=f"{rel(src)}:{r._line}"))
    axs[1].legend(handles=legend_handles([(lab, c, "o", True) for lab, _, _, c in cats]
                                         + [("the four factors", ff.BLACK, "o", False)]),
                  loc="upper left", bbox_to_anchor=(1.03, 1.0), frameon=False, handletextpad=0.2)
    ff.panel_letters(fig, [[(axs[0], "A"), (axs[1], "B")]])
    save(fig, "FigS1", data)


# --------------------------------------------------------------------------- S2
def figS2():
    ff.use_style()
    fig = ff.new_figure(4.5, 2)
    axA = ff.axes_in(fig, 1.05, 2.65, 2.3, 1.45)
    axB = ff.axes_in(fig, 4.75, 2.65, 2.3, 1.45)
    axC = ff.axes_in(fig, 1.05, 0.5, 2.3, 1.45)
    axD = ff.axes_in(fig, 4.75, 0.5, 2.3, 1.45)
    data = []
    src = SENG3 / "d3_offset.csv"
    rows = [("S1", "A"), ("S1", "I"), ("S2A_k20", "A"), ("S2A_k20", "I")]
    words = [f"{SPACE_WORDS[s]}, {'CRA' if l == 'A' else 'CRI'}" for s, l in rows]
    for ax, col, xl in ((axA, "thin_offset_over_a", "Library offset / ||a||"),
                        (axB, "thin_cos_offset_a", "cos(library offset, a)")):
        for i, (sp, lib) in enumerate(rows):
            r = csv_rows(src, variant="half_A", space=sp, library=lib).iloc[0]
            ax.plot([r[col + "_lo"], r[col + "_hi"]], [i, i], color=SPACE_COLOR[sp], lw=ff.LW_DATA)
            dot(ax, r[col], i, SPACE_COLOR[sp], LIB_MARKER["CRA" if lib == "A" else "CRI"], True)
            data.append(dict(panel="A" if ax is axA else "B", space=sp, library=lib, quantity=col,
                             thinned_median=r[col], lo_2_5=r[col + "_lo"], hi_97_5=r[col + "_hi"],
                             source=f"{rel(src)}:{r._line}"))
        ax.set_yticks(range(4), words)
        ax.set_ylim(3.6, -0.6)
        ax.set_xlabel(xl)
        ax.axvline(0 if ax is axB else 1, color=ff.GREY, lw=ff.LW_AXIS, ls=(0, (4, 2.5)), zorder=0)
    axA.set_title("Library offset, PD14 half A, depth-matched", loc="left")
    axB.set_title("Its direction against the aging axis a", loc="left")

    src = SENG3 / "d4_table.csv"
    t = csv_rows(src)
    bins = np.arange(0, t.n_cells.max() + 20, 10)
    for lib, color in (("CRA", ff.VERMILLION), ("CRI", ff.BLUE)):
        v = t[t.modality == lib].n_cells
        axC.hist(v, bins=bins, histtype="step", color=color, lw=ff.LW_DATA, label=f"{lib} ({len(v)})")
    for f in FOUR:
        n = int(t[t.perturbation == f].n_cells.iloc[0])
        axC.axvline(n, color=ff.BLACK, lw=ff.LW_AXIS, ls=(0, (2, 2)), zorder=0)
    for _, r in t.iterrows():
        data.append(dict(panel="C", perturbation=r.perturbation, modality=r.modality, n_cells=r.n_cells,
                         source=f"{rel(src)}:{r._line}"))
    axC.set_xlabel("Cells per perturbation (dotted: the four factors)")
    axC.set_ylabel("Perturbations")
    axC.legend(frameon=False, loc="upper right")
    axC.set_title("Cells per perturbation, by library", loc="left")

    d2 = AXANG / "d2_summary.json"
    cv = AXANG / "correction_values.json"
    cov = []
    for name in ("Southard", "Sengstack"):
        pr, l1 = json_at(d2, name, "present")
        ab, _ = json_at(d2, name, "absent")
        pc, l3 = json_at(d2, name, "pct_gap_sq_in_unmeasured")
        cov.append((name, pr, ab, pc, f"{rel(d2)}:{l1}-{l3}"))
    pr, l1 = json_at(cv, "lu_present")
    ab, _ = json_at(cv, "lu_absent")
    pc, l3 = json_at(cv, "pct_gap_sq_unmeasured_corrected")
    cov.append(("Lu (corrected)", pr, ab, pc, f"{rel(cv)}:{l1}-{l3}"))
    for i, (name, pr, ab, pc, s) in enumerate(cov):
        axD.barh(i, pr, color=ff.GREY, height=0.6)
        axD.barh(i, ab, left=pr, color="white", edgecolor=ff.BLACK, lw=ff.LW_AXIS, height=0.6)
        axD.text(300, i, f"{pr:,} present\n{pc:.1f}% of GTEx age gap unmeasured", va="center",
                 ha="left")
        data.append(dict(panel="D", dataset=name, present=pr, absent=ab, pct_gap_sq_unmeasured=pc, source=s))
    axD.set_yticks(range(3), [c[0] for c in cov])
    axD.set_ylim(2.6, -0.6)
    axD.set_xlim(0, 23485)
    axD.set_xticks([0, 10000, 20000, 23485], ["0", "10,000", "20,000", "23,485"])
    axD.set_xlabel("Ruler genes (grey present, white absent)")
    axD.set_title("Ruler-gene coverage", loc="left")
    axD.set_clip_on(False)
    ff.panel_letters(fig, [[(axA, "A"), (axB, "B")], [(axC, "C"), (axD, "D")]])
    save(fig, "FigS2", data)


# --------------------------------------------------------------------------- S3
def figS3():
    ff.use_style()
    fig = ff.new_figure(2.4, 1)
    ax = ff.axes_in(fig, 1.3, 0.65, 2.6, 1.45)
    prolif = ("C2", "C1_gtex")
    vals = angles_s2a(("A1", "A2", "A3", "A4"), prolif)
    draw_angles(ax, vals, prolif)
    ax.legend(handles=legend_handles([(PROLIF_WORDS[p], *PROLIF_STYLE[p], True) for p in prolif]),
              loc="upper left", bbox_to_anchor=(1.03, 1.0), frameon=False, handletextpad=0.2)
    ax.set_title("Aging axes vs gene-list and bulk-donor cell-division axes", loc="left")
    save(fig, "FigS3", angle_rows(fig, "S3", vals))


def save(fig, stem, data):
    ff.OUT = OUT  # figures_final.save writes to its module-level OUT
    ff.save(fig, stem, data)


# --------------------------------------------------------------------------- angle check (author's Q2)
def angle_check():
    fa = ROOT / "FINDINGS_AXANG.md"
    lines = fa.read_text(encoding="utf-8").splitlines()
    start, _ = md_line(fa, "### Each aging axis against every proliferation construction, S2A_k20")
    end, _ = md_line(fa, "### Aging axis against aging axis, S2A_k20", start)
    a_words = {"A1 GTEx donor age": "A1", "A2 frozen ruler": "A2", "A3 Sengstack passage": "A3",
               "A4 Lu donor pair": "A4"}
    p_words = {"C1 Southard NT cells (14,013)": "C1_south", "C1 Sengstack NT CRA (333)": "C1_seng_CRA",
               "C1 Sengstack NT CRI (739)": "C1_seng_CRI", "C1 Lu day-0 (12,803)": "C1_lu",
               "C2 Tirosh gene list (94 genes)": "C2", "C1 GTEx donors (652, bulk)": "C1_gtex"}
    table = {}
    for i in range(start, end):
        cells = [c.strip() for c in lines[i].split("|")[1:-1]]
        if len(cells) == 7 and cells[0] in a_words and cells[1] in p_words:
            table[(a_words[cells[0]], p_words[cells[1]])] = (float(cells[5].strip("*°")), i + 1)
    co = AXANG / "CORRECTION.md"
    clines = co.read_text(encoding="utf-8").splitlines()
    cstart, _ = md_line(co, "**S2A_k20** (the primary space):")
    corr = {}
    cend, _ = md_line(co, "**Full-gene space**", cstart)  # the S2A_k20 table only, not the full-gene one
    for i in range(cstart, cend - 1):
        m = re.match(r"\|\s*\**(A\d)[^|]*vs C1_lu\**\s*\|[^|]*\|[^|]*\|\s*\**([\d.]+)°", clines[i])
        if m:
            corr[m.group(1)] = (float(m.group(2)), i + 1)
    if sorted(corr) != ["A1", "A2", "A3", "A4"]:
        stop("value_not_stored", f"CORRECTION.md S2A_k20 table: found C1_lu rows for {sorted(corr)}")
    rows = []
    plotted = [("Fig3A", p) for p in ("C1_south", "C1_seng_CRA", "C1_seng_CRI", "C1_lu")] + \
              [("FigS3", p) for p in ("C2", "C1_gtex")]
    for fig, p in plotted:
        for a in ("A1", "A2", "A3", "A4"):
            v, src, _, _ = angles_s2a((a,), (p,))[(a, p)]
            if (a, p) not in table:
                stop("value_not_stored", f"FINDINGS_AXANG.md §5 has no row for {a} vs {p}")
            fv, fl = table[(a, p)]
            row = dict(figure=fig, pair=f"{a} vs {p}", plotted=v, plotted_2dp=f"{v:.2f}", source=src,
                       findings_s5=f"{fv:.2f}", findings_line=f"FINDINGS_AXANG.md:{fl}",
                       match_s5=round(v, 2) == fv, correction="", correction_line="", match_correction="")
            if p == "C1_lu":
                cv, cl = corr[a]
                row.update(correction=f"{cv:.2f}", correction_line=f"results/axang/CORRECTION.md:{cl}",
                           match_correction=round(v, 2) == cv)
            rows.append(row)
    return rows


# --------------------------------------------------------------------------- paper2_numbers.md
def numbers():
    R = []

    def add(label, stated, stored, src, line, d, note="", kind="value"):
        """MATCH if `stated` equals `stored` rounded to `d` decimals (d=None: exact)."""
        if isinstance(stored, bool) or isinstance(stated, str):
            ok = stored == stated
        elif d is None:
            ok = stored == stated
        else:
            ok = round(float(stored), d) == float(stated)
        R.append(dict(label=label, stated=stated, stored=stored, where=f"{rel(src)}:{line}",
                      result="MATCH" if ok else "MISMATCH", note=note))

    m = SENGINV / "addendum_meta.json"
    v, l = json_at(m, "n_cells"); add("cells in the h5ad", 36514, v, m, l, None)
    a = SENGINV / "ADDENDUM_SENGINV.md"
    for g, s in (("`PD14`", 581), ("`PD26`", 1880), ("`PD32`", 4348), ("barcode `_A`", 14333),
                 ("barcode `_I`", 15372)):
        l, txt = md_line(a, f"| {g} |")
        add(f"{g.strip('`').replace('barcode `', 'barcode ')} cells", s, int(txt.split("|")[2]), a, l, None,
            "group table of the SENGINV addendum (a markdown results file)")
    m = SENG2 / "d1_summary.json"
    for k, s in (("NT_CRA", 666), ("NT_CRI", 1479)):
        v, l = json_at(m, "sizes", k); add(f"{k} cells", s, v, m, l, None)
    m = AXANG / "d2_summary.json"
    v, l = json_at(m, "Sengstack", "present"); add("Sengstack ruler genes present", 16781, v, m, l, None)
    v2, l2 = json_at(m, "Sengstack", "absent")
    add("ruler genes in total (present + absent)", 23485, v + v2, m, f"{l}-{l2}", None,
        "sum of two stored counts; also results/axang/d1_summary.json meta.n_ruler")
    o = SENG3 / "d3_offset.csv"
    for lib, s_off, s_cos in (("A", 1.0532, 0.1123), ("I", 0.7298, 0.0708)):
        v, l = csv_at(o, "thin_offset_over_a", variant="half_A", space="S1", library=lib)
        add(f"offset/||a||, half A, full-gene, {lib}", s_off, v, o, l, 4, "thinned median")
        v, l = csv_at(o, "thin_cos_offset_a", variant="half_A", space="S1", library=lib)
        add(f"cos(offset, a), half A, full-gene, {lib}", s_cos, v, o, l, 4, "thinned median")
    for sp, lib, s, d in (("S1", "A", 0.6445, 4), ("S1", "I", 0.7934, 4),
                          ("S2A_k20", "A", 0.811375, 6), ("S2A_k20", "I", 0.887005, 6)):
        v, l = csv_at(o, "thin_cos_t_a", variant="half_A", space=sp, library=lib)
        add(f"cos(target, aging), {SPACE_WORDS[sp]}, {lib}", s, v, o, l, d, "thinned median")
    m = SENG4 / "d1_summary.json"
    v, l = json_at(m, "cos_a_u_full_gene"); add("cos(a, u), full-gene", -0.040273, v, m, l, 6)
    nv = SENG3 / "d1_null_validation.csv"
    t = csv_rows(nv, space="S1", origin="i")
    lo, hi = t.loc[t.eval_median.idxmin()], t.loc[t.eval_median.idxmax()]
    add("full-gene null median, lowest (origin i)", 20.607, lo.eval_median, nv, lo._line, 3,
        f"{lo.library} {lo.n_cells} cells; lowest of the 8 origin (i) S1 rows")
    add("full-gene null median, highest (origin i)", 33.654, hi.eval_median, nv, hi._line, 3,
        f"{hi.library} {hi.n_cells} cells; highest of the 8 origin (i) S1 rows")
    t = csv_rows(nv, space="S1")
    add("full-gene effect-free draws negative, every S1 row", 0, int(t.n_eval_negative.max()), nv,
        f"{t._line.min()}-{t._line.max()} (S1 rows)", None, "largest n_eval_negative over the 16 S1 rows, of 200")
    t = csv_rows(nv, space="S2A_k20", origin="i")
    lo, hi = t.loc[t.firing_rate.idxmin()], t.loc[t.firing_rate.idxmax()]
    add("compressed firing rate, lowest (i)", 0.040, lo.firing_rate, nv, lo._line, 3, f"{lo.library} {lo.n_cells}")
    add("compressed firing rate, highest (i)", 0.085, hi.firing_rate, nv, hi._line, 3, f"{hi.library} {hi.n_cells}")
    m = SENG3 / "d1_summary.json"
    v, l = json_at(m, "f_max_origin_ii"); add("compressed firing rate, highest (ii)", 0.110, v, m, l, 3)
    p = SENG3 / "d2_plant.csv"
    for dose, s_lo, s_hi in ((0.25, 0.565, 0.740), (0.5, 0.940, 1.000)):
        t = csv_rows(p, space="S2A_k20", origin="i", dose=dose)
        lo, hi = t.loc[t.power.idxmin()], t.loc[t.power.idxmax()]
        add(f"plant power at {int(dose * 100)}%, lowest (S2A_k20, i)", s_lo, lo.power, p, lo._line, 3,
            f"{lo.library} {lo.n_cells}")
        add(f"plant power at {int(dose * 100)}%, highest (S2A_k20, i)", s_hi, hi.power, p, hi._line, 3,
            f"{hi.library} {hi.n_cells}")
    m = SENG3 / "d2_summary.json"
    k20, l1 = json_at(m, "a_norms", "half_A", "S2A_k20")
    s1, l2 = json_at(m, "a_norms", "half_A", "S1")
    fl, _ = md_line(ROOT / "FINDINGS_SENG3.md", "45.1% of the aging gap")
    add("% of the aging gap's length kept by S2A_k20", 45.1, 100 * k20 / s1, m, f"{l1}-{l2}", 1,
        f"ratio of two stored norms ({k20:.4f} / {s1:.4f}); the 45.1% itself is stated only at "
        f"FINDINGS_SENG3.md:{fl}", kind="ratio")
    m = SENG2 / "d2_summary.json"
    v, l = json_at(m, "k_origin_i"); add("full-gene: of the four called, origin (i)", 1, v, m, l, None)
    v, l = json_at(m, "k_origin_ii"); add("full-gene: of the four called, origin (ii)", 0, v, m, l, None)
    m = SENG3 / "d4_summary.json"
    v, l = json_at(m, "summary", "delta_i", "called"); add("compressed called of 340, origin (i)", 80, v, m, l, None)
    v, l = json_at(m, "summary", "delta_ii", "called"); add("compressed called of 340, origin (ii)", 73, v, m, l, None)
    d4 = SENG3 / "d4_table.csv"
    for f, s in zip(FOUR, (-19.211, -13.211, -10.334, -19.936)):
        v, l = csv_at(d4, "delta_i", perturbation=f); add(f"{f} delta to PD14, origin (i)", s, v, d4, l, 3)
    c = SENG3 / "d5_confounds.json"
    for dec, s in ((0, 0.056), (9, 0.545)):
        v, l = json_at(c, "decile_table", dec, "rate_i")
        n_c, _ = json_at(c, "decile_table", dec, "called_i")
        n_n, _ = json_at(c, "decile_table", dec, "n")
        # json_at's line walk skips list indices, so locate the decile's own rate_i line directly
        txt = c.read_text(encoding="utf-8").splitlines()
        l = [i + 1 for i, x in enumerate(txt) if '"rate_i":' in x][dec]
        add(f"call rate, cell-count decile {dec + 1} (origin i)", s, v, c, l, 3,
            f"stored to 4 dp; called/n = {n_c}/{n_n} = {n_c / n_n:.4f}")
    v, l = json_at(c, "depth_band_96_157", "called_i"); add("depth band 96-157: called (i)", 24, v, c, l, None)
    v, l = json_at(c, "depth_band_96_157", "n"); add("depth band 96-157: n", 60, v, c, l, None)
    txt = c.read_text(encoding="utf-8").splitlines()
    plines = [i + 1 for i, x in enumerate(txt) if '"pctl_delta_i_in_band":' in x]
    for k, (f, s) in enumerate(zip(FOUR, (3.3, 15.0, 21.7, 1.7))):
        v = json.loads(c.read_text(encoding="utf-8"))["named_in_band_percentiles"][k]
        if v["factor"] != f:
            stop("value_not_stored", f"named_in_band_percentiles[{k}] is {v['factor']}, not {f}")
        add(f"{f} delta percentile in the depth band", s, v["pctl_delta_i_in_band"], c, plines[k], 1)
    v, l = json_at(c, "spearman", "prolif_shift_vs_delta_i", "rho")
    add("rho(prolif. shift, delta i)", -0.8639, v, c, l, 4)
    va = SENG4 / "d1_vs_aging.csv"
    v, l = csv_at(va, "angle_deg", axis="P1", modality="CRA"); add("angle(P1, a), CRA", 19.66, v, va, l, 2)
    v, l = csv_at(va, "angle_deg", axis="P1", modality="CRI"); add("angle(P1, a), CRI", 20.56, v, va, l, 2)
    v, l = csv_at(va, "pct_a_removed", axis="P1", modality="CRA"); add("% of a removed by P1 (CRA)", 66.4, v, va, l, 1)
    v, l = csv_at(va, "angle_deg", axis="P3", modality="CRA"); add("angle(P3, a)", 61.68, v, va, l, 2)
    v, l = csv_at(va, "pct_a_removed", axis="P3", modality="CRA"); add("% of a removed by P3", 12.0, v, va, l, 1)
    m = SENG4 / "d2_summary.json"
    v, l = json_at(m, "seng3_p50"); add("p50 before residualizing (SENG3)", 0.940, v, m, l, 3)
    v, l = json_at(m, "p50"); add("p50 after projecting out P1", 0.010, v, m, l, 3)
    v, l = json_at(m, "by_axis_origin", "NONE_i", "p50"); add("p50 control (nothing removed, reserved-half origin)", 0.895, v, m, l, 3)
    m = SENG4 / "d3_summary.json"
    v, l = json_at(m, "gate"); add("SENG4 D3 gate", "SPLIT", v, m, l, None)
    v, l = json_at(m, "k_named", "P1", "i"); add("SPLIT: of the four called, P1 origin (i)", 2, v, m, l, None)
    v, l = json_at(m, "k_named", "P1", "ii"); add("SPLIT: of the four called, P1 origin (ii)", 0, v, m, l, None)
    an = AXANG / "d1_angles.csv"
    v, l = csv_at(an, "angle_deg", space="S2A_k20", axis_a="A1", axis_b="A3"); add("angle A1-A3, S2A_k20", 101.50, v, an, l, 2)
    v, l = csv_at(an, "angle_deg", space="S2A_k20", axis_a="A2", axis_b="A3"); add("angle A2-A3, S2A_k20", 100.79, v, an, l, 2)
    m = AXANG / "d1_summary.json"
    v, l = json_at(m, "d_angle_magnitude"); add("A1 vs C1_gtex, angle from parallel", 35.83, v, m, l, 2)
    v, l = json_at(m, "d_cos"); add("cos(A1, C1_gtex)", -0.8108, v, m, l, 4)
    p = AXANG / "d3a_plant.csv"
    t = csv_rows(p, dose=0.25)
    lo, hi = t.loc[t.power.idxmin()], t.loc[t.power.idxmax()]
    add("GTEx plant power at 25%, lowest", 0.945, lo.power, p, lo._line, 3, f"{lo.library} {lo.n_cells}")
    add("GTEx plant power at 25%, highest", 0.990, hi.power, p, hi._line, 3, f"{hi.library} {hi.n_cells}")
    p = AXANG / "d3b_table.csv"
    for f, s in zip(FOUR, (11.448, 9.022, 5.327, 15.336)):
        v, l = csv_at(p, "delta", perturbation=f); add(f"{f} delta to GTEx young", s, v, p, l, 3)
    m = AXANG / "d3_summary.json"
    v, l = json_at(m, "summary", "decreased"); add("GTEx test: closer, of 340", 17, v, m, l, None)
    v, l = json_at(m, "summary", "called"); add("GTEx test: called, of 340", 5, v, m, l, None)
    m = AXANG / "d2_summary.json"
    v, l = json_at(m, "PD32_vs_GTEx", "ratio_to_gap_young"); add("PD32 distance to GTEx young / gap", 14.4, v, m, l, 1)
    v, l = json_at(m, "Sengstack", "pct_gap_sq_in_unmeasured")
    add("% of GTEx age gap in genes Sengstack does not measure", 21.97, v, m, l, 2,
        "Sengstack row. The identical Lu row in this file is superseded (CORRECTION.md C1); "
        "Lu corrected = 8.97% in results/axang/correction_values.json")
    m = SENG3 / "d4_summary.json"
    v, l = json_at(m, "ruler", "i", "n_younger"); add("ruler: score younger, of 340", 181, v, m, l, None)
    v, l = json_at(m, "ruler", "i", "n_younger_not_called"); add("ruler: younger, not called (i)", 116, v, m, l, None)
    v, l = json_at(m, "ruler", "ii", "n_younger_not_called"); add("ruler: younger, not called (ii)", 120, v, m, l, None)
    return R


def write_numbers(R, A):
    def fmt(x):
        return f"{x:.6g}" if isinstance(x, float) else str(x)
    L = ["# PAPER2FIGS: manuscript numbers against the stored results", "",
         "Written by `python src/paper2_figs.py`. Each value is read from the file and line given; MATCH",
         "means the stated value equals the stored value rounded to the decimals as stated (integers",
         "exact). Line numbers are 1-based; for CSVs line 1 is the header. Nothing was adjusted to match.",
         "", f"**{sum(r['result'] == 'MATCH' for r in R)} MATCH, "
             f"{sum(r['result'] == 'MISMATCH' for r in R)} MISMATCH, of {len(R)}.**", "",
         "| # | value | stated | stored | file:line | result | note |", "|---|---|---|---|---|---|---|"]
    for i, r in enumerate(R, 1):
        L.append(f"| {i} | {r['label']} | {fmt(r['stated'])} | {fmt(r['stored'])} | `{r['where']}` | "
                 f"**{r['result']}** | {r['note']} |")
    n_bad = sum(not a["match_s5"] for a in A)
    L += ["", "## Angle check (author's instruction, Q2 of PLAN.flag)", "",
          "Every plotted angle from parallel (Fig 3A, Fig S3) against the FINDINGS_AXANG.md §5 table, to 2 dp.",
          "Plotted values for C1_lu come from the corrected CSV, so they are expected to differ from §5,",
          "which holds the superseded C1_lu values; those rows are also checked against the corrected table",
          "in results/axang/CORRECTION.md.", "",
          f"**{len(A) - n_bad} of {len(A)} match §5; {n_bad} do not.**", "",
          "| figure | pair | plotted | §5 | §5 line | match §5 | CORRECTION.md | line | match corrected |",
          "|---|---|---|---|---|---|---|---|---|"]
    for a in A:
        L.append(f"| {a['figure']} | {a['pair']} | {a['plotted_2dp']} | {a['findings_s5']} | `{a['findings_line']}` | "
                 f"{'yes' if a['match_s5'] else '**NO**'} | {a['correction']} | "
                 f"{'`' + a['correction_line'] + '`' if a['correction_line'] else ''} | "
                 f"{'' if a['match_correction'] == '' else ('yes' if a['match_correction'] else '**NO**')} |")
    L += ["", "Values that do not match §5:", ""]
    L += [f"- {a['figure']} {a['pair']}: plotted {a['plotted_2dp']}°, §5 says {a['findings_s5']}° "
          f"(`{a['findings_line']}`)" + (f"; superseded value, corrected {a['correction']}° at "
                                         f"`{a['correction_line']}`, which the plotted value "
                                         f"{'matches' if a['match_correction'] else 'does NOT match'}"
                                         if a['correction'] else "") for a in A if not a["match_s5"]] or ["- none"]
    (OUT / "paper2_numbers.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    pd.DataFrame(A).to_csv(OUT / "angle_check.csv", index=False)
    print(f"wrote {rel(OUT / 'paper2_numbers.md')} and angle_check.csv")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig1()
    fig2()
    fig3()
    fig4()
    figS1()
    figS2()
    figS3()
    write_numbers(numbers(), angle_check())


if __name__ == "__main__":
    main()

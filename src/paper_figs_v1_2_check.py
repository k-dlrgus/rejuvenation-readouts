"""Check that the v1.2 figure layout changed no plotted value.

Runs the figure functions of src/figures_final.py twice, once as committed at a reference revision
(default 77b74e0, the v12 figures in results/figures_final/) and once from the working tree, with
`save` replaced by a recorder. For every figure it records what is plotted, read back from the
matplotlib artists in data coordinates:

  data    every Line2D x/y, every collection offset and segment (scatter points, error bars),
          every patch geometry, axis limits, tick positions
  labels  every text string: tick labels, axis labels, titles, legend entries, annotations,
          panel letters; plus the position of any text placed in data coordinates

and the plotting-data table each figure function hands to `save`. Data and labels must be
identical (floats compared exactly). Positions of text placed in axes or figure coordinates,
font sizes, line widths and marker sizes are layout, and are listed but not gated.

Also checks each PNG in results/paper_figs_v1_2/: width equal to the text width, height within
the one-row or two-row limit, 300 dpi.

Writes results/paper_figs_v1_2/check_data_identity.csv, check_png_size.csv and check_report.txt.

Usage: python src/paper_figs_v1_2_check.py [REF_REVISION]
"""
from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "paper_figs_v1_2"
FIGS = ("Fig1", "Fig2", "Fig3", "Fig4", "Fig5", "FigS1")

# Text width of paper/rejuvenation_readouts_preprint_v1.2_draft.docx: US Letter 12240 twips wide,
# left and right margins 720 twips, so 10800 twips = 7.5 in. Height limits from the task.
TEXT_WIDTH_IN = 7.5
ROWS = dict(Fig1=1, Fig2=1, Fig3=2, Fig4=1, Fig5=1, FigS1=2)
MAX_H_IN = {1: 3.0, 2: 4.5}


def load_module(path, name):
    sys.path.insert(0, str(ROOT / "src"))
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _arr(a):
    return json.dumps(np.asarray(np.ma.filled(np.ma.asarray(a, dtype=float), np.nan), dtype=float).tolist())


def _kind(t, ax, fig):
    tr = t.get_transform()
    if tr == ax.transData:
        return "data"
    if tr == ax.get_yaxis_transform():
        return "yaxis"
    if tr == ax.get_xaxis_transform():
        return "xaxis"
    if tr == ax.transAxes:
        return "axes"
    if fig is not None and tr == fig.transFigure:
        return "figure"
    return "other"


def _is_letter(t):
    # Panel letters: a single bold capital. v12 drew some on the axes and some on the figure, so they
    # are compared as one set per figure rather than by the artist that holds them.
    return len(t.get_text()) == 1 and t.get_text() in "ABCDEF" and t.get_fontweight() in ("bold", 700)


def dump(fig):
    """Plotted content of a figure. Returns (gated, layout): lists of (key, value-string)."""
    fig.canvas.draw()
    gated, layout = [], []
    letters = sorted(t.get_text() for t in fig.texts if _is_letter(t))
    letters += [t.get_text() for ax in fig.axes for t in ax.texts if _is_letter(t)]
    for i, ax in enumerate(fig.axes):
        p = f"ax{i}"
        gated += [(f"{p}.xlim", _arr(ax.get_xlim())), (f"{p}.ylim", _arr(ax.get_ylim())),
                  (f"{p}.xticks", _arr(ax.get_xticks())), (f"{p}.yticks", _arr(ax.get_yticks())),
                  (f"{p}.xticklabels", json.dumps([t.get_text() for t in ax.get_xticklabels()])),
                  (f"{p}.yticklabels", json.dumps([t.get_text() for t in ax.get_yticklabels()])),
                  (f"{p}.xlabel", json.dumps(ax.get_xlabel())), (f"{p}.ylabel", json.dumps(ax.get_ylabel())),
                  (f"{p}.title", json.dumps(ax.get_title()))]
        for j, ln in enumerate(ax.get_lines()):
            gated.append((f"{p}.line{j}.xy", _arr(ln.get_xydata())))
            layout.append((f"{p}.line{j}.style", f"lw={ln.get_linewidth()} ms={ln.get_markersize()} "
                                                 f"mew={ln.get_markeredgewidth()}"))
        for j, c in enumerate(ax.collections):
            gated.append((f"{p}.coll{j}.offsets", _arr(c.get_offsets())))
            if hasattr(c, "get_segments"):
                gated.append((f"{p}.coll{j}.segments", json.dumps([json.loads(_arr(s)) for s in c.get_segments()])))
            layout.append((f"{p}.coll{j}.style", f"lw={list(np.atleast_1d(c.get_linewidths()))} "
                                                 f"sizes={list(np.atleast_1d(c.get_sizes())) if hasattr(c, 'get_sizes') else ''}"))
        for j, pa in enumerate(ax.patches):
            geo = dict(type=type(pa).__name__)
            if hasattr(pa, "center"):
                geo.update(center=list(map(float, pa.center)), radius=float(getattr(pa, "radius", np.nan)))
            gated.append((f"{p}.patch{j}", json.dumps(geo)))
        for j, t in enumerate([t for t in ax.texts if not _is_letter(t)]):
            k = _kind(t, ax, fig)
            gated.append((f"{p}.text{j}.string", json.dumps(t.get_text())))
            x, y = t.get_position()
            if k == "data":
                gated.append((f"{p}.text{j}.xy_data", _arr([x, y])))
            elif k == "yaxis":
                gated.append((f"{p}.text{j}.y_data", _arr([y])))
                layout.append((f"{p}.text{j}.x_axes", _arr([x])))
            else:
                layout.append((f"{p}.text{j}.pos_{k}", _arr([x, y])))
            if hasattr(t, "xy"):  # annotation: the point it refers to is data
                gated.append((f"{p}.text{j}.annot_xy", _arr(t.xy)))
            layout.append((f"{p}.text{j}.fontsize", str(t.get_fontsize())))
        legends = [a for a in ax.artists if isinstance(a, matplotlib.legend.Legend)]
        if ax.get_legend() is not None:
            legends.append(ax.get_legend())
        for j, lg in enumerate(legends):
            gated.append((f"{p}.legend{j}.texts", json.dumps([t.get_text() for t in lg.get_texts()])))
    gated.append(("panel_letters", json.dumps(sorted(letters))))
    gated.append(("fig.texts", json.dumps(sorted(t.get_text() for t in fig.texts if not _is_letter(t)))))
    for j, lg in enumerate(fig.legends):
        gated.append((f"fig.legend{j}.texts", json.dumps([t.get_text() for t in lg.get_texts()])))
    layout.append(("fig.size_in", _arr(fig.get_size_inches())))
    return gated, layout


def style_audit(fig):
    """Font families and sizes of every visible text, and line widths and marker sizes of every line."""
    from matplotlib.lines import Line2D
    from matplotlib.text import Text
    fig.canvas.draw()
    texts = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    letters = [t for t in texts if _is_letter(t)]
    other = [t for t in texts if not _is_letter(t)]
    lines = [ln for ln in fig.findobj(Line2D) if ln.get_visible() and ln.axes is not None]
    data_lines = [ln for ln in lines if ln.get_linestyle() not in ("None", "none", "")]
    markers = [ln for ln in lines if ln.get_marker() not in (None, "None", "none", "", " ")]
    from matplotlib.font_manager import FontProperties, findfont
    return dict(
        font_files=sorted({Path(findfont(FontProperties(family=t.get_fontfamily()))).name for t in texts}),
        text_sizes_pt=sorted({float(t.get_fontsize()) for t in other}),
        letter_sizes_pt=sorted({float(t.get_fontsize()) for t in letters}),
        letters_bold=all(t.get_fontweight() in ("bold", 700) for t in letters),
        line_widths_pt=sorted({round(float(ln.get_linewidth()), 3) for ln in data_lines}),
        marker_sizes_pt=sorted({round(float(ln.get_markersize()), 3) for ln in markers if ln.get_linestyle()
                                in ("None", "none", "") or ln.get_marker() != ""}),
        spine_widths_pt=sorted({round(float(sp.get_linewidth()), 3) for ax in fig.axes
                                for sp in ax.spines.values() if sp.get_visible()}),
    )


def capture(mod):
    rec = {}

    def fake_save(fig, stem, data):
        gated, layout = dump(fig)
        rec[stem] = dict(gated=gated, layout=layout, data=pd.DataFrame(data), style=style_audit(fig))
        plt.close(fig)

    mod.save = fake_save
    mod.fig1()
    mod.fig2()
    d = mod._load_manuscript()
    mod.fig3(d)
    mod.fig4(d)
    mod.fig5()
    mod.figS1()
    return rec


def main():
    ref = sys.argv[1] if len(sys.argv) > 1 else "77b74e0"
    src = subprocess.run(["git", "show", f"{ref}:src/figures_final.py"], cwd=ROOT, check=True,
                         capture_output=True).stdout
    with tempfile.TemporaryDirectory() as tmp:
        ref_path = Path(tmp) / "figures_final_ref.py"
        ref_path.write_bytes(src)
        before = capture(load_module(ref_path, "figures_final_ref"))
    after = capture(load_module(ROOT / "src" / "figures_final.py", "figures_final_new"))

    rows, report = [], [f"Reference: src/figures_final.py at {ref}. Compared with the working tree.", ""]
    ok_all = True
    for stem in FIGS:
        b, a = before[stem], after[stem]
        gb, ga = dict(b["gated"]), dict(a["gated"])
        keys_same = list(gb) == list(ga)
        diff = sorted(k for k in set(gb) | set(ga) if gb.get(k) != ga.get(k))
        n_num = sum(1 for k in gb if not k.endswith(("labels", "label", "title", "string", "texts")))
        tb, ta = b["data"], a["data"]
        tab_same = tb.shape == ta.shape and list(tb.columns) == list(ta.columns) and tb.equals(ta)
        csv_old = (ROOT / "results" / "figures_final" / f"{stem}_data.csv").read_bytes()
        csv_new_path = OUT / f"{stem}_data.csv"
        csv_new = csv_new_path.read_bytes() if csv_new_path.exists() else b""
        csv_same = csv_old == csv_new
        ok = keys_same and not diff and tab_same and csv_same
        ok_all &= ok
        rows.append(dict(figure=stem, artist_items_compared=len(gb), numeric_items=n_num,
                         text_items=len(gb) - n_num, artist_items_differing=len(diff),
                         plotting_table_rows=len(tb), plotting_table_identical=tab_same,
                         data_csv_bytes_identical_to_figures_final=csv_same, identical=ok))
        report.append(f"{stem}: {'IDENTICAL' if ok else 'DIFFERS'} "
                      f"({len(gb)} artist items, {len(tb)} table rows)")
        for k in diff:
            report.append(f"  DIFF {k}\n    before {gb.get(k)}\n    after  {ga.get(k)}")
        lb, la = dict(b["layout"]), dict(a["layout"])
        moved = sorted(k for k in set(lb) | set(la) if lb.get(k) != la.get(k))
        report.append(f"  layout items changed (not gated): {len(moved)}")
        for k in moved:
            report.append(f"    {k}: {lb.get(k)} -> {la.get(k)}")
    report.insert(1, f"ALL FIGURES IDENTICAL IN PLOTTED DATA: {ok_all}")
    pd.DataFrame(rows).to_csv(OUT / "check_data_identity.csv", index=False)

    from PIL import Image
    prow = []
    for stem in FIGS:
        im = Image.open(OUT / f"{stem}.png")
        w, h = im.size
        dpi = im.info.get("dpi", (np.nan, np.nan))
        lim = MAX_H_IN[ROWS[stem]]
        good = (w == round(TEXT_WIDTH_IN * 300) and h <= round(lim * 300) and abs(dpi[0] - 300) < 0.1
                and abs(dpi[1] - 300) < 0.1)
        ok_all &= good
        prow.append(dict(figure=stem, rows=ROWS[stem], width_px=w, height_px=h, dpi=round(dpi[0], 2),
                         width_in=w / 300, height_in=round(h / 300, 4), height_limit_in=lim, passes=good))
    pd.DataFrame(prow).to_csv(OUT / "check_png_size.csv", index=False)
    srows = [dict(figure=stem, **{k: (" ".join(map(str, v)) if isinstance(v, list) else v)
                                  for k, v in after[stem]["style"].items()}) for stem in FIGS]
    pd.DataFrame(srows).to_csv(OUT / "check_style.csv", index=False)
    report += ["", "Style audit of the new figures (every visible text and line; scatter sizes are in "
               "check_report layout lines):"]
    report += [f"  {r['figure']}: " + "; ".join(f"{k}={v}" for k, v in r.items() if k != "figure") for r in srows]
    report += ["", "PNG size (text width 7.5 in = 2250 px at 300 dpi):"]
    report += [f"  {r['figure']}: {r['width_px']} x {r['height_px']} px, {r['dpi']} dpi, "
               f"{r['width_in']:.3f} x {r['height_in']:.3f} in, limit {r['height_limit_in']} in, "
               f"{'pass' if r['passes'] else 'FAIL'}" for r in prow]
    (OUT / "check_report.txt").write_text("\n".join(report) + "\n", encoding="utf-8")
    print("\n".join(report[:2] + [r for r in report if not r.startswith(("    ", "  layout"))][2:]))
    if not ok_all:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

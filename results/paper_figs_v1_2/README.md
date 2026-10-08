# Preprint v1.2 figures: uniform size and style

A layout task, not an analysis. No statistic was recomputed and no plotted value changed.
The v12 figures in `results/figures_final/` are untouched. This directory holds the v1.2 figures.

## Target

| | |
|---|---|
| Text width of `paper/rejuvenation_readouts_preprint_v1.2_draft.docx` | 7.5 in: US Letter, 12240 twips wide, left and right margins 720 twips each, so 10800 twips. All six figures sit in single-column sections. |
| Figure width | 7.5 in = 2250 px at 300 dpi, inserted at 100% |
| Height limit | 3.0 in for one row of panels, 4.5 in for two |
| Font | Arial everywhere. Every label, tick label, title, legend entry and annotation is 7 pt. Panel letters are 9 pt bold. |
| Lines | 0.6 pt for spines, ticks and reference lines (zero lines, dashed guides, the circle). 1.0 pt for data lines and error bars. |
| Markers | 4 pt for every data marker (scatter `s = 16`). Open-marker edges are 0.8 pt. Error-bar caps are 2 pt. |

| Figure | Rows | Size (in) | v12 size (in) |
|---|---|---|---|
| Fig1 | 1 | 7.5 × 2.8 | 4.8 × 3.0 |
| Fig2 | 1 | 7.5 × 3.0 | 7.0 × 3.1 |
| Fig3 | 2 | 7.5 × 4.5 | 7.0 × 5.6 |
| Fig4 | 1 | 7.5 × 2.8 | 7.0 × 3.1 |
| Fig5 | 1 | 7.5 × 3.0 | 7.0 × 3.0 |
| FigS1 | 2 | 7.5 × 4.5 | 7.0 × 5.0 |

Figure 3 keeps its arrangement: A (transfer forest with its n and ρ columns) on top, B and C side
by side below. It was compressed from 5.6 in to 4.5 in. Letters A–C mean what they meant before.

## Files

- `Fig*.png` (300 dpi), `Fig*.pdf` and `Fig*_data.csv`: written by `python src/figures_final.py`.
- `check_report.txt`, `check_data_identity.csv`, `check_png_size.csv`, `check_style.csv`: written by
  `python src/paper_figs_v1_2_check.py`. It reruns every figure function of `src/figures_final.py`
  twice: as committed at 77b74e0 (the v12 figures) and as in this commit. Each run records what
  every figure plots, read back from the matplotlib artists in data coordinates. That covers every
  line's x/y, every scatter offset and error-bar segment, patches, axis limits, tick positions and
  tick labels, axis labels, titles, legend entries, annotation strings and anchors, and panel letters.
  It also records the plotting table each figure hands to `save`. The check then compares the two runs.
  - Result: all six figures are identical in plotted data, and every `_data.csv` is byte-identical
    to its `results/figures_final/` counterpart.
  - The check can see changes. An intermediate version changed Figure 3C's automatic tick positions,
    and the check flagged it; the step was then pinned (see deviation 4 below).
  - Before any edit, the committed script reproduced all six stored `_data.csv` tables byte for
    byte, and two runs compared identical.
- `docx_check.txt`: comparison of `paper/rejuvenation_readouts_preprint_v1.2.docx` with the draft,
  plus facts about the PDF. `paper/` is git-ignored, so this file is the tracked record. The
  docx/PDF helper scripts were not committed.

## Decisions the task did not fix, and deviations

1. **Output directory.** `OUT` in `src/figures_final.py` now points to `results/paper_figs_v1_2/`,
   and the docstring says so. Rerunning the script at this commit no longer rewrites
   `results/figures_final/`; those files stay reproducible from commit 77b74e0. The docx figures
   came from `src/figures_final.py`, not from `src/paper_figs_manuscript.py`, which was not changed.
2. **Marker sizes: three exceptions to "one size".**
   - The target star in the plane panels (Fig1B, Fig2A/B) is 6 pt. A 4 pt star reads smaller than a
     4 pt circle.
   - The 400-point background cloud of random gene sets in Fig2A/B is 2 pt. Its legend symbols are
     4 pt, because 2 pt symbols cannot be read; v12 also enlarged them in the legend.
   - v12 used per-shape sizes in Fig5B (3.8–4.4) to equalise visual area. All shapes are now 4 pt.
   - Points that v12 emphasised by size (Fig1B forward 5.5 pt, Fig2 real gene sets 5.5 pt, Fig3C
     open circles s = 28) are now 4 pt like everything else. They stay distinguished by colour or
     fill.
3. **Text size.** The task fixed axis and tick labels at 7 pt. Every other text (titles, legends,
   in-panel annotations, the table columns in Fig3A, the CI text in Fig4B) was also set to 7 pt.
   In v12 these ranged from 5.8 to 9 pt.
4. **Fig3C y-ticks.** At the new panel height, matplotlib's automatic locator would have switched
   from steps of 2.5 to steps of 5. The step is pinned at 2.5 (`MultipleLocator(2.5)`), which
   reproduces v12's tick positions and labels exactly. The axis limits are unchanged.
5. **FigS1D value labels.** At 7 pt in a shorter panel, the labels of folds 1, 3 and 4 overlapped;
   in v12 they already touched at 6.25 pt. Each label keeps its v12 anchor (x + 0.03, cosine) and
   string. Only the vertical alignment changed: folds 1 and 4 sit above the anchor, fold 3 below.
6. **Panel letters** are now placed by one rule in every figure: at the left edge of the panel's
   labelled area, on one baseline per row of panels. In v12, Fig1, Fig5 and FigS1 drew them on the
   axes and Fig2–4 drew them on the figure. The check compares the set of letters per figure.
7. **Guards added to `save`.** It stops if a figure is not exactly 7.5 in wide, or if anything is
   drawn outside the figure. `new_figure` stops if a height exceeds the row limit.
8. **Not changed:** colours, axis limits, every string, line breaks inside strings, and draw order.

## What the layout change does to the document (`docx_check.txt`)

- The page count goes from 13 to 14. In the PDF, page 3 now ends about one third empty, because
  Figure 2 at its larger size no longer fits below the text and moves to page 4. Fixing this needs
  a text or float change, which this task excluded.
- Word downsamples images when it exports a PDF. In `rejuvenation_readouts_preprint_v1.2.pdf` the
  figures are 200 ppi at 7.5 in, where v12 had 137–192 ppi. The docx holds the 300 dpi PNGs.

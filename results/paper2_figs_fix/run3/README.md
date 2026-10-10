# PAPER2FIGS-fix, run 3: Figs 1–4 fixed, every data file identical to 7a1efb9, row 76 MATCH

A drawing task. No plotted value changed, and no statistic was computed or refit. Row 76 is a count over rows
already plotted.

## The plan, verbatim

```
PLAN written 2026-10-10T05:51:06Z (local 2026-10-10 14:51:06 KST)
Run 3 of PAPER2FIGS-fix. Written before anything is drawn in this run and
before the new count is read. Reproduced verbatim at the top of
results/paper2_figs_fix/run3/README.md. Not edited after this timestamp.

CONTEXT. Run 1 (results/paper2_figs_fix/) STOPPED `reproduction_failed`
(raw-byte comparison; CLAUDE.md key). Run 2 (results/paper2_figs_fix/run2/)
passed reproduction 9/9 under the CRLF->LF rule, then STOPPED
`label_unplaceable` in Fig 4 (its own flag key). Neither is resumed or
edited. This run starts over from step 1.

TASK (author's words, run 1): Fig 1 arrows from the CRA box to both targets;
Fig 4 and Fig 3C factor labels moved so none overlaps another label or a
data point; one new paper2_numbers.md row ("of the 80 called against
early-passage cells (origin i), number called against GTEx young, and number
with delta_gtex_young > 0", from Fig4_data.csv with its source lines; stated
0 and 80). Commit locally as one new commit; do not push.

AUTHOR'S INSTRUCTIONS FOR RUN 3 (after run 2's STOP):
 - Fig 4 labels placed by hand at offsets written here before drawing, each
   in clear space with a leader line to its point.
 - Check: label text overlaps no data point and no other label; leader
   lines may cross points.
 - Label placement is a drawing choice, not a result. A failed label check is
   reported with the image; it is not a STOP.
 - The only hard rule: every data file identical to 7a1efb9 under the
   CRLF->LF comparison.
 - Finish Figs 1-4, add the new row, commit.

Writes only: results/paper2_figs_fix/run3/, src/paper2_figs.py,
src/paper2_figs_fix_run.py.

--- HARD RULE (the only STOP in this run besides CLAUDE.md's) ---

C (as run 2): for a written file F and blob B = `git show 7a1efb9:results/paper2_figs/F`:
  C1 F with b"\r\n" -> b"\n" equals B exactly; C2 count of b"\n" equal;
  C3 no lone b"\r" in F, no b"\r" in B. PASS = C1 and C2 and C3.
"Data files" = the 7 *_data.csv files and angle_check.csv.
S1 Reproduce PAPER2FIGS as committed (FIX off) in a temporary directory inside
   run3/ (deleted afterwards); apply C to the 7 *_data.csv, angle_check.csv
   and paper2_numbers.md. Any FAIL: STOP `reproduction_failed` (CLAUDE.md
   reproduce-first rule).
S2 Re-render Fig1-Fig4 and FigS1-FigS3 with FIX on into run3/.
G1 Apply C to the 8 data files of S2. Any FAIL: STOP `data_changed`.

--- FIG 4 LABELS, fixed here (drawing choice) ---

Each label is placed with its text anchor at the data coordinate given,
ha = "left", va = "center", textcoords = data, with a straight leader
(lw = 0.8 x LW_AXIS, shrinkB = 1.5 pt) from the text to the point. Chosen
from the committed 7a1efb9 Fig4.png (the axes run about x -23..41,
y -2.8..23.5), in regions with no points there:
  CRI_ZFX   point (-20.0, 15.3)  text at (-22.0, 21.0)   above the cluster, left
  CRA_EZH2  point (-13.0,  9.0)  text at (-10.0, 21.0)   above the cluster, right of ZFX's label
  CRA_E2F3  point (-19.2, 11.5)  text at (-22.0,  3.2)   empty area below-left of the cluster
  CRI_STAT3 point (-10.3,  5.3)  text at (-20.0,  0.0)   empty area below-left, under E2F3's label
(The point coordinates are approximate, read off the image. The leader goes to
the stored point.) These positions are not adjusted after drawing in this run.

--- CHECKS (reported, not STOPs) ---

L  For each label in Fig 4 and Fig 3C: the text's window extent vs (a) every
   plotted point's marker disc, radius 2.6 pt + 1.0 pt pad, (b) every other
   label's extent, and (c) for Fig 3C the legend and the rho text. Reported
   per label as counts in label_check.csv. Any non-zero count is reported in
   the README with the image. Not a STOP. Labels may lie outside the axes;
   that is reported too.
Fig 3C labels: placed by the same deterministic search as run 2 (216 fixed
   offsets, shortest leader first, inside the axes). If the search finds no
   clear offset for a label, that label keeps its 7a1efb9 offset and the
   failure is logged. Not a STOP.
G3 Fig 1 connectors and caption vs the boxes, checked numerically as in run 2
   and on the PNG; reported, not a STOP.
G4 paper2_numbers.md vs the blob after CRLF->LF: only the title line, the
   "Written by" line, the MATCH tally line and the appended row 76 may differ.
   Reported, not a STOP (it is not a data file).

--- THE NEW ROW 76 (as run 2's flag) ---

From run3/Fig4_data.csv (G1 must have passed), rows with perturbation != "all".
 R1 subset = called_early_passage_i == True; size stated 80.
 R2 count with called_gtex_young == True in the subset; stated 0.
 R3 count with delta_gtex_young > 0 (strict) in the subset; stated 80.
One row with both stated and both stored values; MATCH only if R2 = 0 and
R3 = 80; R1 is checked in the row's note. Cited as Fig4_data.csv line ranges;
underlying source lines in run3/new_row_sources.csv. Labelled as a count over
plotted rows.
Readings (exactly one is reported):
 - R1 = 80, R2 = 0, R3 = 80: "Of the 80 perturbations called against
   early-passage cells (origin i), 0 are called against GTEx young and all 80
   move away from GTEx young (delta_gtex_young > 0)."
 - otherwise: "MISMATCH: stated 0 and 80 (of 80); stored <R2> and <R3> (of <R1>)."
Chance baseline beside R2: 5th-percentile floor, 4.0 of 80 expected with no
effect. Not covered: the GTEx test transplants a single-cell displacement
onto a bulk GTEx origin (AXANG D3); the cells are a passaged neonatal line;
nothing transfers to aged donor cells.

STOP keys: reproduction_failed; data_changed.

COMMIT. One new commit carrying run 1, run 2 and run 3 directories (run 2's
partial figures stay as written, marked ungated in its README), src/paper2_figs.py
and src/paper2_figs_fix_run.py. CLAUDE.md is untracked and not part of it.
```

## The hard rule: data files unchanged (gated)

`gate_bytes.csv` lists, for each file, the raw byte count and line count of the blob and of the written
file, the SHA-256 of each, and the result of C1–C3.

- **S1, reproduction:** 9 of 9 PASS. This covers the 7 `*_data.csv` files, `angle_check.csv` and
  `paper2_numbers.md`, drawn as committed into a temporary directory that was then deleted.
- **G1, after the fixes:** 8 of 8 PASS. This covers the 7 `*_data.csv` files and `angle_check.csv`. In each
  file, the SHA-256 after converting CRLF to LF equals the 7a1efb9 blob's SHA-256. The line counts are equal
  (for example, Fig4_data.csv has 344 lines in both), every line ending is a CRLF, and no other CR is left
  after conversion.

The raw bytes differ only by one CR per line, because this checkout has `core.autocrlf=true`. That is why
run 1 stopped.

## Row 76 (gated reading)

**Of the 80 perturbations called against early-passage cells (origin i), 0 are called against GTEx young
and all 80 move away from GTEx young (delta_gtex_young > 0).** Stored R1 = 80, R2 = 0, R3 = 80. Row 76 is
the MATCH reading.

- **Chance baseline beside R2:** the GTEx call uses a 5th-percentile floor, so about 4.0 of 80 would be
  called with no effect. The stored 0 is below that.
- **This is a count over plotted rows.** The source is the line ranges of `Fig4_data.csv`, which G1 shows
  equal 7a1efb9. Each row's lines in `results/seng3/d4_table.csv` and `results/axang/d3b_table.csv` are in
  `new_row_sources.csv` (80 rows).
- **What it does not cover:** the GTEx test transplants a single-cell displacement onto a bulk GTEx origin
  (AXANG D3). The cells are a passaged neonatal fibroblast line whose "old" is passage number. Nothing here
  transfers to aged donor cells.
- **paper2_numbers.md (G4):** after converting CRLF to LF, the only differing lines are the allowed ones:
  line 1 (title), line 3 (Written by), line 7 (tally, now 76 MATCH, 0 MISMATCH, of 76), and the appended
  row 76. See `numbers_diff.csv`. Rows 1–75 are byte-identical to PAPER2FIGS's.

## Drawing checks (reported, not gates)

`label_check.csv`:

- **Fig 1 (G3):** no connector segment crosses a box it does not start or end at. The caption
  "displacements go to both targets" overlaps no box or connector. CRA leaves the top of its box and runs
  along the gap between the rows to a vertical bus at x = 4.88 in. CRI joins the bus from its right edge.
  The bus carries arrows to Target 1 and Target 2. The old CRI-only arrows are not drawn.
- **Fig 4 (by hand, at the flag's positions):** every label has 0 point discs, 0 other labels and 0 avoided
  artists overlapped, and every label lies inside the axes.
- **Fig 3C (search):** every label has 0 overlaps and lies inside the axes. The search placed all four, so
  no label fell back to its 7a1efb9 offset.

### Seen on the PNGs, not covered by the check

The check measures points and labels only, not reference lines or how a reader resolves a label. Each of
these is a drawing choice the author may want to change. Positions were fixed before drawing, so they were
not adjusted in this run.

1. **Fig 4, CRI_STAT3:** the text sits on the y = 0 zero line, which runs through the word. The flag put the
   text anchor at y = 0.0.
2. **Fig 3C, CRI_ZFX:** the text sits on the y = 0 zero line. Its leader and CRA_EZH2's cross.
3. **Fig 3C, CRA_E2F3:** the text ends right beside the CRI_ZFX point, so its short leader can be misread as
   pointing to CRI_ZFX.
4. **Fig 4 leaders are long** and cross the orange cluster, as the author allowed.

## Files

- `Fig1`–`Fig4`, `FigS1`–`FigS3` (`.png`, `.pdf`, `_data.csv`): FigS1–S3 were re-rendered unchanged in
  drawing.
- `paper2_numbers.md`, `angle_check.csv`, `gate_bytes.csv`, `label_check.csv`, `numbers_diff.csv` and
  `new_row_sources.csv`.
- Reproduce with `.venv/Scripts/python src/paper2_figs_fix_run.py`. `--run2` repeats run 2 into `run2/`
  (STOPs at Fig 4), and `--run1` repeats run 1's raw-byte comparison (prints only).
  `python src/paper2_figs.py` still draws PAPER2FIGS as committed, because `FIX` defaults to off.

## Deviations and decisions beyond the instruction

1. **Three runs, not one.**
   - Run 1 (`../README.md`) STOPPED with `reproduction_failed`, a CLAUDE.md key, not one from its flag. Its
     raw-byte comparison could not pass under autocrlf.
   - Run 2 (`../run2/README.md`) passed reproduction, then STOPPED with `label_unplaceable`, its own flag's
     key.
   - Neither run was resumed. Their partial outputs are committed as written: run 2's Fig1–3 are ungated.
2. **Directory.** Per CLAUDE.md and the author's answer, nothing in `results/paper2_figs/` was touched. This
   directory's paper2_numbers.md is a full copy (rows 1–75 plus row 76), not an edit of the PAPER2FIGS file.
3. **`src/paper2_figs.py` was modified** (it is PAPER2FIGS's script). The changes are the module switches
   `FIX`, `FIG4_MANUAL` and `LABEL_FALLBACK`, all off by default, plus `place_labels`, `check_labels` and
   the Fig 1 bus. S1 confirms that the default path still reproduces PAPER2FIGS's data files.
4. **Fig 1** gained a caption and a bus layout. The task asked only for CRA arrows "as the CRI box has". The
   CRI box's own two arrows were replaced by the shared bus, so that the CRA path does not cross them.
5. **The Fig 4 label positions were read off the 7a1efb9 PNG by eye** before drawing, as the flag states.
6. **The point-disc radius in the label check is 2.6 pt + 1.0 pt pad,** set in the run 2 flag and chosen
   by me, not the author.

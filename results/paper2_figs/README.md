# PAPER2FIGS: figures for the second paper

A drawing task. No statistic was computed and nothing was refit. Every plotted value is read from a
stored file in `results/senginv/`, `results/seng2/`, `results/seng3/`, `results/seng4/` or `results/axang/`,
and each figure's `_data.csv` names the source file and line of every plotted value.

## The plan, verbatim

From `PLAN.flag`. It was written at 2026-10-10T04:35:18Z, before `src/paper2_figs.py` existed, and has not been edited since.

```
PLAN written 2026-10-10T04:35:18Z (local 2026-10-10 13:35:18 KST)
Written before src/paper2_figs.py exists and before any figure or check is run.
Reproduced verbatim at the top of results/paper2_figs/README.md. Not edited
after this timestamp.

TASK: PAPER2FIGS. Drawing task only: no new statistic, no recomputation, no
refitting. Every plotted value is read from a stored result file; if one is
not stored, STOP `value_not_stored` and list it rather than compute it.

Sources (read only): results/senginv/, results/seng2/, results/seng3/,
results/seng4/, results/axang/. Any C1_lu angle comes from
results/axang/d1_angles_corrected_C1lu.csv, never the superseded d1_angles.csv
rows (results/axang/CORRECTION.md).
Writes only: src/paper2_figs.py, results/paper2_figs/.
Style: the constants and helpers of src/figures_final.py (results/paper_figs_v1_2/),
imported, not copied: 7.5 in wide, <= 3.0 in per row / 4.5 in for two rows,
Arial 7 pt, letters 9 pt bold, same line widths, marker size and palette.

--- ANSWERS FROM THE AUTHOR, before any code (verbatim) ---

Q1 Fig 2B full-gene rows: "seng3 only". Both spaces are drawn from
   results/seng3/d2_plant.csv (S1 = full-gene, S2A_k20) at 96/107/128/157
   cells. results/seng2/d3_plant_summary.csv is not drawn.

Q2 Angle from parallel: "Convert from the CSV: plot min(θ, 180−θ) from
   d1_angles.csv at full precision, and disclose it as a display conversion
   (amendment F10's rule). For C1_lu, read angle_from_parallel_corrected
   straight from the corrected CSV, never the superseded rows. Check every
   plotted value against the FINDINGS_AXANG.md §5 table to 2 decimal places,
   and list any value that doesn't match in paper2_numbers.md."

--- DECISIONS THE TASK DID NOT FIX, fixed here before running ---

P1 Origins. Where a stored file holds both origins, both are drawn: origin
   (i) filled markers, origin (ii) open markers. Neither is privileged in the
   drawing; (i) is the origin the SENG3/SENG4 gates were read on and is said
   so in the README. Fig 3C uses delta_i only, because the stored rho
   (-0.8639) is for delta_i. S1 draws origins (i) and (ii) as two panels.
P2 Fig 2A. y = null delta median with its stored 2.5-97.5% range, by cell
   count, per library (marker shape) and origin (fill), colour = space.
P3 Fig 2B. y = power, at doses 0.25 / 0.50 / 1.00, one dot per
   library x cell count, colour = space, fill = origin. The 0% dose is not
   drawn (it is the noise floor, 0.05 or 0.00 by construction) and is
   written to the data CSV.
P4 Fig 2C. Grid of the four factors x {full-gene (i), full-gene (ii),
   S2A_k20 (i), S2A_k20 (ii)}; filled = called, open = not called. Full-gene
   from results/seng2/d2_table.csv called_delta_i/ii, S2A_k20 from
   results/seng3/d4_table.csv called_delta_i/ii. Each column carries the
   stored screen-wide called count of 340 and the chance expectation 17.0
   (CLAUDE.md: chance baseline beside any count).
P5 Fig 3A / S3 rows and columns. Space S2A_k20 only. 3A: A1-A4 x {C1_south,
   C1_seng_CRA, C1_seng_CRI, C1_lu (corrected)}. S3: A1-A4 x {C2, C1_gtex}
   (the 652-donor primary C1_gtex; the 344-donor anchors variant and the
   duplicate C2_seng4_P3 rows are not drawn).
P6 Fig 3B. Power at the 0.50 dose from results/seng4/d2_plant_residualized.csv
   for axis NONE, P1, P2, P3, one dot per library x cell count, origin (i)
   filled and (ii) open; SENG3's stored unresidualized power
   (column seng3_power_unresidualized, origin (i)) drawn as a reference group.
P7 Fig 3C. x = prolif_shift, y = delta_i, 340 rows of
   results/seng3/d4_table.csv; colour = called_delta_i; the four labelled.
   The stored rho is printed in the panel, read from
   results/seng3/d5_confounds.json, not recomputed.
P8 Fig 4. x = delta_i from results/seng3/d4_table.csv, y = delta from
   results/axang/d3b_table.csv, inner join on `perturbation`; the join must
   give 340 rows with matching modality and n_cells, else STOP
   `join_mismatch`. Marker = called on each test (both / early-passage only /
   GTEx only / neither). Stored counts printed with chance 17.0.
P9 S1. y = ruler_change, x = delta_i (panel A) and delta_ii (panel B), from
   results/seng3/d4_table.csv; colour = ruler_scores_younger x called.
P10 S2. (A) library offsets from results/seng3/d3_offset.csv, variant half_A,
   thinned medians with their stored 2.5-97.5% ranges: offset/||a|| and
   cos(offset, a), per space and library. (B) cells per perturbation by
   library, histogram of n_cells in results/seng3/d4_table.csv (binning is
   display only). (C) ruler-gene coverage present/absent per dataset:
   Sengstack and Southard from results/axang/d2_summary.json, Lu from
   results/axang/correction_values.json (corrected), never the superseded
   Lu row.
P11 Fig 1 numbers are read from results/senginv/addendum_meta.json,
   results/senginv/ADDENDUM_SENGINV.md (group table), results/seng2/d1_summary.json
   (NT 666 / 1,479) and results/axang/stage0_summary.json (GTEx young n).
P12 paper2_numbers.md. For each listed value: the stored file and line it is
   read from (a CSV/JSON result file where one holds it; a findings or
   markdown file only where the value exists nowhere else, said so), the
   stored value, and MATCH if the stated value equals the stored value
   rounded to the decimals as stated (integers exact; percentages compared
   as percentages), else MISMATCH with both values. A stated value that is a
   ratio of stored numbers is checked by that ratio and labelled as such.
   No value is adjusted to match.
P13 Angle check (Q2). Every plotted angle is compared to 2 dp with the
   FINDINGS_AXANG.md §5 table. The four C1_lu rows are expected to differ
   there because that table holds the superseded values; they are listed as
   mismatches and additionally checked against CORRECTION.md's corrected
   table.
P14 Every figure writes FigX_data.csv holding each plotted value with its
   source file and row key.

STOP keys: value_not_stored; join_mismatch.
```

## Files

- `Fig1`–`Fig4`, `FigS1`–`FigS3`: `.png` at 300 dpi, `.pdf`, and `_data.csv` with every plotted value and
  its source file and line. All are 7.5 in wide. One-row figures are at most 3.0 in tall; Fig3 and FigS2
  have two rows and are 4.5 in.
- `paper2_numbers.md`: the 75 manuscript values the task listed, each with the stored file and line it was
  read from. **75 MATCH, 0 MISMATCH.**
- `angle_check.csv` and the second section of `paper2_numbers.md`: every plotted angle from parallel
  checked against the FINDINGS_AXANG.md §5 table to 2 dp (Q2). **20 of 24 match. The 4 that do not are the
  C1_lu rows.** §5 holds the superseded C1_lu values. All 4 plotted C1_lu values match the corrected
  S2A_k20 table in `results/axang/CORRECTION.md` (lines 77–80).
- Reproduce with `.venv/Scripts/python src/paper2_figs.py`. Two runs gave byte-identical CSV and markdown
  outputs.

## Deviations and decisions beyond the plan

1. **No FINDINGS_PAPER2FIGS.md.** This is a drawing task, so the record is this README, following
   PAPERFIGS v1.2 (`results/paper_figs_v1_2/README.md`). CLAUDE.md's layout names a `FINDINGS_<TASK>.md`;
   no gated result exists to put in one.
2. **Style is imported.** `src/paper2_figs.py` imports `src/figures_final.py` and uses its constants,
   `use_style`, `new_figure`, `axes_in`, `panel_letters` and `save`. These enforce the width, height and
   drawing-bounds guards. It redirects `figures_final.OUT` to this directory at save time.
   `src/figures_final.py` is not modified.
3. **Colours carry new meanings here.** Blue is the full-gene space and vermillion is S2A_k20, the
   compressed space. Elsewhere, vermillion also marks "called". In Fig 3A the C1 axes are told apart by
   colour and marker. The v1.2 meanings (vermillion = MD, purple = age lists) do not apply.
4. **Marker sizes.** Each figure has one marker size, with three exceptions. The 340-point scatters
   (Fig3C, Fig4, FigS1) use 3 pt markers. The power dots (Fig2B, Fig3B) are 3.4 pt so overlapping
   points stay visible. The Fig2C grid squares are 5.6 pt.
5. **Fig2B** draws no 0% dose. That dose is the noise floor (0.05 by construction in S2A_k20, 0.00 in
   full-gene). Its rows are in `Fig2_data.csv` with `drawn = False`.
6. **Fig1 numbers** include the 167/173 perturbations per library
   (`results/seng3/d4_summary.json`) and the 231 GTEx donors aged 60–79 that form the origin
   (`results/axang/stage0_summary.json`). The task did not list these; they were added to make the
   schematic complete.
7. **Fig2C** prints the stored called-of-340 count and the chance expectation of 17.0 for each column.
   Fig4 and Fig3C print stored counts the same way. Every count and the rho are read from the summary
   JSONs and are not recounted from the plotted rows.
8. **Value checks that are arithmetic on stored numbers.** Three entries in `paper2_numbers.md` are
   arithmetic on stored values and are labelled as such:
   - "45.1% of gap kept" is the ratio of two stored norms. The 45.1% itself is stated only in
     FINDINGS_SENG3.md:662.
   - "23,485" is present plus absent.
   - The ranges (for example 0.565–0.740) are the minimum and maximum of the stored rows, each cited by line.
9. **A checker bug, found and fixed before the final run.** The first version of the CORRECTION.md
   parser also read the full-gene table and took two full-gene rows for A1/A2 vs C1_lu. It reported
   those two as not matching the corrected values. The parser now stops at the full-gene heading and
   requires exactly the four S2A_k20 rows. No plotted value changed.

## What the figures do not show, and caveats that apply wherever they are used

- **Chance baseline.** Each called count rests on a 5th-percentile floor, so 17.0 of 340 is expected by
  construction. The realised effect-free firing rate in S2A_k20 was 0.040–0.085 under origin (i) and up to
  0.110 under (ii) (Fig2A data). That means up to about 29 of 340 under (i).
- **Sensitivity.** SENG3 D2 fired LIMITED: S2A_k20 detects a 50% approach but not reliably a 25% one
  (lowest powers 0.940 and 0.565, Fig2B). A full-gene zero is not a null: the full-gene test cannot fire,
  because 0 of 200 effect-free draws come out negative (Fig2A). The residualized statistic is blind
  (p50 = 0.010, Fig3B), so no residualized count anywhere is evidence of absence.
- **Cell count.** The four factors sit at the 82nd–99th percentile of cell count. In their depth band,
  24 of 60 are called (40%).
- **Proliferation.** In S2A_k20, delta (i) is rank-correlated −0.8639 with the cell-division shift
  (Fig3C), and in this dataset passage aging is 19.66° from proliferation (Fig3A, A3 row). Fig 3 should be
  read with this in mind.
- **Circularity, disclosed rather than refitted.**
  - In Fig 3A and S3, C1_lu and A4 are both built from Lu day-0 cells and overlap (AXANG F6).
  - In S3, C1_gtex's 652 donors include the 344 anchor donors that define A1 (AXANG F2).
  - C1_sengstack is built from NT half one, which no origin, target or floor uses (SENG4).
  - In Fig3C, the cell-division shift uses the same NT cells that define origin (i) (SENG3 D5).
- **Cross-space numbers are not comparable.** Every angle drawn is in S2A_k20. Fig 2's full-gene and
  S2A_k20 values come from different spaces and are shown side by side, not as one scale.
- **Fig4's GTEx axis** transplants a single-cell displacement measured on one platform onto a bulk GTEx
  origin (AXANG D3).
- **Scope.** The Sengstack cells are a passaged neonatal-line fibroblast culture, and their "old" is
  passage number. Passage aging is near-orthogonal to donor aging (cos(a, u) = −0.040273). Nothing in
  these figures transfers to aged donor cells.
- **Superseded values.** Every C1_lu value is the corrected one. The Lu coverage in FigS2D is
  20,452 / 8.97% from `results/axang/correction_values.json`. The 16,781 / 21.97% Lu row in
  `results/axang/d2_summary.json` is superseded and is not drawn. The 21.97% in `paper2_numbers.md` is the
  Sengstack value, which happens to equal that superseded Lu figure.

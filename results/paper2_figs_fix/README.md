# PAPER2FIGS-fix: STOPPED, `reproduction_failed`

**The run stopped at step 1, before any fixed figure was rendered and before the new row was counted.**
No Fig 1 / Fig 3C / Fig 4 fix exists in this directory, and no value for the new paper2_numbers.md row
was read. Per CLAUDE.md the run is not resumed or restarted with an adjustment.

## The plan, verbatim

From `PLAN.flag`, written at 2026-10-10T04:55:47Z, with decision D1 appended at 04:58Z by the author
before any code change. Not edited otherwise.

```
PLAN written 2026-10-10T04:55:47Z (local 2026-10-10 13:55:47 KST)
Written before any figure is re-rendered and before the new count is read.
Reproduced verbatim at the top of the task README. Not edited after this
timestamp; any later decision is appended below a dated line, not edited in.

TASK: PAPER2FIGS-fix. Cosmetic only. No plotted value may change.

Author's instruction (verbatim):
1. Fig 1: draw arrows from the CRA box to both targets, as the CRI box has.
   Both libraries' displacements go to both targets.
2. Fig 4 and Fig 3C: move the factor labels so none overlaps another label or
   a data point.
3. Add one row to paper2_numbers.md: "of the 80 called against early-passage
   cells (origin i), number called against GTEx young, and number with
   delta_gtex_young > 0", from Fig4_data.csv with its source lines.
   Stated: 0 and 80.
Re-render, confirm every *_data.csv is byte-identical to commit 7a1efb9, and
commit locally as one new commit. Do not push.

--- GATES ---

G1 Byte identity. Every *_data.csv written by the re-render (Fig1, Fig2, Fig3,
   Fig4, FigS1, FigS2, FigS3) is compared byte for byte (SHA-256) with
   `git show 7a1efb9:results/paper2_figs/<name>`. Any difference: STOP
   `data_changed`, naming the files. paper2_numbers.md and angle_check.csv
   are also compared; paper2_numbers.md is expected to differ only by the
   added row(s) and the MATCH tally line, angle_check.csv not at all.
G2 Overlap (item 2). Checked from the rendered figure, not by eye alone: for
   each of the four labels in Fig 4 and in Fig 3C, the label's text bounding
   box (renderer window extent) must not intersect another label's box, nor
   the marker extent of any plotted point in that panel (marker radius in
   display units), nor the panel's legend or in-panel text. Leader lines may
   cross points. A failure is fixed by moving labels and re-checking; it is
   not a STOP. The final placement and the check result are reported.
G3 Fig 1 (item 1). Arrows from the CRA box to Target 1 and to Target 2, not
   crossing the CRI box or any text. Checked on the rendered PNG.

--- THE NEW ROWS (item 3) ---

Read from the re-rendered Fig4_data.csv (which G1 requires to equal the
7a1efb9 copy), rows with perturbation != "all".
 R1 subset = called_early_passage_i == True. Its size must be 80, else
    MISMATCH is reported on a row "called against early-passage (i)", stated 80.
 R2 "of those 80, called against GTEx young": count of called_gtex_young ==
    True in the subset. Stated 0. Exact integer.
 R3 "of those 80, delta_gtex_young > 0": count of delta_gtex_young > 0
    (strictly) in the subset. Stated 80. Exact integer.
Each row cites Fig4_data.csv with the line range of the subset, plus the
underlying source lines those rows carry (results/seng3/d4_table.csv and
results/axang/d3b_table.csv, from the `source` column). Rule P12 of the
PAPER2FIGS plan applies: MATCH only if equal, else MISMATCH with both values;
nothing adjusted. These are counts over plotted rows (arithmetic on stored
values) and are labelled as such. The instruction says "one row"; it names two
numbers, so it is written as one table row holding both stated and both stored
values, with MATCH only if both agree. (R1 is checked inside that row's note.)

Readings (exactly one is reported):
 - both match: "Of the 80 perturbations called against early-passage cells
   (origin i), 0 are called against GTEx young and all 80 move away from GTEx
   young (delta_gtex_young > 0)." [with the stored numbers]
 - otherwise: "MISMATCH: stated 0 and 80; stored <a> and <b>." The figure
   commit still proceeds; the mismatch is reported, not corrected.
Chance baseline beside the count: the GTEx test's called rate under an
effect-free null is 5% by construction (chance expectation 17.0 of 340, so
4.0 of 80). This is stated beside R2.

STOP keys: data_changed; value_not_stored.

--- 2026-10-10T04:58Z, decided by the author before any code change ---
D1 Directory. CLAUDE.md forbids writing into results/paper2_figs/. Author's
   choice: every output of this task goes to results/paper2_figs_fix/,
   including a full paper2_numbers.md (the 75 PAPER2FIGS rows plus the new
   row). results/paper2_figs/ is not touched. G1 compares the new *_data.csv
   files against 7a1efb9:results/paper2_figs/.
```

## What happened

`src/paper2_figs_fix_run.py` step 1 drew PAPER2FIGS exactly as committed (`paper2_figs.FIX = False`) into a
temporary directory under this one (deleted on exit), then compared each output byte for byte with
`git show 7a1efb9:results/paper2_figs/<file>`. These are the plan's G1 terms, and step 1 applies them to the
reproduction. **All 9 files differed** (`gate_bytes.csv`):

| file | bytes at 7a1efb9 (blob) | bytes now | identical |
|---|---|---|---|
| Fig1_data.csv | 838 | 852 | no |
| Fig2_data.csv | 15835 | 16016 | no |
| Fig3_data.csv | 43573 | 44011 | no |
| Fig4_data.csv | 45631 | 45975 | no |
| FigS1_data.csv | 64559 | 65244 | no |
| FigS2_data.csv | 22340 | 22692 | no |
| FigS3_data.csv | 1111 | 1120 | no |
| paper2_numbers.md | 13136 | 13264 | no |
| angle_check.csv | 3776 | 3801 | no |

## Diagnostic (not a gated result)

Read-only checks after the STOP, with nothing recomputed:

- This checkout has `core.autocrlf=true`. `git ls-files --eol` shows every PAPER2FIGS CSV and
  `paper2_numbers.md` stored as `i/lf` in the repository and checked out as `w/crlf`.
- On Windows, the script writes CRLF. Each byte excess equals the file's line count (Fig1_data.csv: 14 lines,
  +14 bytes).
- For Fig1_data.csv, Fig4_data.csv and paper2_numbers.md, the **checked-out** files in `results/paper2_figs/`
  equal the 7a1efb9 blobs exactly once CRLF is converted to LF. The new outputs were not compared
  under that conversion. Doing so would be the adjusted re-run that CLAUDE.md forbids without a decision.

The likely cause is therefore line endings, not content. That is not established for the new outputs.

## Why this is a STOP and not a fix

The task said "byte-identical to commit 7a1efb9". With autocrlf on, that has two readings. One is the
repository blob (LF), which a Windows render can never match byte for byte. The other is the
checked-out file (CRLF), which is what git would write. The flag fixed the first reading. Choosing the
second now, to let the run proceed, is what CLAUDE.md says not to do. The author decides.

## Not done

- Fig 1 arrows, Fig 3C / Fig 4 label placement: code is written (`src/paper2_figs.py` `FIX` branches and
  `place_labels`), never rendered.
- paper2_numbers.md row 76: code written (`new_row` in `src/paper2_figs_fix_run.py`), never run.
- No commit.

---

## Addendum 2026-10-10T05:47Z, appended after the STOP. Nothing above is edited.

- **Where the STOP key came from.** `reproduction_failed` is not one of the flag's STOP keys (`data_changed`,
  `value_not_stored`). It comes from CLAUDE.md's standing rule: reproduce the stored numbers first, and if
  reproduction fails, STOP. The flag's G1 terms (raw bytes against the `git show` blob) were applied to that
  reproduction step.
- **Author's decision after the STOP:** use the first reading. Compare after converting CRLF to LF, as git
  does on commit with autocrlf, and also require an exact line-count match. That rule is pre-registered in a
  new run, `results/paper2_figs_fix/run2/`, with its own flag. This attempt is not resumed; it stays on
  record as it is.

## Addendum 2026-10-10T05:53Z, appended. Nothing above is edited.

- **Run 2** (`run2/`) passed reproduction 9/9 under the CRLF->LF rule, then STOPPED `label_unplaceable` in Fig 4.
- **Run 3** (`run3/`) is the completed run: Figs 1–4 fixed; all 8 data files identical to 7a1efb9 under
  CRLF->LF with exact line counts; row 76 MATCH (0 and 80). Findings: `run3/README.md`.

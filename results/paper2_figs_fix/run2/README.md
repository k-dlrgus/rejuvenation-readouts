# PAPER2FIGS-fix, run 2: STOPPED, `label_unplaceable`

**The run stopped while drawing Fig 4. It had passed reproduction (S1), but no gate after S1 was
evaluated.** G1, G2, G3 and G4 were not run, and row 76 of paper2_numbers.md was not counted. Per CLAUDE.md
the run is not resumed, and the label search is not adjusted and rerun. There is no commit.

The STOP key is one of this run's own flag keys. It is not a CLAUDE.md key.

## The plan, verbatim

```
PLAN written 2026-10-10T05:47:05Z (local 2026-10-10 14:47:05 KST)
Run 2 of PAPER2FIGS-fix. Written before anything is drawn in this run, before
the comparison below is applied to any file, and before the new count is read.
Reproduced verbatim at the top of results/paper2_figs_fix/run2/README.md. Not
edited after this timestamp; any later decision is appended below a dated line.

CONTEXT. Run 1 (results/paper2_figs_fix/, flag 2026-10-10T04:55:47Z) STOPPED at
reproduction, key `reproduction_failed` (from CLAUDE.md's reproduce-first rule,
not from that flag). Its comparison was raw bytes against the `git show`
blob; this checkout has core.autocrlf=true, blobs are LF, files written here
are CRLF. Run 1 is not resumed and its files are not edited (a dated
addendum was appended to its README). This run starts over from step 1.

TASK (unchanged from run 1, author's words):
1. Fig 1: draw arrows from the CRA box to both targets, as the CRI box has.
   Both libraries' displacements go to both targets.
2. Fig 4 and Fig 3C: move the factor labels so none overlaps another label or
   a data point.
3. Add one row to paper2_numbers.md: "of the 80 called against early-passage
   cells (origin i), number called against GTEx young, and number with
   delta_gtex_young > 0", from Fig4_data.csv with its source lines.
   Stated: 0 and 80.
Cosmetic only; no plotted value may change. Commit locally as one new commit.
Do not push.

Writes only: results/paper2_figs_fix/run2/, src/paper2_figs.py (FIX switch,
default off = PAPER2FIGS as committed), src/paper2_figs_fix_run.py.

--- THE COMPARISON RULE (author's decision after run 1) ---

C  For a file F written here and its blob B = `git show 7a1efb9:results/paper2_figs/F`:
   C1 EOL-normalized identity: bytes(F) with every b"\r\n" replaced by b"\n"
      equals bytes(B) exactly. (This is how git stores the file on commit
      with autocrlf.)
   C2 Line count: the number of b"\n" in F equals the number in B, exactly.
   C3 No other CR: F holds no b"\r" that is not followed by b"\n", and B
      holds no b"\r" at all. (So the conversion removed only line-ending CRs.)
   A file PASSES only if C1, C2 and C3 all hold. Each is reported per file
   with byte counts, line counts and SHA-256 of raw and normalized bytes.

--- STEPS AND GATES ---

S1 Reproduce. Draw PAPER2FIGS as committed (FIX = False) into a temporary
   directory inside run2/ (deleted afterwards) and apply C to the 7 *_data.csv
   files, paper2_numbers.md and angle_check.csv. Any FAIL: STOP
   `reproduction_failed`, naming every file and which of C1-C3 failed.
S2 Re-render with FIX = True into run2/.
G1 Apply C to the 7 *_data.csv files of S2 and to angle_check.csv. Any FAIL:
   STOP `data_changed`.
G2 Label overlap, from the rendered figure: for each of the four factor
   labels in Fig 3C and in Fig 4, the label text's window extent must not
   intersect another label's, nor any plotted point's marker disc (radius
   2.6 pt + 1.0 pt pad, covering the 4 pt open highlight marker and the 3 pt
   scatter markers), nor (Fig 3C) the legend or the rho text. Labels must lie
   inside the axes. Leader lines may cross points. Placement is a
   deterministic search over fixed candidate offsets (shortest leader first).
   If no candidate clears a label: STOP `label_unplaceable`. Result per label
   in label_check.csv.
G3 Fig 1: the CRA and CRI connectors reach both targets; no connector or the
   connector caption crosses a box or other text. Checked numerically (caption
   extent vs box extents; connector coordinates vs box extents) and on the
   rendered PNG; reported.
G4 paper2_numbers.md vs the blob, after C1-style EOL normalization: the only
   differing lines may be (a) the title line (PAPER2FIGS -> PAPER2FIGS-fix),
   (b) the "Written by" line (now names src/paper2_figs_fix_run.py and says
   row 76 is new), (c) the MATCH tally line, and (d) the appended row 76.
   Every other line identical, in order. Anything else: STOP `numbers_changed`.

--- THE NEW ROW (item 3), as run 1's flag ---

From run2/Fig4_data.csv (which G1 requires to equal the 7a1efb9 copy), rows
with perturbation != "all".
 R1 subset = called_early_passage_i == True; size stated 80.
 R2 count of called_gtex_young == True in the subset; stated 0.
 R3 count of delta_gtex_young > 0 (strict) in the subset; stated 80.
One table row (row 76) holds both stated and both stored values; MATCH only
if R2 and R3 both equal the stated values exactly; R1 is checked in the
row's note. Cited as Fig4_data.csv line ranges; each subset row's underlying
source lines (results/seng3/d4_table.csv, results/axang/d3b_table.csv) are
listed in run2/new_row_sources.csv. Labelled as a count over plotted rows.

Readings (exactly one is reported):
 - R1 = 80, R2 = 0, R3 = 80: "Of the 80 perturbations called against
   early-passage cells (origin i), 0 are called against GTEx young and all 80
   move away from GTEx young (delta_gtex_young > 0)."
 - otherwise: "MISMATCH: stated 0 and 80 (of 80); stored <R2> and <R3> (of <R1>)."
   Reported, not corrected; the figure commit still proceeds.
Chance baseline beside R2: the GTEx call uses a 5th-percentile floor, so 5%
of effect-free perturbations are called by construction: 4.0 of 80 expected.
What R2 = 0 does not cover: the GTEx test transplants a single-cell
displacement onto a bulk GTEx origin (AXANG D3), and the cells are a passaged
neonatal line; nothing here transfers to aged donor cells.

STOP keys: reproduction_failed; data_changed; label_unplaceable;
numbers_changed; value_not_stored.
```

## What happened

1. **S1 reproduction: 9 of 9 PASS** under rule C (`gate_bytes.csv`). For every file, the LF-normalized output
   equals the 7a1efb9 blob and the line counts are equal. Every LF in the written file is part of a CRLF,
   and no CR is left after normalization. This confirms run 1's diagnosis: run 1's mismatch was line endings
   only.

   | file | lines blob / written | bytes blob / written | C1 | C2 | C3 |
   |---|---|---|---|---|---|
   | Fig1_data.csv | 14 / 14 | 838 / 852 | ✓ | ✓ | ✓ |
   | Fig2_data.csv | 181 / 181 | 15835 / 16016 | ✓ | ✓ | ✓ |
   | Fig3_data.csv | 438 / 438 | 43573 / 44011 | ✓ | ✓ | ✓ |
   | Fig4_data.csv | 344 / 344 | 45631 / 45975 | ✓ | ✓ | ✓ |
   | FigS1_data.csv | 685 / 685 | 64559 / 65244 | ✓ | ✓ | ✓ |
   | FigS2_data.csv | 352 / 352 | 22340 / 22692 | ✓ | ✓ | ✓ |
   | FigS3_data.csv | 9 / 9 | 1111 / 1120 | ✓ | ✓ | ✓ |
   | paper2_numbers.md | 128 / 128 | 13136 / 13264 | ✓ | ✓ | ✓ |
   | angle_check.csv | 25 / 25 | 3776 / 3801 | ✓ | ✓ | ✓ |

2. **S2 re-render** wrote Fig1, Fig2 and Fig3 (`.png`, `.pdf`, `_data.csv`) to this directory. It then stopped
   in Fig4 with `STOP label_unplaceable: Fig4: no candidate offset clears CRA_E2F3`. Fig4, FigS1–FigS3,
   label_check.csv, paper2_numbers.md, angle_check.csv and new_row_sources.csv were never written.

## The files left here are not results

`Fig1.*`, `Fig2.*` and `Fig3.*` in this directory are partial S2 output, **ungated**. G1 was never applied to
their `_data.csv`. The Fig 1 connector check (G3) and the Fig 3C label check (G2) were computed in memory,
but the script exited before writing them. They are not reported and were not inspected after the STOP.

## The likely cause (not a gated result; from the code, not from new computation)

The search tries 216 fixed offsets (9 distances from 8 to 48 pt × 24 angles) and requires the label box to
lie inside the axes. CRA_E2F3 sits at (−19.2, 11.5), inside the dense orange cluster near the left edge of
the axes. The search found no candidate within 48 pt that avoids every point and stays inside the axes. The
flag fixed the search space and made exhaustion a STOP. Run 1's flag instead said a G2 failure is fixed by
moving labels; run 2's flag replaced that with a STOP for an exhausted search.

## Not done

- No Fig 4, no label or connector check reported, no row 76, no commit.

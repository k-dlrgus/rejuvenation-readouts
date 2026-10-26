# FINDINGS_CORR3 — finishing the corrections to both papers

**Task:** CORR3. **Status:** complete. No STOP fired. **Nature:** writing only. No statistic was
computed. Every number written into either manuscript or either notice is copied from
`FINDINGS_CORR2CHK.md`, `FINDINGS_FLOORCHK.md` or their `results/` tables, and the script asserts
the paragraph-142 wording is byte-identical to the blockquote it came from rather than retyped.

## Step 0 — v0.7 is the source of the posted paper-2 PDF

Checked before anything was written. `paper/two_definitions_of_young_preprint_v1.pdf` was
extracted as UTF-8 text in reading order, running headers removed, and compared with the body text
of `paper/paper2_passage_draft_v0.7.docx`.

| check | result |
|---|---|
| v0.7 non-empty body paragraphs | 99 |
| found in the PDF, verbatim | **99 of 99** |
| in the same order | **yes** |
| share of the PDF's text covered by them | 98.2% |
| Table 1 non-empty cells found in the PDF | **52 of 52** |
| paragraphs that differ | **none** |

**No paragraph differs, so v0.7 is the source and no STOP fired.** Two notes on the method, so
the result can be judged:

- **Hyphens had to be normalised out of both sides.** On a first pass 12 of 99 paragraphs appeared
  to diverge. Every one of the 12 was the same artefact: a real hyphen in the docx
  ("non-targeting", "early-passage", "pre-registered" and so on) that falls at a line break in the
  PDF, where `pdftotext` drops it. With hyphens removed from both sides all 99 match. This means
  the comparison could not have detected a change that was only a hyphen; nothing else is
  affected.
- **Table 1 was checked cell by cell, not as a block.** Its 52 cells account for the 1.8% of the
  PDF not covered by body paragraphs. One cell, the header "Cell-division direction(one per row),
  built from", does not appear as a contiguous string because it is one docx cell holding two
  paragraphs that the PDF lays out on separate lines with other header text between them; both
  halves are present separately.

## Step 1 — the grammar fix in paper 2, v0.8

`paper/paper2_passage_draft_v0.8.docx`, paragraph 30: "1 move away less far" → "1 moves away less
far". Nothing else changed.

`results/corr2chk/readings.txt` and `FINDINGS_CORR2CHK.md` were **not** touched. They carry the
original form because they are the record, and the slip is disclosed in them and in
`FINDINGS_CORR2.md` deviation 3.

## Step 2 — `CORRECTIONS_PAPER2.md`

Written in the form of `CORRECTIONS_PAPER1.md`, with one substantive difference stated in its
first paragraph: **paper 2 is not public**, so it is a pre-publication revision and not a
correction of any record. It covers paragraph 31's Southard clause (A1b) and paragraph 30's
young-donor clause (C5), keeps their two failures distinct, states prominently that the paper's
conclusion is unaffected, explains why the Summary's "all four move them away from young donors"
is a different statement that is not affected, gives the corrected wording for both, discloses how
the pre-registration behind the check is weaker than paper 1's and what still constrains it, and
appends all 21 inventoried statements with their verdicts and the three length ratios with no
size-matched floor.

## Step 3 — paper 1, paragraph 142

The replacement wording written as step 4 of `FINDINGS_CORR2CHK.md` is now in
`paper/rejuvenation_readouts_preprint_v1.3.docx`, alongside the Supplementary Note 3 change CORR1
had already made. Both of that wording's reversible judgement calls were accepted on instruction:
the `q = 1.0` is dropped, and the space B floor comparison is added.

`CORRECTIONS_PAPER1.md` changed in seven places:

1. The metadata gained a **Revised** line naming CORR3 and stating that CORR1's two sentences are
   unchanged, plus a note that both corrected passages sit in Supplementary Note 3 and fail for
   different reasons.
2. **What was stated** is now (a) the two sentences and (b) the bullet, each quoted verbatim.
3. **What is wrong** is split into "(a) the two sentences", unchanged, and a new "(b) the bullet",
   which gives the three reasons it carries no information: no perturbation's change in distance
   could have been negative in space C; in space B none of the 1,836 is negative either and these
   two sit inside their own cell-count-matched floors, as 1,088 of 1,836 do; and the quoted
   `q = 1.0` came from a permutation null SOUTH8 found anti-conservative in that space.
4. **What still holds** gained a paragraph stating that the decision the bullet supported still
   stands on the paragraph's other three diagnostics, and the corrected bullet beside the
   corrected sentences, with the two deliberate changes named.
5. The old passage recording paragraph 142 as a detail that "loses its force but is not corrected
   here" is replaced by a parenthesis saying CORR1 recorded it that way and CORR3 corrected it.
   The old wording is quoted there rather than simply removed.
6. **Where the check is recorded** is split by passage, since (b)'s numbers come from
   `FINDINGS_CORR2CHK.md` and `results/coordspace/`, not from FLOORCHK.
7. **Appendix A** gained rows 9 to 12 for the bullet, and **Appendix B** now describes four edits
   rather than three.

### One row of that appendix needs a decision outside this task

`SUPPLEMENTARY_NOTE_COORDSPACE.md:52` carries the superseded bullet verbatim, in a tracked file
whose first line reads "Coordinate-space robustness check: text for the preprint". It was **not
edited**, following CORR1's own convention and CLAUDE.md's rule that superseded text is marked
elsewhere rather than edited in place. But unlike a findings file it presents itself as preprint
copy, so the repository now holds drafted text that contradicts the corrected manuscript. It is
row 11 of Appendix A and is named in `results/verify/WITHDRAWN_READINGS.md`. Whether to revise it
is not this task's call.

## Step 4 — the repository record

- `results/verify/WITHDRAWN_READINGS.md`, section D: the paper-1 entry now covers paragraph 142
  as passage (b), with what it said, what was wrong with it, where it is corrected, why it is not
  counted as a reading, and what was not edited. A new dated paper-2 entry follows it in the same
  shape, and states in its own bullet that paper 2 is not public and that the pre-registration
  behind the check is weaker than FLOORCHK's.
- `README.md`: the paper-1 correction line now says two passages rather than one and cites both
  checks; a new line under the second preprint points at `CORRECTIONS_PAPER2.md`.

## Verification

`src/corr3_run.py` rebuilds both revised files **from their originals** rather than patching the
existing revisions, so every check CORR1 and CORR2 made holds end to end. It imports their edit
strings instead of restating them.

| check | paper 2 (v0.7 → v0.8) | paper 1 (v1.2 → v1.3) |
|---|---|---|
| string edits | 4 | 4 |
| each replaced string occurs exactly once | yes | yes |
| undoing the edits reproduces the original XML exactly | yes | yes |
| other zip members byte-identical | 38 of 38 | 53 of 53 |
| body paragraphs | 117 → 117 | 145 → 146 |
| paragraph-142 wording matches the findings blockquote | — | asserted |

Verified independently of that script, by a separate comparison: the two source files are
untouched, and their sha256 values are **identical to the ones CORR1 and CORR2 recorded in their
own committed diffs** (`v1.2 5dfeaf0e…` at commit `d8d892a`, `v0.7 74f5c339…` at commit
`46a05c1`). Only `word/document.xml` differs in each revision. For paper 1 an order-aware
comparison shows exactly three changes: one inserted paragraph at position 134, v1.2 paragraph 142
replaced, and v1.2 paragraph 144 replaced.

## Deviations and decisions not specified by the task (disclosed)

1. **Two previous tasks' results directories were written into, against CLAUDE.md.** The task
   instructed regenerating `results/corr1/paper1_diff.txt` and `results/corr2/paper2_diff.txt`,
   which belong to CORR1 and CORR2. CLAUDE.md says never to write into a previous task's results
   directory. The instruction was followed, for a reason that is stated here rather than assumed:
   a file whose entire content is "the paragraph-level diff of A against B" misdescribes the B it
   names once B changes, and both notices point readers at these files as the record of the
   revision. Each regenerated diff now carries a provenance header naming CORR3, stating what
   changed, and naming the commit (`d8d892a`, `46a05c1`) where the previous version can be read.
   Nothing is lost: git holds both earlier versions. This is the one place this task departed from
   a standing rule, and it was not put to the author first.
2. **Both revised docx files were rebuilt from their originals, not patched.** The task says to
   change v0.8 and to fold wording into v1.3. Rebuilding from v0.7 and v1.2 achieves the same
   content and makes the reversibility and byte-identity checks meaningful end to end; patching
   the existing revisions would have checked only the increment. The previous v0.8 and v1.3 are
   overwritten; both are git-ignored, so no tracked file was replaced.
3. **CORR1's in-document correction note was rewritten.** v1.3 carries a one-line italic
   "Correction (v1.3)" paragraph under the Supplementary Note 3 heading. CORR1's text named only
   "two sentences … in the last paragraph of this note", which would now be incomplete, so it was
   extended to mention the bullet. CORR1's XML wrapper is reused verbatim by substituting only its
   text, so the paragraph's formatting cannot drift. The task did not ask for this; leaving the
   note inaccurate was the alternative.
4. **`CORRECTIONS_PAPER2.md` is not titled a correction notice.** Paper 2 is unpublished. The
   filename follows the task; the title and first paragraph say plainly that this is a
   pre-publication revision. `CORRECTIONS_PAPER1.md`'s structure is otherwise followed.
5. **The appendix lists 22 clause-level verdicts for 21 statements**, because paragraph 31's
   sentence is one statement and two clauses. The same arithmetic slip is disclosed in
   `FINDINGS_CORR2CHK.md` deviation 8, where the pre-registration's three counts are each labelled
   "of 21" and sum to 22.
6. **`SUPPLEMENTARY_NOTE_COORDSPACE.md` was found to carry the superseded bullet and was left
   alone.** See step 3 above. The task did not mention the file; it was found by the Appendix A
   sweep.
7. **`FINDINGS_CORR3.md` and `PROGRESS_CORR3.md` were written although the task named neither.**
   CLAUDE.md's repository layout calls for them, and every other task here has them.
8. **Nothing was computed, so nothing here confirms any number it carries.** The numbers are
   copied. `FINDINGS_CORR2CHK.md` and `FINDINGS_FLOORCHK.md` are where they are checked, and
   `python -I src/corr2chk_run.py verify` re-derives the paper-2 readings from tracked tables.
9. **Where things were written.** `paper/paper2_passage_draft_v0.8.docx`,
   `paper/rejuvenation_readouts_preprint_v1.3.docx` (both git-ignored),
   `results/corr1/paper1_diff.txt`, `results/corr2/paper2_diff.txt`,
   `results/verify/WITHDRAWN_READINGS.md`, `CORRECTIONS_PAPER1.md`, `CORRECTIONS_PAPER2.md`,
   `README.md`, `src/corr3_run.py`, `FINDINGS_CORR3.md`, `PROGRESS_CORR3.md`. No `FINDINGS_*.md`
   of a previous task was edited, nothing under `results/corr2chk/` or `results/floorchk/` was
   touched, no original docx or PDF was modified, no PDF was produced, nothing was uploaded and
   nothing was pushed. The pre-existing uncommitted change to `CITATION.cff` and the untracked
   `results/paper2_figs_plain/` are not this task's work and are left out of its commit.

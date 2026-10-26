# Correction notice: paper 1, Supplementary Note 3

**Paper:** R. Kwon, *Rejuvenation readouts cannot distinguish younger from less old*, v1.2
(https://doi.org/10.5281/zenodo.23113817; PDF also at `paper/rejuvenation_readouts_preprint_v1.2.pdf`).
**Written:** 2026-10-26, task CORR1. **Revised:** 2026-10-26, task CORR3, which added the
paragraph-142 bullet below to the corrected sentences; CORR1's two sentences are unchanged.
**Basis:** `FINDINGS_FLOORCHK.md` and `results/floorchk/`, and for the bullet's numbers
`FINDINGS_CORR2CHK.md` and `results/coordspace/`. This notice computes nothing new. Every number
below is copied from those files.

Both corrected passages are in Supplementary Note 3. They are treated separately below because
they fail for different reasons: the sentences overstate what was measured, while the bullet
states something that could not have come out otherwise.

## What was stated

**(a)** Supplementary Note 3, last-but-one paragraph (PDF page 13; docx v1.2 body paragraph 144):

> In space C, for the GTEx old-to-young comparison, with the GTEx old centroid (60 to 79 years) as
> the origin and the GTEx young centroid (20 to 39 years) as the target, every one of the 1,836
> perturbations moved more than twice the distance between the two centroids (smallest ratio 2.10;
> per-perturbation values in results/south6/ in the repository), so from that origin no approach to
> that target could have been registered whatever its direction. This holds for that origin and
> target, which sit only 32 units apart, and is not a statement about the screen in general; it is
> the same limitation the main text discloses for the same-platform test.

**(b)** Supplementary Note 3, the fourth bullet of the space A diagnostics (PDF page 13; docx v1.2
body paragraph 142):

> •  In space C both perturbations move away from the young centroid (+106.6 and +113.8, q = 1.0),
> and in space B both end farther from it (+32.9 and +35.1). They have 84 and 82 cells, against a
> median of 106.

## What is wrong

### (a) the two sentences

The sentence compares each perturbation's displacement with the old-to-young gap. It does not
compare that ratio with a noise floor matched to the group size. In the Southard screen a
perturbation is measured on 13 to 341 cells. At those cell counts, non-targeting cells with no
effect already read more than twice the gap at every one of the screen's 212 group sizes. At the
screen's median group size (106 cells) their median is 3.76 gap lengths, against the
perturbations' own median of 3.82. Across group sizes the median runs from 2.06 (341 cells) to
11.9 (13 cells). The measured lengths are therefore mostly sampling noise.

FLOORCHK applied a pre-registered noise correction: sqrt(max(0, ‖d‖² − median‖d₀‖²)) / ‖gap‖, with
the floor drawn at each perturbation's own cell count. The correction passed its own validation
gate. Planted displacements of 1, 2 and 3 gap lengths were recovered within ±20% at every gating
size; the worst case was 0.960 for a planted 1 at 106 cells. After this correction:

- **18 of 1,836** perturbations are estimated to move more than twice the gap. About 5 (5.05)
  would be expected from noise alone. That chance baseline is a diagnostic, not part of the
  pre-registration.
- The median estimate is **0.61** gap lengths (0.607).

The 18 are smaller groups than the screen as a whole (median 54.5 cells against 106), so they are
not a shortlist. The pre-registered reading that fired is `CLAIM1-PARTIAL: 18 of 1,836`.

The closing comparison also fails. The same-platform length statement survives the same check (see
below), but the screen statement rests on sampling noise. They are not the same limitation.

### (b) the bullet

The bullet is offered as a diagnostic: it is one of four reasons given for not treating the two
space A labels as candidates. But neither of its two clauses distinguishes those two
perturbations from the rest of the screen, so it carries no information about them.

- **In space C no perturbation's change in distance could have been negative.** Every one of the
  1,836 measured displacements is longer than twice the 32-unit gap (smallest ratio 2.0966;
  `FINDINGS_FLOORCHK.md`, reproduction check R1), and `n_delta_negative` is 0 of 1,836
  (`results/coordspace/stage3_south6_summary.csv`). "Both perturbations move away" is therefore
  true of all 1,836, and the median non-targeting floor in that space is 2.9 gap lengths
  (`median_floor_over_gap` = 2.9415).
- **In space B none of the 1,836 moves closer either**, and both of these two sit just *inside*
  their own cell-count-matched floors: +32.886 against +33.344 and +35.053 against +35.705
  (`results/coordspace/stage3_south6_factors.csv`), as 1,088 of the 1,836 do
  (`n_below_floor` = 1088). So they move away by less than sampling noise does at their cell
  count, which is the opposite of what the clause implies.
- **The `q = 1.0` should not have been quoted as a gate.** It comes from the permutation null
  that SOUTH8 found anti-conservative in space C, flagging 0.09 to 0.14 of true-null draws at the
  0.05 level (`results/south8/d2_calibration.csv`, S1 rows; verdict in `d2_verdicts.csv`).

The only part of the bullet that distinguished these two perturbations was their cell count, 84
and 82 against a screen median of 106, and that part is kept.

## What still holds

It is still true that the test could not have registered an approach from that origin to that
target. Every measured displacement is longer than twice the gap, so no measured change in distance
could come out negative. The reason is the screen's depth, not the length of the displacements.

**The decision the bullet supported still stands**, on the other three diagnostics the paragraph
names: ten genes carry 94% of the squared old-to-young distance in space A, every endpoint there
has negative coordinates, and the space A ranking agrees poorly with the main analysis. The two
space A labels are still not treated as candidates. What changes is that the space C and space B
comparison is no longer offered as a reason, because it was true of the whole screen.

Corrected sentences (a), as in the revised source, v1.3:

> In space C, for the GTEx old-to-young comparison, with the GTEx old centroid (60 to 79 years) as
> the origin and the GTEx young centroid (20 to 39 years) as the target, every one of the 1,836
> measured displacements is longer than twice the distance between the two centroids (smallest
> ratio 2.10; per-perturbation values in results/south6/ in the repository), so from that origin
> the test could not register an approach to that target whatever the direction. This holds for
> that origin and target, which sit only 32 units apart, and is not a statement about the screen in
> general. It reflects depth, not effect size: non-targeting cells with no effect, sampled at the
> same group sizes, also read more than twice that distance (2.1 to 12.9 at the 95th percentile).
> After the noise is subtracted, 18 perturbations are estimated to move more than twice the
> distance (about 5 expected by chance), and the median estimate is 0.6 times it
> (results/floorchk/ in the repository). The space C null for the screen is therefore a null of an
> instrument that could not see approach at this depth, not evidence that the perturbations
> overshoot, and it is not the same limitation as the one the main text discloses for the
> same-platform test, whose length statement survives the same check.

Corrected bullet (b), as in the revised source, v1.3:

> •  Neither of the other spaces separates these two from the rest of the screen. In space C both
> end farther from the young centroid (+106.6 and +113.8), as do all 1,836 perturbations: there
> every displacement is longer than twice the 32-unit gap, so no change in distance could have
> come out negative, and the median non-targeting floor is 2.9 gap lengths. In space B both also
> end farther (+32.9 and +35.1), again with none of the 1,836 moving closer, and both sit just
> inside their own cell-count-matched floors (+33.3 and +35.7), as 1,088 of the 1,836 do. They
> have 84 and 82 cells, against a median of 106.

Two deliberate changes in that bullet, beyond adding the floor comparison: the `q = 1.0` is
dropped, for the reason given above, and the space B clause now carries its floors, because
"both end farther" was as uninformative there as in space C and the floors are the part that
says something.

## What was rechecked and holds

The main text states that the aged donor's partially reprogrammed cells (5,832 cells) travel more
than twice the line between the two donors, and in the MD genes 2.3 times that line. This
appears in Results (docx paragraph 22), Discussion (paragraph 51) and Limitations (paragraph 66).
FLOORCHK checked it with the same correction. It holds (`CLAIM2-HOLDS`).

| gene space | raw ratio | corrected, floor extrapolated to 5,832 cells | corrected, observed 2,470-cell floor |
|---|---|---|---|
| **all genes (23,485), primary** | 2.087 | **2.046** | **2.022** |
| MD genes (205), secondary | 2.300 | not used (fit failed) | 2.296 |
| age genes (3,089), secondary | 2.084 | not used (fit failed) | 2.047 |

**The all-genes margin is narrow.** Both corrected values are only 1 to 2% above 2.

- **One of the two floors is an extrapolation.** It is read from a 1/n fit at 5,832 cells,
  about 2.4 times (5,832 / 2,470 = 2.36) past the largest size drawn.
- **The floor comes from day-0 fibroblasts, not from partially reprogrammed cells.** If real noise
  were 1.33 times the observed floor, the all-genes reading would fall below 2. There is no
  no-effect pool of partially reprogrammed cells to check this.

In the MD and age genes the 1/n fit failed its 10% check. Only the observed floor was used there,
as pre-registered.

## Effect on the rest of the paper

**No other conclusion of the paper depends on the corrected sentences.** This was checked against
every passage of the v1.2 text that mentions the transcription-factor screen (docx paragraphs 58,
69, 138 to 144) or a length of more than twice a gap (paragraphs 22, 49, 51, 66):

- **The transcription-factor screen does not bear on the paper's question.** The paper says so
  itself (paragraph 69): the screen "does not bear on the question tested here and is not reported
  in this paper". The pre-registered space C result, that no perturbation is labelled as
  approaching, is a count. It stands, and the corrected text says what it means.
- **Every main-text length statement is about the reprogramming cells.** These are paragraphs 22,
  49, 51 and 66, and they rest on claim 2, which holds.
- **The decision not to treat the two space A labels as candidates is unaffected.** It rests on the
  diagnostics the paragraph names: ten genes and negative coordinates. One of the four
  diagnostics, the paragraph-142 bullet, is corrected above and no longer offered as a reason;
  the other three carry the decision.

(When this notice was first written, as task CORR1, the paragraph-142 bullet was recorded here as
a supporting detail that "loses its force but is not corrected here, because it was not among the
sentences this notice covers". Task CORR3 corrected it; it is now section (b) above. Nothing else
in this notice changed.)

## Where the check is recorded

In the public repository (https://github.com/k-dlrgus/rejuvenation-readouts):

For (a), the two sentences:

- Pre-registration: `results/floorchk/PREREG.flag`, written 2026-10-25T23:17:26+09:00 before any
  statistic.
- Findings: `FINDINGS_FLOORCHK.md`.
- Tables: `results/floorchk/`, in particular `claim1_floor.csv`, `claim1_units.csv`,
  `claim2_unit.csv` and `readings.txt`.
- Code: `src/floorchk_run.py`.

For (b), the bullet:

- Its wording was written as step 4 of `FINDINGS_CORR2CHK.md`, which computed nothing for it and
  quotes the source of every number.
- Tables: `results/coordspace/stage3_south6_factors.csv` and `stage3_south6_summary.csv` for the
  deltas, floors and counts; `results/south8/d2_calibration.csv` and `d2_verdicts.csv` for the
  permutation null's calibration; `FINDINGS_FLOORCHK.md` reproduction check R1 for the
  twice-the-gap statement.

For both: revised-source diff `results/corr1/paper1_diff.txt`; code `src/corr3_run.py`, which
rebuilt v1.3 from v1.2 and asserts that the bullet's wording is byte-identical to the blockquote
in `FINDINGS_CORR2CHK.md`.

---

## Appendix A: every place the failed claim appears

The search covered tracked and untracked text, the docx and PDF manuscripts, and the git-ignored
`bundle/`. It looked for: the "every one of the 1,836 perturbations" length statement,
"2.10" or "2.097" used for it, and the "same limitation the main text discloses" sentence.
Numeric matches in data tables that have nothing to do with the claim were dropped.

| # | file:line | tracked | what it says | action |
|---|---|---|---|---|
| 1 | `paper/rejuvenation_readouts_preprint_v1.2.pdf`, page 13 (pdftotext lines 728–734) | yes | The published sentences quoted above | Not edited. The PDF is the published record; a v1.3 PDF was not produced. |
| 2 | `paper/rejuvenation_readouts_preprint_v1.2.docx`, body paragraph 144 | no (`paper/` is git-ignored) | Source of #1 | Not edited. The revised copy is `paper/rejuvenation_readouts_preprint_v1.3.docx` (also git-ignored). |
| 3 | `paper/rejuvenation_readouts_preprint_v1.2_draft.docx`, same paragraph | no | Earlier draft of #2 | Not edited. Superseded by #2. |
| 4 | `FINDINGS_RATIO.md:35–36` | yes | "…exceeded twice the old-to-young gap" holds for SOUTH6 in the full-gene space S1 (1,836 of 1,836, minimum ratio 2.097) | Not edited (findings file). The raw ratio is correct as a measurement; FLOORCHK shows it is mostly sampling noise. |
| 5 | `FINDINGS_RATIO.md:89` | yes | The RATIO pre-registration's premise: "In the Southard screen, every perturbation's displacement exceeded twice the old-to-young gap, so no approach could register whatever its direction." | Not edited (verbatim pre-registration copy). |
| 6 | `results/ratio/PREREG.flag:33` | yes | Same sentence as #5 | Not edited (flag file; CLAUDE.md forbids writing into a previous task's results). |
| 7 | `bundle/RESULTS_SUMMARY.md:77` and `:467`, `bundle/_build/write_report.py:88` | no (`bundle/` is git-ignored, a review package) | "in S1 every factor's ‖d‖ exceeds twice ‖T − O‖ (min ratio 2.097)" | Not edited. Not public. |
| 8 | `FINDINGS_FLOORCHK.md`, `PROGRESS_FLOORCHK.md`, `results/floorchk/*`, `src/floorchk_run.py` | yes | The check itself, quoting the claim | Not edited. |

Rows 9 to 12 were added by CORR3 for the paragraph-142 bullet, (b). The search looked for
"106.6 and +113.8", "both perturbations move away", and `ZNF441` / `ZNF669` in prose.

| # | file:line | tracked | what it says | action |
|---|---|---|---|---|
| 9 | `paper/rejuvenation_readouts_preprint_v1.2.pdf`, page 13 | yes | The published bullet quoted above | Not edited. The PDF is the published record; no v1.3 PDF was produced. |
| 10 | `paper/rejuvenation_readouts_preprint_v1.2.docx`, body paragraph 142, and `_draft.docx` | no (`paper/` is git-ignored) | Source of #9 | Not edited. The revised copy is v1.3, paragraph 143. |
| 11 | **`SUPPLEMENTARY_NOTE_COORDSPACE.md`** | **yes** | The bullet verbatim, in a file whose first line reads "Coordinate-space robustness check: text for the preprint" | **Marked superseded in place; its wording is unedited.** A one-line marker was added directly above the bullet, dated 2026-10-26, naming paper 1 v1.3 paragraph 142 and pointing here to passage (b), and a one-line note under the file's title now states that it is drafted text and that the manuscript is authoritative where the two differ. The bullet itself is untouched, as CLAUDE.md requires of a superseded text. |
| 12 | `FINDINGS_COORDSPACE.md` (Stage 3 tables and the `toward_young` rows), `results/coordspace/stage3_south6_factors.csv` | yes | The underlying deltas +106.609 and +113.750 with their floors | Not edited. The numbers are correct as measurements and are the source this correction quotes; what failed was the inference drawn from them. |

No hit in `README.md`, `STATE.md`, `CITATION.cff`, `REPO_INVENTORY.md`, `FINDINGS_SOUTH6.md`,
`FINDINGS_COORDSPACE.md`, `SUPPLEMENTARY_NOTE_COORDSPACE.md`, `paper/PREPRINT_DRAFT.md`,
`paper_package/`, the earlier preprint drafts (`_3`, `_9`, `_10`, `_11`, `_12`), or the paper-2 drafts
and PDF. There is no blueprint or business document in the repository.

## Appendix B: preprint source

The preprint source is `paper/rejuvenation_readouts_preprint_v1.2.docx`. It is in the working tree
but **not tracked**, because `paper/` is git-ignored ("review packages, not public"). Only the v1.2
PDF is tracked.

The revised source `paper/rejuvenation_readouts_preprint_v1.3.docx` sits beside the original,
which is unchanged. It was first made by `src/corr1_run.py` and was rebuilt from v1.2 by
`src/corr3_run.py`, which imports CORR1's edit strings rather than restating them. It differs
from v1.2 in four ways:

1. A one-line, italic "Correction (v1.3)" paragraph directly under the Supplementary Note 3
   heading, naming both corrected passages.
2. The two sentences of (a), replaced by the corrected text.
3. The bullet of (b), replaced by the corrected text.
4. Nothing else. The script checks that each replaced string occurs exactly once, that undoing
   the four edits reproduces the v1.2 XML exactly, and that the other 53 files in the docx are
   byte-identical. v1.2's sha256 is recorded in the diff and is unchanged from CORR1's run.

The paragraph-level diff is `results/corr1/paper1_diff.txt`, which is tracked. The v1.3 docx
itself is git-ignored like the rest of `paper/`, so it is not committed. No v1.3 PDF was made, and
nothing was uploaded.

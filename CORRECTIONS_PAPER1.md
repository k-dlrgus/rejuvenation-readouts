# Correction notice: paper 1, Supplementary Note 3

**Paper:** R. Kwon, *Rejuvenation readouts cannot distinguish younger from less old*, v1.2
(https://doi.org/10.5281/zenodo.23113817; PDF also at `paper/rejuvenation_readouts_preprint_v1.2.pdf`).
**Written:** 2026-10-26, task CORR1. **Basis:** `FINDINGS_FLOORCHK.md` and `results/floorchk/`.
This notice computes nothing new. Every number below is copied from those files.

## What was stated

Supplementary Note 3, last-but-one paragraph (PDF page 13; docx v1.2 body paragraph 144):

> In space C, for the GTEx old-to-young comparison, with the GTEx old centroid (60 to 79 years) as
> the origin and the GTEx young centroid (20 to 39 years) as the target, every one of the 1,836
> perturbations moved more than twice the distance between the two centroids (smallest ratio 2.10;
> per-perturbation values in results/south6/ in the repository), so from that origin no approach to
> that target could have been registered whatever its direction. This holds for that origin and
> target, which sit only 32 units apart, and is not a statement about the screen in general; it is
> the same limitation the main text discloses for the same-platform test.

## What is wrong

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

## What still holds

It is still true that the test could not have registered an approach from that origin to that
target. Every measured displacement is longer than twice the gap, so no measured change in distance
could come out negative. The reason is the screen's depth, not the length of the displacements.

Corrected sentences (as in the revised source, v1.3):

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
  diagnostics the paragraph names: ten genes and negative coordinates.

One supporting detail loses its force but is not corrected here, because it was not among the
sentences this notice covers. The bullet in paragraph 142 reads: "In space C both perturbations
move away from the young centroid (+106.6 and +113.8, q = 1.0)". In space C no perturbation's
measured change in distance could have been negative, so that clause carries no information about
those two perturbations. The paper's decision does not cite it.

## Where the check is recorded

In the public repository (https://github.com/k-dlrgus/rejuvenation-readouts):

- Pre-registration: `results/floorchk/PREREG.flag`, written 2026-10-25T23:17:26+09:00 before any
  statistic.
- Findings: `FINDINGS_FLOORCHK.md`.
- Tables: `results/floorchk/`, in particular `claim1_floor.csv`, `claim1_units.csv`,
  `claim2_unit.csv` and `readings.txt`.
- Code: `src/floorchk_run.py`.
- Revised-source diff: `results/corr1/paper1_diff.txt`.

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

No hit in `README.md`, `STATE.md`, `CITATION.cff`, `REPO_INVENTORY.md`, `FINDINGS_SOUTH6.md`,
`FINDINGS_COORDSPACE.md`, `SUPPLEMENTARY_NOTE_COORDSPACE.md`, `paper/PREPRINT_DRAFT.md`,
`paper_package/`, the earlier preprint drafts (`_3`, `_9`, `_10`, `_11`, `_12`), or the paper-2 drafts
and PDF. There is no blueprint or business document in the repository.

## Appendix B: preprint source

The preprint source is `paper/rejuvenation_readouts_preprint_v1.2.docx`. It is in the working tree
but **not tracked**, because `paper/` is git-ignored ("review packages, not public"). Only the v1.2
PDF is tracked.

The revised source `paper/rejuvenation_readouts_preprint_v1.3.docx` was made by
`src/corr1_run.py`, beside the original, which is unchanged. It differs from v1.2 in three ways:

1. A one-line, italic "Correction (v1.3)" paragraph directly under the Supplementary Note 3
   heading.
2. The two sentences quoted above, replaced by the corrected text.
3. Nothing else. The script checks that undoing the three edits reproduces the v1.2 XML exactly,
   and that the other 53 files in the docx are byte-identical.

The paragraph-level diff is `results/corr1/paper1_diff.txt`, which is tracked. The v1.3 docx
itself is git-ignored like the rest of `paper/`, so it is not committed. No v1.3 PDF was made, and
nothing was uploaded.

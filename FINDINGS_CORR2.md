# FINDINGS_CORR2 — paper 2, v0.7 → v0.8: the two clauses CORR2CHK found unsupported

**Task:** CORR2. **Status:** complete. **Date:** 2026-10-26.
**Kind:** a writing task. **No statistic was computed.** Every number in the new text is copied
from `FINDINGS_CORR2CHK.md` (sections "Corrected wording (step 3)", "A1b — the Southard clause of
paragraph 31" and "C5 — the young-donor clause of paragraph 30") or from `results/corr2chk/`.
This task read no expression data and ran nothing on any dataset.

**Paper 2 is not public.** Its source is `paper/paper2_passage_draft_v0.7.docx`, whose text matches
`paper/two_definitions_of_young_preprint_v1.pdf`; neither is tracked, and no version of paper 2 has
been deposited or published anywhere. **This is therefore a pre-publication revision of a draft,
not a correction notice.** Nothing here corrects a public record, nothing here is addressed to
readers of a published paper, and no notice paragraph was inserted into the draft. That is the one
structural difference from CORR1, which revised an already-published paper (paper 1, v1.2 → v1.3)
and so carried a public notice, `CORRECTIONS_PAPER1.md`, and an italic in-document correction line.
There is no `CORRECTIONS_PAPER2.md` and there should not be one.

## What changed

`paper/paper2_passage_draft_v0.8.docx`, written beside v0.7 by `src/corr2_run.py`. Three string
edits to `word/document.xml`, in two body paragraphs; nothing else in the file differs.

| # | paragraph | CORR2CHK id | what was replaced |
|---|---|---|---|
| 1 | 30 | C5 | the sentence "None of the 80 … and all 80 move away from them." |
| 2 | 31 | A1b | "Separately, the four genes are also uncalled and move away from young donors in the Southard screen" → "… uncalled in the Southard screen" |
| 3 | 31 | A1b | the remainder of that sentence, after the superscript reference run |

Edits 2 and 3 are one clause-level replacement split in two, because a superscript run (the
Southard citation, "5") sits inside the sentence and is kept exactly where it was. CORR2CHK's
replacement text restates the sentence's closing activation-only clause, so the whole sentence
from "Separately," onward is replaced; nothing of the old sentence is left duplicated or orphaned.

### Paragraph 30 — before

> All four move away (Figure 2). Their change in distance to young donors is +11.4 (E2F3), +9.0
> (EZH2), +5.3 (STAT3) and +15.3 (ZFX), against thresholds of −1.0 to −1.3. … Across the screen,
> only 17 of 340 perturbations end any closer to young donors, and 5 are called, fewer than the 17
> expected by chance. **None of the 80 perturbations called closer to early-passage cells is called
> closer to young donors, and all 80 move away from them.** Figure 2 plots each change relative to
> its own threshold …

### Paragraph 30 — after

> **None of the 80 perturbations called closer to early-passage cells is called closer to young
> donors, and all 80 have a positive change in distance to them. Against no-effect samples of the
> same size, 56 of the 80 move away further than the 95th percentile, 23 lie inside the no-effect
> spread and 1 move away less far than the median no-effect sample, so the departure is established
> for 56 of them and not for the rest.**

The rest of paragraph 30 is unchanged, including its first sentence about the four validated
factors and its closing sentence about the thresholds depending on group size.

### Paragraph 31 — before

> Against this target, distance adds little to direction. Each effect starts from the old-donor
> average, so any effect that points away from the donor-age direction, as all four do, must
> increase the distance to young donors. Distance would add information only for an effect that
> pointed toward young donors, where it could tell an approach from an overshoot. **Separately, the
> four genes are also uncalled and move away from young donors in the Southard screen⁵, a different
> neonatal line; that screen used activation only, so for STAT3 and ZFX it tests the opposite
> perturbation.**

### Paragraph 31 — after

> **Separately, the four genes are also uncalled in the Southard screen⁵, a different neonatal
> line, where the test detects a quarter-gap approach toward young donors at their cell counts with
> power 0.925 to 1.000. Their changes in distance to young donors are positive, but so are those of
> non-targeting cells with no effect at the same cell counts: the no-effect median runs from +1.3 to
> +2.7 against the four factors' +0.5 to +2.4. Three of the four move away less far than the median
> no-effect sample of their own size, and none exceeds its 95th percentile, so that screen shows no
> approach and does not show a departure either. That screen used activation only, so for STAT3 and
> ZFX it tests the opposite perturbation.**

The paragraph's first three sentences — distance adding little to direction, the geometric
consequence of a negative cosine from the old-donor start, and distance adding information only for
an effect pointing toward young donors — are unchanged. The activation-only clause survives, as the
last sentence of the replacement.

## Why each clause failed, with the numbers copied from CORR2CHK

Both clauses fired CORR2CHK's pre-registered key `PAPER2-DEPTH`. Both are of the form "these units
move away from young donors", and CORR2CHK's rule fails such a claim if even one unit is not
distinguishable from no effect, because the published wording is universal. They fail very
differently, and the revision says so differently.

**Paragraph 31 (A1b) — no support at all.** Of the four factors in the Southard screen, in paper 2's
20-component space, from the GTEx old centroid: **0 of 4 exceed the 95th percentile of no-effect
draws at their own cell count (about 0.2 expected by chance), 1 is inside the no-effect spread, and
3 move away less far than the median no-effect draw (about 2.0 expected by chance).** At those
depths a no-effect draw's change in distance is itself positive, median **+1.3 to +2.7**, against
the four factors' **+0.5 to +2.4** — so the positive sign the old clause reported is the baseline's
own sign. The "uncalled" half of the sentence stands and is better supported than the old sentence
claimed: none of the four falls below its own floor, and the test detects a quarter-gap approach at
those cell counts with power **0.925 to 1.000**. That power is why the new text can say the screen
"shows no approach" rather than merely failing to see one. Per-factor values:
`results/corr2chk/a1b_units.csv`.

**Paragraph 30 (C5) — the aggregate holds; only the universal quantifier failed.** Of the 80
perturbations called closer to early-passage cells, against no-effect samples at each one's own
cell count: **56 exceed the 95th percentile (about 4.0 expected by chance), 23 lie inside the
no-effect spread, and 1 falls below its own median (about 40.0 expected by chance).** The sign is
positive for all 80, which is why the new text keeps "all 80 have a positive change in distance",
and the departure is established against a size-matched baseline for 56. The 24 that are not
individually distinguishable are the smaller groups (CORR2CHK: cell counts 32 to 187, median 80;
the no-effect baseline rises as cell count falls), so the failure of "all" is a depth effect, not a
sign that those perturbations behave differently. Per-unit values: `results/corr2chk/c5_units.csv`;
per-size baselines, two seeds: `results/corr2chk/c5_baseline.csv`.

**What neither change does.** Neither says any factor moves *toward* young donors; neither is
evidence of rejuvenation by either definition. Both call-based statements are untouched: paragraph
31's "uncalled", and paragraph 30's "None of the 80 … is called closer to young donors". Paper 2's
conclusion rests on the calls — 80 of 340 called closer to early-passage cells against 5 of 340
called closer to young donors, 17 expected by chance in each case — and CORR2CHK found both counts
gated by size-matched thresholds and unaffected. What the revision removes is a stronger directional
claim laid on top of the calls.

**Honesty obligations carried over from CORR2CHK, which hold for the revised text too.** The chance
baselines are in the paragraph above beside every count; the 95th percentile is exceeded by 5% of
no-effect draws by construction and the median by 50%. The nulls are qualified by power, not left
bare (0.925 to 1.000 for A1b; AXANG's D3a gate, QUALIFIED at p25 = 0.945, for C5). Both baselines
rest on 200 draws, so verdicts for units close to the boundary are not stable at that precision,
and CORR2CHK reports the seed-to-seed movement. Both statements are circular in a way CORR2CHK
discloses: the same control cells define zero and supply every no-effect draw, and for C5 the 80
units are selected by a statistic built on those cells. Neither screen has an aged donor — both are
single neonatal fibroblast lines — so nothing here transfers to aged human cells.

## Where the check is recorded

- Pre-registration: `results/corr2chk/PREREG.flag`, written 2026-10-26T11:21:38+09:00 before any
  statistic of that task existed.
- Findings: `FINDINGS_CORR2CHK.md`; progress: `PROGRESS_CORR2CHK.md`.
- Tables: `results/corr2chk/`, in particular `a1b_units.csv`, `c5_units.csv`, `c5_baseline.csv`,
  `reproduction.csv` and `readings.txt`.
- Code: `src/corr2chk_run.py` (modes `run`, `verify`) and `src/corr2chk_findings.py`.
  `python -I src/corr2chk_run.py verify` re-derives every reading from the tracked CSVs alone.
- This revision: `src/corr2_run.py`; paragraph-level diff `results/corr2/paper2_diff.txt`.

## Verification of the revision

`src/corr2_run.py` performs, and asserts, all three checks CORR1's script performs:

1. **Each old string occurs exactly once** in v0.7's `word/document.xml` before its replacement.
2. **Every other zip member is byte-identical**: 38 of 38 non-`document.xml` members, and the
   written `document.xml` equals the intended string exactly.
3. **The edits are reversible**: undoing the three replacements in reverse order reproduces v0.7's
   `word/document.xml` byte for byte.

It also asserts that the body-paragraph count is unchanged (117 in both), since this revision
inserts no paragraph, and it refuses to overwrite an existing v0.8.

| file | sha256 |
|---|---|
| `paper/paper2_passage_draft_v0.7.docx` | `74f5c3392a652cd7e4034d65ff9174eefb8a2cf79ba21021168dec29e9deaabd` |
| `paper/paper2_passage_draft_v0.8.docx` | `bc751df1937a344872b3d35e606980c12dfbbe34646d0707ec6b2b3f45f9ee7b` |

v0.7's hash was taken before and after the run and did not change; the earlier drafts v0.1 to v0.6
and their PDFs were not opened for writing and their modification times are unchanged.

## Consistency sweep: every place either clause's content appears

The search covered tracked and untracked text: the whole repository (excluding `.git/` and
`data/raw/`), every paper-2 docx draft and PDF, paper 1's docx and PDF, the git-ignored `bundle/`
review package, `paper_package/`, `notebooks/`, and the figure-provenance directories
`results/paper2_figs*/`. It looked for: "all 80 move away", "80 move away from", "move away from
them", "uncalled and move away", "move away from young donors", "move them away", "80 perturbations
called closer", "uncalled", "move away from GTEx young". Matches about paper 1's reprogramming
figures, which have nothing to do with these two clauses, were dropped.

| # | file:where | tracked | what it says | action |
|---|---|---|---|---|
| 1 | `paper/paper2_passage_draft_v0.7.docx`, body paragraph 30 | no (`paper/` is git-ignored) | C5: "and all 80 move away from them" | **Not edited.** The revised copy is `paper/paper2_passage_draft_v0.8.docx` (also git-ignored). v0.7 is byte-identical to before this task. |
| 2 | `paper/paper2_passage_draft_v0.7.docx`, body paragraph 31 | no | A1b: "uncalled and move away from young donors in the Southard screen" | **Not edited**, as #1. |
| 3 | `paper/paper2_passage_draft_v0.8.docx`, paragraphs 30 and 31 | no | the corrected clauses | **Written by this task.** The only file of paper 2 that carries the new wording. |
| 4 | `paper/two_definitions_of_young_preprint_v1.pdf` | no | the PDF of v0.7; carries both clauses verbatim | **Not edited, and no v0.8 PDF was produced.** Paper 2 is unpublished, so this PDF is a draft rendering, not a record of anything; typesetting the revised draft is the author's step, outside this task. |
| 5 | `paper/paper2_passage_draft_v0.6.docx` paragraphs 34, 35, and `v0.6.pdf` | no | both clauses, in near-final wording ("it is the opposite perturbation" for "it tests the opposite perturbation") | **Not edited. Superseded by v0.7.** |
| 6 | `paper/paper2_passage_draft_v0.2` to `v0.5` `.docx` (paragraphs 51–54) and their PDFs | no | both clauses in earlier wording: "called against early-passage cells … and all 80 move away from them (descriptive, from the stored per-perturbation tables)"; "the four gene names are also uncalled and move away from young donors in the Southard screen" | **Not edited. Superseded.** |
| 7 | `paper/paper2_passage_draft_v0.1.docx` paragraph 51 and `v0.1.pdf` | no | C5 only, as "The pattern holds beyond the four: none of the 80 … and all 80 move away from them (descriptive …)". No Southard clause yet. | **Not edited. Superseded.** |
| 8 | `paper/paper2_passage_draft_v0.7.docx` paragraph 7 (SUMMARY) | no | "All four factors the authors validated … move cells closer to early-passage cells, and all four move them away from young donors." | **Not edited, deliberately.** This is the *Sengstack* young-donor test of the four factors (CORR2CHK id **C4**, paragraphs 7, 29, 30), which CORR2CHK ruled `PAPER2-NOT-AFFECTED`: its thresholds are drawn at each perturbation's own cell count, from the same space and starting profile (`src/axang_run.py:516`, floor built at `:476`; plant gate `:486-506`, power 0.945–0.990 at a quarter of the gap). The four factors' +11.4, +9.0, +5.3 and +15.3 are called against thresholds of −1.0 to −1.3. Nothing in A1b or C5 touches it; "correcting" it would be wrong. |
| 9 | `paper/paper2_passage_draft_v0.7.docx` paragraph 28 (section heading) and the first two sentences of paragraph 30 | no | "Against young donors, all four move away"; "+11.4 (E2F3), +9.0 (EZH2), +5.3 (STAT3) and +15.3 (ZFX), against thresholds of −1.0 to −1.3" | **Not edited.** Same statement as #8 (C4), gated and not affected. The heading and the four factors' own numbers stand as written. |
| 10 | `paper/paper2_passage_draft_v0.7.docx` paragraph 34 (Figure 2 caption) | no | "The 80 called closer to early-passage cells, including the four validated factors, all lie above the horizontal line." | **Not edited.** On both axes zero *is* the threshold, so "above the horizontal line" says these 80 are not called closer to young donors — the call-based half of C5, which holds. The caption makes no claim that the departure is established for each of the 80, so it stays consistent with the revised paragraph 30. |
| 11 | `paper/paper2_passage_draft_v0.7.docx` paragraph 57 (DISCUSSION) | no | "The same effects carry cells away from young human donors." | **Not edited.** "The same effects" are the four validated factors of #8 (C4), gated, not affected. |
| 12 | `paper/paper2_passage_draft_v0.7.docx` paragraph 70 (LIMITATIONS) | no | "Against this target, a departure follows from the sign of each effect's cosine." | **Not edited.** A geometric statement (CORR2CHK id C12, no floor needed) that the revision strengthens rather than contradicts: it is exactly why a positive change in distance is not by itself evidence of an effect. |
| 13 | rest of paper 2's DISCUSSION and LIMITATIONS (paragraphs 56, 58–62, 64–72) | no | checked for any repeat of either clause's content | **No hit.** Nothing else in the Discussion asserts a size-unmatched departure, for the four factors or for the 80. |
| 14 | `results/paper2_figs_fix/PLAN.flag:59`, `README.md:71`, `run2/README.md:95-96`, `run3/README.md:98-99,131` | yes | "Of the 80 perturbations called against early-passage cells (origin i), 0 are called against GTEx young and all 80 move away from GTEx young (delta_gtex_young > 0)" | **Not edited**, for two reasons. It is a previous task's results directory, which CLAUDE.md forbids writing into; and the wording already defines itself as the *sign* of `delta_gtex_young`, which is the part that survives and which v0.8 states as "all 80 have a positive change in distance to them". These are figure-provenance records of what was plotted, not prose claims. |
| 15 | `results/paper2_figs/paper2_numbers.md`, `results/paper2_figs_fix/run3/paper2_numbers.md` | yes | searched for both clauses | **Not edited: no hit.** Neither number file carries either clause's wording. |
| 16 | `results/paper2_figs_plain/`, `results/paper2_figs_v06/` (figures, PDFs and `figs_plain.py`) | `paper2_figs_plain/` is untracked; `paper2_figs_v06/` tracked | searched for both clauses | **Not edited: no hit.** `results/paper2_figs_plain/` is a pre-existing untracked working-tree directory that is not this task's work and is excluded from its commit. |
| 17 | `FINDINGS_AXANG.md:869` | yes | "All four are positive-delta and uncalled in Southard too, on a different [line]" — A1b's content, as AXANG recorded it | **Not edited.** CLAUDE.md: a superseded reading is marked superseded in the new file, never edited in the old. CORR2CHK supersedes the directional half of this sentence; the "uncalled" half stands. This file is named here so the record points both ways. |
| 18 | `FINDINGS_CORR2CHK.md`, `PROGRESS_CORR2CHK.md`, `results/corr2chk/PREREG.flag`, `src/corr2chk_run.py`, `src/corr2chk_findings.py` | yes | quote both clauses as the subject of the check, and hold the replacement wording | **Not edited.** The check itself; the flag is a verbatim pre-registration. |
| 19 | `bundle/` (git-ignored paper-1 review package), `paper_package/`, `paper/PREPRINT_DRAFT.md` | no / yes / no | searched for both clauses | **Not edited: no hit on either clause.** The only "move away" hits are paper 1's Figure 1 caption and a paper-1 section heading about partially reprogrammed fibroblasts, which are a different claim about a different dataset. `bundle/` contains no paper-2 Results text. |
| 20 | `README.md`, `STATE.md`, `REPO_INVENTORY.md`, `CITATION.cff`, `notebooks/`, `NUMBERS_NEWSTORY.md` and every other top-level `.md` | yes | searched for both clauses | **Not edited: no hit.** The two `README.md` hits are paper 1's Figure 1 and the repository's one-line summary of the age-score argument. |
| 21 | paper 1: `paper/rejuvenation_readouts_preprint_v1.2.pdf` / `.docx` paragraph 142 bullet, and `v1.3.docx` | PDF tracked; docx not | "In space C both perturbations move away from the young centroid (+106.6 and +113.8, q = 1.0)" | **Not edited; out of scope.** A different claim, in paper 1, about two labels in space C and space B. CORR2CHK's step 4 wrote replacement wording for it, for folding into a paper-1 revision; that is a separate task, and this one touched no paper-1 file. |

**Nothing became inconsistent.** No passage of paper 2 outside paragraphs 30 and 31 asserts either
failed clause's content. Every other "move away" statement in paper 2 is one of three things: the
four factors' gated Sengstack result (#8, #9, #11), a call-based statement (#10), or the geometric
identity the revision relies on (#12).

## Deviations and decisions not specified by the task (disclosed)

1. **No correction notice, and no in-document notice paragraph.** CORR1 wrote both, because paper 1
   is published. Paper 2 is not public, so there is no record to correct and nobody to notify; the
   task brief says the same. The consequence is that v0.8 differs from v0.7 only in the two
   paragraphs' text, with no added paragraph, and that `src/corr2_run.py` asserts the paragraph
   count is unchanged — a check CORR1's script could not make, since it inserted a paragraph.
2. **CORR2CHK's replacement text is used verbatim, with one typographic substitution.** The
   blockquote in `FINDINGS_CORR2CHK.md` writes "the four factors' +0.5 to +2.4" with a straight
   apostrophe; the docx uses the typographic apostrophe (U+2019) everywhere, including in the
   adjacent "screen's" of the same paragraph, so U+2019 was used. No word was changed.
3. **One grammatical slip in the pre-registered text was kept.** The paragraph-30 replacement reads
   "1 move away less far than the median no-effect sample", where "moves" would agree. The same
   form appears in CORR2CHK's own pre-registered conclusion wording ("1 move away less far than the
   median no-effect draw"). It is the pre-registered output of the previous task, so it was copied
   as written rather than improved silently. This note is the disclosure; changing it is the
   author's call, not this session's.
4. **Paragraph 31's replacement is applied as two string edits, not one.** A superscript reference
   run sits mid-sentence, so replacing the sentence as a single string would have deleted the
   Southard citation or moved it. The two edits replace the text either side of that run and leave
   the run itself untouched. The resulting paragraph text is exactly CORR2CHK's wording with the
   citation marker where it was, after "the Southard screen".
5. **The scope of each replacement was a judgement call.** CORR2CHK names the failing clauses as
   "and all 80 move away from them" and "and move away from young donors in the Southard screen",
   but its replacement text spans whole sentences in both cases. Each replacement therefore covers
   the sentence the clause sits in — for paragraph 31, the sentence from "Separately," to the end —
   so that the paragraph reads as written in CORR2CHK and no fragment is duplicated. Every other
   sentence of both paragraphs is byte-identical to v0.7.
6. **Interpreter note.** `lxml` is not in `requirements.txt` and is absent from the repository's
   `.venv`; `src/corr1_run.py` imports it too. Both scripts were run under the machine's Anaconda
   interpreter, which has it. Nothing else about the environment matters: the script reads one docx
   and writes one docx and one text file.
7. **Where things were written.** `paper/paper2_passage_draft_v0.8.docx`, `results/corr2/`,
   `src/corr2_run.py`, `FINDINGS_CORR2.md`, `PROGRESS_CORR2.md`. Nothing else was created or
   modified. No previous task's results directory was written into, no existing `FINDINGS_*.md`,
   `PROGRESS_*.md`, `CORRECTIONS_PAPER1.md`, `SCOPING_*.md`, `README.md`, `STATE.md` or
   `CITATION.cff` was touched, and no earlier draft or PDF of either paper was opened for writing.
   The pre-existing uncommitted change to `CITATION.cff` and the untracked
   `results/paper2_figs_plain/` are not this task's work and are left out of its commit.
8. **What this task did not do.** It computed nothing, so it cannot and does not confirm CORR2CHK's
   numbers; it copied them. It produced no PDF, uploaded nothing, and pushed nothing. It did not
   revisit any of CORR2CHK's 17 `PAPER2-NOT-AFFECTED` verdicts or its three unfloored length
   ratios, and it did not act on CORR2CHK's step-4 wording for paper 1.

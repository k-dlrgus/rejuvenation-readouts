# Correction notice: paper 2, paragraphs 30 and 31

**Paper:** R. Kwon, *Two definitions of young give opposite answers in a fibroblast rejuvenation
screen*, v1 (https://doi.org/10.5281/zenodo.23591623; PDF also at
`paper/two_definitions_of_young_preprint_v1.pdf`). Not to be confused with
https://doi.org/10.5281/zenodo.23591315, which is the paper's code and results archive, not the
paper.

**Written:** 2026-10-26, task CORR3. **Basis:** `FINDINGS_CORR2CHK.md` and `results/corr2chk/`.
This notice computes nothing. Every number below is copied from those files.

**The corrected text is the paper's v2.** Its source is
`paper/paper2_passage_draft_v0.8.docx`, beside the v1 source
`paper/paper2_passage_draft_v0.7.docx`, which is unchanged. The paragraph-level diff is
`results/corr2/paper2_diff.txt`.

**v2 carries no in-document correction note.** A new version carries the corrected text only;
the explanation of the change belongs in the Zenodo version notes and in this notice. No
paragraph is inserted or removed, so v2's paragraph numbering is v1's: the two corrected
paragraphs are 30 and 31 in both.

That the v1 source is `paper2_passage_draft_v0.7.docx` was checked before anything else, by
extracting the posted PDF's text and comparing it with that file's body text: all 99 of its
non-empty body paragraphs appear in the PDF, verbatim and in the same order, covering 98.2% of
its text, and all 52 non-empty Table 1 cells are present. No paragraph differs. The only
normalisation needed was to drop hyphens, because `pdftotext` loses a hyphen that falls at a line
break, so a change that was only a hyphen would not have been detected.

## What was stated

**(a)** Paragraph 31, the last sentence (PDF page 3):

> Separately, the four genes are also uncalled and move away from young donors in the Southard
> screen5, a different neonatal line; that screen used activation only, so for STAT3 and ZFX it
> tests the opposite perturbation.

**(b)** Paragraph 30, the fourth sentence (PDF page 3):

> None of the 80 perturbations called closer to early-passage cells is called closer to young
> donors, and all 80 move away from them.

Both are of the form "these units move away from young donors". Both were measured as the sign of
a change in distance against zero, with no comparison to what a group of the same size reads with
no effect at all. That comparison is what CORR2CHK added, using no-effect samples drawn at each
unit's own cell count, in the same coordinate space and from the same starting profile.

## What is wrong

The two clauses fail very differently, and they are kept apart here for that reason.

### (a) paragraph 31: the clause has no support

The sentence joins two claims with "and". **"Uncalled" holds** and is dealt with under "What
still holds". **"Move away from young donors" does not.**

At the four factors' Southard cell counts, a group of non-targeting cells with no effect at all
already reads a **positive** change in distance to the GTEx young centroid. Its median runs from
+1.3 to +2.7. The four factors read +0.5 to +2.4. So the positive sign the clause reports is the
sign the no-effect baseline already has at that depth.

| gene | cells | change in distance | floor (5th pct) | no-effect median | no-effect 95th pct | power at 25% of gap | verdict |
|---|---|---|---|---|---|---|---|
| E2F3 | 96 | +1.562 | −0.107 | **+2.520** | +5.643 | 0.975 | moves away less than the median no-effect sample |
| EZH2 | 140 | +0.481 | −0.694 | **+1.294** | +3.967 | 1.000 | moves away less than the median no-effect sample |
| STAT3 | 83 | +2.320 | +0.468 | **+2.695** | +5.932 | 0.925 | moves away less than the median no-effect sample |
| ZFX | 129 | +2.441 | −0.995 | **+1.514** | +4.180 | 0.995 | inside the no-effect spread |

- **0 of 4** exceed the 95th percentile of no-effect samples at their own cell count. Under a true
  null of no effect anywhere, about **0.2** of four would.
- **3 of 4** fall below the median no-effect sample of their own size. Under the same null, about
  **2.0** would.

What was observed is what no effect looks like. The pre-registered reading that fired is
`PAPER2-DEPTH: A1b`.

This is the same defect that SOUTH7's criterion amendment was written for, in the opposite
direction: there a factor with a positive change in distance was called *approaching*, because it
moved away less far than noise usually does.

**The clause would also have carried no information in the full gene space**, which the published
sentence does not name. There no factor's change in distance could have been negative at all:
every one of the 1,836 Southard displacements is longer than twice the gap. The no-effect medians
at these four depths are +110.5, +100.6, +82.2 and +78.0 against the four factors' +112.7, +98.2,
+83.8 and +79.2 — again all inside the spread. So the clause fails in the 20-component space and
is vacuous in the full gene space.

### (b) paragraph 30: the direction holds for the set; only the word "all" fails

Here the aggregate claim is strongly supported, and that should be said first.

- **56 of the 80** exceed the 95th percentile of no-effect samples at their own cell count,
  against about **4.0** expected by chance. That is a large, real departure.
- **1 of the 80** falls below its own median, against about **40.0** expected by chance.
- **23** lie between the two, inside the no-effect spread.

So 24 of the 80 are not individually distinguishable from no effect, and the universal "all 80
move away from them" claims more than the measurement supports. The one that moves away by less
than its own median no-effect sample is `CRA_TFCP2`, at 120 cells: +0.927 against a no-effect
median of +1.192.

**The units that fail are the shallow ones**, which is why this is a depth effect and not a
statement about those perturbations:

| verdict | units | cells: min / median / max | median change in distance |
|---|---|---|---|
| beyond the no-effect 95th percentile | 56 | 32 / 82 / 187 | +9.61 |
| inside the no-effect spread | 23 | 33 / 64 / 140 | +5.44 |
| below the no-effect median | 1 | 120 / 120 / 120 | +0.93 |

The no-effect baseline rises as cell count falls, so the units it swallows are the small ones. The
80 span 32 to 187 cells, median 80. The pre-registered reading that fired is `PAPER2-DEPTH: C5`.

## What still holds

**The paper's conclusion is unaffected.** That conclusion is that the two definitions of young
give opposite answers in the same dataset, and it rests on the **calls**, not on either corrected
clause:

- **80 of 340** perturbations are called closer to early-passage cells, against 17 expected by
  chance.
- **5 of 340** are called closer to young donors, against the same 17 expected by chance.
- Both counts are gated by thresholds drawn from control samples at each perturbation's own cell
  count, in the same space and from the same starting profile. Neither is affected.

**"Uncalled" in paragraph 31 stands, and is better supported than the sentence claimed.** None of
the four falls below its own floor, and the test detects a quarter-gap approach toward young
donors at those cell counts with power **0.925 to 1.000**. It is a qualified null, not a blind
one. The floor is a 5th percentile, so 95% of true nulls are "uncalled" by construction; it is
the power, not the count, that makes the null informative.

**The Summary's "all four move them away from young donors" is not affected, and a reader will
ask why.** That sentence is about the *Sengstack* young-donor test, not the Southard screen —
CORR2CHK's statement **C4**, which it ruled `PAPER2-NOT-AFFECTED`. There the four factors read
+11.4, +9.0, +5.3 and +15.3 against thresholds of −1.0 to −1.3, while the no-effect median at
those depths is only +0.81 to +1.50, and the test's gate fired QUALIFIED with power 0.945 at a
quarter of the donor-age gap. The four are three and a half to nineteen times the no-effect median
there. The corrected clause (a) is about a different screen, at a third of the effect size; the
corrected clause (b) is about the other 76 of the 80, not about the four. Nothing in either
correction touches the Summary.

**Of the 21 paper-2 statements inventoried, 17 need no correction** — their thresholds were
already matched to group size. Three length ratios have no size-matched floor anywhere and were
recorded rather than checked. The full list is in the appendix.

**What these corrections do not say.** Neither says any factor moves *toward* young donors.
Nothing here is evidence of rejuvenation by either definition. Clause (a) says the Southard
measurement cannot tell a departure from no effect at those four depths; clause (b) says the
Sengstack measurement can tell one for most of the 80 but not for the smallest groups. No factor
named here is a candidate for anything.

## Corrected wording

As in v2, whose source is `paper/paper2_passage_draft_v0.8.docx`.

**(a)** Paragraph 31, last sentence:

> Separately, the four genes are also uncalled in the Southard screen5, a different neonatal line,
> where the test detects a quarter-gap approach toward young donors at their cell counts with
> power 0.925 to 1.000. Their changes in distance to young donors are positive, but so are those
> of non-targeting cells with no effect at the same cell counts: the no-effect median runs from
> +1.3 to +2.7 against the four factors' +0.5 to +2.4. Three of the four move away less far than
> the median no-effect sample of their own size, and none exceeds its 95th percentile, so that
> screen shows no approach and does not show a departure either. That screen used activation only,
> so for STAT3 and ZFX it tests the opposite perturbation.

**(b)** Paragraph 30, fourth sentence:

> None of the 80 perturbations called closer to early-passage cells is called closer to young
> donors, and all 80 have a positive change in distance to them. Against no-effect samples of the
> same size, 56 of the 80 move away further than the 95th percentile, 23 lie inside the no-effect
> spread and 1 moves away less far than the median no-effect sample, so the departure is
> established for 56 of them and not for the rest.

One wording note: as first drafted by task CORR2, (b) read "1 move away less far", which is the
form carried by CORR2CHK's own pre-registered conclusion text. CORR2 copied it rather than improve
pre-registered wording silently, and disclosed it. Task CORR3 corrected the verb in the v2 source
only. `results/corr2chk/readings.txt` and `FINDINGS_CORR2CHK.md` keep the original form, because
they are the record.

## How strong the pre-registration behind this is, and where it is weaker than paper 1's

This has to be stated plainly, because it is weaker than the check behind `CORRECTIONS_PAPER1.md`
and a reader cannot tell from the outside.

- **The flag was written after the inventory, and the inventory already implied both answers.**
  CORR2CHK's step 1 was an inventory of published statements, and the author asked for it to be
  reported before the flag was written. That inventory read the stored baselines and the stored
  changes in distance. By the time `results/corr2chk/PREREG.flag` was written, at
  2026-10-26T11:21:38+09:00, the session already expected clause (a) to fail and expected at
  least one of (b)'s 80 units to fall below its baseline. Part 0 of the flag records exactly what
  had been read, with the numbers.
- **FLOORCHK, behind paper 1's correction, was stronger**: its flag was written before its
  statistic existed and its outcome was not implied by anything already read, beyond one hand
  calculation it discloses.
- **What still constrains this result.** The exact wording of all four possible conclusions was
  fixed by the author in the task instruction **before** the inventory ran and before any of those
  stored values were read. So the session could not choose, after seeing the numbers, what a
  failure would be called or how it would be phrased. The decision rule — a claim of the form
  "these move away" fails if even one unit is not distinguishable from no effect — was fixed in
  the flag before the comparison was tabulated, and it is a rule about the published wording, not
  about the data. And both corrected clauses are **universal** claims, so the finding does not
  depend on where a threshold sits: for (a) nothing supports it at all, and for (b) a single
  exception is enough and there are 24.
- **What would have made it stronger.** Writing the flag before the inventory, and reporting the
  inventory afterwards. The author's instruction asked for the inventory first so that only
  statements that needed checking would be checked; that was the right call for scope and the
  wrong one for blinding, and the two cannot both be had in one pass.

## Where the check is recorded

In the public repository (https://github.com/k-dlrgus/rejuvenation-readouts):

- Pre-registration: `results/corr2chk/PREREG.flag`, written **2026-10-26T11:21:38+09:00**, before
  any statistic of that task existed. Its Part 0 discloses what was already known; Parts 1 and 2
  reproduce the task instruction and the author's rulings verbatim.
- Findings: `FINDINGS_CORR2CHK.md`, with the flag reproduced verbatim at its top.
- Tables: `results/corr2chk/`, in particular `a1b_units.csv`, `c5_units.csv`, `c5_baseline.csv`
  and `readings.txt`. Reproduction of every stored value used: `reproduction.csv`, 41 gating
  checks and 0 mismatches.
- Code: `src/corr2chk_run.py` (modes `run`, `verify`) and `src/corr2chk_findings.py`.
  `python -I src/corr2chk_run.py verify` re-derives every reading from the tracked tables alone.
- Revised-source diff: `results/corr2/paper2_diff.txt`. Code for the revision: `src/corr2_run.py`
  and `src/corr3_run.py`.

The baseline for clause (a) was not drawn by that task: it already existed, keyed to each
factor's own cell count, in `results/south7/stage_a1_depth_nulls.csv`, built from 200 draws with
the reference being the non-targeting pool minus the draw. For clause (b) AXANG's floor loop was
re-run to recover the medians it had computed and discarded, at AXANG's own seed, and it
reproduced every stored 5th percentile to 2.2e-16.

---

## Appendix: all 21 inventoried statements, with verdicts

Paragraph numbers are docx body-paragraph indices, the same in v1 and v2. Splitting paragraph 31's sentence into
its two clauses gives 22 clause-level verdicts from 21 statements. **k20** is the paper's
20-component space; **full** is the 23,485-gene space.

| id | ¶ | what it rests on | space | verdict |
|---|---|---|---|---|
| **A1a** | 31 | "the four genes are also uncalled" (Southard) | k20 | NOT-AFFECTED — threshold at each factor's own cell count; power 0.925–1.000 |
| **A1b** | 31 | "and move away from young donors in the Southard screen" | k20 | **PAPER2-DEPTH — corrected** |
| A2 | 18, 77 | the 20-component count carried over from the Southard screen; quoted powers | k20, full | NOT-AFFECTED — the quoted numbers are themselves size-matched positive-control powers |
| A3 | Table 1; 48, 80, 101 | angles to the Southard cell-division direction | both | NOT-AFFECTED — no length and no change in distance; labelled descriptive |
| A4 | 16 | the Southard and Lu datasets listed as references | — | NOT-AFFECTED — no statistic |
| B1 | 18, 67 | "keeps 45.1% of the passage gap's length" | both | **no size-matched floor exists** |
| B2 | 20, 102 | the screen-to-PD32 offset, 0.65 to 1.05 times the passage gap | both | **no size-matched floor exists** |
| B3 | 70 | "14.4 and 14.3 times the donor-age gap" | full | **no size-matched floor exists** |
| B4 | 40, 47, 104 | shares of a direction removed by projection | both | NOT-AFFECTED — functions of an angle, not of a length |
| B5 | 79 | the mixing positive control, in gap fractions | k20 | NOT-AFFECTED — drawn at 96 and 157 cells; it *is* the positive control |
| C1 | 7, 26 | 80 of 340 called closer to early-passage cells; the four factors' −19.2 to −19.9 | k20 | NOT-AFFECTED — false-call rate 0.040–0.085 against a pre-set 0.10; no-effect median ≈ 0 |
| C2 | 26 | 73 called from the late-passage start | k20 | NOT-AFFECTED — same construction, that start's own floors |
| C3 | 19, 60 | only ZFX and 3 of 340 called in the full gene space | full | NOT-AFFECTED — the paper already states this is the space's insensitivity |
| C4 | 7, 29, 30 | the four factors' +11.4 to +15.3 against young donors; 5 of 340 called | k20 | NOT-AFFECTED — gate QUALIFIED, power 0.945 at a quarter-gap |
| **C5** | 30 | "and all 80 move away from them" | k20 | **PAPER2-DEPTH — corrected** |
| C6 | 27 | call rate by cell-count decile; the band of 60 | k20 | NOT-AFFECTED — the band *is* the size control |
| C7 | 37 | Spearman ρ = −0.864, distance against cell-cycle score | k20 | NOT-AFFECTED — no threshold of any kind; the paper labels it descriptive |
| C8 | 39, 40, 105, 106 | called counts after removing a cell-division direction | k20 | NOT-AFFECTED — floors recomputed after each removal, reported with the power collapse |
| C9 | 54 | 181 score younger, 116 not called closer | k20 | NOT-AFFECTED — the paper labels it a diagnostic |
| C10 | 67 | false-call rates 0.085 and 0.110 | k20 | NOT-AFFECTED — this *is* the floor validation |
| C11 | 19, 68 | half-gap power ≥ 0.940, quarter-gap ≥ 0.565 | k20 | NOT-AFFECTED — this *is* the plant gate |
| C12 | 31 | "any effect that points away … must increase the distance" | k20 | NOT-AFFECTED — a geometric identity; no length enters it |

### The three statements with no size-matched floor

None was checked. For each, the direction of the sampling-noise bias is known, and in each case it
**strengthens rather than weakens the paper's own caveat**, which is why none was corrected.

- **B1**, "keeps 45.1% of the passage gap's length" — a ratio of one vector's norm in two spaces.
  Noise inflates the full-gene norm far more than the 20-component one, so 45.1% **understates**
  the true retention. The paper uses it to argue that compression discards more noise than signal;
  the bias works against that argument, not for it.
- **B2**, the screen-to-PD32 offset at 0.65 to 1.05 times the passage gap — numerator and
  denominator are both noisy pseudobulk differences, and the denominator's PD14 half A is the
  smallest group in the comparison. No no-effect pool spans the two populations being compared, so
  the method used elsewhere in this check does not transfer. Paragraph 65 calls the offset "an
  upper bound on the technical difference", so noise inflation **strengthens** that disclosure.
- **B3**, "14.4 and 14.3 times the donor-age gap" — a distance from a 4,348-cell pseudobulk to two
  bulk centroids, not a displacement. Noise inflates it, so a correction could only **reduce** it,
  and the sentence's role is to warn that the GTEx-anchored test starts far outside the donor-age
  span.

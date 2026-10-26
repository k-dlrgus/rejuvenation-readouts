# FINDINGS_CORR2CHK — do paper 2's two "move away" clauses survive a size-matched no-effect baseline?

**Task:** CORR2CHK. **Status:** complete. No STOP fired. Reproduction: 41 gating checks, 0 mismatches.

**Primary readings, in the pre-registered wording:**

```
== A1b (paper 2 paragraph 31) ==
CONCLUSION  PAPER2-DEPTH: A1b, k20 (paper 2's 20-component space): the statement was reached partly by sampling noise; 0 of 4 units move away further than the 95th percentile of no-effect draws at their own cell count (about 0.2 expected by chance), 1 are inside the no-effect spread, and 3 move away less far than the median no-effect draw (about 2.0 expected by chance). Corrected wording: see FINDINGS_CORR2CHK.md, section "Corrected wording".
== C5 (paper 2 paragraph 30) ==
CONCLUSION  PAPER2-DEPTH: C5, k20 (paper 2's 20-component space): the statement was reached partly by sampling noise; 56 of 80 units move away further than the 95th percentile of no-effect draws at their own cell count (about 4.0 expected by chance), 23 are inside the no-effect spread, and 1 move away less far than the median no-effect draw (about 40.0 expected by chance). Corrected wording: see FINDINGS_CORR2CHK.md, section "Corrected wording".
== the 17 statements not checked ==
CONCLUSION  PAPER2-NOT-AFFECTED: A1a: its threshold was already matched to group size, so FLOORCHK's correction does not apply.
CONCLUSION  PAPER2-NOT-AFFECTED: A2: its threshold was already matched to group size, so FLOORCHK's correction does not apply.
CONCLUSION  PAPER2-NOT-AFFECTED: A3: its threshold was already matched to group size, so FLOORCHK's correction does not apply.
CONCLUSION  PAPER2-NOT-AFFECTED: A4: its threshold was already matched to group size, so FLOORCHK's correction does not apply.
CONCLUSION  PAPER2-NOT-AFFECTED: B4: its threshold was already matched to group size, so FLOORCHK's correction does not apply.
CONCLUSION  PAPER2-NOT-AFFECTED: B5: its threshold was already matched to group size, so FLOORCHK's correction does not apply.
CONCLUSION  PAPER2-NOT-AFFECTED: C1: its threshold was already matched to group size, so FLOORCHK's correction does not apply.
CONCLUSION  PAPER2-NOT-AFFECTED: C2: its threshold was already matched to group size, so FLOORCHK's correction does not apply.
CONCLUSION  PAPER2-NOT-AFFECTED: C3: its threshold was already matched to group size, so FLOORCHK's correction does not apply.
CONCLUSION  PAPER2-NOT-AFFECTED: C4: its threshold was already matched to group size, so FLOORCHK's correction does not apply.
CONCLUSION  PAPER2-NOT-AFFECTED: C6: its threshold was already matched to group size, so FLOORCHK's correction does not apply.
CONCLUSION  PAPER2-NOT-AFFECTED: C7: its threshold was already matched to group size, so FLOORCHK's correction does not apply.
CONCLUSION  PAPER2-NOT-AFFECTED: C8: its threshold was already matched to group size, so FLOORCHK's correction does not apply.
CONCLUSION  PAPER2-NOT-AFFECTED: C9: its threshold was already matched to group size, so FLOORCHK's correction does not apply.
CONCLUSION  PAPER2-NOT-AFFECTED: C10: its threshold was already matched to group size, so FLOORCHK's correction does not apply.
CONCLUSION  PAPER2-NOT-AFFECTED: C11: its threshold was already matched to group size, so FLOORCHK's correction does not apply.
CONCLUSION  PAPER2-NOT-AFFECTED: C12: its threshold was already matched to group size, so FLOORCHK's correction does not apply.
            (qualification, flag Part 3) C7 has no threshold of any kind; it is not affected because there is nothing to size-match.
            (qualification, flag Part 3) A2's quoted numbers are themselves size-matched positive-control powers.
== B1, B2, B3 ==
            no size-matched floor exists for these three length ratios; per the ruling they are recorded with the direction of their noise bias and no key fires for them.
```

## What this found, in plain words

**The inventory is the main result: 20 of the 22 clauses need no correction.** Paper 2 has 21 statements that depend on the Southard screen, on a displacement length or length ratio, or on a change in distance to a target; splitting paragraph 31's sentence into its two clauses, as the ruling on Q1 requires, makes 22 clause-level verdicts. Paper 2's own thresholds are drawn from control samples at each perturbation's own cell count, in the same space and from the same starting profile as the statement they gate. Paper 2 also nowhere repeats paper 1's corrected "more than twice the gap" claim, and it already applies FLOORCHK's own argument to its own full-gene-space result (paragraphs 19 and 60). Three length ratios have no size-matched floor at all and are recorded as such without being checked.

**Two clauses were checked. Both fire PAPER2-DEPTH, but they fail very differently, and the single key hides that.** Each is of the form "these units move away from young donors", and the pre-registered rule fails such a claim if even one unit is not distinguishable from no effect, because the published wording is universal.

- **A1b has no support at all.** 0 of 4 units exceed the 95th percentile of no-effect samples at their own cell count and 3 fall below the median. Under a true null about 0.2 of four would exceed the 95th percentile and about 2.0 would fall below the median, so what was observed is what no effect looks like. At these depths a no-effect sample's change in distance is itself positive, with a median of +1.3 to +2.7, so the positive sign the clause reports is the baseline's own sign.
- **C5's aggregate claim is strongly supported; only the word "all" fails.** 56 of the 80 units exceed their own 95th percentile, against about 4.0 expected by chance, and only 1 falls below its own median, against about 40.0 expected. That is a large, real departure for most of the 80. What does not hold is the universal quantifier: 24 of the 80 are not individually distinguishable from no effect, and they are the smaller groups.

Reading those two as one finding would be wrong. A1b's clause has nothing behind it. C5's clause is right about the direction and overstated only in claiming it of every one of the 80.

**Paragraph 31's sentence mixes two kinds of claim.** "Uncalled" is gated by a threshold matched to each factor's own cell count and is well powered at those depths (0.925 to 1.000 for a quarter-gap approach), so it stands and is better supported than the sentence claims. "Move away" is not gated by anything: it reports the sign of a measured change in distance against zero, when the matched no-effect baseline at those depths is itself positive. The two clauses sit in one sentence joined by "and", which is why the inventory had to split them.

**What the failures are not.** Neither result says any factor moves *toward* young donors; nothing here is evidence of rejuvenation by either definition. A1b says the Southard measurement cannot tell a departure from no effect at those four depths. C5 says the Sengstack measurement can tell a departure for most of the 80 but not for the smallest groups, so the claim holds of the set and not of every member. The call-based half of each sentence is untouched in both cases, and no factor named here is a candidate for anything.

Files, all in `results/corr2chk/`: `PREREG.flag`, `reproduction.csv`, `a1b_units.csv`, `c5_units.csv`, `c5_baseline.csv`, `readings.txt`, `summary.json`, `corr2chk_run.log`, `run_stdout.txt`, `run_stderr.txt`. Code: `src/corr2chk_run.py` (modes `run`, `verify`) and `src/corr2chk_findings.py`, which writes this file from those tables and computes nothing. No paper source was edited.

## Pre-registration (verbatim copy of `results/corr2chk/PREREG.flag`)

~~~~text
CORR2CHK pre-registration flag
Written 2026-10-26T11:21:38+09:00 (2026-10-26T02:21:38Z), before any statistic of this task
exists. No floor draw, no re-derived median, no per-perturbation comparison table and no
reading has been computed. No new code has been run on any expression data.

================================================================================
PART 0 — WHAT WAS KNOWN BEFORE THIS FILE WAS WRITTEN (disclosed first, because it
bears on the outcome)
================================================================================

This task's step 1 was an inventory of published statements, carried out and reported to the
user in chat BEFORE this flag was written, at the user's instruction ("Option 2. Finish the
step-1 inventory and come back before writing the flag"). The inventory read stored result
files. It therefore established, before this flag existed, what the answer for both checked
statements is very likely to be. That is disclosed here in full rather than discovered below.

[K1] OUTCOME KNOWLEDGE FOR A1b. The inventory read results/south7/stage_a1_depth_nulls.csv,
  which stores, per cell count and space, the 5th percentile, median, sd and 95th percentile
  of the change in distance for 200 no-effect draws. In S2A_k20, at the four factors' own
  Southard cell counts, it read:
     83 cells  floor_p05 +0.467661  median +2.694518  sd 1.743163  p95 +5.932085
     96 cells  floor_p05 -0.106862  median +2.520107  sd 1.787256  p95 +5.643203
    129 cells  floor_p05 -0.995048  median +1.514257  sd 1.573311  p95 +4.179653
    140 cells  floor_p05 -0.694436  median +1.294052  sd 1.415215  p95 +3.967368
  and from results/south8/d3c_all_1836_umi_calls.csv the four factors' stored changes in
  distance in the same space: STAT3 (83 cells) +2.319734, E2F3 (96) +1.562149,
  ZFX (129) +2.441500, EZH2 (140) +0.481121.
  Comparing those by eye — which the inventory did, and reported — three of the four sit
  BELOW the no-effect median at their own cell count, and none of the four reaches its own
  P95. So the session already expects PAPER2-DEPTH to fire for A1b. The computation below
  formalises, tabulates and gates that expectation; it does not discover it.

[K2] OUTCOME KNOWLEDGE FOR C5. The inventory joined results/seng3/d4_table.csv to
  results/axang/d3b_table.csv and read that the 80 perturbations called closer to
  early-passage cells have young-donor changes in distance from +0.927 to +18.348, median
  +8.478, and that exactly one of them (CRA_TFCP2, 120 cells, +0.927422) falls below the
  largest no-effect median stored for that test (+1.50 at 96 cells, CRI). It also read that
  55 of the 80 have fewer than 96 cells and 18 fewer than 60, depths at which the young-donor
  no-effect median is larger and is NOT stored. So the session already expects at least one
  of the 80 to fall below its own baseline, and expects the count to be small. Whether it
  stays at one is not known, because the baseline at the smaller depths has never been
  recorded.

[K3] The exact wording of all four possible conclusions was fixed by the user in the task
  message, before the inventory ran and before anything in [K1] or [K2] was read. The
  per-statement rulings Q1-Q4 reproduced in Part 2 were given after the inventory was
  reported, with [K1] and [K2] in front of the user.

[K4] Everything else read before this file: CLAUDE.md; FINDINGS_FLOORCHK.md (whole);
  CORRECTIONS_PAPER1.md (whole); PROGRESS_FLOORCHK.md; FINDINGS_AXANG.md (lines 1-60,
  225-300, 370-520, 560-700, 720-900); FINDINGS_SOUTH8.md (lines 1-130, 180-430, 540-600);
  FINDINGS_COORDSPACE.md (lines 1-30, 540-600); README.md (lines 15-45); .gitignore;
  paper/paper2_passage_draft_v0.7.docx (all 117 body paragraphs, text extracted read-only to
  a scratch file, not edited) and paper/two_definitions_of_young_preprint_v1.pdf (text
  extracted read-only); paper/rejuvenation_readouts_preprint_v1.2.docx (paragraphs 138-145
  only, read-only). Code read, none modified: src/corr1_run.py (head), src/axang_run.py,
  src/seng3_run.py (parts), src/seng4_run.py (parts), src/south7_run.py (parts),
  src/south8_run.py (parts), src/south6rec_run.py (parts), src/seng2_common.py (seeds).
  Stored tables read: results/south7/stage_a1_depth_nulls.csv;
  results/south8/{d3c_all_1836_umi_calls.csv, d1_plant_curves.csv, d1_plant_curves_all.csv,
  d2_calibration.csv, d2_verdicts.csv}; results/seng3/{d1_null_validation.csv, d2_plant.csv,
  d4_table.csv}; results/seng4/{d2_plant_residualized.csv, d3_deciles.csv};
  results/axang/{d3a_plant.csv, d3b_table.csv};
  results/coordspace/{stage3_south6_factors.csv, stage3_south6_summary.csv};
  results/paper2_figs_v06/table1.csv.

[K5] Values computed by hand or by one-line reads during the inventory, beyond the stored
  columns themselves, and so also known before this flag: the count of the 80 called
  perturbations with fewer than 96 and fewer than 60 cells (55 and 18); the min, 5th
  percentile, median and max of the 80's young-donor changes in distance; that all 80 of the
  passage test's called perturbations lie below their own stored null_median_i; that the
  minima over cell count of SENG4's residualised half-gap powers are 0.895 (nothing removed),
  0.010 (primary axis), 0.040 (PD26 axis) and 0.885 (gene-list axis) from the own-control
  start. None of these is a result of this task. They are reported in the findings as
  inventory readings and are labelled as such.

================================================================================
PART 1 — TASK, VERBATIM AS GIVEN IN CHAT
================================================================================

Read CLAUDE.md first; its rules apply. Read FINDINGS_FLOORCHK.md and CORRECTIONS_PAPER1.md.

Task: `CORR2CHK`. Deliverables: `results/corr2chk/PREREG.flag`, `FINDINGS_CORR2CHK.md`,
`PROGRESS_CORR2CHK.md`. One local commit, no push. Do not modify any existing FINDINGS_*.md,
CORRECTIONS_PAPER1.md, or any file under results/floorchk/, results/corr1/ or results/south6/.

## The question

FLOORCHK showed that in the Southard screen's full-gene space, displacement length is
dominated by sampling noise at the screen's cell counts: non-targeting cells with no effect
read a median 3.76 gap lengths at 106 cells, and after correction the perturbations' median is
0.61.

Paper 2 ("Two definitions of young give opposite answers in a fibroblast rejuvenation screen")
states that the four validated factors are uncalled and move away from young donors in the
Southard screen. This task asks whether that statement, and any other paper-2 statement
resting on Southard displacement length or on a change in distance computed in a space where
the floor exceeds the effect, survives the same check.

## Step 1 — inventory before any statistic

List every paper-2 statement that depends on (a) the Southard screen, (b) a displacement
length or length ratio, or (c) a change in distance to a target, with file, paragraph and
verbatim quote. Include main text, supplementary notes, figures and tables. For each, record
which dataset, which coordinate space, and which starting profile it uses.

Paper 2's own tests are in the Sengstack data, which has its own control-sample floor already
built into its thresholds. Distinguish clearly, for each statement, whether its floor was
already matched to group size (as the Sengstack call thresholds are) or not (as the Southard
full-gene ratio was not). A statement whose threshold was already size-matched does NOT need
the FLOORCHK correction; say so and move on.

## Step 2 — pre-register, then check only what step 1 shows needs checking

Write the flag after step 1 and before any statistic, naming exactly which statements are
being checked, the statistic, the floor construction, the gate and the exact wording of every
conclusion.

For each statement needing a check, apply the FLOORCHK method unchanged: a size-matched floor
from non-targeting cells drawn at the statement's own group size, 200 draws, reference = pool
minus the draw, the same noise-corrected statistic, and the same validation gate with planted
displacements. Reuse the FLOORCHK code and the stored Southard intermediates; recompute
nothing that is already stored, and reproduce every stored number used before computing
anything new. Seeds: fresh, named in the flag.

Conclusions, worded now, per statement:
  "PAPER2-STANDS: <statement id>, <space>: the statement survives the size-matched floor;
   <numbers>."
  "PAPER2-DEPTH: <statement id>, <space>: the statement was reached partly by sampling noise;
   <numbers>. Corrected wording: <text>."
  "PAPER2-NOT-AFFECTED: <statement id>: its threshold was already matched to group size, so
   FLOORCHK's correction does not apply."
  "PAPER2-UNRESOLVED: <statement id>, <space>: the correction failed its gate; only raw
   values and floor are reported."

## Step 3 — if anything fails

Write the corrected wording for each failing statement in FINDINGS_CORR2CHK.md, in the same
plain style as CORRECTIONS_PAPER1.md. Do not edit the paper-2 source and do not write a
correction notice in this task; report what one would have to say.

## Step 4 — paper 1, paragraph 142

Separately, report the exact replacement wording for the paragraph-142 bullet that
CORRECTIONS_PAPER1.md identifies as carrying no information ("In space C both perturbations
move away from the young centroid (+106.6 and +113.8, q = 1.0)"). Do not edit any file; give
the wording only, so it can be folded into the paper-1 revision before publication.

## Hard rules
No conclusion about any named gene or factor beyond what the existing papers already state.
Compute nothing in a space paper 2 does not use.
Disclose every deviation. If anything is ambiguous, STOP and ask.

--- second user message, verbatim ---

Option 2. Finish the step-1 inventory and come back before writing the flag.

For every statement you propose as PAPER2-NOT-AFFECTED, give the evidence, not just the
ruling:
- the file and line where its floor or threshold is defined, and the group size it was drawn
  at;
- whether that floor was built in the same coordinate space and from the same starting profile
  as the statement;
- the file and line of the plant gate it relied on, and its numbers.

Also list any statement where you are not sure, separately, with what would settle it.

Write nothing to results/corr2chk/ yet. Step 4 (the paper-1 paragraph-142 wording) does not
depend on the inventory; include it in the same report.

--- third user message, verbatim ---

Before the rulings: go back over the inventory and show me the full statement list with the
per-statement floor evidence for the 18 you grouped, not just the four classes.

================================================================================
PART 2 — RULINGS BY THE USER, VERBATIM, AFTER THE INVENTORY AND BEFORE THIS FILE
================================================================================

Rulings on Q1-Q4.

Q1 — a statement may take two keys. Split A1 into A1a ("uncalled", NOT-AFFECTED on the
evidence given) and A1b ("move away from young donors", to be checked). Record both clause
verdicts against the same source sentence, and say in the findings that the sentence as
published mixes a size-matched claim with an unmatched one.

Q2 — yes, a change-in-distance comparison against a size-matched no-effect baseline at the
statement's own cell count substitutes for FLOORCHK's quadrature correction, and is preferred
where it is available. Reason, stated in the flag: FLOORCHK's correction exists to remove a
noise contribution to a LENGTH; where the statistic is already a change in distance measured
against no-effect draws of the same size, the noise is subtracted by construction and the
quadrature step would double-count it. Where no such baseline exists, use the FLOORCHK
correction as written. Name in the flag which of the two applies to A1b and to C5.

Q3 — read stored floors where they exist and are keyed to the statement's own cell count;
draw fresh only where nothing stored fits. Reproduce every stored value used before using it.

Q4 — yes, check C5. Re-deriving AXANG's discarded per-size medians is permitted: it recomputes
a quantity that was computed and thrown away, not one that is stored, so there is no stored
value to contradict. Do it by re-running AXANG's existing loop unchanged except for recording
the median alongside the 5th percentile. Fresh seeds, named in the flag. Before using the
re-derived values, reproduce AXANG's stored 5th percentiles from the same re-run and report
the comparison; if they do not reproduce, STOP and report.

Also fixed now, for both A1b and C5: report the result per perturbation at its own cell count,
not only in aggregate, and state how many move away by less than their own no-effect median.

The three unchanged uncertain statements (B1, B2, B3) are recorded as having no size-matched
floor, with the direction of the noise bias stated, and are not checked. Say in the findings
that each strengthens rather than weakens the paper's own caveat.

Proceed to the flag and the run after these rulings.

================================================================================
PART 3 — THE INVENTORY THIS FLAG RESTS ON (step 1, settled before this file)
================================================================================

Source of every quote: paper/paper2_passage_draft_v0.7.docx, whose text matches
paper/two_definitions_of_young_preprint_v1.pdf. Paragraph numbers are docx body-paragraph
indices, as CORR1 used for paper 1. Table 1 is a docx table and has no body-paragraph index;
its values are also stored at results/paper2_figs_v06/table1.csv.

21 statements were found in categories (a) Southard, (b) a displacement length or length
ratio, (c) a change in distance to a target. The full list, with verbatim quotes, dataset,
space, starting profile and per-statement floor evidence, was reported to the user in chat
over two messages and is reproduced in FINDINGS_CORR2CHK.md. Space abbreviations below:
k20 = paper 2's 20-component space = SOUTH6/7/8's S2A_k20 = results/south7/basis.npz key
L_mid_full, the first 20 PCs of 308 GTEx donors aged 40-59, read frozen by src/seng3_run.py:64,
src/axang_run.py:87 and src/south8_run.py:439-441 alike; full = the 23,485-gene space = S1 =
paper 1's "space C".

CHECKED BY THIS TASK (2 of 21):
 A1b  P2 para 31, the clause "and move away from young donors in the Southard screen".
      Southard; k20 (the space is not named in the published sentence; the numbers behind it
      are k20); start = GTEx old centroid, target = GTEx young centroid.
 C5   P2 para 30, the clause "and all 80 move away from them" (them = young donors).
      Sengstack; k20; start = GTEx old centroid, target = GTEx young centroid.

NOT CHECKED, PAPER2-NOT-AFFECTED (17 of 21): A1a (para 31 "uncalled"), A2 (paras 18, 77),
 A3 (Table 1 row 1; paras 48, 80, 101), A4 (para 16), B4 (paras 40, 47, 104), B5 (para 79),
 C1 (paras 7, 26), C2 (para 26), C3 (paras 19, 60), C4 (paras 7, 29, 30), C6 (para 27),
 C7 (para 37), C8 (paras 39, 40, 105, 106), C9 (para 54), C10 (para 67), C11 (paras 19, 68),
 C12 (para 31). The evidence for each — the file and line of its floor, the group size it was
 drawn at, whether the floor shares the statement's space and starting profile, and the file,
 line and numbers of its plant gate — is in FINDINGS_CORR2CHK.md and was reported to the user
 before this flag. Two of the 17 carry a qualification recorded now: C7 has no threshold of
 any kind and is NOT-AFFECTED because there is nothing to size-match, not because anything was
 size-matched; A2's quoted numbers are themselves size-matched positive-control powers.

NOT CHECKED, NO SIZE-MATCHED FLOOR EXISTS (3 of 21): B1 (paras 18, 67, "45.1% of the passage
 gap's length"), B2 (paras 20, 102, the screen-to-PD32 offset as a multiple of the passage
 gap), B3 (para 70, "14.4 and 14.3 times the donor-age gap"). Per the ruling these are
 recorded with the direction of their noise bias and are not checked. No key fires for them.

Three findings of the inventory that bear on scope, recorded now:
 - Paper 2 nowhere repeats paper 1's corrected "more than twice the gap" claim. A search of
   the source for "twice", "2.10", "1,836" and "overshoot" returns only para 13, "apply that
   readout to the Sengstack data twice". CORR1's correction does not propagate into paper-2
   text.
 - Paper 2 already applies FLOORCHK's own argument to its own full-gene-space result, in
   paras 19 and 60. A check there would confirm the published text, not correct it.
 - Paper 2's k20 and the Southard S2A_k20 are the same frozen basis, read not refit, so A1b
   and C4 are in one space and the "compute nothing in a space paper 2 does not use" rule is
   satisfied for both.

================================================================================
PART 4 — THE STATISTIC, FOR BOTH CHECKED STATEMENTS
================================================================================

[D1] Which of the two corrections applies, as Q2 requires it to be named.
  A1b: the change-in-distance comparison. A size-matched no-effect baseline exists at the
    statement's own cell counts (results/south7/stage_a1_depth_nulls.csv, S2A_k20). FLOORCHK's
    quadrature correction is NOT used and would be wrong here: the statistic is already a
    change in distance taken against draws of the same size, so the noise is subtracted by
    construction and a quadrature step would double-count it.
  C5: the change-in-distance comparison. The baseline does not exist in stored form and is
    re-derived under Q4, from the same construction, at each perturbation's own cell count.
    FLOORCHK's quadrature correction is not used, for the same reason.
  Neither checked statement is a length, so the FLOORCHK correction as written is used in
  neither. No statement in this task uses it.

[D2] The decision statistic, both statements.
  For a perturbation with n cells, let d be its stored change in distance to the GTEx young
  centroid from the GTEx old centroid in k20, and let the no-effect draws at n cells give
  F5(n) (5th percentile), M(n) (median) and P95(n). The statement under test is of the form
  "it moves away". The decision statistic is the pair
      (d - M(n),  whether d >= P95(n)).
  A perturbation SUPPORTS "moves away" iff d >= P95(n): it moves away further than 95% of
  no-effect draws of its own size. A perturbation CONTRADICTS "moves away" iff d < M(n): it
  moves away less far than the median no-effect draw of its own size. A perturbation between
  M(n) and P95(n) is INSIDE THE NOISE: its sign is the floor's sign and it supports nothing.
  F5(n) is reported beside these because it is the floor the published calls were made
  against; it gates nothing here.

[D3] Chance baseline, fixed now (CLAUDE.md).
  P95 is exceeded by 5% of no-effect draws by construction, and M(n) by 50%. So under a true
  null of no effect anywhere, 5% of perturbations are expected to SUPPORT and 50% to
  CONTRADICT. Both expectations are reported beside every count.

[D4] Per-perturbation reporting, as the ruling fixes.
  Both statements are reported per perturbation at its own cell count, with n, d, F5(n),
  M(n), P95(n), the verdict of D2, and the count that CONTRADICTS, beside the D3 expectation.
  A1b's units are the four named factors (the sentence is about them). C5's units are the 80
  perturbations the published sentence is about.

================================================================================
PART 5 — WHERE EVERY NUMBER COMES FROM
================================================================================

[D5] A1b. Nothing is drawn. Per Q3 the stored floor is read.
  Units: the four gene names paper 2 names, as they appear in the Southard screen: E2F3,
    EZH2, STAT3, ZFX. Cell counts 96, 140, 83, 129 respectively.
  d: results/south8/d3c_all_1836_umi_calls.csv column delta_S2A_k20.
  F5, M, P95: results/south7/stage_a1_depth_nulls.csv, space S2A_k20, at the unit's own
    n_cells. Built by src/south7_run.py:319-324 from 200 draws (R_NULL, src/south7_run.py:51)
    of non-targeting cells at that n, with the reference being the non-targeting pool minus
    the draw (src/south7_run.py:323), which is FLOORCHK's U3 rule and the task's required
    construction.
  Recovery power at the unit's own depth: results/south8/d1_plant_curves_all.csv column
    power_25, space S2A_k20, built by src/south8_run.py:157-193 with the same reference
    (stated at :160-161, implemented at :176-178). Reported beside every verdict per CLAUDE.md.
  Space and start: S2A_k20, origin = GTEx old centroid, target = GTEx young centroid, from
    results/toward/anchors.npz, read not refit.

[D6] C5. The baseline is re-derived, per Q4.
  Units: the 80 perturbations with called_delta_i true in results/seng3/d4_table.csv. The
    count must be exactly 80 or STOP (S3).
  d: results/axang/d3b_table.csv column delta, which is the young-donor change in distance in
    k20 from the GTEx old centroid.
  F5, M, P95: re-derived by re-running AXANG's own D3b floor loop (src/axang_run.py:508-518,
    with draws at :471-478 and the delta at :481-486) with one change only, recording the
    median and the 95th percentile alongside the 5th percentile. The loop's reference is
    pool minus draw (src/axang_run.py:476), 200 reps, at every distinct cell count in each
    screen.
  How the loop is re-run, disclosed because it is a deviation: axang_run.d3() is a single
    function that writes into results/axang/, so it is not called. Its prologue and floor loop
    are repeated line for line in src/corr2chk_run.py with the helpers imported from
    axang_run, seng3_run and seng2_common unchanged. The proof that the repetition is faithful
    is R5 below: the re-derived 5th percentiles must equal the stored ones.
  Space and start: k20 and the GTEx old centroid, identical to A1b's, from the same anchors
    file and the same frozen basis.

[D7] Seeds, and a tension in the ruling, disclosed.
  The ruling asks for the loop "unchanged except for recording the median", for "fresh seeds",
  and for the stored 5th percentiles to reproduce "from the same re-run". Those three cannot
  all hold at once: a 5th percentile of 200 resampled draws is a random quantity, so a run at
  a fresh seed cannot reproduce the stored value exactly, while a run that reproduces it
  exactly is by definition not at a fresh seed. This is resolved as follows, and the
  resolution is reported prominently in the findings:
   - PRIMARY re-run, seed unchanged at AXANG's own K.SEED = 20260914 (src/seng2_common.py:38),
     which is what "the loop unchanged" means, since the seed is an argument of the loop
     (src/axang_run.py:471). This run supplies M(n) and P95(n) and must reproduce the stored
     F5(n) exactly (R5). It is the run every C5 number comes from.
   - STABILITY replicate at a fresh seed, 20261030, named here as the ruling requires. It
     gates nothing and decides nothing. It is reported so that the Monte Carlo spread of
     F5, M and P95 at 200 draws is visible, and so that the primary run's agreement with the
     stored floor is not mistaken for evidence that 200 draws are precise.
  No other draw is made by this task. A1b draws nothing.

================================================================================
PART 6 — REPRODUCTION BEFORE ANYTHING NEW (CLAUDE.md; gating)
================================================================================

Run first, by src/corr2chk_run.py, into results/corr2chk/reproduction.csv. All gating unless
marked otherwise. Counts exact. Numeric tolerances in D11.

 R1 A1b, the four factors' change in distance. Recompute each from
    results/south6/stage1.npz (key d_f), results/south7/basis.npz (key L_mid_full) and
    results/toward/anchors.npz, by SOUTH7's own formula
      delta = ||(c_old - c_young) L + d L|| - ||(c_old - c_young) L||,
    and compare with the stored delta_S2A_k20. Also compare with the independent stored copy
    in results/south7/stage_a1_observed.csv if that file carries it.
 R2 A1b, the floor. For each of the four, the stored depth_floor_S2A_k20 in
    results/south8/d3c_all_1836_umi_calls.csv must equal floor_p05 in
    results/south7/stage_a1_depth_nulls.csv at the same n_cells and space. This is the check
    that the floor read in D5 is the floor the published call was made against.
 R3 A1b, cell counts. The four factors' n_cells must equal
    results/south6/stage0_cells_per_target.csv. Exact.
 R4 A1b, the four gene names must each appear exactly once in the Southard table. Exact.
 R5 C5, the re-derived floor. The PRIMARY re-run's 5th percentile at every (modality, cell
    count) must equal the stored floor in results/axang/d3b_table.csv. Exact to D11's
    rebuild tolerance. This is the gate the ruling requires; failure is STOP S2.
 R6 C5, the young-donor changes in distance. The re-run recomputes d for all 340
    perturbations as AXANG does; every one must equal the stored delta in
    results/axang/d3b_table.csv.
 R7 C5, the units. results/seng3/d4_table.csv must have exactly 340 rows and exactly 80 with
    called_delta_i true, and must join to results/axang/d3b_table.csv on
    (modality, perturbation, n_cells) with 340 matches and none lost.
 R8 C5, the published count. results/axang/d3b_table.csv must reproduce AXANG's stored
    summary: 5 called of 340, 17 with distance decreased, and the four named factors' changes
    in distance +11.448, +9.022, +5.327, +15.336. Exact to D11.
 R9 Descriptive, not gating. The stored dose-0 medians in results/axang/d3a_plant.csv at
    96, 107, 128 and 157 cells are reported beside the re-derived M(n) at the same depths, as
    a second check that the re-run matches AXANG's own construction.

If any gating check fails: STOP S1, naming every mismatch. A STOP is never resumed and never
restarted with an adjustment.

================================================================================
PART 7 — CONCLUSIONS, WORDED NOW
================================================================================

The four wordings are the user's, fixed in the task message before the inventory ran. They
are used verbatim, with <statement id>, <space> and <numbers> filled in.

  "PAPER2-STANDS: <statement id>, <space>: the statement survives the size-matched floor;
   <numbers>."
  "PAPER2-DEPTH: <statement id>, <space>: the statement was reached partly by sampling noise;
   <numbers>. Corrected wording: <text>."
  "PAPER2-NOT-AFFECTED: <statement id>: its threshold was already matched to group size, so
   FLOORCHK's correction does not apply."
  "PAPER2-UNRESOLVED: <statement id>, <space>: the correction failed its gate; only raw
   values and floor are reported."

[D8] Firing rules, fixed now.
  For a checked statement of the form "these units move away", with m units:
   - PAPER2-STANDS fires iff every one of the m units SUPPORTS (d >= P95(n)).
   - PAPER2-DEPTH fires iff at least one unit CONTRADICTS (d < M(n)) or at least one unit is
     INSIDE THE NOISE. <numbers> gives m, the count SUPPORTING, the count INSIDE THE NOISE,
     the count CONTRADICTING, the D3 chance expectations, and the recovery power at the
     depths in question.
   - The two are exclusive and exhaustive over the m units, so exactly one fires per checked
     statement.
   - PAPER2-UNRESOLVED fires instead, for that statement, if its reproduction gate fails in a
     way Part 6 routes to S2 rather than S1, or if the re-derived baseline cannot be computed.
   - PAPER2-NOT-AFFECTED fires for each of the 17 statements listed in Part 3, on the evidence
     recorded there. It is not a computed result and no number of this task gates it.
  <space> is filled as "k20 (paper 2's 20-component space)" for both checked statements.
  A1b's and C5's readings are both primary. There is no secondary reading in this task.

[D9] If PAPER2-DEPTH fires, step 3 applies: FINDINGS_CORR2CHK.md carries the corrected
  wording for that statement, in CORRECTIONS_PAPER1.md's plain style. The paper-2 source is
  not edited and no correction notice is written, per the task.

[D10] A1a is recorded as PAPER2-NOT-AFFECTED against the same source sentence as A1b, per Q1,
  and the findings state that paragraph 31's sentence as published mixes a size-matched claim
  ("uncalled") with an unmatched one ("move away").

================================================================================
PART 8 — TOLERANCES, STOPS, HONESTY OBLIGATIONS, PROCESS
================================================================================

[D11] Tolerances. Counts exact. Values read from one stored file and compared with the same
  value in another stored file: relative 1e-12. Values recomputed from stored vectors or a
  stored basis: relative 1e-9. Values rebuilt by re-running a seeded resampling loop against
  the stored output of that same loop at the same seed: relative 1e-9, and 0 expected.
  Published rounded values quoted in the paper: compared at the precision printed.

STOP conditions, fixed now:
 [S1] reproduction_failed — any gating check in Part 6 fails. Every mismatch named.
 [S2] axang_floor_not_reproduced — R5 fails. Named separately because the ruling names this
      case and its consequence: C5 is then PAPER2-UNRESOLVED and no re-derived median is used.
 [S3] units_not_identified — R7 fails: the 80 cannot be identified, or the join loses rows.
 [S4] flag_ambiguous — this flag is ambiguous or self-contradictory in a way Parts 2 to 7 do
      not settle. Write the findings with the STOP key, name the ambiguity, ask.
A STOP is never resumed and never restarted with an adjustment.

Honesty obligations (CLAUDE.md), applied:
 - Chance baseline beside every count: D3 fixes it at 5% SUPPORTING and 50% CONTRADICTING
   under a true null, and both are printed beside every count.
 - Recovery power beside any null: the per-unit power_25 at each unit's own depth is reported
   for A1b, and AXANG's D3a gate (QUALIFIED, p25 = 0.945, p50 = 1.000) for C5. A1a's
   "uncalled" verdict is reported with the power that qualifies it, never as a bare count.
 - Shortlist against cell count: both checked statements are reported per unit with n beside
   every verdict, and the relation between n and the verdict is stated. No unit is presented
   as a candidate for anything.
 - Circularity, disclosed where the results appear: for both statements the non-targeting
   pool of the relevant screen defines zero for every change in distance and is also the pool
   every no-effect draw is taken from; and for C5 the 80 are selected by a statistic (the
   passage test) that uses the same control cells as the young-donor statistic being tested.
 - What this does not cover: cell-sampling noise only, not guide-to-guide, batch or
   biological variation. Southard and Sengstack are each one neonatal line with no aged
   donor, so nothing here transfers to aged human cells. Nothing here tests any factor's
   phenotype. The GTEx-anchored construction transplants a displacement measured in single
   cells onto a bulk post-mortem origin, as both papers state.
 - No conclusion about any named gene or factor beyond what the existing papers already
   state. A1b's units carry gene names because the published sentence names them; the result
   is about the sentence, not about the genes.

Process:
 Code: src/corr2chk_run.py (modes: run, verify). Interpreter: Anaconda Python, run with -I.
 Imports, unchanged and not modified: axang_run, seng3_run, seng2_common, seng_run,
   south7_run helpers as needed, toward_run. Importing axang_run and seng3_run creates
   results/axang/ and results/seng3/ if absent via mkdir(exist_ok=True); both already exist
   and no file in them is written.
 Results: results/corr2chk/ only. Not modified: any FINDINGS_*.md, CORRECTIONS_PAPER1.md,
   anything under results/floorchk/, results/corr1/, results/south6/, results/south7/,
   results/south8/, results/axang/, results/seng3/, results/seng4/, results/coordspace/, or
   any paper source.
 Deliverables: this file; FINDINGS_CORR2CHK.md with this file reproduced verbatim at the top;
   PROGRESS_CORR2CHK.md. One local commit, no push.
 The pre-existing uncommitted change to CITATION.cff and the untracked
   results/paper2_figs_plain/ are not this task's work and are left out of its commit.
 Step 4 of the task (the paper-1 paragraph-142 replacement wording) computes nothing. Every
   number in it is read from results/coordspace/ and from FINDINGS_FLOORCHK.md and is quoted
   with its source. It was drafted and reported to the user before this flag was written, and
   is reproduced unchanged in FINDINGS_CORR2CHK.md.
~~~~

## Reproduction of stored numbers, before anything new

`python -I src/corr2chk_run.py run` writes `results/corr2chk/reproduction.csv`. It has **41 gating checks and 0 mismatches**, so no STOP fired. 6 further rows are descriptive and gate nothing.

| check | item | expected | observed | tol | match | gating |
|---|---|---|---|---|---|---|
| R4 Southard table | rows for E2F3 | 1.0 | 1.0 | exact | True | True |
| R4 Southard table | rows for EZH2 | 1.0 | 1.0 | exact | True | True |
| R4 Southard table | rows for STAT3 | 1.0 | 1.0 | exact | True | True |
| R4 Southard table | rows for ZFX | 1.0 | 1.0 | exact | True | True |
| R1 basis | ||O1-Y1|| in S2A_k20 (stored rnorm vs ||r||) | 27.92915445 | 27.92915445 | 1e-09 | True | True |
| R1 E2F3 | change in distance, rebuilt vs stored delta_S2A_k20 | 1.562149403 | 1.562149403 | 1e-09 | True | True |
| R1 E2F3 | stored delta in stage_a1_observed vs d3c table | 1.562149403 | 1.562149403 | 1e-12 | True | True |
| R3 E2F3 | rows in SOUTH6 stage 0 cells-per-target | 1.0 | 1.0 | exact | True | True |
| R3 E2F3 | cells per target vs SOUTH6 stage 0 | 96.0 | 96.0 | exact | True | True |
| R2 E2F3 | depth floor, SOUTH7 nulls vs SOUTH8 table | -0.1068623398 | -0.1068623398 | 1e-12 | True | True |
| R2 E2F3 | draws behind the baseline | 200.0 | 200.0 | exact | True | True |
| R1 EZH2 | change in distance, rebuilt vs stored delta_S2A_k20 | 0.4811212715 | 0.4811212715 | 1e-09 | True | True |
| R1 EZH2 | stored delta in stage_a1_observed vs d3c table | 0.4811212715 | 0.4811212715 | 1e-12 | True | True |
| R3 EZH2 | rows in SOUTH6 stage 0 cells-per-target | 1.0 | 1.0 | exact | True | True |
| R3 EZH2 | cells per target vs SOUTH6 stage 0 | 140.0 | 140.0 | exact | True | True |
| R2 EZH2 | depth floor, SOUTH7 nulls vs SOUTH8 table | -0.6944363312 | -0.6944363312 | 1e-12 | True | True |
| R2 EZH2 | draws behind the baseline | 200.0 | 200.0 | exact | True | True |
| R1 STAT3 | change in distance, rebuilt vs stored delta_S2A_k20 | 2.319733597 | 2.319733597 | 1e-09 | True | True |
| R1 STAT3 | stored delta in stage_a1_observed vs d3c table | 2.319733597 | 2.319733597 | 1e-12 | True | True |
| R3 STAT3 | rows in SOUTH6 stage 0 cells-per-target | 1.0 | 1.0 | exact | True | True |
| R3 STAT3 | cells per target vs SOUTH6 stage 0 | 83.0 | 83.0 | exact | True | True |
| R2 STAT3 | depth floor, SOUTH7 nulls vs SOUTH8 table | 0.467660657 | 0.467660657 | 1e-12 | True | True |
| R2 STAT3 | draws behind the baseline | 200.0 | 200.0 | exact | True | True |
| R1 ZFX | change in distance, rebuilt vs stored delta_S2A_k20 | 2.441499656 | 2.441499656 | 1e-09 | True | True |
| R1 ZFX | stored delta in stage_a1_observed vs d3c table | 2.441499656 | 2.441499656 | 1e-12 | True | True |
| R3 ZFX | rows in SOUTH6 stage 0 cells-per-target | 1.0 | 1.0 | exact | True | True |
| R3 ZFX | cells per target vs SOUTH6 stage 0 | 129.0 | 129.0 | exact | True | True |
| R2 ZFX | depth floor, SOUTH7 nulls vs SOUTH8 table | -0.9950479374 | -0.9950479374 | 1e-12 | True | True |
| R2 ZFX | draws behind the baseline | 200.0 | 200.0 | exact | True | True |
| R7 units | rows in seng3 d4_table | 340.0 | 340.0 | exact | True | True |
| R7 units | rows in axang d3b_table | 340.0 | 340.0 | exact | True | True |
| R7 units | called closer to early-passage cells (own-control start) | 80.0 | 80.0 | exact | True | True |
| R7 units | rows surviving the join | 340.0 | 340.0 | exact | True | True |
| R8 axang | called closer to young donors, of 340 | 5.0 | 5.0 | exact | True | True |
| R8 axang | distance decreased, of 340 | 17.0 | 17.0 | exact | True | True |
| R8 CRA_E2F3 | stored young-donor change in distance (3 dp) | 11.448 | 11.448 | 1e-09 | True | True |
| R8 CRA_EZH2 | stored young-donor change in distance (3 dp) | 9.022 | 9.022 | 1e-09 | True | True |
| R8 CRI_STAT3 | stored young-donor change in distance (3 dp) | 5.327 | 5.327 | 1e-09 | True | True |
| R8 CRI_ZFX | stored young-donor change in distance (3 dp) | 15.336 | 15.336 | 1e-09 | True | True |
| R6 axang | max relative error, 340 recomputed changes in distance | 0.0 | 2.457726482e-16 | 1e-09 | True | True |
| R5 C5 baseline | max relative error, re-derived vs stored 5th percentile | 0.0 | 2.0042975e-16 | 1e-09 | True | True |
| R9 CRA 96 | dose-0 median, re-derived vs AXANG d3a | 1.066771004 | 1.066771004 | 1e-09 | True | False |
| R9 CRI 96 | dose-0 median, re-derived vs AXANG d3a | 1.502393873 | 1.502393873 | 1e-09 | True | False |
| R9 CRA 107 | dose-0 median, re-derived vs AXANG d3a | 1.242412619 | 1.242412619 | 1e-09 | True | False |
| R9 CRI 107 | dose-0 median, re-derived vs AXANG d3a | 1.050074009 | 1.050074009 | 1e-09 | True | False |
| R9 CRA 128 | dose-0 median, re-derived vs AXANG d3a | 1.291713124 | 1.291713124 | 1e-09 | True | False |
| R9 CRI 157 | dose-0 median, re-derived vs AXANG d3a | 0.8507320729 | 0.8507320729 | 1e-09 | True | False |

The two checks that carry the most weight:

- **R2** shows the baseline read for A1b is the same floor the published call was made against: `depth_floor_S2A_k20` in `results/south8/d3c_all_1836_umi_calls.csv` equals `floor_p05` in `results/south7/stage_a1_depth_nulls.csv` at each factor's own cell count, to 1e-12.
- **R5** shows the re-derived C5 baseline is AXANG's own: re-running AXANG's floor loop at its own seed reproduces every stored 5th percentile. That is the proof the repeated loop is faithful, and the ruling made it a STOP condition.

## The inventory (step 1)

Quotes are from `paper/paper2_passage_draft_v0.7.docx`, whose text matches `paper/two_definitions_of_young_preprint_v1.pdf`. Paragraph numbers are docx body-paragraph indices. **k20** is paper 2's 20-component space, which is SOUTH6/7/8's `S2A_k20`: both read `results/south7/basis.npz` key `L_mid_full` frozen, at `src/seng3_run.py:64`, `src/axang_run.py:87` and `src/south8_run.py:439-441`. **full** is the 23,485-gene space, SOUTH6's `S1`, paper 1's "space C".

| id | para | what it says | dataset | space | start | floor defined at | drawn at | same space / start | plant gate | gate numbers |
|---|---|---|---|---|---|---|---|---|---|---|
| A1a | 31 | "the four genes are also uncalled" (Southard) | Southard | k20 | GTEx old centroid | `src/south7_run.py:319-324`, reference = pool − draw at `:323`; stored `results/south7/stage_a1_depth_nulls.csv`; attached per factor at `:537-543,555` | own n = 83, 96, 129, 140; 200 draws (`src/south7_run.py:51`) | yes / yes | `src/south8_run.py:157-193`, reference stated `:160-161`, built `:176-178`; `results/south8/d1_plant_curves_all.csv` | `power_25` = 0.925 (83 cells), 0.975 (96), 0.995 (129), 1.000 (140) |
| A1b | 31 | "and move away from young donors in the Southard screen" | Southard | k20 | GTEx old centroid | **CHECKED BY THIS TASK** | — | yes / yes | as A1a | as A1a |
| A2 | 18, 77 | the 20-component count was carried over from the Southard screen; quoted recovery powers 0.990 / 0.965 / 0.930, full gene space 0.115 | Southard | k20/k50/k100 + full | GTEx old centroid | no floor — the quoted numbers **are** size-matched positive-control powers | 106 cells, the Southard screen's median | yes / yes | SOUTH6REC, re-run as a stage-0 check at `src/south8_run.py:266` | 0.990 / 0.965 / 0.930; full gene space 0.115, which FLOORCHK independently confirms |
| A3 | Table 1 row 1; 48, 80, 101 | angles between the Southard cell-division direction and four definitions of aging | Southard | both, reported separately | none (angles) | none — `arccos|cos|`, defined at paragraph 80; no length, no change in distance | — | n/a — paper 2 states angles are never compared across spaces (paras 66, 80, 104) | none | labelled "Descriptive only" (Table 1); no intervals (para 71) |
| A4 | 16 | the Southard and Lu datasets listed as reference datasets | Southard | — | — | no statistic | — | — | — | — |
| B1 | 18, 67 | "keeps 45.1% of the passage gap's length" | Sengstack | ratio across both | PD32 half A | **none exists** — see "Three length ratios with no size-matched floor" | — | — | — | — |
| B2 | 20, 102 | the screen-to-PD32 offset as a multiple of the passage gap (0.65 to 1.05 times) | Sengstack | both | PD32 | **none exists** — `results/seng3/d3_offset.csv` holds a depth-matched bootstrap interval, not a no-effect floor | — | — | — | — |
| B3 | 70 | "14.4 and 14.3 times the donor-age gap" | Sengstack + GTEx | full | — (a position, not a displacement) | **none exists** | — | — | — | — |
| B4 | 40, 47, 104 | shares of a direction removed by projection (12.0%, 41.5%, 66.4% and 64.9%, 13.6% and 20.3%) | Sengstack / GTEx | both | — | none — `1 − √(1 − cos²)` of an angle, defined at paragraph 80 | — | n/a | none | cosine noise attenuates toward zero, so the near-parallel shares are if anything understated |
| B5 | 79 | the mixing positive control moved control samples toward the target by 0.087–0.105, 0.215–0.248 and 0.433–0.484 of the gap | Sengstack | k20 | own-control | `src/seng3_run.py:591`, reference = pool − draw at `:588` | 96 and 157 cells, both screens, 200 draws | yes / yes | the mixture **is** the positive control; `results/seng3/d2_mixture.csv` | as quoted in the paper |
| C1 | 7, 26 | 80 of 340 called closer to early-passage cells; the four factors' −19.2, −13.2, −10.3, −19.9 against thresholds −5.9 to −7.2 | Sengstack | k20 | own-control | `src/seng3_run.py:591-593`, reference = pool − draw at `:588` | own n, every distinct count per screen (`:578`), 200 draws | yes / yes | `results/seng3/d2_plant.csv`; false-call check `src/seng3_run.py:285-287` → `results/seng3/d1_null_validation.csv` | half-gap power 0.940–1.000, quarter-gap 0.565–0.740; false-call 0.040–0.085 against a limit of 0.10. Per-factor `null_median_i` = +0.42, +0.02, +0.46, +0.55 |
| C2 | 26 | 73 called from the late-passage start, including all four | Sengstack | k20 | PD32 half A | as C1, with `t="ii"` | own n | yes / yes | as C1, origin `ii` | half-gap power 0.815–0.995 |
| C3 | 19, 60 | only ZFX and 3 of 340 called in the full gene space; "the near-empty result there carries no information" | Sengstack | full | own-control | `src/seng3_run.py:591`, space `S1` | own n and the four gate depths | yes / yes | `results/seng3/d2_plant.csv`, `S1` rows | floor +14.5 to +24.1, no-effect median +20.6 to +33.7, firing rate 0.000, half-gap power 0.000 in seven of eight cells. This is FLOORCHK's finding, already in the published text |
| C4 | 7, 29, 30 | the four factors' +11.4, +9.0, +5.3, +15.3 against thresholds −1.0 to −1.3; 5 of 340 called | Sengstack | k20 | GTEx old centroid | `src/axang_run.py:516`, reference = pool − draw at `:476` | own n, every distinct count per screen, 200 reps | yes / yes | `src/axang_run.py:486-506` → `results/axang/d3a_plant.csv` | power 0.945–0.990 at 25% of the gap, 1.000 at 50% and 100%; gate QUALIFIED, p25 = 0.945. Dose-0 medians +0.81 to +1.50 against the four factors' +5.3 to +15.3 |
| C5 | 30 | "and all 80 move away from them" (them = young donors) | Sengstack | k20 | GTEx old centroid | **CHECKED BY THIS TASK** | — | yes / yes | as C4 | as C4 |
| C6 | 27 | call rate 0.056 → 0.545 across cell-count deciles; 24 of the 60 perturbations closest in cell count are called | Sengstack | k20 | own-control | as C1 | own n; the band is 96–157 cells, n = 60 | yes / yes | `results/seng3/d5_confounds.json` keys `decile_table`, `depth_band_96_157` | 0.056 → 0.545; 24 of 60 in the band. The band **is** the size control |
| C7 | 37 | Spearman ρ = −0.864 between change in distance and change in cell-cycle score | Sengstack | k20 | own-control | **none** — no threshold of any kind | — | — | none | paper 2 labels it descriptive in the same sentence; paragraph 37 adds the size-controlled version (3 of 4 at or above the 90th percentile among perturbations of similar size) |
| C8 | 39, 40, 105, 106 | called counts after removing a cell-division direction: 82→27, 78→9, and 75 for the gene-list direction | Sengstack | k20 | both | `src/seng4_run.py:458`, keyed `(modality, n, axis, start)` — floors **recomputed after each removal** | own n | yes / yes | `src/seng4_run.py:407` → `results/seng4/d2_plant_residualized.csv` | half-gap power, own-control start: nothing removed 0.895, primary axis 0.010, PD26 axis 0.040, gene-list axis 0.885; late-passage start 0.985 and 0.760 |
| C9 | 54 | of 181 perturbations that score younger on the ruler, 116 are not called closer (120 from the late-passage start) | Sengstack | k20 | both | the call half uses C1's floors; the ruler-score half has no floor by design | own n | yes / yes | as C1 | paper 2 labels it "a diagnostic, not evidence about any factor" |
| C10 | 67 | false-call rate 0.085 from the own-control start, 0.110 from the late-passage start | Sengstack | k20 | both | this **is** the floor validation, `src/seng3_run.py:285-287` | 96, 107, 128, 157; threshold built on one control half, evaluated on the other | yes / yes | self | 0.040–0.085 (own-control), 0.110 (late-passage). Paragraph 67 already discloses that the 0.110 would have failed its pre-set limit and that the choice was written down after the number was known |
| C11 | 19, 68 | half-gap power at least 0.940, quarter-gap at least 0.565 | Sengstack | k20 | own-control | this **is** the plant gate | 96, 107, 128, 157, both screens | yes / yes | self, `results/seng3/d2_plant.csv` | as quoted in the paper |
| C12 | 31 | "any effect that points away from the donor-age direction … must increase the distance to young donors" | Sengstack | k20 | GTEx old centroid | none needed — a geometric identity | — | — | none | with the start at the old centroid, cos(e, Y−O) < 0 forces the distance to increase whatever ‖e‖ is, so no length enters it |

Three findings of the inventory that bear on scope:

- **Paper 2 does not repeat paper 1's corrected claim.** Searching the source for "twice", "2.10", "1,836" and "overshoot" returns only paragraph 13, "apply that readout to the Sengstack data twice". CORR1's correction does not propagate into paper-2 text.
- **Paper 2 already applies FLOORCHK's argument to its own result**, at paragraphs 19 and 60, for the full gene space. A check there would confirm the published text.
- **A1b and C4/C5 are in one space**, so the rule "compute nothing in a space paper 2 does not use" is satisfied for both checked statements.

## A1b — the Southard clause of paragraph 31

> Separately, the four genes are also **uncalled and move away from young donors in the Southard screen**, a different neonatal line; that screen used activation only, so for STAT3 and ZFX it tests the opposite perturbation.

Nothing was drawn. The size-matched baseline already exists at each factor's own cell count, from 200 non-targeting draws with the reference being the non-targeting pool minus the draw (`src/south7_run.py:319-324`), in k20, from the GTEx old centroid to the GTEx young centroid. Per the ruling on Q3 it was read, and reproduced before use.

| gene | cells | change in distance | floor (P5) | no-effect median | no-effect sd | no-effect P95 | d − median | power at 25% of gap | verdict |
|---|---|---|---|---|---|---|---|---|---|
| E2F3 | 96 | +1.562 | -0.107 | +2.520 | 1.787 | +5.643 | -0.958 | 0.975 | CONTRADICTS |
| EZH2 | 140 | +0.481 | -0.694 | +1.294 | 1.415 | +3.967 | -0.813 | 1.000 | CONTRADICTS |
| STAT3 | 83 | +2.320 | +0.468 | +2.695 | 1.743 | +5.932 | -0.375 | 0.925 | CONTRADICTS |
| ZFX | 129 | +2.441 | -0.995 | +1.514 | 1.573 | +4.180 | +0.927 | 0.995 | INSIDE_NOISE |

**3 of the four move away by less than the median no-effect sample of their own size, 1 is inside the no-effect spread, and 0 exceeds its own 95th percentile.** Under a true null of no effect anywhere, about 0.2 of four would exceed the 95th percentile and about 2.0 would fall below the median, so what was observed is what no effect looks like.

At these depths a no-effect draw's change in distance to the GTEx young centroid is **positive**, with a median of +1.3 to +2.7. The positive sign the published clause reports is the baseline's own sign. This is the same defect SOUTH7's criterion amendment was written for, in the opposite direction: there a factor with a positive change in distance was called *approaching* because it moved away less far than noise usually does.

**The "uncalled" half of the sentence stands.** None of the four falls below its own floor, and the test detects a quarter-gap approach at these cell counts with power 0.925 to 1.000, so it is a qualified null and not a blind one. The floor is a 5th percentile, so 95% of true nulls are "uncalled" by construction; it is the power, not the count, that makes the null informative.

**The clause is also uninformative in the full gene space**, where paper 2 does not say it is measured. There no factor's change in distance could be negative at all: FLOORCHK confirmed that all 1,836 Southard displacements are longer than twice the gap. The stored S1 no-effect medians at these four depths are +110.5, +100.6, +82.2 and +78.0 against the four factors' +112.7, +98.2, +83.8 and +79.2 — again all inside the no-effect spread. So the clause fails in k20 and is vacuous in the full gene space.

## C5 — the young-donor clause of paragraph 30

> None of the 80 perturbations called closer to early-passage cells is called closer to young donors, and **all 80 move away from them**.

AXANG computed a floor at every distinct cell count but stored only the 5th percentile. Per the ruling on Q4 its loop was re-run with one change — the median and the 95th percentile are recorded too — at AXANG's own seed 20260914, which is what reproducing the stored 5th percentiles requires. A replicate at the fresh seed 20261030 is reported beside it and gates nothing. Per-size baselines, both seeds: `c5_baseline.csv` (161 rows).

The 80 units are those with `called_delta_i` true in `results/seng3/d4_table.csv`, joined to their young-donor change in distance in `results/axang/d3b_table.csv`. Their cell counts run from 32 to 187, median 80.

| verdict | units | expected under a true null |
|---|---|---|
| SUPPORTS (d ≥ no-effect P95) | 56 | 4.0 |
| INSIDE THE NOISE | 23 | — |
| CONTRADICTS (d < no-effect median) | 1 | 40.0 |

Every one of the 80, at its own cell count:

| perturbation | cells | change in distance | floor (P5) | no-effect median | no-effect P95 | d − median | verdict |
|---|---|---|---|---|---|---|---|
| CRI_SOX10 | 32 | +17.229 | +0.798 | +4.768 | +12.513 | +12.461 | SUPPORTS |
| CRA_DRGX | 33 | +10.182 | +0.316 | +5.533 | +15.048 | +4.650 | INSIDE_NOISE |
| CRI_RXRB | 36 | +17.849 | +0.414 | +4.722 | +12.187 | +13.127 | SUPPORTS |
| CRA_TCF3 | 37 | +10.380 | +0.112 | +5.277 | +12.049 | +5.103 | INSIDE_NOISE |
| CRA_FOXM1 | 38 | +17.518 | +0.445 | +5.185 | +14.192 | +12.334 | SUPPORTS |
| CRI_HOXA13 | 38 | +11.439 | +0.149 | +4.253 | +11.248 | +7.186 | SUPPORTS |
| CRA_PPARGC1A | 40 | +15.713 | +0.405 | +4.852 | +10.763 | +10.861 | SUPPORTS |
| CRA_STAT4 | 43 | +13.561 | -0.215 | +3.588 | +10.037 | +9.973 | SUPPORTS |
| CRI_NPAS1 | 43 | +9.326 | +0.449 | +4.576 | +9.919 | +4.751 | INSIDE_NOISE |
| CRI_ZBTB7A | 44 | +11.192 | -0.104 | +3.596 | +10.104 | +7.596 | SUPPORTS |
| CRI_CEBPA | 49 | +6.787 | +0.125 | +3.375 | +9.013 | +3.412 | INSIDE_NOISE |
| CRI_MITF | 50 | +9.114 | -0.481 | +3.080 | +7.572 | +6.034 | SUPPORTS |
| CRI_PPARG | 51 | +7.149 | -0.093 | +2.845 | +8.911 | +4.304 | INSIDE_NOISE |
| CRA_ZFX | 51 | +5.442 | -0.507 | +3.141 | +8.887 | +2.302 | INSIDE_NOISE |
| CRA_PLAU | 53 | +6.009 | -0.346 | +2.995 | +9.885 | +3.014 | INSIDE_NOISE |
| CRI_GATA1 | 54 | +18.348 | -0.422 | +2.578 | +7.057 | +15.770 | SUPPORTS |
| CRI_DUXA | 55 | +8.169 | -0.397 | +3.248 | +8.285 | +4.921 | INSIDE_NOISE |
| CRI_EGR2 | 55 | +8.210 | -0.397 | +3.248 | +8.285 | +4.962 | INSIDE_NOISE |
| CRA_TFEB | 61 | +6.993 | -0.903 | +2.618 | +8.189 | +4.375 | INSIDE_NOISE |
| CRA_POU3F2 | 61 | +6.573 | -0.903 | +2.618 | +8.189 | +3.956 | INSIDE_NOISE |
| CRI_THRA | 63 | +9.406 | -0.807 | +2.053 | +6.604 | +7.353 | SUPPORTS |
| CRA_PPARG | 64 | +5.340 | -0.565 | +2.441 | +7.122 | +2.900 | INSIDE_NOISE |
| CRI_FOXP1 | 64 | +9.509 | -0.582 | +2.611 | +7.247 | +6.897 | SUPPORTS |
| CRI_ATG5 | 65 | +12.183 | -0.363 | +2.574 | +6.827 | +9.609 | SUPPORTS |
| CRA_MEF2C | 67 | +10.036 | -0.981 | +2.100 | +6.877 | +7.937 | SUPPORTS |
| CRI_CTCFL | 68 | +16.374 | -0.719 | +2.372 | +7.086 | +14.002 | SUPPORTS |
| CRI_GATA2 | 70 | +7.576 | -0.773 | +2.218 | +7.384 | +5.358 | SUPPORTS |
| CRI_EGR1 | 70 | +16.588 | -0.773 | +2.218 | +7.384 | +14.369 | SUPPORTS |
| CRA_SALL4 | 71 | +7.037 | -0.660 | +2.131 | +7.426 | +4.906 | INSIDE_NOISE |
| CRI_ARID5A | 72 | +10.875 | -0.588 | +1.956 | +6.293 | +8.919 | SUPPORTS |
| CRI_E2F5 | 72 | +9.709 | -0.588 | +1.956 | +6.293 | +7.753 | SUPPORTS |
| CRI_ATF4 | 73 | +14.025 | -0.739 | +2.521 | +7.166 | +11.503 | SUPPORTS |
| CRI_XRCC5 | 73 | +8.596 | -0.739 | +2.521 | +7.166 | +6.074 | SUPPORTS |
| CRI_SOX5 | 73 | +9.809 | -0.739 | +2.521 | +7.166 | +7.288 | SUPPORTS |
| CRI_SOX7 | 73 | +7.550 | -0.739 | +2.521 | +7.166 | +5.029 | SUPPORTS |
| CRI_APOD | 74 | +7.217 | -0.910 | +1.695 | +6.824 | +5.521 | SUPPORTS |
| CRI_MRRF | 75 | +7.804 | -0.452 | +2.081 | +6.980 | +5.722 | SUPPORTS |
| CRI_TFEB | 77 | +10.637 | -0.877 | +1.861 | +6.844 | +8.776 | SUPPORTS |
| CRA_VAX1 | 79 | +11.143 | -0.871 | +2.128 | +5.970 | +9.015 | SUPPORTS |
| CRI_FOXA3 | 79 | +8.282 | -0.562 | +1.994 | +6.051 | +6.287 | SUPPORTS |
| CRI_PARP1 | 81 | +7.471 | -0.869 | +1.994 | +6.642 | +5.476 | SUPPORTS |
| CRA_FOXD2 | 82 | +1.620 | -1.174 | +1.430 | +6.625 | +0.191 | INSIDE_NOISE |
| CRI_GSX2 | 83 | +6.436 | -0.938 | +2.066 | +5.925 | +4.370 | SUPPORTS |
| CRA_GATA3 | 85 | +6.788 | -0.937 | +1.795 | +6.326 | +4.992 | SUPPORTS |
| CRI_ONECUT1 | 86 | +7.356 | -0.974 | +1.592 | +5.568 | +5.764 | SUPPORTS |
| CRI_SIM1 | 87 | +11.921 | -0.920 | +1.530 | +4.969 | +10.391 | SUPPORTS |
| CRI_TFAM | 88 | +8.831 | -1.101 | +1.609 | +6.200 | +7.222 | SUPPORTS |
| CRI_MAZ | 88 | +11.071 | -1.101 | +1.609 | +6.200 | +9.462 | SUPPORTS |
| CRA_IRF4 | 89 | +3.942 | -1.439 | +1.421 | +6.268 | +2.521 | INSIDE_NOISE |
| CRA_DLX6 | 89 | +13.059 | -1.439 | +1.421 | +6.268 | +11.637 | SUPPORTS |
| CRI_PGR | 89 | +11.540 | -0.477 | +1.838 | +5.506 | +9.702 | SUPPORTS |
| CRI_DLX4 | 89 | +5.283 | -0.477 | +1.838 | +5.506 | +3.445 | INSIDE_NOISE |
| CRI_FERD3L | 89 | +13.568 | -0.477 | +1.838 | +5.506 | +11.730 | SUPPORTS |
| CRI_NFE2L2 | 91 | +9.459 | -0.842 | +1.564 | +5.790 | +7.895 | SUPPORTS |
| CRI_HMGB1 | 92 | +7.359 | -0.931 | +1.412 | +4.819 | +5.947 | SUPPORTS |
| CRI_STAT3 | 96 | +5.327 | -1.024 | +1.502 | +5.430 | +3.825 | INSIDE_NOISE |
| CRI_VEZF1 | 97 | +8.646 | -1.083 | +1.504 | +5.532 | +7.142 | SUPPORTS |
| CRI_BSX | 98 | +7.408 | -1.252 | +1.588 | +6.026 | +5.820 | SUPPORTS |
| CRI_LCOR | 99 | +8.361 | -0.840 | +1.517 | +5.334 | +6.844 | SUPPORTS |
| CRA_HOXD3 | 100 | +5.222 | -0.667 | +1.704 | +5.469 | +3.517 | INSIDE_NOISE |
| CRI_FOXB2 | 106 | +6.443 | -1.189 | +1.278 | +5.400 | +5.166 | SUPPORTS |
| CRA_EZH2 | 107 | +9.022 | -1.223 | +1.242 | +5.156 | +7.779 | SUPPORTS |
| CRA_ARX | 112 | +5.026 | -1.434 | +1.141 | +5.887 | +3.885 | INSIDE_NOISE |
| CRI_SOX9 | 113 | +7.734 | -1.102 | +1.270 | +4.158 | +6.464 | SUPPORTS |
| CRA_JUND | 118 | +4.773 | -1.392 | +0.986 | +4.544 | +3.786 | SUPPORTS |
| CRA_ETS2 | 119 | +3.477 | -1.190 | +1.058 | +4.753 | +2.420 | INSIDE_NOISE |
| CRA_TEAD3 | 119 | +4.752 | -1.190 | +1.058 | +4.753 | +3.694 | INSIDE_NOISE |
| CRA_TFCP2 | 120 | +0.927 | -1.382 | +1.192 | +4.946 | -0.265 | CONTRADICTS |
| CRA_HSF2 | 121 | +3.955 | -1.108 | +1.107 | +4.853 | +2.848 | INSIDE_NOISE |
| CRI_ZFHX3 | 121 | +14.052 | -1.007 | +1.195 | +4.297 | +12.857 | SUPPORTS |
| CRI_SOX2 | 125 | +9.714 | -1.317 | +1.359 | +5.602 | +8.355 | SUPPORTS |
| CRI_FOXP2 | 125 | +8.775 | -1.317 | +1.359 | +5.602 | +7.417 | SUPPORTS |
| CRA_E2F3 | 128 | +11.448 | -1.346 | +1.292 | +5.001 | +10.156 | SUPPORTS |
| CRA_TLX2 | 140 | +4.118 | -1.441 | +1.011 | +4.712 | +3.107 | INSIDE_NOISE |
| CRA_SOX15 | 140 | +9.773 | -1.441 | +1.011 | +4.712 | +8.762 | SUPPORTS |
| CRI_CDX1 | 146 | +5.023 | -1.045 | +1.073 | +4.155 | +3.950 | SUPPORTS |
| CRA_FOSL1 | 149 | +7.566 | -0.975 | +1.149 | +4.374 | +6.417 | SUPPORTS |
| CRI_EOMES | 155 | +3.768 | -1.312 | +1.004 | +3.724 | +2.764 | SUPPORTS |
| CRI_ZFX | 157 | +15.336 | -1.202 | +0.851 | +3.891 | +14.486 | SUPPORTS |
| CRA_MYC | 187 | +7.812 | -1.176 | +0.725 | +3.369 | +7.087 | SUPPORTS |

**1 of the 80 move away by less than the median no-effect sample of their own size:** `CRA_TFCP2`. The published "all 80 move away from them" is true of the measured sign for all 80, and established against a size-matched baseline for 56.

**The units that fail are the shallow ones.** By verdict, their cell counts are:

| verdict | units | cells: min / median / max | median change in distance |
|---|---|---|---|
| Supports | 56 | 32 / 82 / 187 | +9.61 |
| Inside Noise | 23 | 33 / 64 / 140 | +5.44 |
| Contradicts | 1 | 120 / 120 / 120 | +0.93 |

So the failure of the universal claim is a depth effect, not a sign that those perturbations behave differently: the no-effect baseline rises as cell count falls, and the units it swallows are the small ones.

**How precise the baseline is.** It rests on 200 draws, so it is itself a random quantity. Between the primary seed and the fresh-seed replicate the no-effect median moves by 0.203 at the median cell count and up to 1.057 at worst, and the 5th percentile by 0.227 and up to 0.761. Units whose `d − median` is smaller than that could change verdict at another seed, so the three counts are approximate at the margin. The per-unit `d − median` column is printed above so the margin is visible for every unit, and both seeds' baselines are in `c5_baseline.csv`.

**Circularity, disclosed here because this is where the result appears.** The 80 are selected by the passage test, which uses the same screen control cells that define zero for the young-donor statistic and that every no-effect draw is taken from. The selection and the quantity tested are not independent.

## Corrected wording (step 3)

The paper-2 source was not edited and no correction notice was written, as the task requires. What one would have to say:

**Paragraph 31, the Southard clause.** The sentence as published joins a size-matched claim to an unmatched one with "and". "Uncalled" holds. "Move away from young donors" does not: it compares a change in distance with zero when the matched no-effect baseline at those cell counts is positive, and three of the four factors fall below that baseline. A corrected clause could read:

> Separately, the four genes are also uncalled in the Southard screen, a different neonatal
> line, where the test detects a quarter-gap approach toward young donors at their cell
> counts with power 0.925 to 1.000. Their changes in distance to young donors are positive,
> but so are those of non-targeting cells with no effect at the same cell counts: the
> no-effect median runs from +1.3 to +2.7 against the four factors' +0.5 to +2.4. Three of
> the four move away less far than the median no-effect sample of their own size, and none
> exceeds its 95th percentile, so that screen shows no approach and does not show a
> departure either. That screen used activation only, so for STAT3 and ZFX it tests the
> opposite perturbation.

**Paragraph 30, the young-donor clause.** "None of the 80 … is called closer to young donors" holds; it rests on the same size-matched thresholds as the rest of the young-donor test, whose gate fired QUALIFIED at p25 = 0.945. The clause "and all 80 move away from them" overstates what a size-matched comparison supports. A corrected clause could read:

> None of the 80 perturbations called closer to early-passage cells is called closer to young donors, and all 80 have a positive change in distance to them. Against no-effect samples of the same size, 56 of the 80 move away further than the 95th percentile, 23 lie inside the no-effect spread and 1 move away less far than the median no-effect sample, so the departure is established for 56 of them and not for the rest.

Neither correction touches paper 2's conclusion. That conclusion is that the two definitions of young give opposite answers, and it rests on the *calls*: 80 of 340 called closer to early-passage cells against 5 of 340 called closer to young donors, with 17 expected by chance in each case. Both counts are gated by size-matched thresholds and neither is affected. What the corrections remove is a stronger directional claim laid on top of the calls.

## Three length ratios with no size-matched floor (B1, B2, B3)

Per the ruling these are recorded, not checked. For each, the direction of the noise bias is stated, and in each case it **strengthens rather than weakens paper 2's own caveat**, which is why none was checked.

- **B1, paragraphs 18 and 67**, "keeps 45.1% of the passage gap's length". A ratio of one vector's norm in two spaces. Sampling noise inflates the full-gene norm far more than the 20-component one, so 45.1% understates the true retention. Paper 2 uses it to argue that compression discards more noise than signal; the bias works against that argument, not for it.
- **B2, paragraphs 20 and 102**, the screen-to-PD32 offset at 0.65 to 1.05 times the passage gap. Numerator and denominator are both noisy pseudobulk differences, and the denominator's PD14 half A is the smallest group in the comparison. `results/seng3/d3_offset.csv` holds a depth-matched bootstrap interval, not a no-effect floor, and FLOORCHK's construction does not transfer: the two populations compared sit in different sample groups, so there is no single no-effect pool to draw both ends from. Paragraph 65 calls the offset "an upper bound on the technical difference", so noise inflation strengthens the disclosure.
- **B3, paragraph 70**, "14.4 and 14.3 times the donor-age gap". A distance from a 4,348-cell pseudobulk to two bulk centroids, not a displacement. Noise inflates it, so correction could only reduce it, and the statement's role is to warn that the GTEx-anchored test starts far outside the donor-age span. For scale, FLOORCHK's stored full-gene floor at the largest Southard group size (341 cells) is a median of 2.06 gap lengths and falls with cell count, while this group has 4,348 cells — but that is a different dataset with a different gene overlap, so the number does not transfer and none was computed.

## Step 4 — paper 1, paragraph 142

Separate from everything above, and computing nothing. `CORRECTIONS_PAPER1.md` identifies a bullet in paragraph 142 of paper 1 as carrying no information:

> • In space C both perturbations move away from the young centroid (+106.6 and +113.8, q = 1.0), and in space B both end farther from it (+32.9 and +35.1). They have 84 and 82 cells, against a median of 106.

Replacement wording, for folding into the paper-1 revision. No file was edited.

> • Neither of the other spaces separates these two from the rest of the screen. In space C
> both end farther from the young centroid (+106.6 and +113.8), as do all 1,836
> perturbations: there every displacement is longer than twice the 32-unit gap, so no change
> in distance could have come out negative, and the median non-targeting floor is 2.9 gap
> lengths. In space B both also end farther (+32.9 and +35.1), again with none of the 1,836
> moving closer, and both sit just inside their own cell-count-matched floors (+33.3 and
> +35.7), as 1,088 of the 1,836 do. They have 84 and 82 cells, against a median of 106.

Every number is read from stored files:

| number | source |
|---|---|
| +106.6, +113.8; +32.9, +35.1; 84 and 82 cells | `results/coordspace/stage3_south6_factors.csv` (106.609, 113.750; 32.886, 35.053) — the paper's own values |
| none of the 1,836 moves closer, in space C or space B | `results/coordspace/stage3_south6_summary.csv`, `n_delta_negative` = 0 for both |
| every space-C displacement longer than twice the 32-unit gap | FLOORCHK reproduction R1: smallest ratio 2.0966, 1,836 of 1,836; gap 32.284 |
| median non-targeting floor 2.9 gap lengths in space C | `median_floor_over_gap` = 2.9415 |
| space-B floors +33.3 and +35.7 | `floor_p05` = 33.344, 35.705 |
| 1,088 of 1,836 below floor in space B | `n_below_floor` = 1088 |

Two judgement calls in that wording, both reversible by the author:

1. **"q = 1.0" is dropped.** That q comes from the permutation null SOUTH8 found **anti-conservative in space C** — it flags 0.09 to 0.14 of true-null draws at the 0.05 level (`results/south8/d2_calibration.csv`, S1 rows; verdict in `d2_verdicts.csv`). Quoting it as though it were a gate is the thing the bullet should stop doing, and in context it adds nothing the +106.6 does not.
2. **The space-B floor comparison is added.** The original bullet's space-B clause has the same weakness as its space-C clause: "both end farther" is true of all 1,836. But space B is not vacuous by geometry, and there the two labels are in fact *below* their own floors, which is the informative fact. Leaving it out would replace one uninformative clause with another.

## Deviations and decisions not in the task (disclosed)

1. **The session knew both outcomes before the flag was written.** The task's step 1 required an inventory, and the user's second message required it to be reported in chat before the flag. That inventory read the stored baselines and the stored changes in distance, so it established that A1b would fail and that at least one of C5's 80 units would fall below its baseline. Flag Part 0 records exactly what was known, with the numbers. The four conclusion wordings were fixed by the user in the task message before any of it was read, which is what protects the result.
2. **Four questions went to the user before the flag, and two of the session's own groupings were wrong.** The inventory first reported C5 as needing no check and placed C7 under a floor construction it does not have. Both were corrected in the second inventory pass, before the flag: C5 moved to the checked list and C7's verdict was re-grounded on having no threshold at all. The rulings Q1 to Q4 are reproduced verbatim in the flag.
3. **The FLOORCHK quadrature correction is used nowhere in this task.** Per the ruling on Q2, a change-in-distance statistic already measured against draws of the same size has its noise subtracted by construction, so a quadrature step would double-count it. Both checked statements are of that kind. This is a substitution of statistic relative to the task's original wording, authorised in writing.
4. **A tension in the ruling on Q4 was resolved by the session, not by the user.** The ruling asked for AXANG's loop "unchanged except for recording the median", for "fresh seeds", and for the stored 5th percentiles to reproduce "from the same re-run". Those cannot all hold: a 5th percentile of 200 draws is random, so a fresh-seed run cannot reproduce the stored value exactly. The resolution, fixed in flag D7 before the run: the primary re-run keeps AXANG's own seed, because the seed is an argument of the loop and "unchanged" therefore includes it, and a fresh-seed replicate is reported beside it and gates nothing. Both are in `c5_baseline.csv`, so nothing is hidden by the choice. This is the one place the session resolved an ambiguity instead of stopping; it affects no conclusion, because both runs are reported.
5. **AXANG's `d3()` was not called.** It writes into `results/axang/`, which CLAUDE.md forbids. Its prologue and floor loop are repeated in `src/corr2chk_run.py` with the helpers imported unchanged from `axang_run`, `seng3_run` and `seng2_common`. R5 and R6 are the proof the repetition is faithful: every stored 5th percentile and all 340 stored changes in distance reproduce. Importing `axang_run` and `seng3_run` calls `mkdir(exist_ok=True)` on their results directories, which already exist; no file in them was written.
6. **A1b draws nothing.** Per the ruling on Q3 the stored baseline was read, because it exists at exactly the four cell counts with the required construction. No seed was needed for it, and the task's "seeds: fresh" therefore applies only to C5's stability replicate.
7. **The 17 NOT-AFFECTED verdicts are inventory judgements, not computed results.** No number of this task gates them. Their evidence is the table above: the file and line of each floor, the group size it was drawn at, whether it shares the statement's space and starting profile, and the file, line and numbers of its plant gate.
8. **An arithmetic slip in the flag, disclosed and not edited.** Flag Part 3 says "21 statements were found" and then gives three counts — 2 checked, 17 NOT-AFFECTED, 3 with no size-matched floor — each labelled "of 21". Those sum to 22, because paragraph 31's sentence is one statement among the 21 and two clause-level entries (A1a and A1b) once the ruling on Q1 splits it. The correct statement is 21 statements and 22 clause-level verdicts. The flag is the pre-registration and is reproduced above exactly as written, slip included; no count used in any result depends on it.
9. **Where things were written.** `results/corr2chk/` only, plus the two tracked source files and the three deliverables. Nothing under `results/floorchk/`, `results/corr1/`, `results/south6/`, `results/south7/`, `results/south8/`, `results/axang/`, `results/seng3/`, `results/seng4/` or `results/coordspace/` was modified, no existing `FINDINGS_*.md` or `CORRECTIONS_PAPER1.md` was touched, and no paper source was edited. The pre-existing uncommitted change to `CITATION.cff` and the untracked `results/paper2_figs_plain/` are not this task's work and are left out of its commit.

## Honesty obligations

- **Chance baseline beside every count.** The 95th percentile is exceeded by 5% of no-effect draws by construction and the median by 50%, so under a true null about 0.2 of A1b's four and 4.0 of C5's 80 would read SUPPORTS, and about 2.0 and 40.0 would read CONTRADICTS. Both are printed beside every count above.
- **Recovery power beside every null.** A1b's per-factor power at 25% of the gap at its own depth is in the result table (0.925 to 1.000). C5's gate is AXANG's D3a, QUALIFIED at p25 = 0.945 and p50 = 1.000. A1a's "uncalled" verdict is never reported as a bare count.
- **Cell count beside every verdict.** Both result tables carry each unit's own cell count. For C5 the units run from 32 to 187 cells, and the no-effect median falls with cell count, so the baseline differs across units by more than the quantity being tested in some cases. No unit is presented as a candidate for anything.
- **Circularity, disclosed and not refitted away.** For both statements the non-targeting pool of the relevant screen defines zero for every change in distance and is also the pool every no-effect draw is taken from. For C5 the 80 units are selected by a statistic that uses those same control cells.
- **Precision of the baseline.** The medians rest on 200 draws. The fresh-seed replicate in `c5_baseline.csv` shows the Monte Carlo spread, and verdicts for units close to the boundary are not stable at that precision.
- **What this does not cover.** Cell-sampling noise only: not guide-to-guide, batch or biological variation. Southard and Sengstack are each one neonatal fibroblast line with no aged donor, so nothing here transfers to aged human cells. Nothing here tests any factor's phenotype. The GTEx-anchored construction transplants a displacement measured in single cells onto a bulk post-mortem origin, as both papers state. Both checked clauses are about what the measurement can resolve, not about what any factor does.
- **No conclusion about any named gene or factor** beyond what the existing papers already state. A1b's units carry gene names because the published sentence names them; the result is about the sentence.

`python -I src/corr2chk_run.py verify` re-derives every reading line from the tracked CSVs alone.

# FINDINGS_SOUTH8 — three diagnostics on the SOUTH7 result, plus the all-1,836 reading

Status: Stage 0 passed (46 checks). D1, D2, D3 and D4 all ran; no STOP fired.

Headline: the permutation null that produced SOUTH7's zero is **not calibrated** in S2A_k20, and the excess of 27 factors below the cell-count-matched floor is **not accounted for by read depth**. FINDINGS_SOUTH7.md is unchanged and remains the record of the primary run.

## Pre-registration (verbatim, `results/south8/PREREG.flag`)

```
TASK: SOUTH8 — three diagnostics on the SOUTH7 result, plus the all-1,836
reading. Written before any SOUTH8 number exists.

Nature: computational. No new data. Reads results/south7/ (including the
completed results/south7/stage_a2_perm_rest.npz), results/south6/,
results/south6rec/ and the Southard singlets file already on disk. Writes only
results/south8/ and FINDINGS_SOUTH8.md. FINDINGS_SOUTH7.md is not regenerated
and no file in results/south7/ is modified; it stays as the record of the
primary run. Nothing is committed or pushed.

--- WHY THIS SESSION EXISTS ---

SOUTH7's gate read zero in S2A_k20 on the published 304-factor set. Three things
about that zero need testing before anything is built on it.

1. The zero has no size attached. SOUTH7 reports that no factor is called, and
   that a planted 25%-of-gap approach would be detected with power 0.980. It
   never reports how large the largest apparent approach actually was. Without
   that, "none" cannot be turned into a bound.

2. Two nulls disagree. In S1, SOUTH7's permutation test passes 51 of 304 factors
   at q <= 0.05 while the depth-matched empirical floor passes 3. In S2A_k20 the
   permutation test passes 0 while 27 fall below the empirical floor. Both
   directions of disagreement cannot be right, and the permutation null — a
   uniform shuffle of gene coordinates — destroys the correlation structure that
   real expression vectors carry, which is a known way for such a null to be
   mis-specified. The permutation test is what produced the zero, so its
   calibration is load-bearing.

3. The depth-matched null is matched on cell count. Sampling noise is set by
   total reads, which is cells multiplied by reads per cell. If the 27 factors
   below the floor simply have richer cells than the control draws of the same
   cell count, the excess is a depth artifact of the same kind that defeated the
   original analysis, in a different variable.

--- STAGE 0: REPRODUCE FIRST (gate, before any new number) ---
Re-derive from the stored SOUTH7 files, with this task's code:
 (a) the Stage A3 row for primary_304 in each of the four spaces: called, below
     depth floor, q <= 0.05, SOUTH6 rule, tied at min p, smallest q;
 (b) that the smallest q in S2A_k20 is 0.0578 and is held by AFF1, EGR3, FOXN3;
 (c) the MYC row in S1 (cells 99, delta 94.0537, floor 94.3758, q 0.0018);
 (d) the S2A_k20 basis gap retention 0.8651 and ||O1-Y1|| 27.9292;
 (e) the R6 recovery rows at 106 cells in all four spaces.
Pass: counts exact, every numeric item within 1e-4 * max(1, |stored|).
Fail -> STOP `prior_results_not_reproduced`, naming every mismatch.

--- CRITERION AMENDMENT (declared now, before any SOUTH8 number) ---
SOUTH7's criterion called MYC APPROACHING in S1 although MYC's change in
distance is positive: it moved away from the young target, only less far than
sampling noise usually does. Comparing each factor with its own null rather than
with zero is the preprint's remedy and is kept, but it is not sufficient on its
own. From here:

  APPROACHING requires BOTH (i) change in distance below the 5th percentile of
  the factor's own depth-matched null, AND (ii) change in distance < 0.

SOUTH7's S1 reading is marked superseded, not edited; its text stays in
results/south7/. Both counts are reported in every space so the effect of the
amendment is visible.

--- D1. The size of the zero (descriptive; no decision rule) ---
For each space, on the 304-factor set and on all 1,836:
 (a) the five most negative changes in distance, with cell count, cell-count
     percentile, depth-matched floor, p, q, and the high_depth_candidate flag;
 (b) each one expressed as a fraction of the young-old gap, by linear
     interpolation against the R6 plant dose-response (0, 25, 50, 100%) rebuilt
     at that factor's own cell count. State the interpolation method used and
     report the value as an interval between the two bracketing doses where it
     falls outside the measured range;
 (c) the resulting sentence, filled in: "The largest apparent approach anywhere
     in the screen corresponds to at most <X>% of the young-old gap, measured
     with a statistic that detects a 25% approach at this depth with power
     0.980."
This is descriptive and carries no decision rule. It is not an effect size in
aged cells: the screen has no aged donor.

--- D2. Is the permutation null calibrated? (gated) ---
The depth-matched null draws contain, by construction, no effect. Run SOUTH7's
permutation test on them exactly as on a real factor.
 - For each of the 200 non-targeting draws at cell counts 50, 106, 180 and 341,
   in each of the four spaces, permute that draw's displacement vector as
   SOUTH7's Stage A2 does and compute its p.
 - 10,000 permutations per draw is sufficient and is permitted here: the
   question is the rate at 0.05, not a tail at 1e-5. State that it was used.
 - Report, per space and per cell count, the fraction of true-null draws with
   p <= 0.05 and with p <= 0.001.
Gate, wording fixed now:
 - CALIBRATED: "the permutation null flags <f> of true-null draws at the 0.05
   level in <space>, consistent with the 0.05 it should flag."  Condition:
   0.03 <= f <= 0.08 at every cell count tested in that space.
 - ANTI-CONSERVATIVE (f > 0.08 anywhere): "the permutation null flags <f> of
   true-null draws at the 0.05 level in <space>, so it cannot be read as
   evidence that a factor's direction is special."
 - CONSERVATIVE (f < 0.03 anywhere): "the permutation null flags <f> of
   true-null draws at the 0.05 level in <space>, so it cannot be used to
   withhold a call."
Consequence, fixed now: if the null is not CALIBRATED in S2A_k20, the SOUTH7
primary call is re-read using the depth-matched floor and the amended criterion
alone, without the permutation condition. That re-reading is reported as the
primary reading of SOUTH8, and SOUTH7's is marked superseded, not edited. Both
are shown side by side.

--- D3. Is the null matched on the right quantity? (gated) ---
 (a) For every factor, total UMI across its cells, from the counts already used
     in Stage 0 (b). Report Spearman of total UMI against change in distance,
     and against cell count, in each space, on the 304 set and on all 1,836.
 (b) Rebuild the null matched on total UMI rather than on cell count: for each
     factor, draw 200 non-targeting pseudobulks whose total UMI is within 2% of
     that factor's, by the SOUTH6REC scheme (disjoint blocks within a
     permutation, reference = complement). Report how many factors can be
     matched within tolerance, the realised tolerance distribution, and the
     cells-per-draw actually used.
 (c) Recount how many of the 304, and of the 1,836, fall below the UMI-matched
     5th percentile in each space, and re-run the amended criterion against it.
Gate, wording fixed now, for S2A_k20 on the 304 set. 15.2 is 5% of 304, and the
95% binomial interval for 5% of 304 is 8 to 23 inclusive:
 - EXPLAINED (k <= 23): "the excess of 27 factors below the cell-count-matched
   floor is accounted for by read depth; matched on total reads, <k> of 304 fall
   below, which is consistent with chance."
 - SURVIVES (k > 23): "the excess survives read-depth matching: <k> of 304 fall
   below a floor matched on total reads, against 15.2 expected by chance."
Both counts are reported regardless of which fires.

--- D4. The all-1,836 reading (report-only) ---
From results/south7/stage_a2_perm_rest.npz, compute q over m = 1,836 and report
SOUTH7's Stage A3 table for the full screen in all four spaces: called under the
amended criterion, below depth floor, q <= 0.05, SOUTH6 rule, tied at min p,
smallest q. Add the two factors of SOUTH6's recurring set that fell outside the
published primary set, CERS3 and REST, to the R3 table. Report whether any
factor called under the amended criterion in any space carries the
high_depth_candidate flag.

--- WHAT DOES NOT COUNT ---
- choosing a space, a criterion, a null or a factor set after seeing a number;
- reading D1's bound as a biological effect size in aged cells; the screen uses
  neonatal foreskin fibroblasts with no aged donor and no matched young target;
- any result from an O2 pair;
- the SOUTH6 recurring 11, or any earlier shortlist, as a prior;
- a factor called under the amended criterion being described as a candidate for
  the pilot panel in this findings file.

--- STOP RULES ---
On any STOP: write the findings with the key and exit. Never compute after a
STOP; never resume a stopped run.

--- OUTPUT ---
FINDINGS_SOUTH8.md: this text verbatim; the Stage 0 reproduction table; D1, D2,
D3 with the gate reading that fired in each; D4; limitations. Tables in
results/south8/. Code src/south8_run.py. Hand over files, not summaries. Do not
regenerate FINDINGS_SOUTH7.md. Do not commit or push.
```

## 1. Stage 0: SOUTH7's reported numbers re-derived (gate)

46 checks, 46 gating, 0 gating failures. Counts exact; numeric items within the stated tolerance. Full table: `results/south8/stage0_reproduction.csv`.

| block | items | max abs diff | mismatches | passed |
|---|---|---|---|---|
| (a) SOUTH7 A3 primary_304 S2A_k20 | 6 | 4.86e-17 | 0 | yes |
| (a) SOUTH7 A3 primary_304 S2A_k50 | 6 | 0 | 0 | yes |
| (a) SOUTH7 A3 primary_304 S2A_k100 | 6 | 0 | 0 | yes |
| (a) SOUTH7 A3 primary_304 S1 | 6 | 3.31e-17 | 0 | yes |
| (b) smallest q in S2A_k20 | 2 | 4.06e-05 | 0 | yes |
| (c) MYC in S1 | 4 | 3.82e-05 | 0 | yes |
| (d) S2A_k20 basis | 4 | 4.55e-05 | 0 | yes |
| (e) R6 at 106 cells S2A_k20 | 3 | 4.44e-16 | 0 | yes |
| (e) R6 at 106 cells S2A_k50 | 3 | 0 | 0 | yes |
| (e) R6 at 106 cells S2A_k100 | 3 | 2.78e-17 | 0 | yes |
| (e) R6 at 106 cells S1 | 3 | 0 | 0 | yes |

Nothing in `results/south7/`, `results/south6/` or `results/south6rec/` was written. SOUTH7's Stage A3 was never called; its per-factor tables were rebuilt here from `stage_a1_observed.csv`, `stage_a1_depth_nulls.csv` and the two `stage_a2_perm_*.npz` files, and compared against its stored `stage_a3_calls.csv`.

## 2. The criterion amendment

Declared in the pre-registration before any SOUTH8 number: APPROACHING now requires both (i) a change in distance below the 5th percentile of the factor's own depth-matched null and (ii) a change in distance below zero. SOUTH7's S1 reading, which called MYC approaching on (i) alone although MYC moved away from the young target, is **superseded, not edited**: its text stands in `FINDINGS_SOUTH7.md` and `results/south7/`.

Both counts, every space, so the effect of the amendment is visible:

| factor set | space | m | SOUTH7 criterion (floor+q) | amended (floor+delta<0) | amended+q | below depth floor | delta < 0 | q <= 0.05 | smallest q |
|---|---|---|---|---|---|---|---|---|---|
| primary_304 | S2A_k20 | 304 | 0 | 24 | 0 | 27 | 35 | 0 | 0.0578 |
| primary_304 | S2A_k50 | 304 | 0 | 4 | 0 | 23 | 4 | 0 | 0.3040 |
| primary_304 | S2A_k100 | 304 | 0 | 1 | 0 | 20 | 1 | 0 | 0.9576 |
| primary_304 | S1 | 304 | 1 | 0 | 0 | 3 | 0 | 51 | 0.0003 |
| all_1836 | S2A_k20 | 1836 | 4 | 125 | 4 | 153 | 214 | 4 | 0.0184 |
| all_1836 | S2A_k50 | 1836 | 0 | 19 | 0 | 130 | 21 | 0 | 0.3213 |
| all_1836 | S2A_k100 | 1836 | 0 | 6 | 0 | 118 | 7 | 0 | 0.5692 |
| all_1836 | S1 | 1836 | 6 | 0 | 0 | 48 | 0 | 210 | 0.0007 |

The amendment removes MYC from the S1 call and adds nothing anywhere: in every space the amended count is the number of factors that are both below their own floor and actually moving toward the target.

## 3. D1: the size of the zero (descriptive, no decision rule)

Dose equivalents come from the planted dose-response rebuilt at each factor's own cell count (doses 0, 25, 50 and 100% of the measured gap, 200 draws each, reference = complement of the draw). **Interpolation method:** the curve of median planted change in distance is monotone decreasing in dose, so a factor's observed change is placed by piecewise-linear interpolation between the two bracketing doses; a value outside the measured range is reported as a bound, not a point.

### primary_304, S2A_k20

| factor | cells | cell pctl in set | high_depth_candidate | delta | depth floor | below floor | p | q | % of gap | bracket | power at 25% at this depth |
|---|---|---|---|---|---|---|---|---|---|---|---|
| EGR3 | 71 | 10.8553 | no | -2.3610 | 0.9528 | yes | 0.0004 | 0.0578 | 27.8152 | between 25% and 50% | 0.8850 |
| FOXN3 | 128 | 61.1842 | no | -2.1206 | -0.8317 | yes | 0.0002 | 0.0578 | 13.9859 | between 0% and 25% | 0.9900 |
| SOX8 | 114 | 50.0000 | no | -1.8734 | -0.7587 | yes | 0.0022 | 0.1313 | 14.6811 | between 0% and 25% | 0.9850 |
| AFF1 | 177 | 85.8553 | no | -1.8316 | -1.2665 | yes | 0.0006 | 0.0578 | 10.5669 | between 0% and 25% | 1.0000 |
| BNC1 | 152 | 74.3421 | no | -1.7147 | -1.1506 | yes | 0.0013 | 0.0996 | 11.9193 | between 0% and 25% | 1.0000 |

### primary_304, S2A_k50

| factor | cells | cell pctl in set | high_depth_candidate | delta | depth floor | below floor | p | q | % of gap | bracket | power at 25% at this depth |
|---|---|---|---|---|---|---|---|---|---|---|---|
| AFF1 | 177 | 85.8553 | no | -1.2575 | -0.3630 | yes | 0.0010 | 0.3040 | 11.6522 | between 0% and 25% | 1.0000 |
| SOX8 | 114 | 50.0000 | no | -0.5731 | 1.1536 | yes | 0.0093 | 1.0000 | 16.9214 | between 0% and 25% | 0.9700 |
| BNC1 | 152 | 74.3421 | no | -0.4569 | 0.0059 | yes | 0.0151 | 1.0000 | 10.7691 | between 0% and 25% | 1.0000 |
| AFF4 | 199 | 91.4474 | yes | -0.3399 | -0.0374 | yes | 0.0338 | 1.0000 | 7.5683 | between 0% and 25% | 1.0000 |
| ARID3A | 242 | 97.6974 | yes | 0.0564 | -0.2958 | no | 0.2183 | 1.0000 | 3.7000 | between 0% and 25% | 1.0000 |

### primary_304, S2A_k100

| factor | cells | cell pctl in set | high_depth_candidate | delta | depth floor | below floor | p | q | % of gap | bracket | power at 25% at this depth |
|---|---|---|---|---|---|---|---|---|---|---|---|
| AFF1 | 177 | 85.8553 | no | -0.5266 | 0.1430 | yes | 0.0031 | 0.9576 | 11.1655 | between 0% and 25% | 1.0000 |
| AFF4 | 199 | 91.4474 | yes | 0.2203 | 0.4985 | yes | 0.0500 | 1.0000 | 7.6113 | between 0% and 25% | 1.0000 |
| BNC1 | 152 | 74.3421 | no | 0.4863 | 0.5933 | yes | 0.0398 | 1.0000 | 9.7655 | between 0% and 25% | 1.0000 |
| ARID3A | 242 | 97.6974 | yes | 0.5327 | 0.0380 | no | 0.3114 | 1.0000 | 3.2380 | between 0% and 25% | 1.0000 |
| SOX8 | 114 | 50.0000 | no | 0.6535 | 2.1237 | yes | 0.0266 | 1.0000 | 17.0810 | between 0% and 25% | 0.9550 |

### primary_304, S1

| factor | cells | cell pctl in set | high_depth_candidate | delta | depth floor | below floor | p | q | % of gap | bracket | power at 25% at this depth |
|---|---|---|---|---|---|---|---|---|---|---|---|
| FOXD2 | 307 | 99.6711 | yes | 47.5927 | 43.2887 | no | 0.4662 | 0.9387 | NA | above the 0% point (44.912); no approach to express | 0.9350 |
| JDP2 | 331 | 100.0000 | yes | 48.2105 | 41.1527 | no | 1.0000 | 1.0000 | NA | above the 0% point (42.663); no approach to express | 0.9400 |
| ELK3 | 284 | 99.3421 | yes | 49.6092 | 45.6637 | no | 0.9992 | 1.0000 | NA | above the 0% point (47.581); no approach to express | 0.7350 |
| HMX1 | 280 | 99.0132 | yes | 53.7484 | 46.4151 | no | 0.9895 | 1.0000 | NA | above the 0% point (47.979); no approach to express | 0.8350 |
| NFIL3 | 255 | 98.0263 | yes | 55.4402 | 49.6076 | no | 0.7562 | 1.0000 | NA | above the 0% point (51.370); no approach to express | 0.8100 |

### all_1836, S2A_k20

| factor | cells | cell pctl in set | high_depth_candidate | delta | depth floor | below floor | p | q | % of gap | bracket | power at 25% at this depth |
|---|---|---|---|---|---|---|---|---|---|---|---|
| DMRTA1 | 183 | 94.0087 | yes | -2.5876 | -0.8752 | yes | 0.0000 | 0.0184 | 13.8422 | between 0% and 25% | 1.0000 |
| RELA | 179 | 92.9194 | yes | -2.3627 | -1.0005 | yes | 0.0000 | 0.0184 | 13.0127 | between 0% and 25% | 1.0000 |
| EGR3 | 71 | 13.0719 | no | -2.3610 | 0.9528 | yes | 0.0004 | 0.0964 | 27.8152 | between 25% and 50% | 0.8850 |
| ZNF473 | 111 | 56.4815 | no | -2.2928 | -0.6335 | yes | 0.0002 | 0.0734 | 16.6681 | between 0% and 25% | 0.9850 |
| GSC2 | 166 | 89.9237 | no | -2.1765 | -0.8813 | yes | 0.0001 | 0.0413 | 12.5008 | between 0% and 25% | 0.9950 |

### all_1836, S2A_k50

| factor | cells | cell pctl in set | high_depth_candidate | delta | depth floor | below floor | p | q | % of gap | bracket | power at 25% at this depth |
|---|---|---|---|---|---|---|---|---|---|---|---|
| GSC2 | 166 | 89.9237 | no | -1.5425 | -0.0419 | yes | 0.0003 | 0.3213 | 13.1833 | between 0% and 25% | 0.9950 |
| RELA | 179 | 92.9194 | yes | -1.4725 | -0.1392 | yes | 0.0003 | 0.3213 | 12.4452 | between 0% and 25% | 1.0000 |
| AFF1 | 177 | 92.3747 | yes | -1.2575 | -0.3630 | yes | 0.0010 | 0.6120 | 11.6522 | between 0% and 25% | 1.0000 |
| CRTC1 | 190 | 95.2070 | yes | -1.0672 | -0.0175 | yes | 0.0016 | 0.7436 | 10.6824 | between 0% and 25% | 1.0000 |
| DMRTA1 | 183 | 94.0087 | yes | -0.9135 | 0.1205 | yes | 0.0051 | 1.0000 | 10.4997 | between 0% and 25% | 1.0000 |

### all_1836, S2A_k100

| factor | cells | cell pctl in set | high_depth_candidate | delta | depth floor | below floor | p | q | % of gap | bracket | power at 25% at this depth |
|---|---|---|---|---|---|---|---|---|---|---|---|
| GSC2 | 166 | 89.9237 | no | -0.9874 | 0.4601 | yes | 0.0006 | 0.5692 | 13.3751 | between 0% and 25% | 0.9950 |
| RELA | 179 | 92.9194 | yes | -0.9867 | 0.3233 | yes | 0.0006 | 0.5692 | 12.8518 | between 0% and 25% | 1.0000 |
| AFF1 | 177 | 92.3747 | yes | -0.5266 | 0.1430 | yes | 0.0031 | 1.0000 | 11.1655 | between 0% and 25% | 1.0000 |
| CRTC1 | 190 | 95.2070 | yes | -0.4256 | 0.6096 | yes | 0.0039 | 1.0000 | 10.3638 | between 0% and 25% | 1.0000 |
| DMRTA1 | 183 | 94.0087 | yes | -0.2674 | 0.6715 | yes | 0.0103 | 1.0000 | 10.2484 | between 0% and 25% | 1.0000 |

### all_1836, S1

| factor | cells | cell pctl in set | high_depth_candidate | delta | depth floor | below floor | p | q | % of gap | bracket | power at 25% at this depth |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ZBTB46 | 341 | 100.0000 | yes | 43.0018 | 40.1237 | no | 0.9129 | 1.0000 | NA | above the 0% point (41.679); no approach to express | 0.9550 |
| REPIN1 | 317 | 99.8911 | yes | 45.6856 | 42.5139 | no | 0.9429 | 1.0000 | NA | above the 0% point (43.951); no approach to express | 0.9400 |
| FOXD2 | 307 | 99.8366 | yes | 47.5927 | 43.2887 | no | 0.4662 | 0.9416 | NA | above the 0% point (44.912); no approach to express | 0.9350 |
| ZNF789 | 287 | 99.7277 | yes | 47.6904 | 45.4511 | no | 0.6685 | 1.0000 | NA | above the 0% point (47.192); no approach to express | 0.8500 |
| JDP2 | 331 | 99.9455 | yes | 48.2105 | 41.1527 | no | 1.0000 | 1.0000 | NA | above the 0% point (42.663); no approach to express | 0.9400 |

**The sentence, filled in.** The largest apparent approach anywhere in the screen corresponds to at most 32.1% of the young-old gap, measured with a statistic that detects a 25% approach at this depth with power 0.980. It is held by RUNX1 (54 cells). The bound is computed for every one of the 1,836 factors, not only the five listed above, restricted to factors whose change in distance is actually negative.

Two qualifications on that number. First, the 0.980 power figure is at 106 cells; at RUNX1's own depth of 54 cells the power at a 25% plant is 0.830, so the statistic is weaker there than the sentence's headline figure suggests.
Second, 694 factors with a **positive** change in distance can also be placed on the planted curve by interpolation, because at low depth the 0%-dose median sits far above zero: the largest is NCL (28 cells, change in distance +9.94), which interpolates to 40.8% of the gap while moving away from the young target. Reading that as an approach is the same category error the criterion amendment removed, so such factors are excluded from the bound above. Had they been included the bound would read 40.8%.

This is descriptive and carries no decision rule. It is not an effect size in aged cells: the screen has no aged donor.

## 4. D2: is the permutation null calibrated? (gated)

SOUTH7's permutation test was run on draws that contain no effect by construction - the depth-matched null draws themselves - at 4 cell counts [50, 106, 180, 341], 200 draws each, in all four spaces. **n_perm = 10,000 here, not SOUTH7's 100,000**: the question is the rate at 0.05, not a tail at 1e-5, and the pre-registration permits it. Each draw's own displacement vector is permuted exactly as Stage A2 permutes a factor's.

| space | cells | draws | fraction p <= 0.05 | fraction p <= 0.001 | median p |
|---|---|---|---|---|---|
| S2A_k20 | 50 | 200 | 0.0000 | 0.0000 | 1.0000 |
| S2A_k20 | 106 | 200 | 0.0350 | 0.0000 | 0.9475 |
| S2A_k20 | 180 | 200 | 0.1100 | 0.0200 | 0.8044 |
| S2A_k20 | 341 | 200 | 0.0950 | 0.0150 | 0.7301 |
| S2A_k50 | 50 | 200 | 0.0000 | 0.0000 | 1.0000 |
| S2A_k50 | 106 | 200 | 0.0000 | 0.0000 | 0.9991 |
| S2A_k50 | 180 | 200 | 0.0400 | 0.0000 | 0.9296 |
| S2A_k50 | 341 | 200 | 0.0700 | 0.0100 | 0.7507 |
| S2A_k100 | 50 | 200 | 0.0000 | 0.0000 | 1.0000 |
| S2A_k100 | 106 | 200 | 0.0000 | 0.0000 | 0.9997 |
| S2A_k100 | 180 | 200 | 0.0300 | 0.0000 | 0.9525 |
| S2A_k100 | 341 | 200 | 0.0650 | 0.0100 | 0.7990 |
| S1 | 50 | 200 | 0.0900 | 0.0050 | 0.6783 |
| S1 | 106 | 200 | 0.1200 | 0.0150 | 0.6474 |
| S1 | 180 | 200 | 0.1400 | 0.0200 | 0.6049 |
| S1 | 341 | 200 | 0.0950 | 0.0150 | 0.5566 |

### Gate reading that fired, per space (pre-registered wording)

- **S2A_k20: NOT CALIBRATED, in both directions.** The permutation null flags between 0.000 and 0.110 of true-null draws at the 0.05 level in S2A_k20 depending on cell count. Where it exceeds 0.08 it cannot be read as evidence that a factor's direction is special; where it falls below 0.03 it cannot be used to withhold a call.
- **S2A_k50: CONSERVATIVE.** The permutation null flags as little as 0.000 of true-null draws at the 0.05 level in S2A_k50, so it cannot be used to withhold a call.
- **S2A_k100: CONSERVATIVE.** The permutation null flags as little as 0.000 of true-null draws at the 0.05 level in S2A_k100, so it cannot be used to withhold a call.
- **S1: ANTI-CONSERVATIVE.** The permutation null flags up to 0.140 of true-null draws at the 0.05 level in S1, so it cannot be read as evidence that a factor's direction is special.

### Consequence, as fixed in the pre-registration

The null is **not** calibrated in S2A_k20. The pre-registration fixed the consequence in advance: the SOUTH7 primary call is re-read using the depth-matched floor and the amended criterion alone, without the permutation condition. **That re-reading is the primary reading of SOUTH8.** SOUTH7's reading is marked superseded, not edited. Both side by side:

| factor set | m | SOUTH7 reading (superseded) | SOUTH8 primary reading (floor + delta<0) | factors called |
|---|---|---|---|---|
| primary_304 | 304 | 0 | 24 | AFF1;AFF4;BNC1;BNC2;CDX2;CEBPA;EGR2;EGR3;FOXC1;FOXN3;HNF1B;HOXA13;KAT2A;KLF12;KLF15;MYRF;NFATC1;NKX3-2;NR2F1;PLAGL1;PTF1A;SOX8;TSHZ1;ZXDA |
| all_1836 | 1836 | 4 | 125 | AFF1;AFF4;ASCL4;ATMIN;BAZ2B;BCLAF1;BLZF1;BNC1;BNC2;BTAF1;CDX2;CEBPA;CRTC1;CRTC3;CTBP1;DAXX;DBX2;DDB1;DMRTA1;DPRX;EGR2;EGR3;ETV2;FOXC1;FOXJ3;FOXN3;FUBP1;GABPB1;GLMP;GRHL2;GSC2;GSX1;GTF3C2;HMGXB3;HNF1B;HOXA10;HOXA13;HOXD1;IFI16;KAT2A;KLF12;KLF15;L3MBTL4;LDB1;LIN28B;MAEL;MAFB;MEF2B;MEF2D;MYRF;NFATC1;NFATC4;NFKB2;NFX1;NHLH1;NHLH2;NKRF;NKX2-5;NKX3-2;NR2F1;OLIG1;OVOL3;PAX4;PBX3;PLAGL1;POGZ;POU5F1;POU6F1;PTF1A;PUF60;RAD21;RCOR2;RELA;RUNX1;SCRT1;SKP2;SMAD3;SMAD9;SNAPC2;SOX8;SP9;SPIB;TFDP1;THAP1;TLE1;TRIM33;TSHZ1;VPS72;WDHD1;WNK1;WT1;XRCC4;YY2;ZBTB24;ZFAT;ZFP64;ZFP90;ZNF184;ZNF189;ZNF200;ZNF222;ZNF224;ZNF263;ZNF267;ZNF268;ZNF275;ZNF286B;ZNF418;ZNF430;ZNF473;ZNF492;ZNF502;ZNF561;ZNF615;ZNF652;ZNF683;ZNF70;ZNF780B;ZNF813;ZNF84;ZNF865;ZSCAN1;ZSCAN23;ZSCAN31;ZXDA |

The permutation condition is what produced SOUTH7's zero. Without it, the count in the primary space is 24 of 304 in primary_304, 125 of 1836 in all_1836.

**Read that count with its chance baseline.** The depth-matched floor is a 5th percentile, so about 5% of factors fall below it with no effect at all: 15.2 of 304, and 91.8 of 1,836. The counts above are therefore not that many discoveries. In the primary set, 24 factors are called against 15.2 expected by chance, an excess of roughly nine to twelve factors spread across the named list; D3 below shows the excess is larger than chance (27 below the floor against a 95% interval of 8 to 23) and is not explained by read depth. No individual factor on the list is established by this, and the amended criterion carries no multiple-testing control of its own now that the permutation condition has been removed.

Why this null can fail: a uniform shuffle of gene coordinates destroys the correlation structure a real expression vector carries, so the shuffled vectors are not exchangeable with the observed one. The depth-matched null does not have this problem - it is built from real non-targeting cells - which is why the two nulls disagree.

## 5. D3: is the null matched on the right quantity? (gated)

### (a) Total UMI against the statistic

| factor set | space | factors | Spearman total UMI vs delta | Spearman cells vs delta | Spearman total UMI vs cells |
|---|---|---|---|---|---|
| primary_304 | S2A_k20 | 304 | -0.1766 | -0.1966 | 0.9889 |
| primary_304 | S2A_k50 | 304 | -0.5629 | -0.5797 | 0.9889 |
| primary_304 | S2A_k100 | 304 | -0.6801 | -0.6936 | 0.9889 |
| primary_304 | S1 | 304 | -0.9929 | -0.9855 | 0.9889 |
| all_1836 | S2A_k20 | 1836 | -0.3435 | -0.3532 | 0.9906 |
| all_1836 | S2A_k50 | 1836 | -0.7138 | -0.7205 | 0.9906 |
| all_1836 | S2A_k100 | 1836 | -0.7978 | -0.8014 | 0.9906 |
| all_1836 | S1 | 1836 | -0.9955 | -0.9869 | 0.9906 |

### (b) A null matched on total reads instead of cell count

Draws are non-targeting pseudobulks whose total UMI is matched to the factor's, by the SOUTH6REC scheme (disjoint blocks within a permutation, reference = complement of the draw). The stopping point is whichever of the two candidate cells lands closer to the target, so draws are centred on it rather than all falling to one side - stopping on first crossing would make every draw systematically thinner than its factor, which would lower the floor and bias this test toward 'explained'. One null is computed per 1% geometric bin of total UMI and shared by the factors in that bin; 200 bins cover all 1,836 factors.

Matching achieved: the realised offset between a factor's total UMI and the total of each of its draws is at most 0.1134 in the worst case and 0.0129 at the median, against the 2% tolerance; 90.5% of factors have every draw inside tolerance.

Cells per draw, against the cell-count-matched null, at a few representative bins:

| target total UMI | draws | median cells/draw | min | max | 5th-pct floor | realised tol (median) | realised tol (max) | mutually disjoint |
|---|---|---|---|---|---|---|---|---|
| 76023.00000 | 200 | 12.00000 | 9 | 15 | 39.70020 | 0.02207 | 0.11344 | 200 |
| 298387.68757 | 200 | 49.00000 | 42 | 56 | 3.35151 | 0.00472 | 0.02461 | 200 |
| 398946.53352 | 200 | 65.00000 | 57 | 74 | 1.30788 | 0.00354 | 0.01427 | 200 |
| 531998.03008 | 200 | 87.50000 | 77 | 97 | -0.16969 | 0.00303 | 0.01369 | 160 |
| 703090.58489 | 200 | 115.00000 | 103 | 124 | -0.75284 | 0.00233 | 0.01074 | 121 |
| 936488.33011 | 200 | 154.00000 | 141 | 164 | -0.95373 | 0.00177 | 0.00783 | 91 |
| 1238451.78164 | 200 | 203.00000 | 191 | 221 | -0.92433 | 0.00114 | 0.00556 | 68 |
| 2052725.00000 | 200 | 338.00000 | 319 | 355 | -1.04678 | 0.00070 | 0.00331 | 41 |

### (c) Recount against the read-matched floor

| factor set | space | m | below cell-count floor | below total-UMI floor | amended call, cell-count floor | amended call, UMI floor |
|---|---|---|---|---|---|---|
| primary_304 | S2A_k20 | 304 | 27 | 27 | 24 | 22 |
| primary_304 | S2A_k50 | 304 | 23 | 24 | 4 | 4 |
| primary_304 | S2A_k100 | 304 | 20 | 22 | 1 | 1 |
| primary_304 | S1 | 304 | 3 | 3 | 0 | 0 |
| all_1836 | S2A_k20 | 1836 | 153 | 150 | 125 | 121 |
| all_1836 | S2A_k50 | 1836 | 130 | 143 | 19 | 20 |
| all_1836 | S2A_k100 | 1836 | 118 | 137 | 6 | 7 |
| all_1836 | S1 | 1836 | 48 | 50 | 0 | 0 |

### Gate reading that fired (pre-registered wording)

- **SURVIVES.** The excess survives read-depth matching: 27 of 304 fall below a floor matched on total reads, against 15.2 expected by chance.

(5% of 304 is 15.2; the 95% binomial interval is 8 to 23 inclusive, as the pre-registration fixed.)

## 6. D4: the all-1,836 reading (report-only)

SOUTH7's permutation run completed for all 1,836 factors (100,000 draws each; the 1,532 outside the published set finished in 617 minutes). q is computed over m = 1,836 here.

| space | m | amended (floor+delta<0) | amended+q | SOUTH7 criterion | below depth floor | q <= 0.05 | delta < 0 | tied at min p | smallest q |
|---|---|---|---|---|---|---|---|---|---|
| S2A_k20 | 1836 | 125 | 4 | 4 | 153 | 4 | 214 | 1 | 0.0184 |
| S2A_k50 | 1836 | 19 | 0 | 0 | 130 | 0 | 21 | 0 | 0.3213 |
| S2A_k100 | 1836 | 6 | 0 | 0 | 118 | 0 | 7 | 0 | 0.5692 |
| S1 | 1836 | 0 | 0 | 6 | 48 | 210 | 0 | 28 | 0.0007 |

SOUTH6's recurring set, now complete - CERS3 and REST included, the two that fell outside the published primary set:

| factor | in published set | cells | cell pctl in screen | high_depth_candidate | delta | depth floor | q | amended call | ruler scores younger |
|---|---|---|---|---|---|---|---|---|---|
| CERS3 | no | 110 | 55.5011 | no | 2.0331 | -0.9122 | 1.0000 | no | no |
| ETV3L | yes | 168 | 90.5773 | yes | 9.3303 | -1.0206 | 1.0000 | no | yes |
| FOSB | yes | 228 | 98.5294 | yes | 9.1936 | -0.9982 | 1.0000 | no | yes |
| FOXF1 | yes | 137 | 77.9956 | no | 3.8218 | -0.7802 | 1.0000 | no | no |
| FOXI3 | yes | 171 | 91.1765 | yes | 4.2143 | -0.7284 | 1.0000 | no | yes |
| GATA1 | yes | 239 | 98.9651 | yes | 7.6969 | -0.9565 | 1.0000 | no | yes |
| HELZ2 | yes | 140 | 79.5752 | no | -0.3020 | -0.6944 | 1.0000 | no | no |
| JDP2 | yes | 331 | 99.9455 | yes | 8.2053 | -1.0294 | 1.0000 | no | yes |
| NEUROD1 | yes | 185 | 94.2266 | yes | 7.5108 | -0.9873 | 1.0000 | no | yes |
| OSR2 | yes | 183 | 94.0087 | yes | 1.1661 | -0.8752 | 1.0000 | no | no |
| REST | no | 227 | 98.3115 | yes | 2.0698 | -1.0508 | 1.0000 | no | yes |

**Factors called under the amended criterion that carry the high_depth_candidate flag** (cell count above the 90th percentile of the set, so reported as such and not as a result): S2A_k20: AFF1, AFF4, CRTC1, DMRTA1, HOXD1, RELA, SMAD3, ZNF652; S2A_k50: AFF1, AFF4, CRTC1, DMRTA1, HOXD1, RELA, SMAD3, TCF25, ZHX3, ZNF652; S2A_k100: AFF1, CRTC1, DMRTA1, HOXD1, RELA.

#### Called under the amended criterion in S2A_k20, all 1,836

| factor | in published set | cells | cell pctl | high_depth_candidate | delta | depth floor | p | q | guide agreement | identity drop | ruler score change |
|---|---|---|---|---|---|---|---|---|---|---|---|
| DMRTA1 | no | 183 | 94.0087 | yes | -2.5876 | -0.8752 | 0.0000 | 0.0184 | 0.1425 | -0.0142 | -0.0561 |
| RELA | no | 179 | 92.9194 | yes | -2.3627 | -1.0005 | 0.0000 | 0.0184 | 0.1444 | -0.2269 | 0.5372 |
| EGR3 | yes | 71 | 13.0719 | no | -2.3610 | 0.9528 | 0.0004 | 0.0964 | 0.2597 | -0.0765 | -0.1754 |
| ZNF473 | no | 111 | 56.4815 | no | -2.2928 | -0.6335 | 0.0002 | 0.0734 | 0.1938 | -0.0871 | 0.6889 |
| GSC2 | no | 166 | 89.9237 | no | -2.1765 | -0.8813 | 0.0001 | 0.0413 | 0.1500 | 0.0028 | -0.5831 |
| CRTC1 | no | 190 | 95.2070 | yes | -2.1713 | -0.8638 | 0.0000 | 0.0245 | 0.1383 | 0.0185 | 0.2843 |
| FOXN3 | yes | 128 | 71.7320 | no | -2.1206 | -0.8317 | 0.0002 | 0.0734 | 0.1635 | -0.0151 | -0.4337 |
| ZNF780B | no | 128 | 71.7320 | no | -2.0914 | -0.8317 | 0.0003 | 0.0734 | 0.1664 | 0.0015 | -0.0600 |
| GTF3C2 | no | 113 | 58.8780 | no | -2.0783 | -0.5914 | 0.0007 | 0.1025 | 0.1671 | 0.0466 | 0.4816 |
| FUBP1 | no | 94 | 36.8736 | no | -1.9470 | 0.0053 | 0.0013 | 0.1426 | 0.2331 | 0.0885 | 1.0310 |
| MEF2B | no | 73 | 15.3595 | no | -1.9137 | 0.2844 | 0.0018 | 0.1720 | 0.2277 | 0.0088 | 0.9094 |
| SOX8 | yes | 114 | 60.0763 | no | -1.8734 | -0.7587 | 0.0022 | 0.1724 | 0.1999 | 0.0244 | -1.8477 |
| BTAF1 | no | 92 | 34.0959 | no | -1.8355 | 0.1596 | 0.0013 | 0.1426 | 0.2182 | -0.0443 | 1.0311 |
| XRCC4 | no | 75 | 17.2113 | no | -1.8319 | 0.7067 | 0.0032 | 0.2131 | 0.2296 | -0.0463 | 0.4060 |
| AFF1 | yes | 177 | 92.3747 | yes | -1.8316 | -1.2665 | 0.0006 | 0.1001 | 0.1441 | -0.0457 | -0.1530 |
| GLMP | no | 97 | 40.3050 | no | -1.8244 | -0.2455 | 0.0010 | 0.1412 | 0.2209 | 0.0974 | 0.6415 |
| ZBTB24 | no | 134 | 75.4902 | no | -1.7922 | -1.1467 | 0.0011 | 0.1426 | 0.1935 | -0.0204 | 1.3231 |
| SPIB | no | 104 | 48.8017 | no | -1.7825 | -0.2712 | 0.0021 | 0.1724 | 0.1818 | 0.0196 | -0.2703 |
| CTBP1 | no | 150 | 84.6405 | no | -1.7689 | -0.8147 | 0.0005 | 0.1001 | 0.1549 | -0.0239 | 0.7878 |
| SP9 | no | 112 | 57.6797 | no | -1.7198 | -0.6785 | 0.0026 | 0.1997 | 0.2239 | 0.0840 | 0.2011 |
| BNC1 | yes | 152 | 85.7843 | no | -1.7147 | -1.1506 | 0.0013 | 0.1426 | 0.1599 | -0.0320 | -0.1231 |
| SMAD9 | no | 84 | 25.4357 | no | -1.7145 | 0.5761 | 0.0031 | 0.2128 | 0.2010 | 0.0352 | -0.2340 |
| RAD21 | no | 140 | 79.5752 | no | -1.6879 | -0.6944 | 0.0016 | 0.1591 | 0.1633 | -0.0730 | 0.8908 |
| ATMIN | no | 103 | 47.3856 | no | -1.6748 | -0.4579 | 0.0031 | 0.2128 | 0.1891 | -0.1060 | 0.3996 |
| PLAGL1 | yes | 118 | 63.5621 | no | -1.6283 | -0.7335 | 0.0045 | 0.2539 | 0.2029 | -0.1005 | -1.2399 |
| TLE1 | no | 149 | 84.2593 | no | -1.5984 | -0.9321 | 0.0021 | 0.1724 | 0.1682 | 0.1137 | 0.5085 |
| SMAD3 | no | 222 | 98.0392 | yes | -1.5833 | -0.9325 | 0.0006 | 0.1001 | 0.1143 | 0.0288 | -0.2127 |
| KLF12 | yes | 117 | 62.7451 | no | -1.5769 | -0.5573 | 0.0044 | 0.2539 | 0.1884 | -0.0353 | 0.3469 |
| TFDP1 | no | 140 | 79.5752 | no | -1.5715 | -0.6944 | 0.0020 | 0.1724 | 0.1489 | -0.0062 | 0.3758 |
| MAFB | no | 117 | 62.7451 | no | -1.5587 | -0.5573 | 0.0047 | 0.2539 | 0.1908 | -0.0630 | -0.4586 |
| OVOL3 | no | 89 | 30.8824 | no | -1.5305 | -0.1242 | 0.0048 | 0.2539 | 0.2069 | 0.0727 | -0.5509 |
| GABPB1 | no | 117 | 62.7451 | no | -1.5207 | -0.5573 | 0.0040 | 0.2423 | 0.2011 | -0.0192 | 0.6784 |
| ZNF70 | no | 127 | 71.0240 | no | -1.4616 | -0.8346 | 0.0058 | 0.2873 | 0.1650 | 0.0002 | -1.0667 |
| ZFP64 | no | 161 | 88.1264 | no | -1.3967 | -1.1324 | 0.0046 | 0.2539 | 0.1746 | -0.0729 | -0.0517 |
| ZNF189 | no | 103 | 47.3856 | no | -1.3804 | -0.4579 | 0.0064 | 0.3087 | 0.2031 | -0.0486 | -0.4432 |
| POU6F1 | no | 135 | 76.4161 | no | -1.3739 | -0.5556 | 0.0057 | 0.2873 | 0.1564 | 0.0076 | 0.2386 |
| HOXD1 | no | 181 | 93.5185 | yes | -1.3622 | -1.1630 | 0.0040 | 0.2423 | 0.1684 | 0.0437 | -0.0917 |
| ZNF502 | no | 107 | 52.1786 | no | -1.3472 | -0.1591 | 0.0087 | 0.3984 | 0.1829 | -0.0760 | 0.7985 |
| AFF4 | yes | 199 | 96.4052 | yes | -1.3419 | -1.0042 | 0.0030 | 0.2128 | 0.1466 | -0.0521 | 0.1591 |
| TRIM33 | no | 119 | 64.3791 | no | -1.3357 | -0.8300 | 0.0092 | 0.4129 | 0.1677 | -0.0257 | 1.2983 |
| FOXJ3 | no | 84 | 25.4357 | no | -1.3227 | 0.5761 | 0.0112 | 0.4590 | 0.2791 | -0.0448 | 0.9728 |
| DPRX | no | 112 | 57.6797 | no | -1.3118 | -0.6785 | 0.0119 | 0.4656 | 0.1666 | -0.0602 | -0.4492 |
| NFKB2 | no | 133 | 74.7277 | no | -1.2842 | -0.6818 | 0.0079 | 0.3719 | 0.1608 | -0.0383 | -0.2844 |
| ZFP90 | no | 141 | 80.2288 | no | -1.2451 | -0.7600 | 0.0095 | 0.4166 | 0.1371 | 0.0510 | 0.3042 |
| ZNF268 | no | 117 | 62.7451 | no | -1.2225 | -0.5573 | 0.0171 | 0.5920 | 0.2055 | -0.1611 | 1.3763 |
| KAT2A | yes | 104 | 48.8017 | no | -1.2149 | -0.2712 | 0.0135 | 0.5179 | 0.1869 | 0.0223 | 0.6764 |
| ZNF275 | no | 73 | 15.3595 | no | -1.2074 | 0.2844 | 0.0154 | 0.5544 | 0.2233 | 0.1387 | 1.1698 |
| BAZ2B | no | 119 | 64.3791 | no | -1.2036 | -0.8300 | 0.0192 | 0.6242 | 0.2070 | -0.0289 | 0.0849 |
| OLIG1 | no | 118 | 63.5621 | no | -1.1742 | -0.7335 | 0.0193 | 0.6242 | 0.1795 | -0.0006 | 0.7658 |
| GRHL2 | no | 135 | 76.4161 | no | -1.1643 | -0.5556 | 0.0116 | 0.4650 | 0.1559 | 0.0289 | 1.3330 |
| HMGXB3 | no | 122 | 66.9390 | no | -1.1604 | -0.6602 | 0.0152 | 0.5544 | 0.1826 | 0.0885 | -0.6063 |
| ZNF652 | no | 178 | 92.5926 | yes | -1.1373 | -1.0544 | 0.0112 | 0.4590 | 0.1424 | -0.0763 | 0.6646 |
| YY2 | no | 102 | 46.1874 | no | -1.1361 | -0.2384 | 0.0166 | 0.5879 | 0.1759 | -0.0243 | -0.3043 |
| CDX2 | yes | 80 | 21.7865 | no | -1.1307 | -0.0351 | 0.0207 | 0.6242 | 0.2287 | -0.2897 | 0.5348 |
| BNC2 | yes | 120 | 65.0871 | no | -1.1147 | -0.5986 | 0.0253 | 0.7135 | 0.1798 | 0.1896 | -0.0548 |
| THAP1 | no | 102 | 46.1874 | no | -1.0966 | -0.2384 | 0.0206 | 0.6242 | 0.1829 | 0.0198 | -0.2689 |
| ZNF184 | no | 99 | 43.1917 | no | -1.0855 | -0.1996 | 0.0145 | 0.5437 | 0.1985 | 0.0222 | 0.5229 |
| DDB1 | no | 122 | 66.9390 | no | -1.0564 | -0.6602 | 0.0182 | 0.6188 | 0.2002 | -0.0059 | 0.1780 |
| RUNX1 | no | 54 | 3.4314 | no | -1.0474 | 2.0352 | 0.0206 | 0.6242 | 0.2663 | -0.0857 | 1.7521 |
| ZNF84 | no | 104 | 48.8017 | no | -1.0053 | -0.2712 | 0.0251 | 0.7135 | 0.1976 | 0.0330 | -0.5977 |
| ZNF200 | no | 125 | 69.7712 | no | -0.9998 | -0.5991 | 0.0207 | 0.6242 | 0.1589 | -0.0626 | 0.1485 |
| KLF15 | yes | 143 | 81.3181 | no | -0.9909 | -0.4996 | 0.0234 | 0.6926 | 0.1693 | -0.0325 | -0.0941 |
| ZSCAN1 | no | 123 | 67.8105 | no | -0.9851 | -0.6664 | 0.0249 | 0.7135 | 0.1736 | 0.0643 | 1.0623 |
| NKRF | no | 127 | 71.0240 | no | -0.9829 | -0.8346 | 0.0261 | 0.7274 | 0.1860 | -0.0323 | 0.5849 |
| ZXDA | yes | 84 | 25.4357 | no | -0.9624 | 0.5761 | 0.0319 | 0.8262 | 0.2671 | -0.0536 | 0.1704 |
| ASCL4 | no | 81 | 22.7669 | no | -0.9423 | 0.2220 | 0.0330 | 0.8262 | 0.2377 | 0.0910 | 0.9263 |
| NFATC1 | yes | 142 | 80.7734 | no | -0.9403 | -0.7434 | 0.0203 | 0.6242 | 0.1807 | 0.0181 | 0.5786 |
| SNAPC2 | no | 115 | 60.8932 | no | -0.9223 | -0.6671 | 0.0367 | 0.8525 | 0.1817 | -0.0215 | -0.2478 |
| MYRF | yes | 130 | 72.9303 | no | -0.9079 | -0.8448 | 0.0333 | 0.8262 | 0.1881 | 0.2207 | -0.2304 |
| ZNF813 | no | 124 | 68.6819 | no | -0.9035 | -0.5545 | 0.0326 | 0.8262 | 0.1769 | -0.0472 | 0.5591 |
| NKX3-2 | yes | 94 | 36.8736 | no | -0.8888 | 0.0053 | 0.0275 | 0.7528 | 0.1968 | 0.1360 | 0.5134 |
| SCRT1 | no | 155 | 86.6558 | no | -0.8792 | -0.5114 | 0.0388 | 0.8690 | 0.1581 | -0.0587 | 0.0026 |
| EGR2 | yes | 133 | 74.7277 | no | -0.8786 | -0.6818 | 0.0343 | 0.8373 | 0.1784 | -0.0865 | -0.2135 |
| ZNF615 | no | 107 | 52.1786 | no | -0.8771 | -0.1591 | 0.0325 | 0.8262 | 0.2586 | 0.0585 | -1.1983 |
| ZNF418 | no | 78 | 20.2070 | no | -0.8717 | -0.1571 | 0.0347 | 0.8373 | 0.2365 | 0.0035 | 1.0399 |
| LDB1 | no | 80 | 21.7865 | no | -0.8621 | -0.0351 | 0.0357 | 0.8512 | 0.2394 | -0.0307 | 0.3776 |
| PTF1A | yes | 125 | 69.7712 | no | -0.8604 | -0.5991 | 0.0378 | 0.8666 | 0.1588 | -0.0140 | 1.3133 |
| NHLH2 | no | 108 | 52.9412 | no | -0.8309 | -0.6407 | 0.0429 | 0.9388 | 0.1848 | -0.0450 | 1.3952 |
| NFATC4 | no | 107 | 52.1786 | no | -0.7909 | -0.1591 | 0.0501 | 1.0000 | 0.2180 | 0.1347 | -0.5743 |
| ZNF224 | no | 91 | 32.7342 | no | -0.7866 | 0.0911 | 0.0386 | 0.8690 | 0.2095 | 0.0494 | -0.4876 |
| HOXA13 | yes | 96 | 39.3246 | no | -0.7382 | -0.1069 | 0.0421 | 0.9310 | 0.2263 | 0.0454 | -0.1781 |
| BLZF1 | no | 104 | 48.8017 | no | -0.7314 | -0.2712 | 0.0506 | 1.0000 | 0.1881 | -0.0478 | 1.5923 |
| IFI16 | no | 77 | 19.3900 | no | -0.7192 | 0.3335 | 0.0462 | 0.9745 | 0.2208 | 0.0085 | -0.1499 |
| POGZ | no | 125 | 69.7712 | no | -0.7102 | -0.5991 | 0.0580 | 1.0000 | 0.1678 | 0.0492 | 0.5706 |
| MAEL | no | 135 | 76.4161 | no | -0.6924 | -0.5556 | 0.0558 | 1.0000 | 0.2079 | -0.0354 | 0.1509 |
| HOXA10 | no | 116 | 61.6013 | no | -0.6388 | -0.5192 | 0.0862 | 1.0000 | 0.1736 | -0.0008 | 0.6127 |
| NR2F1 | yes | 76 | 18.2462 | no | -0.6318 | 0.1834 | 0.0580 | 1.0000 | 0.2276 | -0.0105 | 0.5711 |
| GSX1 | no | 131 | 73.4749 | no | -0.6308 | -0.6278 | 0.0724 | 1.0000 | 0.1559 | 0.0496 | -0.7118 |
| FOXC1 | yes | 120 | 65.0871 | no | -0.6239 | -0.5986 | 0.0713 | 1.0000 | 0.2561 | 0.0112 | 0.6902 |
| MEF2D | no | 100 | 44.0632 | no | -0.6171 | -0.2085 | 0.0549 | 1.0000 | 0.2092 | -0.0236 | 1.1184 |
| DBX2 | no | 77 | 19.3900 | no | -0.6062 | 0.3335 | 0.0726 | 1.0000 | 0.2410 | -0.0426 | -0.0552 |
| PBX3 | no | 76 | 18.2462 | no | -0.6061 | 0.1834 | 0.0713 | 1.0000 | 0.2204 | 0.0261 | -0.7034 |
| ZNF865 | no | 143 | 81.3181 | no | -0.5873 | -0.4996 | 0.0687 | 1.0000 | 0.1907 | 0.0168 | 0.2363 |
| PUF60 | no | 81 | 22.7669 | no | -0.5526 | 0.2220 | 0.0723 | 1.0000 | 0.2099 | 0.0876 | -0.9535 |
| VPS72 | no | 143 | 81.3181 | no | -0.5433 | -0.4996 | 0.0927 | 1.0000 | 0.1412 | 0.0806 | 0.0667 |
| POU5F1 | no | 98 | 41.4488 | no | -0.5406 | -0.2147 | 0.0672 | 1.0000 | 0.1995 | 0.0658 | -0.0964 |
| TSHZ1 | yes | 77 | 19.3900 | no | -0.5359 | 0.3335 | 0.0702 | 1.0000 | 0.2479 | 0.0479 | -0.2600 |
| ZNF267 | no | 106 | 50.8170 | no | -0.4922 | -0.4660 | 0.0791 | 1.0000 | 0.1968 | 0.0136 | -0.2441 |
| ZNF263 | no | 89 | 30.8824 | no | -0.4810 | -0.1242 | 0.0826 | 1.0000 | 0.2121 | -0.0704 | 0.5917 |
| ZNF430 | no | 107 | 52.1786 | no | -0.4666 | -0.1591 | 0.0926 | 1.0000 | 0.2114 | -0.1486 | -0.0605 |
| ZNF561 | no | 90 | 32.0261 | no | -0.4610 | -0.0458 | 0.0790 | 1.0000 | 0.2607 | 0.0153 | -1.0203 |
| NKX2-5 | no | 81 | 22.7669 | no | -0.4539 | 0.2220 | 0.0875 | 1.0000 | 0.2702 | 0.0384 | -0.3161 |
| ZNF222 | no | 102 | 46.1874 | no | -0.4397 | -0.2384 | 0.1231 | 1.0000 | 0.1887 | 0.0159 | -0.1151 |
| NHLH1 | no | 91 | 32.7342 | no | -0.4366 | 0.0911 | 0.0705 | 1.0000 | 0.2191 | 0.0749 | 0.6783 |
| ZNF492 | no | 102 | 46.1874 | no | -0.4275 | -0.2384 | 0.1256 | 1.0000 | 0.1755 | 0.0963 | 0.7418 |
| CEBPA | yes | 74 | 16.2854 | no | -0.3894 | 0.7661 | 0.0765 | 1.0000 | 0.2682 | 0.0056 | 0.6765 |
| L3MBTL4 | no | 78 | 20.2070 | no | -0.3863 | -0.1571 | 0.0837 | 1.0000 | 0.2271 | 0.0616 | -0.2749 |
| WT1 | no | 102 | 46.1874 | no | -0.3759 | -0.2384 | 0.1285 | 1.0000 | 0.2120 | 0.0540 | 1.8201 |
| DAXX | no | 102 | 46.1874 | no | -0.3731 | -0.2384 | 0.0936 | 1.0000 | 0.2079 | 0.0116 | 0.3883 |
| ZNF286B | no | 88 | 30.1198 | no | -0.3727 | -0.1197 | 0.1104 | 1.0000 | 0.2350 | -0.0335 | -0.0191 |
| RCOR2 | no | 78 | 20.2070 | no | -0.3526 | -0.1571 | 0.0997 | 1.0000 | 0.2323 | -0.0689 | -0.1138 |
| LIN28B | no | 82 | 23.7473 | no | -0.3460 | -0.2571 | 0.1057 | 1.0000 | 0.2410 | 0.0713 | -0.4229 |
| NFX1 | no | 74 | 16.2854 | no | -0.3230 | 0.7661 | 0.1032 | 1.0000 | 0.2129 | 0.0682 | -0.4436 |
| ETV2 | no | 71 | 13.0719 | no | -0.2900 | 0.9528 | 0.0898 | 1.0000 | 0.2388 | 0.0587 | 0.4659 |
| SKP2 | no | 86 | 27.7778 | no | -0.2891 | -0.0586 | 0.1236 | 1.0000 | 0.2307 | 0.0390 | 0.5767 |
| PAX4 | no | 107 | 52.1786 | no | -0.2864 | -0.1591 | 0.1457 | 1.0000 | 0.1943 | 0.0071 | -0.1814 |
| CRTC3 | no | 96 | 39.3246 | no | -0.2571 | -0.1069 | 0.1155 | 1.0000 | 0.2507 | 0.0916 | -0.1116 |
| WNK1 | no | 78 | 20.2070 | no | -0.2392 | -0.1571 | 0.1319 | 1.0000 | 0.2357 | -0.0814 | 0.7753 |
| ZSCAN31 | no | 76 | 18.2462 | no | -0.2322 | 0.1834 | 0.1348 | 1.0000 | 0.2180 | 0.0633 | -0.5550 |
| WDHD1 | no | 71 | 13.0719 | no | -0.2082 | 0.9528 | 0.1301 | 1.0000 | 0.2306 | 0.0508 | 0.7333 |
| ZNF683 | no | 96 | 39.3246 | no | -0.1902 | -0.1069 | 0.1161 | 1.0000 | 0.2054 | -0.0113 | 0.4657 |
| ZSCAN23 | no | 94 | 36.8736 | no | -0.1561 | 0.0053 | 0.1517 | 1.0000 | 0.1834 | 0.0138 | -1.3138 |
| ZFAT | no | 94 | 36.8736 | no | -0.1463 | 0.0053 | 0.1680 | 1.0000 | 0.1879 | 0.0353 | 0.1408 |
| HNF1B | yes | 84 | 25.4357 | no | -0.0315 | 0.5761 | 0.1866 | 1.0000 | 0.2537 | 0.0668 | -0.3233 |
| BCLAF1 | no | 85 | 26.6340 | no | -0.0236 | 0.1084 | 0.2066 | 1.0000 | 0.2136 | 0.0692 | 1.6283 |

#### Called under the amended criterion in S2A_k50, all 1,836

| factor | in published set | cells | cell pctl | high_depth_candidate | delta | depth floor | p | q | guide agreement | identity drop | ruler score change |
|---|---|---|---|---|---|---|---|---|---|---|---|
| GSC2 | no | 166 | 89.9237 | no | -1.5425 | -0.0419 | 0.0003 | 0.3213 | 0.1500 | 0.0028 | -0.5831 |
| RELA | no | 179 | 92.9194 | yes | -1.4725 | -0.1392 | 0.0003 | 0.3213 | 0.1444 | -0.2269 | 0.5372 |
| AFF1 | yes | 177 | 92.3747 | yes | -1.2575 | -0.3630 | 0.0010 | 0.6120 | 0.1441 | -0.0457 | -0.1530 |
| CRTC1 | no | 190 | 95.2070 | yes | -1.0672 | -0.0175 | 0.0016 | 0.7436 | 0.1383 | 0.0185 | 0.2843 |
| DMRTA1 | no | 183 | 94.0087 | yes | -0.9135 | 0.1205 | 0.0051 | 1.0000 | 0.1425 | -0.0142 | -0.0561 |
| HOXD1 | no | 181 | 93.5185 | yes | -0.7439 | -0.1925 | 0.0087 | 1.0000 | 0.1684 | 0.0437 | -0.0917 |
| GTF3C2 | no | 113 | 58.8780 | no | -0.5749 | 1.3584 | 0.0086 | 1.0000 | 0.1671 | 0.0466 | 0.4816 |
| SOX8 | yes | 114 | 60.0763 | no | -0.5731 | 1.1536 | 0.0093 | 1.0000 | 0.1999 | 0.0244 | -1.8477 |
| SMAD3 | no | 222 | 98.0392 | yes | -0.5527 | -0.4631 | 0.0179 | 1.0000 | 0.1143 | 0.0288 | -0.2127 |
| BNC1 | yes | 152 | 85.7843 | no | -0.4569 | 0.0059 | 0.0151 | 1.0000 | 0.1599 | -0.0320 | -0.1231 |
| ZMYND11 | no | 164 | 89.2157 | no | -0.4284 | 0.1444 | 0.0544 | 1.0000 | 0.1838 | 0.0833 | -0.2213 |
| AFF4 | yes | 199 | 96.4052 | yes | -0.3399 | -0.0374 | 0.0338 | 1.0000 | 0.1466 | -0.0521 | 0.1591 |
| ZFP90 | no | 141 | 80.2288 | no | -0.2448 | 0.8773 | 0.0480 | 1.0000 | 0.1371 | 0.0510 | 0.3042 |
| MAFA | no | 164 | 89.2157 | no | -0.1931 | 0.1444 | 0.0764 | 1.0000 | 0.1354 | 0.0914 | 0.1834 |
| TCF25 | no | 201 | 96.6776 | yes | -0.1884 | -0.1733 | 0.0574 | 1.0000 | 0.1191 | 0.0649 | 0.3671 |
| ZNF780B | no | 128 | 71.7320 | no | -0.1412 | 1.1862 | 0.0291 | 1.0000 | 0.1664 | 0.0015 | -0.0600 |
| ZNF652 | no | 178 | 92.5926 | yes | -0.1074 | -0.0075 | 0.0759 | 1.0000 | 0.1424 | -0.0763 | 0.6646 |
| TLE1 | no | 149 | 84.2593 | no | -0.0894 | 0.4371 | 0.0631 | 1.0000 | 0.1682 | 0.1137 | 0.5085 |
| ZHX3 | no | 169 | 90.7407 | yes | -0.0707 | 0.1032 | 0.1147 | 1.0000 | 0.1465 | -0.0765 | -0.1227 |

#### Called under the amended criterion in S2A_k100, all 1,836

| factor | in published set | cells | cell pctl | high_depth_candidate | delta | depth floor | p | q | guide agreement | identity drop | ruler score change |
|---|---|---|---|---|---|---|---|---|---|---|---|
| GSC2 | no | 166 | 89.9237 | no | -0.9874 | 0.4601 | 0.0006 | 0.5692 | 0.1500 | 0.0028 | -0.5831 |
| RELA | no | 179 | 92.9194 | yes | -0.9867 | 0.3233 | 0.0006 | 0.5692 | 0.1444 | -0.2269 | 0.5372 |
| AFF1 | yes | 177 | 92.3747 | yes | -0.5266 | 0.1430 | 0.0031 | 1.0000 | 0.1441 | -0.0457 | -0.1530 |
| CRTC1 | no | 190 | 95.2070 | yes | -0.4256 | 0.6096 | 0.0039 | 1.0000 | 0.1383 | 0.0185 | 0.2843 |
| DMRTA1 | no | 183 | 94.0087 | yes | -0.2674 | 0.6715 | 0.0103 | 1.0000 | 0.1425 | -0.0142 | -0.0561 |
| HOXD1 | no | 181 | 93.5185 | yes | -0.1126 | 0.3329 | 0.0174 | 1.0000 | 0.1684 | 0.0437 | -0.0917 |

## 7. Limitations

- D1's bound is **not** a biological effect size in aged cells. The screen uses neonatal foreskin fibroblasts with no aged donor and no matched young target. The bound says how large an apparent approach could be hiding in this screen as measured against a GTEx-anchored young reference.
- No factor named here is a candidate for a pilot panel, and nothing in this file should be read as proposing one.
- D2 tests the permutation null on draws built from the same 14,013 non-targeting cells that define the depth-matched null. A null tested on its own construction's draws is the right test for exchangeability, but it does not test either null against an independent population.
- D3 matches on total UMI over the matched genes. Total reads and cell count are not the only depth variables: per-cell complexity, guide composition and capture efficiency are not matched.
- The UMI-matched null shares one set of draws across the factors in each 1% bin of total UMI rather than redrawing per factor; the realised offsets are reported above so the approximation is visible.
- The amended criterion still compares each factor with a null built by resampling one pool of non-targeting cells from one line, which carries no guide-to-guide variation of the kind a factor's ~6 guides carry.
- O2 pairs were not run. SOUTH6's recurring set came entirely from them.
- 3,221 of the ruler's 23,485 genes are outside the panel and are assumed unmoved.

## 8. Files

`src/south8_run.py`. In `results/south8/`: `PREREG.flag`, `stage0_reproduction.csv`, `calls_amended.csv`, `d1_plant_curves.csv`, `d1_largest_apparent.csv`, `d2_null_vectors.npz`, `d2_draw_meta.csv`, `d2_true_null_p.csv`, `d2_calibration.csv`, `d2_verdicts.csv`, `d3_factor_umi.csv`, `d3a_correlations.csv`, `d3b_umi_floors.csv`, `d3b_factor_tolerance.csv`, `d3c_*_umi_calls.csv`, `d3c_counts.csv`, `d4_recurring11.csv`, `cell_umi.npy`, `south8_run.log`.

`FINDINGS_SOUTH7.md` and every file in `results/south7/` are unchanged.

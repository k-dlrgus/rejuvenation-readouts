# FINDINGS_SAME — same-platform young target

**Status:** DONE. key=`no_approach_same_platform`. Seed `20260914`. boot `20260918`. n_perm=200. n_boot=200. P pass=True. Frozen ruler exists=True.

Does aged partially-reprogrammed tissue move toward the YOUNG donor's untreated day-0 fibroblasts from the same experiment? And is that rejuvenation, or just both donors converging?

Does not modify any existing `FINDINGS_*.md`, `PROGRESS_*.md`, or `FALSIFICATION.md`. No GSE325735. Uncalibrated R² is never a gate. Frozen ruler is a context column only. Donors never pooled. States never averaged.

Reproduced by `src/same_run.py`. Reuses functions from `src/toward_run.py` (not modified). Frozen ruler: `<repo>\results\fibro\frozen_ruler_ridge_raw.npz` exists=True.

## Pre-registration (verbatim, written before any SAME statistic)

Written **before any SAME statistic**, 2026-09-21. No threshold in this block is re-tuned after numbers exist.

Question: does aged partially-reprogrammed tissue move toward the YOUNG donor's untreated day-0 fibroblasts from the same experiment? And is that rejuvenation, or just both donors converging?

Seeds: 20260914 (splits, permutations), 20260918 (bootstrap). n_perm=200, n_boot=200. GSE325735 out of scope. Uncalibrated R² never a gate.

INPUTS (read from disk; refit/re-label/re-cluster/rescore nothing)
- GSE297234 per-cell counts, md2 labels (results/md2/t2_cluster_labels_*.csv, results/md2/louvain_obs_*.csv) exactly as used in TOWARD.
- Frozen-ruler overlap gene set and frozen GTEx mu/sd (same gene space and transform as TOWARD: pseudobulk -> TMM log2-CPM prior.count=2 across this task's full pseudobulk panel -> frozen GTEx mu/sd -> z; missing genes z=0).
- Frozen ruler: context column only. Not the statistic.

SPLIT (frozen before any statistic)
Randomly split each donor's d0 Fibroblast cells 50/50 (seed 20260914):
half A = anchors, half B = held out for the positive control only.
  O  = GM00731 d0 Fibroblast, half A   (aged origin)
  Y  = GM23815 d0 Fibroblast, half A   (young target)
Report n cells in every half. If any half has <500 cells, STOP.

PRIMARY STATISTICS (for a destination state S of one donor)
Forward test (aged -> young):
  v      = z(Y) - z(O)
  d_S    = z(S) - z(O),       S = GM00731 PartialReprog (all timepoints pooled)
  cos_S  = dot(d_S, v) / (||d_S|| * ||v||)
  frac_S = dot(d_S, v) / ||v||^2      (0 = no progress, 1 = reached young
                                        along the axis)
  delta_S = ||z(S) - z(Y)|| - ||z(O) - z(Y)||   (negative = closer to young)

Reverse test (young -> aged), identical formulas with donors swapped:
  v_rev = z(O) - z(Y); origin Y; S_rev = GM23815 PartialReprog (all timepoints
  pooled). Produces cos_rev, frac_rev, delta_rev.

Asymmetry = frac_S - frac_rev.

NULLS / CIs
- N1: 200 permutations of v's entries across genes (seed 20260914).
  p = (#{cos >= observed} + 1)/201. Same for the reverse test with v_rev.
- Bootstrap (seed 20260918, B=200): resample cells with replacement within
  each pseudobulk (O, Y, S, S_rev) independently; rebuild; 95% percentile CIs
  on frac_S, frac_rev, delta_S, delta_rev, and Asymmetry (recomputed on the
  same bootstrap draw).

POSITIVE CONTROL P (run BEFORE any PartialReprog vector)
Using ONLY half-B cells: build synthetic pseudobulks by mixing aged half-B d0
cells with young half-B d0 cells at young fractions f = 0, 0.10, 0.25, 0.50
(cell counts, not weights; seed 20260914; total cells fixed = min available).
Score each mixture as S in the forward test against the half-A anchors.
P PASSES only if: frac is monotonically increasing in f, AND frac at f=0.50
has bootstrap CI excluding 0, AND delta at f=0.50 < 0.
Also report the noise floor: frac and delta for f=0 (should be ~0).
If P fails: STOP. The statistic cannot detect a known approach. Report nothing
downstream.

CONTEXT (reported, never the reading)
- GM00731 Pluripotency as S in the forward test (overshoot reference).
- GM00731 NonReprog as S in the forward test.
- Frozen-ruler score of every pseudobulk.
- cos(v, u_TOWARD) where u_TOWARD is the GTEx young-minus-old axis from
  results/toward/anchors.npz: how much of the same-platform donor gap looks
  like age as GTEx defines it. Report only.

PRE-REGISTERED READING (decision order; only what fired)
1. P fails -> key `pos_control_broken`. STOP.
2. Forward: frac_S CI excludes 0 (positive) AND delta_S < 0 AND N1 p <= 0.05,
   AND Asymmetry CI excludes 0 (positive) -> key `toward_young_same_platform`.
   Aged partially-reprogrammed cells approach the young donor's untreated
   cells more than the reverse. State n_cells and frac_S in the same sentence.
   Also state: one young donor, so "young" is confounded with that donor.
3. Forward meets the conditions in 2 BUT Asymmetry CI includes 0 or is
   negative -> key `convergence_not_rejuvenation`. Both donors drift toward
   each other; reprogramming erases donor differences rather than age.
4. delta_S >= 0 OR frac_S CI includes 0 -> key `no_approach_same_platform`.
   The TOWARD negative holds with the platform gap removed.
5. Anything else -> key `mixed`; report exactly which conditions held.
If GM00731 Pluripotency also has delta < 0 and frac CI > 0, add the flag
`overshoot_also_closer` to whatever key fired and say the forward result may
reflect generic de-differentiation, not approach to young fibroblast.

WHAT DOES NOT COUNT
- The frozen ruler moving. Context only.
- The forward test alone. Without the reverse test, convergence and
  rejuvenation look identical. Asymmetry is required.
- Using half-B cells anywhere except P. Using half-A cells in P.
- Per-timepoint slicing to find a better result. PartialReprog is pooled
  across timepoints; report timepoint composition, do not split it.
- Relaxing the 500-cell split minimum, changing f values, or re-choosing
  anchors after seeing numbers.
- Averaging donors, directions, or states.

Flag: `results/same/PREREG.flag` exists=True.

## STOP / failures

Not substituting columns or repairing rows.

None recorded.

## Split sizes (frozen before any statistic)

Each donor's d0 Fibroblast cells split 50/50 with an independent `Generator(20260914)`. Half A = anchors (O, Y). Half B = positive control P only. STOP if any half has <500 cells.

| donor | half | n_cells | role | used_in |
|---|---|---|---|---|
| GM00731 | A | 2469 | O (aged origin; forward-test anchor) | anchors |
| GM00731 | B | 2470 | held out; P mixtures only | P |
| GM23815 | A | 3851 | Y (young target; forward-test anchor) | anchors |
| GM23815 | B | 3852 | held out; P mixtures only | P |

## Positive control P (run BEFORE any PartialReprog vector)

Mixtures of aged half-B d0 Fibroblast cells with young half-B d0 Fibroblast cells at young fractions f=[0.0, 0.1, 0.25, 0.5]; cell counts, not weights; total cells fixed = min(n_aged_B, n_young_B); seed `20260914`. Scored as S in the forward test against half-A anchors. P PASSES only if frac is monotonically increasing in f, AND frac at f=0.50 has bootstrap CI excluding 0, AND delta at f=0.50 < 0.

- P pass=True. monotonic=True. frac(f=0.50) CI excludes 0 = True [+0.438, +0.463]. delta(f=0.50)=-34.443 < 0 = True.
- Noise floor f=0: frac=+0.115 delta=+0.068 (should be ~0). n_total=2470.

| f | n_aged | n_young | n_total | frac | delta | cos | frac_ci_lo | frac_ci_hi | delta_ci_lo | delta_ci_hi | noise_floor |
|---|---|---|---|---|---|---|---|---|---|---|---|
| +0.000 | 2470 | 0 | 2470 | +0.115 | +0.068 | +0.239 | +0.168 | +0.195 | -2.060 | +1.972 | True |
| +0.100 | 2223 | 247 | 2470 | +0.182 | -7.870 | +0.371 | +0.225 | +0.252 | -9.312 | -5.248 | False |
| +0.250 | 1852 | 618 | 2470 | +0.283 | -18.913 | +0.530 | +0.309 | +0.335 | -18.807 | -14.793 | False |
| +0.500 | 1235 | 1235 | 2470 | +0.441 | -34.443 | +0.694 | +0.438 | +0.463 | -31.045 | -26.909 | False |

## Timepoint composition (PartialReprog pooled; not split)

Labels: md2 argmax of mean AddModuleScore (mmc3 Reprog_cell_state_signatures), exactly as used in TOWARD. Not re-labelled. PartialReprog is pooled across timepoints; composition is reported and is not used to slice a better result. Half-B d0 Fibroblast cells are not in any destination vector.

| cell_line | state | n_cells | n_d0 | n_d3 | n_d7 | n_d10 |
|---|---|---|---|---|---|---|
| GM00731 | Fibroblast | 5988 | 4939 | 5 | 14 | 1030 |
| GM00731 | PartialReprog | 5832 | 3 | 5220 | 578 | 31 |
| GM00731 | Pluripotency | 439 | 0 | 5 | 354 | 80 |
| GM00731 | NonReprog | 5459 | 65 | 348 | 1932 | 3114 |
| GM23815 | Fibroblast | 11355 | 7703 | 46 | 2233 | 1373 |
| GM23815 | PartialReprog | 6481 | 9 | 6113 | 309 | 50 |

## Forward, reverse, and asymmetry

Forward: v = z(Y)−z(O), S = GM00731 PartialReprog (timepoints pooled). Reverse: donors swapped; S_rev = GM23815 PartialReprog (timepoints pooled). Asymmetry = frac_S − frac_rev. N1: 200 permutations of v (resp. v_rev) entries (seed 20260914); p = (#{cos ≥ observed} + 1)/201. Bootstrap B=200 seed 20260918; 95% percentile CIs. Uncalibrated R² is never a gate.

| test | origin | target | S | n_cells_S | cos | frac | delta | p_N1 | frac_ci_lo | frac_ci_hi | delta_ci_lo | delta_ci_hi | d_norm | v_norm |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| forward | GM00731 d0 Fibroblast half A | GM23815 d0 Fibroblast half A | GM00731 PartialReprog pooled | 5832 | +0.323 | +0.674 | +130.894 | 0.0050 | +0.624 | +0.666 | +119.669 | +125.224 | +272.734 | +130.688 |
| reverse | GM23815 d0 Fibroblast half A | GM00731 d0 Fibroblast half A | GM23815 PartialReprog pooled | 6481 | -0.029 | -0.057 | +163.342 | 1.0000 | +0.026 | +0.063 | +152.488 | +156.719 | +259.649 | +130.688 |

- Asymmetry = frac_S − frac_rev = +0.732 CI [+0.562, +0.640].

## Context (reported, never the reading)

GM00731 Pluripotency and NonReprog scored as S in the forward test. Frozen-ruler score of every pseudobulk (Z @ w). cos(v, u_TOWARD) from `results/toward/anchors.npz` — how much of the same-platform donor gap looks like age as GTEx defines it. Report only.

| pseudobulk | cell_line | state | n_cells | role | ruler_score_context | frac | delta | cos | frac_ci_lo | frac_ci_hi | delta_ci_lo | delta_ci_hi |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| O | GM00731 | Fibroblast_d0_halfA | 2469 | origin_forward | +2.647 | NA | NA | NA | NA | NA | NA | NA |
| Y | GM23815 | Fibroblast_d0_halfA | 3851 | target_forward | -1.934 | NA | NA | NA | NA | NA | NA | NA |
| S | GM00731 | PartialReprog | 5832 | forward_S | +0.929 | +0.674 | +130.894 | +0.323 | +0.624 | +0.666 | +119.669 | +125.224 |
| S_rev | GM23815 | PartialReprog | 6481 | reverse_S | -4.350 | -0.057 | +163.342 | -0.029 | +0.026 | +0.063 | +152.488 | +156.719 |
| Pluri | GM00731 | Pluripotency | 439 | overshoot_reference_forward | -10.041 | +1.403 | +382.106 | +0.338 | +1.165 | +1.294 | +379.979 | +406.190 |
| NonReprog | GM00731 | NonReprog | 5459 | off_trajectory_forward | +7.089 | +0.555 | +60.331 | +0.370 | +0.532 | +0.561 | +52.577 | +56.671 |

- cos(v, u_TOWARD)=+0.039 (report only; not a gate).

## Pre-registered reading (only the outcome that fired)

**key: `no_approach_same_platform`.**

The TOWARD negative holds with the platform gap removed. GM00731 PartialReprog n_cells=5832, frac_S=+0.674, delta_S=+130.894 (delta_S >= 0: not closer to the young donor).

## Limitations

- One young donor, so "young" is confounded with that donor.
- Day-0 culture/passage differences between donors are part of v.
- Two donors total. Donors are not averaged.
- State labels are Louvain-cluster argmax of signature scores, not per-cell argmax.
- TMM is among a handful of pseudobulk rows (this task's panel).
- Bootstrap CIs re-TMM the panel, so they are CIs on the TMM-coupled statistic and need not contain the observed point estimate. Not a gate.
- Half-A anchors are a random 50% of d0 Fibroblast cells; the complementary half is used only in P.
- Uncalibrated R² is never a gate and is not computed here.
- GSE325735 is out of scope and is not opened.

## Files touched

- `src/same_run.py`
- `FINDINGS_SAME.md`
- `PROGRESS_SAME.md`
- `results/same/PREREG.flag` (not rewritten if it already existed)
- `results/same/split.json`
- `results/same/split_table.csv`
- `results/same/split_idx.npz`
- `results/same/posctrl.csv`
- `results/same/posctrl_summary.json`
- `results/same/occupancy.csv`
- `results/same/stats.csv`
- `results/same/context.csv`
- `results/same/reading.json`
- `results/same/summary.json`
- `results/same/manifest.json`
- `results/same/run_report.txt`
- Existing `FINDINGS_*.md`, `PROGRESS_*.md`, and `FALSIFICATION.md` were not modified (except FINDINGS_SAME.md and PROGRESS_SAME.md).


# FINDINGS_PLANE — does the age–identity curve bend?

**Status:** Task 1: GM00731 frozen ruler=`mixed`; MD=`mixed`; age-up − age-down=`bend`. GM23815 frozen ruler=`mixed`; MD=`mixed`; age-up − age-down=`bend`. Task 2: report-only extrapolation tabulated; no gate. Task 3: one figure per donor. Seed `20260914`. boot `20260918`. n_perm=200. n_boot=200. n_random=200. Frozen ruler `frozen_ruler_ridge_raw.npz` exists=True.

Does not modify any existing `FINDINGS_*.md` or `FALSIFICATION.md`. No GSE325735. Uncalibrated R² is reported and is never a gate. Nothing averaged across donors or instruments. Refit nothing. d7 is not the primary Stage 2 endpoint. Task 2 is report-only and is not used to discount Task 1.

Reproduced by `src/plane_run.py`. Frozen ruler: `<repo>\results\fibro\frozen_ruler_ridge_raw.npz` exists=True.

Flag: `PREREG_TASK1.flag` exists=True. `PREREG_TASK2.flag` exists=True. `PREREG_TASK3.flag` exists=True. `DECLARED_BEFORE_SCORES.flag` exists=True.

## The question

Place reprogramming cells on the age–identity plane and ask whether the curve bends — whether age can move before identity is lost — or whether age and identity move in lockstep. FINDINGS_MD5.md established that partially reprogrammed cells read younger than their starting fibroblasts on the frozen ruler. The live alternative is that the ruler reads younger because cells are leaving fibroblast identity, not because they are rejuvenating.

- Lockstep (identity loss): the age reading keeps falling as cells move further from fibroblast, monotonically, all the way to pluripotency.
- Bend (rejuvenation separable from dedifferentiation): the age reading falls from Fibroblast to PartialReprog, then stops falling or reverses as cells continue toward EarlyPluripotency and Pluripotency, while identity keeps moving.

A bend is the result the project was built to find. Lockstep is equally reportable and is not argued away. Trajectory order, fixed before scoring: Fibroblast → PartialReprog → EarlyPluripotency → Pluripotency. `NonReprog` is off-trajectory: scored and plotted, excluded from every trajectory statistic. Unit: cell state per donor; timepoints pooled within state. ≥30 cells per state per donor or the state is not scored. x = identity = `pluri_primary` (FINDINGS_FIBRO.md frozen lists; md4 primary, drop_oskm=False). y = age. Three instruments, never combined: frozen ruler; MD; Table S3 age-up minus age-down.

## Failures recorded

None recorded.

## Task 1 pre-registration (verbatim, written before any PLANE state-median, slope, bend, or null p)

Written **before any PLANE state-median, slope, bend statistic, or null p-value**, 2026-09-21. No threshold in this block is re-tuned after numbers exist.

Using the src/md2_*.py Louvain clustering already computed (not re-run, not re-tuned). Cell-state labels are the md2 assignment: argmax of mean AddModuleScore of mmc3 Reprog_cell_state_signatures columns Fibroblast / PartialReprog / EarlyPluripotency / Pluripotency / NonReprog. Labels are read from results/md2/t2_cluster_labels_*.csv. Not re-labelled. AddModuleScore for MD is the per-cell vector already stored in data/processed/md2/louvain_ams_*.npz (the md2 computation). Table S3 age-up minus age-down and the frozen-ruler / pluripotency−fibroblast per-cell scores are the vectors already stored in data/processed/md4/cell_scores_*.npz (the md4 computation). Frozen ruler: results/fibro/frozen_ruler_ridge_raw.npz used as-is. Refit nothing. Do not re-cluster, re-label, re-tune, rescore observed instruments, or re-order the trajectory after seeing scores.

Trajectory order, fixed before scoring: Fibroblast → PartialReprog → EarlyPluri → Pluripotency. NonReprog is off-trajectory: score and plot it, but exclude it from every trajectory statistic.

Unit of analysis: the cell state, per donor. Pool timepoints within each state (the state, not the day, is the trajectory position here). Report the timepoint composition of every state. Require ≥30 cells per state per donor; a state below that is reported as not scored, not dropped silently.

Two coordinates per cell: x = identity = pluripotency−fibroblast score (frozen lists from FINDINGS_FIBRO.md; primary version as stored in md4 cell_scores, drop_oskm=False). y = age. Three age instruments, plotted and tested separately, never combined: frozen ruler; MD; Table S3 age-up minus age-down.

Per donor, per instrument:

1. State medians and bootstrap 95% CIs (B=200, cells resampled within state) on both coordinates.
2. Segment slopes. Δage / Δidentity for segment A (Fibroblast → PartialReprog) and segment B (PartialReprog → Pluripotency; if Pluripotency is not scored, use EarlyPluri and say so).
3. Bend statistic = slope_A − slope_B, with bootstrap CI.
4. Null: for the frozen ruler, the same statistic on 200 permuted-weight random directions (seed 20260914). For each gene-list score (MD, age-up−age-down), 200 size-matched, expression-bin-matched random gene sets. Report the empirical p for the observed bend being at least as large as the null's.

Young direction = destination state median age below origin (Δage < 0). If Δidentity > 0, that is a negative slope; if Δidentity < 0, a positive slope. Slope CI excludes 0 in the young direction iff the bootstrap CI of the slope lies entirely on that side of 0. Segment B is flat if its slope CI includes 0; reversed (old direction) if the slope CI excludes 0 on the age-rises side (Δage > 0).

Empirical p for the bend: p_more_negative = (n_null ≤ obs + 1) / (n_null + 1). That is the one-sided test that the observed bend is at least as large as the null's in the rejuvenation-separable direction (slope_A more negative than slope_B when identity increases). Also report p_more_positive. “Beats its null” = p_more_negative ≤ 0.05. Seed 20260914 for random directions and random gene sets. Bootstrap seed 20260918 with a documented per-(donor, instrument) offset. Donors never pooled. Instruments never combined. NonReprog never in a trajectory statistic.

Sanity, before any new median or bend: (1) reproduce FINDINGS_MD2.md pooled MD means on the aged donor (PartialReprog +0.150, NonReprog +0.439; t2_reproduction.json); (2) occupancy matches results/md4/t1_cell_counts.csv summed within state. If either fails, stop. Cell-state labels have drifted.

**Pre-registered reading (per donor, per instrument; only what fired):**
- **Bend:** age falls in segment A (slope_A CI excludes 0 in the young direction) **and** segment B's age change is flat or reversed (slope_B CI includes 0 or is in the old direction) **and** the bend statistic beats its null at p ≤ 0.05 → age moves before identity is lost; rejuvenation is separable from dedifferentiation on this instrument.
- **Lockstep:** age falls in both segments with overlapping slopes and the bend statistic does not beat its null → this instrument cannot separate rejuvenation from identity loss.
- **Mixed / other:** report exactly what fired.

Instrument comparison, report-only: state which instruments bend and which run in lockstep. FINDINGS_MD3.md Task 3 found MD correlates ρ ≈ −0.94 with the pluripotency−fibroblast score, so MD is expected to run in lockstep; say whether it does. If the frozen ruler bends and MD does not, that is the separation between the instruments, stated as observed, with no winner declared.


## Task 1 — state cell counts and timepoint composition (md2 Louvain labels, not re-labelled)

Assignment basis: argmax of mean AddModuleScore of mmc3 `Reprog_cell_state_signatures` (Fibroblast, PartialReprog, EarlyPluripotency, Pluripotency, NonReprog). Labels from `results/md2/t2_cluster_labels_*.csv`. Not re-labelled. Occupancy matches `results/md4/t1_cell_counts.csv` or this file stops.

Pooled n cells by state (timepoints pooled within state for scoring; composition reported):

| cell_line | label | n_clusters | n_d0 | n_d3 | n_d7 | n_d10 | n_all | ge30 | scored | reason |
|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | Fibroblast | 4 | 4939 | 5 | 14 | 1030 | 5988 | True | True | NA |
| GM00731 | PartialReprog | 6 | 3 | 5220 | 578 | 31 | 5832 | True | True | NA |
| GM00731 | EarlyPluripotency | 4 | 2 | 410 | 3917 | 85 | 4414 | True | True | NA |
| GM00731 | Pluripotency | 1 | 0 | 5 | 354 | 80 | 439 | True | True | NA |
| GM00731 | NonReprog | 4 | 65 | 348 | 1932 | 3114 | 5459 | True | True | NA |
| GM23815 | Fibroblast | 5 | 7703 | 46 | 2233 | 1373 | 11355 | True | True | NA |
| GM23815 | PartialReprog | 4 | 9 | 6113 | 309 | 50 | 6481 | True | True | NA |
| GM23815 | EarlyPluripotency | 8 | 56 | 316 | 8172 | 4417 | 12961 | True | True | NA |
| GM23815 | Pluripotency | 1 | 0 | 91 | 977 | 181 | 1249 | True | True | NA |
| GM23815 | NonReprog | 1 | 0 | 10 | 8 | 1794 | 1812 | True | True | NA |

Timepoint composition of every state:

| cell_line | day | label | n_cells | n_timepoint | n_state | frac_state | frac_timepoint | state_ge30 |
|---|---|---|---|---|---|---|---|---|
| GM00731 | 0 | Fibroblast | 4939 | 5009 | 5988 | +0.825 | +0.986 | True |
| GM00731 | 0 | PartialReprog | 3 | 5009 | 5832 | +0.001 | +0.001 | True |
| GM00731 | 0 | EarlyPluripotency | 2 | 5009 | 4414 | +0.000 | +0.000 | True |
| GM00731 | 0 | Pluripotency | 0 | 5009 | 439 | +0.000 | +0.000 | True |
| GM00731 | 0 | NonReprog | 65 | 5009 | 5459 | +0.012 | +0.013 | True |
| GM00731 | 3 | Fibroblast | 5 | 5988 | 5988 | +0.001 | +0.001 | True |
| GM00731 | 3 | PartialReprog | 5220 | 5988 | 5832 | +0.895 | +0.872 | True |
| GM00731 | 3 | EarlyPluripotency | 410 | 5988 | 4414 | +0.093 | +0.068 | True |
| GM00731 | 3 | Pluripotency | 5 | 5988 | 439 | +0.011 | +0.001 | True |
| GM00731 | 3 | NonReprog | 348 | 5988 | 5459 | +0.064 | +0.058 | True |
| GM00731 | 7 | Fibroblast | 14 | 6795 | 5988 | +0.002 | +0.002 | True |
| GM00731 | 7 | PartialReprog | 578 | 6795 | 5832 | +0.099 | +0.085 | True |
| GM00731 | 7 | EarlyPluripotency | 3917 | 6795 | 4414 | +0.887 | +0.576 | True |
| GM00731 | 7 | Pluripotency | 354 | 6795 | 439 | +0.806 | +0.052 | True |
| GM00731 | 7 | NonReprog | 1932 | 6795 | 5459 | +0.354 | +0.284 | True |
| GM00731 | 10 | Fibroblast | 1030 | 4340 | 5988 | +0.172 | +0.237 | True |
| GM00731 | 10 | PartialReprog | 31 | 4340 | 5832 | +0.005 | +0.007 | True |
| GM00731 | 10 | EarlyPluripotency | 85 | 4340 | 4414 | +0.019 | +0.020 | True |
| GM00731 | 10 | Pluripotency | 80 | 4340 | 439 | +0.182 | +0.018 | True |
| GM00731 | 10 | NonReprog | 3114 | 4340 | 5459 | +0.570 | +0.718 | True |
| GM23815 | 0 | Fibroblast | 7703 | 7768 | 11355 | +0.678 | +0.992 | True |
| GM23815 | 0 | PartialReprog | 9 | 7768 | 6481 | +0.001 | +0.001 | True |
| GM23815 | 0 | EarlyPluripotency | 56 | 7768 | 12961 | +0.004 | +0.007 | True |
| GM23815 | 0 | Pluripotency | 0 | 7768 | 1249 | +0.000 | +0.000 | True |
| GM23815 | 0 | NonReprog | 0 | 7768 | 1812 | +0.000 | +0.000 | True |
| GM23815 | 3 | Fibroblast | 46 | 6576 | 11355 | +0.004 | +0.007 | True |
| GM23815 | 3 | PartialReprog | 6113 | 6576 | 6481 | +0.943 | +0.930 | True |
| GM23815 | 3 | EarlyPluripotency | 316 | 6576 | 12961 | +0.024 | +0.048 | True |
| GM23815 | 3 | Pluripotency | 91 | 6576 | 1249 | +0.073 | +0.014 | True |
| GM23815 | 3 | NonReprog | 10 | 6576 | 1812 | +0.006 | +0.002 | True |
| GM23815 | 7 | Fibroblast | 2233 | 11699 | 11355 | +0.197 | +0.191 | True |
| GM23815 | 7 | PartialReprog | 309 | 11699 | 6481 | +0.048 | +0.026 | True |
| GM23815 | 7 | EarlyPluripotency | 8172 | 11699 | 12961 | +0.631 | +0.699 | True |
| GM23815 | 7 | Pluripotency | 977 | 11699 | 1249 | +0.782 | +0.084 | True |
| GM23815 | 7 | NonReprog | 8 | 11699 | 1812 | +0.004 | +0.001 | True |
| GM23815 | 10 | Fibroblast | 1373 | 7815 | 11355 | +0.121 | +0.176 | True |
| GM23815 | 10 | PartialReprog | 50 | 7815 | 6481 | +0.008 | +0.006 | True |
| GM23815 | 10 | EarlyPluripotency | 4417 | 7815 | 12961 | +0.341 | +0.565 | True |
| GM23815 | 10 | Pluripotency | 181 | 7815 | 1249 | +0.145 | +0.023 | True |
| GM23815 | 10 | NonReprog | 1794 | 7815 | 1812 | +0.990 | +0.230 | True |

Qualify (segment A = Fibroblast → PartialReprog; segment B = PartialReprog → Pluripotency, else EarlyPluripotency and said so):

| cell_line | n_Fibroblast | n_PartialReprog | n_EarlyPluripotency | n_Pluripotency | n_NonReprog | segment_A_ok | segment_B_ok | segment_B_end | segment_B_end_is_fallback | NonReprog_scored | min_cells | reason |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | 5988 | 5832 | 4414 | 439 | 5459 | True | True | Pluripotency | False | True | 30 | NA |
| GM23815 | 11355 | 6481 | 12961 | 1249 | 1812 | True | True | Pluripotency | False | True | 30 | NA |

## Task 1 — sanity check (FINDINGS_MD2.md pooled MD means; occupancy vs md4)

- n_PartialReprog=5832 n_NonReprog=5459 MD_PR=+0.150471 MD_NR=+0.439448 ok=True
- md2 json n_PR=5832 n_NR=5459 MD_PR=+0.150471 MD_NR=+0.439448
- printed 3 d.p. match (+0.150 / +0.439): True

## Task 1 — state medians with bootstrap 95% CIs (B=200, cells within state)

Identity = pluripotency−fibroblast (`pluri_primary`). Age instruments in separate rows. NonReprog is scored and plotted; it is not used in slopes or the bend statistic.

| cell_line | instrument | label | n_cells | scored | on_trajectory | identity_median | identity_ci_lo | identity_ci_hi | age_median | age_ci_lo | age_ci_hi |
|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | frozen_ruler | Fibroblast | 5988 | True | True | +8.281 | +8.266 | +8.295 | +18.229 | +18.216 | +18.251 |
| GM00731 | frozen_ruler | PartialReprog | 5832 | True | True | +12.135 | +12.067 | +12.217 | +17.490 | +17.464 | +17.521 |
| GM00731 | frozen_ruler | EarlyPluripotency | 4414 | True | True | +12.165 | +12.106 | +12.209 | +17.041 | +17.005 | +17.064 |
| GM00731 | frozen_ruler | Pluripotency | 439 | True | True | +12.986 | +12.895 | +13.101 | +15.349 | +15.225 | +15.421 |
| GM00731 | frozen_ruler | NonReprog | 5459 | True | False | +8.908 | +8.884 | +8.929 | +18.386 | +18.360 | +18.412 |
| GM00731 | md_score | Fibroblast | 5988 | True | True | +8.281 | +8.267 | +8.298 | +0.505 | +0.503 | +0.508 |
| GM00731 | md_score | PartialReprog | 5832 | True | True | +12.135 | +12.063 | +12.223 | +0.128 | +0.122 | +0.134 |
| GM00731 | md_score | EarlyPluripotency | 4414 | True | True | +12.165 | +12.091 | +12.209 | +0.126 | +0.120 | +0.133 |
| GM00731 | md_score | Pluripotency | 439 | True | True | +12.986 | +12.888 | +13.101 | +0.004 | -0.003 | +0.018 |
| GM00731 | md_score | NonReprog | 5459 | True | False | +8.908 | +8.887 | +8.931 | +0.442 | +0.439 | +0.445 |
| GM00731 | age_up_minus_age_down | Fibroblast | 5988 | True | True | +8.281 | +8.266 | +8.296 | +0.025 | +0.024 | +0.027 |
| GM00731 | age_up_minus_age_down | PartialReprog | 5832 | True | True | +12.135 | +12.067 | +12.220 | -0.101 | -0.103 | -0.099 |
| GM00731 | age_up_minus_age_down | EarlyPluripotency | 4414 | True | True | +12.165 | +12.095 | +12.206 | -0.000 | -0.002 | +0.002 |
| GM00731 | age_up_minus_age_down | Pluripotency | 439 | True | True | +12.986 | +12.873 | +13.100 | +0.073 | +0.068 | +0.082 |
| GM00731 | age_up_minus_age_down | NonReprog | 5459 | True | False | +8.908 | +8.888 | +8.931 | +0.029 | +0.028 | +0.031 |
| GM23815 | frozen_ruler | Fibroblast | 11355 | True | True | +9.442 | +9.432 | +9.451 | +18.854 | +18.840 | +18.867 |
| GM23815 | frozen_ruler | PartialReprog | 6481 | True | True | +12.087 | +12.040 | +12.125 | +18.198 | +18.176 | +18.214 |
| GM23815 | frozen_ruler | EarlyPluripotency | 12961 | True | True | +12.218 | +12.183 | +12.248 | +18.130 | +18.114 | +18.144 |
| GM23815 | frozen_ruler | Pluripotency | 1249 | True | True | +13.684 | +13.578 | +13.751 | +17.095 | +17.029 | +17.165 |
| GM23815 | frozen_ruler | NonReprog | 1812 | True | False | +9.343 | +9.312 | +9.362 | +19.188 | +19.166 | +19.215 |
| GM23815 | md_score | Fibroblast | 11355 | True | True | +9.442 | +9.432 | +9.452 | +0.459 | +0.457 | +0.461 |
| GM23815 | md_score | PartialReprog | 6481 | True | True | +12.087 | +12.044 | +12.127 | +0.075 | +0.069 | +0.081 |
| GM23815 | md_score | EarlyPluripotency | 12961 | True | True | +12.218 | +12.184 | +12.253 | +0.115 | +0.111 | +0.120 |
| GM23815 | md_score | Pluripotency | 1249 | True | True | +13.684 | +13.559 | +13.751 | -0.062 | -0.066 | -0.059 |
| GM23815 | md_score | NonReprog | 1812 | True | False | +9.343 | +9.322 | +9.362 | +0.422 | +0.416 | +0.426 |
| GM23815 | age_up_minus_age_down | Fibroblast | 11355 | True | True | +9.442 | +9.432 | +9.453 | +0.010 | +0.009 | +0.011 |
| GM23815 | age_up_minus_age_down | PartialReprog | 6481 | True | True | +12.087 | +12.047 | +12.128 | -0.118 | -0.120 | -0.117 |
| GM23815 | age_up_minus_age_down | EarlyPluripotency | 12961 | True | True | +12.218 | +12.186 | +12.254 | -0.029 | -0.030 | -0.027 |
| GM23815 | age_up_minus_age_down | Pluripotency | 1249 | True | True | +13.684 | +13.579 | +13.751 | +0.023 | +0.021 | +0.025 |
| GM23815 | age_up_minus_age_down | NonReprog | 1812 | True | False | +9.343 | +9.322 | +9.364 | +0.024 | +0.022 | +0.025 |

## Task 1 — segment slopes, bend statistic, nulls

slope = Δage / Δidentity. Bend = slope_A − slope_B. Frozen-ruler null = 200 permuted-weight random directions, seed 20260914. Gene-list nulls = 200 size-matched, expression-bin-matched random gene sets. Identity held fixed. p_more_negative = (n_null ≤ obs + 1) / (n_null + 1). Beats null iff p_more_negative ≤ 0.05.

| cell_line | instrument | state_B1 | segment_B_end_is_fallback | ok | delta_identity_A | delta_age_A | slope_A | slope_A_ci_lo | slope_A_ci_hi | delta_identity_B | delta_age_B | slope_B | slope_B_ci_lo | slope_B_ci_hi | bend | bend_ci_lo | bend_ci_hi | slope_A_young | slope_B_young | slope_B_flat | slope_B_old | slopes_overlap |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | frozen_ruler | Pluripotency | False | True | +3.853 | -0.740 | -0.192 | -0.200 | -0.182 | +0.851 | -2.140 | -2.514 | -3.038 | -2.202 | +2.322 | +2.006 | +2.853 | True | True | False | False | False |
| GM00731 | md_score | Pluripotency | False | True | +3.853 | -0.377 | -0.098 | -0.099 | -0.096 | +0.851 | -0.124 | -0.146 | -0.167 | -0.128 | +0.048 | +0.029 | +0.071 | True | True | False | False | False |
| GM00731 | age_up_minus_age_down | Pluripotency | False | True | +3.853 | -0.126 | -0.033 | -0.034 | -0.032 | +0.851 | +0.175 | +0.205 | +0.180 | +0.257 | -0.238 | -0.289 | -0.213 | True | False | False | True | False |
| GM23815 | frozen_ruler | Pluripotency | False | True | +2.646 | -0.657 | -0.248 | -0.259 | -0.240 | +1.597 | -1.103 | -0.691 | -0.755 | -0.630 | +0.443 | +0.377 | +0.504 | True | True | False | False | False |
| GM23815 | md_score | Pluripotency | False | True | +2.646 | -0.384 | -0.145 | -0.147 | -0.144 | +1.597 | -0.137 | -0.086 | -0.091 | -0.081 | -0.060 | -0.065 | -0.053 | True | True | False | False | False |
| GM23815 | age_up_minus_age_down | Pluripotency | False | True | +2.646 | -0.128 | -0.048 | -0.049 | -0.047 | +1.597 | +0.141 | +0.089 | +0.083 | +0.096 | -0.137 | -0.144 | -0.132 | True | False | False | True | False |

| cell_line | instrument | ok | null_kind | n_random | n_null_finite | obs_bend | null_bend_median | null_bend_ci_lo | null_bend_ci_hi | p_more_negative | p_more_positive | segment_B_end | segment_B_end_is_fallback |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | frozen_ruler | True | permuted_ruler_weights | 200 | 200 | +2.322 | +0.056 | -1.957 | +1.641 | +0.985 | +0.020 | Pluripotency | False |
| GM00731 | md_score | True | size_matched_random_gene_sets_mean_expression_bin | 200 | 200 | +0.048 | +0.147 | +0.065 | +0.230 | +0.010 | +0.995 | Pluripotency | False |
| GM00731 | age_up_minus_age_down | True | size_matched_random_gene_sets_mean_expression_bin_up_minus_down | 200 | 200 | -0.238 | +0.012 | -0.015 | +0.036 | +0.005 | +1.000 | Pluripotency | False |
| GM23815 | frozen_ruler | True | permuted_ruler_weights | 200 | 200 | +0.443 | +0.090 | -0.800 | +0.803 | +0.816 | +0.189 | Pluripotency | False |
| GM23815 | md_score | True | size_matched_random_gene_sets_mean_expression_bin | 200 | 200 | -0.060 | +0.086 | +0.042 | +0.133 | +0.005 | +1.000 | Pluripotency | False |
| GM23815 | age_up_minus_age_down | True | size_matched_random_gene_sets_mean_expression_bin_up_minus_down | 200 | 200 | -0.137 | +0.001 | -0.013 | +0.014 | +0.005 | +1.000 | Pluripotency | False |

| cell_line | instrument | key | ok | slope_A | slope_B | bend | p_more_negative | p_more_positive | bend_beats_null | slope_A_young | slope_B_young | slope_B_flat | slope_B_old | slopes_overlap |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | frozen_ruler | mixed | True | -0.192 | -2.514 | +2.322 | +0.985 | +0.020 | False | True | True | False | False | False |
| GM00731 | md_score | mixed | True | -0.098 | -0.146 | +0.048 | +0.010 | +0.995 | True | True | True | False | False | False |
| GM00731 | age_up_minus_age_down | bend | True | -0.033 | +0.205 | -0.238 | +0.005 | +1.000 | True | True | False | False | True | False |
| GM23815 | frozen_ruler | mixed | True | -0.248 | -0.691 | +0.443 | +0.816 | +0.189 | False | True | True | False | False | False |
| GM23815 | md_score | mixed | True | -0.145 | -0.086 | -0.060 | +0.005 | +1.000 | True | True | True | False | False | False |
| GM23815 | age_up_minus_age_down | bend | True | -0.048 | +0.089 | -0.137 | +0.005 | +1.000 | True | True | False | False | True | False |

## Task 1 — reading (only what fired)

### GM00731

**frozen ruler** key=`mixed` segment_B_end=`Pluripotency` fallback=False null=`permuted_ruler_weights`.

slope_A=-0.1920 slope_B=-2.5140 bend=+2.3221 p_more_negative=0.9851 p_more_positive=0.0199.

mixed / other: slope_A_young=True; slope_B_young=True; slope_B_flat=False; slope_B_old=False; slopes_overlap=False; bend_beats_null=False; age_falls_A=True; age_falls_B=True; p_more_negative=0.9850746268656716

flags: slope_A_young=True; slope_B_young=True; slope_B_flat=False; slope_B_old=False; slopes_overlap=False; bend_beats_null=False; age_falls_A=True; age_falls_B=True; ok=True

**MD** key=`mixed` segment_B_end=`Pluripotency` fallback=False null=`size_matched_random_gene_sets_mean_expression_bin`.

slope_A=-0.0979 slope_B=-0.1455 bend=+0.0477 p_more_negative=0.0100 p_more_positive=0.9950.

mixed / other: slope_A_young=True; slope_B_young=True; slope_B_flat=False; slope_B_old=False; slopes_overlap=False; bend_beats_null=True; age_falls_A=True; age_falls_B=True; p_more_negative=0.009950248756218905

flags: slope_A_young=True; slope_B_young=True; slope_B_flat=False; slope_B_old=False; slopes_overlap=False; bend_beats_null=True; age_falls_A=True; age_falls_B=True; ok=True

**age-up − age-down** key=`bend` segment_B_end=`Pluripotency` fallback=False null=`size_matched_random_gene_sets_mean_expression_bin_up_minus_down`.

slope_A=-0.0328 slope_B=+0.2050 bend=-0.2377 p_more_negative=0.0050 p_more_positive=1.0000.

age falls in segment A (slope_A CI excludes 0 in the young direction) and segment B is flat or reversed and the bend statistic beats its null. Age moves before identity is lost; rejuvenation is separable from dedifferentiation on this instrument.

flags: slope_A_young=True; slope_B_young=False; slope_B_flat=False; slope_B_old=True; slopes_overlap=False; bend_beats_null=True; age_falls_A=True; age_falls_B=False; ok=True

### GM23815

**frozen ruler** key=`mixed` segment_B_end=`Pluripotency` fallback=False null=`permuted_ruler_weights`.

slope_A=-0.2483 slope_B=-0.6909 bend=+0.4426 p_more_negative=0.8159 p_more_positive=0.1891.

mixed / other: slope_A_young=True; slope_B_young=True; slope_B_flat=False; slope_B_old=False; slopes_overlap=False; bend_beats_null=False; age_falls_A=True; age_falls_B=True; p_more_negative=0.8159203980099502

flags: slope_A_young=True; slope_B_young=True; slope_B_flat=False; slope_B_old=False; slopes_overlap=False; bend_beats_null=False; age_falls_A=True; age_falls_B=True; ok=True

**MD** key=`mixed` segment_B_end=`Pluripotency` fallback=False null=`size_matched_random_gene_sets_mean_expression_bin`.

slope_A=-0.1453 slope_B=-0.0857 bend=-0.0596 p_more_negative=0.0050 p_more_positive=1.0000.

mixed / other: slope_A_young=True; slope_B_young=True; slope_B_flat=False; slope_B_old=False; slopes_overlap=False; bend_beats_null=True; age_falls_A=True; age_falls_B=True; p_more_negative=0.004975124378109453

flags: slope_A_young=True; slope_B_young=True; slope_B_flat=False; slope_B_old=False; slopes_overlap=False; bend_beats_null=True; age_falls_A=True; age_falls_B=True; ok=True

**age-up − age-down** key=`bend` segment_B_end=`Pluripotency` fallback=False null=`size_matched_random_gene_sets_mean_expression_bin_up_minus_down`.

slope_A=-0.0483 slope_B=+0.0886 bend=-0.1369 p_more_negative=0.0050 p_more_positive=1.0000.

age falls in segment A (slope_A CI excludes 0 in the young direction) and segment B is flat or reversed and the bend statistic beats its null. Age moves before identity is lost; rejuvenation is separable from dedifferentiation on this instrument.

flags: slope_A_young=True; slope_B_young=False; slope_B_flat=False; slope_B_old=True; slopes_overlap=False; bend_beats_null=True; age_falls_A=True; age_falls_B=False; ok=True

t1_summary by_donor={'GM00731': {'frozen_ruler': 'mixed', 'md_score': 'mixed', 'age_up_minus_age_down': 'bend'}, 'GM23815': {'frozen_ruler': 'mixed', 'md_score': 'mixed', 'age_up_minus_age_down': 'bend'}} sanity_ok=True

## Instrument comparison (report-only)

Instruments are not combined. No winner is declared.

FINDINGS_MD3.md Task 3: MD vs pluri_primary ρ=-0.938 (GM00731 n=33) and ρ=-0.96 (GM23815 n=34). MD is expected to run in lockstep. No winner is declared.

- frozen ruler: GM00731=`mixed`; GM23815=`mixed`
- MD: GM00731=`mixed`; GM23815=`mixed`
- age-up − age-down: GM00731=`bend`; GM23815=`bend`

MD vs the lockstep definition: GM00731 key=`mixed` age_falls_A=True age_falls_B=True slopes_overlap=False bend_beats_null=True; GM23815 key=`mixed` age_falls_A=True age_falls_B=True slopes_overlap=False bend_beats_null=True.
MD did not fire `lockstep` on either donor. FINDINGS_MD3.md Task 3 ρ ≈ −0.94 predicted lockstep; that is not the key that fired. No winner is declared.

## Task 2 pre-registration (verbatim, written before any PLANE extrapolation distance)

Written **before any PLANE extrapolation distance**, 2026-09-21. Report-only. No gate.

Per state per donor, report the extrapolation distance to the GTEx fibroblast training distribution (same PCA k=50 nearest-neighbour Euclidean and Mahalanobis measures as FINDINGS_FIBRO2.md / src/fibro2_extrap.py). Sum Louvain-kept UMIs to a pseudobulk on the frozen-ruler gene space (md4 louvain_rulerY aligned). TMM among that donor's scored-state rows, then frozen μ/σ/w. Do not re-TMM among cells. Do not refit PCA. Do not rescore the per-cell instruments of Task 1.

State the direction extrapolation pushes the result, using the evidence already in the repo: FINDINGS_FIBRO3.md Task 2 (the same PCA k=50 NN / Mahalanobis measures as FINDINGS_FIBRO2.md) found the ruler score falls with distance from training (ρ = −0.625); FINDINGS_MD4.md found a weaker version (ρ = −0.280). Pluripotent cells sit furthest from fibroblast training. If extrapolation pushes readings down with distance, it manufactures lockstep, not a bend — so a bend observed in the frozen ruler would be conservative with respect to extrapolation, and lockstep would be uninterpretable. State whichever applies from the numbers, and say so plainly. Do not use this to rescue a lockstep result.


## Task 2 — extrapolation to GTEx fibroblast training (report-only, no gate)

PCA k=50 nearest-neighbour Euclidean and Mahalanobis, same measures as FINDINGS_FIBRO2.md / `src/fibro2_extrap.py`. TMM among that donor's scored-state rows (not the md4 state×timepoint rows). Frozen μ/σ/w. Not a gate.

| cell_line | label | n_cells | scored | on_trajectory | age_score | nn_euclidean | mahalanobis_pca | pca_k | reason |
|---|---|---|---|---|---|---|---|---|---|
| GM00731 | Fibroblast | 5988 | True | True | +4.277 | +451.439 | +23.945 | 50 | NA |
| GM00731 | PartialReprog | 5832 | True | True | +0.593 | +463.976 | +24.384 | 50 | NA |
| GM00731 | EarlyPluripotency | 4414 | True | True | -2.675 | +491.733 | +24.973 | 50 | NA |
| GM00731 | Pluripotency | 439 | True | True | -10.546 | +630.747 | +28.506 | 50 | NA |
| GM00731 | NonReprog | 5459 | True | False | +6.807 | +469.960 | +23.761 | 50 | NA |
| GM23815 | Fibroblast | 11355 | True | True | -2.208 | +448.772 | +23.320 | 50 | NA |
| GM23815 | PartialReprog | 6481 | True | True | -4.695 | +458.154 | +23.181 | 50 | NA |
| GM23815 | EarlyPluripotency | 12961 | True | True | -8.016 | +495.576 | +25.043 | 50 | NA |
| GM23815 | Pluripotency | 1249 | True | True | -15.240 | +629.554 | +28.273 | 50 | NA |
| GM23815 | NonReprog | 1812 | True | False | -0.216 | +490.574 | +24.648 | 50 | NA |

Inherited ρ (not recomputed as a gate): FINDINGS_FIBRO3.md Task 2 nn_euclidean ρ=-0.625 CI=[-0.834, -0.32] n=22. FINDINGS_MD4.md ρ_nn=-0.28.

Direction statement (from the distances in `t2_extrap.csv` and those inherited ρ values):

FINDINGS_FIBRO3.md Task 2 (same PCA k=50 NN / Mahalanobis measures as FINDINGS_FIBRO2.md) ρ(nn_euclidean, frozen ruler)=-0.625 CI=[-0.834, -0.32] n=22. FINDINGS_MD4.md ρ_nn=-0.28 on state×timepoint rows. Negative ρ: the ruler score falls as distance from GTEx fibroblast training rises. GM00731 scored states by nn_euclidean (near→far): Fibroblast nn=+451.439 mahal=+23.945; PartialReprog nn=+463.976 mahal=+24.384; NonReprog nn=+469.960 mahal=+23.761; EarlyPluripotency nn=+491.733 mahal=+24.973; Pluripotency nn=+630.747 mahal=+28.506. GM00731: Pluripotency nn=+630.747 Fibroblast nn=+451.439; Pluripotency farther than Fibroblast=True; farthest scored state=Pluripotency. GM00731 frozen ruler reading=`mixed`. Extrapolation that pushes age readings down with distance manufactures lockstep, not a bend. Not used to reclassify Task 1. GM23815 scored states by nn_euclidean (near→far): Fibroblast nn=+448.772 mahal=+23.320; PartialReprog nn=+458.154 mahal=+23.181; NonReprog nn=+490.574 mahal=+24.648; EarlyPluripotency nn=+495.576 mahal=+25.043; Pluripotency nn=+629.554 mahal=+28.273. GM23815: Pluripotency nn=+629.554 Fibroblast nn=+448.772; Pluripotency farther than Fibroblast=True; farthest scored state=Pluripotency. GM23815 frozen ruler reading=`mixed`. Extrapolation that pushes age readings down with distance manufactures lockstep, not a bend. Not used to reclassify Task 1. Task 2 is report-only. It is not a gate and is not used to explain away a lockstep result.

report_only=True used_to_discount_task1=False

## Task 3 — figure

Written **before any PLANE figure**, 2026-09-21.

One figure per donor, three panels (frozen ruler, MD, age-up−age-down), x = identity, y = age, the four trajectory states as points with 2D bootstrap CIs joined in order, NonReprog as a separate marker, the random-direction or random-set envelope shaded behind. Same identity-axis scaling within a donor. results/plane/figures/plane_<donor>.png.


- GM00731: `<repo>\results\plane\figures\plane_GM00731.png`
- GM23815: `<repo>\results\plane\figures\plane_GM23815.png`

Three panels per donor (frozen ruler, MD, age-up−age-down). x = identity, y = age. Trajectory states as points with 2D bootstrap CIs joined in order. NonReprog as a separate marker. Random-direction / random-set envelope shaded. Same identity-axis scaling within a donor. Y-axis per instrument (instruments are not combined).

## Updated ledger (quoted from FINDINGS_MD5.md; that file is not edited)

Append this file's reading. Quote the FINDINGS_MD5.md ledger; edit nothing in it.

## Task 3 — one-line ledger (quote, do not edit the source files)

Written **before any MD5 contrast p-value**, 2026-09-19. Ledger only. Existing FINDINGS_*.md are not edited.

List, with the file each lives in, which readings currently stand and which are withdrawn:
- FINDINGS_FIBRO.md Stage 2 d0→d10 (stands, all-cell)
- FINDINGS_MD3.md Task 1 (withdrawn in md4)
- FINDINGS_MD3.md Task 2 Δρ (stands)
- FINDINGS_MD4.md Task 1 (stands, as qualified by Task 1 of this file)

Quote each; edit none of them. This is the ledger the write-up will be built from, so it must be exact.


### FINDINGS_FIBRO.md Stage 2 d0→d10 (stands, all-cell)

File: `FINDINGS_FIBRO.md`. Stands as written. Nothing here promotes d7 to the primary Stage 2 result.

> age score does not decline → the ruler reads a static donor property, not a modifiable state. Step 3 is not supported by this data.

### FINDINGS_MD3.md Task 1 (withdrawn in md4)

File: `FINDINGS_MD3.md`. Withdrawn as unsupported in `FINDINGS_MD4.md` (n₀=3 PartialReprog cells at day 0). Not edited.

> Ruler does not decline in the pooled PartialReprog state while MD does, both with their own nulls. The instruments disagree on exactly the population the published claim concerns. No winner is declared.

Status quoted: Task 1: `instruments_disagree_on_claim_population`.

### FINDINGS_MD3.md Task 2 Δρ (stands)

File: `FINDINGS_MD3.md`. Unaffected by md4 or by this file.

> Δρ CI excludes zero in favour of the ruler in 2 of 3 named cohorts. MD and/or their aging signature predict donor age materially worse than the ruler on independent fibroblast cohorts.

### FINDINGS_MD4.md Task 1 (stands, as qualified by Task 1 here)

File: `FINDINGS_MD4.md`. Not edited. Qualification is this file's Task 1 reading, per donor, above.

> Task 1: GM00731 d3: `both_lower`; d7: `both_lower`; GM23815 d10: `both_lower`. Fired: GM00731 `both_lower`; GM23815 `both_lower`.

> Both instruments agree that partially reprogrammed cells read younger. Their claim is supported and our earlier negatives reflected the wrong test (a within-state trajectory with no baseline).

### FINDINGS_PLANE.md Task 1 (this file)

File: `FINDINGS_PLANE.md`. Not a replacement of any prior ledger entry. Not an edit of FINDINGS_FIBRO.md Stage 2 or of any FINDINGS_MD5.md ledger line.

- GM00731 frozen ruler: `mixed`. mixed / other: slope_A_young=True; slope_B_young=True; slope_B_flat=False; slope_B_old=False; slopes_overlap=False; bend_beats_null=False; age_falls_A=True; age_falls_B=True; p_more_negative=0.9850746268656716

- GM00731 MD: `mixed`. mixed / other: slope_A_young=True; slope_B_young=True; slope_B_flat=False; slope_B_old=False; slopes_overlap=False; bend_beats_null=True; age_falls_A=True; age_falls_B=True; p_more_negative=0.009950248756218905

- GM00731 age-up − age-down: `bend`. age falls in segment A (slope_A CI excludes 0 in the young direction) and segment B is flat or reversed and the bend statistic beats its null. Age moves before identity is lost; rejuvenation is separable from dedifferentiation on this instrument.

- GM23815 frozen ruler: `mixed`. mixed / other: slope_A_young=True; slope_B_young=True; slope_B_flat=False; slope_B_old=False; slopes_overlap=False; bend_beats_null=False; age_falls_A=True; age_falls_B=True; p_more_negative=0.8159203980099502

- GM23815 MD: `mixed`. mixed / other: slope_A_young=True; slope_B_young=True; slope_B_flat=False; slope_B_old=False; slopes_overlap=False; bend_beats_null=True; age_falls_A=True; age_falls_B=True; p_more_negative=0.004975124378109453

- GM23815 age-up − age-down: `bend`. age falls in segment A (slope_A CI excludes 0 in the young direction) and segment B is flat or reversed and the bend statistic beats its null. Age moves before identity is lost; rejuvenation is separable from dedifferentiation on this instrument.

## Limitations

1. Two donors in GSE297234.
2. States are inferred from expression, not lineage. PartialReprog / Fibroblast / EarlyPluripotency / Pluripotency / NonReprog labels are argmax of mmc3 Reprog_cell_state_signatures AddModuleScore, not lineage-tracked cells from day-0 fibroblasts.
3. A state ordering is not a measured trajectory. The order Fibroblast → PartialReprog → EarlyPluripotency → Pluripotency was fixed before scoring; it is not a Slingshot or time-resolved path of the same cell.
4. Timepoints are pooled within state. The state, not the day, is the trajectory position. Timepoint composition of every state is reported and is not a substitute for a day-resolved path.
5. The ruler was trained on GTEx V10 cultured fibroblasts only, public `AGE` 10-year bins. Per-cell scores are the md4 log2-CPM (TMM nf=1) vectors; not refit. Missing overlap genes at z=0.
6. Louvain is a Python deviation from their R/sctransform pipeline, as recorded in FINDINGS_MD2.md (LogNormalize + quadratic HVG + percent.mt residualization + PCA 50 + SNN + networkx louvain, resolution=0.8 Seurat default; STAR Methods omit resolution).
7. Our reproduction of their metrics is ours, not theirs. AddModuleScore is a Python reimplementation of the published algorithm (nbin=24 ctrl=100), not Seurat's C++/R object.
8. Cross-platform shift from GTEx bulk polyA (RNASeQCv2.4.2) to 10x 3' scRNA-seq.
9. MIN_CELLS=30. The threshold is not lowered. Donors are not pooled. Instruments are not combined.
10. Seed `20260914`. Bootstrap seed `20260918`. n_perm=200, n_boot=200, n_random=200.
11. GSE325735 not opened. d7 is not the Stage 2 endpoint. FINDINGS_FIBRO.md Stage 2 d0→d10 is not changed.
12. Task 2 distances are on TMM state-pooled pseudobulks; Task 1 ruler scores are per-cell log2-CPM. Those two ruler numbers are not interchangeable. Task 2 is not used to explain away a lockstep result.
13. Frozen pluripotency list ['POU5F1', 'NANOG', 'LIN28A', 'SALL4', 'DPPA4', 'ZFP42', 'DNMT3B']; fibroblast identity ['COL1A1', 'COL1A2', 'THY1', 'S100A4', 'FN1', 'VIM', 'POSTN']; identity version md4 pluri_primary drop_oskm=False.

## Gene-list versions (from disk, never from memory)

- MD n=205 Age up n=1533 Age down n=2007
- used_as_deseq2_rebuild=False used_as_published_lists=True
- identity_version=md4_pluri_primary_drop_oskm_False
- pluripotency endogenous: ['POU5F1', 'NANOG', 'LIN28A', 'SALL4', 'DPPA4', 'ZFP42', 'DNMT3B']
- fibroblast identity: ['COL1A1', 'COL1A2', 'THY1', 'S100A4', 'FN1', 'VIM', 'POSTN']

- GM00731 S3 AMS meta: {'age_up': {'n_requested': 1533, 'n_mapped': 896, 'n_missing': 636, 'missing': ['ADH1B', 'SLC25A52', 'FMO2', 'ENSG00000263244', 'TMSB4XP8', 'ENSG00000290383', 'ENSG00000283907', 'ENSG00000220370', 'RPL36AP21', 'PPIAP82', 'ENSG00000271047', 'EEF1A1P1', 'ENSG00000254732', 'TMSB10P1', 'ENSG00000196656', 'ENSG00000270249', 'RPL12P38', 'LINC01616', 'MT-TA', 'RPL5P8', 'RPSAP53', 'SUMO1P4', 'RPL39P5', 'ENSG00000236015', 'ENSG00000239719', 'PRELID1P6', 'EEF1A1P50', 'ENSG00000289470', 'ENSG00000236768', 'ACKR4P1', 'ENSG00000205037', 'ENSG00000287666', 'ENSG00000286289', 'ENSG00000234287', 'RPL3P4', 'ENSG00000290003', 'PSMC2P1', 'SLC25A51P4', 'RPS4XP22', 'RPS14P8', 'ENSG00000240898', 'HSPE1P3', 'ENSG00000269983', 'LDHAP7', 'PPIAP65', 'ENSG00000285914', 'ENSG00000268193', 'FTH1P7', 'ENSG00000286033', 'EEF1A1P16', 'PPIAP3', 'MRPS18CP7', 'RPL6P20', 'ENSG00000289316', 'PSMC1P1', 'SNX19P2', 'ENSG00000279392', 'ENSG00000272872', 'LINC02908', 'ENSG00000232454', 'MICOS10P2', 'OR7E28P', 'ARMC10P1', 'GAPDHP46', 'ENSG00000290018', 'RPS3AP36', 'ENSG00000286282', 'HSP90AB2P', 'TBC1D3P1-DHX40P1', 'MTND4P12', 'ENSG00000235400', 'MTIF2P1', 'FTH1P15', 'DPP4-DT', 'RPL7P11', 'ENSG00000273017', 'PARP4P2', 'ENSG00000234268', 'ENSG00000254373', 'ENSG00000289860', 'EIF2S2P4', 'ENSG00000277954', 'ENSG00000226965', 'CTAGE11P', 'CCDC54-AS1', 'ENSG00000225170', 'RPS3AP49', 'ENSG00000259623', 'ENSG00000286892', 'ENSG00000289171', 'SMC3P1', 'ENSG00000250049', 'ENSG00000272010', 'RPL36AP37', 'ENSG00000253204', 'ENSG00000280180', 'CRADD-AS1', 'EEF1A1P25', 'RPL7P6', 'ENSG00000225519', 'RPL10P3', 'ENSG00000289045', 'RPL5P29', 'ENSG00000261588', 'GAPDHP35', 'ENSG00000223711', 'ENSG00000260517', 'RPL6P13', 'ENSG00000257879', 'ENSG00000285082', 'ENSG00000279805', 'DUXAP10', 'ADAMDEC1', 'DKK\xa02.00', 'ENSG00000233597', 'ENSG00000239470', 'ENSG00000248583', 'ENSG00000272825', 'ENSG00000279208', 'CKS1BP3', 'RPS7P4', 'ENSG00000278673', 'ENSG00000287723', 'ENSG00000256204', 'DSTNP4', 'ENSG00000272864', 'MRPL30P1', 'APP-DT', 'ENSG00000260772', 'DUXAP9', 'ENSG00000227482', 'EEF1A1P10', 'ENSG00000280339', 'ENSG00000229021', 'ENSG00000267559', 'ENSG00000225790', 'ENSG00000226622', 'ENSG00000286909', 'ENSG00000290124', 'ENSG00000275426', 'ENSG00000286695', 'ENSG00000278601', 'ENSG00000277245', 'ENSG00000269902', 'IL6STP1', 'ENSG00000260197', 'ENSG00000261600', 'ENSG00000273035', 'PIGFP1', 'EIF4BP7', 'FGF7P8', 'RPSAP18', 'MT-TL1', 'ENSG00000272209', 'ENSG00000225498', 'TPT1P5', 'ENSG00000278981', 'ENSG00000291007', 'ENSG00000266968', 'ADAMTS16-DT', 'ENSG00000273012', 'ENSG00000284685', 'ENSG00000219133', 'ENSG00000259828', 'ENSG00000271730', 'HMGN2P17', 'ENSG00000279320', 'ENSG00000238150', 'GAPDHP20', 'ENSG00000236501', 'ENSG00000289968', 'PKD1L2', 'ENSG00000272156', 'ENSG00000242951', 'ENSG00000289880', 'ENSG00000254343', 'EEF1A1P8', 'UQCRFS1P1', 'ENSG00000240634', 'ENSG00000173867', 'ENSG00000289027', 'ENSG00000279108', 'ENSG00000244593', 'ENSG00000241352', 'ENSG00000280106', 'ENSG00000280778', 'ENSG00000286401', 'ENSG00000273368', 'ENSG00000267152', 'ENSG00000261671', 'ENSG00000287278', 'ENSG00000275854', 'ENSG00000272386', 'BZW1P2', 'ENSG00000272420', 'CYP2G1P', 'MACC1-DT', 'ENSG00000242858', 'ENSG00000287367', 'ENSG00000272181', 'ENSG00000276174', 'ENSG00000287008', 'ENSG00000227192', 'PPIAP58', 'ENSG00000286545', 'ENSG00000226945', 'CSNK1A1P1', 'ENSG00000289003', 'ENSG00000291188', 'EEF1GP1', 'ENSG00000236114', 'ENSG00000225279', 'NDUFAF4P3', 'ENSG00000272320', 'ROCK1P1', 'ENSG00000226526', 'ENSG00000286342', 'ENSG00000249664', 'MIR4800', 'ENSG00000259711', 'ENSG00000228509', 'EVX2', 'ENSG00000259648', 'ENSG00000240401', 'ENSG00000273733', 'ENSG00000260971', 'ENSG00000278492', 'ENSG00000272842', 'ENSG00000273077', 'HSPA8P1', 'ENSG00000271553', 'ENSG00000290597', 'ENSG00000259807', 'ENSG00000283413', 'ENSG00000259436', 'ENSG00000288899', 'ENSG00000253968', 'ENSG00000260160', 'ENSG00000280200', 'ENSG00000287109', 'ENSG00000272112', 'RPL7AP27', 'ENSG00000266573', 'ENSG00000289449', 'ENSG00000278989', 'ENSG00000271930', 'MTND1P23', 'RPL15P14', 'ENSG00000278847', 'ENSG00000272862', 'ENSG00000280241', 'ENSG00000273264', 'ENSG00000272279', 'ENSG00000278962', 'ENSG00000260400', 'ENSG00000272777', 'ENSG00000286314', 'ENSG00000226578', 'ENSG00000273149', 'ENSG00000280159', 'RPL5P5', 'ENSG00000284977', 'ENSG00000275580', 'GOT2P2', 'ENSG00000286561', 'ENSG00000272936', 'ENSG00000287090', 'ITGB8-AS1', 'NBEAP2', 'ENSG00000259865', 'CYP2B7P', 'LGR4-AS1', 'YWHAZP6', 'ENSG00000289210', 'ENSG00000287998', 'HNRNPA3P1', 'ENSG00000257496', 'ENSG00000218027', 'MTND5P11', 'ENSG00000230623', 'ENSG00000288770', 'EEF1A1P29', 'EEF1A1P7', 'ENSG00000272650', 'ENSG00000243181', 'ENSG00000277806', 'LRRK2-DT', 'ENSG00000289397', 'ENSG00000288815', 'FTLP2', 'PPIAP9', 'ENSG00000273221', 'EIF4A2P3', 'ENSG00000270096', 'ENSG00000288583', 'ENSG00000280022', 'ENSG00000271228', 'TMSB4XP4', 'ENSG00000276223', 'ENSG00000271420', 'ENSG00000286215', 'GASK1B-AS1', 'ENSG00000278058', 'ENSG00000285569', 'ENSG00000273181', 'ENSG00000281904', 'MTCO2P22', 'NT5C3AP1', 'ENSG00000290616', 'ENSG00000232934', 'ENSG00000279149', 'ENSG00000285831', 'ENSG00000274561', 'ENSG00000259171', 'ENSG00000276417', 'ENSG00000279254', 'ENSG00000276957', 'HNRNPA3P12', 'ENSG00000270482', 'ENSG00000288888', 'SETP20', 'ENSG00000267160', 'LINC02861', 'ENSG00000203327', 'KPNA2P3', 'LY75-CD302', 'ENSG00000228510', 'DSEL-AS1', 'ENSG00000272755', 'ENSG00000272316', 'ENSG00000280310', 'ENSG00000239920', 'ENSG00000286147', 'ENSG00000288025', 'ENSG00000285873', 'ENSG00000288632', 'RPL5P18', 'ENSG00000288543', 'ENSG00000289613', 'ENSG00000273402', 'FKBP9P1', 'ENSG00000254859', 'ENSG00000274281', 'IFITM3P2', 'ENSG00000288007', 'ELAPOR2', 'LINC02986', 'ENSG00000159239', 'ENSG00000279289', 'ENSG00000288049', 'RPL32P29', 'ENSG00000254810', 'EEF1A1P35', 'ENSG00000244063', 'ENSG00000178715', 'UBR5-DT', 'ENSG00000288172', 'ENSG00000289443', 'SOCAR', 'ENSG00000272733', 'ENSG00000276334', 'ENSG00000286828', 'ENSG00000274987', 'FAM156B', 'RPS23P8', 'ENSG00000278991', 'ENSG00000274515', 'CBX3P2', 'ENSG00000289689', 'ENSG00000226992', 'ENSG00000270031', 'ENSG00000271259', 'RPL31P17', 'H4C8', 'ENSG00000272885', 'ENSG00000260578', 'ENSG00000290879', 'ENSG00000278206', 'ENSG00000267546', 'ENSG00000270179', 'RPL5P24', 'ENSG00000284773', 'ENSG00000278385', 'ENSG00000268205', 'ENSG00000273284', 'FTH1P11', 'ENSG00000267711', 'ENSG00000255366', 'PRELID1P1', 'ENSG00000291042', 'ENSG00000260743', 'ENSG00000288924', 'ENSG00000289145', 'ENSG00000288890', 'ENSG00000213939', 'ENSG00000287419', 'ENSG00000279453', 'ENSG00000244538', 'EIF2S2P3', 'ENSG00000268401', 'ENSG00000272379', 'ENSG00000281100', 'ENSG00000270589', 'ENSG00000242861', 'ENSG00000261019', 'ENSG00000267838', 'RPL6P27', 'EIF3FP3', 'FTH1P4', 'ENSG00000279199', 'ENSG00000283175', 'ENSG00000253330', 'ENSG00000287119', 'ENSG00000279118', 'ENSG00000279619', 'ENSG00000273226', 'RPL18AP3', 'OR7E102P', 'ENSG00000289507', 'ZNNT1', 'ENSG00000272444', 'ENSG00000287619', 'ENSG00000272338', 'ENSG00000272732', 'ENSG00000225096', 'RPL7P32', 'THBS1-IT1', 'ENSG00000284669', 'EIF5-DT', 'EIF4BP6', 'ENSG00000230532', 'ENSG00000278959', 'ENSG00000276718', 'ENSG00000272463', 'ENSG00000279722', 'ENSG00000251023', 'RANP4', 'EEF1A1P17', 'ENSG00000276855', 'ENSG00000254539', 'PDC-AS1', 'ENSG00000272144', 'ENSG00000272983', 'ENSG00000234431', 'ENSG00000269892', 'ENSG00000235609', 'ENSG00000273384', 'ENSG00000232692', 'ATF4P4', 'ENSG00000277152', 'ENSG00000279204', 'ENSG00000261542', 'ENSG00000251615', 'ENSG00000278071', 'ENSG00000290046', 'FTH1P16', 'ENSG00000279133', 'ENSG00000277662', 'ENSG00000259940', 'RPL37P6', 'ENSG00000280099', 'MIR221', 'ENSG00000291233', 'ENSG00000203644', 'ENSG00000271882', 'EEF1A1P24', 'RPS3AP6', 'ENSG00000261625', 'ENSG00000286724', 'ENSG00000255089', 'ENSG00000288826', 'ENSG00000259994', 'HEXIM2-AS1', 'ENSG00000288018', 'ENSG00000261487', 'ENSG00000279865', 'ENSG00000279267', 'ENSG00000250041', 'AHI1-DT', 'MIR130AHG', 'ENSG00000290032', 'IRF1-AS1', 'ENSG00000286913', 'ENSG00000273014', 'COSMOC', 'ENSG00000272375', 'ENSG00000275601', 'ENSG00000273156', 'ENSG00000287414', 'ENSG00000273183', 'ENSG00000272182', 'ENSG00000290876', 'ENSG00000274444', 'EEF1A1P3', 'ENSG00000261888', 'DNAI3', 'ENSG00000272970', 'ENSG00000259976', 'ENSG00000288548', 'EEF1A1P6', 'ENSG00000273073', 'ENSG00000275367', 'RPL21P16', 'LINC02941', 'ENSG00000266910', 'ENSG00000276007', 'ENSG00000290578', 'ENSG00000262089', 'ENSG00000283352', 'ENSG00000279089', 'SYT9-AS1', 'RPS4XP16', 'ENSG00000272335', 'ENSG00000272630', 'ENSG00000286366', 'ENSG00000254477', 'ENSG00000272941', 'ENSG00000287262', 'MIRLET7IHG', 'ENSG00000275764', 'FTH1P23', 'ENSG00000255389', 'FTH1P10', 'ENSG00000272345', 'ENSG00000280157', 'ENSG00000279923', 'ENSG00000267481', 'ENSG00000289223', 'ENSG00000260257', 'EEF1A1P4', 'ENSG00000285851', 'ENSG00000260855', 'LINC02977', 'ENSG00000273271', 'ENSG00000290399', 'ENSG00000251867', 'ENSG00000289361', 'ENSG00000278730', 'ENSG00000274315', 'ENSG00000260296', 'ENSG00000288829', 'ENSG00000289007', 'ENSG00000269044', 'LINC00933', 'ENSG00000203392', 'ENSG00000288973', 'ENSG00000271976', 'ENSG00000289574', 'ENSG00000272807', 'SNX25P1', 'ENSG00000291032', 'ENSG00000273680', 'EEF1A1P22', 'ENSG00000272341', 'ENSG00000271993', 'ENSG00000272692', 'ENSG00000261468', 'ENSG00000261799', 'H2AC6', 'AMZ2P1', 'ENSG00000259972', 'ENSG00000233178', 'ENSG00000272977', 'ENSG00000272040', 'RPL17-C18ORF32', 'ENSG00000236263', 'ENSG00000275120', 'ENSG00000272750', 'ENSG00000276724', 'ENSG00000239763', 'ENSG00000279882', 'ENSG00000289463', 'ENSG00000288644', 'ENSG00000267904', 'FGGY-DT', 'ENSG00000280077', 'ENSG00000279041', 'ENSG00000271327', 'DARS1-AS1', 'ENSG00000253636', 'ENSG00000270055', 'ENSG00000243155', 'USP34-DT', 'ENSG00000270091', 'C10ORF95-AS1', 'ENSG00000272129', 'ENSG00000274760', 'H2BC21', 'DPY19L1P1', 'ENSG00000248927', 'FAM151B-DT', 'EEF1A1P12', 'ENSG00000246090', 'ENSG00000261098', 'ENSG00000232611', 'ENSG00000285804', 'ENSG00000214558', 'ENSG00000279696', 'ENSG00000289754', 'ENSG00000273576', 'ENSG00000280047', 'ENSG00000283236', 'ENSG00000258768', 'MFF-DT', 'ENSG00000283228', 'ENSG00000260948', 'ENSG00000271533', 'CUTALP', 'ENSG00000230551', 'ENSG00000262526', 'ENSG00000290537', 'ENSG00000276278', 'ENSG00000289106', 'ENSG00000286964', 'ENSG00000279253', 'ENSG00000272909', 'ENSG00000272604', 'ENSG00000261324', 'HNRNPA1L3', 'ENSG00000235859', 'ENSG00000291081', 'ENSG00000230896', 'ENSG00000272772', 'ENSG00000253106', 'ENSG00000291066', 'WASL-DT', 'CFAP418', 'RIMOC1', 'RP9P', 'ENSG00000259959', 'ENSG00000279059', 'ENSG00000278784', 'CHASERR', 'GTF2IP4', 'ENSG00000285796', 'ENSG00000282393'], 'n_ctrl': 24488, 'nbin': 24, 'ctrl_per_gene': 100, 'seed': 20260914}, 'age_down': {'n_requested': 2007, 'n_mapped': 1586, 'n_missing': 421, 'missing': ['TMEM191A', 'TMEM191A', 'RNA5-8SN3', 'RNVU1-29', 'RNU1-2', 'ENSG00000243679', 'CGB5', 'UBE2SP2', 'ENSG00000288725', 'FAM86DP', 'AGAP11', 'ENSG00000111780', 'ENSG00000235510', 'PSMC1P5', 'CGB8', 'URGCP-MRPS24', 'SLC2A3P1', 'ENSG00000260170', 'TRIM6-TRIM34', 'SLC7A5P2', 'ENSG00000256663', 'H4C11', 'ENSG00000236120', 'ENSG00000277047', 'ENSG00000217231', 'ENSG00000269981', 'ENSG00000264545', 'TCF21', 'PPAN-P2RY11', 'ENSG00000257767', 'MT2P1', 'ENSG00000228981', 'PRKAR1B-AS1', 'ENSG00000288698', 'GOLGA8F', 'ENSG00000268287', 'ENSG00000269706', 'ENSG00000288529', 'CORO7-PAM16', 'ENSG00000260899', 'ENSG00000267314', 'RN7SKP160', 'ENSG00000253163', 'SMIM45', 'ENSG00000259316', 'ENSG00000291210', 'OR7E47P', 'PDCD6-AHRR', 'USP32P1', 'ENSG00000287828', 'RANBP1P1', 'ENSG00000288637', 'PPIAP29', 'ENSG00000263620', 'ZFP91-CNTF', 'ENSG00000272163', 'ENSG00000267335', 'PKD1P6', 'ENSG00000258924', 'ENSG00000285505', 'NOC2LP2', 'ENSG00000285382', 'ENSG00000286001', 'ENSG00000288564', 'ENSG00000288685', 'NUS1P1', 'COMMD3-BMI1', 'SYT15B', 'ENSG00000286048', 'ENSG00000268108', 'ENSG00000250575', 'ENSG00000288604', 'ENSG00000244716', 'ENSG00000222032', 'ENSG00000268655', 'HMGN2P5', 'H2AC13', 'H2BP1', 'ENSG00000286518', 'PKD1P3', 'ENSG00000248710', 'ENSG00000291273', 'ENSG00000222022', 'ENSG00000285976', 'ENSG00000287235', 'H2AC20', 'ADAMTS7P1', 'ENSG00000250657', 'ENSG00000289701', 'ENSG00000290920', 'H3C3', 'H3C8', 'ENSG00000230202', 'UBE2L4', 'SPANXA2', 'ENSG00000268434', 'HNRNPCL2', 'ENSG00000257522', 'UBE2SP1', 'TEKTIP1', 'H2AX', 'ENSG00000278546', 'CFL1P4', 'COX6A1P2', 'ENSG00000291288', 'ENSG00000262413', 'EDIL3-DT', 'NHP2P1', 'MIR3652', 'ENSG00000267059', 'H3P4', 'HSPE1P4', 'ZNF286A-TBC1D26', 'SETP8', 'MBTPS1-DT', 'ENSG00000288728', 'ENSG00000248187', 'ENSG00000275092', 'PGM5P2', 'ENSG00000227080', 'ENSG00000254719', 'MSRB3-AS1', 'HOATZ', 'RNF207-AS1', 'NTF6B', 'DBF4P1', 'HNRNPA1P7', 'ARL6IP1P2', 'ENSG00000272512', 'ENSG00000283217', 'ENSG00000279140', 'ENSG00000214776', 'ENSG00000273420', 'ENSG00000268750', 'MT1XP1', 'ALDOAP1', 'H2BC9', 'ENSG00000268173', 'ENSG00000255054', 'SNX18P7', 'H2AC11', 'TMEM256-PLSCR3', 'MGC15885', 'ATP1B3P1', 'ENSG00000272702', 'ENSG00000223612', 'C10ORF88B', 'ENSG00000276131', 'DBIL5P', 'RNPS1P1', 'ENSG00000255983', 'ENSG00000289297', 'HMGN2P4', 'ENSG00000288844', 'ENSG00000274825', 'ZNF849P', 'ENSG00000287064', 'EIF1AXP1', 'ENSG00000285730', 'CLIC4P1', 'ENSAP2', 'ENSG00000288772', 'RAB43P1', 'ENSG00000269927', 'ENSG00000256802', 'ENSG00000283235', 'ENSG00000279337', 'ACTN4P1', 'ENSG00000290230', 'ENSG00000250038', 'ACTG1P10', 'ENSG00000259522', 'STK25P1', 'MT1P3', 'ENSG00000241634', 'DDX12P', 'ENSG00000272438', 'C16ORF95-DT', 'MADCAM1-AS1', 'ENSG00000279825', 'ENSG00000269210', 'RPSAP15', 'FSCN1P1', 'ENSG00000235605', 'FAM192BP', 'ABTB3', 'ENSG00000289847', 'GLTPP1', 'ENSG00000250264', 'ENSG00000258529', 'ENSG00000236680', 'ENSG00000289635', 'G3BP1P1', 'HNRNPCP1', 'MIR4280HG', 'ENSG00000214593', 'ENSG00000291105', 'GRASLND', 'ENSG00000273734', 'ENSG00000245869', 'ENSG00000273175', 'GREM1-AS1', 'COPS8P2', 'ENSG00000267698', 'SIGLEC10-AS1', 'ENSG00000255872', 'C5ORF34-AS1', 'ENSG00000289151', 'ENSG00000288891', 'ENSG00000291234', 'ATP6V1B1-AS1', 'NOC2LP1', 'PIGAP1', 'ENSG00000286192', 'ENSG00000228477', 'ENSG00000273796', 'ENSG00000277368', 'RPS7P3', 'ENSG00000271709', 'ENSG00000273361', 'ENSG00000232626', 'CSPG4P13', 'ENSG00000232389', 'ENSG00000290886', 'ENSG00000274272', 'RPS21P4', 'RASA4DP', 'ENSG00000251468', 'ENSG00000232499', 'ENSG00000213480', 'ATF6-DT', 'CROCCP3', 'ENSG00000111788', 'MYL6P1', 'ENSG00000284906', 'ENSG00000216285', 'ENSG00000268603', 'RPS4XP3', 'ENSG00000275740', 'GPX1P1', 'HMGB1P5', 'GAPDHP61', 'ENSG00000272540', 'MAGOH3P', 'EIF4A1P4', 'ENSG00000233558', 'PKD1P4', 'ENSG00000275993', 'EMSLR', 'ENSG00000277383', 'PSME2P1', 'CENPS-CORT', 'DGKZP1', 'ICMT-DT', 'ENSG00000187186', 'ENSG00000271737', 'ENSG00000235852', 'SENP3-EIF4A1', 'RANP1', 'ST13P19', 'ENO1P1', 'SNRPGP10', 'LASTR', 'ENSG00000284946', 'RTEL1-TNFRSF6B', 'ENSG00000266677', 'FMC1-LUC7L2', 'RPL12P4', 'EEF1DP1', 'ACTBP2', 'ENSG00000289612', 'ENSG00000280893', 'ENSG00000222000', 'ENSG00000272980', 'ENSG00000289318', 'TFDP1P2', 'ENSG00000267048', 'RNASEK-C17ORF49', 'ENSG00000267228', 'HMGA1P2', 'ENSG00000278817', 'H1-10', 'ENSG00000219023', 'DLSTP1', 'ENSG00000180015', 'HSPA8P4', 'ENSG00000213080', 'ZNF887P', 'ENSG00000259132', 'ENSG00000250697', 'ARHGAP11A-SCG5', 'MPRIPP1', 'ENSG00000267001', 'ENSG00000280287', 'BOP\xa01.00', 'ENSG00000234961', 'ENSG00000279568', 'ENSG00000249738', 'DHFRP1', 'SERBP1P6', 'H2AC19', 'TUBAP2', 'ENSG00000230699', 'ENSG00000279827', 'ATP6V1G2-DDX39B', 'ENSG00000289260', 'ENSG00000280571', 'ENSG00000273760', 'ENSG00000239246', 'SNORD17', 'ENSG00000290386', 'SNRPEP2', 'CCNP', 'DYNC2I2', 'POC1B-GALNT4', 'RPL5P12', 'RBMS1P1', 'RPS15P5', 'ENSG00000261253', 'TRABD-AS1', 'MIR3677HG', 'ENSG00000289006', 'ENSG00000289293', 'ENSG00000283390', 'ENSG00000290394', 'ENSG00000228818', 'HNRNPCP2', 'ENSG00000285446', 'ENSG00000290948', 'HROB', 'ENSG00000288856', 'ENSG00000274248', 'ENSG00000288550', 'ENSG00000289517', 'ENSG00000289486', 'GAPDHP60', 'ENSG00000197332', 'SLX9', 'ENSG00000287356', 'NUP153-AS1', 'RPL41P2', 'ENSG00000271936', 'NME2P1', 'ENSG00000284874', 'ENSG00000272277', 'ENSG00000289653', 'ITGB1P1', 'ENSG00000167774', 'DPY19L2P2', 'BBLN', 'ENSG00000258691', 'ENSG00000273154', 'ENSG00000289123', 'RPL10AP6', 'ENSG00000244313', 'ENSG00000213985', 'H2AZ1', 'SETP14', 'MSNP1', 'ENSG00000234241', 'ENSG00000260342', 'ENSG00000290045', 'ENSG00000285641', 'SRP9P1', 'ENSG00000291224', 'ENSG00000276672', 'ENSG00000266709', 'ENSG00000288663', 'SBDSP1', 'BRME1', 'ENSG00000254910', 'CA5BP1', 'H2AW', 'PTGES3P1', 'ENSG00000283761', 'NAV2-AS6', 'ENSG00000260877', 'POLR1G', 'LSM\xa04.00', 'ENSG00000285043', 'FAUP1', 'ENSG00000241889', 'TOLLIP-DT', 'TUBA3FP', 'RPS2P46', 'ENSG00000279415', 'ZNRF2P1', 'HNRNPLP2', 'ARHGEF17-AS1', 'ENSG00000289579', 'YBX1P1', 'SSBP3P1', 'CALM2P2', 'CSPG4P10', 'LSM\xa07.00', 'ENSG00000287168', 'ENSG00000284292', 'ST20-MTHFS', 'TUBBP1', 'PPIAP22', 'CFAP119', 'ENSG00000280239', 'GAPDHP1', 'ENSG00000273398', 'MARCHF3', 'ENSG00000288983', 'ENSG00000289720', 'ENSG00000290689', 'ENSG00000286112', 'ENSG00000261130', 'ENSG00000218426', 'ENSG00000274425', 'ENSG00000255639', 'H2AJ', 'ZFHX3-AS1', 'ENSG00000259529', 'MORF4L1P1', 'RPLP0P6', 'ENSG00000249209', 'RPL13P5', 'ENSG00000184441'], 'n_ctrl': 24313, 'nbin': 24, 'ctrl_per_gene': 100, 'seed': 20260914}}

- GM23815 S3 AMS meta: {'age_up': {'n_requested': 1533, 'n_mapped': 894, 'n_missing': 638, 'missing': ['ADH1B', 'KCNA4', 'SLC25A52', 'ENSG00000263244', 'TMSB4XP8', 'ENSG00000290383', 'ENSG00000283907', 'ENSG00000220370', 'RPL36AP21', 'PPIAP82', 'ENSG00000271047', 'EEF1A1P1', 'ENSG00000254732', 'TMSB10P1', 'ENSG00000196656', 'ENSG00000270249', 'RPL12P38', 'LINC01616', 'MT-TA', 'RPL5P8', 'RPSAP53', 'SUMO1P4', 'RPL39P5', 'ENSG00000236015', 'ENSG00000239719', 'PRELID1P6', 'EEF1A1P50', 'ENSG00000289470', 'ENSG00000236768', 'ACKR4P1', 'ENSG00000205037', 'ENSG00000287666', 'ENSG00000286289', 'ENSG00000234287', 'RPL3P4', 'ENSG00000290003', 'PSMC2P1', 'SLC25A51P4', 'RPS4XP22', 'RPS14P8', 'ENSG00000240898', 'HSPE1P3', 'ENSG00000269983', 'LDHAP7', 'PPIAP65', 'ENSG00000285914', 'ENSG00000268193', 'FTH1P7', 'ENSG00000286033', 'EEF1A1P16', 'PPIAP3', 'MRPS18CP7', 'RPL6P20', 'ENSG00000289316', 'PSMC1P1', 'SNX19P2', 'ENSG00000279392', 'ENSG00000272872', 'LINC02908', 'ENSG00000232454', 'MICOS10P2', 'OR7E28P', 'ARMC10P1', 'GAPDHP46', 'ENSG00000290018', 'RPS3AP36', 'ENSG00000286282', 'HSP90AB2P', 'TBC1D3P1-DHX40P1', 'MTND4P12', 'ENSG00000235400', 'MTIF2P1', 'FTH1P15', 'DPP4-DT', 'RPL7P11', 'ENSG00000273017', 'PARP4P2', 'ENSG00000234268', 'ENSG00000254373', 'ENSG00000289860', 'EIF2S2P4', 'ENSG00000277954', 'ENSG00000226965', 'CTAGE11P', 'CCDC54-AS1', 'ENSG00000225170', 'RPS3AP49', 'ENSG00000259623', 'ENSG00000286892', 'ENSG00000289171', 'SMC3P1', 'ENSG00000250049', 'ENSG00000272010', 'RPL36AP37', 'ENSG00000253204', 'ENSG00000280180', 'CRADD-AS1', 'EEF1A1P25', 'RPL7P6', 'ENSG00000225519', 'RPL10P3', 'ENSG00000289045', 'RPL5P29', 'ENSG00000261588', 'GAPDHP35', 'ENSG00000223711', 'ENSG00000260517', 'MYH1', 'RPL6P13', 'ENSG00000257879', 'ENSG00000285082', 'ENSG00000279805', 'DUXAP10', 'ADAMDEC1', 'DKK\xa02.00', 'ENSG00000233597', 'ENSG00000239470', 'ENSG00000248583', 'ENSG00000272825', 'ENSG00000279208', 'CKS1BP3', 'RPS7P4', 'ENSG00000278673', 'ENSG00000287723', 'ENSG00000256204', 'DSTNP4', 'ENSG00000272864', 'MRPL30P1', 'APP-DT', 'ENSG00000260772', 'DUXAP9', 'ENSG00000227482', 'EEF1A1P10', 'ENSG00000280339', 'ENSG00000229021', 'ENSG00000267559', 'ENSG00000225790', 'ENSG00000226622', 'ENSG00000286909', 'ENSG00000290124', 'ENSG00000275426', 'ENSG00000286695', 'ENSG00000278601', 'ENSG00000277245', 'ENSG00000269902', 'IL6STP1', 'ENSG00000260197', 'ENSG00000261600', 'ENSG00000273035', 'PIGFP1', 'EIF4BP7', 'FGF7P8', 'RPSAP18', 'MT-TL1', 'ENSG00000272209', 'ENSG00000225498', 'TPT1P5', 'ENSG00000278981', 'ENSG00000291007', 'ENSG00000266968', 'ADAMTS16-DT', 'ENSG00000273012', 'ENSG00000284685', 'ENSG00000219133', 'ENSG00000259828', 'ENSG00000271730', 'HMGN2P17', 'ENSG00000279320', 'ENSG00000238150', 'GAPDHP20', 'ENSG00000236501', 'ENSG00000289968', 'PKD1L2', 'ENSG00000272156', 'ENSG00000242951', 'ENSG00000289880', 'ENSG00000254343', 'EEF1A1P8', 'UQCRFS1P1', 'ENSG00000240634', 'ENSG00000173867', 'ENSG00000289027', 'ENSG00000279108', 'ENSG00000244593', 'ENSG00000241352', 'ENSG00000280106', 'ENSG00000280778', 'ENSG00000286401', 'ENSG00000273368', 'ENSG00000267152', 'ENSG00000261671', 'ENSG00000287278', 'ENSG00000275854', 'ENSG00000272386', 'BZW1P2', 'ENSG00000272420', 'CYP2G1P', 'MACC1-DT', 'ENSG00000242858', 'ENSG00000287367', 'ENSG00000272181', 'ENSG00000276174', 'ENSG00000287008', 'ENSG00000227192', 'PPIAP58', 'ENSG00000286545', 'ENSG00000226945', 'CSNK1A1P1', 'ENSG00000289003', 'ENSG00000291188', 'EEF1GP1', 'ENSG00000236114', 'ENSG00000225279', 'NDUFAF4P3', 'ENSG00000272320', 'ROCK1P1', 'ENSG00000226526', 'ENSG00000286342', 'ENSG00000249664', 'MIR4800', 'ENSG00000259711', 'ENSG00000228509', 'EVX2', 'ENSG00000259648', 'ENSG00000240401', 'ENSG00000273733', 'ENSG00000260971', 'ENSG00000278492', 'ENSG00000272842', 'ENSG00000273077', 'HSPA8P1', 'ENSG00000271553', 'ENSG00000290597', 'ENSG00000259807', 'ENSG00000283413', 'ENSG00000259436', 'ENSG00000288899', 'ENSG00000253968', 'ENSG00000260160', 'ENSG00000280200', 'ENSG00000287109', 'ENSG00000272112', 'RPL7AP27', 'ENSG00000266573', 'ENSG00000289449', 'ENSG00000278989', 'ENSG00000271930', 'MTND1P23', 'RPL15P14', 'ENSG00000278847', 'ENSG00000272862', 'ENSG00000280241', 'ENSG00000273264', 'ENSG00000272279', 'ENSG00000278962', 'ENSG00000260400', 'ENSG00000272777', 'ENSG00000286314', 'ENSG00000226578', 'ENSG00000273149', 'ENSG00000280159', 'RPL5P5', 'ENSG00000284977', 'ENSG00000275580', 'GOT2P2', 'ENSG00000286561', 'ENSG00000272936', 'ENSG00000287090', 'ITGB8-AS1', 'NBEAP2', 'ENSG00000259865', 'CYP2B7P', 'LGR4-AS1', 'YWHAZP6', 'ENSG00000289210', 'ENSG00000287998', 'HNRNPA3P1', 'ENSG00000257496', 'ENSG00000218027', 'MTND5P11', 'ENSG00000230623', 'ENSG00000288770', 'EEF1A1P29', 'EEF1A1P7', 'ENSG00000272650', 'ENSG00000243181', 'ENSG00000277806', 'LRRK2-DT', 'ENSG00000289397', 'ENSG00000288815', 'FTLP2', 'PPIAP9', 'ENSG00000273221', 'EIF4A2P3', 'ENSG00000270096', 'ENSG00000288583', 'ENSG00000280022', 'ENSG00000271228', 'TMSB4XP4', 'ENSG00000276223', 'ENSG00000271420', 'ENSG00000286215', 'GASK1B-AS1', 'ENSG00000278058', 'ENSG00000285569', 'ENSG00000273181', 'ENSG00000281904', 'MTCO2P22', 'NT5C3AP1', 'ENSG00000290616', 'ENSG00000232934', 'ENSG00000279149', 'ENSG00000285831', 'ENSG00000274561', 'ENSG00000259171', 'ENSG00000276417', 'ENSG00000279254', 'ENSG00000276957', 'HNRNPA3P12', 'ENSG00000270482', 'ENSG00000288888', 'SETP20', 'ENSG00000267160', 'LINC02861', 'ENSG00000203327', 'KPNA2P3', 'LY75-CD302', 'ENSG00000228510', 'DSEL-AS1', 'ENSG00000272755', 'ENSG00000272316', 'ENSG00000280310', 'ENSG00000239920', 'ENSG00000286147', 'ENSG00000288025', 'ENSG00000285873', 'ENSG00000288632', 'RPL5P18', 'ENSG00000288543', 'ENSG00000289613', 'ENSG00000273402', 'FKBP9P1', 'ENSG00000254859', 'ENSG00000274281', 'IFITM3P2', 'ENSG00000288007', 'ELAPOR2', 'LINC02986', 'ENSG00000159239', 'ENSG00000279289', 'ENSG00000288049', 'RPL32P29', 'ENSG00000254810', 'EEF1A1P35', 'ENSG00000244063', 'ENSG00000178715', 'UBR5-DT', 'ENSG00000288172', 'ENSG00000289443', 'SOCAR', 'ENSG00000272733', 'ENSG00000276334', 'ENSG00000286828', 'ENSG00000274987', 'FAM156B', 'RPS23P8', 'ENSG00000278991', 'ENSG00000274515', 'CBX3P2', 'ENSG00000289689', 'ENSG00000226992', 'ENSG00000270031', 'ENSG00000271259', 'RPL31P17', 'H4C8', 'ENSG00000272885', 'ENSG00000260578', 'ENSG00000290879', 'ENSG00000278206', 'ENSG00000267546', 'ENSG00000270179', 'RPL5P24', 'ENSG00000284773', 'ENSG00000278385', 'OR2A1', 'ENSG00000268205', 'ENSG00000273284', 'FTH1P11', 'ENSG00000267711', 'ENSG00000255366', 'PRELID1P1', 'ENSG00000291042', 'ENSG00000260743', 'ENSG00000288924', 'ENSG00000289145', 'ENSG00000288890', 'ENSG00000213939', 'ENSG00000287419', 'ENSG00000279453', 'ENSG00000244538', 'EIF2S2P3', 'ENSG00000268401', 'ENSG00000272379', 'ENSG00000281100', 'ENSG00000270589', 'ENSG00000242861', 'ENSG00000261019', 'ENSG00000267838', 'RPL6P27', 'EIF3FP3', 'FTH1P4', 'ENSG00000279199', 'ENSG00000283175', 'ENSG00000253330', 'ENSG00000287119', 'ENSG00000279118', 'ENSG00000279619', 'ENSG00000273226', 'RPL18AP3', 'OR7E102P', 'ENSG00000289507', 'ZNNT1', 'ENSG00000272444', 'ENSG00000287619', 'ENSG00000272338', 'ENSG00000272732', 'ENSG00000225096', 'RPL7P32', 'THBS1-IT1', 'ENSG00000284669', 'EIF5-DT', 'EIF4BP6', 'ENSG00000230532', 'ENSG00000278959', 'ENSG00000276718', 'ENSG00000272463', 'ENSG00000279722', 'ENSG00000251023', 'RANP4', 'EEF1A1P17', 'ENSG00000276855', 'ENSG00000254539', 'PDC-AS1', 'ENSG00000272144', 'ENSG00000272983', 'ENSG00000234431', 'ENSG00000269892', 'ENSG00000235609', 'ENSG00000273384', 'ENSG00000232692', 'ATF4P4', 'ENSG00000277152', 'ENSG00000279204', 'ENSG00000261542', 'ENSG00000251615', 'ENSG00000278071', 'ENSG00000290046', 'FTH1P16', 'ENSG00000279133', 'ENSG00000277662', 'ENSG00000259940', 'RPL37P6', 'ENSG00000280099', 'MIR221', 'ENSG00000291233', 'ENSG00000203644', 'ENSG00000271882', 'EEF1A1P24', 'RPS3AP6', 'ENSG00000261625', 'ENSG00000286724', 'ENSG00000255089', 'ENSG00000288826', 'ENSG00000259994', 'HEXIM2-AS1', 'ENSG00000288018', 'ENSG00000261487', 'ENSG00000279865', 'ENSG00000279267', 'ENSG00000250041', 'AHI1-DT', 'MIR130AHG', 'ENSG00000290032', 'IRF1-AS1', 'ENSG00000286913', 'ENSG00000273014', 'COSMOC', 'ENSG00000272375', 'ENSG00000275601', 'ENSG00000273156', 'ENSG00000287414', 'ENSG00000273183', 'ENSG00000272182', 'ENSG00000290876', 'ENSG00000274444', 'EEF1A1P3', 'ENSG00000261888', 'DNAI3', 'ENSG00000272970', 'ENSG00000259976', 'ENSG00000288548', 'EEF1A1P6', 'ENSG00000273073', 'ENSG00000275367', 'RPL21P16', 'LINC02941', 'ENSG00000266910', 'ENSG00000276007', 'ENSG00000290578', 'ENSG00000262089', 'ENSG00000283352', 'ENSG00000279089', 'SYT9-AS1', 'RPS4XP16', 'ENSG00000272335', 'ENSG00000272630', 'ENSG00000286366', 'ENSG00000254477', 'ENSG00000272941', 'ENSG00000287262', 'MIRLET7IHG', 'ENSG00000275764', 'FTH1P23', 'ENSG00000255389', 'FTH1P10', 'ENSG00000272345', 'ENSG00000280157', 'ENSG00000279923', 'ENSG00000267481', 'ENSG00000289223', 'ENSG00000260257', 'EEF1A1P4', 'ENSG00000285851', 'ENSG00000260855', 'LINC02977', 'ENSG00000273271', 'ENSG00000290399', 'ENSG00000251867', 'ENSG00000289361', 'ENSG00000278730', 'ENSG00000274315', 'ENSG00000260296', 'ENSG00000288829', 'ENSG00000289007', 'ENSG00000269044', 'LINC00933', 'ENSG00000203392', 'ENSG00000288973', 'ENSG00000271976', 'ENSG00000289574', 'ENSG00000272807', 'SNX25P1', 'ENSG00000291032', 'ENSG00000273680', 'EEF1A1P22', 'ENSG00000272341', 'ENSG00000271993', 'ENSG00000272692', 'ENSG00000261468', 'ENSG00000261799', 'H2AC6', 'AMZ2P1', 'ENSG00000259972', 'ENSG00000233178', 'ENSG00000272977', 'ENSG00000272040', 'RPL17-C18ORF32', 'ENSG00000236263', 'ENSG00000275120', 'ENSG00000272750', 'ENSG00000276724', 'ENSG00000239763', 'ENSG00000279882', 'ENSG00000289463', 'ENSG00000288644', 'ENSG00000267904', 'FGGY-DT', 'ENSG00000280077', 'ENSG00000279041', 'ENSG00000271327', 'DARS1-AS1', 'ENSG00000253636', 'ENSG00000270055', 'ENSG00000243155', 'USP34-DT', 'ENSG00000270091', 'C10ORF95-AS1', 'ENSG00000272129', 'ENSG00000274760', 'H2BC21', 'DPY19L1P1', 'ENSG00000248927', 'FAM151B-DT', 'EEF1A1P12', 'ENSG00000246090', 'ENSG00000261098', 'ENSG00000232611', 'ENSG00000285804', 'ENSG00000214558', 'ENSG00000279696', 'ENSG00000289754', 'ENSG00000273576', 'ENSG00000280047', 'ENSG00000283236', 'ENSG00000258768', 'MFF-DT', 'ENSG00000283228', 'ENSG00000260948', 'ENSG00000271533', 'CUTALP', 'ENSG00000230551', 'ENSG00000262526', 'ENSG00000290537', 'ENSG00000276278', 'ENSG00000289106', 'ENSG00000286964', 'ENSG00000279253', 'ENSG00000272909', 'ENSG00000272604', 'ENSG00000261324', 'HNRNPA1L3', 'ENSG00000235859', 'ENSG00000291081', 'ENSG00000230896', 'ENSG00000272772', 'ENSG00000253106', 'ENSG00000291066', 'WASL-DT', 'CFAP418', 'RIMOC1', 'RP9P', 'ENSG00000259959', 'ENSG00000279059', 'ENSG00000278784', 'CHASERR', 'GTF2IP4', 'ENSG00000285796', 'ENSG00000282393'], 'n_ctrl': 24930, 'nbin': 24, 'ctrl_per_gene': 100, 'seed': 20260914}, 'age_down': {'n_requested': 2007, 'n_mapped': 1582, 'n_missing': 425, 'missing': ['TMEM191A', 'TMEM191A', 'RNA5-8SN3', 'RNVU1-29', 'RNU1-2', 'ENSG00000243679', 'CGB5', 'UBE2SP2', 'ENSG00000288725', 'FAM86DP', 'AGAP11', 'ENSG00000111780', 'ENSG00000235510', 'PSMC1P5', 'URGCP-MRPS24', 'SLC2A3P1', 'ENSG00000260170', 'TRIM6-TRIM34', 'SLC7A5P2', 'ENSG00000256663', 'H4C11', 'ENSG00000236120', 'ENSG00000277047', 'ENSG00000217231', 'ENSG00000269981', 'ENSG00000264545', 'PPAN-P2RY11', 'ENSG00000257767', 'MT2P1', 'ENSG00000228981', 'PRKAR1B-AS1', 'CBLN2', 'ENSG00000288698', 'GOLGA8F', 'ENSG00000268287', 'ENSG00000269706', 'CGB2', 'ENSG00000288529', 'CORO7-PAM16', 'ENSG00000260899', 'ENSG00000267314', 'RN7SKP160', 'ENSG00000253163', 'SMIM45', 'ENSG00000259316', 'ENSG00000291210', 'OR7E47P', 'PDCD6-AHRR', 'USP32P1', 'ENSG00000287828', 'RANBP1P1', 'ENSG00000288637', 'PPIAP29', 'ENSG00000263620', 'ZFP91-CNTF', 'ENSG00000272163', 'ENSG00000267335', 'PKD1P6', 'ENSG00000258924', 'ENSG00000285505', 'LCE3D', 'NOC2LP2', 'ENSG00000285382', 'ENSG00000286001', 'ENSG00000288564', 'ENSG00000288685', 'NUS1P1', 'COMMD3-BMI1', 'SYT15B', 'ENSG00000286048', 'ENSG00000268108', 'ENSG00000250575', 'ENSG00000288604', 'ENSG00000244716', 'ENSG00000222032', 'ENSG00000268655', 'HMGN2P5', 'H2AC13', 'H2BP1', 'ENSG00000286518', 'PKD1P3', 'ENSG00000248710', 'ENSG00000291273', 'ENSG00000222022', 'ENSG00000285976', 'ENSG00000287235', 'H2AC20', 'ADAMTS7P1', 'ENSG00000250657', 'ENSG00000289701', 'ENSG00000290920', 'H3C3', 'H3C8', 'ENSG00000230202', 'UBE2L4', 'SPANXA2', 'ENSG00000268434', 'HNRNPCL2', 'ENSG00000257522', 'UBE2SP1', 'TEKTIP1', 'H2AX', 'ENSG00000278546', 'CFL1P4', 'COX6A1P2', 'ENSG00000291288', 'ENSG00000262413', 'EDIL3-DT', 'NHP2P1', 'MIR3652', 'ENSG00000267059', 'H3P4', 'LINC00161', 'HSPE1P4', 'ZNF286A-TBC1D26', 'SETP8', 'MBTPS1-DT', 'LCE2A', 'ENSG00000288728', 'ENSG00000248187', 'ENSG00000275092', 'PGM5P2', 'ENSG00000227080', 'ENSG00000254719', 'MSRB3-AS1', 'HOATZ', 'RNF207-AS1', 'NTF6B', 'DBF4P1', 'HNRNPA1P7', 'ARL6IP1P2', 'ENSG00000272512', 'ENSG00000283217', 'ENSG00000279140', 'ENSG00000214776', 'ENSG00000273420', 'ENSG00000268750', 'MT1XP1', 'ALDOAP1', 'SST', 'H2BC9', 'ENSG00000268173', 'ENSG00000255054', 'SNX18P7', 'H2AC11', 'TMEM256-PLSCR3', 'MGC15885', 'ATP1B3P1', 'ENSG00000272702', 'ENSG00000223612', 'C10ORF88B', 'ENSG00000276131', 'DBIL5P', 'RNPS1P1', 'ENSG00000255983', 'ENSG00000289297', 'HMGN2P4', 'ENSG00000288844', 'ENSG00000274825', 'ZNF849P', 'ENSG00000287064', 'EIF1AXP1', 'ENSG00000285730', 'CLIC4P1', 'ENSAP2', 'ENSG00000288772', 'RAB43P1', 'ENSG00000269927', 'ENSG00000256802', 'ENSG00000283235', 'ENSG00000279337', 'ACTN4P1', 'ENSG00000290230', 'ENSG00000250038', 'ACTG1P10', 'ENSG00000259522', 'STK25P1', 'MT1P3', 'ENSG00000241634', 'DDX12P', 'ENSG00000272438', 'C16ORF95-DT', 'MADCAM1-AS1', 'ENSG00000279825', 'ENSG00000269210', 'RPSAP15', 'FSCN1P1', 'ENSG00000235605', 'FAM192BP', 'ABTB3', 'ENSG00000289847', 'GLTPP1', 'ENSG00000250264', 'ENSG00000258529', 'ENSG00000236680', 'ENSG00000289635', 'G3BP1P1', 'HNRNPCP1', 'MIR4280HG', 'ENSG00000214593', 'ENSG00000291105', 'GRASLND', 'ENSG00000273734', 'ENSG00000245869', 'ENSG00000273175', 'GREM1-AS1', 'COPS8P2', 'ENSG00000267698', 'SIGLEC10-AS1', 'ENSG00000255872', 'C5ORF34-AS1', 'ENSG00000289151', 'ENSG00000288891', 'ENSG00000291234', 'ATP6V1B1-AS1', 'NOC2LP1', 'PIGAP1', 'ENSG00000286192', 'ENSG00000228477', 'ENSG00000273796', 'ENSG00000277368', 'RPS7P3', 'ENSG00000271709', 'ENSG00000273361', 'ENSG00000232626', 'CSPG4P13', 'ENSG00000232389', 'ENSG00000290886', 'ENSG00000274272', 'RPS21P4', 'RASA4DP', 'ENSG00000251468', 'ENSG00000232499', 'ENSG00000213480', 'ATF6-DT', 'CROCCP3', 'ENSG00000111788', 'MYL6P1', 'ENSG00000284906', 'ENSG00000216285', 'ENSG00000268603', 'RPS4XP3', 'ENSG00000275740', 'GPX1P1', 'HMGB1P5', 'GAPDHP61', 'ENSG00000272540', 'MAGOH3P', 'EIF4A1P4', 'ENSG00000233558', 'PKD1P4', 'ENSG00000275993', 'EMSLR', 'ENSG00000277383', 'PSME2P1', 'CENPS-CORT', 'DGKZP1', 'ICMT-DT', 'ENSG00000187186', 'ENSG00000271737', 'ENSG00000235852', 'SENP3-EIF4A1', 'RANP1', 'ST13P19', 'ENO1P1', 'SNRPGP10', 'LASTR', 'ENSG00000284946', 'RTEL1-TNFRSF6B', 'ENSG00000266677', 'FMC1-LUC7L2', 'RPL12P4', 'EEF1DP1', 'ACTBP2', 'ENSG00000289612', 'ENSG00000280893', 'ENSG00000222000', 'ENSG00000272980', 'ENSG00000289318', 'TFDP1P2', 'ENSG00000267048', 'RNASEK-C17ORF49', 'ENSG00000267228', 'HMGA1P2', 'ENSG00000278817', 'H1-10', 'ENSG00000219023', 'DLSTP1', 'ENSG00000180015', 'HSPA8P4', 'ENSG00000213080', 'ZNF887P', 'ENSG00000259132', 'ENSG00000250697', 'ARHGAP11A-SCG5', 'MPRIPP1', 'ENSG00000267001', 'ENSG00000280287', 'BOP\xa01.00', 'ENSG00000234961', 'ENSG00000279568', 'ENSG00000249738', 'DHFRP1', 'SERBP1P6', 'H2AC19', 'TUBAP2', 'ENSG00000230699', 'ENSG00000279827', 'ATP6V1G2-DDX39B', 'ENSG00000289260', 'ENSG00000280571', 'ENSG00000273760', 'ENSG00000239246', 'SNORD17', 'ENSG00000290386', 'SNRPEP2', 'CCNP', 'DYNC2I2', 'POC1B-GALNT4', 'RPL5P12', 'RBMS1P1', 'RPS15P5', 'ENSG00000261253', 'TRABD-AS1', 'MIR3677HG', 'ENSG00000289006', 'ENSG00000289293', 'ENSG00000283390', 'ENSG00000290394', 'ENSG00000228818', 'HNRNPCP2', 'ENSG00000285446', 'ENSG00000290948', 'HROB', 'ENSG00000288856', 'ENSG00000274248', 'ENSG00000288550', 'ENSG00000289517', 'ENSG00000289486', 'GAPDHP60', 'ENSG00000197332', 'SLX9', 'ENSG00000287356', 'NUP153-AS1', 'RPL41P2', 'ENSG00000271936', 'NME2P1', 'ENSG00000284874', 'ENSG00000272277', 'ENSG00000289653', 'ITGB1P1', 'ENSG00000167774', 'DPY19L2P2', 'BBLN', 'ENSG00000258691', 'ENSG00000273154', 'ENSG00000289123', 'RPL10AP6', 'ENSG00000244313', 'ENSG00000213985', 'H2AZ1', 'SETP14', 'MSNP1', 'ENSG00000234241', 'ENSG00000260342', 'ENSG00000290045', 'ENSG00000285641', 'SRP9P1', 'ENSG00000291224', 'ENSG00000276672', 'ENSG00000266709', 'ENSG00000288663', 'SBDSP1', 'BRME1', 'ENSG00000254910', 'CA5BP1', 'H2AW', 'PTGES3P1', 'ENSG00000283761', 'NAV2-AS6', 'ENSG00000260877', 'POLR1G', 'LSM\xa04.00', 'ENSG00000285043', 'FAUP1', 'ENSG00000241889', 'TOLLIP-DT', 'TUBA3FP', 'RPS2P46', 'ENSG00000279415', 'ZNRF2P1', 'HNRNPLP2', 'ARHGEF17-AS1', 'ENSG00000289579', 'YBX1P1', 'SSBP3P1', 'CALM2P2', 'CSPG4P10', 'LSM\xa07.00', 'ENSG00000287168', 'ENSG00000284292', 'ST20-MTHFS', 'TUBBP1', 'PPIAP22', 'CFAP119', 'ENSG00000280239', 'GAPDHP1', 'ENSG00000273398', 'MARCHF3', 'ENSG00000288983', 'ENSG00000289720', 'ENSG00000290689', 'ENSG00000286112', 'ENSG00000261130', 'ENSG00000218426', 'ENSG00000274425', 'ENSG00000255639', 'H2AJ', 'ZFHX3-AS1', 'ENSG00000259529', 'MORF4L1P1', 'RPLP0P6', 'ENSG00000249209', 'RPL13P5', 'ENSG00000184441'], 'n_ctrl': 24081, 'nbin': 24, 'ctrl_per_gene': 100, 'seed': 20260914}}

## ID types actually read

- **louvain_GM00731_symbols** kind=`symbol` n=29058 n_ensembl=0 n_refseq=0 n_symbol_like=29057 source=`<repo>\data\processed\md2\louvain_counts_GM00731.npz:symbols` examples=['AL627309.1', 'AL627309.5', 'AP006222.2', 'AL669831.2', 'LINC01409', 'FAM87B', 'LINC01128', 'LINC00115', 'FAM41C', 'AL645608.6', 'AL645608.2', 'AL645608.4']
- **louvain_GM00731_gene_id** kind=`ensembl` n=29058 n_ensembl=29058 n_refseq=0 n_symbol_like=0 source=`<repo>\data\processed\md2\louvain_counts_GM00731.npz:gene_id` examples=['ENSG00000238009', 'ENSG00000241860', 'ENSG00000286448', 'ENSG00000229905', 'ENSG00000237491', 'ENSG00000177757', 'ENSG00000228794', 'ENSG00000225880', 'ENSG00000230368', 'ENSG00000272438', 'ENSG00000230699', 'ENSG00000241180']
- **louvain_GM23815_symbols** kind=`symbol` n=29154 n_ensembl=0 n_refseq=0 n_symbol_like=29154 source=`<repo>\data\processed\md2\louvain_counts_GM23815.npz:symbols` examples=['AL627309.1', 'AL627309.3', 'AL627309.5', 'AP006222.2', 'LINC01409', 'FAM87B', 'LINC01128', 'LINC00115', 'FAM41C', 'AL645608.6', 'AL645608.2', 'AL645608.4']
- **louvain_GM23815_gene_id** kind=`ensembl` n=29154 n_ensembl=29154 n_refseq=0 n_symbol_like=0 source=`<repo>\data\processed\md2\louvain_counts_GM23815.npz:gene_id` examples=['ENSG00000238009', 'ENSG00000239945', 'ENSG00000241860', 'ENSG00000286448', 'ENSG00000237491', 'ENSG00000177757', 'ENSG00000228794', 'ENSG00000225880', 'ENSG00000230368', 'ENSG00000272438', 'ENSG00000230699', 'ENSG00000241180']

## Columns actually read

- **md2_genesets_json** source=`<repo>\results\md2\genesets.json` columns=['MD', 'TGFB', 'EMT', 'reprog', 'mmc3_MD', 'mmc3_TGFB', 'mmc3_age_up', 'mmc3_age_down', 'msigdb_version']
- **md3_genesets_json** source=`<repo>\results\md3\genesets.json` columns=['MD', 'TGFB', 'age_up', 'age_down', 'reprog', 'mmc3_path', 'msigdb_version']
- **md4_genesets_json** source=`<repo>\results\md4\genesets.json` columns=['MD', 'TGFB', 'age_up', 'age_down']
- **t1_s3_ams_meta_GM00731** source=`<repo>\results\md4\t1_s3_ams_meta_GM00731.json` columns=['age_up', 'age_down']
- **t1_s3_ams_meta_GM23815** source=`<repo>\results\md4\t1_s3_ams_meta_GM23815.json` columns=['age_up', 'age_down']
- **t2_cluster_labels_GM00731** source=`<repo>\results\md2\t2_cluster_labels_GM00731.csv` columns=['cell_line', 'cluster', 'n_cells', 'label', 'mean_Fibroblast', 'mean_PartialReprog', 'mean_EarlyPluripotency', 'mean_Pluripotency', 'mean_NonReprog']
- **louvain_obs_GM00731** source=`<repo>\results\md2\louvain_obs_GM00731.csv` columns=['barcode', 'gsm', 'file', 'cell_line', 'day', 'age_years', 'age_source', 'umi', 'n_genes', 'mito_frac', 'cluster', 'louvain_kept', 'percent_mt']
- **louvain_ams_GM00731** source=`<repo>\data\processed\md2\louvain_ams_GM00731.npz` columns=['md', 'tgfb', 'bins', 'gene_mean', 'cluster', 'day', 'n_genes']
- **t2_cluster_labels_GM23815** source=`<repo>\results\md2\t2_cluster_labels_GM23815.csv` columns=['cell_line', 'cluster', 'n_cells', 'label', 'mean_Fibroblast', 'mean_PartialReprog', 'mean_EarlyPluripotency', 'mean_Pluripotency', 'mean_NonReprog']
- **louvain_obs_GM23815** source=`<repo>\results\md2\louvain_obs_GM23815.csv` columns=['barcode', 'gsm', 'file', 'cell_line', 'day', 'age_years', 'age_source', 'umi', 'n_genes', 'mito_frac', 'cluster', 'louvain_kept', 'percent_mt']
- **louvain_ams_GM23815** source=`<repo>\data\processed\md2\louvain_ams_GM23815.npz` columns=['md', 'tgfb', 'bins', 'gene_mean', 'cluster', 'day']
- **md4_t1_cell_counts** source=`<repo>\results\md4\t1_cell_counts.csv` columns=['cell_line', 'day', 'label', 'n_cells', 'n_timepoint', 'frac_timepoint', 'ge50']
- **t2_reproduction** source=`<repo>\results\md2\t2_reproduction.json` columns=['n_PartialReprog', 'n_NonReprog', 'md_mean_PartialReprog', 'md_mean_NonReprog', 'md_Partial_minus_NonReprog', 'n_cells_by_label', 'min_cells', 'has_both_labels_ge_min', 'md_partial_lower', 'matches_reported_direction', 'criterion']
- **md4_cell_scores_GM00731** source=`<repo>\data\processed\md4\cell_scores_GM00731.npz` columns=['frozen_ruler', 'md_score', 'tgfb_score', 'age_up', 'age_down', 'age_up_minus_age_down', 'pluri_primary', 'cluster', 'day', 'label']
- **md4_cell_scores_GM23815** source=`<repo>\data\processed\md4\cell_scores_GM23815.npz` columns=['frozen_ruler', 'md_score', 'tgfb_score', 'age_up', 'age_down', 'age_up_minus_age_down', 'pluri_primary', 'cluster', 'day', 'label']
- **t1_cell_counts** source=`<repo>\results\plane\t1_cell_counts.csv` columns=['cell_line', 'label', 'n_clusters', 'n_d0', 'n_d3', 'n_d7', 'n_d10', 'n_all', 'ge30', 'scored', 'reason']
- **t2_louvain_obs_GM00731** source=`<repo>\results\md2\louvain_obs_GM00731.csv` columns=['barcode', 'gsm', 'file', 'cell_line', 'day', 'age_years', 'age_source', 'umi', 'n_genes', 'mito_frac', 'cluster', 'louvain_kept', 'percent_mt']
- **t2_labels_GM00731** source=`<repo>\results\md2\t2_cluster_labels_GM00731.csv` columns=['cell_line', 'cluster', 'n_cells', 'label', 'mean_Fibroblast', 'mean_PartialReprog', 'mean_EarlyPluripotency', 'mean_Pluripotency', 'mean_NonReprog']
- **t2_louvain_obs_GM23815** source=`<repo>\results\md2\louvain_obs_GM23815.csv` columns=['barcode', 'gsm', 'file', 'cell_line', 'day', 'age_years', 'age_source', 'umi', 'n_genes', 'mito_frac', 'cluster', 'louvain_kept', 'percent_mt']
- **t2_labels_GM23815** source=`<repo>\results\md2\t2_cluster_labels_GM23815.csv` columns=['cell_line', 'cluster', 'n_cells', 'label', 'mean_Fibroblast', 'mean_PartialReprog', 'mean_EarlyPluripotency', 'mean_Pluripotency', 'mean_NonReprog']

## Gene-set names actually read

- **MD_built_from_md2** name=`MD (STAR Methods built; md2)` n=205 id_type=symbol source=`<repo>\results\md2\genesets.json`
- **mmc3_Aging_signatures_Age up** name=`Aging_signatures:Age up` n=1533 id_type=symbol source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx via md4/genesets.json`
- **mmc3_Aging_signatures_Age down** name=`Aging_signatures:Age down` n=2007 id_type=symbol source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx via md4/genesets.json`
- **PLURI_ENDOGENOUS** name=`FINDINGS_FIBRO.md PLURI_ENDOGENOUS` n=7 id_type=symbol source=`fibro_common.PLURI_ENDOGENOUS`
- **FIBRO_IDENTITY** name=`FINDINGS_FIBRO.md FIBRO_IDENTITY` n=7 id_type=symbol source=`fibro_common.FIBRO_IDENTITY`

## Files

- `src/plane_common.py`, `src/plane_task1.py`, `src/plane_nulls.py`, `src/plane_task2.py`, `src/plane_task3.py`, `src/plane_findings.py`, `src/plane_run.py`
- `results/plane/`
- `results/plane/figures/plane_GM00731.png`, `results/plane/figures/plane_GM23815.png`
- `FINDINGS_PLANE.md`
- `PROGRESS_PLANE.md`

Existing `FINDINGS_*.md` and `FALSIFICATION.md` were not modified.

## Supersession note, appended 2026-10-02: atlas ruler − MD Δρ on 93 pseudobulks

Appended; nothing above is changed. Where this file says that the FINDINGS_MD3.md Task 2 Δρ result stands, its fibroblast-atlas ruler − MD interval no longer stands.

The ruler − MD Δρ computed with each of the 93 CELLxGENE atlas pseudobulks (a19d1667) as a unit (+0.202, 95% CI +0.006 to +0.440, `excludes_zero=True`; `results/md3/t2_delta_rho.csv`, row `cxg_ruler_minus_MD`) is **superseded** by the per-person re-analysis, in which pseudobulk scores are averaged within each of 65 donors (+0.171, 95% CI −0.118 to +0.412, `excludes_zero=False`; `results/paper_figs/cxg_true_donor_results.csv`, `analysis=true_donor`, `contrast=ruler_minus_MD_AddModuleScore`). The manuscript reports the per-person value in Table 2 (`paper/rejuvenation_readouts_preprint_11.docx`). **The per-person re-analysis was not pre-registered:** `results/paper_figs/` has no PREREG flag, and its results were written on 2026-09-23, after the pre-registered 93-pseudobulk result was known. After regrouping, the ruler − MD interval excludes zero only in GSE226189 (+0.467, 95% CI +0.239 to +0.696). The atlas ruler − (age-up − age-down) interval excludes zero in both units (+0.295, CI +0.107 to +0.487 on 93 pseudobulks; +0.372, CI +0.094 to +0.595 on 65 donors). The pre-registration flags are not edited. The reading is listed as B7 in `results/verify/WITHDRAWN_READINGS.md`.

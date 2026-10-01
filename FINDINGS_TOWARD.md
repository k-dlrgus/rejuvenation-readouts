# FINDINGS_TOWARD — do partially-reprogrammed cells move TOWARD young fibroblasts?

**Status:** DONE. Claim donor GM00731 key=`away_from_old_only`. Seed `20260914`. boot `20260918`. n_perm=200. n_boot=200. n_random=200. C3 pass=True. Frozen ruler exists=True.

Does not modify any existing `FINDINGS_*.md` or `FALSIFICATION.md`. No GSE325735. Uncalibrated R² is never a gate. Donors never pooled. Frozen ruler is a context column only.

Reproduced by `src/toward_run.py`. Frozen ruler: `<repo>\results\fibro\frozen_ruler_ridge_raw.npz` exists=True.

## Pre-registration (verbatim, written before any TOWARD statistic)

Written **before any TOWARD statistic**, 2026-09-21. No threshold in this block is re-tuned after numbers exist.

Question: do partially-reprogrammed cells move TOWARD young fibroblasts, or merely AWAY from old ones?

Primary metric is the cosine cos_S defined below. Uncalibrated R² is never a gate. Spearman ρ is not the statistic here. Frozen ruler score is a context column only. GSE325735 is out of scope. Donors GM00731 and GM23815 are never pooled. States are never averaged. YOUNG/OLD bin boundaries are not re-chosen. Gene space is not switched to PCs. MIN_CELLS=30 is not relaxed.

Seeds: 20260914 (N1 permutation of u / N2 donor splits / C3 fold shuffle / n_random). 20260918 (bootstrap). n_perm=200. n_boot=200. n_random=200.

Inputs (read from disk; refit nothing, rescore nothing, re-cluster nothing):
- GTEx V10 cultured-fibroblast Stage 0 matrix already on disk as data/processed/fibro/stage1_fibro.npz plus results/fibro/stage1_obs.csv (SMAFRZE==RNASEQ, SMTSD=="Cells - Cultured fibroblasts"; gene filter ≥6 counts in ≥20% samples; TMM/edgeR log2-CPM prior.count=2). This is the only GTEx matrix used. It is not rebuilt from the GCT.
- Frozen ruler results/fibro/frozen_ruler_ridge_raw.npz — used only (i) to fix the overlap gene set and the GTEx training μ/σ, and (ii) to report a context column (Z @ w). It is NOT the statistic. Do not refit it.
- GSE297234 labels and scores on disk: results/md2/t2_cluster_labels_*.csv, results/md2/louvain_obs_*.csv, data/processed/md4/cell_scores_*.npz, data/processed/md4/louvain_rulerY_*.npz. Labels are the md2 argmax assignment. Do NOT re-label, re-run Louvain, or re-tune AddModuleScore. cell_scores labels must match md2 labels or STOP.
- results/md4/t2_extrap.csv (PCA k=50 nn_euclidean / Mahalanobis). Read as-is. Do not recompute. Do not average across days.

Anchors (frozen to results/toward/anchors.npz before any GSE297234 vector is computed):
- YOUNG = AGE bins 20-29 and 30-39 (n=113 expected; report actual).
- OLD   = AGE bins 60-69 and 70-79 (n=231 expected; report actual).
- Middle 40-49 and 50-59 is excluded from both anchors and is reported, not dropped silently.
- If either anchor has n<50, STOP and report.
- Gene space = frozen-ruler overlap gene set, unchanged. Per-gene z-score using the frozen GTEx training μ/σ only. Missing genes at z=0.
- c_young = mean z-vector of YOUNG donors. c_old = mean z-vector of OLD donors.
- u = (c_young − c_old) / ||c_young − c_old||  (the aging-axis anchor; unit length).

Displacement, per donor (never pooled), states Fibroblast / PartialReprog / EarlyPluripotency / Pluripotency / NonReprog:
- Timepoints pooled within each destination state. Report timepoint composition.
- Fibroblast pseudobulk uses d0 cells only and is the ORIGIN. Fibroblast is not a destination. Non-d0 Fibroblast cells are reported as unused, not a destination.
- If a donor has <30 d0 Fibroblast cells, STOP for that donor.
- A destination state qualifies iff n_cells ≥ 30 in that donor. If it does not qualify, it is reported as skipped with its n, not dropped.
- TMM among that donor's panel (origin + qualifying destinations), then log2-CPM prior.count=2, then frozen μ/σ.
- For each non-origin qualifying state S: d_S = z(S) − z(Fibroblast_d0).

Primary statistic:
  cos_S = dot(d_S, u) / ||d_S||
(u is unit, so this is cosine alignment of the move with the young-minus-old axis.)

Reported alongside, never substituted for cos_S: ||d_S||; dist_young_S, dist_young_origin, delta_young = dist_young_S − dist_young_origin (negative = closer to young); dist_old_S, dist_old_origin, delta_old (positive = farther from old); frozen-ruler score of each pseudobulk (context); nn_euclidean and Mahalanobis for that state from md4 t2_extrap.csv as-is by day (not averaged).

"TOWARD young" for a state REQUIRES ALL of: cos_S > 0, p_N1 ≤ 0.05, p_N2 ≤ 0.05, AND delta_young < 0. Either cosine-without-closer or closer-without-cosine is "AWAY from old only" / not toward-young. Passing one null is not passing. Both N1 and N2 required. delta_old > 0 alone is not evidence of toward-young.

Nulls (not interchangeable), both reported for every qualifying destination state:
- N1 (direction null, primary): 200 random unit vectors u_r by permuting the entries of u across genes (seed 20260914, independent Generator). p = (#{cos(d_S, u_r) ≥ cos_S} + 1) / 201.
- N2 (anchor null): 200 splits of ALL GTEx cultured-fibroblast donors in the Stage 0 matrix into two equal-size pseudo-groups at random IGNORING age, rebuilding u each time (seed 20260914, independent Generator). Remainder of 1 donor if n odd is dropped from both groups that draw. p as in N1.
PASS_NULLS = p_N1 ≤ 0.05 AND p_N2 ≤ 0.05.

Bootstrap (seed 20260918, independent Generator per donor with offset 0 for GM00731 and 1 for GM23815): B=200 cell resamples with replacement within origin and within each qualifying destination independently; re-sum, re-TMM the donor panel, re-apply frozen μ/σ; 95% percentile CI on cos_S and delta_young. Anchors stay frozen. Not a gate.

C3 positive control (GTEx only; run after anchors.npz is written and BEFORE any GSE297234 vector):
- Shuffle all Stage 0 donors with seed 20260914 (independent Generator), assign fold = i % 5.
- Each fold: rebuild YOUNG/OLD/u on the training donors only, same bins, same frozen μ/σ (μ/σ are not refit per fold).
- Origin = mean z of held-out OLD donors. Destination = c_young of the training fold. d = c_young_train − mean(z_old_test). cos_fold = dot(d, u_train) / ||d||.
- Construction minimum, frozen: n_young_train ≥ 20, n_old_train ≥ 20, n_old_test ≥ 5 in every fold. If any fold fails construction, C3 fails.
- C3 pass = all 5 folds finite AND every fold cos > 0.50 (strongly positive). This bar is not re-tuned.
- If C3 fails, the statistic is broken. STOP. Nothing downstream (no GSE vectors, no reading) is reported.

C1 Pluripotency (per donor, after GSE vectors): if Pluripotency qualifies and meets TOWARD young (cos>0, both nulls, delta_young<0), the statistic is reading displacement magnitude / extrapolation, not direction toward young. Then PartialReprog is uninterpretable regardless of how it looks.
C2 NonReprog: off-trajectory. Reported, not a gate.
If Pluripotency is skipped (n<30), C1 cannot be evaluated; do not fire toward-young; key = c1_unevaluable.

Pre-registered reading (only the outcome that fired; per donor). GM23815 is a contrast, never the claim-supporting reading. Decision order:
1. C3 fails → key `c3_broken`. STOP. Nothing downstream is reported.
2. Donor origin STOP → key `origin_stop`.
3. PartialReprog skipped → key `skipped_partial`.
4. PartialReprog has cos_S > 0 and delta_young ≥ 0 → key `away_from_old_only`. This is a real negative finding: movement is away-from-old only. A 1D age score cannot distinguish these two cases. Do not fall back to the frozen-ruler score.
5. PartialReprog meets TOWARD young AND Pluripotency also meets TOWARD young → key `uninterpretable_c1`. Do NOT report a toward-young result. What would fix it: a direction statistic that Pluripotency does not pass, or an in-manifold target that is not the young-fibroblast centroid in the full overlap space.
6. PartialReprog meets TOWARD young AND Pluripotency does not (C1 clean) AND C3 passed → key `toward_young`. Partially-reprogrammed cells move toward the young fibroblast target, not merely away from the old state. The safe-target framing is supported in this dataset. State n_cells in the same sentence.
7. PartialReprog does not pass both nulls → key `no_directional_structure`. If both p>0.05, say so. If exactly one p≤0.05, say passing one null is not passing. Do not fall back to the frozen-ruler score.

What does not count: frozen-ruler movement (label any such sentence context); large ||d_S||; delta_old>0 alone; one of two nulls; averaging donors or states; re-choosing bins; cosine-on-PCs; relaxing MIN_CELLS=30.

Flag: `results/toward/PREREG.flag` exists=True.

## STOP / failures

Not substituting columns or repairing rows.

None recorded.

## Anchors (frozen before any GSE297234 vector)

- YOUNG bins ('20-29', '30-39'): n=113 (expected 113).
- OLD bins ('60-69', '70-79'): n=231 (expected 231).
- Middle bins ('40-49', '50-59'): n=308 excluded from both anchors, not dropped silently.
- Gene space n=23485 (frozen-ruler overlap). ||c_young − c_old||=+32.2839.
- Written to `<repo>\results\toward\anchors.npz`.

| AGE_bin | n | role |
|---|---|---|
| 20-29 | 60 | YOUNG |
| 30-39 | 53 | YOUNG |
| 40-49 | 106 | middle_excluded |
| 50-59 | 202 | middle_excluded |
| 60-69 | 210 | OLD |
| 70-79 | 21 | OLD |
| YOUNG_total | 113 | YOUNG |
| OLD_total | 231 | OLD |
| middle_total | 308 | middle_excluded |

## C3 positive control (GTEx held-out OLD → train YOUNG centroid)

- C3 pass=True (all 5 folds finite and cos > 0.5; construction_ok=True). min_cos=+0.615 mean_cos=+0.775 max_cos=+0.874.
- Bootstrap 95% CI (seed 20260918) resamples held-out OLD donors within fold; not a gate.

| fold | n_train | n_test | n_young_train | n_old_train | n_old_test | cos | d_norm | cos_ci_lo | cos_ci_hi | pass_bar | construction_ok | reason |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 521 | 131 | 94 | 187 | 44 | +0.615 | +31.425 | +0.123 | +0.754 | True | True | NA |
| 1 | 521 | 131 | 92 | 181 | 50 | +0.849 | +44.373 | +0.701 | +0.835 | True | True | NA |
| 2 | 522 | 130 | 86 | 180 | 51 | +0.714 | +35.448 | +0.300 | +0.821 | True | True | NA |
| 3 | 522 | 130 | 87 | 191 | 40 | +0.825 | +38.051 | +0.457 | +0.846 | True | True | NA |
| 4 | 522 | 130 | 93 | 185 | 46 | +0.874 | +41.185 | +0.677 | +0.844 | True | True | NA |

## Timepoint composition per state per donor

Labels: md2 argmax of mean AddModuleScore (mmc3 Reprog_cell_state_signatures). Not re-labelled. Fibroblast origin uses d0 only. Non-d0 Fibroblast cells are unused, not a destination. Destination states pool timepoints. MIN_CELLS=30. No destination state was skipped.

| cell_line | state | n_cells | n_d0 | n_d3 | n_d7 | n_d10 | n_d0_fibroblast | qualifies_dest | origin_ok | unused_fibroblast_nond0 | skipped_reason |
|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | Fibroblast | 5988 | 4939 | 5 | 14 | 1030 | 4939 | False | True | 1049 | — |
| GM00731 | PartialReprog | 5832 | 3 | 5220 | 578 | 31 | — | True | False | 0 | — |
| GM00731 | EarlyPluripotency | 4414 | 2 | 410 | 3917 | 85 | — | True | False | 0 | — |
| GM00731 | Pluripotency | 439 | 0 | 5 | 354 | 80 | — | True | False | 0 | — |
| GM00731 | NonReprog | 5459 | 65 | 348 | 1932 | 3114 | — | True | False | 0 | — |
| GM23815 | Fibroblast | 11355 | 7703 | 46 | 2233 | 1373 | 7703 | False | True | 3652 | — |
| GM23815 | PartialReprog | 6481 | 9 | 6113 | 309 | 50 | — | True | False | 0 | — |
| GM23815 | EarlyPluripotency | 12961 | 56 | 316 | 8172 | 4417 | — | True | False | 0 | — |
| GM23815 | Pluripotency | 1249 | 0 | 91 | 977 | 181 | — | True | False | 0 | — |
| GM23815 | NonReprog | 1812 | 0 | 10 | 8 | 1794 | — | True | False | 0 | — |

## Statistic table (primary = cos_S)

cos_S = dot(d_S, u) / ||d_S||. delta_young < 0 means closer to the YOUNG centroid. Frozen-ruler score is a **context** column only. nn_euclidean and Mahalanobis are read as-is from results/md4/t2_extrap.csv by day and are not averaged.

### Primary and accompanying numbers

| cell_line | state | role | n_cells | n_d0 | n_d3 | n_d7 | n_d10 | skipped | cos_S | d_norm | dist_young | dist_young_origin | delta_young | dist_old | dist_old_origin | delta_old | ruler_score_context | p_N1 | p_N2 | pass_nulls | meets_toward | cos_ci_lo | cos_ci_hi | delta_young_ci_lo | delta_young_ci_hi |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | Fibroblast | origin | 4939 | 4939 | 5 | 14 | 1030 | False | NA | +0.000 | +467.450 | +467.450 | +0.000 | +464.431 | +464.431 | +0.000 | +3.222 | NA | NA | False | False | NA | NA | NA | NA |
| GM00731 | PartialReprog | destination | 5832 | 3 | 5220 | 578 | 31 | False | +0.077 | +271.686 | +476.820 | +467.450 | +9.370 | +475.288 | +464.431 | +10.857 | +0.650 | 0.0050 | 0.3881 | False | False | +0.073 | +0.079 | +7.722 | +10.125 |
| GM00731 | EarlyPluripotency | destination | 4414 | 2 | 410 | 3917 | 85 | False | +0.117 | +268.580 | +497.906 | +467.450 | +30.456 | +497.126 | +464.431 | +32.695 | -2.621 | 0.0050 | 0.0796 | False | False | +0.113 | +0.119 | +28.992 | +31.568 |
| GM00731 | Pluripotency | destination | 439 | 0 | 5 | 354 | 80 | False | +0.126 | +546.420 | +651.500 | +467.450 | +184.049 | +652.744 | +464.431 | +188.313 | -10.402 | 0.0050 | 0.1393 | False | False | +0.113 | +0.120 | +192.650 | +212.674 |
| GM00731 | NonReprog | destination | 5459 | 65 | 348 | 1932 | 3114 | False | +0.063 | +192.764 | +476.382 | +467.450 | +8.932 | +474.249 | +464.431 | +9.818 | +6.855 | 0.0050 | 0.3085 | False | False | +0.056 | +0.065 | +8.058 | +10.075 |
| GM23815 | Fibroblast | origin | 7703 | 7703 | 46 | 2233 | 1373 | False | NA | +0.000 | +457.031 | +457.031 | +0.000 | +454.287 | +454.287 | +0.000 | -1.893 | NA | NA | False | False | NA | NA | NA | NA |
| GM23815 | PartialReprog | destination | 6481 | 9 | 6113 | 309 | 50 | False | +0.052 | +260.471 | +471.220 | +457.031 | +14.189 | +469.497 | +454.287 | +15.210 | -4.659 | 0.0050 | 0.4627 | False | False | +0.049 | +0.055 | +12.780 | +14.675 |
| GM23815 | EarlyPluripotency | destination | 12961 | 56 | 316 | 8172 | 4417 | False | +0.076 | +272.149 | +503.353 | +457.031 | +46.322 | +502.190 | +454.287 | +47.903 | -7.959 | 0.0050 | 0.2139 | False | False | +0.073 | +0.078 | +44.739 | +46.664 |
| GM23815 | Pluripotency | destination | 1249 | 0 | 91 | 977 | 181 | False | +0.112 | +515.120 | +642.989 | +457.031 | +185.957 | +643.944 | +454.287 | +189.656 | -15.095 | 0.0050 | 0.1493 | False | False | +0.102 | +0.110 | +189.733 | +206.295 |
| GM23815 | NonReprog | destination | 1812 | 0 | 10 | 8 | 1794 | False | +0.010 | +219.796 | +499.003 | +457.031 | +41.971 | +496.633 | +454.287 | +42.345 | -0.098 | 0.0398 | 0.4478 | False | False | -0.003 | +0.009 | +48.519 | +52.503 |

### Context columns (frozen ruler; md4 t2_extrap as-is; not substituted for cos_S)

| cell_line | state | ruler_score_context | nn_euclidean_by_day | mahalanobis_by_day | skipped_reason |
|---|---|---|---|---|---|
| GM00731 | Fibroblast | +3.222 | d0:465.617(n=4939); d3:1215.554(n=5); d7:918.913(n=14); d10:474.154(n=1030) | d0:24.411(n=4939); d3:56.410(n=5); d7:38.900(n=14); d10:24.294(n=1030) | origin d0 n=4939; unused non-d0 Fibroblast n=1049 |
| GM00731 | PartialReprog | +0.650 | d0:1342.161(n=3); d3:479.017(n=5220); d7:468.982(n=578); d10:709.005(n=31) | d0:64.592(n=3); d3:23.982(n=5220); d7:23.158(n=578); d10:31.008(n=31) | NA |
| GM00731 | EarlyPluripotency | -2.621 | d0:2192.731(n=2); d3:511.647(n=410); d7:503.349(n=3917); d10:531.941(n=85) | d0:137.899(n=2); d3:24.076(n=410); d7:24.453(n=3917); d10:24.800(n=85) | NA |
| GM00731 | Pluripotency | -10.402 | d3:1553.696(n=5); d7:616.361(n=354); d10:631.384(n=80) | d3:80.698(n=5); d7:24.829(n=354); d10:25.461(n=80) | NA |
| GM00731 | NonReprog | +6.855 | d0:670.751(n=65); d3:476.754(n=348); d7:470.259(n=1932); d10:485.365(n=3114) | d0:26.971(n=65); d3:22.042(n=348); d7:24.100(n=1932); d10:24.012(n=3114) | NA |
| GM23815 | Fibroblast | -1.893 | d0:449.237(n=7703); d3:593.699(n=46); d7:465.875(n=2233); d10:468.409(n=1373) | d0:23.054(n=7703); d3:25.444(n=46); d7:23.737(n=2233); d10:24.716(n=1373) | origin d0 n=7703; unused non-d0 Fibroblast n=3652 |
| GM23815 | PartialReprog | -4.659 | d0:1026.072(n=9); d3:455.379(n=6113); d7:482.403(n=309); d10:566.537(n=50) | d0:47.539(n=9); d3:23.014(n=6113); d7:24.240(n=309); d10:26.410(n=50) | NA |
| GM23815 | EarlyPluripotency | -7.959 | d0:854.835(n=56); d3:514.939(n=316); d7:508.372(n=8172); d10:484.765(n=4417) | d0:35.877(n=56); d3:24.337(n=316); d7:25.420(n=8172); d10:24.461(n=4417) | NA |
| GM23815 | Pluripotency | -15.095 | d3:642.148(n=91); d7:627.900(n=977); d10:682.821(n=181) | d3:24.855(n=91); d7:25.657(n=977); d10:27.762(n=181) | NA |
| GM23815 | NonReprog | -0.098 | d3:1066.480(n=10); d7:1515.850(n=8); d10:480.957(n=1794) | d3:48.430(n=10); d7:80.299(n=8); d10:24.440(n=1794) | NA |

## N1 / N2 p-values

N1: permute entries of frozen u (seed 20260914), n_random=200. N2: equal-size GTEx donor splits ignoring age (seed 20260914), n_perm=200. p = (# null ≥ observed + 1) / 201. PASS_NULLS requires both p≤0.05.

| cell_line | state | n_cells | cos_S | p_N1 | p_N2 | pass_nulls | delta_young | meets_toward |
|---|---|---|---|---|---|---|---|---|
| GM00731 | PartialReprog | 5832 | +0.077 | 0.0050 | 0.3881 | False | +9.370 | False |
| GM00731 | EarlyPluripotency | 4414 | +0.117 | 0.0050 | 0.0796 | False | +30.456 | False |
| GM00731 | Pluripotency | 439 | +0.126 | 0.0050 | 0.1393 | False | +184.049 | False |
| GM00731 | NonReprog | 5459 | +0.063 | 0.0050 | 0.3085 | False | +8.932 | False |
| GM23815 | PartialReprog | 6481 | +0.052 | 0.0050 | 0.4627 | False | +14.189 | False |
| GM23815 | EarlyPluripotency | 12961 | +0.076 | 0.0050 | 0.2139 | False | +46.322 | False |
| GM23815 | Pluripotency | 1249 | +0.112 | 0.0050 | 0.1493 | False | +185.957 | False |
| GM23815 | NonReprog | 1812 | +0.010 | 0.0398 | 0.4478 | False | +41.971 | False |

## C1 Pluripotency / C2 NonReprog

- **GM00731 C1 Pluripotency:** skipped=False meets_toward=False cos_S=+0.126 delta_young=+184.049 p_N1=0.0050 p_N2=0.1393 n_cells=439 d_norm=+546.420.
- **GM00731 C2 NonReprog (not a gate):** skipped=False meets_toward=False cos_S=+0.063 delta_young=+8.932 p_N1=0.0050 p_N2=0.3085 n_cells=5459.
- **GM23815 C1 Pluripotency:** skipped=False meets_toward=False cos_S=+0.112 delta_young=+185.957 p_N1=0.0050 p_N2=0.1493 n_cells=1249 d_norm=+515.120.
- **GM23815 C2 NonReprog (not a gate):** skipped=False meets_toward=False cos_S=+0.010 delta_young=+41.971 p_N1=0.0398 p_N2=0.4478 n_cells=1812.

## Pre-registered reading (only the outcome that fired)

**Claim donor (GM00731) key: `away_from_old_only`.** GM23815 is a contrast, never the reading. Cosines are not averaged.

GM00731: PartialReprog has cos_S=+0.077 but delta_young=+9.370 (n_cells=5832). Movement is away-from-old only. This is a real negative finding: a 1D age score cannot distinguish toward-young from away-from-old. The frozen-ruler score is not used to rescue this.

### GM00731 (claim) key=`away_from_old_only`

GM00731: PartialReprog has cos_S=+0.077 but delta_young=+9.370 (n_cells=5832). Movement is away-from-old only. This is a real negative finding: a 1D age score cannot distinguish toward-young from away-from-old. The frozen-ruler score is not used to rescue this.

### GM23815 (contrast) key=`away_from_old_only`

CONTRAST (young donor GM23815; never the reading): PartialReprog has cos_S=+0.052 but delta_young=+14.189 (n_cells=6481). Movement is away-from-old only. This is a real negative finding: a 1D age score cannot distinguish toward-young from away-from-old. The frozen-ruler score is not used to rescue this.

## Limitations

- GTEx public AGE is a 10-year bin; YOUNG and OLD are bin unions, not chronological years.
- GTEx cultured fibroblasts are not the same system as GSE297234 primary dermal fibroblasts.
- Two donors. GM23815 is a contrast and is not pooled into the claim.
- State labels are Louvain-cluster argmax of signature scores, not per-cell argmax.
- TMM is among a handful of state-pseudobulk rows per donor.
- The aging-axis anchor lives in the full frozen-ruler overlap space; states that leave the fibroblast manifold (Pluripotency) can produce large displacements that are not rejuvenation. C1 exists to catch that.
- nn_euclidean / Mahalanobis are copied from md4 t2_extrap.csv (state × timepoint) and are not recomputed on the pooled-state vector.
- Day-0 fibroblast origins already sit hundreds of z-units from both GTEx centroids (scRNA pseudobulk vs bulk cultured fibroblasts). Subsequent states recede from both anchors; delta_young and delta_old are therefore not a within-cloud slide.
- Bootstrap CIs re-TMM the donor panel, so they are CIs on the TMM-coupled statistic and need not contain the observed point estimate. Not a gate.
- Uncalibrated R² is never a gate and is not computed here.
- GSE325735 is out of scope and is not opened.

## Files touched

- `src/toward_run.py`
- `FINDINGS_TOWARD.md`
- `PROGRESS_TOWARD.md`
- `results/toward/PREREG.flag` (not rewritten if it already existed)
- `results/toward/anchors.npz`
- `results/toward/anchors.json`
- `results/toward/anchor_table.csv`
- `results/toward/null_directions.npz`
- `results/toward/c3_folds.csv`
- `results/toward/c3_summary.json`
- `results/toward/occupancy.csv`
- `results/toward/stats.csv`
- `results/toward/state_z.npz`
- `results/toward/reading.json`
- `results/toward/summary.json`
- `results/toward/manifest.json`
- `results/toward/run_report.txt`
- Existing `FINDINGS_*.md`, `PROGRESS_*.md`, and `FALSIFICATION.md` were not modified (except FINDINGS_TOWARD.md and PROGRESS_TOWARD.md).


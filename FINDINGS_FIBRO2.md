# FINDINGS_FIBRO2 — is the fibroblast ruler the instrument, or is the result real?

**Status:** T-A: not interpretable (`primary_unscored`). T-B: d10 not a material technical difference (`d7_unexplained`). T-C: linearity not the limiting factor (`boosting_below_ridge`). Seed `20260914`. boot `20260918`. n_perm=200. n_boot=200. n_random=200. Frozen ruler `frozen_ruler_ridge_raw.npz` exists=True.

Does not modify any existing `FINDINGS_*.md` or `FALSIFICATION.md`. No GSE325735. Uncalibrated R² is reported and is never a gate. Nothing averaged across cohorts, donors, or regimes. T-A / T-B / T-C cannot overturn the Stage 2 `no_decline` verdict by themselves; they establish whether it is interpretable.

Reproduced by `src/fibro2_run.py`. Frozen ruler: `<repo>\results\fibro\frozen_ruler_ridge_raw.npz` exists=True.

Flag: `results/fibro2/PREREG_TA.flag` exists=True. `PREREG_TB.flag` exists=True. `PREREG_TC.flag` exists=True.

## T-A pre-registration (verbatim, written before any T-A download or ρ)

Written **before any T-A download, projection, or ρ**, 2026-09-18. No threshold in this block is re-tuned after numbers exist.

Take the Stage 1 ridge direction **exactly as frozen in `results/fibro/frozen_ruler_ridge_raw.npz`**. Refit nothing. Recalibrate nothing. Rescale nothing. Project, score, report Spearman ρ vs donor age.

Same normalization as the frozen GTEx spec (FINDINGS_GTEX.md Stage 0 / FINDINGS_FIBRO.md Stage 1): TMM/edgeR log2-CPM (`prior.count=2`). Genes restricted to the overlap with the frozen ruler. Missing overlap genes left at z=0 (FINDINGS_EXTERNAL.md rule). Gene-overlap count reported. z uses the frozen train mean and sd.

Metric: Spearman ρ vs the age the record states. Donor-level permutation null n_perm=200 (scores held fixed; shuffle donor ages; each cell starts a fresh Generator on seed `20260914`). Donor-bootstrap 95% percentile CI, B=200, seed `20260918`. Calibrated R² and uncalibrated R² reported; uncalibrated R² is never a gate. Pearson r reported, not a gate. Nothing averaged across cohorts.

Adult restriction on the primary cell: GSE113957 samples with stated age < 18 are excluded to match the GTEx training range; the excluded count is reported. Do not infer a donor age that the record does not state. A record that does not state per-sample donor age is unusable.

Cohorts, in priority order, each resolved from the record this session:
1. **GSE113957** — primary external cell. Human dermal fibroblast. Ages stated in titles. Bulk. Adults ≥18.
2. **Any single-cell fibroblast/skin aging cohort with ≥10 donors and stated ages** resolved this session (CELLxGENE, GEO, or SCOPING_TISSUE.md `near_miss` rows GSE226189 n=82 bulk, GSE307377 n=9). A 10x cohort is worth more than a larger bulk one because the platform shift is the thing under test. Say which platform each is.
3. If no single-cell fibroblast cohort resolves, say so plainly. Do not substitute a bulk cohort and call the platform question answered.

GSE325735 is out of scope.

**Pre-registered reading** (only the outcome that fired; T-A cannot overturn the Stage 2 `no_decline` verdict by itself):
- ρ > 0 with null ≤ 0.05 in the primary external cell → the ruler is a working instrument outside its training cohort; the Stage 2 `no_decline` verdict is interpretable as a statement about OSKM.
- ρ ≤ 0 or null > 0.05 → the ruler does not validate externally; Stage 2 is uninterpretable and `no_decline` must not be reported as a finding about reprogramming.
- bulk external passes, single-cell external fails (or none resolves) → the platform shift is unresolved; the Stage 2 projection rests on an unvalidated cross-platform step; do not resolve it by argument.


## T-A — cohort table

| accession | role | title | platform_kind | types | platforms | n_gsm | n_with_age | n_missing_age | n_adult | n_child | n_donors_with_age | usable | unusable_reason |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GSE113957 | primary_external | Predicting age from the transcriptome of human dermal fibroblasts | bulk_rnaseq | Expression profiling by high throughput sequencing | GPL16791; GPL18573 | 143 | 135 | +8.000 | +111.000 | +24.000 | NA | False | FPKM not integer counts (GSE113957_fpkm.txt.gz); series matrix 404 |
| a19d1667-a7b5-4556-9e5f-f9bfa690c0f1 | single_cell_external | Human dermal fibroblast atlas | single_cell | single-cell RNA-seq | NA | 108440 | 93 | NA | NA | NA | NA | True | NA |
| GSE226189 | near_miss_bulk | Transcriptomes of human primary skin fibroblasts reveal age-associated mRNAs and long noncoding RNAs | bulk_rnaseq | Expression profiling by high throughput sequencing; Non-coding RNA profiling by high throughput sequencing | GPL24676 | 82 | 82 | +0.000 | NA | NA | +82.000 | True | NA |
| GSE307377 | near_miss_n9 | Bulk RNA-seq analysis of primary human dermal fibroblasts from young and aged donors | bulk_rnaseq | Expression profiling by high throughput sequencing | GPL24676 | 9 | 9 | +0.000 | NA | NA | +1.000 | False | HOMER tags not integer counts; series matrix 0 genes |

## T-A — results table

| accession | role | platform_kind | tag | n_samples | n_donors | n_excluded_no_age | n_excluded_child | n_overlap | n_missing_z0 | rho | rho_null | rho_p | rho_ci_lo | rho_ci_hi | cal_r2 | r2 | r | pass_rho | n_perm | n_boot | used_10x_only | platform_detail |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| a19d1667-a7b5-4556-9e5f-f9bfa690c0f1 | single_cell_external | single_cell | cxg_a19d1667 | 93 | 93 | 0 | 0 | 17047 | 6580 | +0.455 | -0.003 | +0.005 | +0.253 | +0.605 | +0.171 | -9.664 | +0.414 | True | 200 | 200 | True | 10x (restricted) |
| GSE226189 | near_miss_bulk | bulk_rnaseq | GSE226189_near_miss_bulk | 82 | 82 | 0 | 0 | 21351 | 2134 | +0.549 | +0.005 | +0.005 | +0.374 | +0.667 | +0.299 | -9.963 | +0.547 | True | 200 | 200 | NA | NA |

## T-A — reading

Fired key: `primary_unscored`. interpretable=False.

primary external cell (GSE113957) was not scored: GSE113957_fpkm.txt.gz values are not integer counts (cannot apply TMM/log2-CPM). min=0.0 max=249213345.0 mean=1009452.0862567581. Series matrix HEAD 404. n_gsm=143 n_with_stated_age=135 n_missing_age=8 n_stated_age_<18=24 n_adult>=18=111. The pre-registered primary-cell pass/fail on rho therefore cannot fire. A single-cell external cohort (CELLxGENE Human dermal fibroblast atlas a19d1667-a7b5-4556-9e5f-f9bfa690c0f1, 10x-restricted adult donors) passed rho=+0.455 null=-0.003 p=+0.005 CI=[+0.253, +0.605] n=93 overlap=17047. It is not the pre-registered primary cell, so the primary-cell interpretable bullet does not fire. Stage 2 no_decline must not be reported as a finding about reprogramming until the primary external cell is scored and passes. A 10x cohort was scored; this is not the case no single-cell cohort resolved. Do not substitute the single-cell pass or the GSE226189 bulk pass for GSE113957.

Single-cell fibroblast cohort resolved this session: **True**.

## T-B pre-registration (verbatim, written before any T-B QC table or d0→d7 null)

Written **before any T-B QC table or d0→d7 null**, 2026-09-18. No threshold in this block is re-tuned after numbers exist.

From the GSE297234 Cell Ranger `filtered_feature_bc_matrix.h5` matrices already on disk under `data/raw/gse297234/RAW/`, report **per donor × timepoint**: cells recovered, median UMI/cell, median genes/cell, mitochondrial read fraction, and n cells per unsupervised cluster. Table only — this is a QC audit, not a model. Unsupervised clusters are the Stage 2 spec: within each donor × timepoint independently, log1p(CP10k) on genes overlapping the frozen ruler, PCA 20 (or n_cells−1 if smaller), k-means k=3, seed `20260914`. Skip a group if n_cells < 50.

Then rerun the Stage 2 **random-direction null exactly as pre-registered** (200 permuted-weight directions, seed `20260914`) on the **d0→d7 endpoint** in addition to the d0→d10 endpoint already reported, for both donors, both pseudobulk variants (all-cell primary, per-cluster secondary). Cluster labels are assigned independently per donor × timepoint; the secondary pairing is by **within-timepoint cell-count rank** (largest / middle / smallest), not by k-means ID. Do not average cluster ranks into the all-cell primary. Do not average donors.

Endpoint decline = age_score(day 0) − age_score(day T). Positive = younger on the ruler. Empirical p = (n_random ≥ real + 1) / (n_random + 1). Pass bar remains p ≤ 0.05. A post-hoc d7 endpoint is not permitted to become the primary Stage 2 result.

**Pre-registered reading** (caveat on the verdict, not a reversal of it):
- d0→d7 also gives p > 0.05 → the `no_decline` verdict does not depend on d10; d10's behaviour is a side observation.
- d0→d7 gives p ≤ 0.05 while d0→d10 does not, **and** the QC table shows d10 differing materially from d0/d3/d7 on cells recovered, depth, or mitochondrial fraction → the endpoint choice is confounded by a technical difference at d10. Report as a caveat. State what would settle it (a cohort with more donors and a pre-registered endpoint).
- d0→d7 gives p ≤ 0.05 with no technical difference at d10 → report the discrepancy and say it is unexplained.


## T-B — QC table (per donor × timepoint)

| cell_line | day | gsm | age_years | n_cells | median_UMI_per_cell | median_genes_per_cell | mito_read_fraction_library | mito_read_fraction_median_cell | n_mito_genes | n_cluster_0 | n_cluster_1 | n_cluster_2 | chemistry | software |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | 0 | GSM8986586 | 96 | 5021 | +14380.000 | +4078.000 | +0.023 | +0.026 | 13 | 73 | 2535 | 2413 | Single Cell 3' v3 | cellranger-7.1.0 |
| GM00731 | 3 | GSM8986587 | 96 | 6318 | +10199.000 | +3646.000 | +0.039 | +0.036 | 13 | 2159 | 915 | 3244 | Single Cell 3' v3 | cellranger-7.1.0 |
| GM00731 | 7 | GSM8986588 | 96 | 6915 | +10922.000 | +3700.000 | +0.033 | +0.031 | 13 | 3054 | 3011 | 850 | Single Cell 3' v3 | cellranger-7.1.0 |
| GM00731 | 10 | GSM8986589 | 96 | 4399 | +13456.000 | +4187.000 | +0.036 | +0.034 | 13 | 2645 | 708 | 1046 | Single Cell 3' v3 | cellranger-7.1.0 |
| GM23815 | 0 | GSM8986590 | 22 | 7782 | +8140.000 | +3030.000 | +0.027 | +0.029 | 13 | 2896 | 3315 | 1571 | Single Cell 3' v3 | cellranger-7.1.0 |
| GM23815 | 3 | GSM8986591 | 22 | 6738 | +10662.000 | +3881.500 | +0.051 | +0.048 | 13 | 2637 | 3346 | 755 | Single Cell 3' v3 | cellranger-7.1.0 |
| GM23815 | 7 | GSM8986592 | 22 | 12209 | +6108.000 | +2565.000 | +0.038 | +0.035 | 13 | 4407 | 6163 | 1639 | Single Cell 3' v3 | cellranger-7.1.0 |
| GM23815 | 10 | GSM8986593 | 22 | 7935 | +9063.000 | +3183.000 | +0.050 | +0.044 | 13 | 2910 | 793 | 4232 | Single Cell 3' v3 | cellranger-7.1.0 |

Mitochondrial rule: `symbol upper startswith MT-` n_mito_genes=13 symbols=['MT-ATP6', 'MT-ATP8', 'MT-CO1', 'MT-CO2', 'MT-CO3', 'MT-CYB', 'MT-ND1', 'MT-ND2', 'MT-ND3', 'MT-ND4', 'MT-ND4L', 'MT-ND5', 'MT-ND6'].

## T-B — random-direction null table

| variant | cell_line | endpoint | primary | rank | n_cells_d0 | n_cells_end | score_d0 | score_end | decline | p | n_random_ge_real | n_random | pass_p | ok |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all_cell | GM00731 | d0→d10 | True | NA | NA | NA | +2.769 | +6.120 | -3.351 | +0.990 | 198 | 200 | False | True |
| all_cell | GM00731 | d0→d7 | True | NA | NA | NA | +2.769 | -1.416 | +4.185 | +0.005 | 0 | 200 | True | True |
| all_cell | GM23815 | d0→d10 | False | NA | NA | NA | -2.423 | -4.885 | +2.462 | +0.060 | 11 | 200 | False | True |
| all_cell | GM23815 | d0→d7 | False | NA | NA | NA | -2.423 | -9.740 | +7.317 | +0.005 | 0 | 200 | True | True |
| cluster_size_rank | GM00731 | d0→d10 | False | +0.000 | +2535.000 | +2645.000 | +5.464 | +8.757 | -3.293 | +0.990 | 198 | 200 | False | True |
| cluster_size_rank | GM00731 | d0→d10 | False | +1.000 | +2413.000 | +1046.000 | +2.584 | +4.973 | -2.389 | +0.955 | 191 | 200 | False | True |
| cluster_size_rank | GM00731 | d0→d10 | False | +2.000 | +73.000 | +708.000 | -4.814 | +2.248 | -7.062 | +0.960 | 192 | 200 | False | True |
| cluster_size_rank | GM00731 | d0→d7 | False | +0.000 | +2535.000 | +3054.000 | +5.464 | +4.406 | +1.058 | +0.224 | 44 | 200 | False | True |
| cluster_size_rank | GM00731 | d0→d7 | False | +1.000 | +2413.000 | +3011.000 | +2.584 | -3.795 | +6.379 | +0.005 | 0 | 200 | True | True |
| cluster_size_rank | GM00731 | d0→d7 | False | +2.000 | +73.000 | +850.000 | -4.814 | -5.240 | +0.427 | +0.224 | 44 | 200 | False | True |
| cluster_size_rank | GM23815 | d0→d10 | False | +0.000 | +3315.000 | +4232.000 | +0.026 | -4.218 | +4.244 | +0.010 | 1 | 200 | True | True |
| cluster_size_rank | GM23815 | d0→d10 | False | +1.000 | +2896.000 | +2910.000 | -0.514 | +0.965 | -1.479 | +0.871 | 174 | 200 | False | True |
| cluster_size_rank | GM23815 | d0→d10 | False | +2.000 | +1571.000 | +793.000 | -1.839 | -8.148 | +6.308 | +0.050 | 9 | 200 | True | True |
| cluster_size_rank | GM23815 | d0→d7 | False | +0.000 | +3315.000 | +6163.000 | +0.026 | -5.693 | +5.719 | +0.005 | 0 | 200 | True | True |
| cluster_size_rank | GM23815 | d0→d7 | False | +1.000 | +2896.000 | +4407.000 | -0.514 | -10.151 | +9.637 | +0.005 | 0 | 200 | True | True |
| cluster_size_rank | GM23815 | d0→d7 | False | +2.000 | +1571.000 | +1639.000 | -1.839 | -13.728 | +11.889 | +0.005 | 0 | 200 | True | True |

## T-B — reading

Fired key: `d7_unexplained`. d10_technical=False. p_d0d7=0.004975124378109453 p_d0d10=0.9900497512437811.

d0→d7 gives p ≤ 0.05 with no technical difference at d10 → the discrepancy is unexplained.

A post-hoc d7 endpoint is not the primary Stage 2 result.

d10 vs d0/d3/d7 (from disk):

- GM00731 cells recovered: d0/d3/d7 median=6318 d10=4399 ratio=0.696 material=False
- GM00731 median UMI/cell: d0/d3/d7 median=1.092e+04 d10=1.346e+04 ratio=1.232 material=False
- GM00731 mitochondrial read fraction: d0/d3/d7 median=0.03317 d10=0.03588 ratio=1.082 material=False
- GM23815 cells recovered: d0/d3/d7 median=7782 d10=7935 ratio=1.020 material=False
- GM23815 median UMI/cell: d0/d3/d7 median=8140 d10=9063 ratio=1.113 material=False
- GM23815 mitochondrial read fraction: d0/d3/d7 median=0.03839 d10=0.04951 ratio=1.290 material=False

## T-C pre-registration (verbatim, written before any boosting fit)

Written **before any T-C boosting fit**, 2026-09-18. No threshold in this block is re-tuned after numbers exist.

On GTEx cultured fibroblasts only (the Stage 1 cohort and preprocessing, unchanged: `SMAFRZE==RNASEQ`, `SMTSD == "Cells - Cultured fibroblasts"`, gene filter ≥6 counts in ≥20% of samples, TMM/edgeR log2-CPM `prior.count=2`, per-gene z-score on training rows only). D1 n is too small; excluded from B1↔C1 transfer and from within-site CV, stated, not silently dropped.

Fit `sklearn.ensemble.HistGradientBoostingRegressor`. Hyperparameters by inner grouped CV on training folds only (outer folds unchanged: within-site `SMNABTCH`-grouped; transfer is B1→C1 and C1→B1). Inner-CV selection metric is Spearman ρ on the inner held-out batches; the selected hyperparameter tuple is frozen for that train split before any permutation of that split.

Report within-site ρ and **B1→C1 / C1→B1 transfer ρ** with the same donor-level permutation null (n_perm=200, seed `20260914`) and donor-bootstrap 95% CI (B=200, seed `20260918`), side by side with the ridge numbers from `results/fibro/stage1_cv.csv` and `results/fibro/stage1_transfer.csv`. Uncalibrated R² reported, never a gate.

This is a **diagnostic, not a replacement**. A nonlinear model yields a score, not a direction; the project's plane, angles, and "age movement per unit identity loss" criterion all require a direction, so a nonlinear model cannot be substituted into Stage 2 without discarding the geometry. Do not rescore GSE297234 with boosting.

**Pre-registered reading:**
- boosting transfer ρ within the ridge bootstrap CI (both directions) → linearity is not the limiting factor; ridge stays.
- boosting transfer ρ materially above the ridge CI (upper bound exceeded in both directions) → the linear direction is leaving signal on the table. Report the gap. Do not rerun Stage 2 with boosting. A direction-preserving nonlinear extension (e.g. a kernel or autoencoder latent with a linear age axis inside it) would be the next design question.
- boosting transfer ρ below ridge → note it and move on.


## T-C — ridge vs HistGradientBoosting

| cell | site | ridge_rho | ridge_rho_null | ridge_p | ridge_ci_lo | ridge_ci_hi | ridge_cal_r2 | ridge_r2 | hgb_rho | hgb_rho_null | hgb_p | hgb_ci_lo | hgb_ci_hi | hgb_cal_r2 | hgb_r2 | n |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| within_site | B1 | +0.549 | -0.010 | +0.005 | +0.481 | +0.604 | +0.303 | +0.302 | +0.442 | +0.009 | +0.005 | +0.360 | +0.510 | +0.220 | +0.218 | 421 |
| within_site | C1 | +0.355 | -0.021 | +0.005 | +0.211 | +0.464 | +0.113 | +0.101 | +0.184 | +0.004 | +0.005 | +0.067 | +0.312 | +0.037 | -0.001 | 224 |
| transfer | B1→C1 | +0.529 | -0.004 | +0.005 | +0.418 | +0.612 | +0.353 | +0.289 | +0.452 | -0.009 | +0.005 | +0.330 | +0.562 | +0.267 | +0.182 | 224 |
| transfer | C1→B1 | +0.561 | +0.013 | +0.005 | +0.500 | +0.617 | +0.312 | +0.251 | +0.456 | +0.003 | +0.005 | +0.359 | +0.516 | +0.229 | +0.139 | 421 |

## T-C — reading

Fired key: `boosting_below_ridge`. linearity_limiting=False.

boosting transfer ρ below ridge → note it and move on.

HGB null: donor-level shuffle of ages; fitted scores held fixed (n_perm=200, seed `20260914`). Ridge Stage 1 p-values on disk used refit-on-permuted-y. T-C reading compares HGB transfer ρ to the ridge bootstrap CI, not to ridge p-values. D1 n=7 skipped (n<10), stated.

Per-direction comparison (from disk):

- B1→C1 hgb ρ=+0.452 ridge ρ=+0.529 ridge CI=[+0.418, +0.612] inside=True above=False
- C1→B1 hgb ρ=+0.456 ridge ρ=+0.561 ridge CI=[+0.500, +0.617] inside=False above=False

## Direction vs score

A nonlinear model yields a score, not a direction; the project's plane, angles, and 'age movement per unit identity loss' criterion all require a direction, so a nonlinear model cannot be substituted into Stage 2 without discarding the geometry.

## Extrapolation (report-only, no gate)

Mahalanobis is in the leading PCA of the GTEx fibroblast training z-space. Nearest neighbour is Euclidean in that z-space. Missing overlap genes at z=0.

GSE297234 all-cell rows:

| cohort | variant | cell_line | day | n_cells | age_score | nn_euclidean | mahalanobis_pca | pca_k |
|---|---|---|---|---|---|---|---|---|
| GSE297234 | all_cell | GM00731 | +0.000 | +5021.000 | +2.769 | +455.630 | +25.015 | 50 |
| GSE297234 | all_cell | GM00731 | +3.000 | +6318.000 | -0.188 | +470.580 | +25.359 | 50 |
| GSE297234 | all_cell | GM00731 | +7.000 | +6915.000 | -1.416 | +483.361 | +25.067 | 50 |
| GSE297234 | all_cell | GM00731 | +10.000 | +4399.000 | +6.120 | +475.307 | +24.370 | 50 |
| GSE297234 | all_cell | GM23815 | +0.000 | +7782.000 | -2.423 | +452.379 | +23.997 | 50 |
| GSE297234 | all_cell | GM23815 | +3.000 | +6738.000 | -4.992 | +461.629 | +24.212 | 50 |
| GSE297234 | all_cell | GM23815 | +7.000 | +12209.000 | -9.740 | +506.756 | +26.324 | 50 |
| GSE297234 | all_cell | GM23815 | +10.000 | +7935.000 | -4.885 | +480.484 | +25.016 | 50 |

Cohort summaries:

| cohort | variant | n | nn_euclidean_min | nn_euclidean_median | nn_euclidean_max | mahalanobis_pca_min | mahalanobis_pca_median | mahalanobis_pca_max | pca_k |
|---|---|---|---|---|---|---|---|---|---|
| GSE226189_near_miss_bulk | external | 82 | +469.619 | +532.814 | +835.961 | +21.206 | +29.972 | +47.828 | 50 |
| GSE297234 | all_cell | 8 | +452.379 | +472.943 | +506.756 | +23.997 | +25.015 | +26.324 | 50 |
| GSE297234 | cluster | 24 | +433.289 | +473.666 | +792.533 | +21.919 | +24.659 | +34.970 | 50 |
| cxg_a19d1667 | external | 93 | +493.301 | +615.800 | +1651.660 | +24.580 | +30.544 | +99.708 | 50 |

GSE297234 cluster rows:

| cell_line | day | cluster | n_cells | age_score | nn_euclidean | mahalanobis_pca |
|---|---|---|---|---|---|---|
| GM00731 | +0.000 | +0.000 | +73.000 | -4.814 | +792.533 | +34.970 |
| GM00731 | +0.000 | +1.000 | +2535.000 | +5.464 | +449.088 | +25.343 |
| GM00731 | +0.000 | +2.000 | +2413.000 | +2.584 | +453.418 | +23.927 |
| GM00731 | +3.000 | +0.000 | +2159.000 | +1.497 | +495.800 | +26.509 |
| GM00731 | +3.000 | +1.000 | +915.000 | -0.162 | +560.790 | +28.221 |
| GM00731 | +3.000 | +2.000 | +3244.000 | +0.811 | +444.297 | +22.622 |
| GM00731 | +7.000 | +0.000 | +3054.000 | +4.406 | +453.908 | +24.154 |
| GM00731 | +7.000 | +1.000 | +3011.000 | -3.795 | +502.599 | +25.708 |
| GM00731 | +7.000 | +2.000 | +850.000 | -5.240 | +589.803 | +27.863 |
| GM00731 | +10.000 | +0.000 | +2645.000 | +8.757 | +471.954 | +24.391 |
| GM00731 | +10.000 | +1.000 | +708.000 | +2.248 | +580.964 | +27.566 |
| GM00731 | +10.000 | +2.000 | +1046.000 | +4.973 | +459.087 | +24.604 |
| GM23815 | +0.000 | +0.000 | +2896.000 | -0.514 | +455.640 | +23.727 |
| GM23815 | +0.000 | +1.000 | +3315.000 | +0.026 | +440.293 | +24.128 |
| GM23815 | +0.000 | +2.000 | +1571.000 | -1.839 | +445.276 | +23.380 |
| GM23815 | +3.000 | +0.000 | +2637.000 | -2.268 | +433.289 | +21.919 |
| GM23815 | +3.000 | +1.000 | +3346.000 | -5.127 | +474.361 | +23.898 |
| GM23815 | +3.000 | +2.000 | +755.000 | -2.272 | +552.144 | +27.622 |
| GM23815 | +7.000 | +0.000 | +4407.000 | -10.151 | +531.327 | +27.212 |
| GM23815 | +7.000 | +1.000 | +6163.000 | -5.693 | +472.514 | +24.467 |
| GM23815 | +7.000 | +2.000 | +1639.000 | -13.728 | +635.755 | +29.566 |
| GM23815 | +10.000 | +0.000 | +2910.000 | +0.965 | +472.971 | +24.714 |
| GM23815 | +10.000 | +1.000 | +793.000 | -8.148 | +603.696 | +27.961 |
| GM23815 | +10.000 | +2.000 | +4232.000 | -4.218 | +478.730 | +24.575 |

## unusable

- **GSE113957:** no usable count matrix for GSE113957: last error [counts] GSE113957_fpkm.txt.gz: values are not integer counts (cannot apply TMM/log2-CPM). min=0.0 max=249213345.0 mean=1009452.0862567581
- **GSE307377:** no usable count matrix for GSE307377: last error [counts] GSE307377_series_matrix.txt.gz: zero-size array to reduction operation fmin which has no identity

GSE113957 record (from SOFT + files this session): n_gsm=143, n_with_stated_age=135, n_missing_age=8, n_stated_age_<18=24, n_adult≥18=111. Ages from Sample_title (e.g. years in the title). Series_type is high-throughput sequencing (GPL16791, GPL18573). Supplementary `GSE113957_fpkm.txt.gz` is FPKM, not integer counts (min=0 max=249213345 mean=1009452.0862567581). Series matrix HEAD 404. TMM/log2-CPM not applied. Primary external cell not scored.

GSE307377 record (from SOFT + files this session): n_gsm=9, ages stated in `Sample_characteristics_ch1` key `age` (23–72 y). `GSE307377_raw.txt.gz` is HOMER `analyzeRepeats.pl` output: gene IDs are `NM_`, columns include `start`/`end`/`Length`/`Copies` and `tags/*` with half-integer values. Not integer counts; TMM/log2-CPM not applied. `GSE307377_tpm.txt.gz` is TPM. Series matrix table has 0 genes. n=9 is the SCOPING_TISSUE `near_miss` row. SOFT `cell line: Primary` is not a donor id; ages are per GSM. Not scored.

## unverified

None recorded.

## What this changes in FINDINGS_FIBRO.md (quote, not edited)

Stage 2 verdict sentence, quoted:

> age score does not decline → the ruler reads a static donor property, not a modifiable state. Step 3 is not supported by this data.

This file does not edit `FINDINGS_FIBRO.md`. T-A/T-B/T-C do not replace that sentence. They say whether it can be read as a statement about reprogramming.

T-A: Stage 2 `no_decline` must not be reported as a finding about reprogramming (`primary_unscored`).

T-B: d0→d7 gives p ≤ 0.05 with no technical difference at d10 → the discrepancy is unexplained.

T-C: boosting transfer ρ below ridge → note it and move on.

## Limitations

1. Frozen ruler trained on GTEx V10 cultured fibroblasts, public `AGE` 10-year bins; ρ is against a coarse ordinal on GTEx and against stated years on external records.
2. Cross-platform: GTEx bulk polyA RNASeQCv2.4.2 vs 10x 3' scRNA-seq pseudobulk (GSE297234) vs whatever T-A resolved.
3. GSE297234 has two donors. T-B cannot estimate donor-to-donor variance of a technical effect.
4. Cluster labels are independent per donor × timepoint; the secondary null pairs by cell-count rank, not by biological identity of a cluster.
5. Missing overlap genes at z=0. Nothing fitted on GSE297234 or on T-A matrices.
6. HistGradientBoosting is a score, not a direction; T-C is not a Stage 2 replacement.
7. Seed `20260914` (permutation / random directions / HGB). Donor-bootstrap seed `20260918`. n_perm=200, n_boot=200, n_random=200.
8. GSE325735 not opened.
9. Extrapolation distances have no pass/fail.
10. A record that does not state per-sample donor age is unusable; ages are not inferred.
11. GSE113957 supplementary is FPKM (transcript-level IDs), not integer counts; series matrix 404. Primary cell unscored.
12. HGB permutation null holds scores fixed. Ridge Stage 1 null refit on permuted y.
13. CXG atlas scored on `development_stage` (stated year-old stages); `Age` is also in obs and was not the column used.
14. GSE226189 is bulk (GPL24676), not a substitute for the primary cell and not a 10x platform check.

## Files

| path | content |
|---|---|
| `src/fibro2_common.py`, `fibro2_ta.py`, `fibro2_tb.py`, `fibro2_tc.py`, `fibro2_extrap.py`, `fibro2_findings.py`, `fibro2_run.py` | code |
| `results/fibro2/` | manifest, tables, logs, frozen-projection scores |
| `FINDINGS_FIBRO2.md` | this file |
| `PROGRESS_FIBRO2.md` | resume state |


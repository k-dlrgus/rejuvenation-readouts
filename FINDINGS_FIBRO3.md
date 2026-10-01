# FINDINGS_FIBRO3 — what happens between day 7 and day 10?

**Status:** Task 4: superseded primary passed (`instrument_validated`); Stage 2 trajectory interpretable as a statement about the cells. Task 1: `cell_intrinsic`. Task 2: nn_euclidean: age score does not rise with distance ρ=-0.625 CI=[-0.834, -0.320] n=22. Not a gate. Task 3: n_ok=8 / n_rows=12 (d7 not the Stage 2 endpoint). Seed `20260914`. boot `20260918`. n_perm=200. n_boot=200. n_random=200. Frozen ruler `frozen_ruler_ridge_raw.npz` exists=True.

Does not modify any existing `FINDINGS_*.md` or `FALSIFICATION.md`. No GSE325735. Uncalibrated R² is reported and is never a gate. Nothing averaged across donors, clusters, or regimes. Refit nothing. d7 is not the primary Stage 2 endpoint.

Reproduced by `src/fibro3_run.py`. Frozen ruler: `<repo>\results\fibro\frozen_ruler_ridge_raw.npz` exists=True.

Flag: `results/fibro3/PREREG_TASK4.flag` exists=True. `PREREG_TASK1.flag` exists=True. `CLUSTER_SETTINGS.flag` exists=True. `DECLARED_BEFORE_SCORES.flag` exists=True.

## Standing constraint

d7 was not the pre-registered Stage 2 endpoint. Nothing in this prompt may promote d7 to the primary result. The pre-registered d0→d10 verdict stands as written in `FINDINGS_FIBRO.md`. This prompt characterizes the trajectory; it does not re-decide it.

## Task 4 supersession (verbatim, written before any Task 4 ρ)

Written **before any Task 4 ρ**, 2026-09-18. No threshold in this block is re-tuned after numbers exist. Option (a) is not run.

Choice: **(b)**. The primary external cell is superseded from GSE113957 to CELLxGENE dataset `a19d1667-a7b5-4556-9e5f-f9bfa690c0f1` (Human dermal fibroblast atlas; 10x-restricted adult donors).

Ground: FINDINGS_FIBRO2.md T-A pre-registration already stated that a 10x cohort is worth more than a larger bulk one because the platform shift is the thing under test. Stage 2 projects a GTEx bulk polyA (RNASeQCv2.4.2) frozen ruler onto 10x 3' scRNA-seq pseudobulks. GSE113957 is bulk FPKM; an FPKM rank-transform would be a different transform from the GTEx TMM/log2-CPM spec and would remain bulk-to-bulk, so it would not test the platform shift. The CELLxGENE 10x cohort tests that shift directly.

Take the Stage 1 ridge direction exactly as frozen in `results/fibro/frozen_ruler_ridge_raw.npz`. Refit nothing. Recalibrate nothing. Rescale nothing. Score Spearman ρ vs donor age on the superseded primary.

Freeze-and-project path (FINDINGS_EXTERNAL.md / T-A): TMM/edgeR log2-CPM (`prior.count=2`); genes restricted to the overlap with the frozen ruler; missing overlap genes left at z=0; z uses the frozen train mean and sd. Adult restriction ≥18. Ages from the record column actually used. Do not infer a donor age that the record does not state.

Metric: Spearman ρ vs stated age. Donor-level permutation null n_perm=200 (scores held fixed; shuffle donor ages; fresh Generator on seed `20260914`). Donor-bootstrap 95% percentile CI, B=200, seed `20260918`. Calibrated R² and uncalibrated R² reported; uncalibrated R² is never a gate. Pearson r reported, not a gate.

**Pre-registered reading:** if the superseded primary passes (ρ > 0, null ≤ 0.05), the frozen fibroblast ruler is a validated instrument outside GTEx and the Stage 2 trajectory is interpretable as a statement about the cells. If it fails, Stage 2 remains uninterpretable and every number in Tasks 1–3 is descriptive only.

GSE113957 is not scored in this prompt. Option (a) is not run. GSE226189 is not the superseded primary. GSE325735 is out of scope. Do not report both (a) and (b) and pick the better number.


## Task 4 — result

choice=`b` option_a_run=False accession=`a19d1667-a7b5-4556-9e5f-f9bfa690c0f1` role=`superseded_primary_external` platform_detail=`10x (restricted)` used_10x_only=True.

n_samples=93 n_donors=93 n_excluded_child=0 n_overlap=17047 n_missing_z0=6580 age_source=['h5ad_obs:development_stage'] age range 18.0–79.0 n_unique_age=39.

ρ=+0.455 null=-0.003 p=+0.005 CI=[+0.253, +0.605] calR²=+0.171 R²=-9.664 r=+0.414 pass_rho=True n_perm=200 n_boot=200.

Fired key: `instrument_validated`. interpretable=True.

the frozen fibroblast ruler is a validated instrument outside GTEx and the Stage 2 trajectory is interpretable as a statement about the cells.

| accession | role | choice | tag | n_samples | n_donors | n_overlap | n_missing_z0 | rho | rho_null | rho_p | rho_ci_lo | rho_ci_hi | cal_r2 | r2 | r | pass_rho | used_10x_only | platform_detail |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| a19d1667-a7b5-4556-9e5f-f9bfa690c0f1 | superseded_primary_external | b | cxg_a19d1667 | 93 | 93 | 17047 | 6580 | +0.455 | -0.003 | +0.005 | +0.253 | +0.605 | +0.171 | -9.664 | +0.414 | True | True | 10x (restricted) |

## Task 1 pre-registration (verbatim, written before any joint-cluster age score)

Written **before any joint-cluster age score**, 2026-09-18. Not swept after numbers exist.

JOINT CLUSTERING (code-listed; chosen before looking at any age score):
- Scope: all cells of one donor, all timepoints concatenated. Donors clustered separately. Never pooled across donors.
- Genes: frozen-ruler overlap genes with count-sum > 0 (same gene filter as Stage 2 cluster spec in FINDINGS_FIBRO.md).
- Transform: log1p(CP10k) on those genes. Library size = row sum of the kept overlap genes (Stage 2).
- Dimensionality: sklearn.decomposition.PCA(n_components=20, or n_cells−1 if smaller, svd_solver="randomized", random_state=20260914) on float32 log1p(CP10k). Cell order: libraries concatenated in day order 0, 3, 7, 10.
- Cluster: sklearn.cluster.KMeans(n_clusters=3, random_state=20260914, n_init=10).
- k=3 is the Stage 2 frozen k, applied jointly so cluster IDs are comparable across days of the same donor. Resolution is k (not a Leiden resolution). Not swept.
- Seed `20260914`.
- MIN_CELLS = 20 (prompt minimum) to score age / pluri / extrapolation / random-direction null at a cluster×timepoint. QC (n, fraction, median UMI, median genes, mito fraction) is reported at all n, including n<20.

Amendment, still before any age score: the first listed dimensionality was IncrementalPCA(n_components=20, batch_size=2048). A run reached IncrementalPCA partial_fit on 2048×20354 and did not finish in usable time (~3 min/batch). No cluster labels and no age scores were written. Dimensionality was replaced by randomized PCA on the same matrix, same n_components, same KMeans k. This is not a k/resolution sweep.


Written **before any joint-cluster age score**, 2026-09-18. No threshold in this block is re-tuned after numbers exist.

Using the GSE297234 Cell Ranger `filtered_feature_bc_matrix.h5` matrices already on disk. Cluster cells jointly across all timepoints within each donor (not within timepoint). Clustering settings are in `results/fibro3/CLUSTER_SETTINGS.flag` and are listed before any age score.

Per donor × cluster × timepoint report: n cells, fraction of that timepoint's cells, median UMI, median genes, mito fraction, frozen age score, pluripotency−fibroblast score (frozen lists from FINDINGS_FIBRO.md, with and without OSKM-family genes). Frozen ruler `results/fibro/frozen_ruler_ridge_raw.npz` used as-is. Refit nothing. Missing overlap genes left at z=0.

MIN_CELLS = 20. Do not read a d7–d10 trend off a cluster×timepoint with n < 20.

Operational definitions (frozen before any joint-cluster age score):
- C_ok(donor) = clusters with n ≥ 20 at both d7 and d10.
- Cluster reverse = age_score(d10) > age_score(d7) (score increased; older on the ruler).
- Cluster stable_or_decline = age_score(d10) ≤ age_score(d7).
- High-age cluster at d7 = cluster with n ≥ 20 at d7 and age_score(d7) equal to the maximum among those.
- Low-age cluster at d7 = cluster with n ≥ 20 at d7 and age_score(d7) equal to the minimum among those.
- Fraction = n_cells(cluster, day) / n_cells(donor, day). A missing cluster has fraction 0.
- Mix shift (reversal-direction) = frac(high-age, d10) > frac(high-age, d7) OR frac(low-age, d10) < frac(low-age, d7).
- Primary reading is the aged donor GM00731. GM23815 is tabulated as a contrast, never pooled, never the reading.

**Pre-registered reading** (only the outcome that fired; d7 is not the Stage 2 endpoint):
- **Compositional.** Individual clusters' age scores are stable or keep declining from d7 to d10, but the *mix* shifts — a high-age-score cluster grows, or a low-scoring one shrinks or drops out. Then the all-cell reversal is a composition artifact, the same failure mode as the blood result in `FINDINGS_BRAIN_PHASE1.md` §2. Say so in those words.
- **Cell-intrinsic.** The same cells' clusters individually reverse from d7 to d10. Then reprogramming genuinely pushes cells back along the age axis after d7, and that is a biological claim worth its own test.
- **Mixed.** Report which clusters do which. Do not average them into one number.
- **Unresolvable.** Too few cells per cluster × timepoint (threshold 20 cells) to read either way. Say so and stop rather than reading a trend off thin cells.

Firing rule (aged donor only):
- Unresolvable if C_ok is empty.
- Cell-intrinsic if C_ok is non-empty AND every cluster in C_ok reverses.
- Compositional if C_ok is non-empty AND every cluster in C_ok is stable_or_decline AND mix shift is True.
- Mixed if C_ok is non-empty AND at least one cluster reverses AND at least one does not.
- If all C_ok are stable_or_decline and mix shift is False: key=`no_listed_bullet`; report the tables; do not invent a fifth story.

Clusters are not lineage-tracked. "The same cells" is an inference from joint cluster identity, not from barcodes.


## Task 1 — clustering

method=PCA svd_solver=randomized k=3 n_components=20 random_state=`20260914` n_init=10. MIN_CELLS=20.

- GM00731: n_cells=22653 n_genes_kept=20354 npc=20 k=3 sizes=[12039, 2570, 8044] days=[0, 3, 7, 10] mito_rule=`symbol upper startswith MT-` n_mito_genes=13
- GM23815: n_cells=34664 n_genes_kept=20350 npc=20 k=3 sizes=[16115, 3085, 15464] days=[0, 3, 7, 10] mito_rule=`symbol upper startswith MT-` n_mito_genes=13

## Task 1 — per donor × cluster × timepoint

| cell_line | day | cluster | gsm | age_years | n_cells | frac_timepoint | median_UMI | median_genes | mito_frac_median_cell | below_min_cells | age_score | pluri_with_OSKM | pluri_without_OSKM | pluri_primary |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | 0 | 0 | GSM8986586 | 96 | 4941 | +0.984 | +14472.000 | +4096.000 | +0.025 | False | +4.224 | -1.391 | -2.546 | -1.391 |
| GM00731 | 0 | 1 | GSM8986586 | 96 | 76 | +0.015 | +3000.500 | +1438.500 | +0.064 | False | -3.165 | +0.566 | +1.328 | +0.566 |
| GM00731 | 0 | 2 | GSM8986586 | 96 | 4 | +0.001 | +9300.500 | +3400.500 | +0.037 | True | NA | NA | NA | NA |
| GM00731 | 3 | 0 | GSM8986587 | 96 | 948 | +0.150 | +7593.000 | +3190.000 | +0.028 | False | +5.541 | +5.360 | -0.081 | +5.360 |
| GM00731 | 3 | 1 | GSM8986587 | 96 | 946 | +0.150 | +2857.500 | +671.500 | +0.030 | False | -0.053 | +10.803 | +2.960 | +10.803 |
| GM00731 | 3 | 2 | GSM8986587 | 96 | 4424 | +0.700 | +12613.500 | +4110.500 | +0.038 | False | +0.639 | +8.503 | +2.141 | +8.503 |
| GM00731 | 7 | 0 | GSM8986588 | 96 | 2719 | +0.393 | +10222.000 | +3508.000 | +0.034 | False | +4.806 | +4.665 | +1.466 | +4.665 |
| GM00731 | 7 | 1 | GSM8986588 | 96 | 850 | +0.123 | +2776.000 | +1550.000 | +0.126 | False | -5.380 | +10.059 | +7.753 | +10.059 |
| GM00731 | 7 | 2 | GSM8986588 | 96 | 3346 | +0.484 | +13085.000 | +4188.500 | +0.026 | False | -2.904 | +9.832 | +6.377 | +9.832 |
| GM00731 | 10 | 0 | GSM8986589 | 96 | 3431 | +0.780 | +14192.000 | +4309.000 | +0.032 | False | +8.595 | +2.012 | -0.948 | +2.012 |
| GM00731 | 10 | 1 | GSM8986589 | 96 | 698 | +0.159 | +4086.500 | +1997.500 | +0.131 | False | +2.249 | +5.950 | +4.680 | +5.950 |
| GM00731 | 10 | 2 | GSM8986589 | 96 | 270 | +0.061 | +22356.000 | +5510.500 | +0.017 | False | +0.619 | +8.088 | +4.929 | +8.088 |
| GM23815 | 0 | 1 | GSM8986590 | 22 | 68 | +0.009 | +1587.500 | +1015.000 | +0.091 | False | -6.688 | +0.570 | +1.047 | +0.570 |
| GM23815 | 0 | 2 | GSM8986590 | 22 | 7714 | +0.991 | +8169.500 | +3039.000 | +0.029 | False | -0.524 | -0.960 | -1.051 | -0.960 |
| GM23815 | 3 | 0 | GSM8986591 | 22 | 5458 | +0.810 | +11909.500 | +4151.000 | +0.049 | False | -3.878 | +9.202 | +4.484 | +9.202 |
| GM23815 | 3 | 1 | GSM8986591 | 22 | 594 | +0.088 | +2454.500 | +835.500 | +0.047 | False | -5.038 | +12.670 | +5.454 | +12.670 |
| GM23815 | 3 | 2 | GSM8986591 | 22 | 686 | +0.102 | +7542.000 | +3190.000 | +0.039 | False | +0.224 | +5.794 | +0.506 | +5.794 |
| GM23815 | 7 | 0 | GSM8986592 | 22 | 7129 | +0.584 | +7860.000 | +3068.000 | +0.031 | False | -9.197 | +11.560 | +8.686 | +11.560 |
| GM23815 | 7 | 1 | GSM8986592 | 22 | 1662 | +0.136 | +1272.000 | +804.000 | +0.142 | False | -13.366 | +13.220 | +11.581 | +13.220 |
| GM23815 | 7 | 2 | GSM8986592 | 22 | 3418 | +0.280 | +4929.000 | +2148.500 | +0.038 | False | -3.439 | +5.270 | +2.151 | +5.270 |
| GM23815 | 10 | 0 | GSM8986593 | 22 | 3528 | +0.445 | +12567.000 | +3959.000 | +0.046 | False | -4.728 | +9.171 | +6.102 | +9.171 |
| GM23815 | 10 | 1 | GSM8986593 | 22 | 761 | +0.096 | +2446.000 | +1306.000 | +0.200 | False | -8.412 | +9.523 | +8.176 | +9.523 |
| GM23815 | 10 | 2 | GSM8986593 | 22 | 3646 | +0.459 | +6578.000 | +2607.000 | +0.040 | False | +0.388 | +3.686 | +0.384 | +3.686 |

## Task 1 — composition d7 vs d10

| cell_line | cluster | n_d7 | frac_d7 | n_d10 | frac_d10 | delta_frac | status | age_score_d7 | age_score_d10 | below_min_d7 | below_min_d10 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | 0 | 2719 | +0.393 | 3431 | +0.780 | +0.387 | grow | +4.806 | +8.595 | False | False |
| GM00731 | 1 | 850 | +0.123 | 698 | +0.159 | +0.036 | grow | -5.380 | +2.249 | False | False |
| GM00731 | 2 | 3346 | +0.484 | 270 | +0.061 | -0.422 | shrink | -2.904 | +0.619 | False | False |
| GM23815 | 0 | 7129 | +0.584 | 3528 | +0.445 | -0.139 | shrink | -9.197 | -4.728 | False | False |
| GM23815 | 1 | 1662 | +0.136 | 761 | +0.096 | -0.040 | shrink | -13.366 | -8.412 | False | False |
| GM23815 | 2 | 3418 | +0.280 | 3646 | +0.459 | +0.180 | grow | -3.439 | +0.388 | False | False |

## Task 1 — reading

Fired key: `cell_intrinsic`. donor=GM00731 MIN_CELLS=20 C_ok=[0, 1, 2] mix_shift=True.

Cell-intrinsic. The same cells' clusters individually reverse from d7 to d10. Then reprogramming genuinely pushes cells back along the age axis after d7, and that is a biological claim worth its own test. clusters=[0, 1, 2].

C_ok cluster d7→d10 (aged donor; not averaged):

| cluster | n_d7 | n_d10 | age_score_d7 | age_score_d10 | delta_age_d10_minus_d7 | reverse | stable_or_decline |
|---|---|---|---|---|---|---|---|
| 0 | 2719 | 3431 | +4.806 | +8.595 | +3.789 | True | False |
| 1 | 850 | 698 | -5.380 | +2.249 | +7.629 | True | False |
| 2 | 3346 | 270 | -2.904 | +0.619 | +3.523 | True | False |

mix_notes: high-age cluster 0 grows frac 0.3932→0.7799

d7 is not the pre-registered Stage 2 endpoint. This reading characterizes the d7–d10 trajectory. It does not re-decide FINDINGS_FIBRO.md.

## Task 2 pre-registration (verbatim, written before any distance)

Written **before any Task 2 distance or correlation**, 2026-09-18.

Report-only, no gate. Per donor × cluster × timepoint with n ≥ 20, report extrapolation distance to the GTEx fibroblast training distribution: PCA k=50 nearest-neighbour Euclidean and Mahalanobis, same measures as FINDINGS_FIBRO2.md. Missing overlap genes at z=0. Frozen ruler used as-is.

Then Spearman ρ, across all donor × cluster × timepoint rows with n ≥ 20, between each distance and frozen age score. Donor-level permutation is not the unit here (the unit is a cluster×timepoint pseudobulk). Permutation null: shuffle age scores, n_perm=200, seed `20260914`. Bootstrap 95% percentile CI on the rows, B=200, seed `20260918`. Per-donor ρ is reported and is not averaged. Uncalibrated R² reported, never a gate.

If age score rises with distance from the training distribution, the ruler is being read outside its calibrated region and the d10 value is an extrapolation artifact rather than a measurement. State this as an observation with its correlation and CI. It is not a gate and does not overturn any verdict.


## Task 2 — extrapolation (report-only, no gate)

| cell_line | day | cluster | n_cells | below_min_cells | age_score | nn_euclidean | mahalanobis_pca | pca_k |
|---|---|---|---|---|---|---|---|---|
| GM00731 | 0 | 0 | 4941 | False | +4.224 | +453.970 | +24.150 | 50 |
| GM00731 | 0 | 1 | 76 | False | -3.165 | +764.998 | +32.640 | 50 |
| GM00731 | 0 | 2 | 4 | True | NA | +1466.888 | +74.549 | 50 |
| GM00731 | 3 | 0 | 948 | False | +5.541 | +441.622 | +21.887 | 50 |
| GM00731 | 3 | 1 | 946 | False | -0.053 | +554.940 | +27.448 | 50 |
| GM00731 | 3 | 2 | 4424 | False | +0.639 | +469.682 | +24.433 | 50 |
| GM00731 | 7 | 0 | 2719 | False | +4.806 | +458.703 | +23.528 | 50 |
| GM00731 | 7 | 1 | 850 | False | -5.380 | +589.008 | +27.625 | 50 |
| GM00731 | 7 | 2 | 3346 | False | -2.904 | +496.744 | +24.897 | 50 |
| GM00731 | 10 | 0 | 3431 | False | +8.595 | +470.636 | +23.915 | 50 |
| GM00731 | 10 | 1 | 698 | False | +2.249 | +577.495 | +25.575 | 50 |
| GM00731 | 10 | 2 | 270 | False | +0.619 | +481.955 | +24.563 | 50 |
| GM23815 | 0 | 1 | 68 | False | -6.688 | +938.376 | +39.981 | 50 |
| GM23815 | 0 | 2 | 7714 | False | -0.524 | +446.635 | +23.046 | 50 |
| GM23815 | 3 | 0 | 5458 | False | -3.878 | +459.140 | +22.914 | 50 |
| GM23815 | 3 | 1 | 594 | False | -5.038 | +574.634 | +28.069 | 50 |
| GM23815 | 3 | 2 | 686 | False | +0.224 | +441.992 | +21.503 | 50 |
| GM23815 | 7 | 0 | 7129 | False | -9.197 | +513.218 | +25.741 | 50 |
| GM23815 | 7 | 1 | 1662 | False | -13.366 | +628.240 | +28.561 | 50 |
| GM23815 | 7 | 2 | 3418 | False | -3.439 | +466.403 | +23.506 | 50 |
| GM23815 | 10 | 0 | 3528 | False | -4.728 | +484.735 | +24.340 | 50 |
| GM23815 | 10 | 1 | 761 | False | -8.412 | +606.858 | +27.117 | 50 |
| GM23815 | 10 | 2 | 3646 | False | +0.388 | +470.006 | +23.938 | 50 |

## Task 2 — correlation

| scope | cell_line | distance | n | rho | rho_null | rho_p | rho_ci_lo | rho_ci_hi | r | r2 | cal_r2 | age_rises_with_distance |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all_donor_cluster_timepoint_n_ge_20 | NA | nn_euclidean | 22 | -0.625 | +0.014 | +1.000 | -0.834 | -0.320 | -0.482 | -11479.969 | +0.232 | False |
| per_donor_GM00731 | GM00731 | nn_euclidean | 11 | -0.782 | -0.036 | +0.990 | -0.972 | -0.486 | -0.622 | -17404.947 | +0.387 | False |
| per_donor_GM23815 | GM23815 | nn_euclidean | 11 | -0.809 | -0.027 | +1.000 | -0.981 | -0.424 | -0.494 | -19661.603 | +0.244 | False |
| all_donor_cluster_timepoint_n_ge_20 | NA | mahalanobis_pca | 22 | -0.572 | +0.019 | +0.990 | -0.772 | -0.253 | -0.473 | -30.289 | +0.223 | False |
| per_donor_GM00731 | GM00731 | mahalanobis_pca | 11 | -0.891 | -0.018 | +1.000 | -1.000 | -0.654 | -0.679 | -37.425 | +0.462 | False |
| per_donor_GM23815 | GM23815 | mahalanobis_pca | 11 | -0.773 | -0.018 | +1.000 | -0.944 | -0.360 | -0.477 | -61.383 | +0.227 | False |

- nn_euclidean: age score does not rise with distance ρ=-0.625 CI=[-0.834, -0.320] n=22. Not a gate.
- mahalanobis_pca: age score does not rise with distance ρ=-0.572 CI=[-0.772, -0.253] n=22. Not a gate.

## Task 3 pre-registration (verbatim, written before any per-cluster p)

Written **before any per-cluster random-direction p**, 2026-09-18.

For each donor × cluster with ≥20 cells at both endpoints, rerun the pre-registered random-direction null (200 permuted-weight directions, seed `20260914`) on d0→d7 and d0→d10 separately. Endpoint decline = age_score(day 0) − age_score(day T). Positive = younger on the ruler. Empirical p = (n_random ≥ real + 1) / (n_random + 1). Joint cluster IDs pair the same cluster across days (not size-rank). Do not average clusters or donors. d7 is not the Stage 2 endpoint.


## Task 3 — per-cluster random-direction null

| cell_line | cluster | endpoint | primary | n_cells_d0 | n_cells_end | score_d0 | score_end | decline | p | n_random_ge_real | n_random | pass_p | ok | reason |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | 0 | d0→d7 | True | +4941.000 | 2719 | +4.224 | +4.806 | -0.582 | +0.652 | +130.000 | 200 | False | True | NA |
| GM00731 | 0 | d0→d10 | True | +4941.000 | 3431 | +4.224 | +8.595 | -4.371 | +1.000 | +200.000 | 200 | False | True | NA |
| GM00731 | 1 | d0→d7 | True | +76.000 | 850 | -3.165 | -5.380 | +2.215 | +0.144 | +28.000 | 200 | False | True | NA |
| GM00731 | 1 | d0→d10 | True | +76.000 | 698 | -3.165 | +2.249 | -5.414 | +0.881 | +176.000 | 200 | False | True | NA |
| GM00731 | 2 | d0→d7 | True | +4.000 | 3346 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint |
| GM00731 | 2 | d0→d10 | True | +4.000 | 270 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint |
| GM23815 | 0 | d0→d7 | False | NA | 7129 | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint |
| GM23815 | 0 | d0→d10 | False | NA | 3528 | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint |
| GM23815 | 1 | d0→d7 | False | +68.000 | 1662 | -6.688 | -13.366 | +6.678 | +0.020 | +3.000 | 200 | True | True | NA |
| GM23815 | 1 | d0→d10 | False | +68.000 | 761 | -6.688 | -8.412 | +1.724 | +0.174 | +34.000 | 200 | False | True | NA |
| GM23815 | 2 | d0→d7 | False | +7714.000 | 3418 | -0.524 | -3.439 | +2.915 | +0.005 | +0.000 | 200 | True | True | NA |
| GM23815 | 2 | d0→d10 | False | +7714.000 | 3646 | -0.524 | +0.388 | -0.911 | +0.746 | +149.000 | 200 | False | True | NA |

d7 is not the Stage 2 endpoint. Rows are not averaged across clusters or donors.

Aged GM00731 d0→d7 (the all-cell T-B cell had p=+0.005): cluster 0 decline=-0.582 p=+0.652 n_ge=130/200 n_d0=4941 n_d7=2719; cluster 1 decline=+2.215 p=+0.144 n_ge=28/200 n_d0=76 n_d7=850; cluster 2 skipped (n_cells<20 at an endpoint).

## What this changes in FINDINGS_FIBRO.md and FINDINGS_FIBRO2.md (quote, not edited)

This file does not edit `FINDINGS_FIBRO.md` or `FINDINGS_FIBRO2.md`.

FINDINGS_FIBRO.md Stage 2 verdict sentence, quoted:

> age score does not decline → the ruler reads a static donor property, not a modifiable state. Step 3 is not supported by this data.

That d0→d10 all-cell verdict stands. Nothing here promotes d7 to the primary Stage 2 result.

FINDINGS_FIBRO2.md T-A sentence, quoted:

> Stage 2 no_decline must not be reported as a finding about reprogramming until the primary external cell is scored and passes.

Task 4 option (b) supersedes the primary external cell to the CELLxGENE 10x cohort `a19d1667-a7b5-4556-9e5f-f9bfa690c0f1` and that cell passed (ρ > 0, null ≤ 0.05). Under the Task 4 reading, the frozen fibroblast ruler is a validated instrument outside GTEx and the Stage 2 trajectory is interpretable as a statement about the cells. The quoted T-A sentence is not edited in FINDINGS_FIBRO2.md.

FINDINGS_FIBRO2.md T-B sentence, quoted:

> d0→d7 gives p ≤ 0.05 with no technical difference at d10 → the discrepancy is unexplained.

Task 1 locates that discrepancy in the joint-cluster tables above. It does not replace the T-B sentence and does not make d7 the Stage 2 endpoint.

## Limitations

1. Two donors in GSE297234.
2. One timepoint spacing (days 0, 3, 7, 10); nothing is measured between d7 and d10.
3. Clustering choice is frozen (PCA svd_solver=randomized n_components=20, KMeans k=3, seed `20260914`) and was not swept. IncrementalPCA was listed first and did not complete; randomized PCA replaced it before any age score.
4. Clusters are not lineage-tracked, so "the same cells" is an inference from cluster identity, not from barcodes.
5. Frozen ruler trained on GTEx V10 cultured fibroblasts, public `AGE` 10-year bins.
6. Cross-platform: GTEx bulk polyA RNASeQCv2.4.2 vs 10x 3' scRNA-seq pseudobulk vs CELLxGENE 10x.
7. Missing overlap genes at z=0. Nothing fitted on GSE297234 or on the CELLxGENE matrix.
8. MIN_CELLS=20. Cluster×timepoint rows below that threshold are QC-only.
9. Seed `20260914` (permutation / random directions / k-means). Donor-bootstrap seed `20260918`. n_perm=200, n_boot=200, n_random=200.
10. GSE325735 not opened. Option (a) (GSE113957 FPKM rank-transform) was not run.
11. Task 2 distances have no pass/fail.
12. Uncalibrated R² is reported and is never a gate.

## Files

| path | content |
|---|---|
| `src/fibro3_common.py`, `fibro3_task4.py`, `fibro3_task1.py`, `fibro3_task2.py`, `fibro3_task3.py`, `fibro3_findings.py`, `fibro3_run.py` | code |
| `results/fibro3/` | manifest, flags, tables, logs |
| `FINDINGS_FIBRO3.md` | this file |
| `PROGRESS_FIBRO3.md` | resume state |

## Columns actually read

- **t4_scores_csv** source=`<repo>\results\fibro2\ta_cxg_a19d1667_scores.csv` columns=['sample', 'donor', 'age_years', 'age_source', 'n_cells', 'assay', 'age_score']
- **t4_proj_npz** source=`<repo>\data\processed\fibro2\ta_cxg_a19d1667_proj.npz` columns=['logcpm', 'pred', 'y', 'donor', 'missing', 'Z']
- **t4_ta_result_json** source=`<repo>\results\fibro2\ta_cxg_a19d1667_result.json` columns=['accession', 'age_max', 'age_min', 'boot_seed', 'cal_r2', 'n_boot', 'n_donors', 'n_ensembl', 'n_excluded_child', 'n_excluded_no_age', 'n_missing_z0', 'n_overlap', 'n_perm', 'n_ruler', 'n_samples', 'n_symbol_fallback', 'n_unique_age', 'pass_rho', 'platform_detail', 'platform_kind', 'r', 'r2', 'rho', 'rho_ci_hi', 'rho_ci_lo', 'rho_null', 'rho_p', 'role', 'seed', 'tag', 'used_10x_only']
- **t1_stage2_geo_samples** source=`<repo>\results\fibro\stage2_geo_samples.csv` columns=['gsm', 'title', 'source', 'organism', 'molecule', 'description', 'characteristics', 'growth_protocol', 'cell_line', 'cell_type', 'tissue', 'treatment', 'day', 'age_from_characteristics', 'age_from_growth_protocol', 'library_name', 'fields_read', 'age_years', 'age_source']
- **t1_h5_GSM8986586** source=`GSM8986586_GM00731_D0_filtered_feature_bc_matrix.h5` columns=['barcode', 'gene_id', 'symbol', 'feature_type', 'genome']
- **t1_h5_GSM8986587** source=`GSM8986587_GM00731_D3_filtered_feature_bc_matrix.h5` columns=['barcode', 'gene_id', 'symbol', 'feature_type', 'genome']
- **t1_h5_GSM8986588** source=`GSM8986588_GM00731_D7_filtered_feature_bc_matrix.h5` columns=['barcode', 'gene_id', 'symbol', 'feature_type', 'genome']
- **t1_h5_GSM8986589** source=`GSM8986589_GM00731_D10_filtered_feature_bc_matrix.h5` columns=['barcode', 'gene_id', 'symbol', 'feature_type', 'genome']
- **t1_h5_GSM8986590** source=`GSM8986590_GM23815_D0_filtered_feature_bc_matrix.h5` columns=['barcode', 'gene_id', 'symbol', 'feature_type', 'genome']
- **t1_h5_GSM8986591** source=`GSM8986591_GM23815_D3_filtered_feature_bc_matrix.h5` columns=['barcode', 'gene_id', 'symbol', 'feature_type', 'genome']
- **t1_h5_GSM8986592** source=`GSM8986592_GM23815_D7_filtered_feature_bc_matrix.h5` columns=['barcode', 'gene_id', 'symbol', 'feature_type', 'genome']
- **t1_h5_GSM8986593** source=`GSM8986593_GM23815_D10_filtered_feature_bc_matrix.h5` columns=['barcode', 'gene_id', 'symbol', 'feature_type', 'genome']
- **t2_gtex_obs** source=`<repo>\results\fibro\stage1_obs.csv` columns=['SAMPID', 'SMATSSCR', 'SMCENTER', 'SMPTHNTS', 'SMRIN', 'SMTS', 'SMTSD', 'SMUBRID', 'SMTSISCH', 'SMTSPAX', 'SMNABTCH', 'SMNABTCHT', 'SMNABTCHD', 'SMGEBTCH', 'SMGEBTCHD', 'SMGEBTCHT', 'ANALYTE_TYPE', 'SMAFRZE', 'SMGTC', 'SMRDTTL', 'SMALTTL', 'SMALTALG', 'SMSUPALG', 'SMRDLGTH', 'SMVQCFL', 'SMLMAPQ', 'SMUMPRD', 'SMUNPDRD', 'SMMPPD', 'SMMAPRT', 'SMMPPDUN', 'SMUNMPRT', 'SMMPDP', 'SMDPMPRT', 'SMMPPDXG', 'SMMPDPXG', 'SMDPRTXG', 'SMCHMRD', 'SMCHMRT', 'SMMPPDPR', 'SMMPHQRD', 'SMMPHQRT', 'SMMPLQRD', 'SMSPLTRT', 'SME1MPRD', 'SME2MPRD', 'SME1MPRT', 'SME2MPRT', 'SME1MMB', 'SME2MMB', 'SME1TTLB', 'SME2TTLB', 'SME1MMRT', 'SME2MMRT', 'SMTTLMM', 'SMTTLB', 'SMBSMMRT', 'SMESTLBS', 'SMEXNCRD', 'SMEXNCRT', 'SMEXPEFF', 'SMNTRNRD', 'SMNTRNRT', 'SMNTRARD', 'SMNTRART', 'SMNTERRD', 'SMNTERRT', 'SMAMBRD', 'SMAMBRT', 'SMNTEXC', 'SMDSCRT', 'SMEXNCRTHQ', 'SMNTRNRTHQ', 'SMNTRARTHQ', 'SMNTERRTHQ', 'SMAMBRTHQ', 'SME1SNSE', 'SME2SNSE', 'SME1ANTI', 'SME2ANTI', 'SME1PCTS', 'SME2PCTS', 'SMGNSDTC', 'SMRRNARD', 'SMRRNART', 'SMMFLGTH', 'SMSFLGTH', 'SMMDFLGTH', 'SMSMFLGTH', 'SMFGCMN', 'SMFGCSD', 'SMFGCSK', 'SMFGCKT', 'SM3PBMN', 'SM3PBSD', 'SM3PBMD', 'SM3PB25P', 'SM3PB75P', 'SM3PBSDM', 'SM3PBGN', 'SMMDMNCV', 'SMMDCVSD', 'SMMDCVCV', 'SMEXCVMD', 'SMEXCVMAD', 'SMMNCV', 'SMUVCRD', 'SMUVCRT', 'SMSHRTRD', 'SMSHRTRT', 'SMSMRDHQ', 'SMSMRTHQ', 'SMPRERDHQ', 'SMPRERTHQ', 'SMSMGNDT', 'SMPREGNDT', 'SMRDLNMN', 'SMRDLNMD', 'SMRDLNSD', 'donor', 'SUBJID', 'SEX', 'AGE', 'DTHHRDY', 'age_mid', 'age_ordinal']
- **t2_t1_z** source=`<repo>\results\fibro3\t1_cluster_z.npz` columns=['Z', 'age', 'pluri_with', 'pluri_without', 'counts', 'logcpm', 'missing', 'cell_line', 'day', 'cluster', 'n_cells', 'below_min_cells']
- **t2_t1_table** source=`<repo>\results\fibro3\t1_cluster_table.csv` columns=['cell_line', 'day', 'cluster', 'age_years', 'gsm', 'age_source', 'n_cells', 'n_timepoint', 'frac_timepoint', 'median_UMI', 'median_genes', 'mito_frac_median_cell', 'below_min_cells', 'age_score', 'pluri_with_OSKM', 'pluri_without_OSKM', 'pluri_primary', 'pluri_primary_name']
- **t3_t1_z** source=`<repo>\results\fibro3\t1_cluster_z.npz` columns=['Z', 'age', 'pluri_with', 'pluri_without', 'counts', 'logcpm', 'missing', 'cell_line', 'day', 'cluster', 'n_cells', 'below_min_cells']
- **t3_t1_table** source=`<repo>\results\fibro3\t1_cluster_table.csv` columns=['cell_line', 'day', 'cluster', 'age_years', 'gsm', 'age_source', 'n_cells', 'n_timepoint', 'frac_timepoint', 'median_UMI', 'median_genes', 'mito_frac_median_cell', 'below_min_cells', 'age_score', 'pluri_with_OSKM', 'pluri_without_OSKM', 'pluri_primary', 'pluri_primary_name']


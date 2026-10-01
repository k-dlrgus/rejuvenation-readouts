# FINDINGS_GTEX — GTEx cortex: is the age direction stable when logistics are recorded?

**Status:** Stage 1 complete (S2_i). Seed `20260914`. Release GTEx Analysis V10 (RNASeQCv2.4.2; phs000424.v10). Primary transfer metric: held-out Pearson r (n_perm=200).

Does not modify `FINDINGS_LOWDIM.md`, `FINDINGS_GEOMETRY.md`, `FINDINGS_POSCTRL.md`, `FINDINGS_TARGET.md`, `FINDINGS_TRAJECTORY.md`, `FINDINGS_BRAIN_PHASE1.md`, or `FALSIFICATION.md`. No perturbation data. No TF Atlas. No candidate interventions.

Reproduced by `notebooks/gtex_stage0.ipynb`, `gtex_stage1.ipynb`. RNA-seq analysis freeze: `SMAFRZE==RNASEQ`. Six tissues, no more.

## Caveats (GTEx, before Stage 1)

1. Public `AGE` is a 10-year bin. Use bin midpoints; report Spearman with the bin as an ordinal alongside R². Year-level MAE is
   not meaningful. Exact ages need dbGaP (`phs000424.v10`) — out of scope here.
2. Bulk tissue. A stable bulk age direction can be compositional (the blood lesson). Composition is handled explicitly below.
3. Per the GTEx FAQ, `Brain - Cortex` was sampled at the collection site and preserved in PAXgene; `Brain - Frontal Cortex (BA9)`
   and the other brain sub-regions were sampled later at the Miami Brain Endowment Bank and snap-frozen (longer ischemic time).
   So `SMCENTER` is the collection site, not the dissection site, for BA9. The BA9/Cortex pair on the same donors is a
   same-region, same-donor, different-preservation contrast — use it as such.

## Stage 0 — cohort

| tissue | n_samples | n_donors | n_SMCENTER | n_SMNABTCH | n_SMGEBTCH | n_male | n_female |
|---|---:|---:|---:|---:|---:|---:|---:|
| Brain - Frontal Cortex (BA9) | 269 | 269 | 4 | 147 | 108 | 201 | 68 |
| Brain - Cortex | 270 | 270 | 3 | 164 | 141 | 195 | 75 |
| Brain - Anterior cingulate cortex (BA24) | 233 | 233 | 4 | 146 | 105 | 177 | 56 |
| Brain - Hippocampus | 255 | 255 | 4 | 145 | 111 | 190 | 65 |
| Heart - Left Ventricle | 452 | 452 | 3 | 269 | 204 | 308 | 144 |
| Muscle - Skeletal | 818 | 818 | 3 | 338 | 226 | 552 | 266 |

AGE bins (n samples):

| tissue | 20-29 | 30-39 | 40-49 | 50-59 | 60-69 | 70-79 |
|---|---|---|---|---|---|---|
| Brain - Frontal Cortex (BA9) | 5 | 8 | 23 | 83 | 138 | 12 |
| Brain - Cortex | 7 | 8 | 25 | 84 | 129 | 17 |
| Brain - Anterior cingulate cortex (BA24) | 7 | 7 | 25 | 62 | 119 | 13 |
| Brain - Hippocampus | 4 | 8 | 25 | 75 | 128 | 15 |
| Heart - Left Ventricle | 20 | 25 | 67 | 159 | 167 | 14 |
| Muscle - Skeletal | 71 | 64 | 121 | 261 | 273 | 28 |

| tissue | SMRIN median (min–max, miss) | SMTSISCH median (min–max, miss) |
|---|---|---|
| Brain - Frontal Cortex (BA9) | 7.10 (5.50–9.40, miss=0) | 766.0 (223.0–1350.0, miss=0) |
| Brain - Cortex | 6.50 (5.40–8.70, miss=0) | 1004.0 (383.0–1668.0, miss=0) |
| Brain - Anterior cingulate cortex (BA24) | 6.60 (5.50–9.40, miss=0) | 790.0 (223.0–1375.0, miss=0) |
| Brain - Hippocampus | 6.70 (5.50–9.20, miss=0) | 762.0 (223.0–1375.0, miss=0) |
| Heart - Left Ventricle | 7.25 (5.50–9.40, miss=0) | 492.0 (55.0–1737.0, miss=1) |
| Muscle - Skeletal | 8.00 (5.50–9.90, miss=0) | 573.0 (35.0–2075.0, miss=0) |

Full SMCENTER / SMNABTCH / DTHHRDY / SEX counts: `results/gtex/stage0_cohort.md`.

## Stage 0 STOP checks

- BA9 n after freeze filter: **269** (underpowered vs POSCTRL n=150 bar if n<150: False).
- SMCENTER levels with ≥25 BA9 samples: **2** (T2 not testable if <3: True).
- BA9 SMCENTER counts: `{'B1, A1': 149, 'C1, A1': 113, 'D1, A1': 4, 'C1, B1, A1': 3}`. Values are the file strings; BA9/BA24/hippocampus are comma-joined (`B1, A1`); they are not split.

## Stage 0 — covariate audit

OLS R² of age-bin midpoint on each covariate, descriptive, no CV. Categorical covariates one-hot (drop first). Joint = all six together. n is complete cases for that regression; incomplete rows are not dropped from the cohort. L = number of levels (SMNABTCH L is large relative to n; that R² is a saturated-batch fit, not a logistics slope).

| tissue | SMRIN | SMTSISCH | DTHHRDY | SEX | SMCENTER | SMNABTCH | JOINT |
|---|---|---|---|---|---|---|---|
| Brain - Frontal Cortex (BA9) | +0.000 (n=269) | +0.007 (n=269) | +0.093 (n=269, L=5) | +0.003 (n=269, L=2) | +0.010 (n=269, L=4) | +0.601 (n=269, L=147) | +0.640 (n=269, L=156) |
| Brain - Cortex | +0.010 (n=270) | +0.024 (n=270) | +0.129 (n=270, L=5) | +0.002 (n=270, L=2) | +0.007 (n=270, L=3) | +0.641 (n=270, L=164) | +0.688 (n=270, L=172) |
| Brain - Anterior cingulate cortex (BA24) | +0.003 (n=233) | +0.016 (n=233) | +0.138 (n=233, L=5) | +0.007 (n=233, L=2) | +0.017 (n=233, L=4) | +0.671 (n=233, L=146) | +0.716 (n=233, L=155) |
| Brain - Hippocampus | +0.010 (n=255) | +0.016 (n=255) | +0.086 (n=255, L=5) | +0.005 (n=255, L=2) | +0.005 (n=255, L=4) | +0.554 (n=255, L=145) | +0.619 (n=255, L=154) |
| Heart - Left Ventricle | +0.070 (n=452) | +0.040 (n=451) | +0.121 (n=451, L=5) | +0.007 (n=452, L=2) | +0.000 (n=452, L=3) | +0.647 (n=452, L=269) | +0.673 (n=450, L=277) |
| Muscle - Skeletal | +0.001 (n=818) | +0.054 (n=818) | +0.186 (n=808, L=5) | +0.001 (n=818, L=2) | +0.017 (n=818, L=3) | +0.429 (n=818, L=338) | +0.540 (n=808, L=346) |

## Stage 0 — gene filter and DLPFC overlap

| tissue | n_samples | n_genes_raw | n_genes_filter | n_overlap_dlpfc |
|---|---:|---:|---:|---:|
| Brain - Frontal Cortex (BA9) | 269 | 59033 | 26847 | 20115 |
| Brain - Cortex | 270 | 59033 | 27151 | 20195 |
| Brain - Anterior cingulate cortex (BA24) | 233 | 59033 | 26595 | 20059 |
| Brain - Hippocampus | 255 | 59033 | 26657 | 19949 |
| Heart - Left Ventricle | 452 | 59033 | 23818 | 18226 |
| Muscle - Skeletal | 818 | 59033 | 23228 | 17943 |

Filter (frozen): ≥6 counts in ≥20% of samples per tissue. Then TMM (edgeR defaults) → log2-CPM (`prior.count=2`). z-score per gene on training rows only inside every fold (Stage 1).

## Stage 0 — composition proxies (brain tissues)

Neuronal markers (frozen): RBFOX3, SNAP25, SYT1, GAD1, SLC17A7. Glial markers (frozen a priori): GFAP, AQP4, MBP, PLP1, OLIG2, CX3CR1, AIF1. Score = sum of gene-z log-CPM over the set (full-tissue z, Stage 0 descriptive).

| tissue | score | Spearman age_mid | Spearman age ordinal | Spearman SMTSISCH | R² ~ C(SMCENTER) |
|---|---|---:|---:|---:|---:|
| Brain - Frontal Cortex (BA9) | neuronal | -0.193 | -0.193 | +0.012 | +0.014 |
| Brain - Frontal Cortex (BA9) | glial | +0.191 | +0.191 | +0.105 | +0.027 |
| Brain - Cortex | neuronal | -0.251 | -0.251 | -0.009 | +0.012 |
| Brain - Cortex | glial | +0.170 | +0.170 | +0.024 | +0.029 |
| Brain - Anterior cingulate cortex (BA24) | neuronal | -0.177 | -0.177 | +0.122 | +0.006 |
| Brain - Anterior cingulate cortex (BA24) | glial | +0.136 | +0.136 | -0.076 | +0.006 |
| Brain - Hippocampus | neuronal | -0.187 | -0.187 | +0.085 | +0.007 |
| Brain - Hippocampus | glial | +0.207 | +0.207 | -0.117 | +0.017 |

## Stage 1 gate

- FINDINGS_POSCTRL.md exists: True
- C1 ran: True (source=s2_c1_dense_shared_r030.json+s2_gene_ridge.json)
- post-supersession status: S2_i
- open: True
- reason: gate open: S2 outcome (i) — planted and real-age r transfer both ways (gene_ridge r HBCC→MSSM +0.568 / MSSM→HBCC +0.578; uncalibrated R² is a calibration artifact)

S2 transfer used for the gate (Pearson r): real-age gene_ridge HBCC→MSSM +0.568 (null -0.005), MSSM→HBCC +0.578 (null -0.022); uncalibrated R² -0.124 / -0.463 (calibration artifact). Planted dense_shared r +0.865 / +0.768 (null +0.015 / +0.013).

## Pre-registration 2026-09-17 (verbatim, written before any Stage 1 fit)

Written **before any Stage 1 fit**, 2026-09-17. No threshold in this block is re-tuned after numbers exist.

FINDINGS_POSCTRL.md Supersession 2 fired outcome (i): the DLPFC age direction transfers between banks as a correlation
(gene_ridge r +0.568 / +0.578, nulls ≈ 0); uncalibrated transfer R² is a calibration artifact. Stage 1 gate is open.

Primary transfer metric is held-out Pearson r between predicted and true age-bin midpoint, with a donor-level permutation
null within `SMCENTER` (`n_perm=200`). **Pass = r > 0 with null ≤ 0.05.** Also report uncalibrated R², calibrated R²
(intercept and slope refit on the test set), and Spearman ρ. Bootstrap median pairwise angle is **reported** against the
POSCTRL reference values (sex control 41.3°, planted dense_shared at C3 n=233: 41.191° — nearest n to BA9 n=269), not
against 35°. The 35° bar is retired. **Angle is not a pass/fail gate.**

POSCTRL planted `dense_shared` transferred as r at n=150 (10/10 draws) and n=233 (1/1). BA9 n=269 is at or above that n.

The original Prompt B bar (transfer R²) was never evaluated on GTEx numbers; it is superseded by this block before any
T-cell is fit.

## Pre-registered reading (verbatim, r in place of R²)

- BA9 raw passes transfer (held-out Pearson r > 0 with permutation null ≤ 0.05) AND resid does not change it materially
  (sign of r unchanged, Δangle < 5°) → a stable bulk cortex age direction exists and is not driven by recorded logistics;
  the DLPFC instability is then attributable to something GTEx does not share: two-bank structure, snRNA/pseudobulk noise,
  or unrecorded covariates. **Say all three; do not pick one.**
- BA9 raw fails, resid passes → recorded logistics (RIN / ischemic time / Hardy / composition) were hiding a stable
  direction. This makes unrecorded PMI/RIN a plausible, not proven, cause of the DLPFC failure.
- Both fail at BA9 n ≥ the n where POSCTRL's planted direction transferred as r → cortex age directions fail transfer in a
  second, differently-structured cohort. The negative generalizes. Say so.
- Both fail at n below that → uninformative; report as such.
- T3 fails while T2 passes → the preservation pathway defeats transfer even when collection center does not. Report as a
  distinct finding.

## Original Prompt B bars (superseded before any GTEx Stage 1 number)

primary bar is **transfer R² > 0 with permutation null ≤ 0.05**; bootstrap median
pairwise angle is **reported against the POSCTRL reference values** (sex control 41.3°, and the planted dense_shared angle at the
nearest n from the C3 power curve), not against 35°. The 35° bar is retired.

- BA9 raw passes transfer AND resid does not change it materially (sign unchanged, Δangle < 5°) → a stable bulk cortex age direction
  exists and is not driven by recorded logistics; the DLPFC instability is then attributable to something GTEx does not share:
  two-bank structure, snRNA/pseudobulk noise, or unrecorded covariates. **Say all three; do not pick one.**
- BA9 raw fails, resid passes → recorded logistics (RIN / ischemic time / Hardy / composition) were hiding a stable direction. This
  makes unrecorded PMI/RIN a plausible, not proven, cause of the DLPFC failure.
- Both fail at BA9 n ≥ the n where POSCTRL's planted direction transferred → cortex age directions fail transfer in a second,
  differently-structured cohort. The negative generalizes. Say so.
- Both fail at n below that → uninformative; report as such.
- T3 fails while T2 passes → the preservation pathway defeats transfer even when collection center does not. Report as a distinct
  finding.

## T1

| test | tissue | regime | method | boot_median_angle | boot_p05 | boot_p95 | boot_n_ok | within_batch_r_median_over_centers | within_batch_r_null | within_batch_r_p | within_batch_r2_median_over_centers | within_batch_r2_null | within_batch_r2_p | within_batch_rho | within_batch_rho_null | within_batch_cal_r2 | within_batch_cal_r2_null | pass_r | n_perm | sex_angle_ref | planted_angle_n233 | n_centers |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T1 | Brain - Frontal Cortex (BA9) | raw | ridge | +46.129 | +40.960 | +51.312 | 20 | +0.508 | -0.006 | +0.005 | +0.268 | -0.201 | +0.005 | +0.438 | -0.007 | +0.279 | +0.012 | True | 200 | +41.300 | +41.191 | 2 |
| T1 | Brain - Frontal Cortex (BA9) | raw | pls1 | +38.223 | +26.079 | +50.940 | 20 | +0.192 | -0.016 | +0.015 | +0.003 | -0.175 | +0.015 | +0.189 | -0.025 | +0.055 | +0.012 | True | 200 | +41.300 | +41.191 | 2 |
| T1 | Brain - Frontal Cortex (BA9) | resid | ridge | +49.107 | +42.650 | +55.540 | 20 | +0.327 | -0.032 | +0.005 | +0.102 | -0.160 | +0.005 | +0.252 | -0.031 | +0.122 | +0.009 | True | 200 | +41.300 | +41.191 | 2 |
| T1 | Brain - Frontal Cortex (BA9) | resid | pls1 | +50.089 | +39.103 | +59.331 | 20 | +0.143 | -0.025 | +0.030 | -0.141 | -0.289 | +0.070 | +0.096 | -0.029 | +0.053 | +0.011 | True | 200 | +41.300 | +41.191 | 2 |
| T1 | Brain - Cortex | raw | ridge | +46.702 | +39.910 | +53.295 | 20 | +0.421 | -0.027 | +0.005 | +0.184 | -0.209 | +0.005 | +0.391 | -0.020 | +0.217 | +0.010 | True | 200 | +41.300 | +41.191 | 2 |
| T1 | Brain - Cortex | raw | pls1 | +30.387 | +23.721 | +40.209 | 20 | +0.119 | -0.012 | +0.065 | -0.103 | -0.183 | +0.179 | +0.128 | -0.010 | +0.076 | +0.011 | True | 200 | +41.300 | +41.191 | 2 |
| T1 | Brain - Cortex | resid | ridge | +48.072 | +41.711 | +53.783 | 20 | +0.382 | -0.015 | +0.005 | +0.141 | -0.146 | +0.005 | +0.387 | -0.015 | +0.149 | +0.008 | True | 200 | +41.300 | +41.191 | 2 |
| T1 | Brain - Cortex | resid | pls1 | +47.630 | +38.757 | +58.673 | 20 | +0.273 | -0.017 | +0.005 | +0.005 | -0.271 | +0.005 | +0.288 | -0.015 | +0.076 | +0.010 | True | 200 | +41.300 | +41.191 | 2 |
| T1 | Brain - Anterior cingulate cortex (BA24) | raw | ridge | +45.822 | +39.408 | +51.888 | 20 | +0.522 | -0.020 | +0.005 | +0.278 | -0.213 | +0.005 | +0.470 | -0.022 | +0.287 | +0.013 | True | 200 | +41.300 | +41.191 | 2 |
| T1 | Brain - Anterior cingulate cortex (BA24) | raw | pls1 | +36.452 | +27.727 | +47.941 | 20 | +0.279 | -0.014 | +0.005 | +0.037 | -0.181 | +0.005 | +0.272 | -0.023 | +0.080 | +0.012 | True | 200 | +41.300 | +41.191 | 2 |
| T1 | Brain - Anterior cingulate cortex (BA24) | resid | ridge | +48.311 | +42.459 | +54.715 | 20 | +0.397 | -0.035 | +0.005 | +0.154 | -0.164 | +0.005 | +0.356 | -0.036 | +0.157 | +0.012 | True | 200 | +41.300 | +41.191 | 2 |
| T1 | Brain - Anterior cingulate cortex (BA24) | resid | pls1 | +45.784 | +37.938 | +56.689 | 20 | +0.261 | -0.014 | +0.005 | +0.009 | -0.275 | +0.005 | +0.231 | -0.017 | +0.070 | +0.013 | True | 200 | +41.300 | +41.191 | 2 |
| T1 | Brain - Hippocampus | raw | ridge | +48.857 | +42.095 | +55.465 | 20 | +0.454 | -0.016 | +0.005 | +0.199 | -0.242 | +0.005 | +0.470 | -0.015 | +0.212 | +0.010 | True | 200 | +41.300 | +41.191 | 2 |
| T1 | Brain - Hippocampus | raw | pls1 | +34.597 | +27.854 | +47.419 | 20 | +0.281 | -0.018 | +0.005 | +0.047 | -0.205 | +0.005 | +0.276 | -0.018 | +0.087 | +0.012 | True | 200 | +41.300 | +41.191 | 2 |
| T1 | Brain - Hippocampus | resid | ridge | +50.492 | +43.706 | +58.845 | 20 | +0.285 | -0.029 | +0.005 | +0.057 | -0.176 | +0.005 | +0.294 | -0.026 | +0.081 | +0.011 | True | 200 | +41.300 | +41.191 | 2 |
| T1 | Brain - Hippocampus | resid | pls1 | +44.244 | +35.430 | +54.308 | 20 | +0.250 | -0.008 | +0.005 | -0.014 | -0.279 | +0.005 | +0.264 | -0.013 | +0.063 | +0.012 | True | 200 | +41.300 | +41.191 | 2 |
| T1 | Heart - Left Ventricle | raw | ridge | +47.753 | +43.137 | +52.592 | 20 | +0.450 | -0.004 | +0.005 | +0.200 | -0.199 | +0.005 | +0.419 | -0.000 | +0.213 | +0.006 | True | 200 | +41.300 | +41.191 | 2 |
| T1 | Heart - Left Ventricle | raw | pls1 | +24.742 | +20.152 | +32.037 | 20 | +0.239 | -0.003 | +0.005 | +0.036 | -0.127 | +0.005 | +0.228 | -0.007 | +0.058 | +0.006 | True | 200 | +41.300 | +41.191 | 2 |
| T1 | Heart - Left Ventricle | resid | ridge | +47.670 | +42.983 | +51.390 | 20 | +0.220 | -0.019 | +0.005 | +0.041 | -0.125 | +0.005 | +0.209 | -0.023 | +0.073 | +0.007 | True | 200 | +41.300 | +41.191 | 2 |
| T1 | Heart - Left Ventricle | resid | pls1 | +54.728 | +44.534 | +71.943 | 20 | +0.048 | -0.006 | +0.259 | -0.109 | -0.176 | +0.164 | +0.091 | -0.009 | +0.003 | +0.008 | True | 200 | +41.300 | +41.191 | 2 |
| T1 | Muscle - Skeletal | raw | ridge | +45.519 | +42.891 | +48.276 | 20 | +0.708 | -0.010 | +0.005 | +0.503 | -0.211 | +0.005 | +0.662 | -0.009 | +0.505 | +0.004 | True | 200 | +41.300 | +41.191 | 2 |
| T1 | Muscle - Skeletal | raw | pls1 | +18.031 | +15.644 | +21.347 | 20 | +0.407 | -0.008 | +0.005 | +0.162 | -0.111 | +0.005 | +0.408 | -0.011 | +0.166 | +0.003 | True | 200 | +41.300 | +41.191 | 2 |
| T1 | Muscle - Skeletal | resid | ridge | +43.445 | +40.598 | +46.445 | 20 | +0.531 | -0.001 | +0.005 | +0.278 | -0.108 | +0.005 | +0.483 | -0.004 | +0.287 | +0.003 | True | 200 | +41.300 | +41.191 | 2 |
| T1 | Muscle - Skeletal | resid | pls1 | +29.615 | +24.853 | +34.823 | 20 | +0.288 | -0.003 | +0.005 | +0.072 | -0.138 | +0.005 | +0.268 | -0.006 | +0.083 | +0.003 | True | 200 | +41.300 | +41.191 | 2 |

## T2

SKIPPED: t2_not_testable: fewer than 3 SMCENTER with ≥25 BA9 samples

## T3 / T4

| test | label | train_tissue | test_tissue | regime | method | r | r_null | r_p | r2 | r2_null | cal_r2 | cal_r2_null | rho | rho_null | pass_r | n_train | n_test | n_genes | n_overlap_donors_dropped | n_perm |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T3 | Brain - Frontal Cortex (BA9) -> Brain - Cortex | Brain - Frontal Cortex (BA9) | Brain - Cortex | raw | ridge | +0.480 | +0.020 | +0.005 | +0.202 | -0.171 | +0.230 | +0.010 | +0.441 | +0.012 | True | 75 | 270 | 26590 | 194 | 200 |
| T3 | Brain - Frontal Cortex (BA9) -> Brain - Cortex | Brain - Frontal Cortex (BA9) | Brain - Cortex | raw | pls1 | +0.244 | +0.102 | +0.114 | +0.040 | -0.071 | +0.059 | +0.029 | +0.251 | +0.130 | False | 75 | 270 | 26590 | 194 | 200 |
| T3 | Brain - Frontal Cortex (BA9) -> Brain - Cortex | Brain - Frontal Cortex (BA9) | Brain - Cortex | resid | ridge | +0.377 | +0.031 | +0.005 | +0.097 | -0.116 | +0.142 | +0.008 | +0.380 | +0.034 | True | 75 | 270 | 26590 | 194 | 200 |
| T3 | Brain - Frontal Cortex (BA9) -> Brain - Cortex | Brain - Frontal Cortex (BA9) | Brain - Cortex | resid | pls1 | +0.226 | +0.048 | +0.045 | -0.005 | -0.199 | +0.051 | +0.009 | +0.249 | +0.040 | True | 75 | 270 | 26590 | 194 | 200 |
| T3 | Brain - Cortex -> Brain - Frontal Cortex (BA9) | Brain - Cortex | Brain - Frontal Cortex (BA9) | raw | ridge | +0.407 | -0.003 | +0.005 | -0.157 | -0.385 | +0.165 | +0.008 | +0.318 | -0.002 | True | 76 | 269 | 26590 | 194 | 200 |
| T3 | Brain - Cortex -> Brain - Frontal Cortex (BA9) | Brain - Cortex | Brain - Frontal Cortex (BA9) | raw | pls1 | +0.251 | +0.074 | +0.005 | -0.181 | -0.271 | +0.063 | +0.024 | +0.227 | +0.068 | False | 76 | 269 | 26590 | 194 | 200 |
| T3 | Brain - Cortex -> Brain - Frontal Cortex (BA9) | Brain - Cortex | Brain - Frontal Cortex (BA9) | resid | ridge | +0.235 | +0.035 | +0.015 | -0.330 | -0.263 | +0.055 | +0.005 | +0.203 | +0.033 | True | 76 | 269 | 26590 | 194 | 200 |
| T3 | Brain - Cortex -> Brain - Frontal Cortex (BA9) | Brain - Cortex | Brain - Frontal Cortex (BA9) | resid | pls1 | +0.187 | +0.022 | +0.010 | -0.458 | -0.457 | +0.035 | +0.006 | +0.131 | +0.023 | True | 76 | 269 | 26590 | 194 | 200 |
| T4 | Brain - Frontal Cortex (BA9) -> Brain - Anterior cingulate cortex (BA24) | Brain - Frontal Cortex (BA9) | Brain - Anterior cingulate cortex (BA24) | raw | ridge | +0.360 | -0.045 | +0.005 | +0.092 | -0.188 | +0.129 | +0.012 | +0.319 | -0.014 | True | 84 | 233 | 26274 | 185 | 200 |
| T4 | Brain - Frontal Cortex (BA9) -> Brain - Anterior cingulate cortex (BA24) | Brain - Frontal Cortex (BA9) | Brain - Anterior cingulate cortex (BA24) | raw | pls1 | +0.154 | +0.007 | +0.194 | -0.004 | -0.097 | +0.024 | +0.020 | +0.147 | +0.011 | True | 84 | 233 | 26274 | 185 | 200 |
| T4 | Brain - Frontal Cortex (BA9) -> Brain - Anterior cingulate cortex (BA24) | Brain - Frontal Cortex (BA9) | Brain - Anterior cingulate cortex (BA24) | resid | ridge | +0.226 | +0.024 | +0.100 | +0.048 | -0.138 | +0.051 | +0.012 | +0.195 | +0.016 | True | 84 | 233 | 26274 | 185 | 200 |
| T4 | Brain - Frontal Cortex (BA9) -> Brain - Anterior cingulate cortex (BA24) | Brain - Frontal Cortex (BA9) | Brain - Anterior cingulate cortex (BA24) | resid | pls1 | +0.045 | +0.007 | +0.433 | -0.244 | -0.183 | +0.002 | +0.010 | +0.004 | +0.009 | True | 84 | 233 | 26274 | 185 | 200 |
| T4 | Brain - Frontal Cortex (BA9) -> Brain - Hippocampus | Brain - Frontal Cortex (BA9) | Brain - Hippocampus | raw | ridge | +0.334 | -0.010 | +0.005 | +0.097 | -0.427 | +0.111 | +0.007 | +0.286 | -0.003 | True | 67 | 255 | 25728 | 202 | 200 |
| T4 | Brain - Frontal Cortex (BA9) -> Brain - Hippocampus | Brain - Frontal Cortex (BA9) | Brain - Hippocampus | raw | pls1 | +0.251 | +0.057 | +0.080 | -0.018 | -0.222 | +0.063 | +0.031 | +0.214 | +0.038 | False | 67 | 255 | 25728 | 202 | 200 |
| T4 | Brain - Frontal Cortex (BA9) -> Brain - Hippocampus | Brain - Frontal Cortex (BA9) | Brain - Hippocampus | resid | ridge | +0.021 | -0.013 | +0.368 | -0.653 | -0.405 | +0.000 | +0.006 | +0.069 | -0.020 | True | 67 | 255 | 25728 | 202 | 200 |
| T4 | Brain - Frontal Cortex (BA9) -> Brain - Hippocampus | Brain - Frontal Cortex (BA9) | Brain - Hippocampus | resid | pls1 | -0.103 | -0.025 | +0.711 | -0.673 | -0.434 | +0.011 | +0.011 | -0.071 | -0.029 | False | 67 | 255 | 25728 | 202 | 200 |

## T5

```
{
 "n_common_genes": 20945,
 "identity_rank": 5,
 "identity_sv": [
  196.076038519334,
  103.21330954875826,
  62.547256494037455,
  37.517359074661265,
  22.58801759230702
 ],
 "ba9_angle_to_identity": 62.97933952920945,
 "ba9_proj_r2": 0.20639917631617533,
 "null_angle_mean": 70.35209728027063,
 "null_angle_median": 69.87495134545492,
 "p_angle": 0.10945273631840796,
 "pairwise_median_deg": 68.64961381794926,
 "pairwise_null_mean": 84.89764363242153,
 "p_pairwise": 0.004975124378109453,
 "tissues": [
  "Brain - Frontal Cortex (BA9)",
  "Brain - Cortex",
  "Brain - Anterior cingulate cortex (BA24)",
  "Brain - Hippocampus",
  "Heart - Left Ventricle",
  "Muscle - Skeletal"
 ],
 "n_perm": 200,
 "note": "bulk analogue of FINDINGS_GEOMETRY Part C, not a replication"
}
```

## T6

```
{
 "ran": true,
 "n_overlap": 20115,
 "angle_deg": 80.19970901762603,
 "null_mean": 89.67042162311502,
 "null_median": 89.7175332213159,
 "p": 0.000999000999000999,
 "source": "v1_axes.npz",
 "w_key": "W_age_ridge",
 "n_types_collapsed": 20,
 "n_random": 1000,
 "note": "one number with a null; not evidence beyond closer-than-random"
}
```

## Composition-only rerun (BA9 raw passed)

| test | tissue | regime | method | boot_median_angle | boot_p05 | boot_p95 | boot_n_ok | within_batch_r_median_over_centers | within_batch_r_null | within_batch_r_p | within_batch_r2_median_over_centers | within_batch_r2_null | within_batch_r2_p | within_batch_rho | within_batch_rho_null | within_batch_cal_r2 | within_batch_cal_r2_null | pass_r | n_perm | sex_angle_ref | planted_angle_n233 | n_centers |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T1 | Brain - Frontal Cortex (BA9) | comp_only | ridge | +44.537 | +39.691 | +50.721 | 20 | +0.455 | -0.026 | +0.005 | +0.211 | -0.155 | +0.005 | +0.394 | -0.020 | +0.222 | +0.011 | True | 200 | +41.300 | +41.191 | 2 |
| T1 | Brain - Frontal Cortex (BA9) | comp_only | pls1 | +46.230 | +33.983 | +60.597 | 20 | +0.189 | -0.013 | +0.010 | -0.100 | -0.241 | +0.040 | +0.135 | -0.014 | +0.078 | +0.009 | True | 200 | +41.300 | +41.191 | 2 |

## Verdict

T2 not testable (fewer than 3 SMCENTER with ≥25 BA9 samples). BA9 transfer is T3 (BA9 ↔ Cortex). ridge vs pls1 T3 raw disagree (ridge pass=True, pls1 pass=False); do not average. BA9 raw passes transfer (held-out Pearson r > 0 with permutation null ≤ 0.05) AND resid does not change it materially (sign of r unchanged, Δangle=+2.98° < 5°). A stable bulk cortex age direction exists and is not driven by recorded logistics; the DLPFC instability is then attributable to something GTEx does not share: two-bank structure, snRNA/pseudobulk noise, or unrecorded covariates. Say all three; do not pick one. Composition-only T1 still passes; the direction is not only the neuronal/glial mix.

## What this changes in FINDINGS_LOWDIM L4 (quote, not edited)

**The negative result stands.** What would be needed to change it: more donors per cell type (so that even gene space is not p ≫ n, or so that a k ≪ n space is estimated from hundreds of independent brains rather than ~200); **more sites** (two banks are one degree of freedom of transfer, and that transfer failed); or a **different tissue** with less bank structure and more independent donors. Re-fitting more unsupervised variants on this same 233-donor DLPFC matrix will not make an underdetermined direction unique.

## Limitations

1. Public `AGE` is a 10-year bin. Use bin midpoints; report Spearman with the bin as an ordinal alongside R². Year-level MAE is
   not meaningful. Exact ages need dbGaP (`phs000424.v10`) — out of scope here.
2. Bulk tissue. A stable bulk age direction can be compositional (the blood lesson). Composition is handled explicitly below.
3. Per the GTEx FAQ, `Brain - Cortex` was sampled at the collection site and preserved in PAXgene; `Brain - Frontal Cortex (BA9)`
   and the other brain sub-regions were sampled later at the Miami Brain Endowment Bank and snap-frozen (longer ischemic time).
   So `SMCENTER` is the collection site, not the dissection site, for BA9. The BA9/Cortex pair on the same donors is a
   same-region, same-donor, different-preservation contrast — use it as such.

4. Two covariate regimes (`raw`, `resid`) are both primary; `resid` is not selected after seeing numbers.
5. Seed `20260914`. SVD-LOO ridge α grid and PLS-1 are the geometry/target methods called on a different X.
6. Public AGE is a bin; Spearman with the ordinal bin is reported alongside R². Year-level MAE is not used.
7. `SMCENTER` is collection site, not BA9 dissection site (Miami Brain Endowment Bank).
8. Exact ages require dbGaP `phs000424.v10` and are out of scope.

## Files

| path | content |
|---|---|
| `src/gtex_common.py`, `gtex_stage0.py`, `gtex_stage1.py`, `gtex_findings.py` | code |
| `notebooks/gtex_stage0.ipynb`, `gtex_stage1.ipynb` | runnable from a clean checkout |
| `results/gtex/` | manifest, cohort, audit, composition, figures |
| `FINDINGS_GTEX.md` | this file |
| `PROGRESS_GTEX.md` | resume state |

# FINDINGS_TRAJECTORY — geometric scores and the reprogramming trajectory

**Status:** G1 done; G2 done; G3 no_transfer. Seed `20260914`. Dataset ID `4442d412-91cb-4261-acca-8adf5fa04c11` (cached Aging_Cohort h5ad; not invented). Does not modify `FINDINGS_GEOMETRY.md`, `FINDINGS_BRAIN_PHASE1.md`, or `FALSIFICATION.md`. No gene selection: every gene in the Phase-1 matrix enters every model with a weight.

Reproduced by `notebooks/trajectory_g1.ipynb`, `trajectory_g2.ipynb`, `trajectory_g3.ipynb`, `trajectory_g4.ipynb`. Null usability line: site-stratified shuffle R² ≤ 0.05.

## Headline

**G1 age (within-site, site-stratified):** ridge R²=+0.317 (null -0.055); PLS-1 R²=+0.016 (null -0.156). **G1 identity (subspace nearest centroid):** acc=0.9961 (chance 0.0698; shuffle 0.0457). **G2a (within-site PLS-1):** angle=74.6° (null 79.3°, p=0.039). **G2b:** identity-coords R²=+0.026 vs age-dir +0.317; redundant=False. **G2c:** residual ridge R²=+0.300; collapsed=False. **G4 verdict: (d)** The frozen human DLPFC age score does not track chronological age in independent mouse brain aging (Tabula Muris Senis) after 1:1 ortholog mapping (13028/25526 genes mapped, 51.0% of frozen). OLS R² of mouse age ~ score = +0.130 (null +0.068, p=0.235); Spearman ρ=+0.343 (p=0.235). Per-type Spearman signs +2 / -6 (mixed signs; a donor-mean can cancel). The trajectory is not computed. Needed: a human adult cortical/DLPFC scRNA dataset (preferred) in which the frozen score tracks donor age, plus an adult cortical OSK or partial-reprogramming time course with cell-type labels; or a mouse cortical aging atlas on which this frozen score tracks age, plus an adult cortical (not SVZ, not E15.5) OSK time course.

## G1 — geometric scores

Cohort: 233 donors, 20 cell types, 25526 genes, 3323 pseudobulks. Sites: `{'H': 120, 'M': 113}`. Age score = projection onto a within-type supervised direction over **all genes** (PLS-1: covariance direction + 1D OLS to years; ridge: SVD-LOO α, linear predictor). Identity subspace = SVD of train type centroids; classification = nearest centroid in that subspace; identity scalar = distance to own type centroid. Directions and the identity basis are fit on **training folds only**.

**Primary age method: ridge** (the actual age clock). PLS-1 is the geometric companion; both are reported. If they differ by >0.05 in primary R² the difference is called material.

| scheme | model | r2 | mae | pooled_r2 | site_strat_r2 | null_mean | null_p95 | p | null_usable | n_perm |
|---|---|---|---|---|---|---|---|---|---|---|
| loso | age_pls1 | -0.354 | +16.725 | -0.354 | -0.753 | -0.901 | -0.630 | +0.048 | True | 20 |
| loso | age_ridge | +0.101 | +13.844 | +0.101 | -0.174 | -0.822 | -0.698 | +0.048 | True | 20 |
| loso | identity_nearest_centroid | +0.991 | +0.959 | +0.991 | NA | +0.045 | +0.056 | +0.048 | True | 20 |
| within_site | age_pls1 | +0.016 | +12.113 | +0.257 | +0.016 | -0.156 | -0.081 | +0.048 | True | 20 |
| within_site | age_ridge | +0.317 | +9.978 | +0.488 | +0.317 | -0.055 | -0.015 | +0.048 | True | 20 |
| within_site | identity_nearest_centroid | +0.996 | +0.991 | +0.996 | NA | +0.046 | +0.053 | +0.048 | True | 20 |

### within_site

- **age pls1:** primary R²=+0.016, MAE=12.11 y; pooled R²=+0.257, site-strat R²=+0.016; null mean=-0.156 (p95=-0.081, p=0.048, n_perm=20, usable=True)
- **age ridge:** primary R²=+0.317, MAE=9.98 y; pooled R²=+0.488, site-strat R²=+0.317; null mean=-0.055 (p95=-0.015, p=0.048, n_perm=20, usable=True)
- **identity nearest-centroid:** acc=0.9961, macro-F1=0.9908, chance=0.0698; null acc mean=0.0457, p=0.048, n_perm=20; mean dist-to-own-centroid=24.034

### loso

- **age pls1:** primary R²=-0.354, MAE=16.72 y; pooled R²=-0.354, site-strat R²=-0.753; null mean=-0.901 (p95=-0.630, p=0.048, n_perm=20, usable=True)
- **age ridge:** primary R²=+0.101, MAE=13.84 y; pooled R²=+0.101, site-strat R²=-0.174; null mean=-0.822 (p95=-0.698, p=0.048, n_perm=20, usable=True)
- **identity nearest-centroid:** acc=0.9910, macro-F1=0.9594, chance=0.0698; null acc mean=0.0452, p=0.048, n_perm=20; mean dist-to-own-centroid=27.925

PLS-1 vs ridge (absolute primary-R² difference): loso Δ=0.455 material=True; within_site Δ=0.300 material=True.

**STOP G1:** False — age direction produces usable held-out R²; null is not inflated.

## G2 — three orthogonality controls

G2a: angle of the age direction to the identity subspace, both fit on the training fold. G2b: within-type age R² from the identity coordinates vs from the age direction (redundancy if R²_id ≥ 0.8 × R²_age or R²_age − R²_id < 0.1). G2c: refit the age direction on expression with the identity subspace regressed out (collapse if residual R² ≤ 0.05 or residual < 0.5 × R²_age).

| control | scheme | method | value | null | p |
|---|---|---|---|---|---|
| G2a_angle | loso | pls1 | +74.738 | +79.453 | +0.039 |
| G2a_proj_r2 | loso | pls1 | +0.070 | +0.035 | +0.059 |
| G2a_angle | loso | ridge | +83.515 | +85.299 | +0.020 |
| G2a_proj_r2 | loso | ridge | +0.013 | +0.007 | +0.020 |
| G2b_id_coords_age_r2 | loso | ridge_on_Z | -0.258 | -0.834 | +0.020 |
| G2b_age_dir_r2 | loso | pls1 | -0.354 | NA | NA |
| G2c_residual_age_r2 | loso | pls1 | -0.329 | -0.931 | +0.020 |
| G2b_age_dir_r2 | loso | ridge | +0.101 | NA | NA |
| G2c_residual_age_r2 | loso | ridge | +0.063 | -0.886 | +0.020 |
| G2a_angle | within_site | pls1 | +74.585 | +79.266 | +0.039 |
| G2a_proj_r2 | within_site | pls1 | +0.071 | +0.036 | +0.039 |
| G2a_angle | within_site | ridge | +83.154 | +84.983 | +0.020 |
| G2a_proj_r2 | within_site | ridge | +0.014 | +0.008 | +0.020 |
| G2b_id_coords_age_r2 | within_site | ridge_on_Z | +0.026 | -0.049 | +0.020 |
| G2b_age_dir_r2 | within_site | pls1 | +0.016 | NA | NA |
| G2c_residual_age_r2 | within_site | pls1 | +0.025 | -0.153 | +0.020 |
| G2b_age_dir_r2 | within_site | ridge | +0.317 | NA | NA |
| G2c_residual_age_r2 | within_site | ridge | +0.300 | -0.119 | +0.020 |

### within_site

- **G2a pls1:** angle=74.6° (null 79.3°, p=0.039); proj R²=0.071 (null 0.036, p=0.039, n_perm=50)
- **G2a ridge:** angle=83.2° (null 85.0°, p=0.020); proj R²=0.014 (null 0.008, p=0.020, n_perm=50)
- **G2b:** identity-coords age R²=+0.026 vs ridge age-dir R²=+0.317 (ratio=0.082, Δ=+0.291, null=-0.049, p=0.020); redundant=False
- **G2c ridge:** residual age R²=+0.300 (null -0.119, p=0.020); PLS-1 residual R²=+0.025; collapsed=False

### loso

- **G2a pls1:** angle=74.7° (null 79.5°, p=0.039); proj R²=0.070 (null 0.035, p=0.059, n_perm=50)
- **G2a ridge:** angle=83.5° (null 85.3°, p=0.020); proj R²=0.013 (null 0.007, p=0.020, n_perm=50)
- **G2b:** identity-coords age R²=-0.258 vs ridge age-dir R²=+0.101 (ratio=-2.541, Δ=+0.359, null=-0.834, p=0.020); redundant=False
- **G2c ridge:** residual age R²=+0.063 (null -0.886, p=0.020); PLS-1 residual R²=-0.329; collapsed=False

**STOP G2:** False — identity coordinates do not predict age as well as the age direction; residual age R² survives.

## G3 — trajectory

Search seed `20260914`. Candidate accessions re-fetched: ['GSE276656', 'GSE224438', 'GSE271794']. Found: `{'GSE276656': True, 'GSE224438': True, 'GSE271794': True}`. GEO hits=81; CXG mouse-brain RNA datasets=32; CXG human OSK title matches=0.

The frozen human DLPFC age score does not track chronological age in independent mouse brain aging (Tabula Muris Senis) after 1:1 ortholog mapping (13028/25526 genes mapped, 51.0% of frozen). OLS R² of mouse age ~ score = +0.130 (null +0.068, p=0.235); Spearman ρ=+0.343 (p=0.235). Per-type Spearman signs +2 / -6 (mixed signs; a donor-mean can cancel). The trajectory is not computed. Needed: a human adult cortical/DLPFC scRNA dataset (preferred) in which the frozen score tracks donor age, plus an adult cortical OSK or partial-reprogramming time course with cell-type labels; or a mouse cortical aging atlas on which this frozen score tracks age, plus an adult cortical (not SVZ, not E15.5) OSK time course.

Ortholog mapping: human genes=25526, mapped 1:1=14239, used=13028 (0.510 of the frozen human genes).

Independent mouse aging validation: dataset `66ff82b4-9380-469c-bc4b-cfa08eacd325 + c08f8441-4a10-4748-872a-e70c0bcccdba` (Tabula Muris Senis brain non-myeloid + myeloid (CXG)). OLS R² of mouse age ~ frozen score = +0.130 (null +0.068, p=0.235); Spearman ρ=+0.343 (p=0.235). Per-type Spearman signs +2 / -6. Transfer usable=False.
 Cortex-subtissue sensitivity: n_donors=11, OLS R²=+0.002, Spearman=-0.254.

Per-type transfer (Tabula Muris Senis, frozen scores):

| celltype | n_donors | spearman | ols_r2 |
|---|---|---|---|
| GABAergic neuron | 13 | -0.553 | +0.155 |
| astrocyte | 11 | -0.178 | +0.045 |
| endothelial cell | 11 | -0.178 | +0.065 |
| microglial cell | 14 | -0.037 | +0.127 |
| oligodendrocyte | 12 | +0.419 | +0.181 |
| oligodendrocyte precursor cell | 11 | -0.637 | +0.417 |
| pericyte | 10 | -0.657 | +0.490 |
| perivascular macrophage | 10 | +0.560 | +0.117 |

**Dataset that would be needed:** Human adult DLPFC/cortex OSK or partial-reprogramming scRNA time course with numeric ages and cell-type labels; or a mouse cortical aging atlas on which the frozen score tracks age, plus an adult cortical (not SVZ, not E15.5) OSK time course.

## G4 — verdict

**(d)** The frozen human DLPFC age score does not track chronological age in independent mouse brain aging (Tabula Muris Senis) after 1:1 ortholog mapping (13028/25526 genes mapped, 51.0% of frozen). OLS R² of mouse age ~ score = +0.130 (null +0.068, p=0.235); Spearman ρ=+0.343 (p=0.235). Per-type Spearman signs +2 / -6 (mixed signs; a donor-mean can cancel). The trajectory is not computed. Needed: a human adult cortical/DLPFC scRNA dataset (preferred) in which the frozen score tracks donor age, plus an adult cortical OSK or partial-reprogramming time course with cell-type labels; or a mouse cortical aging atlas on which this frozen score tracks age, plus an adult cortical (not SVZ, not E15.5) OSK time course.

The four options were: (a) independent scores and a bending trajectory (safe window measured); (b) independent scores, trajectory does not bend (axes separable in aging, fused under reprogramming); (c) scores not independent (G2 failed); (d) scores independent but do not transfer to available reprogramming data.

## Limitations

1. Two brain banks only (HBCC, MSSM). LOSO is one df of transfer and is reported, not averaged away.
2. Unmeasured 6-plex hashing pools, PMI, RIN — as in FINDINGS_BRAIN_PHASE1.md.
3. Cross-species transfer: scores are fit on human DLPFC; mouse application requires ortholog mapping and loses unmapped genes. Transfer is not assumed — it is tested on an independent mouse aging dataset before any reprogramming plot.
4. Seed `20260914` (`numpy.random.default_rng`). Folds match Phase 1 / geometry.
5. No gene selection. p ≫ n; ridge uses SVD-LOO α. PLS-1 is the covariance direction.
6. Identity classification from a 19-dimensional centroid subspace of 20 types is expected to be near-perfect; that is a geometric fact, not evidence that identity is 'solved'.
7. Per-type age directions were largely type-specific in Part C (pairwise ~78°). A single frozen consensus direction is not assumed; G3 applies per-type weights after mapping mouse types.
8. Cross-species orthologs: NCBI strict 1:1 Ensembl pairs cover 14,239 / 25,526 Phase-1 genes (55.8%). Unmapped genes sit at the human mean (z = 0) in the frozen scores.
9. GSE224438 (mouse SVZ partial reprogramming) is the usable reprogramming accession found; whole-body 10x libraries multiplex young/old/OSK in one lane and RAW has no hashing barcodes. Those libraries, and the SVZ-targeted lanes, were **not scored** because the frozen age score failed the independent mouse aging transfer test.
10. GSE224438 is a 3-condition protocol contrast (young / old / old+OSK), not a multi-day time course. Even if transfer had passed, quadratic vs linear would be unidentified with three points.
11. GSE276656 (mPFC engram OSK multiome, 14.4 GB RAW) was not downloaded. GSE271794 is E15.5 development, not adult aging. No human brain OSK scRNA dataset was found on CELLxGENE (0 title matches).
12. Transfer is tested on Tabula Muris Senis brain (CXG `66ff82b4-…` + `c08f8441-…`) before any reprogramming plot. That test failed (OLS R²=+0.130, p=0.235; cortex subset Spearman=−0.254); no trajectory is reported.

## Files

| path | content |
|---|---|
| `src/trajectory_common.py`, `trajectory_g1.py`, `trajectory_g2.py`, `trajectory_g3_search.py`, `trajectory_g3.py`, `trajectory_findings.py` | code |
| `notebooks/trajectory_g1.ipynb`, `g2`, `g3`, `g4` | runnable from a clean checkout |
| `results/trajectory/` | tables, json, logs, `figures/` |
| `FINDINGS_TRAJECTORY.md` | this file |


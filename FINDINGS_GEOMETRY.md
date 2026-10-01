# FINDINGS_GEOMETRY — leakage audit, gene-count control, geometric reframe

**Status:** Part A resolved; Part B done; Part C done. Seed `20260914`. Dataset ID `4442d412-91cb-4261-acca-8adf5fa04c11` (cached Aging_Cohort h5ad; not invented). Does not modify `FINDINGS_BRAIN_PHASE1.md` or `FALSIFICATION.md`. Does not derive new A/I gene-set definitions.

Reproduced by `notebooks/geometry_partA.ipynb`, `geometry_partB.ipynb`, `geometry_partC.ipynb`. Null usability line: site-stratified shuffle R² ≤ 0.05.

## Headline

**Leakage.** The P4 within-site shuffle R² of +0.195 is real. It is not donor leakage and not a per-pseudobulk shuffle. It is the between-bank age gap (R²(age~site)=+0.264) recovered by a site-specific intercept when OOF R² is pooled across HBCC and MSSM. Site-stratified R² puts that null at -0.232.

**Corrected P4c.** Within-site (site-stratified): R²_A=+0.250 vs size-matched R²_I=+0.413, shuffle=-0.232, verdict **falsified (not separable)**. LOSO (original pooled metric): R²_A=-0.016 vs R²_I=+0.279, shuffle=-1.005, verdict **falsified (not separable)**. P4c change: **STRENGTHEN: soundness-fail cleared; separability still falsified on the cross-test**.

**B3.** A-genes at matched size do **not** beat random (within-site R²_A=+0.250 vs random median +0.254, p95 +0.298). I-genes **do** (R²_I=+0.424). So the A-bucket is empty as an age axis (methods / weak-signal aggregation), while identity genes carry extra age signal (biology). This is not clean (a) or clean (b).

**C3.** Median-over-types angle of the PLS-1 age direction to the identity subspace: **66.2°** (permutation null mean 70.9°, p=0.030). Projection R²=0.163 (null 0.108). Ridge sensitivity: 83.3°. Geometry: **largely orthogonal** (a small excess overlap vs null, p=0.03, is not 'substantially within' the identity subspace). Per-type age directions are **largely type-specific** (pairwise ~78°, only ~12° from orthogonal; PC1 0.40 vs null 0.33 is a weak common component, not a shared axis).

**Claim the evidence now supports:** After correcting the site-pooled R² artifact, the gene-set cross-test is within-site **falsified (not separable)** and LOSO **falsified (not separable)** (STRENGTHEN: soundness-fail cleared; separability still falsified on the cross-test). **Genes cannot be partitioned into clean age vs identity buckets:** A-genes match size-matched random sets; I-genes beat both, so chronological age is carried by identity programmes rather than by a privileged A-set. That is a methods finding about the construct, plus a biological finding that identity genes are age-informative. **The axes are not geometrically fused.** The within-type age direction is largely orthogonal to the identity (centroid) subspace (66° vs null 71°; ridge 83°). Gene-set non-partitionability and geometric fusion are different claims; the evidence now supports the first, not the second. Per-type age directions are **largely type-specific** (pairwise ~78°, only ~12° from orthogonal; PC1 0.40 vs null 0.33 is a weak common component, not a shared axis).

## Part A — Leakage audit

Cohort: 233 donors, 20 cell types, frozen A=1355 I=1876.

### A1. The 0.195 is real

Stored P4 `within_site.shuffle_age_r2` = `0.1953817244152412` (JSON, not a rounding artifact). Independent original-ridge donor-level shuffle: within-site pooled R²=+0.177 (same magnitude as +0.195), site-stratified R²=−0.074. An intercept-only model — predict the training-fold mean age, no genes — already yields pooled within-site R²=+0.217 vs R²(age~site)=+0.264, and site-stratified R²=−0.032. LOSO shuffle remains largely negative, as in P4, because the intercept is the *other* bank's mean age.

### A2. Hypotheses

| id | fires | result |
|---|---|---|
| H1_donor_leakage | False | PASS: no donor in both sides of any fold |
| H2_shuffle_granularity | False | PASS: permute_age_within_site_obs remaps a permuted donor-level age onto every pseudobulk of that donor |
| H3_feature_selection_leakage | True | frozen pooled A R²=+0.435; nested=+0.345; P3 saved=+0.345 |
| H4_site_structure | True | intercept-only within-site pooled R²=+0.217 matches R²(age~site)=+0.264; site-stratified=-0.032 |
| H5_target_centering | False | PASS as a full-cohort-mean bug: Ridge(fit_intercept=True) on train; r2_mae uses the evaluation vector's mean |

- **H1_donor_leakage** (fires=False): asserted on loso and within_site folds
- **H2_shuffle_granularity** (fires=False): within-donor row shuffle is a no-op because age is donor-constant; per-pseudobulk shuffle is NOT what P4 ran
- **H3_feature_selection_leakage** (fires=True): P4c used frozen full-data sets; P3 nested inside folds
- **H4_site_structure** (fires=True): within-site k-fold never mixes banks in a fold, so the intercept equals the training bank's mean age; pooled-across-banks R² recovers the bank gap, which the within-site age permutation PRESERVES
- **H5_target_centering** (fires=False): the eval vector for within-site CV is two banks concatenated — see H4

**Fired:** H3_feature_selection_leakage, H4_site_structure. **Fix:** Primary fix: report site-stratified median-over-types R² for the within-site scheme (mean of the two banks' R² per cell type, then median over types). LOSO is left on the original pooled/test-site metric — intercept mismatch across banks is a real transfer failure, not a scoring artifact. P4c frozen vs nested is reported side-by-side; frozen remains the size-matched cross-test.

Donor age is constant across a donor's pseudobulks (True). Shuffling within a donor's rows is therefore a no-op. P4 already used donor-level within-site permutation (`permute_age_within_site_obs`). That permutation preserves each bank's age mean, which is exactly what the pooled within-site metric reads out via the intercept.

### A3. Reported vs corrected

Corrected within-site number = site-stratified median-over-types R² (mean of HBCC and MSSM R² per cell type, then median over types). Corrected LOSO number = original pooled median-over-types R² (intercept mismatch is a real transfer failure). Every predictive number below carries its permutation null.

**On whether the null is usable.** STOP A asked whether the *positive* leak remained (shuffle still above ~0.05). It does not: site-stratified A-gene shuffle mean is -0.232. The intercept-only site-stratified R² is −0.032, so the *metric* is sound. SVD-LOO ridge under permutation is more negative (−0.23) than original sklearn ridge (−0.07) because p ≫ n overfits; that is conservative (it cannot create a false age clock), not residual donor leak. FALSIFICATION.md fails soundness if shuffle R² > 0.10; that bar is now cleared. Treat observed R² as a contrast against this overfit-negative null, not as a number whose null sits in ±0.05. Original-ridge site-stratified shuffle (−0.074) *is* near zero.

#### within_site

| metric | reported_in_P4c | corrected | permutation_null | null_usable |
|---|---|---|---|---|
| P3 nested A-gene age R² (median types) | +0.345 | +0.123 | -0.244 | False |
| P3 unrestricted age R² (median types) | +0.483 | +0.310 | NA | NA |
| P4c R²_A frozen A-genes | +0.423 | +0.250 | -0.232 | False |
| P4c R²_I size-matched I | +0.559 | +0.413 | -0.262 | False |
| P4c R²_I full I set | +0.571 | +0.424 | -0.262 | False |
| P4c shuffle-age A-genes | +0.195 | -0.232 | -0.232 | False |
| P4c Acc_I | +0.999 | +0.999 | NA | NA |
| P4c Acc_A (identity from A) | +0.994 | +0.994 | NA | NA |
| nested A-gene age R² (this audit) | +0.345 | +0.125 | -0.244 | False |
| nested I-gene age R² (this audit) | NA | +0.462 | -0.262 | False |

#### loso

| metric | reported_in_P4c | corrected | permutation_null | null_usable |
|---|---|---|---|---|
| P3 nested A-gene age R² (median types) | -0.102 | -0.102 | -0.828 | False |
| P3 unrestricted age R² (median types) | +0.096 | +0.096 | NA | NA |
| P4c R²_A frozen A-genes | -0.019 | -0.016 | -1.005 | False |
| P4c R²_I size-matched I | +0.260 | +0.279 | -0.965 | False |
| P4c R²_I full I set | +0.267 | +0.289 | -0.965 | False |
| P4c shuffle-age A-genes | -0.913 | -1.005 | -1.005 | False |
| P4c Acc_I | +0.999 | +0.999 | NA | NA |
| P4c Acc_A (identity from A) | +0.991 | +0.991 | NA | NA |
| nested A-gene age R² (this audit) | -0.102 | -0.091 | -0.828 | False |
| nested I-gene age R² (this audit) | NA | +0.330 | -0.965 | False |

**STOP A:** False — corrected within-site shuffle null is usable

**P4c verdict change.** Reported within-site `soundness-fail` → corrected `falsified (not separable)` (R²_I >= 0.8 × R²_A (1.652); R²_A − R²_I = -0.163 < 0.10; Acc_A >= 0.9 × Acc_I (0.995)). Reported LOSO `falsified (not separable)` → corrected `falsified (not separable)` (R²_A − R²_I = -0.294 < 0.10; Acc_A >= 0.9 × Acc_I (0.992)). STRENGTHEN: soundness-fail cleared; separability still falsified on the cross-test. Ordering within-site: I_beats_A → I_beats_A; LOSO: I_beats_A → I_beats_A.

Feature-selection note: P4c used frozen full-data A/I sets. Nested re-derivation (declared quantile levels + P2 filters inside each training fold) is in the table as `nested A/I`. Frozen remains the size-matched cross-test, as in P4; nested is the leakage-safer clock.

## Part B — Gene-count control

Random sets drawn from all 25526 genes. A/I pools are the frozen P2 sets (A=1355, I=1876), not re-derived. 50 draws at size |A|, each with a paired donor-level shuffle. Size sweep [10, 50, 100, 500, 1000, 1500] with 10 observed draws and one permutation null per (scheme, family, k). Primary metric: site-stratified R² within-site, pooled R² LOSO.

### B1. Size-matched random sets

| scheme | n_A | n_draws | random_median | random_p05 | random_p95 | random_null_median | A_frozen | A_null | I_full | I_null | A_beats_random | I_beats_A | random_matches_A |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| within_site | 1355 | 50 | +0.254 | +0.209 | +0.298 | -0.234 | +0.250 | -0.155 | +0.424 | -0.230 | False | True | True |
| loso | 1355 | 50 | +0.025 | -0.025 | +0.086 | -1.026 | -0.016 | -0.971 | +0.289 | -1.265 | False | True | True |

### B2. Size sweep (median R²)

| scheme | family | k | k_effective | n_draws | r2_median | null_median | null_usable |
|---|---|---|---|---|---|---|---|
| within_site | A | 10 | 10 | 10 | -0.035 | -0.049 | True |
| within_site | I | 10 | 10 | 10 | -0.035 | -0.043 | True |
| within_site | random | 10 | 10 | 10 | -0.039 | -0.028 | True |
| within_site | A | 50 | 50 | 10 | -0.022 | -0.070 | True |
| within_site | I | 50 | 50 | 10 | -0.009 | -0.036 | True |
| within_site | random | 50 | 50 | 10 | -0.038 | -0.044 | True |
| within_site | A | 100 | 100 | 10 | -3.215 | -6.631 | True |
| within_site | I | 100 | 100 | 10 | -2.056 | -4.965 | True |
| within_site | random | 100 | 100 | 10 | -2.962 | -4.266 | True |
| within_site | A | 500 | 500 | 10 | +0.094 | -0.306 | True |
| within_site | I | 500 | 500 | 10 | +0.300 | -0.534 | True |
| within_site | random | 500 | 500 | 10 | +0.122 | -0.308 | True |
| within_site | A | 1000 | 1000 | 10 | +0.225 | -0.243 | True |
| within_site | I | 1000 | 1000 | 10 | +0.391 | -0.465 | True |
| within_site | random | 1000 | 1000 | 10 | +0.226 | -0.249 | True |
| within_site | A | 1500 | 1355 | 1 | +0.250 | -0.254 | True |
| within_site | I | 1500 | 1500 | 10 | +0.418 | -0.321 | True |
| within_site | random | 1500 | 1500 | 10 | +0.258 | -0.134 | True |
| loso | A | 10 | 10 | 10 | -0.653 | -0.831 | True |
| loso | I | 10 | 10 | 10 | -0.660 | -0.819 | True |
| loso | random | 10 | 10 | 10 | -0.730 | -0.784 | True |
| loso | A | 50 | 50 | 10 | -0.554 | -0.767 | True |
| loso | I | 50 | 50 | 10 | -0.387 | -0.992 | True |
| loso | random | 50 | 50 | 10 | -0.575 | -0.786 | True |
| loso | A | 100 | 100 | 10 | -0.949 | -0.868 | True |
| loso | I | 100 | 100 | 10 | -0.623 | -1.504 | True |
| loso | random | 100 | 100 | 10 | -0.845 | -1.783 | True |
| loso | A | 500 | 500 | 10 | -0.138 | -1.218 | True |
| loso | I | 500 | 500 | 10 | +0.194 | -1.488 | True |
| loso | random | 500 | 500 | 10 | -0.129 | -1.689 | True |
| loso | A | 1000 | 1000 | 10 | -0.026 | -1.071 | True |
| loso | I | 1000 | 1000 | 10 | +0.233 | -0.972 | True |
| loso | random | 1000 | 1000 | 10 | -0.002 | -1.407 | True |
| loso | A | 1500 | 1355 | 1 | -0.016 | -0.668 | True |
| loso | I | 1500 | 1500 | 10 | +0.279 | -1.329 | True |
| loso | random | 1500 | 1500 | 10 | +0.058 | -0.899 | True |

Figure: `results/geometry/figures/partB_size_sweep.png`.

### B3. Verdict

Not clean (a) and not clean (b). **A-genes match random at matched size** (within-site 0.250 vs random median 0.254, inside the random 5–95% band 0.209–0.298; LOSO both ~0). The A-bucket is not a privileged age axis — that is the arithmetic/methods result. **I-genes beat both A and random** (within-site 0.424; LOSO 0.289 vs random 0.025). Identity-associated genes carry extra age information that a random 1,355-gene set does not. Size sweep: A tracks random at every k; I pulls away from k ≥ 500. At k=100, SVD-LOO ridge is numerically unstable (hat diagonals near 1 when n≈p); do not interpret that cell. The fusion *gene-set* finding is therefore: genes cannot be partitioned into an age bucket, and age prediction lives in identity programmes.

- within-site B1: A=+0.250 random_med=+0.254 random_p95=+0.298 I=+0.424 null_random=-0.234
- LOSO B1: A=-0.016 random_med=+0.025 I=+0.289 null_random=-1.026

## Part C — Geometric reframe

All 25526 genes. Age direction: PLS-1 (covariance `X'y` within cell type), ridge coefficients as sensitivity. Identity subspace: SVD of the 20 type centroids (19 axes). Angle of a vector `w` to that subspace is $\arccos(\\|P_{\mathrm{id}} w\\|)$; projection R² is `||P_id w||²`. A random unit vector in R^25526 has E[proj R²]=0.0007.

### C3. Age direction vs identity subspace

| quantity | observed | permutation null mean | p |
|---|---:|---:|---:|
| median-over-types angle (PLS-1) | 66.21° | 70.89° | 0.030 |
| median-over-types proj R² | 0.163 | 0.108 | 0.030 |
| consensus age-dir angle | 58.83° | 61.88° | — |
| ridge median angle (sensitivity) | 83.32° | — | — |

Both CV schemes (directions fit on each training fold, z-scored on train, identity basis on train):

| scheme | median angle | proj R² | null angle | p_angle | n_perm |
|---|---:|---:|---:|---:|---:|
| loso | 74.74° | 0.070 | 79.89° | 0.048 | 20 |
| within_site | 74.59° | 0.071 | 78.92° | 0.095 | 20 |

**C3 reading:** largely orthogonal — the age direction is largely ORTHOGONAL to identity directions — the axes are separable in geometry even though genes cannot be partitioned. The excess vs the permutation null is real (p=0.03) but small (66° vs 71°; proj R² 0.163 vs 0.108). Ridge coefficients are closer still to orthogonal (83°). This is not 'substantially within' the identity subspace.

### C4. Shared vs type-specific age direction

Pairwise unsigned angle among the 20 PLS-1 age directions: median 78.3° (null 80.3°). PC1 of the 20 directions explains 0.405 of their variance (null 0.330, p=0.040). Per-type age directions are **largely type-specific** (pairwise ~78°, only ~12° from orthogonal; PC1 0.40 vs null 0.33 is a weak common component, not a shared axis). A rejuvenation intervention that moved one shared axis would not automatically move 20 type-specific 78° neighbours. CV schemes agree (pairwise ~81–83°).

Figures: `partC_angle_vs_null.png`, `partC_pairwise_angles.png`, `partC_per_type_angles.png`, `partC_shared_pc1.png`.

## Which claim

After correcting the site-pooled R² artifact, the gene-set cross-test is within-site **falsified (not separable)** and LOSO **falsified (not separable)** (STRENGTHEN: soundness-fail cleared; separability still falsified on the cross-test). **Genes cannot be partitioned into clean age vs identity buckets:** A-genes match size-matched random sets; I-genes beat both, so chronological age is carried by identity programmes rather than by a privileged A-set. That is a methods finding about the construct, plus a biological finding that identity genes are age-informative. **The axes are not geometrically fused.** The within-type age direction is largely orthogonal to the identity (centroid) subspace (66° vs null 71°; ridge 83°). Gene-set non-partitionability and geometric fusion are different claims; the evidence now supports the first, not the second. Per-type age directions are **largely type-specific** (pairwise ~78°, only ~12° from orthogonal; PC1 0.40 vs null 0.33 is a weak common component, not a shared axis).

These are different findings. Gene-set non-partitionability does not imply geometric fusion, and geometric orthogonality would reverse the project's conclusion even if A/I buckets fail.

## Limitations

1. Two brain banks only. LOSO is one df of transfer and is reported, not averaged away.
2. Unmeasured 6-plex hashing pools, PMI, RIN — as in FINDINGS_BRAIN_PHASE1.md.
3. Corrected within-site metric is site-stratified R², not a new gene set and not a retune of FALSIFICATION.md thresholds.
4. Seed `20260914` (`numpy.random.default_rng`).
5. Part B/C age models use SVD-LOO ridge (same α grid as P3/P4). A1 reproduced the shuffle with the original sklearn inner-CV ridge.
6. Size-sweep k=100 is numerically unstable (n≈p, LOO hat diagonals → 1). Dropped from interpretation.
7. Pairwise angles among per-type age directions (~78°) are close to the permutation null (~80°); do not read PC1=0.40 as a single shared aging axis.

## Files

| path | content |
|---|---|
| `src/geometry_common.py`, `geometry_partA.py`, `geometry_partB.py`, `geometry_partC.py` | code |
| `notebooks/geometry_partA.ipynb`, `partB`, `partC` | runnable from a clean checkout |
| `results/geometry/` | tables, json, logs, `figures/` |
| `FINDINGS_GEOMETRY.md` | this file |


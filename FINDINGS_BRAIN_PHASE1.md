# FINDINGS_BRAIN_PHASE1 — age/identity axes on DLPFC adult controls

**Status: P0, P1 (thresholds declared before scores), P2, P3, P4, P5.** Seed `20260914`. Cohort: 233 adult neurotypical DLPFC donors (ages 20–89), 3,323 pseudobulks, 20 cell types, 2 sites (Source = HBCC/MSSM). Reproduced by `notebooks/brain_phase1_p0.ipynb` … `p5`. Run date 2026-09-14.

**Do not modify** `FINDINGS.md`, `FINDINGS_RERUN.md`, `FINDINGS_STAGEB.md`, `FINDINGS_TISSUE.md`, `FINDINGS_BRAIN.md`, or `FALSIFICATION.md`. `FALSIFICATION.md` is pre-registration and governs Stage P4.

**P4c verdict:** NOT SEPARABLE. In both CV schemes the age model is as good or better on I-genes than on A-genes (within-site R²_A=+0.423 vs R²_I=+0.559; LOSO R²_A=-0.019 vs R²_I=+0.260). FALSIFICATION.md labels disagree (LOSO: falsified (not separable); within-site: soundness-fail) and are both reported; they agree on the comparison.

Blood remains closed (compositional). This file is the Phase 1 payoff run on brain.

LIMITATION: the source paper describes 6-sample hashing pools; they are not in obs and cannot be corrected for. Source (HBCC/MSSM) is a 2-level brain-bank factor (1 df), not equivalent to OneK1K's 75 multiplexed 10x pools. Barcodekey looks like Donor-{n}-{10x_barcode}-{gem_group} but is not used as a pool ID.

## P0 — Brain-specific confounders

Audited every obs column and `uns` on the cached Aging_Cohort h5ad (dataset ID `4442d412-91cb-4261-acca-8adf5fa04c11`, fetched in T1/R0 — not invented). Neuronal fraction and n_genes/n_counts use obs only (no 12 GB matrix). MALAT1/MT fractions use the existing 3,323-row pseudobulk.

| variable | in obs / recoverable? | notes |
|---|---|---|
| post-mortem interval (PMI) | **NO** | ABSENT — real limitation. Known post-mortem confound, often age-correlated. |
| RIN / RNA quality | **NO** | ABSENT — real limitation. |
| tissue pH | **NO** | ABSENT. |
| cause of death | **NO** | ABSENT. |
| SoupX / decontX / ambient scores | **NO** | Proxies computed: MALAT1 UMI fraction, MT- UMI fraction, fraction of nuclei with n_genes < 500. |
| 6-plex hashing pools | **NO** | LIMITATION: the source paper describes 6-sample hashing pools; they are not in obs and cannot be corrected for. Source (HBCC/MSSM) is a 2-level brain-bank factor (1 df), not equivalent to OneK1K's 75 multiplexed 10x pools. Barcodekey looks like Donor-{n}-{10x_barcode}-{gem_group} but is not used as a pool ID. `uns/batch_condition` = `Source` only. |
| neuronal nucleus fraction | **computed** | Donor-level fraction of nuclei with `class` in {EN, IN}. Dissociation-bias proxy. |

**A priori inclusion rule:** a donor-level expression-affecting covariate enters the P1 design if within-site Spearman |r| with age ≥ **0.2** in either site or in the site-residualized pool. log(mean UMI/nucleus) and log(mean genes/nucleus) are **always** in the design (specified), not selected on |r|.

Donors in the P0 table: **233**. Mean neuronal fraction: **0.385**.

Within-site Spearman(age, covariate):

| covariate | H | M | within_site_pooled |
|---|---|---|---|
| cc_umi_frac | 0.049 | 0.263 | 0.147 |
| frac_Immune | -0.292 | -0.032 | -0.159 |
| frac_Oligo | 0.087 | 0.276 | 0.185 |
| frac_low_n_genes | NA | NA | NA |
| frac_neuronal | 0.044 | -0.187 | -0.064 |
| log_mean_genes | -0.034 | -0.143 | -0.086 |
| log_mean_umi | 0.005 | -0.084 | -0.036 |
| malat1_frac | -0.191 | 0.031 | -0.079 |
| mean_n_counts | 0.005 | -0.084 | -0.035 |
| mean_n_genes | -0.034 | -0.143 | -0.085 |
| mt_frac | 0.000 | -0.008 | 0.005 |
| n_nuclei | -0.117 | -0.024 | -0.068 |

**P0 extras added to the P1 model (donor-level, merged as `p0_*`):** ['cc_umi_frac', 'frac_Oligo', 'frac_Immune']

`frac_low_n_genes` is NA because it is identically 0 on this QC'd upload (no nuclei with n_genes < 500 among adult analysis donors).

- `frac_neuronal`: not included: all within-site |r| < 0.2
- `frac_low_n_genes`: not included: all within-site |r| < 0.2
- `malat1_frac`: not included: all within-site |r| < 0.2
- `mt_frac`: not included: all within-site |r| < 0.2
- `mean_n_counts`: not included: all within-site |r| < 0.2
- `mean_n_genes`: not included: all within-site |r| < 0.2
- `n_nuclei`: not included: all within-site |r| < 0.2
- `cc_umi_frac`: within-site |Spearman r| >= 0.2 (M r=+0.263)
- `frac_Oligo`: within-site |Spearman r| >= 0.2 (M r=+0.276)
- `frac_Immune`: within-site |Spearman r| >= 0.2 (H r=-0.292)

Figures: `results/brain_phase1/figures/p0_covariates_vs_age.png`.
Obs columns audited: 33. `uns/batch_condition`=['Source'].

## P1 — Variance decomposition

LIMITATION: the source paper describes 6-sample hashing pools; they are not in obs and cannot be corrected for. Source (HBCC/MSSM) is a 2-level brain-bank factor (1 df), not equivalent to OneK1K's 75 multiplexed 10x pools. Barcodekey looks like Donor-{n}-{10x_barcode}-{gem_group} but is not used as a pool ID.

Model on the 3,323 pseudobulks:

```
expression ~ C(cell_type) + C(Source) + log(mean UMI/nucleus) + log(mean genes/nucleus)
             + [P0 extras] + age
```

`unique_age = R²(N+ct+age) − R²(N+ct)` with N = Source + within-type-centred depth + P0 extras. This is within-site, depth/P0-adjusted age variance. `unique_ct` analogously.

### Declared A/I thresholds (written BEFORE any score)

The blood lesson: `t_age` was an absolute unique_age cutoff calibrated on a distribution later shown to be ~65% pool-depth artifact; after correction that cutoff was absurdly strict (99.9th percentile). **Brain thresholds are QUANTILES of the corrected distribution**, with quantile *levels* declared a priori in `src/brain_phase1_common.py` as `{'t_age_q': 0.9, 'k_ct_q': 0.5, 't_ct_q': 0.9, 'k_age_q': 0.5}` — not chosen by looking at downstream R².

- **A-genes:** `unique_age ≥ Q90(unique_age)` AND `unique_ct ≤ Q50(unique_ct)`
- **I-genes:** `unique_ct ≥ Q90(unique_ct)` AND `unique_age ≤ Q50(unique_age)`
- **MIXED:** `unique_age ≥ Q90` AND `unique_ct ≥ Q90` — excluded from both scores

A sweep over high ∈ {0.80, 0.90, 0.95, 0.99} × low ∈ {0.25, 0.50, 0.75} is reported; the primary cell is the declared (0.90, 0.50) pair. **No cell is selected on score R².**

**Numeric cutoffs on this cohort's corrected distribution:**

| cutoff | quantile | value |
|---|---:|---:|
| t_age | Q0.90 unique_age | **0.002387** |
| k_ct | Q0.50 unique_ct | **0.405209** |
| t_ct | Q0.90 unique_ct | **0.815260** |
| k_age | Q0.50 unique_age | **0.000370** |

**Set sizes at the declared quantiles (pre-P2 filters):** A = **1578**, I = **1876**, MIXED = **37**, OTHER = **22035** (n genes = 25526).

These sizes were written into this file after P1 and before P3 (`results/brain_phase1/p1_DECLARED_BEFORE_SCORES.flag`).

### unique_age / unique_ct distribution vs blood

| | DLPFC corrected | blood within-pool (FINDINGS_RERUN) |
|---|---:|---:|
| unique_age median | 0.00037 | 0.00011 |
| unique_age 99th pct | 0.00666 (0.67%) | **0.0019 (0.19%)** |
| unique_age max | 0.03490 (3.49%) | **0.0126 (1.26%)** |
| unique_ct median | 0.4052 | 0.1248 |
| unique_ct max | 0.9774 | 0.9610 |

P0 extras in this partition: `['p0_cc_umi_frac', 'p0_frac_Oligo', 'p0_frac_Immune']`.

Within-site donor-age permutation null (orientation only; **not** used to pick thresholds):

| t_age | observed | null_mean | null_p95 | empirical_FDR |
|---|---|---|---|---|
| 0.0010 | 7128.0000 | 2736.8900 | 4409.1000 | 0.3840 |
| 0.0020 | 3357.0000 | 869.3450 | 1626.4500 | 0.2590 |
| 0.0050 | 584.0000 | 94.8050 | 173.2000 | 0.1623 |
| 0.0100 | 86.0000 | 15.1800 | 31.0500 | 0.1765 |
| 0.0200 | 11.0000 | 2.0050 | 6.0000 | 0.1823 |

Sweep: n A-genes (rows = high quantile for unique_age / unique_ct; columns = low quantile). Primary cell marked in the csv `is_primary`.

```
low_q   0.25  0.50  0.75
high_q                  
0.80    1565  3076  4482
0.90     792  1578  2309
0.95     395   790  1175
0.99      85   159   242
```

Sweep: n I-genes:

```
low_q   0.25  0.50  0.75
high_q                  
0.80    1858  3407  4532
0.90    1062  1876  2365
0.95     594   997  1214
0.99     127   205   249
```

Tables: `p1_variance_decomposition.csv`, `p1_thresholds.json`, `p1_threshold_sweep_set_sizes.csv`.
Figures: `p1_unique_age_vs_unique_ct.png`, `p1_unique_age_hist_vs_blood.png`.

## P2 — Confound filters on A-genes

LIMITATION: the source paper describes 6-sample hashing pools; they are not in obs and cannot be corrected for. Source (HBCC/MSSM) is a 2-level brain-bank factor (1 df), not equivalent to OneK1K's 75 multiplexed 10x pools. Barcodekey looks like Donor-{n}-{10x_barcode}-{gem_group} but is not used as a pool ID.

Filters (same lists/thresholds as blood 1b, plus P0-correlated genes). Sequential survival of the declared A set:

| filter | flagged_marginal | removed_incremental | surviving |
|---|---|---|---|
| cell-cycle (Tirosh S/G2M + proliferation list) | 5 | 5 | 1573 |
| cell-cycle (data: |r|>=0.3 with cycling covariates) | 94 | 94 | 1479 |
| ribosomal (RPL/RPS/MRP) | 35 | 26 | 1453 |
| mitochondrial (MT-) | 2 | 2 | 1451 |
| sex-linked (chrX/chrY) | 97 | 89 | 1362 |
| sex-DE (|d|>=0.5, p<1e-6, within cell type) | 6 | 6 | 1356 |
| depth-correlated (|r|>=0.3 with mean UMI/genes per nucleus) | 6 | 1 | 1355 |
| P0-covariate-correlated (|r|>=0.3, within type) | 11 | 0 | 1355 |

A-genes in: **1578**. Surviving: **1355**. I-genes (unfiltered by A-filters; MIXED still excluded): **1876**.

STOP P2 (< 30 surviving A-genes): **False**.

**Cell-cycle note (DLPFC is largely post-mitotic):** 99/1578 declared A-genes are cell-cycle flagged (list=5, data=94). DLPFC is largely post-mitotic; OPCs are the main cycling population. Small CC hit, as expected for post-mitotic cortex.

## P3 — The two scores

LIMITATION: the source paper describes 6-sample hashing pools; they are not in obs and cannot be corrected for. Source (HBCC/MSSM) is a 2-level brain-bank factor (1 df), not equivalent to OneK1K's 75 multiplexed 10x pools. Barcodekey looks like Donor-{n}-{10x_barcode}-{gem_group} but is not used as a pool ID.

Two CV schemes, both with **no donor in both folds** (asserted):

1. **LOSO** — leave-one-site-out (2 folds: train HBCC/test MSSM and vice versa). Batch-grouped, weak (n_batches=2).
2. **within-site donor k-fold** — 5-fold grouped by donor, nested *inside each site*.

Gene selection (P1 quantiles + P2 filters) is **re-done inside each outer training fold** (FALSIFICATION.md). Reported performance is nested-CV. Frozen full-data gene sets (after P2) are saved for P5 and for listing, and are not the CV features.

Age model: ridge on A-genes → donor age, within cell type. Identity: multinomial logistic on I-genes → cell type. Unrestricted age model: kernel ridge on all genes, within cell type (R3-like ceiling).

### loso

- Median-over-types age R² (A-genes): **-0.102** (MAE 15.44)
- Median-over-types unrestricted R²: **+0.096**
- Joint (all types) A-gene age R²: **+0.032**; unrestricted joint: **+0.217**
- Identity accuracy (I-genes): **0.999**  macro-F1 **0.996**  chance **0.070**

### within_site

- Median-over-types age R² (A-genes): **+0.345** (MAE 11.08)
- Median-over-types unrestricted R²: **+0.483**
- Joint (all types) A-gene age R²: **+0.465**; unrestricted joint: **+0.588**
- Identity accuracy (I-genes): **0.999**  macro-F1 **0.998**  chance **0.070**

Per-cell-type tables: `p3_age_per_type.csv`, `p3_identity.csv`.

STOP P3: **False**. 

Interpretation vs R3 ceiling (~0.44 joint within-site expression): A-gene construct captures the within-type age signal (A 0.345 vs unrestricted 0.483 median-over-types, within-site). Proceed.

## P4 — Negative controls and the cross-test (the verdict)

Governed by `FALSIFICATION.md` (written before any 1c/1d model; not revised). Size-matched I-gene age models: 20 draws, median reported.

### loso

- P4a shuffled ages, A-gene R² (median types): **-0.913**
- P4b shuffled types, identity acc: **0.058**
- R²_A (proper A-genes) = **-0.019**; R²_I (size-matched I-genes, median of 20) = **+0.260**; R²_I (full I set) = **+0.267**
- Acc_I = **0.999**; Acc_A (identity from A-genes) = **0.991**
- Δ R²_A − R²_I(matched) = **-0.279**; R²_I / R²_A = **-13.492**
- Acc_A / Acc_I = **0.992**
- scheme verdict: **falsified (not separable)**

### within_site

- P4a shuffled ages, A-gene R² (median types): **+0.195**
- P4b shuffled types, identity acc: **0.054**
- R²_A (proper A-genes) = **+0.423**; R²_I (size-matched I-genes, median of 20) = **+0.559**; R²_I (full I set) = **+0.571**
- Acc_I = **0.999**; Acc_A (identity from A-genes) = **0.994**
- Δ R²_A − R²_I(matched) = **-0.135**; R²_I / R²_A = **1.319**
- Acc_A / Acc_I = **0.995**
- scheme verdict: **soundness-fail**

**THE DECISION (plain):** DISAGREEMENT between CV schemes

LOSO verdict = falsified (not separable); within-site verdict = soundness-fail. Hard rule: report the disagreement rather than picking the better number. LOSO: R²_A − R²_I = -0.279 < 0.10; Acc_A >= 0.9 × Acc_I (0.992). Within-site: shuffle-age R²=0.195 > 0.10 (pipeline soundness fail); R²_I >= 0.8 × R²_A (1.319); R²_A − R²_I = -0.135 < 0.10; Acc_A >= 0.9 × Acc_I (0.995).

Plain scientific result (both CV schemes): chronological age is predicted as well or BETTER from I-genes than from A-genes (within-site median-over-types R2_A=+0.423 vs size-matched R2_I=+0.559; LOSO R2_A=-0.019 vs R2_I=+0.260). Identity accuracy from A-genes is 0.99 vs 1.00 from I-genes. Under FALSIFICATION.md this falsifies separability: age is not a distinct axis carried by high-unique_age/low-unique_ct genes. This is the failure mode the pre-registration said to report as a result, not a gene-set bug (age acting through identity programmes; and low per-gene unique_ct does not yield no identity information in a 1000-gene set when unique_ct can be as high as the median 0.40). Shuffles: identity collapsed to chance in both schemes. Age shuffle collapsed in LOSO (R2=-0.91) but within-site OOF R2=+0.195 exceeds the 0.10 soundness line (high-p ridge on ~100 donors/type; not used to rescue separability). Combined with blood (2 surviving A-genes after pool correction), the per-gene unique_age vs unique_ct construction does not yield separable age/identity axes in PBMC or in DLPFC under the declared rules.

Where LOSO and within-site disagree on the *FALSIFICATION.md label*, both labels are reported. They **agree** on the comparison the decision rule is about: in both schemes the age model is as good or better on I-genes as on A-genes, and identity accuracy on A-genes is essentially the same as on I-genes. That is the cross-test result. Thresholds were not retuned.

## P5 — Replication on the frozen sets (white-matter atlas)

Dataset `c05e6940-729c-47bd-a2a6-6ce3730c4919` (cached). **Nothing was re-derived.** Frozen DLPFC A/I gene IDs (Ensembl, version-stripped) were mapped onto the WM pseudobulk.

- A-genes mapped: **893** / 1355
- I-genes mapped: **1802** / 1876
- A-gene age R² (median over matched types, SequencingPool-grouped CV): **-0.518**
- size-matched I-gene age R²: **-0.425**
- cross-test holds? **False**

Frozen A-genes do **not** retain a clear age association on WM (median-over-types R²=-0.518, n_donors=20). Cross-test cannot confirm DLPFC.

**Does not replicate** as a clean copy of the DLPFC verdict. Recommend the retina 104-donor snRNA atlas (`d6505c89-c43d-4c28-8c4f-7351a5fd5528`) as the second replication target.

## Limitations (every stage that depends on them)

1. **Two-batch CV.** Source has 2 levels (HBCC/MSSM). Leave-one-site-out is the only batch-grouped scheme and is weak: the two folds are the two banks, MSSM is older (R2(age~Source)=0.264), and repeats are not independent. Donor-grouped k-fold *within each site* is reported alongside it. Where they disagree, both numbers stand. This is not OneK1K's 75-pool CV.
2. **Unmeasured 6-plex hashing pools.** LIMITATION: the source paper describes 6-sample hashing pools; they are not in obs and cannot be corrected for. Source (HBCC/MSSM) is a 2-level brain-bank factor (1 df), not equivalent to OneK1K's 75 multiplexed 10x pools. Barcodekey looks like Donor-{n}-{10x_barcode}-{gem_group} but is not used as a pool ID. Pool-structured age would inflate `unique_age` the way OneK1K pools did; we cannot test that on this upload.
3. **PMI and RIN are not in obs.** Post-mortem interval and RNA integrity are the standard brain-snRNA confounds and are often age-correlated. They were not uploaded. Neuronal fraction, MALAT1, MT fraction, and low-n_genes debris are the recoverable proxies; they are not substitutes for PMI/RIN.
4. **Single-nucleus, 10x 3′ v3.** Metrics are comparable to blood scRNA; absolute sensitivity is not.
5. **Gene selection nested in CV** uses the same declared *quantile levels* on each training fold's own unique_age / unique_ct distribution. Frozen full-data sets are for P5 and for listing genes, not for the reported scores.
6. Seed `20260914` for every RNG (`numpy.random.default_rng`).
7. **P4 within-site age shuffle** left median-over-types R² = +0.195 (above the 0.10 soundness line). LOSO shuffle collapsed (R² = −0.91). High-p ridge (1,355 genes) on ~100 donors per type can leave noisy positive OOF R² after a label permutation; this is reported, not used to claim or rescue separability. Identity shuffles hit chance in both schemes.
8. P4 cross-test uses the **frozen** full-data A/I sets so A vs I are size-and-selection comparable. P3 nested A-gene R² (within-site 0.345) is the leakage-safer estimate of the A-gene clock; frozen P4 R²_A is 0.423. The I > A ordering holds for both.

## Files

| path | content |
|---|---|
| `src/brain_phase1_common.py`, `brain_phase1_p0.py` … `p5.py`, `brain_phase1_models.py`, `brain_phase1_findings.py`, `brain_phase1_run.py`, `brain_phase1_make_notebooks.py` | code |
| `notebooks/brain_phase1_p0.ipynb` … `p5` | runnable per stage |
| `results/brain_phase1/` | tables, json, logs, `figures/` |
| `FINDINGS_BRAIN_PHASE1.md` | this file |
| `FALSIFICATION.md` | pre-registration (unchanged) |

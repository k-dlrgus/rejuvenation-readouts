# FINDINGS_MD3 — test the ruler where their claim lives; compare instruments against donor age

**Status:** Task 1: `instruments_disagree_on_claim_population`. Task 2: `ruler_stronger_age_in_ge2`. Task 3: report-only MD co-variation tabulated; no gate. Fired: Task 1 `instruments_disagree_on_claim_population`; Task 2 `ruler_stronger_age_in_ge2`. Seed `20260914`. boot `20260918`. n_perm=200. n_boot=200. n_random=200. Frozen ruler `frozen_ruler_ridge_raw.npz` exists=True.

Does not modify any existing `FINDINGS_*.md` or `FALSIFICATION.md`. No GSE325735. Uncalibrated R² is reported and is never a gate. Nothing averaged across donors, clusters, or instruments. Refit nothing. d7 is not the primary Stage 2 endpoint. ≥10-cell threshold is a sensitivity check, not primary. Paired Δρ is the age-instrument difference; overlapping separate CIs are not.

Reproduced by `src/md3_run.py`. Frozen ruler: `<repo>\results\fibro\frozen_ruler_ridge_raw.npz` exists=True.

Flag: `PREREG_TASK1.flag` exists=True. `PREREG_TASK2.flag` exists=True. `PREREG_TASK3.flag` exists=True. `DECLARED_BEFORE_SCORES.flag` exists=True.

## Why Table S3 age-up / age-down lists are exact, not an approximation

FINDINGS_MD2.md did not run the authors' aging signature because GSE113957 is FPKM
and the STAR Methods DESeq2 old-versus-young split (adjusted p < 0.05, log2FC > 0.5,
old/young from PCA separation, DESeq2 v1.40.2) could not be rebuilt from integer
counts. That refusal still stands: this prompt does not rebuild the DE analysis.

mmc3.xlsx Table S3 (sheet `Aging_signatures`, columns `Age up` and `Age down`) is
the authors' own published output of that DESeq2 analysis. Scoring those lists is
using the published result, not reconstructing the test that produced it. That is
therefore permitted here, where approximating the DE from FPKM was not.

Mismatch versus a rebuilt DE, if one were later released, would be reported; it is
not a reason to withhold the published lists.


- mmc3 Age up n=1533 id_type=symbol Age down n=2007 id_type=symbol
- mmc3 equals md2 genesets.json copies: up=True down=True
- used_as_deseq2_rebuild=False used_as_published_lists=True

## Failures recorded

- **addmodulescore:** GSE113957: cut_number did not yield 24 bins (got 22). Not dropping bins.

## Task 1 pre-registration (verbatim, written before any MD3 instrument score)

Written **before any MD3 age / MD / TGF-β / Table S3 / pluripotency score**, 2026-09-19. No threshold in this block is re-tuned after numbers exist.

Using the src/md2_*.py Louvain clustering already computed (19 clusters per donor; not re-run, not re-tuned). Cell-state labels are the md2 assignment: argmax of mean AddModuleScore of mmc3 Reprog_cell_state_signatures columns Fibroblast / PartialReprog / EarlyPluripotency / Pluripotency / NonReprog. Labels are read from results/md2/t2_cluster_labels_*.csv. Not re-labelled.

Three aggregation levels, each reported separately and never averaged together:
- (a) pooled by cell state — all PartialReprog clusters pooled into one pseudobulk per timepoint, likewise each other state. This is the level their Figure 3G operates at and is the **primary cell**. No MIN_CELLS cutoff is applied to the pool (pooling is the reason this prompt exists); n cells at each timepoint is reported. A state with n=0 at an endpoint is not pairable.
- (b) per cluster at ≥20 cells (the FINDINGS_MD2.md threshold, reproduced unchanged).
- (c) per cluster at ≥10 cells, reported as a relaxed-threshold sensitivity check and labelled as such. Not primary.

Per donor, never averaged: frozen ruler score, MD score, TGF-β score, Table S3 age-up score, Table S3 age-down score, and pluripotency−fibroblast score (frozen FINDINGS_FIBRO.md lists, same differentiation_score as md2), at d0/d3/d7/d10, with d0→d7 and d0→d10 endpoint changes.

Nulls (not interchangeable):
- Frozen ruler: 200 permuted-ruler-weight directions, seed 20260914. Empirical p = (n_random ≥ real + 1) / (n_random + 1). Endpoint decline = score(day 0) − score(day T). Positive = younger on the ruler.
- Every gene-list score (MD, TGF-β, Age up, Age down): 200 size-matched, expression-bin-matched random gene sets, AddModuleScore nbin=24 ctrl=100, seed 20260914 with a documented per-list offset so MD uses the same stream as md2. Not permuted ruler weights.
- Pluripotency−fibroblast is the frozen Stage 2 mean-z difference on the same TMM/z as the ruler, not AddModuleScore. Its null is 200 size-matched random gene sets in frozen-ruler gene space on that Z (not AMS bins). The difference is stated.

Donors GM00731 and GM23815 reported separately. d7 is not the Stage 2 endpoint. Frozen ruler results/fibro/frozen_ruler_ridge_raw.npz used as-is. Refit nothing. Louvain resolution, AddModuleScore parameters, and cell filters are inherited from md2 and are not re-tuned.

**Pre-registered reading (only the outcome that fired):**
- Ruler declines in the pooled PartialReprog state with p ≤ 0.05 → the FINDINGS_MD2.md `negative_holds` reading was limited by cell counts, not by the ruler. Say so plainly and quote the sentence superseded.
- Ruler does not decline in the pooled PartialReprog state while MD does, both with their own nulls → the instruments disagree on exactly the population the published claim concerns. This is the sharpest form of the disagreement; report both numbers side by side and declare no winner.
- Neither moves in the pooled state → our reproduction does not recover their Figure 3G at the pooled level; report the reproduction failure as the finding, not as a refutation.
- Pooled and per-cluster (≥20) levels disagree → report both, say which is which, do not pick.

“Declines” / “moves” = endpoint decline d0→d7 or d0→d10 with p ≤ 0.05 on that instrument's own null, on the aged donor GM00731 pooled PartialReprog state for (a), or on at least one aged-donor Louvain cluster for (b). GM23815 is tabulated as a contrast, never pooled, never the reading.


## Task 1 — cluster labels (md2 assignment, reused)

Assignment basis: argmax of mean AddModuleScore of mmc3 `Reprog_cell_state_signatures` (Fibroblast, PartialReprog, EarlyPluripotency, Pluripotency, NonReprog). Labels from `results/md2/t2_cluster_labels_*.csv`. Not re-labelled. Prompt's EarlyPluri is this EarlyPluripotency column.

| cell_line | cluster | label | n_cells_all | n_d0 | n_d3 | n_d7 | n_d10 | mean_Fibroblast | mean_PartialReprog | mean_EarlyPluripotency | mean_Pluripotency | mean_NonReprog |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | 0 | Fibroblast | 4012 | 3998 | 4 | 9 | 1 | +0.815 | +0.105 | +0.372 | -0.143 | +0.197 |
| GM00731 | 1 | NonReprog | 2322 | 0 | 0 | 9 | 2313 | +0.451 | +0.121 | +0.266 | -0.154 | +0.898 |
| GM00731 | 2 | EarlyPluripotency | 2096 | 0 | 0 | 2092 | 4 | +0.296 | +0.266 | +0.840 | +0.122 | -0.018 |
| GM00731 | 3 | PartialReprog | 2039 | 1 | 1589 | 444 | 5 | +0.362 | +0.508 | +0.295 | -0.153 | +0.261 |
| GM00731 | 4 | NonReprog | 2031 | 9 | 259 | 1732 | 31 | +0.500 | +0.157 | +0.342 | -0.122 | +0.578 |
| GM00731 | 5 | PartialReprog | 1594 | 2 | 1523 | 48 | 21 | +0.187 | +0.618 | +0.407 | -0.074 | -0.030 |
| GM00731 | 6 | Fibroblast | 1030 | 0 | 0 | 1 | 1029 | +0.434 | +0.259 | +0.379 | -0.134 | +0.424 |
| GM00731 | 7 | EarlyPluripotency | 1014 | 0 | 108 | 884 | 22 | +0.360 | +0.300 | +0.630 | +0.045 | +0.283 |
| GM00731 | 8 | PartialReprog | 986 | 0 | 973 | 13 | 0 | +0.130 | +0.665 | +0.487 | -0.009 | -0.104 |
| GM00731 | 9 | EarlyPluripotency | 961 | 0 | 0 | 909 | 52 | +0.106 | +0.400 | +0.845 | +0.474 | -0.116 |
| GM00731 | 10 | NonReprog | 863 | 56 | 89 | 191 | 527 | +0.299 | -0.115 | +0.108 | -0.039 | +0.512 |
| GM00731 | 11 | PartialReprog | 626 | 0 | 557 | 64 | 5 | +0.163 | +0.563 | +0.390 | -0.080 | -0.045 |
| GM00731 | 12 | Fibroblast | 570 | 569 | 0 | 1 | 0 | +0.872 | +0.143 | +0.426 | -0.147 | +0.158 |
| GM00731 | 13 | PartialReprog | 445 | 0 | 442 | 3 | 0 | +0.116 | +0.415 | +0.252 | -0.099 | -0.054 |
| GM00731 | 14 | Pluripotency | 439 | 0 | 5 | 354 | 80 | -0.001 | -0.052 | +0.591 | +0.659 | -0.069 |
| GM00731 | 15 | Fibroblast | 376 | 372 | 1 | 3 | 0 | +0.653 | +0.353 | +0.288 | -0.170 | +0.121 |
| GM00731 | 16 | EarlyPluripotency | 343 | 2 | 302 | 32 | 7 | +0.075 | +0.120 | +0.208 | +0.006 | +0.019 |
| GM00731 | 17 | NonReprog | 243 | 0 | 0 | 0 | 243 | +0.334 | +0.370 | +0.253 | -0.170 | +0.673 |
| GM00731 | 18 | PartialReprog | 142 | 0 | 136 | 6 | 0 | +0.213 | +0.526 | +0.322 | -0.118 | +0.044 |
| GM23815 | 0 | Fibroblast | 6489 | 6454 | 19 | 14 | 2 | +0.730 | +0.103 | +0.381 | -0.198 | +0.314 |
| GM23815 | 1 | EarlyPluripotency | 3251 | 0 | 3 | 3240 | 8 | +0.341 | +0.241 | +0.730 | -0.008 | +0.089 |
| GM23815 | 2 | EarlyPluripotency | 2731 | 0 | 0 | 1 | 2730 | +0.268 | +0.423 | +0.685 | +0.013 | +0.051 |
| GM23815 | 3 | EarlyPluripotency | 2570 | 0 | 3 | 2537 | 30 | +0.183 | +0.344 | +0.823 | +0.304 | -0.043 |
| GM23815 | 4 | Fibroblast | 2242 | 2 | 21 | 2210 | 9 | +0.470 | +0.152 | +0.445 | -0.139 | +0.354 |
| GM23815 | 5 | PartialReprog | 2171 | 9 | 2039 | 118 | 5 | +0.327 | +0.554 | +0.310 | -0.222 | +0.199 |
| GM23815 | 6 | EarlyPluripotency | 2149 | 0 | 0 | 1819 | 330 | +0.091 | +0.495 | +0.694 | +0.334 | -0.088 |
| GM23815 | 7 | PartialReprog | 2069 | 0 | 2039 | 20 | 10 | +0.146 | +0.638 | +0.526 | -0.021 | -0.049 |
| GM23815 | 8 | NonReprog | 1812 | 0 | 10 | 8 | 1794 | +0.484 | +0.123 | +0.349 | -0.178 | +0.742 |
| GM23815 | 9 | Fibroblast | 1365 | 0 | 0 | 3 | 1362 | +0.439 | +0.286 | +0.432 | -0.188 | +0.309 |
| GM23815 | 10 | Pluripotency | 1249 | 0 | 91 | 977 | 181 | -0.023 | -0.032 | +0.387 | +0.521 | -0.057 |
| GM23815 | 11 | PartialReprog | 1140 | 0 | 1065 | 65 | 10 | +0.140 | +0.500 | +0.351 | -0.101 | -0.020 |
| GM23815 | 12 | PartialReprog | 1101 | 0 | 970 | 106 | 25 | +0.226 | +0.528 | +0.407 | -0.162 | +0.051 |
| GM23815 | 13 | EarlyPluripotency | 898 | 56 | 246 | 164 | 432 | +0.160 | -0.017 | +0.281 | +0.057 | +0.163 |
| GM23815 | 14 | Fibroblast | 896 | 886 | 5 | 5 | 0 | +0.604 | +0.337 | +0.304 | -0.225 | +0.256 |
| GM23815 | 15 | EarlyPluripotency | 761 | 0 | 0 | 0 | 761 | +0.366 | +0.395 | +0.628 | -0.048 | +0.263 |
| GM23815 | 16 | Fibroblast | 363 | 361 | 1 | 1 | 0 | +0.741 | +0.099 | +0.383 | -0.182 | +0.320 |
| GM23815 | 17 | EarlyPluripotency | 307 | 0 | 64 | 122 | 121 | +0.285 | +0.331 | +0.567 | -0.007 | +0.287 |
| GM23815 | 18 | EarlyPluripotency | 294 | 0 | 0 | 289 | 5 | +0.177 | +0.253 | +0.797 | +0.283 | -0.022 |

Pooled n cells by state (sum of clusters with that label):

| cell_line | label | n_clusters | n_d0 | n_d3 | n_d7 | n_d10 |
|---|---|---|---|---|---|---|
| GM00731 | Fibroblast | 4 | 4939 | 5 | 14 | 1030 |
| GM00731 | PartialReprog | 6 | 3 | 5220 | 578 | 31 |
| GM00731 | EarlyPluripotency | 4 | 2 | 410 | 3917 | 85 |
| GM00731 | Pluripotency | 1 | 0 | 5 | 354 | 80 |
| GM00731 | NonReprog | 4 | 65 | 348 | 1932 | 3114 |
| GM23815 | Fibroblast | 5 | 7703 | 46 | 2233 | 1373 |
| GM23815 | PartialReprog | 4 | 9 | 6113 | 309 | 50 |
| GM23815 | EarlyPluripotency | 8 | 56 | 316 | 8172 | 4417 |
| GM23815 | Pluripotency | 1 | 0 | 91 | 977 | 181 |
| GM23815 | NonReprog | 1 | 0 | 10 | 8 | 1794 |

## Task 1 — three aggregation levels (not averaged)

### Primary cell excerpt — aged GM00731 pooled PartialReprog (not a substitute for the tables below)

| cell_line | label | day | n_cells | age_score | md_score | tgfb_score | age_up | age_down | pluri_primary |
|---|---|---|---|---|---|---|---|---|---|
| GM00731 | PartialReprog | 0 | 3 | +1.397 | +0.291 | -0.028 | +0.009 | +0.032 | -1.044 |
| GM00731 | PartialReprog | 3 | 5220 | +2.337 | +0.135 | -0.047 | -0.018 | +0.081 | +8.231 |
| GM00731 | PartialReprog | 7 | 578 | +4.737 | +0.284 | -0.013 | +0.011 | +0.077 | +6.796 |
| GM00731 | PartialReprog | 10 | 31 | +4.933 | +0.297 | -0.013 | +0.007 | +0.038 | +4.272 |

Aged pooled PartialReprog ruler and MD nulls:

| instrument | endpoint | n_cells_d0 | n_cells_end | score_d0 | score_end | decline | p | n_random_ge_real | pass_p | ok | reason | null_kind |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| frozen_ruler | d0→d7 | +3.000 | +578.000 | +1.397 | +4.737 | -3.339 | +0.403 | +80.000 | False | True | NA | permuted_ruler_weights |
| frozen_ruler | d0→d10 | +3.000 | +31.000 | +1.397 | +4.933 | -3.535 | +0.483 | +96.000 | False | True | NA | permuted_ruler_weights |
| md_score | d0→d7 | +3.000 | +578.000 | +0.291 | +0.284 | +0.007 | +0.015 | +2.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| md_score | d0→d10 | +3.000 | +31.000 | +0.291 | +0.297 | -0.006 | +0.886 | +177.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |

### (a) pooled by cell state — primary cell (Figure 3G level)

| cell_line | label | day | n_cells | age_score | md_score | tgfb_score | age_up | age_down | pluri_primary | below_min_cells |
|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | EarlyPluripotency | 0 | 2 | -4.307 | +0.249 | -0.021 | +0.064 | -0.001 | +9.508 | False |
| GM00731 | Fibroblast | 0 | 4939 | +5.253 | +0.506 | +0.042 | +0.055 | +0.026 | -1.183 | False |
| GM00731 | NonReprog | 0 | 65 | -0.077 | +0.350 | +0.086 | +0.097 | -0.040 | +0.141 | False |
| GM00731 | PartialReprog | 0 | 3 | +1.397 | +0.291 | -0.028 | +0.009 | +0.032 | -1.044 | False |
| GM00731 | EarlyPluripotency | 3 | 410 | +1.927 | +0.050 | -0.024 | +0.007 | +0.015 | +9.800 | False |
| GM00731 | Fibroblast | 3 | 5 | +10.582 | +0.378 | +0.023 | +0.001 | +0.050 | +1.872 | False |
| GM00731 | NonReprog | 3 | 348 | +8.122 | +0.476 | +0.075 | +0.028 | +0.026 | +4.549 | False |
| GM00731 | PartialReprog | 3 | 5220 | +2.337 | +0.135 | -0.047 | -0.018 | +0.081 | +8.231 | False |
| GM00731 | Pluripotency | 3 | 5 | -18.978 | +0.042 | +0.006 | +0.027 | +0.005 | +6.791 | False |
| GM00731 | EarlyPluripotency | 7 | 3917 | -0.739 | +0.152 | -0.026 | +0.033 | +0.034 | +9.059 | False |
| GM00731 | Fibroblast | 7 | 14 | +8.596 | +0.425 | +0.064 | +0.037 | +0.025 | -0.495 | False |
| GM00731 | NonReprog | 7 | 1932 | +7.015 | +0.445 | +0.048 | +0.046 | +0.020 | +3.203 | False |
| GM00731 | PartialReprog | 7 | 578 | +4.737 | +0.284 | -0.013 | +0.011 | +0.077 | +6.796 | False |
| GM00731 | Pluripotency | 7 | 354 | -8.138 | +0.011 | +0.005 | +0.056 | -0.022 | +11.677 | False |
| GM00731 | EarlyPluripotency | 10 | 85 | -0.190 | +0.079 | -0.086 | +0.016 | +0.051 | +10.334 | False |
| GM00731 | Fibroblast | 10 | 1030 | +7.873 | +0.434 | -0.009 | +0.020 | +0.049 | +2.789 | False |
| GM00731 | NonReprog | 10 | 3114 | +9.474 | +0.434 | +0.001 | +0.055 | +0.023 | +1.688 | False |
| GM00731 | PartialReprog | 10 | 31 | +4.933 | +0.297 | -0.013 | +0.007 | +0.038 | +4.272 | False |
| GM00731 | Pluripotency | 10 | 80 | -5.242 | +0.069 | +0.019 | +0.051 | -0.019 | +10.335 | False |
| GM23815 | EarlyPluripotency | 0 | 56 | -5.297 | +0.379 | +0.077 | +0.083 | -0.038 | +0.733 | False |
| GM23815 | Fibroblast | 0 | 7703 | +0.241 | +0.479 | +0.029 | +0.039 | +0.029 | -0.795 | False |
| GM23815 | PartialReprog | 0 | 9 | +2.361 | +0.316 | +0.006 | -0.010 | +0.093 | -2.091 | False |
| GM23815 | EarlyPluripotency | 3 | 316 | -1.037 | +0.059 | +0.017 | +0.005 | +0.021 | +9.668 | False |
| GM23815 | Fibroblast | 3 | 46 | +5.947 | +0.385 | +0.081 | -0.000 | +0.049 | +1.930 | False |
| GM23815 | NonReprog | 3 | 10 | -0.790 | +0.380 | +0.051 | +0.020 | +0.055 | +4.556 | False |
| GM23815 | PartialReprog | 3 | 6113 | -2.952 | +0.085 | -0.036 | -0.021 | +0.092 | +8.918 | False |
| GM23815 | Pluripotency | 3 | 91 | -13.423 | +0.028 | -0.012 | +0.021 | +0.007 | +10.499 | False |
| GM23815 | EarlyPluripotency | 7 | 8172 | -8.203 | +0.077 | -0.058 | +0.017 | +0.044 | +11.120 | False |
| GM23815 | Fibroblast | 7 | 2233 | -1.817 | +0.410 | -0.015 | +0.033 | +0.025 | +4.290 | False |
| GM23815 | NonReprog | 7 | 8 | -5.052 | +0.407 | -0.004 | +0.060 | +0.028 | -2.206 | False |
| GM23815 | PartialReprog | 7 | 309 | -2.809 | +0.174 | -0.058 | +0.009 | +0.056 | +7.891 | False |
| GM23815 | Pluripotency | 7 | 977 | -12.384 | -0.059 | +0.006 | +0.027 | -0.004 | +12.645 | False |
| GM23815 | EarlyPluripotency | 10 | 4417 | -3.487 | +0.166 | -0.023 | +0.009 | +0.050 | +8.721 | False |
| GM23815 | Fibroblast | 10 | 1373 | +0.181 | +0.379 | +0.037 | +0.012 | +0.047 | +2.943 | False |
| GM23815 | NonReprog | 10 | 1794 | +2.866 | +0.423 | +0.035 | +0.046 | +0.026 | +2.402 | False |
| GM23815 | PartialReprog | 10 | 50 | -1.208 | +0.128 | -0.043 | -0.007 | +0.062 | +8.184 | False |
| GM23815 | Pluripotency | 10 | 181 | -10.192 | -0.009 | +0.008 | +0.018 | +0.004 | +11.937 | False |

Pooled endpoint nulls:

| cell_line | label | instrument | endpoint | n_cells_d0 | n_cells_end | score_d0 | score_end | decline | p | n_random_ge_real | pass_p | ok | reason | null_kind |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | NonReprog | frozen_ruler | d0→d7 | +65.000 | +1932.000 | -0.077 | +7.015 | -7.093 | +0.925 | +185.000 | False | True | NA | permuted_ruler_weights |
| GM00731 | NonReprog | frozen_ruler | d0→d10 | +65.000 | +3114.000 | -0.077 | +9.474 | -9.551 | +0.980 | +196.000 | False | True | NA | permuted_ruler_weights |
| GM00731 | Fibroblast | frozen_ruler | d0→d7 | +4939.000 | +14.000 | +5.253 | +8.596 | -3.343 | +0.866 | +173.000 | False | True | NA | permuted_ruler_weights |
| GM00731 | Fibroblast | frozen_ruler | d0→d10 | +4939.000 | +1030.000 | +5.253 | +7.873 | -2.620 | +0.995 | +199.000 | False | True | NA | permuted_ruler_weights |
| GM00731 | Pluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | Pluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | PartialReprog | frozen_ruler | d0→d7 | +3.000 | +578.000 | +1.397 | +4.737 | -3.339 | +0.403 | +80.000 | False | True | NA | permuted_ruler_weights |
| GM00731 | PartialReprog | frozen_ruler | d0→d10 | +3.000 | +31.000 | +1.397 | +4.933 | -3.535 | +0.483 | +96.000 | False | True | NA | permuted_ruler_weights |
| GM00731 | EarlyPluripotency | frozen_ruler | d0→d7 | +2.000 | +3917.000 | -4.307 | -0.739 | -3.567 | +0.194 | +38.000 | False | True | NA | permuted_ruler_weights |
| GM00731 | EarlyPluripotency | frozen_ruler | d0→d10 | +2.000 | +85.000 | -4.307 | -0.190 | -4.117 | +0.204 | +40.000 | False | True | NA | permuted_ruler_weights |
| GM00731 | NonReprog | pluri_primary | d0→d7 | +65.000 | +1932.000 | +0.141 | +3.203 | -3.062 | +0.940 | +188.000 | False | True | NA | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | NonReprog | pluri_primary | d0→d10 | +65.000 | +3114.000 | +0.141 | +1.688 | -1.548 | +0.791 | +158.000 | False | True | NA | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | Fibroblast | pluri_primary | d0→d7 | +4939.000 | +14.000 | -1.183 | -0.495 | -0.688 | +0.542 | +108.000 | False | True | NA | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | Fibroblast | pluri_primary | d0→d10 | +4939.000 | +1030.000 | -1.183 | +2.789 | -3.972 | +1.000 | +200.000 | False | True | NA | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | Pluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | Pluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | PartialReprog | pluri_primary | d0→d7 | +3.000 | +578.000 | -1.044 | +6.796 | -7.840 | +0.896 | +179.000 | False | True | NA | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | PartialReprog | pluri_primary | d0→d10 | +3.000 | +31.000 | -1.044 | +4.272 | -5.316 | +0.811 | +162.000 | False | True | NA | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | EarlyPluripotency | pluri_primary | d0→d7 | +2.000 | +3917.000 | +9.508 | +9.059 | +0.449 | +0.463 | +92.000 | False | True | NA | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | EarlyPluripotency | pluri_primary | d0→d10 | +2.000 | +85.000 | +9.508 | +10.334 | -0.825 | +0.572 | +114.000 | False | True | NA | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | NonReprog | md_score | d0→d7 | +65.000 | +1932.000 | +0.350 | +0.445 | -0.095 | +0.184 | +36.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | NonReprog | md_score | d0→d10 | +65.000 | +3114.000 | +0.350 | +0.434 | -0.084 | +0.070 | +13.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | Fibroblast | md_score | d0→d7 | +4939.000 | +14.000 | +0.506 | +0.425 | +0.081 | +0.159 | +31.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | Fibroblast | md_score | d0→d10 | +4939.000 | +1030.000 | +0.506 | +0.434 | +0.072 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | Pluripotency | md_score | d0→d7 | +0.000 | +354.000 | NA | NA | NA | NA | NA | NA | False | n_cells<1 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | Pluripotency | md_score | d0→d10 | +0.000 | +80.000 | NA | NA | NA | NA | NA | NA | False | n_cells<1 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | PartialReprog | md_score | d0→d7 | +3.000 | +578.000 | +0.291 | +0.284 | +0.007 | +0.015 | +2.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | PartialReprog | md_score | d0→d10 | +3.000 | +31.000 | +0.291 | +0.297 | -0.006 | +0.886 | +177.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | EarlyPluripotency | md_score | d0→d7 | +2.000 | +3917.000 | +0.249 | +0.152 | +0.096 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | EarlyPluripotency | md_score | d0→d10 | +2.000 | +85.000 | +0.249 | +0.079 | +0.170 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | NonReprog | tgfb_score | d0→d7 | +65.000 | +1932.000 | +0.086 | +0.048 | +0.038 | +0.164 | +32.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | NonReprog | tgfb_score | d0→d10 | +65.000 | +3114.000 | +0.086 | +0.001 | +0.085 | +0.030 | +5.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | Fibroblast | tgfb_score | d0→d7 | +4939.000 | +14.000 | +0.042 | +0.064 | -0.022 | +0.811 | +162.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | Fibroblast | tgfb_score | d0→d10 | +4939.000 | +1030.000 | +0.042 | -0.009 | +0.051 | +0.030 | +5.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | Pluripotency | tgfb_score | d0→d7 | +0.000 | +354.000 | NA | NA | NA | NA | NA | NA | False | n_cells<1 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | Pluripotency | tgfb_score | d0→d10 | +0.000 | +80.000 | NA | NA | NA | NA | NA | NA | False | n_cells<1 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | PartialReprog | tgfb_score | d0→d7 | +3.000 | +578.000 | -0.028 | -0.013 | -0.015 | +0.537 | +107.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | PartialReprog | tgfb_score | d0→d10 | +3.000 | +31.000 | -0.028 | -0.013 | -0.015 | +0.692 | +138.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | EarlyPluripotency | tgfb_score | d0→d7 | +2.000 | +3917.000 | -0.021 | -0.026 | +0.005 | +0.338 | +67.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | EarlyPluripotency | tgfb_score | d0→d10 | +2.000 | +85.000 | -0.021 | -0.086 | +0.065 | +0.139 | +27.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | NonReprog | age_up | d0→d7 | +65.000 | +1932.000 | +0.097 | +0.046 | +0.051 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | NonReprog | age_up | d0→d10 | +65.000 | +3114.000 | +0.097 | +0.055 | +0.042 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | Fibroblast | age_up | d0→d7 | +4939.000 | +14.000 | +0.055 | +0.037 | +0.018 | +0.010 | +1.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | Fibroblast | age_up | d0→d10 | +4939.000 | +1030.000 | +0.055 | +0.020 | +0.035 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | Pluripotency | age_up | d0→d7 | +0.000 | +354.000 | NA | NA | NA | NA | NA | NA | False | n_cells<1 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | Pluripotency | age_up | d0→d10 | +0.000 | +80.000 | NA | NA | NA | NA | NA | NA | False | n_cells<1 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | PartialReprog | age_up | d0→d7 | +3.000 | +578.000 | +0.009 | +0.011 | -0.002 | +0.343 | +68.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | PartialReprog | age_up | d0→d10 | +3.000 | +31.000 | +0.009 | +0.007 | +0.002 | +0.896 | +179.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | EarlyPluripotency | age_up | d0→d7 | +2.000 | +3917.000 | +0.064 | +0.033 | +0.031 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | EarlyPluripotency | age_up | d0→d10 | +2.000 | +85.000 | +0.064 | +0.016 | +0.048 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | NonReprog | age_down | d0→d7 | +65.000 | +1932.000 | -0.040 | +0.020 | -0.060 | +1.000 | +200.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | NonReprog | age_down | d0→d10 | +65.000 | +3114.000 | -0.040 | +0.023 | -0.063 | +1.000 | +200.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | Fibroblast | age_down | d0→d7 | +4939.000 | +14.000 | +0.026 | +0.025 | +0.001 | +0.990 | +198.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | Fibroblast | age_down | d0→d10 | +4939.000 | +1030.000 | +0.026 | +0.049 | -0.023 | +1.000 | +200.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | Pluripotency | age_down | d0→d7 | +0.000 | +354.000 | NA | NA | NA | NA | NA | NA | False | n_cells<1 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | Pluripotency | age_down | d0→d10 | +0.000 | +80.000 | NA | NA | NA | NA | NA | NA | False | n_cells<1 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | PartialReprog | age_down | d0→d7 | +3.000 | +578.000 | +0.032 | +0.077 | -0.046 | +1.000 | +200.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | PartialReprog | age_down | d0→d10 | +3.000 | +31.000 | +0.032 | +0.038 | -0.006 | +1.000 | +200.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | EarlyPluripotency | age_down | d0→d7 | +2.000 | +3917.000 | -0.001 | +0.034 | -0.035 | +0.925 | +185.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | EarlyPluripotency | age_down | d0→d10 | +2.000 | +85.000 | -0.001 | +0.051 | -0.052 | +1.000 | +200.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | NonReprog | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | NonReprog | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | Fibroblast | frozen_ruler | d0→d7 | +7703.000 | +2233.000 | +0.241 | -1.817 | +2.058 | +0.020 | +3.000 | True | True | NA | permuted_ruler_weights |
| GM23815 | Fibroblast | frozen_ruler | d0→d10 | +7703.000 | +1373.000 | +0.241 | +0.181 | +0.060 | +0.557 | +111.000 | False | True | NA | permuted_ruler_weights |
| GM23815 | Pluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | Pluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | PartialReprog | frozen_ruler | d0→d7 | +9.000 | +309.000 | +2.361 | -2.809 | +5.170 | +0.080 | +15.000 | False | True | NA | permuted_ruler_weights |
| GM23815 | PartialReprog | frozen_ruler | d0→d10 | +9.000 | +50.000 | +2.361 | -1.208 | +3.569 | +0.159 | +31.000 | False | True | NA | permuted_ruler_weights |
| GM23815 | EarlyPluripotency | frozen_ruler | d0→d7 | +56.000 | +8172.000 | -5.297 | -8.203 | +2.906 | +0.124 | +24.000 | False | True | NA | permuted_ruler_weights |
| GM23815 | EarlyPluripotency | frozen_ruler | d0→d10 | +56.000 | +4417.000 | -5.297 | -3.487 | -1.810 | +0.473 | +94.000 | False | True | NA | permuted_ruler_weights |
| GM23815 | NonReprog | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | NonReprog | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | Fibroblast | pluri_primary | d0→d7 | +7703.000 | +2233.000 | -0.795 | +4.290 | -5.085 | +1.000 | +200.000 | False | True | NA | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | Fibroblast | pluri_primary | d0→d10 | +7703.000 | +1373.000 | -0.795 | +2.943 | -3.738 | +1.000 | +200.000 | False | True | NA | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | Pluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | Pluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | PartialReprog | pluri_primary | d0→d7 | +9.000 | +309.000 | -2.091 | +7.891 | -9.982 | +0.990 | +198.000 | False | True | NA | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | PartialReprog | pluri_primary | d0→d10 | +9.000 | +50.000 | -2.091 | +8.184 | -10.275 | +0.980 | +196.000 | False | True | NA | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | EarlyPluripotency | pluri_primary | d0→d7 | +56.000 | +8172.000 | +0.733 | +11.120 | -10.387 | +0.980 | +196.000 | False | True | NA | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | EarlyPluripotency | pluri_primary | d0→d10 | +56.000 | +4417.000 | +0.733 | +8.721 | -7.988 | +0.960 | +192.000 | False | True | NA | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | NonReprog | md_score | d0→d7 | +0.000 | +8.000 | NA | NA | NA | NA | NA | NA | False | n_cells<1 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | NonReprog | md_score | d0→d10 | +0.000 | +1794.000 | NA | NA | NA | NA | NA | NA | False | n_cells<1 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | Fibroblast | md_score | d0→d7 | +7703.000 | +2233.000 | +0.479 | +0.410 | +0.069 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | Fibroblast | md_score | d0→d10 | +7703.000 | +1373.000 | +0.479 | +0.379 | +0.099 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | Pluripotency | md_score | d0→d7 | +0.000 | +977.000 | NA | NA | NA | NA | NA | NA | False | n_cells<1 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | Pluripotency | md_score | d0→d10 | +0.000 | +181.000 | NA | NA | NA | NA | NA | NA | False | n_cells<1 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | PartialReprog | md_score | d0→d7 | +9.000 | +309.000 | +0.316 | +0.174 | +0.143 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | PartialReprog | md_score | d0→d10 | +9.000 | +50.000 | +0.316 | +0.128 | +0.188 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | EarlyPluripotency | md_score | d0→d7 | +56.000 | +8172.000 | +0.379 | +0.077 | +0.303 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | EarlyPluripotency | md_score | d0→d10 | +56.000 | +4417.000 | +0.379 | +0.166 | +0.213 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | NonReprog | tgfb_score | d0→d7 | +0.000 | +8.000 | NA | NA | NA | NA | NA | NA | False | n_cells<1 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | NonReprog | tgfb_score | d0→d10 | +0.000 | +1794.000 | NA | NA | NA | NA | NA | NA | False | n_cells<1 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | Fibroblast | tgfb_score | d0→d7 | +7703.000 | +2233.000 | +0.029 | -0.015 | +0.044 | +0.020 | +3.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | Fibroblast | tgfb_score | d0→d10 | +7703.000 | +1373.000 | +0.029 | +0.037 | -0.008 | +0.652 | +130.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | Pluripotency | tgfb_score | d0→d7 | +0.000 | +977.000 | NA | NA | NA | NA | NA | NA | False | n_cells<1 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | Pluripotency | tgfb_score | d0→d10 | +0.000 | +181.000 | NA | NA | NA | NA | NA | NA | False | n_cells<1 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | PartialReprog | tgfb_score | d0→d7 | +9.000 | +309.000 | +0.006 | -0.058 | +0.064 | +0.030 | +5.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | PartialReprog | tgfb_score | d0→d10 | +9.000 | +50.000 | +0.006 | -0.043 | +0.049 | +0.124 | +24.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | EarlyPluripotency | tgfb_score | d0→d7 | +56.000 | +8172.000 | +0.077 | -0.058 | +0.136 | +0.020 | +3.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | EarlyPluripotency | tgfb_score | d0→d10 | +56.000 | +4417.000 | +0.077 | -0.023 | +0.100 | +0.060 | +11.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | NonReprog | age_up | d0→d7 | +0.000 | +8.000 | NA | NA | NA | NA | NA | NA | False | n_cells<1 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | NonReprog | age_up | d0→d10 | +0.000 | +1794.000 | NA | NA | NA | NA | NA | NA | False | n_cells<1 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | Fibroblast | age_up | d0→d7 | +7703.000 | +2233.000 | +0.039 | +0.033 | +0.006 | +0.274 | +54.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | Fibroblast | age_up | d0→d10 | +7703.000 | +1373.000 | +0.039 | +0.012 | +0.027 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | Pluripotency | age_up | d0→d7 | +0.000 | +977.000 | NA | NA | NA | NA | NA | NA | False | n_cells<1 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | Pluripotency | age_up | d0→d10 | +0.000 | +181.000 | NA | NA | NA | NA | NA | NA | False | n_cells<1 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | PartialReprog | age_up | d0→d7 | +9.000 | +309.000 | -0.010 | +0.009 | -0.019 | +1.000 | +200.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | PartialReprog | age_up | d0→d10 | +9.000 | +50.000 | -0.010 | -0.007 | -0.003 | +0.915 | +183.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | EarlyPluripotency | age_up | d0→d7 | +56.000 | +8172.000 | +0.083 | +0.017 | +0.066 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | EarlyPluripotency | age_up | d0→d10 | +56.000 | +4417.000 | +0.083 | +0.009 | +0.074 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | NonReprog | age_down | d0→d7 | +0.000 | +8.000 | NA | NA | NA | NA | NA | NA | False | n_cells<1 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | NonReprog | age_down | d0→d10 | +0.000 | +1794.000 | NA | NA | NA | NA | NA | NA | False | n_cells<1 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | Fibroblast | age_down | d0→d7 | +7703.000 | +2233.000 | +0.029 | +0.025 | +0.004 | +0.940 | +188.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | Fibroblast | age_down | d0→d10 | +7703.000 | +1373.000 | +0.029 | +0.047 | -0.018 | +1.000 | +200.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | Pluripotency | age_down | d0→d7 | +0.000 | +977.000 | NA | NA | NA | NA | NA | NA | False | n_cells<1 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | Pluripotency | age_down | d0→d10 | +0.000 | +181.000 | NA | NA | NA | NA | NA | NA | False | n_cells<1 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | PartialReprog | age_down | d0→d7 | +9.000 | +309.000 | +0.093 | +0.056 | +0.038 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | PartialReprog | age_down | d0→d10 | +9.000 | +50.000 | +0.093 | +0.062 | +0.031 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | EarlyPluripotency | age_down | d0→d7 | +56.000 | +8172.000 | -0.038 | +0.044 | -0.082 | +1.000 | +200.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | EarlyPluripotency | age_down | d0→d10 | +56.000 | +4417.000 | -0.038 | +0.050 | -0.087 | +1.000 | +200.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |

### (b) per cluster at ≥20 cells (FINDINGS_MD2.md threshold, unchanged)

| cell_line | cluster | label | day | n_cells | below_min_cells | age_score | md_score | tgfb_score | age_up | age_down | pluri_primary |
|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | 0 | Fibroblast | 0 | 3998 | False | +7.037 | +0.513 | +0.040 | +0.058 | +0.018 | -0.743 |
| GM00731 | 3 | PartialReprog | 0 | 1 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 4 | NonReprog | 0 | 9 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 5 | PartialReprog | 0 | 2 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 10 | NonReprog | 0 | 56 | False | -1.425 | +0.356 | +0.117 | +0.112 | -0.052 | +1.860 |
| GM00731 | 12 | Fibroblast | 0 | 569 | False | +7.865 | +0.510 | +0.049 | +0.062 | +0.032 | -1.198 |
| GM00731 | 15 | Fibroblast | 0 | 372 | False | +4.533 | +0.428 | +0.050 | +0.009 | +0.097 | +1.003 |
| GM00731 | 16 | EarlyPluripotency | 0 | 2 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 0 | Fibroblast | 3 | 4 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 3 | PartialReprog | 3 | 1589 | False | +6.055 | +0.308 | +0.038 | -0.024 | +0.097 | +5.941 |
| GM00731 | 4 | NonReprog | 3 | 259 | False | +9.226 | +0.531 | +0.060 | +0.030 | +0.031 | +3.900 |
| GM00731 | 5 | PartialReprog | 3 | 1523 | False | +4.013 | +0.060 | -0.096 | -0.015 | +0.071 | +9.042 |
| GM00731 | 7 | EarlyPluripotency | 3 | 108 | False | +5.193 | +0.163 | -0.047 | +0.012 | +0.051 | +7.227 |
| GM00731 | 8 | PartialReprog | 3 | 973 | False | -0.379 | +0.020 | -0.098 | -0.013 | +0.081 | +8.942 |
| GM00731 | 10 | NonReprog | 3 | 89 | False | +7.215 | +0.316 | +0.117 | +0.024 | +0.011 | +5.974 |
| GM00731 | 11 | PartialReprog | 3 | 557 | False | +3.280 | +0.080 | -0.089 | -0.016 | +0.068 | +7.831 |
| GM00731 | 13 | PartialReprog | 3 | 442 | False | +4.836 | +0.079 | -0.031 | -0.024 | +0.083 | +9.860 |
| GM00731 | 14 | Pluripotency | 3 | 5 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 15 | Fibroblast | 3 | 1 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 16 | EarlyPluripotency | 3 | 302 | False | +0.353 | +0.009 | -0.015 | +0.006 | +0.002 | +12.562 |
| GM00731 | 18 | PartialReprog | 3 | 136 | False | +3.692 | +0.171 | -0.016 | -0.020 | +0.065 | +7.768 |
| GM00731 | 0 | Fibroblast | 7 | 9 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 1 | NonReprog | 7 | 9 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 2 | EarlyPluripotency | 7 | 2092 | False | +0.098 | +0.167 | -0.004 | +0.033 | +0.029 | +8.439 |
| GM00731 | 3 | PartialReprog | 7 | 444 | False | +6.644 | +0.343 | +0.006 | +0.008 | +0.090 | +4.785 |
| GM00731 | 4 | NonReprog | 7 | 1732 | False | +8.658 | +0.456 | +0.045 | +0.046 | +0.024 | +3.152 |
| GM00731 | 5 | PartialReprog | 7 | 48 | False | +0.508 | +0.052 | -0.088 | +0.030 | +0.020 | +11.007 |
| GM00731 | 6 | Fibroblast | 7 | 1 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 7 | EarlyPluripotency | 7 | 884 | False | +4.745 | +0.299 | +0.001 | +0.043 | +0.042 | +7.337 |
| GM00731 | 8 | PartialReprog | 7 | 13 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 9 | EarlyPluripotency | 7 | 909 | False | -2.000 | -0.024 | -0.106 | +0.023 | +0.041 | +11.825 |
| GM00731 | 10 | NonReprog | 7 | 191 | False | +3.966 | +0.347 | +0.079 | +0.053 | -0.017 | +5.835 |
| GM00731 | 11 | PartialReprog | 7 | 64 | False | +4.110 | +0.117 | -0.064 | +0.017 | +0.043 | +6.920 |
| GM00731 | 12 | Fibroblast | 7 | 1 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 13 | PartialReprog | 7 | 3 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 14 | Pluripotency | 7 | 354 | False | -5.699 | +0.011 | +0.005 | +0.056 | -0.022 | +11.693 |
| GM00731 | 15 | Fibroblast | 7 | 3 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 16 | EarlyPluripotency | 7 | 32 | False | +1.891 | +0.188 | +0.033 | +0.038 | +0.004 | +4.024 |
| GM00731 | 18 | PartialReprog | 7 | 6 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 0 | Fibroblast | 10 | 1 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 1 | NonReprog | 10 | 2313 | False | +11.686 | +0.446 | -0.018 | +0.058 | +0.027 | +1.788 |
| GM00731 | 2 | EarlyPluripotency | 10 | 4 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 3 | PartialReprog | 10 | 5 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 4 | NonReprog | 10 | 31 | False | +13.645 | +0.408 | -0.045 | +0.028 | +0.026 | +2.259 |
| GM00731 | 5 | PartialReprog | 10 | 21 | False | +2.744 | +0.320 | +0.008 | +0.007 | +0.023 | +6.696 |
| GM00731 | 6 | Fibroblast | 10 | 1029 | False | +9.102 | +0.434 | -0.009 | +0.020 | +0.049 | +2.984 |
| GM00731 | 7 | EarlyPluripotency | 10 | 22 | False | +6.501 | +0.289 | -0.071 | +0.052 | +0.031 | +6.230 |
| GM00731 | 9 | EarlyPluripotency | 10 | 52 | False | -0.290 | -0.040 | -0.109 | +0.001 | +0.067 | +12.005 |
| GM00731 | 10 | NonReprog | 10 | 527 | False | +7.934 | +0.417 | +0.102 | +0.063 | -0.023 | +3.313 |
| GM00731 | 11 | PartialReprog | 10 | 5 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 14 | Pluripotency | 10 | 80 | False | -2.372 | +0.069 | +0.019 | +0.051 | -0.019 | +10.355 |
| GM00731 | 16 | EarlyPluripotency | 10 | 7 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 17 | NonReprog | 10 | 243 | False | +8.127 | +0.356 | -0.036 | +0.016 | +0.092 | +1.015 |
| GM23815 | 0 | Fibroblast | 0 | 6454 | False | +3.253 | +0.483 | +0.033 | +0.042 | +0.021 | -0.154 |
| GM23815 | 4 | Fibroblast | 0 | 2 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 5 | PartialReprog | 0 | 9 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 13 | EarlyPluripotency | 0 | 56 | False | -2.160 | +0.379 | +0.077 | +0.083 | -0.038 | +1.778 |
| GM23815 | 14 | Fibroblast | 0 | 886 | False | +1.849 | +0.420 | +0.016 | +0.005 | +0.096 | +0.918 |
| GM23815 | 16 | Fibroblast | 0 | 361 | False | +3.401 | +0.547 | -0.005 | +0.057 | +0.010 | -0.696 |
| GM23815 | 0 | Fibroblast | 3 | 19 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 1 | EarlyPluripotency | 3 | 3 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 3 | EarlyPluripotency | 3 | 3 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 4 | Fibroblast | 3 | 21 | False | +8.218 | +0.391 | +0.045 | +0.006 | +0.045 | +3.558 |
| GM23815 | 5 | PartialReprog | 3 | 2039 | False | +1.227 | +0.219 | +0.029 | -0.030 | +0.112 | +6.379 |
| GM23815 | 7 | PartialReprog | 3 | 2039 | False | -3.618 | -0.011 | -0.084 | -0.017 | +0.091 | +9.994 |
| GM23815 | 8 | NonReprog | 3 | 10 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 10 | Pluripotency | 3 | 91 | False | -8.931 | +0.028 | -0.012 | +0.021 | +0.007 | +10.539 |
| GM23815 | 11 | PartialReprog | 3 | 1065 | False | -0.322 | +0.007 | -0.065 | -0.018 | +0.078 | +10.360 |
| GM23815 | 12 | PartialReprog | 3 | 970 | False | +0.720 | +0.095 | -0.043 | -0.016 | +0.070 | +7.424 |
| GM23815 | 13 | EarlyPluripotency | 3 | 246 | False | +2.602 | +0.034 | +0.023 | +0.006 | +0.008 | +9.709 |
| GM23815 | 14 | Fibroblast | 3 | 5 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 16 | Fibroblast | 3 | 1 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 17 | EarlyPluripotency | 3 | 64 | False | -1.597 | +0.156 | -0.011 | +0.001 | +0.072 | +7.501 |
| GM23815 | 0 | Fibroblast | 7 | 14 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 1 | EarlyPluripotency | 7 | 3240 | False | -3.860 | +0.207 | -0.041 | +0.028 | +0.035 | +7.946 |
| GM23815 | 2 | EarlyPluripotency | 7 | 1 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 3 | EarlyPluripotency | 7 | 2537 | False | -6.390 | +0.003 | -0.069 | +0.017 | +0.042 | +12.127 |
| GM23815 | 4 | Fibroblast | 7 | 2210 | False | +1.091 | +0.410 | -0.015 | +0.033 | +0.024 | +4.365 |
| GM23815 | 5 | PartialReprog | 7 | 118 | False | +1.778 | +0.314 | -0.018 | -0.003 | +0.083 | +4.847 |
| GM23815 | 6 | EarlyPluripotency | 7 | 1819 | False | -5.449 | -0.048 | -0.083 | -0.006 | +0.072 | +12.685 |
| GM23815 | 7 | PartialReprog | 7 | 20 | False | -6.363 | +0.006 | -0.107 | +0.024 | +0.046 | +11.332 |
| GM23815 | 8 | NonReprog | 7 | 8 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 9 | Fibroblast | 7 | 3 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 10 | Pluripotency | 7 | 977 | False | -8.857 | -0.059 | +0.006 | +0.027 | -0.004 | +12.656 |
| GM23815 | 11 | PartialReprog | 7 | 65 | False | -2.332 | +0.023 | -0.101 | +0.021 | +0.037 | +10.838 |
| GM23815 | 12 | PartialReprog | 7 | 106 | False | +0.145 | +0.143 | -0.067 | +0.011 | +0.039 | +7.850 |
| GM23815 | 13 | EarlyPluripotency | 7 | 164 | False | -7.076 | +0.159 | +0.056 | +0.035 | +0.001 | +8.758 |
| GM23815 | 14 | Fibroblast | 7 | 5 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 16 | Fibroblast | 7 | 1 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 17 | EarlyPluripotency | 7 | 122 | False | -2.776 | +0.083 | -0.072 | +0.041 | +0.035 | +9.561 |
| GM23815 | 18 | EarlyPluripotency | 7 | 289 | False | -7.347 | +0.001 | -0.055 | +0.029 | +0.025 | +11.933 |
| GM23815 | 0 | Fibroblast | 10 | 2 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 1 | EarlyPluripotency | 10 | 8 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 2 | EarlyPluripotency | 10 | 2730 | False | -1.570 | +0.155 | -0.037 | +0.005 | +0.052 | +7.637 |
| GM23815 | 3 | EarlyPluripotency | 10 | 30 | False | -3.984 | -0.036 | -0.096 | +0.002 | +0.041 | +11.657 |
| GM23815 | 4 | Fibroblast | 10 | 9 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 5 | PartialReprog | 10 | 5 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 6 | EarlyPluripotency | 10 | 330 | False | -4.946 | -0.115 | -0.078 | -0.024 | +0.081 | +14.675 |
| GM23815 | 7 | PartialReprog | 10 | 10 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 8 | NonReprog | 10 | 1794 | False | +5.408 | +0.423 | +0.035 | +0.046 | +0.026 | +2.530 |
| GM23815 | 9 | Fibroblast | 10 | 1362 | False | +3.204 | +0.380 | +0.036 | +0.012 | +0.047 | +3.122 |
| GM23815 | 10 | Pluripotency | 10 | 181 | False | -6.660 | -0.009 | +0.008 | +0.018 | +0.004 | +11.955 |
| GM23815 | 11 | PartialReprog | 10 | 10 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 12 | PartialReprog | 10 | 25 | False | +3.381 | +0.194 | +0.007 | -0.014 | +0.073 | +4.677 |
| GM23815 | 13 | EarlyPluripotency | 10 | 432 | False | -3.794 | +0.235 | +0.071 | +0.046 | -0.017 | +6.874 |
| GM23815 | 15 | EarlyPluripotency | 10 | 761 | False | +1.379 | +0.278 | -0.002 | +0.016 | +0.066 | +7.499 |
| GM23815 | 17 | EarlyPluripotency | 10 | 121 | False | +1.877 | +0.299 | +0.001 | +0.028 | +0.044 | +6.431 |
| GM23815 | 18 | EarlyPluripotency | 10 | 5 | True | NA | NA | NA | NA | NA | NA |

| cell_line | cluster | label | instrument | endpoint | n_cells_d0 | n_cells_end | score_d0 | score_end | decline | p | n_random_ge_real | pass_p | ok | reason | null_kind |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | +0.000 | Fibroblast | frozen_ruler | d0→d7 | +3998.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| GM00731 | +0.000 | Fibroblast | frozen_ruler | d0→d10 | +3998.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| GM00731 | +1.000 | NonReprog | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +1.000 | NonReprog | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +2.000 | EarlyPluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +2.000 | EarlyPluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +3.000 | PartialReprog | frozen_ruler | d0→d7 | +1.000 | +444.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| GM00731 | +3.000 | PartialReprog | frozen_ruler | d0→d10 | +1.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| GM00731 | +4.000 | NonReprog | frozen_ruler | d0→d7 | +9.000 | +1732.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| GM00731 | +4.000 | NonReprog | frozen_ruler | d0→d10 | +9.000 | +31.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| GM00731 | +5.000 | PartialReprog | frozen_ruler | d0→d7 | +2.000 | +48.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| GM00731 | +5.000 | PartialReprog | frozen_ruler | d0→d10 | +2.000 | +21.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| GM00731 | +6.000 | Fibroblast | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +6.000 | Fibroblast | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +7.000 | EarlyPluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +7.000 | EarlyPluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +8.000 | PartialReprog | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +8.000 | PartialReprog | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +9.000 | EarlyPluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +9.000 | EarlyPluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +10.000 | NonReprog | frozen_ruler | d0→d7 | +56.000 | +191.000 | -1.425 | +3.966 | -5.391 | +0.930 | +186.000 | False | True | NA | permuted_ruler_weights |
| GM00731 | +10.000 | NonReprog | frozen_ruler | d0→d10 | +56.000 | +527.000 | -1.425 | +7.934 | -9.359 | +1.000 | +200.000 | False | True | NA | permuted_ruler_weights |
| GM00731 | +11.000 | PartialReprog | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +11.000 | PartialReprog | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +12.000 | Fibroblast | frozen_ruler | d0→d7 | +569.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| GM00731 | +12.000 | Fibroblast | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +13.000 | PartialReprog | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +13.000 | PartialReprog | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +14.000 | Pluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +14.000 | Pluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +15.000 | Fibroblast | frozen_ruler | d0→d7 | +372.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| GM00731 | +15.000 | Fibroblast | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +16.000 | EarlyPluripotency | frozen_ruler | d0→d7 | +2.000 | +32.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| GM00731 | +16.000 | EarlyPluripotency | frozen_ruler | d0→d10 | +2.000 | +7.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| GM00731 | +17.000 | NonReprog | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +17.000 | NonReprog | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +18.000 | PartialReprog | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +18.000 | PartialReprog | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +0.000 | Fibroblast | pluri_primary | d0→d7 | +3998.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +0.000 | Fibroblast | pluri_primary | d0→d10 | +3998.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +1.000 | NonReprog | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +1.000 | NonReprog | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +2.000 | EarlyPluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +2.000 | EarlyPluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +3.000 | PartialReprog | pluri_primary | d0→d7 | +1.000 | +444.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +3.000 | PartialReprog | pluri_primary | d0→d10 | +1.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +4.000 | NonReprog | pluri_primary | d0→d7 | +9.000 | +1732.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +4.000 | NonReprog | pluri_primary | d0→d10 | +9.000 | +31.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +5.000 | PartialReprog | pluri_primary | d0→d7 | +2.000 | +48.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +5.000 | PartialReprog | pluri_primary | d0→d10 | +2.000 | +21.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +6.000 | Fibroblast | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +6.000 | Fibroblast | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +7.000 | EarlyPluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +7.000 | EarlyPluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +8.000 | PartialReprog | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +8.000 | PartialReprog | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +9.000 | EarlyPluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +9.000 | EarlyPluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +10.000 | NonReprog | pluri_primary | d0→d7 | +56.000 | +191.000 | +1.860 | +5.835 | -3.975 | +0.970 | +194.000 | False | True | NA | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +10.000 | NonReprog | pluri_primary | d0→d10 | +56.000 | +527.000 | +1.860 | +3.313 | -1.453 | +0.831 | +166.000 | False | True | NA | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +11.000 | PartialReprog | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +11.000 | PartialReprog | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +12.000 | Fibroblast | pluri_primary | d0→d7 | +569.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +12.000 | Fibroblast | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +13.000 | PartialReprog | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +13.000 | PartialReprog | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +14.000 | Pluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +14.000 | Pluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +15.000 | Fibroblast | pluri_primary | d0→d7 | +372.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +15.000 | Fibroblast | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +16.000 | EarlyPluripotency | pluri_primary | d0→d7 | +2.000 | +32.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +16.000 | EarlyPluripotency | pluri_primary | d0→d10 | +2.000 | +7.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +17.000 | NonReprog | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +17.000 | NonReprog | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +18.000 | PartialReprog | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +18.000 | PartialReprog | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +0.000 | Fibroblast | md_score | d0→d7 | +3998.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +0.000 | Fibroblast | md_score | d0→d10 | +3998.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +1.000 | NonReprog | md_score | d0→d7 | +0.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +1.000 | NonReprog | md_score | d0→d10 | +0.000 | +2313.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +2.000 | EarlyPluripotency | md_score | d0→d7 | +0.000 | +2092.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +2.000 | EarlyPluripotency | md_score | d0→d10 | +0.000 | +4.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +3.000 | PartialReprog | md_score | d0→d7 | +1.000 | +444.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +3.000 | PartialReprog | md_score | d0→d10 | +1.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +4.000 | NonReprog | md_score | d0→d7 | +9.000 | +1732.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +4.000 | NonReprog | md_score | d0→d10 | +9.000 | +31.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +5.000 | PartialReprog | md_score | d0→d7 | +2.000 | +48.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +5.000 | PartialReprog | md_score | d0→d10 | +2.000 | +21.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +6.000 | Fibroblast | md_score | d0→d7 | +0.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +6.000 | Fibroblast | md_score | d0→d10 | +0.000 | +1029.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +7.000 | EarlyPluripotency | md_score | d0→d7 | +0.000 | +884.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +7.000 | EarlyPluripotency | md_score | d0→d10 | +0.000 | +22.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +8.000 | PartialReprog | md_score | d0→d7 | +0.000 | +13.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +8.000 | PartialReprog | md_score | d0→d10 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +9.000 | EarlyPluripotency | md_score | d0→d7 | +0.000 | +909.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +9.000 | EarlyPluripotency | md_score | d0→d10 | +0.000 | +52.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +10.000 | NonReprog | md_score | d0→d7 | +56.000 | +191.000 | +0.356 | +0.347 | +0.009 | +0.035 | +6.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +10.000 | NonReprog | md_score | d0→d10 | +56.000 | +527.000 | +0.356 | +0.417 | -0.061 | +0.886 | +177.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +11.000 | PartialReprog | md_score | d0→d7 | +0.000 | +64.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +11.000 | PartialReprog | md_score | d0→d10 | +0.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +12.000 | Fibroblast | md_score | d0→d7 | +569.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +12.000 | Fibroblast | md_score | d0→d10 | +569.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +13.000 | PartialReprog | md_score | d0→d7 | +0.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +13.000 | PartialReprog | md_score | d0→d10 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +14.000 | Pluripotency | md_score | d0→d7 | +0.000 | +354.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +14.000 | Pluripotency | md_score | d0→d10 | +0.000 | +80.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +15.000 | Fibroblast | md_score | d0→d7 | +372.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +15.000 | Fibroblast | md_score | d0→d10 | +372.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +16.000 | EarlyPluripotency | md_score | d0→d7 | +2.000 | +32.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +16.000 | EarlyPluripotency | md_score | d0→d10 | +2.000 | +7.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +17.000 | NonReprog | md_score | d0→d7 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +17.000 | NonReprog | md_score | d0→d10 | +0.000 | +243.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +18.000 | PartialReprog | md_score | d0→d7 | +0.000 | +6.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +18.000 | PartialReprog | md_score | d0→d10 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +0.000 | Fibroblast | tgfb_score | d0→d7 | +3998.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +0.000 | Fibroblast | tgfb_score | d0→d10 | +3998.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +1.000 | NonReprog | tgfb_score | d0→d7 | +0.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +1.000 | NonReprog | tgfb_score | d0→d10 | +0.000 | +2313.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +2.000 | EarlyPluripotency | tgfb_score | d0→d7 | +0.000 | +2092.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +2.000 | EarlyPluripotency | tgfb_score | d0→d10 | +0.000 | +4.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +3.000 | PartialReprog | tgfb_score | d0→d7 | +1.000 | +444.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +3.000 | PartialReprog | tgfb_score | d0→d10 | +1.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +4.000 | NonReprog | tgfb_score | d0→d7 | +9.000 | +1732.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +4.000 | NonReprog | tgfb_score | d0→d10 | +9.000 | +31.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +5.000 | PartialReprog | tgfb_score | d0→d7 | +2.000 | +48.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +5.000 | PartialReprog | tgfb_score | d0→d10 | +2.000 | +21.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +6.000 | Fibroblast | tgfb_score | d0→d7 | +0.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +6.000 | Fibroblast | tgfb_score | d0→d10 | +0.000 | +1029.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +7.000 | EarlyPluripotency | tgfb_score | d0→d7 | +0.000 | +884.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +7.000 | EarlyPluripotency | tgfb_score | d0→d10 | +0.000 | +22.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +8.000 | PartialReprog | tgfb_score | d0→d7 | +0.000 | +13.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +8.000 | PartialReprog | tgfb_score | d0→d10 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +9.000 | EarlyPluripotency | tgfb_score | d0→d7 | +0.000 | +909.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +9.000 | EarlyPluripotency | tgfb_score | d0→d10 | +0.000 | +52.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +10.000 | NonReprog | tgfb_score | d0→d7 | +56.000 | +191.000 | +0.117 | +0.079 | +0.038 | +0.119 | +23.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +10.000 | NonReprog | tgfb_score | d0→d10 | +56.000 | +527.000 | +0.117 | +0.102 | +0.015 | +0.294 | +58.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +11.000 | PartialReprog | tgfb_score | d0→d7 | +0.000 | +64.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +11.000 | PartialReprog | tgfb_score | d0→d10 | +0.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +12.000 | Fibroblast | tgfb_score | d0→d7 | +569.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +12.000 | Fibroblast | tgfb_score | d0→d10 | +569.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +13.000 | PartialReprog | tgfb_score | d0→d7 | +0.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +13.000 | PartialReprog | tgfb_score | d0→d10 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +14.000 | Pluripotency | tgfb_score | d0→d7 | +0.000 | +354.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +14.000 | Pluripotency | tgfb_score | d0→d10 | +0.000 | +80.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +15.000 | Fibroblast | tgfb_score | d0→d7 | +372.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +15.000 | Fibroblast | tgfb_score | d0→d10 | +372.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +16.000 | EarlyPluripotency | tgfb_score | d0→d7 | +2.000 | +32.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +16.000 | EarlyPluripotency | tgfb_score | d0→d10 | +2.000 | +7.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +17.000 | NonReprog | tgfb_score | d0→d7 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +17.000 | NonReprog | tgfb_score | d0→d10 | +0.000 | +243.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +18.000 | PartialReprog | tgfb_score | d0→d7 | +0.000 | +6.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +18.000 | PartialReprog | tgfb_score | d0→d10 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +0.000 | Fibroblast | age_up | d0→d7 | +3998.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +0.000 | Fibroblast | age_up | d0→d10 | +3998.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +1.000 | NonReprog | age_up | d0→d7 | +0.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +1.000 | NonReprog | age_up | d0→d10 | +0.000 | +2313.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +2.000 | EarlyPluripotency | age_up | d0→d7 | +0.000 | +2092.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +2.000 | EarlyPluripotency | age_up | d0→d10 | +0.000 | +4.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +3.000 | PartialReprog | age_up | d0→d7 | +1.000 | +444.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +3.000 | PartialReprog | age_up | d0→d10 | +1.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +4.000 | NonReprog | age_up | d0→d7 | +9.000 | +1732.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +4.000 | NonReprog | age_up | d0→d10 | +9.000 | +31.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +5.000 | PartialReprog | age_up | d0→d7 | +2.000 | +48.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +5.000 | PartialReprog | age_up | d0→d10 | +2.000 | +21.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +6.000 | Fibroblast | age_up | d0→d7 | +0.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +6.000 | Fibroblast | age_up | d0→d10 | +0.000 | +1029.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +7.000 | EarlyPluripotency | age_up | d0→d7 | +0.000 | +884.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +7.000 | EarlyPluripotency | age_up | d0→d10 | +0.000 | +22.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +8.000 | PartialReprog | age_up | d0→d7 | +0.000 | +13.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +8.000 | PartialReprog | age_up | d0→d10 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +9.000 | EarlyPluripotency | age_up | d0→d7 | +0.000 | +909.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +9.000 | EarlyPluripotency | age_up | d0→d10 | +0.000 | +52.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +10.000 | NonReprog | age_up | d0→d7 | +56.000 | +191.000 | +0.112 | +0.053 | +0.059 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +10.000 | NonReprog | age_up | d0→d10 | +56.000 | +527.000 | +0.112 | +0.063 | +0.049 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +11.000 | PartialReprog | age_up | d0→d7 | +0.000 | +64.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +11.000 | PartialReprog | age_up | d0→d10 | +0.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +12.000 | Fibroblast | age_up | d0→d7 | +569.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +12.000 | Fibroblast | age_up | d0→d10 | +569.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +13.000 | PartialReprog | age_up | d0→d7 | +0.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +13.000 | PartialReprog | age_up | d0→d10 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +14.000 | Pluripotency | age_up | d0→d7 | +0.000 | +354.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +14.000 | Pluripotency | age_up | d0→d10 | +0.000 | +80.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +15.000 | Fibroblast | age_up | d0→d7 | +372.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +15.000 | Fibroblast | age_up | d0→d10 | +372.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +16.000 | EarlyPluripotency | age_up | d0→d7 | +2.000 | +32.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +16.000 | EarlyPluripotency | age_up | d0→d10 | +2.000 | +7.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +17.000 | NonReprog | age_up | d0→d7 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +17.000 | NonReprog | age_up | d0→d10 | +0.000 | +243.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +18.000 | PartialReprog | age_up | d0→d7 | +0.000 | +6.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +18.000 | PartialReprog | age_up | d0→d10 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +0.000 | Fibroblast | age_down | d0→d7 | +3998.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +0.000 | Fibroblast | age_down | d0→d10 | +3998.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +1.000 | NonReprog | age_down | d0→d7 | +0.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +1.000 | NonReprog | age_down | d0→d10 | +0.000 | +2313.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +2.000 | EarlyPluripotency | age_down | d0→d7 | +0.000 | +2092.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +2.000 | EarlyPluripotency | age_down | d0→d10 | +0.000 | +4.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +3.000 | PartialReprog | age_down | d0→d7 | +1.000 | +444.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +3.000 | PartialReprog | age_down | d0→d10 | +1.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +4.000 | NonReprog | age_down | d0→d7 | +9.000 | +1732.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +4.000 | NonReprog | age_down | d0→d10 | +9.000 | +31.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +5.000 | PartialReprog | age_down | d0→d7 | +2.000 | +48.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +5.000 | PartialReprog | age_down | d0→d10 | +2.000 | +21.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +6.000 | Fibroblast | age_down | d0→d7 | +0.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +6.000 | Fibroblast | age_down | d0→d10 | +0.000 | +1029.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +7.000 | EarlyPluripotency | age_down | d0→d7 | +0.000 | +884.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +7.000 | EarlyPluripotency | age_down | d0→d10 | +0.000 | +22.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +8.000 | PartialReprog | age_down | d0→d7 | +0.000 | +13.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +8.000 | PartialReprog | age_down | d0→d10 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +9.000 | EarlyPluripotency | age_down | d0→d7 | +0.000 | +909.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +9.000 | EarlyPluripotency | age_down | d0→d10 | +0.000 | +52.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +10.000 | NonReprog | age_down | d0→d7 | +56.000 | +191.000 | -0.052 | -0.017 | -0.035 | +1.000 | +200.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +10.000 | NonReprog | age_down | d0→d10 | +56.000 | +527.000 | -0.052 | -0.023 | -0.028 | +1.000 | +200.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +11.000 | PartialReprog | age_down | d0→d7 | +0.000 | +64.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +11.000 | PartialReprog | age_down | d0→d10 | +0.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +12.000 | Fibroblast | age_down | d0→d7 | +569.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +12.000 | Fibroblast | age_down | d0→d10 | +569.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +13.000 | PartialReprog | age_down | d0→d7 | +0.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +13.000 | PartialReprog | age_down | d0→d10 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +14.000 | Pluripotency | age_down | d0→d7 | +0.000 | +354.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +14.000 | Pluripotency | age_down | d0→d10 | +0.000 | +80.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +15.000 | Fibroblast | age_down | d0→d7 | +372.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +15.000 | Fibroblast | age_down | d0→d10 | +372.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +16.000 | EarlyPluripotency | age_down | d0→d7 | +2.000 | +32.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +16.000 | EarlyPluripotency | age_down | d0→d10 | +2.000 | +7.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +17.000 | NonReprog | age_down | d0→d7 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +17.000 | NonReprog | age_down | d0→d10 | +0.000 | +243.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +18.000 | PartialReprog | age_down | d0→d7 | +0.000 | +6.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +18.000 | PartialReprog | age_down | d0→d10 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +0.000 | Fibroblast | frozen_ruler | d0→d7 | +6454.000 | +14.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| GM23815 | +0.000 | Fibroblast | frozen_ruler | d0→d10 | +6454.000 | +2.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| GM23815 | +1.000 | EarlyPluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +1.000 | EarlyPluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +2.000 | EarlyPluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +2.000 | EarlyPluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +3.000 | EarlyPluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +3.000 | EarlyPluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +4.000 | Fibroblast | frozen_ruler | d0→d7 | +2.000 | +2210.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| GM23815 | +4.000 | Fibroblast | frozen_ruler | d0→d10 | +2.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| GM23815 | +5.000 | PartialReprog | frozen_ruler | d0→d7 | +9.000 | +118.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| GM23815 | +5.000 | PartialReprog | frozen_ruler | d0→d10 | +9.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| GM23815 | +6.000 | EarlyPluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +6.000 | EarlyPluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +7.000 | PartialReprog | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +7.000 | PartialReprog | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +8.000 | NonReprog | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +8.000 | NonReprog | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +9.000 | Fibroblast | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +9.000 | Fibroblast | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +10.000 | Pluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +10.000 | Pluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +11.000 | PartialReprog | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +11.000 | PartialReprog | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +12.000 | PartialReprog | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +12.000 | PartialReprog | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +13.000 | EarlyPluripotency | frozen_ruler | d0→d7 | +56.000 | +164.000 | -2.160 | -7.076 | +4.917 | +0.050 | +9.000 | True | True | NA | permuted_ruler_weights |
| GM23815 | +13.000 | EarlyPluripotency | frozen_ruler | d0→d10 | +56.000 | +432.000 | -2.160 | -3.794 | +1.634 | +0.159 | +31.000 | False | True | NA | permuted_ruler_weights |
| GM23815 | +14.000 | Fibroblast | frozen_ruler | d0→d7 | +886.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| GM23815 | +14.000 | Fibroblast | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +15.000 | EarlyPluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +15.000 | EarlyPluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +16.000 | Fibroblast | frozen_ruler | d0→d7 | +361.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| GM23815 | +16.000 | Fibroblast | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +17.000 | EarlyPluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +17.000 | EarlyPluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +18.000 | EarlyPluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +18.000 | EarlyPluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +0.000 | Fibroblast | pluri_primary | d0→d7 | +6454.000 | +14.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +0.000 | Fibroblast | pluri_primary | d0→d10 | +6454.000 | +2.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +1.000 | EarlyPluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +1.000 | EarlyPluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +2.000 | EarlyPluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +2.000 | EarlyPluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +3.000 | EarlyPluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +3.000 | EarlyPluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +4.000 | Fibroblast | pluri_primary | d0→d7 | +2.000 | +2210.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +4.000 | Fibroblast | pluri_primary | d0→d10 | +2.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +5.000 | PartialReprog | pluri_primary | d0→d7 | +9.000 | +118.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +5.000 | PartialReprog | pluri_primary | d0→d10 | +9.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +6.000 | EarlyPluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +6.000 | EarlyPluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +7.000 | PartialReprog | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +7.000 | PartialReprog | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +8.000 | NonReprog | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +8.000 | NonReprog | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +9.000 | Fibroblast | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +9.000 | Fibroblast | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +10.000 | Pluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +10.000 | Pluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +11.000 | PartialReprog | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +11.000 | PartialReprog | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +12.000 | PartialReprog | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +12.000 | PartialReprog | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +13.000 | EarlyPluripotency | pluri_primary | d0→d7 | +56.000 | +164.000 | +1.778 | +8.758 | -6.980 | +0.990 | +198.000 | False | True | NA | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +13.000 | EarlyPluripotency | pluri_primary | d0→d10 | +56.000 | +432.000 | +1.778 | +6.874 | -5.096 | +0.970 | +194.000 | False | True | NA | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +14.000 | Fibroblast | pluri_primary | d0→d7 | +886.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +14.000 | Fibroblast | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +15.000 | EarlyPluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +15.000 | EarlyPluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +16.000 | Fibroblast | pluri_primary | d0→d7 | +361.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +16.000 | Fibroblast | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +17.000 | EarlyPluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +17.000 | EarlyPluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +18.000 | EarlyPluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +18.000 | EarlyPluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +0.000 | Fibroblast | md_score | d0→d7 | +6454.000 | +14.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +0.000 | Fibroblast | md_score | d0→d10 | +6454.000 | +2.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +1.000 | EarlyPluripotency | md_score | d0→d7 | +0.000 | +3240.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +1.000 | EarlyPluripotency | md_score | d0→d10 | +0.000 | +8.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +2.000 | EarlyPluripotency | md_score | d0→d7 | +0.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +2.000 | EarlyPluripotency | md_score | d0→d10 | +0.000 | +2730.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +3.000 | EarlyPluripotency | md_score | d0→d7 | +0.000 | +2537.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +3.000 | EarlyPluripotency | md_score | d0→d10 | +0.000 | +30.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +4.000 | Fibroblast | md_score | d0→d7 | +2.000 | +2210.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +4.000 | Fibroblast | md_score | d0→d10 | +2.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +5.000 | PartialReprog | md_score | d0→d7 | +9.000 | +118.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +5.000 | PartialReprog | md_score | d0→d10 | +9.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +6.000 | EarlyPluripotency | md_score | d0→d7 | +0.000 | +1819.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +6.000 | EarlyPluripotency | md_score | d0→d10 | +0.000 | +330.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +7.000 | PartialReprog | md_score | d0→d7 | +0.000 | +20.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +7.000 | PartialReprog | md_score | d0→d10 | +0.000 | +10.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +8.000 | NonReprog | md_score | d0→d7 | +0.000 | +8.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +8.000 | NonReprog | md_score | d0→d10 | +0.000 | +1794.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +9.000 | Fibroblast | md_score | d0→d7 | +0.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +9.000 | Fibroblast | md_score | d0→d10 | +0.000 | +1362.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +10.000 | Pluripotency | md_score | d0→d7 | +0.000 | +977.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +10.000 | Pluripotency | md_score | d0→d10 | +0.000 | +181.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +11.000 | PartialReprog | md_score | d0→d7 | +0.000 | +65.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +11.000 | PartialReprog | md_score | d0→d10 | +0.000 | +10.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +12.000 | PartialReprog | md_score | d0→d7 | +0.000 | +106.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +12.000 | PartialReprog | md_score | d0→d10 | +0.000 | +25.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +13.000 | EarlyPluripotency | md_score | d0→d7 | +56.000 | +164.000 | +0.379 | +0.159 | +0.220 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +13.000 | EarlyPluripotency | md_score | d0→d10 | +56.000 | +432.000 | +0.379 | +0.235 | +0.145 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +14.000 | Fibroblast | md_score | d0→d7 | +886.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +14.000 | Fibroblast | md_score | d0→d10 | +886.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +15.000 | EarlyPluripotency | md_score | d0→d7 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +15.000 | EarlyPluripotency | md_score | d0→d10 | +0.000 | +761.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +16.000 | Fibroblast | md_score | d0→d7 | +361.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +16.000 | Fibroblast | md_score | d0→d10 | +361.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +17.000 | EarlyPluripotency | md_score | d0→d7 | +0.000 | +122.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +17.000 | EarlyPluripotency | md_score | d0→d10 | +0.000 | +121.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +18.000 | EarlyPluripotency | md_score | d0→d7 | +0.000 | +289.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +18.000 | EarlyPluripotency | md_score | d0→d10 | +0.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +0.000 | Fibroblast | tgfb_score | d0→d7 | +6454.000 | +14.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +0.000 | Fibroblast | tgfb_score | d0→d10 | +6454.000 | +2.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +1.000 | EarlyPluripotency | tgfb_score | d0→d7 | +0.000 | +3240.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +1.000 | EarlyPluripotency | tgfb_score | d0→d10 | +0.000 | +8.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +2.000 | EarlyPluripotency | tgfb_score | d0→d7 | +0.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +2.000 | EarlyPluripotency | tgfb_score | d0→d10 | +0.000 | +2730.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +3.000 | EarlyPluripotency | tgfb_score | d0→d7 | +0.000 | +2537.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +3.000 | EarlyPluripotency | tgfb_score | d0→d10 | +0.000 | +30.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +4.000 | Fibroblast | tgfb_score | d0→d7 | +2.000 | +2210.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +4.000 | Fibroblast | tgfb_score | d0→d10 | +2.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +5.000 | PartialReprog | tgfb_score | d0→d7 | +9.000 | +118.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +5.000 | PartialReprog | tgfb_score | d0→d10 | +9.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +6.000 | EarlyPluripotency | tgfb_score | d0→d7 | +0.000 | +1819.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +6.000 | EarlyPluripotency | tgfb_score | d0→d10 | +0.000 | +330.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +7.000 | PartialReprog | tgfb_score | d0→d7 | +0.000 | +20.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +7.000 | PartialReprog | tgfb_score | d0→d10 | +0.000 | +10.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +8.000 | NonReprog | tgfb_score | d0→d7 | +0.000 | +8.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +8.000 | NonReprog | tgfb_score | d0→d10 | +0.000 | +1794.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +9.000 | Fibroblast | tgfb_score | d0→d7 | +0.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +9.000 | Fibroblast | tgfb_score | d0→d10 | +0.000 | +1362.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +10.000 | Pluripotency | tgfb_score | d0→d7 | +0.000 | +977.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +10.000 | Pluripotency | tgfb_score | d0→d10 | +0.000 | +181.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +11.000 | PartialReprog | tgfb_score | d0→d7 | +0.000 | +65.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +11.000 | PartialReprog | tgfb_score | d0→d10 | +0.000 | +10.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +12.000 | PartialReprog | tgfb_score | d0→d7 | +0.000 | +106.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +12.000 | PartialReprog | tgfb_score | d0→d10 | +0.000 | +25.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +13.000 | EarlyPluripotency | tgfb_score | d0→d7 | +56.000 | +164.000 | +0.077 | +0.056 | +0.022 | +0.373 | +74.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +13.000 | EarlyPluripotency | tgfb_score | d0→d10 | +56.000 | +432.000 | +0.077 | +0.071 | +0.006 | +0.403 | +80.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +14.000 | Fibroblast | tgfb_score | d0→d7 | +886.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +14.000 | Fibroblast | tgfb_score | d0→d10 | +886.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +15.000 | EarlyPluripotency | tgfb_score | d0→d7 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +15.000 | EarlyPluripotency | tgfb_score | d0→d10 | +0.000 | +761.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +16.000 | Fibroblast | tgfb_score | d0→d7 | +361.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +16.000 | Fibroblast | tgfb_score | d0→d10 | +361.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +17.000 | EarlyPluripotency | tgfb_score | d0→d7 | +0.000 | +122.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +17.000 | EarlyPluripotency | tgfb_score | d0→d10 | +0.000 | +121.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +18.000 | EarlyPluripotency | tgfb_score | d0→d7 | +0.000 | +289.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +18.000 | EarlyPluripotency | tgfb_score | d0→d10 | +0.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +0.000 | Fibroblast | age_up | d0→d7 | +6454.000 | +14.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +0.000 | Fibroblast | age_up | d0→d10 | +6454.000 | +2.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +1.000 | EarlyPluripotency | age_up | d0→d7 | +0.000 | +3240.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +1.000 | EarlyPluripotency | age_up | d0→d10 | +0.000 | +8.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +2.000 | EarlyPluripotency | age_up | d0→d7 | +0.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +2.000 | EarlyPluripotency | age_up | d0→d10 | +0.000 | +2730.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +3.000 | EarlyPluripotency | age_up | d0→d7 | +0.000 | +2537.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +3.000 | EarlyPluripotency | age_up | d0→d10 | +0.000 | +30.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +4.000 | Fibroblast | age_up | d0→d7 | +2.000 | +2210.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +4.000 | Fibroblast | age_up | d0→d10 | +2.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +5.000 | PartialReprog | age_up | d0→d7 | +9.000 | +118.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +5.000 | PartialReprog | age_up | d0→d10 | +9.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +6.000 | EarlyPluripotency | age_up | d0→d7 | +0.000 | +1819.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +6.000 | EarlyPluripotency | age_up | d0→d10 | +0.000 | +330.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +7.000 | PartialReprog | age_up | d0→d7 | +0.000 | +20.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +7.000 | PartialReprog | age_up | d0→d10 | +0.000 | +10.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +8.000 | NonReprog | age_up | d0→d7 | +0.000 | +8.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +8.000 | NonReprog | age_up | d0→d10 | +0.000 | +1794.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +9.000 | Fibroblast | age_up | d0→d7 | +0.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +9.000 | Fibroblast | age_up | d0→d10 | +0.000 | +1362.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +10.000 | Pluripotency | age_up | d0→d7 | +0.000 | +977.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +10.000 | Pluripotency | age_up | d0→d10 | +0.000 | +181.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +11.000 | PartialReprog | age_up | d0→d7 | +0.000 | +65.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +11.000 | PartialReprog | age_up | d0→d10 | +0.000 | +10.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +12.000 | PartialReprog | age_up | d0→d7 | +0.000 | +106.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +12.000 | PartialReprog | age_up | d0→d10 | +0.000 | +25.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +13.000 | EarlyPluripotency | age_up | d0→d7 | +56.000 | +164.000 | +0.083 | +0.035 | +0.048 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +13.000 | EarlyPluripotency | age_up | d0→d10 | +56.000 | +432.000 | +0.083 | +0.046 | +0.037 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +14.000 | Fibroblast | age_up | d0→d7 | +886.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +14.000 | Fibroblast | age_up | d0→d10 | +886.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +15.000 | EarlyPluripotency | age_up | d0→d7 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +15.000 | EarlyPluripotency | age_up | d0→d10 | +0.000 | +761.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +16.000 | Fibroblast | age_up | d0→d7 | +361.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +16.000 | Fibroblast | age_up | d0→d10 | +361.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +17.000 | EarlyPluripotency | age_up | d0→d7 | +0.000 | +122.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +17.000 | EarlyPluripotency | age_up | d0→d10 | +0.000 | +121.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +18.000 | EarlyPluripotency | age_up | d0→d7 | +0.000 | +289.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +18.000 | EarlyPluripotency | age_up | d0→d10 | +0.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +0.000 | Fibroblast | age_down | d0→d7 | +6454.000 | +14.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +0.000 | Fibroblast | age_down | d0→d10 | +6454.000 | +2.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +1.000 | EarlyPluripotency | age_down | d0→d7 | +0.000 | +3240.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +1.000 | EarlyPluripotency | age_down | d0→d10 | +0.000 | +8.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +2.000 | EarlyPluripotency | age_down | d0→d7 | +0.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +2.000 | EarlyPluripotency | age_down | d0→d10 | +0.000 | +2730.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +3.000 | EarlyPluripotency | age_down | d0→d7 | +0.000 | +2537.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +3.000 | EarlyPluripotency | age_down | d0→d10 | +0.000 | +30.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +4.000 | Fibroblast | age_down | d0→d7 | +2.000 | +2210.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +4.000 | Fibroblast | age_down | d0→d10 | +2.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +5.000 | PartialReprog | age_down | d0→d7 | +9.000 | +118.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +5.000 | PartialReprog | age_down | d0→d10 | +9.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +6.000 | EarlyPluripotency | age_down | d0→d7 | +0.000 | +1819.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +6.000 | EarlyPluripotency | age_down | d0→d10 | +0.000 | +330.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +7.000 | PartialReprog | age_down | d0→d7 | +0.000 | +20.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +7.000 | PartialReprog | age_down | d0→d10 | +0.000 | +10.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +8.000 | NonReprog | age_down | d0→d7 | +0.000 | +8.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +8.000 | NonReprog | age_down | d0→d10 | +0.000 | +1794.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +9.000 | Fibroblast | age_down | d0→d7 | +0.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +9.000 | Fibroblast | age_down | d0→d10 | +0.000 | +1362.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +10.000 | Pluripotency | age_down | d0→d7 | +0.000 | +977.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +10.000 | Pluripotency | age_down | d0→d10 | +0.000 | +181.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +11.000 | PartialReprog | age_down | d0→d7 | +0.000 | +65.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +11.000 | PartialReprog | age_down | d0→d10 | +0.000 | +10.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +12.000 | PartialReprog | age_down | d0→d7 | +0.000 | +106.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +12.000 | PartialReprog | age_down | d0→d10 | +0.000 | +25.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +13.000 | EarlyPluripotency | age_down | d0→d7 | +56.000 | +164.000 | -0.038 | +0.001 | -0.039 | +1.000 | +200.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +13.000 | EarlyPluripotency | age_down | d0→d10 | +56.000 | +432.000 | -0.038 | -0.017 | -0.020 | +1.000 | +200.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +14.000 | Fibroblast | age_down | d0→d7 | +886.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +14.000 | Fibroblast | age_down | d0→d10 | +886.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +15.000 | EarlyPluripotency | age_down | d0→d7 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +15.000 | EarlyPluripotency | age_down | d0→d10 | +0.000 | +761.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +16.000 | Fibroblast | age_down | d0→d7 | +361.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +16.000 | Fibroblast | age_down | d0→d10 | +361.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +17.000 | EarlyPluripotency | age_down | d0→d7 | +0.000 | +122.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +17.000 | EarlyPluripotency | age_down | d0→d10 | +0.000 | +121.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +18.000 | EarlyPluripotency | age_down | d0→d7 | +0.000 | +289.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +18.000 | EarlyPluripotency | age_down | d0→d10 | +0.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |

### (c) per cluster at ≥10 cells — sensitivity check, not primary

| cell_line | cluster | label | day | n_cells | below_min_cells | age_score | md_score | tgfb_score | age_up | age_down | pluri_primary |
|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | 0 | Fibroblast | 0 | 3998 | False | +7.037 | +0.513 | +0.040 | +0.058 | +0.018 | -0.743 |
| GM00731 | 3 | PartialReprog | 0 | 1 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 4 | NonReprog | 0 | 9 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 5 | PartialReprog | 0 | 2 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 10 | NonReprog | 0 | 56 | False | -1.425 | +0.356 | +0.117 | +0.112 | -0.052 | +1.860 |
| GM00731 | 12 | Fibroblast | 0 | 569 | False | +7.865 | +0.510 | +0.049 | +0.062 | +0.032 | -1.198 |
| GM00731 | 15 | Fibroblast | 0 | 372 | False | +4.533 | +0.428 | +0.050 | +0.009 | +0.097 | +1.003 |
| GM00731 | 16 | EarlyPluripotency | 0 | 2 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 0 | Fibroblast | 3 | 4 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 3 | PartialReprog | 3 | 1589 | False | +6.055 | +0.308 | +0.038 | -0.024 | +0.097 | +5.941 |
| GM00731 | 4 | NonReprog | 3 | 259 | False | +9.226 | +0.531 | +0.060 | +0.030 | +0.031 | +3.900 |
| GM00731 | 5 | PartialReprog | 3 | 1523 | False | +4.013 | +0.060 | -0.096 | -0.015 | +0.071 | +9.042 |
| GM00731 | 7 | EarlyPluripotency | 3 | 108 | False | +5.193 | +0.163 | -0.047 | +0.012 | +0.051 | +7.227 |
| GM00731 | 8 | PartialReprog | 3 | 973 | False | -0.379 | +0.020 | -0.098 | -0.013 | +0.081 | +8.942 |
| GM00731 | 10 | NonReprog | 3 | 89 | False | +7.215 | +0.316 | +0.117 | +0.024 | +0.011 | +5.974 |
| GM00731 | 11 | PartialReprog | 3 | 557 | False | +3.280 | +0.080 | -0.089 | -0.016 | +0.068 | +7.831 |
| GM00731 | 13 | PartialReprog | 3 | 442 | False | +4.836 | +0.079 | -0.031 | -0.024 | +0.083 | +9.860 |
| GM00731 | 14 | Pluripotency | 3 | 5 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 15 | Fibroblast | 3 | 1 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 16 | EarlyPluripotency | 3 | 302 | False | +0.353 | +0.009 | -0.015 | +0.006 | +0.002 | +12.562 |
| GM00731 | 18 | PartialReprog | 3 | 136 | False | +3.692 | +0.171 | -0.016 | -0.020 | +0.065 | +7.768 |
| GM00731 | 0 | Fibroblast | 7 | 9 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 1 | NonReprog | 7 | 9 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 2 | EarlyPluripotency | 7 | 2092 | False | +0.098 | +0.167 | -0.004 | +0.033 | +0.029 | +8.439 |
| GM00731 | 3 | PartialReprog | 7 | 444 | False | +6.644 | +0.343 | +0.006 | +0.008 | +0.090 | +4.785 |
| GM00731 | 4 | NonReprog | 7 | 1732 | False | +8.658 | +0.456 | +0.045 | +0.046 | +0.024 | +3.152 |
| GM00731 | 5 | PartialReprog | 7 | 48 | False | +0.508 | +0.052 | -0.088 | +0.030 | +0.020 | +11.007 |
| GM00731 | 6 | Fibroblast | 7 | 1 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 7 | EarlyPluripotency | 7 | 884 | False | +4.745 | +0.299 | +0.001 | +0.043 | +0.042 | +7.337 |
| GM00731 | 8 | PartialReprog | 7 | 13 | False | -0.156 | +0.071 | -0.099 | +0.002 | +0.075 | +7.183 |
| GM00731 | 9 | EarlyPluripotency | 7 | 909 | False | -2.000 | -0.024 | -0.106 | +0.023 | +0.041 | +11.825 |
| GM00731 | 10 | NonReprog | 7 | 191 | False | +3.966 | +0.347 | +0.079 | +0.053 | -0.017 | +5.835 |
| GM00731 | 11 | PartialReprog | 7 | 64 | False | +4.110 | +0.117 | -0.064 | +0.017 | +0.043 | +6.920 |
| GM00731 | 12 | Fibroblast | 7 | 1 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 13 | PartialReprog | 7 | 3 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 14 | Pluripotency | 7 | 354 | False | -5.699 | +0.011 | +0.005 | +0.056 | -0.022 | +11.693 |
| GM00731 | 15 | Fibroblast | 7 | 3 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 16 | EarlyPluripotency | 7 | 32 | False | +1.891 | +0.188 | +0.033 | +0.038 | +0.004 | +4.024 |
| GM00731 | 18 | PartialReprog | 7 | 6 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 0 | Fibroblast | 10 | 1 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 1 | NonReprog | 10 | 2313 | False | +11.686 | +0.446 | -0.018 | +0.058 | +0.027 | +1.788 |
| GM00731 | 2 | EarlyPluripotency | 10 | 4 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 3 | PartialReprog | 10 | 5 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 4 | NonReprog | 10 | 31 | False | +13.645 | +0.408 | -0.045 | +0.028 | +0.026 | +2.259 |
| GM00731 | 5 | PartialReprog | 10 | 21 | False | +2.744 | +0.320 | +0.008 | +0.007 | +0.023 | +6.696 |
| GM00731 | 6 | Fibroblast | 10 | 1029 | False | +9.102 | +0.434 | -0.009 | +0.020 | +0.049 | +2.984 |
| GM00731 | 7 | EarlyPluripotency | 10 | 22 | False | +6.501 | +0.289 | -0.071 | +0.052 | +0.031 | +6.230 |
| GM00731 | 9 | EarlyPluripotency | 10 | 52 | False | -0.290 | -0.040 | -0.109 | +0.001 | +0.067 | +12.005 |
| GM00731 | 10 | NonReprog | 10 | 527 | False | +7.934 | +0.417 | +0.102 | +0.063 | -0.023 | +3.313 |
| GM00731 | 11 | PartialReprog | 10 | 5 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 14 | Pluripotency | 10 | 80 | False | -2.372 | +0.069 | +0.019 | +0.051 | -0.019 | +10.355 |
| GM00731 | 16 | EarlyPluripotency | 10 | 7 | True | NA | NA | NA | NA | NA | NA |
| GM00731 | 17 | NonReprog | 10 | 243 | False | +8.127 | +0.356 | -0.036 | +0.016 | +0.092 | +1.015 |
| GM23815 | 0 | Fibroblast | 0 | 6454 | False | +3.253 | +0.483 | +0.033 | +0.042 | +0.021 | -0.154 |
| GM23815 | 4 | Fibroblast | 0 | 2 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 5 | PartialReprog | 0 | 9 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 13 | EarlyPluripotency | 0 | 56 | False | -2.160 | +0.379 | +0.077 | +0.083 | -0.038 | +1.778 |
| GM23815 | 14 | Fibroblast | 0 | 886 | False | +1.849 | +0.420 | +0.016 | +0.005 | +0.096 | +0.918 |
| GM23815 | 16 | Fibroblast | 0 | 361 | False | +3.401 | +0.547 | -0.005 | +0.057 | +0.010 | -0.696 |
| GM23815 | 0 | Fibroblast | 3 | 19 | False | +6.573 | +0.400 | +0.126 | -0.008 | +0.051 | +3.086 |
| GM23815 | 1 | EarlyPluripotency | 3 | 3 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 3 | EarlyPluripotency | 3 | 3 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 4 | Fibroblast | 3 | 21 | False | +8.218 | +0.391 | +0.045 | +0.006 | +0.045 | +3.558 |
| GM23815 | 5 | PartialReprog | 3 | 2039 | False | +1.227 | +0.219 | +0.029 | -0.030 | +0.112 | +6.379 |
| GM23815 | 7 | PartialReprog | 3 | 2039 | False | -3.618 | -0.011 | -0.084 | -0.017 | +0.091 | +9.994 |
| GM23815 | 8 | NonReprog | 3 | 10 | False | +2.399 | +0.380 | +0.051 | +0.020 | +0.055 | +5.793 |
| GM23815 | 10 | Pluripotency | 3 | 91 | False | -8.931 | +0.028 | -0.012 | +0.021 | +0.007 | +10.539 |
| GM23815 | 11 | PartialReprog | 3 | 1065 | False | -0.322 | +0.007 | -0.065 | -0.018 | +0.078 | +10.360 |
| GM23815 | 12 | PartialReprog | 3 | 970 | False | +0.720 | +0.095 | -0.043 | -0.016 | +0.070 | +7.424 |
| GM23815 | 13 | EarlyPluripotency | 3 | 246 | False | +2.602 | +0.034 | +0.023 | +0.006 | +0.008 | +9.709 |
| GM23815 | 14 | Fibroblast | 3 | 5 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 16 | Fibroblast | 3 | 1 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 17 | EarlyPluripotency | 3 | 64 | False | -1.597 | +0.156 | -0.011 | +0.001 | +0.072 | +7.501 |
| GM23815 | 0 | Fibroblast | 7 | 14 | False | +1.026 | +0.401 | +0.004 | +0.023 | +0.027 | +1.495 |
| GM23815 | 1 | EarlyPluripotency | 7 | 3240 | False | -3.860 | +0.207 | -0.041 | +0.028 | +0.035 | +7.946 |
| GM23815 | 2 | EarlyPluripotency | 7 | 1 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 3 | EarlyPluripotency | 7 | 2537 | False | -6.390 | +0.003 | -0.069 | +0.017 | +0.042 | +12.127 |
| GM23815 | 4 | Fibroblast | 7 | 2210 | False | +1.091 | +0.410 | -0.015 | +0.033 | +0.024 | +4.365 |
| GM23815 | 5 | PartialReprog | 7 | 118 | False | +1.778 | +0.314 | -0.018 | -0.003 | +0.083 | +4.847 |
| GM23815 | 6 | EarlyPluripotency | 7 | 1819 | False | -5.449 | -0.048 | -0.083 | -0.006 | +0.072 | +12.685 |
| GM23815 | 7 | PartialReprog | 7 | 20 | False | -6.363 | +0.006 | -0.107 | +0.024 | +0.046 | +11.332 |
| GM23815 | 8 | NonReprog | 7 | 8 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 9 | Fibroblast | 7 | 3 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 10 | Pluripotency | 7 | 977 | False | -8.857 | -0.059 | +0.006 | +0.027 | -0.004 | +12.656 |
| GM23815 | 11 | PartialReprog | 7 | 65 | False | -2.332 | +0.023 | -0.101 | +0.021 | +0.037 | +10.838 |
| GM23815 | 12 | PartialReprog | 7 | 106 | False | +0.145 | +0.143 | -0.067 | +0.011 | +0.039 | +7.850 |
| GM23815 | 13 | EarlyPluripotency | 7 | 164 | False | -7.076 | +0.159 | +0.056 | +0.035 | +0.001 | +8.758 |
| GM23815 | 14 | Fibroblast | 7 | 5 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 16 | Fibroblast | 7 | 1 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 17 | EarlyPluripotency | 7 | 122 | False | -2.776 | +0.083 | -0.072 | +0.041 | +0.035 | +9.561 |
| GM23815 | 18 | EarlyPluripotency | 7 | 289 | False | -7.347 | +0.001 | -0.055 | +0.029 | +0.025 | +11.933 |
| GM23815 | 0 | Fibroblast | 10 | 2 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 1 | EarlyPluripotency | 10 | 8 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 2 | EarlyPluripotency | 10 | 2730 | False | -1.570 | +0.155 | -0.037 | +0.005 | +0.052 | +7.637 |
| GM23815 | 3 | EarlyPluripotency | 10 | 30 | False | -3.984 | -0.036 | -0.096 | +0.002 | +0.041 | +11.657 |
| GM23815 | 4 | Fibroblast | 10 | 9 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 5 | PartialReprog | 10 | 5 | True | NA | NA | NA | NA | NA | NA |
| GM23815 | 6 | EarlyPluripotency | 10 | 330 | False | -4.946 | -0.115 | -0.078 | -0.024 | +0.081 | +14.675 |
| GM23815 | 7 | PartialReprog | 10 | 10 | False | +0.171 | -0.040 | -0.129 | -0.012 | +0.068 | +9.305 |
| GM23815 | 8 | NonReprog | 10 | 1794 | False | +5.408 | +0.423 | +0.035 | +0.046 | +0.026 | +2.530 |
| GM23815 | 9 | Fibroblast | 10 | 1362 | False | +3.204 | +0.380 | +0.036 | +0.012 | +0.047 | +3.122 |
| GM23815 | 10 | Pluripotency | 10 | 181 | False | -6.660 | -0.009 | +0.008 | +0.018 | +0.004 | +11.955 |
| GM23815 | 11 | PartialReprog | 10 | 10 | False | -8.761 | +0.042 | -0.077 | +0.026 | +0.012 | +11.799 |
| GM23815 | 12 | PartialReprog | 10 | 25 | False | +3.381 | +0.194 | +0.007 | -0.014 | +0.073 | +4.677 |
| GM23815 | 13 | EarlyPluripotency | 10 | 432 | False | -3.794 | +0.235 | +0.071 | +0.046 | -0.017 | +6.874 |
| GM23815 | 15 | EarlyPluripotency | 10 | 761 | False | +1.379 | +0.278 | -0.002 | +0.016 | +0.066 | +7.499 |
| GM23815 | 17 | EarlyPluripotency | 10 | 121 | False | +1.877 | +0.299 | +0.001 | +0.028 | +0.044 | +6.431 |
| GM23815 | 18 | EarlyPluripotency | 10 | 5 | True | NA | NA | NA | NA | NA | NA |

| cell_line | cluster | label | instrument | endpoint | n_cells_d0 | n_cells_end | score_d0 | score_end | decline | p | n_random_ge_real | pass_p | ok | reason | null_kind |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | +0.000 | Fibroblast | frozen_ruler | d0→d7 | +3998.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | permuted_ruler_weights |
| GM00731 | +0.000 | Fibroblast | frozen_ruler | d0→d10 | +3998.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | permuted_ruler_weights |
| GM00731 | +1.000 | NonReprog | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +1.000 | NonReprog | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +2.000 | EarlyPluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +2.000 | EarlyPluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +3.000 | PartialReprog | frozen_ruler | d0→d7 | +1.000 | +444.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | permuted_ruler_weights |
| GM00731 | +3.000 | PartialReprog | frozen_ruler | d0→d10 | +1.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | permuted_ruler_weights |
| GM00731 | +4.000 | NonReprog | frozen_ruler | d0→d7 | +9.000 | +1732.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | permuted_ruler_weights |
| GM00731 | +4.000 | NonReprog | frozen_ruler | d0→d10 | +9.000 | +31.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | permuted_ruler_weights |
| GM00731 | +5.000 | PartialReprog | frozen_ruler | d0→d7 | +2.000 | +48.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | permuted_ruler_weights |
| GM00731 | +5.000 | PartialReprog | frozen_ruler | d0→d10 | +2.000 | +21.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | permuted_ruler_weights |
| GM00731 | +6.000 | Fibroblast | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +6.000 | Fibroblast | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +7.000 | EarlyPluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +7.000 | EarlyPluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +8.000 | PartialReprog | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +8.000 | PartialReprog | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +9.000 | EarlyPluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +9.000 | EarlyPluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +10.000 | NonReprog | frozen_ruler | d0→d7 | +56.000 | +191.000 | -1.425 | +3.966 | -5.391 | +0.930 | +186.000 | False | True | NA | permuted_ruler_weights |
| GM00731 | +10.000 | NonReprog | frozen_ruler | d0→d10 | +56.000 | +527.000 | -1.425 | +7.934 | -9.359 | +1.000 | +200.000 | False | True | NA | permuted_ruler_weights |
| GM00731 | +11.000 | PartialReprog | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +11.000 | PartialReprog | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +12.000 | Fibroblast | frozen_ruler | d0→d7 | +569.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | permuted_ruler_weights |
| GM00731 | +12.000 | Fibroblast | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +13.000 | PartialReprog | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +13.000 | PartialReprog | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +14.000 | Pluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +14.000 | Pluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +15.000 | Fibroblast | frozen_ruler | d0→d7 | +372.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | permuted_ruler_weights |
| GM00731 | +15.000 | Fibroblast | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +16.000 | EarlyPluripotency | frozen_ruler | d0→d7 | +2.000 | +32.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | permuted_ruler_weights |
| GM00731 | +16.000 | EarlyPluripotency | frozen_ruler | d0→d10 | +2.000 | +7.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | permuted_ruler_weights |
| GM00731 | +17.000 | NonReprog | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +17.000 | NonReprog | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +18.000 | PartialReprog | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +18.000 | PartialReprog | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM00731 | +0.000 | Fibroblast | pluri_primary | d0→d7 | +3998.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +0.000 | Fibroblast | pluri_primary | d0→d10 | +3998.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +1.000 | NonReprog | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +1.000 | NonReprog | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +2.000 | EarlyPluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +2.000 | EarlyPluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +3.000 | PartialReprog | pluri_primary | d0→d7 | +1.000 | +444.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +3.000 | PartialReprog | pluri_primary | d0→d10 | +1.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +4.000 | NonReprog | pluri_primary | d0→d7 | +9.000 | +1732.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +4.000 | NonReprog | pluri_primary | d0→d10 | +9.000 | +31.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +5.000 | PartialReprog | pluri_primary | d0→d7 | +2.000 | +48.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +5.000 | PartialReprog | pluri_primary | d0→d10 | +2.000 | +21.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +6.000 | Fibroblast | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +6.000 | Fibroblast | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +7.000 | EarlyPluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +7.000 | EarlyPluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +8.000 | PartialReprog | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +8.000 | PartialReprog | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +9.000 | EarlyPluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +9.000 | EarlyPluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +10.000 | NonReprog | pluri_primary | d0→d7 | +56.000 | +191.000 | +1.860 | +5.835 | -3.975 | +0.970 | +194.000 | False | True | NA | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +10.000 | NonReprog | pluri_primary | d0→d10 | +56.000 | +527.000 | +1.860 | +3.313 | -1.453 | +0.831 | +166.000 | False | True | NA | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +11.000 | PartialReprog | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +11.000 | PartialReprog | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +12.000 | Fibroblast | pluri_primary | d0→d7 | +569.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +12.000 | Fibroblast | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +13.000 | PartialReprog | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +13.000 | PartialReprog | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +14.000 | Pluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +14.000 | Pluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +15.000 | Fibroblast | pluri_primary | d0→d7 | +372.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +15.000 | Fibroblast | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +16.000 | EarlyPluripotency | pluri_primary | d0→d7 | +2.000 | +32.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +16.000 | EarlyPluripotency | pluri_primary | d0→d10 | +2.000 | +7.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +17.000 | NonReprog | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +17.000 | NonReprog | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +18.000 | PartialReprog | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +18.000 | PartialReprog | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM00731 | +0.000 | Fibroblast | md_score | d0→d7 | +3998.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +0.000 | Fibroblast | md_score | d0→d10 | +3998.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +1.000 | NonReprog | md_score | d0→d7 | +0.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +1.000 | NonReprog | md_score | d0→d10 | +0.000 | +2313.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +2.000 | EarlyPluripotency | md_score | d0→d7 | +0.000 | +2092.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +2.000 | EarlyPluripotency | md_score | d0→d10 | +0.000 | +4.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +3.000 | PartialReprog | md_score | d0→d7 | +1.000 | +444.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +3.000 | PartialReprog | md_score | d0→d10 | +1.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +4.000 | NonReprog | md_score | d0→d7 | +9.000 | +1732.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +4.000 | NonReprog | md_score | d0→d10 | +9.000 | +31.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +5.000 | PartialReprog | md_score | d0→d7 | +2.000 | +48.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +5.000 | PartialReprog | md_score | d0→d10 | +2.000 | +21.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +6.000 | Fibroblast | md_score | d0→d7 | +0.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +6.000 | Fibroblast | md_score | d0→d10 | +0.000 | +1029.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +7.000 | EarlyPluripotency | md_score | d0→d7 | +0.000 | +884.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +7.000 | EarlyPluripotency | md_score | d0→d10 | +0.000 | +22.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +8.000 | PartialReprog | md_score | d0→d7 | +0.000 | +13.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +8.000 | PartialReprog | md_score | d0→d10 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +9.000 | EarlyPluripotency | md_score | d0→d7 | +0.000 | +909.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +9.000 | EarlyPluripotency | md_score | d0→d10 | +0.000 | +52.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +10.000 | NonReprog | md_score | d0→d7 | +56.000 | +191.000 | +0.356 | +0.347 | +0.009 | +0.035 | +6.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +10.000 | NonReprog | md_score | d0→d10 | +56.000 | +527.000 | +0.356 | +0.417 | -0.061 | +0.886 | +177.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +11.000 | PartialReprog | md_score | d0→d7 | +0.000 | +64.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +11.000 | PartialReprog | md_score | d0→d10 | +0.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +12.000 | Fibroblast | md_score | d0→d7 | +569.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +12.000 | Fibroblast | md_score | d0→d10 | +569.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +13.000 | PartialReprog | md_score | d0→d7 | +0.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +13.000 | PartialReprog | md_score | d0→d10 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +14.000 | Pluripotency | md_score | d0→d7 | +0.000 | +354.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +14.000 | Pluripotency | md_score | d0→d10 | +0.000 | +80.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +15.000 | Fibroblast | md_score | d0→d7 | +372.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +15.000 | Fibroblast | md_score | d0→d10 | +372.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +16.000 | EarlyPluripotency | md_score | d0→d7 | +2.000 | +32.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +16.000 | EarlyPluripotency | md_score | d0→d10 | +2.000 | +7.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +17.000 | NonReprog | md_score | d0→d7 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +17.000 | NonReprog | md_score | d0→d10 | +0.000 | +243.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +18.000 | PartialReprog | md_score | d0→d7 | +0.000 | +6.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +18.000 | PartialReprog | md_score | d0→d10 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +0.000 | Fibroblast | tgfb_score | d0→d7 | +3998.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +0.000 | Fibroblast | tgfb_score | d0→d10 | +3998.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +1.000 | NonReprog | tgfb_score | d0→d7 | +0.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +1.000 | NonReprog | tgfb_score | d0→d10 | +0.000 | +2313.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +2.000 | EarlyPluripotency | tgfb_score | d0→d7 | +0.000 | +2092.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +2.000 | EarlyPluripotency | tgfb_score | d0→d10 | +0.000 | +4.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +3.000 | PartialReprog | tgfb_score | d0→d7 | +1.000 | +444.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +3.000 | PartialReprog | tgfb_score | d0→d10 | +1.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +4.000 | NonReprog | tgfb_score | d0→d7 | +9.000 | +1732.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +4.000 | NonReprog | tgfb_score | d0→d10 | +9.000 | +31.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +5.000 | PartialReprog | tgfb_score | d0→d7 | +2.000 | +48.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +5.000 | PartialReprog | tgfb_score | d0→d10 | +2.000 | +21.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +6.000 | Fibroblast | tgfb_score | d0→d7 | +0.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +6.000 | Fibroblast | tgfb_score | d0→d10 | +0.000 | +1029.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +7.000 | EarlyPluripotency | tgfb_score | d0→d7 | +0.000 | +884.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +7.000 | EarlyPluripotency | tgfb_score | d0→d10 | +0.000 | +22.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +8.000 | PartialReprog | tgfb_score | d0→d7 | +0.000 | +13.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +8.000 | PartialReprog | tgfb_score | d0→d10 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +9.000 | EarlyPluripotency | tgfb_score | d0→d7 | +0.000 | +909.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +9.000 | EarlyPluripotency | tgfb_score | d0→d10 | +0.000 | +52.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +10.000 | NonReprog | tgfb_score | d0→d7 | +56.000 | +191.000 | +0.117 | +0.079 | +0.038 | +0.119 | +23.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +10.000 | NonReprog | tgfb_score | d0→d10 | +56.000 | +527.000 | +0.117 | +0.102 | +0.015 | +0.294 | +58.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +11.000 | PartialReprog | tgfb_score | d0→d7 | +0.000 | +64.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +11.000 | PartialReprog | tgfb_score | d0→d10 | +0.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +12.000 | Fibroblast | tgfb_score | d0→d7 | +569.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +12.000 | Fibroblast | tgfb_score | d0→d10 | +569.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +13.000 | PartialReprog | tgfb_score | d0→d7 | +0.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +13.000 | PartialReprog | tgfb_score | d0→d10 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +14.000 | Pluripotency | tgfb_score | d0→d7 | +0.000 | +354.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +14.000 | Pluripotency | tgfb_score | d0→d10 | +0.000 | +80.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +15.000 | Fibroblast | tgfb_score | d0→d7 | +372.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +15.000 | Fibroblast | tgfb_score | d0→d10 | +372.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +16.000 | EarlyPluripotency | tgfb_score | d0→d7 | +2.000 | +32.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +16.000 | EarlyPluripotency | tgfb_score | d0→d10 | +2.000 | +7.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +17.000 | NonReprog | tgfb_score | d0→d7 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +17.000 | NonReprog | tgfb_score | d0→d10 | +0.000 | +243.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +18.000 | PartialReprog | tgfb_score | d0→d7 | +0.000 | +6.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +18.000 | PartialReprog | tgfb_score | d0→d10 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +0.000 | Fibroblast | age_up | d0→d7 | +3998.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +0.000 | Fibroblast | age_up | d0→d10 | +3998.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +1.000 | NonReprog | age_up | d0→d7 | +0.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +1.000 | NonReprog | age_up | d0→d10 | +0.000 | +2313.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +2.000 | EarlyPluripotency | age_up | d0→d7 | +0.000 | +2092.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +2.000 | EarlyPluripotency | age_up | d0→d10 | +0.000 | +4.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +3.000 | PartialReprog | age_up | d0→d7 | +1.000 | +444.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +3.000 | PartialReprog | age_up | d0→d10 | +1.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +4.000 | NonReprog | age_up | d0→d7 | +9.000 | +1732.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +4.000 | NonReprog | age_up | d0→d10 | +9.000 | +31.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +5.000 | PartialReprog | age_up | d0→d7 | +2.000 | +48.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +5.000 | PartialReprog | age_up | d0→d10 | +2.000 | +21.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +6.000 | Fibroblast | age_up | d0→d7 | +0.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +6.000 | Fibroblast | age_up | d0→d10 | +0.000 | +1029.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +7.000 | EarlyPluripotency | age_up | d0→d7 | +0.000 | +884.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +7.000 | EarlyPluripotency | age_up | d0→d10 | +0.000 | +22.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +8.000 | PartialReprog | age_up | d0→d7 | +0.000 | +13.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +8.000 | PartialReprog | age_up | d0→d10 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +9.000 | EarlyPluripotency | age_up | d0→d7 | +0.000 | +909.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +9.000 | EarlyPluripotency | age_up | d0→d10 | +0.000 | +52.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +10.000 | NonReprog | age_up | d0→d7 | +56.000 | +191.000 | +0.112 | +0.053 | +0.059 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +10.000 | NonReprog | age_up | d0→d10 | +56.000 | +527.000 | +0.112 | +0.063 | +0.049 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +11.000 | PartialReprog | age_up | d0→d7 | +0.000 | +64.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +11.000 | PartialReprog | age_up | d0→d10 | +0.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +12.000 | Fibroblast | age_up | d0→d7 | +569.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +12.000 | Fibroblast | age_up | d0→d10 | +569.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +13.000 | PartialReprog | age_up | d0→d7 | +0.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +13.000 | PartialReprog | age_up | d0→d10 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +14.000 | Pluripotency | age_up | d0→d7 | +0.000 | +354.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +14.000 | Pluripotency | age_up | d0→d10 | +0.000 | +80.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +15.000 | Fibroblast | age_up | d0→d7 | +372.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +15.000 | Fibroblast | age_up | d0→d10 | +372.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +16.000 | EarlyPluripotency | age_up | d0→d7 | +2.000 | +32.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +16.000 | EarlyPluripotency | age_up | d0→d10 | +2.000 | +7.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +17.000 | NonReprog | age_up | d0→d7 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +17.000 | NonReprog | age_up | d0→d10 | +0.000 | +243.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +18.000 | PartialReprog | age_up | d0→d7 | +0.000 | +6.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +18.000 | PartialReprog | age_up | d0→d10 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +0.000 | Fibroblast | age_down | d0→d7 | +3998.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +0.000 | Fibroblast | age_down | d0→d10 | +3998.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +1.000 | NonReprog | age_down | d0→d7 | +0.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +1.000 | NonReprog | age_down | d0→d10 | +0.000 | +2313.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +2.000 | EarlyPluripotency | age_down | d0→d7 | +0.000 | +2092.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +2.000 | EarlyPluripotency | age_down | d0→d10 | +0.000 | +4.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +3.000 | PartialReprog | age_down | d0→d7 | +1.000 | +444.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +3.000 | PartialReprog | age_down | d0→d10 | +1.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +4.000 | NonReprog | age_down | d0→d7 | +9.000 | +1732.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +4.000 | NonReprog | age_down | d0→d10 | +9.000 | +31.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +5.000 | PartialReprog | age_down | d0→d7 | +2.000 | +48.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +5.000 | PartialReprog | age_down | d0→d10 | +2.000 | +21.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +6.000 | Fibroblast | age_down | d0→d7 | +0.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +6.000 | Fibroblast | age_down | d0→d10 | +0.000 | +1029.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +7.000 | EarlyPluripotency | age_down | d0→d7 | +0.000 | +884.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +7.000 | EarlyPluripotency | age_down | d0→d10 | +0.000 | +22.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +8.000 | PartialReprog | age_down | d0→d7 | +0.000 | +13.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +8.000 | PartialReprog | age_down | d0→d10 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +9.000 | EarlyPluripotency | age_down | d0→d7 | +0.000 | +909.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +9.000 | EarlyPluripotency | age_down | d0→d10 | +0.000 | +52.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +10.000 | NonReprog | age_down | d0→d7 | +56.000 | +191.000 | -0.052 | -0.017 | -0.035 | +1.000 | +200.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +10.000 | NonReprog | age_down | d0→d10 | +56.000 | +527.000 | -0.052 | -0.023 | -0.028 | +1.000 | +200.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +11.000 | PartialReprog | age_down | d0→d7 | +0.000 | +64.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +11.000 | PartialReprog | age_down | d0→d10 | +0.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +12.000 | Fibroblast | age_down | d0→d7 | +569.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +12.000 | Fibroblast | age_down | d0→d10 | +569.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +13.000 | PartialReprog | age_down | d0→d7 | +0.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +13.000 | PartialReprog | age_down | d0→d10 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +14.000 | Pluripotency | age_down | d0→d7 | +0.000 | +354.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +14.000 | Pluripotency | age_down | d0→d10 | +0.000 | +80.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +15.000 | Fibroblast | age_down | d0→d7 | +372.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +15.000 | Fibroblast | age_down | d0→d10 | +372.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +16.000 | EarlyPluripotency | age_down | d0→d7 | +2.000 | +32.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +16.000 | EarlyPluripotency | age_down | d0→d10 | +2.000 | +7.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +17.000 | NonReprog | age_down | d0→d7 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +17.000 | NonReprog | age_down | d0→d10 | +0.000 | +243.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +18.000 | PartialReprog | age_down | d0→d7 | +0.000 | +6.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM00731 | +18.000 | PartialReprog | age_down | d0→d10 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +0.000 | Fibroblast | frozen_ruler | d0→d7 | +6454.000 | +14.000 | +3.253 | +1.026 | +2.227 | +0.517 | +103.000 | False | True | NA | permuted_ruler_weights |
| GM23815 | +0.000 | Fibroblast | frozen_ruler | d0→d10 | +6454.000 | +2.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | permuted_ruler_weights |
| GM23815 | +1.000 | EarlyPluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +1.000 | EarlyPluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +2.000 | EarlyPluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +2.000 | EarlyPluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +3.000 | EarlyPluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +3.000 | EarlyPluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +4.000 | Fibroblast | frozen_ruler | d0→d7 | +2.000 | +2210.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | permuted_ruler_weights |
| GM23815 | +4.000 | Fibroblast | frozen_ruler | d0→d10 | +2.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | permuted_ruler_weights |
| GM23815 | +5.000 | PartialReprog | frozen_ruler | d0→d7 | +9.000 | +118.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | permuted_ruler_weights |
| GM23815 | +5.000 | PartialReprog | frozen_ruler | d0→d10 | +9.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | permuted_ruler_weights |
| GM23815 | +6.000 | EarlyPluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +6.000 | EarlyPluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +7.000 | PartialReprog | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +7.000 | PartialReprog | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +8.000 | NonReprog | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +8.000 | NonReprog | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +9.000 | Fibroblast | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +9.000 | Fibroblast | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +10.000 | Pluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +10.000 | Pluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +11.000 | PartialReprog | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +11.000 | PartialReprog | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +12.000 | PartialReprog | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +12.000 | PartialReprog | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +13.000 | EarlyPluripotency | frozen_ruler | d0→d7 | +56.000 | +164.000 | -2.160 | -7.076 | +4.917 | +0.050 | +9.000 | True | True | NA | permuted_ruler_weights |
| GM23815 | +13.000 | EarlyPluripotency | frozen_ruler | d0→d10 | +56.000 | +432.000 | -2.160 | -3.794 | +1.634 | +0.159 | +31.000 | False | True | NA | permuted_ruler_weights |
| GM23815 | +14.000 | Fibroblast | frozen_ruler | d0→d7 | +886.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | permuted_ruler_weights |
| GM23815 | +14.000 | Fibroblast | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +15.000 | EarlyPluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +15.000 | EarlyPluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +16.000 | Fibroblast | frozen_ruler | d0→d7 | +361.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | permuted_ruler_weights |
| GM23815 | +16.000 | Fibroblast | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +17.000 | EarlyPluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +17.000 | EarlyPluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +18.000 | EarlyPluripotency | frozen_ruler | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +18.000 | EarlyPluripotency | frozen_ruler | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | permuted_ruler_weights |
| GM23815 | +0.000 | Fibroblast | pluri_primary | d0→d7 | +6454.000 | +14.000 | -0.154 | +1.495 | -1.649 | +0.657 | +131.000 | False | True | NA | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +0.000 | Fibroblast | pluri_primary | d0→d10 | +6454.000 | +2.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +1.000 | EarlyPluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +1.000 | EarlyPluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +2.000 | EarlyPluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +2.000 | EarlyPluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +3.000 | EarlyPluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +3.000 | EarlyPluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +4.000 | Fibroblast | pluri_primary | d0→d7 | +2.000 | +2210.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +4.000 | Fibroblast | pluri_primary | d0→d10 | +2.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +5.000 | PartialReprog | pluri_primary | d0→d7 | +9.000 | +118.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +5.000 | PartialReprog | pluri_primary | d0→d10 | +9.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +6.000 | EarlyPluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +6.000 | EarlyPluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +7.000 | PartialReprog | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +7.000 | PartialReprog | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +8.000 | NonReprog | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +8.000 | NonReprog | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +9.000 | Fibroblast | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +9.000 | Fibroblast | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +10.000 | Pluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +10.000 | Pluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +11.000 | PartialReprog | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +11.000 | PartialReprog | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +12.000 | PartialReprog | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +12.000 | PartialReprog | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +13.000 | EarlyPluripotency | pluri_primary | d0→d7 | +56.000 | +164.000 | +1.778 | +8.758 | -6.980 | +0.990 | +198.000 | False | True | NA | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +13.000 | EarlyPluripotency | pluri_primary | d0→d10 | +56.000 | +432.000 | +1.778 | +6.874 | -5.096 | +0.970 | +194.000 | False | True | NA | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +14.000 | Fibroblast | pluri_primary | d0→d7 | +886.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +14.000 | Fibroblast | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +15.000 | EarlyPluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +15.000 | EarlyPluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +16.000 | Fibroblast | pluri_primary | d0→d7 | +361.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +16.000 | Fibroblast | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +17.000 | EarlyPluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +17.000 | EarlyPluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +18.000 | EarlyPluripotency | pluri_primary | d0→d7 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +18.000 | EarlyPluripotency | pluri_primary | d0→d10 | NA | NA | NA | NA | NA | NA | NA | NA | False | missing group row at an endpoint | size_matched_random_gene_sets_frozen_ruler_Z_not_AMS |
| GM23815 | +0.000 | Fibroblast | md_score | d0→d7 | +6454.000 | +14.000 | +0.483 | +0.401 | +0.082 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +0.000 | Fibroblast | md_score | d0→d10 | +6454.000 | +2.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +1.000 | EarlyPluripotency | md_score | d0→d7 | +0.000 | +3240.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +1.000 | EarlyPluripotency | md_score | d0→d10 | +0.000 | +8.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +2.000 | EarlyPluripotency | md_score | d0→d7 | +0.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +2.000 | EarlyPluripotency | md_score | d0→d10 | +0.000 | +2730.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +3.000 | EarlyPluripotency | md_score | d0→d7 | +0.000 | +2537.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +3.000 | EarlyPluripotency | md_score | d0→d10 | +0.000 | +30.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +4.000 | Fibroblast | md_score | d0→d7 | +2.000 | +2210.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +4.000 | Fibroblast | md_score | d0→d10 | +2.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +5.000 | PartialReprog | md_score | d0→d7 | +9.000 | +118.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +5.000 | PartialReprog | md_score | d0→d10 | +9.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +6.000 | EarlyPluripotency | md_score | d0→d7 | +0.000 | +1819.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +6.000 | EarlyPluripotency | md_score | d0→d10 | +0.000 | +330.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +7.000 | PartialReprog | md_score | d0→d7 | +0.000 | +20.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +7.000 | PartialReprog | md_score | d0→d10 | +0.000 | +10.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +8.000 | NonReprog | md_score | d0→d7 | +0.000 | +8.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +8.000 | NonReprog | md_score | d0→d10 | +0.000 | +1794.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +9.000 | Fibroblast | md_score | d0→d7 | +0.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +9.000 | Fibroblast | md_score | d0→d10 | +0.000 | +1362.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +10.000 | Pluripotency | md_score | d0→d7 | +0.000 | +977.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +10.000 | Pluripotency | md_score | d0→d10 | +0.000 | +181.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +11.000 | PartialReprog | md_score | d0→d7 | +0.000 | +65.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +11.000 | PartialReprog | md_score | d0→d10 | +0.000 | +10.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +12.000 | PartialReprog | md_score | d0→d7 | +0.000 | +106.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +12.000 | PartialReprog | md_score | d0→d10 | +0.000 | +25.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +13.000 | EarlyPluripotency | md_score | d0→d7 | +56.000 | +164.000 | +0.379 | +0.159 | +0.220 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +13.000 | EarlyPluripotency | md_score | d0→d10 | +56.000 | +432.000 | +0.379 | +0.235 | +0.145 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +14.000 | Fibroblast | md_score | d0→d7 | +886.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +14.000 | Fibroblast | md_score | d0→d10 | +886.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +15.000 | EarlyPluripotency | md_score | d0→d7 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +15.000 | EarlyPluripotency | md_score | d0→d10 | +0.000 | +761.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +16.000 | Fibroblast | md_score | d0→d7 | +361.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +16.000 | Fibroblast | md_score | d0→d10 | +361.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +17.000 | EarlyPluripotency | md_score | d0→d7 | +0.000 | +122.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +17.000 | EarlyPluripotency | md_score | d0→d10 | +0.000 | +121.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +18.000 | EarlyPluripotency | md_score | d0→d7 | +0.000 | +289.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +18.000 | EarlyPluripotency | md_score | d0→d10 | +0.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +0.000 | Fibroblast | tgfb_score | d0→d7 | +6454.000 | +14.000 | +0.033 | +0.004 | +0.029 | +0.154 | +30.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +0.000 | Fibroblast | tgfb_score | d0→d10 | +6454.000 | +2.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +1.000 | EarlyPluripotency | tgfb_score | d0→d7 | +0.000 | +3240.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +1.000 | EarlyPluripotency | tgfb_score | d0→d10 | +0.000 | +8.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +2.000 | EarlyPluripotency | tgfb_score | d0→d7 | +0.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +2.000 | EarlyPluripotency | tgfb_score | d0→d10 | +0.000 | +2730.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +3.000 | EarlyPluripotency | tgfb_score | d0→d7 | +0.000 | +2537.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +3.000 | EarlyPluripotency | tgfb_score | d0→d10 | +0.000 | +30.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +4.000 | Fibroblast | tgfb_score | d0→d7 | +2.000 | +2210.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +4.000 | Fibroblast | tgfb_score | d0→d10 | +2.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +5.000 | PartialReprog | tgfb_score | d0→d7 | +9.000 | +118.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +5.000 | PartialReprog | tgfb_score | d0→d10 | +9.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +6.000 | EarlyPluripotency | tgfb_score | d0→d7 | +0.000 | +1819.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +6.000 | EarlyPluripotency | tgfb_score | d0→d10 | +0.000 | +330.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +7.000 | PartialReprog | tgfb_score | d0→d7 | +0.000 | +20.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +7.000 | PartialReprog | tgfb_score | d0→d10 | +0.000 | +10.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +8.000 | NonReprog | tgfb_score | d0→d7 | +0.000 | +8.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +8.000 | NonReprog | tgfb_score | d0→d10 | +0.000 | +1794.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +9.000 | Fibroblast | tgfb_score | d0→d7 | +0.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +9.000 | Fibroblast | tgfb_score | d0→d10 | +0.000 | +1362.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +10.000 | Pluripotency | tgfb_score | d0→d7 | +0.000 | +977.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +10.000 | Pluripotency | tgfb_score | d0→d10 | +0.000 | +181.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +11.000 | PartialReprog | tgfb_score | d0→d7 | +0.000 | +65.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +11.000 | PartialReprog | tgfb_score | d0→d10 | +0.000 | +10.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +12.000 | PartialReprog | tgfb_score | d0→d7 | +0.000 | +106.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +12.000 | PartialReprog | tgfb_score | d0→d10 | +0.000 | +25.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +13.000 | EarlyPluripotency | tgfb_score | d0→d7 | +56.000 | +164.000 | +0.077 | +0.056 | +0.022 | +0.373 | +74.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +13.000 | EarlyPluripotency | tgfb_score | d0→d10 | +56.000 | +432.000 | +0.077 | +0.071 | +0.006 | +0.403 | +80.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +14.000 | Fibroblast | tgfb_score | d0→d7 | +886.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +14.000 | Fibroblast | tgfb_score | d0→d10 | +886.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +15.000 | EarlyPluripotency | tgfb_score | d0→d7 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +15.000 | EarlyPluripotency | tgfb_score | d0→d10 | +0.000 | +761.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +16.000 | Fibroblast | tgfb_score | d0→d7 | +361.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +16.000 | Fibroblast | tgfb_score | d0→d10 | +361.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +17.000 | EarlyPluripotency | tgfb_score | d0→d7 | +0.000 | +122.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +17.000 | EarlyPluripotency | tgfb_score | d0→d10 | +0.000 | +121.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +18.000 | EarlyPluripotency | tgfb_score | d0→d7 | +0.000 | +289.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +18.000 | EarlyPluripotency | tgfb_score | d0→d10 | +0.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +0.000 | Fibroblast | age_up | d0→d7 | +6454.000 | +14.000 | +0.042 | +0.023 | +0.020 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +0.000 | Fibroblast | age_up | d0→d10 | +6454.000 | +2.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +1.000 | EarlyPluripotency | age_up | d0→d7 | +0.000 | +3240.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +1.000 | EarlyPluripotency | age_up | d0→d10 | +0.000 | +8.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +2.000 | EarlyPluripotency | age_up | d0→d7 | +0.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +2.000 | EarlyPluripotency | age_up | d0→d10 | +0.000 | +2730.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +3.000 | EarlyPluripotency | age_up | d0→d7 | +0.000 | +2537.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +3.000 | EarlyPluripotency | age_up | d0→d10 | +0.000 | +30.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +4.000 | Fibroblast | age_up | d0→d7 | +2.000 | +2210.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +4.000 | Fibroblast | age_up | d0→d10 | +2.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +5.000 | PartialReprog | age_up | d0→d7 | +9.000 | +118.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +5.000 | PartialReprog | age_up | d0→d10 | +9.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +6.000 | EarlyPluripotency | age_up | d0→d7 | +0.000 | +1819.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +6.000 | EarlyPluripotency | age_up | d0→d10 | +0.000 | +330.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +7.000 | PartialReprog | age_up | d0→d7 | +0.000 | +20.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +7.000 | PartialReprog | age_up | d0→d10 | +0.000 | +10.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +8.000 | NonReprog | age_up | d0→d7 | +0.000 | +8.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +8.000 | NonReprog | age_up | d0→d10 | +0.000 | +1794.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +9.000 | Fibroblast | age_up | d0→d7 | +0.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +9.000 | Fibroblast | age_up | d0→d10 | +0.000 | +1362.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +10.000 | Pluripotency | age_up | d0→d7 | +0.000 | +977.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +10.000 | Pluripotency | age_up | d0→d10 | +0.000 | +181.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +11.000 | PartialReprog | age_up | d0→d7 | +0.000 | +65.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +11.000 | PartialReprog | age_up | d0→d10 | +0.000 | +10.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +12.000 | PartialReprog | age_up | d0→d7 | +0.000 | +106.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +12.000 | PartialReprog | age_up | d0→d10 | +0.000 | +25.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +13.000 | EarlyPluripotency | age_up | d0→d7 | +56.000 | +164.000 | +0.083 | +0.035 | +0.048 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +13.000 | EarlyPluripotency | age_up | d0→d10 | +56.000 | +432.000 | +0.083 | +0.046 | +0.037 | +0.005 | +0.000 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +14.000 | Fibroblast | age_up | d0→d7 | +886.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +14.000 | Fibroblast | age_up | d0→d10 | +886.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +15.000 | EarlyPluripotency | age_up | d0→d7 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +15.000 | EarlyPluripotency | age_up | d0→d10 | +0.000 | +761.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +16.000 | Fibroblast | age_up | d0→d7 | +361.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +16.000 | Fibroblast | age_up | d0→d10 | +361.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +17.000 | EarlyPluripotency | age_up | d0→d7 | +0.000 | +122.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +17.000 | EarlyPluripotency | age_up | d0→d10 | +0.000 | +121.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +18.000 | EarlyPluripotency | age_up | d0→d7 | +0.000 | +289.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +18.000 | EarlyPluripotency | age_up | d0→d10 | +0.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +0.000 | Fibroblast | age_down | d0→d7 | +6454.000 | +14.000 | +0.021 | +0.027 | -0.007 | +1.000 | +200.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +0.000 | Fibroblast | age_down | d0→d10 | +6454.000 | +2.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +1.000 | EarlyPluripotency | age_down | d0→d7 | +0.000 | +3240.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +1.000 | EarlyPluripotency | age_down | d0→d10 | +0.000 | +8.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +2.000 | EarlyPluripotency | age_down | d0→d7 | +0.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +2.000 | EarlyPluripotency | age_down | d0→d10 | +0.000 | +2730.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +3.000 | EarlyPluripotency | age_down | d0→d7 | +0.000 | +2537.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +3.000 | EarlyPluripotency | age_down | d0→d10 | +0.000 | +30.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +4.000 | Fibroblast | age_down | d0→d7 | +2.000 | +2210.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +4.000 | Fibroblast | age_down | d0→d10 | +2.000 | +9.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +5.000 | PartialReprog | age_down | d0→d7 | +9.000 | +118.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +5.000 | PartialReprog | age_down | d0→d10 | +9.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +6.000 | EarlyPluripotency | age_down | d0→d7 | +0.000 | +1819.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +6.000 | EarlyPluripotency | age_down | d0→d10 | +0.000 | +330.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +7.000 | PartialReprog | age_down | d0→d7 | +0.000 | +20.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +7.000 | PartialReprog | age_down | d0→d10 | +0.000 | +10.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +8.000 | NonReprog | age_down | d0→d7 | +0.000 | +8.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +8.000 | NonReprog | age_down | d0→d10 | +0.000 | +1794.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +9.000 | Fibroblast | age_down | d0→d7 | +0.000 | +3.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +9.000 | Fibroblast | age_down | d0→d10 | +0.000 | +1362.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +10.000 | Pluripotency | age_down | d0→d7 | +0.000 | +977.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +10.000 | Pluripotency | age_down | d0→d10 | +0.000 | +181.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +11.000 | PartialReprog | age_down | d0→d7 | +0.000 | +65.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +11.000 | PartialReprog | age_down | d0→d10 | +0.000 | +10.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +12.000 | PartialReprog | age_down | d0→d7 | +0.000 | +106.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +12.000 | PartialReprog | age_down | d0→d10 | +0.000 | +25.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +13.000 | EarlyPluripotency | age_down | d0→d7 | +56.000 | +164.000 | -0.038 | +0.001 | -0.039 | +1.000 | +200.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +13.000 | EarlyPluripotency | age_down | d0→d10 | +56.000 | +432.000 | -0.038 | -0.017 | -0.020 | +1.000 | +200.000 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +14.000 | Fibroblast | age_down | d0→d7 | +886.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +14.000 | Fibroblast | age_down | d0→d10 | +886.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +15.000 | EarlyPluripotency | age_down | d0→d7 | +0.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +15.000 | EarlyPluripotency | age_down | d0→d10 | +0.000 | +761.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +16.000 | Fibroblast | age_down | d0→d7 | +361.000 | +1.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +16.000 | Fibroblast | age_down | d0→d10 | +361.000 | +0.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +17.000 | EarlyPluripotency | age_down | d0→d7 | +0.000 | +122.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +17.000 | EarlyPluripotency | age_down | d0→d10 | +0.000 | +121.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +18.000 | EarlyPluripotency | age_down | d0→d7 | +0.000 | +289.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |
| GM23815 | +18.000 | EarlyPluripotency | age_down | d0→d10 | +0.000 | +5.000 | NA | NA | NA | NA | NA | NA | False | n_cells<10 at an endpoint | size_matched_random_gene_sets_mean_expression_bin |

MD cluster means vs `results/md2/t2_cluster_table.csv` (check, not a gate):

| cell_line | n_shared | max_abs_md_delta | n_gt_1e_6 |
|---|---|---|---|
| GM00731 | 33 | +0.000 | 0 |
| GM23815 | 34 | +0.000 | 0 |

## Task 1 reading

Fired key: `instruments_disagree_on_claim_population`.

Ruler does not decline in the pooled PartialReprog state while MD does, both with their own nulls. The instruments disagree on exactly the population the published claim concerns. No winner is declared. Pooled PartialReprog ruler: n_ok=2 n_pass=0 hits=[]. Pooled PartialReprog MD: n_ok=2 n_pass=1 hits=[{'cell_line': 'GM00731', 'group': 'PartialReprog', 'label': 'PartialReprog', 'endpoint': 'd0→d7', 'decline': 0.007173997855181413, 'p': 0.014925373134328358, 'n_cells_d0': 3.0, 'n_cells_end': 578.0}].

Aged GM00731 pooled PartialReprog ruler: n_ok=2 n_pass=0 hits=[].
Aged GM00731 pooled PartialReprog MD: n_ok=2 n_pass=1 hits=[{'cell_line': 'GM00731', 'group': 'PartialReprog', 'label': 'PartialReprog', 'endpoint': 'd0→d7', 'decline': 0.007173997855181413, 'p': 0.014925373134328358, 'n_cells_d0': 3.0, 'n_cells_end': 578.0}].

d7 is not the Stage 2 endpoint. Rows are not averaged across donors, clusters, or instruments. ≥10 is not primary.

## Task 2 pre-registration (verbatim, written before any three-way ρ or Δρ)

Written **before any MD3 three-way ρ vs donor age or paired Δρ**, 2026-09-19. No threshold in this block is re-tuned after numbers exist.

Three instruments, all applied unchanged to the same donors: (1) frozen GTEx ruler from results/fibro/frozen_ruler_ridge_raw.npz, used as-is, not refit; (2) MD score, AddModuleScore as in md2; (3) their age-up minus age-down score (mmc3 Table S3 Age up and Age down, AddModuleScore, same implementation as md2). Instrument (3) is AMS(Age up) − AMS(Age down).

Cohorts:
- CELLxGENE a19d1667-a7b5-4556-9e5f-f9bfa690c0f1 (10x, n=93) — already in results/fibro2/ and results/md2/. Frozen-ruler and MD numbers are read from those files and recomputed from the saved per-donor scores as a check. Age-up minus age-down is newly scored on the cached donor-count matrix.
- GSE226189 (bulk, n=82) — same: ruler and MD from disk; age-up minus age-down newly scored. geneCOUNT Tracking_ID is Ensembl; ID type is checked before mapping.
- GSE113957 (bulk, n=143, ages 1–92). FPKM blocks DESeq2 and TMM/log2-CPM but does not block rank-based scoring. Rank-transform genes within sample and score all three instruments on ranks. This is a different transform from the GTEx spec and is therefore a weaker test. Restrict to adults ≥18 to match the ruler's training range; report the excluded count. Transcript IDs are RefSeq; gene symbols are parsed from Annotation/Divergence. Unmappable IDs stop this cohort, recorded, not substituted.

Per instrument × cohort: Spearman ρ vs donor age, permutation null (shuffle ages, hold scores; n_perm=200, seed 20260914), p, bootstrap CI (B=200, seed 20260918), n genes mapped out of the list's total, and the ID type of the matrix. Uncalibrated R² is reported and is never a gate.

Then the statistic that matters: per cohort, the paired bootstrap difference in ρ between the frozen ruler and each of the other two instruments — resample donors, recompute both ρ on the same resampled donors, report Δρ = ρ_ruler − ρ_other with its 95% percentile CI. A Δρ CI that excludes zero is the claim; two separate ρ values with overlapping CIs are not. Report Δρ even where it favours the other instrument.

**Pre-registered reading:**
- Δρ CI excludes zero in favour of the ruler (CI lower bound > 0) in ≥2 of 3 named cohorts → MD and/or their aging signature predict donor age materially worse than the ruler on independent fibroblast cohorts.
- Δρ CI includes zero (in the cohorts that distinguish the instruments under test) → the instruments are not distinguishable on age prediction at these n, and the Task 1 disagreement cannot be attributed to one being a weaker age instrument. Say which.
- Mixed (one cohort excludes zero in favour of the ruler, another in favour of the other instrument, or fewer than two ruler-favouring cohorts) → report per cohort, do not average, do not declare a winner.

An unscored cohort does not count as excluding zero. Do not report two separate ρ values as evidence of a difference without the paired Δρ CI. Do not declare a winner between the instruments.


## Task 2 — instrument × cohort (ρ vs donor age)

| instrument | cohort | accession | n | n_donors | rho | rho_null | rho_p | rho_ci_lo | rho_ci_hi | r | cal_r2 | r2 | platform | id_type | transform | n_MD_mapped | n_age_up_mapped | n_age_down_mapped | n_overlap | reason |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| frozen_GTEx_fibroblast_ruler | CELLxGENE a19d1667-a7b5-4556-9e5f-f9bfa690c0f1 | a19d1667-a7b5-4556-9e5f-f9bfa690c0f1 | +93.000 | +93.000 | +0.455 | -0.003 | +0.005 | +0.253 | +0.605 | +0.414 | +0.171 | -9.664 | 10x (restricted) | symbol (var feature_name); gene_id Ensembl | TMM/log2-CPM frozen spec (read from fibro2; not recomputed) | NA | NA | NA | +17047.000 | NA |
| MD_AddModuleScore | CELLxGENE a19d1667-a7b5-4556-9e5f-f9bfa690c0f1 | a19d1667-a7b5-4556-9e5f-f9bfa690c0f1 | +93.000 | +93.000 | +0.253 | -0.004 | +0.005 | +0.025 | +0.408 | +0.251 | +0.063 | -9.655 | 10x (restricted) | symbol | AddModuleScore on LogNormalize of donor pseudobulk (md2) | +202.000 | NA | NA | NA | NA |
| age_up_minus_age_down_AddModuleScore | CELLxGENE a19d1667-a7b5-4556-9e5f-f9bfa690c0f1 | a19d1667-a7b5-4556-9e5f-f9bfa690c0f1 | +93.000 | +93.000 | +0.161 | +0.008 | +0.085 | -0.003 | +0.325 | +0.052 | +0.003 | -9.821 | 10x (restricted) | symbol | AddModuleScore on LogNormalize of donor pseudobulk | NA | +986.000 | +1595.000 | NA | NA |
| frozen_GTEx_fibroblast_ruler | GSE226189 bulk | GSE226189 | +82.000 | +82.000 | +0.549 | +0.005 | +0.005 | +0.374 | +0.667 | +0.547 | +0.299 | -9.963 | bulk_rnaseq | ensembl Tracking_ID | TMM/log2-CPM frozen spec (read from fibro2) | NA | NA | NA | +21351.000 | NA |
| MD_AddModuleScore | GSE226189 bulk | GSE226189 | +82.000 | +82.000 | +0.082 | +0.011 | +0.259 | -0.132 | +0.291 | +0.051 | +0.003 | -7.603 | bulk_rnaseq | ensembl mapped to symbol | AddModuleScore on LogNormalize (md2) | +204.000 | NA | NA | NA | NA |
| age_up_minus_age_down_AddModuleScore | GSE226189 bulk | GSE226189 | +82.000 | +82.000 | +0.037 | +0.004 | +0.363 | -0.182 | +0.226 | +0.018 | +0.000 | -7.749 | bulk_rnaseq | symbol | AddModuleScore on LogNormalize; Ensembl mapped to symbols | NA | +1040.000 | +1745.000 | NA | NA |
| STOP | GSE113957 bulk FPKM ranks | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | cut_number did not yield 24 bins (got 22). Not dropping bins. |

## Task 2 — paired bootstrap Δρ (ρ_ruler − ρ_other; this is the claim)

| cohort | other | n_donors | rho_ruler | rho_other | delta_rho | delta_ci_lo | delta_ci_hi | excludes_zero | favours | note |
|---|---|---|---|---|---|---|---|---|---|---|
| CELLxGENE a19d1667-a7b5-4556-9e5f-f9bfa690c0f1 | MD_AddModuleScore | 93 | +0.455 | +0.253 | +0.202 | +0.006 | +0.440 | True | ruler | Δρ = ρ_ruler − ρ_other; paired donor bootstrap. Overlapping separate CIs are not the claim. |
| CELLxGENE a19d1667-a7b5-4556-9e5f-f9bfa690c0f1 | age_up_minus_age_down | 93 | +0.455 | +0.161 | +0.295 | +0.107 | +0.487 | True | ruler | Δρ = ρ_ruler − ρ_other; paired donor bootstrap. Overlapping separate CIs are not the claim. |
| GSE226189 bulk | MD_AddModuleScore | 82 | +0.549 | +0.082 | +0.467 | +0.239 | +0.696 | True | ruler | Δρ = ρ_ruler − ρ_other; paired donor bootstrap. Overlapping separate CIs are not the claim. |
| GSE226189 bulk | age_up_minus_age_down | 82 | +0.549 | +0.037 | +0.511 | +0.248 | +0.772 | True | ruler | Δρ = ρ_ruler − ρ_other; paired donor bootstrap. Overlapping separate CIs are not the claim. |

## Task 2 reading

Fired key: `ruler_stronger_age_in_ge2`.

Δρ CI excludes zero in favour of the ruler in 2 of 3 named cohorts. MD and/or their aging signature predict donor age materially worse than the ruler on independent fibroblast cohorts.

## Task 3 pre-registration (verbatim, written before any MD co-variation ρ)

Written **before any MD3 MD co-variation ρ**, 2026-09-19. Report-only. No gate.

Across the GSE297234 Louvain cluster × timepoint rows with ≥20 cells (the FINDINGS_MD2.md unit, so the MD-vs-ruler cell is a check), report Spearman ρ between MD score and each of: the frozen ruler (reproduce the FINDINGS_MD2.md values as a check), the pluripotency−fibroblast score, the TGF-β score, and cell-state label treated as an ordinal along their trajectory.

Ordinal, frozen before this ρ: Fibroblast=0, PartialReprog=1, EarlyPluripotency=2, Pluripotency=3. NonReprog is not on that axis and is excluded from the ordinal ρ; n excluded is reported. This is our coding of their OSKM progression, not theirs.

Permutation null n_perm=200 seed 20260914 (shuffle MD, hold the other variable). Bootstrap CI B=200 seed 20260918. One table, both donors separate, never averaged. This says what MD co-varies with in this dataset without claiming what it “really” measures.


Ordinal used: {'Fibroblast': 0, 'PartialReprog': 1, 'EarlyPluripotency': 2, 'Pluripotency': 3}. NonReprog excluded.

## Task 3 — MD co-variation (report-only, no gate)

| cell_line | other | n | rho | rho_null | rho_p | rho_ci_lo | rho_ci_hi | n_excluded_NonReprog | note |
|---|---|---|---|---|---|---|---|---|---|
| GM00731 | frozen_ruler | 33 | +0.782 | -0.014 | +0.005 | +0.603 | +0.879 | NA | permutation shuffles MD, holds the other variable; unit=cluster×timepoint ≥20 |
| GM00731 | pluri_primary | 33 | -0.938 | +0.018 | +1.000 | -0.969 | -0.850 | NA | permutation shuffles MD, holds the other variable; unit=cluster×timepoint ≥20 |
| GM00731 | tgfb_score | 33 | +0.651 | +0.017 | +0.005 | +0.352 | +0.835 | NA | permutation shuffles MD, holds the other variable; unit=cluster×timepoint ≥20 |
| GM00731 | state_ordinal_Fibroblast0_Partial1_EarlyPluri2_Pluri3 | 24 | -0.601 | +0.012 | +1.000 | -0.817 | -0.156 | +9.000 | permutation shuffles MD, holds the other variable; unit=cluster×timepoint ≥20; NonReprog excluded n=9 |
| GM23815 | frozen_ruler | 34 | +0.766 | +0.013 | +0.005 | +0.577 | +0.878 | NA | permutation shuffles MD, holds the other variable; unit=cluster×timepoint ≥20 |
| GM23815 | pluri_primary | 34 | -0.960 | -0.015 | +1.000 | -0.979 | -0.883 | NA | permutation shuffles MD, holds the other variable; unit=cluster×timepoint ≥20 |
| GM23815 | tgfb_score | 34 | +0.650 | +0.002 | +0.005 | +0.434 | +0.822 | NA | permutation shuffles MD, holds the other variable; unit=cluster×timepoint ≥20 |
| GM23815 | state_ordinal_Fibroblast0_Partial1_EarlyPluri2_Pluri3 | 33 | -0.581 | +0.016 | +1.000 | -0.768 | -0.216 | +1.000 | permutation shuffles MD, holds the other variable; unit=cluster×timepoint ≥20; NonReprog excluded n=1 |

MD vs frozen ruler vs FINDINGS_MD2.md (check):

| cell_line | md3_rho | md2_rho | abs_delta | n_md3 | n_md2 |
|---|---|---|---|---|---|
| GM00731 | +0.782 | +0.782 | +0.000 | 33 | 33 |
| GM23815 | +0.766 | +0.766 | +0.000 | 34 | 34 |

## What this changes in FINDINGS_MD2.md / FINDINGS_FIBRO3.md / FINDINGS_FIBRO.md (quote, not edited)

This file does not edit `FINDINGS_MD2.md`, `FINDINGS_FIBRO.md`, `FINDINGS_FIBRO2.md`, `FINDINGS_FIBRO3.md`, or `FALSIFICATION.md`.

FINDINGS_MD2.md Task 1 reading, quoted:

> No cluster at any resolution shows p ≤ 0.05. The per-cluster negative holds and is not a resolution artifact. The frozen ruler does not move in these cells at any granularity tested.

FINDINGS_MD2.md Task 2 reading, quoted:

> MD moves per cluster, the frozen ruler does not. The instruments disagree on the same cells. What would settle it: an independent fibroblast age instrument scored on these same Louvain clusters, or the authors' Seurat/sctransform object with the frozen ruler projected onto it. No winner is declared.

FINDINGS_MD2.md pairable-cluster sentence, quoted:

> Aged-donor Louvain clusters pairable at both endpoints (ok=True): n=1 (10:NonReprog). PartialReprog among them: n=0.

FINDINGS_FIBRO.md Stage 2 verdict sentence, quoted:

> age score does not decline → the ruler reads a static donor property, not a modifiable state. Step 3 is not supported by this data.

That d0→d10 all-cell verdict stands. Nothing here promotes d7 to the primary Stage 2 result.

FINDINGS_FIBRO3.md per-cluster d0→d7 sentence, quoted:

> Aged GM00731 d0→d7 (the all-cell T-B cell had p=+0.005): cluster 0 decline=-0.582 p=+0.652 n_ge=130/200 n_d0=4941 n_d7=2719; cluster 1 decline=+2.215 p=+0.144 n_ge=28/200 n_d0=76 n_d7=850; cluster 2 skipped (n_cells<20 at an endpoint).

Task 1 of this file tests the pooled PartialReprog state, which those sentences did not. The quoted sentences are not edited there.

## Limitations

1. Two donors in GSE297234.
2. Clusters are not lineage-tracked. PartialReprog / NonReprog labels are argmax of mmc3 Reprog_cell_state_signatures AddModuleScore, not Slingshot trajectories from day-0 fibroblasts. Pooled PartialReprog at d0 is the cells that sit in clusters whose all-timepoint argmax is PartialReprog, not a lineage from a d0 fibroblast cluster.
3. Our reproduction of their metrics is ours, not theirs. AddModuleScore is a Python reimplementation of the published algorithm, not Seurat's C++/R object.
4. Louvain is a Python deviation from their R/sctransform pipeline, as recorded in FINDINGS_MD2.md (LogNormalize + quadratic HVG + percent.mt residualization + PCA 50 + SNN + networkx louvain, resolution=0.8 Seurat default; STAR Methods omit resolution).
5. GSE113957 uses a rank transform. That is a different transform from the GTEx TMM/log2-CPM spec and is a weaker test. Frozen mu/sd (log2-CPM units) are applied as-is to ranks; the ruler is not refit.
6. Frozen ruler trained on GTEx V10 cultured fibroblasts, public `AGE` 10-year bins.
7. Cross-platform shift from GTEx bulk polyA (RNASeQCv2.4.2) to 10x 3' scRNA-seq pseudobulk and to bulk FPKM ranks.
8. Missing overlap genes at z=0. Nothing fitted on GSE297234 or on the external matrices.
9. MIN_CELLS=20 is primary for per-cluster. ≥10 is a labelled sensitivity check. Pooled state has no extra cutoff; n is reported.
10. Seed `20260914`. Bootstrap seed `20260918`. n_perm=200, n_boot=200, n_random=200.
11. GSE325735 not opened. d7 is not the Stage 2 endpoint.
12. Table S3 Age up/down are the published lists, not a rebuilt DESeq2 analysis.
13. Pluripotency−fibroblast null is size-matched in frozen-ruler Z space, not AMS expression-bin matching.

## Open questions

1. Would the authors' sctransform v2 object put more PartialReprog cells at d0, so the pooled d0→d7 ruler test is not a 3-cell versus 578-cell comparison?
2. If GSE113957 integer counts and the PCA old/young split were released, how would a rebuilt DE list compare to Table S3 on these cohorts?
3. Slingshot from day-0 fibroblast clusters was not run.

## ID types actually read

- **mmc3_Age_up** kind=`symbol` n=1533 n_ensembl=437 n_refseq=0 n_symbol_like=1095 source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx:Aging_signatures:Age up` examples=['SCN1A', 'STEAP4', 'ADH1B', 'KCNA4', 'LEP', 'SLC25A52', 'SRD5A2', 'FMO2', 'MAOA', 'ENPP5', 'GOLGA8M', 'KCNJ2']
- **mmc3_Age_down** kind=`symbol` n=2007 n_ensembl=213 n_refseq=0 n_symbol_like=1791 source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx:Aging_signatures:Age down` examples=['CACNG4', 'CACNA1H', 'STRA6', 'SLC35F3', 'RRM2', 'DGKG', 'ALDH1A1', 'TMEM191A', 'NDUFA4L2', 'HS3ST3A1', 'HAS1', 'DHRS2']
- **louvain_GM00731_symbols** kind=`symbol` n=29058 n_ensembl=0 n_refseq=0 n_symbol_like=29057 source=`<repo>\data\processed\md2\louvain_counts_GM00731.npz:symbols` examples=['AL627309.1', 'AL627309.5', 'AP006222.2', 'AL669831.2', 'LINC01409', 'FAM87B', 'LINC01128', 'LINC00115', 'FAM41C', 'AL645608.6', 'AL645608.2', 'AL645608.4']
- **louvain_GM00731_gene_id** kind=`ensembl` n=29058 n_ensembl=29058 n_refseq=0 n_symbol_like=0 source=`<repo>\data\processed\md2\louvain_counts_GM00731.npz:gene_id` examples=['ENSG00000238009', 'ENSG00000241860', 'ENSG00000286448', 'ENSG00000229905', 'ENSG00000237491', 'ENSG00000177757', 'ENSG00000228794', 'ENSG00000225880', 'ENSG00000230368', 'ENSG00000272438', 'ENSG00000230699', 'ENSG00000241180']
- **louvain_GM23815_symbols** kind=`symbol` n=29154 n_ensembl=0 n_refseq=0 n_symbol_like=29154 source=`<repo>\data\processed\md2\louvain_counts_GM23815.npz:symbols` examples=['AL627309.1', 'AL627309.3', 'AL627309.5', 'AP006222.2', 'LINC01409', 'FAM87B', 'LINC01128', 'LINC00115', 'FAM41C', 'AL645608.6', 'AL645608.2', 'AL645608.4']
- **louvain_GM23815_gene_id** kind=`ensembl` n=29154 n_ensembl=29154 n_refseq=0 n_symbol_like=0 source=`<repo>\data\processed\md2\louvain_counts_GM23815.npz:gene_id` examples=['ENSG00000238009', 'ENSG00000239945', 'ENSG00000241860', 'ENSG00000286448', 'ENSG00000237491', 'ENSG00000177757', 'ENSG00000228794', 'ENSG00000225880', 'ENSG00000230368', 'ENSG00000272438', 'ENSG00000230699', 'ENSG00000241180']
- **cxg_gene_id** kind=`ensembl` n=25952 n_ensembl=25952 n_refseq=0 n_symbol_like=0 source=`<repo>\data\processed\md2\cxg_donor_counts.npz:genes` examples=['ENSG00000278457', 'ENSG00000121410', 'ENSG00000268895', 'ENSG00000148584', 'ENSG00000175899', 'ENSG00000245105', 'ENSG00000166535', 'ENSG00000256661', 'ENSG00000256904', 'ENSG00000256069', 'ENSG00000128274', 'ENSG00000118017']
- **cxg_symbols** kind=`symbol` n=25952 n_ensembl=8 n_refseq=0 n_symbol_like=25940 source=`<repo>\data\processed\md2\cxg_donor_counts.npz:symbols` examples=['5S_rRNA', 'A1BG', 'A1BG-AS1', 'A1CF', 'A2M', 'A2M-AS1', 'A2ML1', 'A2ML1-AS1', 'A2ML1-AS2', 'A2MP1', 'A4GALT', 'A4GNT']
- **cxg_a19d1667_S3_symbols** kind=`symbol` n=25952 n_ensembl=8 n_refseq=0 n_symbol_like=25940 source=`cxg_a19d1667_S3` examples=['5S_rRNA', 'A1BG', 'A1BG-AS1', 'A1CF', 'A2M', 'A2M-AS1', 'A2ML1', 'A2ML1-AS1', 'A2ML1-AS2', 'A2MP1', 'A4GALT', 'A4GNT']
- **gse226189_Tracking_ID** kind=`ensembl` n=57773 n_ensembl=57773 n_refseq=0 n_symbol_like=0 source=`<repo>\data\raw\fibro2\GSE226189\GSE226189_RAW.tar:Tracking_ID` examples=['ENSG00000000003', 'ENSG00000000005', 'ENSG00000000419', 'ENSG00000000457', 'ENSG00000000460', 'ENSG00000000938', 'ENSG00000000971', 'ENSG00000001036', 'ENSG00000001084', 'ENSG00000001167', 'ENSG00000001460', 'ENSG00000001461']
- **gse226189_mapped_symbols** kind=`symbol` n=57773 n_ensembl=23884 n_refseq=0 n_symbol_like=33883 source=`Ensembl→symbol via 10x+ruler` examples=['TSPAN6', 'TNMD', 'DPM1', 'SCYL3', 'C1ORF112', 'FGR', 'CFH', 'FUCA2', 'GCLC', 'NFYA', 'STPG1', 'NIPAL3']
- **GSE226189_S3_symbols** kind=`symbol` n=57773 n_ensembl=23884 n_refseq=0 n_symbol_like=33883 source=`GSE226189_S3` examples=['TSPAN6', 'TNMD', 'DPM1', 'SCYL3', 'C1ORF112', 'FGR', 'CFH', 'FUCA2', 'GCLC', 'NFYA', 'STPG1', 'NIPAL3']
- **gse113957_Transcript_ID** kind=`refseq_transcript` n=27142 n_ensembl=0 n_refseq=27142 n_symbol_like=0 source=`<repo>\data\raw\fibro2\GSE113957\GSE113957_fpkm.txt.gz:Transcript ID` examples=['NM_173803', 'NM_014423', 'NM_001103167', 'NR_134623', 'NR_024490', 'NM_018397', 'NM_001037671', 'NR_106859', 'NM_133458', 'NM_001080424', 'NM_001271816', 'NR_027283']
- **gse113957_parsed_symbols** kind=`symbol` n=27142 n_ensembl=0 n_refseq=0 n_symbol_like=27130 source=`<repo>\data\raw\fibro2\GSE113957\GSE113957_fpkm.txt.gz:Annotation/Divergence` examples=['MPV17L', 'AFF4', 'ZGLP1', 'LOC105373021', 'GABPB1-AS1', 'CHDH', 'C12ORF74', 'MIR6801', 'ZFP90', 'KDM6B', 'DCLRE1A', 'LOC440461']
- **gse113957_gene_symbols_collapsed** kind=`symbol` n=27142 n_ensembl=0 n_refseq=0 n_symbol_like=27130 source=`collapsed Annotation/Divergence` examples=['MPV17L', 'AFF4', 'ZGLP1', 'LOC105373021', 'GABPB1-AS1', 'CHDH', 'C12ORF74', 'MIR6801', 'ZFP90', 'KDM6B', 'DCLRE1A', 'LOC440461']

## Columns actually read

- **mmc3_MD_signatures** source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx` columns=['MD score', 'TGFB score']
- **mmc3_Aging_signatures** source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx` columns=['Age up', 'Age down']
- **mmc3_Fibroblast_subtype_signatures** source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx` columns=['PI16_univ', 'LRRC15_myo', 'COL3A1_myo']
- **mmc3_Reprog_cell_state_signatures** source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx` columns=['Fibroblast', 'PartialReprog', 'EarlyPluripotency', 'Pluripotency', 'NonReprog']
- **t2_cluster_labels_GM00731** source=`<repo>\results\md2\t2_cluster_labels_GM00731.csv` columns=['cell_line', 'cluster', 'n_cells', 'label', 'mean_Fibroblast', 'mean_PartialReprog', 'mean_EarlyPluripotency', 'mean_Pluripotency', 'mean_NonReprog']
- **louvain_obs_GM00731** source=`<repo>\results\md2\louvain_obs_GM00731.csv` columns=['barcode', 'gsm', 'file', 'cell_line', 'day', 'age_years', 'age_source', 'umi', 'n_genes', 'mito_frac', 'cluster', 'louvain_kept', 'percent_mt']
- **t2_cluster_labels_GM23815** source=`<repo>\results\md2\t2_cluster_labels_GM23815.csv` columns=['cell_line', 'cluster', 'n_cells', 'label', 'mean_Fibroblast', 'mean_PartialReprog', 'mean_EarlyPluripotency', 'mean_Pluripotency', 'mean_NonReprog']
- **louvain_obs_GM23815** source=`<repo>\results\md2\louvain_obs_GM23815.csv` columns=['barcode', 'gsm', 'file', 'cell_line', 'day', 'age_years', 'age_source', 'umi', 'n_genes', 'mito_frac', 'cluster', 'louvain_kept', 'percent_mt']
- **md2_t2_cluster_table** source=`<repo>\results\md2\t2_cluster_table.csv` columns=['resolution', 'cell_line', 'day', 'cluster', 'n_cells', 'n_timepoint', 'frac_timepoint', 'median_UMI', 'median_genes', 'mito_frac_median_cell', 'below_min_cells', 'md_score', 'tgfb_score', 'state_Fibroblast', 'state_PartialReprog', 'state_EarlyPluripotency', 'state_Pluripotency', 'state_NonReprog', 'age_score', 'pluri_with_OSKM', 'pluri_without_OSKM', 'pluri_primary', 'age_source', 'label']
- **t2_cxg_ruler_scores** source=`<repo>\results\fibro2\ta_cxg_a19d1667_scores.csv` columns=['sample', 'donor', 'age_years', 'age_source', 'n_cells', 'assay', 'age_score']
- **t2_cxg_md_scores** source=`<repo>\results\md2\t3_cxg_a19d1667_MD_scores.csv` columns=['sample', 'donor', 'age_years', 'age_source', 'n_cells', 'assay', 'age_score', 'md_score']
- **t2_gse226189_ruler_scores** source=`<repo>\results\fibro2\ta_GSE226189_near_miss_bulk_scores.csv` columns=['sample', 'donor', 'age_years', 'age_source', 'title', 'gsm', 'age_score']
- **t2_gse226189_md_scores** source=`<repo>\results\md2\t3_GSE226189_MD_scores.csv` columns=['sample', 'donor', 'age_years', 'age_source', 'title', 'gsm', 'age_score', 'md_score']
- **t2_gse226189_10x_map** source=`<repo>\data\processed\md2\allcell_counts_GM00731.npz` columns=['gene_id', 'symbol']
- **t2_gse113957_samples** source=`<repo>\results\fibro2\ta_GSE113957_samples.csv` columns=['gsm', 'title', 'source_name', 'organism', 'molecule', 'description', 'characteristics', 'platform', 'supplementary', 'donor', 'age_years', 'age_source', 'char_keys']
- **t2_gse113957_fpkm** source=`<repo>\data\raw\fibro2\GSE113957\GSE113957_fpkm.txt.gz` columns=['Transcript ID', 'chr', 'start', 'end', 'strand', 'Length', 'Copies', 'Annotation/Divergence', '101_19yr_Female_Caucasian', '102_19yr_Male_Caucasian', '103_19yr_Male_Caucasian', '104_19yr_Male_Caucasian', '105_20yr_Male_Caucasian', '106_20yr_Female_Caucasian', '107_31yr_Female_Caucasian', '108_31_female_Caucasian', '109_32_male_Unknown', '110_32_female_BlackPuertoRican', '111_33yr_Male_Caucasian', '112_33yr_Male_Caucasian']
- **t1_cluster20_scores** source=`<repo>\results\md3\t1_cluster20_scores.csv` columns=['day', 'n_cells', 'n_timepoint', 'frac_timepoint', 'below_min_cells', 'cluster', 'label', 'md_score', 'tgfb_score', 'age_up', 'age_down', 'age_score', 'pluri_with_OSKM', 'pluri_without_OSKM', 'pluri_primary', 'age_source', 'cell_line', 'resolution', 'min_cells', 'aggregation']
- **md2_t2_crosscorr** source=`<repo>\results\md2\t2_crosscorr.csv` columns=['cell_line', 'n', 'min_cells', 'n_perm', 'n_boot', 'seed', 'boot_seed', 'rho', 'rho_null', 'rho_p', 'rho_ci_lo', 'rho_ci_hi', 'note']

## Gene-set names actually read

- **mmc3_Aging_signatures_Age up** name=`Aging_signatures:Age up` n=1533 id_type=symbol source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx`
- **mmc3_Aging_signatures_Age down** name=`Aging_signatures:Age down` n=2007 id_type=symbol source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx`
- **MD_built_from_md2** name=`MD (STAR Methods built; md2)` n=205 id_type=symbol source=`<repo>\results\md2\genesets.json`
- **TGFB_built_from_md2** name=`TGFB (STAR Methods built; md2)` n=54 id_type=symbol source=`<repo>\results\md2\genesets.json`

## Files

- `src/md3_common.py`, `md3_idtype.py`, `md3_genesets.py`, `md3_task1.py`, `md3_task2.py`, `md3_task3.py`, `md3_findings.py`, `md3_run.py`
- `results/md3/`
- `FINDINGS_MD3.md`
- `PROGRESS_MD3.md`

## Supersession note, appended 2026-10-02: atlas ruler − MD Δρ on 93 pseudobulks

Appended; nothing above is changed. The Task 2 reading above (`ruler_stronger_age_in_ge2`, "Δρ CI excludes zero in favour of the ruler in 2 of 3 named cohorts") counted the fibroblast atlas as one of the two cohorts whose ruler − MD interval excluded zero. That atlas interval no longer stands.

The ruler − MD Δρ computed with each of the 93 CELLxGENE atlas pseudobulks (a19d1667) as a unit (+0.202, 95% CI +0.006 to +0.440, `excludes_zero=True`; `results/md3/t2_delta_rho.csv`, row `cxg_ruler_minus_MD`) is **superseded** by the per-person re-analysis, in which pseudobulk scores are averaged within each of 65 donors (+0.171, 95% CI −0.118 to +0.412, `excludes_zero=False`; `results/paper_figs/cxg_true_donor_results.csv`, `analysis=true_donor`, `contrast=ruler_minus_MD_AddModuleScore`). The manuscript reports the per-person value in Table 2 (`paper/rejuvenation_readouts_preprint_11.docx`). **The per-person re-analysis was not pre-registered:** `results/paper_figs/` has no PREREG flag, and its results were written on 2026-09-23, after the pre-registered 93-pseudobulk result was known. After regrouping, the ruler − MD interval excludes zero only in GSE226189 (+0.467, 95% CI +0.239 to +0.696). The atlas ruler − (age-up − age-down) interval excludes zero in both units (+0.295, CI +0.107 to +0.487 on 93 pseudobulks; +0.372, CI +0.094 to +0.595 on 65 donors). The pre-registration flags are not edited. The reading is listed as B7 in `results/verify/WITHDRAWN_READINGS.md`.

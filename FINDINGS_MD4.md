# FINDINGS_MD4 — cell-state contrast at a shared timepoint (Figure 3G)

**Status:** Task 1: GM00731 d3: `both_lower`; d7: `both_lower`; GM23815 d10: `both_lower`. Task 2: report-only extrapolation tabulated; no gate. Task 3: FINDINGS_MD3.md Task 1 `instruments_disagree_on_claim_population` withdrawn as unsupported (n₀=3); superseded here. MD3 Task 2 Δρ stands. Fired: GM00731 `both_lower`; GM23815 `both_lower`. Seed `20260914`. boot `20260918`. n_perm=200. n_boot=200. n_random=200. Frozen ruler `frozen_ruler_ridge_raw.npz` exists=True.

Does not modify any existing `FINDINGS_*.md` or `FALSIFICATION.md`. No GSE325735. Uncalibrated R² is reported and is never a gate. Nothing averaged across donors, states, timepoints, or instruments. Refit nothing. d7 is not the primary Stage 2 endpoint. Task 2 extrapolation is report-only and is not used to discount Task 1.

Reproduced by `src/md4_run.py`. Frozen ruler: `<repo>\results\fibro\frozen_ruler_ridge_raw.npz` exists=True.

Flag: `PREREG_TASK1.flag` exists=True. `PREREG_TASK2.flag` exists=True. `PREREG_TASK3.flag` exists=True. `DECLARED_BEFORE_SCORES.flag` exists=True.

## Why the MD3 test was the wrong one

FINDINGS_MD3.md Task 1 tested a within-state trajectory (d0→d7, d0→d10) on the pooled PartialReprog population. At day 0 that pool has n₀=3 cells on aged GM00731 (FINDINGS_MD3.md occupancy: PartialReprog n_d0=3, n_d3=5220, n_d7=578, n_d10=31). The state essentially does not exist before induction. Every number resting on that baseline (ruler p=+0.403, MD p=+0.015) is uninterpretable, and the recorded status `instruments_disagree_on_claim_population` is not supported by it.

Lu et al. Figure 3G compares cell states against each other at timepoints where both exist: partially reprogrammed vs non-reprogrammed. FINDINGS_MD2.md already had the counts for that contrast on the aged donor (PartialReprog n=5,832, NonReprog n=5,459, pooled across timepoints). This file runs that contrast per donor, per timepoint, with no pooling.

Paper quote (Lu et al., Cell 2025, 188:5895–5911):

> The partially reprogrammed populations showed a robust reversal of transcriptomic aging changes, while non-reprogrammed populations maintained the aging transcriptome (Figure 3G). This pattern was accompanied by downregulation of the MD and TGF-β pathway scores in the partially reprogramming cells and their maintained expression in non-reprogrammed cells (Figure 3G).

## Failures recorded

None recorded.

## Task 1 pre-registration (verbatim, written before any MD4 cell-level instrument score)

Written **before any MD4 cell-level instrument score**, 2026-09-19. No threshold in this block is re-tuned after numbers exist.

Using the src/md2_*.py Louvain clustering already computed (not re-run, not re-tuned). Cell-state labels are the md2 assignment: argmax of mean AddModuleScore of mmc3 Reprog_cell_state_signatures columns Fibroblast / PartialReprog / EarlyPluripotency / Pluripotency / NonReprog. Labels are read from results/md2/t2_cluster_labels_*.csv. Not re-labelled. AddModuleScore for MD and TGF-β is the per-cell vector already stored in data/processed/md2/louvain_ams_*.npz (the md2 computation). Table S3 Age up / Age down are scored with the same AddModuleScore implementation (nbin=24, ctrl=100) on the same LogNormalize Louvain matrix, using the bins stored in that npz. Frozen ruler: results/fibro/frozen_ruler_ridge_raw.npz used as-is. Refit nothing.

Per donor, per timepoint, never pooled across donors or timepoints. A timepoint qualifies for a two-state contrast iff both states have ≥50 Louvain-kept cells at that donor × day. Cell counts and which timepoints qualify are reported before contrasts.

Per cell:
- frozen ruler score = edgeR log2-CPM (prior.count=2) on that cell's ruler-aligned UMI counts with TMM size factor fixed at 1 (library-size only). TMM is not re-estimated among cells, among states, or among timepoints: a TMM among ~20k cells would be a new size-factor fit on GSE297234 and is not the frozen bulk spec. Frozen μ, σ, w applied as-is; overlap genes that are all-zero in that donor's Louvain ruler matrix are left at z=0.
- MD score, from md2 louvain_ams (not recomputed).
- TGF-β score, from md2 louvain_ams (not recomputed).
- Table S3 age-up minus age-down = AMS(Age up) − AMS(Age down).
- pluripotency−fibroblast = mean-z(PLURI_ENDOGENOUS) − mean-z(FIBRO_IDENTITY) on the same per-cell Z as the ruler, drop_oskm=False (md2 pluri_primary).

Primary contrast: PartialReprog vs NonReprog. Context rows (same procedure, not pooled into the primary): Fibroblast vs NonReprog; PartialReprog vs Fibroblast.

Per qualifying timepoint, per instrument: Δmean = mean(state_A) − mean(state_B); Δmedian = median(state_A) − median(state_B). Cell-level permutation null: 200 shuffles of the two-state labels among the cells in those two states at that donor × timepoint (group sizes held; timepoint and donor held, so neither can drive the null). Empirical p_lower = (n_null ≤ Δobs + 1) / (n_null + 1); p_two_sided = (n_|null| ≥ |Δobs| + 1) / (n_null + 1). Donor-cell bootstrap 95% percentile CI on Δmean and Δmedian: resample cells within each state independently, B=200, seed 20260918. Effect size: Cohen's d on the mean difference, pooled sd with n−1, same sign as Δmean. Not rank-biserial. Seed 20260914 for permutations with a documented per-(donor, day, contrast) offset.

Sanity check, before any new score: reproduce FINDINGS_MD2.md pooled MD means on the aged donor (PartialReprog +0.150, NonReprog +0.439; t2_reproduction.json). If they do not reproduce, stop. Cell-state labels have drifted.

**Pre-registered reading (only the outcome that fired; each donor separately, never pooled; each qualifying timepoint reported, never pooled):**
- MD lower in PartialReprog than NonReprog and the frozen ruler lower too, both beating their nulls → both instruments agree that partially reprogrammed cells read younger. Their claim is supported and our earlier negatives reflected the wrong test (a within-state trajectory with no baseline). Say so plainly and quote the sentences superseded.
- MD lower, ruler not lower (or higher) → the instruments disagree on the population the published claim concerns, in a well-powered comparison. Report both numbers side by side, declare no winner, and name what would settle it.
- Neither differs → our reproduction does not recover their Figure 3G contrast; report the reproduction failure as the finding, not as a refutation.
- Ruler lower, MD not → report as an unexplained discrepancy.

“Lower” / “beating their nulls” = Δmean = mean(PartialReprog) − mean(NonReprog) is < 0 and p_lower ≤ 0.05 on that instrument's cell-level permutation null, on that donor × timepoint. GM23815 is tabulated, never pooled with GM00731, never averaged.

If two qualifying timepoints in one donor fire different bullets, report both, do not pick, do not average.


## Task 1 — cell counts by state × timepoint (md2 Louvain labels, not re-labelled)

Assignment basis: argmax of mean AddModuleScore of mmc3 `Reprog_cell_state_signatures` (Fibroblast, PartialReprog, EarlyPluripotency, Pluripotency, NonReprog). Labels from `results/md2/t2_cluster_labels_*.csv`. Not re-labelled.

Pooled n cells by state (sum of clusters with that label; timepoints not pooled for testing):

| cell_line | label | n_clusters | n_d0 | n_d3 | n_d7 | n_d10 | n_all |
|---|---|---|---|---|---|---|---|
| GM00731 | Fibroblast | 4 | 4939 | 5 | 14 | 1030 | 5988 |
| GM00731 | PartialReprog | 6 | 3 | 5220 | 578 | 31 | 5832 |
| GM00731 | EarlyPluripotency | 4 | 2 | 410 | 3917 | 85 | 4414 |
| GM00731 | Pluripotency | 1 | 0 | 5 | 354 | 80 | 439 |
| GM00731 | NonReprog | 4 | 65 | 348 | 1932 | 3114 | 5459 |
| GM23815 | Fibroblast | 5 | 7703 | 46 | 2233 | 1373 | 11355 |
| GM23815 | PartialReprog | 4 | 9 | 6113 | 309 | 50 | 6481 |
| GM23815 | EarlyPluripotency | 8 | 56 | 316 | 8172 | 4417 | 12961 |
| GM23815 | Pluripotency | 1 | 0 | 91 | 977 | 181 | 1249 |
| GM23815 | NonReprog | 1 | 0 | 10 | 8 | 1794 | 1812 |

Per state × timepoint:

| cell_line | day | label | n_cells | n_timepoint | frac_timepoint | ge50 |
|---|---|---|---|---|---|---|
| GM00731 | 0 | Fibroblast | 4939 | 5009 | +0.986 | True |
| GM00731 | 0 | PartialReprog | 3 | 5009 | +0.001 | False |
| GM00731 | 0 | EarlyPluripotency | 2 | 5009 | +0.000 | False |
| GM00731 | 0 | Pluripotency | 0 | 5009 | +0.000 | False |
| GM00731 | 0 | NonReprog | 65 | 5009 | +0.013 | True |
| GM00731 | 3 | Fibroblast | 5 | 5988 | +0.001 | False |
| GM00731 | 3 | PartialReprog | 5220 | 5988 | +0.872 | True |
| GM00731 | 3 | EarlyPluripotency | 410 | 5988 | +0.068 | True |
| GM00731 | 3 | Pluripotency | 5 | 5988 | +0.001 | False |
| GM00731 | 3 | NonReprog | 348 | 5988 | +0.058 | True |
| GM00731 | 7 | Fibroblast | 14 | 6795 | +0.002 | False |
| GM00731 | 7 | PartialReprog | 578 | 6795 | +0.085 | True |
| GM00731 | 7 | EarlyPluripotency | 3917 | 6795 | +0.576 | True |
| GM00731 | 7 | Pluripotency | 354 | 6795 | +0.052 | True |
| GM00731 | 7 | NonReprog | 1932 | 6795 | +0.284 | True |
| GM00731 | 10 | Fibroblast | 1030 | 4340 | +0.237 | True |
| GM00731 | 10 | PartialReprog | 31 | 4340 | +0.007 | False |
| GM00731 | 10 | EarlyPluripotency | 85 | 4340 | +0.020 | True |
| GM00731 | 10 | Pluripotency | 80 | 4340 | +0.018 | True |
| GM00731 | 10 | NonReprog | 3114 | 4340 | +0.718 | True |
| GM23815 | 0 | Fibroblast | 7703 | 7768 | +0.992 | True |
| GM23815 | 0 | PartialReprog | 9 | 7768 | +0.001 | False |
| GM23815 | 0 | EarlyPluripotency | 56 | 7768 | +0.007 | True |
| GM23815 | 0 | Pluripotency | 0 | 7768 | +0.000 | False |
| GM23815 | 0 | NonReprog | 0 | 7768 | +0.000 | False |
| GM23815 | 3 | Fibroblast | 46 | 6576 | +0.007 | False |
| GM23815 | 3 | PartialReprog | 6113 | 6576 | +0.930 | True |
| GM23815 | 3 | EarlyPluripotency | 316 | 6576 | +0.048 | True |
| GM23815 | 3 | Pluripotency | 91 | 6576 | +0.014 | True |
| GM23815 | 3 | NonReprog | 10 | 6576 | +0.002 | False |
| GM23815 | 7 | Fibroblast | 2233 | 11699 | +0.191 | True |
| GM23815 | 7 | PartialReprog | 309 | 11699 | +0.026 | True |
| GM23815 | 7 | EarlyPluripotency | 8172 | 11699 | +0.699 | True |
| GM23815 | 7 | Pluripotency | 977 | 11699 | +0.084 | True |
| GM23815 | 7 | NonReprog | 8 | 11699 | +0.001 | False |
| GM23815 | 10 | Fibroblast | 1373 | 7815 | +0.176 | True |
| GM23815 | 10 | PartialReprog | 50 | 7815 | +0.006 | True |
| GM23815 | 10 | EarlyPluripotency | 4417 | 7815 | +0.565 | True |
| GM23815 | 10 | Pluripotency | 181 | 7815 | +0.023 | True |
| GM23815 | 10 | NonReprog | 1794 | 7815 | +0.230 | True |

A timepoint qualifies for a two-state contrast iff both states have ≥50 cells. Timepoints are not pooled to reach the threshold.

| cell_line | day | state_a | state_b | n_a | n_b | min_cells | qualifies | primary | reason |
|---|---|---|---|---|---|---|---|---|---|
| GM00731 | 0 | PartialReprog | NonReprog | 3 | 65 | 50 | False | True | n<50 in at least one state |
| GM00731 | 0 | Fibroblast | NonReprog | 4939 | 65 | 50 | True | False | NA |
| GM00731 | 0 | PartialReprog | Fibroblast | 3 | 4939 | 50 | False | False | n<50 in at least one state |
| GM00731 | 3 | PartialReprog | NonReprog | 5220 | 348 | 50 | True | True | NA |
| GM00731 | 3 | Fibroblast | NonReprog | 5 | 348 | 50 | False | False | n<50 in at least one state |
| GM00731 | 3 | PartialReprog | Fibroblast | 5220 | 5 | 50 | False | False | n<50 in at least one state |
| GM00731 | 7 | PartialReprog | NonReprog | 578 | 1932 | 50 | True | True | NA |
| GM00731 | 7 | Fibroblast | NonReprog | 14 | 1932 | 50 | False | False | n<50 in at least one state |
| GM00731 | 7 | PartialReprog | Fibroblast | 578 | 14 | 50 | False | False | n<50 in at least one state |
| GM00731 | 10 | PartialReprog | NonReprog | 31 | 3114 | 50 | False | True | n<50 in at least one state |
| GM00731 | 10 | Fibroblast | NonReprog | 1030 | 3114 | 50 | True | False | NA |
| GM00731 | 10 | PartialReprog | Fibroblast | 31 | 1030 | 50 | False | False | n<50 in at least one state |
| GM23815 | 0 | PartialReprog | NonReprog | 9 | 0 | 50 | False | True | n<50 in at least one state |
| GM23815 | 0 | Fibroblast | NonReprog | 7703 | 0 | 50 | False | False | n<50 in at least one state |
| GM23815 | 0 | PartialReprog | Fibroblast | 9 | 7703 | 50 | False | False | n<50 in at least one state |
| GM23815 | 3 | PartialReprog | NonReprog | 6113 | 10 | 50 | False | True | n<50 in at least one state |
| GM23815 | 3 | Fibroblast | NonReprog | 46 | 10 | 50 | False | False | n<50 in at least one state |
| GM23815 | 3 | PartialReprog | Fibroblast | 6113 | 46 | 50 | False | False | n<50 in at least one state |
| GM23815 | 7 | PartialReprog | NonReprog | 309 | 8 | 50 | False | True | n<50 in at least one state |
| GM23815 | 7 | Fibroblast | NonReprog | 2233 | 8 | 50 | False | False | n<50 in at least one state |
| GM23815 | 7 | PartialReprog | Fibroblast | 309 | 2233 | 50 | True | False | NA |
| GM23815 | 10 | PartialReprog | NonReprog | 50 | 1794 | 50 | True | True | NA |
| GM23815 | 10 | Fibroblast | NonReprog | 1373 | 1794 | 50 | True | False | NA |
| GM23815 | 10 | PartialReprog | Fibroblast | 50 | 1373 | 50 | True | False | NA |

## Task 1 — sanity check (FINDINGS_MD2.md pooled MD means)

- n_PartialReprog=5832 (md2 json 5832)
- n_NonReprog=5459 (md2 json 5459)
- md_mean_PartialReprog=+0.150471 (md2 json +0.150471; FINDINGS_MD2.md printed +0.150)
- md_mean_NonReprog=+0.439448 (md2 json +0.439448; FINDINGS_MD2.md printed +0.439)
- n_match=True md_match_1e12=True printed_3dp_match=True ok=True

## Task 1 — cell-level means by state × timepoint

These are means of per-cell scores, not TMM among pseudobulk rows. Frozen ruler is per-cell log2-CPM (TMM nf=1). MD and TGF-β are md2 AddModuleScore.

| cell_line | day | label | n_cells | frozen_ruler_mean | frozen_ruler_median | md_score_mean | md_score_median | tgfb_score_mean | tgfb_score_median | age_up_minus_age_down_mean | age_up_minus_age_down_median | pluri_primary_mean | pluri_primary_median |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | 0 | Fibroblast | 4939 | +18.157 | +18.191 | +0.506 | +0.511 | +0.042 | +0.041 | +0.029 | +0.033 | +8.289 | +8.230 |
| GM00731 | 0 | PartialReprog | 3 | +17.175 | +16.875 | +0.291 | +0.239 | -0.028 | +0.011 | -0.023 | -0.031 | +9.762 | +10.069 |
| GM00731 | 0 | EarlyPluripotency | 2 | +16.674 | +16.674 | +0.249 | +0.249 | -0.021 | -0.021 | +0.065 | +0.065 | +11.542 | +11.542 |
| GM00731 | 0 | Pluripotency | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| GM00731 | 0 | NonReprog | 65 | +16.398 | +16.077 | +0.350 | +0.321 | +0.086 | +0.102 | +0.137 | +0.159 | +10.146 | +10.447 |
| GM00731 | 3 | Fibroblast | 5 | +18.083 | +18.047 | +0.378 | +0.473 | +0.023 | +0.002 | -0.049 | -0.028 | +10.472 | +9.554 |
| GM00731 | 3 | PartialReprog | 5220 | +17.394 | +17.447 | +0.135 | +0.110 | -0.047 | -0.053 | -0.099 | -0.104 | +12.280 | +12.369 |
| GM00731 | 3 | EarlyPluripotency | 410 | +16.744 | +16.858 | +0.050 | +0.015 | -0.024 | -0.027 | -0.007 | -0.008 | +14.842 | +14.947 |
| GM00731 | 3 | Pluripotency | 5 | +15.159 | +14.929 | +0.042 | +0.047 | +0.006 | +0.004 | +0.022 | +0.024 | +13.260 | +13.113 |
| GM00731 | 3 | NonReprog | 348 | +17.855 | +17.898 | +0.476 | +0.516 | +0.075 | +0.075 | +0.003 | +0.002 | +10.131 | +9.813 |
| GM00731 | 7 | Fibroblast | 14 | +17.953 | +17.891 | +0.425 | +0.392 | +0.064 | +0.066 | +0.012 | +0.041 | +9.669 | +9.663 |
| GM00731 | 7 | PartialReprog | 578 | +17.781 | +17.864 | +0.284 | +0.308 | -0.013 | -0.010 | -0.067 | -0.067 | +10.222 | +9.824 |
| GM00731 | 7 | EarlyPluripotency | 3917 | +17.098 | +17.059 | +0.152 | +0.143 | -0.026 | -0.025 | -0.001 | +0.002 | +11.916 | +12.037 |
| GM00731 | 7 | Pluripotency | 354 | +15.278 | +15.281 | +0.011 | -0.003 | +0.005 | +0.012 | +0.078 | +0.074 | +12.988 | +13.094 |
| GM00731 | 7 | NonReprog | 1932 | +18.154 | +18.232 | +0.445 | +0.453 | +0.048 | +0.045 | +0.026 | +0.025 | +9.269 | +9.143 |
| GM00731 | 10 | Fibroblast | 1030 | +18.402 | +18.425 | +0.434 | +0.454 | -0.009 | -0.012 | -0.029 | -0.019 | +8.858 | +8.681 |
| GM00731 | 10 | PartialReprog | 31 | +17.470 | +17.655 | +0.297 | +0.291 | -0.013 | -0.023 | -0.031 | -0.033 | +10.100 | +9.898 |
| GM00731 | 10 | EarlyPluripotency | 85 | +17.159 | +17.068 | +0.079 | +0.035 | -0.086 | -0.097 | -0.036 | -0.037 | +12.239 | +12.535 |
| GM00731 | 10 | Pluripotency | 80 | +15.746 | +15.632 | +0.069 | +0.050 | +0.019 | +0.029 | +0.070 | +0.076 | +12.298 | +12.335 |
| GM00731 | 10 | NonReprog | 3114 | +18.433 | +18.546 | +0.434 | +0.432 | +0.001 | -0.007 | +0.032 | +0.034 | +8.821 | +8.706 |
| GM23815 | 0 | Fibroblast | 7703 | +18.879 | +18.894 | +0.479 | +0.480 | +0.029 | +0.031 | +0.010 | +0.014 | +9.400 | +9.388 |
| GM23815 | 0 | PartialReprog | 9 | +18.585 | +18.602 | +0.316 | +0.301 | +0.006 | +0.014 | -0.103 | -0.126 | +10.116 | +10.178 |
| GM23815 | 0 | EarlyPluripotency | 56 | +17.686 | +17.598 | +0.379 | +0.360 | +0.077 | +0.064 | +0.121 | +0.113 | +10.870 | +11.152 |
| GM23815 | 0 | Pluripotency | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| GM23815 | 0 | NonReprog | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| GM23815 | 3 | Fibroblast | 46 | +18.849 | +18.929 | +0.385 | +0.400 | +0.081 | +0.088 | -0.049 | -0.053 | +10.345 | +9.938 |
| GM23815 | 3 | PartialReprog | 6113 | +18.184 | +18.176 | +0.085 | +0.070 | -0.036 | -0.043 | -0.114 | -0.121 | +12.175 | +12.118 |
| GM23815 | 3 | EarlyPluripotency | 316 | +18.059 | +18.100 | +0.059 | +0.003 | +0.017 | +0.002 | -0.017 | -0.007 | +15.050 | +16.595 |
| GM23815 | 3 | Pluripotency | 91 | +17.200 | +17.335 | +0.028 | +0.005 | -0.012 | -0.016 | +0.015 | +0.009 | +13.306 | +13.318 |
| GM23815 | 3 | NonReprog | 10 | +18.969 | +18.909 | +0.380 | +0.402 | +0.051 | +0.071 | -0.035 | -0.011 | +10.766 | +10.203 |
| GM23815 | 7 | Fibroblast | 2233 | +18.653 | +18.655 | +0.410 | +0.409 | -0.015 | -0.019 | +0.008 | +0.013 | +9.618 | +9.559 |
| GM23815 | 7 | PartialReprog | 309 | +18.544 | +18.559 | +0.174 | +0.193 | -0.058 | -0.065 | -0.047 | -0.047 | +11.409 | +10.775 |
| GM23815 | 7 | EarlyPluripotency | 8172 | +17.978 | +17.987 | +0.077 | +0.077 | -0.058 | -0.063 | -0.027 | -0.024 | +12.719 | +12.656 |
| GM23815 | 7 | Pluripotency | 977 | +17.032 | +17.032 | -0.059 | -0.068 | +0.006 | +0.005 | +0.031 | +0.025 | +13.641 | +13.751 |
| GM23815 | 7 | NonReprog | 8 | +18.615 | +18.550 | +0.407 | +0.399 | -0.004 | -0.029 | +0.031 | +0.043 | +9.334 | +9.405 |
| GM23815 | 10 | Fibroblast | 1373 | +18.894 | +18.912 | +0.379 | +0.382 | +0.037 | +0.036 | -0.036 | -0.025 | +9.680 | +9.652 |
| GM23815 | 10 | PartialReprog | 50 | +18.420 | +18.549 | +0.128 | +0.144 | -0.043 | -0.057 | -0.069 | -0.068 | +12.019 | +11.133 |
| GM23815 | 10 | EarlyPluripotency | 4417 | +18.398 | +18.452 | +0.166 | +0.183 | -0.023 | -0.028 | -0.041 | -0.039 | +11.742 | +11.604 |
| GM23815 | 10 | Pluripotency | 181 | +17.496 | +17.478 | -0.009 | -0.057 | +0.008 | +0.004 | +0.014 | +0.014 | +13.241 | +13.400 |
| GM23815 | 10 | NonReprog | 1794 | +19.189 | +19.192 | +0.423 | +0.423 | +0.035 | +0.035 | +0.020 | +0.024 | +9.368 | +9.340 |

## Task 1 — contrasts (per donor, per timepoint, per instrument)

Δ = mean(state_A) − mean(state_B) and median(state_A) − median(state_B). Permutation: 200 shuffles of the two-state label among cells in those two states at that donor × timepoint. Bootstrap: B=200 cells within state, seed `20260918`. Effect size: Cohen's d (pooled sd, n−1), same sign as Δmean. Primary pair: PartialReprog vs NonReprog. Context rows: Fibroblast vs NonReprog, PartialReprog vs Fibroblast. Rows with qualifies=False are not tested.

| cell_line | day | state_a | state_b | instrument | primary | qualifies | n_a | n_b | mean_a | mean_b | delta_mean | p_mean_lower | p_mean_two_sided | delta_mean_ci_lo | delta_mean_ci_hi | median_a | median_b | delta_median | p_median_lower | cohens_d | beats_null_lower | ok | reason |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | 0 | PartialReprog | NonReprog | frozen_ruler | True | False | 3 | 65 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 0 | PartialReprog | NonReprog | md_score | True | False | 3 | 65 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 0 | PartialReprog | NonReprog | tgfb_score | True | False | 3 | 65 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 0 | PartialReprog | NonReprog | age_up_minus_age_down | True | False | 3 | 65 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 0 | PartialReprog | NonReprog | pluri_primary | True | False | 3 | 65 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 0 | Fibroblast | NonReprog | frozen_ruler | False | True | 4939 | 65 | +18.157 | +16.398 | +1.758 | +1.000 | +0.005 | +1.413 | +2.065 | +18.191 | +16.077 | +2.114 | +1.000 | +2.756 | False | True | NA |
| GM00731 | 0 | Fibroblast | NonReprog | md_score | False | True | 4939 | 65 | +0.506 | +0.350 | +0.156 | +1.000 | +0.005 | +0.123 | +0.196 | +0.511 | +0.321 | +0.189 | +1.000 | +2.218 | False | True | NA |
| GM00731 | 0 | Fibroblast | NonReprog | tgfb_score | False | True | 4939 | 65 | +0.042 | +0.086 | -0.044 | +0.005 | +0.005 | -0.080 | -0.011 | +0.041 | +0.102 | -0.062 | +0.005 | -0.642 | True | True | NA |
| GM00731 | 0 | Fibroblast | NonReprog | age_up_minus_age_down | False | True | 4939 | 65 | +0.029 | +0.137 | -0.108 | +0.005 | +0.005 | -0.127 | -0.085 | +0.033 | +0.159 | -0.126 | +0.005 | -1.876 | True | True | NA |
| GM00731 | 0 | Fibroblast | NonReprog | pluri_primary | False | True | 4939 | 65 | +8.289 | +10.146 | -1.858 | +0.005 | +0.005 | -2.130 | -1.597 | +8.230 | +10.447 | -2.217 | +0.005 | -4.032 | True | True | NA |
| GM00731 | 0 | PartialReprog | Fibroblast | frozen_ruler | False | False | 3 | 4939 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 0 | PartialReprog | Fibroblast | md_score | False | False | 3 | 4939 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 0 | PartialReprog | Fibroblast | tgfb_score | False | False | 3 | 4939 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 0 | PartialReprog | Fibroblast | age_up_minus_age_down | False | False | 3 | 4939 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 0 | PartialReprog | Fibroblast | pluri_primary | False | False | 3 | 4939 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 3 | PartialReprog | NonReprog | frozen_ruler | True | True | 5220 | 348 | +17.394 | +17.855 | -0.462 | +0.005 | +0.005 | -0.571 | -0.362 | +17.447 | +17.898 | -0.451 | +0.005 | -0.577 | True | True | NA |
| GM00731 | 3 | PartialReprog | NonReprog | md_score | True | True | 5220 | 348 | +0.135 | +0.476 | -0.341 | +0.005 | +0.005 | -0.363 | -0.322 | +0.110 | +0.516 | -0.406 | +0.005 | -2.104 | True | True | NA |
| GM00731 | 3 | PartialReprog | NonReprog | tgfb_score | True | True | 5220 | 348 | -0.047 | +0.075 | -0.122 | +0.005 | +0.005 | -0.133 | -0.109 | -0.053 | +0.075 | -0.128 | +0.005 | -1.149 | True | True | NA |
| GM00731 | 3 | PartialReprog | NonReprog | age_up_minus_age_down | True | True | 5220 | 348 | -0.099 | +0.003 | -0.102 | +0.005 | +0.005 | -0.107 | -0.096 | -0.104 | +0.002 | -0.106 | +0.005 | -2.160 | True | True | NA |
| GM00731 | 3 | PartialReprog | NonReprog | pluri_primary | True | True | 5220 | 348 | +12.280 | +10.131 | +2.148 | +1.000 | +0.005 | +1.977 | +2.324 | +12.369 | +9.813 | +2.557 | +1.000 | +1.265 | False | True | NA |
| GM00731 | 3 | Fibroblast | NonReprog | frozen_ruler | False | False | 5 | 348 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 3 | Fibroblast | NonReprog | md_score | False | False | 5 | 348 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 3 | Fibroblast | NonReprog | tgfb_score | False | False | 5 | 348 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 3 | Fibroblast | NonReprog | age_up_minus_age_down | False | False | 5 | 348 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 3 | Fibroblast | NonReprog | pluri_primary | False | False | 5 | 348 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 3 | PartialReprog | Fibroblast | frozen_ruler | False | False | 5220 | 5 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 3 | PartialReprog | Fibroblast | md_score | False | False | 5220 | 5 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 3 | PartialReprog | Fibroblast | tgfb_score | False | False | 5220 | 5 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 3 | PartialReprog | Fibroblast | age_up_minus_age_down | False | False | 5220 | 5 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 3 | PartialReprog | Fibroblast | pluri_primary | False | False | 5220 | 5 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 7 | PartialReprog | NonReprog | frozen_ruler | True | True | 578 | 1932 | +17.781 | +18.154 | -0.373 | +0.005 | +0.005 | -0.447 | -0.306 | +17.864 | +18.232 | -0.369 | +0.005 | -0.424 | True | True | NA |
| GM00731 | 7 | PartialReprog | NonReprog | md_score | True | True | 578 | 1932 | +0.284 | +0.445 | -0.161 | +0.005 | +0.005 | -0.174 | -0.149 | +0.308 | +0.453 | -0.145 | +0.005 | -1.339 | True | True | NA |
| GM00731 | 7 | PartialReprog | NonReprog | tgfb_score | True | True | 578 | 1932 | -0.013 | +0.048 | -0.061 | +0.005 | +0.005 | -0.070 | -0.053 | -0.010 | +0.045 | -0.056 | +0.005 | -0.686 | True | True | NA |
| GM00731 | 7 | PartialReprog | NonReprog | age_up_minus_age_down | True | True | 578 | 1932 | -0.067 | +0.026 | -0.093 | +0.005 | +0.005 | -0.098 | -0.088 | -0.067 | +0.025 | -0.092 | +0.005 | -1.900 | True | True | NA |
| GM00731 | 7 | PartialReprog | NonReprog | pluri_primary | True | True | 578 | 1932 | +10.222 | +9.269 | +0.953 | +1.000 | +0.005 | +0.847 | +1.075 | +9.824 | +9.143 | +0.680 | +1.000 | +0.952 | False | True | NA |
| GM00731 | 7 | Fibroblast | NonReprog | frozen_ruler | False | False | 14 | 1932 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 7 | Fibroblast | NonReprog | md_score | False | False | 14 | 1932 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 7 | Fibroblast | NonReprog | tgfb_score | False | False | 14 | 1932 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 7 | Fibroblast | NonReprog | age_up_minus_age_down | False | False | 14 | 1932 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 7 | Fibroblast | NonReprog | pluri_primary | False | False | 14 | 1932 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 7 | PartialReprog | Fibroblast | frozen_ruler | False | False | 578 | 14 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 7 | PartialReprog | Fibroblast | md_score | False | False | 578 | 14 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 7 | PartialReprog | Fibroblast | tgfb_score | False | False | 578 | 14 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 7 | PartialReprog | Fibroblast | age_up_minus_age_down | False | False | 578 | 14 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 7 | PartialReprog | Fibroblast | pluri_primary | False | False | 578 | 14 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 10 | PartialReprog | NonReprog | frozen_ruler | True | False | 31 | 3114 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 10 | PartialReprog | NonReprog | md_score | True | False | 31 | 3114 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 10 | PartialReprog | NonReprog | tgfb_score | True | False | 31 | 3114 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 10 | PartialReprog | NonReprog | age_up_minus_age_down | True | False | 31 | 3114 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 10 | PartialReprog | NonReprog | pluri_primary | True | False | 31 | 3114 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 10 | Fibroblast | NonReprog | frozen_ruler | False | True | 1030 | 3114 | +18.402 | +18.433 | -0.031 | +0.149 | +0.289 | -0.087 | +0.021 | +18.425 | +18.546 | -0.121 | +0.005 | -0.036 | False | True | NA |
| GM00731 | 10 | Fibroblast | NonReprog | md_score | False | True | 1030 | 3114 | +0.434 | +0.434 | +0.000 | +0.577 | +0.940 | -0.008 | +0.008 | +0.454 | +0.432 | +0.022 | +1.000 | +0.004 | False | True | NA |
| GM00731 | 10 | Fibroblast | NonReprog | tgfb_score | False | True | 1030 | 3114 | -0.009 | +0.001 | -0.010 | +0.005 | +0.005 | -0.015 | -0.005 | -0.012 | -0.007 | -0.005 | +0.030 | -0.117 | True | True | NA |
| GM00731 | 10 | Fibroblast | NonReprog | age_up_minus_age_down | False | True | 1030 | 3114 | -0.029 | +0.032 | -0.060 | +0.005 | +0.005 | -0.064 | -0.057 | -0.019 | +0.034 | -0.054 | +0.005 | -1.071 | True | True | NA |
| GM00731 | 10 | Fibroblast | NonReprog | pluri_primary | False | True | 1030 | 3114 | +8.858 | +8.821 | +0.037 | +0.925 | +0.129 | -0.012 | +0.085 | +8.681 | +8.706 | -0.024 | +0.164 | +0.050 | False | True | NA |
| GM00731 | 10 | PartialReprog | Fibroblast | frozen_ruler | False | False | 31 | 1030 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 10 | PartialReprog | Fibroblast | md_score | False | False | 31 | 1030 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 10 | PartialReprog | Fibroblast | tgfb_score | False | False | 31 | 1030 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 10 | PartialReprog | Fibroblast | age_up_minus_age_down | False | False | 31 | 1030 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM00731 | 10 | PartialReprog | Fibroblast | pluri_primary | False | False | 31 | 1030 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 0 | PartialReprog | NonReprog | frozen_ruler | True | False | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 0 | PartialReprog | NonReprog | md_score | True | False | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 0 | PartialReprog | NonReprog | tgfb_score | True | False | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 0 | PartialReprog | NonReprog | age_up_minus_age_down | True | False | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 0 | PartialReprog | NonReprog | pluri_primary | True | False | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 0 | Fibroblast | NonReprog | frozen_ruler | False | False | 7703 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 0 | Fibroblast | NonReprog | md_score | False | False | 7703 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 0 | Fibroblast | NonReprog | tgfb_score | False | False | 7703 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 0 | Fibroblast | NonReprog | age_up_minus_age_down | False | False | 7703 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 0 | Fibroblast | NonReprog | pluri_primary | False | False | 7703 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 0 | PartialReprog | Fibroblast | frozen_ruler | False | False | 9 | 7703 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 0 | PartialReprog | Fibroblast | md_score | False | False | 9 | 7703 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 0 | PartialReprog | Fibroblast | tgfb_score | False | False | 9 | 7703 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 0 | PartialReprog | Fibroblast | age_up_minus_age_down | False | False | 9 | 7703 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 0 | PartialReprog | Fibroblast | pluri_primary | False | False | 9 | 7703 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 3 | PartialReprog | NonReprog | frozen_ruler | True | False | 6113 | 10 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 3 | PartialReprog | NonReprog | md_score | True | False | 6113 | 10 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 3 | PartialReprog | NonReprog | tgfb_score | True | False | 6113 | 10 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 3 | PartialReprog | NonReprog | age_up_minus_age_down | True | False | 6113 | 10 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 3 | PartialReprog | NonReprog | pluri_primary | True | False | 6113 | 10 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 3 | Fibroblast | NonReprog | frozen_ruler | False | False | 46 | 10 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 3 | Fibroblast | NonReprog | md_score | False | False | 46 | 10 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 3 | Fibroblast | NonReprog | tgfb_score | False | False | 46 | 10 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 3 | Fibroblast | NonReprog | age_up_minus_age_down | False | False | 46 | 10 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 3 | Fibroblast | NonReprog | pluri_primary | False | False | 46 | 10 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 3 | PartialReprog | Fibroblast | frozen_ruler | False | False | 6113 | 46 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 3 | PartialReprog | Fibroblast | md_score | False | False | 6113 | 46 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 3 | PartialReprog | Fibroblast | tgfb_score | False | False | 6113 | 46 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 3 | PartialReprog | Fibroblast | age_up_minus_age_down | False | False | 6113 | 46 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 3 | PartialReprog | Fibroblast | pluri_primary | False | False | 6113 | 46 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 7 | PartialReprog | NonReprog | frozen_ruler | True | False | 309 | 8 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 7 | PartialReprog | NonReprog | md_score | True | False | 309 | 8 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 7 | PartialReprog | NonReprog | tgfb_score | True | False | 309 | 8 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 7 | PartialReprog | NonReprog | age_up_minus_age_down | True | False | 309 | 8 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 7 | PartialReprog | NonReprog | pluri_primary | True | False | 309 | 8 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 7 | Fibroblast | NonReprog | frozen_ruler | False | False | 2233 | 8 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 7 | Fibroblast | NonReprog | md_score | False | False | 2233 | 8 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 7 | Fibroblast | NonReprog | tgfb_score | False | False | 2233 | 8 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 7 | Fibroblast | NonReprog | age_up_minus_age_down | False | False | 2233 | 8 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 7 | Fibroblast | NonReprog | pluri_primary | False | False | 2233 | 8 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<50 in at least one state |
| GM23815 | 7 | PartialReprog | Fibroblast | frozen_ruler | False | True | 309 | 2233 | +18.544 | +18.653 | -0.110 | +0.005 | +0.010 | -0.187 | -0.020 | +18.559 | +18.655 | -0.095 | +0.005 | -0.177 | True | True | NA |
| GM23815 | 7 | PartialReprog | Fibroblast | md_score | False | True | 309 | 2233 | +0.174 | +0.410 | -0.236 | +0.005 | +0.005 | -0.255 | -0.218 | +0.193 | +0.409 | -0.216 | +0.005 | -2.421 | True | True | NA |
| GM23815 | 7 | PartialReprog | Fibroblast | tgfb_score | False | True | 309 | 2233 | -0.058 | -0.015 | -0.044 | +0.005 | +0.005 | -0.053 | -0.031 | -0.065 | -0.019 | -0.046 | +0.005 | -0.527 | True | True | NA |
| GM23815 | 7 | PartialReprog | Fibroblast | age_up_minus_age_down | False | True | 309 | 2233 | -0.047 | +0.008 | -0.056 | +0.005 | +0.005 | -0.062 | -0.049 | -0.047 | +0.013 | -0.061 | +0.005 | -1.171 | True | True | NA |
| GM23815 | 7 | PartialReprog | Fibroblast | pluri_primary | False | True | 309 | 2233 | +11.409 | +9.618 | +1.791 | +1.000 | +0.005 | +1.642 | +1.972 | +10.775 | +9.559 | +1.217 | +1.000 | +2.222 | False | True | NA |
| GM23815 | 10 | PartialReprog | NonReprog | frozen_ruler | True | True | 50 | 1794 | +18.420 | +19.189 | -0.769 | +0.005 | +0.005 | -0.987 | -0.539 | +18.549 | +19.192 | -0.644 | +0.005 | -1.379 | True | True | NA |
| GM23815 | 10 | PartialReprog | NonReprog | md_score | True | True | 50 | 1794 | +0.128 | +0.423 | -0.295 | +0.005 | +0.005 | -0.331 | -0.257 | +0.144 | +0.423 | -0.279 | +0.005 | -3.700 | True | True | NA |
| GM23815 | 10 | PartialReprog | NonReprog | tgfb_score | True | True | 50 | 1794 | -0.043 | +0.035 | -0.078 | +0.005 | +0.005 | -0.105 | -0.053 | -0.057 | +0.035 | -0.092 | +0.005 | -0.924 | True | True | NA |
| GM23815 | 10 | PartialReprog | NonReprog | age_up_minus_age_down | True | True | 50 | 1794 | -0.069 | +0.020 | -0.089 | +0.005 | +0.005 | -0.107 | -0.069 | -0.068 | +0.024 | -0.091 | +0.005 | -2.087 | True | True | NA |
| GM23815 | 10 | PartialReprog | NonReprog | pluri_primary | True | True | 50 | 1794 | +12.019 | +9.368 | +2.651 | +1.000 | +0.005 | +2.201 | +3.116 | +11.133 | +9.340 | +1.793 | +1.000 | +5.464 | False | True | NA |
| GM23815 | 10 | Fibroblast | NonReprog | frozen_ruler | False | True | 1373 | 1794 | +18.894 | +19.189 | -0.294 | +0.005 | +0.005 | -0.326 | -0.263 | +18.912 | +19.192 | -0.280 | +0.005 | -0.567 | True | True | NA |
| GM23815 | 10 | Fibroblast | NonReprog | md_score | False | True | 1373 | 1794 | +0.379 | +0.423 | -0.043 | +0.005 | +0.005 | -0.050 | -0.038 | +0.382 | +0.423 | -0.040 | +0.005 | -0.503 | True | True | NA |
| GM23815 | 10 | Fibroblast | NonReprog | tgfb_score | False | True | 1373 | 1794 | +0.037 | +0.035 | +0.001 | +0.697 | +0.592 | -0.004 | +0.008 | +0.036 | +0.035 | +0.001 | +0.622 | +0.017 | False | True | NA |
| GM23815 | 10 | Fibroblast | NonReprog | age_up_minus_age_down | False | True | 1373 | 1794 | -0.036 | +0.020 | -0.056 | +0.005 | +0.005 | -0.059 | -0.052 | -0.025 | +0.024 | -0.049 | +0.005 | -1.150 | True | True | NA |
| GM23815 | 10 | Fibroblast | NonReprog | pluri_primary | False | True | 1373 | 1794 | +9.680 | +9.368 | +0.312 | +1.000 | +0.005 | +0.279 | +0.344 | +9.652 | +9.340 | +0.311 | +1.000 | +0.674 | False | True | NA |
| GM23815 | 10 | PartialReprog | Fibroblast | frozen_ruler | False | True | 50 | 1373 | +18.420 | +18.894 | -0.474 | +0.005 | +0.005 | -0.739 | -0.239 | +18.549 | +18.912 | -0.363 | +0.005 | -0.955 | True | True | NA |
| GM23815 | 10 | PartialReprog | Fibroblast | md_score | False | True | 50 | 1373 | +0.128 | +0.379 | -0.251 | +0.005 | +0.005 | -0.289 | -0.210 | +0.144 | +0.382 | -0.238 | +0.005 | -2.551 | True | True | NA |
| GM23815 | 10 | PartialReprog | Fibroblast | tgfb_score | False | True | 50 | 1373 | -0.043 | +0.037 | -0.080 | +0.005 | +0.005 | -0.108 | -0.052 | -0.057 | +0.036 | -0.093 | +0.005 | -1.024 | True | True | NA |
| GM23815 | 10 | PartialReprog | Fibroblast | age_up_minus_age_down | False | True | 50 | 1373 | -0.069 | -0.036 | -0.033 | +0.005 | +0.005 | -0.054 | -0.011 | -0.068 | -0.025 | -0.043 | +0.005 | -0.583 | True | True | NA |
| GM23815 | 10 | PartialReprog | Fibroblast | pluri_primary | False | True | 50 | 1373 | +12.019 | +9.680 | +2.338 | +1.000 | +0.005 | +1.855 | +2.837 | +11.133 | +9.652 | +1.482 | +1.000 | +3.868 | False | True | NA |

Primary pair, qualifying rows only (excerpt of the table above):

| cell_line | day | state_a | state_b | instrument | n_a | n_b | mean_a | mean_b | delta_mean | p_mean_lower | p_mean_two_sided | delta_mean_ci_lo | delta_mean_ci_hi | delta_median | p_median_lower | cohens_d | beats_null_lower |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | 3 | PartialReprog | NonReprog | frozen_ruler | 5220 | 348 | +17.394 | +17.855 | -0.462 | +0.005 | +0.005 | -0.571 | -0.362 | -0.451 | +0.005 | -0.577 | True |
| GM00731 | 3 | PartialReprog | NonReprog | md_score | 5220 | 348 | +0.135 | +0.476 | -0.341 | +0.005 | +0.005 | -0.363 | -0.322 | -0.406 | +0.005 | -2.104 | True |
| GM00731 | 3 | PartialReprog | NonReprog | tgfb_score | 5220 | 348 | -0.047 | +0.075 | -0.122 | +0.005 | +0.005 | -0.133 | -0.109 | -0.128 | +0.005 | -1.149 | True |
| GM00731 | 3 | PartialReprog | NonReprog | age_up_minus_age_down | 5220 | 348 | -0.099 | +0.003 | -0.102 | +0.005 | +0.005 | -0.107 | -0.096 | -0.106 | +0.005 | -2.160 | True |
| GM00731 | 3 | PartialReprog | NonReprog | pluri_primary | 5220 | 348 | +12.280 | +10.131 | +2.148 | +1.000 | +0.005 | +1.977 | +2.324 | +2.557 | +1.000 | +1.265 | False |
| GM00731 | 7 | PartialReprog | NonReprog | frozen_ruler | 578 | 1932 | +17.781 | +18.154 | -0.373 | +0.005 | +0.005 | -0.447 | -0.306 | -0.369 | +0.005 | -0.424 | True |
| GM00731 | 7 | PartialReprog | NonReprog | md_score | 578 | 1932 | +0.284 | +0.445 | -0.161 | +0.005 | +0.005 | -0.174 | -0.149 | -0.145 | +0.005 | -1.339 | True |
| GM00731 | 7 | PartialReprog | NonReprog | tgfb_score | 578 | 1932 | -0.013 | +0.048 | -0.061 | +0.005 | +0.005 | -0.070 | -0.053 | -0.056 | +0.005 | -0.686 | True |
| GM00731 | 7 | PartialReprog | NonReprog | age_up_minus_age_down | 578 | 1932 | -0.067 | +0.026 | -0.093 | +0.005 | +0.005 | -0.098 | -0.088 | -0.092 | +0.005 | -1.900 | True |
| GM00731 | 7 | PartialReprog | NonReprog | pluri_primary | 578 | 1932 | +10.222 | +9.269 | +0.953 | +1.000 | +0.005 | +0.847 | +1.075 | +0.680 | +1.000 | +0.952 | False |
| GM23815 | 10 | PartialReprog | NonReprog | frozen_ruler | 50 | 1794 | +18.420 | +19.189 | -0.769 | +0.005 | +0.005 | -0.987 | -0.539 | -0.644 | +0.005 | -1.379 | True |
| GM23815 | 10 | PartialReprog | NonReprog | md_score | 50 | 1794 | +0.128 | +0.423 | -0.295 | +0.005 | +0.005 | -0.331 | -0.257 | -0.279 | +0.005 | -3.700 | True |
| GM23815 | 10 | PartialReprog | NonReprog | tgfb_score | 50 | 1794 | -0.043 | +0.035 | -0.078 | +0.005 | +0.005 | -0.105 | -0.053 | -0.092 | +0.005 | -0.924 | True |
| GM23815 | 10 | PartialReprog | NonReprog | age_up_minus_age_down | 50 | 1794 | -0.069 | +0.020 | -0.089 | +0.005 | +0.005 | -0.107 | -0.069 | -0.091 | +0.005 | -2.087 | True |
| GM23815 | 10 | PartialReprog | NonReprog | pluri_primary | 50 | 1794 | +12.019 | +9.368 | +2.651 | +1.000 | +0.005 | +2.201 | +3.116 | +1.793 | +1.000 | +5.464 | False |

## Task 1 reading

### GM00731

Fired key: `both_lower`. n_qualifying=2.

d3: `both_lower`.

MD lower in PartialReprog than NonReprog and the frozen ruler lower too, both beating their nulls. Both instruments agree that partially reprogrammed cells read younger. Their claim is supported and our earlier negatives reflected the wrong test (a within-state trajectory with no baseline). GM00731 d3 PartialReprog vs NonReprog n_PR=5220 n_NR=348 MD Δmean=-0.34147022679488415 p_lower=0.004975124378109453 Cohen_d=-2.10407688578756 ruler Δmean=-0.4615273199823484 p_lower=0.004975124378109453 Cohen_d=-0.5769177381403222.

d7: `both_lower`.

MD lower in PartialReprog than NonReprog and the frozen ruler lower too, both beating their nulls. Both instruments agree that partially reprogrammed cells read younger. Their claim is supported and our earlier negatives reflected the wrong test (a within-state trajectory with no baseline). GM00731 d7 PartialReprog vs NonReprog n_PR=578 n_NR=1932 MD Δmean=-0.16091675133010636 p_lower=0.004975124378109453 Cohen_d=-1.3385320275469947 ruler Δmean=-0.3729116978678704 p_lower=0.004975124378109453 Cohen_d=-0.4240641770694765.

Sentences superseded (quoted, not edited in the source files):

> Ruler does not decline in the pooled PartialReprog state while MD does, both with their own nulls. The instruments disagree on exactly the population the published claim concerns. No winner is declared.

> MD moves per cluster, the frozen ruler does not. The instruments disagree on the same cells. What would settle it: an independent fibroblast age instrument scored on these same Louvain clusters, or the authors' Seurat/sctransform object with the frozen ruler projected onto it. No winner is declared.

### GM23815

Fired key: `both_lower`. n_qualifying=1.

d10: `both_lower`.

MD lower in PartialReprog than NonReprog and the frozen ruler lower too, both beating their nulls. Both instruments agree that partially reprogrammed cells read younger. Their claim is supported and our earlier negatives reflected the wrong test (a within-state trajectory with no baseline). GM23815 d10 PartialReprog vs NonReprog n_PR=50 n_NR=1794 MD Δmean=-0.2947491712855021 p_lower=0.004975124378109453 Cohen_d=-3.7001463221212694 ruler Δmean=-0.7686614248183865 p_lower=0.004975124378109453 Cohen_d=-1.378699803086083.

Sentences superseded (quoted, not edited in the source files):

> Ruler does not decline in the pooled PartialReprog state while MD does, both with their own nulls. The instruments disagree on exactly the population the published claim concerns. No winner is declared.

> MD moves per cluster, the frozen ruler does not. The instruments disagree on the same cells. What would settle it: an independent fibroblast age instrument scored on these same Louvain clusters, or the authors' Seurat/sctransform object with the frozen ruler projected onto it. No winner is declared.

Δmean < 0 with p_lower ≤ 0.05 is “lower / beating the null”. Donors are not pooled. Timepoints are not pooled. No winner is declared between instruments unless the pre-registered both_lower bullet fired, in which case they agree on direction; that is not a declaration that one instrument is better.

## Task 2 pre-registration (verbatim, written before any MD4 extrapolation distance or ρ)

Written **before any MD4 extrapolation distance or ρ**, 2026-09-19. Report-only. No gate.

For each donor × state × timepoint with n_cells ≥ 1 (n=0 is not a row): sum Louvain-kept UMIs to a pseudobulk on the frozen-ruler gene space (rulerY aligned as in md2/md3). TMM among that donor's state×timepoint rows, then frozen μ/σ/w. Distance to the GTEx fibroblast training z-space: PCA k=50 nearest-neighbour Euclidean and Mahalanobis, same functions as src/fibro2_extrap.py / FINDINGS_FIBRO2.md. Spearman ρ across those rows between each distance and the frozen ruler score, permutation null n_perm=200 seed 20260914 (shuffle the ruler score, hold distance), bootstrap CI B=200 seed 20260918. Donors reported separately and as concatenated rows; never averaged.

FINDINGS_FIBRO3.md Task 2 (the same PCA k=50 NN / Mahalanobis measures as FINDINGS_FIBRO2.md) found ruler score vs nn_euclidean ρ=-0.625. If PartialReprog sits further from training than NonReprog, that correlation would push the Task 1 frozen-ruler contrast toward “PR reads younger”. Reported as direction only. Not a gate. Not used to explain a Task 1 result away.


## Task 2 — extrapolation (report-only, no gate)

PCA k=50. Nearest neighbour is Euclidean in GTEx fibroblast training z-space. Mahalanobis is in the leading PCA of that z-space. Missing overlap genes at z=0. TMM among that donor's state×timepoint rows, then frozen μ/σ/w. Same measures as FINDINGS_FIBRO2.md / `src/fibro2_extrap.py`.

FINDINGS_FIBRO3.md Task 2 (same measures), quoted: `nn_euclidean: age score does not rise with distance ρ=-0.625 CI=[-0.834, -0.320] n=22. Not a gate.`

| cell_line | day | label | n_cells | age_score | nn_euclidean | mahalanobis_pca | pca_k |
|---|---|---|---|---|---|---|---|
| GM00731 | 0 | Fibroblast | 4939 | +5.253 | +465.617 | +24.411 | 50 |
| GM00731 | 0 | PartialReprog | 3 | +1.397 | +1342.161 | +64.592 | 50 |
| GM00731 | 0 | EarlyPluripotency | 2 | -4.307 | +2192.731 | +137.899 | 50 |
| GM00731 | 0 | NonReprog | 65 | -0.077 | +670.751 | +26.971 | 50 |
| GM00731 | 3 | Fibroblast | 5 | +10.582 | +1215.554 | +56.410 | 50 |
| GM00731 | 3 | PartialReprog | 5220 | +2.337 | +479.017 | +23.982 | 50 |
| GM00731 | 3 | EarlyPluripotency | 410 | +1.927 | +511.647 | +24.076 | 50 |
| GM00731 | 3 | Pluripotency | 5 | -18.978 | +1553.696 | +80.698 | 50 |
| GM00731 | 3 | NonReprog | 348 | +8.122 | +476.754 | +22.042 | 50 |
| GM00731 | 7 | Fibroblast | 14 | +8.596 | +918.913 | +38.900 | 50 |
| GM00731 | 7 | PartialReprog | 578 | +4.737 | +468.982 | +23.158 | 50 |
| GM00731 | 7 | EarlyPluripotency | 3917 | -0.739 | +503.349 | +24.453 | 50 |
| GM00731 | 7 | Pluripotency | 354 | -8.138 | +616.361 | +24.829 | 50 |
| GM00731 | 7 | NonReprog | 1932 | +7.015 | +470.259 | +24.100 | 50 |
| GM00731 | 10 | Fibroblast | 1030 | +7.873 | +474.154 | +24.294 | 50 |
| GM00731 | 10 | PartialReprog | 31 | +4.933 | +709.005 | +31.008 | 50 |
| GM00731 | 10 | EarlyPluripotency | 85 | -0.190 | +531.941 | +24.800 | 50 |
| GM00731 | 10 | Pluripotency | 80 | -5.242 | +631.384 | +25.461 | 50 |
| GM00731 | 10 | NonReprog | 3114 | +9.474 | +485.365 | +24.012 | 50 |
| GM23815 | 0 | Fibroblast | 7703 | +0.241 | +449.237 | +23.054 | 50 |
| GM23815 | 0 | PartialReprog | 9 | +2.361 | +1026.072 | +47.539 | 50 |
| GM23815 | 0 | EarlyPluripotency | 56 | -5.297 | +854.835 | +35.877 | 50 |
| GM23815 | 3 | Fibroblast | 46 | +5.947 | +593.699 | +25.444 | 50 |
| GM23815 | 3 | PartialReprog | 6113 | -2.952 | +455.379 | +23.014 | 50 |
| GM23815 | 3 | EarlyPluripotency | 316 | -1.037 | +514.939 | +24.337 | 50 |
| GM23815 | 3 | Pluripotency | 91 | -13.423 | +642.148 | +24.855 | 50 |
| GM23815 | 3 | NonReprog | 10 | -0.790 | +1066.480 | +48.430 | 50 |
| GM23815 | 7 | Fibroblast | 2233 | -1.817 | +465.875 | +23.737 | 50 |
| GM23815 | 7 | PartialReprog | 309 | -2.809 | +482.403 | +24.240 | 50 |
| GM23815 | 7 | EarlyPluripotency | 8172 | -8.203 | +508.372 | +25.420 | 50 |
| GM23815 | 7 | Pluripotency | 977 | -12.384 | +627.900 | +25.657 | 50 |
| GM23815 | 7 | NonReprog | 8 | -5.052 | +1515.850 | +80.299 | 50 |
| GM23815 | 10 | Fibroblast | 1373 | +0.181 | +468.409 | +24.716 | 50 |
| GM23815 | 10 | PartialReprog | 50 | -1.208 | +566.537 | +26.410 | 50 |
| GM23815 | 10 | EarlyPluripotency | 4417 | -3.487 | +484.765 | +24.461 | 50 |
| GM23815 | 10 | Pluripotency | 181 | -10.192 | +682.821 | +27.762 | 50 |
| GM23815 | 10 | NonReprog | 1794 | +2.866 | +480.957 | +24.440 | 50 |

Spearman ρ across state × timepoint rows (shuffle ruler score, hold distance). `rho_p` is p_greater = (n_null ≥ ρ_obs + 1)/(n_null+1), not p_lower. A negative ρ with large p_greater is on the 'ruler falls with distance' side; the CI is the bootstrap on ρ.

| scope | cell_line | distance | n | rho | rho_null | rho_p | rho_ci_lo | rho_ci_hi | note |
|---|---|---|---|---|---|---|---|---|---|
| all_donor_state_timepoint_n_ge_1 | NA | nn_euclidean | 37 | -0.280 | +0.026 | +0.955 | -0.576 | +0.065 | permutation shuffles frozen ruler score, holds distance; unit=state×timepoint |
| per_donor_GM00731 | GM00731 | nn_euclidean | 19 | -0.368 | +0.055 | +0.945 | -0.720 | +0.197 | permutation shuffles frozen ruler score, holds distance; unit=state×timepoint |
| per_donor_GM23815 | GM23815 | nn_euclidean | 18 | -0.273 | -0.013 | +0.871 | -0.682 | +0.184 | permutation shuffles frozen ruler score, holds distance; unit=state×timepoint |
| all_donor_state_timepoint_n_ge_1 | NA | mahalanobis_pca | 37 | -0.292 | +0.014 | +0.970 | -0.625 | -0.018 | permutation shuffles frozen ruler score, holds distance; unit=state×timepoint |
| per_donor_GM00731 | GM00731 | mahalanobis_pca | 19 | -0.370 | +0.018 | +0.950 | -0.759 | +0.209 | permutation shuffles frozen ruler score, holds distance; unit=state×timepoint |
| per_donor_GM23815 | GM23815 | mahalanobis_pca | 18 | -0.168 | -0.007 | +0.756 | -0.550 | +0.304 | permutation shuffles frozen ruler score, holds distance; unit=state×timepoint |

Direction this would push Task 1 (not used to explain Task 1 away):

| cell_line | day | ok | n_PR | n_NR | nn_PR | nn_NR | nn_PR_minus_NR | mahal_PR | mahal_NR | mahal_PR_minus_NR | age_PR | age_NR | age_PR_minus_NR | PR_farther_nn | PR_farther_mahal | would_push_ruler_PR_lower_via_nn | reason |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | 0 | True | +3.000 | +65.000 | +1342.161 | +670.751 | +671.410 | +64.592 | +26.971 | +37.621 | +1.397 | -0.077 | +1.475 | True | True | True | NA |
| GM00731 | 3 | True | +5220.000 | +348.000 | +479.017 | +476.754 | +2.263 | +23.982 | +22.042 | +1.940 | +2.337 | +8.122 | -5.785 | True | True | True | NA |
| GM00731 | 7 | True | +578.000 | +1932.000 | +468.982 | +470.259 | -1.277 | +23.158 | +24.100 | -0.942 | +4.737 | +7.015 | -2.278 | False | False | False | NA |
| GM00731 | 10 | True | +31.000 | +3114.000 | +709.005 | +485.365 | +223.640 | +31.008 | +24.012 | +6.996 | +4.933 | +9.474 | -4.541 | True | True | True | NA |
| GM23815 | 0 | False | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | missing PartialReprog or NonReprog row |
| GM23815 | 3 | True | +6113.000 | +10.000 | +455.379 | +1066.480 | -611.101 | +23.014 | +48.430 | -25.416 | -2.952 | -0.790 | -2.162 | False | False | False | NA |
| GM23815 | 7 | True | +309.000 | +8.000 | +482.403 | +1515.850 | -1033.446 | +24.240 | +80.299 | -56.059 | -2.809 | -5.052 | +2.243 | False | False | False | NA |
| GM23815 | 10 | True | +50.000 | +1794.000 | +566.537 | +480.957 | +85.580 | +26.410 | +24.440 | +1.970 | -1.208 | +2.866 | -4.074 | True | True | True | NA |

At Task 1 qualifying primary timepoints (from `t2_direction.csv`, not reconstructed):

| cell_line | day | n_PR | n_NR | nn_PR | nn_NR | nn_PR_minus_NR | mahal_PR | mahal_NR | mahal_PR_minus_NR | age_PR | age_NR | age_PR_minus_NR | PR_farther_nn | would_push_ruler_PR_lower_via_nn |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | 3 | +5220.000 | +348.000 | +479.017 | +476.754 | +2.263 | +23.982 | +22.042 | +1.940 | +2.337 | +8.122 | -5.785 | True | True |
| GM00731 | 7 | +578.000 | +1932.000 | +468.982 | +470.259 | -1.277 | +23.158 | +24.100 | -0.942 | +4.737 | +7.015 | -2.278 | False | False |
| GM23815 | 10 | +50.000 | +1794.000 | +566.537 | +480.957 | +85.580 | +26.410 | +24.440 | +1.970 | -1.208 | +2.866 | -4.074 | True | True |

- n_rows=37 pca_k=50 ρ_nn_all=-0.280
- Report-only. If ρ(distance, ruler)<0 and PartialReprog sits farther from training than NonReprog, that would push the Task 1 frozen-ruler contrast toward PR reading younger. Not a gate. Not used to discount Task 1.

If PartialReprog sits farther from the GTEx training distribution than NonReprog, and ruler score falls with distance, that would push the Task 1 frozen-ruler contrast toward PartialReprog reading younger. Reported alongside Task 1. Not a gate.

## Task 3 — withdrawal of FINDINGS_MD3.md Task 1 (dated 2026-09-19)

Dated 2026-09-19. FINDINGS_MD3.md Task 1 status `instruments_disagree_on_claim_population` rested on n₀=3 PartialReprog cells at day 0 and is withdrawn as unsupported, superseded by MD4 Task 1. FINDINGS_MD3.md is not edited. Task 2 of MD3 (the Δρ result) is unaffected and stands.


Quoted from FINDINGS_MD3.md, not edited there:

> Ruler does not decline in the pooled PartialReprog state while MD does, both with their own nulls. The instruments disagree on exactly the population the published claim concerns. No winner is declared.

Status line quoted: Task 1: `instruments_disagree_on_claim_population`.

Numbers quoted: Pooled PartialReprog ruler: n_ok=2 n_pass=0. Pooled PartialReprog MD: n_ok=2 n_pass=1 hits d0→d7 decline=+0.007 p=+0.015 n_cells_d0=3 n_cells_end=578. Frozen ruler d0→d7 decline=-3.339 p=+0.403.

Those numbers rest on n₀=3 PartialReprog cells at day 0. The status is withdrawn as unsupported. It is superseded by Task 1 of this file, which compares PartialReprog vs NonReprog at timepoints where both exist with ≥50 cells.

FINDINGS_MD3.md Task 2 (paired Δρ vs donor age on independent fibroblast cohorts) is unaffected and stands. Quoted:

> Δρ CI excludes zero in favour of the ruler in 2 of 3 named cohorts. MD and/or their aging signature predict donor age materially worse than the ruler on independent fibroblast cohorts.

## What this changes in FINDINGS_MD3.md / FINDINGS_MD2.md / FINDINGS_FIBRO.md (quote, not edited)

This file does not edit `FINDINGS_MD3.md`, `FINDINGS_MD2.md`, `FINDINGS_FIBRO.md`, `FINDINGS_FIBRO2.md`, `FINDINGS_FIBRO3.md`, or `FALSIFICATION.md`.

FINDINGS_MD3.md Task 1 status and reading, quoted, withdrawn as unsupported by this file:

> Ruler does not decline in the pooled PartialReprog state while MD does, both with their own nulls. The instruments disagree on exactly the population the published claim concerns. No winner is declared.

FINDINGS_MD3.md Task 2 Δρ reading, quoted, stands:

> Δρ CI excludes zero in favour of the ruler in 2 of 3 named cohorts. MD and/or their aging signature predict donor age materially worse than the ruler on independent fibroblast cohorts.

FINDINGS_MD2.md Task 1 reading, quoted, is a within-cluster d0→dT test and is not this contrast:

> No cluster at any resolution shows p ≤ 0.05. The per-cluster negative holds and is not a resolution artifact. The frozen ruler does not move in these cells at any granularity tested.

FINDINGS_MD2.md Task 2 reading, quoted, was also a within-cluster trajectory test (pairable Louvain clusters at both endpoints; PartialReprog among them: n=0):

> MD moves per cluster, the frozen ruler does not. The instruments disagree on the same cells. What would settle it: an independent fibroblast age instrument scored on these same Louvain clusters, or the authors' Seurat/sctransform object with the frozen ruler projected onto it. No winner is declared.

FINDINGS_FIBRO.md Stage 2 verdict sentence, quoted. That d0→d10 all-cell verdict stands. Nothing here promotes d7 to the primary Stage 2 result:

> age score does not decline → the ruler reads a static donor property, not a modifiable state. Step 3 is not supported by this data.

## Limitations

1. Two donors in GSE297234.
2. Clusters are not lineage-tracked. PartialReprog / NonReprog labels are argmax of mmc3 Reprog_cell_state_signatures AddModuleScore, not Slingshot trajectories from day-0 fibroblasts. States are inferred from expression, not lineage.
3. A state contrast at one timepoint is not a trajectory. This file does not revive the MD3 d0→d7 / d0→d10 within-PartialReprog test.
4. Louvain is a Python deviation from their R/sctransform pipeline, as recorded in FINDINGS_MD2.md (LogNormalize + quadratic HVG + percent.mt residualization + PCA 50 + SNN + networkx louvain, resolution=0.8 Seurat default; STAR Methods omit resolution).
5. Our reproduction of their metrics is ours, not theirs. AddModuleScore is a Python reimplementation of the published algorithm, not Seurat's C++/R object.
6. Frozen ruler trained on GTEx V10 cultured fibroblasts, public `AGE` 10-year bins. Per-cell scores use log2-CPM with TMM nf=1; TMM among 20k cells was not fit. That is a different size-factor treatment from the cluster×timepoint TMM in md2/md3.
7. Cross-platform shift from GTEx bulk polyA (RNASeQCv2.4.2) to 10x 3' scRNA-seq.
8. Missing overlap genes at z=0. Nothing fitted on GSE297234.
9. MIN_CELLS_CONTRAST=50. Timepoints below that on either side are not tested.
10. Seed `20260914`. Bootstrap seed `20260918`. n_perm=200, n_boot=200, n_random=200.
11. GSE325735 not opened. d7 is not the Stage 2 endpoint.
12. Table S3 Age up/down are the published lists, not a rebuilt DESeq2 analysis.
13. Task 2 distances are on TMM state×timepoint pseudobulks; Task 1 ruler scores are per-cell log2-CPM. Those two ruler numbers are not interchangeable.

## Open questions

1. Would the authors' sctransform v2 object change which Louvain clusters are PartialReprog vs NonReprog at these timepoints?
2. If cells were lineage-tracked from day-0 fibroblasts, would the same contrast hold?
3. An independent fibroblast age instrument on these same labelled cells at d3/d7 (aged) and d10 (young) would settle an instrument disagreement if one remains.

## Gene-list versions

- MD n=205 TGFB n=54 Age up n=1533 Age down n=2007
- mmc3 equals md2 copies: up=True down=True
- used_as_deseq2_rebuild=False used_as_published_lists=True
- AddModuleScore nbin=24 ctrl=100

- GM00731 Age up mapped 896/1533 missing=636; Age down mapped 1586/2007 missing=421 (from `t1_s3_ams_meta_GM00731.json`)
- GM23815 Age up mapped 894/1533 missing=638; Age down mapped 1582/2007 missing=425 (from `t1_s3_ams_meta_GM23815.json`)

## ID types actually read

- **louvain_GM00731_symbols** kind=`symbol` n=29058 n_ensembl=0 n_refseq=0 n_symbol_like=29057 source=`<repo>\data\processed\md2\louvain_counts_GM00731.npz:symbols` examples=['AL627309.1', 'AL627309.5', 'AP006222.2', 'AL669831.2', 'LINC01409', 'FAM87B', 'LINC01128', 'LINC00115', 'FAM41C', 'AL645608.6', 'AL645608.2', 'AL645608.4']
- **louvain_GM00731_gene_id** kind=`ensembl` n=29058 n_ensembl=29058 n_refseq=0 n_symbol_like=0 source=`<repo>\data\processed\md2\louvain_counts_GM00731.npz:gene_id` examples=['ENSG00000238009', 'ENSG00000241860', 'ENSG00000286448', 'ENSG00000229905', 'ENSG00000237491', 'ENSG00000177757', 'ENSG00000228794', 'ENSG00000225880', 'ENSG00000230368', 'ENSG00000272438', 'ENSG00000230699', 'ENSG00000241180']
- **louvain_GM23815_symbols** kind=`symbol` n=29154 n_ensembl=0 n_refseq=0 n_symbol_like=29154 source=`<repo>\data\processed\md2\louvain_counts_GM23815.npz:symbols` examples=['AL627309.1', 'AL627309.3', 'AL627309.5', 'AP006222.2', 'LINC01409', 'FAM87B', 'LINC01128', 'LINC00115', 'FAM41C', 'AL645608.6', 'AL645608.2', 'AL645608.4']
- **louvain_GM23815_gene_id** kind=`ensembl` n=29154 n_ensembl=29154 n_refseq=0 n_symbol_like=0 source=`<repo>\data\processed\md2\louvain_counts_GM23815.npz:gene_id` examples=['ENSG00000238009', 'ENSG00000239945', 'ENSG00000241860', 'ENSG00000286448', 'ENSG00000237491', 'ENSG00000177757', 'ENSG00000228794', 'ENSG00000225880', 'ENSG00000230368', 'ENSG00000272438', 'ENSG00000230699', 'ENSG00000241180']

## Columns actually read

- **md2_genesets_json** source=`<repo>\results\md2\genesets.json` columns=['MD', 'TGFB', 'EMT', 'reprog', 'mmc3_MD', 'mmc3_TGFB', 'mmc3_age_up', 'mmc3_age_down', 'msigdb_version']
- **md3_genesets_json** source=`<repo>\results\md3\genesets.json` columns=['MD', 'TGFB', 'age_up', 'age_down', 'reprog', 'mmc3_path', 'msigdb_version']
- **t2_cluster_labels_GM00731** source=`<repo>\results\md2\t2_cluster_labels_GM00731.csv` columns=['cell_line', 'cluster', 'n_cells', 'label', 'mean_Fibroblast', 'mean_PartialReprog', 'mean_EarlyPluripotency', 'mean_Pluripotency', 'mean_NonReprog']
- **louvain_obs_GM00731** source=`<repo>\results\md2\louvain_obs_GM00731.csv` columns=['barcode', 'gsm', 'file', 'cell_line', 'day', 'age_years', 'age_source', 'umi', 'n_genes', 'mito_frac', 'cluster', 'louvain_kept', 'percent_mt']
- **louvain_ams_GM00731** source=`<repo>\data\processed\md2\louvain_ams_GM00731.npz` columns=['md', 'tgfb', 'bins', 'gene_mean', 'cluster', 'day', 'n_genes']
- **t2_cluster_labels_GM23815** source=`<repo>\results\md2\t2_cluster_labels_GM23815.csv` columns=['cell_line', 'cluster', 'n_cells', 'label', 'mean_Fibroblast', 'mean_PartialReprog', 'mean_EarlyPluripotency', 'mean_Pluripotency', 'mean_NonReprog']
- **louvain_obs_GM23815** source=`<repo>\results\md2\louvain_obs_GM23815.csv` columns=['barcode', 'gsm', 'file', 'cell_line', 'day', 'age_years', 'age_source', 'umi', 'n_genes', 'mito_frac', 'cluster', 'louvain_kept', 'percent_mt']
- **louvain_ams_GM23815** source=`<repo>\data\processed\md2\louvain_ams_GM23815.npz` columns=['md', 'tgfb', 'bins', 'gene_mean', 'cluster', 'day']
- **md3_t1_occupancy_pooled** source=`<repo>\results\md3\t1_occupancy_pooled.csv` columns=['cell_line', 'label', 'n_clusters', 'n_d0', 'n_d3', 'n_d7', 'n_d10']
- **t2_reproduction** source=`<repo>\results\md2\t2_reproduction.json` columns=['n_PartialReprog', 'n_NonReprog', 'md_mean_PartialReprog', 'md_mean_NonReprog', 'md_Partial_minus_NonReprog', 'n_cells_by_label', 'min_cells', 'has_both_labels_ge_min', 'md_partial_lower', 'matches_reported_direction', 'criterion']
- **allcell_obs_GM00731** source=`<repo>\results\md2\allcell_obs_GM00731.csv` columns=['barcode', 'gsm', 'file', 'cell_line', 'day', 'age_years', 'age_source', 'umi', 'n_genes', 'mito_frac']
- **allcell_obs_GM23815** source=`<repo>\results\md2\allcell_obs_GM23815.csv` columns=['barcode', 'gsm', 'file', 'cell_line', 'day', 'age_years', 'age_source', 'umi', 'n_genes', 'mito_frac']
- **t2_gtex_obs** source=`<repo>\results\fibro\stage1_obs.csv` columns=['SAMPID', 'SMATSSCR', 'SMCENTER', 'SMPTHNTS', 'SMRIN', 'SMTS', 'SMTSD', 'SMUBRID', 'SMTSISCH', 'SMTSPAX', 'SMNABTCH', 'SMNABTCHT', 'SMNABTCHD', 'SMGEBTCH', 'SMGEBTCHD', 'SMGEBTCHT', 'ANALYTE_TYPE', 'SMAFRZE', 'SMGTC', 'SMRDTTL', 'SMALTTL', 'SMALTALG', 'SMSUPALG', 'SMRDLGTH', 'SMVQCFL', 'SMLMAPQ', 'SMUMPRD', 'SMUNPDRD', 'SMMPPD', 'SMMAPRT', 'SMMPPDUN', 'SMUNMPRT', 'SMMPDP', 'SMDPMPRT', 'SMMPPDXG', 'SMMPDPXG', 'SMDPRTXG', 'SMCHMRD', 'SMCHMRT', 'SMMPPDPR', 'SMMPHQRD', 'SMMPHQRT', 'SMMPLQRD', 'SMSPLTRT', 'SME1MPRD', 'SME2MPRD', 'SME1MPRT', 'SME2MPRT', 'SME1MMB', 'SME2MMB', 'SME1TTLB', 'SME2TTLB', 'SME1MMRT', 'SME2MMRT', 'SMTTLMM', 'SMTTLB', 'SMBSMMRT', 'SMESTLBS', 'SMEXNCRD', 'SMEXNCRT', 'SMEXPEFF', 'SMNTRNRD', 'SMNTRNRT', 'SMNTRARD', 'SMNTRART', 'SMNTERRD', 'SMNTERRT', 'SMAMBRD', 'SMAMBRT', 'SMNTEXC', 'SMDSCRT', 'SMEXNCRTHQ', 'SMNTRNRTHQ', 'SMNTRARTHQ', 'SMNTERRTHQ', 'SMAMBRTHQ', 'SME1SNSE', 'SME2SNSE', 'SME1ANTI', 'SME2ANTI', 'SME1PCTS', 'SME2PCTS', 'SMGNSDTC', 'SMRRNARD', 'SMRRNART', 'SMMFLGTH', 'SMSFLGTH', 'SMMDFLGTH', 'SMSMFLGTH', 'SMFGCMN', 'SMFGCSD', 'SMFGCSK', 'SMFGCKT', 'SM3PBMN', 'SM3PBSD', 'SM3PBMD', 'SM3PB25P', 'SM3PB75P', 'SM3PBSDM', 'SM3PBGN', 'SMMDMNCV', 'SMMDCVSD', 'SMMDCVCV', 'SMEXCVMD', 'SMEXCVMAD', 'SMMNCV', 'SMUVCRD', 'SMUVCRT', 'SMSHRTRD', 'SMSHRTRT', 'SMSMRDHQ', 'SMSMRTHQ', 'SMPRERDHQ', 'SMPRERTHQ', 'SMSMGNDT', 'SMPREGNDT', 'SMRDLNMN', 'SMRDLNMD', 'SMRDLNSD', 'donor', 'SUBJID', 'SEX', 'AGE', 'DTHHRDY', 'age_mid', 'age_ordinal']
- **t2_louvain_obs_GM00731** source=`<repo>\results\md2\louvain_obs_GM00731.csv` columns=['barcode', 'gsm', 'file', 'cell_line', 'day', 'age_years', 'age_source', 'umi', 'n_genes', 'mito_frac', 'cluster', 'louvain_kept', 'percent_mt']
- **t2_labels_GM00731** source=`<repo>\results\md2\t2_cluster_labels_GM00731.csv` columns=['cell_line', 'cluster', 'n_cells', 'label', 'mean_Fibroblast', 'mean_PartialReprog', 'mean_EarlyPluripotency', 'mean_Pluripotency', 'mean_NonReprog']
- **t2_louvain_obs_GM23815** source=`<repo>\results\md2\louvain_obs_GM23815.csv` columns=['barcode', 'gsm', 'file', 'cell_line', 'day', 'age_years', 'age_source', 'umi', 'n_genes', 'mito_frac', 'cluster', 'louvain_kept', 'percent_mt']
- **t2_labels_GM23815** source=`<repo>\results\md2\t2_cluster_labels_GM23815.csv` columns=['cell_line', 'cluster', 'n_cells', 'label', 'mean_Fibroblast', 'mean_PartialReprog', 'mean_EarlyPluripotency', 'mean_Pluripotency', 'mean_NonReprog']

## Gene-set names actually read

- **MD_built_from_md2** name=`MD (STAR Methods built; md2)` n=205 id_type=symbol source=`<repo>\results\md2\genesets.json`
- **TGFB_built_from_md2** name=`TGFB (STAR Methods built; md2)` n=54 id_type=symbol source=`<repo>\results\md2\genesets.json`
- **mmc3_Aging_signatures_Age up** name=`Aging_signatures:Age up` n=1533 id_type=symbol source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx via md3/genesets.json`
- **mmc3_Aging_signatures_Age down** name=`Aging_signatures:Age down` n=2007 id_type=symbol source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx via md3/genesets.json`

## Files

- `src/md4_common.py`, `md4_task1.py`, `md4_task2.py`, `md4_findings.py`, `md4_run.py`
- `results/md4/`
- `FINDINGS_MD4.md`
- `PROGRESS_MD4.md`

## Supersession note, appended 2026-10-02: atlas ruler − MD Δρ on 93 pseudobulks

Appended; nothing above is changed. Where this file says that the FINDINGS_MD3.md Task 2 Δρ result stands, its fibroblast-atlas ruler − MD interval no longer stands.

The ruler − MD Δρ computed with each of the 93 CELLxGENE atlas pseudobulks (a19d1667) as a unit (+0.202, 95% CI +0.006 to +0.440, `excludes_zero=True`; `results/md3/t2_delta_rho.csv`, row `cxg_ruler_minus_MD`) is **superseded** by the per-person re-analysis, in which pseudobulk scores are averaged within each of 65 donors (+0.171, 95% CI −0.118 to +0.412, `excludes_zero=False`; `results/paper_figs/cxg_true_donor_results.csv`, `analysis=true_donor`, `contrast=ruler_minus_MD_AddModuleScore`). The manuscript reports the per-person value in Table 2 (`paper/rejuvenation_readouts_preprint_11.docx`). **The per-person re-analysis was not pre-registered:** `results/paper_figs/` has no PREREG flag, and its results were written on 2026-09-23, after the pre-registered 93-pseudobulk result was known. After regrouping, the ruler − MD interval excludes zero only in GSE226189 (+0.467, 95% CI +0.239 to +0.696). The atlas ruler − (age-up − age-down) interval excludes zero in both units (+0.295, CI +0.107 to +0.487 on 93 pseudobulks; +0.372, CI +0.094 to +0.595 on 65 donors). The pre-registration flags are not edited. The reading is listed as B7 in `results/verify/WITHDRAWN_READINGS.md`.

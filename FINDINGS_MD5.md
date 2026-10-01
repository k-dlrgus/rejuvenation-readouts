# FINDINGS_MD5 — PartialReprog vs Fibroblast (did responding cells get younger than where they started?)

**Status:** Task 1: GM00731 scheme=a d10: `pr_below_fibroblast`; GM23815 scheme=a d3: `pr_below_fibroblast`; d7: `pr_below_fibroblast`; d10: `pr_below_fibroblast`. Task 2: report-only QC / extrap / NonReprog − Fibroblast tabulated; no gate. Task 3: ledger of standing/withdrawn readings. Fired: GM00731 `pr_below_fibroblast` scheme=a; GM23815 `pr_below_fibroblast` scheme=a. Seed `20260914`. boot `20260918`. n_perm=200. n_boot=200. n_random=200. Frozen ruler `frozen_ruler_ridge_raw.npz` exists=True.

Does not modify any existing `FINDINGS_*.md` or `FALSIFICATION.md`. No GSE325735. Uncalibrated R² is reported and is never a gate. Nothing averaged across donors, states, timepoints, or instruments. Refit nothing. d7 is not the primary Stage 2 endpoint. Task 2 is report-only and is not used to discount Task 1. Scheme (b) is the weaker design and is not primary where scheme (a) qualified.

Reproduced by `src/md5_run.py`. Frozen ruler: `<repo>\results\fibro\frozen_ruler_ridge_raw.npz` exists=True.

Flag: `PREREG_TASK1.flag` exists=True. `PREREG_TASK2.flag` exists=True. `PREREG_TASK3.flag` exists=True. `DECLARED_BEFORE_SCORES.flag` exists=True.

## Why this contrast is the decisive one

FINDINGS_MD4.md Task 1 found PartialReprog lower than NonReprog on both the frozen ruler and MD, at every qualifying timepoint, in both donors (p at the permutation floor). Its own context rows raise the alternative named in the prompt: on young GM23815 at d10, Fibroblast is also lower than NonReprog (ruler −0.294, MD −0.043). If starting fibroblasts already read younger than NonReprog, then PartialReprog < NonReprog may reflect NonReprog reading old rather than PartialReprog reading young.

Lu et al., Cell 2025, 188:5895–5911:

> another is to the non-reprogrammed state that retains many of the starting fibroblast marker genes (Figure S6B), but it also expresses distinct markers perhaps due to viral infection and/or different culture conditions (Figure 3E)

The contrast that distinguishes the two readings is PartialReprog vs Fibroblast. On aged GM00731 that pair never reached md4's n≥50 threshold at any timepoint. This file runs it at n≥30 (scheme a) and, as the weaker design, pooled Fibroblast d0+d3 vs PartialReprog d3+d7 with timepoint-stratified permutation (scheme b).

Paper quote (Lu et al., Figure 3G):

> The partially reprogrammed populations showed a robust reversal of transcriptomic aging changes, while non-reprogrammed populations maintained the aging transcriptome (Figure 3G). This pattern was accompanied by downregulation of the MD and TGF-β pathway scores in the partially reprogramming cells and their maintained expression in non-reprogrammed cells (Figure 3G).

## Failures recorded

- **unusable cell_cycle:** md2 louvain_obs / louvain_ams have no phase / S.Score / G2M.Score columns. Cell-cycle phase fractions are not computable from the md2 object. Not substituting Seurat CellCycleScoring or a new Tirosh score.

## Task 1 pre-registration (verbatim, written before any MD5 contrast p-value)

Written **before any MD5 contrast p-value**, 2026-09-19. No threshold in this block is re-tuned after numbers exist.

Using the src/md2_*.py Louvain clustering already computed (not re-run, not re-tuned). Cell-state labels are the md2 assignment: argmax of mean AddModuleScore of mmc3 Reprog_cell_state_signatures columns Fibroblast / PartialReprog / EarlyPluripotency / Pluripotency / NonReprog. Labels are read from results/md2/t2_cluster_labels_*.csv. Not re-labelled. AddModuleScore for MD and TGF-β is the per-cell vector already stored in data/processed/md2/louvain_ams_*.npz (the md2 computation). Table S3 Age up / Age down and the frozen-ruler / pluripotency−fibroblast per-cell scores are the vectors already stored in data/processed/md4/cell_scores_*.npz (the md4 computation). Frozen ruler: results/fibro/frozen_ruler_ridge_raw.npz used as-is. Refit nothing. Do not re-cluster, re-label, re-tune, or rescore.

Report aged donor GM00731 and young donor GM23815 separately, never pooled. First report cell counts for Fibroblast and PartialReprog at every timepoint per donor. Then run the contrast under two pre-declared schemes, both reported, neither substituted for the other:

- (a) relaxed threshold: both states ≥30 cells at the same timepoint, no pooling across timepoints. Primary if it qualifies.
- (b) pooled timepoints within state: Fibroblast cells from d0+d3 against PartialReprog cells from d3+d7, each pooled within donor. This breaks md4's within-timepoint permutation guard, so the permutation must shuffle the state label within timepoint strata and the timepoint composition of both arms must be reported. Label this clearly as the weaker design; use it only where (a) does not qualify. Do not report (b) as primary where (a) qualified. Do not lower the threshold below 30. Do not pool donors.

Instruments, all applied unchanged from md4 cell_scores / md2 louvain_ams: frozen ruler, MD, TGF-β, Table S3 age-up minus age-down, pluripotency−fibroblast. Per contrast report Δmean and Δmedian, cell-level permutation p (200 shuffles, stratified as above for scheme (b); unstratified two-state shuffle among the cells in those two states at that donor × timepoint for scheme (a), which is md4_task1._one_contrast unchanged), bootstrap 95% CI (B=200), and Cohen's d (pooled sd, n−1, same sign as Δmean). Seed 20260914 for permutations. Bootstrap seed 20260918. Scheme (a) permutation/bootstrap offsets reuse md4_common.perm_seed / boot_seed for pair PartialReprog vs Fibroblast so already-published md4 context rows must reproduce. Scheme (b) uses a documented offset that does not collide with md4.

Δmean = mean(PartialReprog) − mean(Fibroblast). “Below” / “beating the null” = Δmean < 0 and p_lower ≤ 0.05 on that instrument's cell-level permutation null. Empirical p_lower = (n_null ≤ Δobs + 1) / (n_null + 1).

Sanity, before any new contrast p-value: (1) reproduce FINDINGS_MD2.md pooled MD means on the aged donor (PartialReprog +0.150, NonReprog +0.439; t2_reproduction.json); (2) occupancy matches results/md4/t1_cell_counts.csv. If either fails, stop. Cell-state labels have drifted.

**Pre-registered reading (only the outcome that fired, per donor):**
- PartialReprog below Fibroblast on the frozen ruler, beating its null → responding cells read younger than their own starting population. The rejuvenation reading holds and is not an artifact of NonReprog reading old. State this plainly.
- PartialReprog not below Fibroblast while md4's PartialReprog < NonReprog stands → the md4 result is driven by NonReprog reading old, not by PartialReprog reading young. The rejuvenation reading is not supported by this data, and md4's reading must be restated in those terms. Say so plainly and quote the md4 sentence at issue.
- Instruments disagree (MD below Fibroblast, ruler not, or the reverse) → report both, declare no winner, name what would settle it.
- Neither scheme qualifies on the aged donor → the decisive test cannot be run in this dataset. Say so; do not lower the threshold further, and do not let the young donor stand in for the aged one.

If two qualifying scheme-(a) timepoints in one donor fire different bullets, report both, do not pick, do not average. Scheme (b) is not averaged with scheme (a). Donors are not pooled.


## Task 1 — cell counts by state × timepoint (md2 Louvain labels, not re-labelled)

Assignment basis: argmax of mean AddModuleScore of mmc3 `Reprog_cell_state_signatures` (Fibroblast, PartialReprog, EarlyPluripotency, Pluripotency, NonReprog). Labels from `results/md2/t2_cluster_labels_*.csv`. Not re-labelled. Occupancy matches `results/md4/t1_cell_counts.csv` or this file stops.

Fibroblast and PartialReprog at every timepoint (the pair under test):

| cell_line | day | n_Fibroblast | n_PartialReprog | n_timepoint | both_ge30 |
|---|---|---|---|---|---|
| GM00731 | 0 | 4939 | 3 | 5009 | False |
| GM00731 | 3 | 5 | 5220 | 5988 | False |
| GM00731 | 7 | 14 | 578 | 6795 | False |
| GM00731 | 10 | 1030 | 31 | 4340 | True |
| GM23815 | 0 | 7703 | 9 | 7768 | False |
| GM23815 | 3 | 46 | 6113 | 6576 | True |
| GM23815 | 7 | 2233 | 309 | 11699 | True |
| GM23815 | 10 | 1373 | 50 | 7815 | True |

Pooled n cells by state (timepoints not pooled for scheme a):

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

| cell_line | day | label | n_cells | n_timepoint | frac_timepoint | ge30 | ge50 |
|---|---|---|---|---|---|---|---|
| GM00731 | 0 | Fibroblast | 4939 | 5009 | +0.986 | True | True |
| GM00731 | 0 | PartialReprog | 3 | 5009 | +0.001 | False | False |
| GM00731 | 0 | EarlyPluripotency | 2 | 5009 | +0.000 | False | False |
| GM00731 | 0 | Pluripotency | 0 | 5009 | +0.000 | False | False |
| GM00731 | 0 | NonReprog | 65 | 5009 | +0.013 | True | True |
| GM00731 | 3 | Fibroblast | 5 | 5988 | +0.001 | False | False |
| GM00731 | 3 | PartialReprog | 5220 | 5988 | +0.872 | True | True |
| GM00731 | 3 | EarlyPluripotency | 410 | 5988 | +0.068 | True | True |
| GM00731 | 3 | Pluripotency | 5 | 5988 | +0.001 | False | False |
| GM00731 | 3 | NonReprog | 348 | 5988 | +0.058 | True | True |
| GM00731 | 7 | Fibroblast | 14 | 6795 | +0.002 | False | False |
| GM00731 | 7 | PartialReprog | 578 | 6795 | +0.085 | True | True |
| GM00731 | 7 | EarlyPluripotency | 3917 | 6795 | +0.576 | True | True |
| GM00731 | 7 | Pluripotency | 354 | 6795 | +0.052 | True | True |
| GM00731 | 7 | NonReprog | 1932 | 6795 | +0.284 | True | True |
| GM00731 | 10 | Fibroblast | 1030 | 4340 | +0.237 | True | True |
| GM00731 | 10 | PartialReprog | 31 | 4340 | +0.007 | True | False |
| GM00731 | 10 | EarlyPluripotency | 85 | 4340 | +0.020 | True | True |
| GM00731 | 10 | Pluripotency | 80 | 4340 | +0.018 | True | True |
| GM00731 | 10 | NonReprog | 3114 | 4340 | +0.718 | True | True |
| GM23815 | 0 | Fibroblast | 7703 | 7768 | +0.992 | True | True |
| GM23815 | 0 | PartialReprog | 9 | 7768 | +0.001 | False | False |
| GM23815 | 0 | EarlyPluripotency | 56 | 7768 | +0.007 | True | True |
| GM23815 | 0 | Pluripotency | 0 | 7768 | +0.000 | False | False |
| GM23815 | 0 | NonReprog | 0 | 7768 | +0.000 | False | False |
| GM23815 | 3 | Fibroblast | 46 | 6576 | +0.007 | True | False |
| GM23815 | 3 | PartialReprog | 6113 | 6576 | +0.930 | True | True |
| GM23815 | 3 | EarlyPluripotency | 316 | 6576 | +0.048 | True | True |
| GM23815 | 3 | Pluripotency | 91 | 6576 | +0.014 | True | True |
| GM23815 | 3 | NonReprog | 10 | 6576 | +0.002 | False | False |
| GM23815 | 7 | Fibroblast | 2233 | 11699 | +0.191 | True | True |
| GM23815 | 7 | PartialReprog | 309 | 11699 | +0.026 | True | True |
| GM23815 | 7 | EarlyPluripotency | 8172 | 11699 | +0.699 | True | True |
| GM23815 | 7 | Pluripotency | 977 | 11699 | +0.084 | True | True |
| GM23815 | 7 | NonReprog | 8 | 11699 | +0.001 | False | False |
| GM23815 | 10 | Fibroblast | 1373 | 7815 | +0.176 | True | True |
| GM23815 | 10 | PartialReprog | 50 | 7815 | +0.006 | True | True |
| GM23815 | 10 | EarlyPluripotency | 4417 | 7815 | +0.565 | True | True |
| GM23815 | 10 | Pluripotency | 181 | 7815 | +0.023 | True | True |
| GM23815 | 10 | NonReprog | 1794 | 7815 | +0.230 | True | True |

Scheme (a): a timepoint qualifies iff both Fibroblast and PartialReprog have ≥30 cells. Timepoints are not pooled to reach the threshold. Primary if it qualifies.

| cell_line | day | state_a | state_b | n_a | n_b | min_cells | qualifies | primary | reason |
|---|---|---|---|---|---|---|---|---|---|
| GM00731 | 0 | PartialReprog | Fibroblast | 3 | 4939 | 30 | False | True | n<30 in at least one state |
| GM00731 | 3 | PartialReprog | Fibroblast | 5220 | 5 | 30 | False | True | n<30 in at least one state |
| GM00731 | 7 | PartialReprog | Fibroblast | 578 | 14 | 30 | False | True | n<30 in at least one state |
| GM00731 | 10 | PartialReprog | Fibroblast | 31 | 1030 | 30 | True | True | NA |
| GM23815 | 0 | PartialReprog | Fibroblast | 9 | 7703 | 30 | False | True | n<30 in at least one state |
| GM23815 | 3 | PartialReprog | Fibroblast | 6113 | 46 | 30 | True | True | NA |
| GM23815 | 7 | PartialReprog | Fibroblast | 309 | 2233 | 30 | True | True | NA |
| GM23815 | 10 | PartialReprog | Fibroblast | 50 | 1373 | 30 | True | True | NA |

Scheme (b), weaker design: Fibroblast d0+d3 vs PartialReprog d3+d7, pooled within donor. Permutation shuffles the state label within timepoint strata. Timepoint composition of both arms:

| cell_line | arm | day | n_cells | n_arm | frac_arm |
|---|---|---|---|---|---|
| GM00731 | Fibroblast | 0 | 4939 | 4944 | +0.999 |
| GM00731 | Fibroblast | 3 | 5 | 4944 | +0.001 |
| GM00731 | PartialReprog | 3 | 5220 | 5798 | +0.900 |
| GM00731 | PartialReprog | 7 | 578 | 5798 | +0.100 |
| GM23815 | Fibroblast | 0 | 7703 | 7749 | +0.994 |
| GM23815 | Fibroblast | 3 | 46 | 7749 | +0.006 |
| GM23815 | PartialReprog | 3 | 6113 | 6422 | +0.952 |
| GM23815 | PartialReprog | 7 | 309 | 6422 | +0.048 |

| cell_line | state_a | state_b | n_a | n_b | min_cells | qualifies | primary | scheme_a_qualifies_on_this_donor | n_Fibroblast_d0 | n_Fibroblast_d3 | n_PartialReprog_d3 | n_PartialReprog_d7 | reason | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | PartialReprog | Fibroblast | 5798 | 4944 | 30 | True | False | True | 4939 | 5 | 5220 | 578 | NA | weaker design; permutation shuffles state label within timepoint strata; not primary where scheme (a) qualified |
| GM23815 | PartialReprog | Fibroblast | 6422 | 7749 | 30 | True | False | True | 7703 | 46 | 6113 | 309 | NA | weaker design; permutation shuffles state label within timepoint strata; not primary where scheme (a) qualified |

## Task 1 — sanity check (FINDINGS_MD2.md pooled MD means; occupancy vs md4)

- n_PartialReprog=5832 (md2 json 5832)
- n_NonReprog=5459 (md2 json 5459)
- md_mean_PartialReprog=+0.150471 (md2 json +0.150471; FINDINGS_MD2.md printed +0.150)
- md_mean_NonReprog=+0.439448 (md2 json +0.439448; FINDINGS_MD2.md printed +0.439)
- n_match=True md_match_1e12=True printed_3dp_match=True ok=True
- occupancy vs `results/md4/t1_cell_counts.csv`: match required or STOP

## Task 1 — scheme (a) contrasts (primary where it qualifies)

Δ = mean(PartialReprog) − mean(Fibroblast). Permutation: md4_task1._one_contrast unchanged — 200 shuffles of the two-state label among cells in those two states at that donor × timepoint. Bootstrap: B=200 cells within state, seed `20260918` with md4 pair-offset for PartialReprog vs Fibroblast. Effect size: Cohen's d (pooled sd, n−1), same sign as Δmean. Rows with qualifies=False are not tested. MIN_CELLS=30.

| cell_line | day | state_a | state_b | instrument | qualifies | n_a | n_b | mean_a | mean_b | delta_mean | p_mean_lower | p_mean_two_sided | delta_mean_ci_lo | delta_mean_ci_hi | median_a | median_b | delta_median | p_median_lower | cohens_d | beats_null_lower | ok | reason |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | 0 | PartialReprog | Fibroblast | frozen_ruler | False | 3 | 4939 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM00731 | 0 | PartialReprog | Fibroblast | md_score | False | 3 | 4939 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM00731 | 0 | PartialReprog | Fibroblast | tgfb_score | False | 3 | 4939 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM00731 | 0 | PartialReprog | Fibroblast | age_up_minus_age_down | False | 3 | 4939 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM00731 | 0 | PartialReprog | Fibroblast | pluri_primary | False | 3 | 4939 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM00731 | 3 | PartialReprog | Fibroblast | frozen_ruler | False | 5220 | 5 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM00731 | 3 | PartialReprog | Fibroblast | md_score | False | 5220 | 5 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM00731 | 3 | PartialReprog | Fibroblast | tgfb_score | False | 5220 | 5 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM00731 | 3 | PartialReprog | Fibroblast | age_up_minus_age_down | False | 5220 | 5 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM00731 | 3 | PartialReprog | Fibroblast | pluri_primary | False | 5220 | 5 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM00731 | 7 | PartialReprog | Fibroblast | frozen_ruler | False | 578 | 14 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM00731 | 7 | PartialReprog | Fibroblast | md_score | False | 578 | 14 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM00731 | 7 | PartialReprog | Fibroblast | tgfb_score | False | 578 | 14 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM00731 | 7 | PartialReprog | Fibroblast | age_up_minus_age_down | False | 578 | 14 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM00731 | 7 | PartialReprog | Fibroblast | pluri_primary | False | 578 | 14 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM00731 | 10 | PartialReprog | Fibroblast | frozen_ruler | True | 31 | 1030 | +17.470 | +18.402 | -0.931 | +0.005 | +0.005 | -1.484 | -0.496 | +17.655 | +18.425 | -0.770 | +0.005 | -1.283 | True | True | NA |
| GM00731 | 10 | PartialReprog | Fibroblast | md_score | True | 31 | 1030 | +0.297 | +0.434 | -0.137 | +0.005 | +0.005 | -0.207 | -0.075 | +0.291 | +0.454 | -0.163 | +0.005 | -1.009 | True | True | NA |
| GM00731 | 10 | PartialReprog | Fibroblast | tgfb_score | True | 31 | 1030 | -0.013 | -0.009 | -0.004 | +0.403 | +0.731 | -0.038 | +0.025 | -0.023 | -0.012 | -0.011 | +0.299 | -0.054 | False | True | NA |
| GM00731 | 10 | PartialReprog | Fibroblast | age_up_minus_age_down | True | 31 | 1030 | -0.031 | -0.029 | -0.003 | +0.413 | +0.801 | -0.013 | +0.009 | -0.033 | -0.019 | -0.014 | +0.149 | -0.044 | False | True | NA |
| GM00731 | 10 | PartialReprog | Fibroblast | pluri_primary | True | 31 | 1030 | +10.100 | +8.858 | +1.241 | +1.000 | +0.005 | +0.673 | +1.889 | +9.898 | +8.681 | +1.217 | +1.000 | +1.351 | False | True | NA |
| GM23815 | 0 | PartialReprog | Fibroblast | frozen_ruler | False | 9 | 7703 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM23815 | 0 | PartialReprog | Fibroblast | md_score | False | 9 | 7703 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM23815 | 0 | PartialReprog | Fibroblast | tgfb_score | False | 9 | 7703 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM23815 | 0 | PartialReprog | Fibroblast | age_up_minus_age_down | False | 9 | 7703 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM23815 | 0 | PartialReprog | Fibroblast | pluri_primary | False | 9 | 7703 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM23815 | 3 | PartialReprog | Fibroblast | frozen_ruler | True | 6113 | 46 | +18.184 | +18.849 | -0.664 | +0.005 | +0.005 | -0.857 | -0.463 | +18.176 | +18.929 | -0.753 | +0.005 | -1.135 | True | True | NA |
| GM23815 | 3 | PartialReprog | Fibroblast | md_score | True | 6113 | 46 | +0.085 | +0.385 | -0.300 | +0.005 | +0.005 | -0.326 | -0.263 | +0.070 | +0.400 | -0.331 | +0.005 | -2.216 | True | True | NA |
| GM23815 | 3 | PartialReprog | Fibroblast | tgfb_score | True | 6113 | 46 | -0.036 | +0.081 | -0.118 | +0.005 | +0.005 | -0.143 | -0.090 | -0.043 | +0.088 | -0.131 | +0.005 | -1.188 | True | True | NA |
| GM23815 | 3 | PartialReprog | Fibroblast | age_up_minus_age_down | True | 6113 | 46 | -0.114 | -0.049 | -0.064 | +0.005 | +0.005 | -0.079 | -0.052 | -0.121 | -0.053 | -0.067 | +0.005 | -1.304 | True | True | NA |
| GM23815 | 3 | PartialReprog | Fibroblast | pluri_primary | True | 6113 | 46 | +12.175 | +10.345 | +1.830 | +1.000 | +0.005 | +1.475 | +2.135 | +12.118 | +9.938 | +2.180 | +1.000 | +1.394 | False | True | NA |
| GM23815 | 7 | PartialReprog | Fibroblast | frozen_ruler | True | 309 | 2233 | +18.544 | +18.653 | -0.110 | +0.005 | +0.010 | -0.187 | -0.020 | +18.559 | +18.655 | -0.095 | +0.005 | -0.177 | True | True | NA |
| GM23815 | 7 | PartialReprog | Fibroblast | md_score | True | 309 | 2233 | +0.174 | +0.410 | -0.236 | +0.005 | +0.005 | -0.255 | -0.218 | +0.193 | +0.409 | -0.216 | +0.005 | -2.421 | True | True | NA |
| GM23815 | 7 | PartialReprog | Fibroblast | tgfb_score | True | 309 | 2233 | -0.058 | -0.015 | -0.044 | +0.005 | +0.005 | -0.053 | -0.031 | -0.065 | -0.019 | -0.046 | +0.005 | -0.527 | True | True | NA |
| GM23815 | 7 | PartialReprog | Fibroblast | age_up_minus_age_down | True | 309 | 2233 | -0.047 | +0.008 | -0.056 | +0.005 | +0.005 | -0.062 | -0.049 | -0.047 | +0.013 | -0.061 | +0.005 | -1.171 | True | True | NA |
| GM23815 | 7 | PartialReprog | Fibroblast | pluri_primary | True | 309 | 2233 | +11.409 | +9.618 | +1.791 | +1.000 | +0.005 | +1.642 | +1.972 | +10.775 | +9.559 | +1.217 | +1.000 | +2.222 | False | True | NA |
| GM23815 | 10 | PartialReprog | Fibroblast | frozen_ruler | True | 50 | 1373 | +18.420 | +18.894 | -0.474 | +0.005 | +0.005 | -0.739 | -0.239 | +18.549 | +18.912 | -0.363 | +0.005 | -0.955 | True | True | NA |
| GM23815 | 10 | PartialReprog | Fibroblast | md_score | True | 50 | 1373 | +0.128 | +0.379 | -0.251 | +0.005 | +0.005 | -0.289 | -0.210 | +0.144 | +0.382 | -0.238 | +0.005 | -2.551 | True | True | NA |
| GM23815 | 10 | PartialReprog | Fibroblast | tgfb_score | True | 50 | 1373 | -0.043 | +0.037 | -0.080 | +0.005 | +0.005 | -0.108 | -0.052 | -0.057 | +0.036 | -0.093 | +0.005 | -1.024 | True | True | NA |
| GM23815 | 10 | PartialReprog | Fibroblast | age_up_minus_age_down | True | 50 | 1373 | -0.069 | -0.036 | -0.033 | +0.005 | +0.005 | -0.054 | -0.011 | -0.068 | -0.025 | -0.043 | +0.005 | -0.583 | True | True | NA |
| GM23815 | 10 | PartialReprog | Fibroblast | pluri_primary | True | 50 | 1373 | +12.019 | +9.680 | +2.338 | +1.000 | +0.005 | +1.855 | +2.837 | +11.133 | +9.652 | +1.482 | +1.000 | +3.868 | False | True | NA |

Qualifying scheme-(a) rows only:

| cell_line | day | instrument | n_a | n_b | mean_a | mean_b | delta_mean | p_mean_lower | p_mean_two_sided | delta_mean_ci_lo | delta_mean_ci_hi | delta_median | p_median_lower | cohens_d | beats_null_lower |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | 10 | frozen_ruler | 31 | 1030 | +17.470 | +18.402 | -0.931 | +0.005 | +0.005 | -1.484 | -0.496 | -0.770 | +0.005 | -1.283 | True |
| GM00731 | 10 | md_score | 31 | 1030 | +0.297 | +0.434 | -0.137 | +0.005 | +0.005 | -0.207 | -0.075 | -0.163 | +0.005 | -1.009 | True |
| GM00731 | 10 | tgfb_score | 31 | 1030 | -0.013 | -0.009 | -0.004 | +0.403 | +0.731 | -0.038 | +0.025 | -0.011 | +0.299 | -0.054 | False |
| GM00731 | 10 | age_up_minus_age_down | 31 | 1030 | -0.031 | -0.029 | -0.003 | +0.413 | +0.801 | -0.013 | +0.009 | -0.014 | +0.149 | -0.044 | False |
| GM00731 | 10 | pluri_primary | 31 | 1030 | +10.100 | +8.858 | +1.241 | +1.000 | +0.005 | +0.673 | +1.889 | +1.217 | +1.000 | +1.351 | False |
| GM23815 | 3 | frozen_ruler | 6113 | 46 | +18.184 | +18.849 | -0.664 | +0.005 | +0.005 | -0.857 | -0.463 | -0.753 | +0.005 | -1.135 | True |
| GM23815 | 3 | md_score | 6113 | 46 | +0.085 | +0.385 | -0.300 | +0.005 | +0.005 | -0.326 | -0.263 | -0.331 | +0.005 | -2.216 | True |
| GM23815 | 3 | tgfb_score | 6113 | 46 | -0.036 | +0.081 | -0.118 | +0.005 | +0.005 | -0.143 | -0.090 | -0.131 | +0.005 | -1.188 | True |
| GM23815 | 3 | age_up_minus_age_down | 6113 | 46 | -0.114 | -0.049 | -0.064 | +0.005 | +0.005 | -0.079 | -0.052 | -0.067 | +0.005 | -1.304 | True |
| GM23815 | 3 | pluri_primary | 6113 | 46 | +12.175 | +10.345 | +1.830 | +1.000 | +0.005 | +1.475 | +2.135 | +2.180 | +1.000 | +1.394 | False |
| GM23815 | 7 | frozen_ruler | 309 | 2233 | +18.544 | +18.653 | -0.110 | +0.005 | +0.010 | -0.187 | -0.020 | -0.095 | +0.005 | -0.177 | True |
| GM23815 | 7 | md_score | 309 | 2233 | +0.174 | +0.410 | -0.236 | +0.005 | +0.005 | -0.255 | -0.218 | -0.216 | +0.005 | -2.421 | True |
| GM23815 | 7 | tgfb_score | 309 | 2233 | -0.058 | -0.015 | -0.044 | +0.005 | +0.005 | -0.053 | -0.031 | -0.046 | +0.005 | -0.527 | True |
| GM23815 | 7 | age_up_minus_age_down | 309 | 2233 | -0.047 | +0.008 | -0.056 | +0.005 | +0.005 | -0.062 | -0.049 | -0.061 | +0.005 | -1.171 | True |
| GM23815 | 7 | pluri_primary | 309 | 2233 | +11.409 | +9.618 | +1.791 | +1.000 | +0.005 | +1.642 | +1.972 | +1.217 | +1.000 | +2.222 | False |
| GM23815 | 10 | frozen_ruler | 50 | 1373 | +18.420 | +18.894 | -0.474 | +0.005 | +0.005 | -0.739 | -0.239 | -0.363 | +0.005 | -0.955 | True |
| GM23815 | 10 | md_score | 50 | 1373 | +0.128 | +0.379 | -0.251 | +0.005 | +0.005 | -0.289 | -0.210 | -0.238 | +0.005 | -2.551 | True |
| GM23815 | 10 | tgfb_score | 50 | 1373 | -0.043 | +0.037 | -0.080 | +0.005 | +0.005 | -0.108 | -0.052 | -0.093 | +0.005 | -1.024 | True |
| GM23815 | 10 | age_up_minus_age_down | 50 | 1373 | -0.069 | -0.036 | -0.033 | +0.005 | +0.005 | -0.054 | -0.011 | -0.043 | +0.005 | -0.583 | True |
| GM23815 | 10 | pluri_primary | 50 | 1373 | +12.019 | +9.680 | +2.338 | +1.000 | +0.005 | +1.855 | +2.837 | +1.482 | +1.000 | +3.868 | False |

## Task 1 — scheme (b) contrasts (weaker design; not primary where scheme (a) qualified)

Fibroblast cells from days [0, 3] vs PartialReprog cells from days [3, 7], pooled within donor. Permutation shuffles the state label **within timepoint strata** (200 shuffles). Bootstrap resamples cells within each state independently (B=200), as in md4. This design lets timepoint composition of the two arms persist in every permutation. It is the weaker design. It is not the primary reading on any donor where scheme (a) qualified.

| cell_line | instrument | qualifies | primary | n_a | n_b | n_Fibroblast_d0 | n_Fibroblast_d3 | n_PartialReprog_d3 | n_PartialReprog_d7 | mean_a | mean_b | delta_mean | p_mean_lower | p_mean_two_sided | delta_mean_ci_lo | delta_mean_ci_hi | delta_median | p_median_lower | cohens_d | beats_null_lower | n_strata_mixed | ok |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | frozen_ruler | True | False | 5798 | 4944 | 4939 | 5 | 5220 | 578 | +17.432 | +18.157 | -0.724 | +0.035 | +0.035 | -0.754 | -0.696 | -0.702 | +0.139 | -1.004 | True | 1 | True |
| GM00731 | md_score | True | False | 5798 | 4944 | 4939 | 5 | 5220 | 578 | +0.150 | +0.506 | -0.356 | +0.005 | +0.005 | -0.362 | -0.352 | -0.383 | +0.045 | -2.715 | True | 1 | True |
| GM00731 | tgfb_score | True | False | 5798 | 4944 | 4939 | 5 | 5220 | 578 | -0.044 | +0.042 | -0.086 | +0.075 | +0.075 | -0.089 | -0.083 | -0.089 | +0.030 | -0.948 | False | 1 | True |
| GM00731 | age_up_minus_age_down | True | False | 5798 | 4944 | 4939 | 5 | 5220 | 578 | -0.096 | +0.029 | -0.125 | +0.020 | +0.020 | -0.127 | -0.123 | -0.135 | +0.194 | -2.352 | True | 1 | True |
| GM00731 | pluri_primary | True | False | 5798 | 4944 | 4939 | 5 | 5220 | 578 | +12.075 | +8.291 | +3.784 | +0.990 | +0.015 | +3.734 | +3.830 | +3.920 | +0.975 | +2.792 | False | 1 | True |
| GM23815 | frozen_ruler | True | False | 6422 | 7749 | 7703 | 46 | 6113 | 309 | +18.202 | +18.879 | -0.677 | +0.005 | +0.005 | -0.698 | -0.656 | -0.699 | +0.005 | -1.200 | True | 1 | True |
| GM23815 | md_score | True | False | 6422 | 7749 | 7703 | 46 | 6113 | 309 | +0.090 | +0.478 | -0.389 | +0.005 | +0.005 | -0.393 | -0.385 | -0.406 | +0.005 | -3.587 | True | 1 | True |
| GM23815 | tgfb_score | True | False | 6422 | 7749 | 7703 | 46 | 6113 | 309 | -0.037 | +0.029 | -0.067 | +0.005 | +0.005 | -0.070 | -0.064 | -0.075 | +0.005 | -0.752 | True | 1 | True |
| GM23815 | age_up_minus_age_down | True | False | 6422 | 7749 | 7703 | 46 | 6113 | 309 | -0.110 | +0.009 | -0.120 | +0.005 | +0.005 | -0.121 | -0.118 | -0.132 | +0.005 | -2.197 | True | 1 | True |
| GM23815 | pluri_primary | True | False | 6422 | 7749 | 7703 | 46 | 6113 | 309 | +12.138 | +9.406 | +2.732 | +1.000 | +0.005 | +2.698 | +2.765 | +2.701 | +1.000 | +2.854 | False | 1 | True |

## Task 1 reading

Scheme (a) is primary where it qualifies. Scheme (b) is tabulated above and is used for the fired bullet only where scheme (a) does not qualify. Donors are not pooled. Timepoints are not averaged. No winner is declared between instruments.

### GM00731

Fired key: `pr_below_fibroblast`. scheme_used=`a`. scheme_a_is_primary=True. n_qualifying_a=1 n_qualifying_b=1.

scheme=a d10: `pr_below_fibroblast`.

PartialReprog below Fibroblast on the frozen ruler, beating its null. Responding cells read younger than their own starting population. The rejuvenation reading holds and is not an artifact of NonReprog reading old. GM00731 scheme=a d10 PartialReprog vs Fibroblast n_PR=31 n_Fib=1030 MD Δmean=-0.13701588483197658 p_lower=0.004975124378109453 Cohen_d=-1.0091524212841425 ruler Δmean=-0.931403319338088 p_lower=0.004975124378109453 Cohen_d=-1.2831818356450844.

Scheme (b) also qualifies on this donor and is tabulated; it is not the primary reading.

### GM23815

Fired key: `pr_below_fibroblast`. scheme_used=`a`. scheme_a_is_primary=True. n_qualifying_a=3 n_qualifying_b=1.

scheme=a d3: `pr_below_fibroblast`.

PartialReprog below Fibroblast on the frozen ruler, beating its null. Responding cells read younger than their own starting population. The rejuvenation reading holds and is not an artifact of NonReprog reading old. GM23815 scheme=a d3 PartialReprog vs Fibroblast n_PR=6113 n_Fib=46 MD Δmean=-0.29987847032668635 p_lower=0.004975124378109453 Cohen_d=-2.2158509134958058 ruler Δmean=-0.6643573111844852 p_lower=0.004975124378109453 Cohen_d=-1.1345691552224035.

scheme=a d7: `pr_below_fibroblast`.

PartialReprog below Fibroblast on the frozen ruler, beating its null. Responding cells read younger than their own starting population. The rejuvenation reading holds and is not an artifact of NonReprog reading old. GM23815 scheme=a d7 PartialReprog vs Fibroblast n_PR=309 n_Fib=2233 MD Δmean=-0.2357721618356915 p_lower=0.004975124378109453 Cohen_d=-2.420975098758663 ruler Δmean=-0.10963599174700889 p_lower=0.004975124378109453 Cohen_d=-0.17682718606618064.

scheme=a d10: `pr_below_fibroblast`.

PartialReprog below Fibroblast on the frozen ruler, beating its null. Responding cells read younger than their own starting population. The rejuvenation reading holds and is not an artifact of NonReprog reading old. GM23815 scheme=a d10 PartialReprog vs Fibroblast n_PR=50 n_Fib=1373 MD Δmean=-0.25128324204493047 p_lower=0.004975124378109453 Cohen_d=-2.5514397494328436 ruler Δmean=-0.4743939721469417 p_lower=0.004975124378109453 Cohen_d=-0.9545810171291685.

Scheme (b) also qualifies on this donor and is tabulated; it is not the primary reading.

Δmean < 0 with p_lower ≤ 0.05 on the frozen ruler is “PartialReprog below Fibroblast / beating the null”. Donors are not pooled. Scheme (a) timepoints are not pooled.

## Task 2 pre-registration (verbatim, written before any MD5 QC / extrap / NR − Fib p-value)

Written **before any MD5 QC median, cell-cycle fraction, extrapolation row, or NonReprog − Fibroblast contrast p-value**, 2026-09-19. Report-only. No gate.

Per donor × timepoint, for each of NonReprog, Fibroblast and PartialReprog with ≥30 Louvain-kept cells: median UMI, median genes, mitochondrial fraction, cell-cycle phase fractions if computable from the md2 object (louvain_obs / louvain_ams / cluster labels). If those objects have no phase / S.Score / G2M.Score columns, record that cell-cycle is not computable from the md2 object and do not substitute Seurat CellCycleScoring or a new Tirosh score.

Extrapolation distance to the GTEx fibroblast training distribution: the same PCA k=50 nearest-neighbour Euclidean and Mahalanobis measures as FINDINGS_FIBRO2.md / src/fibro2_extrap.py, already computed in results/md4/t2_extrap.csv (TMM among that donor's state×timepoint rows, then frozen μ/σ/w). Read those rows. Do not re-TMM, do not refit PCA, do not rescore.

Also report, for each instrument, the NonReprog − Fibroblast contrast at every timepoint where both states have ≥30 cells, calling md4_task1._one_contrast unchanged (unstratified within timepoint). If NonReprog reads consistently older than the starting fibroblasts across donors and timepoints, state that as an observation and note the paper's own remark about viral infection and culture conditions in those cells. Do not use Task 2 to discount Task 1.


## Task 2 — QC (report-only, no gate)

Per donor × timepoint × {NonReprog, Fibroblast, PartialReprog} with n≥30: median UMI, median genes, mitochondrial fraction from `results/md2/louvain_obs_*.csv`. Cell-cycle phase fractions if computable from the md2 object.

| cell_line | day | label | n_cells | ok | median_umi | median_genes | median_mito_frac | cell_cycle_computable | reason |
|---|---|---|---|---|---|---|---|---|---|
| GM00731 | 0 | NonReprog | 65 | True | +3904.000 | +1659.000 | +0.051 | False | NA |
| GM00731 | 0 | Fibroblast | 4939 | True | +14469.000 | +4094.000 | +0.025 | False | NA |
| GM00731 | 0 | PartialReprog | 3 | False | NA | NA | NA | False | n<30 |
| GM00731 | 3 | NonReprog | 348 | True | +4735.000 | +2138.000 | +0.025 | False | NA |
| GM00731 | 3 | Fibroblast | 5 | False | NA | NA | NA | False | n<30 |
| GM00731 | 3 | PartialReprog | 5220 | True | +11444.500 | +3922.000 | +0.037 | False | NA |
| GM00731 | 7 | NonReprog | 1932 | True | +8500.000 | +3058.000 | +0.038 | False | NA |
| GM00731 | 7 | Fibroblast | 14 | False | NA | NA | NA | False | n<30 |
| GM00731 | 7 | PartialReprog | 578 | True | +11502.000 | +4018.500 | +0.036 | False | NA |
| GM00731 | 10 | NonReprog | 3114 | True | +12389.500 | +3977.000 | +0.036 | False | NA |
| GM00731 | 10 | Fibroblast | 1030 | True | +20629.500 | +5162.000 | +0.025 | False | NA |
| GM00731 | 10 | PartialReprog | 31 | True | +1979.000 | +1205.000 | +0.017 | False | NA |
| GM23815 | 0 | NonReprog | 0 | False | NA | NA | NA | False | n<30 |
| GM23815 | 0 | Fibroblast | 7703 | True | +8169.000 | +3038.000 | +0.029 | False | NA |
| GM23815 | 0 | PartialReprog | 9 | False | NA | NA | NA | False | n<30 |
| GM23815 | 3 | NonReprog | 10 | False | NA | NA | NA | False | n<30 |
| GM23815 | 3 | Fibroblast | 46 | True | +6851.500 | +2799.000 | +0.042 | False | NA |
| GM23815 | 3 | PartialReprog | 6113 | True | +11265.000 | +4014.000 | +0.048 | False | NA |
| GM23815 | 7 | NonReprog | 8 | False | NA | NA | NA | False | n<30 |
| GM23815 | 7 | Fibroblast | 2233 | True | +4637.000 | +2041.000 | +0.040 | False | NA |
| GM23815 | 7 | PartialReprog | 309 | True | +7489.000 | +3001.000 | +0.046 | False | NA |
| GM23815 | 10 | NonReprog | 1794 | True | +4912.000 | +2206.000 | +0.037 | False | NA |
| GM23815 | 10 | Fibroblast | 1373 | True | +8846.000 | +3140.000 | +0.042 | False | NA |
| GM23815 | 10 | PartialReprog | 50 | True | +9196.000 | +3437.000 | +0.047 | False | NA |

Cell-cycle from the md2 object:

- GM00731: computable=False obs_hits=[] ams_hits=[]
  - md2 louvain_obs / louvain_ams have no phase / S.Score / G2M.Score columns. Cell-cycle phase fractions are not computable from the md2 object. Not substituting Seurat CellCycleScoring or a new Tirosh score.
- GM23815: computable=False obs_hits=[] ams_hits=[]
  - md2 louvain_obs / louvain_ams have no phase / S.Score / G2M.Score columns. Cell-cycle phase fractions are not computable from the md2 object. Not substituting Seurat CellCycleScoring or a new Tirosh score.

## Task 2 — extrapolation (report-only, no gate; read from md4, not recomputed)

PCA k=50. Nearest neighbour is Euclidean in GTEx fibroblast training z-space. Mahalanobis is in the leading PCA of that z-space. Rows read from `results/md4/t2_extrap.csv` (TMM among that donor's state×timepoint rows, then frozen μ/σ/w). Same measures as FINDINGS_FIBRO2.md / `src/fibro2_extrap.py`. Not refit. Not re-TMM'd.

| cell_line | day | label | n_cells | age_score | nn_euclidean | mahalanobis_pca | pca_k |
|---|---|---|---|---|---|---|---|
| GM00731 | 0 | NonReprog | 65 | -0.077 | +670.751 | +26.971 | 50 |
| GM00731 | 0 | Fibroblast | 4939 | +5.253 | +465.617 | +24.411 | 50 |
| GM00731 | 3 | NonReprog | 348 | +8.122 | +476.754 | +22.042 | 50 |
| GM00731 | 3 | PartialReprog | 5220 | +2.337 | +479.017 | +23.982 | 50 |
| GM00731 | 7 | NonReprog | 1932 | +7.015 | +470.259 | +24.100 | 50 |
| GM00731 | 7 | PartialReprog | 578 | +4.737 | +468.982 | +23.158 | 50 |
| GM00731 | 10 | NonReprog | 3114 | +9.474 | +485.365 | +24.012 | 50 |
| GM00731 | 10 | Fibroblast | 1030 | +7.873 | +474.154 | +24.294 | 50 |
| GM00731 | 10 | PartialReprog | 31 | +4.933 | +709.005 | +31.008 | 50 |
| GM23815 | 0 | Fibroblast | 7703 | +0.241 | +449.237 | +23.054 | 50 |
| GM23815 | 3 | Fibroblast | 46 | +5.947 | +593.699 | +25.444 | 50 |
| GM23815 | 3 | PartialReprog | 6113 | -2.952 | +455.379 | +23.014 | 50 |
| GM23815 | 7 | Fibroblast | 2233 | -1.817 | +465.875 | +23.737 | 50 |
| GM23815 | 7 | PartialReprog | 309 | -2.809 | +482.403 | +24.240 | 50 |
| GM23815 | 10 | NonReprog | 1794 | +2.866 | +480.957 | +24.440 | 50 |
| GM23815 | 10 | Fibroblast | 1373 | +0.181 | +468.409 | +24.716 | 50 |
| GM23815 | 10 | PartialReprog | 50 | -1.208 | +566.537 | +26.410 | 50 |

## Task 2 — NonReprog − Fibroblast contrasts (report-only, no gate)

Δ = mean(NonReprog) − mean(Fibroblast). Same nulls and CIs as Task 1 scheme (a): md4_task1._one_contrast unchanged, within timepoint, n≥30 both states. `nr_older` = Δmean > 0 and p_two_sided ≤ 0.05.

| cell_line | day | instrument | qualifies | n_a | n_b | mean_a | mean_b | delta_mean | p_mean_lower | p_mean_two_sided | delta_mean_ci_lo | delta_mean_ci_hi | delta_median | cohens_d | beats_null_lower | nr_older | ok | reason |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | 0 | frozen_ruler | True | 65 | 4939 | +16.398 | +18.157 | -1.758 | +0.005 | +0.005 | -2.057 | -1.412 | -2.114 | -2.756 | True | False | True | NA |
| GM00731 | 0 | md_score | True | 65 | 4939 | +0.350 | +0.506 | -0.156 | +0.005 | +0.005 | -0.189 | -0.119 | -0.189 | -2.218 | True | False | True | NA |
| GM00731 | 0 | tgfb_score | True | 65 | 4939 | +0.086 | +0.042 | +0.044 | +1.000 | +0.005 | +0.013 | +0.078 | +0.062 | +0.642 | False | True | True | NA |
| GM00731 | 0 | age_up_minus_age_down | True | 65 | 4939 | +0.137 | +0.029 | +0.108 | +1.000 | +0.005 | +0.088 | +0.128 | +0.126 | +1.876 | False | True | True | NA |
| GM00731 | 0 | pluri_primary | True | 65 | 4939 | +10.146 | +8.289 | +1.858 | +1.000 | +0.005 | +1.591 | +2.139 | +2.217 | +4.032 | False | True | True | NA |
| GM00731 | 3 | frozen_ruler | False | 348 | 5 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM00731 | 3 | md_score | False | 348 | 5 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM00731 | 3 | tgfb_score | False | 348 | 5 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM00731 | 3 | age_up_minus_age_down | False | 348 | 5 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM00731 | 3 | pluri_primary | False | 348 | 5 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM00731 | 7 | frozen_ruler | False | 1932 | 14 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM00731 | 7 | md_score | False | 1932 | 14 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM00731 | 7 | tgfb_score | False | 1932 | 14 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM00731 | 7 | age_up_minus_age_down | False | 1932 | 14 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM00731 | 7 | pluri_primary | False | 1932 | 14 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM00731 | 10 | frozen_ruler | True | 3114 | 1030 | +18.433 | +18.402 | +0.031 | +0.861 | +0.279 | -0.022 | +0.086 | +0.121 | +0.036 | False | False | True | NA |
| GM00731 | 10 | md_score | True | 3114 | 1030 | +0.434 | +0.434 | -0.000 | +0.463 | +0.930 | -0.009 | +0.009 | -0.022 | -0.004 | False | False | True | NA |
| GM00731 | 10 | tgfb_score | True | 3114 | 1030 | +0.001 | -0.009 | +0.010 | +0.995 | +0.010 | +0.004 | +0.015 | +0.005 | +0.117 | False | True | True | NA |
| GM00731 | 10 | age_up_minus_age_down | True | 3114 | 1030 | +0.032 | -0.029 | +0.060 | +1.000 | +0.005 | +0.056 | +0.065 | +0.054 | +1.071 | False | True | True | NA |
| GM00731 | 10 | pluri_primary | True | 3114 | 1030 | +8.821 | +8.858 | -0.037 | +0.119 | +0.209 | -0.092 | +0.019 | +0.024 | -0.050 | False | False | True | NA |
| GM23815 | 0 | frozen_ruler | False | 0 | 7703 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM23815 | 0 | md_score | False | 0 | 7703 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM23815 | 0 | tgfb_score | False | 0 | 7703 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM23815 | 0 | age_up_minus_age_down | False | 0 | 7703 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM23815 | 0 | pluri_primary | False | 0 | 7703 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM23815 | 3 | frozen_ruler | False | 10 | 46 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM23815 | 3 | md_score | False | 10 | 46 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM23815 | 3 | tgfb_score | False | 10 | 46 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM23815 | 3 | age_up_minus_age_down | False | 10 | 46 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM23815 | 3 | pluri_primary | False | 10 | 46 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM23815 | 7 | frozen_ruler | False | 8 | 2233 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM23815 | 7 | md_score | False | 8 | 2233 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM23815 | 7 | tgfb_score | False | 8 | 2233 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM23815 | 7 | age_up_minus_age_down | False | 8 | 2233 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM23815 | 7 | pluri_primary | False | 8 | 2233 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | False | n<30 in at least one state |
| GM23815 | 10 | frozen_ruler | True | 1794 | 1373 | +19.189 | +18.894 | +0.294 | +1.000 | +0.005 | +0.261 | +0.322 | +0.280 | +0.567 | False | True | True | NA |
| GM23815 | 10 | md_score | True | 1794 | 1373 | +0.423 | +0.379 | +0.043 | +1.000 | +0.005 | +0.037 | +0.049 | +0.040 | +0.503 | False | True | True | NA |
| GM23815 | 10 | tgfb_score | True | 1794 | 1373 | +0.035 | +0.037 | -0.001 | +0.343 | +0.667 | -0.007 | +0.004 | -0.001 | -0.017 | False | False | True | NA |
| GM23815 | 10 | age_up_minus_age_down | True | 1794 | 1373 | +0.020 | -0.036 | +0.056 | +1.000 | +0.005 | +0.052 | +0.059 | +0.049 | +1.150 | False | True | True | NA |
| GM23815 | 10 | pluri_primary | True | 1794 | 1373 | +9.368 | +9.680 | -0.312 | +0.005 | +0.005 | -0.343 | -0.281 | -0.311 | -0.674 | True | False | True | NA |

| instrument | n_qualifying | n_nr_older | consistently_older |
|---|---|---|---|
| frozen_ruler | 3 | 1 | False |
| md_score | 3 | 1 | False |
| tgfb_score | 3 | 2 | False |
| age_up_minus_age_down | 3 | 3 | True |
| pluri_primary | 3 | 1 | False |

NonReprog does not read consistently older than Fibroblast on the frozen ruler across all qualifying donor × timepoints. Observation only. Not used to discount Task 1. Paper remark (not used as a gate): "another is to the non-reprogrammed state that retains many of the starting fibroblast marker genes (Figure S6B), but it also expresses distinct markers perhaps due to viral infection and/or different culture conditions (Figure 3E)".

- report_only=True used_to_discount_task1=False

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

## What this changes in FINDINGS_MD4.md (quote, do not edit)

This file does not edit `FINDINGS_MD4.md`, `FINDINGS_MD3.md`, `FINDINGS_MD2.md`, `FINDINGS_FIBRO.md`, `FINDINGS_FIBRO2.md`, `FINDINGS_FIBRO3.md`, or `FALSIFICATION.md`.

FINDINGS_MD4.md Task 1 reading, quoted:

> Both instruments agree that partially reprogrammed cells read younger. Their claim is supported and our earlier negatives reflected the wrong test (a within-state trajectory with no baseline).

Qualification from this file, aged donor, primary scheme: PartialReprog below Fibroblast on the frozen ruler, beating its null. The md4 PartialReprog < NonReprog reading is not an artifact of NonReprog reading old on that donor. The quoted md4 sentence stands with that qualification.

Young donor fired key=`pr_below_fibroblast` scheme=`a`. Not pooled with the aged donor. Not a substitute for the aged-donor test.

## Limitations

1. Two donors in GSE297234.
2. Clusters are not lineage-tracked. PartialReprog / Fibroblast / NonReprog labels are argmax of mmc3 Reprog_cell_state_signatures AddModuleScore, not Slingshot trajectories from day-0 fibroblasts. States are inferred from expression rather than lineage.
3. A state contrast is not a trajectory. This file does not revive the MD3 d0→d7 / d0→d10 within-PartialReprog test.
4. Scheme (b) mixes timepoints. Permutation is stratified within timepoint, so timepoint composition of the two arms persists in every shuffle. Scheme (b) is the weaker design.
5. Louvain is a Python deviation from their R/sctransform pipeline, as recorded in FINDINGS_MD2.md (LogNormalize + quadratic HVG + percent.mt residualization + PCA 50 + SNN + networkx louvain, resolution=0.8 Seurat default; STAR Methods omit resolution).
6. Our reproduction of their metrics is ours, not theirs. AddModuleScore is a Python reimplementation of the published algorithm, not Seurat's C++/R object.
7. Frozen ruler trained on GTEx V10 cultured fibroblasts, public `AGE` 10-year bins. Per-cell scores are the md4 log2-CPM (TMM nf=1) vectors; not refit.
8. Cross-platform shift from GTEx bulk polyA (RNASeQCv2.4.2) to 10x 3' scRNA-seq.
9. Missing overlap genes at z=0. Nothing fitted on GSE297234.
10. Scheme (a) MIN_CELLS=30. The threshold is not lowered further. Donors are not pooled.
11. Seed `20260914`. Bootstrap seed `20260918`. n_perm=200, n_boot=200, n_random=200.
12. GSE325735 not opened. d7 is not the Stage 2 endpoint.
13. Task 2 distances are on TMM state×timepoint pseudobulks from md4; Task 1 ruler scores are per-cell log2-CPM. Those two ruler numbers are not interchangeable.
14. Scheme (a) on the aged donor, if it qualifies, is at the timepoint(s) that reach n≥30 in both states. Fibroblast cells at that timepoint are expression-labelled fibroblasts, not lineage-tracked day-0 cells.
15. Task 2 QC and extrapolation are not used to explain away a Task 1 result.

## Open questions

1. Would the authors' sctransform v2 object change which Louvain clusters are PartialReprog vs Fibroblast at these timepoints?
2. If cells were lineage-tracked from day-0 fibroblasts, would the same contrast hold?
3. An independent fibroblast age instrument on these same labelled cells would settle an instrument disagreement if one remains.

## Gene-list versions

- MD n=205 TGFB n=54 Age up n=1533 Age down n=2007
- mmc3 equals md2 copies: up=True down=True
- used_as_deseq2_rebuild=False used_as_published_lists=True
- AddModuleScore nbin=24 ctrl=100

- GM00731 Age up mapped 896/1533 missing=636; Age down mapped 1586/2007 missing=421 (from md4 `t1_s3_ams_meta_GM00731.json`)
- GM23815 Age up mapped 894/1533 missing=638; Age down mapped 1582/2007 missing=425 (from md4 `t1_s3_ams_meta_GM23815.json`)

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
- **md4_t1_contrasts** source=`<repo>\results\md4\t1_contrasts.csv` columns=['cell_line', 'day', 'state_a', 'state_b', 'n_a', 'n_b', 'qualifies', 'primary', 'min_cells', 'instrument', 'ok', 'reason', 'n_perm', 'n_boot', 'perm_seed', 'boot_seed', 'effect_size', 'mean_a', 'mean_b', 'delta_mean', 'median_a', 'median_b', 'delta_median', 'cohens_d', 'p_mean_lower', 'p_mean_two_sided', 'n_perm_mean_le_obs', 'p_median_lower', 'p_median_two_sided', 'null_mean_median', 'beats_null_lower', 'delta_mean_ci_lo', 'delta_mean_ci_hi', 'delta_median_ci_lo', 'delta_median_ci_hi']
- **md4_t2_extrap** source=`<repo>\results\md4\t2_extrap.csv` columns=['cell_line', 'day', 'label', 'n_cells', 'age_score', 'nn_euclidean', 'nn_train_index', 'mahalanobis_pca', 'pca_k', 'age_source', 'n_missing_z0']
- **t2_louvain_ams_GM00731** source=`<repo>\data\processed\md2\louvain_ams_GM00731.npz` columns=['md', 'tgfb', 'bins', 'gene_mean', 'cluster', 'day', 'n_genes']
- **t2_louvain_ams_GM23815** source=`<repo>\data\processed\md2\louvain_ams_GM23815.npz` columns=['md', 'tgfb', 'bins', 'gene_mean', 'cluster', 'day']

## Gene-set names actually read

- **MD_built_from_md2** name=`MD (STAR Methods built; md2)` n=205 id_type=symbol source=`<repo>\results\md2\genesets.json`
- **TGFB_built_from_md2** name=`TGFB (STAR Methods built; md2)` n=54 id_type=symbol source=`<repo>\results\md2\genesets.json`
- **mmc3_Aging_signatures_Age up** name=`Aging_signatures:Age up` n=1533 id_type=symbol source=`<repo>\results\md4\genesets.json`
- **mmc3_Aging_signatures_Age down** name=`Aging_signatures:Age down` n=2007 id_type=symbol source=`<repo>\results\md4\genesets.json`

## Files

- `src/md5_common.py`, `md5_task1.py`, `md5_task2.py`, `md5_findings.py`, `md5_run.py`
- `results/md5/`
- `FINDINGS_MD5.md`
- `PROGRESS_MD5.md`

## Supersession note, appended 2026-10-02: atlas ruler − MD Δρ on 93 pseudobulks

Appended; nothing above is changed. Where this file says that the FINDINGS_MD3.md Task 2 Δρ result stands, its fibroblast-atlas ruler − MD interval no longer stands.

The ruler − MD Δρ computed with each of the 93 CELLxGENE atlas pseudobulks (a19d1667) as a unit (+0.202, 95% CI +0.006 to +0.440, `excludes_zero=True`; `results/md3/t2_delta_rho.csv`, row `cxg_ruler_minus_MD`) is **superseded** by the per-person re-analysis, in which pseudobulk scores are averaged within each of 65 donors (+0.171, 95% CI −0.118 to +0.412, `excludes_zero=False`; `results/paper_figs/cxg_true_donor_results.csv`, `analysis=true_donor`, `contrast=ruler_minus_MD_AddModuleScore`). The manuscript reports the per-person value in Table 2 (`paper/rejuvenation_readouts_preprint_11.docx`). **The per-person re-analysis was not pre-registered:** `results/paper_figs/` has no PREREG flag, and its results were written on 2026-09-23, after the pre-registered 93-pseudobulk result was known. After regrouping, the ruler − MD interval excludes zero only in GSE226189 (+0.467, 95% CI +0.239 to +0.696). The atlas ruler − (age-up − age-down) interval excludes zero in both units (+0.295, CI +0.107 to +0.487 on 93 pseudobulks; +0.372, CI +0.094 to +0.595 on 65 donors). The pre-registration flags are not edited. The reading is listed as B7 in `results/verify/WITHDRAWN_READINGS.md`.

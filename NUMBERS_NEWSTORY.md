# NUMBERS_NEWSTORY

Numbers for the restructured paper. Produced by `python src/newstory.py all`. Writes only `results/newstory/` and this file. paper/ and paper_package/ were not read.

## Part A. The MD score per state (GSE297234)

MD = `results/md3/genesets.json` "MD" (205 genes, all mapped). Seurat AddModuleScore from `src/md2_score.py` (LogNormalize scale 10,000; nbin 24; ctrl 100; seed 20260914). Labels: md2 Louvain clusters (`results/md2/louvain_obs_*.csv` + `t2_cluster_labels_*.csv`). Day-0 Fibroblast = label Fibroblast at day 0; other states pooled over days 0/3/7/10. Donors never pooled.

**Joint scoring.** Both donors' Louvain count matrices were put on one gene axis (union by Ensembl id: 30,000 genes; 28,212 shared, 846 aged-only, 942 young-only, zero-filled where a donor lacks the gene, which leaves every cell's library size and LogNormalize values unchanged). Gene means were taken over all 55,990 cells of both donors; the 24 bins and the control genes (9,483 unique) were drawn once from those means and used for every cell. Because the score of a cell depends only on its own row, this equals scoring one stacked matrix; it was run donor by donor for memory (3.8 GB free). MD genes absent from a donor's own Louvain matrix (fewer than 5 cells; zero in the joint matrix for that donor): GM00731: none; GM23815: PRSS2.

**How md4 was scored: per donor, not jointly.** md_score_mean in results/md4/t1_state_timepoint_means.csv is the md2 louvain_ams vector, computed by add_module_score on each donor's own Louvain matrix separately (src/md2_task1.py for GM00731, src/md2_task2.py for GM23815); gene means, bins and control genes differ between donors. Rerunning that per-donor scoring reproduces the stored vectors: max |diff| 0.0e+00 (GM00731), 0.0e+00 (GM23815); bins identical (True, True); per-donor control sets 9,120 and 9,690 genes.

### A1. MD score per state

Mean over cells; 95% percentile bootstrap CI over cells (B = 200, one Generator(20260918) used in order aged then young, states in table order). `per-donor` = md2 stored score (each donor scored alone). `md4` = n-weighted mean of `md_score_mean` over days in `results/md4/t1_state_timepoint_means.csv`.

Source: `results/newstory/partA_md_by_state.csv`

| donor | cell_line | state | n_cells | joint_mean | joint_ci_lo | joint_ci_hi | perdonor_mean | md4_nweighted_mean | n_cells_md4 |
|---|---|---|---|---|---|---|---|---|---|
| aged | GM00731 | Fibroblast_d0 | 4939 | 0.5151 | 0.5135 | 0.5169 | 0.5062 | 0.5062 | 4939 |
| aged | GM00731 | PartialReprog | 5832 | 0.1575 | 0.1531 | 0.1617 | 0.1505 | 0.1505 | 5832 |
| aged | GM00731 | EarlyPluripotency | 4414 | 0.1467 | 0.1414 | 0.1515 | 0.1415 | 0.1415 | 4414 |
| aged | GM00731 | Pluripotency | 439 | 0.0219 | 0.0127 | 0.0309 | 0.0216 | 0.0216 | 439 |
| aged | GM00731 | NonReprog | 5459 | 0.4479 | 0.4451 | 0.4506 | 0.4394 | 0.4394 | 5459 |
| young | GM23815 | Fibroblast_d0 | 7703 | 0.4688 | 0.4672 | 0.4706 | 0.4788 | 0.4788 | 7703 |
| young | GM23815 | PartialReprog | 6481 | 0.0825 | 0.0791 | 0.0858 | 0.0904 | 0.0904 | 6481 |
| young | GM23815 | EarlyPluripotency | 12961 | 0.1014 | 0.0991 | 0.1037 | 0.1081 | 0.1081 | 12961 |
| young | GM23815 | Pluripotency | 1249 | -0.0469 | -0.0513 | -0.0414 | -0.0451 | -0.0451 | 1249 |
| young | GM23815 | NonReprog | 1812 | 0.4133 | 0.4097 | 0.4165 | 0.4226 | 0.4226 | 1812 |

### A2. Old–young gap, drop, drop ÷ gap

gap = aged day-0 − young day-0; drop = aged day-0 − aged PartialReprog. CIs (joint only) from the same bootstrap draws, states resampled independently.

Source: `results/newstory/partA_gap_drop.csv`

| scoring | aged_d0 | young_d0 | aged_PartialReprog | gap | gap_CI | drop | drop_CI | drop_over_gap | ratio_CI |
|---|---|---|---|---|---|---|---|---|---|
| joint | 0.5151 | 0.4688 | 0.1575 | 0.0463 | [0.0434, 0.0487] | 0.3576 | [0.3530, 0.3628] | 7.72 | [7.33, 8.21] |
| per-donor (md2 stored) | 0.5062 | 0.4788 | 0.1505 | 0.0274 |  | 0.3557 |  | 12.99 |  |
| md4 n-weighted | 0.5062 | 0.4788 | 0.1505 | 0.0274 |  | 0.3557 |  | 12.99 |  |

## Part B. Extra geometry numbers (point estimates)

Recomputed with `src/genespace_run.py` functions (`map_gene_spaces`, `same_point_block`, `gtex_point_block`, `draw_baseline_sets`) on the same panels. Every recomputed quantity that is also stored in `results/genespace/` matches it: max |diff| 5.7e-14 (`results/newstory/partB_repro_checks.csv`). Distances are Euclidean in frozen-ruler z units within each gene space.

### B1. Same-platform test: starting and final distance

Forward: start = aged day-0 half A (O), target = young day-0 half A (Y), cells = aged PartialReprog. Reverse: start = Y, target = O, cells = young PartialReprog. Starting distance = |O − Y| (same for both tests). Pluripotency and NonReprog (aged, forward frame) are context rows.

Source: `results/newstory/partB_same_distances.csv`

| space | n_genes | test | n_cells | dist_start | dist_final | final_over_start | progress | perp | delta |
|---|---|---|---|---|---|---|---|---|---|
| FULL | 23485 | forward | 5832 | 130.688 | 261.582 | 2.0016 | 0.6744 | 1.9749 | 130.894 |
| FULL | 23485 | reverse | 6481 | 130.688 | 294.030 | 2.2499 | -0.0573 | 1.9860 | 163.342 |
| FULL | 23485 | Pluri | 439 | 130.688 | 512.794 | 3.9238 | 1.4032 | 3.9030 | 382.106 |
| FULL | 23485 | NonReprog | 5459 | 130.688 | 191.019 | 1.4616 | 0.5549 | 1.3922 | 60.331 |
| MD | 205 | forward | 5832 | 14.598 | 31.686 | 2.1706 | 0.7901 | 2.1604 | 17.088 |
| MD | 205 | reverse | 6481 | 14.598 | 39.479 | 2.7045 | -0.2295 | 2.4089 | 24.881 |
| MD | 205 | Pluri | 439 | 14.598 | 56.885 | 3.8969 | 1.3273 | 3.8831 | 42.288 |
| MD | 205 | NonReprog | 5459 | 14.598 | 20.683 | 1.4169 | 0.6330 | 1.3685 | 6.086 |
| AGE | 3089 | forward | 5832 | 48.338 | 92.944 | 1.9228 | 0.8221 | 1.9145 | 44.605 |
| AGE | 3089 | reverse | 6481 | 48.338 | 111.381 | 2.3042 | -0.2304 | 1.9482 | 63.042 |
| AGE | 3089 | Pluri | 439 | 48.338 | 174.723 | 3.6146 | 1.3118 | 3.6011 | 126.384 |
| AGE | 3089 | NonReprog | 5459 | 48.338 | 71.141 | 1.4717 | 0.6016 | 1.4168 | 22.802 |

Mixture control (aged half B + young half B cells, f = young fraction; start O, target Y):

Source: `results/newstory/partB_mixture.csv`

| space | f | n_aged | n_young | dist_start | dist_final | delta | rel_delta | progress | perp |
|---|---|---|---|---|---|---|---|---|---|
| FULL | 0.00 | 2470 | 0 | 124.240 | 124.308 | 0.068 | 0.0005 | 0.1146 | 0.4660 |
| FULL | 0.10 | 2223 | 247 | 124.240 | 116.370 | -7.870 | -0.0633 | 0.1824 | 0.4571 |
| FULL | 0.25 | 1852 | 618 | 124.240 | 105.327 | -18.913 | -0.1522 | 0.2829 | 0.4522 |
| FULL | 0.50 | 1235 | 1235 | 124.240 | 89.797 | -34.443 | -0.2772 | 0.4408 | 0.4579 |
| MD | 0.00 | 2470 | 0 | 14.434 | 14.427 | -0.007 | -0.0005 | 0.0104 | 0.1403 |
| MD | 0.10 | 2223 | 247 | 14.434 | 13.239 | -1.195 | -0.0828 | 0.0955 | 0.1525 |
| MD | 0.25 | 1852 | 618 | 14.434 | 11.758 | -2.676 | -0.1854 | 0.2074 | 0.1882 |
| MD | 0.50 | 1235 | 1235 | 14.434 | 9.291 | -5.143 | -0.3563 | 0.3883 | 0.2005 |
| AGE | 0.00 | 2470 | 0 | 46.551 | 46.982 | 0.431 | 0.0093 | 0.0576 | 0.3613 |
| AGE | 0.10 | 2223 | 247 | 46.551 | 43.078 | -3.473 | -0.0746 | 0.1477 | 0.3606 |
| AGE | 0.25 | 1852 | 618 | 46.551 | 38.377 | -8.174 | -0.1756 | 0.2640 | 0.3715 |
| AGE | 0.50 | 1235 | 1235 | 46.551 | 31.668 | -14.882 | -0.3197 | 0.4367 | 0.3814 |

### B2. GTEx test: distances to the GTEx centroids

GTEx young centroid = mean z of GTEx fibroblast donors aged 20–39 (n = 113); old = 60–79 (n = 231). Day-0 pseudobulk = all day-0 Fibroblast cells of the donor, TMM with that donor's state rows (the GTEx-test panel, not the same-platform half A).

Source: `results/newstory/partB_gtex_distances.csv`

| space | n_genes | cell_line | dist_d0_to_gtex_young | dist_d0_to_gtex_old | gtex_young_to_old | young_to_old_over_d0_to_young |
|---|---|---|---|---|---|---|
| FULL | 23485 | GM00731 | 467.450 | 464.431 | 32.284 | 0.0691 |
| FULL | 23485 | GM23815 | 457.031 | 454.287 | 32.284 | 0.0706 |
| MD | 205 | GM00731 | 41.988 | 41.565 | 3.288 | 0.0783 |
| MD | 205 | GM23815 | 40.405 | 40.175 | 3.288 | 0.0814 |
| AGE | 3089 | GM00731 | 164.036 | 163.624 | 12.296 | 0.0750 |
| AGE | 3089 | GM23815 | 158.695 | 158.115 | 12.296 | 0.0775 |

Change in distance from day 0, aged donor (positive = farther; Figure N1-C uses the MD rows):

| space | state | n_cells | delta_young | rel_delta_young | delta_old | rel_delta_old |
|---|---|---|---|---|---|---|
| FULL | PartialReprog | 5832 | 9.370 | 0.0200 | 10.857 | 0.0234 |
| FULL | EarlyPluripotency | 4414 | 30.456 | 0.0652 | 32.695 | 0.0704 |
| FULL | Pluripotency | 439 | 184.049 | 0.3937 | 188.313 | 0.4055 |
| FULL | NonReprog | 5459 | 8.932 | 0.0191 | 9.818 | 0.0211 |
| MD | PartialReprog | 5832 | 6.461 | 0.1539 | 7.002 | 0.1684 |
| MD | EarlyPluripotency | 4414 | 7.447 | 0.1774 | 8.025 | 0.1931 |
| MD | Pluripotency | 439 | 28.216 | 0.6720 | 29.236 | 0.7034 |
| MD | NonReprog | 5459 | -1.400 | -0.0333 | -1.320 | -0.0318 |
| AGE | PartialReprog | 5832 | -5.439 | -0.0332 | -5.925 | -0.0362 |
| AGE | EarlyPluripotency | 4414 | 1.686 | 0.0103 | 2.093 | 0.0128 |
| AGE | Pluripotency | 439 | 44.720 | 0.2726 | 45.931 | 0.2807 |
| AGE | NonReprog | 5459 | -1.015 | -0.0062 | -0.742 | -0.0045 |

### B3. Random gene sets (Amendment 1: matched on 24 GTEx-mean bins × zero/nonzero; 200 per size)

Per-set values were not stored, so they were recomputed from `results/genespace/size_baseline_sets.npz` (`A1_mu_zero__MD`, `A1_mu_zero__AGE`; redrawing with `draw_baseline_sets` gives identical sets). They reproduce the stored percentiles and positions in `results/genespace/size_baseline.csv`: max |diff| 7.1e-15.

Per-set values: `results/newstory/partB_random_sets_per_set.csv`. Percentile check: `results/newstory/partB_random_sets_repro.csv`.

| space_size | test | state | stat | real | p2_5 | p50 | p97_5 | position | max_abs_diff |
|---|---|---|---|---|---|---|---|---|---|
| MD | same | forward | progress | 0.7901 | 0.6552 | 0.9800 | 1.3870 | 0.170 | 1.1e-16 |
| MD | same | forward | perp | 2.1604 | 1.8911 | 2.2914 | 2.6973 | 0.255 | 0.0e+00 |
| MD | same | forward | delta | 17.0880 | 12.2020 | 15.0116 | 18.2738 | 0.860 | 3.6e-15 |
| MD | same | reverse | progress | -0.2295 | -0.8672 | -0.4686 | -0.1034 | 0.910 | 5.6e-17 |
| MD | same | reverse | perp | 2.4089 | 1.9947 | 2.4255 | 2.8288 | 0.465 | 2.2e-16 |
| MD | same | reverse | delta | 24.8813 | 18.0005 | 21.3186 | 25.0194 | 0.965 | 0.0e+00 |
| MD | mixture | f=0.50 | delta | -5.1430 | -5.0297 | -3.9240 | -2.9148 | 0.015 | 4.4e-16 |
| AGE | same | forward | progress | 0.8221 | 0.6112 | 0.6852 | 0.7600 | 1.000 | 0.0e+00 |
| AGE | same | forward | perp | 1.9145 | 1.8834 | 1.9663 | 2.0655 | 0.140 | 0.0e+00 |
| AGE | same | forward | delta | 44.6052 | 44.8707 | 47.5853 | 51.0443 | 0.020 | 0.0e+00 |
| AGE | same | reverse | progress | -0.2304 | -0.1651 | -0.0956 | -0.0262 | 0.000 | 8.3e-17 |
| AGE | same | reverse | perp | 1.9482 | 1.9025 | 1.9947 | 2.0691 | 0.195 | 0.0e+00 |
| AGE | same | reverse | delta | 63.0425 | 58.2049 | 61.1386 | 63.7664 | 0.910 | 7.1e-15 |
| AGE | mixture | f=0.50 | delta | -14.8824 | -15.0122 | -13.7070 | -12.3806 | 0.050 | 0.0e+00 |

How many random sets end closer to the target (change in distance < 0):

Source: `results/newstory/partB_random_sets_counts.csv`

| space_size | n_sets | fwd_closer | rev_closer | mix50_closer | fwd_min_final_over_start | real_fwd_final_over_start |
|---|---|---|---|---|---|---|
| MD | 200 | 0 | 0 | 200 | 1.8142 | 2.1706 |
| AGE | 200 | 0 | 0 | 200 | 1.8719 | 1.9228 |

Range of per-set values (min / median / max):

| space_size | fwd_progress_min | fwd_progress_median | fwd_progress_max | fwd_perp_min | fwd_perp_median | fwd_perp_max | rev_progress_min | rev_progress_median | rev_progress_max | rev_perp_min | rev_perp_median | rev_perp_max | mix50_delta_min | mix50_delta_median | mix50_delta_max | mix50_rel_delta_min | mix50_rel_delta_median | mix50_rel_delta_max |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| AGE | 0.5841 | 0.6852 | 0.7837 | 1.8443 | 1.9663 | 2.0828 | -0.2056 | -0.0956 | -0.0053 | 1.8728 | 1.9947 | 2.0959 | -15.6031 | -13.7070 | -12.1161 | -0.3289 | -0.2975 | -0.2670 |
| MD | 0.4847 | 0.9800 | 1.5122 | 1.8141 | 2.2914 | 2.9476 | -0.9605 | -0.4686 | -0.0697 | 1.8764 | 2.4255 | 3.1362 | -6.1937 | -3.9240 | -2.3636 | -0.4765 | -0.3410 | -0.2372 |

## Figures

All in `results/newstory/` (PNG 300 dpi, 7 in wide; PDF; plotting data CSV).

- Figure N1, "MD's own genes": `fig_N1_md_own_genes.png`, `fig_N1_md_own_genes.pdf`, `fig_N1_md_own_genes_data.csv`
  - A: MD score (joint scoring) for aged day-0, PartialReprog, EarlyPluripotency, Pluripotency; dashed line = young day-0.
  - B: same-platform plane in the MD genes (x = progress, y = sideways, axis lengths).
  - C: change in distance to GTEx young / old in the MD genes, aged donor.
- Figure N2, "Any genes": `fig_N2_any_genes.png`, `fig_N2_any_genes.pdf`, `fig_N2_any_genes_data.csv`
  - A / B: forward / reverse plane; real MD, age, all genes; 200 random sets per size.
  - C: mixture control, change in distance ÷ starting distance vs f, three spaces.

Other outputs: `partA_cell_scores.npz` (per-cell joint and per-donor MD scores, labels, shared control columns), `partA_meta.json`, `partB_meta.json`, `run_log.txt`.

## Integrity

Files outside `results/newstory/` (results/, src/, data/processed/md2, data/processed/md4, root *.md) compared by size and mtime before and after the run: changed = none.


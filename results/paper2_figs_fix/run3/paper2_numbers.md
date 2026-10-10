# PAPER2FIGS-fix: manuscript numbers against the stored results

Written by `python src/paper2_figs_fix_run.py`. Rows 1-75 are PAPER2FIGS's; row 76 is new. Each value is read from the file and line given; MATCH
means the stated value equals the stored value rounded to the decimals as stated (integers
exact). Line numbers are 1-based; for CSVs line 1 is the header. Nothing was adjusted to match.

**76 MATCH, 0 MISMATCH, of 76.**

| # | value | stated | stored | file:line | result | note |
|---|---|---|---|---|---|---|
| 1 | cells in the h5ad | 36514 | 36514 | `results/senginv/addendum_meta.json:4` | **MATCH** |  |
| 2 | PD14 cells | 581 | 581 | `results/senginv/ADDENDUM_SENGINV.md:32` | **MATCH** | group table of the SENGINV addendum (a markdown results file) |
| 3 | PD26 cells | 1880 | 1880 | `results/senginv/ADDENDUM_SENGINV.md:33` | **MATCH** | group table of the SENGINV addendum (a markdown results file) |
| 4 | PD32 cells | 4348 | 4348 | `results/senginv/ADDENDUM_SENGINV.md:34` | **MATCH** | group table of the SENGINV addendum (a markdown results file) |
| 5 | barcode _A cells | 14333 | 14333 | `results/senginv/ADDENDUM_SENGINV.md:30` | **MATCH** | group table of the SENGINV addendum (a markdown results file) |
| 6 | barcode _I cells | 15372 | 15372 | `results/senginv/ADDENDUM_SENGINV.md:31` | **MATCH** | group table of the SENGINV addendum (a markdown results file) |
| 7 | NT_CRA cells | 666 | 666 | `results/seng2/d1_summary.json:74` | **MATCH** |  |
| 8 | NT_CRI cells | 1479 | 1479 | `results/seng2/d1_summary.json:75` | **MATCH** |  |
| 9 | Sengstack ruler genes present | 16781 | 16781 | `results/axang/d2_summary.json:3` | **MATCH** |  |
| 10 | ruler genes in total (present + absent) | 23485 | 23485 | `results/axang/d2_summary.json:3-4` | **MATCH** | sum of two stored counts; also results/axang/d1_summary.json meta.n_ruler |
| 11 | offset/||a||, half A, full-gene, A | 1.0532 | 1.05316 | `results/seng3/d3_offset.csv:4` | **MATCH** | thinned median |
| 12 | cos(offset, a), half A, full-gene, A | 0.1123 | 0.112275 | `results/seng3/d3_offset.csv:4` | **MATCH** | thinned median |
| 13 | offset/||a||, half A, full-gene, I | 0.7298 | 0.729838 | `results/seng3/d3_offset.csv:5` | **MATCH** | thinned median |
| 14 | cos(offset, a), half A, full-gene, I | 0.0708 | 0.0708182 | `results/seng3/d3_offset.csv:5` | **MATCH** | thinned median |
| 15 | cos(target, aging), Full-gene, A | 0.6445 | 0.644549 | `results/seng3/d3_offset.csv:4` | **MATCH** | thinned median |
| 16 | cos(target, aging), Full-gene, I | 0.7934 | 0.793429 | `results/seng3/d3_offset.csv:5` | **MATCH** | thinned median |
| 17 | cos(target, aging), S2A_k20, A | 0.811375 | 0.811375 | `results/seng3/d3_offset.csv:2` | **MATCH** | thinned median |
| 18 | cos(target, aging), S2A_k20, I | 0.887005 | 0.887005 | `results/seng3/d3_offset.csv:3` | **MATCH** | thinned median |
| 19 | cos(a, u), full-gene | -0.040273 | -0.0402729 | `results/seng4/d1_summary.json:51` | **MATCH** |  |
| 20 | full-gene null median, lowest (origin i) | 20.607 | 20.6068 | `results/seng3/d1_null_validation.csv:31` | **MATCH** | CRI 157 cells; lowest of the 8 origin (i) S1 rows |
| 21 | full-gene null median, highest (origin i) | 33.654 | 33.6543 | `results/seng3/d1_null_validation.csv:3` | **MATCH** | CRA 96 cells; highest of the 8 origin (i) S1 rows |
| 22 | full-gene effect-free draws negative, every S1 row | 0 | 0 | `results/seng3/d1_null_validation.csv:3-33 (S1 rows)` | **MATCH** | largest n_eval_negative over the 16 S1 rows, of 200 |
| 23 | compressed firing rate, lowest (i) | 0.04 | 0.04 | `results/seng3/d1_null_validation.csv:10` | **MATCH** | CRA 128 |
| 24 | compressed firing rate, highest (i) | 0.085 | 0.085 | `results/seng3/d1_null_validation.csv:6` | **MATCH** | CRA 107 |
| 25 | compressed firing rate, highest (ii) | 0.11 | 0.11 | `results/seng3/d1_summary.json:8` | **MATCH** |  |
| 26 | plant power at 25%, lowest (S2A_k20, i) | 0.565 | 0.565 | `results/seng3/d2_plant.csv:3` | **MATCH** | CRA 96 |
| 27 | plant power at 25%, highest (S2A_k20, i) | 0.74 | 0.74 | `results/seng3/d2_plant.csv:115` | **MATCH** | CRI 157 |
| 28 | plant power at 50%, lowest (S2A_k20, i) | 0.94 | 0.94 | `results/seng3/d2_plant.csv:20` | **MATCH** | CRA 107 |
| 29 | plant power at 50%, highest (S2A_k20, i) | 1 | 1 | `results/seng3/d2_plant.csv:100` | **MATCH** | CRI 128 |
| 30 | % of the aging gap's length kept by S2A_k20 | 45.1 | 45.0785 | `results/seng3/d2_summary.json:11-12` | **MATCH** | ratio of two stored norms (40.6380 / 90.1493); the 45.1% itself is stated only at FINDINGS_SENG3.md:662 |
| 31 | full-gene: of the four called, origin (i) | 1 | 1 | `results/seng2/d2_summary.json:226` | **MATCH** |  |
| 32 | full-gene: of the four called, origin (ii) | 0 | 0 | `results/seng2/d2_summary.json:227` | **MATCH** |  |
| 33 | compressed called of 340, origin (i) | 80 | 80 | `results/seng3/d4_summary.json:7` | **MATCH** |  |
| 34 | compressed called of 340, origin (ii) | 73 | 73 | `results/seng3/d4_summary.json:57` | **MATCH** |  |
| 35 | CRA_E2F3 delta to PD14, origin (i) | -19.211 | -19.211 | `results/seng3/d4_table.csv:32` | **MATCH** |  |
| 36 | CRA_EZH2 delta to PD14, origin (i) | -13.211 | -13.2105 | `results/seng3/d4_table.csv:40` | **MATCH** |  |
| 37 | CRI_STAT3 delta to PD14, origin (i) | -10.334 | -10.3342 | `results/seng3/d4_table.csv:308` | **MATCH** |  |
| 38 | CRI_ZFX delta to PD14, origin (i) | -19.936 | -19.9362 | `results/seng3/d4_table.csv:332` | **MATCH** |  |
| 39 | call rate, cell-count decile 1 (origin i) | 0.056 | 0.0556 | `results/seng3/d5_confounds.json:5` | **MATCH** | stored to 4 dp; called/n = 2/36 = 0.0556 |
| 40 | call rate, cell-count decile 10 (origin i) | 0.545 | 0.5455 | `results/seng3/d5_confounds.json:132` | **MATCH** | stored to 4 dp; called/n = 18/33 = 0.5455 |
| 41 | depth band 96-157: called (i) | 24 | 24 | `results/seng3/d5_confounds.json:4` | **MATCH** |  |
| 42 | depth band 96-157: n | 60 | 60 | `results/seng3/d5_confounds.json:3` | **MATCH** |  |
| 43 | CRA_E2F3 delta percentile in the depth band | 3.3 | 3.3 | `results/seng3/d5_confounds.json:15` | **MATCH** |  |
| 44 | CRA_EZH2 delta percentile in the depth band | 15 | 15 | `results/seng3/d5_confounds.json:25` | **MATCH** |  |
| 45 | CRI_STAT3 delta percentile in the depth band | 21.7 | 21.7 | `results/seng3/d5_confounds.json:35` | **MATCH** |  |
| 46 | CRI_ZFX delta percentile in the depth band | 1.7 | 1.7 | `results/seng3/d5_confounds.json:45` | **MATCH** |  |
| 47 | rho(prolif. shift, delta i) | -0.8639 | -0.863943 | `results/seng3/d5_confounds.json:158` | **MATCH** |  |
| 48 | angle(P1, a), CRA | 19.66 | 19.6596 | `results/seng4/d1_vs_aging.csv:2` | **MATCH** |  |
| 49 | angle(P1, a), CRI | 20.56 | 20.5629 | `results/seng4/d1_vs_aging.csv:3` | **MATCH** |  |
| 50 | % of a removed by P1 (CRA) | 66.4 | 66.3568 | `results/seng4/d1_vs_aging.csv:2` | **MATCH** |  |
| 51 | angle(P3, a) | 61.68 | 61.6807 | `results/seng4/d1_vs_aging.csv:6` | **MATCH** |  |
| 52 | % of a removed by P3 | 12 | 11.9683 | `results/seng4/d1_vs_aging.csv:6` | **MATCH** |  |
| 53 | p50 before residualizing (SENG3) | 0.94 | 0.94 | `results/seng4/d2_summary.json:47` | **MATCH** |  |
| 54 | p50 after projecting out P1 | 0.01 | 0.01 | `results/seng4/d2_summary.json:3` | **MATCH** |  |
| 55 | p50 control (nothing removed, reserved-half origin) | 0.895 | 0.895 | `results/seng4/d2_summary.json:8` | **MATCH** |  |
| 56 | SENG4 D3 gate | SPLIT | SPLIT | `results/seng4/d3_summary.json:2` | **MATCH** |  |
| 57 | SPLIT: of the four called, P1 origin (i) | 2 | 2 | `results/seng4/d3_summary.json:15` | **MATCH** |  |
| 58 | SPLIT: of the four called, P1 origin (ii) | 0 | 0 | `results/seng4/d3_summary.json:16` | **MATCH** |  |
| 59 | angle A1-A3, S2A_k20 | 101.5 | 101.496 | `results/axang/d1_angles.csv:69` | **MATCH** |  |
| 60 | angle A2-A3, S2A_k20 | 100.79 | 100.791 | `results/axang/d1_angles.csv:79` | **MATCH** |  |
| 61 | A1 vs C1_gtex, angle from parallel | 35.83 | 35.825 | `results/axang/d1_summary.json:84` | **MATCH** |  |
| 62 | cos(A1, C1_gtex) | -0.8108 | -0.810808 | `results/axang/d1_summary.json:4` | **MATCH** |  |
| 63 | GTEx plant power at 25%, lowest | 0.945 | 0.945 | `results/axang/d3a_plant.csv:3` | **MATCH** | CRA 96 |
| 64 | GTEx plant power at 25%, highest | 0.99 | 0.99 | `results/axang/d3a_plant.csv:27` | **MATCH** | CRA 157 |
| 65 | CRA_E2F3 delta to GTEx young | 11.448 | 11.4478 | `results/axang/d3b_table.csv:32` | **MATCH** |  |
| 66 | CRA_EZH2 delta to GTEx young | 9.022 | 9.02184 | `results/axang/d3b_table.csv:40` | **MATCH** |  |
| 67 | CRI_STAT3 delta to GTEx young | 5.327 | 5.32743 | `results/axang/d3b_table.csv:308` | **MATCH** |  |
| 68 | CRI_ZFX delta to GTEx young | 15.336 | 15.3364 | `results/axang/d3b_table.csv:332` | **MATCH** |  |
| 69 | GTEx test: closer, of 340 | 17 | 17 | `results/axang/d3_summary.json:10` | **MATCH** |  |
| 70 | GTEx test: called, of 340 | 5 | 5 | `results/axang/d3_summary.json:11` | **MATCH** |  |
| 71 | PD32 distance to GTEx young / gap | 14.4 | 14.4048 | `results/axang/d2_summary.json:21` | **MATCH** |  |
| 72 | % of GTEx age gap in genes Sengstack does not measure | 21.97 | 21.9724 | `results/axang/d2_summary.json:5` | **MATCH** | Sengstack row. The identical Lu row in this file is superseded (CORRECTION.md C1); Lu corrected = 8.97% in results/axang/correction_values.json |
| 73 | ruler: score younger, of 340 | 181 | 181 | `results/seng3/d4_summary.json:259` | **MATCH** |  |
| 74 | ruler: younger, not called (i) | 116 | 116 | `results/seng3/d4_summary.json:260` | **MATCH** |  |
| 75 | ruler: younger, not called (ii) | 120 | 120 | `results/seng3/d4_summary.json:271` | **MATCH** |  |
| 76 | of the 80 called against early-passage cells (origin i): number called against GTEx young, and number with delta_gtex_young > 0 | 0 and 80 | 0 and 80 | `results/paper2_figs_fix/run3/Fig4_data.csv:5,28,30,32,38,40,43,48,51,57,73,75,78,81,89,95,115,117,119-120,124,130,137,139-140,143,146,148,150,157,170-171,174,179,181-182,184,189-190,194,198,200-202,208,211-212,218-219,222-223,226-227,231,249,252,254-255,261,269,275-276,281,287,290,292,295,299-302,308,313,318,320,324,326,328,331-332` | **MATCH** | Count over plotted rows (arithmetic on stored values). Subset = called_early_passage_i True: 80 rows (stated 80, MATCH). Chance baseline for the GTEx count: 5th-percentile floor, 4.0 of 80 expected with no effect. Each row's lines in results/seng3/d4_table.csv and results/axang/d3b_table.csv: results/paper2_figs_fix/run3/new_row_sources.csv. |

## Angle check (author's instruction, Q2 of PLAN.flag)

Every plotted angle from parallel (Fig 3A, Fig S3) against the FINDINGS_AXANG.md §5 table, to 2 dp.
Plotted values for C1_lu come from the corrected CSV, so they are expected to differ from §5,
which holds the superseded C1_lu values; those rows are also checked against the corrected table
in results/axang/CORRECTION.md.

**20 of 24 match §5; 4 do not.**

| figure | pair | plotted | §5 | §5 line | match §5 | CORRECTION.md | line | match corrected |
|---|---|---|---|---|---|---|---|---|
| Fig3A | A1 vs C1_south | 70.22 | 70.22 | `FINDINGS_AXANG.md:574` | yes |  |  |  |
| Fig3A | A2 vs C1_south | 80.49 | 80.49 | `FINDINGS_AXANG.md:581` | yes |  |  |  |
| Fig3A | A3 vs C1_south | 27.30 | 27.30 | `FINDINGS_AXANG.md:588` | yes |  |  |  |
| Fig3A | A4 vs C1_south | 64.48 | 64.48 | `FINDINGS_AXANG.md:595` | yes |  |  |  |
| Fig3A | A1 vs C1_seng_CRA | 73.19 | 73.19 | `FINDINGS_AXANG.md:575` | yes |  |  |  |
| Fig3A | A2 vs C1_seng_CRA | 76.44 | 76.44 | `FINDINGS_AXANG.md:582` | yes |  |  |  |
| Fig3A | A3 vs C1_seng_CRA | 19.66 | 19.66 | `FINDINGS_AXANG.md:589` | yes |  |  |  |
| Fig3A | A4 vs C1_seng_CRA | 65.40 | 65.40 | `FINDINGS_AXANG.md:596` | yes |  |  |  |
| Fig3A | A1 vs C1_seng_CRI | 74.40 | 74.40 | `FINDINGS_AXANG.md:576` | yes |  |  |  |
| Fig3A | A2 vs C1_seng_CRI | 78.05 | 78.05 | `FINDINGS_AXANG.md:583` | yes |  |  |  |
| Fig3A | A3 vs C1_seng_CRI | 20.56 | 20.56 | `FINDINGS_AXANG.md:590` | yes |  |  |  |
| Fig3A | A4 vs C1_seng_CRI | 71.83 | 71.83 | `FINDINGS_AXANG.md:597` | yes |  |  |  |
| Fig3A | A1 vs C1_lu | 70.36 | 70.95 | `FINDINGS_AXANG.md:577` | **NO** | 70.36 | `results/axang/CORRECTION.md:77` | yes |
| Fig3A | A2 vs C1_lu | 82.45 | 81.98 | `FINDINGS_AXANG.md:584` | **NO** | 82.45 | `results/axang/CORRECTION.md:78` | yes |
| Fig3A | A3 vs C1_lu | 31.94 | 32.63 | `FINDINGS_AXANG.md:591` | **NO** | 31.94 | `results/axang/CORRECTION.md:79` | yes |
| Fig3A | A4 vs C1_lu | 57.97 | 57.63 | `FINDINGS_AXANG.md:598` | **NO** | 57.97 | `results/axang/CORRECTION.md:80` | yes |
| FigS3 | A1 vs C2 | 48.79 | 48.79 | `FINDINGS_AXANG.md:578` | yes |  |  |  |
| FigS3 | A2 vs C2 | 84.60 | 84.60 | `FINDINGS_AXANG.md:585` | yes |  |  |  |
| FigS3 | A3 vs C2 | 61.68 | 61.68 | `FINDINGS_AXANG.md:592` | yes |  |  |  |
| FigS3 | A4 vs C2 | 77.06 | 77.06 | `FINDINGS_AXANG.md:599` | yes |  |  |  |
| FigS3 | A1 vs C1_gtex | 35.83 | 35.83 | `FINDINGS_AXANG.md:572` | yes |  |  |  |
| FigS3 | A2 vs C1_gtex | 85.51 | 85.51 | `FINDINGS_AXANG.md:579` | yes |  |  |  |
| FigS3 | A3 vs C1_gtex | 55.38 | 55.38 | `FINDINGS_AXANG.md:586` | yes |  |  |  |
| FigS3 | A4 vs C1_gtex | 83.68 | 83.68 | `FINDINGS_AXANG.md:593` | yes |  |  |  |

Values that do not match §5:

- Fig3A A1 vs C1_lu: plotted 70.36°, §5 says 70.95° (`FINDINGS_AXANG.md:577`); superseded value, corrected 70.36° at `results/axang/CORRECTION.md:77`, which the plotted value matches
- Fig3A A2 vs C1_lu: plotted 82.45°, §5 says 81.98° (`FINDINGS_AXANG.md:584`); superseded value, corrected 82.45° at `results/axang/CORRECTION.md:78`, which the plotted value matches
- Fig3A A3 vs C1_lu: plotted 31.94°, §5 says 32.63° (`FINDINGS_AXANG.md:591`); superseded value, corrected 31.94° at `results/axang/CORRECTION.md:79`, which the plotted value matches
- Fig3A A4 vs C1_lu: plotted 57.97°, §5 says 57.63° (`FINDINGS_AXANG.md:598`); superseded value, corrected 57.97° at `results/axang/CORRECTION.md:80`, which the plotted value matches

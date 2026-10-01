# REVIEW_CHECK

Read-only check, 2026-10-02. Values computed from stored CSV/TXT files only; no analysis was rerun and no other file was modified.

## 1. `results/genespace/same_stats.csv`: rho = d_norm / v_norm vs 2*cos

`rho_lt_2cos` = whether rho < 2*cos. Cross-check: sqrt(progress² + perp²) equals d_norm / v_norm to 6 decimals in all 12 rows.

| space | test | rho | two_cos | rho < two_cos | delta (stored) |
|---|---|---|---|---|---|
| FULL | forward | 2.086903 | 0.646356 | False | 130.893565 |
| FULL | reverse | 1.986781 | -0.057656 | False | 163.341673 |
| FULL | Pluri | 4.147585 | 0.676619 | False | 382.105569 |
| FULL | NonReprog | 1.498723 | 0.740487 | False | 60.330758 |
| MD | forward | 2.300376 | 0.686941 | False | 17.087969 |
| MD | reverse | 2.419771 | -0.189649 | False | 24.881302 |
| MD | Pluri | 4.103719 | 0.646900 | False | 42.287779 |
| MD | NonReprog | 1.507838 | 0.839606 | False | 6.085638 |
| AGE | forward | 2.083558 | 0.789121 | False | 44.605176 |
| AGE | reverse | 1.961772 | -0.234857 | False | 63.042484 |
| AGE | Pluri | 3.832582 | 0.684548 | False | 126.384207 |
| AGE | NonReprog | 1.539223 | 0.781719 | False | 22.802353 |

- Rows with rho < two_cos: 0 of 12.
- rho > 2 for the three forward rows: FULL 2.086903 (True), MD 2.300376 (True), AGE 2.083558 (True). Confirmed: 3 of 3.

## 2. `results/newstory/partB_random_sets_per_set.csv`

fwd_rho = sqrt(fwd_progress² + fwd_perp²); fwd_cos = fwd_progress / fwd_rho.

| subset | n sets | fwd_rho < 2*fwd_cos | fwd_rho < 2 | fwd_rho min | fwd_rho max |
|---|---|---|---|---|---|
| all | 400 | 0 | 8 | 1.952238 | 3.312845 |
| space_size = AGE | 200 | 0 | 7 | 1.965617 | 2.202471 |
| space_size = MD | 200 | 0 | 1 | 1.952238 | 3.312845 |

Also from the stored columns: fwd_delta < 0 in 0 of 400; fwd_final_over_start < 1 in 0 of 400.

## 3. `results/genespace/same_posctrl.csv`: rho per mixture fraction

This file has no d_norm / v_norm columns, so rho = sqrt(progress² + perp²) (the identity cross-checked in item 1).

| space | f | rho | two_cos | rho < two_cos | delta (stored) |
|---|---|---|---|---|---|
| FULL | 0.00 | 0.479865 | 0.477582 | False | 0.068051 |
| FULL | 0.10 | 0.492131 | 0.741410 | True | -7.870024 |
| FULL | 0.25 | 0.533455 | 1.060740 | True | -18.912864 |
| FULL | 0.50 | 0.635638 | 1.387012 | True | -34.443006 |
| MD | 0.00 | 0.140634 | 0.147393 | True | -0.006862 |
| MD | 0.10 | 0.179931 | 1.062047 | True | -1.194962 |
| MD | 0.25 | 0.280112 | 1.481187 | True | -2.676180 |
| MD | 0.50 | 0.437010 | 1.777156 | True | -5.142994 |
| AGE | 0.00 | 0.365874 | 0.315013 | False | 0.431130 |
| AGE | 0.10 | 0.389689 | 0.758277 | True | -3.472675 |
| AGE | 0.25 | 0.455783 | 1.158636 | True | -8.173837 |
| AGE | 0.50 | 0.579823 | 1.506306 | True | -14.882392 |

- Positive-control rho range, all rows: 0.140634 to 0.635638.
- By space: FULL 0.479865 to 0.635638; MD 0.140634 to 0.437010; AGE 0.365874 to 0.579823.
- Forward-test rho (item 1), for comparison: FULL 2.086903, MD 2.300376, AGE 2.083558.

## 4. `PREREG*.flag` per directory under `results/`

40 top-level directories. Of these, 23 contain at least one `PREREG*.flag` and 17 contain none. None of the 25 subdirectories (for example `*/figures`, `lowdim/cache`, `survey/opened/*`) contains one.

| directory | PREREG*.flag files |
|---|---|
| brain | none |
| brain_phase1 | none |
| external | PREREG_20260917.flag |
| fibro | PREREG_STAGE1.flag, PREREG_STAGE2.flag |
| fibro2 | PREREG_TA.flag, PREREG_TB.flag, PREREG_TC.flag |
| fibro3 | PREREG_TASK1.flag, PREREG_TASK2.flag, PREREG_TASK3.flag, PREREG_TASK4.flag |
| figures | none |
| figures_final | none |
| genespace | PREREG.flag |
| genespace_diag | none |
| geometry | none |
| gtex | PREREG_20260917.flag |
| lowdim | none |
| md | PREREG_STAGE2.flag |
| md2 | PREREG_TASK1.flag, PREREG_TASK2.flag, PREREG_TASK3.flag |
| md3 | PREREG_TASK1.flag, PREREG_TASK2.flag, PREREG_TASK3.flag |
| md4 | PREREG_TASK1.flag, PREREG_TASK2.flag, PREREG_TASK3.flag |
| md5 | PREREG_TASK1.flag, PREREG_TASK2.flag, PREREG_TASK3.flag |
| newstory | none |
| paper_figs | none |
| plane | PREREG_TASK1.flag, PREREG_TASK2.flag, PREREG_TASK3.flag |
| posctrl | none |
| relabel_check | PREREG.flag |
| rerun | none |
| same | PREREG.flag |
| scoping | none |
| seng | PREREG.flag |
| seng_why | none |
| south | PREREG.flag |
| south2 | PREREG.flag |
| south3 | PREREG.flag, PREREG_SOUTH2_inherited.flag |
| south5 | PREREG.flag |
| south6 | PREREG.flag |
| survey | PREREG.flag |
| tables | none |
| target | none |
| tissue | none |
| toward | PREREG.flag |
| trajectory | none |
| verify | PREREG.flag |

- `results/brain/`: no `PREREG*.flag` (no `.flag` file of any name; 16 files plus 4 in `figures/`).
- `results/brain_phase1/`: no `PREREG*.flag`. It contains one other `.flag` file, `p1_DECLARED_BEFORE_SCORES.flag`, which does not match the pattern.

## 5. `results/genespace/run_report.txt` gene-list counts

| list | entries | unique | duplicate symbol |
|---|---|---|---|
| age_up | 1533 | 1532 | GOLGA8M |
| age_down | 2007 | 2006 | TMEM191A |

- unique_up + unique_down = 1532 + 2006 = 3538. Confirmed (the report's checks line also gives `'age_unique': 3538`).
- entries_up + entries_down = 1533 + 2007 = 3540. Confirmed.
- The checks line reports `'up_down_overlap_lists': 0`.

## `*_ci_valid` columns equal to False

### `same_stats.csv` (key: space, test)

| quantity | n False | rows |
|---|---|---|
| cos | 3 | FULL Pluri; FULL NonReprog; AGE Pluri |
| progress | 0 | — |
| delta | 4 | FULL Pluri; FULL NonReprog; AGE Pluri; AGE NonReprog |
| perp | 4 | FULL Pluri; FULL NonReprog; AGE Pluri; AGE NonReprog |
| rel_delta | 4 | FULL Pluri; FULL NonReprog; AGE Pluri; AGE NonReprog |

`asymmetry_ci_valid` (only filled for forward rows): 0 False. No forward or reverse row has any False `*_ci_valid`.

### `same_posctrl.csv` (key: space, f)

| quantity | n False | rows |
|---|---|---|
| cos | 7 | FULL 0.00, 0.10, 0.25, 0.50; AGE 0.10, 0.25, 0.50 |
| progress | 0 | — |
| delta | 8 | FULL 0.00, 0.10, 0.25, 0.50; AGE 0.00, 0.10, 0.25, 0.50 |
| perp | 10 | FULL 0.00, 0.10, 0.25, 0.50; MD 0.00, 0.50; AGE 0.00, 0.10, 0.25, 0.50 |
| rel_delta | 8 | FULL 0.00, 0.10, 0.25, 0.50; AGE 0.00, 0.10, 0.25, 0.50 |

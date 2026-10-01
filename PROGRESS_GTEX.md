# PROGRESS_GTEX

**STOP status:** Stage 1 complete.
**Next action:** none

## Stage gate

- `FINDINGS_POSCTRL.md` exists: True
- status: `S2_i`  open=True
- source: `s2_c1_dense_shared_r030.json+s2_gene_ridge.json`
- reason: gate open: S2 outcome (i) — planted and real-age r transfer both ways (gene_ridge r HBCC→MSSM +0.568 / MSSM→HBCC +0.578; uncalibrated R² is a calibration artifact)
- primary transfer metric: held-out Pearson r, n_perm=200, pass = r>0 with null ≤ 0.05
- bootstrap angle reported only (sex ref 41.3°, planted n=233 41.191°)
- pre-registration 2026-09-17 written before any T-cell: True

## Pre-registered bars (verbatim, 2026-09-17)

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

## Pre-registered reading (verbatim, r not R²)

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

## Seeds

- main: `20260914`
- T-cell RNGs: {'T1': 20260915, 'T2': 20260916, 'T3': 20260917, 'T4': 20260918, 'T5': 20260919, 'T6': 20260920, 'COMP': 20260921}

## Frozen preprocessing spec (Stage 0)

- Tissues (SMTSD), six only, `SMAFRZE == RNASEQ`
- Gene filter: ≥6 counts in ≥20% of samples per tissue
- TMM (edgeR defaults) → log2-CPM `prior.count=2`
- z-score per gene on training rows only inside every fold (Stage 1)
- Map to DLPFC gene universe, Ensembl version-stripped

## Frozen marker lists

- Neuronal: RBFOX3, SNAP25, SYT1, GAD1, SLC17A7
- Glial: GFAP, AQP4, MBP, PLP1, OLIG2, CX3CR1, AIF1

## Stage 0 STOP checks

- BA9 n=269  t2_not_testable=True  n_SMCENTER_ge25=2
- BA9 SMCENTER: `{'B1, A1': 149, 'C1, A1': 113, 'D1, A1': 4, 'C1, B1, A1': 3}`

## Finished Stage 1 files

- `results/gtex/t1.json`
- `results/gtex/t2.json`
- `results/gtex/t3.json`
- `results/gtex/t4.json`
- `results/gtex/t5.json`
- `results/gtex/t6.json`
- `results/gtex/t1_comp.json`
- `results/gtex/stage1_summary.json`

## Finished cell numbers (from disk)

### T1
- Brain - Frontal Cortex (BA9) | raw | ridge | angle=+46.129° | r=+0.508 (null -0.006) | R²=+0.268
- Brain - Frontal Cortex (BA9) | raw | pls1 | angle=+38.223° | r=+0.192 (null -0.016) | R²=+0.003
- Brain - Frontal Cortex (BA9) | resid | ridge | angle=+49.107° | r=+0.327 (null -0.032) | R²=+0.102
- Brain - Frontal Cortex (BA9) | resid | pls1 | angle=+50.089° | r=+0.143 (null -0.025) | R²=-0.141
- Brain - Cortex | raw | ridge | angle=+46.702° | r=+0.421 (null -0.027) | R²=+0.184
- Brain - Cortex | raw | pls1 | angle=+30.387° | r=+0.119 (null -0.012) | R²=-0.103
- Brain - Cortex | resid | ridge | angle=+48.072° | r=+0.382 (null -0.015) | R²=+0.141
- Brain - Cortex | resid | pls1 | angle=+47.630° | r=+0.273 (null -0.017) | R²=+0.005
- Brain - Anterior cingulate cortex (BA24) | raw | ridge | angle=+45.822° | r=+0.522 (null -0.020) | R²=+0.278
- Brain - Anterior cingulate cortex (BA24) | raw | pls1 | angle=+36.452° | r=+0.279 (null -0.014) | R²=+0.037
- Brain - Anterior cingulate cortex (BA24) | resid | ridge | angle=+48.311° | r=+0.397 (null -0.035) | R²=+0.154
- Brain - Anterior cingulate cortex (BA24) | resid | pls1 | angle=+45.784° | r=+0.261 (null -0.014) | R²=+0.009
- Brain - Hippocampus | raw | ridge | angle=+48.857° | r=+0.454 (null -0.016) | R²=+0.199
- Brain - Hippocampus | raw | pls1 | angle=+34.597° | r=+0.281 (null -0.018) | R²=+0.047
- Brain - Hippocampus | resid | ridge | angle=+50.492° | r=+0.285 (null -0.029) | R²=+0.057
- Brain - Hippocampus | resid | pls1 | angle=+44.244° | r=+0.250 (null -0.008) | R²=-0.014
- Heart - Left Ventricle | raw | ridge | angle=+47.753° | r=+0.450 (null -0.004) | R²=+0.200
- Heart - Left Ventricle | raw | pls1 | angle=+24.742° | r=+0.239 (null -0.003) | R²=+0.036
- Heart - Left Ventricle | resid | ridge | angle=+47.670° | r=+0.220 (null -0.019) | R²=+0.041
- Heart - Left Ventricle | resid | pls1 | angle=+54.728° | r=+0.048 (null -0.006) | R²=-0.109
- Muscle - Skeletal | raw | ridge | angle=+45.519° | r=+0.708 (null -0.010) | R²=+0.503
- Muscle - Skeletal | raw | pls1 | angle=+18.031° | r=+0.407 (null -0.008) | R²=+0.162
- Muscle - Skeletal | resid | ridge | angle=+43.445° | r=+0.531 (null -0.001) | R²=+0.278
- Muscle - Skeletal | resid | pls1 | angle=+29.615° | r=+0.288 (null -0.003) | R²=+0.072

### T2
- SKIPPED: t2_not_testable: fewer than 3 SMCENTER with ≥25 BA9 samples

### T3 / T4
- T3 Brain - Frontal Cortex (BA9) -> Brain - Cortex raw ridge r=+0.480 (null +0.020) R²=+0.202 calR²=+0.230 ρ=+0.441 pass=True
- T3 Brain - Frontal Cortex (BA9) -> Brain - Cortex raw pls1 r=+0.244 (null +0.102) R²=+0.040 calR²=+0.059 ρ=+0.251 pass=False
- T3 Brain - Frontal Cortex (BA9) -> Brain - Cortex resid ridge r=+0.377 (null +0.031) R²=+0.097 calR²=+0.142 ρ=+0.380 pass=True
- T3 Brain - Frontal Cortex (BA9) -> Brain - Cortex resid pls1 r=+0.226 (null +0.048) R²=-0.005 calR²=+0.051 ρ=+0.249 pass=True
- T3 Brain - Cortex -> Brain - Frontal Cortex (BA9) raw ridge r=+0.407 (null -0.003) R²=-0.157 calR²=+0.165 ρ=+0.318 pass=True
- T3 Brain - Cortex -> Brain - Frontal Cortex (BA9) raw pls1 r=+0.251 (null +0.074) R²=-0.181 calR²=+0.063 ρ=+0.227 pass=False
- T3 Brain - Cortex -> Brain - Frontal Cortex (BA9) resid ridge r=+0.235 (null +0.035) R²=-0.330 calR²=+0.055 ρ=+0.203 pass=True
- T3 Brain - Cortex -> Brain - Frontal Cortex (BA9) resid pls1 r=+0.187 (null +0.022) R²=-0.458 calR²=+0.035 ρ=+0.131 pass=True
- T4 Brain - Frontal Cortex (BA9) -> Brain - Anterior cingulate cortex (BA24) raw ridge r=+0.360 (null -0.045) R²=+0.092 calR²=+0.129 ρ=+0.319 pass=True
- T4 Brain - Frontal Cortex (BA9) -> Brain - Anterior cingulate cortex (BA24) raw pls1 r=+0.154 (null +0.007) R²=-0.004 calR²=+0.024 ρ=+0.147 pass=True
- T4 Brain - Frontal Cortex (BA9) -> Brain - Anterior cingulate cortex (BA24) resid ridge r=+0.226 (null +0.024) R²=+0.048 calR²=+0.051 ρ=+0.195 pass=True
- T4 Brain - Frontal Cortex (BA9) -> Brain - Anterior cingulate cortex (BA24) resid pls1 r=+0.045 (null +0.007) R²=-0.244 calR²=+0.002 ρ=+0.004 pass=True
- T4 Brain - Frontal Cortex (BA9) -> Brain - Hippocampus raw ridge r=+0.334 (null -0.010) R²=+0.097 calR²=+0.111 ρ=+0.286 pass=True
- T4 Brain - Frontal Cortex (BA9) -> Brain - Hippocampus raw pls1 r=+0.251 (null +0.057) R²=-0.018 calR²=+0.063 ρ=+0.214 pass=False
- T4 Brain - Frontal Cortex (BA9) -> Brain - Hippocampus resid ridge r=+0.021 (null -0.013) R²=-0.653 calR²=+0.000 ρ=+0.069 pass=True
- T4 Brain - Frontal Cortex (BA9) -> Brain - Hippocampus resid pls1 r=-0.103 (null -0.025) R²=-0.673 calR²=+0.011 ρ=-0.071 pass=False

### T5
- {'n_common_genes': 20945, 'identity_rank': 5, 'identity_sv': [196.076038519334, 103.21330954875826, 62.547256494037455, 37.517359074661265, 22.58801759230702], 'ba9_angle_to_identity': 62.97933952920945, 'ba9_proj_r2': 0.20639917631617533, 'null_angle_mean': 70.35209728027063, 'null_angle_median': 69.87495134545492, 'p_angle': 0.10945273631840796, 'pairwise_median_deg': 68.64961381794926, 'pairwise_null_mean': 84.89764363242153, 'p_pairwise': 0.004975124378109453, 'tissues': ['Brain - Frontal Cortex (BA9)', 'Brain - Cortex', 'Brain - Anterior cingulate cortex (BA24)', 'Brain - Hippocampus', 'Heart - Left Ventricle', 'Muscle - Skeletal'], 'n_perm': 200, 'note': 'bulk analogue of FINDINGS_GEOMETRY Part C, not a replication'}

### T6
- {'ran': True, 'n_overlap': 20115, 'angle_deg': 80.19970901762603, 'null_mean': 89.67042162311502, 'null_median': 89.7175332213159, 'p': 0.000999000999000999, 'source': 'v1_axes.npz', 'w_key': 'W_age_ridge', 'n_types_collapsed': 20, 'n_random': 1000, 'note': 'one number with a null; not evidence beyond closer-than-random'}

## Note

FINDINGS_GTEX.md written from T1-T6 + COMP

## Files

- `src/gtex_common.py`, `gtex_stage0.py`, `gtex_stage1.py`, `gtex_findings.py`
- `notebooks/gtex_stage0.ipynb`, `gtex_stage1.ipynb`
- `FINDINGS_GTEX.md`
- `PROGRESS_GTEX.md`


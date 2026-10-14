# FINDINGS_DEPTH2 — the DEPTH test with references at the experiment's own depth

**Status: DONE.** Readings: n_total=400: `USABLE`, n_total=1200: `USABLE`.

Reproduced by `python src/depth2_run.py`; Stage 0 comparison, tables and readings re-checked from tracked files by `python src/depth2_run.py verify`. DEPTH's readings (400 `INSUFFICIENT`, 1,200 `MARGINAL`, FINDINGS_DEPTH.md) stand on record and are not revised here.

## Pre-registration (verbatim from `results/depth2/PREREG.flag`)

```
TASK: DEPTH2 — the DEPTH test with references at the experiment's own depth.
Written 2026-10-14T05:02Z (2026-10-14 14:02 +0900), before any DEPTH2 number
exists. No threshold or wording in this block is changed after numbers exist.
The only stored numbers read so far are prior results: DEPTH's flag, findings
and stored tables (results/depth/PREREG.flag, FINDINGS_DEPTH.md,
detection.csv, readings.json; commit 5255f87) and the sizes of the SAME split
index arrays in results/same/split_idx.npz (O 2,469, Y 3,851, aged_B 2,470,
young_B 3,852), which DEPTH's flag already states.

--- WHY THIS TASK EXISTS ---
DEPTH (FINDINGS_DEPTH.md, commit 5255f87) subsampled only the sample S and
kept the anchors O and Y at half-A depth (2,469 and 3,851 cells). Its
diagnostic section shows the failures at 400 cells came from the delta < 0
clause: a shallow S sits farther from Y than a deep O through sampling noise.
In the planned pilot every sample, including the aged 0% reference, is
sequenced at the same depth. DEPTH2 asks the pilot's question: does the
mixing control work when O, Y and S are all at n cells? DEPTH's readings
stand on record and are not revised by this task.

--- FIXED UNCHANGED FROM DEPTH (read from results/depth/PREREG.flag and
src/depth_run.py; reused, not modified) ---
- Per-cell counts data/processed/md4/louvain_rulerY_<line>.npz and md2
  labels, loaded by toward_run.load_donor; cell index arrays read from
  results/same/split_idx.npz (O = GM00731 half A, Y = GM23815 half A,
  aged_B = GM00731 half B, young_B = GM23815 half B).
- Frozen ruler gene space and frozen GTEx mu/sd (results/fibro/
  frozen_ruler_ridge_raw.npz) via fibro2_common.load_frozen_ruler. Not refit.
- Transform: pseudobulk (sum of cell counts) -> TMM log2-CPM, prior.count=2,
  across the panel -> frozen z; genes with zero counts in every panel row
  z=0. Called as same_run.panel_z.
- Statistics: same_run.axis_stats(z_O, z_Y, z_S) -> frac = dot(d,v)/||v||^2,
  delta = ||z_S - z_Y|| - ||z_O - z_Y||, with v = z_Y - z_O, d = z_S - z_O.
- f grid 0, 0.10, 0.25, 0.50 (f_index 0..3); n_young = round(f * n_total)
  (numpy round), n_aged = n_total - n_young.
  400:  aged/young = 400/0, 360/40, 300/100, 200/200
  1200: aged/young = 1200/0, 1080/120, 900/300, 600/600
- 200 draws per (n_total, f).
- Detection rule: per depth, floor_p5 = numpy.percentile(that depth's 200
  f = 0 deltas, 5) (numpy default linear interpolation). A draw is DETECTED
  iff delta < floor_p5 AND delta < 0, both strict. A non-finite delta is not
  detected. Detection rate = n_detected / 200, compared as >= 0.80.
- Chance baseline beside every rate: the same rule applied to the f = 0
  draws themselves (in-sample; at most 5% by construction).
- Monotonic: for each (n_total, draw_index), frac at f = 0, 0.10, 0.25, 0.50
  passed to same_run.monotonic_increasing (non-decreasing; any non-finite ->
  not monotonic). Monotonic share = share of the 200 draw indices.
- Readings, nesting and wording exactly as DEPTH (below).

--- THE ONLY CHANGE ---
In every draw, O, Y and S are all built at n_total cells, from disjoint cell
pools:
- O: n_total cells drawn without replacement from GM00731 half A (2,469).
- Y: n_total cells drawn without replacement from GM23815 half A (3,851).
- S: the mixture, drawn from the half-B pools exactly as DEPTH did: aged
  cells = rng.choice(aged_B, n_aged, replace=False) (skipped when n_aged = 0),
  then young cells = rng.choice(young_B, n_young, replace=False) (skipped
  when n_young = 0), both from the S generator.
O and Y are redrawn in every draw, so the f = 0 floor includes reference
noise. v = z(Y) - z(O) is recomputed per draw.
Seeds: each role in draw (n_index, f_index, draw_index) has its own generator
  rng = numpy.random.default_rng([20260915, n_index, f_index, draw_index, role])
  role 0 = O, role 1 = Y, role 2 = S.
  n_index 0 = 400, n_index 1 = 1200 (DEPTH's mapping kept); f_index 0..3 as
  above; draw_index 0..199.
Depths: 400 and 1,200, each read independently. No 2,470 rung.

Decisions the task text did not fix, fixed here before computing:
- Scoring panel. Because the seed key includes f_index for every role, O and
  Y are independent per (n_total, f, draw); there is no shared O, Y to put
  into P's 6-row panel [O, Y, four mixes]. Each (n_total, f, draw) is
  therefore scored as its own 3-row panel [O, Y, S], transformed jointly by
  same_run.panel_z, and axis_stats(Z[0], Z[1], Z[2]). This differs from
  DEPTH's 6-row panel. The f = 0 floor and the spike-in draws are both under
  the same 3-row shape, so the rule is internally consistent; but a
  difference between DEPTH and DEPTH2 numbers is not attributable to anchor
  depth alone, since TMM's normalization context also changed from 6 rows
  to 3.
- Monotonicity under independent references. Within a draw index the four f
  values now have different, independently drawn O and Y, so monotonicity
  compares fracs measured against different references. That is the pilot's
  situation and is kept as the definition.
- Disjointness is checked, not assumed: O pool vs aged_B and Y pool vs
  young_B must have empty intersection (gated, under Stage 0).

Known properties of this design, stated before any number:
- At f = 0, O and S are both GM00731 samples of the same size from disjoint
  halves, so the f = 0 delta is no longer biased positive by an S-only depth
  mismatch; its sign is not forced either way. This is what the task is
  designed to test, and is not used to set any threshold.
- Half-A pools are reused across draws (400 or 1,200 of 2,469 aged; of
  3,851 young), so draws overlap in cells.

--- STAGE 0 (before anything new) ---
With this task's code, run in DEPTH's configuration (deep anchors):
O = all 2,469 half-A aged cells, Y = all 3,851 half-A young cells, the same
in every draw; S cells from depth_run.draw_cells (default_rng([20260914,
n_index, f_index, draw_index]), n_index 0 = 400, 1 = 1200); one 6-row panel
[O, Y, mix f=0, 0.10, 0.25, 0.50] per (n_total, draw_index), as DEPTH. This
task's detection-table code is applied to those draws and compared with
results/depth/detection.csv rows n_total = 400 and 1,200:
  counts exact: n_aged, n_young, n_draws, n_detected; role string exact;
  numerics within 1e-4 * max(1, |stored|): every other column
  (detection_rate, floor_p5_delta, median_frac, median_delta and all diag_
  columns).
Checks added beyond the task text, also gated:
  (a) split: recomputed d0 Fibroblast split (same_run.d0_fibroblast_idx,
      split_half, seed 20260914) equals results/same/split_idx.npz exactly;
  (b) disjointness of O vs aged_B and Y vs young_B;
  (c) per-draw frac and delta for every 400 / 1,200 row of
      results/depth/draws.csv within the same tolerance (1,600 rows);
  (d) results/depth/monotonic.csv n_monotonic exact for 400 and 1,200;
  (e) results/depth/readings.json rule_outcome exact for 400 and 1,200
      (INSUFFICIENT, MARGINAL).
Any gated mismatch: STOP `prior_results_not_reproduced`, naming every one.
Stage 0 output goes to results/depth2/stage0/; nothing is written into
results/depth/ or results/same/.

--- READINGS (per depth; 400 and 1200 each get their own; wording verbatim
from DEPTH, the only substitution is the depth N) ---
Let det10, det25 = detection rates at f = 0.10, 0.25; mono = monotonic share.
The four keys are mutually exclusive and exhaustive (DEPTH's nesting).

USABLE: det10 >= 0.80 AND det25 >= 0.80 AND mono >= 0.80.
  "At N cells per sample, this depth preserves the sensitivity the mixing
  control was validated at."

MARGINAL: det25 >= 0.80 AND det10 < 0.80.
  "At N cells per sample, the mixing control would detect a quarter-gap
  approach (a 25% young spike-in) but not a tenth-gap one (a 10% young
  spike-in). An effect smaller than a quarter of the donor gap would be
  undetectable by construction at this depth."

INSUFFICIENT: det25 < 0.80.
  "At N cells per sample, the mixing control does not work: even a 25% young
  spike-in is not detected in 80% of draws. A null result from an experiment
  at this depth could not be read."

NO_READING: det10 >= 0.80 AND det25 >= 0.80 AND mono < 0.80.
  "At N cells per sample, no pre-registered reading fired: both spike-in rungs
  were detected in at least 80% of draws, but frac was monotonic in f in fewer
  than 80% of draws. No conclusion about this depth's sensitivity is drawn."

A depth that fails here cannot be rescued by analysis. The reading is about what the experiment could see, not about biology.

Only the reading whose conditions fired is reported as the conclusion for
each depth. Everything else is a diagnostic, labelled as such.

--- DIAGNOSTICS (not gated, labelled) ---
- The clause breakdown DEPTH reported: per depth and f, share of draws below
  the floor, share with delta < 0.
- frac separation: the f = 0 draws' 95th percentile of frac, and the share
  of draws at each f with frac above it.
- Spread: 5th/95th percentiles of frac and delta, median cos, median v_norm
  (||v|| now varies by draw); per-draw cos, d_norm, v_norm, n_missing_z0 in
  draws.csv.
- DEPTH's stored detection rates printed beside DEPTH2's, descriptive only,
  with the caveat that panel shape also differs (above).

--- STOP KEYS ---
prior_results_not_reproduced: any gated Stage 0 mismatch (every one named).
pool_too_small: a pool cannot supply a requested draw (O: n_total > 2,469;
  Y: n_total > 3,851; S: n_aged > 2,470 or n_young > 3,852). Checked for the
  whole grid before any DEPTH2 draw.
After a STOP: write the findings with the STOP key and exit. No computation
after a STOP, no restart with an adjustment.

--- WHAT THIS DOES NOT COVER ---
- Two donors, one young (GM23815) and one aged (GM00731), from a
  reprogramming dataset (GSE297234), not the six donors of the planned pilot.
- Reads per cell are not modelled: cells are subsampled from a deeply
  sequenced run, and a CRO's recovered cells may be shallower per cell. The
  result is optimistic in that respect.
- The spike-in is genuinely young cells, the most favourable possible
  signal. A partial approach to young by an intervention is weaker.
- Nothing here about donor-to-donor variation. Each reference is a
  resample of one donor line, so the f = 0 floor contains cell-sampling
  noise only, not the between-donor spread a real aged reference would carry.
- Circularity, disclosed: the f = 0 floor and the spike-in draws come from
  the same pools and share cells across draws; the floor is the 5th
  percentile of the same draws the chance baseline is computed on.
- Day-0 fibroblasts of two cell lines only; nothing about other cell types,
  partially reprogrammed cells, or donors without an aged counterpart.

--- OUTPUTS ---
src/depth2_run.py; results/depth2/ (PREREG.flag, stage0/ (DEPTH-configuration
draws and tables), stage0_reproduction.csv, draws.csv, detection.csv,
monotonic.csv, readings.json, run.log); FINDINGS_DEPTH2.md with this block
verbatim at the top. One local commit, not pushed. results/transfer/
NOT_A_RESULT_gate_failed/ is not opened.
```

## Stage 0 — reproduction of DEPTH's stored numbers (DEPTH configuration, deep anchors)

This task's code, run with DEPTH's deep anchors, 6-row panels and seeds, compared with `results/depth/` for n_total = 400 and 1,200. Counts exact; numerics within 1e-4 x max(1, |stored|). All rows are gated; rows marked "extra check" were added beyond the task text. For the two per-draw rows, `stored` is the required number of failures (0), `reproduced` the number found, and `abs_diff` the largest absolute difference over the 1,600 draws. Full table: `results/depth2/stage0_reproduction.csv`.

| item | kind | gated | stored | reproduced | abs_diff | ok |
|---|---|---|---|---|---|---|
| split_idx.npz O indices identical (extra check) | exact | True | True | True |  | True |
| split_idx.npz aged_B indices identical (extra check) | exact | True | True | True |  | True |
| split_idx.npz Y indices identical (extra check) | exact | True | True | True |  | True |
| split_idx.npz young_B indices identical (extra check) | exact | True | True | True |  | True |
| O pool (GM00731 half A) and aged_B disjoint (extra check) | exact | True | 0 | 0 |  | True |
| Y pool (GM23815 half A) and young_B disjoint (extra check) | exact | True | 0 | 0 |  | True |
| detection n=400 f=0.00 n_aged | exact | True | 400 | 400 |  | True |
| detection n=400 f=0.00 n_young | exact | True | 0 | 0 |  | True |
| detection n=400 f=0.00 n_draws | exact | True | 200 | 200 |  | True |
| detection n=400 f=0.00 n_detected | exact | True | 0 | 0 |  | True |
| detection n=400 f=0.00 detection_rate | tol | True | 0.0 | 0.0 | 0.00e+00 | True |
| detection n=400 f=0.00 role | exact | True | chance baseline (f=0 against its own floor) | chance baseline (f=0 against its own floor) |  | True |
| detection n=400 f=0.00 floor_p5_delta | tol | True | 23.90194974342085 | 23.90194974342085 | 0.00e+00 | True |
| detection n=400 f=0.00 median_frac | tol | True | 0.0683745970437417 | 0.06837459704374177 | 6.94e-17 | True |
| detection n=400 f=0.00 median_delta | tol | True | 26.060753216995955 | 26.060753216995955 | 0.00e+00 | True |
| detection n=400 f=0.00 diag_p5_frac | tol | True | 0.0485516400443739 | 0.04855164004437392 | 2.08e-17 | True |
| detection n=400 f=0.00 diag_p95_frac | tol | True | 0.0874302437259712 | 0.08743024372597123 | 2.78e-17 | True |
| detection n=400 f=0.00 diag_p5_delta | tol | True | 23.90194974342085 | 23.90194974342085 | 0.00e+00 | True |
| detection n=400 f=0.00 diag_p95_delta | tol | True | 28.492105721942355 | 28.492105721942355 | 0.00e+00 | True |
| detection n=400 f=0.00 diag_median_cos | tol | True | 0.0848110119793038 | 0.0848110119793038 | 0.00e+00 | True |
| detection n=400 f=0.00 diag_share_delta_lt_0 | tol | True | 0.0 | 0.0 | 0.00e+00 | True |
| detection n=400 f=0.00 diag_share_below_floor | tol | True | 0.05 | 0.05 | 0.00e+00 | True |
| detection n=400 f=0.00 diag_p95_frac_f0 | tol | True | 0.0874302437259712 | 0.08743024372597123 | 2.78e-17 | True |
| detection n=400 f=0.00 diag_share_frac_gt_f0_p95 | tol | True | 0.05 | 0.05 | 0.00e+00 | True |
| detection n=400 f=0.10 n_aged | exact | True | 360 | 360 |  | True |
| detection n=400 f=0.10 n_young | exact | True | 40 | 40 |  | True |
| detection n=400 f=0.10 n_draws | exact | True | 200 | 200 |  | True |
| detection n=400 f=0.10 n_detected | exact | True | 0 | 0 |  | True |
| detection n=400 f=0.10 detection_rate | tol | True | 0.0 | 0.0 | 0.00e+00 | True |
| detection n=400 f=0.10 role | exact | True | spike-in | spike-in |  | True |
| detection n=400 f=0.10 floor_p5_delta | tol | True | 23.90194974342085 | 23.90194974342085 | 0.00e+00 | True |
| detection n=400 f=0.10 median_frac | tol | True | 0.1399461782342503 | 0.13994617823425032 | 2.78e-17 | True |
| detection n=400 f=0.10 median_delta | tol | True | 20.44850017103548 | 20.448500171035477 | 3.55e-15 | True |
| detection n=400 f=0.10 diag_p5_frac | tol | True | 0.1153966293678792 | 0.11539662936787923 | 2.78e-17 | True |
| detection n=400 f=0.10 diag_p95_frac | tol | True | 0.1648768424187326 | 0.16487684241873268 | 8.33e-17 | True |
| detection n=400 f=0.10 diag_p5_delta | tol | True | 17.83905160953385 | 17.83905160953385 | 0.00e+00 | True |
| detection n=400 f=0.10 diag_p95_delta | tol | True | 23.00707934224254 | 23.00707934224254 | 0.00e+00 | True |
| detection n=400 f=0.10 diag_median_cos | tol | True | 0.1708112984903896 | 0.17081129849038967 | 8.33e-17 | True |
| detection n=400 f=0.10 diag_share_delta_lt_0 | tol | True | 0.0 | 0.0 | 0.00e+00 | True |
| detection n=400 f=0.10 diag_share_below_floor | tol | True | 0.98 | 0.98 | 0.00e+00 | True |
| detection n=400 f=0.10 diag_p95_frac_f0 | tol | True | 0.0874302437259712 | 0.08743024372597123 | 2.78e-17 | True |
| detection n=400 f=0.10 diag_share_frac_gt_f0_p95 | tol | True | 1.0 | 1.0 | 0.00e+00 | True |
| detection n=400 f=0.25 n_aged | exact | True | 300 | 300 |  | True |
| detection n=400 f=0.25 n_young | exact | True | 100 | 100 |  | True |
| detection n=400 f=0.25 n_draws | exact | True | 200 | 200 |  | True |
| detection n=400 f=0.25 n_detected | exact | True | 0 | 0 |  | True |
| detection n=400 f=0.25 detection_rate | tol | True | 0.0 | 0.0 | 0.00e+00 | True |
| detection n=400 f=0.25 role | exact | True | spike-in | spike-in |  | True |
| detection n=400 f=0.25 floor_p5_delta | tol | True | 23.90194974342085 | 23.90194974342085 | 0.00e+00 | True |
| detection n=400 f=0.25 median_frac | tol | True | 0.239145862999553 | 0.23914586299955307 | 5.55e-17 | True |
| detection n=400 f=0.25 median_delta | tol | True | 13.469688319744598 | 13.469688319744598 | 0.00e+00 | True |
| detection n=400 f=0.25 diag_p5_frac | tol | True | 0.2163946365409283 | 0.21639463654092836 | 5.55e-17 | True |
| detection n=400 f=0.25 diag_p95_frac | tol | True | 0.2645696098557411 | 0.26456960985574113 | 5.55e-17 | True |
| detection n=400 f=0.25 diag_p5_delta | tol | True | 10.90635674772605 | 10.90635674772605 | 0.00e+00 | True |
| detection n=400 f=0.25 diag_p95_delta | tol | True | 15.99828525305128 | 15.998285253051279 | 1.78e-15 | True |
| detection n=400 f=0.25 diag_median_cos | tol | True | 0.2809507543255569 | 0.2809507543255569 | 0.00e+00 | True |
| detection n=400 f=0.25 diag_share_delta_lt_0 | tol | True | 0.0 | 0.0 | 0.00e+00 | True |
| detection n=400 f=0.25 diag_share_below_floor | tol | True | 1.0 | 1.0 | 0.00e+00 | True |
| detection n=400 f=0.25 diag_p95_frac_f0 | tol | True | 0.0874302437259712 | 0.08743024372597123 | 2.78e-17 | True |
| detection n=400 f=0.25 diag_share_frac_gt_f0_p95 | tol | True | 1.0 | 1.0 | 0.00e+00 | True |
| detection n=400 f=0.50 n_aged | exact | True | 200 | 200 |  | True |
| detection n=400 f=0.50 n_young | exact | True | 200 | 200 |  | True |
| detection n=400 f=0.50 n_draws | exact | True | 200 | 200 |  | True |
| detection n=400 f=0.50 n_detected | exact | True | 1 | 1 |  | True |
| detection n=400 f=0.50 detection_rate | tol | True | 0.005 | 0.005 | 0.00e+00 | True |
| detection n=400 f=0.50 role | exact | True | spike-in | spike-in |  | True |
| detection n=400 f=0.50 floor_p5_delta | tol | True | 23.90194974342085 | 23.90194974342085 | 0.00e+00 | True |
| detection n=400 f=0.50 median_frac | tol | True | 0.4161071358735982 | 0.41610713587359827 | 5.55e-17 | True |
| detection n=400 f=0.50 median_delta | tol | True | 3.886322991980684 | 3.886322991980684 | 0.00e+00 | True |
| detection n=400 f=0.50 diag_p5_frac | tol | True | 0.3840621179735534 | 0.3840621179735534 | 0.00e+00 | True |
| detection n=400 f=0.50 diag_p95_frac | tol | True | 0.4427490066713603 | 0.4427490066713603 | 0.00e+00 | True |
| detection n=400 f=0.50 diag_p5_delta | tol | True | 1.4411515075400034 | 1.4411515075400034 | 0.00e+00 | True |
| detection n=400 f=0.50 diag_p95_delta | tol | True | 6.98453864293967 | 6.98453864293967 | 0.00e+00 | True |
| detection n=400 f=0.50 diag_median_cos | tol | True | 0.4381556583403169 | 0.4381556583403169 | 0.00e+00 | True |
| detection n=400 f=0.50 diag_share_delta_lt_0 | tol | True | 0.005 | 0.005 | 0.00e+00 | True |
| detection n=400 f=0.50 diag_share_below_floor | tol | True | 1.0 | 1.0 | 0.00e+00 | True |
| detection n=400 f=0.50 diag_p95_frac_f0 | tol | True | 0.0874302437259712 | 0.08743024372597123 | 2.78e-17 | True |
| detection n=400 f=0.50 diag_share_frac_gt_f0_p95 | tol | True | 1.0 | 1.0 | 0.00e+00 | True |
| detection n=1200 f=0.00 n_aged | exact | True | 1200 | 1200 |  | True |
| detection n=1200 f=0.00 n_young | exact | True | 0 | 0 |  | True |
| detection n=1200 f=0.00 n_draws | exact | True | 200 | 200 |  | True |
| detection n=1200 f=0.00 n_detected | exact | True | 0 | 0 |  | True |
| detection n=1200 f=0.00 detection_rate | tol | True | 0.0 | 0.0 | 0.00e+00 | True |
| detection n=1200 f=0.00 role | exact | True | chance baseline (f=0 against its own floor) | chance baseline (f=0 against its own floor) |  | True |
| detection n=1200 f=0.00 floor_p5_delta | tol | True | 5.559254499106801 | 5.559254499106801 | 0.00e+00 | True |
| detection n=1200 f=0.00 median_frac | tol | True | 0.0945561021159389 | 0.09455610211593898 | 8.33e-17 | True |
| detection n=1200 f=0.00 median_delta | tol | True | 6.760856243241072 | 6.760856243241072 | 0.00e+00 | True |
| detection n=1200 f=0.00 diag_p5_frac | tol | True | 0.0848587978396841 | 0.08485879783968411 | 1.39e-17 | True |
| detection n=1200 f=0.00 diag_p95_frac | tol | True | 0.1050001960571969 | 0.10500019605719699 | 8.33e-17 | True |
| detection n=1200 f=0.00 diag_p5_delta | tol | True | 5.559254499106801 | 5.559254499106801 | 0.00e+00 | True |
| detection n=1200 f=0.00 diag_p95_delta | tol | True | 7.796826575631127 | 7.796826575631127 | 0.00e+00 | True |
| detection n=1200 f=0.00 diag_median_cos | tol | True | 0.1704303513775017 | 0.17043035137750173 | 2.78e-17 | True |
| detection n=1200 f=0.00 diag_share_delta_lt_0 | tol | True | 0.0 | 0.0 | 0.00e+00 | True |
| detection n=1200 f=0.00 diag_share_below_floor | tol | True | 0.05 | 0.05 | 0.00e+00 | True |
| detection n=1200 f=0.00 diag_p95_frac_f0 | tol | True | 0.1050001960571969 | 0.10500019605719699 | 8.33e-17 | True |
| detection n=1200 f=0.00 diag_share_frac_gt_f0_p95 | tol | True | 0.05 | 0.05 | 0.00e+00 | True |
| detection n=1200 f=0.10 n_aged | exact | True | 1080 | 1080 |  | True |
| detection n=1200 f=0.10 n_young | exact | True | 120 | 120 |  | True |
| detection n=1200 f=0.10 n_draws | exact | True | 200 | 200 |  | True |
| detection n=1200 f=0.10 n_detected | exact | True | 146 | 146 |  | True |
| detection n=1200 f=0.10 detection_rate | tol | True | 0.73 | 0.73 | 0.00e+00 | True |
| detection n=1200 f=0.10 role | exact | True | spike-in | spike-in |  | True |
| detection n=1200 f=0.10 floor_p5_delta | tol | True | 5.559254499106801 | 5.559254499106801 | 0.00e+00 | True |
| detection n=1200 f=0.10 median_frac | tol | True | 0.1650593455132629 | 0.16505934551326296 | 5.55e-17 | True |
| detection n=1200 f=0.10 median_delta | tol | True | -0.4620107591254907 | -0.4620107591254907 | 0.00e+00 | True |
| detection n=1200 f=0.10 diag_p5_frac | tol | True | 0.1560073005281476 | 0.15600730052814765 | 5.55e-17 | True |
| detection n=1200 f=0.10 diag_p95_frac | tol | True | 0.176152470734779 | 0.17615247073477905 | 5.55e-17 | True |
| detection n=1200 f=0.10 diag_p5_delta | tol | True | -1.7229765254069724 | -1.7229765254069724 | 0.00e+00 | True |
| detection n=1200 f=0.10 diag_p95_delta | tol | True | 0.6990867757828686 | 0.6990867757828686 | 0.00e+00 | True |
| detection n=1200 f=0.10 diag_median_cos | tol | True | 0.2904429139101483 | 0.2904429139101483 | 0.00e+00 | True |
| detection n=1200 f=0.10 diag_share_delta_lt_0 | tol | True | 0.73 | 0.73 | 0.00e+00 | True |
| detection n=1200 f=0.10 diag_share_below_floor | tol | True | 1.0 | 1.0 | 0.00e+00 | True |
| detection n=1200 f=0.10 diag_p95_frac_f0 | tol | True | 0.1050001960571969 | 0.10500019605719699 | 8.33e-17 | True |
| detection n=1200 f=0.10 diag_share_frac_gt_f0_p95 | tol | True | 1.0 | 1.0 | 0.00e+00 | True |
| detection n=1200 f=0.25 n_aged | exact | True | 900 | 900 |  | True |
| detection n=1200 f=0.25 n_young | exact | True | 300 | 300 |  | True |
| detection n=1200 f=0.25 n_draws | exact | True | 200 | 200 |  | True |
| detection n=1200 f=0.25 n_detected | exact | True | 200 | 200 |  | True |
| detection n=1200 f=0.25 detection_rate | tol | True | 1.0 | 1.0 | 0.00e+00 | True |
| detection n=1200 f=0.25 role | exact | True | spike-in | spike-in |  | True |
| detection n=1200 f=0.25 floor_p5_delta | tol | True | 5.559254499106801 | 5.559254499106801 | 0.00e+00 | True |
| detection n=1200 f=0.25 median_frac | tol | True | 0.2621620638447243 | 0.26216206384472435 | 5.55e-17 | True |
| detection n=1200 f=0.25 median_delta | tol | True | -9.848023888356408 | -9.848023888356408 | 0.00e+00 | True |
| detection n=1200 f=0.25 diag_p5_frac | tol | True | 0.2491775315797108 | 0.2491775315797108 | 0.00e+00 | True |
| detection n=1200 f=0.25 diag_p95_frac | tol | True | 0.2761345180354062 | 0.2761345180354062 | 0.00e+00 | True |
| detection n=1200 f=0.25 diag_p5_delta | tol | True | -11.315703829286685 | -11.315703829286685 | 0.00e+00 | True |
| detection n=1200 f=0.25 diag_p95_delta | tol | True | -8.417181266150676 | -8.417181266150676 | 0.00e+00 | True |
| detection n=1200 f=0.25 diag_median_cos | tol | True | 0.435350181033115 | 0.43535018103311507 | 5.55e-17 | True |
| detection n=1200 f=0.25 diag_share_delta_lt_0 | tol | True | 1.0 | 1.0 | 0.00e+00 | True |
| detection n=1200 f=0.25 diag_share_below_floor | tol | True | 1.0 | 1.0 | 0.00e+00 | True |
| detection n=1200 f=0.25 diag_p95_frac_f0 | tol | True | 0.1050001960571969 | 0.10500019605719699 | 8.33e-17 | True |
| detection n=1200 f=0.25 diag_share_frac_gt_f0_p95 | tol | True | 1.0 | 1.0 | 0.00e+00 | True |
| detection n=1200 f=0.50 n_aged | exact | True | 600 | 600 |  | True |
| detection n=1200 f=0.50 n_young | exact | True | 600 | 600 |  | True |
| detection n=1200 f=0.50 n_draws | exact | True | 200 | 200 |  | True |
| detection n=1200 f=0.50 n_detected | exact | True | 200 | 200 |  | True |
| detection n=1200 f=0.50 detection_rate | tol | True | 1.0 | 1.0 | 0.00e+00 | True |
| detection n=1200 f=0.50 role | exact | True | spike-in | spike-in |  | True |
| detection n=1200 f=0.50 floor_p5_delta | tol | True | 5.559254499106801 | 5.559254499106801 | 0.00e+00 | True |
| detection n=1200 f=0.50 median_frac | tol | True | 0.4303392249012456 | 0.43033922490124565 | 5.55e-17 | True |
| detection n=1200 f=0.50 median_delta | tol | True | -24.1968195195136 | -24.196819519513596 | 3.55e-15 | True |
| detection n=1200 f=0.50 diag_p5_frac | tol | True | 0.416044013380861 | 0.416044013380861 | 0.00e+00 | True |
| detection n=1200 f=0.50 diag_p95_frac | tol | True | 0.445558034853455 | 0.445558034853455 | 0.00e+00 | True |
| detection n=1200 f=0.50 diag_p5_delta | tol | True | -25.734150940036788 | -25.734150940036788 | 0.00e+00 | True |
| detection n=1200 f=0.50 diag_p95_delta | tol | True | -22.78703766454259 | -22.78703766454259 | 0.00e+00 | True |
| detection n=1200 f=0.50 diag_median_cos | tol | True | 0.6121822846404189 | 0.6121822846404189 | 0.00e+00 | True |
| detection n=1200 f=0.50 diag_share_delta_lt_0 | tol | True | 1.0 | 1.0 | 0.00e+00 | True |
| detection n=1200 f=0.50 diag_share_below_floor | tol | True | 1.0 | 1.0 | 0.00e+00 | True |
| detection n=1200 f=0.50 diag_p95_frac_f0 | tol | True | 0.1050001960571969 | 0.10500019605719699 | 8.33e-17 | True |
| detection n=1200 f=0.50 diag_share_frac_gt_f0_p95 | tol | True | 1.0 | 1.0 | 0.00e+00 | True |
| draws.csv 400/1200 rows matched (extra check) | exact | True | 1600 | 1600 |  | True |
| draws.csv per-draw n_aged all equal (extra check) | exact | True | True | True |  | True |
| draws.csv per-draw n_young all equal (extra check) | exact | True | True | True |  | True |
| draws.csv per-draw frac within tol, n_failed (extra check) | exact | True | 0 | 0 | 9.71e-17 | True |
| draws.csv per-draw delta within tol, n_failed (extra check) | exact | True | 0 | 0 | 3.55e-15 | True |
| monotonic n=400 n_monotonic (extra check) | exact | True | 200 | 200 |  | True |
| readings n=400 rule_outcome (extra check) | exact | True | INSUFFICIENT | INSUFFICIENT |  | True |
| monotonic n=1200 n_monotonic (extra check) | exact | True | 200 | 200 |  | True |
| readings n=1200 rule_outcome (extra check) | exact | True | MARGINAL | MARGINAL |  | True |

159 checks, 159 gated, 0 gated failures.

## Pre-registered readings (only the reading that fired, per depth)

### n_total = 400: `USABLE`

> At 400 cells per sample, this depth preserves the sensitivity the mixing control was validated at.
>
> A depth that fails here cannot be rescued by analysis. The reading is about what the experiment could see, not about biology.

Conditions: detection at f=0.10 = 86.5% (>= 80%: True); detection at f=0.25 = 100.0% (>= 80%: True); monotonic share = 99.0% (>= 80%: True). Chance baseline (f=0 draws against their own floor, in-sample, at most 5% by construction): 5.0%.

### n_total = 1,200: `USABLE`

> At 1,200 cells per sample, this depth preserves the sensitivity the mixing control was validated at.
>
> A depth that fails here cannot be rescued by analysis. The reading is about what the experiment could see, not about biology.

Conditions: detection at f=0.10 = 100.0% (>= 80%: True); detection at f=0.25 = 100.0% (>= 80%: True); monotonic share = 100.0% (>= 80%: True). Chance baseline (f=0 draws against their own floor, in-sample, at most 5% by construction): 5.0%.

## Detection table

Per depth and f, over 200 draws; O, Y and S all at n_total cells, each (n_total, f, draw) its own 3-row panel. `floor_p5_delta` is the 5th percentile of that depth's f=0 deltas. A draw is detected iff delta < floor_p5_delta and delta < 0. The f=0 row is the chance baseline: the same rule applied to the floor's own draws. File: `results/depth2/detection.csv`.

| n_total | f | n_aged | n_young | n_detected | detection_rate | floor_p5_delta | median_frac | median_delta |
|---|---|---|---|---|---|---|---|---|
| 400 | 0.00 | 400 | 0 | 10 | 5.0% | -2.749 | +0.208 | -0.015 |
| 400 | 0.10 | 360 | 40 | 173 | 86.5% | -2.749 | +0.244 | -4.541 |
| 400 | 0.25 | 300 | 100 | 200 | 100.0% | -2.749 | +0.302 | -10.601 |
| 400 | 0.50 | 200 | 200 | 200 | 100.0% | -2.749 | +0.404 | -19.828 |
| 1200 | 0.00 | 1200 | 0 | 10 | 5.0% | -1.424 | +0.139 | +0.094 |
| 1200 | 0.10 | 1080 | 120 | 200 | 100.0% | -1.424 | +0.194 | -6.429 |
| 1200 | 0.25 | 900 | 300 | 200 | 100.0% | -1.424 | +0.272 | -15.102 |
| 1200 | 0.50 | 600 | 600 | 200 | 100.0% | -1.424 | +0.409 | -28.227 |

Monotonicity: share of the 200 draw indices in which frac(f=0) <= frac(0.10) <= frac(0.25) <= frac(0.50). Within a draw index each f has its own independently drawn O, Y and S. File: `results/depth2/monotonic.csv`.

| n_total | n_draws | n_monotonic | monotonic_share |
|---|---|---|---|
| 400 | 200 | 198 | 99.0% |
| 1200 | 200 | 200 | 100.0% |

## Notes on the readings (not gated; added after the first run's numbers existed)

This section changes no reading. It records what the reading keys alone do not say.

- Margin at n_total=400, f=0.10: 173 of 200 detected, 160 needed for the 80% bar (binomial standard error of the rate at this value over 200 draws: 2.4 points).
- Margin at n_total=1,200, f=0.10: 200 of 200 detected, 160 needed for the 80% bar (rate at the boundary, so a binomial standard error is not informative; the rule-of-three 95% bound is 98.5%).
- `USABLE` is the pre-registered key for meeting the 80% bars. Its wording says the depth "preserves the sensitivity the mixing control was validated at"; it does not mean detection is unchanged. At 400 cells the 10% rung is detected in 86.5% of draws, against 100% in DEPTH's 2,470 reference rung (deep anchors).
- frac at f=0 is not zero: its median is +0.208 at 400, +0.139 at 1,200. O's sampling noise enters both d = z_S - z_O and v = z_Y - z_O, which adds a positive term to dot(d, v) and grows as depth falls; this is a property of the statistic, not tested separately here. frac is therefore only readable against the same depth's f=0 draws, never as an absolute number across depths.
- The panel-shape change (6 rows in DEPTH, 3 here) is confounded with the anchor-depth change. This task cannot say how much of the difference from DEPTH's 400 `INSUFFICIENT` is due to each.
- The f=0 floor contains cell-sampling noise of one aged line only. A pilot's aged reference from different donors would add between-donor spread, widening the floor; nothing here measures that.

## Diagnostics (not gated)

### Which clause of the detection rule binds

The same breakdown DEPTH reported: share of draws below the f=0 floor, share with delta < 0, and a frac-based separation that the rule does not use (share of draws whose frac exceeds the f=0 draws' 95th percentile). Not gated.

| n_total | f | detection_rate | diag_share_below_floor | diag_share_delta_lt_0 | diag_p95_frac_f0 | diag_share_frac_gt_f0_p95 |
|---|---|---|---|---|---|---|
| 400 | 0.00 | 5.0% | 5.0% | 50.5% | +0.224 | 5.0% |
| 400 | 0.10 | 86.5% | 86.5% | 99.5% | +0.224 | 97.5% |
| 400 | 0.25 | 100.0% | 100.0% | 100.0% | +0.224 | 100.0% |
| 400 | 0.50 | 100.0% | 100.0% | 100.0% | +0.224 | 100.0% |
| 1200 | 0.00 | 5.0% | 5.0% | 46.5% | +0.152 | 5.0% |
| 1200 | 0.10 | 100.0% | 100.0% | 100.0% | +0.152 | 100.0% |
| 1200 | 0.25 | 100.0% | 100.0% | 100.0% | +0.152 | 100.0% |
| 1200 | 0.50 | 100.0% | 100.0% | 100.0% | +0.152 | 100.0% |

- n_total=400: the f=0 floor's 5th-percentile delta is -2.749, and 50.5% of f=0 draws have delta < 0. At f=0.10, 86.5% of draws fall below the floor and 99.5% have delta < 0. At f=0.25 the shares are 100.0% and 100.0%. frac exceeds the f=0 draws' 95th percentile (+0.224) in 97.5% of f=0.10 draws and 100.0% of f=0.25 draws.
- n_total=1,200: the f=0 floor's 5th-percentile delta is -1.424, and 46.5% of f=0 draws have delta < 0. At f=0.10, 100.0% of draws fall below the floor and 100.0% have delta < 0. At f=0.25 the shares are 100.0% and 100.0%. frac exceeds the f=0 draws' 95th percentile (+0.152) in 100.0% of f=0.10 draws and 100.0% of f=0.25 draws.

### Spread

5th/95th percentiles of frac and delta over draws, median cos, median ||v|| (v = z(Y) - z(O) is recomputed per draw). Not gated.

| n_total | f | diag_p5_frac | diag_p95_frac | diag_p5_delta | diag_p95_delta | diag_median_cos | diag_median_v_norm |
|---|---|---|---|---|---|---|---|
| 400 | 0.00 | +0.186 | +0.224 | -2.749 | +2.880 | +0.323 | +136.388 |
| 400 | 0.10 | +0.227 | +0.265 | -7.695 | -1.975 | +0.375 | +136.009 |
| 400 | 0.25 | +0.280 | +0.321 | -13.163 | -7.797 | +0.449 | +135.788 |
| 400 | 0.50 | +0.378 | +0.423 | -22.871 | -16.700 | +0.553 | +135.224 |
| 1200 | 0.00 | +0.129 | +0.152 | -1.424 | +1.528 | +0.263 | +128.795 |
| 1200 | 0.10 | +0.182 | +0.206 | -8.078 | -4.678 | +0.360 | +128.582 |
| 1200 | 0.25 | +0.260 | +0.284 | -16.839 | -13.403 | +0.479 | +128.236 |
| 1200 | 0.50 | +0.396 | +0.424 | -30.022 | -26.453 | +0.628 | +127.756 |

### DEPTH beside DEPTH2 (descriptive)

DEPTH's stored numbers (deep anchors, one 6-row panel per draw index) beside DEPTH2's (anchors at n_total, one 3-row panel per draw). Two things changed between the columns: anchor depth, and the TMM panel shape. A difference is therefore not attributable to anchor depth alone. Descriptive only; DEPTH's readings are not revised.

| n_total | f | DEPTH_detection | DEPTH2_detection | DEPTH_floor_p5 | DEPTH2_floor_p5 | DEPTH_median_frac | DEPTH2_median_frac |
|---|---|---|---|---|---|---|---|
| 400 | 0.00 | 0.0% | 5.0% | +23.902 | -2.749 | +0.068 | +0.208 |
| 400 | 0.10 | 0.0% | 86.5% | +23.902 | -2.749 | +0.140 | +0.244 |
| 400 | 0.25 | 0.0% | 100.0% | +23.902 | -2.749 | +0.239 | +0.302 |
| 400 | 0.50 | 0.5% | 100.0% | +23.902 | -2.749 | +0.416 | +0.404 |
| 1200 | 0.00 | 0.0% | 5.0% | +5.559 | -1.424 | +0.095 | +0.139 |
| 1200 | 0.10 | 73.0% | 100.0% | +5.559 | -1.424 | +0.165 | +0.194 |
| 1200 | 0.25 | 100.0% | 100.0% | +5.559 | -1.424 | +0.262 | +0.272 |
| 1200 | 0.50 | 100.0% | 100.0% | +5.559 | -1.424 | +0.430 | +0.409 |

Monotonic share, DEPTH vs DEPTH2: 400: 100.0% vs 99.0%; 1,200: 100.0% vs 100.0%.

## What this does not cover

- Two donors, one young (GM23815) and one aged (GM00731), from a reprogramming dataset (GSE297234), not the six donors of the planned pilot.
- Reads per cell are not modelled. Cells are subsampled from a deeply sequenced run; a CRO's recovered cells may be shallower per cell. The result is optimistic in that respect.
- The spike-in is genuinely young cells, the most favourable possible signal. A partial approach to young by an intervention is weaker and less clean.
- Nothing here about donor-to-donor variation. Each reference is a resample of one donor line, so the f=0 floor contains cell-sampling noise only, not the between-donor spread a real aged reference would carry.
- Circularity, disclosed: the f=0 floor and the spike-in draws come from the same pools and share cells across draws, and the floor is the 5th percentile of the same draws the chance baseline is computed on.
- Day-0 fibroblasts of two cell lines only. Nothing here speaks to other cell types, to partially reprogrammed cells, or to donors without an aged counterpart.

## Deviations and decisions beyond the task text

- **Scoring panel** (not fixed by the task text; fixed in the flag before computing): the seed key includes f_index for O and Y, so the references are independent per (n_total, f, draw) and DEPTH's shared 6-row panel cannot be formed. Each draw is its own 3-row TMM panel [O, Y, S]. This changes TMM's normalization context relative to DEPTH, so DEPTH-vs-DEPTH2 differences are not attributable to anchor depth alone.
- **Monotonicity under independent references:** within a draw index, the four f values have different O and Y. Kept as the definition, stated in the flag.
- **n_index mapping** 0 = 400, 1 = 1200 (DEPTH's), stated in the flag; the task text gave the key but not the mapping.
- **Detection-table code** is this task's own function (`tables_from_draws(df, depths)`), because DEPTH's version iterates over its module-level 3-depth grid. It applies DEPTH's rule and produces DEPTH's columns; Stage 0 checks it against every column of DEPTH's stored table. One diagnostic column was added: `diag_median_v_norm`. Readings text, the rescue sentence, the floor percentile, the 0.80 bar, the tolerance and `n_split` are imported from `src/depth_run.py`, and the Stage 0 S draws use `depth_run.draw_cells`.
- **Extra Stage 0 checks** (all gated): split indices recomputed and equal to `results/same/split_idx.npz`; disjointness of the O pool vs aged_B and the Y pool vs young_B; all 1,600 per-draw frac and delta values of `results/depth/draws.csv` (400, 1,200); `monotonic.csv` counts; `readings.json` rule outcomes.
- **Stage 0 outputs** are written to `results/depth2/stage0/` (draws_deep.csv, detection_deep.csv, monotonic_deep.csv, readings_deep.json). Nothing was written into `results/depth/` or `results/same/`.
- **Import side effects:** importing `same_run` and `toward_run` (directly and via `depth_run`) runs `mkdir(exist_ok=True)` on `results/same/` and `results/toward/`. Both already existed; no file in them was written. Importing `depth_run` writes nothing.
- **Added diagnostics** beyond the task text: spread percentiles of frac and delta, median cos, median ||v||, and DEPTH's stored rates printed beside DEPTH2's. All declared in the flag, none gated.
- **Added after the readings existed:** the section "Notes on the readings". It was added after the first complete run had produced the readings. Then the whole script was re-run from scratch, so the committed code produces every committed file. The readings and every gated number of the re-run are identical to the first run.
- **Not opened:** `results/transfer/NOT_A_RESULT_gate_failed/`.

## Files

- `src/depth2_run.py`
- `results/depth2/PREREG.flag`
- `results/depth2/stage0_reproduction.csv`
- `results/depth2/stage0/` (DEPTH-configuration draws_deep.csv, detection_deep.csv, monotonic_deep.csv, readings_deep.json)
- `results/depth2/run.log`
- `results/depth2/draws.csv` (one row per draw: frac, delta, cos, d_norm, v_norm, ...)
- `results/depth2/detection.csv`
- `results/depth2/monotonic.csv`
- `results/depth2/readings.json`
- `FINDINGS_DEPTH2.md`


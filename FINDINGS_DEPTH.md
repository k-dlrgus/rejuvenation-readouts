# FINDINGS_DEPTH — cell depth for the mixing-series positive control

**Status: DONE.** Readings: n_total=400: `INSUFFICIENT`, n_total=1200: `MARGINAL`.

Reproduced by `python src/depth_run.py`; gates re-checked from tracked files by `python src/depth_run.py verify`.

## Pre-registration (verbatim from `results/depth/PREREG.flag`)

```
TASK: DEPTH — at what per-sample cell count does the mixing-series positive
control still detect a young spike-in?
Written 2026-10-14T04:22:14Z (2026-10-14 13:22 +0900), before any DEPTH number
exists. No threshold or wording in this block is changed after numbers exist.
The only stored numbers read so far are the SAME split sizes and P table
(results/same/split_table.csv, split.json, posctrl.csv, posctrl_summary.json),
which are prior results, not DEPTH numbers.

QUESTION. FINDINGS_SAME's positive control P mixed aged (GM00731) and young
(GM23815) half-B day-0 Fibroblast cells at young fractions f = 0, 0.10, 0.25,
0.50 with n_total = 2,470, and passed. A CRO quote offers ~400 or ~1,200
recovered cells per sample. Does P still fire at those depths?

--- INSTRUMENTS (read, not refit) ---
Everything P used, unchanged:
- Half-A anchors O (GM00731, 2,469 cells) and Y (GM23815, 3,851 cells) and the
  half-B pools (GM00731 2,470 cells, GM23815 3,852 cells), cell indices read
  from results/same/split_idx.npz.
- Per-cell counts data/processed/md4/louvain_rulerY_<line>.npz and md2 labels,
  loaded by toward_run.load_donor exactly as SAME did.
- Frozen ruler gene space and frozen GTEx mu/sd (results/fibro/
  frozen_ruler_ridge_raw.npz) via fibro2_common.load_frozen_ruler.
- Transform: pseudobulk (sum of cell counts) -> TMM log2-CPM, prior.count=2,
  across the panel -> frozen z; genes with zero counts in every panel row
  z=0. Called as same_run.panel_z.
- Statistics: same_run.axis_stats(z_O, z_Y, z_S) -> frac = dot(d,v)/||v||^2,
  delta = ||z_S - z_Y|| - ||z_O - z_Y||, cos (diagnostic), with v = z_Y - z_O,
  d = z_S - z_O.
No anchor, axis, basis, ruler or mu/sd is refit. src/same_run.py and
src/toward_run.py are imported, not modified.

--- STAGE 0 (before anything new) ---
With this task's code:
(a) Recompute each donor's d0 Fibroblast set and 50/50 split
    (same_run.d0_fibroblast_idx, same_run.split_half, seed 20260914). Half sizes
    must equal results/same/split_table.csv exactly (2,469 / 2,470 / 3,851 /
    3,852). The recomputed index arrays are also checked for exact equality
    with results/same/split_idx.npz (a check added beyond the task text).
(b) Re-run P itself: same_run.run_posctrl on the split_idx.npz indices, with
    same_run's module-level SAME_DIR redirected to results/depth/stage0/ so
    that nothing is written into results/same/. Compare to
    results/same/posctrl.csv and results/same/posctrl_summary.json:
      counts exact: f, n_aged, n_young, n_total at all four f;
      numerics within 1e-4 * max(1, |stored|): frac, delta, cos at all four f;
      pass flags exact: pass, monotonic, frac50_ci_excludes_0, delta50_neg.
    frac50_ci_excludes_0 is a bootstrap flag, so Stage 0 re-runs P's own
    bootstrap (B=200, seed 20260918) for this comparison only. No DEPTH
    reading uses a bootstrap. The CI bounds themselves are compared as a
    diagnostic and are not gated.
Any gated mismatch: STOP `prior_results_not_reproduced`, naming every one.

--- DESIGN ---
Depth grid n_total: 400 (n_index 0), 1200 (n_index 1), 2470 (n_index 2).
f grid: 0, 0.10, 0.25, 0.50 (f_index 0..3), by cell counts as P did:
  n_young = round(f * n_total) (numpy round), n_aged = n_total - n_young.
  400:  aged/young = 400/0, 360/40, 300/100, 200/200
  1200: aged/young = 1200/0, 1080/120, 900/300, 600/600
  2470: aged/young = 2470/0, 2223/247, 1852/618, 1235/1235
200 independent draws per (n_total, f). Draw (n_index, f_index, draw_index)
uses its own generator rng = numpy.random.default_rng([20260914, n_index,
f_index, draw_index]); aged cells = rng.choice(half-B aged pool, n_aged,
replace=False) (skipped when n_aged = 0), then young cells = rng.choice(half-B
young pool, n_young, replace=False) (skipped when n_young = 0). Aged and young
cells come from their own pools.
Cells are NOT shared across f within a draw index: each f has its own
generator and its own independent sample.
Scoring panel (a decision the task text did not fix): for each (n_total,
draw_index) one 6-row panel [O, Y, mix f=0, mix f=0.10, mix f=0.25,
mix f=0.50] is built and transformed jointly, exactly the panel shape P used.
Each mix row is scored as S in P's forward test against the O, Y rows of that
same panel. O and Y are the same half-A count vectors in every panel; their z
moves slightly between panels only through the joint TMM, as in P.
Per draw, report frac and delta (and cos, d_norm, v_norm as diagnostics).
Bootstrap CIs are not used; the spread over draws is the uncertainty.

Known properties of this design, stated before any number:
- At n_total = 2470, f = 0 needs all 2,470 cells of the aged half-B pool, so
  every f = 0 draw at that depth is the same cell set; its spread comes only
  from TMM coupling to the other panel rows, not from sampling. f = 0.10 at
  2470 takes 2,223 of 2,470 aged cells. The 2470 rung is therefore a
  reproduction / reference rung: its table is reported, but it gets no gated
  reading and its "floor" is not a sampling noise floor.
- Only S is subsampled. O and Y stay at half-A depth (2,469 and 3,851 cells).
  Sampling noise in a shallow S enters delta with a positive sign
  (||z_S - z_Y|| grows with noise). This is part of what is being measured and
  is not corrected.

--- DETECTION RULE ---
At each depth, the 200 f = 0 deltas are the noise floor. floor_p5 =
numpy.percentile(f=0 deltas, 5) (numpy default linear interpolation).
A draw at f > 0 is DETECTED iff delta < floor_p5 AND delta < 0 (both strict).
A non-finite delta is not detected. Detection rate = n_detected / 200.
Chance baseline, reported beside every rate: the same rule applied to the
f = 0 draws themselves (at most 5% by construction, since the floor is their
own 5th percentile; this self-application is in-sample).
Monotonic: for each (n_total, draw_index), frac at f = 0, 0.10, 0.25, 0.50
is passed to same_run.monotonic_increasing (non-decreasing; any non-finite
value -> not monotonic). Report the share of the 200 draw indices that are
monotonic.
Per depth and f report: detection rate, median frac, median delta, floor_p5.
All rates compare with >= 0.80 (i.e. >= 160 of 200).

--- READINGS (per depth; 400 and 1200 each get their own; wording verbatim,
the only substitution is the depth N) ---
The task text's three readings could overlap (10% rung and monotonicity pass,
25% rung fails: USABLE and INSUFFICIENT both true) or leave a gap (both rungs
pass, monotonicity fails: nothing true). Asked in chat before writing this
flag; the author chose (2026-10-14): USABLE also requires the 25% rung, and
the gap gets NO_READING. The four keys below are mutually exclusive and
exhaustive.

Let det10, det25 = detection rates at f = 0.10, 0.25; mono = monotonic share.

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
each depth. The 2470 rung's table, cos, d_norm, the f = 0.50 rung's rate and
any spread summaries are diagnostics, labelled as such.

--- STOP KEYS ---
prior_results_not_reproduced: any gated Stage 0 mismatch (every one named).
pool_too_small: a half-B pool cannot supply a requested draw (n_aged > 2,470
  or n_young > 3,852). Checked for the whole grid before any draw.
After a STOP: write the findings with the STOP key and exit. No computation
after a STOP, no restart with an adjustment.

--- WHAT THIS DOES NOT COVER ---
- Two donors, one young and one aged, from a reprogramming dataset
  (GSE297234: GM00731 aged, GM23815 young), not the six donors of the planned
  pilot. This measures depth sensitivity only; it says nothing about
  donor-to-donor variation, which more donors and not more cells would
  address.
- Cells here are subsampled from a deeply sequenced run. A CRO's 400
  recovered cells may also differ in reads per cell, which this cannot model.
  The result is therefore optimistic for the low-depth rungs.
- The spike-in is cells of a genuinely younger donor, which is the most
  favourable possible positive signal.
- The anchors O and Y stay at half-A depth (~2,500 and ~3,900 cells). In a
  real experiment at 400 or 1,200 cells per sample the references would also
  be shallow. This is a second reason the low-depth rungs are optimistic.
- The f = 0 floor and the spike-in draws come from the same half-B pools and
  overlap in cells across draws; the floor is estimated from the same cells
  the rule is applied to (disclosed, not refitted away).
- Day-0 fibroblasts only, from one aged and one young donor line; nothing here
  speaks to other cell types or to partially reprogrammed cells.

--- OUTPUTS ---
src/depth_run.py; results/depth/ (PREREG.flag, stage0/, stage0_reproduction.csv,
draws.csv, detection.csv, monotonic.csv, readings.json, run.log);
FINDINGS_DEPTH.md with this block verbatim at the top. One local commit, not
pushed.
```

## Stage 0 — reproduction of SAME's stored numbers

Counts exact; numerics within 1e-4 x max(1, |stored|). CI-bound rows are diagnostic (not gated). Full table: `results/depth/stage0_reproduction.csv`.

| item | kind | gated | stored | reproduced | abs_diff | ok |
|---|---|---|---|---|---|---|
| split_table GM00731 half A n_cells | exact | True | 2469 | 2469 |  | True |
| split_table GM00731 half B n_cells | exact | True | 2470 | 2470 |  | True |
| split_table GM23815 half A n_cells | exact | True | 3851 | 3851 |  | True |
| split_table GM23815 half B n_cells | exact | True | 3852 | 3852 |  | True |
| split_idx.npz O indices identical (extra check) | exact | True | True | True |  | True |
| split_idx.npz aged_B indices identical (extra check) | exact | True | True | True |  | True |
| split_idx.npz Y indices identical (extra check) | exact | True | True | True |  | True |
| split_idx.npz young_B indices identical (extra check) | exact | True | True | True |  | True |
| posctrl f=0.00 n_aged | exact | True | 2470 | 2470 |  | True |
| posctrl f=0.00 n_young | exact | True | 0 | 0 |  | True |
| posctrl f=0.00 n_total | exact | True | 2470 | 2470 |  | True |
| posctrl f=0.00 frac | tol | True | 0.1145874459296361 | 0.11458744592963613 | 2.78e-17 | True |
| posctrl f=0.00 delta | tol | True | 0.0680510968813763 | 0.06805109688137634 | 4.16e-17 | True |
| posctrl f=0.00 cos | tol | True | 0.2387908767452672 | 0.23879087674526725 | 5.55e-17 | True |
| posctrl f=0.00 frac_ci_lo (diagnostic) | tol | False | 0.0987166334829383 | 0.09871663348293837 | 6.94e-17 | True |
| posctrl f=0.00 frac_ci_hi (diagnostic) | tol | False | 0.1209518273591894 | 0.12095182735918947 | 6.94e-17 | True |
| posctrl f=0.00 delta_ci_lo (diagnostic) | tol | False | 5.8070464703982125 | 5.8070464703982125 | 0.00e+00 | True |
| posctrl f=0.00 delta_ci_hi (diagnostic) | tol | False | 8.334360668817 | 8.334360668817 | 0.00e+00 | True |
| posctrl f=0.10 n_aged | exact | True | 2223 | 2223 |  | True |
| posctrl f=0.10 n_young | exact | True | 247 | 247 |  | True |
| posctrl f=0.10 n_total | exact | True | 2470 | 2470 |  | True |
| posctrl f=0.10 frac | tol | True | 0.1824353709141252 | 0.18243537091412526 | 5.55e-17 | True |
| posctrl f=0.10 delta | tol | True | -7.870024259533409 | -7.870024259533409 | 0.00e+00 | True |
| posctrl f=0.10 cos | tol | True | 0.3707048047536838 | 0.37070480475368384 | 5.55e-17 | True |
| posctrl f=0.10 frac_ci_lo (diagnostic) | tol | False | 0.1673149980838796 | 0.16731499808387965 | 5.55e-17 | True |
| posctrl f=0.10 frac_ci_hi (diagnostic) | tol | False | 0.1906980196151144 | 0.19069801961511446 | 5.55e-17 | True |
| posctrl f=0.10 delta_ci_lo (diagnostic) | tol | False | -1.8997635252250715 | -1.8997635252250713 | 2.22e-16 | True |
| posctrl f=0.10 delta_ci_hi (diagnostic) | tol | False | 1.0069515019350677 | 1.0069515019350677 | 0.00e+00 | True |
| posctrl f=0.25 n_aged | exact | True | 1852 | 1852 |  | True |
| posctrl f=0.25 n_young | exact | True | 618 | 618 |  | True |
| posctrl f=0.25 n_total | exact | True | 2470 | 2470 |  | True |
| posctrl f=0.25 frac | tol | True | 0.2829286266887166 | 0.2829286266887166 | 0.00e+00 | True |
| posctrl f=0.25 delta | tol | True | -18.91286402632018 | -18.91286402632018 | 0.00e+00 | True |
| posctrl f=0.25 cos | tol | True | 0.5303698817642954 | 0.5303698817642954 | 0.00e+00 | True |
| posctrl f=0.25 frac_ci_lo (diagnostic) | tol | False | 0.2671343628713917 | 0.26713436287139175 | 5.55e-17 | True |
| posctrl f=0.25 frac_ci_hi (diagnostic) | tol | False | 0.2913917403292013 | 0.29139174032920134 | 5.55e-17 | True |
| posctrl f=0.25 delta_ci_lo (diagnostic) | tol | False | -11.793133013896414 | -11.793133013896414 | 0.00e+00 | True |
| posctrl f=0.25 delta_ci_hi (diagnostic) | tol | False | -9.372419797939743 | -9.372419797939745 | 1.78e-15 | True |
| posctrl f=0.50 n_aged | exact | True | 1235 | 1235 |  | True |
| posctrl f=0.50 n_young | exact | True | 1235 | 1235 |  | True |
| posctrl f=0.50 n_total | exact | True | 2470 | 2470 |  | True |
| posctrl f=0.50 frac | tol | True | 0.4408188614368522 | 0.4408188614368522 | 0.00e+00 | True |
| posctrl f=0.50 delta | tol | True | -34.443006230983826 | -34.443006230983826 | 0.00e+00 | True |
| posctrl f=0.50 cos | tol | True | 0.6935059545570107 | 0.6935059545570107 | 0.00e+00 | True |
| posctrl f=0.50 frac_ci_lo (diagnostic) | tol | False | 0.424343381063624 | 0.42434338106362407 | 5.55e-17 | True |
| posctrl f=0.50 frac_ci_hi (diagnostic) | tol | False | 0.4474586408649893 | 0.4474586408649893 | 0.00e+00 | True |
| posctrl f=0.50 delta_ci_lo (diagnostic) | tol | False | -25.403148658470336 | -25.403148658470336 | 0.00e+00 | True |
| posctrl f=0.50 delta_ci_hi (diagnostic) | tol | False | -22.41415003731634 | -22.41415003731634 | 0.00e+00 | True |
| posctrl_summary pass | exact | True | True | True |  | True |
| posctrl_summary monotonic | exact | True | True | True |  | True |
| posctrl_summary frac50_ci_excludes_0 | exact | True | True | True |  | True |
| posctrl_summary delta50_neg | exact | True | True | True |  | True |
| posctrl_summary n_total | exact | True | 2470 | 2470 |  | True |

53 checks, 37 gated, 0 gated failures, 0 diagnostic failures.

## Pre-registered readings (only the reading that fired, per depth)

### n_total = 400: `INSUFFICIENT`

> At 400 cells per sample, the mixing control does not work: even a 25% young spike-in is not detected in 80% of draws. A null result from an experiment at this depth could not be read.
>
> A depth that fails here cannot be rescued by analysis. The reading is about what the experiment could see, not about biology.

Conditions: detection at f=0.10 = 0.0% (>= 80%: False); detection at f=0.25 = 0.0% (>= 80%: False); monotonic share = 100.0% (>= 80%: True). Chance baseline (f=0 draws against their own floor): 0.0%.

### n_total = 1,200: `MARGINAL`

> At 1,200 cells per sample, the mixing control would detect a quarter-gap approach (a 25% young spike-in) but not a tenth-gap one (a 10% young spike-in). An effect smaller than a quarter of the donor gap would be undetectable by construction at this depth.
>
> A depth that fails here cannot be rescued by analysis. The reading is about what the experiment could see, not about biology.

Conditions: detection at f=0.10 = 73.0% (>= 80%: False); detection at f=0.25 = 100.0% (>= 80%: True); monotonic share = 100.0% (>= 80%: True). Chance baseline (f=0 draws against their own floor): 0.0%.

## Detection table

Per depth and f, over 200 draws. `floor_p5_delta` is the 5th percentile of that depth's f=0 deltas. A draw at f>0 is detected iff delta < floor_p5_delta and delta < 0. The f=0 row is the chance baseline: the same rule applied to the floor's own draws (in-sample, at most 5% by construction). The 2470 rung is a reference rung with no gated reading; its f=0 draws are all the same 2,470 cells, so its floor is not a sampling floor. File: `results/depth/detection.csv`.

| n_total | f | n_aged | n_young | n_detected | detection_rate | floor_p5_delta | median_frac | median_delta |
|---|---|---|---|---|---|---|---|---|
| 400 | 0.00 | 400 | 0 | 0 | 0.0% | +23.902 | +0.068 | +26.061 |
| 400 | 0.10 | 360 | 40 | 0 | 0.0% | +23.902 | +0.140 | +20.449 |
| 400 | 0.25 | 300 | 100 | 0 | 0.0% | +23.902 | +0.239 | +13.470 |
| 400 | 0.50 | 200 | 200 | 1 | 0.5% | +23.902 | +0.416 | +3.886 |
| 1200 | 0.00 | 1200 | 0 | 0 | 0.0% | +5.559 | +0.095 | +6.761 |
| 1200 | 0.10 | 1080 | 120 | 146 | 73.0% | +5.559 | +0.165 | -0.462 |
| 1200 | 0.25 | 900 | 300 | 200 | 100.0% | +5.559 | +0.262 | -9.848 |
| 1200 | 0.50 | 600 | 600 | 200 | 100.0% | +5.559 | +0.430 | -24.197 |
| 2470 | 0.00 | 2470 | 0 | 0 | 0.0% | +0.044 | +0.115 | +0.064 |
| 2470 | 0.10 | 2223 | 247 | 200 | 100.0% | +0.044 | +0.184 | -8.001 |
| 2470 | 0.25 | 1852 | 618 | 200 | 100.0% | +0.044 | +0.279 | -18.387 |
| 2470 | 0.50 | 1235 | 1235 | 200 | 100.0% | +0.044 | +0.441 | -34.557 |

Monotonicity: share of the 200 draw indices in which frac(f=0) <= frac(0.10) <= frac(0.25) <= frac(0.50) (`same_run.monotonic_increasing`). Cells are not shared across f within a draw index. File: `results/depth/monotonic.csv`.

| n_total | n_draws | n_monotonic | monotonic_share |
|---|---|---|---|
| 400 | 200 | 200 | 100.0% |
| 1200 | 200 | 200 | 100.0% |
| 2470 | 200 | 200 | 100.0% |

## Diagnostics (not gated)

- 2470 reference rung: the rule, if applied, would give `USABLE` (det10=100.0%, det25=100.0%, mono=100.0%). Not a reading: the f=0 floor at this depth has no sampling spread (every draw is the full aged pool), so the floor is narrower than a real replicate would be.
- Spread of frac and delta over draws, median cos and the share of draws with delta < 0, per depth and f:

| n_total | f | diag_p5_frac | diag_p95_frac | diag_p5_delta | diag_p95_delta | diag_median_cos | diag_share_delta_lt_0 |
|---|---|---|---|---|---|---|---|
| 400 | 0.00 | +0.049 | +0.087 | +23.902 | +28.492 | +0.085 | +0.000 |
| 400 | 0.10 | +0.115 | +0.165 | +17.839 | +23.007 | +0.171 | +0.000 |
| 400 | 0.25 | +0.216 | +0.265 | +10.906 | +15.998 | +0.281 | +0.000 |
| 400 | 0.50 | +0.384 | +0.443 | +1.441 | +6.985 | +0.438 | +0.005 |
| 1200 | 0.00 | +0.085 | +0.105 | +5.559 | +7.797 | +0.170 | +0.000 |
| 1200 | 0.10 | +0.156 | +0.176 | -1.723 | +0.699 | +0.290 | +0.730 |
| 1200 | 0.25 | +0.249 | +0.276 | -11.316 | -8.417 | +0.435 | +1.000 |
| 1200 | 0.50 | +0.416 | +0.446 | -25.734 | -22.787 | +0.612 | +1.000 |
| 2470 | 0.00 | +0.115 | +0.115 | +0.044 | +0.083 | +0.239 | +0.000 |
| 2470 | 0.10 | +0.179 | +0.189 | -8.525 | -7.393 | +0.373 | +1.000 |
| 2470 | 0.25 | +0.272 | +0.286 | -19.281 | -17.541 | +0.524 | +1.000 |
| 2470 | 0.50 | +0.432 | +0.450 | -35.633 | -33.537 | +0.695 | +1.000 |

## Which clause of the detection rule binds (diagnostic; added after the readings were computed)

This section was added after the first run's readings existed. It changes no reading and is not a gate. It splits the detection rule into its two clauses: below the f=0 floor, and delta < 0. It also gives a frac-based separation that the rule does not use.

| n_total | f | detection_rate | diag_share_below_floor | diag_share_delta_lt_0 | diag_p95_frac_f0 | diag_share_frac_gt_f0_p95 |
|---|---|---|---|---|---|---|
| 400 | 0.10 | 0.0% | 98.0% | 0.0% | +0.087 | 100.0% |
| 400 | 0.25 | 0.0% | 100.0% | 0.0% | +0.087 | 100.0% |
| 400 | 0.50 | 0.5% | 100.0% | 0.5% | +0.087 | 100.0% |
| 1200 | 0.10 | 73.0% | 100.0% | 73.0% | +0.105 | 100.0% |
| 1200 | 0.25 | 100.0% | 100.0% | 100.0% | +0.105 | 100.0% |
| 1200 | 0.50 | 100.0% | 100.0% | 100.0% | +0.105 | 100.0% |
| 2470 | 0.10 | 100.0% | 100.0% | 100.0% | +0.115 | 100.0% |
| 2470 | 0.25 | 100.0% | 100.0% | 100.0% | +0.115 | 100.0% |
| 2470 | 0.50 | 100.0% | 100.0% | 100.0% | +0.115 | 100.0% |

- n_total=400: the f=0 floor's 5th-percentile delta is +23.902, and 0.0% of f=0 draws have delta < 0. At f=0.10, 98.0% of draws fall below the floor and 0.0% have delta < 0. At f=0.25 the shares are 100.0% and 0.0%. frac exceeds the f=0 draws' 95th percentile (+0.087) in 100.0% of f=0.10 draws.
- n_total=1,200: the f=0 floor's 5th-percentile delta is +5.559, and 0.0% of f=0 draws have delta < 0. At f=0.10, 100.0% of draws fall below the floor and 73.0% have delta < 0. At f=0.25 the shares are 100.0% and 100.0%. frac exceeds the f=0 draws' 95th percentile (+0.105) in 100.0% of f=0.10 draws.

Plainly: wherever the floor clause holds but detection fails, the failure is the delta < 0 clause. A shallow S is farther from Y than the deep half-A origin O is, through sampling noise alone. The spike-in moves S toward Y, but not past O's distance. frac is a projection, so noise does not add to it the way it adds to a distance. But frac is not depth-invariant either: its f=0 median is +0.068 at 400, +0.095 at 1,200, +0.115 at 2,470. This describes how the pre-registered rule behaves. It is not a re-reading. Swapping to a frac-based rule now would be the after-the-fact rescue the flag forbids. A future flag that wanted a depth-robust rule would have to fix it in advance, and would have to model anchors at the experiment's own depth, which this design keeps deep.

## What this does not cover

- Two donors, one young (GM23815) and one aged (GM00731), from a reprogramming dataset, not the six donors of the planned pilot. This measures depth sensitivity only; it says nothing about donor-to-donor variation, which more donors and not more cells would address.
- Cells here are subsampled from a deeply sequenced run. A CRO's 400 recovered cells may also differ in reads per cell, which this cannot model. The result is therefore optimistic for the low-depth rungs.
- The spike-in is cells of a genuinely younger donor, which is the most favourable possible positive signal. A partial approach to young by an intervention is a weaker and less clean signal than mixing in real young cells.
- Only the sample S is subsampled. The anchors O and Y stay at half-A depth (2,469 and 3,851 cells). In a real experiment at 400 or 1,200 cells per sample the references would also be shallow, so this is a second reason the low-depth rungs are optimistic.
- Circularity, disclosed: the f=0 floor and the spike-in draws come from the same half-B pools and share cells across draws, and the floor is the 5th percentile of the same draws the chance baseline is computed on.
- Day-0 fibroblasts of two cell lines only. Nothing here speaks to other cell types, to partially reprogrammed cells, or to donors without an aged counterpart.

## Deviations and decisions beyond the task text

- **Reading overlap and gap.** The task text's USABLE and INSUFFICIENT could both hold (10% rung and monotonicity pass, 25% rung fails), and no reading held when both rungs passed but monotonicity failed. Asked in chat before the flag was written. The author chose to make USABLE also require the 25% rung and to add `NO_READING` for the gap. The USABLE condition is therefore stricter than the task text.
- **Panel composition** (not fixed by the task text): each (depth, draw index) is one 6-row TMM panel [O, Y, four mixes], the shape P used. TMM coupling means the O and Y z-scores move slightly between panels, as in P.
- **Cells across f:** not shared. Each (depth, f, draw) has its own generator.
- **Stage 0 re-runs P's bootstrap** (B=200, seed 20260918), to reproduce the stored `frac50_ci_excludes_0` flag. No DEPTH reading uses a bootstrap.
- **Output redirect:** `same_run.SAME_DIR` was set to `results/depth/stage0/` for the Stage 0 call of `same_run.run_posctrl`, which writes `posctrl.csv`, `posctrl_summary.json` and `posctrl_z.npz` there. Nothing was written into `results/same/`.
- **Import side effects:** importing `same_run` and `toward_run` runs `mkdir(exist_ok=True)` on `results/same/` and `results/toward/`. Both already existed, so nothing was created and no file in them was written.
- **Extra Stage 0 check:** the recomputed split index arrays were compared for exact equality with `results/same/split_idx.npz`, beyond the half sizes.
- **Fixed in the flag, not in the task text:** numpy-default (linear) percentile for the floor; strict inequalities; a non-finite delta counts as not detected; monotonic means non-decreasing (P's `monotonic_increasing`); rates compared as >= 0.80; the 2470 rung gets no gated reading.
- **Added diagnostics:** cos, d_norm, v_norm, n_missing_z0 per draw; 5th/95th percentiles of frac and delta; share of draws with delta < 0. None of these are gated.
- **Added after the readings existed:** the columns `diag_share_below_floor`, `diag_p95_frac_f0` and `diag_share_frac_gt_f0_p95`, and the section on which clause binds. They were added after the first complete run had produced the readings. Then the whole script was re-run from scratch, so the committed code produces every committed file. The readings and every gated number of the re-run are identical to the first run.
- **Not opened:** `results/transfer/NOT_A_RESULT_gate_failed/`.

## Files

- `src/depth_run.py`
- `results/depth/PREREG.flag`
- `results/depth/stage0_reproduction.csv`
- `results/depth/stage0/` (P re-run: posctrl.csv, posctrl_summary.json, posctrl_z.npz)
- `results/depth/run.log`
- `results/depth/draws.csv` (one row per draw: frac, delta, cos, ...)
- `results/depth/detection.csv`
- `results/depth/monotonic.csv`
- `results/depth/readings.json`
- `FINDINGS_DEPTH.md`


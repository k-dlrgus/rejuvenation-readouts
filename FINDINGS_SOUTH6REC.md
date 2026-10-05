# FINDINGS_SOUTH6REC — recovery control for the SOUTH6 readout, depth grid

Status: completed. Stage 0 reproduced SOUTH6; Stage 1 and the report-only readings ran.

## Pre-registration (verbatim, `results/south6rec/PREREG.flag`)

```
TASK: SOUTH6REC — recovery control for the SOUTH6 screen readout, on a depth grid.
Written 2026-10-04 09:45 +0900, before any SOUTH6REC number exists.

Nature of this task: computational control on data SOUTH6 already analyzed.
No lab work, no new data, no download. Nothing in results/south6/ is modified;
every file there is read as stored. The frozen ruler and the TOWARD anchors are
read, never refit. src/south6_run.py and every other existing src/ file are
imported unchanged.

QUESTION. SOUTH6 found toward_young = 0 of 1,836 factors in every primary
(O1 -> Y1) space. Could the readout have recovered a real approach to young at
the depth the screen had? Plant a displacement of known size along the
old-to-young direction into non-targeting (NT) pseudobulks of fixed size and
measure whether the change in distance recovers it.

--- DATA ---
The SOUTH6 file as analyzed: Southard Hs27 per-cell singlets,
data/raw/southard/fibroblast_CRISPRa_final_pop_singlets_normalized_log1p.h5ad,
md5 2c926c7cae2eaab1ed2519b06ca21db5, 8,048,055,098 bytes. The md5 is
confirmed; STOP `file_mismatch` if it differs. Counts come from SOUTH6's
matched-gene CSR cache (results/south6/matched_data.f64, matched_indices.i32,
matched_indptr.npy), which holds layers/counts restricted to the 20,264 matched
columns. obs is read from the h5ad for guide labels only.
NT cells: guide_target == 'non', the SOUTH6 definition (14,013 cells, 78
guides). STOP `nt_cells_mismatch` if the obs-derived NT cell set differs from
results/south6/stage1.npz nt_cells.

--- SPACES, ORIGIN, TARGET ---
Four spaces, the SOUTH6 primary settings, origin O1 = GTEx 60-79 centroid
(c_old), target Y1 = GTEx 20-39 centroid (c_young), from
results/toward/anchors.npz:
  S1       gene space: the 23,485 ruler genes, unmatched genes at z = 0.
  S2_k20   S2_k50   S2_k100: projection onto the first k PCs of the SOUTH6
           basis, rebuilt exactly as src/south6_run.py build_spaces does
           (lowdim_common.fit_pca on the 652 GTEx donors' frozen-z matrix,
           centered on the donor mean). No other basis enters the gate.
Gap g = c_young - c_old (S1 coordinates). r = O1 - Y1 in each space (projected
for S2). delta = ||O1 + d - Y1|| - ||O1 - Y1|| in that space, with d projected
for S2. This is block_delta from src/south6_run.py.

--- STAGE 0: REPRODUCE SOUTH6 FIRST (gate, before any new number) ---
Rebuild with this task's code, from the CSR cache and obs, and compare with the
stored files:
  (a) per-factor counts C_f (1,836 x 20,264) from stage1.npz cell lists, per-NT-
      guide counts C_nt (78 rows, guides mapped from obs guide_identity), and
      C_pool: against stage1.npz. Exact (integer sums), tol 0.
  (b) TMM / log2 CPM (prior 2) / frozen z on the rebuilt panel, d_f = z_f -
      z_pool: against stage1.npz d_f.
  (c) ||O1 - Y1|| in S1, S2_k20, S2_k50, S2_k100: against
      origin_target_distances.csv.
  (d) per factor, S1 and each S2_k, O1_Y1: delta, cos, frac, d_norm against
      stage3_full_table.csv.
  (e) Stage 2B random-split floor and split median, all78 control set, O1_Y1,
      four spaces: 200 half-splits of the 14,013 NT cells with
      default_rng(20260914), as split_effects_of does, against
      stage2b_control_floor.csv.
  (f) Stage 2C dose-response deltas, O1_Y1, four spaces, against
      stage2c_dose_response.csv.
Pass: every numeric item within |new - stored| <= 1e-4 * max(1, |stored|)
(the SOUTH6 tolerance COORDSPACE used), counts exact. Fail -> STOP
`south6_not_reproduced`, naming every mismatch; no new number is computed.
Not recomputed: Stage 2A permutation floors, p and q. This task computes no
permutation; COORDSPACE reproduced those on 2026-10-03.

--- STAGE 1: RECOVERY CONTROL (the gated experiment) ---
Depths n = 106 (the observed median cells per factor, stage0_cells_per_target),
390, 2,500 cells. Doses f = 0, 0.25, 0.50, 1.00 of the gap. R = 200 draws per
(depth, dose).
Plant. A real perturbation can only show on genes the panel measures. The
planted vector is f * g_m, where g_m is g with every coordinate the SOUTH6
pipeline holds at z = 0 set to 0 (ruler genes outside the 20,264 matched
columns, and the matched columns with zero counts across the SOUTH6 panel).
The plant is added in z space after the transform: d = e + f * g_m, where e is
the null displacement of the draw. ||g_m|| / ||g|| is reported per space.
Null displacement e (PRIMARY): e = z(draw) - z(complement), where the draw is n
NT cells and the complement is the other 14,013 - n NT cells, so the draw's
cells are never in its own reference, as a real factor's cells are never in the
NT pool. z() is z_of from src/south6_run.py with the TMM cache of the stored
SOUTH6 panel.
Draw scheme, per depth, stream default_rng([20260914, n]): repeatedly permute
the 14,013 NT cells and cut the permutation into m = floor(14013 / n)
consecutive disjoint blocks of n; assign blocks round-robin to doses in the
order 0, 0.25, 0.50, 1.00, using 4 * floor(m / 4) blocks per permutation;
stop when every dose holds 200 draws. All draws from one permutation are
mutually disjoint, across doses too. Draws from different permutations
overlap at random.
Readout per draw: delta in each of the four spaces.
Per (space, depth, dose) report: median delta with a 95% bootstrap CI (2,000
resamples of the draws, seed 20260918, one generator per condition), mean, SD,
5th and 95th percentiles, and the noise-free delta (e = 0).
Monotone: median delta strictly decreasing across f = 0 < 0.25 < 0.50 < 1.00.
Separable at dose f: power(f) >= 0.80, where power(f) is the fraction of the
200 dose-f draws whose delta is below the 5th percentile (numpy default linear
interpolation) of the 200 dose-0 draws' delta, same space and depth. Power gets
a 95% bootstrap CI that resamples both the dose-0 and the dose-f draws (2,000
resamples, seed 20260918); the CI is reported, not gated on.
Smallest separable effect: the smallest f in {0.25, 0.50, 1.00} that is
separable; "none" if none is.

--- GATE (declared now, per space and depth) ---
A space QUALIFIES at depth n if median delta falls monotonically across the four
doses AND the 25% dose is separable from 0%. The existing-depth answer is
n = 106; 390 and 2,500 are read with the same rule, each on its own row.
Never pooled across spaces or depths.
Wording, fixed now:
  QUALIFIES: "<space>, <n> cells: recovered. Median delta falls monotonically
  with dose and a 25% plant is separated from the 0% draws (power <p>)."
  FAILS: "<space>, <n> cells: not recovered. <monotone: yes/no>; 25% plant
  power <p> against the 0.80 bar. At this depth and in this space, the readout
  cannot detect an approach of 25% of the gap."

--- OVERLAP AND VARIANCE, reported explicitly (never a gate) ---
Per depth: m, disjoint blocks used per permutation, permutations used, and how
many of each dose's 200 draws are mutually disjoint within the condition (the
largest set drawn from one permutation). Per (depth, dose): mean and maximum
pairwise cell overlap between draws in the condition, as a fraction of n,
against the n/14013 expected for independent draws.
State plainly: (1) every draw resamples the same 14,013 cells, so the replicate
spread is the spread of resampling this pool, not of new cells; (2) random
cells mix all 78 guides, so guide-to-guide variation, which a factor's ~6
guides carry, is absent from e.

--- REPORT-ONLY READINGS (no decision rule; never a gate) ---
R1 Guide-structured null. Per draw: guides in random order (stream
   default_rng([20260914, n, 1])), taken until at least 6 guides are chosen and
   their cells total >= n; n cells sampled from those guides' cells; reference
   = all cells of the unchosen guides. 200 draws per (depth, dose), independent
   draws. Same statistics, same separability rule, reported as report-only.
R2 Pool reference. e = z(draw) - z_pool, where z_pool is the stored SOUTH6
   all-NT reference and contains the draw's own cells (n/14013 of it). Same
   draws as the primary. Reported to show how much that overlap shrinks e.
R3 Full-gap plant. f * g including coordinates the panel cannot measure, as the
   stored Stage 2C did. Same draws as the primary.
R4 Mann-Whitney one-sided p and AUC (P(delta_dose < delta_0)) per condition.
R5 Basis circularity. The S2 basis is fit on all 652 GTEx donors, which include
   the 113 young and 231 old donors whose means are Y1 and O1. Report the gap
   retention ||proj_k g|| / ||g|| in the SOUTH6 basis and in a basis refit
   (same fit_pca, same centering rule) on the 308 donors aged 40-59 only, which
   excludes every anchor donor. Diagnostic only; the gate uses the SOUTH6
   basis.

--- WHAT DOES NOT COUNT ---
- any report-only reading as a pass or fail of the gate;
- a separability call from the bootstrap CI instead of the point power;
- choosing depth, dose, space, reference or plant after seeing results;
- pooling spaces or depths;
- reading recovery here as evidence about any factor's effect. A pass says
  the readout could see a planted effect of that size; it says nothing about
  whether any factor has one.

--- STOP RULES ---
On any STOP: write the findings with the key and exit. Never compute after a
STOP; never resume a stopped run.

--- OUTPUT ---
FINDINGS_SOUTH6REC.md: this text verbatim; the reconnaissance of the SOUTH6
space labels, origins, targets and basis; the Stage 0 table; Stage 1 per space
and depth; the gate reading that fired per row; overlap; report-only readings
labeled as such; limitations. Tables in results/south6rec/. Code
src/south6rec_run.py, run from the repository root.
```

## 1. Reconnaissance: what the SOUTH6 settings are (read-only)

> **The anchors entered the S2 basis.** The PCA behind S2_k20, S2_k50 and S2_k100 was fitted on all 652 GTEx fibroblast donors. That set includes the 113 donors aged 20-39 whose mean is Y1 and the 231 donors aged 60-79 whose mean is O1: 344 of the 652 rows (52.8%) are anchor donors. Only the 308 donors aged 40-59 are not. Gap retention in the S2 spaces is therefore partly circular. The basis was built from rows that contain the young-old difference it is later asked to retain. R5 below measures how much retention survives a basis that excludes every anchor donor.

Sources: `src/south6_run.py` `build_spaces` (lines 1293-1346) and `make_settings` (lines 1349-1402); `src/toward_run.py` `load_gtex_z` and `freeze_anchors`; `src/lowdim_common.py` `fit_pca`; `results/south6/south6_run.log`.

| label | what it is |
|---|---|
| S1 | Gene space. The 23,485 genes of the frozen ruler, in frozen-GTEx z units (log2 CPM minus the ruler's mu, over its sd). 20,265 ruler genes match the Southard panel (20,264 distinct source columns). The other 3,220, and 164 matched genes with zero counts across the panel, sit at z = 0 in every screen vector. |
| S2_k20, S2_k50, S2_k100 | Projection of S1 vectors onto the first k principal components of the GTEx donor matrix. Fit: `fit_pca(Zg - Zg.mean(0), 100)` on Zg = all 652 GTEx fibroblast donors x 23,485 genes in frozen z, all age bins (20-29: 60, 30-39: 53, 40-49: 106, 50-59: 202, 60-69: 210, 70-79: 21). Gram eigendecomposition, loadings scaled to unit norm. One fit at k = 100; k = 20 and 50 are its leading columns. No Southard cell enters the fit. Origins, targets and screen vectors are projected as (x - Zg.mean(0)) @ L[:, :k]. |
| O1 (primary origin) | `c_old` in `results/toward/anchors.npz`: the mean frozen-z profile of the 231 GTEx donors aged 60-79 (60-69: 210, 70-79: 21). Frozen, not refit. |
| Y1 (primary target) | `c_young` in the same file: the mean of the 113 donors aged 20-39 (20-29: 60, 30-39: 53). Frozen, not refit. |
| O2 (secondary origin) | The Lu aged donor GM00731 day-0 pseudobulk, row 'O' of `results/same/panel_z.npz`. |
| Y2 (secondary target) | A young region in S2 only. Its mean and covariance come from 56 of the 113 young donors (default_rng(20260914) permutation). Its radius is the 95th-percentile Mahalanobis distance of the other 57. Defined at k = 20 and 50; skipped at k = 100 (covariance rank 55 < 100). O1 already sits inside Y2 at both defined k (Mahalanobis 1.06 against radius 11.88 at k = 20; 6.01 against 59.30 at k = 50). |

Gap retention in the stored settings: ||O1 - Y1|| = 32.284 in S1, 28.645 in S2_k20 (88.7%), 30.864 in S2_k50 (95.6%), 31.434 in S2_k100 (97.4%) (`origin_target_distances.csv`).

Also shared with the anchors, and not part of the basis question: the frozen mu and sd that define z in every space, S1 included, come from the same 652 GTEx donors (`frozen_ruler_ridge_raw.npz`, note 'trained on all GTEx fibroblast donors').

The stored Stage 2C dose response, for contrast with Stage 1 here, is noise-free. It adds f x (median factor norm) along the unit gap direction to O1, with no sampling noise and including gap coordinates the panel cannot measure. Its MDA of 0.05 therefore says nothing about recovery from real pseudobulks.

## 2. Stage 0: SOUTH6 reproduced with this task's code (gate)

42 checks, 42 gating, 0 gating failures. Tolerance |new - stored| <= 1e-4 x max(1, |stored|); counts exact. The full table is `results/south6rec/stage0_reproduction.csv`.

| block | item | n | max_abs_diff | n_mismatch | passed |
|---|---|---|---|---|---|
| obs | NT cells (guide_target == 'non') equal stage1.npz nt_cells | 14013 | NA | 0 | True |
| (a) counts | C_f, 1836 factors x 20264 genes | 37204704 | 0.000 | 0 | True |
| (a) counts | C_nt, 78 guides x 20264 genes | 1580592 | 0.000 | 0 | True |
| (a) counts | C_pool | 20264 | 0.000 | 0 | True |
| (a) counts | C_pool equals sum of C_nt | 20264 | 0.000 | 0 | True |
| (b) transform | missing-gene mask | 20264 | NA | 0 | True |
| (b) transform | z_pool (matched genes) | 20264 | 0.000 | 0 | True |
| (b) transform | z_pool unmatched genes are 0 | 3221 | 0.000 | 0 | True |
| (b) transform | d_f, matched genes | 37204704 | 0.000 | 0 | True |
| (b) transform | d_f, unmatched genes are 0 | 5913756 | 0.000 | 0 | True |
| (c) spaces | ||O1 - Y1|| S1 | 1 | 0.000 | 0 | True |
| (c) spaces | ||O1 - Y1|| S2_k20 | 1 | 0.000 | 0 | True |
| (c) spaces | ||O1 - Y1|| S2_k50 | 1 | 0.000 | 0 | True |
| (c) spaces | ||O1 - Y1|| S2_k100 | 1 | 0.000 | 0 | True |
| (d) stage3 S1 O1_Y1 | delta | 1836 | 0.000 | 0 | True |
| (d) stage3 S1 O1_Y1 | cos | 1836 | 0.000 | 0 | True |
| (d) stage3 S1 O1_Y1 | frac | 1836 | 0.000 | 0 | True |
| (d) stage3 S1 O1_Y1 | d_norm | 1836 | 0.000 | 0 | True |
| (d) stage3 S2_k20 O1_Y1 | delta | 1836 | 0.000 | 0 | True |
| (d) stage3 S2_k20 O1_Y1 | cos | 1836 | 0.000 | 0 | True |
| (d) stage3 S2_k20 O1_Y1 | frac | 1836 | 0.000 | 0 | True |
| (d) stage3 S2_k20 O1_Y1 | d_norm | 1836 | 0.000 | 0 | True |
| (d) stage3 S2_k50 O1_Y1 | delta | 1836 | 0.000 | 0 | True |
| (d) stage3 S2_k50 O1_Y1 | cos | 1836 | 0.000 | 0 | True |
| (d) stage3 S2_k50 O1_Y1 | frac | 1836 | 0.000 | 0 | True |
| (d) stage3 S2_k50 O1_Y1 | d_norm | 1836 | 0.000 | 0 | True |
| (d) stage3 S2_k100 O1_Y1 | delta | 1836 | 0.000 | 0 | True |
| (d) stage3 S2_k100 O1_Y1 | cos | 1836 | 0.000 | 0 | True |
| (d) stage3 S2_k100 O1_Y1 | frac | 1836 | 0.000 | 0 | True |
| (d) stage3 S2_k100 O1_Y1 | d_norm | 1836 | 0.000 | 0 | True |
| (e) stage2b S1 O1_Y1 all78 | split_floor_p05 | 1 | 0.000 | 0 | True |
| (e) stage2b S1 O1_Y1 all78 | split_median | 1 | 0.000 | 0 | True |
| (e) stage2b S2_k20 O1_Y1 all78 | split_floor_p05 | 1 | 0.000 | 0 | True |
| (e) stage2b S2_k20 O1_Y1 all78 | split_median | 1 | 0.000 | 0 | True |
| (e) stage2b S2_k50 O1_Y1 all78 | split_floor_p05 | 1 | 0.000 | 0 | True |
| (e) stage2b S2_k50 O1_Y1 all78 | split_median | 1 | 0.000 | 0 | True |
| (e) stage2b S2_k100 O1_Y1 all78 | split_floor_p05 | 1 | 0.000 | 0 | True |
| (e) stage2b S2_k100 O1_Y1 all78 | split_median | 1 | 0.000 | 0 | True |
| (f) stage2c S1 O1_Y1 | delta at f = 0..1 | 6 | 0.000 | 0 | True |
| (f) stage2c S2_k20 O1_Y1 | delta at f = 0..1 | 6 | 0.000 | 0 | True |
| (f) stage2c S2_k50 O1_Y1 | delta at f = 0..1 | 6 | 0.000 | 0 | True |
| (f) stage2c S2_k100 O1_Y1 | delta at f = 0..1 | 6 | 0.000 | 0 | True |

## 3. Stage 1: the gate, per space and depth (primary null, measured-gene plant)

Plant f x g_m, where g_m is the young-minus-old gap restricted to the coordinates the panel measures. Its norm relative to the full gap in each space:

| space | ||g_m|| | ||g|| | ratio |
|---|---|---|---|
| S1 | 30.515 | 32.284 | 0.945 |
| S2_k20 | 26.019 | 28.645 | 0.908 |
| S2_k50 | 27.939 | 30.864 | 0.905 |
| S2_k100 | 28.362 | 31.434 | 0.902 |

| space | cells | monotone | power 25% | power 50% | power 100% | smallest separable | QUALIFIES |
|---|---|---|---|---|---|---|---|
| S1 | 106 | True | 0.115 | 0.340 | 0.470 | none | False |
| S1 | 390 | True | 0.955 | 1.000 | 1.000 | 25% | True |
| S1 | 2500 | True | 1.000 | 1.000 | 1.000 | 25% | True |
| S2_k20 | 106 | True | 0.990 | 1.000 | 1.000 | 25% | True |
| S2_k20 | 390 | True | 1.000 | 1.000 | 1.000 | 25% | True |
| S2_k20 | 2500 | True | 1.000 | 1.000 | 1.000 | 25% | True |
| S2_k50 | 106 | True | 0.965 | 1.000 | 1.000 | 25% | True |
| S2_k50 | 390 | True | 1.000 | 1.000 | 1.000 | 25% | True |
| S2_k50 | 2500 | True | 1.000 | 1.000 | 1.000 | 25% | True |
| S2_k100 | 106 | True | 0.930 | 1.000 | 1.000 | 25% | True |
| S2_k100 | 390 | True | 1.000 | 1.000 | 1.000 | 25% | True |
| S2_k100 | 2500 | True | 1.000 | 1.000 | 1.000 | 25% | True |

Power = fraction of the 200 dose draws whose delta is below the 5th percentile of the 200 dose-0 draws (bar 0.80).

### Readings that fired (pre-registered wording)

- S1, 106 cells: not recovered. Monotone: yes; 25% plant power 0.115 against the 0.80 bar. At this depth and in this space, the readout cannot detect an approach of 25% of the gap.
- S1, 390 cells: recovered. Median delta falls monotonically with dose and a 25% plant is separated from the 0% draws (power 0.955).
- S1, 2500 cells: recovered. Median delta falls monotonically with dose and a 25% plant is separated from the 0% draws (power 1.000).
- S2_k20, 106 cells: recovered. Median delta falls monotonically with dose and a 25% plant is separated from the 0% draws (power 0.990).
- S2_k20, 390 cells: recovered. Median delta falls monotonically with dose and a 25% plant is separated from the 0% draws (power 1.000).
- S2_k20, 2500 cells: recovered. Median delta falls monotonically with dose and a 25% plant is separated from the 0% draws (power 1.000).
- S2_k50, 106 cells: recovered. Median delta falls monotonically with dose and a 25% plant is separated from the 0% draws (power 0.965).
- S2_k50, 390 cells: recovered. Median delta falls monotonically with dose and a 25% plant is separated from the 0% draws (power 1.000).
- S2_k50, 2500 cells: recovered. Median delta falls monotonically with dose and a 25% plant is separated from the 0% draws (power 1.000).
- S2_k100, 106 cells: recovered. Median delta falls monotonically with dose and a 25% plant is separated from the 0% draws (power 0.930).
- S2_k100, 390 cells: recovered. Median delta falls monotonically with dose and a 25% plant is separated from the 0% draws (power 1.000).
- S2_k100, 2500 cells: recovered. Median delta falls monotonically with dose and a 25% plant is separated from the 0% draws (power 1.000).

## 4. Measured change in distance per dose (primary)

Median delta over 200 draws with its 95% bootstrap CI. Noise-free = delta with no sampling noise (e = 0). Units are each space's own; compare within a row block, not across spaces.

### S1

| cells | dose | median delta [95% CI] | mean | SD | 5th-95th pct | noise-free | power [95% CI] | MWU p (R4) | AUC (R4) |
|---|---|---|---|---|---|---|---|---|---|
| 106 | 0% | 93.603 [93.085, 94.103] | 93.711 | 2.735 | 89.244 to 98.235 | 0.000 | — | — | — |
| 106 | 25% | 92.152 [91.693, 92.642] | 92.157 | 2.658 | 88.208 to 96.250 | -7.087 | 0.115 [0.060, 0.230] | 2.72e-08 | 0.343 |
| 106 | 50% | 90.693 [89.963, 91.221] | 90.835 | 3.025 | 86.262 to 95.980 | -13.740 | 0.340 [0.225, 0.460] | 1.59e-19 | 0.241 |
| 106 | 100% | 89.535 [88.903, 89.981] | 89.668 | 2.757 | 85.600 to 94.331 | -21.743 | 0.470 [0.380, 0.605] | 5.05e-34 | 0.150 |
| 390 | 0% | 37.674 [37.503, 37.907] | 37.772 | 0.906 | 36.398 to 39.210 | 0.000 | — | — | — |
| 390 | 25% | 34.798 [34.680, 34.928] | 34.862 | 0.844 | 33.543 to 36.294 | -7.087 | 0.955 [0.905, 0.990] | 4.55e-65 | 0.009 |
| 390 | 50% | 32.710 [32.563, 32.832] | 32.682 | 0.788 | 31.490 to 33.921 | -13.740 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |
| 390 | 100% | 30.811 [30.674, 30.984] | 30.847 | 0.826 | 29.641 to 32.138 | -21.743 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |
| 2500 | 0% | 9.475 [9.440, 9.512] | 9.469 | 0.247 | 9.059 to 9.851 | 0.000 | — | — | — |
| 2500 | 25% | 4.242 [4.182, 4.280] | 4.234 | 0.245 | 3.802 to 4.606 | -7.087 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |
| 2500 | 50% | 0.014 [-0.038, 0.042] | 0.012 | 0.268 | -0.412 to 0.444 | -13.740 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |
| 2500 | 100% | -3.833 [-3.872, -3.803] | -3.838 | 0.245 | -4.231 to -3.439 | -21.743 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |

### S2_k20

| cells | dose | median delta [95% CI] | mean | SD | 5th-95th pct | noise-free | power [95% CI] | MWU p (R4) | AUC (R4) |
|---|---|---|---|---|---|---|---|---|---|
| 106 | 0% | 2.518 [2.264, 2.806] | 2.646 | 1.567 | 0.131 to 5.148 | 0.000 | — | — | — |
| 106 | 25% | -3.478 [-3.692, -3.142] | -3.407 | 1.560 | -5.656 to -0.809 | -6.502 | 0.990 [0.960, 1.000] | 1.72e-66 | 0.003 |
| 106 | 50% | -9.482 [-9.712, -9.238] | -9.377 | 1.530 | -11.729 to -6.481 | -13.001 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |
| 106 | 100% | -18.343 [-18.538, -18.035] | -18.175 | 1.623 | -20.427 to -15.122 | -25.915 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |
| 390 | 0% | 0.395 [0.314, 0.557] | 0.433 | 0.737 | -0.855 to 1.554 | 0.000 | — | — | — |
| 390 | 25% | -5.920 [-6.069, -5.843] | -5.917 | 0.710 | -7.198 to -4.843 | -6.502 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |
| 390 | 50% | -12.311 [-12.471, -12.136] | -12.262 | 0.770 | -13.563 to -10.938 | -13.001 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |
| 390 | 100% | -24.016 [-24.104, -23.953] | -23.963 | 0.801 | -25.140 to -22.737 | -25.915 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |
| 2500 | 0% | 0.063 [0.013, 0.133] | 0.098 | 0.322 | -0.377 to 0.666 | 0.000 | — | — | — |
| 2500 | 25% | -6.459 [-6.508, -6.395] | -6.447 | 0.333 | -7.009 to -5.906 | -6.502 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |
| 2500 | 50% | -12.911 [-12.975, -12.838] | -12.914 | 0.341 | -13.478 to -12.381 | -13.001 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |
| 2500 | 100% | -25.613 [-25.666, -25.567] | -25.608 | 0.298 | -26.047 to -25.103 | -25.915 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |

### S2_k50

| cells | dose | median delta [95% CI] | mean | SD | 5th-95th pct | noise-free | power [95% CI] | MWU p (R4) | AUC (R4) |
|---|---|---|---|---|---|---|---|---|---|
| 106 | 0% | 4.209 [3.818, 4.599] | 4.353 | 1.729 | 1.678 to 7.296 | 0.000 | — | — | — |
| 106 | 25% | -1.768 [-2.070, -1.390] | -1.634 | 1.706 | -4.260 to 1.171 | -6.974 | 0.965 [0.935, 0.990] | 3.81e-65 | 0.009 |
| 106 | 50% | -7.437 [-7.860, -6.911] | -7.317 | 1.859 | -10.021 to -3.751 | -13.939 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |
| 106 | 100% | -14.792 [-15.076, -14.559] | -14.597 | 2.074 | -17.775 to -10.606 | -27.609 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |
| 390 | 0% | 0.524 [0.399, 0.694] | 0.554 | 0.659 | -0.532 to 1.586 | 0.000 | — | — | — |
| 390 | 25% | -6.274 [-6.360, -6.162] | -6.221 | 0.663 | -7.309 to -5.142 | -6.974 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |
| 390 | 50% | -12.938 [-13.038, -12.803] | -12.874 | 0.750 | -14.013 to -11.513 | -13.939 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |
| 390 | 100% | -24.605 [-24.721, -24.474] | -24.524 | 0.803 | -25.709 to -23.208 | -27.609 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |
| 2500 | 0% | 0.098 [0.051, 0.152] | 0.106 | 0.279 | -0.364 to 0.573 | 0.000 | — | — | — |
| 2500 | 25% | -6.934 [-6.971, -6.886] | -6.898 | 0.281 | -7.334 to -6.389 | -6.974 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |
| 2500 | 50% | -13.819 [-13.883, -13.761] | -13.838 | 0.293 | -14.306 to -13.385 | -13.939 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |
| 2500 | 100% | -27.185 [-27.228, -27.143] | -27.161 | 0.252 | -27.529 to -26.687 | -27.609 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |

### S2_k100

| cells | dose | median delta [95% CI] | mean | SD | 5th-95th pct | noise-free | power [95% CI] | MWU p (R4) | AUC (R4) |
|---|---|---|---|---|---|---|---|---|---|
| 106 | 0% | 6.075 [5.677, 6.333] | 6.182 | 1.930 | 3.158 to 9.559 | 0.000 | — | — | — |
| 106 | 25% | 0.221 [0.011, 0.674] | 0.388 | 1.935 | -2.662 to 3.392 | -7.075 | 0.930 [0.890, 0.975] | 9.16e-63 | 0.018 |
| 106 | 50% | -5.025 [-5.521, -4.538] | -4.926 | 2.226 | -8.118 to -0.970 | -14.136 | 1.000 [1.000, 1.000] | 2.64e-67 | 0.000 |
| 106 | 100% | -11.745 [-11.970, -11.318] | -11.461 | 2.337 | -14.886 to -7.065 | -27.889 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |
| 390 | 0% | 0.765 [0.665, 0.939] | 0.802 | 0.664 | -0.251 to 1.824 | 0.000 | — | — | — |
| 390 | 25% | -6.025 [-6.181, -5.892] | -5.997 | 0.670 | -7.106 to -4.893 | -7.075 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |
| 390 | 50% | -12.724 [-12.837, -12.586] | -12.635 | 0.768 | -13.744 to -11.174 | -14.136 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |
| 390 | 100% | -23.821 [-23.994, -23.695] | -23.775 | 0.815 | -24.942 to -22.439 | -27.889 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |
| 2500 | 0% | 0.142 [0.077, 0.180] | 0.134 | 0.275 | -0.295 to 0.590 | 0.000 | — | — | — |
| 2500 | 25% | -6.984 [-7.032, -6.944] | -6.968 | 0.280 | -7.399 to -6.441 | -7.075 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |
| 2500 | 50% | -13.967 [-14.031, -13.920] | -13.987 | 0.287 | -14.443 to -13.530 | -14.136 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |
| 2500 | 100% | -27.307 [-27.344, -27.261] | -27.277 | 0.242 | -27.643 to -26.767 | -27.889 | 1.000 [1.000, 1.000] | 2.42e-67 | 0.000 |

## 5. Overlap and variance (stated explicitly)

Draw design per depth (14,013 NT cells):

| cells | disjoint blocks per permutation | blocks used | draws per dose per permutation | permutations |
|---|---|---|---|---|
| 106 | 132 | 132 | 33 | 7 |
| 390 | 35 | 32 | 8 | 25 |
| 2500 | 5 | 4 | 1 | 200 |

| cells | dose | draws | mutually disjoint (largest one-permutation set) | permutations spanned | mean pairwise overlap / n | max pairwise overlap / n | share of pairs disjoint | n / 14013 |
|---|---|---|---|---|---|---|---|---|
| 106 | 0.000 | 200 | 33.000 | 7.000 | 0.006 | 0.057 | 0.528 | 0.008 |
| 106 | 0.250 | 200 | 33.000 | 7.000 | 0.006 | 0.057 | 0.533 | 0.008 |
| 106 | 0.500 | 200 | 33.000 | 7.000 | 0.006 | 0.057 | 0.537 | 0.008 |
| 106 | 1.000 | 200 | 33.000 | 7.000 | 0.006 | 0.057 | 0.529 | 0.008 |
| 390 | 0.000 | 200 | 8.000 | 25.000 | 0.027 | 0.069 | 0.035 | 0.028 |
| 390 | 0.250 | 200 | 8.000 | 25.000 | 0.027 | 0.064 | 0.035 | 0.028 |
| 390 | 0.500 | 200 | 8.000 | 25.000 | 0.027 | 0.064 | 0.035 | 0.028 |
| 390 | 1.000 | 200 | 8.000 | 25.000 | 0.027 | 0.064 | 0.035 | 0.028 |
| 2500 | 0.000 | 200 | 1.000 | 200.000 | 0.178 | 0.208 | 0.000 | 0.178 |
| 2500 | 0.250 | 200 | 1.000 | 200.000 | 0.178 | 0.207 | 0.000 | 0.178 |
| 2500 | 0.500 | 200 | 1.000 | 200.000 | 0.178 | 0.207 | 0.000 | 0.178 |
| 2500 | 1.000 | 200 | 1.000 | 200.000 | 0.178 | 0.204 | 0.000 | 0.178 |

Draws within one permutation share no cell, across doses too. Draws in the same condition from different permutations overlap at random, at about n/14,013 of their cells per pair. At 2,500 cells only one draw per dose fits in a permutation. Each of the 200 draws per dose there comes from its own permutation, and any two share about 18% of their cells. The primary reference is the complement of the draw, so a draw never shares a cell with its own reference.

What the replicate spread does and does not include:

1. Every draw resamples the same 14,013 cells. The spread is the spread of resampling this pool, not of new cells from the line. Given the pool, draws are independent, so overlap does not shrink the spread around the pool mean. But no draw can show variation the pool lacks. At 2,500 cells, 200 draws cover each cell about 36 times per dose.
2. Random cells mix all 78 guides, so the primary null carries no guide-to-guide variation. A factor's cells come from about 6 guides. R1 below restores that structure; compare its spread with the primary.
3. A draw's reference is the other NT cells. A real factor is compared with all 14,013. The primary reference is smaller by n cells, which adds a little noise at 2,500 cells (11,513 reference cells).

## 6. Report-only readings (pre-registered, no decision rule, never a gate)

### R1 guide-structured null

| space | cells | monotone | power 25% | power 50% | power 100% | smallest separable |
|---|---|---|---|---|---|---|
| S1 | 106 | True | 0.240 | 0.415 | 0.505 | none |
| S1 | 390 | True | 0.975 | 1.000 | 1.000 | 25% |
| S1 | 2500 | True | 1.000 | 1.000 | 1.000 | 25% |
| S2_k20 | 106 | True | 0.990 | 1.000 | 1.000 | 25% |
| S2_k20 | 390 | True | 1.000 | 1.000 | 1.000 | 25% |
| S2_k20 | 2500 | True | 1.000 | 1.000 | 1.000 | 25% |
| S2_k50 | 106 | True | 0.990 | 1.000 | 1.000 | 25% |
| S2_k50 | 390 | True | 1.000 | 1.000 | 1.000 | 25% |
| S2_k50 | 2500 | True | 1.000 | 1.000 | 1.000 | 25% |
| S2_k100 | 106 | True | 0.960 | 1.000 | 1.000 | 25% |
| S2_k100 | 390 | True | 1.000 | 1.000 | 1.000 | 25% |
| S2_k100 | 2500 | True | 1.000 | 1.000 | 1.000 | 25% |

Dose-0 spread against the primary:

| space | cells | median delta at 0% | SD at 0% | primary SD at 0% | SD ratio |
|---|---|---|---|---|---|
| S1 | 106 | 93.798 | 2.765 | 2.735 | 1.011 |
| S1 | 390 | 37.897 | 0.806 | 0.906 | 0.890 |
| S1 | 2500 | 9.567 | 0.318 | 0.247 | 1.287 |
| S2_k20 | 106 | 2.913 | 1.568 | 1.567 | 1.000 |
| S2_k20 | 390 | 0.471 | 0.928 | 0.737 | 1.259 |
| S2_k20 | 2500 | 0.090 | 0.447 | 0.322 | 1.389 |
| S2_k50 | 106 | 4.650 | 1.580 | 1.729 | 0.913 |
| S2_k50 | 390 | 0.539 | 0.847 | 0.659 | 1.286 |
| S2_k50 | 2500 | 0.105 | 0.379 | 0.279 | 1.356 |
| S2_k100 | 106 | 6.359 | 1.738 | 1.930 | 0.901 |
| S2_k100 | 390 | 0.845 | 0.867 | 0.664 | 1.305 |
| S2_k100 | 2500 | 0.136 | 0.375 | 0.275 | 1.362 |

| cells | dose | guides per draw (median) | min | max | mean pairwise overlap / n | share of pairs disjoint |
|---|---|---|---|---|---|---|
| 106 | 0.000 | 6.000 | 6.000 | 6.000 | 0.008 | 0.672 |
| 106 | 0.250 | 6.000 | 6.000 | 6.000 | 0.007 | 0.683 |
| 106 | 0.500 | 6.000 | 6.000 | 6.000 | 0.008 | 0.677 |
| 106 | 1.000 | 6.000 | 6.000 | 6.000 | 0.008 | 0.680 |
| 390 | 0.000 | 6.000 | 6.000 | 6.000 | 0.028 | 0.613 |
| 390 | 0.250 | 6.000 | 6.000 | 6.000 | 0.028 | 0.607 |
| 390 | 0.500 | 6.000 | 6.000 | 6.000 | 0.028 | 0.604 |
| 390 | 1.000 | 6.000 | 6.000 | 6.000 | 0.029 | 0.609 |
| 2500 | 0.000 | 14.000 | 11.000 | 18.000 | 0.177 | 0.042 |
| 2500 | 0.250 | 15.000 | 10.000 | 19.000 | 0.179 | 0.040 |
| 2500 | 0.500 | 14.000 | 10.000 | 18.000 | 0.180 | 0.041 |
| 2500 | 1.000 | 14.500 | 11.000 | 19.000 | 0.178 | 0.039 |

### R2 pool reference (draw inside its own reference)

| space | cells | monotone | power 25% | power 50% | power 100% | smallest separable |
|---|---|---|---|---|---|---|
| S1 | 106 | True | 0.115 | 0.340 | 0.470 | none |
| S1 | 390 | True | 0.965 | 1.000 | 1.000 | 25% |
| S1 | 2500 | True | 1.000 | 1.000 | 1.000 | 25% |
| S2_k20 | 106 | True | 0.990 | 1.000 | 1.000 | 25% |
| S2_k20 | 390 | True | 1.000 | 1.000 | 1.000 | 25% |
| S2_k20 | 2500 | True | 1.000 | 1.000 | 1.000 | 25% |
| S2_k50 | 106 | True | 0.965 | 1.000 | 1.000 | 25% |
| S2_k50 | 390 | True | 1.000 | 1.000 | 1.000 | 25% |
| S2_k50 | 2500 | True | 1.000 | 1.000 | 1.000 | 25% |
| S2_k100 | 106 | True | 0.940 | 1.000 | 1.000 | 25% |
| S2_k100 | 390 | True | 1.000 | 1.000 | 1.000 | 25% |
| S2_k100 | 2500 | True | 1.000 | 1.000 | 1.000 | 25% |

Dose-0 spread against the primary:

| space | cells | median delta at 0% | SD at 0% | primary SD at 0% | SD ratio |
|---|---|---|---|---|---|
| S1 | 106 | 92.773 | 2.746 | 2.735 | 1.004 |
| S1 | 390 | 36.184 | 0.899 | 0.906 | 0.992 |
| S1 | 2500 | 6.678 | 0.207 | 0.247 | 0.839 |
| S2_k20 | 106 | 2.510 | 1.555 | 1.567 | 0.993 |
| S2_k20 | 390 | 0.390 | 0.717 | 0.737 | 0.973 |
| S2_k20 | 2500 | 0.056 | 0.264 | 0.322 | 0.819 |
| S2_k50 | 106 | 4.199 | 1.717 | 1.729 | 0.993 |
| S2_k50 | 390 | 0.516 | 0.641 | 0.659 | 0.974 |
| S2_k50 | 2500 | 0.079 | 0.229 | 0.279 | 0.819 |
| S2_k100 | 106 | 6.077 | 1.917 | 1.930 | 0.993 |
| S2_k100 | 390 | 0.750 | 0.647 | 0.664 | 0.974 |
| S2_k100 | 2500 | 0.111 | 0.226 | 0.275 | 0.819 |

### R3 full-gap plant (includes unmeasurable coordinates)

| space | cells | monotone | power 25% | power 50% | power 100% | smallest separable |
|---|---|---|---|---|---|---|
| S1 | 106 | True | 0.135 | 0.400 | 0.535 | none |
| S1 | 390 | True | 0.985 | 1.000 | 1.000 | 25% |
| S1 | 2500 | True | 1.000 | 1.000 | 1.000 | 25% |
| S2_k20 | 106 | True | 0.990 | 1.000 | 1.000 | 25% |
| S2_k20 | 390 | True | 1.000 | 1.000 | 1.000 | 25% |
| S2_k20 | 2500 | True | 1.000 | 1.000 | 1.000 | 25% |
| S2_k50 | 106 | True | 0.975 | 1.000 | 1.000 | 25% |
| S2_k50 | 390 | True | 1.000 | 1.000 | 1.000 | 25% |
| S2_k50 | 2500 | True | 1.000 | 1.000 | 1.000 | 25% |
| S2_k100 | 106 | True | 0.965 | 1.000 | 1.000 | 25% |
| S2_k100 | 390 | True | 1.000 | 1.000 | 1.000 | 25% |
| S2_k100 | 2500 | True | 1.000 | 1.000 | 1.000 | 25% |

Dose-0 spread against the primary:

| space | cells | median delta at 0% | SD at 0% | primary SD at 0% | SD ratio |
|---|---|---|---|---|---|
| S1 | 106 | 93.603 | 2.735 | 2.735 | 1.000 |
| S1 | 390 | 37.674 | 0.906 | 0.906 | 1.000 |
| S1 | 2500 | 9.475 | 0.247 | 0.247 | 1.000 |
| S2_k20 | 106 | 2.518 | 1.567 | 1.567 | 1.000 |
| S2_k20 | 390 | 0.395 | 0.737 | 0.737 | 1.000 |
| S2_k20 | 2500 | 0.063 | 0.322 | 0.322 | 1.000 |
| S2_k50 | 106 | 4.209 | 1.729 | 1.729 | 1.000 |
| S2_k50 | 390 | 0.524 | 0.659 | 0.659 | 1.000 |
| S2_k50 | 2500 | 0.098 | 0.279 | 0.279 | 1.000 |
| S2_k100 | 106 | 6.075 | 1.930 | 1.930 | 1.000 |
| S2_k100 | 390 | 0.765 | 0.664 | 0.664 | 1.000 |
| S2_k100 | 2500 | 0.142 | 0.275 | 0.275 | 1.000 |

### R4 Mann-Whitney

Per condition in the Section 4 tables (one-sided p that dose draws sit below dose-0 draws, and AUC). With 200 against 200 draws, this tests a shift of the whole distribution, not whether one factor would be called.

### R5 Basis circularity: gap retention ||proj_k g|| / ||g||

| k | donors in SOUTH6 basis | retention, SOUTH6 basis | donors in 40-59 basis | retention, 40-59-only basis |
|---|---|---|---|---|
| 20 | 652 | 0.887 | 308 | 0.865 |
| 50 | 652 | 0.956 | 308 | 0.897 |
| 100 | 652 | 0.974 | 308 | 0.905 |

## 7. Limitations

- The plant is added in z after the transform. A real effect would pass through counts, TMM and the log, and would add its own cell-to-cell variance. Here only the null carries sampling noise.
- The plant points exactly along the gap. A real factor's displacement also carries components off the gap axis. In S1, the screen's median ||d|| (123.4) is about 3.8 times the gap (32.3), so off-axis movement dominates real vectors in a way this control does not include.
- The plant is restricted to measured genes. At 100% it closes only the measured part of the gap, so even a noise-free 100% plant leaves the distance above zero (see 'noise-free').
- One NT pool of 14,013 cells from one neonatal line, with no aged cells. Every delta is a counterfactual against GTEx anchors from tissue donors.
- The S2 basis shares donors with the anchors (Section 1, R5). The gate uses that basis as stored.
- Separability is called against the 5th percentile of 200 dose-0 draws. That percentile is itself estimated, and its uncertainty is carried only in the bootstrap CI on power, which is reported but not gated on.

## 8. Files

`src/south6rec_run.py`; `results/south6rec/`: `PREREG.flag`, `stage0_reproduction.csv`, `draw_design.csv`, `overlap.csv`, `stage1_draws.csv` (every draw x space x variant), `stage1_summary.csv`, `stage1_gate.csv`, `r5_basis_retention.csv`, `south6rec_run.log`.

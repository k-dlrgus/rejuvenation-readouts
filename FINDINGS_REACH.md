# FINDINGS_REACH — the reachability bound

**Status: STOP `solver_not_converged`.** set D target 11: a solve failed the duality-gap test (max gap 4.268e-06).

Per CLAUDE.md, nothing was computed after the STOP and the run was not restarted.

Reproduced by `python src/reach_run.py` then `python src/reach_run.py findings`. `python src/reach_run.py verify` re-solves sets A-D from tracked `results/reach/inputs_k20.npz` and applies the same convergence test, so the STOP is re-checkable from tracked files.

## Pre-registration (verbatim from `results/reach/PREREG.flag`)

```
TASK: REACH — the reachability bound. How close to the young anchor can any
weighted sum of Southard's measured factor displacements get?
Written 2026-10-17T22:22Z (2026-10-18 07:22 +0900), before any REACH number
exists. No threshold or wording in this block is changed after numbers exist.
Stored prior results read so far, all as written by earlier tasks:
results/south6/origin_target_distances.csv, stage2a_floor_summary.csv,
stage0_cells_per_target.csv (factor count 1,836; cell counts 13-341, 212
distinct values, largest multiplicity of one count 32), results/south7/
basis_retention.csv and the first rows of stage_a1_depth_nulls.csv, and the
code of src/south6_run.py, south6rec_run.py, south7_run.py, transfer_run.py,
seng_run.py. Solver behaviour was checked on synthetic random matrices only
(no project data) to choose the convergence test below.
results/transfer/NOT_A_RESULT_gate_failed/ is not opened.

--- WHAT THIS IS ---
The TF search model blueprint calls this the reachability bound and labels it
"sizing information, never evidence". That label is binding here. The
calculation assumes Southard's directions are measured rather than
estimated, and that they transfer to aged cells. Neither is conservative. A
positive result is an upper bound under optimistic assumptions; a negative
result does not show unreachability in aged cells.

No reading of this task licenses a claim about aged cells. A positive bound sizes an experiment; it does not support one.

--- INSTRUMENTS (read, never refit) ---
- d_f: results/south6/stage1.npz key `d_f` (1,836 x 23,485 frozen-z ruler
  coordinates; 0 off the 20,264 matched genes), labels `keep`.
- Space S2A_k20: L20 = results/south7/basis.npz `L_mid_full[:, :20]`.
  Projection is applied last and to differences only (SENG3 amendment B6):
  b = (c_young - c_old) @ L20; D20 = d_f @ L20; a null displacement E
  (a difference) is projected the same way. center_mid never enters.
- c_old, c_young: results/south7/basis.npz (the O1/Y1 GTEx anchor centroids
  of SOUTH6/SOUTH7).
- Per-factor cell counts: diff of stage1.npz `cell_indptr`, checked against
  results/south6/stage0_cells_per_target.csv.
- Control pool: stage1.npz `nt_cells` (14,013 non-targeting cells), counts
  read from the SOUTH6 matched-gene CSR cache (results/south6/matched_*);
  the panel TMM cache rebuilt from stored stage1.npz `C_f`, `C_nt`, `C_pool`
  exactly as south6rec_run.stage0 / transfer_run.stage0 do; z against the
  frozen GTEx mu/sd (fibro2_common.load_frozen_ruler), absent genes z = 0.
- Seed 20260914 everywhere (generators named below).

--- THE OPTIMISATION ---
g = ||b|| (S2A_k20). For w in [0, 1]^1836 (1 = the dose Southard delivered):
  frac_closed(w) = (g - ||D20^T w - b||) / g
(equivalently the predicted state c_old + sum_f w_f d_f, all projected,
against c_young). w = 0 gives 0, so every optimum is >= 0.
Sparsity caps k = 1, 3, 5, 10, 50 and unlimited; k caps the number of
non-zero weights (nnz <= k). No penalty other than the dose cap and k.
- Unlimited: bounded least squares, scipy.optimize.lsq_linear(method="bvls",
  bounds (0, 1), tol 1e-12, max_iter 10000). Exact.
- Greedy (finite k): start w = 0, active set S empty. Each step, for every
  j not in S, t_j = clip(-(a_j . e)/||a_j||^2, 0, 1) with e = A w - b
  (A = D20^T, columns a_j), gain_j = ||e||^2 - ||e + t_j a_j||^2. Add the j
  with the largest gain (lowest index on ties); stop early if no gain > 0.
  Re-optimise all weights in S by BVLS (same settings). Repeat until nnz(w)
  reaches k or no candidate gains. At k = 1 this is exhaustive and exact.
- L1 relaxation (finite k): P(lam): minimise 0.5||A w - b||^2 + lam*sum(w)
  over [0, 1]^n (sum(w) is the L1 norm since w >= 0), solved by L-BFGS-B
  (scipy.optimize.minimize, analytic gradient, ftol 1e-16, gtol 1e-12,
  maxiter 20000, maxcor 30), warm-started. lam_max = max_j (A^T b)_j; if
  lam_max <= 0 the support is empty. lam is bisected in log space on
  [lam_max * 1e-10, lam_max] for 60 steps to find nnz == k, nnz counting
  w_j > 1e-6. If no step gives exactly k, the step with the largest
  nnz <= k is used (smallest lam among ties); the nnz actually hit is
  reported. The support's weights are then refit by BVLS with no penalty
  (relaxed L1), and frac is computed from the refit.
- Gated value at finite k: max(greedy, L1-refit); both are feasible, so the
  better is the tighter lower estimate of the k-sparse optimum. Both are
  reported, with |greedy - L1| and the Jaccard overlap of their supports.
- Convergence (every solve, BVLS and L1): certified by duality gap. For
  primal w and y = A w - b, P = 0.5||y||^2 + lam*sum(w), dual
  D = -0.5||y||^2 - b.y + sum_j min(0, a_j.y + lam) (lam = 0 for BVLS;
  restricted to the active columns for greedy re-fits). Converged iff
  P - D <= 1e-8 * max(1, 0.5||b||^2), and for BVLS also status >= 1.
  L-BFGS-B is restarted from its last point up to 20 times to meet this.

--- FOUR VECTOR SETS, identical optimisation ---
(A) RAW. D20 as stored.
(B) SHRUNK. column f scaled by s_f = max(0, 1 - sigma_f^2 / ||D20_f||^2).
    sigma_f^2: results/south6/stage2a_floors.npz holds gene-label
    permutation delta floors (keys floors, les, pooled, medians), not indexed
    by cell count and not squared norms; it does not supply sigma_f^2. So
    sigma_f^2 is built with this task's code (deviation, stated): for each
    distinct cell count n, the 200 depth-matched non-targeting draws of
    south7_run.null_draws(n, 14013) (SOUTH6REC primary scheme: disjoint
    n-cell blocks of permutations from default_rng([20260914, n])), each
    turned into E = z(draw counts) - z(C_pool - draw counts) with the SOUTH6
    panel TMM cache; sigma^2(n) = mean over the 200 draws of ||E @ L20||^2.
    sigma_f^2 = sigma^2(n_f), exact per n, no interpolation.
(C) NULL FACTOR BANK. 1,836 displacements from non-targeting cells only,
    one per factor at that factor's cell count, by the same construction as
    (B)'s draws. "The SOUTH6 control-draw construction" names no
    depth-matched draw in src/south6_run.py (its control draws are
    7,006/7,007-cell half splits and leave-one-guide-out); the construction
    used is the one SOUTH6REC built for the SOUTH6 readout and SOUTH7 used
    for its depth-matched floor (fixed here, stated). Factor i with count n
    takes draw index m_i = number of earlier factors (keep order) with the
    same n (max 32 < 200). These are the same draws that enter sigma^2(n);
    sigma enters only (B). (C) is the null analogue of (A): its vectors are
    NOT shrunk. R0 compares unshrunk (C) with shrunk (B) as the task text
    writes; this is the stricter comparison for the bound.
(D) RANDOM TARGETS. (B)'s vectors against 200 targets u_t * g, u_t = z/||z||,
    z ~ N(0, I_20) from default_rng([20260914, 4]) (one call,
    standard_normal((200, 20))). Every k, both finite-k methods, same
    gated value rule.

--- STAGE 0 (before anything new) ---
Counts exact; numerics within 1e-4 * max(1, |stored|). Any gated mismatch:
STOP prior_results_not_reproduced, naming every one.
 (a) factor count: d_f rows = len(keep) = csv factor rows = 1,836; d_f has
     23,485 columns; cells per factor equal stage0_cells_per_target.csv.
 (b) basis: max|L20^T L20 - I| <= 1e-8, else STOP basis_not_orthonormal
     (the threshold TRANSFER used on the same basis).
 (c) ||O1 - Y1||: S1 (all 23,485 coordinates) against
     results/south6/origin_target_distances.csv (S1, O1, Y1). S2A_k20:
     that file has no S2A row (its k20 row is SOUTH6's S2_k20 basis), so the
     S2A_k20 distance is checked against results/south7/basis_retention.csv
     (k = 20, origin_target_distance) — deviation from the task text, stated.
     Extra: basis.npz `rnorm` (S2A_k20 entry) and vector `r_S2A_k20`.
 (d) 20 per-factor norms: factors default_rng([20260914, 0]).choice(1836,
     20, replace=False); ||d_f|| (S1) against results/south6/
     stage3_full_table.csv d_norm (space S1, pair O1_Y1) and ||d_f @ L20||
     against results/south7/stage_a1_observed.csv dnorm_S2A_k20.
Extra checks, also gated:
 (e) the same two norms for all 1,836 factors;
 (f) panel rebuilt from stored counts: missing-gene mask identical; rebuilt
     d_f on matched genes within tol of stored; X_nt has 14,013 rows and
     X_nt.sum(0) equals stored C_pool exactly;
 (g) the null construction: for all 212 distinct n, this task's draws
     reproduce results/south7/stage_a1_depth_nulls.csv (floor_p05, median,
     sd, p95; spaces S2A_k20 and S1), and per-draw deltas in
     results/south7/stage_a1_null_draws.npz (keys "n|S2A_k20", "n|S1").
     (g) is checked after the draws are built and before any sigma^2,
     null bank vector or frac is written or used.

--- SOLVER CHECK (added; the task defines no positive control) ---
The bound is an optimisation, not a test, so recovery power in the usual
sense does not apply; (C) and (D) are its comparisons. To show the optimiser
finds a direction that is there: (B)'s matrix with one extra column equal
to b. Greedy and L1 at k = 1 must each return frac >= 1 - 1e-6 and select
that column. Failure is STOP solver_not_converged.

--- READINGS (evaluated in this order; only the one that fires is the
conclusion) ---
Let B(k), C(k) be the gated values for sets (B), (C) on the young target at
k in {1, 3, 5, 10, 50, unlimited}; Dp5(k), Dp95(k) the 5th and 95th
percentiles (numpy default linear) of the 200 random-target values at k.

(R0) UNINTERPRETABLE fires if, at any k, C(k) >= B(k) - 1e-6 (a tie,
  including both reaching the target exactly, counts as "at least as far"),
  or if, at any k, Dp5(k) - 1e-6 <= B(k) <= Dp95(k) + 1e-6. ("At any k" is
  applied to both clauses; this is the stricter reading.)
  Wording: "UNINTERPRETABLE. The bound cannot distinguish the measured
  factors from noise vectors of the same depth, so it carries no
  information about whether young is reachable, and no number from (A) or
  (B) may be quoted as a bound. Fired because: {list of clauses and k at
  which each fired}."

(R1) BOUND REPORTED fires otherwise.
  Wording: "BOUND REPORTED. In S2A_k20, with Southard's factor
  displacements shrunk for depth noise (set B), the largest fraction of the
  O1-Y1 gap closed by any combination with doses between 0 and the
  delivered dose is {B(k) for k = 1, 3, 5, 10, 50, unlimited}. Beside it:
  the null factor bank (set C) closes {C(k), same k}; against 200 random
  targets of the same length (set D) the 5th-95th percentile range is
  {Dp5(k)-Dp95(k), same k} and the young target sits {above / below} that
  range at every k. This is an upper bound computed under assumptions known
  to be optimistic: that these directions are measured rather than
  estimated, that they add, and that they transfer to aged cells. It is
  sizing information for the combination screen. It is not evidence that
  any combination rejuvenates aged cells, and it must not be cited as
  such."

No reading of this task licenses a claim about aged cells. A positive bound sizes an experiment; it does not support one.

--- DIAGNOSTICS (not gated, labelled) ---
- frac for (A) beside (B), per k and method: what shrinkage costs.
- Which factors carry the weight at each k (sets A and B, both methods):
  factor, weight, cell count, cell-count percentile among the 1,836.
  Spearman (scipy.stats.spearmanr) between weight and cell count over the
  selected (non-zero) factors where at least 3 are selected, and over all
  1,836 (zeros included) at unlimited k. The project's shortlist-versus-
  depth check.
- s_f: count of factors shrunk to 0 in (B), with the chance baseline: the
  count of (C)'s null vectors that the same rule (their own n) shrinks to 0.
- (C) shrunk by the same rule, and (A) vs (C): not gated.
- Full gene space S1 (all 23,485 coordinates, b = c_young - c_old): sets A,
  B (sigma^2 from the same draws, full-space squared norms) and C, same
  algorithms, run in Gram-reduced form (G = A^T A = U L U^T, A' = L^1/2 U^T
  over eigenvalues > 1e-10 * max, b' = L^-1/2 U^T A^T b, constant
  ||b||^2 - ||b'||^2 added back). D is not run in S1. Caveat: S1 includes
  the 3,221 ruler genes the factors cannot move (z = 0 in d_f) and many
  noise dimensions; its numbers are not comparable to S2A_k20 and are not
  gated. A solver failure here is recorded, not a STOP.

--- STOP KEYS ---
prior_results_not_reproduced: any gated Stage 0 mismatch (every one named).
basis_not_orthonormal: Stage 0 (b).
solver_not_converged: any gated solve (sets A-D in S2A_k20, the solver
  check) fails the convergence test.
After a STOP: write the findings with the STOP key and exit. No computation
after a STOP, no restart with an adjustment.

--- WHAT THIS DOES NOT COVER ---
- Southard is Hs27 neonatal foreskin fibroblasts, not aged donors. Nothing
  here transfers to aged cells.
- The transfer check (FINDINGS_TRANSFER.md, commit 7230271) could not
  measure whether these directions carry to another fibroblast context.
- Doses above the delivered one are not modelled; w <= 1 is the cap, and
  any w in (0, 1) assumes a partial dose acts as a scaled displacement,
  itself an assumption.
- Additivity of displacements is assumed, not tested.
- Nothing here models toxicity, identity loss or proliferation.
- Circularity, disclosed: the null draws, sigma^2 and (C) all come from the
  same 14,013 non-targeting cells that form d_f's reference pool; (C)'s
  vectors are among the draws that estimate sigma^2.

--- OUTPUTS ---
src/reach_run.py; results/reach/ (PREREG.flag, stage0_reproduction.csv,
null_depths.csv, factors.csv, inputs_k20.npz, frac_table.csv,
random_targets.csv, weights.csv, spearman.csv, fullspace.csv,
readings.json, run.log, STOP.json if a STOP fires); FINDINGS_REACH.md with
this block verbatim at the top. One local commit, not pushed.
```

## Deviations and decisions, disclosed

1. **S2A_k20 O1-Y1 distance compared against `results/south7/basis_retention.csv`, not `results/south6/origin_target_distances.csv`.** The task named the SOUTH6 file for both S1 and S2A_k20; that file has no S2A row (its k20 row is SOUTH6's S2_k20 basis, 28.645). Fixed in the flag before any number; S1 is compared against the SOUTH6 file as asked.
2. **sigma_f^2 built with this task's code.** `results/south6/stage2a_floors.npz` holds gene-label permutation delta floors (keys `floors`, `les`, `pooled`, `medians`), not indexed by cell count and not squared norms, so it does not supply sigma_f^2. sigma^2(n) is the mean of ||E L20||^2 over the 200 depth-matched non-targeting draws of `south7_run.null_draws(n, 14013)` at each factor's own n; no interpolation.
3. **"The SOUTH6 control-draw construction" read as the SOUTH6REC/SOUTH7 depth-matched construction.** `src/south6_run.py` has no draw at a factor's cell count (its control draws are 7,006/7,007-cell half splits and leave-one-guide-out). The construction used is SOUTH6REC's primary scheme, which SOUTH7 used for its depth-matched floor: E = z(draw) - z(C_pool - draw) with the SOUTH6 panel TMM cache. Its reproduction of SOUTH7's stored null deltas is gated in Stage 0 (g).
4. **(C) is not shrunk.** The task did not say; the flag fixed (C) as the null analogue of (A), so R0 compares unshrunk noise with shrunk real factors (the stricter comparison for the bound). (C) shrunk by the same rule is a labelled diagnostic.
5. **(C) reuses the sigma^2 draws.** Factor i with count n takes draw index m_i (its rank among factors with that n, keep order; max 32 of 200). The same draws estimate sigma^2(n), which enters only (B).
6. **Decisions the task did not specify, fixed in the flag:** gated value at finite k = max(greedy, L1-refit); L1 support refit by BVLS without penalty; L1 bisection rule and fallback when nnz never equals k; nnz counts weights > 1e-6 (weights at or below 1e-6 are zeroed after every BVLS fit); convergence certified by duality gap <= 1e-8 * max(1, 0.5 g^2); R0 tie tolerance 1e-6 with a tie (including both sets reaching the target exactly) counting as "at least as far"; "at any k" applied to both R0 clauses; the D clause inclusive with 1e-6 slack; random targets from default_rng([20260914, 4]); Stage 0 sample from default_rng([20260914, 0]).
7. **Added beyond the task text:** the planted-column solver check (the task defines no positive control); gated extra Stage 0 checks (e)-(g) (all-factor norms, panel rebuild from stored counts, reproduction of SOUTH7's stored null deltas); `stage0_norm_sample.csv`, `null_depths.csv`, `factors.csv`, `inputs_k20.npz`, `solver_check.json`, `diagnostics.json`, `verify` mode.
8. **Imported unchanged:** `south6_run` (block_delta, open_matched), `south6rec_run.indicator`, `south7_run.null_draws`, `seng_run` (tmm_cache_from_panel, z_from_counts), `toward_run.tmm_logcpm_quiet`, `fibro2_common.load_frozen_ruler`. Importing `south7_run` runs its `OUT.mkdir` on `results/south7/` (already present); nothing is written there.
9. **Where written:** only `results/reach/`, `src/reach_run.py`, `FINDINGS_REACH.md`. `results/transfer/NOT_A_RESULT_gate_failed/` was not opened. The pool for set (D) runs 16 single-threaded workers; BLAS thread count can change floating-point summation order, not results beyond ~1e-12.
10. **Script edited while the run was in progress.** The findings renderer (`render`, `_stop_section`, the lists in this section) was added to `src/reach_run.py` after the run started, at Stage 0 (f)/(g). The main process had already loaded the earlier file; set (D)'s spawned workers re-import the script from disk and so loaded the edited file. No solver, draw or Stage 0 function was changed by those edits; the committed file is the edited one.
11. **Flag wording slip, not edited:** the flag says S1 includes "the 3,221 ruler genes the factors cannot move". 3,221 is the unmatched count (23,485 - 20,264); genes in the missing-gene mask (zero in every panel row) are also z = 0 in every d_f, so the untouchable set is at least that large.

## Stage 0 — reproduction of stored numbers

Counts exact; numerics within 1e-4 x max(1, |stored|). Every row is gated; rows marked `extra` were added beyond the task text. For multi-value rows `stored`/`reproduced` are blank and `n_mismatch` counts values outside tolerance. File: `results/reach/stage0_reproduction.csv`.

| item | kind | extra | n | stored | reproduced | max_scaled_diff | n_mismatch | ok |
|---|---|---|---|---|---|---|---|---|
| (a) factor count: d_f rows | exact |  | 1 | 1836.0000 | 1836.0000 |  | 0 | True |
| (a) factor count: len(keep) | exact |  | 1 | 1836.0000 | 1836.0000 |  | 0 | True |
| (a) factor count: stage0_cells_per_target.csv factor rows | exact |  | 1 | 1836.0000 | 1836.0000 |  | 0 | True |
| (a) d_f columns (ruler coordinates) | exact |  | 1 | 23485.0000 | 23485.0000 |  | 0 | True |
| (a) cells per factor equal stage0_cells_per_target.csv | exact |  | 1836 |  |  |  | 0 | True |
| (b) basis orthonormality max|L20^T L20 - I| <= 1e-8 | threshold |  | 400 | 1.00e-08 | 3.55e-15 |  | 0 | True |
| (c) ||O1 - Y1|| S1 vs south6/origin_target_distances.csv | tol |  | 1 | 32.2839 | 32.2839 | 0.0000 | 0 | True |
| (c) ||O1 - Y1|| S2A_k20 vs south7/basis_retention.csv (deviation: not in the south6 file) | tol |  | 1 | 27.9292 | 27.9292 | 0.0000 | 0 | True |
| (c) ||O1 - Y1|| S2A_k20 vs south7/basis.npz rnorm | tol | extra | 1 | 27.9292 | 27.9292 | 0.0000 | 0 | True |
| (c) r_S2A_k20 vector vs south7/basis.npz | tol | extra | 20 |  |  | 0.0000 | 0 | True |
| (c) r_S1 on matched genes vs south7/basis.npz | tol | extra | 20264 |  |  | 0.0000 | 0 | True |
| (d) ||d_f|| S1, 20 sampled factors, vs south6/stage3_full_table.csv d_norm | tol |  | 20 |  |  | 1.84e-15 | 0 | True |
| (d) ||d_f L20|| S2A_k20, 20 sampled factors, vs south7/stage_a1_observed.csv | tol |  | 20 |  |  | 5.76e-16 | 0 | True |
| (e) ||d_f|| S1, all 1836 factors, vs south6/stage3_full_table.csv | tol | extra | 1836 |  |  | 2.58e-15 | 0 | True |
| (e) ||d_f L20|| S2A_k20, all 1836 factors, vs south7/stage_a1_observed.csv | tol | extra | 1836 |  |  | 8.64e-16 | 0 | True |
| (f) stage1.npz matched_ruler equals basis.npz matched_ruler | exact | extra | 20264 |  |  |  | 0 | True |
| (f) missing-gene mask identical to stage1.npz | exact | extra | 20264 |  |  |  | 0 | True |
| (f) rebuilt d_f on matched genes vs stage1.npz d_f | tol | extra | 37204704 |  |  | 0.0000 | 0 | True |
| (f) X_nt rows (non-targeting cells) | exact | extra | 1 | 14013.0000 | 14013.0000 |  | 0 | True |
| (f) X_nt column sums equal stage1.npz C_pool | exact | extra | 20264 |  |  |  | 0 | True |
| (g) SOUTH7 depth nulls S2A_k20 floor_p05, 212 depths vs south7/stage_a1_depth_nulls.csv | tol | extra | 212 |  |  | 2.22e-16 | 0 | True |
| (g) SOUTH7 depth nulls S2A_k20 median, 212 depths vs south7/stage_a1_depth_nulls.csv | tol | extra | 212 |  |  | 1.96e-16 | 0 | True |
| (g) SOUTH7 depth nulls S2A_k20 sd, 212 depths vs south7/stage_a1_depth_nulls.csv | tol | extra | 212 |  |  | 2.16e-16 | 0 | True |
| (g) SOUTH7 depth nulls S2A_k20 p95, 212 depths vs south7/stage_a1_depth_nulls.csv | tol | extra | 212 |  |  | 2.43e-16 | 0 | True |
| (g) SOUTH7 per-draw deltas S2A_k20, 212 x 200 draws vs south7/stage_a1_null_draws.npz (n outside tol) | tol | extra | 42400 | 0.0000 | 0.0000 | 0.0000 | 0 | True |
| (g) SOUTH7 depth nulls S1 floor_p05, 212 depths vs south7/stage_a1_depth_nulls.csv | tol | extra | 212 |  |  | 2.19e-16 | 0 | True |
| (g) SOUTH7 depth nulls S1 median, 212 depths vs south7/stage_a1_depth_nulls.csv | tol | extra | 212 |  |  | 2.19e-16 | 0 | True |
| (g) SOUTH7 depth nulls S1 sd, 212 depths vs south7/stage_a1_depth_nulls.csv | tol | extra | 212 |  |  | 2.25e-16 | 0 | True |
| (g) SOUTH7 depth nulls S1 p95, 212 depths vs south7/stage_a1_depth_nulls.csv | tol | extra | 212 |  |  | 2.17e-16 | 0 | True |
| (g) SOUTH7 per-draw deltas S1, 212 x 200 draws vs south7/stage_a1_null_draws.npz (n outside tol) | tol | extra | 42400 | 0.0000 | 0.0000 | 0.0000 | 0 | True |

30 checks, 30 gated, 0 gated failures.

The 20 sampled factors for (d) (default_rng([20260914, 0])): AATF (79), GTF2A1 (198), HMGB1 (91), JDP2 (331), KLF2 (152), LHX6 (129), MAML1 (113), PIAS3 (133), SREBF2 (88), TAF7 (106), THAP10 (87), TOX2 (100), TTF1 (71), VAX2 (109), ZFHX3 (50), ZNF202 (141), ZNF394 (132), ZNF418 (78), ZNF614 (136), ZNF671 (124) — name (cells). Values: `results/reach/stage0_norm_sample.csv`.

## What happened, in order

1. Stage 0 (a)-(g): 30 of 30 gated checks passed (table above), including the reproduction of SOUTH7's stored depth-matched null deltas for all 212 cell counts.
2. Inputs were built and written before any solve: `null_depths.csv` (sigma^2(n)), `factors.csv` (per-factor ||D20_f||, sigma^2, s_f, null-bank draw index and norm), `inputs_k20.npz`.
3. Planted-column solver check passed (greedy and L1 at k=1 both selected the planted column; passed = True). `solver_check.json`.
4. Sets (A), (B), (C) on the young target were solved and each passed the convergence test. Their values were printed to `run.log` as the code logs every set. **They are not results.** No reading was evaluated: both R0 clauses and R1 need set (D), and only a pre-registered reading whose conditions fired is a conclusion. They are not tabulated here and are not a bound.
5. Set (D): 200 random targets in a pool of 16 workers, results consumed in target order. Target 11 (0-based) returned a failed convergence test: the largest duality gap over that target's solves was 4.268e-06 against the pre-registered tolerance 1e-8 x max(1, 0.5 g^2) = 3.900e-06 (ratio 1.09). The pool was terminated and STOP `solver_not_converged` written. Which of that target's solves failed (a greedy BVLS refit, an L1 L-BFGS-B solve after 20 restarts, the L1 refit, or the unlimited BVLS) is not recorded: the worker returns only the maximum gap. Targets 0-10 had completed; results of targets in flight on other workers were discarded unrecorded; `random_targets.csv` was not written.
6. After the STOP nothing was computed: no rerun, no change of tolerance or solver, no re-solve of target 11 to find which step failed.

## What the STOP does and does not say

- It is a failure of the convergence certificate fixed in the flag, which was chosen from solver behaviour on synthetic matrices. The gap exceeded it by 9%. It is not evidence about reachability in either direction, and it leaves no reading on record.
- No frac_closed table across the four sets is reported, because (D) was not completed and the readings were never evaluated.
- Whether a new task should re-pose the bound with a different, pre-registered convergence rule is a decision for outside this session. This run is not to be resumed or adjusted.
- The STOP is reproducible from tracked files: `python src/reach_run.py verify` re-solves from `inputs_k20.npz` and runs the same convergence test (not run in this session, as it would compute after the STOP).

## What this does not cover

- Southard is Hs27 neonatal foreskin fibroblasts, not aged donors. Nothing here transfers to aged cells.
- The transfer check (FINDINGS_TRANSFER.md, commit 7230271) could not measure whether these directions carry to another fibroblast context.
- Doses above the delivered one are not modelled; any w in (0, 1) assumes a partial dose acts as a scaled displacement, itself an assumption. Doses above 1 would be extrapolation.
- Additivity of displacements is assumed, not tested.
- Nothing here models toxicity, identity loss or proliferation.
- Circularity, disclosed: the null draws, sigma^2 and (C) all come from the same 14,013 non-targeting cells that form d_f's reference pool; (C)'s vectors are among the draws that estimate sigma^2. The S2A basis was fit on GTEx donors aged 40-59 with no anchor donor in the fit (SOUTH7), so the space is not fit to the factors or the anchors.

## Files

`src/reach_run.py`; `results/reach/`: `PREREG.flag`, `STOP.json`, `factors.csv`, `inputs_k20.npz`, `null_depths.csv`, `run.log`, `run_stderr.txt`, `run_stdout.txt`, `solver_check.json`, `stage0_norm_sample.csv`, `stage0_reproduction.csv`. Not written, because the STOP fired first: `random_targets.csv`, `readings.json`, `frac_table.csv`, `weights.csv`, `spearman.csv`, `diagnostics.json`.

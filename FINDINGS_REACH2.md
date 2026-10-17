# FINDINGS_REACH2 — the reachability bound, rerun with a relative convergence rule

**Status: DONE.** Reading: `UNINTERPRETABLE` (R0).

This is sizing information, never evidence (TF search model blueprint). No reading of this task licenses a claim about aged cells. A positive bound sizes an experiment; it does not support one.

**Prior partial output, disclosed.** REACH solved sets A, B and C against the young target before its STOP and printed those values to `results/reach/run.log`. They are not a result. This session did not open `run.log` or `run_stdout.txt` (the same log stream), and this task's code does not read them. No number from REACH's partial solve reached this session by any route. The REACH numbers this session saw are its Stage 0 table, the STOP's gap (4.268e-06) and tolerance (3.900e-06), and g = 27.929154.

Reproduced by `python src/reach2_run.py` then `python src/reach2_run.py findings`. `python src/reach2_run.py verify` re-solves sets A-D from tracked `results/reach2/inputs_k20.npz` under the same convergence rule.

## Pre-registration (verbatim from `results/reach2/PREREG.flag`)

```
TASK: REACH2 — rerun the REACH reachability bound with a workable
convergence rule. Same question, same instruments, same readings.
Written 2026-10-17T22:59Z (2026-10-18 07:59 +0900), before any REACH2 number exists. No
threshold or wording in this block is changed after numbers exist.

--- WHY THIS TASK EXISTS ---
REACH (results/reach/PREREG.flag, FINDINGS_REACH.md, commit ddec105) stopped
with solver_not_converged: one solve among the 200 random-target controls
(set D, target 11) left a duality gap of 4.268e-06 against REACH's absolute
tolerance 1e-8 * max(1, 0.5 g^2) = 3.900e-06, 9% above it. That tolerance
was chosen by the previous session from synthetic tests. The science was not
in question; the tolerance was. REACH's readings stand unchanged and are
carried over below verbatim. REACH is not resumed: this is a new task with
its own directory, code and pre-registration; nothing in results/reach/ is
written.

--- WHAT THIS SESSION HAS READ, AND WHAT IT HAS NOT ---
Read before this flag: CLAUDE.md; results/reach/PREREG.flag;
FINDINGS_REACH.md; results/reach/stage0_reproduction.csv; src/reach_run.py
(code only); the listing of results/reach/ (file names and sizes); git log
and .gitignore; the Python environment (.venv: Python 3.12.7, numpy 2.5.3,
scipy 1.18.1, pandas 3.0.5, 18 logical CPUs).
Not opened: results/reach/run.log; results/reach/run_stdout.txt (same size
as run.log and, per src/reach_run.py, the same log lines printed to stdout,
so treated as run.log); results/reach/STOP.json, solver_check.json,
run_stderr.txt; results/transfer/NOT_A_RESULT_gate_failed/.
results/reach/inputs_k20.npz, factors.csv and null_depths.csv are read only
by this task's code, in Stage 0 (h) below. Per src/reach_run.py they were
written before REACH's first solve and hold inputs (D20, b20, g, sigma^2,
s_f, C20, the 200 targets, per-factor norms), not solve output.
PRIOR PARTIAL OUTPUT, disclosed: REACH solved sets A, B and C against the
young target before its STOP and printed those values to results/reach/
run.log. They are not a result: no reading fired. No number from that
partial solve has reached this session as of this flag. The only REACH
numbers seen are its Stage 0 table, the STOP's gap and tolerance, and
g = 27.929154 (Stage 0 row (c)).

=== CARRIED OVER VERBATIM: results/reach/PREREG.flag lines 15-217 ===
(copied by program, byte for byte. Parts superseded for this task are named
under OVERRIDES below; the carried text is not edited.)

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

=== END OF CARRIED-OVER TEXT ===

--- OVERRIDES: THE ONLY CHANGES (task items 1-3) ---
O1. Convergence. Replaces, in THE OPTIMISATION's "Convergence" bullet, the
  test "P - D <= 1e-8 * max(1, 0.5||b||^2)", and, in STOP KEYS, the
  solver_not_converged line (see O2). A solve converges iff
      P - D <= 1e-6 * max(P, 1.0)
  with P and D exactly as in the carried bullet (P = 0.5||y||^2 +
  lam*sum(w), lam = 0 for BVLS, restricted to the active columns for greedy
  re-fits) and, for BVLS, also status >= 1 (carried). "Its own objective
  value" is P, the primal objective of that solve at the w it returns (fixed
  here; the task text does not say primal or dual). In the S1 Gram-reduced
  diagnostic, P includes the constant 0.5*c0 that the reduction removes, so
  P is the full-space objective (fixed here). L-BFGS-B is still restarted
  from its last point up to 20 times to meet the test (carried).
  Arithmetic, stated before any number: REACH's tolerance was 3.900e-06 for
  every S2A_k20 solve. The new tolerance is larger when P > 3.900 and
  smaller when P < 3.900; for P <= 1 it is 1e-6 absolute. A BVLS solve that
  passes has frac within 0.5e-6 * (residual / g) of that problem's optimum
  when P >= 1, and within sqrt(2e-6) / g = 5.1e-5 when P < 1. The readings'
  tie tolerance is 1e-6.

O2. One pre-registered retry. "A solve" is one solver call: one BVLS call (a
  greedy re-fit at one step, an L1 support refit at one k, or the unlimited
  fit) or one L-BFGS-B solve at one lam (with its up to 20 restarts). If it
  fails O1 it is re-solved once, from the same starting point, with the
  strict setting; the strict result replaces the first if it passes.
  Base setting (carried from REACH): BVLS lsq_linear(method="bvls", bounds
  (0, 1), tol 1e-12, max_iter 10000); L-BFGS-B ftol 1e-16, gtol 1e-12,
  maxiter 20000, maxcor 30, with maxfun 15000 and maxls 20 (scipy defaults:
  REACH's flag did not name them and its code used them, so they stand).
  Strict setting: BVLS tol 1e-15, max_iter 100000. L-BFGS-B ftol 1e-20,
  gtol 1e-15, maxiter 200000, maxfun 200000 (maxcor 30, maxls 20 and up to
  20 restarts unchanged).
  A second failure of the same solve in a gated problem (sets A, B, C on the
  young target; the 200 set-D targets; the solver check) is STOP
  solver_not_converged. In the not-gated diagnostics ((C) shrunk; S1) a
  second failure is recorded, the row is marked not converged and the
  result with the smaller gap/tolerance ratio is kept, as REACH fixed for S1
  ("recorded, not a STOP").
  A retry is part of the solve it belongs to; it is not a resumption of a
  stopped run. Nothing is computed after a STOP fires.

O3. Failures must be diagnosable. Written whatever the outcome:
  results/reach2/solver_diagnostics.csv, one row per (space, set, target
  index, k, method). Target index is "young" for A, B, C, 0-199 for D,
  "planted" for the solver check; method is greedy or l1 at finite k, bvls
  at unlimited. Columns: number of solves in the row; the final gap (the
  largest gap any of its solves ended with, after any retry) with that
  solve's objective P and tolerance; the largest gap/tolerance ratio;
  whether a retry was used and how many; whether the row converged.
  A solve belongs to the row being computed when it ran: a greedy re-fit at
  step s to the smallest k in {1, 3, 5, 10, 50} with s <= k (s > 50: k =
  50); an L-BFGS-B solve to the k whose bisection ran it (a solve reused
  from an earlier k's cache is not counted again); the L1 refit for k to k;
  the unlimited BVLS to unlimited. Because greedy follows one path to
  k = 50 and L1 reuses solves and warm starts across k, a row's value also
  depends on earlier rows' solves. The full record is
  results/reach2/solver_calls.csv.gz: every solve with its row keys, kind,
  step or bisection index, lam, P, gap, tolerance, restarts, status, retry
  used, gap before retry, passed. The 20 largest gaps and the 20 largest
  gap/tolerance ratios among gated solves are reported in the findings even
  on a clean run (worst_gaps.csv).
  Added diagnostic (not gated): the count of gated solves that pass O1 but
  would have failed REACH's test (gap > 3.900e-06), and the count that fail
  O1 at the first attempt but would have passed REACH's.

--- STAGE 0 FOR REACH2 (REACH's (a)-(g) and solver check carried; added) ---
- (a)-(g) are re-run in full by this task's script from the stored inputs
  REACH named, and gated as REACH gated them.
- Row for row against results/reach/stage0_reproduction.csv (30 rows):
  item, kind, gated, extra and n identical; stored, reproduced,
  max_abs_diff and max_scaled_diff within 1e-4 * max(1, |REACH value|)
  where REACH's cell is a number, blank where REACH's is blank; n_mismatch
  and ok identical. Any difference: STOP prior_results_not_reproduced,
  naming every one. Written to stage0_vs_reach.csv.
- (h) added, gated: the inputs this task builds equal REACH's stored inputs
  (results/reach/inputs_k20.npz: keep and n_cells exact; D20, b20, g,
  sigma2, s_f, C20, targets within tolerance; results/reach/null_depths.csv
  and factors.csv every column, counts exact, numerics within tolerance).
- Then the planted-column solver check (carried text above), under O1-O2.
  It is the planted problem whose answer is known (w = 1 on the planted
  column, frac = 1).

--- STOP KEYS FOR REACH2 ---
prior_results_not_reproduced: any gated Stage 0 mismatch, (a)-(h) or the
  row-for-row comparison (every one named).
basis_not_orthonormal: Stage 0 (b).
solver_not_converged: a second failure of the same solve (O2) in a gated
  problem, or the solver check not selecting the planted column at
  frac >= 1 - 1e-6.
After a STOP: write the findings with the STOP key and the solver
diagnostics of what completed, and exit. No frac value of any set is
written after a STOP (fixed here). No computation after a STOP, no restart
with an adjustment.

--- READINGS: carried unchanged ---
R0 and R1, their wording, the evaluation order and the sentence "No reading
of this task licenses a claim about aged cells. A positive bound sizes an
experiment; it does not support one." are the carried text above. The R1
slot "{above / below}" when B sits above the D range at some k and below it
at others is filled as REACH's code (reach_run.reading_text) fills it:
"above that range at k = ... and below it at k = ..." (REACH's code, not
its flag; stated here).

--- EXECUTION (housekeeping, disclosed) ---
- src/reach2_run.py imports unchanged from src/reach_run.py: stage0_static,
  stage0_panel, null_draws_all and check (none writes a file); readings_of
  and reading_text (the R0/R1 rule and wording); frac_rows, weight_rows,
  spearman_rows, jaccard, fmt, md_table. Through them REACH's imports stand
  (south6_run, south6rec_run.indicator, south7_run.null_draws, seng_run,
  toward_run, fibro2_common). REACH's solvers (greedy, L1, unlimited BVLS,
  Gram reduction) are copied into reach2_run.py with only O1-O3 changed.
- The 203 gated problems (A, B, C young; D targets 0-199) run in one pool of
  16 spawned single-threaded workers (REACH ran A-C in the main process and
  D in the pool). Results are taken as they complete; on a STOP the pool is
  terminated and problems in flight are recorded as not completed.
- No frac value is logged or written before all 203 gated problems have
  converged.
- Written only to results/reach2/, src/reach2_run.py, FINDINGS_REACH2.md,
  and one .gitignore line if solver_calls.csv.gz is 1 MB or more.
- verify mode re-solves sets A-D from results/reach2/inputs_k20.npz under
  the same rules. If no STOP fires it is run once after the findings are
  written and its agreement is reported; a mismatch is disclosed and does
  not change the reading.

--- OUTPUTS ---
src/reach2_run.py; results/reach2/ (PREREG.flag, stage0_reproduction.csv,
stage0_vs_reach.csv, stage0_norm_sample.csv, null_depths.csv, factors.csv,
inputs_k20.npz, solver_check.json, solver_diagnostics.csv,
solver_calls.csv.gz, worst_gaps.csv, frac_table.csv, random_targets.csv,
weights.csv, spearman.csv, readings.json, diagnostics.json, run.log,
STOP.json if a STOP fires, verify.log); FINDINGS_REACH2.md with this block
verbatim at the top. One local commit, not pushed.
```

## Carried over from REACH, item by item

| item | source |
|---|---|
| THE OPTIMISATION: origin c_old, target c_young, b = (c_young - c_old) @ L20, g = ||b||, w in [0, 1]^1836 (dose cap), frac_closed(w) = (g - ||D20^T w - b||) / g; sparsity caps k = 1, 3, 5, 10, 50 and unlimited; unlimited by BVLS; greedy and L1 (L-BFGS-B, log-space bisection on lam, BVLS refit of the support) both run at finite k; gated value max(greedy, L1-refit) | results/reach/PREREG.flag lines 54-81 (flag text, copied) |
| The duality-gap definition (P, D, lam = 0 for BVLS, restricted to the active columns for greedy re-fits), BVLS status >= 1, and up to 20 L-BFGS-B restarts from the last point | results/reach/PREREG.flag lines 82-87; only the threshold is replaced (O1) |
| INSTRUMENTS: d_f and keep from results/south6/stage1.npz, L20 = south7 basis L_mid_full[:, :20], c_old and c_young from south7 basis.npz, projection of differences only, cell counts from cell_indptr, the non-targeting pool and the panel TMM cache rebuilt from stored counts | results/reach/PREREG.flag lines 36-52 |
| Set (A) raw D20 | results/reach/PREREG.flag line 90 |
| Set (B) shrinkage s_f = max(0, 1 - sigma_f^2 / ||D20_f||^2), sigma^2(n) the mean of ||E L20||^2 over the 200 depth-matched draws of south7_run.null_draws(n, 14013), exact per n | results/reach/PREREG.flag lines 91-101 |
| Set (C) null factor bank: the depth-matched construction (SOUTH6REC primary scheme, as SOUTH7), draw index m_i, not shrunk | results/reach/PREREG.flag lines 102-113 |
| Set (D) 200 random targets u_t * g from default_rng([20260914, 4]).standard_normal((200, 20)), against (B)'s vectors, every k, both methods | results/reach/PREREG.flag lines 114-117 |
| Stage 0 (a)-(g), tolerance 1e-4 * max(1, |stored|), Stage 0 sample default_rng([20260914, 0]) | results/reach/PREREG.flag lines 119-146 |
| The planted-column solver check (B's matrix plus a column equal to b; greedy and L1 at k = 1 must reach frac >= 1 - 1e-6 and select that column) | results/reach/PREREG.flag lines 148-153 |
| Readings R0 (UNINTERPRETABLE) and R1 (BOUND REPORTED), their wording, their evaluation order, the tie tolerance 1e-6, "at any k" on both R0 clauses, numpy default percentiles | results/reach/PREREG.flag lines 155-184 |
| The sentence "No reading of this task licenses a claim about aged cells. A positive bound sizes an experiment; it does not support one." | results/reach/PREREG.flag lines 34 and 186 |
| The R1 slot when B sits above the D range at some k and below at others | src/reach_run.py reading_text (REACH's code, not its flag); stated in the REACH2 flag |
| DIAGNOSTICS: (A) beside (B); weights and cell counts with Spearman; s_f zero count with the null-bank baseline; (C) shrunk; the full gene space S1 in Gram-reduced form (D not run there; failures recorded, not a STOP) | results/reach/PREREG.flag lines 188-206 |
| STOP keys prior_results_not_reproduced and basis_not_orthonormal; "no computation after a STOP, no restart with an adjustment" | results/reach/PREREG.flag lines 208-214 |
| WHAT THIS DOES NOT COVER (all six items) | results/reach/PREREG.flag lines 216-228 |
| Seed 20260914 and every generator named | results/reach/PREREG.flag line 52 and the set definitions |
| L-BFGS-B maxfun 15000 and maxls 20 (scipy defaults) | src/reach_run.py l1_solve (REACH's flag did not name them) |
| 16 single-threaded spawned workers for the pool; weights at or below 1e-6 zeroed after every BVLS fit | src/reach_run.py (N_WORKERS, _clean) |

## Deviations and decisions, disclosed

1. **Only the three pre-registered changes alter the computation** (O1 relative convergence test, O2 one retry with a stricter setting, O3 the per-solve record). Everything else is REACH's, listed under "Carried over".
2. **"Its own objective value" fixed as the primal objective P** of that solve at the w it returns; in the S1 Gram-reduced diagnostic P includes the constant 0.5 c0 the reduction removes. Fixed in the flag; the task text did not say primal or dual.
3. **"A solve" fixed as one solver call** (one BVLS call, or one L-BFGS-B solve at one lam with its up to 20 restarts); the retry re-solves that call from the same starting point. The per-row table (set, target, k, method) aggregates the calls of a row; the rule assigning calls to rows is in the flag (O3).
4. **Strict solver setting chosen by this session** (the task asked only for a tighter internal tolerance and more iterations): BVLS tol 1e-15, max_iter 100000; L-BFGS-B ftol 1e-20, gtol 1e-15, maxiter 200000, maxfun 200000. Not tuned on any data.
5. **Added beyond the task text:** Stage 0 (h), this task's inputs against REACH's stored `inputs_k20.npz`, `null_depths.csv` and `factors.csv` (gated); `stage0_vs_reach.csv`; `solver_calls.csv.gz`, `worst_gaps.csv`, `solver_summary.json`; the old-rule counterfactual count (not gated); `selftest` mode; `verify` mode as in REACH.
6. **Files of REACH's directory read by this task's code:** `stage0_reproduction.csv` (row-for-row), `inputs_k20.npz`, `null_depths.csv`, `factors.csv` (Stage 0 (h)). Per `src/reach_run.py` these were all written before REACH's first solve and hold no solve output. Not opened by this session or its code: `run.log`, `run_stdout.txt`, `run_stderr.txt`, `STOP.json`, `solver_check.json`.
7. **Execution changed, not the computation:** sets A, B, C now run in the same pool of 16 spawned single-threaded workers as set D (REACH ran them in the main process with multi-threaded BLAS), results are taken as they complete, and no frac value is logged before all 203 gated problems have converged. BLAS thread count can change floating-point summation order, not results beyond ~1e-12.
8. **Imported unchanged from `src/reach_run.py`:** `stage0_static`, `stage0_panel`, `null_draws_all`, `check`, `readings_of`, `reading_text`, `frac_rows`, `weight_rows`, `spearman_rows`, `jaccard`, `fmt`, `md_table`, `Stop`, `Log`, `_js` and its constants. Through them REACH's imports stand (`south6_run`, `south6rec_run.indicator`, `south7_run.null_draws`, `seng_run`, `toward_run`, `fibro2_common`); importing `south7_run`, `seng_run`, `toward_run` and `south6_run` runs their `mkdir` on directories that already exist; nothing is written there. REACH's `gate_stage0` is not called (it writes into `results/reach/`); this task has its own gate. REACH's solvers and Gram reduction are copied into `src/reach2_run.py` with only O1-O3 changed.
9. **Where written:** `results/reach2/`, `src/reach2_run.py`, `FINDINGS_REACH2.md`. Nothing in `results/reach/` or any other earlier results directory. `results/transfer/NOT_A_RESULT_gate_failed/` was not opened.
10. **The first flag-writing attempt failed** (a shell here-document was cut off before the flag file was created; no file was written and nothing had been computed). The flag was then assembled by a script that copies REACH's lines 15-217 byte for byte and asserts the copy.
11. **Synthetic tests before the real run, no project data.** `selftest` (added) runs the solver check, pool, readings and findings on random matrices (120 vectors, 24 targets) into the session scratchpad. Four runs: clean (2,357 solves, no retry); every first attempt forced to fail (all 2,357 retried, all passed); tolerance 1e-30 (STOP at the solver check); tolerance 1e-11 (STOP inside the pool, diagnostics table written). The forcing switches are environment variables (`REACH2_SELFTEST_FAIL_FIRST`, `REACH2_SELFTEST_GAP_REL`, `REACH2_OUT`, `REACH2_FINDINGS`); `run`, `verify` and `findings` refuse to start when any is set. Seen before the real run and recorded here: in the 1e-11 test the strict L-BFGS-B retry ended at exactly the gap of the first attempt, both ending in a line-search stop (scipy status 2). Where L-BFGS-B stops on its line search, a tighter tolerance and more iterations do not help, so O2's retry mainly helps BVLS and iteration-limited solves. The flag was not changed.
12. **A stray file, removed.** Testing that `verify` refuses an empty directory created `results/reach2/verify.log` with only its header line (no number). It was deleted before the real run.
13. **Script edited after the run, renderer only.** The run and the verify re-solve used `src/reach2_run.py` with SHA-256 a46a5d2b082bee8b27e3a42f3661ddee9fd20ffa88af6120f98d533f387ad154 (hashed at launch, unchanged at the end of the run). After verify finished, two edits were made: the hand-written notes (`NOTES`, labelled as notes in the findings) and string formatting of the counts in the set-D greedy/L1 table. No solver, Stage 0, reading or output-writing code was changed. The findings were then rewritten by `findings` mode, which computes nothing. The committed file is the edited one.
14. **Console output kept outside the repository.** The run's stdout (the same lines as `run.log`) and stderr went to the session scratchpad and are not committed. stderr was empty. `solver_calls.csv.gz` is 553 KB, under the 1 MB line, so it is tracked and no `.gitignore` line was added.

## Stage 0 — reproduction of stored numbers

Counts exact; numerics within 1e-4 x max(1, |stored|). Every row is gated; rows marked `extra` were added beyond the task text (REACH's or this task's). Rows (a)-(g) are REACH's 30 checks re-run by this task's script; rows (h) are added by REACH2. File: `results/reach2/stage0_reproduction.csv`.

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
| (h) inputs_k20.npz keep | exact | extra | 1836 |  |  |  | 0 | True |
| (h) inputs_k20.npz n_cells | exact | extra | 1836 |  |  |  | 0 | True |
| (h) inputs_k20.npz D20 | tol | extra | 36720 |  |  | 0.0000 | 0 | True |
| (h) inputs_k20.npz b20 | tol | extra | 20 |  |  | 0.0000 | 0 | True |
| (h) inputs_k20.npz g | tol | extra | 1 | 27.9292 | 27.9292 | 0.0000 | 0 | True |
| (h) inputs_k20.npz sigma2 | tol | extra | 1836 |  |  | 0.0000 | 0 | True |
| (h) inputs_k20.npz s_f | tol | extra | 1836 |  |  | 0.0000 | 0 | True |
| (h) inputs_k20.npz C20 | tol | extra | 36720 |  |  | 0.0000 | 0 | True |
| (h) inputs_k20.npz targets | tol | extra | 4000 |  |  | 0.0000 | 0 | True |
| (h) null_depths.csv n_cells | exact | extra | 212 |  |  |  | 0 | True |
| (h) null_depths.csv sigma2_S2A_k20 | tol | extra | 212 |  |  | 2.13e-16 | 0 | True |
| (h) null_depths.csv sigma2_S1 | tol | extra | 212 |  |  | 1.98e-16 | 0 | True |
| (h) null_depths.csv R | exact | extra | 212 |  |  |  | 0 | True |
| (h) factors.csv factor | exact | extra | 1836 |  |  |  | 0 | True |
| (h) factors.csv n_cells | exact | extra | 1836 |  |  |  | 0 | True |
| (h) factors.csv n_cells_pctile | tol | extra | 1836 |  |  | 1.89e-16 | 0 | True |
| (h) factors.csv dnorm_S2A_k20 | tol | extra | 1836 |  |  | 1.97e-16 | 0 | True |
| (h) factors.csv sigma2_S2A_k20 | tol | extra | 1836 |  |  | 2.13e-16 | 0 | True |
| (h) factors.csv s_f | tol | extra | 1836 |  |  | 1.11e-16 | 0 | True |
| (h) factors.csv null_draw_index | exact | extra | 1836 |  |  |  | 0 | True |
| (h) factors.csv null_dnorm_S2A_k20 | tol | extra | 1836 |  |  | 1.97e-16 | 0 | True |
| (h) factors.csv null_s | tol | extra | 1836 |  |  | 9.93e-17 | 0 | True |

52 checks, 52 gated, 0 gated failures.

### Row for row against `results/reach/stage0_reproduction.csv`

30 rows compared, 330 gated cells (item, kind, gated, extra, n, n_mismatch, ok identical; stored, reproduced, max_abs_diff, max_scaled_diff within tolerance or blank where REACH's is blank): **0 mismatches**. The `note` column (not gated) is identical in 30 of 30 rows. File: `results/reach2/stage0_vs_reach.csv`.

## Solver check (carried; the planted problem with a known answer)

Set (B)'s 1,836 vectors plus one planted column equal to the gap b (index 1836). Greedy k=1: frac 1.000000, selected [1836]; L1 k=1: frac 1.000000, selected [1836]. Passed: True. This shows the optimiser finds a direction that is present; the bound is an optimisation, not a test, so it has no recovery power in the usual sense. Sets (C) and (D) are its comparisons.

## Pre-registered reading: `UNINTERPRETABLE` (R0)

> UNINTERPRETABLE. The bound cannot distinguish the measured factors from noise vectors of the same depth, so it carries no information about whether young is reachable, and no number from (A) or (B) may be quoted as a bound. Fired because: the null bank (C) reaches at least as far as the shrunken real factors (B) at k = 1, 3, 5, 10, 50, unlimited; the young target sits inside the middle 90% of the random-target distribution (D) at k = 1, 3, 5, 10, 50, unlimited.

No reading of this task licenses a claim about aged cells. A positive bound sizes an experiment; it does not support one.

Conditions, per k (gated values; finite k = max(greedy, L1-refit), unlimited = exact BVLS):

| k | B shrunk | C null bank | C >= B - 1e-6 | D 5th pct | D 95th pct | B inside D middle 90% |
|---|---|---|---|---|---|---|
| 1 | 0.0670 | 0.0938 | True | 0.0197 | 0.1755 | True |
| 3 | 0.1487 | 0.1861 | True | 0.0274 | 0.2902 | True |
| 5 | 0.1669 | 0.1957 | True | 0.0299 | 0.3341 | True |
| 10 | 0.1738 | 0.1958 | True | 0.0299 | 0.3871 | True |
| 50 | 0.1738 | 0.1958 | True | 0.0299 | 0.3946 | True |
| unlimited | 0.1738 | 0.1958 | True | 0.0299 | 0.3946 | True |

Under R0 the numbers from (A) and (B) are not a bound and must not be quoted as one. They appear below only because the hand-off requires the table across sets; they are labelled as such.

## frac_closed across k, all four sets (S2A_k20)

g = ||c_old - c_young|| in S2A_k20 = 27.929154. Gated value per set and k. (D) is (B)'s vectors against 200 random targets of length g; its columns summarise those 200 values. Under the R0 reading that fired, none of these numbers is a bound.

| k | A raw | B shrunk | C null bank | D median | D 5th | D 95th | D min | D max |
|---|---|---|---|---|---|---|---|---|
| 1 | 0.0926 | 0.0670 | 0.0938 | 0.0777 | 0.0197 | 0.1755 | 0.0065 | 0.2644 |
| 3 | 0.1771 | 0.1487 | 0.1861 | 0.1309 | 0.0274 | 0.2902 | 0.0113 | 0.4061 |
| 5 | 0.1901 | 0.1669 | 0.1957 | 0.1510 | 0.0299 | 0.3341 | 0.0128 | 0.4570 |
| 10 | 0.1933 | 0.1738 | 0.1958 | 0.1676 | 0.0299 | 0.3871 | 0.0128 | 0.5154 |
| 50 | 0.1933 | 0.1738 | 0.1958 | 0.1679 | 0.0299 | 0.3946 | 0.0128 | 0.5534 |
| unlimited | 0.1933 | 0.1738 | 0.1958 | 0.1679 | 0.0299 | 0.3946 | 0.0128 | 0.5534 |

Files: `results/reach2/frac_table.csv` (every set, k and method), `results/reach2/random_targets.csv` (every target, k and method).

**Note (not a pre-registered reading):** Both R0 clauses fired at every k, and either alone fires R0. The D clause does not depend on the flag's choice to compare unshrunk (C) with shrunk (B). At unlimited k, B (0.174) sits next to the D median (0.168), inside a middle 90% of 0.030-0.395. The C clause, by contrast, rests on that choice: the diagnostic (C) shrunk by B's rule closes less than B at every k (0.096 at unlimited, against 0.174); that comparison is not gated and does not change the reading.

**Note (not a pre-registered reading):** R0 is not a finding that young is unreachable. It says this bound carries no information either way. The planted-column check shows the optimiser finds a direction that is present, so the bound could register such a direction. A positive R1 would still have been sizing information only.

**Note (not a pre-registered reading):** For sets A, B and C the unlimited optimum is sparse: 9 factors for A, 10 for B and 6 for C (`frac_table.csv`). So k = 10, 50 and unlimited give the same value. In 20 dimensions, a handful of vectors is enough to reach the best point of the cone the dose cap allows.

**Note (not a pre-registered reading):** Shortlist against cell count (this project's known trap). Nine of the 10 factors in B's unlimited support sit above the 75th cell-count percentile (chance about 2.5 of 10); CDX2 at p21 is the exception. The support's median cell count is 158 against 106 overall, and Spearman rho between weight and cell count over all 1,836 factors is 0.080 (p = 0.0006). Part of this is mechanical: shrinkage cuts low-count factors hardest. Under R0 no factor list here is a shortlist.

**Note (not a pre-registered reading):** What the convergence change did. 531 gated solves end with a gap above REACH's 3.900e-06 and pass O1 (largest gap/tolerance 0.199). All 531 are intermediate L-BFGS-B solves of the L1 bisection: 526 in set D (122 of the 200 targets), 3 in B, 2 in A. Under REACH's rule each would have been restarted. They decide which support the bisection finds; they do not decide the reported frac. Every reported frac comes from a BVLS fit (greedy re-fit, L1 refit or unlimited), and the largest BVLS gap relative to max(P, 1) is 9.6e-15. No solve needed a retry. Three L-BFGS-B solves used a restart. On set D target 11, where REACH stopped, the largest gap is 1.01e-05: above REACH's tolerance, and 0.045 of that solve's O1 tolerance.

**Note (not a pre-registered reading):** Because REACH's partial values for sets A, B and C were never seen, this run's A-C values were not compared with them, deliberately. This run is not a check on REACH's partial solve.

## Solver diagnostics (O3)

Every solve was recorded whatever the outcome. A solve is one solver call; a row of `solver_diagnostics.csv` is one (space, set, target, k, method) and aggregates the solves assigned to it by the flag's rule. `final gap` is the largest gap any solve of the row ended with (after any retry); `gap/tol` is gap over that solve's own O1 tolerance 1e-6 x max(P, 1).

- Gated solves (solver check, A, B, C young, D targets): 16371. Failed O1 at the first attempt: 0. Retries used: 0. Second failures: 0.
- Largest gated gap: 5.088e-05; largest gap/tolerance ratio: 0.199.
- Counterfactual against REACH's rule (not gated; evaluated on this run's gaps, not a replay of REACH, whose restart loop stopped at different points): 531 gated solves pass O1 but have a gap above REACH's 3.900e-06; 0 failed O1 at the first attempt but would have passed REACH's test.

Gated solves by kind:

| kind | solves | max gap | max gap/tol | retries | median objective P | max L-BFGS-B restarts |
|---|---|---|---|---|---|---|
| bvls_unlimited | 203 | 1.02e-12 | 3.84e-09 | 0 | 269.0010 | 0 |
| greedy_refit | 2393 | 6.11e-13 | 7.11e-09 | 0 | 267.1761 | 0 |
| l1_refit | 1020 | 2.16e-12 | 9.62e-09 | 0 | 300.3789 | 0 |
| lbfgsb | 12755 | 5.09e-05 | 0.1988 | 0 | 279.0130 | 1 |

The 20 largest gated gaps, ranked by gap:

| set | target | k | method | kind | step | lam | P | gap | tol | ratio | retry_used | passed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | young | 1 | l1 | lbfgsb | 2 | 0.7285 | 255.8746 | 5.09e-05 | 0.0003 | 0.1988 | False | True |
| D | 21 | 50 | l1 | lbfgsb | 3 | 2.81e-07 | 243.5192 | 1.93e-05 | 0.0002 | 0.0792 | False | True |
| D | 151 | 50 | l1 | lbfgsb | 3 | 2.73e-07 | 242.0552 | 1.88e-05 | 0.0002 | 0.0779 | False | True |
| D | 103 | 50 | l1 | lbfgsb | 3 | 4.16e-07 | 223.5147 | 1.80e-05 | 0.0002 | 0.0806 | False | True |
| A | young | 1 | l1 | lbfgsb | 1 | 0.0023 | 253.8045 | 1.67e-05 | 0.0003 | 0.0660 | False | True |
| D | 118 | 50 | l1 | lbfgsb | 4 | 5.69e-08 | 286.8722 | 1.67e-05 | 0.0003 | 0.0580 | False | True |
| D | 192 | 50 | l1 | lbfgsb | 3 | 3.16e-07 | 269.0010 | 1.56e-05 | 0.0003 | 0.0578 | False | True |
| D | 105 | 10 | l1 | lbfgsb | 5 | 26.6258 | 269.0484 | 1.54e-05 | 0.0003 | 0.0571 | False | True |
| D | 1 | 50 | l1 | lbfgsb | 4 | 6.37e-08 | 211.6016 | 1.46e-05 | 0.0002 | 0.0689 | False | True |
| D | 21 | 50 | l1 | lbfgsb | 5 | 3.24e-08 | 243.5192 | 1.45e-05 | 0.0002 | 0.0594 | False | True |
| D | 192 | 50 | l1 | lbfgsb | 4 | 7.49e-08 | 269.0010 | 1.44e-05 | 0.0003 | 0.0534 | False | True |
| D | 109 | 50 | l1 | lbfgsb | 6 | 2.53e-08 | 190.2216 | 1.40e-05 | 0.0002 | 0.0736 | False | True |
| D | 26 | 50 | l1 | lbfgsb | 3 | 1.99e-07 | 284.6739 | 1.37e-05 | 0.0003 | 0.0482 | False | True |
| D | 56 | 50 | l1 | lbfgsb | 3 | 3.00e-07 | 280.2055 | 1.35e-05 | 0.0003 | 0.0480 | False | True |
| D | 63 | 50 | l1 | lbfgsb | 2 | 3.13e-06 | 272.0090 | 1.34e-05 | 0.0003 | 0.0493 | False | True |
| D | 118 | 50 | l1 | lbfgsb | 3 | 2.40e-07 | 286.8722 | 1.32e-05 | 0.0003 | 0.0462 | False | True |
| D | 1 | 50 | l1 | lbfgsb | 3 | 2.69e-07 | 211.6016 | 1.28e-05 | 0.0002 | 0.0607 | False | True |
| D | 98 | 50 | l1 | lbfgsb | 3 | 4.11e-07 | 134.7817 | 1.23e-05 | 0.0001 | 0.0911 | False | True |
| D | 152 | 50 | l1 | lbfgsb | 2 | 6.42e-06 | 139.0511 | 1.22e-05 | 0.0001 | 0.0878 | False | True |
| D | 100 | 50 | l1 | lbfgsb | 4 | 7.96e-08 | 236.6531 | 1.20e-05 | 0.0002 | 0.0509 | False | True |

The 20 largest gated gaps, ranked by gap/tol:

| set | target | k | method | kind | step | lam | P | gap | tol | ratio | retry_used | passed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | young | 1 | l1 | lbfgsb | 2 | 0.7285 | 255.8746 | 5.09e-05 | 0.0003 | 0.1988 | False | True |
| D | 123 | 50 | l1 | lbfgsb | 6 | 3.86e-08 | 77.7832 | 7.27e-06 | 7.78e-05 | 0.0935 | False | True |
| D | 98 | 50 | l1 | lbfgsb | 3 | 4.11e-07 | 134.7817 | 1.23e-05 | 0.0001 | 0.0911 | False | True |
| D | 152 | 50 | l1 | lbfgsb | 2 | 6.42e-06 | 139.0511 | 1.22e-05 | 0.0001 | 0.0878 | False | True |
| D | 103 | 50 | l1 | lbfgsb | 3 | 4.16e-07 | 223.5147 | 1.80e-05 | 0.0002 | 0.0806 | False | True |
| D | 21 | 50 | l1 | lbfgsb | 3 | 2.81e-07 | 243.5192 | 1.93e-05 | 0.0002 | 0.0792 | False | True |
| D | 151 | 50 | l1 | lbfgsb | 3 | 2.73e-07 | 242.0552 | 1.88e-05 | 0.0002 | 0.0779 | False | True |
| D | 109 | 50 | l1 | lbfgsb | 6 | 2.53e-08 | 190.2216 | 1.40e-05 | 0.0002 | 0.0736 | False | True |
| D | 123 | 50 | l1 | lbfgsb | 5 | 5.53e-08 | 77.7832 | 5.48e-06 | 7.78e-05 | 0.0704 | False | True |
| D | 137 | 50 | l1 | lbfgsb | 3 | 5.27e-07 | 168.5111 | 1.17e-05 | 0.0002 | 0.0693 | False | True |
| D | 1 | 50 | l1 | lbfgsb | 4 | 6.37e-08 | 211.6016 | 1.46e-05 | 0.0002 | 0.0689 | False | True |
| D | 123 | 50 | l1 | lbfgsb | 3 | 4.79e-07 | 77.7832 | 5.26e-06 | 7.78e-05 | 0.0676 | False | True |
| D | 123 | 50 | l1 | lbfgsb | 4 | 1.13e-07 | 77.7832 | 5.25e-06 | 7.78e-05 | 0.0675 | False | True |
| D | 175 | 10 | l1 | lbfgsb | 4 | 3.4631 | 140.1419 | 9.38e-06 | 0.0001 | 0.0669 | False | True |
| D | 77 | 1 | l1 | lbfgsb | 2 | 0.6927 | 136.0902 | 9.00e-06 | 0.0001 | 0.0661 | False | True |
| A | young | 1 | l1 | lbfgsb | 1 | 0.0023 | 253.8045 | 1.67e-05 | 0.0003 | 0.0660 | False | True |
| D | 71 | 50 | l1 | lbfgsb | 4 | 9.15e-08 | 191.0951 | 1.20e-05 | 0.0002 | 0.0629 | False | True |
| D | 193 | 50 | l1 | lbfgsb | 5 | 4.18e-08 | 185.1326 | 1.15e-05 | 0.0002 | 0.0623 | False | True |
| D | 94 | 10 | l1 | lbfgsb | 5 | 1.0555 | 194.5146 | 1.18e-05 | 0.0002 | 0.0608 | False | True |
| D | 1 | 50 | l1 | lbfgsb | 3 | 2.69e-07 | 211.6016 | 1.28e-05 | 0.0002 | 0.0607 | False | True |

Per (set, k, method), over targets (set D has 200 rows per k and method; full per-target table: `results/reach2/solver_diagnostics.csv`):

| set | k | method | rows | solves | max final gap | max gap/tol | retries | rows not converged | rows not reached |
|---|---|---|---|---|---|---|---|---|---|
| B+planted | 1 | greedy | 1 | 1 | 9.21e-30 | 9.21e-24 | 0 | 0 | 0 |
| B+planted | 1 | l1 | 1 | 2 | 9.21e-30 | 9.21e-24 | 0 | 0 | 0 |
| B+planted | 3 | greedy | 1 | 0 |  |  | 0 | 0 | 0 |
| B+planted | 3 | l1 | 1 | 54 | 4.68e-08 | 0.0468 | 0 | 0 | 0 |
| B+planted | 5 | greedy | 1 | 0 |  |  | 0 | 0 | 0 |
| B+planted | 5 | l1 | 1 | 1 | 9.21e-30 | 9.21e-24 | 0 | 0 | 0 |
| B+planted | 10 | greedy | 1 | 0 |  |  | 0 | 0 | 0 |
| B+planted | 10 | l1 | 1 | 1 | 9.21e-30 | 9.21e-24 | 0 | 0 | 0 |
| B+planted | 50 | greedy | 1 | 0 |  |  | 0 | 0 | 0 |
| B+planted | 50 | l1 | 1 | 1 | 9.21e-30 | 9.21e-24 | 0 | 0 | 0 |
| A | 1 | greedy | 1 | 1 | 1.14e-13 | 3.54e-10 | 0 | 0 | 0 |
| A | 1 | l1 | 1 | 7 | 5.09e-05 | 0.1988 | 0 | 0 | 0 |
| A | 3 | greedy | 1 | 2 | 5.68e-14 | 2.15e-10 | 0 | 0 | 0 |
| A | 3 | l1 | 1 | 1 | 5.68e-14 | 1.82e-10 | 0 | 0 | 0 |
| A | 5 | greedy | 1 | 2 | 1.14e-13 | 4.44e-10 | 0 | 0 | 0 |
| A | 5 | l1 | 1 | 1 | 0.0000 | 0.0000 | 0 | 0 | 0 |
| A | 10 | greedy | 1 | 4 | 3.13e-13 | 1.23e-09 | 0 | 0 | 0 |
| A | 10 | l1 | 1 | 54 | 1.51e-06 | 0.0060 | 0 | 0 | 0 |
| A | 50 | greedy | 1 | 0 |  |  | 0 | 0 | 0 |
| A | 50 | l1 | 1 | 1 | 8.53e-14 | 3.36e-10 | 0 | 0 | 0 |
| A | unlimited | bvls | 1 | 1 | 8.53e-14 | 3.36e-10 | 0 | 0 | 0 |
| B | 1 | greedy | 1 | 1 | 1.14e-13 | 3.35e-10 | 0 | 0 | 0 |
| B | 1 | l1 | 1 | 8 | 3.84e-06 | 0.0143 | 0 | 0 | 0 |
| B | 3 | greedy | 1 | 2 | 5.68e-14 | 1.91e-10 | 0 | 0 | 0 |
| B | 3 | l1 | 1 | 3 | 2.60e-07 | 0.0007 | 0 | 0 | 0 |
| B | 5 | greedy | 1 | 2 | 0.0000 | 0.0000 | 0 | 0 | 0 |
| B | 5 | l1 | 1 | 1 | 0.0000 | 0.0000 | 0 | 0 | 0 |
| B | 10 | greedy | 1 | 5 | 2.27e-13 | 8.53e-10 | 0 | 0 | 0 |
| B | 10 | l1 | 1 | 1 | 8.53e-13 | 3.20e-09 | 0 | 0 | 0 |
| B | 50 | greedy | 1 | 1 | 2.27e-13 | 8.54e-10 | 0 | 0 | 0 |
| B | 50 | l1 | 1 | 52 | 7.28e-06 | 0.0273 | 0 | 0 | 0 |
| B | unlimited | bvls | 1 | 1 | 1.02e-12 | 3.84e-09 | 0 | 0 | 0 |
| C | 1 | greedy | 1 | 1 | -1.71e-13 | -5.32e-10 | 0 | 0 | 0 |
| C | 1 | l1 | 1 | 7 | 6.02e-08 | 0.0002 | 0 | 0 | 0 |
| C | 3 | greedy | 1 | 2 | 5.68e-14 | 2.20e-10 | 0 | 0 | 0 |
| C | 3 | l1 | 1 | 1 | 5.68e-14 | 1.75e-10 | 0 | 0 | 0 |
| C | 5 | greedy | 1 | 2 | -5.68e-14 | -2.25e-10 | 0 | 0 | 0 |
| C | 5 | l1 | 1 | 4 | 9.27e-12 | 3.08e-08 | 0 | 0 | 0 |
| C | 10 | greedy | 1 | 1 | 2.84e-14 | 1.13e-10 | 0 | 0 | 0 |
| C | 10 | l1 | 1 | 54 | 1.98e-06 | 0.0079 | 0 | 0 | 0 |
| C | 50 | greedy | 1 | 0 |  |  | 0 | 0 | 0 |
| C | 50 | l1 | 1 | 1 | 2.84e-14 | 1.13e-10 | 0 | 0 | 0 |
| C | unlimited | bvls | 1 | 1 | 1.42e-13 | 5.63e-10 | 0 | 0 | 0 |
| D | 1 | greedy | 200 | 200 | 2.84e-13 | 9.08e-10 | 0 | 0 | 0 |
| D | 1 | l1 | 200 | 1537 | 1.17e-05 | 0.0661 | 0 | 0 | 0 |
| D | 3 | greedy | 200 | 398 | 3.41e-13 | 1.03e-09 | 0 | 0 | 0 |
| D | 3 | l1 | 200 | 470 | 4.13e-06 | 0.0125 | 0 | 0 | 0 |
| D | 5 | greedy | 200 | 389 | 5.68e-13 | 1.89e-09 | 0 | 0 | 0 |
| D | 5 | l1 | 200 | 800 | 5.26e-06 | 0.0154 | 0 | 0 | 0 |
| D | 10 | greedy | 200 | 813 | 4.55e-13 | 1.91e-09 | 0 | 0 | 0 |
| D | 10 | l1 | 200 | 3570 | 1.54e-05 | 0.0669 | 0 | 0 | 0 |
| D | 50 | greedy | 200 | 566 | 6.11e-13 | 7.11e-09 | 0 | 0 | 0 |
| D | 50 | l1 | 200 | 7143 | 1.93e-05 | 0.0935 | 0 | 0 | 0 |
| D | unlimited | bvls | 200 | 200 | 5.97e-13 | 3.41e-09 | 0 | 0 | 0 |
| Cs | 1 | greedy | 1 | 1 | 5.68e-14 | 1.59e-10 | 0 | 0 | 0 |
| Cs | 1 | l1 | 1 | 10 | 1.73e-07 | 0.0005 | 0 | 0 | 0 |
| Cs | 3 | greedy | 1 | 2 | 1.71e-13 | 5.05e-10 | 0 | 0 | 0 |
| Cs | 3 | l1 | 1 | 1 | 2.84e-13 | 8.29e-10 | 0 | 0 | 0 |
| Cs | 5 | greedy | 1 | 2 | 1.14e-13 | 3.50e-10 | 0 | 0 | 0 |
| Cs | 5 | l1 | 1 | 4 | 1.12e-07 | 0.0003 | 0 | 0 | 0 |
| Cs | 10 | greedy | 1 | 5 | 1.71e-13 | 5.33e-10 | 0 | 0 | 0 |
| Cs | 10 | l1 | 1 | 4 | 1.39e-06 | 0.0042 | 0 | 0 | 0 |
| Cs | 50 | greedy | 1 | 6 | 2.27e-13 | 7.14e-10 | 0 | 0 | 0 |
| Cs | 50 | l1 | 1 | 52 | 1.95e-06 | 0.0061 | 0 | 0 | 0 |
| Cs | unlimited | bvls | 1 | 1 | 2.27e-13 | 7.14e-10 | 0 | 0 | 0 |

Full gene space S1 (diagnostic, not gated):

| set | solves | max_final_gap | max_ratio | retries | not_converged |
|---|---|---|---|---|---|
| A | 74 | 3.02e-05 | 0.0580 | 0 | 0 |
| B | 77 | 6.30e-06 | 0.0122 | 0 | 0 |
| C | 79 | 3.79e-05 | 0.0734 | 0 | 0 |

## Greedy versus L1 at finite k

Both methods are feasible solutions; neither is guaranteed optimal for k > 1 (greedy is exhaustive at k = 1). `L1 nnz hit` is the support size the bisection reached before the BVLS refit; "not exact" means no lam in 60 steps gave exactly k and the fallback rule was used.

| set | k | greedy frac | greedy nnz | L1 frac | L1 nnz | L1 nnz hit | abs diff | support Jaccard |
|---|---|---|---|---|---|---|---|---|
| A raw | 1 | 0.0926 | 1 | 0.0670 | 1 | 1 | 0.0256 | 0.0000 |
| A raw | 3 | 0.1771 | 3 | 0.1043 | 3 | 3 | 0.0728 | 0.2000 |
| A raw | 5 | 0.1901 | 5 | 0.1267 | 4 | 5 | 0.0634 | 0.1250 |
| A raw | 10 | 0.1933 | 9 | 0.1933 | 9 | 9 (not exact) | 0.0000 | 1.0000 |
| A raw | 50 | 0.1933 | 9 | 0.1933 | 9 | 9 (not exact) | 0.0000 | 1.0000 |
| B shrunk | 1 | 0.0670 | 1 | 0.0589 | 1 | 1 | 0.0081 | 0.0000 |
| B shrunk | 3 | 0.1487 | 3 | 0.1173 | 3 | 3 | 0.0314 | 0.2000 |
| B shrunk | 5 | 0.1669 | 5 | 0.1658 | 5 | 5 | 0.0011 | 0.6667 |
| B shrunk | 10 | 0.1738 | 10 | 0.1738 | 10 | 10 | 0.0000 | 1.0000 |
| B shrunk | 50 | 0.1738 | 10 | 0.1738 | 10 | 10 (not exact) | 0.0000 | 1.0000 |
| C null bank | 1 | 0.0938 | 1 | 0.0141 | 1 | 1 | 0.0797 | 0.0000 |
| C null bank | 3 | 0.1861 | 3 | 0.0871 | 3 | 3 | 0.0990 | 0.0000 |
| C null bank | 5 | 0.1957 | 5 | 0.1889 | 4 | 5 | 0.0068 | 0.5000 |
| C null bank | 10 | 0.1958 | 6 | 0.1958 | 6 | 6 (not exact) | 0.0000 | 1.0000 |
| C null bank | 50 | 0.1958 | 6 | 0.1958 | 6 | 6 (not exact) | 0.0000 | 1.0000 |
| C null bank, shrunk (diagnostic) | 1 | 0.0415 | 1 | 0.0246 | 1 | 1 | 0.0169 | 0.0000 |
| C null bank, shrunk (diagnostic) | 3 | 0.0796 | 3 | 0.0624 | 3 | 3 | 0.0172 | 0.5000 |
| C null bank, shrunk (diagnostic) | 5 | 0.0879 | 5 | 0.0864 | 5 | 5 | 0.0015 | 0.6667 |
| C null bank, shrunk (diagnostic) | 10 | 0.0951 | 10 | 0.0937 | 9 | 10 | 0.0013 | 0.5833 |
| C null bank, shrunk (diagnostic) | 50 | 0.0962 | 15 | 0.0962 | 15 | 15 (not exact) | 0.0000 | 1.0000 |

Set (D), over the 200 targets:

| k | median abs diff | max abs diff | L1 better (targets) | greedy better (targets) | L1 did not hit k exactly |
|---|---|---|---|---|---|
| 1 | 0.0099 | 0.0928 | 0 | 134 | 0 |
| 3 | 0.0180 | 0.1720 | 4 | 176 | 2 |
| 5 | 0.0090 | 0.1565 | 12 | 160 | 7 |
| 10 | 1.94e-16 | 0.0627 | 12 | 82 | 67 |
| 50 | 0.0000 | 2.91e-16 | 0 | 0 | 200 |

## Diagnostics (not gated)

### What shrinkage costs: (A) beside (B)

| k | A raw | B shrunk |
|---|---|---|
| 1 | 0.0926 | 0.0670 |
| 3 | 0.1771 | 0.1487 |
| 5 | 0.1901 | 0.1669 |
| 10 | 0.1933 | 0.1738 |
| 50 | 0.1933 | 0.1738 |
| unlimited | 0.1933 | 0.1738 |

Factors shrunk to zero in (B) (||D20_f||^2 <= sigma^2(n_f)): 575 of 1836. Chance baseline: the same rule shrinks 1011 of the 1836 null-bank vectors (C) to zero. Median s_f 0.130; median ||D20_f|| 11.456 against median sigma(n_f) 9.699.

### (C) shrunk by the same rule

| k | C shrunk (diag) | C null bank | B shrunk |
|---|---|---|---|
| 1 | 0.0415 | 0.0938 | 0.0670 |
| 3 | 0.0796 | 0.1861 | 0.1487 |
| 5 | 0.0879 | 0.1957 | 0.1669 |
| 10 | 0.0951 | 0.1958 | 0.1738 |
| 50 | 0.0962 | 0.1958 | 0.1738 |
| unlimited | 0.0962 | 0.1958 | 0.1738 |

### Which factors carry the weight, and cell count

This project's shortlist-versus-depth check. Cell-count percentile is among the 1,836 factors (median cell count 106).

- B shrunk, k=1, greedy: DMRT1 (w 0.90, 62 cells, p6)
- B shrunk, k=1, l1: HELZ2 (w 0.69, 140 cells, p79)
- B shrunk, k=3, greedy: BNC1 (w 1.00, 152 cells, p85), HELZ2 (w 0.71, 140 cells, p79), GFI1 (w 0.69, 223 cells, p98)
- B shrunk, k=3, l1: MYC (w 1.00, 99 cells, p42), HELZ2 (w 0.74, 140 cells, p79), DMRT1 (w 0.06, 62 cells, p6)
- B shrunk, k=5, greedy: RELA (w 1.00, 179 cells, p93), MAFA (w 1.00, 164 cells, p89), BNC1 (w 0.92, 152 cells, p85), HELZ2 (w 0.69, 140 cells, p79), GFI1 (w 0.44, 223 cells, p98)
- B shrunk, k=5, l1: BNC1 (w 1.00, 152 cells, p85), RELA (w 1.00, 179 cells, p93), HELZ2 (w 0.71, 140 cells, p79), GFI1 (w 0.54, 223 cells, p98), FOXF1 (w 0.25, 137 cells, p78)
- B shrunk, k=10, greedy: MAFA (w 1.00, 164 cells, p89), DMRTA1 (w 1.00, 183 cells, p94), RELA (w 0.99, 179 cells, p93), HELZ2 (w 0.66, 140 cells, p79), TLE1 (w 0.60, 149 cells, p84), BNC1 (w 0.53, 152 cells, p85), GSC2 (w 0.43, 166 cells, p90), GFI1 (w 0.36, 223 cells, p98), FOXF1 (w 0.29, 137 cells, p78), CDX2 (w 0.12, 80 cells, p21)
- B shrunk, k=10, l1: MAFA (w 1.00, 164 cells, p89), DMRTA1 (w 1.00, 183 cells, p94), RELA (w 0.99, 179 cells, p93), HELZ2 (w 0.66, 140 cells, p79), TLE1 (w 0.60, 149 cells, p84), BNC1 (w 0.53, 152 cells, p85), GSC2 (w 0.43, 166 cells, p90), GFI1 (w 0.36, 223 cells, p98), FOXF1 (w 0.29, 137 cells, p78), CDX2 (w 0.12, 80 cells, p21)
- A raw, k=1, greedy: DMRTA1 (w 1.00, 183 cells, p94)
- A raw, k=1, l1: DMRT1 (w 0.44, 62 cells, p6)
- A raw, k=3, greedy: DMRTA1 (w 1.00, 183 cells, p94), FOXN3 (w 0.96, 128 cells, p71), HELZ2 (w 0.37, 140 cells, p79)
- A raw, k=3, l1: RUNX1 (w 0.32, 54 cells, p3), HELZ2 (w 0.29, 140 cells, p79), DMRT1 (w 0.19, 62 cells, p6)
- A raw, k=5, greedy: DMRTA1 (w 1.00, 183 cells, p94), FOXN3 (w 0.76, 128 cells, p71), TLE1 (w 0.42, 149 cells, p84), HELZ2 (w 0.39, 140 cells, p79), GFI1 (w 0.17, 223 cells, p98)
- A raw, k=5, l1: HELZ2 (w 0.41, 140 cells, p79), FUBP1 (w 0.38, 94 cells, p36), CDX2 (w 0.32, 80 cells, p21), RUNX1 (w 0.02, 54 cells, p3)
- A raw, k=10, greedy: DMRTA1 (w 1.00, 183 cells, p94), FOXN3 (w 0.58, 128 cells, p71), HELZ2 (w 0.42, 140 cells, p79), TLE1 (w 0.26, 149 cells, p84), RELA (w 0.25, 179 cells, p93), GFI1 (w 0.15, 223 cells, p98), FOXF1 (w 0.11, 137 cells, p78), MAFA (w 0.07, 164 cells, p89), GSC2 (w 0.02, 166 cells, p90)
- A raw, k=10, l1: DMRTA1 (w 1.00, 183 cells, p94), FOXN3 (w 0.58, 128 cells, p71), HELZ2 (w 0.42, 140 cells, p79), TLE1 (w 0.26, 149 cells, p84), RELA (w 0.25, 179 cells, p93), GFI1 (w 0.15, 223 cells, p98), FOXF1 (w 0.11, 137 cells, p78), MAFA (w 0.07, 164 cells, p89), GSC2 (w 0.02, 166 cells, p90)
- B shrunk, k=50, greedy: 10 factors; median cell count 158 (all factors 106); 2 at the dose cap w = 1.
- B shrunk, k=50, l1: 10 factors; median cell count 158 (all factors 106); 2 at the dose cap w = 1.
- B shrunk, k=unlimited, bvls: 10 factors; median cell count 158 (all factors 106); 2 at the dose cap w = 1.
- A raw, k=50, greedy: 9 factors; median cell count 164 (all factors 106); 1 at the dose cap w = 1.
- A raw, k=50, l1: 9 factors; median cell count 164 (all factors 106); 1 at the dose cap w = 1.
- A raw, k=unlimited, bvls: 9 factors; median cell count 164 (all factors 106); 1 at the dose cap w = 1.

Full lists: `results/reach2/weights.csv`.

| set | k | method | over | n | Spearman rho (weight vs cells) | p | median cells, selected |
|---|---|---|---|---|---|---|---|
| A | 1 | greedy | selected | 1 |  |  | 183.0000 |
| A | 1 | l1 | selected | 1 |  |  | 62.0000 |
| A | 3 | greedy | selected | 3 | 0.5000 | 0.6667 | 140.0000 |
| A | 3 | l1 | selected | 3 | -0.5000 | 0.6667 | 62.0000 |
| A | 5 | greedy | selected | 5 | -0.3000 | 0.6238 | 149.0000 |
| A | 5 | l1 | selected | 4 | 1.0000 | 0.0000 | 87.0000 |
| A | 10 | greedy | selected | 9 | -0.1167 | 0.7650 | 164.0000 |
| A | 10 | l1 | selected | 9 | -0.1167 | 0.7650 | 164.0000 |
| A | 50 | greedy | selected | 9 | -0.1167 | 0.7650 | 164.0000 |
| A | 50 | l1 | selected | 9 | -0.1167 | 0.7650 | 164.0000 |
| A | unlimited | bvls | selected | 9 | -0.1167 | 0.7650 | 164.0000 |
| A | unlimited | bvls | all 1836 (zeros included) | 1836 | 0.0880 | 0.0002 |  |
| B | 1 | greedy | selected | 1 |  |  | 62.0000 |
| B | 1 | l1 | selected | 1 |  |  | 140.0000 |
| B | 3 | greedy | selected | 3 | -0.5000 | 0.6667 | 152.0000 |
| B | 3 | l1 | selected | 3 | 0.5000 | 0.6667 | 99.0000 |
| B | 5 | greedy | selected | 5 | -0.0513 | 0.9347 | 164.0000 |
| B | 5 | l1 | selected | 5 | 0.3591 | 0.5528 | 152.0000 |
| B | 10 | greedy | selected | 10 | 0.4499 | 0.1921 | 158.0000 |
| B | 10 | l1 | selected | 10 | 0.4499 | 0.1921 | 158.0000 |
| B | 50 | greedy | selected | 10 | 0.4499 | 0.1921 | 158.0000 |
| B | 50 | l1 | selected | 10 | 0.4499 | 0.1921 | 158.0000 |
| B | unlimited | bvls | selected | 10 | 0.4499 | 0.1921 | 158.0000 |
| B | unlimited | bvls | all 1836 (zeros included) | 1836 | 0.0799 | 0.0006 |  |

### Full gene space S1 (caveat)

All 23,485 ruler coordinates, b = c_young - c_old, Gram-reduced form. S1 includes coordinates no factor can move (unmatched and missing-mask genes, z = 0 in every d_f) and many noise dimensions; the numbers are not comparable to S2A_k20 and are not gated. (D) is not run here. (B) uses sigma^2 from the same draws as full-space squared norms; 537 of 1,836 factors shrink to zero there.

| k | A raw | B shrunk | C null bank |
|---|---|---|---|
| 1 | 0.0008 | 0.0008 | 0.0005 |
| 3 | 0.0020 | 0.0020 | 0.0011 |
| 5 | 0.0029 | 0.0029 | 0.0016 |
| 10 | 0.0043 | 0.0042 | 0.0025 |
| 50 | 0.0061 | 0.0059 | 0.0040 |
| unlimited | 0.0061 | 0.0059 | 0.0041 |

Solver status in S1: A: converged True, Gram rank 1836; B: converged True, Gram rank 1299; C: converged True, Gram rank 1836.

## Verify (re-solve from tracked inputs)

`python src/reach2_run.py verify` re-solved sets A-D from `results/reach2/inputs_k20.npz`: reading R0 (stored R0); mismatches beyond 1e-4: none; largest |difference| over the D table 9.71e-17; 16311 solves, 0 retries, max gap/tol 0.199.

## What this does not cover

- Southard is Hs27 neonatal foreskin fibroblasts, not aged donors. Nothing here transfers to aged cells.
- The transfer check (FINDINGS_TRANSFER.md, commit 7230271) could not measure whether these directions carry to another fibroblast context.
- Doses above the delivered one are not modelled; any w in (0, 1) assumes a partial dose acts as a scaled displacement, itself an assumption. Doses above 1 would be extrapolation.
- Additivity of displacements is assumed, not tested.
- Nothing here models toxicity, identity loss or proliferation.
- Circularity, disclosed: the null draws, sigma^2 and (C) all come from the same 14,013 non-targeting cells that form d_f's reference pool; (C)'s vectors are among the draws that estimate sigma^2. The S2A basis was fit on GTEx donors aged 40-59 with no anchor donor in the fit (SOUTH7), so the space is not fit to the factors or the anchors.
- The relative convergence test certifies each solve against its own optimum to the stated accuracy. It does not make greedy or L1 optimal for k > 1: both remain lower estimates of the k-sparse optimum, as REACH said.

## Files

`src/reach2_run.py`; `results/reach2/`: `PREREG.flag`, `diagnostics.json`, `factors.csv`, `frac_table.csv`, `inputs_k20.npz`, `null_depths.csv`, `random_targets.csv`, `readings.json`, `run.log`, `solver_calls.csv.gz`, `solver_check.json`, `solver_diagnostics.csv`, `solver_summary.json`, `spearman.csv`, `stage0_norm_sample.csv`, `stage0_reproduction.csv`, `stage0_vs_reach.csv`, `verify.json`, `verify.log`, `weights.csv`, `worst_gaps.csv`.

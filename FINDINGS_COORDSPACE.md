# FINDINGS_COORDSPACE — coordinate-space robustness check (SECONDARY)

**Status:** DONE. Decision-table rows fired: [1].

Reproduced by `src/coordspace_run.py`. Imports `src/same_run.py`, `src/south6_run.py`, `src/toward_run.py`, `src/seng_run.py` and `src/gtex_common.py` without modifying them. The frozen ruler is not refit. Space C is the pre-registered primary and is not replaced by anything below.

## Decision-table outcome

Rows fired: 1. Surviving (qualified) spaces: A, B, C. Disqualified: none.

- Row 1: In every coordinate space that detected known young-cell mixtures (A, B, C), aged partially reprogrammed cells end farther from the young donor than they started (A: R_S = 2.054, top-10 genes carry 77.3% of the squared final distance; B: R_S = 1.879, top-10 genes carry 3.6% of the squared final distance; C: R_S = 2.002, top-10 genes carry 1.8% of the squared final distance). The same-platform result does not depend on the log or the z-score.
- Screen, space A: CHANGED. toward_young = 2 (ZNF441;ZNF669) against 0 in space C. Reported immediately; SOUTH6 space C remains the screen's primary.
- Screen, space B: UNCHANGED. toward_young = 0, the same set as space C.

## Pre-registration (verbatim, written before any COORDSPACE statistic)

```
TASK: COORDSPACE — coordinate-space robustness check. SECONDARY analysis.

Written **before any COORDSPACE statistic**, 2026-10-03T10:07:01+09:00. No threshold in this block is re-tuned after numbers exist.

Repository state when written: commit 78cfee67c4c5996e8ccb45e6339135aaf4b7b6d7
("state before coordinate-space task"). No number in space A or space B exists,
and space C has not been recomputed by this task. Every space C value quoted
below is an existing SAME or SOUTH6 result, read from results/same/ or
results/south6/ as stored.

Nature of this task: computational reanalysis of public data already on disk.
No new data, no lab work.

--- QUESTION ---
Every distance in this project is computed in one coordinate system:
TMM -> log2 CPM (prior.count=2) -> z-score against the frozen GTEx donor mean
and SD. Two of those steps are untested choices. The log compresses a
million-fold abundance range to near parity. The z-score divides by each
gene's donor-to-donor spread, so a gene that barely varies between people can
give a large score from a tiny absolute change. Do the conclusions of SAME and
SOUTH6 survive the choice of coordinates?

Which spaces count is decided by controls, not by argument: a space counts only
if, in that space, the mixture control detects known young-cell fractions. No
view on which coordinates are philosophically right enters any reading.

--- STATUS ---
- SECONDARY. Space C is the pre-registered primary of SAME and of SOUTH6 and
  stays primary whatever this task finds. Nothing here replaces a SAME or
  SOUTH6 number.
- The preprint is submitted to bioRxiv and will not be withdrawn. The only
  consequences this task can have are those in the DECISION TABLE below.
- The frozen ruler is NOT refit in any space. It stays exactly as it is
  (results/fibro/frozen_ruler_ridge_raw.npz). Ridge needs z, and it is frozen.
  Only the distance test moves: centroids, displacement, distance. That is
  geometry and needs no regression.

--- EXECUTION RULES (binding) ---
- src/coordspace_run.py runs standalone (`python src/coordspace_run.py`) and
  refuses to compute if this flag is missing.
- On a STOP, write the findings and exit. Never compute after a STOP; never
  resume a stopped run.
- No checkpoint/resume logic unless resumed random draws are verified
  bit-identical to an uninterrupted stream, and that check is reported.
- FINDINGS and PROGRESS report the same stop status; every failure printed
  verbatim. Do not substitute columns or repair rows.
- Write FINDINGS_COORDSPACE.md, PROGRESS_COORDSPACE.md, results/coordspace/*,
  src/coordspace_run.py. Modify no other FINDINGS_*, PROGRESS_* or src/ file.
  Import from src/same_run.py, src/south6_run.py, src/toward_run.py and
  src/gtex_common.py without changing them.
- This file is reproduced verbatim at the top of FINDINGS_COORDSPACE.md.
- Seeds: 20260914 (splits, permutations, mixture draws), 20260918
  (bootstraps), used exactly as SAME and SOUTH6 use them.

--- THE THREE SPACES (everything else held fixed) ---
Only the last step of the transform changes. Everything upstream is shared and
built once: the pseudobulk panels, the summed counts, the TMM factors
(estimated on each analysis's own panel exactly as now), the SAME half-split,
the mixture draws, and the bootstrap and permutation index streams.

For panel row i with summed counts C_ij over gene j, TMM factor nf_i, and
effective library size L_i = nf_i * sum_j C_ij:

  A. CPM, no log, no z      x_ij = 1e6 * C_ij / L_i
                            No prior count: it exists only to stabilise the
                            log. Magnitude dominates.
  B. log2 CPM, no z         x_ij = log2_cpm_edger(C, nf, prior=2), the exact
                            matrix the current pipeline z-scores.
                            Fold-change space.
  C. log2 CPM + frozen z    x_ij = (B_ij - mu_j) / sd_j, sd_j < 1e-12 set to 1.
                            The current pipeline, unchanged. Pre-registered
                            primary.

Genes with no counts in any row of a panel, and ruler genes absent from a
dataset, are set to one constant per gene across every row of that panel: 0 in
A, mu_j in B (the B image of z = 0), 0 in C (the current rule). Their
contribution to every within-panel difference is therefore exactly 0 in every
space, as it is now.

Same gene universe in every space: the 23,485 frozen-ruler genes (SAME), and
the same 20,265 matched genes inside them (SOUTH6, mapping
best_after_md3_idtype, unchanged).

No clipping, flooring or re-centring in any space. In A, origin + displacement
can have negative entries (a CPM displacement measured in Hs27 added to a GTEx
CPM origin). They are kept as they are and counted in the findings.

--- GTEX CENTROIDS PER SPACE (SOUTH6 O1 and Y1) ---
Same donors in every space: the YOUNG bins (20-29, 30-39) and OLD bins
(60-69, 70-79) frozen in results/toward/anchors.npz. A centroid is the
arithmetic mean of the member donors' coordinates in that space, as in C.
  C: c_young and c_old from results/toward/anchors.npz, used as stored.
  B: mean of the same donors' columns of the GTEx pack log2-CPM
     (data/processed/fibro/stage1_fibro.npz, key `logcpm`, TMM over all 652
     donors as stored, not re-estimated).
  A: mean of the same donors' plain TMM-CPM, from the pack's `counts` and
     `tmm_factors`: 1e6 * count / (tmm_factor * column sum). Not
     back-transformed from B: the mean of 2^x is not 2 to the mean of x.
SAME uses no GTEx centroid in its statistics. Its context column
cos(v, u_TOWARD) uses the young-minus-old unit vector of that space.

--- STAGE 0: REPRODUCE SPACE C FIRST (gate, before any A or B number) ---
Rebuild space C through this task's code and compare it with the stored
results:
  (a) SAME: results/same/posctrl.csv (every row), results/same/stats.csv
      (forward and reverse rows), results/same/summary.json (asymmetry and its
      CI), results/same/context.csv, and the key `no_approach_same_platform`
      with delta_S = +130.894. The comparison is against results/same/ as
      stored, after the DEBUG_SAME bootstrap fix, not against the pre-fix
      intervals printed in FINDINGS_SAME.md.
  (b) SOUTH6, setting S1 O1_Y1 only: ||O1 - Y1|| (origin_target_distances.csv),
      the per-factor Stage 2A floors (stage2a_floors.npz), the Stage 2C MDA
      (stage2c_mda.csv), and per factor delta, cos, frac, p, q and the
      toward_young label (stage3_full_table.csv).
  (c) GTEx centroids: the B centroids against mu + sd * (C centroids).
Pass requires every categorical outcome identical (keys, P pass, and per factor
the below-floor, q <= 0.05 and toward_young indicators) and every compared
number within |new - stored| <= tol * max(1, |stored|), with tol = 1e-6 for
SAME and the centroids and tol = 1e-4 for SOUTH6 (its permutation pool is
stored in float32).
A pass also shows that A and B will see the same permutation and bootstrap
index streams that produced the stored space C results.
Fail -> STOP `space_c_not_reproduced`, naming every mismatch. No A or B number
is computed.

--- STAGE 1: MIXTURE CONTROL IN EVERY SPACE (before any PartialReprog vector in A or B) ---
SAME positive control P, procedure unchanged, in A and B (C comes from
Stage 0):
  half-B d0 Fibroblast cells only, split as stored in
  results/same/split_idx.npz (not re-drawn); young fractions
  f = 0, 0.10, 0.25, 0.50; cell counts, not weights; total cells fixed = min
  available; seed 20260914; each mixture scored as S in the forward test
  against the half-A anchors O and Y built in that space.
  P PASSES in a space only if frac is monotonically increasing in f, AND frac
  at f = 0.50 has a bootstrap CI excluding 0, AND delta at f = 0.50 < 0.
  Report the noise floor (frac and delta at f = 0) in every space. Never
  subtract it.
Interval check, added here and applied to every space: if the frac CI at
f = 0.50 does not contain its own point estimate, that condition is
`ci_invalid` and counts as not met. A control that cannot be verified does not
qualify a space. (In space C as stored this interval contains its estimate.)

QUALIFICATION: a space is QUALIFIED if and only if P passes in it. Otherwise
it is DISQUALIFIED (key `disqualified_mixture`). Disqualification applies to
the whole space, SAME and SOUTH6 halves alike.

Deviation from SAME, stated now: in SAME a failed P stopped everything
downstream. Here a disqualified space is still run to the end, because the
decision table needs to know whether a disqualified space flips. Its numbers
go in a separate block headed DISQUALIFIED, and no reading or conclusion is
drawn from them.

After Stage 1, before any PartialReprog vector in A or B, write
results/coordspace/QUALIFIED_BEFORE_PRIMARY.flag. It holds this flag verbatim
plus, per space, P pass or fail and each of its three conditions. It is not
rewritten afterwards.

--- STAGE 2: SAME IN EVERY SPACE ---
The SAME pre-registration (results/same/PREREG.flag) applied with x in place of
z, using src/same_run.py as it stands at commit 78cfee6 (after the DEBUG_SAME
bootstrap fix):
  v = x(Y) - x(O); d_S = x(S) - x(O), S = GM00731 PartialReprog pooled over
  timepoints; cos_S; frac_S = dot(d_S, v) / ||v||^2;
  delta_S = ||x(S) - x(Y)|| - ||x(O) - x(Y)||.
Reversed-roles control, run in every space including disqualified ones: the
SAME reverse test, v_rev = x(O) - x(Y), origin Y, S_rev = GM23815
PartialReprog pooled, giving cos_rev, frac_rev and delta_rev.
Asymmetry = frac_S - frac_rev, with its CI taken from the same bootstrap draw.
The reversed-roles control does not qualify or disqualify a space. It enters
the SAME reading exactly as in SAME: it is what separates
`toward_young_same_platform` from `convergence_not_rejuvenation`.
Nulls and intervals: N1 = 200 permutations of v (and v_rev) across genes,
seed 20260914, p = (#{cos >= observed} + 1) / 201. Bootstrap B = 200, seed
20260918.
Context, reported only: GM00731 Pluripotency and GM00731 NonReprog as S;
cos(v, u_TOWARD) with u from that space's GTEx centroids. The frozen-ruler
score is reported in space C only, because the ruler needs z.
Unit-free summaries in every space: frac_S, and
R_S = ||x(S) - x(Y)|| / ||x(O) - x(Y)||.
Reading per space: SAME's decision order 1-5 and its `overshoot_also_closer`
flag, verbatim. Interval check as in Stage 1: if a gating interval (frac_S,
Asymmetry, Pluripotency frac) does not contain its point estimate, its
condition is `ci_invalid`, and any key whose firing depends on it becomes
`mixed`, naming the interval.

SIGN AND OUTCOME PER SPACE (SAME):
  sign = sign of the delta_S point estimate. SAME's own reading uses the point
  estimate; the delta interval is not a gate.
  HOLDS: delta_S >= 0 (SAME key `no_approach_same_platform`, as in C).
  FLIPS: delta_S < 0, whatever key fires. The key is reported with it.
Space C as stored: delta_S = +130.894, HOLDS (re-confirmed in Stage 0).

--- STAGE 3: SOUTH6 IN EVERY SPACE ---
Same factors (all 1,836 targets, none dropped), same control set (all 78
non-targeting guides, guide_target == 'non'), same panel (1,836 factor rows +
78 guide rows + 1 pooled non-targeting row, TMM over that panel), same gates,
and the same pre-registered primary comparison: O1 -> Y1 in S1 gene space.
  d_g = x(target g) - x(all non-targeting cells pooled), in that space.
  Stage 2A floor: each factor's own vector permuted across the 20,265 matched
  coordinates, n_perm = 5000, seed 20260914, as the same uninterrupted index
  stream SOUTH6 used. delta = ||(O1 + d_perm) - Y1|| - ||O1 - Y1||; the
  per-factor floor is the 5th percentile; the pooled floor is taken at that
  space's median factor norm. Every space recomputes its floors by the same
  procedure, because floors are in that space's units. A space C floor is
  never applied to A or B.
  Stage 2B control floor (78 guides, and random splits of the non-targeting
  cells): reported alongside, as in SOUTH6.
  Stage 2C dose-response: that space's Y1 - O1 added to O1 at
  f = 0, 0.05, 0.1, 0.25, 0.5, 1.0, scaled to that space's median factor
  norm; MDA = smallest f clearing the primary floor. If no f <= 1.0 clears it
  in a space, that space gets key `no_detection_power` and no ranking is read
  out there. This applies to the screen half only and does not disqualify the
  space for SAME.
  Stage 3: per factor cos, frac, delta; one-sided permutation p from Stage 2A;
  BH across the 1,836 factors within the space.
  toward_young only if delta < 0 AND delta is below that factor's own Stage 2A
  floor in that space AND q <= 0.05 AND not identity_loss AND not
  guides_disagree.
  identity_loss and guides_disagree are carried unchanged from SOUTH6 space C
  (stage3_full_table.csv, S1 O1_Y1) and are not recomputed. They are
  factor-level quality gates, not part of the distance test, and the 0.5
  identity_loss threshold is defined in frozen-z units with no counterpart in A
  or B.
  Report per space: ||O1 - Y1||; per-factor floor min, median and max; the
  pooled floor; MDA; the counts of delta < 0, below floor, q <= 0.05 and
  toward_young; and each toward_young factor by name with its delta, floor, q
  and n_cells. Unit-free: median per-factor floor / ||O1 - Y1||.
  Rank agreement, reported only: Spearman rho of per-factor delta across the
  1,836 factors for A-C, B-C and A-B; overlap of the top 20 by delta for each
  pair; overlap of the toward_young sets.
  Not run in A or B:
  - the S2 PCA settings: their basis was fit on z-scored GTEx data, and
    refitting one in A or B would be a new fit, which this task excludes;
  - O2 and Y2 (secondary in SOUTH6);
  - proliferation-residualised cos and delta;
  - Stage 3 bootstrap CIs (never a gate in SOUTH6);
  - Stages 4 and 5;
  - the 61-guide report-only floor;
  - the frozen-ruler score change.
SCREEN OUTCOME PER QUALIFIED SPACE: UNCHANGED if its toward_young set equals
space C's (space C as stored has 0 factors in S1 O1_Y1); CHANGED otherwise.

--- STAGE 4: DIAGNOSTIC, TOP-10 GENE SHARE (report only, never a gate) ---
Euclidean distance adds up over genes only in its square, so the share is
taken on squared distance. For a difference vector x:
  share10(x) = (sum of the 10 largest x_j^2) / (sum over all j of x_j^2)
In each space, report share10 and the 10 gene symbols for:
  SAME: x(O) - x(Y) (start distance), x(S) - x(Y) (final distance) and d_S;
        for the reverse test, x(S_rev) - x(O) and d_rev.
  SOUTH6: O1 - Y1, with symbols; and, across factors, the median and the 5th
        and 95th percentiles of share10 for d_g and for (O1 + d_g) - Y1.
Wording rule: wherever a space's SAME result is stated, the share10 of its
final distance is stated in the same sentence.

--- EXPECTATIONS, stated before running (not gates) ---
Context established before this file and not re-derived here: in the screen,
effect lengths reflect sequencing depth, not biology. ||d|| against n_cells
has Spearman ~ -0.99, and ||d|| ~ 1342 / sqrt(n_cells), which matches pure
sampling noise.
Expected: the two halves behave differently.
  - Screen: A weights abundant genes, which are the ones measured most
    reliably at low cell counts, so A may help the screen. Its floors may fall
    relative to ||O1 - Y1||, and its ranking may differ from C's more than
    B's does.
  - SAME: cell counts are in the thousands, and depth is not the problem. A
    lets a few highly expressed genes dominate the whole distance, which is a
    risk here. share10 is there to show whether A's SAME result rests on a
    handful of genes.
  - B: no directional expectation.
These expectations change no threshold and no reading. If the halves differ as
expected, that is not evidence for or against either half.

--- CROSS-SPACE RULES ---
Raw distances, deltas and floors are in different units in A, B and C, and
they are never compared across spaces. Cross-space statements use only: the
sign of delta_S, the SAME key, frac_S, R_S, the toward_young sets and counts,
rank correlations and top-20 overlap, and floor / ||O1 - Y1||. Spaces are
never pooled or averaged.

--- DECISION TABLE (written before any number exists) ---
It governs the preprint's result, which is the SAME half.

| Outcome                                   | Response                                               |
|-------------------------------------------|--------------------------------------------------------|
| Sign holds in all surviving spaces        | Strongest version of the paper; one supplementary note |
| A space flips but fails its mixture control | Report as disqualified; no change to conclusions     |
| A space flips AND passes its controls     | Revise the preprint (not retract); report immediately  |

Definitions:
  - surviving space = QUALIFIED (P passes in it), C included.
  - flips = delta_S < 0 in that space (Stage 2).
  - passes its controls = QUALIFIED. The reversed-roles control does not
    qualify a space. It is reported in every space and decides the SAME key
    inside it.
  - Rows 1 and 2 can both apply (for example, A disqualified and flipped, B
    qualified and holding); report both. Row 3 overrides rows 1 and 2.
  - If neither A nor B qualifies, row 1 holds only trivially (C alone). The
    supplementary note must then say that no alternative space could detect a
    known approach, so coordinate robustness was not tested.

Pre-declared wording:
  Row 1: "In every coordinate space that detected known young-cell mixtures
    ([spaces]), aged partially reprogrammed cells end farther from the young
    donor than they started ([R_S and share10 per space]). The same-platform
    result does not depend on the log or the z-score." Then list any
    disqualified space.
  Row 2: "Space [X] ends with delta_S < 0 but failed its mixture control. Its
    numbers are reported as disqualified and support no conclusion."
  Row 3, by the key that fired in the flipping space:
    `toward_young_same_platform`: "In space [X], which detected known
      young-cell mixtures, aged partially reprogrammed cells end closer to the
      young donor than they started and approach it more than the reverse.
      The same-platform result depends on the choice of coordinates."
    `convergence_not_rejuvenation`: "In space [X], which detected known
      young-cell mixtures, aged partially reprogrammed cells end closer to the
      young donor than they started, but the reverse test does not separate
      this from mutual convergence of the two donors. The statement that the
      cells end farther from young depends on the choice of coordinates."
    `no_approach_same_platform` or `mixed`, with delta_S < 0: "In space [X],
      which detected known young-cell mixtures, the distance to the young
      donor falls (delta_S < 0), but [the condition that failed]. The sign of
      the change in distance depends on the choice of coordinates."
  If row 3 fires:
  - write results/coordspace/REVISION_TRIGGERED.flag naming the space and the
    key, and put it at the top of FINDINGS_COORDSPACE.md;
  - stop once the findings are written, and report to the investigator before
    anything else is run.
  Space C stays primary. The revision adds the space dependence; it does not
  replace the C result.

SOUTH6 half: the preprint names the screen only as a direct next test and
does not report it. Each qualified space's screen outcome (UNCHANGED or
CHANGED) is reported against SOUTH6's own recorded reading. A CHANGED space is
reported immediately with its factor list, and SOUTH6 space C remains the
screen's primary. The screen half does not by itself trigger any row of the
decision table.

--- WHAT DOES NOT COUNT ---
- Promoting A or B to primary, or choosing a space after seeing results.
- Arguing which coordinates are right. Qualification is P and nothing else.
- Changing f, the P pass rule, the split, seeds, n_perm, n_boot, the gates,
  the factors, the control set, or the origin and target after seeing numbers.
- Refitting the ruler, or fitting any PCA or other basis in A or B.
- Applying a space C floor in A or B, or comparing raw magnitudes across
  spaces.
- Clipping, flooring or re-centring coordinates, or dropping the genes that
  dominate share10.
- Adding a space (square root, arcsinh, z without log, ...) after seeing these
  results. Any further space needs its own dated pre-registration.
- Pooling or averaging spaces, or letting one half override the other.
- Reading anything from a DISQUALIFIED block.

--- OUT OF SCOPE ---
TOWARD's GTEx-anchored test; GENESPACE's MD and age gene spaces; SOUTH6's S2
PCA settings, O2 and Y2, and Stages 4 and 5; every other analysis. GSE325735
is not opened.

--- OUTPUT ---
FINDINGS_COORDSPACE.md contains:
- this flag, verbatim;
- the Stage 0 reproduction table;
- the per-space P table, with the noise floor and qualification;
- SAME per space: forward, reverse, asymmetry, context, key, sign, R_S and
  share10;
- SOUTH6 per space: floors, MDA, counts, the toward_young list and rank
  agreement;
- the diagnostics;
- the decision-table row or rows that fired, in the pre-declared wording;
- DISQUALIFIED blocks, kept separate;
- limitations.
results/coordspace/ contains:
- PREREG.flag (this file);
- QUALIFIED_BEFORE_PRIMARY.flag;
- REVISION_TRIGGERED.flag, if row 3 fires;
- per-stage CSVs, summary.json, manifest.json and the run log.
PROGRESS_COORDSPACE.md gives the stop status, the terminal command and the
runtime.

--- LIMITATIONS, declared now ---
- SAME rests on two donors. No choice of coordinates changes that.
- B and C share the prior count and A has none, so A and B differ in the log
  and the prior together.
- TMM factors are shared across spaces and are themselves estimated from
  log-ratios. This task does not test TMM.
- In A, the screen's displacements are Hs27 CPM differences added to a GTEx CPM
  origin, and a gene's absolute scale in Hs27 and in GTEx need not match.
  Negative coordinates are kept and counted.
- identity_loss and guides_disagree are not recomputed in A or B.
- Passing P shows that a space can see a known approach between these two
  donors' day-0 cells. It does not show that the space can see one across the
  bulk-to-single-cell gap, which is where the screen's GTEx origin and target
  sit. Stage 2C is the screen's own check on that.
```

## Notes on the pre-registration text

The flag is a timestamped record and was not edited. Two statements in it are factually wrong about SOUTH6. Neither changes any computation.

- **matched coordinates.** The flag says: 20,265 matched genes; permutation across the 20,265 matched coordinates. Actual: 20,264. SOUTH6's coverage table counts 20,265 ruler genes matched by best_after_md3_idtype; two of them map to one source column, which SOUTH6 collapses, so its panel, its d_g and its Stage 2A permutation run over 20,264 coordinates (results/south6/stage1.npz, matched_ruler). Effect: None on any computation. The pre-registered rule is 'the same matched genes, unchanged', and this task uses SOUTH6's own matched set.
- **identity_loss units.** The flag says: the 0.5 identity_loss threshold is defined in frozen-z units. Actual: In src/south6_run.py identity_drop is computed from the panel's log2-CPM (logcpm_f, logcpm_pool), so the 0.5 threshold is in log2-CPM units, the units of space B, not z units. Effect: None on any computation. The flags are carried unchanged from SOUTH6 space C as pre-registered. In space B a recomputed identity_drop would equal the carried one exactly; space A still has no counterpart.

## Stage 0 — reproduction of space C

Gate: every categorical outcome identical, and every compared number within |new − stored| ≤ tol × max(1, |stored|), tol = 1e-06 for SAME and the GTEx centroids, 0.0001 for SOUTH6. Checks marked diagnostic are not in the pre-registered list and do not gate.

**Result: PASS.** Gating checks passed: 72 of 72. Diagnostic checks failed: 0.

| block | item | kind | n | max_abs_diff | max_scaled_diff | tol | n_mismatch | pass | gating |
|---|---|---|---|---|---|---|---|---|---|
| SAME | stored split equals a fresh split_half draw | categorical | 1 | NA | NA | NA | 0 | True | False |
| SAME posctrl.csv | f | numeric | 4 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME posctrl.csv | n_aged | numeric | 4 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME posctrl.csv | n_young | numeric | 4 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME posctrl.csv | n_total | numeric | 4 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME posctrl.csv | frac | numeric | 4 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME posctrl.csv | delta | numeric | 4 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME posctrl.csv | cos | numeric | 4 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME posctrl.csv | d_norm | numeric | 4 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME posctrl.csv | v_norm | numeric | 4 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME posctrl.csv | frac_ci_lo | numeric | 4 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME posctrl.csv | frac_ci_hi | numeric | 4 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME posctrl.csv | delta_ci_lo | numeric | 4 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME posctrl.csv | delta_ci_hi | numeric | 4 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME posctrl.csv | ruler_score_context | numeric | 4 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME posctrl_summary.json | P pass | categorical | 1 | NA | NA | NA | 0 | True | True |
| SAME stats.csv | n_cells_S | numeric | 2 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME stats.csv | cos | numeric | 2 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME stats.csv | frac | numeric | 2 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME stats.csv | delta | numeric | 2 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME stats.csv | p_N1 | numeric | 2 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME stats.csv | frac_ci_lo | numeric | 2 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME stats.csv | frac_ci_hi | numeric | 2 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME stats.csv | delta_ci_lo | numeric | 2 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME stats.csv | delta_ci_hi | numeric | 2 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME stats.csv | d_norm | numeric | 2 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME stats.csv | v_norm | numeric | 2 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME stats.csv | ruler_score_S | numeric | 2 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME summary.json | asymmetry | numeric | 1 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME summary.json | asymmetry_ci_lo | numeric | 1 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME summary.json | asymmetry_ci_hi | numeric | 1 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME context.csv | n_cells | numeric | 4 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME context.csv | ruler_score_context | numeric | 4 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME context.csv | frac | numeric | 4 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME context.csv | delta | numeric | 4 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME context.csv | cos | numeric | 4 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME context.csv | frac_ci_lo | numeric | 4 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME context.csv | frac_ci_hi | numeric | 4 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME context.csv | delta_ci_lo | numeric | 4 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME context.csv | delta_ci_hi | numeric | 4 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME context.csv | ruler_score_context (O, Y) | numeric | 2 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SAME reading | key | categorical | 1 | NA | NA | NA | 0 | True | True |
| SAME reading | key is no_approach_same_platform | categorical | 1 | NA | NA | NA | 0 | True | True |
| SAME reading | delta_S (+130.894) | numeric | 1 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| GTEx centroids | n_young (113), n_old (231) | numeric | 2 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| GTEx centroids | B young vs mu + sd * C young | numeric | 23485 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| GTEx centroids | B old vs mu + sd * C old | numeric | 23485 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| GTEx centroids | C young recomputed vs anchors.npz (diagnostic) | numeric | 23485 | +0.000 | +0.000 | +0.000 | 0 | True | False |
| GTEx centroids | C old recomputed vs anchors.npz (diagnostic) | numeric | 23485 | +0.000 | +0.000 | +0.000 | 0 | True | False |
| SOUTH6 stage1.npz | panel missing-gene mask | categorical | 20264 | NA | NA | NA | 0 | True | False |
| SOUTH6 stage1.npz | d_g, all factors and genes (diagnostic) | numeric | 43118460 | +0.000 | +0.000 | +0.000 | 0 | True | False |
| SOUTH6 origin_target_distances.csv | ||O1 - Y1|| | numeric | 1 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SOUTH6 stage2a_floors.npz | per-factor floor | numeric | 1836 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SOUTH6 stage2a_floors.npz | per-factor permutation count le | numeric | 1836 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SOUTH6 stage2a_floors.npz | pooled floor | numeric | 1 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SOUTH6 stage2a_floors.npz | median ||d|| | numeric | 1 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SOUTH6 stage2c_mda.csv | mda_primary | numeric | 1 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SOUTH6 stage2c_mda.csv | delta_at_mda_primary | numeric | 1 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SOUTH6 stage2c_mda.csv | pooled_floor | numeric | 1 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SOUTH6 stage2c_mda.csv | mda_stricter | numeric | 1 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SOUTH6 stage2c_mda.csv | delta_at_mda_stricter | numeric | 1 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SOUTH6 stage2c_mda.csv | stricter_floor | numeric | 1 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SOUTH6 stage2c_mda.csv | mda_margin | numeric | 1 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SOUTH6 stage2c_mda.csv | median_norm | numeric | 1 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SOUTH6 stage2c_mda.csv | nt_guide_floor | numeric | 1 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SOUTH6 stage2c_mda.csv | split_floor | numeric | 1 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SOUTH6 stage2c_mda.csv | stricter_source | categorical | 1 | NA | NA | NA | 0 | True | True |
| SOUTH6 stage3_full_table.csv (S1 O1_Y1) | delta | numeric | 1836 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SOUTH6 stage3_full_table.csv (S1 O1_Y1) | cos | numeric | 1836 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SOUTH6 stage3_full_table.csv (S1 O1_Y1) | frac | numeric | 1836 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SOUTH6 stage3_full_table.csv (S1 O1_Y1) | d_norm | numeric | 1836 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SOUTH6 stage3_full_table.csv (S1 O1_Y1) | floor_p05 | numeric | 1836 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SOUTH6 stage3_full_table.csv (S1 O1_Y1) | p_delta | numeric | 1836 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SOUTH6 stage3_full_table.csv (S1 O1_Y1) | q_delta | numeric | 1836 | +0.000 | +0.000 | +0.000 | 0 | True | True |
| SOUTH6 stage3_full_table.csv (S1 O1_Y1) | below floor (per factor) | categorical | 1836 | NA | NA | NA | 0 | True | True |
| SOUTH6 stage3_full_table.csv (S1 O1_Y1) | q <= 0.05 (per factor) | categorical | 1836 | NA | NA | NA | 0 | True | True |
| SOUTH6 stage3_full_table.csv (S1 O1_Y1) | toward_young (per factor) | categorical | 1836 | NA | NA | NA | 0 | True | True |

## Stage 1 — mixture control in every space

SAME positive control P, unchanged, in each space. A space is QUALIFIED only if frac rises monotonically with f, the frac CI at f = 0.50 excludes 0 and contains its own point estimate, and delta at f = 0.50 < 0. The f = 0 row is the noise floor and is not subtracted.

| space | qualified | monotonic | frac50_ci_excludes_0 | frac50_ci_valid | delta50_lt_0 | pass_same_rule |
|---|---|---|---|---|---|---|
| A | True | True | True | True | True | True |
| B | True | True | True | True | True | True |
| C | True | True | True | True | True | True |

| space | f | n_total | frac | frac_ci_lo | frac_ci_hi | delta | delta_ci_lo | delta_ci_hi | cos | frac_ci_covers_point | delta_ci_covers_point | noise_floor |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | +0.000 | 2470 | +0.039 | -0.009 | +0.085 | -497.940 | -1059.048 | +138.209 | +0.738 | True | True | True |
| A | +0.100 | 2470 | +0.104 | +0.064 | +0.155 | -1338.440 | -1968.745 | -813.686 | +0.948 | True | True | False |
| A | +0.250 | 2470 | +0.218 | +0.167 | +0.263 | -2820.358 | -3342.202 | -2134.769 | +0.989 | True | True | False |
| A | +0.500 | 2470 | +0.417 | +0.371 | +0.459 | -5402.161 | -5833.482 | -4741.512 | +0.994 | True | True | False |
| B | +0.000 | 2470 | +0.106 | +0.091 | +0.108 | +0.342 | +5.078 | +6.805 | +0.227 | True | False | True |
| B | +0.100 | 2470 | +0.195 | +0.181 | +0.197 | -8.335 | -3.232 | -1.366 | +0.400 | True | False | False |
| B | +0.250 | 2470 | +0.305 | +0.290 | +0.310 | -18.303 | -12.505 | -10.611 | +0.564 | True | False | False |
| B | +0.500 | 2470 | +0.461 | +0.446 | +0.466 | -30.690 | -23.478 | -21.398 | +0.708 | True | False | False |
| C | +0.000 | 2470 | +0.115 | +0.099 | +0.121 | +0.068 | +5.807 | +8.334 | +0.239 | True | False | True |
| C | +0.100 | 2470 | +0.182 | +0.167 | +0.191 | -7.870 | -1.900 | +1.007 | +0.371 | True | False | False |
| C | +0.250 | 2470 | +0.283 | +0.267 | +0.291 | -18.913 | -11.793 | -9.372 | +0.530 | True | False | False |
| C | +0.500 | 2470 | +0.441 | +0.424 | +0.447 | -34.443 | -25.403 | -22.414 | +0.694 | True | False | False |

## Stage 2 — SAME in every space

Raw distances are in different units in each space and are not compared across spaces. Cross-space statements use the sign of delta_S, the key, frac_S and R_S = final / start distance.

### Qualified spaces

| space | key | sign | flags | ci_invalid | n_cells_S | delta_S | R_S | frac_S | frac_S_ci_lo | frac_S_ci_hi | p_N1 | frac_rev | asymmetry | asymmetry_ci_lo | asymmetry_ci_hi | pluri_delta | cos_v_u_toward |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | no_approach_same_platform | HOLDS | NA | NA | 5832 | +17020.885 | +2.054 | +2.045 | +2.026 | +2.052 | +0.005 | -1.208 | +3.253 | +3.224 | +3.261 | +26676.177 | +0.421 |
| B | no_approach_same_platform | HOLDS | NA | NA | 5832 | +99.407 | +1.879 | +0.581 | +0.574 | +0.583 | +0.005 | +0.074 | +0.507 | +0.498 | +0.511 | +196.392 | +0.008 |
| C | no_approach_same_platform | HOLDS | NA | NA | 5832 | +130.894 | +2.002 | +0.674 | +0.667 | +0.680 | +0.005 | -0.057 | +0.732 | +0.724 | +0.740 | +382.106 | +0.039 |

- Space A: GM00731 PartialReprog (n_cells = 5832) ends farther from the young donor than it started (delta_S = +1.702e+04 in space A units, R_S = 2.054); the top 10 genes carry 77.3% of the squared final distance (MALAT1, KLF4, VIM, LGALS1, POU5F1, S100A6, FTL, TMSB10, COL1A2, TMSB4X). Key `no_approach_same_platform` (delta_S >= 0).
- Space B: GM00731 PartialReprog (n_cells = 5832) ends farther from the young donor than it started (delta_S = +99.41 in space B units, R_S = 1.879); the top 10 genes carry 3.6% of the squared final distance (SOX2, POU5F1, HBA2, HBA1, LY6D, SPRR3, KRT6A, KLF4, PI3, POU5F1B). Key `no_approach_same_platform` (delta_S >= 0).
- Space C: GM00731 PartialReprog (n_cells = 5832) ends farther from the young donor than it started (delta_S = +130.9 in space C units, R_S = 2.002); the top 10 genes carry 1.8% of the squared final distance (POU5F1, SOX2, KLF4, LY6D, SLC9A3, CALY, OXR1, GPX2, ELP5, CALB1). Key `no_approach_same_platform` (delta_S >= 0).

### DISQUALIFIED spaces (no reading is drawn from these numbers)

None.

Context, reported only (frozen-ruler score in space C only):

| space | pseudobulk | n_cells | ruler_score_context | frac | delta | cos | frac_ci_lo | frac_ci_hi | delta_ci_lo | delta_ci_hi |
|---|---|---|---|---|---|---|---|---|---|---|
| A | S | 5832 | NA | +2.045 | +17020.885 | +0.756 | +2.026 | +2.052 | +16273.713 | +17617.548 |
| A | S_rev | 6481 | NA | -1.208 | +26863.619 | -0.629 | -1.213 | -1.195 | +26648.532 | +27040.084 |
| A | Pluri | 439 | NA | +2.570 | +26676.177 | +0.769 | +2.530 | +2.598 | +25983.482 | +27150.064 |
| A | NonReprog | 5459 | NA | +1.199 | +3085.300 | +0.715 | +1.182 | +1.222 | +2632.228 | +3253.851 |
| B | S | 5832 | NA | +0.581 | +99.407 | +0.303 | +0.574 | +0.583 | +99.273 | +101.579 |
| B | S_rev | 6481 | NA | +0.074 | +114.095 | +0.042 | +0.070 | +0.078 | +113.945 | +116.202 |
| B | Pluri | 439 | NA | +0.773 | +196.392 | +0.273 | +0.741 | +0.778 | +206.429 | +217.632 |
| B | NonReprog | 5459 | NA | +0.490 | +38.476 | +0.368 | +0.479 | +0.493 | +40.663 | +42.109 |
| C | S | 5832 | +0.929 | +0.674 | +130.894 | +0.323 | +0.667 | +0.680 | +130.473 | +133.918 |
| C | S_rev | 6481 | -4.350 | -0.057 | +163.342 | -0.029 | -0.063 | -0.054 | +163.312 | +166.330 |
| C | Pluri | 439 | -10.041 | +1.403 | +382.106 | +0.338 | +1.365 | +1.442 | +391.791 | +418.643 |
| C | NonReprog | 5459 | +7.089 | +0.555 | +60.331 | +0.370 | +0.543 | +0.561 | +62.804 | +64.575 |

## Stage 3 — SOUTH6 S1 O1→Y1 in every space

identity_loss and guides_disagree are carried from SOUTH6 space C. Floors are recomputed in each space. Raw floors and deltas are in each space's units; the unit-free comparison is median floor / ||O1 − Y1||.

### Qualified spaces

| space | qualified | outcome | gap | median_d_norm | floor_min | floor_median | floor_max | pooled_floor | median_floor_over_gap | nt_guide_floor | split_floor | mda_primary | mda_stricter | n_delta_negative | n_below_floor | n_q | n_toward_young | n_tie_at_min_p | n_negative_coords_A | n_factors_with_negative_coords_A |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | True | CHANGED | +8225.098 | +3563.583 | +106.023 | +719.249 | +10030.506 | +719.939 | +0.087 | +23.544 | -14.178 | +0.050 | +0.050 | 2 | 1292 | 1272 | 2 | 44 | +3497133.000 | +1836.000 |
| B | True | UNCHANGED | +29.906 | +51.655 | +11.737 | +29.653 | +115.299 | +29.642 | +0.992 | +11.593 | +0.858 | +0.050 | +0.050 | 0 | 1088 | 1063 | 0 | 838 | NA | NA |
| C | True | reference | +32.284 | +123.350 | +42.416 | +94.964 | +363.574 | +94.966 | +2.942 | +41.846 | +5.483 | +0.050 | +0.050 | 0 | 377 | 202 | 0 | 63 | NA | NA |

toward_young factors, every space:

| space | factor | delta | floor_p05 | q_delta | n_cells |
|---|---|---|---|---|---|
| A | ZNF441 | -4.525 | +374.613 | +0.005 | 84 |
| A | ZNF669 | -69.489 | +302.000 | +0.005 | 82 |

### DISQUALIFIED spaces (no conclusion is drawn from these numbers)

None.

Rank agreement, reported only (all 1,836 factors, S1 O1→Y1):

| pair | spearman_delta | top20_overlap | toward_young_overlap | toward_young_a | toward_young_b |
|---|---|---|---|---|---|
| A-C | +0.241 | 1 | 0 | 2 | 0 |
| B-C | +0.990 | 14 | 0 | 0 | 0 |
| A-B | +0.283 | 1 | 0 | 2 | 0 |

## Stage 4 — diagnostic: share of squared distance in the top 10 genes

share10 = (sum of the 10 largest squared coordinates) / (squared length). Report only, never a gate.

| space | half | vector | share10 | top10_genes | p05 | p95 |
|---|---|---|---|---|---|---|
| A | SAME | start distance x(O) - x(Y) | +0.727 | FN1, FTL, MALAT1, IGFBP3, FTH1, TMSB4X, RPL41, COL1A2, RPS19, NEAT1 | NA | NA |
| A | SAME | final distance x(S) - x(Y) | +0.773 | MALAT1, KLF4, VIM, LGALS1, POU5F1, S100A6, FTL, TMSB10, COL1A2, TMSB4X | NA | NA |
| A | SAME | displacement d_S | +0.720 | MALAT1, VIM, KLF4, FTL, FN1, LGALS1, IGFBP3, POU5F1, FTH1, S100A6 | NA | NA |
| A | SAME | reverse final distance x(S_rev) - x(O) | +0.715 | MALAT1, VIM, FTL, FN1, FTH1, LGALS1, IGFBP3, S100A6, TMSB10, TMSB4X | NA | NA |
| A | SAME | reverse displacement d_rev | +0.720 | MALAT1, VIM, LGALS1, FTL, S100A6, TMSB10, FTH1, COL1A2, TMSB4X, MT-CO3 | NA | NA |
| B | SAME | start distance x(O) - x(Y) | +0.034 | F2RL2, MEOX2, PI16, RSPO4, RASIP1, PSG4, PIEZO2, EFHD1, CMKLR1, MIR924HG | NA | NA |
| B | SAME | final distance x(S) - x(Y) | +0.036 | SOX2, POU5F1, HBA2, HBA1, LY6D, SPRR3, KRT6A, KLF4, PI3, POU5F1B | NA | NA |
| B | SAME | displacement d_S | +0.033 | SOX2, HBA2, POU5F1, LY6D, HBA1, SPRR3, KRT6A, SPINK1, KLK10, PI3 | NA | NA |
| B | SAME | reverse final distance x(S_rev) - x(O) | +0.025 | SOX2, HBA2, LY6D, HBA1, POU5F1, KRT6A, CALB1, SPRR3, IGFL2-AS1, MT1H | NA | NA |
| B | SAME | reverse displacement d_rev | +0.033 | HBA2, SOX2, HBA1, POU5F1, LY6D, KRT6A, SPRR3, MAL, CALB1, KLK1 | NA | NA |
| C | SAME | start distance x(O) - x(Y) | +0.022 | ELP5, ABCC9, KCNJ8, ENSG00000268081, HOXC10, TDRP, MTRNR2L8, PLCE1, MGST1, DNAH17 | NA | NA |
| C | SAME | final distance x(S) - x(Y) | +0.018 | POU5F1, SOX2, KLF4, LY6D, SLC9A3, CALY, OXR1, GPX2, ELP5, CALB1 | NA | NA |
| C | SAME | displacement d_S | +0.016 | POU5F1, SOX2, SLC9A3, KLF4, LY6D, CALB1, NRARP, PARP12, MYC, PITX1-AS1 | NA | NA |
| C | SAME | reverse final distance x(S_rev) - x(O) | +0.013 | POU5F1, MIR497HG, SOX2, PPP1R13B, MT1H, CALB1, MT1G, SLC9A3, LY6D, HS3ST5 | NA | NA |
| C | SAME | reverse displacement d_rev | +0.015 | POU5F1, LY6D, SOX2, CALB1, CALY, SLC9A3, MT1H, KLF4, TPPP, SEMA4D | NA | NA |
| A | SOUTH6 | O1 - Y1 (gap) | +0.937 | FN1, SERPINE1, COL1A2, TGFBI, FTL, MT-ND2, COL1A1, MT2A, COL3A1, ANXA2 | NA | NA |
| A | SOUTH6 | d_g across factors | +0.677 | NA | +0.401 | +0.903 |
| A | SOUTH6 | (O1 + d_g) - Y1 across factors | +0.896 | NA | +0.852 | +0.921 |
| B | SOUTH6 | O1 - Y1 (gap) | +0.021 | CHI3L1, CADM3-AS1, CADM3, DUXAP8, ACKR1, PRSS2, ENSG00000228560, DUXAP10, PRSS1, KIF20A | NA | NA |
| B | SOUTH6 | d_g across factors | +0.016 | NA | +0.014 | +0.026 |
| B | SOUTH6 | (O1 + d_g) - Y1 across factors | +0.014 | NA | +0.012 | +0.020 |
| C | SOUTH6 | O1 - Y1 (gap) | +0.006 | DUXAP8, ENSG00000177788, DUXAP10, ENSG00000255310, RFPL1S, XKR6, ZNF518B, ENSG00000228560, DUXAP9, MARVELD2 | NA | NA |
| C | SOUTH6 | d_g across factors | +0.038 | NA | +0.030 | +0.048 |
| C | SOUTH6 | (O1 + d_g) - Y1 across factors | +0.036 | NA | +0.028 | +0.046 |

## Limitations

- SAME rests on two donors. No choice of coordinates changes that.
- B and C share the prior count and A has none, so A and B differ in the log and the prior together.
- TMM factors are shared across spaces and are themselves estimated from log-ratios. This task does not test TMM.
- In A, the screen's displacements are Hs27 CPM differences added to a GTEx CPM origin; negative coordinates are kept and counted above.
- identity_loss and guides_disagree are not recomputed in A or B.
- Passing P shows that a space can see a known approach between these two donors' day-0 cells, not across the bulk-to-single-cell gap where the screen's GTEx origin and target sit. Stage 2C is the screen's own check on that.
- SOUTH6 Stage 1 count sums are read from results/south6/stage1.npz and results/south6/matched_* as SOUTH6 wrote them; Stage 0 shows they reproduce SOUTH6's space C vectors.

## STOP / failures

None recorded.

## Files

- `src/coordspace_run.py`
- `results/coordspace/PREREG.flag` (not rewritten)
- `results/coordspace/coordspace_run.log`
- `results/coordspace/decision.json`
- `results/coordspace/full_stderr.txt`
- `results/coordspace/full_stdout.txt`
- `results/coordspace/manifest.json`
- `results/coordspace/QUALIFIED_BEFORE_PRIMARY.flag`
- `results/coordspace/stage0_reproduction.csv`
- `results/coordspace/stage0_stderr.txt`
- `results/coordspace/stage0_stdout.txt`
- `results/coordspace/stage0_summary.json`
- `results/coordspace/stage1_posctrl.csv`
- `results/coordspace/stage1_qualification.csv`
- `results/coordspace/stage2_same_context.csv`
- `results/coordspace/stage2_same_stats.csv`
- `results/coordspace/stage2_same_summary.csv`
- `results/coordspace/stage3_rank_agreement.csv`
- `results/coordspace/stage3_south6_factors.csv`
- `results/coordspace/stage3_south6_summary.csv`
- `results/coordspace/stage3_toward_young.csv`
- `results/coordspace/stage4_share10.csv`
- `results/coordspace/summary.json`
- `FINDINGS_COORDSPACE.md`, `PROGRESS_COORDSPACE.md`


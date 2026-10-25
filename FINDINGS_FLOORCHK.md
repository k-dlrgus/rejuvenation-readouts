# FINDINGS_FLOORCHK — do paper 1's two length claims survive a size-matched noise floor?

**Task:** FLOORCHK. **Status:** complete. No STOP fired. Reproduction: 30 gating checks, 0 mismatches.
Every gate passed: claim 1 in S1, S2_k20, S2_k50 and S2_k100; claim 2 in all, MD and age genes.

**Primary readings, in pre-registered wording:**

```
CONCLUSION  CLAIM1-PARTIAL: 18 of 1,836
CONCLUSION  CLAIM2-HOLDS
```

**Claim 1 does not hold as written.** In S1 (paper 1's "space C"), random non-targeting cells with no
effect already read more than twice the GTEx old-to-young distance at every group size in the screen.
Their median is 2.06 at 341 cells, 3.76 at the median 106 cells and 11.9 at 13 cells. The published
"every one of the 1,836 perturbations moved more than twice the distance (smallest ratio 2.10)" is
mostly sampling noise. After the pre-registered noise correction, 18 of 1,836 are estimated to move
more than twice the distance, against about 5 expected from noise alone (a diagnostic, not
pre-registered). The median estimate is 0.61. The 18 skew toward small groups (median 54.5 cells,
against 106 for the screen) and are not a shortlist. That no approach *could be registered* from that
origin is still true of the measurement. But it reflects the screen's depth, not perturbations
overshooting. The sentence that fails, and a corrected wording, are in "What no longer holds in paper 1"
below.

**Claim 2 holds, narrowly in all genes.** The aged donor's 5,832 partially reprogrammed cells give a
corrected ratio of 2.046 with the floor extrapolated to 5,832 cells, and 2.022 with the observed
2,470-cell floor, against a raw 2.087. As secondary readings, the MD genes give 2.296 (raw 2.300; the
main-text "2.3") and the age genes 2.047 (raw 2.084). Both used the observed floor alone, because the
1/n fit failed its 10% check there.

The parts that undercut claim 2:
- The 5,832-cell floor is an extrapolation about 2.4 times past the largest size drawn.
- The floor comes from day-0 fibroblasts, not partially reprogrammed cells. A 33% excess of real noise
  over the observed floor would put the all-genes reading below 2.
- The uncorrected gap biases the corrected ratio downward, which works against the claim, not for it.

**Secondary readings.** S2_k20/k50/k100 fired FLOOR-REACHES-2 at 1, 3 and 4 of 212 sizes, all small,
and LENGTH-DEPTH for the 1, 2 and 2 factors whose raw ratio was ≥ 2. No factor reaches a corrected
ratio of 2 there; the maximum is 1.67. Claim 2 MD and age fired LENGTH-STANDS. All readings are listed
below and in `results/floorchk/readings.txt`.

Files are all in `results/floorchk/`:
- `PREREG.flag`, `reproduction.csv`
- Floors: `claim1_floor.csv`, `claim2_floor.csv`, `claim2_fit.csv`
- Gates: `claim1_gate.csv`, `claim1_gate_summary.csv`, `claim2_gate.csv`, `claim2_gate_summary.csv`
- Units: `claim1_units.csv`, `claim1_summary.csv`, `claim2_unit.csv`
- Gap and diagnostic: `claim2_gap_spread.csv`, `claim1_diagnostic.csv` (not pre-registered)
- Readings and logs: `readings.txt`, `summary.json`, `floorchk_run.log`, `run_stdout.txt`, `run_stderr.txt`
- Draws: `claim1_floor_draws.npz` (git-ignored, regenerable), `claim2_floor_draws.npz`

Code: `src/floorchk_run.py`, with modes `run`, `verify` and `diagnostic`. No preprint text was edited.

## Pre-registration (verbatim copy of `results/floorchk/PREREG.flag`)

~~~~text
FLOORCHK pre-registration flag
Written 2026-10-25T23:17:26+09:00 (2026-10-25T14:17:26Z), before any statistic of this task
exists. No floor draw, gate draw, corrected ratio, re-split gap, fit or new pseudobulk has been
computed. No new code has been run on any expression data.

What had been read or seen before this file was written (disclosed):
 - Repository files: CLAUDE.md; FINDINGS_SOUTH6.md (lines 1-552 and 743-771); FINDINGS_RATIO.md
   (lines 1-53 and 344-600); results/ratio/PREREG.flag (whole); FINDINGS_SAME.md (whole);
   FINDINGS_GENESPACE.md (lines 1-300, which include the stored per-space same-platform table and
   the per-space f = 0 mixture rows); code in src/south6_run.py, src/same_run.py,
   src/seng_run.py (tmm_cache_from_panel, z_from_counts), src/toward_run.py (tmm_logcpm_quiet,
   zscore_frozen, sum_rows, load_donor), src/genespace_run.py (function list); .gitignore.
 - Paper 1 (paper/rejuvenation_readouts_preprint_v1.2.docx, text extracted to a scratch file, not
   edited): the passages of the main text (paragraph on the 205 MD genes), the Discussion, the
   Limitations bullet on the same-platform test, and Supplementary Note 3.
 - Stored values read from files (these ARE the stored statistics the task asks to confirm; read,
   not recomputed by new code):
     results/south6/stage3_full_table.csv, O1_Y1, column d_over_v:
       S1 1,836 rows, min 2.0966, 1,836 with ratio >= 2, median 3.8208;
       S2_k20 min 0.1312, 1 >= 2; S2_k50 min 0.1940, 2 >= 2; S2_k100 min 0.2384, 2 >= 2.
     results/genespace/same_stats.csv, forward rows, d_norm / v_norm:
       FULL 272.7340 / 130.6884 = 2.0869; MD 33.5800 / 14.5976 = 2.3004;
       AGE 100.7160 / 48.3385 = 2.0836.
     results/south6/stage0_cells_per_target.csv (targets other than 'non'): 1,836 targets, cells
       per target 13..341, 212 distinct counts; numpy-linear percentiles P5/P25/P50/P95 =
       59 / 84 / 106 / 189.
     results/same/split.json: aged half A (O) 2,469, aged half B 2,470, young half A (Y) 3,851,
       young half B 3,852, aged d0 4,939, young d0 7,703. GM00731 PartialReprog 5,832 cells.
     results/south6/stage1.npz key list and shapes; results/same/panel_z.npz key list.
 - Outcome-relevant arithmetic done by hand before this file, disclosed because it bears on
   claim 2: the stored GENESPACE f = 0 mixture row (2,470 aged half-B cells vs O) has progress
   0.115 and perp 0.466 in FULL, i.e. a noise length of about 0.48 x ||gap||. Squared and
   subtracted from the stored ||d||, that puts the all-genes corrected ratio within a few percent
   of 2 (roughly 2.00 to 2.05 depending on the floor size used). For MD the same arithmetic gives
   about 2.30; for AGE about 2.05. This was stated to the user in chat before the rulings below.
   It is an estimate from stored numbers, not a result of this task.
 - Three questions were put to the user in chat; their rulings are reproduced verbatim below
   (U1-U3).

== Task (verbatim, as given in chat) ==

Read CLAUDE.md first; its rules apply. Then read FINDINGS_SOUTH6.md, FINDINGS_RATIO.md, results/ratio/PREREG.flag, and the paper-1 preprint passages cited below.

Task: `FLOORCHK`. Deliverables: `results/floorchk/PREREG.flag`, `FINDINGS_FLOORCHK.md`, `PROGRESS_FLOORCHK.md`. One local commit, no push. Modify no existing FINDINGS_*.md, no SCOPING_*.md, no file under results/south6/ or results/ratio/, and no preprint source.

## Why this task exists

RATIO (GSE278572) found that the displacement-length ratio ‖d‖/‖gap‖ carries a large sampling-noise component at low cell counts: random control-guide cells with no effect read 1.5–1.8 at 30 cells, and a planted full-gap displacement read 2.4–2.6 at 30 cells. SOUTH6 computed the same kind of ratio and never compared it with a size-matched control floor.

Two published length claims rest on such ratios. This task checks both against a noise floor. It changes no published text; it records whether the claims survive.

Claim 1 (paper 1, Supplementary Note 3; SOUTH6): in space C, for the GTEx old-to-young comparison with origin O1 and target Y1, every one of the 1,836 Southard perturbations moved more than twice the distance between the two centroids (smallest ratio 2.10), "so from that origin no approach to that target could have been registered whatever its direction". Locate the stored per-perturbation values in results/south6/ and confirm the 1,836 and the 2.10 from the files before anything else; if they do not reproduce, STOP and report.

Claim 2 (paper 1, main text and Limitations; the same-platform test on GSE297234): the aged donor's partially reprogrammed cells travel 2.3 times the line between the two donors' day-0 cells, "a move longer than twice that line cannot end closer to the target whichever way it points, so these cells overshoot on length alone". Locate the stored values for MD genes, age genes and all genes, and confirm them from the files before anything else; if they do not reproduce, STOP and report.

## Inputs — no download

Southard data and code are local: data/raw/southard/, the SOUTH6 run script and its stored intermediates under results/south6/. GSE297234 is local under data/raw/gse297234/. Reuse the existing pipelines' functions for normalisation, gene matching and coordinates exactly as they ran — no refit, no new basis, no new gene filter. If an intermediate needed to rebuild a pseudobulk is missing, rebuild it with the original code and show that the stored statistic reproduces to within rounding before using it. Download nothing.

## Pre-register before any statistic

Write the flag before computing anything new, including the statistic, the controls, the gates and the exact wording of every conclusion.

### Floor
For each perturbation (claim 1) or state (claim 2), draw 200 random samples of control cells at that unit's own cell count, with no effect, using the seeds in CLAUDE.md. Claim 1 controls: Southard non-targeting cells, as SOUTH6 defined them. Claim 2 controls: the aged donor's day-0 fibroblasts, drawn from the half NOT used to build the origin, matching the existing half-split. Compute ‖d0‖ for each draw exactly as ‖d‖ was computed. Report per space: median and P95 of ‖d0‖/‖gap‖ at each group size, beside the observed ratio.

For claim 2 the gap itself is a difference of two single-cell pseudobulks, so it also carries noise: estimate it by repeating the half-split for both day-0 arms (200 draws) and report the spread of ‖gap‖. For claim 1 the gap is between two GTEx bulk centroids; state that it is treated as fixed and why.

### Corrected ratio
Sampling noise adds approximately in quadrature to an independent effect. Define the noise-corrected ratio as sqrt(max(0, ‖d‖² − median(‖d0‖²))) / ‖gap‖. This is the pre-registered decision statistic. RATIO computed a similar estimate without pre-registration; this task pre-registers it.

### Gate — run before reading any real unit
Validate the correction with planted effects. Take control samples at four group sizes spanning each claim's range (claim 1: the 5th, 25th, 50th and 95th percentile of Southard cells per perturbation; claim 2: the partially reprogrammed population's own cell count), add a displacement of known length 0.5, 1, 2 and 3 times ‖gap‖ in a random direction (200 draws each), and compute the corrected ratio.
The gate passes in a space if the median corrected ratio is within ±20% of the planted length at every planted length >= 1, at every group size tested. If the gate fails in a space, the corrected ratio is NOT used there; report only the raw ratio beside the floor, and the conclusions below that require the corrected ratio do not fire in that space.

Spaces: claim 1, every space SOUTH6 reported for O1→Y1 (S1 and the compressed S2_k20/50/100), with S1 primary because it is the space paper 1 reports. Claim 2, the three gene spaces of paper 1 (MD 205, age 3,089, all 23,485), all-genes primary.

### Conclusions, worded now
Per claim, per space, every reading whose conditions fired:

  "FLOOR-REACHES-2: in <claim>, <space>, the control floor's P95 is >= 2 at <n> of <N> group sizes; at those depths a ratio of 2 is reached by sampling noise alone."
  "LENGTH-STANDS: in <claim>, <space>, <n> of <N> units have corrected ratio >= 2; the published 'more than twice' holds for them after removing the noise contribution."
  "LENGTH-DEPTH: in <claim>, <space>, <n> of <N> units have raw ratio >= 2 but corrected ratio < 2; for them the published 'more than twice' was reached partly by sampling noise."
  "LENGTH-UNRESOLVED: in <claim>, <space>, the correction failed its gate; only raw ratio and floor are reported and no statement about the published claim is made."

For claim 1 additionally, worded now, fires only in S1 if the gate passes there:
  "CLAIM1-HOLDS" if LENGTH-STANDS counts all 1,836 units; "CLAIM1-PARTIAL: <n> of 1,836" if it counts some; "CLAIM1-FAILS" if it counts none.
For claim 2 additionally, fires only in all-genes if the gate passes there:
  "CLAIM2-HOLDS" if the partially reprogrammed corrected ratio is >= 2; "CLAIM2-FAILS" if it is < 2. Report the MD and age spaces beside it, labelled secondary.

If CLAIM1-PARTIAL, CLAIM1-FAILS or CLAIM2-FAILS fires, write in FINDINGS_FLOORCHK.md, in plain words, which sentence of paper 1 no longer holds as written and what the corrected statement would be. Do not edit the preprint.

## Hard rules
No conclusion about any named factor or gene; this is about the published length claims only.
Reproduce every stored statistic used before computing anything new, and report the reproduction table.
Disclose every deviation, including anything decided that the flag did not specify. If a bar or definition is ambiguous, STOP and ask.

== Rulings by the user in chat, before this file was written ==

Questions are paraphrased; the text after "Ruling:" is the user's answer, verbatim.

[U1] Claim 2: the unit (GM00731 PartialReprog) has 5,832 cells but the aged half-B control pool
has 2,470, so neither the floor nor the gate can be drawn at the unit's own size without
replacement. How should they be drawn?
Ruling: Ladder + 1/n, as in the recommended option: draw without replacement from half B at 154, 309, 617, 1,235 and 2,470 cells (200 draws each), fit median ||d0||^2 = a + b/n, read it at 5,832, run the gate at each ladder size, and disclose the model.

Additions, fixed now, because 5,832 is beyond the largest ladder point:
- Report the fitted a and b, the fit's residuals at all five ladder points, and the floor read at 5,832 beside the floor observed at 2,470. State in the file that the 5,832 value is an extrapolation of roughly 2.4x past the largest size drawn.
- Compute the corrected ratio twice: once with the extrapolated floor at 5,832, and once with the observed floor at 2,470 used unchanged, which overstates the noise and so biases the corrected ratio against the published claim. Report both.
- CLAIM2-HOLDS fires only if the corrected ratio is >= 2 under BOTH. If the two disagree, fire LENGTH-UNRESOLVED for claim 2 instead and say that the verdict depends on the extrapolation.
- If the 1/n fit does not describe the ladder (any ladder point's median off the fitted curve by more than 10%), do not extrapolate: use the observed 2,470 floor alone and say so.

[U2] CLAUDE.md names no seeds. Which seeds?
Ruling: Fresh 20261027-29 (Recommended) — the option as offered: 20261027 for floors, 20261028
for the gate, 20261029 for the claim-2 gap re-splits.

[U3] Claim 1: what is each control draw's ||d0|| measured against?
Ruling: Pool minus the draw (Recommended) — the option as offered: the reference is the 14,013
label-defined non-targeting cells minus the drawn cells, as RATIO did, so the draw is independent
of its reference as a real perturbation is.

== Definitions fixed by the session now, because the task and rulings do not specify them
   (disclosed). Stated to the user in chat before the rulings, without objection: the shared
   per-size draws (D2), the single claim-2 unit (D1), the uncorrected claim-2 gap (D6), the S2
   rebuild (D9). Everything else below was NOT shown to the user before this file. ==

[D1] Units, spaces, gaps.
  Claim 1. Units: the 1,836 ranked factors of SOUTH6 (results/south6/stage1.npz). Unit cell count
    n_i = the factor's number of cells (stage1.npz cell_indptr; must equal stage0_cells_per_target).
    Spaces: S1 (full ruler coordinates, 23,485, unmatched ruler genes 0) and S2_k20/k50/k100
    (d projected on the first k GTEx-donor PCA loadings, as SOUTH6 project_D). Pair O1->Y1 only.
    ||d_i|| = SOUTH6's d_f row (S1) or its projection (S2). ||gap|| = ||c_old - c_young|| from
    results/toward/anchors.npz (S1), or ||(c_old - c_young) L_k|| (S2), i.e. SOUTH6's v_norm.
    The gap is treated as FIXED. Reason: c_old and c_young are means over 231 and 113 GTEx donor
    bulk profiles, frozen anchors read and not refit; they are not single-cell pseudobulks and
    carry no cell-sampling noise of the kind this task measures. Their donor-sampling uncertainty
    is a different quantity and is not estimated here.
  Claim 2. One unit: GM00731 PartialReprog (5,832 cells), the population claim 2 is about. No
    other state is a unit (Pluripotency, NonReprog and the reverse test are not tested).
    Spaces: FULL (all 23,485 ruler columns), MD (205 columns, results/genespace/gene_map.csv
    in_MD), AGE (3,089 columns, in_AGE), each sliced from the full-gene z, exactly as GENESPACE.
    d = z(S) - z(O), gap = z(Y) - z(O), from the SAME main panel [O, Y, S, S_rev, Pluri,
    NonReprog], TMM across that panel, frozen GTEx mu/sd, panel-missing genes z = 0.
    <space> in the reading wording is filled as "all genes (23,485)", "MD genes (205)",
    "age genes (3,089)" for claim 2 and "S1", "S2_k20", "S2_k50", "S2_k100" for claim 1;
    <claim> as "claim 1" / "claim 2".

[D2] Claim-1 floor (implements the task's Floor and U3).
  Per-cell matched counts: results/south6/matched_{data,indices,indptr} (SOUTH6's own per-cell
  matrix), summed by south6_run.gather_sum. Pool: the 14,013 label-defined non-targeting cells
  (stage1.npz nt_cells). For a draw of n cells without replacement from the pool:
    C_draw = gather_sum(draw); C_ref = C_pool - C_draw;
    z(.) = SOUTH6's normalisation of one extra pseudobulk: seng_run.z_from_counts against the
    cached TMM reference of SOUTH6's panel [C_f; C_nt; C_pool] (seng_run.tmm_cache_from_panel),
    frozen GTEx mu/sd of the matched ruler genes, SOUTH6's missing mask, scattered into the
    23,485 ruler coordinates with unmatched genes 0 (as south6_run z_of);
    d0 = z(C_draw) - z(C_ref); S2: d0 L_k.
  Floor ratio = ||d0|| / ||gap|| per space. Draws are made once per distinct n among the 1,836
  factors (212 sizes; plus the four gate sizes if any is not among them), 200 draws per size,
  shared by all factors of that size (a factor's floor depends only on n). Generator per size:
  numpy default_rng([20261027, 1, n]), so computation order cannot change a number.
  Reported per space and size: median and P95 (numpy linear) of the floor ratio, and
  m(n) = median over draws of ||d0||^2.

[D3] Claim-2 floor (implements U1).
  Pool: aged half B (results/same/split_idx.npz aged_B, 2,470 cells), never half A.
  Ladder n in {154, 309, 617, 1235, 2470}; 200 draws without replacement per n; generator
  default_rng([20261027, 2, n]). At n = 2,470 every draw is the whole of half B, so the 200 draws
  are identical; this is a property of the ruling and is reported, not changed.
  C_draw = toward_run.sum_rows(Y_aged, draw); z(C_draw) by seng_run.z_from_counts against the
  cached TMM reference of the SAME main panel, frozen mu/sd of all 23,485 ruler genes, the
  panel's missing mask; d0 = z(C_draw) - z(O) with z(O) the panel's O row; sliced per space.
  Fit, per space separately: ordinary least squares of m(n) on (1, 1/n) over the five ladder
  points, unweighted. Residual check: |m(n) - fit(n)| / fit(n) <= 0.10 at all five points, else
  "fit_fails" in that space and no extrapolation there (U1). F_ext = a + b/5832 (if fit passes);
  F_obs = m(2470).

[D4] Corrected ratio (the task's definition).
  Claim 1: CR_i = sqrt(max(0, ||d_i||^2 - m(n_i))) / ||gap||, per space.
  Claim 2: CR_ext = sqrt(max(0, ||d||^2 - F_ext)) / ||gap||; CR_obs likewise with F_obs, per space.
  Raw ratio = ||d|| / ||gap|| as stored.

[D5] Gate (implements the task's Gate and U1).
  Claim 1: sizes 59, 84, 106, 189 (numpy-linear P5/P25/P50/P95 of the 1,836 factor cell counts,
  recomputed by the script; a different value is a reproduction failure). Claim 2: the five
  ladder sizes. Per space and size: 200 fresh control draws (claim 1: default_rng([20261028, 1,
  n]); claim 2: default_rng([20261028, 2, n])), each d0_j built as in D2/D3. For each draw j and
  each L in {0.5, 1, 2, 3}: a direction u drawn uniformly on the unit sphere (normalised standard
  normal) of the space's coordinates that can carry a measurement (S1: matched ruler genes not in
  the missing mask; S2: R^k; claim 2: that space's columns not in the panel missing mask), from
  the same generator, after that draw's cells; x = d0_j + L ||gap|| u;
  CR = sqrt(max(0, ||x||^2 - m(n))) / ||gap||, with m(n) the FLOOR's (D2/D3) median at that n,
  not the gate draws' own.
  Gate passes in a space iff, for every L in {1, 2, 3} and every size tested,
  |median_j CR - L| <= 0.20 L. L = 0.5 is reported and does not gate.
  For claim 2 the gate validates the subtraction at the ladder sizes. It cannot test the
  extrapolation to 5,832; U1's two-floor rule is what guards that.
  Order: the gate is evaluated, and its pass/fail written to file, before any real unit's
  corrected ratio is computed. Real units' RAW ratios are already known (they are the stored
  statistics the task requires reproducing first); this is unavoidable and disclosed.

[D6] Claim-2 gap spread (descriptive; gates nothing).
  200 re-splits, one generator default_rng(20261029); per draw: permute the aged d0 cells
  (4,939) then the young d0 cells (7,703), half A' = first n//2 (2,469; 3,851) as in
  same_run.split_half; gap_j = z(Y_A') - z(O_A'), both normalised by z_from_counts against the
  SAME main panel cache; ||gap_j|| per space. Reported: median, P5, P95, beside the observed
  ||gap||. Not shown to the user: the spread across re-splits measures how much ||gap|| varies
  between splits; it does not measure how much cell-sampling noise inflates ||gap|| itself (every
  split carries that inflation). The findings say so. The corrected ratio uses the observed
  ||gap||, uncorrected, as the task defines it; noise in ||gap|| makes the corrected ratio biased
  low, and the findings state that direction.

[D7] Firing rules (implements the task's Conclusions and U1).
  FLOOR-REACHES-2 fires in a claim x space when the floor P95 is >= 2 at >= 1 size.
    <N> = number of floor sizes: claim 1, the distinct factor sizes (212); claim 2, the 5 ladder
    sizes. Gate-only sizes are not counted. It needs no gate.
  Where the gate passes:
    claim 1: LENGTH-STANDS fires if n >= 1 factors have CR >= 2; LENGTH-DEPTH fires if n >= 1 have
      raw >= 2 and CR < 2. <N> = 1,836. Both may fire. Counts of 0 are written as descriptive lines.
    claim 2 (N = 1): if fit_fails in the space, CR_obs alone decides. Otherwise the unit counts for
      LENGTH-STANDS only if CR_ext >= 2 AND CR_obs >= 2, and for LENGTH-DEPTH only if raw >= 2 AND
      CR_ext < 2 AND CR_obs < 2. If CR_ext and CR_obs fall on opposite sides of 2, neither fires
      and LENGTH-UNRESOLVED fires in that space with this wording, fixed now (not shown to the
      user; the task's fixed LENGTH-UNRESOLVED sentence says "failed its gate", which would be
      untrue here, so U1's "fire LENGTH-UNRESOLVED ... and say that the verdict depends on the
      extrapolation" is implemented as):
      "LENGTH-UNRESOLVED: in claim 2, <space>, the corrected ratio is <x> with the floor
      extrapolated to 5,832 cells and <y> with the observed 2,470-cell floor, on opposite sides of
      2; the verdict depends on the extrapolation and no statement about the published claim is
      made."
  Where the gate fails: LENGTH-UNRESOLVED fires with the task's wording; LENGTH-STANDS and
    LENGTH-DEPTH do not fire; CR is not reported as a result there (computed values, if any, are
    not written as conclusions).
  CLAIM1-*: only in S1, only if the S1 gate passes: HOLDS if the LENGTH-STANDS count is 1,836;
    PARTIAL with <n> if 1 <= n < 1,836; FAILS if 0.
  CLAIM2-*: only in all genes, only if its gate passes: HOLDS iff the unit counts for
    LENGTH-STANDS there; FAILS iff CR < 2 under both floors (or under CR_obs alone if fit_fails);
    if they disagree, neither CLAIM2 key fires and the all-genes LENGTH-UNRESOLVED above stands
    in its place. MD and AGE readings are reported beside it, labelled secondary.
  The claim 1 S1 and claim 2 all-genes readings are primary; all others are secondary.
  Note, not shown to the user: the "2.3" of the main text is the MD-space value; all-genes is
  2.09. CLAIM2 is decided in all genes as the task fixes. If CLAIM2-FAILS fires, the findings
  name which sentences fail, including whether the MD "2.3" sentence survives in its own space.

[D8] Tolerances (the task names none; fixed now).
  Counts exact. Recomputation from stored vectors: relative 1e-9. Rebuilt pseudobulks or rebuilt
  bases (z rows, PCA projections): max abs difference 1e-8 on z, relative 1e-6 on norms and
  ratios. Cache normalisation self-check (z_from_counts of a panel row vs the panel's z):
  max abs 1e-6.

[D9] Rebuilding the S2 basis (prominent: this re-runs a fit).
  SOUTH6 fit its S2 PCA at run time and did not store it. Per the task ("rebuild it with the
  original code and show that the stored statistic reproduces"), it is rebuilt with SOUTH6's
  own steps: toward_run.load_gtex_z, center = column mean, lowdim_common.fit_pca(Zg - center,
  KMAX) with south6_run's KMAX. Inputs to that fit: the 652 GTEx fibroblast donors' frozen-z
  matrix only, no Southard or GSE297234 cell. The rebuilt basis is used only if every
  ||d_f L_k|| / ||(c_old - c_young) L_k|| reproduces the stored S2 d_over_v (D8). south6_run
  .build_spaces is not called, because its Y2 branch writes to results/south6/.

== STOP conditions, fixed now ==

 [S1] `reproduction_failed` — any gating check below fails. Every mismatch named.
 [S2] `flag_ambiguous` — the flag is ambiguous or self-contradictory in a way U1-U3 and D1-D9 do
      not settle. Write the findings with the STOP key, name the ambiguity, ask.
A STOP is never resumed and never restarted with an adjustment.

== Reproduction of stored numbers before anything new (CLAUDE.md) ==

Run first, by src/floorchk_run.py, into results/floorchk/reproduction.csv. All gating.
 R1 Claim 1, S1. From stage1.npz d_f and anchors.npz: factor count 1,836 (exact); per-factor
    ||d_f|| / ||c_old - c_young|| vs stored d_over_v (S1, O1_Y1), rel 1e-9; ||gap|| = stored
    v_norm 32.2839...; min 2.0966 (which rounds to the published 2.10) and 1,836 of 1,836 >= 2,
    both from the recomputed values.
 R2 Claim 1, panel and pool. tmm_logcpm_quiet on [C_f; C_nt; C_pool] from stage1.npz, z as
    SOUTH6: reproduces stored z_f, z_pool, d_f (max abs 1e-8). gather_sum over the 14,013
    nt_cells equals C_pool exactly; gather_sum over each factor's cells equals its C_f row
    exactly; cells per factor equal stage0_cells_per_target (exact). Cache self-check
    z_from_counts(C_pool) vs z_pool (1e-6).
 R3 Claim 1, S2. Rebuilt basis (D9): per factor, ratio vs stored d_over_v for S2_k20/50/100
    (rel 1e-6), and ||gap|| vs stored v_norm (rel 1e-6).
 R4 Claim 2. load_donor for both lines, split_idx.npz halves, state_idx: n O 2,469, Y 3,851,
    aged_B 2,470, S 5,832 (exact); SAME main panel rebuilt with same_run.panel_z: Z vs stored
    results/same/panel_z.npz Z (max abs 1e-8); MD 205 and AGE 3,089 columns (exact); per space
    ||d||, ||gap|| vs results/genespace/same_stats.csv forward d_norm, v_norm (rel 1e-9); raw
    ratios 2.300 (MD), 2.084 (AGE), 2.087 (FULL). Cache self-check on rows O and S (1e-6).

== Honesty obligations (CLAUDE.md), applied ==

- Chance baseline beside counts: the floor P95 is by construction exceeded by 5% of control
  draws; FLOOR-REACHES-2 is about the floor itself. Each LENGTH count is shown beside the floor
  median ratio at the unit's own size.
- Recovery power beside any null: the gate is a planted-effect recovery check; its full table
  (every L, every size, every space) is reported, including L = 0.5.
- Circularity, disclosed where results appear: claim 1, the non-targeting pool defines zero for
  every d and is the pool every floor draw comes from. Claim 2, O is one end of the gap, the
  origin of d and the reference of every floor draw, so O's own sampling noise enters d, d0 and
  the gap alike.
- Cell-count check: Spearman of corrected ratio and of raw ratio with n over the 1,836 factors.
- What this does not cover: the floor measures cell-sampling noise only, not guide-to-guide,
  batch or biological variation; the quadrature subtraction assumes noise is independent of and
  roughly orthogonal to the effect; Hs27 is one neonatal line with no aged cells; claim 2 is two
  donors, one each. Nothing here transfers to other screens or to aged cells.
- No conclusion about any named factor or gene.

== Process ==

Code: src/floorchk_run.py (modes: run, verify). Interpreter: anaconda python, run with -I.
Imports, unchanged: toward_run (tmm_logcpm_quiet, sum_rows, load_donor, load_gtex_z,
ANCHORS_PATH), seng_run (tmm_cache_from_panel, z_from_counts), south6_run (gather_sum,
open_matched, KMAX, S2_KS), lowdim_common.fit_pca, same_run (d0_fibroblast_idx, state_idx,
panel_z), fibro2_common.load_frozen_ruler. Results: results/floorchk/ only. Deliverables: this
file; FINDINGS_FLOORCHK.md with this file reproduced verbatim at the top; PROGRESS_FLOORCHK.md.
One local commit, no push. Not modified: any FINDINGS_*.md, any SCOPING_*.md, anything under
results/south6/ or results/ratio/, any preprint source.
~~~~

## Reproduction of stored numbers, before anything new

`python -I src/floorchk_run.py run`, first stage, writes `results/floorchk/reproduction.csv`. It has
**30 gating checks and 0 mismatches**, so no STOP fired.

| check | expected | observed | match |
|---|---|---|---|
| R1 claim-1 factor count (stage1.npz vs stage3 S1 O1_Y1) | 1836 | 1836 / 1836 | yes |
| R1 S1 per-factor ratio vs stored `d_over_v`, max rel err | ≤1e-9 | 2.6e-15 | yes |
| R1 S1 ‖gap‖ (anchors.npz, read not refit) vs stored v_norm | 32.2839490499 | 32.2839490499 | yes |
| R1 S1 smallest ratio (published "2.10") | 2.0966 | 2.0966 | yes |
| R1 S1 factors with ratio ≥ 2 (published "every one of the 1,836") | 1836 | 1836 | yes |
| R2 SOUTH6 missing mask | equal | 164 vs 164 | yes |
| R2 panel z_f / z_pool / d_f rebuilt vs stored, max abs | ≤1e-8 | 0 / 0 / 0 | yes |
| R2 cache self-check z_from_counts(C_pool) | ≤1e-6 | 0 | yes |
| R2 sum over the 14,013 non-targeting cells = C_pool | exact | exact | yes |
| R2 cells per factor = stage0_cells_per_target | 1836 equal | 1836 equal | yes |
| R2 per-factor cell sums = C_f rows | 0 mismatches | 0 | yes |
| R2 gate sizes = P5/P25/P50/P95 of cells per factor | 59, 84, 106, 189 | 59, 84, 106, 189 | yes |
| R3 S2_k20 / k50 / k100 per-factor ratio, rebuilt basis vs stored, max rel err | ≤1e-6 | 6.5e-16 / 6.1e-16 / 7.8e-16 | yes |
| R3 S2_k20 / k50 / k100 ‖gap‖ vs stored v_norm | 28.645166 / 30.864266 / 31.434263 | same | yes |
| R4 cells O / Y / aged half B / S / aged d0 / young d0 | 2469, 3851, 2470, 5832, 4939, 7703 | same | yes |
| R4 halves disjoint and covering aged d0 | yes | yes | yes |
| R4 SAME main panel z rebuilt vs stored `panel_z.npz`, max abs | ≤1e-8 | 0 | yes |
| R4 space sizes FULL / MD / AGE | 23485 / 205 / 3089 | same | yes |
| R4 FULL ‖d‖, ‖gap‖ vs genespace `same_stats.csv` | ≤1e-9 | 0 (272.7340, 130.6884; ratio 2.0869) | yes |
| R4 MD ‖d‖, ‖gap‖ | ≤1e-9 | 0 (33.5800, 14.5976; ratio 2.3004) | yes |
| R4 AGE ‖d‖, ‖gap‖ | ≤1e-9 | 1.4e-16 (100.7160, 48.3385; ratio 2.0836) | yes |
| R4 published ratios: MD "2.3"; AGE and FULL "> 2" | 2.30 / >2 / >2 | 2.300 / 2.084 / 2.087 | yes |
| R4 cache self-check on rows O and S | ≤1e-6 | 0 / 0 | yes |

The stored minima in the compressed spaces were also read; they are descriptive and not gating. S2_k20 has min 0.1312, with 1 factor ≥ 2; S2_k50 min 0.1940, 2 factors; S2_k100 min 0.2384, 2 factors.

**Claim 1 is about S1.** Paper 1's Supplementary Note 3 calls this space "space C", the
pre-registered primary space of the screen check (log2 CPM, z-scored with the frozen GTEx mu/sd).
SOUTH6 calls it S1. The two are identified by their numbers: 1,836 perturbations, smallest ratio 2.10,
and an origin and target "32 units apart" (‖gap‖ = 32.28). These match S1 O1_Y1 and no other stored
setting.

**The rebuilt S2 basis re-runs a fit.** SOUTH6 never stored its S2 PCA basis. It was rebuilt with
SOUTH6's own steps: `fit_pca` on the 652 GTEx fibroblast donors' frozen-z matrix, centered, with no
Southard or GSE297234 cell. The rebuilt basis reproduces every stored S2 ratio to about 1e-15, so it
is the same instrument. It is still a re-run of a fit, and is stated here as one.

## Claim 1 — Southard, O1 → Y1 (S1 primary)

### Floor: non-targeting cells at each perturbation's own size, with no effect

There are 212 distinct factor sizes, with 200 draws each. Every draw is compared with the
non-targeting pool minus the drawn cells (U3). Tables: `claim1_floor.csv`, all 212 sizes × 4 spaces.

| cells | factors at this size | S1 floor median | S1 floor P95 | S2_k20 median / P95 | S2_k50 median / P95 | S2_k100 median / P95 |
|---|---|---|---|---|---|---|
| 13 (min) | 1 | 11.906 | 12.927 | 2.438 / 2.895 | 3.541 / 4.074 | 4.215 / 4.846 |
| 59 (P5) | 10 | 5.183 | 5.417 | 0.606 / 0.717 | 0.927 / 1.065 | 1.113 / 1.256 |
| 84 (P25) | 18 | 4.279 | 4.433 | 0.461 / 0.559 | 0.689 / 0.798 | 0.832 / 0.940 |
| 106 (P50) | 19 | 3.764 | 3.922 | 0.345 / 0.454 | 0.510 / 0.636 | 0.613 / 0.746 |
| 189 (P95) | 3 | 2.777 | 2.850 | 0.219 / 0.294 | 0.316 / 0.390 | 0.381 / 0.457 |
| 341 (max) | 1 | 2.057 | 2.109 | 0.133 / 0.184 | 0.182 / 0.229 | 0.218 / 0.267 |

**In S1, random non-targeting cells already read more than twice the old-to-young gap at every group
size in the screen.** The floor median runs from 2.06 at 341 cells to 11.9 at 13 cells. At the median
perturbation size (106 cells) it is 3.76, against the perturbations' stored median of 3.82. A ratio of
2 in S1 therefore says nothing about any perturbation. Only 406 of 1,836 factors (22%) exceed the P95 of
their own-size floor, against 91.8 expected by chance (5%). The rest are inside the noise. Raw ratio
tracks depth almost perfectly (Spearman with cells −0.987). After the correction that dependence is
gone (−0.091).

In the compressed spaces the floor is far lower. It reaches 2 only at the smallest sizes: 13 cells in
S2_k20, 13 to 28 cells in S2_k50, and 13 to 29 cells in S2_k100. Factors exceeding their own-size floor
P95 are 367, 301 and 288 of 1,836, against 91.8 expected.

### Gate: planted displacements recovered by the correction

The gate was written to disk before any real unit's corrected ratio was computed. It passes in all four
spaces, with 0 failing cells of 12 in each. The worst gating cell is S1 at 106 cells with planted length
1: median 0.960, or −4.0%. Every other gating cell is within ±2%. The full table is in
`claim1_gate.csv`; the S1 rows are below.

| S1, cells | planted 0.5 (not gating) | 1 | 2 | 3 |
|---|---|---|---|---|
| 59 | 0.452 [IQR 0.000–1.168] | 0.988 [0.385–1.462] | 1.975 | 3.001 |
| 84 | 0.462 [0.000–0.839] | 0.982 [0.650–1.211] | 1.975 | 2.985 |
| 106 | 0.391 [0.000–0.860] | 0.960 [0.738–1.197] | 1.973 | 2.981 |
| 189 | 0.507 [0.290–0.642] | 1.000 [0.899–1.077] | 2.001 | 2.999 |

What the gate does and does not show:

- In S1 the median is right, but single values are wide at small depth. At 59 cells, a planted
  length of 1 reads anywhere from 0.39 to 1.46 across the middle half of draws.
- Planted 0.5 reads low (−22% at 106 cells), because the max(0, ·) clips many draws to 0. Short
  effects are under-read.
- The gate plants effects in uniformly random directions. A real effect that is not orthogonal to
  the noise is not covered by it.

### Corrected ratio of the 1,836 perturbations

| space | raw ratio: min / median | raw ≥ 2 | corrected: median / max | corrected ≥ 2 | corrected = 0 |
|---|---|---|---|---|---|
| **S1** | 2.097 / 3.821 | 1836 | **0.607 / 2.961** | **18** | 551 |
| S2_k20 | 0.131 / 0.412 | 1 | 0.149 / 0.947 | 0 | — |
| S2_k50 | 0.194 / 0.585 | 2 | 0.165 / 1.569 | 0 | — |
| S2_k100 | 0.238 / 0.692 | 2 | 0.182 / 1.670 | 0 | — |

The S1 corrected-ratio quantiles are: P10 0.000, P25 0.000, P50 0.607, P75 0.981, P90 1.276, P99 1.99.
Per factor, per space: `claim1_units.csv`. Summary: `claim1_summary.csv`.

**The chance baseline for the 18 is not pre-registered and gates nothing.** It comes from
`python -I src/floorchk_run.py diagnostic` and is written to `claim1_diagnostic.csv`. Each floor draw is
treated as a no-effect unit of its size and given the same correction. If all 1,836 factors were pure
noise at their own sizes, about **5.05** would reach a corrected ratio ≥ 2 in S1 (S2: 0.005, 0.055 and
0.16). The 18 that do sit beside that 5. Seventeen of them exceed the P95 of no-effect corrected ratios
at their own size.

**Cell-count check (CLAUDE.md).** The 18 are smaller groups than the screen as a whole: their median is
54.5 cells, ranging from 13 to 98, against 106 for all factors. Small groups are where the corrected
ratio is least precise. The no-effect corrected ratio has P95 1.58 at 59 cells, and P95 ≥ 2 at 16 of
the 212 sizes, all of 47 cells or fewer. The 18 should not be read as a shortlist. This task names no factor.

### Readings, claim 1 (pre-registered wording; `results/floorchk/readings.txt`)

```
== claim 1 / S1 (primary) ==
CONCLUSION  FLOOR-REACHES-2: in claim 1, S1, the control floor's P95 is >= 2 at 212 of 212 group sizes; at those depths a ratio of 2 is reached by sampling noise alone.
CONCLUSION  LENGTH-STANDS: in claim 1, S1, 18 of 1836 units have corrected ratio >= 2; the published 'more than twice' holds for them after removing the noise contribution.
CONCLUSION  LENGTH-DEPTH: in claim 1, S1, 1818 of 1836 units have raw ratio >= 2 but corrected ratio < 2; for them the published 'more than twice' was reached partly by sampling noise.
CONCLUSION  CLAIM1-PARTIAL: 18 of 1,836
== claim 1 / S2_k20 (secondary) ==
CONCLUSION  FLOOR-REACHES-2: in claim 1, S2_k20, the control floor's P95 is >= 2 at 1 of 212 group sizes; at those depths a ratio of 2 is reached by sampling noise alone.
            (descriptive) 0 of 1836 units have corrected ratio >= 2; LENGTH-STANDS did not fire.
CONCLUSION  LENGTH-DEPTH: in claim 1, S2_k20, 1 of 1836 units have raw ratio >= 2 but corrected ratio < 2; for them the published 'more than twice' was reached partly by sampling noise.
== claim 1 / S2_k50 (secondary) ==
CONCLUSION  FLOOR-REACHES-2: in claim 1, S2_k50, the control floor's P95 is >= 2 at 3 of 212 group sizes; at those depths a ratio of 2 is reached by sampling noise alone.
            (descriptive) 0 of 1836 units have corrected ratio >= 2; LENGTH-STANDS did not fire.
CONCLUSION  LENGTH-DEPTH: in claim 1, S2_k50, 2 of 1836 units have raw ratio >= 2 but corrected ratio < 2; for them the published 'more than twice' was reached partly by sampling noise.
== claim 1 / S2_k100 (secondary) ==
CONCLUSION  FLOOR-REACHES-2: in claim 1, S2_k100, the control floor's P95 is >= 2 at 4 of 212 group sizes; at those depths a ratio of 2 is reached by sampling noise alone.
            (descriptive) 0 of 1836 units have corrected ratio >= 2; LENGTH-STANDS did not fire.
CONCLUSION  LENGTH-DEPTH: in claim 1, S2_k100, 2 of 1836 units have raw ratio >= 2 but corrected ratio < 2; for them the published 'more than twice' was reached partly by sampling noise.
```

The claim-1 gap is fixed. c_old and c_young are means over 231 and 113 GTEx donor bulk profiles: frozen
anchors, read and not refit. They are not single-cell pseudobulks and carry no cell-sampling noise.
Their donor-sampling uncertainty is a different quantity and was not estimated.

## Claim 2 — GSE297234 same-platform test (all genes primary)

### Floor: the aged donor's half-B day-0 cells, on the ladder (U1)

The ladder draws 200 samples at each size, without replacement, from the 2,470 half-B cells. At 2,470
cells every draw is the whole of half B, so that rung is a single realisation (`n_distinct_draws` = 1).
Each draw's d0 = z(draw) − z(O), with O the stored half A. Table: `claim2_floor.csv`.

| cells | all genes: floor median / P95 | MD: median / P95 | age: median / P95 |
|---|---|---|---|
| 154 | 1.621 / 1.683 | 0.542 / 0.822 | 1.385 / 1.484 |
| 309 | 1.172 / 1.203 | 0.365 / 0.499 | 0.967 / 1.026 |
| 617 | 0.853 / 0.868 | 0.251 / 0.375 | 0.680 / 0.715 |
| 1,235 | 0.640 / 0.648 | 0.187 / 0.231 | 0.494 / 0.513 |
| 2,470 (one realisation) | 0.518 | 0.146 | 0.388 |

No ladder size reaches a floor P95 of 2 in any space. At 154 cells the all-genes floor is already 1.62.

### The 1/n fit and the extrapolation to 5,832 cells (U1)

| space | a | b | relative residuals at 154 / 309 / 617 / 1,235 / 2,470 | fit passes (≤10%) | floor² read at 5,832 (F_ext) | floor² observed at 2,470 (F_obs) |
|---|---|---|---|---|---|---|
| **all genes** | 1756.98 | 6,651,156 | −0.001 / +0.007 / −0.009 / −0.019 / +0.028 | **yes** | 2897.44 (floor ratio 0.412) | 4575.25 (0.518) |
| MD | −0.866 | 9,595 | +0.018 / −0.062 / −0.086 / +0.076 / **+0.502** | no | (not used) | 4.533 (0.146) |
| age | 18.04 | 682,877 | +0.006 / −0.018 / −0.039 / −0.003 / **+0.197** | no | (not used) | 352.67 (0.388) |

**The all-genes value at 5,832 cells is an extrapolation,** about 2.4 times (5,832 / 2,470 = 2.36) past
the largest size drawn. In MD and age the fit fails its 10% check, so nothing was extrapolated there.
The observed 2,470-cell floor is used alone, as U1 requires. In both, the failure is at the 2,470-cell
rung, which is a single fixed set of cells rather than a median over 200 draws. That rung sits above the
curve, so using it overstates the noise and works against the claim. Fit details: `claim2_fit.csv`.

### Gate

The gate passes in all three spaces, with 0 failing cells of 15 in each. Every gating cell is within
±0.7% of the planted length, and planted 0.5 within ±1.9%. Table: `claim2_gate.csv`. The gate checks
the subtraction at the ladder sizes. It cannot check the extrapolation to 5,832; U1's rule of requiring
both floors is the guard there.

### Gap spread (descriptive; D6)

There were 200 re-splits of both donors' day-0 cells (seed 20261029). Table: `claim2_gap_spread.csv`.

| space | observed ‖gap‖ | re-split median [P5, P95] | (P95 − P5) / observed |
|---|---|---|---|
| all genes | 130.688 | 130.667 [129.418, 131.760] | 1.8% |
| MD | 14.598 | 14.776 [14.493, 15.061] | 3.9% |
| age | 48.338 | 48.756 [48.115, 49.378] | 2.6% |

This spread says how much ‖gap‖ moves between splits. It does **not** say how much cell-sampling noise
inflates ‖gap‖ itself, because every split carries that inflation. The corrected ratio divides by the
observed, uncorrected ‖gap‖, as the task defines it. Any noise inflation of the gap therefore biases the
corrected ratio **downward**, against the claim.

### Corrected ratio of the partially reprogrammed population (5,832 cells)

| space | raw ratio | corrected, extrapolated floor | corrected, observed 2,470 floor | ≥ 2 under the rule |
|---|---|---|---|---|
| **all genes (23,485)** | 2.0869 | **2.0459** | **2.0217** | **yes, under both** |
| MD genes (205), secondary | 2.3004 | — (fit failed) | 2.2957 | yes |
| age genes (3,089), secondary | 2.0836 | — (fit failed) | 2.0470 | yes |

The all-genes margin is thin. This is arithmetic on `claim2_unit.csv`, not a pre-registered test. For
the all-genes corrected ratio to fall below 2, the noise in d would have to be 6066 in squared length.
That is 2.09 times the extrapolated floor and 1.33 times the observed 2,470-cell floor. In the age
genes the noise would need to be 2.26 times the observed floor, and in the MD genes 61 times.

**The floor comes from day-0 fibroblasts, not from partially reprogrammed cells.** If partially
reprogrammed cells vary more from cell to cell, their real sampling noise at 5,832 cells is larger than
this floor. A 33% excess over the observed 2,470-cell floor in all genes would be enough to reverse the
reading. This task cannot measure it: there is no no-effect pool of partially reprogrammed cells.

### Readings, claim 2 (pre-registered wording)

```
== claim 2 / all genes (23,485) (primary) ==
            (descriptive) floor P95 >= 2 at 0 of 5 group sizes; FLOOR-REACHES-2 did not fire.
            corrected ratio 2.0459 (floor extrapolated to 5,832) and 2.0217 (observed 2,470 floor); raw ratio 2.0869
CONCLUSION  LENGTH-STANDS: in claim 2, all genes (23,485), 1 of 1 units have corrected ratio >= 2; the published 'more than twice' holds for them after removing the noise contribution.
CONCLUSION  CLAIM2-HOLDS
== claim 2 / MD genes (205) (secondary) ==
            (descriptive) floor P95 >= 2 at 0 of 5 group sizes; FLOOR-REACHES-2 did not fire.
            corrected ratio 2.2957 (observed 2,470 floor alone: the 1/n fit failed its 10% residual check, so no extrapolation); raw ratio 2.3004
CONCLUSION  LENGTH-STANDS: in claim 2, MD genes (205), 1 of 1 units have corrected ratio >= 2; the published 'more than twice' holds for them after removing the noise contribution.
== claim 2 / age genes (3,089) (secondary) ==
            (descriptive) floor P95 >= 2 at 0 of 5 group sizes; FLOOR-REACHES-2 did not fire.
            corrected ratio 2.0470 (observed 2,470 floor alone: the 1/n fit failed its 10% residual check, so no extrapolation); raw ratio 2.0836
CONCLUSION  LENGTH-STANDS: in claim 2, age genes (3,089), 1 of 1 units have corrected ratio >= 2; the published 'more than twice' holds for them after removing the noise contribution.
```

`python -I src/floorchk_run.py verify` re-derives every reading line and both gate summaries from the
tracked CSVs alone. It reports 30 reading lines and 0 mismatches.

## What no longer holds in paper 1 (CLAIM1-PARTIAL fired)

The preprint is not edited. The affected sentence is in Supplementary Note 3, in the paragraph beginning
"Space C remains the screen's pre-registered primary". Paraphrased: in space C, from the GTEx old
centroid to the GTEx young centroid, every one of the 1,836 perturbations moved more than twice the
distance between them (smallest ratio 2.10), so no approach could have been registered whatever its
direction. The paragraph adds that this is the same limitation the main text discloses for the
same-platform test.

**What fails as written.** That the perturbations *moved* more than twice the distance. The measured
displacement vectors are indeed all longer than twice the gap. Inside this geometry, a measured delta
from that origin therefore could not have come out negative; that part is still literally true of the
measurement. But non-targeting cells with no effect, drawn at each perturbation's own group size,
already read 2.06 to 11.9 times the gap. That is above 2 at every one of the screen's 212 group sizes.
After the noise contribution is removed in quadrature, 18 of 1,836 perturbations are estimated to have
moved more than twice the gap, about 5 would be expected from noise alone, and the median estimate is
0.61 gap lengths. The test was blind because of the screen's depth (13 to 341 cells per perturbation in
a 20,264-gene z space), not because the perturbations moved too far. The closing comparison also fails
as written. The same-platform length statement survives this check (CLAIM2-HOLDS), but the screen
statement rests on sampling noise, so the two are not the same limitation.

**A corrected statement could read:** "In space C, every one of the 1,836 measured displacements is
longer than twice the distance between the GTEx old and young centroids (smallest 2.10), so from that
origin the test could not register an approach whatever the direction. This reflects depth, not effect
size: non-targeting cells with no effect, sampled at the same group sizes, also read more than twice
that distance (2.1 to 12.9 at the 95th percentile). After the noise is subtracted, 18 perturbations are
estimated to move more than twice the distance (about 5 expected by chance), and the median estimate is
0.6 times it. The space C null for the screen is therefore a null of an instrument that could not see
approach at this depth, not evidence that the perturbations overshoot."

**What still holds.** Claim 2, the main-text "2.3 times" sentence (MD genes, corrected 2.30), and the
Limitations bullet that the same-platform displacement exceeds twice the donor gap (all genes, corrected
2.05 and 2.02). The all-genes margin is about 1 to 2%, as set out above.

## Deviations and decisions not in the task (disclosed)

1. **Three questions went to the user before the flag was written.** They covered the claim-2 pool
   being smaller than the unit, the missing seeds, and the claim-1 reference. The rulings (U1 to U3)
   are reproduced verbatim in the flag. CLAUDE.md contains no seeds, contrary to the task's wording.
   The seeds used are the ones the user ruled on: 20261027, 20261028 and 20261029.
2. **Outcome knowledge before the flag.** I had seen the stored ratios. From the stored GENESPACE f = 0
   mixture row, a hand calculation put the all-genes corrected ratio "within a few percent of 2". This
   was said to the user before U1 was ruled on, and it is recorded in the flag. U1's choice was made
   knowing the floor-size choice could decide CLAIM2.
3. **Definitions D1 to D9 were fixed by the session.** Four were shown to the user before the rulings;
   the rest were not. The flag marks which. One wording was added in D7: a LENGTH-UNRESOLVED variant
   for floor disagreement in claim 2, because the task's own sentence says "failed its gate". It did
   not fire.
4. **Claim 2 has one unit**, GM00731 PartialReprog. The task says "state (claim 2)". No other state was
   tested.
5. **The S2 basis was rebuilt** (D9). This is a re-run of a fit, and it reproduces the stored ratios to
   about 1e-15. `south6_run.build_spaces` was not called, because its Y2 branch writes into
   results/south6/. Its PCA lines were repeated with the imported `fit_pca`, `KMAX` and `load_gtex_z`.
6. **Claim-1 floor draws are shared by size.** All 212 sizes got 200 draws each, used by every factor
   of that size. The four gate sizes were all factor sizes, so no extra floor size was drawn.
7. **Claim-2 fit.** In MD and age the fit failed at the single-realisation 2,470 rung, so only the
   observed floor was used there (U1). In all genes both floors were used, and they agree.
8. **A diagnostic was added after the readings and is not pre-registered.** `diagnostic` mode, writing
   `claim1_diagnostic.csv`, gives the chance baseline beside the LENGTH-STANDS count (CLAUDE.md). It
   gates nothing. The noise-margin arithmetic for claim 2, and the exceed-floor counts (406, 367, 301,
   288; derived from `claim1_summary.csv` `n_raw_within_floor_p95`), are descriptive and were computed
   from tracked CSVs.
9. **Interpreter and smoke checks.** The interpreter is Anaconda Python 3.12 (`numpy` 2.x), run with
   `-I`. The default `python` on PATH has no numpy. Before the run, the module was imported once
   (`python -c "import floorchk_run"` from `src/`, without `-I`) to check that the imports resolved.
   Nothing was computed. The run completed in one pass, about 11 minutes, with no crash and no rerun.
10. **Where things were written.** Results went to `results/floorchk/` only. `run_stdout.txt`
    duplicates `floorchk_run.log`. One line was added to `.gitignore` for
    `results/floorchk/claim1_floor_draws.npz` (1.26 MB, regenerable from seeds). `claim1_units.csv`
    (1.06 MB) is a result CSV and is tracked. Nothing under results/south6/, results/ratio/,
    results/same/ or results/genespace/ was modified. No FINDINGS_*, SCOPING_* or preprint file was
    modified. The pre-existing uncommitted change to `CITATION.cff` and the untracked
    `results/paper2_figs_plain/` are not this task's work, and are left out of its commit.

## Honesty obligations

- **Chance baseline beside counts.** 18 corrected ratios ≥ 2 sit beside about 5.05 expected from
  noise. 406 floor-exceeding factors in S1 sit beside 91.8. FLOOR-REACHES-2 is a statement about the
  floor itself.
- **Recovery power.** The gate is the planted-effect check, and it is reported in full. In S1 at small
  depth, single-unit corrected ratios are wide: at 59 cells, a planted length of 1 gives an IQR of 0.39
  to 1.46. Lengths of 0.5 are under-read.
- **Circularity, disclosed and not refitted.** In claim 1, the non-targeting pool defines zero for every
  d and supplies every floor draw. In claim 2, O (aged half A) is the origin of d, one end of the gap
  and the reference of every floor draw, so O's sampling noise enters all three.
- **Assumption.** The quadrature subtraction assumes noise is independent of, and roughly orthogonal
  to, the effect. The gate supports this for random directions only.
- **What the floor covers.** Cell-sampling noise only. It does not cover guide-to-guide, batch or
  biological variation. In claim 2 it does not cover any excess heterogeneity of partially
  reprogrammed cells over day-0 fibroblasts.
- **What this does not transfer to.** Hs27 is one neonatal line with no aged cells. Claim 2 is one aged
  and one young donor. Nothing here transfers to aged cells in general, to other screens, or to other
  origins and targets.
- **No factor or gene is named or discussed.** The per-factor CSV carries labels because it is data.

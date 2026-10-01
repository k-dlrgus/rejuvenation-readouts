# FINDINGS_GENESPACE — the SAME and TOWARD tests inside smaller gene spaces

**Status:** DONE. Run valid (protected files unchanged): True. Seeds 20260914 / 20260918. n_perm=200, n_random=200, n_boot=200, n_baseline=200.

Reproduced by `src/genespace_run.py`. Only the section-7 reading at the end is the conclusion; every other number is **diagnostic**. No manuscript is edited; paper/ and paper_package/ are not read.

## Pre-registration (verbatim copy of results/genespace/PREREG.flag)

~~~~text
Written **before any GENESPACE statistic**, 2026-09-30. No threshold in this block is re-tuned after numbers exist.

Question: do the paper's two displacement tests (SAME: young target = GM23815 day-0 cells from the same experiment; TOWARD: young/old = GTEx centroids) give the same answer inside smaller gene spaces: (A) the MD genes, (B) the published age-up ∪ age-down genes? The full frozen-ruler space (23,485 genes) is rerun only as a reproduction check.

Seeds: 20260914 (split check, mixtures, gene permutations, donor splits, C3 folds, size baseline). 20260918 (bootstrap). n_perm=200, n_random=200, n_boot=200, n_baseline=200. GSE325735 out of scope. Uncalibrated R² never a gate. No manuscript is edited. paper/ and paper_package/ are not read for any number; reference numbers come only from results/same/ and results/toward/.

STATUS AT WRITING
Read: FINDINGS_SAME.md, FINDINGS_TOWARD.md, src/same_run.py, src/toward_run.py, both PREREG.flag files, results/md3/genesets.json (+ src/md3_genesets.py, genesets_summary.json, fibro2_common.align_counts_to_ruler for the paper's symbol-matching convention). Computed: only the symbol -> ruler-column counts below, from the frozen ruler's `symbol`/`ensembl` arrays. Not loaded: any expression matrix, GTEx mu/sd, ruler weights, results/same/summary.json, results/toward/stats.csv. No statistic exists in any gene space.

INPUTS (read from disk; refit/re-label/re-cluster/rescore nothing)
- Everything under INPUTS in results/same/PREREG.flag and results/toward/PREREG.flag, unchanged (GSE297234 counts and md2 labels; GTEx Stage 0 pack; frozen ruler mu/sd; md4 files).
- results/same/split_idx.npz (O, Y, aged_B, young_B): the paper's half-A / half-B cells.
- results/md3/genesets.json keys "MD", "age_up", "age_down". Keys TGFB and reprog are not used.
- Stop-rule references: results/same/summary.json, results/toward/stats.csv. Read by the script at the check step only.

1. GENE SPACES
Mapping rule (frozen): strip and upper-case each list entry; keep repeated entries once (first occurrence); an entry maps iff it exactly equals an upper-cased entry of the frozen ruler's `symbol` array (results/fibro/frozen_ruler_ridge_raw.npz). No alias table, no HGNC update, no Ensembl lookup, no fuzzy match. Unmapped entries are listed in gene_map.csv and dropped, not rescued. If a symbol matched more than one ruler column, all its columns would be included and counted (13 ruler symbols span >1 column; none is in these lists).
Counts at writing (n mapped / n in list):
  MD        : 205 / 205 entries (205 unique; 0 unmapped)        -> 205 columns
  age_up    : 1312 / 1533 entries (1532 unique; GOLGA8M listed twice; 220 unmapped)
  age_down  : 1777 / 2007 entries (2006 unique; TMEM191A listed twice; 229 unmapped)
  age_up ∩ age_down: 0 in the lists, 0 among mapped. No overlap rule is needed.
  AGE = age_up ∪ age_down: 3089 / 3538 unique                   -> 3089 columns
  MD ∩ AGE (mapped): 43 genes (18 age_up, 25 age_down). Spaces overlap; nothing is removed to make them disjoint.
  FULL: all 23,485 ruler columns.
The script recomputes these counts; if any differs, STOP.
A space is an index set of ruler columns, kept in ruler order (ascending column index). The space's matrix = those columns of the SAME frozen-z matrix the paper builds in the full space: TMM, log2-CPM prior.count=2, and z with frozen GTEx mu/sd run on all 23,485 genes of the paper's panel exactly as in the paper; columns are taken afterwards. TMM and z-scoring are never rerun on a subset. Genes set to z=0 (zero counts in every panel row) stay 0; report n such genes per space.
Space order everywhere: FULL, MD, AGE.

2. SAME-PLATFORM TEST (per space)
Cells: O, Y, aged_B, young_B from results/same/split_idx.npz. The script also reruns split_half(d0, 20260914) per donor and STOPs if it does not equal the file. Mixtures: mix_indices(aged_B, young_B, 20260914), f = 0, 0.10, 0.25, 0.50, total = min(n_aged_B, n_young_B). PartialReprog, Pluripotency, NonReprog pooled across timepoints, md2 labels, as in the paper.
Panels as in src/same_run.py (full-gene TMM): mixture panel [O, Y, mix f=0 … 0.50]; main panel [O, Y, S, S_rev, Pluri, NonReprog].
Inside each space X (z restricted to X):
  forward : v = z(Y) − z(O); S = GM00731 PartialReprog; d = z(S) − z(O).
  reverse : v_rev = z(O) − z(Y); S_rev = GM23815 PartialReprog; d = z(S_rev) − z(Y).
  mixture control: each mixture as S in the forward frame of X.
  context : GM00731 Pluripotency and GM00731 NonReprog as S in the forward frame of X.
  Asymmetry = progress_forward − progress_reverse (reported; not part of the section-7 reading).
  N1 (reported, not a gate here): 200 permutations of v (resp. v_rev) within X, fresh Generator(20260914) per call; p = (#{cos ≥ observed} + 1)/201.
Anchors, control and reverse are built from X's own columns in every space; nothing is copied between spaces. If ‖v‖ < 1e-12 in X, STOP for X.

3. GTEx-ANCHORED TEST (per space; secondary)
Anchors as in src/toward_run.py: GTEx Stage 0 matrix, frozen mu/sd z; YOUNG = AGE 20-29 ∪ 30-39, OLD = 60-69 ∪ 70-79 (report actual n; middle reported, not dropped). c_young(X), c_old(X) = centroids restricted to X; w(X) = c_young(X) − c_old(X); u(X) = w(X)/‖w(X)‖, re-normalized inside X.
Per donor (GM00731 and GM23815 separately, never pooled), for every destination state (PartialReprog, EarlyPluripotency, Pluripotency, NonReprog; MIN_CELLS=30), with the paper's per-donor full-gene panel TMM: d = z(S) − z(d0 Fibroblast), restricted to X. Report cos_S = d·u(X)/‖d‖, dist_young, delta_young, dist_old, delta_old, and the section-4 statistics against w(X).
Rebuilt inside X:
  N1: 200 permutations of u(X) entries, fresh Generator(20260914) per space.
  N2: 200 equal-size random GTEx donor splits ignoring age, fresh Generator(20260914) per space (same splits as the paper); each split's mean difference restricted to X, then unit-normalized.
  C3 (five-fold positive control): same folds (Generator(20260914) shuffle, fold = i % 5); per fold rebuild c_young, c_old, u on training donors restricted to X; cos_fold = d·u_train/‖d‖ with d = c_young_train − mean z(held-out OLD), in X. Pass = all 5 folds finite and every cos_fold > 0.50; construction minimums unchanged (20/20/5).
If C3 fails in X, the GTEx test for X is `c3_broken`: nothing about GSE297234 is reported for X in this test.
If C3 passes, the decision order of results/toward/PREREG.flag is applied unchanged per donor in X as a secondary key (GM23815 contrast only). It never changes the section-7 reading.

4. PER-STATE STATISTICS (every state, every space, both tests)
  cos       = d·v / (‖d‖ ‖v‖)
  progress  = d·v / ‖v‖²
  delta     = ‖z(S) − target‖ − ‖z(origin) − target‖   (negative = closer to target)
  perp      = ‖d − progress·v‖ / ‖v‖                   (sideways part, in units of ‖v‖)
  rel_delta = delta / ‖z(origin) − target‖             (dimensionless)
Same-platform test: v = target − origin, so ‖z(origin) − target‖ = ‖v‖ and
  final distance² / ‖v‖² = (1 − progress)² + perp²,
  hence delta < 0 iff (1 − progress)² + perp² < 1.
The script checks this identity for every state, space and bootstrap draw; relative mismatch > 1e-8 is a code bug: STOP.
GTEx test: v = w(X) = c_young − c_old, but the origin is the donor's d0 pseudobulk, not c_old, so the identity does NOT hold. progress and perp are reported against w(X); delta_young and delta_old are computed directly, never derived from progress/perp.

5. READOUT-VECTOR COSINE (context only; never a gate or reading)
  r_MD(X) : −1 on each column of X whose gene is in MD, 0 elsewhere. Signed so "MD goes down" is a positive cosine.
  r_AGE(X): −1 on age_up columns, +1 on age_down columns, 0 elsewhere (a gene in both lists would get 0; there are none). Signed so "younger" is a positive cosine.
Report cos(d, r) for every state in both tests, plus cos(v, r) and cos(u(X), r), in every space where r has ≥1 nonzero entry, with its n nonzero (FULL: r_MD 205, r_AGE 3089; MD: r_MD 205, r_AGE 43; AGE: r_MD 43, r_AGE 3089).
This approximates the published AddModuleScore and does not reproduce it. AddModuleScore is per cell on log-normalized counts, subtracts expression-bin-matched control genes, and is summarized per cluster; here it is an unweighted sign vector dotted with frozen-GTEx z of TMM pseudobulk displacements.

6. SIZE BASELINE (context; not a gate; never changes a reading)
Cosines and distances are not compared across spaces of different size. For MD (205) and AGE (3089) separately:
  - Bins: 24 equal-count bins of frozen GTEx mu over all 23,485 ruler genes; edges = quantiles of mu at 0, 1/24, …, 1; bin = searchsorted(edges, mu, side="right") − 1, clipped to 0..23.
  - 200 draws. In each draw, for bins 0..23 in order, draw (n space genes in that bin) genes without replacement from all ruler genes in that bin. Space genes are not excluded; report the mean overlap of random sets with the real space.
  - One Generator(20260914) per space, MD before AGE. Drawn column indices are saved.
  - Each random set gets point estimates only (no bootstrap, no nulls), sliced from the same stored full-gene z matrices: every section-4 statistic in both tests, mixture-control progress and delta, and the section-7 geometry label.
  - For the real space report position = (#{random < real} + 0.5·#{random = real}) / 200, the random 2.5 / 50 / 97.5 percentiles, and the count of random sets per geometry label.
Cross-space comparisons use progress, perp, rel_delta, cos and baseline position; never raw distances.

7. PRE-REGISTERED READING (per space, MD and AGE; forward same-platform test; only what fired)
Gate = that space's mixture control passes: progress non-decreasing in f (f=0 ≤ 0.10 ≤ 0.25 ≤ 0.50), AND the progress CI at f=0.50 excludes 0 (lo > 0 or hi < 0, as in results/same/PREREG.flag) AND is valid (section 9), AND delta at f=0.50 < 0.
Decision order:
1. Stop rule (section 8) failed -> key `reproduction_failed`. No space is read.
2. Gate failed -> NO READING.
3. delta < 0 -> APPROACH.
4. delta > 0 AND progress > 1 AND (progress − 1)² > perp² -> OVERSHOOT.
5. delta > 0 AND not OVERSHOOT -> SIDEWAYS.
6. delta = 0 exactly -> NO READING (tie).
Point estimates decide the key. The delta and progress CIs, with valid/invalid flags, are printed after the sentence and do not change the key.
Exact sentences (fill braces; no other wording):
APPROACH: "In the {space} gene space ({n_genes} ruler genes), aged partially-reprogrammed cells (GM00731 PartialReprog, n_cells={n_S}) end closer to the young donor's day-0 cells than they started: change in distance {delta} ({rel_delta} of the starting gap), progress {progress}, sideways {perp}. The reverse test in the same space gives progress {progress_rev} and change in distance {delta_rev}; approach alone does not separate rejuvenation from donor convergence. One young donor, so 'young' is confounded with that donor."
OVERSHOOT: "In the {space} gene space ({n_genes} ruler genes), aged partially-reprogrammed cells (GM00731 PartialReprog, n_cells={n_S}) move past the young donor's day-0 cells along the aged-to-young axis (progress {progress} > 1) and end farther from them than they started (change in distance {delta}); the overshoot term (progress − 1)² = {a} is larger than the sideways term perp² = {b}. This is not approach to the young donor."
SIDEWAYS: "In the {space} gene space ({n_genes} ruler genes), aged partially-reprogrammed cells (GM00731 PartialReprog, n_cells={n_S}) end farther from the young donor's day-0 cells than they started (change in distance {delta}), mainly through movement off the aged-to-young axis: sideways term perp² = {b}, overshoot term (progress − 1)² = {a}, progress {progress}. This is not approach to the young donor."
NO READING: "In the {space} gene space ({n_genes} ruler genes), the mixture positive control failed ({failed_conditions}); the statistic cannot detect a known approach in this space, and its forward numbers are not interpreted."
NO READING (tie): "In the {space} gene space ({n_genes} ruler genes), change in distance is exactly 0; no reading."
The same geometry label (APPROACH / OVERSHOOT / SIDEWAYS) is printed for reverse, Pluripotency and NonReprog in every space, and for forward in FULL, as context only. FULL is a reproduction check and gets no reading.
Spaces are not averaged, ranked, or combined. There is no "best space". MD and AGE readings are reported side by side.

8. STOP RULE (before any MD or AGE statistic)
FULL is run first. Compare to disk, |new − reference| ≤ 1e-6 absolute:
  results/same/summary.json : frac_S, delta_S, frac_rev, p.fracs (all four f).
  results/toward/stats.csv  : row cell_line=GM00731, state=PartialReprog: cos_S, delta_young.
Any miss -> STOP: key `reproduction_failed`; write every compared value, reference and difference to repro_check.json; compute no MD or AGE statistic. Nothing is re-tuned to make it match.

9. INTERVALS
Bootstrap B=200, seed 20260918, with the paper's procedure and generator order, so resampled cells/donors are identical in every space:
  - same-platform mixture: fresh Generator(20260918); O and Y frozen; mixture cells resampled.
  - same-platform main   : fresh Generator(20260918); O and Y frozen; S, S_rev, Pluri, NonReprog resampled.
  - GTEx per donor       : Generator(20260918 + 0) GM00731, Generator(20260918 + 1) GM23815; origin and destinations resampled.
  - C3 fold CIs          : fresh Generator(20260918) per space; held-out OLD donors resampled.
Each draw re-TMMs the full-gene panel and z-scores with frozen mu/sd once; every space is sliced from that same draw. 95% percentile CIs.
Validity: a CI is valid iff lo ≤ point ≤ hi. Otherwise it is flagged `invalid` in every table and text. An invalid CI cannot satisfy any condition (it can only fail the section-7 gate). FINDINGS_SAME.md and FINDINGS_TOWARD.md already note that the re-TMM bootstrap need not contain the point; the procedure is not changed here to make CIs valid.

10. OUTPUTS / FILES
New files only: results/genespace/ (this flag; manifest.json; gene_map.csv; spaces.json; repro_check.json; same_stats.csv; same_posctrl.csv; gtex_stats.csv; gtex_c3.csv; readout_cos.csv; size_baseline.csv; size_baseline_sets.npz; boot.npz; reading.json; summary.json; run_report.txt), src/genespace_run.py, FINDINGS_GENESPACE.md, PROGRESS_GENESPACE.md.
No existing file is modified. src/genespace_run.py may import pure functions from src/same_run.py and src/toward_run.py; it must not call any function that writes (freeze_anchors, build_null_directions, run_c3, run_posctrl, write_findings, write_progress, save_manifest, record_failure of those modules). SHA-256 of every file under results/same/, results/toward/, results/md3/, plus results/fibro/frozen_ruler_ridge_raw.npz, src/same_run.py, src/toward_run.py, FINDINGS_SAME.md, FINDINGS_TOWARD.md is recorded before and after the run; any change -> run flagged invalid.
This flag is not rewritten by the script.

WHAT DOES NOT COUNT
- Rerunning TMM or z-scoring on a subset. Refitting mu/sd, the ruler, labels, clusters or AddModuleScore.
- Adding, removing or re-mapping genes, or choosing among spaces, after numbers exist.
- Comparing raw cosines or distances across spaces of different size. Use progress, perp, rel_delta and the size baseline.
- Using the size baseline or the readout cosine to change a reading.
- The forward test alone as rejuvenation: an APPROACH sentence must carry the reverse numbers.
- Any MD or AGE number if the FULL rerun misses the stop rule.
- An invalid CI passing a gate.
- Using half-B cells anywhere except the mixture control; using half-A cells in it; changing f values.
- Per-timepoint slicing. Averaging donors, directions, states or spaces.
- Frozen-ruler movement. Context only.
- Any number from paper/ or paper_package/.

AMENDMENT 1 — written 2026-09-30, **before any GENESPACE statistic**. Adds to sections 6, 2/7 and 8. Nothing above is changed.

A1. SIZE BASELINE MATCHED ON MEASUREMENT (adds to section 6)
Reason (owner): the age lists hold many pseudogenes and non-coding genes the reprogramming data do not measure (results/md3/t1_ams_meta_*.json: ~896 of 1,533 age-up, ~1,586 of 2,007 age-down found there; cited as context, not an input). Random genes are measured more often, so unmatched random sets would hold more genes that can move.
- Zero status: a ruler gene is "zero" iff its summed counts are 0 in every row of the full-gene main SAME panel [O, Y, S, S_rev, Pluri, NonReprog] (section 2); otherwise "nonzero". This is the same mask same_run.py sets to z=0. The script fixes it from the panel counts before computing any statistic and saves it per gene in gene_map.csv.
- Strata = 24 GTEx-mu bins (section 6) × {zero, nonzero} = 48 strata. In each draw, for each stratum in order (bin 0..23; zero before nonzero within a bin), draw (n space genes in that stratum) genes without replacement from all ruler genes in that stratum. Space genes are not excluded. Everything else in section 6 is unchanged (200 draws, one Generator(20260914) per space, MD before AGE, point estimates only, reporting).
- Report per space (FULL, MD, AGE) and per list (age_up, age_down): n mapped to the ruler, n nonzero in the main SAME panel, n zero, and the per-stratum counts.

A2. MIXTURE NOISE FLOOR (adds to sections 2 and 7)
In every space, the f = 0 mixture's progress, delta, perp and cos are reported as the noise floor, as in the paper. They are not subtracted from any mixture, forward, reverse or context statistic. The section-7 gate and reading use uncorrected values.

A3. FULL BOOTSTRAP-CI COMPARISON (adds to section 8; reported only, NOT a stop rule)
After the section-8 point check, compare the FULL rerun's 95% CIs with the paper's:
  results/same/summary.json : frac_S_ci_lo/hi, delta_S_ci_lo/hi, frac_rev_ci_lo/hi, and mixture progress CI at f=0.50 (p.frac50_ci_lo/hi).
  The mixture progress CIs at f = 0, 0.10, 0.25 are not in summary.json; they are compared with results/same/posctrl.csv (frac_ci_lo/hi by f).
Write each rerun bound, reference bound and difference to repro_check.json and FINDINGS_GENESPACE.md. A difference never stops the run, never replaces the rerun's CIs, and never changes a gate or reading.
~~~~

## Results

### Run integrity

- No STOP.
- SHA-256 of 82 protected files before and after: unchanged = True.

### Gene spaces (section 1, Amendment 1)

Counts recomputed by the script match PREREG.flag (otherwise the run would have stopped). Zero = summed counts 0 in every row of the full-gene main SAME panel [O, Y, S, S_rev, Pluri, NonReprog].

| set | n_mapped_to_ruler | n_nonzero_main_SAME | n_zero_main_SAME |
|---|---|---|---|
| space FULL | 23485 | 20353 | 3132 |
| space MD | 205 | 205 | 0 |
| space AGE | 3089 | 2756 | 333 |
| list age_up | 1312 | 1127 | 185 |
| list age_down | 1777 | 1629 | 148 |
| list MD | 205 | 205 | 0 |

Genes at z=0 per space and panel (zero counts in every row of that panel):

| space | main_same | mixture | gtex_GM00731 | gtex_GM23815 |
|---|---|---|---|---|
| FULL | 3132 | 3507 | 3139 | 3138 |
| MD | 0 | 0 | 0 | 0 |
| AGE | 333 | 342 | 332 | 333 |

Per-stratum counts (24 GTEx-mu bins × zero/nonzero) are in `results/genespace/spaces.json`.

### Reproduction check, FULL space (section 8; stop rule |new − reference| ≤ 1e-6)

Passed: **True**.

| quantity | rerun | reference | diff | pass |
|---|---|---|---|---|
| same/summary.json frac_S | 0.674440748 | 0.674440748 | 0.00e+00 | True |
| same/summary.json delta_S | 130.893565 | 130.893565 | 0.00e+00 | True |
| same/summary.json frac_rev | -0.0572752268 | -0.0572752268 | 0.00e+00 | True |
| same/summary.json p.fracs[f=0.00] | 0.114587446 | 0.114587446 | 0.00e+00 | True |
| same/summary.json p.fracs[f=0.10] | 0.182435371 | 0.182435371 | 0.00e+00 | True |
| same/summary.json p.fracs[f=0.25] | 0.282928627 | 0.282928627 | 0.00e+00 | True |
| same/summary.json p.fracs[f=0.50] | 0.440818861 | 0.440818861 | 0.00e+00 | True |
| toward/stats.csv GM00731 PartialReprog cos_S | 0.0772462933 | 0.0772462933 | 2.78e-17 | True |
| toward/stats.csv GM00731 PartialReprog delta_young | 9.36967585 | 9.36967585 | 0.00e+00 | True |

### FULL bootstrap CI comparison (Amendment 1, A3; reported only, not a stop rule)

| bound | rerun | reference | diff |
|---|---|---|---|
| frac_S_ci_lo | 0.667391996 | 0.667391996 | 0.00e+00 |
| frac_S_ci_hi | 0.680092116 | 0.680092116 | 0.00e+00 |
| delta_S_ci_lo | 130.473493 | 130.473493 | 0.00e+00 |
| delta_S_ci_hi | 133.918108 | 133.918108 | 0.00e+00 |
| frac_rev_ci_lo | -0.0628257203 | -0.0628257203 | 0.00e+00 |
| frac_rev_ci_hi | -0.0538346533 | -0.0538346533 | 0.00e+00 |
| p.frac50_ci_lo | 0.424343381 | 0.424343381 | 0.00e+00 |
| p.frac50_ci_hi | 0.447458641 | 0.447458641 | 0.00e+00 |
| posctrl.csv f=0.00 frac_ci_lo | 0.0987166335 | 0.0987166335 | 6.94e-17 |
| posctrl.csv f=0.00 frac_ci_hi | 0.120951827 | 0.120951827 | 6.94e-17 |
| posctrl.csv f=0.10 frac_ci_lo | 0.167314998 | 0.167314998 | 5.55e-17 |
| posctrl.csv f=0.10 frac_ci_hi | 0.19069802 | 0.19069802 | 5.55e-17 |
| posctrl.csv f=0.25 frac_ci_lo | 0.267134363 | 0.267134363 | 5.55e-17 |
| posctrl.csv f=0.25 frac_ci_hi | 0.29139174 | 0.29139174 | 5.55e-17 |

### Same-platform test per space (diagnostic except where used by the section-7 reading)

Forward: v = z(Y) − z(O), S = GM00731 PartialReprog. Reverse: v_rev = z(O) − z(Y), S_rev = GM23815 PartialReprog. Pluri / NonReprog: GM00731, forward frame (context). N1 p for forward and reverse. Geometry labels are context except MD/AGE forward, which the section-7 reading uses only after the gate. Raw delta is not compared across spaces.

| space | n_genes | test | n_cells | cos | progress | delta | perp | rel_delta | geometry_label | p_N1 | asymmetry |
|---|---|---|---|---|---|---|---|---|---|---|---|
| FULL | 23485 | forward | 5832 | +0.323 | +0.674 | +130.894 | +1.975 | +1.002 | SIDEWAYS | +0.005 | +0.732 |
| FULL | 23485 | reverse | 6481 | -0.029 | -0.057 | +163.342 | +1.986 | +1.250 | SIDEWAYS | +1.000 | NA |
| FULL | 23485 | Pluri | 439 | +0.338 | +1.403 | +382.106 | +3.903 | +2.924 | SIDEWAYS | NA | NA |
| FULL | 23485 | NonReprog | 5459 | +0.370 | +0.555 | +60.331 | +1.392 | +0.462 | SIDEWAYS | NA | NA |
| MD | 205 | forward | 5832 | +0.343 | +0.790 | +17.088 | +2.160 | +1.171 | SIDEWAYS | +0.005 | +1.020 |
| MD | 205 | reverse | 6481 | -0.095 | -0.229 | +24.881 | +2.409 | +1.704 | SIDEWAYS | +0.289 | NA |
| MD | 205 | Pluri | 439 | +0.323 | +1.327 | +42.288 | +3.883 | +2.897 | SIDEWAYS | NA | NA |
| MD | 205 | NonReprog | 5459 | +0.420 | +0.633 | +6.086 | +1.369 | +0.417 | SIDEWAYS | NA | NA |
| AGE | 3089 | forward | 5832 | +0.395 | +0.822 | +44.605 | +1.915 | +0.923 | SIDEWAYS | +0.005 | +1.052 |
| AGE | 3089 | reverse | 6481 | -0.117 | -0.230 | +63.042 | +1.948 | +1.304 | SIDEWAYS | +1.000 | NA |
| AGE | 3089 | Pluri | 439 | +0.342 | +1.312 | +126.384 | +3.601 | +2.615 | SIDEWAYS | NA | NA |
| AGE | 3089 | NonReprog | 5459 | +0.391 | +0.602 | +22.802 | +1.417 | +0.472 | SIDEWAYS | NA | NA |

95% percentile bootstrap CIs (B=200; O and Y frozen; valid iff lo ≤ point ≤ hi):

| space | test | progress_CI | delta_CI | perp_CI | cos_CI | rel_delta_CI | asymmetry_CI |
|---|---|---|---|---|---|---|---|
| FULL | forward | [+0.667, +0.680] valid | [+130.473, +133.918] valid | [+1.972, +1.998] valid | [+0.318, +0.325] valid | [+0.998, +1.025] valid | [+0.724, +0.740] valid |
| FULL | reverse | [-0.063, -0.054] valid | [+163.312, +166.330] valid | [+1.984, +2.009] valid | [-0.032, -0.027] valid | [+1.249, +1.272] valid |  |
| FULL | Pluri | [+1.365, +1.442] valid | [+391.791, +418.643] invalid | [+3.979, +4.178] invalid | [+0.320, +0.331] invalid | [+2.996, +3.201] invalid |  |
| FULL | NonReprog | [+0.543, +0.561] valid | [+62.804, +64.575] invalid | [+1.411, +1.426] invalid | [+0.358, +0.368] invalid | [+0.480, +0.494] invalid |  |
| MD | forward | [+0.780, +0.811] valid | [+16.733, +17.419] valid | [+2.135, +2.184] valid | [+0.339, +0.352] valid | [+1.146, +1.192] valid | [+1.009, +1.055] valid |
| MD | reverse | [-0.251, -0.224] valid | [+24.598, +25.295] valid | [+2.379, +2.437] valid | [-0.103, -0.092] valid | [+1.682, +1.732] valid |  |
| MD | Pluri | [+1.251, +1.431] valid | [+40.188, +46.044] valid | [+3.734, +4.134] valid | [+0.303, +0.343] valid | [+2.747, +3.152] valid |  |
| MD | NonReprog | [+0.621, +0.648] valid | [+5.926, +6.241] valid | [+1.358, +1.379] valid | [+0.412, +0.428] valid | [+0.405, +0.427] valid |  |
| AGE | forward | [+0.812, +0.830] valid | [+44.293, +45.573] valid | [+1.908, +1.934] valid | [+0.389, +0.397] valid | [+0.916, +0.942] valid | [+1.040, +1.065] valid |
| AGE | reverse | [-0.240, -0.222] valid | [+62.812, +63.948] valid | [+1.941, +1.968] valid | [-0.122, -0.112] valid | [+1.299, +1.322] valid |  |
| AGE | Pluri | [+1.264, +1.374] valid | [+129.051, +140.346] invalid | [+3.654, +3.883] invalid | [+0.317, +0.342] invalid | [+2.667, +2.900] invalid |  |
| AGE | NonReprog | [+0.588, +0.613] valid | [+23.232, +24.173] invalid | [+1.425, +1.444] invalid | [+0.378, +0.394] valid | [+0.480, +0.500] invalid |  |

### Mixture positive control per space (f = 0 row = noise floor, Amendment 1 A2; not subtracted)

| space | f | n_aged | n_young | progress | delta | perp | cos | progress_CI | delta_CI | geometry_label |
|---|---|---|---|---|---|---|---|---|---|---|
| FULL | +0.000 | 2470 | 0 | +0.115 | +0.068 | +0.466 | +0.239 | [+0.099, +0.121] valid | [+5.807, +8.334] invalid | SIDEWAYS |
| FULL | +0.100 | 2223 | 247 | +0.182 | -7.870 | +0.457 | +0.371 | [+0.167, +0.191] valid | [-1.900, +1.007] invalid | APPROACH |
| FULL | +0.250 | 1852 | 618 | +0.283 | -18.913 | +0.452 | +0.530 | [+0.267, +0.291] valid | [-11.793, -9.372] invalid | APPROACH |
| FULL | +0.500 | 1235 | 1235 | +0.441 | -34.443 | +0.458 | +0.694 | [+0.424, +0.447] valid | [-25.403, -22.414] invalid | APPROACH |
| MD | +0.000 | 2470 | 0 | +0.010 | -0.007 | +0.140 | +0.074 | [-0.011, +0.028] valid | [-0.197, +0.376] valid | APPROACH |
| MD | +0.100 | 2223 | 247 | +0.096 | -1.195 | +0.152 | +0.531 | [+0.073, +0.116] valid | [-1.361, -0.797] valid | APPROACH |
| MD | +0.250 | 1852 | 618 | +0.207 | -2.676 | +0.188 | +0.741 | [+0.183, +0.227] valid | [-2.848, -2.134] valid | APPROACH |
| MD | +0.500 | 1235 | 1235 | +0.388 | -5.143 | +0.200 | +0.889 | [+0.367, +0.407] valid | [-5.286, -4.708] valid | APPROACH |
| AGE | +0.000 | 2470 | 0 | +0.058 | +0.431 | +0.361 | +0.158 | [+0.037, +0.074] valid | [+1.384, +2.932] invalid | SIDEWAYS |
| AGE | +0.100 | 2223 | 247 | +0.148 | -3.473 | +0.361 | +0.379 | [+0.129, +0.165] valid | [-2.395, -0.785] invalid | APPROACH |
| AGE | +0.250 | 1852 | 618 | +0.264 | -8.174 | +0.372 | +0.579 | [+0.245, +0.279] valid | [-6.942, -5.253] invalid | APPROACH |
| AGE | +0.500 | 1235 | 1235 | +0.437 | -14.882 | +0.381 | +0.753 | [+0.414, +0.452] valid | [-13.180, -11.300] invalid | APPROACH |

Section-7 gate per space:

| space | gate_pass | non_decreasing | ci50_excludes_0 | ci50_valid | delta50_lt_0 | use |
|---|---|---|---|---|---|---|
| FULL | True | True | True | True | True | reproduction check only |
| MD | True | True | True | True | True | section-7 gate |
| AGE | True | True | True | True | True | section-7 gate |

### GTEx-anchored test per space (secondary; diagnostic)

Anchors: YOUNG (20-29 ∪ 30-39) n=113, OLD (60-69 ∪ 70-79) n=231, middle (40-49, 50-59) n=308 (reported, not in anchors), all donors n=652.

C3 five-fold positive control (pass = all folds finite and cos_fold > 0.50):

| space | fold | n_young_train | n_old_train | n_old_test | cos | cos_CI | pass_bar | construction_ok |
|---|---|---|---|---|---|---|---|---|
| FULL | 0 | 94 | 187 | 44 | +0.615 | [+0.123, +0.754] valid | True | True |
| FULL | 1 | 92 | 181 | 50 | +0.849 | [+0.701, +0.835] invalid | True | True |
| FULL | 2 | 86 | 180 | 51 | +0.714 | [+0.300, +0.821] valid | True | True |
| FULL | 3 | 87 | 191 | 40 | +0.825 | [+0.457, +0.846] valid | True | True |
| FULL | 4 | 93 | 185 | 46 | +0.874 | [+0.677, +0.844] invalid | True | True |
| MD | 0 | 94 | 187 | 44 | +0.733 | [+0.308, +0.844] valid | True | True |
| MD | 1 | 92 | 181 | 50 | +0.864 | [+0.716, +0.868] valid | True | True |
| MD | 2 | 86 | 180 | 51 | +0.809 | [+0.313, +0.908] valid | True | True |
| MD | 3 | 87 | 191 | 40 | +0.828 | [+0.434, +0.870] valid | True | True |
| MD | 4 | 93 | 185 | 46 | +0.895 | [+0.658, +0.877] invalid | True | True |
| AGE | 0 | 94 | 187 | 44 | +0.681 | [+0.117, +0.799] valid | True | True |
| AGE | 1 | 92 | 181 | 50 | +0.866 | [+0.754, +0.867] valid | True | True |
| AGE | 2 | 86 | 180 | 51 | +0.721 | [+0.181, +0.877] valid | True | True |
| AGE | 3 | 87 | 191 | 40 | +0.842 | [+0.328, +0.891] valid | True | True |
| AGE | 4 | 93 | 185 | 46 | +0.895 | [+0.664, +0.871] invalid | True | True |

C3 pass per space: FULL=True, MD=True, AGE=True.

State statistics (d = z(S) − z(d0 Fibroblast); v = w(X) = c_young(X) − c_old(X); delta_young/delta_old computed directly):

| space | cell_line | state | n_cells | cos_S | progress | perp | delta_young | rel_delta_young | delta_old | rel_delta_old | p_N1 | p_N2 | meets_toward |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| FULL | GM00731 | PartialReprog | 5832 | +0.077 | +0.650 | +8.390 | +9.370 | +0.020 | +10.857 | +0.023 | +0.005 | +0.388 | False |
| FULL | GM00731 | EarlyPluripotency | 4414 | +0.117 | +0.977 | +8.262 | +30.456 | +0.065 | +32.695 | +0.070 | +0.005 | +0.080 | False |
| FULL | GM00731 | Pluripotency | 439 | +0.126 | +2.128 | +16.791 | +184.049 | +0.394 | +188.313 | +0.405 | +0.005 | +0.139 | False |
| FULL | GM00731 | NonReprog | 5459 | +0.063 | +0.377 | +5.959 | +8.932 | +0.019 | +9.818 | +0.021 | +0.005 | +0.308 | False |
| FULL | GM23815 | PartialReprog | 6481 | +0.052 | +0.422 | +8.057 | +14.189 | +0.031 | +15.210 | +0.033 | +0.005 | +0.463 | False |
| FULL | GM23815 | EarlyPluripotency | 12961 | +0.076 | +0.639 | +8.406 | +46.322 | +0.101 | +47.903 | +0.105 | +0.005 | +0.214 | False |
| FULL | GM23815 | Pluripotency | 1249 | +0.112 | +1.789 | +15.855 | +185.957 | +0.407 | +189.656 | +0.417 | +0.005 | +0.149 | False |
| FULL | GM23815 | NonReprog | 1812 | +0.010 | +0.068 | +6.808 | +41.971 | +0.092 | +42.345 | +0.093 | +0.040 | +0.448 | False |
| MD | GM00731 | PartialReprog | 5832 | +0.209 | +2.162 | +10.136 | +6.461 | +0.154 | +7.002 | +0.168 | +0.080 | +0.279 | False |
| MD | GM00731 | EarlyPluripotency | 4414 | +0.225 | +2.344 | +10.143 | +7.447 | +0.177 | +8.025 | +0.193 | +0.050 | +0.199 | False |
| MD | GM00731 | Pluripotency | 439 | +0.300 | +5.529 | +17.611 | +28.216 | +0.672 | +29.236 | +0.703 | +0.005 | +0.109 | False |
| MD | GM00731 | NonReprog | 5459 | +0.053 | +0.351 | +6.646 | -1.400 | -0.033 | -1.320 | -0.032 | +0.259 | +0.413 | False |
| MD | GM23815 | PartialReprog | 6481 | +0.210 | +2.278 | +10.610 | +11.798 | +0.292 | +12.322 | +0.307 | +0.100 | +0.274 | False |
| MD | GM23815 | EarlyPluripotency | 12961 | +0.201 | +2.187 | +10.646 | +11.686 | +0.289 | +12.191 | +0.303 | +0.095 | +0.184 | False |
| MD | GM23815 | Pluripotency | 1249 | +0.273 | +5.137 | +18.066 | +33.473 | +0.828 | +34.326 | +0.854 | +0.020 | +0.134 | False |
| MD | GM23815 | NonReprog | 1812 | +0.017 | +0.116 | +6.718 | +5.685 | +0.141 | +5.741 | +0.143 | +0.532 | +0.453 | False |
| AGE | GM00731 | PartialReprog | 5832 | -0.060 | -0.493 | +8.135 | -5.439 | -0.033 | -5.925 | -0.036 | +0.995 | +0.632 | False |
| AGE | GM00731 | EarlyPluripotency | 4414 | +0.057 | +0.441 | +7.744 | +1.686 | +0.010 | +2.093 | +0.013 | +0.005 | +0.338 | False |
| AGE | GM00731 | Pluripotency | 439 | +0.102 | +1.552 | +15.154 | +44.720 | +0.273 | +45.931 | +0.281 | +0.005 | +0.234 | False |
| AGE | GM00731 | NonReprog | 5459 | +0.049 | +0.297 | +5.999 | -1.015 | -0.006 | -0.742 | -0.005 | +0.020 | +0.303 | False |
| AGE | GM23815 | PartialReprog | 6481 | -0.065 | -0.504 | +7.719 | -0.054 | -0.000 | -0.537 | -0.003 | +1.000 | +0.637 | False |
| AGE | GM23815 | EarlyPluripotency | 12961 | +0.014 | +0.111 | +7.760 | +10.716 | +0.068 | +10.852 | +0.069 | +0.209 | +0.547 | False |
| AGE | GM23815 | Pluripotency | 1249 | +0.073 | +1.032 | +14.049 | +46.968 | +0.296 | +47.860 | +0.303 | +0.025 | +0.279 | False |
| AGE | GM23815 | NonReprog | 1812 | +0.049 | +0.322 | +6.545 | +15.948 | +0.100 | +16.281 | +0.103 | +0.040 | +0.199 | False |

Bootstrap CIs (Generator(20260918 + donor offset); origin and destinations resampled):

| space | cell_line | state | cos_S_CI | progress_CI | delta_young_CI | delta_old_CI |
|---|---|---|---|---|---|---|
| FULL | GM00731 | PartialReprog | [+0.073, +0.079] valid | [+0.623, +0.677] valid | [+7.722, +10.125] valid | [+9.187, +11.620] valid |
| FULL | GM00731 | EarlyPluripotency | [+0.113, +0.119] valid | [+0.953, +1.005] valid | [+28.992, +31.568] valid | [+31.211, +33.778] valid |
| FULL | GM00731 | Pluripotency | [+0.113, +0.120] invalid | [+1.982, +2.142] valid | [+192.650, +212.674] invalid | [+196.728, +216.867] invalid |
| FULL | GM00731 | NonReprog | [+0.056, +0.065] valid | [+0.348, +0.404] valid | [+8.058, +10.075] valid | [+8.930, +10.962] valid |
| FULL | GM23815 | PartialReprog | [+0.049, +0.055] valid | [+0.399, +0.449] valid | [+12.780, +14.675] valid | [+13.772, +15.731] valid |
| FULL | GM23815 | EarlyPluripotency | [+0.073, +0.078] valid | [+0.625, +0.663] valid | [+44.739, +46.664] valid | [+46.289, +48.236] valid |
| FULL | GM23815 | Pluripotency | [+0.102, +0.110] invalid | [+1.672, +1.810] valid | [+189.733, +206.295] invalid | [+193.234, +209.935] invalid |
| FULL | GM23815 | NonReprog | [-0.003, +0.009] invalid | [-0.021, +0.072] valid | [+48.519, +52.503] invalid | [+48.858, +52.722] invalid |
| MD | GM00731 | PartialReprog | [+0.201, +0.214] valid | [+2.068, +2.216] valid | [+6.045, +6.691] valid | [+6.579, +7.241] valid |
| MD | GM00731 | EarlyPluripotency | [+0.217, +0.231] valid | [+2.246, +2.413] valid | [+6.957, +7.801] valid | [+7.522, +8.384] valid |
| MD | GM00731 | Pluripotency | [+0.285, +0.311] valid | [+5.188, +5.971] valid | [+26.174, +31.647] valid | [+27.183, +32.667] valid |
| MD | GM00731 | NonReprog | [+0.041, +0.063] valid | [+0.271, +0.425] valid | [-1.650, -1.128] valid | [-1.583, -1.041] valid |
| MD | GM23815 | PartialReprog | [+0.204, +0.215] valid | [+2.194, +2.341] valid | [+11.242, +12.019] valid | [+11.752, +12.541] valid |
| MD | GM23815 | EarlyPluripotency | [+0.196, +0.206] valid | [+2.122, +2.249] valid | [+11.259, +11.912] valid | [+11.754, +12.423] valid |
| MD | GM23815 | Pluripotency | [+0.262, +0.284] valid | [+4.883, +5.467] valid | [+31.517, +36.216] valid | [+32.370, +37.041] valid |
| MD | GM23815 | NonReprog | [-0.007, +0.054] valid | [-0.047, +0.401] valid | [+5.174, +7.540] valid | [+5.218, +7.591] valid |
| AGE | GM00731 | PartialReprog | [-0.067, -0.055] valid | [-0.554, -0.453] valid | [-6.281, -4.798] valid | [-6.769, -5.265] valid |
| AGE | GM00731 | EarlyPluripotency | [+0.048, +0.063] valid | [+0.376, +0.495] valid | [+0.957, +2.367] valid | [+1.383, +2.766] valid |
| AGE | GM00731 | Pluripotency | [+0.090, +0.108] valid | [+1.427, +1.698] valid | [+46.237, +56.614] invalid | [+47.494, +57.798] invalid |
| AGE | GM00731 | NonReprog | [+0.037, +0.058] valid | [+0.225, +0.353] valid | [-1.709, -0.092] valid | [-1.483, +0.192] valid |
| AGE | GM23815 | PartialReprog | [-0.071, -0.059] valid | [-0.555, -0.463] valid | [-0.991, +0.526] valid | [-1.463, +0.051] valid |
| AGE | GM23815 | EarlyPluripotency | [+0.008, +0.019] valid | [+0.064, +0.147] valid | [+9.702, +11.138] valid | [+9.819, +11.303] valid |
| AGE | GM23815 | Pluripotency | [+0.062, +0.078] valid | [+0.906, +1.147] valid | [+48.508, +54.844] invalid | [+49.327, +55.770] invalid |
| AGE | GM23815 | NonReprog | [+0.033, +0.061] valid | [+0.233, +0.429] valid | [+16.787, +20.185] invalid | [+17.204, +20.545] invalid |

Secondary key per donor (results/toward/PREREG.flag decision order in X; GM23815 is a contrast; never changes the section-7 reading):

- FULL GM00731: `away_from_old_only` — GM00731: PartialReprog has cos_S=+0.077 but delta_young=+9.370 (n_cells=5832). Movement is away-from-old only. This is a real negative finding: a 1D age score cannot distinguish toward-young from away-from-old. The frozen-ruler score is not used to rescue this.
- FULL GM23815: `away_from_old_only` — CONTRAST (young donor GM23815; never the reading): PartialReprog has cos_S=+0.052 but delta_young=+14.189 (n_cells=6481). Movement is away-from-old only. This is a real negative finding: a 1D age score cannot distinguish toward-young from away-from-old. The frozen-ruler score is not used to rescue this.
- MD GM00731: `away_from_old_only` — GM00731: PartialReprog has cos_S=+0.209 but delta_young=+6.461 (n_cells=5832). Movement is away-from-old only. This is a real negative finding: a 1D age score cannot distinguish toward-young from away-from-old. The frozen-ruler score is not used to rescue this.
- MD GM23815: `away_from_old_only` — CONTRAST (young donor GM23815; never the reading): PartialReprog has cos_S=+0.210 but delta_young=+11.798 (n_cells=6481). Movement is away-from-old only. This is a real negative finding: a 1D age score cannot distinguish toward-young from away-from-old. The frozen-ruler score is not used to rescue this.
- AGE GM00731: `no_directional_structure` — GM00731: PartialReprog (n_cells=5832, cos_S=-0.060, p_N1=0.9950, p_N2=0.6318, delta_young=-5.439): Both nulls fail (p>0.05). No directional structure detectable at this n. The frozen-ruler score is not used to rescue this.
- AGE GM23815: `no_directional_structure` — CONTRAST (young donor GM23815; never the reading): PartialReprog (n_cells=6481, cos_S=-0.065, p_N1=1.0000, p_N2=0.6368, delta_young=-0.054): Both nulls fail (p>0.05). No directional structure detectable at this n. The frozen-ruler score is not used to rescue this.

### Readout-vector cosine (section 5; context only)

Unweighted sign vectors (r_MD: −1 on MD genes; r_AGE: −1 age_up, +1 age_down) dotted with frozen-GTEx z of TMM pseudobulk displacements. This approximates the published AddModuleScore and does not reproduce it (AddModuleScore is per cell on log-normalized counts, subtracts expression-bin-matched control genes, and is summarized per cluster).

| space | readout | n_nonzero | test | vector | cell_line | state | cos |
|---|---|---|---|---|---|---|---|
| FULL | r_MD | 205 | same | forward | GM00731 | PartialReprog | +0.066 |
| FULL | r_MD | 205 | same | reverse | GM23815 | PartialReprog | +0.083 |
| FULL | r_MD | 205 | same | Pluri | GM00731 | Pluripotency | +0.073 |
| FULL | r_MD | 205 | same | NonReprog | GM00731 | NonReprog | -0.001 |
| FULL | r_MD | 205 | same | mixture |  | f=0.00 | +0.000 |
| FULL | r_MD | 205 | same | mixture |  | f=0.10 | -0.003 |
| FULL | r_MD | 205 | same | mixture |  | f=0.25 | -0.001 |
| FULL | r_MD | 205 | same | mixture |  | f=0.50 | +0.001 |
| FULL | r_MD | 205 | same | axis |  | v = z(Y) - z(O) | +0.023 |
| FULL | r_MD | 205 | gtex | axis |  | u(X) | +0.022 |
| FULL | r_MD | 205 | gtex | d0_to_state | GM00731 | PartialReprog | +0.069 |
| FULL | r_MD | 205 | gtex | d0_to_state | GM00731 | EarlyPluripotency | +0.064 |
| FULL | r_MD | 205 | gtex | d0_to_state | GM00731 | Pluripotency | +0.074 |
| FULL | r_MD | 205 | gtex | d0_to_state | GM00731 | NonReprog | +0.001 |
| FULL | r_MD | 205 | gtex | d0_to_state | GM23815 | PartialReprog | +0.085 |
| FULL | r_MD | 205 | gtex | d0_to_state | GM23815 | EarlyPluripotency | +0.070 |
| FULL | r_MD | 205 | gtex | d0_to_state | GM23815 | Pluripotency | +0.083 |
| FULL | r_MD | 205 | gtex | d0_to_state | GM23815 | NonReprog | +0.010 |
| FULL | r_AGE | 3089 | same | forward | GM00731 | PartialReprog | +0.147 |
| FULL | r_AGE | 3089 | same | reverse | GM23815 | PartialReprog | +0.136 |
| FULL | r_AGE | 3089 | same | Pluri | GM00731 | Pluripotency | -0.037 |
| FULL | r_AGE | 3089 | same | NonReprog | GM00731 | NonReprog | +0.025 |
| FULL | r_AGE | 3089 | same | mixture |  | f=0.00 | +0.003 |
| FULL | r_AGE | 3089 | same | mixture |  | f=0.10 | +0.017 |
| FULL | r_AGE | 3089 | same | mixture |  | f=0.25 | +0.034 |
| FULL | r_AGE | 3089 | same | mixture |  | f=0.50 | +0.049 |
| FULL | r_AGE | 3089 | same | axis |  | v = z(Y) - z(O) | +0.067 |
| FULL | r_AGE | 3089 | gtex | axis |  | u(X) | -0.142 |
| FULL | r_AGE | 3089 | gtex | d0_to_state | GM00731 | PartialReprog | +0.144 |
| FULL | r_AGE | 3089 | gtex | d0_to_state | GM00731 | EarlyPluripotency | +0.048 |
| FULL | r_AGE | 3089 | gtex | d0_to_state | GM00731 | Pluripotency | -0.039 |
| FULL | r_AGE | 3089 | gtex | d0_to_state | GM00731 | NonReprog | +0.023 |
| FULL | r_AGE | 3089 | gtex | d0_to_state | GM23815 | PartialReprog | +0.135 |
| FULL | r_AGE | 3089 | gtex | d0_to_state | GM23815 | EarlyPluripotency | +0.060 |
| FULL | r_AGE | 3089 | gtex | d0_to_state | GM23815 | Pluripotency | -0.024 |
| FULL | r_AGE | 3089 | gtex | d0_to_state | GM23815 | NonReprog | +0.014 |
| MD | r_MD | 205 | same | forward | GM00731 | PartialReprog | +0.534 |
| MD | r_MD | 205 | same | reverse | GM23815 | PartialReprog | +0.608 |
| MD | r_MD | 205 | same | Pluri | GM00731 | Pluripotency | +0.662 |
| MD | r_MD | 205 | same | NonReprog | GM00731 | NonReprog | -0.008 |
| MD | r_MD | 205 | same | mixture |  | f=0.00 | +0.006 |
| MD | r_MD | 205 | same | mixture |  | f=0.10 | -0.065 |
| MD | r_MD | 205 | same | mixture |  | f=0.25 | -0.015 |
| MD | r_MD | 205 | same | mixture |  | f=0.50 | +0.017 |
| MD | r_MD | 205 | same | axis |  | v = z(Y) - z(O) | +0.202 |
| MD | r_MD | 205 | gtex | axis |  | u(X) | +0.215 |
| MD | r_MD | 205 | gtex | d0_to_state | GM00731 | PartialReprog | +0.551 |
| MD | r_MD | 205 | gtex | d0_to_state | GM00731 | EarlyPluripotency | +0.500 |
| MD | r_MD | 205 | gtex | d0_to_state | GM00731 | Pluripotency | +0.669 |
| MD | r_MD | 205 | gtex | d0_to_state | GM00731 | NonReprog | +0.008 |
| MD | r_MD | 205 | gtex | d0_to_state | GM23815 | PartialReprog | +0.617 |
| MD | r_MD | 205 | gtex | d0_to_state | GM23815 | EarlyPluripotency | +0.529 |
| MD | r_MD | 205 | gtex | d0_to_state | GM23815 | Pluripotency | +0.689 |
| MD | r_MD | 205 | gtex | d0_to_state | GM23815 | NonReprog | +0.100 |
| MD | r_AGE | 43 | same | forward | GM00731 | PartialReprog | +0.160 |
| MD | r_AGE | 43 | same | reverse | GM23815 | PartialReprog | +0.103 |
| MD | r_AGE | 43 | same | Pluri | GM00731 | Pluripotency | -0.004 |
| MD | r_AGE | 43 | same | NonReprog | GM00731 | NonReprog | +0.115 |
| MD | r_AGE | 43 | same | mixture |  | f=0.00 | +0.093 |
| MD | r_AGE | 43 | same | mixture |  | f=0.10 | +0.107 |
| MD | r_AGE | 43 | same | mixture |  | f=0.25 | +0.125 |
| MD | r_AGE | 43 | same | mixture |  | f=0.50 | +0.104 |
| MD | r_AGE | 43 | same | axis |  | v = z(Y) - z(O) | +0.093 |
| MD | r_AGE | 43 | gtex | axis |  | u(X) | -0.213 |
| MD | r_AGE | 43 | gtex | d0_to_state | GM00731 | PartialReprog | +0.153 |
| MD | r_AGE | 43 | gtex | d0_to_state | GM00731 | EarlyPluripotency | +0.065 |
| MD | r_AGE | 43 | gtex | d0_to_state | GM00731 | Pluripotency | -0.007 |
| MD | r_AGE | 43 | gtex | d0_to_state | GM00731 | NonReprog | +0.110 |
| MD | r_AGE | 43 | gtex | d0_to_state | GM23815 | PartialReprog | +0.103 |
| MD | r_AGE | 43 | gtex | d0_to_state | GM23815 | EarlyPluripotency | +0.054 |
| MD | r_AGE | 43 | gtex | d0_to_state | GM23815 | Pluripotency | -0.027 |
| MD | r_AGE | 43 | gtex | d0_to_state | GM23815 | NonReprog | +0.145 |
| AGE | r_MD | 43 | same | forward | GM00731 | PartialReprog | +0.059 |
| AGE | r_MD | 43 | same | reverse | GM23815 | PartialReprog | +0.073 |
| AGE | r_MD | 43 | same | Pluri | GM00731 | Pluripotency | +0.084 |
| AGE | r_MD | 43 | same | NonReprog | GM00731 | NonReprog | -0.022 |
| AGE | r_MD | 43 | same | mixture |  | f=0.00 | -0.005 |
| AGE | r_MD | 43 | same | mixture |  | f=0.10 | -0.009 |
| AGE | r_MD | 43 | same | mixture |  | f=0.25 | -0.009 |
| AGE | r_MD | 43 | same | mixture |  | f=0.50 | -0.005 |
| AGE | r_MD | 43 | same | axis |  | v = z(Y) - z(O) | +0.015 |
| AGE | r_MD | 43 | gtex | axis |  | u(X) | +0.039 |
| AGE | r_MD | 43 | gtex | d0_to_state | GM00731 | PartialReprog | +0.063 |
| AGE | r_MD | 43 | gtex | d0_to_state | GM00731 | EarlyPluripotency | +0.068 |
| AGE | r_MD | 43 | gtex | d0_to_state | GM00731 | Pluripotency | +0.085 |
| AGE | r_MD | 43 | gtex | d0_to_state | GM00731 | NonReprog | -0.020 |
| AGE | r_MD | 43 | gtex | d0_to_state | GM23815 | PartialReprog | +0.074 |
| AGE | r_MD | 43 | gtex | d0_to_state | GM23815 | EarlyPluripotency | +0.062 |
| AGE | r_MD | 43 | gtex | d0_to_state | GM23815 | Pluripotency | +0.098 |
| AGE | r_MD | 43 | gtex | d0_to_state | GM23815 | NonReprog | -0.020 |
| AGE | r_AGE | 3089 | same | forward | GM00731 | PartialReprog | +0.397 |
| AGE | r_AGE | 3089 | same | reverse | GM23815 | PartialReprog | +0.372 |
| AGE | r_AGE | 3089 | same | Pluri | GM00731 | Pluripotency | -0.109 |
| AGE | r_AGE | 3089 | same | NonReprog | GM00731 | NonReprog | +0.066 |
| AGE | r_AGE | 3089 | same | mixture |  | f=0.00 | +0.011 |
| AGE | r_AGE | 3089 | same | mixture |  | f=0.10 | +0.057 |
| AGE | r_AGE | 3089 | same | mixture |  | f=0.25 | +0.108 |
| AGE | r_AGE | 3089 | same | mixture |  | f=0.50 | +0.143 |
| AGE | r_AGE | 3089 | same | axis |  | v = z(Y) - z(O) | +0.182 |
| AGE | r_AGE | 3089 | gtex | axis |  | u(X) | -0.373 |
| AGE | r_AGE | 3089 | gtex | d0_to_state | GM00731 | PartialReprog | +0.391 |
| AGE | r_AGE | 3089 | gtex | d0_to_state | GM00731 | EarlyPluripotency | +0.136 |
| AGE | r_AGE | 3089 | gtex | d0_to_state | GM00731 | Pluripotency | -0.113 |
| AGE | r_AGE | 3089 | gtex | d0_to_state | GM00731 | NonReprog | +0.061 |
| AGE | r_AGE | 3089 | gtex | d0_to_state | GM23815 | PartialReprog | +0.371 |
| AGE | r_AGE | 3089 | gtex | d0_to_state | GM23815 | EarlyPluripotency | +0.171 |
| AGE | r_AGE | 3089 | gtex | d0_to_state | GM23815 | Pluripotency | -0.073 |
| AGE | r_AGE | 3089 | gtex | d0_to_state | GM23815 | NonReprog | +0.038 |

### Size baselines (section 6 and Amendment 1 A1; context, never changes a reading)

`sec6_mu`: section 6 as written (24 equal-count GTEx-mu bins). `A1_mu_zero`: Amendment 1 (24 bins × zero/nonzero in the main SAME panel). Position = (#{random < real} + 0.5·#{random = real}) / 200. Selected rows; the full table is `results/genespace/size_baseline.csv`.

| baseline | space | test | cell_line | state | stat | real | position | rand_p2_5 | rand_p50 | rand_p97_5 |
|---|---|---|---|---|---|---|---|---|---|---|
| sec6_mu | MD | same |  | forward | cos | +0.343 | +0.190 | +0.270 | +0.402 | +0.510 |
| sec6_mu | MD | same |  | forward | progress | +0.790 | +0.140 | +0.674 | +1.025 | +1.414 |
| sec6_mu | MD | same |  | forward | perp | +2.160 | +0.190 | +1.908 | +2.317 | +2.724 |
| sec6_mu | MD | same |  | forward | rel_delta | +1.171 | +0.210 | +0.911 | +1.326 | +1.757 |
| sec6_mu | MD | same |  | reverse | cos | -0.095 | +0.930 | -0.333 | -0.214 | -0.057 |
| sec6_mu | MD | same |  | reverse | progress | -0.229 | +0.925 | -0.917 | -0.520 | -0.145 |
| sec6_mu | MD | same |  | reverse | perp | +2.409 | +0.420 | +2.007 | +2.441 | +2.882 |
| sec6_mu | MD | same |  | reverse | rel_delta | +1.704 | +0.225 | +1.398 | +1.893 | +2.343 |
| sec6_mu | MD | mixture |  | f=0.50 | progress | +0.388 | +0.365 | +0.332 | +0.399 | +0.468 |
| sec6_mu | MD | mixture |  | f=0.50 | delta | -5.143 | +0.010 | -4.979 | -3.831 | -2.874 |
| sec6_mu | MD | gtex | GM00731 | PartialReprog | cos_S | +0.209 | +0.845 | +0.020 | +0.148 | +0.254 |
| sec6_mu | MD | gtex | GM00731 | PartialReprog | progress | +2.162 | +0.930 | +0.187 | +1.320 | +2.379 |
| sec6_mu | MD | gtex | GM00731 | PartialReprog | perp | +10.136 | +0.985 | +7.928 | +8.911 | +9.962 |
| sec6_mu | MD | gtex | GM00731 | PartialReprog | rel_delta_young | +0.154 | +0.950 | -0.010 | +0.078 | +0.170 |
| sec6_mu | AGE | same |  | forward | cos | +0.395 | +1.000 | +0.298 | +0.328 | +0.361 |
| sec6_mu | AGE | same |  | forward | progress | +0.822 | +1.000 | +0.618 | +0.683 | +0.764 |
| sec6_mu | AGE | same |  | forward | perp | +1.915 | +0.185 | +1.867 | +1.972 | +2.077 |
| sec6_mu | AGE | same |  | forward | rel_delta | +0.923 | +0.075 | +0.898 | +0.996 | +1.104 |
| sec6_mu | AGE | same |  | reverse | cos | -0.117 | +0.000 | -0.085 | -0.047 | -0.019 |
| sec6_mu | AGE | same |  | reverse | progress | -0.230 | +0.000 | -0.170 | -0.096 | -0.037 |
| sec6_mu | AGE | same |  | reverse | perp | +1.948 | +0.225 | +1.897 | +1.993 | +2.108 |
| sec6_mu | AGE | same |  | reverse | rel_delta | +1.304 | +0.690 | +1.173 | +1.275 | +1.378 |
| sec6_mu | AGE | mixture |  | f=0.50 | progress | +0.437 | +0.345 | +0.421 | +0.441 | +0.463 |
| sec6_mu | AGE | mixture |  | f=0.50 | delta | -14.882 | +0.035 | -14.992 | -13.671 | -12.360 |
| sec6_mu | AGE | gtex | GM00731 | PartialReprog | cos_S | -0.060 | +0.000 | +0.037 | +0.067 | +0.098 |
| sec6_mu | AGE | gtex | GM00731 | PartialReprog | progress | -0.493 | +0.000 | +0.318 | +0.560 | +0.815 |
| sec6_mu | AGE | gtex | GM00731 | PartialReprog | perp | +8.135 | +0.115 | +7.997 | +8.338 | +8.668 |
| sec6_mu | AGE | gtex | GM00731 | PartialReprog | rel_delta_young | -0.033 | +0.000 | -0.022 | -0.001 | +0.018 |
| A1_mu_zero | MD | same |  | forward | cos | +0.343 | +0.200 | +0.259 | +0.393 | +0.514 |
| A1_mu_zero | MD | same |  | forward | progress | +0.790 | +0.170 | +0.655 | +0.980 | +1.387 |
| A1_mu_zero | MD | same |  | forward | perp | +2.160 | +0.255 | +1.891 | +2.291 | +2.697 |
| A1_mu_zero | MD | same |  | forward | rel_delta | +1.171 | +0.255 | +0.926 | +1.301 | +1.702 |
| A1_mu_zero | MD | same |  | reverse | cos | -0.095 | +0.900 | -0.336 | -0.192 | -0.046 |
| A1_mu_zero | MD | same |  | reverse | progress | -0.229 | +0.910 | -0.867 | -0.469 | -0.103 |
| A1_mu_zero | MD | same |  | reverse | perp | +2.409 | +0.465 | +1.995 | +2.425 | +2.829 |
| A1_mu_zero | MD | same |  | reverse | rel_delta | +1.704 | +0.275 | +1.367 | +1.843 | +2.259 |
| A1_mu_zero | MD | mixture |  | f=0.50 | progress | +0.388 | +0.350 | +0.331 | +0.404 | +0.482 |
| A1_mu_zero | MD | mixture |  | f=0.50 | delta | -5.143 | +0.015 | -5.030 | -3.924 | -2.915 |
| A1_mu_zero | MD | gtex | GM00731 | PartialReprog | cos_S | +0.209 | +0.835 | +0.017 | +0.146 | +0.261 |
| A1_mu_zero | MD | gtex | GM00731 | PartialReprog | progress | +2.162 | +0.900 | +0.152 | +1.312 | +2.381 |
| A1_mu_zero | MD | gtex | GM00731 | PartialReprog | perp | +10.136 | +0.945 | +7.975 | +8.976 | +10.424 |
| A1_mu_zero | MD | gtex | GM00731 | PartialReprog | rel_delta_young | +0.154 | +0.975 | -0.005 | +0.083 | +0.153 |
| A1_mu_zero | AGE | same |  | forward | cos | +0.395 | +1.000 | +0.293 | +0.329 | +0.364 |
| A1_mu_zero | AGE | same |  | forward | progress | +0.822 | +1.000 | +0.611 | +0.685 | +0.760 |
| A1_mu_zero | AGE | same |  | forward | perp | +1.915 | +0.140 | +1.883 | +1.966 | +2.066 |
| A1_mu_zero | AGE | same |  | forward | rel_delta | +0.923 | +0.060 | +0.913 | +0.991 | +1.087 |
| A1_mu_zero | AGE | same |  | reverse | cos | -0.117 | +0.000 | -0.084 | -0.048 | -0.013 |
| A1_mu_zero | AGE | same |  | reverse | progress | -0.230 | +0.000 | -0.165 | -0.096 | -0.026 |
| A1_mu_zero | AGE | same |  | reverse | perp | +1.948 | +0.195 | +1.902 | +1.995 | +2.069 |
| A1_mu_zero | AGE | same |  | reverse | rel_delta | +1.304 | +0.705 | +1.188 | +1.276 | +1.349 |
| A1_mu_zero | AGE | mixture |  | f=0.50 | progress | +0.437 | +0.325 | +0.420 | +0.440 | +0.464 |
| A1_mu_zero | AGE | mixture |  | f=0.50 | delta | -14.882 | +0.050 | -15.012 | -13.707 | -12.381 |
| A1_mu_zero | AGE | gtex | GM00731 | PartialReprog | cos_S | -0.060 | +0.000 | +0.038 | +0.067 | +0.100 |
| A1_mu_zero | AGE | gtex | GM00731 | PartialReprog | progress | -0.493 | +0.000 | +0.320 | +0.563 | +0.824 |
| A1_mu_zero | AGE | gtex | GM00731 | PartialReprog | perp | +8.135 | +0.110 | +8.048 | +8.316 | +8.564 |
| A1_mu_zero | AGE | gtex | GM00731 | PartialReprog | rel_delta_young | -0.033 | +0.010 | -0.027 | -0.002 | +0.021 |

Forward geometry label: real space vs counts among 200 random sets (no gate applied to random sets):

| baseline | space | real | label_counts | mean_overlap | mean_overlap_frac |
|---|---|---|---|---|---|
| sec6_mu | MD | SIDEWAYS | {"APPROACH": 0, "OVERSHOOT": 0, "SIDEWAYS": 200, "TIE": 0, "NA": 0} | +10.030 | +0.049 |
| sec6_mu | AGE | SIDEWAYS | {"APPROACH": 0, "OVERSHOOT": 0, "SIDEWAYS": 200, "TIE": 0, "NA": 0} | +447.325 | +0.145 |
| A1_mu_zero | MD | SIDEWAYS | {"APPROACH": 0, "OVERSHOOT": 0, "SIDEWAYS": 200, "TIE": 0, "NA": 0} | +10.485 | +0.051 |
| A1_mu_zero | AGE | SIDEWAYS | {"APPROACH": 0, "OVERSHOOT": 0, "SIDEWAYS": 200, "TIE": 0, "NA": 0} | +453.985 | +0.147 |

## Pre-registered reading (section 7; only what fired)

### MD: `SIDEWAYS`

In the MD gene space (205 ruler genes), aged partially-reprogrammed cells (GM00731 PartialReprog, n_cells=5832) end farther from the young donor's day-0 cells than they started (change in distance +17.088), mainly through movement off the aged-to-young axis: sideways term perp² = +4.667, overshoot term (progress − 1)² = +0.044, progress +0.790. This is not approach to the young donor.

delta CI [+16.733, +17.419] (valid); progress CI [+0.780, +0.811] (valid).

### AGE: `SIDEWAYS`

In the AGE gene space (3089 ruler genes), aged partially-reprogrammed cells (GM00731 PartialReprog, n_cells=5832) end farther from the young donor's day-0 cells than they started (change in distance +44.605), mainly through movement off the aged-to-young axis: sideways term perp² = +3.665, overshoot term (progress − 1)² = +0.032, progress +0.822. This is not approach to the young donor.

delta CI [+44.293, +45.573] (valid); progress CI [+0.812, +0.830] (valid).

MD and AGE are reported side by side; they are not averaged, ranked or combined. FULL is a reproduction check and gets no reading.

Context geometry labels (never a reading): FULL: reverse=SIDEWAYS, Pluri=SIDEWAYS, NonReprog=SIDEWAYS, forward=SIDEWAYS; MD: reverse=SIDEWAYS, Pluri=SIDEWAYS, NonReprog=SIDEWAYS; AGE: reverse=SIDEWAYS, Pluri=SIDEWAYS, NonReprog=SIDEWAYS.

## Files

- `src/genespace_run.py`
- `results/genespace/` (manifest.json, gene_map.csv, spaces.json, repro_check.json, same_stats.csv, same_posctrl.csv, gtex_stats.csv, gtex_c3.csv, readout_cos.csv, size_baseline.csv, size_baseline_sets.npz, boot.npz [bootstrap draws and null cosines], reading.json, summary.json, run_report.txt)
- `FINDINGS_GENESPACE.md`
- `PROGRESS_GENESPACE.md`


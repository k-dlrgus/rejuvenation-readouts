# FINDINGS_SOUTH7 — re-rank of the Southard CRISPRa screen with a qualified statistic

Status: Stage 0 passed. Stage A complete for the primary set; the all-1,836 secondary run (R2) is reported where available.
Stage B: STOP `fulldata_infeasible` at the B0 feasibility probe. Stage A is unaffected.

## Pre-registration (verbatim, `results/south7/PREREG.flag`)

```
TASK: SOUTH7 — re-rank the Southard CRISPRa screen with a qualified statistic,
a depth-matched floor, an anchor-free basis, and a permutation floor that does
not cap the answer.
Written before any SOUTH7 number exists.

Nature of this task: computational. No lab work, no new data beyond the Southard
files already on disk. Nothing in results/south6/ or results/south6rec/ is
modified; every file there is read as stored. The frozen ruler and the TOWARD
anchors are read, never refit. Existing src/ files are imported unchanged.

Two stages. Stage A is the result and must complete. Stage B is upside and may
STOP without affecting Stage A.

--- WHY EACH CHANGE IS MADE (fixed now, so no choice is made after seeing data) ---

1. SPACE. SOUTH6REC planted a known approach into real non-targeting pseudobulks
   at the screen's own depth and measured recovery. At 106 cells, S1 recovered a
   25%-of-gap approach with power 0.115 and a 100%-of-gap approach with power
   0.470. S2_k20, S2_k50 and S2_k100 recovered 25% with power 0.990, 0.965 and
   0.930. S1 at the screen's depth is blind; the compressed spaces are not.
   SOUTH6's primary was S1. That is the single reason its ranking carried no
   information.

2. BASIS. The SOUTH6 S2 basis was fit on all 652 GTEx donors, of which 344
   (52.8%) are the anchor donors averaged into O1 and Y1. SOUTH6REC R5 showed
   gap retention at k=20 falls only from 0.887 to 0.865 when the basis is refit
   on the 308 donors aged 40-59, which contains no anchor donor. The cost of
   removing the circularity is small, so it is removed.

3. FLOOR. SOUTH6's control floor (stage2b) was built from half-splits of the
   14,013 NT cells, about 7,006 cells per half, against factors with a median of
   106 cells. Displacement scales as 1/sqrt(cells), so that floor sits at roughly
   1/8 of the noise the factors carry. A floor at the wrong depth is not a floor.
   The same arithmetic explains the stored per-guide NT figure of 69 against the
   factors' 95: the 78 NT guides hold about 180 cells each, and 1273/sqrt(180)
   against 1273/sqrt(106) reproduces that difference with no biology in it.
   SOUTH7 builds the floor at each factor's own cell count.

4. PERMUTATIONS. SOUTH6 used 5,000 permutations, so the smallest attainable
   p is 2e-4. With Benjamini-Hochberg at m = 1,836 a factor can only pass if at
   least 8 factors are tied at that floor. S2_k20 on the primary pair had 4. That
   setting could not have reported a hit whatever the biology. SOUTH7 raises the
   permutation count so the floor cannot bind.

5. CRITERION. SOUTH6 called a factor approaching only if its change in distance
   was negative outright. The preprint's own remedy list says to report the
   zero-dose point as the noise floor rather than assuming it is zero. SOUTH7
   adopts the preprint's rule: a factor is judged against a depth-matched
   zero-effect null, not against zero. The SOUTH6 rule is reported alongside for
   every factor so the two can be compared directly.

--- DATA ---
The SOUTH6 file as analyzed, same md5 check as SOUTH6REC. Counts from SOUTH6's
matched-gene CSR cache. NT cells: guide_target == 'non' (14,013 cells, 78
guides). STOP `nt_cells_mismatch` if the obs-derived NT set differs from
results/south6/stage1.npz nt_cells.

--- SPACES ---
PRIMARY: S2A_k20. Projection onto the first 20 principal components of a basis
fit exactly as src/south6_run.py build_spaces does, except on the 308 GTEx
donors aged 40-59 only. No anchor donor enters the fit. Report the gap retention
of this basis.
SECONDARY, reported, never primary: S2A_k50, S2A_k100, and S1 (gene space).
Origin O1 = GTEx 60-79 centroid, target Y1 = GTEx 20-39 centroid, frozen, as
stored. The O2 pairs are not run. Their axis is the 468-unit origin-to-anchor
platform gap, not an age gap, and every SOUTH6 hit came from them.

--- FACTOR SET ---
PRIMARY: Southard's confirmed-effect factors — the set reported as showing
strong transcriptional effects in Hs27 cells, and the additional single-guide
set, taken from the Southard et al. supplementary tables. Record the table, the
column, the inclusion rule and the resulting count. Do not substitute a
self-derived list; if the published list cannot be sourced, STOP
`factor_list_unavailable`, report what was searched, and run nothing further.
SECONDARY: all 1,836 factors, so the effect of m on multiple testing is visible
rather than assumed.

--- STAGE 0: REPRODUCE FIRST (gate, before any new number) ---
Rebuild with this task's code and compare against stored files:
 (a) the SOUTH6 Stage 0 items that SOUTH6REC already checked (counts, transform,
     d_f, origin-target distances, per-factor delta/cos/frac/d_norm in S1 and
     each S2_k on O1_Y1);
 (b) the SOUTH6REC Stage 1 gate table: the 12 space-by-depth rows, their power
     at 25%, and the dose-0 median delta in each space at 106 cells.
Pass: counts exact, every numeric item within 1e-4 * max(1, |stored|).
Fail -> STOP `prior_results_not_reproduced`, naming every mismatch. No new
number is computed.

--- STAGE A: THE RE-RANK (gated, must complete) ---

A1. Depth-matched null. For each factor of n cells, draw R = 200 NT pseudobulks
of exactly n cells with the SOUTH6REC primary scheme: disjoint blocks within a
permutation, reference = the complement of the draw, never the pool containing
it. Record the 5th percentile of the null change in distance at that n, in each
space. For the primary factor set this is computed exactly, per factor. For the
all-1,836 secondary set, a grid of cell counts with interpolation of the 5th
percentile in 1/sqrt(n) is permitted; if used, validate it against at least 5
exactly computed grid points spanning the range and report the maximum
interpolation error. State which was done.

A2. Permutation null. Unchanged in construction from SOUTH6 — state in the
findings exactly what is permuted and confirm it is the same null — but with at
least 100,000 draws, so the attainable minimum p is at most 1e-5. Report the
attainable minimum p and the number of factors tied at it in every setting.

A3. Call. A factor is APPROACHING if, in the primary space:
      its change in distance is below the 5th percentile of its own
      depth-matched null, AND its Benjamini-Hochberg q over the factor set
      is <= 0.05.
Report for every factor, alongside: the SOUTH6 rule (change in distance < 0),
cell count, percentile rank of cell count, guide agreement, identity drop,
cosine, and the call in S2A_k50, S2A_k100 and S1.

A4. Flags, declared now. Any factor called approaching whose cell count sits
above the 90th percentile of the factor set is flagged `high_depth_candidate`
and reported as such, not as a result. SOUTH6's recurring 11 had a median of 183
cells against the screen's 106, with 8 of 11 above the 90th percentile; that set
was an artifact of how many cells each guide happened to capture. Any new
shortlist must be checked the same way before it is believed.

--- GATE WORDING, fixed now ---
PASS (k = 0): "<set>, S2A_k20: of <m> factors, none is called approaching. Each
factor is judged against a null built at its own cell count, and SOUTH6REC
showed this statistic recovers a 25%-of-gap approach at this depth with power
0.990. The screen contains no transcription factor that moves the reference
cells a quarter of the way toward young."
PASS (k > 0): "<set>, S2A_k20: <k> of <m> factors are called approaching:
<list>. <j> of them carry the high_depth_candidate flag." Each is reported with
its full row and its call in the three secondary spaces. No factor is described
as a candidate for the pilot panel in this findings file.

--- STAGE B: FULL DATASET (upside; STOP here does not touch Stage A) ---
B0. Feasibility probe, before any computation. Confirm the non-singlet Southard
file is present and readable; confirm Southard's droplet-overloading
normalization is available or implementable from the published method; report
the cells per factor achievable and how that compares with the 106 of the
singlets file. If any of these fails -> STOP `fulldata_infeasible`, write the
findings with Stage A complete and this stage marked not attempted.
B1. Measure the information content of recovered cells; do not assume it.
Compare the null spread of a pseudobulk of n singlets against the null spread of
a pseudobulk of n cells drawn at the realised singlet-to-recovered ratio, at
n = 106 and n = 390. Report the ratio, and from it an "effective cells per
recovered cell" figure. The project has been assuming 1.0; this measures it.
B2. Re-run A1 to A4 at the achieved depth, in S2A_k20 and in S1. S1 is the point
of this stage: SOUTH6REC showed S1 recovers a 25% approach at 390 cells with
power 0.955, so at that depth the finding can be stated in gene space, with no
basis to defend. Note in the findings that even at 390 cells a full-gap plant
leaves the S1 change in distance positive, so in S1 the depth-matched criterion
is the only one that can call an approach; the SOUTH6 rule cannot, at any
effect size, at that depth.

--- REPORT-ONLY (no decision rule, never a gate) ---
R1. The SOUTH6 rule (change in distance < 0) for every factor, both sets.
R2. The all-1,836 run with m = 1,836, to show the multiple-testing effect.
R3. The SOUTH6 recurring 11 under the new setup, with cell-count percentile.
R4. Spearman of the new statistic against cell count, and displacement ×
    sqrt(cells) across cell-count bins, in the new primary space. In S1 that
    product was flat at 1,273 within 3%, which is the signature of pure sampling
    error. Report whether it is flat in S2A_k20.
R5. For every factor, the frozen ruler score change beside the distance call.
    In SOUTH6, 761 of 1,836 factors scored younger on the ruler while none ended
    closer in S1. Report the same two columns here. This is the preprint's thesis
    in a second dataset and is reported as a diagnostic, not as evidence about
    any factor.

--- WHAT DOES NOT COUNT ---
- the SOUTH6 recurring 11, or any earlier shortlist, as a prior;
- any result from an O2 pair;
- choosing space, criterion, factor set or depth after seeing results;
- Stage B failing as any reflection on Stage A;
- reading a zero in Stage A as evidence that no transcription factor can
  rejuvenate aged cells. The screen uses neonatal foreskin fibroblasts with no
  aged donor and no matched young target. A zero here says the screen contains
  no such factor as measured against a GTEx-anchored young reference; it does not
  test the question in aged cells.

--- STOP RULES ---
On any STOP: write the findings with the key and exit. Never compute after a
STOP; never resume a stopped run.

--- OUTPUT ---
FINDINGS_SOUTH7.md: this text verbatim; the Stage 0 reproduction table; the
sourcing record for the factor set; Stage A per space and set; the gate reading
that fired; Stage B or its STOP; report-only readings labeled as such;
limitations. Tables in results/south7/. Code src/south7_run.py, run from the
repository root. Hand over files, not summaries.
```

## 1. Factor set: sourcing record

| item | value |
|---|---|
| Paper | Southard, Ardy et al., *Comprehensive transcription factor perturbations recapitulate fibroblast transcriptional states*, Nature Genetics 57:2323-2334 (2025), doi:10.1038/s41588-025-02284-1 |
| File | `41588_2025_2284_MOESM4_ESM.xlsx` (Supplementary Tables 1-19), version of record, 16,946,878 bytes, md5 `031d1831dc8fdf26b855f1ba0ce5708b`, downloaded 2026-10-04 from static-content.springer.com |
| Table | sheet `Table S4 Guide Activity`; column semantics from sheet `Table S17 Column descriptions` |
| Columns | `active_fibro` - "guide exhibited phenotype driven clustering (i.e. clustering with guides to same target and/or related genes)"; `expanded_active_fibro` - "guide clustered with active guides" |
| Inclusion rule | a gene is in the primary set if at least one of its guides has `active_fibro` TRUE or `expanded_active_fibro` TRUE. The published classification is used as published; no threshold was applied here and no list was re-derived. |
| Rows read | 10,708 guides |
| Result | `active_fibro`: 659 guides, **212 genes**. `expanded_active_fibro`: 106 guides, **102 genes**. Union: **304 genes**. |

Both published counts reproduce exactly: the paper's text states "212 and 178 transcription factors with strong effects in Hs27 and RPE-1 cells, respectively, plus 102 and 82 additional genes that had single guides that clustered strongly with these gold standards".

Discrepancy, reported and not adjusted. The two flags are not disjoint: 10 genes carry both. The paper's wording ("plus 102 additional genes") implies 212 + 102 = 314 distinct genes; the table yields 304. Each of the 10 has 6 guides, of which 2-4 are `active_fibro` and 1 further guide is `expanded_active_fibro`, so for these the expanded flag marks an additional guide of an already-active gene rather than an additional gene. Genes added by the expanded flag: 92. The 10 are CEBPA, FOXI3, HOXC13, KLF15, NKX1-2, PHOX2A, RREB1, SOX8, TRIM28, ZNF711. Computation was halted at this point and the author chose the union (304) as the primary set; nothing was altered to make the counts match.

All 304 are among the screen's 1,836 targets. Cell counts in the set: median 114 against the screen's 106.

## 2. Stage 0: prior results reproduced with this task's code (gate)

80 checks, 80 gating, 0 gating failures. Counts exact; every other item within |new - stored| <= 1e-4 x max(1, |stored|). Full table: `results/south7/stage0_reproduction.csv`.

| block | items | max abs diff | mismatches | passed |
|---|---|---|---|---|
| obs | 1 | n/a | 0 | yes |
| (a) counts | 6 | 0 | 0 | yes |
| (b) transform | 5 | 0 | 0 | yes |
| (c) spaces | 4 | 1.42e-14 | 0 | yes |
| (d) stage3 S1 O1_Y1 | 4 | 2.27e-13 | 0 | yes |
| (d) stage3 S2_k20 O1_Y1 | 4 | 1.42e-14 | 0 | yes |
| (d) stage3 S2_k50 O1_Y1 | 4 | 1.14e-13 | 0 | yes |
| (d) stage3 S2_k100 O1_Y1 | 4 | 7.82e-14 | 0 | yes |
| (e) stage2b S1 O1_Y1 all78 | 2 | 0 | 0 | yes |
| (e) stage2b S2_k20 O1_Y1 all78 | 2 | 5.55e-17 | 0 | yes |
| (e) stage2b S2_k50 O1_Y1 all78 | 2 | 3.22e-15 | 0 | yes |
| (e) stage2b S2_k100 O1_Y1 all78 | 2 | 5.33e-15 | 0 | yes |
| (f) stage2c S1 O1_Y1 | 1 | 7.82e-14 | 0 | yes |
| (f) stage2c S2_k20 O1_Y1 | 1 | 3.55e-15 | 0 | yes |
| (f) stage2c S2_k50 O1_Y1 | 1 | 7.11e-15 | 0 | yes |
| (f) stage2c S2_k100 O1_Y1 | 1 | 1.78e-15 | 0 | yes |
| (b) SOUTH6REC gate S2_k20 n=106 | 3 | 0 | 0 | yes |
| (b) SOUTH6REC gate S2_k50 n=106 | 3 | 0 | 0 | yes |
| (b) SOUTH6REC gate S2_k100 n=106 | 3 | 0 | 0 | yes |
| (b) SOUTH6REC gate S1 n=106 | 3 | 0 | 0 | yes |
| (b) SOUTH6REC gate S2_k20 n=390 | 3 | 0 | 0 | yes |
| (b) SOUTH6REC gate S2_k50 n=390 | 3 | 0 | 0 | yes |
| (b) SOUTH6REC gate S2_k100 n=390 | 3 | 0 | 0 | yes |
| (b) SOUTH6REC gate S1 n=390 | 3 | 0 | 0 | yes |
| (b) SOUTH6REC gate S2_k20 n=2500 | 3 | 1.39e-17 | 0 | yes |
| (b) SOUTH6REC gate S2_k50 n=2500 | 3 | 8.33e-17 | 0 | yes |
| (b) SOUTH6REC gate S2_k100 n=2500 | 3 | 2.78e-17 | 0 | yes |
| (b) SOUTH6REC gate S1 n=2500 | 3 | 0 | 0 | yes |

(a) covers the SOUTH6 items SOUTH6REC already checked - per-factor and per-guide counts, the TMM/log2-CPM/frozen-z transform, d_f, the origin-target distances, and per-factor delta, cos, frac and ||d|| in S1 and each S2_k on O1_Y1 - plus two additions: `results/south6/perm_D.npy` equals `stage1.npz` d_f on the matched genes exactly, and the per-factor cell counts equal `stage0_cells_per_target.csv` exactly. (b) covers the SOUTH6REC Stage 1 gate: all 12 space-by-depth rows, their power at 25%, their dose-0 median change in distance, and their monotonicity flag.

Nothing in `results/south6/` or `results/south6rec/` was written. `src/south6rec_run.py` was imported unchanged and its Stage 0 reused; its module-level output directory was redirected to `results/south7/` for the duration of that call so no file in `results/south6rec/` was touched.

## 3. The anchor-free basis (S2A)

Fitted exactly as `src/south6_run.py` `build_spaces` fits its basis - `fit_pca` on the frozen-z donor matrix, centred on that matrix's own mean - but on the 308 GTEx fibroblast donors aged 40-59 only. **0 of the 308 is an anchor donor**, against 344 of 652 (52.8%) in the SOUTH6 basis. O1 and Y1 themselves are the frozen TOWARD anchors, unchanged.

| k | gap retention, S2A (anchor-free) | gap retention, SOUTH6 basis | ||O1-Y1|| in S2A | ||O1-Y1|| in SOUTH6 basis |
|---|---|---|---|---|
| 20 | 0.8651 | 0.8873 | 27.9292 | 28.6452 |
| 50 | 0.8973 | 0.9560 | 28.9696 | 30.8643 |
| 100 | 0.9054 | 0.9737 | 29.2292 | 31.4343 |

Removing the circularity costs little, as SOUTH6REC's R5 predicted: retention at k=20 falls from 0.8873 to 0.8651.

## 4. Stage A1: the depth-matched null

The null change in distance depends only on how many cells a pseudobulk holds, so it was computed once per distinct cell count: 212 distinct counts spanning 13-341 cells, R = 200 draws each, 42400 draws per space. **Every factor's floor is exact at its own cell count; the interpolation the pre-registration permits was not used, so there is no interpolation error to report.**

Draws follow the SOUTH6REC primary scheme: a permutation of the 14,013 non-targeting cells cut into disjoint blocks of n, each draw's reference being the complement of that draw, so no draw is inside its own reference.

| space | cells | 5th-pct floor | null median | null SD | mutually disjoint draws | permutations |
|---|---|---|---|---|---|---|
| S1 | 13 | 316.7843 | 353.0495 | 22.4220 | 200 | 1 |
| S1 | 106 | 89.6697 | 93.4181 | 2.7670 | 132 | 2 |
| S1 | 341 | 40.1237 | 41.6270 | 0.8624 | 41 | 5 |
| S2A_k100 | 13 | 75.9453 | 93.1609 | 11.3505 | 200 | 1 |
| S2A_k100 | 106 | 2.8845 | 5.3354 | 1.6543 | 132 | 2 |
| S2A_k100 | 341 | -0.3288 | 0.8543 | 0.7939 | 41 | 5 |
| S2A_k20 | 13 | 33.0091 | 44.5950 | 7.2687 | 200 | 1 |
| S2A_k20 | 106 | -0.4660 | 1.8803 | 1.5012 | 132 | 2 |
| S2A_k20 | 341 | -0.9816 | 0.4011 | 0.8426 | 41 | 5 |
| S2A_k50 | 13 | 62.7044 | 78.2973 | 9.9929 | 200 | 1 |
| S2A_k50 | 106 | 1.6297 | 4.0504 | 1.5536 | 132 | 2 |
| S2A_k50 | 341 | -0.5264 | 0.6205 | 0.7965 | 41 | 5 |

Full table: `results/south7/stage_a1_depth_nulls.csv` (every depth, every space, with overlap figures). The floor rises steeply as cells fall - the reason a floor built from ~7,006-cell half-splits, as SOUTH6's was, is not a floor for a 106-cell factor.

## 5. Stage A2: the permutation null

n_perm = 100,000. Attainable minimum p = 1/(n_perm+1) = 1e-05, which is <= 1e-5 as the pre-registration requires.

What is permuted, stated exactly: for each factor, its own displacement vector d over the 20,264 matched gene coordinates is permuted uniformly at random; the 3,221 ruler coordinates the panel does not measure stay at 0 in every draw, so a draw cannot invent a measurement; ||d|| is preserved in gene space. Draws are generated in chunks of 1,000 from `default_rng([20260914, factor_index])`. This is SOUTH6's construction unchanged - the same generator, seeded the same way, consuming the stream the same way.

Confirmation that it is the same null, not merely a similar one: the first 5,000 draws of each factor's stream were compared against SOUTH6's stored per-factor floors and <=-counts (`results/south6/stage2a_floors.npz`, rows 0-3, its four primary O1_Y1 settings).

| setting | SOUTH6 row | factors compared | max abs floor diff | floors bit-identical | counts identical | count mismatches |
|---|---|---|---|---|---|---|
| S1 | 0 | 304 | 0.000000 | yes | yes | 0 |
| south6basis_S2_k20 | 1 | 304 | 0.000000 | no | yes | 0 |
| south6basis_S2_k50 | 2 | 304 | 0.000000 | no | yes | 0 |
| south6basis_S2_k100 | 3 | 304 | 0.000000 | no | yes | 0 |

Because each factor is seeded independently of every other, the draws do not depend on how many worker processes ran or in what order the factors were taken, and no checkpoint or resume logic was used: each factor's stream runs start to finish in one process.

Multiple testing no longer caps the answer. With m = 304 the smallest attainable q is 0.0030, and with m = 1,836 it is 0.0184; both are below 0.05, so a single factor can clear Benjamini-Hochberg on its own. Under SOUTH6's 5,000 draws, 8 factors had to tie at the minimum p before any could pass, and S2_k20 on the primary pair had 4.

## 6. Stage A3: the call

A factor is APPROACHING in the primary space if its change in distance is below the 5th percentile of its own depth-matched null AND its Benjamini-Hochberg q over the factor set is <= 0.05.

| factor set | m | space | primary | called approaching | below depth floor | q <= 0.05 | SOUTH6 rule (delta<0) | tied at min p | smallest q |
|---|---|---|---|---|---|---|---|---|---|
| primary_304 | 304 | S2A_k20 | yes | 0 | 27 | 0 | 35 | 0 | 0.0578 |
| primary_304 | 304 | S2A_k50 | no | 0 | 23 | 0 | 4 | 0 | 0.3040 |
| primary_304 | 304 | S2A_k100 | no | 0 | 20 | 0 | 1 | 0 | 0.9576 |
| primary_304 | 304 | S1 | no | 1 | 3 | 51 | 0 | 12 | 0.0003 |

### Gate reading that fired (pre-registered wording)

- **the published primary set, S2A_k20: of 304 factors, none is called approaching.** Each factor is judged against a null built at its own cell count, and SOUTH6REC showed this statistic recovers a 25%-of-gap approach at this depth with power 0.990. The screen contains no transcription factor that moves the reference cells a quarter of the way toward young.

Qualification on that wording. The 0.990 figure is SOUTH6REC's, measured in S2_k20 on the 652-donor basis, not in S2A_k20. Because the pre-registered sentence rests on it, the same plant-recovery measurement was repeated in the SOUTH7 primary space (report-only, R6 below): in S2A_k20 at 106 cells the power at a 25% plant is 0.980. The premise holds in the space actually used.

How close it came in primary_304: 27 of 304 factors do fall below their own depth-matched floor, and 35 have a negative change in distance, but none reaches q <= 0.05. The smallest q is 0.0578, held by AFF1, EGR3, FOXN3. The margin is therefore thin, not wide: the permutation test, not the depth floor, is what none of them clears - their displacement direction is not special given its magnitude.

#### S1 (secondary, never primary): 1 called in primary_304

| factor | cells | delta | depth floor | p | q | delta < 0 | high_depth_candidate | guide agreement | identity drop | ruler score change |
|---|---|---|---|---|---|---|---|---|---|---|
| MYC | 99 | 94.0537 | 94.3758 | 0.0001 | 0.0018 | no | no | 0.2147 | 0.0047 | 0.5502 |

**Read this one carefully.** MYC's change in distance in S1 is +94.05, in a space where the depth-matched null at 99 cells has median +98.43 and 5th percentile +94.38. MYC moved **away** from the young target, only less far than sampling noise usually does. The pre-registered criterion labels that APPROACHING because it compares each factor with its own null rather than with zero, which is the remedy the preprint asks for; the label is not a claim that the cells moved toward young. In S1, 0 of 304 factors have a negative change in distance at all, and a full-gap plant at this depth still leaves it positive (R6).

A4 depth flag: the 90th percentile of cell count within the primary set is 192 cells. Every factor's flag is in the per-factor table, whether or not it was called.

Per-factor tables, every factor and every column the pre-registration asks for (change in distance, depth-matched floor, permutation floor, p, q, the call, the SOUTH6 rule, cell count and its percentile, guide agreement, identity drop, cosine, and the call in each secondary space): `results/south7/stage_a3_primary_factors.csv`.

## 7. Report-only readings (no decision rule, never a gate)

### R1 The SOUTH6 rule (change in distance < 0), every factor

| factor set | space | m | delta < 0 (SOUTH6 rule) | below depth floor | called approaching |
|---|---|---|---|---|---|
| primary_304 | S2A_k20 | 304 | 35 | 27 | 0 |
| primary_304 | S2A_k50 | 304 | 4 | 23 | 0 |
| primary_304 | S2A_k100 | 304 | 1 | 20 | 0 |
| primary_304 | S1 | 304 | 0 | 3 | 1 |

### R2 Multiple testing at m = 1,836

The all-1,836 permutation run was not complete when this file was written; the primary set is unaffected, since its q values are computed within the primary set.

### R3 SOUTH6's recurring set under the new setup

The recurring set is the 11 factors labelled toward_young in at least 3 of the 5 SOUTH6 settings that produced any (all of them O2 pairs): CERS3, ETV3L, FOSB, FOXF1, FOXI3, GATA1, HELZ2, JDP2, NEUROD1, OSR2, REST. Median cell count 183 against the screen's median of 106, and 8 of 11 sit above the screen's 90th percentile (167 cells) - reproducing the figures the pre-registration states (183 and 8). 9 of them are in the published primary set.

The table below covers the 9 of the 11 that are in the factor set reported here; CERS3, REST are outside the published primary set and appear only once the all-1,836 run is complete.

| factor | cells | cell pctl in screen | high_depth_candidate | delta | depth floor | q | called approaching | SOUTH6 rule | ruler scores younger |
|---|---|---|---|---|---|---|---|---|---|
| ETV3L | 168 | 90.5773 | no | 9.3303 | -1.0206 | 1.0000 | no | no | yes |
| FOSB | 228 | 98.5294 | yes | 9.1936 | -0.9982 | 1.0000 | no | no | yes |
| FOXF1 | 137 | 77.9956 | no | 3.8218 | -0.7802 | 1.0000 | no | no | no |
| FOXI3 | 171 | 91.1765 | no | 4.2143 | -0.7284 | 1.0000 | no | no | yes |
| GATA1 | 239 | 98.9651 | yes | 7.6969 | -0.9565 | 1.0000 | no | no | yes |
| HELZ2 | 140 | 79.5752 | no | -0.3020 | -0.6944 | 1.0000 | no | yes | no |
| JDP2 | 331 | 99.9455 | yes | 8.2053 | -1.0294 | 1.0000 | no | no | yes |
| NEUROD1 | 185 | 94.2266 | no | 7.5108 | -0.9873 | 1.0000 | no | no | yes |
| OSR2 | 183 | 94.0087 | no | 1.1661 | -0.8752 | 1.0000 | no | no | no |

### R4 Depth dependence of the statistic

| factor set | space | factors | Spearman delta vs cells | Spearman ||d|| vs cells | median ||d||x sqrt(cells) | IQR as % of median |
|---|---|---|---|---|---|---|
| primary_304 | S2A_k20 | 304 | -0.1966 | -0.7128 | 125.0021 | 31.8456 |
| primary_304 | S2A_k50 | 304 | -0.5797 | -0.8512 | 183.4552 | 26.1067 |
| primary_304 | S2A_k100 | 304 | -0.6936 | -0.8805 | 211.2571 | 26.3575 |
| primary_304 | S1 | 304 | -0.9855 | -0.9856 | 1286.7974 | 5.0307 |

By quartile of cell count:

| factor set | space | quartile | factors | median cells | median ||d|| | median ||d||x sqrt(cells) |
|---|---|---|---|---|---|---|
| primary_304 | S2A_k20 | Q1 | 75 | 73.0000 | 16.5447 | 142.2803 |
| primary_304 | S2A_k20 | Q2 | 77 | 99.0000 | 13.3440 | 129.2791 |
| primary_304 | S2A_k20 | Q3 | 76 | 130.0000 | 9.6019 | 110.6581 |
| primary_304 | S2A_k20 | Q4 | 76 | 182.0000 | 8.4396 | 117.8862 |
| primary_304 | S2A_k50 | Q1 | 75 | 73.0000 | 24.1126 | 205.6622 |
| primary_304 | S2A_k50 | Q2 | 77 | 99.0000 | 18.9925 | 188.7917 |
| primary_304 | S2A_k50 | Q3 | 76 | 130.0000 | 14.4680 | 164.0290 |
| primary_304 | S2A_k50 | Q4 | 76 | 182.0000 | 11.4269 | 156.9576 |
| primary_304 | S2A_k100 | Q1 | 75 | 73.0000 | 28.0213 | 237.2582 |
| primary_304 | S2A_k100 | Q2 | 77 | 99.0000 | 21.6202 | 215.4927 |
| primary_304 | S2A_k100 | Q3 | 76 | 130.0000 | 16.7075 | 190.8808 |
| primary_304 | S2A_k100 | Q4 | 76 | 182.0000 | 12.8557 | 180.7462 |
| primary_304 | S1 | Q1 | 75 | 73.0000 | 152.8701 | 1313.5929 |
| primary_304 | S1 | Q2 | 77 | 99.0000 | 130.0276 | 1287.2074 |
| primary_304 | S1 | Q3 | 76 | 130.0000 | 112.2919 | 1281.6193 |
| primary_304 | S1 | Q4 | 76 | 182.0000 | 94.1974 | 1269.6317 |

The answer, stated rather than left to the table. In S1 the product is flat: median 1287, moving only from 1314 to 1270 between the lowest and highest cell-count quartile, with an interquartile range of 5.0% of the median. That reproduces the ~1,273-within-3% signature the pre-registration describes, and it is what displacement looks like when it is pure sampling error. In S2A_k20 it is **not** flat: 142 to 118 across the same quartiles, an interquartile range of 31.8% of the median. So the compressed space is not merely a rescaling of the sampling error, though depth still drives much of it: Spearman of ||d|| against cell count is -0.71 there, against -0.99 in S1. Any shortlist from either space has to be checked against cell count before it is believed, which is what the A4 flag is for.

### R5 The frozen ruler beside the distance call

- primary_304: 126 of 304 factors score younger on the frozen ruler (ruler score change < 0), while 35 end closer to young in S2A_k20 by the SOUTH6 rule and 0 are called approaching. In SOUTH6, 761 of 1,836 scored younger on the ruler while none ended closer in S1.

This is the preprint's thesis in a second dataset: a one-dimensional age score falls for many factors that do not move closer to a young target of the same cell type. It is reported as a diagnostic and is not evidence about any factor. The ruler was not refit.

### R6 Plant recovery in the SOUTH7 spaces (an addition, not in the pre-registration's report-only list)

The SOUTH6REC recovery control, repeated in the anchor-free spaces so the gate's PASS wording does not rest on a figure measured in a different basis. Same construction as SOUTH6REC: a known fraction of the measured gap planted into real non-targeting pseudobulks, 200 draws per dose, power against the 5th percentile of the 0% draws.

| space | cells | monotone | median delta 0% | median delta 25% | power 25% | power 50% | power 100% |
|---|---|---|---|---|---|---|---|
| S2A_k20 | 106 | yes | 1.8302 | -3.9401 | 0.9800 | 1.0000 | 1.0000 |
| S2A_k50 | 106 | yes | 3.7580 | -1.7101 | 0.9550 | 1.0000 | 1.0000 |
| S2A_k100 | 106 | yes | 5.1828 | -0.2244 | 0.9300 | 1.0000 | 1.0000 |
| S1 | 106 | yes | 93.6030 | 92.1517 | 0.1150 | 0.3400 | 0.4700 |
| S2A_k20 | 390 | yes | 0.2267 | -5.9227 | 1.0000 | 1.0000 | 1.0000 |
| S2A_k50 | 390 | yes | 0.4957 | -5.9059 | 1.0000 | 1.0000 | 1.0000 |
| S2A_k100 | 390 | yes | 0.6827 | -5.7418 | 1.0000 | 1.0000 | 1.0000 |
| S1 | 390 | yes | 37.6739 | 34.7980 | 0.9550 | 1.0000 | 1.0000 |

This carries no decision rule. It states what the statistic could have detected at the screen's depth in the space the gate uses.

## 8. Stage B: full dataset

**STOP `fulldata_infeasible` at B0, the feasibility probe, before any computation.** Stage A is complete and unaffected.

- The non-singlet file `data/raw/southard/fibroblast_CRISPRa_final_pop.h5ad` is **not present** on disk. Only the singlets file is: `fibroblast_CRISPRa_final_pop_singlets_normalized_log1p.h5ad`, 8,048,055,098 bytes, md5 confirmed.
- It exists in the same Zenodo record (10.5281/zenodo.15200179) as `fibroblast_CRISPRa_final_pop.h5ad`, listed at 9.7 GB. The pre-registration admits no new expression data beyond the Southard files already on disk, so it was not downloaded.
- Droplet-overloading normalization was therefore not assessed against data. The singlets file does carry the per-cell columns the published method produces (`pairwise_coef`, `pairwise_UMI_count`, `pairwise_singlet_equivalent_UMI_count`, `final_pairwise_UMI_count`), so the route is visible, but it cannot be exercised without the non-singlet matrix.
- Cells per factor achievable: not measurable without the file. The singlets file gives a median of 106 cells per factor.

B1 and B2 were not attempted. Nothing in Stage A depends on them, and no Stage B number exists.

## 9. Limitations

- **A zero in Stage A is not evidence that no transcription factor can rejuvenate aged cells.** The screen uses neonatal foreskin fibroblasts with no aged donor and no matched young target. A zero says the screen contains no such factor as measured against a GTEx-anchored young reference; it does not test the question in aged cells.
- Recovery of a planted effect says what the statistic could have seen, not that any factor has such an effect. A plant points exactly along the gap; real displacement vectors are longer than the gap and mostly point off it.
- The depth-matched null is built by resampling one pool of 14,013 non-targeting cells from one line. Its spread is the spread of resampling that pool, and random cells mix all 78 non-targeting guides, so it carries no guide-to-guide variation of the kind a factor's ~6 guides carry.
- The plant and the null both enter in z, after the transform. A real effect would pass through counts, TMM and the log and would carry its own cell-to-cell variance.
- S2A removes the anchor donors from the basis but not from the z-scoring: the frozen mu and sd come from all 652 GTEx donors, the anchors included. That is the frozen ruler and was not refit.
- O1 and Y1 are GTEx fibroblast tissue donors, not cells of this line; every change in distance is a counterfactual that assumes the effect transfers. The O2 pairs were not run: their axis is a 468-unit platform gap, not an age gap.
- Unmeasured genes are assumed unmoved. 3,221 of the ruler's 23,485 genes are outside the panel, and 164 matched genes have no counts.
- Guide agreement and identity drop are carried from SOUTH6's stored Stage 1, which Stage 0 reproduced; they are gene-space quantities and do not depend on the basis.

## 10. Files

`src/south7_run.py`. In `results/south7/`: `PREREG.flag`, `southard_table_s4_hs27.csv` (the extracted published columns), `factor_set.csv`, `stage0_reproduction.csv`, `basis.npz`, `basis_retention.csv`, `stage_a1_observed.csv`, `stage_a1_depth_nulls.csv`, `stage_a1_null_draws.npz`, `stage_a2_perm_*.npz`, `stage_a2_same_null_check.csv`, `stage_a3_primary_factors.csv`, `stage_a3_calls.csv`, `r3_recurring11.csv`, `r4_depth_relation.csv`, `r6_recovery_primary_space.csv`, `stage_b0_probe.csv`, `south7_run.log`.

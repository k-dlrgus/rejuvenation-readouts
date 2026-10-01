**SUPERSEDED BY SOUTH3 (see FINDINGS_SOUTH3.md): SOUTH3 re-runs this analysis on the same data and the same pre-registration with three changes — n_perm = 20000 for the Stage 2A per-factor permutation null so BH across 1,836 factors can resolve q, guide agreement computed on V_all vectors always and applied to both variants, and a report-only section on ||d||, ||d|| / ||target - origin|| and the small-displacement identity. Nothing below is changed.**

# FINDINGS_SOUTH2 — rank Southard Hs27 CRISPRa transcription factors by approach to a young target, from the per-guide regression output

**Status:** Ran to completion. No STOP key fired. Report-only flags: `scale_unknown`, `coverage_low`.

**Result:** `toward_young` is 0 in 19 of the 28 settings and non-zero in 9, ranging 189-231 of 1836 factors. Every non-zero setting is V_all; V_sig returns 0 everywhere. The counts are not pooled and no single setting is primary. Read the "What produced these counts" section before any count: several settings are zero for geometric or resolution reasons rather than because factors were tested and found null, and the O2 settings are dominated by a cross-dataset offset.

SOUTH2 supersedes SOUTH. SOUTH halted at Stage 0 with `no_raw_counts` because its pre-registration described this file as a count matrix; the file holds per-guide regression output. `FINDINGS_SOUTH.md` was not edited except for one line at its top marking it superseded.

The pre-registration below was written verbatim to `FINDINGS_SOUTH2.md` and `results/south2/PREREG.flag`, and `results/south2/PREREG.flag` was created, before any statistic was computed. Nothing was re-tuned after numbers existed.

## Pre-registration (verbatim, written before any SOUTH2 statistic)

```
TASK: SOUTH2

Nature of this task: a purely computational reanalysis of a published,
publicly released human cell dataset already on disk in this repository.
No laboratory work, no protocol, no biological material, no new data.
The output is a ranked table and a set of statistics.

SOUTH halted at Stage 0 with key `no_raw_counts` because its pre-registration
described the input as a count matrix and the pinned file holds per-guide
regression output instead. SOUTH2 supersedes SOUTH with an input description
that matches the file. Do not edit FINDINGS_SOUTH.md; add one line at its top
marking it superseded by SOUTH2, and leave everything else in it unchanged.

Write FINDINGS_SOUTH2.md, PROGRESS_SOUTH2.md, results/south2/*, src/south2_run.py.
Do not modify any other existing FINDINGS_*.md, PROGRESS_*.md, or src/ file.
Import from src/toward_run.py, src/same_run.py, src/seng_run.py,
src/lowdim_common.py, src/md3_idtype.py; do not change them.
Do not refit the frozen ruler.
Seeds: 20260914 (splits, nulls, permutations), 20260918 (bootstrap).
n_null=200, n_perm=200, n_boot=200. Uncalibrated R2 is never a gate.

WRITE THE PRE-REGISTRATION FIRST. Copy this whole file verbatim into
FINDINGS_SOUTH2.md and touch results/south2/PREREG.flag before computing any
statistic. Stage 0 may run after the flag exists. Nothing is re-tuned after
numbers exist. Every STOP key halts the task and is reported as written.

--- DATA (described as the file actually is) ---
results/survey/opened/southard/fibroblast_CRISPRa_mean_pop.h5ad
1,720,196,368 bytes, md5 ba44c7813903bb5df900348d6b0d589a.
Rows are per-guide perturbations; X holds per-gene effect estimates from the
source's own regression, with p and adj_p layers and a masked layer.
This task uses those effect estimates directly as perturbation vectors.
It does not use L4A_matrix.h5 and does not download anything.

--- STAGE 0: INVENTORY (descriptive only, no statistic) ---
Report:
- exact shape; what one row is; what the column index holds (symbols,
  Ensembl ids, versioned Ensembl ids, or mixed).
- the units and log base of X, from uns, var, obs, or any bundled
  documentation. If undeterminable, record key `scale_unknown` as a
  REPORT-ONLY flag, not a stop; note that the Stage 2 floor is
  magnitude-matched and therefore invariant to a global scale factor,
  while the Stage 2 MDA is not.
- what the masked layer marks, and how many entries it covers. Masked
  entries are set to 0 and the count is reported.
- ruler-gene coverage computed THREE ways: by symbol, by Ensembl id with
  version suffixes stripped, and by the better of the two after
  src/md3_idtype.py normalisation. For each, report the number of the
  23,485 ruler genes matched AND the fraction of the ruler's total
  absolute weight matched. Use the mapping with the highest weight
  fraction for all later stages and name it explicitly.
  Record flag `coverage_low` if that weight fraction is below 0.25.
- number of distinct transcription factors, guides per factor
  (min, median, max), and the count with at least 2 guides.
- the untargeted control rows, identified by the `off-target` label and by
  empty target_gene_id / NaN target_expr / active=False. Report their count
  and, for each, sequence_driven status, de_genes count, and cell count.
  These are used only as a secondary floor and are flagged as not inert.
STOP, key `not_an_effect_matrix`, if X is not per-gene effect estimates.
STOP, key `too_few_factors`, if fewer than 500 factors have 2+ guides.

--- STAGE 1: VECTORS AND SPACES ---
Per-factor effect vector: mean of that factor's guide rows.
Conversion to frozen-ruler coordinates: divide each gene's effect by that
gene's frozen training sd, giving a displacement in frozen-z units. Genes
absent from the table are set to 0, meaning "this factor does not move this
gene"; state that this is a weaker assumption than setting an expression
level to zero, and report how many coordinates it fixes.
Two gene variants, both pre-registered, both reported, never combined:
  V_all  = every matched gene.
  V_sig  = matched genes with adj_p <= 0.05 kept, the rest set to 0.
Two spaces, both reported:
  S1 gene space, frozen-ruler coordinates as above.
  S2 low-dim, PCA fit ONLY on the GTEx fibroblast donor matrix, never on
  this dataset or the targets, k in (20, 50, 100), applied unchanged.
  Report every k; do not choose one after seeing results.

--- TARGETS AND ORIGINS (fixed before anything is scored) ---
Young target, reported separately, never combined:
  Y1 = GTEx fibroblast donors aged 20-39 (existing TOWARD anchor).
  Y2 = the same donors as a region: mean and covariance, S2 only,
       Mahalanobis distance, membership radius at the 95th percentile of
       held-out young donors' own distances. Report the radius.
Aged origin, reported separately, never combined:
  O1 = GTEx fibroblast donors aged 60-79 (existing anchor).
  O2 = the Lu et al. aged donor GM00731 untreated day-0 fibroblast
       pseudobulk, as built in SAME.
Report each origin-to-target distance before scoring anything.

--- STAGE 2: DETECTION LIMIT (before any factor is ranked) ---
A. Primary floor, magnitude-matched random direction. For each factor,
   permute its own vector's entries across genes 200 times and compute
   delta = ||(origin + d_perm) - target|| - ||origin - target||.
   The per-factor floor is the 5th percentile of that distribution.
   Also report a pooled floor built from a vector at the median factor norm.
B. Secondary floor, controls. Same delta using the untargeted control rows
   as pseudo-effects. Report it, flag it as built on non-inert controls,
   and do not use it as the gate.
C. Dose-response. Add the true young-minus-old direction to the aged origin
   at f = 0, 0.05, 0.1, 0.25, 0.5, 1.0, scaled to the median factor norm.
   MDA = the smallest f whose delta falls below the primary floor.
   STOP, key `no_detection_power`, if no f <= 1.0 clears the floor in any
   space or variant. Report Stages 0-2 only and state plainly that this
   data cannot support a ranking.
Never subtract a floor from a later number. It is a bar to clear.

--- STAGE 3: RANKING (only if Stage 2 passed) ---
For each factor with 2+ guides, each origin/target pair, each space, each
gene variant, report:
  cos    = dot(d, v) / (||d|| ||v||), v = target - origin
  frac   = dot(d, v) / ||v||^2
  delta  = ||(origin + d) - target|| - ||origin - target||
  for Y2 also: change in Mahalanobis distance, and whether origin + d lands
  inside the membership radius.
Guide reliability: compute each guide's own vector and report the median
pairwise cosine between guides of the same factor. Flag `guides_disagree`
if that median is <= 0. Factors so flagged are reported but never labelled
a candidate.
Also per factor: the identity check on COL1A1, COL1A2, FN1, LUM, PDGFRA,
PDGFRB, POSTN, PRRX1, SERPINH1, VIM, using only the genes present (POSTN is
absent; report the panel actually used and its size), flagging
`identity_loss` if the drop exceeds 0.5; and the frozen ruler score change
as a context column only.
Nulls: the Stage 2A per-factor permutations, one-sided empirical p,
BH-FDR across factors within each space and variant.
Bootstrap over guides where a factor has 3+; flag any interval that does not
contain its point estimate as INVALID and never use it as a gate.
Label a factor `toward_young` only if delta is negative AND below its own
Stage 2A floor AND q <= 0.05 AND identity_loss is not flagged AND
guides_disagree is not flagged.
Report the count in every space, variant and origin/target pair. Zero is a
valid and reportable answer.
If `coverage_low` fired in Stage 0, the candidate table is titled
"exploratory" and every candidate carries that word.

--- STAGE 4: DISAGREEMENT WITH DIRECTIONAL READOUTS ---
Rank all factors three ways: by delta, by cos, and by drop in the frozen
ruler score. Report Spearman correlation between each pair and the overlap
of the top 20. State how many of the top 20 by cosine are not
`toward_young` by delta.

--- WHAT DOES NOT COUNT ---
- cosine or ruler-score drop alone as evidence of approach.
- any delta that does not clear the measured primary floor.
- subtracting a floor instead of clearing it.
- choosing k, gene variant, mapping, origin or target after seeing results.
- using the control-row floor as the gate.
- dropping or relabelling factors because they look odd; report them.
- pooling spaces, variants, origins or targets into one number.
- substituting a different file for the pinned one.
- claims about the source paper's own conclusions; this task does not test them.

--- OUTPUT ---
FINDINGS_SOUTH2.md: this file verbatim; Stage 0 inventory including all three
coverage numbers; conversion and its assumption; Stage 2 floor, control floor
and MDA tables; Stage 3 full table plus a top-20 table with factor, delta,
margin over floor, q, guide agreement, identity flag; Stage 4 disagreement
table; fired keys; limitations, which must include: inputs are regression
estimates, not expression, so the frozen count pipeline was not applied;
effects are measured in one neonatal fibroblast line with no aged cells, so
every delta is a counterfactual assuming the effect transfers to an aged cell;
unmeasured genes are assumed unmoved; guide counts are small; the log base of
the input may be unresolved; the identity panel is incomplete.
PROGRESS_SOUTH2.md: stop status, next action.
Record every failure verbatim. Do not substitute columns or repair rows.
```

## Stage 0 inventory (descriptive only; no statistic)

- File: `<repo>\results\survey\opened\southard\fibroblast_CRISPRa_mean_pop.h5ad`
- Size: 1,720,196,368 bytes (pinned 1,720,196,368), md5 `ba44c7813903bb5df900348d6b0d589a` (pinned `ba44c7813903bb5df900348d6b0d589a`). Both match; no file was substituted. `L4A_matrix.h5` was not opened and nothing was downloaded.
- Exact shape: **10,916 rows x 4,914 genes**, dtype `float64`, dense. Layers: `adj_p`, `masked`, `p`, each the same shape as X. `uns`: empty. `obsm`: empty. `varm`: empty.

### What one row is

**per-guide perturbation row (the source's own mean-population regression output).** `guide_identity` is unique on every row (10,916 values for 10,916 rows), so a row is one guide, not a per-factor mean. There is no cell-barcode column (`cell_barcode_column=False`), only a per-row `cell_count` (min 1, median 59, max 1369, total 762,352), so a row is not a cell either. Each row is the guide-level output of the source's own regression over the cells carrying that guide.

### What the column index holds

- unversioned Ensembl gene ids in var/gene_id (4,914/4,914 match ENSG, 0 carry a version suffix) with HGNC symbols in var/gene_name (4,914/4,914 symbol-shaped). The h5ad row index of var is a positional integer, not an identifier.

### Units and log base of X

- **Undeterminable from the file.** `uns` is empty; `X` carries only `{'encoding-type': 'array', 'encoding-version': '0.2.0'}` as attributes, which state an encoding and no unit; `var` holds only `mean`, `std`, `cv`, `fano`, `excess_cv` and related columns with no unit statement; `obs` states no unit; and the only documentation bundled in the folder is three plain gene-symbol lists (`hs27_lambert_tfs.txt`, `hs27_targets.txt`, `mean_pop_targets.txt`), which carry no units.
- Descriptive diagnostics that were tried and did **not** settle it: Spearman of the per-gene mean |effect| against the source's own `var/mean` = +0.645, against `var/std` = +0.684, and of mean |effect| / `var/std` against `var/mean` = -0.991. |effect| grows with expression but sub-linearly, which is consistent with more than one scale. No claim about a log base is made.
- Key `scale_unknown` is recorded as a **REPORT-ONLY flag, not a stop**, exactly as pre-registered. The Stage 2A primary floor and the Stage 2B control floor are magnitude-matched — they permute or re-use the same vector, so a global scale factor on X multiplies the factor vector and its null identically and the comparison is invariant to it. **The Stage 2C MDA is not invariant**: it scales the true young-minus-old direction to the median factor norm, and that norm moves with any global rescaling of X, so the MDA would move too.

### The masked layer

- `layers/masked` is X with masked entries written as NaN. It covers **2,247 entries** spread over 2,247 rows (between 1 and 1 per masked row).
- What it marks: **each row's own targeted gene**. Of the 2,247 masked entries, 2,247 sit exactly at the column of that row's `target_gene_id` and 0 sit anywhere else. The count equals the number of rows whose `target_gene_id` is non-empty, i.e. the rows whose targeted gene is itself one of the measured genes; the source masks the guide's direct effect on its own target.
- Masked entries are **set to 0** before anything is derived, and the count is reported above. Nothing else in X is altered.

### Ruler-gene coverage, computed three ways

All three are against the 23,485 frozen-ruler genes. "Weight fraction" is the share of the ruler's total absolute weight sum(|w|) carried by the matched genes.

| mapping | genes_matched | of_ruler_genes | gene_fraction | weight_matched | weight_total | weight_fraction |
|---|---|---|---|---|---|---|
| symbol | 4,834 | 23,485 | 0.2058 | 13.323 | 110.264 | 0.1208 |
| ensembl_version_stripped | 4,906 | 23,485 | 0.2089 | 13.534 | 110.264 | 0.1227 |
| best_of_two_after_md3_idtype | 4,907 | 23,485 | 0.2089 | 13.536 | 110.264 | 0.1228 |

- **Mapping used for all later stages: `best_of_two_after_md3_idtype`**, named explicitly because it has the highest weight fraction (0.1228; ties would go to more genes matched).
- `best_of_two_after_md3_idtype` is Ensembl-id-first with a symbol fallback (`seng_run.ruler_column_index`), applied after `src/md3_idtype.py` confirmed that `var/gene_id` is an Ensembl column (kind `ensembl`) and `var/gene_name` a symbol column (kind `symbol`).
- Flag `coverage_low` (weight fraction < 0.25): **FIRED**.

### Transcription factors and guides

- Distinct transcription factors: **1,836** (the `off-target` label is the untargeted control arm and is not counted as a factor).
- Guides per factor: min **2**, median **6**, max **6**.
- Factors with at least 2 guides: **1,836** (the `too_few_factors` threshold is 500).

### The untargeted control rows

- Count: **17** rows labelled `off-target`. All 17 of them also carry an empty `target_gene_id`, a NaN `target_expr`, and `active=False`.
- That evidence triple is **necessary but not sufficient**: 8,156 targeted guide rows meet it too. `target_gene_id` is filled only when the targeted gene is one of the 4,914 measured genes, `target_expr` likewise, and most guides are `active=False`, so an inactive guide against an unmeasured target looks the same on those three columns. The control arm is therefore identified by the `off-target` label, with the triple as a consistency check that every labelled row passes. No targeted row was reassigned to the controls.
- They are **flagged as not inert** and are used only as the Stage 2B secondary floor, never as the gate.

| guide_identity | sequence_driven | de_genes | cell_count |
|---|---|---|---|
| off-target_GAAGCAGCATGAATACGCCG | True | 304 | 1,369 |
| off-target_GACGAATGAAGCGTCGATAA | False | 1 | 684 |
| off-target_GACTGAAAGCCGATATCGGG | True | 375 | 341 |
| off-target_GAGATATGAGGCGACGATAT | False | 0 | 516 |
| off-target_GCCTAAATACTATTCGCGGA | True | 13 | 479 |
| off-target_GCTGTTTCGACCCGTCGAAT | False | 0 | 373 |
| off-target_GGAATCCAGCTCGACGACCA | False | 236 | 659 |
| off-target_GGAATCTTCGCGTAACGAGC | False | 0 | 433 |
| off-target_GGATTCTGGAAAGTCGGTCG | False | 0 | 444 |
| off-target_GGCTTCTACACCGCGATGAC | False | 94 | 611 |
| off-target_GGGAACAGGGGCGGTCCGTA | False | 20 | 823 |
| off-target_GTAAATTCTCGCGTAACGTT | False | 0 | 238 |
| off-target_GTCCATAGGGTCTAGCGCCG | False | 64 | 653 |
| off-target_GTCGTTATCTCGCTATTTCG | False | 4 | 411 |
| off-target_GTGGATTATCTGACGCGAAT | False | 5 | 363 |
| off-target_GTGGTTATACCCGACTAGAC | False | 0 | 255 |
| off-target_GTTCCTGTTGGGTCGCGAAT | False | 0 | 397 |

3 of 17 control rows are flagged `sequence_driven` and 10 have `de_genes` > 0, so they move the transcriptome and are not an inert baseline.

### The two Stage 0 STOP tests

- `not_an_effect_matrix`: **not fired**. X is per-gene effect estimates. Evidence: the file carries `p` and `adj_p` layers of exactly X's shape, which only exist for a per-entry test; 30,370,039 of 53,641,224 X entries are negative (56.6%, impossible for counts); only 0.00% of entries are integral; X spans -9.156 to 184.805; and the source's own `obs/de_genes` column is reproduced exactly by counting `adj_p <= 0.05` per row on **10,916 of 10,916 rows** (all of them), which ties X's layers to the source's own differential-expression call.
- `too_few_factors`: **not fired**. 1,836 factors carry >= 2 guides, above the threshold of 500.

## Stage 1 conversion, gene variants, and spaces

### Conversion to frozen-ruler coordinates, and its assumption

- Per-factor effect vector = the **mean of that factor's guide rows** (no weighting by cell count, no re-estimation).
- Each gene's effect is divided by that gene's **frozen training sd** from `frozen_ruler_ridge_raw.npz`, giving a displacement in frozen-z units. The frozen ruler is loaded, never refit. sd was floored at 1e-12 on 0 ruler genes.
- Genes absent from the table are set to **0**, meaning "this factor does not move this gene". That fixes **18,578 of the 23,485 coordinates** (79.1%), leaving 4,907 coordinates free to move.
- This is a **weaker assumption than setting an expression level to zero**: it says the perturbation leaves the unmeasured gene where it already was, so the unmeasured coordinates contribute an identical, perturbation-independent offset to the origin and to every candidate end point. Setting an expression level to zero would instead assert that the gene is not expressed, which would move the point itself and would corrupt the origin as well as the displacement.

### The two gene variants (both reported, never combined)

- **V_all**: every matched gene (4,907 coordinates).
- **V_sig**: matched genes with `adj_p <= 0.05` kept, the rest set to 0. The threshold is applied **per guide row** and the surviving effects are then averaged, so each guide contributes only its own significant effects and the per-guide vectors used for guide agreement and for the bootstrap are the same objects that build the factor vector.

### The two spaces (both reported)

- **S1**: gene space, the 23,485 frozen-ruler coordinates in frozen-z units as above.
- **S2**: low-dim. PCA fit **only** on the GTEx fibroblast donor matrix (652 donors x 23,485 frozen-ruler genes, in frozen-z, centred on the GTEx donor mean), then applied unchanged. Loadings are orthonormal. Every k in (20, 50, 100) is reported; no k is chosen after seeing results.

| k | cumulative_fraction_of_GTEx_donor_variance |
|---|---|
| 20 | 0.6337 |
| 50 | 0.7098 |
| 100 | 0.7621 |

Stated plainly because it bears on how S2 should be read: Y1 and Y2 are subsets of the very GTEx donor matrix the PCA is fit on, so the S2 basis is **not** independent of the young target. It is independent of the Southard data and of O2, which is what the pre-registration forbids fitting on. No alternative basis was tried and none was chosen after seeing results.

## Targets, origins, and the distances between them

- **Y1** = GTEx fibroblast donors aged 20-39, the existing TOWARD anchor (n=113; bins ['20-29', '30-39']).
- **Y2** = the same donors as a region in S2 only: mean and covariance fit on 56 young donors, membership radius at the 95th percentile of the 57 **held-out** young donors' own Mahalanobis distances (split seed 20260914).
- **O1** = GTEx fibroblast donors aged 60-79, the existing anchor (bins ['60-69', '70-79']).
- **O2** = the Lu et al. aged donor GM00731 untreated day-0 fibroblast pseudobulk, as built in SAME (`results/same/panel_z.npz`, row `O`); 3,132 of its coordinates are exactly 0 because those genes are absent from that dataset and were zeroed by SAME's own missing-gene rule.

Membership radii, and where Y2 is undefined:

| k | defined | n_fit_donors | covariance_rank | condition_number | membership_radius |
|---|---|---|---|---|---|
| 20 | True | 56 | 20 | 2.207e+02 | 11.8805 |
| 50 | True | 56 | 50 | 1.597e+04 | 59.2987 |
| 100 | False | 56 | 55 | NA | NA |

At k=50 the covariance is non-singular but rests on 56 donors in 50 dimensions (condition number 1.597e+04); read Y2 at this k as fragile.

At k=100 the covariance of 56 fit donors has rank 55, so the Mahalanobis metric is singular. **Y2 is reported as undefined at k=100.** No pseudo-inverse and no ridge term were substituted, and no k was dropped to make a number appear. Y1 at k=100 is unaffected and is reported.

Every origin-to-target distance, reported before anything was scored:

| space | origin | target | metric | origin_to_target_distance | origin_to_target_mahalanobis | membership_radius | origin_inside_radius |
|---|---|---|---|---|---|---|---|
| S1 | O1 | Y1 | euclidean_frozen_z | 32.2839 | NA | NA | NA |
| S1 | O2 | Y1 | euclidean_frozen_z | 468.3878 | NA | NA | NA |
| S2_k20 | O1 | Y1 | euclidean_pc | 28.6452 | NA | NA | NA |
| S2_k20 | O1 | Y2 | euclidean_pc+mahalanobis | 33.2100 | 1.0554 | 11.8805 | True |
| S2_k20 | O2 | Y1 | euclidean_pc | 196.4504 | NA | NA | NA |
| S2_k20 | O2 | Y2 | euclidean_pc+mahalanobis | 204.3757 | 14.2055 | 11.8805 | False |
| S2_k50 | O1 | Y1 | euclidean_pc | 30.8643 | NA | NA | NA |
| S2_k50 | O1 | Y2 | euclidean_pc+mahalanobis | 35.1680 | 6.0114 | 59.2987 | True |
| S2_k50 | O2 | Y1 | euclidean_pc | 259.3335 | NA | NA | NA |
| S2_k50 | O2 | Y2 | euclidean_pc+mahalanobis | 265.1887 | 127.0348 | 59.2987 | False |
| S2_k100 | O1 | Y1 | euclidean_pc | 31.4343 | NA | NA | NA |
| S2_k100 | O1 | Y2 | euclidean_pc+mahalanobis | 35.8353 | NA | NA | NA |
| S2_k100 | O2 | Y1 | euclidean_pc | 272.0445 | NA | NA | NA |
| S2_k100 | O2 | Y2 | euclidean_pc+mahalanobis | 277.4758 | NA | NA | NA |

The true young-minus-old direction has ||v|| = 32.2839 in S1; restricted to the 4,907 measured coordinates its norm is 15.7557, i.e. 48.8% of the full norm. Both are used in the Stage 2C dose-response and both are reported.

## Stage 2 detection limit (computed before any factor was ranked)

### A. Primary floor — magnitude-matched random direction

For each factor, its own vector's entries are permuted across the 4,907 matched gene coordinates 200 times (seed 20260914, an independent generator per factor), and `delta = ||(origin + d_perm) - target|| - ||origin - target||` is recomputed each time. The per-factor floor is the 5th percentile of that distribution. Permuting preserves ||d|| exactly, so the floor is magnitude-matched by construction.

The pooled floor pools every factor's permuted vectors after rescaling each to the median factor norm, and takes the 5th percentile of that pooled distribution. It is the bar the Stage 2C dose-response has to clear.

**Implementation decision, stated because the pre-registration says only "across genes":** the permutation is over the 4,907 measured coordinates, not over all 23,485 ruler coordinates. The other 18,578 are structurally fixed at 0 by the conversion rule, so permuting effects into them would place the null in a subspace the experiment cannot report on and would no longer be a like-for-like random direction for the real effect. This was fixed before any floor was computed and was not revisited.

Distribution of the per-factor floor across factors:

| variant | space | pair | floor_min | floor_p25 | floor_median | floor_p75 | floor_max | n_factors |
|---|---|---|---|---|---|---|---|---|
| V_all | S1 | O1_Y1 | +1.33906 | +4.18368 | +5.16862 | +6.74920 | +88.62504 | 1,836 |
| V_all | S1 | O2_Y1 | -1.23924 | -0.49288 | -0.31466 | -0.13458 | +11.66674 | 1,836 |
| V_all | S2_k20 | O1_Y1 | -2.37348 | -0.96756 | -0.76244 | -0.56149 | +0.27454 | 1,836 |
| V_all | S2_k20 | O1_Y2 | -1.89105 | -0.62225 | -0.51696 | -0.41097 | +0.05721 | 1,836 |
| V_all | S2_k20 | O2_Y1 | -1.52881 | -0.35516 | -0.30347 | -0.26425 | -0.12227 | 1,836 |
| V_all | S2_k20 | O2_Y2 | -1.43213 | -0.31590 | -0.27253 | -0.23984 | -0.13753 | 1,836 |
| V_all | S2_k50 | O1_Y1 | -1.93954 | -0.86213 | -0.69254 | -0.52165 | +0.24777 | 1,836 |
| V_all | S2_k50 | O1_Y2 | -1.44457 | -0.56192 | -0.47352 | -0.38372 | +0.04976 | 1,836 |
| V_all | S2_k50 | O2_Y1 | -1.64296 | -0.42128 | -0.35417 | -0.30091 | -0.10962 | 1,836 |
| V_all | S2_k50 | O2_Y2 | -1.62913 | -0.38768 | -0.32974 | -0.28244 | -0.13736 | 1,836 |
| V_all | S2_k100 | O1_Y1 | -1.29123 | -0.77523 | -0.63770 | -0.48835 | +0.24943 | 1,836 |
| V_all | S2_k100 | O1_Y2 | -1.00931 | -0.49510 | -0.43540 | -0.36338 | +0.07461 | 1,836 |
| V_all | S2_k100 | O2_Y1 | -1.98441 | -0.59315 | -0.48706 | -0.38824 | -0.01194 | 1,836 |
| V_all | S2_k100 | O2_Y2 | -1.87086 | -0.56086 | -0.46211 | -0.37199 | -0.02915 | 1,836 |
| V_sig | S1 | O1_Y1 | +0.00396 | +0.09209 | +0.16433 | +0.30606 | +13.07451 | 1,836 |
| V_sig | S1 | O2_Y1 | -0.15819 | -0.00367 | +0.00416 | +0.01695 | +0.90306 | 1,836 |
| V_sig | S2_k20 | O1_Y1 | -0.50337 | -0.00670 | -0.00036 | +0.00907 | +0.40992 | 1,836 |
| V_sig | S2_k20 | O1_Y2 | -0.42742 | -0.02320 | -0.01631 | -0.01213 | +0.13174 | 1,836 |
| V_sig | S2_k20 | O2_Y1 | -0.31318 | -0.04480 | -0.03348 | -0.02580 | -0.00747 | 1,836 |
| V_sig | S2_k20 | O2_Y2 | -0.30907 | -0.04840 | -0.03588 | -0.02754 | -0.00773 | 1,836 |
| V_sig | S2_k50 | O1_Y1 | -0.44310 | -0.00837 | -0.00277 | +0.00503 | +0.37153 | 1,836 |
| V_sig | S2_k50 | O1_Y2 | -0.37233 | -0.02335 | -0.01674 | -0.01232 | +0.12032 | 1,836 |
| V_sig | S2_k50 | O2_Y1 | -0.35396 | -0.04197 | -0.03152 | -0.02460 | -0.00683 | 1,836 |
| V_sig | S2_k50 | O2_Y2 | -0.35292 | -0.04471 | -0.03322 | -0.02576 | -0.00699 | 1,836 |
| V_sig | S2_k100 | O1_Y1 | -0.42154 | -0.00807 | -0.00283 | +0.00480 | +0.39950 | 1,836 |
| V_sig | S2_k100 | O1_Y2 | -0.34073 | -0.02241 | -0.01617 | -0.01157 | +0.15144 | 1,836 |
| V_sig | S2_k100 | O2_Y1 | -0.38253 | -0.03024 | -0.02210 | -0.01707 | +0.01565 | 1,836 |
| V_sig | S2_k100 | O2_Y2 | -0.37234 | -0.03241 | -0.02400 | -0.01844 | +0.00208 | 1,836 |

Pooled floor, built from a vector at the median factor norm:

| variant | space | pair | median_factor_norm | pooled_floor_p05 | null_median | null_min | null_max |
|---|---|---|---|---|---|---|---|
| V_all | S1 | O1_Y1 | 20.2543 | +4.97725 | +5.39102 | +4.34768 | +6.52897 |
| V_all | S1 | O2_Y1 | 20.2543 | -0.58195 | -0.08031 | -1.41165 | +1.28820 |
| V_all | S2_k20 | O1_Y1 | 20.2543 | -0.98568 | -0.53447 | -1.56151 | +0.91052 |
| V_all | S2_k20 | O1_Y2 | 20.2543 | -0.57592 | -0.27524 | -1.04446 | +0.79928 |
| V_all | S2_k20 | O2_Y1 | 20.2543 | -0.31272 | -0.07202 | -0.68256 | +0.67399 |
| V_all | S2_k20 | O2_Y2 | 20.2543 | -0.27754 | -0.03988 | -0.67737 | +0.65304 |
| V_all | S2_k50 | O1_Y1 | 20.2543 | -0.85007 | -0.46145 | -1.34830 | +0.93857 |
| V_all | S2_k50 | O1_Y2 | 20.2543 | -0.51112 | -0.23501 | -0.95454 | +0.74806 |
| V_all | S2_k50 | O2_Y1 | 20.2543 | -0.37139 | -0.11393 | -0.83925 | +0.68014 |
| V_all | S2_k50 | O2_Y2 | 20.2543 | -0.34115 | -0.08958 | -0.78546 | +0.66883 |
| V_all | S2_k100 | O1_Y1 | 20.2543 | -0.74080 | -0.40784 | -1.22598 | +0.93370 |
| V_all | S2_k100 | O1_Y2 | 20.2543 | -0.45058 | -0.19631 | -0.89998 | +0.72730 |
| V_all | S2_k100 | O2_Y1 | 20.2543 | -0.55156 | -0.24206 | -1.10014 | +0.66755 |
| V_all | S2_k100 | O2_Y2 | 20.2543 | -0.51497 | -0.21751 | -1.04517 | +0.66663 |
| V_sig | S1 | O1_Y1 | 3.2889 | +0.15732 | +0.19966 | +0.05119 | +0.31996 |
| V_sig | S1 | O2_Y1 | 3.2889 | -0.00202 | +0.04410 | -0.10550 | +0.18455 |
| V_sig | S2_k20 | O1_Y1 | 3.2889 | -0.00776 | +0.03784 | -0.11908 | +0.16468 |
| V_sig | S2_k20 | O1_Y2 | 3.2889 | -0.02147 | +0.02117 | -0.10634 | +0.13908 |
| V_sig | S2_k20 | O2_Y1 | 3.2889 | -0.03454 | +0.00551 | -0.10527 | +0.10635 |
| V_sig | S2_k20 | O2_Y2 | 3.2889 | -0.03663 | +0.00344 | -0.10244 | +0.10219 |
| V_sig | S2_k50 | O1_Y1 | 3.2889 | -0.00950 | +0.03491 | -0.11569 | +0.15993 |
| V_sig | S2_k50 | O1_Y2 | 3.2889 | -0.02144 | +0.02020 | -0.09661 | +0.13044 |
| V_sig | S2_k50 | O2_Y1 | 3.2889 | -0.03300 | +0.00802 | -0.11460 | +0.11814 |
| V_sig | S2_k50 | O2_Y2 | 3.2889 | -0.03441 | +0.00649 | -0.11267 | +0.11533 |
| V_sig | S2_k100 | O1_Y1 | 3.2889 | -0.00937 | +0.03411 | -0.11129 | +0.15851 |
| V_sig | S2_k100 | O1_Y2 | 3.2889 | -0.02086 | +0.02020 | -0.10334 | +0.13476 |
| V_sig | S2_k100 | O2_Y1 | 3.2889 | -0.02569 | +0.01633 | -0.10490 | +0.13128 |
| V_sig | S2_k100 | O2_Y2 | 3.2889 | -0.02698 | +0.01479 | -0.10351 | +0.12803 |

### B. Secondary floor — the untargeted control rows

The same delta, using the 17 untargeted `off-target` rows as pseudo-effects, converted exactly as the factors were.

**This is reported and is never the gate.** It is built on controls that are not inert: 3 of 17 are flagged `sequence_driven` and 10 carry `de_genes` > 0 (Stage 0). A floor built from rows that themselves move the transcriptome is not a noise floor.

| variant | space | pair | n_controls | control_floor_p05 | control_min | control_median | control_max | mean_control_norm | is_gate |
|---|---|---|---|---|---|---|---|---|---|
| V_all | S1 | O1_Y1 | 17 | +1.67856 | +1.52113 | +3.44063 | +6.66310 | 17.2928 | False |
| V_all | S1 | O2_Y1 | 17 | -1.00609 | -1.11283 | -0.66421 | +0.09477 | 17.2928 | False |
| V_all | S2_k20 | O1_Y1 | 17 | -1.66758 | -1.98266 | -0.77381 | -0.06219 | 17.2928 | False |
| V_all | S2_k20 | O1_Y2 | 17 | -0.91047 | -0.95697 | -0.34069 | +0.00283 | 17.2928 | False |
| V_all | S2_k20 | O2_Y1 | 17 | -0.63028 | -0.65299 | -0.24418 | +0.36333 | 17.2928 | False |
| V_all | S2_k20 | O2_Y2 | 17 | -0.59387 | -0.63898 | -0.17814 | +0.43913 | 17.2928 | False |
| V_all | S2_k50 | O1_Y1 | 17 | -1.34013 | -1.69142 | -0.63998 | -0.06249 | 17.2928 | False |
| V_all | S2_k50 | O1_Y2 | 17 | -0.75910 | -0.84322 | -0.31063 | +0.03071 | 17.2928 | False |
| V_all | S2_k50 | O2_Y1 | 17 | -0.80889 | -0.82317 | -0.37091 | +0.39407 | 17.2928 | False |
| V_all | S2_k50 | O2_Y2 | 17 | -0.79715 | -0.81261 | -0.29745 | +0.48700 | 17.2928 | False |
| V_all | S2_k100 | O1_Y1 | 17 | -1.09925 | -1.37771 | -0.52451 | +0.08679 | 17.2928 | False |
| V_all | S2_k100 | O1_Y2 | 17 | -0.58097 | -0.77055 | -0.20781 | +0.11973 | 17.2928 | False |
| V_all | S2_k100 | O2_Y1 | 17 | -0.98210 | -1.01248 | -0.53684 | +0.17495 | 17.2928 | False |
| V_all | S2_k100 | O2_Y2 | 17 | -0.95118 | -0.99938 | -0.49203 | +0.25968 | 17.2928 | False |
| V_sig | S1 | O1_Y1 | 17 | +0.00000 | +0.00000 | +0.02841 | +2.21314 | 2.6506 | False |
| V_sig | S1 | O2_Y1 | 17 | -0.25736 | -0.28532 | +0.00000 | +0.02646 | 2.6506 | False |
| V_sig | S2_k20 | O1_Y1 | 17 | -0.34574 | -0.73164 | -0.00285 | +0.00905 | 2.6506 | False |
| V_sig | S2_k20 | O1_Y2 | 17 | -0.28360 | -0.53171 | -0.00322 | +0.00217 | 2.6506 | False |
| V_sig | S2_k20 | O2_Y1 | 17 | -0.22432 | -0.27723 | +0.00000 | +0.00698 | 2.6506 | False |
| V_sig | S2_k20 | O2_Y2 | 17 | -0.21845 | -0.25092 | +0.00000 | +0.00837 | 2.6506 | False |
| V_sig | S2_k50 | O1_Y1 | 17 | -0.28263 | -0.63977 | -0.00554 | +0.01356 | 2.6506 | False |
| V_sig | S2_k50 | O1_Y2 | 17 | -0.23115 | -0.44983 | -0.00236 | +0.00203 | 2.6506 | False |
| V_sig | S2_k50 | O2_Y1 | 17 | -0.20743 | -0.27620 | -0.00302 | +0.01092 | 2.6506 | False |
| V_sig | S2_k50 | O2_Y2 | 17 | -0.20409 | -0.25566 | -0.00147 | +0.01191 | 2.6506 | False |
| V_sig | S2_k100 | O1_Y1 | 17 | -0.25608 | -0.58508 | -0.00398 | +0.01269 | 2.6506 | False |
| V_sig | S2_k100 | O1_Y2 | 17 | -0.21053 | -0.40740 | -0.00199 | +0.00143 | 2.6506 | False |
| V_sig | S2_k100 | O2_Y1 | 17 | -0.24879 | -0.32631 | -0.00503 | +0.01112 | 2.6506 | False |
| V_sig | S2_k100 | O2_Y2 | 17 | -0.24171 | -0.30655 | -0.00474 | +0.01198 | 2.6506 | False |

### C. Dose-response and the minimum detectable approach (MDA)

The true young-minus-old direction is added to the aged origin at f = 0.0, 0.05, 0.1, 0.25, 0.5, 1.0, scaled so that f = 1.0 has the median factor norm. MDA is the smallest f whose delta falls below the **primary** (pooled) floor. `support = full_v` is the direction as it is; `support = measured_coords_only` is the same direction restricted to the coordinates a perturbation in this file could actually move, reported alongside because a real effect cannot move the unmeasured genes.

| variant | space | pair | support | median_factor_norm | pooled_floor_p05 | MDA | any_f_clears |
|---|---|---|---|---|---|---|---|
| V_all | S1 | O1_Y1 | full_v | 20.2543 | +4.97725 | 0.050 | True |
| V_all | S1 | O2_Y1 | full_v | 20.2543 | -0.58195 | 0.250 | True |
| V_all | S2_k20 | O1_Y1 | full_v | 20.2543 | -0.98568 | 0.100 | True |
| V_all | S2_k20 | O1_Y2 | full_v | 20.2543 | -0.57592 | 0.050 | True |
| V_all | S2_k20 | O2_Y1 | full_v | 20.2543 | -0.31272 | 0.100 | True |
| V_all | S2_k20 | O2_Y2 | full_v | 20.2543 | -0.27754 | 0.100 | True |
| V_all | S2_k50 | O1_Y1 | full_v | 20.2543 | -0.85007 | 0.050 | True |
| V_all | S2_k50 | O1_Y2 | full_v | 20.2543 | -0.51112 | 0.050 | True |
| V_all | S2_k50 | O2_Y1 | full_v | 20.2543 | -0.37139 | 0.100 | True |
| V_all | S2_k50 | O2_Y2 | full_v | 20.2543 | -0.34115 | 0.100 | True |
| V_all | S2_k100 | O1_Y1 | full_v | 20.2543 | -0.74080 | 0.050 | True |
| V_all | S2_k100 | O1_Y2 | full_v | 20.2543 | -0.45058 | 0.050 | True |
| V_all | S2_k100 | O2_Y1 | full_v | 20.2543 | -0.55156 | 0.250 | True |
| V_all | S2_k100 | O2_Y2 | full_v | 20.2543 | -0.51497 | 0.250 | True |
| V_all | S1 | O1_Y1 | measured_coords_only | 20.2543 | +4.97725 | 0.050 | True |
| V_all | S1 | O2_Y1 | measured_coords_only | 20.2543 | -0.58195 | 0.500 | True |
| V_all | S2_k20 | O1_Y1 | measured_coords_only | 20.2543 | -0.98568 | 0.250 | True |
| V_all | S2_k20 | O1_Y2 | measured_coords_only | 20.2543 | -0.57592 | 0.100 | True |
| V_all | S2_k20 | O2_Y1 | measured_coords_only | 20.2543 | -0.31272 | 0.100 | True |
| V_all | S2_k20 | O2_Y2 | measured_coords_only | 20.2543 | -0.27754 | 0.100 | True |
| V_all | S2_k50 | O1_Y1 | measured_coords_only | 20.2543 | -0.85007 | 0.100 | True |
| V_all | S2_k50 | O1_Y2 | measured_coords_only | 20.2543 | -0.51112 | 0.100 | True |
| V_all | S2_k50 | O2_Y1 | measured_coords_only | 20.2543 | -0.37139 | 0.250 | True |
| V_all | S2_k50 | O2_Y2 | measured_coords_only | 20.2543 | -0.34115 | 0.100 | True |
| V_all | S2_k100 | O1_Y1 | measured_coords_only | 20.2543 | -0.74080 | 0.100 | True |
| V_all | S2_k100 | O1_Y2 | measured_coords_only | 20.2543 | -0.45058 | 0.100 | True |
| V_all | S2_k100 | O2_Y1 | measured_coords_only | 20.2543 | -0.55156 | 0.250 | True |
| V_all | S2_k100 | O2_Y2 | measured_coords_only | 20.2543 | -0.51497 | 0.250 | True |
| V_sig | S1 | O1_Y1 | full_v | 3.2889 | +0.15732 | 0.050 | True |
| V_sig | S1 | O2_Y1 | full_v | 3.2889 | -0.00202 | 0.050 | True |
| V_sig | S2_k20 | O1_Y1 | full_v | 3.2889 | -0.00776 | 0.050 | True |
| V_sig | S2_k20 | O1_Y2 | full_v | 3.2889 | -0.02147 | 0.050 | True |
| V_sig | S2_k20 | O2_Y1 | full_v | 3.2889 | -0.03454 | 0.050 | True |
| V_sig | S2_k20 | O2_Y2 | full_v | 3.2889 | -0.03663 | 0.050 | True |
| V_sig | S2_k50 | O1_Y1 | full_v | 3.2889 | -0.00950 | 0.050 | True |
| V_sig | S2_k50 | O1_Y2 | full_v | 3.2889 | -0.02144 | 0.050 | True |
| V_sig | S2_k50 | O2_Y1 | full_v | 3.2889 | -0.03300 | 0.050 | True |
| V_sig | S2_k50 | O2_Y2 | full_v | 3.2889 | -0.03441 | 0.050 | True |
| V_sig | S2_k100 | O1_Y1 | full_v | 3.2889 | -0.00937 | 0.050 | True |
| V_sig | S2_k100 | O1_Y2 | full_v | 3.2889 | -0.02086 | 0.050 | True |
| V_sig | S2_k100 | O2_Y1 | full_v | 3.2889 | -0.02569 | 0.050 | True |
| V_sig | S2_k100 | O2_Y2 | full_v | 3.2889 | -0.02698 | 0.050 | True |
| V_sig | S1 | O1_Y1 | measured_coords_only | 3.2889 | +0.15732 | 0.050 | True |
| V_sig | S1 | O2_Y1 | measured_coords_only | 3.2889 | -0.00202 | 0.050 | True |
| V_sig | S2_k20 | O1_Y1 | measured_coords_only | 3.2889 | -0.00776 | 0.050 | True |
| V_sig | S2_k20 | O1_Y2 | measured_coords_only | 3.2889 | -0.02147 | 0.050 | True |
| V_sig | S2_k20 | O2_Y1 | measured_coords_only | 3.2889 | -0.03454 | 0.100 | True |
| V_sig | S2_k20 | O2_Y2 | measured_coords_only | 3.2889 | -0.03663 | 0.100 | True |
| V_sig | S2_k50 | O1_Y1 | measured_coords_only | 3.2889 | -0.00950 | 0.050 | True |
| V_sig | S2_k50 | O1_Y2 | measured_coords_only | 3.2889 | -0.02144 | 0.050 | True |
| V_sig | S2_k50 | O2_Y1 | measured_coords_only | 3.2889 | -0.03300 | 0.100 | True |
| V_sig | S2_k50 | O2_Y2 | measured_coords_only | 3.2889 | -0.03441 | 0.100 | True |
| V_sig | S2_k100 | O1_Y1 | measured_coords_only | 3.2889 | -0.00937 | 0.050 | True |
| V_sig | S2_k100 | O1_Y2 | measured_coords_only | 3.2889 | -0.02086 | 0.050 | True |
| V_sig | S2_k100 | O2_Y1 | measured_coords_only | 3.2889 | -0.02569 | 0.050 | True |
| V_sig | S2_k100 | O2_Y2 | measured_coords_only | 3.2889 | -0.02698 | 0.050 | True |

Full dose-response deltas at every f are in `results/south2/stage2c_dose_response.csv`.

- `no_detection_power`: **not fired**. Some f <= 1.0 clears the primary floor, so Stage 3 runs.
- No floor is ever subtracted from a later number. It is a bar to clear.

## Stage 3 ranking

The full table is `results/south2/stage3_full_table.csv`: **51,408 rows** = 1,836 factors x 2 gene variants x 14 space/origin/target settings, with cos, frac, delta, ||d||, ||v||, the factor's own Stage 2A floor, margin over that floor, permutation p, BH q, bootstrap intervals, the Mahalanobis change and membership test for Y2, guide agreement, the identity panel, and the frozen ruler score change. Nothing is pooled across spaces, variants, origins, or targets.

- **Guide reliability**: each guide's own vector is built and the median pairwise cosine between guides of the same factor is reported per factor and per variant. `guides_disagree` fires when that median is <= 0. Those factors are reported and are never labelled a candidate.
- **Identity check**: the panel actually used is COL1A1, COL1A2, FN1, LUM, PDGFRA, PDGFRB, PRRX1, SERPINH1, VIM — **9 of the 10 pre-registered genes**. Absent: POSTN. `identity_drop` is the negated mean displacement over that panel **in frozen-z units** (the analysis's own coordinates), and `identity_loss` fires above 0.5.
- **Frozen ruler score change** is a context column only. It gates nothing.
- **Nulls**: the Stage 2A per-factor permutations, one-sided empirical p (P(perm delta <= observed delta)), then BH-FDR across factors **within each space and variant** (and within each origin/target pair, since a pair is part of the setting).
- **Bootstrap**: 200 resamples over guides (seed 20260918) for factors with >= 3 guides. Intervals that do not contain their own point estimate are flagged INVALID and are never a gate: **618** of 51,380 computed delta intervals are so flagged.

A factor is labelled `toward_young` only if delta is negative AND below its own Stage 2A floor AND q <= 0.05 AND `identity_loss` is not flagged AND `guides_disagree` is not flagged.

### Candidate counts in every space, variant, and origin/target pair

Zero is a valid and reportable answer.

| variant | space | pair | n_factors | delta_negative | below_own_floor | q_le_005 | guides_disagree | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|
| V_all | S1 | O1_Y1 | 1,836 | 0 | 582 | 407 | 299 | 0 | 0 |
| V_all | S1 | O2_Y1 | 1,836 | 918 | 477 | 280 | 299 | 0 | 231 |
| V_all | S2_k20 | O1_Y1 | 1,836 | 1,645 | 502 | 258 | 299 | 0 | 224 |
| V_all | S2_k20 | O1_Y2 | 1,836 | 1,209 | 194 | 0 | 299 | 0 | 0 |
| V_all | S2_k20 | O2_Y1 | 1,836 | 849 | 383 | 220 | 299 | 0 | 189 |
| V_all | S2_k20 | O2_Y2 | 1,836 | 775 | 375 | 221 | 299 | 0 | 190 |
| V_all | S2_k50 | O1_Y1 | 1,836 | 1,595 | 451 | 223 | 299 | 0 | 194 |
| V_all | S2_k50 | O1_Y2 | 1,836 | 1,111 | 203 | 0 | 299 | 0 | 0 |
| V_all | S2_k50 | O2_Y1 | 1,836 | 913 | 408 | 250 | 299 | 0 | 207 |
| V_all | S2_k50 | O2_Y2 | 1,836 | 840 | 397 | 238 | 299 | 0 | 194 |
| V_all | S2_k100 | O1_Y1 | 1,836 | 1,545 | 377 | 0 | 299 | 0 | 0 |
| V_all | S2_k100 | O1_Y2 | 1,836 | 957 | 158 | 0 | 299 | 0 | 0 |
| V_all | S2_k100 | O2_Y1 | 1,836 | 1,124 | 430 | 262 | 299 | 0 | 216 |
| V_all | S2_k100 | O2_Y2 | 1,836 | 1,057 | 411 | 256 | 299 | 0 | 208 |
| V_sig | S1 | O1_Y1 | 1,836 | 0 | 324 | 0 | 1,692 | 0 | 0 |
| V_sig | S1 | O2_Y1 | 1,836 | 394 | 436 | 0 | 1,692 | 0 | 0 |
| V_sig | S2_k20 | O1_Y1 | 1,836 | 398 | 380 | 0 | 1,692 | 0 | 0 |
| V_sig | S2_k20 | O1_Y2 | 1,836 | 571 | 252 | 0 | 1,692 | 0 | 0 |
| V_sig | S2_k20 | O2_Y1 | 1,836 | 903 | 266 | 0 | 1,692 | 0 | 0 |
| V_sig | S2_k20 | O2_Y2 | 1,836 | 918 | 260 | 0 | 1,692 | 0 | 0 |
| V_sig | S2_k50 | O1_Y1 | 1,836 | 411 | 348 | 0 | 1,692 | 0 | 0 |
| V_sig | S2_k50 | O1_Y2 | 1,836 | 552 | 235 | 0 | 1,692 | 0 | 0 |
| V_sig | S2_k50 | O2_Y1 | 1,836 | 972 | 341 | 0 | 1,692 | 0 | 0 |
| V_sig | S2_k50 | O2_Y2 | 1,836 | 983 | 334 | 0 | 1,692 | 0 | 0 |
| V_sig | S2_k100 | O1_Y1 | 1,836 | 396 | 308 | 0 | 1,692 | 0 | 0 |
| V_sig | S2_k100 | O1_Y2 | 1,836 | 530 | 225 | 0 | 1,692 | 0 | 0 |
| V_sig | S2_k100 | O2_Y1 | 1,836 | 786 | 335 | 0 | 1,692 | 0 | 0 |
| V_sig | S2_k100 | O2_Y2 | 1,836 | 806 | 335 | 0 | 1,692 | 0 | 0 |

Y2 Mahalanobis columns, reported and never a gate. The Y2 `delta` above is the pre-registered Euclidean delta to the Y2 region's centre; the Mahalanobis change has its own floor from the same Stage 2A permutations. Where Y2 is undefined (see the radius table) these counts are 0 because nothing could be computed, not because nothing moved.

| variant | space | pair | md_change_below_own_floor | md_q_le_005 | origin_plus_d_inside_radius |
|---|---|---|---|---|---|
| V_all | S2_k20 | O1_Y2 | 108 | 0 | 1,836 |
| V_all | S2_k20 | O2_Y2 | 295 | 0 | 0 |
| V_all | S2_k50 | O1_Y2 | 229 | 0 | 1,836 |
| V_all | S2_k50 | O2_Y2 | 554 | 376 | 0 |
| V_all | S2_k100 | O1_Y2 | 0 | 0 | 0 |
| V_all | S2_k100 | O2_Y2 | 0 | 0 | 0 |
| V_sig | S2_k20 | O1_Y2 | 186 | 0 | 1,836 |
| V_sig | S2_k20 | O2_Y2 | 219 | 0 | 0 |
| V_sig | S2_k50 | O1_Y2 | 231 | 0 | 1,836 |
| V_sig | S2_k50 | O2_Y2 | 369 | 0 | 0 |
| V_sig | S2_k100 | O1_Y2 | 0 | 0 | 0 |
| V_sig | S2_k100 | O2_Y2 | 0 | 0 | 0 |

### Top 20 by delta in every setting (exploratory candidate tables)

`coverage_low` fired in Stage 0, so these tables are **exploratory** and every factor labelled `toward_young` carries that word in its name below.

**V_all / S1 / O1_Y1** — 0 of these 20 are `toward_young`; 0 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| REPIN1 | 6 | +1.53031 | +1.44036 | -0.08995 | 0.6667 | +0.034 | False | +0.019 | False | False |
| ZBTB46 | 6 | +1.53040 | +1.33906 | -0.19135 | 1.0000 | +0.019 | False | +0.015 | False | False |
| LHX3 | 6 | +1.77375 | +1.43557 | -0.33818 | 1.0000 | +0.033 | False | +0.009 | False | False |
| RXRA | 6 | +1.80255 | +1.72471 | -0.07784 | 0.4665 | +0.040 | False | +0.016 | False | False |
| CAPN15 | 6 | +2.06517 | +1.75257 | -0.31260 | 1.0000 | +0.001 | False | +0.035 | False | False |
| MAFG | 6 | +2.08524 | +2.50288 | +0.41764 | 0.0270 | +0.039 | False | +0.015 | False | False |
| ZIC1 | 6 | +2.15891 | +1.95406 | -0.20485 | 1.0000 | +0.035 | False | +0.085 | False | False |
| ZNF770 | 6 | +2.17700 | +2.01411 | -0.16289 | 1.0000 | -0.002 | True | +0.023 | False | False |
| MIER2 | 6 | +2.31056 | +2.33189 | +0.02133 | 0.1216 | +0.015 | False | -0.029 | False | False |
| HEYL | 6 | +2.35534 | +2.33584 | -0.01950 | 0.2675 | +0.021 | False | -0.020 | False | False |
| STAT5B | 6 | +2.36037 | +2.35677 | -0.00360 | 0.1830 | +0.006 | False | +0.052 | False | False |
| ZBTB8B | 6 | +2.37583 | +2.22815 | -0.14768 | 0.9246 | +0.002 | False | +0.014 | False | False |
| GZF1 | 6 | +2.38620 | +2.25738 | -0.12882 | 0.8597 | +0.020 | False | +0.011 | False | False |
| ZNF211 | 6 | +2.44489 | +2.60811 | +0.16322 | 0.0270 | +0.001 | False | +0.030 | False | False |
| LZTR1 | 6 | +2.46532 | +2.27009 | -0.19524 | 1.0000 | +0.015 | False | +0.054 | False | False |
| MAPK8IP1 | 6 | +2.51469 | +2.35687 | -0.15782 | 0.8597 | +0.008 | False | +0.066 | False | False |
| ZMIZ2 | 6 | +2.53936 | +2.37899 | -0.16037 | 0.9340 | +0.006 | False | -0.006 | False | False |
| TONSL | 6 | +2.55780 | +2.59201 | +0.03420 | 0.0614 | -0.007 | True | -0.000 | False | False |
| ZNF496 | 6 | +2.56473 | +2.40373 | -0.16100 | 1.0000 | +0.004 | False | -0.044 | False | False |
| TADA2B | 6 | +2.56560 | +2.73849 | +0.17289 | 0.0270 | +0.020 | False | +0.046 | False | False |

**V_all / S1 / O2_Y1** — 18 of these 20 are `toward_young`; 231 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -2.65651 | +0.61099 | +3.26750 | 0.0326 | +0.324 | False | -0.036 | False | True |
| PRDM16 (exploratory) | 6 | -1.81109 | -1.08076 | +0.73033 | 0.0326 | +0.068 | False | +0.012 | False | True |
| XRCC5 (exploratory) | 6 | -1.73822 | -0.66334 | +1.07488 | 0.0326 | +0.022 | False | +0.221 | False | True |
| ZNF208 (exploratory) | 4 | -1.67677 | +0.11962 | +1.79639 | 0.0326 | +0.019 | False | -0.109 | False | True |
| PLAG1 (exploratory) | 6 | -1.57243 | -0.42505 | +1.14738 | 0.0326 | +0.057 | False | -0.160 | False | True |
| HAND1 (exploratory) | 6 | -1.51582 | -0.23806 | +1.27775 | 0.0326 | +0.108 | False | +0.054 | False | True |
| HDAC1 (exploratory) | 6 | -1.49956 | -0.67432 | +0.82524 | 0.0326 | +0.019 | False | +0.059 | False | True |
| KLF5 (exploratory) | 6 | -1.49952 | -0.80814 | +0.69138 | 0.0326 | +0.166 | False | +0.266 | False | True |
| SPDEF (exploratory) | 6 | -1.40700 | -0.52702 | +0.87998 | 0.0326 | +0.030 | False | +0.107 | False | True |
| ZNF470 (exploratory) | 6 | -1.32587 | -0.51740 | +0.80846 | 0.0326 | +0.006 | False | +0.054 | False | True |
| ZNF615 (exploratory) | 6 | -1.31116 | -0.74120 | +0.56997 | 0.0326 | +0.026 | False | +0.124 | False | True |
| ETV1 (exploratory) | 6 | -1.31019 | -0.95298 | +0.35722 | 0.0326 | +0.061 | False | -0.086 | False | True |
| TFAP4 (exploratory) | 6 | -1.30706 | -0.59492 | +0.71214 | 0.0326 | +0.049 | False | +0.051 | False | True |
| TSHZ1 (exploratory) | 6 | -1.29290 | -0.72486 | +0.56804 | 0.0326 | +0.055 | False | +0.095 | False | True |
| LBX1 (exploratory) | 6 | -1.25818 | -0.94646 | +0.31172 | 0.0326 | +0.044 | False | +0.142 | False | True |
| KLF2 (exploratory) | 6 | -1.25818 | -0.80803 | +0.45015 | 0.0326 | +0.044 | False | +0.090 | False | True |
| ZNF230 | 6 | -1.25248 | -1.06857 | +0.18391 | 0.0562 | +0.028 | False | +0.142 | False | False |
| ENO1 (exploratory) | 6 | -1.25080 | -0.53634 | +0.71446 | 0.0326 | +0.021 | False | -0.120 | False | True |
| ZSCAN20 (exploratory) | 6 | -1.24636 | -0.64070 | +0.60566 | 0.0326 | +0.024 | False | +0.039 | False | True |
| HHEX | 6 | -1.24307 | -1.01919 | +0.22388 | 0.0562 | +0.044 | False | +0.282 | False | False |

**V_all / S2_k20 / O1_Y1** — 19 of these 20 are `toward_young`; 224 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| NCL (exploratory) | 6 | -2.73849 | -1.47827 | +1.26022 | 0.0354 | +0.020 | False | +0.404 | False | True |
| PRDM1 (exploratory) | 6 | -2.59509 | -1.66979 | +0.92531 | 0.0354 | +0.099 | False | +0.414 | False | True |
| DMRT1 (exploratory) | 6 | -2.40814 | -1.60754 | +0.80060 | 0.0354 | +0.128 | False | +0.328 | False | True |
| HOXA13 (exploratory) | 6 | -2.33688 | -1.31125 | +1.02563 | 0.0354 | +0.068 | False | +0.296 | False | True |
| LBX1 (exploratory) | 6 | -2.20883 | -1.53253 | +0.67630 | 0.0354 | +0.044 | False | +0.142 | False | True |
| PLAG1 (exploratory) | 6 | -2.14772 | -1.51380 | +0.63391 | 0.0354 | +0.057 | False | -0.160 | False | True |
| CEBPB (exploratory) | 6 | -2.14483 | -1.45686 | +0.68798 | 0.0354 | +0.091 | False | -0.233 | False | True |
| NFX1 (exploratory) | 6 | -2.10116 | -1.46185 | +0.63931 | 0.0354 | +0.022 | False | +0.180 | False | True |
| ARID2 | 6 | -2.06375 | -1.99781 | +0.06594 | 0.1284 | +0.045 | False | +0.288 | False | False |
| ZNF850 (exploratory) | 6 | -2.02629 | -1.62104 | +0.40524 | 0.0354 | +0.029 | False | +0.260 | False | True |
| PSMD14 (exploratory) | 6 | -1.99307 | -1.34070 | +0.65238 | 0.0354 | +0.031 | False | +0.125 | False | True |
| CEBPD (exploratory) | 6 | -1.99152 | -1.14880 | +0.84271 | 0.0354 | +0.043 | False | -0.064 | False | True |
| HOXB13 (exploratory) | 6 | -1.95009 | -1.20205 | +0.74804 | 0.0354 | +0.034 | False | +0.084 | False | True |
| ARID1A (exploratory) | 6 | -1.94249 | -1.22199 | +0.72050 | 0.0354 | +0.059 | False | +0.077 | False | True |
| IRF4 (exploratory) | 6 | -1.92825 | -1.53168 | +0.39657 | 0.0354 | +0.180 | False | +0.379 | False | True |
| ZBTB7A (exploratory) | 6 | -1.91444 | -1.53814 | +0.37630 | 0.0354 | +0.127 | False | -0.047 | False | True |
| ZNF544 (exploratory) | 6 | -1.90356 | -1.51098 | +0.39259 | 0.0354 | +0.031 | False | +0.103 | False | True |
| FOXJ3 (exploratory) | 6 | -1.90299 | -1.65064 | +0.25234 | 0.0354 | +0.006 | False | +0.058 | False | True |
| FUBP1 (exploratory) | 6 | -1.80866 | -1.33038 | +0.47829 | 0.0354 | +0.032 | False | +0.183 | False | True |
| ZFP64 (exploratory) | 6 | -1.77533 | -1.18509 | +0.59024 | 0.0354 | +0.045 | False | -0.065 | False | True |

**V_all / S2_k20 / O1_Y2** — 0 of these 20 are `toward_young`; 0 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 | 6 | -2.85365 | -0.70726 | +2.14639 | 0.1324 | +0.324 | False | -0.036 | False | False |
| NCL | 6 | -1.69696 | -0.93170 | +0.76526 | 0.1324 | +0.020 | False | +0.404 | False | False |
| PRDM1 | 6 | -1.52551 | -1.01716 | +0.50835 | 0.1324 | +0.099 | False | +0.414 | False | False |
| PIAS1 | 6 | -1.35891 | -0.47203 | +0.88688 | 0.1324 | +0.008 | False | +0.059 | False | False |
| ZNF850 | 6 | -1.28385 | -0.97550 | +0.30835 | 0.1324 | +0.029 | False | +0.260 | False | False |
| XRCC5 | 6 | -1.09735 | -0.91264 | +0.18472 | 0.1324 | +0.022 | False | +0.221 | False | False |
| OSR2 | 6 | -1.03599 | -0.52472 | +0.51127 | 0.1324 | +0.106 | False | -0.223 | False | False |
| TBX20 | 6 | -1.02785 | -0.51881 | +0.50904 | 0.1324 | +0.071 | False | +0.156 | False | False |
| ZNF275 | 6 | -1.02304 | -0.55644 | +0.46660 | 0.1324 | +0.015 | False | +0.073 | False | False |
| MEOX1 | 6 | -1.01666 | -0.57814 | +0.43851 | 0.1324 | -0.001 | True | -0.019 | False | False |
| ZNF286B | 6 | -1.00927 | -0.63915 | +0.37012 | 0.1324 | +0.006 | False | +0.158 | False | False |
| ZNF322 | 6 | -0.98257 | -0.44680 | +0.53577 | 0.1324 | +0.011 | False | +0.096 | False | False |
| HOXD1 | 6 | -0.97273 | -0.69620 | +0.27653 | 0.1324 | +0.029 | False | -0.053 | False | False |
| ZNF420 | 6 | -0.96477 | -0.60826 | +0.35651 | 0.2076 | +0.015 | False | +0.030 | False | False |
| ZNF845 | 5 | -0.96044 | -0.61885 | +0.34159 | 0.1324 | -0.002 | True | -0.172 | False | False |
| MKX | 6 | -0.95355 | -0.70161 | +0.25195 | 0.1324 | +0.063 | False | +0.032 | False | False |
| PPARGC1A | 6 | -0.94249 | -0.83088 | +0.11161 | 0.2995 | +0.042 | False | +0.224 | False | False |
| GFI1 | 6 | -0.93186 | -0.70153 | +0.23033 | 0.1324 | +0.231 | False | +0.058 | False | False |
| RUNX1T1 | 6 | -0.91839 | -0.94401 | -0.02562 | 0.6115 | +0.043 | False | +0.260 | False | False |
| ZFP64 | 6 | -0.90023 | -0.73200 | +0.16823 | 0.2076 | +0.045 | False | -0.065 | False | False |

**V_all / S2_k20 / O2_Y1** — 17 of these 20 are `toward_young`; 189 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -3.66590 | -0.49669 | +3.16921 | 0.0415 | +0.324 | False | -0.036 | False | True |
| TOX3 | 6 | -2.16853 | -1.52881 | +0.63972 | 0.1210 | +0.014 | False | +0.155 | False | False |
| ZEB2 (exploratory) | 6 | -1.71618 | -0.40652 | +1.30966 | 0.0415 | +0.067 | False | -0.313 | False | True |
| PIAS1 (exploratory) | 6 | -1.40572 | -0.37600 | +1.02971 | 0.0415 | +0.008 | False | +0.059 | False | True |
| MEOX1 | 6 | -1.37588 | -0.41152 | +0.96436 | 0.0415 | -0.001 | True | -0.019 | False | False |
| ZNF845 | 5 | -1.35292 | -0.39387 | +0.95905 | 0.0415 | -0.002 | True | -0.172 | False | False |
| OSR2 (exploratory) | 6 | -1.33947 | -0.40775 | +0.93172 | 0.0415 | +0.106 | False | -0.223 | False | True |
| ZNF208 (exploratory) | 4 | -1.21165 | -0.34993 | +0.86172 | 0.0415 | +0.019 | False | -0.109 | False | True |
| KLF5 (exploratory) | 6 | -1.16186 | -0.45764 | +0.70422 | 0.0415 | +0.166 | False | +0.266 | False | True |
| JDP2 (exploratory) | 6 | -1.13892 | -0.22133 | +0.91759 | 0.0415 | +0.091 | False | +0.153 | False | True |
| HOXD1 (exploratory) | 6 | -1.12963 | -0.39750 | +0.73213 | 0.0415 | +0.029 | False | -0.053 | False | True |
| PRRX2 (exploratory) | 6 | -1.07871 | -0.26045 | +0.81825 | 0.0415 | +0.028 | False | +0.084 | False | True |
| ZNF470 (exploratory) | 6 | -1.06156 | -0.29807 | +0.76348 | 0.0415 | +0.006 | False | +0.054 | False | True |
| TFAP4 (exploratory) | 6 | -1.04961 | -0.24824 | +0.80136 | 0.0415 | +0.049 | False | +0.051 | False | True |
| SPDEF (exploratory) | 6 | -1.03284 | -0.28721 | +0.74563 | 0.0415 | +0.030 | False | +0.107 | False | True |
| HAND1 (exploratory) | 6 | -1.02582 | -0.33530 | +0.69052 | 0.0415 | +0.108 | False | +0.054 | False | True |
| TWIST2 (exploratory) | 6 | -1.02426 | -0.28687 | +0.73739 | 0.0415 | +0.036 | False | -0.044 | False | True |
| ZNF841 (exploratory) | 6 | -1.01587 | -0.37790 | +0.63798 | 0.0415 | +0.011 | False | -0.058 | False | True |
| DMRTB1 (exploratory) | 6 | -1.00533 | -0.29268 | +0.71266 | 0.0415 | +0.020 | False | +0.024 | False | True |
| PRDM16 (exploratory) | 6 | -1.00300 | -0.53902 | +0.46398 | 0.0415 | +0.068 | False | +0.012 | False | True |

**V_all / S2_k20 / O2_Y2** — 18 of these 20 are `toward_young`; 190 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -3.74991 | -0.45447 | +3.29544 | 0.0413 | +0.324 | False | -0.036 | False | True |
| TOX3 (exploratory) | 6 | -2.32832 | -1.43213 | +0.89618 | 0.0413 | +0.014 | False | +0.155 | False | True |
| ZEB2 (exploratory) | 6 | -1.74136 | -0.38128 | +1.36008 | 0.0413 | +0.067 | False | -0.313 | False | True |
| PIAS1 (exploratory) | 6 | -1.46475 | -0.38316 | +1.08159 | 0.0413 | +0.008 | False | +0.059 | False | True |
| MEOX1 | 6 | -1.39784 | -0.37811 | +1.01972 | 0.0413 | -0.001 | True | -0.019 | False | False |
| ZNF845 | 5 | -1.29759 | -0.34693 | +0.95066 | 0.0413 | -0.002 | True | -0.172 | False | False |
| OSR2 (exploratory) | 6 | -1.29683 | -0.38125 | +0.91558 | 0.0413 | +0.106 | False | -0.223 | False | True |
| JDP2 (exploratory) | 6 | -1.29296 | -0.23144 | +1.06152 | 0.0413 | +0.091 | False | +0.153 | False | True |
| KLF5 (exploratory) | 6 | -1.17133 | -0.38105 | +0.79027 | 0.0413 | +0.166 | False | +0.266 | False | True |
| ZNF208 (exploratory) | 4 | -1.14237 | -0.33182 | +0.81055 | 0.0413 | +0.019 | False | -0.109 | False | True |
| PRRX2 (exploratory) | 6 | -1.11114 | -0.23065 | +0.88050 | 0.0413 | +0.028 | False | +0.084 | False | True |
| TWIST2 (exploratory) | 6 | -1.07150 | -0.25925 | +0.81225 | 0.0413 | +0.036 | False | -0.044 | False | True |
| SPDEF (exploratory) | 6 | -1.06323 | -0.25305 | +0.81018 | 0.0413 | +0.030 | False | +0.107 | False | True |
| ZNF470 (exploratory) | 6 | -1.06202 | -0.26591 | +0.79611 | 0.0413 | +0.006 | False | +0.054 | False | True |
| TFAP4 (exploratory) | 6 | -1.04625 | -0.20796 | +0.83830 | 0.0413 | +0.049 | False | +0.051 | False | True |
| DMRTB1 (exploratory) | 6 | -1.03716 | -0.26171 | +0.77545 | 0.0413 | +0.020 | False | +0.024 | False | True |
| ZNF841 (exploratory) | 6 | -1.02495 | -0.33716 | +0.68779 | 0.0413 | +0.011 | False | -0.058 | False | True |
| TBX20 (exploratory) | 6 | -1.02208 | -0.25014 | +0.77194 | 0.0413 | +0.071 | False | +0.156 | False | True |
| HAND1 (exploratory) | 6 | -1.00500 | -0.29978 | +0.70522 | 0.0413 | +0.108 | False | +0.054 | False | True |
| HOXD1 (exploratory) | 6 | -1.00169 | -0.35315 | +0.64854 | 0.0413 | +0.029 | False | -0.053 | False | True |

**V_all / S2_k50 / O1_Y1** — 19 of these 20 are `toward_young`; 194 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| PRDM1 (exploratory) | 6 | -2.38834 | -1.37295 | +1.01539 | 0.0410 | +0.099 | False | +0.414 | False | True |
| NCL (exploratory) | 6 | -2.13134 | -1.29317 | +0.83817 | 0.0410 | +0.020 | False | +0.404 | False | True |
| HOXA13 (exploratory) | 6 | -2.10348 | -1.13369 | +0.96980 | 0.0410 | +0.068 | False | +0.296 | False | True |
| PLAG1 (exploratory) | 6 | -2.06176 | -1.26366 | +0.79810 | 0.0410 | +0.057 | False | -0.160 | False | True |
| DMRT1 (exploratory) | 6 | -1.83980 | -1.33040 | +0.50940 | 0.0410 | +0.128 | False | +0.328 | False | True |
| ZNF850 (exploratory) | 6 | -1.81767 | -1.35105 | +0.46663 | 0.0410 | +0.029 | False | +0.260 | False | True |
| HOXD1 (exploratory) | 6 | -1.81471 | -0.93884 | +0.87587 | 0.0410 | +0.029 | False | -0.053 | False | True |
| LBX1 (exploratory) | 6 | -1.79339 | -1.31463 | +0.47876 | 0.0410 | +0.044 | False | +0.142 | False | True |
| HOXB13 (exploratory) | 6 | -1.78387 | -1.03732 | +0.74655 | 0.0410 | +0.034 | False | +0.084 | False | True |
| ARID1A (exploratory) | 6 | -1.71644 | -1.04450 | +0.67194 | 0.0410 | +0.059 | False | +0.077 | False | True |
| CEBPB (exploratory) | 6 | -1.69760 | -1.25123 | +0.44636 | 0.0410 | +0.091 | False | -0.233 | False | True |
| CEBPD (exploratory) | 6 | -1.65247 | -1.02399 | +0.62848 | 0.0410 | +0.043 | False | -0.064 | False | True |
| SMAD1 (exploratory) | 6 | -1.62847 | -1.01905 | +0.60943 | 0.0410 | +0.098 | False | +0.366 | False | True |
| TSHZ1 (exploratory) | 6 | -1.62371 | -1.05808 | +0.56563 | 0.0410 | +0.055 | False | +0.095 | False | True |
| RUNX1T1 (exploratory) | 6 | -1.60719 | -1.37062 | +0.23657 | 0.0410 | +0.043 | False | +0.260 | False | True |
| SMAD9 | 6 | -1.59781 | -1.51140 | +0.08641 | 0.1571 | +0.011 | False | +0.252 | False | False |
| ZFP64 (exploratory) | 6 | -1.59269 | -0.99077 | +0.60192 | 0.0410 | +0.045 | False | -0.065 | False | True |
| PSMD14 (exploratory) | 6 | -1.58880 | -1.14612 | +0.44268 | 0.0410 | +0.031 | False | +0.125 | False | True |
| NFX1 (exploratory) | 6 | -1.53231 | -1.29020 | +0.24211 | 0.0410 | +0.022 | False | +0.180 | False | True |
| FUBP1 (exploratory) | 6 | -1.53019 | -1.14249 | +0.38770 | 0.0410 | +0.032 | False | +0.183 | False | True |

**V_all / S2_k50 / O1_Y2** — 0 of these 20 are `toward_young`; 0 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 | 6 | -2.17852 | -0.62375 | +1.55476 | 0.1202 | +0.324 | False | -0.036 | False | False |
| PRDM1 | 6 | -1.31278 | -0.82140 | +0.49138 | 0.1809 | +0.099 | False | +0.414 | False | False |
| NCL | 6 | -1.24995 | -0.91524 | +0.33471 | 0.1202 | +0.020 | False | +0.404 | False | False |
| PIAS1 | 6 | -1.17333 | -0.42720 | +0.74612 | 0.1202 | +0.008 | False | +0.059 | False | False |
| ZNF850 | 6 | -1.11269 | -0.83404 | +0.27865 | 0.1809 | +0.029 | False | +0.260 | False | False |
| TBX20 | 6 | -1.04806 | -0.48107 | +0.56699 | 0.1202 | +0.071 | False | +0.156 | False | False |
| HOXD1 | 6 | -0.97492 | -0.62375 | +0.35117 | 0.1202 | +0.029 | False | -0.053 | False | False |
| NFIX | 6 | -0.96107 | -0.62522 | +0.33585 | 0.1202 | +0.004 | False | +0.014 | False | False |
| MEOX1 | 6 | -0.95621 | -0.56269 | +0.39352 | 0.1202 | -0.001 | True | -0.019 | False | False |
| ZNF322 | 6 | -0.94021 | -0.40203 | +0.53818 | 0.1202 | +0.011 | False | +0.096 | False | False |
| ZNF845 | 5 | -0.88267 | -0.59139 | +0.29129 | 0.1202 | -0.002 | True | -0.172 | False | False |
| ZNF286B | 6 | -0.86885 | -0.60036 | +0.26849 | 0.1202 | +0.006 | False | +0.158 | False | False |
| ZNF281 | 6 | -0.86818 | -0.40418 | +0.46400 | 0.1202 | +0.028 | False | -0.015 | False | False |
| RBMS1 | 6 | -0.86536 | -0.46611 | +0.39925 | 0.1202 | +0.016 | False | +0.281 | False | False |
| MKX | 6 | -0.86508 | -0.58843 | +0.27666 | 0.1202 | +0.063 | False | +0.032 | False | False |
| ZNF420 | 6 | -0.85666 | -0.57428 | +0.28238 | 0.1202 | +0.015 | False | +0.030 | False | False |
| ZNF814 | 6 | -0.85037 | -0.36056 | +0.48982 | 0.1202 | +0.018 | False | +0.043 | False | False |
| ILF3 | 6 | -0.84569 | -0.54245 | +0.30324 | 0.1202 | +0.020 | False | -0.030 | False | False |
| XRCC5 | 6 | -0.83388 | -0.73451 | +0.09937 | 0.2322 | +0.022 | False | +0.221 | False | False |
| MBD3 | 6 | -0.82689 | -0.53587 | +0.29103 | 0.1202 | +0.005 | False | -0.036 | False | False |

**V_all / S2_k50 / O2_Y1** — 18 of these 20 are `toward_young`; 207 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -3.72251 | -0.50230 | +3.22021 | 0.0365 | +0.324 | False | -0.036 | False | True |
| TOX3 (exploratory) | 6 | -3.15294 | -1.64296 | +1.50998 | 0.0365 | +0.014 | False | +0.155 | False | True |
| PRDM1 (exploratory) | 6 | -1.73425 | -0.78900 | +0.94525 | 0.0365 | +0.099 | False | +0.414 | False | True |
| ZNF208 (exploratory) | 4 | -1.72839 | -0.45184 | +1.27655 | 0.0365 | +0.019 | False | -0.109 | False | True |
| PIAS1 (exploratory) | 6 | -1.56676 | -0.36963 | +1.19714 | 0.0365 | +0.008 | False | +0.059 | False | True |
| ZNF845 | 5 | -1.49186 | -0.44184 | +1.05002 | 0.0365 | -0.002 | True | -0.172 | False | False |
| MEOX1 | 6 | -1.41899 | -0.46435 | +0.95464 | 0.0365 | -0.001 | True | -0.019 | False | False |
| HOXD1 (exploratory) | 6 | -1.35789 | -0.49013 | +0.86775 | 0.0365 | +0.029 | False | -0.053 | False | True |
| HAND1 (exploratory) | 6 | -1.32447 | -0.35546 | +0.96901 | 0.0365 | +0.108 | False | +0.054 | False | True |
| KLF5 (exploratory) | 6 | -1.31821 | -0.57800 | +0.74021 | 0.0365 | +0.166 | False | +0.266 | False | True |
| XRCC5 (exploratory) | 6 | -1.28777 | -0.68729 | +0.60047 | 0.0365 | +0.022 | False | +0.221 | False | True |
| NCOA1 (exploratory) | 6 | -1.25190 | -0.36710 | +0.88481 | 0.0365 | +0.012 | False | -0.072 | False | True |
| ZNF615 (exploratory) | 6 | -1.24714 | -0.50465 | +0.74249 | 0.0365 | +0.026 | False | +0.124 | False | True |
| ZEB2 (exploratory) | 6 | -1.21148 | -0.44895 | +0.76253 | 0.0365 | +0.067 | False | -0.313 | False | True |
| PLAG1 (exploratory) | 6 | -1.18013 | -0.62866 | +0.55147 | 0.0365 | +0.057 | False | -0.160 | False | True |
| ZNF322 (exploratory) | 6 | -1.16868 | -0.34669 | +0.82199 | 0.0365 | +0.011 | False | +0.096 | False | True |
| ZNF850 (exploratory) | 6 | -1.13320 | -0.67057 | +0.46262 | 0.0365 | +0.029 | False | +0.260 | False | True |
| CDCA7L (exploratory) | 6 | -1.12732 | -0.70894 | +0.41838 | 0.0365 | +0.030 | False | -0.111 | False | True |
| SHOX (exploratory) | 4 | -1.12564 | -0.46654 | +0.65910 | 0.0365 | +0.074 | False | -0.021 | False | True |
| TCF3 (exploratory) | 6 | -1.11758 | -0.37660 | +0.74099 | 0.0365 | +0.010 | False | +0.156 | False | True |

**V_all / S2_k50 / O2_Y2** — 18 of these 20 are `toward_young`; 194 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -3.79218 | -0.48486 | +3.30732 | 0.0384 | +0.324 | False | -0.036 | False | True |
| TOX3 (exploratory) | 6 | -3.24535 | -1.62913 | +1.61622 | 0.0384 | +0.014 | False | +0.155 | False | True |
| ZNF208 (exploratory) | 4 | -1.66424 | -0.43343 | +1.23081 | 0.0384 | +0.019 | False | -0.109 | False | True |
| PIAS1 (exploratory) | 6 | -1.60232 | -0.34562 | +1.25671 | 0.0384 | +0.008 | False | +0.059 | False | True |
| PRDM1 (exploratory) | 6 | -1.59874 | -0.71836 | +0.88038 | 0.0384 | +0.099 | False | +0.414 | False | True |
| ZNF845 | 5 | -1.45213 | -0.41485 | +1.03728 | 0.0384 | -0.002 | True | -0.172 | False | False |
| MEOX1 | 6 | -1.43961 | -0.44720 | +0.99241 | 0.0384 | -0.001 | True | -0.019 | False | False |
| KLF5 (exploratory) | 6 | -1.33401 | -0.51315 | +0.82086 | 0.0384 | +0.166 | False | +0.266 | False | True |
| HAND1 (exploratory) | 6 | -1.30501 | -0.34231 | +0.96270 | 0.0384 | +0.108 | False | +0.054 | False | True |
| ZEB2 (exploratory) | 6 | -1.25240 | -0.42839 | +0.82400 | 0.0384 | +0.067 | False | -0.313 | False | True |
| HOXD1 (exploratory) | 6 | -1.24988 | -0.46056 | +0.78932 | 0.0384 | +0.029 | False | -0.053 | False | True |
| XRCC5 (exploratory) | 6 | -1.24069 | -0.62218 | +0.61851 | 0.0384 | +0.022 | False | +0.221 | False | True |
| NCOA1 (exploratory) | 6 | -1.19030 | -0.36892 | +0.82138 | 0.0384 | +0.012 | False | -0.072 | False | True |
| ZNF322 (exploratory) | 6 | -1.16571 | -0.32610 | +0.83961 | 0.0384 | +0.011 | False | +0.096 | False | True |
| ZNF615 (exploratory) | 6 | -1.16296 | -0.46080 | +0.70215 | 0.0384 | +0.026 | False | +0.124 | False | True |
| ZNF470 (exploratory) | 6 | -1.09580 | -0.31568 | +0.78012 | 0.0384 | +0.006 | False | +0.054 | False | True |
| TFAP4 (exploratory) | 6 | -1.09455 | -0.28870 | +0.80585 | 0.0384 | +0.049 | False | +0.051 | False | True |
| SHOX (exploratory) | 4 | -1.07585 | -0.43073 | +0.64512 | 0.0384 | +0.074 | False | -0.021 | False | True |
| TCF3 (exploratory) | 6 | -1.06378 | -0.34692 | +0.71686 | 0.0384 | +0.010 | False | +0.156 | False | True |
| ZNF589 (exploratory) | 6 | -1.06377 | -0.41891 | +0.64486 | 0.0384 | +0.002 | False | +0.139 | False | True |

**V_all / S2_k100 / O1_Y1** — 0 of these 20 are `toward_young`; 0 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| PRDM1 | 6 | -1.96796 | -1.11915 | +0.84882 | 0.0534 | +0.099 | False | +0.414 | False | False |
| HOXA13 | 6 | -1.87016 | -0.97827 | +0.89188 | 0.0534 | +0.068 | False | +0.296 | False | False |
| NCL | 6 | -1.84777 | -1.10834 | +0.73942 | 0.0534 | +0.020 | False | +0.404 | False | False |
| PLAG1 | 6 | -1.68989 | -1.03581 | +0.65409 | 0.0534 | +0.057 | False | -0.160 | False | False |
| HOXD1 | 6 | -1.58945 | -0.82860 | +0.76085 | 0.0534 | +0.029 | False | -0.053 | False | False |
| HOXB13 | 6 | -1.51522 | -0.91627 | +0.59895 | 0.0534 | +0.034 | False | +0.084 | False | False |
| DMRT1 | 6 | -1.48518 | -1.05758 | +0.42760 | 0.0534 | +0.128 | False | +0.328 | False | False |
| ARID1A | 6 | -1.45248 | -0.90260 | +0.54988 | 0.0534 | +0.059 | False | +0.077 | False | False |
| ZNF850 | 6 | -1.43945 | -1.10236 | +0.33709 | 0.0858 | +0.029 | False | +0.260 | False | False |
| SMAD1 | 6 | -1.42093 | -0.89174 | +0.52919 | 0.0534 | +0.098 | False | +0.366 | False | False |
| BCL11B | 6 | -1.38915 | -0.56987 | +0.81928 | 0.0534 | +0.030 | False | -0.056 | False | False |
| PSMD14 | 6 | -1.35540 | -0.96915 | +0.38626 | 0.0534 | +0.031 | False | +0.125 | False | False |
| ZNF583 | 6 | -1.35528 | -0.83265 | +0.52262 | 0.0534 | +0.021 | False | +0.006 | False | False |
| TSHZ1 | 6 | -1.34612 | -0.91423 | +0.43189 | 0.0534 | +0.055 | False | +0.095 | False | False |
| CEBPD | 6 | -1.33768 | -0.89643 | +0.44125 | 0.0534 | +0.043 | False | -0.064 | False | False |
| ZFP64 | 6 | -1.33561 | -0.86153 | +0.47409 | 0.0534 | +0.045 | False | -0.065 | False | False |
| LBX1 | 6 | -1.33441 | -1.05406 | +0.28035 | 0.0534 | +0.044 | False | +0.142 | False | False |
| ZNF705B | 4 | -1.28701 | -0.71937 | +0.56764 | 0.0534 | +0.005 | False | -0.015 | False | False |
| NFIA | 6 | -1.28114 | -0.83418 | +0.44696 | 0.0534 | +0.053 | False | +0.188 | False | False |
| CEBPB | 6 | -1.26746 | -1.00654 | +0.26092 | 0.0858 | +0.091 | False | -0.233 | False | False |

**V_all / S2_k100 / O1_Y2** — 0 of these 20 are `toward_young`; 0 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 | 6 | -1.55534 | -0.59003 | +0.96531 | 0.1631 | +0.324 | False | -0.036 | False | False |
| NCL | 6 | -1.02231 | -0.75541 | +0.26690 | 0.1631 | +0.020 | False | +0.404 | False | False |
| PIAS1 | 6 | -1.00180 | -0.41543 | +0.58638 | 0.1631 | +0.008 | False | +0.059 | False | False |
| PRDM1 | 6 | -0.95290 | -0.61635 | +0.33655 | 0.2312 | +0.099 | False | +0.414 | False | False |
| TBX20 | 6 | -0.90707 | -0.43786 | +0.46921 | 0.1631 | +0.071 | False | +0.156 | False | False |
| ZNF322 | 6 | -0.86984 | -0.37354 | +0.49630 | 0.1631 | +0.011 | False | +0.096 | False | False |
| NFIX | 6 | -0.83850 | -0.56719 | +0.27131 | 0.1631 | +0.004 | False | +0.014 | False | False |
| MEOX1 | 6 | -0.82433 | -0.52181 | +0.30252 | 0.1631 | -0.001 | True | -0.019 | False | False |
| ZNF850 | 6 | -0.81710 | -0.65277 | +0.16433 | 0.3582 | +0.029 | False | +0.260 | False | False |
| ZNF814 | 6 | -0.80966 | -0.33917 | +0.47049 | 0.1631 | +0.018 | False | +0.043 | False | False |
| ZNF281 | 6 | -0.79291 | -0.38718 | +0.40573 | 0.1631 | +0.028 | False | -0.015 | False | False |
| RBMS1 | 6 | -0.79258 | -0.40828 | +0.38430 | 0.1631 | +0.016 | False | +0.281 | False | False |
| HOXD1 | 6 | -0.79220 | -0.56289 | +0.22932 | 0.1631 | +0.029 | False | -0.053 | False | False |
| HMG20A | 6 | -0.78659 | -0.30477 | +0.48182 | 0.1631 | +0.013 | False | +0.021 | False | False |
| ILF3 | 6 | -0.77577 | -0.46980 | +0.30598 | 0.1631 | +0.020 | False | -0.030 | False | False |
| ZNF845 | 5 | -0.77576 | -0.53663 | +0.23913 | 0.1631 | -0.002 | True | -0.172 | False | False |
| ZNF420 | 6 | -0.76648 | -0.54232 | +0.22416 | 0.3186 | +0.015 | False | +0.030 | False | False |
| ZNF805 | 6 | -0.75177 | -0.37905 | +0.37272 | 0.1631 | -0.003 | True | +0.031 | False | False |
| MAFF | 6 | -0.74021 | -0.43522 | +0.30499 | 0.1631 | +0.026 | False | +0.095 | False | False |
| ZNF275 | 6 | -0.73788 | -0.45208 | +0.28579 | 0.1631 | +0.015 | False | +0.073 | False | False |

**V_all / S2_k100 / O2_Y1** — 18 of these 20 are `toward_young`; 216 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -3.75179 | -0.63727 | +3.11452 | 0.0349 | +0.324 | False | -0.036 | False | True |
| TOX3 (exploratory) | 6 | -3.62221 | -1.98441 | +1.63781 | 0.0349 | +0.014 | False | +0.155 | False | True |
| PRDM1 (exploratory) | 6 | -2.20173 | -1.11125 | +1.09048 | 0.0349 | +0.099 | False | +0.414 | False | True |
| ZNF208 (exploratory) | 4 | -1.87202 | -0.50330 | +1.36872 | 0.0349 | +0.019 | False | -0.109 | False | True |
| PIAS1 (exploratory) | 6 | -1.72373 | -0.37062 | +1.35311 | 0.0349 | +0.008 | False | +0.059 | False | True |
| XRCC5 (exploratory) | 6 | -1.71284 | -0.99301 | +0.71983 | 0.0349 | +0.022 | False | +0.221 | False | True |
| KLF5 (exploratory) | 6 | -1.71241 | -0.85648 | +0.85593 | 0.0349 | +0.166 | False | +0.266 | False | True |
| CDCA7L (exploratory) | 6 | -1.66618 | -1.03100 | +0.63518 | 0.0349 | +0.030 | False | -0.111 | False | True |
| MEOX1 | 6 | -1.61731 | -0.53592 | +1.08138 | 0.0349 | -0.001 | True | -0.019 | False | False |
| ZNF845 | 5 | -1.61410 | -0.57856 | +1.03553 | 0.0349 | -0.002 | True | -0.172 | False | False |
| HOXD1 (exploratory) | 6 | -1.59306 | -0.67072 | +0.92234 | 0.0349 | +0.029 | False | -0.053 | False | True |
| ZNF615 (exploratory) | 6 | -1.53075 | -0.73784 | +0.79291 | 0.0349 | +0.026 | False | +0.124 | False | True |
| PRDM16 (exploratory) | 6 | -1.50178 | -1.11222 | +0.38956 | 0.0349 | +0.068 | False | +0.012 | False | True |
| ZNF850 (exploratory) | 6 | -1.48319 | -1.02031 | +0.46288 | 0.0349 | +0.029 | False | +0.260 | False | True |
| PLAG1 (exploratory) | 6 | -1.47870 | -0.90625 | +0.57244 | 0.0349 | +0.057 | False | -0.160 | False | True |
| HAND1 (exploratory) | 6 | -1.44323 | -0.49282 | +0.95042 | 0.0349 | +0.108 | False | +0.054 | False | True |
| ZEB2 (exploratory) | 6 | -1.35886 | -0.49281 | +0.86605 | 0.0349 | +0.067 | False | -0.313 | False | True |
| SHOX (exploratory) | 4 | -1.35530 | -0.66526 | +0.69005 | 0.0349 | +0.074 | False | -0.021 | False | True |
| ZIK1 (exploratory) | 6 | -1.32059 | -0.57564 | +0.74494 | 0.0349 | +0.023 | False | +0.091 | False | True |
| SULT2A1 (exploratory) | 5 | -1.30770 | -0.77474 | +0.53296 | 0.0349 | +0.008 | False | +0.046 | False | True |

**V_all / S2_k100 / O2_Y2** — 18 of these 20 are `toward_young`; 208 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -3.82159 | -0.61442 | +3.20718 | 0.0357 | +0.324 | False | -0.036 | False | True |
| TOX3 (exploratory) | 6 | -3.71650 | -1.87086 | +1.84564 | 0.0357 | +0.014 | False | +0.155 | False | True |
| PRDM1 (exploratory) | 6 | -2.06303 | -1.05116 | +1.01187 | 0.0357 | +0.099 | False | +0.414 | False | True |
| ZNF208 (exploratory) | 4 | -1.81112 | -0.47230 | +1.33882 | 0.0357 | +0.019 | False | -0.109 | False | True |
| PIAS1 (exploratory) | 6 | -1.75483 | -0.38054 | +1.37429 | 0.0357 | +0.008 | False | +0.059 | False | True |
| KLF5 (exploratory) | 6 | -1.72088 | -0.78976 | +0.93112 | 0.0357 | +0.166 | False | +0.266 | False | True |
| XRCC5 (exploratory) | 6 | -1.66605 | -0.93278 | +0.73327 | 0.0357 | +0.022 | False | +0.221 | False | True |
| MEOX1 | 6 | -1.63400 | -0.50324 | +1.13076 | 0.0357 | -0.001 | True | -0.019 | False | False |
| ZNF845 | 5 | -1.57941 | -0.55526 | +1.02415 | 0.0357 | -0.002 | True | -0.172 | False | False |
| CDCA7L (exploratory) | 6 | -1.55871 | -0.96161 | +0.59709 | 0.0357 | +0.030 | False | -0.111 | False | True |
| HOXD1 (exploratory) | 6 | -1.48696 | -0.64998 | +0.83698 | 0.0357 | +0.029 | False | -0.053 | False | True |
| ZNF615 (exploratory) | 6 | -1.44890 | -0.69627 | +0.75263 | 0.0357 | +0.026 | False | +0.124 | False | True |
| HAND1 (exploratory) | 6 | -1.42558 | -0.46292 | +0.96266 | 0.0357 | +0.108 | False | +0.054 | False | True |
| ZNF850 (exploratory) | 6 | -1.39871 | -0.94187 | +0.45683 | 0.0357 | +0.029 | False | +0.260 | False | True |
| PRDM16 (exploratory) | 6 | -1.39641 | -1.02510 | +0.37131 | 0.0357 | +0.068 | False | +0.012 | False | True |
| ZEB2 (exploratory) | 6 | -1.39343 | -0.46901 | +0.92442 | 0.0357 | +0.067 | False | -0.313 | False | True |
| PLAG1 (exploratory) | 6 | -1.32929 | -0.84307 | +0.48622 | 0.0357 | +0.057 | False | -0.160 | False | True |
| ZIK1 (exploratory) | 6 | -1.31572 | -0.54158 | +0.77414 | 0.0357 | +0.023 | False | +0.091 | False | True |
| SULT2A1 (exploratory) | 5 | -1.31079 | -0.75462 | +0.55617 | 0.0357 | +0.008 | False | +0.046 | False | True |
| SHOX (exploratory) | 4 | -1.30688 | -0.61741 | +0.68947 | 0.0357 | +0.074 | False | -0.021 | False | True |

**V_sig / S1 / O1_Y1** — 0 of these 20 are `toward_young`; 0 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| ZNF175 | 6 | +0.00492 | +0.00670 | +0.00178 | 0.2514 | +0.000 | True | -0.000 | False | False |
| TCF25 | 6 | +0.00846 | +0.00886 | +0.00040 | 0.2687 | +0.000 | True | -0.000 | False | False |
| NR5A1 | 6 | +0.01182 | +0.00763 | -0.00418 | 0.4881 | +0.000 | True | -0.000 | False | False |
| SKI | 6 | +0.01245 | +0.00513 | -0.00732 | 0.6742 | +0.000 | True | -0.000 | False | False |
| PKNOX1 | 6 | +0.01297 | +0.00396 | -0.00901 | 0.8433 | +0.000 | True | -0.000 | False | False |
| EN1 | 6 | +0.01331 | +0.00756 | -0.00575 | 0.5656 | +0.000 | True | -0.000 | False | False |
| CREBL2 | 6 | +0.01381 | +0.00836 | -0.00545 | 0.5238 | +0.000 | True | -0.000 | False | False |
| POU5F1B | 6 | +0.01463 | +0.00558 | -0.00905 | 0.7186 | +0.000 | True | -0.000 | False | False |
| ZBTB8B | 6 | +0.01539 | +0.01536 | -0.00002 | 0.3114 | +0.000 | True | -0.000 | False | False |
| ZNF341 | 6 | +0.01557 | +0.00752 | -0.00805 | 0.7161 | +0.000 | True | -0.000 | False | False |
| MESP2 | 6 | +0.01622 | +0.00713 | -0.00909 | 0.6691 | +0.000 | True | -0.000 | False | False |
| XBP1 | 6 | +0.01970 | +0.01217 | -0.00752 | 0.5753 | +0.000 | True | -0.000 | False | False |
| NFXL1 | 6 | +0.01983 | +0.00580 | -0.01402 | 0.9048 | +0.000 | True | -0.000 | False | False |
| ZNF232 | 6 | +0.01984 | +0.00928 | -0.01055 | 0.8487 | +0.000 | True | -0.012 | False | False |
| ZNF853 | 6 | +0.01988 | +0.03397 | +0.01408 | 0.1522 | +0.000 | True | -0.000 | False | False |
| ZNF146 | 6 | +0.01990 | +0.01055 | -0.00935 | 0.5510 | +0.000 | True | -0.000 | False | False |
| CCNT2 | 6 | +0.01998 | +0.01949 | -0.00049 | 0.3114 | +0.000 | True | -0.000 | False | False |
| HMGXB4 | 6 | +0.02205 | +0.01827 | -0.00377 | 0.3724 | +0.000 | True | -0.000 | False | False |
| ISX | 6 | +0.02444 | +0.02442 | -0.00002 | 0.3114 | +0.000 | True | -0.000 | False | False |
| CTBP1 | 6 | +0.02444 | +0.01870 | -0.00574 | 0.5197 | +0.000 | True | -0.000 | False | False |

**V_sig / S1 / O2_Y1** — 0 of these 20 are `toward_young`; 0 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 | 6 | -1.47611 | +0.18582 | +1.66194 | 0.0593 | +0.579 | False | -0.081 | False | False |
| BSX | 6 | -0.48456 | -0.15819 | +0.32637 | 0.0593 | +0.003 | False | +0.021 | False | False |
| KLF5 | 6 | -0.41908 | -0.03417 | +0.38491 | 0.0593 | +0.157 | False | +0.126 | False | False |
| GFI1 | 6 | -0.38395 | -0.06203 | +0.32192 | 0.0593 | +0.370 | False | +0.061 | False | False |
| HHEX | 6 | -0.31065 | -0.10147 | +0.20917 | 0.0593 | +0.002 | False | +0.053 | False | False |
| NKX3-1 | 6 | -0.23803 | -0.02031 | +0.21772 | 0.0593 | +0.205 | False | +0.061 | False | False |
| ZFP2 | 6 | -0.20805 | -0.07640 | +0.13166 | 0.0593 | +0.000 | True | +0.042 | False | False |
| ELOC | 6 | -0.19802 | -0.04453 | +0.15349 | 0.0593 | +0.000 | True | -0.020 | False | False |
| POU3F3 | 6 | -0.19064 | -0.04688 | +0.14376 | 0.0593 | +0.000 | True | -0.040 | False | False |
| GATA1 | 6 | -0.19002 | +0.03166 | +0.22168 | 0.0593 | +0.000 | True | +0.019 | False | False |
| ZNF557 | 6 | -0.18286 | -0.06913 | +0.11374 | 0.0593 | +0.000 | True | -0.020 | False | False |
| ZNF574 | 6 | -0.17789 | -0.02828 | +0.14962 | 0.0593 | +0.000 | True | +0.027 | False | False |
| TFCP2 | 6 | -0.17013 | -0.04297 | +0.12716 | 0.0593 | +0.000 | True | +0.059 | False | False |
| RBMS1 | 6 | -0.15807 | -0.04615 | +0.11192 | 0.0593 | +0.000 | True | +0.037 | False | False |
| PRRX2 | 6 | -0.15292 | -0.02902 | +0.12390 | 0.0593 | +0.000 | True | +0.034 | False | False |
| MNX1 | 6 | -0.14675 | -0.05270 | +0.09405 | 0.0593 | +0.018 | False | +0.015 | False | False |
| TBX20 | 6 | -0.14334 | -0.03205 | +0.11129 | 0.0593 | +0.012 | False | +0.077 | False | False |
| MYOG | 6 | -0.14008 | +0.01389 | +0.15397 | 0.0593 | +0.000 | True | -0.049 | False | False |
| MXI1 | 6 | -0.13846 | +0.04476 | +0.18322 | 0.0593 | +0.012 | False | -0.000 | False | False |
| CLOCK | 6 | -0.13785 | -0.03589 | +0.10196 | 0.0593 | +0.000 | True | -0.004 | False | False |

**V_sig / S2_k20 / O1_Y1** — 0 of these 20 are `toward_young`; 0 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 | 6 | -0.66365 | -0.50337 | +0.16027 | 0.1202 | +0.579 | False | -0.081 | False | False |
| BSX | 6 | -0.53870 | -0.28753 | +0.25117 | 0.0932 | +0.003 | False | +0.021 | False | False |
| GFI1 | 6 | -0.50501 | -0.18194 | +0.32307 | 0.0932 | +0.370 | False | +0.061 | False | False |
| HHEX | 6 | -0.29910 | -0.18704 | +0.11206 | 0.0932 | +0.002 | False | +0.053 | False | False |
| IRF1 | 6 | -0.22983 | +0.14856 | +0.37839 | 0.0932 | +0.537 | False | +0.025 | False | False |
| OLIG1 | 6 | -0.19257 | -0.06215 | +0.13042 | 0.0932 | +0.000 | True | +0.026 | False | False |
| PRDM16 | 6 | -0.18864 | -0.03266 | +0.15598 | 0.0932 | +0.006 | False | -0.111 | False | False |
| KLF5 | 6 | -0.18758 | -0.21301 | -0.02544 | 0.3116 | +0.157 | False | +0.126 | False | False |
| PRDM1 | 6 | -0.17827 | +0.11815 | +0.29642 | 0.0932 | +0.031 | False | +0.073 | False | False |
| ZNF763 | 6 | -0.16801 | -0.00393 | +0.16408 | 0.0932 | +0.000 | True | +0.053 | False | False |
| ETS1 | 6 | -0.16546 | -0.04319 | +0.12226 | 0.0932 | +0.000 | True | +0.041 | False | False |
| SOX11 | 5 | -0.16103 | -0.01485 | +0.14618 | 0.0932 | +0.001 | False | +0.010 | False | False |
| RBMS1 | 6 | -0.15626 | -0.07105 | +0.08521 | 0.0932 | +0.000 | True | +0.037 | False | False |
| LHX4 | 6 | -0.15404 | -0.01097 | +0.14307 | 0.0932 | +0.021 | False | -0.039 | False | False |
| MYOD1 | 6 | -0.15278 | -0.06475 | +0.08803 | 0.0932 | +0.000 | True | +0.013 | False | False |
| CDX2 | 6 | -0.15203 | -0.04098 | +0.11105 | 0.0932 | +0.215 | False | +0.041 | False | False |
| HNRNPAB | 6 | -0.13949 | +0.18254 | +0.32204 | 0.0932 | +0.000 | True | -0.000 | False | False |
| POU3F3 | 6 | -0.13350 | -0.05738 | +0.07612 | 0.0932 | +0.000 | True | -0.040 | False | False |
| NKX3-1 | 6 | -0.13251 | -0.18501 | -0.05250 | 0.4675 | +0.205 | False | +0.061 | False | False |
| SOX5 | 6 | -0.13147 | +0.01801 | +0.14948 | 0.0932 | +0.000 | True | +0.012 | False | False |

**V_sig / S2_k20 / O1_Y2** — 0 of these 20 are `toward_young`; 0 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 | 6 | -1.76899 | -0.42742 | +1.34157 | 0.1305 | +0.579 | False | -0.081 | False | False |
| BSX | 6 | -0.59579 | -0.20902 | +0.38677 | 0.1305 | +0.003 | False | +0.021 | False | False |
| NKX3-1 | 6 | -0.50802 | -0.19439 | +0.31363 | 0.1305 | +0.205 | False | +0.061 | False | False |
| KLF5 | 6 | -0.50050 | -0.20024 | +0.30027 | 0.1305 | +0.157 | False | +0.126 | False | False |
| GFI1 | 6 | -0.46155 | -0.16178 | +0.29977 | 0.1305 | +0.370 | False | +0.061 | False | False |
| HHEX | 6 | -0.38523 | -0.15156 | +0.23367 | 0.1305 | +0.002 | False | +0.053 | False | False |
| HNRNPAB | 6 | -0.36118 | -0.04677 | +0.31441 | 0.1305 | +0.000 | True | -0.000 | False | False |
| MXI1 | 6 | -0.27466 | -0.04839 | +0.22627 | 0.1305 | +0.012 | False | -0.000 | False | False |
| POU3F3 | 6 | -0.24870 | -0.05253 | +0.19617 | 0.1305 | +0.000 | True | -0.040 | False | False |
| MNX1 | 6 | -0.23618 | -0.10965 | +0.12653 | 0.1305 | +0.018 | False | +0.015 | False | False |
| PRDM1 | 6 | -0.23410 | +0.00211 | +0.23622 | 0.1305 | +0.031 | False | +0.073 | False | False |
| PIAS1 | 6 | -0.21816 | -0.00952 | +0.20864 | 0.1305 | +0.000 | True | +0.015 | False | False |
| IRF1 | 6 | -0.21396 | -0.04845 | +0.16551 | 0.1305 | +0.537 | False | +0.025 | False | False |
| ETS1 | 6 | -0.19745 | -0.06426 | +0.13319 | 0.1305 | +0.000 | True | +0.041 | False | False |
| ZFP2 | 6 | -0.18757 | -0.08197 | +0.10559 | 0.1305 | +0.000 | True | +0.042 | False | False |
| TBX20 | 6 | -0.17291 | -0.06226 | +0.11066 | 0.1305 | +0.012 | False | +0.077 | False | False |
| ZNF763 | 6 | -0.15757 | -0.03345 | +0.12412 | 0.1305 | +0.000 | True | +0.053 | False | False |
| FOXD2 | 6 | -0.15525 | -0.05963 | +0.09562 | 0.1305 | +0.000 | True | +0.146 | False | False |
| SP140 | 6 | -0.15082 | +0.04215 | +0.19296 | 0.1305 | +0.000 | True | -0.061 | False | False |
| RHOXF2B | 4 | -0.14824 | -0.02770 | +0.12054 | 0.1305 | +0.000 | True | -0.000 | False | False |

**V_sig / S2_k20 / O2_Y1** — 0 of these 20 are `toward_young`; 0 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 | 6 | -1.56271 | -0.31318 | +1.24953 | 0.1202 | +0.579 | False | -0.081 | False | False |
| KLF5 | 6 | -0.54361 | -0.17495 | +0.36866 | 0.1202 | +0.157 | False | +0.126 | False | False |
| ZEB2 | 6 | -0.35116 | -0.09612 | +0.25504 | 0.1202 | +0.037 | False | -0.208 | False | False |
| GFI1 | 6 | -0.34893 | -0.13092 | +0.21801 | 0.1202 | +0.370 | False | +0.061 | False | False |
| BSX | 6 | -0.34840 | -0.14874 | +0.19966 | 0.1202 | +0.003 | False | +0.021 | False | False |
| CREB5 | 6 | -0.33051 | -0.22218 | +0.10833 | 0.1943 | +0.000 | True | -0.012 | False | False |
| NCOA1 | 6 | -0.31245 | -0.13258 | +0.17986 | 0.1202 | +0.000 | True | +0.007 | False | False |
| NKX3-1 | 6 | -0.25801 | -0.16138 | +0.09663 | 0.1202 | +0.205 | False | +0.061 | False | False |
| MXI1 | 6 | -0.25336 | -0.08568 | +0.16768 | 0.1202 | +0.012 | False | -0.000 | False | False |
| JDP2 | 6 | -0.24437 | -0.07094 | +0.17342 | 0.1202 | +0.121 | False | +0.098 | False | False |
| AFF4 | 6 | -0.23121 | -0.11843 | +0.11278 | 0.1589 | +0.000 | True | -0.000 | False | False |
| MYOG | 6 | -0.22093 | -0.05405 | +0.16688 | 0.1202 | +0.000 | True | -0.049 | False | False |
| GATA1 | 6 | -0.20039 | -0.05230 | +0.14809 | 0.1202 | +0.000 | True | +0.019 | False | False |
| PIAS1 | 6 | -0.19780 | -0.08739 | +0.11042 | 0.1589 | +0.000 | True | +0.015 | False | False |
| HHEX | 6 | -0.19725 | -0.12594 | +0.07131 | 0.1943 | +0.002 | False | +0.053 | False | False |
| PAX2 | 6 | -0.18867 | -0.07985 | +0.10881 | 0.1202 | +0.119 | False | -0.149 | False | False |
| POU3F3 | 6 | -0.18418 | -0.04586 | +0.13832 | 0.1202 | +0.000 | True | -0.040 | False | False |
| MNX1 | 6 | -0.17668 | -0.09479 | +0.08189 | 0.1202 | +0.018 | False | +0.015 | False | False |
| ZFP2 | 6 | -0.17415 | -0.06415 | +0.11000 | 0.1202 | +0.000 | True | +0.042 | False | False |
| ZNF420 | 6 | -0.16158 | -0.08312 | +0.07846 | 0.1202 | +0.000 | True | +0.008 | False | False |

**V_sig / S2_k20 / O2_Y2** — 0 of these 20 are `toward_young`; 0 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 | 6 | -1.69098 | -0.30907 | +1.38191 | 0.1171 | +0.579 | False | -0.081 | False | False |
| KLF5 | 6 | -0.57713 | -0.16950 | +0.40763 | 0.1171 | +0.157 | False | +0.126 | False | False |
| ZEB2 | 6 | -0.38152 | -0.10800 | +0.27352 | 0.1171 | +0.037 | False | -0.208 | False | False |
| BSX | 6 | -0.35605 | -0.14855 | +0.20751 | 0.1171 | +0.003 | False | +0.021 | False | False |
| CREB5 | 6 | -0.34534 | -0.25508 | +0.09026 | 0.2076 | +0.000 | True | -0.012 | False | False |
| GFI1 | 6 | -0.33971 | -0.12801 | +0.21169 | 0.1171 | +0.370 | False | +0.061 | False | False |
| NCOA1 | 6 | -0.31160 | -0.13870 | +0.17290 | 0.1171 | +0.000 | True | +0.007 | False | False |
| NKX3-1 | 6 | -0.31147 | -0.15495 | +0.15652 | 0.1171 | +0.205 | False | +0.061 | False | False |
| JDP2 | 6 | -0.27962 | -0.07820 | +0.20142 | 0.1171 | +0.121 | False | +0.098 | False | False |
| MXI1 | 6 | -0.27237 | -0.09038 | +0.18199 | 0.1171 | +0.012 | False | -0.000 | False | False |
| AFF4 | 6 | -0.25003 | -0.13582 | +0.11420 | 0.1631 | +0.000 | True | -0.000 | False | False |
| MYOG | 6 | -0.24490 | -0.05706 | +0.18784 | 0.1171 | +0.000 | True | -0.049 | False | False |
| GATA1 | 6 | -0.23109 | -0.05751 | +0.17358 | 0.1171 | +0.000 | True | +0.019 | False | False |
| PIAS1 | 6 | -0.22840 | -0.10633 | +0.12207 | 0.1171 | +0.000 | True | +0.015 | False | False |
| PAX2 | 6 | -0.21040 | -0.08762 | +0.12278 | 0.1171 | +0.119 | False | -0.149 | False | False |
| HHEX | 6 | -0.21015 | -0.12042 | +0.08973 | 0.1631 | +0.002 | False | +0.053 | False | False |
| MNX1 | 6 | -0.20178 | -0.09580 | +0.10597 | 0.1171 | +0.018 | False | +0.015 | False | False |
| POU3F3 | 6 | -0.19865 | -0.04473 | +0.15391 | 0.1171 | +0.000 | True | -0.040 | False | False |
| HNRNPAB | 6 | -0.19291 | -0.21464 | -0.02173 | 0.4291 | +0.000 | True | -0.000 | False | False |
| ZFP2 | 6 | -0.18434 | -0.06244 | +0.12190 | 0.1171 | +0.000 | True | +0.042 | False | False |

**V_sig / S2_k50 / O1_Y1** — 0 of these 20 are `toward_young`; 0 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 | 6 | -0.50325 | -0.44310 | +0.06016 | 0.1883 | +0.579 | False | -0.081 | False | False |
| BSX | 6 | -0.45314 | -0.25761 | +0.19553 | 0.0972 | +0.003 | False | +0.021 | False | False |
| GFI1 | 6 | -0.43713 | -0.17014 | +0.26700 | 0.0972 | +0.370 | False | +0.061 | False | False |
| HHEX | 6 | -0.28901 | -0.17829 | +0.11071 | 0.0972 | +0.002 | False | +0.053 | False | False |
| OLIG1 | 6 | -0.20966 | -0.06143 | +0.14823 | 0.0972 | +0.000 | True | +0.026 | False | False |
| KLF5 | 6 | -0.20653 | -0.19530 | +0.01123 | 0.2887 | +0.157 | False | +0.126 | False | False |
| PRDM1 | 6 | -0.20250 | +0.11423 | +0.31673 | 0.0972 | +0.031 | False | +0.073 | False | False |
| PRDM16 | 6 | -0.18431 | -0.02879 | +0.15552 | 0.0972 | +0.006 | False | -0.111 | False | False |
| NKX3-1 | 6 | -0.17257 | -0.16885 | +0.00372 | 0.2887 | +0.205 | False | +0.061 | False | False |
| CDX2 | 6 | -0.16898 | -0.03453 | +0.13445 | 0.0972 | +0.215 | False | +0.041 | False | False |
| SOX11 | 5 | -0.16510 | -0.02393 | +0.14117 | 0.0972 | +0.001 | False | +0.010 | False | False |
| RBMS1 | 6 | -0.16196 | -0.06544 | +0.09652 | 0.0972 | +0.000 | True | +0.037 | False | False |
| ZNF714 | 6 | -0.16075 | -0.01600 | +0.14475 | 0.1427 | +0.000 | True | -0.000 | False | False |
| ZNF763 | 6 | -0.16054 | -0.00635 | +0.15419 | 0.0972 | +0.000 | True | +0.053 | False | False |
| MYOD1 | 6 | -0.15921 | -0.05664 | +0.10256 | 0.0972 | +0.000 | True | +0.013 | False | False |
| ETS1 | 6 | -0.15916 | -0.04789 | +0.11127 | 0.0972 | +0.000 | True | +0.041 | False | False |
| LHX4 | 6 | -0.15258 | -0.02062 | +0.13196 | 0.0972 | +0.021 | False | -0.039 | False | False |
| POU3F3 | 6 | -0.13171 | -0.06040 | +0.07131 | 0.0972 | +0.000 | True | -0.040 | False | False |
| MXI1 | 6 | -0.12489 | -0.01525 | +0.10963 | 0.0972 | +0.012 | False | -0.000 | False | False |
| SOX5 | 6 | -0.11923 | +0.02070 | +0.13994 | 0.0972 | +0.000 | True | +0.012 | False | False |

**V_sig / S2_k50 / O1_Y2** — 0 of these 20 are `toward_young`; 0 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 | 6 | -1.47333 | -0.37233 | +1.10100 | 0.1363 | +0.579 | False | -0.081 | False | False |
| BSX | 6 | -0.50456 | -0.19680 | +0.30776 | 0.1363 | +0.003 | False | +0.021 | False | False |
| NKX3-1 | 6 | -0.50036 | -0.16011 | +0.34025 | 0.1363 | +0.205 | False | +0.061 | False | False |
| KLF5 | 6 | -0.49372 | -0.17334 | +0.32038 | 0.1363 | +0.157 | False | +0.126 | False | False |
| GFI1 | 6 | -0.39506 | -0.14415 | +0.25092 | 0.1363 | +0.370 | False | +0.061 | False | False |
| HHEX | 6 | -0.35536 | -0.14557 | +0.20979 | 0.1363 | +0.002 | False | +0.053 | False | False |
| HNRNPAB | 6 | -0.31781 | -0.01365 | +0.30416 | 0.1363 | +0.000 | True | -0.000 | False | False |
| PRDM1 | 6 | -0.26395 | +0.00401 | +0.26796 | 0.1363 | +0.031 | False | +0.073 | False | False |
| MXI1 | 6 | -0.25680 | -0.05130 | +0.20550 | 0.1363 | +0.012 | False | -0.000 | False | False |
| POU3F3 | 6 | -0.23279 | -0.05083 | +0.18196 | 0.1363 | +0.000 | True | -0.040 | False | False |
| ETS1 | 6 | -0.17890 | -0.06180 | +0.11710 | 0.1363 | +0.000 | True | +0.041 | False | False |
| ZFP2 | 6 | -0.17601 | -0.08027 | +0.09574 | 0.1363 | +0.000 | True | +0.042 | False | False |
| MNX1 | 6 | -0.16814 | -0.10418 | +0.06397 | 0.1363 | +0.018 | False | +0.015 | False | False |
| TBX20 | 6 | -0.16453 | -0.06239 | +0.10214 | 0.1363 | +0.012 | False | +0.077 | False | False |
| PIAS1 | 6 | -0.15752 | -0.01729 | +0.14022 | 0.1363 | +0.000 | True | +0.015 | False | False |
| ELOC | 6 | -0.15741 | -0.05269 | +0.10472 | 0.1363 | +0.000 | True | -0.020 | False | False |
| ZNF763 | 6 | -0.14628 | -0.03520 | +0.11108 | 0.1363 | +0.000 | True | +0.053 | False | False |
| OLIG1 | 6 | -0.14117 | -0.05329 | +0.08788 | 0.1363 | +0.000 | True | +0.026 | False | False |
| AFF4 | 6 | -0.13512 | -0.04218 | +0.09294 | 0.2265 | +0.000 | True | -0.000 | False | False |
| PRDM16 | 6 | -0.12826 | -0.05285 | +0.07541 | 0.1363 | +0.006 | False | -0.111 | False | False |

**V_sig / S2_k50 / O2_Y1** — 0 of these 20 are `toward_young`; 0 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 | 6 | -1.59107 | -0.33135 | +1.25972 | 0.0808 | +0.579 | False | -0.081 | False | False |
| KLF5 | 6 | -0.52448 | -0.18130 | +0.34318 | 0.0808 | +0.157 | False | +0.126 | False | False |
| CREB5 | 6 | -0.41561 | -0.14022 | +0.27539 | 0.0808 | +0.000 | True | -0.012 | False | False |
| GFI1 | 6 | -0.39884 | -0.12993 | +0.26891 | 0.0808 | +0.370 | False | +0.061 | False | False |
| BSX | 6 | -0.36578 | -0.14295 | +0.22283 | 0.0808 | +0.003 | False | +0.021 | False | False |
| NCOA1 | 6 | -0.33414 | -0.10201 | +0.23212 | 0.0808 | +0.000 | True | +0.007 | False | False |
| SP140 | 6 | -0.29254 | -0.11820 | +0.17434 | 0.0808 | +0.000 | True | -0.061 | False | False |
| NKX3-1 | 6 | -0.28238 | -0.17113 | +0.11125 | 0.0808 | +0.205 | False | +0.061 | False | False |
| ZEB2 | 6 | -0.27964 | -0.07948 | +0.20016 | 0.0808 | +0.037 | False | -0.208 | False | False |
| PRDM1 | 6 | -0.27917 | -0.11932 | +0.15986 | 0.0808 | +0.031 | False | +0.073 | False | False |
| HHEX | 6 | -0.24886 | -0.14099 | +0.10787 | 0.0808 | +0.002 | False | +0.053 | False | False |
| BCL6B | 6 | -0.24752 | -0.09535 | +0.15218 | 0.0808 | +0.000 | True | -0.000 | False | False |
| MYOG | 6 | -0.24270 | -0.04865 | +0.19405 | 0.0808 | +0.000 | True | -0.049 | False | False |
| GATA1 | 6 | -0.23835 | -0.04706 | +0.19129 | 0.0808 | +0.000 | True | +0.019 | False | False |
| AFF4 | 6 | -0.22533 | -0.10233 | +0.12300 | 0.0808 | +0.000 | True | -0.000 | False | False |
| MXI1 | 6 | -0.21831 | -0.07931 | +0.13900 | 0.0808 | +0.012 | False | -0.000 | False | False |
| PAX2 | 6 | -0.21684 | -0.08700 | +0.12983 | 0.0808 | +0.119 | False | -0.149 | False | False |
| ZFP2 | 6 | -0.21130 | -0.06662 | +0.14468 | 0.0808 | +0.000 | True | +0.042 | False | False |
| PROX1 | 6 | -0.19275 | -0.06792 | +0.12483 | 0.0808 | +0.138 | False | -0.141 | False | False |
| MNX1 | 6 | -0.19087 | -0.10125 | +0.08962 | 0.0808 | +0.018 | False | +0.015 | False | False |

**V_sig / S2_k50 / O2_Y2** — 0 of these 20 are `toward_young`; 0 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 | 6 | -1.68975 | -0.32374 | +1.36601 | 0.0794 | +0.579 | False | -0.081 | False | False |
| KLF5 | 6 | -0.55402 | -0.18134 | +0.37268 | 0.0794 | +0.157 | False | +0.126 | False | False |
| CREB5 | 6 | -0.41949 | -0.16221 | +0.25728 | 0.0794 | +0.000 | True | -0.012 | False | False |
| GFI1 | 6 | -0.39160 | -0.12956 | +0.26204 | 0.0794 | +0.370 | False | +0.061 | False | False |
| BSX | 6 | -0.37179 | -0.13368 | +0.23811 | 0.0794 | +0.003 | False | +0.021 | False | False |
| NCOA1 | 6 | -0.33330 | -0.11009 | +0.22320 | 0.0794 | +0.000 | True | +0.007 | False | False |
| NKX3-1 | 6 | -0.32204 | -0.17472 | +0.14732 | 0.0794 | +0.205 | False | +0.061 | False | False |
| ZEB2 | 6 | -0.30597 | -0.08848 | +0.21749 | 0.0794 | +0.037 | False | -0.208 | False | False |
| SP140 | 6 | -0.30251 | -0.12273 | +0.17977 | 0.0794 | +0.000 | True | -0.061 | False | False |
| PRDM1 | 6 | -0.28440 | -0.12617 | +0.15822 | 0.0794 | +0.031 | False | +0.073 | False | False |
| GATA1 | 6 | -0.26186 | -0.05041 | +0.21145 | 0.0794 | +0.000 | True | +0.019 | False | False |
| MYOG | 6 | -0.26154 | -0.05126 | +0.21029 | 0.0794 | +0.000 | True | -0.049 | False | False |
| HHEX | 6 | -0.25678 | -0.14054 | +0.11624 | 0.0794 | +0.002 | False | +0.053 | False | False |
| BCL6B | 6 | -0.24917 | -0.09702 | +0.15215 | 0.0794 | +0.000 | True | -0.000 | False | False |
| AFF4 | 6 | -0.23833 | -0.11554 | +0.12279 | 0.0794 | +0.000 | True | -0.000 | False | False |
| PAX2 | 6 | -0.23295 | -0.09084 | +0.14211 | 0.0794 | +0.119 | False | -0.149 | False | False |
| MXI1 | 6 | -0.23293 | -0.08366 | +0.14926 | 0.0794 | +0.012 | False | -0.000 | False | False |
| ZFP2 | 6 | -0.21860 | -0.06337 | +0.15523 | 0.0794 | +0.000 | True | +0.042 | False | False |
| MNX1 | 6 | -0.21039 | -0.09528 | +0.11511 | 0.0794 | +0.018 | False | +0.015 | False | False |
| PROX1 | 6 | -0.19136 | -0.06798 | +0.12337 | 0.0794 | +0.138 | False | -0.141 | False | False |

**V_sig / S2_k100 / O1_Y1** — 0 of these 20 are `toward_young`; 0 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| BSX | 6 | -0.43806 | -0.25260 | +0.18546 | 0.1142 | +0.003 | False | +0.021 | False | False |
| GFI1 | 6 | -0.37277 | -0.16400 | +0.20877 | 0.1142 | +0.370 | False | +0.061 | False | False |
| HHEX | 6 | -0.24377 | -0.17615 | +0.06762 | 0.1646 | +0.002 | False | +0.053 | False | False |
| OLIG1 | 6 | -0.19165 | -0.05991 | +0.13174 | 0.1142 | +0.000 | True | +0.026 | False | False |
| ZNF714 | 6 | -0.19024 | -0.01835 | +0.17189 | 0.1142 | +0.000 | True | -0.000 | False | False |
| PRDM16 | 6 | -0.18220 | -0.03701 | +0.14519 | 0.1142 | +0.006 | False | -0.111 | False | False |
| CDX2 | 6 | -0.17600 | -0.03012 | +0.14588 | 0.1142 | +0.215 | False | +0.041 | False | False |
| PRDM1 | 6 | -0.17451 | +0.11832 | +0.29283 | 0.1142 | +0.031 | False | +0.073 | False | False |
| KLF5 | 6 | -0.17384 | -0.18208 | -0.00825 | 0.3422 | +0.157 | False | +0.126 | False | False |
| SOX11 | 5 | -0.17159 | -0.01788 | +0.15372 | 0.1142 | +0.001 | False | +0.010 | False | False |
| MYOD1 | 6 | -0.16035 | -0.05580 | +0.10455 | 0.1142 | +0.000 | True | +0.013 | False | False |
| RBMS1 | 6 | -0.15518 | -0.06253 | +0.09265 | 0.1142 | +0.000 | True | +0.037 | False | False |
| ZNF763 | 6 | -0.14534 | -0.00766 | +0.13768 | 0.1142 | +0.000 | True | +0.053 | False | False |
| ETS1 | 6 | -0.14136 | -0.05178 | +0.08958 | 0.1142 | +0.000 | True | +0.041 | False | False |
| NKX3-1 | 6 | -0.13391 | -0.16313 | -0.02922 | 0.4350 | +0.205 | False | +0.061 | False | False |
| HELZ2 | 6 | -0.13202 | -0.42154 | -0.28951 | 0.7484 | +0.579 | False | -0.081 | False | False |
| LHX4 | 6 | -0.12712 | -0.02311 | +0.10401 | 0.1142 | +0.021 | False | -0.039 | False | False |
| POU3F3 | 6 | -0.12679 | -0.05715 | +0.06964 | 0.1142 | +0.000 | True | -0.040 | False | False |
| MXI1 | 6 | -0.11846 | -0.00512 | +0.11334 | 0.1142 | +0.012 | False | -0.000 | False | False |
| ARID1A | 6 | -0.10880 | +0.01041 | +0.11921 | 0.1142 | +0.000 | True | +0.057 | False | False |

**V_sig / S2_k100 / O1_Y2** — 0 of these 20 are `toward_young`; 0 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 | 6 | -1.10573 | -0.34073 | +0.76500 | 0.1497 | +0.579 | False | -0.081 | False | False |
| BSX | 6 | -0.48223 | -0.18524 | +0.29698 | 0.1497 | +0.003 | False | +0.021 | False | False |
| NKX3-1 | 6 | -0.45610 | -0.16014 | +0.29596 | 0.1497 | +0.205 | False | +0.061 | False | False |
| KLF5 | 6 | -0.45305 | -0.17234 | +0.28071 | 0.1497 | +0.157 | False | +0.126 | False | False |
| GFI1 | 6 | -0.32678 | -0.14113 | +0.18565 | 0.1497 | +0.370 | False | +0.061 | False | False |
| HHEX | 6 | -0.30507 | -0.14756 | +0.15751 | 0.1497 | +0.002 | False | +0.053 | False | False |
| HNRNPAB | 6 | -0.29884 | -0.00509 | +0.29375 | 0.1497 | +0.000 | True | -0.000 | False | False |
| MXI1 | 6 | -0.24589 | -0.04628 | +0.19961 | 0.1497 | +0.012 | False | -0.000 | False | False |
| PRDM1 | 6 | -0.22807 | +0.01026 | +0.23833 | 0.1497 | +0.031 | False | +0.073 | False | False |
| POU3F3 | 6 | -0.22803 | -0.04970 | +0.17834 | 0.1497 | +0.000 | True | -0.040 | False | False |
| MNX1 | 6 | -0.16785 | -0.09808 | +0.06977 | 0.1497 | +0.018 | False | +0.015 | False | False |
| ZFP2 | 6 | -0.16380 | -0.07423 | +0.08957 | 0.1497 | +0.000 | True | +0.042 | False | False |
| ETS1 | 6 | -0.16028 | -0.06015 | +0.10013 | 0.1497 | +0.000 | True | +0.041 | False | False |
| ELOC | 6 | -0.15866 | -0.05089 | +0.10777 | 0.1497 | +0.000 | True | -0.020 | False | False |
| AFF4 | 6 | -0.13956 | -0.03667 | +0.10289 | 0.2149 | +0.000 | True | -0.000 | False | False |
| ZNF714 | 6 | -0.13492 | -0.09852 | +0.03640 | 0.2854 | +0.000 | True | -0.000 | False | False |
| ZNF763 | 6 | -0.13422 | -0.03058 | +0.10364 | 0.1497 | +0.000 | True | +0.053 | False | False |
| PIAS1 | 6 | -0.13088 | +0.00054 | +0.13142 | 0.1497 | +0.000 | True | +0.015 | False | False |
| PRDM16 | 6 | -0.12990 | -0.05678 | +0.07313 | 0.1497 | +0.006 | False | -0.111 | False | False |
| TBX20 | 6 | -0.12974 | -0.05978 | +0.06996 | 0.1497 | +0.012 | False | +0.077 | False | False |

**V_sig / S2_k100 / O2_Y1** — 0 of these 20 are `toward_young`; 0 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 | 6 | -1.58452 | -0.38253 | +1.20199 | 0.0801 | +0.579 | False | -0.081 | False | False |
| KLF5 | 6 | -0.59052 | -0.19036 | +0.40016 | 0.0801 | +0.157 | False | +0.126 | False | False |
| GFI1 | 6 | -0.46090 | -0.14921 | +0.31170 | 0.0801 | +0.370 | False | +0.061 | False | False |
| BSX | 6 | -0.42553 | -0.17681 | +0.24872 | 0.0801 | +0.003 | False | +0.021 | False | False |
| HHEX | 6 | -0.30767 | -0.16644 | +0.14122 | 0.0801 | +0.002 | False | +0.053 | False | False |
| NCOA1 | 6 | -0.28039 | -0.07243 | +0.20796 | 0.0801 | +0.000 | True | +0.007 | False | False |
| CREB5 | 6 | -0.26621 | -0.00091 | +0.26530 | 0.0801 | +0.000 | True | -0.012 | False | False |
| BCL6B | 6 | -0.25634 | -0.07546 | +0.18088 | 0.0801 | +0.000 | True | -0.000 | False | False |
| NKX3-1 | 6 | -0.25453 | -0.16896 | +0.08556 | 0.0801 | +0.205 | False | +0.061 | False | False |
| ZEB2 | 6 | -0.24408 | -0.02383 | +0.22024 | 0.0801 | +0.037 | False | -0.208 | False | False |
| GATA1 | 6 | -0.24381 | -0.02888 | +0.21493 | 0.0801 | +0.000 | True | +0.019 | False | False |
| PRDM1 | 6 | -0.24265 | -0.06339 | +0.17926 | 0.0801 | +0.031 | False | +0.073 | False | False |
| MYOG | 6 | -0.24037 | -0.03639 | +0.20398 | 0.0801 | +0.000 | True | -0.049 | False | False |
| ZFP2 | 6 | -0.23952 | -0.07374 | +0.16578 | 0.0801 | +0.000 | True | +0.042 | False | False |
| PAX2 | 6 | -0.23450 | -0.05917 | +0.17533 | 0.0801 | +0.119 | False | -0.149 | False | False |
| SP140 | 6 | -0.21428 | -0.03791 | +0.17637 | 0.0801 | +0.000 | True | -0.061 | False | False |
| MNX1 | 6 | -0.19978 | -0.10050 | +0.09927 | 0.0801 | +0.018 | False | +0.015 | False | False |
| FOXD2 | 6 | -0.18703 | -0.06731 | +0.11972 | 0.0801 | +0.000 | True | +0.146 | False | False |
| AFF4 | 6 | -0.18525 | -0.05576 | +0.12949 | 0.1473 | +0.000 | True | -0.000 | False | False |
| ZNF574 | 6 | -0.18191 | -0.04481 | +0.13710 | 0.0801 | +0.000 | True | +0.027 | False | False |

**V_sig / S2_k100 / O2_Y2** — 0 of these 20 are `toward_young`; 0 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 | 6 | -1.67974 | -0.37234 | +1.30740 | 0.0801 | +0.579 | False | -0.081 | False | False |
| KLF5 | 6 | -0.61752 | -0.18346 | +0.43406 | 0.0801 | +0.157 | False | +0.126 | False | False |
| GFI1 | 6 | -0.45190 | -0.14469 | +0.30721 | 0.0801 | +0.370 | False | +0.061 | False | False |
| BSX | 6 | -0.42978 | -0.17185 | +0.25794 | 0.0801 | +0.003 | False | +0.021 | False | False |
| HHEX | 6 | -0.31337 | -0.16668 | +0.14669 | 0.0801 | +0.002 | False | +0.053 | False | False |
| NKX3-1 | 6 | -0.29298 | -0.16791 | +0.12507 | 0.0801 | +0.205 | False | +0.061 | False | False |
| NCOA1 | 6 | -0.28311 | -0.07946 | +0.20365 | 0.0801 | +0.000 | True | +0.007 | False | False |
| CREB5 | 6 | -0.27269 | -0.03086 | +0.24184 | 0.0801 | +0.000 | True | -0.012 | False | False |
| ZEB2 | 6 | -0.27004 | -0.03392 | +0.23612 | 0.0801 | +0.037 | False | -0.208 | False | False |
| GATA1 | 6 | -0.26574 | -0.03598 | +0.22976 | 0.0801 | +0.000 | True | +0.019 | False | False |
| MYOG | 6 | -0.25887 | -0.04072 | +0.21815 | 0.0801 | +0.000 | True | -0.049 | False | False |
| BCL6B | 6 | -0.25682 | -0.07766 | +0.17916 | 0.0801 | +0.000 | True | -0.000 | False | False |
| PAX2 | 6 | -0.24912 | -0.06553 | +0.18359 | 0.0801 | +0.119 | False | -0.149 | False | False |
| PRDM1 | 6 | -0.24755 | -0.07773 | +0.16982 | 0.0801 | +0.031 | False | +0.073 | False | False |
| ZFP2 | 6 | -0.24627 | -0.06964 | +0.17662 | 0.0801 | +0.000 | True | +0.042 | False | False |
| SP140 | 6 | -0.22462 | -0.04526 | +0.17937 | 0.0801 | +0.000 | True | -0.061 | False | False |
| MNX1 | 6 | -0.21809 | -0.10035 | +0.11774 | 0.0801 | +0.018 | False | +0.015 | False | False |
| AFF4 | 6 | -0.19777 | -0.06325 | +0.13452 | 0.1522 | +0.000 | True | -0.000 | False | False |
| FOXD2 | 6 | -0.18805 | -0.06797 | +0.12008 | 0.0801 | +0.000 | True | +0.146 | False | False |
| CLOCK | 6 | -0.18461 | -0.05286 | +0.13175 | 0.0801 | +0.000 | True | -0.004 | False | False |

### What produced these counts

Several of the counts above read backwards without their mechanism, so each is stated.

**Why each of the 19 zero settings is zero.** For **18 of 19** the binding constraint is the q-value gate: no factor reaches q <= 0.05 there. That is a resolution limit of the pre-registered n_perm = 200 (next paragraph), not a finding that those factors were tested and found null. No zero setting is zero because factors cleared every statistical bar and were then removed by a flag.

| variant | space | pair | binding_constraint |
|---|---|---|---|
| V_all | S1 | O1_Y1 | no factor has a negative delta |
| V_all | S2_k20 | O1_Y2 | no factor reaches q <= 0.05 |
| V_all | S2_k50 | O1_Y2 | no factor reaches q <= 0.05 |
| V_all | S2_k100 | O1_Y1 | no factor reaches q <= 0.05 |
| V_all | S2_k100 | O1_Y2 | no factor reaches q <= 0.05 |
| V_sig | S1 | O1_Y1 | no factor has a negative delta; no factor reaches q <= 0.05 |
| V_sig | S1 | O2_Y1 | no factor reaches q <= 0.05 |
| V_sig | S2_k20 | O1_Y1 | no factor reaches q <= 0.05 |
| V_sig | S2_k20 | O1_Y2 | no factor reaches q <= 0.05 |
| V_sig | S2_k20 | O2_Y1 | no factor reaches q <= 0.05 |
| V_sig | S2_k20 | O2_Y2 | no factor reaches q <= 0.05 |
| V_sig | S2_k50 | O1_Y1 | no factor reaches q <= 0.05 |
| V_sig | S2_k50 | O1_Y2 | no factor reaches q <= 0.05 |
| V_sig | S2_k50 | O2_Y1 | no factor reaches q <= 0.05 |
| V_sig | S2_k50 | O2_Y2 | no factor reaches q <= 0.05 |
| V_sig | S2_k100 | O1_Y1 | no factor reaches q <= 0.05 |
| V_sig | S2_k100 | O1_Y2 | no factor reaches q <= 0.05 |
| V_sig | S2_k100 | O2_Y1 | no factor reaches q <= 0.05 |
| V_sig | S2_k100 | O2_Y2 | no factor reaches q <= 0.05 |

**In 2 of 28 settings not one factor has a negative delta, so `toward_young` = 0 there regardless of significance.** Adding a displacement of norm ||d|| to a point at distance ||v|| from the target reduces that distance only if cos(d, v) > ||d|| / (2||v||); a large step needs to be well aimed or it overshoots. Below is that required cosine against the cosine actually achieved. Meeting one's own requirement is the same event as delta < 0, so the count column restates the result; the informative comparison is how far the best achieved cosine sits from the required one:

| variant | space | pair | origin_to_target | median_factor_norm | median_required_cos | smallest_required_cos | observed_max_cos | n_factors_meeting_own_requirement | smallest_delta |
|---|---|---|---|---|---|---|---|---|---|
| V_all | S1 | O1_Y1 | 32.284 | 20.254 | 0.3137 | 0.1709 | 0.0974 | 0 | +1.5303 |
| V_sig | S1 | O1_Y1 | 32.284 | 3.289 | 0.0509 | 0.0115 | 0.0535 | 0 | +0.0049 |

In those settings the `below_own_floor` count is still non-zero. That is consistent: a factor can move *less far away* than a random direction of its own size would, which clears its floor while still increasing the distance. It is not an approach, and the pre-registered `toward_young` rule requires delta < 0, so none of them is labelled a candidate.

**The BH-FDR gate is resolution-limited by the pre-registered n_perm = 200.** The smallest p an empirical one-sided permutation test can return is 1/(200+1) = 0.00498. With BH across 1,836 factors, q <= 0.05 is unreachable unless at least **183 factors are tied at that minimum p** in the same setting. So a setting with `q_le_005` = 0 is reporting that fewer than that many factors hit the p floor; it is not reporting that the remaining factors were tested and found null. Observed ties at the p floor, per setting, are in `stage3_full_table.csv`.

**`guides_disagree` in V_sig is almost always an exact zero, not a disagreement.** 1,692 of 1,836 factors have a median pairwise guide cosine of **exactly 0** in V_sig, 0 are negative, and 143 are positive. A cosine of exactly 0 here means two guides for the same factor share **no significant gene at all**, so their masked vectors are orthogonal by construction. The pre-registered rule flags `median <= 0`, so these fire. They are reported as fired and are not relabelled. The finding is that per-guide `adj_p <= 0.05` masking leaves most guides of the same factor with disjoint support. This is **not the sole reason** V_sig has no candidate: 14 of 14 V_sig settings also have no factor at q <= 0.05, so the q-value gate alone already returns zero there, even for the 143 factors whose guides do agree.

**The identity check was live and never fired.** The largest identity-panel drop over all factors and variants is +0.4952 in frozen-z units against a threshold of 0.5, and 3 factor-variants exceed 0.4. So `identity_loss` = 0 everywhere is a real result at this threshold, not a check that could never have fired; it came within 0.0048 of firing.

**The Y2 membership test has no discriminating power for O1, because the aged origin is already inside the young region.** At k=20, k=50, O1's own Mahalanobis distance to the Y2 region is below the membership radius before any perturbation is added (see the distance table in Stage 1), so `origin + d` lands inside the radius for every one of the 1,836 factors. That is not evidence that any factor rejuvenates anything; it is a statement about the GTEx young and old centroids being close relative to the spread among young donors. For O2 the opposite holds: O2 starts outside the radius and no factor brings it inside.

**The size of the best approach, stated so the deltas are not read as large.** The single most negative delta in each setting, with the fraction of the origin-to-target gap it closes (`frac`):

| variant | space | pair | factor | delta | origin_to_target | frac_of_gap_closed | cos | toward_young |
|---|---|---|---|---|---|---|---|---|
| V_all | S2_k100 | O2_Y2 | HELZ2 | -3.8216 | 277.476 | +0.01496 | +0.2956 | True |
| V_all | S2_k50 | O2_Y2 | HELZ2 | -3.7922 | 265.189 | +0.01531 | +0.3248 | True |
| V_all | S2_k100 | O2_Y1 | HELZ2 | -3.7518 | 272.045 | +0.01503 | +0.2912 | True |
| V_all | S2_k20 | O2_Y2 | HELZ2 | -3.7499 | 204.376 | +0.01950 | +0.3791 | True |
| V_all | S2_k50 | O2_Y1 | HELZ2 | -3.7225 | 259.334 | +0.01541 | +0.3197 | True |
| V_all | S2_k20 | O2_Y1 | HELZ2 | -3.6659 | 196.450 | +0.01992 | +0.3721 | True |
| V_all | S2_k20 | O1_Y2 | HELZ2 | -2.8536 | 33.210 | +0.13236 | +0.4180 | False |
| V_all | S2_k20 | O1_Y1 | NCL | -2.7385 | 28.645 | +0.11444 | +0.5289 | True |
| V_all | S1 | O2_Y1 | HELZ2 | -2.6565 | 468.388 | +0.00903 | +0.1099 | True |
| V_all | S2_k50 | O1_Y1 | PRDM1 | -2.3883 | 30.864 | +0.12562 | +0.3924 | True |
| V_all | S2_k50 | O1_Y2 | HELZ2 | -2.1785 | 35.168 | +0.12320 | +0.3466 | False |
| V_all | S2_k100 | O1_Y1 | PRDM1 | -1.9680 | 31.434 | +0.12149 | +0.3483 | False |
| V_sig | S2_k20 | O1_Y2 | HELZ2 | -1.7690 | 33.210 | +0.06971 | +0.3688 | False |
| V_sig | S2_k20 | O2_Y2 | HELZ2 | -1.6910 | 204.376 | +0.00871 | +0.2836 | False |
| V_sig | S2_k50 | O2_Y2 | HELZ2 | -1.6897 | 265.189 | +0.00678 | +0.2321 | False |
| V_sig | S2_k100 | O2_Y2 | HELZ2 | -1.6797 | 277.476 | +0.00658 | +0.1998 | False |
| V_sig | S2_k50 | O2_Y1 | HELZ2 | -1.5911 | 259.334 | +0.00656 | +0.2197 | False |
| V_sig | S2_k100 | O2_Y1 | HELZ2 | -1.5845 | 272.045 | +0.00637 | +0.1898 | False |
| V_sig | S2_k20 | O2_Y1 | HELZ2 | -1.5627 | 196.450 | +0.00843 | +0.2639 | False |
| V_all | S2_k100 | O1_Y2 | HELZ2 | -1.5553 | 35.835 | +0.11921 | +0.3043 | False |
| V_sig | S1 | O2_Y1 | HELZ2 | -1.4761 | 468.388 | +0.00462 | +0.0851 | False |
| V_sig | S2_k50 | O1_Y2 | HELZ2 | -1.4733 | 35.168 | +0.06526 | +0.2964 | False |
| V_sig | S2_k100 | O1_Y2 | HELZ2 | -1.1057 | 35.835 | +0.06285 | +0.2466 | False |
| V_sig | S2_k20 | O1_Y1 | HELZ2 | -0.6636 | 28.645 | +0.04691 | +0.2141 | False |
| V_sig | S2_k50 | O1_Y1 | HELZ2 | -0.5033 | 30.864 | +0.04765 | +0.1899 | False |
| V_sig | S2_k100 | O1_Y1 | BSX | -0.4381 | 31.434 | +0.01608 | +0.2403 | False |
| V_sig | S1 | O1_Y1 | ZNF175 | +0.0049 | 32.284 | +0.00028 | +0.0095 | False |
| V_all | S1 | O1_Y1 | REPIN1 | +1.5303 | 32.284 | +0.01592 | +0.0443 | False |

## Stage 4 disagreement with directional readouts

All factors are ranked three ways in every setting: by delta (most negative first), by cos (highest first), and by drop in the frozen ruler score (largest drop first). Spearman correlations are between those rank vectors, so positive means agreement.

| variant | space | pair | rho_delta_cos | rho_delta_ruler | rho_cos_ruler | top20_overlap_delta_cos | top20_overlap_delta_ruler | top20_overlap_cos_ruler | top20_by_cos_not_toward_young |
|---|---|---|---|---|---|---|---|---|---|
| V_all | S1 | O1_Y1 | -0.008 | +0.094 | +0.293 | 1 | 0 | 6 | 20 |
| V_all | S1 | O2_Y1 | +0.974 | +0.321 | +0.318 | 14 | 3 | 2 | 1 |
| V_all | S2_k20 | O1_Y1 | +0.919 | +0.166 | +0.155 | 1 | 4 | 0 | 2 |
| V_all | S2_k20 | O1_Y2 | +0.974 | +0.325 | +0.310 | 7 | 3 | 2 | 20 |
| V_all | S2_k20 | O2_Y1 | +0.970 | +0.323 | +0.315 | 4 | 2 | 0 | 8 |
| V_all | S2_k20 | O2_Y2 | +0.967 | +0.309 | +0.308 | 5 | 2 | 0 | 6 |
| V_all | S2_k50 | O1_Y1 | +0.951 | +0.235 | +0.229 | 3 | 5 | 1 | 3 |
| V_all | S2_k50 | O1_Y2 | +0.972 | +0.352 | +0.348 | 6 | 2 | 1 | 20 |
| V_all | S2_k50 | O2_Y1 | +0.974 | +0.333 | +0.328 | 5 | 4 | 0 | 9 |
| V_all | S2_k50 | O2_Y2 | +0.973 | +0.326 | +0.324 | 5 | 3 | 0 | 9 |
| V_all | S2_k100 | O1_Y1 | +0.972 | +0.259 | +0.259 | 6 | 5 | 2 | 20 |
| V_all | S2_k100 | O1_Y2 | +0.943 | +0.360 | +0.372 | 10 | 2 | 1 | 20 |
| V_all | S2_k100 | O2_Y1 | +0.969 | +0.341 | +0.327 | 3 | 4 | 0 | 7 |
| V_all | S2_k100 | O2_Y2 | +0.973 | +0.336 | +0.324 | 3 | 4 | 0 | 8 |
| V_sig | S1 | O1_Y1 | +0.157 | +0.096 | +0.248 | 0 | 0 | 1 | 20 |
| V_sig | S1 | O2_Y1 | +0.881 | +0.097 | +0.101 | 15 | 0 | 0 | 20 |
| V_sig | S2_k20 | O1_Y1 | +0.871 | +0.073 | +0.067 | 9 | 2 | 1 | 20 |
| V_sig | S2_k20 | O1_Y2 | +0.923 | +0.084 | +0.074 | 7 | 1 | 0 | 20 |
| V_sig | S2_k20 | O2_Y1 | +0.947 | +0.103 | +0.098 | 1 | 2 | 0 | 20 |
| V_sig | S2_k20 | O2_Y2 | +0.948 | +0.097 | +0.094 | 0 | 3 | 0 | 20 |
| V_sig | S2_k50 | O1_Y1 | +0.880 | +0.157 | +0.153 | 9 | 2 | 1 | 20 |
| V_sig | S2_k50 | O1_Y2 | +0.922 | +0.143 | +0.140 | 9 | 3 | 1 | 20 |
| V_sig | S2_k50 | O2_Y1 | +0.952 | +0.129 | +0.133 | 4 | 2 | 0 | 20 |
| V_sig | S2_k50 | O2_Y2 | +0.952 | +0.125 | +0.128 | 3 | 2 | 0 | 20 |
| V_sig | S2_k100 | O1_Y1 | +0.874 | +0.190 | +0.190 | 9 | 2 | 2 | 20 |
| V_sig | S2_k100 | O1_Y2 | +0.922 | +0.173 | +0.172 | 9 | 4 | 1 | 20 |
| V_sig | S2_k100 | O2_Y1 | +0.947 | +0.117 | +0.119 | 5 | 2 | 0 | 20 |
| V_sig | S2_k100 | O2_Y2 | +0.949 | +0.114 | +0.116 | 6 | 2 | 0 | 20 |

The last column is the number asked for directly: of the 20 best factors by cosine, how many are **not** `toward_young` by delta.

## Fired keys

STOP keys fired: **none**.
STOP keys checked and not fired: `not_an_effect_matrix`, `too_few_factors`, `no_detection_power`.

REPORT-ONLY flags fired: `scale_unknown`, `coverage_low`.

- **scale_unknown** (report-only, not a stop): The units and log base of X are undeterminable from the file: uns is empty, X carries no unit attribute, var holds only mean/std/cv/fano/excess_cv with no unit statement, obs holds no unit statement, and the only bundled documentation in the folder is three plain gene-symbol lists. This is a REPORT-ONLY flag, not a stop. The Stage 2A/2B floors are magnitude-matched (they permute or re-use the same vector, so they are invariant to a global scale factor on X); the Stage 2C MDA is NOT invariant, because it scales the true young-minus-old direction to the median factor norm, and that norm moves with any global rescaling of X.

- **coverage_low** (report-only, not a stop): The chosen mapping (best_of_two_after_md3_idtype) matches 0.1228 of the frozen ruler's total absolute weight, which is below 0.25. The Stage 3 candidate table is titled exploratory and every candidate carries that word.

Recorded non-fatal failures and implementation decisions, verbatim:

- **Y2_singular:** Y2 at k=100: the young-donor covariance built on 56 fit donors has rank 55 < k=100, so the Mahalanobis metric is singular and undefined. No pseudo-inverse or ridge was substituted. Y2 is reported as undefined at this k; Y1 at the same k is unaffected.
- **stage2C_clearing_rule:** Implementation reading fixed before any floor existed: a dose 'clears the floor' when its delta is below the primary floor AND negative, mirroring the Stage 3 toward_young rule. Without the second clause, f = 0 (no displacement, delta exactly 0) would count as clearing whenever a floor happens to be positive, which would report detection power for adding nothing.

## Limitations

- **The inputs are regression estimates, not expression.** X holds the source's own per-gene effect estimates, so the frozen count pipeline (sum -> TMM log2-CPM with prior.count=2 -> frozen mu/sd -> z) was **not** applied and could not be. Only the frozen sd was used, to turn an effect into a displacement in frozen-z units. No library-size normalisation of this file was performed, because there are no counts here to normalise.
- **The effects are measured in one neonatal fibroblast line with no aged cells.** Every delta in this document is therefore a counterfactual: it assumes the perturbation effect measured in Hs27 transfers unchanged to an aged cell. Nothing in this data can test that assumption. The origin is an aged donor and the effect is not.
- **Unmeasured genes are assumed unmoved.** 18,578 of the 23,485 ruler coordinates (79.1%) are fixed at 0 for every factor. The chosen mapping covers 0.1228 of the ruler's absolute weight, which fired `coverage_low`; the candidate tables are titled exploratory for that reason. A factor that moves a gene this file does not measure is invisible here.
- **Guide counts are small.** Guides per factor run 2 to 6 (median 6), so the per-factor mean averages a handful of rows and the bootstrap over guides has very few distinct resamples. Factors with 2 guides get no bootstrap interval at all.
- **The log base of the input may be unresolved — and here it is.** `scale_unknown` fired. The magnitude-matched floors are invariant to a global scale factor on X, but the MDA is not, and neither is the absolute size of any delta or identity drop in frozen-z units. Read the deltas as relative to their own floor, not as absolute distances.
- **The identity panel is incomplete.** 9 of the 10 pre-registered identity genes are present (POSTN absent), so `identity_loss` is a weaker check than intended.
- **The S2 basis is not independent of the young target.** The PCA is fit on the GTEx fibroblast donor matrix, and Y1/Y2 are subsets of exactly that matrix. It is independent of the Southard data and of O2.
- **Y2 is undefined at k=100** (covariance rank 55 on 56 fit donors), and fragile wherever the fit donors number fewer than twice k. Reported as such; nothing was substituted.
- **O1 and O2 are not on the same footing.** O1 is a mean over GTEx donors in the same frozen-z space that the ruler was built in; O2 comes through SAME's own TMM panel of six pseudobulk rows and carries its own missing-gene zeroing (3,132 coordinates). Their distances to Y1 are reported separately and are never pooled.
- **O2's distance to Y1 is dominated by a cross-dataset offset, not by age.** In S1, ||O2 - Y1|| = 468.4 against ||O1 - Y1|| = 32.3, a ratio of 14.5x. At that distance almost any displacement with a slightly favourable projection lowers delta, so the O2 settings are the easiest place to earn a `toward_young` label and the least informative about ageing. The O2 candidate counts should be read with that in mind.
- **The q-value gate is resolution-limited.** With n_perm = 200, the smallest attainable p is 1/201, and BH over 1,836 factors needs at least 183 factors tied at that minimum before any q can reach 0.05. Settings with no q <= 0.05 are therefore not evidence of absence, and the candidate counts depend partly on how many factors tie at the p floor. More permutations would change them; that would be a new pre-registration.
- **The control-row floor is not a noise floor** and was never used as the gate: those rows are untargeted but not inert.
- This task makes **no claim about the source paper's own conclusions**. It did not test them.

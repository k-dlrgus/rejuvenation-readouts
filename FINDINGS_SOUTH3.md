**SUPERSEDED BY SOUTH4 (see FINDINGS_SOUTH4.md): SOUTH4 re-runs the ranking on the full measured per-cell singlets file named in SURVEY_TF.md before SOUTH2, not on the 4,914-gene mean-population regression subset. Nothing below is changed.**

# FINDINGS_SOUTH3 — rank Southard Hs27 CRISPRa transcription factors by approach to a young target; SOUTH2 re-run with a resolvable null, one guide-agreement rule, and report-only displacement geometry

**Status:** Ran to completion. No STOP key fired. Report-only flags: `scale_unknown`, `coverage_low`.

**Result:** `toward_young` is 0 in 2 of the 28 settings and non-zero in 26, ranging 40-277 of 1836 factors. Non-zero settings are V_all, V_sig. The counts are not pooled and no single setting is primary. **No.** Across the 26 settings that have survivors, the Spearman correlation between the delta ranking and the cosine ranking over all 1,836 factors is at least +0.871, and within every quartile of ||d|| / ||target - origin|| in those settings it is at least +0.878. Taking the top-n factors by cosine alone, n being that setting's number of survivors, recovers 2,329 of the 3,326 survivors (70.0%) — and taking the top-n by delta itself, the statistic the label is built on, recovers 2,334 (70.2%). The two agree to within a percentage point, so the ~30% that neither recovers is the conjunction of the floor and q gates, not information the cosine lacks. ||d|| alone recovers only 542 (16.3%), so the label is not merely picking the biggest displacements. The mechanism is the identity: in every setting that has a survivor the median residual of delta = -||d||cos(theta) is at most 5.6% of ||d||, so delta there is the cosine re-weighted by the displacement's own length and carries nothing beyond it. No survivor in any setting is distinguishable on direction from what a plain cosine ranking would give. What the delta machinery does buy is not a different ordering but a bar: it says which of those cosines are larger than the factor's own reshuffled null, and a cosine ranking on its own never says that. Read the "What produced these counts" section before any count: some settings are zero for geometric reasons rather than because factors were tested and found null, and the O2 settings are dominated by a cross-dataset offset.

**SOUTH3 supersedes SOUTH2.** It is a patch, not a new analysis: same pinned file, same pre-registration, same spaces, variants, origins, targets, floors, labelling rule, what-does-not-count list, exploratory titling under `coverage_low`, and same seeds. Exactly three things change: (1) `n_perm` = 20,000 for the Stage 2A per-factor permutation null instead of 200, so BH across 1,836 factors can resolve q; (2) guide agreement is computed on the V_all vectors always and that one value is applied to both gene variants; (3) a report-only section on ||d||, ||d|| / ||target - origin||, rank agreement between delta and cosine as a function of that ratio, and the small-displacement identity. `FINDINGS_SOUTH2.md` was not edited except for one line at its top marking it superseded.

SOUTH2 itself superseded SOUTH, which halted at Stage 0 with `no_raw_counts` because its pre-registration described this file as a count matrix; the file holds per-guide regression output.

The SOUTH3 patch below was written verbatim to `results/south3/PREREG.flag`, and the inherited SOUTH2 block to `results/south3/PREREG_SOUTH2_inherited.flag`, before any SOUTH3 statistic was computed. Nothing was re-tuned after numbers existed.

## Pre-registration — the SOUTH3 patch (verbatim, written before any SOUTH3 statistic)

```
TASK: SOUTH3 — patch of SOUTH2. Same data, same prereg, three changes.

Write FINDINGS_SOUTH3.md, PROGRESS_SOUTH3.md, results/south3/*, src/south3_run.py.
Add one line to the top of FINDINGS_SOUTH2.md marking it superseded; change
nothing else in it. Touch results/south3/PREREG.flag with this file verbatim
before any statistic. Seeds unchanged.

1. n_perm = 20000 for the Stage 2A per-factor permutation null, so BH across
   1,836 factors can resolve q. Report the attainable minimum p. Everything
   else about the null is unchanged.
2. Guide agreement is computed on V_all vectors ALWAYS, and that one value is
   applied to both variants. Report how many factors change label versus SOUTH2
   because of this.
3. Report-only, no gate: for each factor report ||d|| and the ratio
   ||d|| / ||target - origin||, and the Spearman correlation between delta rank
   and cosine rank as a function of that ratio (quartiles). State the small-
   displacement identity delta ≈ -||d||cos(theta) and whether the data follow it.

Nothing else changes: same spaces, variants, origins, targets, floors, labels,
what-does-not-count list, and exploratory titling under coverage_low.
Report the final toward_young counts per setting, and state plainly whether any
survivor is distinguishable from the ranking a plain cosine would give.
```

## Pre-registration — the inherited SOUTH2 block (verbatim, unchanged)

Everything in this block still governs SOUTH3 except the three patched items above.

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

**CHANGE 1 lives here and nowhere else.** `n_perm` is 20,000 instead of SOUTH2's 200. The permuted coordinate set, the per-factor generator (`np.random.default_rng([20260914, factor_index])`), the delta definition, the 5th percentile rule, and the pooled-floor construction are all unchanged. The smallest p an empirical one-sided permutation test can now return is 1/(20,000+1) = **5.000e-05**, against 1/201 = 0.00498 in SOUTH2. With BH across 1,836 factors, q <= 0.05 now needs only 2 factors tied at that minimum p, against 183 in SOUTH2. That is the whole point of the change: the BH gate can resolve q instead of being pinned above it by the resolution of the null.

For each factor, its own vector's entries are permuted across the 4,907 matched gene coordinates 20,000 times (seed 20260914, an independent generator per factor), and `delta = ||(origin + d_perm) - target|| - ||origin - target||` is recomputed each time. The per-factor floor is the 5th percentile of that distribution. Permuting preserves ||d|| exactly, so the floor is magnitude-matched by construction.

The pooled floor pools every factor's permuted vectors after rescaling each to the median factor norm, and takes the 5th percentile of that pooled distribution. It is the bar the Stage 2C dose-response has to clear.

**Implementation decision, stated because the pre-registration says only "across genes":** the permutation is over the 4,907 measured coordinates, not over all 23,485 ruler coordinates. The other 18,578 are structurally fixed at 0 by the conversion rule, so permuting effects into them would place the null in a subspace the experiment cannot report on and would no longer be a like-for-like random direction for the real effect. This was fixed before any floor was computed and was not revisited.

Distribution of the per-factor floor across factors:

| variant | space | pair | floor_min | floor_p25 | floor_median | floor_p75 | floor_max | n_factors |
|---|---|---|---|---|---|---|---|---|
| V_all | S1 | O1_Y1 | +1.32162 | +4.18527 | +5.16423 | +6.73916 | +88.63016 | 1,836 |
| V_all | S1 | O2_Y1 | -1.22246 | -0.49689 | -0.31389 | -0.13702 | +11.67097 | 1,836 |
| V_all | S2_k20 | O1_Y1 | -2.36318 | -0.97435 | -0.76671 | -0.56858 | +0.25672 | 1,836 |
| V_all | S2_k20 | O1_Y2 | -1.79525 | -0.62472 | -0.52019 | -0.41699 | +0.05300 | 1,836 |
| V_all | S2_k20 | O2_Y1 | -1.51797 | -0.35791 | -0.30976 | -0.27048 | -0.14679 | 1,836 |
| V_all | S2_k20 | O2_Y2 | -1.43114 | -0.31649 | -0.27742 | -0.24624 | -0.15275 | 1,836 |
| V_all | S2_k50 | O1_Y1 | -1.88489 | -0.86446 | -0.69175 | -0.52351 | +0.23206 | 1,836 |
| V_all | S2_k50 | O1_Y2 | -1.40810 | -0.56231 | -0.47844 | -0.38981 | +0.05271 | 1,836 |
| V_all | S2_k50 | O2_Y1 | -1.65061 | -0.42387 | -0.35886 | -0.30814 | -0.12285 | 1,836 |
| V_all | S2_k50 | O2_Y2 | -1.58758 | -0.38991 | -0.33433 | -0.29118 | -0.14102 | 1,836 |
| V_all | S2_k100 | O1_Y1 | -1.35299 | -0.77796 | -0.63822 | -0.49123 | +0.24408 | 1,836 |
| V_all | S2_k100 | O1_Y2 | -0.97725 | -0.49954 | -0.43863 | -0.36712 | +0.07063 | 1,836 |
| V_all | S2_k100 | O2_Y1 | -1.96119 | -0.59455 | -0.48908 | -0.39280 | -0.02424 | 1,836 |
| V_all | S2_k100 | O2_Y2 | -1.89058 | -0.56285 | -0.46339 | -0.37720 | -0.04214 | 1,836 |
| V_sig | S1 | O1_Y1 | +0.00398 | +0.09196 | +0.16418 | +0.30468 | +13.08838 | 1,836 |
| V_sig | S1 | O2_Y1 | -0.17129 | -0.00388 | +0.00331 | +0.01496 | +0.90424 | 1,836 |
| V_sig | S2_k20 | O1_Y1 | -0.47855 | -0.00665 | -0.00088 | +0.00761 | +0.41364 | 1,836 |
| V_sig | S2_k20 | O1_Y2 | -0.40656 | -0.02267 | -0.01705 | -0.01307 | +0.14641 | 1,836 |
| V_sig | S2_k20 | O2_Y1 | -0.35337 | -0.04634 | -0.03427 | -0.02653 | -0.00839 | 1,836 |
| V_sig | S2_k20 | O2_Y2 | -0.36741 | -0.04935 | -0.03647 | -0.02817 | -0.00862 | 1,836 |
| V_sig | S2_k50 | O1_Y1 | -0.44820 | -0.00829 | -0.00333 | +0.00405 | +0.38455 | 1,836 |
| V_sig | S2_k50 | O1_Y2 | -0.38214 | -0.02275 | -0.01733 | -0.01348 | +0.14827 | 1,836 |
| V_sig | S2_k50 | O2_Y1 | -0.35121 | -0.04337 | -0.03236 | -0.02536 | -0.00826 | 1,836 |
| V_sig | S2_k50 | O2_Y2 | -0.35243 | -0.04569 | -0.03408 | -0.02653 | -0.00843 | 1,836 |
| V_sig | S2_k100 | O1_Y1 | -0.42116 | -0.00820 | -0.00339 | +0.00388 | +0.40208 | 1,836 |
| V_sig | S2_k100 | O1_Y2 | -0.36169 | -0.02191 | -0.01671 | -0.01288 | +0.17402 | 1,836 |
| V_sig | S2_k100 | O2_Y1 | -0.39715 | -0.03011 | -0.02292 | -0.01823 | +0.03998 | 1,836 |
| V_sig | S2_k100 | O2_Y2 | -0.38865 | -0.03246 | -0.02479 | -0.01960 | +0.01583 | 1,836 |

Pooled floor, built from a vector at the median factor norm:

| variant | space | pair | median_factor_norm | pooled_floor_p05 | null_median | null_min | null_max |
|---|---|---|---|---|---|---|---|
| V_all | S1 | O1_Y1 | 20.2543 | +4.97706 | +5.39106 | +4.21771 | +6.59353 |
| V_all | S1 | O2_Y1 | 20.2543 | -0.58202 | -0.08058 | -1.53530 | +1.49957 |
| V_all | S2_k20 | O1_Y1 | 20.2543 | -0.98582 | -0.53456 | -1.73709 | +1.04338 |
| V_all | S2_k20 | O1_Y2 | 20.2543 | -0.57580 | -0.27516 | -1.20730 | +0.90887 |
| V_all | S2_k20 | O2_Y1 | 20.2543 | -0.31319 | -0.07205 | -0.84563 | +0.72217 |
| V_all | S2_k20 | O2_Y2 | 20.2543 | -0.27774 | -0.03990 | -0.80231 | +0.73220 |
| V_all | S2_k50 | O1_Y1 | 20.2543 | -0.85011 | -0.46158 | -1.50234 | +1.02122 |
| V_all | S2_k50 | O1_Y2 | 20.2543 | -0.51073 | -0.23537 | -1.12189 | +0.90406 |
| V_all | S2_k50 | O2_Y1 | 20.2543 | -0.37197 | -0.11420 | -0.94616 | +0.78549 |
| V_all | S2_k50 | O2_Y2 | 20.2543 | -0.34164 | -0.08996 | -0.91821 | +0.78710 |
| V_all | S2_k100 | O1_Y1 | 20.2543 | -0.74075 | -0.40793 | -1.34560 | +1.01847 |
| V_all | S2_k100 | O1_Y2 | 20.2543 | -0.44990 | -0.19625 | -1.00840 | +0.89786 |
| V_all | S2_k100 | O2_Y1 | 20.2543 | -0.55137 | -0.24265 | -1.20092 | +0.84773 |
| V_all | S2_k100 | O2_Y2 | 20.2543 | -0.51470 | -0.21814 | -1.14655 | +0.84901 |
| V_sig | S1 | O1_Y1 | 3.2889 | +0.15738 | +0.19972 | +0.03514 | +0.33772 |
| V_sig | S1 | O2_Y1 | 3.2889 | -0.00201 | +0.04408 | -0.14123 | +0.21747 |
| V_sig | S2_k20 | O1_Y1 | 3.2889 | -0.00776 | +0.03789 | -0.14538 | +0.19632 |
| V_sig | S2_k20 | O1_Y2 | 3.2889 | -0.02133 | +0.02124 | -0.14984 | +0.16595 |
| V_sig | S2_k20 | O2_Y1 | 3.2889 | -0.03460 | +0.00550 | -0.12854 | +0.13262 |
| V_sig | S2_k20 | O2_Y2 | 3.2889 | -0.03664 | +0.00344 | -0.13299 | +0.13058 |
| V_sig | S2_k50 | O1_Y1 | 3.2889 | -0.00952 | +0.03493 | -0.14034 | +0.18158 |
| V_sig | S2_k50 | O1_Y2 | 3.2889 | -0.02133 | +0.02025 | -0.14594 | +0.16507 |
| V_sig | S2_k50 | O2_Y1 | 3.2889 | -0.03305 | +0.00803 | -0.12981 | +0.15375 |
| V_sig | S2_k50 | O2_Y2 | 3.2889 | -0.03443 | +0.00649 | -0.13359 | +0.15245 |
| V_sig | S2_k100 | O1_Y1 | 3.2889 | -0.00937 | +0.03416 | -0.13341 | +0.17792 |
| V_sig | S2_k100 | O1_Y2 | 3.2889 | -0.02072 | +0.02024 | -0.13691 | +0.16288 |
| V_sig | S2_k100 | O2_Y1 | 3.2889 | -0.02572 | +0.01634 | -0.13354 | +0.16577 |
| V_sig | S2_k100 | O2_Y2 | 3.2889 | -0.02701 | +0.01478 | -0.13222 | +0.16108 |

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
| V_all | S1 | O1_Y1 | full_v | 20.2543 | +4.97706 | 0.050 | True |
| V_all | S1 | O2_Y1 | full_v | 20.2543 | -0.58202 | 0.250 | True |
| V_all | S2_k20 | O1_Y1 | full_v | 20.2543 | -0.98582 | 0.100 | True |
| V_all | S2_k20 | O1_Y2 | full_v | 20.2543 | -0.57580 | 0.050 | True |
| V_all | S2_k20 | O2_Y1 | full_v | 20.2543 | -0.31319 | 0.100 | True |
| V_all | S2_k20 | O2_Y2 | full_v | 20.2543 | -0.27774 | 0.100 | True |
| V_all | S2_k50 | O1_Y1 | full_v | 20.2543 | -0.85011 | 0.050 | True |
| V_all | S2_k50 | O1_Y2 | full_v | 20.2543 | -0.51073 | 0.050 | True |
| V_all | S2_k50 | O2_Y1 | full_v | 20.2543 | -0.37197 | 0.100 | True |
| V_all | S2_k50 | O2_Y2 | full_v | 20.2543 | -0.34164 | 0.100 | True |
| V_all | S2_k100 | O1_Y1 | full_v | 20.2543 | -0.74075 | 0.050 | True |
| V_all | S2_k100 | O1_Y2 | full_v | 20.2543 | -0.44990 | 0.050 | True |
| V_all | S2_k100 | O2_Y1 | full_v | 20.2543 | -0.55137 | 0.250 | True |
| V_all | S2_k100 | O2_Y2 | full_v | 20.2543 | -0.51470 | 0.250 | True |
| V_all | S1 | O1_Y1 | measured_coords_only | 20.2543 | +4.97706 | 0.050 | True |
| V_all | S1 | O2_Y1 | measured_coords_only | 20.2543 | -0.58202 | 0.500 | True |
| V_all | S2_k20 | O1_Y1 | measured_coords_only | 20.2543 | -0.98582 | 0.250 | True |
| V_all | S2_k20 | O1_Y2 | measured_coords_only | 20.2543 | -0.57580 | 0.100 | True |
| V_all | S2_k20 | O2_Y1 | measured_coords_only | 20.2543 | -0.31319 | 0.100 | True |
| V_all | S2_k20 | O2_Y2 | measured_coords_only | 20.2543 | -0.27774 | 0.100 | True |
| V_all | S2_k50 | O1_Y1 | measured_coords_only | 20.2543 | -0.85011 | 0.100 | True |
| V_all | S2_k50 | O1_Y2 | measured_coords_only | 20.2543 | -0.51073 | 0.100 | True |
| V_all | S2_k50 | O2_Y1 | measured_coords_only | 20.2543 | -0.37197 | 0.250 | True |
| V_all | S2_k50 | O2_Y2 | measured_coords_only | 20.2543 | -0.34164 | 0.100 | True |
| V_all | S2_k100 | O1_Y1 | measured_coords_only | 20.2543 | -0.74075 | 0.100 | True |
| V_all | S2_k100 | O1_Y2 | measured_coords_only | 20.2543 | -0.44990 | 0.100 | True |
| V_all | S2_k100 | O2_Y1 | measured_coords_only | 20.2543 | -0.55137 | 0.250 | True |
| V_all | S2_k100 | O2_Y2 | measured_coords_only | 20.2543 | -0.51470 | 0.250 | True |
| V_sig | S1 | O1_Y1 | full_v | 3.2889 | +0.15738 | 0.050 | True |
| V_sig | S1 | O2_Y1 | full_v | 3.2889 | -0.00201 | 0.050 | True |
| V_sig | S2_k20 | O1_Y1 | full_v | 3.2889 | -0.00776 | 0.050 | True |
| V_sig | S2_k20 | O1_Y2 | full_v | 3.2889 | -0.02133 | 0.050 | True |
| V_sig | S2_k20 | O2_Y1 | full_v | 3.2889 | -0.03460 | 0.050 | True |
| V_sig | S2_k20 | O2_Y2 | full_v | 3.2889 | -0.03664 | 0.050 | True |
| V_sig | S2_k50 | O1_Y1 | full_v | 3.2889 | -0.00952 | 0.050 | True |
| V_sig | S2_k50 | O1_Y2 | full_v | 3.2889 | -0.02133 | 0.050 | True |
| V_sig | S2_k50 | O2_Y1 | full_v | 3.2889 | -0.03305 | 0.050 | True |
| V_sig | S2_k50 | O2_Y2 | full_v | 3.2889 | -0.03443 | 0.050 | True |
| V_sig | S2_k100 | O1_Y1 | full_v | 3.2889 | -0.00937 | 0.050 | True |
| V_sig | S2_k100 | O1_Y2 | full_v | 3.2889 | -0.02072 | 0.050 | True |
| V_sig | S2_k100 | O2_Y1 | full_v | 3.2889 | -0.02572 | 0.050 | True |
| V_sig | S2_k100 | O2_Y2 | full_v | 3.2889 | -0.02701 | 0.050 | True |
| V_sig | S1 | O1_Y1 | measured_coords_only | 3.2889 | +0.15738 | 0.050 | True |
| V_sig | S1 | O2_Y1 | measured_coords_only | 3.2889 | -0.00201 | 0.050 | True |
| V_sig | S2_k20 | O1_Y1 | measured_coords_only | 3.2889 | -0.00776 | 0.050 | True |
| V_sig | S2_k20 | O1_Y2 | measured_coords_only | 3.2889 | -0.02133 | 0.050 | True |
| V_sig | S2_k20 | O2_Y1 | measured_coords_only | 3.2889 | -0.03460 | 0.100 | True |
| V_sig | S2_k20 | O2_Y2 | measured_coords_only | 3.2889 | -0.03664 | 0.100 | True |
| V_sig | S2_k50 | O1_Y1 | measured_coords_only | 3.2889 | -0.00952 | 0.050 | True |
| V_sig | S2_k50 | O1_Y2 | measured_coords_only | 3.2889 | -0.02133 | 0.050 | True |
| V_sig | S2_k50 | O2_Y1 | measured_coords_only | 3.2889 | -0.03305 | 0.100 | True |
| V_sig | S2_k50 | O2_Y2 | measured_coords_only | 3.2889 | -0.03443 | 0.100 | True |
| V_sig | S2_k100 | O1_Y1 | measured_coords_only | 3.2889 | -0.00937 | 0.050 | True |
| V_sig | S2_k100 | O1_Y2 | measured_coords_only | 3.2889 | -0.02072 | 0.050 | True |
| V_sig | S2_k100 | O2_Y1 | measured_coords_only | 3.2889 | -0.02572 | 0.050 | True |
| V_sig | S2_k100 | O2_Y2 | measured_coords_only | 3.2889 | -0.02701 | 0.050 | True |

Full dose-response deltas at every f are in `results/south3/stage2c_dose_response.csv`.

- `no_detection_power`: **not fired**. Some f <= 1.0 clears the primary floor, so Stage 3 runs.
- No floor is ever subtracted from a later number. It is a bar to clear.

## Stage 3 ranking

The full table is `results/south3/stage3_full_table.csv`: **51,408 rows** = 1,836 factors x 2 gene variants x 14 space/origin/target settings, with cos, frac, delta, ||d||, ||v||, the factor's own Stage 2A floor, margin over that floor, permutation p, BH q, bootstrap intervals, the Mahalanobis change and membership test for Y2, guide agreement, the identity panel, and the frozen ruler score change. Nothing is pooled across spaces, variants, origins, or targets.

- **Guide reliability — CHANGE 2 lives here.** Each guide's own vector is built and the median pairwise cosine between guides of the same factor is computed **on the V_all vectors always**, and that single number per factor is applied to both gene variants. `guides_disagree` fires when it is <= 0. Those factors are reported and are never labelled a candidate. The V_sig median pairwise cosine is still computed and carried in the full table as `guide_median_cos_vsig_reported_only`; it gates nothing. SOUTH2 gated each variant on its own value, which is what this change replaces.
- **Identity check**: the panel actually used is COL1A1, COL1A2, FN1, LUM, PDGFRA, PDGFRB, PRRX1, SERPINH1, VIM — **9 of the 10 pre-registered genes**. Absent: POSTN. `identity_drop` is the negated mean displacement over that panel **in frozen-z units** (the analysis's own coordinates), and `identity_loss` fires above 0.5.
- **Frozen ruler score change** is a context column only. It gates nothing.
- **Nulls**: the Stage 2A per-factor permutations, 20,000 per factor and setting, one-sided empirical p = (count(perm delta <= observed delta) + 1) / (20,000 + 1), so the attainable minimum is **p = 5.000e-05**; then BH-FDR across factors **within each space and variant** (and within each origin/target pair, since a pair is part of the setting).
- **Report-only geometry (CHANGE 3)**: the full table also carries `d_norm`, `d_over_v` = ||d|| / ||target - origin||, and `delta_small_disp_approx` = -||d||cos(theta). None of the three gates anything; see the geometry section below.
- **Bootstrap**: 200 resamples over guides (seed 20260918) for factors with >= 3 guides. Intervals that do not contain their own point estimate are flagged INVALID and are never a gate: **618** of 51,380 computed delta intervals are so flagged.

A factor is labelled `toward_young` only if delta is negative AND below its own Stage 2A floor AND q <= 0.05 AND `identity_loss` is not flagged AND `guides_disagree` is not flagged.

### Candidate counts in every space, variant, and origin/target pair

Zero is a valid and reportable answer.

| variant | space | pair | n_factors | delta_negative | below_own_floor | at_p_floor | q_le_005 | guides_disagree | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| V_all | S1 | O1_Y1 | 1,836 | 0 | 578 | 140 | 436 | 299 | 0 | 0 |
| V_all | S1 | O2_Y1 | 1,836 | 918 | 477 | 139 | 325 | 299 | 0 | 266 |
| V_all | S2_k20 | O1_Y1 | 1,836 | 1,645 | 498 | 93 | 322 | 299 | 0 | 277 |
| V_all | S2_k20 | O1_Y2 | 1,836 | 1,209 | 194 | 17 | 50 | 299 | 0 | 44 |
| V_all | S2_k20 | O2_Y1 | 1,836 | 849 | 384 | 117 | 260 | 299 | 0 | 219 |
| V_all | S2_k20 | O2_Y2 | 1,836 | 775 | 372 | 114 | 246 | 299 | 0 | 208 |
| V_all | S2_k50 | O1_Y1 | 1,836 | 1,595 | 449 | 73 | 283 | 299 | 0 | 247 |
| V_all | S2_k50 | O1_Y2 | 1,836 | 1,111 | 194 | 18 | 61 | 299 | 0 | 53 |
| V_all | S2_k50 | O2_Y1 | 1,836 | 913 | 398 | 109 | 288 | 299 | 0 | 236 |
| V_all | S2_k50 | O2_Y2 | 1,836 | 840 | 383 | 113 | 280 | 299 | 0 | 230 |
| V_all | S2_k100 | O1_Y1 | 1,836 | 1,545 | 366 | 50 | 198 | 299 | 0 | 173 |
| V_all | S2_k100 | O1_Y2 | 1,836 | 957 | 160 | 12 | 46 | 299 | 0 | 40 |
| V_all | S2_k100 | O2_Y1 | 1,836 | 1,124 | 421 | 121 | 301 | 299 | 0 | 248 |
| V_all | S2_k100 | O2_Y2 | 1,836 | 1,057 | 405 | 123 | 293 | 299 | 0 | 243 |
| V_sig | S1 | O1_Y1 | 1,836 | 0 | 317 | 22 | 56 | 299 | 0 | 0 |
| V_sig | S1 | O2_Y1 | 1,836 | 394 | 406 | 45 | 159 | 299 | 0 | 128 |
| V_sig | S2_k20 | O1_Y1 | 1,836 | 398 | 373 | 21 | 74 | 299 | 0 | 67 |
| V_sig | S2_k20 | O1_Y2 | 1,836 | 571 | 238 | 21 | 49 | 299 | 0 | 45 |
| V_sig | S2_k20 | O2_Y1 | 1,836 | 903 | 257 | 16 | 44 | 299 | 0 | 40 |
| V_sig | S2_k20 | O2_Y2 | 1,836 | 918 | 252 | 19 | 50 | 299 | 0 | 45 |
| V_sig | S2_k50 | O1_Y1 | 1,836 | 411 | 332 | 20 | 69 | 299 | 0 | 62 |
| V_sig | S2_k50 | O1_Y2 | 1,836 | 552 | 225 | 20 | 49 | 299 | 0 | 43 |
| V_sig | S2_k50 | O2_Y1 | 1,836 | 972 | 327 | 22 | 86 | 299 | 0 | 79 |
| V_sig | S2_k50 | O2_Y2 | 1,836 | 983 | 317 | 25 | 92 | 299 | 0 | 80 |
| V_sig | S2_k100 | O1_Y1 | 1,836 | 396 | 310 | 19 | 52 | 299 | 0 | 48 |
| V_sig | S2_k100 | O1_Y2 | 1,836 | 530 | 217 | 15 | 45 | 299 | 0 | 40 |
| V_sig | S2_k100 | O2_Y1 | 1,836 | 786 | 326 | 27 | 92 | 299 | 0 | 80 |
| V_sig | S2_k100 | O2_Y2 | 1,836 | 806 | 325 | 28 | 99 | 299 | 0 | 85 |

Y2 Mahalanobis columns, reported and never a gate. The Y2 `delta` above is the pre-registered Euclidean delta to the Y2 region's centre; the Mahalanobis change has its own floor from the same Stage 2A permutations. Where Y2 is undefined (see the radius table) these counts are 0 because nothing could be computed, not because nothing moved.

| variant | space | pair | md_change_below_own_floor | md_q_le_005 | origin_plus_d_inside_radius |
|---|---|---|---|---|---|
| V_all | S2_k20 | O1_Y2 | 104 | 22 | 1,836 |
| V_all | S2_k20 | O2_Y2 | 295 | 132 | 0 |
| V_all | S2_k50 | O1_Y2 | 221 | 53 | 1,836 |
| V_all | S2_k50 | O2_Y2 | 544 | 405 | 0 |
| V_all | S2_k100 | O1_Y2 | 0 | 0 | 0 |
| V_all | S2_k100 | O2_Y2 | 0 | 0 | 0 |
| V_sig | S2_k20 | O1_Y2 | 181 | 44 | 1,836 |
| V_sig | S2_k20 | O2_Y2 | 211 | 16 | 0 |
| V_sig | S2_k50 | O1_Y2 | 224 | 48 | 1,836 |
| V_sig | S2_k50 | O2_Y2 | 357 | 121 | 0 |
| V_sig | S2_k100 | O1_Y2 | 0 | 0 | 0 |
| V_sig | S2_k100 | O2_Y2 | 0 | 0 | 0 |

### Top 20 by delta in every setting (exploratory candidate tables)

`coverage_low` fired in Stage 0, so these tables are **exploratory** and every factor labelled `toward_young` carries that word in its name below.

**V_all / S1 / O1_Y1** — 0 of these 20 are `toward_young`; 0 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| REPIN1 | 6 | +1.53031 | +1.43477 | -0.09554 | 0.7239 | +0.034 | False | +0.019 | False | False |
| ZBTB46 | 6 | +1.53040 | +1.32162 | -0.20879 | 1.0000 | +0.019 | False | +0.015 | False | False |
| LHX3 | 6 | +1.77375 | +1.43623 | -0.33752 | 1.0000 | +0.033 | False | +0.009 | False | False |
| RXRA | 6 | +1.80255 | +1.74227 | -0.06029 | 0.4366 | +0.040 | False | +0.016 | False | False |
| CAPN15 | 6 | +2.06517 | +1.74485 | -0.32032 | 1.0000 | +0.001 | False | +0.035 | False | False |
| MAFG | 6 | +2.08524 | +2.51233 | +0.42710 | 0.0007 | +0.039 | False | +0.015 | False | False |
| ZIC1 | 6 | +2.15891 | +1.95673 | -0.20219 | 1.0000 | +0.035 | False | +0.085 | False | False |
| ZNF770 | 6 | +2.17700 | +2.00386 | -0.17314 | 1.0000 | -0.002 | True | +0.023 | False | False |
| MIER2 | 6 | +2.31056 | +2.32301 | +0.01245 | 0.1232 | +0.015 | False | -0.029 | False | False |
| HEYL | 6 | +2.35534 | +2.31627 | -0.03906 | 0.3015 | +0.021 | False | -0.020 | False | False |
| STAT5B | 6 | +2.36037 | +2.35198 | -0.00839 | 0.1833 | +0.006 | False | +0.052 | False | False |
| ZBTB8B | 6 | +2.37583 | +2.21854 | -0.15729 | 0.9479 | +0.002 | False | +0.014 | False | False |
| GZF1 | 6 | +2.38620 | +2.23672 | -0.14948 | 0.9095 | +0.020 | False | +0.011 | False | False |
| ZNF211 | 6 | +2.44489 | +2.60922 | +0.16433 | 0.0035 | +0.001 | False | +0.030 | False | False |
| LZTR1 | 6 | +2.46532 | +2.26685 | -0.19848 | 1.0000 | +0.015 | False | +0.054 | False | False |
| MAPK8IP1 | 6 | +2.51469 | +2.37512 | -0.13958 | 0.8421 | +0.008 | False | +0.066 | False | False |
| ZMIZ2 | 6 | +2.53936 | +2.38102 | -0.15834 | 0.9414 | +0.006 | False | -0.006 | False | False |
| TONSL | 6 | +2.55780 | +2.58346 | +0.02566 | 0.0938 | -0.007 | True | -0.000 | False | False |
| ZNF496 | 6 | +2.56473 | +2.39967 | -0.16506 | 0.9587 | +0.004 | False | -0.044 | False | False |
| TADA2B | 6 | +2.56560 | +2.71096 | +0.14537 | 0.0045 | +0.020 | False | +0.046 | False | False |

**V_all / S1 / O2_Y1** — 20 of these 20 are `toward_young`; 266 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -2.65651 | +0.58696 | +3.24347 | 0.0007 | +0.324 | False | -0.036 | False | True |
| PRDM16 (exploratory) | 6 | -1.81109 | -1.08341 | +0.72768 | 0.0007 | +0.068 | False | +0.012 | False | True |
| XRCC5 (exploratory) | 6 | -1.73822 | -0.69599 | +1.04224 | 0.0007 | +0.022 | False | +0.221 | False | True |
| ZNF208 (exploratory) | 4 | -1.67677 | +0.14702 | +1.82379 | 0.0007 | +0.019 | False | -0.109 | False | True |
| PLAG1 (exploratory) | 6 | -1.57243 | -0.48367 | +1.08876 | 0.0007 | +0.057 | False | -0.160 | False | True |
| HAND1 (exploratory) | 6 | -1.51582 | -0.24406 | +1.27176 | 0.0007 | +0.108 | False | +0.054 | False | True |
| HDAC1 (exploratory) | 6 | -1.49956 | -0.66710 | +0.83246 | 0.0007 | +0.019 | False | +0.059 | False | True |
| KLF5 (exploratory) | 6 | -1.49952 | -0.83436 | +0.66515 | 0.0007 | +0.166 | False | +0.266 | False | True |
| SPDEF (exploratory) | 6 | -1.40700 | -0.52156 | +0.88544 | 0.0007 | +0.030 | False | +0.107 | False | True |
| ZNF470 (exploratory) | 6 | -1.32587 | -0.52373 | +0.80214 | 0.0007 | +0.006 | False | +0.054 | False | True |
| ZNF615 (exploratory) | 6 | -1.31116 | -0.75086 | +0.56030 | 0.0007 | +0.026 | False | +0.124 | False | True |
| ETV1 (exploratory) | 6 | -1.31019 | -0.98048 | +0.32972 | 0.0028 | +0.061 | False | -0.086 | False | True |
| TFAP4 (exploratory) | 6 | -1.30706 | -0.60621 | +0.70085 | 0.0007 | +0.049 | False | +0.051 | False | True |
| TSHZ1 (exploratory) | 6 | -1.29290 | -0.72768 | +0.56522 | 0.0007 | +0.055 | False | +0.095 | False | True |
| LBX1 (exploratory) | 6 | -1.25818 | -0.95088 | +0.30730 | 0.0074 | +0.044 | False | +0.142 | False | True |
| KLF2 (exploratory) | 6 | -1.25818 | -0.79033 | +0.46785 | 0.0007 | +0.044 | False | +0.090 | False | True |
| ZNF230 (exploratory) | 6 | -1.25248 | -1.02169 | +0.23078 | 0.0087 | +0.028 | False | +0.142 | False | True |
| ENO1 (exploratory) | 6 | -1.25080 | -0.58738 | +0.66342 | 0.0007 | +0.021 | False | -0.120 | False | True |
| ZSCAN20 (exploratory) | 6 | -1.24636 | -0.62056 | +0.62580 | 0.0007 | +0.024 | False | +0.039 | False | True |
| HHEX (exploratory) | 6 | -1.24307 | -1.04506 | +0.19802 | 0.0363 | +0.044 | False | +0.282 | False | True |

**V_all / S2_k20 / O1_Y1** — 19 of these 20 are `toward_young`; 277 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| NCL (exploratory) | 6 | -2.73849 | -1.49381 | +1.24468 | 0.0010 | +0.020 | False | +0.404 | False | True |
| PRDM1 (exploratory) | 6 | -2.59509 | -1.68217 | +0.91292 | 0.0010 | +0.099 | False | +0.414 | False | True |
| DMRT1 (exploratory) | 6 | -2.40814 | -1.57270 | +0.83544 | 0.0010 | +0.128 | False | +0.328 | False | True |
| HOXA13 (exploratory) | 6 | -2.33688 | -1.32427 | +1.01261 | 0.0010 | +0.068 | False | +0.296 | False | True |
| LBX1 (exploratory) | 6 | -2.20883 | -1.55975 | +0.64908 | 0.0010 | +0.044 | False | +0.142 | False | True |
| PLAG1 (exploratory) | 6 | -2.14772 | -1.52957 | +0.61815 | 0.0010 | +0.057 | False | -0.160 | False | True |
| CEBPB (exploratory) | 6 | -2.14483 | -1.48673 | +0.65811 | 0.0010 | +0.091 | False | -0.233 | False | True |
| NFX1 (exploratory) | 6 | -2.10116 | -1.50186 | +0.59930 | 0.0010 | +0.022 | False | +0.180 | False | True |
| ARID2 | 6 | -2.06375 | -1.98450 | +0.07925 | 0.1005 | +0.045 | False | +0.288 | False | False |
| ZNF850 (exploratory) | 6 | -2.02629 | -1.55666 | +0.46963 | 0.0016 | +0.029 | False | +0.260 | False | True |
| PSMD14 (exploratory) | 6 | -1.99307 | -1.34212 | +0.65096 | 0.0010 | +0.031 | False | +0.125 | False | True |
| CEBPD (exploratory) | 6 | -1.99152 | -1.18994 | +0.80157 | 0.0010 | +0.043 | False | -0.064 | False | True |
| HOXB13 (exploratory) | 6 | -1.95009 | -1.19562 | +0.75447 | 0.0010 | +0.034 | False | +0.084 | False | True |
| ARID1A (exploratory) | 6 | -1.94249 | -1.26190 | +0.68059 | 0.0010 | +0.059 | False | +0.077 | False | True |
| IRF4 (exploratory) | 6 | -1.92825 | -1.53392 | +0.39433 | 0.0016 | +0.180 | False | +0.379 | False | True |
| ZBTB7A (exploratory) | 6 | -1.91444 | -1.54880 | +0.36564 | 0.0041 | +0.127 | False | -0.047 | False | True |
| ZNF544 (exploratory) | 6 | -1.90356 | -1.46744 | +0.43612 | 0.0010 | +0.031 | False | +0.103 | False | True |
| FOXJ3 (exploratory) | 6 | -1.90299 | -1.62565 | +0.27734 | 0.0107 | +0.006 | False | +0.058 | False | True |
| FUBP1 (exploratory) | 6 | -1.80866 | -1.26409 | +0.54457 | 0.0010 | +0.032 | False | +0.183 | False | True |
| ZFP64 (exploratory) | 6 | -1.77533 | -1.16921 | +0.60613 | 0.0010 | +0.045 | False | -0.065 | False | True |

**V_all / S2_k20 / O1_Y2** — 12 of these 20 are `toward_young`; 44 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -2.85365 | -0.72344 | +2.13021 | 0.0054 | +0.324 | False | -0.036 | False | True |
| NCL (exploratory) | 6 | -1.69696 | -0.99601 | +0.70095 | 0.0054 | +0.020 | False | +0.404 | False | True |
| PRDM1 (exploratory) | 6 | -1.52551 | -1.03449 | +0.49102 | 0.0235 | +0.099 | False | +0.414 | False | True |
| PIAS1 (exploratory) | 6 | -1.35891 | -0.40980 | +0.94911 | 0.0054 | +0.008 | False | +0.059 | False | True |
| ZNF850 (exploratory) | 6 | -1.28385 | -0.95140 | +0.33245 | 0.0459 | +0.029 | False | +0.260 | False | True |
| XRCC5 | 6 | -1.09735 | -0.91345 | +0.18391 | 0.1567 | +0.022 | False | +0.221 | False | False |
| OSR2 (exploratory) | 6 | -1.03599 | -0.54241 | +0.49358 | 0.0054 | +0.106 | False | -0.223 | False | True |
| TBX20 (exploratory) | 6 | -1.02785 | -0.50061 | +0.52724 | 0.0054 | +0.071 | False | +0.156 | False | True |
| ZNF275 (exploratory) | 6 | -1.02304 | -0.57120 | +0.45184 | 0.0136 | +0.015 | False | +0.073 | False | True |
| MEOX1 | 6 | -1.01666 | -0.58556 | +0.43110 | 0.0110 | -0.001 | True | -0.019 | False | False |
| ZNF286B (exploratory) | 6 | -1.00927 | -0.66650 | +0.34277 | 0.0083 | +0.006 | False | +0.158 | False | True |
| ZNF322 (exploratory) | 6 | -0.98257 | -0.44397 | +0.53860 | 0.0054 | +0.011 | False | +0.096 | False | True |
| HOXD1 | 6 | -0.97273 | -0.72694 | +0.24579 | 0.0692 | +0.029 | False | -0.053 | False | False |
| ZNF420 (exploratory) | 6 | -0.96477 | -0.58952 | +0.37525 | 0.0235 | +0.015 | False | +0.030 | False | True |
| ZNF845 | 5 | -0.96044 | -0.64371 | +0.31673 | 0.0326 | -0.002 | True | -0.172 | False | False |
| MKX (exploratory) | 6 | -0.95355 | -0.67984 | +0.27371 | 0.0178 | +0.063 | False | +0.032 | False | True |
| PPARGC1A | 6 | -0.94249 | -0.88940 | +0.05309 | 0.3406 | +0.042 | False | +0.224 | False | False |
| GFI1 | 6 | -0.93186 | -0.71972 | +0.21214 | 0.0510 | +0.231 | False | +0.058 | False | False |
| RUNX1T1 | 6 | -0.91839 | -0.93767 | -0.01927 | 0.5277 | +0.043 | False | +0.260 | False | False |
| ZFP64 | 6 | -0.90023 | -0.72107 | +0.17916 | 0.1017 | +0.045 | False | -0.065 | False | False |

**V_all / S2_k20 / O2_Y1** — 17 of these 20 are `toward_young`; 219 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -3.66590 | -0.53340 | +3.13250 | 0.0008 | +0.324 | False | -0.036 | False | True |
| TOX3 | 6 | -2.16853 | -1.51797 | +0.65056 | 0.0569 | +0.014 | False | +0.155 | False | False |
| ZEB2 (exploratory) | 6 | -1.71618 | -0.41544 | +1.30074 | 0.0008 | +0.067 | False | -0.313 | False | True |
| PIAS1 (exploratory) | 6 | -1.40572 | -0.37869 | +1.02703 | 0.0008 | +0.008 | False | +0.059 | False | True |
| MEOX1 | 6 | -1.37588 | -0.41626 | +0.95962 | 0.0008 | -0.001 | True | -0.019 | False | False |
| ZNF845 | 5 | -1.35292 | -0.41146 | +0.94146 | 0.0008 | -0.002 | True | -0.172 | False | False |
| OSR2 (exploratory) | 6 | -1.33947 | -0.36246 | +0.97701 | 0.0008 | +0.106 | False | -0.223 | False | True |
| ZNF208 (exploratory) | 4 | -1.21165 | -0.35473 | +0.85692 | 0.0008 | +0.019 | False | -0.109 | False | True |
| KLF5 (exploratory) | 6 | -1.16186 | -0.43361 | +0.72825 | 0.0008 | +0.166 | False | +0.266 | False | True |
| JDP2 (exploratory) | 6 | -1.13892 | -0.20688 | +0.93204 | 0.0008 | +0.091 | False | +0.153 | False | True |
| HOXD1 (exploratory) | 6 | -1.12963 | -0.44624 | +0.68339 | 0.0008 | +0.029 | False | -0.053 | False | True |
| PRRX2 (exploratory) | 6 | -1.07871 | -0.25172 | +0.82699 | 0.0008 | +0.028 | False | +0.084 | False | True |
| ZNF470 (exploratory) | 6 | -1.06156 | -0.29662 | +0.76493 | 0.0008 | +0.006 | False | +0.054 | False | True |
| TFAP4 (exploratory) | 6 | -1.04961 | -0.26912 | +0.78049 | 0.0008 | +0.049 | False | +0.051 | False | True |
| SPDEF (exploratory) | 6 | -1.03284 | -0.30266 | +0.73018 | 0.0008 | +0.030 | False | +0.107 | False | True |
| HAND1 (exploratory) | 6 | -1.02582 | -0.33380 | +0.69202 | 0.0008 | +0.108 | False | +0.054 | False | True |
| TWIST2 (exploratory) | 6 | -1.02426 | -0.28607 | +0.73820 | 0.0008 | +0.036 | False | -0.044 | False | True |
| ZNF841 (exploratory) | 6 | -1.01587 | -0.34135 | +0.67452 | 0.0008 | +0.011 | False | -0.058 | False | True |
| DMRTB1 (exploratory) | 6 | -1.00533 | -0.30920 | +0.69613 | 0.0008 | +0.020 | False | +0.024 | False | True |
| PRDM16 (exploratory) | 6 | -1.00300 | -0.52316 | +0.47984 | 0.0008 | +0.068 | False | +0.012 | False | True |

**V_all / S2_k20 / O2_Y2** — 18 of these 20 are `toward_young`; 208 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -3.74991 | -0.50425 | +3.24566 | 0.0008 | +0.324 | False | -0.036 | False | True |
| TOX3 (exploratory) | 6 | -2.32832 | -1.43114 | +0.89718 | 0.0279 | +0.014 | False | +0.155 | False | True |
| ZEB2 (exploratory) | 6 | -1.74136 | -0.39773 | +1.34363 | 0.0008 | +0.067 | False | -0.313 | False | True |
| PIAS1 (exploratory) | 6 | -1.46475 | -0.37657 | +1.08818 | 0.0008 | +0.008 | False | +0.059 | False | True |
| MEOX1 | 6 | -1.39784 | -0.39145 | +1.00638 | 0.0008 | -0.001 | True | -0.019 | False | False |
| ZNF845 | 5 | -1.29759 | -0.37850 | +0.91909 | 0.0008 | -0.002 | True | -0.172 | False | False |
| OSR2 (exploratory) | 6 | -1.29683 | -0.33835 | +0.95848 | 0.0008 | +0.106 | False | -0.223 | False | True |
| JDP2 (exploratory) | 6 | -1.29296 | -0.21535 | +1.07762 | 0.0008 | +0.091 | False | +0.153 | False | True |
| KLF5 (exploratory) | 6 | -1.17133 | -0.36105 | +0.81027 | 0.0008 | +0.166 | False | +0.266 | False | True |
| ZNF208 (exploratory) | 4 | -1.14237 | -0.33938 | +0.80298 | 0.0008 | +0.019 | False | -0.109 | False | True |
| PRRX2 (exploratory) | 6 | -1.11114 | -0.22633 | +0.88481 | 0.0008 | +0.028 | False | +0.084 | False | True |
| TWIST2 (exploratory) | 6 | -1.07150 | -0.25796 | +0.81354 | 0.0008 | +0.036 | False | -0.044 | False | True |
| SPDEF (exploratory) | 6 | -1.06323 | -0.26342 | +0.79981 | 0.0008 | +0.030 | False | +0.107 | False | True |
| ZNF470 (exploratory) | 6 | -1.06202 | -0.25685 | +0.80517 | 0.0008 | +0.006 | False | +0.054 | False | True |
| TFAP4 (exploratory) | 6 | -1.04625 | -0.22794 | +0.81832 | 0.0008 | +0.049 | False | +0.051 | False | True |
| DMRTB1 (exploratory) | 6 | -1.03716 | -0.27508 | +0.76208 | 0.0008 | +0.020 | False | +0.024 | False | True |
| ZNF841 (exploratory) | 6 | -1.02495 | -0.31093 | +0.71401 | 0.0008 | +0.011 | False | -0.058 | False | True |
| TBX20 (exploratory) | 6 | -1.02208 | -0.24988 | +0.77219 | 0.0008 | +0.071 | False | +0.156 | False | True |
| HAND1 (exploratory) | 6 | -1.00500 | -0.30413 | +0.70087 | 0.0008 | +0.108 | False | +0.054 | False | True |
| HOXD1 (exploratory) | 6 | -1.00169 | -0.39831 | +0.60337 | 0.0008 | +0.029 | False | -0.053 | False | True |

**V_all / S2_k50 / O1_Y1** — 19 of these 20 are `toward_young`; 247 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| PRDM1 (exploratory) | 6 | -2.38834 | -1.41456 | +0.97377 | 0.0013 | +0.099 | False | +0.414 | False | True |
| NCL (exploratory) | 6 | -2.13134 | -1.29282 | +0.83852 | 0.0013 | +0.020 | False | +0.404 | False | True |
| HOXA13 (exploratory) | 6 | -2.10348 | -1.15604 | +0.94744 | 0.0013 | +0.068 | False | +0.296 | False | True |
| PLAG1 (exploratory) | 6 | -2.06176 | -1.29806 | +0.76370 | 0.0013 | +0.057 | False | -0.160 | False | True |
| DMRT1 (exploratory) | 6 | -1.83980 | -1.31365 | +0.52615 | 0.0021 | +0.128 | False | +0.328 | False | True |
| ZNF850 (exploratory) | 6 | -1.81767 | -1.31785 | +0.49983 | 0.0013 | +0.029 | False | +0.260 | False | True |
| HOXD1 (exploratory) | 6 | -1.81471 | -0.97419 | +0.84051 | 0.0013 | +0.029 | False | -0.053 | False | True |
| LBX1 (exploratory) | 6 | -1.79339 | -1.30013 | +0.49326 | 0.0013 | +0.044 | False | +0.142 | False | True |
| HOXB13 (exploratory) | 6 | -1.78387 | -1.04192 | +0.74196 | 0.0013 | +0.034 | False | +0.084 | False | True |
| ARID1A (exploratory) | 6 | -1.71644 | -1.09625 | +0.62019 | 0.0013 | +0.059 | False | +0.077 | False | True |
| CEBPB (exploratory) | 6 | -1.69760 | -1.25949 | +0.43811 | 0.0026 | +0.091 | False | -0.233 | False | True |
| CEBPD (exploratory) | 6 | -1.65247 | -1.03991 | +0.61256 | 0.0013 | +0.043 | False | -0.064 | False | True |
| SMAD1 (exploratory) | 6 | -1.62847 | -0.99386 | +0.63461 | 0.0013 | +0.098 | False | +0.366 | False | True |
| TSHZ1 (exploratory) | 6 | -1.62371 | -1.03630 | +0.58741 | 0.0013 | +0.055 | False | +0.095 | False | True |
| RUNX1T1 (exploratory) | 6 | -1.60719 | -1.33849 | +0.26870 | 0.0163 | +0.043 | False | +0.260 | False | True |
| SMAD9 | 6 | -1.59781 | -1.51634 | +0.08147 | 0.1326 | +0.011 | False | +0.252 | False | False |
| ZFP64 (exploratory) | 6 | -1.59269 | -1.02001 | +0.57268 | 0.0013 | +0.045 | False | -0.065 | False | True |
| PSMD14 (exploratory) | 6 | -1.58880 | -1.16000 | +0.42880 | 0.0021 | +0.031 | False | +0.125 | False | True |
| NFX1 (exploratory) | 6 | -1.53231 | -1.29008 | +0.24223 | 0.0330 | +0.022 | False | +0.180 | False | True |
| FUBP1 (exploratory) | 6 | -1.53019 | -1.09417 | +0.43602 | 0.0013 | +0.032 | False | +0.183 | False | True |

**V_all / S2_k50 / O1_Y2** — 16 of these 20 are `toward_young`; 53 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -2.17852 | -0.66445 | +1.51407 | 0.0051 | +0.324 | False | -0.036 | False | True |
| PRDM1 (exploratory) | 6 | -1.31278 | -0.86493 | +0.44785 | 0.0239 | +0.099 | False | +0.414 | False | True |
| NCL | 6 | -1.24995 | -0.86724 | +0.38271 | 0.0661 | +0.020 | False | +0.404 | False | False |
| PIAS1 (exploratory) | 6 | -1.17333 | -0.38658 | +0.78675 | 0.0051 | +0.008 | False | +0.059 | False | True |
| ZNF850 (exploratory) | 6 | -1.11269 | -0.80622 | +0.30647 | 0.0467 | +0.029 | False | +0.260 | False | True |
| TBX20 (exploratory) | 6 | -1.04806 | -0.46141 | +0.58665 | 0.0051 | +0.071 | False | +0.156 | False | True |
| HOXD1 (exploratory) | 6 | -0.97492 | -0.65526 | +0.31966 | 0.0268 | +0.029 | False | -0.053 | False | True |
| NFIX (exploratory) | 6 | -0.96107 | -0.61947 | +0.34160 | 0.0106 | +0.004 | False | +0.014 | False | True |
| MEOX1 | 6 | -0.95621 | -0.54365 | +0.41256 | 0.0083 | -0.001 | True | -0.019 | False | False |
| ZNF322 (exploratory) | 6 | -0.94021 | -0.41562 | +0.52460 | 0.0051 | +0.011 | False | +0.096 | False | True |
| ZNF845 | 5 | -0.88267 | -0.59033 | +0.29234 | 0.0318 | -0.002 | True | -0.172 | False | False |
| ZNF286B (exploratory) | 6 | -0.86885 | -0.58828 | +0.28057 | 0.0157 | +0.006 | False | +0.158 | False | True |
| ZNF281 (exploratory) | 6 | -0.86818 | -0.42921 | +0.43896 | 0.0051 | +0.028 | False | -0.015 | False | True |
| RBMS1 (exploratory) | 6 | -0.86536 | -0.48518 | +0.38018 | 0.0051 | +0.016 | False | +0.281 | False | True |
| MKX (exploratory) | 6 | -0.86508 | -0.59930 | +0.26578 | 0.0141 | +0.063 | False | +0.032 | False | True |
| ZNF420 (exploratory) | 6 | -0.85666 | -0.54802 | +0.30864 | 0.0367 | +0.015 | False | +0.030 | False | True |
| ZNF814 (exploratory) | 6 | -0.85037 | -0.38584 | +0.46454 | 0.0051 | +0.018 | False | +0.043 | False | True |
| ILF3 (exploratory) | 6 | -0.84569 | -0.53560 | +0.31009 | 0.0051 | +0.020 | False | -0.030 | False | True |
| XRCC5 | 6 | -0.83388 | -0.76218 | +0.07170 | 0.2948 | +0.022 | False | +0.221 | False | False |
| MBD3 (exploratory) | 6 | -0.82689 | -0.52038 | +0.30651 | 0.0141 | +0.005 | False | -0.036 | False | True |

**V_all / S2_k50 / O2_Y1** — 18 of these 20 are `toward_young`; 236 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -3.72251 | -0.59008 | +3.13243 | 0.0008 | +0.324 | False | -0.036 | False | True |
| TOX3 (exploratory) | 6 | -3.15294 | -1.65061 | +1.50233 | 0.0045 | +0.014 | False | +0.155 | False | True |
| PRDM1 (exploratory) | 6 | -1.73425 | -0.76554 | +0.96872 | 0.0008 | +0.099 | False | +0.414 | False | True |
| ZNF208 (exploratory) | 4 | -1.72839 | -0.39194 | +1.33645 | 0.0008 | +0.019 | False | -0.109 | False | True |
| PIAS1 (exploratory) | 6 | -1.56676 | -0.39640 | +1.17036 | 0.0008 | +0.008 | False | +0.059 | False | True |
| ZNF845 | 5 | -1.49186 | -0.47241 | +1.01945 | 0.0008 | -0.002 | True | -0.172 | False | False |
| MEOX1 | 6 | -1.41899 | -0.46015 | +0.95884 | 0.0008 | -0.001 | True | -0.019 | False | False |
| HOXD1 (exploratory) | 6 | -1.35789 | -0.51253 | +0.84536 | 0.0008 | +0.029 | False | -0.053 | False | True |
| HAND1 (exploratory) | 6 | -1.32447 | -0.38200 | +0.94247 | 0.0008 | +0.108 | False | +0.054 | False | True |
| KLF5 (exploratory) | 6 | -1.31821 | -0.53728 | +0.78092 | 0.0008 | +0.166 | False | +0.266 | False | True |
| XRCC5 (exploratory) | 6 | -1.28777 | -0.64024 | +0.64753 | 0.0008 | +0.022 | False | +0.221 | False | True |
| NCOA1 (exploratory) | 6 | -1.25190 | -0.37724 | +0.87467 | 0.0008 | +0.012 | False | -0.072 | False | True |
| ZNF615 (exploratory) | 6 | -1.24714 | -0.49437 | +0.75277 | 0.0008 | +0.026 | False | +0.124 | False | True |
| ZEB2 (exploratory) | 6 | -1.21148 | -0.44469 | +0.76679 | 0.0008 | +0.067 | False | -0.313 | False | True |
| PLAG1 (exploratory) | 6 | -1.18013 | -0.64850 | +0.53163 | 0.0015 | +0.057 | False | -0.160 | False | True |
| ZNF322 (exploratory) | 6 | -1.16868 | -0.32517 | +0.84351 | 0.0008 | +0.011 | False | +0.096 | False | True |
| ZNF850 (exploratory) | 6 | -1.13320 | -0.67129 | +0.46191 | 0.0015 | +0.029 | False | +0.260 | False | True |
| CDCA7L (exploratory) | 6 | -1.12732 | -0.66707 | +0.46025 | 0.0031 | +0.030 | False | -0.111 | False | True |
| SHOX (exploratory) | 4 | -1.12564 | -0.43687 | +0.68877 | 0.0008 | +0.074 | False | -0.021 | False | True |
| TCF3 (exploratory) | 6 | -1.11758 | -0.35856 | +0.75903 | 0.0008 | +0.010 | False | +0.156 | False | True |

**V_all / S2_k50 / O2_Y2** — 18 of these 20 are `toward_young`; 230 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -3.79218 | -0.56848 | +3.22371 | 0.0008 | +0.324 | False | -0.036 | False | True |
| TOX3 (exploratory) | 6 | -3.24535 | -1.58758 | +1.65777 | 0.0036 | +0.014 | False | +0.155 | False | True |
| ZNF208 (exploratory) | 4 | -1.66424 | -0.37789 | +1.28635 | 0.0008 | +0.019 | False | -0.109 | False | True |
| PIAS1 (exploratory) | 6 | -1.60232 | -0.39157 | +1.21075 | 0.0008 | +0.008 | False | +0.059 | False | True |
| PRDM1 (exploratory) | 6 | -1.59874 | -0.70028 | +0.89846 | 0.0008 | +0.099 | False | +0.414 | False | True |
| ZNF845 | 5 | -1.45213 | -0.44469 | +1.00744 | 0.0008 | -0.002 | True | -0.172 | False | False |
| MEOX1 | 6 | -1.43961 | -0.43947 | +1.00014 | 0.0008 | -0.001 | True | -0.019 | False | False |
| KLF5 (exploratory) | 6 | -1.33401 | -0.48059 | +0.85342 | 0.0008 | +0.166 | False | +0.266 | False | True |
| HAND1 (exploratory) | 6 | -1.30501 | -0.35892 | +0.94609 | 0.0008 | +0.108 | False | +0.054 | False | True |
| ZEB2 (exploratory) | 6 | -1.25240 | -0.43072 | +0.82168 | 0.0008 | +0.067 | False | -0.313 | False | True |
| HOXD1 (exploratory) | 6 | -1.24988 | -0.47709 | +0.77279 | 0.0008 | +0.029 | False | -0.053 | False | True |
| XRCC5 (exploratory) | 6 | -1.24069 | -0.57839 | +0.66230 | 0.0008 | +0.022 | False | +0.221 | False | True |
| NCOA1 (exploratory) | 6 | -1.19030 | -0.36979 | +0.82051 | 0.0008 | +0.012 | False | -0.072 | False | True |
| ZNF322 (exploratory) | 6 | -1.16571 | -0.30716 | +0.85855 | 0.0008 | +0.011 | False | +0.096 | False | True |
| ZNF615 (exploratory) | 6 | -1.16296 | -0.44399 | +0.71897 | 0.0008 | +0.026 | False | +0.124 | False | True |
| ZNF470 (exploratory) | 6 | -1.09580 | -0.32469 | +0.77111 | 0.0008 | +0.006 | False | +0.054 | False | True |
| TFAP4 (exploratory) | 6 | -1.09455 | -0.29965 | +0.79491 | 0.0008 | +0.049 | False | +0.051 | False | True |
| SHOX (exploratory) | 4 | -1.07585 | -0.40008 | +0.67578 | 0.0008 | +0.074 | False | -0.021 | False | True |
| TCF3 (exploratory) | 6 | -1.06378 | -0.33420 | +0.72958 | 0.0008 | +0.010 | False | +0.156 | False | True |
| ZNF589 (exploratory) | 6 | -1.06377 | -0.42568 | +0.63809 | 0.0008 | +0.002 | False | +0.139 | False | True |

**V_all / S2_k100 / O1_Y1** — 20 of these 20 are `toward_young`; 173 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| PRDM1 (exploratory) | 6 | -1.96796 | -1.14968 | +0.81828 | 0.0018 | +0.099 | False | +0.414 | False | True |
| HOXA13 (exploratory) | 6 | -1.87016 | -0.99942 | +0.87074 | 0.0018 | +0.068 | False | +0.296 | False | True |
| NCL (exploratory) | 6 | -1.84777 | -1.09962 | +0.74815 | 0.0018 | +0.020 | False | +0.404 | False | True |
| PLAG1 (exploratory) | 6 | -1.68989 | -1.07447 | +0.61542 | 0.0018 | +0.057 | False | -0.160 | False | True |
| HOXD1 (exploratory) | 6 | -1.58945 | -0.86840 | +0.72105 | 0.0018 | +0.029 | False | -0.053 | False | True |
| HOXB13 (exploratory) | 6 | -1.51522 | -0.90662 | +0.60861 | 0.0018 | +0.034 | False | +0.084 | False | True |
| DMRT1 (exploratory) | 6 | -1.48518 | -1.05473 | +0.43045 | 0.0032 | +0.128 | False | +0.328 | False | True |
| ARID1A (exploratory) | 6 | -1.45248 | -0.94617 | +0.50631 | 0.0018 | +0.059 | False | +0.077 | False | True |
| ZNF850 (exploratory) | 6 | -1.43945 | -1.08968 | +0.34976 | 0.0132 | +0.029 | False | +0.260 | False | True |
| SMAD1 (exploratory) | 6 | -1.42093 | -0.87814 | +0.54278 | 0.0018 | +0.098 | False | +0.366 | False | True |
| BCL11B (exploratory) | 6 | -1.38915 | -0.56593 | +0.82322 | 0.0018 | +0.030 | False | -0.056 | False | True |
| PSMD14 (exploratory) | 6 | -1.35540 | -0.98988 | +0.36553 | 0.0060 | +0.031 | False | +0.125 | False | True |
| ZNF583 (exploratory) | 6 | -1.35528 | -0.85883 | +0.49644 | 0.0018 | +0.021 | False | +0.006 | False | True |
| TSHZ1 (exploratory) | 6 | -1.34612 | -0.89738 | +0.44874 | 0.0018 | +0.055 | False | +0.095 | False | True |
| CEBPD (exploratory) | 6 | -1.33768 | -0.90863 | +0.42905 | 0.0018 | +0.043 | False | -0.064 | False | True |
| ZFP64 (exploratory) | 6 | -1.33561 | -0.89051 | +0.44510 | 0.0018 | +0.045 | False | -0.065 | False | True |
| LBX1 (exploratory) | 6 | -1.33441 | -1.04114 | +0.29327 | 0.0060 | +0.044 | False | +0.142 | False | True |
| ZNF705B (exploratory) | 4 | -1.28701 | -0.68110 | +0.60591 | 0.0018 | +0.005 | False | -0.015 | False | True |
| NFIA (exploratory) | 6 | -1.28114 | -0.86396 | +0.41718 | 0.0018 | +0.053 | False | +0.188 | False | True |
| CEBPB (exploratory) | 6 | -1.26746 | -1.03987 | +0.22759 | 0.0288 | +0.091 | False | -0.233 | False | True |

**V_all / S2_k100 / O1_Y2** — 11 of these 20 are `toward_young`; 40 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -1.55534 | -0.60648 | +0.94886 | 0.0076 | +0.324 | False | -0.036 | False | True |
| NCL | 6 | -1.02231 | -0.72320 | +0.29911 | 0.1471 | +0.020 | False | +0.404 | False | False |
| PIAS1 (exploratory) | 6 | -1.00180 | -0.36371 | +0.63809 | 0.0076 | +0.008 | False | +0.059 | False | True |
| PRDM1 | 6 | -0.95290 | -0.66286 | +0.29004 | 0.1122 | +0.099 | False | +0.414 | False | False |
| TBX20 (exploratory) | 6 | -0.90707 | -0.42519 | +0.48188 | 0.0076 | +0.071 | False | +0.156 | False | True |
| ZNF322 (exploratory) | 6 | -0.86984 | -0.38983 | +0.48000 | 0.0076 | +0.011 | False | +0.096 | False | True |
| NFIX (exploratory) | 6 | -0.83850 | -0.54649 | +0.29201 | 0.0212 | +0.004 | False | +0.014 | False | True |
| MEOX1 | 6 | -0.82433 | -0.50516 | +0.31917 | 0.0357 | -0.001 | True | -0.019 | False | False |
| ZNF850 | 6 | -0.81710 | -0.63391 | +0.18319 | 0.1801 | +0.029 | False | +0.260 | False | False |
| ZNF814 (exploratory) | 6 | -0.80966 | -0.36329 | +0.44637 | 0.0076 | +0.018 | False | +0.043 | False | True |
| ZNF281 (exploratory) | 6 | -0.79291 | -0.40359 | +0.38932 | 0.0076 | +0.028 | False | -0.015 | False | True |
| RBMS1 (exploratory) | 6 | -0.79258 | -0.43736 | +0.35522 | 0.0076 | +0.016 | False | +0.281 | False | True |
| HOXD1 | 6 | -0.79220 | -0.58120 | +0.21101 | 0.0933 | +0.029 | False | -0.053 | False | False |
| HMG20A (exploratory) | 6 | -0.78659 | -0.30813 | +0.47847 | 0.0076 | +0.013 | False | +0.021 | False | True |
| ILF3 (exploratory) | 6 | -0.77577 | -0.47925 | +0.29652 | 0.0108 | +0.020 | False | -0.030 | False | True |
| ZNF845 | 5 | -0.77576 | -0.53925 | +0.23652 | 0.0793 | -0.002 | True | -0.172 | False | False |
| ZNF420 | 6 | -0.76648 | -0.50631 | +0.26017 | 0.0688 | +0.015 | False | +0.030 | False | False |
| ZNF805 | 6 | -0.75177 | -0.37139 | +0.38038 | 0.0125 | -0.003 | True | +0.031 | False | False |
| MAFF (exploratory) | 6 | -0.74021 | -0.45173 | +0.28848 | 0.0108 | +0.026 | False | +0.095 | False | True |
| ZNF275 | 6 | -0.73788 | -0.49171 | +0.24616 | 0.0666 | +0.015 | False | +0.073 | False | False |

**V_all / S2_k100 / O2_Y1** — 18 of these 20 are `toward_young`; 248 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -3.75179 | -0.71113 | +3.04066 | 0.0008 | +0.324 | False | -0.036 | False | True |
| TOX3 (exploratory) | 6 | -3.62221 | -1.96119 | +1.66103 | 0.0032 | +0.014 | False | +0.155 | False | True |
| PRDM1 (exploratory) | 6 | -2.20173 | -1.09597 | +1.10576 | 0.0008 | +0.099 | False | +0.414 | False | True |
| ZNF208 (exploratory) | 4 | -1.87202 | -0.46272 | +1.40930 | 0.0008 | +0.019 | False | -0.109 | False | True |
| PIAS1 (exploratory) | 6 | -1.72373 | -0.41263 | +1.31111 | 0.0008 | +0.008 | False | +0.059 | False | True |
| XRCC5 (exploratory) | 6 | -1.71284 | -0.96494 | +0.74790 | 0.0008 | +0.022 | False | +0.221 | False | True |
| KLF5 (exploratory) | 6 | -1.71241 | -0.83090 | +0.88151 | 0.0008 | +0.166 | False | +0.266 | False | True |
| CDCA7L (exploratory) | 6 | -1.66618 | -0.98012 | +0.68606 | 0.0008 | +0.030 | False | -0.111 | False | True |
| MEOX1 | 6 | -1.61731 | -0.56105 | +1.05626 | 0.0008 | -0.001 | True | -0.019 | False | False |
| ZNF845 | 5 | -1.61410 | -0.61817 | +0.99593 | 0.0008 | -0.002 | True | -0.172 | False | False |
| HOXD1 (exploratory) | 6 | -1.59306 | -0.70231 | +0.89075 | 0.0008 | +0.029 | False | -0.053 | False | True |
| ZNF615 (exploratory) | 6 | -1.53075 | -0.75169 | +0.77907 | 0.0008 | +0.026 | False | +0.124 | False | True |
| PRDM16 (exploratory) | 6 | -1.50178 | -1.09151 | +0.41027 | 0.0019 | +0.068 | False | +0.012 | False | True |
| ZNF850 (exploratory) | 6 | -1.48319 | -0.98167 | +0.50152 | 0.0008 | +0.029 | False | +0.260 | False | True |
| PLAG1 (exploratory) | 6 | -1.47870 | -0.95568 | +0.52302 | 0.0008 | +0.057 | False | -0.160 | False | True |
| HAND1 (exploratory) | 6 | -1.44323 | -0.50582 | +0.93741 | 0.0008 | +0.108 | False | +0.054 | False | True |
| ZEB2 (exploratory) | 6 | -1.35886 | -0.51866 | +0.84020 | 0.0008 | +0.067 | False | -0.313 | False | True |
| SHOX (exploratory) | 4 | -1.35530 | -0.62667 | +0.72863 | 0.0008 | +0.074 | False | -0.021 | False | True |
| ZIK1 (exploratory) | 6 | -1.32059 | -0.61095 | +0.70963 | 0.0008 | +0.023 | False | +0.091 | False | True |
| SULT2A1 (exploratory) | 5 | -1.30770 | -0.77856 | +0.52914 | 0.0008 | +0.008 | False | +0.046 | False | True |

**V_all / S2_k100 / O2_Y2** — 18 of these 20 are `toward_young`; 243 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -3.82159 | -0.68701 | +3.13458 | 0.0007 | +0.324 | False | -0.036 | False | True |
| TOX3 (exploratory) | 6 | -3.71650 | -1.89058 | +1.82592 | 0.0023 | +0.014 | False | +0.155 | False | True |
| PRDM1 (exploratory) | 6 | -2.06303 | -1.03152 | +1.03151 | 0.0007 | +0.099 | False | +0.414 | False | True |
| ZNF208 (exploratory) | 4 | -1.81112 | -0.44969 | +1.36143 | 0.0007 | +0.019 | False | -0.109 | False | True |
| PIAS1 (exploratory) | 6 | -1.75483 | -0.40951 | +1.34531 | 0.0007 | +0.008 | False | +0.059 | False | True |
| KLF5 (exploratory) | 6 | -1.72088 | -0.77352 | +0.94736 | 0.0007 | +0.166 | False | +0.266 | False | True |
| XRCC5 (exploratory) | 6 | -1.66605 | -0.90016 | +0.76589 | 0.0007 | +0.022 | False | +0.221 | False | True |
| MEOX1 | 6 | -1.63400 | -0.54239 | +1.09161 | 0.0007 | -0.001 | True | -0.019 | False | False |
| ZNF845 | 5 | -1.57941 | -0.59019 | +0.98921 | 0.0007 | -0.002 | True | -0.172 | False | False |
| CDCA7L (exploratory) | 6 | -1.55871 | -0.91979 | +0.63892 | 0.0007 | +0.030 | False | -0.111 | False | True |
| HOXD1 (exploratory) | 6 | -1.48696 | -0.66639 | +0.82057 | 0.0007 | +0.029 | False | -0.053 | False | True |
| ZNF615 (exploratory) | 6 | -1.44890 | -0.70029 | +0.74861 | 0.0007 | +0.026 | False | +0.124 | False | True |
| HAND1 (exploratory) | 6 | -1.42558 | -0.48143 | +0.94415 | 0.0007 | +0.108 | False | +0.054 | False | True |
| ZNF850 (exploratory) | 6 | -1.39871 | -0.92113 | +0.47758 | 0.0018 | +0.029 | False | +0.260 | False | True |
| PRDM16 (exploratory) | 6 | -1.39641 | -1.00820 | +0.38821 | 0.0018 | +0.068 | False | +0.012 | False | True |
| ZEB2 (exploratory) | 6 | -1.39343 | -0.50383 | +0.88960 | 0.0007 | +0.067 | False | -0.313 | False | True |
| PLAG1 (exploratory) | 6 | -1.32929 | -0.89800 | +0.43129 | 0.0041 | +0.057 | False | -0.160 | False | True |
| ZIK1 (exploratory) | 6 | -1.31572 | -0.57558 | +0.74014 | 0.0007 | +0.023 | False | +0.091 | False | True |
| SULT2A1 (exploratory) | 5 | -1.31079 | -0.73437 | +0.57642 | 0.0007 | +0.008 | False | +0.046 | False | True |
| SHOX (exploratory) | 4 | -1.30688 | -0.59007 | +0.71681 | 0.0007 | +0.074 | False | -0.021 | False | True |

**V_sig / S1 / O1_Y1** — 0 of these 20 are `toward_young`; 0 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| ZNF175 | 6 | +0.00492 | +0.00731 | +0.00238 | 0.2005 | +0.014 | False | -0.000 | False | False |
| TCF25 | 6 | +0.00846 | +0.00779 | -0.00067 | 0.3114 | +0.014 | False | -0.000 | False | False |
| NR5A1 | 6 | +0.01182 | +0.00709 | -0.00472 | 0.4974 | +0.008 | False | -0.000 | False | False |
| SKI | 6 | +0.01245 | +0.00542 | -0.00703 | 0.6422 | +0.043 | False | -0.000 | False | False |
| PKNOX1 | 6 | +0.01297 | +0.00398 | -0.00899 | 0.8149 | +0.005 | False | -0.000 | False | False |
| EN1 | 6 | +0.01331 | +0.00741 | -0.00590 | 0.5373 | +0.028 | False | -0.000 | False | False |
| CREBL2 | 6 | +0.01381 | +0.00867 | -0.00514 | 0.5124 | +0.006 | False | -0.000 | False | False |
| POU5F1B | 6 | +0.01463 | +0.00641 | -0.00822 | 0.7001 | +0.006 | False | -0.000 | False | False |
| ZBTB8B | 6 | +0.01539 | +0.01637 | +0.00098 | 0.2651 | +0.002 | False | -0.000 | False | False |
| ZNF341 | 6 | +0.01557 | +0.00732 | -0.00824 | 0.6617 | +0.043 | False | -0.000 | False | False |
| MESP2 | 6 | +0.01622 | +0.00750 | -0.00873 | 0.6708 | +0.001 | False | -0.000 | False | False |
| XBP1 | 6 | +0.01970 | +0.01235 | -0.00734 | 0.5663 | +0.033 | False | -0.000 | False | False |
| NFXL1 | 6 | +0.01983 | +0.00739 | -0.01243 | 0.8764 | -0.003 | True | -0.000 | False | False |
| ZNF232 | 6 | +0.01984 | +0.00831 | -0.01153 | 0.8357 | +0.004 | False | -0.012 | False | False |
| ZNF853 | 6 | +0.01988 | +0.03372 | +0.01383 | 0.0673 | +0.022 | False | -0.000 | False | False |
| ZNF146 | 6 | +0.01990 | +0.01260 | -0.00730 | 0.5198 | +0.007 | False | -0.000 | False | False |
| CCNT2 | 6 | +0.01998 | +0.02024 | +0.00027 | 0.2826 | +0.022 | False | -0.000 | False | False |
| HMGXB4 | 6 | +0.02205 | +0.02000 | -0.00205 | 0.3385 | +0.024 | False | -0.000 | False | False |
| ISX | 6 | +0.02444 | +0.02446 | +0.00003 | 0.2884 | +0.004 | False | -0.000 | False | False |
| CTBP1 | 6 | +0.02444 | +0.01649 | -0.00795 | 0.5634 | +0.025 | False | -0.000 | False | False |

**V_sig / S1 / O2_Y1** — 18 of these 20 are `toward_young`; 128 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -1.47611 | +0.18771 | +1.66382 | 0.0020 | +0.324 | False | -0.081 | False | True |
| BSX (exploratory) | 6 | -0.48456 | -0.17129 | +0.31328 | 0.0020 | +0.048 | False | +0.021 | False | True |
| KLF5 (exploratory) | 6 | -0.41908 | -0.05090 | +0.36818 | 0.0020 | +0.166 | False | +0.126 | False | True |
| GFI1 (exploratory) | 6 | -0.38395 | -0.06931 | +0.31464 | 0.0020 | +0.231 | False | +0.061 | False | True |
| HHEX (exploratory) | 6 | -0.31065 | -0.09016 | +0.22049 | 0.0020 | +0.044 | False | +0.053 | False | True |
| NKX3-1 (exploratory) | 6 | -0.23803 | -0.03437 | +0.20366 | 0.0063 | +0.163 | False | +0.061 | False | True |
| ZFP2 (exploratory) | 6 | -0.20805 | -0.08029 | +0.12776 | 0.0020 | +0.028 | False | +0.042 | False | True |
| ELOC (exploratory) | 6 | -0.19802 | -0.04268 | +0.15534 | 0.0020 | +0.018 | False | -0.020 | False | True |
| POU3F3 | 6 | -0.19064 | -0.04895 | +0.14170 | 0.0020 | -0.004 | True | -0.040 | False | False |
| GATA1 (exploratory) | 6 | -0.19002 | +0.03142 | +0.22144 | 0.0020 | +0.047 | False | +0.019 | False | True |
| ZNF557 (exploratory) | 6 | -0.18286 | -0.07239 | +0.11048 | 0.0020 | +0.014 | False | -0.020 | False | True |
| ZNF574 | 6 | -0.17789 | -0.02735 | +0.15054 | 0.0020 | -0.018 | True | +0.027 | False | False |
| TFCP2 (exploratory) | 6 | -0.17013 | -0.04520 | +0.12493 | 0.0020 | +0.002 | False | +0.059 | False | True |
| RBMS1 (exploratory) | 6 | -0.15807 | -0.05116 | +0.10691 | 0.0020 | +0.016 | False | +0.037 | False | True |
| PRRX2 (exploratory) | 6 | -0.15292 | -0.02834 | +0.12458 | 0.0020 | +0.028 | False | +0.034 | False | True |
| MNX1 (exploratory) | 6 | -0.14675 | -0.04984 | +0.09691 | 0.0083 | +0.087 | False | +0.015 | False | True |
| TBX20 (exploratory) | 6 | -0.14334 | -0.03776 | +0.10559 | 0.0052 | +0.071 | False | +0.077 | False | True |
| MYOG (exploratory) | 6 | -0.14008 | +0.01716 | +0.15724 | 0.0020 | +0.008 | False | -0.049 | False | True |
| MXI1 (exploratory) | 6 | -0.13846 | +0.04990 | +0.18836 | 0.0020 | +0.075 | False | -0.000 | False | True |
| CLOCK (exploratory) | 6 | -0.13785 | -0.04141 | +0.09644 | 0.0020 | +0.034 | False | -0.004 | False | True |

**V_sig / S2_k20 / O1_Y1** — 15 of these 20 are `toward_young`; 67 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 | 6 | -0.66365 | -0.47855 | +0.18509 | 0.0758 | +0.324 | False | -0.081 | False | False |
| BSX (exploratory) | 6 | -0.53870 | -0.27988 | +0.25882 | 0.0044 | +0.048 | False | +0.021 | False | True |
| GFI1 (exploratory) | 6 | -0.50501 | -0.18291 | +0.32210 | 0.0044 | +0.231 | False | +0.061 | False | True |
| HHEX (exploratory) | 6 | -0.29910 | -0.18231 | +0.11679 | 0.0146 | +0.044 | False | +0.053 | False | True |
| IRF1 (exploratory) | 6 | -0.22983 | +0.11291 | +0.34274 | 0.0095 | +0.200 | False | +0.025 | False | True |
| OLIG1 (exploratory) | 6 | -0.19257 | -0.06272 | +0.12984 | 0.0044 | +0.018 | False | +0.026 | False | True |
| PRDM16 (exploratory) | 6 | -0.18864 | -0.03443 | +0.15420 | 0.0044 | +0.068 | False | -0.111 | False | True |
| KLF5 | 6 | -0.18758 | -0.19985 | -0.01227 | 0.2868 | +0.166 | False | +0.126 | False | False |
| PRDM1 (exploratory) | 6 | -0.17827 | +0.11120 | +0.28947 | 0.0044 | +0.099 | False | +0.073 | False | True |
| ZNF763 (exploratory) | 6 | -0.16801 | -0.01416 | +0.15385 | 0.0044 | +0.008 | False | +0.053 | False | True |
| ETS1 (exploratory) | 6 | -0.16546 | -0.04443 | +0.12103 | 0.0044 | +0.044 | False | +0.041 | False | True |
| SOX11 (exploratory) | 5 | -0.16103 | -0.01336 | +0.14767 | 0.0076 | +0.045 | False | +0.010 | False | True |
| RBMS1 (exploratory) | 6 | -0.15626 | -0.06510 | +0.09116 | 0.0044 | +0.016 | False | +0.037 | False | True |
| LHX4 (exploratory) | 6 | -0.15404 | -0.01977 | +0.13427 | 0.0044 | +0.084 | False | -0.039 | False | True |
| MYOD1 (exploratory) | 6 | -0.15278 | -0.06492 | +0.08785 | 0.0044 | +0.065 | False | +0.013 | False | True |
| CDX2 (exploratory) | 6 | -0.15203 | -0.02499 | +0.12704 | 0.0128 | +0.093 | False | +0.041 | False | True |
| HNRNPAB (exploratory) | 6 | -0.13949 | +0.18672 | +0.32622 | 0.0044 | +0.029 | False | -0.000 | False | True |
| POU3F3 | 6 | -0.13350 | -0.06159 | +0.07191 | 0.0044 | -0.004 | True | -0.040 | False | False |
| NKX3-1 | 6 | -0.13251 | -0.18325 | -0.05073 | 0.4169 | +0.163 | False | +0.061 | False | False |
| SOX5 | 6 | -0.13147 | +0.01148 | +0.14295 | 0.0044 | -0.002 | True | +0.012 | False | False |

**V_sig / S2_k20 / O1_Y2** — 18 of these 20 are `toward_young`; 45 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -1.76899 | -0.40656 | +1.36243 | 0.0044 | +0.324 | False | -0.081 | False | True |
| BSX (exploratory) | 6 | -0.59579 | -0.21137 | +0.38442 | 0.0044 | +0.048 | False | +0.021 | False | True |
| NKX3-1 (exploratory) | 6 | -0.50802 | -0.17322 | +0.33480 | 0.0044 | +0.163 | False | +0.061 | False | True |
| KLF5 (exploratory) | 6 | -0.50050 | -0.18143 | +0.31908 | 0.0044 | +0.166 | False | +0.126 | False | True |
| GFI1 (exploratory) | 6 | -0.46155 | -0.15925 | +0.30231 | 0.0044 | +0.231 | False | +0.061 | False | True |
| HHEX (exploratory) | 6 | -0.38523 | -0.15546 | +0.22977 | 0.0044 | +0.044 | False | +0.053 | False | True |
| HNRNPAB (exploratory) | 6 | -0.36118 | -0.00568 | +0.35550 | 0.0044 | +0.029 | False | -0.000 | False | True |
| MXI1 (exploratory) | 6 | -0.27466 | -0.04909 | +0.22557 | 0.0044 | +0.075 | False | -0.000 | False | True |
| POU3F3 | 6 | -0.24870 | -0.05628 | +0.19242 | 0.0044 | -0.004 | True | -0.040 | False | False |
| MNX1 (exploratory) | 6 | -0.23618 | -0.09366 | +0.14252 | 0.0044 | +0.087 | False | +0.015 | False | True |
| PRDM1 (exploratory) | 6 | -0.23410 | -0.00441 | +0.22969 | 0.0089 | +0.099 | False | +0.073 | False | True |
| PIAS1 (exploratory) | 6 | -0.21816 | -0.00134 | +0.21681 | 0.0063 | +0.008 | False | +0.015 | False | True |
| IRF1 | 6 | -0.21396 | -0.06170 | +0.15226 | 0.1250 | +0.200 | False | +0.025 | False | False |
| ETS1 (exploratory) | 6 | -0.19745 | -0.05827 | +0.13918 | 0.0044 | +0.044 | False | +0.041 | False | True |
| ZFP2 (exploratory) | 6 | -0.18757 | -0.08562 | +0.10195 | 0.0044 | +0.028 | False | +0.042 | False | True |
| TBX20 (exploratory) | 6 | -0.17291 | -0.07000 | +0.10292 | 0.0044 | +0.071 | False | +0.077 | False | True |
| ZNF763 (exploratory) | 6 | -0.15757 | -0.04102 | +0.11655 | 0.0044 | +0.008 | False | +0.053 | False | True |
| FOXD2 (exploratory) | 6 | -0.15525 | -0.06357 | +0.09168 | 0.0089 | +0.037 | False | +0.146 | False | True |
| SP140 (exploratory) | 6 | -0.15082 | +0.02870 | +0.17952 | 0.0299 | +0.006 | False | -0.061 | False | True |
| RHOXF2B (exploratory) | 4 | -0.14824 | -0.03095 | +0.11729 | 0.0063 | +0.003 | False | -0.000 | False | True |

**V_sig / S2_k20 / O2_Y1** — 13 of these 20 are `toward_young`; 40 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -1.56271 | -0.32989 | +1.23281 | 0.0057 | +0.324 | False | -0.081 | False | True |
| KLF5 (exploratory) | 6 | -0.54361 | -0.15648 | +0.38712 | 0.0057 | +0.166 | False | +0.126 | False | True |
| ZEB2 (exploratory) | 6 | -0.35116 | -0.09190 | +0.25926 | 0.0057 | +0.067 | False | -0.208 | False | True |
| GFI1 (exploratory) | 6 | -0.34893 | -0.13585 | +0.21308 | 0.0057 | +0.231 | False | +0.061 | False | True |
| BSX (exploratory) | 6 | -0.34840 | -0.14343 | +0.20497 | 0.0083 | +0.048 | False | +0.021 | False | True |
| CREB5 | 6 | -0.33051 | -0.17088 | +0.15962 | 0.1015 | +0.013 | False | -0.012 | False | False |
| NCOA1 (exploratory) | 6 | -0.31245 | -0.11932 | +0.19313 | 0.0115 | +0.012 | False | +0.007 | False | True |
| NKX3-1 | 6 | -0.25801 | -0.15349 | +0.10452 | 0.0707 | +0.163 | False | +0.061 | False | False |
| MXI1 (exploratory) | 6 | -0.25336 | -0.08771 | +0.16566 | 0.0083 | +0.075 | False | -0.000 | False | True |
| JDP2 (exploratory) | 6 | -0.24437 | -0.06985 | +0.17451 | 0.0083 | +0.091 | False | +0.098 | False | True |
| AFF4 | 6 | -0.23121 | -0.12912 | +0.10209 | 0.1015 | +0.023 | False | -0.000 | False | False |
| MYOG (exploratory) | 6 | -0.22093 | -0.05301 | +0.16793 | 0.0057 | +0.008 | False | -0.049 | False | True |
| GATA1 (exploratory) | 6 | -0.20039 | -0.05675 | +0.14364 | 0.0057 | +0.047 | False | +0.019 | False | True |
| PIAS1 | 6 | -0.19780 | -0.11246 | +0.08534 | 0.1142 | +0.008 | False | +0.015 | False | False |
| HHEX | 6 | -0.19725 | -0.12481 | +0.07245 | 0.1015 | +0.044 | False | +0.053 | False | False |
| PAX2 (exploratory) | 6 | -0.18867 | -0.09064 | +0.09803 | 0.0387 | +0.114 | False | -0.149 | False | True |
| POU3F3 | 6 | -0.18418 | -0.04892 | +0.13526 | 0.0057 | -0.004 | True | -0.040 | False | False |
| MNX1 (exploratory) | 6 | -0.17668 | -0.08583 | +0.09085 | 0.0326 | +0.087 | False | +0.015 | False | True |
| ZFP2 (exploratory) | 6 | -0.17415 | -0.06546 | +0.10869 | 0.0057 | +0.028 | False | +0.042 | False | True |
| ZNF420 | 6 | -0.16158 | -0.09524 | +0.06634 | 0.1061 | +0.015 | False | +0.008 | False | False |

**V_sig / S2_k20 / O2_Y2** — 14 of these 20 are `toward_young`; 45 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -1.69098 | -0.32049 | +1.37049 | 0.0048 | +0.324 | False | -0.081 | False | True |
| KLF5 (exploratory) | 6 | -0.57713 | -0.15372 | +0.42341 | 0.0048 | +0.166 | False | +0.126 | False | True |
| ZEB2 (exploratory) | 6 | -0.38152 | -0.10657 | +0.27495 | 0.0048 | +0.067 | False | -0.208 | False | True |
| BSX (exploratory) | 6 | -0.35605 | -0.13441 | +0.22164 | 0.0048 | +0.048 | False | +0.021 | False | True |
| CREB5 | 6 | -0.34534 | -0.20306 | +0.14228 | 0.1217 | +0.013 | False | -0.012 | False | False |
| GFI1 (exploratory) | 6 | -0.33971 | -0.13251 | +0.20719 | 0.0048 | +0.231 | False | +0.061 | False | True |
| NCOA1 (exploratory) | 6 | -0.31160 | -0.12790 | +0.18371 | 0.0158 | +0.012 | False | +0.007 | False | True |
| NKX3-1 (exploratory) | 6 | -0.31147 | -0.15234 | +0.15913 | 0.0106 | +0.163 | False | +0.061 | False | True |
| JDP2 (exploratory) | 6 | -0.27962 | -0.07821 | +0.20141 | 0.0076 | +0.091 | False | +0.098 | False | True |
| MXI1 (exploratory) | 6 | -0.27237 | -0.09283 | +0.17954 | 0.0076 | +0.075 | False | -0.000 | False | True |
| AFF4 | 6 | -0.25003 | -0.14045 | +0.10958 | 0.0837 | +0.023 | False | -0.000 | False | False |
| MYOG (exploratory) | 6 | -0.24490 | -0.05637 | +0.18853 | 0.0048 | +0.008 | False | -0.049 | False | True |
| GATA1 (exploratory) | 6 | -0.23109 | -0.06113 | +0.16996 | 0.0048 | +0.047 | False | +0.019 | False | True |
| PIAS1 | 6 | -0.22840 | -0.12612 | +0.10228 | 0.0724 | +0.008 | False | +0.015 | False | False |
| PAX2 (exploratory) | 6 | -0.21040 | -0.09662 | +0.11377 | 0.0204 | +0.114 | False | -0.149 | False | True |
| HHEX | 6 | -0.21015 | -0.12111 | +0.08904 | 0.0554 | +0.044 | False | +0.053 | False | False |
| MNX1 (exploratory) | 6 | -0.20178 | -0.08476 | +0.11702 | 0.0048 | +0.087 | False | +0.015 | False | True |
| POU3F3 | 6 | -0.19865 | -0.04807 | +0.15058 | 0.0048 | -0.004 | True | -0.040 | False | False |
| HNRNPAB | 6 | -0.19291 | -0.20632 | -0.01341 | 0.3997 | +0.029 | False | -0.000 | False | False |
| ZFP2 (exploratory) | 6 | -0.18434 | -0.06286 | +0.12148 | 0.0048 | +0.028 | False | +0.042 | False | True |

**V_sig / S2_k50 / O1_Y1** — 14 of these 20 are `toward_young`; 62 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 | 6 | -0.50325 | -0.44820 | +0.05505 | 0.2048 | +0.324 | False | -0.081 | False | False |
| BSX (exploratory) | 6 | -0.45314 | -0.26191 | +0.19122 | 0.0046 | +0.048 | False | +0.021 | False | True |
| GFI1 (exploratory) | 6 | -0.43713 | -0.17427 | +0.26286 | 0.0046 | +0.231 | False | +0.061 | False | True |
| HHEX (exploratory) | 6 | -0.28901 | -0.17260 | +0.11641 | 0.0139 | +0.044 | False | +0.053 | False | True |
| OLIG1 (exploratory) | 6 | -0.20966 | -0.06091 | +0.14876 | 0.0046 | +0.018 | False | +0.026 | False | True |
| KLF5 | 6 | -0.20653 | -0.19144 | +0.01509 | 0.2326 | +0.166 | False | +0.126 | False | False |
| PRDM1 (exploratory) | 6 | -0.20250 | +0.09639 | +0.29888 | 0.0046 | +0.099 | False | +0.073 | False | True |
| PRDM16 (exploratory) | 6 | -0.18431 | -0.03622 | +0.14809 | 0.0046 | +0.068 | False | -0.111 | False | True |
| NKX3-1 | 6 | -0.17257 | -0.17670 | -0.00412 | 0.2858 | +0.163 | False | +0.061 | False | False |
| CDX2 (exploratory) | 6 | -0.16898 | -0.02940 | +0.13958 | 0.0046 | +0.093 | False | +0.041 | False | True |
| SOX11 (exploratory) | 5 | -0.16510 | -0.01787 | +0.14723 | 0.0076 | +0.045 | False | +0.010 | False | True |
| RBMS1 (exploratory) | 6 | -0.16196 | -0.06228 | +0.09969 | 0.0046 | +0.016 | False | +0.037 | False | True |
| ZNF714 | 6 | -0.16075 | -0.06089 | +0.09986 | 0.0856 | +0.016 | False | -0.000 | False | False |
| ZNF763 (exploratory) | 6 | -0.16054 | -0.01770 | +0.14284 | 0.0046 | +0.008 | False | +0.053 | False | True |
| MYOD1 (exploratory) | 6 | -0.15921 | -0.06386 | +0.09535 | 0.0046 | +0.065 | False | +0.013 | False | True |
| ETS1 (exploratory) | 6 | -0.15916 | -0.04585 | +0.11331 | 0.0046 | +0.044 | False | +0.041 | False | True |
| LHX4 (exploratory) | 6 | -0.15258 | -0.02345 | +0.12913 | 0.0046 | +0.084 | False | -0.039 | False | True |
| POU3F3 | 6 | -0.13171 | -0.05958 | +0.07213 | 0.0046 | -0.004 | True | -0.040 | False | False |
| MXI1 (exploratory) | 6 | -0.12489 | -0.00920 | +0.11568 | 0.0095 | +0.075 | False | -0.000 | False | True |
| SOX5 | 6 | -0.11923 | +0.00617 | +0.12540 | 0.0046 | -0.002 | True | +0.012 | False | False |

**V_sig / S2_k50 / O1_Y2** — 17 of these 20 are `toward_young`; 43 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -1.47333 | -0.38214 | +1.09119 | 0.0046 | +0.324 | False | -0.081 | False | True |
| BSX (exploratory) | 6 | -0.50456 | -0.20086 | +0.30370 | 0.0046 | +0.048 | False | +0.021 | False | True |
| NKX3-1 (exploratory) | 6 | -0.50036 | -0.16608 | +0.33428 | 0.0046 | +0.163 | False | +0.061 | False | True |
| KLF5 (exploratory) | 6 | -0.49372 | -0.17447 | +0.31925 | 0.0046 | +0.166 | False | +0.126 | False | True |
| GFI1 (exploratory) | 6 | -0.39506 | -0.15285 | +0.24221 | 0.0046 | +0.231 | False | +0.061 | False | True |
| HHEX (exploratory) | 6 | -0.35536 | -0.14912 | +0.20624 | 0.0046 | +0.044 | False | +0.053 | False | True |
| HNRNPAB (exploratory) | 6 | -0.31781 | -0.00015 | +0.31766 | 0.0046 | +0.029 | False | -0.000 | False | True |
| PRDM1 (exploratory) | 6 | -0.26395 | -0.00519 | +0.25876 | 0.0071 | +0.099 | False | +0.073 | False | True |
| MXI1 (exploratory) | 6 | -0.25680 | -0.04867 | +0.20813 | 0.0046 | +0.075 | False | -0.000 | False | True |
| POU3F3 | 6 | -0.23279 | -0.05456 | +0.17824 | 0.0046 | -0.004 | True | -0.040 | False | False |
| ETS1 (exploratory) | 6 | -0.17890 | -0.05709 | +0.12181 | 0.0046 | +0.044 | False | +0.041 | False | True |
| ZFP2 (exploratory) | 6 | -0.17601 | -0.08290 | +0.09311 | 0.0046 | +0.028 | False | +0.042 | False | True |
| MNX1 | 6 | -0.16814 | -0.09068 | +0.07747 | 0.0514 | +0.087 | False | +0.015 | False | False |
| TBX20 (exploratory) | 6 | -0.16453 | -0.06766 | +0.09687 | 0.0046 | +0.071 | False | +0.077 | False | True |
| PIAS1 (exploratory) | 6 | -0.15752 | -0.00309 | +0.15443 | 0.0148 | +0.008 | False | +0.015 | False | True |
| ELOC (exploratory) | 6 | -0.15741 | -0.05058 | +0.10684 | 0.0046 | +0.018 | False | -0.020 | False | True |
| ZNF763 (exploratory) | 6 | -0.14628 | -0.04154 | +0.10473 | 0.0046 | +0.008 | False | +0.053 | False | True |
| OLIG1 (exploratory) | 6 | -0.14117 | -0.05788 | +0.08328 | 0.0046 | +0.018 | False | +0.026 | False | True |
| AFF4 | 6 | -0.13512 | -0.04025 | +0.09487 | 0.1006 | +0.023 | False | -0.000 | False | False |
| PRDM16 (exploratory) | 6 | -0.12826 | -0.05436 | +0.07390 | 0.0270 | +0.068 | False | -0.111 | False | True |

**V_sig / S2_k50 / O2_Y1** — 18 of these 20 are `toward_young`; 79 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -1.59107 | -0.35121 | +1.23986 | 0.0042 | +0.324 | False | -0.081 | False | True |
| KLF5 (exploratory) | 6 | -0.52448 | -0.16446 | +0.36002 | 0.0042 | +0.166 | False | +0.126 | False | True |
| CREB5 (exploratory) | 6 | -0.41561 | -0.13728 | +0.27833 | 0.0252 | +0.013 | False | -0.012 | False | True |
| GFI1 (exploratory) | 6 | -0.39884 | -0.14394 | +0.25490 | 0.0042 | +0.231 | False | +0.061 | False | True |
| BSX (exploratory) | 6 | -0.36578 | -0.15827 | +0.20750 | 0.0042 | +0.048 | False | +0.021 | False | True |
| NCOA1 (exploratory) | 6 | -0.33414 | -0.11070 | +0.22343 | 0.0063 | +0.012 | False | +0.007 | False | True |
| SP140 (exploratory) | 6 | -0.29254 | -0.11032 | +0.18222 | 0.0262 | +0.006 | False | -0.061 | False | True |
| NKX3-1 (exploratory) | 6 | -0.28238 | -0.16229 | +0.12009 | 0.0411 | +0.163 | False | +0.061 | False | True |
| ZEB2 (exploratory) | 6 | -0.27964 | -0.07655 | +0.20309 | 0.0042 | +0.067 | False | -0.208 | False | True |
| PRDM1 (exploratory) | 6 | -0.27917 | -0.10136 | +0.17782 | 0.0124 | +0.099 | False | +0.073 | False | True |
| HHEX (exploratory) | 6 | -0.24886 | -0.13606 | +0.11280 | 0.0262 | +0.044 | False | +0.053 | False | True |
| BCL6B | 6 | -0.24752 | -0.10333 | +0.14419 | 0.0145 | -0.001 | True | -0.000 | False | False |
| MYOG (exploratory) | 6 | -0.24270 | -0.05040 | +0.19230 | 0.0042 | +0.008 | False | -0.049 | False | True |
| GATA1 (exploratory) | 6 | -0.23835 | -0.05346 | +0.18490 | 0.0042 | +0.047 | False | +0.019 | False | True |
| AFF4 | 6 | -0.22533 | -0.12022 | +0.10511 | 0.0673 | +0.023 | False | -0.000 | False | False |
| MXI1 (exploratory) | 6 | -0.21831 | -0.08452 | +0.13379 | 0.0063 | +0.075 | False | -0.000 | False | True |
| PAX2 (exploratory) | 6 | -0.21684 | -0.08575 | +0.13108 | 0.0083 | +0.114 | False | -0.149 | False | True |
| ZFP2 (exploratory) | 6 | -0.21130 | -0.07055 | +0.14074 | 0.0042 | +0.028 | False | +0.042 | False | True |
| PROX1 (exploratory) | 6 | -0.19275 | -0.06168 | +0.13107 | 0.0042 | +0.120 | False | -0.141 | False | True |
| MNX1 (exploratory) | 6 | -0.19087 | -0.08998 | +0.10089 | 0.0175 | +0.087 | False | +0.015 | False | True |

**V_sig / S2_k50 / O2_Y2** — 18 of these 20 are `toward_young`; 80 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -1.68975 | -0.34299 | +1.34676 | 0.0037 | +0.324 | False | -0.081 | False | True |
| KLF5 (exploratory) | 6 | -0.55402 | -0.16182 | +0.39220 | 0.0037 | +0.166 | False | +0.126 | False | True |
| CREB5 (exploratory) | 6 | -0.41949 | -0.16147 | +0.25802 | 0.0258 | +0.013 | False | -0.012 | False | True |
| GFI1 (exploratory) | 6 | -0.39160 | -0.14164 | +0.24996 | 0.0037 | +0.231 | False | +0.061 | False | True |
| BSX (exploratory) | 6 | -0.37179 | -0.15084 | +0.22095 | 0.0037 | +0.048 | False | +0.021 | False | True |
| NCOA1 (exploratory) | 6 | -0.33330 | -0.11670 | +0.21660 | 0.0059 | +0.012 | False | +0.007 | False | True |
| NKX3-1 (exploratory) | 6 | -0.32204 | -0.16027 | +0.16177 | 0.0141 | +0.163 | False | +0.061 | False | True |
| ZEB2 (exploratory) | 6 | -0.30597 | -0.08786 | +0.21811 | 0.0037 | +0.067 | False | -0.208 | False | True |
| SP140 (exploratory) | 6 | -0.30251 | -0.12447 | +0.17803 | 0.0277 | +0.006 | False | -0.061 | False | True |
| PRDM1 (exploratory) | 6 | -0.28440 | -0.11092 | +0.17347 | 0.0131 | +0.099 | False | +0.073 | False | True |
| GATA1 (exploratory) | 6 | -0.26186 | -0.05675 | +0.20511 | 0.0037 | +0.047 | False | +0.019 | False | True |
| MYOG (exploratory) | 6 | -0.26154 | -0.05260 | +0.20894 | 0.0037 | +0.008 | False | -0.049 | False | True |
| HHEX (exploratory) | 6 | -0.25678 | -0.13266 | +0.12412 | 0.0217 | +0.044 | False | +0.053 | False | True |
| BCL6B | 6 | -0.24917 | -0.10559 | +0.14358 | 0.0112 | -0.001 | True | -0.000 | False | False |
| AFF4 | 6 | -0.23833 | -0.12824 | +0.11009 | 0.0594 | +0.023 | False | -0.000 | False | False |
| PAX2 (exploratory) | 6 | -0.23295 | -0.09041 | +0.14254 | 0.0059 | +0.114 | False | -0.149 | False | True |
| MXI1 (exploratory) | 6 | -0.23293 | -0.08829 | +0.14464 | 0.0059 | +0.075 | False | -0.000 | False | True |
| ZFP2 (exploratory) | 6 | -0.21860 | -0.06858 | +0.15002 | 0.0037 | +0.028 | False | +0.042 | False | True |
| MNX1 (exploratory) | 6 | -0.21039 | -0.08931 | +0.12108 | 0.0099 | +0.087 | False | +0.015 | False | True |
| PROX1 (exploratory) | 6 | -0.19136 | -0.06360 | +0.12776 | 0.0037 | +0.120 | False | -0.141 | False | True |

**V_sig / S2_k100 / O1_Y1** — 15 of these 20 are `toward_young`; 48 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| BSX (exploratory) | 6 | -0.43806 | -0.25158 | +0.18648 | 0.0048 | +0.048 | False | +0.021 | False | True |
| GFI1 (exploratory) | 6 | -0.37277 | -0.16798 | +0.20479 | 0.0048 | +0.231 | False | +0.061 | False | True |
| HHEX | 6 | -0.24377 | -0.16621 | +0.07756 | 0.0607 | +0.044 | False | +0.053 | False | False |
| OLIG1 (exploratory) | 6 | -0.19165 | -0.05922 | +0.13243 | 0.0048 | +0.018 | False | +0.026 | False | True |
| ZNF714 (exploratory) | 6 | -0.19024 | -0.05233 | +0.13791 | 0.0408 | +0.016 | False | -0.000 | False | True |
| PRDM16 (exploratory) | 6 | -0.18220 | -0.03511 | +0.14710 | 0.0048 | +0.068 | False | -0.111 | False | True |
| CDX2 (exploratory) | 6 | -0.17600 | -0.02632 | +0.14968 | 0.0048 | +0.093 | False | +0.041 | False | True |
| PRDM1 (exploratory) | 6 | -0.17451 | +0.09941 | +0.27392 | 0.0048 | +0.099 | False | +0.073 | False | True |
| KLF5 | 6 | -0.17384 | -0.18380 | -0.00996 | 0.3324 | +0.166 | False | +0.126 | False | False |
| SOX11 (exploratory) | 5 | -0.17159 | -0.01656 | +0.15504 | 0.0048 | +0.045 | False | +0.010 | False | True |
| MYOD1 (exploratory) | 6 | -0.16035 | -0.06242 | +0.09793 | 0.0048 | +0.065 | False | +0.013 | False | True |
| RBMS1 (exploratory) | 6 | -0.15518 | -0.06093 | +0.09425 | 0.0048 | +0.016 | False | +0.037 | False | True |
| ZNF763 (exploratory) | 6 | -0.14534 | -0.01761 | +0.12773 | 0.0048 | +0.008 | False | +0.053 | False | True |
| ETS1 (exploratory) | 6 | -0.14136 | -0.04448 | +0.09688 | 0.0048 | +0.044 | False | +0.041 | False | True |
| NKX3-1 | 6 | -0.13391 | -0.17029 | -0.03638 | 0.4161 | +0.163 | False | +0.061 | False | False |
| HELZ2 | 6 | -0.13202 | -0.42116 | -0.28914 | 0.7989 | +0.324 | False | -0.081 | False | False |
| LHX4 (exploratory) | 6 | -0.12712 | -0.02160 | +0.10552 | 0.0076 | +0.084 | False | -0.039 | False | True |
| POU3F3 | 6 | -0.12679 | -0.05845 | +0.06834 | 0.0048 | -0.004 | True | -0.040 | False | False |
| MXI1 (exploratory) | 6 | -0.11846 | -0.00851 | +0.10995 | 0.0098 | +0.075 | False | -0.000 | False | True |
| ARID1A (exploratory) | 6 | -0.10880 | +0.00507 | +0.11387 | 0.0048 | +0.059 | False | +0.057 | False | True |

**V_sig / S2_k100 / O1_Y2** — 16 of these 20 are `toward_young`; 40 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -1.10573 | -0.36169 | +0.74404 | 0.0061 | +0.324 | False | -0.081 | False | True |
| BSX (exploratory) | 6 | -0.48223 | -0.19428 | +0.28794 | 0.0061 | +0.048 | False | +0.021 | False | True |
| NKX3-1 (exploratory) | 6 | -0.45610 | -0.16043 | +0.29568 | 0.0061 | +0.163 | False | +0.061 | False | True |
| KLF5 (exploratory) | 6 | -0.45305 | -0.16750 | +0.28556 | 0.0061 | +0.166 | False | +0.126 | False | True |
| GFI1 (exploratory) | 6 | -0.32678 | -0.14878 | +0.17800 | 0.0061 | +0.231 | False | +0.061 | False | True |
| HHEX (exploratory) | 6 | -0.30507 | -0.14449 | +0.16058 | 0.0061 | +0.044 | False | +0.053 | False | True |
| HNRNPAB (exploratory) | 6 | -0.29884 | +0.01778 | +0.31663 | 0.0061 | +0.029 | False | -0.000 | False | True |
| MXI1 (exploratory) | 6 | -0.24589 | -0.04637 | +0.19952 | 0.0061 | +0.075 | False | -0.000 | False | True |
| PRDM1 (exploratory) | 6 | -0.22807 | +0.00132 | +0.22938 | 0.0092 | +0.099 | False | +0.073 | False | True |
| POU3F3 | 6 | -0.22803 | -0.05348 | +0.17455 | 0.0061 | -0.004 | True | -0.040 | False | False |
| MNX1 | 6 | -0.16785 | -0.08797 | +0.07988 | 0.0559 | +0.087 | False | +0.015 | False | False |
| ZFP2 (exploratory) | 6 | -0.16380 | -0.08042 | +0.08338 | 0.0061 | +0.028 | False | +0.042 | False | True |
| ETS1 (exploratory) | 6 | -0.16028 | -0.05554 | +0.10473 | 0.0061 | +0.044 | False | +0.041 | False | True |
| ELOC (exploratory) | 6 | -0.15866 | -0.04963 | +0.10903 | 0.0061 | +0.018 | False | -0.020 | False | True |
| AFF4 | 6 | -0.13956 | -0.03541 | +0.10415 | 0.0820 | +0.023 | False | -0.000 | False | False |
| ZNF714 | 6 | -0.13492 | -0.11275 | +0.02216 | 0.3360 | +0.016 | False | -0.000 | False | False |
| ZNF763 (exploratory) | 6 | -0.13422 | -0.04007 | +0.09415 | 0.0092 | +0.008 | False | +0.053 | False | True |
| PIAS1 (exploratory) | 6 | -0.13088 | +0.00385 | +0.13473 | 0.0372 | +0.008 | False | +0.015 | False | True |
| PRDM16 (exploratory) | 6 | -0.12990 | -0.05286 | +0.07704 | 0.0275 | +0.068 | False | -0.111 | False | True |
| TBX20 (exploratory) | 6 | -0.12974 | -0.06653 | +0.06321 | 0.0306 | +0.071 | False | +0.077 | False | True |

**V_sig / S2_k100 / O2_Y1** — 16 of these 20 are `toward_young`; 80 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -1.58452 | -0.39715 | +1.18737 | 0.0034 | +0.324 | False | -0.081 | False | True |
| KLF5 (exploratory) | 6 | -0.59052 | -0.17761 | +0.41291 | 0.0034 | +0.166 | False | +0.126 | False | True |
| GFI1 (exploratory) | 6 | -0.46090 | -0.15832 | +0.30258 | 0.0034 | +0.231 | False | +0.061 | False | True |
| BSX (exploratory) | 6 | -0.42553 | -0.19625 | +0.22928 | 0.0034 | +0.048 | False | +0.021 | False | True |
| HHEX (exploratory) | 6 | -0.30767 | -0.15218 | +0.15548 | 0.0090 | +0.044 | False | +0.053 | False | True |
| NCOA1 (exploratory) | 6 | -0.28039 | -0.07541 | +0.20498 | 0.0074 | +0.012 | False | +0.007 | False | True |
| CREB5 (exploratory) | 6 | -0.26621 | -0.00306 | +0.26315 | 0.0260 | +0.013 | False | -0.012 | False | True |
| BCL6B | 6 | -0.25634 | -0.08814 | +0.16820 | 0.0090 | -0.001 | True | -0.000 | False | False |
| NKX3-1 | 6 | -0.25453 | -0.17247 | +0.08206 | 0.0792 | +0.163 | False | +0.061 | False | False |
| ZEB2 (exploratory) | 6 | -0.24408 | -0.01969 | +0.22439 | 0.0034 | +0.067 | False | -0.208 | False | True |
| GATA1 (exploratory) | 6 | -0.24381 | -0.03574 | +0.20807 | 0.0034 | +0.047 | False | +0.019 | False | True |
| PRDM1 (exploratory) | 6 | -0.24265 | -0.04474 | +0.19791 | 0.0074 | +0.099 | False | +0.073 | False | True |
| MYOG (exploratory) | 6 | -0.24037 | -0.03610 | +0.20427 | 0.0034 | +0.008 | False | -0.049 | False | True |
| ZFP2 (exploratory) | 6 | -0.23952 | -0.08192 | +0.15759 | 0.0034 | +0.028 | False | +0.042 | False | True |
| PAX2 (exploratory) | 6 | -0.23450 | -0.06154 | +0.17296 | 0.0034 | +0.114 | False | -0.149 | False | True |
| SP140 (exploratory) | 6 | -0.21428 | -0.02887 | +0.18541 | 0.0233 | +0.006 | False | -0.061 | False | True |
| MNX1 (exploratory) | 6 | -0.19978 | -0.09322 | +0.10656 | 0.0128 | +0.087 | False | +0.015 | False | True |
| FOXD2 (exploratory) | 6 | -0.18703 | -0.06548 | +0.12155 | 0.0034 | +0.037 | False | +0.146 | False | True |
| AFF4 | 6 | -0.18525 | -0.07553 | +0.10972 | 0.0615 | +0.023 | False | -0.000 | False | False |
| ZNF574 | 6 | -0.18191 | -0.04613 | +0.13579 | 0.0034 | -0.018 | True | +0.027 | False | False |

**V_sig / S2_k100 / O2_Y2** — 18 of these 20 are `toward_young`; 85 in this setting overall.

| factor | n_guides | delta | floor_p05 | margin_over_floor | q_delta | guide_median_cos | guides_disagree | identity_drop | identity_loss | toward_young |
|---|---|---|---|---|---|---|---|---|---|---|
| HELZ2 (exploratory) | 6 | -1.67974 | -0.38865 | +1.29109 | 0.0033 | +0.324 | False | -0.081 | False | True |
| KLF5 (exploratory) | 6 | -0.61752 | -0.17580 | +0.44172 | 0.0033 | +0.166 | False | +0.126 | False | True |
| GFI1 (exploratory) | 6 | -0.45190 | -0.15509 | +0.29681 | 0.0033 | +0.231 | False | +0.061 | False | True |
| BSX (exploratory) | 6 | -0.42978 | -0.18891 | +0.24087 | 0.0033 | +0.048 | False | +0.021 | False | True |
| HHEX (exploratory) | 6 | -0.31337 | -0.14949 | +0.16388 | 0.0087 | +0.044 | False | +0.053 | False | True |
| NKX3-1 (exploratory) | 6 | -0.29298 | -0.17097 | +0.12201 | 0.0367 | +0.163 | False | +0.061 | False | True |
| NCOA1 (exploratory) | 6 | -0.28311 | -0.08259 | +0.20052 | 0.0087 | +0.012 | False | +0.007 | False | True |
| CREB5 (exploratory) | 6 | -0.27269 | -0.02773 | +0.24496 | 0.0302 | +0.013 | False | -0.012 | False | True |
| ZEB2 (exploratory) | 6 | -0.27004 | -0.03025 | +0.23979 | 0.0033 | +0.067 | False | -0.208 | False | True |
| GATA1 (exploratory) | 6 | -0.26574 | -0.03888 | +0.22686 | 0.0033 | +0.047 | False | +0.019 | False | True |
| MYOG (exploratory) | 6 | -0.25887 | -0.03898 | +0.21989 | 0.0033 | +0.008 | False | -0.049 | False | True |
| BCL6B | 6 | -0.25682 | -0.09049 | +0.16633 | 0.0087 | -0.001 | True | -0.000 | False | False |
| PAX2 (exploratory) | 6 | -0.24912 | -0.06617 | +0.18295 | 0.0033 | +0.114 | False | -0.149 | False | True |
| PRDM1 (exploratory) | 6 | -0.24755 | -0.05548 | +0.19207 | 0.0096 | +0.099 | False | +0.073 | False | True |
| ZFP2 (exploratory) | 6 | -0.24627 | -0.07948 | +0.16679 | 0.0033 | +0.028 | False | +0.042 | False | True |
| SP140 (exploratory) | 6 | -0.22462 | -0.04431 | +0.18031 | 0.0236 | +0.006 | False | -0.061 | False | True |
| MNX1 (exploratory) | 6 | -0.21809 | -0.09235 | +0.12574 | 0.0074 | +0.087 | False | +0.015 | False | True |
| AFF4 | 6 | -0.19777 | -0.08348 | +0.11429 | 0.0538 | +0.023 | False | -0.000 | False | False |
| FOXD2 (exploratory) | 6 | -0.18805 | -0.06529 | +0.12276 | 0.0033 | +0.037 | False | +0.146 | False | True |
| CLOCK (exploratory) | 6 | -0.18461 | -0.05473 | +0.12988 | 0.0033 | +0.034 | False | -0.004 | False | True |

### What produced these counts

Several of the counts above read backwards without their mechanism, so each is stated.

**Why each of the 2 zero settings is zero.** **In neither is the q-value gate the binding constraint.** Both are zero because no factor has a negative delta at all, which is a geometric fact about step length and aim and is settled before significance is consulted. Under SOUTH2's n_perm = 200 the q gate was the binding constraint in 18 of its 19 zero settings; with n_perm = 20,000 it binds in none. That is change 1 doing exactly what it was for. Observed ties at the p floor per setting are in the `at_p_floor` column of the count table and in `stage3_full_table.csv`.

| variant | space | pair | binding_constraint |
|---|---|---|---|
| V_all | S1 | O1_Y1 | no factor has a negative delta |
| V_sig | S1 | O1_Y1 | no factor has a negative delta |

**In 2 of 28 settings not one factor has a negative delta, so `toward_young` = 0 there regardless of significance.** Adding a displacement of norm ||d|| to a point at distance ||v|| from the target reduces that distance only if cos(d, v) > ||d|| / (2||v||); a large step needs to be well aimed or it overshoots. Below is that required cosine against the cosine actually achieved. Meeting one's own requirement is the same event as delta < 0, so the count column restates the result; the informative comparison is how far the best achieved cosine sits from the required one:

| variant | space | pair | origin_to_target | median_factor_norm | median_required_cos | smallest_required_cos | observed_max_cos | n_factors_meeting_own_requirement | smallest_delta |
|---|---|---|---|---|---|---|---|---|---|
| V_all | S1 | O1_Y1 | 32.284 | 20.254 | 0.3137 | 0.1709 | 0.0974 | 0 | +1.5303 |
| V_sig | S1 | O1_Y1 | 32.284 | 3.289 | 0.0509 | 0.0115 | 0.0535 | 0 | +0.0049 |

In those settings the `below_own_floor` count is still non-zero. That is consistent: a factor can move *less far away* than a random direction of its own size would, which clears its floor while still increasing the distance. It is not an approach, and the pre-registered `toward_young` rule requires delta < 0, so none of them is labelled a candidate.

**The BH-FDR gate is no longer resolution-limited — that is CHANGE 1.** The smallest p an empirical one-sided permutation test can return is 1/(20,000+1) = **5.000e-05**. With BH across 1,836 factors, q <= 0.05 needs only **2 factors tied at that minimum p** in the same setting, against 183 under SOUTH2's n_perm = 200. The most any setting shows at the p floor is 140 factors. So a `q_le_005` of 0 in SOUTH3 is a statement about the data, not about the null's resolution — which is exactly the sentence SOUTH2 could not write.

**The V_sig guide-agreement collapse is now reported and not gated on — that is CHANGE 2.** 1,692 of 1,836 factors have a median pairwise guide cosine of **exactly 0** when it is computed on V_sig vectors, 0 are negative, and 143 are positive. A cosine of exactly 0 there means two guides for the same factor share **no significant gene at all**, so their masked vectors are orthogonal by construction — an artefact of per-guide `adj_p <= 0.05` masking, not a disagreement about direction. SOUTH2 gated V_sig on that number and so refused 1,692 factors for a reason that was not about reliability. SOUTH3 computes the agreement on the V_all vectors always and applies that one value to both variants, so `guides_disagree` fires on 299 factors in V_sig as well, the same set as in V_all. The V_sig number is kept as a reported-only column.

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
| V_all | S2_k20 | O1_Y2 | HELZ2 | -2.8536 | 33.210 | +0.13236 | +0.4180 | True |
| V_all | S2_k20 | O1_Y1 | NCL | -2.7385 | 28.645 | +0.11444 | +0.5289 | True |
| V_all | S1 | O2_Y1 | HELZ2 | -2.6565 | 468.388 | +0.00903 | +0.1099 | True |
| V_all | S2_k50 | O1_Y1 | PRDM1 | -2.3883 | 30.864 | +0.12562 | +0.3924 | True |
| V_all | S2_k50 | O1_Y2 | HELZ2 | -2.1785 | 35.168 | +0.12320 | +0.3466 | True |
| V_all | S2_k100 | O1_Y1 | PRDM1 | -1.9680 | 31.434 | +0.12149 | +0.3483 | True |
| V_sig | S2_k20 | O1_Y2 | HELZ2 | -1.7690 | 33.210 | +0.06971 | +0.3688 | True |
| V_sig | S2_k20 | O2_Y2 | HELZ2 | -1.6910 | 204.376 | +0.00871 | +0.2836 | True |
| V_sig | S2_k50 | O2_Y2 | HELZ2 | -1.6897 | 265.189 | +0.00678 | +0.2321 | True |
| V_sig | S2_k100 | O2_Y2 | HELZ2 | -1.6797 | 277.476 | +0.00658 | +0.1998 | True |
| V_sig | S2_k50 | O2_Y1 | HELZ2 | -1.5911 | 259.334 | +0.00656 | +0.2197 | True |
| V_sig | S2_k100 | O2_Y1 | HELZ2 | -1.5845 | 272.045 | +0.00637 | +0.1898 | True |
| V_sig | S2_k20 | O2_Y1 | HELZ2 | -1.5627 | 196.450 | +0.00843 | +0.2639 | True |
| V_all | S2_k100 | O1_Y2 | HELZ2 | -1.5553 | 35.835 | +0.11921 | +0.3043 | True |
| V_sig | S1 | O2_Y1 | HELZ2 | -1.4761 | 468.388 | +0.00462 | +0.0851 | True |
| V_sig | S2_k50 | O1_Y2 | HELZ2 | -1.4733 | 35.168 | +0.06526 | +0.2964 | True |
| V_sig | S2_k100 | O1_Y2 | HELZ2 | -1.1057 | 35.835 | +0.06285 | +0.2466 | True |
| V_sig | S2_k20 | O1_Y1 | HELZ2 | -0.6636 | 28.645 | +0.04691 | +0.2141 | False |
| V_sig | S2_k50 | O1_Y1 | HELZ2 | -0.5033 | 30.864 | +0.04765 | +0.1899 | False |
| V_sig | S2_k100 | O1_Y1 | BSX | -0.4381 | 31.434 | +0.01608 | +0.2403 | True |
| V_sig | S1 | O1_Y1 | ZNF175 | +0.0049 | 32.284 | +0.00028 | +0.0095 | False |
| V_all | S1 | O1_Y1 | REPIN1 | +1.5303 | 32.284 | +0.01592 | +0.0443 | False |

## Stage 4 disagreement with directional readouts

All factors are ranked three ways in every setting: by delta (most negative first), by cos (highest first), and by drop in the frozen ruler score (largest drop first). Spearman correlations are between those rank vectors, so positive means agreement.

| variant | space | pair | rho_delta_cos | rho_delta_ruler | rho_cos_ruler | top20_overlap_delta_cos | top20_overlap_delta_ruler | top20_overlap_cos_ruler | top20_by_cos_not_toward_young |
|---|---|---|---|---|---|---|---|---|---|
| V_all | S1 | O1_Y1 | -0.008 | +0.094 | +0.293 | 1 | 0 | 6 | 20 |
| V_all | S1 | O2_Y1 | +0.974 | +0.321 | +0.318 | 14 | 3 | 2 | 1 |
| V_all | S2_k20 | O1_Y1 | +0.919 | +0.166 | +0.155 | 1 | 4 | 0 | 2 |
| V_all | S2_k20 | O1_Y2 | +0.974 | +0.325 | +0.310 | 7 | 3 | 2 | 4 |
| V_all | S2_k20 | O2_Y1 | +0.970 | +0.323 | +0.315 | 4 | 2 | 0 | 8 |
| V_all | S2_k20 | O2_Y2 | +0.967 | +0.309 | +0.308 | 5 | 2 | 0 | 6 |
| V_all | S2_k50 | O1_Y1 | +0.951 | +0.235 | +0.229 | 3 | 5 | 1 | 3 |
| V_all | S2_k50 | O1_Y2 | +0.972 | +0.352 | +0.348 | 6 | 2 | 1 | 4 |
| V_all | S2_k50 | O2_Y1 | +0.974 | +0.333 | +0.328 | 5 | 4 | 0 | 9 |
| V_all | S2_k50 | O2_Y2 | +0.973 | +0.326 | +0.324 | 5 | 3 | 0 | 9 |
| V_all | S2_k100 | O1_Y1 | +0.972 | +0.259 | +0.259 | 6 | 5 | 2 | 4 |
| V_all | S2_k100 | O1_Y2 | +0.943 | +0.360 | +0.372 | 10 | 2 | 1 | 5 |
| V_all | S2_k100 | O2_Y1 | +0.969 | +0.341 | +0.327 | 3 | 4 | 0 | 7 |
| V_all | S2_k100 | O2_Y2 | +0.973 | +0.336 | +0.324 | 3 | 4 | 0 | 8 |
| V_sig | S1 | O1_Y1 | +0.157 | +0.096 | +0.248 | 0 | 0 | 1 | 20 |
| V_sig | S1 | O2_Y1 | +0.881 | +0.097 | +0.101 | 15 | 0 | 0 | 2 |
| V_sig | S2_k20 | O1_Y1 | +0.871 | +0.073 | +0.067 | 9 | 2 | 1 | 2 |
| V_sig | S2_k20 | O1_Y2 | +0.923 | +0.084 | +0.074 | 7 | 1 | 0 | 5 |
| V_sig | S2_k20 | O2_Y1 | +0.947 | +0.103 | +0.098 | 1 | 2 | 0 | 13 |
| V_sig | S2_k20 | O2_Y2 | +0.948 | +0.097 | +0.094 | 0 | 3 | 0 | 12 |
| V_sig | S2_k50 | O1_Y1 | +0.880 | +0.157 | +0.153 | 9 | 2 | 1 | 3 |
| V_sig | S2_k50 | O1_Y2 | +0.922 | +0.143 | +0.140 | 9 | 3 | 1 | 3 |
| V_sig | S2_k50 | O2_Y1 | +0.952 | +0.129 | +0.133 | 4 | 2 | 0 | 4 |
| V_sig | S2_k50 | O2_Y2 | +0.952 | +0.125 | +0.128 | 3 | 2 | 0 | 4 |
| V_sig | S2_k100 | O1_Y1 | +0.874 | +0.190 | +0.190 | 9 | 2 | 2 | 2 |
| V_sig | S2_k100 | O1_Y2 | +0.922 | +0.173 | +0.172 | 9 | 4 | 1 | 3 |
| V_sig | S2_k100 | O2_Y1 | +0.947 | +0.117 | +0.119 | 5 | 2 | 0 | 2 |
| V_sig | S2_k100 | O2_Y2 | +0.949 | +0.114 | +0.116 | 6 | 2 | 0 | 2 |

The last column is the number asked for directly: of the 20 best factors by cosine, how many are **not** `toward_young` by delta.

## CHANGE 2 — one guide-agreement number, computed on V_all, applied to both variants

SOUTH2 computed the median pairwise guide cosine separately inside each gene variant and gated each variant on its own value. In V_sig that number is **exactly 0** for 1,692 of 1,836 factors, because per-guide `adj_p <= 0.05` masking often leaves two guides of the same factor with no significant gene in common, so their vectors are orthogonal by construction rather than by disagreement. SOUTH3 computes the agreement on the V_all vectors always and applies that one number to both variants.

**How many factors change label because of this.** Two isolations are given, because only both together are honest: change 2 evaluated inside SOUTH2's own delta / floor / q, and change 2 evaluated inside SOUTH3's. Both count factor-settings whose `toward_young` label flips when only the guide rule is swapped.

| variant | space | pair | guides_disagree_S2 | guides_disagree_S3 | flag_changed | toward_young_S2 | toward_young_S3 | label_changed_total | label_changed_by_change2_in_S2 | label_changed_by_change2_in_S3 |
|---|---|---|---|---|---|---|---|---|---|---|
| V_all | S1 | O1_Y1 | 299 | 299 | 0 | 0 | 0 | 0 | 0 | 0 |
| V_all | S1 | O2_Y1 | 299 | 299 | 0 | 231 | 266 | 37 | 0 | 0 |
| V_all | S2_k100 | O1_Y1 | 299 | 299 | 0 | 0 | 173 | 173 | 0 | 0 |
| V_all | S2_k100 | O1_Y2 | 299 | 299 | 0 | 0 | 40 | 40 | 0 | 0 |
| V_all | S2_k100 | O2_Y1 | 299 | 299 | 0 | 216 | 248 | 42 | 0 | 0 |
| V_all | S2_k100 | O2_Y2 | 299 | 299 | 0 | 208 | 243 | 39 | 0 | 0 |
| V_all | S2_k20 | O1_Y1 | 299 | 299 | 0 | 224 | 277 | 69 | 0 | 0 |
| V_all | S2_k20 | O1_Y2 | 299 | 299 | 0 | 0 | 44 | 44 | 0 | 0 |
| V_all | S2_k20 | O2_Y1 | 299 | 299 | 0 | 189 | 219 | 36 | 0 | 0 |
| V_all | S2_k20 | O2_Y2 | 299 | 299 | 0 | 190 | 208 | 28 | 0 | 0 |
| V_all | S2_k50 | O1_Y1 | 299 | 299 | 0 | 194 | 247 | 57 | 0 | 0 |
| V_all | S2_k50 | O1_Y2 | 299 | 299 | 0 | 0 | 53 | 53 | 0 | 0 |
| V_all | S2_k50 | O2_Y1 | 299 | 299 | 0 | 207 | 236 | 37 | 0 | 0 |
| V_all | S2_k50 | O2_Y2 | 299 | 299 | 0 | 194 | 230 | 42 | 0 | 0 |
| V_sig | S1 | O1_Y1 | 1,692 | 299 | 1,395 | 0 | 0 | 0 | 0 | 0 |
| V_sig | S1 | O2_Y1 | 1,692 | 299 | 1,395 | 0 | 128 | 128 | 0 | 107 |
| V_sig | S2_k100 | O1_Y1 | 1,692 | 299 | 1,395 | 0 | 48 | 48 | 0 | 36 |
| V_sig | S2_k100 | O1_Y2 | 1,692 | 299 | 1,395 | 0 | 40 | 40 | 0 | 29 |
| V_sig | S2_k100 | O2_Y1 | 1,692 | 299 | 1,395 | 0 | 80 | 80 | 0 | 65 |
| V_sig | S2_k100 | O2_Y2 | 1,692 | 299 | 1,395 | 0 | 85 | 85 | 0 | 68 |
| V_sig | S2_k20 | O1_Y1 | 1,692 | 299 | 1,395 | 0 | 67 | 67 | 0 | 47 |
| V_sig | S2_k20 | O1_Y2 | 1,692 | 299 | 1,395 | 0 | 45 | 45 | 0 | 31 |
| V_sig | S2_k20 | O2_Y1 | 1,692 | 299 | 1,395 | 0 | 40 | 40 | 0 | 28 |
| V_sig | S2_k20 | O2_Y2 | 1,692 | 299 | 1,395 | 0 | 45 | 45 | 0 | 33 |
| V_sig | S2_k50 | O1_Y1 | 1,692 | 299 | 1,395 | 0 | 62 | 62 | 0 | 44 |
| V_sig | S2_k50 | O1_Y2 | 1,692 | 299 | 1,395 | 0 | 43 | 43 | 0 | 31 |
| V_sig | S2_k50 | O2_Y1 | 1,692 | 299 | 1,395 | 0 | 79 | 79 | 0 | 63 |
| V_sig | S2_k50 | O2_Y2 | 1,692 | 299 | 1,395 | 0 | 80 | 80 | 0 | 63 |

**The direct answer.** The `guides_disagree` flag changes on **19,530** of 51,408 factor-settings — every one of them in V_sig, since the V_all flag is by definition unchanged. The number of factors whose `toward_young` label changes **because of this rule alone** is **0** when evaluated against SOUTH2's own numbers and **645** when evaluated against SOUTH3's. The total SOUTH2 -> SOUTH3 label change, from all three patches together, is **1,539** factor-settings.

Read the first of those two numbers with its cause: under SOUTH2's n_perm = 200 no V_sig setting had a single factor at q <= 0.05, so the V_sig guide flag was never the binding constraint there and relaxing it could not by itself create or destroy a label. **0 is therefore not evidence that the rule does not matter** — it is evidence that SOUTH2's null was too coarse for the rule to matter. The rule becomes load-bearing exactly once the null can resolve q, which is what change 1 does, and that is why the two isolations differ so sharply.

**Change 2 is what makes V_sig reportable at all.** V_sig has 842 `toward_young` factor-settings in SOUTH3, and 645 of them would not be labelled under SOUTH2's per-variant guide rule — they were being refused for sharing no significant gene between guides, which is an artefact of the masking and not a disagreement about direction. In V_all the change alters 0 labels, as it must: the V_all agreement value is the same number in both analyses.

Per-factor rows for this comparison are in `results/south3/change2_label_comparison_per_factor.csv`.

## CHANGE 3 — report-only displacement geometry (no gate)

Nothing in this section gates, labels, or filters anything. It is here because the size of a displacement relative to the gap it has to close determines whether delta and cosine can disagree at all.

### The identity

Write r = origin - target, so ||r|| = ||target - origin||; d for the factor's displacement; rho = ||d|| / ||r||; and theta for the angle between d and -r, so that cos(theta) = -d.r / (||d|| ||r||) is exactly the `cos` column. Then, with no approximation,

```
delta = ||r|| * ( sqrt(1 - 2*rho*cos(theta) + rho^2) - 1 )
```

Expanding in rho:

```
delta = -||d||*cos(theta)  +  ||d||^2 * (1 - cos^2(theta)) / (2*||r||)  +  O(rho^3)
```

So **delta ~= -||d|| cos(theta)** is the small-displacement identity. Three things follow and are worth stating before the numbers: it is exact only in the limit rho -> 0; the leading correction is **positive** whenever |cos| < 1, so the identity **overstates** how much ground a real displacement covers; and because the correction depends on ||d|| as well as on the angle, ranking by delta and ranking by cosine must coincide as rho -> 0 and can only come apart at larger rho.

The residual is therefore quoted as a fraction of **||d||**, not of |delta|. Dividing by |delta| looks natural and is useless here: delta passes through zero, so the relative error is unbounded near it and its large values say nothing about the identity. Divided by ||d||, the expansion predicts a residual of exactly rho*(1 - cos^2)/2, which is bounded above by rho/2.

### Do the data follow the identity?

**Yes, to the leading term and with the predicted sign.** The residual delta - (-||d||cos(theta)) is non-negative on **all 51,408 factor-settings**, as the expansion requires. Dividing the residual by the leading term the expansion predicts for it, ||d||*rho*(1 - cos^2)/2, gives per-setting medians of 0.9317 to 1.0223 — that is, the second-order term is not merely the right order of magnitude, it is the right number.

**Whether the identity is a usable approximation, though, varies by two orders of magnitude across settings, and the ratio says which.** Where rho stays small the residual is negligible: across the 14 O2 settings rho never exceeds 0.250 and the median residual is under 2.2% of ||d||, because the cross-dataset offset makes ||target - origin|| 6 to 15 times larger than it is for O1 in the same space. One setting sits outside the small-displacement regime entirely: **V_all / S1 / O1_Y1**, where the median factor moves 0.63 of the whole origin-to-target gap and the largest moves 3.63 times it, leaving a median residual of 29% of ||d||. That is the same fact as the Stage 3 result that **no factor has a negative delta there**: when the step is as long as the journey the quadratic term dominates, and no achievable angle is enough to end up closer.

The other zero setting, **V_sig / S1 / O1_Y1**, is zero for a different reason and the ratio is the wrong explanation for it: its rho median is 0.102, no larger than several settings that do produce survivors. What it has instead is cosines of essentially zero: significance-masked displacements in gene space are close to orthogonal to the within-GTEx age direction, and the required cosine rho/2 is not met by any factor. Both zero settings are S1 / O1_Y1 — the one origin-target pair whose gap is a real within-GTEx age difference rather than a cross-dataset offset — but one fails on step length and the other on aim, and they should not be given the same explanation.

### ||d|| and the ratio ||d|| / ||target - origin||, per setting

| variant | space | pair | origin_to_target | d_norm_min | d_norm_median | d_norm_max | ratio_min | ratio_median | ratio_max |
|---|---|---|---|---|---|---|---|---|---|
| V_all | S1 | O1_Y1 | 32.284 | 11.0365 | 20.2543 | 117.2566 | 0.34186 | 0.62738 | 3.63204 |
| V_all | S1 | O2_Y1 | 468.388 | 11.0365 | 20.2543 | 117.2566 | 0.02356 | 0.04324 | 0.25034 |
| V_all | S2_k20 | O1_Y1 | 28.645 | 0.6620 | 2.6017 | 10.5153 | 0.02311 | 0.09083 | 0.36709 |
| V_all | S2_k20 | O1_Y2 | 33.210 | 0.6620 | 2.6017 | 10.5153 | 0.01993 | 0.07834 | 0.31663 |
| V_all | S2_k20 | O2_Y1 | 196.450 | 0.6620 | 2.6017 | 10.5153 | 0.00337 | 0.01324 | 0.05353 |
| V_all | S2_k20 | O2_Y2 | 204.376 | 0.6620 | 2.6017 | 10.5153 | 0.00324 | 0.01273 | 0.05145 |
| V_all | S2_k50 | O1_Y1 | 30.864 | 0.8844 | 3.0911 | 12.8166 | 0.02865 | 0.10015 | 0.41526 |
| V_all | S2_k50 | O1_Y2 | 35.168 | 0.8844 | 3.0911 | 12.8166 | 0.02515 | 0.08790 | 0.36444 |
| V_all | S2_k50 | O2_Y1 | 259.334 | 0.8844 | 3.0911 | 12.8166 | 0.00341 | 0.01192 | 0.04942 |
| V_all | S2_k50 | O2_Y2 | 265.189 | 0.8844 | 3.0911 | 12.8166 | 0.00334 | 0.01166 | 0.04833 |
| V_all | S2_k100 | O1_Y1 | 31.434 | 1.1146 | 3.5825 | 14.9721 | 0.03546 | 0.11397 | 0.47630 |
| V_all | S2_k100 | O1_Y2 | 35.835 | 1.1146 | 3.5825 | 14.9721 | 0.03110 | 0.09997 | 0.41780 |
| V_all | S2_k100 | O2_Y1 | 272.045 | 1.1146 | 3.5825 | 14.9721 | 0.00410 | 0.01317 | 0.05504 |
| V_all | S2_k100 | O2_Y2 | 277.476 | 1.1146 | 3.5825 | 14.9721 | 0.00402 | 0.01291 | 0.05396 |
| V_sig | S1 | O1_Y1 | 32.284 | 0.7435 | 3.2889 | 32.0358 | 0.02303 | 0.10187 | 0.99231 |
| V_sig | S1 | O2_Y1 | 468.388 | 0.7435 | 3.2889 | 32.0358 | 0.00159 | 0.00702 | 0.06840 |
| V_sig | S2_k20 | O1_Y1 | 28.645 | 0.0212 | 0.1626 | 6.2771 | 0.00074 | 0.00568 | 0.21913 |
| V_sig | S2_k20 | O1_Y2 | 33.210 | 0.0212 | 0.1626 | 6.2771 | 0.00064 | 0.00490 | 0.18901 |
| V_sig | S2_k20 | O2_Y1 | 196.450 | 0.0212 | 0.1626 | 6.2771 | 0.00011 | 0.00083 | 0.03195 |
| V_sig | S2_k20 | O2_Y2 | 204.376 | 0.0212 | 0.1626 | 6.2771 | 0.00010 | 0.00080 | 0.03071 |
| V_sig | S2_k50 | O1_Y1 | 30.864 | 0.0329 | 0.2233 | 7.7447 | 0.00107 | 0.00724 | 0.25093 |
| V_sig | S2_k50 | O1_Y2 | 35.168 | 0.0329 | 0.2233 | 7.7447 | 0.00093 | 0.00635 | 0.22022 |
| V_sig | S2_k50 | O2_Y1 | 259.334 | 0.0329 | 0.2233 | 7.7447 | 0.00013 | 0.00086 | 0.02986 |
| V_sig | S2_k50 | O2_Y2 | 265.189 | 0.0329 | 0.2233 | 7.7447 | 0.00012 | 0.00084 | 0.02920 |
| V_sig | S2_k100 | O1_Y1 | 31.434 | 0.0439 | 0.2751 | 9.1317 | 0.00140 | 0.00875 | 0.29050 |
| V_sig | S2_k100 | O1_Y2 | 35.835 | 0.0439 | 0.2751 | 9.1317 | 0.00122 | 0.00768 | 0.25482 |
| V_sig | S2_k100 | O2_Y1 | 272.045 | 0.0439 | 0.2751 | 9.1317 | 0.00016 | 0.00101 | 0.03357 |
| V_sig | S2_k100 | O2_Y2 | 277.476 | 0.0439 | 0.2751 | 9.1317 | 0.00016 | 0.00099 | 0.03291 |

Per-factor ||d|| and the ratio are the `d_norm` and `d_over_v` columns of `results/south3/stage3_full_table.csv`; this table is the per-setting summary in `results/south3/stage3_dnorm_ratio_summary.csv`.

### Rank agreement between delta and cosine, by quartile of the ratio

Within each setting the factors are split at their own 25th, 50th and 75th percentiles of ||d|| / ||target - origin||. `rho_delta_cos` is the Spearman correlation between the rank by delta (most negative first) and the rank by cosine (highest first), so **+1 means the two orderings agree**. The last three columns test the identity inside that quartile: the median residual as a fraction of ||d||, the bound rho/2 that the expansion puts on it, and the residual divided by the leading term predicted for it.

| variant | space | pair | quartile | n | ratio_median | rho_delta_cos | err_over_dnorm_median | err_bound_half_ratio | resid_over_leading_term |
|---|---|---|---|---|---|---|---|---|---|
| V_all | S1 | O1_Y1 | Q1 | 459 | 0.51995 | +0.498 | 0.24628 | 0.28026 | 0.9513 |
| V_all | S1 | O1_Y1 | Q2 | 459 | 0.59209 | +0.652 | 0.27717 | 0.31369 | 0.9373 |
| V_all | S1 | O1_Y1 | Q3 | 459 | 0.66520 | +0.579 | 0.30747 | 0.35837 | 0.9249 |
| V_all | S1 | O1_Y1 | Q4 | 459 | 0.80162 | +0.221 | 0.35876 | 1.81602 | 0.8933 |
| V_all | S1 | O2_Y1 | Q1 | 459 | 0.03584 | +0.999 | 0.01790 | 0.01932 | 1.0005 |
| V_all | S1 | O2_Y1 | Q2 | 459 | 0.04081 | +0.999 | 0.02039 | 0.02162 | 1.0005 |
| V_all | S1 | O2_Y1 | Q3 | 459 | 0.04585 | +0.999 | 0.02291 | 0.02470 | 1.0004 |
| V_all | S1 | O2_Y1 | Q4 | 459 | 0.05525 | +0.974 | 0.02760 | 0.12517 | 1.0005 |
| V_all | S2_k20 | O1_Y1 | Q1 | 459 | 0.05405 | +0.989 | 0.02573 | 0.03411 | 1.0072 |
| V_all | S2_k20 | O1_Y1 | Q2 | 459 | 0.07907 | +0.989 | 0.03722 | 0.04539 | 1.0192 |
| V_all | S2_k20 | O1_Y1 | Q3 | 459 | 0.10318 | +0.986 | 0.04818 | 0.05922 | 1.0280 |
| V_all | S2_k20 | O1_Y1 | Q4 | 459 | 0.13979 | +0.932 | 0.06583 | 0.18354 | 1.0436 |
| V_all | S2_k20 | O1_Y2 | Q1 | 459 | 0.04662 | +0.991 | 0.02255 | 0.02942 | 1.0011 |
| V_all | S2_k20 | O1_Y2 | Q2 | 459 | 0.06820 | +0.998 | 0.03350 | 0.03915 | 1.0038 |
| V_all | S2_k20 | O1_Y2 | Q3 | 459 | 0.08900 | +0.999 | 0.04390 | 0.05108 | 1.0077 |
| V_all | S2_k20 | O1_Y2 | Q4 | 459 | 0.12057 | +0.993 | 0.06000 | 0.15831 | 1.0123 |
| V_all | S2_k20 | O2_Y1 | Q1 | 459 | 0.00788 | +0.990 | 0.00374 | 0.00497 | 1.0000 |
| V_all | S2_k20 | O2_Y1 | Q2 | 459 | 0.01153 | +0.998 | 0.00558 | 0.00662 | 0.9998 |
| V_all | S2_k20 | O2_Y1 | Q3 | 459 | 0.01505 | +0.999 | 0.00735 | 0.00864 | 0.9999 |
| V_all | S2_k20 | O2_Y1 | Q4 | 459 | 0.02038 | +0.992 | 0.01003 | 0.02676 | 0.9997 |
| V_all | S2_k20 | O2_Y2 | Q1 | 459 | 0.00758 | +0.991 | 0.00359 | 0.00478 | 1.0000 |
| V_all | S2_k20 | O2_Y2 | Q2 | 459 | 0.01108 | +0.998 | 0.00535 | 0.00636 | 0.9995 |
| V_all | S2_k20 | O2_Y2 | Q3 | 459 | 0.01446 | +0.998 | 0.00705 | 0.00830 | 0.9996 |
| V_all | S2_k20 | O2_Y2 | Q4 | 459 | 0.01959 | +0.990 | 0.00962 | 0.02573 | 0.9991 |
| V_all | S2_k50 | O1_Y1 | Q1 | 459 | 0.06263 | +0.993 | 0.03035 | 0.03818 | 1.0059 |
| V_all | S2_k50 | O1_Y1 | Q2 | 459 | 0.08766 | +0.994 | 0.04213 | 0.05005 | 1.0163 |
| V_all | S2_k50 | O1_Y1 | Q3 | 459 | 0.11318 | +0.993 | 0.05459 | 0.06441 | 1.0235 |
| V_all | S2_k50 | O1_Y1 | Q4 | 459 | 0.15339 | +0.976 | 0.07351 | 0.20763 | 1.0343 |
| V_all | S2_k50 | O1_Y2 | Q1 | 459 | 0.05497 | +0.994 | 0.02712 | 0.03351 | 1.0011 |
| V_all | S2_k50 | O1_Y2 | Q2 | 459 | 0.07693 | +0.999 | 0.03797 | 0.04393 | 1.0037 |
| V_all | S2_k50 | O1_Y2 | Q3 | 459 | 0.09933 | +0.999 | 0.04924 | 0.05653 | 1.0068 |
| V_all | S2_k50 | O1_Y2 | Q4 | 459 | 0.13462 | +0.983 | 0.06673 | 0.18222 | 1.0096 |
| V_all | S2_k50 | O2_Y1 | Q1 | 459 | 0.00745 | +0.993 | 0.00360 | 0.00454 | 1.0000 |
| V_all | S2_k50 | O2_Y1 | Q2 | 459 | 0.01043 | +0.999 | 0.00506 | 0.00596 | 0.9999 |
| V_all | S2_k50 | O2_Y1 | Q3 | 459 | 0.01347 | +0.999 | 0.00660 | 0.00767 | 1.0000 |
| V_all | S2_k50 | O2_Y1 | Q4 | 459 | 0.01826 | +0.994 | 0.00903 | 0.02471 | 1.0002 |
| V_all | S2_k50 | O2_Y2 | Q1 | 459 | 0.00729 | +0.993 | 0.00350 | 0.00444 | 0.9999 |
| V_all | S2_k50 | O2_Y2 | Q2 | 459 | 0.01020 | +0.998 | 0.00495 | 0.00583 | 0.9998 |
| V_all | S2_k50 | O2_Y2 | Q3 | 459 | 0.01317 | +0.999 | 0.00645 | 0.00750 | 0.9999 |
| V_all | S2_k50 | O2_Y2 | Q4 | 459 | 0.01785 | +0.993 | 0.00882 | 0.02416 | 0.9999 |
| V_all | S2_k100 | O1_Y1 | Q1 | 459 | 0.07239 | +0.996 | 0.03558 | 0.04363 | 1.0051 |
| V_all | S2_k100 | O1_Y1 | Q2 | 459 | 0.09916 | +0.997 | 0.04853 | 0.05697 | 1.0135 |
| V_all | S2_k100 | O1_Y1 | Q3 | 459 | 0.12908 | +0.997 | 0.06321 | 0.07409 | 1.0205 |
| V_all | S2_k100 | O1_Y1 | Q4 | 459 | 0.17380 | +0.986 | 0.08520 | 0.23815 | 1.0290 |
| V_all | S2_k100 | O1_Y2 | Q1 | 459 | 0.06350 | +0.992 | 0.03123 | 0.03827 | 1.0006 |
| V_all | S2_k100 | O1_Y2 | Q2 | 459 | 0.08698 | +0.998 | 0.04300 | 0.04997 | 1.0029 |
| V_all | S2_k100 | O1_Y2 | Q3 | 459 | 0.11323 | +0.998 | 0.05623 | 0.06499 | 1.0054 |
| V_all | S2_k100 | O1_Y2 | Q4 | 459 | 0.15245 | +0.942 | 0.07584 | 0.20890 | 1.0069 |
| V_all | S2_k100 | O2_Y1 | Q1 | 459 | 0.00836 | +0.993 | 0.00406 | 0.00504 | 1.0002 |
| V_all | S2_k100 | O2_Y1 | Q2 | 459 | 0.01146 | +0.999 | 0.00560 | 0.00658 | 1.0004 |
| V_all | S2_k100 | O2_Y1 | Q3 | 459 | 0.01492 | +0.999 | 0.00737 | 0.00856 | 1.0006 |
| V_all | S2_k100 | O2_Y1 | Q4 | 459 | 0.02008 | +0.990 | 0.00997 | 0.02752 | 1.0011 |
| V_all | S2_k100 | O2_Y2 | Q1 | 459 | 0.00820 | +0.993 | 0.00398 | 0.00494 | 1.0001 |
| V_all | S2_k100 | O2_Y2 | Q2 | 459 | 0.01123 | +0.999 | 0.00549 | 0.00645 | 1.0003 |
| V_all | S2_k100 | O2_Y2 | Q3 | 459 | 0.01462 | +0.999 | 0.00722 | 0.00839 | 1.0005 |
| V_all | S2_k100 | O2_Y2 | Q4 | 459 | 0.01969 | +0.992 | 0.00976 | 0.02698 | 1.0009 |
| V_sig | S1 | O1_Y1 | Q1 | 459 | 0.06308 | +0.496 | 0.03150 | 0.03883 | 0.9986 |
| V_sig | S1 | O1_Y1 | Q2 | 459 | 0.09037 | +0.743 | 0.04508 | 0.05094 | 0.9974 |
| V_sig | S1 | O1_Y1 | Q3 | 459 | 0.11642 | +0.613 | 0.05798 | 0.06904 | 0.9958 |
| V_sig | S1 | O1_Y1 | Q4 | 459 | 0.17705 | +0.216 | 0.08764 | 0.49616 | 0.9909 |
| V_sig | S1 | O2_Y1 | Q1 | 459 | 0.00435 | +0.976 | 0.00217 | 0.00268 | 1.0000 |
| V_sig | S1 | O2_Y1 | Q2 | 459 | 0.00623 | +0.995 | 0.00311 | 0.00351 | 1.0000 |
| V_sig | S1 | O2_Y1 | Q3 | 459 | 0.00802 | +0.992 | 0.00401 | 0.00476 | 0.9999 |
| V_sig | S1 | O2_Y1 | Q4 | 459 | 0.01220 | +0.878 | 0.00610 | 0.03420 | 0.9999 |
| V_sig | S2_k20 | O1_Y1 | Q1 | 459 | 0.00291 | +0.947 | 0.00134 | 0.00193 | 0.9996 |
| V_sig | S2_k20 | O1_Y1 | Q2 | 459 | 0.00475 | +0.992 | 0.00224 | 0.00284 | 0.9992 |
| V_sig | S2_k20 | O1_Y1 | Q3 | 459 | 0.00684 | +0.993 | 0.00322 | 0.00433 | 0.9990 |
| V_sig | S2_k20 | O1_Y1 | Q4 | 459 | 0.01280 | +0.939 | 0.00582 | 0.10957 | 0.9984 |
| V_sig | S2_k20 | O1_Y2 | Q1 | 459 | 0.00251 | +0.968 | 0.00119 | 0.00167 | 0.9998 |
| V_sig | S2_k20 | O1_Y2 | Q2 | 459 | 0.00410 | +0.995 | 0.00196 | 0.00245 | 0.9996 |
| V_sig | S2_k20 | O1_Y2 | Q3 | 459 | 0.00590 | +0.996 | 0.00282 | 0.00373 | 0.9995 |
| V_sig | S2_k20 | O1_Y2 | Q4 | 459 | 0.01104 | +0.961 | 0.00515 | 0.09451 | 0.9991 |
| V_sig | S2_k20 | O2_Y1 | Q1 | 459 | 0.00042 | +0.984 | 0.00020 | 0.00028 | 1.0000 |
| V_sig | S2_k20 | O2_Y1 | Q2 | 459 | 0.00069 | +0.998 | 0.00033 | 0.00041 | 1.0000 |
| V_sig | S2_k20 | O2_Y1 | Q3 | 459 | 0.00100 | +0.997 | 0.00048 | 0.00063 | 1.0000 |
| V_sig | S2_k20 | O2_Y1 | Q4 | 459 | 0.00187 | +0.968 | 0.00090 | 0.01598 | 1.0000 |
| V_sig | S2_k20 | O2_Y2 | Q1 | 459 | 0.00041 | +0.984 | 0.00019 | 0.00027 | 1.0000 |
| V_sig | S2_k20 | O2_Y2 | Q2 | 459 | 0.00067 | +0.998 | 0.00032 | 0.00040 | 1.0000 |
| V_sig | S2_k20 | O2_Y2 | Q3 | 459 | 0.00096 | +0.997 | 0.00046 | 0.00061 | 1.0000 |
| V_sig | S2_k20 | O2_Y2 | Q4 | 459 | 0.00179 | +0.970 | 0.00086 | 0.01536 | 1.0000 |
| V_sig | S2_k50 | O1_Y1 | Q1 | 459 | 0.00387 | +0.954 | 0.00188 | 0.00254 | 0.9996 |
| V_sig | S2_k50 | O1_Y1 | Q2 | 459 | 0.00609 | +0.995 | 0.00299 | 0.00362 | 0.9994 |
| V_sig | S2_k50 | O1_Y1 | Q3 | 459 | 0.00872 | +0.991 | 0.00426 | 0.00542 | 0.9991 |
| V_sig | S2_k50 | O1_Y1 | Q4 | 459 | 0.01537 | +0.949 | 0.00740 | 0.12546 | 0.9987 |
| V_sig | S2_k50 | O1_Y2 | Q1 | 459 | 0.00339 | +0.969 | 0.00165 | 0.00223 | 0.9998 |
| V_sig | S2_k50 | O1_Y2 | Q2 | 459 | 0.00535 | +0.996 | 0.00262 | 0.00317 | 0.9996 |
| V_sig | S2_k50 | O1_Y2 | Q3 | 459 | 0.00766 | +0.995 | 0.00375 | 0.00476 | 0.9994 |
| V_sig | S2_k50 | O1_Y2 | Q4 | 459 | 0.01349 | +0.961 | 0.00659 | 0.11011 | 0.9992 |
| V_sig | S2_k50 | O2_Y1 | Q1 | 459 | 0.00046 | +0.983 | 0.00022 | 0.00030 | 1.0000 |
| V_sig | S2_k50 | O2_Y1 | Q2 | 459 | 0.00073 | +0.998 | 0.00036 | 0.00043 | 1.0000 |
| V_sig | S2_k50 | O2_Y1 | Q3 | 459 | 0.00104 | +0.997 | 0.00051 | 0.00065 | 1.0000 |
| V_sig | S2_k50 | O2_Y1 | Q4 | 459 | 0.00183 | +0.975 | 0.00090 | 0.01493 | 1.0000 |
| V_sig | S2_k50 | O2_Y2 | Q1 | 459 | 0.00045 | +0.982 | 0.00022 | 0.00030 | 1.0000 |
| V_sig | S2_k50 | O2_Y2 | Q2 | 459 | 0.00071 | +0.998 | 0.00035 | 0.00042 | 1.0000 |
| V_sig | S2_k50 | O2_Y2 | Q3 | 459 | 0.00102 | +0.997 | 0.00050 | 0.00063 | 1.0000 |
| V_sig | S2_k50 | O2_Y2 | Q4 | 459 | 0.00179 | +0.974 | 0.00088 | 0.01460 | 1.0000 |
| V_sig | S2_k100 | O1_Y1 | Q1 | 459 | 0.00476 | +0.952 | 0.00233 | 0.00313 | 0.9996 |
| V_sig | S2_k100 | O1_Y1 | Q2 | 459 | 0.00751 | +0.994 | 0.00370 | 0.00437 | 0.9993 |
| V_sig | S2_k100 | O1_Y1 | Q3 | 459 | 0.01062 | +0.991 | 0.00521 | 0.00661 | 0.9991 |
| V_sig | S2_k100 | O1_Y1 | Q4 | 459 | 0.01832 | +0.942 | 0.00898 | 0.14525 | 0.9986 |
| V_sig | S2_k100 | O1_Y2 | Q1 | 459 | 0.00418 | +0.967 | 0.00207 | 0.00275 | 0.9998 |
| V_sig | S2_k100 | O1_Y2 | Q2 | 459 | 0.00658 | +0.996 | 0.00324 | 0.00384 | 0.9996 |
| V_sig | S2_k100 | O1_Y2 | Q3 | 459 | 0.00932 | +0.995 | 0.00458 | 0.00580 | 0.9995 |
| V_sig | S2_k100 | O1_Y2 | Q4 | 459 | 0.01607 | +0.956 | 0.00793 | 0.12741 | 0.9991 |
| V_sig | S2_k100 | O2_Y1 | Q1 | 459 | 0.00055 | +0.982 | 0.00027 | 0.00036 | 1.0000 |
| V_sig | S2_k100 | O2_Y1 | Q2 | 459 | 0.00087 | +0.998 | 0.00043 | 0.00051 | 1.0000 |
| V_sig | S2_k100 | O2_Y1 | Q3 | 459 | 0.00123 | +0.997 | 0.00061 | 0.00076 | 1.0000 |
| V_sig | S2_k100 | O2_Y1 | Q4 | 459 | 0.00212 | +0.969 | 0.00104 | 0.01678 | 0.9999 |
| V_sig | S2_k100 | O2_Y2 | Q1 | 459 | 0.00054 | +0.982 | 0.00026 | 0.00035 | 1.0000 |
| V_sig | S2_k100 | O2_Y2 | Q2 | 459 | 0.00085 | +0.998 | 0.00042 | 0.00050 | 1.0000 |
| V_sig | S2_k100 | O2_Y2 | Q3 | 459 | 0.00120 | +0.996 | 0.00059 | 0.00075 | 1.0000 |
| V_sig | S2_k100 | O2_Y2 | Q4 | 459 | 0.00208 | +0.970 | 0.00102 | 0.01645 | 0.9999 |

Read down the `rho_delta_cos` column within any one setting and the pattern is the same everywhere: the two rankings agree almost perfectly in the middle quartiles and loosen at both ends — at the bottom because the deltas there are so small that ties and rounding dominate the ordering, at the top because that is where the quadratic term is large enough to re-order factors that the angle alone would have ranked differently. The only settings where the agreement is poor in every quartile are the two S1 / O1_Y1 settings (+0.216 to +0.743 across their eight quartiles), which are also the two with no survivor; every quartile of every other setting is at least +0.878.

Full per-quartile numbers, including the error relative to |delta| for completeness, are in `results/south3/stage3_geometry_quartiles.csv`.

## Is any survivor distinguishable from the ranking a plain cosine would give?

Stated plainly, because it is the question the pre-registration ends on. For every setting with at least one `toward_young` survivor, the survivor set is compared with the top-n factors under four one-number rankings, n being the number of survivors in that setting: cosine alone, the small-displacement identity -||d||cos(theta), delta itself, and ||d|| alone.

**`recovered_by_delta` is the control that makes this table readable.** The `toward_young` label is a conjunction — delta < 0 AND below the factor's own magnitude-matched floor AND q <= 0.05 AND neither flag — not a ranking. So even delta, the statistic the label is built on, cannot reproduce the survivor set exactly by taking its own top n. Cosine only fails in a way that means anything if it fails by more than delta does.

| variant | space | pair | n_survivors | recovered_by_cos | recovered_by_identity | recovered_by_delta | recovered_by_d_norm | cos_pct | delta_pct | worst_survivor_cos_rank | rho_delta_cos_all_factors |
|---|---|---|---|---|---|---|---|---|---|---|---|
| V_all | S1 | O2_Y1 | 266 | 184 | 165 | 180 | 36 | 69.2 | 67.7 | 886 | +0.974 |
| V_all | S2_k20 | O1_Y1 | 277 | 208 | 148 | 169 | 79 | 75.1 | 61.0 | 1,186 | +0.919 |
| V_all | S2_k20 | O1_Y2 | 44 | 32 | 12 | 21 | 4 | 72.7 | 47.7 | 262 | +0.974 |
| V_all | S2_k20 | O2_Y1 | 219 | 155 | 178 | 180 | 25 | 70.8 | 82.2 | 439 | +0.970 |
| V_all | S2_k20 | O2_Y2 | 208 | 144 | 165 | 167 | 21 | 69.2 | 80.3 | 417 | +0.967 |
| V_all | S2_k50 | O1_Y1 | 247 | 194 | 123 | 150 | 57 | 78.5 | 60.7 | 1,515 | +0.951 |
| V_all | S2_k50 | O1_Y2 | 53 | 35 | 17 | 31 | 3 | 66.0 | 58.5 | 352 | +0.972 |
| V_all | S2_k50 | O2_Y1 | 236 | 165 | 189 | 191 | 33 | 69.9 | 80.9 | 426 | +0.974 |
| V_all | S2_k50 | O2_Y2 | 230 | 159 | 185 | 185 | 30 | 69.1 | 80.4 | 468 | +0.973 |
| V_all | S2_k100 | O1_Y1 | 173 | 143 | 73 | 109 | 25 | 82.7 | 63.0 | 933 | +0.972 |
| V_all | S2_k100 | O1_Y2 | 40 | 27 | 6 | 22 | 1 | 67.5 | 55.0 | 167 | +0.943 |
| V_all | S2_k100 | O2_Y1 | 248 | 183 | 182 | 187 | 34 | 73.8 | 75.4 | 499 | +0.969 |
| V_all | S2_k100 | O2_Y2 | 243 | 178 | 185 | 191 | 34 | 73.3 | 78.6 | 790 | +0.973 |
| V_sig | S1 | O2_Y1 | 128 | 97 | 90 | 93 | 16 | 75.8 | 72.7 | 379 | +0.881 |
| V_sig | S2_k20 | O1_Y1 | 67 | 35 | 42 | 42 | 11 | 52.2 | 62.7 | 349 | +0.871 |
| V_sig | S2_k20 | O1_Y2 | 45 | 28 | 34 | 34 | 13 | 62.2 | 75.6 | 293 | +0.923 |
| V_sig | S2_k20 | O2_Y1 | 40 | 13 | 23 | 23 | 7 | 32.5 | 57.5 | 245 | +0.947 |
| V_sig | S2_k20 | O2_Y2 | 45 | 17 | 26 | 26 | 9 | 37.8 | 57.8 | 382 | +0.948 |
| V_sig | S2_k50 | O1_Y1 | 62 | 37 | 38 | 38 | 10 | 59.7 | 61.3 | 313 | +0.880 |
| V_sig | S2_k50 | O1_Y2 | 43 | 29 | 26 | 27 | 12 | 67.4 | 62.8 | 269 | +0.922 |
| V_sig | S2_k50 | O2_Y1 | 79 | 48 | 49 | 49 | 17 | 60.8 | 62.0 | 437 | +0.952 |
| V_sig | S2_k50 | O2_Y2 | 80 | 51 | 49 | 49 | 17 | 63.8 | 61.2 | 387 | +0.952 |
| V_sig | S2_k100 | O1_Y1 | 48 | 29 | 30 | 30 | 6 | 60.4 | 62.5 | 285 | +0.874 |
| V_sig | S2_k100 | O1_Y2 | 40 | 28 | 26 | 26 | 12 | 70.0 | 65.0 | 226 | +0.922 |
| V_sig | S2_k100 | O2_Y1 | 80 | 53 | 56 | 56 | 14 | 66.2 | 70.0 | 322 | +0.947 |
| V_sig | S2_k100 | O2_Y2 | 85 | 57 | 58 | 58 | 16 | 67.1 | 68.2 | 343 | +0.949 |

Two things to read off it. First, **neither direction statistic dominates**: cosine recovers more of the survivor set than delta does in 11 of the 26 settings and fewer in 15, and in aggregate they tie (70.0% against 70.2%). Per setting they differ by as much as 25 percentage points in either direction, which is the size of the disagreement you would expect between two near-identical orderings cut at an arbitrary rank; what does not happen anywhere is delta recovering the survivor set while cosine fails to. The shortfall in both columns is the conjunction of gates, not the choice of direction statistic. Second, `recovered_by_d_norm` is far lower than either, so the survivor set is not simply the largest displacements: the label does use the angle. What it does not use is anything the angle does not already carry.

**Verdict, stated plainly as the pre-registration asks.** **No.** Across the 26 settings that have survivors, the Spearman correlation between the delta ranking and the cosine ranking over all 1,836 factors is at least +0.871, and within every quartile of ||d|| / ||target - origin|| in those settings it is at least +0.878. Taking the top-n factors by cosine alone, n being that setting's number of survivors, recovers 2,329 of the 3,326 survivors (70.0%) — and taking the top-n by delta itself, the statistic the label is built on, recovers 2,334 (70.2%). The two agree to within a percentage point, so the ~30% that neither recovers is the conjunction of the floor and q gates, not information the cosine lacks. ||d|| alone recovers only 542 (16.3%), so the label is not merely picking the biggest displacements. The mechanism is the identity: in every setting that has a survivor the median residual of delta = -||d||cos(theta) is at most 5.6% of ||d||, so delta there is the cosine re-weighted by the displacement's own length and carries nothing beyond it. No survivor in any setting is distinguishable on direction from what a plain cosine ranking would give. What the delta machinery does buy is not a different ordering but a bar: it says which of those cosines are larger than the factor's own reshuffled null, and a cosine ranking on its own never says that.

## Fired keys

STOP keys fired: **none**.
STOP keys checked and not fired: `not_an_effect_matrix`, `too_few_factors`, `no_detection_power`.

REPORT-ONLY flags fired: `scale_unknown`, `coverage_low`.

- **scale_unknown** (report-only, not a stop): The units and log base of X are undeterminable from the file: uns is empty, X carries no unit attribute, var holds only mean/std/cv/fano/excess_cv with no unit statement, obs holds no unit statement, and the only bundled documentation in the folder is three plain gene-symbol lists. This is a REPORT-ONLY flag, not a stop. The Stage 2A/2B floors are magnitude-matched (they permute or re-use the same vector, so they are invariant to a global scale factor on X); the Stage 2C MDA is NOT invariant, because it scales the true young-minus-old direction to the median factor norm, and that norm moves with any global rescaling of X.

- **coverage_low** (report-only, not a stop): The chosen mapping (best_of_two_after_md3_idtype) matches 0.1228 of the frozen ruler's total absolute weight, which is below 0.25. The Stage 3 candidate table is titled exploratory and every candidate carries that word.

Recorded non-fatal failures and implementation decisions, verbatim:

- **Y2_singular:** Y2 at k=100: the young-donor covariance built on 56 fit donors has rank 55 < k=100, so the Mahalanobis metric is singular and undefined. No pseudo-inverse or ridge was substituted. Y2 is reported as undefined at this k; Y1 at the same k is unaffected.
- **stage2A_streaming:** Implementation decision forced by n_perm = 20,000 and fixed before any floor existed. SOUTH2 held every permuted delta in RAM (200 per factor and setting) and derived the floor and the permutation p from that array afterwards. 20,000 permuted deltas for 1,836 factors across 28 factor-settings do not fit, so each factor's permutations are drawn in row chunks of 1,000 and reduced as they are produced: the 5th percentile (the floor), the count of permuted deltas at or below the already-computed observed delta (the numerator of the one-sided empirical p), and the rescaled draw that feeds the pooled floor. The per-factor generator, np.random.default_rng([SEED, i]), the permuted coordinate set, and the delta definition are all unchanged. Two consequences are stated rather than hidden: (a) the observed delta is now computed before the permutations instead of after, which changes no value because it does not depend on them; (b) the pooled draw is evaluated from the algebraic form sqrt(r2 + 2*s*dot + s^2*qn2) - ||r|| instead of by rescaling a 20000 x 4907 block, which is the same quantity up to float64 rounding. Both were checked against the SOUTH2 code path before the run: drawing the permutations in chunks of 1,000 from np.random.default_rng([SEED, i]) gives a bit-identical block to drawing them in one call, so SOUTH3's null is SOUTH2's stream continued — the first 200 of SOUTH3's 20,000 draws are literally SOUTH2's 200 draws — and the streamed delta matches the block_stats delta exactly, the pooled draw to 5e-14 absolute. The reduced null — the per-factor floors and the observed-value counts, which is everything any later stage reads — is kept at `results/south3/stage2a_null_reduced.npz` so the report can be rebuilt without redrawing the permutations; it is keyed on the seed, n_perm, chunk size, factor count and setting list, and is ignored if any of them differs.
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
- **The q-value gate is no longer resolution-limited, but it is still a permutation gate.** With n_perm = 20,000 the smallest attainable p is 1/20,001 = 5.000e-05 and BH over 1,836 factors needs only 2 factors tied at that minimum before q can reach 0.05, so SOUTH2's resolution ceiling is gone. What remains is that the null is a within-factor coordinate permutation: it asks whether this factor's displacement is better aimed than a reshuffling of its own entries, and nothing else. It does not test measurement error in the effect estimates, it does not test transfer to an aged cell, and a small q is not a small effect size.
- **The permutations were drawn in chunks and reduced on the fly.** 20,000 permuted deltas per factor and setting cannot be held in memory, so the floor, the permutation count and the pooled draw are accumulated as the chunks are produced (see the `stage2A_streaming` note). The permuted values themselves are not retained, so the full null distributions are not recoverable from `results/south3/` — only their 5th percentiles, the observed-value counts, and the pooled floor.
- **The control-row floor is not a noise floor** and was never used as the gate: those rows are untargeted but not inert.
- This task makes **no claim about the source paper's own conclusions**. It did not test them.

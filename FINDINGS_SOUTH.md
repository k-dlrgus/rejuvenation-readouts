**SUPERSEDED BY SOUTH2 (see FINDINGS_SOUTH2.md): SOUTH halted at Stage 0 with `no_raw_counts` because its pre-registration described the pinned file as a count matrix; SOUTH2 re-registers the same file as the per-guide regression output it actually is. Nothing below is changed.**

# FINDINGS_SOUTH — rank Southard Hs27 CRISPRa transcription factors by approach to a young target

**Status:** STOP at Stage 0. Keys fired: `no_raw_counts`. No statistic was computed; Stages 1-4 did not run.

The pre-registration below was written verbatim to `FINDINGS_SOUTH.md` and `results/south/PREREG.flag` before any Southard row was opened and before any statistic.

## Pre-registration (verbatim, written before any SOUTH statistic)

```
TASK: SOUTH

Nature of this task: a purely computational reanalysis of published,
publicly released human cell gene-expression tables that are already
downloaded in this repository. No laboratory work, no protocol, no
biological material, no new data generation. The output is a ranked
table and a set of statistics.

Write FINDINGS_SOUTH.md, PROGRESS_SOUTH.md, results/south/*, src/south_run.py.
Do not modify any existing FINDINGS_*.md, PROGRESS_*.md, or existing src/ files.
Import from src/toward_run.py, src/same_run.py, src/seng_run.py,
src/lowdim_common.py; do not change them. Do not refit the frozen ruler.
Seeds: 20260914 (splits, nulls, permutations), 20260918 (bootstrap).
n_null=200, n_perm=200, n_boot=200. Uncalibrated R2 is never a gate.

WRITE THE PRE-REGISTRATION FIRST. Copy this whole file verbatim into
FINDINGS_SOUTH.md and touch results/south/PREREG.flag before computing any
statistic. Stage 0 may run after the flag exists. Nothing is re-tuned after
numbers exist. Every STOP key halts the task and is reported as written.

--- DATA ---
The already-downloaded expression table at
results/survey/opened/southard/fibroblast_CRISPRa_mean_pop.h5ad
(Southard et al. 2025, transcription-factor perturbations in the Hs27
human fibroblast line, public). Record file size and md5.

--- STAGE 0: INVENTORY (descriptive only, no statistic) ---
From obs/var, determine and report:
- what one row is: a single cell, a replicate, or a per-factor mean.
- how many transcription factors, how many replicate rows per factor,
  how many untargeted control rows, and their labels.
- whether the matrix holds raw integer counts; if not, whether a raw
  layer exists. If neither, STOP, key `no_raw_counts`.
- how many of the 23,485 frozen-ruler genes are present.
- what unit is available for resampling (cells, replicates, or none).
  If none, STOP, key `no_resampling_unit`, and report the inventory only.
STOP, key `no_controls`, if fewer than 10 control rows exist.

--- SPACES (every statistic is reported in both) ---
S1, gene space: frozen-ruler overlap genes; sum, TMM log2-CPM with
prior.count=2 across this task's own panel, frozen GTEx mu/sd, z-score;
missing genes z=0. Identical to TOWARD/SAME/SENG.
S2, low-dim space: PCA fit only on the GTEx fibroblast donor matrix,
never on this dataset and never on the targets, k in (20, 50, 100),
applied unchanged to everything else. Report every k. Do not choose a k
after seeing results.

--- TARGETS AND ORIGINS (fixed before anything is scored) ---
Young target, two definitions, reported separately, never combined:
  Y1 = GTEx fibroblast donors aged 20-39 (existing TOWARD anchor).
  Y2 = the same donors as a region: mean and covariance, in S2 only,
       distance measured as Mahalanobis distance, membership radius set
       at the 95th percentile of held-out young donors' own distances.
       Report the radius.
Aged origin, two definitions, reported separately, never combined:
  O1 = GTEx fibroblast donors aged 60-79 (existing anchor).
  O2 = the Lu et al. aged donor GM00731 untreated day-0 fibroblast
       pseudobulk, as already built in SAME.
Report each origin-to-target distance before anything else is scored.

--- STAGE 1: DETECTION LIMIT (before any factor is scored) ---
A. Control floor. Split the control rows into two disjoint halves,
   200 times, taking d_null = z(half1) - z(half2), matched in row count
   to the median real factor. For each origin/target pair report the
   distribution of
   delta = ||(origin + d_null) - target|| - ||origin - target||.
   The floor is the 5th percentile: the most negative delta that noise
   alone can reach. Report median, 5th percentile, and whether each
   bootstrap interval contains its point estimate.
B. Dose-response. Add the true young-minus-old direction to the aged
   origin at fractions f = 0, 0.05, 0.1, 0.25, 0.5, 1.0, with magnitude
   matched to the median real factor's vector norm. Report delta at each f.
   MDA = the smallest f whose delta falls below the control floor.
C. STOP, key `no_detection_power`, if no f <= 1.0 clears the floor in any
   space. Report Stages 0 and 1 only and state plainly that this data
   cannot support a ranking.
Never subtract the floor from a later number. It is a bar to clear.

--- STAGE 2: SINGLE-FACTOR RANKING (only if Stage 1 passed) ---
For each transcription factor g with at least 2 replicate rows, the effect
vector is d_g = z(g rows) - z(control rows). For each origin/target pair and
each space report:
  cos_g   = dot(d_g, v) / (||d_g|| ||v||), where v = target - origin
  frac_g  = dot(d_g, v) / ||v||^2
  delta_g = ||(origin + d_g) - target|| - ||origin - target||
  for Y2 also: whether origin + d_g falls inside the membership radius,
  and the change in Mahalanobis distance.
Also per factor: proliferation-residualized cos and delta (project out the
proliferation axis exactly as in seng_run.py); the identity check on
COL1A1, COL1A2, FN1, LUM, PDGFRA, PDGFRB, POSTN, PRRX1, SERPINH1, VIM,
flagging `identity_loss` if the drop exceeds 0.5; and the frozen ruler
score as a context column only.
Nulls: the Stage 1A control draws, size-matched to g, one-sided empirical p,
BH-FDR across all factors. Bootstrap intervals by the Stage 0 resampling
unit; flag any interval that does not contain its point estimate as INVALID
and never use it as a gate.
Label a factor `toward_young` only if delta_g is negative AND below the
Stage 1A floor AND q <= 0.05 AND the residualized delta is negative with
p <= 0.05 AND identity_loss is not flagged.
Report the count in every space and every origin/target pair. Zero is a
valid and reportable answer.

--- STAGE 3: DISAGREEMENT WITH DIRECTIONAL READOUTS ---
Rank all factors three ways: by delta_g, by cos_g, and by drop in the frozen
ruler score. Report Spearman correlation between each pair of rankings and
the overlap of the top 20. State how many of the top 20 by cosine are not
`toward_young` by delta.

--- STAGE 4: COMBINATIONS (assumption-laden; label as such) ---
Only if Stage 1 passed. State at the top of the section that this stage
assumes effects add linearly and that this data cannot test that assumption.
Greedy forward search: start empty, add at each step the factor whose sum
most reduces distance to the target, stop when no addition reduces it or when
the state enters the Y2 membership region, maximum size 5. Repeat from 50
random starts (seed 20260914) and report how often the same set is found.
For each reported set give delta, margin over the Stage 1A floor, identity
flag, and the delta with each member removed in turn.
Report no set as a candidate unless its delta clears the floor by at least
the Stage 1B MDA margin.

--- WHAT DOES NOT COUNT ---
- cosine or ruler-score drop alone as evidence of approach.
- any delta that does not clear the measured control floor.
- subtracting the floor instead of clearing it.
- choosing k, origin, target, or proliferation cut after seeing results.
- dropping or relabelling factors because they look odd; report them.
- pooling the two origins, two targets, or two spaces into one number.
- claims about the source paper's own conclusions; this task does not test them.

--- OUTPUT ---
FINDINGS_SOUTH.md: this file verbatim; Stage 0 inventory; both spaces and
their sizes; Stage 1 floor and MDA table; Stage 2 full table plus a top-20
table with factor, delta, margin over floor, q, identity flag, proliferation
flag; Stage 3 disagreement table; Stage 4 with its additivity warning; fired
keys; limitations, which must include: Hs27 is one neonatal line with no aged
cells, so every delta is a counterfactual assuming the effect transfers;
effects are measured in a different cell line and platform from the origin
and target; combinations assume additivity; per-factor mean rows may not
permit cell-level resampling.
PROGRESS_SOUTH.md: stop status, next action.
Record every failure verbatim. Do not substitute columns or repair rows.
```

## Stage 0 inventory (descriptive only; no statistic)

- File: `<repo>\results\survey\opened\southard\fibroblast_CRISPRa_mean_pop.h5ad`
- Size: 1,720,196,368 bytes. md5 `ba44c7813903bb5df900348d6b0d589a`.
- X: dense 10,916 rows x 4,914 genes, dtype `float64`. Layers: `adj_p`, `masked`, `p`. `raw` group: False. `uns` keys: none.
- obs columns (26): `protospacer`, `target_gene`, `target_expr`, `active`, `cluster`, `cluster_name`, `cluster_description`, `gene_driven`, `sequence_driven`, `sequence_suffix`, `bulk_expr`, `expressed`, `p_expressed`, `de_genes`, `strength`, `cell_count`, `target_gene_id`, `masked_active`, `masked_cluster`, `masked_cluster_name`, `masked_cluster_description`, `masked_cluster_size`, `masked_cluster_num_genes`, `stricter_masked_cluster`, `stricter_masked_cluster_description`, `guide_identity`.
- var columns (13): `feature_types`, `genome`, `mean`, `in_matrix`, `std`, `cv`, `fano`, `pairwise_p_cost`, `mean_adjusted_pairwise_p_cost`, `excess_cv`, `pairwise_chosen`, `gene_id`, `gene_name`.

### What one row is

**replicate (per-guide mean-population row).** `guide_identity` is unique on every row (10,916 values for 10,916 rows), while a targeted transcription factor has between 2 and 6 rows (the 17-row `off-target` group is the untargeted control arm, not a factor). So a row is not a per-factor mean (multiple rows per factor) and not a single cell (no barcode column; only a per-row `cell_count`: min 1, median 59, max 1369, total 762,352). Each row is the guide-level 'mean population' produced by the source's regression model, i.e. a per-guide replicate summary, not raw expression.

### Transcription factors, replicate rows, controls

- Distinct `target_gene` labels: 1,837. Of these, 1,836 are targeted transcription factors and 1 is an untargeted control label.
- Untargeted / non-targeting control rows: **17** (threshold for `no_controls` is 10). Control labels: `off-target`. Matched on `target_gene` or `guide_identity` against `(?i)(non[-_ ]?targeting|non[-_ ]?target|^off[-_ ]?target|scrambled?|scrambl|intergenic|safe[-_ ]?harbor|^nt[_-]|^ntc$|^ntc[_-]|^control$|^ctrl$|^ctrl[_-]|negative[-_ ]?control)`.
- Evidence those 17 rows are genuinely untargeted: 17/17 carry an empty `target_gene_id`, 17/17 a NaN `target_expr` (there is no intended target whose activation could be measured), 0/17 are `active`, and 0/17 are `gene_driven`. The label is also absent from the source's own 1,836-entry TF target list (`hs27_targets.txt`), which equals the 1,837 labels here minus this one.
- Those control rows are **not inert**: 3 of 17 are flagged `sequence_driven` and 10 have `de_genes` > 0. They also carry far more cells per row than the targeted guides (median 444, range 238-1369, against a median of 59 over all rows). A control floor built from them would therefore absorb real sequence-driven effects, which would raise the floor, not lower it.
- Hyphenated labels (18), audited explicitly because a hyphenated label passes a gene-symbol shape test: `BORCS8-MEF2B`, `CCDC169-SOHLH2`, `NKX1-1`, `NKX1-2`, `NKX2-1`, `NKX2-2`, `NKX2-3`, `NKX2-4`, `NKX2-5`, `NKX2-6`, `NKX2-8`, `NKX3-1`, `NKX3-2`, `NKX6-1`, `NKX6-2`, `NKX6-3`, `ZNF559-ZNF177`, `off-target`. 17 of them are standard HGNC symbols kept as transcription factors; the remainder are the control label above.

Replicate rows per transcription factor:

| n_replicate_rows | n_factors |
|---|---|
| 2 | 1 |
| 3 | 3 |
| 4 | 37 |
| 5 | 30 |
| 6 | 1765 |

1,836 factors have at least 2 replicate rows.

### Raw counts

- X holds raw integer counts: **False**. A raw layer with integer counts exists: **False**.

Audit of every matrix checked (same rule as `seng_run.locate_count_matrix`):

| where | dtype | n_seen | n_fractional | n_negative | n_nonfinite | min | max | raw_integer_counts |
|---|---|---|---|---|---|---|---|---|
| X | float64 | 53,641,224 | 53,640,987 | 30,370,037 | 0 | -9.156 | 184.8 | False |
| layers/adj_p | float64 | 53,641,224 | 53,619,366 | 0 | 0 | 0 | 1 | False |
| layers/masked | float64 | 53,641,224 | 53,638,740 | 30,369,689 | 2,247 | -9.156 | 184.8 | False |
| layers/p | float64 | 53,641,224 | 53,593,911 | 0 | 0 | 0 | 1 | False |

The other already-downloaded Southard files in the same folder were checked, so that what is and is not available locally is a verified fact rather than an assumption. None of them was substituted for the pinned table:

- `fibroblast_CRISPRa_aggr_total_guide_umis.h5` (58,166,462 bytes): `guide_umis/index_label0` (4799297,) int32; `guide_umis/index_label1` (4799297,) int16; `guide_umis/index_level0` (514841,) |S19; `guide_umis/index_level1` (10979,) |S35; `guide_umis/values` (4799297,) float32.
- `L4A_matrix.h5` (76,636,358 bytes): CellRanger matrix, 47,580 features x 18,197 barcodes (CRISPR Guide Capture 10,979, Gene Expression 36,601); `matrix/data` dtype `int32`, integer dtype True, range 1-13,784.

`fibroblast_CRISPRa_aggr_total_guide_umis.h5` cannot supply counts: it is a sparse COO table whose two index levels are cell barcodes and guides, i.e. guide-assignment UMIs, with no gene axis.

**`L4A_matrix.h5` does hold raw integer gene counts** (36,601 Gene Expression features x 18,197 cells, `int32`, minimum 1). Stating this plainly because it is the one fact that bears on what could be done next, and because it would be wrong to report that no counts for this screen exist locally. It was **not** used, for two reasons that are separate:

1. The pre-registration pins the analysed table to `fibroblast_CRISPRa_mean_pop.h5ad` under DATA, and applies the `no_raw_counts` test to that table. Swapping in a different file after the pinned one failed its own test would be choosing the data after seeing the result, which this task forbids.
2. It is a single lane. Its 18,197 cells are 2.4% of the 762,352 cells the pinned table summarises, spread over 10,979 guide features - on the order of 1.7 cells per guide. Even as a future path it would not by itself support a 1,836-factor ranking; the remaining lanes would have to be obtained.

### Frozen-ruler genes present

- 4,907 of 23,485 frozen-ruler genes (expected 23,485) are present among the file's 4,914 genes (Ensembl id first, symbol fallback; `seng_run.ruler_column_index`). The ruler is not refit.
- Identity genes present: COL1A1, COL1A2, FN1, LUM, PDGFRA, PDGFRB, PRRX1, SERPINH1, VIM. Missing: POSTN.

### Resampling unit

- **replicates.** Cells are not in this file (only a per-row cell count), so cell-level resampling is impossible; 1,836 factors carry >= 2 replicate rows, so replicate-level resampling would be the available unit.

### Recorded failure of this task's own first Stage 0 pass

The first pass of SOUTH reported `0` untargeted control rows and `1,837` transcription factors, and fired `no_controls` on that basis. That was wrong. Verbatim, it reported:

> - Distinct `target_gene` labels: 1,837, all targeted transcription factors (1,837 after removing control-matched labels, of which there are none).
> - Untargeted / non-targeting control rows: **0**.

The error was in the control-label pattern, not in the data: it tested for `non-targeting`, `scrambled`, `intergenic`, `safe-harbor`, `NTC`, `ctrl` and similar, but not for `off-target`, which is the name this library gives its untargeted arm. The 17 `off-target` rows were therefore counted as a 1,837th transcription factor — the one factor the first pass reported as having 17 replicate rows, where every real factor has 2 to 6.

Corrected here: 1,836 transcription factors and 17 control rows, so `no_controls` does **not** fire. The `off-target` rows were not dropped or relabelled to make anything look better; they are reported with the evidence above.

This correction does not change the task's outcome. `no_raw_counts` fires independently and halts SOUTH at Stage 0 either way. It changes which keys are reported as fired, and it removes a false claim that this dataset carries no controls.

## Spaces and their sizes

- S1 (gene space): the 23,485 frozen-ruler overlap genes, of which 4,907 are present in this file; missing genes are z=0 by construction. S1 was **not built**: its pipeline (sum -> TMM log2-CPM with prior.count=2 -> frozen GTEx mu/sd -> z) requires raw counts, and this file has none.
- S2 (low-dim space): PCA fit only on the GTEx fibroblast donor matrix at k in (20, 50, 100), applied unchanged elsewhere. **Not fitted**: the task halted at Stage 0 before any space was built.

## Targets and origins

Not built. Y1 (GTEx fibroblast donors 20-39), Y2 (the same donors as a Mahalanobis region in S2, with a held-out 95th-percentile radius), O1 (GTEx fibroblast donors 60-79), O2 (the Lu et al. GM00731 aged day-0 pseudobulk as built in SAME), and every origin-to-target distance are downstream of the Stage 0 STOP and were not computed.

## Stage 1 detection limit (control floor, dose-response, MDA)

Not run. Stage 0 fired `no_raw_counts`, which halts the task before any statistic. The row count needed for Stage 1A exists (17 control rows, enough to split into two disjoint halves 200 times), but the quantity being split is z, and z is defined here as sum -> TMM log2-CPM with prior.count=2 -> frozen GTEx mu/sd. That pipeline starts from raw counts, which this file does not contain. No floor, no dose-response, no MDA, and no `no_detection_power` determination were computed.

## Stage 2 single-factor ranking

Not run. Stage 1 did not pass (it did not run), and the per-factor effect vector d_g = z(g rows) - z(control rows) cannot be formed, because z cannot be formed without raw counts. The count of `toward_young` factors is therefore undefined, not zero.

## Stage 3 disagreement with directional readouts

Not run. It ranks factors by delta, cosine, and ruler-score drop, all of which require Stage 2 outputs.

## Stage 4 combinations

**Additivity warning (stated as required):** this stage would assume that perturbation effects add linearly, and this data cannot test that assumption. Not run, because Stage 1 did not pass.

## Fired STOP keys

Fired, in pre-registration order: `no_raw_counts`.
Checked and not fired: `no_resampling_unit`, `no_controls`.

Failures recorded verbatim:

- **no_raw_counts:** X is not raw integer counts and no raw layer exists. audits=[{'dtype': 'float64', 'integer_dtype': False, 'n_seen': 53641224, 'n_fractional': 53640987, 'n_negative': 30370037, 'n_nonfinite': 0, 'min': -9.156079292297363, 'max': 184.80467224121094, 'raw_integer_counts': False, 'where': 'X'}, {'dtype': 'float64', 'integer_dtype': False, 'n_seen': 53641224, 'n_fractional': 53619366, 'n_negative': 0, 'n_nonfinite': 0, 'min': 0.0, 'max': 0.9999999985175676, 'raw_integer_counts': False, 'where': 'layers/adj_p'}, {'dtype': 'float64', 'integer_dtype': False, 'n_seen': 53641224, 'n_fractional': 53638740, 'n_negative': 30369689, 'n_nonfinite': 2247, 'min': -9.156079292297363, 'max': 184.80467224121094, 'raw_integer_counts': False, 'where': 'layers/masked'}, {'dtype': 'float64', 'integer_dtype': False, 'n_seen': 53641224, 'n_fractional': 53593911, 'n_negative': 0, 'n_nonfinite': 0, 'min': 0.0, 'max': 0.9999999985175676, 'raw_integer_counts': False, 'where': 'layers/p'}] layers=['adj_p', 'masked', 'p'] raw_group=False.

## Limitations

- Hs27 is one neonatal fibroblast line with no aged cells, so every delta this task could have produced is a counterfactual that assumes the perturbation effect transfers to aged donors.
- Effects here are measured in a different cell line and platform (Hs27 CRISPRa screen) from the origin and target (GTEx / Lu et al. donor fibroblasts), so any comparison crosses a batch and platform boundary.
- Combinations (Stage 4) assume additivity, which this data cannot test.
- The rows are per-guide 'mean population' summaries from a regression model, not raw counts and not single cells, so cell-level resampling is not possible and the frozen z-pipeline cannot be applied at all.
- Decisive blocker this run: the matrix is not raw integer counts and no raw layer exists (`no_raw_counts`). That halts the ranking on its own; no delta, floor, or MDA exists to report, and none was invented.
- Even had counts existed, only 4,907 of the 23,485 frozen-ruler genes (21%) are present in this file, so 79% of the S1 coordinates would have been set to z=0 by the missing-gene rule. Every S1 distance would then be dominated by a fixed, perturbation-independent offset. This is recorded as a limitation, not used as a gate.
- The 17 control rows are untargeted but not inert (3 are flagged `sequence_driven` and 10 have `de_genes` > 0), and they carry about 8x more cells per row than the median row. A floor built from them would be both noisier and less exchangeable with a real factor's rows than the pre-registration's row-count matching alone assumes.

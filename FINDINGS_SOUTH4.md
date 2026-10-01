> **SUPERSEDED AND WITHDRAWN by SOUTH5 (FINDINGS_SOUTH5.md).** The `guide_assignment` stop was a parser limitation, not a property of the data; the obs columns are named explicitly in SOUTH5. `results/south4/` has been deleted as untrustworthy. Nothing below is changed.

# FINDINGS_SOUTH4 — rerun the Southard ranking on the full measured data

**Status:** STOP `guide_assignment`.

Guide assignment is not stored in a recognised target/guide column. obs columns=['dataset', 'UMI_count', 'guide_identity', 'guide_umi_count', 'thresholded_features', 'thresholded_guide_umi', 'num_cells', 'guide_target', 'gem_group', 'mt_frac', 'keep', 'single_cell', 'control', 'pairwise_coef', 'pairwise_UMI_count', 'pairwise_singlet_equivalent_UMI_count', 'final_pairwise_UMI_count', 'protospacer']. Not inventing a parser. target_col=None guide_col=guide_identity.

The fenced pre-registration was written to this file and to `results/south4/PREREG.flag` before any statistic. The fence is the flag file, not a paraphrase.

## Pre-registration (verbatim, written before any statistic)

```
TASK: SOUTH4 — rerun the Southard ranking on the full measured data.

Nature of this task: computational reanalysis of a public human cell dataset.
No lab work, no protocol, no new data generation. Output is a ranked table.

SOUTH2 and SOUTH3 used results/survey/opened/southard/fibroblast_CRISPRa_mean_pop.h5ad,
which holds the authors' 4,914-gene analysis subset as regression output. That
covered 12% of the frozen ruler's absolute weight and its units were unstated.
SURVEY_TF.md recorded, before any of those runs, that the full measured data
covers 63% of ruler weight across 36,601 genes as raw counts. SOUTH4 uses the
full data. This is not a post-hoc dataset change: the target file was named in
the survey before SOUTH2 ran.

Write FINDINGS_SOUTH4.md, PROGRESS_SOUTH4.md, results/south4/*, src/south4_run.py.
Add one line at the top of FINDINGS_SOUTH3.md marking it superseded; change
nothing else in it. Do not modify any other FINDINGS_*, PROGRESS_*, or src/ file.
Import from src/toward_run.py, src/same_run.py, src/seng_run.py,
src/lowdim_common.py; do not change them. Do not refit the frozen ruler.
Seeds: 20260914, 20260918. n_perm=20000, n_boot=200.

WRITE THE PRE-REGISTRATION FIRST: copy this file verbatim into
FINDINGS_SOUTH4.md and write it to results/south4/PREREG.flag before any
statistic. Nothing is re-tuned after numbers exist.

--- DATA ---
Southard et al. 2025 Hs27 CRISPRa Perturb-seq, Zenodo 10.5281/zenodo.15200179.
Use the per-cell singlets file (~8.05 GB). Download it to data/raw/southard/.
Record size and md5. If the download cannot complete, STOP with key
`download_failed` and report how far it got. Do not substitute the
mean-population file or the single L4A lane.

--- STAGE 0: INVENTORY (descriptive only) ---
Report: matrix shape; whether X is raw integer counts (STOP `no_raw_counts`
if not and no raw layer exists); how guide assignment is stored; number of
distinct target genes; guides per target (min/median/max); number of
non-targeting guides and cells; cells per target gene.
Ruler coverage three ways (symbol, version-stripped Ensembl, best after
src/md3_idtype.py), reporting for each both the gene count and the fraction
of the ruler's total absolute weight. Use the highest-weight mapping and name
it. Record `coverage_low` if best weight fraction < 0.25.
STOP `too_few_factors` if fewer than 500 targets have 2+ guides.
STOP `no_controls` if fewer than 20 non-targeting guides.

--- STAGE 1: VECTORS AND SPACES ---
Pseudobulk per target gene (all its cells) and per non-targeting guide.
Process exactly as TOWARD/SAME: sum counts, TMM, log2 CPM prior.count=2
across this task's own panel, frozen GTEx mu/sd z-score, missing genes z=0.
Effect vector d_g = z(target g) - z(all non-targeting cells pooled).
One gene variant only: all matched genes. No significance masking; SOUTH3
showed masking makes guides orthogonal by construction.
Spaces: S1 gene space; S2 low-dim PCA fit ONLY on the GTEx fibroblast donor
matrix, k in (20, 50, 100), applied unchanged. Report every k.

--- ORIGINS AND TARGETS ---
PRIMARY, declared now and reported first: O1 = GTEx donors aged 60-79,
Y1 = GTEx donors aged 20-39, in S1 and each S2 k.
SECONDARY, reported but never primary: O2 = Lu aged donor GM00731 day-0
pseudobulk; Y2 = young donors as a region (mean, covariance, Mahalanobis,
radius at the 95th percentile of held-out young donors), S2 only, skipped
at any k where the covariance is singular.
Reason for the primary choice, stated now: SOUTH3 found ||O2-Y1|| is 7-15x
||O1-Y1||, so O2 makes any favourable projection lower delta, and the aged
GTEx origin already sits inside the Y2 region, so Y2 membership does not
discriminate at O1. This choice is made before any SOUTH4 number exists.
Report each origin-to-target distance before scoring anything.

--- STAGE 2: DETECTION LIMIT (before any ranking) ---
A. Primary floor: permute each factor's own vector across genes 20,000 times;
   delta = ||(origin + d_perm) - target|| - ||origin - target||; the floor is
   the 5th percentile. Report per-factor and pooled-at-median-norm floors.
B. Control floor: same delta using each non-targeting guide as a pseudo-effect,
   and using random splits of the non-targeting cells. Report both. With 78
   guides this is now a real floor; report it alongside A and say which is
   stricter.
C. Dose-response: add the true young-minus-old direction to the origin at
   f = 0, 0.05, 0.1, 0.25, 0.5, 1.0, scaled to the median factor norm.
   MDA = smallest f clearing the primary floor.
   STOP `no_detection_power` if no f <= 1.0 clears it anywhere.
Never subtract a floor. It is a bar to clear.

--- STAGE 3: RANKING ---
Per factor, per setting: cos, frac, delta, and for Y2 the Mahalanobis change
and membership. Guide agreement = median pairwise cosine between that factor's
own guide-level vectors, computed on unmasked vectors; flag `guides_disagree`
if <= 0. Identity check on COL1A1, COL1A2, FN1, LUM, PDGFRA, PDGFRB, POSTN,
PRRX1, SERPINH1, VIM (report which are present); flag `identity_loss` if the
drop exceeds 0.5. Proliferation-residualized cos and delta as in seng_run.py.
Frozen ruler score change as context only.
Nulls: Stage 2A permutations, one-sided p, BH across factors within setting.
Bootstrap over cells; flag any interval not containing its point estimate as
INVALID and never gate on it.
Label `toward_young` only if delta < 0 AND below the factor's own Stage 2A
floor AND q <= 0.05 AND not identity_loss AND not guides_disagree.
Report counts per setting; never pool. Zero is a valid answer.

--- STAGE 4: DOES DELTA ADD ANYTHING? ---
Repeat SOUTH3's check at this coverage. Report ‖d‖/‖target-origin‖ per factor;
Spearman between delta and cosine rankings overall and by quartile of that
ratio; the residual of delta ≈ -‖d‖cos(theta); and what fraction of survivors
a plain top-n-by-cosine recovers. State plainly whether delta separates from
cosine at this coverage, and in which ratio quartile if anywhere.

--- STAGE 5: COMBINATIONS ---
State at the top of the section that this assumes effects add linearly and
that no public fibroblast screen can test that assumption.
Greedy forward search on summed vectors: add the factor that most reduces
distance; stop when nothing reduces it, or on entering the Y2 region, or at
size 5. Repeat from 50 random starts (seed 20260914); report how often the
same set recurs. Per set: delta, margin over floor, ‖d‖/‖gap‖ ratio, identity
flag, and delta with each member removed.
Report no set unless its delta clears the floor by at least the Stage 2B MDA
margin. Report whether combinations reach a ratio where Stage 4 showed delta
and cosine separating.

--- WHAT DOES NOT COUNT ---
- cosine or ruler-score drop alone as evidence of approach.
- any delta not clearing the measured primary floor.
- subtracting a floor instead of clearing it.
- choosing k, mapping, origin or target after seeing results.
- substituting a smaller file for the pinned one.
- dropping or relabelling factors because they look odd; report them.
- pooling settings into one number.
- claims about the source paper's own conclusions.

--- OUTPUT ---
FINDINGS_SOUTH4.md: prereg verbatim; Stage 0 inventory with all three coverage
numbers; Stage 2 floors and MDA; Stage 3 counts plus a top-20 candidate table
(factor, delta, margin over floor, q, guide agreement, identity flag,
proliferation flag) written to results/south4/top20.csv; Stage 4 verdict;
Stage 5 combinations; fired keys; limitations, which must include: Hs27 is one
neonatal line with no aged cells, so every delta is a counterfactual assuming
the effect transfers; combinations assume additivity, untestable here; unmeasured
genes assumed unmoved.
PROGRESS_SOUTH4.md: stop status, next action.
Record every failure verbatim. Do not substitute columns or repair rows.
```

## Limitations

- Hs27 is one neonatal line with no aged cells, so every delta is a counterfactual assuming the effect transfers.
- Combinations assume additivity, untestable here. No public fibroblast screen can test that assumption.
- Unmeasured genes are assumed unmoved (z = 0) and are not permuted into.
- The primary origin and target are GTEx fibroblast donors, not cells from this line. O2 is one aged donor's day-0 pseudobulk from a different study.
- Non-targeting guides are the control the experiment measured. They are not a proof that the guides are inert.
- The frozen ruler was not refit. Its score is context only and is not evidence of approach.

## Failures (verbatim)

- **guide_assignment:** Guide assignment is not stored in a recognised target/guide column. obs columns=['dataset', 'UMI_count', 'guide_identity', 'guide_umi_count', 'thresholded_features', 'thresholded_guide_umi', 'num_cells', 'guide_target', 'gem_group', 'mt_frac', 'keep', 'single_cell', 'control', 'pairwise_coef', 'pairwise_UMI_count', 'pairwise_singlet_equivalent_UMI_count', 'final_pairwise_UMI_count', 'protospacer']. Not inventing a parser. target_col=None guide_col=guide_identity.

## Flags

None recorded.

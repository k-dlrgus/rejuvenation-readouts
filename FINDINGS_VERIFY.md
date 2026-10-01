# FINDINGS_VERIFY

**Status:** Completed. No stop key fired.

The fenced pre-registration was written to this file and to `results/verify/PREREG.flag` before the draft was opened. The fence is the flag file, not a paraphrase.

## Draft identity

Audited file: `paper/rejuvenation_readouts_preprint_3.docx`

- size: 1,102,439 bytes
- md5: `6585db8d02191542dc867abe015ddcf5`
- modification time: 2026-09-28 03:39:46

Other draft versions in the repository, not audited:

- `paper/PREPRINT_DRAFT.md` — 25,348 bytes, modification time 2026-09-28 03:40:52
- `paper_package/paper/PREPRINT_DRAFT.md` — 24,562 bytes, modification time 2026-09-23 18:42:52

No other `rejuvenation_readouts_preprint_*.docx` is in the repository.

## Stage 1

Claim rows in `results/verify/claims.csv`: 409.

## Stage 2

Verdict counts: match 274; rounding 5; mismatch 0; untraceable 36; not_a_claim 94.

### Rounding

These agree with the stored value except in the last printed digit.

- `T2R3` number `+0.413`. Stored `0.41248891146806693` in `results/paper_figs/cxg_true_donor_results.csv` field `true_donor ruler_minus_MD delta_ci_hi`.

> Fibroblast atlas, 65 donors | 0.433 | 0.26 | 0.06 | +0.171 [−0.118, +0.413] | +0.372 [+0.094, +0.595]

- `P39` number `+0.413`. Stored `0.41248891146806693` in `results/paper_figs/cxg_true_donor_results.csv` field `MD delta_ci_hi`.

> The evidence is therefore uneven, and we state it that way. In GSE226189 the ruler beats both curated instruments clearly, and neither MD nor the age-up and age-down lists order donors better than chance there (MD permutation p = 0.26). In the fibroblast atlas the ruler still beats the age-up and age-down lists after the unit correction (+0.372, 95% CI +0.094 to +0.595), but its advantage over MD (+0.171, CI −0.118 to +0.413) is no longer distinguishable from zero. The original analysis had reported that interval as excluding zero; that result depended on counting repeated samples from the same person as independent, and we withdraw it. The claim that survives is that the ruler tracks donor age better than MD in one clean cohort and at least as well in a second. MD was also not applied exactly as published: the gene list is the authors' own, but the score is a Python reimplementation of the Seurat AddModuleScore algorithm on log-normalised pseudobulks, because the published environment could not be run here. A third pre-registered cohort, GSE113957, could not be scored at all, because the deposit provides FPKM values rather than counts and our normalisation requires counts; it is reported as not scored rather than as a negative result.

- `P8` number `two thirds`. Stored `0.674440748216413` in `results/same/stats.csv` field `forward frac`.

> Applying the distance-to-target framework to OSKM partial reprogramming of fibroblasts (Lu et al.), the aged donor's partially reprogrammed cells end farther from young fibroblasts than the cells they came from, whether young is defined by GTEx donors (+9.4 z units) or by a 22-year-old donor assayed in the same experiment (+131), even though against the latter target they travel two thirds of the way along the aged-to-young axis. A dose-response control on the same statistic detects a known approach, so the null result is not simply insensitivity, and a reverse test rules out mutual convergence of the two donors.

- `P51` number `two thirds`. Stored `0.674440748216413` in `results/same/stats.csv` field `forward frac`.

> Figure 5. Same-platform test against a young target measured in the same experiment. (A) The forward test drawn in the plane spanned by the origin-to-target axis and the orthogonal component of each displacement, in units of the axis length. Points inside the dashed circle are closer to the target than the origin was. The aged donor's partially reprogrammed cells project two thirds of the way along the axis and still land far outside the circle. Grey points are the positive-control mixtures. (B) The reverse test, with donors swapped. (C) Progress and (D) change in distance for the mixture control; open orange circles mark the zero-mixture noise floor.

- `P53` number `two thirds`. Stored `0.674440748216413` in `results/same/stats.csv` field `forward frac`.

> In the dataset where the distinction can be tested, what is scored as rejuvenation is well explained by departure from the aged state, and is not accompanied by a measurable approach to a young cell of the same type. The aged donor's partially reprogrammed cells project two thirds of the way along the aged-to-young axis while ending twice as far from the young target as they started, and the same sign holds against a GTEx-defined young centroid. A dose-response control shows the statistic would have seen an approach had one occurred, and a reverse test shows the movement is not mutual convergence of two donors.

### Mismatches

None.

### Untraceable

#### 1. `P7` — `95%`

Searched: searched results/fibro/manifest.json, results/md3/t2_delta_rho.csv and results/toward/summary.json for a stored percentile 95, 2.5 or 97.5; the draft states the interval level, and those files store n_boot=200 without the percentile

> We built a transcriptomic age ruler in cultured human dermal fibroblasts (652 GTEx donors) by supervised ridge regression on donor age, and froze it. It orders donors by age in cohorts that differ in laboratory, collection site and platform (Spearman ρ 0.43 to 0.56 across four transfers), but it does not predict age in calibrated years. In an independent bulk cohort of 82 donors, the frozen ruler tracks donor age better than the curated mesenchymal-drift (MD) signature used as an age readout by Lu et al. (paired bootstrap Δρ +0.47, 95% CI +0.24 to +0.70). In a second, single-cell cohort the difference points the same way, but once its samples are grouped into the 65 people they came from, the interval includes zero (+0.17, −0.12 to +0.41).

#### 2. `P27` — `95%`

Searched: searched results/fibro/manifest.json, results/md3/t2_delta_rho.csv and results/toward/summary.json for a stored percentile 95, 2.5 or 97.5; the draft states the interval level, and those files store n_boot=200 without the percentile

> One of these transfers is weaker than it looks. The dermal fibroblast atlas is an integration of published studies rather than a single cohort, and the pre-registered analysis treated each of its 93 donor identifiers as a person. They are not. Tracing each identifier to its source study shows that the 93 pseudobulks come from 65 people: one study contributes 23 sequencing accessions from five donors, whose recorded age of 42 is a placeholder for an author-reported range of 25 to 60 years, and several other studies contribute repeated samples from the same donor. With the correct unit the transfer holds (ρ = +0.433, 95% CI +0.174 to +0.610, permutation p = 0.005, against +0.455 as run), but the effective sample is smaller and some ages are imprecise. We report both values and interpret the corrected one. The bulk cohort GSE226189, in which 82 samples carry 82 distinct donor identifiers, is the cleaner of the two external tests.

#### 3. `P29` — `95%`

Searched: searched results/fibro/manifest.json, results/md3/t2_delta_rho.csv and results/toward/summary.json for a stored percentile 95, 2.5 or 97.5; the draft states the interval level, and those files store n_boot=200 without the percentile

> Figure 1. Transfer of the frozen fibroblast ruler. (A) Spearman ρ with donor age for every transfer in Table 1, with 95% percentile bootstrap intervals (B = 200). Filled circles are pre-registered tests; the open circle is the same fibroblast atlas test with the analysis unit corrected from pseudobulks to people. (B) GSE226189 bulk cohort, 82 donors. (C) Dermal fibroblast atlas, 93 pseudobulks from 65 people; the 23 pseudobulks from five donors of one constituent study, whose age is recorded as 42 but reported by the source authors only as 25 to 60 years, are drawn as open circles.

#### 4. `P30` — `95%`

Searched: searched results/fibro/manifest.json, results/md3/t2_delta_rho.csv and results/toward/summary.json for a stored percentile 95, 2.5 or 97.5; the draft states the interval level, and those files store n_boot=200 without the percentile

> One property of the fibroblast ruler matters for everything that follows. Its score is anticorrelated with distance from its own training distribution: across 22 donor × cluster × timepoint pseudobulks in the reprogramming dataset, Spearman ρ between nearest-neighbour distance to the GTEx training donors and the frozen score is −0.625 (95% CI −0.834 to −0.320; Figure 2), and the same sign holds within each donor separately (−0.78 and −0.81). Cells far outside the training manifold are scored as young. Because distance and reprogramming state move together here, this correlation cannot by itself separate an extrapolation artefact from genuine change; pluripotent cells, for example, are reset to near-zero age by DNA-methylation clocks. What it does show is that a low ruler score on these cells cannot be taken at face value, which is a plausible mechanism by which a one-dimensional score mistakes departure for rejuvenation, and it applies to our ruler exactly as it applies to the readouts we examine below. This analysis was pre-registered as report-only, and it had anticipated the opposite sign, a score that rises with distance; the observation is that it falls.

#### 5. `T2R1` — `95%`

Searched: searched results/fibro/manifest.json, results/md3/t2_delta_rho.csv and results/toward/summary.json for a stored percentile 95, 2.5 or 97.5; the draft states the interval level, and those files store n_boot=200 without the percentile

> Cohort | Ruler ρ | MD ρ | Age-up − age-down ρ | Δρ vs MD [95% CI] | Δρ vs age-up/down [95% CI]

#### 6. `T2R1` — `95%`

Searched: searched results/fibro/manifest.json, results/md3/t2_delta_rho.csv and results/toward/summary.json for a stored percentile 95, 2.5 or 97.5; the draft states the interval level, and those files store n_boot=200 without the percentile

> Cohort | Ruler ρ | MD ρ | Age-up − age-down ρ | Δρ vs MD [95% CI] | Δρ vs age-up/down [95% CI]

#### 7. `P38` — `95%`

Searched: searched results/fibro/manifest.json, results/md3/t2_delta_rho.csv and results/toward/summary.json for a stored percentile 95, 2.5 or 97.5; the draft states the interval level, and those files store n_boot=200 without the percentile

> Figure 3. The frozen ruler against two curated instruments. (A) Spearman ρ with donor age for each instrument in GSE226189, with 95% percentile bootstrap intervals. (B) Paired Δρ (ruler minus the other instrument). Grey open points are the fibroblast atlas as originally run, with pseudobulks treated as donors; black points use real donors. Against MD, the atlas interval includes zero once the unit is corrected. A third pre-registered cohort, GSE113957, is not shown because it could not be scored.

#### 8. `P39` — `95%`

Searched: searched results/fibro/manifest.json, results/md3/t2_delta_rho.csv and results/toward/summary.json for a stored percentile 95, 2.5 or 97.5; the draft states the interval level, and those files store n_boot=200 without the percentile

> The evidence is therefore uneven, and we state it that way. In GSE226189 the ruler beats both curated instruments clearly, and neither MD nor the age-up and age-down lists order donors better than chance there (MD permutation p = 0.26). In the fibroblast atlas the ruler still beats the age-up and age-down lists after the unit correction (+0.372, 95% CI +0.094 to +0.595), but its advantage over MD (+0.171, CI −0.118 to +0.413) is no longer distinguishable from zero. The original analysis had reported that interval as excluding zero; that result depended on counting repeated samples from the same person as independent, and we withdraw it. The claim that survives is that the ruler tracks donor age better than MD in one clean cohort and at least as well in a second. MD was also not applied exactly as published: the gene list is the authors' own, but the score is a Python reimplementation of the Seurat AddModuleScore algorithm on log-normalised pseudobulks, because the published environment could not be run here. A third pre-registered cohort, GSE113957, could not be scored at all, because the deposit provides FPKM values rather than counts and our normalisation requires counts; it is reported as not scored rather than as a negative result.

#### 9. `P47` — `95%`

Searched: searched results/fibro/manifest.json, results/md3/t2_delta_rho.csv and results/toward/summary.json for a stored percentile 95, 2.5 or 97.5; the draft states the interval level, and those files store n_boot=200 without the percentile

> Young defined within the experiment. To remove the bulk-versus-single-cell platform gap, we repeated the test with the 22-year-old donor's untreated day-0 cells as the target, splitting each donor's day-0 fibroblasts in half so that the anchors and the positive control use disjoint cells. The result is the clearest illustration of the paper's argument. The aged donor's partially reprogrammed cells travel 67% of the way along the aged-to-young donor axis (progress 0.674, 95% CI 0.667 to 0.680), which on a projection-based readout would look like substantial rejuvenation, and yet their distance to the young target increases by 131 z units (Figure 5A). The displacement is long and largely orthogonal: its cosine with the axis is 0.323, so the cells move far enough sideways that projecting a large distance along the axis is compatible with ending much farther away.

#### 10. `P49` — `95%`

Searched: searched results/fibro/manifest.json, results/md3/t2_delta_rho.csv and results/toward/summary.json for a stored percentile 95, 2.5 or 97.5; the draft states the interval level, and those files store n_boot=200 without the percentile

> Convergence versus rejuvenation. A treatment that merely erased donor-specific expression would move both donors toward each other and mimic rejuvenation. Reversing the test, the young donor's partially reprogrammed cells do not move toward the aged donor at all (progress −0.057, 95% CI −0.063 to −0.054), giving a clear asymmetry (+0.732, CI +0.724 to +0.740). The forward movement is therefore not symmetric mutual convergence. But it is also not obviously age: the aged-to-young donor axis is itself nearly orthogonal to the GTEx donor-age axis (cosine 0.039), so with one donor at each age we cannot separate age from everything else that differs between two cell lines, including proliferative state and culture history.

#### 11. `P74` — `0.01`

Searched: searched results/fibro/freeze.json and results/fibro/manifest.json for the alpha grid; only alpha=0.01 is stored, not 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000 or 10000

> The ruler is a ridge regression of donor chronological age on z-scored gene expression, fit through the n × n Gram matrix, and frozen with its coefficients and its training mean and standard deviation before any evaluation. The regularisation parameter was selected over the grid 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000 by minimising a leave-one-out residual sum of squares computed from the hat-matrix diagonal, not by grouped cross-validation. The frozen fibroblast ruler stores alpha = 0.01, the smallest value on the grid, so a smaller optimum cannot be excluded. The frozen fibroblast vector is the ridge coefficient vector, unit-normalised and sign-aligned so that the score increases with age. Transfers within GTEx use the full prediction (training mean age plus the ridge coefficients applied to the z-scored profile), which is in years. The external cohorts and all reprogramming analyses use the frozen unit-length score, which orders samples identically but has no intercept or scale and is not in years. A PLS-1 model was fit in parallel for comparison and is neither the frozen ruler nor any gate. Within-site fibroblast cross-validation, grouped by sequencing batch with up to five folds, reports ρ but does not choose alpha. The ruler was never refit on any evaluation cohort.

#### 12. `P74` — `0.03`

Searched: alpha grid value is not stored; results/fibro/freeze.json has alpha 0.01 only

> The ruler is a ridge regression of donor chronological age on z-scored gene expression, fit through the n × n Gram matrix, and frozen with its coefficients and its training mean and standard deviation before any evaluation. The regularisation parameter was selected over the grid 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000 by minimising a leave-one-out residual sum of squares computed from the hat-matrix diagonal, not by grouped cross-validation. The frozen fibroblast ruler stores alpha = 0.01, the smallest value on the grid, so a smaller optimum cannot be excluded. The frozen fibroblast vector is the ridge coefficient vector, unit-normalised and sign-aligned so that the score increases with age. Transfers within GTEx use the full prediction (training mean age plus the ridge coefficients applied to the z-scored profile), which is in years. The external cohorts and all reprogramming analyses use the frozen unit-length score, which orders samples identically but has no intercept or scale and is not in years. A PLS-1 model was fit in parallel for comparison and is neither the frozen ruler nor any gate. Within-site fibroblast cross-validation, grouped by sequencing batch with up to five folds, reports ρ but does not choose alpha. The ruler was never refit on any evaluation cohort.

#### 13. `P74` — `0.1`

Searched: alpha grid value is not stored; results/fibro/freeze.json has alpha 0.01 only

> The ruler is a ridge regression of donor chronological age on z-scored gene expression, fit through the n × n Gram matrix, and frozen with its coefficients and its training mean and standard deviation before any evaluation. The regularisation parameter was selected over the grid 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000 by minimising a leave-one-out residual sum of squares computed from the hat-matrix diagonal, not by grouped cross-validation. The frozen fibroblast ruler stores alpha = 0.01, the smallest value on the grid, so a smaller optimum cannot be excluded. The frozen fibroblast vector is the ridge coefficient vector, unit-normalised and sign-aligned so that the score increases with age. Transfers within GTEx use the full prediction (training mean age plus the ridge coefficients applied to the z-scored profile), which is in years. The external cohorts and all reprogramming analyses use the frozen unit-length score, which orders samples identically but has no intercept or scale and is not in years. A PLS-1 model was fit in parallel for comparison and is neither the frozen ruler nor any gate. Within-site fibroblast cross-validation, grouped by sequencing batch with up to five folds, reports ρ but does not choose alpha. The ruler was never refit on any evaluation cohort.

#### 14. `P74` — `0.3`

Searched: alpha grid value is not stored; results/fibro/freeze.json has alpha 0.01 only

> The ruler is a ridge regression of donor chronological age on z-scored gene expression, fit through the n × n Gram matrix, and frozen with its coefficients and its training mean and standard deviation before any evaluation. The regularisation parameter was selected over the grid 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000 by minimising a leave-one-out residual sum of squares computed from the hat-matrix diagonal, not by grouped cross-validation. The frozen fibroblast ruler stores alpha = 0.01, the smallest value on the grid, so a smaller optimum cannot be excluded. The frozen fibroblast vector is the ridge coefficient vector, unit-normalised and sign-aligned so that the score increases with age. Transfers within GTEx use the full prediction (training mean age plus the ridge coefficients applied to the z-scored profile), which is in years. The external cohorts and all reprogramming analyses use the frozen unit-length score, which orders samples identically but has no intercept or scale and is not in years. A PLS-1 model was fit in parallel for comparison and is neither the frozen ruler nor any gate. Within-site fibroblast cross-validation, grouped by sequencing batch with up to five folds, reports ρ but does not choose alpha. The ruler was never refit on any evaluation cohort.

#### 15. `P74` — `1`

Searched: alpha grid value is not stored; results/fibro/freeze.json has alpha 0.01 only

> The ruler is a ridge regression of donor chronological age on z-scored gene expression, fit through the n × n Gram matrix, and frozen with its coefficients and its training mean and standard deviation before any evaluation. The regularisation parameter was selected over the grid 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000 by minimising a leave-one-out residual sum of squares computed from the hat-matrix diagonal, not by grouped cross-validation. The frozen fibroblast ruler stores alpha = 0.01, the smallest value on the grid, so a smaller optimum cannot be excluded. The frozen fibroblast vector is the ridge coefficient vector, unit-normalised and sign-aligned so that the score increases with age. Transfers within GTEx use the full prediction (training mean age plus the ridge coefficients applied to the z-scored profile), which is in years. The external cohorts and all reprogramming analyses use the frozen unit-length score, which orders samples identically but has no intercept or scale and is not in years. A PLS-1 model was fit in parallel for comparison and is neither the frozen ruler nor any gate. Within-site fibroblast cross-validation, grouped by sequencing batch with up to five folds, reports ρ but does not choose alpha. The ruler was never refit on any evaluation cohort.

#### 16. `P74` — `3`

Searched: alpha grid value is not stored; results/fibro/freeze.json has alpha 0.01 only

> The ruler is a ridge regression of donor chronological age on z-scored gene expression, fit through the n × n Gram matrix, and frozen with its coefficients and its training mean and standard deviation before any evaluation. The regularisation parameter was selected over the grid 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000 by minimising a leave-one-out residual sum of squares computed from the hat-matrix diagonal, not by grouped cross-validation. The frozen fibroblast ruler stores alpha = 0.01, the smallest value on the grid, so a smaller optimum cannot be excluded. The frozen fibroblast vector is the ridge coefficient vector, unit-normalised and sign-aligned so that the score increases with age. Transfers within GTEx use the full prediction (training mean age plus the ridge coefficients applied to the z-scored profile), which is in years. The external cohorts and all reprogramming analyses use the frozen unit-length score, which orders samples identically but has no intercept or scale and is not in years. A PLS-1 model was fit in parallel for comparison and is neither the frozen ruler nor any gate. Within-site fibroblast cross-validation, grouped by sequencing batch with up to five folds, reports ρ but does not choose alpha. The ruler was never refit on any evaluation cohort.

#### 17. `P74` — `10`

Searched: alpha grid value is not stored; results/fibro/freeze.json has alpha 0.01 only

> The ruler is a ridge regression of donor chronological age on z-scored gene expression, fit through the n × n Gram matrix, and frozen with its coefficients and its training mean and standard deviation before any evaluation. The regularisation parameter was selected over the grid 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000 by minimising a leave-one-out residual sum of squares computed from the hat-matrix diagonal, not by grouped cross-validation. The frozen fibroblast ruler stores alpha = 0.01, the smallest value on the grid, so a smaller optimum cannot be excluded. The frozen fibroblast vector is the ridge coefficient vector, unit-normalised and sign-aligned so that the score increases with age. Transfers within GTEx use the full prediction (training mean age plus the ridge coefficients applied to the z-scored profile), which is in years. The external cohorts and all reprogramming analyses use the frozen unit-length score, which orders samples identically but has no intercept or scale and is not in years. A PLS-1 model was fit in parallel for comparison and is neither the frozen ruler nor any gate. Within-site fibroblast cross-validation, grouped by sequencing batch with up to five folds, reports ρ but does not choose alpha. The ruler was never refit on any evaluation cohort.

#### 18. `P74` — `30`

Searched: alpha grid value is not stored; results/fibro/freeze.json has alpha 0.01 only

> The ruler is a ridge regression of donor chronological age on z-scored gene expression, fit through the n × n Gram matrix, and frozen with its coefficients and its training mean and standard deviation before any evaluation. The regularisation parameter was selected over the grid 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000 by minimising a leave-one-out residual sum of squares computed from the hat-matrix diagonal, not by grouped cross-validation. The frozen fibroblast ruler stores alpha = 0.01, the smallest value on the grid, so a smaller optimum cannot be excluded. The frozen fibroblast vector is the ridge coefficient vector, unit-normalised and sign-aligned so that the score increases with age. Transfers within GTEx use the full prediction (training mean age plus the ridge coefficients applied to the z-scored profile), which is in years. The external cohorts and all reprogramming analyses use the frozen unit-length score, which orders samples identically but has no intercept or scale and is not in years. A PLS-1 model was fit in parallel for comparison and is neither the frozen ruler nor any gate. Within-site fibroblast cross-validation, grouped by sequencing batch with up to five folds, reports ρ but does not choose alpha. The ruler was never refit on any evaluation cohort.

#### 19. `P74` — `100`

Searched: alpha grid value is not stored; results/fibro/freeze.json has alpha 0.01 only

> The ruler is a ridge regression of donor chronological age on z-scored gene expression, fit through the n × n Gram matrix, and frozen with its coefficients and its training mean and standard deviation before any evaluation. The regularisation parameter was selected over the grid 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000 by minimising a leave-one-out residual sum of squares computed from the hat-matrix diagonal, not by grouped cross-validation. The frozen fibroblast ruler stores alpha = 0.01, the smallest value on the grid, so a smaller optimum cannot be excluded. The frozen fibroblast vector is the ridge coefficient vector, unit-normalised and sign-aligned so that the score increases with age. Transfers within GTEx use the full prediction (training mean age plus the ridge coefficients applied to the z-scored profile), which is in years. The external cohorts and all reprogramming analyses use the frozen unit-length score, which orders samples identically but has no intercept or scale and is not in years. A PLS-1 model was fit in parallel for comparison and is neither the frozen ruler nor any gate. Within-site fibroblast cross-validation, grouped by sequencing batch with up to five folds, reports ρ but does not choose alpha. The ruler was never refit on any evaluation cohort.

#### 20. `P74` — `300`

Searched: alpha grid value is not stored; results/fibro/freeze.json has alpha 0.01 only

> The ruler is a ridge regression of donor chronological age on z-scored gene expression, fit through the n × n Gram matrix, and frozen with its coefficients and its training mean and standard deviation before any evaluation. The regularisation parameter was selected over the grid 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000 by minimising a leave-one-out residual sum of squares computed from the hat-matrix diagonal, not by grouped cross-validation. The frozen fibroblast ruler stores alpha = 0.01, the smallest value on the grid, so a smaller optimum cannot be excluded. The frozen fibroblast vector is the ridge coefficient vector, unit-normalised and sign-aligned so that the score increases with age. Transfers within GTEx use the full prediction (training mean age plus the ridge coefficients applied to the z-scored profile), which is in years. The external cohorts and all reprogramming analyses use the frozen unit-length score, which orders samples identically but has no intercept or scale and is not in years. A PLS-1 model was fit in parallel for comparison and is neither the frozen ruler nor any gate. Within-site fibroblast cross-validation, grouped by sequencing batch with up to five folds, reports ρ but does not choose alpha. The ruler was never refit on any evaluation cohort.

#### 21. `P74` — `1000`

Searched: alpha grid value is not stored; results/fibro/freeze.json has alpha 0.01 only

> The ruler is a ridge regression of donor chronological age on z-scored gene expression, fit through the n × n Gram matrix, and frozen with its coefficients and its training mean and standard deviation before any evaluation. The regularisation parameter was selected over the grid 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000 by minimising a leave-one-out residual sum of squares computed from the hat-matrix diagonal, not by grouped cross-validation. The frozen fibroblast ruler stores alpha = 0.01, the smallest value on the grid, so a smaller optimum cannot be excluded. The frozen fibroblast vector is the ridge coefficient vector, unit-normalised and sign-aligned so that the score increases with age. Transfers within GTEx use the full prediction (training mean age plus the ridge coefficients applied to the z-scored profile), which is in years. The external cohorts and all reprogramming analyses use the frozen unit-length score, which orders samples identically but has no intercept or scale and is not in years. A PLS-1 model was fit in parallel for comparison and is neither the frozen ruler nor any gate. Within-site fibroblast cross-validation, grouped by sequencing batch with up to five folds, reports ρ but does not choose alpha. The ruler was never refit on any evaluation cohort.

#### 22. `P74` — `3000`

Searched: alpha grid value is not stored; results/fibro/freeze.json has alpha 0.01 only

> The ruler is a ridge regression of donor chronological age on z-scored gene expression, fit through the n × n Gram matrix, and frozen with its coefficients and its training mean and standard deviation before any evaluation. The regularisation parameter was selected over the grid 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000 by minimising a leave-one-out residual sum of squares computed from the hat-matrix diagonal, not by grouped cross-validation. The frozen fibroblast ruler stores alpha = 0.01, the smallest value on the grid, so a smaller optimum cannot be excluded. The frozen fibroblast vector is the ridge coefficient vector, unit-normalised and sign-aligned so that the score increases with age. Transfers within GTEx use the full prediction (training mean age plus the ridge coefficients applied to the z-scored profile), which is in years. The external cohorts and all reprogramming analyses use the frozen unit-length score, which orders samples identically but has no intercept or scale and is not in years. A PLS-1 model was fit in parallel for comparison and is neither the frozen ruler nor any gate. Within-site fibroblast cross-validation, grouped by sequencing batch with up to five folds, reports ρ but does not choose alpha. The ruler was never refit on any evaluation cohort.

#### 23. `P74` — `10000`

Searched: alpha grid value is not stored; results/fibro/freeze.json has alpha 0.01 only

> The ruler is a ridge regression of donor chronological age on z-scored gene expression, fit through the n × n Gram matrix, and frozen with its coefficients and its training mean and standard deviation before any evaluation. The regularisation parameter was selected over the grid 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000 by minimising a leave-one-out residual sum of squares computed from the hat-matrix diagonal, not by grouped cross-validation. The frozen fibroblast ruler stores alpha = 0.01, the smallest value on the grid, so a smaller optimum cannot be excluded. The frozen fibroblast vector is the ridge coefficient vector, unit-normalised and sign-aligned so that the score increases with age. Transfers within GTEx use the full prediction (training mean age plus the ridge coefficients applied to the z-scored profile), which is in years. The external cohorts and all reprogramming analyses use the frozen unit-length score, which orders samples identically but has no intercept or scale and is not in years. A PLS-1 model was fit in parallel for comparison and is neither the frozen ruler nor any gate. Within-site fibroblast cross-validation, grouped by sequencing batch with up to five folds, reports ρ but does not choose alpha. The ruler was never refit on any evaluation cohort.

#### 24. `P75` — `0.30`

Searched: searched results/fibro/manifest.json key normalization; stored text is 'TMM then log2-CPM prior.count=2' and does not contain trim 0.30 or 0.05

> Training data. GTEx v10 cultured fibroblasts, 652 donors including the seven donors at the small collection site D1, 23,485 genes after filtering. Genes were retained at a count of at least 6 in at least 20% of samples. Counts were TMM-normalised (edgeR-style, log-ratio trim 0.30, absolute-expression trim 0.05, factors scaled to geometric mean 1), converted to log2 counts per million with prior count 2, and z-scored with the training mean and standard deviation; standard deviations below 1e-12 were set to 1. Site D1 is included in the freeze and excluded from the B1 to C1 transfer and from within-site cross-validation, which is stated rather than silently applied.

#### 25. `P75` — `0.05`

Searched: searched results/fibro/manifest.json key normalization; stored text is 'TMM then log2-CPM prior.count=2' and does not contain trim 0.30 or 0.05

> Training data. GTEx v10 cultured fibroblasts, 652 donors including the seven donors at the small collection site D1, 23,485 genes after filtering. Genes were retained at a count of at least 6 in at least 20% of samples. Counts were TMM-normalised (edgeR-style, log-ratio trim 0.30, absolute-expression trim 0.05, factors scaled to geometric mean 1), converted to log2 counts per million with prior count 2, and z-scored with the training mean and standard deviation; standard deviations below 1e-12 were set to 1. Site D1 is included in the freeze and excluded from the B1 to C1 transfer and from within-site cross-validation, which is stated rather than silently applied.

#### 26. `P75` — `1`

Searched: searched results/fibro/manifest.json normalization string for geometric mean; it is not stored there

> Training data. GTEx v10 cultured fibroblasts, 652 donors including the seven donors at the small collection site D1, 23,485 genes after filtering. Genes were retained at a count of at least 6 in at least 20% of samples. Counts were TMM-normalised (edgeR-style, log-ratio trim 0.30, absolute-expression trim 0.05, factors scaled to geometric mean 1), converted to log2 counts per million with prior count 2, and z-scored with the training mean and standard deviation; standard deviations below 1e-12 were set to 1. Site D1 is included in the freeze and excluded from the B1 to C1 transfer and from within-site cross-validation, which is stated rather than silently applied.

#### 27. `P75` — `SCI:1e-12`

Searched: searched results/fibro/manifest.json normalization and gene_filter strings for 1e-12; the stored strings are 'TMM then log2-CPM prior.count=2' and '>=6 counts in >=20% of samples'

> Training data. GTEx v10 cultured fibroblasts, 652 donors including the seven donors at the small collection site D1, 23,485 genes after filtering. Genes were retained at a count of at least 6 in at least 20% of samples. Counts were TMM-normalised (edgeR-style, log-ratio trim 0.30, absolute-expression trim 0.05, factors scaled to geometric mean 1), converted to log2 counts per million with prior count 2, and z-scored with the training mean and standard deviation; standard deviations below 1e-12 were set to 1. Site D1 is included in the freeze and excluded from the B1 to C1 transfer and from within-site cross-validation, which is stated rather than silently applied.

#### 28. `P75` — `1`

Searched: searched results/fibro/manifest.json for the standard-deviation floor replacement; 1e-12 and the replacement value 1 are not in that file

> Training data. GTEx v10 cultured fibroblasts, 652 donors including the seven donors at the small collection site D1, 23,485 genes after filtering. Genes were retained at a count of at least 6 in at least 20% of samples. Counts were TMM-normalised (edgeR-style, log-ratio trim 0.30, absolute-expression trim 0.05, factors scaled to geometric mean 1), converted to log2 counts per million with prior count 2, and z-scored with the training mean and standard deviation; standard deviations below 1e-12 were set to 1. Site D1 is included in the freeze and excluded from the B1 to C1 transfer and from within-site cross-validation, which is stated rather than silently applied.

#### 29. `P85` — `1`

Searched: searched results/fibro/manifest.json and results/toward/stats.csv for the string 201 or the formula (exceedances+1)/201; they store n_perm=200 and p=0.004975124378109453

> All nulls use 200 draws and an empirical p-value of (exceedances + 1) / 201. Three constructions appear, and they are not interchangeable. A donor-label permutation shuffles donor age within the training bank, holding the expression matrix fixed, and is the null for every transfer statistic. A gene-permutation null shuffles the entries of the frozen target direction across genes, and tests whether a cosine is larger than for a direction with the same coefficient distribution. An anchor null rebuilds the young-minus-old direction from equal-size random splits of all GTEx fibroblast donors, ignoring age, and tests whether an alignment is specifically with age rather than with any donor contrast.

#### 30. `P85` — `201`

Searched: searched results/fibro/manifest.json and results/toward/stats.csv for the string 201 or the formula (exceedances+1)/201; they store n_perm=200 and p=0.004975124378109453

> All nulls use 200 draws and an empirical p-value of (exceedances + 1) / 201. Three constructions appear, and they are not interchangeable. A donor-label permutation shuffles donor age within the training bank, holding the expression matrix fixed, and is the null for every transfer statistic. A gene-permutation null shuffles the entries of the frozen target direction across genes, and tests whether a cosine is larger than for a direction with the same coefficient distribution. An anchor null rebuilds the young-minus-old direction from equal-size random splits of all GTEx fibroblast donors, ignoring age, and tests whether an alignment is specifically with age rather than with any donor contrast.

#### 31. `P86` — `2.5`

Searched: searched results/fibro/manifest.json and results/md3/t2_delta_rho.csv for 2.5 or 97.5; not present

> Confidence intervals are 200-draw percentile bootstraps at the 2.5th and 97.5th percentiles. The resampling unit is the donor for transfer statistics and the cell for displacement statistics. Paired comparisons recompute both quantities on the same resampled donors, so that the interval is on the difference rather than on two separate intervals. Any interval that does not contain its own point estimate is flagged invalid and is not used as a gate; this happens for distance intervals in the high-dimensional displacement tests and is reported at each occurrence.

#### 32. `P86` — `97.5`

Searched: searched results/fibro/manifest.json and results/md3/t2_delta_rho.csv for 2.5 or 97.5; not present

> Confidence intervals are 200-draw percentile bootstraps at the 2.5th and 97.5th percentiles. The resampling unit is the donor for transfer statistics and the cell for displacement statistics. Paired comparisons recompute both quantities on the same resampled donors, so that the interval is on the difference rather than on two separate intervals. Any interval that does not contain its own point estimate is flagged invalid and is not used as a gate; this happens for distance intervals in the high-dimensional displacement tests and is reported at each occurrence.

#### 33. `P7` — `four`

Searched: searched results/fibro/stage1_transfer.csv, results/fibro2/ta_summary.json and results/paper_figs/cxg_true_donor_results.csv for a stored transfer count of 4; no such field. Table 1 prints five Spearman rows.

> We built a transcriptomic age ruler in cultured human dermal fibroblasts (652 GTEx donors) by supervised ridge regression on donor age, and froze it. It orders donors by age in cohorts that differ in laboratory, collection site and platform (Spearman ρ 0.43 to 0.56 across four transfers), but it does not predict age in calibrated years. In an independent bulk cohort of 82 donors, the frozen ruler tracks donor age better than the curated mesenchymal-drift (MD) signature used as an age readout by Lu et al. (paired bootstrap Δρ +0.47, 95% CI +0.24 to +0.70). In a second, single-cell cohort the difference points the same way, but once its samples are grouped into the 65 people they came from, the interval includes zero (+0.17, −0.12 to +0.41).

#### 34. `P12` — `twenty thousand`

Searched: searched results/fibro/freeze.json for a gene count of 20000; the stored n_genes is 23485, and the sentence does not say it is that count

> This form has a structural weakness: it records the direction a cell moves and not where it ends up. In a space of twenty thousand genes almost every direction leads away from the old state, and only a narrow cone of them leads toward young cells of the same type. A score anchored on not-old falls across that whole space. The condition for genuine approach is quantitative: a displacement ends closer to the young target only when its length is less than twice the target distance times the cosine of the angle between them. A directional score carries the angle and discards the length, so it cannot evaluate that condition even in principle. Three outcomes therefore become hard to tell apart:

#### 35. `P53` — `twice`

Searched: searched results/same/stats.csv and results/same/summary.json for a stored ratio of 2. Forward row has delta=130.89356478446803 and v_norm=130.68843692001553. No field stores twice or a ratio.

> In the dataset where the distinction can be tested, what is scored as rejuvenation is well explained by departure from the aged state, and is not accompanied by a measurable approach to a young cell of the same type. The aged donor's partially reprogrammed cells project two thirds of the way along the aged-to-young axis while ending twice as far from the young target as they started, and the same sign holds against a GTEx-defined young centroid. A dose-response control shows the statistic would have seen an approach had one occurred, and a reverse test shows the movement is not mutual convergence of two donors.

#### 36. `P72` — `three`

Searched: searched results/**/*.json for withdrawn; the only hit is results/md4/t1_reading.json key md3_status_withdrawn. No file stores a count of three withdrawn readings.

> Each analysis was specified in writing before any number was computed, including its statistics, nulls, positive controls, decision thresholds and the exact wording of each possible conclusion. A flag file was written before execution and the pre-registration text is reproduced verbatim in the results file for each analysis. Only the pre-registered reading whose conditions fired is reported as the conclusion of that analysis. Where a later analysis superseded an earlier reading, the earlier text was retained and marked withdrawn rather than edited; three readings were withdrawn in this way. Analyses run as diagnostics after a gate had fired are labelled as such and are not treated as gated results. One later correction is reported alongside the original result rather than in place of it: the regrouping of the fibroblast atlas into real donors, which was done after the pre-registered result was known and is reported as a correction of the analysis unit. Analyses of brain tissue and of a published transcription-factor screen were carried out in the same project under the same protocols; they do not bear on the question tested here and are not reported in this paper. Their pre-registrations, code and results, including a transfer that failed its gate, are in the archive.

## Stage 3

### Placeholders

- `P96`: `[GITHUB REPOSITORY URL]`
- `P96`: `[ZENODO DOI]`

Sentence: All input data are public and are listed in Table 3. Analysis code, the frozen instruments, the pre-registration documents with their timestamp flags, and every result table used in this paper are available at [GITHUB REPOSITORY URL] and archived at [ZENODO DOI]. The frozen fibroblast ruler is a single file containing the unit weight vector, the un-normalised coefficients and the training mean and standard deviation over 23,485 genes; the young and old GTEx anchors are provided as a separate file. The repository README lists the entry-point script that reproduces each table and figure. Raw input data are not redistributed and must be obtained from the sources in Table 3.

No other bracketed placeholder is in the draft.

### Named accessions, URLs, DOIs, and results files

The draft names no path under `results/` and no `http` URL. The two bracketed availability strings are placeholders, not resolvable links.

| Name in the draft | On disk |
| --- | --- |
| GSE226189 | yes — `results/fibro2/ta_GSE226189_near_miss_bulk_result.json` |
| GSE113957 | yes as a recorded non-score — `results/fibro2/ta_summary.json` text says the FPKM file is not integer counts |
| GSE297234 | yes — `results/survey/datasets.csv` and `results/toward/stats.csv` |
| GM00731, GM23815 | yes — `results/toward/stats.csv` column `cell_line` |
| CELLxGENE `a19d1667-a7b5-4556-9e5f-f9bfa690c0f1` | yes — `results/fibro2/ta_cxg_a19d1667_result.json` |
| dbGaP `phs000424` | yes — `results/fibro/stage1_report.txt` header `phs000424.v10` |
| bioRxiv `2024.09.05.611379` | no file in `results/` contains that id |
| `[ZENODO DOI]` | placeholder; nothing fetched |
| `[GITHUB REPOSITORY URL]` | placeholder; nothing fetched |

### Citations

Author-year citations in the running text are Lu et al. and Fleischer et al. Both are in the reference list.

Datasets used as data whose papers are not named in the running text:

- GSE226189 is used as a scored cohort. The reference list has no paper for that accession. Table row: Human skin fibroblast bulk RNA-seq | GSE226189 | External fibroblast transfer, 82 donors.
- The dermal fibroblast atlas is used. Ascensión & Izeta 2024 is in the reference list and is not named in the running text.
- The HALLMARK EMT gene set is used. Liberzon et al. 2015 is in the reference list and is not named in the running text.
- The sentence that states 1,836 transcription factors does not name Southard et al. Southard et al. 2025 is in the reference list. The count itself is the string in `results/survey/opened/bioproject_esummary.json`.

GTEx is named in the acknowledgements and GTEx Consortium 2020 is in the reference list. GSE297234 is introduced as Lu et al. GSE113957 is named as not scored; Fleischer et al. 2018 is cited for the age-up and age-down lists and is in the reference list.

The methods say "edgeR-style". Robinson & Oshlack 2010 is in the reference list. The methods say "Seurat". Hao et al. 2024 is in the reference list.

## Stage 4

### 1. Numerical discrepancies

1. **resolved.** The current abstract reads: "It orders donors by age in cohorts that differ in laboratory, collection site and platform (Spearman ρ 0.43 to 0.56 across four transfers), but it does not predict age in calibrated years."
2. **stale.** The GTEx bulk-cortex row and the cell "0.45" were cut. Table 1 is fibroblast transfers only.
3. **resolved.** "Within GTEx, where the ruler outputs a predicted age in years, uncalibrated R² is positive but modest (+0.29 and +0.25)."
4. **stale.** The sentence "approximately 83°" is not in this draft. The age-to-identity angle paragraph was cut.
5. **resolved.** Table 2 prints `+0.29 [+0.11, +0.49]` for the 93-pseudobulk age-up/down contrast.
6. **resolved.** "A five-fold positive control, in which held-out old GTEx donors are displaced toward the young centroid, gave fold cosines of 0.615, 0.714, 0.825, 0.849 and 0.874, clearing the pre-registered bar of 0.5."
7. **resolved.** "We note that the distance travelled is not monotonic across these states (271.7, 268.6, 546.4 z units), so cosine order and displacement order are not the same ordering."
8. **resolved.** "This distinction matters numerically: the change in distance at pluripotency is +184 on the GTEx axis and +382 within the experiment, and those are two different measurements of two different things."
9. **stale.** The CRISPRi ranking paragraph, including ZFX, was cut.
10. **stale.** The late-passage proliferation paragraph was cut.

The further "All transfers exceed donor-label permutation nulls" item is **stale**. That sentence, and the cortex and SEA-AD transfers it depended on, are not in this draft. Each fibroblast row now prints `p = 0.005`.

### 2. Placeholders still in the draft

- **resolved.** "Roger Kwon, independent researcher, Seoul, Republic of Korea. Correspondence: kih7920@gmail.com"
- **stale.** The version-and-date line is not in this draft.
- **resolved.** "The regularisation parameter was selected over the grid 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000 by minimising a leave-one-out residual sum of squares computed from the hat-matrix diagonal, not by grouped cross-validation. The frozen fibroblast ruler stores alpha = 0.01, the smallest value on the grid, so a smaller optimum cannot be excluded." The cortex preprocessing sentence was cut.
- **open.** "Accession-level detail and per-cohort sample tables are in the repository." Table 3 now names the accessions. The per-cohort sample tables are still not in the draft.
- **resolved.** "Analyses were implemented in Python 3.12.7 with numpy 2.5.3, scipy 1.18.1, pandas 3.0.5, scikit-learn 1.9.1, anndata 0.13.3.post0 and h5py 3.16.0, and figures were drawn with matplotlib 3.11.2." Those strings match `results/paper_figs/environment_versions.txt`.
- **open.** "available at [GITHUB REPOSITORY URL] and archived at [ZENODO DOI]."
- **resolved.** "R.K. is the sole author. He conceived the study, specified the analyses and their pre-registered decision rules, directed their execution, and is responsible for the manuscript."
- **resolved.** "The author thanks the investigators and consortia who generated and publicly released the datasets analysed here, including the GTEx Consortium, the authors of the studies integrated in the dermal fibroblast atlas, and the authors of GSE226189 and GSE297234."
- **resolved.** "The author declares no competing interests."

### 3. Methods that are still missing or wrong

- **resolved.** "Counts were TMM-normalised (edgeR-style, log-ratio trim 0.30, absolute-expression trim 0.05, factors scaled to geometric mean 1), converted to log2 counts per million with prior count 2." The draft no longer describes a cortex ruler in the same paragraph.
- **stale.** The ≥65 SEA-AD parenthetical was cut.
- **stale.** The R_rej reproduction paragraph was cut with the screen.
- **resolved.** "The two cosines in this paper are not interchangeable and are never combined." The same paragraph states +184 on the GTEx axis and +382 within the experiment.
- **stale.** The Sengstack +20.6 zero-mixture result was cut. The screen is not in this draft.
- **resolved.** "And the bootstrap intervals on the distance values do not contain their point estimates and are flagged invalid, so the monotonicity claim rests on the point estimates alone, as stated."
- **resolved.** The package versions in the software paragraph match `results/paper_figs/environment_versions.txt`.

### 4. Figures that need computation

- **stale.** Figure 1 no longer has brain-bank, BA9→Cortex, or B1↔C1-as-brain panels. The caption is the fibroblast transfer figure.
- **stale.** The SEA-AD donor panel was cut.
- **stale.** The passage-versus-donor-age null panel was cut. Figure 5 is the same-platform reprogramming test.
- **resolved.** The extrapolation scatter is Figure 2. "Each point is one donor × cluster × timepoint pseudobulk from the reprogramming dataset with at least 20 cells (n = 22). Distance is the nearest-neighbour Euclidean distance to the GTEx training donors in a 50-component PCA space, read as computed and not recomputed here."

### 5. References cited in the text and not in the reference list

- **resolved.** GTEx Consortium 2020 is reference item P108.
- **stale.** SEA-AD is not in this draft.
- **resolved.** Hao et al. 2024, the Seurat paper, is in the reference list.
- **stale.** The cell-cycle gene-set citation belonged to the screen section, which was cut.
- **resolved.** Horvath 2013, Hannum et al. 2013, and Kabacik et al. 2022 are in the reference list.
- **resolved.** Takahashi & Yamanaka 2006, Ocampo et al. 2016, and Gill et al. 2022 are in the reference list.
- **stale.** The ROSMAP / DLPFC cohort was cut.
- **resolved.** "Fleischer, J. G. et al. Predicting age from the transcriptome of human dermal fibroblasts. Genome Biol. 19, 221 (2018)."
- **open.** GSE226189 is used and is not given a paper in the reference list. "Human skin fibroblast bulk RNA-seq | GSE226189 | External fibroblast transfer, 82 donors."
- **resolved.** GSE113957 is named in Table 3, and Fleischer et al. 2018 is in the reference list.
- **resolved.** "Liberzon, A. et al. The Molecular Signatures Database hallmark gene set collection. Cell Syst. 1, 417–425 (2015)."
- **open.** Southard et al. 2025 is in the reference list, and the 1,836 sentence does not name it. "Genome-scale CRISPRa Perturb-seq in primary human fibroblasts now covers 1,836 transcription factors with a transcriptome-wide readout, enough to rank candidate interventions by approach to a young target rather than by departure from an old one, with a dose-response control built in from the start."

## Fired keys

None.

## Pre-registration (verbatim, written before any check)

```
TASK: VERIFY — audit every number in the preprint draft against its source file.

Nature of this task: a read-only audit of files already in this repository.
No statistic is computed, no analysis is rerun, and nothing is edited —
including the preprint itself.

Write FINDINGS_VERIFY.md, PROGRESS_VERIFY.md, results/verify/*.
Create no other files. Modify no existing file. Do not refit the frozen ruler.
Do not rerun any analysis to regenerate a number: every check reads a stored
result file as it is on disk.

Write this file verbatim into FINDINGS_VERIFY.md and results/verify/PREREG.flag
before any check.

--- INPUT ---
The current preprint draft: the most recent rejuvenation_readouts_preprint_3.docx
in the repository paper/. Report its filename, size, md5 and modification time, and
list every other draft version found so the author can confirm the right one
was audited. STOP `no_draft` if none exists.

--- STAGE 1: EXTRACT ---
Extract every numeric claim in the draft: correlations, confidence intervals,
distances, counts, sample sizes, fractions, angles, p and q values, thresholds,
accessions, file sizes. Include the abstract, tables, figure captions and
methods. Write them to results/verify/claims.csv with the number, the sentence
containing it, and its location in the document. Report the total.

--- STAGE 2: TRACE ---
For each claim, find the file under results/ that produced it. Record the path
and the exact field or row. Read the number; never compute it. Mark each claim:
  match        — draft and source agree at the precision printed
  rounding     — differ only in the last printed digit
  mismatch     — both exist and differ; record both values verbatim
  untraceable  — no source found; record what was searched
  not_a_claim  — a date, page number or similar
Report the count of each. List every mismatch and every untraceable claim in
full, with the draft sentence quoted exactly. Do not summarise them.

--- STAGE 3: PLACEHOLDERS AND REFERENCES ---
Report every bracketed placeholder still in the draft, verbatim with location.
Report every results file, accession, URL and DOI the draft names, and whether
it exists on disk. Do not fetch anything requiring credentials.
Report every reference cited in the text but absent from the reference list,
and every dataset used as data but not cited.

--- STAGE 4: GAPS.md TRIAGE ---
GAPS.md was written against an earlier draft and may be out of date. For each
item in all five of its sections report exactly one verdict:
  resolved — the current draft carries the corrected value; quote it
  open     — the discrepancy or gap is still present; quote the draft sentence
  stale    — the item refers to text no longer in the draft; say what was cut
Do not edit GAPS.md.

--- WHAT DOES NOT COUNT ---
- editing the preprint, GAPS.md, or any result file.
- recomputing a number instead of reading it.
- judging a mismatch acceptable; report it and let the author decide.
- guessing a source file; if the trace is uncertain, mark it untraceable.
- collapsing several mismatches into one line.

--- OUTPUT ---
FINDINGS_VERIFY.md: this file verbatim; draft identity and other versions found;
Stage 1 count; Stage 2 verdict counts plus the full mismatch and untraceable
lists; Stage 3 placeholders, references and broken links; Stage 4 triage of all
five GAPS sections; fired keys.
PROGRESS_VERIFY.md: stop status, and a short ordered list of what must be fixed
before posting, hardest first.
results/verify/claims.csv: every claim with its verdict and source path.
Record every failure verbatim.
```

superseded by the docx, kept for history.



# Rejuvenation readouts cannot distinguish younger from different

*A one-dimensional age score cannot separate rejuvenation from loss of cell identity. Distance to a young target of the same cell type can.*

**Author:** name, independent researcher, Seoul, Republic of Korea. Correspondence: email.

**Preprint draft, not peer reviewed. Version v1, date.**

## Abstract

Interventions that claim to rejuvenate cells are judged by scores that fall when cells "get younger". Because aging is itself a drift in cell identity, such a score also falls when a cell stops being the cell type it was. Rejuvenation and loss of identity are therefore indistinguishable on a one-dimensional readout. We formalise the alternative: the safe target is a point, not a direction, and progress is the distance to young cells of the same type, with departure from that type counted as failure.

We built transcriptomic age rulers in human dorsolateral prefrontal cortex (233 donors) and cultured fibroblasts (652 GTEx donors) by supervised regression on donor age. Both transfer to independent cohorts, laboratories and platforms as a direction (Spearman ρ 0.45–0.57) but not as a calibrated age. In two independent fibroblast cohorts, the learned ruler predicts donor age better than the curated mesenchymal-drift (MD) signature used as the age readout by Lu et al. (paired bootstrap Δρ +0.20 and +0.47, both excluding zero).

Applying the distance-to-target framework to two flagship datasets gives a consistent result. In OSKM partial reprogramming of fibroblasts (Lu et al.), partially reprogrammed cells end up farther from young fibroblasts than the aged cells they came from, whether "young" is defined by GTEx donors or by a 22-year-old donor assayed in the same experiment; their apparent rejuvenation is movement away from the old state, not toward the young one. In a Perturb-seq screen for rejuvenating transcription factors (Sengstack et al.), whose readout is transcriptome-wide and which controls for dedifferentiation, we reproduce the published ranking under the authors' own statistic (Spearman −0.98 between their score and ours); but that statistic is a correlation, which registers orientation without destination, the screen-wide score correlates with a proliferation signature at ρ 0.83–0.86, and the replicative-passage axis used to define aging is indistinguishable from random with respect to a donor-age axis (cos −0.04, p = 1.00).

These results suggest that current rejuvenation readouts largely measure departure from the aged state and proliferative capacity rather than approach to a young cell of the same type. We provide frozen rulers, young-target anchors and a pre-registered test harness so that any intervention can be scored on this basis.

## Introduction

Partial reprogramming, transcription-factor perturbation and other rejuvenation interventions are evaluated by asking whether a cell's molecular profile has moved in the direction opposite to aging. Epigenetic clocks, curated aging signatures and correlation-based reversal scores all share this form: a single number that falls when the cell is judged younger.

This form has a structural weakness. Loss of cellular identity is itself among the most reproducible features of aging, and the same is true of aging at the tissue level, where the composition of cell states drifts. An intervention that pushes a cell out of its differentiated state therefore reduces an aging score for a reason that has nothing to do with restoring youthful function. Two very different outcomes become indistinguishable:

- **Approach to young.** The cell becomes more similar to a young cell of the same type. This is the intended outcome.
- **Departure from old.** The cell becomes less similar to an old cell by becoming something else, for example a less differentiated or more proliferative state. This is overshoot, and in the limit it is dedifferentiation and loss of tissue function.

The distinction is not hypothetical for reprogramming specifically. Complete OSKM reprogramming proceeds through a mesenchymal-to-epithelial transition and ends in pluripotency; the entire trajectory moves away from the aged fibroblast state, and any score anchored only on "not old" will register the whole of it as rejuvenation.

We therefore replace the direction with a target. For a given cell type there is a point in expression space that the intervention should approach: this cell type, young. Movement toward that point is progress; movement past it is failure. This converts identity from a competing axis into a stopping condition and converts age from a direction into a distance.

Operationalising this requires three things: a ruler that orders cells by donor age and demonstrably transfers to new data; an anchor for the young target of the same cell type, measured on the same platform where possible; and a statistic that separates "closer to young" from "farther from old", together with controls that detect the confounds that make the two look alike. We build these, validate them against the existing readout, and apply them to two recent datasets that report rejuvenation in human fibroblasts.

## Results



### A learned age ruler transfers across cohorts as a direction

We trained ridge-regression rulers on donor chronological age from transcriptome-wide expression, separately in human cortex and in cultured dermal fibroblasts, and froze them before any evaluation. Each ruler was then applied unchanged to cohorts differing in laboratory, tissue bank, dissociation protocol, preservation pathway and sequencing platform.


| Tissue         | Trained on                        | Held-out cohort                                   | ρ vs donor age                  |
| -------------- | --------------------------------- | ------------------------------------------------- | ------------------------------- |
| Cortex (DLPFC) | snRNA-seq, 233 donors             | Independent brain bank                            | 0.57                            |
| Cortex         | DLPFC                             | GTEx bulk cortex (different preservation pathway) | 0.45                            |
| Cortex         | DLPFC                             | SEA-AD middle temporal gyrus                      | 0.48 (attenuated in donors ≥65) |
| Fibroblast     | GTEx collection site B1 (n = 421) | GTEx site C1 (n = 224)                            | 0.53                            |
| Fibroblast     | GTEx site C1                      | GTEx site B1                                      | 0.56                            |
| Fibroblast     | GTEx, 652 donors                  | CELLxGENE 10x cohort (n = 93)                     | 0.46                            |
| Fibroblast     | GTEx, 652 donors                  | GSE226189 bulk (n = 82)                           | 0.55                            |


All transfers exceed donor-label permutation nulls. Transfer is directional only: predicted age is not calibrated in years, and the coefficient of determination remains negative in every external cohort. We therefore treat the ruler as an ordering instrument and never as an age estimate, and every statistic below uses it as a direction or a distance rather than as a score in years.

In cortex we additionally fit a cell-identity ruler per cell type. The age and identity axes are close to orthogonal (approximately 83°), confirming that they are separable directions within a single expression space rather than two names for the same drift.

One property of the ruler matters for everything that follows. Its score is correlated with distance from its own training distribution (ρ = −0.63 between nearest-neighbour distance to the training data and predicted age in the reprogramming dataset). Cells far outside the training manifold are scored as young irrespective of biology. This extrapolation artefact is the concrete mechanism by which a one-dimensional score mistakes departure for rejuvenation, and it applies to our ruler exactly as it applies to the readouts we examine below.

### The ruler outperforms a curated mesenchymal-drift signature

Lu et al. use a mesenchymal-drift (MD) signature — the HALLMARK epithelial-mesenchymal transition gene set augmented with SNAI1, ZEB1, ZEB2, TWIST1 and TWIST2 — as their readout of aging, and separately report age-up and age-down gene lists. If MD is to serve as an age readout, it should predict donor age. We compared all three instruments, applied exactly as published, on two fibroblast cohorts that none of them was fitted on.


| Cohort                  | Ruler ρ | MD ρ  | Age-up minus age-down ρ | Δρ vs MD 95% CI    | Δρ vs age-up/down 95% CI |
| ----------------------- | ------- | ----- | ----------------------- | ------------------ | ------------------------ |
| CELLxGENE 10x (n = 93)  | 0.455   | 0.253 | 0.161                   | +0.20 +0.01, +0.44 | +0.30 +0.11, +0.49       |
| GSE226189 bulk (n = 82) | 0.549   | 0.082 | 0.037                   | +0.47 +0.24, +0.70 | +0.51 +0.25, +0.77       |


Differences are paired: Δρ is recomputed on the same resampled donors (B = 200), so the interval reflects the difference directly rather than the overlap of two separate intervals. A third pre-registered cohort (GSE113957) could not be scored because only FPKM values were available, and is reported as not scored rather than as a result.

MD nonetheless behaves as Lu et al. describe within the reprogramming trajectory: it falls as cells leave the fibroblast state. Across clusters in their own OSKM dataset, the MD score is almost perfectly anticorrelated with a pluripotency-minus-fibroblast score (ρ = −0.94 and −0.96 in the aged and young donor respectively). MD genes are largely mesenchymal identity genes, so a falling MD score is primarily a statement about how fibroblast-like a cell remains. This is consistent with MD being a good marker of the reprogramming trajectory and a weak marker of donor age.

### Partially reprogrammed fibroblasts move away from old, not toward young

We analysed GSE297234 (Lu et al.), an OSKM reprogramming time course in dermal fibroblasts from a 96-year-old and a 22-year-old donor sampled at days 0, 3, 7 and 10, using the published cell-state assignments. Displacement was measured from the aged donor's day-0 fibroblast pseudobulk to its partially reprogrammed pseudobulk (5,832 cells), in the frozen ruler's gene space.

**Young defined by GTEx donors.** Anchors were the mean expression of GTEx fibroblast donors aged 20–39 (n = 113) and 60–79 (n = 231), with 40–59 excluded and reported. The displacement is almost orthogonal to the young-minus-old axis (cosine +0.077). That small alignment survives a gene-permutation null (p = 0.005) but not a null in which GTEx donors are split at random ignoring age (p = 0.39), so it is not specifically an age alignment. The cells move farther from the young centroid (+9.4) and farther from the old centroid (+10.9): they recede from both. A five-fold positive control in which held-out old GTEx donors are displaced toward young gave cosines of at least 0.62, confirming the statistic detects a genuine approach.

**Young defined within the experiment.** To remove the bulk-versus-single-cell platform gap, we repeated the test with the 22-year-old donor's untreated day-0 cells as the target. Aged partially reprogrammed cells again end farther from young (+131). A dose-response positive control, mixing young into aged day-0 cells at 10%, 25% and 50%, produced monotonically increasing progress and decreasing distance, so the statistic could have detected an approach had one occurred.

**Convergence versus rejuvenation.** A treatment that merely erased donor-specific expression would move both donors toward each other and mimic rejuvenation. Reversing the test, the young donor's partially reprogrammed cells do not move toward the aged donor (progress −0.06), giving a clear asymmetry. However, the aged-to-young donor axis is itself nearly orthogonal to the GTEx donor-age axis (cosine 0.039), so the partial alignment along it is better explained by differences between these two cell lines, including proliferative state, than by age.

The trajectory as a whole illustrates the central point. Cosine alignment with the youth axis increases monotonically from partially reprogrammed to early pluripotent to pluripotent cells, in the same order as the total distance travelled, while distance to young fibroblasts increases throughout (pluripotency: +382). The farther cells travel from the aged fibroblast state, the more rejuvenated they appear on a directional score and the less they resemble a young fibroblast.

### A published TF rejuvenation screen tracks proliferation, and its aging axis is not donor aging

Sengstack et al. screened 200 transcription factors by CRISPRa and CRISPRi in replicatively aged human dermal fibroblasts and scored rejuvenation as Rrej, the Pearson correlation between the late-versus-early passage fold changes and the perturbation-versus-non-targeting fold changes. Their readout is transcriptome-wide rather than a curated signature, and they explicitly guard against the dedifferentiation confound: the four validated factors leave fibroblast identity genes, telomere length, DNA damage and methylation age largely unchanged, show no resemblance to oncogenic transformation, and partially reverse the MD signature of Lu et al. Their screen is therefore not vulnerable to the criticism in the preceding section, and we treat it as the stronger of the two readouts.

A correlation nonetheless measures orientation, not destination. Rrej is invariant to the magnitude of the perturbation's effect and to any component orthogonal to the aging contrast, so a displacement can have a strongly negative Rrej while ending farther from the early-passage state than it began — the situation we document above for OSKM. Identity checks exclude one specific failure mode, dedifferentiation toward stemness, but they do not establish approach to a target. Two questions therefore remain open in their data: whether the perturbations approach the young state, and whether their aging axis corresponds to donor aging.

**The published ranking reproduces; the effect sizes do not.** Implementing Rrej as defined in their methods and comparing it with our own cosine against the passage-reversal direction, the two statistics order the screen almost identically (Spearman −0.984 across 340 perturbations; the sign reflects opposite conventions). E2F3 ranks first among CRISPRa perturbations and ZFX first among CRISPRi perturbations under both. Absolute values are roughly half those published (E2F3 −0.28 versus −0.53; ZFX −0.36 versus −0.51), and STAT3, one of the four factors selected for extensive validation, ranks 49 of 173 in our recomputation (−0.155 versus a published −0.41). The papers' methods do not state which genes enter the correlation, and we did not search for a gene filter that recovers the published table, since fitting a cutoff to a target result is not a reproduction.

**The screen tracks proliferation.** Scoring every non-targeting cell with S and G2M cell-cycle gene sets and defining a proliferation axis from the top and bottom quartiles, the passage-reversal direction is itself substantially proliferative (cosine 0.567 and 0.634 in the two modalities). Across all perturbations, the rejuvenation cosine correlates with the perturbation's shift in cell-cycle score at ρ = 0.828 (CRISPRa) and 0.862 (CRISPRi). After projecting out the proliferation axis, the named hits retain little signal (E2F3 +0.248 → +0.068; ZFX +0.316 → +0.039; STAT3 +0.155 → −0.002). Both correlations are unchanged under either choice of late-passage reference.

**Replicative passage is not donor aging.** The passage-reversal direction is indistinguishable from random with respect to the GTEx donor-age axis (cosine −0.040, permutation p = 1.00; −0.006 and p = 0.74 under the alternative reference). The cells are neonatal foreskin fibroblasts aged by serial passaging. Whatever the screen reverses, it is not the transcriptional signature that distinguishes old from young human donors in our data.

We do not report per-factor toward-young classifications for this dataset. Our distance statistic carries a positive offset of approximately +20 in this matrix, established by the zero-mixture point of the positive control, and the named hits' distances (+13 to +28) fall within that offset. Whether these perturbations approach the early-passage state therefore remains untested here; it is not answered by Rrej, and answering it would require a distance statistic calibrated on this matrix. The proliferation and donor-age results above do not depend on that statistic.

We note that the authors reach compatible conclusions from their own data. They report that most top perturbations upregulate cell-cycle TF modules and shift cells out of G1, that cell types sharing the rejuvenation signature in parabiosis are enriched for mitotically active populations, and that cell-cycle slowdown is a feature of replicative aging that "may not be an important feature of in vivo aging". They address the confound by recomputing Rrej within a single cell-cycle phase and find residual signal for several factors, including EZH2. Our contribution is to quantify the confound across the whole screen rather than for selected factors, and to test the aging axis itself against donor age.

## Discussion

Two independent datasets, two different interventions and two different readouts produce the same pattern: what is scored as rejuvenation is well explained by departure from the aged state and by proliferative capacity, and is not accompanied by measurable approach to a young cell of the same type.

This is not a claim that partial reprogramming or transcription-factor perturbation fails to rejuvenate. Both papers report functional improvements — proliferation, proteostasis, mitochondrial activity, and in aged mouse liver a reduction in steatosis and fibrosis and improved glucose tolerance — that transcriptome distance cannot adjudicate, and the in vivo results in particular lie outside what these data can test. Nor is it a claim that the authors were unaware of these confounds; Sengstack et al. in particular test for dedifferentiation, control for cell-cycle phase, and note the limits of the replicative model themselves. The claim is narrower: the readouts used to rank candidate interventions cannot distinguish approach from departure, the checks that are performed exclude only specific failure modes rather than establishing approach, and in both datasets examined the positive evidence for approach is absent where it should be strongest.

The mechanism differs by readout and is the same in kind. Directional scores learned by regression are confounded with distance from the training distribution, so cells that leave the manifold score young. Curated identity-gene signatures such as MD fall whenever the cell type is abandoned. Correlation-based reversal scores such as Rrej are invariant to magnitude and to orthogonal displacement, so they register orientation without destination; this is why a transcriptome-wide readout does not by itself solve the problem, and why adding genes to a signature does not either. Replicative-senescence models make proliferation and youth nearly the same axis, so any proliferative perturbation scores as rejuvenating. Each is a specific instance of scoring a direction where a target is required.

The remedy we propose is cheap to adopt. Define the young target on the same platform as the intervention, report distance to it rather than progress along a direction, include a dose-response positive control so that a null result can be distinguished from an insensitive assay, and report an identity check and a proliferation-residualised statistic alongside the headline number. None of this requires new data; all four analyses here were performed on public data.

The obvious next step is to apply the framework prospectively rather than retrospectively. Genome-scale CRISPRa Perturb-seq in primary human fibroblasts now covers 1,836 transcription factors with transcriptome-wide readout, which is sufficient to rank every human transcription factor by approach to a young target rather than by departure from an old one. Whether any single factor produces measurable approach, and whether the factors that do so overlap with those nominated by directional scores, is an open and directly testable question.

## Methods

**Pre-registration.** Every analysis was specified in writing, including its statistics, nulls, positive controls, decision thresholds and the exact wording of each possible conclusion, before any number was computed. Pre-registration text is reproduced verbatim in each results file and timestamped by a flag written before execution. Superseded conclusions are retained and marked, not edited.

**Rulers.** Ridge regression on donor chronological age from transcriptome-wide log-expression, fit once per tissue and frozen, with coefficients and the training mean and standard deviation stored. No ruler was refit on any evaluation cohort. Add: regularisation selection, exact preprocessing, cell-type handling in cortex.

**Expression processing.** Pseudobulk sums per group, TMM-normalised log2 counts per million with prior count 2 computed across the analysis panel, then z-scored using the frozen training mean and standard deviation; genes absent from a dataset are set to zero.

**Statistics.** For a displacement *d* from an origin pseudobulk and a target direction *v*, we report cosine alignment, progress along the axis (*d·v* / ‖*v*‖²) and change in Euclidean distance to the target. Approach requires both positive alignment beating its nulls and decreasing distance; alignment alone is reported as departure only.

**Nulls.** Gene-permuted target directions; random donor splits ignoring age; and, for perturbation screens, empirical nulls constructed from size-matched random subsets of non-targeting cells, with Benjamini–Hochberg control within modality. Confidence intervals are percentile bootstrap (B = 200); any interval not containing its point estimate is flagged invalid and is not used as a gate.

**Positive controls.** Each test must first recover a known effect — held-out old donors displaced toward young, or synthetic mixtures of young into old cells at defined fractions — before any experimental result is interpreted. The zero-mixture point is reported as the noise floor.

**Datasets.** GTEx v10 cultured fibroblasts and cortex; ROSMAP-derived DLPFC snRNA-seq; SEA-AD MTG; CELLxGENE fibroblast collections; GSE226189; GSE297234 (Lu et al.); figshare 10.6084/m9.figshare.30898748 (Sengstack et al.). Add: accession-level detail and per-cohort sample tables.

**Software.** Add: language, package versions, random seeds, entry-point scripts.

## Limitations

The reprogramming analysis rests on two donors in one experiment, so donor-level variance cannot be estimated and, in the same-platform test, age is confounded with individual. The transcription-factor screen uses neonatal foreskin fibroblasts aged by serial passaging in one cell line, and the passage of the perturbed cells is not recorded in the deposited matrix, so the late-passage reference was inferred; both screen-level results are reported under either choice. Cell-state assignments in the reprogramming data are cluster-level signature labels, not lineage tracing. Our distance statistic in the perturbation matrix carries a positive offset that precludes per-factor classification there. The rulers transfer as orderings, not as calibrated ages. Transcriptome distance is a molecular proxy and cannot adjudicate functional rejuvenation, in vivo results, or non-transcriptional mechanisms. Finally, our evidence is that approach is not detectable in these datasets at these sample sizes, which is weaker than evidence that approach does not occur.

## Data and code availability

All data are public. Add: repository URL, frozen ruler and anchor files, pre-registration documents, entry-point scripts for each analysis, and archived DOI.

## Author contributions and competing interests

[Add. The author declares no competing interests. Analyses were implemented with the assistance of AI coding agents operating under the written pre-registration protocols described in Methods; all pre-registrations, decision rules and interpretations are the author's.

## Acknowledgements

[Add.

## References

1. Lu, J. Y. *et al.* Prevalent mesenchymal drift in aging and disease is reversed by partial reprogramming. *Cell* **188**, 5895–5911 (2025).
2. Sengstack, J. *et al.* Systematic identification of single transcription factor perturbations that drive cellular and tissue rejuvenation. *Proc. Natl Acad. Sci. USA* **123**, e2515183123 (2026).
3. Southard, K. M. *et al.* Comprehensive transcription factor perturbations recapitulate fibroblast transcriptional states. *Nat. Genet.* **57**, 2323–2334 (2025).

[Add: GTEx, SEA-AD, Seurat, cell-cycle gene sets, epigenetic clock and prior reprogramming references cited in the text.
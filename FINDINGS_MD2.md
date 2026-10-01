# FINDINGS_MD2 — Lu et al. cell states, then both instruments on the same cells

**Status:** Task 1: `negative_holds`. Task 2: `instruments_disagree`. Task 3: aging signature not run (GSE113957 not DESeq2-rebuildable); MD and frozen ruler tabulated. Fired: Task 1 `negative_holds`; Task 2 `instruments_disagree`. Seed `20260914`. boot `20260918`. n_perm=200. n_boot=200. n_random=200. Frozen ruler `frozen_ruler_ridge_raw.npz` exists=True.

Does not modify any existing `FINDINGS_*.md` or `FALSIFICATION.md`. No GSE325735. Uncalibrated R² is reported and is never a gate. Nothing averaged across donors, clusters, or resolutions. Refit nothing. d7 is not the primary Stage 2 endpoint.

Reproduced by `src/md2_run.py`. Frozen ruler: `<repo>\results\fibro\frozen_ruler_ridge_raw.npz` exists=True.

Flag: `PREREG_TASK1.flag` exists=True. `PREREG_TASK2.flag` exists=True. `PREREG_TASK3.flag` exists=True. `CLUSTER_SETTINGS.flag` exists=True. `DECLARED_BEFORE_SCORES.flag` exists=True.

## Methods as published

Lu et al., Cell 2025, 188:5895–5911; DOI 10.1016/j.cell.2025.07.031. STAR Methods quotes:

> We obtained the gene lists from the hallmark gene sets: HALLMARK EPITHELIAL MESENCHYMAL TRANSITION and HALLMARK TGF BETA SIGNALING (https://www.gsea-msigdb.org/). We added SNAI1, ZEB1, ZEB2, TWIST1, and TWIST2 to HALLMARK EPITHELIAL MESENCHYMAL TRANSITION to compose the MD signature based on literature review (Table S3). We derived human fibroblast aging signatures using RNA-seq data from Fleischer et al. (2018) by determining significant differentially expressed genes (adjusted p-value < 0.05, log2(fold change) > 0.5) in old cells versus young cells as classified based on PCA separation using DESeq2 v1.40.2 (Table S3).

> Count matrices were loaded to create Seurat objects using Seurat v5.0.0. Cells with more than 500 genes and genes detected in at least 5 cells were kept for further processing. Separate samples were merged and normalized using sctransform v2, with regression for percent mitochondrial genes. Gene lists were used to calculate pathway score using the Seurat AddModuleScore function. Trajectory inference was performed using Slingshot v2.12.0, designating day 0 fibroblast populations as starting clusters.

> All the module expressions of the MD gene were based on log-scaled normalized UMI counts and calculated by the function AddModuleScore in Seurat package.

> To normalize sequencing reads for each gene, we utilized the "NormalizeData" function with the parameter "normalization.method = "LogNormalize". This involved dividing the UMI counts of each gene by the total UMIs of the cell, multiplying by the median of the total UMIs, and transforming the result using the natural logarithm.

> The partially reprogrammed populations showed a robust reversal of transcriptomic aging changes, while non-reprogrammed populations maintained the aging transcriptome (Figure 3G). This pattern was accompanied by downregulation of the MD and TGF-β pathway scores in the partially reprogramming cells and their maintained expression in non-reprogrammed cells (Figure 3G).

- MSigDB version used: `https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt`
- HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION n=200
- TFs added to EMT: ['SNAI1', 'ZEB1', 'ZEB2', 'TWIST1', 'TWIST2'] already in EMT: []
- MD built n=205
- HALLMARK_TGF_BETA_SIGNALING n=54
- mmc3 `MD_signatures` column `MD score` n=205 equal_to_built=True
- mmc3 only_built=[] only_supplement=[]
- mmc3 `TGFB score` n=54 equal_to_built=True
- mmc3 TGFB only_built=[] only_supplement=[]
- Mismatch is reported, not reconciled. Scoring uses the STAR Methods built set.
- AddModuleScore: Python implementation of Seurat's control-gene-bin algorithm (nbin=24, ctrl=100, quantile bins with 1e-30 Gaussian noise). LogNormalize scale.factor=10000 (Seurat NormalizeData default). Paper public-data prose multiplies by median total UMIs; that prose is not used.
- Clustering Louvain: STAR Methods specify Seurat v5 sctransform v2. This environment has no R/sctransform. Deviation (not a silent substitute): LogNormalize + quadratic HVG n=3000 + ScaleData-style percent.mt residualization + PCA 50 + SNN k.param=20 prune.SNN=1/15 + networkx louvain_communities resolution=0.8 (Seurat FindClusters default; STAR Methods omit resolution).
- Louvain n clusters GM00731=19 GM23815=19. Not retuned.
- Aging signature (Fleischer DESeq2): built=False reason=GSE113957_fpkm.txt.gz values are not integer counts (cannot apply DESeq2). FINDINGS_FIBRO2.md: series matrix HEAD 404. STAR Methods require DESeq2 old-vs-young with PCA-based old/young split. Not approximating from FPKM and not substituting mmc3 Age up (n=1533) / Age down (n=2007). Aging-signature comparison reported as not run.
- mmc3 Age up n=1533 Age down n=2007 used_as_substitute=False

## Failures recorded

- **addmodulescore:** GSE226189_MD:MD: 0 feature genes mapped of 205. Not scoring. Resolved on retry: geneCOUNT Tracking_ID is Ensembl; mapped to symbols via GSE297234 10x features + frozen ruler; MD then scored (see Task 3 table).

## Task 1 pre-registration (verbatim, written before any age score)

Written **before any age score**, 2026-09-19. All four resolutions are pre-registered. Not chosen after seeing which one gives a cleaner age trajectory.

Scope: aged donor GM00731, all timepoints concatenated, jointly clustered. Task 2 Louvain is repeated independently on GM23815 (never pooled across donors).

k = 3, 8, 15 (FIBRO3 setting, only k changes):
- Genes: frozen-ruler overlap genes with count-sum > 0 (FINDINGS_FIBRO3.md Task 1 / FINDINGS_FIBRO.md Stage 2 cluster spec).
- Transform: log1p(CP10k) on those genes. Library size = row sum of the kept overlap genes.
- Dimensionality: sklearn.decomposition.PCA(n_components=20, or n_cells−1 if smaller, svd_solver="randomized", random_state=20260914) on float32 log1p(CP10k). Cell order: libraries concatenated in day order 0, 3, 7, 10.
- Cluster: sklearn.cluster.KMeans(n_clusters=k, random_state=20260914, n_init=10).
- k=3 reuses the frozen labels in results/fibro3/t1_cell_labels.csv for GM00731 (the FIBRO3 setting). k=8 and k=15 are fit on the same PCA recipe, new KMeans.
- No n_genes>500 filter (FIBRO3 did not apply it). MIN_CELLS=20 to score a cluster×timepoint.

Louvain (their pipeline, STAR Methods of Lu et al. Cell 2025 188:5895–5911, “Single-cell RNA-seq analysis of newly generated data”):
- Cells with more than 500 genes; genes detected in at least 5 cells.
- STAR Methods specify Seurat v5.0.0, sctransform v2, regression on percent mitochondrial, Louvain, AddModuleScore. This environment has no R/Seurat/sctransform. Deviation, not a silent substitute: LogNormalize (Seurat NormalizeData default: ln(1 + 10000 * count / tot); paper public-data prose also describes multiplying by the median total UMIs — we use scale.factor=10000 and report the prose/default discrepancy), then ScaleData-style residualization of percent.mt (symbol upper startswith MT-) on the top 3000 highly variable genes (Seurat SCTransform default variable.features.n; STAR Methods do not state n HVGs), PCA 50 PCs (Seurat RunPCA default; STAR Methods do not state n PCs), shared-nearest-neighbour graph (k.param=20, prune.SNN=1/15, Seurat FindNeighbors defaults; STAR Methods do not state these), Louvain via networkx.community.louvain_communities, resolution=0.8 (Seurat FindClusters default; STAR Methods do not state the resolution). Seed 20260914.
- n clusters obtained is reported, not retuned.
- Slingshot from day-0 fibroblast clusters is optional and report-only; not required for scoring.

MIN_CELLS=20. Seed 20260914. Frozen ruler used as-is. Refit nothing.


Written **before any age score**, 2026-09-19. No threshold in this block is re-tuned after numbers exist.

Cluster the aged donor GM00731 jointly across timepoints at k = 3 (the FIBRO3 setting), 8, 15, and the Louvain pipeline in CLUSTER_SETTINGS. Fix all settings before any age score is computed.

For each resolution, per cluster × timepoint (≥20 cells): n cells, fraction of timepoint, frozen age score, MD score, TGF-β score, pluripotency−fibroblast score (frozen lists from FINDINGS_FIBRO.md), and the pre-registered 200-direction random null on d0→d7 and d0→d10 separately.

Frozen ruler: results/fibro/frozen_ruler_ridge_raw.npz used as-is. Age null = 200 permuted-ruler-weight directions (seed 20260914), empirical p = (n_random ≥ real + 1) / (n_random + 1). Endpoint decline = age_score(day 0) − age_score(day T). Positive = younger on the ruler. MD/TGF-β are Seurat AddModuleScore on LogNormalize counts (control-gene bins matched on average expression; nbin=24, ctrl=100). Do not average across donors, clusters, or resolutions. d7 is not the Stage 2 endpoint.

**Pre-registered reading (only the outcome that fired):**
- At finer resolution some cluster shows a d0→d7 or d0→d10 age decline with p ≤ 0.05 → FINDINGS_FIBRO3.md's per-cluster negative was a resolution artifact. Report which clusters, at which resolutions, and how stable the result is across resolutions. State plainly that this supersedes the FIBRO3 per-cluster reading, quoting the sentence it supersedes.
- No cluster at any resolution shows p ≤ 0.05 → the per-cluster negative holds and is not a resolution artifact. The frozen ruler does not move in these cells at any granularity tested.
- Results appear at k=15 but not at Louvain, or vice versa → report per resolution, do not average, and say the finding is resolution-dependent.


## Task 1 — resolution sweep tables

| resolution | cell_line | day | cluster | n_cells | frac_timepoint | below_min_cells | age_score | md_score | tgfb_score | pluri_primary | pluri_with_OSKM |
|---|---|---|---|---|---|---|---|---|---|---|---|
| k3 | GM00731 | 0 | 0 | 4941 | +0.984 | False | +4.224 | +0.508 | +0.049 | -1.391 | -1.391 |
| k3 | GM00731 | 0 | 1 | 76 | +0.015 | False | -3.165 | +0.343 | +0.101 | +0.566 | +0.566 |
| k3 | GM00731 | 0 | 2 | 4 | +0.001 | True | NA | NA | NA | NA | NA |
| k3 | GM00731 | 3 | 0 | 948 | +0.150 | False | +5.541 | +0.448 | +0.072 | +5.360 | +5.360 |
| k3 | GM00731 | 3 | 1 | 946 | +0.150 | False | -0.053 | +0.073 | -0.004 | +10.803 | +10.803 |
| k3 | GM00731 | 3 | 2 | 4424 | +0.700 | False | +0.639 | +0.095 | -0.057 | +8.503 | +8.503 |
| k3 | GM00731 | 7 | 0 | 2719 | +0.393 | False | +4.806 | +0.432 | +0.049 | +4.665 | +4.665 |
| k3 | GM00731 | 7 | 1 | 850 | +0.123 | False | -5.380 | +0.140 | +0.034 | +10.059 | +10.059 |
| k3 | GM00731 | 7 | 2 | 3346 | +0.484 | False | -2.904 | +0.106 | -0.033 | +9.832 | +9.832 |
| k3 | GM00731 | 10 | 0 | 3431 | +0.780 | False | +8.595 | +0.449 | -0.002 | +2.012 | +2.012 |
| k3 | GM00731 | 10 | 1 | 698 | +0.159 | False | +2.249 | +0.358 | +0.091 | +5.950 | +5.950 |
| k3 | GM00731 | 10 | 2 | 270 | +0.061 | False | +0.619 | +0.174 | -0.071 | +8.088 | +8.088 |
| k8 | GM00731 | 0 | 0 | 79 | +0.016 | False | +6.130 | +0.459 | +0.001 | -2.277 | -2.277 |
| k8 | GM00731 | 0 | 2 | 268 | +0.053 | False | +2.848 | +0.388 | +0.057 | +1.289 | +1.289 |
| k8 | GM00731 | 0 | 3 | 14 | +0.003 | True | NA | NA | NA | NA | NA |
| k8 | GM00731 | 0 | 4 | 4599 | +0.916 | False | +5.485 | +0.515 | +0.049 | -1.176 | -1.176 |
| k8 | GM00731 | 0 | 6 | 61 | +0.012 | False | -2.334 | +0.366 | +0.117 | +1.538 | +1.538 |
| k8 | GM00731 | 3 | 0 | 287 | +0.045 | False | +6.676 | +0.460 | +0.037 | +5.567 | +5.567 |
| k8 | GM00731 | 3 | 1 | 2433 | +0.385 | False | +2.173 | +0.040 | -0.096 | +9.163 | +9.163 |
| k8 | GM00731 | 3 | 2 | 2345 | +0.371 | False | +2.209 | +0.227 | +0.017 | +7.040 | +7.040 |
| k8 | GM00731 | 3 | 3 | 744 | +0.118 | False | +3.069 | +0.049 | -0.022 | +11.473 | +11.473 |
| k8 | GM00731 | 3 | 4 | 12 | +0.002 | True | NA | NA | NA | NA | NA |
| k8 | GM00731 | 3 | 5 | 138 | +0.022 | False | -1.424 | -0.012 | -0.089 | +9.031 | +9.031 |
| k8 | GM00731 | 3 | 6 | 224 | +0.035 | False | -0.770 | +0.156 | +0.050 | +7.613 | +7.613 |
| k8 | GM00731 | 3 | 7 | 135 | +0.021 | False | +7.651 | +0.550 | +0.104 | +3.473 | +3.473 |
| k8 | GM00731 | 7 | 0 | 1995 | +0.289 | False | +4.835 | +0.395 | +0.031 | +5.364 | +5.364 |
| k8 | GM00731 | 7 | 1 | 101 | +0.015 | False | +0.449 | +0.016 | -0.130 | +10.955 | +10.955 |
| k8 | GM00731 | 7 | 2 | 379 | +0.055 | False | +3.072 | +0.292 | +0.010 | +5.882 | +5.882 |
| k8 | GM00731 | 7 | 3 | 211 | +0.031 | False | -1.229 | +0.105 | +0.014 | +10.365 | +10.365 |
| k8 | GM00731 | 7 | 4 | 56 | +0.008 | False | +5.535 | +0.467 | +0.077 | +3.244 | +3.244 |
| k8 | GM00731 | 7 | 5 | 2753 | +0.398 | False | -2.693 | +0.084 | -0.033 | +10.243 | +10.243 |
| k8 | GM00731 | 7 | 6 | 638 | +0.092 | False | -4.046 | +0.154 | +0.041 | +9.451 | +9.451 |
| k8 | GM00731 | 7 | 7 | 782 | +0.113 | False | +7.326 | +0.473 | +0.072 | +4.067 | +4.067 |
| k8 | GM00731 | 10 | 0 | 439 | +0.100 | False | +9.145 | +0.454 | -0.020 | +2.124 | +2.124 |
| k8 | GM00731 | 10 | 1 | 14 | +0.003 | True | NA | NA | NA | NA | NA |
| k8 | GM00731 | 10 | 2 | 282 | +0.064 | False | +4.861 | +0.307 | -0.024 | +2.966 | +2.966 |
| k8 | GM00731 | 10 | 3 | 100 | +0.023 | False | +1.390 | +0.253 | +0.022 | +6.017 | +6.017 |
| k8 | GM00731 | 10 | 4 | 10 | +0.002 | True | NA | NA | NA | NA | NA |
| k8 | GM00731 | 10 | 5 | 89 | +0.020 | False | -1.005 | +0.063 | -0.088 | +10.500 | +10.500 |
| k8 | GM00731 | 10 | 6 | 604 | +0.137 | False | +3.430 | +0.371 | +0.102 | +5.973 | +5.973 |
| k8 | GM00731 | 10 | 7 | 2861 | +0.650 | False | +9.537 | +0.451 | -0.001 | +1.830 | +1.830 |
| k15 | GM00731 | 0 | 3 | 11 | +0.002 | True | NA | NA | NA | NA | NA |
| k15 | GM00731 | 0 | 6 | 14 | +0.003 | True | NA | NA | NA | NA | NA |
| k15 | GM00731 | 0 | 7 | 2276 | +0.453 | False | +7.841 | +0.524 | +0.029 | -1.391 | -1.391 |
| k15 | GM00731 | 0 | 8 | 235 | +0.047 | False | +4.057 | +0.383 | +0.055 | +1.722 | +1.722 |
| k15 | GM00731 | 0 | 10 | 60 | +0.012 | False | -1.904 | +0.363 | +0.116 | +2.286 | +2.286 |
| k15 | GM00731 | 0 | 12 | 30 | +0.006 | False | +3.565 | +0.407 | +0.044 | -1.665 | -1.665 |
| k15 | GM00731 | 0 | 13 | 2395 | +0.477 | False | +5.430 | +0.506 | +0.068 | -0.352 | -0.352 |
| k15 | GM00731 | 3 | 0 | 28 | +0.004 | False | +10.147 | +0.547 | +0.122 | +3.244 | +3.244 |
| k15 | GM00731 | 3 | 1 | 1 | +0.000 | True | NA | NA | NA | NA | NA |
| k15 | GM00731 | 3 | 2 | 127 | +0.020 | False | -3.029 | +0.026 | -0.008 | +9.434 | +9.434 |
| k15 | GM00731 | 3 | 3 | 132 | +0.021 | False | +6.537 | +0.368 | +0.014 | +5.614 | +5.614 |
| k15 | GM00731 | 3 | 4 | 1352 | +0.214 | False | +4.311 | +0.068 | -0.095 | +9.028 | +9.028 |
| k15 | GM00731 | 3 | 5 | 12 | +0.002 | True | NA | NA | NA | NA | NA |
| k15 | GM00731 | 3 | 6 | 708 | +0.112 | False | +3.770 | +0.044 | -0.022 | +11.506 | +11.506 |
| k15 | GM00731 | 3 | 7 | 5 | +0.001 | True | NA | NA | NA | NA | NA |
| k15 | GM00731 | 3 | 8 | 1639 | +0.259 | False | +5.227 | +0.296 | +0.048 | +6.724 | +6.724 |
| k15 | GM00731 | 3 | 9 | 430 | +0.068 | False | +3.190 | +0.045 | -0.050 | +10.261 | +10.261 |
| k15 | GM00731 | 3 | 10 | 88 | +0.014 | False | +3.949 | +0.363 | +0.137 | +5.552 | +5.552 |
| k15 | GM00731 | 3 | 11 | 166 | +0.026 | False | +8.479 | +0.537 | +0.054 | +5.668 | +5.668 |
| k15 | GM00731 | 3 | 12 | 123 | +0.019 | False | +7.402 | +0.458 | +0.098 | +4.324 | +4.324 |
| k15 | GM00731 | 3 | 13 | 7 | +0.001 | True | NA | NA | NA | NA | NA |
| k15 | GM00731 | 3 | 14 | 1500 | +0.237 | False | +0.366 | +0.022 | -0.092 | +8.807 | +8.807 |
| k15 | GM00731 | 7 | 0 | 192 | +0.028 | False | +9.114 | +0.443 | +0.074 | +4.375 | +4.375 |
| k15 | GM00731 | 7 | 1 | 891 | +0.129 | False | -2.638 | -0.039 | -0.106 | +12.224 | +12.224 |
| k15 | GM00731 | 7 | 2 | 437 | +0.063 | False | -5.817 | +0.025 | +0.014 | +11.321 | +11.321 |
| k15 | GM00731 | 7 | 3 | 1110 | +0.161 | False | +5.323 | +0.343 | +0.005 | +6.487 | +6.487 |
| k15 | GM00731 | 7 | 4 | 42 | +0.006 | False | +3.990 | +0.032 | -0.125 | +10.965 | +10.965 |
| k15 | GM00731 | 7 | 5 | 1746 | +0.252 | False | -0.743 | +0.137 | +0.002 | +9.047 | +9.047 |
| k15 | GM00731 | 7 | 6 | 206 | +0.030 | False | +0.660 | +0.111 | +0.014 | +9.592 | +9.592 |
| k15 | GM00731 | 7 | 7 | 8 | +0.001 | True | NA | NA | NA | NA | NA |
| k15 | GM00731 | 7 | 8 | 307 | +0.044 | False | +5.013 | +0.293 | +0.005 | +5.197 | +5.197 |
| k15 | GM00731 | 7 | 9 | 32 | +0.005 | False | -0.328 | +0.010 | -0.133 | +7.876 | +7.876 |
| k15 | GM00731 | 7 | 10 | 210 | +0.030 | False | +4.391 | +0.404 | +0.093 | +5.511 | +5.511 |
| k15 | GM00731 | 7 | 11 | 635 | +0.092 | False | +10.068 | +0.486 | +0.054 | +2.878 | +2.878 |
| k15 | GM00731 | 7 | 12 | 1059 | +0.153 | False | +5.129 | +0.423 | +0.065 | +5.089 | +5.089 |
| k15 | GM00731 | 7 | 13 | 6 | +0.001 | True | NA | NA | NA | NA | NA |
| k15 | GM00731 | 7 | 14 | 34 | +0.005 | False | -0.881 | +0.060 | -0.075 | +8.979 | +8.979 |
| k15 | GM00731 | 10 | 0 | 1942 | +0.441 | False | +9.864 | +0.443 | +0.011 | +2.233 | +2.233 |
| k15 | GM00731 | 10 | 1 | 47 | +0.011 | False | -1.106 | -0.067 | -0.113 | +12.481 | +12.481 |
| k15 | GM00731 | 10 | 2 | 83 | +0.019 | False | -4.519 | +0.068 | +0.033 | +10.531 | +10.531 |
| k15 | GM00731 | 10 | 3 | 301 | +0.068 | False | +9.435 | +0.414 | -0.030 | +3.059 | +3.059 |
| k15 | GM00731 | 10 | 4 | 2 | +0.000 | True | NA | NA | NA | NA | NA |
| k15 | GM00731 | 10 | 5 | 39 | +0.009 | False | +2.309 | +0.172 | -0.061 | +8.505 | +8.505 |
| k15 | GM00731 | 10 | 6 | 99 | +0.023 | False | +2.376 | +0.252 | +0.021 | +6.063 | +6.063 |
| k15 | GM00731 | 10 | 7 | 3 | +0.001 | True | NA | NA | NA | NA | NA |
| k15 | GM00731 | 10 | 8 | 259 | +0.059 | False | +6.397 | +0.305 | -0.029 | +2.634 | +2.634 |
| k15 | GM00731 | 10 | 9 | 5 | +0.001 | True | NA | NA | NA | NA | NA |
| k15 | GM00731 | 10 | 10 | 532 | +0.121 | False | +6.522 | +0.423 | +0.112 | +3.381 | +3.381 |
| k15 | GM00731 | 10 | 11 | 954 | +0.217 | False | +12.561 | +0.472 | -0.032 | +0.470 | +0.470 |
| k15 | GM00731 | 10 | 12 | 132 | +0.030 | False | +8.526 | +0.472 | +0.053 | +2.631 | +2.631 |
| k15 | GM00731 | 10 | 14 | 1 | +0.000 | True | NA | NA | NA | NA | NA |
| louvain | GM00731 | 0 | 0 | 3998 | +0.798 | False | +7.037 | +0.513 | +0.040 | -0.743 | -0.743 |
| louvain | GM00731 | 0 | 3 | 1 | +0.000 | True | NA | NA | NA | NA | NA |
| louvain | GM00731 | 0 | 4 | 9 | +0.002 | True | NA | NA | NA | NA | NA |
| louvain | GM00731 | 0 | 5 | 2 | +0.000 | True | NA | NA | NA | NA | NA |
| louvain | GM00731 | 0 | 10 | 56 | +0.011 | False | -1.425 | +0.356 | +0.117 | +1.860 | +1.860 |
| louvain | GM00731 | 0 | 12 | 569 | +0.114 | False | +7.865 | +0.510 | +0.049 | -1.198 | -1.198 |
| louvain | GM00731 | 0 | 15 | 372 | +0.074 | False | +4.533 | +0.428 | +0.050 | +1.003 | +1.003 |
| louvain | GM00731 | 0 | 16 | 2 | +0.000 | True | NA | NA | NA | NA | NA |
| louvain | GM00731 | 3 | 0 | 4 | +0.001 | True | NA | NA | NA | NA | NA |
| louvain | GM00731 | 3 | 3 | 1589 | +0.265 | False | +6.055 | +0.308 | +0.038 | +5.941 | +5.941 |
| louvain | GM00731 | 3 | 4 | 259 | +0.043 | False | +9.226 | +0.531 | +0.060 | +3.900 | +3.900 |
| louvain | GM00731 | 3 | 5 | 1523 | +0.254 | False | +4.013 | +0.060 | -0.096 | +9.042 | +9.042 |
| louvain | GM00731 | 3 | 7 | 108 | +0.018 | False | +5.193 | +0.163 | -0.047 | +7.227 | +7.227 |
| louvain | GM00731 | 3 | 8 | 973 | +0.162 | False | -0.379 | +0.020 | -0.098 | +8.942 | +8.942 |
| louvain | GM00731 | 3 | 10 | 89 | +0.015 | False | +7.215 | +0.316 | +0.117 | +5.974 | +5.974 |
| louvain | GM00731 | 3 | 11 | 557 | +0.093 | False | +3.280 | +0.080 | -0.089 | +7.831 | +7.831 |
| louvain | GM00731 | 3 | 13 | 442 | +0.074 | False | +4.836 | +0.079 | -0.031 | +9.860 | +9.860 |
| louvain | GM00731 | 3 | 14 | 5 | +0.001 | True | NA | NA | NA | NA | NA |
| louvain | GM00731 | 3 | 15 | 1 | +0.000 | True | NA | NA | NA | NA | NA |
| louvain | GM00731 | 3 | 16 | 302 | +0.050 | False | +0.353 | +0.009 | -0.015 | +12.562 | +12.562 |
| louvain | GM00731 | 3 | 18 | 136 | +0.023 | False | +3.692 | +0.171 | -0.016 | +7.768 | +7.768 |
| louvain | GM00731 | 7 | 0 | 9 | +0.001 | True | NA | NA | NA | NA | NA |
| louvain | GM00731 | 7 | 1 | 9 | +0.001 | True | NA | NA | NA | NA | NA |
| louvain | GM00731 | 7 | 2 | 2092 | +0.308 | False | +0.098 | +0.167 | -0.004 | +8.439 | +8.439 |
| louvain | GM00731 | 7 | 3 | 444 | +0.065 | False | +6.644 | +0.343 | +0.006 | +4.785 | +4.785 |
| louvain | GM00731 | 7 | 4 | 1732 | +0.255 | False | +8.658 | +0.456 | +0.045 | +3.152 | +3.152 |
| louvain | GM00731 | 7 | 5 | 48 | +0.007 | False | +0.508 | +0.052 | -0.088 | +11.007 | +11.007 |
| louvain | GM00731 | 7 | 6 | 1 | +0.000 | True | NA | NA | NA | NA | NA |
| louvain | GM00731 | 7 | 7 | 884 | +0.130 | False | +4.745 | +0.299 | +0.001 | +7.337 | +7.337 |
| louvain | GM00731 | 7 | 8 | 13 | +0.002 | True | NA | NA | NA | NA | NA |
| louvain | GM00731 | 7 | 9 | 909 | +0.134 | False | -2.000 | -0.024 | -0.106 | +11.825 | +11.825 |
| louvain | GM00731 | 7 | 10 | 191 | +0.028 | False | +3.966 | +0.347 | +0.079 | +5.835 | +5.835 |
| louvain | GM00731 | 7 | 11 | 64 | +0.009 | False | +4.110 | +0.117 | -0.064 | +6.920 | +6.920 |
| louvain | GM00731 | 7 | 12 | 1 | +0.000 | True | NA | NA | NA | NA | NA |
| louvain | GM00731 | 7 | 13 | 3 | +0.000 | True | NA | NA | NA | NA | NA |
| louvain | GM00731 | 7 | 14 | 354 | +0.052 | False | -5.699 | +0.011 | +0.005 | +11.693 | +11.693 |
| louvain | GM00731 | 7 | 15 | 3 | +0.000 | True | NA | NA | NA | NA | NA |
| louvain | GM00731 | 7 | 16 | 32 | +0.005 | False | +1.891 | +0.188 | +0.033 | +4.024 | +4.024 |
| louvain | GM00731 | 7 | 18 | 6 | +0.001 | True | NA | NA | NA | NA | NA |
| louvain | GM00731 | 10 | 0 | 1 | +0.000 | True | NA | NA | NA | NA | NA |
| louvain | GM00731 | 10 | 1 | 2313 | +0.533 | False | +11.686 | +0.446 | -0.018 | +1.788 | +1.788 |
| louvain | GM00731 | 10 | 2 | 4 | +0.001 | True | NA | NA | NA | NA | NA |
| louvain | GM00731 | 10 | 3 | 5 | +0.001 | True | NA | NA | NA | NA | NA |
| louvain | GM00731 | 10 | 4 | 31 | +0.007 | False | +13.645 | +0.408 | -0.045 | +2.259 | +2.259 |
| louvain | GM00731 | 10 | 5 | 21 | +0.005 | False | +2.744 | +0.320 | +0.008 | +6.696 | +6.696 |
| louvain | GM00731 | 10 | 6 | 1029 | +0.237 | False | +9.102 | +0.434 | -0.009 | +2.984 | +2.984 |
| louvain | GM00731 | 10 | 7 | 22 | +0.005 | False | +6.501 | +0.289 | -0.071 | +6.230 | +6.230 |
| louvain | GM00731 | 10 | 9 | 52 | +0.012 | False | -0.290 | -0.040 | -0.109 | +12.005 | +12.005 |
| louvain | GM00731 | 10 | 10 | 527 | +0.121 | False | +7.934 | +0.417 | +0.102 | +3.313 | +3.313 |
| louvain | GM00731 | 10 | 11 | 5 | +0.001 | True | NA | NA | NA | NA | NA |
| louvain | GM00731 | 10 | 14 | 80 | +0.018 | False | -2.372 | +0.069 | +0.019 | +10.355 | +10.355 |
| louvain | GM00731 | 10 | 16 | 7 | +0.002 | True | NA | NA | NA | NA | NA |
| louvain | GM00731 | 10 | 17 | 243 | +0.056 | False | +8.127 | +0.356 | -0.036 | +1.015 | +1.015 |

## Task 1 — frozen-ruler random-direction null (d0→d7 and d0→d10; not averaged)

| resolution | cluster | endpoint | n_cells_d0 | n_cells_end | score_d0 | score_end | decline | p | n_random_ge_real | n_random | pass_p | ok | reason | null_kind |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| k3 | 0 | d0→d7 | +4941.000 | +2719.000 | +4.224 | +4.806 | -0.582 | +0.652 | +130.000 | 200 | False | True | NA | permuted_ruler_weights_reused_from_FIBRO3 |
| k3 | 0 | d0→d10 | +4941.000 | +3431.000 | +4.224 | +8.595 | -4.371 | +1.000 | +200.000 | 200 | False | True | NA | permuted_ruler_weights_reused_from_FIBRO3 |
| k3 | 1 | d0→d7 | +76.000 | +850.000 | -3.165 | -5.380 | +2.215 | +0.144 | +28.000 | 200 | False | True | NA | permuted_ruler_weights_reused_from_FIBRO3 |
| k3 | 1 | d0→d10 | +76.000 | +698.000 | -3.165 | +2.249 | -5.414 | +0.881 | +176.000 | 200 | False | True | NA | permuted_ruler_weights_reused_from_FIBRO3 |
| k3 | 2 | d0→d7 | +4.000 | +3346.000 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights_reused_from_FIBRO3 |
| k3 | 2 | d0→d10 | +4.000 | +270.000 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights_reused_from_FIBRO3 |
| k8 | 0 | d0→d7 | +79.000 | +1995.000 | +6.130 | +4.835 | +1.295 | +0.139 | +27.000 | 200 | False | True | NA | permuted_ruler_weights |
| k8 | 0 | d0→d10 | +79.000 | +439.000 | +6.130 | +9.145 | -3.015 | +0.771 | +154.000 | 200 | False | True | NA | permuted_ruler_weights |
| k8 | 1 | d0→d7 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| k8 | 1 | d0→d10 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| k8 | 2 | d0→d7 | +268.000 | +379.000 | +2.848 | +3.072 | -0.224 | +0.488 | +97.000 | 200 | False | True | NA | permuted_ruler_weights |
| k8 | 2 | d0→d10 | +268.000 | +282.000 | +2.848 | +4.861 | -2.013 | +0.925 | +185.000 | 200 | False | True | NA | permuted_ruler_weights |
| k8 | 3 | d0→d7 | +14.000 | +211.000 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| k8 | 3 | d0→d10 | +14.000 | +100.000 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| k8 | 4 | d0→d7 | +4599.000 | +56.000 | +5.485 | +5.535 | -0.050 | +0.657 | +131.000 | 200 | False | True | NA | permuted_ruler_weights |
| k8 | 4 | d0→d10 | +4599.000 | +10.000 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| k8 | 5 | d0→d7 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| k8 | 5 | d0→d10 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| k8 | 6 | d0→d7 | +61.000 | +638.000 | -2.334 | -4.046 | +1.712 | +0.100 | +19.000 | 200 | False | True | NA | permuted_ruler_weights |
| k8 | 6 | d0→d10 | +61.000 | +604.000 | -2.334 | +3.430 | -5.764 | +0.920 | +184.000 | 200 | False | True | NA | permuted_ruler_weights |
| k8 | 7 | d0→d7 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| k8 | 7 | d0→d10 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| k15 | 0 | d0→d7 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| k15 | 0 | d0→d10 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| k15 | 1 | d0→d7 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| k15 | 1 | d0→d10 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| k15 | 2 | d0→d7 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| k15 | 2 | d0→d10 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| k15 | 3 | d0→d7 | +11.000 | +1110.000 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| k15 | 3 | d0→d10 | +11.000 | +301.000 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| k15 | 4 | d0→d7 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| k15 | 4 | d0→d10 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| k15 | 5 | d0→d7 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| k15 | 5 | d0→d10 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| k15 | 6 | d0→d7 | +14.000 | +206.000 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| k15 | 6 | d0→d10 | +14.000 | +99.000 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| k15 | 7 | d0→d7 | +2276.000 | +8.000 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| k15 | 7 | d0→d10 | +2276.000 | +3.000 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| k15 | 8 | d0→d7 | +235.000 | +307.000 | +4.057 | +5.013 | -0.956 | +0.796 | +159.000 | 200 | False | True | NA | permuted_ruler_weights |
| k15 | 8 | d0→d10 | +235.000 | +259.000 | +4.057 | +6.397 | -2.340 | +0.970 | +194.000 | 200 | False | True | NA | permuted_ruler_weights |
| k15 | 9 | d0→d7 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| k15 | 9 | d0→d10 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| k15 | 10 | d0→d7 | +60.000 | +210.000 | -1.904 | +4.391 | -6.295 | +0.980 | +196.000 | 200 | False | True | NA | permuted_ruler_weights |
| k15 | 10 | d0→d10 | +60.000 | +532.000 | -1.904 | +6.522 | -8.426 | +0.995 | +199.000 | 200 | False | True | NA | permuted_ruler_weights |
| k15 | 11 | d0→d7 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| k15 | 11 | d0→d10 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| k15 | 12 | d0→d7 | +30.000 | +1059.000 | +3.565 | +5.129 | -1.565 | +0.577 | +115.000 | 200 | False | True | NA | permuted_ruler_weights |
| k15 | 12 | d0→d10 | +30.000 | +132.000 | +3.565 | +8.526 | -4.962 | +0.950 | +190.000 | 200 | False | True | NA | permuted_ruler_weights |
| k15 | 13 | d0→d7 | +2395.000 | +6.000 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| k15 | 13 | d0→d10 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| k15 | 14 | d0→d7 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| k15 | 14 | d0→d10 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| louvain | 0 | d0→d7 | +3998.000 | +9.000 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| louvain | 0 | d0→d10 | +3998.000 | +1.000 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| louvain | 1 | d0→d7 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| louvain | 1 | d0→d10 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| louvain | 2 | d0→d7 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| louvain | 2 | d0→d10 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| louvain | 3 | d0→d7 | +1.000 | +444.000 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| louvain | 3 | d0→d10 | +1.000 | +5.000 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| louvain | 4 | d0→d7 | +9.000 | +1732.000 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| louvain | 4 | d0→d10 | +9.000 | +31.000 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| louvain | 5 | d0→d7 | +2.000 | +48.000 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| louvain | 5 | d0→d10 | +2.000 | +21.000 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| louvain | 6 | d0→d7 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| louvain | 6 | d0→d10 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| louvain | 7 | d0→d7 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| louvain | 7 | d0→d10 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| louvain | 8 | d0→d7 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| louvain | 8 | d0→d10 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| louvain | 9 | d0→d7 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| louvain | 9 | d0→d10 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| louvain | 10 | d0→d7 | +56.000 | +191.000 | -1.425 | +3.966 | -5.391 | +0.930 | +186.000 | 200 | False | True | NA | permuted_ruler_weights |
| louvain | 10 | d0→d10 | +56.000 | +527.000 | -1.425 | +7.934 | -9.359 | +1.000 | +200.000 | 200 | False | True | NA | permuted_ruler_weights |
| louvain | 11 | d0→d7 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| louvain | 11 | d0→d10 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| louvain | 12 | d0→d7 | +569.000 | +1.000 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| louvain | 12 | d0→d10 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| louvain | 13 | d0→d7 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| louvain | 13 | d0→d10 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| louvain | 14 | d0→d7 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| louvain | 14 | d0→d10 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| louvain | 15 | d0→d7 | +372.000 | +3.000 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| louvain | 15 | d0→d10 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| louvain | 16 | d0→d7 | +2.000 | +32.000 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| louvain | 16 | d0→d10 | +2.000 | +7.000 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | permuted_ruler_weights |
| louvain | 17 | d0→d7 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| louvain | 17 | d0→d10 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| louvain | 18 | d0→d7 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |
| louvain | 18 | d0→d10 | NA | NA | NA | NA | NA | NA | NA | 200 | NA | False | missing cluster row at an endpoint | permuted_ruler_weights |

k=3 age scores and nulls reused from `results/fibro3/` (not recomputed; TMM among both donors in FIBRO3). k=8, k=15, Louvain TMM among GM00731 cluster×timepoint rows of that resolution only.

## Task 1 reading

Fired key: `negative_holds`.

No cluster at any resolution shows p ≤ 0.05. The per-cluster negative holds and is not a resolution artifact. The frozen ruler does not move in these cells at any granularity tested.

Louvain occupancy: most Louvain clusters are occupied at one timepoint. Aged GM00731 tests with ≥20 cells at both endpoints (from t1_reading.json): k3 n_ok=4 n_pass=0; k8 n_ok=7 n_pass=0; k15 n_ok=6 n_pass=0; louvain n_ok=2 n_pass=0.

d7 is not the Stage 2 endpoint. Rows are not averaged across clusters or resolutions.

## Task 2 pre-registration (verbatim, written before any MD vs ruler comparison)

Written **before any MD vs frozen-ruler comparison**, 2026-09-19. No threshold in this block is re-tuned after numbers exist.

On the Louvain clustering (their pipeline), per cluster × timepoint, both donors reported separately:

1. Reproduce their Figure 3G result. Label each Louvain cluster by the mmc3.xlsx sheet Reprog_cell_state_signatures column with the highest mean AddModuleScore (Fibroblast, PartialReprog, EarlyPluripotency, Pluripotency, NonReprog). Reproduction matches their reported direction iff (a) at least one PartialReprog-assigned cluster and one NonReprog-assigned cluster exist with ≥20 cells in the aged donor, and (b) mean per-cell MD score in PartialReprog-assigned clusters is lower than in NonReprog-assigned clusters (aged donor). If it does not match, the comparison below is not valid — say so and stop Task 2.
2. Score both instruments on the same cells: frozen GTEx fibroblast ruler and MD score (and their aging signature if buildable).
3. Cross-correlate. Across cluster × timepoint rows with ≥20 cells, Spearman ρ between frozen age score and MD score, with permutation null (shuffle MD, hold age; n_perm=200, seed 20260914) and bootstrap CI (B=200, seed 20260918). Donors reported separately; not averaged.
4. Composition test on their metric. Per cluster, d0→d7 and d0→d10 MD decline with a 200-draw size-matched random gene-set null (same n genes, matched on mean-expression bin; AddModuleScore; seed 20260914). Not permuted ruler weights. State the difference.

**Pre-registered reading (only the outcome that fired):**
- MD moves per cluster and the frozen ruler does too (Task 1 outcome 1) → both instruments agree; their bend claim is supported and our earlier negative was our clustering. Report it that way, without hedging.
- MD moves per cluster, the frozen ruler does not → the instruments disagree on the same cells. Describe both, name what would settle it, and declare no winner.
- Neither moves per cluster → our reproduction does not recover their published result; report the reproduction failure as the finding and do not present it as a refutation.
- MD moves but fails its size-matched random-gene-set null → their metric's movement is not specific to mesenchymal genes in this dataset. Report the null value; this is a strong claim and needs the number stated next to it.

“MD moves per cluster” = at least one aged-donor Louvain cluster with ≥20 cells at both endpoints shows MD decline d0→d7 or d0→d10 with p ≤ 0.05 on the size-matched gene-set null.


## Task 2 — Figure 3G reproduction check

matches_reported_direction=True
n_PartialReprog=5832 n_NonReprog=5459
md_mean_PartialReprog=+0.150 md_mean_NonReprog=+0.439 delta=-0.289
n_cells_by_label={'EarlyPluripotency': 4414, 'Fibroblast': 5988, 'NonReprog': 5459, 'PartialReprog': 5832, 'Pluripotency': 439}
criterion=(a) ≥1 PartialReprog-assigned cluster and ≥1 NonReprog-assigned cluster with ≥20 cells in the aged donor; (b) mean per-cell MD in PartialReprog-assigned clusters < mean per-cell MD in NonReprog-assigned clusters.

Cluster state labels (argmax mean AddModuleScore of mmc3 Reprog_cell_state_signatures):

| cell_line | cluster | n_cells | label | mean_Fibroblast | mean_PartialReprog | mean_EarlyPluripotency | mean_Pluripotency | mean_NonReprog |
|---|---|---|---|---|---|---|---|---|
| GM00731 | 0 | 4012 | Fibroblast | +0.815 | +0.105 | +0.372 | -0.143 | +0.197 |
| GM00731 | 1 | 2322 | NonReprog | +0.451 | +0.121 | +0.266 | -0.154 | +0.898 |
| GM00731 | 2 | 2096 | EarlyPluripotency | +0.296 | +0.266 | +0.840 | +0.122 | -0.018 |
| GM00731 | 3 | 2039 | PartialReprog | +0.362 | +0.508 | +0.295 | -0.153 | +0.261 |
| GM00731 | 4 | 2031 | NonReprog | +0.500 | +0.157 | +0.342 | -0.122 | +0.578 |
| GM00731 | 5 | 1594 | PartialReprog | +0.187 | +0.618 | +0.407 | -0.074 | -0.030 |
| GM00731 | 6 | 1030 | Fibroblast | +0.434 | +0.259 | +0.379 | -0.134 | +0.424 |
| GM00731 | 7 | 1014 | EarlyPluripotency | +0.360 | +0.300 | +0.630 | +0.045 | +0.283 |
| GM00731 | 8 | 986 | PartialReprog | +0.130 | +0.665 | +0.487 | -0.009 | -0.104 |
| GM00731 | 9 | 961 | EarlyPluripotency | +0.106 | +0.400 | +0.845 | +0.474 | -0.116 |
| GM00731 | 10 | 863 | NonReprog | +0.299 | -0.115 | +0.108 | -0.039 | +0.512 |
| GM00731 | 11 | 626 | PartialReprog | +0.163 | +0.563 | +0.390 | -0.080 | -0.045 |
| GM00731 | 12 | 570 | Fibroblast | +0.872 | +0.143 | +0.426 | -0.147 | +0.158 |
| GM00731 | 13 | 445 | PartialReprog | +0.116 | +0.415 | +0.252 | -0.099 | -0.054 |
| GM00731 | 14 | 439 | Pluripotency | -0.001 | -0.052 | +0.591 | +0.659 | -0.069 |
| GM00731 | 15 | 376 | Fibroblast | +0.653 | +0.353 | +0.288 | -0.170 | +0.121 |
| GM00731 | 16 | 343 | EarlyPluripotency | +0.075 | +0.120 | +0.208 | +0.006 | +0.019 |
| GM00731 | 17 | 243 | NonReprog | +0.334 | +0.370 | +0.253 | -0.170 | +0.673 |
| GM00731 | 18 | 142 | PartialReprog | +0.213 | +0.526 | +0.322 | -0.118 | +0.044 |

| cell_line | cluster | n_cells | label | mean_Fibroblast | mean_PartialReprog | mean_EarlyPluripotency | mean_Pluripotency | mean_NonReprog |
|---|---|---|---|---|---|---|---|---|
| GM23815 | 0 | 6489 | Fibroblast | +0.730 | +0.103 | +0.381 | -0.198 | +0.314 |
| GM23815 | 1 | 3251 | EarlyPluripotency | +0.341 | +0.241 | +0.730 | -0.008 | +0.089 |
| GM23815 | 2 | 2731 | EarlyPluripotency | +0.268 | +0.423 | +0.685 | +0.013 | +0.051 |
| GM23815 | 3 | 2570 | EarlyPluripotency | +0.183 | +0.344 | +0.823 | +0.304 | -0.043 |
| GM23815 | 4 | 2242 | Fibroblast | +0.470 | +0.152 | +0.445 | -0.139 | +0.354 |
| GM23815 | 5 | 2171 | PartialReprog | +0.327 | +0.554 | +0.310 | -0.222 | +0.199 |
| GM23815 | 6 | 2149 | EarlyPluripotency | +0.091 | +0.495 | +0.694 | +0.334 | -0.088 |
| GM23815 | 7 | 2069 | PartialReprog | +0.146 | +0.638 | +0.526 | -0.021 | -0.049 |
| GM23815 | 8 | 1812 | NonReprog | +0.484 | +0.123 | +0.349 | -0.178 | +0.742 |
| GM23815 | 9 | 1365 | Fibroblast | +0.439 | +0.286 | +0.432 | -0.188 | +0.309 |
| GM23815 | 10 | 1249 | Pluripotency | -0.023 | -0.032 | +0.387 | +0.521 | -0.057 |
| GM23815 | 11 | 1140 | PartialReprog | +0.140 | +0.500 | +0.351 | -0.101 | -0.020 |
| GM23815 | 12 | 1101 | PartialReprog | +0.226 | +0.528 | +0.407 | -0.162 | +0.051 |
| GM23815 | 13 | 898 | EarlyPluripotency | +0.160 | -0.017 | +0.281 | +0.057 | +0.163 |
| GM23815 | 14 | 896 | Fibroblast | +0.604 | +0.337 | +0.304 | -0.225 | +0.256 |
| GM23815 | 15 | 761 | EarlyPluripotency | +0.366 | +0.395 | +0.628 | -0.048 | +0.263 |
| GM23815 | 16 | 363 | Fibroblast | +0.741 | +0.099 | +0.383 | -0.182 | +0.320 |
| GM23815 | 17 | 307 | EarlyPluripotency | +0.285 | +0.331 | +0.567 | -0.007 | +0.287 |
| GM23815 | 18 | 294 | EarlyPluripotency | +0.177 | +0.253 | +0.797 | +0.283 | -0.022 |

## Task 2 — both instruments on the same Louvain cluster × timepoint rows

| resolution | cell_line | day | cluster | label | n_cells | below_min_cells | age_score | md_score | tgfb_score | pluri_primary |
|---|---|---|---|---|---|---|---|---|---|---|
| louvain | GM00731 | 0 | 0 | Fibroblast | 3998 | False | +7.037 | +0.513 | +0.040 | -0.743 |
| louvain | GM00731 | 0 | 3 | PartialReprog | 1 | True | NA | NA | NA | NA |
| louvain | GM00731 | 0 | 4 | NonReprog | 9 | True | NA | NA | NA | NA |
| louvain | GM00731 | 0 | 5 | PartialReprog | 2 | True | NA | NA | NA | NA |
| louvain | GM00731 | 0 | 10 | NonReprog | 56 | False | -1.425 | +0.356 | +0.117 | +1.860 |
| louvain | GM00731 | 0 | 12 | Fibroblast | 569 | False | +7.865 | +0.510 | +0.049 | -1.198 |
| louvain | GM00731 | 0 | 15 | Fibroblast | 372 | False | +4.533 | +0.428 | +0.050 | +1.003 |
| louvain | GM00731 | 0 | 16 | EarlyPluripotency | 2 | True | NA | NA | NA | NA |
| louvain | GM00731 | 3 | 0 | Fibroblast | 4 | True | NA | NA | NA | NA |
| louvain | GM00731 | 3 | 3 | PartialReprog | 1589 | False | +6.055 | +0.308 | +0.038 | +5.941 |
| louvain | GM00731 | 3 | 4 | NonReprog | 259 | False | +9.226 | +0.531 | +0.060 | +3.900 |
| louvain | GM00731 | 3 | 5 | PartialReprog | 1523 | False | +4.013 | +0.060 | -0.096 | +9.042 |
| louvain | GM00731 | 3 | 7 | EarlyPluripotency | 108 | False | +5.193 | +0.163 | -0.047 | +7.227 |
| louvain | GM00731 | 3 | 8 | PartialReprog | 973 | False | -0.379 | +0.020 | -0.098 | +8.942 |
| louvain | GM00731 | 3 | 10 | NonReprog | 89 | False | +7.215 | +0.316 | +0.117 | +5.974 |
| louvain | GM00731 | 3 | 11 | PartialReprog | 557 | False | +3.280 | +0.080 | -0.089 | +7.831 |
| louvain | GM00731 | 3 | 13 | PartialReprog | 442 | False | +4.836 | +0.079 | -0.031 | +9.860 |
| louvain | GM00731 | 3 | 14 | Pluripotency | 5 | True | NA | NA | NA | NA |
| louvain | GM00731 | 3 | 15 | Fibroblast | 1 | True | NA | NA | NA | NA |
| louvain | GM00731 | 3 | 16 | EarlyPluripotency | 302 | False | +0.353 | +0.009 | -0.015 | +12.562 |
| louvain | GM00731 | 3 | 18 | PartialReprog | 136 | False | +3.692 | +0.171 | -0.016 | +7.768 |
| louvain | GM00731 | 7 | 0 | Fibroblast | 9 | True | NA | NA | NA | NA |
| louvain | GM00731 | 7 | 1 | NonReprog | 9 | True | NA | NA | NA | NA |
| louvain | GM00731 | 7 | 2 | EarlyPluripotency | 2092 | False | +0.098 | +0.167 | -0.004 | +8.439 |
| louvain | GM00731 | 7 | 3 | PartialReprog | 444 | False | +6.644 | +0.343 | +0.006 | +4.785 |
| louvain | GM00731 | 7 | 4 | NonReprog | 1732 | False | +8.658 | +0.456 | +0.045 | +3.152 |
| louvain | GM00731 | 7 | 5 | PartialReprog | 48 | False | +0.508 | +0.052 | -0.088 | +11.007 |
| louvain | GM00731 | 7 | 6 | Fibroblast | 1 | True | NA | NA | NA | NA |
| louvain | GM00731 | 7 | 7 | EarlyPluripotency | 884 | False | +4.745 | +0.299 | +0.001 | +7.337 |
| louvain | GM00731 | 7 | 8 | PartialReprog | 13 | True | NA | NA | NA | NA |
| louvain | GM00731 | 7 | 9 | EarlyPluripotency | 909 | False | -2.000 | -0.024 | -0.106 | +11.825 |
| louvain | GM00731 | 7 | 10 | NonReprog | 191 | False | +3.966 | +0.347 | +0.079 | +5.835 |
| louvain | GM00731 | 7 | 11 | PartialReprog | 64 | False | +4.110 | +0.117 | -0.064 | +6.920 |
| louvain | GM00731 | 7 | 12 | Fibroblast | 1 | True | NA | NA | NA | NA |
| louvain | GM00731 | 7 | 13 | PartialReprog | 3 | True | NA | NA | NA | NA |
| louvain | GM00731 | 7 | 14 | Pluripotency | 354 | False | -5.699 | +0.011 | +0.005 | +11.693 |
| louvain | GM00731 | 7 | 15 | Fibroblast | 3 | True | NA | NA | NA | NA |
| louvain | GM00731 | 7 | 16 | EarlyPluripotency | 32 | False | +1.891 | +0.188 | +0.033 | +4.024 |
| louvain | GM00731 | 7 | 18 | PartialReprog | 6 | True | NA | NA | NA | NA |
| louvain | GM00731 | 10 | 0 | Fibroblast | 1 | True | NA | NA | NA | NA |
| louvain | GM00731 | 10 | 1 | NonReprog | 2313 | False | +11.686 | +0.446 | -0.018 | +1.788 |
| louvain | GM00731 | 10 | 2 | EarlyPluripotency | 4 | True | NA | NA | NA | NA |
| louvain | GM00731 | 10 | 3 | PartialReprog | 5 | True | NA | NA | NA | NA |
| louvain | GM00731 | 10 | 4 | NonReprog | 31 | False | +13.645 | +0.408 | -0.045 | +2.259 |
| louvain | GM00731 | 10 | 5 | PartialReprog | 21 | False | +2.744 | +0.320 | +0.008 | +6.696 |
| louvain | GM00731 | 10 | 6 | Fibroblast | 1029 | False | +9.102 | +0.434 | -0.009 | +2.984 |
| louvain | GM00731 | 10 | 7 | EarlyPluripotency | 22 | False | +6.501 | +0.289 | -0.071 | +6.230 |
| louvain | GM00731 | 10 | 9 | EarlyPluripotency | 52 | False | -0.290 | -0.040 | -0.109 | +12.005 |
| louvain | GM00731 | 10 | 10 | NonReprog | 527 | False | +7.934 | +0.417 | +0.102 | +3.313 |
| louvain | GM00731 | 10 | 11 | PartialReprog | 5 | True | NA | NA | NA | NA |
| louvain | GM00731 | 10 | 14 | Pluripotency | 80 | False | -2.372 | +0.069 | +0.019 | +10.355 |
| louvain | GM00731 | 10 | 16 | EarlyPluripotency | 7 | True | NA | NA | NA | NA |
| louvain | GM00731 | 10 | 17 | NonReprog | 243 | False | +8.127 | +0.356 | -0.036 | +1.015 |
| louvain | GM23815 | 0 | 0 | Fibroblast | 6454 | False | +3.253 | +0.483 | +0.033 | -0.154 |
| louvain | GM23815 | 0 | 4 | Fibroblast | 2 | True | NA | NA | NA | NA |
| louvain | GM23815 | 0 | 5 | PartialReprog | 9 | True | NA | NA | NA | NA |
| louvain | GM23815 | 0 | 13 | EarlyPluripotency | 56 | False | -2.160 | +0.379 | +0.077 | +1.778 |
| louvain | GM23815 | 0 | 14 | Fibroblast | 886 | False | +1.849 | +0.420 | +0.016 | +0.918 |
| louvain | GM23815 | 0 | 16 | Fibroblast | 361 | False | +3.401 | +0.547 | -0.005 | -0.696 |
| louvain | GM23815 | 3 | 0 | Fibroblast | 19 | True | NA | NA | NA | NA |
| louvain | GM23815 | 3 | 1 | EarlyPluripotency | 3 | True | NA | NA | NA | NA |
| louvain | GM23815 | 3 | 3 | EarlyPluripotency | 3 | True | NA | NA | NA | NA |
| louvain | GM23815 | 3 | 4 | Fibroblast | 21 | False | +8.218 | +0.391 | +0.045 | +3.558 |
| louvain | GM23815 | 3 | 5 | PartialReprog | 2039 | False | +1.227 | +0.219 | +0.029 | +6.379 |
| louvain | GM23815 | 3 | 7 | PartialReprog | 2039 | False | -3.618 | -0.011 | -0.084 | +9.994 |
| louvain | GM23815 | 3 | 8 | NonReprog | 10 | True | NA | NA | NA | NA |
| louvain | GM23815 | 3 | 10 | Pluripotency | 91 | False | -8.931 | +0.028 | -0.012 | +10.539 |
| louvain | GM23815 | 3 | 11 | PartialReprog | 1065 | False | -0.322 | +0.007 | -0.065 | +10.360 |
| louvain | GM23815 | 3 | 12 | PartialReprog | 970 | False | +0.720 | +0.095 | -0.043 | +7.424 |
| louvain | GM23815 | 3 | 13 | EarlyPluripotency | 246 | False | +2.602 | +0.034 | +0.023 | +9.709 |
| louvain | GM23815 | 3 | 14 | Fibroblast | 5 | True | NA | NA | NA | NA |
| louvain | GM23815 | 3 | 16 | Fibroblast | 1 | True | NA | NA | NA | NA |
| louvain | GM23815 | 3 | 17 | EarlyPluripotency | 64 | False | -1.597 | +0.156 | -0.011 | +7.501 |
| louvain | GM23815 | 7 | 0 | Fibroblast | 14 | True | NA | NA | NA | NA |
| louvain | GM23815 | 7 | 1 | EarlyPluripotency | 3240 | False | -3.860 | +0.207 | -0.041 | +7.946 |
| louvain | GM23815 | 7 | 2 | EarlyPluripotency | 1 | True | NA | NA | NA | NA |
| louvain | GM23815 | 7 | 3 | EarlyPluripotency | 2537 | False | -6.390 | +0.003 | -0.069 | +12.127 |
| louvain | GM23815 | 7 | 4 | Fibroblast | 2210 | False | +1.091 | +0.410 | -0.015 | +4.365 |
| louvain | GM23815 | 7 | 5 | PartialReprog | 118 | False | +1.778 | +0.314 | -0.018 | +4.847 |
| louvain | GM23815 | 7 | 6 | EarlyPluripotency | 1819 | False | -5.449 | -0.048 | -0.083 | +12.685 |
| louvain | GM23815 | 7 | 7 | PartialReprog | 20 | False | -6.363 | +0.006 | -0.107 | +11.332 |
| louvain | GM23815 | 7 | 8 | NonReprog | 8 | True | NA | NA | NA | NA |
| louvain | GM23815 | 7 | 9 | Fibroblast | 3 | True | NA | NA | NA | NA |
| louvain | GM23815 | 7 | 10 | Pluripotency | 977 | False | -8.857 | -0.059 | +0.006 | +12.656 |
| louvain | GM23815 | 7 | 11 | PartialReprog | 65 | False | -2.332 | +0.023 | -0.101 | +10.838 |
| louvain | GM23815 | 7 | 12 | PartialReprog | 106 | False | +0.145 | +0.143 | -0.067 | +7.850 |
| louvain | GM23815 | 7 | 13 | EarlyPluripotency | 164 | False | -7.076 | +0.159 | +0.056 | +8.758 |
| louvain | GM23815 | 7 | 14 | Fibroblast | 5 | True | NA | NA | NA | NA |
| louvain | GM23815 | 7 | 16 | Fibroblast | 1 | True | NA | NA | NA | NA |
| louvain | GM23815 | 7 | 17 | EarlyPluripotency | 122 | False | -2.776 | +0.083 | -0.072 | +9.561 |
| louvain | GM23815 | 7 | 18 | EarlyPluripotency | 289 | False | -7.347 | +0.001 | -0.055 | +11.933 |
| louvain | GM23815 | 10 | 0 | Fibroblast | 2 | True | NA | NA | NA | NA |
| louvain | GM23815 | 10 | 1 | EarlyPluripotency | 8 | True | NA | NA | NA | NA |
| louvain | GM23815 | 10 | 2 | EarlyPluripotency | 2730 | False | -1.570 | +0.155 | -0.037 | +7.637 |
| louvain | GM23815 | 10 | 3 | EarlyPluripotency | 30 | False | -3.984 | -0.036 | -0.096 | +11.657 |
| louvain | GM23815 | 10 | 4 | Fibroblast | 9 | True | NA | NA | NA | NA |
| louvain | GM23815 | 10 | 5 | PartialReprog | 5 | True | NA | NA | NA | NA |
| louvain | GM23815 | 10 | 6 | EarlyPluripotency | 330 | False | -4.946 | -0.115 | -0.078 | +14.675 |
| louvain | GM23815 | 10 | 7 | PartialReprog | 10 | True | NA | NA | NA | NA |
| louvain | GM23815 | 10 | 8 | NonReprog | 1794 | False | +5.408 | +0.423 | +0.035 | +2.530 |
| louvain | GM23815 | 10 | 9 | Fibroblast | 1362 | False | +3.204 | +0.380 | +0.036 | +3.122 |
| louvain | GM23815 | 10 | 10 | Pluripotency | 181 | False | -6.660 | -0.009 | +0.008 | +11.955 |
| louvain | GM23815 | 10 | 11 | PartialReprog | 10 | True | NA | NA | NA | NA |
| louvain | GM23815 | 10 | 12 | PartialReprog | 25 | False | +3.381 | +0.194 | +0.007 | +4.677 |
| louvain | GM23815 | 10 | 13 | EarlyPluripotency | 432 | False | -3.794 | +0.235 | +0.071 | +6.874 |
| louvain | GM23815 | 10 | 15 | EarlyPluripotency | 761 | False | +1.379 | +0.278 | -0.002 | +7.499 |
| louvain | GM23815 | 10 | 17 | EarlyPluripotency | 121 | False | +1.877 | +0.299 | +0.001 | +6.431 |
| louvain | GM23815 | 10 | 18 | EarlyPluripotency | 5 | True | NA | NA | NA | NA |

## Task 2 — cross-correlation (Spearman ρ, age vs MD; permutation shuffles MD)

| cell_line | n | min_cells | n_perm | n_boot | seed | boot_seed | rho | rho_null | rho_p | rho_ci_lo | rho_ci_hi | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | 33 | 20 | 200 | 200 | 20260914 | 20260918 | +0.782 | -0.014 | +0.005 | +0.603 | +0.879 | permutation shuffles MD, holds age; unit=cluster×timepoint |
| GM23815 | 34 | 20 | 200 | 200 | 20260914 | 20260918 | +0.766 | +0.013 | +0.005 | +0.577 | +0.878 | permutation shuffles MD, holds age; unit=cluster×timepoint |

## Task 2 — MD composition test (size-matched random gene-set null, not permuted weights)

| cell_line | cluster | endpoint | n_cells_d0 | n_cells_end | md_d0 | md_end | decline | p | n_random_ge_real | n_random | pass_p | ok | reason | null_kind | n_MD_genes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM00731 | 0 | d0→d7 | 3998 | 9 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 0 | d0→d10 | 3998 | 1 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 1 | d0→d7 | 0 | 9 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 1 | d0→d10 | 0 | 2313 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 2 | d0→d7 | 0 | 2092 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 2 | d0→d10 | 0 | 4 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 3 | d0→d7 | 1 | 444 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 3 | d0→d10 | 1 | 5 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 4 | d0→d7 | 9 | 1732 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 4 | d0→d10 | 9 | 31 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 5 | d0→d7 | 2 | 48 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 5 | d0→d10 | 2 | 21 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 6 | d0→d7 | 0 | 1 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 6 | d0→d10 | 0 | 1029 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 7 | d0→d7 | 0 | 884 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 7 | d0→d10 | 0 | 22 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 8 | d0→d7 | 0 | 13 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 8 | d0→d10 | 0 | 0 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 9 | d0→d7 | 0 | 909 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 9 | d0→d10 | 0 | 52 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 10 | d0→d7 | 56 | 191 | +0.356 | +0.347 | +0.009 | +0.035 | +6.000 | 200 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 10 | d0→d10 | 56 | 527 | +0.356 | +0.417 | -0.061 | +0.886 | +177.000 | 200 | False | True | NA | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 11 | d0→d7 | 0 | 64 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 11 | d0→d10 | 0 | 5 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 12 | d0→d7 | 569 | 1 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 12 | d0→d10 | 569 | 0 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 13 | d0→d7 | 0 | 3 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 13 | d0→d10 | 0 | 0 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 14 | d0→d7 | 0 | 354 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 14 | d0→d10 | 0 | 80 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 15 | d0→d7 | 372 | 3 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 15 | d0→d10 | 372 | 0 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 16 | d0→d7 | 2 | 32 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 16 | d0→d10 | 2 | 7 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 17 | d0→d7 | 0 | 0 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 17 | d0→d10 | 0 | 243 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 18 | d0→d7 | 0 | 6 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM00731 | 18 | d0→d10 | 0 | 0 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 205 |
| GM23815 | 0 | d0→d7 | 6454 | 14 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 0 | d0→d10 | 6454 | 2 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 1 | d0→d7 | 0 | 3240 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 1 | d0→d10 | 0 | 8 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 2 | d0→d7 | 0 | 1 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 2 | d0→d10 | 0 | 2730 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 3 | d0→d7 | 0 | 2537 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 3 | d0→d10 | 0 | 30 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 4 | d0→d7 | 2 | 2210 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 4 | d0→d10 | 2 | 9 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 5 | d0→d7 | 9 | 118 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 5 | d0→d10 | 9 | 5 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 6 | d0→d7 | 0 | 1819 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 6 | d0→d10 | 0 | 330 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 7 | d0→d7 | 0 | 20 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 7 | d0→d10 | 0 | 10 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 8 | d0→d7 | 0 | 8 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 8 | d0→d10 | 0 | 1794 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 9 | d0→d7 | 0 | 3 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 9 | d0→d10 | 0 | 1362 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 10 | d0→d7 | 0 | 977 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 10 | d0→d10 | 0 | 181 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 11 | d0→d7 | 0 | 65 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 11 | d0→d10 | 0 | 10 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 12 | d0→d7 | 0 | 106 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 12 | d0→d10 | 0 | 25 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 13 | d0→d7 | 56 | 164 | +0.379 | +0.159 | +0.220 | +0.005 | +0.000 | 200 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 13 | d0→d10 | 56 | 432 | +0.379 | +0.235 | +0.145 | +0.005 | +0.000 | 200 | True | True | NA | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 14 | d0→d7 | 886 | 5 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 14 | d0→d10 | 886 | 0 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 15 | d0→d7 | 0 | 0 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 15 | d0→d10 | 0 | 761 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 16 | d0→d7 | 361 | 1 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 16 | d0→d10 | 361 | 0 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 17 | d0→d7 | 0 | 122 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 17 | d0→d10 | 0 | 121 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 18 | d0→d7 | 0 | 289 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |
| GM23815 | 18 | d0→d10 | 0 | 5 | NA | NA | NA | NA | NA | 200 | NA | False | n_cells<20 at an endpoint | size_matched_random_gene_sets_mean_expression_bin | 204 |

Null difference: MD is a fixed gene list, so the null is 200 draws of the same n genes matched on mean-expression bin, scored with AddModuleScore. The frozen ruler's null is 200 permutations of its weights. Those nulls are not interchangeable.

## Task 2 reading

Fired key: `instruments_disagree`. comparison_valid=True.

MD moves per cluster, the frozen ruler does not. The instruments disagree on the same cells. What would settle it: an independent fibroblast age instrument scored on these same Louvain clusters, or the authors' Seurat/sctransform object with the frozen ruler projected onto it. No winner is declared.

Aged-donor MD null passes (size-matched gene-set null): GM00731 cluster 10 label=NonReprog d0→d7 MD decline=+0.009 p=+0.035 n_ge=6/200 n_d0=56 n_end=191 md_d0=+0.356 md_end=+0.347.

Aged-donor Louvain clusters pairable at both endpoints (ok=True): n=1 (10:NonReprog). PartialReprog among them: n=0.

Cross-correlation frozen age vs MD (cluster×timepoint ≥20 cells): GM00731 n=33 ρ=+0.782 null=-0.014 p=+0.005 CI=[+0.603, +0.879]; GM23815 n=34 ρ=+0.766 null=+0.013 p=+0.005 CI=[+0.577, +0.878].

## Task 3 pre-registration (verbatim, written before any external-cohort MD ρ)

Written **before any external-cohort MD ρ**, 2026-09-19. Report-only. No gate.

Their aging signature comes from Fleischer et al. 2018 (~100 dermal fibroblasts); our ruler comes from GTEx cultured fibroblasts (n=652). If the Fleischer DE sets are buildable from GSE113957 exactly as published (DESeq2 old-vs-young, adjusted p < 0.05 and log2FC > 0.5, old/young from PCA separation), report each instrument's Spearman ρ against donor age in the external cohorts already validated in FINDINGS_FIBRO2.md (CELLxGENE a19d1667-a7b5-4556-9e5f-f9bfa690c0f1, n=93, 10x; GSE226189, n=82, bulk), with nulls and CIs. One table. No gate.

If GSE113957's format prevents exact reconstruction, say so and report the aging-signature comparison as not run rather than approximating. MD and the frozen ruler are still scored. Frozen ruler numbers are read from results/fibro2/; the ruler is not refit.


## Task 3 — instrument validity vs donor age (report-only, no gate)

| instrument | cohort | accession | n | n_donors | rho | rho_null | rho_p | rho_ci_lo | rho_ci_hi | r | cal_r2 | r2 | platform | n_MD_mapped | n_MD_missing | note | reason |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| frozen_GTEx_fibroblast_ruler | CELLxGENE a19d1667-a7b5-4556-9e5f-f9bfa690c0f1 | a19d1667-a7b5-4556-9e5f-f9bfa690c0f1 | +93.000 | +93.000 | +0.455 | -0.003 | +0.005 | +0.253 | +0.605 | +0.414 | +0.171 | -9.664 | 10x (restricted) | NA | NA | NA | NA |
| frozen_GTEx_fibroblast_ruler | GSE226189 bulk | GSE226189 | +82.000 | +82.000 | +0.549 | +0.005 | +0.005 | +0.374 | +0.667 | +0.547 | +0.299 | -9.963 | NA | NA | NA | NA | NA |
| MD_AddModuleScore | CELLxGENE a19d1667-a7b5-4556-9e5f-f9bfa690c0f1 | a19d1667-a7b5-4556-9e5f-f9bfa690c0f1 | +93.000 | +93.000 | +0.253 | -0.004 | +0.005 | +0.025 | +0.408 | +0.251 | +0.063 | -9.655 | 10x (restricted) | +202.000 | +3.000 | NA | NA |
| MD_AddModuleScore | GSE226189 bulk | GSE226189 | +82.000 | +82.000 | +0.082 | +0.011 | +0.259 | -0.132 | +0.291 | +0.051 | +0.003 | -7.603 | bulk_rnaseq | +204.000 | +1.000 | NA | NA |
| Fleischer_DESeq2_aging_signature | CELLxGENE a19d1667-a7b5-4556-9e5f-f9bfa690c0f1 | a19d1667-a7b5-4556-9e5f-f9bfa690c0f1 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | not run | GSE113957_fpkm.txt.gz values are not integer counts (cannot apply DESeq2). FINDINGS_FIBRO2.md: series matrix HEAD 404. STAR Methods require DESeq2 old-vs-young with PCA-based old/young split. Not approximating from FPKM and not substituting mmc3 Age up (n=1533) / Age down (n=2007). Aging-signature comparison reported as not run. |
| Fleischer_DESeq2_aging_signature | GSE226189 bulk | GSE226189 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | not run | GSE113957_fpkm.txt.gz values are not integer counts (cannot apply DESeq2). FINDINGS_FIBRO2.md: series matrix HEAD 404. STAR Methods require DESeq2 old-vs-young with PCA-based old/young split. Not approximating from FPKM and not substituting mmc3 Age up (n=1533) / Age down (n=2007). Aging-signature comparison reported as not run. |

Aging signature: not run. GSE113957_fpkm.txt.gz values are not integer counts (cannot apply DESeq2). FINDINGS_FIBRO2.md: series matrix HEAD 404. STAR Methods require DESeq2 old-vs-young with PCA-based old/young split. Not approximating from FPKM and not substituting mmc3 Age up (n=1533) / Age down (n=2007). Aging-signature comparison reported as not run.

Frozen ruler numbers are read from `results/fibro2/`. The ruler is not refit.

## What this changes in FINDINGS_FIBRO.md and FINDINGS_FIBRO3.md (quote, not edited)

This file does not edit `FINDINGS_FIBRO.md`, `FINDINGS_FIBRO2.md`, `FINDINGS_FIBRO3.md`, `FINDINGS_MD.md`, or `FALSIFICATION.md`.

FINDINGS_FIBRO.md Stage 2 verdict sentence, quoted:

> age score does not decline → the ruler reads a static donor property, not a modifiable state. Step 3 is not supported by this data.

That d0→d10 all-cell verdict stands. Nothing here promotes d7 to the primary Stage 2 result.

FINDINGS_FIBRO3.md per-cluster d0→d7 sentence, quoted:

> Aged GM00731 d0→d7 (the all-cell T-B cell had p=+0.005): cluster 0 decline=-0.582 p=+0.652 n_ge=130/200 n_d0=4941 n_d7=2719; cluster 1 decline=+2.215 p=+0.144 n_ge=28/200 n_d0=76 n_d7=850; cluster 2 skipped (n_cells<20 at an endpoint).

Task 1 fired `negative_holds`. The quoted FIBRO3 per-cluster sentence stands. The frozen ruler does not move in these cells at k=3, 8, 15, or Louvain.

FINDINGS_FIBRO3.md Task 1 cell-intrinsic sentence, quoted (d7–d10 composition, a different claim):

> Cell-intrinsic. The same cells' clusters individually reverse from d7 to d10. Then reprogramming genuinely pushes cells back along the age axis after d7, and that is a biological claim worth its own test. clusters=[0, 1, 2].

That sentence is about d7→d10 reversal within joint k=3 clusters. It is not the per-cluster d0→d7/d0→d10 random-direction negative this prompt re-tests. It is not edited.

FINDINGS_MD.md closed the gate because the MD score's computation was not publicly retrievable. This file opens that computation from STAR Methods and runs the instrument comparison. FINDINGS_MD.md is not edited.

## Limitations

1. Two donors in GSE297234.
2. Our reproduction of their metric is ours, not theirs. Seurat v5 / sctransform v2 was not run; LogNormalize + residualized PCA Louvain is a reported deviation. AddModuleScore is a Python reimplementation of the published algorithm, not Seurat's C++/R object.
3. Clusters are not lineage-tracked. PartialReprog / NonReprog labels are argmax of mmc3 Reprog_cell_state_signatures AddModuleScore, not Slingshot trajectories from day-0 fibroblasts.
4. Cross-platform shift from GTEx bulk polyA (RNASeQCv2.4.2) to 10x 3' scRNA-seq pseudobulk.
5. Frozen ruler trained on GTEx V10 cultured fibroblasts, public `AGE` 10-year bins.
6. Missing overlap genes at z=0. Nothing fitted on GSE297234 or on the external matrices.
7. GSE113957 is FPKM; the Fleischer DESeq2 aging signature was not rebuilt and was not substituted from mmc3 Age up/down.
8. MIN_CELLS=20. Cluster×timepoint rows below that threshold are not scored.
9. Seed `20260914` (permutation / random directions / clustering / AddModuleScore sampling). Bootstrap seed `20260918`. n_perm=200, n_boot=200, n_random=200.
10. GSE325735 not opened. d7 is not the Stage 2 endpoint.
11. Louvain resolution 0.8 is the Seurat default; STAR Methods do not state the resolution used for Figure 3.
12. Louvain cluster occupancy is time-specific: on aged GM00731 only cluster 10 has ≥20 cells at both d0 and d7 (n_d0=56, n_d7=191). The MD per-cluster null therefore rests on that one pairable cluster.
13. GSE226189 geneCOUNT `Tracking_ID` is Ensembl. Symbols were mapped from the GSE297234 Cell Ranger feature table plus the frozen ruler; n_MD_mapped=204 n_MD_missing=1 (PRSS2; t3_report.txt).
14. CELLxGENE a19d1667 MD: n_MD_mapped=202 n_MD_missing=3 (GREM1, PRSS2, TWIST2; t3_report.txt).

## Open questions

1. Would the authors' sctransform v2 object change which Louvain clusters are PartialReprog vs NonReprog?
2. If the Fleischer integer-count matrix and PCA old/young split were released, how does that aging signature correlate with the frozen GTEx ruler on CELLxGENE a19d1667 and GSE226189?
3. Slingshot from day-0 fibroblast clusters was not run (optional, report-only).

## Columns actually read

- **hallmark_local_set_names** source=`<repo>\results\lowdim\genesets\MSigDB_Hallmark_2020.txt` columns=['Adipogenesis', 'Allograft Rejection', 'Androgen Response', 'Angiogenesis', 'Apical Junction', 'Apical Surface', 'Apoptosis', 'Bile Acid Metabolism', 'Cholesterol Homeostasis', 'Coagulation', 'Complement', 'DNA Repair', 'E2F Targets', 'Epithelial Mesenchymal Transition', 'Estrogen Response Early', 'Estrogen Response Late', 'Fatty Acid Metabolism', 'G2-M Checkpoint', 'Glycolysis', 'Hedgehog Signaling', 'Hypoxia', 'IL-2/STAT5 Signaling', 'IL-6/JAK/STAT3 Signaling', 'Inflammatory Response', 'Interferon Alpha Response', 'Interferon Gamma Response', 'KRAS Signaling Dn', 'KRAS Signaling Up', 'Mitotic Spindle', 'Myc Targets V1', 'Myc Targets V2', 'Myogenesis', 'Notch Signaling', 'Oxidative Phosphorylation', 'PI3K/AKT/mTOR  Signaling', 'Pancreas Beta Cells', 'Pperoxisome', 'Protein Secretion', 'Reactive Oxygen Species Pathway', 'Spermatogenesis', 'TGF-beta Signaling', 'TNF-alpha Signaling via NF-kB', 'UV Response Dn', 'UV Response Up', 'Unfolded Protein Response', 'Wnt-beta Catenin Signaling', 'Xenobiotic Metabolism', 'heme Metabolism', 'mTORC1 Signaling', 'p53 Pathway']
- **msigdb_gmt_set_names** source=`https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt` columns=['HALLMARK_ADIPOGENESIS', 'HALLMARK_ALLOGRAFT_REJECTION', 'HALLMARK_ANDROGEN_RESPONSE', 'HALLMARK_ANGIOGENESIS', 'HALLMARK_APICAL_JUNCTION', 'HALLMARK_APICAL_SURFACE', 'HALLMARK_APOPTOSIS', 'HALLMARK_BILE_ACID_METABOLISM', 'HALLMARK_CHOLESTEROL_HOMEOSTASIS', 'HALLMARK_COAGULATION', 'HALLMARK_COMPLEMENT', 'HALLMARK_DNA_REPAIR', 'HALLMARK_E2F_TARGETS', 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', 'HALLMARK_ESTROGEN_RESPONSE_EARLY', 'HALLMARK_ESTROGEN_RESPONSE_LATE', 'HALLMARK_FATTY_ACID_METABOLISM', 'HALLMARK_G2M_CHECKPOINT', 'HALLMARK_GLYCOLYSIS', 'HALLMARK_HEDGEHOG_SIGNALING', 'HALLMARK_HEME_METABOLISM', 'HALLMARK_HYPOXIA', 'HALLMARK_IL2_STAT5_SIGNALING', 'HALLMARK_IL6_JAK_STAT3_SIGNALING', 'HALLMARK_INFLAMMATORY_RESPONSE', 'HALLMARK_INTERFERON_ALPHA_RESPONSE', 'HALLMARK_INTERFERON_GAMMA_RESPONSE', 'HALLMARK_KRAS_SIGNALING_DN', 'HALLMARK_KRAS_SIGNALING_UP', 'HALLMARK_MITOTIC_SPINDLE', 'HALLMARK_MTORC1_SIGNALING', 'HALLMARK_MYC_TARGETS_V1', 'HALLMARK_MYC_TARGETS_V2', 'HALLMARK_MYOGENESIS', 'HALLMARK_NOTCH_SIGNALING', 'HALLMARK_OXIDATIVE_PHOSPHORYLATION', 'HALLMARK_P53_PATHWAY', 'HALLMARK_PANCREAS_BETA_CELLS', 'HALLMARK_PEROXISOME', 'HALLMARK_PI3K_AKT_MTOR_SIGNALING', 'HALLMARK_PROTEIN_SECRETION', 'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY', 'HALLMARK_SPERMATOGENESIS', 'HALLMARK_TGF_BETA_SIGNALING', 'HALLMARK_TNFA_SIGNALING_VIA_NFKB', 'HALLMARK_UNFOLDED_PROTEIN_RESPONSE', 'HALLMARK_UV_RESPONSE_DN', 'HALLMARK_UV_RESPONSE_UP', 'HALLMARK_WNT_BETA_CATENIN_SIGNALING', 'HALLMARK_XENOBIOTIC_METABOLISM']
- **mmc3_MD_signatures** source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx` columns=['MD score', 'TGFB score']
- **mmc3_Aging_signatures** source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx` columns=['Age up', 'Age down']
- **mmc3_Fibroblast_subtype_signatures** source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx` columns=['PI16_univ', 'LRRC15_myo', 'COL3A1_myo']
- **mmc3_Reprog_cell_state_signatures** source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx` columns=['Fibroblast', 'PartialReprog', 'EarlyPluripotency', 'Pluripotency', 'NonReprog']
- **GSE113957_fpkm_header** source=`<repo>\data\raw\fibro2\GSE113957\GSE113957_fpkm.txt.gz` columns=['Transcript ID', 'chr', 'start', 'end', 'strand', 'Length', 'Copies', 'Annotation/Divergence', '101_19yr_Female_Caucasian', '102_19yr_Male_Caucasian', '103_19yr_Male_Caucasian', '104_19yr_Male_Caucasian', '105_20yr_Male_Caucasian', '106_20yr_Female_Caucasian', '107_31yr_Female_Caucasian', '108_31_female_Caucasian', '109_32_male_Unknown', '110_32_female_BlackPuertoRican', '111_33yr_Male_Caucasian', '112_33yr_Male_Caucasian']
- **cluster_stage2_geo_samples** source=`<repo>\results\fibro\stage2_geo_samples.csv` columns=['gsm', 'title', 'source', 'organism', 'molecule', 'description', 'characteristics', 'growth_protocol', 'cell_line', 'cell_type', 'tissue', 'treatment', 'day', 'age_from_characteristics', 'age_from_growth_protocol', 'library_name', 'fields_read', 'age_years', 'age_source']
- **fibro3_t1_cell_labels** source=`<repo>\results\fibro3\t1_cell_labels.csv` columns=['barcode', 'gsm', 'file', 'cell_line', 'day', 'age_years', 'age_source', 'umi', 'n_genes', 'mito_frac', 'cluster']
- **fibro3_t1_cluster_table** source=`<repo>\results\fibro3\t1_cluster_table.csv` columns=['cell_line', 'day', 'cluster', 'age_years', 'gsm', 'age_source', 'n_cells', 'n_timepoint', 'frac_timepoint', 'median_UMI', 'median_genes', 'mito_frac_median_cell', 'below_min_cells', 'age_score', 'pluri_with_OSKM', 'pluri_without_OSKM', 'pluri_primary', 'pluri_primary_name']
- **fibro3_t3_null** source=`<repo>\results\fibro3\t3_null.csv` columns=['cell_line', 'cluster', 'endpoint', 'd_end', 'primary', 'n_random', 'seed', 'min_cells', 'n_cells_d0', 'n_cells_end', 'ok', 'score_d0', 'score_end', 'decline', 'p', 'n_random_ge_real', 'pass_p', 'reason']
- **t1_cluster_table** source=`<repo>\results\md2\t1_cluster_table.csv` columns=['resolution', 'cell_line', 'day', 'cluster', 'n_cells', 'n_timepoint', 'frac_timepoint', 'median_UMI', 'median_genes', 'mito_frac_median_cell', 'below_min_cells', 'md_score', 'tgfb_score', 'state_Fibroblast', 'state_PartialReprog', 'state_EarlyPluripotency', 'state_Pluripotency', 'state_NonReprog', 'age_score', 'pluri_with_OSKM', 'pluri_without_OSKM', 'pluri_primary', 'age_source']
- **t3_cxg_scores** source=`<repo>\results\fibro2\ta_cxg_a19d1667_scores.csv` columns=['sample', 'donor', 'age_years', 'age_source', 'n_cells', 'assay', 'age_score']
- **t3_cxg_h5ad_obs** source=`<repo>\data\raw\fibro2\cxg\8ae4471a-2e70-47b5-bc42-4f0df26a9996.h5ad` columns=['cluster_robust', 'donor_id', 'development_stage_ontology_term_id', 'sex_ontology_term_id', 'self_reported_ethnicity_ontology_term_id', 'disease_ontology_term_id', 'tissue_ontology_term_id', 'cell_type_ontology_term_id', 'assay_ontology_term_id', 'suspension_type', 'Author', 'Accession (Sample)', 'Aligner', 'Genome', 'Donor identifier', 'Sample identifier', 'Internal sample identifier', 'Sequencer', 'Age', 'Race', 'Ethnicity', 'Sample location', 'Condition', 'Condition (other)', 'is_primary_data', 'tissue_type', 'cell_type', 'assay', 'disease', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage', 'observation_joinid']
- **t3_gse226189_scores** source=`<repo>\results\fibro2\ta_GSE226189_near_miss_bulk_scores.csv` columns=['sample', 'donor', 'age_years', 'age_source', 'title', 'gsm', 'age_score']
- **t3_gse226189_gene_index** source=`<repo>\data\raw\fibro2\GSE226189\GSE226189_RAW.tar:Tracking_ID` columns=['ENSG00000000003', 'ENSG00000000005', 'ENSG00000000419', 'ENSG00000000457', 'ENSG00000000460', 'ENSG00000000938', 'ENSG00000000971', 'ENSG00000001036', 'ENSG00000001084', 'ENSG00000001167', 'ENSG00000001460', 'ENSG00000001461', 'ENSG00000001497', 'ENSG00000001561', 'ENSG00000001617', 'ENSG00000001626', 'ENSG00000001629', 'ENSG00000001630', 'ENSG00000001631', 'ENSG00000002016']
- **t3_gse226189_10x_map** source=`<repo>\data\processed\md2\allcell_counts_GM00731.npz` columns=['gene_id', 'symbol']

## Gene-set names actually read

- **HALLMARK_EMT_local** name=`Epithelial Mesenchymal Transition` n=200 source=`<repo>\results\lowdim\genesets\MSigDB_Hallmark_2020.txt`
- **HALLMARK_TGFB_local** name=`TGF-beta Signaling` n=54 source=`<repo>\results\lowdim\genesets\MSigDB_Hallmark_2020.txt`
- **HALLMARK_EMT_gmt** name=`HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION` n=200 source=`https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt`
- **HALLMARK_TGFB_gmt** name=`HALLMARK_TGF_BETA_SIGNALING` n=54 source=`https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt`
- **MD_built** name=`MD = HALLMARK_EMT + SNAI1,ZEB1,ZEB2,TWIST1,TWIST2` n=205 source=`https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt`
- **TGFB_built** name=`HALLMARK_TGF_BETA_SIGNALING` n=54 source=`https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt`
- **mmc3_MD_signatures_MD score** name=`MD_signatures:MD score` n=205 source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx`
- **mmc3_MD_signatures_TGFB score** name=`MD_signatures:TGFB score` n=54 source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx`
- **mmc3_Aging_signatures_Age up** name=`Aging_signatures:Age up` n=1533 source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx`
- **mmc3_Aging_signatures_Age down** name=`Aging_signatures:Age down` n=2007 source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx`
- **mmc3_Fibroblast_subtype_signatures_PI16_univ** name=`Fibroblast_subtype_signatures:PI16_univ` n=778 source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx`
- **mmc3_Fibroblast_subtype_signatures_LRRC15_myo** name=`Fibroblast_subtype_signatures:LRRC15_myo` n=472 source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx`
- **mmc3_Fibroblast_subtype_signatures_COL3A1_myo** name=`Fibroblast_subtype_signatures:COL3A1_myo` n=278 source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx`
- **mmc3_Reprog_cell_state_signatures_Fibroblast** name=`Reprog_cell_state_signatures:Fibroblast` n=200 source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx`
- **mmc3_Reprog_cell_state_signatures_PartialReprog** name=`Reprog_cell_state_signatures:PartialReprog` n=200 source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx`
- **mmc3_Reprog_cell_state_signatures_EarlyPluripotency** name=`Reprog_cell_state_signatures:EarlyPluripotency` n=200 source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx`
- **mmc3_Reprog_cell_state_signatures_Pluripotency** name=`Reprog_cell_state_signatures:Pluripotency` n=200 source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx`
- **mmc3_Reprog_cell_state_signatures_NonReprog** name=`Reprog_cell_state_signatures:NonReprog` n=200 source=`<repo>\results\md\raw\suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx`

## Files

- `src/md2_common.py`, `src/md2_genesets.py`, `src/md2_score.py`, `src/md2_cluster.py`, `src/md2_task1.py`, `src/md2_task2.py`, `src/md2_task3.py`, `src/md2_findings.py`, `src/md2_run.py`
- `results/md2/`
- `FINDINGS_MD2.md`
- `PROGRESS_MD2.md`


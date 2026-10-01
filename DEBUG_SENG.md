# DEBUG_SENG — why the four named hits missed our cosine gate

Diagnostic only. `FINDINGS_SENG.md` is unchanged. Its fired key `pipeline_disagrees` stands. No gate in that file was re-run, re-tuned, or overturned.

Seeds: `20260914` (splits, N_NT nulls, permutations of `u`). The PD32 cosine below matches `results/seng/screen.csv` to a maximum absolute difference of `9.7e-17`.

## The fact to explain

Named hits `CRA_E2F3`, `CRA_EZH2`, `CRI_STAT3`, `CRI_ZFX`. In the SENG run only `CRI_ZFX` had `cos_g > 0` at BH `q <= 0.05` (`E2F3` q=0.092, `EZH2` q=0.092, `STAT3` q=0.100, `ZFX` q=0.034). All four had positive cosine and positive delta. The deltas were small (`+13` to `+28`).

## H1 — wrong or partial file

**Verdict: resolved. `file_explains` does not fire.** This h5ad is the perturb-seq screen. No second file on either figshare record cited by the paper carries a passage label for the perturbed cells. No additional file was downloaded.

Records opened:

- Perturb-seq: https://api.figshare.com/v2/articles/30898748 (DOI https://doi.org/10.6084/m9.figshare.30898748, HTML https://figshare.com/articles/dataset/Single_Cell_Perturb-seq_data_in_h5ad_format_For_the_paper_published_in_PNAS_/30898748). Version 2. Description: the combined 10x count matrix and Perturb-seq guide UMIs.
- Mouse liver, the other DOI in the data statement: https://api.figshare.com/v2/articles/30899018 (DOI https://doi.org/10.6084/m9.figshare.30899018). Description: bulk RNA-seq counts for young/aged liver, GFP vs EZH2. Not the fibroblast screen.
- Paper: https://www.pnas.org/doi/10.1073/pnas.2515183123 and https://pmc.ncbi.nlm.nih.gov/articles/PMC12799168/

Data availability statement, from the paper:

> All data needed to understand and assess the conclusions of this study are included in the text, figures, and supplementary materials. The full set of perturb-seq data and mouse liver gene expression data are available from figshare (https://doi.org/10.6084/m9.figshare.30898748, and https://doi.org/10.6084/m9.figshare.30899018) (47, 48).

Files as listed:

| record | file | size (bytes) |
|---|---|---|
| 30898748 | `exp3_merged_not_normalized.h5ad` | 1426871504 |
| 30899018 | `mouse_rnaseq_sample_name.csv` | 233 |
| 30899018 | `mouse_liver_raw_count_figshare.csv` | 1669439 |

The perturb-seq record contains one file. It is an h5ad. The mouse record contains two csv files and no h5ad. The SI Appendix says the WT passage stages were pooled and identified with cell-membrane barcodes, and that the CRA/CRI Perturb-seq was a separate 5' library. Those WT passage labels are the boolean columns `PD14`, `PD26`, `PD32` inside this h5ad. They are not a second file.

## H2 — wrong LATE reference

**Verdict: resolved. `reference_explains` does not fire.** No obs column states the passage of the perturbed cells. Recomputing the cosine with LATE = PD26 makes every named hit smaller, and none of the four has `cos_g > 0` and BH `q <= 0.05` under that choice. `CRI_ZFX`, the one hit that cleared `q <= 0.05` at PD32, does not clear it at PD26 (q=0.057).

`obs` has 814 columns. The full list, with dtype and up to 10 unique values, is `results/seng_why/obs_columns.csv`. Seven columns are low-cardinality and are not guide UMI counts:

| column | dtype | unique values | on perturbed cells |
|---|---|---|---|
| PD14 | bool | true, false | false only (0 cells flagged) |
| PD26 | bool | true, false | false only (0 cells flagged) |
| PD32 | bool | true, false | false only (0 cells flagged) |
| CRA | int64 | 0, 1 | 0 and 1 |
| CRI | int64 | 0, 1 | 0 and 1 |
| WT | int64 | 0, 1 | 0 only |
| SampleName | int8 | 0, 1, 2, 3 | 3 only |

The other 807 columns are `n_counts` (float32), `log_counts` (float32), `n_genes` (int64), and 804 guide-UMI columns (float64, names `{TF}_{CRA|CRI}{k}` or `NT_{k}`). There is no `batch`, `library`, `lane`, hashing, or `orig.ident` column. Guide columns are zero on every WT cell and are UMI counts, not a passage label.

`SampleName` is a library id, not a passage:

| SampleName | cells |
|---|---|
| 0 | WT PD14 only (581) |
| 1 | WT PD26 only (1880) |
| 2 | WT PD32 only (4348) |
| 3 | every non-WT cell (CRA 12264, CRI 12346, NT 2145, unassigned 2950) |

Sample 3 is the perturb-seq library. It has no PD value. The paper's protocol calls the perturbed cultures late passage, and its figures name WT PD14 as early and WT PD32 as late. That text is not an obs column. Neither LATE is chosen here.

Cosine recomputed with the same split seed, the same TMM/z panel, and the same NT null (200 draws, seed `20260914`, BH within modality). PD32 `q` is the value already stored in `FINDINGS_SENG.md`, confirmed by an exact cosine match. PD26 `q` is new and is not a re-run of that gate.

| perturbation | n | cos PD32 | q PD32 | cos PD26 | p PD26 | q PD26 |
|---|---|---|---|---|---|---|
| CRA_E2F3 | 128 | +0.248 | 0.092 | +0.143 | 0.0050 | 0.138 |
| CRA_EZH2 | 107 | +0.197 | 0.092 | +0.101 | 0.0249 | 0.231 |
| CRI_STAT3 | 96 | +0.155 | 0.100 | +0.099 | 0.0149 | 0.092 |
| CRI_ZFX | 157 | +0.316 | 0.034 | +0.137 | 0.0050 | 0.057 |

Full screen: `results/seng_why/cos_by_late.csv`.

## H3 — metric difference

**Verdict: the test that the paper's wording supports does not fire `metric_explains`.** All-gene Pearson `R_rej` and our `cos_g` rank the screen almost the same way (Spearman −0.984). The four named hits are not high under `R_rej` and low under `cos_g`. Numerical recovery of Table 1 is unresolved, because the paper states no gene filter and none was searched for.

### Quoted definition

Main text (PNAS / PMC):

> We quantify the degree of reversal by R_rej, defined as the correlation between log fold change of gene expression of late passage cells vs. early passage cells and log fold change of late passage cells with a TF perturbation vs. late passage cells with a nontargeting control. TF perturbations with a significant negative R_rej indicated the perturbation reversed the gene expression changes due to replicative aging. We then ranked R_rej and selected the top hits (TFs with r =< −0.3 Table 1).

SI Appendix, Materials and Methods, "Single-cell RNA Sequencing (scRNA-seq) Analysis":

> 10x Genomics Cell Ranger and Scanpy (7) were computational packages used to analyze scRNA-seq data. The potential rejuvenation effect of the TF perturbations was measured by how well the gene expression profile in perturbed cells mimicked the gene expression profile in the early passage cells, compared to the late passage cells. We first computed the gene expression fold changes (log2) in the late passage cells compared to the early passage cells. Then, for each TF perturbation (CRA or CRI) we computed the gene expression fold changes (log2) by comparing the cells with the guides targeting the TF and the cells with the non-targeting guides (NT). We then computed the Pearson correlation of the gene expression changes between late passage and early passage cells and those between CRA or CRI targeted cells and the NT control (R_rej). The TF perturbations with the strongest negative R_rej had the most significant change in gene expression towards being like earlier passage cells.

### What the paper fixes, and what was assumed

Fixed by the paper:

- The two vectors are log2 fold changes: late vs early, and perturbation vs NT of the same modality.
- The correlation is Pearson (SI). The main text says "correlation" and does not name Spearman.
- Rejuvenating is the negative direction.

Assumptions, not facts:

- Gene selection. Neither paragraph states a significance cut or a differentially-expressed subset. The correlation was computed on all 36,601 genes in the file. No cutoff was tried.
- The log2 fold change is a pseudobulk library-size CPM ratio, `log2((CPM_group+1)/(CPM_reference+1))`, with no TMM. The paper does not state pseudobulk versus per-cell mean, the library target, or the pseudocount.
- The aging arm uses all WT PD32 cells versus all WT PD14 cells, the two passages the paper's figures name. It does not use the half-A split. The perturbation arm uses every cell of that guide pool versus every NT cell of the same modality.
- Rank 1 means the most negative `R_rej`, and rank 1 under `cos_g` means the most positive cosine.

| perturbation | R_rej | published Table 1 | rank R_rej | our cos_g | rank cos_g | n in modality |
|---|---|---|---|---|---|---|
| CRA_E2F3 | −0.281 | −0.53 | 1 / 167 | +0.248 | 1 / 167 | 167 |
| CRA_EZH2 | −0.203 | −0.36 | 9 / 167 | +0.197 | 6 / 167 | 167 |
| CRI_STAT3 | −0.155 | −0.41 | 49 / 173 | +0.155 | 37 / 173 | 173 |
| CRI_ZFX | −0.355 | −0.51 | 1 / 173 | +0.316 | 2 / 173 | 173 |

Spearman correlation between `R_rej` and `cos_g` across all 340 perturbations with n ≥ 30: **−0.984** (CRA −0.983, n=167; CRI −0.986, n=173). The negative sign is the two conventions: their rejuvenating direction is negative, ours is positive. The orderings agree.

`E2F3`, `EZH2`, and `ZFX` are near the top under both scores. `STAT3` is mid-pack under both. That is the opposite of "`R_rej` ranks them at the top and `cos_g` does not."

The absolute values are weaker than Table 1. Nine of the 15 CRA names in Table 1, and seven of the 15 CRI names, fall in our top 15 by `R_rej`; those same names are also high by `cos_g`. `STAT3` is not among them (our `R_rej` −0.15 versus their published −0.41). Closing that gap would require the gene list they actually correlated. Choosing a cutoff until Table 1 reappears would be tuning, and it was not done.

Per-perturbation table: `results/seng_why/rrej_vs_cos.csv`.

## H4 — gene-space restriction

**Verdict: resolved. `gene_space_explains` does not fire.** The cosine over every gene in the file is smaller than the ruler-restricted cosine for all four hits, and the ranks stay in the same part of the list. The extra genes do not create a hit the ruler cosine missed.

The ruler cosine uses 16,781 of 23,485 frozen-ruler genes present in the file (the PD32 values above, identical to the SENG run). The all-gene cosine uses the same groups and the same TMM log2-CPM (`prior.count=2`). Frozen GTEx mean and sd are applied where the gene is on the ruler. The file has 19,830 genes with no frozen mean or sd. Those were z-scored with the mean and sd of that gene's TMM log2-CPM across the analysis panel. That scale is an assumption. Genes with no counts in the panel (10,078) are z=0, as in the ruler path. Ten file columns are the symbol match for two ruler genes; the ruler cosine follows the original map (and matches `screen.csv`); the all-gene cosine keeps each file column once.

| perturbation | cos, ruler (16,781) | rank | cos, all genes | rank | \|\|d\|\| on mapped genes | \|\|d\|\| on unmapped genes |
|---|---|---|---|---|---|---|
| CRA_E2F3 | +0.248 | 1 | +0.066 | 3 | 92.9 | 133.3 |
| CRA_EZH2 | +0.197 | 6 | +0.058 | 7 | 91.6 | 127.5 |
| CRI_STAT3 | +0.155 | 37 | +0.052 | 50 | 92.2 | 206.4 |
| CRI_ZFX | +0.316 | 2 | +0.077 | 6 | 86.9 | 198.3 |

The unmapped block is a large part of `||d||`, so the numeric value depends on the panel-sd assumption. The direction of that dependence is dilution: the all-gene cosine is about a third of the ruler cosine, and `STAT3` moves from rank 37 to rank 50. A larger cosine on the extra genes is not what this assumption produced.

Full screen: `results/seng_why/cos_all_genes.csv`.

## Delta noise floor

`delta` in `FINDINGS_SENG.md` carries a positive offset of about +20.6 (positive control f=0 gave delta=+20.6), and the named-hit deltas (+13 to +28) sit within it, so they do not support a direction.

## Do the proliferation and donor-age results depend on LATE?

No. Both LATE choices leave `rho_prolif` near 0.83 / 0.86 and leave `cos(v, u)` indistinguishable from the gene-permutation null. The S+G2M shift itself does not use LATE; only `cos_g` and `v` change.

| LATE | rho_prolif CRA | rho_prolif CRI | cos(v, u) | N1 p |
|---|---|---|---|---|
| PD32 | +0.828 | +0.862 | −0.040 | 1.000 |
| PD26 | +0.819 | +0.859 | −0.006 | 0.736 |

## Fired reading

**`unexplained`.**

- `metric_explains` does not fire. Under the stated Pearson definition with no gene filter, the named hits do not rank at the top under `R_rej` and away from the top under `cos_g`. The two scores agree (Spearman −0.984). `E2F3` is rank 1 under both and still has q=0.092.
- `reference_explains` does not fire. PD26 does not put any of the four across `q <= 0.05`.
- `gene_space_explains` does not fire. The all-gene cosine is smaller and does not reproduce them.
- `file_explains` does not fire. This file is the screen.

## What is still unknown

The passage of the perturbed cells is not in this file and not in the mouse-liver file. The paper describes those cultures as late passage. That sentence was not treated as a metadata column.

The gene set inside their correlation is not stated. All-gene `R_rej` is systematically weaker than Table 1, and it leaves `CRI_STAT3` at rank 49 rather than inside their `r <= −0.3` list. Their gene list, or the code that produced Table 1, is what would be required to test a restricted correlation. A cutoff was not fit to Table 1.

Nothing in these four checks replaces `pipeline_disagrees`.

# SURVEY_TF — usable transcription-factor perturbation data

Pre-registration written before any dataset page was opened and before any search.
`results/survey/PREREG.flag` was touched in the same step. Tiers are not changed after results.

## TIERS (pre-registered)

Tier A (screen-ready): human fibroblast; >=50 Lambert TFs overexpressed or
  CRISPRa; >=10,000 genes measured; matched controls in the same batch;
  count-level or processed data downloadable.
Tier B (usable with a named caveat): meets Tier A except exactly ONE of:
  non-fibroblast human SOMATIC cell (not stem, not cancer); knockdown only;
  reduced readout (<10,000 measured genes, e.g., L1000); 10-49 TFs.
Tier C: human stem-cell or cancer-line only, or mouse only, or fails two or
  more Tier A criteria.
Not usable: no controls, or no downloadable expression data.
Each dataset gets exactly one tier with one sentence explaining why.

## PRE-REGISTERED READING (report which fired)

- >=1 Tier A -> key `screen_ready_fibroblast`. An in-silico fibroblast TF
  screen is feasible with existing data. Name the dataset(s).
- No Tier A, >=1 Tier B -> key `screen_with_caveat`. Name the best dataset and
  its one caveat in the same sentence.
- Only Tier C or nothing -> key `data_generation_needed`. The honest pitch is
  "fund generating the data," not "fund validating candidates."

## WHAT DOES NOT COUNT

- A dataset mentioned in a review or abstract but not opened.
- Counts taken from an abstract when the data page says otherwise.
- Calling a dataset "fibroblast" because the paper mentions fibroblasts
  somewhere; the perturbed cells themselves must be fibroblasts.
- Counting non-TF genes (signalling genes, chromatin modifiers not on the
  Lambert list) as TFs.
- Mouse data presented as human, or inferred L1000 genes counted as measured.
- Promoting a dataset up a tier because it "almost" qualifies.

---

The sections below are results. They are not part of the pre-registration.
The flat table is `results/survey/datasets.csv`.

## Lambert list

Human TFs v1.01 (Lambert et al. 2018, PMID 29425488) was downloaded from
https://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt
and saved as `data/reference/TF_names_v_1.01.txt`. The file has 1,639 unique
symbols. Only those symbols are counted as TFs below.

## Main table (verified only)

| Dataset | Paper | Accession | Cells | Perturbation | Lambert TFs | Readout | Genes measured | Controls | Tier |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Hs27 CRISPRa Perturb-seq | Southard 2025, Nat Genet | Zenodo 10.5281/zenodo.15200179; PRJNA1108254 | primary Hs27 fibroblasts | CRISPRa | 1,488 | 10x scRNA | 36,601 | Y | A |
| Fibroblast CRISPRa/CRISPRi Perturb-seq | Sengstack 2026, PNAS | figshare 10.6084/m9.figshare.30898748; GEO not found | human fibroblasts | CRISPRa and CRISPRi | 179 | 10x scRNA | shape 36,514 x 36,601 | Y | A |
| RPE-1 CRISPRa Perturb-seq | Southard 2025, Nat Genet | Zenodo 10.5281/zenodo.15213619 | RPE-1, not fibroblast | CRISPRa | not counted | 10x files listed, not opened | not counted | not counted | C |
| LINCS L1000 Phase I/II | Subramanian 2017, Cell | GSE92742; GSE70138 | many lines; MCH58 is the skin fibroblast line | ORF overexpression and shRNA | 645 OE; MCH58 has 1 knockdown | L1000 | 978 landmarks | Y at cell_id level | C |
| TF atlas, 90 ORF library | Joung 2023, Cell | GSE216595 | hESC | ORF overexpression | 46 | 10x scRNA | 33,694 feature rows | Y | C |
| Human TFome | Ng 2021, Nat Biotechnol | GSE159786 | hiPSC | ORF; pooled arm is barcode PCR | 4 in the RNA-seq arm | barcode PCR, plus bulk RNA-seq of 4 factors | not counted | Y for the bulk arm | C |
| Norman Perturb-seq | Norman 2019, Science | GSE133344 | K562 | overexpression (sample page also says dCas9-KRAB) | 44 | 10x scRNA | 33,694 feature rows | Y | C |
| SEUSS | Parekh 2018, Cell Systems | GSE107185 | hPSC | ORF overexpression | 60 | scRNA-seq | not counted | Y | C |
| Fibroblast core TF knockdown / monocyte TF overexpression | Li 2016, PLoS One | GSE80676 | dermal fibroblasts | OE and knockdown, few genes | 7 | microarray | 48,107 probes | Y | C |
| ENCODE K562-shX | ENCODE / HudsonAlpha | GSE33816 | K562 | shRNA | 3 | bulk RNA-seq | not counted | Y | C |
| ENCODE HepG2 TFIP11 shRNA | ENCODE | GSE88370 | HepG2 | shRNA of TFIP11 | 0 | bulk RNA-seq | not counted | not read | C |
| Genome-scale Perturb-seq | BioProject title; cell line not on that page | PRJNA831566 | not stated | CRISPRi of expressed genes | not counted | raw reads only | not counted | not read | Not usable |

Tier sentences:

- Hs27 is Tier A because the perturbed cells are primary human fibroblasts, 1,488 Lambert TFs are CRISPRa, 36,601 genes are measured, 78 non-targeting guides sit in the same guide matrix, and the count data are downloadable.
- Sengstack is Tier A because the perturbed cells are human fibroblasts, 179 Lambert TFs are CRISPRa (CRISPRi of the same genes is also present, so this is not knockdown-only), both axes of the matrix shape are above 10,000, non-targeting labels are in the same file, and the h5ad is downloadable.
- RPE-1 is Tier C because the cells are not fibroblasts and this record's own TF and gene counts were not read, so it is not scored as a one-caveat somatic screen.
- LINCS is Tier C because MCH58, the only skin-fibroblast line, has a single Lambert knockdown, the large overexpression sets are in cancer lines or immortalized kidney cells, and only 978 genes are measured landmarks.
- Joung GSE216595 is Tier C because the perturbed cells are human embryonic stem cells.
- Ng is Tier C because the perturbed cells are human pluripotent stem cells and the only transcriptome is four factors.
- Norman is Tier C because the perturbed cells are the K562 cancer line.
- Parekh is Tier C because the perturbed cells are human pluripotent stem cells.
- GSE80676 is Tier C because seven Lambert TFs is below the Tier B band of 10 to 49.
- GSE33816 is Tier C because the cells are K562 and only three Lambert TFs are knocked down.
- GSE88370 is Tier C because HepG2 is a cancer line and TFIP11 is not on the Lambert list.
- Replogle PRJNA831566 is Not usable because the opened page lists raw reads and no expression matrix was opened.

No dataset is Tier B. Nothing was moved up a tier for almost qualifying.

### How the two Tier A counts were read

Hs27. Zenodo 10.5281/zenodo.15200179 says the files are CRISPRa perturbations in Hs27 cells (fibroblasts). BioProject PRJNA1108254 says primary fibroblasts and retinal pigment epithelial cells and "1,836 transcription factors." The guide HDF5 has 10,979 guides and 1,836 target symbols, of which 1,488 are on Lambert v1.01 and 348 are not (examples: AATF, ABL1, ARID1A). Non-targeting guides: 78, in that same matrix. The Cell Ranger feature list was read from lane L4A inside the 21.6 GB zip, without downloading the zip: 36,601 Gene Expression features and 10,979 guide-capture features. The guide-capture count matches the guide HDF5, so this is the shared feature reference rather than a lane-specific detected subset. The 1.72 GB mean-population h5ad was also opened; its matrix is 10,916 by 4,914 genes. Those 4,914 genes are the analysis subset, not the measured set. The 8.05 GB singlets file was not downloaded.

Sengstack. Figshare 10.6084/m9.figshare.30898748 is one 1.43 GB h5ad, public, CC BY 4.0. The description says it is the 10x Cell Ranger matrix plus Perturb-seq guide UMIs, and that each TF was targeted by two guides. The obs name heap, read by HTTP range, contains 200 gene symbols, each with CRA and CRI labels, plus NT_1 through NT_6. 179 of the 200 symbols are on Lambert v1.01. The other 21 include ATG3, BRCA1, EZH2, HMGB1, PARP1, and SIRT1. Dataset X carries a shape attribute of 36,514 by 36,601. Gene symbols for that matrix were not extracted, so ruler coverage was not computed. Donor age in years is not on the figshare record. Obs labels include PD14, PD26, and PD32.

LINCS detail that the table compresses. Phase I sig_info has 473,647 signatures. Landmark genes in gene_info: 978. Non-landmark rows: 11,350, not counted as measured. Overexpression Lambert counts by cell_id include HEK293T 523, A375 316, PC3 316, HA1E 284, HCC515 284, HEPG2 284, HT29 284, MCF7 284, VCAP 237, A549 226. shRNA Lambert counts are highest in PC3 (555), HA1E (531), and A375 (530). MCH58 signatures are controls plus EIF2AK3, LRPPRC, TCF7L2, WFS1, and tunicamycin. Only TCF7L2 is on Lambert v1.01. Phase II cell_info has the same 98 cell ids. Phase II signatures were not counted. HEK293T is embryonal kidney, not a fibroblast screen.

Joung detail. The GEO design says 90 ORFs. The barcode map collapses to 67 symbols, 46 of them Lambert. GSE217215 (198 ORFs, SHARE-seq in hESCs) was opened; its ORF list was not intersected. GSE218789 (3,548 ORFs) is flow cytometry, not a transcriptome.

Norman detail. The series design says overexpression. GSM3906020 says K562 cells expressing dCas9-KRAB. The tier is the cancer line either way. The identity file yields 44 Lambert symbols, and NegCtrl barcodes are present.

## Secondary table: partial reprogramming in human cells

These are not TF screens. They are listed because a multi-donor partial-reprogramming transcriptome could repeat a young-versus-same comparison.

| Dataset | Paper | Accession | Donors and ages | Factors | Readout | URL opened |
| --- | --- | --- | --- | --- | --- | --- |
| OSKM Sendai, young and aged fibroblasts | Lu 2025, Cell | GSE297234 | 2 donors: GM23815 age 22, GM00731 age 96 | OSKM | 10x scRNA | GEO quick SOFT |
| Maturation-phase transient reprogramming | Gill 2022, eLife. GEO lists pubmed 5390271; a PubMed title search returns 35390271 | GSE165177 | 3 fibroblast donors: age 38, 53, and 53 | lentiviral hOKMS | bulk RNA-seq | GEO quick SOFT and GSM brief |
| Chemical reprogramming | Schoenfeldt 2025, EMBO Mol Med | GSE297984 | age labels 22, 23, 56, 83; not a proven unique-donor count | chemical 2c and 7c, not TFs | bulk RNA-seq | GEO quick SOFT and GSM brief |

GSE226189 was opened and left out of this table. It is RNA-seq of primary skin fibroblasts from 82 donors ages 22 to 89, with no perturbation.

Gill is the partial-reprogramming set with the most fibroblast donors among the pages opened here (three donors, ages 38, 53, and 53), with negative-control samples in the same series.

## UNVERIFIED leads

These were not given main-table counts.

- Replogle cell lines. The opened BioProject summary does not name K562 or RPE1. A remembered assignment of those lines, and any Lambert-TF count, stays here.
- A genome-wide ENCODE TF-knockdown catalog in K562 and HepG2. The ENCODE search endpoint returned a bot check, so that catalog was not opened. Only GSE33816 and GSE88370 are in the main table.
- GSE297233 (a possible OSK bulk RNA-seq series) was not opened.
- Mouse partial-reprogramming datasets were not opened.
- The Addgene MORF plasmid library is a reagent, not an expression matrix. Joung's 3,548-ORF experiment on GEO is the flow series GSE218789, already noted above.

## Ruler coverage

Frozen ruler: `results/fibro/frozen_ruler_ridge_raw.npz`, weights `w` and symbols `symbol`. Not refit. Duplicate ruler symbols were summed on absolute weight, leaving 23,467 symbols. Coverage is in `results/survey/coverage.csv`.

| Dataset | Tier | Measured symbols used | Ruler genes hit | Fraction of ruler genes | Fraction of total absolute weight | Note |
| --- | --- | --- | --- | --- | --- | --- |
| Southard Hs27 | A | 36,591 unique Gene Expression symbols from lane L4A | 16,771 / 23,467 | 0.715 | 0.630 | Measured genes only. Guide-capture features excluded. The 4,914-gene mean-population matrix was not used. |
| Sengstack | A | not read | coverage not computed | coverage not computed | coverage not computed | Gene symbols were not extracted. The h5ad is 1.43 GB. The shape attribute is 36,514 x 36,601. |

No other dataset is Tier A or Tier B, so no other coverage row was computed.

## Fired key

`screen_ready_fibroblast`. An in-silico fibroblast TF screen is feasible with existing data. The datasets are Southard et al. 2025 Hs27 CRISPRa Perturb-seq (Zenodo 10.5281/zenodo.15200179; 1,488 Lambert TFs; 36,601 measured genes; ruler coverage 0.715 of genes and 0.630 of absolute weight) and Sengstack et al. 2026 fibroblast CRISPRa/CRISPRi Perturb-seq (figshare 10.6084/m9.figshare.30898748; 179 Lambert TFs overexpressed; matrix shape 36,514 by 36,601).

Neither Tier A page states donor age in years. That does not change the tier. It does limit a direct young-versus-old contrast inside the screen itself.

## Limitations

- The search was not exhaustive.
- Tiers depend on what each data page, and the files that were actually read, report.
- The ENCODE portal search was blocked by a bot check, so a larger knockdown catalog may exist and was not scored.
- L1000 inferred genes were not counted as measured.
- Files larger than 5 GB were not downloaded. Hs27's measured-gene count comes from one Cell Ranger lane matrix pulled out of the zip, not from the 8.05 GB singlets file.
- Sengstack gene symbols were not extracted, so its ruler coverage is blank even though the shape attribute clears the 10,000-gene bar.
- Lambert v1.01 has 1,639 names. Genes off that list were not counted, including 348 Hs27 targets and 21 Sengstack targets.
- RPE-1 numbers were not copied from the Hs27 files.

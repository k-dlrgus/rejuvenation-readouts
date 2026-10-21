# SRSFINV - data inventory, GEO series GSE325735

Inventory only. Nothing here is a pre-registered result and nothing is
interpreted. No statistic, distance or fitted axis was computed.
Ambiguities are recorded as ambiguities, not resolved. Every count below
is a count of GEO sample records (GSMs) unless it says otherwise.

## 0. Provenance

Files read, all untracked under `data/raw/srsfinv/` (`.gitignore` rule `data/raw/`):

| file | source | bytes | md5 | fetched this run |
|---|---|---|---|---|
| `GSE325735_family.soft.gz` | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE325nnn/GSE325735/soft/GSE325735_family.soft.gz | 49,917 | `6257c37dbec6f49a5cff6f110ca641a4` | False |
| `sra_runinfo.csv` | https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA1441554 | 1,518,357 | `e4e068f26c23c3c19af21340cad3c0d8` | False |
| `human_counts_head.gz.part` | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE325nnn/GSE325735/suppl/GSE325735_human_raw_counts.txt.gz (bytes 0-262143) | 262,144 | `128cbcf604f90530bfd9415ffad3d59f` | False |

Supplementary files on the GEO FTP (sizes from HTTP HEAD, not downloaded):

- `GSE325735_human_raw_counts.txt.gz` 26,007,361 B, Last-Modified Mon, 23 Mar 2026 16:14:19 GMT
- `GSE325735_mouse_raw_counts.txt.gz` 1,773,429 B, Last-Modified Mon, 23 Mar 2026 16:14:20 GMT

Download status, stated plainly: **no expression matrix was downloaded.**
The first 262,144 bytes of the human counts file were read with an HTTP
range request (1.0% of the file) to see its header and identifier
type. The mouse counts file was not opened. The SRA run table is
metadata, fetched for per-sample read depth, which GEO does not carry.

Series record fields:

- `Series_title`: Transcriptomic reprogramming screen identifies SRSF1 as a regulator of aging
- `Series_status`: Public on Aug 15 2026
- `Series_submission_date`: Mar 23 2026
- `Series_last_update_date`: Oct 14 2026
- `Series_web_link`: https://www.nature.com/articles/s43587-026-01234-8
- `Series_type`: Expression profiling by high throughput sequencing
- `Series_overall_design`: RNA-seq profiling of wild type normal human dermal fibroblast cells or iPSCs and their derivatives overexpressing different transgenes.
- `Series_relation`: BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1441554
- `Series_sample_id`: 839 GSMs; family SOFT holds 839 `^SAMPLE` blocks
- `Series_supplementary_file`: ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE325nnn/GSE325735/suppl/GSE325735_human_raw_counts.txt.gz; ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE325nnn/GSE325735/suppl/GSE325735_mouse_raw_counts.txt.gz

The GEO title today is "Transcriptomic reprogramming screen identifies SRSF1 as a regulator of aging". The task names the series
"Transcriptomic reprogramming screen identifies SRSF1 as a central regulator of aging and lifespan". They differ; the GEO record
was last updated Oct 14 2026, after `SCOPING_TISSUE.md` (2026-09-17).
This file does not know which title was current on which date.

## 0b. Comparison with the stored scoping parse

`results/scoping/gse325735_parse.json` was built from GSM titles by
`results/scoping/_parse_325735.py`. The same title rules, re-applied here to
the titles in the family SOFT:

- quantities compared: 124; exact matches: 124; mismatches: 0
- full comparison: `results/srsfinv/scoping_parse_comparison.csv`

The scoping parse reads titles. Sections 2-4 below read the structured
`!Sample_characteristics_ch1` fields, which give different counts for
some labels (WT, NTg, OSKM, SRSF1); section 3c tabulates the difference.

## 1. Assay per sample

Fields that say what kind of assay each GSM is, all 839 GSMs:

- `sample_type`: SRA (839)
- `library_strategy`: RNA-Seq (839)
- `library_source`: transcriptomic (839)
- `library_selection`: cDNA (839)

| organism | platform | instrument | molecule | n GSMs | platform title |
|---|---|---|---|---|---|
| Homo sapiens | GPL18573 | NextSeq 500 | polyA RNA | 15 | Illumina NextSeq 500 (Homo sapiens) |
| Homo sapiens | GPL18573 | NextSeq 500 | total RNA | 40 | Illumina NextSeq 500 (Homo sapiens) |
| Homo sapiens | GPL24676 | Illumina NovaSeq 6000 | polyA RNA | 716 | Illumina NovaSeq 6000 (Homo sapiens) |
| Homo sapiens | GPL34284 | Illumina NovaSeq X Plus | polyA RNA | 36 | Illumina NovaSeq X Plus (Homo sapiens) |
| Mus musculus | GPL34290 | Illumina NovaSeq X Plus | polyA RNA | 32 | Illumina NovaSeq X Plus (Mus musculus) |

`total RNA` GSMs (40): batches LIB049469 16, LIB047784 8, LIB048130 8, LIB050518 8; groups Homo sapiens 40; ages 29/31/36/41/52/55/60/66/68/69/73/79/86/92/94.
The extract protocol below names a polyA mRNA workflow for the NEBNext libraries
and does not mention total RNA (ambiguity A12).

Whole-word search of the family SOFT text for single-cell terms (`single-cell`,
`single cell`, `scRNA`, `10x`, `Chromium`, `droplet`, `snRNA`, `cell barcode`, `UMI`):
no hits.
`Series_type` is "Expression profiling by high throughput sequencing". Every GSM is `RNA-Seq`, `cDNA`, one
library per GSM, and the processing text below describes STAR + per-gene
read counting of bulk libraries. **No single-cell assay is recorded.**
Counts below are therefore not split by assay; they are split by organism.

`!Sample_data_processing`, identical in all 839 GSMs: True

> To minimize variability arising from differences in sequencing depth, FASTQ files for all libraries (except LIB058979, LIB092645, and LIB083102) were downsampled to 20 million reads using Seqtk-1.3.
> Raw RNA-seq reads were aligned to the GRCh38/GRCm39 reference genome using STAR v2.5.2b and checked for quality control using FastQC v0.11.5
> Alignment files were indexed using SAMtools v1.3.154 and gene-level counts were obtained using featureCounts from the Subread v2.0.1 package
> Assembly: GRCh38, GRCm39
> Supplementary files format and content: tab-delimited text file includes raw counts for each sample

`!Sample_extract_protocol_ch1` (identical in all GSMs: True):

> Cells were lysed using the Monarch DNA/RNA Protection Reagent (NEB T2011) and stored at -80°C until RNA extraction. Total RNA was extrqacted using the Monarch Total RNA Miniprep kit (NEB T2010). Mouse tissues were processed by Novogene.
> Libraries were generated using the NEBNext Ultra II Directional RNA Library Prep Kit for Illumina using the polyA mRNA workflow according to manufacturer's instructions. LIB058979 was genereated using the Tecan MagicPrepRevelo mRNA-seq kit acoording to manufacturare's instructions. LIB083102 was generated by Novogene using a directional mRNA library preparation method.

### 1b. Read depth

GEO carries no per-sample read-count field (`!Sample_data_row_count` is
0 in every GSM). Depth below is from the SRA run table.
- SRA runs: 2871; distinct `SampleName` (= GSM): 839; GSMs with no run: 0; runs whose `SampleName` is not a GSM of this series: 0
- `LibraryLayout`: PAIRED (2871 runs)
- `Sex` column in the SRA table: '<empty>' (2871)
- per-GSM depth table: `results/srsfinv/sra_depth_per_gsm.csv`. `spots` are
  read pairs for PAIRED runs, summed over a GSM's runs. These are depths
  **before** the downsampling the processing text describes.

| organism platform | count | min | 5% | 50% | 95% | max |
|---|---|---|---|---|---|---|
| Homo sapiens GPL18573 | 55 | 6,164,291 | 9,045,880 | 15,326,488 | 35,921,447 | 40,118,491 |
| Homo sapiens GPL24676 | 716 | 10,367,411 | 14,303,033 | 21,041,086 | 35,870,382 | 115,185,263 |
| Homo sapiens GPL34284 | 36 | 30,190,589 | 30,592,903 | 34,826,980 | 41,585,725 | 44,632,112 |
| Mus musculus GPL34290 | 32 | 20,393,923 | 21,050,492 | 24,810,818 | 40,835,552 | 51,653,529 |

GSMs with fewer than 20,000,000 spots in SRA: 368 (Homo sapiens 368).

GSMs in the three libraries the processing text exempts from 20 M downsampling
(`ch_batch` value): LIB058979 96, LIB092645 36, LIB083102 32; all other GSMs 675.

`!Sample_description` free-text flags beyond donor age and passage (human GSMs):

| description_flags | n GSMs | batches |
|---|---|---|
| <none> | 707 | LIB047784, LIB048130, LIB049469, LIB050518, LIB051945, LIB052643, LIB052667, LIB053449, LIB053907, LIB054839, LIB092645 |
| Batch 9 | 8 | LIB058979 |
| Batch 1 | 8 | LIB058979 |
| Batch 10 | 8 | LIB058979 |
| Batch 2 | 8 | LIB058979 |
| Batch 7 | 8 | LIB058979 |
| Batch 6 | 8 | LIB058979 |
| Batch 3 | 8 | LIB058979 |
| Batch 4 | 8 | LIB058979 |
| Batch 11 | 8 | LIB058979 |
| Batch 12 | 8 | LIB058979 |
| Batch 5 | 8 | LIB058979 |
| Batch 8 | 8 | LIB058979 |
| <10M reads | 3 | LIB052667 |
| High GC | 1 | LIB052667 |

- `<10M reads`: GSM9611695 (M65 line with Dox-inducible OGT transgene, Dox_1000ng/mL, P10, Rep2, 14,887,523 spots), GSM9611749 (M65 line with Dox-inducible UBA52 transgene, Dox_1000ng/mL, P10, Rep2, 12,232,670 spots), GSM9611766 (M79 line with Dox-inducible ACTG1 transgene, Dox_1000ng/mL, P13, Rep1, 16,643,422 spots)
- `High GC`: GSM9611762 (M65 line with Dox-inducible mTagBFP2 transgene, Dox_1000ng/mL, P10, Rep1, 24,801,777 spots)
- titles containing `redo`: GSM9611954 (M65 line with Dox-inducible mTagBFP2 transgene, Dox_1000ng/mL, P10, Rep1, redo), GSM9611955 (M65 line with Dox-inducible MDM2 transgene, Dox_1000ng/mL, P10, Rep1, redo), GSM9611956 (M65 line with Dox-inducible UBA52 transgene, Dox_1000ng/mL, P10, Rep2, redo)

## 2. Donor table (human)

Where each field comes from:

- donor line: `!Sample_characteristics_ch1` `cell line:` (equals
  `!Sample_source_name_ch1` in 807 of 807 human GSMs)
- donor age: `!Sample_description` text `Donor age N`. Not in
  `!Sample_characteristics_ch1` (no `age` key exists there).
- title M-code: the `M<n>` token in `!Sample_title`
- passage: `!Sample_description` `Passage N`; titles also carry a `P<n>` token
- sex: **not stated** anywhere in the GEO record (no `sex`/`gender` key in
  characteristics; see 1b for the SRA `Sex` column). Not looked up elsewhere.

Characteristic keys present, with GSM counts: `genotype` 839, `treatment` 839, `batch` 839, `cell line` 807, `cell type` 807, `tissue` 32
Whole-SOFT search for `sex|gender|male|female`: no hits.

Age provenance per human GSM:

- age stated in `!Sample_description`: 807 of 807
- title M-code present: 799; of these M-code equals described age: 799; disagrees: 0
- age stated in description but no M-code in title: 8 (lines B6, F3)
- age only inferable from title (no described age): 0
- no age from either source: 0

So no human GSM's age is inferred from the title alone within this record.
What the record does not say is the meaning of the `M` prefix; it equals
the described age wherever both exist (see ambiguity A2).

Distinct `cell line` values among human GSMs: 19. With at least one `Fibroblast`
GSM: 17. With at least one `iPSC` GSM: 3 (F3, B6, AG07141; AG07141 has both).
Distinct described ages among human fibroblast GSMs: 16. (`SCOPING_TISSUE.md`
wrote "16 M-coded NHDF lines"; 16 is the number of distinct M-codes, i.e. ages,
not of lines.)

| cell line | cell type | n Fibroblast GSMs | n iPSC GSMs | described age(s) | title M-code(s) | described passage(s) | title P token(s) | sex | n GSMs | n batches |
|---|---|---|---|---|---|---|---|---|---|---|
| AG04054B | Fibroblast | 3 | 0 | 29 | 29 | 3/4 | - | not stated | 3 | 2 |
| F3 | iPSC | 0 | 5 | 31 | - | - | - | not stated | 5 | 2 |
| B6 | iPSC | 0 | 3 | 36 | - | - | P+10/P+2/P0 | not stated | 3 | 1 |
| AG11557A | Fibroblast | 4 | 0 | 36 | 36 | 5 | - | not stated | 4 | 1 |
| AG08906 | Fibroblast | 1 | 0 | 41 | 41 | 5 | - | not stated | 1 | 1 |
| AG05186 | Fibroblast | 1 | 0 | 52 | 52 | 4 | - | not stated | 1 | 1 |
| AG05192B | Fibroblast | 4 | 0 | 55 | 55 | 5 | - | not stated | 4 | 1 |
| AG12657 | Fibroblast | 273 | 0 | 55 | 55 | 6/9/10/11/12/13/14 | P10/P11/P12/P13/P14/P9 | not stated | 273 | 7 |
| AG10883 | Fibroblast | 1 | 0 | 60 | 60 | 3 | - | not stated | 1 | 1 |
| AG06294A | Fibroblast | 147 | 0 | 65 | 65 | 5/7/8/9/10/11/12/13 | P10/P11/P12/P13/P5/P7/P8/P9 | not stated | 147 | 5 |
| AG07141 | Fibroblast/iPSC | 1 | 1 | 66 | 66 | 2/4 | - | not stated | 2 | 1 |
| AG09697 | Fibroblast | 22 | 0 | 67 | 67 | 5/8 | P5/P8 | not stated | 22 | 1 |
| AG13369 | Fibroblast | 35 | 0 | 68 | 68 | 9/11/12 | P11/P12/P9 | not stated | 35 | 3 |
| AG12587 | Fibroblast | 23 | 0 | 69 | 69 | 8/9/12 | P12/P9 | not stated | 23 | 2 |
| AG05271B | Fibroblast | 4 | 0 | 73 | 73 | 16 | - | not stated | 4 | 1 |
| AG16359 | Fibroblast | 271 | 0 | 79 | 79 | 7/9/10/11/12/13/14/15 | P10/P11/P12/P13/P14/P15/P9 | not stated | 271 | 7 |
| AG13349A | Fibroblast | 1 | 0 | 86 | 86 | 12 | - | not stated | 1 | 1 |
| AG04064B | Fibroblast | 4 | 0 | 92 | 92 | 5 | - | not stated | 4 | 1 |
| AG08433C | Fibroblast | 3 | 0 | 94 | 94 | 3/4 | - | not stated | 3 | 2 |

The 9 iPSC GSMs, verbatim titles: GSM9611426 `F3 line with Dox-inducible GFP transgene, Dox_250ng/mL` (line F3, described age 31); GSM9611427 `F3 line with Dox-inducible GFP transgene, Dox_2000ng/mL` (line F3, described age 31); GSM9611428 `F3 line with Dox-inducible OTX2 transgene` (line F3, described age 31); GSM9611429 `F3 line with Dox-inducible ZGLP1 transgene` (line F3, described age 31); GSM9611430 `B6 line, P0` (line B6, described age 36); GSM9611431 `B6 line, P+2` (line B6, described age 36); GSM9611432 `B6 line, P+10` (line B6, described age 36); GSM9611433 `F3 line with Dox-inducible ZGLP1 transgene, Dox_2000ng/mL` (line F3, described age 31); GSM9611461 `M66 WT iPSC line (isogenic with S39)` (line AG07141, described age 66).
Titles containing `isogenic` (8): GSM9611454 `M41 WT line (isogenic with S34)` (library S33); GSM9611455 `M55 WT line (isogenic with S33)` (library S34); GSM9611456 `M52 WT line (isogenic with S36)` (library S35); GSM9611457 `M68 WT line (isogenic with S35)` (library S36); GSM9611458 `M60 WT line (isogenic with S38)` (library S37); GSM9611459 `M79 WT line (isogenic with S37)` (library S38); GSM9611460 `M66 WT line (isogenic with S40)` (library S39); GSM9611461 `M66 WT iPSC line (isogenic with S39)` (library S40).

Ages carried by more than one `cell line` value: 36: AG11557A, B6; 55: AG05192B, AG12657.
Cell lines with more than one described age: none.
Whether two `cell line` values are the same human donor is not stated; the
table counts `cell line` values, not verified donors.

## 3. Perturbation table

Source: `!Sample_characteristics_ch1` `genotype:`. Every value is listed as
written. Which ones are controls is section 4.

- distinct `genotype` values, all GSMs: 102
- distinct `genotype` values in human fibroblasts: 99; in human iPSC: 4; in mouse: 2
- human-fibroblast `genotype` values other than `WT`, `NTg`, `mTagBFP2`, `GFP`: 96

| genotype | human fibroblast | human iPSC | mouse | total | n fibroblast donor ages | n fibroblast cell lines |
|---|---|---|---|---|---|---|
| SRSF1 | 72 | 0 | 16 | 88 | 6 | 6 |
| mTagBFP2 | 73 | 0 | 0 | 73 | 6 | 6 |
| WT | 53 | 4 | 0 | 57 | 16 | 17 |
| OSKM | 38 | 0 | 0 | 38 | 6 | 6 |
| GFP | 0 | 2 | 16 | 18 | 0 | 0 |
| SLC2A13 | 18 | 0 | 0 | 18 | 6 | 6 |
| SIRT1 | 18 | 0 | 0 | 18 | 6 | 6 |
| PPP1CB | 18 | 0 | 0 | 18 | 6 | 6 |
| OGT | 18 | 0 | 0 | 18 | 6 | 6 |
| WDTC1 | 18 | 0 | 0 | 18 | 6 | 6 |
| RNASEL | 16 | 0 | 0 | 16 | 6 | 6 |
| HDAC1 | 16 | 0 | 0 | 16 | 6 | 6 |
| SUMO1 | 10 | 0 | 0 | 10 | 3 | 3 |
| NTg | 8 | 0 | 0 | 8 | 4 | 4 |
| MDM2 | 7 | 0 | 0 | 7 | 3 | 3 |
| UBA52 | 7 | 0 | 0 | 7 | 3 | 3 |
| HES1 | 6 | 0 | 0 | 6 | 3 | 3 |
| GSK3B | 6 | 0 | 0 | 6 | 3 | 3 |
| PPP4C | 6 | 0 | 0 | 6 | 3 | 3 |
| PSMC6 | 6 | 0 | 0 | 6 | 3 | 3 |
| PTEN | 6 | 0 | 0 | 6 | 3 | 3 |
| RAD23B | 6 | 0 | 0 | 6 | 3 | 3 |
| MAPK7 | 6 | 0 | 0 | 6 | 3 | 3 |
| P4HB | 6 | 0 | 0 | 6 | 3 | 3 |
| CDC20 | 6 | 0 | 0 | 6 | 3 | 3 |
| CDK1 | 6 | 0 | 0 | 6 | 3 | 3 |
| EEF2 | 6 | 0 | 0 | 6 | 3 | 3 |
| ESR1 | 6 | 0 | 0 | 6 | 3 | 3 |
| CHEK2 | 6 | 0 | 0 | 6 | 3 | 3 |
| CDKN2A | 6 | 0 | 0 | 6 | 3 | 3 |
| CYCS | 6 | 0 | 0 | 6 | 3 | 3 |
| CXCR2 | 6 | 0 | 0 | 6 | 3 | 3 |
| KAT7 | 6 | 0 | 0 | 6 | 3 | 3 |
| ISG15 | 6 | 0 | 0 | 6 | 3 | 3 |
| CCNB1 | 6 | 0 | 0 | 6 | 3 | 3 |
| CCNA2 | 6 | 0 | 0 | 6 | 3 | 3 |
| ACTL6A | 6 | 0 | 0 | 6 | 3 | 3 |
| ATG5 | 6 | 0 | 0 | 6 | 3 | 3 |
| YWHAQ | 6 | 0 | 0 | 6 | 3 | 3 |
| YWHAB | 6 | 0 | 0 | 6 | 3 | 3 |
| UBE2I | 6 | 0 | 0 | 6 | 3 | 3 |
| UBE2S | 6 | 0 | 0 | 6 | 3 | 3 |
| SLC2A10 | 6 | 0 | 0 | 6 | 3 | 3 |
| SOD2 | 6 | 0 | 0 | 6 | 3 | 3 |
| RPLP0 | 6 | 0 | 0 | 6 | 3 | 3 |
| RPS27A | 6 | 0 | 0 | 6 | 3 | 3 |
| UBQLN1 | 6 | 0 | 0 | 6 | 3 | 3 |
| PPP3CA | 6 | 0 | 0 | 6 | 3 | 3 |
| PPP1CA | 6 | 0 | 0 | 6 | 3 | 3 |
| PPP1CC | 6 | 0 | 0 | 6 | 3 | 3 |
| PIK3CA | 6 | 0 | 0 | 6 | 3 | 3 |
| PPIA | 6 | 0 | 0 | 6 | 3 | 3 |
| NFKB1 | 6 | 0 | 0 | 6 | 3 | 3 |
| PIK3C2A | 6 | 0 | 0 | 6 | 3 | 3 |
| PARK7 | 6 | 0 | 0 | 6 | 3 | 3 |
| UBB | 6 | 0 | 0 | 6 | 3 | 3 |
| SOD1 | 6 | 0 | 0 | 6 | 3 | 3 |
| TXN | 6 | 0 | 0 | 6 | 3 | 3 |
| PKM | 6 | 0 | 0 | 6 | 3 | 3 |
| CHEK1 | 4 | 0 | 0 | 4 | 2 | 2 |
| ACTG1 | 4 | 0 | 0 | 4 | 2 | 2 |
| AKT1 | 4 | 0 | 0 | 4 | 2 | 2 |
| CDK4 | 4 | 0 | 0 | 4 | 2 | 2 |
| GART | 4 | 0 | 0 | 4 | 2 | 2 |
| CDC34 | 4 | 0 | 0 | 4 | 2 | 2 |
| CDK2 | 4 | 0 | 0 | 4 | 2 | 2 |
| DDX3X | 4 | 0 | 0 | 4 | 2 | 2 |
| CTNNB1 | 4 | 0 | 0 | 4 | 2 | 2 |
| FAU | 4 | 0 | 0 | 4 | 2 | 2 |
| HSP90AB1 | 4 | 0 | 0 | 4 | 2 | 2 |
| XPO1 | 4 | 0 | 0 | 4 | 2 | 2 |
| LRRK2 | 4 | 0 | 0 | 4 | 2 | 2 |
| HSPA8 | 4 | 0 | 0 | 4 | 2 | 2 |
| KDM6A | 4 | 0 | 0 | 4 | 2 | 2 |
| HSPA5 | 4 | 0 | 0 | 4 | 2 | 2 |
| HSP90B1 | 4 | 0 | 0 | 4 | 2 | 2 |
| HNRNPA2B1 | 4 | 0 | 0 | 4 | 2 | 2 |
| HSPD1 | 4 | 0 | 0 | 4 | 2 | 2 |
| HDAC6 | 4 | 0 | 0 | 4 | 2 | 2 |
| NOTCH1 | 4 | 0 | 0 | 4 | 2 | 2 |
| MCM7 | 4 | 0 | 0 | 4 | 2 | 2 |
| MAPK14 | 4 | 0 | 0 | 4 | 2 | 2 |
| NPM1 | 4 | 0 | 0 | 4 | 2 | 2 |
| PSMD4 | 4 | 0 | 0 | 4 | 2 | 2 |
| PLK1 | 4 | 0 | 0 | 4 | 2 | 2 |
| PPP2CA | 4 | 0 | 0 | 4 | 2 | 2 |
| RUVBL2 | 4 | 0 | 0 | 4 | 2 | 2 |
| VCP | 4 | 0 | 0 | 4 | 2 | 2 |
| UBE2N | 4 | 0 | 0 | 4 | 2 | 2 |
| RAD51 | 4 | 0 | 0 | 4 | 2 | 2 |
| UBE2A | 4 | 0 | 0 | 4 | 2 | 2 |
| RRM1 | 4 | 0 | 0 | 4 | 2 | 2 |
| TP53 (DN) | 4 | 0 | 0 | 4 | 2 | 2 |
| DICER1 | 2 | 0 | 0 | 2 | 1 | 1 |
| TP53 | 2 | 0 | 0 | 2 | 1 | 1 |
| PIK3C3 | 2 | 0 | 0 | 2 | 1 | 1 |
| MTOR | 2 | 0 | 0 | 2 | 1 | 1 |
| SMARCA4 | 2 | 0 | 0 | 2 | 1 | 1 |
| UBXN7 | 2 | 0 | 0 | 2 | 1 | 1 |
| SOX2 | 2 | 0 | 0 | 2 | 1 | 1 |
| ZGLP1 | 0 | 2 | 0 | 2 | 0 | 0 |
| OTX2 | 0 | 1 | 0 | 1 | 0 | 0 |

### 3b. Genotype x described donor age, human fibroblasts

|  | age 29 | age 36 | age 41 | age 52 | age 55 | age 60 | age 65 | age 66 | age 67 | age 68 | age 69 | age 73 | age 79 | age 86 | age 92 | age 94 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SRSF1 | 0 | 0 | 0 | 0 | 28 | 0 | 4 | 0 | 2 | 8 | 2 | 0 | 28 | 0 | 0 | 0 |
| mTagBFP2 | 0 | 0 | 0 | 0 | 26 | 0 | 5 | 0 | 2 | 8 | 2 | 0 | 30 | 0 | 0 | 0 |
| WT | 3 | 4 | 1 | 1 | 11 | 1 | 6 | 1 | 2 | 1 | 3 | 4 | 7 | 1 | 4 | 3 |
| OSKM | 0 | 0 | 0 | 0 | 18 | 0 | 2 | 0 | 2 | 2 | 2 | 0 | 12 | 0 | 0 | 0 |
| SLC2A13 | 0 | 0 | 0 | 0 | 4 | 0 | 4 | 0 | 2 | 2 | 2 | 0 | 4 | 0 | 0 | 0 |
| SIRT1 | 0 | 0 | 0 | 0 | 4 | 0 | 4 | 0 | 2 | 2 | 2 | 0 | 4 | 0 | 0 | 0 |
| PPP1CB | 0 | 0 | 0 | 0 | 4 | 0 | 4 | 0 | 2 | 2 | 2 | 0 | 4 | 0 | 0 | 0 |
| OGT | 0 | 0 | 0 | 0 | 4 | 0 | 4 | 0 | 2 | 2 | 2 | 0 | 4 | 0 | 0 | 0 |
| WDTC1 | 0 | 0 | 0 | 0 | 4 | 0 | 4 | 0 | 2 | 2 | 2 | 0 | 4 | 0 | 0 | 0 |
| RNASEL | 0 | 0 | 0 | 0 | 4 | 0 | 2 | 0 | 2 | 2 | 2 | 0 | 4 | 0 | 0 | 0 |
| HDAC1 | 0 | 0 | 0 | 0 | 4 | 0 | 2 | 0 | 2 | 2 | 2 | 0 | 4 | 0 | 0 | 0 |
| SUMO1 | 0 | 0 | 0 | 0 | 2 | 0 | 6 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| NTg | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 2 | 0 | 0 | 2 | 0 | 0 | 0 |
| MDM2 | 0 | 0 | 0 | 0 | 2 | 0 | 3 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| UBA52 | 0 | 0 | 0 | 0 | 2 | 0 | 3 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| HES1 | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| GSK3B | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| PPP4C | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| PSMC6 | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| PTEN | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| RAD23B | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| MAPK7 | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| P4HB | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| CDC20 | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| CDK1 | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| EEF2 | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| ESR1 | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| CHEK2 | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| CDKN2A | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| CYCS | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| CXCR2 | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| KAT7 | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| ISG15 | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| CCNB1 | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| CCNA2 | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| ACTL6A | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| ATG5 | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| YWHAQ | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| YWHAB | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| UBE2I | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| UBE2S | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| SLC2A10 | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| SOD2 | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| RPLP0 | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| RPS27A | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| UBQLN1 | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| PPP3CA | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| PPP1CA | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| PPP1CC | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| PIK3CA | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| PPIA | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| NFKB1 | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| PIK3C2A | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| PARK7 | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| UBB | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| SOD1 | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| TXN | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| PKM | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| CHEK1 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| ACTG1 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| AKT1 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| CDK4 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| GART | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| CDC34 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| CDK2 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| DDX3X | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| CTNNB1 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| FAU | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| HSP90AB1 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| XPO1 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| LRRK2 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| HSPA8 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| KDM6A | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| HSPA5 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| HSP90B1 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| HNRNPA2B1 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| HSPD1 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| HDAC6 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| NOTCH1 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| MCM7 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| MAPK14 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| NPM1 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| PSMD4 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| PLK1 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| PPP2CA | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| RUVBL2 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| VCP | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| UBE2N | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| RAD51 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| UBE2A | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| RRM1 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| TP53 (DN) | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| DICER1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| TP53 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| PIK3C3 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| MTOR | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| SMARCA4 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| UBXN7 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| SOX2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |

**Genotypes applied to more than one donor age in human fibroblasts (92):**

- `SRSF1`: ages 55, 65, 67, 68, 69, 79
- `mTagBFP2`: ages 55, 65, 67, 68, 69, 79
- `WT`: ages 29, 36, 41, 52, 55, 60, 65, 66, 67, 68, 69, 73, 79, 86, 92, 94
- `OSKM`: ages 55, 65, 67, 68, 69, 79
- `SLC2A13`: ages 55, 65, 67, 68, 69, 79
- `SIRT1`: ages 55, 65, 67, 68, 69, 79
- `PPP1CB`: ages 55, 65, 67, 68, 69, 79
- `OGT`: ages 55, 65, 67, 68, 69, 79
- `WDTC1`: ages 55, 65, 67, 68, 69, 79
- `RNASEL`: ages 55, 65, 67, 68, 69, 79
- `HDAC1`: ages 55, 65, 67, 68, 69, 79
- `SUMO1`: ages 55, 65, 79
- `NTg`: ages 55, 65, 68, 79
- `MDM2`: ages 55, 65, 79
- `UBA52`: ages 55, 65, 79
- `HES1`: ages 55, 65, 79
- `GSK3B`: ages 55, 65, 79
- `PPP4C`: ages 55, 65, 79
- `PSMC6`: ages 55, 65, 79
- `PTEN`: ages 55, 65, 79
- `RAD23B`: ages 55, 65, 79
- `MAPK7`: ages 55, 65, 79
- `P4HB`: ages 55, 65, 79
- `CDC20`: ages 55, 65, 79
- `CDK1`: ages 55, 65, 79
- `EEF2`: ages 55, 65, 79
- `ESR1`: ages 55, 65, 79
- `CHEK2`: ages 55, 65, 79
- `CDKN2A`: ages 55, 65, 79
- `CYCS`: ages 55, 65, 79
- `CXCR2`: ages 55, 65, 79
- `KAT7`: ages 55, 65, 79
- `ISG15`: ages 55, 65, 79
- `CCNB1`: ages 55, 65, 79
- `CCNA2`: ages 55, 65, 79
- `ACTL6A`: ages 55, 65, 79
- `ATG5`: ages 55, 65, 79
- `YWHAQ`: ages 55, 65, 79
- `YWHAB`: ages 55, 65, 79
- `UBE2I`: ages 55, 65, 79
- `UBE2S`: ages 55, 65, 79
- `SLC2A10`: ages 55, 65, 79
- `SOD2`: ages 55, 65, 79
- `RPLP0`: ages 55, 65, 79
- `RPS27A`: ages 55, 65, 79
- `UBQLN1`: ages 55, 65, 79
- `PPP3CA`: ages 55, 65, 79
- `PPP1CA`: ages 55, 65, 79
- `PPP1CC`: ages 55, 65, 79
- `PIK3CA`: ages 55, 65, 79
- `PPIA`: ages 55, 65, 79
- `NFKB1`: ages 55, 65, 79
- `PIK3C2A`: ages 55, 65, 79
- `PARK7`: ages 55, 65, 79
- `UBB`: ages 55, 65, 79
- `SOD1`: ages 55, 65, 79
- `TXN`: ages 55, 65, 79
- `PKM`: ages 55, 65, 79
- `CHEK1`: ages 55, 79
- `ACTG1`: ages 55, 79
- `AKT1`: ages 55, 79
- `CDK4`: ages 55, 79
- `GART`: ages 55, 79
- `CDC34`: ages 55, 79
- `CDK2`: ages 55, 79
- `DDX3X`: ages 55, 79
- `CTNNB1`: ages 55, 79
- `FAU`: ages 55, 79
- `HSP90AB1`: ages 55, 79
- `XPO1`: ages 55, 79
- `LRRK2`: ages 65, 79
- `HSPA8`: ages 55, 79
- `KDM6A`: ages 55, 79
- `HSPA5`: ages 55, 79
- `HSP90B1`: ages 55, 79
- `HNRNPA2B1`: ages 55, 79
- `HSPD1`: ages 55, 79
- `HDAC6`: ages 55, 79
- `NOTCH1`: ages 55, 79
- `MCM7`: ages 55, 79
- `MAPK14`: ages 55, 79
- `NPM1`: ages 55, 79
- `PSMD4`: ages 55, 79
- `PLK1`: ages 55, 79
- `PPP2CA`: ages 55, 79
- `RUVBL2`: ages 55, 79
- `VCP`: ages 55, 79
- `UBE2N`: ages 55, 79
- `RAD51`: ages 55, 79
- `UBE2A`: ages 55, 79
- `RRM1`: ages 55, 79
- `TP53 (DN)`: ages 55, 65

Of these, 89 are not `WT`, `NTg` or a reporter (`mTagBFP2`, `GFP`).
Human fibroblast GSMs carrying a genotype other than `WT`/`NTg`/reporter: 664,
from cell lines AG06294A, AG09697, AG12587, AG12657, AG13369, AG16359, at described ages 55, 65, 67, 68, 69, 79. **No such GSM comes from a donor
younger than 55.**

Genotypes seen at exactly one fibroblast donor age: 7, by that age: 55 (3), 65 (1), 79 (3).
Genotypes with no human fibroblast GSM: GFP, ZGLP1, OTX2.

### 3c. Title token versus `genotype` field

Title token = the scoping rule (`WT line` -> WT, else `Dox-inducible X transgene` -> X).
Human GSMs where it differs from `genotype`:

| title_token | ch_genotype | n GSMs |
|---|---|---|
| <none> | WT | 3 |
| NTg | WT | 22 |
| OSKM | SOX2 | 2 |
| TP53 (DN) | TP53 | 2 |

Examples (first GSM of each pair):

- GSM9611430 title `B6 line, P0`; genotype `WT`; treatment `NT`
- GSM9611654 title `M55 line with Dox-inducible NTg transgene, Dox_1000ng/mL, P9, Rep1`; genotype `WT`; treatment `Dox_1000ng/mL, 3days`
- GSM9611856 title `M79 line with Dox-inducible OSKM transgene, Dox_1000ng/mL, P13, Rep1`; genotype `SOX2`; treatment `Dox_1000ng/mL, 3days`
- GSM9611918 title `M79 line with Dox-inducible TP53 (DN) transgene, Dox_1000ng/mL, P13, Rep1`; genotype `TP53`; treatment `Dox_1000ng/mL, 3days`

## 4. Controls

The record never uses the words control, untreated or vector. The
candidate control conditions below are the stated field values; which of
them counts as "untreated" is not stated (ambiguity A4).

Values of `treatment:` among human GSMs: `Dox_1000ng/mL, 3days` 501, `Dox_2000ng/mL, 3days` 120, `NT` 39, `Dox_1000ng/mL, 3d-ON` 36, `Dox_1000ng/mL, 1d-ON` 12, `Dox_1000ng/mL, 1d-ON/8d-OFF` 12, `Dox_1000ng/mL, 1d-ON/16d-OFF` 12, `Dox_1000ng/mL, 2d-ON` 12, `Dox_1000ng/mL, 2d-ON/8d-OFF` 12, `Dox_1000ng/mL, 5d-ON/2d-OFF` 12, `Dox_1000ng/mL, 3d-ON/16d-OFF` 8, `Dox_1000ng/mL, 3d-ON/8d-OFF` 8, `Dox_1000ng/mL, 2d-ON/16d-OFF` 8, `Rapa_100ng/mL_24h` 4, `Rapa_100ng/mL_48h` 4, `Rapa_100ng/mL_72h` 4, `Rapa_100ng/mL, 3days` 2, `Dox_250ng/mL, 3days` 1

| candidate / group | age 29 | age 31 | age 36 | age 41 | age 52 | age 55 | age 60 | age 65 | age 66 | age 67 | age 68 | age 69 | age 73 | age 79 | age 86 | age 92 | age 94 | total |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| C1a genotype WT, treatment NT, title says WT/plain line / human fibroblast | 2 | 0 | 1 | 1 | 1 | 2 | 1 | 0 | 1 | 0 | 1 | 1 | 1 | 1 | 1 | 1 | 2 | 17 |
| C1a genotype WT, treatment NT, title says WT/plain line / human iPSC | 0 | 0 | 3 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 4 |
| C1b genotype WT, treatment NT, title says NTg transgene / human fibroblast | 0 | 0 | 0 | 0 | 0 | 4 | 0 | 4 | 0 | 2 | 0 | 0 | 0 | 4 | 0 | 0 | 0 | 14 |
| C2 genotype WT, rapamycin / human fibroblast | 1 | 0 | 3 | 0 | 0 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 3 | 0 | 0 | 3 | 1 | 14 |
| C3 genotype WT, dox (title says NTg transgene) / human fibroblast | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 0 | 8 |
| C4 genotype NTg, dox / human fibroblast | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 0 | 0 | 2 | 0 | 0 | 2 | 0 | 0 | 0 | 8 |
| C5 reporter (GFP), dox / human iPSC | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 |
| C5 reporter (mTagBFP2), dox / human fibroblast | 0 | 0 | 0 | 0 | 0 | 26 | 0 | 5 | 0 | 2 | 8 | 2 | 0 | 30 | 0 | 0 | 0 | 73 |
| C6 other transgene, treatment NT (not induced) / human fibroblast | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 |
| C6 other transgene, treatment NT (not induced) / human iPSC | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 |

Labels are this file's shorthand for field combinations, not the authors'.
Primed labels (C3', C4', C5') would be rows with zero GSMs if absent above.
- C1a genotype WT, treatment NT, title says WT/plain line: treatments `NT` 21; batches LIB050518 8, LIB048130 5, LIB047784 4, LIB049469 4
- C1b genotype WT, treatment NT, title says NTg transgene: treatments `NT` 14; batches LIB054839 8, LIB053907 6
- C2 genotype WT, rapamycin: treatments `Rapa_100ng/mL_24h` 4, `Rapa_100ng/mL_48h` 4, `Rapa_100ng/mL_72h` 4, `Rapa_100ng/mL, 3days` 2; batches LIB049469 12, LIB048130 2
- C3 genotype WT, dox (title says NTg transgene): treatments `Dox_1000ng/mL, 3days` 6, `Dox_2000ng/mL, 3days` 2; batches LIB052667 6, LIB054839 2
- C4 genotype NTg, dox: treatments `Dox_1000ng/mL, 3days` 6, `Dox_2000ng/mL, 3days` 2; batches LIB054839 8
- C5 reporter (GFP), dox: treatments `Dox_250ng/mL, 3days` 1, `Dox_2000ng/mL, 3days` 1; batches LIB047784 2
- C5 reporter (mTagBFP2), dox: treatments `Dox_2000ng/mL, 3days` 12, `Dox_1000ng/mL, 3d-ON` 12, `Dox_1000ng/mL, 3days` 11, `Dox_1000ng/mL, 5d-ON/2d-OFF` 6, `Dox_1000ng/mL, 3d-ON/16d-OFF` 4, `Dox_1000ng/mL, 3d-ON/8d-OFF` 4, `Dox_1000ng/mL, 1d-ON` 4, `Dox_1000ng/mL, 1d-ON/8d-OFF` 4, `Dox_1000ng/mL, 2d-ON/8d-OFF` 4, `Dox_1000ng/mL, 2d-ON/16d-OFF` 4, `Dox_1000ng/mL, 1d-ON/16d-OFF` 4, `Dox_1000ng/mL, 2d-ON` 4; batches LIB058979 36, LIB092645 14, LIB052667 6, LIB053907 6, LIB054839 6, LIB052643 4, LIB053449 1
- C6 other transgene, treatment NT (not induced): treatments `NT` 4; batches LIB047784 2, LIB051945 2

Human GSMs with `genotype: WT` and a dox treatment: 8; of these, titles say
"Dox-inducible NTg transgene": 8 (ambiguity A3).
Of the 39 human `treatment: NT` GSMs, 14 have `Dox_0ng/mL` in the title (all C1b);
the rest have no dox token in the title. That is the only place the record
hints at what `NT` means.

`molecule` by candidate: C1a genotype WT, treatment NT, title says WT/plain line: total RNA 21; C1b genotype WT, treatment NT, title says NTg transgene: polyA RNA 14; C2 genotype WT, rapamycin: total RNA 14; C3 genotype WT, dox (title says NTg transgene): polyA RNA 8; C4 genotype NTg, dox: polyA RNA 8; C5 reporter (GFP), dox: total RNA 2; C5 reporter (mTagBFP2), dox: polyA RNA 73; C6 other transgene, treatment NT (not induced): total RNA 2, polyA RNA 2; transgene GSMs: polyA RNA 664.

### 4b. Every human GSM from a donor under 40

Human GSMs with described age < 40: 15 (fibroblast 7, iPSC 8).

| gsm | library_name | title | cell_line | cell_type | age_description | passage_description | genotype | treatment | batch | control_candidate |
|---|---|---|---|---|---|---|---|---|---|---|
| GSM9611422 | S1 | M29 WT line | AG04054B | Fibroblast | 29 | 3 | WT | NT | LIB047784 | C1a genotype WT, treatment NT, title says WT/plain line |
| GSM9611434 | S13 | M29 WT line, No Rapa | AG04054B | Fibroblast | 29 | 4 | WT | NT | LIB048130 | C1a genotype WT, treatment NT, title says WT/plain line |
| GSM9611435 | S14 | M29 WT line, Rapa_100, 3 days | AG04054B | Fibroblast | 29 | 4 | WT | Rapa_100ng/mL, 3days | LIB048130 | C2 genotype WT, rapamycin |
| GSM9611426 | S5 | F3 line with Dox-inducible GFP transgene, Dox_250ng/mL | F3 | iPSC | 31 |  | GFP | Dox_250ng/mL, 3days | LIB047784 | C5 reporter (GFP), dox |
| GSM9611427 | S6 | F3 line with Dox-inducible GFP transgene, Dox_2000ng/mL | F3 | iPSC | 31 |  | GFP | Dox_2000ng/mL, 3days | LIB047784 | C5 reporter (GFP), dox |
| GSM9611428 | S7 | F3 line with Dox-inducible OTX2 transgene | F3 | iPSC | 31 |  | OTX2 | NT | LIB047784 | C6 other transgene, treatment NT (not induced) |
| GSM9611429 | S8 | F3 line with Dox-inducible ZGLP1 transgene | F3 | iPSC | 31 |  | ZGLP1 | NT | LIB047784 | C6 other transgene, treatment NT (not induced) |
| GSM9611433 | S12 | F3 line with Dox-inducible ZGLP1 transgene, Dox_2000ng/mL | F3 | iPSC | 31 |  | ZGLP1 | Dox_2000ng/mL, 3days | LIB048130 | - |
| GSM9611438 | S17 | M36 WT line, No Rapa | AG11557A | Fibroblast | 36 | 5 | WT | NT | LIB049469 | C1a genotype WT, treatment NT, title says WT/plain line |
| GSM9611439 | S18 | M36 WT line, Rapa_100ng/mL, 24h | AG11557A | Fibroblast | 36 | 5 | WT | Rapa_100ng/mL_24h | LIB049469 | C2 genotype WT, rapamycin |
| GSM9611440 | S19 | M36 WT line, Rapa_100ng/mL, 48h | AG11557A | Fibroblast | 36 | 5 | WT | Rapa_100ng/mL_48h | LIB049469 | C2 genotype WT, rapamycin |
| GSM9611441 | S20 | M36 WT line, Rapa_100ng/mL, 72h | AG11557A | Fibroblast | 36 | 5 | WT | Rapa_100ng/mL_72h | LIB049469 | C2 genotype WT, rapamycin |
| GSM9611430 | S9 | B6 line, P0 | B6 | iPSC | 36 |  | WT | NT | LIB048130 | C1a genotype WT, treatment NT, title says WT/plain line |
| GSM9611431 | S10 | B6 line, P+2 | B6 | iPSC | 36 |  | WT | NT | LIB048130 | C1a genotype WT, treatment NT, title says WT/plain line |
| GSM9611432 | S11 | B6 line, P+10 | B6 | iPSC | 36 |  | WT | NT | LIB048130 | C1a genotype WT, treatment NT, title says WT/plain line |

**Untreated (`genotype WT`, `treatment NT`) human fibroblast GSMs from donors under 40: 3.**
Their batches: LIB047784 1, LIB048130 1, LIB049469 1.

### 4c. Batch co-measurement of the under-40 fibroblast GSMs

For each batch holding an under-40 human fibroblast GSM, all human GSMs in
that batch by group and age (descriptive; the record does not say batch is a
sequencing run, a plate or a culture round):

| batch / group | age 29 | age 31 | age 36 | age 55 | age 69 | age 73 | age 86 | age 92 | age 94 |
|---|---|---|---|---|---|---|---|---|---|
| LIB047784 / human fibroblast | 1 | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 1 |
| LIB047784 / human iPSC | 0 | 4 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| LIB048130 / human fibroblast | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 |
| LIB048130 / human iPSC | 0 | 1 | 3 | 0 | 0 | 0 | 0 | 0 | 0 |
| LIB049469 / human fibroblast | 0 | 0 | 4 | 4 | 0 | 4 | 0 | 4 | 0 |

All batches ("transgene" = genotype other than `WT`/`NTg`/reporter):

| ch_batch | human fibroblast | human iPSC | mouse | described ages | genotypes | molecule | platform | fibroblast transgene GSMs | fibroblast WT+NT GSMs | fibroblast WT+NT ages |
|---|---|---|---|---|---|---|---|---|---|---|
| LIB047784 | 4 | 4 | 0 | 29/31/69/86/94 | GFP, OTX2, WT, ZGLP1 | total RNA | GPL18573 | 0 | 4 | 29/69/86/94 |
| LIB048130 | 4 | 4 | 0 | 29/31/36/94 | WT, ZGLP1 | total RNA | GPL18573 | 0 | 2 | 29/94 |
| LIB049469 | 16 | 0 | 0 | 36/55/73/92 | WT | total RNA | GPL18573 | 0 | 4 | 36/55/73/92 |
| LIB050518 | 7 | 1 | 0 | 41/52/55/60/66/68/79 | WT | total RNA | GPL18573 | 0 | 7 | 41/52/55/60/66/68/79 |
| LIB051945 | 4 | 0 | 0 | 65 | SUMO1 | polyA RNA | GPL18573 | 4 | 0 | - |
| LIB052643 | 4 | 0 | 0 | 55/79 | mTagBFP2 | polyA RNA | GPL18573 | 0 | 0 | - |
| LIB052667 | 480 | 0 | 0 | 55/65/79 | 98 values | polyA RNA | GPL24676 | 468 | 0 | - |
| LIB053449 | 7 | 0 | 0 | 65 | CYCS, MDM2, UBA52, UBE2I, mTagBFP2 | polyA RNA | GPL18573 | 6 | 0 | - |
| LIB053907 | 60 | 0 | 0 | 55/65/79 | 10 values | polyA RNA | GPL24676 | 48 | 6 | 55/65/79 |
| LIB054839 | 80 | 0 | 0 | 55/65/67/68/69/79 | 13 values | polyA RNA | GPL24676 | 56 | 8 | 55/65/67/79 |
| LIB058979 | 96 | 0 | 0 | 55/79 | OSKM, SRSF1, mTagBFP2 | polyA RNA | GPL24676 | 60 | 0 | - |
| LIB083102 | 0 | 0 | 32 | - | GFP, SRSF1 | polyA RNA | GPL34290 | 0 | 0 | - |
| LIB092645 | 36 | 0 | 0 | 55/68/79 | OSKM, SRSF1, mTagBFP2 | polyA RNA | GPL34284 | 22 | 0 | - |

Batches holding both an under-40 human fibroblast GSM and a fibroblast
transgene GSM: **none**.

## 5. Timepoints, doses, induction conditions

Source: `!Sample_characteristics_ch1` `treatment:`. Treatment protocol text
(identical in all GSMs: True):

> Transgene was induced using doxycycline at the specified concentration for the specified amount of time. Similarly, rapamycin was added to the media at the 100nM for the specified amount of time.

The protocol text says rapamycin at 100 nM; the field values say
`Rapa_100ng/mL`. Recorded, not resolved (ambiguity A6).

| treatment | human fibroblast | human iPSC | mouse | total | n genotypes | genotypes (if <= 6) | fibroblast ages |
|---|---|---|---|---|---|---|---|
| Dox_1000ng/mL, 3days | 501 | 0 | 0 | 501 | 99 |  | 55/65/79 |
| Dox_2000ng/mL, 3days | 118 | 2 | 0 | 120 | 14 |  | 55/65/67/68/69/79 |
| NT | 33 | 6 | 0 | 39 | 4 | OTX2, SUMO1, WT, ZGLP1 | 29/36/41/52/55/60/65/66/67/68/69/73/79/86/92/94 |
| Dox_1000ng/mL, 3d-ON | 36 | 0 | 0 | 36 | 3 | OSKM, SRSF1, mTagBFP2 | 55/68/79 |
| Dox_200mg/kg, 1wk-ON/1wk-OFF | 0 | 0 | 32 | 32 | 2 | GFP, SRSF1 | - |
| Dox_1000ng/mL, 1d-ON | 12 | 0 | 0 | 12 | 3 | OSKM, SRSF1, mTagBFP2 | 55/79 |
| Dox_1000ng/mL, 1d-ON/16d-OFF | 12 | 0 | 0 | 12 | 3 | OSKM, SRSF1, mTagBFP2 | 55/79 |
| Dox_1000ng/mL, 1d-ON/8d-OFF | 12 | 0 | 0 | 12 | 3 | OSKM, SRSF1, mTagBFP2 | 55/79 |
| Dox_1000ng/mL, 5d-ON/2d-OFF | 12 | 0 | 0 | 12 | 2 | SRSF1, mTagBFP2 | 55/68/79 |
| Dox_1000ng/mL, 2d-ON/8d-OFF | 12 | 0 | 0 | 12 | 3 | OSKM, SRSF1, mTagBFP2 | 55/79 |
| Dox_1000ng/mL, 2d-ON | 12 | 0 | 0 | 12 | 3 | OSKM, SRSF1, mTagBFP2 | 55/79 |
| Dox_1000ng/mL, 3d-ON/16d-OFF | 8 | 0 | 0 | 8 | 2 | SRSF1, mTagBFP2 | 55/79 |
| Dox_1000ng/mL, 3d-ON/8d-OFF | 8 | 0 | 0 | 8 | 2 | SRSF1, mTagBFP2 | 55/79 |
| Dox_1000ng/mL, 2d-ON/16d-OFF | 8 | 0 | 0 | 8 | 2 | SRSF1, mTagBFP2 | 55/79 |
| Rapa_100ng/mL_24h | 4 | 0 | 0 | 4 | 1 | WT | 36/55/73/92 |
| Rapa_100ng/mL_72h | 4 | 0 | 0 | 4 | 1 | WT | 36/55/73/92 |
| Rapa_100ng/mL_48h | 4 | 0 | 0 | 4 | 1 | WT | 36/55/73/92 |
| Rapa_100ng/mL, 3days | 2 | 0 | 0 | 2 | 1 | WT | 29/94 |
| Dox_250ng/mL, 3days | 0 | 1 | 0 | 1 | 1 | GFP | - |

Reading of the codes, as written and not interpreted further: `Dox_<c>` is
doxycycline concentration; `<n>d-ON/<m>d-OFF` reads as n days on dox then m
days off; `3days` and `3d-ON` are both written and may or may not be the same
schedule (ambiguity A5); `NT` is not defined in the record.

Title tokens: `P<n>` / `P+<n>` (766 human GSMs) and `Rep<n>` (767 human GSMs).
Where both a title `P<n>` and a described `Passage N` exist (763 GSMs), they are equal in 761.
- unequal: GSM9611955 title `M65 line with Dox-inducible MDM2 transgene, Dox_1000ng/mL, P10, Rep1, redo`; description `Donor age 65, Passage 9`
- unequal: GSM9611956 title `M65 line with Dox-inducible UBA52 transgene, Dox_1000ng/mL, P10, Rep2, redo`; description `Donor age 65, Passage 9`
Described passage present for 799 of 807 human GSMs.

## 6. Human versus mouse

| organism | n GSMs | platforms | sample kind |
|---|---|---|---|
| Homo sapiens | 807 | GPL18573, GPL24676, GPL34284 | `cell line` + `cell type` (Fibroblast 798, iPSC 9) |
| Mus musculus | 32 | GPL34290 | `tissue` (Brain 8, Heart 8, Liver 8, Spleen 8) |

Human and mouse are separated by `!Sample_organism_ch1`; no GSM lacks it.
The human counts file is a separate supplementary file from the mouse one;
its header is compared with the human GSMs in section 7.

Mouse GSMs (all AAV8, dox 200 mg/kg; `treatment` `Dox_200mg/kg, 1wk-ON/1wk-OFF`):

| genotype / mouse | tissues | age_at_death_days | RIN | batch |
|---|---|---|---|---|
| GFP / Mouse#5 | Brain/Heart/Liver/Spleen | 785 | 7.0/7.3/8.1/7.2 | LIB083102 |
| GFP / Mouse#7 | Brain/Heart/Liver/Spleen | 778 | 6.8/5.7/5.8/5.8 | LIB083102 |
| GFP / Mouse#8 | Brain/Heart/Liver/Spleen | 865 | 6.7/6.2/8.1/4.1 | LIB083102 |
| GFP / Mouse#9 | Brain/Heart/Liver/Spleen | 806 | 6.8/5.5/7.4/3.7 | LIB083102 |
| SRSF1 / Mouse#1 | Brain/Heart/Liver/Spleen | 901 | 7.3/6.0/4.4/5.8 | LIB083102 |
| SRSF1 / Mouse#10 | Brain/Heart/Liver/Spleen | 872 | 7.4/5.5/8.5/8.4 | LIB083102 |
| SRSF1 / Mouse#15 | Brain/Heart/Liver/Spleen | 872 | 7.6/6.9/7.7/8.0 | LIB083102 |
| SRSF1 / Mouse#19 | Brain/Heart/Liver/Spleen | 818 | 7.2/6.1/7.8/7.7 | LIB083102 |

Mouse: 8 animals (GFP 4, SRSF1 4), 4 tissues each. The mouse ages are "Age at death",
in days, from `!Sample_description`. No young mouse group is recorded.

## 7. Gene identifiers

From the 262,144-byte range read of `GSE325735_human_raw_counts.txt.gz`:

- first header field: `Geneid`; sample columns: 807; line endings CRLF: True
- complete data rows decoded: 494 (the file's total row count is not known
  without reading it all)
- identifier pattern of those rows: `ENSG###########.v` 494
- first 5 identifiers: `ENSG00000223972.5`, `ENSG00000227232.5`, `ENSG00000278267.1`, `ENSG00000243485.5`, `ENSG00000284332.1`
- identifier type: versioned Ensembl gene IDs (`ENSG...`.version); no
  symbol column in the rows read
- header columns vs human GSM library names (`Library name:` in
  `!Sample_description`): equal sets: True; in header only: 0; in GEO only: 0

**Frozen-ruler overlap: not counted.** The ruler (`results/fibro/frozen_ruler_ridge_raw.npz`,
n = 23485 genes, keys `gene_id` versioned and `ensembl` unversioned) can only be
matched against the full identifier column, which needs the whole
26,007,361-byte matrix. That is a matrix download, so it was not done.

Key-format check only, not an overlap estimate. Among the 494 identifiers read,
223 match a ruler `ensembl` once the version suffix is dropped; 62 of those also
match the ruler's versioned `gene_id` exactly. So the Ensembl versions in this
file and in the ruler's annotation differ for some genes, and any later match
would have to be made on unversioned IDs (a choice for that task, not made
here). These rows are the first in file order, not a sample; the numbers say
nothing about how many of the 23,485 ruler genes are present.

The mouse counts file (1,773,429 B) was not opened.

## 8. Ambiguities, recorded and not resolved

- **A1.** Series title. GEO now reads "Transcriptomic reprogramming screen identifies SRSF1 as a regulator of aging"; the task names a different title. The record was updated Oct 14 2026.
- **A2.** Donor identity and age. Age is free text in `!Sample_description`, not a structured characteristic. The title `M<n>` equals it wherever both exist, but the record never says what `M` means. Sex is not stated. Whether each `cell line` value is a distinct donor is not stated; ages 36 and 55 each belong to two lines. Coriell catalog entries for the AG IDs were not consulted.
- **A3.** `NTg`. 8 human GSMs have `genotype: NTg`; 22 more have `genotype: WT` but a title reading "line with Dox-inducible NTg transgene" (8 with dox, 14 with `NT`/`Dox_0ng/mL`). The record does not define NTg, nor say whether these lines carry a construct. So `genotype: WT` covers both plain "WT line" GSMs and NTg-titled GSMs.
- **A4.** What is "untreated". The record does not use the word, and `NT` is undefined (NTg-titled `NT` rows say `Dox_0ng/mL`). Candidates: C1a plain WT line, no dox; C1b NTg-titled, 0 dox; C3/C4 NTg with dox; C5 `mTagBFP2` reporter with dox. Below age 55 only C1a (and rapamycin-treated C2) exist.
- **A5.** Schedules. `Dox_1000ng/mL, 3days` and `Dox_1000ng/mL, 3d-ON` are both used; whether they are the same exposure, and what precedes harvest after a `d-OFF` period, is not stated beyond the code.
- **A6.** Rapamycin dose: protocol text says 100 nM, field values say 100 ng/mL.
- **A7.** Label disagreements between title and `genotype`: 2 GSMs titled `OSKM` carry `genotype: SOX2`; 2 GSMs titled `TP53 (DN)` carry `genotype: TP53`. (`TP53 (DN)` is the genotype of 4 other GSMs.)
- **A8.** Cell type. The design mentions "iPSCs and their derivatives"; `cell type` holds only `Fibroblast` and `iPSC`, so no derivative is labelled. The iPSC GSMs carry a `Donor age`, presumably of the source donor; the record does not say so. AG07141 contributes one Fibroblast and one iPSC GSM.
- **A9.** Batch. The `batch` field holds `LIB...` identifiers; the description adds `Batch 1`-`Batch 12` for 96 GSMs. Neither is defined as culture, plate or sequencing batch.
- **A10.** Depth. SRA spots are pre-downsampling. The processing text says libraries were downsampled to 20 million reads before counting (three LIB IDs exempt) and does not say reads or read pairs. GEO counts were not opened, so post-downsampling depth is unknown here.
- **A11.** Passage. Two GSMs have a title `P<n>` that disagrees with the described passage (listed in section 5). B6 titles use `P0`, `P+2`, `P+10`.
- **A12.** Molecule. 40 GSMs are `total RNA`, 799 `polyA RNA`; the extract protocol describes a polyA workflow for NEBNext libraries and never mentions total-RNA libraries. Every under-40 fibroblast GSM is `total RNA`.
- **A13.** Gene overlap. Not counted (section 7). Ensembl versions differ between the file and the ruler annotation for some genes.

## 9. Plain statement (not a result)

What a distance-to-young readout would need from this series:

**(a) Perturbations in donor-aged human fibroblasts: present.** 664 human
fibroblast GSMs carry 96 distinct genotype values other than `WT`/`NTg`/reporter,
from 6 lines at stated donor ages 55, 65, 67, 68, 69, 79. 89 of those values occur at more than one
donor age. The ages are stated in `!Sample_description`, not inferred.

**(b) Untreated young human fibroblasts of the same cell type, same series:
present but thin, and ambiguous.** Fibroblast GSMs with `genotype: WT`,
`treatment: NT` from donors under 40: 3 (GSM9611422 age 29, line AG04054B, passage 3, batch LIB047784, GPL18573 total RNA, 14,208,426 SRA spots; GSM9611434 age 29, line AG04054B, passage 4, batch LIB048130, GPL18573 total RNA, 15,000,551 SRA spots; GSM9611438 age 36, line AG11557A, passage 5, batch LIB049469, GPL18573 total RNA, 11,327,576 SRA spots).
They are labelled `cell type: Fibroblast`, as the perturbed GSMs are. What
makes them ambiguous as a reference, stated without resolving it:

- no under-40 fibroblast GSM shares a `batch` value with any transgene GSM (section 4c);
- none is a reporter (`mTagBFP2`) or dox-exposed control: those exist only at ages 55, 65, 67, 68, 69, 79;
- their described passages are 3, 4, 5; transgene GSMs' passages span 8-15;
- all 3 are `total RNA` on GPL18573; the transgene GSMs are GPL18573 polyA RNA 10, GPL24676 polyA RNA 632, GPL34284 polyA RNA 22;
- untreated fibroblasts aged 40 to 54: age 41 1 GSM(s), age 52 1 GSM(s); the youngest transgene line is 55.
- the under-40 fibroblast batches (LIB047784, LIB048130, LIB049469) hold only plain "WT line" GSMs, iPSC GSMs
  and rapamycin GSMs (section 4c).
- read together: the plain "WT line" GSMs (C1a, C2; 31 fibroblast GSMs, 15 lines,
  ages 29-94) are all `total RNA`, in batches LIB047784, LIB048130, LIB049469, LIB050518. The transgene GSMs come from
  ages 55, 65, 67, 68, 69, 79; the untreated rows inside the screen's batches (C1b, polyA) come
  only from ages 55, 65, 67, 79. The record
  does not say whether the WT-line set and the screen are one experiment.

The 8 iPSC GSMs from donors under 40 are a different cell type and are not
counted under (b).

## 10. Disclosures

- No flag file: the task is an inventory with no pre-registered statistic.
  No count here is a gated result, so the chance-baseline, cell-count and
  recovery-power obligations of `CLAUDE.md` have nothing to attach to.
- Raw metadata stored under `data/raw/srsfinv/` (git-ignored); the script
  refetches a missing file. md5s are in section 0 and `inventory_meta.json`.
- Beyond the GEO record: the SRA run table (for depth; GEO has none) and a
  262,144-byte HTTP range read of the human counts file (for identifier type).
  Neither is an expression matrix download.
- The frozen ruler npz was opened for its `gene_id`/`ensembl` arrays only; no
  weight, mean or sd was read or used.
- Added beyond the brief: section 0b (scoping-parse comparison), section 3c
  (title vs genotype), section 4c (batch co-measurement), the C1-C6 shorthand.
- Decided here because the brief did not specify: "donor" is counted as a
  distinct `cell line` value; "transgene" excludes `WT`, `NTg`, `mTagBFP2`, `GFP`
  only where a sentence says so, and the full genotype table lists all values;
  "under 40" means described age < 40.
- Not covered: anything about donors younger than 55 receiving a transgene
  (there are none); the mouse arm beyond its sample table; expression values.

Machine-readable tables in `results/srsfinv/`: `samples.csv`,
`samples_human_annotated.csv`, `donors.csv`, `genotype_counts.csv`,
`genotype_by_age_fibroblast.csv`, `control_candidates_by_age.csv`,
`under40_samples.csv`, `batches.csv`, `treatments.csv`, `sra_depth_per_gsm.csv`,
`scoping_parse_comparison.csv`, `inventory_meta.json`.

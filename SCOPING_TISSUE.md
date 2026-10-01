# SCOPING_TISSUE — is there a tissue where steps 2 and 3 are actually runnable?

**Status:** search complete 2026-09-17. Fit nothing. Download nothing large. Does not modify any `FINDINGS_*.md`.

**Verdict: STRONG GO, skin / dermal fibroblast.** (A) GTEx V10 `Cells - Cultured fibroblasts` n=652, AGE bins 20–29 to 70–79, SMRIN miss=0. (B) GSE297234 Sendai OSKM days 0/3/7/10, donors 22 yr and 96 yr. (C) GSE325735 NHDF cDNA overexpression, 94 distinct transgenes, sample_description donor ages include 79 and 94.

## Pre-registered go/no-go (verbatim, written 2026-09-17 before any accession search)

Written **before any GEO/SRA/CXG/GTEx sample-table fetch**, 2026-09-17. No bar in this block is re-tuned after numbers exist.

Two rulers from FINDINGS_GEOMETRY.md (unchanged):
- **age ruler** — within-type supervised direction over genes (ridge / PLS-1).
- **identity ruler** — SVD of type centroids. In a culture with one type this axis degenerates.

FINDINGS_TRAJECTORY.md: a 3-condition protocol contrast cannot identify a bend. That is why (B) requires ≥4 distinct timepoints.

FINDINGS_BRAIN_PHASE1.md §2 / FINDINGS.md: blood age signal is compositional. Blood is scored on (B) and (C) only. (A) is not re-litigated.

Per tissue:
- **(A) passes** if a public cohort has ≥150 donors, adult age span ≥40 years, and per-sample RIN or an equivalent quality field. Bulk counts; single-cell is better and is noted.
- **(B) passes** if a public dataset has ≥4 distinct timepoints of partial reprogramming or another rejuvenation perturbation, in primary cells from donors with stated age, **at least one donor ≥60**, with transcriptomic readout. ≥4 timepoints is the bar because a curve with 3 points cannot show a bend (FINDINGS_TRAJECTORY's original objection).
- **(C) passes** if a public dataset has ≥20 distinct perturbations applied to primary cells of that tissue from donors with stated age, at least one ≥60, with transcriptomic readout.
- **GO** = a single tissue passes (A) and (B). **STRONG GO** = the same tissue passes (A), (B) and (C).
- **NO-GO** = no tissue passes (A)+(B). Then say so plainly and list, for the closest near-miss, the single missing ingredient.

Do not relax a bar after seeing what exists. If a dataset misses a bar by a little, it goes in a `near_miss` block with the exact shortfall (e.g. "3 timepoints, needs 4"; "ages not stated in record").

Tissues checked, in this order: (1) skin / dermal fibroblast; (2) blood / PBMC — (B)+(C) only; (3) muscle, adipose, lung; (4) iPSC-derived neurons from aged donors (reprogramming-to-iPSC erases age markers; stated next to any hit). Also: Joung 2023 TF Atlas non-hESC arm.

Hard rule: every dataset in the tables was resolved from the record page in this session.

## 1. Skin / dermal fibroblast

| criterion | best dataset | accession | n donors | age span | timepoints or perturbations | access tier | passes |
|---|---|---|---|---|---|---|---|
| (A) | GTEx V10 RNASEQ `Cells - Cultured fibroblasts` | phs000424.v10 / GCS v10 annotations | 652 | public AGE bins 20–29 to 70–79 (span ≥50 y) | n/a (cohort) | annotations + gene-reads open; exact year age dbGaP controlled | True |
| (B) | Sendai OSKM scRNA time course, GM00731 + GM23815 | GSE297234 | 2 | 22 yr and 96 yr | days 0, 3, 7, 10 (4) | open (GEO + SRA PRJNA1263211) | True |
| (C) | NHDF dox-inducible cDNA overexpression screen | GSE325735 | 16 M-coded NHDF lines (GSM-verified ages include 55, 65, 68, 79, 94) | 29–94 as M-codes in titles; verified `Donor age` in `!Sample_description` | 94 distinct transgenes (93 parsed gene/OSKM names + TP53-DN) | open (GEO; pdat 2026/08/15) | True |

(A) also passes on GTEx bulk `Skin - Sun Exposed (Lower leg)` n=754 / `Skin - Not Sun Exposed (Suprapubic)` n=651, same AGE bins, SMRIN miss=0. Those are mixed-cell bulk, not the cultured fibroblast used in (B)/(C). Single-cell skin atlases on CELLxGENE were listed (96 title/tissue hits); per-dataset donor counts were not recovered (dataset GET 404) and are not used as (A).

## 2. Blood / PBMC

(A) is closed (compositional; FINDINGS_BRAIN_PHASE1 §2). Not re-scored.

| criterion | best dataset | accession | n donors | age span | timepoints or perturbations | access tier | passes |
|---|---|---|---|---|---|---|---|
| (B) | none resolved | — | — | — | — | — | False |
| (C) | none resolved | — | — | — | — | — | False |

Closest resolved blood-adjacent hit: GSE254389 NEUROD1 conversion of T cells to neurons, n=18, two induction methods, ages not in series design. Not a ≥4-timepoint rejuvenation course and not ≥20 perturbations.

## 3. Muscle, adipose, lung

| tissue | criterion | best dataset | accession | n donors | age span | timepoints or perturbations | access tier | passes |
|---|---|---|---|---|---|---|---|---|
| muscle | (A) | GTEx V10 RNASEQ `Muscle - Skeletal` | phs000424.v10 | 818 | bins 20–29 to 70–79 | n/a | annotations + gene-reads open; exact year dbGaP | True |
| muscle | (B) | none resolved in human primary skeletal muscle | — | — | — | — | — | False |
| muscle | (C) | none resolved | — | — | — | — | — | False |
| adipose | (A) | GTEx V10 RNASEQ `Adipose - Subcutaneous` | phs000424.v10 | 714 | bins 20–29 to 70–79 | n/a | same | True |
| adipose | (A) alt | GTEx V10 RNASEQ `Adipose - Visceral (Omentum)` | phs000424.v10 | 587 | bins 20–29 to 70–79 | n/a | same | True |
| adipose | (B) | none resolved | — | — | — | — | — | False |
| adipose | (C) | none resolved | — | — | — | — | — | False |
| lung | (A) | GTEx V10 RNASEQ `Lung` | phs000424.v10 | 604 | bins 20–29 to 70–79 | n/a | same | True |
| lung | (B) | none resolved in aged primary lung | — | — | — | — | — | False |
| lung | (C) | none resolved in aged primary lung | — | — | — | — | — | False |

Muscle (B) search returned GSE201710 (LAKI mouse + IMR90) and GSE190665 (4F mouse methylation). IMR90 is a fetal lung fibroblast line, not aged primary skeletal muscle. Adipose (B) search returned no human adipose OSK time course among resolved GSE records. Lung (B) hits are IMR90 reprogramming/senescence screens (e.g. GSE103938, GSE95021 in the hit list), not donors ≥60.

## 4. iPSC-derived neurons from aged donors

Caveat (stated before search, restated next to hits): reprogramming to iPSC erases age markers, so a derived neuron may not be "aged" at all.

| criterion | best dataset | accession | n donors | age span | timepoints or perturbations | access tier | passes |
|---|---|---|---|---|---|---|---|
| (A) | none ≥150 donors | — | — | — | — | — | False |
| (B) | none | — | — | — | — | — | False |
| (C) | GSE246052 CRISPRi iPSC-neurons (search hit) | GSE246052 | not used; n_samples=4 in esummary | not resolved as aged-donor ≥60 | 4 samples | GEO (esummary only) | False |

GSE255982 SAD iPSC RNA-seq n=52 is below the (A) n=150 bar and is iPSC, not aged primary neurons. GSE117720 blood-derived iNSC "epigenetic rejuvenation" n=19, no numeric ≥60 in the titles returned.

## Joung 2023 TF Atlas (non-hESC arm?)

No aged-primary arm. SuperSeries GSE216481, 15 SubSeries. Resolved designs:

- GSE217460: "hESCs were transduced with a library of 3,550 transcription factor ORFs"
- GSE216595: "hESCs were transduced with a library of 90 transcription factor ORFs"
- GSE216479: "introduced into pluripotent stem cells"
- GSE219000 / GSE216601: hESC → cardiomyocyte
- GSE219058 / GSE216463 / GSE216602: hESC → neural progenitor / astrocyte screen

The atlas supplies a differentiation axis from hESC, not aged cells of any tissue in this survey.

## near_miss

| accession | tissue | aimed | shortfall |
|---|---|---|---|
| GSE165177 / GSE165180 | skin fibroblast | (B) | ≥4 timepoints (10/13/15/17 d) and RNA-seq; `donor age (years)` on GSM is 38 and 53. Needs ≥60. |
| GSE297984 | skin fibroblast | (B) | RNA-seq; titles `22y/23y/56y/83y`; 83≥60; timepoints D6 and D14 only. Needs 4. |
| GSE297233 | skin fibroblast | (B) | aged donor 96 yr in overall_design; OSK vs O4YRSK at 4 days only. Needs 4 timepoints. |
| GSE142439 | skin fibroblast + vein EC | (B) | methylation not transcriptomic; 2 conditions (Normal/Treated); GSM `age: 61`. RNA is SRA PRJNA598923 (18 runs: Fib/EC × young/aged/treated). Needs 4 timepoints and a transcriptomic GEO sample table with years. |
| GSE237269 | skin fibroblast | (B) | design states GM09503 10 yr and GM08401 75 yr; 3 conditions (Young / Old / day-8 rejuvenated). Needs 4 timepoints. |
| GSE226189 | skin fibroblast | (A) | 82 donors, titles `SKIN_AGE22`–`AGE89`, `age (years): 72` on GSM. Needs n≥150. |
| GSE113957 | skin fibroblast | (A) | 143 samples; ages in titles include 1 yr–92 yr. Needs n≥150; no RIN field in the series text fetched. |
| GSE307377 | skin fibroblast | (A) | 9 donors, ages 23–72 in titles and `age:` characteristics. Needs n≥150. |
| GSE201710 | muscle-adjacent | (B) | human arm is IMR90 (fetal lung fibroblast), not primary skeletal muscle; mouse in vivo OSK. |
| GSE246052 | iPSC-neuron | (C) | n_samples=4. Needs ≥20 perturbations and a stated donor ≥60. |

## unverified

Accessions encountered but not given a usable sample table in this session:

- PRJNA1108254 — esearch n=1090 SRA runs (Southard/Sengstack fibroblast Perturb-seq, from paper/GitHub); esummary HTTP 429. Donor age not read from a record.
- Figshare 10.6084/m9.figshare.30898748 — cited as the h5ad for that Perturb-seq; page not fetched.
- CELLxGENE `a19d1667-a7b5-4556-9e5f-f9bfa690c0f1` "Human dermal fibroblast atlas" (108440 cells in the list payload) — dataset GET 404; n_donors not read.
- LINCS L1000 fibroblast GEO (e.g. GSE92742) — not fetched; the query `LINCS AND fibroblast AND Homo sapiens` returned only GSE169045 (iPSC-EC, n=36).

## Winning tissue

**Skin / dermal fibroblast.** Age-ruler training cohort: GTEx Analysis V10 RNASEQ freeze (`SMAFRZE==RNASEQ`) `SMTSD==Cells - Cultured fibroblasts`. n=652 donors, one sample per donor. Public `AGE` bins 20–29 (n=60), 30–39 (53), 40–49 (106), 50–59 (202), 60–69 (210), 70–79 (21). SMRIN miss=0, median 9.9 (6.7–10.0). `SMCENTER` donor counts: B1=421, C1=224, D1=7. Two-site transfer (B1 vs C1) is testable; a three-site test is not (D1=7). Bulk `Skin - Sun Exposed` n=754 is an alternate (A) matrix (SMCENTER B1=509, C1=241, D1=4) but is mixed-cell bulk, not the cultured fibroblast used in GSE297234 / GSE325735.

Identity(type) ruler: not meaningful in a fibroblast culture. GTEx cultured fibroblasts and the (B)/(C) plates are one cell type. The 20-type DLPFC identity subspace has no analogue. Safety against identity loss would have to use a differentiation / fibroblast-state axis (the Joung atlas is hESC, not aged NHDF). That is the hidden risk of this pivot.

Gene-universe overlap with the DLPFC 25526-gene matrix was not computed for skin in this survey. GTEx V10 `gene_reads` raw n_genes=59033 (FINDINGS_GTEX Stage 0, same GCT). After the frozen GTEx filter, muscle overlap was 17943 and heart 18226; skin/fibroblast overlap is unknown until that filter is run.

## Limitations

1. GTEx public AGE is a 10-year bin. Exact years need dbGaP `phs000424.v10`.
2. GSE297234 has two donors. The time course can show a bend; it cannot estimate donor-to-donor variance.
3. GSE325735 mixes human NHDF with mouse AAV samples (32 mouse-like titles of 839). Use the human fibroblast rows only.
4. GSE325735 ages are in `!Sample_description` ("Donor age 79"), not in `!Sample_characteristics_ch1`. Five GSMs were opened; remaining M-codes were not all GSM-fetched.
5. Joung has no aged-cell arm; it cannot supply the missing identity ruler for NHDF.
6. Blood (B)/(C) and muscle/adipose/lung (B)/(C) were searched via NCBI GDS esearch + esummary; no resolved human primary hit passed the bars. A dataset that exists but was not returned by those queries is outside this file.

## Next

Freeze GTEx V10 cultured-fibroblast samples as the age-ruler cohort and score GSE297234 days 0/3/7/10 on that frozen direction; measure gene overlap with DLPFC before opening GSE325735.

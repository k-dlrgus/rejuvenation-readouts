# PROGRESS_SCOPING

**STOP status:** complete. Verdict **STRONG GO, skin / dermal fibroblast**.
**Next action:** freeze GTEx V10 cultured-fibroblast samples as the age-ruler cohort and score GSE297234 days 0/3/7/10 on that frozen direction; measure gene overlap with DLPFC before opening GSE325735.

## Pre-registered bars (verbatim, 2026-09-17, before any accession search)

Per tissue:
- **(A) passes** if a public cohort has ≥150 donors, adult age span ≥40 years, and per-sample RIN or an equivalent quality field. Bulk counts; single-cell is better and is noted.
- **(B) passes** if a public dataset has ≥4 distinct timepoints of partial reprogramming or another rejuvenation perturbation, in primary cells from donors with stated age, **at least one donor ≥60**, with transcriptomic readout. ≥4 timepoints is the bar because a curve with 3 points cannot show a bend (FINDINGS_TRAJECTORY's original objection).
- **(C) passes** if a public dataset has ≥20 distinct perturbations applied to primary cells of that tissue from donors with stated age, at least one ≥60, with transcriptomic readout.
- **GO** = a single tissue passes (A) and (B). **STRONG GO** = the same tissue passes (A), (B) and (C).
- **NO-GO** = no tissue passes (A)+(B). Then say so plainly and list, for the closest near-miss, the single missing ingredient.

Do not relax a bar after seeing what exists.

## Verified candidates so far

| accession | criterion | tissue | pass |
|---|---|---|---|
| GTEx V10 `Cells - Cultured fibroblasts` n=652 | A | skin fibroblast | True |
| GTEx V10 Skin sun-exposed n=754 / not-sun n=651 | A | bulk skin | True |
| GTEx V10 Muscle n=818, Adipose-Sub n=714, Adipose-Vis n=587, Lung n=604 | A | those tissues | True |
| GSE297234 days 0/3/7/10, donors 22 and 96 yr | B | skin fibroblast | True |
| GSE325735 94 transgenes, Donor age 79/94 in GSM description | C | skin fibroblast | True |
| GSE216481 / GSE217460 Joung | TF atlas | hESC only | False (no aged arm) |

Blood (A) closed, not re-scored. Blood (B)/(C) fail. Muscle/adipose/lung (B)/(C) fail. iPSC-neurons fail (A)(B)(C); iPSC-erasure caveat stands.

## Searches already run

Do not repeat:

1. GTEx V10 annotations counted (`SMAFRZE==RNASEQ`) and HTTP-confirmed on GCS.
2. NCBI GDS esearch queries in `results/scoping/ncbi_searches.json` (partial/transient reprogramming, OSKM/OSK fibroblast, CRISPRa/CRISPR screen fibroblast, tissue OSKM, iPSC-neuron, LINCS, Joung).
3. GEO series + GSM text in `geo_records.json`, `geo_records_round2.json`, `gsm_key.json`.
4. Joung subseries designs in `round3.json`.
5. SRA PRJNA598923, PRJNA1263211; PRJNA1108254 count only (429 on summary).
6. CXG list 2226 datasets; individual dataset GET 404.

## Single next action

Fit nothing in this prompt. Next prompt: train the age ruler on GTEx V10 cultured fibroblasts (two-site B1 vs C1) and project GSE297234.

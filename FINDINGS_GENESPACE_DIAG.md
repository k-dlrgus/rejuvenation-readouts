# FINDINGS_GENESPACE_DIAG — do the published age lists point the right way in GTEx?

**Label: diagnostic, not pre-registered.** Nothing here was in `results/genespace/PREREG.flag`. No manuscript or existing file was edited. Script: `src/genespace_diag.py`. Outputs: `results/genespace_diag/`.

## 1. Loading check

- `src/md3_genesets.py` reads sheet `Aging_signatures` and indexes columns by header name (`aging["Age up"]`, `aging["Age down"]`), with a StopStep if either header is missing.
- Independent read (stdlib, sheet resolved through `workbook.xml.rels` → `xl/worksheets/sheet2.xml`), headers: ['Age up', 'Age down'].
- `genesets.json` `age_up` equals the xlsx `Age up` column (stripped/uppercased, same order): **True** (n=1533). `age_down` equals `Age down`: **True** (n=2007). `age_up` equals the `Age down` column: False.

| list | first 5 in xlsx column | first 5 in genesets.json |
|---|---|---|
| Age up / age_up | SCN1A, STEAP4, ADH1B, KCNA4, LEP | SCN1A, STEAP4, ADH1B, KCNA4, LEP |
| Age down / age_down | CACNG4, CACNA1H, STRA6, SLC35F3, RRM2 | CACNG4, CACNA1H, STRA6, SLC35F3, RRM2 |

## 2–3. Per-gene Spearman ρ with GTEx age bin

Frozen Stage 0 fibroblast pack (TMM log2-CPM, prior 2; `toward_run.load_gtex_z`), 652 donors, age = ordinal 10-year bin (20-29 … 70-79; counts {'20-29': 60, '30-39': 53, '40-49': 106, '50-59': 202, '60-69': 210, '70-79': 21}). Ties use average ranks. Genes are the mapped ruler columns (same symbol mapping as GENESPACE).

| set | n genes | share ρ > 0 | share ρ < 0 | median ρ | IQR |
|---|---|---|---|---|---|
| age_up | 1312 | 35.1% | 64.9% | -0.0330 | [-0.085, +0.021] |
| age_down | 1777 | 69.2% | 30.8% | +0.0535 | [-0.018, +0.110] |
| all_ruler | 23485 | 49.0% | 51.0% | -0.0026 | [-0.067, +0.061] |
| ruler_not_in_lists | 20396 | 48.2% | 51.8% | -0.0047 | [-0.069, +0.058] |

Mann–Whitney age_up ρ vs age_down ρ: p = 3.31e-107.

## Consistency with the GENESPACE readout

- Recomputed cos(GTEx c_young − c_old, r_AGE) in AGE space: -0.3732 (r_AGE = −1 on age_up, +1 on age_down; FINDINGS_GENESPACE reports −0.37).
- Share of age_up genes higher in young (20–39, n=113) than old (60–79, n=231) frozen-z: 63.8%; median young−old z +0.0692.
- Share of age_down genes higher in young: 31.1%; median young−old z -0.1337.
- Age coding check: GTEx `age_ordinal` equals the bin order used here: True. Reference genes (ρ with age; MKI67 and MMP3 are in age_down): CDKN2A +0.081, CDKN1A -0.032, MKI67 +0.113, MMP3 +0.107, XIST -0.088, RPS4Y1 +0.080. XIST negative and RPS4Y1 positive mean the older bins are more male; sex is not adjusted here.
- Per-gene effects are small (IQR of ρ about ±0.1); the signal is in the consistent direction across thousands of genes, not in individual genes.

## Plain summary

- Age-up genes go up with age in GTEx: **no** (35.1% have ρ > 0 vs 49.0% of all ruler genes; median ρ -0.0330).
- Age-down genes go down with age in GTEx: **no** (30.8% have ρ < 0 vs 51.0% of all ruler genes; median ρ +0.0535).
- Loading error: **no sign of one**. Columns are read by header name and the JSON matches the xlsx columns exactly.

_diagnostic, not pre-registered_

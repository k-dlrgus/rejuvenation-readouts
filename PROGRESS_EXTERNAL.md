# PROGRESS_EXTERNAL

**STOP status:** none.
**Next action:** done

## Prompt (re-read after any context summary)

Goal: does the frozen DLPFC age direction predict age in an independent cohort it was never trained on? Cohort: SEA-AD single-nucleus MTG (Allen Institute), open access.

Stage 0: download donor metadata and per-nucleus cell-type labels; map SEA-AD subclass → 20 DLPFC types; report unmapped types and n donors per type; report age range, sex, and AD pathology fields. Log every column name read into `results/external/manifest.json`. Stop if not downloadable or fewer than 40 donors map.

Pre-register before any fit. Primary = no/low AD pathology. Secondary = all donors, pathology covariate. Pass = r > 0 with null ≤ 0.05 in the primary cell. Reading, only the outcome that fired: pass → third bank/region/lab; fail + pathology-covariate pass → masked by pathology; both fail → does not transfer beyond the two DLPFC banks.

## Pre-registered bars

- pass: r > 0 with permutation-null r ≤ 0.05 in the primary cell
- n_perm=200  n_boot=200  perm seed `20260914`  boot seed `20260918`
- directions: raw ridge and identity-residualized ridge (0.317→0.300 fit), trained on all 233 DLPFC donors
- aggregation: median over mapped types with n≥8 donors

## Seeds

- main / permutation: `20260914`
- donor-bootstrap: `20260918`
- SEA-AD dataset: `c2876b1b-06d8-4d96-a56b-5304f815b99a`
- DLPFC dataset: `4442d412-91cb-4261-acca-8adf5fa04c11`

## Finished cells

- `freeze.json`
- `gene_overlap.json`
- `primary_idresid.json`
- `primary_raw.json`
- `project_summary.json`
- `secondary_idresid.json`
- `secondary_raw.json`
- `stage0_summary.json`

## Paths written

- `results/external/manifest.json`
- `results/external/PREREG_20260917.flag`
- `PROGRESS_EXTERNAL.md`
- `FINDINGS_EXTERNAL.md`

PASS: frozen DLPFC direction transfers to SEA-AD MTG. primary raw r=+0.877 (null +0.004, p=0.005, CI [+0.639,+0.942]); idresid r=+0.873. Same reading both directions. 89 donors mapped; ADNC; age 29-89; 26 no/low.


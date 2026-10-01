# PROGRESS_FIBRO

**STOP status:** STAGE2_DONE
**Next action:** done. GSE325735 out of scope

## Stage gate

- Stage 1 prereg flag: `PREREG_STAGE1.flag` exists=True
- Stage 2 prereg flag: `PREREG_STAGE2.flag` exists=True
- Stage 1 pass (ridge raw both transfer directions): True
- Stage 1 reason: Stage 1 pass: ridge raw ρ>0 both ways with null ≤ 0.05 (B1→C1 ρ=+0.529 null=-0.004 p=+0.005; C1→B1 ρ=+0.561 null=+0.013 p=+0.005).
- Stage 2 opened: True
- frozen ruler: `frozen_ruler_ridge_raw.npz` exists=True
- seed `20260914`  boot seed `20260918`  n_perm=200  n_boot=200
- primary metric: Spearman ρ; uncalibrated R² is never a gate
- D1 exclusion: D1 n is too small for a transfer cell; excluded from B1↔C1 transfer and from within-site CV, stated, not silently dropped from the cohort table.

## Frozen preprocessing spec (FINDINGS_GTEX.md Stage 0, verbatim)

- `SMAFRZE==RNASEQ` and `SMTSD == 'Cells - Cultured fibroblasts'`
- Gene filter: ≥6 counts in ≥20% of samples
- TMM (edgeR defaults) → log2-CPM `prior.count=2`
- z-score per gene on training rows only inside every fold (Stage 1)
- Map to DLPFC gene universe, Ensembl version-stripped

## Frozen gene lists (Stage 2; not scored until Stage 1 passes)

- Pluripotency endogenous: POU5F1, NANOG, LIN28A, SALL4, DPPA4, ZFP42, DNMT3B
- Fibroblast identity: COL1A1, COL1A2, THY1, S100A4, FN1, VIM, POSTN
- OSKM-family genes dropped for the without-OSKM pluripotency score: POU5F1, OCT4, OCT3, SOX2, KLF4, MYC, MYCL, MYCN

## Frozen ruler path

- `<repo>\results\fibro\frozen_ruler_ridge_raw.npz` exists=True

## Stage 1 pre-registration (verbatim)

Written **before any Stage 1 fit**, 2026-09-17. No threshold in this block is re-tuned after numbers exist.

Cohort: GTEx Analysis V10 RNASEQ freeze (`SMAFRZE==RNASEQ`), `SMTSD == "Cells - Cultured fibroblasts"`.
One sample per donor. Target = public `AGE` bin midpoint. Model = SVD-LOO ridge (same fitter as
POSCTRL gene_ridge / TARGET age companion). PLS-1 is fit and reported; it is **not** averaged with
ridge (GTEx T3: PLS-1 has an inflated null).

Preprocessing (FINDINGS_GTEX.md Stage 0, verbatim): gene filter ≥6 counts in ≥20% of samples;
TMM/edgeR log2-CPM (`prior.count=2`); per-gene z-score on training rows only.

Covariate regimes `raw` and `resid`. `resid` columns: `SMRIN`, `SMTSISCH`, `DTHHRDY` one-hot, `SEX`,
residualized on training rows only. `resid` is reported alongside `raw`; it is not a gate.

Within-site CV grouped by `SMNABTCH` (batch-held-out). Donor-held-out is asserted in the same folds
(one sample per donor). Metrics: Spearman ρ, donor-level permutation null (n_perm=200, seed `20260914`),
donor-bootstrap 95% percentile CI (B=200, seed `20260918`), calibrated R², uncalibrated R².
Pearson r is reported and is not a gate.

Transfer: B1→C1 and C1→B1 (`SMCENTER`). D1 is too small and is excluded, stated, not silently dropped.

**Stage 1 pass = ρ > 0 in both transfer directions with permutation-null ρ ≤ 0.05, in `raw`, ridge.**
If Stage 1 fails, STOP. Do not open GSE297234. The fibroblast ruler does not transfer between GTEx
sites, so a reprogramming projection would be uninterpretable.

GTEx public `AGE` is a 10-year bin, so ρ is against a coarse ordinal target and year-level error is
not meaningful. Gene-universe overlap with the DLPFC 25526-gene matrix (Ensembl, version-stripped)
is reported before Stage 2.


## Finished cells (from disk)

### Cohort
- SMTSD='Cells - Cultured fibroblasts' n_samples=652 n_donors=652 SMCENTER={'B1': 421, 'C1': 224, 'D1': 7} n_SMNABTCH=72

### Gene overlap
- n_genes_raw=59033 n_genes_filter=23485 n_overlap_dlpfc=17791 dlpfc_n=25526

### Within-site CV (SMNABTCH-grouped)
- site=B1 regime=raw method=ridge ρ=+0.549 (null -0.010 p=+0.005) CI=[+0.481, +0.604] calR²=+0.303 R²=+0.302 n=421 n_perm=200
- site=B1 regime=raw method=pls1 ρ=+0.183 (null +0.005 p=+0.005) CI=[+0.085, +0.271] calR²=+0.039 R²=+0.032 n=421 n_perm=200
- site=B1 regime=resid method=ridge ρ=+0.246 (null -0.016 p=+0.005) CI=[+0.162, +0.327] calR²=+0.072 R²=+0.071 n=421 n_perm=200
- site=B1 regime=resid method=pls1 ρ=-0.002 (null +0.003 p=+0.537) CI=[-0.088, +0.093] calR²=+0.003 R²=-0.030 n=421 n_perm=200
- site=C1 regime=raw method=ridge ρ=+0.355 (null -0.021 p=+0.005) CI=[+0.211, +0.464] calR²=+0.113 R²=+0.101 n=224 n_perm=200
- site=C1 regime=raw method=pls1 ρ=+0.139 (null +0.005 p=+0.040) CI=[-0.014, +0.257] calR²=+0.018 R²=-0.006 n=224 n_perm=200
- site=C1 regime=resid method=ridge ρ=+0.212 (null -0.024 p=+0.005) CI=[+0.090, +0.334] calR²=+0.031 R²=+0.017 n=224 n_perm=200
- site=C1 regime=resid method=pls1 ρ=+0.060 (null -0.001 p=+0.229) CI=[-0.093, +0.192] calR²=+0.010 R²=-0.018 n=224 n_perm=200

### Transfer B1↔C1
- B1→C1 regime=raw method=ridge ρ=+0.529 (null -0.004 p=+0.005) CI=[+0.418, +0.612] calR²=+0.353 R²=+0.289 pass_rho=True n_train=421 n_test=224 n_perm=200
- C1→B1 regime=raw method=ridge ρ=+0.561 (null +0.013 p=+0.005) CI=[+0.500, +0.617] calR²=+0.312 R²=+0.251 pass_rho=True n_train=224 n_test=421 n_perm=200
- B1→C1 regime=raw method=pls1 ρ=+0.218 (null +0.009 p=+0.005) CI=[+0.097, +0.324] calR²=+0.040 R²=-0.091 pass_rho=True n_train=421 n_test=224 n_perm=200
- C1→B1 regime=raw method=pls1 ρ=+0.236 (null -0.001 p=+0.005) CI=[+0.135, +0.317] calR²=+0.056 R²=-0.040 pass_rho=True n_train=224 n_test=421 n_perm=200
- B1→C1 regime=resid method=ridge ρ=+0.394 (null +0.016 p=+0.005) CI=[+0.281, +0.484] calR²=+0.187 R²=+0.064 pass_rho=True n_train=421 n_test=224 n_perm=200
- C1→B1 regime=resid method=ridge ρ=+0.430 (null +0.002 p=+0.005) CI=[+0.362, +0.498] calR²=+0.197 R²=+0.068 pass_rho=True n_train=224 n_test=421 n_perm=200
- B1→C1 regime=resid method=pls1 ρ=+0.146 (null +0.019 p=+0.030) CI=[+0.007, +0.276] calR²=+0.026 R²=-0.113 pass_rho=True n_train=421 n_test=224 n_perm=200
- C1→B1 regime=resid method=pls1 ρ=+0.066 (null -0.009 p=+0.179) CI=[-0.033, +0.152] calR²=+0.025 R²=-0.098 pass_rho=True n_train=224 n_test=421 n_perm=200

### Stage 1 gate
- pass=True method=ridge regime=raw B1→C1 ρ=+0.529 (null -0.004) C1→B1 ρ=+0.561 (null +0.013) reason=Stage 1 pass: ridge raw ρ>0 both ways with null ≤ 0.05 (B1→C1 ρ=+0.529 null=-0.004 p=+0.005; C1→B1 ρ=+0.561 null=+0.013 p=+0.005).

### Stage 2 (from disk)
- {"sanity": {"day0_aged": 2.768658769832732, "day0_young": -2.4229074198440266, "aged_gt_young": true}, "aged": {"cell_line": "GM00731", "days": [0, 3, 7, 10], "age_score": [2.768658769832732, -0.18832033320412486, -1.4159321966519975, 6.119542463860377], "pluri": [-1.5979352973543126, 8.174864864441993, 8.059927154669232, 3.845824200884512], "decline_0_to_10": -3.3508836940276447, "spearman_vs_day": 0.19999999999999998, "p_endpoint": 0.9900497512437811, "p_monotonicity": 0.5920398009950248, "n_random": 200, "n_null_ge_real": 198, "n_null_le_real_mono": 118}, "young": {"cell_line": "GM23815", "days": [0, 3, 7, 10], "age_score": [-2.4229074198440266, -4.992245093656465, -9.739773328221705, -4.884594873620896], "pluri": [-1.2646001062377195, 9.280109220766006, 10.14182083977496, 7.669076027518511], "decline_0_to_10": 2.46168745377687, "spearman_vs_day": -0.39999999999999997}, "pluri_primary": "with_OSKM", "pluri_disagree": false, "sendai": {"n_hits": 0, "names": [], "n_features_scanned": 73194, "n_features_first_library": 36601, "n_ensg": 36601, "n_non_ensg": 0, "n_libraries": 8, "genomes": ["GRCh38"], "feature_types": ["Gene Expression"], "columns_actually_read": ["gene_id", "symbol", "feature_type", "genome"], "pattern": "sendai|\\bsev\\b|\\bsev[-_]|kozak|transgene|cytotune|cyto.?tune|vector|orfeome|\\boskm\\b|sevoskm", "note": "Sendai/vector features reported; not used in the frozen GTEx ruler. Pluripotency is scored with and without OSKM-family endogenous genes. Zero hits means the Cell Ranger matrix has no labelled Sendai/vector features; vector reads can still map onto endogenous OSKM gene models."}, "verify": {"n_samples": 8, "days": [0, 3, 7, 10], "cell_lines": ["GM00731", "GM23815"], "age_map": {"GM00731": {"age_years": 96, "source": "Series_overall_design"}, "GM23815": {"age_years": 22, "source": "Series_overall_design"}}, "discrepancies": [], "used_record": true, "overall_design": "10x Genomics scRNA-seq of human fibroblasts from a young (GM23815, 22yr) and aged donor (GM00731, 96yr) treated with Sendai virus OSKM for up to 10 days."}, "outcome": {"key": "no_decline", "text": "age score does not decline → the ruler reads a static donor property, not a modifiable state. Step 3 is not supported by this data.", "decline": -3.3508836940276447, "monotonic": false, "drop_time": 3, "rise_time": 3, "beats_null": false, "declined": false, "p_endpoint": 0.9900497512437811, "sanity_ok": true}, "figure": "<repo>\\results\\fibro\\figures\\plane.png", "matrix_format": "10x_h5_filtered_feature_bc_matrix", "n_allcell": 8, "n_cluster_rows": 24, "cluster_skipped": [], "missing_genes_z0": 3086, "n_ruler_genes": 23485, "seed": 20260914, "n_perm_random": 200, "note": "nothing fitted on GSE297234", "spearman_note": "trajectory Spearman uses n>=3 (4 days); not gtex spearman_safe n>=8"}

## Note

<repo>\FINDINGS_FIBRO.md

## Failures (manifest)

- **download:** https://ftp.ncbi.nlm.nih.gov/geo/series/GSE297nnn/GSE297234/suppl/GSE297234_GM00731_SEVOSKM.rds: RuntimeError: failed to download https://ftp.ncbi.nlm.nih.gov/geo/series/GSE297nnn/GSE297234/suppl/GSE297234_GM00731_SEVOSKM.rds
- **download:** https://ftp.ncbi.nlm.nih.gov/geo/series/GSE297nnn/GSE297234/suppl/GSE297234_HFIB_COMBINED_SEVOSKM.rds: RuntimeError: failed to download https://ftp.ncbi.nlm.nih.gov/geo/series/GSE297nnn/GSE297234/suppl/GSE297234_HFIB_COMBINED_SEVOSKM.rds
- **10x:** no 10x matrix/features/barcodes triplets in GSE297234_RAW.tar. RDS present but R/pyreadr are not available; not substituting a fitted object.

## Files

- `src/fibro_common.py`, `fibro_stage1.py`, `fibro_stage2.py`, `fibro_findings.py`, `fibro_run.py`
- `results/fibro/`
- `FINDINGS_FIBRO.md`
- `PROGRESS_FIBRO.md`


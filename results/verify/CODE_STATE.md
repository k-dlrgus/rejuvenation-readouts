# CODE_STATE — which version of `src/pseudobulk.py` produced the stored results

Written 2026-10-02 from code, run logs, timestamps and stored outputs only. Nothing was re-run.

## Decision

**Mixed. Nothing was changed.** `src/pseudobulk.py` keeps its uncommitted modifications. No patch was written to `wip/` and `git checkout -- src/pseudobulk.py` was not run.

- The OneK1K pseudobulk, and everything built from it, was produced by the **committed** version.
- The tissue pseudobulks (brain white matter, skeletal muscle, retina) and the brain Aging_Cohort pseudobulk, and everything built from them, were produced by the **modified** version. Those two scripts cannot run on the committed version at all.

`src/download.py` and `src/rerun_common.py` were left as they are (their uncommitted changes only handle console encoding).

## The two versions

| | committed (`HEAD`, commit `cb31762`, 2026-09-14 03:03:49 +0900) | modified (working tree, last written 2026-09-14 10:05:17 +0900) |
| --- | --- | --- |
| `build()` signature | `..., celltype_map=None` | `..., celltype_map=None, require_primary=False, keep_mask=None` |
| caller cell mask | none | `keep_mask` applied first; logged as `caller keep_mask (analysis cohort)` |
| age parsing | `parse_age(development_stage)` | same, or `pd.to_numeric` when `age_col == "age"`; then cells with a non-numeric age are filled from `donor_age` if that column exists, printing `[age] filled N cells from donor_age ...` |
| primary-data filter | none | if `require_primary`, keep `is_primary_data == True`; logged as `is_primary_data==True` |

`cb31762` is the only commit that touches `src/pseudobulk.py`. The diff is +23 / −2 lines (`bundle/_build/uncommitted.diff`).

## Evidence, by output

### OneK1K (blood): committed version

- Caller: `src/run_pseudobulk_onek1k.py` calls `build(..., age_col="development_stage", ..., celltype_map=CELLTYPE_MAP)`. It uses no argument the committed version lacks.
- `data/processed/onek1k_pseudobulk.h5ad`, `results/phase1_pseudobulk_build.log` and `results/tables/onek1k_filter_log.csv` all carry the timestamp 2026-09-14 02:31:02. That is 33 minutes before commit `cb31762`. The commit added the committed `src/pseudobulk.py` together with that build log and filter log, and neither log has changed since (`git diff HEAD` is empty for both).
- The build log has no `[age] filled` line and no `caller keep_mask` or `is_primary_data==True` step. Its obs column list has `age` and `development_stage` but no `donor_age`, and the numeric-age filter removed 0 cells.
- From reading the code (not re-run): for this call the modified version takes the same path. `keep_mask` is None, `require_primary` is False, `age_col` is not `"age"`, there is no `donor_age` column and no age is missing. So it would write the same pseudobulk. Provenance is still the committed version.

Depends on this pseudobulk (through `data/processed/onek1k_logcpm.npy` and `onek1k_pseudobulk_obs.csv`, written 04:20:45):

- `src/run_1a_decomp.py`, `src/run_1b_confounds.py`, `src/run_1b_batch_verify.py`, `src/run_1b_depth_mediator.py`: `results/phase1a_*.txt`, `results/phase1b_*.txt`, `results/tables/phase1a_*`, `results/tables/phase1b_*`, `results/tables/onek1k_*`
- `src/rerun_stageA.py`, `src/rerun_stageB.py` (via `rerun_common.load_data`): `results/rerun/*`, `results/tables/rerun_stage*`
- `notebooks/phase1_scores.ipynb`, `notebooks/rerun_stageA_age_axis.ipynb`
- written up in `FINDINGS.md`, `FINDINGS_RERUN.md`, `FINDINGS_STAGEB.md`

### Tissue pseudobulks: modified version (an intermediate edit of it)

- Caller: `src/tissue_analyze.py` line 340 calls `build(..., require_primary=require_primary)`. The committed `build()` has no `require_primary` parameter and would raise `TypeError`, so every completed tissue build ran on a modified file.
- `results/tables/tissue_muscle_filter_log.csv` and `results/tables/tissue_brain_wm_filter_log.csv` contain the step `is_primary_data==True`, written by the modified `build()`. It comes right after `cells with numeric donor age`, the same order as in the current file. Muscle: 195,843 → 165,109 cells. White matter: 45,528 → 45,528.
- `results/tables/tissue_retina_sc_filter_log.csv` has no step unique to the modified file, because retina ran with `require_primary=False` (`results/tissue/retina_sc_analyze.txt` line 32: "NOT applied"). It still needed the modified signature.
- Timestamps: `tissue_brain_wm_pseudobulk.h5ad` 07:51:42, `tissue_muscle_pseudobulk.h5ad` 07:57:51, `tissue_retina_sc_pseudobulk.h5ad` 08:03:51. All three are **earlier** than the file's last write (10:05:17). So they ran on an intermediate edit that already had `require_primary`. Whether the `keep_mask` and `donor_age` hunks existed then cannot be told from the record. Neither would change these three outputs:
  - `tissue_analyze.py` passes no `keep_mask`.
  - Muscle has no `donor_age` column (obs column list in `results/tissue/muscle_analyze.txt`).
  - White matter dropped 0 cells for non-numeric age.
  - For retina, `tissue_analyze.py` ran the same `donor_age` fill on the same obs and logged "filled 0 cells" (`retina_sc_analyze.txt` line 33).

Depends on these: `results/tissue/*`, `results/tables/tissue_{brain_wm,muscle,retina_sc}_filter_log.csv`, `FINDINGS_TISSUE.md`, and `src/brain_phase1_p5.py` (reads `tissue_brain_wm_pseudobulk.h5ad`; `results/brain_phase1/p5_*`).

### Brain Aging_Cohort pseudobulk: modified version (the current file)

- Caller: `src/brain_analyze.py` lines 298–306 call `build(..., require_primary=False, keep_mask=keep_mask)`. This is impossible on the committed version.
- `results/tables/tissue_brain_aging_filter_log.csv` and `results/brain/filter_log.csv` begin with `caller keep_mask (analysis cohort),cells,1332155,1054845,...`. That step name exists only in the modified file. `results/brain/r1_r3_analyze.txt` line 40 logs the same step.
- `data/processed/tissue_brain_aging_pseudobulk.h5ad` 10:53:59, after the file's last write (10:05:17).

Depends on this:

- `results/brain/*` (`FINDINGS_BRAIN.md`)
- `src/brain_phase1_common.py` reads it (`PB_PATH`) to write `data/processed/brain_aging_phase1_logcpm*.{npy,csv}` (12:30:26). That feeds `results/brain_phase1/*`, `notebooks/brain_phase1_p*.ipynb` and `FINDINGS_BRAIN_PHASE1.md`, and the analyses that load the Phase 1 matrix through `geometry_common.load_phase1_matrix` and related loaders.
- The DLPFC gene list `brain_aging_phase1_logcpm_genes.csv` is used by `src/gtex_stage1.py` (BA9 transfer gene space), `src/gtex_stage0.py` (overlap counts) and `src/fibro_stage1.py`. In `fibro_stage1.py` it is used **only** for the overlap count written to `results/fibro/stage1_gene_overlap.json`. It does not restrict the genes of the frozen fibroblast ruler.

### Not dependent on either version's differences

`src/external_project.py`, `src/fibro2_ta.py` and `src/md2_task3.py` import only `find_counts_path` and `read_counts_chunk`. `src/profile_dataset.py` imports only `parse_age`. None of those functions changed. `paper_package/src/pseudobulk.py` is byte-identical to the modified working-tree file (same SHA-256), so the package snapshot already carries the modified version.

## What this means for a reader

- Restoring the committed file would break `src/tissue_analyze.py` and `src/brain_analyze.py` (`TypeError`), so those results could no longer be regenerated. By code reading, the modified file reproduces the OneK1K pseudobulk unchanged. That is why it stays in place.
- For a clean history, commit the modified `src/pseudobulk.py` together with the brain/tissue results, with a message noting that the OneK1K pseudobulk predates it. That decision is outside this task (nothing was committed).
- No preprint_10 table or figure source named in `README.md` sits under `results/rerun/`, `results/tissue/`, `results/brain/` or `results/brain_phase1/`. The preprint's fibroblast results do not depend on either version of `build()`.

# Rejuvenation readouts cannot distinguish younger from less old

Code, frozen instruments, pre-registrations and result tables for the preprint
*Rejuvenation readouts cannot distinguish younger from less old* (R. Kwon).

A one-dimensional age score falls for any move away from the old state. This repository
contains the analyses behind the claim that distance to a young target of the same cell
type tells approach from departure, and a one-dimensional score cannot.

Preprint: https://doi.org/10.5281/zenodo.23113818 (PDF also at [paper/rejuvenation_readouts_preprint_12.pdf](paper/rejuvenation_readouts_preprint_12.pdf))

Code and results archive (v1.0): https://doi.org/10.5281/zenodo.23092319
---

## Quick start

```bash
pip install -r requirements.txt
python src/download.py          # fetch the public inputs listed in Table 3
python src/figures_final.py     # re-plot Fig1-Fig5 and FigS1 from stored tables
```

Raw input data is not redistributed here. Every input is public; see Table 3 of the
preprint and `src/config.py` for the dataset registry. `src/figures_final.py` re-plots
from stored result tables only and refits nothing, so it runs without the raw data.

Scripts are run from the repository root and put the root on `sys.path` themselves via
`src/config.py`. Seeds are fixed: `20260914` for splits, permutations and null directions,
`20260918` for bootstraps.

---

## Layout

| Path | What is in it |
|---|---|
| `src/` | All analysis code. One entry point per analysis (table below). |
| `results/<name>/` | Output of the analysis named `<name>`, including its pre-registration flag. |
| `results/figures_final/` | The figures as they appear in the preprint, plus the plotting data for each. |
| `paper/` | paper/ holds the preprint PDF, rejuvenation_readouts_preprint_12.pdf. |
| `notebooks/` | Exploratory notebooks, outputs stripped. |
| `FINDINGS_*.md` | The written reading of each analysis, including withdrawn and superseded readings, marked as such. |
| `PROGRESS_*.md` | Running logs for each analysis. |

---

## Figures

Every figure in the preprint is produced by one script, `src/figures_final.py`, which
re-plots from stored tables and recomputes nothing. Each figure ships with the exact
table it was drawn from.

| Preprint | File | Plotting data |
|---|---|---|
| Figure 1, a published aging score falls while the cells move away from young cells | `results/figures_final/Fig1.png` | `Fig1_data.csv` |
| Figure 2, the result does not depend on which genes are read | `results/figures_final/Fig2.png` | `Fig2_data.csv` |
| Figure 3, transfer of the frozen fibroblast ruler | `results/figures_final/Fig3.png` | `Fig3_data.csv` |
| Figure 4, the frozen ruler against two curated instruments | `results/figures_final/Fig4.png` | `Fig4_data.csv` |
| Figure 5, a better age score makes the same mistake | `results/figures_final/Fig5.png` | `Fig5_data.csv` |
| Supplementary Figure S1, reprogramming states against a GTEx-defined youth axis, all genes | `results/figures_final/FigS1.png` | `FigS1_data.csv` |

```bash
python src/figures_final.py
```

## Tables

No script emits Tables 1 to 3 or Supplementary Table S1 as formatted tables. Each was
assembled by hand from the files below, and every cell traces to one of them.

| Preprint | Source file | Script |
|---|---|---|
| Table 1, rows 1-2 (GTEx site B1 to C1, and C1 to B1) | `results/fibro/stage1_transfer.csv` (`regime=raw`, `method=ridge`) | `src/fibro_run.py` |
| Table 1, row 3 (GSE226189) | `results/fibro2/ta_GSE226189_near_miss_bulk_result.json` | `src/fibro2_run.py` |
| Table 1, row 4 (atlas, 65 donors) | `results/paper_figs/cxg_true_donor_results.csv` (`analysis=true_donor`, `contrast=ruler_spearman`) | `src/paper_figs_cxg_donor.py` |
| Table 1, "652 donors (frozen)" | `results/fibro/freeze.json`, field `n_donors` | `src/fibro_run.py` |
| Table 2, row 1 (atlas, 65 donors), all columns | `results/paper_figs/cxg_true_donor_results.csv` (`analysis=true_donor`) | `src/paper_figs_cxg_donor.py` |
| Table 2, row 2 (GSE226189), instrument rho columns | `results/md3/t2_instruments.csv` (`gse226189_*` and `GSE226189_S3` rows) | `src/md3_run.py` |
| Table 2, row 2 (GSE226189), paired delta-rho | `results/md3/t2_delta_rho.csv` (`gse226189_*` rows) | `src/md3_run.py` |
| Table 3, dataset list | `src/config.py` registry, `src/gtex_common.py` URLs, `results/fibro2/ta_*_downloads.json`, `results/scoping/candidates.json`, `results/survey/coverage.csv` | several |
| Table 3, donor counts | `results/fibro/freeze.json` (652), `results/toward/anchors.json` (113 and 231), `results/paper_figs/cxg_true_donor_results.csv` (65), `results/fibro2/ta_GSE226189_near_miss_bulk_result.json` (82) | as above |
| Supplementary Table S1, changes in distance | `results/newstory/partB_gtex_distances.csv` (`cell_line=GM00731`; `delta_young`, `delta_old`) | `src/newstory.py` |
| Supplementary Table S1, lowest fold cosine | `results/genespace/gtex_c3.csv` (minimum `cos` per `space`) | `src/genespace_run.py` |

The atlas rows of `results/md3/t2_instruments.csv` and `results/md3/t2_delta_rho.csv` hold the
pre-registered run on 93 pseudobulks (ruler rho +0.455), which the preprint does not tabulate.
The 93-pseudobulk ruler-minus-MD delta-rho is superseded by the 65-donor value in Table 2. The
65-donor grouping was not pre-registered (Methods, "Analysis unit in the fibroblast atlas").

---

## Frozen instruments

Frozen before any evaluation and applied unchanged thereafter.

| Instrument | File | What it is |
|---|---|---|
| Fibroblast age ruler | `results/fibro/frozen_ruler_ridge_raw.npz` | Unit weight vector, un-normalised coefficients, training mean and standard deviation over 23,485 genes, from 652 GTEx fibroblast donors |
| GTEx young and old anchors | `results/toward/anchors.npz` | Mean expression of donors aged 20-39 and 60-79 |

The ruler orders donors. It is not a calibrated age in years, and is never used as an age
estimate.

---

## Entry points

Each analysis writes into `results/<name>/` and has one entry point.

### Analyses reported in the preprint

| Analysis | Command | What it does |
|---|---|---|
| FIBRO | `python src/fibro_run.py` | The frozen fibroblast age ruler and its within-GTEx site transfers |
| FIBRO2 | `python src/fibro2_run.py` | Transfer of the frozen ruler to GSE226189 and the CELLxGENE atlas (93 pseudobulks as run). GSE113957 was examined but not scored: its deposit is FPKM, not counts |
| FIBRO3 | `python src/fibro3_run.py` | Extrapolation check: ruler score against distance from the training manifold |
| SAME | `python src/same_run.py` | The same-platform test, with its forward, reverse and mixture controls |
| GENESPACE | `python src/genespace_run.py` | The same test repeated in the MD and age gene spaces |
| TOWARD | `python src/toward_run.py` | Do partially reprogrammed cells move toward a GTEx-defined young target |
| MD, MD2-MD5 | `python src/md_run.py`, `python src/md2_run.py` … | Reconstructing the mesenchymal-drift signature; instrument comparisons; the paired delta-rho bootstrap; cell-state signatures |
| NEWSTORY | `python src/newstory.py` | Plotting-data tables behind the narrative figures |
| PAPER_FIGS | `python src/paper_figs_cxg_donor.py`, `python src/paper_figs_manuscript.py` | The 65-donor regrouping of the atlas (not pre-registered), then the manuscript figure inputs |
| FINAL FIGURES | `python src/figures_final.py` | Fig1-Fig5 and FigS1 |

### Analyses in the repository but not in the preprint

Carried out in the same project and pre-registered in the same way, with any deviations
recorded alongside them. They do not bear on the question tested in the preprint and are
included so the record is complete. Several are declared negatives.

| Analysis | Command | Note |
|---|---|---|
| POSCTRL | `python src/posctrl_run.py` | Planted-effect positive controls for the brain (DLPFC) identifiability pipeline |
| SOUTH6 | `python src/south6_run.py` | The TF ranking: 1,836 transcription factors scored in 12 space by origin by target settings |
| SOUTH, SOUTH2-SOUTH5 | `python src/south_run.py` … | Earlier passes at the same screen. SOUTH halted at Stage 0; SOUTH5 stopped on a pre-registered STOP key |
| SENG, SENG_WHY | `python src/seng_run.py`, `python src/seng_why.py` | Senescence screen. Its reproduction gate failed 3 of 4 |
| LOWDIM | `python src/lowdim_run.py` | Stopped with `winner = NONE` |
| BRAIN, BRAIN_PHASE1, EXTERNAL | `python src/brain_r0.py`, `python src/brain_phase1_run.py`, `python src/external_run.py` | Cortex and DLPFC arm |
| GTEX, GEOMETRY, TARGET, TRAJECTORY, PLANE, RERUN, TISSUE, SCOPING, SURVEY, RELABEL_CHECK | see `src/` | Supporting and exploratory analyses |

Run `python src/<name>_run.py` for any analysis under `results/` that has a `src/<name>_run.py`.

---

## Pre-registration

Each gated analysis was specified in writing before any number was computed, including
its statistics, nulls, positive controls, decision thresholds and the exact wording of
each possible conclusion. A flag file was written before execution, and the
pre-registration text is reproduced verbatim in the results file for that analysis. Only
the pre-registered reading whose conditions fired is reported.

Where a later analysis superseded an earlier reading, the earlier text was kept and
marked withdrawn or superseded rather than edited. Nine readings across the project are
marked this way, three of which concern analyses reported in the preprint; all are listed in
`results/verify/WITHDRAWN_READINGS.md`.

These are self-administered protocols recorded here, not deposits with an external
registry. Diagnostics and descriptive summaries are not gated results.

Flags: `results/*/PREREG*.flag` and `results/*/PREREG_written_at.txt`.

---

## Environment

Python 3.12.7, numpy 2.5.3, scipy 1.18.1, pandas 3.0.5, scikit-learn 1.9.1,
anndata 0.13.3.post0, h5py 3.16.0, matplotlib 3.11.2. Full listing in
`requirements.txt`.

Seurat v5 and sctransform were not available in this environment. The `AddModuleScore`
reimplementation is a deviation from the published MD pipeline rather than a
substitution; see the Methods.

---

## Data

All inputs are public and listed in Table 3 of the preprint. Raw data is not
redistributed. `src/download.py` fetches it; `src/pseudobulk.py` rebuilds the derived
matrices.

GTEx is used from the open-access bucket only. No controlled-access tier is used
anywhere, and GTEx donor age is handled only as the public decade bin or its midpoint.

---

## Citation

R. Kwon. *Rejuvenation readouts cannot distinguish younger from less old.* Preprint.
https://doi.org/10.5281/zenodo.23092319

## Contact

kih7920@gmail.com

## License

See `LICENSE`.

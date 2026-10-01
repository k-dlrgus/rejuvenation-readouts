"""Shared paths / helpers for the TISSUE SELECTION stage (FINDINGS_TISSUE.md).

Everything is namespaced under results/tissue/ and data/raw/tissue/<name>/, data/processed/tissue_<name>_*; nothing from the
OneK1K runs is overwritten. Reused OneK1K code: src/pseudobulk.py (chunked donor x cell-type pseudobulk with a filter log),
src/decomp.normalize_logcpm (log-CPM + gene filter), src/rerun_stageB.py (kernel-ridge batch-grouped CV, within-batch R^2, CLR
composition, pool-cluster bootstrap), src/download.py (resumable downloader), src/search_cxg.py / src/geo_meta.py (metadata).
"""
import sys
from pathlib import Path

import numpy as np

# anndata 0.13 + xarray on this machine imports np.unicode_, removed in NumPy 2. Local workaround only.
if not hasattr(np, "unicode_"):
    np.unicode_ = np.str_  # type: ignore[attr-defined]

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import ROOT, DATA_RAW, DATA_PROC, RESULTS, SEED  # noqa: E402
from rerun_common import Logger  # noqa: E402,F401  (same console+file logger as the rerun)

TISSUE_DIR = RESULTS / "tissue"
TISSUE_FIG = TISSUE_DIR / "figures"
TISSUE_TAB = TISSUE_DIR  # tables live next to the T1 logs (results/tissue/*.csv)
TISSUE_RAW = DATA_RAW / "tissue"
TISSUE_PROC = DATA_PROC
for p in (TISSUE_DIR, TISSUE_FIG, TISSUE_RAW):
    p.mkdir(parents=True, exist_ok=True)

TISSUE_SEED = SEED  # 20260914 — same project seed; every rng is numpy.random.default_rng([TISSUE_SEED, ...])

# Blood Stage B reference numbers (FINDINGS_STAGEB.md). Do not recompute.
BLOOD_COMP_R2 = 0.493
BLOOD_EXPR_R2 = 0.482
BLOOD_RATIO = BLOOD_EXPR_R2 / BLOOD_COMP_R2  # 0.978 — expression does not beat composition

PRIORITY_TISSUES = ["retina", "brain", "skeletal muscle", "heart", "kidney", "liver"]

# Columns that, if present in obs, can identify a technical batch that is not the donor.
BATCH_COL_CANDIDATES = (
    "pool_number", "pool", "batch", "batch_id", "library_id", "library", "10x_lane",
    "suspension_uuid", "sample_id", "sample", "specimen", "specimen_id",
    "dataset", "dataset_id", "source", "study", "study_id", "geo_accession",
    "assay", "assay_ontology_term_id", "suspension_type",
    "tissue", "tissue_ontology_term_id",
    "development_stage",  # not a batch; listed so we do not confuse it
    "donor_id",
)

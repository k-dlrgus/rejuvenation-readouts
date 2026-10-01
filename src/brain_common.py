"""Shared paths / constants for the DLPFC brain-replication stage (FINDINGS_BRAIN.md).

Namespaced under results/brain/ and data/raw/tissue/brain_aging/. Does not overwrite
results/tissue/ or any prior FINDINGS_*.md. Reuses Stage B / T3 metric code in
src/tissue_analyze.py (kernel ridge, batch-grouped CV, within-batch R2, CLR, bootstrap).
"""
import sys
from pathlib import Path

import numpy as np

if not hasattr(np, "unicode_"):
    np.unicode_ = np.str_  # type: ignore[attr-defined]

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import ROOT, DATA_RAW, DATA_PROC, RESULTS, SEED  # noqa: E402
from rerun_common import Logger  # noqa: E402,F401
from tissue_common import BLOOD_COMP_R2, BLOOD_EXPR_R2  # noqa: E402
from tissue_download_top import REGISTRY  # noqa: E402

BRAIN_DIR = RESULTS / "brain"
BRAIN_FIG = BRAIN_DIR / "figures"
BRAIN_RAW = DATA_RAW / "tissue" / "brain_aging"
BRAIN_RETINA_RAW = DATA_RAW / "tissue" / "retina_sn"
for p in (BRAIN_DIR, BRAIN_FIG):
    p.mkdir(parents=True, exist_ok=True)

BRAIN_SEED = SEED  # 20260914

# Blood Stage B (FINDINGS_STAGEB.md) and white-matter T3 (FINDINGS_TISSUE.md). Do not recompute.
WM_COMP_R2 = -0.290
WM_EXPR_R2 = 0.194
WM_N_DONORS = 20
WM_N_BATCHES = 7

# Registry IDs were taken from results/tissue/t1_cxg_primary_candidates.csv (fetched, not invented).
DLPFC_NAME = "brain_aging"
DLPFC_DATASET_ID = REGISTRY[DLPFC_NAME]  # 4442d412-91cb-4261-acca-8adf5fa04c11
RETINA_SN_NAME = "retina_sn"
RETINA_SN_DATASET_ID = REGISTRY[RETINA_SN_NAME]  # d6505c89-c43d-4c28-8c4f-7351a5fd5528

# Adult cutoff: paper treats 12–19 as adolescence (developmental) and 20–39 as young adulthood;
# transcriptomic changes "stabilized after age 20". Developmental samples are not aging.
ADULT_MIN_AGE = 20.0
STOP_R1_MIN_DONORS = 40
STOP_R1_MIN_SPAN = 40.0

# Keywords that look like a sequencing pool / run / library (R0 decision).
TECH_BATCH_KEYS = (
    "pool", "batch", "library", "lane", "channel", "hash", "gem", "round",
    "10x", "chromium", "run_id", "seq_run", "sequencing", "flowcell", "chip",
)
# Shared-donor factors that are NOT a sequencing pool (still scored, not preferred).
SOFT_BATCH_KEYS = (
    "assay", "suspension", "tissue", "source", "dataset", "study", "institute",
    "brain_bank", "brainbank", "cohort", "site", "preservation", "platform",
)

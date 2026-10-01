"""Project-wide constants: paths, seeds, dataset registry.

All accessions/IDs below were fetched from GEO / CELLxGENE Discover APIs in Phase 0
(see results/phase0_*.txt); none were assumed.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA_RAW = ROOT / "data" / "raw"
DATA_PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"
FIG = RESULTS / "figures"
TAB = RESULTS / "tables"
for p in (DATA_RAW, DATA_PROC, FIG, TAB):
    p.mkdir(parents=True, exist_ok=True)

SEED = 20260914  # fixed project seed (YYYYMMDD of run start); recorded in FINDINGS.md

# --- Phase 0 registry -------------------------------------------------------
TF_ATLAS = dict(
    geo_superseries="GSE216481",
    geo_subseries="GSE217460",  # SHAREseq_210322_TFAtlas: the ~670k-cell hESC TF-overexpression atlas
    pubmed="36608654",
    raw_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217460/suppl/GSE217460_210322_TFAtlas_raw.h5ad.gz",
    small_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217460/suppl/GSE217460_210322_TFAtlas_differentiated_raw.h5ad.gz",
    local_dir=DATA_RAW / "GSE217460",
)

PRIMARY = dict(
    name="OneK1K",
    citation="Yazar et al. 2022 Science, doi:10.1126/science.abf3041",
    geo="GSE196830 (SuperSeries) / GSE196735 (scRNA-seq)",
    cxg_collection="dde06e0f-ab3b-46be-96a2-a8082383c4a1",
    cxg_dataset="3faad104-2ab8-4434-816d-474d8d2641db",
    h5ad_url="https://datasets.cellxgene.cziscience.com/1e44db10-b572-46cc-adae-dcc7acd44ca6.h5ad",
    local=DATA_RAW / "onek1k" / "1e44db10-b572-46cc-adae-dcc7acd44ca6.h5ad",
)

REPLICATION = dict(  # held out; NOT touched in Phase 1
    name="AIDA Phase 1 v2",
    citation="Asian Immune Diversity Atlas, doi:10.1016/j.cell.2025.02.017",
    cxg_collection="ced320a1-29f3-47c1-a735-513c7084d508",
    cxg_dataset="c838aec3-03ef-4398-b882-0e3912abfff0",
    h5ad_url="https://datasets.cellxgene.cziscience.com/f89a12c2-7a3b-415b-ab87-bbc550fe17f4.h5ad",
    local=DATA_RAW / "aida" / "f89a12c2-7a3b-415b-ab87-bbc550fe17f4.h5ad",
)

# --- Phase 1 parameters -----------------------------------------------------
MIN_CELLS_PER_PSEUDOBULK = 20     # (donor, cell type) groups with fewer cells are dropped (logged)
MIN_DONORS_PER_CELLTYPE = 100     # cell types present in fewer donors are dropped (logged)
N_OUTER_FOLDS = 5                 # donor-grouped outer CV

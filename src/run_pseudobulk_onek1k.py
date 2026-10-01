"""Build OneK1K donor x cell-type pseudobulks (Phase 1 input). Cached: skips if output exists."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import PRIMARY, DATA_PROC  # noqa
from pseudobulk import build  # noqa

# Azimuth L2 label -> analysis cell type. None = dropped (logged with counts).
# Rationale: keep the authors' fine-grained immune labels; merge 'NK Proliferating' into NK as the CL
# ontology does; drop non-PBMC/rare/ambiguous classes (platelet, erythrocyte, HSPC, doublet) and the two tiny
# proliferating T classes whose parent subtype is undefined (773 + 305 cells).
CELLTYPE_MAP = {
    "CD4 TCM": "CD4_TCM", "CD4 Naive": "CD4_Naive", "CD4 TEM": "CD4_TEM", "CD4 CTL": "CD4_CTL", "Treg": "Treg",
    "CD8 TEM": "CD8_TEM", "CD8 Naive": "CD8_Naive", "CD8 TCM": "CD8_TCM",
    "NK": "NK", "NK Proliferating": "NK", "NK_CD56bright": "NK_CD56bright",
    "B naive": "B_naive", "B memory": "B_memory", "B intermediate": "B_intermediate", "Plasmablast": "Plasmablast",
    "CD14 Mono": "CD14_Mono", "CD16 Mono": "CD16_Mono",
    "cDC2": "cDC", "cDC1": "cDC", "ASDC": "cDC", "pDC": "pDC",
    "gdT": "gdT", "MAIT": "MAIT", "dnT": "dnT", "ILC": "ILC",
    "CD4 Proliferating": None, "CD8 Proliferating": None,
    "HSPC": None, "Platelet": None, "Eryth": None, "Doublet": None,
}

if __name__ == "__main__":
    out = DATA_PROC / "onek1k_pseudobulk.h5ad"
    if out.exists():
        print(f"[cached] {out} exists; skipping build")
    else:
        build(PRIMARY["local"], "onek1k", donor_col="donor_id", ct_col="predicted.celltype.l2",
              age_col="development_stage", sex_col="sex", extra_obs=("pool_number", "self_reported_ethnicity"),
              celltype_map=CELLTYPE_MAP)

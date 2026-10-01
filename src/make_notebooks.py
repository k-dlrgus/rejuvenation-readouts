"""Generate the per-phase notebooks (notebooks/*.ipynb). Each notebook drives the src/ modules so that a clean
checkout reproduces every number in FINDINGS.md (downloads are cached in data/raw and skipped when complete)."""
import nbformat as nbf
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NB = ROOT / "notebooks"
NB.mkdir(exist_ok=True)

SETUP = """import os, sys, warnings
warnings.filterwarnings("ignore")
ROOT = os.path.abspath(os.path.join(os.getcwd(), "..")) if os.path.basename(os.getcwd()) == "notebooks" else os.getcwd()
os.chdir(ROOT); sys.path.insert(0, os.path.join(ROOT, "src"))
os.environ["PYTHONUTF8"] = "1"
from config import *
import subprocess
def sh(*args):
    \"\"\"run a src/ script in a subprocess and embed its stdout in the notebook (child stdout bypasses the kernel otherwise)\"\"\"
    r = subprocess.run([sys.executable, *args], capture_output=True, text=True, encoding="utf-8", errors="replace")
    print(r.stdout)
    if r.returncode:
        print(r.stderr); raise RuntimeError(f"{args[0]} exited with {r.returncode}")
print("ROOT:", ROOT, "| SEED:", SEED)"""


def md(s):
    return nbf.v4.new_markdown_cell(s)


def code(s):
    return nbf.v4.new_code_cell(s)


# ----------------------------------------------------------------------------- Phase 0
nb0 = nbf.v4.new_notebook()
nb0.cells = [
    md("""# Phase 0 — Data acquisition

**Hard rules enforced here**: no accession is assumed; every accession's metadata is fetched from GEO / CELLxGENE Discover and
printed before use. Downloads are cached in `data/raw/` and never repeated once complete.

Outputs: `results/phase0_*.txt`, `results/tables/phase0_*.csv`, `results/figures/phase0_*.png`."""),
    code(SETUP),
    md("""## 0.1 Confirmed accession: GSE216481 (Joung et al. 2023, TF Atlas)
GSE216481 is a **SuperSeries**; the ~670k-cell hESC TF-overexpression atlas is SubSeries **GSE217460** (`SHAREseq_210322_TFAtlas`).
We fetch metadata for the SuperSeries and all SubSeries, then download only the processed raw-count h5ad of GSE217460.
This dataset is **not analysed** in Phases 0–1 (reserved for Phase 2)."""),
    code("""from geo_meta import summarize
meta = summarize("GSE216481")
import requests
t = requests.get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216481&targ=self&form=text&view=full", timeout=60).text
subs = [l.split(": ")[-1] for l in t.splitlines() if l.startswith("!Series_relation = SuperSeries of")]
print("SubSeries:", subs)
_ = summarize("GSE217460")"""),
    code("""# download (cached). The full raw h5ad is 4.69 GB gz; NCBI throttles per connection, so a segmented downloader is used.
sh("src/download.py", TF_ATLAS["small_url"], "--out", str(TF_ATLAS["local_dir"]))
sh("src/download_parallel.py", TF_ATLAS["raw_url"], "--out", str(TF_ATLAS["local_dir"]), "--n", "24")"""),
    code("""# verify the small companion file loads (backed mode) — X is stored DENSE float32 in these files
sh("src/verify_h5ad.py", str(TF_ATLAS["local_dir"] / "GSE217460_210322_TFAtlas_differentiated_raw.h5ad.gz"))"""),
    code("""# verify the full raw file. Its X is DENSE float32: 1,145,823 cells x 37,528 genes = 172 GB decompressed (HDF5 superblock EOF
# address read from the gzip stream), which does not fit on disk. It is therefore converted in a single pass directly from the gzip
# stream into a sparse CSR h5ad (src/stream_convert_tf_atlas.py; ~11 min, <1 GB RAM). Prints shape and obs columns.
sh("src/stream_convert_tf_atlas.py")"""),
    md("""## 0.2 Search for human single-cell aging atlases
Criteria: per-donor chronological age, ≥10 donors, age span ≥40 y (ideally 20–80), ≥3 cell types, scRNA-seq, PBMC/immune preferred.
Two independent sources are queried and everything they return is recorded:
1. **CELLxGENE Discover** curation API (structured donor / development-stage metadata for 2,226 datasets).
2. **NCBI GEO** (`gds` esearch) with aging + single-cell + PBMC terms."""),
    code("""from search_cxg import main as cxg_search
cand = cxg_search()"""),
    code("""sh("src/search_geo.py")"""),
    code("""# top immune candidates with collection DOI and file size
import pandas as pd, requests
c = pd.read_csv(TAB / "phase0_cxg_candidates.csv")
ids = ["3faad104-2ab8-4434-816d-474d8d2641db","c838aec3-03ef-4398-b882-0e3912abfff0","ff4235bb-1cdf-4dcf-a13a-dc4f32244400",
       "ca7d95ac-53aa-4a94-8c76-300c608ce6a0","789ad837-9bf5-436c-9b60-9f1153546e4c","d86edd6a-4b5d-437a-ad80-1fd976a5e23a",
       "1b350d0a-4535-4879-beb6-1142f3f94947","53d208b0-2cfd-4366-9866-c3c6114081bc"]
for _, r in c[c.dataset_id.isin(ids)].iterrows():
    col = requests.get(f"https://api.cellxgene.cziscience.com/curation/v1/collections/{r.collection_id}", timeout=60).json()
    print(f"{r.dataset_id} | {col.get('name')[:70]} | DOI {col.get('doi')}\\n   cells={r.cell_count:,} donors={r.n_donors} cell_types={r.n_cell_types} "
          f"age={r.age_min:.0f}-{r.age_max:.0f} ({r.n_stages_numeric_years} distinct years) | {r.assays} | {r.tissues[:40]} | {r.diseases[:50]} | h5ad {(r.h5ad_bytes or 0)/1e9:.1f} GB")"""),
    md("""### Decision
| role | dataset | why |
|---|---|---|
| **Primary** | **OneK1K** (Yazar et al. 2022, *Science*; GEO GSE196830/GSE196735; CELLxGENE `3faad104…`) | 981 healthy donors, ages 19–97 (78 distinct years → continuous), 29 cell types, 1.25 M PBMCs, 10x 3' v2, 4.4 GB |
| **Replication (held out, untouched in Phase 1)** | **AIDA Phase 1 v2** (*Cell* 2025; CELLxGENE `c838aec3…`) | 625 healthy donors, 19–77, 32 cell types, 10x 5' v2, independent cohort/ancestry/chemistry |

Also seen but not chosen: Immunobiology of Aging Cohort (234 donors, 40–89, 40 GB, CMV status mixed), Indonesia PBMC (199 donors, 18–60 — span 42 y borderline),
Tabula Sapiens (24 donors), Immune-aging project (24 donors), GEO GSE157007 (Luo et al.: 3 young / 6 old / 5 frail → **binned**, would fire STOP 0)."""),
    code("""sh("src/download.py", PRIMARY["h5ad_url"], "--out", str(PRIMARY["local"].parent))
sh("src/verify_h5ad.py", str(PRIMARY["local"]))"""),
    code("""from profile_dataset import main as profile
donors = profile(str(PRIMARY["local"]), "onek1k")"""),
    code("""from IPython.display import Image, display
display(Image(str(FIG / "phase0_onek1k_donor_age_distribution.png")))"""),
    code("""# replication set: download only (NOT opened in Phase 1)
sh("src/download.py", REPLICATION["h5ad_url"], "--out", str(REPLICATION["local"].parent))"""),
]
nbf.write(nb0, NB / "phase0_data_acquisition.ipynb")

# ----------------------------------------------------------------------------- Phase 1
nb1 = nbf.v4.new_notebook()
nb1.cells = [
    md("""# Phase 1 — The two scores (1a variance decomposition, 1b confound removal)

Input: OneK1K (Phase 0). Unit of analysis: **pseudobulk per donor × cell type** (summed raw UMIs, log2 CPM).
Every filtering step is logged with counts (`results/tables/onek1k_filter_log.csv`). Seed: `config.SEED`.

**Outcome of this run: STOP CONDITION 1b fired — see `FINDINGS.md`.** 1c/1d were therefore not run."""),
    code(SETUP),
    md("""## 1.0 Pseudobulk construction
Cells → (donor, cell type) groups using the authors' Azimuth L2 labels (mapping in `src/run_pseudobulk_onek1k.py`).
Groups with <20 cells and cell types present in <100 donors are dropped (logged). Per-pseudobulk covariates recorded:
n_cells, mean UMI/cell, mean genes/cell, cell-cycle UMI fraction, fraction of cells expressing proliferation markers, pool, sex."""),
    code("""sh("src/run_pseudobulk_onek1k.py")
import pandas as pd
pd.read_csv(TAB / "onek1k_filter_log.csv")"""),
    md("""## 1a. Variance decomposition
Model (as specified): `expression ~ C(cell_type) + age`, OLS per gene on log2-CPM pseudobulks. Commonality partition:
`unique_ct = R²(full) − R²(age)`, `unique_age = R²(full) − R²(ct)`, `shared = R²(ct)+R²(age)−R²(full)`, `resid = 1 − R²(full)`.
Extras: interaction (`C(cell_type)*age`) and `partial_age` (fraction of within-cell-type variance explained by age).

**Why this model** (see FINDINGS.md §1a): pseudobulking removes cell-level Poisson noise and cell-level confounders, makes the
design almost balanced (every donor contributes most cell types, so `shared ≈ 0` and the partition is unambiguous), and puts
age and cell type on the same footing as *donor-level* vs *type-level* sources of variance. A per-cell model would let cell-type
composition and cell-level depth dominate the age term.

Threshold sweep: the grid was recalibrated after the first look at the distribution (initial a-priori grid gave empty A sets);
this happened before any downstream model existed and is recorded in `src/decomp.py` and `results/phase1a_decomp_initial.txt`."""),
    code("""from run_1a_decomp import main as run_1a
df = run_1a()"""),
    code("""from IPython.display import Image, display
display(Image(str(FIG / "phase1a_joint_variance_distribution.png")))
display(Image(str(FIG / "phase1a_mean_variance_partition.png")))"""),
    md("""## 1b. Confound removal
Flags: Tirosh S/G2M + proliferation list; data-driven cycling (|r| ≥ 0.3 with cell-cycle UMI fraction / proliferating-cell
fraction, within cell type); ribosomal; mitochondrial; chrX/chrY (Ensembl) + within-type sex DE; depth (|r| ≥ 0.3 with mean
UMI/cell or genes/cell, within cell type). Plus a **batch check**: age variance retained after adding the 10x `pool` to the model."""),
    code("""from run_1b_confounds import main as run_1b
F = run_1b()"""),
    code("""from run_1b_batch_verify import main as verify
verify()"""),
    code("""from run_1b_depth_mediator import main as mediator
mediator()"""),
    code("""display(Image(str(FIG / "phase1b_Agene_confound_filters.png")))
display(Image(str(FIG / "phase1b_batch_adjustment_age_signal.png")))
display(Image(str(FIG / "phase1b_pool_age_depth_confound.png")))"""),
    md("""## STOP CONDITION 1b — verdict
By the literal pre-specified flags 29/140 (20.7%) of primary A-genes are cell-cycle or technical. But the batch analysis shows the
A-gene set is dominated by a **technical artifact**: 136/140 (97%) retain <50% of their age variance after adjusting for 10x pool
(a real effect retains ≥ ~77%, confirmed by simulation), depth covariates alone do the same for 118/140 (84%), and within-pool the
A-genes' age correlation is below the permutation null. Mechanism: donors were not age-randomized across pools; per-cell UMI depth is a
pool property (R² = 0.86); depth correlates with age between pools (r = +0.43) and not within (r = −0.03).

Per the operating mode, Phase 1 halts here. 1c (score models) and 1d (negative controls) were **not run**. See `FINDINGS.md`."""),
    code("""from run_1b_verdict import main as verdict
fired = verdict()
assert fired, "verdict changed — FINDINGS.md must be revisited"
print("\\nPhase 1 halted at 1b. 1c/1d not executed.")"""),
]
nbf.write(nb1, NB / "phase1_scores.ipynb")
print("wrote", NB / "phase0_data_acquisition.ipynb", "and", NB / "phase1_scores.ipynb")

"""Generate notebooks/tissue_selection.ipynb — runnable end to end.

Drives the T1–T4 src/ modules. Downloads are cached. T2/T3 run on the three shortlisted
atlases (muscle, brain_wm, retina_sc) which are the smallest downloadable full atlases
that still meet the protocol cuts; the 104-donor retina snRNA and 291-donor DLPFC Aging_Cohort
are listed in T1 and can be swapped in by changing ANALYZE.
"""
import nbformat as nbf
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NB = ROOT / "notebooks"
NB.mkdir(exist_ok=True)

SETUP = r'''import os, sys, warnings
warnings.filterwarnings("ignore")
ROOT = os.path.abspath(os.path.join(os.getcwd(), "..")) if os.path.basename(os.getcwd()) == "notebooks" else os.getcwd()
os.chdir(ROOT); sys.path.insert(0, os.path.join(ROOT, "src"))
os.environ["PYTHONUTF8"] = "1"
from tissue_common import TISSUE_DIR, TISSUE_FIG, TISSUE_SEED, BLOOD_COMP_R2, BLOOD_EXPR_R2
import subprocess
def sh(*args):
    r = subprocess.run([sys.executable, *args], capture_output=True, text=True, encoding="utf-8", errors="replace")
    print(r.stdout)
    if r.returncode:
        print(r.stderr); raise RuntimeError(f"{args[0]} exited with {r.returncode}")
print("ROOT:", ROOT, "| TISSUE_SEED:", TISSUE_SEED)
print("blood Stage B reference: composition R2 =", BLOOD_COMP_R2, "expression R2 =", BLOOD_EXPR_R2)'''

ANALYZE = ["muscle", "brain_wm", "retina_sc"]


def md(s):
    return nbf.v4.new_markdown_cell(s)


def code(s):
    return nbf.v4.new_code_cell(s)


nb = nbf.v4.new_notebook()
nb.cells = [
    md("""# Tissue selection — find a low-turnover tissue where age is within-cell

Blood (OneK1K) is closed. Age in PBMC is at least half compositional (composition R² = 0.493 vs
expression R² = 0.482, `FINDINGS_STAGEB.md`). This notebook searches CELLxGENE Discover and GEO
for human single-cell / single-nucleus atlases in low-turnover tissues, audits batch recoverability,
and runs the Stage B comparison (composition vs within-type expression, batch-grouped CV) on the
three smallest atlases that still meet the protocol cuts.

**Hard rules.** No accession is invented — every ID is printed from a fetched record. CV is
batch-grouped. Filters are logged with counts. Seed `20260914`.

**Do not run Stage C/D on OneK1K. Do not modify prior FINDINGS files.**"""),
    code(SETUP),
    md("""## T1 — CELLxGENE Discover
Fetch the public curation listing (or reuse the cached `t1_cxg_all_human_datasets.csv` if present)
and recurate tissue groups so that kidney `cortex of kidney` is not scored as brain."""),
    code("""from pathlib import Path
if not (TISSUE_DIR / "t1_cxg_all_human_datasets.csv").exists():
    sh("src/tissue_search_cxg.py")
sh("src/tissue_curate.py")
import pandas as pd
c = pd.read_csv(TISSUE_DIR / "t1_cxg_primary_candidates.csv")
print(c.groupby("tissue_group").agg(n=("dataset_id","size"), full=("is_subset", lambda s: (~s).sum()),
                                    donors=("n_donors","max")).loc[["retina","brain","skeletal muscle","heart","kidney","liver"]])"""),
    md("""## T1 — GEO aging + reprogramming
NCBI E-utilities. Every accession below was returned by GEO; SOFT metadata is fetched for the
reprogramming-query hits."""),
    code("""sh("src/tissue_search_geo.py")
sh("src/tissue_search_reprogramming.py")"""),
    md("""## T1 shortlist
Primary (non-organoid, non-cancer, tissue-majority) full atlases, one row per tissue. The T2/T3
trio is the three smallest that still have ≥10 adult donors and a ≥40 y span: skeletal muscle
(HCA-harmonized, 32 donors), white-matter glia (20 donors), retina scRNA (20 donors). Larger
atlases (retina snRNA 104 donors / 37 GB; DLPFC Aging_Cohort 291 donors / 12 GB) remain the
preferred Phase-2 resources if download capacity allows — swap them into `ANALYZE`."""),
    code("""pd.set_option("display.width", 220); pd.set_option("display.max_colwidth", 72)
full = c[(~c.is_subset)].sort_values(["priority","n_donors"], ascending=[True, False])
print(full[["tissue_group","dataset_id","title","n_donors","cell_count","age_min","age_max",
            "assays","suspension_type","is_integrated","h5ad_bytes"]].to_string(index=False))"""),
    md("""## Download the T2/T3 trio (cached)"""),
    code(f"""ANALYZE = {ANALYZE!r}
for name in ANALYZE:
    sh("src/tissue_download_top.py", name)"""),
    md("""## T2 + T3 — batch audit and the turnover test
For each atlas: recover a shared-donor batch column from obs; R²(age ~ batch), R²(depth ~ batch),
age–depth between vs within batch; then the Stage B comparison (CLR composition vs all-type
expression, batch-grouped kernel ridge, within-batch R²). Blood numbers are plotted as reference."""),
    code("""for name in ANALYZE:
    sh("src/tissue_analyze.py", name)
sh("src/tissue_t4.py")
import pandas as pd
print("\\n=== T2 ===")
print(pd.read_csv(TISSUE_DIR / "t2_batch_audit.csv").to_string(index=False))
print("\\n=== T3 ===")
print(pd.read_csv(TISSUE_DIR / "t3_turnover.csv").to_string(index=False))
print("\\n=== T4 ===")
print(pd.read_csv(TISSUE_DIR / "t4_comparison.csv").to_string(index=False))"""),
    md("""## T4 — comparison table and recommendation
See `FINDINGS_TISSUE.md`. Rank by whether expression beats composition on the same within-batch
metric as blood (0.482 vs 0.493)."""),
    code("""from IPython.display import Image, display
from pathlib import Path
for p in sorted(TISSUE_FIG.glob("t3_*_comp_vs_expr.png")):
    print(p.name); display(Image(filename=str(p)))"""),
]

(NB / "tissue_selection.ipynb").write_text(nbf.writes(nb), encoding="utf-8")
print("wrote", NB / "tissue_selection.ipynb")


if __name__ == "__main__":
    pass

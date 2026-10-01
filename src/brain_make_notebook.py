"""Generate notebooks/brain_replication.ipynb — runnable end to end.

Drives src/brain_r0.py (obs peek, no matrix download) then src/tissue_download_top.py
and src/brain_analyze.py (R1–R3). Does not modify prior FINDINGS files.
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
from brain_common import BRAIN_DIR, BRAIN_FIG, BRAIN_SEED, BLOOD_COMP_R2, BLOOD_EXPR_R2, WM_COMP_R2, WM_EXPR_R2, ADULT_MIN_AGE
import subprocess, pandas as pd
from IPython.display import Image, display
def sh(*args):
    r = subprocess.run([sys.executable, *args], capture_output=True, text=True, encoding="utf-8", errors="replace")
    print(r.stdout)
    if r.returncode:
        print(r.stderr); raise RuntimeError(f"{args[0]} exited with {r.returncode}")
print("ROOT:", ROOT, "| BRAIN_SEED:", BRAIN_SEED, "| adult cutoff:", ADULT_MIN_AGE)
print("blood Stage B: composition", BLOOD_COMP_R2, "expression", BLOOD_EXPR_R2)
print("white matter T3: composition", WM_COMP_R2, "expression", WM_EXPR_R2)'''


def md(s):
    return nbf.v4.new_markdown_cell(s)


def code(s):
    return nbf.v4.new_code_cell(s)


nb = nbf.v4.new_notebook()
nb.cells = [
    md("""# Brain replication — DLPFC Aging_Cohort at scale

Blood (OneK1K) is closed: age is compositional (comp **0.493** vs expr **0.482**). White-matter glia
(20 donors) flipped that (comp **−0.290** vs expr **+0.194**). This notebook tests whether the flip
survives at scale on CELLxGENE `4442d412-91cb-4261-acca-8adf5fa04c11` (Aging_Cohort DLPFC, 291 donors,
1.33 M nuclei), using the same Stage B / T3 metric (batch-grouped kernel ridge, within-batch R²).

**Hard rules.** No accession is invented — the ID is read from `t1_cxg_primary_candidates.csv` and
re-verified against the live curation API. CV is batch-grouped. Filters are logged with counts.
Seed `20260914`. Adult cutoff **≥ 20** (adolescence is developmental).

**Do not modify prior FINDINGS files or FALSIFICATION.md.**"""),
    code(SETUP),
    md("""## R0 — batch peek (obs only, no 12 GB download)

Open the published h5ad with HTTP range requests. A usable batch is a sequencing pool / run / library
shared by multiple donors. If every candidate aliases donor 1:1, halt and peek the retina 104-donor
snRNA atlas instead."""),
    code("""sh("src/brain_r0.py", "brain_aging")
import json
from pathlib import Path
s = json.loads((BRAIN_DIR / "r0_summary.json").read_text(encoding="utf-8"))[0]
print("stop_r0:", s["stop_r0"], "batch_column:", s["batch_column"], "kind:", s["batch_kind"])
print("n_donors:", s["n_donors"], "n_cells:", s["n_cells"])
print(pd.read_csv(BRAIN_DIR / "r0_brain_aging_batch_columns.csv").to_string(index=False))"""),
    md("""## Download (cached) and R1–R3

R1 filters to adult (≥20) nuclei. This file is already the neurotypical Aging_Cohort subset of the
cross-disorder PsychAD atlas (psychiatric/neurodegeneration flags are constant). R2 audits
`Source` (HBCC vs MSSM) — there is no sequencing-pool column. R3 is the turnover test."""),
    code("""sh("src/tissue_download_top.py", "brain_aging")
sh("src/brain_analyze.py")
print("\\n=== R1 cohort ===")
print(pd.read_csv(BRAIN_DIR / "r1_cohort.csv").to_string(index=False))
print("\\n=== R2 batch audit ===")
print(pd.read_csv(BRAIN_DIR / "r2_batch_audit.csv").to_string(index=False))
print("\\n=== R3 summary ===")
print(pd.read_csv(BRAIN_DIR / "r3_summary.csv").to_string(index=False))
print("\\n=== filter log ===")
print(pd.read_csv(BRAIN_DIR / "filter_log.csv").to_string(index=False))"""),
    md("""## Comparison table and figures"""),
    code("""print("=== R4 comparison ===")
print(pd.read_csv(BRAIN_DIR / "r4_comparison.csv").to_string(index=False))
print("\\n=== per-type expression R2 ===")
print(pd.read_csv(BRAIN_DIR / "r3_per_type.csv").to_string(index=False))
for p in sorted(BRAIN_FIG.glob("*.png")):
    print(p)
    display(Image(filename=str(p)))"""),
]

out = NB / "brain_replication.ipynb"
out.write_text(nbf.writes(nb), encoding="utf-8")
print("wrote notebooks/brain_replication.ipynb")


if __name__ == "__main__":
    pass

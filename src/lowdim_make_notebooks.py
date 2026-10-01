"""Generate notebooks/lowdim_l1.ipynb … l4."""
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
from lowdim_common import LD_DIR, LD_FIG, LD_SEED, DATASET_ID, ANGLE_MATERIAL, KS
import subprocess, pandas as pd
from pathlib import Path
from IPython.display import Image, display
def sh(*args):
    r = subprocess.run([sys.executable, *args], capture_output=True, text=True, encoding="utf-8", errors="replace")
    print(r.stdout)
    if r.returncode:
        print(r.stderr); raise RuntimeError(f"{args[0]} exited with {r.returncode}")
print("ROOT:", ROOT, "| LD_SEED:", LD_SEED, "| dataset:", DATASET_ID)
print("k:", list(KS), "| material angle bar:", ANGLE_MATERIAL)'''


def md(s):
    return nbf.v4.new_markdown_cell(s)


def code(s):
    return nbf.v4.new_code_cell(s)


def write(name, cells):
    nb = nbf.v4.new_notebook()
    nb.cells = cells
    path = NB / name
    nbf.write(nb, path)
    print("wrote", path.name)


write("lowdim_l1.ipynb", [
    md("""# Low-dim L1 — build the spaces and flag technical components

Three variants (PCA, curated Hallmark/GO/Reactome modules, WGCNA-style co-expression),
k ∈ {20, 50, 100, 150}, on the 233-donor DLPFC pseudobulks.
Flag any component with |r| > 0.5 vs site, log depth, or neuronal fraction.
L1 fits are descriptive; L2 refits inside training folds.
Does not modify prior `FINDINGS*.md` or `FALSIFICATION.md`."""),
    code(SETUP),
    md("## Run L1"),
    code("""sh("src/lowdim_l1.py")
print((LD_DIR / "l1_report.txt").read_text(encoding="utf-8")[-3000:])
if (LD_DIR / "l1_overview.csv").exists():
    print(pd.read_csv(LD_DIR / "l1_overview.csv").to_string(index=False))"""),
    code("""for f in ("l1_frac_flagged.png", "l1_tech_r_hist_k50.png", "l1_per_type_frac_k50.png"):
    p = LD_FIG / f
    if p.exists():
        display(Image(str(p)))"""),
])

write("lowdim_l2.ipynb", [
    md("""# Low-dim L2 — refit the age direction and re-test identifiability

L2a donor-bootstrap angles (gene space 45.1°).
L2b held-out prediction, both CV schemes, permutation null.
L2c site transfer HBCC↔MSSM (gene space −0.213 / −0.499) — the single most important number.
L2d age-range split with bank control (gene space 86.1°).

Pick the best space by L2a and L2c together, not by R².
STOP if no variant/k is both more stable than gene space and has positive transfer
in at least one direction."""),
    code(SETUP),
    md("## Run L2"),
    code("""sh("src/lowdim_l2.py")
print((LD_DIR / "l2_STOP.txt").read_text(encoding="utf-8") if (LD_DIR / "l2_STOP.txt").exists() else "no STOP file")
if (LD_DIR / "l2_comparison.csv").exists():
    print(pd.read_csv(LD_DIR / "l2_comparison.csv").to_string(index=False))"""),
    code("""for f in ("l2a_bootstrap_vs_k.png", "l2b_within_site_r2.png", "l2c_transfer_vs_k.png", "l2d_agerange_vs_k.png"):
    p = LD_FIG / f
    if p.exists():
        display(Image(str(p)))"""),
])

write("lowdim_l3.ipynb", [
    md("""# Low-dim L3 — TARGET in the winning space

Identity-residual TARGET, V3 pairwise-angle retest vs null, consensus vs own,
gene-mapped top-50 weights and Enrichr. Skipped if STOP L2."""),
    code(SETUP),
    md("## Run L3"),
    code("""sh("src/lowdim_l3.py")
print((LD_DIR / "l3_STOP.txt").read_text(encoding="utf-8") if (LD_DIR / "l3_STOP.txt").exists() else "no skip file")
if (LD_DIR / "l3_geom_full.csv").exists():
    print(pd.read_csv(LD_DIR / "l3_geom_full.csv").head().to_string(index=False))
if (LD_DIR / "l3_top_genes.csv").exists():
    print(pd.read_csv(LD_DIR / "l3_top_genes.csv").query("source=='consensus'").head(20).to_string(index=False))"""),
    code("""for f in ("l3_pairwise_vs_null.png", "l3_top_weights.png"):
    p = LD_FIG / f
    if p.exists():
        display(Image(str(p)))"""),
])

write("lowdim_l4.ipynb", [
    md("""# Low-dim L4 — verdict

Writes `FINDINGS_LOWDIM.md` from L1–L3. Does not modify prior FINDINGS files
or `FALSIFICATION.md`."""),
    code(SETUP),
    md("## Write findings"),
    code("""sh("src/lowdim_findings.py")
print(Path("FINDINGS_LOWDIM.md").read_text(encoding="utf-8")[:6000])"""),
])

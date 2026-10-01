"""Generate notebooks/posctrl_c1.ipynb … posctrl_verdict.ipynb."""
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
from posctrl_common import PC_DIR, PC_FIG, PC_SEED, DATASET_ID, ANGLE_BAR, PREREG_BLOCK
import subprocess, pandas as pd
from pathlib import Path
from IPython.display import Image, display
def sh(*args):
    r = subprocess.run([sys.executable, *args], capture_output=True, text=True, encoding="utf-8", errors="replace")
    print(r.stdout)
    if r.returncode:
        print(r.stderr); raise RuntimeError(f"{args[0]} exited with {r.returncode}")
print("ROOT:", ROOT, "| PC_SEED:", PC_SEED, "| dataset:", DATASET_ID, "| angle bar:", ANGLE_BAR)
print("pre-registered block is in FINDINGS_POSCTRL.md (written before any model fit)")'''


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


write("posctrl_c2.ipynb", [
    md("""# POSCTRL C2 — sex control (gate)

C2a: all genes, target = sex ±1. Must pass bootstrap ≤35° and transfer >0.
If it fails: STOP, do not run C1/C3.
C2b: chrX, chrY, PAR removed.
Does not modify prior `FINDINGS*.md` or `FALSIFICATION.md`."""),
    code(SETUP),
    md("## Pre-registered interpretation (verbatim, written before fits)"),
    code("print(PREREG_BLOCK)"),
    md("## Run C2a then C2b"),
    code("""sh("src/posctrl_c2.py")
for name in ("c2a", "c2b"):
    p = PC_DIR / f"{name}.json"
    print(name, "exists" if p.exists() else "MISSING", p)
if (PC_DIR / "STOP_P1.txt").exists():
    print((PC_DIR / "STOP_P1.txt").read_text(encoding="utf-8"))"""),
])

write("posctrl_c1.ipynb", [
    md("""# POSCTRL C1 — planted-signal control

Permute age within bank, plant a unit direction at within-site ridge R² ≈ 0.30
(dense_shared / dense_typespec / sparse_shared / sparse_typespec), run the
gene-space L2 battery. dense_shared also at 0.15 and 0.50 (B/2).
Skipped if STOP P1."""),
    code(SETUP),
    md("## Run C1"),
    code("""sh("src/posctrl_c1.py")
print("beta:", (PC_DIR / "frozen_beta.json").read_text(encoding="utf-8")[:2000] if (PC_DIR / "frozen_beta.json").exists() else "no beta")
rows = []
for p in sorted(PC_DIR.glob("c1_*.json")):
    if p.name.endswith("_battery.json") or p.name.endswith("_bootstrap.json"):
        continue
    print(p.name)"""),
])

write("posctrl_c3.ipynb", [
    md("""# POSCTRL C3 — power curve

Subsample n ∈ {60, 100, 150, 233}, bank-stratified, 10 draws.
Real age (full B) and planted dense_shared R²≈0.30 (B/2).
Bootstrap angle and bank transfer vs n.
Skipped if STOP P1."""),
    code(SETUP),
    md("## Run C3"),
    code("""sh("src/posctrl_c3.py")
if (PC_DIR / "c3_summary.csv").exists():
    print(pd.read_csv(PC_DIR / "c3_summary.csv").to_string(index=False))
if (PC_DIR / "c3_fit.json").exists():
    print((PC_DIR / "c3_fit.json").read_text(encoding="utf-8")[:2000])"""),
    code("""for f in ("c3_angle_vs_n.png", "c3_transfer_vs_n.png"):
    p = PC_FIG / f
    if p.exists():
        display(Image(str(p)))"""),
])

write("posctrl_verdict.ipynb", [
    md("""# POSCTRL verdict

Writes `FINDINGS_POSCTRL.md` from `results/posctrl/*.json`.
Does not modify prior FINDINGS files or `FALSIFICATION.md`.
The pre-registered block is copied verbatim."""),
    code(SETUP),
    md("## Write findings"),
    code("""sh("src/posctrl_findings.py")
print(Path("FINDINGS_POSCTRL.md").read_text(encoding="utf-8")[:8000])"""),
])

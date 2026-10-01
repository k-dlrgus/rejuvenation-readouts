"""Generate notebooks/geometry_partA.ipynb, partB, partC."""
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
from geometry_common import GEO_DIR, GEO_FIG, GEO_SEED, NULL_USABLE, DATASET_ID
import subprocess, pandas as pd
from pathlib import Path
from IPython.display import Image, display, Markdown
def sh(*args):
    r = subprocess.run([sys.executable, *args], capture_output=True, text=True, encoding="utf-8", errors="replace")
    print(r.stdout)
    if r.returncode:
        print(r.stderr); raise RuntimeError(f"{args[0]} exited with {r.returncode}")
print("ROOT:", ROOT, "| GEO_SEED:", GEO_SEED, "| dataset:", DATASET_ID)
print("null usability line:", NULL_USABLE)'''


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


write("geometry_partA.ipynb", [
    md("""# Geometry Part A — leakage audit

Reproduce the within-site shuffle R², test every leakage hypothesis, fix what fires, re-score P3/P4c.
Halt if the corrected shuffle null is still above ~0.05. Does not modify `FALSIFICATION.md` or prior `FINDINGS*.md`."""),
    code(SETUP),
    md("## Run Part A"),
    code("""sh("src/geometry_partA.py")
print((GEO_DIR / "partA_STOP.txt").read_text(encoding="utf-8"))
print(pd.read_csv(GEO_DIR / "partA_hypotheses.csv").to_string(index=False))
print(pd.read_csv(GEO_DIR / "partA_side_by_side.csv").to_string(index=False))"""),
    code("""for f in ("partA_reported_vs_corrected.png", "partA_shuffle_null.png"):
    p = GEO_FIG / f
    if p.exists():
        display(Image(str(p)))"""),
])

write("geometry_partB.ipynb", [
    md("""# Geometry Part B — gene-count control

Random gene sets of size |A|, then a size sweep for A / I / random.
Skip if Part A stopped. Does not derive new A/I definitions."""),
    code(SETUP),
    md("## Run Part B"),
    code("""sh("src/geometry_partB.py")
print((GEO_DIR / "partB_B3.txt").read_text(encoding="utf-8") if (GEO_DIR / "partB_B3.txt").exists() else "skipped")
if (GEO_DIR / "partB_b1_summary.csv").exists():
    print(pd.read_csv(GEO_DIR / "partB_b1_summary.csv").to_string(index=False))"""),
    code("""for f in ("partB_random_nA.png", "partB_size_sweep.png"):
    p = GEO_FIG / f
    if p.exists():
        display(Image(str(p)))"""),
])

write("geometry_partC.ipynb", [
    md("""# Geometry Part C — age vs identity as directions

Abandon gene bucketing. Angle of the within-type age direction to the identity (centroid) subspace,
with a donor-age permutation null. Skip if Part A stopped."""),
    code(SETUP),
    md("## Run Part C"),
    code("""sh("src/geometry_partC.py")
print((GEO_DIR / "partC_claim.txt").read_text(encoding="utf-8") if (GEO_DIR / "partC_claim.txt").exists() else "skipped")
if (GEO_DIR / "partC_per_type_full_pls.csv").exists():
    print(pd.read_csv(GEO_DIR / "partC_per_type_full_pls.csv").to_string(index=False))"""),
    code("""for f in ("partC_angle_vs_null.png", "partC_per_type_angles.png",
          "partC_pairwise_angles.png", "partC_shared_pc1.png"):
    p = GEO_FIG / f
    if p.exists():
        display(Image(str(p)))"""),
])

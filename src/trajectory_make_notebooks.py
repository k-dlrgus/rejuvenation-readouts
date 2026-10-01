"""Generate notebooks/trajectory_g1.ipynb … g4."""
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
from trajectory_common import TRAJ_DIR, TRAJ_FIG, TRAJ_SEED, DATASET_ID, STOP_NULL
import subprocess, pandas as pd
from pathlib import Path
from IPython.display import Image, display, Markdown
def sh(*args):
    r = subprocess.run([sys.executable, *args], capture_output=True, text=True, encoding="utf-8", errors="replace")
    print(r.stdout)
    if r.returncode:
        print(r.stderr); raise RuntimeError(f"{args[0]} exited with {r.returncode}")
print("ROOT:", ROOT, "| TRAJ_SEED:", TRAJ_SEED, "| dataset:", DATASET_ID)
print("null usability line:", STOP_NULL)'''


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


write("trajectory_g1.ipynb", [
    md("""# Trajectory G1 — geometric age and identity scores

No gene selection. Age = projection onto a within-type supervised direction (PLS-1 and ridge).
Identity = position in the type-centroid subspace plus distance to own centroid.
Both CV schemes; permutation nulls. Halt if the corrected shuffle null is unusable or age R² collapses.
Does not modify `FALSIFICATION.md` or prior `FINDINGS*.md`."""),
    code(SETUP),
    md("## Run G1"),
    code("""sh("src/trajectory_g1.py")
print((TRAJ_DIR / "g1_STOP.txt").read_text(encoding="utf-8"))
print(pd.read_csv(TRAJ_DIR / "g1_summary_table.csv").to_string(index=False))"""),
    code("""for f in ("g1_age_vs_null.png", "g1_age_per_type.png", "g1_identity_acc.png"):
    p = TRAJ_FIG / f
    if p.exists():
        display(Image(str(p)))"""),
])

write("trajectory_g2.ipynb", [
    md("""# Trajectory G2 — orthogonality controls

G2a: angle with directions fit on training folds only.
G2b: can identity coordinates predict age?
G2c: age direction on identity-residual expression.
Halt before G3 if G2b is redundant or G2c collapses."""),
    code(SETUP),
    md("## Run G2"),
    code("""sh("src/trajectory_g2.py")
print((TRAJ_DIR / "g2_STOP.txt").read_text(encoding="utf-8") if (TRAJ_DIR / "g2_STOP.txt").exists() else "no STOP file")
if (TRAJ_DIR / "g2_summary_table.csv").exists():
    print(pd.read_csv(TRAJ_DIR / "g2_summary_table.csv").to_string(index=False))"""),
    code("""for f in ("g2a_angle_vs_null.png", "g2b_g2c_r2.png"):
    p = TRAJ_FIG / f
    if p.exists():
        display(Image(str(p)))"""),
])

write("trajectory_g3.ipynb", [
    md("""# Trajectory G3 — frozen scores on reprogramming data

Search and verify accessions first. Cross-species transfer is validated on an independent
mouse aging dataset before any reprogramming trajectory is plotted. Scores are frozen from G1;
nothing is refit."""),
    code(SETUP),
    md("## Search"),
    code("""sh("src/trajectory_g3_search.py")
print((TRAJ_DIR / "g3_search.txt").read_text(encoding="utf-8")[-4000:])"""),
    md("## Transfer + trajectory (skipped if G1/G2 halted or transfer fails)"),
    code("""sh("src/trajectory_g3.py")
for p in sorted(TRAJ_DIR.glob("g3_*.txt")):
    print("====", p.name, "====")
    print(p.read_text(encoding="utf-8")[-2500:])"""),
    code("""for f in ("g3_mouse_aging_validation.png", "g3_trajectory.png", "g3_didt.png"):
    p = TRAJ_FIG / f
    if p.exists():
        display(Image(str(p)))"""),
])

write("trajectory_g4.ipynb", [
    md("""# Trajectory G4 — verdict

Writes `FINDINGS_TRAJECTORY.md` from G1–G3 results. Does not modify prior FINDINGS files
or `FALSIFICATION.md`."""),
    code(SETUP),
    md("## Write findings"),
    code("""sh("src/trajectory_findings.py")
print(Path("FINDINGS_TRAJECTORY.md").read_text(encoding="utf-8")[:4000])"""),
])

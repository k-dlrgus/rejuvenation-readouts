"""Generate notebooks/target_v1.ipynb … v5."""
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
from target_common import TARGET_DIR, TARGET_FIG, TARGET_SEED, DATASET_ID, STOP_NULL, STOP_BOOT_ANGLE
import subprocess, pandas as pd
from pathlib import Path
from IPython.display import Image, display
def sh(*args):
    r = subprocess.run([sys.executable, *args], capture_output=True, text=True, encoding="utf-8", errors="replace")
    print(r.stdout)
    if r.returncode:
        print(r.stderr); raise RuntimeError(f"{args[0]} exited with {r.returncode}")
print("ROOT:", ROOT, "| TARGET_SEED:", TARGET_SEED, "| dataset:", DATASET_ID)
print("null usability line:", STOP_NULL, "| STOP V1 bootstrap angle:", STOP_BOOT_ANGLE)'''


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


write("target_v1.ipynb", [
    md("""# Target V1 — recover the axes and bootstrap stability

Recompute, per cell type, the ridge (and PLS-1) age direction and the identity
centroid subspace on the adult-control DLPFC cohort. Bootstrap donors; halt if
the median pairwise angle between ridge replicates is above ~60°.
No gene selection. Does not modify `FALSIFICATION.md` or prior `FINDINGS*.md`."""),
    code(SETUP),
    md("## Run V1"),
    code("""sh("src/target_v1.py")
print((TARGET_DIR / "v1_STOP.txt").read_text(encoding="utf-8"))
print(pd.read_csv(TARGET_DIR / "v1_bootstrap_summary.csv").to_string(index=False))"""),
    code("""for f in ("v1_bootstrap_per_type.png", "v1_bootstrap_hist.png"):
    p = TARGET_FIG / f
    if p.exists():
        display(Image(str(p)))"""),
])

write("target_v2.ipynb", [
    md("""# Target V2 — the target vector

TARGET = the age direction with its identity-subspace component removed, then
renormalized. Report norm retained, angle to the raw age direction, and whether
TARGET still predicts age held out (both CV schemes, permutation null).
Halt if a majority of cell types lose age predictivity vs null."""),
    code(SETUP),
    md("## Run V2"),
    code("""sh("src/target_v2.py")
print((TARGET_DIR / "v2_STOP.txt").read_text(encoding="utf-8") if (TARGET_DIR / "v2_STOP.txt").exists() else "no STOP file")
if (TARGET_DIR / "v2_geom_full.csv").exists():
    print(pd.read_csv(TARGET_DIR / "v2_geom_full.csv").query("method=='ridge'").to_string(index=False))"""),
    code("""for f in ("v2_norm_retained.png", "v2_per_type_r2.png", "v2_r2_vs_age.png"):
    p = TARGET_FIG / f
    if p.exists():
        display(Image(str(p)))"""),
])

write("target_v3.ipynb", [
    md("""# Target V3 — one target, or twenty?

Pairwise angles of per-type TARGETs vs a permuted-age null. Consensus = PC1 of
those vectors. Held-out: consensus vs own-type TARGET R².
(a) one shared target, or (b) type-specific — stated plainly."""),
    code(SETUP),
    md("## Run V3"),
    code("""sh("src/target_v3.py")
print((TARGET_DIR / "v3_verdict.txt").read_text(encoding="utf-8") if (TARGET_DIR / "v3_verdict.txt").exists() else "no verdict")
if (TARGET_DIR / "v3c_within_site_per_type.csv").exists():
    print(pd.read_csv(TARGET_DIR / "v3c_within_site_per_type.csv").to_string(index=False))"""),
    code("""for f in ("v3_pairwise.png", "v3_pairwise_pc1_vs_null.png", "v3_consensus_vs_own.png"):
    p = TARGET_FIG / f
    if p.exists():
        display(Image(str(p)))"""),
])

write("target_v4.ipynb", [
    md("""# Target V4 — how trustworthy is it?

Site transfer (HBCC ↔ MSSM), age-range split (20–50 vs 50–89), and biological
sanity of the largest weights (Enrichr; interpretation only, not gene selection)."""),
    code(SETUP),
    md("## Run V4"),
    code("""sh("src/target_v4.py")
if (TARGET_DIR / "v4a_summary.csv").exists():
    print(pd.read_csv(TARGET_DIR / "v4a_summary.csv").to_string(index=False))
if (TARGET_DIR / "v4b_per_type_angles.csv").exists():
    print(pd.read_csv(TARGET_DIR / "v4b_per_type_angles.csv").to_string(index=False))"""),
    code("""for f in ("v4_site_transfer.png", "v4_age_range_angles.png", "v4_enrichr.png"):
    p = TARGET_FIG / f
    if p.exists():
        display(Image(str(p)))"""),
])

write("target_v5.ipynb", [
    md("""# Target V5 — verdict

Writes `FINDINGS_TARGET.md` from V1–V4. Does not modify prior FINDINGS files
or `FALSIFICATION.md`."""),
    code(SETUP),
    md("## Write findings"),
    code("""sh("src/target_findings.py")
print(Path("FINDINGS_TARGET.md").read_text(encoding="utf-8")[:5000])"""),
])

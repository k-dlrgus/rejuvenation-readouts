"""Generate notebooks/brain_phase1_p0.ipynb … p5 — each drives the corresponding src script."""
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
from brain_phase1_common import PHASE1_DIR, PHASE1_FIG, PHASE1_SEED, PRIMARY_QUANTILES, P0_AGE_CORR_ABS
import subprocess, pandas as pd
from pathlib import Path
from IPython.display import Image, display, Markdown
def sh(*args):
    r = subprocess.run([sys.executable, *args], capture_output=True, text=True, encoding="utf-8", errors="replace")
    print(r.stdout)
    if r.returncode:
        print(r.stderr); raise RuntimeError(f"{args[0]} exited with {r.returncode}")
print("ROOT:", ROOT, "| PHASE1_SEED:", PHASE1_SEED)
print("declared quantiles:", PRIMARY_QUANTILES, "| P0 |r| threshold:", P0_AGE_CORR_ABS)'''


def md(s):
    return nbf.v4.new_markdown_cell(s)


def code(s):
    return nbf.v4.new_code_cell(s)


def write(name, cells):
    nb = nbf.v4.new_notebook()
    nb.cells = cells
    path = NB / name
    nbf.write(nb, path)
    print("wrote", path)


write("brain_phase1_p0.ipynb", [
    md("""# Brain Phase 1 — P0 confounders (DLPFC adult controls)

Audit PMI, RIN, ambient-RNA proxies, neuronal nucleus fraction, cause of death, tissue pH.
Hashing pools are not in obs. Do not modify `FALSIFICATION.md` or prior `FINDINGS*.md`."""),
    code(SETUP),
    md("## Run P0"),
    code("""sh("src/brain_phase1_p0.py")
print(pd.read_csv(PHASE1_DIR / "p0_covariate_age_correlations.csv").to_string(index=False))
print("\\nincluded:", (PHASE1_DIR / "p0_included_covariates.json").read_text(encoding="utf-8")[:1500])"""),
    code("""p = PHASE1_FIG / "p0_covariates_vs_age.png"
if p.exists():
    display(Image(str(p)))"""),
])

write("brain_phase1_p1.ipynb", [
    md("""# Brain Phase 1 — P1 variance decomposition

Quantile thresholds are declared **before any score**. Primary: A = Q90 unique_age & ≤ Q50 unique_ct;
I = Q90 unique_ct & ≤ Q50 unique_age; MIXED = both high."""),
    code(SETUP),
    code("""sh("src/brain_phase1_p1.py")
import json
print(json.loads((PHASE1_DIR / "p1_thresholds.json").read_text(encoding="utf-8")))
print(pd.read_csv(PHASE1_DIR / "p1_threshold_sweep_set_sizes.csv").to_string(index=False))"""),
    code("""for f in ("p1_unique_age_vs_unique_ct.png", "p1_unique_age_hist_vs_blood.png"):
    p = PHASE1_FIG / f
    if p.exists():
        display(Image(str(p)))"""),
])

write("brain_phase1_p2.ipynb", [
    md("""# Brain Phase 1 — P2 confound filters

STOP if fewer than 30 A-genes survive. Cell-cycle is expected to bite less in post-mitotic DLPFC."""),
    code(SETUP),
    code("""sh("src/brain_phase1_p2.py")
print(pd.read_csv(PHASE1_DIR / "p2_Agene_survival.csv").to_string(index=False))
print((PHASE1_DIR / "p2_STOP_verdict.txt").read_text(encoding="utf-8"))"""),
    code("""p = PHASE1_FIG / "p2_Agene_survival.png"
if p.exists():
    display(Image(str(p)))"""),
])

write("brain_phase1_p3.ipynb", [
    md("""# Brain Phase 1 — P3 scores (both CV schemes)

Leave-one-site-out **and** donor-grouped k-fold within each site. Nested gene selection.
Unrestricted within-type kernel ridge is the ceiling."""),
    code(SETUP),
    code("""sh("src/brain_phase1_p3.py")
print(pd.read_csv(PHASE1_DIR / "p3_identity.csv").to_string(index=False))
print(pd.read_csv(PHASE1_DIR / "p3_age_per_type.csv").head(40).to_string(index=False))
print((PHASE1_DIR / "p3_STOP_verdict.txt").read_text(encoding="utf-8"))"""),
    code("""p = PHASE1_FIG / "p3_age_A_vs_unrestricted.png"
if p.exists():
    display(Image(str(p)))"""),
])

write("brain_phase1_p4.ipynb", [
    md("""# Brain Phase 1 — P4 negatives and cross-test

Governed by `FALSIFICATION.md` (unchanged). Either outcome is the result; do not retune."""),
    code(SETUP),
    code("""sh("src/brain_phase1_p4.py")
print((PHASE1_DIR / "p4_verdict.txt").read_text(encoding="utf-8"))"""),
    code("""for f in ("p4_cross_test_age.png", "p4_cross_test_identity.png"):
    p = PHASE1_FIG / f
    if p.exists():
        display(Image(str(p)))"""),
])

write("brain_phase1_p5.ipynb", [
    md("""# Brain Phase 1 — P5 frozen sets on white-matter atlas

CELLxGENE `c05e6940-729c-47bd-a2a6-6ce3730c4919`. Do not re-derive genes."""),
    code(SETUP),
    code("""sh("src/brain_phase1_p5.py")
import json
print(json.dumps(json.loads((PHASE1_DIR / "p5_summary.json").read_text(encoding="utf-8")), indent=2)[:4000])"""),
    code("""p = PHASE1_FIG / "p5_wm_cross_test.png"
if p.exists():
    display(Image(str(p)))"""),
])

if __name__ == "__main__":
    print("done")

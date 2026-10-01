"""Generate notebooks/gtex_stage0.ipynb and gtex_stage1.ipynb."""
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
from gtex_common import GTEX_DIR, GTEX_FIG, GTEX_SEED, RELEASE, stage1_gate
import subprocess, pandas as pd
from pathlib import Path
from IPython.display import Image, display, Markdown
def sh(*args):
    r = subprocess.run([sys.executable, *args], capture_output=True, text=True, encoding="utf-8", errors="replace")
    print(r.stdout)
    if r.returncode:
        print(r.stderr); raise RuntimeError(f"{args[0]} exited with {r.returncode}")
print("ROOT:", ROOT, "| GTEX_SEED:", GTEX_SEED, "| release:", RELEASE)
print("Stage 1 gate:", stage1_gate())'''


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


write("gtex_stage0.ipynb", [
    md("""# GTEx Stage 0 — download, QC, covariate audit

GTEx v10 open-access: gene read counts GCT (stream and subset), sample attributes,
subject phenotypes. Six tissues, RNA-seq analysis freeze only.
Does not modify prior `FINDINGS*.md` or `FALSIFICATION.md`.
Stage 1 is not run from this notebook."""),
    code(SETUP),
    md("## Frozen tissues and markers"),
    code("""from gtex_common import TISSUES, NEURONAL_MARKERS, GLIAL_MARKERS, SMAFRZE_RNASEQ
print("freeze:", SMAFRZE_RNASEQ)
print("tissues:")
for t in TISSUES:
    print(" ", t)
print("neuronal:", NEURONAL_MARKERS)
print("glial (frozen a priori):", GLIAL_MARKERS)"""),
    md("## Run Stage 0"),
    code("""sh("src/gtex_stage0.py")
print((GTEX_DIR / "stage0_report.txt").read_text(encoding="utf-8")[-4000:])"""),
    md("## Cohort and covariate audit"),
    code("""print((GTEX_DIR / "stage0_cohort.md").read_text(encoding="utf-8")[:6000] if (GTEX_DIR / "stage0_cohort.md").exists() else "no cohort md")
if (GTEX_DIR / "stage0_covariate_audit.csv").exists():
    print(pd.read_csv(GTEX_DIR / "stage0_covariate_audit.csv").to_string(index=False))
if (GTEX_DIR / "stage0_gene_overlap.csv").exists():
    print(pd.read_csv(GTEX_DIR / "stage0_gene_overlap.csv").to_string(index=False))
if (GTEX_DIR / "stage0_composition.csv").exists():
    print(pd.read_csv(GTEX_DIR / "stage0_composition.csv").to_string(index=False))"""),
    code("""for f in ("stage0_age_bins.png", "stage0_rin_vs_age.png", "stage0_isch_vs_age.png", "stage0_composition_vs_age.png"):
    p = GTEX_FIG / f
    if p.exists():
        display(Image(str(p)))"""),
    md("## Manifest and findings"),
    code("""sh("src/gtex_findings.py")
import json
print(json.dumps(json.loads((GTEX_DIR / "manifest.json").read_text(encoding="utf-8")), indent=1)[:4000])
print(Path("FINDINGS_GTEX.md").read_text(encoding="utf-8")[:5000])"""),
])

write("gtex_stage1.ipynb", [
    md("""# GTEx Stage 1 — identifiability battery (gated)

Runs only if FINDINGS_POSCTRL.md Supersession 2 fired outcome (i) or (ii)
(planted dense_shared transfers as Pearson r). Primary transfer metric is held-out
Pearson r (`n_perm=200`). Pass = r > 0 with null ≤ 0.05. Bootstrap angle is reported
only. Pre-registration 2026-09-17 is written before any T-cell is fit.
Does not modify prior `FINDINGS*.md` or `FALSIFICATION.md`."""),
    code(SETUP),
    md("## Gate and pre-registration (before any fit)"),
    code("""g = stage1_gate()
print(g)
from gtex_common import PREREG_FLAG, PREREG_BARS, PREREG_READING
print("PREREG flag:", PREREG_FLAG.exists(), PREREG_FLAG)
print(PREREG_BARS)
print(PREREG_READING)
if not g["open"]:
    print("STAGE 1 not opened. Do not fit T-cells.")"""),
    md("## Write FINDINGS with frozen pre-registration (no T-cell fit)"),
    code("""sh("src/gtex_findings.py")
print("PREREG flag now:", (GTEX_DIR / "PREREG_20260917.flag").exists())
print(Path("FINDINGS_GTEX.md").read_text(encoding="utf-8")[:4000])"""),
    md("## T1 — bootstrap angles + within-batch scores"),
    code("""sh("src/gtex_stage1.py", "--cell", "T1")
print(Path("PROGRESS_GTEX.md").read_text(encoding="utf-8")[-2500:])"""),
    md("## T2 — center transfer (skipped if <3 SMCENTER with ≥25 BA9)"),
    code("""sh("src/gtex_stage1.py", "--cell", "T2")
print(Path("PROGRESS_GTEX.md").read_text(encoding="utf-8")[-2000:])"""),
    md("## T3 — BA9 ↔ Cortex (donors disjoint)"),
    code("""sh("src/gtex_stage1.py", "--cell", "T3")
print(Path("PROGRESS_GTEX.md").read_text(encoding="utf-8")[-2000:])"""),
    md("## T4 — BA9 → BA24 / hippocampus"),
    code("""sh("src/gtex_stage1.py", "--cell", "T4")
print(Path("PROGRESS_GTEX.md").read_text(encoding="utf-8")[-2000:])"""),
    md("## T5 — BA9 age direction vs tissue-identity subspace"),
    code("""sh("src/gtex_stage1.py", "--cell", "T5")
print(Path("PROGRESS_GTEX.md").read_text(encoding="utf-8")[-2000:])"""),
    md("## T6 — GTEx BA9 vs DLPFC age direction in overlap"),
    code("""sh("src/gtex_stage1.py", "--cell", "T6")
print(Path("PROGRESS_GTEX.md").read_text(encoding="utf-8")[-2000:])"""),
    md("## Composition-only T1 if BA9 raw T3 passed; verdict"),
    code("""sh("src/gtex_stage1.py", "--cell", "COMP")
print(Path("PROGRESS_GTEX.md").read_text(encoding="utf-8")[-2500:])"""),
    md("## Write findings"),
    code("""sh("src/gtex_findings.py")
print(Path("FINDINGS_GTEX.md").read_text(encoding="utf-8")[-6000:])"""),
])

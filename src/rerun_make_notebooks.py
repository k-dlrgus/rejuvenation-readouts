"""Generate the RERUN notebooks (notebooks/rerun_*.ipynb). Each notebook drives src/rerun_*.py so that a clean checkout
reproduces every number in FINDINGS_RERUN.md. Only the stages that were actually reached are generated: the rerun halted at
STOP A, so there is one notebook (Stage A). Stage B-D notebooks are deliberately absent — see FINDINGS_RERUN.md."""
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
from rerun_common import RERUN_DIR, RERUN_SEED
import subprocess
def sh(*args):
    \"\"\"run a src/ script in a subprocess and embed its stdout in the notebook\"\"\"
    r = subprocess.run([sys.executable, *args], capture_output=True, text=True, encoding="utf-8", errors="replace")
    print(r.stdout)
    if r.returncode:
        print(r.stderr); raise RuntimeError(f"{args[0]} exited with {r.returncode}")
print("ROOT:", ROOT, "| RERUN_SEED:", RERUN_SEED)"""


def md(s):
    return nbf.v4.new_markdown_cell(s)


def code(s):
    return nbf.v4.new_code_cell(s)


nbA = nbf.v4.new_notebook()
nbA.cells = [
    md("""# Phase 1 RERUN — Stage A: re-derive the age axis within batch

**Context.** The prior run (`FINDINGS.md`) halted at STOP 1b: under `expression ~ cell_type + age`, OneK1K's A-genes were
depth-artifact genes, because donors were not age-randomised across the 75 multiplexed 10x pools and per-cell depth is a pool
property (R² = 0.86). This rerun discards the prior A-gene set and re-derives Phase 1 from 1a under the corrected model

`expression ~ C(cell_type) + C(pool) + log(mean UMI/cell) + log(mean genes/cell) + age`

Reused from the prior run: cached raw data (`data/raw/`), the pseudobulk construction (`src/pseudobulk.py`,
`src/run_pseudobulk_onek1k.py`), the log-CPM normalisation (`decomp.normalize_logcpm`) and the confound gene lists / flags
(`src/gene_lists.py`, `src/confounds.py`). Thresholds (`decomp.PRIMARY_THRESH`, `decomp.SWEEP`) are **unchanged**.
All outputs are namespaced: `results/rerun/`, `results/tables/rerun_*.csv`, `results/figures/rerun_*.png`.

**Outcome of this notebook: STOP CONDITION A fired (2 A-genes survive; stop line = 10). Stages B–D were not run.**
See `FINDINGS_RERUN.md`."""),
    code(SETUP),
    md("""## Input (reused): pseudobulk and log-CPM
From a clean checkout, `rerun_common.ensure_logcpm()` (called inside `rerun_stageA.py`) rebuilds the pseudobulk from the cached OneK1K
h5ad via `src/run_pseudobulk_onek1k.py` and regenerates the log-CPM matrix with the same normalisation as the prior run. Filter log
(unchanged): `results/tables/onek1k_filter_log.csv`."""),
    code("""sh("src/run_pseudobulk_onek1k.py")
import pandas as pd
pd.read_csv(TAB / "onek1k_filter_log.csv")"""),
    md("""## A1–A4 and STOP A
`src/rerun_stageA.py` does, in order:

* **A0** design audit: donors nest in pools (asserted), `age ~ C(pool)` R², depth ~ pool R².
* **A1** per-gene commonality partition under the corrected model. Nuisance N = pool + within-type-centred log depth.
  `unique_age = R²(N+ct+age) − R²(N+ct)` is the **within-pool, depth-adjusted** age variance. The prior pooled model
  (`ct + age`) and a pool-only model (`ct + pool + age`) are computed alongside for comparison. Because the pseudobulks of one
  donor share its age, the nominal F test is anti-conservative; the primary inference is a **within-pool, donor-level permutation
  null** (1,000 draws, seed 20260914): empirical FDR per threshold, per-gene q-values, and a family-wise threshold on the maximum.
* **A2** A / I / MIXED sets with the unchanged primary thresholds and the unchanged sweep grid; comparison with the prior sets.
* **A3** known positives CDKN2A (p16, expected up) and LRRN3 (expected down): rank, slope sign, class, per-cell-type within-pool r.
* **A4** the prior run's confound filters (same lists, same fixed thresholds) applied to the new A set; survival at the primary
  thresholds and across the sweep.
* **STOP A** verdict → `results/rerun/stageA_STOP_verdict.txt`."""),
    code("""sh("src/rerun_stageA.py")"""),
    code("""from IPython.display import Image, display
for f in ("rerun_stageA_pooled_vs_corrected_unique_age.png", "rerun_stageA_joint_variance_distribution.png",
          "rerun_stageA_permutation_null.png", "rerun_stageA_Agene_survival.png"):
    display(Image(str(FIG / f)))"""),
    md("""## Verdict
STOP CONDITION A: **fewer than 10 A-genes survive** → the rerun halts here. Thresholds were not relaxed; Stages B, C, D were not run.
The cell below re-reads the verdict written by the script and asserts it, so a re-execution that changes the outcome fails loudly."""),
    code("""import json
s = json.load(open(RERUN_DIR / "stageA_summary.json"))
print(open(RERUN_DIR / "stageA_STOP_verdict.txt", encoding="utf-8").read())
assert s["fired"], "STOP A verdict changed on re-execution — FINDINGS_RERUN.md must be revisited"
print("\\nRerun halted at STOP A. Stages B-D not executed.")"""),
]
nbf.write(nbA, NB / "rerun_stageA_age_axis.ipynb")
print("wrote", NB / "rerun_stageA_age_axis.ipynb")

# ----------------------------------------------------------------------------- Stage B (run on the protocol owner's instruction after STOP A)
nbB = nbf.v4.new_notebook()
nbB.cells = [
    md("""# Phase 1 RERUN — Stage B: the compositional question

**Context.** Stage A halted at STOP A (2 A-genes). Stage B was then run **on the protocol owner's explicit instruction** as a measurement:
no stop condition, gene sets and thresholds untouched (`t_age` unchanged; the Stage A classes are read from
`results/tables/rerun_stageA_variance_decomposition.csv`). Write-up: `FINDINGS_STAGEB.md`.

* **B1** per-donor cell-type composition from **all** cells of the raw OneK1K h5ad (Azimuth L2 labels mapped as in the pseudobulk build;
  Doublet/Platelet/Eryth/HSPC excluded and logged; 16 analysis types + `other_PBMC`) → centred log-ratio → `age ~ C(pool) + CLR`.
  The number reported is the fraction of **within-pool** donor-age variance explained, with a within-pool donor-level permutation null.
* **B2** pool-grouped 5-fold CV, 10 repeats, inner pool-grouped alpha selection, HARD RULE 1 asserted in every fold: age from
  (a) composition alone, (b) within-cell-type expression, **all 13,904 genes, unrestricted** kernel ridge — per cell type and all 16
  types jointly, (c) a 2-coefficient recombination of the out-of-fold predictions of (a) and (b). Primary metric: within-pool R²
  (features pool-centred, target pool-centred, held-out predictions re-centred within each held-out pool). Secondary: naive R² on raw age.
  CIs: pool-cluster bootstrap (2,000 draws) plus the spread over the 10 repeats. One within-pool age shuffle per model as a pipeline check.
* **B3** share of the summed within-pool per-gene age variance (Stage A `unique_age`) in A vs MIXED vs I vs OTHER and by identity bin."""),
    code(SETUP),
    code("""sh("src/rerun_stageB.py")"""),
    code("""from IPython.display import Image, display
for f in ("rerun_stageB_composition_vs_age.png", "rerun_stageB_cv_composition_vs_expression.png",
          "rerun_stageB_oof_predictions.png", "rerun_stageB_age_signal_by_class.png"):
    display(Image(str(FIG / f)))"""),
    code("""import json
s = json.load(open(RERUN_DIR / "stageB_summary.json"))
print("B1  composition explains", f"{s['B1']['frac_within_pool_insample']:.1%}", "of within-pool age variance (in-sample; null", f"{s['B1']['null_increment_mean']/(1-s['B1']['r2_pool']):.1%})")
print("B2  pool-grouped CV within-pool R2: composition", f"{s['B2']['composition_r2_within']:.3f}", s['B2']['composition_ci'],
      "| expression all types", f"{s['B2']['expression_alltype_r2_within']:.3f}", s['B2']['expression_alltype_ci'],
      "| per-type median", f"{s['B2']['expression_per_type_r2_within_median']:.3f}", "| combined 2-df", f"{s['B2']['combined_2df_r2']:.3f}")
print("B3  share of summed within-pool age variance: A", f"{s['B3']['A_share_sum_unique_age']:.1%}", "| MIXED+I", f"{s['B3']['MIXED_plus_I_share']:.1%}",
      "| unique_ct>=0.5", f"{s['B3']['unique_ct_ge_0p5_share']:.1%}")"""),
]
nbf.write(nbB, NB / "rerun_stageB_composition.ipynb")
print("wrote", NB / "rerun_stageB_composition.ipynb")

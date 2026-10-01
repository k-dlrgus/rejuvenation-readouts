"""Prompt I — reproduce Lu et al. cell states, then score both instruments.

Namespaced under results/md2/. Does not modify any existing FINDINGS_*.md or FALSIFICATION.md.
Does not open GSE325735. Refit nothing. Frozen ruler used as-is.
Seed 20260914. Bootstrap seed 20260918. Spearman ρ is the reported statistic.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import RESULTS, DATA_RAW, DATA_PROC, ROOT, SEED  # noqa: E402
from brain_phase1_common import Logger, dump_json, load_json  # noqa: E402
from gtex_common import StopStep, spearman_safe, pearson_safe  # noqa: E402
from fibro_common import (  # noqa: E402
    FROZEN_RULER, FIBRO_SEED, BOOT_SEED, N_PERM, N_BOOT, N_RANDOM_DIR,
    jsonable, fmt_ci, bootstrap_rho_ci, PLURI_ENDOGENOUS, FIBRO_IDENTITY, OSKM_FAMILY,
)
from fibro3_common import MIN_CELLS, AGED_LINE, YOUNG_LINE, DAYS  # noqa: E402
from fibro2_common import load_frozen_ruler, rho_null_donor  # noqa: E402

MD2_DIR = RESULTS / "md2"
MD2_FIG = MD2_DIR / "figures"
MD2_PROC = DATA_PROC / "md2"
for _p in (MD2_DIR, MD2_FIG, MD2_PROC):
    _p.mkdir(parents=True, exist_ok=True)

MD2_SEED = int(SEED)  # 20260914
MD2_BOOT = int(BOOT_SEED)  # 20260918
PREREG_DATE = "2026-09-19"
MMC3_XLSX = RESULTS / "md" / "raw" / "suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx"
HALLMARK_LOCAL = RESULTS / "lowdim" / "genesets" / "MSigDB_Hallmark_2020.txt"
FIBRO3_DIR = RESULTS / "fibro3"
FIBRO2_DIR = RESULTS / "fibro2"
CXG_H5AD = DATA_RAW / "fibro2" / "cxg" / "8ae4471a-2e70-47b5-bc42-4f0df26a9996.h5ad"
CXG_DATASET_ID = "a19d1667-a7b5-4556-9e5f-f9bfa690c0f1"
GSE226189_TAR = DATA_RAW / "fibro2" / "GSE226189" / "GSE226189_RAW.tar"
GSE113957_FPKM = DATA_RAW / "fibro2" / "GSE113957" / "GSE113957_fpkm.txt.gz"

MANIFEST_PATH = MD2_DIR / "manifest.json"
PROGRESS_PATH = ROOT / "PROGRESS_MD2.md"
FINDINGS_PATH = ROOT / "FINDINGS_MD2.md"
PREREG_TASK1_FLAG = MD2_DIR / "PREREG_TASK1.flag"
PREREG_TASK2_FLAG = MD2_DIR / "PREREG_TASK2.flag"
PREREG_TASK3_FLAG = MD2_DIR / "PREREG_TASK3.flag"
CLUSTER_SETTINGS_FLAG = MD2_DIR / "CLUSTER_SETTINGS.flag"
DECLARED_BEFORE_SCORES_FLAG = MD2_DIR / "DECLARED_BEFORE_SCORES.flag"

KMEANS_KS = (3, 8, 15)
CLUSTER_NPC = 20
CLUSTER_N_INIT = 10
LOUVAIN_RES = 0.8
LOUVAIN_K_PARAM = 20
LOUVAIN_PRUNE = 1.0 / 15.0
LOUVAIN_NPC = 50
LOUVAIN_N_HVG = 3000
LOUVAIN_MIN_GENES = 500
LOUVAIN_MIN_CELLS_PER_GENE = 5
AMS_NBIN = 24
AMS_CTRL = 100
SEURAT_LOGNORMALIZE_SCALE = 10000
MD_TF_ADD = ("SNAI1", "ZEB1", "ZEB2", "TWIST1", "TWIST2")
HALLMARK_EMT_NAME = "Epithelial Mesenchymal Transition"
HALLMARK_TGFB_NAME = "TGF-beta Signaling"
HALLMARK_EMT_MSIGDB = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
HALLMARK_TGFB_MSIGDB = "HALLMARK_TGF_BETA_SIGNALING"

STAGE2_VERDICT_FIBRO = (
    "age score does not decline → the ruler reads a static donor property, not a "
    "modifiable state. Step 3 is not supported by this data."
)
FIBRO3_TASK3_AGED_D0D7 = (
    "Aged GM00731 d0→d7 (the all-cell T-B cell had p=+0.005): cluster 0 decline=-0.582 "
    "p=+0.652 n_ge=130/200 n_d0=4941 n_d7=2719; cluster 1 decline=+2.215 p=+0.144 "
    "n_ge=28/200 n_d0=76 n_d7=850; cluster 2 skipped (n_cells<20 at an endpoint)."
)
FIBRO3_TASK1_SENTENCE = (
    "Cell-intrinsic. The same cells' clusters individually reverse from d7 to d10. "
    "Then reprogramming genuinely pushes cells back along the age axis after d7, "
    "and that is a biological claim worth its own test. clusters=[0, 1, 2]."
)

CLUSTER_SETTINGS = """Written **before any age score**, 2026-09-19. All four resolutions are pre-registered. Not chosen after seeing which one gives a cleaner age trajectory.

Scope: aged donor GM00731, all timepoints concatenated, jointly clustered. Task 2 Louvain is repeated independently on GM23815 (never pooled across donors).

k = 3, 8, 15 (FIBRO3 setting, only k changes):
- Genes: frozen-ruler overlap genes with count-sum > 0 (FINDINGS_FIBRO3.md Task 1 / FINDINGS_FIBRO.md Stage 2 cluster spec).
- Transform: log1p(CP10k) on those genes. Library size = row sum of the kept overlap genes.
- Dimensionality: sklearn.decomposition.PCA(n_components=20, or n_cells−1 if smaller, svd_solver="randomized", random_state=20260914) on float32 log1p(CP10k). Cell order: libraries concatenated in day order 0, 3, 7, 10.
- Cluster: sklearn.cluster.KMeans(n_clusters=k, random_state=20260914, n_init=10).
- k=3 reuses the frozen labels in results/fibro3/t1_cell_labels.csv for GM00731 (the FIBRO3 setting). k=8 and k=15 are fit on the same PCA recipe, new KMeans.
- No n_genes>500 filter (FIBRO3 did not apply it). MIN_CELLS=20 to score a cluster×timepoint.

Louvain (their pipeline, STAR Methods of Lu et al. Cell 2025 188:5895–5911, “Single-cell RNA-seq analysis of newly generated data”):
- Cells with more than 500 genes; genes detected in at least 5 cells.
- STAR Methods specify Seurat v5.0.0, sctransform v2, regression on percent mitochondrial, Louvain, AddModuleScore. This environment has no R/Seurat/sctransform. Deviation, not a silent substitute: LogNormalize (Seurat NormalizeData default: ln(1 + 10000 * count / tot); paper public-data prose also describes multiplying by the median total UMIs — we use scale.factor=10000 and report the prose/default discrepancy), then ScaleData-style residualization of percent.mt (symbol upper startswith MT-) on the top 3000 highly variable genes (Seurat SCTransform default variable.features.n; STAR Methods do not state n HVGs), PCA 50 PCs (Seurat RunPCA default; STAR Methods do not state n PCs), shared-nearest-neighbour graph (k.param=20, prune.SNN=1/15, Seurat FindNeighbors defaults; STAR Methods do not state these), Louvain via networkx.community.louvain_communities, resolution=0.8 (Seurat FindClusters default; STAR Methods do not state the resolution). Seed 20260914.
- n clusters obtained is reported, not retuned.
- Slingshot from day-0 fibroblast clusters is optional and report-only; not required for scoring.

MIN_CELLS=20. Seed 20260914. Frozen ruler used as-is. Refit nothing.
"""

PREREG_TASK1 = """Written **before any age score**, 2026-09-19. No threshold in this block is re-tuned after numbers exist.

Cluster the aged donor GM00731 jointly across timepoints at k = 3 (the FIBRO3 setting), 8, 15, and the Louvain pipeline in CLUSTER_SETTINGS. Fix all settings before any age score is computed.

For each resolution, per cluster × timepoint (≥20 cells): n cells, fraction of timepoint, frozen age score, MD score, TGF-β score, pluripotency−fibroblast score (frozen lists from FINDINGS_FIBRO.md), and the pre-registered 200-direction random null on d0→d7 and d0→d10 separately.

Frozen ruler: results/fibro/frozen_ruler_ridge_raw.npz used as-is. Age null = 200 permuted-ruler-weight directions (seed 20260914), empirical p = (n_random ≥ real + 1) / (n_random + 1). Endpoint decline = age_score(day 0) − age_score(day T). Positive = younger on the ruler. MD/TGF-β are Seurat AddModuleScore on LogNormalize counts (control-gene bins matched on average expression; nbin=24, ctrl=100). Do not average across donors, clusters, or resolutions. d7 is not the Stage 2 endpoint.

**Pre-registered reading (only the outcome that fired):**
- At finer resolution some cluster shows a d0→d7 or d0→d10 age decline with p ≤ 0.05 → FINDINGS_FIBRO3.md's per-cluster negative was a resolution artifact. Report which clusters, at which resolutions, and how stable the result is across resolutions. State plainly that this supersedes the FIBRO3 per-cluster reading, quoting the sentence it supersedes.
- No cluster at any resolution shows p ≤ 0.05 → the per-cluster negative holds and is not a resolution artifact. The frozen ruler does not move in these cells at any granularity tested.
- Results appear at k=15 but not at Louvain, or vice versa → report per resolution, do not average, and say the finding is resolution-dependent.
"""

PREREG_TASK2 = """Written **before any MD vs frozen-ruler comparison**, 2026-09-19. No threshold in this block is re-tuned after numbers exist.

On the Louvain clustering (their pipeline), per cluster × timepoint, both donors reported separately:

1. Reproduce their Figure 3G result. Label each Louvain cluster by the mmc3.xlsx sheet Reprog_cell_state_signatures column with the highest mean AddModuleScore (Fibroblast, PartialReprog, EarlyPluripotency, Pluripotency, NonReprog). Reproduction matches their reported direction iff (a) at least one PartialReprog-assigned cluster and one NonReprog-assigned cluster exist with ≥20 cells in the aged donor, and (b) mean per-cell MD score in PartialReprog-assigned clusters is lower than in NonReprog-assigned clusters (aged donor). If it does not match, the comparison below is not valid — say so and stop Task 2.
2. Score both instruments on the same cells: frozen GTEx fibroblast ruler and MD score (and their aging signature if buildable).
3. Cross-correlate. Across cluster × timepoint rows with ≥20 cells, Spearman ρ between frozen age score and MD score, with permutation null (shuffle MD, hold age; n_perm=200, seed 20260914) and bootstrap CI (B=200, seed 20260918). Donors reported separately; not averaged.
4. Composition test on their metric. Per cluster, d0→d7 and d0→d10 MD decline with a 200-draw size-matched random gene-set null (same n genes, matched on mean-expression bin; AddModuleScore; seed 20260914). Not permuted ruler weights. State the difference.

**Pre-registered reading (only the outcome that fired):**
- MD moves per cluster and the frozen ruler does too (Task 1 outcome 1) → both instruments agree; their bend claim is supported and our earlier negative was our clustering. Report it that way, without hedging.
- MD moves per cluster, the frozen ruler does not → the instruments disagree on the same cells. Describe both, name what would settle it, and declare no winner.
- Neither moves per cluster → our reproduction does not recover their published result; report the reproduction failure as the finding and do not present it as a refutation.
- MD moves but fails its size-matched random-gene-set null → their metric's movement is not specific to mesenchymal genes in this dataset. Report the null value; this is a strong claim and needs the number stated next to it.

“MD moves per cluster” = at least one aged-donor Louvain cluster with ≥20 cells at both endpoints shows MD decline d0→d7 or d0→d10 with p ≤ 0.05 on the size-matched gene-set null.
"""

PREREG_TASK3 = """Written **before any external-cohort MD ρ**, 2026-09-19. Report-only. No gate.

Their aging signature comes from Fleischer et al. 2018 (~100 dermal fibroblasts); our ruler comes from GTEx cultured fibroblasts (n=652). If the Fleischer DE sets are buildable from GSE113957 exactly as published (DESeq2 old-vs-young, adjusted p < 0.05 and log2FC > 0.5, old/young from PCA separation), report each instrument's Spearman ρ against donor age in the external cohorts already validated in FINDINGS_FIBRO2.md (CELLxGENE a19d1667-a7b5-4556-9e5f-f9bfa690c0f1, n=93, 10x; GSE226189, n=82, bulk), with nulls and CIs. One table. No gate.

If GSE113957's format prevents exact reconstruction, say so and report the aging-signature comparison as not run rather than approximating. MD and the frozen ruler are still scored. Frozen ruler numbers are read from results/fibro2/; the ruler is not refit.
"""

STAR_METHODS_QUOTES = {
    "gene_sets": (
        "We obtained the gene lists from the hallmark gene sets: HALLMARK EPITHELIAL "
        "MESENCHYMAL TRANSITION and HALLMARK TGF BETA SIGNALING (https://www.gsea-msigdb.org/). "
        "We added SNAI1, ZEB1, ZEB2, TWIST1, and TWIST2 to HALLMARK EPITHELIAL MESENCHYMAL "
        "TRANSITION to compose the MD signature based on literature review (Table S3). "
        "We derived human fibroblast aging signatures using RNA-seq data from Fleischer et al. "
        "(2018) by determining significant differentially expressed genes (adjusted p-value < 0.05, "
        "log2(fold change) > 0.5) in old cells versus young cells as classified based on PCA "
        "separation using DESeq2 v1.40.2 (Table S3)."
    ),
    "new_scrna": (
        "Count matrices were loaded to create Seurat objects using Seurat v5.0.0. Cells with more "
        "than 500 genes and genes detected in at least 5 cells were kept for further processing. "
        "Separate samples were merged and normalized using sctransform v2, with regression for "
        "percent mitochondrial genes. Gene lists were used to calculate pathway score using the "
        "Seurat AddModuleScore function. Trajectory inference was performed using Slingshot v2.12.0, "
        "designating day 0 fibroblast populations as starting clusters."
    ),
    "public_scrna_score": (
        "All the module expressions of the MD gene were based on log-scaled normalized UMI counts "
        "and calculated by the function AddModuleScore in Seurat package."
    ),
    "lognormalize_prose": (
        "To normalize sequencing reads for each gene, we utilized the \"NormalizeData\" function "
        "with the parameter \"normalization.method = \"LogNormalize\". This involved dividing the "
        "UMI counts of each gene by the total UMIs of the cell, multiplying by the median of the "
        "total UMIs, and transforming the result using the natural logarithm."
    ),
    "fig3g": (
        "The partially reprogrammed populations showed a robust reversal of transcriptomic aging "
        "changes, while non-reprogrammed populations maintained the aging transcriptome (Figure 3G). "
        "This pattern was accompanied by downregulation of the MD and TGF-β pathway scores in the "
        "partially reprogramming cells and their maintained expression in non-reprogrammed cells "
        "(Figure 3G)."
    ),
}


def md2_log_banner(log, stage: str):
    log("=" * 100)
    log(f"MD2 {stage}  seed={MD2_SEED}  boot_seed={MD2_BOOT}  "
        f"n_perm={N_PERM}  n_boot={N_BOOT}  n_random={N_RANDOM_DIR}")
    log("=" * 100)
    log("FALSIFICATION.md and existing FINDINGS*.md are not modified.")
    log(f"Frozen ruler: {FROZEN_RULER} exists={FROZEN_RULER.exists()}")
    log("Primary metric: Spearman ρ (n_perm=200); uncalibrated R² is never a gate.")
    log("Refit nothing. Do not open GSE325735. Do not promote d7 to the Stage 2 endpoint.")
    log("Fail loudly: no substitute columns, no inferred ages, no silent coerce.")


def load_manifest():
    if MANIFEST_PATH.exists():
        return load_json(MANIFEST_PATH)
    return dict(
        seed=MD2_SEED, boot_seed=MD2_BOOT, files={}, failures=[],
        columns={}, unusable=[], unverified=[], gene_sets_read=[],
        status="INIT",
    )


def save_manifest(man):
    dump_json(MANIFEST_PATH, jsonable(man))
    return MANIFEST_PATH


def _clip_msg(message, n=400):
    s = str(message)
    return s if len(s) <= n else s[:n] + "…"


def record_failure(step, message, details=None):
    man = load_manifest()
    rec = dict(step=step, message=_clip_msg(message), details=jsonable(details or {}))
    fails = man.setdefault("failures", [])
    if not any(f.get("step") == step and f.get("message") == rec["message"] for f in fails):
        fails.append(rec)
    man["status"] = "STOP"
    save_manifest(man)
    return rec


def log_columns(man_key, columns, source):
    man = load_manifest()
    man.setdefault("columns", {})[man_key] = dict(
        source=str(source), columns=[str(c) for c in list(columns)],
    )
    save_manifest(man)


def log_geneset(man_key, name, genes, source):
    man = load_manifest()
    man.setdefault("gene_sets_read", [])
    rec = dict(
        key=man_key, name=str(name), n=int(len(genes)),
        source=str(source), genes=[str(g) for g in list(genes)],
    )
    existing = man["gene_sets_read"]
    if not any(x.get("key") == man_key and x.get("name") == rec["name"] for x in existing):
        existing.append(rec)
    save_manifest(man)


def write_prereg_flag(path: Path, text: str):
    if path.exists():
        existing = path.read_text(encoding="utf-8")
        if existing.strip() != text.strip():
            raise StopStep(
                "prereg",
                f"{path.name} already exists and does not match the frozen block. Not rewriting.",
            )
        return path
    path.write_text(text, encoding="utf-8")
    return path


def write_all_prereg_flags():
    write_prereg_flag(PREREG_TASK1_FLAG, PREREG_TASK1)
    write_prereg_flag(PREREG_TASK2_FLAG, PREREG_TASK2)
    write_prereg_flag(PREREG_TASK3_FLAG, PREREG_TASK3)
    write_prereg_flag(CLUSTER_SETTINGS_FLAG, CLUSTER_SETTINGS)
    return PREREG_TASK1_FLAG, PREREG_TASK2_FLAG, PREREG_TASK3_FLAG, CLUSTER_SETTINGS_FLAG


def _fmt_num(x, d=3):
    if x is None:
        return "NA"
    if isinstance(x, (bool, np.bool_)):
        return str(bool(x))
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    if not np.isfinite(v):
        return "NA"
    return f"{v:+.{d}f}"


def _read_csv_if(path: Path):
    if path.exists():
        return pd.read_csv(path)
    return None


def _disk_row_lines(title, path, max_cols=16, max_rows=80):
    df = _read_csv_if(path)
    out = [f"### {title} (`{path.name}`)"]
    if df is None or not len(df):
        out.append("- missing or empty")
        out.append("")
        return out
    for i, (_, r) in enumerate(df.iterrows()):
        if i >= max_rows:
            out.append(f"- … {len(df) - max_rows} more rows")
            break
        bits = [f"{c}={r[c]}" for c in list(df.columns)[:max_cols]]
        out.append("- " + " ".join(bits))
    out.append("")
    return out


def progress_snapshot(next_action, stop=None, extra=""):
    man = load_manifest() if MANIFEST_PATH.exists() else {}
    stop_s = stop if stop else man.get("status", "running")
    lines = [
        "# PROGRESS_MD2",
        "",
        f"**STOP status:** {stop_s}",
        f"**Next action:** {next_action}",
        "",
        "## Seeds / gates",
        "",
        f"- seed `{MD2_SEED}`  boot `{MD2_BOOT}`  n_perm={N_PERM}  n_boot={N_BOOT}  n_random={N_RANDOM_DIR}",
        "- primary metric: Spearman ρ; uncalibrated R² is never a gate",
        f"- frozen ruler: `{FROZEN_RULER}` exists={FROZEN_RULER.exists()}",
        f"- MIN_CELLS={MIN_CELLS}  k_means k={list(KMEANS_KS)}  Louvain resolution={LOUVAIN_RES} (Seurat default; not stated in STAR Methods)",
        f"- PREREG_TASK1.flag exists={PREREG_TASK1_FLAG.exists()}",
        f"- PREREG_TASK2.flag exists={PREREG_TASK2_FLAG.exists()}",
        f"- PREREG_TASK3.flag exists={PREREG_TASK3_FLAG.exists()}",
        f"- CLUSTER_SETTINGS.flag exists={CLUSTER_SETTINGS_FLAG.exists()}",
        f"- DECLARED_BEFORE_SCORES.flag exists={DECLARED_BEFORE_SCORES_FLAG.exists()}",
        "- GSE325735: out of scope",
        "- d7 is not the Stage 2 endpoint",
        "- Refit nothing",
        "",
        "## Frozen ruler path",
        "",
        f"- `{FROZEN_RULER}` exists={FROZEN_RULER.exists()}",
        "",
        "## Clustering settings (verbatim, frozen before any age score)",
        "",
        CLUSTER_SETTINGS,
        "",
        "## Pre-registered blocks (verbatim)",
        "",
        "### Task 1",
        "",
        PREREG_TASK1,
        "",
        "### Task 2",
        "",
        PREREG_TASK2,
        "",
        "### Task 3",
        "",
        PREREG_TASK3,
        "",
        "## Gene-set versions and counts (from disk; never from memory)",
        "",
    ]
    gs_path = MD2_DIR / "genesets_summary.json"
    if gs_path.exists():
        rec = load_json(gs_path)
        lines.append(f"- {json.dumps(jsonable(rec), ensure_ascii=False)[:8000]}")
        lines.append("")
    else:
        lines.append("- genesets_summary.json not yet written")
        lines.append("")
    lines += ["## Finished cells (from disk)", ""]
    for name, path in (
        ("Task 1 cluster table", MD2_DIR / "t1_cluster_table.csv"),
        ("Task 1 age null", MD2_DIR / "t1_age_null.csv"),
        ("Task 1 reading", MD2_DIR / "t1_reading.json"),
        ("Task 2 Louvain table", MD2_DIR / "t2_cluster_table.csv"),
        ("Task 2 reproduction", MD2_DIR / "t2_reproduction.json"),
        ("Task 2 cross-correlation", MD2_DIR / "t2_crosscorr.csv"),
        ("Task 2 MD composition null", MD2_DIR / "t2_md_null.csv"),
        ("Task 2 reading", MD2_DIR / "t2_reading.json"),
        ("Task 3 instruments", MD2_DIR / "t3_instruments.csv"),
    ):
        if str(path).endswith(".json"):
            if path.exists():
                rec = load_json(path)
                lines.append(f"### {name} (`{path.name}`)")
                lines.append(f"- {json.dumps(jsonable(rec), ensure_ascii=False)[:4000]}")
                lines.append("")
            else:
                lines.append(f"### {name} (`{path.name}`)")
                lines.append("- missing")
                lines.append("")
        else:
            lines += _disk_row_lines(name, path)
    for js_name in (
        "genesets_summary.json", "t1_summary.json", "t2_summary.json", "t3_summary.json",
        "cluster_meta.json",
    ):
        p = MD2_DIR / js_name
        if p.exists():
            rec = load_json(p)
            lines.append(f"### {js_name}")
            lines.append(f"- {json.dumps(jsonable(rec), ensure_ascii=False)[:4000]}")
            lines.append("")
    cols = man.get("columns") or {}
    if cols:
        lines += ["## Columns actually read", ""]
        for k, v in cols.items():
            lines.append(f"- **{k}** source=`{v.get('source')}` columns={v.get('columns')}")
        lines.append("")
    gsets = man.get("gene_sets_read") or []
    if gsets:
        lines += ["## Gene-set names actually read", ""]
        for g in gsets:
            lines.append(
                f"- **{g.get('key')}** name=`{g.get('name')}` n={g.get('n')} source=`{g.get('source')}`"
            )
        lines.append("")
    if extra:
        lines += ["## Note", "", extra, ""]
    if man.get("failures"):
        lines += ["## Failures (manifest)", ""]
        for f in man["failures"]:
            lines.append(f"- **{f.get('step')}:** {_clip_msg(f.get('message'), 240)}")
        lines.append("")
    lines += [
        "## Files",
        "",
        "- `src/md2_common.py`, `md2_genesets.py`, `md2_score.py`, `md2_cluster.py`, "
        "`md2_task1.py`, `md2_task2.py`, `md2_task3.py`, `md2_findings.py`, `md2_run.py`",
        "- `results/md2/`",
        "- `FINDINGS_MD2.md`",
        "- `PROGRESS_MD2.md`",
        "",
    ]
    PROGRESS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return PROGRESS_PATH

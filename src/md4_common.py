"""Prompt K — compare cell states at a shared timepoint (Figure 3G contrast).

Namespaced under results/md4/. Does not modify any existing FINDINGS_*.md or FALSIFICATION.md.
Does not open GSE325735. Refit nothing. Frozen ruler used as-is.
Seed 20260914. Bootstrap seed 20260918. Spearman ρ is the reported statistic.
Reuse src/md2_*.py Louvain labels + AddModuleScore exactly; do not re-cluster, re-label, or re-tune.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import RESULTS, DATA_PROC, ROOT, SEED  # noqa: E402
from brain_phase1_common import Logger, dump_json, load_json  # noqa: E402
from gtex_common import StopStep, spearman_safe  # noqa: E402
from fibro_common import (  # noqa: E402
    FROZEN_RULER, BOOT_SEED, N_PERM, N_BOOT, N_RANDOM_DIR, jsonable,
    PLURI_ENDOGENOUS, FIBRO_IDENTITY, OSKM_FAMILY,
)
from fibro3_common import AGED_LINE, YOUNG_LINE, DAYS  # noqa: E402
from fibro2_common import load_frozen_ruler, PCA_MAHAL_K  # noqa: E402
from md2_common import (  # noqa: E402
    MD2_DIR, MD2_PROC, MMC3_XLSX, AMS_NBIN, AMS_CTRL, SEURAT_LOGNORMALIZE_SCALE,
    LOUVAIN_RES, CLUSTER_SETTINGS,
)
from md3_common import MD3_DIR  # noqa: E402

MD4_DIR = RESULTS / "md4"
MD4_FIG = MD4_DIR / "figures"
MD4_PROC = DATA_PROC / "md4"
for _p in (MD4_DIR, MD4_FIG, MD4_PROC):
    _p.mkdir(parents=True, exist_ok=True)

MD4_SEED = int(SEED)  # 20260914
MD4_BOOT = int(BOOT_SEED)  # 20260918
PREREG_DATE = "2026-09-19"
FIBRO2_DIR = RESULTS / "fibro2"
FIBRO3_DIR = RESULTS / "fibro3"

MANIFEST_PATH = MD4_DIR / "manifest.json"
PROGRESS_PATH = ROOT / "PROGRESS_MD4.md"
FINDINGS_PATH = ROOT / "FINDINGS_MD4.md"
PREREG_TASK1_FLAG = MD4_DIR / "PREREG_TASK1.flag"
PREREG_TASK2_FLAG = MD4_DIR / "PREREG_TASK2.flag"
PREREG_TASK3_FLAG = MD4_DIR / "PREREG_TASK3.flag"
DECLARED_BEFORE_SCORES_FLAG = MD4_DIR / "DECLARED_BEFORE_SCORES.flag"

MIN_CELLS_CONTRAST = 50
STATE_NAMES = ("Fibroblast", "PartialReprog", "EarlyPluripotency", "Pluripotency", "NonReprog")
PRIMARY_PAIR = ("PartialReprog", "NonReprog")
CONTEXT_PAIRS = (
    ("Fibroblast", "NonReprog"),
    ("PartialReprog", "Fibroblast"),
)
ALL_PAIRS = (PRIMARY_PAIR,) + CONTEXT_PAIRS
INSTRUMENTS = (
    "frozen_ruler",
    "md_score",
    "tgfb_score",
    "age_up_minus_age_down",
    "pluri_primary",
)
# Permutation seed offsets: MD4_SEED + 1000*donor_i + 10*day + pair_i
# Bootstrap seed offsets: MD4_BOOT + same.
PAIR_INDEX = {PRIMARY_PAIR: 0, CONTEXT_PAIRS[0]: 1, CONTEXT_PAIRS[1]: 2}
DONOR_INDEX = {AGED_LINE: 0, YOUNG_LINE: 1}

MD2_POOLED_MD_PR = 0.15047059685308203
MD2_POOLED_MD_NR = 0.4394484253245684
MD2_N_PR = 5832
MD2_N_NR = 5459

MD3_TASK1_STATUS = "instruments_disagree_on_claim_population"
MD3_TASK1_READING = (
    "Ruler does not decline in the pooled PartialReprog state while MD does, both with "
    "their own nulls. The instruments disagree on exactly the population the published "
    "claim concerns. No winner is declared."
)
MD3_TASK1_NUMBERS = (
    "Pooled PartialReprog ruler: n_ok=2 n_pass=0. Pooled PartialReprog MD: n_ok=2 n_pass=1 "
    "hits d0→d7 decline=+0.007 p=+0.015 n_cells_d0=3 n_cells_end=578. "
    "Frozen ruler d0→d7 decline=-3.339 p=+0.403."
)
MD3_TASK2_READING = (
    "Δρ CI excludes zero in favour of the ruler in 2 of 3 named cohorts. MD and/or their "
    "aging signature predict donor age materially worse than the ruler on independent "
    "fibroblast cohorts."
)
MD2_NEGATIVE_HOLDS = (
    "No cluster at any resolution shows p ≤ 0.05. The per-cluster negative holds and is "
    "not a resolution artifact. The frozen ruler does not move in these cells at any "
    "granularity tested."
)
MD2_INSTRUMENTS_DISAGREE = (
    "MD moves per cluster, the frozen ruler does not. The instruments disagree on the "
    "same cells. What would settle it: an independent fibroblast age instrument scored "
    "on these same Louvain clusters, or the authors' Seurat/sctransform object with the "
    "frozen ruler projected onto it. No winner is declared."
)
STAGE2_VERDICT_FIBRO = (
    "age score does not decline → the ruler reads a static donor property, not a "
    "modifiable state. Step 3 is not supported by this data."
)
FIBRO3_EXTRAP_RHO = (
    "nn_euclidean: age score does not rise with distance ρ=-0.625 CI=[-0.834, -0.320] n=22. "
    "Not a gate."
)
WHAT_WOULD_SETTLE = (
    "An independent fibroblast age instrument scored on these same Louvain-labelled cells "
    "at the same donor × timepoint, or the authors' Seurat v5 / sctransform v2 object with "
    "the frozen ruler projected onto the same PartialReprog vs NonReprog split."
)

PREREG_TASK1 = """Written **before any MD4 cell-level instrument score**, 2026-09-19. No threshold in this block is re-tuned after numbers exist.

Using the src/md2_*.py Louvain clustering already computed (not re-run, not re-tuned). Cell-state labels are the md2 assignment: argmax of mean AddModuleScore of mmc3 Reprog_cell_state_signatures columns Fibroblast / PartialReprog / EarlyPluripotency / Pluripotency / NonReprog. Labels are read from results/md2/t2_cluster_labels_*.csv. Not re-labelled. AddModuleScore for MD and TGF-β is the per-cell vector already stored in data/processed/md2/louvain_ams_*.npz (the md2 computation). Table S3 Age up / Age down are scored with the same AddModuleScore implementation (nbin=24, ctrl=100) on the same LogNormalize Louvain matrix, using the bins stored in that npz. Frozen ruler: results/fibro/frozen_ruler_ridge_raw.npz used as-is. Refit nothing.

Per donor, per timepoint, never pooled across donors or timepoints. A timepoint qualifies for a two-state contrast iff both states have ≥50 Louvain-kept cells at that donor × day. Cell counts and which timepoints qualify are reported before contrasts.

Per cell:
- frozen ruler score = edgeR log2-CPM (prior.count=2) on that cell's ruler-aligned UMI counts with TMM size factor fixed at 1 (library-size only). TMM is not re-estimated among cells, among states, or among timepoints: a TMM among ~20k cells would be a new size-factor fit on GSE297234 and is not the frozen bulk spec. Frozen μ, σ, w applied as-is; overlap genes that are all-zero in that donor's Louvain ruler matrix are left at z=0.
- MD score, from md2 louvain_ams (not recomputed).
- TGF-β score, from md2 louvain_ams (not recomputed).
- Table S3 age-up minus age-down = AMS(Age up) − AMS(Age down).
- pluripotency−fibroblast = mean-z(PLURI_ENDOGENOUS) − mean-z(FIBRO_IDENTITY) on the same per-cell Z as the ruler, drop_oskm=False (md2 pluri_primary).

Primary contrast: PartialReprog vs NonReprog. Context rows (same procedure, not pooled into the primary): Fibroblast vs NonReprog; PartialReprog vs Fibroblast.

Per qualifying timepoint, per instrument: Δmean = mean(state_A) − mean(state_B); Δmedian = median(state_A) − median(state_B). Cell-level permutation null: 200 shuffles of the two-state labels among the cells in those two states at that donor × timepoint (group sizes held; timepoint and donor held, so neither can drive the null). Empirical p_lower = (n_null ≤ Δobs + 1) / (n_null + 1); p_two_sided = (n_|null| ≥ |Δobs| + 1) / (n_null + 1). Donor-cell bootstrap 95% percentile CI on Δmean and Δmedian: resample cells within each state independently, B=200, seed 20260918. Effect size: Cohen's d on the mean difference, pooled sd with n−1, same sign as Δmean. Not rank-biserial. Seed 20260914 for permutations with a documented per-(donor, day, contrast) offset.

Sanity check, before any new score: reproduce FINDINGS_MD2.md pooled MD means on the aged donor (PartialReprog +0.150, NonReprog +0.439; t2_reproduction.json). If they do not reproduce, stop. Cell-state labels have drifted.

**Pre-registered reading (only the outcome that fired; each donor separately, never pooled; each qualifying timepoint reported, never pooled):**
- MD lower in PartialReprog than NonReprog and the frozen ruler lower too, both beating their nulls → both instruments agree that partially reprogrammed cells read younger. Their claim is supported and our earlier negatives reflected the wrong test (a within-state trajectory with no baseline). Say so plainly and quote the sentences superseded.
- MD lower, ruler not lower (or higher) → the instruments disagree on the population the published claim concerns, in a well-powered comparison. Report both numbers side by side, declare no winner, and name what would settle it.
- Neither differs → our reproduction does not recover their Figure 3G contrast; report the reproduction failure as the finding, not as a refutation.
- Ruler lower, MD not → report as an unexplained discrepancy.

“Lower” / “beating their nulls” = Δmean = mean(PartialReprog) − mean(NonReprog) is < 0 and p_lower ≤ 0.05 on that instrument's cell-level permutation null, on that donor × timepoint. GM23815 is tabulated, never pooled with GM00731, never averaged.

If two qualifying timepoints in one donor fire different bullets, report both, do not pick, do not average.
"""

PREREG_TASK2 = """Written **before any MD4 extrapolation distance or ρ**, 2026-09-19. Report-only. No gate.

For each donor × state × timepoint with n_cells ≥ 1 (n=0 is not a row): sum Louvain-kept UMIs to a pseudobulk on the frozen-ruler gene space (rulerY aligned as in md2/md3). TMM among that donor's state×timepoint rows, then frozen μ/σ/w. Distance to the GTEx fibroblast training z-space: PCA k=50 nearest-neighbour Euclidean and Mahalanobis, same functions as src/fibro2_extrap.py / FINDINGS_FIBRO2.md. Spearman ρ across those rows between each distance and the frozen ruler score, permutation null n_perm=200 seed 20260914 (shuffle the ruler score, hold distance), bootstrap CI B=200 seed 20260918. Donors reported separately and as concatenated rows; never averaged.

FINDINGS_FIBRO3.md Task 2 (the same PCA k=50 NN / Mahalanobis measures as FINDINGS_FIBRO2.md) found ruler score vs nn_euclidean ρ=-0.625. If PartialReprog sits further from training than NonReprog, that correlation would push the Task 1 frozen-ruler contrast toward “PR reads younger”. Reported as direction only. Not a gate. Not used to explain a Task 1 result away.
"""

PREREG_TASK3 = """Dated 2026-09-19. FINDINGS_MD3.md Task 1 status `instruments_disagree_on_claim_population` rested on n₀=3 PartialReprog cells at day 0 and is withdrawn as unsupported, superseded by MD4 Task 1. FINDINGS_MD3.md is not edited. Task 2 of MD3 (the Δρ result) is unaffected and stands.
"""

FIG3G_QUOTE = (
    "The partially reprogrammed populations showed a robust reversal of transcriptomic aging "
    "changes, while non-reprogrammed populations maintained the aging transcriptome (Figure 3G). "
    "This pattern was accompanied by downregulation of the MD and TGF-β pathway scores in the "
    "partially reprogramming cells and their maintained expression in non-reprogrammed cells "
    "(Figure 3G)."
)


def md4_log_banner(log, stage: str):
    log("=" * 100)
    log(f"MD4 {stage}  seed={MD4_SEED}  boot_seed={MD4_BOOT}  "
        f"n_perm={N_PERM}  n_boot={N_BOOT}  n_random={N_RANDOM_DIR}")
    log("=" * 100)
    log("FALSIFICATION.md and existing FINDINGS*.md are not modified.")
    log(f"Frozen ruler: {FROZEN_RULER} exists={FROZEN_RULER.exists()}")
    log("Primary metric: Spearman ρ (n_perm=200); uncalibrated R² is never a gate.")
    log("Refit nothing. Do not open GSE325735. Do not promote d7 to the Stage 2 endpoint.")
    log("Fail loudly: no substitute columns, no inferred ages, no silent coerce.")
    log("Louvain / labels / AddModuleScore inherited from md2; not re-tuned.")
    log(f"MIN_CELLS_CONTRAST={MIN_CELLS_CONTRAST}. Donors never pooled. Timepoints never pooled.")


def load_manifest():
    if MANIFEST_PATH.exists():
        return load_json(MANIFEST_PATH)
    return dict(
        seed=MD4_SEED, boot_seed=MD4_BOOT, files={}, failures=[],
        columns={}, unusable=[], unverified=[], gene_sets_read=[],
        id_types={}, status="INIT",
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
    save_manifest(man)
    return rec


def log_columns(man_key, columns, source):
    man = load_manifest()
    man.setdefault("columns", {})[man_key] = dict(
        source=str(source), columns=[str(c) for c in list(columns)],
    )
    save_manifest(man)


def log_geneset(man_key, name, genes, source, id_type=None):
    man = load_manifest()
    man.setdefault("gene_sets_read", [])
    rec = dict(
        key=man_key, name=str(name), n=int(len(genes)),
        source=str(source), genes=[str(g) for g in list(genes)],
        id_type=id_type,
    )
    existing = man["gene_sets_read"]
    if not any(x.get("key") == man_key and x.get("name") == rec["name"] for x in existing):
        existing.append(rec)
    save_manifest(man)


def log_id_type(man_key, rec):
    man = load_manifest()
    man.setdefault("id_types", {})[man_key] = jsonable(rec)
    save_manifest(man)


_ENS_RE = re.compile(r"^ENS[A-Z]*[GT]\d+", re.I)
_NM_RE = re.compile(r"^(NM_|NR_|XM_|XR_)\d+", re.I)
_SYM_RE = re.compile(r"^[A-Z][A-Z0-9\-.]{1,18}$")


def detect_id_type(ids, log, tag, source):
    """Copy of md3_idtype.detect_id_type; writes to the md4 manifest, not md3."""
    ids = [str(x) for x in np.asarray(ids).astype(str)]
    n = int(len(ids))
    if n == 0:
        rec = dict(tag=tag, source=str(source), n=0, kind="empty",
                   n_ensembl=0, n_refseq=0, n_symbol_like=0, examples=[])
        log_id_type(tag, rec)
        if log is not None:
            log(f"[id_type {tag}] EMPTY source={source}")
        return rec
    n_ens = int(sum(1 for x in ids if _ENS_RE.match(x.split(".")[0])))
    n_ref = int(sum(1 for x in ids if _NM_RE.match(x)))
    n_sym = int(sum(1 for x in ids if _SYM_RE.match(x.upper()) and not _ENS_RE.match(x) and not _NM_RE.match(x)))
    frac_ens = n_ens / n
    frac_ref = n_ref / n
    frac_sym = n_sym / n
    if frac_ens >= 0.8:
        kind = "ensembl"
    elif frac_ref >= 0.5:
        kind = "refseq_transcript"
    elif frac_sym >= 0.5:
        kind = "symbol"
    else:
        kind = "unknown"
    rec = dict(
        tag=tag, source=str(source), n=n, kind=kind,
        n_ensembl=n_ens, n_refseq=n_ref, n_symbol_like=n_sym,
        frac_ensembl=frac_ens, frac_refseq=frac_ref, frac_symbol_like=frac_sym,
        examples=ids[:12],
    )
    log_id_type(tag, rec)
    if log is not None:
        log(f"[id_type {tag}] kind={kind} n={n} n_ensembl={n_ens} n_refseq={n_ref} "
            f"n_symbol_like={n_sym} source={source} examples={ids[:8]}")
    return rec


def require_mappable(rec, allowed, log, step):
    if rec.get("kind") not in allowed:
        msg = (
            f"{rec.get('tag')}: ID type {rec.get('kind')!r} not in {allowed}. "
            f"n={rec.get('n')} n_ensembl={rec.get('n_ensembl')} n_refseq={rec.get('n_refseq')} "
            f"n_symbol_like={rec.get('n_symbol_like')} examples={rec.get('examples')}. "
            "Not mapping. This step stops."
        )
        if log is not None:
            log(f"[id_type STOP] {msg}")
        raise StopStep(step, msg, details=rec)
    return rec


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
    return PREREG_TASK1_FLAG, PREREG_TASK2_FLAG, PREREG_TASK3_FLAG


def perm_seed(line, day, pair):
    return int(MD4_SEED + 1000 * DONOR_INDEX[line] + 10 * int(day) + PAIR_INDEX[pair])


def boot_seed(line, day, pair):
    return int(MD4_BOOT + 1000 * DONOR_INDEX[line] + 10 * int(day) + PAIR_INDEX[pair])


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


def _disk_row_lines(title, path, max_cols=18, max_rows=80):
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
        "# PROGRESS_MD4",
        "",
        f"**STOP status:** {stop_s}",
        f"**Next action:** {next_action}",
        "",
        "## Seeds / gates",
        "",
        f"- seed `{MD4_SEED}`  boot `{MD4_BOOT}`  n_perm={N_PERM}  n_boot={N_BOOT}  n_random={N_RANDOM_DIR}",
        "- primary metric: Spearman ρ; uncalibrated R² is never a gate",
        f"- frozen ruler: `{FROZEN_RULER}` exists={FROZEN_RULER.exists()}  **not refit**",
        f"- MIN_CELLS_CONTRAST={MIN_CELLS_CONTRAST} (both states, per donor × timepoint; never pooled)",
        f"- Louvain resolution={LOUVAIN_RES} inherited from md2; AddModuleScore nbin={AMS_NBIN} ctrl={AMS_CTRL}",
        f"- PCA_MAHAL_K={PCA_MAHAL_K} (FINDINGS_FIBRO2.md measures)",
        f"- PREREG_TASK1.flag exists={PREREG_TASK1_FLAG.exists()}",
        f"- PREREG_TASK2.flag exists={PREREG_TASK2_FLAG.exists()}",
        f"- PREREG_TASK3.flag exists={PREREG_TASK3_FLAG.exists()}",
        f"- DECLARED_BEFORE_SCORES.flag exists={DECLARED_BEFORE_SCORES_FLAG.exists()}",
        "- GSE325735: out of scope",
        "- d7 is not the Stage 2 endpoint",
        "- Refit nothing. Do not re-cluster, re-label, or re-tune Louvain or AddModuleScore.",
        "- Donors never pooled. Timepoints never pooled.",
        "",
        "## Frozen ruler path",
        "",
        f"- `{FROZEN_RULER}` exists={FROZEN_RULER.exists()}",
        "",
        "## Clustering and scoring settings inherited from md2 (not re-tuned)",
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
        "## Gene-list versions and mapped counts (from disk; never from memory)",
        "",
    ]
    gs_path = MD4_DIR / "genesets_summary.json"
    if gs_path.exists():
        rec = load_json(gs_path)
        lines.append(f"- {json.dumps(jsonable(rec), ensure_ascii=False)[:8000]}")
        lines.append("")
    else:
        lines.append("- genesets_summary.json not yet written")
        lines.append("")
    idt = man.get("id_types") or {}
    if idt:
        lines += ["## ID type per matrix (from disk)", ""]
        for k, v in idt.items():
            lines.append(f"- **{k}** {json.dumps(jsonable(v), ensure_ascii=False)[:1500]}")
        lines.append("")
    lines += ["## Finished cells (from disk)", ""]
    for name, path in (
        ("Task 1 cell counts", MD4_DIR / "t1_cell_counts.csv"),
        ("Task 1 qualify", MD4_DIR / "t1_qualify.csv"),
        ("Task 1 sanity", MD4_DIR / "t1_sanity.json"),
        ("Task 1 state×timepoint means", MD4_DIR / "t1_state_timepoint_means.csv"),
        ("Task 1 contrasts", MD4_DIR / "t1_contrasts.csv"),
        ("Task 1 reading", MD4_DIR / "t1_reading.json"),
        ("Task 2 extrap", MD4_DIR / "t2_extrap.csv"),
        ("Task 2 correlation", MD4_DIR / "t2_correlation.csv"),
        ("Task 2 summary", MD4_DIR / "t2_summary.json"),
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
    for js_name in ("genesets_summary.json", "t1_summary.json", "t2_summary.json"):
        p = MD4_DIR / js_name
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
                f"- **{g.get('key')}** name=`{g.get('name')}` n={g.get('n')} "
                f"id_type={g.get('id_type')} source=`{g.get('source')}`"
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
        "- `src/md4_common.py`, `md4_task1.py`, `md4_task2.py`, `md4_findings.py`, `md4_run.py`",
        "- `results/md4/`",
        "- `FINDINGS_MD4.md`",
        "- `PROGRESS_MD4.md`",
        "",
    ]
    PROGRESS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return PROGRESS_PATH

"""Prompt L — PartialReprog vs Fibroblast (did responding cells get younger than where they started?).

Namespaced under results/md5/. Does not modify any existing FINDINGS_*.md or FALSIFICATION.md.
Does not open GSE325735. Refit nothing. Frozen ruler used as-is.
Seed 20260914. Bootstrap seed 20260918. Spearman ρ is the reported statistic.
Reuse src/md2_*.py Louvain labels + AddModuleScore exactly; do not re-cluster, re-label, or re-tune.
Reuse src/md4_*.py contrast code unchanged, called on different state pairs. Do not rescore cells.
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
from gtex_common import StopStep  # noqa: E402
from fibro_common import (  # noqa: E402
    FROZEN_RULER, BOOT_SEED, N_PERM, N_BOOT, N_RANDOM_DIR, jsonable,
)
from fibro3_common import AGED_LINE, YOUNG_LINE, DAYS  # noqa: E402
from fibro2_common import PCA_MAHAL_K  # noqa: E402
from md2_common import (  # noqa: E402
    MD2_DIR, MD2_PROC, MMC3_XLSX, AMS_NBIN, AMS_CTRL,
    LOUVAIN_RES, CLUSTER_SETTINGS,
)
from md3_common import MD3_DIR  # noqa: E402
from md4_common import (  # noqa: E402
    MD4_DIR, MD4_PROC, INSTRUMENTS, STATE_NAMES, PAIR_INDEX, DONOR_INDEX,
    MD2_POOLED_MD_PR, MD2_POOLED_MD_NR, MD2_N_PR, MD2_N_NR,
    MD3_TASK1_STATUS, MD3_TASK1_READING, MD3_TASK2_READING,
    MD2_NEGATIVE_HOLDS, STAGE2_VERDICT_FIBRO, FIG3G_QUOTE,
    perm_seed as md4_perm_seed, boot_seed as md4_boot_seed,
)

MD5_DIR = RESULTS / "md5"
MD5_FIG = MD5_DIR / "figures"
MD5_PROC = DATA_PROC / "md5"
for _p in (MD5_DIR, MD5_FIG, MD5_PROC):
    _p.mkdir(parents=True, exist_ok=True)

MD5_SEED = int(SEED)  # 20260914
MD5_BOOT = int(BOOT_SEED)  # 20260918
PREREG_DATE = "2026-09-19"

MANIFEST_PATH = MD5_DIR / "manifest.json"
PROGRESS_PATH = ROOT / "PROGRESS_MD5.md"
FINDINGS_PATH = ROOT / "FINDINGS_MD5.md"
PREREG_TASK1_FLAG = MD5_DIR / "PREREG_TASK1.flag"
PREREG_TASK2_FLAG = MD5_DIR / "PREREG_TASK2.flag"
PREREG_TASK3_FLAG = MD5_DIR / "PREREG_TASK3.flag"
DECLARED_BEFORE_SCORES_FLAG = MD5_DIR / "DECLARED_BEFORE_SCORES.flag"

MIN_CELLS_A = 30
SCHEME_B_FIB_DAYS = (0, 3)
SCHEME_B_PR_DAYS = (3, 7)
PRIMARY_PAIR = ("PartialReprog", "Fibroblast")
TASK2_PAIR = ("NonReprog", "Fibroblast")
TASK2_STATES = ("NonReprog", "Fibroblast", "PartialReprog")
QC_COLS = ("umi", "n_genes", "mito_frac")
CELL_CYCLE_COLS = ("phase", "Phase", "S.Score", "G2M.Score", "S_score", "G2M_score")

MD4_TASK1_YOUNGER_CLAIM = (
    "Both instruments agree that partially reprogrammed cells read younger. "
    "Their claim is supported and our earlier negatives reflected the wrong test "
    "(a within-state trajectory with no baseline)."
)
MD4_TASK1_STATUS_FIRED = (
    "Task 1: GM00731 d3: `both_lower`; d7: `both_lower`; GM23815 d10: `both_lower`. "
    "Fired: GM00731 `both_lower`; GM23815 `both_lower`."
)
MD4_CONTEXT_GM23815_D10_FIB_NR = (
    "on young GM23815 at d10, Fibroblast is also lower than NonReprog "
    "(ruler −0.294, MD −0.043)"
)
PAPER_NONREPROG_QUOTE = (
    "another is to the non-reprogrammed state that retains many of the starting "
    "fibroblast marker genes (Figure S6B), but it also expresses distinct markers "
    "perhaps due to viral infection and/or different culture conditions (Figure 3E)"
)
WHAT_WOULD_SETTLE = (
    "An independent fibroblast age instrument scored on these same Louvain-labelled "
    "PartialReprog vs Fibroblast cells at the same donor × timepoint, or the authors' "
    "Seurat v5 / sctransform v2 object with the frozen ruler projected onto that split."
)

PREREG_TASK1 = """Written **before any MD5 contrast p-value**, 2026-09-19. No threshold in this block is re-tuned after numbers exist.

Using the src/md2_*.py Louvain clustering already computed (not re-run, not re-tuned). Cell-state labels are the md2 assignment: argmax of mean AddModuleScore of mmc3 Reprog_cell_state_signatures columns Fibroblast / PartialReprog / EarlyPluripotency / Pluripotency / NonReprog. Labels are read from results/md2/t2_cluster_labels_*.csv. Not re-labelled. AddModuleScore for MD and TGF-β is the per-cell vector already stored in data/processed/md2/louvain_ams_*.npz (the md2 computation). Table S3 Age up / Age down and the frozen-ruler / pluripotency−fibroblast per-cell scores are the vectors already stored in data/processed/md4/cell_scores_*.npz (the md4 computation). Frozen ruler: results/fibro/frozen_ruler_ridge_raw.npz used as-is. Refit nothing. Do not re-cluster, re-label, re-tune, or rescore.

Report aged donor GM00731 and young donor GM23815 separately, never pooled. First report cell counts for Fibroblast and PartialReprog at every timepoint per donor. Then run the contrast under two pre-declared schemes, both reported, neither substituted for the other:

- (a) relaxed threshold: both states ≥30 cells at the same timepoint, no pooling across timepoints. Primary if it qualifies.
- (b) pooled timepoints within state: Fibroblast cells from d0+d3 against PartialReprog cells from d3+d7, each pooled within donor. This breaks md4's within-timepoint permutation guard, so the permutation must shuffle the state label within timepoint strata and the timepoint composition of both arms must be reported. Label this clearly as the weaker design; use it only where (a) does not qualify. Do not report (b) as primary where (a) qualified. Do not lower the threshold below 30. Do not pool donors.

Instruments, all applied unchanged from md4 cell_scores / md2 louvain_ams: frozen ruler, MD, TGF-β, Table S3 age-up minus age-down, pluripotency−fibroblast. Per contrast report Δmean and Δmedian, cell-level permutation p (200 shuffles, stratified as above for scheme (b); unstratified two-state shuffle among the cells in those two states at that donor × timepoint for scheme (a), which is md4_task1._one_contrast unchanged), bootstrap 95% CI (B=200), and Cohen's d (pooled sd, n−1, same sign as Δmean). Seed 20260914 for permutations. Bootstrap seed 20260918. Scheme (a) permutation/bootstrap offsets reuse md4_common.perm_seed / boot_seed for pair PartialReprog vs Fibroblast so already-published md4 context rows must reproduce. Scheme (b) uses a documented offset that does not collide with md4.

Δmean = mean(PartialReprog) − mean(Fibroblast). “Below” / “beating the null” = Δmean < 0 and p_lower ≤ 0.05 on that instrument's cell-level permutation null. Empirical p_lower = (n_null ≤ Δobs + 1) / (n_null + 1).

Sanity, before any new contrast p-value: (1) reproduce FINDINGS_MD2.md pooled MD means on the aged donor (PartialReprog +0.150, NonReprog +0.439; t2_reproduction.json); (2) occupancy matches results/md4/t1_cell_counts.csv. If either fails, stop. Cell-state labels have drifted.

**Pre-registered reading (only the outcome that fired, per donor):**
- PartialReprog below Fibroblast on the frozen ruler, beating its null → responding cells read younger than their own starting population. The rejuvenation reading holds and is not an artifact of NonReprog reading old. State this plainly.
- PartialReprog not below Fibroblast while md4's PartialReprog < NonReprog stands → the md4 result is driven by NonReprog reading old, not by PartialReprog reading young. The rejuvenation reading is not supported by this data, and md4's reading must be restated in those terms. Say so plainly and quote the md4 sentence at issue.
- Instruments disagree (MD below Fibroblast, ruler not, or the reverse) → report both, declare no winner, name what would settle it.
- Neither scheme qualifies on the aged donor → the decisive test cannot be run in this dataset. Say so; do not lower the threshold further, and do not let the young donor stand in for the aged one.

If two qualifying scheme-(a) timepoints in one donor fire different bullets, report both, do not pick, do not average. Scheme (b) is not averaged with scheme (a). Donors are not pooled.
"""

PREREG_TASK2 = """Written **before any MD5 QC median, cell-cycle fraction, extrapolation row, or NonReprog − Fibroblast contrast p-value**, 2026-09-19. Report-only. No gate.

Per donor × timepoint, for each of NonReprog, Fibroblast and PartialReprog with ≥30 Louvain-kept cells: median UMI, median genes, mitochondrial fraction, cell-cycle phase fractions if computable from the md2 object (louvain_obs / louvain_ams / cluster labels). If those objects have no phase / S.Score / G2M.Score columns, record that cell-cycle is not computable from the md2 object and do not substitute Seurat CellCycleScoring or a new Tirosh score.

Extrapolation distance to the GTEx fibroblast training distribution: the same PCA k=50 nearest-neighbour Euclidean and Mahalanobis measures as FINDINGS_FIBRO2.md / src/fibro2_extrap.py, already computed in results/md4/t2_extrap.csv (TMM among that donor's state×timepoint rows, then frozen μ/σ/w). Read those rows. Do not re-TMM, do not refit PCA, do not rescore.

Also report, for each instrument, the NonReprog − Fibroblast contrast at every timepoint where both states have ≥30 cells, calling md4_task1._one_contrast unchanged (unstratified within timepoint). If NonReprog reads consistently older than the starting fibroblasts across donors and timepoints, state that as an observation and note the paper's own remark about viral infection and culture conditions in those cells. Do not use Task 2 to discount Task 1.
"""

PREREG_TASK3 = """Written **before any MD5 contrast p-value**, 2026-09-19. Ledger only. Existing FINDINGS_*.md are not edited.

List, with the file each lives in, which readings currently stand and which are withdrawn:
- FINDINGS_FIBRO.md Stage 2 d0→d10 (stands, all-cell)
- FINDINGS_MD3.md Task 1 (withdrawn in md4)
- FINDINGS_MD3.md Task 2 Δρ (stands)
- FINDINGS_MD4.md Task 1 (stands, as qualified by Task 1 of this file)

Quote each; edit none of them. This is the ledger the write-up will be built from, so it must be exact.
"""


def md5_log_banner(log, stage: str):
    log("=" * 100)
    log(f"MD5 {stage}  seed={MD5_SEED}  boot_seed={MD5_BOOT}  "
        f"n_perm={N_PERM}  n_boot={N_BOOT}  n_random={N_RANDOM_DIR}")
    log("=" * 100)
    log("FALSIFICATION.md and existing FINDINGS*.md are not modified.")
    log(f"Frozen ruler: {FROZEN_RULER} exists={FROZEN_RULER.exists()}")
    log("Primary metric: Spearman ρ (n_perm=200); uncalibrated R² is never a gate.")
    log("Refit nothing. Do not open GSE325735. Do not promote d7 to the Stage 2 endpoint.")
    log("Fail loudly: no substitute columns, no inferred ages, no silent coerce.")
    log("Louvain / labels / AddModuleScore inherited from md2; not re-tuned.")
    log("md4 cell_scores reused; cells are not rescored. md4 _one_contrast reused.")
    log(f"MIN_CELLS_A={MIN_CELLS_A}. Scheme (b) Fib days={SCHEME_B_FIB_DAYS} "
        f"PR days={SCHEME_B_PR_DAYS}. Donors never pooled.")


def load_manifest():
    if MANIFEST_PATH.exists():
        return load_json(MANIFEST_PATH)
    return dict(
        seed=MD5_SEED, boot_seed=MD5_BOOT, files={}, failures=[],
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


def record_unusable(step, message, details=None):
    man = load_manifest()
    rec = dict(step=step, message=_clip_msg(message), details=jsonable(details or {}))
    rows = man.setdefault("unusable", [])
    if not any(f.get("step") == step and f.get("message") == rec["message"] for f in rows):
        rows.append(rec)
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


def perm_seed_b(line):
    """Scheme (b) offset 200: does not collide with md4 (10*day + pair ≤ 102)."""
    return int(MD5_SEED + 1000 * DONOR_INDEX[line] + 200)


def boot_seed_b(line):
    return int(MD5_BOOT + 1000 * DONOR_INDEX[line] + 200)


def perm_seed_t2(line, day):
    """NonReprog vs Fibroblast. Offset pair_i=20: does not collide with md4 pair 0..2."""
    return int(MD5_SEED + 1000 * DONOR_INDEX[line] + 10 * int(day) + 20)


def boot_seed_t2(line, day):
    return int(MD5_BOOT + 1000 * DONOR_INDEX[line] + 10 * int(day) + 20)


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
        "# PROGRESS_MD5",
        "",
        f"**STOP status:** {stop_s}",
        f"**Next action:** {next_action}",
        "",
        "## Seeds / gates",
        "",
        f"- seed `{MD5_SEED}`  boot `{MD5_BOOT}`  n_perm={N_PERM}  n_boot={N_BOOT}  n_random={N_RANDOM_DIR}",
        "- primary metric: Spearman ρ; uncalibrated R² is never a gate",
        f"- frozen ruler: `{FROZEN_RULER}` exists={FROZEN_RULER.exists()}  **not refit**",
        f"- MIN_CELLS_A={MIN_CELLS_A} (scheme a; both states, per donor × timepoint; never pooled)",
        f"- scheme (b) Fibroblast days={list(SCHEME_B_FIB_DAYS)} PartialReprog days={list(SCHEME_B_PR_DAYS)}; "
        "permutation stratified within timepoint; weaker design; not primary where (a) qualified",
        f"- Louvain resolution={LOUVAIN_RES} inherited from md2; AddModuleScore nbin={AMS_NBIN} ctrl={AMS_CTRL}",
        f"- PCA_MAHAL_K={PCA_MAHAL_K} (FINDINGS_FIBRO2.md measures; read from md4 t2_extrap.csv, not refit)",
        f"- PREREG_TASK1.flag exists={PREREG_TASK1_FLAG.exists()}",
        f"- PREREG_TASK2.flag exists={PREREG_TASK2_FLAG.exists()}",
        f"- PREREG_TASK3.flag exists={PREREG_TASK3_FLAG.exists()}",
        f"- DECLARED_BEFORE_SCORES.flag exists={DECLARED_BEFORE_SCORES_FLAG.exists()}",
        "- GSE325735: out of scope",
        "- d7 is not the Stage 2 endpoint",
        "- Refit nothing. Do not re-cluster, re-label, or re-tune Louvain or AddModuleScore.",
        "- Do not rescore cells. md4 cell_scores / md4 _one_contrast reused.",
        "- Donors never pooled. Scheme (a) timepoints never pooled.",
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
    gs_path = MD5_DIR / "genesets_summary.json"
    if gs_path.exists():
        rec = load_json(gs_path)
        lines.append(f"- {json.dumps(jsonable(rec), ensure_ascii=False)[:8000]}")
        lines.append("")
    else:
        md4_gs = MD4_DIR / "genesets_summary.json"
        if md4_gs.exists():
            rec = load_json(md4_gs)
            lines.append(f"- inherited from md4 `{md4_gs}`: "
                         f"{json.dumps(jsonable(rec), ensure_ascii=False)[:8000]}")
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
        ("Task 1 cell counts", MD5_DIR / "t1_cell_counts.csv"),
        ("Task 1 Fibroblast / PartialReprog counts", MD5_DIR / "t1_fib_pr_counts.csv"),
        ("Task 1 scheme (a) qualify", MD5_DIR / "t1_qualify_a.csv"),
        ("Task 1 scheme (b) composition", MD5_DIR / "t1_scheme_b_composition.csv"),
        ("Task 1 scheme (b) qualify", MD5_DIR / "t1_qualify_b.csv"),
        ("Task 1 sanity", MD5_DIR / "t1_sanity.json"),
        ("Task 1 scheme (a) contrasts", MD5_DIR / "t1_contrasts_a.csv"),
        ("Task 1 scheme (b) contrasts", MD5_DIR / "t1_contrasts_b.csv"),
        ("Task 1 reading", MD5_DIR / "t1_reading.json"),
        ("Task 2 QC", MD5_DIR / "t2_qc.csv"),
        ("Task 2 extrap", MD5_DIR / "t2_extrap.csv"),
        ("Task 2 NR − Fib contrasts", MD5_DIR / "t2_contrasts_nr_fib.csv"),
        ("Task 2 summary", MD5_DIR / "t2_summary.json"),
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
    if man.get("unusable"):
        lines += ["## Unusable / not computable (manifest)", ""]
        for f in man["unusable"]:
            lines.append(f"- **{f.get('step')}:** {_clip_msg(f.get('message'), 240)}")
        lines.append("")
    lines += [
        "## Files",
        "",
        "- `src/md5_common.py`, `md5_task1.py`, `md5_task2.py`, `md5_findings.py`, `md5_run.py`",
        "- `results/md5/`",
        "- `FINDINGS_MD5.md`",
        "- `PROGRESS_MD5.md`",
        "",
    ]
    PROGRESS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return PROGRESS_PATH

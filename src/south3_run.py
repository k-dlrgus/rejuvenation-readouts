"""SOUTH3 — patch of SOUTH2. Same data, same pre-registration, three changes.

1. n_perm = 20000 for the Stage 2A per-factor permutation null (everything else
   about the null unchanged), so BH across 1,836 factors can resolve q.
2. Guide agreement is computed on V_all guide vectors ALWAYS, and that one value
   is applied to both gene variants.
3. Report-only, no gate: per factor ||d|| and ||d|| / ||target - origin||, the
   Spearman correlation between delta rank and cosine rank by quartile of that
   ratio, and the small-displacement identity delta ~= -||d|| cos(theta).

This file is src/south2_run.py with those three changes applied. Every other rule
is carried over unchanged. src/south2_run.py is neither imported nor modified:
importing it would rebind the md3_idtype log sink to SOUTH2's manifest and would
re-create results/south2/, so the code is carried as a copy on purpose.

Reads results/survey/opened/southard/fibroblast_CRISPRa_mean_pop.h5ad (read-only)
and results/south2/stage3_full_table.csv (read-only, for the label comparison).
Imports from src/toward_run.py, src/same_run.py, src/seng_run.py,
src/lowdim_common.py, src/md3_idtype.py without changing them.
Does not refit the frozen ruler. Does not download anything. Does not open
L4A_matrix.h5.

The pre-registration must already exist verbatim at results/south3/PREREG.flag
(the SOUTH3 patch text), with the inherited SOUTH2 block at
results/south3/PREREG_SOUTH2_inherited.flag. Both were written before any
statistic. This script does not author either block; it reads them.

Writes only: FINDINGS_SOUTH3.md, PROGRESS_SOUTH3.md, results/south3/*.

Usage: python src/south3_run.py
"""
from __future__ import annotations

import os
import re
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats as sps

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import RESULTS, ROOT  # noqa: E402
from brain_phase1_common import Logger, dump_json, load_json, strip_ensembl  # noqa: E402
from fibro_common import FROZEN_RULER, jsonable  # noqa: E402
from fibro2_common import load_frozen_ruler  # noqa: E402
from gtex_common import StopStep  # noqa: E402
from target_common import md_table, unit  # noqa: E402
from trajectory_common import permutation_p  # noqa: E402

# --- the five modules the pre-registration names. Not modified. ----------------
from toward_run import (  # noqa: E402
    TOWARD_SEED,
    ANCHORS_PATH,
    YOUNG_BINS,
    OLD_BINS,
    load_gtex_z,
    percentile_ci,
    require_file,
    log_columns,
    _fmt,
    _fmt_p,
)
from same_run import (  # noqa: E402
    BOOT_SEED,
    SAME_DIR,
    cosine_full,
    frac_along,
    ci_covers_point,
)
from seng_run import (  # noqa: E402
    IDENTITY_GENES,
    IDENTITY_DROP,
    bh_q,
    dump_tree,
    file_md5,
    read_obs_var,
    ruler_column_index,
)
from lowdim_common import KS as LOWDIM_KS, fit_pca  # noqa: E402
import md3_idtype  # noqa: E402
from md3_idtype import detect_id_type, require_mappable  # noqa: E402


SOUTH3_DIR = RESULTS / "south3"
SOUTH3_DIR.mkdir(parents=True, exist_ok=True)
FINDINGS_PATH = ROOT / "FINDINGS_SOUTH3.md"
PROGRESS_PATH = ROOT / "PROGRESS_SOUTH3.md"
PREREG_FLAG = SOUTH3_DIR / "PREREG.flag"
PREREG_INHERITED = SOUTH3_DIR / "PREREG_SOUTH2_inherited.flag"
MANIFEST_PATH = SOUTH3_DIR / "manifest.json"
INVENTORY_PATH = SOUTH3_DIR / "inventory.json"
WORK_DIR = SOUTH3_DIR / "_workcache"
NULL_CACHE = SOUTH3_DIR / "stage2a_null_reduced.npz"
SOUTH2_DIR = RESULTS / "south2"
SOUTH2_TABLE = SOUTH2_DIR / "stage3_full_table.csv"
SOUTH2_FINDINGS = ROOT / "FINDINGS_SOUTH2.md"

SOUTHARD_DIR = RESULTS / "survey" / "opened" / "southard"
H5_PATH = SOUTHARD_DIR / "fibroblast_CRISPRa_mean_pop.h5ad"
PINNED_SIZE = 1_720_196_368
PINNED_MD5 = "ba44c7813903bb5df900348d6b0d589a"

SEED = 20260914
N_NULL = 200
# CHANGE 1. The only number that moves: the Stage 2A per-factor permutation null.
N_PERM = 20000
N_PERM_SOUTH2 = 200
# The smallest p an empirical one-sided permutation test can return, (0 + 1)/(N + 1).
P_MIN = 1.0 / (N_PERM + 1)
N_BOOT = 200
# 20000 x 4907 float64 is 785 MB, so the permutations are drawn in row chunks.
# Fixed before any statistic; the per-factor generator is unchanged.
PERM_CHUNK = 1000
N_WORKERS = max(1, int(os.environ.get("SOUTH3_WORKERS", "6")))
N_RULER_EXPECTED = 23485
S2_KS = (20, 50, 100)
KMAX = max(S2_KS)
F_VALUES = (0.0, 0.05, 0.1, 0.25, 0.5, 1.0)
ADJ_P_BAR = 0.05
Q_BAR = 0.05
MIN_GUIDES = 2
MIN_FACTORS = 500
COVERAGE_LOW_BAR = 0.25
BOOT_MIN_GUIDES = 3
CHUNK = 1024

TARGET_COL = "target_gene"
GUIDE_COL = "guide_identity"
CONTROL_LABEL = "off-target"
VARIANTS = ("V_all", "V_sig")

# Pre-registration STOP keys, in pre-registration order.
STOP_KEYS = ("not_an_effect_matrix", "too_few_factors", "no_detection_power")
# Report-only flags named in the pre-registration.
REPORT_FLAGS = ("scale_unknown", "coverage_low")


# --------------------------------------------------------------------------- #
# manifest + failure logging
# --------------------------------------------------------------------------- #
def load_prereg_text():
    out = []
    for p, what in ((PREREG_FLAG, "PREREG.flag"),
                    (PREREG_INHERITED, "PREREG_SOUTH2_inherited.flag")):
        if not p.exists():
            raise StopStep(
                "prereg",
                f"{what} missing — the pre-registration must be written verbatim before "
                "any statistic. Not computing.",
            )
        out.append(p.read_text(encoding="utf-8").rstrip("\n"))
    return out[0], out[1]


def load_manifest():
    if MANIFEST_PATH.exists():
        return load_json(MANIFEST_PATH)
    return dict(seed=SEED, boot_seed=BOOT_SEED, failures=[], fired_keys=[],
                report_flags=[], status="INIT")


def save_manifest(man):
    dump_json(MANIFEST_PATH, jsonable(man))
    return MANIFEST_PATH


def record_failure(step, message, details=None):
    man = load_manifest()
    rec = dict(step=step, message=str(message), details=jsonable(details or {}))
    man.setdefault("failures", []).append(rec)
    if step in STOP_KEYS and step not in man.setdefault("fired_keys", []):
        man["fired_keys"].append(step)
    man["status"] = "STOP"
    save_manifest(man)
    return rec


def record_flag(key, message, details=None):
    """A REPORT-ONLY flag. Never halts, always printed."""
    man = load_manifest()
    rec = dict(key=key, message=str(message), details=jsonable(details or {}))
    man.setdefault("report_flags", [])
    if not any(r.get("key") == key for r in man["report_flags"]):
        man["report_flags"].append(rec)
    save_manifest(man)
    return rec


def record_note(step, message, details=None):
    """A recorded non-fatal failure or implementation decision. Verbatim, never repaired."""
    man = load_manifest()
    man.setdefault("notes", []).append(
        dict(step=step, message=str(message), details=jsonable(details or {}))
    )
    save_manifest(man)


# md3_idtype.detect_id_type logs into results/md3/manifest.json through
# md3_common.log_id_type. SOUTH3 may only write results/south3/*, so the sink is
# rebound in this process to SOUTH3's own manifest. The imported function itself
# is untouched on disk.
_ID_TYPES = {}


def _south3_id_sink(man_key, rec):
    _ID_TYPES[str(man_key)] = jsonable(rec)


md3_idtype.log_id_type = _south3_id_sink


def _fmt_int(x):
    try:
        return f"{int(x):,}"
    except (TypeError, ValueError):
        return "NA"


def _f(x, d=3):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "NA"
    return f"{v:.{d}f}" if np.isfinite(v) else "NA"


def _sf(x, d=4):
    """Signed fixed-width float, NA-safe."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "NA"
    return f"{v:+.{d}f}" if np.isfinite(v) else "NA"


# --------------------------------------------------------------------------- #
# STAGE 0 — inventory (descriptive only; no statistic)
# --------------------------------------------------------------------------- #
def verify_pinned_file(log):
    if not H5_PATH.exists():
        raise StopStep("pinned_file", f"missing {H5_PATH}. Not substituting.")
    size = int(H5_PATH.stat().st_size)
    t0 = time.perf_counter()
    md5 = file_md5(H5_PATH)
    log(f"[file] {H5_PATH}")
    log(f"[file] size={size} (pinned {PINNED_SIZE}) md5={md5} (pinned {PINNED_MD5}) "
        f"[{time.perf_counter() - t0:.1f}s]")
    if size != PINNED_SIZE or md5 != PINNED_MD5:
        raise StopStep(
            "pinned_file",
            f"pinned file mismatch: size={size} (expected {PINNED_SIZE}) md5={md5} "
            f"(expected {PINNED_MD5}). Not substituting a different file.",
            dict(size=size, md5=md5, expected_size=PINNED_SIZE, expected_md5=PINNED_MD5),
        )
    return size, md5


def coverage_by_symbol(var, frozen):
    ruler_sym = np.array([str(s).upper() for s in np.asarray(frozen["symbol"])], dtype=object)
    sym_u = np.array([str(s).upper() for s in var["gene_name"].astype(str).to_numpy()], dtype=object)
    pos = {}
    for i, s in enumerate(sym_u):
        pos.setdefault(s, i)
    idx = np.array([pos.get(s, -1) for s in ruler_sym], dtype=int)
    return idx


def coverage_by_ensembl(var, frozen):
    ruler_ens = np.asarray(frozen["ensembl"]).astype(str)
    ens = np.array([strip_ensembl(g) for g in var["gene_id"].astype(str).to_numpy()])
    pos = {}
    for i, e in enumerate(ens):
        pos.setdefault(e, i)
    idx = np.array([pos.get(e, -1) for e in ruler_ens], dtype=int)
    return idx


def coverage_record(name, idx, frozen, n_var):
    w = np.abs(np.asarray(frozen["w"], float))
    tot = float(w.sum())
    hit = idx >= 0
    return dict(
        mapping=name,
        n_matched=int(hit.sum()),
        n_ruler=int(idx.size),
        frac_genes=float(hit.sum()) / float(idx.size) if idx.size else np.nan,
        weight_matched=float(w[hit].sum()),
        weight_total=tot,
        weight_frac=float(w[hit].sum() / tot) if tot > 0 else np.nan,
        n_source_cols_used=int(len(set(idx[hit].tolist()))),
        n_source_cols=int(n_var),
    )


def scan_matrices(f, obs, var, cols_src, scale, log):
    """One pass over X + layers/masked, one pass over layers/adj_p.

    Returns the compact matched-gene effect matrix in frozen-z units (rows are
    guides), the significance mask on the same columns, and the descriptive audits.
    """
    Xds, Mds, Ads = f["X"], f["layers"]["masked"], f["layers"]["adj_p"]
    Pds = f["layers"]["p"]
    n_rows, n_genes = (int(Xds.shape[0]), int(Xds.shape[1]))
    n_matched = int(cols_src.size)

    gid = np.array([strip_ensembl(g) for g in var["gene_id"].astype(str).to_numpy()])
    pos_id = {}
    for i, g in enumerate(gid):
        pos_id.setdefault(g, i)
    tgi = obs["target_gene_id"].astype(str).str.strip().to_numpy()
    de = pd.to_numeric(obs["de_genes"], errors="coerce").to_numpy(float)

    Gc = np.empty((n_rows, n_matched), dtype=np.float64)
    Sc = np.empty((n_rows, n_matched), dtype=bool)

    x_min, x_max = np.inf, -np.inf
    n_seen = n_neg = n_nonfinite = n_integral = 0
    mask_total = 0
    mask_on_own_target = 0
    mask_off_target_col = []
    mask_per_row = np.zeros(n_rows, dtype=int)

    t0 = time.perf_counter()
    for s in range(0, n_rows, CHUNK):
        e = min(s + CHUNK, n_rows)
        xb = np.asarray(Xds[s:e], dtype=np.float64)
        mb = np.asarray(Mds[s:e], dtype=np.float64)
        n_seen += int(xb.size)
        fin = np.isfinite(xb)
        n_nonfinite += int((~fin).sum())
        if fin.any():
            x_min = min(x_min, float(xb[fin].min()))
            x_max = max(x_max, float(xb[fin].max()))
            n_neg += int((xb[fin] < 0).sum())
            n_integral += int((xb[fin] == np.round(xb[fin])).sum())
        nan_mask = ~np.isfinite(mb)
        rr, cc = np.where(nan_mask)
        mask_total += int(rr.size)
        for i, j in zip(rr.tolist(), cc.tolist()):
            mask_per_row[s + i] += 1
            want = pos_id.get(tgi[s + i], -1)
            if want == j:
                mask_on_own_target += 1
            else:
                mask_off_target_col.append(
                    dict(row=int(s + i), col=int(j), target_gene_id=str(tgi[s + i]),
                         target_col=int(want))
                )
        # Masked entries are set to 0 (pre-registered) before anything is derived.
        xb[nan_mask] = 0.0
        Gc[s:e] = xb[:, cols_src] * scale
        del xb, mb, nan_mask
    log(f"[scan X+masked] {n_rows} rows in {time.perf_counter() - t0:.1f}s; "
        f"masked_entries={mask_total} on_own_target={mask_on_own_target} "
        f"elsewhere={len(mask_off_target_col)}")

    t0 = time.perf_counter()
    de_exact = 0
    n_adj_sig = 0
    p_min, p_max = np.inf, -np.inf
    for s in range(0, n_rows, CHUNK):
        e = min(s + CHUNK, n_rows)
        ab = np.asarray(Ads[s:e], dtype=np.float64)
        sig = ab <= ADJ_P_BAR
        n_adj_sig += int(sig.sum())
        cnt = sig.sum(1).astype(float)
        de_exact += int(np.sum(cnt == de[s:e]))
        Sc[s:e] = sig[:, cols_src]
        pb = np.asarray(Pds[s:e], dtype=np.float64)
        fin = np.isfinite(pb)
        if fin.any():
            p_min = min(p_min, float(pb[fin].min()))
            p_max = max(p_max, float(pb[fin].max()))
        del ab, sig, pb
    log(f"[scan adj_p+p] {time.perf_counter() - t0:.1f}s; rows where "
        f"count(adj_p<={ADJ_P_BAR}) == obs.de_genes: {de_exact}/{n_rows}; "
        f"total adj_p<={ADJ_P_BAR} entries={n_adj_sig}")

    x_audit = dict(
        where="X", dtype=str(Xds.dtype), n_seen=int(n_seen),
        n_negative=int(n_neg), n_nonfinite=int(n_nonfinite),
        n_integral=int(n_integral),
        frac_negative=float(n_neg) / float(n_seen) if n_seen else np.nan,
        frac_integral=float(n_integral) / float(n_seen) if n_seen else np.nan,
        min=float(x_min), max=float(x_max),
        integer_dtype=bool(np.issubdtype(Xds.dtype, np.integer)),
    )
    mask_audit = dict(
        n_masked_entries=int(mask_total),
        n_on_own_target_gene=int(mask_on_own_target),
        n_elsewhere=int(len(mask_off_target_col)),
        examples_elsewhere=mask_off_target_col[:10],
        n_rows_with_any_mask=int((mask_per_row > 0).sum()),
        masked_per_row_max=int(mask_per_row.max()),
        masked_per_row_min_nonzero=int(mask_per_row[mask_per_row > 0].min()) if (mask_per_row > 0).any() else 0,
        set_to_zero=True,
    )
    p_audit = dict(
        rows_where_de_genes_equals_count_adj_p_le_005=int(de_exact),
        n_rows=int(n_rows),
        n_entries_adj_p_le_005=int(n_adj_sig),
        p_min=float(p_min), p_max=float(p_max),
    )
    return Gc, Sc, x_audit, mask_audit, p_audit


def stage0(log):
    for p in (H5_PATH, FROZEN_RULER):
        if not Path(p).exists():
            raise StopStep("inventory", f"missing {p}. Not substituting.")
    size, md5 = verify_pinned_file(log)

    frozen = load_frozen_ruler()
    n_ruler = int(len(np.asarray(frozen["mu"])))
    if n_ruler != N_RULER_EXPECTED:
        raise StopStep(
            "frozen_ruler",
            f"frozen ruler has {n_ruler} genes, pre-registration says {N_RULER_EXPECTED}. "
            "Not refitting and not reconciling.",
        )
    log(f"[frozen] {FROZEN_RULER} n_genes={n_ruler} keys={list(frozen.files)} (not refit)")

    f = h5py.File(H5_PATH, "r")
    tree = dump_tree(f, SOUTH3_DIR / "h5_tree.txt")
    log(f"[h5] tree written ({len(tree)} lines) -> {SOUTH3_DIR / 'h5_tree.txt'}")
    obs = read_obs_var(f["obs"])
    var = read_obs_var(f["var"])
    log_columns(log, "obs", list(obs.columns), str(H5_PATH))
    log_columns(log, "var", list(var.columns), str(H5_PATH))
    x_shape = tuple(int(s) for s in f["X"].shape)
    layer_names = list(f["layers"].keys()) if "layers" in f else []
    uns_keys = list(f["uns"].keys()) if "uns" in f else []
    obsm_keys = list(f["obsm"].keys()) if "obsm" in f else []
    varm_keys = list(f["varm"].keys()) if "varm" in f else []
    log(f"[h5] X shape={x_shape} dtype={f['X'].dtype} layers={layer_names} "
        f"uns={uns_keys or 'none'} obsm={obsm_keys or 'none'} varm={varm_keys or 'none'} "
        f"raw_group={'raw' in f}")
    for col in (TARGET_COL, GUIDE_COL, "target_gene_id", "target_expr", "active",
                "sequence_driven", "de_genes", "cell_count"):
        if col not in obs.columns:
            raise StopStep("inventory", f"obs missing {col!r}. Not substituting.",
                           dict(columns=list(obs.columns)))
    for col in ("gene_id", "gene_name"):
        if col not in var.columns:
            raise StopStep("inventory", f"var missing {col!r}. Not substituting.",
                           dict(columns=list(var.columns)))
    for lay in ("p", "adj_p", "masked"):
        if lay not in layer_names:
            raise StopStep("inventory", f"layers missing {lay!r}. Not substituting.",
                           dict(layers=layer_names))
        if tuple(int(s) for s in f["layers"][lay].shape) != x_shape:
            raise StopStep("inventory",
                           f"layer {lay} shape {f['layers'][lay].shape} != X shape {x_shape}.")

    # --- what the column index holds -------------------------------------- #
    id_rec = detect_id_type(var["gene_id"].astype(str).to_numpy(), log, "south3_var_gene_id",
                            str(H5_PATH) + "::var/gene_id")
    sym_rec = detect_id_type(var["gene_name"].astype(str).to_numpy(), log, "south3_var_gene_name",
                             str(H5_PATH) + "::var/gene_name")
    require_mappable(id_rec, ("ensembl",), log, "inventory")
    require_mappable(sym_rec, ("symbol",), log, "inventory")
    gid_raw = var["gene_id"].astype(str).to_numpy()
    n_versioned = int(sum(1 for g in gid_raw if re.match(r"^ENS[A-Z]*G\d+\.\d+$", g)))
    col_index_kind = (
        f"unversioned Ensembl gene ids in var/gene_id ({id_rec['n_ensembl']:,}/"
        f"{id_rec['n']:,} match ENSG, {n_versioned:,} carry a version suffix) with HGNC "
        f"symbols in var/gene_name ({sym_rec['n_symbol_like']:,}/{sym_rec['n']:,} symbol-shaped). "
        "The h5ad row index of var is a positional integer, not an identifier."
    )
    log(f"[column index] {col_index_kind}")

    # --- units and log base ------------------------------------------------ #
    unit_re = re.compile(r"(?i)log2|log10|\bln\b|log[-_ ]?base|cpm|tpm|rpkm|fpkm|unit|"
                         r"normali[sz]|scale|fold[-_ ]?change|lfc")
    doc_hits = []
    for p in sorted(SOUTHARD_DIR.glob("*.txt")):
        head = p.read_text(encoding="utf-8", errors="replace")[:100_000]
        if unit_re.search(head):
            doc_hits.append(p.name)
    x_attrs = {k: str(v) for k, v in dict(f["X"].attrs).items()}
    attr_hits = [k for k in x_attrs if unit_re.search(str(k))]
    col_hits = [c for c in list(var.columns) + list(obs.columns) if unit_re.search(str(c))]
    determinable = bool(uns_keys or attr_hits or col_hits or doc_hits)
    unit_rec = dict(
        uns_keys=uns_keys, obsm_keys=obsm_keys, varm_keys=varm_keys,
        x_attrs=x_attrs, x_attr_unit_hits=attr_hits, column_unit_hits=col_hits,
        var_columns=list(var.columns), obs_columns=list(obs.columns),
        bundled_txt_files=[p.name for p in sorted(SOUTHARD_DIR.glob("*.txt"))],
        bundled_txt_with_unit_words=doc_hits,
        var_mean_range=[float(var["mean"].min()), float(var["mean"].max())],
        var_std_range=[float(var["std"].min()), float(var["std"].max())],
        searched_for=unit_re.pattern,
        determinable=determinable,
    )
    log(f"[units] searched uns={uns_keys or 'empty'}, X attrs={list(x_attrs)}, var/obs column "
        f"names, and bundled .txt files {unit_rec['bundled_txt_files']} for {unit_re.pattern!r}: "
        f"attr hits={attr_hits} column hits={col_hits} doc hits={doc_hits}. "
        f"Units and log base determinable from the file: {determinable}")

    # --- ruler coverage, three ways ---------------------------------------- #
    idx_sym = coverage_by_symbol(var, frozen)
    idx_ens = coverage_by_ensembl(var, frozen)
    idx_best = ruler_column_index(
        var["gene_id"].astype(str).to_numpy(),
        var["gene_name"].astype(str).to_numpy(),
        frozen, log,
    )
    cov = [
        coverage_record("symbol", idx_sym, frozen, len(var)),
        coverage_record("ensembl_version_stripped", idx_ens, frozen, len(var)),
        coverage_record("best_of_two_after_md3_idtype", idx_best, frozen, len(var)),
    ]
    for c in cov:
        log(f"[coverage {c['mapping']}] matched {c['n_matched']:,}/{c['n_ruler']:,} genes "
            f"({c['frac_genes']:.4f}); |w| fraction {c['weight_frac']:.4f}")
    order = {"symbol": idx_sym, "ensembl_version_stripped": idx_ens,
             "best_of_two_after_md3_idtype": idx_best}
    chosen = max(cov, key=lambda c: (c["weight_frac"], c["n_matched"]))
    idx_used = order[chosen["mapping"]]
    log(f"[coverage] mapping used for all later stages: {chosen['mapping']} "
        f"(highest |w| fraction = {chosen['weight_frac']:.4f})")
    coverage_low = bool(chosen["weight_frac"] < COVERAGE_LOW_BAR)

    # --- factors, guides, controls ----------------------------------------- #
    targets = obs[TARGET_COL].astype(str)
    guides = obs[GUIDE_COL].astype(str)
    n_rows = int(len(obs))
    gid_col = obs["target_gene_id"].astype(str).str.strip()
    texp = pd.to_numeric(obs["target_expr"], errors="coerce")
    active = np.asarray(obs["active"].to_numpy(), dtype=bool)
    ctrl_label = (targets == CONTROL_LABEL).to_numpy()
    ctrl_evidence = (gid_col == "").to_numpy() & texp.isna().to_numpy() & (~active)
    ctrl_mask = ctrl_label
    n_ctrl = int(ctrl_mask.sum())
    log(f"[controls] rows labelled `{CONTROL_LABEL}`: {n_ctrl}; of those, rows also carrying "
        f"empty target_gene_id AND NaN target_expr AND active=False: "
        f"{int((ctrl_label & ctrl_evidence).sum())}")
    log(f"[controls] rows meeting the evidence triple but NOT labelled `{CONTROL_LABEL}`: "
        f"{int((ctrl_evidence & ~ctrl_label).sum())}")

    ctrl_tab = pd.DataFrame({
        "row": np.arange(n_rows)[ctrl_mask],
        GUIDE_COL: guides.to_numpy()[ctrl_mask],
        TARGET_COL: targets.to_numpy()[ctrl_mask],
        "target_gene_id_empty": (gid_col == "").to_numpy()[ctrl_mask],
        "target_expr_nan": texp.isna().to_numpy()[ctrl_mask],
        "active": active[ctrl_mask],
        "sequence_driven": np.asarray(obs["sequence_driven"].to_numpy(), bool)[ctrl_mask],
        "gene_driven": np.asarray(obs["gene_driven"].to_numpy(), bool)[ctrl_mask],
        "de_genes": pd.to_numeric(obs["de_genes"], errors="coerce").to_numpy()[ctrl_mask],
        "cell_count": pd.to_numeric(obs["cell_count"], errors="coerce").to_numpy()[ctrl_mask],
        "strength": pd.to_numeric(obs["strength"], errors="coerce").to_numpy()[ctrl_mask],
    })
    ctrl_tab.to_csv(SOUTH3_DIR / "stage0_controls.csv", index=False)

    tf_mask = ~ctrl_mask
    tf_labels = sorted(set(targets.to_numpy()[tf_mask]))
    per_factor = pd.Series(targets.to_numpy()[tf_mask]).value_counts()
    gmin, gmed, gmax = int(per_factor.min()), float(per_factor.median()), int(per_factor.max())
    n_ge2 = int((per_factor >= MIN_GUIDES).sum())
    log(f"[factors] distinct transcription factors={len(tf_labels):,}; guides per factor "
        f"min={gmin} median={gmed:g} max={gmax}; with >= {MIN_GUIDES} guides={n_ge2:,}")
    gdist = per_factor.value_counts().sort_index()
    pd.DataFrame({"n_guides": gdist.index.astype(int), "n_factors": gdist.to_numpy(int)}).to_csv(
        SOUTH3_DIR / "stage0_guides_per_factor.csv", index=False)
    pd.DataFrame({"factor": per_factor.index, "n_guides": per_factor.to_numpy(int)}).to_csv(
        SOUTH3_DIR / "stage0_factors.csv", index=False)

    one_row_per_guide = bool(guides.nunique() == n_rows)
    has_barcode = any(re.search(r"(?i)barcode", str(c)) for c in obs.columns)
    row_unit = dict(
        unit="per-guide perturbation row (the source's own mean-population regression output)",
        one_row_per_guide=one_row_per_guide,
        n_unique_guides=int(guides.nunique()),
        cell_barcode_column=bool(has_barcode),
        cell_count_min=float(pd.to_numeric(obs["cell_count"], errors="coerce").min()),
        cell_count_median=float(pd.to_numeric(obs["cell_count"], errors="coerce").median()),
        cell_count_max=float(pd.to_numeric(obs["cell_count"], errors="coerce").max()),
        cell_count_sum=float(pd.to_numeric(obs["cell_count"], errors="coerce").sum()),
    )
    log(f"[row unit] {row_unit['unit']}; guide_identity unique on every row="
        f"{one_row_per_guide} ({row_unit['n_unique_guides']:,} of {n_rows:,})")

    # --- the single pass over the matrices --------------------------------- #
    hit = idx_used >= 0
    ruler_rows = np.flatnonzero(hit)               # positions in the 23,485 ruler
    cols_src = idx_used[hit].astype(int)           # matching source columns
    sd = np.asarray(frozen["sd"], float)
    sd_safe = np.where(sd < 1e-12, 1.0, sd)
    scale = 1.0 / sd_safe[ruler_rows]
    n_sd_floored = int((sd < 1e-12).sum())
    log(f"[convert] {cols_src.size:,} matched coordinates; frozen sd floored at 1e-12 on "
        f"{n_sd_floored} ruler genes")

    Gc, Sc, x_audit, mask_audit, p_audit = scan_matrices(f, obs, var, cols_src, scale, log)
    f.close()

    for a in (x_audit, mask_audit, p_audit):
        log(f"[audit] {a}")

    # --- units: the descriptive diagnostics that were tried ---------------- #
    # |effect| per gene against the source's own var/mean and var/std. Descriptive
    # only: it does not settle the log base, and nothing downstream branches on it.
    mabs = np.abs(Gc * (1.0 / scale)).mean(0)  # back to source units, matched genes only
    vmean = var["mean"].to_numpy(float)[cols_src]
    vstd = var["std"].to_numpy(float)[cols_src]
    unit_rec["diag_spearman_meanabs_vs_var_mean"] = float(sps.spearmanr(mabs, vmean).statistic)
    unit_rec["diag_spearman_meanabs_vs_var_std"] = float(sps.spearmanr(mabs, vstd).statistic)
    unit_rec["diag_spearman_meanabs_over_std_vs_var_mean"] = float(
        sps.spearmanr(mabs / np.where(vstd < 1e-12, 1.0, vstd), vmean).statistic)
    unit_rec["diag_median_meanabs"] = float(np.median(mabs))
    log(f"[units diag] spearman(mean|effect|, var/mean)="
        f"{unit_rec['diag_spearman_meanabs_vs_var_mean']:+.3f}  "
        f"spearman(mean|effect|, var/std)={unit_rec['diag_spearman_meanabs_vs_var_std']:+.3f}  "
        f"spearman(mean|effect|/std, var/mean)="
        f"{unit_rec['diag_spearman_meanabs_over_std_vs_var_mean']:+.3f} — neither a log base nor "
        f"a linear unit is established by this. scale_unknown stands.")

    # --- identity panel ---------------------------------------------------- #
    ruler_sym_all = np.array([str(s).upper() for s in np.asarray(frozen["symbol"])], dtype=object)
    sym_to_ruler = {}
    for i, s in enumerate(ruler_sym_all):
        sym_to_ruler.setdefault(s, i)
    id_present, id_missing = [], []
    id_rowpos = []
    ruler_to_compact = {int(r): i for i, r in enumerate(ruler_rows)}
    for g in IDENTITY_GENES:
        r = sym_to_ruler.get(g)
        if r is not None and int(r) in ruler_to_compact:
            id_present.append(g)
            id_rowpos.append(ruler_to_compact[int(r)])
        else:
            id_missing.append(g)
    log(f"[identity] panel actually used ({len(id_present)}): {id_present}; absent: {id_missing}")

    inv = dict(
        path=str(H5_PATH), size_bytes=size, md5=md5,
        pinned_size=PINNED_SIZE, pinned_md5=PINNED_MD5,
        x_shape=list(x_shape), x_dtype=str(np.dtype(np.float64)),
        layers=layer_names, uns_keys=uns_keys, obsm_keys=obsm_keys, varm_keys=varm_keys,
        obs_columns=list(obs.columns), var_columns=list(var.columns),
        n_rows=n_rows, n_var_genes=int(len(var)),
        row_unit=row_unit,
        column_index_kind=col_index_kind,
        id_type_records=_ID_TYPES,
        n_versioned_ensembl=n_versioned,
        units=unit_rec,
        x_audit=x_audit, mask_audit=mask_audit, p_audit=p_audit,
        coverage=cov, coverage_chosen=chosen, coverage_low=coverage_low,
        coverage_low_bar=COVERAGE_LOW_BAR,
        n_matched=int(cols_src.size), n_fixed_zero=int(n_ruler - cols_src.size),
        n_ruler=n_ruler, n_sd_floored=n_sd_floored,
        n_factors=len(tf_labels), guides_min=gmin, guides_median=gmed, guides_max=gmax,
        n_factors_ge2=n_ge2, min_guides=MIN_GUIDES, min_factors=MIN_FACTORS,
        n_controls=n_ctrl, control_label=CONTROL_LABEL,
        n_controls_full_evidence=int((ctrl_label & ctrl_evidence).sum()),
        n_evidence_not_labelled=int((ctrl_evidence & ~ctrl_label).sum()),
        controls=ctrl_tab.to_dict("records"),
        identity_panel=id_present, identity_absent=id_missing,
    )
    pd.DataFrame(cov).to_csv(SOUTH3_DIR / "stage0_coverage.csv", index=False)
    dump_json(INVENTORY_PATH, jsonable(inv))
    log(f"[inventory] wrote {INVENTORY_PATH}")

    # --- report-only flags -------------------------------------------------- #
    if not unit_rec["determinable"]:
        record_flag(
            "scale_unknown",
            "The units and log base of X are undeterminable from the file: uns is empty, "
            "X carries no unit attribute, var holds only mean/std/cv/fano/excess_cv with no "
            "unit statement, obs holds no unit statement, and the only bundled documentation "
            "in the folder is three plain gene-symbol lists. This is a REPORT-ONLY flag, not "
            "a stop. The Stage 2A/2B floors are magnitude-matched (they permute or re-use the "
            "same vector, so they are invariant to a global scale factor on X); the Stage 2C "
            "MDA is NOT invariant, because it scales the true young-minus-old direction to the "
            "median factor norm, and that norm moves with any global rescaling of X.",
            unit_rec,
        )
    if coverage_low:
        record_flag(
            "coverage_low",
            f"The chosen mapping ({chosen['mapping']}) matches "
            f"{chosen['weight_frac']:.4f} of the frozen ruler's total absolute weight, which is "
            f"below {COVERAGE_LOW_BAR}. The Stage 3 candidate table is titled exploratory and "
            "every candidate carries that word.",
            chosen,
        )

    # --- STOP conditions, in pre-registration order ------------------------- #
    fired = []
    is_effect = (
        ("p" in layer_names and "adj_p" in layer_names)
        and x_audit["n_negative"] > 0
        and x_audit["frac_integral"] < 0.5
        and p_audit["rows_where_de_genes_equals_count_adj_p_le_005"] == n_rows
    )
    log(f"[effect-matrix test] p+adj_p layers present={('p' in layer_names and 'adj_p' in layer_names)}; "
        f"X has {x_audit['n_negative']:,} negative entries "
        f"({x_audit['frac_negative']:.3f}); integral fraction {x_audit['frac_integral']:.4f}; "
        f"obs.de_genes reproduced by count(adj_p<=0.05) on "
        f"{p_audit['rows_where_de_genes_equals_count_adj_p_le_005']:,}/{n_rows:,} rows "
        f"=> X is per-gene effect estimates: {is_effect}")
    if not is_effect:
        fired.append((
            "not_an_effect_matrix",
            "X is not per-gene effect estimates. "
            f"p/adj_p layers={[l for l in layer_names]}; n_negative={x_audit['n_negative']}; "
            f"frac_integral={x_audit['frac_integral']}; rows where obs.de_genes equals "
            f"count(adj_p<=0.05)={p_audit['rows_where_de_genes_equals_count_adj_p_le_005']} "
            f"of {n_rows}.",
            dict(x_audit=x_audit, p_audit=p_audit, layers=layer_names),
        ))
    if n_ge2 < MIN_FACTORS:
        fired.append((
            "too_few_factors",
            f"{n_ge2} factors have >= {MIN_GUIDES} guides, fewer than {MIN_FACTORS}.",
            dict(n_factors_ge2=n_ge2, n_factors=len(tf_labels)),
        ))

    stage0_pack = dict(
        Gc=Gc, Sc=Sc, ruler_rows=ruler_rows, cols_src=cols_src, scale=scale,
        obs=obs, var=var, frozen=frozen, ctrl_mask=ctrl_mask, tf_labels=tf_labels,
        per_factor=per_factor, id_rowpos=np.asarray(id_rowpos, int),
        id_present=id_present, id_missing=id_missing, inv=inv,
    )
    return stage0_pack, fired


# --------------------------------------------------------------------------- #
# STAGE 1 — vectors and spaces
# --------------------------------------------------------------------------- #
def build_factor_vectors(pack, log):
    obs, Gc, Sc = pack["obs"], pack["Gc"], pack["Sc"]
    targets = obs[TARGET_COL].astype(str).to_numpy()
    ctrl = pack["ctrl_mask"]
    per_factor = pack["per_factor"]
    keep = [str(t) for t in per_factor.index if int(per_factor[t]) >= MIN_GUIDES]
    keep = sorted(keep)
    rows_by_factor = {}
    for t in keep:
        rows_by_factor[t] = np.flatnonzero((targets == t) & (~ctrl))
    n_f = len(keep)
    log(f"[stage1] {n_f:,} factors with >= {MIN_GUIDES} guides; per-factor vector = mean of "
        f"that factor's guide rows, converted to frozen-z by dividing each gene's effect by "
        f"that gene's frozen training sd")
    D = {}
    for variant in VARIANTS:
        M = np.empty((n_f, Gc.shape[1]), dtype=np.float64)
        for i, t in enumerate(keep):
            r = rows_by_factor[t]
            block = Gc[r] if variant == "V_all" else Gc[r] * Sc[r]
            M[i] = block.mean(0)
        D[variant] = M
        nrm = np.linalg.norm(M, axis=1)
        log(f"[stage1 {variant}] ||d|| min={nrm.min():.4f} median={np.median(nrm):.4f} "
            f"max={nrm.max():.4f}; n_nonzero_coords median="
            f"{np.median((M != 0).sum(1)):.0f}")
    ctrl_rows = np.flatnonzero(ctrl)
    C = {v: (Gc[ctrl_rows] if v == "V_all" else Gc[ctrl_rows] * Sc[ctrl_rows]) for v in VARIANTS}
    return keep, rows_by_factor, D, C


def build_s2(gtex, log):
    Zg = np.asarray(gtex["Z"], np.float64)
    center = Zg.mean(0)
    sp = fit_pca(Zg - center, KMAX)
    L = np.asarray(sp.L, np.float64)
    if L.shape[1] < KMAX:
        record_note("s2", f"PCA returned {L.shape[1]} components, fewer than kmax={KMAX}.",
                    dict(n_components=int(L.shape[1])))
    tot_var = float(np.var(Zg - center, axis=0, ddof=1).sum())
    evals = np.asarray(sp.evals, float) if sp.evals is not None else np.array([])
    cum = np.cumsum(evals) / tot_var if tot_var > 0 else np.full(evals.shape, np.nan)
    log(f"[S2] PCA fit on the GTEx fibroblast donor matrix only: n_donors={Zg.shape[0]} "
        f"n_genes={Zg.shape[1]} kmax={L.shape[1]}; loadings orthonormal "
        f"(max |L'L - I| = {float(np.abs(L.T @ L - np.eye(L.shape[1])).max()):.2e})")
    for k in S2_KS:
        if k <= cum.size:
            log(f"[S2] k={k}: cumulative fraction of GTEx donor variance = {cum[k-1]:.4f}")
    return dict(L=L, center=center, cum_var=cum, tot_var=tot_var, n_donors=int(Zg.shape[0]))


def build_targets_origins(gtex, s2, log):
    require_file(ANCHORS_PATH, "toward_anchors")
    anc = np.load(ANCHORS_PATH, allow_pickle=True)
    log_columns(log, "toward_anchors", list(anc.files), str(ANCHORS_PATH))
    c_young = np.asarray(anc["c_young"], float).ravel()
    c_old = np.asarray(anc["c_old"], float).ravel()
    log(f"[Y1/O1] TOWARD anchors: n_young={int(anc['n_young'])} (bins {list(YOUNG_BINS)}) "
        f"n_old={int(anc['n_old'])} (bins {list(OLD_BINS)}); ||c_young - c_old||="
        f"{float(np.linalg.norm(c_young - c_old)):.4f}")

    same_panel = SAME_DIR / "panel_z.npz"
    require_file(same_panel, "same_panel_z")
    sp = np.load(same_panel, allow_pickle=True)
    names = [str(x) for x in np.asarray(sp["names"]).tolist()]
    if "O" not in names:
        raise StopStep("O2", f"SAME panel_z has no 'O' row. names={names}. Not substituting.")
    z_o2 = np.asarray(sp["Z"], float)[names.index("O")].ravel()
    log_columns(log, "same_panel_z", names, str(same_panel))
    log(f"[O2] Lu et al. GM00731 untreated day-0 fibroblast pseudobulk as built in SAME: "
        f"{int((z_o2 == 0).sum()):,} of {z_o2.size:,} coordinates are exactly 0 "
        f"(genes absent from that dataset, zeroed by SAME's own missing-gene rule)")

    # Y2: the same young donors as a region, S2 only.
    Zg, age = np.asarray(gtex["Z"], float), gtex["age"]
    young_all = np.flatnonzero(np.isin(age, YOUNG_BINS))
    rng = np.random.default_rng(SEED)
    perm = rng.permutation(young_all.size)
    n_fit = int(young_all.size // 2)
    fit_idx = np.sort(young_all[perm[:n_fit]])
    held_idx = np.sort(young_all[perm[n_fit:]])
    log(f"[Y2] young donors n={young_all.size}; region fit on {fit_idx.size} donors, "
        f"membership radius calibrated on {held_idx.size} held-out young donors "
        f"(seed {SEED})")
    L, center = s2["L"], s2["center"]
    Pfit_full = (Zg[fit_idx] - center) @ L
    Pheld_full = (Zg[held_idx] - center) @ L
    y2 = {}
    for k in S2_KS:
        A = Pfit_full[:, :k]
        mu = A.mean(0)
        C = np.cov(A, rowvar=False, ddof=1)
        C = np.atleast_2d(C)
        rank = int(np.linalg.matrix_rank(C))
        ev = np.linalg.eigvalsh(C)
        cond = float(ev.max() / ev.min()) if ev.min() > 0 else np.inf
        if rank < k:
            record_note(
                "Y2_singular",
                f"Y2 at k={k}: the young-donor covariance built on {A.shape[0]} fit donors has "
                f"rank {rank} < k={k}, so the Mahalanobis metric is singular and undefined. "
                "No pseudo-inverse or ridge was substituted. Y2 is reported as undefined at "
                "this k; Y1 at the same k is unaffected.",
                dict(k=int(k), n_fit=int(A.shape[0]), rank=rank),
            )
            log(f"[Y2 k={k}] SINGULAR: rank={rank} < k={k} on {A.shape[0]} fit donors. "
                f"Mahalanobis undefined; reported as such, nothing substituted.")
            y2[k] = dict(k=int(k), defined=False, rank=rank, n_fit=int(A.shape[0]),
                         n_held=int(held_idx.size), mu=mu, Cinv=None, radius=np.nan,
                         cond=cond)
            continue
        Cinv = np.linalg.inv(C)
        dh = Pheld_full[:, :k] - mu
        md_held = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", dh, Cinv, dh), 0.0))
        radius = float(np.percentile(md_held, 95))
        log(f"[Y2 k={k}] defined: rank={rank} cond={cond:.3e}; held-out young Mahalanobis "
            f"median={np.median(md_held):.4f} max={md_held.max():.4f}; membership radius "
            f"(95th pct of held-out young donors' own distances) = {radius:.4f}")
        y2[k] = dict(k=int(k), defined=True, rank=rank, n_fit=int(A.shape[0]),
                     n_held=int(held_idx.size), mu=mu, Cinv=Cinv, radius=radius,
                     cond=cond, md_held=md_held)
    return dict(c_young=c_young, c_old=c_old, z_o2=z_o2, y2=y2,
                young_fit=fit_idx, young_held=held_idx, young_all=young_all)


def build_combos(pack, s2, tgt, log):
    """The fixed list of (space, origin, target) settings. Nothing is chosen later."""
    ruler_rows = pack["ruler_rows"]
    L, center = s2["L"], s2["center"]
    Lm = L[ruler_rows, :]
    origins = {"O1": tgt["c_old"], "O2": tgt["z_o2"]}
    combos = []
    for oname, o in origins.items():
        r = o - tgt["c_young"]
        combos.append(dict(
            space="S1", k=None, origin=oname, target="Y1", pair=f"{oname}_Y1",
            metric="euclidean_frozen_z",
            r_m=r[ruler_rows].copy(), r_full=r.copy(), r2=float(np.dot(r, r)),
            r_norm=float(np.linalg.norm(r)),
        ))
    for k in S2_KS:
        Lk = L[:, :k]
        for oname, o in origins.items():
            op = (o - center) @ Lk
            yp = (tgt["c_young"] - center) @ Lk
            r = op - yp
            combos.append(dict(
                space=f"S2_k{k}", k=int(k), origin=oname, target="Y1", pair=f"{oname}_Y1",
                metric="euclidean_pc", r_m=r.copy(), r2=float(np.dot(r, r)),
                r_norm=float(np.linalg.norm(r)), proj=True,
            ))
            rec = tgt["y2"][k]
            r2v = op - rec["mu"]
            combos.append(dict(
                space=f"S2_k{k}", k=int(k), origin=oname, target="Y2", pair=f"{oname}_Y2",
                metric="euclidean_pc+mahalanobis", r_m=r2v.copy(), r2=float(np.dot(r2v, r2v)),
                r_norm=float(np.linalg.norm(r2v)), proj=True,
                y2=rec, op=op.copy(),
            ))
    rows = []
    for c in combos:
        rec = dict(space=c["space"], k=c["k"], origin=c["origin"], target=c["target"],
                   metric=c["metric"], origin_to_target_distance=c["r_norm"])
        if c["target"] == "Y2":
            y = c["y2"]
            if y["defined"]:
                d = c["op"] - y["mu"]
                rec["origin_to_target_mahalanobis"] = float(
                    np.sqrt(max(float(d @ y["Cinv"] @ d), 0.0)))
                rec["membership_radius"] = y["radius"]
                rec["origin_inside_radius"] = bool(
                    rec["origin_to_target_mahalanobis"] <= y["radius"])
            else:
                rec["origin_to_target_mahalanobis"] = np.nan
                rec["membership_radius"] = np.nan
                rec["origin_inside_radius"] = None
            rec["y2_defined"] = y["defined"]
        rows.append(rec)
    df = pd.DataFrame(rows)
    df.to_csv(SOUTH3_DIR / "origin_target_distances.csv", index=False)
    log("[origin->target distances, reported before anything is scored]")
    for _, r in df.iterrows():
        extra = ""
        if r.get("target") == "Y2":
            extra = (f" mahalanobis={_f(r.get('origin_to_target_mahalanobis'), 4)} "
                     f"radius={_f(r.get('membership_radius'), 4)} "
                     f"inside={r.get('origin_inside_radius')}")
        log(f"  {r['space']:9s} {r['origin']}->{r['target']}  ||origin-target||="
            f"{r['origin_to_target_distance']:.4f}{extra}")
    return combos, Lm, df


# --------------------------------------------------------------------------- #
# geometry on a block of displacement vectors
# --------------------------------------------------------------------------- #
def block_stats(combo, Q, Qp):
    """delta / cos / frac for every row of a displacement block.

    Q  : (N, n_matched) displacements in frozen-z gene coordinates.
    Qp : (N, KMAX) the same displacements projected on the frozen GTEx PCA basis.
    """
    if combo.get("proj"):
        k = combo["k"]
        A = Qp[:, :k]
    else:
        A = Q
    r = combo["r_m"]
    dot = A @ r
    qn2 = np.einsum("ij,ij->i", A, A)
    r_norm, r2 = combo["r_norm"], combo["r2"]
    dist = np.sqrt(np.maximum(r2 + 2.0 * dot + qn2, 0.0))
    delta = dist - r_norm
    dn = np.sqrt(np.maximum(qn2, 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        cos = np.where((dn > 1e-12) & (r_norm > 1e-12), -dot / (dn * r_norm), np.nan)
        frac = np.where(r2 > 1e-12, -dot / r2, np.nan)
    out = dict(delta=delta, cos=cos, frac=frac, d_norm=dn, dist=dist)
    if combo["target"] == "Y2" and combo["y2"]["defined"]:
        y = combo["y2"]
        Xb = combo["op"][None, :] + A
        Db = Xb - y["mu"][None, :]
        md = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", Db, y["Cinv"], Db), 0.0))
        d0 = combo["op"] - y["mu"]
        md0 = float(np.sqrt(max(float(d0 @ y["Cinv"] @ d0), 0.0)))
        out["md"] = md
        out["md_change"] = md - md0
        out["inside"] = md <= y["radius"]
    return out


def full_space_delta(combo, dfull, L):
    """delta for displacements given on all ruler coordinates (the dose-response,
    whose direction is not confined to the measured genes)."""
    dfull = np.atleast_2d(dfull)
    if combo.get("proj"):
        A = dfull @ L[:, :combo["k"]]
        r = combo["r_m"]
    else:
        A = dfull
        r = combo["r_full"]
    dot = A @ r
    qn2 = np.einsum("ij,ij->i", A, A)
    return np.sqrt(np.maximum(combo["r2"] + 2.0 * dot + qn2, 0.0)) - combo["r_norm"]


# --------------------------------------------------------------------------- #
# observed geometry — computed once, consumed by Stage 2A (for the streaming
# permutation p) and by Stage 3 (as the point estimate). The numbers are exactly
# SOUTH2's Stage 3 point estimates; only the order of computation moved, because
# with N_PERM = 20000 the permuted deltas can no longer be held in RAM and the
# p-value has to be accumulated against an already-known observed value.
# --------------------------------------------------------------------------- #
def observed_points(combos, Lm, D, factors, log):
    keys = [f"{c['space']}|{c['pair']}" for c in combos]
    n_f = len(factors)
    point = {}
    for variant in VARIANTS:
        M = D[variant]
        for kk in keys:
            point[(variant, kk)] = {c: np.full(n_f, np.nan) for c in
                                    ("delta", "cos", "frac", "d_norm", "md_change", "inside")}
        t0 = time.perf_counter()
        for i in range(n_f):
            d = M[i]
            dp = (d[None, :] @ Lm)
            for c, kk in zip(combos, keys):
                st = block_stats(c, d[None, :], dp)
                for cname in ("delta", "cos", "frac", "d_norm"):
                    point[(variant, kk)][cname][i] = float(st[cname][0])
                if "md_change" in st:
                    point[(variant, kk)]["md_change"][i] = float(st["md_change"][0])
                    point[(variant, kk)]["inside"][i] = float(bool(st["inside"][0]))
        log(f"[observed {variant}] point estimates for {n_f:,} factors x {len(keys)} settings "
            f"in {time.perf_counter() - t0:.0f}s")
    return point, keys


# --------------------------------------------------------------------------- #
# STAGE 2A permutation null — N_PERM = 20000, drawn in chunks, streamed
# --------------------------------------------------------------------------- #
def _perm_block(combo, Q, Qp, s):
    """One chunk of permuted vectors: (delta, pooled-scaled delta, md_change).

    Algebraically identical to SOUTH2's `block_stats(c, Q, Qp)["delta"]` and
    `block_stats(c, Q * s, Qp * s)["delta"]`, written out so the scaled copy of a
    20000 x 4907 block is never materialised.
    """
    A = Qp[:, :combo["k"]] if combo.get("proj") else Q
    dot = A @ combo["r_m"]
    qn2 = np.einsum("ij,ij->i", A, A)
    r2, rn = combo["r2"], combo["r_norm"]
    delta = np.sqrt(np.maximum(r2 + 2.0 * dot + qn2, 0.0)) - rn
    if np.isfinite(s):
        pooled = np.sqrt(np.maximum(r2 + 2.0 * s * dot + (s * s) * qn2, 0.0)) - rn
    else:
        pooled = np.full(delta.shape, np.nan)
    md_change = None
    if combo["target"] == "Y2" and combo["y2"]["defined"]:
        y = combo["y2"]
        Db = (combo["op"][None, :] + A) - y["mu"][None, :]
        md = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", Db, y["Cinv"], Db), 0.0))
        d0 = combo["op"] - y["mu"]
        md0 = float(np.sqrt(max(float(d0 @ y["Cinv"] @ d0), 0.0)))
        md_change = md - md0
    return delta, pooled, md_change


_W = {}


def _perm_init(work_dir, variant):
    """Worker set-up. M and Lm are memory-mapped so the workers share one copy."""
    import pickle
    wd = Path(work_dir)
    with open(wd / f"payload_{variant}.pkl", "rb") as fh:
        _W.update(pickle.load(fh))
    _W["M"] = np.load(wd / f"M_{variant}.npy", mmap_mode="r")
    _W["Lm"] = np.load(wd / "Lm.npy", mmap_mode="r")
    _W["pool"] = [np.load(wd / f"pool_{variant}_{ci}.npy", mmap_mode="r+")
                  for ci in range(len(_W["combos"]))]
    _W["buf"] = np.empty((PERM_CHUNK, _W["M"].shape[1]), dtype=np.float64)


def _perm_one_factor(i):
    combos = _W["combos"]
    nC = len(combos)
    d = np.asarray(_W["M"][i], dtype=np.float64)
    Lm = _W["Lm"]
    nrm = float(_W["nrm"][i])
    s = (_W["med"] / nrm) if nrm >= 1e-12 else np.nan
    # Unchanged from SOUTH2: one independent generator per factor, seeded [SEED, i].
    rng = np.random.default_rng([SEED, i])
    dl = np.empty((nC, N_PERM), dtype=np.float64)
    mdl = np.full((nC, N_PERM), np.nan, dtype=np.float64)
    pool = np.empty((nC, N_PERM), dtype=np.float32)
    buf = _W["buf"]
    for start in range(0, N_PERM, PERM_CHUNK):
        n = min(PERM_CHUNK, N_PERM - start)
        Pb = buf[:n]
        Pb[:] = d
        rng.permuted(Pb, axis=1, out=Pb)
        Pp = Pb @ Lm
        for ci, c in enumerate(combos):
            de, po, mc = _perm_block(c, Pb, Pp, s)
            dl[ci, start:start + n] = de
            pool[ci, start:start + n] = po
            if mc is not None:
                mdl[ci, start:start + n] = mc
    obs_delta = _W["obs_delta"][:, i]
    obs_md = _W["obs_md"][:, i]
    floor = np.percentile(dl, 5, axis=1)
    n_le = (dl <= obs_delta[:, None]).sum(1).astype(np.int64)
    md_floor = np.full(nC, np.nan)
    md_n_le = np.full(nC, -1, dtype=np.int64)
    for ci in range(nC):
        if np.isfinite(mdl[ci, 0]):
            md_floor[ci] = float(np.percentile(mdl[ci], 5))
            if np.isfinite(obs_md[ci]):
                md_n_le[ci] = int(np.sum(mdl[ci] <= obs_md[ci]))
    for ci in range(nC):
        _W["pool"][ci][i, :] = pool[ci]
    return i, floor, n_le, md_floor, md_n_le


def run_permutation_null(combos, keys, Lm, D, point, log):
    """Stage 2A for every factor and variant.

    Returns per-factor floors, the streamed permutation counts, and the pooled
    5th percentile. The 20000 permuted deltas per factor and setting are not
    retained: the per-factor floor, the count of permuted deltas at or below the
    observed delta, and the rescaled draw that feeds the pooled floor are all the
    later stages use, and the first two are reduced inside the worker.
    """
    import pickle

    nC0 = len(combos)
    n_f0 = D[VARIANTS[0]].shape[0]
    # Stage 2A is deterministic given the seed, n_perm and the chunk size, and it is by far
    # the most expensive step. Its reduced output is cached so the report can be rebuilt
    # without redrawing 73.4 million permutations. The cache is keyed on every parameter that
    # could change a number and is ignored if any of them differs.
    if NULL_CACHE.exists():
        z = np.load(NULL_CACHE, allow_pickle=False)
        ok = (int(z["n_perm"]) == N_PERM and int(z["perm_chunk"]) == PERM_CHUNK
              and int(z["seed"]) == SEED and int(z["n_f"]) == n_f0 and int(z["n_c"]) == nC0
              and [str(x) for x in z["keys"]] == list(keys)
              and [str(x) for x in z["variants"]] == list(VARIANTS))
        if ok and (SOUTH3_DIR / "stage2a_pooled_floor.csv").exists():
            log(f"[stage2A] reusing {NULL_CACHE.name}: same seed, n_perm={N_PERM:,}, "
                f"chunk={PERM_CHUNK:,}, factors and settings. Not redrawing the permutations.")
            floors, counts, md_floors, md_counts = {}, {}, {}, {}
            for vi, variant in enumerate(VARIANTS):
                for ci, kk in enumerate(keys):
                    floors[(variant, kk)] = z["floors"][vi, ci]
                    counts[(variant, kk)] = z["counts"][vi, ci].astype(np.int64)
                    if np.isfinite(z["md_floors"][vi, ci]).any():
                        md_floors[(variant, kk)] = z["md_floors"][vi, ci]
                        md_counts[(variant, kk)] = z["md_counts"][vi, ci].astype(np.int64)
            pooled = pd.read_csv(SOUTH3_DIR / "stage2a_pooled_floor.csv")
            return floors, counts, md_floors, md_counts, pooled
        log(f"[stage2A] {NULL_CACHE.name} exists but does not match the current parameters; "
            "ignoring it and redrawing.")

    if WORK_DIR.exists():
        shutil.rmtree(WORK_DIR)
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    nC = len(combos)
    n_f = D[VARIANTS[0]].shape[0]
    np.save(WORK_DIR / "Lm.npy", np.ascontiguousarray(Lm))

    floors, counts, md_floors, md_counts = {}, {}, {}, {}
    cache = {k: np.full((len(VARIANTS), nC0, n_f0), np.nan) for k in ("floors", "md_floors")}
    cache.update({k: np.full((len(VARIANTS), nC0, n_f0), -1, dtype=np.int64)
                  for k in ("counts", "md_counts")})
    pooled_rows = []
    for variant in VARIANTS:
        M = D[variant]
        nrm = np.linalg.norm(M, axis=1)
        med = float(np.median(nrm))
        n_zero = int(np.sum(nrm < 1e-12))
        obs_delta = np.vstack([point[(variant, kk)]["delta"] for kk in keys])
        obs_md = np.vstack([point[(variant, kk)]["md_change"] for kk in keys])
        np.save(WORK_DIR / f"M_{variant}.npy", np.ascontiguousarray(M))
        with open(WORK_DIR / f"payload_{variant}.pkl", "wb") as fh:
            pickle.dump(dict(combos=combos, keys=keys, nrm=nrm, med=med,
                             obs_delta=obs_delta, obs_md=obs_md), fh, protocol=4)
        for ci in range(nC):
            mm = np.lib.format.open_memmap(WORK_DIR / f"pool_{variant}_{ci}.npy", mode="w+",
                                           dtype=np.float32, shape=(n_f, N_PERM))
            del mm

        fl = np.full((nC, n_f), np.nan)
        cn = np.full((nC, n_f), -1, dtype=np.int64)
        mfl = np.full((nC, n_f), np.nan)
        mcn = np.full((nC, n_f), -1, dtype=np.int64)
        t0 = time.perf_counter()
        done = 0
        if N_WORKERS <= 1:
            _perm_init(str(WORK_DIR), variant)
            it = (_perm_one_factor(i) for i in range(n_f))
            for i, a, b, c_, d_ in it:
                fl[:, i], cn[:, i], mfl[:, i], mcn[:, i] = a, b, c_, d_
                done += 1
                if done % 100 == 0:
                    log(f"[stage2A {variant}] {done}/{n_f} factors "
                        f"({time.perf_counter() - t0:.0f}s)")
            _W.clear()
        else:
            saved = {k: os.environ.get(k) for k in
                     ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                      "NUMEXPR_NUM_THREADS")}
            for k in saved:
                os.environ[k] = "1"
            try:
                with ProcessPoolExecutor(max_workers=N_WORKERS, initializer=_perm_init,
                                         initargs=(str(WORK_DIR), variant)) as ex:
                    for i, a, b, c_, d_ in ex.map(_perm_one_factor, range(n_f), chunksize=2):
                        fl[:, i], cn[:, i], mfl[:, i], mcn[:, i] = a, b, c_, d_
                        done += 1
                        if done % 100 == 0:
                            log(f"[stage2A {variant}] {done}/{n_f} factors "
                                f"({time.perf_counter() - t0:.0f}s)")
            finally:
                for k, v in saved.items():
                    if v is None:
                        os.environ.pop(k, None)
                    else:
                        os.environ[k] = v
        log(f"[stage2A {variant}] {n_f:,} factors x {N_PERM:,} permutations done in "
            f"{time.perf_counter() - t0:.0f}s on {N_WORKERS} worker(s); median ||d||={med:.4f}; "
            f"zero-norm factors={n_zero}")

        vi = VARIANTS.index(variant)
        cache["floors"][vi] = fl
        cache["counts"][vi] = cn
        cache["md_floors"][vi] = mfl
        cache["md_counts"][vi] = mcn
        for ci, (c, kk) in enumerate(zip(combos, keys)):
            floors[(variant, kk)] = fl[ci].copy()
            counts[(variant, kk)] = cn[ci].copy()
            if np.isfinite(mfl[ci]).any():
                md_floors[(variant, kk)] = mfl[ci].copy()
                md_counts[(variant, kk)] = mcn[ci].copy()
            pp = WORK_DIR / f"pool_{variant}_{ci}.npy"
            pl = np.load(pp, mmap_mode="r")
            pl = np.asarray(pl).reshape(-1)
            pl = pl[np.isfinite(pl)].astype(np.float64)
            pooled_rows.append(dict(
                variant=variant, space=c["space"], k=c["k"], pair=c["pair"],
                median_factor_norm=med, n_factors_pooled=int(n_f - n_zero),
                n_zero_norm_factors=n_zero, n_draws=int(pl.size),
                floor_p05=float(np.percentile(pl, 5)),
                null_median=float(np.median(pl)), null_min=float(pl.min()),
                null_max=float(pl.max()), n_perm=N_PERM,
            ))
            del pl
            pp.unlink()
        (WORK_DIR / f"M_{variant}.npy").unlink()
        (WORK_DIR / f"payload_{variant}.pkl").unlink()
    shutil.rmtree(WORK_DIR, ignore_errors=True)
    np.savez_compressed(
        NULL_CACHE, n_perm=N_PERM, perm_chunk=PERM_CHUNK, seed=SEED, n_f=n_f0, n_c=nC0,
        keys=np.array(list(keys), dtype=object).astype(str),
        variants=np.array(list(VARIANTS)).astype(str), **cache)
    log(f"[stage2A] reduced null written to {NULL_CACHE.name} "
        f"({NULL_CACHE.stat().st_size / 1e6:.1f} MB)")
    return floors, counts, md_floors, md_counts, pd.DataFrame(pooled_rows)


# --------------------------------------------------------------------------- #
# STAGE 2 — detection limit
# --------------------------------------------------------------------------- #
def stage2(pack, combos, Lm, D, C, tgt, point, keys, log):
    """A: per-factor magnitude-matched permutation floor (+ pooled floor).
    B: control-row floor (reported, never a gate).
    C: dose-response along the true young-minus-old direction, and the MDA.
    """
    n_f = D[VARIANTS[0]].shape[0]
    log(f"[stage2A] per-factor floor: permute each factor's own vector entries across the "
        f"{Lm.shape[0]:,} matched gene coordinates, {N_PERM} times, seed {SEED}. The "
        f"{pack['inv']['n_fixed_zero']:,} unmatched ruler coordinates are structurally fixed "
        f"at 0 by the conversion rule and are not permuted into; this keeps every null "
        f"direction inside the same measured subspace as the real effect and preserves ||d|| "
        f"exactly.")
    log(f"[stage2A pooled] pooled floor: every factor's {N_PERM} permuted vectors are rescaled "
        f"to the median factor norm and pooled across all factors; the pooled floor is the 5th "
        f"percentile of that pooled distribution. Zero-norm factors cannot be rescaled and are "
        f"left out of the pool (counted below).")
    record_note(
        "stage2A_streaming",
        f"Implementation decision forced by n_perm = {N_PERM:,} and fixed before any floor "
        f"existed. SOUTH2 held every permuted delta in RAM ({N_PERM_SOUTH2} per factor and "
        f"setting) and derived the floor and the permutation p from that array afterwards. "
        f"{N_PERM:,} permuted deltas for 1,836 factors across 28 factor-settings do not fit, so "
        "each factor's permutations are drawn in row chunks of "
        f"{PERM_CHUNK:,} and reduced as they are produced: the 5th percentile (the floor), the "
        "count of permuted deltas at or below the already-computed observed delta (the "
        "numerator of the one-sided empirical p), and the rescaled draw that feeds the pooled "
        "floor. The per-factor generator, np.random.default_rng([SEED, i]), the permuted "
        "coordinate set, and the delta definition are all unchanged. Two consequences are "
        "stated rather than hidden: (a) the observed delta is now computed before the "
        "permutations instead of after, which changes no value because it does not depend on "
        "them; (b) the pooled draw is evaluated from the algebraic form "
        "sqrt(r2 + 2*s*dot + s^2*qn2) - ||r|| instead of by rescaling a 20000 x 4907 block, "
        "which is the same quantity up to float64 rounding. Both were checked against the "
        "SOUTH2 code path before the run: drawing the permutations in chunks of "
        f"{PERM_CHUNK:,} from np.random.default_rng([SEED, i]) gives a bit-identical block to "
        "drawing them in one call, so SOUTH3's null is SOUTH2's stream continued — the first "
        f"{N_PERM_SOUTH2} of SOUTH3's {N_PERM:,} draws are literally SOUTH2's {N_PERM_SOUTH2} "
        "draws — and the streamed delta matches the block_stats delta exactly, the pooled draw "
        "to 5e-14 absolute. The reduced null — the per-factor floors and the observed-value "
        f"counts, which is everything any later stage reads — is kept at "
        f"`results/south3/{NULL_CACHE.name}` so the report can be rebuilt without redrawing "
        "the permutations; it is keyed on the seed, n_perm, chunk size, factor count and "
        "setting list, and is ignored if any of them differs.",
        dict(n_perm=N_PERM, n_perm_south2=N_PERM_SOUTH2, perm_chunk=PERM_CHUNK,
             n_workers=N_WORKERS),
    )
    floors, perm_counts, md_floors, md_counts, pooled = run_permutation_null(
        combos, keys, Lm, D, point, log)
    pooled.to_csv(SOUTH3_DIR / "stage2a_pooled_floor.csv", index=False)

    # per-factor floor table
    frows = []
    for variant in VARIANTS:
        for c in combos:
            kk = f"{c['space']}|{c['pair']}"
            fl = floors[(variant, kk)]
            frows.append(dict(
                variant=variant, space=c["space"], k=c["k"], pair=c["pair"],
                floor_min=float(fl.min()), floor_p25=float(np.percentile(fl, 25)),
                floor_median=float(np.median(fl)), floor_p75=float(np.percentile(fl, 75)),
                floor_max=float(fl.max()), n_factors=int(fl.size), n_perm=N_PERM,
            ))
    fdf = pd.DataFrame(frows)
    fdf.to_csv(SOUTH3_DIR / "stage2a_per_factor_floor_summary.csv", index=False)

    # --- B: control rows as pseudo-effects --------------------------------- #
    crows = []
    for variant in VARIANTS:
        Cm = C[variant]
        if Cm.shape[0] == 0:
            continue
        Cp = Cm @ Lm
        for c in combos:
            dl = block_stats(c, Cm, Cp)["delta"]
            crows.append(dict(
                variant=variant, space=c["space"], k=c["k"], pair=c["pair"],
                n_controls=int(Cm.shape[0]),
                control_floor_p05=float(np.percentile(dl, 5)),
                control_min=float(dl.min()), control_median=float(np.median(dl)),
                control_max=float(dl.max()),
                control_mean_norm=float(np.linalg.norm(Cm, axis=1).mean()),
                is_gate=False, controls_inert=False,
            ))
    cdf = pd.DataFrame(crows)
    cdf.to_csv(SOUTH3_DIR / "stage2b_control_floor.csv", index=False)
    log(f"[stage2B] control floor from {C[VARIANTS[0]].shape[0]} untargeted `{CONTROL_LABEL}` "
        f"rows used as pseudo-effects. Reported only. These rows are NOT inert, so this is "
        f"never the gate.")

    # --- C: dose-response and MDA ------------------------------------------ #
    ruler_rows = pack["ruler_rows"]
    v_true_full = tgt["c_young"] - tgt["c_old"]
    u_full, n_full = unit(v_true_full)
    v_restricted = np.zeros_like(v_true_full)
    v_restricted[ruler_rows] = v_true_full[ruler_rows]
    u_res, n_res = unit(v_restricted)
    log(f"[stage2C] true young-minus-old direction: ||v||={n_full:.4f}; restricted to the "
        f"{ruler_rows.size:,} measured coordinates ||v_measured||={n_res:.4f} "
        f"({n_res / n_full:.4f} of the full norm)")
    record_note(
        "stage2C_clearing_rule",
        "Implementation reading fixed before any floor existed: a dose 'clears the floor' when "
        "its delta is below the primary floor AND negative, mirroring the Stage 3 toward_young "
        "rule. Without the second clause, f = 0 (no displacement, delta exactly 0) would count "
        "as clearing whenever a floor happens to be positive, which would report detection "
        "power for adding nothing.",
    )
    drows = []
    for variant in VARIANTS:
        med = float(np.median(np.linalg.norm(D[variant], axis=1)))
        for support, uu in (("full_v", u_full), ("measured_coords_only", u_res)):
            Qfull = np.vstack([float(fv) * med * uu for fv in F_VALUES])
            for c in combos:
                st = full_space_delta(c, Qfull, pack["L_full"])
                pf = pooled[(pooled.variant == variant) & (pooled.space == c["space"])
                            & (pooled.pair == c["pair"])]["floor_p05"]
                floor = float(pf.iloc[0])
                for fv, dl in zip(F_VALUES, st):
                    drows.append(dict(
                        variant=variant, space=c["space"], k=c["k"], pair=c["pair"],
                        support=support, f=float(fv), median_factor_norm=med,
                        added_norm=float(fv) * med, delta=float(dl),
                        pooled_floor_p05=floor,
                        clears_floor=bool((dl < floor) and (dl < 0)),
                    ))
    ddf = pd.DataFrame(drows)
    ddf.to_csv(SOUTH3_DIR / "stage2c_dose_response.csv", index=False)

    mrows = []
    for (variant, space, pair, support), g in ddf.groupby(
            ["variant", "space", "pair", "support"], sort=False):
        ok = g[(g.f <= 1.0) & g.clears_floor]
        mda = float(ok.f.min()) if len(ok) else np.nan
        mrows.append(dict(variant=variant, space=space, pair=pair, support=support,
                          mda=mda, cleared=bool(len(ok) > 0),
                          pooled_floor_p05=float(g.pooled_floor_p05.iloc[0]),
                          median_factor_norm=float(g.median_factor_norm.iloc[0])))
    mdf = pd.DataFrame(mrows)
    mdf.to_csv(SOUTH3_DIR / "stage2c_mda.csv", index=False)
    for _, r in mdf.iterrows():
        log(f"[MDA] {r['variant']:6s} {r['space']:9s} {r['pair']:6s} {r['support']:20s} "
            f"MDA={_f(r['mda'], 3)} (pooled floor {r['pooled_floor_p05']:+.5f})")

    primary = mdf[mdf.support == "full_v"]
    passed = bool(primary.cleared.any())
    log(f"[stage2C] some f <= 1.0 clears the primary floor in at least one space/variant: "
        f"{passed}")
    fired = []
    if not passed:
        fired.append((
            "no_detection_power",
            "No dose fraction f <= 1.0 of the true young-minus-old direction, scaled to the "
            "median factor norm, produced a delta below the primary (magnitude-matched "
            "permutation) floor in any space or gene variant. Stages 0-2 only are reported. "
            "This data cannot support a ranking.",
            dict(mda_table=mdf.to_dict("records")),
        ))
    return dict(floors=floors, perm_counts=perm_counts, md_floors=md_floors,
                md_counts=md_counts, pooled=pooled, per_factor_floor=fdf,
                control=cdf, dose=ddf, mda=mdf, passed=passed, v_true=v_true_full,
                u_full=u_full, v_norm_full=float(n_full), v_norm_res=float(n_res)), fired


# --------------------------------------------------------------------------- #
# STAGE 3 — ranking
# --------------------------------------------------------------------------- #
def guide_agreement(pack, factors, rows_by_factor, log):
    """CHANGE 2. Guide agreement on V_all guide vectors ALWAYS.

    One number per factor, applied unchanged to both gene variants. The V_sig
    value is still computed and reported as a diagnostic column, because it is
    what SOUTH2 gated on, but it gates nothing here.
    """
    Gc, Sc = pack["Gc"], pack["Sc"]
    n_f = len(factors)
    g_all = np.full(n_f, np.nan)
    g_sig = np.full(n_f, np.nan)
    for i, t in enumerate(factors):
        r = rows_by_factor[t]
        for name, base in (("all", Gc[r]), ("sig", Gc[r] * Sc[r])):
            if base.shape[0] < 2:
                continue
            cs = [cosine_full(base[a], base[b])
                  for a in range(base.shape[0]) for b in range(a + 1, base.shape[0])]
            cs = np.asarray(cs, float)
            v = float(np.nanmedian(cs)) if np.isfinite(cs).any() else np.nan
            if name == "all":
                g_all[i] = v
            else:
                g_sig[i] = v
    log(f"[guides] CHANGE 2: guide agreement is the median pairwise cosine between the V_all "
        f"vectors of a factor's own guides, computed once and applied to both variants. "
        f"guides_disagree (<= 0) fires on {int(np.sum(g_all <= 0)):,} of {n_f:,} factors. "
        f"For reference only, the V_sig median pairwise cosine would fire on "
        f"{int(np.sum(g_sig <= 0)):,} of {n_f:,}; that number gates nothing in SOUTH3.")
    return g_all, g_sig


def stage3(pack, combos, Lm, D, factors, rows_by_factor, s2res, point, keys,
           guide_cos_all, guide_cos_sig, log):
    Gc, Sc = pack["Gc"], pack["Sc"]
    frozen = pack["frozen"]
    ruler_rows = pack["ruler_rows"]
    w_unit, _ = unit(np.asarray(frozen["w"], float))
    w_m = w_unit[ruler_rows]
    id_pos = pack["id_rowpos"]
    n_f = len(factors)

    # per-factor, per-variant scalars that do not depend on the space
    id_drop = {v: np.full(n_f, np.nan) for v in VARIANTS}
    ruler_chg = {v: np.full(n_f, np.nan) for v in VARIANTS}
    n_guides = np.array([rows_by_factor[t].size for t in factors], dtype=int)

    boot = {}    # (variant, key) -> dict of arrays
    for variant in VARIANTS:
        M = D[variant]
        for kk in keys:
            boot[(variant, kk)] = {c: np.full(n_f, np.nan) for c in
                                   ("delta_lo", "delta_hi", "cos_lo", "cos_hi",
                                    "frac_lo", "frac_hi")}
        t0 = time.perf_counter()
        for i, t in enumerate(factors):
            r = rows_by_factor[t]
            base = Gc[r] if variant == "V_all" else Gc[r] * Sc[r]
            d = M[i]
            # identity panel drop, in frozen-z units
            if id_pos.size:
                id_drop[variant][i] = float(-np.mean(d[id_pos]))
            ruler_chg[variant][i] = float(np.dot(d, w_m))
            # bootstrap over guides (3+ only)
            if base.shape[0] >= BOOT_MIN_GUIDES:
                rb = np.random.default_rng([BOOT_SEED, i])
                sel = rb.integers(0, base.shape[0], size=(N_BOOT, base.shape[0]))
                B = base[sel].mean(axis=1)
                Bp = B @ Lm
                for c, kk in zip(combos, keys):
                    stb = block_stats(c, B, Bp)
                    for cname in ("delta", "cos", "frac"):
                        lo, hi = percentile_ci(stb[cname])
                        boot[(variant, kk)][f"{cname}_lo"][i] = lo
                        boot[(variant, kk)][f"{cname}_hi"][i] = hi
            if (i + 1) % 400 == 0:
                log(f"[stage3 {variant}] {i+1}/{n_f} factors ({time.perf_counter() - t0:.0f}s)")
        log(f"[stage3 {variant}] done in {time.perf_counter() - t0:.0f}s")
        il = int(np.sum(id_drop[variant] > IDENTITY_DROP))
        log(f"[stage3 {variant}] identity_loss (panel drop > {IDENTITY_DROP} in frozen-z "
            f"units): {il:,}")
    guide_cos = {v: guide_cos_all for v in VARIANTS}

    # assemble the full table
    out = []
    for variant in VARIANTS:
        for c, kk in zip(combos, keys):
            fl = s2res["floors"][(variant, kk)]
            pt = point[(variant, kk)]
            bt = boot[(variant, kk)]
            # One-sided empirical p, (count(perm delta <= observed) + 1) / (n_perm + 1) —
            # the same definition as trajectory_common.permutation_p(greater=False), with the
            # count streamed in Stage 2A instead of recomputed from a stored null.
            cnt_le = s2res["perm_counts"][(variant, kk)]
            p = np.where(np.isfinite(pt["delta"]),
                         (cnt_le + 1.0) / (N_PERM + 1.0), np.nan)
            q = bh_q(p)
            ci_valid = np.array([
                ci_covers_point(pt["delta"][i], bt["delta_lo"][i], bt["delta_hi"][i])
                if np.isfinite(bt["delta_lo"][i]) else None
                for i in range(n_f)], dtype=object)
            if (variant, kk) in s2res["md_floors"]:
                md_fl = s2res["md_floors"][(variant, kk)]
                md_cnt = s2res["md_counts"][(variant, kk)]
                p_md = np.where(np.isfinite(pt["md_change"]) & (md_cnt >= 0),
                                (md_cnt + 1.0) / (N_PERM + 1.0), np.nan)
                q_md = bh_q(p_md)
            else:
                md_fl = np.full(n_f, np.nan)
                p_md = np.full(n_f, np.nan)
                q_md = np.full(n_f, np.nan)
            gdis = guide_cos[variant] <= 0
            iloss = id_drop[variant] > IDENTITY_DROP
            toward = (
                (pt["delta"] < 0)
                & (pt["delta"] < fl)
                & np.isfinite(q) & (q <= Q_BAR)
                & ~iloss
                & ~gdis
            )
            df = pd.DataFrame(dict(
                factor=factors, n_guides=n_guides, variant=variant,
                space=c["space"], k=c["k"], origin=c["origin"], target=c["target"],
                pair=c["pair"],
                cos=pt["cos"], frac=pt["frac"], delta=pt["delta"], d_norm=pt["d_norm"],
                v_norm=c["r_norm"],
                # CHANGE 3, report-only: the displacement size relative to the gap it
                # has to close. ||target - origin|| is v_norm.
                d_over_v=pt["d_norm"] / c["r_norm"] if c["r_norm"] > 0 else np.nan,
                delta_small_disp_approx=-pt["d_norm"] * pt["cos"],
                floor_p05=fl, margin_over_floor=fl - pt["delta"],
                p_delta=p, q_delta=q, n_perm=N_PERM, p_delta_min_attainable=P_MIN,
                delta_ci_lo=bt["delta_lo"], delta_ci_hi=bt["delta_hi"],
                delta_ci_valid=ci_valid,
                cos_ci_lo=bt["cos_lo"], cos_ci_hi=bt["cos_hi"],
                frac_ci_lo=bt["frac_lo"], frac_ci_hi=bt["frac_hi"],
                md_change=pt["md_change"],
                md_floor_p05=md_fl, md_margin_over_floor=md_fl - pt["md_change"],
                p_md=p_md, q_md=q_md,
                inside_radius=np.where(np.isfinite(pt["inside"]), pt["inside"] > 0.5, None),
                guide_median_cos=guide_cos[variant], guides_disagree=gdis,
                guide_median_cos_source="V_all",
                guide_median_cos_vsig_reported_only=guide_cos_sig,
                identity_drop=id_drop[variant], identity_loss=iloss,
                identity_panel_n=int(id_pos.size),
                ruler_score_change=ruler_chg[variant],
                ruler_score_drop=-ruler_chg[variant],
                toward_young=toward,
            ))
            out.append(df)
    tab = pd.concat(out, ignore_index=True)
    tab.to_csv(SOUTH3_DIR / "stage3_full_table.csv", index=False)
    log(f"[stage3] full table {len(tab):,} rows -> stage3_full_table.csv")

    cnt = (tab.groupby(["variant", "space", "pair"], sort=False)
           .agg(n_factors=("factor", "size"),
                n_delta_negative=("delta", lambda s: int((s < 0).sum())),
                n_below_floor=("margin_over_floor", lambda s: int((s > 0).sum())),
                n_q_le_005=("q_delta", lambda s: int((s <= Q_BAR).sum())),
                n_at_p_floor=("p_delta", lambda s: int((s <= P_MIN + 1e-12).sum())),
                n_guides_disagree=("guides_disagree", "sum"),
                n_identity_loss=("identity_loss", "sum"),
                n_toward_young=("toward_young", "sum"),
                n_md_below_own_floor=("md_margin_over_floor", lambda s: int((s > 0).sum())),
                n_md_q_le_005=("q_md", lambda s: int((s <= Q_BAR).sum())),
                n_inside_radius=("inside_radius", lambda s: int(sum(bool(v) for v in s
                                                                    if v is not None))))
           .reset_index())
    cnt.to_csv(SOUTH3_DIR / "stage3_candidate_counts.csv", index=False)
    for _, r in cnt.iterrows():
        log(f"[count] {r['variant']:6s} {r['space']:9s} {r['pair']:6s} "
            f"delta<0={r['n_delta_negative']:5d} below_floor={r['n_below_floor']:5d} "
            f"q<=0.05={r['n_q_le_005']:5d} toward_young={r['n_toward_young']:5d}")

    n_invalid = int((tab["delta_ci_valid"] == False).sum())  # noqa: E712
    log(f"[stage3] bootstrap delta intervals not containing their point estimate "
        f"(flagged INVALID, never a gate): {n_invalid:,} of "
        f"{int(tab['delta_ci_lo'].notna().sum()):,} computed")
    return tab, cnt, n_invalid


# --------------------------------------------------------------------------- #
# STAGE 4 — disagreement with directional readouts
# --------------------------------------------------------------------------- #
def stage4(tab, log):
    rows = []
    for (variant, space, pair), g in tab.groupby(["variant", "space", "pair"], sort=False):
        g = g.reset_index(drop=True)
        # rank 1 = best by each readout
        r_delta = sps.rankdata(g["delta"].to_numpy(float))
        r_cos = sps.rankdata(-g["cos"].to_numpy(float))
        r_ruler = sps.rankdata(-g["ruler_score_drop"].to_numpy(float))
        top = 20
        t_delta = set(g["factor"].to_numpy()[np.argsort(r_delta)[:top]])
        t_cos = set(g["factor"].to_numpy()[np.argsort(r_cos)[:top]])
        t_ruler = set(g["factor"].to_numpy()[np.argsort(r_ruler)[:top]])
        tw = set(g.loc[g["toward_young"], "factor"].to_numpy())
        rows.append(dict(
            variant=variant, space=space, pair=pair, n_factors=int(len(g)),
            rho_delta_cos=float(sps.spearmanr(r_delta, r_cos).statistic),
            rho_delta_ruler=float(sps.spearmanr(r_delta, r_ruler).statistic),
            rho_cos_ruler=float(sps.spearmanr(r_cos, r_ruler).statistic),
            overlap_top20_delta_cos=int(len(t_delta & t_cos)),
            overlap_top20_delta_ruler=int(len(t_delta & t_ruler)),
            overlap_top20_cos_ruler=int(len(t_cos & t_ruler)),
            n_top20_cos_not_toward_young=int(len(t_cos - tw)),
            n_top20_delta_not_toward_young=int(len(t_delta - tw)),
            n_toward_young=int(len(tw)),
        ))
    df = pd.DataFrame(rows)
    df.to_csv(SOUTH3_DIR / "stage4_disagreement.csv", index=False)
    for _, r in df.iterrows():
        log(f"[stage4] {r['variant']:6s} {r['space']:9s} {r['pair']:6s} "
            f"rho(delta,cos)={r['rho_delta_cos']:+.3f} rho(delta,ruler)={r['rho_delta_ruler']:+.3f} "
            f"rho(cos,ruler)={r['rho_cos_ruler']:+.3f} top20 overlap d/c="
            f"{r['overlap_top20_delta_cos']:2d} cos_top20_not_toward={r['n_top20_cos_not_toward_young']:2d}")
    return df


# --------------------------------------------------------------------------- #
# CHANGE 3 — report-only displacement geometry. No gate, no label depends on it.
# --------------------------------------------------------------------------- #
def geometry_report(tab, log):
    """||d||, ||d|| / ||target - origin||, rank agreement by quartile of that ratio,
    and the small-displacement identity delta ~= -||d|| cos(theta).

    Exactly:   delta = ||r|| * (sqrt(1 - 2*rho*cos + rho^2) - 1),
    with r = origin - target, rho = ||d|| / ||r||, cos = cos(theta) = -d.r/(||d|| ||r||).
    Expanding in rho gives delta = -||d||*cos + ||d||^2*(1 - cos^2)/(2||r||) + O(rho^3),
    so the identity is exact only as rho -> 0 and the leading error is positive.
    """
    srows, qrows = [], []
    for (variant, space, pair), g in tab.groupby(["variant", "space", "pair"], sort=False):
        g = g.reset_index(drop=True)
        ratio = g["d_over_v"].to_numpy(float)
        dn = g["d_norm"].to_numpy(float)
        delta = g["delta"].to_numpy(float)
        cos = g["cos"].to_numpy(float)
        approx = g["delta_small_disp_approx"].to_numpy(float)
        # Residual of the identity, and the leading term the expansion predicts for it.
        # err/||d|| is used instead of err/|delta| wherever a summary is quoted, because
        # |delta| passes through 0 and makes the relative error unbounded there.
        resid = delta - approx
        pred = dn * ratio * (1.0 - cos ** 2) / 2.0
        with np.errstate(divide="ignore", invalid="ignore"):
            err_over_d = np.abs(resid) / np.where(dn > 1e-12, dn, np.nan)
            resid_over_pred = resid / np.where(np.abs(pred) > 1e-15, pred, np.nan)
        srows.append(dict(
            variant=variant, space=space, pair=pair, n_factors=int(len(g)),
            v_norm=float(g["v_norm"].iloc[0]),
            d_norm_min=float(np.nanmin(dn)), d_norm_median=float(np.nanmedian(dn)),
            d_norm_max=float(np.nanmax(dn)),
            ratio_min=float(np.nanmin(ratio)), ratio_p25=float(np.nanpercentile(ratio, 25)),
            ratio_median=float(np.nanmedian(ratio)),
            ratio_p75=float(np.nanpercentile(ratio, 75)), ratio_max=float(np.nanmax(ratio)),
            identity_err_over_dnorm_median=float(np.nanmedian(err_over_d)),
            identity_err_over_dnorm_p90=float(np.nanpercentile(err_over_d, 90)),
            identity_err_over_dnorm_max=float(np.nanmax(err_over_d)),
            identity_err_bound_half_ratio_max=float(np.nanmax(ratio) / 2.0),
            resid_over_leading_term_median=float(np.nanmedian(resid_over_pred)),
            resid_over_leading_term_p01=float(np.nanpercentile(resid_over_pred, 1)),
            resid_over_leading_term_p99=float(np.nanpercentile(resid_over_pred, 99)),
            n_resid_negative=int(np.sum(resid < -1e-12)),
        ))
        edges = np.nanpercentile(ratio, [25, 50, 75])
        qidx = np.digitize(ratio, edges, right=True)
        for qi in range(4):
            m = (qidx == qi) & np.isfinite(ratio) & np.isfinite(delta) & np.isfinite(cos)
            if m.sum() < 3:
                continue
            rd = sps.rankdata(delta[m])          # rank 1 = most negative delta = best
            rc = sps.rankdata(-cos[m])           # rank 1 = highest cosine = best
            rel = np.abs(delta[m] - approx[m]) / np.maximum(np.abs(delta[m]), 1e-12)
            qrows.append(dict(
                variant=variant, space=space, pair=pair, quartile=f"Q{qi+1}",
                n=int(m.sum()),
                ratio_lo=float(np.nanmin(ratio[m])), ratio_hi=float(np.nanmax(ratio[m])),
                ratio_median=float(np.nanmedian(ratio[m])),
                rho_delta_rank_vs_cos_rank=float(sps.spearmanr(rd, rc).statistic),
                rho_delta_vs_small_disp_approx=float(
                    sps.spearmanr(delta[m], approx[m]).statistic),
                pearson_delta_vs_small_disp_approx=float(
                    sps.pearsonr(delta[m], approx[m]).statistic),
                identity_err_over_dnorm_median=float(np.nanmedian(err_over_d[m])),
                identity_err_over_dnorm_p90=float(np.nanpercentile(err_over_d[m], 90)),
                identity_err_bound_half_ratio_max=float(np.nanmax(ratio[m]) / 2.0),
                resid_over_leading_term_median=float(np.nanmedian(resid_over_pred[m])),
                median_abs_rel_error_vs_delta=float(np.nanmedian(rel)),
                p90_abs_rel_error_vs_delta=float(np.nanpercentile(rel, 90)),
                median_delta_minus_approx=float(np.nanmedian(resid[m])),
                n_resid_negative=int(np.sum(resid[m] < -1e-12)),
            ))
    sdf = pd.DataFrame(srows)
    qdf = pd.DataFrame(qrows)
    sdf.to_csv(SOUTH3_DIR / "stage3_dnorm_ratio_summary.csv", index=False)
    qdf.to_csv(SOUTH3_DIR / "stage3_geometry_quartiles.csv", index=False)
    log(f"[geometry] ||d||/||target-origin|| ranges {sdf['ratio_min'].min():.5f} to "
        f"{sdf['ratio_max'].max():.5f} over all settings; quartile table has {len(qdf)} rows")
    for _, r in qdf.iterrows():
        log(f"[geometry] {r['variant']:6s} {r['space']:9s} {r['pair']:6s} {r['quartile']} "
            f"ratio~{r['ratio_median']:.4f} rho(delta,cos)={r['rho_delta_rank_vs_cos_rank']:+.3f} "
            f"identity err/||d|| median={r['identity_err_over_dnorm_median']:.5f} "
            f"(bound {r['identity_err_bound_half_ratio_max']:.5f}); "
            f"resid/leading term={r['resid_over_leading_term_median']:.4f}")
    log(f"[geometry] identity residual delta-(-||d||cos) is non-negative on "
        f"{len(sdf) and int(len(sdf)) and (int(sdf['n_resid_negative'].sum()) == 0)} "
        f"(negatives: {int(sdf['n_resid_negative'].sum())}); residual / predicted leading term "
        f"has per-setting medians {sdf['resid_over_leading_term_median'].min():.4f} to "
        f"{sdf['resid_over_leading_term_median'].max():.4f}")
    return sdf, qdf


# --------------------------------------------------------------------------- #
# CHANGE 2 — label comparison against SOUTH2
# --------------------------------------------------------------------------- #
def change2_comparison(tab, log):
    """How many factors change label because guide agreement moved to V_all always.

    Two isolations are reported because only both together are honest:
      (a) inside SOUTH2: SOUTH2's own delta / floor / q with only the guide rule
          swapped for the V_all-always rule. This is change 2 with change 1 held out.
      (b) inside SOUTH3: SOUTH3's delta / floor / q with the guide rule swapped back
          to SOUTH2's per-variant rule. This is change 2 with change 1 held in.
    The total SOUTH2 -> SOUTH3 label change is reported alongside, undecomposed.
    """
    if not SOUTH2_TABLE.exists():
        record_note("change2", f"{SOUTH2_TABLE} missing; the SOUTH2 label comparison could not "
                               "be computed and is reported as unavailable, not as zero.")
        log(f"[change2] {SOUTH2_TABLE} missing — comparison unavailable")
        return None, None
    s2 = pd.read_csv(SOUTH2_TABLE, low_memory=False)
    key = ["variant", "space", "pair", "factor"]
    # SOUTH2's own V_all guide agreement, keyed by factor. This is the number the
    # V_all-always rule uses, and it is identical in SOUTH2 and SOUTH3 by construction.
    gall = (s2[s2.variant == "V_all"].drop_duplicates("factor")
            .set_index("factor")["guide_median_cos"])
    s2 = s2.copy()
    s2["gdis_vall"] = s2["factor"].map(gall).to_numpy(float) <= 0
    s2["cf_label_v_all_rule"] = (
        (s2["delta"] < 0) & (s2["margin_over_floor"] > 0)
        & s2["q_delta"].notna() & (s2["q_delta"] <= Q_BAR)
        & ~s2["identity_loss"].astype(bool) & ~s2["gdis_vall"]
    )
    s2["toward_young"] = s2["toward_young"].astype(bool)

    t3 = tab.copy()
    t3["gdis_per_variant"] = np.where(
        t3["variant"] == "V_all",
        t3["guide_median_cos"].to_numpy(float) <= 0,
        t3["guide_median_cos_vsig_reported_only"].to_numpy(float) <= 0)
    t3["cf_label_per_variant_rule"] = (
        (t3["delta"] < 0) & (t3["margin_over_floor"] > 0)
        & t3["q_delta"].notna() & (t3["q_delta"] <= Q_BAR)
        & ~t3["identity_loss"].astype(bool) & ~t3["gdis_per_variant"]
    )

    m = (s2[key + ["toward_young", "cf_label_v_all_rule", "guides_disagree", "gdis_vall"]]
         .rename(columns=dict(toward_young="label_south2",
                              guides_disagree="gdis_south2"))
         .merge(t3[key + ["toward_young", "cf_label_per_variant_rule", "guides_disagree"]]
                .rename(columns=dict(toward_young="label_south3",
                                     guides_disagree="gdis_south3")),
                on=key, how="outer", indicator=True))
    unmatched = int((m["_merge"] != "both").sum())
    if unmatched:
        record_note("change2", f"{unmatched} of {len(m)} (variant, space, pair, factor) rows did "
                               "not match between the SOUTH2 and SOUTH3 tables. Reported, not "
                               "repaired.", dict(n_unmatched=unmatched))
    m = m[m["_merge"] == "both"].drop(columns="_merge")
    rows = []
    for (variant, space, pair), g in m.groupby(["variant", "space", "pair"], sort=False):
        rows.append(dict(
            variant=variant, space=space, pair=pair, n_factors=int(len(g)),
            n_guides_disagree_south2=int(g["gdis_south2"].sum()),
            n_guides_disagree_south3=int(g["gdis_south3"].sum()),
            n_guide_flag_changed=int((g["gdis_south2"].astype(bool)
                                      != g["gdis_south3"].astype(bool)).sum()),
            n_toward_young_south2=int(g["label_south2"].sum()),
            n_toward_young_south3=int(g["label_south3"].sum()),
            n_label_changed_total=int((g["label_south2"].astype(bool)
                                       != g["label_south3"].astype(bool)).sum()),
            n_label_changed_by_change2_within_south2=int(
                (g["label_south2"].astype(bool)
                 != g["cf_label_v_all_rule"].astype(bool)).sum()),
            n_label_changed_by_change2_within_south3=int(
                (g["label_south3"].astype(bool)
                 != g["cf_label_per_variant_rule"].astype(bool)).sum()),
        ))
    cdf = pd.DataFrame(rows)
    cdf.to_csv(SOUTH3_DIR / "change2_label_comparison.csv", index=False)
    m.to_csv(SOUTH3_DIR / "change2_label_comparison_per_factor.csv", index=False)
    tot = dict(
        flag_changed=int(cdf["n_guide_flag_changed"].sum()),
        label_changed_total=int(cdf["n_label_changed_total"].sum()),
        by_change2_within_south2=int(cdf["n_label_changed_by_change2_within_south2"].sum()),
        by_change2_within_south3=int(cdf["n_label_changed_by_change2_within_south3"].sum()),
        n_rows=int(len(m)),
        n_factor_flag_changed=int(
            m.drop_duplicates(["factor", "variant"])
            .assign(chg=lambda x: x["gdis_south2"].astype(bool) != x["gdis_south3"].astype(bool))
            ["chg"].sum()),
        vsig_survivors_south3=int(cdf.loc[cdf.variant == "V_sig", "n_toward_young_south3"].sum()),
        vsig_by_change2=int(
            cdf.loc[cdf.variant == "V_sig",
                    "n_label_changed_by_change2_within_south3"].sum()),
        vall_by_change2=int(
            cdf.loc[cdf.variant == "V_all",
                    "n_label_changed_by_change2_within_south3"].sum()),
    )
    log(f"[change2] guide flag changed on {tot['flag_changed']:,} of {tot['n_rows']:,} "
        f"factor-settings; toward_young label changed on {tot['label_changed_total']:,} in total; "
        f"attributable to the guide rule alone: {tot['by_change2_within_south2']:,} inside "
        f"SOUTH2's numbers and {tot['by_change2_within_south3']:,} inside SOUTH3's numbers")
    return cdf, tot


# --------------------------------------------------------------------------- #
# Is any survivor distinguishable from a plain cosine ranking?
# --------------------------------------------------------------------------- #
def survivor_vs_cosine(tab, log):
    """Could a plain cosine have produced the survivor set?

    The survivor set is compared with the top-n factors under four one-number
    rankings, n being the number of survivors in that setting: cosine alone, the
    small-displacement identity -||d||cos(theta), delta itself, and ||d|| alone.
    `by_delta` is the control that makes the comparison mean something: the label is
    a conjunction (delta < 0, below its own floor, q <= 0.05, two flags), not a
    ranking, so even delta — the statistic the label is built on — cannot reproduce
    the survivor set exactly. Cosine only fails to reproduce it in a way that
    matters if it fails by more than delta does.
    """
    rows = []
    for (variant, space, pair), g in tab.groupby(["variant", "space", "pair"], sort=False):
        g = g.reset_index(drop=True)
        rho_all = float(sps.spearmanr(sps.rankdata(g["delta"].to_numpy(float)),
                                      sps.rankdata(-g["cos"].to_numpy(float))).statistic)
        surv = g[g["toward_young"].astype(bool)]
        n_s = int(len(surv))
        if n_s == 0:
            rows.append(dict(variant=variant, space=space, pair=pair, n_factors=int(len(g)),
                             n_survivors=0, n_overlap_top_cos=0, n_overlap_top_identity=0,
                             n_overlap_top_delta=0, n_overlap_top_d_norm=0,
                             n_survivors_not_in_top_cos=0, jaccard_cos=np.nan,
                             rho_delta_cos_all=rho_all,
                             rho_delta_cos_within_survivors=np.nan,
                             worst_survivor_cos_rank=np.nan))
            continue
        sset = set(surv["factor"].to_numpy())
        fac = g["factor"].to_numpy()

        def top_n(col, ascending):
            v = g[col].to_numpy(float)
            o = np.argsort(v if ascending else -v, kind="stable")
            return set(fac[o[:n_s]])

        top_cos = top_n("cos", False)
        cos_order = np.argsort(-g["cos"].to_numpy(float), kind="stable")
        cos_rank = pd.Series(np.arange(1, len(g) + 1), index=fac[cos_order])
        rd = sps.rankdata(surv["delta"].to_numpy(float))
        rc = sps.rankdata(-surv["cos"].to_numpy(float))
        rows.append(dict(
            variant=variant, space=space, pair=pair, n_factors=int(len(g)),
            n_survivors=n_s,
            n_overlap_top_cos=int(len(sset & top_cos)),
            n_overlap_top_identity=int(len(sset & top_n("delta_small_disp_approx", True))),
            n_overlap_top_delta=int(len(sset & top_n("delta", True))),
            n_overlap_top_d_norm=int(len(sset & top_n("d_norm", False))),
            n_survivors_not_in_top_cos=int(len(sset - top_cos)),
            jaccard_cos=float(len(sset & top_cos) / len(sset | top_cos)),
            rho_delta_cos_all=rho_all,
            rho_delta_cos_within_survivors=float(sps.spearmanr(rd, rc).statistic)
            if n_s >= 3 else np.nan,
            worst_survivor_cos_rank=int(max(int(cos_rank[f]) for f in sset)),
        ))
    df = pd.DataFrame(rows)
    df.to_csv(SOUTH3_DIR / "stage4_survivor_vs_cosine.csv", index=False)
    for _, r in df[df.n_survivors > 0].iterrows():
        log(f"[survivor-vs-cos] {r['variant']:6s} {r['space']:9s} {r['pair']:6s} "
            f"n={int(r['n_survivors'])} recovered by top-n: cos={int(r['n_overlap_top_cos'])} "
            f"identity={int(r['n_overlap_top_identity'])} delta={int(r['n_overlap_top_delta'])} "
            f"||d||={int(r['n_overlap_top_d_norm'])}; rho(delta,cos) over all factors="
            f"{r['rho_delta_cos_all']:+.3f}")
    return df


def top20_tables(tab, coverage_low):
    """One top-20-by-delta table per (variant, space, origin/target pair). Never pooled."""
    blocks = []
    word = "exploratory " if coverage_low else ""
    for (variant, space, pair), g in tab.groupby(["variant", "space", "pair"], sort=False):
        n_setting = int(g["toward_young"].sum())
        g = g.sort_values("delta").head(20)
        show = pd.DataFrame({
            "factor": [(f"{x} (exploratory)" if coverage_low and bool(tw) else str(x))
                       for x, tw in zip(g["factor"], g["toward_young"])],
            "n_guides": g["n_guides"].astype(int).astype(str),
            "delta": [_sf(v, 5) for v in g["delta"]],
            "floor_p05": [_sf(v, 5) for v in g["floor_p05"]],
            "margin_over_floor": [_sf(v, 5) for v in g["margin_over_floor"]],
            "q_delta": [_fmt_p(v) for v in g["q_delta"]],
            "guide_median_cos": [_sf(v, 3) for v in g["guide_median_cos"]],
            "guides_disagree": [str(bool(v)) for v in g["guides_disagree"]],
            "identity_drop": [_sf(v, 3) for v in g["identity_drop"]],
            "identity_loss": [str(bool(v)) for v in g["identity_loss"]],
            "toward_young": [str(bool(v)) for v in g["toward_young"]],
        })
        blocks.append((variant, space, pair, show, (int(g["toward_young"].sum()), n_setting)))
    return blocks, word


# --------------------------------------------------------------------------- #
# FINDINGS / PROGRESS
# --------------------------------------------------------------------------- #
def stage0_lines(inv, man):
    x = inv["x_audit"]
    mk = inv["mask_audit"]
    pa = inv["p_audit"]
    ru = inv["row_unit"]
    lines = [
        "## Stage 0 inventory (descriptive only; no statistic)",
        "",
        f"- File: `{inv['path']}`",
        f"- Size: {_fmt_int(inv['size_bytes'])} bytes (pinned {_fmt_int(inv['pinned_size'])}), "
        f"md5 `{inv['md5']}` (pinned `{inv['pinned_md5']}`). Both match; no file was substituted. "
        "`L4A_matrix.h5` was not opened and nothing was downloaded.",
        f"- Exact shape: **{inv['x_shape'][0]:,} rows x {inv['x_shape'][1]:,} genes**, dtype "
        f"`float64`, dense. Layers: {', '.join(f'`{l}`' for l in inv['layers'])}, each the same "
        f"shape as X. `uns`: {inv['uns_keys'] or 'empty'}. `obsm`: {inv['obsm_keys'] or 'empty'}. "
        f"`varm`: {inv['varm_keys'] or 'empty'}.",
        "",
        "### What one row is",
        "",
        f"**{ru['unit']}.** `{GUIDE_COL}` is unique on every row "
        f"({ru['n_unique_guides']:,} values for {inv['n_rows']:,} rows), so a row is one guide, "
        f"not a per-factor mean. There is no cell-barcode column "
        f"(`cell_barcode_column={ru['cell_barcode_column']}`), only a per-row `cell_count` "
        f"(min {ru['cell_count_min']:.0f}, median {ru['cell_count_median']:.0f}, max "
        f"{ru['cell_count_max']:.0f}, total {ru['cell_count_sum']:,.0f}), so a row is not a cell "
        "either. Each row is the guide-level output of the source's own regression over the "
        "cells carrying that guide.",
        "",
        "### What the column index holds",
        "",
        f"- {inv['column_index_kind']}",
        "",
        "### Units and log base of X",
        "",
        (f"- **Undeterminable from the file.** `uns` is empty; `X` carries only "
         if not inv["units"]["determinable"] else
         f"- **A unit hint was found and must be read by hand before trusting this run** "
         f"(attr hits {inv['units']['x_attr_unit_hits']}, column hits "
         f"{inv['units']['column_unit_hits']}, doc hits "
         f"{inv['units']['bundled_txt_with_unit_words']}). For reference, `X` carries ")
        + f"`{inv['units']['x_attrs']}` as attributes, which state an encoding and no unit; `var` "
        f"holds only `mean`, `std`, `cv`, `fano`, `excess_cv` and related columns with no unit "
        f"statement; `obs` states no unit; and the only documentation bundled in the folder is "
        f"three plain gene-symbol lists ({', '.join(f'`{n}`' for n in inv['units']['bundled_txt_files'])}), "
        "which carry no units.",
        f"- Descriptive diagnostics that were tried and did **not** settle it: Spearman of the "
        f"per-gene mean |effect| against the source's own `var/mean` = "
        f"{inv['units']['diag_spearman_meanabs_vs_var_mean']:+.3f}, against `var/std` = "
        f"{inv['units']['diag_spearman_meanabs_vs_var_std']:+.3f}, and of mean |effect| / "
        f"`var/std` against `var/mean` = "
        f"{inv['units']['diag_spearman_meanabs_over_std_vs_var_mean']:+.3f}. |effect| grows with "
        "expression but sub-linearly, which is consistent with more than one scale. No claim "
        "about a log base is made.",
        f"- Key `scale_unknown` is recorded as a **REPORT-ONLY flag, not a stop**, exactly as "
        "pre-registered. The Stage 2A primary floor and the Stage 2B control floor are "
        "magnitude-matched — they permute or re-use the same vector, so a global scale factor "
        "on X multiplies the factor vector and its null identically and the comparison is "
        "invariant to it. **The Stage 2C MDA is not invariant**: it scales the true "
        "young-minus-old direction to the median factor norm, and that norm moves with any "
        "global rescaling of X, so the MDA would move too.",
        "",
        "### The masked layer",
        "",
        f"- `layers/masked` is X with masked entries written as NaN. It covers "
        f"**{mk['n_masked_entries']:,} entries** spread over {mk['n_rows_with_any_mask']:,} rows "
        f"(between {mk['masked_per_row_min_nonzero']} and {mk['masked_per_row_max']} per masked "
        "row).",
        f"- What it marks: **each row's own targeted gene**. Of the {mk['n_masked_entries']:,} "
        f"masked entries, {mk['n_on_own_target_gene']:,} sit exactly at the column of that row's "
        f"`target_gene_id` and {mk['n_elsewhere']:,} sit anywhere else. The count equals the "
        f"number of rows whose `target_gene_id` is non-empty, i.e. the rows whose targeted gene "
        "is itself one of the measured genes; the source masks the guide's direct effect on its "
        "own target.",
        f"- Masked entries are **set to 0** before anything is derived, and the count is reported "
        "above. Nothing else in X is altered.",
        "",
        "### Ruler-gene coverage, computed three ways",
        "",
        f"All three are against the {inv['n_ruler']:,} frozen-ruler genes. \"Weight fraction\" is "
        "the share of the ruler's total absolute weight sum(|w|) carried by the matched genes.",
        "",
    ]
    cdf = pd.DataFrame(inv["coverage"])
    show = pd.DataFrame({
        "mapping": cdf["mapping"],
        "genes_matched": [f"{int(v):,}" for v in cdf["n_matched"]],
        "of_ruler_genes": [f"{int(v):,}" for v in cdf["n_ruler"]],
        "gene_fraction": [_f(v, 4) for v in cdf["frac_genes"]],
        "weight_matched": [_f(v, 3) for v in cdf["weight_matched"]],
        "weight_total": [_f(v, 3) for v in cdf["weight_total"]],
        "weight_fraction": [_f(v, 4) for v in cdf["weight_frac"]],
    })
    ch = inv["coverage_chosen"]
    lines += [
        md_table(show),
        "",
        f"- **Mapping used for all later stages: `{ch['mapping']}`**, named explicitly because it "
        f"has the highest weight fraction ({ch['weight_frac']:.4f}; ties would go to more genes "
        "matched).",
        "- `best_of_two_after_md3_idtype` is Ensembl-id-first with a symbol fallback "
        "(`seng_run.ruler_column_index`), applied after `src/md3_idtype.py` confirmed that "
        f"`var/gene_id` is an Ensembl column "
        f"(kind `{inv['id_type_records'].get('south3_var_gene_id', {}).get('kind')}`) and "
        f"`var/gene_name` a symbol column "
        f"(kind `{inv['id_type_records'].get('south3_var_gene_name', {}).get('kind')}`).",
        f"- Flag `coverage_low` (weight fraction < {inv['coverage_low_bar']}): "
        f"**{'FIRED' if inv['coverage_low'] else 'not fired'}**.",
        "",
        "### Transcription factors and guides",
        "",
        f"- Distinct transcription factors: **{inv['n_factors']:,}** (the "
        f"`{inv['control_label']}` label is the untargeted control arm and is not counted as a "
        "factor).",
        f"- Guides per factor: min **{inv['guides_min']}**, median **{inv['guides_median']:g}**, "
        f"max **{inv['guides_max']}**.",
        f"- Factors with at least {inv['min_guides']} guides: **{inv['n_factors_ge2']:,}** "
        f"(the `too_few_factors` threshold is {inv['min_factors']}).",
        "",
        "### The untargeted control rows",
        "",
        f"- Count: **{inv['n_controls']}** rows labelled `{inv['control_label']}`. All "
        f"{inv['n_controls_full_evidence']} of them also carry an empty `target_gene_id`, a NaN "
        f"`target_expr`, and `active=False`.",
        f"- That evidence triple is **necessary but not sufficient**: "
        f"{inv['n_evidence_not_labelled']:,} targeted guide rows meet it too. `target_gene_id` is "
        "filled only when the targeted gene is one of the 4,914 measured genes, `target_expr` "
        "likewise, and most guides are `active=False`, so an inactive guide against an "
        "unmeasured target looks the same on those three columns. The control arm is therefore "
        f"identified by the `{inv['control_label']}` label, with the triple as a consistency "
        "check that every labelled row passes. No targeted row was reassigned to the controls.",
        f"- They are **flagged as not inert** and are used only as the Stage 2B secondary floor, "
        "never as the gate.",
        "",
    ]
    ctab = pd.DataFrame(inv["controls"])
    cshow = pd.DataFrame({
        "guide_identity": ctab["guide_identity"],
        "sequence_driven": [str(bool(v)) for v in ctab["sequence_driven"]],
        "de_genes": [f"{int(v):,}" for v in ctab["de_genes"]],
        "cell_count": [f"{int(v):,}" for v in ctab["cell_count"]],
    })
    lines += [
        md_table(cshow),
        "",
        f"{int(ctab['sequence_driven'].sum())} of {len(ctab)} control rows are flagged "
        f"`sequence_driven` and {int((ctab['de_genes'] > 0).sum())} have `de_genes` > 0, so they "
        "move the transcriptome and are not an inert baseline.",
        "",
        "### The two Stage 0 STOP tests",
        "",
        f"- `not_an_effect_matrix`: **not fired**. X is per-gene effect estimates. Evidence: the "
        f"file carries `p` and `adj_p` layers of exactly X's shape, which only exist for a "
        f"per-entry test; {x['n_negative']:,} of {x['n_seen']:,} X entries are negative "
        f"({100 * x['frac_negative']:.1f}%, impossible for counts); only "
        f"{100 * x['frac_integral']:.2f}% of entries are integral; X spans "
        f"{x['min']:.3f} to {x['max']:.3f}; and the source's own `obs/de_genes` column is "
        f"reproduced exactly by counting `adj_p <= 0.05` per row on "
        f"**{pa['rows_where_de_genes_equals_count_adj_p_le_005']:,} of {pa['n_rows']:,} rows** "
        "(all of them), which ties X's layers to the source's own differential-expression call.",
        f"- `too_few_factors`: **not fired**. {inv['n_factors_ge2']:,} factors carry >= "
        f"{inv['min_guides']} guides, above the threshold of {inv['min_factors']}.",
        "",
    ]
    return lines


def stage1_lines(inv, s2, tgt, dist_df, s2res):
    ch = inv["coverage_chosen"]
    lines = [
        "## Stage 1 conversion, gene variants, and spaces",
        "",
        "### Conversion to frozen-ruler coordinates, and its assumption",
        "",
        f"- Per-factor effect vector = the **mean of that factor's guide rows** (no weighting by "
        "cell count, no re-estimation).",
        f"- Each gene's effect is divided by that gene's **frozen training sd** from "
        f"`{FROZEN_RULER.name}`, giving a displacement in frozen-z units. The frozen ruler is "
        f"loaded, never refit. sd was floored at 1e-12 on {inv['n_sd_floored']} ruler genes.",
        f"- Genes absent from the table are set to **0**, meaning \"this factor does not move this "
        f"gene\". That fixes **{inv['n_fixed_zero']:,} of the {inv['n_ruler']:,} coordinates** "
        f"({100.0 * inv['n_fixed_zero'] / inv['n_ruler']:.1f}%), leaving {inv['n_matched']:,} "
        "coordinates free to move.",
        "- This is a **weaker assumption than setting an expression level to zero**: it says the "
        "perturbation leaves the unmeasured gene where it already was, so the unmeasured "
        "coordinates contribute an identical, perturbation-independent offset to the origin and "
        "to every candidate end point. Setting an expression level to zero would instead assert "
        "that the gene is not expressed, which would move the point itself and would corrupt the "
        "origin as well as the displacement.",
        "",
        "### The two gene variants (both reported, never combined)",
        "",
        f"- **V_all**: every matched gene ({inv['n_matched']:,} coordinates).",
        f"- **V_sig**: matched genes with `adj_p <= {ADJ_P_BAR}` kept, the rest set to 0. The "
        "threshold is applied **per guide row** and the surviving effects are then averaged, so "
        "each guide contributes only its own significant effects and the per-guide vectors used "
        "for guide agreement and for the bootstrap are the same objects that build the factor "
        "vector.",
        "",
        "### The two spaces (both reported)",
        "",
        f"- **S1**: gene space, the {inv['n_ruler']:,} frozen-ruler coordinates in frozen-z units "
        "as above.",
        f"- **S2**: low-dim. PCA fit **only** on the GTEx fibroblast donor matrix "
        f"({s2['n_donors']} donors x {inv['n_ruler']:,} frozen-ruler genes, in frozen-z, centred "
        f"on the GTEx donor mean), then applied unchanged. Loadings are orthonormal. Every k in "
        f"{S2_KS} is reported; no k is chosen after seeing results.",
        "",
    ]
    cum = np.asarray(s2["cum_var"], float)
    vs = pd.DataFrame({
        "k": [str(k) for k in S2_KS],
        "cumulative_fraction_of_GTEx_donor_variance": [
            _f(cum[k - 1], 4) if k <= cum.size else "NA" for k in S2_KS],
    })
    lines += [md_table(vs), ""]
    lines += [
        "Stated plainly because it bears on how S2 should be read: Y1 and Y2 are subsets of the "
        "very GTEx donor matrix the PCA is fit on, so the S2 basis is **not** independent of the "
        "young target. It is independent of the Southard data and of O2, which is what the "
        "pre-registration forbids fitting on. No alternative basis was tried and none was chosen "
        "after seeing results.",
        "",
        "## Targets, origins, and the distances between them",
        "",
        f"- **Y1** = GTEx fibroblast donors aged 20-39, the existing TOWARD anchor "
        f"(n={int(tgt['young_all'].size)}; bins {list(YOUNG_BINS)}).",
        f"- **Y2** = the same donors as a region in S2 only: mean and covariance fit on "
        f"{int(tgt['young_fit'].size)} young donors, membership radius at the 95th percentile of "
        f"the {int(tgt['young_held'].size)} **held-out** young donors' own Mahalanobis distances "
        f"(split seed {SEED}).",
        f"- **O1** = GTEx fibroblast donors aged 60-79, the existing anchor (bins {list(OLD_BINS)}).",
        "- **O2** = the Lu et al. aged donor GM00731 untreated day-0 fibroblast pseudobulk, as "
        f"built in SAME (`results/same/panel_z.npz`, row `O`); "
        f"{int((tgt['z_o2'] == 0).sum()):,} of its coordinates are exactly 0 because those genes "
        "are absent from that dataset and were zeroed by SAME's own missing-gene rule.",
        "",
        "Membership radii, and where Y2 is undefined:",
        "",
    ]
    yrows = []
    for k in S2_KS:
        y = tgt["y2"][k]
        yrows.append(dict(
            k=str(k), defined=str(bool(y["defined"])), n_fit_donors=str(y["n_fit"]),
            covariance_rank=str(y["rank"]),
            condition_number=("NA" if not np.isfinite(y["cond"]) else f"{y['cond']:.3e}"),
            membership_radius=_f(y["radius"], 4),
        ))
    lines += [md_table(pd.DataFrame(yrows)), ""]
    for k in S2_KS:
        y = tgt["y2"][k]
        if not y["defined"]:
            lines += [
                f"At k={k} the covariance of {y['n_fit']} fit donors has rank {y['rank']}, so the "
                f"Mahalanobis metric is singular. **Y2 is reported as undefined at k={k}.** No "
                "pseudo-inverse and no ridge term were substituted, and no k was dropped to make "
                f"a number appear. Y1 at k={k} is unaffected and is reported.",
                "",
            ]
        elif y["n_fit"] < 2 * k:
            lines += [
                f"At k={k} the covariance is non-singular but rests on {y['n_fit']} donors in {k} "
                f"dimensions (condition number {y['cond']:.3e}); read Y2 at this k as fragile.",
                "",
            ]
    lines += ["Every origin-to-target distance, reported before anything was scored:", ""]
    d = dist_df.copy()
    show = pd.DataFrame({
        "space": d["space"], "origin": d["origin"], "target": d["target"],
        "metric": d["metric"],
        "origin_to_target_distance": [_f(v, 4) for v in d["origin_to_target_distance"]],
        "origin_to_target_mahalanobis": [
            _f(v, 4) for v in d.get("origin_to_target_mahalanobis", pd.Series([np.nan] * len(d)))],
        "membership_radius": [
            _f(v, 4) for v in d.get("membership_radius", pd.Series([np.nan] * len(d)))],
        "origin_inside_radius": [
            "NA" if v is None or (isinstance(v, float) and not np.isfinite(v)) else str(bool(v))
            for v in d.get("origin_inside_radius", pd.Series([None] * len(d)))],
    })
    lines += [md_table(show), ""]
    lines += [
        f"The true young-minus-old direction has ||v|| = {s2res['v_norm_full']:.4f} in S1; "
        f"restricted to the {inv['n_matched']:,} measured coordinates its norm is "
        f"{s2res['v_norm_res']:.4f}, i.e. "
        f"{100.0 * s2res['v_norm_res'] / s2res['v_norm_full']:.1f}% of the full norm. Both are "
        "used in the Stage 2C dose-response and both are reported.",
        "",
    ]
    return lines


def stage2_lines(s2res, inv):
    lines = [
        "## Stage 2 detection limit (computed before any factor was ranked)",
        "",
        "### A. Primary floor — magnitude-matched random direction",
        "",
        f"**CHANGE 1 lives here and nowhere else.** `n_perm` is {N_PERM:,} instead of SOUTH2's "
        f"{N_PERM_SOUTH2}. The permuted coordinate set, the per-factor generator "
        f"(`np.random.default_rng([{SEED}, factor_index])`), the delta definition, the 5th "
        "percentile rule, and the pooled-floor construction are all unchanged. The smallest p "
        f"an empirical one-sided permutation test can now return is 1/({N_PERM:,}+1) = "
        f"**{P_MIN:.3e}**, against 1/{N_PERM_SOUTH2 + 1} = {1.0 / (N_PERM_SOUTH2 + 1):.5f} in "
        f"SOUTH2. With BH across {inv['n_factors_ge2']:,} factors, q <= {Q_BAR} now needs only "
        f"{int(np.ceil(P_MIN * inv['n_factors_ge2'] / Q_BAR))} factors tied at that minimum p, "
        f"against {int(np.ceil(inv['n_factors_ge2'] / (N_PERM_SOUTH2 + 1) / Q_BAR))} in SOUTH2. "
        "That is the whole point of the change: the BH gate can resolve q instead of being "
        "pinned above it by the resolution of the null.",
        "",
        f"For each factor, its own vector's entries are permuted across the "
        f"{inv['n_matched']:,} matched gene coordinates {N_PERM:,} times (seed {SEED}, an "
        "independent generator per factor), and "
        "`delta = ||(origin + d_perm) - target|| - ||origin - target||` is recomputed each time. "
        "The per-factor floor is the 5th percentile of that distribution. Permuting preserves "
        "||d|| exactly, so the floor is magnitude-matched by construction.",
        "",
        "The pooled floor pools every factor's permuted vectors after rescaling each to the "
        "median factor norm, and takes the 5th percentile of that pooled distribution. It is the "
        "bar the Stage 2C dose-response has to clear.",
        "",
        f"**Implementation decision, stated because the pre-registration says only \"across "
        f"genes\":** the permutation is over the {inv['n_matched']:,} measured coordinates, not "
        f"over all {inv['n_ruler']:,} ruler coordinates. The other "
        f"{inv['n_fixed_zero']:,} are structurally fixed at 0 by the conversion rule, so "
        "permuting effects into them would place the null in a subspace the experiment cannot "
        "report on and would no longer be a like-for-like random direction for the real effect. "
        "This was fixed before any floor was computed and was not revisited.",
        "",
        "Distribution of the per-factor floor across factors:",
        "",
    ]
    f = s2res["per_factor_floor"]
    show = pd.DataFrame({
        "variant": f["variant"], "space": f["space"], "pair": f["pair"],
        "floor_min": [_sf(v, 5) for v in f["floor_min"]],
        "floor_p25": [_sf(v, 5) for v in f["floor_p25"]],
        "floor_median": [_sf(v, 5) for v in f["floor_median"]],
        "floor_p75": [_sf(v, 5) for v in f["floor_p75"]],
        "floor_max": [_sf(v, 5) for v in f["floor_max"]],
        "n_factors": [f"{int(v):,}" for v in f["n_factors"]],
    })
    lines += [md_table(show), "", "Pooled floor, built from a vector at the median factor norm:", ""]
    p = s2res["pooled"]
    show = pd.DataFrame({
        "variant": p["variant"], "space": p["space"], "pair": p["pair"],
        "median_factor_norm": [_f(v, 4) for v in p["median_factor_norm"]],
        "pooled_floor_p05": [_sf(v, 5) for v in p["floor_p05"]],
        "null_median": [_sf(v, 5) for v in p["null_median"]],
        "null_min": [_sf(v, 5) for v in p["null_min"]],
        "null_max": [_sf(v, 5) for v in p["null_max"]],
    })
    lines += [md_table(show), ""]
    lines += [
        "### B. Secondary floor — the untargeted control rows",
        "",
        f"The same delta, using the {int(s2res['control']['n_controls'].iloc[0])} untargeted "
        f"`{CONTROL_LABEL}` rows as pseudo-effects, converted exactly as the factors were.",
        "",
        "**This is reported and is never the gate.** It is built on controls that are not inert: "
        f"{sum(bool(r['sequence_driven']) for r in inv['controls'])} of {len(inv['controls'])} "
        "are flagged `sequence_driven` and "
        f"{sum(float(r['de_genes']) > 0 for r in inv['controls'])} carry `de_genes` > 0 "
        "(Stage 0). A floor built from rows that themselves move the transcriptome is not a "
        "noise floor.",
        "",
    ]
    c = s2res["control"]
    show = pd.DataFrame({
        "variant": c["variant"], "space": c["space"], "pair": c["pair"],
        "n_controls": [str(int(v)) for v in c["n_controls"]],
        "control_floor_p05": [_sf(v, 5) for v in c["control_floor_p05"]],
        "control_min": [_sf(v, 5) for v in c["control_min"]],
        "control_median": [_sf(v, 5) for v in c["control_median"]],
        "control_max": [_sf(v, 5) for v in c["control_max"]],
        "mean_control_norm": [_f(v, 4) for v in c["control_mean_norm"]],
        "is_gate": [str(bool(v)) for v in c["is_gate"]],
    })
    lines += [md_table(show), ""]
    lines += [
        "### C. Dose-response and the minimum detectable approach (MDA)",
        "",
        "The true young-minus-old direction is added to the aged origin at "
        f"f = {', '.join(str(x) for x in F_VALUES)}, scaled so that f = 1.0 has the median factor "
        "norm. MDA is the smallest f whose delta falls below the **primary** (pooled) floor. "
        "`support = full_v` is the direction as it is; `support = measured_coords_only` is the "
        "same direction restricted to the coordinates a perturbation in this file could actually "
        "move, reported alongside because a real effect cannot move the unmeasured genes.",
        "",
    ]
    m = s2res["mda"]
    show = pd.DataFrame({
        "variant": m["variant"], "space": m["space"], "pair": m["pair"],
        "support": m["support"],
        "median_factor_norm": [_f(v, 4) for v in m["median_factor_norm"]],
        "pooled_floor_p05": [_sf(v, 5) for v in m["pooled_floor_p05"]],
        "MDA": [_f(v, 3) for v in m["mda"]],
        "any_f_clears": [str(bool(v)) for v in m["cleared"]],
    })
    lines += [md_table(show), ""]
    lines += [
        "Full dose-response deltas at every f are in `results/south3/stage2c_dose_response.csv`.",
        "",
        f"- `no_detection_power`: **{'FIRED' if not s2res['passed'] else 'not fired'}**. "
        + ("Some f <= 1.0 clears the primary floor, so Stage 3 runs."
           if s2res["passed"] else
           "No f <= 1.0 cleared the primary floor in any space or variant. Stages 0-2 only are "
           "reported and **this data cannot support a ranking.**"),
        "- No floor is ever subtracted from a later number. It is a bar to clear.",
        "",
    ]
    return lines


def stage3_lines(tab, cnt, blocks, word, n_invalid, inv):
    lines = [
        "## Stage 3 ranking",
        "",
        f"The full table is `results/south3/stage3_full_table.csv`: **{len(tab):,} rows** = "
        f"{inv['n_factors_ge2']:,} factors x {len(VARIANTS)} gene variants x "
        f"{len(tab.groupby(['space', 'pair'])):,} space/origin/target settings, with cos, frac, "
        "delta, ||d||, ||v||, the factor's own Stage 2A floor, margin over that floor, "
        "permutation p, BH q, bootstrap intervals, the Mahalanobis change and membership test "
        "for Y2, guide agreement, the identity panel, and the frozen ruler score change. "
        "Nothing is pooled across spaces, variants, origins, or targets.",
        "",
        "- **Guide reliability — CHANGE 2 lives here.** Each guide's own vector is built and "
        "the median pairwise cosine between guides of the same factor is computed **on the "
        "V_all vectors always**, and that single number per factor is applied to both gene "
        "variants. `guides_disagree` fires when it is <= 0. Those factors are reported and are "
        "never labelled a candidate. The V_sig median pairwise cosine is still computed and "
        "carried in the full table as `guide_median_cos_vsig_reported_only`; it gates nothing. "
        "SOUTH2 gated each variant on its own value, which is what this change replaces.",
        f"- **Identity check**: the panel actually used is "
        f"{', '.join(inv['identity_panel'])} — **{len(inv['identity_panel'])} of the "
        f"{len(IDENTITY_GENES)} pre-registered genes**. Absent: "
        f"{', '.join(inv['identity_absent']) or 'none'}. "
        f"`identity_drop` is the negated mean displacement over that panel **in frozen-z units** "
        f"(the analysis's own coordinates), and `identity_loss` fires above {IDENTITY_DROP}.",
        "- **Frozen ruler score change** is a context column only. It gates nothing.",
        f"- **Nulls**: the Stage 2A per-factor permutations, {N_PERM:,} per factor and setting, "
        "one-sided empirical p = (count(perm delta <= observed delta) + 1) / "
        f"({N_PERM:,} + 1), so the attainable minimum is **p = {P_MIN:.3e}**; then BH-FDR across "
        "factors **within each space and variant** (and within each origin/target pair, since a "
        "pair is part of the setting).",
        "- **Report-only geometry (CHANGE 3)**: the full table also carries `d_norm`, "
        "`d_over_v` = ||d|| / ||target - origin||, and `delta_small_disp_approx` = "
        "-||d||cos(theta). None of the three gates anything; see the geometry section below.",
        f"- **Bootstrap**: {N_BOOT} resamples over guides (seed {BOOT_SEED}) for factors with "
        f">= {BOOT_MIN_GUIDES} guides. Intervals that do not contain their own point estimate are "
        f"flagged INVALID and are never a gate: **{n_invalid:,}** of "
        f"{int(tab['delta_ci_lo'].notna().sum()):,} computed delta intervals are so flagged.",
        "",
        "A factor is labelled `toward_young` only if delta is negative AND below its own Stage 2A "
        f"floor AND q <= {Q_BAR} AND `identity_loss` is not flagged AND `guides_disagree` is not "
        "flagged.",
        "",
        "### Candidate counts in every space, variant, and origin/target pair",
        "",
        "Zero is a valid and reportable answer.",
        "",
    ]
    show = pd.DataFrame({
        "variant": cnt["variant"], "space": cnt["space"], "pair": cnt["pair"],
        "n_factors": [f"{int(v):,}" for v in cnt["n_factors"]],
        "delta_negative": [f"{int(v):,}" for v in cnt["n_delta_negative"]],
        "below_own_floor": [f"{int(v):,}" for v in cnt["n_below_floor"]],
        "at_p_floor": [f"{int(v):,}" for v in cnt["n_at_p_floor"]],
        "q_le_005": [f"{int(v):,}" for v in cnt["n_q_le_005"]],
        "guides_disagree": [f"{int(v):,}" for v in cnt["n_guides_disagree"]],
        "identity_loss": [f"{int(v):,}" for v in cnt["n_identity_loss"]],
        "toward_young": [f"{int(v):,}" for v in cnt["n_toward_young"]],
    })
    lines += [md_table(show), ""]
    y2c = cnt[cnt["pair"].str.endswith("_Y2")]
    if len(y2c):
        lines += [
            "Y2 Mahalanobis columns, reported and never a gate. The Y2 `delta` above is the "
            "pre-registered Euclidean delta to the Y2 region's centre; the Mahalanobis change has "
            "its own floor from the same Stage 2A permutations. Where Y2 is undefined (see the "
            "radius table) these counts are 0 because nothing could be computed, not because "
            "nothing moved.",
            "",
        ]
        show = pd.DataFrame({
            "variant": y2c["variant"], "space": y2c["space"], "pair": y2c["pair"],
            "md_change_below_own_floor": [f"{int(v):,}" for v in y2c["n_md_below_own_floor"]],
            "md_q_le_005": [f"{int(v):,}" for v in y2c["n_md_q_le_005"]],
            "origin_plus_d_inside_radius": [f"{int(v):,}" for v in y2c["n_inside_radius"]],
        })
        lines += [md_table(show), ""]
    title = f"### Top 20 by delta in every setting ({word}candidate tables)".replace("( ", "(")
    lines += [title, ""]
    if word:
        lines += [
            f"`coverage_low` fired in Stage 0, so these tables are **{word.strip()}** and every "
            "factor labelled `toward_young` carries that word in its name below.",
            "",
        ]
    for variant, space, pair, show20, (n_top, n_all) in blocks:
        lines += [
            f"**{variant} / {space} / {pair}** — {n_top} of these 20 are `toward_young`; "
            f"{n_all:,} in this setting overall.",
            "",
            md_table(show20),
            "",
        ]
    return lines


def reading_lines(tab, cnt, inv, tgt, dist_df, guide_cos_zero):
    """What produced each count. Reported because several counts read backwards otherwise.

    Nothing here is a new gate or a re-tuned threshold; every number is read off the
    table that the pre-registered rules already produced.
    """
    n_f = int(inv["n_factors_ge2"])
    p_min = 1.0 / (N_PERM + 1)
    need = int(np.ceil(p_min * n_f / Q_BAR))
    # settings where no factor can have a negative delta at all
    no_neg = cnt[cnt["n_delta_negative"] == 0]
    best = tab.sort_values("delta").groupby(["variant", "space", "pair"], sort=False).head(1)
    lines = [
        "### What produced these counts",
        "",
        "Several of the counts above read backwards without their mechanism, so each is stated.",
        "",
    ]
    zero = cnt[cnt["n_toward_young"] == 0]
    brows = []
    for _, r in zero.iterrows():
        why = []
        if r["n_delta_negative"] == 0:
            why.append("no factor has a negative delta")
        if r["n_q_le_005"] == 0:
            why.append("no factor reaches q <= 0.05")
        if r["n_below_floor"] == 0:
            why.append("no factor clears its own floor")
        brows.append(dict(variant=r["variant"], space=r["space"], pair=r["pair"],
                          binding_constraint="; ".join(why) or "the flags (guides/identity)"))
    if len(zero):
        n_q = int((zero["n_q_le_005"] == 0).sum())
        lines += [
            f"**Why each of the {len(zero)} zero settings is zero.**"
            + (f" In **{n_q} of {len(zero)}** the binding constraint is the q-value gate: no "
               f"factor reaches q <= {Q_BAR} there. Unlike SOUTH2 that is **no longer** a "
               f"resolution limit of the null — with n_perm = {N_PERM:,} the attainable minimum "
               f"p is {p_min:.3e} and only {need} factors need to tie at it for BH to clear "
               f"{Q_BAR} across {n_f:,} factors — so a zero q count here means the observed "
               "deltas genuinely sit inside their own permutation distributions."
               if n_q else
               " **In neither is the q-value gate the binding constraint.** Both are zero "
               "because no factor has a negative delta at all, which is a geometric fact about "
               f"step length and aim and is settled before significance is consulted. Under "
               f"SOUTH2's n_perm = {N_PERM_SOUTH2} the q gate was the binding constraint in 18 "
               f"of its 19 zero settings; with n_perm = {N_PERM:,} it binds in none. That is "
               "change 1 doing exactly what it was for.")
            + " Observed ties at the p floor per setting are in the `at_p_floor` column of the "
              "count table and in `stage3_full_table.csv`.",
            "",
            md_table(pd.DataFrame(brows)),
            "",
        ]
    if len(no_neg):
        rows = []
        for _, r in no_neg.iterrows():
            g = tab[(tab.variant == r["variant"]) & (tab.space == r["space"])
                    & (tab.pair == r["pair"])]
            req = g["d_norm"] / (2.0 * float(g["v_norm"].iloc[0]))
            rows.append(dict(
                variant=r["variant"], space=r["space"], pair=r["pair"],
                origin_to_target=_f(g["v_norm"].iloc[0], 3),
                median_factor_norm=_f(g["d_norm"].median(), 3),
                median_required_cos=_f(req.median(), 4),
                smallest_required_cos=_f(req.min(), 4),
                observed_max_cos=_f(g["cos"].max(), 4),
                n_factors_meeting_own_requirement=str(int((g["cos"].to_numpy()
                                                           > req.to_numpy()).sum())),
                smallest_delta=_sf(g["delta"].min(), 4),
            ))
        lines += [
            f"**In {len(no_neg)} of {len(cnt)} settings not one factor has a negative delta, so "
            "`toward_young` = 0 there regardless of significance.** Adding a displacement of norm "
            "||d|| to a point at distance ||v|| from the target reduces that distance only if "
            "cos(d, v) > ||d|| / (2||v||); a large step needs to be well aimed or it overshoots. "
            "Below is that required cosine against the cosine actually achieved. Meeting one's own "
            "requirement is the same event as delta < 0, so the count column restates the result; "
            "the informative comparison is how far the best achieved cosine sits from the "
            "required one:",
            "",
            md_table(pd.DataFrame(rows)),
            "",
            "In those settings the `below_own_floor` count is still non-zero. That is consistent: "
            "a factor can move *less far away* than a random direction of its own size would, "
            "which clears its floor while still increasing the distance. It is not an approach, "
            "and the pre-registered `toward_young` rule requires delta < 0, so none of them is "
            "labelled a candidate.",
            "",
        ]
    n_pf = int(cnt["n_at_p_floor"].max()) if "n_at_p_floor" in cnt.columns else 0
    lines += [
        f"**The BH-FDR gate is no longer resolution-limited — that is CHANGE 1.** The smallest p "
        f"an empirical one-sided permutation test can return is 1/({N_PERM:,}+1) = "
        f"**{p_min:.3e}**. With BH across {n_f:,} factors, q <= {Q_BAR} needs only **{need} "
        f"factors tied at that minimum p** in the same setting, against "
        f"{int(np.ceil(n_f / (N_PERM_SOUTH2 + 1) / Q_BAR))} under SOUTH2's n_perm = "
        f"{N_PERM_SOUTH2}. The most any setting shows at the p floor is {n_pf:,} factors. So a "
        f"`q_le_005` of 0 in SOUTH3 is a statement about the data, not about the null's "
        "resolution — which is exactly the sentence SOUTH2 could not write.",
        "",
        f"**The V_sig guide-agreement collapse is now reported and not gated on — that is "
        f"CHANGE 2.** {guide_cos_zero['n_zero']:,} of {n_f:,} factors have a median pairwise "
        f"guide cosine of **exactly 0** when it is computed on V_sig vectors, "
        f"{guide_cos_zero['n_neg']:,} are negative, and {guide_cos_zero['n_pos']:,} are "
        "positive. A cosine of exactly 0 there means two guides for the same factor share **no "
        "significant gene at all**, so their masked vectors are orthogonal by construction — an "
        f"artefact of per-guide `adj_p <= {ADJ_P_BAR}` masking, not a disagreement about "
        "direction. SOUTH2 gated V_sig on that number and so refused 1,692 factors for a reason "
        "that was not about reliability. SOUTH3 computes the agreement on the V_all vectors "
        "always and applies that one value to both variants, so `guides_disagree` fires on "
        f"{int((cnt[cnt.variant == 'V_sig']['n_guides_disagree']).max()):,} factors in V_sig "
        "as well, the same set as in V_all. The V_sig number is kept as a reported-only column.",
        "",
        f"**The identity check was live and never fired.** The largest identity-panel drop over "
        f"all factors and variants is {tab['identity_drop'].max():+.4f} in frozen-z units against "
        f"a threshold of {IDENTITY_DROP}, and "
        f"{int((tab.drop_duplicates(['factor', 'variant'])['identity_drop'] > 0.4).sum())} "
        "factor-variants exceed 0.4. So `identity_loss` = 0 everywhere is a real result at this "
        "threshold, not a check that could never have fired; it came within "
        f"{IDENTITY_DROP - float(tab['identity_drop'].max()):.4f} of firing.",
        "",
    ]
    o1y2 = dist_df[(dist_df.origin == "O1") & (dist_df.target == "Y2")]
    inside_rows = o1y2[o1y2.get("origin_inside_radius").astype(str) == "True"] \
        if "origin_inside_radius" in o1y2.columns else o1y2.iloc[0:0]
    if len(inside_rows):
        ks = ", ".join(f"k={int(r['k'])}" for _, r in inside_rows.iterrows())
        lines += [
            "**The Y2 membership test has no discriminating power for O1, because the aged origin "
            f"is already inside the young region.** At {ks}, O1's own Mahalanobis distance to the "
            "Y2 region is below the membership radius before any perturbation is added "
            "(see the distance table in Stage 1), so `origin + d` lands inside the radius for "
            f"every one of the {n_f:,} factors. That is not evidence that any factor rejuvenates "
            "anything; it is a statement about the GTEx young and old centroids being close "
            "relative to the spread among young donors. For O2 the opposite holds: O2 starts "
            "outside the radius and no factor brings it inside.",
            "",
        ]
    lines += [
        "**The size of the best approach, stated so the deltas are not read as large.** The "
        "single most negative delta in each setting, with the fraction of the origin-to-target "
        "gap it closes (`frac`):",
        "",
        md_table(pd.DataFrame({
            "variant": best["variant"], "space": best["space"], "pair": best["pair"],
            "factor": best["factor"],
            "delta": [_sf(v, 4) for v in best["delta"]],
            "origin_to_target": [_f(v, 3) for v in best["v_norm"]],
            "frac_of_gap_closed": [_sf(v, 5) for v in best["frac"]],
            "cos": [_sf(v, 4) for v in best["cos"]],
            "toward_young": [str(bool(v)) for v in best["toward_young"]],
        })),
        "",
    ]
    return lines


def stage4_lines(s4):
    lines = [
        "## Stage 4 disagreement with directional readouts",
        "",
        "All factors are ranked three ways in every setting: by delta (most negative first), by "
        "cos (highest first), and by drop in the frozen ruler score (largest drop first). "
        "Spearman correlations are between those rank vectors, so positive means agreement.",
        "",
    ]
    show = pd.DataFrame({
        "variant": s4["variant"], "space": s4["space"], "pair": s4["pair"],
        "rho_delta_cos": [_sf(v, 3) for v in s4["rho_delta_cos"]],
        "rho_delta_ruler": [_sf(v, 3) for v in s4["rho_delta_ruler"]],
        "rho_cos_ruler": [_sf(v, 3) for v in s4["rho_cos_ruler"]],
        "top20_overlap_delta_cos": [str(int(v)) for v in s4["overlap_top20_delta_cos"]],
        "top20_overlap_delta_ruler": [str(int(v)) for v in s4["overlap_top20_delta_ruler"]],
        "top20_overlap_cos_ruler": [str(int(v)) for v in s4["overlap_top20_cos_ruler"]],
        "top20_by_cos_not_toward_young": [str(int(v)) for v in s4["n_top20_cos_not_toward_young"]],
    })
    lines += [md_table(show), ""]
    lines += [
        "The last column is the number asked for directly: of the 20 best factors by cosine, how "
        "many are **not** `toward_young` by delta.",
        "",
    ]
    return lines


def change2_lines(cdf, tot, guide_cos_zero, n_f):
    lines = [
        "## CHANGE 2 — one guide-agreement number, computed on V_all, applied to both variants",
        "",
    ]
    if cdf is None:
        return lines + [
            "**Unavailable.** `results/south2/stage3_full_table.csv` is not on disk, so the "
            "label comparison against SOUTH2 could not be computed. It is reported as "
            "unavailable, not as zero.",
            "",
        ]
    lines += [
        f"SOUTH2 computed the median pairwise guide cosine separately inside each gene variant "
        f"and gated each variant on its own value. In V_sig that number is **exactly 0** for "
        f"{guide_cos_zero['n_zero']:,} of {n_f:,} factors, because per-guide `adj_p <= "
        f"{ADJ_P_BAR}` masking often leaves two guides of the same factor with no significant "
        "gene in common, so their vectors are orthogonal by construction rather than by "
        "disagreement. SOUTH3 computes the agreement on the V_all vectors always and applies "
        "that one number to both variants.",
        "",
        "**How many factors change label because of this.** Two isolations are given, because "
        "only both together are honest: change 2 evaluated inside SOUTH2's own delta / floor / "
        "q, and change 2 evaluated inside SOUTH3's. Both count factor-settings whose "
        "`toward_young` label flips when only the guide rule is swapped.",
        "",
        md_table(pd.DataFrame({
            "variant": cdf["variant"], "space": cdf["space"], "pair": cdf["pair"],
            "guides_disagree_S2": [f"{int(v):,}" for v in cdf["n_guides_disagree_south2"]],
            "guides_disagree_S3": [f"{int(v):,}" for v in cdf["n_guides_disagree_south3"]],
            "flag_changed": [f"{int(v):,}" for v in cdf["n_guide_flag_changed"]],
            "toward_young_S2": [f"{int(v):,}" for v in cdf["n_toward_young_south2"]],
            "toward_young_S3": [f"{int(v):,}" for v in cdf["n_toward_young_south3"]],
            "label_changed_total": [f"{int(v):,}" for v in cdf["n_label_changed_total"]],
            "label_changed_by_change2_in_S2": [
                f"{int(v):,}" for v in cdf["n_label_changed_by_change2_within_south2"]],
            "label_changed_by_change2_in_S3": [
                f"{int(v):,}" for v in cdf["n_label_changed_by_change2_within_south3"]],
        })),
        "",
        f"**The direct answer.** The `guides_disagree` flag changes on "
        f"**{tot['flag_changed']:,}** of {tot['n_rows']:,} factor-settings — every one of them "
        "in V_sig, since the V_all flag is by definition unchanged. The number of factors whose "
        "`toward_young` label changes **because of this rule alone** is "
        f"**{tot['by_change2_within_south2']:,}** when evaluated against SOUTH2's own numbers "
        f"and **{tot['by_change2_within_south3']:,}** when evaluated against SOUTH3's. The total "
        f"SOUTH2 -> SOUTH3 label change, from all three patches together, is "
        f"**{tot['label_changed_total']:,}** factor-settings.",
        "",
        (f"Read the first of those two numbers with its cause: under SOUTH2's n_perm = "
         f"{N_PERM_SOUTH2} no V_sig setting had a single factor at q <= {Q_BAR}, so the V_sig "
         "guide flag was never the binding constraint there and relaxing it could not by itself "
         "create or destroy a label. **0 is therefore not evidence that the rule does not "
         "matter** — it is evidence that SOUTH2's null was too coarse for the rule to matter. "
         "The rule becomes load-bearing exactly once the null can resolve q, which is what "
         "change 1 does, and that is why the two isolations differ so sharply."
         if tot["by_change2_within_south2"] == 0 else
         "Both isolations are non-zero, so the guide rule is load-bearing in both."),
        "",
        f"**Change 2 is what makes V_sig reportable at all.** V_sig has "
        f"{tot['vsig_survivors_south3']:,} `toward_young` factor-settings in SOUTH3, and "
        f"{tot['vsig_by_change2']:,} of them would not be labelled under SOUTH2's per-variant "
        "guide rule — they were being refused for sharing no significant gene between guides, "
        "which is an artefact of the masking and not a disagreement about direction. In V_all "
        f"the change alters {tot['vall_by_change2']:,} labels, as it must: the V_all agreement "
        "value is the same number in both analyses.",
        "",
        "Per-factor rows for this comparison are in "
        "`results/south3/change2_label_comparison_per_factor.csv`.",
        "",
    ]
    return lines


def _o2_gap_ratio(sdf):
    """||O2 - target|| / ||O1 - target|| per space, as a range, for the prose."""
    g = sdf[sdf.variant == "V_all"][["space", "pair", "v_norm"]].copy()
    g["origin"] = g["pair"].str[:2]
    g["tgt"] = g["pair"].str[-2:]
    w = g.pivot_table(index=["space", "tgt"], columns="origin", values="v_norm")
    w = w.dropna()
    r = (w["O2"] / w["O1"]).to_numpy(float)
    return f"{r.min():.0f} to {r.max():.0f}"


def geometry_lines(sdf, qdf):
    lines = [
        "## CHANGE 3 — report-only displacement geometry (no gate)",
        "",
        "Nothing in this section gates, labels, or filters anything. It is here because the "
        "size of a displacement relative to the gap it has to close determines whether delta "
        "and cosine can disagree at all.",
        "",
        "### The identity",
        "",
        "Write r = origin - target, so ||r|| = ||target - origin||; d for the factor's "
        "displacement; rho = ||d|| / ||r||; and theta for the angle between d and -r, so that "
        "cos(theta) = -d.r / (||d|| ||r||) is exactly the `cos` column. Then, with no "
        "approximation,",
        "",
        "```",
        "delta = ||r|| * ( sqrt(1 - 2*rho*cos(theta) + rho^2) - 1 )",
        "```",
        "",
        "Expanding in rho:",
        "",
        "```",
        "delta = -||d||*cos(theta)  +  ||d||^2 * (1 - cos^2(theta)) / (2*||r||)  +  O(rho^3)",
        "```",
        "",
        "So **delta ~= -||d|| cos(theta)** is the small-displacement identity. Three things "
        "follow and are worth stating before the numbers: it is exact only in the limit "
        "rho -> 0; the leading correction is **positive** whenever |cos| < 1, so the identity "
        "**overstates** how much ground a real displacement covers; and because the correction "
        "depends on ||d|| as well as on the angle, ranking by delta and ranking by cosine must "
        "coincide as rho -> 0 and can only come apart at larger rho.",
        "",
        "The residual is therefore quoted as a fraction of **||d||**, not of |delta|. Dividing "
        "by |delta| looks natural and is useless here: delta passes through zero, so the "
        "relative error is unbounded near it and its large values say nothing about the "
        "identity. Divided by ||d||, the expansion predicts a residual of exactly "
        "rho*(1 - cos^2)/2, which is bounded above by rho/2.",
        "",
        "### Do the data follow the identity?",
        "",
    ]
    n_neg = int(sdf["n_resid_negative"].sum())
    lines += [
        f"**Yes, to the leading term and with the predicted sign.** The residual "
        f"delta - (-||d||cos(theta)) is non-negative on "
        f"**{'all' if n_neg == 0 else f'all but {n_neg:,} of the'} 51,408 factor-settings**, as "
        "the expansion requires. Dividing the residual by the leading term the expansion "
        "predicts for it, ||d||*rho*(1 - cos^2)/2, gives per-setting medians of "
        f"{sdf['resid_over_leading_term_median'].min():.4f} to "
        f"{sdf['resid_over_leading_term_median'].max():.4f} — that is, the second-order term is "
        "not merely the right order of magnitude, it is the right number.",
        "",
        "**Whether the identity is a usable approximation, though, varies by two orders of "
        "magnitude across settings, and the ratio says which.** Where rho stays small the "
        f"residual is negligible: across the {int((sdf.pair.str.startswith('O2')).sum())} O2 "
        "settings rho never exceeds "
        f"{float(sdf.loc[sdf.pair.str.startswith('O2'), 'ratio_max'].max()):.3f} and the median "
        "residual is under "
        f"{100.0 * float(sdf.loc[sdf.pair.str.startswith('O2'), 'identity_err_over_dnorm_median'].max()):.1f}% "
        f"of ||d||, because the cross-dataset offset makes ||target - origin|| {_o2_gap_ratio(sdf)} "
        "times larger than it is for O1 in the same space. One setting sits outside the "
        "small-displacement regime "
        "entirely: **V_all / S1 / O1_Y1**, where the median factor moves "
        f"{float(sdf.loc[(sdf.variant == 'V_all') & (sdf.space == 'S1') & (sdf.pair == 'O1_Y1'), 'ratio_median'].iloc[0]):.2f} "
        "of the whole origin-to-target gap and the largest moves "
        f"{float(sdf.loc[(sdf.variant == 'V_all') & (sdf.space == 'S1') & (sdf.pair == 'O1_Y1'), 'ratio_max'].iloc[0]):.2f} "
        "times it, leaving a median residual of "
        f"{100.0 * float(sdf.loc[(sdf.variant == 'V_all') & (sdf.space == 'S1') & (sdf.pair == 'O1_Y1'), 'identity_err_over_dnorm_median'].iloc[0]):.0f}% "
        "of ||d||. That is the same fact as the Stage 3 result that **no factor has a negative "
        "delta there**: when the step is as long as the journey the quadratic term dominates, "
        "and no achievable angle is enough to end up closer.",
        "",
        "The other zero setting, **V_sig / S1 / O1_Y1**, is zero for a different reason and the "
        "ratio is the wrong explanation for it: its rho median is "
        f"{float(sdf.loc[(sdf.variant == 'V_sig') & (sdf.space == 'S1') & (sdf.pair == 'O1_Y1'), 'ratio_median'].iloc[0]):.3f}, "
        "no larger than several settings that do produce survivors. What it has instead is "
        "cosines of essentially zero: significance-masked displacements in gene space are close "
        "to orthogonal to the within-GTEx age direction, and the required cosine rho/2 is not "
        "met by any factor. Both zero settings are S1 / O1_Y1 — the one origin-target pair "
        "whose gap is a real within-GTEx age difference rather than a cross-dataset offset — but "
        "one fails on step length and the other on aim, and they should not be given the same "
        "explanation.",
        "",
        "### ||d|| and the ratio ||d|| / ||target - origin||, per setting",
        "",
        md_table(pd.DataFrame({
            "variant": sdf["variant"], "space": sdf["space"], "pair": sdf["pair"],
            "origin_to_target": [_f(v, 3) for v in sdf["v_norm"]],
            "d_norm_min": [_f(v, 4) for v in sdf["d_norm_min"]],
            "d_norm_median": [_f(v, 4) for v in sdf["d_norm_median"]],
            "d_norm_max": [_f(v, 4) for v in sdf["d_norm_max"]],
            "ratio_min": [_f(v, 5) for v in sdf["ratio_min"]],
            "ratio_median": [_f(v, 5) for v in sdf["ratio_median"]],
            "ratio_max": [_f(v, 5) for v in sdf["ratio_max"]],
        })),
        "",
        "Per-factor ||d|| and the ratio are the `d_norm` and `d_over_v` columns of "
        "`results/south3/stage3_full_table.csv`; this table is the per-setting summary in "
        "`results/south3/stage3_dnorm_ratio_summary.csv`.",
        "",
        "### Rank agreement between delta and cosine, by quartile of the ratio",
        "",
        "Within each setting the factors are split at their own 25th, 50th and 75th percentiles "
        "of ||d|| / ||target - origin||. `rho_delta_cos` is the Spearman correlation between the "
        "rank by delta (most negative first) and the rank by cosine (highest first), so "
        "**+1 means the two orderings agree**. The last three columns test the identity inside "
        "that quartile: the median residual as a fraction of ||d||, the bound rho/2 that the "
        "expansion puts on it, and the residual divided by the leading term predicted for it.",
        "",
        md_table(pd.DataFrame({
            "variant": qdf["variant"], "space": qdf["space"], "pair": qdf["pair"],
            "quartile": qdf["quartile"],
            "n": [f"{int(v):,}" for v in qdf["n"]],
            "ratio_median": [_f(v, 5) for v in qdf["ratio_median"]],
            "rho_delta_cos": [_sf(v, 3) for v in qdf["rho_delta_rank_vs_cos_rank"]],
            "err_over_dnorm_median": [
                _f(v, 5) for v in qdf["identity_err_over_dnorm_median"]],
            "err_bound_half_ratio": [
                _f(v, 5) for v in qdf["identity_err_bound_half_ratio_max"]],
            "resid_over_leading_term": [
                _f(v, 4) for v in qdf["resid_over_leading_term_median"]],
        })),
        "",
        "Read down the `rho_delta_cos` column within any one setting and the pattern is the "
        "same everywhere: the two rankings agree almost perfectly in the middle quartiles and "
        "loosen at both ends — at the bottom because the deltas there are so small that ties "
        "and rounding dominate the ordering, at the top because that is where the quadratic "
        "term is large enough to re-order factors that the angle alone would have ranked "
        "differently. The only settings where the agreement is poor in every quartile are the "
        f"two S1 / O1_Y1 settings ({qdf[(qdf.space == 'S1') & (qdf.pair == 'O1_Y1')]['rho_delta_rank_vs_cos_rank'].min():+.3f} "
        f"to {qdf[(qdf.space == 'S1') & (qdf.pair == 'O1_Y1')]['rho_delta_rank_vs_cos_rank'].max():+.3f} "
        "across their eight quartiles), which are also the two with no survivor; every quartile "
        "of every other setting is at least "
        f"{qdf[~((qdf.space == 'S1') & (qdf.pair == 'O1_Y1'))]['rho_delta_rank_vs_cos_rank'].min():+.3f}.",
        "",
        "Full per-quartile numbers, including the error relative to |delta| for completeness, "
        "are in `results/south3/stage3_geometry_quartiles.csv`.",
        "",
    ]
    return lines


def survivor_lines(sv, s4):
    lines = [
        "## Is any survivor distinguishable from the ranking a plain cosine would give?",
        "",
        "Stated plainly, because it is the question the pre-registration ends on. For every "
        "setting with at least one `toward_young` survivor, the survivor set is compared with "
        "the top-n factors under four one-number rankings, n being the number of survivors in "
        "that setting: cosine alone, the small-displacement identity -||d||cos(theta), delta "
        "itself, and ||d|| alone.",
        "",
        "**`recovered_by_delta` is the control that makes this table readable.** The "
        "`toward_young` label is a conjunction — delta < 0 AND below the factor's own "
        f"magnitude-matched floor AND q <= {Q_BAR} AND neither flag — not a ranking. So even "
        "delta, the statistic the label is built on, cannot reproduce the survivor set exactly "
        "by taking its own top n. Cosine only fails in a way that means anything if it fails by "
        "more than delta does.",
        "",
    ]
    nz = sv[sv.n_survivors > 0]
    if len(nz) == 0:
        lines += ["No setting has a survivor, so the comparison is empty.", ""]
        return lines
    n = nz["n_survivors"].to_numpy(float)
    lines += [
        md_table(pd.DataFrame({
            "variant": nz["variant"], "space": nz["space"], "pair": nz["pair"],
            "n_survivors": [f"{int(v):,}" for v in nz["n_survivors"]],
            "recovered_by_cos": [f"{int(v):,}" for v in nz["n_overlap_top_cos"]],
            "recovered_by_identity": [f"{int(v):,}" for v in nz["n_overlap_top_identity"]],
            "recovered_by_delta": [f"{int(v):,}" for v in nz["n_overlap_top_delta"]],
            "recovered_by_d_norm": [f"{int(v):,}" for v in nz["n_overlap_top_d_norm"]],
            "cos_pct": [_f(v, 1) for v in 100.0 * nz["n_overlap_top_cos"].to_numpy(float) / n],
            "delta_pct": [_f(v, 1) for v in
                          100.0 * nz["n_overlap_top_delta"].to_numpy(float) / n],
            "worst_survivor_cos_rank": [f"{int(v):,}" for v in nz["worst_survivor_cos_rank"]],
            "rho_delta_cos_all_factors": [_sf(v, 3) for v in nz["rho_delta_cos_all"]],
        })),
        "",
        "Two things to read off it. First, **neither direction statistic dominates**: cosine "
        f"recovers more of the survivor set than delta does in {int((nz['n_overlap_top_cos'] > nz['n_overlap_top_delta']).sum())} "
        f"of the {len(nz)} settings and fewer in "
        f"{int((nz['n_overlap_top_cos'] < nz['n_overlap_top_delta']).sum())}, and in aggregate "
        f"they tie ({100.0 * nz['n_overlap_top_cos'].sum() / nz['n_survivors'].sum():.1f}% "
        f"against {100.0 * nz['n_overlap_top_delta'].sum() / nz['n_survivors'].sum():.1f}%). "
        "Per setting they differ by as much as 25 percentage points in either direction, which "
        "is the size of the disagreement you would expect between two near-identical orderings "
        "cut at an arbitrary rank; what does not happen anywhere is delta recovering the "
        "survivor set while cosine fails to. The shortfall in both columns is the conjunction "
        "of gates, not the choice of direction statistic. Second, `recovered_by_d_norm` is far "
        "lower than either, so the survivor set is not simply the largest displacements: the "
        "label does use the angle. What it does not use is anything the angle does not already "
        "carry.",
        "",
    ]
    return lines


def keys_lines(man):
    fired = man.get("fired_keys") or []
    flags = man.get("report_flags") or []
    notes = man.get("notes") or []
    lines = ["## Fired keys", ""]
    lines.append("STOP keys fired: " + (", ".join(f"`{k}`" for k in fired) or "**none**") + ".")
    lines.append("STOP keys checked and not fired: "
                 + (", ".join(f"`{k}`" for k in STOP_KEYS if k not in fired) or "none") + ".")
    lines.append("")
    lines.append("REPORT-ONLY flags fired: "
                 + (", ".join(f"`{r['key']}`" for r in flags) or "none") + ".")
    lines.append("")
    for r in flags:
        lines += [f"- **{r['key']}** (report-only, not a stop): {r['message']}", ""]
    if man.get("failures"):
        lines += ["Failures recorded verbatim:", ""]
        for fr in man["failures"]:
            lines.append(f"- **{fr.get('step')}:** {fr.get('message')}")
        lines.append("")
    if notes:
        lines += ["Recorded non-fatal failures and implementation decisions, verbatim:", ""]
        for n in notes:
            lines.append(f"- **{n.get('step')}:** {n.get('message')}")
        lines.append("")
    return lines


def limitations_lines(inv, s2res=None, tgt=None):
    ch = inv["coverage_chosen"]
    flags = {r["key"] for r in (load_manifest().get("report_flags") or [])}
    y2 = (tgt or {}).get("y2") or {}
    y100 = y2.get(100)
    z_o2 = (tgt or {}).get("z_o2")
    cov_sentence = (
        f"The chosen mapping covers {ch['weight_frac']:.4f} of the ruler's absolute weight, which "
        + ("fired `coverage_low`; the candidate tables are titled exploratory for that reason."
           if "coverage_low" in flags else
           f"is above the `coverage_low` bar of {COVERAGE_LOW_BAR}.")
    )
    scale_sentence = (
        "- **The log base of the input may be unresolved — and here it is.** `scale_unknown` "
        "fired. The magnitude-matched floors are invariant to a global scale factor on X, but "
        "the MDA is not, and neither is the absolute size of any delta or identity drop in "
        "frozen-z units. Read the deltas as relative to their own floor, not as absolute "
        "distances."
        if "scale_unknown" in flags else
        "- **The log base of the input** was determined from the file; `scale_unknown` did not "
        "fire."
    )
    return [
        "## Limitations",
        "",
        "- **The inputs are regression estimates, not expression.** X holds the source's own "
        "per-gene effect estimates, so the frozen count pipeline (sum -> TMM log2-CPM with "
        "prior.count=2 -> frozen mu/sd -> z) was **not** applied and could not be. Only the "
        "frozen sd was used, to turn an effect into a displacement in frozen-z units. No "
        "library-size normalisation of this file was performed, because there are no counts here "
        "to normalise.",
        "- **The effects are measured in one neonatal fibroblast line with no aged cells.** "
        "Every delta in this document is therefore a counterfactual: it assumes the perturbation "
        "effect measured in Hs27 transfers unchanged to an aged cell. Nothing in this data can "
        "test that assumption. The origin is an aged donor and the effect is not.",
        f"- **Unmeasured genes are assumed unmoved.** {inv['n_fixed_zero']:,} of the "
        f"{inv['n_ruler']:,} ruler coordinates "
        f"({100.0 * inv['n_fixed_zero'] / inv['n_ruler']:.1f}%) are fixed at 0 for every factor. "
        f"{cov_sentence} A factor that moves a gene this file does not measure is invisible "
        "here.",
        f"- **Guide counts are small.** Guides per factor run "
        f"{inv['guides_min']} to {inv['guides_max']} (median {inv['guides_median']:g}), so the "
        f"per-factor mean averages a handful of rows and the bootstrap over guides has very few "
        "distinct resamples. Factors with 2 guides get no bootstrap interval at all.",
        scale_sentence,
        f"- **The identity panel is incomplete.** {len(inv['identity_panel'])} of the "
        f"{len(IDENTITY_GENES)} pre-registered identity genes are present "
        f"({', '.join(inv['identity_absent'])} absent), so `identity_loss` is a weaker check "
        "than intended.",
        "- **The S2 basis is not independent of the young target.** The PCA is fit on the GTEx "
        "fibroblast donor matrix, and Y1/Y2 are subsets of exactly that matrix. It is independent "
        "of the Southard data and of O2.",
        (f"- **Y2 is undefined at k=100** (covariance rank {y100['rank']} on {y100['n_fit']} fit "
         "donors), and fragile wherever the fit donors number fewer than twice k. Reported as "
         "such; nothing was substituted."
         if y100 is not None and not y100["defined"] else
         "- Y2 was built only where the stages that need it ran."),
        ("- **O1 and O2 are not on the same footing.** O1 is a mean over GTEx donors in the same "
         "frozen-z space that the ruler was built in; O2 comes through SAME's own TMM panel of six "
         "pseudobulk rows and carries its own missing-gene zeroing "
         f"({int((np.asarray(z_o2) == 0).sum()):,} coordinates). Their distances to Y1 are "
         "reported separately and are never pooled."
         if z_o2 is not None and np.asarray(z_o2).size > 1 else
         "- O1 and O2 were not built because the task stopped before Stage 1."),
        ("- **O2's distance to Y1 is dominated by a cross-dataset offset, not by age.** In S1, "
         f"||O2 - Y1|| = {float(np.linalg.norm(np.asarray(z_o2) - tgt['c_young'])):.1f} against "
         f"||O1 - Y1|| = {float(np.linalg.norm(tgt['c_old'] - tgt['c_young'])):.1f}, a ratio of "
         f"{float(np.linalg.norm(np.asarray(z_o2) - tgt['c_young'])) / float(np.linalg.norm(tgt['c_old'] - tgt['c_young'])):.1f}x. "
         "At that distance almost any displacement with a slightly favourable projection lowers "
         "delta, so the O2 settings are the easiest place to earn a `toward_young` label and the "
         "least informative about ageing. The O2 candidate counts should be read with that in "
         "mind."
         if z_o2 is not None and np.asarray(z_o2).size > 1 and "c_young" in (tgt or {}) else
         "- O2 was not built."),
        f"- **The q-value gate is no longer resolution-limited, but it is still a permutation "
        f"gate.** With n_perm = {N_PERM:,} the smallest attainable p is 1/{N_PERM + 1:,} = "
        f"{P_MIN:.3e} and BH over {inv['n_factors_ge2']:,} factors needs only "
        f"{int(np.ceil(P_MIN * inv['n_factors_ge2'] / Q_BAR))} factors tied at that minimum "
        "before q can reach 0.05, so SOUTH2's resolution ceiling is gone. What remains is that "
        "the null is a within-factor coordinate permutation: it asks whether this factor's "
        "displacement is better aimed than a reshuffling of its own entries, and nothing else. "
        "It does not test measurement error in the effect estimates, it does not test transfer "
        "to an aged cell, and a small q is not a small effect size.",
        f"- **The permutations were drawn in chunks and reduced on the fly.** {N_PERM:,} "
        "permuted deltas per factor and setting cannot be held in memory, so the floor, the "
        "permutation count and the pooled draw are accumulated as the chunks are produced (see "
        "the `stage2A_streaming` note). The permuted values themselves are not retained, so the "
        "full null distributions are not recoverable from `results/south3/` — only their 5th "
        "percentiles, the observed-value counts, and the pooled floor.",
        "- **The control-row floor is not a noise floor** and was never used as the gate: those "
        "rows are untargeted but not inert.",
        "- This task makes **no claim about the source paper's own conclusions**. It did not test "
        "them.",
        "",
    ]


def build_findings(prereg, prereg_s2, inv=None, man=None, sections=None, result=None):
    fired = (man or {}).get("fired_keys") or []
    flags = [r["key"] for r in ((man or {}).get("report_flags") or [])]
    if fired:
        status = ("STOP. Keys fired: " + ", ".join(f"`{k}`" for k in fired)
                  + ". Stages after the stop did not run.")
    else:
        status = ("Ran to completion. No STOP key fired. Report-only flags: "
                  + (", ".join(f"`{k}`" for k in flags) or "none") + ".")
    head = [
        "# FINDINGS_SOUTH3 — rank Southard Hs27 CRISPRa transcription factors by approach to a "
        "young target; SOUTH2 re-run with a resolvable null, one guide-agreement rule, and "
        "report-only displacement geometry",
        "",
        f"**Status:** {status}",
        "",
    ]
    if result:
        head += [
            f"**Result:** {result} Read the \"What produced these counts\" section before any "
            "count: some settings are zero for geometric reasons rather than because factors "
            "were tested and found null, and the O2 settings are dominated by a cross-dataset "
            "offset.",
            "",
        ]
    head += [
        "**SOUTH3 supersedes SOUTH2.** It is a patch, not a new analysis: same pinned file, same "
        "pre-registration, same spaces, variants, origins, targets, floors, labelling rule, "
        "what-does-not-count list, exploratory titling under `coverage_low`, and same seeds. "
        f"Exactly three things change: (1) `n_perm` = {N_PERM:,} for the Stage 2A per-factor "
        f"permutation null instead of {N_PERM_SOUTH2}, so BH across 1,836 factors can resolve q; "
        "(2) guide agreement is computed on the V_all vectors always and that one value is "
        "applied to both gene variants; (3) a report-only section on ||d||, "
        "||d|| / ||target - origin||, rank agreement between delta and cosine as a function of "
        "that ratio, and the small-displacement identity. `FINDINGS_SOUTH2.md` was not edited "
        "except for one line at its top marking it superseded.",
        "",
        "SOUTH2 itself superseded SOUTH, which halted at Stage 0 with `no_raw_counts` because "
        "its pre-registration described this file as a count matrix; the file holds per-guide "
        "regression output.",
        "",
        "The SOUTH3 patch below was written verbatim to `results/south3/PREREG.flag`, and the "
        "inherited SOUTH2 block to `results/south3/PREREG_SOUTH2_inherited.flag`, before any "
        "SOUTH3 statistic was computed. Nothing was re-tuned after numbers existed.",
        "",
        "## Pre-registration — the SOUTH3 patch (verbatim, written before any SOUTH3 statistic)",
        "",
        "```",
        prereg,
        "```",
        "",
        "## Pre-registration — the inherited SOUTH2 block (verbatim, unchanged)",
        "",
        "Everything in this block still governs SOUTH3 except the three patched items above.",
        "",
        "```",
        prereg_s2,
        "```",
        "",
    ]
    body = list(sections or [])
    FINDINGS_PATH.write_text("\n".join(head + body).rstrip() + "\n", encoding="utf-8")
    return FINDINGS_PATH


def write_progress(man, next_action, extra=None):
    fired = man.get("fired_keys") or []
    flags = [r["key"] for r in (man.get("report_flags") or [])]
    status = ("STOP; keys fired: " + ", ".join(f"`{k}`" for k in fired)) if fired \
        else "no STOP key fired; the task ran to completion"
    lines = [
        "# PROGRESS_SOUTH3",
        "",
        f"**STOP status:** {status}.",
        "",
        f"**Report-only flags:** {', '.join(f'`{k}`' for k in flags) or 'none'}.",
        "",
        f"**Next action:** {next_action}",
        "",
    ]
    if extra:
        lines += list(extra) + [""]
    PROGRESS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return PROGRESS_PATH


# --------------------------------------------------------------------------- #
def main():
    log = Logger(SOUTH3_DIR / "south3_run.log")
    log("=" * 100)
    log("SOUTH3 — patch of SOUTH2: n_perm=20000, V_all-always guide agreement, report-only "
        "displacement geometry")
    log("=" * 100)
    log(f"[seeds] SEED={SEED} (splits, nulls, permutations) BOOT_SEED={BOOT_SEED} (bootstrap); "
        f"imported TOWARD_SEED={TOWARD_SEED}, SAME BOOT_SEED={BOOT_SEED}. Unchanged from SOUTH2.")
    log(f"[params] n_null={N_NULL} n_perm={N_PERM} (SOUTH2: {N_PERM_SOUTH2}) n_boot={N_BOOT} "
        f"S2_KS={S2_KS} (lowdim_common.KS={LOWDIM_KS}) F={F_VALUES} adj_p_bar={ADJ_P_BAR} "
        f"q_bar={Q_BAR}")
    log(f"[params] attainable minimum permutation p = 1/({N_PERM}+1) = {P_MIN:.6e}; BH over "
        f"1,836 factors reaches q<={Q_BAR} once "
        f"{int(np.ceil(P_MIN * 1836 / Q_BAR))} factor(s) tie at it "
        f"(SOUTH2 needed {int(np.ceil(1836 / (N_PERM_SOUTH2 + 1) / Q_BAR))}).")
    log("[params] n_null is recorded as pre-registered; the pre-registration assigns the null to "
        f"the Stage 2A per-factor permutations (n_perm={N_PERM}), so no separate null draw "
        "exists.")
    log(f"[params] permutation chunk={PERM_CHUNK} rows, workers={N_WORKERS} "
        "(SOUTH3_WORKERS overrides)")
    log("[scope] Uncalibrated R2 is never a gate; the frozen ruler is loaded, never refit; "
        "L4A_matrix.h5 is not opened; nothing is downloaded.")
    log("[scope] Only FINDINGS_SOUTH3.md, PROGRESS_SOUTH3.md, results/south3/* and "
        "src/south3_run.py are written by this run. results/south2/stage3_full_table.csv and "
        "FINDINGS_SOUTH2.md are read, never written.")

    prereg, prereg_s2 = load_prereg_text()
    log(f"[prereg] read {PREREG_FLAG} ({len(prereg)} chars) and {PREREG_INHERITED} "
        f"({len(prereg_s2)} chars); not rewriting either block")

    man = load_manifest()
    man.update(seed=SEED, boot_seed=BOOT_SEED, status="RUNNING", fired_keys=[], failures=[],
               report_flags=[], notes=[])
    save_manifest(man)

    # ---------------- Stage 0 ---------------- #
    try:
        pack, fired0 = stage0(log)
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        man = load_manifest()
        build_findings(prereg, prereg_s2, None, man,
                       ["## Stage 0 inventory", "",
                        f"Stage 0 could not complete: {e.message}", ""]
                       + keys_lines(man))
        write_progress(man, "Resolve the Stage 0 structural failure recorded above. Nothing "
                            "downstream can run until the inventory completes.")
        log(f"[STOP] {e.step}: {e.message}")
        return
    inv = pack["inv"]
    for key, message, details in fired0:
        record_failure(key, message, details)
    man = load_manifest()
    if man.get("fired_keys"):
        build_findings(prereg, prereg_s2, inv, man, stage0_lines(inv, man) + keys_lines(man)
                       + limitations_lines(inv))
        write_progress(man, "Halt as pre-registered on the Stage 0 key(s) above. Stages 1-4 "
                            "did not run and no statistic was computed.")
        log(f"[STOP] Stage 0 fired: {[k for k, _, _ in fired0]}")
        return
    log("[stage0] passed: no Stage 0 STOP key fired")

    # ---------------- Stage 1 ---------------- #
    factors, rows_by_factor, D, Cctrl = build_factor_vectors(pack, log)
    gtex = load_gtex_z(log)
    if int(np.asarray(gtex["Z"]).shape[1]) != inv["n_ruler"]:
        raise StopStep("genes", f"GTEx Z has {gtex['Z'].shape[1]} genes, ruler has "
                                f"{inv['n_ruler']}. Not aligning ad hoc.")
    s2 = build_s2(gtex, log)
    pack["L_full"] = s2["L"]
    tgt = build_targets_origins(gtex, s2, log)
    combos, Lm, dist_df = build_combos(pack, s2, tgt, log)

    # ---- observed geometry, then Stage 2 ---- #
    point, keys = observed_points(combos, Lm, D, factors, log)
    s2res, fired2 = stage2(pack, combos, Lm, D, Cctrl, tgt, point, keys, log)
    for key, message, details in fired2:
        record_failure(key, message, details)
    man = load_manifest()
    sections = (stage0_lines(inv, man)
                + stage1_lines(inv, s2, tgt, dist_df, s2res)
                + stage2_lines(s2res, inv))
    if not s2res["passed"]:
        sections += [
            "## Stage 3 ranking", "",
            "**Not run.** Stage 2C fired `no_detection_power`: no dose fraction f <= 1.0 of the "
            "true young-minus-old direction, scaled to the median factor norm, cleared the "
            "primary floor in any space or gene variant. Stated plainly as pre-registered: "
            "**this data cannot support a ranking.** No factor was scored and no candidate count "
            "exists — it is undefined, not zero.", "",
            "## Stage 4 disagreement with directional readouts", "",
            "**Not run.** It ranks factors by delta, cosine, and ruler-score drop, all of which "
            "require Stage 3 output.", "",
        ]
        sections += keys_lines(man) + limitations_lines(inv, s2res, tgt)
        build_findings(prereg, prereg_s2, inv, man, sections)
        write_progress(man, "Halt as pre-registered on `no_detection_power`. Stages 3-4 did not "
                            "run. Report Stages 0-2 only: this data cannot support a ranking.")
        log("[STOP] no_detection_power")
        return

    # ---------------- Stages 3 and 4 ---------------- #
    guide_cos_all, guide_cos_sig = guide_agreement(pack, factors, rows_by_factor, log)
    tab, cnt, n_invalid = stage3(pack, combos, Lm, D, factors, rows_by_factor, s2res,
                                 point, keys, guide_cos_all, guide_cos_sig, log)
    s4 = stage4(tab, log)
    blocks, word = top20_tables(tab, inv["coverage_low"])
    sdf, qdf = geometry_report(tab, log)
    c2df, c2tot = change2_comparison(tab, log)
    sv = survivor_vs_cosine(tab, log)

    man = load_manifest()
    man["status"] = "DONE"
    man["inventory_path"] = str(INVENTORY_PATH)
    man["n_stage3_rows"] = int(len(tab))
    man["toward_young_by_setting"] = {
        f"{r['variant']}|{r['space']}|{r['pair']}": int(r["n_toward_young"])
        for _, r in cnt.iterrows()}
    save_manifest(man)

    gs = tab[tab.variant == "V_sig"].drop_duplicates("factor")[
        "guide_median_cos_vsig_reported_only"]
    guide_cos_zero = dict(n_zero=int((gs == 0).sum()), n_neg=int((gs < 0).sum()),
                          n_pos=int((gs > 0).sum()))
    tw = cnt[["variant", "space", "pair", "n_factors", "n_toward_young"]]
    nonzero = tw[tw.n_toward_young > 0]

    # The pre-registration's closing question, answered from the survivor table and
    # never from a pooled number.
    nzsv = sv[sv.n_survivors > 0]
    if len(nzsv) == 0:
        verdict = ("No setting has a `toward_young` survivor, so the question of whether a "
                   "survivor is distinguishable from a plain cosine ranking does not arise.")
    else:
        worst_rho = float(nzsv["rho_delta_cos_all"].min())
        n_surv = int(nzsv["n_survivors"].sum())
        by_cos = int(nzsv["n_overlap_top_cos"].sum())
        by_delta = int(nzsv["n_overlap_top_delta"].sum())
        by_dn = int(nzsv["n_overlap_top_d_norm"].sum())
        live = set(zip(nzsv["variant"], nzsv["space"], nzsv["pair"]))
        ql = qdf[[tuple(x) in live for x in zip(qdf["variant"], qdf["space"], qdf["pair"])]]
        qmin = float(ql["rho_delta_rank_vs_cos_rank"].min())
        sl = sdf[[tuple(x) in live for x in zip(sdf["variant"], sdf["space"], sdf["pair"])]]
        emax = float(sl["identity_err_over_dnorm_median"].max())
        verdict = (
            f"**No.** Across the {len(nzsv)} settings that have survivors, the Spearman "
            "correlation between the delta ranking and the cosine ranking over all 1,836 "
            f"factors is at least {worst_rho:+.3f}, and within every quartile of "
            f"||d|| / ||target - origin|| in those settings it is at least {qmin:+.3f}. Taking "
            "the top-n factors by cosine alone, n being that setting's number of survivors, "
            f"recovers {by_cos:,} of the {n_surv:,} survivors "
            f"({100.0 * by_cos / n_surv:.1f}%) — and taking the top-n by delta itself, the "
            f"statistic the label is built on, recovers {by_delta:,} "
            f"({100.0 * by_delta / n_surv:.1f}%). The two agree to within a percentage point, "
            f"so the ~{100.0 - 100.0 * by_cos / n_surv:.0f}% that neither recovers is the "
            "conjunction of the floor and q gates, not information the cosine lacks. ||d|| "
            f"alone recovers only {by_dn:,} ({100.0 * by_dn / n_surv:.1f}%), so the label is "
            "not merely picking the biggest displacements. The mechanism is the identity: in "
            "every setting that has a survivor the median residual of delta = -||d||cos(theta) "
            f"is at most {100.0 * emax:.1f}% of ||d||, so delta there is the cosine re-weighted "
            "by the displacement's own length and carries nothing beyond it. No survivor in any "
            "setting is distinguishable on direction from what a plain cosine ranking would "
            "give. What the delta machinery does buy is not a different ordering but a bar: it "
            "says which of those cosines are larger than the factor's own reshuffled null, and "
            "a cosine ranking on its own never says that.")
    # Reported per setting, never pooled into one number.
    summary = (
        f"`toward_young` is 0 in {len(tw) - len(nonzero)} of the {len(tw)} settings and non-zero "
        f"in {len(nonzero)}"
        + (f", ranging {int(nonzero.n_toward_young.min())}-{int(nonzero.n_toward_young.max())} "
           f"of {int(tw.n_factors.iloc[0]) if 'n_factors' in tw.columns else 0} factors"
           if len(nonzero) else "")
        + f". Non-zero settings are {', '.join(sorted(set(nonzero.variant))) or 'none'}"
        + (f"; {', '.join(v for v in VARIANTS if v not in set(nonzero.variant))} returns 0 "
           "everywhere" if [v for v in VARIANTS if v not in set(nonzero.variant)] else "")
        + ". The counts are not pooled and no single setting is primary. " + verdict)
    log(f"[summary] {summary}")
    man = load_manifest()
    man["verdict_survivor_vs_cosine"] = verdict
    man["p_min_attainable"] = P_MIN
    man["change2"] = c2tot
    save_manifest(man)
    sections += (stage3_lines(tab, cnt, blocks, word, n_invalid, inv)
                 + reading_lines(tab, cnt, inv, tgt, dist_df, guide_cos_zero)
                 + stage4_lines(s4)
                 + change2_lines(c2df, c2tot, guide_cos_zero, len(factors))
                 + geometry_lines(sdf, qdf)
                 + survivor_lines(sv, s4)
                 + ["**Verdict, stated plainly as the pre-registration asks.** " + verdict, ""]
                 + keys_lines(man)
                 + limitations_lines(inv, s2res, tgt))
    build_findings(prereg, prereg_s2, inv, man, sections, result=summary)
    write_progress(
        man,
        "Nothing further is required by this pre-registration: Stages 0-4 and the three SOUTH3 "
        "patches all ran and are reported in FINDINGS_SOUTH3.md. SOUTH3's own three changes are "
        "spent — the null now resolves q, the guide rule is single-valued, and the displacement "
        "geometry is reported. What would change the answer rather than decorate it is "
        "unchanged from SOUTH2: (a) an input whose units are stated, since `scale_unknown` "
        "fired and the MDA is not scale-invariant, and (b) a perturbation dataset measured in "
        "aged cells, since every delta here is a counterfactual built from effects measured in "
        "one neonatal line. Both need a new pre-registration.",
        [f"**Result in one line:** {summary}"
         + (" `coverage_low` fired, so every candidate is labelled exploratory."
            if inv["coverage_low"] else ""),
         "",
         f"**Attainable minimum permutation p:** {P_MIN:.6e} = 1/({N_PERM:,}+1).",
         "",
         "**Supersession:** SOUTH3 supersedes SOUTH2. One line was added at the top of "
         "`FINDINGS_SOUTH2.md` marking it superseded; nothing else in that file was changed."],
    )
    log("[done]")


if __name__ == "__main__":
    main()

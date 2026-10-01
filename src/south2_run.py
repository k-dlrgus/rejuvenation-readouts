"""SOUTH2 — rank the Southard Hs27 CRISPRa transcription factors by approach to a
young fibroblast target, using the pinned file as what it actually is: per-guide
regression output (per-gene effect estimates with p / adj_p / masked layers).

SOUTH halted at Stage 0 with `no_raw_counts` because its pre-registration described
the same file as a count matrix. SOUTH2 supersedes it with a matching input
description. FINDINGS_SOUTH.md is not edited here beyond the one superseded line
that was added by hand before this script ran.

Reads results/survey/opened/southard/fibroblast_CRISPRa_mean_pop.h5ad (read-only).
Imports from src/toward_run.py, src/same_run.py, src/seng_run.py,
src/lowdim_common.py, src/md3_idtype.py without changing them.
Does not refit the frozen ruler. Does not download anything. Does not open
L4A_matrix.h5.

The pre-registration must already exist verbatim at results/south2/PREREG.flag,
written before any statistic. This script does not author that block; it reads it.

Writes only: FINDINGS_SOUTH2.md, PROGRESS_SOUTH2.md, results/south2/*.

Usage: python src/south2_run.py
"""
from __future__ import annotations

import re
import sys
import time
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


SOUTH2_DIR = RESULTS / "south2"
SOUTH2_DIR.mkdir(parents=True, exist_ok=True)
FINDINGS_PATH = ROOT / "FINDINGS_SOUTH2.md"
PROGRESS_PATH = ROOT / "PROGRESS_SOUTH2.md"
PREREG_FLAG = SOUTH2_DIR / "PREREG.flag"
MANIFEST_PATH = SOUTH2_DIR / "manifest.json"
INVENTORY_PATH = SOUTH2_DIR / "inventory.json"

SOUTHARD_DIR = RESULTS / "survey" / "opened" / "southard"
H5_PATH = SOUTHARD_DIR / "fibroblast_CRISPRa_mean_pop.h5ad"
PINNED_SIZE = 1_720_196_368
PINNED_MD5 = "ba44c7813903bb5df900348d6b0d589a"

SEED = 20260914
N_NULL = 200
N_PERM = 200
N_BOOT = 200
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
    if not PREREG_FLAG.exists():
        raise StopStep(
            "prereg",
            "PREREG.flag missing — the pre-registration must be written verbatim before "
            "any statistic. Not computing.",
        )
    return PREREG_FLAG.read_text(encoding="utf-8").rstrip("\n")


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
# md3_common.log_id_type. SOUTH2 may only write results/south2/*, so the sink is
# rebound in this process to SOUTH2's own manifest. The imported function itself
# is untouched on disk.
_ID_TYPES = {}


def _south2_id_sink(man_key, rec):
    _ID_TYPES[str(man_key)] = jsonable(rec)


md3_idtype.log_id_type = _south2_id_sink


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
    tree = dump_tree(f, SOUTH2_DIR / "h5_tree.txt")
    log(f"[h5] tree written ({len(tree)} lines) -> {SOUTH2_DIR / 'h5_tree.txt'}")
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
    id_rec = detect_id_type(var["gene_id"].astype(str).to_numpy(), log, "south2_var_gene_id",
                            str(H5_PATH) + "::var/gene_id")
    sym_rec = detect_id_type(var["gene_name"].astype(str).to_numpy(), log, "south2_var_gene_name",
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
    ctrl_tab.to_csv(SOUTH2_DIR / "stage0_controls.csv", index=False)

    tf_mask = ~ctrl_mask
    tf_labels = sorted(set(targets.to_numpy()[tf_mask]))
    per_factor = pd.Series(targets.to_numpy()[tf_mask]).value_counts()
    gmin, gmed, gmax = int(per_factor.min()), float(per_factor.median()), int(per_factor.max())
    n_ge2 = int((per_factor >= MIN_GUIDES).sum())
    log(f"[factors] distinct transcription factors={len(tf_labels):,}; guides per factor "
        f"min={gmin} median={gmed:g} max={gmax}; with >= {MIN_GUIDES} guides={n_ge2:,}")
    gdist = per_factor.value_counts().sort_index()
    pd.DataFrame({"n_guides": gdist.index.astype(int), "n_factors": gdist.to_numpy(int)}).to_csv(
        SOUTH2_DIR / "stage0_guides_per_factor.csv", index=False)
    pd.DataFrame({"factor": per_factor.index, "n_guides": per_factor.to_numpy(int)}).to_csv(
        SOUTH2_DIR / "stage0_factors.csv", index=False)

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
    pd.DataFrame(cov).to_csv(SOUTH2_DIR / "stage0_coverage.csv", index=False)
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
    df.to_csv(SOUTH2_DIR / "origin_target_distances.csv", index=False)
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
# STAGE 2 — detection limit
# --------------------------------------------------------------------------- #
def stage2(pack, combos, Lm, D, C, tgt, log):
    """A: per-factor magnitude-matched permutation floor (+ pooled floor).
    B: control-row floor (reported, never a gate).
    C: dose-response along the true young-minus-old direction, and the MDA.
    """
    n_f = D[VARIANTS[0]].shape[0]
    floors = {}          # (variant, pair_space) -> (n_f,) 5th percentile
    perm_store = {}      # (variant, pair_space) -> (n_f, N_PERM) deltas (kept in RAM)
    md_floors = {}       # Y2 only: the same permutations, Mahalanobis change
    md_perm_store = {}
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
    pooled_rows = []
    for variant in VARIANTS:
        M = D[variant]
        nrm = np.linalg.norm(M, axis=1)
        med = float(np.median(nrm))
        n_zero = int(np.sum(nrm < 1e-12))
        keys = [f"{c['space']}|{c['pair']}" for c in combos]
        acc = {kk: np.empty((n_f, N_PERM), dtype=np.float64) for kk in keys}
        acc_md = {kk: np.full((n_f, N_PERM), np.nan, dtype=np.float64)
                  for c, kk in zip(combos, keys) if c["target"] == "Y2" and c["y2"]["defined"]}
        pool = {kk: np.full((n_f, N_PERM), np.nan, dtype=np.float32) for kk in keys}
        t0 = time.perf_counter()
        for i in range(n_f):
            rng = np.random.default_rng([SEED, i])
            P = np.broadcast_to(M[i], (N_PERM, M.shape[1])).copy()
            P = rng.permuted(P, axis=1)
            Pp = P @ Lm
            s = (med / nrm[i]) if nrm[i] >= 1e-12 else np.nan
            for c, kk in zip(combos, keys):
                st = block_stats(c, P, Pp)
                acc[kk][i] = st["delta"]
                if kk in acc_md:
                    acc_md[kk][i] = st["md_change"]
                if np.isfinite(s):
                    pool[kk][i] = block_stats(c, P * s, Pp * s)["delta"]
            if (i + 1) % 400 == 0:
                log(f"[stage2A {variant}] {i+1}/{n_f} factors "
                    f"({time.perf_counter() - t0:.0f}s)")
        for c, kk in zip(combos, keys):
            perm_store[(variant, kk)] = acc[kk]
            floors[(variant, kk)] = np.percentile(acc[kk], 5, axis=1)
            if kk in acc_md:
                md_perm_store[(variant, kk)] = acc_md[kk]
                md_floors[(variant, kk)] = np.percentile(acc_md[kk], 5, axis=1)
            pl = pool[kk][np.isfinite(pool[kk])].astype(np.float64)
            pooled_rows.append(dict(
                variant=variant, space=c["space"], k=c["k"], pair=c["pair"],
                median_factor_norm=med, n_factors_pooled=int(n_f - n_zero),
                n_zero_norm_factors=n_zero, n_draws=int(pl.size),
                floor_p05=float(np.percentile(pl, 5)),
                null_median=float(np.median(pl)), null_min=float(pl.min()),
                null_max=float(pl.max()), n_perm=N_PERM,
            ))
        log(f"[stage2A {variant}] done in {time.perf_counter() - t0:.0f}s; median ||d||="
            f"{med:.4f}; zero-norm factors={n_zero}")
    pooled = pd.DataFrame(pooled_rows)
    pooled.to_csv(SOUTH2_DIR / "stage2a_pooled_floor.csv", index=False)

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
    fdf.to_csv(SOUTH2_DIR / "stage2a_per_factor_floor_summary.csv", index=False)

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
    cdf.to_csv(SOUTH2_DIR / "stage2b_control_floor.csv", index=False)
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
    ddf.to_csv(SOUTH2_DIR / "stage2c_dose_response.csv", index=False)

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
    mdf.to_csv(SOUTH2_DIR / "stage2c_mda.csv", index=False)
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
    return dict(floors=floors, perm_store=perm_store, md_floors=md_floors,
                md_perm_store=md_perm_store, pooled=pooled, per_factor_floor=fdf,
                control=cdf, dose=ddf, mda=mdf, passed=passed, v_true=v_true_full,
                u_full=u_full, v_norm_full=float(n_full), v_norm_res=float(n_res)), fired


# --------------------------------------------------------------------------- #
# STAGE 3 — ranking
# --------------------------------------------------------------------------- #
def stage3(pack, combos, Lm, D, factors, rows_by_factor, s2res, log):
    Gc, Sc = pack["Gc"], pack["Sc"]
    frozen = pack["frozen"]
    ruler_rows = pack["ruler_rows"]
    w_unit, _ = unit(np.asarray(frozen["w"], float))
    w_m = w_unit[ruler_rows]
    id_pos = pack["id_rowpos"]
    n_f = len(factors)
    keys = [f"{c['space']}|{c['pair']}" for c in combos]

    # per-factor, per-variant scalars that do not depend on the space
    guide_cos = {v: np.full(n_f, np.nan) for v in VARIANTS}
    id_drop = {v: np.full(n_f, np.nan) for v in VARIANTS}
    ruler_chg = {v: np.full(n_f, np.nan) for v in VARIANTS}
    n_guides = np.array([rows_by_factor[t].size for t in factors], dtype=int)

    point = {}   # (variant, key) -> dict of arrays
    boot = {}    # (variant, key) -> dict of arrays
    for variant in VARIANTS:
        M = D[variant]
        for kk in keys:
            point[(variant, kk)] = {c: np.full(n_f, np.nan) for c in
                                    ("delta", "cos", "frac", "d_norm", "md_change", "inside")}
            boot[(variant, kk)] = {c: np.full(n_f, np.nan) for c in
                                   ("delta_lo", "delta_hi", "cos_lo", "cos_hi",
                                    "frac_lo", "frac_hi")}
        t0 = time.perf_counter()
        for i, t in enumerate(factors):
            r = rows_by_factor[t]
            base = Gc[r] if variant == "V_all" else Gc[r] * Sc[r]
            d = M[i]
            dp = (d[None, :] @ Lm)
            st1 = {kk: block_stats(c, d[None, :], dp) for c, kk in zip(combos, keys)}
            for kk in keys:
                for cname in ("delta", "cos", "frac", "d_norm"):
                    point[(variant, kk)][cname][i] = float(st1[kk][cname][0])
                if "md_change" in st1[kk]:
                    point[(variant, kk)]["md_change"][i] = float(st1[kk]["md_change"][0])
                    point[(variant, kk)]["inside"][i] = float(bool(st1[kk]["inside"][0]))
            # guide agreement
            if base.shape[0] >= 2:
                cs = []
                for a in range(base.shape[0]):
                    for b in range(a + 1, base.shape[0]):
                        cs.append(cosine_full(base[a], base[b]))
                cs = np.asarray(cs, float)
                guide_cos[variant][i] = float(np.nanmedian(cs)) if np.isfinite(cs).any() else np.nan
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
        gd = int(np.sum(guide_cos[variant] <= 0))
        il = int(np.sum(id_drop[variant] > IDENTITY_DROP))
        log(f"[stage3 {variant}] guides_disagree (median pairwise guide cosine <= 0): {gd:,} of "
            f"{n_f:,}; identity_loss (panel drop > {IDENTITY_DROP} in frozen-z units): {il:,}")

    # assemble the full table
    out = []
    for variant in VARIANTS:
        for c, kk in zip(combos, keys):
            fl = s2res["floors"][(variant, kk)]
            perms = s2res["perm_store"][(variant, kk)]
            pt = point[(variant, kk)]
            bt = boot[(variant, kk)]
            p = np.array([permutation_p(pt["delta"][i], perms[i], greater=False)
                          for i in range(n_f)], dtype=float)
            q = bh_q(p)
            ci_valid = np.array([
                ci_covers_point(pt["delta"][i], bt["delta_lo"][i], bt["delta_hi"][i])
                if np.isfinite(bt["delta_lo"][i]) else None
                for i in range(n_f)], dtype=object)
            if (variant, kk) in s2res["md_floors"]:
                md_fl = s2res["md_floors"][(variant, kk)]
                md_perm = s2res["md_perm_store"][(variant, kk)]
                p_md = np.array([permutation_p(pt["md_change"][i], md_perm[i], greater=False)
                                 for i in range(n_f)], dtype=float)
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
                floor_p05=fl, margin_over_floor=fl - pt["delta"],
                p_delta=p, q_delta=q,
                delta_ci_lo=bt["delta_lo"], delta_ci_hi=bt["delta_hi"],
                delta_ci_valid=ci_valid,
                cos_ci_lo=bt["cos_lo"], cos_ci_hi=bt["cos_hi"],
                frac_ci_lo=bt["frac_lo"], frac_ci_hi=bt["frac_hi"],
                md_change=pt["md_change"],
                md_floor_p05=md_fl, md_margin_over_floor=md_fl - pt["md_change"],
                p_md=p_md, q_md=q_md,
                inside_radius=np.where(np.isfinite(pt["inside"]), pt["inside"] > 0.5, None),
                guide_median_cos=guide_cos[variant], guides_disagree=gdis,
                identity_drop=id_drop[variant], identity_loss=iloss,
                identity_panel_n=int(id_pos.size),
                ruler_score_change=ruler_chg[variant],
                ruler_score_drop=-ruler_chg[variant],
                toward_young=toward,
            ))
            out.append(df)
    tab = pd.concat(out, ignore_index=True)
    tab.to_csv(SOUTH2_DIR / "stage3_full_table.csv", index=False)
    log(f"[stage3] full table {len(tab):,} rows -> stage3_full_table.csv")

    cnt = (tab.groupby(["variant", "space", "pair"], sort=False)
           .agg(n_factors=("factor", "size"),
                n_delta_negative=("delta", lambda s: int((s < 0).sum())),
                n_below_floor=("margin_over_floor", lambda s: int((s > 0).sum())),
                n_q_le_005=("q_delta", lambda s: int((s <= Q_BAR).sum())),
                n_guides_disagree=("guides_disagree", "sum"),
                n_identity_loss=("identity_loss", "sum"),
                n_toward_young=("toward_young", "sum"),
                n_md_below_own_floor=("md_margin_over_floor", lambda s: int((s > 0).sum())),
                n_md_q_le_005=("q_md", lambda s: int((s <= Q_BAR).sum())),
                n_inside_radius=("inside_radius", lambda s: int(sum(bool(v) for v in s
                                                                    if v is not None))))
           .reset_index())
    cnt.to_csv(SOUTH2_DIR / "stage3_candidate_counts.csv", index=False)
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
    df.to_csv(SOUTH2_DIR / "stage4_disagreement.csv", index=False)
    for _, r in df.iterrows():
        log(f"[stage4] {r['variant']:6s} {r['space']:9s} {r['pair']:6s} "
            f"rho(delta,cos)={r['rho_delta_cos']:+.3f} rho(delta,ruler)={r['rho_delta_ruler']:+.3f} "
            f"rho(cos,ruler)={r['rho_cos_ruler']:+.3f} top20 overlap d/c="
            f"{r['overlap_top20_delta_cos']:2d} cos_top20_not_toward={r['n_top20_cos_not_toward_young']:2d}")
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
        f"(kind `{inv['id_type_records'].get('south2_var_gene_id', {}).get('kind')}`) and "
        f"`var/gene_name` a symbol column "
        f"(kind `{inv['id_type_records'].get('south2_var_gene_name', {}).get('kind')}`).",
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
        f"For each factor, its own vector's entries are permuted across the "
        f"{inv['n_matched']:,} matched gene coordinates {N_PERM} times (seed {SEED}, an "
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
        "Full dose-response deltas at every f are in `results/south2/stage2c_dose_response.csv`.",
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
        f"The full table is `results/south2/stage3_full_table.csv`: **{len(tab):,} rows** = "
        f"{inv['n_factors_ge2']:,} factors x {len(VARIANTS)} gene variants x "
        f"{len(tab.groupby(['space', 'pair'])):,} space/origin/target settings, with cos, frac, "
        "delta, ||d||, ||v||, the factor's own Stage 2A floor, margin over that floor, "
        "permutation p, BH q, bootstrap intervals, the Mahalanobis change and membership test "
        "for Y2, guide agreement, the identity panel, and the frozen ruler score change. "
        "Nothing is pooled across spaces, variants, origins, or targets.",
        "",
        f"- **Guide reliability**: each guide's own vector is built and the median pairwise "
        "cosine between guides of the same factor is reported per factor and per variant. "
        "`guides_disagree` fires when that median is <= 0. Those factors are reported and are "
        "never labelled a candidate.",
        f"- **Identity check**: the panel actually used is "
        f"{', '.join(inv['identity_panel'])} — **{len(inv['identity_panel'])} of the "
        f"{len(IDENTITY_GENES)} pre-registered genes**. Absent: "
        f"{', '.join(inv['identity_absent']) or 'none'}. "
        f"`identity_drop` is the negated mean displacement over that panel **in frozen-z units** "
        f"(the analysis's own coordinates), and `identity_loss` fires above {IDENTITY_DROP}.",
        "- **Frozen ruler score change** is a context column only. It gates nothing.",
        f"- **Nulls**: the Stage 2A per-factor permutations, one-sided empirical p "
        f"(P(perm delta <= observed delta)), then BH-FDR across factors **within each space and "
        "variant** (and within each origin/target pair, since a pair is part of the setting).",
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
            f"**Why each of the {len(zero)} zero settings is zero.** For "
            f"**{n_q} of {len(zero)}** the binding constraint is the q-value gate: no factor "
            "reaches q <= 0.05 there. That is a resolution limit of the pre-registered "
            f"n_perm = {N_PERM} (next paragraph), not a finding that those factors were tested "
            "and found null. No zero setting is zero because factors cleared every statistical "
            "bar and were then removed by a flag.",
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
    lines += [
        f"**The BH-FDR gate is resolution-limited by the pre-registered n_perm = {N_PERM}.** The "
        f"smallest p an empirical one-sided permutation test can return is 1/({N_PERM}+1) = "
        f"{p_min:.5f}. With BH across {n_f:,} factors, q <= {Q_BAR} is unreachable unless at "
        f"least **{need} factors are tied at that minimum p** in the same setting. So a setting "
        "with `q_le_005` = 0 is reporting that fewer than that many factors hit the p floor; it "
        "is not reporting that the remaining factors were tested and found null. Observed ties "
        "at the p floor, per setting, are in `stage3_full_table.csv`.",
        "",
        f"**`guides_disagree` in V_sig is almost always an exact zero, not a disagreement.** "
        f"{guide_cos_zero['n_zero']:,} of {n_f:,} factors have a median pairwise guide cosine of "
        f"**exactly 0** in V_sig, {guide_cos_zero['n_neg']:,} are negative, and "
        f"{guide_cos_zero['n_pos']:,} are positive. A cosine of exactly 0 here means two guides "
        "for the same factor share **no significant gene at all**, so their masked vectors are "
        "orthogonal by construction. The pre-registered rule flags `median <= 0`, so these fire. "
        "They are reported as fired and are not relabelled. The finding is that per-guide "
        f"`adj_p <= {ADJ_P_BAR}` masking leaves most guides of the same factor with disjoint "
        "support. This is **not the sole reason** V_sig has no candidate: "
        f"{int((cnt[cnt.variant == 'V_sig']['n_q_le_005'] == 0).sum())} of "
        f"{int((cnt.variant == 'V_sig').sum())} V_sig settings also have no factor at "
        "q <= 0.05, so the q-value gate alone already returns zero there, even for the "
        f"{guide_cos_zero['n_pos']:,} factors whose guides do agree.",
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
        f"- **The q-value gate is resolution-limited.** With n_perm = {N_PERM}, the smallest "
        f"attainable p is 1/{N_PERM + 1}, and BH over {inv['n_factors_ge2']:,} factors needs at "
        f"least {int(np.ceil(inv['n_factors_ge2'] / (N_PERM + 1) / Q_BAR))} factors tied at that "
        "minimum before any q can reach 0.05. Settings with no q <= 0.05 are therefore not "
        "evidence of absence, and the candidate counts depend partly on how many factors tie at "
        "the p floor. More permutations would change them; that would be a new pre-registration.",
        "- **The control-row floor is not a noise floor** and was never used as the gate: those "
        "rows are untargeted but not inert.",
        "- This task makes **no claim about the source paper's own conclusions**. It did not test "
        "them.",
        "",
    ]


def build_findings(prereg, inv=None, man=None, sections=None, result=None):
    fired = (man or {}).get("fired_keys") or []
    flags = [r["key"] for r in ((man or {}).get("report_flags") or [])]
    if fired:
        status = ("STOP. Keys fired: " + ", ".join(f"`{k}`" for k in fired)
                  + ". Stages after the stop did not run.")
    else:
        status = ("Ran to completion. No STOP key fired. Report-only flags: "
                  + (", ".join(f"`{k}`" for k in flags) or "none") + ".")
    head = [
        "# FINDINGS_SOUTH2 — rank Southard Hs27 CRISPRa transcription factors by approach to a "
        "young target, from the per-guide regression output",
        "",
        f"**Status:** {status}",
        "",
    ]
    if result:
        head += [
            f"**Result:** {result} Read the \"What produced these counts\" section before any "
            "count: several settings are zero for geometric or resolution reasons rather than "
            "because factors were tested and found null, and the O2 settings are dominated by a "
            "cross-dataset offset.",
            "",
        ]
    head += [
        "SOUTH2 supersedes SOUTH. SOUTH halted at Stage 0 with `no_raw_counts` because its "
        "pre-registration described this file as a count matrix; the file holds per-guide "
        "regression output. `FINDINGS_SOUTH.md` was not edited except for one line at its top "
        "marking it superseded.",
        "",
        "The pre-registration below was written verbatim to `FINDINGS_SOUTH2.md` and "
        "`results/south2/PREREG.flag`, and `results/south2/PREREG.flag` was created, before any "
        "statistic was computed. Nothing was re-tuned after numbers existed.",
        "",
        "## Pre-registration (verbatim, written before any SOUTH2 statistic)",
        "",
        "```",
        prereg,
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
        "# PROGRESS_SOUTH2",
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
    log = Logger(SOUTH2_DIR / "south2_run.log")
    log("=" * 100)
    log("SOUTH2 — Southard Hs27 CRISPRa TFs ranked by approach to a young fibroblast target")
    log("=" * 100)
    log(f"[seeds] SEED={SEED} (splits, nulls, permutations) BOOT_SEED={BOOT_SEED} (bootstrap); "
        f"imported TOWARD_SEED={TOWARD_SEED}, SAME BOOT_SEED={BOOT_SEED}")
    log(f"[params] n_null={N_NULL} n_perm={N_PERM} n_boot={N_BOOT} S2_KS={S2_KS} "
        f"(lowdim_common.KS={LOWDIM_KS}) F={F_VALUES} adj_p_bar={ADJ_P_BAR} q_bar={Q_BAR}")
    log("[params] n_null is recorded as pre-registered; the pre-registration assigns the null to "
        "the Stage 2A per-factor permutations (n_perm=200), so no separate null draw exists.")
    log("[scope] Uncalibrated R2 is never a gate; the frozen ruler is loaded, never refit; "
        "L4A_matrix.h5 is not opened; nothing is downloaded.")
    log("[scope] Only FINDINGS_SOUTH2.md, PROGRESS_SOUTH2.md, results/south2/* and "
        "src/south2_run.py are written by this run.")

    prereg = load_prereg_text()
    log(f"[prereg] read {PREREG_FLAG} ({len(prereg)} chars); not rewriting the block")

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
        build_findings(prereg, None, man, ["## Stage 0 inventory", "",
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
        build_findings(prereg, inv, man, stage0_lines(inv, man) + keys_lines(man)
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

    # ---------------- Stage 2 ---------------- #
    s2res, fired2 = stage2(pack, combos, Lm, D, Cctrl, tgt, log)
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
        build_findings(prereg, inv, man, sections)
        write_progress(man, "Halt as pre-registered on `no_detection_power`. Stages 3-4 did not "
                            "run. Report Stages 0-2 only: this data cannot support a ranking.")
        log("[STOP] no_detection_power")
        return

    # ---------------- Stages 3 and 4 ---------------- #
    tab, cnt, n_invalid = stage3(pack, combos, Lm, D, factors, rows_by_factor, s2res, log)
    s4 = stage4(tab, log)
    blocks, word = top20_tables(tab, inv["coverage_low"])

    man = load_manifest()
    man["status"] = "DONE"
    man["inventory_path"] = str(INVENTORY_PATH)
    man["n_stage3_rows"] = int(len(tab))
    man["toward_young_by_setting"] = {
        f"{r['variant']}|{r['space']}|{r['pair']}": int(r["n_toward_young"])
        for _, r in cnt.iterrows()}
    save_manifest(man)

    gs = tab[tab.variant == "V_sig"].drop_duplicates("factor")["guide_median_cos"]
    guide_cos_zero = dict(n_zero=int((gs == 0).sum()), n_neg=int((gs < 0).sum()),
                          n_pos=int((gs > 0).sum()))
    tw = cnt[["variant", "space", "pair", "n_factors", "n_toward_young"]]
    nonzero = tw[tw.n_toward_young > 0]
    # Reported per setting, never pooled into one number.
    summary = (
        f"`toward_young` is 0 in {len(tw) - len(nonzero)} of the {len(tw)} settings and non-zero "
        f"in {len(nonzero)}"
        + (f", ranging {int(nonzero.n_toward_young.min())}-{int(nonzero.n_toward_young.max())} "
           f"of {int(tw.n_factors.iloc[0]) if 'n_factors' in tw.columns else 0} factors"
           if len(nonzero) else "")
        + f". Every non-zero setting is {', '.join(sorted(set(nonzero.variant))) or 'none'}; "
        f"{', '.join(v for v in VARIANTS if v not in set(nonzero.variant)) or 'no variant'} "
        "returns 0 everywhere. The counts are not pooled and no single setting is primary.")
    log(f"[summary] {summary}")
    sections += (stage3_lines(tab, cnt, blocks, word, n_invalid, inv)
                 + reading_lines(tab, cnt, inv, tgt, dist_df, guide_cos_zero)
                 + stage4_lines(s4)
                 + keys_lines(man)
                 + limitations_lines(inv, s2res, tgt))
    build_findings(prereg, inv, man, sections, result=summary)
    write_progress(
        man,
        "Nothing further is required by this pre-registration: Stages 0-4 all ran and are "
        "reported in FINDINGS_SOUTH2.md. If this is taken further, the two things that would "
        "change the answer rather than decorate it are (a) an input whose units are stated, "
        "since `scale_unknown` fired and the MDA is not scale-invariant, and (b) a "
        "perturbation dataset measured in aged cells, since every delta here is a "
        "counterfactual built from effects measured in one neonatal line. Both need a new "
        "pre-registration.",
        [f"**Result in one line:** {summary}"
         + (" `coverage_low` fired, so every candidate is labelled exploratory."
            if inv["coverage_low"] else "")],
    )
    log("[done]")


if __name__ == "__main__":
    main()

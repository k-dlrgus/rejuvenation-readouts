"""SOUTH — rank the Southard Hs27 CRISPRa transcription factors by approach to a
young target, after first measuring how small an approach this data could detect.

Reads results/survey/opened/southard/fibroblast_CRISPRa_mean_pop.h5ad (read-only).
Reuses src/toward_run.py, src/same_run.py, src/seng_run.py, src/lowdim_common.py
without changing them. Does not refit the frozen ruler.

The pre-registration must already exist verbatim at results/south/PREREG.flag,
written before any statistic. This script does not author that block; it reads it.

Stage 0 (inventory, descriptive only) runs first. Every Stage 0 STOP key halts the
task: Stages 1-4 do not run and no statistic is computed.

Usage: python src/south_run.py
"""
from __future__ import annotations

import re
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import RESULTS, ROOT  # noqa: E402
from brain_phase1_common import Logger, dump_json, load_json  # noqa: E402
from fibro_common import FROZEN_RULER, jsonable  # noqa: E402
from fibro2_common import load_frozen_ruler  # noqa: E402
from gtex_common import StopStep  # noqa: E402
from target_common import md_table  # noqa: E402

# Required imports from the four named modules (do not change them). Stages 1-4
# would use these; the seeds and space definitions are pinned here for provenance.
from toward_run import TOWARD_SEED  # noqa: E402
from same_run import BOOT_SEED as SAME_BOOT_SEED  # noqa: E402
from lowdim_common import KS as LOWDIM_KS, SpaceFit  # noqa: E402
from seng_run import (  # noqa: E402
    IDENTITY_GENES,
    dump_tree,
    file_md5,
    locate_count_matrix,
    read_obs_var,
    ruler_column_index,
)

SOUTH_DIR = RESULTS / "south"
SOUTH_DIR.mkdir(parents=True, exist_ok=True)
FINDINGS_PATH = ROOT / "FINDINGS_SOUTH.md"
PROGRESS_PATH = ROOT / "PROGRESS_SOUTH.md"
PREREG_FLAG = SOUTH_DIR / "PREREG.flag"
MANIFEST_PATH = SOUTH_DIR / "manifest.json"
INVENTORY_PATH = SOUTH_DIR / "inventory.json"

SOUTHARD_DIR = RESULTS / "survey" / "opened" / "southard"
H5_PATH = SOUTHARD_DIR / "fibroblast_CRISPRa_mean_pop.h5ad"

SEED = 20260914
BOOT_SEED = 20260918
N_NULL = 200
N_PERM = 200
N_BOOT = 200
N_RULER_EXPECTED = 23485
MIN_CONTROLS = 10
MIN_REPS = 2
S2_KS = (20, 50, 100)
F_VALUES = (0.0, 0.05, 0.1, 0.25, 0.5, 1.0)

TARGET_COL = "target_gene"
GUIDE_COL = "guide_identity"
CELL_COUNT_COL = "cell_count"
# Anchored so it flags genuine non-targeting/control labels only, not real TF
# genes that merely contain the substrings (ARNT, NONO, MNT, VENTX, SCRT1/2, ...).
# `off-target` is this library's own untargeted arm: those rows carry an empty
# target_gene_id, a NaN target_expr, active=False, and the label is absent from
# the source's own 1,836-entry TF target list (hs27_targets.txt). It is a control
# label, not a transcription factor, and is counted as such.
CONTROL_RE = re.compile(
    r"(?i)(non[-_ ]?targeting|non[-_ ]?target|^off[-_ ]?target|scrambled?|scrambl|intergenic|"
    r"safe[-_ ]?harbor|^nt[_-]|^ntc$|^ntc[_-]|^control$|^ctrl$|^ctrl[_-]|negative[-_ ]?control)"
)
# Columns that evidence whether a control-labelled row really is untargeted.
CONTROL_EVIDENCE_COLS = (
    "target_gene_id", "target_expr", "active", "masked_active", "gene_driven",
    "sequence_driven", "de_genes", "strength", "cell_count",
)
# Prereg-order Stage 0 STOP keys.
STAGE0_KEYS = ("no_raw_counts", "no_resampling_unit", "no_controls")


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
    return dict(seed=SEED, boot_seed=BOOT_SEED, failures=[], fired_keys=[], status="INIT")


def save_manifest(man):
    dump_json(MANIFEST_PATH, jsonable(man))
    return MANIFEST_PATH


def record_failure(step, message, details=None):
    man = load_manifest()
    rec = dict(step=step, message=str(message), details=jsonable(details or {}))
    man.setdefault("failures", []).append(rec)
    if step in STAGE0_KEYS and step not in man.setdefault("fired_keys", []):
        man["fired_keys"].append(step)
    man["status"] = "STOP"
    save_manifest(man)
    return rec


# --------------------------------------------------------------------------- #
# Stage 0 inventory (descriptive only; no statistic)
# --------------------------------------------------------------------------- #
def row_unit(obs, n_rows):
    guides = obs[GUIDE_COL].astype(str)
    targets = obs[TARGET_COL].astype(str)
    one_row_per_guide = bool(guides.nunique() == n_rows)
    rows_per_target = targets.value_counts()
    per_tf_mean = bool(rows_per_target.max() == 1)
    has_barcodes = any(re.search(r"(?i)barcode", str(c)) for c in obs.columns)
    if per_tf_mean:
        unit = "per-factor mean"
    elif has_barcodes:
        unit = "single cell"
    else:
        unit = "replicate (per-guide mean-population row)"
    return dict(
        unit=unit,
        one_row_per_guide=one_row_per_guide,
        n_unique_guides=int(guides.nunique()),
        max_rows_per_target=int(rows_per_target.max()),
        min_rows_per_target=int(rows_per_target.min()),
        cell_barcode_column=bool(has_barcodes),
        cell_count_column=CELL_COUNT_COL in obs.columns,
    )


def sibling_audit(log):
    """Record what the other already-downloaded Southard files are, so the reason a
    raw-count matrix cannot be obtained from them is a checked fact, not an assertion.

    The pre-registration pins the analysed table; nothing here is substituted for it.
    """
    out = []
    for p in sorted(SOUTHARD_DIR.glob("*")):
        if p.name == H5_PATH.name or p.suffix in (".txt", ".bin", ".deflate"):
            continue
        rec = dict(name=p.name, size_bytes=int(p.stat().st_size), datasets=[])
        if p.suffix in (".h5", ".h5ad"):
            try:
                with h5py.File(p, "r") as f:
                    def walk(g, prefix=""):
                        for k in g.keys():
                            o = g[k]
                            if isinstance(o, h5py.Dataset):
                                rec["datasets"].append(
                                    dict(path=f"{prefix}{k}", shape=list(o.shape), dtype=str(o.dtype)))
                            elif len(prefix.split("/")) < 4:
                                walk(o, f"{prefix}{k}/")
                    walk(f)
                    # CellRanger-style matrix: does it hold raw integer gene counts?
                    if "matrix" in f and "data" in f.get("matrix", {}):
                        m = f["matrix"]
                        d = m["data"]
                        shape = [int(x) for x in m["shape"][:]] if "shape" in m else None
                        ftypes = {}
                        if "features" in m and "feature_type" in m["features"]:
                            vals, cnts = np.unique(m["features"]["feature_type"][:], return_counts=True)
                            ftypes = {v.decode() if isinstance(v, bytes) else str(v): int(c)
                                      for v, c in zip(vals, cnts)}
                        rec["cellranger"] = dict(
                            shape_features_barcodes=shape,
                            n_barcodes=int(m["barcodes"].shape[0]) if "barcodes" in m else None,
                            feature_types=ftypes,
                            data_dtype=str(d.dtype),
                            integer_dtype=bool(np.issubdtype(d.dtype, np.integer)),
                            data_min=int(d[:].min()), data_max=int(d[:].max()),
                            n_gene_expression=int(ftypes.get("Gene Expression", 0)),
                        )
            except OSError as e:
                rec["error"] = str(e)
        out.append(rec)
        log(f"[sibling] {rec}")
    return out


def stage0(log):
    for p in (H5_PATH, FROZEN_RULER):
        if not Path(p).exists():
            raise StopStep("inventory", f"missing {p}. Not substituting.")

    t0 = time.perf_counter()
    size = int(H5_PATH.stat().st_size)
    md5 = file_md5(H5_PATH)
    log(f"[file] {H5_PATH} size={size} md5={md5} ({time.perf_counter() - t0:.1f}s)")

    with h5py.File(H5_PATH, "r") as f:
        tree = dump_tree(f, SOUTH_DIR / "h5_tree.txt")
        log(f"[h5] tree written ({len(tree)} lines) -> {SOUTH_DIR / 'h5_tree.txt'}")
        obs = read_obs_var(f["obs"])
        var = read_obs_var(f["var"])
        x_shape = tuple(int(s) for s in f["X"].shape) if isinstance(f["X"], h5py.Dataset) else None
        x_dtype = str(f["X"].dtype) if isinstance(f["X"], h5py.Dataset) else None
        layer_names = list(f["layers"].keys()) if "layers" in f else []
        uns_keys = list(f["uns"].keys()) if "uns" in f else []
        raw_group = "raw" in f
        log(f"[columns obs] {list(obs.columns)}")
        log(f"[columns var] {list(var.columns)}")
        log(f"[h5] X shape={x_shape} dtype={x_dtype} layers={layer_names} "
            f"raw_group={raw_group} uns={uns_keys}")
        for col in (TARGET_COL, GUIDE_COL):
            if col not in obs.columns:
                raise StopStep("inventory", f"obs missing {col!r}. Not substituting.",
                               dict(columns=list(obs.columns)))

        # Raw-count audit (same rule as seng_run.locate_count_matrix).
        raw_stop = None
        try:
            _, _, source, audits = locate_count_matrix(f)
            raw_rec = dict(raw_integer_counts=True, source=source, audits=audits)
        except StopStep as e:
            if e.step != "no_raw_counts":
                raise
            raw_stop = e
            raw_rec = dict(raw_integer_counts=False, source=None,
                           audits=(e.details or {}).get("audits", []))
    log(f"[raw] raw_integer_counts={raw_rec['raw_integer_counts']} source={raw_rec['source']}")
    for a in raw_rec["audits"]:
        log(f"[raw audit] {a}")
    if raw_rec["audits"]:
        pd.DataFrame(raw_rec["audits"]).to_csv(SOUTH_DIR / "inventory_raw_audit.csv", index=False)

    n_rows = int(len(obs))
    unit = row_unit(obs, n_rows)
    log(f"[unit] {unit}")

    targets = obs[TARGET_COL].astype(str)
    guides = obs[GUIDE_COL].astype(str)
    cc = (obs[CELL_COUNT_COL].astype(float) if CELL_COUNT_COL in obs.columns
          else pd.Series(np.nan, index=obs.index))
    cc_all = cc.to_numpy(float) if CELL_COUNT_COL in obs.columns else None
    ctrl_target = targets.map(lambda s: bool(CONTROL_RE.search(s))).to_numpy()
    ctrl_guide = guides.map(lambda s: bool(CONTROL_RE.search(s))).to_numpy()
    ctrl_mask = ctrl_target | ctrl_guide
    n_controls = int(ctrl_mask.sum())
    control_labels = sorted(set(targets[ctrl_mask]))
    log(f"[controls] rows whose {TARGET_COL} or {GUIDE_COL} matches {CONTROL_RE.pattern!r}: "
        f"{n_controls} labels={control_labels}")
    ctrl_tab = pd.DataFrame({
        "row": np.arange(n_rows)[ctrl_mask],
        TARGET_COL: targets[ctrl_mask].to_numpy(),
        GUIDE_COL: guides[ctrl_mask].to_numpy(),
    })
    for c in CONTROL_EVIDENCE_COLS:
        if c in obs.columns:
            ctrl_tab[c] = obs[c].to_numpy()[ctrl_mask]
    ctrl_tab.to_csv(SOUTH_DIR / "inventory_controls.csv", index=False)

    # Evidence that the control-labelled rows are genuinely untargeted, and that
    # they are not inert: some off-target guides still move the transcriptome.
    ctrl_evidence = {}
    if n_controls:
        gid = obs["target_gene_id"].astype(str).str.strip().to_numpy()[ctrl_mask] \
            if "target_gene_id" in obs.columns else np.array([], dtype=object)
        te = pd.to_numeric(obs["target_expr"], errors="coerce").to_numpy(float)[ctrl_mask] \
            if "target_expr" in obs.columns else np.array([np.nan])
        ctrl_evidence = dict(
            n_empty_target_gene_id=int((gid == "").sum()) if gid.size else None,
            n_nan_target_expr=int(np.isnan(te).sum()),
            n_active=int(np.asarray(obs["active"].to_numpy()[ctrl_mask], dtype=bool).sum())
            if "active" in obs.columns else None,
            n_gene_driven=int(np.asarray(obs["gene_driven"].to_numpy()[ctrl_mask], dtype=bool).sum())
            if "gene_driven" in obs.columns else None,
            n_sequence_driven=int(np.asarray(obs["sequence_driven"].to_numpy()[ctrl_mask], dtype=bool).sum())
            if "sequence_driven" in obs.columns else None,
            n_with_de_genes=int((pd.to_numeric(obs["de_genes"], errors="coerce")
                                 .to_numpy(float)[ctrl_mask] > 0).sum())
            if "de_genes" in obs.columns else None,
            cell_count_min=float(np.nanmin(cc_all[ctrl_mask])) if cc_all is not None else None,
            cell_count_median=float(np.nanmedian(cc_all[ctrl_mask])) if cc_all is not None else None,
            cell_count_max=float(np.nanmax(cc_all[ctrl_mask])) if cc_all is not None else None,
        )
        log(f"[controls evidence] {ctrl_evidence}")

    tf_labels = sorted(set(targets) - set(control_labels))
    rows_tab = pd.DataFrame({
        "row": np.arange(n_rows),
        GUIDE_COL: guides.to_numpy(),
        TARGET_COL: targets.to_numpy(),
        CELL_COUNT_COL: cc.to_numpy(),
        "control_label_match": ctrl_mask,
    })
    rows_tab.to_csv(SOUTH_DIR / "inventory_rows.csv", index=False)

    tgt = (
        rows_tab.groupby(TARGET_COL)
        .agg(n_reps=(GUIDE_COL, "size"), cell_count_sum=(CELL_COUNT_COL, "sum"))
        .reset_index()
        .sort_values("n_reps", ascending=False)
    )
    tgt.to_csv(SOUTH_DIR / "inventory_targets.csv", index=False)
    reps_per = tgt.loc[tgt[TARGET_COL].isin(tf_labels), "n_reps"]
    rp_dist = reps_per.value_counts().sort_index()
    rp_df = pd.DataFrame({"n_reps": rp_dist.index.astype(int), "n_factors": rp_dist.to_numpy(int)})
    rp_df.to_csv(SOUTH_DIR / "inventory_reps_per_factor.csv", index=False)
    log(f"[reps per factor] {dict(zip(rp_df.n_reps, rp_df.n_factors))}")
    n_factors_ge_min = int((reps_per >= MIN_REPS).sum())

    non_symbol = [t for t in tf_labels if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9.\-]*", t)]
    log(f"[targets] non-gene-symbol-shaped labels: {non_symbol}")
    # Hyphenated labels are audited explicitly because the only non-gene label in
    # this file (`off-target`) is hyphenated and would otherwise pass as a symbol.
    hyphen_labels = sorted(t for t in set(targets) if "-" in t)
    hyphen_tf_labels = [t for t in hyphen_labels if t in set(tf_labels)]
    log(f"[targets] hyphenated labels ({len(hyphen_labels)}): {hyphen_labels}")

    # Frozen-ruler overlap (Ensembl first, symbol fallback; ruler is not refit).
    frozen = load_frozen_ruler()
    n_ruler = int(len(np.asarray(frozen["mu"])))
    gid_col = "gene_id" if "gene_id" in var.columns else None
    sym_col = "gene_name" if "gene_name" in var.columns else None
    if gid_col is None and sym_col is None:
        raise StopStep("inventory",
                       f"var has neither gene_id nor gene_name. columns={list(var.columns)}")
    idx = ruler_column_index(
        var[gid_col].astype(str).to_numpy() if gid_col else var[sym_col].astype(str).to_numpy(),
        var[sym_col].astype(str).to_numpy() if sym_col else var[gid_col].astype(str).to_numpy(),
        frozen, log,
    )
    n_overlap = int((idx >= 0).sum())
    ruler_sym = np.array([str(s).upper() for s in np.asarray(frozen["symbol"])], dtype=object)
    present_sym = set(ruler_sym[idx >= 0])
    id_present = [g for g in IDENTITY_GENES if g in present_sym]
    id_missing = [g for g in IDENTITY_GENES if g not in present_sym]
    log(f"[ruler] frozen-ruler genes={n_ruler} (expected {N_RULER_EXPECTED}) present={n_overlap} "
        f"of {len(var)} var genes; identity present={id_present} missing={id_missing}")

    # Resampling unit.
    if unit["cell_barcode_column"]:
        rep_unit = "cells"
    elif n_factors_ge_min > 0:
        rep_unit = "replicates"
    else:
        rep_unit = "none"
    log(f"[resampling] unit={rep_unit} factors with >= {MIN_REPS} replicate rows={n_factors_ge_min}")

    inv = dict(
        path=str(H5_PATH), size_bytes=size, md5=md5,
        x_shape=list(x_shape) if x_shape else None, x_dtype=x_dtype,
        layers=layer_names, raw_group=raw_group, uns_keys=uns_keys,
        obs_columns=list(obs.columns), var_columns=list(var.columns),
        n_rows=n_rows, n_var_genes=int(len(var)),
        unit=unit,
        n_target_labels=int(targets.nunique()), n_factors=len(tf_labels),
        non_symbol_labels=non_symbol,
        hyphen_labels=hyphen_labels, hyphen_tf_labels=hyphen_tf_labels,
        reps_per_factor={int(k): int(v) for k, v in zip(rp_df.n_reps, rp_df.n_factors)},
        n_factors_ge_min_reps=n_factors_ge_min, min_reps=MIN_REPS,
        n_controls=n_controls, control_labels=control_labels, control_regex=CONTROL_RE.pattern,
        min_controls=MIN_CONTROLS, control_evidence=ctrl_evidence,
        cell_count_min=float(np.nanmin(cc)) if CELL_COUNT_COL in obs.columns else None,
        cell_count_median=float(np.nanmedian(cc)) if CELL_COUNT_COL in obs.columns else None,
        cell_count_max=float(np.nanmax(cc)) if CELL_COUNT_COL in obs.columns else None,
        cell_count_sum=float(np.nansum(cc)) if CELL_COUNT_COL in obs.columns else None,
        raw=raw_rec,
        n_ruler=n_ruler, n_ruler_expected=N_RULER_EXPECTED, n_ruler_present=n_overlap,
        identity_present=id_present, identity_missing=id_missing,
        resampling_unit=rep_unit,
        siblings=sibling_audit(log),
    )
    dump_json(INVENTORY_PATH, jsonable(inv))
    log(f"[inventory] wrote {INVENTORY_PATH}")

    # Evaluate STOP conditions in pre-registration order.
    fired = []
    if raw_stop is not None:
        fired.append(("no_raw_counts", raw_stop.message, dict(audits=raw_rec["audits"])))
    if rep_unit == "none":
        fired.append((
            "no_resampling_unit",
            "No replication unit (cells or replicate rows) is available for resampling. "
            "Reporting the inventory only.",
            dict(unit=unit, n_factors_ge_min_reps=n_factors_ge_min),
        ))
    if n_controls < MIN_CONTROLS:
        fired.append((
            "no_controls",
            f"{n_controls} untargeted-control rows exist (< {MIN_CONTROLS}). No {TARGET_COL} or "
            f"{GUIDE_COL} value matches {CONTROL_RE.pattern!r}; every one of the {len(tf_labels)} "
            "labels is a targeted transcription factor. They are reported, not relabelled as controls.",
            dict(n_controls=n_controls, control_labels=control_labels),
        ))
    return inv, fired


# --------------------------------------------------------------------------- #
# FINDINGS rendering
# --------------------------------------------------------------------------- #
def _fmt_int(x):
    try:
        return f"{int(x):,}"
    except (TypeError, ValueError):
        return "NA"


def inventory_lines(inv):
    unit = inv["unit"]
    raw = inv["raw"]
    lines = [
        "## Stage 0 inventory (descriptive only; no statistic)",
        "",
        f"- File: `{inv['path']}`",
        f"- Size: {_fmt_int(inv['size_bytes'])} bytes. md5 `{inv['md5']}`.",
        f"- X: dense {inv['x_shape'][0]:,} rows x {inv['x_shape'][1]:,} genes, dtype "
        f"`{inv['x_dtype']}`. Layers: {', '.join(f'`{l}`' for l in inv['layers']) or 'none'}. "
        f"`raw` group: {inv['raw_group']}. `uns` keys: {inv['uns_keys'] or 'none'}.",
        f"- obs columns ({len(inv['obs_columns'])}): {', '.join(f'`{c}`' for c in inv['obs_columns'])}.",
        f"- var columns ({len(inv['var_columns'])}): {', '.join(f'`{c}`' for c in inv['var_columns'])}.",
        "",
        "### What one row is",
        "",
        f"**{unit['unit']}.** `{GUIDE_COL}` is unique on every row "
        f"({unit['n_unique_guides']:,} values for {inv['n_rows']:,} rows), while a targeted "
        f"transcription factor has between {min(inv['reps_per_factor'])} and "
        f"{max(inv['reps_per_factor'])} rows (the {inv['n_controls']}-row `off-target` group is the "
        "untargeted control arm, not a factor). So a row is not a "
        f"per-factor mean (multiple rows per factor) and not a single cell (no barcode column; only a "
        f"per-row `{CELL_COUNT_COL}`: min {inv['cell_count_min']:.0f}, median "
        f"{inv['cell_count_median']:.0f}, max {inv['cell_count_max']:.0f}, total "
        f"{inv['cell_count_sum']:,.0f}). Each row is the guide-level 'mean population' produced by the "
        "source's regression model, i.e. a per-guide replicate summary, not raw expression.",
        "",
        "### Transcription factors, replicate rows, controls",
        "",
        f"- Distinct `{TARGET_COL}` labels: {inv['n_target_labels']:,}. Of these, "
        f"{inv['n_factors']:,} are targeted transcription factors and "
        f"{len(inv['control_labels'])} is an untargeted control label.",
        f"- Untargeted / non-targeting control rows: **{inv['n_controls']}** "
        f"(threshold for `no_controls` is {inv['min_controls']}). Control labels: "
        + (", ".join(f"`{x}`" for x in inv["control_labels"]) or "none")
        + f". Matched on `{TARGET_COL}` or `{GUIDE_COL}` against `{inv['control_regex']}`.",
    ]
    ev = inv.get("control_evidence") or {}
    if ev:
        lines.append(
            f"- Evidence those {inv['n_controls']} rows are genuinely untargeted: "
            f"{ev.get('n_empty_target_gene_id')}/{inv['n_controls']} carry an empty `target_gene_id`, "
            f"{ev.get('n_nan_target_expr')}/{inv['n_controls']} a NaN `target_expr` (there is no intended "
            f"target whose activation could be measured), {ev.get('n_active')}/{inv['n_controls']} are "
            f"`active`, and {ev.get('n_gene_driven')}/{inv['n_controls']} are `gene_driven`. The label is "
            "also absent from the source's own 1,836-entry TF target list "
            "(`hs27_targets.txt`), which equals the 1,837 labels here minus this one."
        )
        lines.append(
            f"- Those control rows are **not inert**: {ev.get('n_sequence_driven')} of "
            f"{inv['n_controls']} are flagged `sequence_driven` and {ev.get('n_with_de_genes')} have "
            f"`de_genes` > 0. They also carry far more cells per row than the targeted guides "
            f"(median {ev.get('cell_count_median'):.0f}, range {ev.get('cell_count_min'):.0f}-"
            f"{ev.get('cell_count_max'):.0f}, against a median of "
            f"{inv['cell_count_median']:.0f} over all rows). A control floor built from them would "
            "therefore absorb real sequence-driven effects, which would raise the floor, not lower it."
        )
    if inv["non_symbol_labels"]:
        lines.append(f"- Labels that are not gene-symbol-shaped (reported, not relabelled): "
                     f"{', '.join(f'`{x}`' for x in inv['non_symbol_labels'])}.")
    if inv.get("hyphen_labels"):
        lines.append(
            f"- Hyphenated labels ({len(inv['hyphen_labels'])}), audited explicitly because a hyphenated "
            f"label passes a gene-symbol shape test: "
            + ", ".join(f"`{x}`" for x in inv["hyphen_labels"])
            + f". {len(inv.get('hyphen_tf_labels') or [])} of them are standard HGNC symbols kept as "
            "transcription factors; the remainder are the control label above."
        )
    lines += [
        "",
        "Replicate rows per transcription factor:",
        "",
        md_table(pd.DataFrame({
            "n_replicate_rows": list(inv["reps_per_factor"].keys()),
            "n_factors": list(inv["reps_per_factor"].values()),
        }).astype(str)),
        "",
        f"{inv['n_factors_ge_min_reps']:,} factors have at least {inv['min_reps']} replicate rows.",
        "",
        "### Raw counts",
        "",
        f"- X holds raw integer counts: **{raw['raw_integer_counts']}**. A raw layer with integer "
        f"counts exists: **{bool(raw['source'])}**.",
        "",
        "Audit of every matrix checked (same rule as `seng_run.locate_count_matrix`):",
        "",
    ]
    if raw["audits"]:
        adf = pd.DataFrame(raw["audits"])
        keep = [c for c in ("where", "dtype", "n_seen", "n_fractional", "n_negative",
                            "n_nonfinite", "min", "max", "raw_integer_counts") if c in adf.columns]
        adf = adf[keep].copy()
        for c in ("n_seen", "n_fractional", "n_negative", "n_nonfinite"):
            if c in adf.columns:
                adf[c] = [_fmt_int(v) for v in adf[c]]
        for c in ("min", "max"):
            if c in adf.columns:
                adf[c] = [f"{float(v):.4g}" if v is not None and np.isfinite(float(v)) else "NA"
                          for v in adf[c]]
        lines += [md_table(adf.astype(str)), ""]
    if inv.get("siblings"):
        lines += [
            "The other already-downloaded Southard files in the same folder were checked, so that what "
            "is and is not available locally is a verified fact rather than an assumption. None of them "
            "was substituted for the pinned table:",
            "",
        ]
        for s in inv["siblings"]:
            cr = s.get("cellranger")
            if cr:
                ft = ", ".join(f"{k} {v:,}" for k, v in sorted(cr["feature_types"].items()))
                desc = (f"CellRanger matrix, {cr['shape_features_barcodes'][0]:,} features x "
                        f"{cr['shape_features_barcodes'][1]:,} barcodes ({ft}); `matrix/data` dtype "
                        f"`{cr['data_dtype']}`, integer dtype {cr['integer_dtype']}, range "
                        f"{cr['data_min']}-{cr['data_max']:,}")
            else:
                desc = "; ".join(f"`{d['path']}` {tuple(d['shape'])} {d['dtype']}"
                                 for d in s.get("datasets", [])) or "no datasets readable"
            lines.append(f"- `{s['name']}` ({_fmt_int(s['size_bytes'])} bytes): {desc}."
                         + (f" Read error: {s['error']}" if s.get("error") else ""))
        lines += [
            "",
            "`fibroblast_CRISPRa_aggr_total_guide_umis.h5` cannot supply counts: it is a sparse COO "
            "table whose two index levels are cell barcodes and guides, i.e. guide-assignment UMIs, "
            "with no gene axis.",
            "",
        ]
        raw_sib = [s for s in inv["siblings"]
                   if (s.get("cellranger") or {}).get("integer_dtype")
                   and (s.get("cellranger") or {}).get("n_gene_expression")]
        for s in raw_sib:
            cr = s["cellranger"]
            pct = (100.0 * cr["n_barcodes"] / inv["cell_count_sum"]) if inv.get("cell_count_sum") else None
            lines += [
                f"**`{s['name']}` does hold raw integer gene counts** "
                f"({cr['n_gene_expression']:,} Gene Expression features x {cr['n_barcodes']:,} cells, "
                f"`{cr['data_dtype']}`, minimum {cr['data_min']}). Stating this plainly because it is the "
                "one fact that bears on what could be done next, and because it would be wrong to report "
                "that no counts for this screen exist locally. It was **not** used, for two reasons that "
                "are separate:",
                "",
                "1. The pre-registration pins the analysed table to "
                f"`{Path(inv['path']).name}` under DATA, and applies the `no_raw_counts` test to that "
                "table. Swapping in a different file after the pinned one failed its own test would be "
                "choosing the data after seeing the result, which this task forbids.",
                (f"2. It is a single lane. Its {cr['n_barcodes']:,} cells are "
                 f"{pct:.1f}% of the {inv['cell_count_sum']:,.0f} cells the pinned table summarises, "
                 f"spread over {cr['feature_types'].get('CRISPR Guide Capture', 0):,} guide features - "
                 f"on the order of {cr['n_barcodes'] / max(cr['feature_types'].get('CRISPR Guide Capture', 1), 1):.1f} "
                 "cells per guide. Even as a future path it would not by itself support a "
                 f"{inv['n_factors']:,}-factor ranking; the remaining lanes would have to be obtained."
                 if pct else
                 "2. It is a single lane of the screen and would not by itself support a full ranking."),
                "",
            ]
    lines += [
        "### Frozen-ruler genes present",
        "",
        f"- {inv['n_ruler_present']:,} of {inv['n_ruler']:,} frozen-ruler genes (expected "
        f"{inv['n_ruler_expected']:,}) are present among the file's {inv['n_var_genes']:,} genes "
        "(Ensembl id first, symbol fallback; `seng_run.ruler_column_index`). The ruler is not refit.",
        f"- Identity genes present: {', '.join(inv['identity_present']) or 'none'}. "
        f"Missing: {', '.join(inv['identity_missing']) or 'none'}.",
        "",
        "### Resampling unit",
        "",
        f"- **{inv['resampling_unit']}.** Cells are not in this file (only a per-row cell count), so "
        f"cell-level resampling is impossible; {inv['n_factors_ge_min_reps']:,} factors carry "
        f">= {inv['min_reps']} replicate rows, so replicate-level resampling would be the available unit.",
        "",
    ]
    return lines


def correction_lines(inv):
    """Record this task's own first-pass Stage 0 error verbatim, as required."""
    return [
        "### Recorded failure of this task's own first Stage 0 pass",
        "",
        "The first pass of SOUTH reported `0` untargeted control rows and `1,837` transcription "
        "factors, and fired `no_controls` on that basis. That was wrong. Verbatim, it reported:",
        "",
        "> - Distinct `target_gene` labels: 1,837, all targeted transcription factors (1,837 after "
        "removing control-matched labels, of which there are none).",
        "> - Untargeted / non-targeting control rows: **0**.",
        "",
        "The error was in the control-label pattern, not in the data: it tested for "
        "`non-targeting`, `scrambled`, `intergenic`, `safe-harbor`, `NTC`, `ctrl` and similar, but not "
        "for `off-target`, which is the name this library gives its untargeted arm. The 17 "
        "`off-target` rows were therefore counted as a 1,837th transcription factor — the one factor "
        "the first pass reported as having 17 replicate rows, where every real factor has 2 to 6.",
        "",
        f"Corrected here: {inv['n_factors']:,} transcription factors and {inv['n_controls']} control "
        "rows, so `no_controls` does **not** fire. The `off-target` rows were not dropped or "
        "relabelled to make anything look better; they are reported with the evidence above.",
        "",
        "This correction does not change the task's outcome. `no_raw_counts` fires independently and "
        "halts SOUTH at Stage 0 either way. It changes which keys are reported as fired, and it "
        "removes a false claim that this dataset carries no controls.",
        "",
    ]


def spaces_lines(inv):
    ks = ", ".join(str(k) for k in S2_KS)
    return [
        "## Spaces and their sizes",
        "",
        f"- S1 (gene space): the {inv['n_ruler']:,} frozen-ruler overlap genes, of which "
        f"{inv['n_ruler_present']:,} are present in this file; missing genes are z=0 by construction. "
        "S1 was **not built**: its pipeline (sum -> TMM log2-CPM with prior.count=2 -> frozen GTEx "
        "mu/sd -> z) requires raw counts, and this file has none.",
        f"- S2 (low-dim space): PCA fit only on the GTEx fibroblast donor matrix at k in ({ks}), applied "
        "unchanged elsewhere. **Not fitted**: the task halted at Stage 0 before any space was built.",
        "",
    ]


def targets_lines():
    return [
        "## Targets and origins",
        "",
        "Not built. Y1 (GTEx fibroblast donors 20-39), Y2 (the same donors as a Mahalanobis region in "
        "S2, with a held-out 95th-percentile radius), O1 (GTEx fibroblast donors 60-79), O2 (the Lu "
        "et al. GM00731 aged day-0 pseudobulk as built in SAME), and every origin-to-target distance "
        "are downstream of the Stage 0 STOP and were not computed.",
        "",
    ]


def stages_not_run_lines():
    return [
        "## Stage 1 detection limit (control floor, dose-response, MDA)",
        "",
        "Not run. Stage 0 fired `no_raw_counts`, which halts the task before any statistic. The row "
        "count needed for Stage 1A exists (17 control rows, enough to split into two disjoint halves "
        "200 times), but the quantity being split is z, and z is defined here as sum -> TMM log2-CPM "
        "with prior.count=2 -> frozen GTEx mu/sd. That pipeline starts from raw counts, which this file "
        "does not contain. No floor, no dose-response, no MDA, and no `no_detection_power` "
        "determination were computed.",
        "",
        "## Stage 2 single-factor ranking",
        "",
        "Not run. Stage 1 did not pass (it did not run), and the per-factor effect vector "
        "d_g = z(g rows) - z(control rows) cannot be formed, because z cannot be formed without raw "
        "counts. The count of `toward_young` factors is therefore undefined, not zero.",
        "",
        "## Stage 3 disagreement with directional readouts",
        "",
        "Not run. It ranks factors by delta, cosine, and ruler-score drop, all of which require "
        "Stage 2 outputs.",
        "",
        "## Stage 4 combinations",
        "",
        "**Additivity warning (stated as required):** this stage would assume that perturbation effects "
        "add linearly, and this data cannot test that assumption. Not run, because Stage 1 did not pass.",
        "",
    ]


def stop_lines(man):
    lines = ["## Fired STOP keys", ""]
    fired = man.get("fired_keys") or []
    if fired:
        lines.append("Fired, in pre-registration order: " + ", ".join(f"`{k}`" for k in fired) + ".")
        not_fired = [k for k in STAGE0_KEYS if k not in fired]
        if not_fired:
            lines.append("Checked and not fired: " + ", ".join(f"`{k}`" for k in not_fired) + ".")
        lines += ["", "Failures recorded verbatim:", ""]
        for fr in man.get("failures") or []:
            lines.append(f"- **{fr.get('step')}:** {fr.get('message')}")
    else:
        lines.append("None.")
    lines.append("")
    return lines


def limitations_lines(inv=None):
    ev = (inv or {}).get("control_evidence") or {}
    n_ctrl = (inv or {}).get("n_controls")
    ruler_pct = (100.0 * inv["n_ruler_present"] / inv["n_ruler"]) if inv else None
    cc_ratio = (ev.get("cell_count_median") / inv["cell_count_median"]
                if ev.get("cell_count_median") and inv and inv.get("cell_count_median") else None)
    return [
        "## Limitations",
        "",
        "- Hs27 is one neonatal fibroblast line with no aged cells, so every delta this task could have "
        "produced is a counterfactual that assumes the perturbation effect transfers to aged donors.",
        "- Effects here are measured in a different cell line and platform (Hs27 CRISPRa screen) from the "
        "origin and target (GTEx / Lu et al. donor fibroblasts), so any comparison crosses a batch and "
        "platform boundary.",
        "- Combinations (Stage 4) assume additivity, which this data cannot test.",
        "- The rows are per-guide 'mean population' summaries from a regression model, not raw counts "
        "and not single cells, so cell-level resampling is not possible and the frozen z-pipeline "
        "cannot be applied at all.",
        "- Decisive blocker this run: the matrix is not raw integer counts and no raw layer exists "
        "(`no_raw_counts`). That halts the ranking on its own; no delta, floor, or MDA exists to "
        "report, and none was invented.",
        (f"- Even had counts existed, only {inv['n_ruler_present']:,} of the {inv['n_ruler']:,} "
         f"frozen-ruler genes ({ruler_pct:.0f}%) are present in this file, so {100 - ruler_pct:.0f}% of "
         "the S1 coordinates would have been set to z=0 by the missing-gene rule. Every S1 distance "
         "would then be dominated by a fixed, perturbation-independent offset. This is recorded as a "
         "limitation, not used as a gate."
         if inv else
         "- Frozen-ruler coverage in this file is partial; missing genes would enter S1 as z=0."),
        (f"- The {n_ctrl} control rows are untargeted but not inert "
         f"({ev.get('n_sequence_driven')} are flagged `sequence_driven` and "
         f"{ev.get('n_with_de_genes')} have `de_genes` > 0), and they carry about "
         f"{cc_ratio:.0f}x more cells per row than the median row. A floor built from them would be "
         "both noisier and less exchangeable with a real factor's rows than the pre-registration's "
         "row-count matching alone assumes."
         if ev and cc_ratio else
         "- Control-row exchangeability with real factor rows was not established."),
        "",
    ]


def build_findings(prereg, inv, man):
    head = [
        "# FINDINGS_SOUTH — rank Southard Hs27 CRISPRa transcription factors by approach to a young target",
        "",
        "**Status:** STOP at Stage 0. Keys fired: "
        + (", ".join(f"`{k}`" for k in (man.get('fired_keys') or [])) or "none")
        + ". No statistic was computed; Stages 1-4 did not run.",
        "",
        "The pre-registration below was written verbatim to `FINDINGS_SOUTH.md` and "
        "`results/south/PREREG.flag` before any Southard row was opened and before any statistic.",
        "",
        "## Pre-registration (verbatim, written before any SOUTH statistic)",
        "",
        "```",
        prereg,
        "```",
        "",
    ]
    body = []
    if inv is not None:
        body += inventory_lines(inv)
        body += correction_lines(inv)
        body += spaces_lines(inv)
        body += targets_lines()
        body += stages_not_run_lines()
    body += stop_lines(man)
    body += limitations_lines(inv)
    FINDINGS_PATH.write_text("\n".join(head + body).rstrip() + "\n", encoding="utf-8")
    return FINDINGS_PATH


def write_progress(man, next_action):
    fired = man.get("fired_keys") or []
    status = ("STOP at Stage 0; keys fired: " + ", ".join(f"`{k}`" for k in fired)) if fired \
        else "running"
    text = "\n".join([
        "# PROGRESS_SOUTH",
        "",
        f"**STOP status:** {status}.",
        "",
        f"**Next action:** {next_action}",
        "",
    ]) + "\n"
    PROGRESS_PATH.write_text(text, encoding="utf-8")
    return PROGRESS_PATH


# --------------------------------------------------------------------------- #
def main():
    log = Logger(SOUTH_DIR / "south_run.log")
    log(f"[seeds] SEED={SEED} BOOT_SEED={BOOT_SEED} "
        f"(toward TOWARD_SEED={TOWARD_SEED}, same BOOT_SEED={SAME_BOOT_SEED}, lowdim KS={LOWDIM_KS})")
    log(f"[params] n_null={N_NULL} n_perm={N_PERM} n_boot={N_BOOT} S2_KS={S2_KS} F={F_VALUES}")

    prereg = load_prereg_text()
    log(f"[prereg] read {PREREG_FLAG} ({len(prereg)} chars)")

    # A run describes itself: clear any keys and failures left by an earlier run so
    # nothing stale is reported as having fired this time.
    man = load_manifest()
    man.update(seed=SEED, boot_seed=BOOT_SEED, status="INIT", fired_keys=[], failures=[])
    save_manifest(man)

    try:
        inv, fired = stage0(log)
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        man = load_manifest()
        build_findings(prereg, None, man)
        write_progress(man, "Resolve the Stage 0 structural failure above; nothing else can run until "
                            "the inventory completes.")
        log(f"[STOP] {e.step}: {e.message}")
        return

    for key, message, details in fired:
        record_failure(key, message, details)
    man = load_manifest()
    man["inventory"] = jsonable(inv)
    save_manifest(man)

    build_findings(prereg, inv, man)
    x_audit = next((a for a in inv["raw"]["audits"] if a.get("where") == "X"), {})
    x_neg_pct = (100.0 * x_audit.get("n_negative", 0) / x_audit["n_seen"]) if x_audit.get("n_seen") else float("nan")
    x_min = float(x_audit.get("min", float("nan")))
    if man.get("fired_keys"):
        write_progress(
            man,
            "Halt as pre-registered on `no_raw_counts`. This file holds per-guide 'mean population' "
            f"regression output (X is {inv['x_dtype']}, {x_neg_pct:.0f}% of its entries negative, min "
            f"{x_min:.2f}) with no raw layer and no counts layer, so the frozen z-pipeline cannot start "
            "and neither the floor nor any delta can be formed. The "
            f"{inv['n_controls']} `off-target` control rows do exist, so `no_controls` does not "
            "fire and a floor would be constructible from a counts table. Reported rather than acted "
            "on: `L4A_matrix.h5`, in the same folder, does hold raw integer gene counts (36,601 genes x "
            "18,197 cells), but it is not the table this pre-registration pins under DATA, and it is one "
            "lane holding 2.4% of the screen's cells at roughly 1.7 cells per guide. Swapping it in now "
            "would be choosing the data after seeing the pinned table fail. Stages 1-4 would need a new "
            "pre-registration naming a count-level table and the remaining lanes; that is a new-data "
            "decision outside this task.",
        )
        log(f"[STOP] Stage 0 fired: {[k for k, _, _ in fired]}")
    else:
        write_progress(man, "Stage 0 passed; proceed to Stage 1 detection limit.")
        log("[stage0] passed; no STOP fired")


if __name__ == "__main__":
    main()

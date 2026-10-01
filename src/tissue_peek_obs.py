"""Inspect obs (and, if present, per-cell n_counts) of a CELLxGENE h5ad WITHOUT loading X.

Used to decide whether batch/donor structure is recoverable before committing to a multi-GB
download. Prefers a local file; otherwise opens the published HTTPS URL with fsspec range
requests (the CXG CDN supports HTTP Range). Every field printed is read from the file.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tissue_common import TISSUE_DIR, Logger, BATCH_COL_CANDIDATES  # noqa: E402
from search_cxg import parse_age_years  # noqa: E402


def _open_h5(path_or_url: str):
    import h5py
    if path_or_url.startswith("http://") or path_or_url.startswith("https://"):
        import fsspec
        of = fsspec.open(path_or_url, "rb", block_size=2 * 1024 * 1024)
        fobj = of.open()
        return h5py.File(fobj, "r"), of
    return __import__("h5py").File(path_or_url, "r"), None


def _decode(v):
    if isinstance(v, bytes):
        return v.decode("utf-8", "replace")
    return v


def read_obs_column(h, name):
    """Read one obs column from anndata 0.8+ HDF5 (categorical codes + categories, or a raw dataset)."""
    g = h["obs"]
    if name not in g:
        return None
    node = g[name]
    if hasattr(node, "dtype"):
        arr = node[()]
        if getattr(arr, "dtype", None) is not None and arr.dtype.kind in ("S", "O"):
            return np.array([_decode(x) for x in arr], dtype=object)
        return arr
    # group: categorical
    if "codes" in node and "categories" in node:
        codes = node["codes"][()]
        cats = [_decode(x) for x in node["categories"][()]]
        out = np.empty(len(codes), dtype=object)
        for i, c in enumerate(codes):
            out[i] = cats[int(c)] if int(c) >= 0 else None
        return out
    return None


def obs_column_names(h):
    cols = []
    if "obs" not in h:
        return cols
    g = h["obs"]
    # anndata stores _index plus one dataset/group per column; also a column-order attribute
    order = g.attrs.get("column-order")
    if order is not None:
        return [_decode(x) for x in order]
    skip = {"_index", "index"}
    for k in g.keys():
        if k not in skip:
            cols.append(k)
    return cols


def peek(path_or_url: str, name: str, log, max_levels=12):
    log(f"\n{'=' * 90}\nPEEK {name}\n  source: {path_or_url}\n{'=' * 90}")
    h, handle = _open_h5(path_or_url)
    try:
        shape = None
        for p in ("X", "raw/X"):
            if p in h:
                sh = h[p].attrs.get("shape")
                if sh is not None:
                    shape = tuple(int(x) for x in sh)
                    log(f"  matrix {p} shape={shape}")
                    break
        cols = obs_column_names(h)
        log(f"  obs columns ({len(cols)}): {cols}")
        n = None
        # prefer donor_id length
        wanted = [c for c in (
            "donor_id", "development_stage", "sex", "cell_type", "assay", "tissue",
            "suspension_type", "disease", "is_primary_data",
            "n_counts", "nCount_RNA", "n_genes", "nFeature_RNA", "total_counts",
            "n_genes_by_counts",
        ) + BATCH_COL_CANDIDATES if c in cols]
        # unique-preserve
        seen, uniq = set(), []
        for c in wanted:
            if c not in seen:
                seen.add(c); uniq.append(c)
        summary = {"name": name, "source": path_or_url, "n_cells": None if shape is None else shape[0],
                   "n_genes": None if shape is None else shape[1], "obs_columns": cols}
        for c in uniq:
            arr = read_obs_column(h, c)
            if arr is None:
                continue
            n = len(arr)
            vc = pd.Series(arr).value_counts(dropna=False)
            summary[f"nunique_{c}"] = int(vc.shape[0])
            log(f"  {c}: n={n:,}  nunique={vc.shape[0]}  top={list(vc.head(max_levels).items())}")
        # age distribution at donor level
        if "donor_id" in cols and "development_stage" in cols:
            donor = read_obs_column(h, "donor_id")
            stage = read_obs_column(h, "development_stage")
            ages = pd.Series([parse_age_years(s) for s in stage], dtype=float)
            d = pd.DataFrame({"donor": donor, "age": ages})
            if "assay" in cols:
                d["assay"] = read_obs_column(h, "assay")
            if "suspension_type" in cols:
                d["suspension_type"] = read_obs_column(h, "suspension_type")
            if "tissue" in cols:
                d["tissue"] = read_obs_column(h, "tissue")
            don = d.groupby("donor", dropna=False).agg(age=("age", "first"), n_cells=("age", "size"))
            log(f"  donors: {len(don)}  age min={don.age.min()} max={don.age.max()} "
                f"nunique_age={don.age.nunique()}  cells/donor median={don.n_cells.median():.0f} "
                f"(min {don.n_cells.min()}, max {don.n_cells.max()})")
            summary.update(n_donors=int(len(don)), age_min=None if don.age.isna().all() else float(don.age.min()),
                           age_max=None if don.age.isna().all() else float(don.age.max()),
                           nunique_age=int(don.age.nunique(dropna=True)))
        # batch recoverability
        batch_present = [c for c in BATCH_COL_CANDIDATES if c in cols and c not in
                         ("donor_id", "development_stage")]
        log(f"  candidate batch-like columns present: {batch_present or '(NONE)'}")
        summary["batch_columns_present"] = batch_present
        summary["batch_recoverable"] = bool(batch_present)
        return summary
    finally:
        h.close()
        if handle is not None:
            handle.close()


def main(targets=None):
    """targets: list of dicts with keys name, url or path, optional dataset_id."""
    log = Logger(TISSUE_DIR / "t1_obs_peek.txt")
    if targets is None:
        spec = TISSUE_DIR / "t1_peek_targets.json"
        if not spec.exists():
            log("no t1_peek_targets.json — nothing to peek")
            log.close()
            return []
        targets = json.loads(spec.read_text(encoding="utf-8"))
    summaries = []
    for t in targets:
        src = t.get("path") or t.get("url")
        try:
            summaries.append(peek(src, t["name"], log))
        except Exception as e:
            log(f"\nPEEK FAILED {t.get('name')}: {type(e).__name__}: {e}")
            summaries.append(dict(name=t.get("name"), error=str(e)))
    pd.DataFrame(summaries).to_json(TISSUE_DIR / "t1_obs_peek.json", orient="records", indent=1)
    # also a flat csv of the scalar fields
    flat = []
    for s in summaries:
        row = {k: v for k, v in s.items() if not isinstance(v, (list, dict))}
        row["obs_columns"] = ";".join(s.get("obs_columns") or [])
        row["batch_columns_present"] = ";".join(s.get("batch_columns_present") or [])
        flat.append(row)
    pd.DataFrame(flat).to_csv(TISSUE_DIR / "t1_obs_peek.csv", index=False)
    log.close()
    return summaries


if __name__ == "__main__":
    main()

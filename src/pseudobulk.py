"""Build donor x cell-type pseudobulks from a large CELLxGENE h5ad without loading X into RAM.

Every filtering decision is logged with counts (HARD RULE 3). Output:
  data/processed/<name>_pseudobulk.h5ad   (X = summed raw UMI counts, float32)
  results/tables/<name>_filter_log.csv
"""
import re
import sys
import time
import numpy as np
import pandas as pd
import scipy.sparse as sp
import h5py
import anndata as ad

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from config import DATA_PROC, TAB, MIN_CELLS_PER_PSEUDOBULK, MIN_DONORS_PER_CELLTYPE  # noqa: E402
from gene_lists import CELL_CYCLE_ALL  # noqa: E402

FILTER_LOG = []


def log_filter(step, before, after, unit, note=""):
    removed = before - after
    FILTER_LOG.append(dict(step=step, unit=unit, before=before, after=after, removed=removed, note=note))
    print(f"[filter] {step:<55} {unit:<10} {before:>10,} -> {after:>10,}  (removed {removed:,}) {note}", flush=True)


def parse_age(stage: str):
    m = re.match(r"^(\d+)[- ]year[- ]old", str(stage).lower())
    return float(m.group(1)) if m else np.nan


def read_counts_chunk(f: h5py.File, path: str, r0: int, r1: int, n_var: int) -> sp.csr_matrix:
    """Read rows [r0, r1) of a CSR-encoded h5ad matrix group."""
    g = f[path]
    indptr = g["indptr"][r0:r1 + 1]
    s, e = indptr[0], indptr[-1]
    data = g["data"][s:e]
    indices = g["indices"][s:e]
    return sp.csr_matrix((data, indices, indptr - s), shape=(r1 - r0, n_var))


def find_counts_path(f: h5py.File) -> str:
    """Prefer raw/X (CELLxGENE convention: raw counts), else X. Check integer-ness on a sample."""
    for p in ("raw/X", "X"):
        if p in f and isinstance(f[p], h5py.Group):
            enc = f[p].attrs.get("encoding-type", b"")
            enc = enc.decode() if isinstance(enc, bytes) else str(enc)
            if enc != "csr_matrix":
                continue
            d = f[p]["data"][:200000]
            frac_int = np.mean(np.abs(d - np.round(d)) < 1e-6)
            print(f"[counts] candidate {p}: encoding={enc} dtype={d.dtype} integer-fraction(sample)={frac_int:.4f} max={d.max():.1f}")
            if frac_int > 0.999:
                return p
    raise RuntimeError("no integer CSR count matrix found in raw/X or X")


def build(h5ad_path, name, donor_col="donor_id", ct_col="cell_type", age_col="development_stage", sex_col="sex",
          extra_obs=(), chunk=40000, min_cells=MIN_CELLS_PER_PSEUDOBULK, min_donors=MIN_DONORS_PER_CELLTYPE,
          celltype_map=None, require_primary=False, keep_mask=None):
    t0 = time.time()
    a = ad.read_h5ad(h5ad_path, backed="r")
    obs = a.obs.copy()
    obs["_pos"] = np.arange(len(obs))  # integer row position; barcodes may repeat across pools
    var = a.raw.var.copy() if a.raw is not None else a.var.copy()
    n_cells0, n_genes0 = obs.shape[0], var.shape[0]
    print(f"[load] {name}: {n_cells0:,} cells x {n_genes0:,} genes; obs index unique={obs.index.is_unique}; "
          f"obs cols={list(obs.columns)}", flush=True)

    if keep_mask is not None:
        keep_mask = np.asarray(keep_mask, dtype=bool)
        if keep_mask.shape[0] != n_cells0:
            raise ValueError(f"keep_mask length {keep_mask.shape[0]} != n_cells {n_cells0}")
        n = len(obs)
        obs = obs[keep_mask]
        log_filter("caller keep_mask (analysis cohort)", n, len(obs), "cells")

    if age_col == "age":
        obs["age"] = pd.to_numeric(obs[age_col], errors="coerce")
    else:
        obs["age"] = obs[age_col].map(parse_age).astype(float)
    if obs["age"].isna().any() and "donor_age" in obs.columns:
        alt = pd.to_numeric(obs["donor_age"], errors="coerce")
        fill = obs["age"].isna() & alt.notna()
        obs.loc[fill, "age"] = alt[fill]
        print(f"[age] filled {int(fill.sum()):,} cells from donor_age where {age_col} was non-numeric", flush=True)
    obs["donor"] = obs[donor_col].astype(str)
    obs["celltype_raw"] = obs[ct_col].astype(str)
    obs["celltype"] = obs["celltype_raw"].map(celltype_map) if celltype_map else obs["celltype_raw"]
    obs["sex"] = obs[sex_col].astype(str) if sex_col in obs else "unknown"

    # ---- cell-level filters (logged) ----
    n = len(obs)
    keep = obs["age"].notna()
    log_filter("cells with numeric donor age", n, int(keep.sum()), "cells", "non-numeric development_stage dropped")
    obs = obs[keep]
    if require_primary and "is_primary_data" in obs.columns:
        n = len(obs)
        keep = obs["is_primary_data"].astype(bool)
        log_filter("is_primary_data==True", n, int(keep.sum()), "cells")
        obs = obs[keep]
    if celltype_map:
        n = len(obs)
        keep = obs["celltype"].notna()
        log_filter("cells whose cell type is in the analysis map", n, int(keep.sum()), "cells",
                   f"dropped raw labels: {sorted(set(obs.loc[~keep, 'celltype_raw']))[:15]}")
        obs = obs[keep]

    # (donor, celltype) group sizes
    grp = obs.groupby(["donor", "celltype"], observed=True).size().rename("n_cells").reset_index()
    n_grp0 = len(grp)
    grp_keep = grp[grp.n_cells >= min_cells]
    log_filter(f"pseudobulk groups with >= {min_cells} cells", n_grp0, len(grp_keep), "groups")
    cells_before = len(obs)
    key = obs["donor"] + "||" + obs["celltype"]
    keep_keys = set(grp_keep.donor + "||" + grp_keep.celltype)
    obs = obs[key.isin(keep_keys)]
    log_filter("cells in retained pseudobulk groups", cells_before, len(obs), "cells")

    # cell types present in enough donors
    ct_donors = grp_keep.groupby("celltype").donor.nunique()
    ct_keep = ct_donors[ct_donors >= min_donors].index
    log_filter(f"cell types present in >= {min_donors} donors", len(ct_donors), len(ct_keep), "cell types",
               f"dropped: {dict(ct_donors[ct_donors < min_donors])}")
    cells_before = len(obs)
    obs = obs[obs.celltype.isin(ct_keep)]
    log_filter("cells in retained cell types", cells_before, len(obs), "cells")
    grp_keep = grp_keep[grp_keep.celltype.isin(ct_keep)].reset_index(drop=True)
    grp_keep["key"] = grp_keep.donor + "||" + grp_keep.celltype
    gidx = {k: i for i, k in enumerate(grp_keep.key)}
    n_grp = len(grp_keep)

    # ---- aggregate counts in chunks ----
    obs = obs.sort_values("_pos")
    row_pos = obs["_pos"].to_numpy()  # positions in the full matrix (sorted ascending)
    cell_key = (obs["donor"] + "||" + obs["celltype"]).map(gidx).to_numpy()
    assert len(row_pos) == len(obs) and np.all(np.diff(row_pos) > 0)
    sym = var["feature_name"].astype(str).to_numpy() if "feature_name" in var else var.index.to_numpy()
    cc_mask = np.isin(sym, CELL_CYCLE_ALL)
    prolif_markers = np.isin(sym, ["MKI67", "TOP2A", "CENPF", "NUSAP1", "STMN1"])

    X_sum = np.zeros((n_grp, n_genes0), dtype=np.float64)
    tot_counts = np.zeros(n_grp)
    tot_genes = np.zeros(n_grp)
    tot_cc = np.zeros(n_grp)
    n_prolif = np.zeros(n_grp)
    with h5py.File(h5ad_path, "r") as f:
        cpath = find_counts_path(f)
        n_var = f[cpath].attrs["shape"][1]
        assert n_var == n_genes0
        # iterate over contiguous blocks of the full matrix, subset rows we keep
        pos_ptr = 0
        for r0 in range(0, n_cells0, chunk):
            r1 = min(r0 + chunk, n_cells0)
            # rows of interest in this block
            j0 = pos_ptr
            while pos_ptr < len(row_pos) and row_pos[pos_ptr] < r1:
                pos_ptr += 1
            if pos_ptr == j0:
                continue
            sel = row_pos[j0:pos_ptr] - r0
            Xc = read_counts_chunk(f, cpath, r0, r1, n_var)[sel]
            gk = cell_key[j0:pos_ptr]
            M = sp.csr_matrix((np.ones(len(gk)), (gk, np.arange(len(gk)))), shape=(n_grp, len(gk)))
            X_sum += (M @ Xc).toarray()
            per_cell_tot = np.asarray(Xc.sum(1)).ravel()
            per_cell_ng = np.diff(Xc.indptr)
            per_cell_cc = np.asarray(Xc[:, cc_mask].sum(1)).ravel()
            per_cell_pm = (np.asarray(Xc[:, prolif_markers].sum(1)).ravel() > 0).astype(float)
            tot_counts += M @ per_cell_tot
            tot_genes += M @ per_cell_ng
            tot_cc += M @ per_cell_cc
            n_prolif += M @ per_cell_pm
            if (r0 // chunk) % 5 == 0:
                print(f"   rows {r1:,}/{n_cells0:,}  ({time.time()-t0:.0f}s)", flush=True)

    grp_keep["mean_counts_per_cell"] = tot_counts / grp_keep.n_cells
    grp_keep["mean_genes_per_cell"] = tot_genes / grp_keep.n_cells
    grp_keep["cc_umi_frac"] = tot_cc / tot_counts
    grp_keep["prolif_marker_cell_frac"] = n_prolif / grp_keep.n_cells
    donor_meta = obs.groupby("donor", observed=True).agg(age=("age", "first"), sex=("sex", "first"),
                                                          **{c: (c, "first") for c in extra_obs if c in obs})
    grp_keep = grp_keep.merge(donor_meta, left_on="donor", right_index=True, how="left")
    grp_keep.index = grp_keep.key.values

    pb = ad.AnnData(X=X_sum.astype(np.float32), obs=grp_keep.drop(columns=["key"]), var=var)
    pb.uns["source_h5ad"] = str(h5ad_path)
    pb.uns["n_cells_source"] = int(n_cells0)
    pb.uns["n_cells_used"] = int(len(obs))
    out = DATA_PROC / f"{name}_pseudobulk.h5ad"
    pb.write_h5ad(out, compression="gzip")
    pd.DataFrame(FILTER_LOG).to_csv(TAB / f"{name}_filter_log.csv", index=False)
    print(f"[done] pseudobulk {pb.shape} -> {out}  ({time.time()-t0:.0f}s)", flush=True)
    return pb

"""Project frozen DLPFC directions onto SEA-AD donor pseudobulks.

Primary: no/low AD pathology donors. Secondary: all donors, pathology covariate.
Does not run unless frozen_directions.npz and PREREG_20260917.flag exist.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import MIN_CELLS_PER_PSEUDOBULK  # noqa: E402
from external_common import (  # noqa: E402
    EXT_DIR, EXT_PROC, EXT_RAW, EXT_SEED, BOOT_SEED, N_PERM, N_BOOT,
    TRANSFER_NULL_BAR, MIN_CELLS_PB, MIN_N_R, PREREG_FLAG,
    StopStep, Logger, dump_json, load_json, jsonable, ext_log_banner,
    strip_ensembl, pearson_safe, spearman_safe, pred_scores, pass_r_bar,
    median_over_types, refresh_progress, finite, resolve_column,
)
from brain_phase1_common import r2_mae  # noqa: E402
from trajectory_common import residualize_identity, permutation_p, summarize_null_col  # noqa: E402
from tissue_peek_obs import _open_h5, obs_column_names, read_obs_column  # noqa: E402
from pseudobulk import find_counts_path, read_counts_chunk  # noqa: E402


def _norm_logcpm(counts: np.ndarray) -> np.ndarray:
    """Same transform as decomp.normalize_logcpm: log2(CPM+1) from summed UMI. No gene filter."""
    X = np.asarray(counts, np.float64)
    lib = X.sum(1, keepdims=True)
    if np.any(lib <= 0):
        bad = int((lib.ravel() <= 0).sum())
        raise StopStep("logcpm", f"{bad} pseudobulks have library size 0")
    cpm = X / lib * 1e6
    return np.log2(cpm + 1.0)


def residualize_1d(y, z):
    y = np.asarray(y, float)
    z = np.asarray(z, float)
    m = np.isfinite(y) & np.isfinite(z)
    out = np.full_like(y, np.nan, dtype=float)
    if int(m.sum()) < 3 or float(np.nanstd(z[m])) < 1e-12:
        return out
    A = np.column_stack([np.ones(int(m.sum())), z[m]])
    coef, *_ = np.linalg.lstsq(A, y[m], rcond=None)
    out[m] = y[m] - (coef[0] + coef[1] * z[m])
    return out


def read_var_ids(h, group: str) -> np.ndarray:
    if group not in h:
        raise StopStep("var", f"h5ad has no {group}/")
    g = h[group]
    node = g["_index"] if "_index" in g else None
    if node is None:
        raise StopStep("var", f"{group} has no _index. keys={list(g.keys())[:20]}")
    arr = node[()]
    if getattr(arr, "dtype", None) is not None and arr.dtype.kind in ("S", "O"):
        return np.array([x.decode("utf-8", "replace") if isinstance(x, (bytes, np.bytes_)) else str(x) for x in arr], dtype=object)
    return np.array([str(x) for x in arr], dtype=object)


def local_h5ad_if_complete(url: str, expected_bytes) -> str | None:
    fname = url.rstrip("/").rsplit("/", 1)[-1]
    path = EXT_RAW / fname
    if path.exists() and expected_bytes is not None and path.stat().st_size == int(expected_bytes):
        return str(path)
    if path.exists() and expected_bytes is None:
        return str(path)
    return None


def build_seaad_pseudobulk(h5ad_src: str, mapping: pd.DataFrame, log, chunk=40000):
    """Donor × mapped DLPFC-type summed UMI. min_cells = project constant. No min-donors drop."""
    if int(MIN_CELLS_PB) != int(MIN_CELLS_PER_PSEUDOBULK):
        raise StopStep(
            "pseudobulk",
            f"MIN_CELLS_PB={MIN_CELLS_PB} != config MIN_CELLS_PER_PSEUDOBULK={MIN_CELLS_PER_PSEUDOBULK}",
        )
    sub_to = dict(
        zip(mapping.loc[mapping.mapped, "seaad_subclass"].astype(str),
            mapping.loc[mapping.mapped, "dlpfc_type"].astype(str))
    )
    if not sub_to:
        raise StopStep("pseudobulk", "mapping has zero mapped subclasses")
    log(f"[pb] opening {h5ad_src}")
    h, handle = _open_h5(h5ad_src)
    t0 = time.time()
    try:
        cols = obs_column_names(h)
        if "donor_id" not in cols:
            raise StopStep("pseudobulk", f"obs missing donor_id. columns={cols}")
        subclass_col = resolve_column(
            cols, ("subclass", "Subclass"),
            field="subclass", step="pseudobulk", file_tag="seaad_h5ad_obs",
        )
        donor = read_obs_column(h, "donor_id")
        subclass = np.asarray(read_obs_column(h, subclass_col)).astype(str)
        n_cells0 = len(donor)
        dlpfc = pd.Series(subclass).map(sub_to)
        keep = dlpfc.notna().to_numpy()
        log(f"[pb] nuclei mapped by subclass: {int(keep.sum()):,} / {n_cells0:,}")
        obs = pd.DataFrame({
            "donor": np.asarray(donor).astype(str),
            "celltype": dlpfc.astype(object),
            "_pos": np.arange(n_cells0),
        })
        obs = obs.loc[keep].copy()
        grp = obs.groupby(["donor", "celltype"], observed=True).size().rename("n_cells").reset_index()
        n0 = len(grp)
        grp_keep = grp[grp.n_cells >= int(MIN_CELLS_PB)].copy()
        log(f"[pb] groups with ≥{MIN_CELLS_PB} cells: {n0} → {len(grp_keep)}")
        if grp_keep.empty:
            raise StopStep("pseudobulk", f"zero donor×type groups with ≥{MIN_CELLS_PB} cells")
        key = obs["donor"] + "||" + obs["celltype"].astype(str)
        keep_keys = set(grp_keep.donor.astype(str) + "||" + grp_keep.celltype.astype(str))
        obs = obs[key.isin(keep_keys)].copy()
        grp_keep = grp_keep.reset_index(drop=True)
        grp_keep["key"] = grp_keep.donor.astype(str) + "||" + grp_keep.celltype.astype(str)
        gidx = {k: i for i, k in enumerate(grp_keep.key)}
        n_grp = len(grp_keep)
        obs = obs.sort_values("_pos")
        row_pos = obs["_pos"].to_numpy()
        cell_key = (obs["donor"].astype(str) + "||" + obs["celltype"].astype(str)).map(gidx).to_numpy()
        if not np.all(np.diff(row_pos) > 0):
            raise StopStep("pseudobulk", "row positions not strictly increasing after sort")

        try:
            cpath = find_counts_path(h)
        except RuntimeError as e:
            raise StopStep("pseudobulk", str(e))
        n_var = int(h[cpath].attrs["shape"][1])
        var_group = "raw/var" if cpath.startswith("raw/") else "var"
        gene_id = read_var_ids(h, var_group)
        if len(gene_id) != n_var:
            raise StopStep(
                "pseudobulk",
                f"var length {len(gene_id)} != matrix n_var {n_var} (path={cpath}, var={var_group})",
            )
        log(f"[pb] counts path={cpath}  n_var={n_var:,}  groups={n_grp}")
        X_sum = np.zeros((n_grp, n_var), dtype=np.float64)
        pos_ptr = 0
        for r0 in range(0, n_cells0, chunk):
            r1 = min(r0 + chunk, n_cells0)
            j0 = pos_ptr
            while pos_ptr < len(row_pos) and row_pos[pos_ptr] < r1:
                pos_ptr += 1
            if pos_ptr == j0:
                continue
            sel = row_pos[j0:pos_ptr] - r0
            Xc = read_counts_chunk(h, cpath, r0, r1, n_var)[sel]
            gk = cell_key[j0:pos_ptr]
            M = sp.csr_matrix((np.ones(len(gk)), (gk, np.arange(len(gk)))), shape=(n_grp, len(gk)))
            X_sum += (M @ Xc).toarray()
            log(f"   rows {r1:,}/{n_cells0:,}  ({time.time()-t0:.0f}s)")
    finally:
        h.close()
        if handle is not None:
            handle.close()
    log(f"[pb] summed in {time.time()-t0:.0f}s  shape={X_sum.shape}")
    return grp_keep.drop(columns=["key"]).reset_index(drop=True), X_sum, gene_id


def align_genes(Y_sea, sea_ids, frozen, log):
    dlpfc_ids = [strip_ensembl(g) for g in np.asarray(frozen["gene_id"]).astype(str)]
    sea_ids_st = [strip_ensembl(g) for g in np.asarray(sea_ids).astype(str)]
    sea_pos = {}
    for i, g in enumerate(sea_ids_st):
        sea_pos.setdefault(g, i)
    overlap_idx_d = []
    overlap_idx_s = []
    missing = 0
    for j, g in enumerate(dlpfc_ids):
        i = sea_pos.get(g)
        if i is None:
            missing += 1
            continue
        overlap_idx_d.append(j)
        overlap_idx_s.append(i)
    n_ov = len(overlap_idx_d)
    log(f"[genes] DLPFC={len(dlpfc_ids)}  SEA-AD={len(sea_ids_st)}  overlap={n_ov}  "
        f"DLPFC-missing-in-SEA-AD={missing}  frac={n_ov / max(len(dlpfc_ids), 1):.4f}")
    if n_ov < 1:
        raise StopStep(
            "genes",
            "gene overlap is 0 after stripping Ensembl versions. Not switching to symbols.",
            dict(n_dlpfc=len(dlpfc_ids), n_sea=len(sea_ids_st)),
        )
    Y_al = np.zeros((Y_sea.shape[0], len(dlpfc_ids)), dtype=np.float64)
    Y_al[:, overlap_idx_d] = np.asarray(Y_sea, np.float64)[:, overlap_idx_s]
    rec = dict(
        n_dlpfc=len(dlpfc_ids), n_seaad=len(sea_ids_st), n_overlap=n_ov,
        n_dlpfc_missing=missing, frac_overlap=float(n_ov / len(dlpfc_ids)),
    )
    return Y_al, rec


def load_frozen():
    path = EXT_DIR / "frozen_directions.npz"
    if not path.exists():
        raise StopStep("project", f"missing {path}")
    return np.load(path, allow_pickle=True)


def scores_for_direction(Xz, W, intercept, slope, types, obs, kind: str):
    """DLPFC 1-D calibrated predictions per pseudobulk."""
    ct = obs.celltype.astype(str).to_numpy()
    pred = np.full(len(obs), np.nan)
    score = np.full(len(obs), np.nan)
    rows = []
    for j, t in enumerate(types):
        w = np.asarray(W[:, j], float)
        m = ct == t
        n = int(m.sum())
        ok = n >= 1 and np.isfinite(w).all() and float(np.linalg.norm(w)) > 1e-12
        if not ok:
            rows.append(dict(kind=kind, celltype=t, n=n, ok=False))
            continue
        s = Xz[m] @ w
        score[m] = s
        a, b = float(intercept[j]), float(slope[j])
        if np.isfinite(a) and np.isfinite(b):
            pred[m] = a + b * s
        else:
            pred[m] = s
        rows.append(dict(kind=kind, celltype=t, n=n, ok=True, intercept=a, slope=b))
    return pred, score, pd.DataFrame(rows)


def _type_vectors(obs, pred, donors_keep=None, path_score=None):
    """Per-type donor-level (pred, age, pathology). One row per donor (age is donor-constant)."""
    df = obs.copy()
    df["pred"] = pred
    if donors_keep is not None:
        keep = set(map(str, donors_keep))
        df = df[df.donor.astype(str).isin(keep)]
    out = {}
    for t, g in df.groupby(df.celltype.astype(str)):
        g = g[np.isfinite(g.pred.to_numpy(float)) & np.isfinite(g.age.to_numpy(float))].copy()
        if path_score is not None and "pathology_score" not in g.columns:
            g = g.merge(path_score, on="donor", how="left")
        agg = dict(pred=("pred", "first"), age=("age", "first"))
        if "pathology_score" in g.columns:
            agg["pathology_score"] = ("pathology_score", "first")
        d = g.groupby("donor", as_index=False).agg(**agg)
        out[t] = d
    return out


def per_type_metrics(type_dfs, residualize_on=None):
    rows = []
    for t, d in type_dfs.items():
        y = d.age.to_numpy(float)
        p = d.pred.to_numpy(float)
        if residualize_on is not None:
            if residualize_on not in d.columns:
                raise StopStep("secondary", f"missing {residualize_on} on type {t}")
            z = d[residualize_on].to_numpy(float)
            y = residualize_1d(y, z)
            p = residualize_1d(p, z)
        n = int((np.isfinite(y) & np.isfinite(p)).sum())
        if n < MIN_N_R:
            rows.append(dict(celltype=t, n=n, r=np.nan, rho=np.nan, r2=np.nan, cal_r2=np.nan, ok=False))
            continue
        sc = pred_scores(y, p)
        rows.append(dict(celltype=t, n=n, r=sc["r"], rho=sc["rho"], r2=sc["r2"], cal_r2=sc["cal_r2"], ok=True))
    tab = pd.DataFrame(rows)
    med = dict(
        r=median_over_types(tab.loc[tab.ok, "r"] if "ok" in tab else tab.r),
        rho=median_over_types(tab.loc[tab.ok, "rho"] if "ok" in tab else tab.rho),
        r2=median_over_types(tab.loc[tab.ok, "r2"] if "ok" in tab else tab.r2),
        cal_r2=median_over_types(tab.loc[tab.ok, "cal_r2"] if "ok" in tab else tab.cal_r2),
        n_types_ok=int(tab.ok.sum()) if len(tab) and "ok" in tab else 0,
        n_types=int(len(tab)),
    )
    return tab, med


def permute_median_r(type_dfs, rng, n_perm, residualize_on=None):
    """Donor-level age permutation; X/pred held fixed. Same shuffle applied across types."""
    donors = sorted(set().union(*[set(d.donor.astype(str)) for d in type_dfs.values()]))
    age_map0 = {}
    for d in type_dfs.values():
        for _, row in d.iterrows():
            age_map0[str(row.donor)] = float(row.age)
    ages0 = np.array([age_map0[d] for d in donors], float)
    null = []
    for i in range(n_perm):
        shuf = rng.permutation(ages0)
        amap = dict(zip(donors, shuf))
        tdfs = {}
        for t, d in type_dfs.items():
            dd = d.copy()
            dd["age"] = dd.donor.astype(str).map(amap).astype(float)
            tdfs[t] = dd
        _, med = per_type_metrics(tdfs, residualize_on=residualize_on)
        null.append(med["r"])
    return np.asarray(null, float)


def bootstrap_median_r(type_dfs, rng, n_boot, residualize_on=None):
    donors = np.array(sorted(set().union(*[set(d.donor.astype(str)) for d in type_dfs.values()])), dtype=object)
    by = {t: d.set_index(d.donor.astype(str)) for t, d in type_dfs.items()}
    vals = []
    for _ in range(n_boot):
        draw = rng.choice(donors, size=len(donors), replace=True)
        tdfs = {}
        for t, d in by.items():
            rows = []
            for k, don in enumerate(draw):
                if don in d.index:
                    rec = d.loc[don]
                    if isinstance(rec, pd.DataFrame):
                        rec = rec.iloc[0]
                    rec = rec.copy()
                    rec["donor"] = f"{don}#{k}"
                    rows.append(rec)
            tdfs[t] = pd.DataFrame(rows).reset_index(drop=True) if rows else pd.DataFrame(columns=d.reset_index().columns)
        _, med = per_type_metrics(tdfs, residualize_on=residualize_on)
        vals.append(med["r"])
    v = np.asarray(vals, float)
    v = v[np.isfinite(v)]
    if v.size < 2:
        return dict(p025=np.nan, p975=np.nan, n=int(v.size))
    return dict(p025=float(np.percentile(v, 2.5)), p975=float(np.percentile(v, 97.5)), n=int(v.size))


def run_cell(name, type_dfs, log, residualize_on=None):
    tab, med = per_type_metrics(type_dfs, residualize_on=residualize_on)
    log(f"[{name}] median r={med['r']:+.3f}  rho={med['rho']:+.3f}  R²={med['r2']:+.3f}  "
        f"calR²={med['cal_r2']:+.3f}  types_ok={med['n_types_ok']}/{med['n_types']}")
    rng_p = np.random.default_rng(EXT_SEED)
    log(f"[{name}] permutation null n_perm={N_PERM} ...")
    null = permute_median_r(type_dfs, rng_p, N_PERM, residualize_on=residualize_on)
    ns = summarize_null_col(null)
    pval = permutation_p(med["r"], null, greater=True)
    log(f"[{name}] r_null mean={ns['mean']:+.3f}  p95={ns['p95']:+.3f}  p={pval:.3f}")
    rng_b = np.random.default_rng(BOOT_SEED)
    log(f"[{name}] donor bootstrap B={N_BOOT} ...")
    ci = bootstrap_median_r(type_dfs, rng_b, N_BOOT, residualize_on=residualize_on)
    log(f"[{name}] r 95% CI [{ci['p025']:+.3f}, {ci['p975']:+.3f}]")
    passed = pass_r_bar(med["r"], ns["mean"], bar=TRANSFER_NULL_BAR)
    rec = dict(
        cell=name,
        r=med["r"], rho=med["rho"], r2=med["r2"], cal_r2=med["cal_r2"],
        r_null=ns["mean"], r_null_p95=ns["p95"], r_p=pval,
        n_perm=N_PERM, n_boot=N_BOOT,
        r_ci_lo=ci["p025"], r_ci_hi=ci["p975"],
        pass_r=bool(passed),
        n_types_ok=med["n_types_ok"], n_types=med["n_types"],
        residualize_on=residualize_on,
        n_donors=int(len(set().union(*[set(d.donor.astype(str)) for d in type_dfs.values()]))),
    )
    tab.to_csv(EXT_DIR / f"{name}_per_type.csv", index=False)
    dump_json(EXT_DIR / f"{name}.json", jsonable(rec))
    return rec, tab


def project():
    if not PREREG_FLAG.exists():
        raise StopStep("prereg", "pre-registration flag missing; refusing to score SEA-AD")
    if not (EXT_DIR / "frozen_directions.npz").exists():
        raise StopStep("project", "frozen_directions.npz missing — run freeze first")
    log = Logger(EXT_DIR / "project_report.txt")
    try:
        return _project(log)
    finally:
        log.close()


def _project(log):
    ext_log_banner(log, "project frozen DLPFC directions onto SEA-AD MTG")
    summ = load_json(EXT_DIR / "stage0_summary.json")
    mapping = pd.read_csv(EXT_DIR / "stage0_mapping.csv")
    donors = pd.read_csv(EXT_DIR / "stage0_donor_meta.csv")
    donors["donor"] = donors.donor.astype(str)
    donors["age"] = pd.to_numeric(donors["age"], errors="coerce")
    donors["pathology_score"] = pd.to_numeric(donors["pathology_score"], errors="coerce")
    pl = donors["pathology_low"]
    if pl.dtype != bool:
        donors["pathology_low"] = pl.map(
            lambda x: str(x).strip().lower() in ("true", "1", "yes")
        )
    frozen = load_frozen()
    types = [str(t) for t in frozen["types"]]
    url = summ["h5ad_url"]
    src = local_h5ad_if_complete(url, summ.get("h5ad_bytes")) or url
    log(f"[project] expression source: {src}")

    cache = EXT_DIR / "stage1_pseudobulk.npz"
    if cache.exists():
        log(f"[project] loading cached pseudobulk {cache}")
        z = np.load(cache, allow_pickle=True)
        pb_obs = pd.DataFrame({k: z[k] for k in ("donor", "celltype", "n_cells")})
        pb_obs["donor"] = pb_obs.donor.astype(str)
        pb_obs["celltype"] = pb_obs.celltype.astype(str)
        X_sum = z["X"]
        gene_id = z["gene_id"]
    else:
        pb_obs, X_sum, gene_id = build_seaad_pseudobulk(src, mapping, log)
        np.savez_compressed(
            cache,
            X=X_sum.astype(np.float32),
            gene_id=np.asarray(gene_id).astype(str),
            donor=pb_obs.donor.astype(str).to_numpy(),
            celltype=pb_obs.celltype.astype(str).to_numpy(),
            n_cells=pb_obs.n_cells.to_numpy(),
        )
        log(f"[project] wrote {cache}")

    Y = _norm_logcpm(X_sum)
    Y_al, gene_rec = align_genes(Y, gene_id, frozen, log)
    dump_json(EXT_DIR / "gene_overlap.json", jsonable(gene_rec))
    mu = np.asarray(frozen["mu"], float)
    sd = np.asarray(frozen["sd"], float)
    sd = np.where(sd < 1e-12, 1.0, sd)
    Xz = (Y_al - mu) / sd
    basis = np.asarray(frozen["identity_basis"], float)
    R = residualize_identity(Xz, basis)

    pb_obs = pb_obs.merge(donors[["donor", "age", "sex", "pathology_score", "pathology_low"]], on="donor", how="left")
    if pb_obs.age.isna().any():
        n = int(pb_obs.age.isna().sum())
        raise StopStep("project", f"{n} pseudobulks missing age after donor merge")

    path_don = donors[["donor", "pathology_score"]].copy()
    results = {}
    for kind, Xuse, a_key, b_key, Wkey in (
        ("raw", Xz, "cal_raw_intercept", "cal_raw_slope", "W_raw"),
        ("idresid", R, "cal_id_intercept", "cal_id_slope", "W_idresid"),
    ):
        log(f"\n=== direction {kind} ===")
        pred, score, meta = scores_for_direction(
            Xuse, frozen[Wkey], frozen[a_key], frozen[b_key], types, pb_obs, kind)
        meta.to_csv(EXT_DIR / f"project_{kind}_apply.csv", index=False)
        pb_obs[f"pred_{kind}"] = pred
        pb_obs[f"score_{kind}"] = score

        prim_don = donors.loc[donors.pathology_low.astype(bool) & donors.age.notna(), "donor"].astype(str)
        log(f"[primary] no/low pathology donors={prim_don.nunique()}  "
            f"(field={summ.get('pathology_field')})")
        if int(prim_don.nunique()) < MIN_N_R:
            raise StopStep(
                "primary",
                f"primary cell has {int(prim_don.nunique())} no/low-pathology donors; need ≥{MIN_N_R}. Not widening the set.",
                dict(n=int(prim_don.nunique())),
            )
        t_prim = _type_vectors(pb_obs.assign(pred=pred), pred, donors_keep=set(prim_don))
        rec_p, _ = run_cell(f"primary_{kind}", t_prim, log, residualize_on=None)
        rec_p["n_donors_low"] = int(prim_don.nunique())
        rec_p["pathology_field"] = summ.get("pathology_field")
        dump_json(EXT_DIR / f"primary_{kind}.json", jsonable(rec_p))

        sec_don = donors.loc[donors.pathology_score.notna() & donors.age.notna(), "donor"].astype(str)
        log(f"[secondary] donors with pathology score={sec_don.nunique()}")
        t_sec = _type_vectors(
            pb_obs.assign(pred=pred), pred, donors_keep=set(sec_don), path_score=path_don)
        rec_s, _ = run_cell(f"secondary_{kind}", t_sec, log, residualize_on="pathology_score")
        rec_s["n_donors_all"] = int(sec_don.nunique())
        rec_s["pathology_field"] = summ.get("pathology_field")
        dump_json(EXT_DIR / f"secondary_{kind}.json", jsonable(rec_s))
        results[kind] = dict(primary=rec_p, secondary=rec_s)

    pb_obs.to_csv(EXT_DIR / "project_pseudobulk_obs.csv", index=False)
    dump_json(EXT_DIR / "project_summary.json", jsonable(dict(
        gene_overlap=gene_rec, results=results, n_pseudobulks=int(len(pb_obs)),
        n_donors_pb=int(pb_obs.donor.nunique()),
        pathology_field=summ.get("pathology_field"),
        age_min=summ.get("age_min"), age_max=summ.get("age_max"),
        tissue=summ.get("tissue"),
    )))
    refresh_progress(next_action="write FINDINGS_EXTERNAL.md")
    return results


if __name__ == "__main__":
    project()

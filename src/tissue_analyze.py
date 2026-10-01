"""T2 batch audit + T3 turnover test for one downloaded tissue atlas.

Reuses OneK1K code:
  - src/pseudobulk.py (chunked donor x cell-type sums, filter log)
  - src/decomp.normalize_logcpm
  - src/rerun_stageB.py (kernel ridge, batch-grouped CV, within-batch R2, CLR, bootstrap)

HARD RULES: batch-grouped CV (not donor-grouped); every filter logged with counts; seed recorded.
Usage: python src/tissue_analyze.py muscle [brain_aging retina_sn ...]
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tissue_common import (  # noqa: E402
    TISSUE_DIR, TISSUE_FIG, TISSUE_RAW, TISSUE_PROC, TISSUE_SEED, Logger,
    BLOOD_COMP_R2, BLOOD_EXPR_R2, BATCH_COL_CANDIDATES,
)
from tissue_download_top import resolve, local_path  # noqa: E402
from tissue_peek_obs import peek  # noqa: E402
from pseudobulk import build, parse_age, FILTER_LOG  # noqa: E402
from decomp import normalize_logcpm  # noqa: E402

# Reused Stage B helpers (copied so we do not import rerun_stageB.py, which pulls a
# statsmodels/scipy combo that is broken on this machine). Logic is unchanged.
N_REPEATS = 10
N_OUTER = 5
N_INNER = 3
N_BOOT = 2000
MIN_CELLS = 20
ADULT_MIN_AGE = 18.0
ALPHAS = np.array([0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000])


def pool_center(M, pool):
    M = np.asarray(M, float)
    out = M.copy()
    for p in np.unique(pool):
        m = pool == p
        out[m] -= M[m].mean(0, keepdims=True)
    return out


def within_pool_r2_mae(y_c, pred, pool):
    pred_c = pool_center(pred[:, None], pool)[:, 0]
    ss_tot = (y_c ** 2).sum()
    if ss_tot <= 0:
        return np.nan, np.abs(y_c - pred_c).mean()
    return 1 - ((y_c - pred_c) ** 2).sum() / ss_tot, np.abs(y_c - pred_c).mean()


def raw_r2_mae(y, pred):
    ss_tot = ((y - y.mean()) ** 2).sum()
    if ss_tot <= 0:
        return np.nan, np.abs(y - pred).mean()
    return 1 - ((y - pred) ** 2).sum() / ss_tot, np.abs(y - pred).mean()


def make_pool_folds(pools_all, n_folds, rng):
    pools = np.array(sorted(pools_all))
    rng.shuffle(pools)
    return {p: i % n_folds for i, p in enumerate(pools)}


def kernel_ridge_oof(K, y, pool, donor, fold_of_pool, alphas, rng, n_inner=N_INNER, obs_for_assert=None):
    n = len(y)
    fold = np.array([fold_of_pool[p] for p in pool])
    pred = np.full(n, np.nan)
    chosen = []
    for k in np.unique(fold):
        te = np.flatnonzero(fold == k)
        tr = np.flatnonzero(fold != k)
        assert not (set(donor[tr]) & set(donor[te])), "donor leakage"
        assert not (set(pool[tr]) & set(pool[te])), "pool leakage"
        tr_pools = np.array(sorted(set(pool[tr])))
        rng.shuffle(tr_pools)
        n_in = max(2, min(n_inner, len(tr_pools)))
        inner_of_pool = {p: i % n_in for i, p in enumerate(tr_pools)}
        inner_fold = np.array([inner_of_pool[p] for p in pool[tr]])
        score = np.zeros(len(alphas))
        for ki in range(n_in):
            ite = tr[inner_fold == ki]
            itr = tr[inner_fold != ki]
            if len(itr) < 3 or len(ite) < 1:
                continue
            mu_i = y[itr].mean()
            S, U = np.linalg.eigh(K[np.ix_(itr, itr)])
            Uy = U.T @ (y[itr] - mu_i)
            Kte = K[np.ix_(ite, itr)]
            for a_i, a in enumerate(alphas):
                dual = U @ (Uy / (S + a))
                p_ = mu_i + Kte @ dual
                score[a_i] += -((y[ite] - p_) ** 2).sum()
        a_best = alphas[int(np.argmax(score))]
        chosen.append(float(a_best))
        mu = y[tr].mean()
        S, U = np.linalg.eigh(K[np.ix_(tr, tr)])
        dual = U @ ((U.T @ (y[tr] - mu)) / (S + a_best))
        pred[te] = mu + K[np.ix_(te, tr)] @ dual
    assert not np.isnan(pred).any()
    return pred, chosen


def bootstrap_ci(y, pred, pool, rng, n_boot=N_BOOT, within=True):
    pools = np.unique(pool)
    idx_by_pool = [np.flatnonzero(pool == p) for p in pools]
    pred_use = pool_center(pred[:, None], pool)[:, 0] if within else pred
    vals = np.empty(n_boot)
    for b in range(n_boot):
        sel = rng.integers(0, len(pools), len(pools))
        idx = np.concatenate([idx_by_pool[s] for s in sel])
        yy, pp = y[idx], pred_use[idx]
        if within:
            ss = (yy ** 2).sum()
            vals[b] = 1 - ((yy - pp) ** 2).sum() / ss if ss > 0 else np.nan
        else:
            ss = ((yy - yy.mean()) ** 2).sum()
            vals[b] = 1 - ((yy - pp) ** 2).sum() / ss if ss > 0 else np.nan
    return float(np.nanpercentile(vals, 2.5)), float(np.nanpercentile(vals, 97.5))


def _safe_ols_r2(y, group):
    """R^2 of y ~ C(group) via QR. No statsmodels (broken vs current scipy on this machine)."""
    y = np.asarray(y, float)
    g = pd.Series(group).astype(str).to_numpy()
    ok = np.isfinite(y)
    y, g = y[ok], g[ok]
    if len(y) < 4 or pd.Series(g).nunique() < 2:
        return np.nan
    X = pd.get_dummies(pd.Categorical(g), dtype=float).to_numpy()
    yc = y - y.mean()
    Q, _ = np.linalg.qr(X)
    ss_res = ((yc - Q @ (Q.T @ yc)) ** 2).sum()
    ss_tot = (yc ** 2).sum()
    if ss_tot <= 0:
        return np.nan
    return float(1 - ss_res / ss_tot)


def pick_batch(donor_df, log):
    """Score donor-level categorical columns. A usable batch groups >=2 donors in >=2 levels
    and is not a 1-1 alias of donor."""
    n_don = len(donor_df)
    rows = []
    skip = {"donor", "age", "sex", "n_cells", "mean_counts_per_cell", "mean_genes_per_cell"}
    for col in list(donor_df.columns):
        if col in skip:
            continue
        s = donor_df[col].astype(str)
        n_lev = s.nunique(dropna=True)
        sizes = s.value_counts()
        n_shared = int((sizes >= 2).sum())
        frac = float((s.map(sizes) >= 2).mean())
        aliased = n_lev == n_don
        usable = (n_lev >= 2) and (n_shared >= 2) and (not aliased)
        # prefer columns that look like technical batches
        tech = col.lower() in {c.lower() for c in BATCH_COL_CANDIDATES} or any(
            k in col.lower() for k in ("batch", "pool", "library", "dataset", "assay", "suspension", "sample")
        )
        rows.append(dict(column=col, n_levels=n_lev, n_levels_ge2_donors=n_shared,
                         frac_donors_in_shared=frac, aliased_with_donor=aliased,
                         looks_technical=tech, usable=usable))
    sc = pd.DataFrame(rows).sort_values(
        ["usable", "looks_technical", "frac_donors_in_shared", "n_levels_ge2_donors"],
        ascending=[False, False, False, False],
    )
    log("\n[T2] donor-level batch-column audit:")
    if len(sc):
        log(sc.to_string(index=False))
    usable = sc[sc.usable]
    if usable.empty:
        log("[T2] NO usable batch column (no shared-donor technical factor). "
            "Within-batch age R2 is not defined. Dataset flagged unusable for the OneK1K correction.")
        return None, sc
    pick = usable.iloc[0].column
    log(f"[T2] selected batch column: {pick}")
    return pick, sc


def build_composition(obs, donors, log, min_mean_prop=0.005, min_donors_type=8):
    """Per-donor CLR composition from ALL cells (not the >=20-cell pseudobulk groups)."""
    n0 = len(obs)
    o = obs[obs.donor.isin(donors.index)].copy()
    log(f"[filter] cells of analysis donors                              cells {n0:,} -> {len(o):,} (removed {n0-len(o):,})")
    vc = o.celltype.value_counts()
    n_don_ct = o.groupby("celltype").donor.nunique()
    mean_prop = o.groupby("donor").celltype.value_counts(normalize=True).unstack(fill_value=0).mean()
    keep_ct = [c for c in vc.index
               if mean_prop.get(c, 0) >= min_mean_prop and n_don_ct.get(c, 0) >= min_donors_type]
    o["cat"] = o.celltype.where(o.celltype.isin(keep_ct), "other")
    dropped = sorted(set(vc.index) - set(keep_ct))
    n_other = int((o.cat == "other").sum())
    log(f"[filter] cell types kept as composition categories "
        f"(mean prop>={min_mean_prop}, donors>={min_donors_type}): {len(keep_ct)} kept, "
        f"{len(dropped)} collapsed to 'other' ({n_other:,} cells, {100*n_other/len(o):.2f}%) "
        f"dropped={dropped[:20]}")
    cats = keep_ct + (["other"] if n_other else [])
    counts = o.groupby(["donor", "cat"], observed=True).size().unstack(fill_value=0).reindex(donors.index).fillna(0)
    counts = counts.reindex(columns=cats, fill_value=0)
    Ncell = counts.sum(1)
    prop = counts.div(Ncell.replace(0, np.nan), axis=0).fillna(0)
    K = len(cats)
    clr = np.log((counts + 0.5).div(Ncell + 0.5 * K, axis=0))
    clr = clr.sub(clr.mean(1), axis=0)
    log(f"cells/donor for composition: median {Ncell.median():.0f} (min {Ncell.min():.0f}, max {Ncell.max():.0f}); "
        f"CLR dim={clr.shape[1]}")
    return counts, prop, clr, cats, Ncell


def run_one(name, log):
    rng = np.random.default_rng([TISSUE_SEED, abs(hash(name)) % (2**31)])
    info = resolve(name)
    path = local_path(name)
    if not path.exists():
        raise FileNotFoundError(f"{path} missing — run tissue_download_top.py {name}")
    log("=" * 100)
    log(f"TISSUE ANALYZE  name={name}  seed={TISSUE_SEED}  dataset_id={info['dataset_id']}")
    log(f"  title={info['title']}  tissue={info['tissue']}")
    log(f"  local={path}  size={path.stat().st_size:,}")
    log("=" * 100)

    peek_sum = peek(str(path), name, log)

    # ---- obs-level filters (backed, no X) ----
    import anndata as ad
    a = ad.read_h5ad(path, backed="r")
    obs = a.obs.copy()
    obs["_pos"] = np.arange(len(obs))
    n0 = len(obs)
    log(f"[load] {n0:,} cells; obs cols={list(obs.columns)}")

    require_primary = True
    if "is_primary_data" in obs:
        n_pri = int(obs.is_primary_data.astype(bool).sum())
        frac_pri = n_pri / max(n0, 1)
        # CELLxGENE marks cells non-primary when they also appear in a subset file of the
        # same collection. Dropping them can empty an "all cells" atlas (retina scRNA:
        # 266k -> 6k). Only enforce the flag when it keeps a majority of cells.
        if frac_pri < 0.5:
            require_primary = False
            log(f"[filter] is_primary_data==True would keep {n_pri:,}/{n0:,} ({100*frac_pri:.1f}%) "
                f"— NOT applied (cells are duplicated into subset files of the same collection)")
        else:
            keep = obs.is_primary_data.astype(bool)
            log(f"[filter] is_primary_data==True                               cells {n0:,} -> {int(keep.sum()):,} "
                f"(removed {n0-int(keep.sum()):,})")
            obs = obs[keep]
            n0 = len(obs)

    age_col = "development_stage" if "development_stage" in obs.columns else ("age" if "age" in obs.columns else None)
    if age_col is None:
        raise RuntimeError("no development_stage or age column")
    if age_col == "development_stage":
        obs["age"] = obs[age_col].map(parse_age).astype(float)
    else:
        obs["age"] = pd.to_numeric(obs[age_col], errors="coerce")
    if "donor_age" in obs.columns:
        alt = pd.to_numeric(obs["donor_age"], errors="coerce")
        filled = obs["age"].isna() & alt.notna()
        obs.loc[filled, "age"] = alt[filled]
        log(f"[age] filled {int(filled.sum()):,} cells from donor_age where development_stage was non-numeric")
    n = len(obs)
    obs = obs[obs.age.notna()]
    log(f"[filter] cells with numeric donor age                          cells {n:,} -> {len(obs):,} (removed {n-len(obs):,})")
    n = len(obs)
    obs = obs[obs.age >= ADULT_MIN_AGE]
    log(f"[filter] adult cells (age >= {ADULT_MIN_AGE:.0f})                         cells {n:,} -> {len(obs):,} "
        f"(removed {n-len(obs):,})")

    obs["donor"] = obs["donor_id"].astype(str) if "donor_id" in obs else obs.index.astype(str)
    ct_col = "cell_type" if "cell_type" in obs.columns else None
    if ct_col is None:
        raise RuntimeError("no cell_type column")
    obs["celltype"] = obs[ct_col].astype(str)
    if "sex" in obs.columns:
        obs["sex"] = obs["sex"].astype(str)

    # donor table + candidate batch columns from obs
    extra = [c for c in (
        "assay", "suspension_type", "tissue", "disease", "batch", "library_id", "library_uuid",
        "sample_id", "sample_uuid", "dataset", "dataset_id", "self_reported_ethnicity",
        "10XBatch", "SequencingPool", "sequencing_platform", "sample_preservation_method",
        "tissue_handling_interval", "AgeGroup", "Dataset", "dataset",
        "study_name", "institute", "sample_collection_year",
    ) if c in obs.columns]
    agg = dict(age=("age", "first"), n_cells=("age", "size"), age_nunique=("age", "nunique"))
    if "sex" in obs:
        agg["sex"] = ("sex", "first")
    for c in extra:
        def _mode(s):
            s = s.dropna().astype(str)
            if s.empty:
                return "NA"
            return s.value_counts().index[0]
        agg[c] = (c, _mode)
    donors = obs.groupby("donor", observed=True).agg(**agg)
    n_bad_age = int((donors.age_nunique != 1).sum())
    if n_bad_age:
        log(f"[warn] {n_bad_age} donors have >1 age value; using first")
    log(f"[donors] {len(donors)} adult donors; age {donors.age.min():.0f}–{donors.age.max():.0f} "
        f"(nunique={donors.age.nunique()}, span={donors.age.max()-donors.age.min():.0f} y); "
        f"cells/donor median {donors.n_cells.median():.0f}")
    if len(donors) < 10 or (donors.age.max() - donors.age.min()) < 40:
        log("[STOP] after adult filter: <10 donors or span<40 y. Continuing the audit but T3 will be flagged.")

    batch_col, batch_scores = pick_batch(donors, log)
    batch_usable = batch_col is not None
    if batch_usable:
        donors["batch"] = donors[batch_col].astype(str)
        r2_age_batch = _safe_ols_r2(donors.age, donors.batch)
        log(f"[T2] R2(age ~ C(batch={batch_col})) = {r2_age_batch:.3f}   "
            f"n_batches={donors.batch.nunique()}  donors/batch: "
            f"median {donors.groupby('batch').size().median():.1f} "
            f"(min {donors.groupby('batch').size().min()}, max {donors.groupby('batch').size().max()})")
    else:
        donors["batch"] = donors.index.astype(str)  # placeholder; T3 within-batch undefined
        r2_age_batch = 1.0
        log("[T2] R2(age ~ batch) := 1.0 by construction (batch aliases donor)")

    # ---- pseudobulk (needs X; writes processed file) ----
    FILTER_LOG.clear()
    extra_obs = extra + ([] if batch_col is None or batch_col in extra else [batch_col])
    min_don_ct = max(8, min(10, len(donors) // 3))
    age_col_pb = "development_stage" if "development_stage" in obs.columns else "age"
    sex_col_pb = "sex" if "sex" in obs.columns else "sex"
    try:
        a.file.close()
    except Exception:
        pass
    del a
    log(f"\n[pseudobulk] min_cells={MIN_CELLS}  min_donors_per_type={min_don_ct}")
    pb = build(
        str(path), f"tissue_{name}",
        donor_col="donor_id", ct_col="cell_type", age_col=age_col_pb,
        sex_col=sex_col_pb,
        extra_obs=tuple(extra_obs),
        min_cells=MIN_CELLS, min_donors=min_don_ct,
        require_primary=require_primary,
    )
    # restrict to adult donors that survived composition filters
    n_pb = len(pb)
    keep_d = set(donors.index)
    pb = pb[pb.obs.donor.isin(keep_d)].copy()
    log(f"[filter] pseudobulks of adult analysis donors                  groups {n_pb:,} -> {pb.n_obs:,} "
        f"(removed {n_pb-pb.n_obs:,})")
    # drop cell types that now have too few adult donors
    ct_n = pb.obs.groupby("celltype").donor.nunique()
    keep_ct = ct_n[ct_n >= min_don_ct].index
    n1 = pb.n_obs
    pb = pb[pb.obs.celltype.isin(keep_ct)].copy()
    log(f"[filter] cell types with >= {min_don_ct} adult donors            groups {n1:,} -> {pb.n_obs:,}  "
        f"types {len(ct_n)} -> {len(keep_ct)}  dropped={dict(ct_n[ct_n < min_don_ct])}")

    Y, gid, sym = normalize_logcpm(pb, log=log)
    pb.obs["donor"] = pb.obs["donor"].astype(str)
    pb.obs["celltype"] = pb.obs["celltype"].astype(str)
    if batch_usable:
        # map donor -> batch
        pb.obs["batch"] = pb.obs.donor.map(donors.batch)
    else:
        pb.obs["batch"] = pb.obs.donor

    # donor-level depth from the (adult) obs we already have? compute from pb
    depth_d = pb.obs.groupby("donor").apply(
        lambda g: np.average(g.mean_counts_per_cell, weights=g.n_cells), include_groups=False
    )
    genes_d = pb.obs.groupby("donor").apply(
        lambda g: np.average(g.mean_genes_per_cell, weights=g.n_cells), include_groups=False
    )
    if pb.n_obs == 0 or len(depth_d) == 0:
        log("[STOP] no surviving pseudobulks after filters — T3 not run")
        t2 = dict(name=name, dataset_id=info["dataset_id"], tissue=info["tissue"],
                  n_donors=int(len(donors)), batch_usable=bool(batch_usable),
                  batch_column=batch_col, error="no_pseudobulks")
        pd.DataFrame([t2]).to_csv(TISSUE_DIR / f"t2_{name}_batch_audit.csv", index=False)
        return dict(t2=t2, t3=dict(name=name, error="no_pseudobulks"))
    donors = donors.join(depth_d.rename("mean_counts_per_cell"), how="inner")
    donors = donors.join(genes_d.rename("mean_genes_per_cell"), how="left")
    log(f"[filter] donors with at least one retained pseudobulk          donors {len(keep_d):,} -> {len(donors):,}")
    # composition must use the same donor set as T2/T3
    counts, prop, clr, cats, Ncell = build_composition(obs, donors, log, min_donors_type=min_don_ct)
    donors = donors.join(Ncell.rename("n_cells_comp"), how="left")

    # ---- T2 depth ----
    if batch_usable:
        r2_depth_batch = _safe_ols_r2(np.log(donors.mean_counts_per_cell), donors.batch)
        log(f"[T2] R2(log mean UMI/cell ~ C(batch)) = {r2_depth_batch:.3f}")
        donors["age_b"] = donors.groupby("batch").age.transform("mean")
        donors["dep_b"] = donors.groupby("batch").mean_counts_per_cell.transform("mean")
        bm = donors.groupby("batch").agg(age=("age", "mean"), depth=("mean_counts_per_cell", "mean"))
        r_between = float(np.corrcoef(bm.age, bm.depth)[0, 1]) if len(bm) > 2 else np.nan
        r_within = float(np.corrcoef(donors.age - donors.age_b, donors.mean_counts_per_cell - donors.dep_b)[0, 1])
        r_pooled = float(np.corrcoef(donors.age, donors.mean_counts_per_cell)[0, 1])
        log(f"[T2] corr(age, depth): pooled r={r_pooled:+.3f}; between-batch r={r_between:+.3f} "
            f"(n={len(bm)} batches); within-batch r={r_within:+.3f}")
    else:
        r2_depth_batch = np.nan
        r_between = np.nan
        r_within = np.nan
        r_pooled = float(np.corrcoef(donors.age, donors.mean_counts_per_cell)[0, 1])
        log(f"[T2] corr(age, depth) pooled r={r_pooled:+.3f}  (no shared batch; within/between undefined)")

    t2 = dict(name=name, dataset_id=info["dataset_id"], tissue=info["tissue"],
              n_donors=int(len(donors)), n_batches=int(donors.batch.nunique()) if batch_usable else None,
              batch_column=batch_col, batch_usable=bool(batch_usable),
              r2_age_batch=None if not batch_usable else float(r2_age_batch),
              r2_depth_batch=None if not np.isfinite(r2_depth_batch) else float(r2_depth_batch),
              r_age_depth_pooled=float(r_pooled),
              r_age_depth_between=None if not np.isfinite(r_between) else float(r_between),
              r_age_depth_within=None if not np.isfinite(r_within) else float(r_within),
              age_min=float(donors.age.min()), age_max=float(donors.age.max()),
              n_cell_types=int(pb.obs.celltype.nunique()), n_pseudobulks=int(pb.n_obs))
    pd.DataFrame([t2]).to_csv(TISSUE_DIR / f"t2_{name}_batch_audit.csv", index=False)
    batch_scores.to_csv(TISSUE_DIR / f"t2_{name}_batch_columns.csv", index=False)

    # ---- T3 ----
    log("\n" + "=" * 100)
    log("[T3] batch-grouped CV: composition vs within-cell-type expression")
    log("=" * 100)
    if not batch_usable:
        log("[T3] batch is not recoverable as a shared-donor factor. Within-batch R2 is UNDEFINED. "
            "A donor-grouped naive R2 is reported as a secondary number and is NOT comparable to blood 0.493 / 0.482.")
    pool = donors.batch.to_numpy()
    n_batch = donors.batch.nunique()
    n_outer = int(min(N_OUTER, n_batch)) if batch_usable else 5
    n_rep = N_REPEATS if batch_usable and n_batch >= 5 else (max(3, n_batch) if batch_usable else N_REPEATS)
    log(f"  n_donors={len(donors)} n_batches={n_batch} n_outer={n_outer} n_repeats={n_rep}  "
        f"within-batch metric {'ON' if batch_usable else 'OFF (naive only)'}")

    # kernels
    types = sorted(pb.obs.celltype.unique())
    donor_idx = pd.Index(donors.index)
    obs_pb = pb.obs
    ct = obs_pb.celltype.to_numpy()
    pool_pb = obs_pb.batch.to_numpy()
    donor_pb = obs_pb.donor.to_numpy()
    kern = {}
    log("building per-type kernels:")
    for t in types:
        rows_t = np.flatnonzero(ct == t)
        X = Y[rows_t]
        if batch_usable:
            Xc = pool_center(X, pool_pb[rows_t])
        else:
            Xc = X - X.mean(0)
        sd = Xc.std(0); sd[sd == 0] = 1; Xc = Xc / sd
        Xn = X - X.mean(0); sdn = Xn.std(0); sdn[sdn == 0] = 1; Xn = Xn / sdn
        kern[t] = (rows_t, Xc @ Xc.T / X.shape[1], Xn @ Xn.T / X.shape[1])
        log(f"   {t:<40} n={len(rows_t):>4} batches={len(np.unique(pool_pb[rows_t]))}")
    nD = len(donors)
    K_all_w = np.zeros((nD, nD)); K_all_n = np.zeros((nD, nD))
    for t in types:
        rows_t, Kw, Kn = kern[t]
        di = donor_idx.get_indexer(donor_pb[rows_t])
        K_all_w[np.ix_(di, di)] += Kw
        K_all_n[np.ix_(di, di)] += Kn
    K_all_w /= max(len(types), 1); K_all_n /= max(len(types), 1)
    pool_d = donors.batch.to_numpy()
    if batch_usable:
        Cw = pool_center(clr.to_numpy(), pool_d)
    else:
        Cw = clr.to_numpy() - clr.to_numpy().mean(0)
    sd = Cw.std(0); sd[sd == 0] = 1; Cw = Cw / sd
    Cn = clr.to_numpy() - clr.to_numpy().mean(0); sdn = Cn.std(0); sdn[sdn == 0] = 1; Cn = Cn / sdn
    K_comp_w = Cw @ Cw.T / max(Cw.shape[1], 1)
    K_comp_n = Cn @ Cn.T / max(Cn.shape[1], 1)
    y_d = donors.age.to_numpy(float)
    if batch_usable:
        donors["age_c"] = donors.age - donors.groupby("batch").age.transform("mean")
        yc_d = donors.age_c.to_numpy(float)
    else:
        donors["age_c"] = donors.age - donors.age.mean()
        yc_d = donors.age_c.to_numpy(float)

    results = []
    oof_store = {}
    for rep in range(n_rep):
        rrng = np.random.default_rng([TISSUE_SEED, abs(hash(name)) % (2**31), rep])
        fold_of_pool = make_pool_folds(donors.batch.unique(), n_outer, rrng)
        for label, Kw_, Kn_ in (
            (f"composition (CLR, {clr.shape[1]} categories)", K_comp_w, K_comp_n),
            (f"expression: all {len(types)} cell types jointly", K_all_w, K_all_n),
        ):
            pw, aw = kernel_ridge_oof(Kw_, yc_d, pool_d, donors.index.to_numpy(), fold_of_pool, ALPHAS, rrng, n_inner=min(N_INNER, n_outer))
            pn, an = kernel_ridge_oof(Kn_, y_d, pool_d, donors.index.to_numpy(), fold_of_pool, ALPHAS, rrng, n_inner=min(N_INNER, n_outer))
            r2w, maew = within_pool_r2_mae(yc_d, pw, pool_d)
            r2n, maen = raw_r2_mae(y_d, pn)
            results.append(dict(repeat=rep, model=label, n=nD, r2_within_batch=r2w, mae_within=maew,
                                r2_naive=r2n, mae_naive=maen))
            if rep == 0:
                oof_store[label] = (pw, pn)
        log(f"   repeat {rep}: comp R2w={results[-2]['r2_within_batch']:+.3f}  "
            f"expr R2w={results[-1]['r2_within_batch']:+.3f}  "
            f"comp naive={results[-2]['r2_naive']:+.3f}  expr naive={results[-1]['r2_naive']:+.3f}")
    res = pd.DataFrame(results)
    res.to_csv(TISSUE_DIR / f"t3_{name}_cv_repeats.csv", index=False)

    def _summ(sub, pw, pn, y, yc, pool):
        lo, hi = (np.nan, np.nan)
        lon, hin = (np.nan, np.nan)
        if batch_usable:
            lo, hi = bootstrap_ci(yc, pw, pool, np.random.default_rng([TISSUE_SEED, 99]), within=True, n_boot=N_BOOT)
        lon, hin = bootstrap_ci(y, pn, pool, np.random.default_rng([TISSUE_SEED, 98]), within=False, n_boot=N_BOOT)
        return dict(r2_within_mean=float(sub.r2_within_batch.mean()),
                    r2_within_min=float(sub.r2_within_batch.min()),
                    r2_within_max=float(sub.r2_within_batch.max()),
                    r2_within_lo=None if not np.isfinite(lo) else float(lo),
                    r2_within_hi=None if not np.isfinite(hi) else float(hi),
                    r2_naive_mean=float(sub.r2_naive.mean()),
                    r2_naive_lo=float(lon), r2_naive_hi=float(hin),
                    mae_within=float(sub.mae_within.mean()), mae_naive=float(sub.mae_naive.mean()))

    labels = list(dict.fromkeys(res.model))
    summ_rows = []
    for lab in labels:
        sub = res[res.model == lab]
        pw, pn = oof_store[lab]
        s = _summ(sub, pw, pn, y_d, yc_d, pool_d)
        s.update(model=lab, name=name)
        summ_rows.append(s)
        log(f"  {lab}: within-batch R2={s['r2_within_mean']:+.3f} "
            f"[{s['r2_within_lo']}, {s['r2_within_hi']}] {{{s['r2_within_min']:+.3f},{s['r2_within_max']:+.3f}}}  "
            f"naive R2={s['r2_naive_mean']:+.3f} [{s['r2_naive_lo']:+.3f},{s['r2_naive_hi']:+.3f}]")
    summ = pd.DataFrame(summ_rows)
    summ.to_csv(TISSUE_DIR / f"t3_{name}_cv_summary.csv", index=False)

    comp = summ.iloc[0]; expr = summ.iloc[1]
    # ratio uses the same metric as blood when batch is usable; else naive (flagged)
    if batch_usable:
        ratio = expr.r2_within_mean / comp.r2_within_mean if abs(comp.r2_within_mean) > 1e-6 else np.nan
        metric = "within_batch"
        c_r2, e_r2 = comp.r2_within_mean, expr.r2_within_mean
    else:
        ratio = expr.r2_naive_mean / comp.r2_naive_mean if abs(comp.r2_naive_mean) > 1e-6 else np.nan
        metric = "naive_donor_grouped_NOT_comparable_to_blood"
        c_r2, e_r2 = comp.r2_naive_mean, expr.r2_naive_mean
    log(f"\n[T3] {name}: composition R2={c_r2:+.3f}  expression R2={e_r2:+.3f}  "
        f"expr/comp={ratio:+.3f}  metric={metric}")
    log(f"     blood reference: composition {BLOOD_COMP_R2:.3f}  expression {BLOOD_EXPR_R2:.3f}  "
        f"expr/comp={BLOOD_EXPR_R2/BLOOD_COMP_R2:.3f}")

    t3 = dict(name=name, dataset_id=info["dataset_id"], tissue=info["tissue"],
              n_donors=int(nD), n_cell_types=len(types), n_comp_categories=int(clr.shape[1]),
              batch_usable=bool(batch_usable), batch_column=batch_col, metric=metric,
              composition_r2=float(c_r2), expression_r2=float(e_r2),
              ratio_expr_over_comp=None if not np.isfinite(ratio) else float(ratio),
              composition_r2_within=float(comp.r2_within_mean),
              expression_r2_within=float(expr.r2_within_mean),
              composition_r2_naive=float(comp.r2_naive_mean),
              expression_r2_naive=float(expr.r2_naive_mean),
              expression_beats_composition=bool(e_r2 > c_r2),
              blood_composition_r2=BLOOD_COMP_R2, blood_expression_r2=BLOOD_EXPR_R2)
    json.dump(dict(t2=t2, t3=t3, peek={k: v for k, v in peek_sum.items() if not isinstance(v, (list, dict))},
                   seed=TISSUE_SEED, types=types, composition_categories=cats),
              open(TISSUE_DIR / f"{name}_summary.json", "w"), indent=1)
    pd.DataFrame([t3]).to_csv(TISSUE_DIR / f"t3_{name}_summary.csv", index=False)

    # figure
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    labs = ["composition", "expression\n(all types)", "blood\ncomposition", "blood\nexpression"]
    vals = [c_r2, e_r2, BLOOD_COMP_R2, BLOOD_EXPR_R2]
    cols = ["tab:orange", "tab:blue", "0.7", "0.55"]
    ax.bar(range(4), vals, color=cols)
    ax.set_xticks(range(4)); ax.set_xticklabels(labs, fontsize=8)
    ax.set_ylabel(f"held-out R² ({metric})")
    ax.axhline(0, c="k", lw=0.5)
    ax.set_title(f"{name}: composition vs expression  (expr/comp={ratio:+.2f})" if np.isfinite(ratio)
                 else f"{name}: composition vs expression")
    fig.tight_layout()
    fig.savefig(TISSUE_FIG / f"t3_{name}_comp_vs_expr.png", dpi=140)
    plt.close(fig)

    # donor age histogram
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    ax.hist(donors.age, bins=20, color="steelblue")
    ax.set_xlabel("donor age (years)"); ax.set_ylabel("n donors")
    ax.set_title(f"{name}: {len(donors)} adult donors, {donors.age.min():.0f}–{donors.age.max():.0f} y")
    fig.tight_layout(); fig.savefig(TISSUE_FIG / f"t1_{name}_age_hist.png", dpi=130); plt.close(fig)

    log(f"\n[{name}] done.")
    return dict(t2=t2, t3=t3)


def main(names=None):
    names = names or sys.argv[1:]
    if not names:
        raise SystemExit("usage: python src/tissue_analyze.py NAME [NAME ...]")
    out = []
    for name in names:
        log = Logger(TISSUE_DIR / f"{name}_analyze.txt")
        try:
            out.append(run_one(name, log))
        finally:
            log.close()
    # combined T2 / T3 tables (merge with any previously written rows for other tissues)
    t2 = pd.DataFrame([o["t2"] for o in out])
    t3 = pd.DataFrame([o["t3"] for o in out])
    for fname, add in (("t2_batch_audit.csv", t2), ("t3_turnover.csv", t3)):
        p = TISSUE_DIR / fname
        if p.exists():
            old = pd.read_csv(p)
            add = pd.concat([old[~old.name.isin(add.name)], add], ignore_index=True)
        add.to_csv(p, index=False)
    return out


if __name__ == "__main__":
    main()

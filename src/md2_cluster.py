"""Cluster GM00731 at k=3/8/15 and Louvain; Louvain independently on GM23815. No age scores."""
from __future__ import annotations

import gc
import sys
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md2_common import (  # noqa: E402
    MD2_DIR, MD2_PROC, MD2_SEED, KMEANS_KS, CLUSTER_NPC, CLUSTER_N_INIT,
    LOUVAIN_RES, LOUVAIN_K_PARAM, LOUVAIN_PRUNE, LOUVAIN_NPC, LOUVAIN_N_HVG,
    LOUVAIN_MIN_GENES, LOUVAIN_MIN_CELLS_PER_GENE,
    CLUSTER_SETTINGS, CLUSTER_SETTINGS_FLAG, DECLARED_BEFORE_SCORES_FLAG,
    FIBRO3_DIR, AGED_LINE, YOUNG_LINE, DAYS,
    StopStep, Logger, dump_json, jsonable, md2_log_banner,
    load_manifest, save_manifest, record_failure, log_columns,
    progress_snapshot,
)
from fibro_common import FIBRO_DIR, FIBRO_RAW  # noqa: E402
from fibro2_common import load_frozen_ruler  # noqa: E402
from fibro2_tb import mito_mask  # noqa: E402
from fibro_stage2 import (  # noqa: E402
    find_10x_h5, load_one_h5, attach_meta, align_counts_to_ruler, _row_sum,
)
from md2_score import lognormalize_csr  # noqa: E402


def _require_cluster_settings():
    if not CLUSTER_SETTINGS_FLAG.exists():
        raise StopStep("prereg", "CLUSTER_SETTINGS.flag missing — settings must be frozen before clustering")
    existing = CLUSTER_SETTINGS_FLAG.read_text(encoding="utf-8")
    if existing.strip() != CLUSTER_SETTINGS.strip():
        raise StopStep("prereg", "CLUSTER_SETTINGS.flag does not match the frozen block. Not rewriting.")


def load_libs(log, cell_line=None):
    extract = FIBRO_RAW / "RAW"
    if not extract.exists():
        raise StopStep("cluster", f"missing GSE297234 extract {extract}")
    samp_path = FIBRO_DIR / "stage2_geo_samples.csv"
    if not samp_path.exists():
        raise StopStep("cluster", f"missing {samp_path} — Stage 2 sample table is required and is not re-fetched")
    samples = pd.read_csv(samp_path)
    log(f"[cluster] Stage 2 sample table columns actually read: {list(samples.columns)} path={samp_path}")
    log_columns("cluster_stage2_geo_samples", list(samples.columns), str(samp_path))
    need = {"gsm", "cell_line", "day", "age_years", "age_source"}
    missing = sorted(need - set(samples.columns))
    if missing:
        raise StopStep("cluster", f"{samp_path} missing columns {missing}. Not substituting.")
    h5s = find_10x_h5(extract, log)
    if not h5s:
        raise StopStep("cluster", "no filtered_feature_bc_matrix.h5 under GSE297234 RAW")
    libs = [load_one_h5(h, log) for h in h5s]
    libs = attach_meta(libs, samples, log)
    if cell_line is not None:
        libs = [L for L in libs if str(L["cell_line"]).upper() == str(cell_line).upper()]
        if not libs:
            raise StopStep("cluster", f"no libraries for {cell_line}")
    days = sorted(int(L["day"]) for L in libs)
    log(f"[cluster] loaded {len(libs)} libraries days={days} line_filter={cell_line}")
    return libs, samples


def _qc_cells(lib):
    X = lib["X"].tocsr()
    umi = np.asarray(X.sum(axis=1)).ravel()
    ngenes = np.asarray((X > 0).sum(axis=1)).ravel()
    mito, mito_meta = mito_mask(lib["var"])
    if mito.any():
        mito_umi = np.asarray(X[:, mito].sum(axis=1)).ravel()
        frac = mito_umi / np.where(umi > 0, umi, np.nan)
    else:
        frac = np.full(X.shape[0], np.nan)
    return umi, ngenes, frac, mito_meta


def concat_donor(libs_d, log):
    """Concatenate libraries of one donor in day order 0,3,7,10. Keep full X and ruler Y."""
    libs_d = sorted(libs_d, key=lambda L: int(L["day"]))
    line = str(libs_d[0]["cell_line"]).upper()
    days = [int(L["day"]) for L in libs_d]
    if days != list(DAYS):
        raise StopStep("cluster", f"{line}: days={days} expected {list(DAYS)}")
    frozen = load_frozen_ruler()
    X_parts, Y_parts, obs_parts = [], [], []
    mito_meta = None
    var0 = None
    for lib in libs_d:
        X = lib["X"].tocsr()
        if var0 is None:
            var0 = lib["var"].copy()
        else:
            if len(lib["var"]) != len(var0):
                raise StopStep("cluster", f"{line} d{lib['day']}: var length mismatch")
        Y, rec, _idx = align_counts_to_ruler(X, lib["var"], frozen, log)
        umi, ngenes, mito_frac, mito_meta = _qc_cells(lib)
        obs = lib["obs"].copy()
        obs["cell_line"] = line
        obs["day"] = int(lib["day"])
        obs["age_years"] = int(lib["age_years"])
        obs["age_source"] = str(lib["obs"]["age_source"].iloc[0])
        obs["gsm"] = str(lib["obs"]["gsm"].iloc[0])
        obs["umi"] = umi
        obs["n_genes"] = ngenes
        obs["mito_frac"] = mito_frac
        if len(obs) != X.shape[0]:
            raise StopStep("cluster", f"{line} d{lib['day']}: obs n={len(obs)} X n={X.shape[0]}")
        X_parts.append(X)
        Y_parts.append(Y)
        obs_parts.append(obs)
        log(f"[concat] {line} d{lib['day']} n_cells={X.shape[0]} ruler_overlap={rec.get('n_overlap')}")
    Xall = sparse.vstack(X_parts, format="csr")
    Yall = sparse.vstack(Y_parts, format="csr")
    obs = pd.concat(obs_parts, axis=0, ignore_index=True)
    log(f"[concat] {line} joint n_cells={len(obs)} n_genes={Xall.shape[1]}")
    return dict(
        cell_line=line, X=Xall, Y=Yall, obs=obs, var=var0,
        mito_meta=mito_meta or {}, frozen=frozen,
    )


def reuse_k3_labels(pack, log):
    path = FIBRO3_DIR / "t1_cell_labels.csv"
    if not path.exists():
        raise StopStep("k3", f"missing FIBRO3 labels {path}")
    lab = pd.read_csv(path)
    log(f"[k3] FIBRO3 t1_cell_labels columns actually read: {list(lab.columns)}")
    log_columns("fibro3_t1_cell_labels", list(lab.columns), str(path))
    need = {"barcode", "gsm", "cell_line", "day", "cluster"}
    missing = sorted(need - set(lab.columns))
    if missing:
        raise StopStep("k3", f"{path.name} missing {missing}")
    sub = lab[lab.cell_line.astype(str).str.upper() == pack["cell_line"]].copy()
    if sub.empty:
        raise StopStep("k3", f"no {pack['cell_line']} rows in {path}")
    obs = pack["obs"]
    key_obs = obs["gsm"].astype(str) + "||" + obs["barcode"].astype(str)
    key_lab = sub["gsm"].astype(str) + "||" + sub["barcode"].astype(str)
    pos = {k: int(c) for k, c in zip(key_lab, sub["cluster"].astype(int))}
    missing_keys = [k for k in key_obs if k not in pos]
    if missing_keys:
        raise StopStep(
            "k3",
            f"{len(missing_keys)} cells in loaded {pack['cell_line']} not in FIBRO3 labels "
            f"e.g. {missing_keys[:5]}. Not matching by order.",
        )
    extra = int(len(pos) - len(set(key_obs)))
    clusters = np.array([pos[k] for k in key_obs], dtype=int)
    log(f"[k3] reused FIBRO3 labels n={len(clusters)} n_clusters={int(np.unique(clusters).size)} "
        f"extra_in_fibro3={extra}")
    return clusters, dict(
        k=3, n_cells=int(len(clusters)), method="reuse_FIBRO3_t1_cell_labels",
        source=str(path), sizes=np.bincount(clusters).tolist(),
        extra_fibro3_rows_not_in_loaded=extra,
    )


def _gene_expression_mask(var):
    if "feature_type" not in var.columns:
        return np.ones(len(var), dtype=bool)
    ft = var["feature_type"].astype(str)
    return (ft.eq("Gene Expression") | ft.str.lower().eq("gene expression")).to_numpy()


def _hvg_quadratic(X, n_top, log):
    """Quadratic log10(mean)–log10(var) residual ranking. Deviation from Seurat vst loess."""
    X = X.tocsr()
    mean = np.asarray(X.mean(axis=0)).ravel()
    X2 = X.copy()
    X2.data = np.asarray(X2.data, np.float64) ** 2
    mean_sq = np.asarray(X2.mean(axis=0)).ravel()
    var = np.maximum(mean_sq - mean ** 2, 0.0)
    m = (mean > 0) & np.isfinite(var) & (var > 0)
    idx = np.flatnonzero(m)
    if idx.size < int(n_top):
        raise StopStep("hvg", f"genes with mean>0 and var>0 = {idx.size} < n_hvg={n_top}")
    lx = np.log10(mean[idx])
    ly = np.log10(var[idx])
    coef = np.polyfit(lx, ly, 2)
    resid = ly - np.polyval(coef, lx)
    order = idx[np.argsort(resid)[::-1]]
    take = order[:int(n_top)]
    log(f"[hvg] quadratic log-mean/log-var residual n_top={n_top} "
        f"n_candidates={idx.size} coef={coef.tolist()} "
        "(Seurat vst uses loess span=0.3; this is the reported deviation)")
    return take, dict(method="quadratic_logmean_logvar_residual", n_top=int(n_top),
                      n_candidates=int(idx.size), coef=list(map(float, coef)))


def _scale_residualize_mt(logX_hvg, pct_mt):
    """ScaleData-style: residualize each HVG on percent.mt, then z-score."""
    A = np.column_stack([np.ones(len(pct_mt)), np.asarray(pct_mt, float)])
    if not np.isfinite(A).all():
        raise StopStep("scale", "percent.mt has non-finite values")
    beta, *_ = np.linalg.lstsq(A, logX_hvg, rcond=None)
    resid = logX_hvg - A @ beta
    mu = resid.mean(axis=0)
    sd = resid.std(axis=0)
    sd = np.where(sd < 1e-12, 1.0, sd)
    z = (resid - mu) / sd
    return np.asarray(z, np.float32)


def louvain_one_donor(pack, log):
    """LogNormalize + HVG + percent.mt residualization + PCA + SNN + Louvain. Reported SCT deviation."""
    obs = pack["obs"].copy()
    var = pack["var"]
    X = pack["X"].tocsr()
    ge = _gene_expression_mask(var)
    X = X[:, np.asarray(ge, dtype=bool)]
    var = var.iloc[np.flatnonzero(np.asarray(ge, dtype=bool))].reset_index(drop=True)
    n_genes_cell = np.asarray((X > 0).sum(axis=1)).ravel()
    keep_cells = n_genes_cell > int(LOUVAIN_MIN_GENES)
    n_drop_cells = int((~keep_cells).sum())
    X = X[keep_cells]
    obs = obs.iloc[np.flatnonzero(keep_cells)].reset_index(drop=True)
    n_cells_gene = np.asarray((X > 0).sum(axis=0)).ravel()
    keep_genes = n_cells_gene >= int(LOUVAIN_MIN_CELLS_PER_GENE)
    n_drop_genes = int((~keep_genes).sum())
    X = X[:, keep_genes]
    var = var.iloc[np.flatnonzero(keep_genes)].reset_index(drop=True)
    log(f"[louvain-filter] {pack['cell_line']} drop_cells n_genes<={LOUVAIN_MIN_GENES}: {n_drop_cells} "
        f"kept_cells={int(keep_cells.sum())} drop_genes in<{LOUVAIN_MIN_CELLS_PER_GENE}_cells: {n_drop_genes} "
        f"kept_genes={int(keep_genes.sum())}")
    mito, mito_meta = mito_mask(var)
    umi = np.asarray(X.sum(axis=1)).ravel()
    if mito.any():
        mito_umi = np.asarray(X[:, mito].sum(axis=1)).ravel()
        pct_mt = 100.0 * mito_umi / np.where(umi > 0, umi, np.nan)
    else:
        raise StopStep("louvain", f"{pack['cell_line']}: zero MT- genes; cannot regress percent.mt")
    if not np.isfinite(pct_mt).all():
        raise StopStep("louvain", f"{pack['cell_line']}: non-finite percent.mt")
    logX = lognormalize_csr(X)
    hvg_idx, hvg_meta = _hvg_quadratic(X, LOUVAIN_N_HVG, log)
    log_hvg = logX[:, hvg_idx].astype(np.float32).toarray()
    Z = _scale_residualize_mt(log_hvg, pct_mt)
    del log_hvg
    gc.collect()
    npc = int(min(LOUVAIN_NPC, Z.shape[0] - 1, Z.shape[1]))
    pca = PCA(n_components=npc, svd_solver="randomized", random_state=int(MD2_SEED))
    P = pca.fit_transform(Z).astype(np.float32)
    evr = float(np.sum(pca.explained_variance_ratio_))
    del Z
    gc.collect()
    log(f"[louvain-pca] {pack['cell_line']} npc={npc} evr_sum={evr:.4f} P={P.shape}")
    nn = NearestNeighbors(n_neighbors=int(LOUVAIN_K_PARAM), algorithm="auto", metric="euclidean")
    nn.fit(P)
    _, idx = nn.kneighbors(P)
    n = int(P.shape[0])
    k = int(LOUVAIN_K_PARAM)
    rows = np.repeat(np.arange(n), k)
    cols = idx.ravel()
    M = sparse.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)), shape=(n, n))
    inter = (M @ M.T).tocoo()
    union = (2.0 * k) - inter.data
    jacc = inter.data / np.where(union > 0, union, 1.0)
    keep_e = (jacc >= float(LOUVAIN_PRUNE)) & (inter.row < inter.col)
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_weighted_edges_from(
        zip(inter.row[keep_e].tolist(), inter.col[keep_e].tolist(), jacc[keep_e].tolist())
    )
    log(f"[louvain-snn] {pack['cell_line']} n={n} edges={G.number_of_edges()} "
        f"k.param={k} prune={LOUVAIN_PRUNE} resolution={LOUVAIN_RES}")
    comms = nx.community.louvain_communities(
        G, weight="weight", resolution=float(LOUVAIN_RES), seed=int(MD2_SEED),
    )
    lab = np.full(n, -1, dtype=int)
    comms_sorted = sorted(comms, key=lambda s: (-len(s), min(s) if s else 0))
    for i, members in enumerate(comms_sorted):
        for v in members:
            lab[int(v)] = i
    if (lab < 0).any():
        raise StopStep("louvain", f"{pack['cell_line']}: {int((lab < 0).sum())} cells unassigned")
    n_cl = int(lab.max()) + 1
    sizes = np.bincount(lab, minlength=n_cl).tolist()
    log(f"[louvain] {pack['cell_line']} n_clusters={n_cl} sizes={sizes}")
    obs = obs.copy()
    obs["cluster"] = lab
    obs["louvain_kept"] = True
    obs["percent_mt"] = pct_mt
    meta = dict(
        method="LogNormalize+HVG_quadratic+ScaleData_percent.mt+PCA+SNN+networkx_louvain",
        deviation_from_STAR=(
            "STAR Methods specify Seurat v5.0.0 sctransform v2. This environment has no R/"
            "sctransform. LogNormalize scale.factor=10000 (Seurat NormalizeData default; "
            "paper public-data prose multiplies by median total UMIs — that prose is not used). "
            "HVG: quadratic log-mean/log-var residual ranking, n=3000 (Seurat SCTransform "
            "variable.features.n default; STAR Methods do not state n HVGs; Seurat vst uses loess). "
            f"PCA {npc} PCs (Seurat RunPCA default). SNN k.param={LOUVAIN_K_PARAM} "
            f"prune.SNN={LOUVAIN_PRUNE} (Seurat FindNeighbors defaults). "
            f"Louvain resolution={LOUVAIN_RES} (Seurat FindClusters default; STAR Methods omit resolution)."
        ),
        n_cells_in=int(pack["X"].shape[0]),
        n_cells_kept=int(n),
        n_drop_cells=n_drop_cells,
        n_genes_kept=int(X.shape[1]),
        n_drop_genes=n_drop_genes,
        n_clusters=n_cl,
        sizes=sizes,
        resolution=float(LOUVAIN_RES),
        k_param=int(LOUVAIN_K_PARAM),
        prune_snn=float(LOUVAIN_PRUNE),
        npc=npc,
        n_hvg=int(LOUVAIN_N_HVG),
        pca_explained_variance_ratio_sum=evr,
        random_state=MD2_SEED,
        mito_rule=mito_meta.get("rule"),
        n_mito_genes=mito_meta.get("n_mito"),
        hvg=hvg_meta,
        min_genes=int(LOUVAIN_MIN_GENES),
        min_cells_per_gene=int(LOUVAIN_MIN_CELLS_PER_GENE),
    )
    return obs, lab, meta, dict(X=X, var=var, logX=logX, obs=obs, pct_mt=pct_mt)


def _save_sparse(path, X, var):
    X = X.tocsr()
    np.savez_compressed(
        path,
        data=X.data, indices=X.indices, indptr=X.indptr, shape=np.array(X.shape),
        symbols=var["symbol"].astype(str).to_numpy(),
        gene_id=var["gene_id"].astype(str).to_numpy() if "gene_id" in var.columns else np.array([], dtype=object),
        feature_type=var["feature_type"].astype(str).to_numpy() if "feature_type" in var.columns else np.array([], dtype=object),
    )


def run_clustering(log=None):
    _require_cluster_settings()
    close_log = False
    if log is None:
        log = Logger(MD2_DIR / "cluster_report.txt")
        close_log = True
    md2_log_banner(log, "CLUSTER")
    log(CLUSTER_SETTINGS)
    if DECLARED_BEFORE_SCORES_FLAG.exists():
        log(f"[prereg] {DECLARED_BEFORE_SCORES_FLAG.name} already exists; clustering may be a rerun")

    libs, samples = load_libs(log, cell_line=None)
    by_line = {}
    for lib in libs:
        by_line.setdefault(str(lib["cell_line"]).upper(), []).append(lib)

    meta_all = {}
    if AGED_LINE not in by_line:
        raise StopStep("cluster", f"donor {AGED_LINE} missing")
    pack_a = concat_donor(by_line[AGED_LINE], log)
    labels_k3, meta_k3 = reuse_k3_labels(pack_a, log)
    obs_k3 = pack_a["obs"].copy()
    obs_k3["cluster"] = labels_k3
    obs_k3["resolution"] = "k3"
    obs_k3.to_csv(MD2_DIR / "cluster_labels_k3_GM00731.csv", index=False)
    meta_all["k3_GM00731"] = meta_k3

    Yall = pack_a["Y"]
    n = int(Yall.shape[0])
    colsum = _row_sum(Yall)
    keep_idx = np.flatnonzero(colsum > 0)
    n_keep = int(keep_idx.size)
    npc = int(min(CLUSTER_NPC, n - 1, n_keep))
    log(f"[kmeans-shared-pca] {AGED_LINE} n_cells={n} n_genes_kept={n_keep} npc={npc}")
    Yk = Yall.tocsc()[:, keep_idx].tocsr().astype(np.float32)
    Z = Yk.toarray()
    del Yk
    gc.collect()
    libsz = Z.sum(1, keepdims=True)
    np.maximum(libsz, 1.0, out=libsz)
    Z /= libsz
    Z *= 1e4
    np.log1p(Z, out=Z)
    pca = PCA(n_components=npc, svd_solver="randomized", random_state=int(MD2_SEED))
    P = pca.fit_transform(Z).astype(np.float32)
    evr = float(np.sum(pca.explained_variance_ratio_))
    del Z
    gc.collect()
    log(f"[kmeans-shared-pca] done evr_sum={evr:.4f}")
    np.savez_compressed(MD2_PROC / "pca_kmeans_GM00731.npz", P=P, evr=np.array([evr]))
    for k in (8, 15):
        km = KMeans(n_clusters=int(k), random_state=int(MD2_SEED), n_init=CLUSTER_N_INIT)
        lab = km.fit_predict(np.asarray(P, np.float32)).astype(int)
        sizes = np.bincount(lab, minlength=int(k)).tolist()
        log(f"[kmeans] {AGED_LINE} k={k} sizes={sizes} inertia={float(km.inertia_):.4g}")
        obs_k = pack_a["obs"].copy()
        obs_k["cluster"] = lab
        obs_k["resolution"] = f"k{k}"
        obs_k.to_csv(MD2_DIR / f"cluster_labels_k{k}_GM00731.csv", index=False)
        meta_all[f"k{k}_GM00731"] = dict(
            k=int(k), n_cells=n, n_genes_kept=n_keep, npc=npc, sizes=sizes,
            method="PCA_randomized+KMeans", pca_explained_variance_ratio_sum=evr,
            random_state=MD2_SEED, n_init=CLUSTER_N_INIT,
        )
    del P
    gc.collect()

    obs_lv_a, lab_lv_a, meta_lv_a, packed_lv_a = louvain_one_donor(pack_a, log)
    obs_lv_a.to_csv(MD2_DIR / "cluster_labels_louvain_GM00731.csv", index=False)
    meta_all["louvain_GM00731"] = meta_lv_a
    _save_sparse(MD2_PROC / "louvain_counts_GM00731.npz", packed_lv_a["X"], packed_lv_a["var"])
    packed_lv_a["obs"].to_csv(MD2_DIR / "louvain_obs_GM00731.csv", index=False)

    if YOUNG_LINE not in by_line:
        raise StopStep("cluster", f"donor {YOUNG_LINE} missing")
    pack_y = concat_donor(by_line[YOUNG_LINE], log)
    obs_lv_y, lab_lv_y, meta_lv_y, packed_lv_y = louvain_one_donor(pack_y, log)
    obs_lv_y.to_csv(MD2_DIR / "cluster_labels_louvain_GM23815.csv", index=False)
    meta_all["louvain_GM23815"] = meta_lv_y
    _save_sparse(MD2_PROC / "louvain_counts_GM23815.npz", packed_lv_y["X"], packed_lv_y["var"])
    packed_lv_y["obs"].to_csv(MD2_DIR / "louvain_obs_GM23815.csv", index=False)

    # Full (unfiltered) counts for k-means AddModuleScore on aged donor.
    ge = _gene_expression_mask(pack_a["var"])
    _save_sparse(
        MD2_PROC / "allcell_counts_GM00731.npz",
        pack_a["X"][:, ge],
        pack_a["var"].iloc[np.flatnonzero(ge)].reset_index(drop=True),
    )
    pack_a["obs"].to_csv(MD2_DIR / "allcell_obs_GM00731.csv", index=False)
    # ruler-aligned Y for frozen age pseudobulk
    Ya = pack_a["Y"].tocsr()
    np.savez_compressed(
        MD2_PROC / "rulerY_GM00731.npz",
        data=Ya.data, indices=Ya.indices, indptr=Ya.indptr, shape=np.array(Ya.shape),
    )
    Yy = pack_y["Y"].tocsr()
    np.savez_compressed(
        MD2_PROC / "rulerY_GM23815.npz",
        data=Yy.data, indices=Yy.indices, indptr=Yy.indptr, shape=np.array(Yy.shape),
    )
    ge_y = _gene_expression_mask(pack_y["var"])
    _save_sparse(
        MD2_PROC / "allcell_counts_GM23815.npz",
        pack_y["X"][:, ge_y],
        pack_y["var"].iloc[np.flatnonzero(ge_y)].reset_index(drop=True),
    )
    pack_y["obs"].to_csv(MD2_DIR / "allcell_obs_GM23815.csv", index=False)

    dump_json(MD2_DIR / "cluster_meta.json", jsonable(meta_all))

    ncl = {k: (v.get("n_clusters") if v.get("n_clusters") is not None else v.get("k")) for k, v in meta_all.items()}
    if not DECLARED_BEFORE_SCORES_FLAG.exists():
        DECLARED_BEFORE_SCORES_FLAG.write_text(
            CLUSTER_SETTINGS
            + "\n\nThis flag was written after clustering (k=3 reuse, k=8, k=15, Louvain GM00731, "
              "Louvain GM23815) and before any age score or MD score.\n"
            + "\nn_clusters: " + str(ncl) + "\n",
            encoding="utf-8",
        )
        log(f"[prereg] wrote {DECLARED_BEFORE_SCORES_FLAG} after clustering, before scoring")
    else:
        log(f"[prereg] {DECLARED_BEFORE_SCORES_FLAG.name} already exists; not rewriting")

    man = load_manifest()
    man["status"] = "CLUSTER_DONE"
    man["n_clusters"] = ncl
    save_manifest(man)
    progress_snapshot(
        "score MD/TGFB/age on frozen clusters (DECLARED_BEFORE_SCORES written)",
        stop="CLUSTER_DONE",
    )
    if close_log:
        log.close()
    return meta_all


if __name__ == "__main__":
    try:
        run_clustering()
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise

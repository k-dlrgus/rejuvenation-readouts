"""Extrapolation check: distance of each scored matrix to GTEx fibroblast training z-space.

Report-only. No pass/fail.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fibro2_common import (  # noqa: E402
    FIBRO2_DIR, FIBRO2_PROC, FIBRO2_SEED, PCA_MAHAL_K,
    StopStep, Logger, dump_json, jsonable, fibro2_log_banner,
    load_manifest, save_manifest, record_failure, load_frozen_ruler, progress_snapshot,
)
from fibro_stage1 import load_pack  # noqa: E402
from target_common import zscore_train  # noqa: E402


def _pca_mahal_fit(Xz, k, seed):
    Xz = np.asarray(Xz, np.float64)
    n, p = Xz.shape
    k = int(min(k, n - 1, p))
    mu = Xz.mean(0)
    Xc = Xz - mu
    # economy SVD
    rng = np.random.default_rng(int(seed))
    # deterministic SVD
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    V = Vt[:k].T
    pcs = Xc @ V
    var = pcs.var(0)
    var = np.where(var < 1e-12, 1.0, var)
    return dict(mu=mu, V=V, var=var, k=k)


def _nn_and_mahal(Zq, Xz, pca):
    Zq = np.asarray(Zq, np.float64)
    Xz = np.asarray(Xz, np.float64)
    q2 = (Zq * Zq).sum(1, keepdims=True)
    x2 = (Xz * Xz).sum(1, keepdims=True).T
    d2 = np.maximum(q2 + x2 - 2.0 * (Zq @ Xz.T), 0.0)
    nn_idx = d2.argmin(1)
    nn = np.sqrt(d2[np.arange(len(Zq)), nn_idx])
    Zc = Zq - pca["mu"]
    pcs = Zc @ pca["V"]
    mahal = np.sqrt(((pcs ** 2) / pca["var"]).sum(1))
    return nn, nn_idx, mahal


def _gtex_z(pack, frozen):
    X = np.asarray(pack["X"], np.float64)
    mu = np.asarray(frozen["mu"], float)
    sd = np.asarray(frozen["sd"], float)
    sd = np.where(sd < 1e-12, 1.0, sd)
    return (X - mu) / sd


def run_extrap(log=None):
    close_log = False
    if log is None:
        log = Logger(FIBRO2_DIR / "extrap_report.txt")
        close_log = True
    fibro2_log_banner(log, "EXTRAP")
    log("Report-only. No pass/fail. Mahalanobis is in the leading PCA of GTEx fibroblast z-space. "
        "Nearest neighbour is Euclidean in the same z-space (missing genes at 0).")
    frozen = load_frozen_ruler()
    pack = load_pack()
    Xz = _gtex_z(pack, frozen)
    log(f"[extrap] GTEx train z: n={Xz.shape[0]} p={Xz.shape[1]}")
    pca = _pca_mahal_fit(Xz, PCA_MAHAL_K, FIBRO2_SEED)
    log(f"[extrap] PCA k={pca['k']}")
    rows = []

    def add_block(Z, meta_df, cohort, variant):
        nn, nn_idx, mahal = _nn_and_mahal(Z, Xz, pca)
        for i in range(len(meta_df)):
            rec = dict(
                cohort=cohort, variant=variant,
                nn_euclidean=float(nn[i]), nn_train_index=int(nn_idx[i]),
                mahalanobis_pca=float(mahal[i]), pca_k=int(pca["k"]),
            )
            for c in meta_df.columns:
                rec[c] = meta_df.iloc[i][c]
            rows.append(rec)

    # GSE297234 all-cell
    p_a = FIBRO2_DIR / "tb_allcell_z.npz"
    if p_a.exists():
        z = np.load(p_a, allow_pickle=True)
        meta = pd.read_csv(FIBRO2_DIR / "tb_allcell_scores.csv")
        add_block(z["Z"], meta, "GSE297234", "all_cell")
        log(f"[extrap] GSE297234 all-cell n={len(meta)}")
    else:
        log("[extrap] tb_allcell_z.npz missing — T-B not run")

    p_b = FIBRO2_DIR / "tb_cluster_z.npz"
    if p_b.exists() and (FIBRO2_DIR / "tb_cluster_scores.csv").exists():
        z = np.load(p_b, allow_pickle=True)
        meta = pd.read_csv(FIBRO2_DIR / "tb_cluster_scores.csv")
        add_block(z["Z"], meta, "GSE297234", "cluster")
        log(f"[extrap] GSE297234 cluster n={len(meta)}")

    # T-A projections
    for npz in sorted(FIBRO2_PROC.glob("ta_*_proj.npz")):
        tag = npz.name.replace("ta_", "").replace("_proj.npz", "")
        scores_p = FIBRO2_DIR / f"ta_{tag}_scores.csv"
        if not scores_p.exists():
            log(f"[extrap] skip {tag}: no scores csv")
            continue
        z = np.load(npz, allow_pickle=True)
        meta = pd.read_csv(scores_p)
        add_block(z["Z"], meta, tag, "external")
        log(f"[extrap] {tag} n={len(meta)}")

    if not rows:
        log("[extrap] no scored matrices yet")
        dump_json(FIBRO2_DIR / "extrap_summary.json", jsonable(dict(n=0, note="no scored matrices")))
        if close_log:
            log.close()
        return dict(n=0)

    df = pd.DataFrame(rows)
    df.to_csv(FIBRO2_DIR / "extrap_distances.csv", index=False)
    # cohort summary
    summ_rows = []
    for (cohort, variant), g in df.groupby(["cohort", "variant"]):
        summ_rows.append(dict(
            cohort=cohort, variant=variant, n=int(len(g)),
            nn_euclidean_min=float(g.nn_euclidean.min()),
            nn_euclidean_median=float(g.nn_euclidean.median()),
            nn_euclidean_max=float(g.nn_euclidean.max()),
            mahalanobis_pca_min=float(g.mahalanobis_pca.min()),
            mahalanobis_pca_median=float(g.mahalanobis_pca.median()),
            mahalanobis_pca_max=float(g.mahalanobis_pca.max()),
            pca_k=int(pca["k"]),
        ))
    summ = pd.DataFrame(summ_rows)
    summ.to_csv(FIBRO2_DIR / "extrap_summary_by_cohort.csv", index=False)
    dump_json(FIBRO2_DIR / "extrap_summary.json", jsonable(dict(
        n_rows=int(len(df)), pca_k=int(pca["k"]),
        gtex_n=int(Xz.shape[0]), gtex_p=int(Xz.shape[1]),
        note="report-only; no pass/fail",
        by_cohort=summ.to_dict(orient="records"),
    )))
    man = load_manifest()
    man["status"] = "EXTRAP_DONE"
    save_manifest(man)
    progress_snapshot("write FINDINGS_FIBRO2.md from disk", stop="EXTRAP_DONE")
    if close_log:
        log.close()
    return dict(n_rows=len(df))


if __name__ == "__main__":
    try:
        run_extrap()
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        raise

"""Task 1 — resolution sweep on GM00731. Scores only after DECLARED_BEFORE_SCORES."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md2_common import (  # noqa: E402
    MD2_DIR, MD2_PROC, MD2_SEED, MIN_CELLS, AGED_LINE, N_RANDOM_DIR,
    PREREG_TASK1, PREREG_TASK1_FLAG, DECLARED_BEFORE_SCORES_FLAG, CLUSTER_SETTINGS,
    FIBRO3_DIR, FIBRO3_TASK3_AGED_D0D7,
    StopStep, Logger, dump_json, jsonable, load_json, md2_log_banner,
    load_manifest, save_manifest, record_failure, log_columns,
    progress_snapshot, FROZEN_RULER,
)
from fibro2_common import load_frozen_ruler  # noqa: E402
from fibro_stage2 import random_directions  # noqa: E402
from fibro_common import N_RANDOM_DIR as _N  # noqa: E402
from trajectory_common import permutation_p  # noqa: E402
from md2_score import (  # noqa: E402
    lognormalize_csr, add_module_score, score_frozen_pseudobulk, map_features,
)


def _csr_from_npz(path):
    z = np.load(path, allow_pickle=True)
    X = sparse.csr_matrix(
        (z["data"], z["indices"], z["indptr"]),
        shape=tuple(int(x) for x in z["shape"]),
    )
    symbols = np.asarray(z["symbols"]).astype(str)
    return X, symbols, z


def _require():
    if not PREREG_TASK1_FLAG.exists():
        raise StopStep("prereg", "PREREG_TASK1.flag missing")
    if not DECLARED_BEFORE_SCORES_FLAG.exists():
        raise StopStep("prereg", "DECLARED_BEFORE_SCORES.flag missing — refusing to score")
    gs = MD2_DIR / "genesets.json"
    if not gs.exists():
        raise StopStep("task1", "genesets.json missing")
    if not FROZEN_RULER.exists():
        raise StopStep("frozen_ruler", f"missing {FROZEN_RULER}")


def _load_sets():
    rec = load_json(MD2_DIR / "genesets.json")
    return rec


def _pb_table(obs, Y, frozen, log, tag, ams_by_cell, min_cells=MIN_CELLS, reuse_age=None):
    """Cluster × day rows. Age via TMM+frozen ruler unless reuse_age DataFrame given."""
    rows, mats, idx_keep = [], [], []
    for (day, cl), g in obs.groupby(["day", "cluster"], sort=True):
        idx = g.index.to_numpy()
        n = int(len(idx))
        n_tp = int((obs.day == day).sum())
        rec = dict(
            cell_line=AGED_LINE, day=int(day), cluster=int(cl),
            n_cells=n, n_timepoint=n_tp,
            frac_timepoint=float(n / n_tp) if n_tp else np.nan,
            median_UMI=float(np.median(g.umi.to_numpy(float))) if n and "umi" in g.columns else np.nan,
            median_genes=float(np.median(g.n_genes.to_numpy(float))) if n and "n_genes" in g.columns else np.nan,
            mito_frac_median_cell=float(np.nanmedian(g.mito_frac.to_numpy(float))) if n and "mito_frac" in g.columns else np.nan,
            below_min_cells=bool(n < min_cells),
        )
        for key, arr in ams_by_cell.items():
            rec[key] = float(np.mean(arr[idx])) if n else np.nan
        rows.append(rec)
        if Y is not None:
            pb = np.asarray(Y[idx].sum(axis=0), dtype=np.float64).ravel()
            mats.append(pb)
            idx_keep.append(len(rows) - 1)
    df = pd.DataFrame(rows)
    scored = None
    if reuse_age is not None:
        m = reuse_age.copy()
        m["day"] = m["day"].astype(int)
        m["cluster"] = m["cluster"].astype(int)
        df = df.merge(
            m[["day", "cluster", "age_score", "pluri_with_OSKM", "pluri_without_OSKM", "pluri_primary"]],
            on=["day", "cluster"], how="left",
        )
        df["age_source"] = "FIBRO3_t1_cluster_table"
    elif mats:
        C = np.vstack(mats)
        scored = score_frozen_pseudobulk(C, frozen, log, tag=tag)
        df["age_score"] = scored["age"]
        df["pluri_with_OSKM"] = scored["pluri_with"]
        df["pluri_without_OSKM"] = scored["pluri_without"]
        df["pluri_primary"] = scored["pluri_with"]
        df["age_source"] = "TMM_log2CPM_frozen_ruler"
        df.loc[df.below_min_cells, ["age_score", "pluri_with_OSKM", "pluri_without_OSKM", "pluri_primary"]] = np.nan
        for k in list(ams_by_cell):
            df.loc[df.below_min_cells, k] = np.nan
    if reuse_age is not None:
        for k in list(ams_by_cell):
            df.loc[df.below_min_cells, k] = np.nan
    return df, scored


def _age_null_rows(df, scored, frozen, log, resolution):
    rows = []
    if scored is None:
        return rows
    Z = scored["Z"]
    age = np.asarray(scored["age"], float)
    if len(df) != len(Z):
        raise StopStep("task1", f"{resolution}: table n={len(df)} Z n={len(Z)}")
    null = random_directions(frozen["w"], Z, seed=MD2_SEED, n=N_RANDOM_DIR)
    meta = df.reset_index(drop=True)

    def pick(cl, day):
        m = (meta.cluster.astype(int) == int(cl)) & (meta.day.astype(int) == int(day))
        idx = np.flatnonzero(m.to_numpy())
        if len(idx) == 0:
            return None
        if len(idx) != 1:
            raise StopStep("task1", f"{resolution} cluster={cl} day={day} n_rows={len(idx)}")
        return int(idx[0])

    clusters = sorted(meta.cluster.astype(int).unique().tolist())
    for cl in clusters:
        for dend in (7, 10):
            i0 = pick(cl, 0)
            i1 = pick(cl, dend)
            rec = dict(
                resolution=resolution, cell_line=AGED_LINE, cluster=int(cl),
                endpoint=f"d0→d{dend}", d_end=int(dend),
                n_random=N_RANDOM_DIR, seed=MD2_SEED, min_cells=MIN_CELLS,
                instrument="frozen_ruler", null_kind="permuted_ruler_weights",
            )
            if i0 is None or i1 is None:
                rec.update(ok=False, reason="missing cluster row at an endpoint")
                rows.append(rec)
                continue
            n0 = int(meta.n_cells.iloc[i0])
            n1 = int(meta.n_cells.iloc[i1])
            rec["n_cells_d0"] = n0
            rec["n_cells_end"] = n1
            if n0 < MIN_CELLS or n1 < MIN_CELLS:
                rec.update(ok=False, reason=f"n_cells<{MIN_CELLS} at an endpoint")
                rows.append(rec)
                continue
            if bool(meta.below_min_cells.iloc[i0]) or bool(meta.below_min_cells.iloc[i1]):
                rec.update(ok=False, reason="below_min_cells at an endpoint")
                rows.append(rec)
                continue
            real = float(age[i0] - age[i1])
            null_d = np.asarray(null)[i0] - np.asarray(null)[i1]
            p = permutation_p(real, null_d, greater=True)
            n_ge = int(np.sum(null_d >= real))
            rec.update(
                ok=True, score_d0=float(age[i0]), score_end=float(age[i1]),
                decline=real, p=p, n_random_ge_real=n_ge,
                pass_p=bool(np.isfinite(p) and p <= 0.05),
            )
            rows.append(rec)
            log(f"[t1 null {resolution}] c{cl} d0→d{dend} decline={real:+.4f} p={p:+.3f} "
                f"n_ge={n_ge}/{N_RANDOM_DIR} n0={n0} n1={n1}")
    return rows


def _copy_k3_age_null(log):
    path = FIBRO3_DIR / "t3_null.csv"
    if not path.exists():
        raise StopStep("task1", f"missing {path} — k=3 age nulls are reused, not recomputed")
    t3 = pd.read_csv(path)
    log(f"[t1 k3] FIBRO3 t3_null columns actually read: {list(t3.columns)}")
    log_columns("fibro3_t3_null", list(t3.columns), str(path))
    sub = t3[t3.cell_line.astype(str).str.upper() == AGED_LINE].copy()
    if sub.empty:
        raise StopStep("task1", "FIBRO3 t3_null has no GM00731 rows")
    rows = []
    for _, r in sub.iterrows():
        rec = dict(
            resolution="k3", cell_line=AGED_LINE, cluster=int(r.cluster),
            endpoint=str(r.endpoint), d_end=int(r.d_end),
            n_random=int(r.n_random) if pd.notna(r.n_random) else N_RANDOM_DIR,
            seed=MD2_SEED, min_cells=MIN_CELLS,
            instrument="frozen_ruler", null_kind="permuted_ruler_weights_reused_from_FIBRO3",
            ok=bool(r.ok),
            n_cells_d0=float(r.n_cells_d0) if pd.notna(r.n_cells_d0) else np.nan,
            n_cells_end=float(r.n_cells_end) if pd.notna(r.n_cells_end) else np.nan,
        )
        if bool(r.ok):
            rec.update(
                score_d0=float(r.score_d0), score_end=float(r.score_end),
                decline=float(r.decline), p=float(r.p),
                n_random_ge_real=float(r.n_random_ge_real),
                pass_p=bool(r.pass_p),
            )
        else:
            rec["reason"] = str(r.reason) if pd.notna(r.reason) else "not ok"
        rows.append(rec)
        log(f"[t1 null k3-reuse] c{int(r.cluster)} {r.endpoint} ok={r.ok} "
            f"decline={r.decline if bool(r.ok) else 'NA'} p={r.p if bool(r.ok) else 'NA'}")
    return rows


def _reading(null_df):
    """Only the outcome that fired. Aged donor, all four resolutions."""
    g = null_df[(null_df.cell_line == AGED_LINE) & (null_df.instrument == "frozen_ruler")].copy()
    hits = g[(g.ok == True) & (g.pass_p == True)].copy()  # noqa: E712
    by_res = {}
    for res, sub in g.groupby("resolution"):
        ok = sub[sub.ok == True]  # noqa: E712
        by_res[str(res)] = dict(
            n_ok=int(len(ok)),
            n_pass=int((ok.pass_p == True).sum()) if len(ok) else 0,  # noqa: E712
            pass_rows=ok[ok.pass_p == True][  # noqa: E712
                ["cluster", "endpoint", "decline", "p", "n_cells_d0", "n_cells_end"]
            ].to_dict(orient="records") if len(ok) else [],
        )
    k3_pass = by_res.get("k3", {}).get("n_pass", 0)
    k8_pass = by_res.get("k8", {}).get("n_pass", 0)
    k15_pass = by_res.get("k15", {}).get("n_pass", 0)
    lv_pass = by_res.get("louvain", {}).get("n_pass", 0)
    finer = (k8_pass + k15_pass + lv_pass) > 0
    k3_any = k3_pass > 0
    finer_only_k15 = k15_pass > 0 and lv_pass == 0
    finer_only_lv = lv_pass > 0 and k15_pass == 0
    if finer_only_k15 or finer_only_lv:
        key = "resolution_dependent"
        text = (
            "Results appear at k=15 but not at Louvain, or vice versa. Report per resolution, "
            "do not average, and say the finding is resolution-dependent. "
            f"k15_n_pass={k15_pass} louvain_n_pass={lv_pass} k8_n_pass={k8_pass} k3_n_pass={k3_pass}."
        )
    elif finer and not k3_any:
        key = "resolution_artifact"
        text = (
            "At finer resolution some cluster shows a d0→d7 or d0→d10 age decline with p ≤ 0.05. "
            "FINDINGS_FIBRO3.md's per-cluster negative was a resolution artifact. "
            f"This supersedes the FIBRO3 per-cluster reading: \"{FIBRO3_TASK3_AGED_D0D7}\" "
            f"Hits: {hits[['resolution', 'cluster', 'endpoint', 'decline', 'p']].to_dict(orient='records') if len(hits) else []}."
        )
    elif (not finer) and (not k3_any):
        key = "negative_holds"
        text = (
            "No cluster at any resolution shows p ≤ 0.05. The per-cluster negative holds and is "
            "not a resolution artifact. The frozen ruler does not move in these cells at any "
            "granularity tested."
        )
    elif k3_any and finer:
        key = "resolution_artifact"
        text = (
            "At finer resolution some cluster shows a d0→d7 or d0→d10 age decline with p ≤ 0.05. "
            "FINDINGS_FIBRO3.md's per-cluster negative was a resolution artifact. "
            f"This supersedes the FIBRO3 per-cluster reading: \"{FIBRO3_TASK3_AGED_D0D7}\" "
            f"k3 also has pass rows. Hits: "
            f"{hits[['resolution', 'cluster', 'endpoint', 'decline', 'p']].to_dict(orient='records') if len(hits) else []}."
        )
    else:
        # k3 pass but finer none — still report per resolution, not averaged.
        key = "resolution_dependent"
        text = (
            "k=3 has p≤0.05 and finer resolutions do not, or the pattern is not the listed "
            "artifact/hold/split. Report per resolution, do not average. "
            f"k3_n_pass={k3_pass} k8_n_pass={k8_pass} k15_n_pass={k15_pass} louvain_n_pass={lv_pass}."
        )
    return dict(
        key=key, text=text, by_resolution=by_res,
        n_pass_total=int(len(hits)),
        k3_n_pass=int(k3_pass), k8_n_pass=int(k8_pass),
        k15_n_pass=int(k15_pass), louvain_n_pass=int(lv_pass),
        fibro3_sentence_at_issue=FIBRO3_TASK3_AGED_D0D7,
    )


def run_task1(log=None):
    _require()
    close_log = False
    if log is None:
        log = Logger(MD2_DIR / "t1_report.txt")
        close_log = True
    md2_log_banner(log, "TASK1")
    log(PREREG_TASK1)
    log(CLUSTER_SETTINGS)
    sets = _load_sets()
    frozen = load_frozen_ruler()
    gene_sets_ams = {
        "md_score": list(sets["MD"]),
        "tgfb_score": list(sets["TGFB"]),
    }
    for st, genes in (sets.get("reprog") or {}).items():
        gene_sets_ams[f"state_{st}"] = list(genes)

    X, symbols, _z = _csr_from_npz(MD2_PROC / "allcell_counts_GM00731.npz")
    obs = pd.read_csv(MD2_DIR / "allcell_obs_GM00731.csv")
    if len(obs) != X.shape[0]:
        raise StopStep("task1", f"allcell obs n={len(obs)} X n={X.shape[0]}")
    log(f"[t1] allcell GM00731 n_cells={len(obs)} n_genes={X.shape[1]}")
    logX = lognormalize_csr(X)
    ams, ams_meta = add_module_score(logX, symbols, gene_sets_ams, log=log, tag="GM00731_allcell")
    dump_json(MD2_DIR / "t1_ams_meta_allcell.json", jsonable({
        k: {kk: vv for kk, vv in v.items() if kk not in ("mapped",)}  # mapped is long; n_mapped kept
        if isinstance(v, dict) else v
        for k, v in (ams_meta.get("per_set") or {}).items()
    }))
    yz = np.load(MD2_PROC / "rulerY_GM00731.npz")
    Y = sparse.csr_matrix(
        (yz["data"], yz["indices"], yz["indptr"]),
        shape=tuple(int(x) for x in yz["shape"]),
    )
    if Y.shape[0] != len(obs):
        raise StopStep("task1", f"rulerY n={Y.shape[0]} obs n={len(obs)}")

    fibro3_tbl = pd.read_csv(FIBRO3_DIR / "t1_cluster_table.csv")
    log_columns("fibro3_t1_cluster_table", list(fibro3_tbl.columns), str(FIBRO3_DIR / "t1_cluster_table.csv"))
    fibro3_aged = fibro3_tbl[fibro3_tbl.cell_line.astype(str).str.upper() == AGED_LINE].copy()

    tables = []
    null_rows = []
    scored_by_res = {}

    for res, lab_path in (
        ("k3", MD2_DIR / "cluster_labels_k3_GM00731.csv"),
        ("k8", MD2_DIR / "cluster_labels_k8_GM00731.csv"),
        ("k15", MD2_DIR / "cluster_labels_k15_GM00731.csv"),
    ):
        lab = pd.read_csv(lab_path)
        log(f"[t1] {res} label columns actually read: {list(lab.columns)}")
        key_obs = obs["gsm"].astype(str) + "||" + obs["barcode"].astype(str)
        key_lab = lab["gsm"].astype(str) + "||" + lab["barcode"].astype(str)
        pos = {k: int(c) for k, c in zip(key_lab, lab["cluster"].astype(int))}
        miss = [k for k in key_obs if k not in pos]
        if miss:
            raise StopStep("task1", f"{res}: {len(miss)} cells missing from labels")
        obs_r = obs.copy()
        obs_r["cluster"] = [pos[k] for k in key_obs]
        reuse = fibro3_aged if res == "k3" else None
        df, scored = _pb_table(
            obs_r, Y, frozen, log, tag=f"t1_{res}", ams_by_cell=ams, reuse_age=reuse,
        )
        df.insert(0, "resolution", res)
        tables.append(df)
        if res == "k3":
            null_rows.extend(_copy_k3_age_null(log))
        else:
            null_rows.extend(_age_null_rows(df, scored, frozen, log, res))
            scored_by_res[res] = scored

    # Louvain: filtered cells, own AMS (gene means on the Louvain object).
    Xl, symbols_l, _ = _csr_from_npz(MD2_PROC / "louvain_counts_GM00731.npz")
    obs_l = pd.read_csv(MD2_DIR / "louvain_obs_GM00731.csv")
    if len(obs_l) != Xl.shape[0]:
        raise StopStep("task1", f"louvain obs n={len(obs_l)} X n={Xl.shape[0]}")
    logXl = lognormalize_csr(Xl)
    ams_l, ams_l_meta = add_module_score(logXl, symbols_l, gene_sets_ams, log=log, tag="GM00731_louvain")
    dump_json(MD2_DIR / "t1_ams_meta_louvain.json", jsonable({
        k: {kk: vv for kk, vv in v.items() if kk != "mapped"}
        for k, v in (ams_l_meta.get("per_set") or {}).items()
    }))
    # ruler Y restricted to louvain-kept cells
    key_all = obs["gsm"].astype(str) + "||" + obs["barcode"].astype(str)
    pos_all = {k: i for i, k in enumerate(key_all)}
    key_l = obs_l["gsm"].astype(str) + "||" + obs_l["barcode"].astype(str)
    take = np.array([pos_all[k] for k in key_l], dtype=int)
    Yl = Y[take]
    df_l, scored_l = _pb_table(obs_l, Yl, frozen, log, tag="t1_louvain", ams_by_cell=ams_l)
    df_l.insert(0, "resolution", "louvain")
    tables.append(df_l)
    null_rows.extend(_age_null_rows(df_l, scored_l, frozen, log, "louvain"))
    scored_by_res["louvain"] = scored_l

    # persist louvain AMS bins for Task 2 MD null
    np.savez_compressed(
        MD2_PROC / "louvain_ams_GM00731.npz",
        md=ams_l["md_score"], tgfb=ams_l["tgfb_score"],
        bins=ams_l_meta["bins"], gene_mean=ams_l_meta["gene_mean"],
        cluster=obs_l["cluster"].to_numpy(int),
        day=obs_l["day"].to_numpy(int),
        n_genes=obs_l["n_genes"].to_numpy(float) if "n_genes" in obs_l.columns else np.full(len(obs_l), np.nan),
    )
    # also state scores
    state_mat = {k: ams_l[k] for k in ams_l if k.startswith("state_")}
    np.savez_compressed(MD2_PROC / "louvain_state_GM00731.npz", **state_mat)

    cluster_df = pd.concat(tables, axis=0, ignore_index=True)
    cluster_df.to_csv(MD2_DIR / "t1_cluster_table.csv", index=False)
    null_df = pd.DataFrame(null_rows)
    null_df.to_csv(MD2_DIR / "t1_age_null.csv", index=False)
    reading = _reading(null_df)
    dump_json(MD2_DIR / "t1_reading.json", jsonable(reading))
    log(f"[t1 reading] key={reading['key']}")
    log(reading["text"])
    summary = dict(
        reading=reading,
        n_cluster_rows=int(len(cluster_df)),
        n_null_rows=int(len(null_df)),
        n_null_ok=int(null_df.ok.sum()) if "ok" in null_df.columns else 0,
        min_cells=MIN_CELLS,
        resolutions=["k3", "k8", "k15", "louvain"],
        k3_age_reused_from_fibro3=True,
        ams_algorithm=ams_meta.get("algorithm"),
        nbin=ams_meta.get("nbin"), ctrl=ams_meta.get("ctrl"),
        seed=MD2_SEED,
    )
    dump_json(MD2_DIR / "t1_summary.json", jsonable(summary))
    man = load_manifest()
    man["status"] = "TASK1_DONE"
    man["t1_reading"] = reading.get("key")
    save_manifest(man)
    progress_snapshot("Task 2 Louvain instrument comparison", stop="TASK1_DONE")
    if close_log:
        log.close()
    return summary


if __name__ == "__main__":
    try:
        run_task1()
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise

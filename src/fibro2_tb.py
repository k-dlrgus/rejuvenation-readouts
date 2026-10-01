"""T-B — QC audit of GSE297234 and d0→d7 random-direction null."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fibro2_common import (  # noqa: E402
    FIBRO2_DIR, FIBRO2_SEED, N_RANDOM_DIR, PREREG_TB, PREREG_TB_FLAG,
    StopStep, Logger, dump_json, jsonable, fibro2_log_banner,
    load_manifest, save_manifest, record_failure, log_columns,
    load_frozen_ruler, progress_snapshot, score_frozen, tmm_logcpm_rows,
)
from fibro_common import FIBRO_DIR  # noqa: E402
from fibro_stage2 import (  # noqa: E402
    find_10x_h5, load_one_h5, attach_meta, align_counts_to_ruler,
    random_directions, AGED_LINE, YOUNG_LINE, CLUSTER_K, CLUSTER_NPC, CLUSTER_MIN_CELLS,
    _log1p_cp10k, _row_sum, FIBRO_RAW,
)
from trajectory_common import permutation_p  # noqa: E402
from fibro_common import unit  # noqa: E402
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA


def mito_mask(var: pd.DataFrame):
    """MT- symbol prefix on GRCh38. Log the rule; do not guess other prefixes."""
    if "symbol" not in var.columns:
        raise StopStep("qc", f"var missing symbol. columns={list(var.columns)}")
    sym = var["symbol"].astype(str)
    mask = sym.str.upper().str.startswith("MT-")
    return mask.to_numpy(), dict(
        rule="symbol upper startswith MT-",
        n_mito=int(mask.sum()),
        mito_symbols=sorted(sym[mask].unique().tolist())[:40],
        columns_used=["symbol"],
    )


def qc_one(lib, cluster_sizes, log):
    X = lib["X"].tocsr()
    n_cells = int(X.shape[0])
    umi = np.asarray(X.sum(axis=1)).ravel()
    ngenes = np.asarray((X > 0).sum(axis=1)).ravel()
    mito, mito_meta = mito_mask(lib["var"])
    if mito.any():
        mito_umi = np.asarray(X[:, mito].sum(axis=1)).ravel()
        frac_cell = mito_umi / np.where(umi > 0, umi, np.nan)
        lib_frac = float(mito_umi.sum() / umi.sum()) if umi.sum() > 0 else np.nan
        med_frac = float(np.nanmedian(frac_cell))
    else:
        lib_frac = np.nan
        med_frac = np.nan
        log(f"[qc] {lib['cell_line']} d{lib['day']}: zero MT- genes under rule {mito_meta['rule']}")
    rec = dict(
        cell_line=lib["cell_line"], day=int(lib["day"]), gsm=str(lib["obs"]["gsm"].iloc[0]),
        age_years=int(lib["age_years"]), n_cells=n_cells,
        median_UMI_per_cell=float(np.median(umi)) if n_cells else np.nan,
        median_genes_per_cell=float(np.median(ngenes)) if n_cells else np.nan,
        mito_read_fraction_library=lib_frac,
        mito_read_fraction_median_cell=med_frac,
        n_mito_genes=mito_meta["n_mito"],
        mito_rule=mito_meta["rule"],
        chemistry=lib.get("chemistry"), software=lib.get("software"),
    )
    for k in range(CLUSTER_K):
        rec[f"n_cluster_{k}"] = int(cluster_sizes[k]) if cluster_sizes is not None and k < len(cluster_sizes) else np.nan
    log(f"[qc] {rec['cell_line']} d{rec['day']} n_cells={n_cells} "
        f"medUMI={rec['median_UMI_per_cell']:.0f} medGenes={rec['median_genes_per_cell']:.0f} "
        f"mito_lib={lib_frac if np.isfinite(lib_frac) else 'NA'} clusters={cluster_sizes}")
    log_columns(f"qc_{rec['gsm']}", ["barcode", "gene_id", "symbol", "feature_type", "genome"],
                lib["obs"]["file"].iloc[0] if "file" in lib["obs"].columns else "h5")
    return rec, mito_meta


def cluster_sizes_one(lib, frozen, log):
    Y, rec, idx = align_counts_to_ruler(lib["X"], lib["var"], frozen, log)
    n = int(Y.shape[0])
    if n < CLUSTER_MIN_CELLS:
        return None, dict(skipped=True, n_cells=n, reason=f"n_cells<{CLUSTER_MIN_CELLS}")
    colsum = _row_sum(Y)
    keep = colsum > 0
    Ykeep = Y[:, keep]
    if sparse.issparse(Ykeep):
        Ykeep = Ykeep.astype(np.float32).toarray()
    Z = _log1p_cp10k(np.asarray(Ykeep, np.float32))
    npc = int(min(CLUSTER_NPC, n - 1, Z.shape[1]))
    pca = PCA(n_components=npc, random_state=int(FIBRO2_SEED))
    P = pca.fit_transform(Z)
    km = KMeans(n_clusters=CLUSTER_K, random_state=int(FIBRO2_SEED), n_init=10)
    lab = km.fit_predict(P)
    sizes = np.bincount(lab, minlength=CLUSTER_K).tolist()
    return sizes, dict(skipped=False, labels=lab, n_cells=n)


def decline_p(obs, scores, null_scores, line, d_end):
    m0 = (obs.cell_line.astype(str).str.upper() == line) & (obs.day.astype(int) == 0)
    m1 = (obs.cell_line.astype(str).str.upper() == line) & (obs.day.astype(int) == int(d_end))
    if int(m0.sum()) != 1 or int(m1.sum()) != 1:
        return dict(ok=False, n0=int(m0.sum()), n1=int(m1.sum()),
                    reason=f"{line} d0 n={int(m0.sum())} d{d_end} n={int(m1.sum())} (want 1)")
    i0 = int(np.flatnonzero(m0.to_numpy())[0])
    i1 = int(np.flatnonzero(m1.to_numpy())[0])
    real = float(scores[i0] - scores[i1])
    null_d = np.asarray(null_scores)[i0] - np.asarray(null_scores)[i1]
    p = permutation_p(real, null_d, greater=True)
    n_ge = int(np.sum(null_d >= real))
    return dict(
        ok=True, cell_line=line, endpoint=f"d0→d{d_end}", d_end=int(d_end),
        score_d0=float(scores[i0]), score_end=float(scores[i1]),
        decline=real, p=p, n_random_ge_real=n_ge, n_random=int(len(null_d)),
        pass_p=bool(np.isfinite(p) and p <= 0.05),
    )


def cluster_rank_pairs(obs_b, scores_b, null_b, line, d_end):
    """Match clusters by within-timepoint n_cells rank (0=largest)."""
    rows = []
    sub0 = obs_b[(obs_b.cell_line.astype(str).str.upper() == line) & (obs_b.day.astype(int) == 0)].copy()
    sub1 = obs_b[(obs_b.cell_line.astype(str).str.upper() == line) & (obs_b.day.astype(int) == int(d_end))].copy()
    if sub0.empty or sub1.empty:
        return [dict(ok=False, cell_line=line, endpoint=f"d0→d{d_end}",
                     reason="missing cluster rows at an endpoint")]
    sub0 = sub0.sort_values("n_cells", ascending=False).reset_index()
    sub1 = sub1.sort_values("n_cells", ascending=False).reset_index()
    n = min(len(sub0), len(sub1))
    for rnk in range(n):
        i0 = int(sub0.loc[rnk, "index"])
        i1 = int(sub1.loc[rnk, "index"])
        # obs_b index may be the dataframe index; scores aligned to obs_b row order
        pos = {idx: j for j, idx in enumerate(obs_b.index)}
        j0, j1 = pos[i0], pos[i1]
        real = float(scores_b[j0] - scores_b[j1])
        null_d = np.asarray(null_b)[j0] - np.asarray(null_b)[j1]
        p = permutation_p(real, null_d, greater=True)
        rows.append(dict(
            ok=True, cell_line=line, endpoint=f"d0→d{d_end}", d_end=int(d_end),
            variant="cluster_size_rank", rank=int(rnk),
            n_cells_d0=int(sub0.loc[rnk, "n_cells"]),
            n_cells_end=int(sub1.loc[rnk, "n_cells"]),
            cluster_id_d0=int(sub0.loc[rnk, "cluster"]) if "cluster" in sub0.columns else np.nan,
            cluster_id_end=int(sub1.loc[rnk, "cluster"]) if "cluster" in sub1.columns else np.nan,
            score_d0=float(scores_b[j0]), score_end=float(scores_b[j1]),
            decline=real, p=p, n_random_ge_real=int(np.sum(null_d >= real)),
            n_random=int(len(null_d)),
            pass_p=bool(np.isfinite(p) and p <= 0.05),
        ))
    return rows


def d10_differs_materially(qc: pd.DataFrame, log):
    """Compare d10 vs d0/d3/d7 on cells recovered, median UMI, mito fraction. Per donor."""
    notes = []
    flags = []
    for line, g in qc.groupby("cell_line"):
        g = g.sort_values("day")
        pre = g[g.day.isin([0, 3, 7])]
        d10 = g[g.day == 10]
        if pre.empty or d10.empty:
            notes.append(f"{line}: missing pre or d10 row")
            continue
        for col, label in (
            ("n_cells", "cells recovered"),
            ("median_UMI_per_cell", "median UMI/cell"),
            ("mito_read_fraction_library", "mitochondrial read fraction"),
        ):
            pre_med = float(np.nanmedian(pre[col].to_numpy(float)))
            v10 = float(d10[col].iloc[0])
            if not np.isfinite(pre_med) or not np.isfinite(v10) or pre_med == 0:
                notes.append(f"{line} {label}: pre_med={pre_med} d10={v10} (cannot ratio)")
                continue
            ratio = v10 / pre_med
            # material: |log2 ratio| ≥ 1 (2-fold) or mito absolute delta ≥ 0.05
            material = abs(np.log2(max(ratio, 1e-12))) >= 1.0
            if col.startswith("mito") and abs(v10 - pre_med) >= 0.05:
                material = True
            flags.append(dict(
                cell_line=line, metric=col, label=label,
                pre_median=pre_med, d10=v10, ratio=ratio, material=bool(material),
            ))
            notes.append(
                f"{line} {label}: d0/d3/d7 median={pre_med:.4g} d10={v10:.4g} "
                f"ratio={ratio:.3f} material={material}"
            )
    any_mat = any(f["material"] for f in flags)
    log("[qc] d10 vs d0/d3/d7: " + " | ".join(notes))
    return dict(any_material=bool(any_mat), flags=flags, notes=notes)


def tb_reading(null_df, qc_diff, log):
    """All-cell primary, aged donor. Young and cluster are reported, not averaged in."""
    aged_a = null_df[
        (null_df.variant == "all_cell") & (null_df.cell_line == AGED_LINE) & (null_df.ok == True)
    ]
    r7 = aged_a[aged_a.d_end == 7]
    r10 = aged_a[aged_a.d_end == 10]
    if r7.empty or r10.empty:
        text = "aged-donor all-cell d0→d7 or d0→d10 row missing. Cannot fire a T-B bullet."
        return dict(key="missing", text=text, d10_technical=None)
    p7 = float(r7.iloc[0]["p"])
    p10 = float(r10.iloc[0]["p"])
    d10_tech = bool(qc_diff.get("any_material"))
    if p7 > 0.05:
        text = (
            "d0→d7 also gives p > 0.05 → the `no_decline` verdict does not depend on d10; "
            "d10's behaviour is a side observation."
        )
        return dict(key="no_decline_independent_of_d10", text=text, d10_technical=d10_tech,
                    p7=p7, p10=p10)
    if p7 <= 0.05 and not (p10 <= 0.05) and d10_tech:
        text = (
            "d0→d7 gives p ≤ 0.05 while d0→d10 does not, and the QC table shows d10 "
            "differing materially from d0/d3/d7 on cells recovered, depth, or mitochondrial "
            "fraction → the endpoint choice is confounded by a technical difference at d10. "
            "Caveat on the verdict, not a reversal of it: a post-hoc endpoint that was not "
            "pre-registered cannot become the primary result. What would settle it: a cohort "
            "with more donors and a pre-registered endpoint."
        )
        return dict(key="d10_technical_caveat", text=text, d10_technical=True, p7=p7, p10=p10)
    if p7 <= 0.05 and not (p10 <= 0.05) and not d10_tech:
        text = (
            "d0→d7 gives p ≤ 0.05 with no technical difference at d10 → the discrepancy "
            "is unexplained."
        )
        return dict(key="d7_unexplained", text=text, d10_technical=False, p7=p7, p10=p10)
    if p7 <= 0.05 and p10 <= 0.05:
        text = (
            f"d0→d7 p={p7:.3f} and d0→d10 p={p10:.3f} both ≤ 0.05. "
            "That is not a listed T-B bullet (the Stage 2 d0→d10 all-cell p was > 0.05). "
            "See the null table. Not averaged into a headline. Not a reversal of Stage 2."
        )
        return dict(key="both_le_005", text=text, d10_technical=d10_tech, p7=p7, p10=p10)
    text = "T-B pattern did not match a pre-registered bullet. See tables."
    return dict(key="unclassified", text=text, d10_technical=d10_tech, p7=p7, p10=p10)


def _local_allcell(libs, frozen, log):
    rows, mats = [], []
    for lib in libs:
        Y, rec, idx = align_counts_to_ruler(lib["X"], lib["var"], frozen, log)
        pb = _row_sum(Y)
        rows.append(dict(
            cell_line=lib["cell_line"], day=lib["day"], age_years=lib["age_years"],
            gsm=lib["obs"]["gsm"].iloc[0], n_cells=int(Y.shape[0]),
            variant="all_cell", cluster=-1,
            age_source=lib["obs"]["age_source"].iloc[0],
        ))
        mats.append(pb)
        log(f"[allcell] {lib['cell_line']} d{lib['day']} overlap={rec.get('n_overlap')}")
    C = np.vstack(mats)
    logcpm, nf = tmm_logcpm_rows(C, log, tag="tb_allcell")
    return pd.DataFrame(rows), logcpm, C


def _local_cluster(libs, frozen, log):
    rows, mats, skipped = [], [], []
    for lib in libs:
        Y, rec, idx = align_counts_to_ruler(lib["X"], lib["var"], frozen, log)
        n = int(Y.shape[0])
        if n < CLUSTER_MIN_CELLS:
            skipped.append(dict(cell_line=lib["cell_line"], day=lib["day"], n_cells=n,
                                reason=f"n_cells<{CLUSTER_MIN_CELLS}"))
            continue
        colsum = _row_sum(Y)
        keep = colsum > 0
        Ykeep = Y[:, keep]
        if sparse.issparse(Ykeep):
            Ykeep = Ykeep.astype(np.float32).toarray()
        Z = _log1p_cp10k(np.asarray(Ykeep, np.float32))
        npc = int(min(CLUSTER_NPC, n - 1, Z.shape[1]))
        pca = PCA(n_components=npc, random_state=int(FIBRO2_SEED))
        P = pca.fit_transform(Z)
        km = KMeans(n_clusters=CLUSTER_K, random_state=int(FIBRO2_SEED), n_init=10)
        lab = km.fit_predict(P)
        for k in range(CLUSTER_K):
            m = lab == k
            if int(m.sum()) == 0:
                continue
            pb = _row_sum(Y, m)
            rows.append(dict(
                cell_line=lib["cell_line"], day=lib["day"], age_years=lib["age_years"],
                gsm=lib["obs"]["gsm"].iloc[0], n_cells=int(m.sum()),
                variant="cluster", cluster=int(k),
                age_source=lib["obs"]["age_source"].iloc[0],
            ))
            mats.append(pb)
        log(f"[cluster] {lib['cell_line']} d{lib['day']} sizes={np.bincount(lab, minlength=CLUSTER_K).tolist()}")
    if not mats:
        return pd.DataFrame(), np.zeros((0, int(np.asarray(frozen['w']).shape[0]))), skipped, np.zeros((0, 1))
    C = np.vstack(mats)
    logcpm, nf = tmm_logcpm_rows(C, log, tag="tb_cluster")
    return pd.DataFrame(rows), logcpm, skipped, C


def run_tb(log=None):
    if not PREREG_TB_FLAG.exists():
        raise StopStep("prereg", "PREREG_TB.flag missing")
    close_log = False
    if log is None:
        log = Logger(FIBRO2_DIR / "tb_report.txt")
        close_log = True
    fibro2_log_banner(log, "T-B")
    log(PREREG_TB)
    frozen = load_frozen_ruler()
    extract = FIBRO_RAW / "RAW"
    if not extract.exists():
        raise StopStep("tb", f"missing GSE297234 extract {extract}")
    samp_path = FIBRO_DIR / "stage2_geo_samples.csv"
    if not samp_path.exists():
        raise StopStep("tb", f"missing {samp_path} — Stage 2 sample table is required and is not re-fetched here")
    samples = pd.read_csv(samp_path)
    log(f"[tb] loaded Stage 2 sample table {samp_path} n={len(samples)} "
        f"columns actually read: {list(samples.columns)}")
    log_columns("tb_stage2_geo_samples", list(samples.columns), str(samp_path))
    need = {"gsm", "cell_line", "day", "age_years", "age_source"}
    missing_cols = sorted(need - set(samples.columns))
    if missing_cols:
        raise StopStep("tb", f"{samp_path} missing columns {missing_cols}. Not substituting.")
    h5s = find_10x_h5(extract, log)
    if not h5s:
        raise StopStep("tb", "no filtered_feature_bc_matrix.h5 under GSE297234 RAW")
    libs = [load_one_h5(h, log) for h in h5s]
    libs = attach_meta(libs, samples, log)
    log(f"[tb] loaded {len(libs)} libraries")
    verify = dict(source=str(samp_path), n_samples=int(len(samples)),
                  cell_lines=sorted(samples.cell_line.astype(str).str.upper().unique().tolist()),
                  days=sorted(samples.day.dropna().astype(int).unique().tolist()))

    qc_rows = []
    mito_rules = []
    size_map = {}
    for lib in libs:
        sizes, meta = cluster_sizes_one(lib, frozen, log)
        size_map[(lib["cell_line"], int(lib["day"]))] = sizes
        rec, mito_meta = qc_one(lib, sizes, log)
        qc_rows.append(rec)
        mito_rules.append(mito_meta)
    qc = pd.DataFrame(qc_rows).sort_values(["cell_line", "day"])
    qc.to_csv(FIBRO2_DIR / "tb_qc.csv", index=False)
    dump_json(FIBRO2_DIR / "tb_mito_rule.json", jsonable(mito_rules[0] if mito_rules else {}))

    # scores + random directions (same spec as Stage 2; writes only under results/fibro2/)
    obs_a, logcpm_a, counts_a = _local_allcell(libs, frozen, log)
    age_a, Z_a, missing_a = score_frozen(logcpm_a, frozen, counts=counts_a)
    null_a = random_directions(frozen["w"], Z_a, seed=FIBRO2_SEED, n=N_RANDOM_DIR)
    log(f"[tb] all-cell missing z=0 columns={int(missing_a.sum())} n_random={N_RANDOM_DIR}")

    obs_b, logcpm_b, skipped_b, counts_b = _local_cluster(libs, frozen, log)
    age_b = null_b = None
    if len(obs_b):
        age_b, Z_b, missing_b = score_frozen(logcpm_b, frozen, counts=counts_b)
        null_b = random_directions(frozen["w"], Z_b, seed=FIBRO2_SEED, n=N_RANDOM_DIR)
        log(f"[tb] cluster missing z=0 columns={int(missing_b.sum())}")
        obs_b = obs_b.copy()
        obs_b["age_score"] = age_b
        obs_b.to_csv(FIBRO2_DIR / "tb_cluster_scores.csv", index=False)

    obs_a = obs_a.copy()
    obs_a["age_score"] = age_a
    obs_a.to_csv(FIBRO2_DIR / "tb_allcell_scores.csv", index=False)
    np.savez_compressed(
        FIBRO2_DIR / "tb_allcell_z.npz",
        Z=np.asarray(Z_a, np.float32), age=np.asarray(age_a, np.float64),
        null_age=null_a, day=obs_a.day.to_numpy(),
        cell_line=obs_a.cell_line.astype(str).to_numpy(),
        variant=obs_a.variant.astype(str).to_numpy() if "variant" in obs_a.columns else np.array(["all_cell"] * len(obs_a)),
        n_cells=obs_a.n_cells.to_numpy(),
    )
    if age_b is not None:
        np.savez_compressed(
            FIBRO2_DIR / "tb_cluster_z.npz",
            Z=np.asarray(Z_b, np.float32), age=np.asarray(age_b, np.float64),
            null_age=null_b, day=obs_b.day.to_numpy(),
            cell_line=obs_b.cell_line.astype(str).to_numpy(),
            cluster=obs_b.cluster.to_numpy(),
            n_cells=obs_b.n_cells.to_numpy(),
        )

    null_rows = []
    for line in (AGED_LINE, YOUNG_LINE):
        for dend in (10, 7):
            rec = decline_p(obs_a, age_a, null_a, line, dend)
            rec["variant"] = "all_cell"
            rec["primary"] = bool(line == AGED_LINE)
            null_rows.append(rec)
            if rec.get("ok"):
                log(f"[null] all_cell {line} d0→d{dend} decline={rec['decline']:+.4f} "
                    f"p={rec['p']:+.3f} n_ge={rec['n_random_ge_real']}/{rec['n_random']}")
    if age_b is not None:
        for line in (AGED_LINE, YOUNG_LINE):
            for dend in (10, 7):
                for rec in cluster_rank_pairs(obs_b, age_b, null_b, line, dend):
                    rec["primary"] = False
                    null_rows.append(rec)
                    if rec.get("ok"):
                        log(f"[null] cluster_rank{rec.get('rank')} {line} d0→d{dend} "
                            f"decline={rec['decline']:+.4f} p={rec['p']:+.3f}")
    null_df = pd.DataFrame(null_rows)
    null_df.to_csv(FIBRO2_DIR / "tb_null.csv", index=False)

    qc_diff = d10_differs_materially(qc, log)
    dump_json(FIBRO2_DIR / "tb_d10_qc_diff.json", jsonable(qc_diff))
    reading = tb_reading(null_df, qc_diff, log)
    log(f"[tb reading] key={reading['key']} {reading['text']}")
    summary = dict(
        reading=reading, qc_diff=qc_diff, n_qc=int(len(qc)), n_null=int(len(null_df)),
        seed=FIBRO2_SEED, n_random=N_RANDOM_DIR,
        skipped_clusters=skipped_b if isinstance(skipped_b, list) else [],
        verify=verify,
    )
    dump_json(FIBRO2_DIR / "tb_summary.json", jsonable(summary))
    man = load_manifest()
    man["status"] = "TB_DONE"
    man["tb_reading"] = reading
    save_manifest(man)
    progress_snapshot("T-C HistGradientBoosting vs ridge", stop="TB_DONE")
    if close_log:
        log.close()
    return summary


if __name__ == "__main__":
    try:
        run_tb()
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise

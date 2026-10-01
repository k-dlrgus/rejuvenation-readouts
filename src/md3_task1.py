"""Task 1 — score the partially reprogrammed population at pooled and per-cluster levels."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md3_common import (  # noqa: E402
    MD3_DIR, MD3_PROC, MD2_DIR, MD2_PROC, MD3_SEED, MIN_CELLS_20, MIN_CELLS_10,
    MIN_CELLS_POOLED, AGED_LINE, YOUNG_LINE, N_RANDOM_DIR, AMS_CTRL,
    PREREG_TASK1, PREREG_TASK1_FLAG, DECLARED_BEFORE_SCORES_FLAG,
    STATE_NAMES, GENELIST_NULL_SEEDS, FROZEN_RULER, MD2_NEGATIVE_HOLDS,
    StopStep, Logger, dump_json, jsonable, load_json, md3_log_banner,
    load_manifest, save_manifest, record_failure, log_columns,
    progress_snapshot, write_prereg_flag,
)
from md3_idtype import detect_id_type, require_mappable, csr_from_npz  # noqa: E402
from fibro2_common import load_frozen_ruler  # noqa: E402
from fibro_stage2 import random_directions  # noqa: E402
from fibro_common import PLURI_ENDOGENOUS, FIBRO_IDENTITY, unit  # noqa: E402
from trajectory_common import permutation_p  # noqa: E402
from md2_score import (  # noqa: E402
    lognormalize_csr, add_module_score, score_frozen_pseudobulk, map_features,
    draw_size_matched_indices, ams_from_indices,
)
from fibro_stage2 import gene_index_by_symbol  # noqa: E402


def _require():
    if not PREREG_TASK1_FLAG.exists():
        raise StopStep("prereg", "PREREG_TASK1.flag missing")
    if not FROZEN_RULER.exists():
        raise StopStep("frozen_ruler", f"missing {FROZEN_RULER}")
    if not (MD3_DIR / "genesets.json").exists():
        raise StopStep("task1", "genesets.json missing")
    for p in (
        MD2_DIR / "t2_cluster_labels_GM00731.csv",
        MD2_DIR / "t2_cluster_labels_GM23815.csv",
        MD2_DIR / "louvain_obs_GM00731.csv",
        MD2_DIR / "louvain_obs_GM23815.csv",
        MD2_PROC / "louvain_counts_GM00731.npz",
        MD2_PROC / "louvain_counts_GM23815.npz",
        MD2_PROC / "rulerY_GM00731.npz",
        MD2_PROC / "rulerY_GM23815.npz",
        MD2_DIR / "allcell_obs_GM00731.csv",
        MD2_DIR / "allcell_obs_GM23815.csv",
    ):
        if not p.exists():
            raise StopStep("task1", f"missing md2 artifact {p}. Not substituting.")


def _load_ruler_y(path):
    z = np.load(path)
    return sparse.csr_matrix(
        (z["data"], z["indices"], z["indptr"]),
        shape=tuple(int(x) for x in z["shape"]),
    )


def _align_louvain_y(obs_l, obs_all, Yall, line):
    key_all = obs_all["gsm"].astype(str) + "||" + obs_all["barcode"].astype(str)
    pos = {k: i for i, k in enumerate(key_all)}
    key = obs_l["gsm"].astype(str) + "||" + obs_l["barcode"].astype(str)
    miss = [k for k in key if k not in pos]
    if miss:
        raise StopStep("task1", f"{line} louvain {len(miss)} cells not in allcell obs")
    take = np.array([pos[k] for k in key], dtype=int)
    if Yall.shape[0] != len(obs_all):
        raise StopStep("task1", f"{line} rulerY n={Yall.shape[0]} allcell n={len(obs_all)}")
    return Yall[take]


def occupancy_table(log):
    rows = []
    for line in (AGED_LINE, YOUNG_LINE):
        lab = pd.read_csv(MD2_DIR / f"t2_cluster_labels_{line}.csv")
        log(f"[t1 occupancy] {line} label columns actually read: {list(lab.columns)}")
        log_columns(f"t2_cluster_labels_{line}", list(lab.columns),
                    str(MD2_DIR / f"t2_cluster_labels_{line}.csv"))
        need = {"cluster", "label", "n_cells"}
        missing = sorted(need - set(lab.columns))
        if missing:
            raise StopStep("task1", f"t2_cluster_labels_{line}.csv missing {missing}")
        obs = pd.read_csv(MD2_DIR / f"louvain_obs_{line}.csv")
        log_columns(f"louvain_obs_{line}", list(obs.columns), str(MD2_DIR / f"louvain_obs_{line}.csv"))
        if "cluster" not in obs.columns or "day" not in obs.columns:
            raise StopStep("task1", f"louvain_obs_{line} missing cluster/day")
        lab_map = {int(r.cluster): str(r.label) for _, r in lab.iterrows()}
        obs = obs.copy()
        obs["label"] = obs["cluster"].map(lambda c: lab_map.get(int(c)))
        if obs["label"].isna().any():
            raise StopStep("task1", f"{line}: clusters without md2 labels: "
                           f"{sorted(obs.loc[obs.label.isna(), 'cluster'].unique())}")
        for _, r in lab.iterrows():
            cl = int(r.cluster)
            rec = dict(
                cell_line=line, cluster=cl, label=str(r.label),
                n_cells_all=int(r.n_cells),
                assignment_basis=(
                    "argmax of mean AddModuleScore of mmc3 Reprog_cell_state_signatures "
                    "(Fibroblast, PartialReprog, EarlyPluripotency, Pluripotency, NonReprog); "
                    "computed in md2, reused, not re-labelled"
                ),
            )
            for st in STATE_NAMES:
                col = f"mean_{st}"
                rec[col] = float(r[col]) if col in lab.columns else np.nan
            for d in (0, 3, 7, 10):
                rec[f"n_d{d}"] = int(((obs.cluster.astype(int) == cl) & (obs.day.astype(int) == d)).sum())
            rows.append(rec)
            log(f"[t1 occupancy] {line} c{cl} label={r.label} n0={rec['n_d0']} n3={rec['n_d3']} "
                f"n7={rec['n_d7']} n10={rec['n_d10']}")
    df = pd.DataFrame(rows)
    df.to_csv(MD3_DIR / "t1_occupancy.csv", index=False)
    pooled = []
    for line in (AGED_LINE, YOUNG_LINE):
        sub = df[df.cell_line == line]
        for st in STATE_NAMES:
            rec = dict(cell_line=line, label=st, n_clusters=int((sub.label == st).sum()))
            for d in (0, 3, 7, 10):
                rec[f"n_d{d}"] = int(sub.loc[sub.label == st, f"n_d{d}"].sum())
            pooled.append(rec)
            log(f"[t1 occupancy pooled] {line} {st} n0={rec['n_d0']} n3={rec['n_d3']} "
                f"n7={rec['n_d7']} n10={rec['n_d10']}")
    pdf = pd.DataFrame(pooled)
    pdf.to_csv(MD3_DIR / "t1_occupancy_pooled.csv", index=False)
    return df, pdf


def _pb_table(obs, Y, frozen, log, tag, ams, by, min_cells):
    """One row per (day, group). TMM+ruler on the same row order as AMS means."""
    rows, mats = [], []
    for (day, grp), g in obs.groupby(["day", by], sort=True):
        idx = g.index.to_numpy()
        n = int(len(idx))
        n_tp = int((obs.day == day).sum())
        rec = dict(
            day=int(day), n_cells=n, n_timepoint=n_tp,
            frac_timepoint=float(n / n_tp) if n_tp else np.nan,
            below_min_cells=bool(n < min_cells),
        )
        if by == "cluster":
            rec["cluster"] = int(grp)
            rec["label"] = str(g["label"].iloc[0]) if "label" in g.columns else None
        else:
            rec["label"] = str(grp)
        for k, arr in ams.items():
            rec[k] = float(np.mean(arr[idx])) if n else np.nan
        rows.append(rec)
        mats.append(np.asarray(Y[idx].sum(axis=0), dtype=np.float64).ravel())
    df = pd.DataFrame(rows)
    if not mats:
        return df, None
    C = np.vstack(mats)
    scored = score_frozen_pseudobulk(C, frozen, log, tag=tag)
    df["age_score"] = scored["age"]
    df["pluri_with_OSKM"] = scored["pluri_with"]
    df["pluri_without_OSKM"] = scored["pluri_without"]
    df["pluri_primary"] = scored["pluri_with"]
    df["age_source"] = "TMM_log2CPM_frozen_ruler"
    hide = list(ams) + ["age_score", "pluri_with_OSKM", "pluri_without_OSKM", "pluri_primary"]
    for col in hide:
        if col in df.columns:
            df.loc[df.below_min_cells, col] = np.nan
    return df, scored


def _ruler_null_rows(df, scored, key_col, line, aggregation, min_cells, log):
    rows = []
    if scored is None or df is None or not len(df):
        return rows
    Z = scored["Z"]
    age = np.asarray(scored["age"], float)
    if len(df) != len(Z):
        raise StopStep("task1", f"{line} {aggregation}: table n={len(df)} Z n={len(Z)}")
    frozen = load_frozen_ruler()
    null = random_directions(frozen["w"], Z, seed=MD3_SEED, n=N_RANDOM_DIR)
    meta = df.reset_index(drop=True)

    def pick(grp, day):
        m = (meta[key_col].astype(str) == str(grp)) & (meta.day.astype(int) == int(day))
        idx = np.flatnonzero(m.to_numpy())
        if len(idx) == 0:
            return None
        if len(idx) != 1:
            raise StopStep("task1", f"{line} {aggregation} {key_col}={grp} day={day} n_rows={len(idx)}")
        return int(idx[0])

    groups = sorted(meta[key_col].astype(str).unique().tolist(), key=lambda x: (len(x), x))
    for grp in groups:
        for dend in (7, 10):
            i0 = pick(grp, 0)
            i1 = pick(grp, dend)
            rec = dict(
                aggregation=aggregation, min_cells=int(min_cells), cell_line=line,
                group=str(grp), endpoint=f"d0→d{dend}", d_end=int(dend),
                n_random=N_RANDOM_DIR, seed=MD3_SEED,
                instrument="frozen_ruler", null_kind="permuted_ruler_weights",
            )
            if key_col == "cluster":
                rec["cluster"] = int(grp) if str(grp).lstrip("-").isdigit() else grp
                lab = meta.loc[meta[key_col].astype(str) == str(grp), "label"]
                rec["label"] = str(lab.iloc[0]) if len(lab) else None
            else:
                rec["label"] = str(grp)
            if i0 is None or i1 is None:
                rec.update(ok=False, reason="missing group row at an endpoint")
                rows.append(rec)
                continue
            n0 = int(meta.n_cells.iloc[i0])
            n1 = int(meta.n_cells.iloc[i1])
            rec["n_cells_d0"] = n0
            rec["n_cells_end"] = n1
            if n0 < min_cells or n1 < min_cells:
                rec.update(ok=False, reason=f"n_cells<{min_cells} at an endpoint")
                rows.append(rec)
                continue
            if bool(meta.below_min_cells.iloc[i0]) or bool(meta.below_min_cells.iloc[i1]):
                rec.update(ok=False, reason="below_min_cells at an endpoint")
                rows.append(rec)
                continue
            a0, a1 = float(age[i0]), float(age[i1])
            if not np.isfinite(a0) or not np.isfinite(a1):
                rec.update(ok=False, reason="age_score not finite at an endpoint")
                rows.append(rec)
                continue
            real = float(a0 - a1)
            null_d = np.asarray(null)[i0] - np.asarray(null)[i1]
            p = permutation_p(real, null_d, greater=True)
            n_ge = int(np.sum(null_d >= real))
            rec.update(
                ok=True, score_d0=a0, score_end=a1, decline=real, p=p,
                n_random_ge_real=n_ge, pass_p=bool(np.isfinite(p) and p <= 0.05),
            )
            rows.append(rec)
            log(f"[t1 ruler-null {aggregation} min={min_cells}] {line} {grp} d0→d{dend} "
                f"decline={real:+.4f} p={p:+.3f} n_ge={n_ge}/{N_RANDOM_DIR} n0={n0} n1={n1}")
    return rows


def _pluri_null_rows(df, scored, key_col, line, aggregation, min_cells, log):
    """Size-matched random gene sets in frozen-ruler Z space (not AMS)."""
    rows = []
    if scored is None or df is None or not len(df):
        return rows
    Z = np.asarray(scored["Z"], float)
    pluri = np.asarray(scored["pluri_with"], float)
    frozen = load_frozen_ruler()
    pos, _, _ = gene_index_by_symbol(frozen)
    p_idx = [pos[g] for g in PLURI_ENDOGENOUS if str(g).upper() in pos]
    f_idx = [pos[g] for g in FIBRO_IDENTITY if str(g).upper() in pos]
    missing_z = np.asarray(scored["missing"], bool) if "missing" in scored else np.zeros(Z.shape[1], bool)
    pool = np.flatnonzero(~missing_z)
    rec_base = dict(
        aggregation=aggregation, min_cells=int(min_cells), cell_line=line,
        n_random=N_RANDOM_DIR, seed=GENELIST_NULL_SEEDS["pluri_primary"],
        instrument="pluri_primary",
        null_kind="size_matched_random_gene_sets_frozen_ruler_Z_not_AMS",
        n_pluri_mapped=len(p_idx), n_fibro_mapped=len(f_idx),
    )
    if len(p_idx) == 0 or len(f_idx) == 0:
        rec_base.update(ok=False, reason="pluri or fibro list mapped 0 genes on frozen ruler")
        rows.append(rec_base)
        log(f"[t1 pluri-null] {line} STOP this instrument: 0 mapped genes")
        return rows
    if pool.size < (len(p_idx) + len(f_idx) + 1):
        rec_base.update(ok=False, reason=f"frozen-ruler non-missing pool {pool.size} too small")
        rows.append(rec_base)
        return rows
    rng = np.random.default_rng(GENELIST_NULL_SEEDS["pluri_primary"])
    n_draw = N_RANDOM_DIR
    null_scores = np.zeros((Z.shape[0], n_draw), dtype=np.float64)
    for d in range(n_draw):
        take = rng.choice(pool, size=len(p_idx) + len(f_idx), replace=False)
        pi, fi = take[:len(p_idx)], take[len(p_idx):]
        null_scores[:, d] = Z[:, pi].mean(1) - Z[:, fi].mean(1)
    meta = df.reset_index(drop=True)

    def pick(grp, day):
        m = (meta[key_col].astype(str) == str(grp)) & (meta.day.astype(int) == int(day))
        idx = np.flatnonzero(m.to_numpy())
        if len(idx) == 0:
            return None
        if len(idx) != 1:
            raise StopStep("task1", f"pluri-null {line} {grp} day={day} n={len(idx)}")
        return int(idx[0])

    groups = sorted(meta[key_col].astype(str).unique().tolist(), key=lambda x: (len(x), x))
    for grp in groups:
        for dend in (7, 10):
            rec = dict(rec_base)
            rec.update(group=str(grp), endpoint=f"d0→d{dend}", d_end=int(dend))
            if key_col == "cluster":
                rec["cluster"] = int(grp) if str(grp).lstrip("-").isdigit() else grp
                lab = meta.loc[meta[key_col].astype(str) == str(grp), "label"]
                rec["label"] = str(lab.iloc[0]) if len(lab) else None
            else:
                rec["label"] = str(grp)
            i0, i1 = pick(grp, 0), pick(grp, dend)
            if i0 is None or i1 is None:
                rec.update(ok=False, reason="missing group row at an endpoint")
                rows.append(rec)
                continue
            n0, n1 = int(meta.n_cells.iloc[i0]), int(meta.n_cells.iloc[i1])
            rec["n_cells_d0"], rec["n_cells_end"] = n0, n1
            if n0 < min_cells or n1 < min_cells or bool(meta.below_min_cells.iloc[i0]) or bool(meta.below_min_cells.iloc[i1]):
                rec.update(ok=False, reason=f"n_cells<{min_cells} at an endpoint")
                rows.append(rec)
                continue
            s0, s1 = float(pluri[i0]), float(pluri[i1])
            if not np.isfinite(s0) or not np.isfinite(s1):
                rec.update(ok=False, reason="pluri not finite")
                rows.append(rec)
                continue
            real = float(s0 - s1)
            null_d = null_scores[i0] - null_scores[i1]
            p = permutation_p(real, null_d, greater=True)
            n_ge = int(np.sum(null_d >= real))
            rec.update(
                ok=True, score_d0=s0, score_end=s1, decline=real, p=p,
                n_random_ge_real=n_ge, pass_p=bool(np.isfinite(p) and p <= 0.05),
            )
            rows.append(rec)
    return rows


def _geneset_null_pack(logX, symbols, genes, bins, name, seed, line, log):
    feat_idx, mapped, missing = map_features(symbols, genes, log, f"{line}:{name}_null")
    rec_map = dict(
        n_requested=int(len(list(genes))), n_mapped=int(feat_idx.size),
        n_missing=int(len(missing)), missing=list(missing)[:40],
    )
    if feat_idx.size == 0:
        log(f"[t1 AMS-null] {line} {name}: 0 genes mapped of {len(list(genes))}. Not scoring this list.")
        record_failure(
            "addmodulescore",
            f"{line}:{name}: 0 feature genes mapped of {len(list(genes))}. Not scoring.",
            details=rec_map,
        )
        return None, rec_map
    rng = np.random.default_rng(int(seed))
    draws = draw_size_matched_indices(feat_idx, bins, rng, N_RANDOM_DIR)
    return dict(draws=draws, feat_idx=feat_idx, rec_map=rec_map), rec_map


def _ams_from_indices_matvec(X, feat_idx, bins, rng, ctrl=AMS_CTRL):
    """Same control sampling and AMS formula as md2_score.ams_from_indices.

    Means are CSR/CSC matvecs instead of column-submatrix copies. Score is
    mean(features) − mean(unique controls) per cell.
    """
    feat_idx = np.asarray(feat_idx, dtype=int)
    ctrl_pool = []
    for j in feat_idx:
        b = int(bins[j])
        pool = np.flatnonzero(bins == b)
        if pool.size < int(ctrl):
            raise StopStep("addmodulescore", f"bin {b} has {pool.size} < ctrl={ctrl}")
        ctrl_pool.extend(rng.choice(pool, size=int(ctrl), replace=False).tolist())
    ctrl_idx = np.unique(np.asarray(ctrl_pool, dtype=int))
    n_genes = int(X.shape[1])
    wf = np.zeros(n_genes, dtype=np.float64)
    wf[feat_idx] = 1.0 / float(feat_idx.size)
    wc = np.zeros(n_genes, dtype=np.float64)
    wc[ctrl_idx] = 1.0 / float(ctrl_idx.size)
    return np.asarray(X.dot(wf) - X.dot(wc), dtype=np.float64).ravel()


def _group_declines_multi(obs, cell, draws, bins, logX, seed, line, name, n_genes, log):
    """One AMS draw fills pooled + both cluster thresholds. Real scores are the existing per-cell AMS."""
    kind = "size_matched_random_gene_sets_mean_expression_bin"
    specs = (
        ("cluster", "per_cluster", MIN_CELLS_20),
        ("cluster", "per_cluster_sensitivity_ge10", MIN_CELLS_10),
        ("label", "pooled_by_state", MIN_CELLS_POOLED),
    )
    obs = obs.reset_index(drop=True)
    cell = np.asarray(cell, float)
    day = obs.day.to_numpy(int)
    n_draws = int(draws.shape[0])
    spec_state = []
    for key_col, aggregation, min_cells in specs:
        key_arr = obs[key_col].astype(str).to_numpy()
        groups = sorted(pd.unique(key_arr).tolist(), key=lambda x: (len(str(x)), str(x)))
        idx_map, n_map = {}, {}
        for grp in groups:
            for d in (0, 3, 7, 10):
                m = (key_arr == str(grp)) & (day == int(d))
                n = int(m.sum())
                n_map[(str(grp), d)] = n
                idx_map[(str(grp), d)] = np.flatnonzero(m) if n >= min_cells else None
        spec_state.append(dict(
            key_col=key_col, aggregation=aggregation, min_cells=min_cells,
            groups=groups, idx_map=idx_map, n_map=n_map,
            null_mean={(str(g), d): np.full(n_draws, np.nan) for g in groups for d in (0, 3, 7, 10)},
        ))

    def mean_at(scores, st, grp, d):
        idx = st["idx_map"].get((str(grp), int(d)))
        n = st["n_map"].get((str(grp), int(d)), 0)
        if idx is None:
            return np.nan, n
        return float(np.mean(scores[idx])), n

    rng_ams = np.random.default_rng(int(seed) + 17)
    for d in range(n_draws):
        sc = _ams_from_indices_matvec(logX, draws[d], bins, rng_ams, ctrl=AMS_CTRL)
        for st in spec_state:
            for (grp, dy), idx in st["idx_map"].items():
                if idx is None:
                    continue
                st["null_mean"][(grp, dy)][d] = float(np.mean(sc[idx]))
        if d % 20 == 0:
            log(f"[t1 {name}-null] {line} draw {d}/{n_draws}")

    rows = []
    for st in spec_state:
        key_col, aggregation, min_cells = st["key_col"], st["aggregation"], st["min_cells"]
        for grp in st["groups"]:
            for dend in (7, 10):
                m0, n0 = mean_at(cell, st, grp, 0)
                m1, n1 = mean_at(cell, st, grp, dend)
                rec = dict(
                    aggregation=aggregation, min_cells=int(min_cells), cell_line=line,
                    group=str(grp), endpoint=f"d0→d{dend}", d_end=int(dend),
                    n_cells_d0=n0, n_cells_end=n1, n_random=n_draws, seed=int(seed),
                    instrument=name, null_kind=kind, n_genes=int(n_genes),
                )
                if key_col == "cluster":
                    rec["cluster"] = int(grp) if str(grp).lstrip("-").isdigit() else grp
                    lab = obs.loc[obs[key_col].astype(str) == str(grp), "label"]
                    rec["label"] = str(lab.iloc[0]) if len(lab) else None
                else:
                    rec["label"] = str(grp)
                if not np.isfinite(m0) or not np.isfinite(m1):
                    rec.update(ok=False, reason=f"n_cells<{min_cells} at an endpoint")
                    rows.append(rec)
                    continue
                real = float(m0 - m1)
                null_d = st["null_mean"][(str(grp), 0)] - st["null_mean"][(str(grp), int(dend))]
                if not np.isfinite(null_d).all():
                    rec.update(ok=False, reason="null not finite at an endpoint")
                    rows.append(rec)
                    continue
                p = permutation_p(real, null_d, greater=True)
                n_ge = int(np.sum(null_d >= real))
                rec.update(
                    ok=True, score_d0=m0, score_end=m1, decline=real, p=p,
                    n_random_ge_real=n_ge, pass_p=bool(np.isfinite(p) and p <= 0.05),
                )
                rows.append(rec)
                log(f"[t1 {name}-null {aggregation} min={min_cells}] {line} {grp} d0→d{dend} "
                    f"decline={real:+.4f} p={p:+.3f} n_ge={n_ge}/{n_draws} n0={n0} n1={n1}")
    return rows


def _score_one_donor(line, sets, frozen, log):
    Xl, symbols, gene_id, z = csr_from_npz(MD2_PROC / f"louvain_counts_{line}.npz")
    obs = pd.read_csv(MD2_DIR / f"louvain_obs_{line}.csv")
    if len(obs) != Xl.shape[0]:
        raise StopStep("task1", f"{line} louvain obs n={len(obs)} X n={Xl.shape[0]}")
    id_sym = detect_id_type(symbols, log, f"louvain_{line}_symbols",
                            str(MD2_PROC / f"louvain_counts_{line}.npz") + ":symbols")
    require_mappable(id_sym, {"symbol"}, log, "task1")
    if gene_id is not None and len(gene_id):
        detect_id_type(gene_id, log, f"louvain_{line}_gene_id",
                       str(MD2_PROC / f"louvain_counts_{line}.npz") + ":gene_id")
    lab = pd.read_csv(MD2_DIR / f"t2_cluster_labels_{line}.csv")
    lab_map = {int(r.cluster): str(r.label) for _, r in lab.iterrows()}
    obs = obs.reset_index(drop=True)
    obs["label"] = obs["cluster"].map(lambda c: lab_map.get(int(c)))
    if obs["label"].isna().any():
        raise StopStep("task1", f"{line}: unlabeled clusters")

    obs_all = pd.read_csv(MD2_DIR / f"allcell_obs_{line}.csv")
    Yall = _load_ruler_y(MD2_PROC / f"rulerY_{line}.npz")
    Y = _align_louvain_y(obs, obs_all, Yall, line)

    logX = lognormalize_csr(Xl)
    gene_sets_ams = {
        "md_score": list(sets["MD"]),
        "tgfb_score": list(sets["TGFB"]),
    }
    for st, genes in (sets.get("reprog") or {}).items():
        gene_sets_ams[f"state_{st}"] = list(genes)
    ams, ams_meta = add_module_score(logX, symbols, gene_sets_ams, log=log, tag=f"{line}_louvain")
    bins = np.asarray(ams_meta["bins"])
    age_sets = {"age_up": list(sets["age_up"]), "age_down": list(sets["age_down"])}
    try:
        ams_age, ams_age_meta = add_module_score(
            logX, symbols, age_sets, log=log, tag=f"{line}_louvain_S3", bins=bins,
        )
        ams.update(ams_age)
        mapped_s3 = {k: v for k, v in (ams_age_meta.get("per_set") or {}).items()}
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        log(f"[t1] {line} Table S3 AMS STOP [{e.step}] {e.message}")
        mapped_s3 = {}
        ams_age_meta = {}

    dump_json(MD3_DIR / f"t1_ams_meta_{line}.json", jsonable({
        "md_tgfb_states": {
            k: {kk: vv for kk, vv in v.items() if kk != "mapped"}
            for k, v in (ams_meta.get("per_set") or {}).items()
        },
        "table_s3": {
            k: {kk: vv for kk, vv in v.items() if kk != "mapped"}
            for k, v in mapped_s3.items()
        },
        "id_type_symbols": id_sym,
    }))

    ams_report = {k: ams[k] for k in ("md_score", "tgfb_score") if k in ams}
    if "age_up" in ams:
        ams_report["age_up"] = ams["age_up"]
    if "age_down" in ams:
        ams_report["age_down"] = ams["age_down"]

    # TMM among all cluster×day rows (n≥1). Thresholds mask reported scores, matching md2.
    df_c, scored_c = _pb_table(obs, Y, frozen, log, tag=f"t1_cluster_{line}",
                               ams=ams_report, by="cluster", min_cells=1)
    df_c["cell_line"] = line
    df_c["resolution"] = "louvain"
    hide = list(ams_report) + ["age_score", "pluri_with_OSKM", "pluri_without_OSKM", "pluri_primary"]
    df_c20 = df_c.copy()
    df_c20["below_min_cells"] = df_c20.n_cells < MIN_CELLS_20
    df_c20["min_cells"] = MIN_CELLS_20
    df_c20["aggregation"] = "per_cluster"
    for col in hide:
        if col in df_c20.columns:
            df_c20.loc[df_c20.below_min_cells, col] = np.nan
    df_c10 = df_c.copy()
    df_c10["below_min_cells"] = df_c10.n_cells < MIN_CELLS_10
    df_c10["min_cells"] = MIN_CELLS_10
    df_c10["aggregation"] = "per_cluster_sensitivity_ge10"
    for col in hide:
        if col in df_c10.columns:
            df_c10.loc[df_c10.below_min_cells, col] = np.nan

    df_p, scored_p = _pb_table(obs, Y, frozen, log, tag=f"t1_pooled_{line}",
                               ams=ams_report, by="label", min_cells=MIN_CELLS_POOLED)
    df_p["cell_line"] = line
    df_p["resolution"] = "louvain"
    df_p["aggregation"] = "pooled_by_state"
    df_p["min_cells"] = MIN_CELLS_POOLED

    # Ruler/pluri nulls need finite scores at both endpoints; min_cells is applied inside.
    df_c_for_null = df_c.copy()
    df_c_for_null["below_min_cells"] = False
    null_rows = []
    null_rows += _ruler_null_rows(df_c_for_null, scored_c, "cluster", line, "per_cluster", MIN_CELLS_20, log)
    null_rows += _ruler_null_rows(df_c_for_null, scored_c, "cluster", line, "per_cluster_sensitivity_ge10", MIN_CELLS_10, log)
    null_rows += _ruler_null_rows(df_p, scored_p, "label", line, "pooled_by_state", MIN_CELLS_POOLED, log)
    null_rows += _pluri_null_rows(df_c_for_null, scored_c, "cluster", line, "per_cluster", MIN_CELLS_20, log)
    null_rows += _pluri_null_rows(df_c_for_null, scored_c, "cluster", line, "per_cluster_sensitivity_ge10", MIN_CELLS_10, log)
    null_rows += _pluri_null_rows(df_p, scored_p, "label", line, "pooled_by_state", MIN_CELLS_POOLED, log)

    lists = [
        ("md_score", sets["MD"], GENELIST_NULL_SEEDS["md_score"]),
        ("tgfb_score", sets["TGFB"], GENELIST_NULL_SEEDS["tgfb_score"]),
    ]
    if "age_up" in ams:
        lists.append(("age_up", sets["age_up"], GENELIST_NULL_SEEDS["age_up"]))
    if "age_down" in ams:
        lists.append(("age_down", sets["age_down"], GENELIST_NULL_SEEDS["age_down"]))

    for name, genes, seed in lists:
        pack, rec_map = _geneset_null_pack(logX, symbols, genes, bins, name, seed, line, log)
        if pack is None:
            continue
        cell = np.asarray(ams[name], float)
        null_rows += _group_declines_multi(
            obs, cell, pack["draws"], bins, logX, seed, line, name, rec_map["n_mapped"], log,
        )

    return dict(
        cluster20=df_c20, cluster10=df_c10, pooled=df_p,
        null_rows=null_rows, ams_meta=ams_meta,
        n_mapped_s3={k: mapped_s3.get(k, {}).get("n_mapped") for k in ("age_up", "age_down")},
    )


def _reading(null_df):
    def _moves(sub, instrument, label=None):
        g = sub[(sub.instrument == instrument) & (sub.ok == True)].copy()  # noqa: E712
        if label is not None:
            g = g[g.label.astype(str) == str(label)]
        hits = g[g.pass_p == True]  # noqa: E712
        return dict(
            n_ok=int(len(g)), n_pass=int(len(hits)),
            hits=hits[["cell_line", "group", "label", "endpoint", "decline", "p",
                       "n_cells_d0", "n_cells_end"]].to_dict(orient="records") if len(hits) else [],
        )

    aged = null_df[null_df.cell_line.astype(str) == AGED_LINE].copy()
    pooled = aged[aged.aggregation == "pooled_by_state"]
    cl20 = aged[aged.aggregation == "per_cluster"]
    pr_ruler = _moves(pooled, "frozen_ruler", "PartialReprog")
    pr_md = _moves(pooled, "md_score", "PartialReprog")
    cl_ruler = _moves(cl20, "frozen_ruler")
    cl_md = _moves(cl20, "md_score")
    ruler_pooled = pr_ruler["n_pass"] > 0
    md_pooled = pr_md["n_pass"] > 0
    ruler_cl = cl_ruler["n_pass"] > 0
    md_cl = cl_md["n_pass"] > 0

    if ruler_pooled:
        key = "ruler_declines_pooled_PR"
        text = (
            "Ruler declines in the pooled PartialReprog state with p ≤ 0.05. "
            "The FINDINGS_MD2.md `negative_holds` reading was limited by cell counts, not by the ruler. "
            f"Supersedes: \"{MD2_NEGATIVE_HOLDS}\" "
            f"Pooled PartialReprog ruler hits: {pr_ruler['hits']}."
        )
    elif (not ruler_pooled) and md_pooled:
        key = "instruments_disagree_on_claim_population"
        text = (
            "Ruler does not decline in the pooled PartialReprog state while MD does, both with their own nulls. "
            "The instruments disagree on exactly the population the published claim concerns. "
            "No winner is declared. "
            f"Pooled PartialReprog ruler: n_ok={pr_ruler['n_ok']} n_pass={pr_ruler['n_pass']} hits={pr_ruler['hits']}. "
            f"Pooled PartialReprog MD: n_ok={pr_md['n_ok']} n_pass={pr_md['n_pass']} hits={pr_md['hits']}."
        )
    else:
        key = "reproduction_failure_pooled"
        text = (
            "Neither moves in the pooled PartialReprog state. Our reproduction does not recover "
            "their Figure 3G at the pooled level; this is a reproduction failure, not a refutation. "
            f"Pooled PartialReprog ruler n_ok={pr_ruler['n_ok']} n_pass={pr_ruler['n_pass']}; "
            f"MD n_ok={pr_md['n_ok']} n_pass={pr_md['n_pass']}."
        )
    levels_disagree = bool((ruler_pooled != ruler_cl) or (md_pooled != md_cl))
    extra = None
    if levels_disagree:
        extra = (
            "Pooled and per-cluster (≥20) levels disagree. Report both; do not pick. "
            f"pooled PartialReprog ruler_moves={ruler_pooled} md_moves={md_pooled}; "
            f"per-cluster≥20 any-cluster ruler_moves={ruler_cl} md_moves={md_cl}. "
            f"per-cluster ruler hits={cl_ruler['hits']} md hits={cl_md['hits']}."
        )
    return dict(
        key=key, text=text, levels_disagree=levels_disagree, extra=extra,
        pooled_PR_ruler=pr_ruler, pooled_PR_md=pr_md,
        cluster20_ruler=cl_ruler, cluster20_md=cl_md,
        md2_sentence_at_issue=MD2_NEGATIVE_HOLDS,
    )


def run_task1(log=None):
    _require()
    close_log = False
    if log is None:
        log = Logger(MD3_DIR / "t1_report.txt")
        close_log = True
    md3_log_banner(log, "TASK1")
    log(PREREG_TASK1)
    occ, occ_p = occupancy_table(log)
    if not DECLARED_BEFORE_SCORES_FLAG.exists():
        DECLARED_BEFORE_SCORES_FLAG.write_text(
            PREREG_TASK1
            + "\n\nThis flag was written after cluster occupancy (from md2 labels + louvain obs) "
              "and before any MD3 frozen-ruler / MD / TGF-β / Table S3 / pluripotency score.\n"
            + occ_p.to_csv(index=False)
            + "\n",
            encoding="utf-8",
        )
        log(f"[prereg] wrote {DECLARED_BEFORE_SCORES_FLAG}")
    else:
        log(f"[prereg] {DECLARED_BEFORE_SCORES_FLAG.name} already exists; not rewriting")

    sets = load_json(MD3_DIR / "genesets.json")
    frozen = load_frozen_ruler()
    parts20, parts10, parts_p, null_rows = [], [], [], []
    mapped = {}
    for line in (AGED_LINE, YOUNG_LINE):
        rec = _score_one_donor(line, sets, frozen, log)
        parts20.append(rec["cluster20"])
        parts10.append(rec["cluster10"])
        parts_p.append(rec["pooled"])
        null_rows.extend(rec["null_rows"])
        mapped[line] = rec["n_mapped_s3"]

    df20 = pd.concat(parts20, axis=0, ignore_index=True)
    df10 = pd.concat(parts10, axis=0, ignore_index=True)
    dfp = pd.concat(parts_p, axis=0, ignore_index=True)
    df20.to_csv(MD3_DIR / "t1_cluster20_scores.csv", index=False)
    df10.to_csv(MD3_DIR / "t1_cluster10_scores.csv", index=False)
    dfp.to_csv(MD3_DIR / "t1_pooled_scores.csv", index=False)
    null_df = pd.DataFrame(null_rows)
    null_df.to_csv(MD3_DIR / "t1_null.csv", index=False)
    null_df[null_df.aggregation == "pooled_by_state"].to_csv(MD3_DIR / "t1_pooled_null.csv", index=False)
    null_df[null_df.aggregation == "per_cluster"].to_csv(MD3_DIR / "t1_cluster20_null.csv", index=False)
    null_df[null_df.aggregation == "per_cluster_sensitivity_ge10"].to_csv(
        MD3_DIR / "t1_cluster10_null.csv", index=False,
    )

    # Check MD cluster means vs md2 (aged louvain ≥20).
    md2 = pd.read_csv(MD2_DIR / "t2_cluster_table.csv")
    log_columns("md2_t2_cluster_table", list(md2.columns), str(MD2_DIR / "t2_cluster_table.csv"))
    chk_rows = []
    for line in (AGED_LINE, YOUNG_LINE):
        a = df20[(df20.cell_line == line) & (df20.below_min_cells == False)]  # noqa: E712
        b = md2[(md2.cell_line == line) & (md2.resolution.astype(str) == "louvain")
                & (md2.below_min_cells == False)]  # noqa: E712
        m = a.merge(b[["day", "cluster", "md_score"]], on=["day", "cluster"], how="inner", suffixes=("_md3", "_md2"))
        if not len(m):
            continue
        delta = (m["md_score_md3"] - m["md_score_md2"]).abs()
        chk_rows.append(dict(
            cell_line=line, n_shared=int(len(m)),
            max_abs_md_delta=float(delta.max()) if len(delta) else np.nan,
            n_gt_1e_6=int((delta > 1e-6).sum()),
        ))
        log(f"[t1 md-check] {line} n_shared={len(m)} max|ΔMD|={float(delta.max()):.3e}")
    pd.DataFrame(chk_rows).to_csv(MD3_DIR / "t1_md_check_vs_md2.csv", index=False)

    reading = _reading(null_df)
    dump_json(MD3_DIR / "t1_reading.json", jsonable(reading))
    log(f"[t1 reading] key={reading['key']}")
    log(reading["text"])
    if reading.get("extra"):
        log(reading["extra"])
    summary = dict(
        reading=reading, n_cluster20=int(len(df20)), n_cluster10=int(len(df10)),
        n_pooled=int(len(dfp)), n_null=int(len(null_df)),
        min_cells_20=MIN_CELLS_20, min_cells_10=MIN_CELLS_10,
        min_cells_pooled=MIN_CELLS_POOLED, mapped_s3=mapped, seed=MD3_SEED,
    )
    dump_json(MD3_DIR / "t1_summary.json", jsonable(summary))
    man = load_manifest()
    man["status"] = "TASK1_DONE"
    man["t1_reading"] = reading.get("key")
    save_manifest(man)
    progress_snapshot("Task 2 three-way instrument vs donor age + paired Δρ", stop="TASK1_DONE")
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

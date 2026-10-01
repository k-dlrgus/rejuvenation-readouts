"""Task 1 nulls — permuted-ruler-weight directions and size-matched random gene sets."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plane_common import (  # noqa: E402
    PLANE_DIR, MD2_PROC, MD4_PROC, PLANE_SEED, AGED_LINE, YOUNG_LINE,
    STATE_NAMES, TRAJECTORY, SEGMENT_A, IDENTITY_KEY, AGE_INSTRUMENTS,
    MIN_CELLS, N_RANDOM_DIR, AMS_CTRL, AMS_NBIN, AMS_CTRL_STREAM_OFFSET,
    RULER_NULL_CHUNK, GENELIST_NULL_SEEDS, DECLARED_BEFORE_SCORES_FLAG,
    StopStep, dump_json, jsonable, load_frozen_ruler, _slope, segment_B_end,
)
from fibro_common import unit  # noqa: E402
from gtex_common import EDGE_R_PRIOR  # noqa: E402
from md3_idtype import csr_from_npz  # noqa: E402
from md2_score import lognormalize_csr, map_features, draw_size_matched_indices  # noqa: E402
from trajectory_common import permutation_p  # noqa: E402


def _load_Y(line):
    path = MD4_PROC / f"louvain_rulerY_{line}.npz"
    z = np.load(path)
    return sparse.csr_matrix(
        (z["data"], z["indices"], z["indptr"]),
        shape=tuple(int(x) for x in z["shape"]),
    )


def _permuted_weight_matrix(w, seed, n):
    rng = np.random.default_rng(int(seed))
    w = np.asarray(w, float).ravel()
    W = np.zeros((len(w), int(n)), dtype=np.float64)
    for i in range(int(n)):
        wp = rng.permutation(w)
        wp, _ = unit(wp)
        W[:, i] = wp
    return W


def score_random_directions_cells(Y, frozen, log, tag, n_random=N_RANDOM_DIR,
                                  seed=PLANE_SEED, chunk=RULER_NULL_CHUNK):
    """Per-cell log2-CPM (nf=1) then 200 permuted-weight directions. Same Z as md4."""
    Y = Y.tocsr()
    n, p = Y.shape
    mu = np.asarray(frozen["mu"], float)
    sd = np.asarray(frozen["sd"], float)
    sd = np.where(sd < 1e-12, 1.0, sd)
    w = np.asarray(frozen["w"], float)
    w_u, _ = unit(w)
    if len(mu) != p or len(w) != p:
        raise StopStep("nulls", f"{tag}: rulerY p={p} frozen p_mu={len(mu)} p_w={len(w)}")
    lib = np.asarray(Y.sum(axis=1)).ravel().astype(np.float64)
    missing = np.asarray((Y.getnnz(axis=0) == 0))
    finite_lib = lib[np.isfinite(lib) & (lib > 0)]
    ave = float(np.mean(finite_lib)) if finite_lib.size else 1.0
    if ave <= 0:
        ave = 1.0
    W = _permuted_weight_matrix(w, seed, n_random)
    obs = np.full(n, np.nan, dtype=np.float64)
    null = np.full((n, int(n_random)), np.nan, dtype=np.float64)
    log(f"[null ruler {tag}] n={n} p={p} n_missing_allzero={int(missing.sum())} "
        f"n_random={n_random} seed={seed} tmm_nf=1 prior={EDGE_R_PRIOR:g}")
    for start in range(0, n, int(chunk)):
        sl = slice(start, min(start + int(chunk), n))
        dense = np.asarray(Y[sl].todense(), dtype=np.float64)
        lib_s = lib[sl]
        adj_prior = float(EDGE_R_PRIOR) * lib_s / ave
        den = lib_s + 2.0 * adj_prior
        den = np.where(den <= 0, 1.0, den)
        with np.errstate(divide="ignore", invalid="ignore"):
            cpm = (dense + adj_prior[:, None]) / den[:, None] * 1e6
            logcpm = np.log2(np.clip(cpm, 1e-12, None))
        bad = ~(np.isfinite(lib_s) & (lib_s > 0))
        Z = (logcpm - mu) / sd
        Z[:, missing] = 0.0
        age_s = Z @ w_u
        age_s[bad] = np.nan
        obs[sl] = age_s
        ns = Z @ W
        ns[bad, :] = np.nan
        null[sl, :] = ns
        if start == 0 or start % (int(chunk) * 10) == 0:
            log(f"[null ruler {tag}] {start}/{n}")
    return obs, null


def _ams_from_indices_matvec(X, feat_idx, bins, rng, ctrl=AMS_CTRL):
    """Same control sampling and AMS formula as md2_score.ams_from_indices (md3 matvec)."""
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


def _state_medians_1d(labels, scores, states):
    """scores: n_cells or n_cells × n_draw. Returns dict state -> scalar or (n_draw,)"""
    labels = np.asarray(labels).astype(str)
    scores = np.asarray(scores, float)
    out = {}
    nd = None if scores.ndim == 1 else scores.shape[1]
    for st in states:
        m = labels == st
        n = int(m.sum())
        if n < MIN_CELLS:
            out[st] = np.nan if nd is None else np.full(nd, np.nan)
            continue
        block = scores[m]
        if nd is None:
            block = block[np.isfinite(block)]
            out[st] = float(np.median(block)) if len(block) >= MIN_CELLS else np.nan
        else:
            med = np.full(nd, np.nan)
            for d in range(nd):
                col = block[:, d]
                col = col[np.isfinite(col)]
                if len(col) >= MIN_CELLS:
                    med[d] = float(np.median(col))
            out[st] = med
    return out


def _bends_from_state_y(y_by_state, x_by_state, a0, a1, end):
    y0, y1, y2 = y_by_state[a0], y_by_state[a1], y_by_state[end]
    x0, x1, x2 = x_by_state[a0], x_by_state[a1], x_by_state[end]
    y0 = np.asarray(y0, float).ravel()
    y1 = np.asarray(y1, float).ravel()
    y2 = np.asarray(y2, float).ravel()
    n = int(y0.size)
    sA = np.full(n, np.nan)
    sB = np.full(n, np.nan)
    bend = np.full(n, np.nan)
    dxA = float(x1 - x0)
    dxB = float(x2 - x1)
    for i in range(n):
        sA[i] = _slope(float(y1[i] - y0[i]), dxA)
        sB[i] = _slope(float(y2[i] - y1[i]), dxB)
        if np.isfinite(sA[i]) and np.isfinite(sB[i]):
            bend[i] = sA[i] - sB[i]
    return sA, sB, bend


def _scalar_x(x_med, st):
    v = x_med[st]
    if np.ndim(v) != 0:
        return np.nan
    return float(v)


def _gene_list_null_state_y(logX, symbols, bins, genes, seed, labels, log, tag):
    feat_idx, mapped, missing = map_features(symbols, genes, log, tag)
    if feat_idx.size == 0:
        raise StopStep("nulls", f"{tag}: 0 feature genes mapped of {len(list(genes))}. Not scoring.")
    rng_draw = np.random.default_rng(int(seed))
    draws = draw_size_matched_indices(feat_idx, bins, rng_draw, N_RANDOM_DIR)
    rng_ams = np.random.default_rng(int(seed) + AMS_CTRL_STREAM_OFFSET)
    y_by = {st: np.full(N_RANDOM_DIR, np.nan) for st in STATE_NAMES}
    for d in range(N_RANDOM_DIR):
        sc = _ams_from_indices_matvec(logX, draws[d], bins, rng_ams, ctrl=AMS_CTRL)
        meds = _state_medians_1d(labels, sc, STATE_NAMES)
        for st in STATE_NAMES:
            y_by[st][d] = meds[st]
        if d % 20 == 0:
            log(f"[null AMS {tag}] draw {d}/{N_RANDOM_DIR} n_mapped={feat_idx.size}")
    return y_by, dict(
        n_requested=int(len(list(genes))), n_mapped=int(feat_idx.size),
        n_missing=int(len(missing)), missing=list(missing)[:40],
        seed=int(seed), n_random=int(N_RANDOM_DIR), mapped=list(mapped)[:40],
    )


def _pair_null_state_y(logX, symbols, bins, genes_up, seed_up, genes_dn, seed_dn,
                       labels, log, tag):
    feat_up, mapped_up, miss_up = map_features(symbols, genes_up, log, f"{tag}:up")
    feat_dn, mapped_dn, miss_dn = map_features(symbols, genes_dn, log, f"{tag}:down")
    if feat_up.size == 0 or feat_dn.size == 0:
        raise StopStep(
            "nulls",
            f"{tag}: age-up mapped={feat_up.size} age-down mapped={feat_dn.size}. Not scoring.",
        )
    rng_du = np.random.default_rng(int(seed_up))
    rng_dd = np.random.default_rng(int(seed_dn))
    draws_up = draw_size_matched_indices(feat_up, bins, rng_du, N_RANDOM_DIR)
    draws_dn = draw_size_matched_indices(feat_dn, bins, rng_dd, N_RANDOM_DIR)
    rng_au = np.random.default_rng(int(seed_up) + AMS_CTRL_STREAM_OFFSET)
    rng_ad = np.random.default_rng(int(seed_dn) + AMS_CTRL_STREAM_OFFSET)
    y_by = {st: np.full(N_RANDOM_DIR, np.nan) for st in STATE_NAMES}
    for d in range(N_RANDOM_DIR):
        up = _ams_from_indices_matvec(logX, draws_up[d], bins, rng_au, ctrl=AMS_CTRL)
        dn = _ams_from_indices_matvec(logX, draws_dn[d], bins, rng_ad, ctrl=AMS_CTRL)
        sc = up - dn
        meds = _state_medians_1d(labels, sc, STATE_NAMES)
        for st in STATE_NAMES:
            y_by[st][d] = meds[st]
        if d % 20 == 0:
            log(f"[null AMS {tag} up-down] draw {d}/{N_RANDOM_DIR}")
    return y_by, dict(
        age_up=dict(n_requested=int(len(list(genes_up))), n_mapped=int(feat_up.size),
                    n_missing=int(len(miss_up)), missing=list(miss_up)[:40], seed=int(seed_up)),
        age_down=dict(n_requested=int(len(list(genes_dn))), n_mapped=int(feat_dn.size),
                      n_missing=int(len(miss_dn)), missing=list(miss_dn)[:40], seed=int(seed_dn)),
        n_random=int(N_RANDOM_DIR),
    )


def run_nulls(log, obs_by, ams_by, scores_by, sets, qdf, sdf):
    if not DECLARED_BEFORE_SCORES_FLAG.exists():
        raise StopStep("nulls", "DECLARED_BEFORE_SCORES.flag missing. Not scoring nulls.")
    frozen = load_frozen_ruler()
    draw_rows = []
    bend_rows = []
    env_rows = []
    sum_rows = []
    meta = {}

    for line in (AGED_LINE, YOUNG_LINE):
        obs = obs_by[line]
        labels = obs["label"].to_numpy().astype(str)
        x = np.asarray(scores_by[line][IDENTITY_KEY], float)
        x_med = _state_medians_1d(labels, x, STATE_NAMES)
        qrow = qdf[qdf.cell_line == line].iloc[0].to_dict()
        end, fallback = segment_B_end(qrow)
        a0, a1 = SEGMENT_A
        if not (qrow.get("segment_A_ok") and qrow.get("segment_B_ok") and end):
            for inst in AGE_INSTRUMENTS:
                sum_rows.append(dict(
                    cell_line=line, instrument=inst, ok=False,
                    reason="trajectory state n<30",
                    n_random=N_RANDOM_DIR, n_null_finite=0,
                    p_more_negative=np.nan, p_more_positive=np.nan,
                ))
            continue

        Y = _load_Y(line)
        if Y.shape[0] != len(obs):
            raise StopStep("nulls", f"{line} louvain_rulerY n={Y.shape[0]} obs n={len(obs)}")
        obs_age, null_ruler = score_random_directions_cells(Y, frozen, log, tag=line)
        disk_age = np.asarray(scores_by[line]["frozen_ruler"], float)
        delta = np.nanmax(np.abs(obs_age - disk_age))
        if not np.isfinite(delta) or delta > 1e-8:
            raise StopStep(
                "nulls",
                f"{line}: recomputed frozen-ruler scores do not match md4 cell_scores. "
                f"max_abs_delta={delta}. Not substituting.",
            )
        log(f"[null ruler] {line} observed scores match md4 cell_scores max_abs_delta={delta}")
        meta[f"{line}_ruler"] = dict(
            max_abs_delta=float(delta), n_random=N_RANDOM_DIR, seed=PLANE_SEED,
        )
        y_ruler = _state_medians_1d(labels, null_ruler, STATE_NAMES)
        del null_ruler

        Xl, symbols, gene_id, _ = csr_from_npz(MD2_PROC / f"louvain_counts_{line}.npz")
        if Xl.shape[0] != len(obs):
            raise StopStep("nulls", f"{line} counts n={Xl.shape[0]} obs n={len(obs)}")
        logX = lognormalize_csr(Xl)
        bins = np.asarray(ams_by[line]["bins"], int)
        if bins.size != logX.shape[1]:
            raise StopStep(
                "nulls",
                f"{line}: louvain_ams bins n={bins.size} != logX p={logX.shape[1]}. Not substituting.",
            )
        n_unique_bins = int(np.unique(bins).size)
        if n_unique_bins != int(AMS_NBIN):
            raise StopStep("nulls", f"{line}: bins unique={n_unique_bins} != nbin={AMS_NBIN}")
        y_md, md_meta = _gene_list_null_state_y(
            logX, symbols, bins, sets["MD"], GENELIST_NULL_SEEDS["md_score"],
            labels, log, f"{line}:MD",
        )
        y_s3, s3_meta = _pair_null_state_y(
            logX, symbols, bins,
            sets["age_up"], GENELIST_NULL_SEEDS["age_up"],
            sets["age_down"], GENELIST_NULL_SEEDS["age_down"],
            labels, log, f"{line}:S3",
        )
        meta[f"{line}_md"] = md_meta
        meta[f"{line}_age_up_minus_age_down"] = s3_meta

        null_map = {
            "frozen_ruler": (y_ruler, "permuted_ruler_weights", PLANE_SEED),
            "md_score": (y_md, "size_matched_random_gene_sets_mean_expression_bin",
                         GENELIST_NULL_SEEDS["md_score"]),
            "age_up_minus_age_down": (
                y_s3, "size_matched_random_gene_sets_mean_expression_bin_up_minus_down",
                GENELIST_NULL_SEEDS["age_up"],
            ),
        }

        for inst in AGE_INSTRUMENTS:
            y_null, kind, nseed = null_map[inst]
            srow = sdf[(sdf.cell_line == line) & (sdf.instrument == inst)]
            if not len(srow):
                raise StopStep("nulls", f"missing slope row {line} {inst}")
            obs_bend = srow.iloc[0].get("bend")
            sA_n, sB_n, bend_n = _bends_from_state_y(y_null, x_med, a0, a1, end)
            finite = np.isfinite(bend_n)
            n_fin = int(finite.sum())
            if n_fin == 0:
                p_neg = np.nan
                p_pos = np.nan
                null_med = np.nan
                lo = np.nan
                hi = np.nan
            else:
                p_neg = permutation_p(obs_bend, bend_n[finite], greater=False)
                p_pos = permutation_p(obs_bend, bend_n[finite], greater=True)
                null_med = float(np.median(bend_n[finite]))
                lo = float(np.percentile(bend_n[finite], 2.5))
                hi = float(np.percentile(bend_n[finite], 97.5))
            sum_rows.append(dict(
                cell_line=line, instrument=inst, ok=True,
                null_kind=kind, null_seed=int(nseed),
                n_random=int(N_RANDOM_DIR), n_null_finite=n_fin,
                obs_bend=float(obs_bend) if np.isfinite(obs_bend) else np.nan,
                null_bend_median=null_med,
                null_bend_ci_lo=lo, null_bend_ci_hi=hi,
                p_more_negative=p_neg, p_more_positive=p_pos,
                segment_B_end=end, segment_B_end_is_fallback=fallback,
            ))
            log(f"[null] {line} {inst} n_fin={n_fin} obs_bend={obs_bend} "
                f"null_med={null_med} p_neg={p_neg} p_pos={p_pos} kind={kind}")
            for d in range(N_RANDOM_DIR):
                bend_rows.append(dict(
                    cell_line=line, instrument=inst, draw=int(d),
                    slope_A=float(sA_n[d]) if np.isfinite(sA_n[d]) else np.nan,
                    slope_B=float(sB_n[d]) if np.isfinite(sB_n[d]) else np.nan,
                    bend=float(bend_n[d]) if np.isfinite(bend_n[d]) else np.nan,
                ))
                for st in STATE_NAMES:
                    yv = y_null[st]
                    yv = float(yv[d]) if np.ndim(yv) else float(yv)
                    draw_rows.append(dict(
                        cell_line=line, instrument=inst, draw=int(d), label=st,
                        identity_median=_scalar_x(x_med, st),
                        age_median=yv,
                        on_trajectory=bool(st in TRAJECTORY),
                    ))
            for st in STATE_NAMES:
                yv = np.asarray(y_null[st], float)
                yv = yv[np.isfinite(yv)]
                rec = dict(
                    cell_line=line, instrument=inst, label=st,
                    n_null_finite=int(len(yv)),
                    identity_median=_scalar_x(x_med, st),
                    on_trajectory=bool(st in TRAJECTORY),
                )
                if len(yv):
                    rec.update(
                        age_p025=float(np.percentile(yv, 2.5)),
                        age_p50=float(np.percentile(yv, 50)),
                        age_p975=float(np.percentile(yv, 97.5)),
                    )
                else:
                    rec.update(age_p025=np.nan, age_p50=np.nan, age_p975=np.nan)
                env_rows.append(rec)

        dump_json(PLANE_DIR / f"t1_null_meta_{line}.json", jsonable({
            k: meta[k] for k in meta if k.startswith(line)
        }))

    sdf_sum = pd.DataFrame(sum_rows)
    sdf_sum.to_csv(PLANE_DIR / "t1_null_summary.csv", index=False)
    pd.DataFrame(bend_rows).to_csv(PLANE_DIR / "t1_null_bend.csv", index=False)
    pd.DataFrame(draw_rows).to_csv(PLANE_DIR / "t1_null_state_medians.csv", index=False)
    env = pd.DataFrame(env_rows)
    env.to_csv(PLANE_DIR / "t1_null_envelope.csv", index=False)
    dump_json(PLANE_DIR / "t1_null_meta.json", jsonable(meta))
    return sdf_sum, env, pd.DataFrame(bend_rows)

"""STAGE L2 — refit the age direction in each low-dim space; re-test identifiability.

L2a donor-bootstrap pairwise angle (gene space was 45.1°).
L2b held-out age prediction, both CV schemes, permutation null.
L2c site transfer HBCC→MSSM and MSSM→HBCC (gene space −0.213 / −0.499).
L2d age-range split with bank control (gene space 86.1°, confounded with bank).

Primary selection: L2a and L2c together, not prediction R².
STOP L2 if no variant/k achieves both (i) bootstrap angle materially below
gene space (≤ 35°) and (ii) positive site transfer in at least one direction.

Dimensionality reduction is fit on the training rows only.

Usage: python src/lowdim_l2.py
"""
from __future__ import annotations

import time
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lowdim_common import (  # noqa: E402
    LD_DIR, LD_FIG, LD_SEED, DATASET_ID, KS, KMAX, VARIANTS, VARIANT_LABEL,
    N_BOOT, N_PERM, MIN_N_TYPE, ANGLE_MATERIAL, SITE_NAME, STOP_NULL,
    Logger, dump_json, ld_log_banner, load_bundle, load_gene_sets,
    zscore_train, make_folds, fit_spaces_by_type, apply_slice, age_direction,
    gene_direction_from_pack, direction_and_pred, donor_bootstrap_index,
    median_pairwise_unit_cols, permute_age_within_site_obs, score_age_predictions,
    permutation_p, summarize_null_col, config_id, gene_space_baseline,
    unsigned_angle_deg, r2_mae, FilterLog, fmt, fit_pack, TypeProjector,
)


def _configs():
    for variant in VARIANTS:
        for k in KS:
            for clean in (False, True):
                yield variant, int(k), bool(clean)


def _fit_pack(Y, obs, tech, idx, types, variant, sets):
    return fit_pack(Y, obs, tech, idx, types, variant, sets)


def _projector(pack, t, Xte_t=None):
    sp = pack["spaces"].get(t)
    tr_m = pack["ct"] == t
    n_tr = int(tr_m.sum())
    if sp is None or n_tr < MIN_N_TYPE:
        return None
    return TypeProjector(sp, pack["Xz"][tr_m], pack["tech"][tr_m], Xte=Xte_t), pack["y"][tr_m], n_tr, sp.p


def _gene_dir_type(pack, t, k, clean, proj=None):
    if proj is None:
        got = _projector(pack, t)
        if got is None:
            return None
        proj, ytr, _n, p = got
    else:
        ytr = pack["y"][pack["ct"] == t]
        p = pack["Xz"].shape[1]
    sl = proj.pack(k, clean)
    if sl["k_keep"] < 1:
        return None
    w, rec = age_direction(sl["Ztr"], ytr, method="ridge")
    if not rec.get("ok"):
        return None
    v = gene_direction_from_pack(sl, w, p)
    if not np.isfinite(v).all() or float(np.linalg.norm(v)) < 1e-12:
        return None
    return v / np.linalg.norm(v)


def _predict_split(pack, Y, obs, te_idx, types, k, clean, projs=None):
    te_idx = np.asarray(te_idx)
    Xte = (Y[te_idx] - pack["mu"]) / pack["sd"]
    ct_te = obs.celltype.astype(str).to_numpy()[te_idx]
    yte = obs.age.to_numpy(float)[te_idx]
    pred = np.full(len(te_idx), np.nan)
    rows = []
    if projs is None:
        projs = {}
        for t in types:
            te_m = ct_te == t
            projs[t] = _projector(pack, t, Xte[te_m])
    for t in types:
        te_m = ct_te == t
        if projs.get(t) is None:
            rows.append(dict(celltype=t, r2=np.nan, mae=np.nan, n_te=int(te_m.sum()), n_tr=0, k_keep=0))
            continue
        proj, ytr, n_tr, _p = projs[t]
        # ytr from projector construction; refresh from pack (perms)
        ytr = pack["y"][pack["ct"] == t]
        sl = proj.pack(k, clean)
        if sl["k_keep"] < 1 or int(te_m.sum()) < 3:
            rows.append(dict(celltype=t, r2=np.nan, mae=np.nan, n_te=int(te_m.sum()),
                             n_tr=n_tr, k_keep=int(sl["k_keep"])))
            continue
        _w, rec, pte = direction_and_pred(sl["Ztr"], ytr, sl["Zte"])
        pred[te_m] = pte
        met = r2_mae(yte[te_m], pte)
        rows.append(dict(celltype=t, r2=met["r2"], mae=met["mae"], n_te=met["n"], n_tr=n_tr,
                         k_keep=int(sl["k_keep"]), reverted=bool(sl["reverted"]),
                         alpha=rec.get("alpha", np.nan)))
    tab = pd.DataFrame(rows)
    med = float(tab.r2.median()) if tab.r2.notna().any() else np.nan
    return pred, yte, tab, med


def _l2a(Y, obs, tech, types, sets, rng, log):
    log(f"\n[L2a] donor bootstrap n_boot={N_BOOT}  seed={LD_SEED}")
    p = Y.shape[1]
    ct = obs.celltype.astype(str).to_numpy()
    rows = []
    t0 = time.time()
    for variant in VARIANTS:
        log(f"   variant={variant}")
        store = {(k, clean, t): [] for k, clean in [(k, c) for k in KS for c in (False, True)] for t in types}
        for b in range(N_BOOT):
            idx = donor_bootstrap_index(obs, rng)
            pack = _fit_pack(Y, obs, tech, idx, types, variant, sets)
            for t in types:
                got = _projector(pack, t)
                if got is None:
                    continue
                proj, ytr, _n, p = got
                for k in KS:
                    for clean in (False, True):
                        sl = proj.pack(k, clean)
                        if sl["k_keep"] < 1:
                            continue
                        w, rec = age_direction(sl["Ztr"], ytr, method="ridge")
                        if not rec.get("ok"):
                            continue
                        v = gene_direction_from_pack(sl, w, p)
                        if not np.isfinite(v).all() or float(np.linalg.norm(v)) < 1e-12:
                            continue
                        store[(k, clean, t)].append(v / np.linalg.norm(v))
            if (b + 1) % 5 == 0 or b == 0:
                log(f"      {variant} boot {b+1}/{N_BOOT}  ({time.time()-t0:.0f}s)")
        for k in KS:
            for clean in (False, True):
                per_t = []
                for t in types:
                    cols = store[(k, clean, t)]
                    med, n_pairs = median_pairwise_unit_cols(cols)
                    per_t.append(dict(
                        variant=variant, k=int(k), clean=bool(clean),
                        config=config_id(variant, k, clean), celltype=t,
                        n_boot_ok=len(cols), n_boot=N_BOOT,
                        median_pairwise_deg=med, n_pairs=n_pairs, n_genes=int(p),
                    ))
                    rows.append(per_t[-1])
                meds = [r["median_pairwise_deg"] for r in per_t if np.isfinite(r["median_pairwise_deg"])]
                log(f"      {config_id(variant, k, clean)} median-over-types angle="
                    f"{float(np.median(meds)) if meds else np.nan:.1f}°  n_types={len(meds)}")
    boot_df = pd.DataFrame(rows)
    boot_df.to_csv(LD_DIR / "l2a_bootstrap_per_type.csv", index=False)
    g = (boot_df.groupby(["variant", "k", "clean", "config"], as_index=False)
         .agg(median_pairwise_deg=("median_pairwise_deg", "median"),
              n_types=("median_pairwise_deg", "count"),
              p05=("median_pairwise_deg", lambda s: float(np.nanpercentile(s, 5)) if s.notna().any() else np.nan),
              p95=("median_pairwise_deg", lambda s: float(np.nanpercentile(s, 95)) if s.notna().any() else np.nan)))
    g.to_csv(LD_DIR / "l2a_bootstrap_summary.csv", index=False)
    log(f"[L2a] done in {time.time()-t0:.0f}s")
    return boot_df, g


def _scheme_primary(sc):
    """Hard rule: site-stratified R², never pooled."""
    return sc["site_strat_median_type_r2"]


def _l2b(Y, obs, tech, types, sets, rng, log):
    log("\n[L2b] held-out prediction (train-fold spaces only)")
    folds_all = make_folds(obs, np.random.default_rng(LD_SEED))
    y = obs.age.to_numpy(float)
    ct = obs.celltype.astype(str).to_numpy()
    score_rows, null_rows, pt_rows = [], [], []
    t0 = time.time()
    for scheme, folds in (("within_site", folds_all["within_site"]), ("loso", folds_all["loso"])):
        log(f"   scheme={scheme}  n_folds={len(folds)}")
        for variant in VARIANTS:
            log(f"      variant={variant}  fitting fold spaces...")
            fold_pack = []
            fold_slice = []  # list of {t: { (k,clean): sl } }
            for fi, fold in enumerate(folds):
                tr = fold[0]
                te = fold[1]
                pack = _fit_pack(Y, obs, tech, tr, types, variant, sets)
                Xte = (Y[te] - pack["mu"]) / pack["sd"]
                ct_te = ct[te]
                proj_map = {}
                for t in types:
                    te_m = ct_te == t
                    got = _projector(pack, t, Xte[te_m])
                    proj_map[t] = got
                fold_pack.append((pack, te, tr))
                fold_slice.append(proj_map)
                log(f"         fold {fi+1}/{len(folds)}  ({time.time()-t0:.0f}s)")

            def _oof(y_use):
                pred = { (k, clean): np.full(len(obs), np.nan) for k in KS for clean in (False, True) }
                for fi, fold in enumerate(folds):
                    te = fold[1]
                    ct_te = ct[te]
                    pack = fold_pack[fi][0]
                    proj_map = fold_slice[fi]
                    ytr_all = y_use[pack["idx"]]
                    for t in types:
                        tr_m = pack["ct"] == t
                        te_m = ct_te == t
                        ytr = ytr_all[tr_m]
                        got = proj_map.get(t)
                        if got is None or int(te_m.sum()) < 1:
                            continue
                        proj, _ytr_orig, _n, _p = got
                        for k in KS:
                            for clean in (False, True):
                                sl = proj.pack(k, clean)
                                if sl["k_keep"] < 1:
                                    continue
                                _w, rec, pte = direction_and_pred(sl["Ztr"], ytr, sl["Zte"])
                                pred[(k, clean)][te[te_m]] = pte
                return pred

            pred_obs = _oof(y)
            log(f"      {variant} permutation null ({N_PERM})...")
            null_pred = []
            for i in range(N_PERM):
                obs_p = permute_age_within_site_obs(obs, rng)
                yp = obs_p.age.to_numpy(float)
                null_pred.append((yp, _oof(yp)))
                if (i + 1) % 5 == 0:
                    log(f"         perm {i+1}/{N_PERM}  ({time.time()-t0:.0f}s)")

            for k in KS:
                for clean in (False, True):
                    pred = pred_obs[(k, clean)]
                    sc = score_age_predictions(obs, pred)
                    prim = _scheme_primary(sc)
                    pooled = sc["pooled_median_type_r2"]
                    null_prim = []
                    for yp, pmap in null_pred:
                        scn = score_age_predictions(obs, pmap[(k, clean)], y=yp)
                        null_prim.append(_scheme_primary(scn))
                        null_rows.append(dict(
                            scheme=scheme, variant=variant, k=int(k), clean=bool(clean),
                            config=config_id(variant, k, clean),
                            perm=len(null_prim) - 1, primary=_scheme_primary(scn),
                            pooled=scn["pooled_median_type_r2"],
                            site_strat=scn["site_strat_median_type_r2"],
                        ))
                    pval = permutation_p(prim, null_prim, greater=True)
                    ns = summarize_null_col(null_prim)
                    cfg = config_id(variant, k, clean)
                    log(f"      {scheme} {cfg}: site-strat R²={prim:+.3f}  pooled={pooled:+.3f}  "
                        f"null={ns['mean']:+.3f}  p={pval:.3f}")
                    score_rows.append(dict(
                        scheme=scheme, variant=variant, k=int(k), clean=bool(clean),
                        config=cfg, site_strat_r2=prim, pooled_r2=pooled,
                        null_mean=ns["mean"], null_p95=ns["p95"], p=pval, n_perm=N_PERM,
                        null_usable=bool(np.isfinite(ns["mean"]) and ns["mean"] <= STOP_NULL),
                    ))
                    pt = sc["per_type"]
                    sub = pt[pt.site == "mean_of_sites"] if scheme == "within_site" else pt[pt.site == "mean_of_sites"]
                    # always site-stratified per type
                    for _, row in pt[pt.site == "mean_of_sites"].iterrows():
                        pt_rows.append(dict(
                            scheme=scheme, variant=variant, k=int(k), clean=bool(clean),
                            config=cfg, celltype=row.celltype, r2=row.r2, mae=row.mae, n=row.n,
                        ))
    sc_df = pd.DataFrame(score_rows)
    sc_df.to_csv(LD_DIR / "l2b_oof_summary.csv", index=False)
    pd.DataFrame(null_rows).to_csv(LD_DIR / "l2b_oof_null.csv", index=False)
    pd.DataFrame(pt_rows).to_csv(LD_DIR / "l2b_oof_per_type.csv", index=False)
    log(f"[L2b] done in {time.time()-t0:.0f}s")
    return sc_df


def _l2c(Y, obs, tech, types, sets, rng, log):
    log("\n[L2c] site transfer (fit space + direction on one bank, score the other)")
    site = obs.Source.astype(str).to_numpy()
    sites = sorted(pd.unique(site))
    rows, pt_rows, null_rows = [], [], []
    t0 = time.time()
    for tr_s in sites:
        te_s = [s for s in sites if s != tr_s][0]
        tr_idx = np.flatnonzero(site == tr_s)
        te_idx = np.flatnonzero(site == te_s)
        log(f"   train {SITE_NAME.get(tr_s, tr_s)} n_donors="
            f"{obs.iloc[tr_idx].donor.nunique()} → test {SITE_NAME.get(te_s, te_s)} "
            f"n_donors={obs.iloc[te_idx].donor.nunique()}")
        for variant in VARIANTS:
            pack = _fit_pack(Y, obs, tech, tr_idx, types, variant, sets)
            obs_tr = obs.iloc[tr_idx].reset_index(drop=True)
            Xte = (Y[te_idx] - pack["mu"]) / pack["sd"]
            ct_te = obs.celltype.astype(str).to_numpy()[te_idx]
            projs = {t: _projector(pack, t, Xte[ct_te == t]) for t in types}
            log(f"      {variant} observed + {N_PERM} train-age perms...")
            obs_cache = {}
            for k in KS:
                for clean in (False, True):
                    pred, yte, tab, med = _predict_split(pack, Y, obs, te_idx, types, k, clean, projs=projs)
                    tab = tab.assign(train_site=tr_s, test_site=te_s, variant=variant,
                                     k=int(k), clean=bool(clean), config=config_id(variant, k, clean))
                    pt_rows.append(tab)
                    obs_cache[(k, clean)] = med
            null_med = {(k, c): [] for k in KS for c in (False, True)}
            for i in range(N_PERM):
                obs_p = permute_age_within_site_obs(obs_tr, rng)
                pack_p = dict(pack)
                pack_p["y"] = obs_p.age.to_numpy(float)
                for k in KS:
                    for clean in (False, True):
                        _pred, _yte, _tab, mnull = _predict_split(
                            pack_p, Y, obs, te_idx, types, k, clean, projs=projs)
                        null_med[(k, clean)].append(mnull)
                        null_rows.append(dict(
                            train_site=tr_s, test_site=te_s, variant=variant, k=int(k),
                            clean=bool(clean), config=config_id(variant, k, clean),
                            perm=i, median_r2=mnull,
                        ))
            for k in KS:
                for clean in (False, True):
                    med = obs_cache[(k, clean)]
                    nm = null_med[(k, clean)]
                    pval = permutation_p(med, nm, greater=True)
                    ns = summarize_null_col(nm)
                    cfg = config_id(variant, k, clean)
                    log(f"      {SITE_NAME.get(tr_s, tr_s)}→{SITE_NAME.get(te_s, te_s)} {cfg}: "
                        f"R²={med:+.3f}  null={ns['mean']:+.3f}  p={pval:.3f}")
                    rows.append(dict(
                        train_site=tr_s, test_site=te_s, variant=variant, k=int(k),
                        clean=bool(clean), config=cfg, median_r2=med,
                        null_mean=ns["mean"], null_p95=ns["p95"], p=pval, n_perm=N_PERM,
                        n_train_donors=int(obs.iloc[tr_idx].donor.nunique()),
                        n_test_donors=int(obs.iloc[te_idx].donor.nunique()),
                    ))
            log(f"      {variant} {tr_s}→{te_s} done ({time.time()-t0:.0f}s)")
    rec = pd.DataFrame(rows)
    rec.to_csv(LD_DIR / "l2c_transfer_summary.csv", index=False)
    pd.concat(pt_rows, ignore_index=True).to_csv(LD_DIR / "l2c_transfer_per_type.csv", index=False)
    pd.DataFrame(null_rows).to_csv(LD_DIR / "l2c_transfer_null.csv", index=False)
    log(f"[L2c] done in {time.time()-t0:.0f}s")
    return rec


def _age_masks(obs, log, rng):
    don = obs.groupby("donor").agg(age=("age", "first"), site=("Source", "first"))
    donors = obs.donor.astype(str).to_numpy()
    site = obs.Source.astype(str).to_numpy()
    age_d = obs.age.to_numpy(float)
    out = {}
    for s in sorted(pd.unique(site)):
        young = set(don.index[(don.site.astype(str) == s) & (don.age < 50)].astype(str))
        old = set(don.index[(don.site.astype(str) == s) & (don.age >= 50)].astype(str))
        m_y = np.array([d in young for d in donors])
        m_o = np.array([d in old for d in donors])
        log(f"[L2d] within {SITE_NAME.get(s, s)}: young donors={len(young)} rows={int(m_y.sum())}  "
            f"old donors={len(old)} rows={int(m_o.sum())}")
        out[f"within_{s}"] = (m_y, m_o, dict(n_young=len(young), n_old=len(old), site=s))
    # balanced: equal young/old within each bank, then pool
    y_keep, o_keep = [], []
    for s in sorted(pd.unique(don.site.astype(str))):
        y = don.index[(don.site.astype(str) == s) & (don.age < 50)].astype(str).to_numpy()
        o = don.index[(don.site.astype(str) == s) & (don.age >= 50)].astype(str).to_numpy()
        n = int(min(len(y), len(o)))
        log(f"[L2d] balance {SITE_NAME.get(s, s)}: min(young,old)={n} (young={len(y)} old={len(o)})")
        if n < 5:
            continue
        y_keep.extend(rng.choice(y, size=n, replace=False).tolist())
        o_keep.extend(rng.choice(o, size=n, replace=False).tolist())
    y_keep, o_keep = set(y_keep), set(o_keep)
    m_y = np.array([d in y_keep for d in donors])
    m_o = np.array([d in o_keep for d in donors])
    y_sites = don.loc[list(y_keep)].groupby("site").size().to_dict() if y_keep else {}
    o_sites = don.loc[list(o_keep)].groupby("site").size().to_dict() if o_keep else {}
    log(f"[L2d] balanced pool: young donors={len(y_keep)} sites={y_sites}  "
        f"old donors={len(o_keep)} sites={o_sites}")
    out["balanced"] = (m_y, m_o, dict(n_young=len(y_keep), n_old=len(o_keep),
                                      young_sites=y_sites, old_sites=o_sites))
    return out


def _l2d(Y, obs, tech, types, sets, rng, log):
    log("\n[L2d] age-range split, bank-controlled")
    masks = _age_masks(obs, log, rng)
    rows = []
    t0 = time.time()
    for split, (m_y, m_o, meta) in masks.items():
        log(f"   split={split}")
        iy, io = np.flatnonzero(m_y), np.flatnonzero(m_o)
        if len(iy) < 20 or len(io) < 20:
            log(f"      SKIP too few rows young={len(iy)} old={len(io)}")
            continue
        for variant in VARIANTS:
            pack_y = _fit_pack(Y, obs, tech, iy, types, variant, sets)
            pack_o = _fit_pack(Y, obs, tech, io, types, variant, sets)
            py = {t: _projector(pack_y, t) for t in types}
            po = {t: _projector(pack_o, t) for t in types}
            for k in KS:
                for clean in (False, True):
                    angs = []
                    for t in types:
                        vy = _gene_dir_type(pack_y, t, k, clean, proj=None if py[t] is None else py[t][0])
                        vo = _gene_dir_type(pack_o, t, k, clean, proj=None if po[t] is None else po[t][0])
                        n_y = int((pack_y["ct"] == t).sum())
                        n_o = int((pack_o["ct"] == t).sum())
                        ok = vy is not None and vo is not None
                        ang = unsigned_angle_deg(vy, vo) if ok else np.nan
                        if ok:
                            angs.append(ang)
                        rows.append(dict(
                            split=split, variant=variant, k=int(k), clean=bool(clean),
                            config=config_id(variant, k, clean), celltype=t,
                            n_young=n_y, n_old=n_o, angle_deg=ang, ok=bool(ok),
                            n_donors_young=meta.get("n_young"), n_donors_old=meta.get("n_old"),
                        ))
                    med = float(np.median(angs)) if angs else np.nan
                    log(f"      {split} {config_id(variant, k, clean)} median angle={med:.1f}°  n_types={len(angs)}")
            log(f"      {split} {variant} done ({time.time()-t0:.0f}s)")
    df = pd.DataFrame(rows)
    df.to_csv(LD_DIR / "l2d_agerange_per_type.csv", index=False)
    g = (df.groupby(["split", "variant", "k", "clean", "config"], as_index=False)
         .agg(median_angle_deg=("angle_deg", "median"),
              n_types=("ok", "sum")))
    g.to_csv(LD_DIR / "l2d_agerange_summary.csv", index=False)
    log(f"[L2d] done in {time.time()-t0:.0f}s")
    return df, g


def _comparison(boot_sum, oof, trans, age_sum, base, log):
    rows = []
    # gene-space reference
    rows.append(dict(
        config="gene_ridge", variant="gene", k=25526, clean=False,
        boot_median_angle=base["boot_median_angle"],
        within_site_r2=base["within_site_r2_v2"],
        within_site_r2_g1=base["within_site_r2_g1"],
        loso_r2=np.nan, loso_null=np.nan,
        delta_r2_vs_v2=0.0, delta_r2_vs_g1=base["within_site_r2_v2"] - base["within_site_r2_g1"],
        transfer_H_to_M=base["transfer_H_to_M"], transfer_M_to_H=base["transfer_M_to_H"],
        transfer_any_positive=bool(base["transfer_H_to_M"] > 0 or base["transfer_M_to_H"] > 0),
        young_old_within_H=np.nan, young_old_within_M=np.nan, young_old_balanced=base["young_old_angle"],
        pass_angle=False, pass_transfer=False, pass_both=False,
        note="gene-space TARGET/V1 from FINDINGS_TARGET (not refit)",
    ))
    oof_ws = oof[oof.scheme == "within_site"].drop_duplicates("config").set_index("config")
    oof_lo = oof[oof.scheme == "loso"].drop_duplicates("config").set_index("config")
    boot_sum = boot_sum.copy()
    if "config" not in boot_sum.columns:
        boot_sum["config"] = [config_id(v, k, c) for v, k, c in zip(boot_sum.variant, boot_sum.k, boot_sum.clean)]
    bmap = boot_sum.drop_duplicates("config").set_index("config")
    tH = trans[(trans.train_site == "H") & (trans.test_site == "M")].drop_duplicates("config").set_index("config")
    tM = trans[(trans.train_site == "M") & (trans.test_site == "H")].drop_duplicates("config").set_index("config")

    def _age(split, cfg):
        sub = age_sum[(age_sum.split == split) & (age_sum.config == cfg)]
        return float(sub.median_angle_deg.iloc[0]) if len(sub) else np.nan

    for variant, k, clean in _configs():
        cfg = config_id(variant, k, clean)
        ang = float(bmap.at[cfg, "median_pairwise_deg"]) if cfg in bmap.index else np.nan
        ws = float(oof_ws.at[cfg, "site_strat_r2"]) if cfg in oof_ws.index else np.nan
        lo = float(oof_lo.at[cfg, "site_strat_r2"]) if cfg in oof_lo.index else np.nan
        lo_null = float(oof_lo.at[cfg, "null_mean"]) if cfg in oof_lo.index else np.nan
        ws_null = float(oof_ws.at[cfg, "null_mean"]) if cfg in oof_ws.index else np.nan
        hm = float(tH.at[cfg, "median_r2"]) if cfg in tH.index else np.nan
        mh = float(tM.at[cfg, "median_r2"]) if cfg in tM.index else np.nan
        pass_a = bool(np.isfinite(ang) and ang <= ANGLE_MATERIAL)
        pass_t = bool((np.isfinite(hm) and hm > 0) or (np.isfinite(mh) and mh > 0))
        rows.append(dict(
            config=cfg, variant=variant, k=int(k), clean=bool(clean),
            boot_median_angle=ang,
            within_site_r2=ws, within_site_null=ws_null,
            within_site_r2_g1=base["within_site_r2_g1"],
            loso_r2=lo, loso_null=lo_null,
            delta_r2_vs_v2=(ws - base["within_site_r2_v2"]) if np.isfinite(ws) else np.nan,
            delta_r2_vs_g1=(ws - base["within_site_r2_g1"]) if np.isfinite(ws) else np.nan,
            transfer_H_to_M=hm, transfer_M_to_H=mh,
            transfer_any_positive=pass_t,
            young_old_within_H=_age("within_H", cfg),
            young_old_within_M=_age("within_M", cfg),
            young_old_balanced=_age("balanced", cfg),
            pass_angle=pass_a, pass_transfer=pass_t, pass_both=bool(pass_a and pass_t),
            note="",
        ))
    tab = pd.DataFrame(rows)
    tab.to_csv(LD_DIR / "l2_comparison.csv", index=False)

    cand = tab[tab.pass_both & (tab.config != "gene_ridge")]
    winner = None
    reason = ""
    if len(cand):
        # rank by bootstrap angle (lower better) and best one-way transfer (higher better)
        c = cand.copy()
        c["best_transfer"] = c[["transfer_H_to_M", "transfer_M_to_H"]].max(axis=1)
        c["rank_angle"] = c.boot_median_angle.rank(method="min")
        c["rank_transfer"] = c.best_transfer.rank(ascending=False, method="min")
        c["score"] = c.rank_angle + c.rank_transfer
        c = c.sort_values(["score", "boot_median_angle"])
        winner = c.iloc[0].to_dict()
        reason = ("lowest combined rank of bootstrap angle and best one-way site transfer "
                  "among configs that pass both bars")
        # prefer not-cleaned if tied-ish? already in score
    stop = winner is None
    reason = reason if winner else (
        "no variant/k achieved bootstrap angle materially below gene space AND "
        "positive site transfer in at least one direction — direction is not identifiable even at low dimension"
    )
    (LD_DIR / "l2_STOP.txt").write_text(
        f"stop_L2={stop}\n"
        f"angle_bar_deg={ANGLE_MATERIAL}\n"
        f"gene_boot_angle={base['boot_median_angle']}\n"
        f"n_pass_angle={int(tab.pass_angle.sum())}\n"
        f"n_pass_transfer={int(tab.pass_transfer.sum())}\n"
        f"n_pass_both={int(tab.pass_both.sum())}\n"
        f"winner={winner['config'] if winner else 'NONE'}\n"
        f"reason={'no variant/k achieved bootstrap angle materially below gene space AND positive site transfer in at least one direction — direction is not identifiable even at low dimension' if stop else reason}\n",
        encoding="utf-8",
    )
    log(f"[STOP L2] {stop}  winner={winner['config'] if winner else 'NONE'}")
    if winner:
        log(f"   boot={winner['boot_median_angle']:.1f}°  "
            f"H→M={winner['transfer_H_to_M']:+.3f}  M→H={winner['transfer_M_to_H']:+.3f}  "
            f"within-site R²={winner['within_site_r2']:+.3f}  "
            f"(delta vs V2 {winner['delta_r2_vs_v2']:+.3f})")
        if winner["within_site_r2"] < base["within_site_r2_v2"]:
            log("   TRADEOFF: winning space predicts worse than gene-space V2 TARGET; "
                "stability/transfer are the selection criteria, not R².")
    dump_json(LD_DIR / "l2_winner.json", dict(stop_L2=stop, winner=winner, reason=reason, bar_deg=ANGLE_MATERIAL))
    return tab, winner, stop


def _figures(boot_sum, oof, trans, age_sum, tab):
    colors = dict(pca="tab:blue", curated="tab:green", wgcna="tab:orange")
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.0), sharey=True)
    for ax, clean, title in zip(axes, (False, True), ("all components", "technical components dropped")):
        for variant, c in colors.items():
            sub = boot_sum[(boot_sum.variant == variant) & (boot_sum.clean == clean)].sort_values("k")
            if len(sub):
                ax.plot(sub.k, sub.median_pairwise_deg, "-o", color=c, label=VARIANT_LABEL[variant])
        ax.axhline(45.1, c="k", ls="--", lw=1, label="gene space 45.1°")
        ax.axhline(ANGLE_MATERIAL, c="0.4", ls=":", lw=1, label=f"material bar {ANGLE_MATERIAL:.0f}°")
        ax.set_xlabel("k")
        ax.set_title(title)
        ax.set_ylabel("median bootstrap pairwise angle (deg)")
        ax.legend(fontsize=6)
    fig.suptitle("L2a  donor-bootstrap stability", fontsize=11)
    fig.tight_layout()
    fig.savefig(LD_FIG / "l2a_bootstrap_vs_k.png", dpi=140)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.0), sharey=True)
    for ax, direction, title in zip(
        axes,
        ((trans.train_site == "H") & (trans.test_site == "M"),
         (trans.train_site == "M") & (trans.test_site == "H")),
        ("HBCC → MSSM", "MSSM → HBCC"),
    ):
        sub0 = trans[direction]
        for variant, c in colors.items():
            for clean, ls in ((False, "-"), (True, "--")):
                sub = sub0[(sub0.variant == variant) & (sub0.clean == clean)].sort_values("k")
                if len(sub):
                    ax.plot(sub.k, sub.median_r2, ls + "o", color=c,
                            label=f"{VARIANT_LABEL[variant]}{'_clean' if clean else ''}")
        ax.axhline(0, c="0.6", lw=0.8)
        ax.axhline(-0.213 if "HBCC" in title else -0.499, c="k", ls="--", lw=1, label="gene TARGET")
        ax.set_xlabel("k")
        ax.set_title(title)
        ax.set_ylabel("median-over-types R²")
        ax.legend(fontsize=6)
    fig.suptitle("L2c  site transfer", fontsize=11)
    fig.tight_layout()
    fig.savefig(LD_FIG / "l2c_transfer_vs_k.png", dpi=140)
    plt.close(fig)

    ws = oof[oof.scheme == "within_site"]
    fig, ax = plt.subplots(figsize=(6.8, 4.0))
    for variant, c in colors.items():
        for clean, ls in ((False, "-"), (True, "--")):
            sub = ws[(ws.variant == variant) & (ws.clean == clean)].sort_values("k")
            if len(sub):
                ax.plot(sub.k, sub.site_strat_r2, ls + "o", color=c,
                        label=f"{VARIANT_LABEL[variant]}{'_clean' if clean else ''}")
    ax.axhline(0.231, c="k", ls="--", lw=1, label="gene TARGET 0.231")
    ax.axhline(0.317, c="0.4", ls=":", lw=1, label="G1 ridge 0.317")
    ax.set_xlabel("k")
    ax.set_ylabel("within-site site-stratified R²")
    ax.set_title("L2b  held-out age prediction (within-site)")
    ax.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(LD_FIG / "l2b_within_site_r2.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    for split, lab, col in (("within_H", "within HBCC", "tab:blue"),
                            ("within_M", "within MSSM", "tab:orange"),
                            ("balanced", "balanced banks", "tab:green")):
        sub = age_sum[age_sum.split == split]
        if not len(sub):
            continue
        g = sub[sub.clean == False].groupby(["variant", "k"]).median_angle_deg.mean().reset_index()
        for variant, c in colors.items():
            gg = g[g.variant == variant].sort_values("k")
            if len(gg):
                ax.plot(gg.k, gg.median_angle_deg, "-o", color=c, alpha=0.9 if split == "balanced" else 0.45)
    ax.axhline(86.1, c="k", ls="--", lw=1, label="gene space 86.1° (bank-confounded)")
    ax.set_xlabel("k")
    ax.set_ylabel("median young vs old angle (deg)")
    ax.set_title("L2d  age-range split (faint=within-bank, solid-ish=see legend in report)")
    fig.tight_layout()
    fig.savefig(LD_FIG / "l2d_agerange_vs_k.png", dpi=140)
    plt.close(fig)


def _run(log, rng):
    ld_log_banner(log, "L2 — identifiability in low-dimensional spaces")
    flog = FilterLog()
    data = load_bundle(log)
    Y, genes, obs, tech, types = data["Y"], data["genes"], data["obs"], data["tech"], data["types"]
    sets, _ = load_gene_sets(genes, log)
    flog.add("L2", "pseudobulks", len(obs), n_genes=int(Y.shape[1]), n_donors=int(data["n_donors"]))
    base = gene_space_baseline()
    log(f"[L2] gene-space baseline: boot={base['boot_median_angle']:.1f}°  "
        f"V2 within-site={base['within_site_r2_v2']:+.3f}  G1={base['within_site_r2_g1']:+.3f}  "
        f"H→M={base['transfer_H_to_M']:+.3f}  M→H={base['transfer_M_to_H']:+.3f}  "
        f"source={base['source']}")
    log(f"[L2] STOP bars: bootstrap angle ≤ {ANGLE_MATERIAL:.0f}° AND at least one transfer R² > 0")

    boot_df, boot_sum = _l2a(Y, obs, tech, types, sets, rng, log)
    oof = _l2b(Y, obs, tech, types, sets, rng, log)
    trans = _l2c(Y, obs, tech, types, sets, rng, log)
    age_df, age_sum = _l2d(Y, obs, tech, types, sets, rng, log)
    tab, winner, stop = _comparison(boot_sum, oof, trans, age_sum, base, log)
    _figures(boot_sum, oof, trans, age_sum, tab)
    flog.to_csv(LD_DIR / "l2_filter_log.csv")
    summary = dict(
        seed=LD_SEED, dataset_id=DATASET_ID, n_boot=N_BOOT, n_perm=N_PERM,
        angle_bar_deg=ANGLE_MATERIAL, baseline=base,
        stop_L2=stop, winner=winner,
        n_pass_both=int(tab.pass_both.sum()),
    )
    dump_json(LD_DIR / "l2_summary.json", summary)
    log("[L2] wrote results/lowdim/l2_*")
    log("[L2] done.")
    return summary


def run():
    rng = np.random.default_rng(LD_SEED)
    log = Logger(LD_DIR / "l2_report.txt")
    try:
        return _run(log, rng)
    finally:
        log.close()


if __name__ == "__main__":
    run()

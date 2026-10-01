"""Stage 1 — GTEx identifiability battery.

Gated on FINDINGS_POSCTRL.md Supersession 2 (Pearson r). Primary transfer metric
is held-out Pearson r with donor-level permutation null (n_perm=200).
Pass = r > 0 with null ≤ 0.05. Uncalibrated R², calibrated R², and ρ are reported.
Bootstrap angle is reported only, not a gate.

Usage: python src/gtex_stage1.py [--cell T1|T2|T3|T4|T5|T6|COMP|all]
Does not modify prior FINDINGS*.md or FALSIFICATION.md.
Fails loudly on held-out assertions. Does not repair folds.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gtex_common import (  # noqa: E402
    GTEX_DIR, GTEX_PROC, GTEX_SEED, TISSUES, PRIMARY_TISSUE, PAIRED_TISSUE,
    CROSS_REGION, TISSUE_SLUG, BRAIN_TISSUES, N_BOOT, N_PERM_GTEX, N_FOLDS,
    SEX_ANGLE_REF, PLANTED_ANGLE_N233, T6_MIN_OVERLAP, CELL_SEEDS, PREREG_FLAG,
    PREREG_DATE, Logger, dump_json, load_json, StopStep, gtex_log_banner,
    stage1_gate, onehot_train, apply_onehot, residualize_train_test, zscore_train,
    age_direction, unsigned_angle_deg, pairwise_angles_deg, permutation_p,
    strip_ensembl, TARGET_DIR, identity_basis, angle_to_subspace_deg, proj_r2,
    assert_disjoint, grouped_kfold, pred_scores, pass_r_bar, unit, align_age_sign,
    ridge_prestd, progress_snapshot, summarize_null_col,
)
from brain_phase1_common import r2_mae  # noqa: E402
from trajectory_common import pls1_predict  # noqa: E402
from target_common import donor_bootstrap_index  # noqa: E402


REGIMES = ("raw", "resid")
METHODS = ("ridge", "pls1")


def _gate_or_stop(log):
    gate = stage1_gate()
    dump_json(GTEX_DIR / "stage1_gate.json", gate)
    log(f"[gate] open={gate['open']} status={gate.get('status')} reason={gate['reason']}")
    if not gate["open"]:
        (GTEX_DIR / "STAGE1_NOT_RUN.txt").write_text(
            f"Stage 1 not run.\nstatus={gate.get('status')}\nreason={gate['reason']}\n",
            encoding="utf-8",
        )
        log("STAGE 1 not run (gate closed).")
        return None
    stale = GTEX_DIR / "STAGE1_NOT_RUN.txt"
    if stale.exists():
        stale.unlink()
    if not PREREG_FLAG.exists():
        raise StopStep(
            "prereg",
            f"pre-registration flag missing ({PREREG_FLAG.name}). "
            "Write FINDINGS_GTEX.md pre-registration before any T-cell fit.",
        )
    return gate


def permute_age_within_center(obs, rng, y=None):
    """Donor-level age permutation within SMCENTER. Age is donor-constant."""
    y = obs.age_mid.to_numpy(float).copy() if y is None else np.asarray(y, float).copy()
    donor = obs.donor.astype(str).to_numpy()
    center = obs.SMCENTER.astype(str).to_numpy()
    for c in pd.unique(center):
        mc = center == c
        d_u, inv = np.unique(donor[mc], return_inverse=True)
        age_d = np.empty(len(d_u), dtype=np.float64)
        y_c = y[mc]
        d_c = donor[mc]
        for i, d in enumerate(d_u):
            vals = y_c[d_c == d]
            if float(np.max(vals) - np.min(vals)) > 1e-9:
                raise StopStep(
                    "age_perm",
                    f"donor {d} has non-constant age_mid within center {c}",
                )
            age_d[i] = float(vals[0])
        rng.shuffle(age_d)
        y[mc] = age_d[inv]
    return y


def _fill_train_median(tr, te):
    tr = np.asarray(tr, float)
    te = np.asarray(te, float)
    med = float(np.nanmedian(tr)) if np.isfinite(tr).any() else 0.0
    n_tr = int((~np.isfinite(tr)).sum())
    n_te = int((~np.isfinite(te)).sum())
    tr = np.where(np.isfinite(tr), tr, med)
    te = np.where(np.isfinite(te), te, med)
    return tr, te, dict(n_filled_train=n_tr, n_filled_test=n_te, fill=med)


def build_resid_Z(obs, tr, te, comp_only=False):
    """Train-only one-hots. resid columns: SMRIN, SMTSISCH, DTHHRDY, SEX, neuronal, glial.
    Composition columns omitted on non-brain (all-NaN). Missing numeric filled from train median.
    Dummy NaNs (missing DTHHRDY/SEX) set to 0. Rows are not dropped.
    """
    dtr, dte = obs.iloc[tr], obs.iloc[te]
    if comp_only:
        if "neuronal_score" not in obs.columns or "glial_score" not in obs.columns:
            raise StopStep("resid", "composition scores not on obs (Stage 0 must write them)")
        ntr, nte, nf = _fill_train_median(dtr.neuronal_score, dte.neuronal_score)
        gtr, gte, gf = _fill_train_median(dtr.glial_score, dte.glial_score)
        Ztr = np.column_stack([ntr, gtr])
        Zte = np.column_stack([nte, gte])
        return Ztr, Zte, dict(comp_only=True, neuronal_fill=nf, glial_fill=gf)

    need = ["SMRIN", "SMTSISCH", "DTHHRDY", "SEX"]
    missing = [c for c in need if c not in obs.columns]
    if missing:
        raise StopStep("resid", f"obs missing {missing}", dict(columns=list(obs.columns)))
    rin_tr, rin_te, rin_f = _fill_train_median(dtr.SMRIN, dte.SMRIN)
    isch_tr, isch_te, isch_f = _fill_train_median(dtr.SMTSISCH, dte.SMTSISCH)
    dth_tr, lv_dth = onehot_train(dtr.DTHHRDY, drop_first=True)
    dth_te, unseen_dth = apply_onehot(dte.DTHHRDY, lv_dth)
    sex_tr, lv_sex = onehot_train(dtr.SEX, drop_first=True)
    sex_te, unseen_sex = apply_onehot(dte.SEX, lv_sex)
    dth_tr = np.nan_to_num(dth_tr, nan=0.0)
    dth_te = np.nan_to_num(dth_te, nan=0.0)
    sex_tr = np.nan_to_num(sex_tr, nan=0.0)
    sex_te = np.nan_to_num(sex_te, nan=0.0)
    parts_tr = [rin_tr, isch_tr, dth_tr, sex_tr]
    parts_te = [rin_te, isch_te, dth_te, sex_te]
    use_comp = False
    tissue = str(dtr.SMTSD.iloc[0]) if "SMTSD" in dtr.columns else ""
    if "neuronal_score" in obs.columns and "glial_score" in obs.columns:
        n_ok = int(np.isfinite(dtr.neuronal_score.to_numpy(float)).sum())
        g_ok = int(np.isfinite(dtr.glial_score.to_numpy(float)).sum())
        use_comp = n_ok >= 8 and g_ok >= 8
        if tissue in BRAIN_TISSUES and not use_comp:
            raise StopStep("resid", f"{tissue}: composition scores missing on train (n_ok={n_ok} g_ok={g_ok})")
        if use_comp:
            ntr, nte, nf = _fill_train_median(dtr.neuronal_score, dte.neuronal_score)
            gtr, gte, gf = _fill_train_median(dtr.glial_score, dte.glial_score)
            parts_tr += [ntr, gtr]
            parts_te += [nte, gte]
    Ztr = np.column_stack(parts_tr)
    Zte = np.column_stack(parts_te)
    return Ztr, Zte, dict(
        dth_levels=lv_dth, sex_levels=lv_sex, unseen_dth=unseen_dth, unseen_sex=unseen_sex,
        rin_fill=rin_f, isch_fill=isch_f, use_comp=use_comp,
    )


def prepare_fold_X(X, obs, tr, te, regime, method="ridge"):
    X = np.asarray(X, np.float64)
    y = obs.age_mid.to_numpy(float)
    Xtr, Xte = X[tr], X[te]
    ytr, yte = y[tr], y[te]
    meta = dict(regime=regime)
    if regime == "resid":
        Ztr, Zte, zmeta = build_resid_Z(obs, tr, te, comp_only=False)
        Xtr, Xte, rmeta = residualize_train_test(Xtr, Xte, Ztr, Zte)
        meta.update(zmeta)
        meta.update(rmeta)
    elif regime == "comp_only":
        Ztr, Zte, zmeta = build_resid_Z(obs, tr, te, comp_only=True)
        Xtr, Xte, rmeta = residualize_train_test(Xtr, Xte, Ztr, Zte)
        meta.update(zmeta)
        meta.update(rmeta)
    elif regime != "raw":
        raise StopStep("regime", f"unknown regime {regime!r}")
    Xtr_z, mu, sd = zscore_train(Xtr)
    Xte_z = (np.asarray(Xte, np.float64) - mu) / sd
    return Xtr_z, Xte_z, ytr, yte, meta


def ridge_predict_prestd(Xtr_z, ytr, Xte_z, svd=None):
    coef, ymu, alpha, svd_out = ridge_prestd(Xtr_z, ytr, svd=svd)
    pred = ymu + np.asarray(Xte_z, np.float64) @ coef
    w, _ = unit(coef)
    w = align_age_sign(w, Xtr_z, ytr)
    w, _ = unit(w)
    return pred, w, dict(ok=True, method="ridge", alpha=alpha, intercept=float(ymu)), svd_out


def predict_fold(Xtr_z, ytr, Xte_z, method, svd=None):
    if method == "ridge":
        return ridge_predict_prestd(Xtr_z, ytr, Xte_z, svd=svd)
    pred, w, _, cal = pls1_predict(Xtr_z, ytr, Xte_z)
    return pred, w, dict(ok=True, method="pls1", cal=cal), None


def load_tissue_pack(slug):
    path = GTEX_PROC / f"stage0_{slug}.npz"
    if not path.exists():
        raise StopStep("stage1_data", f"missing {path} — run Stage 0")
    z = np.load(path, allow_pickle=True)
    obs = pd.read_csv(GTEX_DIR / "stage0_obs.csv")
    obs = obs[obs.SMTSD.astype(str) == str(z["tissue"])].reset_index(drop=True)
    ids = z["sample_ids"].astype(str)
    pos = {s: i for i, s in enumerate(obs.SAMPID.astype(str))}
    missing = [s for s in ids if s not in pos]
    if missing:
        raise StopStep("stage1_data", f"{slug}: obs missing {len(missing)} sample ids")
    obs = obs.iloc[[pos[s] for s in ids]].reset_index(drop=True)
    if list(obs.SAMPID.astype(str)) != list(ids):
        raise StopStep("stage1_data", f"{slug}: obs/matrix SAMPID order mismatch")
    X = np.asarray(z["logcpm"], np.float64).T  # samples × genes
    return dict(X=X, obs=obs, genes=z["genes"], symbols=z["symbols"],
                ensembl=z["ensembl"].astype(str), tissue=str(z["tissue"]), slug=slug)


def build_oof_cache(pack, regime, method, log):
    """Prepare fold X (independent of y) and ridge Gram cache."""
    X, obs = pack["X"], pack["obs"]
    centers = obs.SMCENTER.astype(str).to_numpy()
    batches = obs.SMNABTCH.astype(str).to_numpy()
    donors = obs.donor.astype(str).to_numpy()
    folds = []
    for c in sorted(pd.unique(centers)):
        mc = np.flatnonzero(centers == c)
        if len(mc) < 10:
            log(f"   skip center {c}: n={len(mc)}")
            continue
        n_lv = int(pd.Series(batches[mc]).nunique())
        kf = grouped_kfold(batches[mc], n_splits=min(N_FOLDS, n_lv))
        for fi, (tr_l, te_l, te_lv) in enumerate(kf):
            tr, te = mc[tr_l], mc[te_l]
            assert_disjoint(donors[tr], donors[te], what="donor",
                            where=f"{pack['slug']} {regime} {c} fold{fi}")
            assert_disjoint(batches[tr], batches[te], what="SMNABTCH",
                            where=f"{pack['slug']} {regime} {c} fold{fi}")
            Xtr, Xte, _, _, meta = prepare_fold_X(X, obs, tr, te, regime, method)
            svd = None
            if method == "ridge":
                # Gram of X (U, S2) does not depend on y; dummy y only satisfies ridge_prestd's API.
                svd = ridge_prestd(Xtr, np.zeros(len(tr), dtype=float))[3]
            folds.append(dict(
                center=c, fold=fi, tr=tr, te=te, Xtr=Xtr, Xte=Xte, svd=svd,
                test_batches=",".join(te_lv), meta=meta,
            ))
    return folds


def oof_from_cache(pack, cache, y, method):
    obs = pack["obs"]
    y = np.asarray(y, float)
    pred = np.full(len(obs), np.nan)
    for f in cache:
        p, w, rec, _ = predict_fold(f["Xtr"], y[f["tr"]], f["Xte"], method, svd=f["svd"])
        pred[f["te"]] = p
    centers = obs.SMCENTER.astype(str).to_numpy()
    per_c = []
    used = sorted({f["center"] for f in cache})
    for c in used:
        m = centers == c
        sc = pred_scores(y[m], pred[m])
        sc["n"] = int(m.sum())
        per_c.append(dict(center=c, **sc))
    return dict(pred=pred, per_center=per_c)


def _median_stat(per_c, key):
    vals = [r[key] for r in per_c if np.isfinite(r.get(key, np.nan))]
    return float(np.median(vals)) if vals else np.nan


def bootstrap_angle(pack, regime, method, rng, log):
    X, obs = pack["X"], pack["obs"]
    y = obs.age_mid.to_numpy(float)
    ws = []
    for b in range(N_BOOT):
        idx = donor_bootstrap_index(obs, rng)
        dummy_te = idx[:1]
        Xtr_z, _, ytr, _, _ = prepare_fold_X(X, obs, idx, dummy_te, regime, method)
        w, rec = age_direction(Xtr_z, ytr, method=method)
        if rec.get("ok") and np.isfinite(w).all() and float(np.linalg.norm(w)) > 1e-12:
            ws.append(w / np.linalg.norm(w))
        if (b + 1) % 5 == 0:
            log(f"   bootstrap {b+1}/{N_BOOT} {pack['slug']} {regime} {method}")
    if len(ws) < 2:
        return dict(median_pairwise_deg=np.nan, n_ok=len(ws), n_boot=N_BOOT, n_pairs=0)
    Wc = np.column_stack(ws)
    ang = pairwise_angles_deg(Wc)
    iu = np.triu_indices(ang.shape[0], 1)
    vals = ang[iu]
    return dict(median_pairwise_deg=float(np.median(vals)), mean_pairwise_deg=float(np.mean(vals)),
                p05=float(np.percentile(vals, 5)), p95=float(np.percentile(vals, 95)),
                n_ok=len(ws), n_boot=N_BOOT, n_pairs=int(len(vals)))


def t1(packs, log, tissues=None, regimes=None, tag="t1"):
    rng = np.random.default_rng(CELL_SEEDS["T1"] if tag == "t1" else CELL_SEEDS["COMP"])
    tissues = tissues or TISSUES
    regimes = regimes or REGIMES
    rows = []
    t0 = time.time()
    for t in tissues:
        slug = TISSUE_SLUG[t]
        pack = packs[slug]
        for regime in regimes:
            for method in METHODS:
                log(f"[T1] {t} {regime} {method}")
                boot = bootstrap_angle(pack, regime, method, rng, log)
                cache = build_oof_cache(pack, regime, method, log)
                y0 = pack["obs"].age_mid.to_numpy(float)
                oof = oof_from_cache(pack, cache, y0, method)
                null_r, null_r2, null_rho, null_cal = [], [], [], []
                for p in range(N_PERM_GTEX):
                    yp = permute_age_within_center(pack["obs"], rng)
                    oof_p = oof_from_cache(pack, cache, yp, method)
                    null_r.append(_median_stat(oof_p["per_center"], "r"))
                    null_r2.append(_median_stat(oof_p["per_center"], "r2"))
                    null_rho.append(_median_stat(oof_p["per_center"], "rho"))
                    null_cal.append(_median_stat(oof_p["per_center"], "cal_r2"))
                    if (p + 1) % 50 == 0:
                        log(f"   perm {p+1}/{N_PERM_GTEX} ({time.time()-t0:.0f}s)")
                r_med = _median_stat(oof["per_center"], "r")
                r2_med = _median_stat(oof["per_center"], "r2")
                rho_med = _median_stat(oof["per_center"], "rho")
                cal_med = _median_stat(oof["per_center"], "cal_r2")
                r_null = float(np.nanmedian(null_r))
                rows.append(dict(
                    test="T1", tissue=t, regime=regime, method=method,
                    boot_median_angle=boot["median_pairwise_deg"],
                    boot_p05=boot.get("p05"), boot_p95=boot.get("p95"),
                    boot_n_ok=boot["n_ok"],
                    within_batch_r_median_over_centers=r_med,
                    within_batch_r_null=r_null,
                    within_batch_r_p=permutation_p(r_med, np.asarray(null_r, float), greater=True),
                    within_batch_r2_median_over_centers=r2_med,
                    within_batch_r2_null=float(np.nanmedian(null_r2)),
                    within_batch_r2_p=permutation_p(r2_med, np.asarray(null_r2, float), greater=True),
                    within_batch_rho=rho_med,
                    within_batch_rho_null=float(np.nanmedian(null_rho)),
                    within_batch_cal_r2=cal_med,
                    within_batch_cal_r2_null=float(np.nanmedian(null_cal)),
                    pass_r=pass_r_bar(r_med, r_null),
                    n_perm=N_PERM_GTEX,
                    sex_angle_ref=SEX_ANGLE_REF,
                    planted_angle_n233=PLANTED_ANGLE_N233,
                    n_centers=len(oof["per_center"]),
                ))
    df = pd.DataFrame(rows)
    df.to_csv(GTEX_DIR / f"{tag}.csv", index=False)
    dump_json(GTEX_DIR / f"{tag}.json", rows)
    log(f"[{tag}] wrote {len(rows)} rows in {time.time()-t0:.0f}s")
    return rows


def _transfer_align(pack_tr, pack_te, donor_hold="disjoint"):
    e_tr = pack_tr["ensembl"].astype(str)
    e_te = pack_te["ensembl"].astype(str)
    common = [e for e in e_tr if e in set(e_te)]
    if len(common) < 100:
        raise StopStep("transfer", f"common genes {len(common)}")
    i_tr = {e: i for i, e in enumerate(e_tr)}
    i_te = {e: i for i, e in enumerate(e_te)}
    idx_tr = np.array([i_tr[e] for e in common], dtype=int)
    idx_te = np.array([i_te[e] for e in common], dtype=int)
    otr, ote = pack_tr["obs"], pack_te["obs"]
    d_tr = otr.donor.astype(str).to_numpy()
    d_te = ote.donor.astype(str).to_numpy()
    tr = np.arange(len(otr))
    te = np.arange(len(ote))
    n_overlap = int(len(set(d_tr) & set(d_te)))
    if donor_hold == "disjoint":
        overlap = set(d_tr) & set(d_te)
        tr = np.flatnonzero(~np.isin(d_tr, list(overlap)))
        assert_disjoint(d_tr[tr], d_te[te], what="donor",
                        where=f"{pack_tr['slug']}->{pack_te['slug']}")
    return dict(idx_tr=idx_tr, idx_te=idx_te, tr=tr, te=te, n_genes=len(common),
                n_overlap_donors_dropped=n_overlap)


def transfer_prepare(pack_tr, pack_te, regime, method, align):
    Xtr = pack_tr["X"][np.ix_(align["tr"], align["idx_tr"])]
    Xte = pack_te["X"][np.ix_(align["te"], align["idx_te"])]
    dummy_tr = pack_tr["obs"].iloc[align["tr"]].reset_index(drop=True)
    dummy_te = pack_te["obs"].iloc[align["te"]].reset_index(drop=True)
    ntr, nte = len(dummy_tr), len(dummy_te)
    Xcat = np.vstack([Xtr, Xte])
    obst = pd.concat([dummy_tr, dummy_te], ignore_index=True)
    tr_i = np.arange(ntr)
    te_i = np.arange(ntr, ntr + nte)
    Xtr_z, Xte_z, ytr, yte, meta = prepare_fold_X(Xcat, obst, tr_i, te_i, regime, method)
    svd = None
    if method == "ridge":
        svd = ridge_prestd(Xtr_z, ytr)[3]
    return dict(Xtr_z=Xtr_z, Xte_z=Xte_z, ytr=ytr, yte=yte, svd=svd, meta=meta,
                n_train=ntr, n_test=nte, obs_tr=dummy_tr)


def transfer_score(prep, method, ytr=None, yte=None):
    ytr = prep["ytr"] if ytr is None else np.asarray(ytr, float)
    yte = prep["yte"] if yte is None else np.asarray(yte, float)
    pred, w, rec, _ = predict_fold(prep["Xtr_z"], ytr, prep["Xte_z"], method, svd=prep["svd"])
    sc = pred_scores(yte, pred)
    sc.update(dict(w=w, pred=pred, rec=rec))
    return sc


def run_transfer_pair(src, dst, regime, method, rng, log, label, test_id):
    align = _transfer_align(src, dst)
    prep = transfer_prepare(src, dst, regime, method, align)
    obs = transfer_score(prep, method)
    nulls = {k: [] for k in ("r", "r2", "rho", "cal_r2")}
    for p in range(N_PERM_GTEX):
        # train-bank donor permutation within center; test y held fixed (S2)
        ytr_p = permute_age_within_center(prep["obs_tr"], rng)
        rec_p = transfer_score(prep, method, ytr=ytr_p, yte=prep["yte"])
        for k in nulls:
            nulls[k].append(rec_p[k])
        if (p + 1) % 50 == 0:
            log(f"   perm {p+1}/{N_PERM_GTEX} {label} {regime} {method}")
    r_null = float(np.nanmedian(nulls["r"]))
    return dict(
        test=test_id, label=label, train_tissue=src["tissue"], test_tissue=dst["tissue"],
        regime=regime, method=method,
        r=obs["r"], r_null=r_null, r_p=permutation_p(obs["r"], np.asarray(nulls["r"], float), greater=True),
        r2=obs["r2"], r2_null=float(np.nanmedian(nulls["r2"])),
        cal_r2=obs["cal_r2"], cal_r2_null=float(np.nanmedian(nulls["cal_r2"])),
        rho=obs["rho"], rho_null=float(np.nanmedian(nulls["rho"])),
        pass_r=pass_r_bar(obs["r"], r_null),
        n_train=prep["n_train"], n_test=prep["n_test"],
        n_genes=align["n_genes"], n_overlap_donors_dropped=align["n_overlap_donors_dropped"],
        n_perm=N_PERM_GTEX,
    )


def t2(packs, log):
    flags = load_json(GTEX_DIR / "stage0_summary.json").get("stop_flags", {})
    if flags.get("t2_not_testable"):
        rec = dict(skipped=True, reason="t2_not_testable: fewer than 3 SMCENTER with ≥25 BA9 samples",
                   n_SMCENTER_ge25=flags.get("n_SMCENTER_ge25"),
                   center_counts=flags.get("center_counts"))
        dump_json(GTEX_DIR / "t2.json", rec)
        pd.DataFrame([rec]).to_csv(GTEX_DIR / "t2.csv", index=False)
        log("[T2] SKIP: fewer than 3 SMCENTER with ≥25 BA9 samples")
        return rec
    rng = np.random.default_rng(CELL_SEEDS["T2"])
    pack = packs[TISSUE_SLUG[PRIMARY_TISSUE]]
    obs = pack["obs"]
    centers = sorted(obs.SMCENTER.astype(str).unique().tolist())
    X, donors = pack["X"], obs.donor.astype(str).to_numpy()
    cvec = obs.SMCENTER.astype(str).to_numpy()
    rows = []
    for regime in REGIMES:
        for method in METHODS:
            for held in centers:
                te = np.flatnonzero(cvec == held)
                tr = np.flatnonzero(cvec != held)
                if len(tr) < 8 or len(te) < 8:
                    log(f"[T2] skip LOCO {held} n_tr={len(tr)} n_te={len(te)}")
                    continue
                assert_disjoint(donors[tr], donors[te], what="donor", where=f"LOCO {held} {regime}")
                log(f"[T2] LOCO hold {held} {regime} {method} n_tr={len(tr)} n_te={len(te)}")
                dummy_src = dict(pack)
                dummy_src["X"] = X
                dummy_src["obs"] = obs
                # Build a fake two-pack split by index via prepare_fold_X
                Xtr, Xte, ytr, yte, _ = prepare_fold_X(X, obs, tr, te, regime, method)
                svd = ridge_prestd(Xtr, ytr)[3] if method == "ridge" else None
                pred, w, rec, _ = predict_fold(Xtr, ytr, Xte, method, svd=svd)
                sc = pred_scores(yte, pred)
                nulls = {k: [] for k in ("r", "r2", "rho", "cal_r2")}
                obs_tr = obs.iloc[tr].reset_index(drop=True)
                for p in range(N_PERM_GTEX):
                    ytr_p = permute_age_within_center(obs_tr, rng)
                    pred_p, _, _, _ = predict_fold(Xtr, ytr_p, Xte, method, svd=svd)
                    scp = pred_scores(yte, pred_p)
                    for k in nulls:
                        nulls[k].append(scp[k])
                r_null = float(np.nanmedian(nulls["r"]))
                rows.append(dict(
                    test="T2_LOCO", tissue=PRIMARY_TISSUE, hold_center=held,
                    train_centers=",".join([c for c in centers if c != held]),
                    regime=regime, method=method,
                    r=sc["r"], r_null=r_null,
                    r_p=permutation_p(sc["r"], np.asarray(nulls["r"], float), greater=True),
                    r2=sc["r2"], r2_null=float(np.nanmedian(nulls["r2"])),
                    cal_r2=sc["cal_r2"], cal_r2_null=float(np.nanmedian(nulls["cal_r2"])),
                    rho=sc["rho"], rho_null=float(np.nanmedian(nulls["rho"])),
                    pass_r=pass_r_bar(sc["r"], r_null),
                    n_train=int(len(tr)), n_test=int(len(te)), n_perm=N_PERM_GTEX,
                ))
            for a in centers:
                for b in centers:
                    if a == b:
                        continue
                    tr = np.flatnonzero(cvec == a)
                    te = np.flatnonzero(cvec == b)
                    if len(tr) < 8 or len(te) < 8:
                        continue
                    assert_disjoint(donors[tr], donors[te], what="donor",
                                    where=f"pair {a}->{b} {regime}")
                    log(f"[T2] pair {a}->{b} {regime} {method}")
                    Xtr, Xte, ytr, yte, _ = prepare_fold_X(X, obs, tr, te, regime, method)
                    svd = ridge_prestd(Xtr, ytr)[3] if method == "ridge" else None
                    pred, w, rec, _ = predict_fold(Xtr, ytr, Xte, method, svd=svd)
                    sc = pred_scores(yte, pred)
                    nulls = {k: [] for k in ("r", "r2", "rho", "cal_r2")}
                    obs_tr = obs.iloc[tr].reset_index(drop=True)
                    for p in range(N_PERM_GTEX):
                        ytr_p = permute_age_within_center(obs_tr, rng)
                        pred_p, _, _, _ = predict_fold(Xtr, ytr_p, Xte, method, svd=svd)
                        scp = pred_scores(yte, pred_p)
                        for k in nulls:
                            nulls[k].append(scp[k])
                    r_null = float(np.nanmedian(nulls["r"]))
                    rows.append(dict(
                        test="T2_pair", tissue=PRIMARY_TISSUE,
                        train_center=a, test_center=b,
                        regime=regime, method=method,
                        r=sc["r"], r_null=r_null,
                        r_p=permutation_p(sc["r"], np.asarray(nulls["r"], float), greater=True),
                        r2=sc["r2"], r2_null=float(np.nanmedian(nulls["r2"])),
                        cal_r2=sc["cal_r2"], cal_r2_null=float(np.nanmedian(nulls["cal_r2"])),
                        rho=sc["rho"], rho_null=float(np.nanmedian(nulls["rho"])),
                        pass_r=pass_r_bar(sc["r"], r_null),
                        n_train=int(len(tr)), n_test=int(len(te)), n_perm=N_PERM_GTEX,
                    ))
    pd.DataFrame(rows).to_csv(GTEX_DIR / "t2.csv", index=False)
    dump_json(GTEX_DIR / "t2.json", rows)
    return rows


def t3(packs, log):
    rng = np.random.default_rng(CELL_SEEDS["T3"])
    ba9 = packs[TISSUE_SLUG[PRIMARY_TISSUE]]
    other = packs[TISSUE_SLUG[PAIRED_TISSUE]]
    rows = []
    directions = [
        (ba9, other, f"{PRIMARY_TISSUE} -> {PAIRED_TISSUE}"),
        (other, ba9, f"{PAIRED_TISSUE} -> {PRIMARY_TISSUE}"),
    ]
    for src, dst, lab in directions:
        for regime in REGIMES:
            for method in METHODS:
                log(f"[T3] {lab} {regime} {method}")
                rows.append(run_transfer_pair(src, dst, regime, method, rng, log, lab, "T3"))
    pd.DataFrame(rows).to_csv(GTEX_DIR / "t3.csv", index=False)
    dump_json(GTEX_DIR / "t3.json", rows)
    _write_t34()
    return rows


def t4(packs, log):
    rng = np.random.default_rng(CELL_SEEDS["T4"])
    ba9 = packs[TISSUE_SLUG[PRIMARY_TISSUE]]
    rows = []
    for tgt in CROSS_REGION:
        other = packs[TISSUE_SLUG[tgt]]
        lab = f"{PRIMARY_TISSUE} -> {tgt}"
        for regime in REGIMES:
            for method in METHODS:
                log(f"[T4] {lab} {regime} {method}")
                rows.append(run_transfer_pair(ba9, other, regime, method, rng, log, lab, "T4"))
    pd.DataFrame(rows).to_csv(GTEX_DIR / "t4.csv", index=False)
    dump_json(GTEX_DIR / "t4.json", rows)
    _write_t34()
    return rows


def _write_t34():
    frames = []
    for name in ("t3.csv", "t4.csv"):
        p = GTEX_DIR / name
        if p.exists():
            frames.append(pd.read_csv(p))
    if frames:
        pd.concat(frames, ignore_index=True).to_csv(GTEX_DIR / "t3_t4.csv", index=False)


def t5(packs, log):
    rng = np.random.default_rng(CELL_SEEDS["T5"])
    sets = [set(packs[TISSUE_SLUG[t]]["ensembl"].astype(str)) for t in TISSUES]
    common = sets[0]
    for s in sets[1:]:
        common &= s
    common = sorted(common)
    log(f"[T5] common filtered Ensembl across 6 tissues: {len(common):,}")
    if len(common) < 100:
        raise StopStep("T5", f"common gene space too small: {len(common)}")
    mats, obs_list = [], []
    for t in TISSUES:
        pack = packs[TISSUE_SLUG[t]]
        idx = {e: i for i, e in enumerate(pack["ensembl"].astype(str))}
        take = np.array([idx[e] for e in common], dtype=int)
        mats.append(pack["X"][:, take])
        obs_list.append(pack["obs"])
    Xall = np.vstack(mats)
    mu, sd = Xall.mean(0), Xall.std(0)
    sd = np.where(sd < 1e-12, 1.0, sd)
    cents = [((M - mu) / sd).mean(0) for M in mats]
    C = np.vstack(cents)
    C = C - C.mean(0, keepdims=True)
    _, S, Vt = np.linalg.svd(C, full_matrices=False)
    rank = int(np.sum(S > 1e-8))
    B = Vt[:rank].T
    log(f"[T5] identity subspace rank={rank} (expect 5)")
    ba9 = packs[TISSUE_SLUG[PRIMARY_TISSUE]]
    idx = {e: i for i, e in enumerate(ba9["ensembl"].astype(str))}
    take = np.array([idx[e] for e in common], dtype=int)
    Xb = (ba9["X"][:, take] - mu) / sd
    y = ba9["obs"].age_mid.to_numpy(float)
    w, rec = age_direction(Xb, y, method="pls1")
    ang = angle_to_subspace_deg(w, B)
    pr = proj_r2(w, B)
    W = []
    for M, o in zip(mats, obs_list):
        Z = (M - mu) / sd
        ww, _ = age_direction(Z, o.age_mid.to_numpy(float), method="pls1")
        W.append(ww)
    units = [u / np.linalg.norm(u) for u in W if np.isfinite(u).all() and np.linalg.norm(u) > 1e-12]
    pang = pairwise_angles_deg(np.column_stack(units))
    iu = np.triu_indices(pang.shape[0], 1)
    pair_med = float(np.median(pang[iu]))
    null_ang, null_pair = [], []
    for p in range(N_PERM_GTEX):
        yp = permute_age_within_center(ba9["obs"], rng)
        wp, _ = age_direction(Xb, yp, method="pls1")
        null_ang.append(angle_to_subspace_deg(wp, B))
        Wp = []
        for M, o in zip(mats, obs_list):
            Z = (M - mu) / sd
            yy = permute_age_within_center(o, rng)
            ww, _ = age_direction(Z, yy, method="pls1")
            Wp.append(ww)
        pa = pairwise_angles_deg(np.column_stack([u / max(np.linalg.norm(u), 1e-12) for u in Wp]))
        null_pair.append(float(np.median(pa[np.triu_indices(pa.shape[0], 1)])))
        if (p + 1) % 50 == 0:
            log(f"   T5 perm {p+1}/{N_PERM_GTEX}")
    out = dict(
        n_common_genes=len(common), identity_rank=rank, identity_sv=S[:rank].tolist(),
        ba9_angle_to_identity=float(ang), ba9_proj_r2=float(pr),
        null_angle_mean=float(np.nanmean(null_ang)),
        null_angle_median=float(np.nanmedian(null_ang)),
        p_angle=permutation_p(ang, np.asarray(null_ang, float), greater=False),
        pairwise_median_deg=pair_med,
        pairwise_null_mean=float(np.nanmean(null_pair)),
        p_pairwise=permutation_p(pair_med, np.asarray(null_pair, float), greater=False),
        tissues=list(TISSUES), n_perm=N_PERM_GTEX,
        note="bulk analogue of FINDINGS_GEOMETRY Part C, not a replication",
    )
    dump_json(GTEX_DIR / "t5.json", out)
    return out


def t6(packs, log):
    rng = np.random.default_rng(CELL_SEEDS["T6"])
    ba9 = packs[TISSUE_SLUG[PRIMARY_TISSUE]]
    gtex_ens = set(ba9["ensembl"].astype(str))
    from gtex_common import load_dlpfc_gene_universe
    universe, gdf = load_dlpfc_gene_universe()
    overlap = sorted(gtex_ens & universe)
    log(f"[T6] BA9 filtered ∩ DLPFC universe = {len(overlap):,}")
    if len(overlap) < T6_MIN_OVERLAP:
        rec = dict(ran=False, n_overlap=len(overlap), min_required=T6_MIN_OVERLAP,
                   reason="overlap < 15000")
        dump_json(GTEX_DIR / "t6.json", rec)
        return rec
    axes = TARGET_DIR / "v1_axes.npz"
    tw = TARGET_DIR / "target_W.npz"
    path = axes if axes.exists() else tw
    if not path.exists():
        raise StopStep("T6", f"neither {axes.name} nor {tw.name} in results/target/")
    z = np.load(path, allow_pickle=True)
    log(f"[T6] loaded {path.name} keys={list(z.files)}")
    w_key = None
    for key in ("W_age_ridge", "W_age", "W_ridge"):
        if key in z.files:
            w_key = key
            W = np.asarray(z[key], float)
            break
    if w_key is None:
        raise StopStep("T6", f"{path.name} has no W_age_ridge/W_age/W_ridge. keys={list(z.files)}")
    log(f"[T6] using {w_key} shape={W.shape}")
    cols = []
    for j in range(W.shape[1]):
        w = W[:, j]
        if np.isfinite(w).all() and float(np.linalg.norm(w)) > 1e-12:
            cols.append(w / np.linalg.norm(w))
    if not cols:
        raise StopStep("T6", "no finite per-type DLPFC age directions")
    w_dlpfc = np.mean(np.column_stack(cols), 1)
    w_dlpfc = w_dlpfc / np.linalg.norm(w_dlpfc)
    if "gene_id" in z.files:
        dlpfc_ids = np.array([strip_ensembl(x) for x in np.asarray(z["gene_id"]).astype(str)])
    else:
        dlpfc_ids = gdf.gene_id.astype(str).map(strip_ensembl).to_numpy()
    if W.shape[0] != len(dlpfc_ids):
        raise StopStep("T6", f"W rows {W.shape[0]} != DLPFC genes {len(dlpfc_ids)}")
    pos_d = {e: i for i, e in enumerate(dlpfc_ids)}
    pos_g = {e: i for i, e in enumerate(ba9["ensembl"].astype(str))}
    take_d = np.array([pos_d[e] for e in overlap], dtype=int)
    take_g = np.array([pos_g[e] for e in overlap], dtype=int)
    wd = w_dlpfc[take_d]
    X = ba9["X"][:, take_g]
    Xz, _, _ = zscore_train(X)
    wg, rec = age_direction(Xz, ba9["obs"].age_mid.to_numpy(float), method="ridge")
    ang = unsigned_angle_deg(wg, wd)
    pdim = len(overlap)
    null = []
    n_rand = max(N_PERM_GTEX * 5, 100)
    for i in range(n_rand):
        v = rng.normal(size=pdim)
        v = v / np.linalg.norm(v)
        null.append(unsigned_angle_deg(wg, v))
    out = dict(
        ran=True, n_overlap=len(overlap), angle_deg=float(ang),
        null_mean=float(np.mean(null)), null_median=float(np.median(null)),
        p=permutation_p(ang, np.asarray(null, float), greater=False),
        source=path.name, w_key=w_key, n_types_collapsed=len(cols), n_random=n_rand,
        note="one number with a null; not evidence beyond closer-than-random",
    )
    dump_json(GTEX_DIR / "t6.json", out)
    return out


def composition_check(packs, log):
    """If BA9 raw ridge T3 both directions pass, rerun T1 with composition proxies alone."""
    t3p = GTEX_DIR / "t3.csv"
    if not t3p.exists():
        rec = dict(ran=False, reason="t3.csv missing")
        dump_json(GTEX_DIR / "t1_comp.json", rec)
        return rec
    t3 = pd.read_csv(t3p)
    sub = t3[(t3.regime == "raw") & (t3.method == "ridge")]
    if len(sub) < 2 or not bool(sub.pass_r.all()):
        rec = dict(ran=False, reason="BA9 raw ridge T3 did not pass both directions; composition-only not triggered")
        dump_json(GTEX_DIR / "t1_comp.json", rec)
        pd.DataFrame([rec]).to_csv(GTEX_DIR / "t1_comp.csv", index=False)
        log("[COMP] not triggered")
        return rec
    log("[COMP] BA9 raw ridge T3 passed both ways; rerunning T1 BA9 with composition proxies alone")
    rows = t1(packs, log, tissues=(PRIMARY_TISSUE,), regimes=("comp_only",), tag="t1_comp")
    return rows


def _t3_pass(regime, method="ridge"):
    p = GTEX_DIR / "t3.csv"
    if not p.exists():
        return None
    df = pd.read_csv(p)
    sub = df[(df.regime == regime) & (df.method == method)]
    if len(sub) < 2:
        return None
    return bool(sub.pass_r.all())


def write_verdict(log):
    t2 = load_json(GTEX_DIR / "t2.json") if (GTEX_DIR / "t2.json").exists() else {}
    t2_skipped = isinstance(t2, dict) and t2.get("skipped")
    raw_pass = _t3_pass("raw", "ridge")
    resid_pass = _t3_pass("resid", "ridge")
    raw_pls = _t3_pass("raw", "pls1")
    resid_pls = _t3_pass("resid", "pls1")
    # angles
    d_ang = np.nan
    if (GTEX_DIR / "t1.csv").exists():
        t1 = pd.read_csv(GTEX_DIR / "t1.csv")
        ba9 = t1[(t1.tissue == PRIMARY_TISSUE) & (t1.method == "ridge")]
        raw_a = ba9[ba9.regime == "raw"]
        res_a = ba9[ba9.regime == "resid"]
        if len(raw_a) and len(res_a):
            d_ang = float(res_a.boot_median_angle.iloc[0] - raw_a.boot_median_angle.iloc[0])
    # r signs T3 ridge
    sign_unchanged = None
    if (GTEX_DIR / "t3.csv").exists():
        t3 = pd.read_csv(GTEX_DIR / "t3.csv")
        rr = t3[(t3.method == "ridge")]
        raw = rr[rr.regime == "raw"].sort_values("label")
        resid = rr[rr.regime == "resid"].sort_values("label")
        if len(raw) == len(resid) and len(raw):
            sign_unchanged = bool(np.all(np.sign(raw.r.to_numpy()) == np.sign(resid.r.to_numpy())))
    ba9_n = (load_json(GTEX_DIR / "stage0_summary.json").get("stop_flags") or {}).get("ba9_n")
    planted_n = 150  # C3 planted dense_shared transferred as r 10/10 at n=150
    lines = []
    if t2_skipped:
        lines.append("T2 not testable (fewer than 3 SMCENTER with ≥25 BA9 samples). BA9 transfer is T3 (BA9 ↔ Cortex).")
    if raw_pls is not None and raw_pls != raw_pass:
        lines.append(f"ridge vs pls1 T3 raw disagree (ridge pass={raw_pass}, pls1 pass={raw_pls}); do not average.")
    material = bool(sign_unchanged) and (np.isfinite(d_ang) and abs(d_ang) < 5)
    if raw_pass and resid_pass and material:
        fired = (
            "BA9 raw passes transfer (held-out Pearson r > 0 with permutation null ≤ 0.05) AND resid does not change it "
            f"materially (sign of r unchanged, Δangle={d_ang:+.2f}° < 5°). A stable bulk cortex age direction exists and "
            "is not driven by recorded logistics; the DLPFC instability is then attributable to something GTEx does not "
            "share: two-bank structure, snRNA/pseudobulk noise, or unrecorded covariates. Say all three; do not pick one."
        )
    elif (not raw_pass) and resid_pass:
        fired = (
            "BA9 raw fails, resid passes. Recorded logistics (RIN / ischemic time / Hardy / composition) were hiding a "
            "stable direction. This makes unrecorded PMI/RIN a plausible, not proven, cause of the DLPFC failure."
        )
    elif (not raw_pass) and (not resid_pass):
        if ba9_n is not None and int(ba9_n) >= planted_n:
            fired = (
                f"Both fail at BA9 n={ba9_n} ≥ the n where POSCTRL's planted direction transferred as r (n=150). "
                "Cortex age directions fail transfer in a second, differently-structured cohort. The negative generalizes."
            )
        else:
            fired = (
                f"Both fail at BA9 n={ba9_n} below the n where POSCTRL's planted direction transferred as r. Uninformative."
            )
    elif raw_pass and resid_pass and not material:
        fired = (
            f"BA9 raw and resid both pass transfer as r, but resid changes it materially "
            f"(sign_unchanged={sign_unchanged}, Δangle={d_ang:+.2f}°). Report both; do not pick resid as primary."
        )
    elif raw_pass and not resid_pass:
        fired = (
            "BA9 raw passes transfer as r and resid does not. Recorded logistics are not required for the direction; "
            "residualizing them removes it. Report both; resid is not selected because it looks worse."
        )
    else:
        fired = "T3 pass pattern did not match a single pre-registered bullet. See T3 table. Do not average."
    if lines:
        fired = " ".join(lines) + " " + fired
    # T3 fail T2 pass
    if (not t2_skipped) and isinstance(t2, list) and resid_pass is not None:
        t2_df = pd.DataFrame(t2)
        t2_raw = t2_df[(t2_df.get("regime") == "raw") & (t2_df.get("method") == "ridge")] if len(t2_df) else t2_df
        t2_ok = bool(t2_raw["pass_r"].all()) if len(t2_raw) and "pass_r" in t2_raw.columns else False
        if (not raw_pass) and t2_ok:
            fired += " T3 fails while T2 passes: the preservation pathway defeats transfer even when collection center does not."
    # composition
    if (GTEX_DIR / "t1_comp.csv").exists():
        tc = pd.read_csv(GTEX_DIR / "t1_comp.csv")
        if "pass_r" in tc.columns and "reason" not in tc.columns:
            if not bool(tc.pass_r.any()):
                fired += (
                    " Composition-only T1: passing disappears. The bulk cortex ages by its mix, parallel to the blood "
                    "finding — not a stable cellular age direction."
                )
            else:
                fired += " Composition-only T1 still passes; the direction is not only the neuronal/glial mix."
    text = fired
    (GTEX_DIR / "stage1_verdict.txt").write_text(text + "\n", encoding="utf-8")
    log("[verdict] " + text)
    return text


def load_all_packs(log):
    packs = {}
    for t in TISSUES:
        slug = TISSUE_SLUG[t]
        log(f"[load] {slug}")
        packs[slug] = load_tissue_pack(slug)
    return packs


def run_cell(name, packs, log):
    t0 = time.time()
    log(f"===== {name} start =====")
    if name == "T1":
        t1(packs, log)
        progress_snapshot("T2 (or SKIP if t2_not_testable)", extra=f"T1 done in {time.time()-t0:.0f}s")
    elif name == "T2":
        t2(packs, log)
        progress_snapshot("T3 BA9 ↔ Cortex", extra=f"T2 done in {time.time()-t0:.0f}s")
    elif name == "T3":
        t3(packs, log)
        progress_snapshot("T4 BA9 → BA24 / hippocampus", extra=f"T3 done in {time.time()-t0:.0f}s")
    elif name == "T4":
        t4(packs, log)
        progress_snapshot("T5 identity subspace", extra=f"T4 done in {time.time()-t0:.0f}s")
    elif name == "T5":
        t5(packs, log)
        progress_snapshot("T6 DLPFC overlap angle", extra=f"T5 done in {time.time()-t0:.0f}s")
    elif name == "T6":
        t6(packs, log)
        progress_snapshot("COMP if T3 raw ridge passed; then findings", extra=f"T6 done in {time.time()-t0:.0f}s")
    elif name == "COMP":
        composition_check(packs, log)
        write_verdict(log)
        dump_json(GTEX_DIR / "stage1_summary.json", dict(
            seed=GTEX_SEED, n_perm=N_PERM_GTEX, prereg=PREREG_DATE,
            cells=["T1", "T2", "T3", "T4", "T5", "T6", "COMP"],
        ))
        progress_snapshot("write FINDINGS_GTEX.md", stop="Stage 1 cells complete",
                          extra=f"COMP/verdict done in {time.time()-t0:.0f}s")
    else:
        raise StopStep("cell", f"unknown cell {name}")
    log(f"===== {name} done {time.time()-t0:.0f}s =====")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", default="all",
                    help="T1 T2 T3 T4 T5 T6 COMP all")
    args = ap.parse_args()
    log = Logger(GTEX_DIR / "stage1_report.txt")
    gtex_log_banner(log, "STAGE 1")
    try:
        gate = _gate_or_stop(log)
        if gate is None:
            progress_snapshot("none (gate closed)", stop="Stage 1 gate closed")
            return
        log(f"[prereg] {PREREG_FLAG} exists; n_perm={N_PERM_GTEX}")
        packs = load_all_packs(log)
        order = ["T1", "T2", "T3", "T4", "T5", "T6", "COMP"]
        cells = order if args.cell == "all" else [args.cell]
        for c in cells:
            run_cell(c, packs, log)
        if args.cell == "all":
            log("STAGE 1 complete.")
    except StopStep as e:
        log(f"STOP [{e.step}] {e.message}")
        dump_json(GTEX_DIR / "stage1_STOP.json", dict(step=e.step, message=e.message, details=e.details))
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise
    except Exception as e:
        log(f"FAIL {type(e).__name__}: {e}")
        progress_snapshot(f"fix {type(e).__name__}", stop=f"FAIL {type(e).__name__}: {e}")
        raise
    finally:
        log.close()


if __name__ == "__main__":
    main()

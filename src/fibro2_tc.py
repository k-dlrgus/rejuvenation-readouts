"""T-C — HistGradientBoosting diagnostic on GTEx cultured fibroblasts. Not a Stage 2 replacement."""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fibro2_common import (  # noqa: E402
    FIBRO2_DIR, FIBRO2_SEED, FIBRO2_BOOT, N_PERM, N_BOOT, HGB_GRID,
    PREREG_TC, PREREG_TC_FLAG, TRANSFER_NULL_BAR,
    StopStep, Logger, dump_json, jsonable, fibro2_log_banner,
    load_manifest, save_manifest, record_failure, progress_snapshot, pass_rho_bar,
)
from fibro_common import (  # noqa: E402
    FIBRO_DIR, TRANSFER_SITES, EXCLUDE_SITE, EXCLUDE_SITE_REASON, SITE_CV_MIN_N,
    TRANSFER_MIN_N, N_FOLDS, bootstrap_rho_ci,
)
from fibro_stage1 import load_pack, _site_folds  # noqa: E402
from gtex_common import spearman_safe, pearson_safe, pred_scores, grouped_kfold, assert_disjoint  # noqa: E402
from gtex_stage1 import permute_age_within_center, prepare_fold_X  # noqa: E402
from target_common import zscore_train  # noqa: E402
from trajectory_common import permutation_p  # noqa: E402


def _hgb(**params):
    kw = dict(random_state=int(FIBRO2_SEED), early_stopping=False)
    kw.update(params)
    return HistGradientBoostingRegressor(**kw)


def _select_params(Xtr, ytr, groups_tr, log, tag):
    """Inner grouped CV on training rows only. Spearman ρ; freeze the winner for this split."""
    groups_tr = np.asarray(groups_tr).astype(str)
    n_lv = int(pd.Series(groups_tr).nunique())
    if n_lv < 2:
        raise StopStep("hgb_inner", f"{tag}: inner CV needs ≥2 SMNABTCH levels, got {n_lv}")
    n_splits = int(min(3, n_lv))
    folds = grouped_kfold(groups_tr, n_splits=n_splits)
    best, best_score, rows = None, -np.inf, []
    Xtr = np.asarray(Xtr, np.float32)
    ytr = np.asarray(ytr, np.float64)
    for params in HGB_GRID:
        rhos = []
        t0 = time.time()
        for tr, te, te_lv in folds:
            est = _hgb(**params)
            est.fit(Xtr[tr], ytr[tr])
            pred = est.predict(Xtr[te])
            r = spearman_safe(ytr[te], pred)
            if np.isfinite(r):
                rhos.append(float(r))
        med = float(np.nanmedian(rhos)) if rhos else np.nan
        rows.append(dict(tag=tag, **params, inner_rho_median=med, n_inner_folds=len(folds)))
        score = med if np.isfinite(med) else -np.inf
        if score > best_score:
            best_score, best = score, dict(params)
        log(f"   inner {tag} {params} ρ_med={med if np.isfinite(med) else 'NA'} ({time.time()-t0:.0f}s)")
    if best is None:
        raise StopStep("hgb_inner", f"{tag}: no finite inner ρ")
    log(f"[hgb] {tag} selected {best} inner_ρ={best_score:+.3f}")
    return best, pd.DataFrame(rows)


def _fit_predict(Xtr, ytr, Xte, params):
    est = _hgb(**params)
    est.fit(np.asarray(Xtr, np.float32), np.asarray(ytr, np.float64))
    return est.predict(np.asarray(Xte, np.float32))


def _rho_null_fixed(obs, y, pred, n_perm=N_PERM, seed=FIBRO2_SEED):
    """Donor-level permutation of ages; HGB scores held fixed (house-rule null)."""
    rng = np.random.default_rng(int(seed))
    pred = np.asarray(pred, float)
    nulls = []
    for _ in range(int(n_perm)):
        yp = permute_age_within_center(obs, rng)
        nulls.append(spearman_safe(yp, pred))
    nulls = np.asarray(nulls, float)
    rho_null = float(np.nanmedian(nulls)) if np.isfinite(nulls).any() else np.nan
    sc_rho = spearman_safe(np.asarray(y, float), pred)
    p = permutation_p(sc_rho, nulls, greater=True)
    return rho_null, p, nulls


def run_cv_hgb(pack, log):
    rows, skipped, inner_all = [], [], []
    oof = {}
    t0 = time.time()
    for site in TRANSFER_SITES + (EXCLUDE_SITE,):
        site_pack, skip = _site_folds(pack, site, log)
        if skip is not None:
            skipped.append(skip)
            log(f"[cv] skip site {site}: {skip}")
            continue
        X, obs, folds = site_pack["X"], site_pack["obs"], site_pack["folds"]
        y = obs.age_mid.to_numpy(float)
        batches = obs.SMNABTCH.astype(str).to_numpy()
        pred = np.full(len(obs), np.nan)
        fold_params = []
        for fi, (tr, te, te_lv) in enumerate(folds):
            Xtr, Xte, ytr, yte, meta = prepare_fold_X(X, obs, tr, te, "raw", "ridge")
            params, inner_df = _select_params(Xtr, ytr, batches[tr], log, f"{site}_fold{fi}")
            inner_df["site"] = site
            inner_df["fold"] = fi
            inner_all.append(inner_df)
            pred[te] = _fit_predict(Xtr, ytr, Xte, params)
            fold_params.append(dict(fold=fi, params=params, test_batches=",".join(te_lv)))
        sc = pred_scores(y, pred)
        rho_null, p, null_rho = _rho_null_fixed(obs, y, pred)
        rng_b = np.random.default_rng(FIBRO2_BOOT)
        ci = bootstrap_rho_ci(y, pred, rng_b, n_boot=N_BOOT)
        rec = dict(
            site=site, regime="raw", method="hgb",
            rho=sc["rho"], rho_null=rho_null,
            rho_p=p,
            rho_ci_lo=ci["p025"], rho_ci_hi=ci["p975"],
            r=sc["r"], r2=sc["r2"], cal_r2=sc["cal_r2"],
            n=int(len(obs)), n_folds=len(folds), n_perm=N_PERM, n_boot=N_BOOT,
            pass_rho=pass_rho_bar(sc["rho"], rho_null),
            selected_params=str(fold_params),
            null_procedure="donor-level age permutation; HGB OOF scores held fixed",
        )
        rows.append(rec)
        oof[site] = dict(y=y, pred=pred, obs=obs)
        log(f"[cv] {site} hgb ρ={sc['rho']:+.3f} null={rho_null:+.3f} p={rec['rho_p']:+.3f} "
            f"CI=[{ci['p025']:+.3f}, {ci['p975']:+.3f}] n={len(obs)} ({time.time()-t0:.0f}s)")
        dump_json(FIBRO2_DIR / "tc_cv_partial.json", jsonable(dict(rows=rows, skipped=skipped)))
        pd.DataFrame(rows).to_csv(FIBRO2_DIR / "tc_cv.csv", index=False)
        progress_snapshot(f"T-C within-site CV after {site}", extra=f"hgb {site} done")
    if inner_all:
        pd.concat(inner_all, ignore_index=True).to_csv(FIBRO2_DIR / "tc_inner_cv.csv", index=False)
    return rows, skipped, oof


def run_transfer_hgb(pack, log):
    obs, X = pack["obs"], pack["X"]
    cvec = obs.SMCENTER.astype(str).to_numpy()
    donors = obs.donor.astype(str).to_numpy()
    batches = obs.SMNABTCH.astype(str).to_numpy()
    rows = []
    t0 = time.time()
    for train_site, test_site in (("B1", "C1"), ("C1", "B1")):
        tr = np.flatnonzero(cvec == train_site)
        te = np.flatnonzero(cvec == test_site)
        if len(tr) < TRANSFER_MIN_N or len(te) < TRANSFER_MIN_N:
            raise StopStep(
                "transfer",
                f"{train_site}→{test_site}: n_train={len(tr)} n_test={len(te)}; "
                f"bar is n≥{TRANSFER_MIN_N}. Not lowering the bar.",
            )
        assert_disjoint(donors[tr], donors[te], what="donor",
                        where=f"T-C transfer {train_site}->{test_site}")
        Xtr, Xte, ytr, yte, meta = prepare_fold_X(X, obs, tr, te, "raw", "ridge")
        params, inner_df = _select_params(Xtr, ytr, batches[tr], log, f"{train_site}->{test_site}")
        inner_df["direction"] = f"{train_site}→{test_site}"
        inner_df.to_csv(FIBRO2_DIR / f"tc_inner_{train_site}_{test_site}.csv", index=False)
        pred = _fit_predict(Xtr, ytr, Xte, params)
        sc = pred_scores(yte, pred)
        obs_te = obs.iloc[te].reset_index(drop=True)
        rho_null, p, nulls = _rho_null_fixed(obs_te, yte, pred)
        rng_b = np.random.default_rng(FIBRO2_BOOT)
        ci = bootstrap_rho_ci(yte, pred, rng_b, n_boot=N_BOOT)
        rec = dict(
            direction=f"{train_site}→{test_site}", train_site=train_site, test_site=test_site,
            regime="raw", method="hgb",
            rho=sc["rho"], rho_null=rho_null,
            rho_p=p,
            rho_ci_lo=ci["p025"], rho_ci_hi=ci["p975"],
            r=sc["r"], r2=sc["r2"], cal_r2=sc["cal_r2"],
            pass_rho=pass_rho_bar(sc["rho"], rho_null),
            n_train=int(len(tr)), n_test=int(len(te)),
            n_perm=N_PERM, n_boot=N_BOOT, params=str(params),
            excluded_site=EXCLUDE_SITE, excluded_reason=EXCLUDE_SITE_REASON,
            null_procedure="donor-level age permutation on the test site; transferred HGB scores held fixed",
        )
        rows.append(rec)
        pd.DataFrame(rows).to_csv(FIBRO2_DIR / "tc_transfer.csv", index=False)
        log(f"[transfer] {train_site}→{test_site} hgb ρ={sc['rho']:+.3f} null={rho_null:+.3f} "
            f"p={rec['rho_p']:+.3f} CI=[{ci['p025']:+.3f}, {ci['p975']:+.3f}] "
            f"({time.time()-t0:.0f}s)")
        progress_snapshot(f"T-C transfer after {train_site}→{test_site}")
    return rows


def load_ridge_side_by_side():
    cv_p = FIBRO_DIR / "stage1_cv.csv"
    tr_p = FIBRO_DIR / "stage1_transfer.csv"
    if not cv_p.exists() or not tr_p.exists():
        raise StopStep("tc", "missing results/fibro/stage1_cv.csv or stage1_transfer.csv")
    cv = pd.read_csv(cv_p)
    tr = pd.read_csv(tr_p)
    ridge_cv = cv[(cv.regime == "raw") & (cv.method == "ridge")].copy()
    ridge_tr = tr[(tr.regime == "raw") & (tr.method == "ridge")].copy()
    return ridge_cv, ridge_tr


def tc_reading(hgb_tr, ridge_tr):
    """Compare boosting transfer ρ to ridge bootstrap CI, both directions."""
    def row(df, direction):
        sub = df[df.direction.astype(str) == direction]
        if sub.empty:
            return None
        return sub.iloc[0]

    comps = []
    both_in = True
    both_above = True
    both_below = True
    for d in ("B1→C1", "C1→B1"):
        h, r = row(hgb_tr, d), row(ridge_tr, d)
        if h is None or r is None:
            return dict(key="missing", text=f"missing transfer row {d}", linearity_limiting=None)
        rho_h = float(h.rho)
        lo, hi = float(r.rho_ci_lo), float(r.rho_ci_hi)
        inside = bool(lo <= rho_h <= hi)
        above = bool(rho_h > hi)
        below = bool(rho_h < float(r.rho))
        both_in = both_in and inside
        both_above = both_above and above
        both_below = both_below and below
        comps.append(dict(
            direction=d, hgb_rho=rho_h, ridge_rho=float(r.rho),
            ridge_ci_lo=lo, ridge_ci_hi=hi, inside=inside, above=above, below_point=below,
        ))
    if both_in:
        text = (
            "boosting transfer ρ within the ridge bootstrap CI → linearity is not the "
            "limiting factor; ridge stays."
        )
        return dict(key="linearity_not_limiting", text=text, linearity_limiting=False, comps=comps)
    if both_above:
        text = (
            "boosting transfer ρ materially above the ridge CI (upper bound exceeded in "
            "both directions) → the linear direction is leaving signal on the table. "
            "Do not rerun Stage 2 with boosting. A direction-preserving nonlinear extension "
            "(e.g. a kernel or autoencoder latent with a linear age axis inside it) would "
            "be the next design question."
        )
        return dict(key="linear_leaves_signal", text=text, linearity_limiting=True, comps=comps)
    # "below ridge" — point estimate below ridge point, or mixed
    h_b1 = float(row(hgb_tr, "B1→C1").rho)
    h_c1 = float(row(hgb_tr, "C1→B1").rho)
    r_b1 = float(row(ridge_tr, "B1→C1").rho)
    r_c1 = float(row(ridge_tr, "C1→B1").rho)
    if h_b1 < r_b1 and h_c1 < r_c1:
        text = "boosting transfer ρ below ridge → note it and move on."
        return dict(key="boosting_below_ridge", text=text, linearity_limiting=False, comps=comps)
    text = (
        "boosting vs ridge transfer ρ is mixed across directions (not both inside CI, "
        "not both above CI, not both below ridge). See the table. Not averaged. Ridge stays "
        "as the Stage 2 instrument because a nonlinear model is not a direction."
    )
    return dict(key="mixed", text=text, linearity_limiting=None, comps=comps)


def run_tc(log=None):
    if not PREREG_TC_FLAG.exists():
        raise StopStep("prereg", "PREREG_TC.flag missing")
    close_log = False
    if log is None:
        log = Logger(FIBRO2_DIR / "tc_report.txt")
        close_log = True
    fibro2_log_banner(log, "T-C")
    log(PREREG_TC)
    log("A nonlinear model yields a score, not a direction; the project's plane, angles, "
        "and 'age movement per unit identity loss' criterion all require a direction, so a "
        "nonlinear model cannot be substituted into Stage 2 without discarding the geometry.")
    log("HGB permutation null: donor-level shuffle of ages; fitted scores held fixed "
        "(n_perm=200, seed 20260914). Ridge Stage 1 numbers on disk used a refit-on-permuted-y "
        "null because ridge is cheap; T-C reading compares HGB transfer ρ to the ridge "
        "bootstrap CI, not to ridge p-values. Do not rescore GSE297234 with boosting.")
    pack = load_pack()
    log(f"[tc] pack n={len(pack['obs'])} genes={pack['X'].shape[1]} "
        f"D1 n={int((pack['obs'].SMCENTER.astype(str)==EXCLUDE_SITE).sum())}. {EXCLUDE_SITE_REASON}")
    ridge_cv, ridge_tr = load_ridge_side_by_side()
    cv_rows, skipped, oof = run_cv_hgb(pack, log)
    pd.DataFrame(cv_rows).to_csv(FIBRO2_DIR / "tc_cv.csv", index=False)
    dump_json(FIBRO2_DIR / "tc_cv_skipped.json", jsonable(skipped))
    tr_rows = run_transfer_hgb(pack, log)
    hgb_tr = pd.DataFrame(tr_rows)
    hgb_tr.to_csv(FIBRO2_DIR / "tc_transfer.csv", index=False)
    reading = tc_reading(hgb_tr, ridge_tr)
    # side-by-side table
    side = []
    for _, r in ridge_cv.iterrows():
        h = next((x for x in cv_rows if x["site"] == r.site), None)
        side.append(dict(
            cell="within_site", site=r.site, ridge_rho=r.rho, ridge_rho_null=r.rho_null,
            ridge_p=r.rho_p, ridge_ci_lo=r.rho_ci_lo, ridge_ci_hi=r.rho_ci_hi,
            ridge_cal_r2=r.cal_r2, ridge_r2=r.r2, n=r.n,
            hgb_rho=None if h is None else h["rho"],
            hgb_rho_null=None if h is None else h["rho_null"],
            hgb_p=None if h is None else h["rho_p"],
            hgb_ci_lo=None if h is None else h["rho_ci_lo"],
            hgb_ci_hi=None if h is None else h["rho_ci_hi"],
            hgb_cal_r2=None if h is None else h["cal_r2"],
            hgb_r2=None if h is None else h["r2"],
        ))
    for _, r in ridge_tr.iterrows():
        h = hgb_tr[hgb_tr.direction == r.direction]
        h = None if h.empty else h.iloc[0]
        side.append(dict(
            cell="transfer", site=r.direction, ridge_rho=r.rho, ridge_rho_null=r.rho_null,
            ridge_p=r.rho_p, ridge_ci_lo=r.rho_ci_lo, ridge_ci_hi=r.rho_ci_hi,
            ridge_cal_r2=r.cal_r2, ridge_r2=r.r2, n=r.n_test,
            hgb_rho=None if h is None else float(h.rho),
            hgb_rho_null=None if h is None else float(h.rho_null),
            hgb_p=None if h is None else float(h.rho_p),
            hgb_ci_lo=None if h is None else float(h.rho_ci_lo),
            hgb_ci_hi=None if h is None else float(h.rho_ci_hi),
            hgb_cal_r2=None if h is None else float(h.cal_r2),
            hgb_r2=None if h is None else float(h.r2),
        ))
    pd.DataFrame(side).to_csv(FIBRO2_DIR / "tc_side_by_side.csv", index=False)
    summary = dict(
        reading=reading, n_cv=len(cv_rows), n_transfer=len(tr_rows),
        skipped=skipped, seed=FIBRO2_SEED, boot_seed=FIBRO2_BOOT,
        direction_vs_score=(
            "A nonlinear model yields a score, not a direction; the project's plane, "
            "angles, and 'age movement per unit identity loss' criterion all require a "
            "direction, so a nonlinear model cannot be substituted into Stage 2 without "
            "discarding the geometry."
        ),
    )
    dump_json(FIBRO2_DIR / "tc_summary.json", jsonable(summary))
    man = load_manifest()
    man["status"] = "TC_DONE"
    man["tc_reading"] = reading
    save_manifest(man)
    progress_snapshot("extrapolation distances (report-only)", stop="TC_DONE")
    if close_log:
        log.close()
    return summary


if __name__ == "__main__":
    try:
        run_tc()
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise

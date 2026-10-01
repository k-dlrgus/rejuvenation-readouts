"""PART A — Leakage audit (blocking).

A1. Reproduce the within-site shuffle R² ≈ 0.195.
A2. Test donor leak, shuffle granularity, feature-selection leak, site structure, target centering.
A3. Fix whatever fired; re-run P3 / P4c under the corrected metric; side-by-side table.

STOP A if the corrected age-shuffle null is still above ~0.05.

Usage: python src/geometry_partA.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geometry_common import (  # noqa: E402
    GEO_DIR, GEO_FIG, GEO_SEED, N_PERM_A, N_MATCH_DRAWS, NULL_USABLE, DATASET_ID,
    Logger, dump_json, load_json, PHASE1_DIR,
    geo_log_banner, load_phase1_matrix, make_folds, assert_donor_folds,
    permute_age_per_pseudobulk, permute_age_within_donor_rows, donor_age_is_constant,
    intercept_only_oof, age_oof_svd, score_age_predictions,
    r2_age_on_site, centering_audit, load_p3_oof, permute_age_within_site_obs,
    age_oof_for_genes, select_genes_on_train, load_chrom_series,
)
from brain_phase1_models import identity_oof_for_genes  # noqa: E402


def _scores_to_json(sc):
    skip = {"per_type", "donor_pred", "pred", "y"}
    out = {}
    for k, v in sc.items():
        if k in skip:
            continue
        if isinstance(v, (np.floating, float)):
            out[k] = float(v) if np.isfinite(v) else None
        else:
            out[k] = v
    return out


def _run_oof(Y, obs, gene_idx, folds, method="svd", rng=None):
    if method == "original":
        if rng is None:
            raise ValueError("original ridge needs rng")
        out = age_oof_for_genes(Y, obs, gene_idx, folds, rng, method="ridge")
        pred, y = out["pred"], obs.age.to_numpy(float)
    else:
        out = age_oof_svd(Y, obs, gene_idx, folds)
        pred, y = out["pred"], out["y"]
    sc = score_age_predictions(obs, pred, y=y)
    sc["pred"] = pred
    sc["y"] = y
    sc["method"] = method
    return sc


def _shuffle_null(Y, obs, gene_idx, folds, rng, n_perm, method="svd", how="donor_within_site"):
    rows = []
    for i in range(n_perm):
        if how == "donor_within_site":
            obs_p = permute_age_within_site_obs(obs, rng)
        elif how == "per_pseudobulk":
            obs_p = permute_age_per_pseudobulk(obs, rng)
        elif how == "within_donor_rows":
            obs_p = permute_age_within_donor_rows(obs, rng)
        else:
            raise ValueError(how)
        y_p = obs_p.age.to_numpy(float)
        if method == "original":
            pred = age_oof_for_genes(Y, obs_p, gene_idx, folds, rng, method="ridge")["pred"]
        else:
            pred = age_oof_svd(Y, obs_p, gene_idx, folds, y_override=y_p)["pred"]
        sc = score_age_predictions(obs_p, pred, y=y_p)
        rows.append(dict(
            perm=i, how=how, method=method,
            pooled=sc["pooled_median_type_r2"],
            site_strat=sc["site_strat_median_type_r2"],
            joint_pooled=sc["joint_pooled_r2"],
            joint_site_strat=sc["joint_site_strat_r2"],
        ))
    return pd.DataFrame(rows)


def _summarize_null(tab):
    if tab is None or len(tab) == 0:
        return {}
    out = {}
    for col in ("pooled", "site_strat"):
        v = tab[col].to_numpy(float)
        out[f"{col}_mean"] = float(np.nanmean(v))
        out[f"{col}_median"] = float(np.nanmedian(v))
        out[f"{col}_p95"] = float(np.nanpercentile(v, 95))
        out[f"{col}_min"] = float(np.nanmin(v))
        out[f"{col}_max"] = float(np.nanmax(v))
        out[f"{col}_n"] = int(np.isfinite(v).sum())
    return out


def run():
    rng = np.random.default_rng(GEO_SEED)
    log = Logger(GEO_DIR / "partA_report.txt")
    try:
        return _run(log, rng)
    finally:
        log.close()


def _run(log, rng):
    geo_log_banner(log, "PART A — leakage audit")
    data = load_phase1_matrix(log)
    Y, genes, obs = data["Y"], data["genes"], data["obs"]
    A_idx, I_idx = data["A_idx"], data["I_idx"]
    extra = data["extra"]

    rng_folds = np.random.default_rng(GEO_SEED)
    schemes = make_folds(obs, rng_folds)
    log(f"[folds] constructed with seed {GEO_SEED}: loso={len(schemes['loso'])} "
        f"within_site={len(schemes['within_site'])}")

    # ------------------------------------------------------------------ A2.1 donor leakage
    log("\n" + "-" * 80)
    log("[A2] H1 donor leakage: is any donor's pseudobulks split across folds?")
    leak_rec = {}
    for name, folds in schemes.items():
        rec = assert_donor_folds(obs, folds, log, name)
        leak_rec[name] = rec
        log(f"   {name}: n_folds={rec['n_folds']} leaked_folds={rec['n_leaked_folds']}")
    donor_leak_fires = any(v["n_leaked_folds"] > 0 for v in leak_rec.values())
    log(f"   FIRES? {donor_leak_fires}")

    # ------------------------------------------------------------------ A2.2 age constancy / shuffle granularity
    log("\n" + "-" * 80)
    log("[A2] H2 shuffle granularity")
    age_const = donor_age_is_constant(obs)
    n_multi = int((obs.groupby("donor")["age"].nunique() > 1).sum())
    log(f"   donor age is constant across a donor's pseudobulks: {age_const} "
        f"(donors with >1 distinct age: {n_multi})")
    obs_sh = permute_age_within_site_obs(obs, np.random.default_rng(GEO_SEED))
    don_true = obs.groupby("donor").age.first()
    don_sh = obs_sh.groupby("donor").age.first()
    n_changed = int((don_true != don_sh).sum())
    within_donor_nunique_after = int((obs_sh.groupby("donor").age.nunique() > 1).sum())
    site_mean_true = obs.groupby(["Source", "donor"]).age.first().groupby("Source").mean()
    site_mean_sh = obs_sh.groupby(["Source", "donor"]).age.first().groupby("Source").mean()
    log(f"   permute_age_within_site_obs: donors with changed age={n_changed}/{len(don_true)}; "
        f"donors with non-constant shuffled age={within_donor_nunique_after}")
    log(f"   site mean ages TRUE: { {k: round(float(v), 3) for k, v in site_mean_true.items()} }")
    log(f"   site mean ages SHUF: { {k: round(float(v), 3) for k, v in site_mean_sh.items()} }")
    granularity_is_donor = (within_donor_nunique_after == 0) and (n_changed > 0)
    log(f"   current helper IS donor-level within-site permutation: {granularity_is_donor}")
    log(f"   FIRES (wrong granularity in current code)? {not granularity_is_donor}")

    r2_site = r2_age_on_site(obs)
    log(f"   R²(donor age ~ site) = {r2_site:+.3f}  (between-bank age gap)")

    # ------------------------------------------------------------------ A1 reproduce original shuffle scoring
    log("\n" + "-" * 80)
    log("[A1] Reproduce within-site shuffle R² under ORIGINAL pooled scoring")
    reported = load_json(PHASE1_DIR / "p4_summary.json")
    log(f"   stored P4 within_site shuffle_age_r2 = {reported['within_site']['shuffle_age_r2']:+.6f}")
    log(f"   stored P4 loso shuffle_age_r2        = {reported['loso']['shuffle_age_r2']:+.6f}")

    a1_rows = []
    for name, folds in schemes.items():
        pred_int = intercept_only_oof(obs, folds)
        sc_int = score_age_predictions(obs, pred_int)
        log(f"   [{name}] intercept-only  pooled={sc_int['pooled_median_type_r2']:+.3f}  "
            f"site_strat={sc_int['site_strat_median_type_r2']:+.3f}")
        a1_rows.append(dict(scheme=name, model="intercept_only_true_age",
                            pooled=sc_int["pooled_median_type_r2"],
                            site_strat=sc_int["site_strat_median_type_r2"],
                            joint_pooled=sc_int["joint_pooled_r2"],
                            joint_site_strat=sc_int["joint_site_strat_r2"]))
        obs_sh1 = permute_age_within_site_obs(obs, np.random.default_rng(GEO_SEED + 1))
        pred_int_sh = intercept_only_oof(obs_sh1, folds)
        sc_int_sh = score_age_predictions(obs_sh1, pred_int_sh, y=obs_sh1.age.to_numpy(float))
        log(f"   [{name}] intercept-only SHUFFLE pooled={sc_int_sh['pooled_median_type_r2']:+.3f}  "
            f"site_strat={sc_int_sh['site_strat_median_type_r2']:+.3f}")
        a1_rows.append(dict(scheme=name, model="intercept_only_shuffled_age",
                            pooled=sc_int_sh["pooled_median_type_r2"],
                            site_strat=sc_int_sh["site_strat_median_type_r2"],
                            joint_pooled=sc_int_sh["joint_pooled_r2"],
                            joint_site_strat=sc_int_sh["joint_site_strat_r2"]))

    log("   original-ridge A-gene shuffle (one draw, donor-level within site)...")
    rng_a1 = np.random.default_rng(GEO_SEED + 17)
    for name, folds in schemes.items():
        obs_sh1 = permute_age_within_site_obs(obs, rng_a1)
        out = age_oof_for_genes(Y, obs_sh1, A_idx, folds, rng_a1, method="ridge")
        sc = score_age_predictions(obs_sh1, out["pred"], y=obs_sh1.age.to_numpy(float))
        log(f"   [{name}] original-ridge shuffle  pooled={sc['pooled_median_type_r2']:+.3f}  "
            f"site_strat={sc['site_strat_median_type_r2']:+.3f}  "
            f"(P4 reported pooled={reported[name]['shuffle_age_r2']:+.3f})")
        a1_rows.append(dict(scheme=name, model="original_ridge_A_shuffle_one_draw",
                            pooled=sc["pooled_median_type_r2"],
                            site_strat=sc["site_strat_median_type_r2"],
                            joint_pooled=sc["joint_pooled_r2"],
                            joint_site_strat=sc["joint_site_strat_r2"]))
        out["per_type"].assign(scheme=name).to_csv(
            GEO_DIR / f"partA_a1_original_shuffle_per_type_{name}.csv", index=False)

    pd.DataFrame(a1_rows).to_csv(GEO_DIR / "partA_a1_reproduce.csv", index=False)

    # ------------------------------------------------------------------ A2.2 empirical granularity comparison
    log("\n   empirical shuffle-granularity comparison (SVD-LOO ridge, A-genes)...")
    ws = schemes["within_site"]
    gran_tabs = {}
    for how in ("donor_within_site", "per_pseudobulk", "within_donor_rows"):
        n = 5 if how != "donor_within_site" else N_PERM_A
        tab = _shuffle_null(Y, obs, A_idx, ws, rng, n_perm=n, method="svd", how=how)
        gran_tabs[how] = tab
        sm = _summarize_null(tab)
        log(f"   how={how:22s}  pooled mean={sm['pooled_mean']:+.3f}  "
            f"site_strat mean={sm['site_strat_mean']:+.3f}  n={sm['pooled_n']}")
        tab.to_csv(GEO_DIR / f"partA_shuffle_{how}_within_site.csv", index=False)
    loso_null = _shuffle_null(Y, obs, A_idx, schemes["loso"], rng, n_perm=N_PERM_A,
                              method="svd", how="donor_within_site")
    loso_null.to_csv(GEO_DIR / "partA_shuffle_donor_within_site_loso.csv", index=False)
    log(f"   LOSO donor-level null  pooled mean={_summarize_null(loso_null)['pooled_mean']:+.3f}  "
        f"site_strat mean={_summarize_null(loso_null)['site_strat_mean']:+.3f}")

    # ------------------------------------------------------------------ A2.4 site structure
    log("\n" + "-" * 80)
    log("[A2] H4 site structure: within-site k-fold + pooled R² recovers the bank age gap")
    ws_int = [r for r in a1_rows if r["scheme"] == "within_site"
              and r["model"] == "intercept_only_true_age"][0]
    site_structure_fires = abs(ws_int["pooled"] - r2_site) < 0.15 and ws_int["pooled"] > 0.10
    log(f"   intercept-only within-site pooled R²={ws_int['pooled']:+.3f} vs R²(age~site)={r2_site:+.3f}")
    log(f"   intercept-only within-site site_strat R²={ws_int['site_strat']:+.3f} (should be ~0)")
    log(f"   FIRES? {site_structure_fires}")

    # ------------------------------------------------------------------ A2.5 target centering
    log("\n" + "-" * 80)
    log("[A2] H5 target centering: ridge uses training-fold intercept; R² uses eval-set mean")
    p3_ws = load_p3_oof("within_site")
    assert len(p3_ws) == len(obs), (len(p3_ws), len(obs))
    cent = centering_audit(obs, p3_ws.pred_A.to_numpy(float), y=p3_ws.age.to_numpy(float))
    log(f"   P3 A-gene OOF pooled (eval-set mean in ss_tot) R²={cent['r2_eval_mean']:+.3f}")
    log(f"   same SS_res but ss_tot vs full-cohort mean: R²={cent['r2_full_cohort_mean_ss_tot']:+.3f}")
    log("   ridge_predict: StandardScaler+Ridge(fit_intercept=True) fit on TRAIN only — "
        "age is NOT centered with the full-cohort mean inside the model.")
    log("   r2_mae uses the evaluation vector's own mean. The leak is that the evaluation "
        "vector for within-site CV concatenates TWO sites, each with its own fitted intercept.")
    centering_fires = False
    log(f"   FIRES as 'full-cohort mean centering'? {centering_fires}  "
        "(the real metric bug is H4, not a global-mean subtraction in the estimator)")

    # ------------------------------------------------------------------ A2.3 feature-selection leakage
    log("\n" + "-" * 80)
    log("[A2] H3 feature-selection leakage: P4c used FROZEN full-data A/I sets")
    log("   P3 nested selection inside each training fold (declared quantile LEVELS + P2 filters).")
    log("   Testing nested vs frozen under BOTH scoring rules.")
    chrom_s = load_chrom_series(genes.gene_id.to_numpy())
    nested = {}
    for name, folds in schemes.items():
        log(f"   nested gene selection on {name} ({len(folds)} folds)...")
        infos = []
        for i, fold in enumerate(folds):
            info = select_genes_on_train(Y, obs, genes, extra, fold[0], chrom_s=chrom_s)
            infos.append(info)
            tag = fold[2] if len(fold) > 2 else i
            log(f"      fold {i} {tag}: nested A={info['n_A']} I={info['n_I']}")
        nested[name] = infos
        dump_json(GEO_DIR / f"partA_nested_gene_counts_{name}.json",
                  [dict(n_A=i["n_A"], n_I=i["n_I"]) for i in infos])

    # ------------------------------------------------------------------ A3 refit
    log("\n" + "-" * 80)
    log("[A3] Re-run age models (frozen A, frozen I, nested A, nested I, size-matched I)")
    log("     SVD-LOO ridge; original sklearn ridge used above for A1 reproduction only.")

    rec = {}
    for name, folds in schemes.items():
        log(f"\n   === scheme {name} ===")
        rec[name] = {}
        scA = _run_oof(Y, obs, A_idx, folds, method="svd")
        rec[name]["frozen_A"] = _scores_to_json(scA)
        scA["per_type"].assign(scheme=name, model="frozen_A").to_csv(
            GEO_DIR / f"partA_per_type_{name}_frozen_A.csv", index=False)
        pd.DataFrame({"donor": obs.donor, "celltype": obs.celltype, "Source": obs.Source,
                      "age": obs.age, "pred": scA["pred"]}).to_csv(
            GEO_DIR / f"partA_oof_{name}_frozen_A.csv", index=False)
        log(f"   frozen A  pooled={scA['pooled_median_type_r2']:+.3f}  "
            f"site_strat={scA['site_strat_median_type_r2']:+.3f}")

        scI = _run_oof(Y, obs, I_idx, folds, method="svd")
        rec[name]["frozen_I"] = _scores_to_json(scI)
        pd.DataFrame({"donor": obs.donor, "celltype": obs.celltype, "Source": obs.Source,
                      "age": obs.age, "pred": scI["pred"]}).to_csv(
            GEO_DIR / f"partA_oof_{name}_frozen_I.csv", index=False)
        log(f"   frozen I  pooled={scI['pooled_median_type_r2']:+.3f}  "
            f"site_strat={scI['site_strat_median_type_r2']:+.3f}")

        A_list = [i["A_idx"] for i in nested[name]]
        I_list = [i["I_idx"] for i in nested[name]]
        scAn = _run_oof(Y, obs, A_list, folds, method="svd")
        scIn = _run_oof(Y, obs, I_list, folds, method="svd")
        rec[name]["nested_A"] = _scores_to_json(scAn)
        rec[name]["nested_I"] = _scores_to_json(scIn)
        rec[name]["nested_n_A_median"] = float(np.median([i["n_A"] for i in nested[name]]))
        rec[name]["nested_n_I_median"] = float(np.median([i["n_I"] for i in nested[name]]))
        pd.DataFrame({"donor": obs.donor, "celltype": obs.celltype, "Source": obs.Source,
                      "age": obs.age, "pred_nested_A": scAn["pred"],
                      "pred_nested_I": scIn["pred"]}).to_csv(
            GEO_DIR / f"partA_oof_{name}_nested.csv", index=False)
        log(f"   nested A  pooled={scAn['pooled_median_type_r2']:+.3f}  "
            f"site_strat={scAn['site_strat_median_type_r2']:+.3f}  "
            f"median n_A={rec[name]['nested_n_A_median']:.0f}")
        log(f"   nested I  pooled={scIn['pooled_median_type_r2']:+.3f}  "
            f"site_strat={scIn['site_strat_median_type_r2']:+.3f}  "
            f"median n_I={rec[name]['nested_n_I_median']:.0f}")

        log(f"   size-matched I ({N_MATCH_DRAWS} draws, k={len(A_idx)})...")
        match_rows = []
        for d in range(N_MATCH_DRAWS):
            sub = rng.choice(I_idx, size=min(len(A_idx), len(I_idx)), replace=False)
            sc = _run_oof(Y, obs, sub, folds, method="svd")
            match_rows.append(dict(draw=d,
                                   pooled=sc["pooled_median_type_r2"],
                                   site_strat=sc["site_strat_median_type_r2"],
                                   joint_pooled=sc["joint_pooled_r2"],
                                   joint_site_strat=sc["joint_site_strat_r2"]))
        match_tab = pd.DataFrame(match_rows)
        match_tab.to_csv(GEO_DIR / f"partA_matched_I_{name}.csv", index=False)
        rec[name]["matched_I"] = dict(
            pooled=float(match_tab.pooled.median()),
            site_strat=float(match_tab.site_strat.median()),
            joint_pooled=float(match_tab.joint_pooled.median()),
            joint_site_strat=float(match_tab.joint_site_strat.median()),
            n_draws=N_MATCH_DRAWS, k=int(min(len(A_idx), len(I_idx))),
        )
        log(f"   matched I median pooled={rec[name]['matched_I']['pooled']:+.3f}  "
            f"site_strat={rec[name]['matched_I']['site_strat']:+.3f}")

        log("   identity frozen I / A...")
        id_I = identity_oof_for_genes(Y, obs, I_idx, folds, rng)
        id_A = identity_oof_for_genes(Y, obs, A_idx, folds, rng)
        rec[name]["id_I_acc"] = float(id_I["acc"])
        rec[name]["id_A_acc"] = float(id_A["acc"])
        rec[name]["id_I_f1"] = float(id_I["macro_f1"])
        rec[name]["id_chance"] = float(id_I["chance"])
        log(f"   Acc_I={id_I['acc']:.3f}  Acc_A={id_A['acc']:.3f}  chance={id_I['chance']:.3f}")
        rec[name]["reported_p4"] = dict(
            r2_A=reported[name]["r2_A"],
            r2_I_matched=reported[name]["r2_I_matched"],
            r2_I_full=reported[name]["r2_I_full"],
            shuffle_age_r2=reported[name]["shuffle_age_r2"],
            acc_I=reported[name]["acc_I"],
            acc_A=reported[name]["acc_A"],
            verdict=reported[name]["verdict"],
        )

        p3 = load_p3_oof(name)
        sc_p3A = score_age_predictions(obs, p3.pred_A.to_numpy(float), y=p3.age.to_numpy(float))
        sc_p3U = score_age_predictions(obs, p3.pred_unrest.to_numpy(float), y=p3.age.to_numpy(float))
        rec[name]["p3_saved_A"] = _scores_to_json(sc_p3A)
        rec[name]["p3_saved_unrest"] = _scores_to_json(sc_p3U)
        log(f"   P3 saved nested-A OOF  pooled={sc_p3A['pooled_median_type_r2']:+.3f}  "
            f"site_strat={sc_p3A['site_strat_median_type_r2']:+.3f}")
        log(f"   P3 saved unrest OOF    pooled={sc_p3U['pooled_median_type_r2']:+.3f}  "
            f"site_strat={sc_p3U['site_strat_median_type_r2']:+.3f}")

    rec["within_site"]["shuffle_A_null"] = _summarize_null(gran_tabs["donor_within_site"])
    rec["loso"]["shuffle_A_null"] = _summarize_null(loso_null)
    rec["within_site"]["shuffle_A_null_per_pseudobulk"] = _summarize_null(gran_tabs["per_pseudobulk"])
    rec["within_site"]["shuffle_A_null_within_donor_rows"] = _summarize_null(
        gran_tabs["within_donor_rows"])

    log("\n   I-gene donor-level shuffle null (10 draws / scheme)...")
    for name, folds in schemes.items():
        tabI = _shuffle_null(Y, obs, I_idx, folds, rng, n_perm=10, method="svd",
                             how="donor_within_site")
        tabI.to_csv(GEO_DIR / f"partA_shuffle_I_{name}.csv", index=False)
        rec[name]["shuffle_I_null"] = _summarize_null(tabI)
        log(f"   [{name}] I-shuffle pooled mean={rec[name]['shuffle_I_null']['pooled_mean']:+.3f}  "
            f"site_strat mean={rec[name]['shuffle_I_null']['site_strat_mean']:+.3f}")

    log("   nested-A donor-level shuffle (5 draws / scheme)...")
    for name, folds in schemes.items():
        A_list = [i["A_idx"] for i in nested[name]]
        tabN = _shuffle_null(Y, obs, A_list, folds, rng, n_perm=5, method="svd",
                             how="donor_within_site")
        tabN.to_csv(GEO_DIR / f"partA_shuffle_nestedA_{name}.csv", index=False)
        rec[name]["shuffle_nestedA_null"] = _summarize_null(tabN)
        log(f"   [{name}] nestedA-shuffle pooled mean="
            f"{rec[name]['shuffle_nestedA_null']['pooled_mean']:+.3f}  "
            f"site_strat mean={rec[name]['shuffle_nestedA_null']['site_strat_mean']:+.3f}")

    fs_frozen = rec["within_site"]["frozen_A"]["pooled_median_type_r2"]
    fs_nested = rec["within_site"]["nested_A"]["pooled_median_type_r2"]
    fs_p3 = rec["within_site"]["p3_saved_A"]["pooled_median_type_r2"]
    fs_leak_fires = (fs_frozen - fs_nested) > 0.05 or (fs_frozen - fs_p3) > 0.05
    log(f"\n   frozen A pooled={fs_frozen:+.3f} vs nested SVD={fs_nested:+.3f} vs P3 saved={fs_p3:+.3f}")
    log(f"   FIRES (selection saw test folds, inflating frozen vs nested by >0.05)? {fs_leak_fires}")
    log("   (P4c used frozen sets by design so A vs I are size-and-selection comparable; "
        "P3 is the leakage-safer A-gene clock.)")

    hypotheses = [
        dict(id="H1_donor_leakage", fires=bool(donor_leak_fires),
             result="PASS: no donor in both sides of any fold" if not donor_leak_fires else "FAIL",
             note="asserted on loso and within_site folds"),
        dict(id="H2_shuffle_granularity", fires=bool(not granularity_is_donor),
             result=("PASS: permute_age_within_site_obs remaps a permuted donor-level age "
                     "onto every pseudobulk of that donor") if granularity_is_donor else "FAIL",
             note=("within-donor row shuffle is a no-op because age is donor-constant; "
                   "per-pseudobulk shuffle is NOT what P4 ran")),
        dict(id="H3_feature_selection_leakage", fires=bool(fs_leak_fires),
             result=(f"frozen pooled A R²={fs_frozen:+.3f}; nested={fs_nested:+.3f}; "
                     f"P3 saved={fs_p3:+.3f}"),
             note="P4c used frozen full-data sets; P3 nested inside folds"),
        dict(id="H4_site_structure", fires=bool(site_structure_fires),
             result=(f"intercept-only within-site pooled R²={ws_int['pooled']:+.3f} "
                     f"matches R²(age~site)={r2_site:+.3f}; site-stratified={ws_int['site_strat']:+.3f}"),
             note=("within-site k-fold never mixes banks in a fold, so the intercept equals the "
                   "training bank's mean age; pooled-across-banks R² recovers the bank gap, "
                   "which the within-site age permutation PRESERVES")),
        dict(id="H5_target_centering", fires=bool(centering_fires),
             result=("PASS as a full-cohort-mean bug: Ridge(fit_intercept=True) on train; "
                     "r2_mae uses the evaluation vector's mean"),
             note="the eval vector for within-site CV is two banks concatenated — see H4"),
    ]
    hyp_tab = pd.DataFrame(hypotheses)
    hyp_tab.to_csv(GEO_DIR / "partA_hypotheses.csv", index=False)
    log("\n[A2] hypothesis summary:")
    for h in hypotheses:
        log(f"   {h['id']}: FIRES={h['fires']}  {h['result']}")

    fired = [h["id"] for h in hypotheses if h["fires"]]
    log(f"\n   FIRED: {fired if fired else 'none'}")
    fix = (
        "Primary fix: report site-stratified median-over-types R² for the within-site scheme "
        "(mean of the two banks' R² per cell type, then median over types). "
        "LOSO is left on the original pooled/test-site metric — intercept mismatch across "
        "banks is a real transfer failure, not a scoring artifact. "
        "P4c frozen vs nested is reported side-by-side; frozen remains the size-matched cross-test."
    )
    log("   FIX: " + fix)

    # ------------------------------------------------------------------ side-by-side table
    log("\n" + "-" * 80)
    log("[A3] side-by-side: reported-in-P4c vs corrected")
    side = []
    p3sum = load_json(PHASE1_DIR / "p3_summary.json")
    for name in ("within_site", "loso"):
        r = rec[name]
        null = r["shuffle_A_null"]
        rows = [
            ("P3 nested A-gene age R² (median types)", p3sum[name]["age_A_median_r2"],
             r["p3_saved_A"]["site_strat_median_type_r2"] if name == "within_site"
             else r["p3_saved_A"]["pooled_median_type_r2"],
             r["shuffle_nestedA_null"]["site_strat_mean"] if name == "within_site"
             else r["shuffle_nestedA_null"]["pooled_mean"]),
            ("P3 unrestricted age R² (median types)", p3sum[name]["age_unrest_median_r2"],
             r["p3_saved_unrest"]["site_strat_median_type_r2"] if name == "within_site"
             else r["p3_saved_unrest"]["pooled_median_type_r2"], np.nan),
            ("P4c R²_A frozen A-genes", r["reported_p4"]["r2_A"],
             r["frozen_A"]["site_strat_median_type_r2"] if name == "within_site"
             else r["frozen_A"]["pooled_median_type_r2"],
             null["site_strat_mean"] if name == "within_site" else null["pooled_mean"]),
            ("P4c R²_I size-matched I", r["reported_p4"]["r2_I_matched"],
             r["matched_I"]["site_strat"] if name == "within_site" else r["matched_I"]["pooled"],
             r["shuffle_I_null"]["site_strat_mean"] if name == "within_site"
             else r["shuffle_I_null"]["pooled_mean"]),
            ("P4c R²_I full I set", r["reported_p4"]["r2_I_full"],
             r["frozen_I"]["site_strat_median_type_r2"] if name == "within_site"
             else r["frozen_I"]["pooled_median_type_r2"],
             r["shuffle_I_null"]["site_strat_mean"] if name == "within_site"
             else r["shuffle_I_null"]["pooled_mean"]),
            ("P4c shuffle-age A-genes", r["reported_p4"]["shuffle_age_r2"],
             null["site_strat_mean"] if name == "within_site" else null["pooled_mean"],
             null["site_strat_mean"] if name == "within_site" else null["pooled_mean"]),
            ("P4c Acc_I", r["reported_p4"]["acc_I"], r["id_I_acc"], np.nan),
            ("P4c Acc_A (identity from A)", r["reported_p4"]["acc_A"], r["id_A_acc"], np.nan),
            ("nested A-gene age R² (this audit)", p3sum[name]["age_A_median_r2"],
             r["nested_A"]["site_strat_median_type_r2"] if name == "within_site"
             else r["nested_A"]["pooled_median_type_r2"],
             r["shuffle_nestedA_null"]["site_strat_mean"] if name == "within_site"
             else r["shuffle_nestedA_null"]["pooled_mean"]),
            ("nested I-gene age R² (this audit)", np.nan,
             r["nested_I"]["site_strat_median_type_r2"] if name == "within_site"
             else r["nested_I"]["pooled_median_type_r2"],
             r["shuffle_I_null"]["site_strat_mean"] if name == "within_site"
             else r["shuffle_I_null"]["pooled_mean"]),
        ]
        for metric, reported_v, corrected_v, null_v in rows:
            rv = float(reported_v) if np.isfinite(reported_v) else np.nan
            cv = float(corrected_v) if np.isfinite(corrected_v) else np.nan
            nv = float(null_v) if np.isfinite(null_v) else np.nan
            side.append(dict(
                scheme=name, metric=metric,
                reported_in_P4c=None if not np.isfinite(rv) else rv,
                corrected=None if not np.isfinite(cv) else cv,
                permutation_null=None if not np.isfinite(nv) else nv,
                null_usable=(bool(np.isfinite(nv) and abs(nv) <= NULL_USABLE)
                             if np.isfinite(nv) else None),
            ))
    side_tab = pd.DataFrame(side)
    side_tab.to_csv(GEO_DIR / "partA_side_by_side.csv", index=False)
    log(side_tab.to_string(index=False, float_format=lambda x: f"{x:+.3f}"))

    # ------------------------------------------------------------------ STOP A / verdict change
    ws_null = rec["within_site"]["shuffle_A_null"]["site_strat_mean"]
    loso_null_m = rec["loso"]["shuffle_A_null"]["pooled_mean"]
    stop_A = bool(np.isfinite(ws_null) and ws_null > NULL_USABLE)
    log("\n" + "-" * 80)
    log(f"[STOP A] corrected within-site shuffle null (site-strat mean) = {ws_null:+.3f}  "
        f"threshold {NULL_USABLE}")
    log(f"         LOSO shuffle null (pooled mean) = {loso_null_m:+.3f}")
    log(f"         STOP A (do not proceed to B/C)? {stop_A}")

    def _p4c_compare(name):
        r = rec[name]
        if name == "within_site":
            r2_A = r["frozen_A"]["site_strat_median_type_r2"]
            r2_I = r["matched_I"]["site_strat"]
            sh = r["shuffle_A_null"]["site_strat_mean"]
        else:
            r2_A = r["frozen_A"]["pooled_median_type_r2"]
            r2_I = r["matched_I"]["pooled"]
            sh = r["shuffle_A_null"]["pooled_mean"]
        acc_I, acc_A = r["id_I_acc"], r["id_A_acc"]
        delta = r2_A - r2_I
        ratio = r2_I / r2_A if (np.isfinite(r2_A) and r2_A != 0) else np.nan
        notes = []
        sound = True
        if np.isfinite(sh) and sh > 0.10:
            sound = False
            notes.append(f"shuffle-age R²={sh:.3f} > 0.10 (pipeline soundness fail)")
        elif np.isfinite(sh) and sh > NULL_USABLE:
            notes.append(f"shuffle-age R²={sh:.3f} above {NULL_USABLE} usability line")
        falsified = False
        if np.isfinite(ratio) and ratio >= 0.8:
            falsified = True
            notes.append(f"R²_I >= 0.8 × R²_A ({ratio:.3f})")
        if np.isfinite(delta) and delta < 0.10:
            falsified = True
            notes.append(f"R²_A − R²_I = {delta:.3f} < 0.10")
        if acc_A / acc_I >= 0.9:
            falsified = True
            notes.append(f"Acc_A >= 0.9 × Acc_I ({acc_A / acc_I:.3f})")
        if not sound:
            label = "soundness-fail"
        elif falsified:
            label = "falsified (not separable)"
        else:
            label = "inconclusive_or_supported"
        return dict(r2_A=r2_A, r2_I=r2_I, shuffle=sh, delta=delta, ratio=ratio,
                    acc_I=acc_I, acc_A=acc_A, verdict=label, detail="; ".join(notes),
                    sound=sound)

    v_ws = _p4c_compare("within_site")
    v_lo = _p4c_compare("loso")
    old_ws, old_lo = reported["within_site"]["verdict"], reported["loso"]["verdict"]
    log(f"\n[P4c verdict] within-site reported={old_ws}  corrected={v_ws['verdict']}")
    log(f"              {v_ws['detail']}")
    log(f"[P4c verdict] LOSO reported={old_lo}  corrected={v_lo['verdict']}")
    log(f"              {v_lo['detail']}")

    def _order(r2_A, r2_I):
        if not np.isfinite(r2_A) or not np.isfinite(r2_I):
            return "undefined"
        if r2_I > r2_A + 0.02:
            return "I_beats_A"
        if r2_A > r2_I + 0.02:
            return "A_beats_I"
        return "tie"

    old_order_ws = _order(reported["within_site"]["r2_A"], reported["within_site"]["r2_I_matched"])
    new_order_ws = _order(v_ws["r2_A"], v_ws["r2_I"])
    old_order_lo = _order(reported["loso"]["r2_A"], reported["loso"]["r2_I_matched"])
    new_order_lo = _order(v_lo["r2_A"], v_lo["r2_I"])
    if (new_order_ws != old_order_ws) or (new_order_lo != old_order_lo):
        change = "REVERSE or shift in A vs I ordering"
    elif old_ws != v_ws["verdict"] or old_lo != v_lo["verdict"]:
        if "soundness-fail" in (old_ws, old_lo) and "soundness-fail" not in (
                v_ws["verdict"], v_lo["verdict"]):
            change = ("STRENGTHEN: soundness-fail cleared; separability still falsified "
                      "on the cross-test")
        else:
            change = "verdict label changed; see details"
    else:
        change = "UNCHANGED label and A vs I ordering"
    log(f"[P4c change] {change}")
    log(f"   within-site order reported={old_order_ws} corrected={new_order_ws} "
        f"(R²_A {v_ws['r2_A']:+.3f} vs R²_I {v_ws['r2_I']:+.3f})")
    log(f"   LOSO order reported={old_order_lo} corrected={new_order_lo} "
        f"(R²_A {v_lo['r2_A']:+.3f} vs R²_I {v_lo['r2_I']:+.3f})")

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.0), sharey=False)
    for ax, name in zip(axes, ("within_site", "loso")):
        r = rec[name]
        labs = ["A frozen", "I matched", "I full", "shuffle null"]
        reported_v = [r["reported_p4"]["r2_A"], r["reported_p4"]["r2_I_matched"],
                      r["reported_p4"]["r2_I_full"], r["reported_p4"]["shuffle_age_r2"]]
        if name == "within_site":
            corr_v = [r["frozen_A"]["site_strat_median_type_r2"], r["matched_I"]["site_strat"],
                      r["frozen_I"]["site_strat_median_type_r2"],
                      r["shuffle_A_null"]["site_strat_mean"]]
            ylab = "median-over-types R²\n(corrected: site-stratified)"
        else:
            corr_v = [r["frozen_A"]["pooled_median_type_r2"], r["matched_I"]["pooled"],
                      r["frozen_I"]["pooled_median_type_r2"],
                      r["shuffle_A_null"]["pooled_mean"]]
            ylab = "median-over-types R²\n(LOSO pooled / original metric)"
        x = np.arange(len(labs))
        ax.bar(x - 0.2, reported_v, 0.4, color="0.7", label="reported P4c (pooled)")
        ax.bar(x + 0.2, corr_v, 0.4, color=["tab:red", "tab:blue", "tab:cyan", "0.35"],
               label="corrected")
        ax.axhline(0, c="k", lw=0.5)
        ax.axhline(NULL_USABLE, c="0.5", ls=":", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(labs, rotation=20, ha="right", fontsize=8)
        ax.set_title(name)
        ax.set_ylabel(ylab, fontsize=8)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(GEO_FIG / "partA_reported_vs_corrected.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    ws_tab = gran_tabs["donor_within_site"]
    ax.hist(ws_tab.pooled, bins=np.linspace(-0.2, 0.5, 15), color="0.6", alpha=0.7,
            label="pooled (original metric)")
    ax.hist(ws_tab.site_strat, bins=np.linspace(-0.2, 0.5, 15), color="tab:red", alpha=0.7,
            label="site-stratified (corrected)")
    ax.axvline(reported["within_site"]["shuffle_age_r2"], c="k", ls="--",
               label=f"P4 reported {reported['within_site']['shuffle_age_r2']:+.3f}")
    ax.axvline(0, c="k", lw=0.5)
    ax.set_xlabel("within-site A-gene shuffle R²")
    ax.set_ylabel("permutation draws")
    ax.set_title("A1/A3 shuffle null, donor-level within site")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(GEO_FIG / "partA_shuffle_null.png", dpi=140)
    plt.close(fig)

    summary = dict(
        seed=GEO_SEED, dataset_id=DATASET_ID,
        n_donors=data["n_donors"], n_types=data["n_types"],
        n_A=int(len(A_idx)), n_I=int(len(I_idx)),
        r2_age_on_site=float(r2_site),
        donor_age_constant=bool(age_const),
        hypotheses=hypotheses, fired=fired, fix=fix,
        stop_A=bool(stop_A),
        stop_A_reason=(f"corrected within-site shuffle null {ws_null:+.3f} > {NULL_USABLE}"
                       if stop_A else "corrected within-site shuffle null is usable"),
        null_usable_threshold=NULL_USABLE,
        n_perm_A=N_PERM_A,
        schemes=rec,
        verdict_within_site=v_ws, verdict_loso=v_lo,
        p4c_change=change,
        old_order=dict(within_site=old_order_ws, loso=old_order_lo),
        new_order=dict(within_site=new_order_ws, loso=new_order_lo),
        proceed_to_B=not bool(stop_A),
    )
    dump_json(GEO_DIR / "partA_summary.json", summary)
    (GEO_DIR / "partA_STOP.txt").write_text(
        ("STOP A\n" + summary["stop_A_reason"] + "\n") if stop_A
        else ("PASS A\n" + summary["stop_A_reason"] + "\n" + change + "\n"),
        encoding="utf-8",
    )
    log("\n[A] wrote results/geometry/partA_*")
    log("[A] done.")
    return summary


if __name__ == "__main__":
    run()

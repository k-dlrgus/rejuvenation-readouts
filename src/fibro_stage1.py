"""Stage 1 — GTEx V10 cultured-fibroblast age ruler.

Cohort, TMM/log-CPM, within-site SMNABTCH CV, B1↔C1 transfer, Stage 1 gate.
Does not open GSE297234. Does not modify existing FINDINGS*.md or FALSIFICATION.md.

Usage: python src/fibro_stage1.py [--cell cohort|preprocess|cv|transfer|gate|freeze|all]
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fibro_common import (  # noqa: E402
    FIBRO_DIR, FIBRO_FIG, FIBRO_PROC, FIBRO_SEED, BOOT_SEED, N_PERM, N_BOOT,
    N_FOLDS, TRANSFER_NULL_BAR, SMTSD_FIBRO, TRANSFER_SITES, EXCLUDE_SITE,
    EXCLUDE_SITE_REASON, SITE_CV_MIN_N, TRANSFER_MIN_N, MANIFEST_PATH,
    PREREG_STAGE1_FLAG, FROZEN_RULER, RELEASE, RNASEQ_PIPELINE, DBGAP,
    SMAFRZE_RNASEQ, SAMPLE_REQUIRED, PHENO_REQUIRED, EXPECTED_SEX,
    EXPECTED_DTHHRDY, EXPECTED_AGE_BINS, AGE_MID, MIN_COUNT, MIN_FRAC,
    EDGE_R_PRIOR, GTEX_RAW, FILES, DLPFC_GENES_CSV, CAVEAT_AGE_BIN,
    PREREG_STAGE1, StopStep, Logger, dump_json, load_json, strip_ensembl,
    sampid_to_subj, require_columns, require_levels, age_midpoint,
    tmm_norm_factors, log2_cpm_edger, gene_filter_mask, load_dlpfc_gene_universe,
    assert_disjoint, grouped_kfold, pred_scores, spearman_safe, pearson_safe,
    pass_rho_bar, fmt, fmt_ci, jsonable, load_manifest, save_manifest,
    record_failure, bootstrap_rho_ci, write_prereg_stage1, progress_snapshot,
    fibro_log_banner, dist_summary,
)
from gtex_stage0 import stream_gct_subset, read_tsv_log_columns  # noqa: E402
from gtex_stage1 import prepare_fold_X, predict_fold, permute_age_within_center  # noqa: E402
from target_common import ridge_prestd, age_direction, zscore_train, unit, align_age_sign  # noqa: E402
from trajectory_common import permutation_p, summarize_null_col  # noqa: E402
from gtex_common import file_md5  # noqa: E402


REGIMES = ("raw", "resid")
METHODS = ("ridge", "pls1")


def _require_prereg():
    if not PREREG_STAGE1_FLAG.exists():
        raise StopStep(
            "prereg",
            f"{PREREG_STAGE1_FLAG.name} missing — write Stage 1 pre-registration before any fit.",
        )


def _update_man(**kwargs):
    man = load_manifest()
    man.update(kwargs)
    save_manifest(man)
    return man


def build_cohort(log):
    samp_path = GTEX_RAW / FILES["sample_attributes"]["name"]
    pheno_path = GTEX_RAW / FILES["subject_phenotypes"]["name"]
    if not samp_path.exists():
        raise StopStep("cohort", f"missing sample attributes: {samp_path}")
    if not pheno_path.exists():
        raise StopStep("cohort", f"missing subject phenotypes: {pheno_path}")
    samp = read_tsv_log_columns(samp_path, "sample_attributes", log)
    pheno = read_tsv_log_columns(pheno_path, "subject_phenotypes", log)
    samp_cols = require_columns(samp, SAMPLE_REQUIRED, file_tag="sample_attributes", step="cohort")
    pheno_cols = require_columns(pheno, PHENO_REQUIRED, file_tag="subject_phenotypes", step="cohort")
    freeze_levels = sorted(samp["SMAFRZE"].dropna().astype(str).unique().tolist())
    log(f"[cohort] SMAFRZE levels observed: {freeze_levels}")
    if SMAFRZE_RNASEQ not in freeze_levels:
        raise StopStep("cohort", f"SMAFRZE does not contain {SMAFRZE_RNASEQ!r}. observed={freeze_levels}")
    smtsd_levels = sorted(samp["SMTSD"].dropna().astype(str).unique().tolist())
    if SMTSD_FIBRO not in set(smtsd_levels):
        nearby = [t for t in smtsd_levels if "fibro" in t.lower() or "skin" in t.lower() or "cell" in t.lower()]
        raise StopStep(
            "cohort",
            f"SMTSD missing exact {SMTSD_FIBRO!r}. Not substituting. Nearby labels: {nearby[:40]}",
            dict(n_smtsd=len(smtsd_levels)),
        )
    keep = samp["SMAFRZE"].astype(str).eq(SMAFRZE_RNASEQ) & samp["SMTSD"].astype(str).eq(SMTSD_FIBRO)
    sub = samp.loc[keep].copy()
    log(f"[filter] SMAFRZE=={SMAFRZE_RNASEQ} and SMTSD=={SMTSD_FIBRO!r}: n={int(keep.sum())}")
    if sub.empty:
        raise StopStep("cohort", "zero samples after freeze × tissue filter")
    sub["donor"] = sub["SAMPID"].map(sampid_to_subj)
    n_per = sub.groupby("donor").size()
    if int(n_per.max()) != 1:
        raise StopStep(
            "cohort",
            f"expected one sample per donor; max samples/donor={int(n_per.max())}. Not dropping extras.",
            dict(n_multi=int((n_per > 1).sum())),
        )
    pheno = pheno.copy()
    pheno["SUBJID"] = pheno["SUBJID"].astype(str)
    missing_d = sorted(set(sub["donor"].astype(str)) - set(pheno["SUBJID"].astype(str)))
    if missing_d:
        raise StopStep("cohort", f"{len(missing_d)} donors have no phenotype row. example={missing_d[:10]}")
    merged = sub.merge(pheno, left_on="donor", right_on="SUBJID", how="left", validate="many_to_one")
    if merged["AGE"].isna().any():
        raise StopStep("cohort", f"{int(merged['AGE'].isna().sum())} samples have missing AGE after merge")
    sex_info = require_levels(merged["SEX"], EXPECTED_SEX, field="SEX", step="cohort")
    dth_info = require_levels(merged["DTHHRDY"], EXPECTED_DTHHRDY, field="DTHHRDY", step="cohort")
    age_info = require_levels(merged["AGE"], frozenset(EXPECTED_AGE_BINS), field="AGE", step="cohort")
    log(f"[cohort] SEX levels: {sex_info}")
    log(f"[cohort] DTHHRDY levels: {dth_info}")
    log(f"[cohort] AGE bins: {age_info}")
    merged["age_mid"] = merged["AGE"].map(lambda a: age_midpoint(str(a)))
    merged["age_ordinal"] = merged["AGE"].map({b: i for i, b in enumerate(EXPECTED_AGE_BINS)})
    merged["SMRIN"] = pd.to_numeric(merged["SMRIN"], errors="coerce")
    merged["SMTSISCH"] = pd.to_numeric(merged["SMTSISCH"], errors="coerce")
    for col, raw in (("SMRIN", sub["SMRIN"]), ("SMTSISCH", sub["SMTSISCH"])):
        raw_s = raw.astype("string")
        blank = raw_s.isna() | raw_s.str.strip().isin(["", "NA", "NaN", "."])
        bad = (~blank.to_numpy()) & merged[col].isna().to_numpy()
        if bool(np.any(bad)):
            examples = raw_s.to_numpy()[bad][:10].tolist()
            raise StopStep("cohort", f"{col}: non-numeric values cannot be coerced. examples={examples}")
    centers = merged["SMCENTER"].astype(str)
    vc = centers.value_counts(dropna=False)
    log(f"[cohort] SMCENTER counts (file strings, not recoded): {vc.to_dict()}")
    observed = sorted(centers.dropna().unique().tolist())
    expected = set(TRANSFER_SITES + (EXCLUDE_SITE,))
    extra = sorted(set(observed) - expected)
    if extra:
        raise StopStep(
            "cohort",
            f"SMCENTER unexpected levels {extra}. observed={observed}. "
            "Not recoding comma-joined values or substituting.",
            dict(observed=observed),
        )
    missing_exp = sorted(expected - set(observed))
    if missing_exp:
        raise StopStep("cohort", f"SMCENTER missing expected levels {missing_exp}. observed={observed}")
    n_d1 = int((centers == EXCLUDE_SITE).sum())
    log(f"[cohort] {EXCLUDE_SITE} n={n_d1}. {EXCLUDE_SITE_REASON}")
    rin = dist_summary(merged.SMRIN, numeric=True)
    isch = dist_summary(merged.SMTSISCH, numeric=True)
    age_vc = merged.AGE.astype(str).value_counts().reindex(list(EXPECTED_AGE_BINS)).fillna(0).astype(int)
    rec = dict(
        tissue=SMTSD_FIBRO, n_samples=int(len(merged)), n_donors=int(merged.donor.nunique()),
        n_SMCENTER=int(merged.SMCENTER.nunique(dropna=True)),
        n_SMNABTCH=int(merged.SMNABTCH.nunique(dropna=True)),
        n_SMGEBTCH=int(merged.SMGEBTCH.nunique(dropna=True)),
        n_male=int((merged.SEX.astype(str) == "1").sum()),
        n_female=int((merged.SEX.astype(str) == "2").sum()),
        SMRIN_n_missing=rin["n_missing"], SMRIN_median=rin["median"],
        SMRIN_min=rin["min"], SMRIN_max=rin["max"],
        SMTSISCH_n_missing=isch["n_missing"], SMTSISCH_median=isch["median"],
        SMTSISCH_min=isch["min"], SMTSISCH_max=isch["max"],
        DTHHRDY_n_missing=int(merged.DTHHRDY.isna().sum()),
        SMCENTER_counts=str({str(k): int(v) for k, v in vc.items()}),
        D1_n=n_d1, D1_excluded_from_transfer=True, D1_reason=EXCLUDE_SITE_REASON,
    )
    for b in EXPECTED_AGE_BINS:
        rec[f"n_AGE_{b}"] = int(age_vc.get(b, 0))
    pd.DataFrame([rec]).to_csv(FIBRO_DIR / "stage1_cohort.csv", index=False)
    merged.to_csv(FIBRO_DIR / "stage1_obs.csv", index=False)
    dump_json(FIBRO_DIR / "stage1_cohort.json", jsonable(dict(
        rec=rec, sex=sex_info, dthhrdy=dth_info, age=age_info,
        SMCENTER_counts={str(k): int(v) for k, v in vc.items()},
        sample_columns_read=samp_cols, phenotype_columns_read=pheno_cols,
        SMAFRZE_levels=freeze_levels, caveat_age_bin=CAVEAT_AGE_BIN,
    )))
    man = load_manifest()
    man["columns"] = dict(sample_attributes=samp_cols, subject_phenotypes=pheno_cols)
    man["cohort"] = jsonable(rec)
    man["status"] = "COHORT_DONE"
    save_manifest(man)
    log(f"[write] stage1_obs.csv n={len(merged)} donors={merged.donor.nunique()}")
    return merged


def preprocess(log):
    obs_path = FIBRO_DIR / "stage1_obs.csv"
    if not obs_path.exists():
        raise StopStep("preprocess", "stage1_obs.csv missing — run cohort")
    obs = pd.read_csv(obs_path)
    gct_path = GTEX_RAW / FILES["gene_reads"]["name"]
    if not gct_path.exists():
        raise StopStep("preprocess", f"missing GCT: {gct_path}")
    keep_ids = set(obs.SAMPID.astype(str))
    counts, genes, symbols, samples, gct_meta = stream_gct_subset(gct_path, keep_ids, log)
    dump_json(FIBRO_DIR / "stage1_gct_meta.json", jsonable(
        {k: v for k, v in gct_meta.items() if k != "columns_read"}
    ))
    (FIBRO_DIR / "stage1_gct_columns_read.txt").write_text(
        "Name\nDescription\n" + "\n".join(map(str, gct_meta.get("kept_sample_ids", []))) + "\n",
        encoding="utf-8",
    )
    pos = {str(s): i for i, s in enumerate(samples)}
    ids = obs.SAMPID.astype(str).tolist()
    missing = [s for s in ids if s not in pos]
    if missing:
        raise StopStep("preprocess", f"{len(missing)} SAMPID not in GCT subset. example={missing[:10]}")
    col = np.array([pos[s] for s in ids], dtype=int)
    C = counts[:, col]
    keep = gene_filter_mask(C)
    n0, n1 = int(C.shape[0]), int(keep.sum())
    log(f"[filter] genes ≥{MIN_COUNT} counts in ≥{MIN_FRAC:.0%} samples  {n0:,} -> {n1:,}  "
        f"(removed {n0 - n1:,})  n_samples={C.shape[1]}")
    Cf = C[keep]
    g_keep = genes[keep]
    s_keep = symbols[keep]
    nf = tmm_norm_factors(Cf.astype(np.float64))
    logcpm = log2_cpm_edger(Cf.astype(np.float64), nf)
    ens = np.array([strip_ensembl(g) for g in g_keep])
    universe, dlpfc_genes = load_dlpfc_gene_universe()
    n_overlap = int(sum(1 for e in ens if e in universe))
    ov = dict(
        n_samples=int(C.shape[1]), n_genes_raw=n0, n_genes_filter=n1,
        n_overlap_dlpfc=n_overlap, frac_overlap=n_overlap / max(n1, 1),
        n_dlpfc=len(universe), dlpfc_path=str(DLPFC_GENES_CSV),
        dlpfc_columns_read=list(dlpfc_genes.columns),
        gene_filter=f">={MIN_COUNT} counts in >={MIN_FRAC:.0%} of samples",
        normalization=f"TMM then log2-CPM prior.count={EDGE_R_PRIOR:g}",
    )
    dump_json(FIBRO_DIR / "stage1_gene_overlap.json", jsonable(ov))
    log(f"[overlap] filtered genes {n1:,} ∩ DLPFC universe {n_overlap:,} / {len(universe):,}")
    out = FIBRO_PROC / "stage1_fibro.npz"
    np.savez_compressed(
        out,
        counts=Cf.astype(np.float32),
        logcpm=logcpm.astype(np.float32),
        tmm_factors=nf.astype(np.float64),
        sample_ids=np.array(ids, dtype=object),
        genes=np.asarray(g_keep, dtype=object),
        symbols=np.asarray(s_keep, dtype=object),
        ensembl=ens.astype(object),
        tissue=np.array(SMTSD_FIBRO),
    )
    log(f"[write] {out} genes={n1} samples={C.shape[1]}")
    man = load_manifest()
    man["gene_overlap"] = jsonable(ov)
    man["processed"] = str(out)
    man["status"] = "PREPROCESS_DONE"
    man["gct"] = dict(
        n_genes_declared=gct_meta.get("n_genes_declared"),
        n_samples_declared=gct_meta.get("n_samples_declared"),
        n_kept=gct_meta.get("n_kept"),
        gct_version=gct_meta.get("gct_version"),
        columns_fixed=["Name", "Description"],
        n_sample_columns_in_header=gct_meta.get("n_sample_columns_in_header"),
    )
    save_manifest(man)
    return ov


def load_pack():
    path = FIBRO_PROC / "stage1_fibro.npz"
    obs_path = FIBRO_DIR / "stage1_obs.csv"
    if not path.exists():
        raise StopStep("stage1_data", f"missing {path} — run preprocess")
    if not obs_path.exists():
        raise StopStep("stage1_data", f"missing {obs_path}")
    z = np.load(path, allow_pickle=True)
    obs = pd.read_csv(obs_path)
    ids = z["sample_ids"].astype(str)
    pos = {s: i for i, s in enumerate(obs.SAMPID.astype(str))}
    missing = [s for s in ids if s not in pos]
    if missing:
        raise StopStep("stage1_data", f"obs missing {len(missing)} sample ids")
    obs = obs.iloc[[pos[s] for s in ids]].reset_index(drop=True)
    if list(obs.SAMPID.astype(str)) != list(ids):
        raise StopStep("stage1_data", "obs/matrix SAMPID order mismatch")
    X = np.asarray(z["logcpm"], np.float64).T
    return dict(
        X=X, obs=obs, genes=z["genes"], symbols=z["symbols"],
        ensembl=z["ensembl"].astype(str), counts=z["counts"],
        tmm_factors=z["tmm_factors"],
    )


def _site_folds(pack, site, log):
    obs, X = pack["obs"], pack["X"]
    m = np.flatnonzero(obs.SMCENTER.astype(str).to_numpy() == site)
    n = int(len(m))
    if n < SITE_CV_MIN_N:
        log(f"[cv] skip site {site}: n={n} < {SITE_CV_MIN_N}")
        return None, dict(skipped=True, site=site, n=n, reason=f"n<{SITE_CV_MIN_N}")
    obs_s = obs.iloc[m].reset_index(drop=True)
    X_s = X[m]
    batches = obs_s.SMNABTCH.astype(str).to_numpy()
    donors = obs_s.donor.astype(str).to_numpy()
    n_lv = int(pd.Series(batches).nunique())
    if n_lv < 2:
        raise StopStep("batch_held_out", f"site {site}: SMNABTCH has {n_lv} level(s); cannot hold out batch")
    kf = grouped_kfold(batches, n_splits=min(N_FOLDS, n_lv))
    folds = []
    for fi, (tr, te, te_lv) in enumerate(kf):
        assert_disjoint(donors[tr], donors[te], what="donor",
                        where=f"within-site {site} SMNABTCH fold{fi}")
        assert_disjoint(batches[tr], batches[te], what="SMNABTCH",
                        where=f"within-site {site} SMNABTCH fold{fi}")
        folds.append((tr, te, te_lv))
    return dict(X=X_s, obs=obs_s, folds=folds, site=site, n=n, n_batches=n_lv), None


def _cv_one(site_pack, regime, method, log):
    X, obs, folds = site_pack["X"], site_pack["obs"], site_pack["folds"]
    cache = []
    for fi, (tr, te, te_lv) in enumerate(folds):
        Xtr, Xte, ytr, yte, meta = prepare_fold_X(X, obs, tr, te, regime, method)
        svd = ridge_prestd(Xtr, np.zeros(len(tr), dtype=float))[3] if method == "ridge" else None
        cache.append(dict(tr=tr, te=te, Xtr=Xtr, Xte=Xte, svd=svd, test_batches=",".join(te_lv)))
    y0 = obs.age_mid.to_numpy(float)
    pred = np.full(len(obs), np.nan)
    for f in cache:
        p, w, rec, _ = predict_fold(f["Xtr"], y0[f["tr"]], f["Xte"], method, svd=f["svd"])
        pred[f["te"]] = p
    sc = pred_scores(y0, pred)
    rng_p = np.random.default_rng(FIBRO_SEED)
    null_rho, null_r, null_r2, null_cal = [], [], [], []
    for i in range(N_PERM):
        yp = permute_age_within_center(obs, rng_p)
        pred_p = np.full(len(obs), np.nan)
        for f in cache:
            p, _, _, _ = predict_fold(f["Xtr"], yp[f["tr"]], f["Xte"], method, svd=f["svd"])
            pred_p[f["te"]] = p
        scp = pred_scores(yp, pred_p)
        null_rho.append(scp["rho"])
        null_r.append(scp["r"])
        null_r2.append(scp["r2"])
        null_cal.append(scp["cal_r2"])
        if (i + 1) % 50 == 0:
            log(f"   {site_pack['site']} {regime} {method} perm {i+1}/{N_PERM}")
    rng_b = np.random.default_rng(BOOT_SEED)
    ci = bootstrap_rho_ci(y0, pred, rng_b, n_boot=N_BOOT)
    rho_null = float(np.nanmedian(null_rho))
    return dict(
        site=site_pack["site"], regime=regime, method=method, cv="SMNABTCH_grouped",
        rho=sc["rho"], rho_null=rho_null,
        rho_p=permutation_p(sc["rho"], np.asarray(null_rho, float), greater=True),
        rho_null_mean=float(np.nanmean(null_rho)),
        rho_ci_lo=ci["p025"], rho_ci_hi=ci["p975"],
        r=sc["r"], r_null=float(np.nanmedian(null_r)),
        r2=sc["r2"], r2_null=float(np.nanmedian(null_r2)),
        cal_r2=sc["cal_r2"], cal_r2_null=float(np.nanmedian(null_cal)),
        n=int(len(obs)), n_folds=len(folds), n_batches=site_pack["n_batches"],
        n_perm=N_PERM, n_boot=N_BOOT,
        pass_rho=pass_rho_bar(sc["rho"], rho_null),
        batch_held_out=True, donor_held_out=True,
    )


def run_cv(log):
    _require_prereg()
    pack = load_pack()
    rows, skipped = [], []
    t0 = time.time()
    for site in TRANSFER_SITES + (EXCLUDE_SITE,):
        site_pack, skip = _site_folds(pack, site, log)
        if skip is not None:
            skipped.append(skip)
            continue
        for regime in REGIMES:
            for method in METHODS:
                log(f"[cv] site={site} regime={regime} method={method} n={site_pack['n']}")
                rows.append(_cv_one(site_pack, regime, method, log))
    df = pd.DataFrame(rows)
    df.to_csv(FIBRO_DIR / "stage1_cv.csv", index=False)
    dump_json(FIBRO_DIR / "stage1_cv.json", jsonable(dict(rows=rows, skipped=skipped)))
    man = load_manifest()
    man["cv_skipped"] = jsonable(skipped)
    man["status"] = "CV_DONE"
    save_manifest(man)
    log(f"[cv] wrote {len(rows)} rows in {time.time()-t0:.0f}s  skipped={skipped}")
    return rows


def _transfer_one(pack, train_site, test_site, regime, method, log):
    obs, X = pack["obs"], pack["X"]
    cvec = obs.SMCENTER.astype(str).to_numpy()
    donors = obs.donor.astype(str).to_numpy()
    tr = np.flatnonzero(cvec == train_site)
    te = np.flatnonzero(cvec == test_site)
    if len(tr) < TRANSFER_MIN_N or len(te) < TRANSFER_MIN_N:
        raise StopStep(
            "transfer",
            f"{train_site}→{test_site}: n_train={len(tr)} n_test={len(te)}; "
            f"bar is n≥{TRANSFER_MIN_N}. Not lowering the bar.",
        )
    assert_disjoint(donors[tr], donors[te], what="donor",
                    where=f"transfer {train_site}->{test_site} {regime} {method}")
    log(f"[transfer] {train_site}→{test_site} {regime} {method} n_tr={len(tr)} n_te={len(te)}")
    Xtr, Xte, ytr, yte, meta = prepare_fold_X(X, obs, tr, te, regime, method)
    svd = ridge_prestd(Xtr, ytr)[3] if method == "ridge" else None
    pred, w, rec, _ = predict_fold(Xtr, ytr, Xte, method, svd=svd)
    sc = pred_scores(yte, pred)
    obs_tr = obs.iloc[tr].reset_index(drop=True)
    rng_p = np.random.default_rng(FIBRO_SEED)
    nulls = {k: [] for k in ("r", "r2", "rho", "cal_r2")}
    for i in range(N_PERM):
        ytr_p = permute_age_within_center(obs_tr, rng_p)
        pred_p, _, _, _ = predict_fold(Xtr, ytr_p, Xte, method, svd=svd)
        scp = pred_scores(yte, pred_p)
        for k in nulls:
            nulls[k].append(scp[k])
        if (i + 1) % 50 == 0:
            log(f"   perm {i+1}/{N_PERM} {train_site}→{test_site} {regime} {method}")
    rng_b = np.random.default_rng(BOOT_SEED)
    ci = bootstrap_rho_ci(yte, pred, rng_b, n_boot=N_BOOT)
    rho_null = float(np.nanmedian(nulls["rho"]))
    return dict(
        direction=f"{train_site}→{test_site}",
        train_site=train_site, test_site=test_site,
        regime=regime, method=method,
        rho=sc["rho"], rho_null=rho_null,
        rho_p=permutation_p(sc["rho"], np.asarray(nulls["rho"], float), greater=True),
        rho_null_mean=float(np.nanmean(nulls["rho"])),
        rho_ci_lo=ci["p025"], rho_ci_hi=ci["p975"],
        r=sc["r"], r_null=float(np.nanmedian(nulls["r"])),
        r2=sc["r2"], r2_null=float(np.nanmedian(nulls["r2"])),
        cal_r2=sc["cal_r2"], cal_r2_null=float(np.nanmedian(nulls["cal_r2"])),
        pass_rho=pass_rho_bar(sc["rho"], rho_null),
        n_train=int(len(tr)), n_test=int(len(te)),
        n_perm=N_PERM, n_boot=N_BOOT, n_genes=int(X.shape[1]),
        excluded_site=EXCLUDE_SITE, excluded_reason=EXCLUDE_SITE_REASON,
    )


def run_transfer(log):
    _require_prereg()
    pack = load_pack()
    n_excl = int((pack["obs"].SMCENTER.astype(str) == EXCLUDE_SITE).sum())
    log(f"[transfer] excluding {EXCLUDE_SITE} n={n_excl}. {EXCLUDE_SITE_REASON}")
    rows = []
    t0 = time.time()
    for regime in REGIMES:
        for method in METHODS:
            rows.append(_transfer_one(pack, "B1", "C1", regime, method, log))
            rows.append(_transfer_one(pack, "C1", "B1", regime, method, log))
    pd.DataFrame(rows).to_csv(FIBRO_DIR / "stage1_transfer.csv", index=False)
    dump_json(FIBRO_DIR / "stage1_transfer.json", jsonable(rows))
    man = load_manifest()
    man["transfer_excluded"] = dict(site=EXCLUDE_SITE, n=n_excl, reason=EXCLUDE_SITE_REASON)
    man["status"] = "TRANSFER_DONE"
    save_manifest(man)
    log(f"[transfer] wrote {len(rows)} rows in {time.time()-t0:.0f}s")
    return rows


def run_gate(log):
    path = FIBRO_DIR / "stage1_transfer.csv"
    if not path.exists():
        raise StopStep("gate", "stage1_transfer.csv missing — run transfer")
    df = pd.read_csv(path)
    sub = df[(df.regime == "raw") & (df.method == "ridge")].copy()
    if len(sub) != 2:
        raise StopStep("gate", f"expected 2 raw ridge transfer rows, got {len(sub)}")
    by = {str(r.direction): r for _, r in sub.iterrows()}
    need = ("B1→C1", "C1→B1")
    missing = [d for d in need if d not in by]
    if missing:
        raise StopStep("gate", f"missing transfer directions {missing}. have={list(by)}")
    a, b = by["B1→C1"], by["C1→B1"]
    pass_a = pass_rho_bar(a.rho, a.rho_null)
    pass_b = pass_rho_bar(b.rho, b.rho_null)
    passed = bool(pass_a and pass_b)
    if passed:
        reason = (
            f"Stage 1 pass: ridge raw ρ>0 both ways with null ≤ {TRANSFER_NULL_BAR} "
            f"(B1→C1 ρ={float(a.rho):+.3f} null={float(a.rho_null):+.3f} p={float(a.rho_p):+.3f}; "
            f"C1→B1 ρ={float(b.rho):+.3f} null={float(b.rho_null):+.3f} p={float(b.rho_p):+.3f})."
        )
    else:
        reason = (
            "Stage 1 fail: ridge raw transfer did not pass both directions. "
            "The fibroblast ruler does not transfer between GTEx sites so the "
            "reprogramming projection would be uninterpretable. GSE297234 not opened."
            f" B1→C1 ρ={float(a.rho):+.3f} null={float(a.rho_null):+.3f} pass={pass_a};"
            f" C1→B1 ρ={float(b.rho):+.3f} null={float(b.rho_null):+.3f} pass={pass_b}."
        )
    rec = {
        "pass": passed, "method": "ridge", "regime": "raw",
        "rho_B1_to_C1": float(a.rho), "rho_null_B1_to_C1": float(a.rho_null),
        "rho_p_B1_to_C1": float(a.rho_p), "pass_B1_to_C1": pass_a,
        "rho_ci_B1_to_C1": [float(a.rho_ci_lo), float(a.rho_ci_hi)],
        "rho_C1_to_B1": float(b.rho), "rho_null_C1_to_B1": float(b.rho_null),
        "rho_p_C1_to_B1": float(b.rho_p), "pass_C1_to_B1": pass_b,
        "rho_ci_C1_to_B1": [float(b.rho_ci_lo), float(b.rho_ci_hi)],
        "r2_B1_to_C1": float(a.r2), "r2_C1_to_B1": float(b.r2),
        "cal_r2_B1_to_C1": float(a.cal_r2), "cal_r2_C1_to_B1": float(b.cal_r2),
        "n_perm": N_PERM, "n_boot": N_BOOT, "seed": FIBRO_SEED, "boot_seed": BOOT_SEED,
        "excluded_site": EXCLUDE_SITE, "reason": reason,
        "pls1_not_used_for_gate": True,
        "uncalibrated_R2_is_not_a_gate": True,
    }
    dump_json(FIBRO_DIR / "stage1_gate.json", jsonable(rec))
    (FIBRO_DIR / "stage1_gate.txt").write_text(reason + "\n", encoding="utf-8")
    man = load_manifest()
    man["stage1_gate"] = jsonable(rec)
    man["stage2_opened"] = False
    man["status"] = "STAGE1_PASS" if passed else "STAGE1_FAIL"
    save_manifest(man)
    log(f"[gate] pass={passed}  {reason}")
    dump_json(FIBRO_DIR / "stage1_summary.json", jsonable(dict(
        seed=FIBRO_SEED, boot_seed=BOOT_SEED, n_perm=N_PERM, n_boot=N_BOOT,
        prereg=str(PREREG_STAGE1_FLAG), gate=rec,
        caveat_age_bin=CAVEAT_AGE_BIN,
    )))
    return rec


def freeze_ruler(log):
    """Fit raw ridge on ALL GTEx fibroblast donors. Called only after Stage 1 pass,
    and only after PREREG_STAGE2.flag exists (Stage 2 owns that flag)."""
    gate_p = FIBRO_DIR / "stage1_gate.json"
    if not gate_p.exists():
        raise StopStep("freeze", "stage1_gate.json missing")
    gate = load_json(gate_p)
    if not gate.get("pass"):
        raise StopStep("freeze", "Stage 1 did not pass; refusing to freeze a non-transferring ruler")
    from fibro_common import PREREG_STAGE2_FLAG
    if not PREREG_STAGE2_FLAG.exists():
        raise StopStep("freeze", "PREREG_STAGE2.flag missing — write it before freezing the ruler")
    pack = load_pack()
    X, obs = pack["X"], pack["obs"]
    y = obs.age_mid.to_numpy(float)
    log(f"[freeze] raw ridge on all donors n={len(obs)} genes={X.shape[1]} "
        f"(includes {EXCLUDE_SITE} n={int((obs.SMCENTER.astype(str)==EXCLUDE_SITE).sum())})")
    Xz_all, mu, sd = zscore_train(X)
    w, rec = age_direction(Xz_all, y, method="ridge")
    if not rec.get("ok"):
        raise StopStep("freeze", f"ridge failed on all donors: {rec}")
    w, _ = unit(w)
    w = align_age_sign(w, Xz_all, y)
    w, _ = unit(w)
    coef = rec.get("coef")
    intercept = rec.get("intercept")
    alpha = rec.get("alpha")
    np.savez_compressed(
        FROZEN_RULER,
        w=np.asarray(w, np.float64),
        coef=np.asarray(coef if coef is not None else w, np.float64),
        mu=np.asarray(mu, np.float64),
        sd=np.asarray(sd, np.float64),
        intercept=np.array(float(intercept) if intercept is not None else np.nan),
        alpha=np.array(float(alpha) if alpha is not None else np.nan),
        gene_id=pack["genes"].astype(str),
        symbol=pack["symbols"].astype(str),
        ensembl=pack["ensembl"].astype(str),
        sample_ids=obs.SAMPID.astype(str).to_numpy(),
        donors=obs.donor.astype(str).to_numpy(),
        age_mid=y,
        seed=np.array(FIBRO_SEED),
        n_donors=np.array(int(obs.donor.nunique())),
        n_genes=np.array(int(X.shape[1])),
        regime=np.array("raw"),
        method=np.array("ridge"),
        note=np.array("trained on all GTEx fibroblast donors; not fit on GSE297234"),
    )
    rec_out = dict(
        path=str(FROZEN_RULER), n_donors=int(obs.donor.nunique()), n_genes=int(X.shape[1]),
        alpha=float(alpha) if alpha is not None else None,
        intercept=float(intercept) if intercept is not None else None,
        n_finite_w=int(np.isfinite(w).sum()),
        seed=FIBRO_SEED, regime="raw", method="ridge",
        includes_D1=True,
        note="trained on all GTEx fibroblast donors; not fit on GSE297234",
    )
    dump_json(FIBRO_DIR / "freeze.json", jsonable(rec_out))
    log(f"[freeze] wrote {FROZEN_RULER}  n_donors={rec_out['n_donors']} n_genes={rec_out['n_genes']} "
        f"alpha={rec_out['alpha']}")
    man = load_manifest()
    man["frozen_ruler"] = jsonable(rec_out)
    save_manifest(man)
    return rec_out


def run_cell(name, log):
    t0 = time.time()
    log(f"===== {name} start =====")
    if name == "cohort":
        write_prereg_stage1()
        build_cohort(log)
        progress_snapshot("preprocess GCT subset + TMM/log2-CPM", extra=f"cohort done in {time.time()-t0:.0f}s")
    elif name == "preprocess":
        write_prereg_stage1()
        preprocess(log)
        progress_snapshot("within-site SMNABTCH CV", extra=f"preprocess done in {time.time()-t0:.0f}s")
    elif name == "cv":
        run_cv(log)
        progress_snapshot("B1↔C1 transfer", extra=f"cv done in {time.time()-t0:.0f}s")
    elif name == "transfer":
        run_transfer(log)
        progress_snapshot("Stage 1 gate (ridge raw both directions)", extra=f"transfer done in {time.time()-t0:.0f}s")
    elif name == "gate":
        rec = run_gate(log)
        if rec["pass"]:
            progress_snapshot(
                "write PREREG_STAGE2.flag then freeze ruler then GSE297234",
                stop="STAGE1_PASS", extra=f"gate done in {time.time()-t0:.0f}s",
            )
        else:
            progress_snapshot(
                "write FINDINGS_FIBRO.md Stage 1 fail; do not open GSE297234",
                stop="STAGE1_FAIL", extra=f"gate done in {time.time()-t0:.0f}s",
            )
    elif name == "freeze":
        freeze_ruler(log)
        progress_snapshot("project GSE297234 (Stage 2)", extra=f"freeze done in {time.time()-t0:.0f}s")
    else:
        raise StopStep("cell", f"unknown cell {name}")
    log(f"===== {name} done {time.time()-t0:.0f}s =====")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", default="all",
                    help="cohort preprocess cv transfer gate freeze all")
    args, _unknown = ap.parse_known_args()
    write_prereg_stage1()
    log = Logger(FIBRO_DIR / "stage1_report.txt")
    fibro_log_banner(log, "STAGE 1")
    log(f"[prereg] wrote/exists {PREREG_STAGE1_FLAG}")
    log(CAVEAT_AGE_BIN)
    try:
        order = ["cohort", "preprocess", "cv", "transfer", "gate"]
        cells = order if args.cell == "all" else [args.cell]
        if args.cell == "all":
            progress_snapshot("cohort from GTEx annotations", stop="Stage 1 starting")
        for c in cells:
            run_cell(c, log)
        if args.cell == "all":
            log("STAGE 1 complete (gate scored). Freeze/Stage 2 only if pass.")
    except StopStep as e:
        log(f"STOP [{e.step}] {e.message}")
        record_failure(e.step, e.message, e.details)
        dump_json(FIBRO_DIR / "stage1_STOP.json", dict(step=e.step, message=e.message, details=e.details))
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise
    except Exception as e:
        log(f"FAIL {type(e).__name__}: {e}")
        record_failure("fail", f"{type(e).__name__}: {e}")
        progress_snapshot(f"fix {type(e).__name__}", stop=f"FAIL {type(e).__name__}: {e}")
        raise
    finally:
        log.close()


if __name__ == "__main__":
    main()

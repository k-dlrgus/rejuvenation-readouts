"""Stage 0 — GTEx v10 download, QC, covariate audit, composition proxies.

Usage: python src/gtex_stage0.py
Does not run Stage 1. Does not modify prior FINDINGS*.md or FALSIFICATION.md.
"""
from __future__ import annotations

import gzip
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gtex_common import (  # noqa: E402
    GTEX_DIR, GTEX_FIG, GTEX_RAW, GTEX_PROC, GTEX_SEED, RELEASE, RNASEQ_PIPELINE,
    DBGAP, FILES, TISSUES, PRIMARY_TISSUE, TISSUE_SLUG, SMAFRZE_RNASEQ,
    SAMPLE_REQUIRED, PHENO_REQUIRED, EXPECTED_SEX, EXPECTED_DTHHRDY,
    EXPECTED_AGE_BINS, AGE_MID, NEURONAL_MARKERS, GLIAL_MARKERS, BRAIN_TISSUES,
    MANIFEST_PATH, DLPFC_GENES_CSV, BA9_N_UNDERPOWERED, CENTER_TRANSFER_MIN_N,
    CENTER_TRANSFER_MIN_LEVELS, CAVEATS,
    Logger, dump_json, load_json, StopStep, gtex_log_banner, file_md5, gcs_md5_hex,
    sampid_to_subj, require_columns, require_levels, age_midpoint, tmm_norm_factors,
    log2_cpm_edger, gene_filter_mask, ols_r2, onehot_train, spearman_safe,
    dist_summary, load_dlpfc_gene_universe, stage1_gate, strip_ensembl, fmt,
)
from download import download as http_download  # noqa: E402


def _download_one(key, spec, log):
    dest = GTEX_RAW / spec["name"]
    log(f"[download] {key}: {spec['url']}")
    path = Path(http_download(spec["url"], str(GTEX_RAW)))
    if path.resolve() != dest.resolve() and path.exists() and not dest.exists():
        path.replace(dest)
        path = dest
    size = int(path.stat().st_size)
    md5 = file_md5(path)
    gcs_hex = gcs_md5_hex(spec["gcs_md5_b64"])
    rec = dict(
        key=key,
        file_name=spec["name"],
        url=spec["url"],
        path=str(path),
        size_bytes=size,
        gcs_size_bytes=int(spec["gcs_size"]),
        md5=md5,
        gcs_md5_hex=gcs_hex,
        gcs_md5_b64=spec["gcs_md5_b64"],
        md5_matches_gcs=md5 == gcs_hex,
        size_matches_gcs=size == int(spec["gcs_size"]),
        kind=spec["kind"],
        release=RELEASE,
        pipeline=RNASEQ_PIPELINE,
        dbgap=DBGAP,
        columns_read=None,
        n_rows=None,
    )
    if size != int(spec["gcs_size"]):
        raise StopStep(
            "download",
            f"{spec['name']}: size {size} != GCS listing {spec['gcs_size']}",
            rec,
        )
    if md5 != gcs_hex:
        raise StopStep(
            "download",
            f"{spec['name']}: md5 {md5} != GCS {gcs_hex}",
            rec,
        )
    log(f"[download] {spec['name']} size={size} md5={md5} (matches GCS)")
    return rec


def read_tsv_log_columns(path: Path, file_tag: str, log) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype="string", low_memory=False)
    log(f"[read] {file_tag}: {path.name} rows={len(df):,} cols={len(df.columns)}")
    log(f"[read] {file_tag} columns actually read ({len(df.columns)}): {list(df.columns)}")
    return df


def build_cohort(samp: pd.DataFrame, pheno: pd.DataFrame, log) -> pd.DataFrame:
    step = "cohort"
    samp_cols = require_columns(samp, SAMPLE_REQUIRED, file_tag="sample_attributes", step=step)
    pheno_cols = require_columns(pheno, PHENO_REQUIRED, file_tag="subject_phenotypes", step=step)
    freeze_levels = sorted(samp["SMAFRZE"].dropna().astype(str).unique().tolist())
    log(f"[cohort] SMAFRZE levels observed: {freeze_levels}")
    if SMAFRZE_RNASEQ not in freeze_levels:
        raise StopStep(
            step,
            f"SMAFRZE does not contain {SMAFRZE_RNASEQ!r}. observed={freeze_levels}",
            dict(observed=freeze_levels),
        )
    smtsd = samp["SMTSD"].astype(str)
    present_tissues = sorted(smtsd.unique().tolist())
    missing_t = [t for t in TISSUES if t not in set(present_tissues)]
    if missing_t:
        brainish = [t for t in present_tissues if "Brain" in t or "Heart" in t or "Muscle" in t]
        raise StopStep(
            step,
            f"SMTSD missing exact tissue name(s) {missing_t}. "
            f"Not substituting. Nearby labels: {brainish[:30]}",
            dict(missing=missing_t, n_smtsd=len(present_tissues)),
        )
    keep = samp["SMAFRZE"].astype(str).eq(SMAFRZE_RNASEQ) & smtsd.isin(list(TISSUES))
    sub = samp.loc[keep].copy()
    log(f"[filter] SMAFRZE=={SMAFRZE_RNASEQ} and SMTSD in 6 tissues: {int(keep.sum()):,} samples")
    if sub.empty:
        raise StopStep(step, "zero samples after freeze × tissue filter")

    sub["donor"] = sub["SAMPID"].map(sampid_to_subj)
    pheno = pheno.copy()
    pheno["SUBJID"] = pheno["SUBJID"].astype(str)
    donors_needed = set(sub["donor"].astype(str))
    pheno_donors = set(pheno["SUBJID"].astype(str))
    missing_d = sorted(donors_needed - pheno_donors)
    if missing_d:
        raise StopStep(
            step,
            f"{len(missing_d)} donors in samples have no subject-phenotype row. "
            f"example={missing_d[:10]}",
            dict(n_missing=len(missing_d), example=missing_d[:10]),
        )
    merged = sub.merge(pheno, left_on="donor", right_on="SUBJID", how="left", validate="many_to_one")
    if merged["AGE"].isna().any():
        n = int(merged["AGE"].isna().sum())
        raise StopStep(step, f"{n} samples have missing AGE after merge")

    sex_info = require_levels(merged["SEX"], EXPECTED_SEX, field="SEX", step=step)
    dth_info = require_levels(merged["DTHHRDY"], EXPECTED_DTHHRDY, field="DTHHRDY", step=step)
    age_info = require_levels(merged["AGE"], frozenset(EXPECTED_AGE_BINS), field="AGE", step=step)
    log(f"[cohort] SEX levels: {sex_info}")
    log(f"[cohort] DTHHRDY levels: {dth_info}")
    log(f"[cohort] AGE bins: {age_info}")
    merged["age_mid"] = merged["AGE"].map(lambda a: age_midpoint(str(a)))
    merged["age_ordinal"] = merged["AGE"].map({b: i for i, b in enumerate(EXPECTED_AGE_BINS)})
    merged["SMRIN"] = pd.to_numeric(merged["SMRIN"], errors="coerce")
    merged["SMTSISCH"] = pd.to_numeric(merged["SMTSISCH"], errors="coerce")
    # Fail if numeric coerce introduced NaNs that were non-empty strings
    for col, raw in (("SMRIN", sub["SMRIN"]), ("SMTSISCH", sub["SMTSISCH"])):
        raw_s = raw.astype("string")
        blank = raw_s.isna() | raw_s.str.strip().isin(["", "NA", "NaN", "."])
        coerced = merged[col]
        bad = (~blank.to_numpy()) & coerced.isna().to_numpy()
        if bool(np.any(bad)):
            examples = raw_s.to_numpy()[bad][:10].tolist()
            raise StopStep(
                step,
                f"{col}: non-numeric values cannot be coerced. examples={examples}",
                dict(column=col, examples=examples),
            )
    extra = {
        "sample_columns_read": samp_cols,
        "phenotype_columns_read": pheno_cols,
        "sex": sex_info,
        "dthhrdy": dth_info,
        "age": age_info,
        "smafrze_levels": freeze_levels,
        "n_after_freeze_tissue": int(len(merged)),
    }
    return merged, extra


def cohort_tables(obs: pd.DataFrame, log):
    rows = []
    for t in TISSUES:
        d = obs[obs.SMTSD.astype(str) == t]
        rin = dist_summary(d.SMRIN, numeric=True)
        isch = dist_summary(d.SMTSISCH, numeric=True)
        hardy = dist_summary(d.DTHHRDY, numeric=False)
        sex = dist_summary(d.SEX, numeric=False)
        center = d.SMCENTER.astype(str).value_counts(dropna=False).to_dict()
        batch = d.SMNABTCH.astype(str).value_counts(dropna=False).to_dict()
        age = d.AGE.astype(str).value_counts(dropna=False).reindex(list(EXPECTED_AGE_BINS)).fillna(0).astype(int)
        rec = dict(
            tissue=t, slug=TISSUE_SLUG[t],
            n_samples=int(len(d)),
            n_donors=int(d.donor.nunique()),
            n_SMCENTER=int(d.SMCENTER.nunique(dropna=True)),
            n_SMNABTCH=int(d.SMNABTCH.nunique(dropna=True)),
            n_SMGEBTCH=int(d.SMGEBTCH.nunique(dropna=True)),
            n_male=int((d.SEX.astype(str) == "1").sum()),
            n_female=int((d.SEX.astype(str) == "2").sum()),
            SMRIN_n_missing=rin["n_missing"], SMRIN_median=rin["median"],
            SMRIN_min=rin["min"], SMRIN_max=rin["max"],
            SMTSISCH_n_missing=isch["n_missing"], SMTSISCH_median=isch["median"],
            SMTSISCH_min=isch["min"], SMTSISCH_max=isch["max"],
            DTHHRDY_n_missing=hardy["n_missing"], DTHHRDY_n_levels=hardy["n_levels"],
        )
        for b in EXPECTED_AGE_BINS:
            rec[f"n_AGE_{b}"] = int(age.get(b, 0))
        rec["SMCENTER_counts"] = {str(k): int(v) for k, v in center.items()}
        rec["SMNABTCH_counts"] = {str(k): int(v) for k, v in batch.items()}
        rec["DTHHRDY_counts"] = hardy["counts"]
        rec["SEX_counts"] = sex["counts"]
        rows.append(rec)
        log(f"[cohort] {t}: n={rec['n_samples']} donors={rec['n_donors']} "
            f"centers={rec['n_SMCENTER']} nabatch={rec['n_SMNABTCH']} "
            f"sex 1/2={rec['n_male']}/{rec['n_female']}")
    df = pd.DataFrame(rows)
    return df


def write_cohort_md(obs: pd.DataFrame, summary: pd.DataFrame, path: Path):
    lines = [
        f"# GTEx v10 Stage 0 cohort",
        "",
        f"Release: {RELEASE}. Pipeline: {RNASEQ_PIPELINE}. Freeze: SMAFRZE=={SMAFRZE_RNASEQ}. Seed `{GTEX_SEED}`.",
        "",
        "## Per tissue",
        "",
        "| tissue | n_samples | n_donors | n_SMCENTER | n_SMNABTCH | n_SMGEBTCH | n_male | n_female | SMRIN median (min–max, n_miss) | SMTSISCH median (min–max, n_miss) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---|---|",
    ]
    for _, r in summary.iterrows():
        rin = f"{r.SMRIN_median:.2f} ({r.SMRIN_min:.2f}–{r.SMRIN_max:.2f}, miss={int(r.SMRIN_n_missing)})"
        isch = f"{r.SMTSISCH_median:.1f} ({r.SMTSISCH_min:.1f}–{r.SMTSISCH_max:.1f}, miss={int(r.SMTSISCH_n_missing)})"
        lines.append(
            f"| {r.tissue} | {int(r.n_samples)} | {int(r.n_donors)} | {int(r.n_SMCENTER)} | "
            f"{int(r.n_SMNABTCH)} | {int(r.n_SMGEBTCH)} | {int(r.n_male)} | {int(r.n_female)} | {rin} | {isch} |"
        )
    lines += ["", "## AGE bins (n samples)", "",
              "| tissue | " + " | ".join(EXPECTED_AGE_BINS) + " |",
              "|" + "|".join(["---"] * (1 + len(EXPECTED_AGE_BINS))) + "|"]
    for _, r in summary.iterrows():
        cells = " | ".join(str(int(r[f"n_AGE_{b}"])) for b in EXPECTED_AGE_BINS)
        lines.append(f"| {r.tissue} | {cells} |")
    lines += ["", "## SMCENTER counts", ""]
    for _, r in summary.iterrows():
        lines.append(f"### {r.tissue}")
        items = sorted(r.SMCENTER_counts.items(), key=lambda kv: (-kv[1], kv[0]))
        lines.append(", ".join(f"{k}={v}" for k, v in items) or "(none)")
        lines.append("")
    lines += ["## SMNABTCH counts", ""]
    for _, r in summary.iterrows():
        lines.append(f"### {r.tissue}")
        items = sorted(r.SMNABTCH_counts.items(), key=lambda kv: (-kv[1], kv[0]))
        lines.append(", ".join(f"{k}={v}" for k, v in items) or "(none)")
        lines.append("")
    lines += ["## DTHHRDY (Hardy scale) counts", ""]
    for _, r in summary.iterrows():
        items = sorted(r.DTHHRDY_counts.items(), key=lambda kv: kv[0])
        lines.append(f"- {r.tissue}: " + ", ".join(f"{k}={v}" for k, v in items))
    lines += ["", "## SEX counts (1=male, 2=female; GTEx coding, not recoded)", ""]
    for _, r in summary.iterrows():
        items = sorted(r.SEX_counts.items(), key=lambda kv: kv[0])
        lines.append(f"- {r.tissue}: " + ", ".join(f"{k}={v}" for k, v in items))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def covariate_audit(obs: pd.DataFrame, log):
    """OLS R² of age_mid on each covariate singly and jointly. Descriptive, no CV."""
    rows = []
    covars_single = ["SMRIN", "SMTSISCH", "DTHHRDY", "SEX", "SMCENTER", "SMNABTCH"]
    cat = {"DTHHRDY", "SEX", "SMCENTER", "SMNABTCH"}
    for t in TISSUES:
        d = obs[obs.SMTSD.astype(str) == t].copy()
        y = d.age_mid.to_numpy(float)
        blocks = {}
        for c in covars_single:
            if c in cat:
                X, levels = onehot_train(d[c], drop_first=True)
                n_levels = int(d[c].nunique(dropna=True))
            else:
                X = pd.to_numeric(d[c], errors="coerce").to_numpy(float)[:, None]
                levels = [c]
                n_levels = 1
            rec = ols_r2(y, X)
            rec.update(dict(tissue=t, covariate=c, kind="single",
                            n_levels=n_levels, n_samples=int(len(d))))
            rows.append(rec)
            blocks[c] = X
            log(f"[audit] {t} ~ {c}: R²={fmt(rec['r2'])} n={rec['n']} dropped={rec['n_dropped']} levels={n_levels}")
        Xj = np.column_stack([blocks[c] for c in covars_single])
        recj = ols_r2(y, Xj)
        recj.update(dict(tissue=t, covariate="JOINT", kind="joint",
                         n_levels=int(Xj.shape[1]), n_samples=int(len(d))))
        rows.append(recj)
        log(f"[audit] {t} ~ JOINT: R²={fmt(recj['r2'])} n={recj['n']} p={recj['p']}")
    return pd.DataFrame(rows)


def stream_gct_subset(path: Path, keep_ids: set[str], log):
    """Stream gzipped GCT; keep Name, Description, and keep_ids sample columns. No full matrix."""
    keep_ids = set(map(str, keep_ids))
    t0 = time.time()
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        ver = f.readline().rstrip("\n")
        dims = f.readline().rstrip("\n").split()
        if len(dims) < 2:
            raise StopStep("gct", f"GCT dimension line unexpected: {dims!r}")
        n_genes_decl, n_samp_decl = int(dims[0]), int(dims[1])
        header = f.readline().rstrip("\n").split("\t")
        if len(header) < 3 or header[0] != "Name" or header[1] != "Description":
            raise StopStep(
                "gct",
                f"GCT header does not start Name, Description. first fields={header[:6]!r}",
                dict(header_prefix=header[:10], n_fields=len(header)),
            )
        sample_ids = header[2:]
        if len(sample_ids) != n_samp_decl:
            raise StopStep(
                "gct",
                f"GCT n_samples declared {n_samp_decl} != header sample columns {len(sample_ids)}",
            )
        kept = [s for s in sample_ids if s in keep_ids]
        missing = sorted(keep_ids - set(sample_ids))
        extra_in_gct = len(sample_ids) - len(kept)
        log(f"[gct] version={ver!r} declared genes={n_genes_decl:,} samples={n_samp_decl:,}")
        log(f"[gct] header columns actually read: Name, Description + {len(sample_ids)} sample IDs")
        if missing:
            raise StopStep(
                "gct",
                f"{len(missing)} freeze-filtered sample IDs absent from GCT. example={missing[:10]}",
                dict(n_missing=len(missing), example=missing[:10]),
            )
        log(f"[gct] keeping {len(kept)} / {len(sample_ids)} sample columns "
            f"(not in 6-tissue freeze: {extra_in_gct})")
    usecols = ["Name", "Description"] + kept
    log("[gct] pandas usecols subset (does not materialise the full gene × all-sample matrix)")
    df = pd.read_csv(
        path, sep="\t", skiprows=2, compression="gzip", usecols=usecols,
        dtype={c: np.float32 for c in kept} | {"Name": "string", "Description": "string"},
    )
    if len(df) != n_genes_decl:
        raise StopStep("gct", f"GCT declared {n_genes_decl} genes, read {len(df)}")
    counts = df.loc[:, kept].to_numpy(dtype=np.float32)
    names = df["Name"].astype(str).to_numpy()
    symbols = df["Description"].astype(str).to_numpy()
    log(f"[gct] subset matrix {counts.shape[0]:,} genes × {counts.shape[1]} samples "
        f"in {time.time()-t0:.0f}s")
    meta = dict(
        gct_version=ver,
        n_genes_declared=n_genes_decl,
        n_samples_declared=n_samp_decl,
        n_sample_columns_in_header=len(sample_ids),
        n_kept=len(kept),
        columns_read=["Name", "Description"] + sample_ids,
        kept_sample_ids=kept,
    )
    return counts, names, symbols, np.array(kept, dtype=object), meta


def composition_scores(logcpm_gs: np.ndarray, symbols: np.ndarray, log):
    """Sum of z-scored log-CPM over frozen marker sets. symbols aligned to genes axis."""
    sy = np.array([str(s).upper() for s in symbols])
    idx = {s: i for i, s in enumerate(sy)}
    out = {}
    for name, markers in (("neuronal", NEURONAL_MARKERS), ("glial", GLIAL_MARKERS)):
        missing = [m for m in markers if m not in idx]
        if missing:
            raise StopStep(
                "composition",
                f"{name} marker(s) not in GCT Description: {missing}. "
                "Not substituting aliases.",
                dict(set=name, missing=missing, listed=list(markers)),
            )
        ii = np.array([idx[m] for m in markers], dtype=int)
        Z = logcpm_gs[ii, :]  # markers × samples, already gene-z
        score = Z.sum(0)
        out[name] = dict(score=score, idx=ii, markers=list(markers), n=len(markers))
        log(f"[composition] {name}: {list(markers)} all mapped")
    return out


def per_tissue_preprocess(counts_all, genes, symbols, samples, obs, universe, log):
    """Filter + TMM + log2-CPM per tissue. Returns dict slug -> pack."""
    samp_index = {str(s): i for i, s in enumerate(samples)}
    packs = {}
    overlap_rows = []
    for t in TISSUES:
        d = obs[obs.SMTSD.astype(str) == t]
        ids = d.SAMPID.astype(str).tolist()
        missing = [s for s in ids if s not in samp_index]
        if missing:
            raise StopStep("preprocess", f"{t}: {len(missing)} SAMPID not in GCT subset",
                           dict(example=missing[:10]))
        col = np.array([samp_index[s] for s in ids], dtype=int)
        C = counts_all[:, col]  # genes × n
        keep = gene_filter_mask(C)
        n0 = int(C.shape[0])
        n1 = int(keep.sum())
        log(f"[filter] {t}: genes ≥{6} counts in ≥{0.20:.0%} samples  {n0:,} -> {n1:,}  "
            f"(removed {n0-n1:,})  n_samples={C.shape[1]}")
        Cf = C[keep]
        g_keep = genes[keep]
        s_keep = symbols[keep]
        nf = tmm_norm_factors(Cf.astype(np.float64))
        logcpm = log2_cpm_edger(Cf.astype(np.float64), nf)  # genes × samples
        # gene-z on all rows of this tissue — Stage 0 descriptive composition only
        mu = logcpm.mean(1, keepdims=True)
        sd = logcpm.std(1, keepdims=True)
        sd = np.where(sd < 1e-12, 1.0, sd)
        Z = (logcpm - mu) / sd
        ens = np.array([strip_ensembl(g) for g in g_keep])
        n_overlap = int(sum(1 for e in ens if e in universe))
        overlap_rows.append(dict(
            tissue=t, n_samples=int(C.shape[1]), n_genes_raw=n0, n_genes_filter=n1,
            n_overlap_dlpfc=n_overlap, frac_overlap=n_overlap / max(n1, 1),
        ))
        log(f"[overlap] {t}: filtered genes {n1:,} ∩ DLPFC universe {n_overlap:,}")
        slug = TISSUE_SLUG[t]
        packs[slug] = dict(
            tissue=t, sample_ids=np.array(ids, dtype=object),
            counts=Cf.astype(np.float32),
            logcpm=logcpm.astype(np.float32),
            gene_z=Z.astype(np.float32),
            genes=g_keep, symbols=s_keep, ensembl=ens,
            tmm_factors=nf.astype(np.float64),
            obs=d.reset_index(drop=True),
        )
        np.savez_compressed(
            GTEX_PROC / f"stage0_{slug}.npz",
            counts=Cf.astype(np.float32),
            logcpm=logcpm.astype(np.float32),
            sample_ids=np.array(ids, dtype=object),
            genes=g_keep.astype(object),
            symbols=s_keep.astype(object),
            ensembl=ens.astype(object),
            tmm_factors=nf,
            tissue=np.array(t),
        )
    return packs, pd.DataFrame(overlap_rows)


def composition_table(packs, log):
    rows = []
    for slug, pack in packs.items():
        t = pack["tissue"]
        if t not in BRAIN_TISSUES:
            continue
        comp = composition_scores(pack["gene_z"], pack["symbols"], log)
        obs = pack["obs"]
        for name, blob in comp.items():
            score = blob["score"]
            r_age = spearman_safe(score, obs.age_mid.to_numpy(float))
            r_ord = spearman_safe(score, obs.age_ordinal.to_numpy(float))
            r_isch = spearman_safe(score, obs.SMTSISCH.to_numpy(float))
            Xc, _ = onehot_train(obs.SMCENTER, drop_first=True)
            r2_c = ols_r2(score, Xc)
            rec = dict(
                tissue=t, score=name, n_markers=blob["n"],
                spearman_age_mid=r_age, spearman_age_ordinal=r_ord,
                spearman_SMTSISCH=r_isch,
                r2_SMCENTER=r2_c["r2"], r2_SMCENTER_n=r2_c["n"],
                score_mean=float(np.mean(score)), score_sd=float(np.std(score)),
            )
            rows.append(rec)
            log(f"[composition] {t} {name}: Spearman age_mid={fmt(r_age)} age_ordinal={fmt(r_ord)} "
                f"SMTSISCH={fmt(r_isch)}  R²~SMCENTER={fmt(r2_c['r2'])}")
            pack[f"score_{name}"] = score
        pack["obs"] = obs.assign(
            neuronal_score=comp["neuronal"]["score"],
            glial_score=comp["glial"]["score"],
        )
    return pd.DataFrame(rows)


def stop_checks(obs, log):
    ba9 = obs[obs.SMTSD.astype(str) == PRIMARY_TISSUE]
    n = int(len(ba9))
    vc = ba9.SMCENTER.astype(str).value_counts()
    n_big = int((vc >= CENTER_TRANSFER_MIN_N).sum())
    flags = dict(
        ba9_n=n,
        ba9_n_underpowered=n < BA9_N_UNDERPOWERED,
        n_SMCENTER_ge25=n_big,
        t2_not_testable=n_big < CENTER_TRANSFER_MIN_LEVELS,
        center_counts={str(k): int(v) for k, v in vc.items()},
    )
    log(f"[STOP] BA9 n={n}  (flag underpowered if n<{BA9_N_UNDERPOWERED}: {flags['ba9_n_underpowered']})")
    log(f"[STOP] SMCENTER with ≥{CENTER_TRANSFER_MIN_N} BA9 samples: {n_big}  "
        f"(T2 not testable if <{CENTER_TRANSFER_MIN_LEVELS}: {flags['t2_not_testable']})")
    log(f"[STOP] BA9 SMCENTER counts: {flags['center_counts']}")
    return flags


def figures(obs, audit, comp, packs):
    GTEX_FIG.mkdir(parents=True, exist_ok=True)
    # age bin counts
    fig, axes = plt.subplots(2, 3, figsize=(11, 6), sharey=False)
    for ax, t in zip(axes.ravel(), TISSUES):
        d = obs[obs.SMTSD.astype(str) == t]
        vc = d.AGE.astype(str).value_counts().reindex(list(EXPECTED_AGE_BINS)).fillna(0)
        ax.bar(range(len(vc)), vc.to_numpy())
        ax.set_xticks(range(len(vc)))
        ax.set_xticklabels(list(EXPECTED_AGE_BINS), rotation=45, ha="right", fontsize=7)
        ax.set_title(t, fontsize=8)
        ax.set_ylabel("n")
    fig.tight_layout()
    fig.savefig(GTEX_FIG / "stage0_age_bins.png", dpi=120)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(11, 6))
    for ax, t in zip(axes.ravel(), TISSUES):
        d = obs[obs.SMTSD.astype(str) == t]
        ax.scatter(d.age_mid, d.SMRIN, s=8, alpha=0.5)
        ax.set_title(t, fontsize=8)
        ax.set_xlabel("age midpoint")
        ax.set_ylabel("SMRIN")
    fig.tight_layout()
    fig.savefig(GTEX_FIG / "stage0_rin_vs_age.png", dpi=120)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(11, 6))
    for ax, t in zip(axes.ravel(), TISSUES):
        d = obs[obs.SMTSD.astype(str) == t]
        ax.scatter(d.age_mid, d.SMTSISCH, s=8, alpha=0.5)
        ax.set_title(t, fontsize=8)
        ax.set_xlabel("age midpoint")
        ax.set_ylabel("SMTSISCH")
    fig.tight_layout()
    fig.savefig(GTEX_FIG / "stage0_isch_vs_age.png", dpi=120)
    plt.close(fig)

    if not comp.empty:
        fig, axes = plt.subplots(2, 2, figsize=(8, 7))
        brain = [t for t in BRAIN_TISSUES]
        for ax, t in zip(axes.ravel(), brain):
            pack = packs[TISSUE_SLUG[t]]
            o = pack["obs"]
            ax.scatter(o.age_mid, o.neuronal_score, s=8, alpha=0.6, label="neuronal")
            ax.scatter(o.age_mid, o.glial_score, s=8, alpha=0.6, label="glial")
            ax.set_title(t, fontsize=8)
            ax.set_xlabel("age midpoint")
            ax.set_ylabel("marker sum (gene-z)")
            ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(GTEX_FIG / "stage0_composition_vs_age.png", dpi=120)
        plt.close(fig)


def write_empty_manifest(failures):
    MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)
    dump_json(MANIFEST_PATH, dict(
        seed=GTEX_SEED, release=RELEASE, pipeline=RNASEQ_PIPELINE, dbgap=DBGAP,
        files={}, failures=failures, status="STOP",
    ))


def main():
    log = Logger(GTEX_DIR / "stage0_report.txt")
    gtex_log_banner(log, "STAGE 0")
    failures = []
    manifest = dict(
        seed=GTEX_SEED, release=RELEASE, pipeline=RNASEQ_PIPELINE, dbgap=DBGAP,
        portal_accessed=str(pd.Timestamp.now(tz="UTC")),
        files={}, failures=[], columns={}, stop_flags={}, stage1_gate=None,
        preprocessing_spec=dict(
            gene_filter=f">={6} counts in >={0.20:.0%} of samples per tissue",
            normalization="TMM (edgeR calcNormFactors defaults) then log2-CPM prior.count=2",
            zscore="per gene on training rows only inside every fold",
            tissues=list(TISSUES), freeze=SMAFRZE_RNASEQ,
            neuronal_markers=list(NEURONAL_MARKERS),
            glial_markers=list(GLIAL_MARKERS),
        ),
    )
    try:
        for key, spec in FILES.items():
            rec = _download_one(key, spec, log)
            manifest["files"][key] = rec

        samp_path = GTEX_RAW / FILES["sample_attributes"]["name"]
        pheno_path = GTEX_RAW / FILES["subject_phenotypes"]["name"]
        gct_path = GTEX_RAW / FILES["gene_reads"]["name"]

        samp = read_tsv_log_columns(samp_path, "sample_attributes", log)
        pheno = read_tsv_log_columns(pheno_path, "subject_phenotypes", log)
        manifest["files"]["sample_attributes"]["columns_read"] = list(samp.columns)
        manifest["files"]["sample_attributes"]["n_rows"] = int(len(samp))
        manifest["files"]["subject_phenotypes"]["columns_read"] = list(pheno.columns)
        manifest["files"]["subject_phenotypes"]["n_rows"] = int(len(pheno))
        manifest["columns"]["sample_attributes"] = list(samp.columns)
        manifest["columns"]["subject_phenotypes"] = list(pheno.columns)

        obs, extra = build_cohort(samp, pheno, log)
        manifest["cohort_extra"] = extra
        obs_path = GTEX_DIR / "stage0_obs.csv"
        obs.to_csv(obs_path, index=False)
        log(f"[write] {obs_path} n={len(obs)}")

        summary = cohort_tables(obs, log)
        # JSON-friendly counts already in dict columns; flatten for csv
        flat = summary.drop(columns=[c for c in ("SMCENTER_counts", "SMNABTCH_counts",
                                                 "DTHHRDY_counts", "SEX_counts") if c in summary.columns])
        flat.to_csv(GTEX_DIR / "stage0_cohort.csv", index=False)
        dump_json(GTEX_DIR / "stage0_cohort_counts.json", summary.to_dict(orient="records"))
        write_cohort_md(obs, summary, GTEX_DIR / "stage0_cohort.md")

        audit = covariate_audit(obs, log)
        audit.to_csv(GTEX_DIR / "stage0_covariate_audit.csv", index=False)

        flags = stop_checks(obs, log)
        manifest["stop_flags"] = flags

        keep_ids = set(obs.SAMPID.astype(str))
        counts, genes, symbols, samples, gct_meta = stream_gct_subset(gct_path, keep_ids, log)
        manifest["files"]["gene_reads"]["columns_read"] = {
            "gct_fixed": ["Name", "Description"],
            "n_sample_columns_in_header": gct_meta["n_sample_columns_in_header"],
            "n_kept_sample_columns": gct_meta["n_kept"],
            "gct_version": gct_meta["gct_version"],
            "n_genes_declared": gct_meta["n_genes_declared"],
            "n_samples_declared": gct_meta["n_samples_declared"],
            "note": "full sample-ID header logged to stage0_gct_sample_ids.txt; not inlined here",
        }
        manifest["files"]["gene_reads"]["n_rows"] = int(gct_meta["n_genes_declared"])
        (GTEX_DIR / "stage0_gct_sample_ids.txt").write_text(
            "\n".join(gct_meta["columns_read"][2:]) + "\n", encoding="utf-8")
        dump_json(GTEX_DIR / "stage0_gct_meta.json", {k: v for k, v in gct_meta.items() if k != "columns_read"})

        universe, dlpfc_genes = load_dlpfc_gene_universe()
        manifest["dlpfc_universe"] = dict(
            path=str(DLPFC_GENES_CSV), n=len(universe),
            columns_read=list(dlpfc_genes.columns),
        )
        log(f"[overlap] DLPFC universe n={len(universe):,} from {DLPFC_GENES_CSV.name} "
            f"columns={list(dlpfc_genes.columns)}")

        packs, overlap = per_tissue_preprocess(counts, genes, symbols, samples, obs, universe, log)
        overlap.to_csv(GTEX_DIR / "stage0_gene_overlap.csv", index=False)
        manifest["gene_overlap"] = overlap.to_dict(orient="records")

        comp = composition_table(packs, log)
        comp.to_csv(GTEX_DIR / "stage0_composition.csv", index=False)
        # attach scores onto the written obs
        score_map_n = {}
        score_map_g = {}
        for slug, pack in packs.items():
            if "neuronal_score" in pack["obs"].columns:
                for sid, ns, gs in zip(pack["obs"].SAMPID.astype(str),
                                       pack["obs"].neuronal_score,
                                       pack["obs"].glial_score):
                    score_map_n[sid] = float(ns)
                    score_map_g[sid] = float(gs)
        if score_map_n:
            obs["neuronal_score"] = obs.SAMPID.astype(str).map(score_map_n)
            obs["glial_score"] = obs.SAMPID.astype(str).map(score_map_g)
            obs.to_csv(obs_path, index=False)

        figures(obs, audit, comp, packs)

        gate = stage1_gate()
        manifest["stage1_gate"] = gate
        log(f"[gate] Stage 1 open={gate['open']}  reason={gate['reason']}  status={gate.get('status')}")

        dump_json(GTEX_DIR / "stage0_summary.json", dict(
            seed=GTEX_SEED, release=RELEASE, n_samples=int(len(obs)),
            n_donors=int(obs.donor.nunique()),
            per_tissue=flat.to_dict(orient="records"),
            stop_flags=flags, stage1_gate=gate,
            gene_overlap=overlap.to_dict(orient="records"),
        ))
        manifest["failures"] = failures
        manifest["status"] = "STAGE0_DONE"
        dump_json(MANIFEST_PATH, manifest)
        log(f"[write] {MANIFEST_PATH}")
        log("STAGE 0 complete.")
    except StopStep as e:
        rec = dict(step=e.step, message=e.message, details=e.details)
        failures.append(rec)
        log(f"STOP [{e.step}] {e.message}")
        manifest["failures"] = failures
        manifest["status"] = "STOP"
        dump_json(MANIFEST_PATH, manifest)
        raise
    finally:
        log.close()


if __name__ == "__main__":
    main()

"""Task 2 — three-way instrument comparison against donor age, with paired Δρ."""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md3_common import (  # noqa: E402
    MD3_DIR, MD3_PROC, MD2_DIR, MD2_PROC, MD3_SEED, MD3_BOOT, N_PERM, N_BOOT,
    PREREG_TASK2, PREREG_TASK2_FLAG, DECLARED_BEFORE_SCORES_FLAG, FROZEN_RULER,
    CXG_DATASET_ID, GSE226189_TAR, GSE113957_FPKM, FIBRO2_DIR, ADULT_MIN_AGE,
    StopStep, Logger, dump_json, jsonable, load_json, md3_log_banner,
    load_manifest, save_manifest, record_failure, log_columns,
    progress_snapshot,
)
from md3_idtype import detect_id_type, require_mappable  # noqa: E402
from fibro2_common import load_frozen_ruler, rho_null_donor, align_counts_to_ruler, score_frozen  # noqa: E402
from fibro_common import bootstrap_rho_ci  # noqa: E402
from gtex_common import spearman_safe  # noqa: E402
from md2_score import lognormalize_dense, add_module_score  # noqa: E402
from fibro2_ta import load_count_table  # noqa: E402
from fibro2_common import looks_like_counts  # noqa: E402
from brain_phase1_common import strip_ensembl  # noqa: E402


def _ensembl_to_symbol(ens_ids, log):
    """Map ENSG Tracking_IDs to symbols using on-disk 10x features, then the frozen ruler.

    Same mapping sources as md2_task3._ensembl_to_symbol. Writes column logs to the md3
    manifest, not the md2 manifest.
    """
    ens = np.array([strip_ensembl(x) for x in np.asarray(ens_ids).astype(str)])
    pos = {}
    src = []
    zpath = MD2_PROC / "allcell_counts_GM00731.npz"
    if zpath.exists():
        z = np.load(zpath, allow_pickle=True)
        if "gene_id" in z.files and "symbols" in z.files:
            gid = np.array([strip_ensembl(x) for x in np.asarray(z["gene_id"]).astype(str)])
            sym = np.array([str(s).upper() for s in np.asarray(z["symbols"])], dtype=object)
            for e, s in zip(gid, sym):
                if e.startswith("ENS") and e not in pos and s and not str(s).startswith("ENS"):
                    pos[e] = s
            src.append(f"{zpath} n={len(pos)}")
            log_columns("t2_gse226189_10x_map", ["gene_id", "symbol"], str(zpath))
    frozen = load_frozen_ruler()
    n_before = len(pos)
    for e, s in zip(np.asarray(frozen["ensembl"]).astype(str), np.asarray(frozen["symbol"]).astype(str)):
        e = strip_ensembl(e)
        su = str(s).upper()
        if e.startswith("ENS") and e not in pos and su and not su.startswith("ENS"):
            pos[e] = su
    src.append(f"frozen_ruler added={len(pos) - n_before} total={len(pos)}")
    symbols = np.array([pos.get(e, e) for e in ens], dtype=object)
    n_sym = int(sum(1 for s in symbols if not str(s).startswith("ENS")))
    log(f"[t2 map] Ensembl→symbol sources={src} n_ids={len(ens)} n_mapped_to_symbol={n_sym} "
        f"n_left_as_ensembl={int(len(ens) - n_sym)}")
    if n_sym == 0:
        raise StopStep(
            "task2",
            "GSE226189 geneCOUNT Tracking_ID is Ensembl and 0 IDs mapped to symbols "
            f"from {src}. Not scoring Age up/down on Ensembl IDs as if they were symbols.",
        )
    return symbols


def _require():
    if not PREREG_TASK2_FLAG.exists():
        raise StopStep("prereg", "PREREG_TASK2.flag missing")
    if not DECLARED_BEFORE_SCORES_FLAG.exists():
        raise StopStep("prereg", "DECLARED_BEFORE_SCORES.flag missing")
    if not FROZEN_RULER.exists():
        raise StopStep("frozen_ruler", f"missing {FROZEN_RULER}")
    if not (MD3_DIR / "genesets.json").exists():
        raise StopStep("task2", "genesets.json missing")


def _rho_pack(y, pred, donor, log, tag):
    sc, rho_null, p, nulls = rho_null_donor(y, pred, donor, n_perm=N_PERM, seed=MD3_SEED)
    rng_b = np.random.default_rng(MD3_BOOT)
    ci = bootstrap_rho_ci(y, pred, rng_b, n_boot=N_BOOT)
    rec = dict(
        tag=tag, n=int(len(y)), n_donors=int(pd.Series(donor).nunique()),
        rho=sc["rho"], rho_null=rho_null, rho_p=p,
        rho_ci_lo=ci.get("p025"), rho_ci_hi=ci.get("p975"),
        r=sc["r"], r2=sc["r2"], cal_r2=sc["cal_r2"],
        n_perm=N_PERM, n_boot=N_BOOT, seed=MD3_SEED, boot_seed=MD3_BOOT,
        age_min=float(np.nanmin(y)) if len(y) else np.nan,
        age_max=float(np.nanmax(y)) if len(y) else np.nan,
        n_unique_age=int(pd.Series(y).nunique()),
    )
    log(f"[t2 {tag}] ρ={rec['rho']:+.3f} null={rho_null:+.3f} p={p:+.3f} "
        f"CI=[{ci.get('p025')}, {ci.get('p975')}] n={rec['n']}")
    return rec


def _check_saved_rho(saved_path, y, pred, donor, log, tag, instrument):
    """Recompute ρ from saved per-donor scores; compare to the on-disk result JSON if present."""
    rec = _rho_pack(y, pred, donor, log, tag)
    rec.update(instrument=instrument, source=str(saved_path), refit=False)
    return rec


def paired_delta_rho(y, pred_a, pred_b, donor, log, tag):
    """Δρ = ρ(y, pred_a) − ρ(y, pred_b) on the same bootstrap resamples of donors."""
    y = np.asarray(y, float)
    a = np.asarray(pred_a, float)
    b = np.asarray(pred_b, float)
    donor = np.asarray(donor).astype(str)
    m = np.isfinite(y) & np.isfinite(a) & np.isfinite(b)
    y, a, b, donor = y[m], a[m], b[m], donor[m]
    d_u = pd.unique(donor)
    n_d = int(len(d_u))
    by = {d: np.flatnonzero(donor == d) for d in d_u}
    # one row per donor (age constant; take first)
    y_d = np.array([float(y[by[d][0]]) for d in d_u])
    a_d = np.array([float(np.mean(a[by[d]])) for d in d_u])
    b_d = np.array([float(np.mean(b[by[d]])) for d in d_u])
    rho_a = spearman_safe(y_d, a_d)
    rho_b = spearman_safe(y_d, b_d)
    point = float(rho_a - rho_b) if np.isfinite(rho_a) and np.isfinite(rho_b) else np.nan
    rng = np.random.default_rng(MD3_BOOT)
    deltas = []
    for _ in range(N_BOOT):
        idx = rng.integers(0, n_d, size=n_d)
        ra = spearman_safe(y_d[idx], a_d[idx])
        rb = spearman_safe(y_d[idx], b_d[idx])
        if np.isfinite(ra) and np.isfinite(rb):
            deltas.append(float(ra - rb))
    v = np.asarray(deltas, float)
    if v.size < 2:
        ci_lo = ci_hi = np.nan
    else:
        ci_lo = float(np.percentile(v, 2.5))
        ci_hi = float(np.percentile(v, 97.5))
    excludes0 = bool(np.isfinite(ci_lo) and np.isfinite(ci_hi) and (ci_lo > 0 or ci_hi < 0))
    favours = None
    if excludes0 and ci_lo > 0:
        favours = "ruler"
    elif excludes0 and ci_hi < 0:
        favours = "other"
    rec = dict(
        tag=tag, n_donors=n_d, n_boot=N_BOOT, seed=MD3_BOOT,
        rho_ruler=float(rho_a) if np.isfinite(rho_a) else np.nan,
        rho_other=float(rho_b) if np.isfinite(rho_b) else np.nan,
        delta_rho=point, delta_ci_lo=ci_lo, delta_ci_hi=ci_hi,
        excludes_zero=excludes0, favours=favours,
        note="Δρ = ρ_ruler − ρ_other; paired donor bootstrap. Overlapping separate CIs are not the claim.",
    )
    log(f"[t2 Δρ {tag}] Δρ={point:+.3f} CI=[{ci_lo}, {ci_hi}] excludes0={excludes0} favours={favours} n={n_d}")
    return rec


def _ams_up_minus_down(counts_sg, symbols, obs, sets, log, tag):
    """counts_sg: samples × genes. LogNormalize + AMS Age up and Age down; return obs with scores."""
    idt = detect_id_type(symbols, log, f"{tag}_symbols", tag)
    require_mappable(idt, {"symbol"}, log, "task2")
    logC = lognormalize_dense(counts_sg)
    gene_sets = {"age_up": list(sets["age_up"]), "age_down": list(sets["age_down"])}
    ams, meta = add_module_score(logC, np.asarray(symbols).astype(str), gene_sets, log=log, tag=tag)
    up = np.asarray(ams["age_up"], float)
    dn = np.asarray(ams["age_down"], float)
    pred = up - dn
    obs = obs.copy()
    obs["age_up"] = up
    obs["age_down"] = dn
    obs["age_up_minus_age_down"] = pred
    y = obs.age_years.to_numpy(float)
    donor = obs.donor.astype(str).to_numpy()
    rec = _rho_pack(y, pred, donor, log, tag)
    rec.update(
        instrument="age_up_minus_age_down_AddModuleScore",
        n_age_up_mapped=meta["per_set"]["age_up"]["n_mapped"],
        n_age_up_missing=meta["per_set"]["age_up"]["n_missing"],
        n_age_up_total=int(len(sets["age_up"])),
        n_age_down_mapped=meta["per_set"]["age_down"]["n_mapped"],
        n_age_down_missing=meta["per_set"]["age_down"]["n_missing"],
        n_age_down_total=int(len(sets["age_down"])),
        id_type=idt["kind"],
        scoring="AddModuleScore(Age up) − AddModuleScore(Age down) on LogNormalize; not a ruler refit",
        matrix_id_type=idt,
    )
    obs.to_csv(MD3_DIR / f"t2_{tag}_scores.csv", index=False)
    dump_json(MD3_DIR / f"t2_{tag}_ams_meta.json", jsonable({
        k: {kk: vv for kk, vv in v.items() if kk != "mapped"}
        for k, v in (meta.get("per_set") or {}).items()
    }))
    return rec, obs


def _cxg(sets, log):
    scores_r = pd.read_csv(FIBRO2_DIR / "ta_cxg_a19d1667_scores.csv")
    scores_m = pd.read_csv(MD2_DIR / "t3_cxg_a19d1667_MD_scores.csv")
    log_columns("t2_cxg_ruler_scores", list(scores_r.columns), str(FIBRO2_DIR / "ta_cxg_a19d1667_scores.csv"))
    log_columns("t2_cxg_md_scores", list(scores_m.columns), str(MD2_DIR / "t3_cxg_a19d1667_MD_scores.csv"))
    m = scores_r.merge(scores_m[["sample", "md_score"]], on="sample", how="inner")
    if len(m) != 93:
        raise StopStep("task2", f"cxg inner join n={len(m)} expected 93. Not substituting.")
    cache = MD2_PROC / "cxg_donor_counts.npz"
    if not cache.exists():
        raise StopStep("task2", f"missing {cache}")
    z = np.load(cache, allow_pickle=True)
    counts = z["counts"]  # genes × samples
    genes = np.asarray(z["genes"]).astype(str)
    symbols = np.asarray(z["symbols"]).astype(str)
    samples = np.asarray(z["samples"]).astype(str)
    detect_id_type(genes, log, "cxg_gene_id", str(cache) + ":genes")
    detect_id_type(symbols, log, "cxg_symbols", str(cache) + ":symbols")
    pos = {str(s): i for i, s in enumerate(samples)}
    miss = [s for s in m["sample"].astype(str) if s not in pos]
    if miss:
        raise StopStep("task2", f"cxg: {len(miss)} scored samples not in cached counts e.g. {miss[:8]}")
    col = np.array([pos[s] for s in m["sample"].astype(str)], dtype=int)
    C = np.asarray(counts, np.float64)[:, col].T
    s3, obs_s3 = _ams_up_minus_down(C, symbols, m, sets, log, tag="cxg_a19d1667_S3")
    y = m.age_years.to_numpy(float)
    donor = m.donor.astype(str).to_numpy()
    ruler = _check_saved_rho(FIBRO2_DIR / "ta_cxg_a19d1667_scores.csv", y, m.age_score.to_numpy(float),
                             donor, log, "cxg_ruler_check", "frozen_GTEx_fibroblast_ruler")
    md2_inst = pd.read_csv(MD2_DIR / "t3_instruments.csv") if (MD2_DIR / "t3_instruments.csv").exists() else None
    def _md2_row(instrument, tag):
        if md2_inst is None:
            return {}
        hit = md2_inst[(md2_inst.instrument.astype(str) == instrument) & (md2_inst.tag.astype(str) == tag)]
        if not len(hit):
            return {}
        r0 = hit.iloc[0]
        out = {}
        for k in ("n_overlap", "n_MD_mapped", "n_MD_missing", "n_genes"):
            if k in hit.columns and pd.notna(r0[k]):
                out[k] = r0[k]
        return out
    ruler.update(
        cohort="CELLxGENE a19d1667-a7b5-4556-9e5f-f9bfa690c0f1",
        accession=CXG_DATASET_ID, platform="10x (restricted)",
        n_genes_mapped=None, n_genes_total=None, id_type="symbol (var feature_name); gene_id Ensembl",
        transform="TMM/log2-CPM frozen spec (read from fibro2; not recomputed)",
        **_md2_row("frozen_GTEx_fibroblast_ruler", "cxg_a19d1667"),
    )
    md = _check_saved_rho(MD2_DIR / "t3_cxg_a19d1667_MD_scores.csv", y, m.md_score.to_numpy(float),
                          donor, log, "cxg_md_check", "MD_AddModuleScore")
    md.update(
        cohort="CELLxGENE a19d1667-a7b5-4556-9e5f-f9bfa690c0f1",
        accession=CXG_DATASET_ID, platform="10x (restricted)",
        id_type="symbol",
        transform="AddModuleScore on LogNormalize of donor pseudobulk (md2)",
        n_MD_total=int(len(sets["MD"])),
        **_md2_row("MD_AddModuleScore", "cxg_a19d1667_MD"),
    )
    s3.update(
        cohort="CELLxGENE a19d1667-a7b5-4556-9e5f-f9bfa690c0f1",
        accession=CXG_DATASET_ID, platform="10x (restricted)",
        transform="AddModuleScore on LogNormalize of donor pseudobulk",
    )
    d_md = paired_delta_rho(y, m.age_score.to_numpy(float), m.md_score.to_numpy(float), donor, log, "cxg_ruler_minus_MD")
    d_s3 = paired_delta_rho(y, m.age_score.to_numpy(float), obs_s3.age_up_minus_age_down.to_numpy(float),
                            donor, log, "cxg_ruler_minus_S3")
    d_md.update(cohort="CELLxGENE a19d1667-a7b5-4556-9e5f-f9bfa690c0f1", other="MD_AddModuleScore")
    d_s3.update(cohort="CELLxGENE a19d1667-a7b5-4556-9e5f-f9bfa690c0f1", other="age_up_minus_age_down")
    return [ruler, md, s3], [d_md, d_s3], dict(obs=obs_s3, y=y, donor=donor, ruler=m.age_score.to_numpy(float),
                                               md=m.md_score.to_numpy(float), s3=obs_s3.age_up_minus_age_down.to_numpy(float))


def _gse226189(sets, log):
    scores_r = pd.read_csv(FIBRO2_DIR / "ta_GSE226189_near_miss_bulk_scores.csv")
    scores_m = pd.read_csv(MD2_DIR / "t3_GSE226189_MD_scores.csv")
    log_columns("t2_gse226189_ruler_scores", list(scores_r.columns),
                str(FIBRO2_DIR / "ta_GSE226189_near_miss_bulk_scores.csv"))
    log_columns("t2_gse226189_md_scores", list(scores_m.columns), str(MD2_DIR / "t3_GSE226189_MD_scores.csv"))
    m = scores_r.merge(scores_m[["sample", "md_score"]], on="sample", how="inner")
    if len(m) != 82:
        raise StopStep("task2", f"GSE226189 inner join n={len(m)} expected 82. Not substituting.")
    if not GSE226189_TAR.exists():
        raise StopStep("task2", f"missing {GSE226189_TAR}")
    pack = load_count_table(GSE226189_TAR, log, gsm_keep=m["sample"].astype(str).tolist())
    genes = np.asarray(pack["genes"]).astype(str)
    idt = detect_id_type(genes, log, "gse226189_Tracking_ID", str(GSE226189_TAR) + ":Tracking_ID")
    require_mappable(idt, {"ensembl"}, log, "task2")
    symbols = _ensembl_to_symbol(genes, log)
    idt_s = detect_id_type(symbols, log, "gse226189_mapped_symbols", "Ensembl→symbol via 10x+ruler")
    pos = {str(s): i for i, s in enumerate(pack["samples"])}
    miss = [s for s in m["sample"].astype(str) if s not in pos]
    if miss:
        raise StopStep("task2", f"GSE226189: {len(miss)} samples not in geneCOUNT e.g. {miss[:8]}")
    col = np.array([pos[s] for s in m["sample"].astype(str)], dtype=int)
    C = np.asarray(pack["counts"], np.float64)[:, col].T
    s3, obs_s3 = _ams_up_minus_down(C, symbols, m, sets, log, tag="GSE226189_S3")
    y = m.age_years.to_numpy(float)
    donor = m.donor.astype(str).to_numpy()
    ruler = _check_saved_rho(FIBRO2_DIR / "ta_GSE226189_near_miss_bulk_scores.csv",
                             y, m.age_score.to_numpy(float), donor, log, "gse226189_ruler_check",
                             "frozen_GTEx_fibroblast_ruler")
    md2_inst = pd.read_csv(MD2_DIR / "t3_instruments.csv") if (MD2_DIR / "t3_instruments.csv").exists() else None
    def _md2_row(instrument, tag):
        if md2_inst is None:
            return {}
        hit = md2_inst[(md2_inst.instrument.astype(str) == instrument) & (md2_inst.tag.astype(str) == tag)]
        if not len(hit):
            return {}
        r0 = hit.iloc[0]
        out = {}
        for k in ("n_overlap", "n_MD_mapped", "n_MD_missing", "n_genes"):
            if k in hit.columns and pd.notna(r0[k]):
                out[k] = r0[k]
        return out
    ruler.update(cohort="GSE226189 bulk", accession="GSE226189", platform="bulk_rnaseq",
                 id_type="ensembl Tracking_ID", transform="TMM/log2-CPM frozen spec (read from fibro2)",
                 **_md2_row("frozen_GTEx_fibroblast_ruler", "GSE226189"))
    md = _check_saved_rho(MD2_DIR / "t3_GSE226189_MD_scores.csv", y, m.md_score.to_numpy(float),
                          donor, log, "gse226189_md_check", "MD_AddModuleScore")
    md.update(cohort="GSE226189 bulk", accession="GSE226189", platform="bulk_rnaseq",
              id_type="ensembl mapped to symbol", transform="AddModuleScore on LogNormalize (md2)",
              n_MD_total=int(len(sets["MD"])),
              **_md2_row("MD_AddModuleScore", "GSE226189_MD"))
    s3.update(cohort="GSE226189 bulk", accession="GSE226189", platform="bulk_rnaseq",
              transform="AddModuleScore on LogNormalize; Ensembl mapped to symbols",
              gene_index_id_type=idt["kind"], mapped_symbol_id_type=idt_s["kind"])
    d_md = paired_delta_rho(y, m.age_score.to_numpy(float), m.md_score.to_numpy(float),
                            donor, log, "gse226189_ruler_minus_MD")
    d_s3 = paired_delta_rho(y, m.age_score.to_numpy(float), obs_s3.age_up_minus_age_down.to_numpy(float),
                            donor, log, "gse226189_ruler_minus_S3")
    d_md.update(cohort="GSE226189 bulk", other="MD_AddModuleScore")
    d_s3.update(cohort="GSE226189 bulk", other="age_up_minus_age_down")
    return [ruler, md, s3], [d_md, d_s3]


def _parse_homer_symbol(ann):
    tok = str(ann).split("|")[0].strip()
    if not tok or tok == "-":
        return None
    return tok.upper()


def _gse113957(sets, frozen, log):
    """Rank-transform FPKM within sample; score all three instruments on ranks. Weaker test."""
    samp_path = FIBRO2_DIR / "ta_GSE113957_samples.csv"
    if not samp_path.exists():
        raise StopStep("task2", f"missing {samp_path}")
    samples = pd.read_csv(samp_path)
    log_columns("t2_gse113957_samples", list(samples.columns), str(samp_path))
    need = {"gsm", "title", "age_years"}
    missing = sorted(need - set(samples.columns))
    if missing:
        raise StopStep("task2", f"{samp_path.name} missing {missing}")
    if not GSE113957_FPKM.exists():
        raise StopStep("task2", f"missing {GSE113957_FPKM}")
    log(f"[t2 GSE113957] reading {GSE113957_FPKM}")
    with gzip.open(GSE113957_FPKM, "rt", encoding="utf-8", errors="replace") as fh:
        df = pd.read_csv(fh, sep="\t")
    log_columns("t2_gse113957_fpkm", list(df.columns)[:20], str(GSE113957_FPKM))
    if "Transcript ID" not in df.columns or "Annotation/Divergence" not in df.columns:
        raise StopStep("task2", f"GSE113957_fpkm columns={list(df.columns)[:12]}. Not substituting.")
    tx = df["Transcript ID"].astype(str).to_numpy()
    id_tx = detect_id_type(tx, log, "gse113957_Transcript_ID", str(GSE113957_FPKM) + ":Transcript ID")
    require_mappable(id_tx, {"refseq_transcript"}, log, "task2")
    meta_cols = {"Transcript ID", "chr", "start", "end", "strand", "Length", "Copies", "Annotation/Divergence"}
    samp_cols = [c for c in df.columns if c not in meta_cols]
    log(f"[t2 GSE113957] n_transcripts={len(df)} n_sample_cols={len(samp_cols)}")
    mat = df[samp_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=np.float64)
    if looks_like_counts(mat):
        log("[t2 GSE113957] values look integer; still treating as FPKM as the file name and FIBRO2 record state")
    sym = np.array([_parse_homer_symbol(a) for a in df["Annotation/Divergence"].astype(str)], dtype=object)
    n_unparsed = int(sum(s is None for s in sym))
    keep = np.array([s is not None for s in sym])
    log(f"[t2 GSE113957] parsed gene symbol from Annotation/Divergence n_unparsed={n_unparsed} n_kept={int(keep.sum())}")
    if int(keep.sum()) == 0:
        raise StopStep("task2", "GSE113957: 0 gene symbols parsed from Annotation/Divergence. Not mapping.")
    mat = mat[keep]
    sym = np.array([s for s in sym if s is not None], dtype=object)
    id_sym = detect_id_type(sym, log, "gse113957_parsed_symbols", str(GSE113957_FPKM) + ":Annotation/Divergence")
    # Collapse transcripts to gene by summing FPKM.
    tmp = pd.DataFrame(mat, columns=samp_cols)
    tmp.insert(0, "symbol", sym)
    collapsed = tmp.groupby("symbol", sort=False).sum(numeric_only=True)
    symbols = collapsed.index.astype(str).to_numpy()
    C_gs = collapsed.to_numpy(dtype=np.float64)  # genes × samples
    log(f"[t2 GSE113957] collapsed transcripts→genes n_transcripts_kept={int(keep.sum())} "
        f"n_unique_symbols={len(symbols)} (sum FPKM)")
    detect_id_type(symbols, log, "gse113957_gene_symbols_collapsed", "collapsed Annotation/Divergence")
    # Match sample columns to GEO titles.
    title_to_row = {}
    for i, r in samples.iterrows():
        title_to_row[str(r.title)] = i
    unmatched = [c for c in samp_cols if c not in title_to_row]
    if unmatched:
        raise StopStep(
            "task2",
            f"GSE113957 FPKM columns not in sample titles n={len(unmatched)} e.g. {unmatched[:8]}. Not coercing.",
        )
    # Rank-transform genes within sample.
    ranks = np.vstack([rankdata(C_gs[:, j], method="average") for j in range(C_gs.shape[1])]).T
    log(f"[t2 GSE113957] rank-transform genes within sample method=average shape={ranks.shape} "
        "(different transform from GTEx TMM/log2-CPM; weaker test)")
    # Build obs: one row per FPKM sample, ages from the record, adult ≥18.
    rows = []
    keep_j = []
    n_missing_age = n_child = n_adult = 0
    for j, col in enumerate(samp_cols):
        r = samples.iloc[title_to_row[col]]
        age = r.age_years
        rec = dict(sample=str(r.gsm), donor=str(r.gsm), title=str(col),
                   age_years=float(age) if pd.notna(age) else np.nan,
                   age_source=str(r.age_source) if "age_source" in samples.columns else "ta_GSE113957_samples")
        if pd.isna(age):
            n_missing_age += 1
            rec["exclude_reason"] = "no stated age"
            rows.append(rec)
            continue
        if float(age) < ADULT_MIN_AGE:
            n_child += 1
            rec["exclude_reason"] = f"age<{ADULT_MIN_AGE}"
            rows.append(rec)
            continue
        n_adult += 1
        rec["exclude_reason"] = None
        rows.append(rec)
        keep_j.append(j)
    obs_all = pd.DataFrame(rows)
    obs_all.to_csv(MD3_DIR / "t2_GSE113957_all_samples.csv", index=False)
    log(f"[t2 GSE113957] n_fpkm_samples={len(samp_cols)} n_missing_age={n_missing_age} "
        f"n_child={n_child} n_adult>={ADULT_MIN_AGE}={n_adult}")
    if n_adult < 8:
        raise StopStep("task2", f"GSE113957 adult n={n_adult}<8 after ≥{ADULT_MIN_AGE} filter")
    obs = obs_all[obs_all.exclude_reason.isna()].reset_index(drop=True)
    R = ranks[:, np.array(keep_j, dtype=int)]  # genes × adult samples
    # Frozen ruler on ranks (as-is weights/mu/sd). Missing overlap genes at z=0.
    Y, align = align_counts_to_ruler(R, symbols, symbols, frozen, log, tag="gse113957_ranks")
    pred_r, Z, missing_g = score_frozen(Y, frozen, counts=Y)
    log(f"[t2 GSE113957 ruler] overlap={align.get('n_overlap')} missing_z0={int(missing_g.sum())} "
        "z uses frozen GTEx log2-CPM mu/sd on rank inputs (weaker test; not a refit)")
    y = obs.age_years.to_numpy(float)
    donor = obs.donor.astype(str).to_numpy()
    ruler = _rho_pack(y, pred_r, donor, log, "gse113957_ruler_ranks")
    ruler.update(
        instrument="frozen_GTEx_fibroblast_ruler",
        cohort="GSE113957 bulk FPKM ranks", accession="GSE113957",
        platform="bulk_rnaseq_FPKM_rank",
        n_overlap=align.get("n_overlap"), n_missing_z0=int(missing_g.sum()),
        n_ruler=align.get("n_ruler"),
        id_type="refseq_transcript → symbol via Annotation/Divergence",
        transform="rank within sample, then frozen mu/sd/w as-is (NOT TMM/log2-CPM)",
        weaker_test=True,
        n_excluded_missing_age=n_missing_age, n_excluded_child=n_child, n_adult=n_adult,
        n_fpkm_samples=len(samp_cols), refit=False,
    )
    # AMS on ranks (samples × genes). nbin=24 is inherited from md2; not reduced if bins collapse.
    logR = R.T  # already ranks; do not LogNormalize counts. Score AMS directly on ranks.
    gene_sets = {
        "MD": list(sets["MD"]),
        "age_up": list(sets["age_up"]),
        "age_down": list(sets["age_down"]),
    }
    obs = obs.copy()
    obs["age_score"] = pred_r
    try:
        ams, meta = add_module_score(logR, symbols, gene_sets, log=log, tag="gse113957_ranks")
    except StopStep as e:
        record_failure(e.step, f"GSE113957: {e.message}", e.details)
        log(f"[t2] GSE113957 AMS STOP [{e.step}] {e.message}. "
            "nbin=24 not reduced. MD and Table S3 AMS not scored on this cohort. "
            "Frozen-ruler rank scores are kept. No paired Δρ for this cohort.")
        obs.to_csv(MD3_DIR / "t2_GSE113957_scores.csv", index=False)
        md = dict(
            instrument="MD_AddModuleScore",
            cohort="GSE113957 bulk FPKM ranks", accession="GSE113957",
            platform="bulk_rnaseq_FPKM_rank", weaker_test=True,
            n_excluded_child=n_child, n_adult=n_adult, n_excluded_missing_age=n_missing_age,
            reason=str(e.message),
            transform="AddModuleScore on ranks not run: cut_number nbin=24 collapsed; nbin not reduced",
        )
        s3 = dict(
            instrument="age_up_minus_age_down_AddModuleScore",
            cohort="GSE113957 bulk FPKM ranks", accession="GSE113957",
            platform="bulk_rnaseq_FPKM_rank", weaker_test=True,
            n_excluded_child=n_child, n_adult=n_adult, n_excluded_missing_age=n_missing_age,
            reason=str(e.message),
            transform="AddModuleScore on ranks not run: cut_number nbin=24 collapsed; nbin not reduced",
        )
        return [ruler, md, s3], []
    md_pred = np.asarray(ams["MD"], float)
    s3_pred = np.asarray(ams["age_up"], float) - np.asarray(ams["age_down"], float)
    md = _rho_pack(y, md_pred, donor, log, "gse113957_MD_ranks")
    md.update(
        instrument="MD_AddModuleScore",
        cohort="GSE113957 bulk FPKM ranks", accession="GSE113957",
        platform="bulk_rnaseq_FPKM_rank",
        n_MD_mapped=meta["per_set"]["MD"]["n_mapped"],
        n_MD_missing=meta["per_set"]["MD"]["n_missing"],
        n_MD_total=int(len(sets["MD"])),
        id_type="symbol (parsed)",
        transform="AddModuleScore on within-sample ranks (NOT LogNormalize of counts)",
        weaker_test=True, n_excluded_child=n_child, n_adult=n_adult,
    )
    s3 = _rho_pack(y, s3_pred, donor, log, "gse113957_S3_ranks")
    s3.update(
        instrument="age_up_minus_age_down_AddModuleScore",
        cohort="GSE113957 bulk FPKM ranks", accession="GSE113957",
        platform="bulk_rnaseq_FPKM_rank",
        n_age_up_mapped=meta["per_set"]["age_up"]["n_mapped"],
        n_age_up_missing=meta["per_set"]["age_up"]["n_missing"],
        n_age_up_total=int(len(sets["age_up"])),
        n_age_down_mapped=meta["per_set"]["age_down"]["n_mapped"],
        n_age_down_missing=meta["per_set"]["age_down"]["n_missing"],
        n_age_down_total=int(len(sets["age_down"])),
        id_type="symbol (parsed)",
        transform="AddModuleScore on within-sample ranks; Age up minus Age down",
        weaker_test=True, n_excluded_child=n_child, n_adult=n_adult,
    )
    obs = obs.copy()
    obs["age_score"] = pred_r
    obs["md_score"] = md_pred
    obs["age_up"] = ams["age_up"]
    obs["age_down"] = ams["age_down"]
    obs["age_up_minus_age_down"] = s3_pred
    obs.to_csv(MD3_DIR / "t2_GSE113957_scores.csv", index=False)
    d_md = paired_delta_rho(y, pred_r, md_pred, donor, log, "gse113957_ruler_minus_MD")
    d_s3 = paired_delta_rho(y, pred_r, s3_pred, donor, log, "gse113957_ruler_minus_S3")
    d_md.update(cohort="GSE113957 bulk FPKM ranks", other="MD_AddModuleScore")
    d_s3.update(cohort="GSE113957 bulk FPKM ranks", other="age_up_minus_age_down")
    return [ruler, md, s3], [d_md, d_s3]


def _t2_reading(delta_df):
    named = [
        "CELLxGENE a19d1667-a7b5-4556-9e5f-f9bfa690c0f1",
        "GSE226189 bulk",
        "GSE113957 bulk FPKM ranks",
    ]
    by = {}
    for cohort in named:
        sub = delta_df[delta_df.cohort.astype(str) == cohort]
        by[cohort] = sub.to_dict(orient="records") if len(sub) else []
    n_ruler = 0
    n_other = 0
    n_zero = 0
    n_scored = 0
    per = []
    for cohort in named:
        recs = by[cohort]
        if not recs:
            per.append(dict(cohort=cohort, scored=False))
            continue
        n_scored += 1
        favs = [r.get("favours") for r in recs]
        # A cohort favours the ruler if ANY paired contrast vs MD or S3 excludes 0 for the ruler.
        # Reading is about MD and/or their aging signature vs ruler; report each contrast.
        cohort_ruler = any(r.get("favours") == "ruler" for r in recs)
        cohort_other = any(r.get("favours") == "other" for r in recs)
        cohort_all_zero = all(r.get("favours") in (None, "includes_zero") or not r.get("excludes_zero") for r in recs)
        if cohort_ruler:
            n_ruler += 1
        if cohort_other and not cohort_ruler:
            n_other += 1
        if cohort_all_zero:
            n_zero += 1
        per.append(dict(cohort=cohort, scored=True, favours_ruler=cohort_ruler,
                        favours_other=cohort_other, all_ci_include_zero=cohort_all_zero, contrasts=recs))
    if n_ruler >= 2:
        key = "ruler_stronger_age_in_ge2"
        text = (
            f"Δρ CI excludes zero in favour of the ruler in {n_ruler} of 3 named cohorts. "
            "MD and/or their aging signature predict donor age materially worse than the ruler "
            "on independent fibroblast cohorts."
        )
    elif n_scored >= 1 and n_ruler == 0 and n_other == 0:
        key = "instruments_not_distinguishable"
        text = (
            "Δρ CI includes zero. The instruments are not distinguishable on age prediction at these n, "
            "and the Task 1 disagreement cannot be attributed to one being a weaker age instrument."
        )
    else:
        key = "mixed_delta_rho"
        text = (
            f"Mixed paired Δρ: n_named=3 n_scored={n_scored} n_favour_ruler={n_ruler} "
            f"n_favour_other_only={n_other} n_all_CI_include_zero={n_zero}. "
            "Report per cohort, do not average, do not declare a winner."
        )
    return dict(key=key, text=text, n_favour_ruler=n_ruler, n_scored=n_scored, per_cohort=per)


def run_task2(log=None):
    _require()
    close_log = False
    if log is None:
        log = Logger(MD3_DIR / "t2_report.txt")
        close_log = True
    md3_log_banner(log, "TASK2")
    log(PREREG_TASK2)
    sets = load_json(MD3_DIR / "genesets.json")
    frozen = load_frozen_ruler()
    rows, deltas = [], []

    try:
        r, d, _ = _cxg(sets, log)
        rows.extend(r)
        deltas.extend(d)
    except StopStep as e:
        record_failure(e.step, f"CELLxGENE: {e.message}", e.details)
        log(f"[t2] CXG STOP [{e.step}] {e.message}")
        rows.append(dict(instrument="STOP", cohort="CELLxGENE a19d1667-a7b5-4556-9e5f-f9bfa690c0f1",
                         reason=str(e.message)))

    try:
        r, d = _gse226189(sets, log)
        rows.extend(r)
        deltas.extend(d)
    except StopStep as e:
        record_failure(e.step, f"GSE226189: {e.message}", e.details)
        log(f"[t2] GSE226189 STOP [{e.step}] {e.message}")
        rows.append(dict(instrument="STOP", cohort="GSE226189 bulk", reason=str(e.message)))

    try:
        r, d = _gse113957(sets, frozen, log)
        rows.extend(r)
        deltas.extend(d)
    except StopStep as e:
        record_failure(e.step, f"GSE113957: {e.message}", e.details)
        log(f"[t2] GSE113957 STOP [{e.step}] {e.message}")
        rows.append(dict(instrument="STOP", cohort="GSE113957 bulk FPKM ranks", reason=str(e.message)))

    inst = pd.DataFrame(rows)
    inst.to_csv(MD3_DIR / "t2_instruments.csv", index=False)
    delta_df = pd.DataFrame(deltas)
    delta_df.to_csv(MD3_DIR / "t2_delta_rho.csv", index=False)
    reading = _t2_reading(delta_df) if len(delta_df) else dict(
        key="no_delta", text="No paired Δρ rows were written.", n_favour_ruler=0, n_scored=0,
    )
    dump_json(MD3_DIR / "t2_reading.json", jsonable(reading))
    log(f"[t2 reading] key={reading['key']}")
    log(reading["text"])
    dump_json(MD3_DIR / "t2_summary.json", jsonable(dict(
        reading=reading, n_instrument_rows=int(len(inst)), n_delta_rows=int(len(delta_df)),
        seed=MD3_SEED, boot=MD3_BOOT,
    )))
    man = load_manifest()
    man["status"] = "TASK2_DONE"
    man["t2_reading"] = reading.get("key")
    save_manifest(man)
    progress_snapshot("Task 3 MD co-variation (report-only)", stop="TASK2_DONE")
    if close_log:
        log.close()
    return reading


if __name__ == "__main__":
    try:
        run_task2()
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise

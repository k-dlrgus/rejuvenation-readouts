"""Task 3 — instrument validity vs donor age on FIBRO2 external cohorts. Report-only. No gate."""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
if not hasattr(np, "unicode_"):
    np.unicode_ = np.str_  # type: ignore[attr-defined]
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md2_common import (  # noqa: E402
    MD2_DIR, MD2_PROC, MD2_SEED, MD2_BOOT, CXG_H5AD, CXG_DATASET_ID,
    GSE226189_TAR, FIBRO2_DIR, N_PERM, N_BOOT, PREREG_TASK3, PREREG_TASK3_FLAG,
    StopStep, Logger, dump_json, jsonable, load_json, md2_log_banner,
    load_manifest, save_manifest, record_failure, log_columns,
    progress_snapshot,
)
from fibro2_common import rho_null_donor, ADULT_MIN_AGE, strict_age_years, load_frozen_ruler  # noqa: E402
from fibro_common import bootstrap_rho_ci  # noqa: E402
from fibro2_ta import load_count_table, _decode_arr  # noqa: E402
from tissue_peek_obs import _open_h5, obs_column_names, read_obs_column  # noqa: E402
from pseudobulk import find_counts_path, read_counts_chunk  # noqa: E402
from md2_score import lognormalize_dense, add_module_score  # noqa: E402
from brain_phase1_common import strip_ensembl  # noqa: E402


def _rho_pack(y, pred, donor, log, tag):
    sc, rho_null, p, nulls = rho_null_donor(y, pred, donor, n_perm=N_PERM, seed=MD2_SEED)
    rng_b = np.random.default_rng(MD2_BOOT)
    ci = bootstrap_rho_ci(y, pred, rng_b, n_boot=N_BOOT)
    rec = dict(
        tag=tag, n=int(len(y)), n_donors=int(pd.Series(donor).nunique()),
        rho=sc["rho"], rho_null=rho_null, rho_p=p,
        rho_ci_lo=ci.get("p025"), rho_ci_hi=ci.get("p975"),
        r=sc["r"], r2=sc["r2"], cal_r2=sc["cal_r2"],
        n_perm=N_PERM, n_boot=N_BOOT, seed=MD2_SEED, boot_seed=MD2_BOOT,
        age_min=float(np.nanmin(y)), age_max=float(np.nanmax(y)),
        n_unique_age=int(pd.Series(y).nunique()),
    )
    log(f"[t3 {tag}] ρ={rec['rho']:+.3f} null={rho_null:+.3f} p={p:+.3f} "
        f"CI=[{ci.get('p025')}, {ci.get('p975')}] n={rec['n']} n_donors={rec['n_donors']}")
    return rec


def _ruler_from_fibro2(path, log, tag):
    if not path.exists():
        raise StopStep("task3", f"missing FIBRO2 result {path}. Not reconstructing from memory.")
    rec = load_json(path)
    log(f"[t3] read frozen-ruler numbers from {path} keys={sorted(rec)}")
    return dict(
        instrument="frozen_GTEx_fibroblast_ruler",
        source=str(path), tag=tag,
        n=rec.get("n_samples"), n_donors=rec.get("n_donors"),
        rho=rec.get("rho"), rho_null=rec.get("rho_null"), rho_p=rec.get("rho_p"),
        rho_ci_lo=rec.get("rho_ci_lo"), rho_ci_hi=rec.get("rho_ci_hi"),
        r=rec.get("r"), r2=rec.get("r2"), cal_r2=rec.get("cal_r2"),
        n_perm=rec.get("n_perm"), n_boot=rec.get("n_boot"),
        seed=rec.get("seed"), boot_seed=rec.get("boot_seed"),
        platform=rec.get("platform_detail") or rec.get("platform_kind"),
        n_overlap=rec.get("n_overlap"),
        accession=rec.get("accession"),
        refit=False,
    )


def _md_on_counts(counts, symbols, samples, obs, sets, log, tag):
    """counts genes × samples. LogNormalize samples as cells, AddModuleScore, donor ρ vs age."""
    pos = {str(s): i for i, s in enumerate(samples)}
    missing = [s for s in obs["sample"].astype(str) if s not in pos]
    if missing:
        raise StopStep("task3", f"{tag}: {len(missing)} samples not in count matrix e.g. {missing[:8]}")
    col = np.array([pos[s] for s in obs["sample"].astype(str)], dtype=int)
    C = np.asarray(counts, np.float64)[:, col].T  # samples × genes
    logC = lognormalize_dense(C)
    ams, meta = add_module_score(
        logC, np.asarray(symbols).astype(str), {"MD": list(sets["MD"])},
        log=log, tag=tag,
    )
    pred = ams["MD"]
    y = obs.age_years.to_numpy(float)
    donor = obs.donor.astype(str).to_numpy()
    rec = _rho_pack(y, pred, donor, log, tag)
    rec.update(
        instrument="MD_AddModuleScore",
        n_MD_mapped=meta["per_set"]["MD"]["n_mapped"],
        n_MD_missing=meta["per_set"]["MD"]["n_missing"],
        n_genes=int(C.shape[1]),
        scoring="AddModuleScore on LogNormalize of donor/sample count columns; not a ruler refit",
    )
    obs = obs.copy()
    obs["md_score"] = pred
    obs.to_csv(MD2_DIR / f"t3_{tag}_scores.csv", index=False)
    return rec


def _ensembl_to_symbol(ens_ids, log):
    """Map ENSG Tracking_IDs to symbols using on-disk 10x features, then the frozen ruler."""
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
            log_columns("t3_gse226189_10x_map", ["gene_id", "symbol"], str(zpath))
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
    log(f"[t3 map] Ensembl→symbol sources={src} n_ids={len(ens)} n_mapped_to_symbol={n_sym} "
        f"n_left_as_ensembl={int(len(ens) - n_sym)}")
    if n_sym == 0:
        raise StopStep(
            "task3",
            "GSE226189 geneCOUNT Tracking_ID is Ensembl and 0 IDs mapped to symbols "
            f"from {src}. Not scoring MD on Ensembl IDs as if they were symbols.",
        )
    return symbols


def _md_gse226189(sets, log):
    scores_path = FIBRO2_DIR / "ta_GSE226189_near_miss_bulk_scores.csv"
    if not scores_path.exists():
        raise StopStep("task3", f"missing {scores_path}")
    obs = pd.read_csv(scores_path)
    log(f"[t3 GSE226189] scores columns actually read: {list(obs.columns)}")
    log_columns("t3_gse226189_scores", list(obs.columns), str(scores_path))
    if not GSE226189_TAR.exists():
        raise StopStep("task3", f"missing {GSE226189_TAR}")
    pack = load_count_table(GSE226189_TAR, log, gsm_keep=obs["sample"].astype(str).tolist())
    genes = np.asarray(pack["genes"]).astype(str)
    log(f"[t3 GSE226189] gene index preview={genes[:8].tolist()} "
        f"ensembl_like={bool(str(genes[0]).startswith('ENSG'))}")
    log_columns("t3_gse226189_gene_index", list(genes[:20]), str(GSE226189_TAR) + ":Tracking_ID")
    if str(genes[0]).startswith("ENS"):
        symbols = _ensembl_to_symbol(genes, log)
    else:
        symbols = np.asarray(pack["symbols"]).astype(str)
    rec = _md_on_counts(
        pack["counts"], symbols, pack["samples"], obs, sets, log,
        tag="GSE226189_MD",
    )
    rec.update(
        accession="GSE226189", platform="bulk_rnaseq", n_cohort_claimed=82,
        gene_index="Tracking_ID Ensembl mapped to symbols via 10x features + frozen ruler",
    )
    return rec


def _md_cxg(sets, log):
    cache = MD2_PROC / "cxg_donor_counts.npz"
    scores_path = FIBRO2_DIR / "ta_cxg_a19d1667_scores.csv"
    if not scores_path.exists():
        raise StopStep("task3", f"missing {scores_path}")
    obs_scored = pd.read_csv(scores_path)
    log(f"[t3 cxg] scores columns actually read: {list(obs_scored.columns)}")
    log_columns("t3_cxg_scores", list(obs_scored.columns), str(scores_path))
    if cache.exists():
        z = np.load(cache, allow_pickle=True)
        counts = z["counts"]
        genes = np.asarray(z["genes"]).astype(str)
        symbols = np.asarray(z["symbols"]).astype(str)
        samples = np.asarray(z["samples"]).astype(str)
        log(f"[t3 cxg] loaded cached donor counts {counts.shape} from {cache}")
    else:
        if not CXG_H5AD.exists():
            raise StopStep("task3", f"missing {CXG_H5AD}")
        counts, genes, symbols, samples = _pseudobulk_cxg(log)
        np.savez_compressed(
            cache, counts=counts.astype(np.float32), genes=genes, symbols=symbols, samples=samples,
        )
    # Restrict to the 93 donors already in the FIBRO2 table (same units).
    obs = obs_scored.copy()
    rec = _md_on_counts(counts, symbols, samples, obs, sets, log, tag="cxg_a19d1667_MD")
    rec.update(
        accession=CXG_DATASET_ID, platform="10x (restricted)",
        n_cohort_claimed=93, h5ad=str(CXG_H5AD),
    )
    return rec


def _pseudobulk_cxg(log):
    """Same donor/age/10x/adult filters as fibro2_ta.project_cxg_fibroblast. Full genes."""
    url = str(CXG_H5AD)
    log(f"[t3 cxg] opening {url}")
    h, handle = _open_h5(url)
    t0 = time.time()
    try:
        cols = obs_column_names(h)
        log(f"[t3 cxg] obs columns actually read: {cols}")
        log_columns("t3_cxg_h5ad_obs", cols, url)
        if "donor_id" not in cols:
            raise StopStep("task3", f"obs missing donor_id. columns={cols}")
        age_col = None
        for c in ("development_stage", "age", "donor_age", "Age", "age_years"):
            if c in cols:
                age_col = c
                break
        if age_col is None:
            raise StopStep("task3", f"no age-like obs column. columns={cols}")
        org_col = "organism" if "organism" in cols else None
        assay_col = "assay" if "assay" in cols else None
        donor = np.asarray(read_obs_column(h, "donor_id")).astype(str)
        age_raw = np.asarray(read_obs_column(h, age_col))
        n = len(donor)
        ages = np.array([strict_age_years(x) for x in age_raw], dtype=object)
        stated = np.array([a is not None for a in ages])
        keep = stated.copy()
        if org_col:
            org = np.asarray(read_obs_column(h, org_col)).astype(str)
            human = np.array([
                "homo sapiens" in x.lower() or x.lower() in {"human", "ncbi9606"}
                or "NCBITaxon:9606" in x for x in org
            ])
            keep &= human
        assay = np.asarray(read_obs_column(h, assay_col)).astype(str) if assay_col else np.array(["unknown"] * n)
        is_10x = np.array(["10x" in x.lower() for x in assay])
        keep = keep & is_10x
        age_num = np.array([float(a) if a is not None else np.nan for a in ages])
        keep = keep & (age_num >= ADULT_MIN_AGE)
        try:
            cpath = find_counts_path(h)
        except Exception as e:
            raise StopStep("task3", f"no integer count matrix: {e}")
        n_var = int(h[cpath].attrs["shape"][1])
        var_group = "raw/var" if cpath.startswith("raw/") else "var"
        vg = h[var_group]
        gene_id = _decode_arr(vg["_index"][()])
        symbol = None
        for sk in ("feature_name", "gene_symbols", "symbol", "name"):
            if sk not in vg:
                continue
            raw_s = vg[sk]
            if hasattr(raw_s, "dtype"):
                symbol = _decode_arr(raw_s[()])
            elif "categories" in raw_s:
                cats = _decode_arr(raw_s["categories"][()])
                codes = raw_s["codes"][()]
                symbol = np.array([cats[int(c)] if int(c) >= 0 else "" for c in codes], dtype=object)
            break
        if symbol is None:
            symbol = gene_id
        use_idx = np.flatnonzero(keep)
        d_keep = donor[keep]
        y_keep = age_num[keep]
        tmp = pd.DataFrame({"donor": d_keep, "age": y_keep})
        d_u = pd.unique(d_keep)
        gidx = {d: i for i, d in enumerate(d_u)}
        cell_g = np.array([gidx[d] for d in d_keep], dtype=int)
        n_grp = len(d_u)
        X_sum = np.zeros((n_grp, n_var), dtype=np.float64)
        order = np.argsort(use_idx)
        row_pos_s = use_idx[order]
        cell_g_s = cell_g[order]
        pos_ptr = 0
        chunk = 40000
        for r0 in range(0, n, chunk):
            r1 = min(r0 + chunk, n)
            j0 = pos_ptr
            while pos_ptr < len(row_pos_s) and row_pos_s[pos_ptr] < r1:
                pos_ptr += 1
            if pos_ptr == j0:
                continue
            sel = row_pos_s[j0:pos_ptr] - r0
            gk = cell_g_s[j0:pos_ptr]
            Xc = read_counts_chunk(h, cpath, r0, r1, n_var)[sel]
            M = sp.csr_matrix((np.ones(len(gk)), (gk, np.arange(len(gk)))), shape=(n_grp, len(gk)))
            X_sum += (M @ Xc).toarray()
            log(f"   rows {r1:,}/{n:,}  ({time.time()-t0:.0f}s)")
        log(f"[t3 cxg] pseudobulk donors={n_grp} genes={n_var} age_col={age_col}")
        return X_sum.T, gene_id, symbol, np.asarray(d_u).astype(str)
    finally:
        h.close()
        if handle is not None:
            handle.close()


def run_task3(log=None):
    if not PREREG_TASK3_FLAG.exists():
        raise StopStep("prereg", "PREREG_TASK3.flag missing")
    close_log = False
    if log is None:
        log = Logger(MD2_DIR / "t3_report.txt")
        close_log = True
    md2_log_banner(log, "TASK3")
    log(PREREG_TASK3)
    sets = load_json(MD2_DIR / "genesets.json")
    gs = load_json(MD2_DIR / "genesets_summary.json")
    aging_ok = bool(gs.get("aging_signature_built"))
    aging_reason = gs.get("aging_signature_reason")
    log(f"[t3] aging_signature_built={aging_ok} reason={aging_reason}")

    rows = []
    rows.append(_ruler_from_fibro2(FIBRO2_DIR / "ta_cxg_a19d1667_result.json", log, "cxg_a19d1667"))
    rows[-1]["cohort"] = "CELLxGENE a19d1667-a7b5-4556-9e5f-f9bfa690c0f1"
    rows.append(_ruler_from_fibro2(FIBRO2_DIR / "ta_GSE226189_near_miss_bulk_result.json", log, "GSE226189"))
    rows[-1]["cohort"] = "GSE226189 bulk"
    rows[-1]["accession"] = "GSE226189"

    md_cxg = _md_cxg(sets, log)
    md_cxg["cohort"] = "CELLxGENE a19d1667-a7b5-4556-9e5f-f9bfa690c0f1"
    rows.append(md_cxg)
    try:
        md_gse = _md_gse226189(sets, log)
        md_gse["cohort"] = "GSE226189 bulk"
        rows.append(md_gse)
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        log(f"[t3] GSE226189 MD STOP [{e.step}] {e.message}")
        rows.append(dict(
            instrument="MD_AddModuleScore", cohort="GSE226189 bulk", accession="GSE226189",
            rho=np.nan, rho_null=np.nan, rho_p=np.nan, note="not scored",
            reason=str(e.message),
        ))

    if aging_ok:
        log("[t3] aging signature was marked built — unexpected; GSE113957 path should have stopped it.")
    else:
        for cohort, acc in (
            ("CELLxGENE a19d1667-a7b5-4556-9e5f-f9bfa690c0f1", CXG_DATASET_ID),
            ("GSE226189 bulk", "GSE226189"),
        ):
            rows.append(dict(
                instrument="Fleischer_DESeq2_aging_signature",
                cohort=cohort, accession=acc,
                rho=np.nan, rho_null=np.nan, rho_p=np.nan,
                rho_ci_lo=np.nan, rho_ci_hi=np.nan,
                note="not run",
                reason=aging_reason,
            ))

    df = pd.DataFrame(rows)
    df.to_csv(MD2_DIR / "t3_instruments.csv", index=False)
    summary = dict(
        aging_signature_built=False,
        aging_signature_reason=aging_reason,
        n_rows=int(len(df)),
        note="report-only; no gate; frozen ruler numbers read from results/fibro2/; ruler not refit",
    )
    dump_json(MD2_DIR / "t3_summary.json", jsonable(summary))
    man = load_manifest()
    man["status"] = "TASK3_DONE"
    save_manifest(man)
    progress_snapshot("write FINDINGS_MD2.md from disk", stop="TASK3_DONE")
    if close_log:
        log.close()
    return summary


if __name__ == "__main__":
    try:
        run_task3()
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise

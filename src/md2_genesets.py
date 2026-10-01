"""Build MD / TGF-β / aging / reprogramming gene sets exactly as published. Fail loudly."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md2_common import (  # noqa: E402
    MD2_DIR, MMC3_XLSX, HALLMARK_LOCAL, GSE113957_FPKM, MD2_SEED,
    HALLMARK_EMT_NAME, HALLMARK_TGFB_NAME, HALLMARK_EMT_MSIGDB, HALLMARK_TGFB_MSIGDB,
    MD_TF_ADD, STAR_METHODS_QUOTES,
    StopStep, Logger, dump_json, jsonable, load_manifest, save_manifest,
    log_columns, log_geneset, md2_log_banner, progress_snapshot, record_failure,
)
from md_stage1 import _xlsx_columns  # noqa: E402
from fibro2_common import looks_like_counts  # noqa: E402

MSIGDB_GMT_CANDIDATES = (
    "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2020.1.Hs/h.all.v2020.1.Hs.symbols.gmt",
    "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt",
    "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt",
)


def _parse_enrichr_hallmark(path: Path) -> dict:
    sets = {}
    text = path.read_text(encoding="utf-8")
    for line in text.splitlines():
        if not line.strip():
            continue
        parts = [p for p in line.split("\t") if p != ""]
        if not parts:
            continue
        name = parts[0].strip()
        genes = []
        for p in parts[1:]:
            g = p.strip()
            if not g or g.lower().startswith("http"):
                continue
            genes.append(g.upper())
        # unique, preserve order
        seen = set()
        out = []
        for g in genes:
            if g not in seen:
                seen.add(g)
                out.append(g)
        sets[name] = out
    return sets


def _try_download_gmt(log) -> dict:
    dest = MD2_DIR / "h.all.downloaded.gmt"
    headers = {"User-Agent": "AgeIdentitySeparabilityBenchmark/md2"}
    errors = []
    for url in MSIGDB_GMT_CANDIDATES:
        try:
            r = requests.get(url, headers=headers, timeout=30)
            log(f"[msigdb] GET {url} status={r.status_code} nbytes={len(r.content)}")
            if r.status_code != 200:
                errors.append(dict(url=url, status=r.status_code))
                continue
            text = r.text
            if "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" not in text:
                errors.append(dict(url=url, status=r.status_code, note="no HALLMARK_EMT line"))
                continue
            dest.write_text(text, encoding="utf-8")
            sets = {}
            for line in text.splitlines():
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                name = parts[0]
                genes = [g.strip().upper() for g in parts[2:] if g.strip()]
                sets[name] = genes
            return dict(ok=True, url=url, path=str(dest), sets=sets, errors=errors)
        except Exception as e:
            errors.append(dict(url=url, error=str(e)))
            log(f"[msigdb] {url} failed: {e}")
    return dict(ok=False, url=None, path=None, sets={}, errors=errors)


def _set_compare(a, b):
    a = [str(x).upper() for x in a]
    b = [str(x).upper() for x in b]
    sa, sb = set(a), set(b)
    return dict(
        n_a=len(a), n_b=len(b), n_shared=len(sa & sb),
        only_a=sorted(sa - sb), only_b=sorted(sb - sa),
        equal=sa == sb,
    )


def build_gene_sets(log=None):
    close_log = False
    if log is None:
        log = Logger(MD2_DIR / "genesets_report.txt")
        close_log = True
    md2_log_banner(log, "GENESETS")
    dump_json(MD2_DIR / "star_methods_quotes.json", STAR_METHODS_QUOTES)
    log("[star] quotes written from paper extract; STAR Methods gene-set sentence:")
    log(STAR_METHODS_QUOTES["gene_sets"])

    if not HALLMARK_LOCAL.exists():
        raise StopStep("genesets", f"missing local Hallmark file {HALLMARK_LOCAL}")
    local_sets = _parse_enrichr_hallmark(HALLMARK_LOCAL)
    log(f"[hallmark-local] path={HALLMARK_LOCAL} n_sets={len(local_sets)} names={sorted(local_sets)}")
    log_columns("hallmark_local_set_names", sorted(local_sets), str(HALLMARK_LOCAL))
    if HALLMARK_EMT_NAME not in local_sets:
        raise StopStep(
            "genesets",
            f"local Hallmark file has no set {HALLMARK_EMT_NAME!r}. names={sorted(local_sets)}",
        )
    if HALLMARK_TGFB_NAME not in local_sets:
        raise StopStep(
            "genesets",
            f"local Hallmark file has no set {HALLMARK_TGFB_NAME!r}. names={sorted(local_sets)}",
        )
    emt_local = local_sets[HALLMARK_EMT_NAME]
    tgfb_local = local_sets[HALLMARK_TGFB_NAME]
    log_geneset("HALLMARK_EMT_local", HALLMARK_EMT_NAME, emt_local, str(HALLMARK_LOCAL))
    log_geneset("HALLMARK_TGFB_local", HALLMARK_TGFB_NAME, tgfb_local, str(HALLMARK_LOCAL))
    log(f"[hallmark-local] EMT n={len(emt_local)} TGFB n={len(tgfb_local)}")

    gmt = _try_download_gmt(log)
    emt_used = list(emt_local)
    tgfb_used = list(tgfb_local)
    msigdb_version = "MSigDB_Hallmark_2020.txt (Enrichr dump at results/lowdim/genesets/)"
    if gmt.get("ok"):
        sets = gmt["sets"]
        log(f"[msigdb-gmt] url={gmt['url']} n_sets={len(sets)}")
        log_columns("msigdb_gmt_set_names", sorted(sets), gmt["url"])
        if HALLMARK_EMT_MSIGDB not in sets:
            raise StopStep(
                "genesets",
                f"downloaded GMT missing {HALLMARK_EMT_MSIGDB}. names sample={sorted(sets)[:10]}",
            )
        if HALLMARK_TGFB_MSIGDB not in sets:
            raise StopStep(
                "genesets",
                f"downloaded GMT missing {HALLMARK_TGFB_MSIGDB}.",
            )
        emt_gmt = sets[HALLMARK_EMT_MSIGDB]
        tgfb_gmt = sets[HALLMARK_TGFB_MSIGDB]
        log_geneset("HALLMARK_EMT_gmt", HALLMARK_EMT_MSIGDB, emt_gmt, gmt["url"])
        log_geneset("HALLMARK_TGFB_gmt", HALLMARK_TGFB_MSIGDB, tgfb_gmt, gmt["url"])
        cmp_e = _set_compare(emt_local, emt_gmt)
        cmp_t = _set_compare(tgfb_local, tgfb_gmt)
        log(f"[msigdb-gmt] EMT local-vs-gmt equal={cmp_e['equal']} n_local={cmp_e['n_a']} "
            f"n_gmt={cmp_e['n_b']} only_local={cmp_e['only_a']} only_gmt={cmp_e['only_b']}")
        log(f"[msigdb-gmt] TGFB local-vs-gmt equal={cmp_t['equal']} n_local={cmp_t['n_a']} "
            f"n_gmt={cmp_t['n_b']} only_local={cmp_t['only_a']} only_gmt={cmp_t['only_b']}")
        # Official MSigDB symbols are the published source named in STAR Methods.
        emt_used = list(emt_gmt)
        tgfb_used = list(tgfb_gmt)
        msigdb_version = gmt["url"]
        gmt_cmp = dict(emt=cmp_e, tgfb=cmp_t)
    else:
        log(f"[msigdb-gmt] NOT RETRIEVED. using local Enrichr dump. errors={gmt.get('errors')}")
        gmt_cmp = None

    md_built = list(emt_used)
    already = {g.upper() for g in md_built}
    added, already_in = [], []
    for g in MD_TF_ADD:
        gu = g.upper()
        if gu in already:
            already_in.append(gu)
        else:
            md_built.append(gu)
            already.add(gu)
            added.append(gu)
    log(f"[MD] EMT n={len(emt_used)} plus TFs={list(MD_TF_ADD)} added={added} "
        f"already_in_EMT={already_in} final n={len(md_built)}")
    log_geneset("MD_built", "MD = HALLMARK_EMT + SNAI1,ZEB1,ZEB2,TWIST1,TWIST2",
                md_built, msigdb_version)
    log_geneset("TGFB_built", HALLMARK_TGFB_MSIGDB, tgfb_used, msigdb_version)

    if not MMC3_XLSX.exists():
        raise StopStep("genesets", f"missing mmc3.xlsx {MMC3_XLSX}")
    sheets = _xlsx_columns(MMC3_XLSX)
    log(f"[mmc3] path={MMC3_XLSX} sheets={list(sheets)} "
        f"columns_per_sheet={ {s: list(cols) for s, cols in sheets.items()} }")
    for sname, cols in sheets.items():
        log_columns(f"mmc3_{sname}", list(cols), str(MMC3_XLSX))
        for cname, vals in cols.items():
            log_geneset(f"mmc3_{sname}_{cname}", f"{sname}:{cname}", vals, str(MMC3_XLSX))

    md_sheet = sheets.get("MD_signatures") or {}
    if "MD score" not in md_sheet:
        raise StopStep("genesets", f"mmc3 MD_signatures missing column 'MD score'. cols={list(md_sheet)}")
    if "TGFB score" not in md_sheet:
        raise StopStep("genesets", f"mmc3 MD_signatures missing column 'TGFB score'. cols={list(md_sheet)}")
    mmc3_md = [g.strip().upper() for g in md_sheet["MD score"] if g.strip()]
    mmc3_tgfb = [g.strip().upper() for g in md_sheet["TGFB score"] if g.strip()]
    cmp_md = _set_compare(md_built, mmc3_md)
    cmp_tgfb = _set_compare(tgfb_used, mmc3_tgfb)
    log(f"[mmc3-crosscheck] MD built-vs-mmc3 equal={cmp_md['equal']} n_built={cmp_md['n_a']} "
        f"n_mmc3={cmp_md['n_b']} only_built={cmp_md['only_a']} only_mmc3={cmp_md['only_b']}")
    log(f"[mmc3-crosscheck] TGFB built-vs-mmc3 equal={cmp_tgfb['equal']} n_built={cmp_tgfb['n_a']} "
        f"n_mmc3={cmp_tgfb['n_b']} only_built={cmp_tgfb['only_a']} only_mmc3={cmp_tgfb['only_b']}")
    if not cmp_md["equal"]:
        log("[mmc3-crosscheck] MD mismatch recorded. Not reconciling silently. Scoring uses STAR Methods built set.")
    if not cmp_tgfb["equal"]:
        log("[mmc3-crosscheck] TGFB mismatch recorded. Not reconciling silently. Scoring uses STAR Methods built set.")

    reprog = sheets.get("Reprog_cell_state_signatures") or {}
    needed_states = ("Fibroblast", "PartialReprog", "EarlyPluripotency", "Pluripotency", "NonReprog")
    missing_st = [s for s in needed_states if s not in reprog]
    if missing_st:
        raise StopStep(
            "genesets",
            f"mmc3 Reprog_cell_state_signatures missing {missing_st}. cols={list(reprog)}",
        )
    state_sets = {s: [g.strip().upper() for g in reprog[s] if g.strip()] for s in needed_states}
    for s, genes in state_sets.items():
        log(f"[reprog] {s} n={len(genes)}")

    aging = sheets.get("Aging_signatures") or {}
    age_up_mmc3 = [g.strip().upper() for g in (aging.get("Age up") or []) if g.strip()]
    age_down_mmc3 = [g.strip().upper() for g in (aging.get("Age down") or []) if g.strip()]
    log(f"[mmc3-aging] Age up n={len(age_up_mmc3)} Age down n={len(age_down_mmc3)} "
        "(logged; not used as a substitute DESeq2 reconstruction)")

    # Fleischer DESeq2 reconstruction: GSE113957 is FPKM.
    aging_built = dict(ok=False, reason=None, n_up=None, n_down=None)
    if not GSE113957_FPKM.exists():
        aging_built["reason"] = (
            f"GSE113957 FPKM file missing at {GSE113957_FPKM}. "
            "FINDINGS_FIBRO2.md already established FPKM with no series matrix. Not approximating."
        )
        log(f"[aging-DE] STOP this step: {aging_built['reason']}")
    else:
        # Peek only: confirm not integer counts. Do not run DESeq2 on FPKM.
        import gzip
        with gzip.open(GSE113957_FPKM, "rt", encoding="utf-8", errors="replace") as fh:
            header = fh.readline()
            n_preview = 0
            vals = []
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                for x in parts[1:12]:
                    try:
                        vals.append(float(x))
                    except ValueError:
                        pass
                n_preview += 1
                if n_preview >= 200:
                    break
        arr = np.asarray(vals, float)
        log_columns("GSE113957_fpkm_header", header.strip().split("\t")[:20], str(GSE113957_FPKM))
        is_counts = bool(len(arr) and looks_like_counts(arr.reshape(1, -1)))
        log(f"[aging-DE] GSE113957_fpkm preview n_values={len(arr)} min={np.nanmin(arr) if len(arr) else 'NA'} "
            f"max={np.nanmax(arr) if len(arr) else 'NA'} looks_like_counts={is_counts}")
        if not is_counts:
            aging_built["reason"] = (
                "GSE113957_fpkm.txt.gz values are not integer counts (cannot apply DESeq2). "
                "FINDINGS_FIBRO2.md: series matrix HEAD 404. STAR Methods require DESeq2 old-vs-young "
                "with PCA-based old/young split. Not approximating from FPKM and not substituting "
                f"mmc3 Age up (n={len(age_up_mmc3)}) / Age down (n={len(age_down_mmc3)}). "
                "Aging-signature comparison reported as not run."
            )
            log(f"[aging-DE] STOP this step: {aging_built['reason']}")
        else:
            aging_built["reason"] = (
                "GSE113957 preview looked like counts but the published matrix is FPKM; "
                "refusing to run DESeq2 without the integer count matrix and the authors' PCA split."
            )
            log(f"[aging-DE] STOP this step: {aging_built['reason']}")

    summary = dict(
        msigdb_version=msigdb_version,
        msigdb_gmt_ok=bool(gmt.get("ok")),
        msigdb_gmt_url=gmt.get("url"),
        msigdb_gmt_errors=gmt.get("errors"),
        n_EMT=len(emt_used),
        n_TGFB=len(tgfb_used),
        n_MD_built=len(md_built),
        MD_TFs_added=added,
        MD_TFs_already_in_EMT=already_in,
        mmc3_n_MD=len(mmc3_md),
        mmc3_n_TGFB=len(mmc3_tgfb),
        mmc3_MD_equal=bool(cmp_md["equal"]),
        mmc3_TGFB_equal=bool(cmp_tgfb["equal"]),
        mmc3_MD_only_built=cmp_md["only_a"],
        mmc3_MD_only_supplement=cmp_md["only_b"],
        mmc3_TGFB_only_built=cmp_tgfb["only_a"],
        mmc3_TGFB_only_supplement=cmp_tgfb["only_b"],
        gmt_vs_local=gmt_cmp,
        aging_signature_built=False,
        aging_signature_reason=aging_built["reason"],
        mmc3_age_up_n=len(age_up_mmc3),
        mmc3_age_down_n=len(age_down_mmc3),
        mmc3_age_used_as_substitute=False,
        reprog_n={k: len(v) for k, v in state_sets.items()},
        hallmark_local_path=str(HALLMARK_LOCAL),
        mmc3_path=str(MMC3_XLSX),
        columns_read={
            "MD_signatures": list((sheets.get("MD_signatures") or {})),
            "Aging_signatures": list((sheets.get("Aging_signatures") or {})),
            "Reprog_cell_state_signatures": list((sheets.get("Reprog_cell_state_signatures") or {})),
            "Fibroblast_subtype_signatures": list((sheets.get("Fibroblast_subtype_signatures") or {})),
        },
    )
    dump_json(MD2_DIR / "genesets_summary.json", jsonable(summary))
    dump_json(MD2_DIR / "genesets.json", jsonable(dict(
        MD=md_built, TGFB=tgfb_used, EMT=emt_used,
        reprog=state_sets,
        mmc3_MD=mmc3_md, mmc3_TGFB=mmc3_tgfb,
        mmc3_age_up=age_up_mmc3, mmc3_age_down=age_down_mmc3,
        msigdb_version=msigdb_version,
    )))
    man = load_manifest()
    man["status"] = "GENESETS_DONE"
    man["aging_signature_built"] = False
    man["aging_signature_reason"] = aging_built["reason"]
    save_manifest(man)
    progress_snapshot(
        "write clustering settings flag, then cluster k=3/8/15/Louvain (no age scores yet)",
        stop="GENESETS_DONE",
    )
    if close_log:
        log.close()
    return summary


if __name__ == "__main__":
    try:
        build_gene_sets()
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise

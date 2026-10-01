"""T-A — external validation of the frozen fibroblast ruler. Refit nothing."""
from __future__ import annotations

import gzip
import io
import re
import sys
import tarfile
import time
from pathlib import Path

import numpy as np
if not hasattr(np, "unicode_"):
    np.unicode_ = np.str_  # type: ignore[attr-defined]
import pandas as pd
import requests
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fibro2_common import (  # noqa: E402
    FIBRO2_DIR, FIBRO2_RAW, FIBRO2_PROC, FIBRO2_SEED, FIBRO2_BOOT, N_PERM, N_BOOT,
    ADULT_MIN_AGE, GSE_PRIMARY, GSE_NEAR_BULK, GSE_NEAR_N9,
    CXG_FIBRO_ATLAS_ID, PREREG_TA, PREREG_TA_FLAG, FROZEN_RULER,
    StopStep, Logger, dump_json, jsonable, fibro2_log_banner,
    load_manifest, save_manifest, load_json, record_failure, record_unusable, record_unverified,
    log_columns, load_frozen_ruler, geo_nnn, strict_age_years,
    align_counts_to_ruler, tmm_logcpm_rows, score_frozen, rho_null_donor,
    looks_like_counts, progress_snapshot, pass_rho_bar,
)
from gtex_common import spearman_safe  # noqa: E402
from geo_meta import fetch_soft, parse_series  # noqa: E402
from download import download as http_download  # noqa: E402
from search_cxg import API, parse_age_years  # noqa: E402
from tissue_peek_obs import _open_h5, obs_column_names, read_obs_column  # noqa: E402
from pseudobulk import find_counts_path, read_counts_chunk  # noqa: E402

CXG_LISTING = f"{API}/datasets"


def _parse_gsm_blocks(text: str, log) -> pd.DataFrame:
    recs, cur = [], {}
    for line in text.splitlines():
        if line.startswith("^SAMPLE"):
            if cur:
                recs.append(cur)
            cur = dict(gsm=line.split("=")[-1].strip())
            continue
        if not line.startswith("!Sample_"):
            continue
        k, _, v = line[1:].partition(" = ")
        k, v = k.strip(), v.strip()
        cur.setdefault(k, [])
        cur[k].append(v)
    if cur:
        recs.append(cur)
    if not recs:
        raise StopStep("geo", "no GSM sample blocks in SOFT")
    rows = []
    fields_union = set()
    for rec in recs:
        fields_union.update(rec.keys())
        chars = rec.get("Sample_characteristics_ch1", [])
        char_map = {}
        for c in chars:
            if ":" in c:
                a, b = c.split(":", 1)
                char_map[a.strip().lower()] = b.strip()
            else:
                char_map.setdefault("_unparsed", []).append(c)
        title = " ".join(rec.get("Sample_title", []))
        desc = " ".join(rec.get("Sample_description", []))
        source = " ".join(rec.get("Sample_source_name_ch1", []))
        growth = " ".join(rec.get("Sample_growth_protocol", []))
        age_from_char = None
        age_char_key = None
        for key in ("age", "age (years)", "donor age", "donor age (years)",
                    "age_years", "age (yrs)", "age(years)"):
            if key in char_map:
                age_from_char = strict_age_years(char_map[key])
                if age_from_char is None:
                    age_from_char = strict_age_years(f"age: {char_map[key]}")
                age_char_key = key
                break
        age, age_src = None, None
        for val, src in (
            (age_from_char, f"Sample_characteristics_ch1:{age_char_key}"),
            (strict_age_years(title), "Sample_title"),
            (strict_age_years(source), "Sample_source_name_ch1"),
            (strict_age_years(desc), "Sample_description"),
            (strict_age_years(growth), "Sample_growth_protocol"),
        ):
            if val is not None:
                age, age_src = val, src
                break
        donor = (
            char_map.get("donor") or char_map.get("donor id") or char_map.get("donor_id")
            or char_map.get("individual") or char_map.get("subject") or char_map.get("cell line")
            or char_map.get("cell_line") or rec.get("gsm")
        )
        rows.append(dict(
            gsm=rec.get("gsm"), title=title, source_name=source,
            organism="; ".join(rec.get("Sample_organism_ch1", [])),
            molecule="; ".join(rec.get("Sample_molecule_ch1", [])),
            description=desc, characteristics="; ".join(chars),
            platform="; ".join(rec.get("Sample_platform_id", [])),
            supplementary="; ".join(rec.get("Sample_supplementary_file", [])),
            donor=str(donor) if donor is not None else rec.get("gsm"),
            age_years=age, age_source=age_src,
            char_keys=";".join(sorted(char_map.keys())),
        ))
    df = pd.DataFrame(rows)
    log(f"[geo] parsed {len(df)} GSM. columns actually read (union of Sample_* keys): "
        f"{sorted(fields_union)}")
    log_columns("geo_gsm", sorted(fields_union), "SOFT GSM")
    return df


def fetch_gse_record(acc: str, log):
    log(f"[geo] fetching SOFT for {acc}")
    time.sleep(0.5)
    text = fetch_soft(acc, targ_gsm=False)
    if "!Series_title" not in text:
        raise StopStep("geo", f"{acc} SOFT has no Series_title. head={text[:400]!r}")
    series = parse_series(text)
    log(f"[geo] {acc} Series_title: {series.get('Series_title', [''])[0]}")
    log(f"[geo] {acc} Series_type: {series.get('Series_type', [])}")
    log(f"[geo] {acc} Series_platform_id: {series.get('Series_platform_id', [])}")
    log(f"[geo] {acc} n Sample_id: {len(series.get('Series_sample_id', []))}")
    log(f"[geo] {acc} Series_supplementary_file: {series.get('Series_supplementary_file', [])}")
    log(f"[geo] {acc} series fields actually read: {sorted(series.keys())}")
    log_columns(f"geo_series_{acc}", sorted(series.keys()), f"{acc} SOFT series")
    time.sleep(0.5)
    gsm_text = fetch_soft(acc, targ_gsm=True)
    samples = _parse_gsm_blocks(gsm_text, log)
    rec = dict(
        accession=acc,
        title=(series.get("Series_title") or [""])[0],
        types=list(series.get("Series_type") or []),
        platforms=list(series.get("Series_platform_id") or []),
        n_sample_id=len(series.get("Series_sample_id") or []),
        supplementary=list(series.get("Series_supplementary_file") or []),
        overall_design=" ".join(series.get("Series_overall_design") or []),
        summary=" ".join(series.get("Series_summary") or [])[:2000],
        series_fields=sorted(series.keys()),
        n_gsm=int(len(samples)),
        n_with_age=int(samples.age_years.notna().sum()),
        n_missing_age=int(samples.age_years.isna().sum()),
        organisms=sorted(set(samples.organism.dropna().astype(str))) if len(samples) else [],
        molecules=sorted(set(samples.molecule.dropna().astype(str))) if len(samples) else [],
    )
    samples.to_csv(FIBRO2_DIR / f"ta_{acc}_samples.csv", index=False)
    dump_json(FIBRO2_DIR / f"ta_{acc}_series.json", jsonable(rec))
    (FIBRO2_DIR / f"ta_{acc}_columns.txt").write_text(
        "series_fields:\n" + "\n".join(rec["series_fields"]) + "\n\nsample_columns:\n"
        + "\n".join(map(str, samples.columns)) + "\n",
        encoding="utf-8",
    )
    return series, samples, rec


def classify_platform(series_rec, samples) -> dict:
    types = [t.lower() for t in series_rec.get("types") or []]
    platforms = [p.upper() for p in series_rec.get("platforms") or []]
    is_array = any("array" in t for t in types)
    is_seq = any("high throughput sequencing" in t or "rna-seq" in t or "rnaseq" in t for t in types)
    is_sc = any("single cell" in t or "single-cell" in t for t in types)
    if not is_sc:
        titles = " ".join(samples.title.astype(str).head(20).tolist()).lower() if len(samples) else ""
        design = (series_rec.get("overall_design") or "").lower()
        is_sc = ("10x" in design or "single-cell" in design or "scrna" in design
                 or "single cell" in design or "10x" in titles)
    kind = "unknown"
    if is_array and not is_seq:
        kind = "microarray"
    elif is_sc:
        kind = "single_cell"
    elif is_seq:
        kind = "bulk_rnaseq"
    return dict(kind=kind, is_array=bool(is_array), is_seq=bool(is_seq), is_sc=bool(is_sc),
                types=series_rec.get("types"), platforms=platforms,
                molecules=series_rec.get("molecules"))


def resolve_gse113957(log):
    series, samples, rec = fetch_gse_record(GSE_PRIMARY, log)
    plat = classify_platform(rec, samples)
    n_age = int(samples.age_years.notna().sum())
    n_miss = int(samples.age_years.isna().sum())
    n_child = int(((samples.age_years.notna()) & (samples.age_years < ADULT_MIN_AGE)).sum())
    n_adult = int(((samples.age_years.notna()) & (samples.age_years >= ADULT_MIN_AGE)).sum())
    log(f"[GSE113957] n_gsm={len(samples)} n_with_age={n_age} missing_age={n_miss} "
        f"age<{ADULT_MIN_AGE}={n_child} adult={n_adult} platform={plat}")
    rec_out = dict(rec)
    rec_out.update(plat)
    rec_out.update(
        role="primary_external", n_with_age=n_age, n_missing_age=n_miss,
        n_child=n_child, n_adult=n_adult, adult_min=ADULT_MIN_AGE,
    )
    if n_age == 0:
        reason = "no per-sample donor age stated in the record"
        record_unusable(GSE_PRIMARY, reason, rec_out)
        rec_out.update(usable=False, unusable_reason=reason)
        return rec_out, samples
    if plat["kind"] == "microarray" or plat["is_array"]:
        reason = (
            f"Series_type={rec.get('types')} platform={plat.get('platforms')}: "
            "microarray, not counts. Frozen GTEx spec is TMM/edgeR log2-CPM on counts. "
            "Not substituting RMA/MAS5 intensities."
        )
        record_unusable(GSE_PRIMARY, reason, rec_out)
        rec_out.update(usable=False, unusable_reason=reason)
        return rec_out, samples
    if plat["kind"] not in {"bulk_rnaseq", "single_cell"} and not plat["is_seq"]:
        reason = f"cannot apply TMM/log2-CPM: Series_type={rec.get('types')} kind={plat['kind']}"
        record_unusable(GSE_PRIMARY, reason, rec_out)
        rec_out.update(usable=False, unusable_reason=reason)
        return rec_out, samples
    rec_out.update(usable=True, unusable_reason=None)
    return rec_out, samples


def resolve_gse_generic(acc, log, role, sc_min_donors=None):
    series, samples, rec = fetch_gse_record(acc, log)
    plat = classify_platform(rec, samples)
    n_age = int(samples.age_years.notna().sum())
    n_miss = int(samples.age_years.isna().sum())
    n_donors = int(samples.loc[samples.age_years.notna(), "donor"].nunique()) if n_age else 0
    rec_out = dict(rec)
    rec_out.update(plat)
    rec_out.update(role=role, n_with_age=n_age, n_missing_age=n_miss, n_donors_with_age=n_donors)
    if n_age == 0:
        reason = "no per-sample donor age stated in the record"
        record_unusable(acc, reason, rec_out)
        rec_out.update(usable=False, unusable_reason=reason)
        return rec_out, samples
    if plat["kind"] == "microarray" or plat["is_array"]:
        reason = (
            f"Series_type={rec.get('types')} platform={plat.get('platforms')}: "
            "microarray, not counts. Frozen spec is TMM/log2-CPM."
        )
        record_unusable(acc, reason, rec_out)
        rec_out.update(usable=False, unusable_reason=reason)
        return rec_out, samples
    if sc_min_donors is not None and plat.get("is_sc") and n_donors < int(sc_min_donors):
        reason = f"single-cell n_donors_with_age={n_donors} < {sc_min_donors}"
        record_unusable(acc, reason, rec_out)
        rec_out.update(usable=False, unusable_reason=reason)
        return rec_out, samples
    if plat["kind"] not in {"bulk_rnaseq", "single_cell"} and not plat["is_seq"]:
        reason = f"cannot apply TMM/log2-CPM: Series_type={rec.get('types')} kind={plat['kind']}"
        record_unusable(acc, reason, rec_out)
        rec_out.update(usable=False, unusable_reason=reason)
        return rec_out, samples
    rec_out.update(usable=True, unusable_reason=None)
    return rec_out, samples


def search_cxg_fibroblast(log):
    log(f"[cxg] GET {CXG_LISTING}")
    r = requests.get(CXG_LISTING, timeout=600)
    r.raise_for_status()
    ds = r.json()
    log(f"[cxg] {len(ds)} datasets returned")
    rows = []
    for d in ds:
        orgs = [o.get("label") for o in d.get("organism", [])]
        if "Homo sapiens" not in orgs:
            continue
        title = str(d.get("title") or "")
        tissues = sorted({t.get("label") for t in d.get("tissue", [])})
        cts = sorted({c.get("label") for c in d.get("cell_type", [])})
        blob = " ".join([title, " ".join(tissues), " ".join(cts)]).lower()
        if not re.search(r"fibroblast|dermal|skin of", blob):
            continue
        assays = sorted({a.get("label") for a in d.get("assay", [])})
        stages = [s.get("label") for s in d.get("development_stage", [])]
        donors = d.get("donor_id", []) or []
        ages = sorted(a for a in (parse_age_years(s) for s in stages) if a is not None)
        h5ad_url = next((a.get("url") for a in d.get("assets", []) if a.get("filetype") == "H5AD"), None)
        h5ad_bytes = next((a.get("filesize") for a in d.get("assets", []) if a.get("filetype") == "H5AD"), None)
        rows.append(dict(
            dataset_id=d.get("dataset_id"), collection_id=d.get("collection_id"),
            title=title, cell_count=d.get("cell_count"), n_donors=len(donors),
            n_cell_types=len(cts), assays="; ".join(assays),
            tissues="; ".join(tissues)[:400], cell_types="; ".join(cts)[:400],
            n_dev_stages=len(stages), n_stages_numeric_years=len(ages),
            age_min=ages[0] if ages else None, age_max=ages[-1] if ages else None,
            h5ad_url=h5ad_url, h5ad_bytes=h5ad_bytes,
            is_10x=bool(re.search(r"10x", "; ".join(assays), re.I)),
            fibroblast_in_title=bool(re.search(r"fibroblast", title, re.I)),
            dermal_in_title=bool(re.search(r"dermal fibroblast", title, re.I)),
        ))
    df = pd.DataFrame(rows)
    df.to_csv(FIBRO2_DIR / "ta_cxg_fibroblast_hits.csv", index=False)
    log(f"[cxg] fibroblast/skin title-or-tissue hits: {len(df)}")
    cand = df.copy()
    if len(cand):
        cand = cand[(cand.n_donors >= 10) & (cand.n_stages_numeric_years >= 8)]
    log(f"[cxg] ≥10 donors and ≥8 numeric development_stage years: {len(cand)}")
    if len(cand):
        cand = cand.sort_values(
            ["dermal_in_title", "fibroblast_in_title", "is_10x", "n_donors"],
            ascending=[False, False, False, False],
        )
        cand.to_csv(FIBRO2_DIR / "ta_cxg_candidates.csv", index=False)
        log("[cxg] top candidates:\n" + cand.head(15).to_string(index=False)[:4000])
    else:
        cand = pd.DataFrame()
    dump_json(FIBRO2_DIR / "ta_cxg_search.json", jsonable(dict(
        n_listing=len(ds), n_fibro_hits=int(len(df)), n_candidates=int(len(cand)),
        atlas_id=CXG_FIBRO_ATLAS_ID,
        atlas_in_listing=bool(len(df) and (df.dataset_id.astype(str) == CXG_FIBRO_ATLAS_ID).any()),
    )))
    return df, cand


def peek_cxg_obs(url, log, name):
    log(f"[cxg] peek obs {name} {url}")
    h, handle = _open_h5(url)
    try:
        cols = obs_column_names(h)
        log(f"[cxg] obs columns actually read ({len(cols)}): {cols}")
        log_columns(f"cxg_obs_{name}", cols, url)
        shape = None
        for p in ("raw/X", "X"):
            if p in h:
                sh = h[p].attrs.get("shape")
                if sh is not None:
                    shape = tuple(int(x) for x in sh)
                    log(f"[cxg] matrix {p} shape={shape}")
        try:
            cpath = find_counts_path(h)
            counts_note = cpath
            log(f"[cxg] integer counts path={cpath}")
        except Exception as e:
            counts_note = f"NO_INTEGER_COUNTS: {e}"
            log(f"[cxg] {counts_note}")
        summary = dict(obs_columns=cols, shape=shape, counts_path=counts_note, source=url)
        wanted = [c for c in (
            "donor_id", "development_stage", "sex", "cell_type", "assay", "tissue",
            "organism", "disease", "suspension_type", "is_primary_data", "age",
        ) if c in cols]
        for c in wanted:
            arr = read_obs_column(h, c)
            vc = pd.Series(arr).value_counts(dropna=False)
            log(f"[cxg] {c}: n={len(arr):,} nunique={vc.shape[0]} top={list(vc.head(12).items())}")
            summary[f"nunique_{c}"] = int(vc.shape[0])
        dump_json(FIBRO2_DIR / f"ta_cxg_{name}_peek.json", jsonable(summary))
        return summary, cols
    finally:
        h.close()
        if handle is not None:
            handle.close()


def _http_url(url: str) -> str:
    u = str(url).strip()
    if u.startswith("ftp://"):
        u = "https://" + u[len("ftp://"):]
    return u


def download_suppl_counts(acc, series_rec, log):
    urls = list(series_rec.get("supplementary") or [])
    dest = FIBRO2_RAW / acc
    dest.mkdir(parents=True, exist_ok=True)
    if not urls:
        record_unusable(acc, "no Series_supplementary_file in SOFT", dict(accession=acc))
        return None
    got, count_like, other = [], [], []
    for url in urls:
        url = _http_url(url)
        if not url:
            continue
        low = url.lower()
        if any(x in low for x in (".cel", ".cel.gz", ".idat", ".chp", ".txt.md5", ".tar.md5")):
            log(f"[download] skip (not expression table): {url}")
            continue
        log(f"[download] {url}")
        try:
            path = Path(http_download(url, str(dest)))
        except Exception as e:
            log(f"[download] FAIL {url}: {type(e).__name__}: {e}")
            record_failure("download", f"{url}: {type(e).__name__}: {e}")
            continue
        rec = dict(url=url, path=str(path), size=int(path.stat().st_size) if path.exists() else 0)
        got.append(rec)
        if any(k in low for k in ("count", "raw", "gene_reads", "counts")):
            count_like.append(rec)
        else:
            other.append(rec)
        log(f"[download] {path.name} size={rec['size']}")
    dump_json(FIBRO2_DIR / f"ta_{acc}_downloads.json", jsonable(dict(got=got, count_like=count_like, other=other)))
    ordered = count_like + other
    nnn = geo_nnn(acc)
    sm_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{nnn}/{acc}/matrix/{acc}_series_matrix.txt.gz"
    log(f"[download] series matrix HEAD {sm_url}")
    try:
        hr = requests.head(sm_url, allow_redirects=True, timeout=30)
        if hr.status_code == 404:
            log(f"[download] series matrix 404 — skip (not retrying)")
        else:
            path = Path(http_download(sm_url, str(dest), retries=2))
            ordered.append(dict(url=sm_url, path=str(path), size=int(path.stat().st_size), kind="series_matrix"))
    except Exception as e:
        log(f"[download] series matrix FAIL: {type(e).__name__}: {e}")
    if not ordered:
        record_unusable(acc, "no supplementary or series-matrix file downloaded", dict(urls=urls))
        return None
    return ordered


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def load_series_matrix_counts(path: Path, log, gsm_keep=None):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()
    begin = end = None
    for i, line in enumerate(lines):
        if line.startswith("!series_matrix_table_begin"):
            begin = i + 1
        if line.startswith("!series_matrix_table_end"):
            end = i
    if begin is None or end is None:
        raise StopStep("series_matrix", f"{path.name}: no series_matrix_table_begin/end")
    df = pd.read_csv(io.StringIO("".join(lines[begin:end])), sep="\t", index_col=0)
    log(f"[series_matrix] genes={df.shape[0]} samples={df.shape[1]} "
        f"columns actually read (first 12)={list(df.columns)[:12]}")
    log_columns(f"series_matrix_{path.stem}", list(df.columns), str(path))
    if gsm_keep is not None:
        keep = [c for c in df.columns if str(c) in set(gsm_keep)]
        if not keep:
            raise StopStep("series_matrix", f"{path.name}: none of requested GSM columns present")
        df = df[keep]
    mat = df.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=np.float64)
    if not looks_like_counts(mat):
        raise StopStep(
            "series_matrix",
            f"{path.name}: table is not integer counts (likely processed intensities). "
            f"min={np.nanmin(mat)} max={np.nanmax(mat)}. Not substituting.",
        )
    return dict(counts=mat, genes=df.index.astype(str).to_numpy(),
                symbols=df.index.astype(str).to_numpy(),
                samples=[str(c) for c in df.columns], path=str(path))


def load_gene_count_tar(path: Path, log, gsm_keep=None):
    """Merge per-GSM *geneCOUNT*.txt.gz members. Sample id = GSM prefix of the member name."""
    path = Path(path)
    keep = set(map(str, gsm_keep)) if gsm_keep is not None else None
    with tarfile.open(path) as tar:
        members = [m for m in tar.getmembers() if m.isfile()]
        names = [Path(m.name).name for m in members]
        gc = [m for m in members if "genecount" in Path(m.name).name.lower()]
        log(f"[counts] {path.name} tar files={len(names)} geneCOUNT={len(gc)}")
        if not gc:
            raise StopStep(
                "counts",
                f"{path.name}: tar has no *geneCOUNT* members. first names={names[:8]}",
            )
        genes0 = None
        cols = []
        samples = []
        for m in sorted(gc, key=lambda x: Path(x.name).name):
            bn = Path(m.name).name
            gsm = bn.split("_", 1)[0]
            if keep is not None and gsm not in keep:
                continue
            fh = tar.extractfile(m)
            if fh is None:
                raise StopStep("counts", f"{path.name}: cannot extract {m.name}")
            raw = gzip.GzipFile(fileobj=fh)
            df = pd.read_csv(raw, sep="\t", index_col=0)
            num = df.select_dtypes(include=[np.number])
            if num.shape[1] == 0:
                num = df.apply(pd.to_numeric, errors="coerce")
                num = num.dropna(axis=1, how="all")
            if num.shape[1] != 1:
                raise StopStep(
                    "counts",
                    f"{bn}: expected one numeric count column, columns actually read={list(df.columns)}",
                )
            genes = df.index.astype(str)
            vec = num.iloc[:, 0].to_numpy(dtype=np.float64)
            if genes0 is None:
                genes0 = genes
            elif not genes.equals(genes0):
                if set(genes) != set(genes0) or len(genes) != len(genes0):
                    raise StopStep(
                        "counts",
                        f"{bn}: gene index differs from first geneCOUNT "
                        f"(n={len(genes)} vs {len(genes0)}). Not aligning by name silently "
                        "beyond a same-set reindex.",
                    )
                vec = pd.Series(vec, index=genes).reindex(genes0).to_numpy(dtype=np.float64)
            cols.append(vec)
            samples.append(gsm)
            if len(samples) == 1:
                log_columns("geneCOUNT_first", list(df.columns), f"{path.name}:{bn}")
        if not cols:
            raise StopStep("counts", f"{path.name}: no geneCOUNT members matched gsm_keep")
        mat = np.column_stack(cols)
        log(f"[counts] {path.name} merged geneCOUNT genes={len(genes0)} samples={len(samples)} "
            f"columns actually read (GSM, first 12)={samples[:12]}")
        if not looks_like_counts(mat):
            raise StopStep(
                "counts",
                f"{path.name} geneCOUNT merge: values are not integer counts. "
                f"min={np.nanmin(mat)} max={np.nanmax(mat)} mean={np.nanmean(mat)}",
            )
        return dict(
            counts=mat, genes=genes0.to_numpy(), symbols=genes0.to_numpy(),
            samples=samples, path=str(path),
        )


def load_count_table(path: Path, log, gsm_keep=None):
    path = Path(path)
    log(f"[counts] reading {path.name}")
    if "series_matrix" in path.name.lower():
        return load_series_matrix_counts(path, log, gsm_keep=gsm_keep)
    low = path.name.lower()
    if low.endswith(".tar") or low.endswith(".tar.gz"):
        return load_gene_count_tar(path, log, gsm_keep=gsm_keep)
    with _open_text(path) as fh:
        head = fh.read(8192)
    sep = "\t" if head.count("\t") >= head.count(",") else ","
    try:
        df = pd.read_csv(path, sep=sep, index_col=0, comment="#", low_memory=False)
    except ValueError as e:
        raise StopStep("counts", f"{path.name}: {e}") from e
    if gsm_keep is not None:
        gsm_cols = [c for c in df.columns if str(c) in set(gsm_keep)]
        if gsm_cols:
            sub = df[gsm_cols]
        else:
            num = df.select_dtypes(include=[np.number])
            if num.shape[1] == 0:
                raise StopStep("counts", f"{path.name}: no numeric sample columns matching GSM")
            sub = num
    else:
        sub = df.select_dtypes(include=[np.number])
        if sub.shape[1] == 0:
            raise StopStep("counts", f"{path.name}: no numeric sample columns")
    mat = sub.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=np.float64)
    samples = [str(c) for c in sub.columns]
    genes = sub.index.astype(str).to_numpy()
    log(f"[counts] {path.name} genes={len(genes)} samples={len(samples)} "
        f"columns actually read (first 12)={samples[:12]}")
    if not looks_like_counts(mat):
        raise StopStep(
            "counts",
            f"{path.name}: values are not integer counts (cannot apply TMM/log2-CPM). "
            f"min={np.nanmin(mat)} max={np.nanmax(mat)} mean={np.nanmean(mat)}",
        )
    return dict(counts=mat, genes=genes, symbols=genes, samples=samples, path=str(path))


def project_cohort(tag, counts, genes, symbols, samples, obs, frozen, log, adult_only=False):
    obs = obs.copy()
    pos = {str(s): i for i, s in enumerate(samples)}
    missing = [s for s in obs["sample"].astype(str) if s not in pos]
    if missing:
        raise StopStep("project", f"{tag}: {len(missing)} samples not in count matrix e.g. {missing[:8]}")
    col = np.array([pos[s] for s in obs["sample"].astype(str)], dtype=int)
    C = np.asarray(counts, np.float64)[:, col]
    n_excl_age = int(obs.age_years.isna().sum())
    if n_excl_age:
        log(f"[{tag}] excluding {n_excl_age} samples with no stated age (not inferred)")
        keep = obs.age_years.notna().to_numpy()
        obs = obs.loc[keep].reset_index(drop=True)
        C = C[:, keep]
    n_child = 0
    if adult_only:
        child = (obs.age_years.astype(float) < ADULT_MIN_AGE).to_numpy()
        n_child = int(child.sum())
        log(f"[{tag}] excluding {n_child} samples with stated age < {ADULT_MIN_AGE}")
        keep = ~child
        obs = obs.loc[keep].reset_index(drop=True)
        C = C[:, keep]
    if len(obs) < 8:
        raise StopStep(
            "project",
            f"{tag}: n={len(obs)} after age filters; spearman_safe requires n≥8. Not lowering the bar.",
        )
    Y, align = align_counts_to_ruler(C, genes, symbols, frozen, log, tag=tag)
    logcpm, nf = tmm_logcpm_rows(Y, log, tag=tag)
    pred, Z, missing_g = score_frozen(logcpm, frozen, counts=Y)
    log(f"[{tag}] missing-gene columns set to z=0: {int(missing_g.sum())}")
    y = obs.age_years.to_numpy(float)
    donor = obs.donor.astype(str).to_numpy()
    sc, rho_null, p, nulls = rho_null_donor(y, pred, donor, n_perm=N_PERM, seed=FIBRO2_SEED)
    rng_b = np.random.default_rng(FIBRO2_BOOT)
    d_u = pd.unique(donor)
    n_d = int(len(d_u))
    by = {d: np.flatnonzero(donor == d) for d in d_u}
    boot_vals = []
    for _ in range(N_BOOT):
        take_d = rng_b.choice(d_u, size=n_d, replace=True)
        yy, pp = [], []
        for d in take_d:
            rows = by[d]
            yy.extend(y[rows].tolist())
            pp.extend(np.asarray(pred)[rows].tolist())
        r = spearman_safe(np.asarray(yy, float), np.asarray(pp, float))
        if np.isfinite(r):
            boot_vals.append(float(r))
    if len(boot_vals) >= 2:
        ci_lo = float(np.percentile(boot_vals, 2.5))
        ci_hi = float(np.percentile(boot_vals, 97.5))
    else:
        ci_lo = ci_hi = np.nan
    passed = pass_rho_bar(sc["rho"], rho_null)
    out = dict(
        tag=tag, n_samples=int(len(obs)), n_donors=int(pd.Series(donor).nunique()),
        n_excluded_no_age=n_excl_age, n_excluded_child=n_child,
        n_overlap=align["n_overlap"], n_missing_z0=int(missing_g.sum()),
        n_ruler=align["n_ruler"], n_ensembl=align["n_ensembl"],
        n_symbol_fallback=align["n_symbol_fallback"],
        rho=sc["rho"], rho_null=rho_null, rho_p=p, rho_ci_lo=ci_lo, rho_ci_hi=ci_hi,
        r=sc["r"], r2=sc["r2"], cal_r2=sc["cal_r2"], pass_rho=passed,
        n_perm=N_PERM, n_boot=N_BOOT, seed=FIBRO2_SEED, boot_seed=FIBRO2_BOOT,
        age_min=float(np.nanmin(y)), age_max=float(np.nanmax(y)),
        n_unique_age=int(pd.Series(y).nunique()),
    )
    obs = obs.copy()
    obs["age_score"] = pred
    obs.to_csv(FIBRO2_DIR / f"ta_{tag}_scores.csv", index=False)
    np.savez_compressed(
        FIBRO2_PROC / f"ta_{tag}_proj.npz",
        logcpm=logcpm.astype(np.float32), pred=np.asarray(pred, np.float64),
        y=y, donor=donor.astype(object), missing=missing_g, Z=Z.astype(np.float32),
    )
    dump_json(FIBRO2_DIR / f"ta_{tag}_align.json", jsonable(align))
    dump_json(FIBRO2_DIR / f"ta_{tag}_result.json", jsonable(out))
    log(f"[{tag}] ρ={out['rho']:+.3f} null={rho_null:+.3f} p={p:+.3f} "
        f"CI=[{ci_lo:+.3f}, {ci_hi:+.3f}] calR²={sc['cal_r2']:+.3f} R²={sc['r2']:+.3f} "
        f"pass={passed} n={out['n_samples']} n_donors={out['n_donors']} overlap={out['n_overlap']}")
    return out, obs, Y, logcpm


def build_obs_from_samples(samples, gsm_in_matrix):
    sub = samples[samples.gsm.astype(str).isin(set(map(str, gsm_in_matrix)))].copy()
    sub["sample"] = sub.gsm.astype(str)
    if "donor" not in sub.columns:
        sub["donor"] = sub["sample"]
    return sub[["sample", "donor", "age_years", "age_source", "title", "gsm"]]


def try_project_gse(acc, rec, samples, frozen, log, adult_only=False, role=""):
    if not rec.get("usable"):
        log(f"[{acc}] not usable: {rec.get('unusable_reason')}")
        return None
    files = download_suppl_counts(acc, rec, log)
    if not files:
        return None
    last_err = None
    gsm_keep = samples.gsm.astype(str).tolist()
    for f in files:
        path = Path(f["path"])
        try:
            pack = load_count_table(path, log, gsm_keep=gsm_keep)
        except StopStep as e:
            last_err = e
            log(f"[{acc}] {path.name} STOP [{e.step}] {e.message}")
            continue
        except ValueError as e:
            last_err = StopStep("counts", f"{path.name}: {e}")
            log(f"[{acc}] {path.name} STOP [counts] {e}")
            continue
        obs = build_obs_from_samples(samples, pack["samples"])
        tag = f"{acc}_{role or rec.get('role') or 'cohort'}".replace(" ", "_")
        try:
            out, obs_sc, Y, logcpm = project_cohort(
                tag, pack["counts"], pack["genes"], pack["symbols"], pack["samples"],
                obs, frozen, log, adult_only=adult_only,
            )
            out["accession"] = acc
            out["role"] = rec.get("role")
            out["platform_kind"] = rec.get("kind")
            out["types"] = rec.get("types")
            out["count_file"] = str(path)
            return dict(result=out, obs=obs_sc, Y=Y, logcpm=logcpm, rec=rec)
        except StopStep as e:
            last_err = e
            log(f"[{acc}] project STOP [{e.step}] {e.message}")
            continue
    reason = f"no usable count matrix for {acc}"
    if last_err is not None:
        reason = f"{reason}: last error [{last_err.step}] {last_err.message}"
    record_unusable(acc, reason, dict(files=files))
    return None


def _decode_arr(arr):
    return np.array([
        x.decode("utf-8", "replace") if isinstance(x, (bytes, np.bytes_)) else str(x)
        for x in arr
    ], dtype=object)


def project_cxg_fibroblast(dataset_id, url, frozen, log, name="cxg_fibro"):
    log(f"[cxg] opening {url}")
    h, handle = _open_h5(str(url) if Path(str(url)).exists() else url)
    t0 = time.time()
    try:
        cols = obs_column_names(h)
        log(f"[cxg] obs columns actually read: {cols}")
        log_columns(f"cxg_{name}", cols, url)
        if "donor_id" not in cols:
            raise StopStep("cxg", f"obs missing donor_id. columns={cols}")
        age_col = None
        for c in ("development_stage", "age", "donor_age", "Age", "age_years"):
            if c in cols:
                age_col = c
                break
        if age_col is None:
            raise StopStep("cxg", f"no age-like obs column. columns={cols}")
        org_col = "organism" if "organism" in cols else None
        assay_col = "assay" if "assay" in cols else None
        donor = np.asarray(read_obs_column(h, "donor_id")).astype(str)
        age_raw = np.asarray(read_obs_column(h, age_col))
        n = len(donor)
        log(f"[cxg] n_cells={n:,} age_col={age_col}")
        ages = np.array([strict_age_years(x) for x in age_raw], dtype=object)
        stated = np.array([a is not None for a in ages])
        keep = stated.copy()
        if org_col:
            org = np.asarray(read_obs_column(h, org_col)).astype(str)
            human = np.array([
                "homo sapiens" in x.lower() or x.lower() in {"human", "ncbi9606"}
                or "NCBITaxon:9606" in x for x in org
            ])
            log(f"[cxg] organism col={org_col} human={int(human.sum()):,}/{n:,}")
            keep &= human
        assay = np.asarray(read_obs_column(h, assay_col)).astype(str) if assay_col else np.array(["unknown"] * n)
        is_10x = np.array(["10x" in x.lower() for x in assay])
        log(f"[cxg] 10x cells={int(is_10x.sum()):,}/{n:,}")
        keep_10x = keep & is_10x
        n_d_10x = int(pd.Series(donor[keep_10x]).nunique()) if keep_10x.any() else 0
        used_10x = n_d_10x >= 10
        if used_10x:
            keep = keep_10x
            platform = "10x (restricted)"
        else:
            platform = f"mixed/non-10x; 10x donors with stated age={n_d_10x}"
            log(f"[cxg] 10x-only donors with stated age={n_d_10x} < 10; not calling this a 10x validation")
        n_d = int(pd.Series(donor[keep]).nunique()) if keep.any() else 0
        if n_d < 10:
            raise StopStep(
                "cxg",
                f"{name}: n_donors with stated age (and organism filter)={n_d} < 10. Not lowering the bar.",
                dict(n_donors=n_d, age_col=age_col, platform=platform),
            )
        age_num = np.array([float(a) if a is not None else np.nan for a in ages])
        adult = keep & (age_num >= ADULT_MIN_AGE)
        n_d_ad = int(pd.Series(donor[adult]).nunique()) if adult.any() else 0
        n_excl_child_donors = n_d - n_d_ad
        if n_d_ad >= 10:
            use = adult
            adult_applied = True
        else:
            log(f"[cxg] adult≥{ADULT_MIN_AGE} donors={n_d_ad} < 10; scoring stated-age donors including <18")
            use = keep
            adult_applied = False
        try:
            cpath = find_counts_path(h)
        except Exception as e:
            raise StopStep("cxg", f"no integer count matrix: {e}")
        n_var = int(h[cpath].attrs["shape"][1])
        var_group = "raw/var" if cpath.startswith("raw/") else "var"
        if var_group not in h:
            raise StopStep("cxg", f"missing {var_group}")
        vg = h[var_group]
        if "_index" not in vg:
            raise StopStep("cxg", f"{var_group} has no _index. keys={list(vg.keys())[:20]}")
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
        log(f"[cxg] counts={cpath} n_var={n_var} gene_id n={len(gene_id)}")
        use_idx = np.flatnonzero(use)
        d_keep = donor[use]
        y_keep = age_num[use]
        tmp = pd.DataFrame({"donor": d_keep, "age": y_keep})
        bad = tmp.groupby("donor")["age"].nunique()
        bad = bad[bad > 1]
        if len(bad):
            raise StopStep("cxg", f"donor age not constant. examples={list(bad.head(10).index)}")
        d_u = pd.unique(d_keep)
        gidx = {d: i for i, d in enumerate(d_u)}
        cell_g = np.array([gidx[d] for d in d_keep], dtype=int)
        n_grp = len(d_u)
        X_sum = np.zeros((n_grp, n_var), dtype=np.float64)
        n_cells = np.zeros(n_grp, dtype=int)
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
            for g in gk:
                n_cells[g] += 1
            log(f"   rows {r1:,}/{n:,}  ({time.time()-t0:.0f}s)")
        age_d = np.array([float(tmp.loc[tmp.donor == d, "age"].iloc[0]) for d in d_u])
        assay_d = []
        for d in d_u:
            m = (donor == d) & use
            assay_d.append(pd.Series(assay[m]).value_counts().index[0])
        obs = pd.DataFrame({
            "sample": d_u, "donor": d_u, "age_years": age_d,
            "age_source": f"h5ad_obs:{age_col}", "n_cells": n_cells, "assay": assay_d,
        })
        log(f"[cxg] pseudobulk donors={n_grp} adult_applied={adult_applied} "
            f"n_excl_child_donors={n_excl_child_donors} platform={platform}")
        out, obs_sc, Y, logcpm = project_cohort(
            name, X_sum.T, gene_id, symbol if symbol is not None else gene_id,
            list(d_u), obs, frozen, log, adult_only=False,
        )
        out.update(
            accession=dataset_id, role="single_cell_external", platform_kind="single_cell",
            platform_detail=platform, used_10x_only=bool(used_10x),
            adult_applied=bool(adult_applied), n_excluded_child_donors=int(n_excl_child_donors),
            age_col=age_col, counts_path=cpath, h5ad_url=url,
        )
        return dict(result=out, obs=obs_sc, Y=Y, logcpm=logcpm)
    finally:
        h.close()
        if handle is not None:
            handle.close()


def _cohort_row(rec, samples):
    return dict(
        accession=rec.get("accession"), role=rec.get("role"), title=rec.get("title"),
        platform_kind=rec.get("kind"), types="; ".join(rec.get("types") or []),
        platforms="; ".join(rec.get("platforms") or []),
        n_gsm=rec.get("n_gsm"), n_with_age=rec.get("n_with_age"),
        n_missing_age=rec.get("n_missing_age"), n_adult=rec.get("n_adult"),
        n_child=rec.get("n_child"), n_donors_with_age=rec.get("n_donors_with_age"),
        usable=rec.get("usable"), unusable_reason=rec.get("unusable_reason"),
        molecules="; ".join(rec.get("molecules") or []),
    )


def _append_resolved(rec):
    man = load_manifest()
    row = jsonable({
        k: rec.get(k) for k in (
            "accession", "role", "title", "kind", "usable", "unusable_reason",
            "n_gsm", "n_with_age", "n_adult", "n_child", "types", "platforms",
        )
    })
    existing = man.setdefault("resolved", [])
    key = (row.get("accession"), row.get("role"))
    if not any((r.get("accession"), r.get("role")) == key for r in existing):
        existing.append(row)
    save_manifest(man)


def _ta_reading(primary, sc_row, sc_resolved, rec_p):
    if primary is None:
        man = load_manifest()
        u = next((x for x in (man.get("unusable") or []) if x.get("accession") == GSE_PRIMARY), None)
        primary_fail_reason = (
            (rec_p or {}).get("unusable_reason")
            or (u or {}).get("reason")
            or "primary external cell was not scored"
        )
        reading = (
            f"primary external cell (GSE113957) was not scored ({primary_fail_reason}). "
            "The pre-registered primary-cell pass/fail on ρ therefore cannot fire. "
        )
        if sc_resolved and sc_row and pass_rho_bar(sc_row.get("rho"), sc_row.get("rho_null")):
            reading += (
                "A single-cell external cohort passed ρ>0 with null≤0.05; see the results table. "
                "It is not the pre-registered primary cell, so the primary-cell interpretable "
                "bullet does not fire. Stage 2 `no_decline` must not be reported as a finding "
                "about reprogramming until the primary external cell is scored and passes. "
                "A 10x cohort was scored; this is not the case 'no single-cell cohort resolved'. "
                "Do not substitute the single-cell pass for GSE113957."
            )
        elif not sc_resolved:
            reading += (
                "No single-cell fibroblast cohort resolved. "
                "A bulk cohort cannot close the platform question by substitution. "
                "The Stage 2 projection rests on an unvalidated cross-platform step; "
                "do not resolve it by argument. `no_decline` must not be reported as a "
                "finding about reprogramming until a working external instrument is shown."
            )
        else:
            reading += "Single-cell external did not pass. Platform shift unresolved."
        return dict(key="primary_unscored", text=reading, interpretable=False)
    p_ok = bool(primary.get("pass_rho")) and float(primary.get("rho") or 0) > 0
    if not p_ok:
        text = (
            "ρ ≤ 0 or null > 0.05 on the primary external cell → the ruler does not "
            "validate externally; Stage 2 is uninterpretable and `no_decline` must not "
            "be reported as a finding about reprogramming."
        )
        return dict(key="primary_fail", text=text, interpretable=False)
    sc_ok = bool(sc_resolved and sc_row and sc_row.get("pass_rho")
                 and float(sc_row.get("rho") or 0) > 0)
    if p_ok and sc_ok:
        text = (
            "ρ > 0 with null ≤ 0.05 in the primary external cell → the ruler is a working "
            "instrument outside its training cohort; the Stage 2 `no_decline` verdict is "
            "interpretable as a statement about OSKM. Single-cell external also passed; "
            "the platform shift is not the remaining objection from T-A."
        )
        return dict(key="primary_pass_sc_pass", text=text, interpretable=True)
    text = (
        "bulk external passes, single-cell external fails (or none resolves) → "
        "the platform shift is unresolved; the Stage 2 projection rests on an "
        "unvalidated cross-platform step; do not resolve it by argument."
    )
    return dict(key="bulk_pass_sc_unresolved", text=text, interpretable=False, bulk_ok=True)


def run_ta(log=None):
    if not PREREG_TA_FLAG.exists():
        raise StopStep("prereg", "PREREG_TA.flag missing — write T-A pre-registration before any projection")
    close_log = False
    if log is None:
        log = Logger(FIBRO2_DIR / "ta_report.txt")
        close_log = True
    fibro2_log_banner(log, "T-A")
    log(PREREG_TA)
    frozen = load_frozen_ruler()
    log(f"[frozen] n_genes={len(frozen['w'])} n_donors={int(frozen['n_donors'])} "
        f"method={frozen['method']} regime={frozen['regime']}")
    man = load_manifest()
    man["status"] = "TA_RUNNING"
    man["frozen_ruler"] = str(FROZEN_RULER)
    save_manifest(man)
    progress_snapshot("T-A resolve GSE113957", stop="TA_RUNNING")

    cohort_rows, results, scored = [], [], []

    rec_p, samp_p = resolve_gse113957(log)
    cohort_rows.append(_cohort_row(rec_p, samp_p))
    _append_resolved(rec_p)
    progress_snapshot("T-A project GSE113957 if usable", extra=f"GSE113957 usable={rec_p.get('usable')}")
    got = try_project_gse(GSE_PRIMARY, rec_p, samp_p, frozen, log, adult_only=True, role="primary")
    if got:
        results.append(got["result"])
        scored.append(got)

    progress_snapshot("T-A search CELLxGENE fibroblast/skin")
    sc_resolved, sc_reason = False, None
    cxg_cached = FIBRO2_DIR / "ta_cxg_a19d1667_result.json"
    if cxg_cached.exists():
        out = load_json(cxg_cached)
        out.setdefault("accession", CXG_FIBRO_ATLAS_ID)
        out.setdefault("role", "single_cell_external")
        out.setdefault("platform_kind", "single_cell")
        out.setdefault("used_10x_only", True)
        out.setdefault("platform_detail", "10x (restricted)")
        dump_json(cxg_cached, jsonable(out))
        results.append(out)
        scored.append(dict(result=out))
        sc_resolved = True
        peek = {}
        peek_p = FIBRO2_DIR / "ta_cxg_a19d1667_peek.json"
        if peek_p.exists():
            peek = load_json(peek_p)
        shape0 = (peek.get("shape") or [None])[0]
        cohort_rows.append(dict(
            accession=CXG_FIBRO_ATLAS_ID, role="single_cell_external",
            title="Human dermal fibroblast atlas",
            platform_kind="single_cell",
            assays="10x 3' v1/v2/v3 + 10x 5' v1 + MARS-seq + Seq-Well S3 (scored 10x-restricted)",
            n_donors=out.get("n_donors"), usable=True, n_gsm=shape0,
            n_with_age=out.get("n_donors"), types="single-cell RNA-seq",
            h5ad_url=peek.get("source"),
        ))
        _append_resolved(dict(
            accession=CXG_FIBRO_ATLAS_ID, role="single_cell_external", usable=True,
            title="Human dermal fibroblast atlas", kind="single_cell",
        ))
        log(f"[cxg] reuse cached result ρ={out.get('rho'):+.3f} p={out.get('rho_p'):+.3f} "
            f"n_donors={out.get('n_donors')} overlap={out.get('n_overlap')}")

    if not sc_resolved:
        hits, cand = search_cxg_fibroblast(log)
        if cand is None or (hasattr(cand, "__len__") and len(cand) == 0):
            sc_reason = (
                "CELLxGENE listing: no fibroblast/skin dataset with ≥10 donors and "
                "≥8 numeric development_stage years"
            )
            record_unverified("CELLxGENE_fibroblast_search", sc_reason)
            cand = None
    else:
        cand = None
    if cand is not None and len(cand) > 0:
        order = cand.copy()
        if (order.dataset_id.astype(str) == CXG_FIBRO_ATLAS_ID).any():
            top_ids = [CXG_FIBRO_ATLAS_ID] + [
                i for i in order.dataset_id.astype(str) if i != CXG_FIBRO_ATLAS_ID
            ]
        else:
            top_ids = order.dataset_id.astype(str).tolist()
        dest = FIBRO2_RAW / "cxg"
        dest.mkdir(parents=True, exist_ok=True)
        for did in top_ids[:5]:
            row = order.loc[order.dataset_id.astype(str) == did].iloc[0]
            url = row.h5ad_url
            if not url or (isinstance(url, float) and not np.isfinite(url)):
                record_unverified(did, "listing has no h5ad_url")
                continue
            log(f"[cxg] candidate {did} title={row.title!r} n_donors={row.n_donors} "
                f"assays={row.assays} bytes={row.h5ad_bytes}")
            try:
                peek_cxg_obs(url, log, did[:8])
            except Exception as e:
                log(f"[cxg] peek FAIL {did}: {type(e).__name__}: {e}")
                record_unverified(did, f"obs peek failed: {type(e).__name__}: {e}")
                continue
            try:
                local = Path(http_download(url, str(dest)))
                got_sc = project_cxg_fibroblast(did, str(local), frozen, log, name=f"cxg_{did[:8]}")
                results.append(got_sc["result"])
                scored.append(got_sc)
                sc_resolved = True
                cohort_rows.append(dict(
                    accession=did, role="single_cell_external", title=row.title,
                    platform_kind="single_cell", assays=row.assays, n_donors=row.n_donors,
                    usable=True, n_gsm=row.cell_count,
                    n_with_age=got_sc["result"].get("n_donors"),
                    types="single-cell RNA-seq", h5ad_url=url,
                ))
                _append_resolved(dict(
                    accession=did, role="single_cell_external", usable=True,
                    title=row.title, kind="single_cell",
                ))
                break
            except StopStep as e:
                log(f"[cxg] STOP [{e.step}] {e.message}")
                record_unusable(did, f"[{e.step}] {e.message}", e.details)
                continue
            except Exception as e:
                log(f"[cxg] FAIL {type(e).__name__}: {e}")
                record_unverified(did, f"{type(e).__name__}: {e}")
                continue
        if not sc_resolved:
            sc_reason = "no CELLxGENE fibroblast candidate yielded a countable donor-age matrix with ≥10 donors"
            record_unverified("CELLxGENE_fibroblast_project", sc_reason)

    for acc, role in ((GSE_NEAR_BULK, "near_miss_bulk"), (GSE_NEAR_N9, "near_miss_n9")):
        progress_snapshot(f"T-A resolve {acc}")
        try:
            rec, samp = resolve_gse_generic(acc, log, role, sc_min_donors=10)
            cohort_rows.append(_cohort_row(rec, samp))
            _append_resolved(rec)
            if rec.get("usable") and rec.get("kind") == "single_cell" and rec.get("n_donors_with_age", 0) >= 10:
                got = try_project_gse(acc, rec, samp, frozen, log, adult_only=True, role=role)
                if got and not sc_resolved:
                    got["result"]["role"] = "single_cell_external"
                    results.append(got["result"])
                    scored.append(got)
                    sc_resolved = True
            elif rec.get("usable"):
                got = try_project_gse(acc, rec, samp, frozen, log, adult_only=True, role=role)
                if got:
                    results.append(got["result"])
                    scored.append(got)
        except StopStep as e:
            log(f"[{acc}] STOP [{e.step}] {e.message}")
            record_unusable(acc, f"[{e.step}] {e.message}", e.details)
        except Exception as e:
            log(f"[{acc}] FAIL {type(e).__name__}: {e}")
            record_unusable(acc, f"{type(e).__name__}: {e}")

    if not sc_resolved:
        log("[T-A] no single-cell fibroblast cohort resolved this session. "
            "Not substituting a bulk cohort as the platform answer.")

    pd.DataFrame(cohort_rows).to_csv(FIBRO2_DIR / "ta_cohorts.csv", index=False)
    if results:
        pd.DataFrame(results).to_csv(FIBRO2_DIR / "ta_results.csv", index=False)
    primary = next((r for r in results if r.get("accession") == GSE_PRIMARY), None)
    sc_row = next((r for r in results if r.get("role") == "single_cell_external"), None)
    reading = _ta_reading(primary, sc_row, sc_resolved, rec_p)
    summary = dict(
        primary=primary, single_cell=sc_row, sc_resolved=sc_resolved,
        sc_reason=sc_reason, reading=reading, n_results=len(results),
        n_cohort_rows=len(cohort_rows), seed=FIBRO2_SEED, boot_seed=FIBRO2_BOOT,
    )
    dump_json(FIBRO2_DIR / "ta_summary.json", jsonable(summary))
    man = load_manifest()
    man["status"] = "TA_DONE"
    man["ta_reading"] = reading
    save_manifest(man)
    progress_snapshot("T-B QC + d0→d7 null", stop="TA_DONE")
    if close_log:
        log.close()
    return summary


if __name__ == "__main__":
    try:
        run_ta()
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise

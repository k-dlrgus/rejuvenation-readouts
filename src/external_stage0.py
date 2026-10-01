"""Stage 0 — SEA-AD MTG donor metadata, nucleus labels, subclass→DLPFC mapping.

Usage: python src/external_stage0.py
Does not freeze directions or project. Does not modify prior FINDINGS*.md.
"""
from __future__ import annotations

import re
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from external_common import (  # noqa: E402
    EXT_DIR, EXT_RAW, EXT_SEED, DATASET_ID, COLLECTION_ID, TITLE_REQUIRED,
    TISSUE_REQUIRED, OUR_DATA_URL, CXG_TABLE, CXG_API, MIN_DONORS_MAP,
    MANIFEST_PATH, STOP_PATH, PREREG_DATE,
    Logger, dump_json, load_json, StopStep, jsonable, ext_log_banner,
    load_dlpfc_types, require_columns, resolve_column, pick_pathology_field, encode_pathology,
    refresh_progress, save_manifest, pathology_columns,
)
from tissue_peek_obs import _open_h5, obs_column_names, read_obs_column  # noqa: E402
from search_cxg import parse_age_years  # noqa: E402
from download import download as http_download  # noqa: E402
from tissue_download_top import local_path as dlpfc_local_path  # noqa: E402


def _record_failure(man, step, message, details=None):
    rec = dict(step=step, message=message, details=jsonable(details or {}))
    man.setdefault("failures", []).append(rec)
    man["status"] = "STOP"
    save_manifest(man)
    STOP_PATH.write_text(f"stop=True\nstep={step}\nmessage={message}\n", encoding="utf-8")
    return rec


def verify_cxg(log) -> dict:
    """Live collection record. Dataset id must match the fetched T1 table."""
    if not CXG_TABLE.exists():
        raise StopStep(
            "cxg_verify",
            f"fetched CELLxGENE table missing: {CXG_TABLE}. Not inventing a URL.",
            dict(path=str(CXG_TABLE)),
        )
    tab = pd.read_csv(CXG_TABLE)
    if "dataset_id" not in tab.columns:
        raise StopStep(
            "cxg_verify",
            f"{CXG_TABLE.name} has no dataset_id. columns={list(tab.columns)}",
            dict(columns=list(tab.columns)),
        )
    hit = tab[tab.dataset_id.astype(str) == DATASET_ID]
    if hit.empty:
        raise StopStep(
            "cxg_verify",
            f"dataset_id {DATASET_ID} not in fetched {CXG_TABLE.name}. Not substituting another SEA-AD object.",
            dict(n_rows=int(len(tab))),
        )
    row = hit.iloc[0]
    table_url = str(row.h5ad_url) if "h5ad_url" in hit.columns else None
    table_title = str(row.title) if "title" in hit.columns else None
    log(f"[cxg] fetched table title={table_title}")
    log(f"[cxg] fetched table url={table_url}")

    url = f"{CXG_API}/collections/{COLLECTION_ID}"
    log(f"[cxg] GET {url}")
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    col = r.json()
    hits = [x for x in (col.get("datasets") or []) if x.get("dataset_id") == DATASET_ID]
    if not hits:
        raise StopStep(
            "cxg_verify",
            f"dataset {DATASET_ID} not in live collection {COLLECTION_ID}",
            dict(n_datasets=len(col.get("datasets") or [])),
        )
    d = hits[0]
    h5ad = next((a for a in d.get("assets") or [] if a.get("filetype") == "H5AD"), None)
    tissues = [t.get("label") for t in (d.get("tissue") or [])]
    title = d.get("title")
    if title != TITLE_REQUIRED:
        raise StopStep(
            "cxg_verify",
            f"live title {title!r} != required {TITLE_REQUIRED!r}. Not substituting.",
            dict(live=title, required=TITLE_REQUIRED),
        )
    if TISSUE_REQUIRED not in tissues:
        raise StopStep(
            "cxg_verify",
            f"required tissue {TISSUE_REQUIRED!r} not in live tissues {tissues}. Not substituting.",
            dict(tissues=tissues),
        )
    if h5ad is None or not h5ad.get("url"):
        raise StopStep("cxg_verify", "live record has no H5AD asset", dict(assets=d.get("assets")))
    rec = dict(
        dataset_id=d.get("dataset_id"),
        collection_id=COLLECTION_ID,
        collection_name=col.get("name"),
        collection_doi=col.get("doi") or col.get("collection_doi"),
        title=title,
        cell_count=d.get("cell_count"),
        n_donors_api=len(d.get("donor_id") or []),
        tissues=tissues,
        diseases=[x.get("label") for x in (d.get("disease") or [])],
        cell_types_api=[c.get("label") for c in (d.get("cell_type") or [])],
        assays=[a.get("label") for a in (d.get("assay") or [])],
        h5ad_url=h5ad.get("url"),
        h5ad_bytes=h5ad.get("filesize"),
        fetched_table_url=table_url,
        fetched_table_title=table_title,
        citation=d.get("citation"),
    )
    log(f"[cxg] title={rec['title']}")
    log(f"[cxg] cells={rec['cell_count']}  donors_api={rec['n_donors_api']}  types_api={len(rec['cell_types_api'])}")
    log(f"[cxg] tissues={rec['tissues']}  diseases={rec['diseases']}")
    log(f"[cxg] h5ad_url={rec['h5ad_url']}")
    log(f"[cxg] h5ad_bytes={rec['h5ad_bytes']}")
    if table_url and table_url != rec["h5ad_url"]:
        log("[cxg] live H5AD URL differs from fetched table; using LIVE URL (verified this fetch)")
    return rec


def find_donor_metadata_xlsx(log) -> dict:
    """The Our Data page lists the donor-metadata xlsx. Require exactly one such link."""
    log(f"[meta] GET {OUR_DATA_URL}")
    r = requests.get(OUR_DATA_URL, timeout=120)
    r.raise_for_status()
    hrefs = re.findall(r"""href=["']([^"']+)["']""", r.text, flags=re.I)
    hits = []
    for u in hrefs:
        lu = u.lower()
        if "donor" in lu and "metadata" in lu and lu.endswith((".xlsx", ".xls")):
            hits.append(u)
        elif "donor_metadata" in lu and lu.endswith((".xlsx", ".xls")):
            hits.append(u)
    hits = list(dict.fromkeys(hits))
    log(f"[meta] donor-metadata xlsx hrefs: {hits}")
    if len(hits) != 1:
        raise StopStep(
            "donor_metadata",
            f"Our Data page must have exactly one donor-metadata xlsx link, found {len(hits)}: {hits}. "
            "Not guessing among files.",
            dict(hits=hits, n_hrefs=len(hrefs)),
        )
    return dict(url=hits[0], page=OUR_DATA_URL, n_hrefs=len(hrefs))


def download_xlsx(url: str, log) -> Path:
    log(f"[download] donor metadata: {url}")
    path = Path(http_download(url, str(EXT_RAW)))
    log(f"[download] wrote {path} size={path.stat().st_size}")
    return path


def read_excel_log(path: Path, log) -> tuple[pd.DataFrame, list[str]]:
    df = pd.read_excel(path, dtype="object")
    cols = [str(c) for c in df.columns]
    log(f"[read] {path.name}: rows={len(df):,} cols={len(cols)}")
    log(f"[read] {path.name} columns actually read ({len(cols)}): {cols}")
    return df, cols


def join_donor_key(cxg_donors: set[str], meta: pd.DataFrame, log) -> str:
    """Find the unique metadata column that covers the most CELLxGENE donor_ids exactly."""
    n = len(cxg_donors)
    scored = []
    for c in meta.columns:
        vals = set(meta[c].dropna().astype(str).str.strip().tolist())
        vals.discard("")
        ov = len(cxg_donors & vals)
        if ov:
            scored.append((ov, c))
    scored.sort(reverse=True)
    log(f"[join] CXG donors={n}  metadata columns by exact overlap: {scored[:12]}")
    if not scored:
        raise StopStep(
            "donor_join",
            "no metadata column shares any exact-string donor_id with CELLxGENE. "
            "Not stripping prefixes or fuzzy-matching. "
            f"example_cxg={sorted(cxg_donors)[:8]}",
            dict(n_cxg=n, example_cxg=sorted(cxg_donors)[:8]),
        )
    best_n = scored[0][0]
    covering = [c for ov, c in scored if ov == best_n]
    if len(covering) > 1:
        raise StopStep(
            "donor_join",
            f"multiple metadata columns share the same max exact overlap ({best_n}): {covering}. Not picking one.",
            dict(covering=covering, overlap=best_n),
        )
    return covering[0]


def dlpfc_nucleus_subclass_celltype(log) -> pd.DataFrame:
    """Per-nucleus subclass × cell_type from the local DLPFC h5ad. Fail if the file is missing."""
    path = dlpfc_local_path("brain_aging")
    if not path.exists():
        raise StopStep(
            "dlpfc_nucleus",
            f"local DLPFC h5ad missing: {path}. Needed for nucleus-level subclass×cell_type (not the pseudobulk first()).",
            dict(path=str(path)),
        )
    log(f"[dlpfc] peeking nucleus subclass×cell_type from {path}")
    h, handle = _open_h5(str(path))
    try:
        cols = obs_column_names(h)
        require_columns(cols, ["subclass", "cell_type"], file_tag="dlpfc_h5ad", step="dlpfc_nucleus")
        sub = read_obs_column(h, "subclass")
        ct = read_obs_column(h, "cell_type")
    finally:
        h.close()
        if handle is not None:
            handle.close()
    tab = pd.DataFrame({"subclass": sub, "cell_type": ct})
    g = tab.groupby(["subclass", "cell_type"], dropna=False).size().rename("n_nuclei").reset_index()
    log(f"[dlpfc] nucleus subclass×cell_type pairs={len(g)}  subclasses={tab.subclass.nunique()}  cell_types={tab.cell_type.nunique()}")
    return g


def build_mapping(sea_pairs: pd.DataFrame, dlpfc_types: list[str], dlpfc_pairs: pd.DataFrame, log) -> tuple[pd.DataFrame, pd.DataFrame]:
    """SEA-AD subclass → DLPFC type by exact strings only.

    A subclass maps iff it has exactly one SEA-AD cell_type AND that cell_type is
    one of the 20 DLPFC types. A 1:1 DLPFC nucleus subclass string match is
    reported as a second key but does not override a cell_type mismatch.
    Multiple cell_types in one subclass → STOP (ambiguous). Never synonymize.
    """
    dlpfc_set = set(dlpfc_types)
    dlpfc_sub_to_ct = {}
    ambiguous_dlpfc_sub = []
    for sub, g in dlpfc_pairs.groupby("subclass"):
        cts = sorted(g.cell_type.astype(str).unique().tolist())
        in20 = [c for c in cts if c in dlpfc_set]
        if len(in20) == 1:
            dlpfc_sub_to_ct[str(sub)] = in20[0]
        elif len(in20) > 1:
            ambiguous_dlpfc_sub.append((str(sub), in20))
    if ambiguous_dlpfc_sub:
        log(f"[map] DLPFC nucleus subclasses with >1 of the 20 types (not used as a key): {ambiguous_dlpfc_sub[:12]}")

    rows, unmapped = [], []
    for sub, g in sea_pairs.groupby("subclass"):
        cts = sorted(g.cell_type.astype(str).unique().tolist())
        n_nuc = int(g.n_nuclei.sum()) if "n_nuclei" in g.columns else int(g.n.sum()) if "n" in g.columns else None
        n_col = "n_nuclei" if "n_nuclei" in g.columns else ("n" if "n" in g.columns else None)
        if len(cts) != 1:
            raise StopStep(
                "mapping",
                f"SEA-AD subclass {sub!r} has {len(cts)} cell_type values {cts}. "
                "Ambiguous; not majority-voting.",
                dict(subclass=str(sub), cell_types=cts),
            )
        ct = cts[0]
        key = None
        dlpfc_type = None
        if ct in dlpfc_set:
            dlpfc_type = ct
            key = "cell_type_exact"
        elif str(sub) in dlpfc_sub_to_ct:
            dlpfc_type = dlpfc_sub_to_ct[str(sub)]
            key = "subclass_exact_to_dlpfc_subclass"
        rec = dict(
            seaad_subclass=str(sub),
            seaad_cell_type=ct,
            dlpfc_type=dlpfc_type,
            mapping_key=key,
            n_nuclei=int(g[n_col].sum()) if n_col else n_nuc,
            mapped=bool(dlpfc_type is not None),
        )
        if dlpfc_type is None:
            unmapped.append(rec)
        rows.append(rec)
        log(f"[map] subclass={sub!r} cell_type={ct!r} → {dlpfc_type!r} ({key or 'UNMAPPED'})")
    mapped = pd.DataFrame(rows)
    um = pd.DataFrame(unmapped) if unmapped else pd.DataFrame(columns=mapped.columns)
    n_map = int(mapped.mapped.sum()) if len(mapped) else 0
    log(f"[map] subclasses mapped={n_map}/{len(mapped)}  unmapped={len(um)}")
    if n_map < 1:
        raise StopStep(
            "mapping",
            "zero SEA-AD subclasses mapped onto the 20 DLPFC types by exact cell_type "
            "or exact subclass string. Not inventing a synonym table.",
            dict(seaad_subclasses=mapped.seaad_subclass.tolist() if len(mapped) else []),
        )
    return mapped, um


def stage0():
    rng_note = EXT_SEED  # recorded; Stage 0 does not draw
    log = Logger(EXT_DIR / "stage0_report.txt")
    man = dict(
        status="RUNNING", seed=EXT_SEED, dataset_id=DATASET_ID, collection_id=COLLECTION_ID,
        columns_read={}, files=[], failures=[], started=time.strftime("%Y-%m-%dT%H:%M:%S"),
    )
    save_manifest(man)
    refresh_progress(next_action="Stage 0 running", stop="none")
    try:
        return _stage0(log, man)
    except StopStep as e:
        log(f"[STOP] {e.step}: {e.message}")
        _record_failure(man, e.step, e.message, e.details)
        refresh_progress(next_action="STOP — write FINDINGS_EXTERNAL.md", stop=f"STOP {e.step}")
        raise
    finally:
        log.close()


def _stage0(log, man):
    if STOP_PATH.exists():
        STOP_PATH.unlink()
    ext_log_banner(log, "Stage 0 — SEA-AD MTG metadata and subclass mapping")
    dlpfc_types, dlpfc_obs = load_dlpfc_types(log)
    man["dlpfc_types"] = dlpfc_types
    man["n_dlpfc_types"] = len(dlpfc_types)
    man["n_dlpfc_donors"] = int(dlpfc_obs.donor.nunique())

    cxg = verify_cxg(log)
    man["cxg"] = cxg
    man["files"].append(dict(kind="cxg_h5ad", url=cxg["h5ad_url"], bytes=cxg["h5ad_bytes"]))

    meta_loc = find_donor_metadata_xlsx(log)
    xlsx_path = download_xlsx(meta_loc["url"], log)
    meta, meta_cols = read_excel_log(xlsx_path, log)
    man["columns_read"]["donor_metadata_xlsx"] = meta_cols
    man["files"].append(dict(kind="donor_metadata_xlsx", url=meta_loc["url"], path=str(xlsx_path),
                             n_rows=int(len(meta)), columns=meta_cols))

    dlpfc_pairs = dlpfc_nucleus_subclass_celltype(log)
    dlpfc_pairs.to_csv(EXT_DIR / "stage0_dlpfc_nucleus_subclass_celltype.csv", index=False)

    log(f"[obs] opening SEA-AD h5ad (obs only): {cxg['h5ad_url']}")
    h, handle = _open_h5(cxg["h5ad_url"])
    try:
        cols = obs_column_names(h)
        log(f"[obs] SEA-AD obs columns ({len(cols)}): {cols}")
        man["columns_read"]["seaad_h5ad_obs"] = cols
        require_columns(
            cols,
            ["donor_id", "cell_type", "sex", "development_stage"],
            file_tag="seaad_h5ad_obs",
            step="obs",
        )
        subclass_col = resolve_column(
            cols, ("subclass", "Subclass"),
            field="subclass", step="obs", file_tag="seaad_h5ad_obs",
        )
        log(f"[obs] subclass column in this file: {subclass_col!r}")
        extra = [c for c in (
            "disease", "tissue", "assay", "suspension_type", "is_primary_data",
            "cell_type_ontology_term_id", "self_reported_ethnicity",
            "Neurotypical reference", "Class", "Supertype", "Age at death",
        ) if c in cols]
        path_in_obs = pathology_columns(cols)
        wanted = ["donor_id", subclass_col, "cell_type", "sex", "development_stage"] + extra + path_in_obs
        wanted = list(dict.fromkeys(wanted))
        arrays = {}
        for c in wanted:
            log(f"[obs] read {c} ...")
            arrays[c] = read_obs_column(h, c)
        n_nuc = len(arrays["donor_id"])
        log(f"[obs] n_nuclei={n_nuc:,}")
        if cxg["cell_count"] is not None and int(cxg["cell_count"]) != n_nuc:
            raise StopStep(
                "obs",
                f"obs length {n_nuc} != live cell_count {cxg['cell_count']}",
                dict(n_obs=n_nuc, cell_count=cxg["cell_count"]),
            )
    finally:
        h.close()
        if handle is not None:
            handle.close()

    donor = np.asarray(arrays["donor_id"]).astype(str)
    subclass = np.asarray(arrays[subclass_col]).astype(str)
    cell_type = np.asarray(arrays["cell_type"]).astype(str)
    sex = np.asarray(arrays["sex"]).astype(str)
    stage = np.asarray(arrays["development_stage"])
    age = np.array([parse_age_years(s) for s in stage], dtype=float)
    man["subclass_column"] = subclass_col

    # subclass × cell_type
    pair = pd.DataFrame({"subclass": subclass, "cell_type": cell_type})
    sea_pairs = pair.groupby(["subclass", "cell_type"], dropna=False).size().rename("n_nuclei").reset_index()
    sea_pairs.to_csv(EXT_DIR / "stage0_seaad_subclass_celltype.csv", index=False)
    log(f"[obs] SEA-AD subclasses={sea_pairs.subclass.nunique()}  cell_types={sea_pairs.cell_type.nunique()}")

    mapped, unmapped = build_mapping(sea_pairs, dlpfc_types, dlpfc_pairs, log)
    mapped.to_csv(EXT_DIR / "stage0_mapping.csv", index=False)
    unmapped.to_csv(EXT_DIR / "stage0_unmapped.csv", index=False)
    sub_to_dlpfc = dict(
        zip(mapped.loc[mapped.mapped, "seaad_subclass"], mapped.loc[mapped.mapped, "dlpfc_type"])
    )

    dlpfc_label = pd.Series(subclass).map(sub_to_dlpfc)
    mapped_mask = dlpfc_label.notna().to_numpy()
    n_mapped_nuc = int(mapped_mask.sum())
    log(f"[map] nuclei in mapped subclasses: {n_mapped_nuc:,} / {len(donor):,}")

    # donor-level from nuclei
    nuc = pd.DataFrame({
        "donor": donor,
        "subclass": subclass,
        "cell_type": cell_type,
        "dlpfc_type": dlpfc_label.to_numpy(),
        "sex": sex,
        "development_stage": np.asarray(stage).astype(str),
        "age": age,
    })
    if "disease" in arrays:
        nuc["disease"] = np.asarray(arrays["disease"]).astype(str)
    if "tissue" in arrays:
        nuc["tissue"] = np.asarray(arrays["tissue"]).astype(str)
    if "is_primary_data" in arrays:
        nuc["is_primary_data"] = arrays["is_primary_data"]
    if "Neurotypical reference" in arrays:
        nuc["Neurotypical reference"] = arrays["Neurotypical reference"]
    for c in pathology_columns(list(arrays.keys())):
        nuc[c] = np.asarray(arrays[c]).astype(str)
    if "Age at death" in arrays:
        nuc["Age at death"] = pd.to_numeric(pd.Series(arrays["Age at death"]), errors="coerce")

    # within-donor consistency
    check_const = ["sex", "age", "development_stage"]
    for c in (
        "disease", "tissue", "Neurotypical reference", "Age at death",
        *pathology_columns(list(nuc.columns)),
    ):
        if c in nuc.columns:
            check_const.append(c)
    for c in check_const:
        nun = nuc.groupby("donor")[c].nunique(dropna=True)
        bad = nun[nun > 1]
        if len(bad):
            raise StopStep(
                "donor_consistency",
                f"{c} is not constant within donor for {len(bad)} donors. Not repairing. example={list(bad.head().index)}",
                dict(field=c, n_bad=int(len(bad))),
            )

    extra_agg = {}
    for c in ("disease", "tissue", "Neurotypical reference", "Age at death"):
        if c in nuc.columns:
            extra_agg[c] = (c, "first")
    for c in pathology_columns(list(nuc.columns)):
        extra_agg[c] = (c, "first")

    don = nuc.groupby("donor", dropna=False).agg(
        age=("age", "first"),
        sex=("sex", "first"),
        development_stage=("development_stage", "first"),
        n_nuclei=("donor", "size"),
        n_nuclei_mapped=("dlpfc_type", lambda s: int(s.notna().sum())),
        n_mapped_types=("dlpfc_type", lambda s: int(s.dropna().nunique())),
        **extra_agg,
    ).reset_index()
    mapped_donors = don.loc[don.n_nuclei_mapped > 0, "donor"].astype(str)
    n_map_don = int(mapped_donors.nunique())
    log(f"[cohort] donors total={len(don)}  donors with ≥1 mapped nucleus={n_map_don}")

    per_type = (
        nuc.dropna(subset=["dlpfc_type"])
        .groupby("dlpfc_type")
        .agg(n_nuclei=("donor", "size"), n_donors=("donor", "nunique"), n_subclasses=("subclass", "nunique"))
        .reset_index()
        .rename(columns={"dlpfc_type": "dlpfc_type"})
    )
    # include unmapped DLPFC types with zeros
    missing_types = [t for t in dlpfc_types if t not in set(per_type.dlpfc_type)]
    if missing_types:
        per_type = pd.concat([
            per_type,
            pd.DataFrame({"dlpfc_type": missing_types, "n_nuclei": 0, "n_donors": 0, "n_subclasses": 0}),
        ], ignore_index=True)
    per_type = per_type.sort_values("n_donors", ascending=False)
    per_type.to_csv(EXT_DIR / "stage0_n_donors_per_type.csv", index=False)
    log("[cohort] n donors per DLPFC type (mapped subclasses):")
    log(per_type.to_string(index=False))

    # join metadata
    cxg_donors = set(don.donor.astype(str))
    key = join_donor_key(cxg_donors, meta, log)
    log(f"[join] using metadata column {key!r}")
    meta = meta.copy()
    meta["_join"] = meta[key].astype(str).str.strip()
    dup = meta["_join"].duplicated(keep=False)
    if dup.any():
        raise StopStep(
            "donor_join",
            f"metadata column {key!r} has duplicate donor ids ({int(dup.sum())} rows). Not collapsing.",
            dict(key=key, n_dup=int(dup.sum())),
        )
    merged = don.merge(meta, left_on="donor", right_on="_join", how="left", suffixes=("", "_meta"))
    not_in_meta = sorted(cxg_donors - set(meta["_join"]))
    extra_meta = sorted(set(meta["_join"]) - cxg_donors)
    log(f"[join] CXG donors not in metadata={len(not_in_meta)}  metadata donors not in CXG={len(extra_meta)}")
    if not_in_meta:
        path_on_don = pathology_columns(list(don.columns))
        leftover = don[don.donor.astype(str).isin(not_in_meta)]
        covered = False
        n_uncovered = len(leftover)
        if path_on_don:
            miss_path = leftover[path_on_don].apply(
                lambda s: s.isna() | s.astype(str).str.strip().isin(("", "nan", "<NA>"))
            )
            n_uncovered = int(miss_path.all(axis=1).sum())
            covered = n_uncovered == 0
            log(f"[join] leftover donors with no h5ad pathology in {path_on_don}: {n_uncovered}")
        if not covered:
            raise StopStep(
                "donor_join",
                f"{len(not_in_meta)} CELLxGENE donors have no metadata row under {key!r} "
                "and do not have ADNC/Braak/CERAD on the h5ad. Not dropping them silently.",
                dict(n=len(not_in_meta), example=not_in_meta[:10], key=key,
                     path_on_don=path_on_don, n_uncovered=n_uncovered),
            )
        log(f"[join] leftover {len(not_in_meta)} donors kept: h5ad already has {path_on_don}. "
            f"example={not_in_meta[:8]}")
        man["cxg_donors_not_in_xlsx"] = not_in_meta
    man["columns_read"]["donor_metadata_join_key"] = key
    man["n_metadata_donors_not_in_cxg"] = len(extra_meta)

    path_field, path_family, path_hits = pick_pathology_field(list(merged.columns), log)
    man["pathology_hits"] = path_hits
    man["pathology_field"] = path_field
    man["pathology_family"] = path_family
    obs_cols = list(man["columns_read"].get("seaad_h5ad_obs") or [])
    xlsx_cols = list(man["columns_read"].get("donor_metadata_xlsx") or [])
    if path_field in obs_cols:
        path_src, path_src_file = "seaad_h5ad_obs", cxg["h5ad_url"]
    elif path_field in xlsx_cols:
        path_src, path_src_file = "donor_metadata_xlsx", str(xlsx_path)
    else:
        raise StopStep(
            "pathology",
            f"chosen field {path_field!r} is not in the h5ad obs columns or the xlsx columns "
            "(possible merge suffix). Not guessing.",
            dict(path_field=path_field, obs_hits=[c for c in obs_cols if c in path_hits],
                 xlsx_hits=[c for c in xlsx_cols if c in path_hits]),
        )
    score, is_low, enc = encode_pathology(merged[path_field], path_family)
    merged["pathology_score"] = score.to_numpy()
    merged["pathology_low"] = is_low.to_numpy()
    log(f"[pathology] field={path_field!r} family={path_family}  n_low={enc['n_low']}  "
        f"n_missing={enc['n_missing']}  observed={enc['observed']}")

    n_age_na = int(merged.age.isna().sum())
    if n_age_na:
        bad = sorted(merged.loc[merged.age.isna(), "development_stage"].astype(str).unique().tolist())
        raise StopStep(
            "age",
            f"{n_age_na} donors have unparsed development_stage. Not substituting Age at death. "
            f"unparsed_levels={bad}",
            dict(n=n_age_na, unparsed=bad),
        )
    ages = merged.age.dropna()
    if ages.empty:
        raise StopStep("age", "no numeric ages after parse_age_years on development_stage")
    sex_counts = merged.sex.astype(str).value_counts(dropna=False).to_dict()
    log(f"[cohort] age min={float(ages.min())} max={float(ages.max())} nunique={int(ages.nunique())} "
        f"missing={int(merged.age.isna().sum())}")
    log(f"[cohort] sex counts: {sex_counts}")

    if n_map_don < MIN_DONORS_MAP:
        raise StopStep(
            "n_donors_map",
            f"only {n_map_don} donors have ≥1 nucleus in a mapped subclass "
            f"(bar ≥ {MIN_DONORS_MAP}). Not lowering the bar.",
            dict(n_mapped_donors=n_map_don, bar=MIN_DONORS_MAP),
        )

    merged.to_csv(EXT_DIR / "stage0_donor_meta.csv", index=False)
    nuc_map_counts = (
        nuc.dropna(subset=["dlpfc_type"])
        .groupby(["donor", "dlpfc_type", "subclass"])
        .size().rename("n_nuclei").reset_index()
    )
    nuc_map_counts.to_csv(EXT_DIR / "stage0_donor_type_nucleus_counts.csv", index=False)

    summary = dict(
        dataset_id=DATASET_ID,
        title=TITLE_REQUIRED,
        tissue=TISSUE_REQUIRED,
        n_nuclei=int(len(donor)),
        n_donors=int(len(merged)),
        n_donors_mapped=n_map_don,
        n_subclasses=int(sea_pairs.subclass.nunique()),
        n_subclasses_mapped=int(mapped.mapped.sum()),
        n_subclasses_unmapped=int((~mapped.mapped).sum()) if "mapped" in mapped else int(len(unmapped)),
        unmapped_subclasses=unmapped.seaad_subclass.tolist() if len(unmapped) else [],
        age_min=float(ages.min()),
        age_max=float(ages.max()),
        nunique_age=int(ages.nunique()),
        sex_counts={str(k): int(v) for k, v in sex_counts.items()},
        pathology_field=path_field,
        pathology_family=path_family,
        pathology_source=path_src,
        pathology_source_file=path_src_file,
        pathology_observed=enc["observed"],
        n_pathology_low=int(enc["n_low"]),
        n_pathology_missing=int(enc["n_missing"]),
        pathology_low_levels=enc["low_levels_used"],
        h5ad_url=cxg["h5ad_url"],
        h5ad_bytes=cxg["h5ad_bytes"],
        donor_metadata_path=str(xlsx_path),
        donor_join_key=key,
        min_donors_map_bar=MIN_DONORS_MAP,
        stop=False,
    )
    dump_json(EXT_DIR / "stage0_summary.json", jsonable(summary))
    man.update(dict(status="STAGE0_OK", summary=summary, columns_read=man["columns_read"]))
    save_manifest(man)
    log("[stage0] complete. Pre-register before any fit.")
    refresh_progress(next_action="write pre-registration flag (before any fit)", stop="none",
                     extra=f"Stage 0: n_donors_mapped={n_map_don}  pathology={path_field!r}  "
                           f"age {summary['age_min']:.0f}–{summary['age_max']:.0f}  "
                           f"unmapped_subclasses={summary['unmapped_subclasses']}")
    return summary


if __name__ == "__main__":
    try:
        stage0()
    except StopStep:
        sys.exit(2)

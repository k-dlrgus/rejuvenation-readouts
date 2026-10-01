"""STAGE R0 — peek DLPFC Aging_Cohort obs WITHOUT downloading the 12 GB matrix.

Verifies the accession/URL against the fetched CELLxGENE curation API and
results/tissue/t1_cxg_primary_candidates.csv, then opens the published HTTPS
h5ad with fsspec range requests and scores every obs column for:

  - 1:1 alias of donor (within-batch metric undefined if ALL candidates alias)
  - shared-donor technical batch (sequencing pool / run / library)

STOP R0: if every candidate batch column aliases donor 1:1, halt and peek the
retina 104-donor snRNA atlas instead.

Usage: python src/brain_r0.py [brain_aging|retina_sn]
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from brain_common import (  # noqa: E402
    BRAIN_DIR, BRAIN_SEED, Logger, DLPFC_NAME, DLPFC_DATASET_ID,
    RETINA_SN_NAME, RETINA_SN_DATASET_ID, TECH_BATCH_KEYS, SOFT_BATCH_KEYS,
)
from tissue_download_top import resolve  # noqa: E402
from tissue_peek_obs import _open_h5, obs_column_names, read_obs_column  # noqa: E402
from search_cxg import parse_age_years, API  # noqa: E402

CXG_API = API


def _decode(v):
    if isinstance(v, bytes):
        return v.decode("utf-8", "replace")
    return v


def verify_cxg(dataset_id: str, expected_url: str, collection_id: str | None, log) -> dict:
    """Fetch the live curation record. Nothing is assumed from the filename.

    Individual /datasets/{id} 404s on this API; collections/{id} and the /datasets
    listing are the endpoints this project already used in T1.
    """
    d = None
    if collection_id:
        url = f"{CXG_API}/collections/{collection_id}"
        log(f"[R0] GET {url}")
        r = requests.get(url, timeout=120)
        r.raise_for_status()
        col = r.json()
        hits = [x for x in (col.get("datasets") or []) if x.get("dataset_id") == dataset_id]
        if not hits:
            raise RuntimeError(f"dataset {dataset_id} not in collection {collection_id}")
        d = hits[0]
        d.setdefault("collection_id", collection_id)
        d.setdefault("collection_doi", col.get("doi") or col.get("collection_doi"))
        d.setdefault("collection_name", col.get("name"))
        log(f"  collection name={col.get('name')}  doi={col.get('doi')}")
    if d is None:
        url = f"{CXG_API}/datasets"
        log(f"[R0] GET {url} (fallback listing)")
        r = requests.get(url, timeout=300)
        r.raise_for_status()
        hits = [x for x in r.json() if x.get("dataset_id") == dataset_id]
        if not hits:
            raise RuntimeError(f"dataset {dataset_id} not in /datasets listing")
        d = hits[0]
    def _labels(xs):
        out = []
        for x in xs or []:
            if isinstance(x, dict):
                out.append(x.get("label") or x.get("name"))
            else:
                out.append(x)
        return out

    h5ad = next((a for a in d.get("assets", []) if a.get("filetype") == "H5AD"), None)
    rec = dict(
        dataset_id=d.get("dataset_id"),
        collection_id=d.get("collection_id"),
        title=d.get("title"),
        cell_count=d.get("cell_count"),
        n_donors=len(d.get("donor_id") or []),
        n_cell_types=len(_labels(d.get("cell_type"))),
        assays=_labels(d.get("assay")),
        tissues=_labels(d.get("tissue")),
        diseases=_labels(d.get("disease")),
        suspension=_labels(d.get("suspension_type")),
        doi=(d.get("citation") or ""),
        h5ad_url=None if h5ad is None else h5ad.get("url"),
        h5ad_bytes=None if h5ad is None else h5ad.get("filesize"),
        collection_doi=d.get("collection_doi") or d.get("doi"),
        collection_name=d.get("collection_name"),
    )
    log(f"  title={rec['title']}")
    log(f"  dataset_id={rec['dataset_id']}  collection_id={rec['collection_id']}")
    log(f"  cells={rec['cell_count']}  donors={rec['n_donors']}  types={rec['n_cell_types']}")
    log(f"  assays={rec['assays']}  suspension={rec['suspension']}  tissues={rec['tissues']}")
    log(f"  diseases={rec['diseases']}")
    log(f"  h5ad_url={rec['h5ad_url']}")
    log(f"  h5ad_bytes={rec['h5ad_bytes']}")
    if rec["dataset_id"] != dataset_id:
        raise RuntimeError(f"API dataset_id {rec['dataset_id']} != requested {dataset_id}")
    if rec["h5ad_url"] != expected_url:
        log(f"  [warn] live h5ad URL differs from fetched candidate table:\n"
            f"         table={expected_url}\n         live ={rec['h5ad_url']}\n"
            f"         using the LIVE URL (verified this fetch)")
    return rec


def _nunique_from_node(h, name):
    """Cheap nunique: len(categories) for categoricals, else None (need to read)."""
    g = h["obs"]
    if name not in g:
        return None
    node = g[name]
    if "categories" in getattr(node, "keys", lambda: [])():
        return int(len(node["categories"]))
    return None


def _looks_tech(col: str) -> bool:
    c = col.lower()
    return any(k in c for k in TECH_BATCH_KEYS)


def _looks_soft(col: str) -> bool:
    c = col.lower()
    return any(k in c for k in SOFT_BATCH_KEYS)


def peek_batch_recoverable(path_or_url: str, name: str, log) -> dict:
    log(f"\n{'=' * 90}\nR0 PEEK {name}\n  source: {path_or_url}\n{'=' * 90}")
    h, handle = _open_h5(path_or_url)
    try:
        shape = None
        for p in ("X", "raw/X"):
            if p in h:
                sh = h[p].attrs.get("shape")
                if sh is not None:
                    shape = tuple(int(x) for x in sh)
                    log(f"  matrix {p} shape={shape}")
                    break
        cols = obs_column_names(h)
        log(f"  obs columns ({len(cols)}): {cols}")

        cheap = []
        for c in cols:
            nu = _nunique_from_node(h, c)
            cheap.append(dict(column=c, nunique_categories=nu, looks_technical=_looks_tech(c),
                              looks_soft_batch=_looks_soft(c)))
        cheap_df = pd.DataFrame(cheap)
        log("  categorical nunique (from categories, no codes read):")
        log(cheap_df.to_string(index=False))

        if "donor_id" not in cols:
            log("[R0] NO donor_id column — cannot score shared-donor batches.")
            return dict(name=name, n_cells=None if shape is None else shape[0],
                        obs_columns=cols, batch_recoverable=False, stop_r0=True,
                        reason="no donor_id")

        log("  reading donor_id + candidate columns (codes only)...")
        donor = read_obs_column(h, "donor_id")
        n = len(donor)
        n_donors = int(pd.Series(donor).nunique(dropna=True))
        log(f"  n_cells={n:,}  n_donors={n_donors}")

        # Age / disease / cell type / assay — needed for R0 report and R1 planning.
        extras_to_read = [c for c in (
            "development_stage", "disease", "cell_type", "assay", "suspension_type",
            "tissue", "sex", "is_primary_data", "self_reported_ethnicity",
        ) if c in cols]
        # Every remaining column that could be a batch or a diagnosis (not unique per cell).
        skip_always = {"donor_id", "observation_joinid", "n_genes", "n_counts"}
        for _, row in cheap_df.iterrows():
            c = row.column
            if c in skip_always:
                continue
            nu = row.nunique_categories
            if nu is None or (isinstance(nu, float) and not np.isfinite(nu)):
                extras_to_read.append(c)
            elif 2 <= int(nu) <= max(n_donors * 3, 50) or row.looks_technical or row.looks_soft_batch:
                extras_to_read.append(c)
            elif c.lower() in ("diagnosis", "dx", "disorder", "condition", "phenotype",
                               "primary_diagnosis", "brain_bank", "source", "subclass",
                               "class", "n_counts", "n_genes", "nCount_RNA", "nFeature_RNA",
                               "total_counts", "n_genes_by_counts"):
                extras_to_read.append(c)
        # unique-preserve
        seen, uniq = set(), []
        for c in extras_to_read:
            if c not in seen and c in cols:
                seen.add(c)
                uniq.append(c)

        arrays = {"donor_id": donor}
        for c in uniq:
            log(f"    read {c} ...")
            arrays[c] = read_obs_column(h, c)

        # value counts for key fields
        for c in uniq:
            arr = arrays[c]
            vc = pd.Series(arr).value_counts(dropna=False)
            log(f"  {c}: nunique={vc.shape[0]}  top={list(vc.head(15).items())}")

        # donor-level table
        d = pd.DataFrame({"donor": donor})
        for c in uniq:
            d[c] = arrays[c]
        if "development_stage" in d.columns:
            d["age"] = [parse_age_years(s) for s in d["development_stage"]]
        donor_first = {}
        for c in d.columns:
            if c == "donor":
                continue
            if c == "age":
                donor_first[c] = ("age", "first")
            else:
                donor_first[c] = (c, "first")
        don = d.groupby("donor", dropna=False).agg(**donor_first, n_cells=("donor", "size"))
        log(f"  donors={len(don)}")
        if "age" in don.columns:
            ages = don.age.dropna()
            log(f"  numeric age: n={len(ages)} min={ages.min()} max={ages.max()} "
                f"nunique={ages.nunique()}  missing={int(don.age.isna().sum())}")
            if "development_stage" in don.columns:
                miss = don[don.age.isna()]["development_stage"].astype(str).value_counts().head(12)
                log(f"  non-numeric development_stage (donor-level): {list(miss.items())}")

        # Score every donor-level column for shared-donor batch recoverability.
        skip_score = {
            "age", "n_cells", "n_genes", "n_counts", "observation_joinid",
            "class", "subclass", "subtype", "cell_type", "cell_type_ontology_term_id",
        }
        rows = []
        for c in don.columns:
            if c in skip_score:
                continue
            s = don[c].astype(str)
            n_lev = int(s.nunique(dropna=True))
            sizes = s.value_counts()
            n_shared = int((sizes >= 2).sum())
            frac = float((s.map(sizes) >= 2).mean())
            aliased = n_lev == len(don)
            usable = (n_lev >= 2) and (n_shared >= 2) and (not aliased)
            tech = _looks_tech(c)
            soft = _looks_soft(c)
            rows.append(dict(
                column=c, n_levels=n_lev, n_levels_ge2_donors=n_shared,
                frac_donors_in_shared=frac, aliased_with_donor=bool(aliased),
                looks_technical=bool(tech), looks_soft_batch=bool(soft),
                usable_shared_donor=bool(usable),
                median_donors_per_level=float(sizes.median()) if len(sizes) else np.nan,
                min_donors_per_level=int(sizes.min()) if len(sizes) else 0,
                max_donors_per_level=int(sizes.max()) if len(sizes) else 0,
            ))
        sc = pd.DataFrame(rows).sort_values(
            ["usable_shared_donor", "looks_technical", "frac_donors_in_shared", "n_levels_ge2_donors"],
            ascending=[False, False, False, False],
        )
        log("\n[R0] donor-level batch-column audit (every obs column scored):")
        log(sc.to_string(index=False))

        # Biological covariates are shared-donor but are NOT a technical batch.
        not_batch = {
            "sex", "sex_ontology_term_id", "genetic_ancestry", "self_reported_ethnicity",
            "self_reported_ethnicity_ontology_term_id", "disease", "disease_ontology_term_id",
            "development_stage", "development_stage_ontology_term_id",
        }
        tech_usable = sc[(sc.usable_shared_donor) & (sc.looks_technical)]
        site_usable = sc[(sc.usable_shared_donor) & (sc.looks_soft_batch) & (~sc.column.isin(not_batch))]
        any_usable = sc[(sc.usable_shared_donor) & (~sc.column.isin(not_batch))]
        stop = bool(any_usable.empty)
        preferred = None
        if not tech_usable.empty:
            preferred = str(tech_usable.iloc[0].column)
            kind = "technical_pool"
        elif not site_usable.empty:
            preferred = str(site_usable.iloc[0].column)
            kind = "shared_donor_site_or_assay_not_sequencing_pool"
        elif not any_usable.empty:
            preferred = str(any_usable.iloc[0].column)
            kind = "shared_donor_but_not_sequencing_pool"
        else:
            kind = "none"

        log(f"\n[R0] preferred batch column: {preferred}  kind={kind}")
        log(f"[R0] STOP_R0={stop}  (True means within-batch metric is undefined)")
        if stop:
            log("[R0] every candidate aliases donor 1:1 or has <2 shared levels. "
                "This atlas cannot be corrected the way OneK1K and white matter were.")

        out = dict(
            name=name,
            source=path_or_url,
            n_cells=int(n),
            n_genes=None if shape is None else int(shape[1]),
            n_donors=int(n_donors),
            obs_columns=cols,
            batch_column=preferred,
            batch_kind=kind,
            batch_recoverable=not stop,
            stop_r0=stop,
            n_technical_usable=int(len(tech_usable)),
            n_any_usable=int(len(any_usable)),
        )
        if "age" in don.columns and don.age.notna().any():
            out.update(age_min=float(don.age.min()), age_max=float(don.age.max()),
                       nunique_age=int(don.age.nunique(dropna=True)))
        sc.to_csv(BRAIN_DIR / f"r0_{name}_batch_columns.csv", index=False)
        don.reset_index().to_csv(BRAIN_DIR / f"r0_{name}_donor_meta.csv", index=False)
        return out
    finally:
        h.close()
        if handle is not None:
            handle.close()


def run_one(name: str, log) -> dict:
    info = resolve(name)
    # collection_id lives in the fetched candidate table (not invented)
    from tissue_common import TISSUE_DIR
    tab = pd.read_csv(TISSUE_DIR / "t1_cxg_primary_candidates.csv")
    hit = tab[tab.dataset_id == info["dataset_id"]]
    collection_id = None if hit.empty else str(hit.iloc[0].collection_id)
    log("=" * 100)
    log(f"R0  name={name}  seed={BRAIN_SEED}  dataset_id={info['dataset_id']}")
    log(f"  title={info['title']}  collection_id={collection_id}")
    log(f"  url_from_fetched_table={info['url']}")
    log("=" * 100)
    rec = verify_cxg(info["dataset_id"], info["url"], collection_id, log)
    url = rec["h5ad_url"] or info["url"]
    peek = peek_batch_recoverable(url, name, log)
    peek["cxg"] = rec
    peek["dataset_id"] = info["dataset_id"]
    peek["title"] = info["title"]
    peek["url_used"] = url
    return peek


def main(names=None):
    names = names or sys.argv[1:] or [DLPFC_NAME]
    log = Logger(BRAIN_DIR / "r0_peek.txt")
    results = []
    try:
        log(f"BRAIN R0  seed={BRAIN_SEED}")
        for name in names:
            if name not in (DLPFC_NAME, RETINA_SN_NAME):
                raise SystemExit(f"unknown name {name}; expected {DLPFC_NAME} or {RETINA_SN_NAME}")
            results.append(run_one(name, log))
        dlpfc = next((r for r in results if r["name"] == DLPFC_NAME), None)
        if dlpfc is not None and dlpfc.get("stop_r0"):
            log("\n[STOP R0] DLPFC has no shared-donor technical batch. Peeking retina snRNA as specified.")
            if RETINA_SN_NAME not in names:
                results.append(run_one(RETINA_SN_NAME, log))
        # JSON-safe
        dump = []
        for r in results:
            row = {k: v for k, v in r.items() if k != "cxg"}
            row["cxg"] = {k: v for k, v in (r.get("cxg") or {}).items()
                          if not isinstance(v, (list, dict)) or k in ("assays", "tissues", "diseases", "suspension")}
            dump.append(row)
        (BRAIN_DIR / "r0_summary.json").write_text(json.dumps(dump, indent=1, default=str), encoding="utf-8")
        log(f"\nwrote {BRAIN_DIR / 'r0_summary.json'}")
        return results
    finally:
        log.close()


if __name__ == "__main__":
    main()

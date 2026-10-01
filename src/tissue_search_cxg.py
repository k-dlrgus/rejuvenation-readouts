"""TISSUE SELECTION — Stage T1a: search CELLxGENE Discover for human single-cell / single-nucleus AGING datasets in
LOW-TURNOVER tissues (retina, brain/cortex, skeletal muscle, heart, kidney, liver).

Requirements (from the protocol): per-donor numeric age, span >= 40 y, >= 10 donors, >= 3 cell types, Homo sapiens.
Nothing is assumed: the full /datasets listing is fetched fresh from the curation API and every field printed comes from it.
Outputs: results/tissue/t1_cxg_all_human_datasets.csv (everything), results/tissue/t1_cxg_tissue_candidates.csv (filtered),
         results/tissue/t1_cxg_search.txt (log).
"""
import sys
import re
import requests
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from search_cxg import API, parse_age_years  # noqa: E402  (reused Phase 0 code)
from tissue_common import TISSUE_DIR, Logger  # noqa: E402

TISSUE_PATTERNS = {                  # priority order from the protocol
    "retina": r"\bretina|macula|fovea|retinal",
    "brain": r"brain|cortex|cortical|hippocamp|prefrontal|frontal lobe|temporal lobe|cerebell|striatum|substantia nigra|entorhinal|thalam|hypothalam|midbrain|white matter|pons|medulla oblongata|amygdala|putamen|caudate|dentate gyrus|basal ganglia",
    "skeletal muscle": r"skeletal muscle|vastus|quadriceps|gastrocnemius|rectus abdominis|muscle of",
    "heart": r"\bheart|cardiac|ventricle|atrium|atrial|myocardium|septum",
    "kidney": r"kidney|renal|nephron",
    "liver": r"\bliver|hepatic",
}
EXCLUDE_TISSUE = r"blood|PBMC|bone marrow|spleen|lymph node|thymus|tonsil|cell culture|organoid|embryo"


def tissue_group(tissues: str):
    t = tissues.lower()
    for g, pat in TISSUE_PATTERNS.items():
        if re.search(pat, t):
            return g
    return None


def main():
    log = Logger(TISSUE_DIR / "t1_cxg_search.txt")
    log("Querying CELLxGENE Discover curation API:", f"{API}/datasets")
    r = requests.get(f"{API}/datasets", timeout=600)
    r.raise_for_status()
    ds = r.json()
    log(f"  {len(ds)} datasets returned (all organisms)")
    rows = []
    for d in ds:
        orgs = [o.get("label") for o in d.get("organism", [])]
        if "Homo sapiens" not in orgs:
            continue
        assays = sorted({a.get("label") for a in d.get("assay", [])})
        tissues = sorted({t.get("label") for t in d.get("tissue", [])})
        diseases = sorted({x.get("label") for x in d.get("disease", [])})
        cts = sorted({c.get("label") for c in d.get("cell_type", [])})
        stages = [s.get("label") for s in d.get("development_stage", [])]
        donors = d.get("donor_id", []) or []
        ages = sorted(a for a in (parse_age_years(s) for s in stages) if a is not None)
        suspension = sorted({s for s in (d.get("suspension_type") or []) if s})
        rows.append(dict(
            dataset_id=d.get("dataset_id"), collection_id=d.get("collection_id"), dataset_version_id=d.get("dataset_version_id"),
            title=d.get("title"), collection_name=d.get("collection_name"), cell_count=d.get("cell_count"),
            primary_cell_count=d.get("primary_cell_count"), n_donors=len(donors), n_cell_types=len(cts),
            assays="; ".join(assays), suspension_type="; ".join(suspension), tissues="; ".join(tissues)[:300],
            diseases="; ".join(diseases)[:200], n_dev_stages=len(stages), n_stages_numeric_years=len(ages),
            age_min=ages[0] if ages else None, age_max=ages[-1] if ages else None,
            age_span=(ages[-1] - ages[0]) if ages else None, ages_numeric=";".join(f"{a:g}" for a in ages),
            dev_stages_sample="; ".join(stages[:15]), cell_types_sample="; ".join(cts[:25]),
            h5ad_url=next((a.get("url") for a in d.get("assets", []) if a.get("filetype") == "H5AD"), None),
            h5ad_bytes=next((a.get("filesize") for a in d.get("assets", []) if a.get("filetype") == "H5AD"), None),
            collection_doi=d.get("collection_doi"), citation=d.get("citation"),
        ))
    df = pd.DataFrame(rows)
    df["tissue_group"] = df.tissues.map(tissue_group)
    df.to_csv(TISSUE_DIR / "t1_cxg_all_human_datasets.csv", index=False)
    log(f"  {len(df)} human datasets; {df.tissue_group.notna().sum()} mention a low-turnover target tissue")

    is_rna = df.assays.str.contains("10x|Seq-Well|Drop-seq|Smart-seq|sci-RNA|inDrop|CEL-seq|BD Rhapsody|snRNA|scRNA|microwell", case=False, regex=True)
    not_excluded = ~df.tissues.str.contains(EXCLUDE_TISSUE, case=False, regex=True) | df.tissue_group.notna()
    cand = df[(df.tissue_group.notna()) & is_rna & not_excluded].copy()
    log("\n[filter] human & RNA assay & target tissue                      "
        f"datasets {len(df):,} -> {len(cand):,}")
    n = len(cand); cand = cand[cand.n_donors >= 10]
    log(f"[filter] >= 10 donors                                             datasets {n:,} -> {len(cand):,}")
    n = len(cand); cand = cand[cand.age_span.fillna(0) >= 40]
    log(f"[filter] numeric development_stage ages spanning >= 40 y         datasets {n:,} -> {len(cand):,}")
    n = len(cand); cand = cand[cand.n_cell_types >= 3]
    log(f"[filter] >= 3 cell types                                          datasets {n:,} -> {len(cand):,}")
    order = {g: i for i, g in enumerate(TISSUE_PATTERNS)}
    cand["priority"] = cand.tissue_group.map(order)
    cand["has_normal"] = cand.diseases.str.contains("normal", case=False)
    cand = cand.sort_values(["priority", "n_donors", "cell_count"], ascending=[True, False, False])
    cand.to_csv(TISSUE_DIR / "t1_cxg_tissue_candidates.csv", index=False)
    pd.set_option("display.width", 300); pd.set_option("display.max_colwidth", 60); pd.set_option("display.max_rows", 500)
    log(f"\n=== CELLxGENE candidates in low-turnover tissues (n = {len(cand)}) ===")
    cols = ["tissue_group", "dataset_id", "title", "cell_count", "n_donors", "n_cell_types", "age_min", "age_max", "n_stages_numeric_years",
            "assays", "suspension_type", "diseases", "h5ad_bytes"]
    for g in TISSUE_PATTERNS:
        sub = cand[cand.tissue_group == g]
        log(f"\n--- {g}: {len(sub)} datasets ---")
        if len(sub):
            log(sub[cols].to_string(index=False))
    log.close()
    return cand


if __name__ == "__main__":
    main()

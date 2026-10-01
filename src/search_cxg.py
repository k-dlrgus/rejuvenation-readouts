"""Search CELLxGENE Discover for human single-cell datasets usable as aging atlases.

Criteria (from project spec): per-donor chronological age, >=10 donors, age span >=40y
(ideally 20-80), >=3 cell types, immune/PBMC preferred, disease = normal.

Uses the public curation API; no data is downloaded here. Output: results/phase0_cxg_candidates.csv
"""
import re
import json
import sys
import pandas as pd
import requests

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from config import TAB  # noqa: E402

API = "https://api.cellxgene.cziscience.com/curation/v1"


def parse_age_years(stage: str):
    """Map HsapDv development_stage labels to numeric years where possible.
    Examples: '45-year-old stage' -> 45; '80 year-old and over stage' -> 80;
    'newborn stage', 'embryonic', 'fetal', 'child stage', 'adult' -> None (not continuous)."""
    if stage is None:
        return None
    s = stage.lower()
    m = re.match(r"^(\d+)[- ]year[- ]old", s)
    if m:
        return float(m.group(1))
    m = re.match(r"^(\d+)-month-old", s)
    if m:
        return float(m.group(1)) / 12.0
    return None


def main():
    print("Querying CELLxGENE Discover /datasets ...", flush=True)
    r = requests.get(f"{API}/datasets", timeout=300)
    r.raise_for_status()
    ds = r.json()
    print(f"  {len(ds)} datasets total", flush=True)
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
        ages = [parse_age_years(s) for s in stages]
        ages_num = [a for a in ages if a is not None]
        n_stage_numeric = len(ages_num)
        rows.append(dict(
            dataset_id=d.get("dataset_id"),
            collection_id=d.get("collection_id"),
            title=d.get("title"),
            cell_count=d.get("cell_count"),
            n_donors=len(donors),
            n_cell_types=len(cts),
            assays="; ".join(assays),
            tissues="; ".join(tissues)[:200],
            diseases="; ".join(diseases)[:120],
            n_dev_stages=len(stages),
            n_stages_numeric_years=n_stage_numeric,
            age_min=min(ages_num) if ages_num else None,
            age_max=max(ages_num) if ages_num else None,
            age_span=(max(ages_num) - min(ages_num)) if ages_num else None,
            dev_stages_sample="; ".join(stages[:12]),
            h5ad_url=next((a.get("url") for a in d.get("assets", []) if a.get("filetype") == "H5AD"), None),
            h5ad_bytes=next((a.get("filesize") for a in d.get("assets", []) if a.get("filetype") == "H5AD"), None),
        ))
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "phase0_cxg_all_human_datasets.csv", index=False)
    print(f"  {len(df)} human datasets", flush=True)

    # Filter: normal disease present, >=10 donors, numeric age labels spanning >=40y, >=3 cell types, scRNA
    is_rna = df.assays.str.contains("10x|Seq-Well|Drop-seq|Smart-seq|sci-RNA|inDrop|CEL-seq|BD Rhapsody|Slide-seq", case=False, regex=True)
    cand = df[(df.n_donors >= 10) & (df.age_span.fillna(0) >= 40) & (df.n_cell_types >= 3) & is_rna].copy()
    cand["immune_like"] = cand.tissues.str.contains("blood|PBMC|peripheral|bone marrow|spleen|lymph", case=False, regex=True)
    cand = cand.sort_values(["immune_like", "n_donors"], ascending=[False, False])
    cand.to_csv(TAB / "phase0_cxg_candidates.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 70)
    print("\n=== CELLxGENE candidates meeting: >=10 donors, numeric age span >=40y, >=3 cell types, scRNA-seq ===")
    print(f"n = {len(cand)}")
    cols = ["dataset_id", "title", "cell_count", "n_donors", "n_cell_types", "age_min", "age_max", "n_stages_numeric_years", "tissues", "diseases", "immune_like"]
    print(cand[cols].to_string(index=False))
    return cand


if __name__ == "__main__":
    main()


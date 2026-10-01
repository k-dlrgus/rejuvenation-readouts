import json
import urllib.request
from pathlib import Path

UA = "AgeIdentitySeparabilityBenchmark/scoping (metadata only)"
ids = [
    "a19d1667-a7b5-4556-9e5f-f9bfa690c0f1",  # Human dermal fibroblast atlas
    "b0ef440b-e303-4ad2-ada8-ce13336280ba",  # Multi-modal skin atlas
    "98b68224-956d-45dd-9eca-538747578af8",  # CTCL
]


def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode())


out = {}
for i in ids:
    print("CXG", i, flush=True)
    d = get(f"https://api.cellxgene.cziscience.com/curation/v1/datasets/{i}")
    donors = d.get("donor_id") or []
    stages = d.get("development_stage") or []
    out[i] = {
        "title": d.get("title"),
        "cell_count": d.get("cell_count"),
        "n_donor_id": len(donors) if isinstance(donors, list) else str(type(donors)),
        "donor_id_head": donors[:20] if isinstance(donors, list) else donors,
        "n_development_stage": len(stages) if isinstance(stages, list) else None,
        "development_stage_head": stages[:15] if isinstance(stages, list) else stages,
        "tissue": d.get("tissue"),
        "assay": d.get("assay"),
        "organism": d.get("organism"),
        "collection_name": d.get("collection_name"),
        "citation": (d.get("citation") or "")[:300],
        "explorer_url": d.get("explorer_url"),
        "is_primary_data": d.get("is_primary_data"),
    }
    print(" ", out[i]["title"], "cells", out[i]["cell_count"], "donors", out[i]["n_donor_id"], flush=True)
    print("  stages", out[i]["development_stage_head"], flush=True)

Path("results/scoping/cxg_skin_detail.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
print("wrote cxg_skin_detail.json")

"""Fetch chromosome for Ensembl gene IDs via the Ensembl REST API (POST /lookup/id), cached.

Used to flag sex-linked (chrX / chrY) genes in Phase 1b. Result is cached at
data/processed/ensembl_gene_chrom.csv and reused on rerun.
"""
import json
import time
from pathlib import Path
import pandas as pd
import requests

from config import DATA_PROC

CACHE = DATA_PROC / "ensembl_gene_chrom.csv"
URL = "https://rest.ensembl.org/lookup/id"


def fetch_chromosomes(ensembl_ids, batch=1000, verbose=True) -> pd.DataFrame:
    ids = [i.split(".")[0] for i in ensembl_ids]
    done = {}
    if CACHE.exists():
        c = pd.read_csv(CACHE)
        done = dict(zip(c.gene_id, c.chrom.astype(str)))
    todo = [i for i in ids if i not in done]
    if verbose:
        print(f"[ensembl] {len(ids)} ids; cached={len(ids)-len(todo)}; to fetch={len(todo)}", flush=True)
    for k in range(0, len(todo), batch):
        chunk = todo[k:k + batch]
        for attempt in range(6):
            try:
                r = requests.post(URL, headers={"Content-Type": "application/json", "Accept": "application/json"},
                                  data=json.dumps({"ids": chunk}), timeout=120)
                if r.status_code == 429:
                    time.sleep(float(r.headers.get("Retry-After", 5)))
                    continue
                r.raise_for_status()
                res = r.json()
                break
            except requests.RequestException as e:
                if verbose:
                    print(f"   retry {attempt}: {e}", flush=True)
                time.sleep(5 * (attempt + 1))
        else:
            raise RuntimeError("Ensembl lookup failed repeatedly")
        for gid in chunk:
            rec = res.get(gid)
            done[gid] = str(rec.get("seq_region_name")) if rec else "NA"
        if verbose:
            print(f"   fetched {min(k+batch, len(todo))}/{len(todo)}", flush=True)
        pd.DataFrame({"gene_id": list(done.keys()), "chrom": list(done.values())}).to_csv(CACHE, index=False)
        time.sleep(0.2)
    out = pd.DataFrame({"gene_id": ids, "chrom": [done.get(i, "NA") for i in ids]})
    return out

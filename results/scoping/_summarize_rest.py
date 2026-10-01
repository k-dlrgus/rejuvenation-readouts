import json
import sys
from pathlib import Path
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
d = json.loads(Path("results/scoping/ncbi_searches.json").read_text(encoding="utf-8"))
for k in ["muscle_OSKM_human", "adipose_OSKM_human", "lung_OSKM_human", "PBMC_OSKM_human", "ipsc_neuron_aged", "CRISPRa_fibroblast_human", "CRISPR_screen_fibroblast_human", "LINCS_fibroblast", "TF_atlas_joung"]:
    v = d.get(k) or {}
    print("\n====", k, "count", v.get("count"), "err", v.get("error"))
    n = 0
    for s in (v.get("summaries") or []):
        acc = s.get("accession") or ""
        if not str(acc).startswith("GSE"):
            continue
        n += 1
        if n > 20:
            print("  ...")
            break
        title = (s.get("title") or "")[:130]
        print(f"  {acc:12} n={s.get('n_samples')} {s.get('taxon')} {s.get('gds_type')} {title}")

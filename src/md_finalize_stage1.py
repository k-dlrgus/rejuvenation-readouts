"""Re-evaluate Stage 1 gate from files already on disk; rewrite FINDINGS_MD.md."""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md_common import (  # noqa: E402
    MD_DIR, MD_RAW, RETRIEVAL_DATE, SOURCES_CSV, Logger, dump_json, jsonable,
    load_json, load_manifest, save_manifest,
)
from md_stage1 import evaluate_gate  # noqa: E402
from md_findings import write_findings  # noqa: E402
import pandas as pd


def main():
    log = Logger(MD_DIR / "s1_gate_reval.txt")
    gate = evaluate_gate(log, [])
    log(
        "GATE "
        + json.dumps(
            {k: gate.get(k) for k in (
                "open", "gene_list_recovered", "computation_recovered",
                "n_MD_genes", "reason",
            )},
            ensure_ascii=False,
        )[:2000]
    )
    src = MD_RAW / "eutils_gds_esummary_two_series.json"
    if src.exists():
        js = json.loads(src.read_text(encoding="utf-8"))
        recs = []
        for uid in (js.get("result") or {}).get("uids", []):
            r = js["result"][uid]
            recs.append(dict(
                accession=r.get("accession"),
                title=r.get("title"),
                n_samples=r.get("n_samples"),
                pubmedids=r.get("pubmedids"),
                relations=r.get("relations"),
                extrelations=r.get("extrelations"),
                bioproject=r.get("bioproject"),
                suppfile=r.get("suppfile"),
                pdat=r.get("pdat"),
            ))
        dump_json(MD_DIR / "geo_gds_two_series.json", jsonable(recs))
        log(f"GDS series {recs}")

    extra_path = MD_DIR / "extra_probes.json"
    extra = load_json(extra_path) if extra_path.exists() else []
    df = pd.read_csv(SOURCES_CSV)
    have = set(df["url"].astype(str))
    new = []
    for e in extra:
        url = e.get("url")
        if not url or url in have:
            continue
        excerpt = str(e.get("excerpt") or "")[:180]
        new.append(dict(
            source="extra probe",
            url=url,
            retrieval_date=RETRIEVAL_DATE,
            found=f"status={e.get('status')} nbytes={e.get('nbytes')} excerpt={excerpt}",
            http_status=e.get("status"),
            nbytes=e.get("nbytes"),
            path=None,
            error=e.get("error"),
            note="follow-up",
        ))
    if new:
        df = pd.concat([df, pd.DataFrame(new)], ignore_index=True)
        df.to_csv(SOURCES_CSV, index=False)
        dump_json(MD_DIR / "sources.json", jsonable(df.to_dict(orient="records")))
        log(f"appended {len(new)} extra sources; n={len(df)}")
    else:
        log(f"no extra sources to append; n={len(df)}")

    man = load_manifest()
    man["status"] = "STAGE1_DONE_GATE_CLOSED"
    man["s1_summary"] = dict(
        n_sources=int(len(df)),
        gate_open=False,
        gene_list_recovered=bool(gate.get("gene_list_recovered")),
        n_MD_genes=gate.get("n_MD_genes"),
        computation_recovered=False,
        gate_reason=gate.get("reason"),
    )
    save_manifest(man)
    p = write_findings()
    log(f"FINDINGS {p}")
    log.close()
    print(p)


if __name__ == "__main__":
    main()

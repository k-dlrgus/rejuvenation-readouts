"""Load Table S3 Age up / Age down from mmc3.xlsx. Exact published lists, not a DE rebuild."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md3_common import (  # noqa: E402
    MD3_DIR, MD2_DIR, MMC3_XLSX, MD3_SEED,
    StopStep, Logger, dump_json, jsonable, load_json, load_manifest, save_manifest,
    log_columns, log_geneset, md3_log_banner, progress_snapshot, record_failure,
    TABLE_S3_REASONING,
)
from md3_idtype import detect_id_type, require_mappable  # noqa: E402
from md_stage1 import _xlsx_columns  # noqa: E402


def _set_compare(a, b):
    a = [str(x).upper() for x in a]
    b = [str(x).upper() for x in b]
    sa, sb = set(a), set(b)
    return dict(
        n_a=len(a), n_b=len(b), n_shared=len(sa & sb),
        only_a=sorted(sa - sb), only_b=sorted(sb - sa),
        equal=sa == sb,
    )


def load_table_s3(log=None):
    close_log = False
    if log is None:
        log = Logger(MD3_DIR / "genesets_report.txt")
        close_log = True
    md3_log_banner(log, "GENESETS")
    log(TABLE_S3_REASONING)
    if not MMC3_XLSX.exists():
        raise StopStep("genesets", f"missing mmc3.xlsx {MMC3_XLSX}")
    sheets = _xlsx_columns(MMC3_XLSX)
    log(f"[mmc3] path={MMC3_XLSX} sheets={list(sheets)}")
    for sname, cols in sheets.items():
        log_columns(f"mmc3_{sname}", list(cols), str(MMC3_XLSX))

    aging = sheets.get("Aging_signatures") or {}
    if "Age up" not in aging:
        raise StopStep("genesets", f"mmc3 Aging_signatures missing 'Age up'. cols={list(aging)}")
    if "Age down" not in aging:
        raise StopStep("genesets", f"mmc3 Aging_signatures missing 'Age down'. cols={list(aging)}")
    age_up = [g.strip().upper() for g in aging["Age up"] if str(g).strip()]
    age_down = [g.strip().upper() for g in aging["Age down"] if str(g).strip()]
    id_up = detect_id_type(age_up, log, "mmc3_Age_up", str(MMC3_XLSX) + ":Aging_signatures:Age up")
    id_dn = detect_id_type(age_down, log, "mmc3_Age_down", str(MMC3_XLSX) + ":Aging_signatures:Age down")
    require_mappable(id_up, {"symbol"}, log, "genesets")
    require_mappable(id_dn, {"symbol"}, log, "genesets")
    log_geneset("mmc3_Aging_signatures_Age up", "Aging_signatures:Age up", age_up, str(MMC3_XLSX),
                id_type=id_up["kind"])
    log_geneset("mmc3_Aging_signatures_Age down", "Aging_signatures:Age down", age_down, str(MMC3_XLSX),
                id_type=id_dn["kind"])
    log(f"[mmc3-aging] Age up n={len(age_up)} Age down n={len(age_down)} "
        f"id_up={id_up['kind']} id_down={id_dn['kind']}")

    md2_gs = MD2_DIR / "genesets.json"
    if not md2_gs.exists():
        raise StopStep("genesets", f"missing md2 genesets {md2_gs}")
    md2 = load_json(md2_gs)
    md = [str(g).upper() for g in md2["MD"]]
    tgfb = [str(g).upper() for g in md2["TGFB"]]
    md2_up = [str(g).upper() for g in (md2.get("mmc3_age_up") or [])]
    md2_dn = [str(g).upper() for g in (md2.get("mmc3_age_down") or [])]
    cmp_up = _set_compare(age_up, md2_up)
    cmp_dn = _set_compare(age_down, md2_dn)
    log(f"[mmc3 vs md2] Age up equal={cmp_up['equal']} n_mmc3={cmp_up['n_a']} n_md2={cmp_up['n_b']} "
        f"only_mmc3={cmp_up['only_a'][:12]} only_md2={cmp_up['only_b'][:12]}")
    log(f"[mmc3 vs md2] Age down equal={cmp_dn['equal']} n_mmc3={cmp_dn['n_a']} n_md2={cmp_dn['n_b']}")
    if not cmp_up["equal"] or not cmp_dn["equal"]:
        raise StopStep(
            "genesets",
            "mmc3 Age up/down do not equal md2 genesets.json copies. Not reconciling. "
            f"up_equal={cmp_up['equal']} down_equal={cmp_dn['equal']}",
        )
    log_geneset("MD_built_from_md2", "MD (STAR Methods built; md2)", md, str(md2_gs), id_type="symbol")
    log_geneset("TGFB_built_from_md2", "TGFB (STAR Methods built; md2)", tgfb, str(md2_gs), id_type="symbol")
    reprog = {k: [str(g).upper() for g in v] for k, v in (md2.get("reprog") or {}).items()}

    summary = dict(
        mmc3_path=str(MMC3_XLSX),
        n_age_up=len(age_up), n_age_down=len(age_down),
        age_up_id_type=id_up["kind"], age_down_id_type=id_dn["kind"],
        n_MD=len(md), n_TGFB=len(tgfb),
        mmc3_equals_md2_age_up=bool(cmp_up["equal"]),
        mmc3_equals_md2_age_down=bool(cmp_dn["equal"]),
        used_as_deseq2_rebuild=False,
        used_as_published_lists=True,
        seed=MD3_SEED,
        columns_read={"Aging_signatures": list(aging)},
    )
    dump_json(MD3_DIR / "genesets_summary.json", jsonable(summary))
    dump_json(MD3_DIR / "genesets.json", jsonable(dict(
        MD=md, TGFB=tgfb,
        age_up=age_up, age_down=age_down,
        reprog=reprog,
        mmc3_path=str(MMC3_XLSX),
        msigdb_version=md2.get("msigdb_version"),
    )))
    man = load_manifest()
    man["status"] = "GENESETS_DONE"
    save_manifest(man)
    progress_snapshot(
        "write PREREG flags if needed, then Task 1 occupancy from md2 labels (no new scores), then score",
        stop="GENESETS_DONE",
    )
    if close_log:
        log.close()
    return summary


if __name__ == "__main__":
    try:
        load_table_s3()
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise

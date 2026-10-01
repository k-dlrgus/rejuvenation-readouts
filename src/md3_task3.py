"""Task 3 — what MD co-varies with. Report-only. No gate."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md3_common import (  # noqa: E402
    MD3_DIR, MD2_DIR, MD3_SEED, MD3_BOOT, N_PERM, N_BOOT, MIN_CELLS_20,
    PREREG_TASK3, PREREG_TASK3_FLAG, STATE_ORDINAL, AGED_LINE, YOUNG_LINE,
    StopStep, Logger, dump_json, jsonable, load_json, md3_log_banner,
    load_manifest, save_manifest, record_failure, log_columns, progress_snapshot,
)
from fibro_common import bootstrap_rho_ci  # noqa: E402
from gtex_common import spearman_safe  # noqa: E402
from trajectory_common import permutation_p  # noqa: E402


def _require():
    if not PREREG_TASK3_FLAG.exists():
        raise StopStep("prereg", "PREREG_TASK3.flag missing")
    p = MD3_DIR / "t1_cluster20_scores.csv"
    if not p.exists():
        raise StopStep("task3", f"missing {p} — run Task 1")
    p2 = MD2_DIR / "t2_crosscorr.csv"
    if not p2.exists():
        raise StopStep("task3", f"missing {p2} — cannot check FINDINGS_MD2.md MD-vs-ruler ρ")


def _rho_one(x, y, log, tag):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = int(len(x))
    rec = dict(tag=tag, n=n, n_perm=N_PERM, n_boot=N_BOOT, seed=MD3_SEED, boot_seed=MD3_BOOT)
    if n < 8:
        rec.update(rho=np.nan, rho_null=np.nan, rho_p=np.nan, rho_ci_lo=np.nan, rho_ci_hi=np.nan,
                   note="n<8")
        return rec
    rho = spearman_safe(x, y)
    rng = np.random.default_rng(MD3_SEED)
    nulls = []
    for _ in range(N_PERM):
        xp = x.copy()
        rng.shuffle(xp)
        nulls.append(spearman_safe(xp, y))
    nulls = np.asarray(nulls, float)
    p = permutation_p(rho, nulls, greater=True)
    rng_b = np.random.default_rng(MD3_BOOT)
    ci = bootstrap_rho_ci(x, y, rng_b, n_boot=N_BOOT)
    rec.update(
        rho=float(rho) if np.isfinite(rho) else np.nan,
        rho_null=float(np.nanmedian(nulls)) if np.isfinite(nulls).any() else np.nan,
        rho_p=float(p) if np.isfinite(p) else np.nan,
        rho_ci_lo=ci.get("p025"), rho_ci_hi=ci.get("p975"),
        note="permutation shuffles MD, holds the other variable; unit=cluster×timepoint ≥20",
    )
    log(f"[t3 {tag}] n={n} ρ={rho:+.3f} null={rec['rho_null']:+.3f} p={p:+.3f} "
        f"CI=[{ci.get('p025')}, {ci.get('p975')}]")
    return rec


def run_task3(log=None):
    _require()
    close_log = False
    if log is None:
        log = Logger(MD3_DIR / "t3_report.txt")
        close_log = True
    md3_log_banner(log, "TASK3")
    log(PREREG_TASK3)
    df = pd.read_csv(MD3_DIR / "t1_cluster20_scores.csv")
    log_columns("t1_cluster20_scores", list(df.columns), str(MD3_DIR / "t1_cluster20_scores.csv"))
    md2cc = pd.read_csv(MD2_DIR / "t2_crosscorr.csv")
    log_columns("md2_t2_crosscorr", list(md2cc.columns), str(MD2_DIR / "t2_crosscorr.csv"))

    rows = []
    checks = []
    for line in (AGED_LINE, YOUNG_LINE):
        g = df[(df.cell_line.astype(str) == line) & (df.below_min_cells == False)].copy()  # noqa: E712
        g = g.dropna(subset=["md_score"])
        log(f"[t3] {line} n_rows≥20 with MD={len(g)}")
        pairs = [
            ("frozen_ruler", "age_score"),
            ("pluri_primary", "pluri_primary"),
            ("tgfb_score", "tgfb_score"),
        ]
        for other_name, col in pairs:
            if col not in g.columns:
                rows.append(dict(cell_line=line, other=other_name, rho=np.nan, reason=f"missing column {col}"))
                continue
            rec = _rho_one(g.md_score.to_numpy(float), g[col].to_numpy(float), log, f"{line}_MD_vs_{other_name}")
            rec.update(cell_line=line, other=other_name, min_cells=MIN_CELLS_20)
            rows.append(rec)
            if other_name == "frozen_ruler":
                hit = md2cc[md2cc.cell_line.astype(str) == line]
                if len(hit):
                    md2_rho = float(hit.rho.iloc[0])
                    delta = rec["rho"] - md2_rho if np.isfinite(rec["rho"]) else np.nan
                    checks.append(dict(
                        cell_line=line, md3_rho=rec["rho"], md2_rho=md2_rho,
                        abs_delta=abs(delta) if np.isfinite(delta) else np.nan,
                        n_md3=rec["n"], n_md2=int(hit.n.iloc[0]) if "n" in hit.columns else None,
                    ))
                    log(f"[t3 check] {line} MD vs ruler md3={rec['rho']:+.3f} md2={md2_rho:+.3f} |Δ|={delta}")
        # ordinal
        g2 = g.copy()
        g2["state_ordinal"] = g2["label"].map(STATE_ORDINAL)
        n_excl = int(g2.state_ordinal.isna().sum())
        g2 = g2.dropna(subset=["state_ordinal"])
        rec = _rho_one(g2.md_score.to_numpy(float), g2.state_ordinal.to_numpy(float),
                       log, f"{line}_MD_vs_state_ordinal")
        rec.update(
            cell_line=line, other="state_ordinal_Fibroblast0_Partial1_EarlyPluri2_Pluri3",
            min_cells=MIN_CELLS_20, n_excluded_NonReprog=n_excl,
            ordinal=STATE_ORDINAL,
            note=(rec.get("note") or "") + f"; NonReprog excluded n={n_excl}",
        )
        rows.append(rec)

    out = pd.DataFrame(rows)
    out.to_csv(MD3_DIR / "t3_md_covariation.csv", index=False)
    chk = pd.DataFrame(checks)
    chk.to_csv(MD3_DIR / "t3_md_vs_ruler_check.csv", index=False)
    summary = dict(
        n_rows=int(len(out)), checks=chk.to_dict(orient="records"),
        ordinal=STATE_ORDINAL, min_cells=MIN_CELLS_20, report_only=True, no_gate=True,
    )
    dump_json(MD3_DIR / "t3_summary.json", jsonable(summary))
    man = load_manifest()
    man["status"] = "TASK3_DONE"
    save_manifest(man)
    progress_snapshot("write FINDINGS_MD3.md from disk", stop="TASK3_DONE")
    if close_log:
        log.close()
    return summary


if __name__ == "__main__":
    try:
        run_task3()
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise

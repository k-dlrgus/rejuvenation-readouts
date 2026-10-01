"""Write FINDINGS_TRAJECTORY.md from results/trajectory/*.json.

Does not modify any prior FINDINGS file or FALSIFICATION.md.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from trajectory_common import (  # noqa: E402
    TRAJ_DIR, TRAJ_SEED, DATASET_ID, STOP_NULL, G2B_RATIO, G2B_DELTA, G2C_RATIO,
    load_json, fmt, fmt_u,
)
from config import ROOT  # noqa: E402


def _md_table(df, cols=None):
    if cols is not None:
        df = df[cols]
    lines = ["| " + " | ".join(df.columns) + " |", "|" + "|".join("---" for _ in df.columns) + "|"]
    for _, row in df.iterrows():
        cells = []
        for c in df.columns:
            v = row[c]
            if isinstance(v, (float, np.floating)):
                cells.append("NA" if not np.isfinite(v) else f"{float(v):+.3f}")
            elif isinstance(v, (bool, np.bool_)):
                cells.append(str(bool(v)))
            elif v is None:
                cells.append("NA")
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def write_findings():
    g1_path = TRAJ_DIR / "g1_summary.json"
    if not g1_path.exists():
        raise FileNotFoundError("results/trajectory/g1_summary.json missing — run G1 first")
    G1 = load_json(g1_path)
    G2 = load_json(TRAJ_DIR / "g2_summary.json") if (TRAJ_DIR / "g2_summary.json").exists() else None
    G3 = load_json(TRAJ_DIR / "g3_summary.json") if (TRAJ_DIR / "g3_summary.json").exists() else None
    search = load_json(TRAJ_DIR / "g3_search_summary.json") if (TRAJ_DIR / "g3_search_summary.json").exists() else None

    stop_g1 = bool(G1.get("stop_G1"))
    stop_g2 = bool(G2.get("stop_G2")) if G2 else False
    skipped_g2 = bool(G2.get("skipped")) if G2 else False

    lines = []
    lines.append("# FINDINGS_TRAJECTORY — geometric scores and the reprogramming trajectory")
    lines.append("")
    g2_status = "not run" if G2 is None else ("skipped" if skipped_g2 else ("STOP" if stop_g2 else "done"))
    g3_status = "not run" if G3 is None else G3.get("status", "done")
    lines.append(
        f"**Status:** G1 {'STOP' if stop_g1 else 'done'}; G2 {g2_status}; G3 {g3_status}. "
        f"Seed `{TRAJ_SEED}`. Dataset ID `{DATASET_ID}` (cached Aging_Cohort h5ad; not invented). "
        f"Does not modify `FINDINGS_GEOMETRY.md`, `FINDINGS_BRAIN_PHASE1.md`, or `FALSIFICATION.md`. "
        f"No gene selection: every gene in the Phase-1 matrix enters every model with a weight."
    )
    lines.append("")
    lines.append(
        "Reproduced by `notebooks/trajectory_g1.ipynb`, `trajectory_g2.ipynb`, "
        "`trajectory_g3.ipynb`, `trajectory_g4.ipynb`. "
        f"Null usability line: site-stratified shuffle R² ≤ {STOP_NULL}."
    )
    lines.append("")

    # Headline / G4 verdict
    verdict, verdict_text = _verdict(G1, G2, G3)
    lines.append("## Headline")
    lines.append("")
    lines.append(_headline(G1, G2, G3, verdict, verdict_text))
    lines.append("")

    lines.append("## G1 — geometric scores")
    lines.append("")
    lines.append(
        f"Cohort: {G1.get('n_donors')} donors, {G1.get('n_types')} cell types, "
        f"{G1.get('n_genes')} genes, {G1.get('n_rows')} pseudobulks. "
        f"Sites: `{G1.get('site_donor_counts')}`. "
        "Age score = projection onto a within-type supervised direction over **all genes** "
        "(PLS-1: covariance direction + 1D OLS to years; ridge: SVD-LOO α, linear predictor). "
        "Identity subspace = SVD of train type centroids; classification = nearest centroid "
        "in that subspace; identity scalar = distance to own type centroid. "
        "Directions and the identity basis are fit on **training folds only**."
    )
    lines.append("")
    lines.append(
        "**Primary age method: ridge** (the actual age clock). PLS-1 is the geometric companion; "
        "both are reported. If they differ by >0.05 in primary R² the difference is called material."
    )
    lines.append("")
    if (TRAJ_DIR / "g1_summary_table.csv").exists():
        tab = pd.read_csv(TRAJ_DIR / "g1_summary_table.csv")
        lines.append(_md_table(tab))
        lines.append("")
    for scheme in ("within_site", "loso"):
        rec = G1.get("schemes", {}).get(scheme, {})
        if not rec:
            continue
        lines.append(f"### {scheme}")
        lines.append("")
        for method in ("pls1", "ridge"):
            m = rec.get(method, {})
            lines.append(
                f"- **age {method}:** primary R²={fmt(m.get('r2'))}, MAE={fmt_u(m.get('mae'), 2)} y; "
                f"pooled R²={fmt(m.get('pooled_r2'))}, site-strat R²={fmt(m.get('site_strat_r2'))}; "
                f"null mean={fmt(m.get('null', {}).get('mean'))} (p95={fmt(m.get('null', {}).get('p95'))}, "
                f"p={fmt_u(m.get('p_r2'))}, n_perm={m.get('n_perm')}, usable={m.get('null_usable')})"
            )
        idm = rec.get("identity", {})
        lines.append(
            f"- **identity nearest-centroid:** acc={fmt_u(idm.get('acc'), 4)}, "
            f"macro-F1={fmt_u(idm.get('macro_f1'), 4)}, chance={fmt_u(idm.get('chance'), 4)}; "
            f"null acc mean={fmt_u(idm.get('null_acc', {}).get('mean'), 4)}, "
            f"p={fmt_u(idm.get('p_acc'))}, n_perm={idm.get('n_perm')}; "
            f"mean dist-to-own-centroid={fmt_u(idm.get('mean_dist_own'), 3)}"
        )
        lines.append("")
    diffs = G1.get("pls_vs_ridge", {})
    if diffs:
        lines.append(
            "PLS-1 vs ridge (absolute primary-R² difference): "
            + "; ".join(
                f"{k} Δ={fmt_u(v.get('abs_delta'))} material={v.get('material')}"
                for k, v in diffs.items()
            )
            + "."
        )
        lines.append("")
    lines.append(
        f"**STOP G1:** {G1.get('stop_G1')} "
        + (f"— {G1.get('stop_reason')}" if G1.get("stop_reason") else "— age direction produces usable held-out R²; null is not inflated.")
    )
    lines.append("")

    lines.append("## G2 — three orthogonality controls")
    lines.append("")
    if G2 is None:
        lines.append("Not run.")
        lines.append("")
    elif G2.get("skipped"):
        lines.append(f"Skipped ({G2.get('reason')}).")
        lines.append("")
    else:
        lines.append(
            "G2a: angle of the age direction to the identity subspace, both fit on the training fold. "
            f"G2b: within-type age R² from the identity coordinates vs from the age direction "
            f"(redundancy if R²_id ≥ {G2B_RATIO} × R²_age or R²_age − R²_id < {G2B_DELTA}). "
            f"G2c: refit the age direction on expression with the identity subspace regressed out "
            f"(collapse if residual R² ≤ {STOP_NULL} or residual < {G2C_RATIO} × R²_age)."
        )
        lines.append("")
        if (TRAJ_DIR / "g2_summary_table.csv").exists():
            tab = pd.read_csv(TRAJ_DIR / "g2_summary_table.csv")
            lines.append(_md_table(tab))
            lines.append("")
        for scheme in ("within_site", "loso"):
            a = G2.get("g2a", {}).get(scheme, {})
            b = G2.get("g2b", {}).get(scheme, {})
            c = G2.get("g2c", {}).get(scheme, {})
            if not a:
                continue
            lines.append(f"### {scheme}")
            lines.append("")
            for method in ("pls1", "ridge"):
                am = a.get(method, {})
                lines.append(
                    f"- **G2a {method}:** angle={fmt_u(am.get('median_angle_deg'), 1)}° "
                    f"(null {fmt_u(am.get('null_angle_mean'), 1)}°, p={fmt_u(am.get('p_angle'))}); "
                    f"proj R²={fmt_u(am.get('median_proj_r2'))} "
                    f"(null {fmt_u(am.get('null_proj_r2_mean'))}, p={fmt_u(am.get('p_proj'))}, "
                    f"n_perm={am.get('n_perm')})"
                )
            lines.append(
                f"- **G2b:** identity-coords age R²={fmt(b.get('r2_identity_coords'))} "
                f"vs ridge age-dir R²={fmt(b.get('age_dir', {}).get('ridge', {}).get('r2'))} "
                f"(ratio={fmt_u(b.get('ratio'))}, Δ={fmt(b.get('delta'))}, "
                f"null={fmt(b.get('null', {}).get('mean'))}, p={fmt_u(b.get('p'))}); "
                f"redundant={b.get('redundant')}"
            )
            rr = c.get("residual", {}).get("ridge", {})
            rp = c.get("residual", {}).get("pls1", {})
            lines.append(
                f"- **G2c ridge:** residual age R²={fmt(rr.get('r2'))} "
                f"(null {fmt(rr.get('null', {}).get('mean'))}, p={fmt_u(rr.get('p'))}); "
                f"PLS-1 residual R²={fmt(rp.get('r2'))}; collapsed={c.get('collapsed')}"
            )
            lines.append("")
        lines.append(
            f"**STOP G2:** {G2.get('stop_G2')} "
            + (f"— {G2.get('stop_reason')}" if G2.get("stop_reason") else
               "— identity coordinates do not predict age as well as the age direction; residual age R² survives.")
        )
        lines.append("")

    lines.append("## G3 — trajectory")
    lines.append("")
    if search:
        lines.append(
            f"Search seed `{TRAJ_SEED}`. Candidate accessions re-fetched: "
            f"{search.get('candidates')}. Found: `{search.get('candidate_found')}`. "
            f"GEO hits={search.get('n_geo_hits')}; CXG mouse-brain RNA datasets={search.get('n_cxg_mouse_brain')}; "
            f"CXG human OSK title matches={search.get('n_cxg_human_osk')}."
        )
        lines.append("")
    if G3 is None:
        lines.append("Trajectory not computed yet.")
        lines.append("")
    else:
        lines.append(G3.get("narrative") or G3.get("report") or str(G3.get("status")))
        lines.append("")
        if G3.get("orthologs"):
            o = G3["orthologs"]
            lines.append(
                f"Ortholog mapping: human genes={o.get('n_human')}, "
                f"mapped 1:1={o.get('n_1to1')}, used={o.get('n_used')} "
                f"({fmt_u(o.get('frac_used'))} of the frozen human genes)."
            )
            lines.append("")
        if G3.get("mouse_aging"):
            ma = G3["mouse_aging"]
            lines.append(
                f"Independent mouse aging validation: dataset `{ma.get('id')}` "
                f"({ma.get('title', '')[:160]}). "
                f"OLS R² of mouse age ~ frozen score = {fmt(ma.get('r2'))} "
                f"(null {fmt(ma.get('null_mean'))}, p={fmt_u(ma.get('p'))}); "
                f"Spearman ρ={fmt(ma.get('spearman'))} (p={fmt_u(ma.get('p_spearman'))}). "
                f"Per-type Spearman signs +{ma.get('n_types_spearman_pos')} / -{ma.get('n_types_spearman_neg')}. "
                f"Transfer usable={ma.get('transfer_ok')}."
            )
            cx = ma.get("cortex_sensitivity") or {}
            if cx:
                lines.append(
                    f" Cortex-subtissue sensitivity: n_donors={cx.get('n_donors')}, "
                    f"OLS R²={fmt(cx.get('r2'))}, Spearman={fmt(cx.get('spearman'))}."
                )
            lines.append("")
            if (TRAJ_DIR / "g3_tms_per_type_transfer.csv").exists():
                lines.append("Per-type transfer (Tabula Muris Senis, frozen scores):")
                lines.append("")
                lines.append(_md_table(pd.read_csv(TRAJ_DIR / "g3_tms_per_type_transfer.csv")))
                lines.append("")
        if G3.get("needed") or G3.get("needed_if_d"):
            lines.append(f"**Dataset that would be needed:** {G3.get('needed') or G3.get('needed_if_d')}")
            lines.append("")
        if G3.get("trajectory"):
            t = G3["trajectory"]
            lines.append(t.get("text", ""))
            lines.append("")
            if t.get("window"):
                w = t["window"]
                lines.append(
                    f"Candidate window: {w.get('range')}; age reduction={fmt(w.get('d_age'))}; "
                    f"identity retained={fmt_u(w.get('id_retained'))} "
                    f"(CIs: age {w.get('age_ci')}, identity {w.get('id_ci')})."
                )
                lines.append("")
        if G3.get("figure"):
            lines.append(f"Figure: `{G3['figure']}`.")
            lines.append("")

    lines.append("## G4 — verdict")
    lines.append("")
    lines.append(f"**({verdict})** {verdict_text}")
    lines.append("")
    lines.append("The four options were: (a) independent scores and a bending trajectory (safe window measured); "
                 "(b) independent scores, trajectory does not bend (axes separable in aging, fused under reprogramming); "
                 "(c) scores not independent (G2 failed); "
                 "(d) scores independent but do not transfer to available reprogramming data.")
    lines.append("")

    lines.append("## Limitations")
    lines.append("")
    lines.append("1. Two brain banks only (HBCC, MSSM). LOSO is one df of transfer and is reported, not averaged away.")
    lines.append("2. Unmeasured 6-plex hashing pools, PMI, RIN — as in FINDINGS_BRAIN_PHASE1.md.")
    lines.append("3. Cross-species transfer: scores are fit on human DLPFC; mouse application requires ortholog mapping and loses unmapped genes. Transfer is not assumed — it is tested on an independent mouse aging dataset before any reprogramming plot.")
    lines.append(f"4. Seed `{TRAJ_SEED}` (`numpy.random.default_rng`). Folds match Phase 1 / geometry.")
    lines.append("5. No gene selection. p ≫ n; ridge uses SVD-LOO α. PLS-1 is the covariance direction.")
    lines.append("6. Identity classification from a 19-dimensional centroid subspace of 20 types is expected to be near-perfect; that is a geometric fact, not evidence that identity is 'solved'.")
    lines.append("7. Per-type age directions were largely type-specific in Part C (pairwise ~78°). A single frozen consensus direction is not assumed; G3 applies per-type weights after mapping mouse types.")
    lines.append("8. Cross-species orthologs: NCBI strict 1:1 Ensembl pairs cover 14,239 / 25,526 Phase-1 genes (55.8%). Unmapped genes sit at the human mean (z = 0) in the frozen scores.")
    lines.append("9. GSE224438 (mouse SVZ partial reprogramming) is the usable reprogramming accession found; whole-body 10x libraries multiplex young/old/OSK in one lane and RAW has no hashing barcodes. Those libraries, and the SVZ-targeted lanes, were **not scored** because the frozen age score failed the independent mouse aging transfer test.")
    lines.append("10. GSE224438 is a 3-condition protocol contrast (young / old / old+OSK), not a multi-day time course. Even if transfer had passed, quadratic vs linear would be unidentified with three points.")
    lines.append("11. GSE276656 (mPFC engram OSK multiome, 14.4 GB RAW) was not downloaded. GSE271794 is E15.5 development, not adult aging. No human brain OSK scRNA dataset was found on CELLxGENE (0 title matches).")
    lines.append("12. Transfer is tested on Tabula Muris Senis brain (CXG `66ff82b4-…` + `c08f8441-…`) before any reprogramming plot. That test failed (OLS R²=+0.130, p=0.235; cortex subset Spearman=−0.254); no trajectory is reported.")
    lines.append("")

    lines.append("## Files")
    lines.append("")
    lines.append("| path | content |")
    lines.append("|---|---|")
    lines.append("| `src/trajectory_common.py`, `trajectory_g1.py`, `trajectory_g2.py`, `trajectory_g3_search.py`, `trajectory_g3.py`, `trajectory_findings.py` | code |")
    lines.append("| `notebooks/trajectory_g1.ipynb`, `g2`, `g3`, `g4` | runnable from a clean checkout |")
    lines.append("| `results/trajectory/` | tables, json, logs, `figures/` |")
    lines.append("| `FINDINGS_TRAJECTORY.md` | this file |")
    lines.append("")

    path = ROOT / "FINDINGS_TRAJECTORY.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _headline(G1, G2, G3, verdict, verdict_text):
    bits = []
    rec = G1.get("schemes", {}).get("within_site", {})
    ridge = rec.get("ridge", {})
    pls = rec.get("pls1", {})
    bits.append(
        f"**G1 age (within-site, site-stratified):** ridge R²={fmt(ridge.get('r2'))} "
        f"(null {fmt(ridge.get('null', {}).get('mean'))}); "
        f"PLS-1 R²={fmt(pls.get('r2'))} (null {fmt(pls.get('null', {}).get('mean'))})."
    )
    idm = rec.get("identity", {})
    bits.append(
        f"**G1 identity (subspace nearest centroid):** acc={fmt_u(idm.get('acc'), 4)} "
        f"(chance {fmt_u(idm.get('chance'), 4)}; shuffle {fmt_u(idm.get('null_acc', {}).get('mean'), 4)})."
    )
    if G2 and not G2.get("skipped"):
        b = G2.get("g2b", {}).get("within_site", {})
        c = G2.get("g2c", {}).get("within_site", {})
        a = G2.get("g2a", {}).get("within_site", {}).get("pls1", {})
        bits.append(
            f"**G2a (within-site PLS-1):** angle={fmt_u(a.get('median_angle_deg'), 1)}° "
            f"(null {fmt_u(a.get('null_angle_mean'), 1)}°, p={fmt_u(a.get('p_angle'))})."
        )
        bits.append(
            f"**G2b:** identity-coords R²={fmt(b.get('r2_identity_coords'))} vs age-dir "
            f"{fmt(b.get('age_dir', {}).get('ridge', {}).get('r2'))}; redundant={b.get('redundant')}."
        )
        bits.append(
            f"**G2c:** residual ridge R²={fmt(c.get('residual', {}).get('ridge', {}).get('r2'))}; "
            f"collapsed={c.get('collapsed')}."
        )
    bits.append(f"**G4 verdict: ({verdict})** {verdict_text}")
    return " ".join(bits)


def _verdict(G1, G2, G3):
    if G1.get("stop_G1"):
        return "c", (
            "The direction method fails at G1 (null unusable or age R² collapsed under both CV schemes). "
            "The project's central construct does not hold on this cohort; report as a negative result."
        )
    if G2 is None:
        return "pending", "G2 has not been run."
    if G2.get("skipped"):
        return "c", f"G2 skipped because G1 stopped ({G2.get('reason')})."
    if G2.get("stop_G2"):
        return "c", (
            "The scores are not independent (G2 failed): "
            f"{G2.get('stop_reason')}. The project's central construct does not hold; report as a negative result."
        )
    if G3 is None:
        return "pending", "Scores are independent on G2; G3 has not been run."
    v = G3.get("verdict")
    if v in ("a", "b", "c", "d"):
        return v, G3.get("verdict_text") or G3.get("narrative") or ""
    if G3.get("status") == "no_transfer":
        return "d", G3.get("verdict_text") or (
            "Scores are independent in human DLPFC but do not transfer to the available reprogramming data."
        )
    if G3.get("status") == "no_bend":
        return "b", G3.get("verdict_text") or (
            "Two independent scores exist, but the trajectory does not bend."
        )
    if G3.get("status") == "window":
        return "a", G3.get("verdict_text") or (
            "Two independent scores exist, and the trajectory bends."
        )
    return "pending", str(G3.get("status"))


if __name__ == "__main__":
    p = write_findings()
    print("wrote", p)

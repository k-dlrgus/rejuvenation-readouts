"""Write FINDINGS_TARGET.md from results/target/*.json.

Does not modify any prior FINDINGS file or FALSIFICATION.md.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from target_common import (  # noqa: E402
    TARGET_DIR, TARGET_SEED, DATASET_ID, STOP_NULL, STOP_BOOT_ANGLE,
    V3_RATIO_NEAR, V3_LOSS_NEAR, load_json, fmt, fmt_u, md_table,
)
from config import ROOT  # noqa: E402


def _site_lab(s):
    s = str(s)
    return {"H": "HBCC", "M": "MSSM"}.get(s, s)


def _junk_symbol(sym):
    import re
    s = str(sym)
    if s.startswith("ENSG") or s.startswith("LINC") or s.startswith("OR"):
        return True
    if "-AS" in s or s.endswith("-DT") or re.search(r"P\d+$", s):
        return True
    return False


def _cv_row(rec, scheme, key):
    r = rec[scheme][key]
    return dict(
        scheme=scheme, model=key,
        r2=r["primary_r2"], mae=r["mae"],
        pooled_r2=r["pooled_r2"], site_strat_r2=r["site_strat_r2"],
        null_mean=r["null_mean"], null_p95=r["null_p95"], p=r["p"],
        null_usable=r["null_usable"], n_perm=r["n_perm"],
        types_lose=f"{r['n_types_lose']}/{r['n_types']}",
    )


def write_findings():
    v1_path = TARGET_DIR / "v1_summary.json"
    if not v1_path.exists():
        raise FileNotFoundError("results/target/v1_summary.json missing — run V1 first")
    V1 = load_json(v1_path)
    V2 = load_json(TARGET_DIR / "v2_summary.json") if (TARGET_DIR / "v2_summary.json").exists() else None
    V3 = load_json(TARGET_DIR / "v3_summary.json") if (TARGET_DIR / "v3_summary.json").exists() else None
    V4 = load_json(TARGET_DIR / "v4_summary.json") if (TARGET_DIR / "v4_summary.json").exists() else None

    stop_v1 = bool(V1.get("stop_V1"))
    skip_v2 = bool(V2 and V2.get("skipped"))
    stop_v2 = bool(V2 and V2.get("stop_V2"))
    skip_v3 = bool(V3 and V3.get("skipped"))
    skip_v4 = bool(V4 and V4.get("skipped"))

    def status():
        s1 = "STOP" if stop_v1 else "done"
        s2 = "not run" if V2 is None else ("skipped" if skip_v2 else ("STOP" if stop_v2 else "done"))
        s3 = "not run" if V3 is None else ("skipped" if skip_v3 else "done")
        s4 = "not run" if V4 is None else ("skipped" if skip_v4 else "done")
        return s1, s2, s3, s4

    s1, s2, s3, s4 = status()
    lines = []
    lines.append("# FINDINGS_TARGET — the target vector (expression space only)")
    lines.append("")
    lines.append(
        f"**Status:** V1 {s1}; V2 {s2}; V3 {s3}; V4 {s4}. "
        f"Seed `{TARGET_SEED}`. Dataset ID `{DATASET_ID}` (cached Aging_Cohort h5ad; not invented). "
        f"Does not modify `FINDINGS_TRAJECTORY.md`, `FINDINGS_GEOMETRY.md`, "
        f"`FINDINGS_BRAIN_PHASE1.md`, or `FALSIFICATION.md`. "
        f"No gene selection: every gene in the Phase-1 matrix enters every model with a weight. "
        f"No perturbation data. No TF Atlas. Expression space only."
    )
    lines.append("")
    lines.append(
        "Reproduced by `notebooks/target_v1.ipynb`, `target_v2.ipynb`, `target_v3.ipynb`, "
        f"`target_v4.ipynb`, `target_v5.ipynb`. Null usability line: site-stratified shuffle R² ≤ {STOP_NULL}."
    )
    lines.append("")
    lines.append("## Headline")
    lines.append("")
    lines.append(_headline(V1, V2, V3, V4, stop_v1, stop_v2))
    lines.append("")

    # ---- V1
    lines.append("## V1 — recover the axes and bootstrap stability")
    lines.append("")
    lines.append(
        f"Cohort: {V1['n_donors']} donors, {V1['n_types']} cell types, {V1['n_genes']} genes. "
        f"Identity subspace: {V1['n_identity_axes']} axes from type-centroid SVD. "
        f"Working space: global z-score; ridge coefficients live in that space. "
        f"Donor bootstrap with replacement (ridge n_boot={V1['n_boot_ridge']}, "
        f"PLS-1 n_boot={V1['n_boot_pls']}). Angles are unsigned."
    )
    lines.append("")
    lines.append(
        f"**Ridge median-over-types pairwise bootstrap angle:** {V1['ridge_median_pairwise_deg']:.1f}° "
        f"(to full-data direction {V1.get('ridge_median_to_full_deg', float('nan')):.1f}°). "
        f"**PLS-1:** {V1['pls_median_pairwise_deg']:.1f}°. "
        f"Types above {STOP_BOOT_ANGLE:.0f}°: {V1['n_types_unstable_ridge']}/{V1['n_types_ok_ridge']}."
    )
    if V1.get("frozen_median_angle_deg") is not None:
        lines.append(
            f"Median angle to G1 frozen ridge directions: {V1['frozen_median_angle_deg']:.2f}° "
            "(spaces differ by G1's second per-type z-score)."
        )
    lines.append("")
    boot = pd.DataFrame(V1["per_type_ridge"])
    cols = [c for c in ("celltype", "n_boot_ok", "median_pairwise_deg", "p05_pairwise_deg",
                        "p95_pairwise_deg", "median_to_full_deg") if c in boot.columns]
    lines.append(md_table(boot, cols))
    lines.append("")
    lines.append(
        f"**STOP V1:** {stop_v1} — "
        + ("median bootstrap angle between ridge replicates > ~60°; the target would be noise."
           if stop_v1 else
           "the ridge age direction is stable enough under donor resampling to define a target.")
    )
    lines.append("")

    if V2 is None or skip_v2:
        lines.append("## V2 — skipped")
        lines.append("")
        lines.append("V2 was not run because V1 halted." if stop_v1 else "V2 has not been run.")
        lines.append("")
        _limitations(lines, V1, V2, V3, V4)
        _files(lines)
        path = ROOT / "FINDINGS_TARGET.md"
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return path

    # ---- V2
    lines.append("## V2 — the target vector")
    lines.append("")
    lines.append(
        "TARGET = the age direction with its identity-subspace component removed, then renormalized. "
        "Positive TARGET is the older direction; a rejuvenating move would be −TARGET. "
        "Held-out scores are 1D OLS of age on the projection, fit on the training fold only."
    )
    lines.append("")
    lines.append(
        f"Full-cohort ridge: median norm retained = {V2['ridge_median_norm_retained']:.3f}; "
        f"median ∠(TARGET, age) = {V2['ridge_median_angle_to_age_deg']:.2f}°; "
        f"median ∠(age, identity) = {V2['ridge_median_angle_to_id_deg']:.1f}°. "
        f"PLS-1 companion: retained {V2['pls_median_norm_retained']:.3f}, "
        f"∠(TARGET, age) {V2['pls_median_angle_to_age_deg']:.2f}°."
    )
    lines.append("")
    geom_path = TARGET_DIR / "v2_geom_full.csv"
    if geom_path.exists():
        g = pd.read_csv(geom_path)
        g = g[g.method == "ridge"]
        lines.append(md_table(g, [c for c in ("celltype", "n", "norm_retained", "angle_to_age_deg",
                                              "angle_to_id_deg", "proj_r2_id", "alpha") if c in g.columns]))
        lines.append("")

    cv_rows = []
    for sch in ("within_site", "loso"):
        if sch not in V2.get("cv", {}):
            continue
        for key in ("ridge_target", "ridge_age", "pls1_target"):
            if key in V2["cv"][sch]:
                cv_rows.append(_cv_row(V2["cv"], sch, key))
    if cv_rows:
        lines.append(md_table(pd.DataFrame(cv_rows)))
        lines.append("")
        ws = V2["cv"]["within_site"]["ridge_target"]
        lines.append(
            f"**Within-site ridge TARGET:** primary R²={fmt(ws['primary_r2'])} "
            f"(MAE {ws['mae']:.2f} y; null {fmt(ws['null_mean'])}, p={ws['p']:.3f}). "
            f"Raw ridge age companion R²={fmt(V2['cv']['within_site']['ridge_age']['primary_r2'])}. "
            f"Types at or below null: {ws['n_types_lose']}/{ws['n_types']}. "
            "G1's within-site ridge R²=+0.317 used a second per-type z-score and the un-normalized "
            "ridge predictor; this number is 1D OLS in the global-z space that TARGET and P_id share. "
            "TARGET ≈ the age direction (6°), so the gap vs G1 is the metric, not identity projection."
        )
        if "loso" in V2["cv"]:
            lo = V2["cv"]["loso"]["ridge_target"]
            lines.append(
                f"**LOSO ridge TARGET:** primary R²={fmt(lo['primary_r2'])} "
                f"(null {fmt(lo['null_mean'])}, p={lo['p']:.3f})."
            )
        lines.append("")
    lines.append(
        f"**STOP V2:** {stop_v2} — "
        + ("TARGET loses held-out age predictivity vs null in a majority of cell types; "
           "the age signal would live entirely in the identity subspace, contradicting G2c."
           if stop_v2 else
           "TARGET retains held-out age predictivity in a majority of cell types "
           "(consistent with G2c: residual age R² survived identity residualization).")
    )
    lines.append("")

    if V3 is None or skip_v3:
        lines.append("## V3 — skipped")
        lines.append("")
        lines.append("V3 was not run because V2 halted." if stop_v2 else "V3 has not been run.")
        lines.append("")
        _v5(lines, V1, V2, None, None, stop_v1, stop_v2)
        _limitations(lines, V1, V2, V3, V4)
        _files(lines)
        path = ROOT / "FINDINGS_TARGET.md"
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return path

    # ---- V3
    lines.append("## V3 — one target, or twenty?")
    lines.append("")
    lines.append(
        f"**V3a.** Pairwise unsigned angle among per-type ridge TARGETs: median "
        f"{V3['pairwise_median_deg']:.1f}° (permutation null {V3['pairwise_null_mean']:.1f}°, "
        f"p={V3['p_pairwise']:.3f}, n_perm={V3['n_perm_pair']}). "
        f"Only ~{90 - V3['pairwise_median_deg']:.0f}° from orthogonal. "
        f"The pairwise angle is **not below the permutation null** — the 20 TARGETs are as "
        f"orthogonal as chance. (Part C's PLS-1 age directions were 78° / PC1 0.40; ridge TARGET "
        f"is more type-specific still.)"
    )
    lines.append("")
    lines.append(
        f"**V3b.** PC1 of the 20 (ok={V3['n_types_ok']}) per-type TARGETs explains "
        f"{V3['shared_pc1_frac']:.3f} of their variance (null {V3['pc1_null_mean']:.3f}, "
        f"p={V3['p_pc1']:.3f}). That is a statistically detectable but tiny excess — not a shared "
        f"axis. Jointly-fit within-type ridge TARGET: median angle to per-type TARGETs "
        f"{V3['joint_median_angle_to_types']:.1f}° "
        f"(norm retained after identity projection {V3['joint_norm_retained']:.3f})."
    )
    lines.append("")
    ws = V3["v3c"]["within_site"]
    lo = V3["v3c"].get("loso")
    lines.append(
        f"**V3c within-site.** Own TARGET R²={fmt(ws['own_primary_r2'])}; "
        f"consensus R²={fmt(ws['cons_primary_r2'])} (null {fmt(ws['cons_null_mean'])}, "
        f"p={ws['cons_p']:.3f}); median per-type loss={fmt(ws['median_loss'])}; "
        f"median R²_cons/R²_own={fmt(ws['median_ratio'])}."
    )
    if lo:
        lines.append(
            f"**V3c LOSO.** Own R²={fmt(lo['own_primary_r2'])}; "
            f"consensus R²={fmt(lo['cons_primary_r2'])} "
            f"(null {fmt(lo['cons_null_mean'])}, p={lo['cons_p']:.3f}); "
            f"median loss={fmt(lo['median_loss'])}."
        )
    lines.append("")
    pt_path = TARGET_DIR / "v3c_within_site_per_type.csv"
    if pt_path.exists():
        ptab = pd.read_csv(pt_path)
        lines.append(md_table(ptab))
        lines.append("")
        glia = ("astrocyte", "oligodendrocyte", "oligodendrocyte precursor cell", "microglial cell")
        g = ptab[ptab.celltype.isin(glia)]
        if len(g):
            lines.append(
                "Glial types lose the consensus (astrocyte / oligodendrocyte / OPC / microglia "
                f"own R² med={g.r2_own.median():+.3f} → consensus {g.r2_consensus.median():+.3f}). "
                "The weak shared component is a neuronal compromise, not a cortex-wide axis."
            )
            lines.append("")
    lines.append(
        f"**V3 verdict: ({V3['verdict']})** {V3['verdict_text']}"
    )
    lines.append("")

    # ---- V4
    if V4 is None or skip_v4:
        lines.append("## V4 — skipped")
        lines.append("")
        lines.append("V4 has not been run.")
        lines.append("")
    else:
        lines.append("## V4 — how trustworthy is it?")
        lines.append("")
        lines.append("### V4a. Site transfer")
        lines.append("")
        lines.append(
            "TARGET fit on one bank, 1D OLS calibrated on that bank, scored on the other "
            "(H=HBCC, M=MSSM). This is the known weak point of the age clock "
            "(G1 within-site 0.317 → LOSO 0.101). Median-over-types R² is **negative** in both "
            "directions: the vector does not predict chronological age in the other bank. "
            "It still beats a severely negative permutation null (p=0.048), which only says it is "
            "not pure overfit noise — it does not make it a usable transfer clock. "
            "L2/3 IT is the exception (positive R² both ways); glia and most interneurons fail."
        )
        lines.append("")
        v4a = pd.DataFrame(V4["v4a"]).copy()
        v4a["train_site"] = v4a.train_site.map(_site_lab)
        v4a["test_site"] = v4a.test_site.map(_site_lab)
        lines.append(md_table(v4a))
        lines.append("")
        pt4 = TARGET_DIR / "v4a_per_type.csv"
        if pt4.exists():
            t = pd.read_csv(pt4)
            pos = t[t.r2 > 0][["train_site", "test_site", "celltype", "r2", "mae", "n_tr", "n_te"]].copy()
            if len(pos):
                pos["train_site"] = pos.train_site.map(_site_lab)
                pos["test_site"] = pos.test_site.map(_site_lab)
                lines.append("Types with R²>0 on transfer (all other type×direction cells are in `v4a_per_type.csv`):")
                lines.append("")
                lines.append(md_table(pos))
                lines.append("")

        lines.append("### V4b. Age-range validity")
        lines.append("")
        b = V4["v4b"]
        lines.append(
            f"Split: young age<50 (n_donors={b['n_donors_young']}) vs old age≥50 "
            f"(n_donors={b['n_donors_old']}). Median per-type angle between the two TARGETs: "
            f"**{b['median_angle_deg']:.1f}°**. Consensus-to-consensus angle: "
            f"{b.get('consensus_angle_deg') if b.get('consensus_angle_deg') is not None else float('nan'):.1f}°. "
            f"Young→old transfer R²={fmt(b['young_to_old_median_r2'])}; "
            f"old→young R²={fmt(b['old_to_young_median_r2'])}."
        )
        lines.append(
            "This split is **confounded with bank**: young donors are 75 HBCC / 27 MSSM; "
            "old donors are 45 HBCC / 86 MSSM. The 86° young–old angle therefore mixes "
            "lifespan non-stationarity with the site-transfer failure in V4a. Either way, "
            "the direction is not a single line that can be extrapolated across the adult span."
        )
        if b["median_angle_deg"] > 60:
            lines.append(
                "The two half-lifespan directions differ sharply (median > 60°): the target is "
                "**not one line across the adult lifespan** and cannot be extrapolated from young to old "
                "or from old to young as if it were a single axis."
            )
        else:
            lines.append(
                "The two half-lifespan directions are closer than the 60° instability bar; "
                "a single line is not ruled out, but the angle is still reported and is not assumed to be zero."
            )
        lines.append("")
        bpt = TARGET_DIR / "v4b_per_type_angles.csv"
        if bpt.exists():
            lines.append(md_table(pd.read_csv(bpt)))
            lines.append("")

        lines.append("### V4c. Biological sanity (interpretation only — not gene selection)")
        lines.append("")
        c = V4["v4c"]
        tech = c.get("tech") or {}
        def _tech_line(name, d):
            if not d:
                return f"- {name}: NA"
            return (f"- {name}: n={d.get('n')}  mito={d.get('n_mito')}  "
                    f"ribo={d.get('n_ribo')}  cell-cycle={d.get('n_cell_cycle')}")
        lines.append("Technical composition of the top-50 weight lists:")
        lines.append(_tech_line("consensus +", tech.get("consensus_pos")))
        lines.append(_tech_line("consensus −", tech.get("consensus_neg")))
        lines.append(_tech_line("mean-per-type +", tech.get("mean_pos")))
        lines.append(_tech_line("mean-per-type −", tech.get("mean_neg")))
        lines.append("")
        lines.append(
            f"Enrichr terms returned: {c.get('n_enrichr', 0)}; "
            f"aging-token terms with p<0.05: {c.get('n_aging_token_p05', 0)}."
        )
        enr = TARGET_DIR / "v4c_enrichr.csv"
        if enr.exists():
            e = pd.read_csv(enr)
            cons_e = e[e.list.str.startswith("consensus")].sort_values("p").head(10)
            aging_e = e[e.aging_token & (e.p < 0.05)].sort_values("p").head(15)
            keep = [c0 for c0 in ("list", "library", "term", "p", "adj_p", "aging_token") if c0 in e.columns]
            if len(cons_e):
                lines.append("")
                lines.append("Enrichr, consensus lists (top 10 by p):")
                lines.append("")
                lines.append(md_table(cons_e, keep))
            if len(aging_e):
                lines.append("")
                lines.append("Aging-token terms with p<0.05 (top 15; mostly per-type, not consensus):")
                lines.append("")
                lines.append(md_table(aging_e, keep))
        topg = TARGET_DIR / "v4c_top_genes.csv"
        if topg.exists():
            tg = pd.read_csv(topg)
            cons = tg[tg.source == "consensus"]
            pos = cons[cons.sign == "+"].symbol.head(15).astype(str).tolist()
            neg = cons[cons.sign == "-"].symbol.head(15).astype(str).tolist()
            lines.append("")
            lines.append(f"Consensus top + weights (15 of 50): {', '.join(pos)}")
            lines.append(f"Consensus top − weights (15 of 50): {', '.join(neg)}")
        lines.append("")
        lines.append(_v4c_reading(V4))
        lines.append("")

    _v5(lines, V1, V2, V3, V4, stop_v1, stop_v2)
    _limitations(lines, V1, V2, V3, V4)
    _files(lines)
    path = ROOT / "FINDINGS_TARGET.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _headline(V1, V2, V3, V4, stop_v1, stop_v2):
    bits = [
        f"**V1 ridge bootstrap:** median pairwise angle {V1['ridge_median_pairwise_deg']:.1f}° "
        f"(STOP={stop_v1}; bar {STOP_BOOT_ANGLE:.0f}°)."
    ]
    if V2 and not V2.get("skipped"):
        ws = V2["cv"]["within_site"]["ridge_target"]
        bits.append(
            f"**V2 ridge TARGET** retained {V2['ridge_median_norm_retained']:.3f} of the age direction "
            f"(∠ {V2['ridge_median_angle_to_age_deg']:.1f}°); within-site R²={fmt(ws['primary_r2'])} "
            f"(null {fmt(ws['null_mean'])}); types≤null {ws['n_types_lose']}/{ws['n_types']}; STOP V2={stop_v2}."
        )
    if V3 and not V3.get("skipped"):
        ws = V3["v3c"]["within_site"]
        bits.append(
            f"**V3 ({V3['verdict']}):** pairwise TARGET {V3['pairwise_median_deg']:.1f}° "
            f"(null {V3['pairwise_null_mean']:.1f}°); PC1 {V3['shared_pc1_frac']:.3f} "
            f"(null {V3['pc1_null_mean']:.3f}); consensus R²={fmt(ws['cons_primary_r2'])} vs own "
            f"{fmt(ws['own_primary_r2'])} (loss {fmt(ws['median_loss'])}). {V3['verdict_text']}"
        )
    if V4 and not V4.get("skipped"):
        a = V4["v4a"]
        trans = "; ".join(
            f"{_site_lab(r['train_site'])}→{_site_lab(r['test_site'])} R²={fmt(r['median_r2'])} "
            f"(null {fmt(r['null_mean'])})"
            for r in a
        )
        b = V4["v4b"]
        bits.append(
            f"**V4a site transfer:** {trans}. **V4b age-range:** young vs old TARGET median angle "
            f"{b['median_angle_deg']:.1f}°."
        )
    return " ".join(bits)


def _v4c_reading(V4):
    topg = TARGET_DIR / "v4c_top_genes.csv"
    n_junk = n_cons = 0
    if topg.exists():
        tg = pd.read_csv(topg)
        cons = tg[tg.source == "consensus"]
        n_cons = int(len(cons))
        n_junk = int(cons.symbol.astype(str).map(_junk_symbol).sum()) if n_cons else 0
    enr = TARGET_DIR / "v4c_enrichr.csv"
    n_micro_infl = 0
    n_cons_aging = 0
    if enr.exists():
        e = pd.read_csv(enr)
        n_micro_infl = int(((e.list.str.contains("microglial")) & e.aging_token & (e.p < 0.05)).sum())
        n_cons_aging = int((e.list.str.startswith("consensus") & e.aging_token & (e.p < 0.05)).sum())
    bits = []
    if n_cons and n_junk >= 0.3 * n_cons:
        bits.append(
            f"The consensus top-weight tails are dominated by sparse/p≫n artifacts "
            f"({n_junk}/{n_cons} symbols are LINCs, pseudogenes, olfactory receptors, or Ensembl IDs). "
            "That is not a readable aging program."
        )
    if n_micro_infl:
        bits.append(
            f"The **microglial** TARGET is the exception: {n_micro_infl} aging-token Enrichr hits "
            "(inflammatory response / inflammasome). That is type-specific biology, which is the V3 point."
        )
    if n_cons_aging == 0:
        bits.append(
            "Consensus Enrichr does not recover a clean inflammation / proteostasis / mitochondrial / "
            "senescence signature. DNA-repair terms appear but sit on a junk-heavy gene list."
        )
    bits.append(
        "Do not select genes from these tails. A perturbation mapped onto TARGET inherits a "
        "correlational linear direction, not a CRISPR menu."
    )
    return "**V4c reading:** " + " ".join(bits)


def _v5(lines, V1, V2, V3, V4, stop_v1, stop_v2):
    lines.append("## V5 — verdict")
    lines.append("")
    usable = (not stop_v1) and (V2 is not None) and (not V2.get("skipped")) and (not stop_v2)
    if not usable:
        lines.append(
            "**There is not a usable target vector** at this step: "
            + ("the age direction is unstable under donor bootstrap (V1)." if stop_v1
               else "TARGET does not retain held-out age predictivity (V2).")
        )
        lines.append("")
        lines.append(_assumptions(V1, V2, V3, V4))
        lines.append("")
        return

    one_or_many = "unknown"
    if V3 and not V3.get("skipped"):
        one_or_many = {"a": "one shared target", "b": "twenty type-specific targets",
                       "partial": "partially shared / not a substitute for per-type targets"}.get(
            V3["verdict"], V3["verdict"])
        lines.append(f"**Is there a usable target vector?** Within one bank, within one cell type, "
                     f"as a correlational linear score: yes (V1 stable, V2 predicts age). "
                     f"**One or many?** **{one_or_many}.** "
                     f"It is not a portable rejuvenation axis (V4a negative site transfer; V4b 86° young vs old).")
        lines.append("")
        lines.append(V3["verdict_text"])
    else:
        lines.append("**Is there a usable target vector?** V2 says the per-type TARGET predicts age. "
                     "V3 (one vs twenty) was not completed.")
    lines.append("")

    trust = []
    if V2 and "cv" in V2 and "loso" in V2["cv"]:
        ws = V2["cv"]["within_site"]["ridge_target"]
        lo = V2["cv"]["loso"]["ridge_target"]
        trust.append(
            f"Within-site R²={fmt(ws['primary_r2'])} vs LOSO R²={fmt(lo['primary_r2'])} "
            "(the known site-transfer drop)."
        )
    if V4 and not V4.get("skipped"):
        for r in V4["v4a"]:
            trust.append(
                f"Fit {_site_lab(r['train_site'])} → test {_site_lab(r['test_site'])}: "
                f"R²={fmt(r['median_r2'])} (null {fmt(r['null_mean'])}, p={r['p']:.3f})."
            )
        b = V4["v4b"]
        trust.append(
            f"Young vs old TARGET median angle {b['median_angle_deg']:.1f}° "
            + ("— not one line across the adult lifespan." if b["median_angle_deg"] > 60
               else "— below the 60° bar, but not assumed identical.")
        )
    if trust:
        lines.append("**How far can it be trusted?**")
        lines.append("")
        for t in trust:
            lines.append(f"- {t}")
        lines.append("")
    lines.append(_assumptions(V1, V2, V3, V4))
    lines.append("")


def _assumptions(V1, V2, V3, V4):
    return (
        "**What would have to be true for a perturbation mapped onto this vector to actually "
        "rejuvenate a cell?** A later step inherits all of the following, none of which this "
        "analysis tests:\n"
        "\n"
        "1. **Correlational, not causal.** The direction is fit as a supervised association of "
        "expression with chronological age. Moving a cell along −TARGET changes the *score*; "
        "it does not by itself mean the cell is younger in any biological sense.\n"
        "2. **Linear.** TARGET is a single unit vector in log-CPM space. Real interventions are "
        "nonlinear, dose-dependent, and constrained to a manifold the linear model never sees.\n"
        "3. **Slow natural aging, not fast intervention.** Weights are estimated from adult "
        "cross-sectional aging (decades), not from a perturbation time course. A two-week OSK "
        "pulse need not travel this line.\n"
        "4. **Identity is the 20-type centroid subspace**, not the rest of what a biologist "
        "means by identity (subtype, state, donor, activity). Orthogonal-to-centroids is not "
        "orthogonal-to-identity in the full sense.\n"
        "5. **The right TARGET is the per-type one** if V3 is (b); a consensus move can age "
        "one type while missing another. A perturbation applied to mixed cortex would have to "
        "be type-resolved, or it is not the vector tested here.\n"
        "6. **Site and age-range transport.** A vector that fails HBCC↔MSSM or young↔old is "
        "a description of this cohort's covariance, not a portable rejuvenation axis. "
        "Any mapped perturbation inherits that non-transport.\n"
        "7. **Held-out prediction is necessary and not sufficient.** V2 only asks whether "
        "projection onto TARGET still tracks chronological age. It does not ask whether "
        "moving along −TARGET reverses function, pathology, or epigenetic age.\n"
        "8. **p ≫ n ridge.** Every gene has a weight; many weights are shrinkage artifacts. "
        "The top-50 list is interpretation, not a CRISPR menu."
    )


def _limitations(lines, V1, V2, V3, V4):
    lines.append("## Limitations")
    lines.append("")
    lines.append("1. Two brain banks only (HBCC, MSSM). LOSO is one df of transfer and is reported, not averaged away.")
    lines.append("2. Unmeasured 6-plex hashing pools, PMI, RIN — as in FINDINGS_BRAIN_PHASE1.md.")
    lines.append(f"3. Seed `{TARGET_SEED}` (`numpy.random.default_rng`). Folds match Phase 1 / geometry / trajectory.")
    lines.append("4. No gene selection. Ridge uses SVD-LOO α on globally z-scored genes (no second per-type z-score, so TARGET and P_id share a metric). PLS-1 is the covariance companion.")
    lines.append("5. Donor bootstrap of a p ≫ n direction is a hard test; pairwise unsigned angles near 90° are the random baseline.")
    lines.append("6. Smooth-muscle and perivascular-macrophage pseudobulks are small-n; their directions are expected to be noisier and are not dropped.")
    lines.append("7. V4c Enrichr is over-representation of the largest-weight tails, not a competitive gene-set test, and is not used for selection.")
    lines.append("8. TARGET is aligned with aging (positive = older). An intervention that 'opposes aging' would move against it. Sign is conventional.")
    lines.append("9. Identity subspace rank is n_types−1. Residualizing it does not residualize donor, site, or unmeasured batch.")
    lines.append("10. The 20–50 vs 50–89 split is confounded with bank (young is HBCC-heavy, old is MSSM-heavy); the 86° young–old angle is not a pure lifespan effect.")
    lines.append("11. Within-site TARGET R²=+0.231 is not G1's +0.317: global-z 1D OLS vs per-type-z ridge predictor. TARGET≈age (6°), so identity projection is not the cause of the gap.")
    lines.append("")


def _files(lines):
    lines.append("## Files")
    lines.append("")
    lines.append("| path | content |")
    lines.append("|---|---|")
    lines.append("| `src/target_common.py`, `target_v1.py` … `target_v4.py`, `target_findings.py` | code |")
    lines.append("| `notebooks/target_v1.ipynb` … `target_v5.ipynb` | runnable from a clean checkout |")
    lines.append("| `results/target/` | tables, json, logs, `figures/`, `target_W.npz` |")
    lines.append("| `FINDINGS_TARGET.md` | this file |")
    lines.append("")


if __name__ == "__main__":
    p = write_findings()
    print("wrote", p)

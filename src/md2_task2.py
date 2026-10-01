"""Task 2 — Figure 3G reproduction and instrument comparison on Louvain clusters."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md2_common import (  # noqa: E402
    MD2_DIR, MD2_PROC, MD2_SEED, MD2_BOOT, MIN_CELLS, AGED_LINE, YOUNG_LINE,
    N_PERM, N_BOOT, N_RANDOM_DIR, PREREG_TASK2, PREREG_TASK2_FLAG,
    DECLARED_BEFORE_SCORES_FLAG, AMS_CTRL,
    StopStep, Logger, dump_json, jsonable, load_json, md2_log_banner,
    load_manifest, save_manifest, record_failure, log_columns,
    progress_snapshot, FROZEN_RULER,
)
from fibro2_common import load_frozen_ruler  # noqa: E402
from fibro_common import bootstrap_rho_ci  # noqa: E402
from gtex_common import spearman_safe  # noqa: E402
from trajectory_common import permutation_p  # noqa: E402
from md2_score import (  # noqa: E402
    lognormalize_csr, add_module_score, score_frozen_pseudobulk,
    map_features, draw_size_matched_indices, ams_from_indices,
)
from md2_task1 import _csr_from_npz  # noqa: E402


STATE_NAMES = ("Fibroblast", "PartialReprog", "EarlyPluripotency", "Pluripotency", "NonReprog")


def _require():
    if not PREREG_TASK2_FLAG.exists():
        raise StopStep("prereg", "PREREG_TASK2.flag missing")
    if not DECLARED_BEFORE_SCORES_FLAG.exists():
        raise StopStep("prereg", "DECLARED_BEFORE_SCORES.flag missing")
    if not (MD2_DIR / "t1_cluster_table.csv").exists():
        raise StopStep("task2", "t1_cluster_table.csv missing — run Task 1")
    if not FROZEN_RULER.exists():
        raise StopStep("frozen_ruler", f"missing {FROZEN_RULER}")


def _score_young_louvain(log, sets, frozen):
    Xl, symbols, _ = _csr_from_npz(MD2_PROC / "louvain_counts_GM23815.npz")
    obs = pd.read_csv(MD2_DIR / "louvain_obs_GM23815.csv")
    if len(obs) != Xl.shape[0]:
        raise StopStep("task2", f"young louvain obs n={len(obs)} X n={Xl.shape[0]}")
    gene_sets = {"md_score": list(sets["MD"]), "tgfb_score": list(sets["TGFB"])}
    for st, genes in (sets.get("reprog") or {}).items():
        gene_sets[f"state_{st}"] = list(genes)
    logX = lognormalize_csr(Xl)
    ams, ams_meta = add_module_score(logX, symbols, gene_sets, log=log, tag="GM23815_louvain")
    yz = np.load(MD2_PROC / "rulerY_GM23815.npz")
    Yall = sparse.csr_matrix(
        (yz["data"], yz["indices"], yz["indptr"]),
        shape=tuple(int(x) for x in yz["shape"]),
    )
    obs_all = pd.read_csv(MD2_DIR / "allcell_obs_GM23815.csv")
    key_all = obs_all["gsm"].astype(str) + "||" + obs_all["barcode"].astype(str)
    pos = {k: i for i, k in enumerate(key_all)}
    key = obs["gsm"].astype(str) + "||" + obs["barcode"].astype(str)
    miss = [k for k in key if k not in pos]
    if miss:
        raise StopStep("task2", f"young louvain {len(miss)} cells not in allcell obs")
    take = np.array([pos[k] for k in key], dtype=int)
    Y = Yall[take]
    rows, mats = [], []
    for (day, cl), g in obs.groupby(["day", "cluster"], sort=True):
        idx = g.index.to_numpy()
        n = int(len(idx))
        n_tp = int((obs.day == day).sum())
        rec = dict(
            resolution="louvain", cell_line=YOUNG_LINE, day=int(day), cluster=int(cl),
            n_cells=n, n_timepoint=n_tp,
            frac_timepoint=float(n / n_tp) if n_tp else np.nan,
            median_UMI=float(np.median(g.umi.to_numpy(float))) if n and "umi" in g.columns else np.nan,
            median_genes=float(np.median(g.n_genes.to_numpy(float))) if n and "n_genes" in g.columns else np.nan,
            mito_frac_median_cell=float(np.nanmedian(g.mito_frac.to_numpy(float))) if n and "mito_frac" in g.columns else np.nan,
            below_min_cells=bool(n < MIN_CELLS),
        )
        for k, arr in ams.items():
            rec[k] = float(np.mean(arr[idx])) if n else np.nan
        rows.append(rec)
        mats.append(np.asarray(Y[idx].sum(axis=0), dtype=np.float64).ravel())
    df = pd.DataFrame(rows)
    C = np.vstack(mats)
    scored = score_frozen_pseudobulk(C, frozen, log, tag="t2_louvain_GM23815")
    df["age_score"] = scored["age"]
    df["pluri_with_OSKM"] = scored["pluri_with"]
    df["pluri_without_OSKM"] = scored["pluri_without"]
    df["pluri_primary"] = scored["pluri_with"]
    df["age_source"] = "TMM_log2CPM_frozen_ruler"
    df.loc[df.below_min_cells, ["age_score", "pluri_with_OSKM", "pluri_without_OSKM", "pluri_primary"]] = np.nan
    for k in list(ams):
        df.loc[df.below_min_cells, k] = np.nan
    np.savez_compressed(
        MD2_PROC / "louvain_ams_GM23815.npz",
        md=ams["md_score"], tgfb=ams["tgfb_score"],
        bins=ams_meta["bins"], gene_mean=ams_meta["gene_mean"],
        cluster=obs["cluster"].to_numpy(int), day=obs["day"].to_numpy(int),
    )
    return df, obs, ams, ams_meta, logX, symbols


def _label_clusters(obs, ams, log, line):
    """Argmax of mean state AddModuleScore per cluster (all timepoints)."""
    recs = []
    for cl, g in obs.groupby("cluster"):
        idx = g.index.to_numpy()
        means = {}
        for st in STATE_NAMES:
            key = f"state_{st}"
            if key not in ams:
                raise StopStep("task2", f"missing {key} in AMS")
            means[st] = float(np.mean(ams[key][idx]))
        label = max(means, key=means.get)
        recs.append(dict(
            cell_line=line, cluster=int(cl), n_cells=int(len(idx)),
            label=label, **{f"mean_{k}": v for k, v in means.items()},
        ))
        log(f"[t2 label] {line} c{int(cl)} n={len(idx)} label={label} means={ {k: round(v,4) for k,v in means.items()} }")
    return pd.DataFrame(recs)


def _reproduction(obs, ams, labels, log):
    """Aged donor. Matches iff PartialReprog and NonReprog exist with ≥20 cells and mean MD Partial < NonReprog."""
    lab_map = {int(r.cluster): str(r.label) for _, r in labels.iterrows()}
    obs = obs.copy()
    obs["state_label"] = obs["cluster"].map(lab_map)
    n_partial = int((obs.state_label == "PartialReprog").sum())
    n_non = int((obs.state_label == "NonReprog").sum())
    md = np.asarray(ams["md_score"], float)
    md_p = float(np.mean(md[obs.state_label.eq("PartialReprog").to_numpy()])) if n_partial else np.nan
    md_n = float(np.mean(md[obs.state_label.eq("NonReprog").to_numpy()])) if n_non else np.nan
    has_both = bool(n_partial >= MIN_CELLS and n_non >= MIN_CELLS)
    md_lower = bool(has_both and np.isfinite(md_p) and np.isfinite(md_n) and md_p < md_n)
    match = bool(has_both and md_lower)
    n_by = obs.groupby("state_label").size().to_dict()
    rec = dict(
        n_PartialReprog=n_partial, n_NonReprog=n_non,
        md_mean_PartialReprog=md_p, md_mean_NonReprog=md_n,
        md_Partial_minus_NonReprog=float(md_p - md_n) if np.isfinite(md_p) and np.isfinite(md_n) else np.nan,
        n_cells_by_label={str(k): int(v) for k, v in n_by.items()},
        min_cells=MIN_CELLS,
        has_both_labels_ge_min=has_both,
        md_partial_lower=md_lower,
        matches_reported_direction=match,
        criterion=(
            "(a) ≥1 PartialReprog-assigned cluster and ≥1 NonReprog-assigned cluster "
            f"with ≥{MIN_CELLS} cells in the aged donor; (b) mean per-cell MD in "
            "PartialReprog-assigned clusters < mean per-cell MD in NonReprog-assigned clusters."
        ),
    )
    log(f"[t2 reproduction] match={match} n_Partial={n_partial} n_Non={n_non} "
        f"MD_P={md_p:+.4f} MD_N={md_n:+.4f}")
    return rec, obs


def _crosscorr(table, line, log):
    g = table[(table.cell_line == line) & (table.below_min_cells == False)].copy()  # noqa: E712
    g = g.dropna(subset=["age_score", "md_score"])
    n = int(len(g))
    rec = dict(cell_line=line, n=n, min_cells=MIN_CELLS, n_perm=N_PERM, n_boot=N_BOOT,
               seed=MD2_SEED, boot_seed=MD2_BOOT)
    if n < 8:
        rec.update(rho=np.nan, rho_null=np.nan, rho_p=np.nan, rho_ci_lo=np.nan, rho_ci_hi=np.nan,
                   note="n<8; spearman_safe not applied")
        log(f"[t2 crosscorr] {line} n={n}<8 skip")
        return rec
    age = g.age_score.to_numpy(float)
    md = g.md_score.to_numpy(float)
    rho = spearman_safe(age, md)
    rng = np.random.default_rng(MD2_SEED)
    nulls = []
    for _ in range(N_PERM):
        md_p = md.copy()
        rng.shuffle(md_p)
        nulls.append(spearman_safe(age, md_p))
    nulls = np.asarray(nulls, float)
    p = permutation_p(rho, nulls, greater=True)
    rho_null = float(np.nanmedian(nulls)) if np.isfinite(nulls).any() else np.nan
    rng_b = np.random.default_rng(MD2_BOOT)
    ci = bootstrap_rho_ci(age, md, rng_b, n_boot=N_BOOT)
    rec.update(
        rho=float(rho) if np.isfinite(rho) else np.nan,
        rho_null=rho_null, rho_p=float(p) if np.isfinite(p) else np.nan,
        rho_ci_lo=ci.get("p025"), rho_ci_hi=ci.get("p975"),
        note="permutation shuffles MD, holds age; unit=cluster×timepoint",
    )
    log(f"[t2 crosscorr] {line} n={n} ρ={rho:+.3f} null={rho_null:+.3f} p={p:+.3f} "
        f"CI=[{ci.get('p025')}, {ci.get('p975')}]")
    return rec


def _md_composition_null(obs, logX, symbols, sets, bins, log, line):
    """Per cluster d0→d7 / d0→d10 MD decline vs 200 size-matched random gene-set AMS."""
    feat_idx, mapped, missing = map_features(symbols, sets["MD"], log, f"{line}:MD_null")
    rng = np.random.default_rng(MD2_SEED)
    draws = draw_size_matched_indices(feat_idx, bins, rng, N_RANDOM_DIR)
    # real MD already in obs via cluster table; recompute per-cell from same AMS object
    # We need per-cell real MD: re-run AMS for MD only using stored bins? Use logX.
    real_ams, _ = add_module_score(
        logX, symbols, {"MD": list(sets["MD"])}, log=log, tag=f"{line}_MDreal",
        bins=bins,
    )
    md_cell = real_ams["MD"]
    rng_ams = np.random.default_rng(MD2_SEED + 17)
    # Precompute null cell scores: n_cells × n_draws would be large (22k × 200).
    # Compute cluster×day means on the fly.
    clusters = sorted(obs.cluster.astype(int).unique().tolist())
    days_u = sorted(obs.day.astype(int).unique().tolist())
    keys = [(int(cl), int(d)) for cl in clusters for d in days_u]
    def mean_at(scores, cl, day, min_n=MIN_CELLS):
        m = (obs.cluster.astype(int) == int(cl)) & (obs.day.astype(int) == int(day))
        n = int(m.sum())
        if n < min_n:
            return np.nan, n
        return float(np.mean(scores[m.to_numpy()])), n

    # null means: draws × keys
    null_mean = {k: np.full(N_RANDOM_DIR, np.nan) for k in keys}
    for d in range(N_RANDOM_DIR):
        sc = ams_from_indices(logX, draws[d], bins, rng_ams, ctrl=AMS_CTRL)
        for cl, day in keys:
            val, n = mean_at(sc, cl, day)
            null_mean[(cl, day)][d] = val
        if d % 20 == 0:
            log(f"[t2 MD-null] {line} draw {d}/{N_RANDOM_DIR}")

    rows = []
    n_move = 0
    for cl in clusters:
        for dend in (7, 10):
            m0, n0 = mean_at(md_cell, cl, 0)
            m1, n1 = mean_at(md_cell, cl, dend)
            rec = dict(
                cell_line=line, cluster=int(cl), endpoint=f"d0→d{dend}", d_end=int(dend),
                n_cells_d0=n0, n_cells_end=n1, n_random=N_RANDOM_DIR, seed=MD2_SEED,
                instrument="MD", null_kind="size_matched_random_gene_sets_mean_expression_bin",
                n_MD_genes=int(feat_idx.size), min_cells=MIN_CELLS,
            )
            if not np.isfinite(m0) or not np.isfinite(m1):
                rec.update(ok=False, reason=f"n_cells<{MIN_CELLS} at an endpoint")
                rows.append(rec)
                continue
            real = float(m0 - m1)  # positive = MD fell
            null_d = null_mean[(cl, 0)] - null_mean[(cl, dend)]
            if not np.isfinite(null_d).all():
                rec.update(ok=False, reason="null not finite at an endpoint")
                rows.append(rec)
                continue
            p = permutation_p(real, null_d, greater=True)
            n_ge = int(np.sum(null_d >= real))
            pass_p = bool(np.isfinite(p) and p <= 0.05)
            rec.update(
                ok=True, md_d0=m0, md_end=m1, decline=real, p=p,
                n_random_ge_real=n_ge, pass_p=pass_p,
            )
            rows.append(rec)
            if pass_p:
                n_move += 1
            log(f"[t2 MD-null] {line} c{cl} d0→d{dend} decline={real:+.4f} p={p:+.3f} "
                f"n_ge={n_ge}/{N_RANDOM_DIR} n0={n0} n1={n1}")
    return pd.DataFrame(rows), int(n_move)


def _t2_reading(repro, t1_key, md_null_aged, ruler_moved):
    if md_null_aged is None or not len(md_null_aged):
        md_any = False
        md_n_pass = 0
    else:
        ok = md_null_aged[md_null_aged.ok == True]  # noqa: E712
        md_n_pass = int((ok.pass_p == True).sum()) if len(ok) else 0  # noqa: E712
        md_any = md_n_pass > 0
    if not repro.get("matches_reported_direction"):
        return dict(
            key="reproduction_failure",
            text=(
                "Our reproduction does not recover their published Figure 3G direction "
                "(PartialReprog mean MD < NonReprog mean MD on the aged donor, both labels ≥20 cells). "
                "The instrument comparison is not valid and was not run as a claim. "
                "This is a reproduction failure, not a refutation."
            ),
            comparison_valid=False,
            md_moves_per_cluster=md_any, ruler_moved_task1=ruler_moved,
            md_n_pass=md_n_pass, t1_key=t1_key,
        )
    # extra reading: MD moves but fails null — that's md_any False while raw MD decline exists.
    # "fails its size-matched random-gene-set null" = movement in the point estimate that is not
    # significant vs the gene-set null. The pre-reg defines "MD moves per cluster" as p≤0.05 on that null.
    if md_any and ruler_moved:
        key = "both_agree"
        text = (
            "MD moves per cluster and the frozen ruler does too (Task 1 outcome 1). "
            "Both instruments agree; their bend claim is supported and our earlier negative "
            "was our clustering."
        )
    elif md_any and not ruler_moved:
        key = "instruments_disagree"
        text = (
            "MD moves per cluster, the frozen ruler does not. The instruments disagree on the "
            "same cells. What would settle it: an independent fibroblast age instrument scored "
            "on these same Louvain clusters, or the authors' Seurat/sctransform object with the "
            "frozen ruler projected onto it. No winner is declared."
        )
    elif (not md_any) and (not ruler_moved):
        key = "reproduction_neither_moves"
        text = (
            "Neither moves per cluster. Our reproduction does not recover their published result "
            "on the MD instrument at the per-cluster null; report the reproduction failure as the "
            "finding and do not present it as a refutation."
        )
    else:
        # ruler moved, MD did not — not a listed bullet; report both, no winner.
        key = "instruments_disagree"
        text = (
            "The frozen ruler moved per cluster (Task 1 outcome 1) and MD did not at p≤0.05 on "
            "the size-matched gene-set null. The instruments disagree on the same cells. No winner "
            "is declared."
        )
    # If MD point estimates decline but null fails: additional strong claim.
    extra = None
    if md_null_aged is not None and len(md_null_aged):
        ok = md_null_aged[md_null_aged.ok == True]  # noqa: E712
        declined = ok[ok.decline > 0] if len(ok) and "decline" in ok.columns else ok.iloc[0:0]
        if len(declined) and md_n_pass == 0:
            extra = (
                "MD point-estimate decline exists in some cluster×endpoint but fails its "
                "size-matched random-gene-set null (no p≤0.05). Their metric's movement is not "
                "specific to mesenchymal genes in this dataset. See t2_md_null.csv for the numbers."
            )
            key = "md_fails_geneset_null"
            text = extra
    return dict(
        key=key, text=text, comparison_valid=True,
        md_moves_per_cluster=md_any, ruler_moved_task1=ruler_moved,
        md_n_pass=md_n_pass, t1_key=t1_key, extra=extra,
    )


def run_task2(log=None):
    _require()
    close_log = False
    if log is None:
        log = Logger(MD2_DIR / "t2_report.txt")
        close_log = True
    md2_log_banner(log, "TASK2")
    log(PREREG_TASK2)
    sets = load_json(MD2_DIR / "genesets.json")
    frozen = load_frozen_ruler()
    t1 = pd.read_csv(MD2_DIR / "t1_cluster_table.csv")
    log_columns("t1_cluster_table", list(t1.columns), str(MD2_DIR / "t1_cluster_table.csv"))
    aged_lv = t1[(t1.resolution.astype(str) == "louvain") & (t1.cell_line.astype(str) == AGED_LINE)].copy()
    t1_read = load_json(MD2_DIR / "t1_reading.json") if (MD2_DIR / "t1_reading.json").exists() else {}
    ruler_moved = t1_read.get("key") == "resolution_artifact"

    obs_a = pd.read_csv(MD2_DIR / "louvain_obs_GM00731.csv")
    ams_a = np.load(MD2_PROC / "louvain_ams_GM00731.npz", allow_pickle=True)
    st = np.load(MD2_PROC / "louvain_state_GM00731.npz", allow_pickle=True)
    ams_dict = dict(md_score=np.asarray(ams_a["md"]), tgfb_score=np.asarray(ams_a["tgfb"]))
    for k in st.files:
        ams_dict[k] = np.asarray(st[k])
    labels_a = _label_clusters(obs_a, ams_dict, log, AGED_LINE)
    labels_a.to_csv(MD2_DIR / "t2_cluster_labels_GM00731.csv", index=False)
    repro, obs_lab = _reproduction(obs_a, ams_dict, labels_a, log)
    dump_json(MD2_DIR / "t2_reproduction.json", jsonable(repro))

    young_df, obs_y, ams_y, ams_y_meta, logX_y, symbols_y = _score_young_louvain(log, sets, frozen)
    labels_y = _label_clusters(obs_y, ams_y, log, YOUNG_LINE)
    labels_y.to_csv(MD2_DIR / "t2_cluster_labels_GM23815.csv", index=False)

    # attach labels to cluster tables
    aged_lv = aged_lv.merge(labels_a[["cluster", "label"]], on="cluster", how="left")
    young_df = young_df.merge(labels_y[["cluster", "label"]], on="cluster", how="left")
    t2tab = pd.concat([aged_lv, young_df], axis=0, ignore_index=True)
    t2tab.to_csv(MD2_DIR / "t2_cluster_table.csv", index=False)

    if not repro.get("matches_reported_direction"):
        reading = _t2_reading(repro, t1_read.get("key"), None, ruler_moved)
        dump_json(MD2_DIR / "t2_reading.json", jsonable(reading))
        dump_json(MD2_DIR / "t2_summary.json", jsonable(dict(
            reproduction=repro, reading=reading, comparison_valid=False,
            note="Task 2 stopped after reproduction check.",
        )))
        log("[t2] reproduction does not match. Stopping Task 2 comparison.")
        man = load_manifest()
        man["status"] = "TASK2_STOPPED_REPRODUCTION"
        save_manifest(man)
        progress_snapshot("Task 3 instrument validity (report-only)", stop="TASK2_STOPPED_REPRODUCTION")
        if close_log:
            log.close()
        return reading

    # Cross-correlate per donor
    cc_rows = [_crosscorr(t2tab, AGED_LINE, log), _crosscorr(t2tab, YOUNG_LINE, log)]
    pd.DataFrame(cc_rows).to_csv(MD2_DIR / "t2_crosscorr.csv", index=False)

    # MD composition null on aged Louvain object
    Xl, symbols_a, _ = _csr_from_npz(MD2_PROC / "louvain_counts_GM00731.npz")
    logX_a = lognormalize_csr(Xl)
    bins_a = np.asarray(ams_a["bins"])
    md_null_a, n_move_a = _md_composition_null(obs_a, logX_a, symbols_a, sets, bins_a, log, AGED_LINE)
    md_null_y, n_move_y = _md_composition_null(obs_y, logX_y, symbols_y, sets, ams_y_meta["bins"], log, YOUNG_LINE)
    md_null = pd.concat([md_null_a, md_null_y], axis=0, ignore_index=True)
    md_null.to_csv(MD2_DIR / "t2_md_null.csv", index=False)

    reading = _t2_reading(repro, t1_read.get("key"), md_null_a, ruler_moved)
    dump_json(MD2_DIR / "t2_reading.json", jsonable(reading))
    dump_json(MD2_DIR / "t2_summary.json", jsonable(dict(
        reproduction=repro, reading=reading, comparison_valid=True,
        md_n_pass_aged=int(n_move_a), md_n_pass_young=int(n_move_y),
        ruler_moved_task1=ruler_moved, t1_key=t1_read.get("key"),
        null_difference=(
            "MD null = 200 size-matched random gene sets (same n, same mean-expression bin, "
            "AddModuleScore). Frozen-ruler null = 200 permuted ruler weights. "
            "MD is a fixed gene list, so permuting weights is not the matching null."
        ),
    )))
    man = load_manifest()
    man["status"] = "TASK2_DONE"
    man["t2_reading"] = reading.get("key")
    save_manifest(man)
    progress_snapshot("Task 3 instrument validity (report-only)", stop="TASK2_DONE")
    if close_log:
        log.close()
    return reading


if __name__ == "__main__":
    try:
        run_task2()
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise

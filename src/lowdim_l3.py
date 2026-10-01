"""STAGE L3 — TARGET in the winning low-dimensional space.

Using the L2 winner:
  - Remove the identity (type-centroid) subspace from the gene-mapped age
    direction, renormalize → TARGET.
  - Re-run the V3 test: pairwise angles among per-type TARGETs vs a
    permuted-age null (gene space was 86.0° vs null 85.8°, p=0.73).
  - If a consensus exists, report how much per-type signal it retains (OOF).
  - Map TARGET back to gene weights; top 50 ±; pathway enrichment; junk check.

Skipped if STOP L2.

Usage: python src/lowdim_l3.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lowdim_common import (  # noqa: E402
    LD_DIR, LD_FIG, LD_SEED, DATASET_ID, KMAX, N_PERM, N_PERM_PAIR, MIN_N_TYPE,
    STOP_NULL, SITE_NAME,
    Logger, dump_json, load_json, ld_log_banner, load_bundle, load_gene_sets,
    zscore_train, make_folds, fit_spaces_by_type, apply_slice, age_direction,
    gene_direction_from_pack, direction_and_pred, target_from_age, identity_basis,
    permute_age_within_site_obs, score_age_predictions, permutation_p,
    summarize_null_col, consensus_direction, pairwise_angles_deg, pairwise_median_deg,
    unsigned_angle_deg, unit, junk_symbol, calibrate_1d, config_id, fit_pack,
)
from gene_lists import is_ribo, is_mito, CELL_CYCLE_ALL  # noqa: E402

ENRICHR_ADD = "https://maayanlab.cloud/Enrichr/addList"
ENRICHR_ENRICH = "https://maayanlab.cloud/Enrichr/enrich"
ENRICHR_LIBS = (
    "MSigDB_Hallmark_2020",
    "GO_Biological_Process_2023",
    "KEGG_2021_Human",
    "Reactome_2022",
)
AGING_TOKENS = (
    "inflamm", "interferon", "complement", "cytokine", "nf-kappa", "tnf",
    "senescen", "p53", "apoptos", "proteostasis", "unfolded protein", "ubiquitin",
    "proteasome", "autophagy", "lysosom", "mitochond", "oxidative phosphoryl",
    "electron transport", "reactive oxygen", "glutathione", "glycolysis",
    "mtorc", "insulin", "igf", "collagen", "extracellular matrix", "synapse",
    "myelination", "oligodendro",
)


def _winner_or_none(log):
    path = LD_DIR / "l2_winner.json"
    stop_path = LD_DIR / "l2_STOP.txt"
    if stop_path.exists() and "stop_L2=True" in stop_path.read_text(encoding="utf-8"):
        log("[L3] STOP L2 is set — not constructing a TARGET from a non-identifiable direction.")
        return None
    if not path.exists():
        log("[L3] l2_winner.json missing")
        return None
    rec = load_json(path)
    if rec.get("stop_L2") or rec.get("winner") is None:
        return None
    return rec["winner"]


def _age_gene_dirs(pack, types, k, clean, y_override=None):
    """Gene-mapped unit age directions, one per type, from a fitted pack."""
    y_orig = pack["y"]
    if y_override is not None:
        pack = dict(pack)
        pack["y"] = np.asarray(y_override, float)
    out = {}
    meta = []
    p = pack["Xz"].shape[1]
    for t in types:
        sp = pack["spaces"].get(t)
        tr_m = pack["ct"] == t
        n = int(tr_m.sum())
        if sp is None or n < MIN_N_TYPE:
            meta.append(dict(celltype=t, n=n, ok=False))
            continue
        sl = apply_slice(sp, pack["Xz"][tr_m], pack["Xz"][:0], k, pack["tech"][tr_m], clean)
        if sl["k_keep"] < 1:
            meta.append(dict(celltype=t, n=n, ok=False, k_keep=0))
            continue
        w, rec = age_direction(sl["Ztr"], pack["y"][tr_m], method="ridge")
        v = gene_direction_from_pack(sl, w, p)
        ok = rec.get("ok") and np.isfinite(v).all() and float(np.linalg.norm(v)) > 1e-12
        if ok:
            v = v / np.linalg.norm(v)
            out[t] = v
        meta.append(dict(celltype=t, n=n, ok=bool(ok), k_keep=int(sl["k_keep"]),
                         alpha=rec.get("alpha", np.nan), reverted=bool(sl["reverted"])))
    pack["y"] = y_orig
    return out, pd.DataFrame(meta)


def _targets_from_dirs(dirs, types, basis):
    p = basis.shape[0] if basis.size else 0
    W = np.full((p, len(types)), np.nan)
    rows = []
    for j, t in enumerate(types):
        v = dirs.get(t)
        if v is None or p == 0:
            rows.append(dict(celltype=t, ok=False, norm_retained=np.nan,
                             angle_to_age_deg=np.nan, proj_r2_id=np.nan))
            continue
        tgt, retained, ang, pr = target_from_age(v, basis)
        W[:, j] = tgt
        rows.append(dict(celltype=t, ok=True, norm_retained=retained,
                         angle_to_age_deg=ang, proj_r2_id=pr))
    return W, pd.DataFrame(rows)


def _pairwise_stats(W):
    med, ang = pairwise_median_deg(W)
    ok = [j for j in range(W.shape[1])
          if np.isfinite(W[:, j]).all() and float(np.linalg.norm(W[:, j])) > 1e-12]
    pc1 = np.nan
    w_cons = None
    if len(ok) >= 2:
        Wc = np.column_stack([W[:, j] / np.linalg.norm(W[:, j]) for j in ok])
        w_cons, pc1 = consensus_direction(Wc)
        dots = Wc.T @ w_cons
        if float(np.nanmean(dots)) < 0:
            w_cons = -w_cons
    return dict(median_pairwise_deg=med, shared_pc1_frac=pc1, n_ok=len(ok),
                w_cons=w_cons, ang=ang)


def _top_genes(w, symbols, gene_ids, k=50):
    w = np.asarray(w, float)
    m = np.isfinite(w)
    order = np.argsort(w)
    neg_i = order[m[order]][:k]
    pos_i = order[m[order]][-k:][::-1]

    def _rows(idx, sign):
        return pd.DataFrame(dict(
            rank=np.arange(1, len(idx) + 1), sign=sign,
            symbol=symbols[idx], gene_id=gene_ids[idx], weight=w[idx],
        ))
    return pd.concat([_rows(pos_i, "+"), _rows(neg_i, "-")], ignore_index=True)


def _enrichr(genes, desc, log):
    import requests
    genes = [g for g in genes if isinstance(g, str) and g and g != "nan"]
    if len(genes) < 5:
        return []
    try:
        r = requests.post(ENRICHR_ADD, files={
            "list": (None, "\n".join(genes)),
            "description": (None, desc),
        }, timeout=60)
        r.raise_for_status()
        uid = r.json()["userListId"]
    except Exception as e:
        log(f"   Enrichr addList failed ({desc}): {e}")
        return []
    hits = []
    for lib in ENRICHR_LIBS:
        try:
            er = requests.get(ENRICHR_ENRICH, params={"userListId": uid, "backgroundType": lib},
                              timeout=60)
            er.raise_for_status()
            data = er.json()
            rows = data.get(lib, [])
            for row in rows[:25]:
                term = row[1]
                p = float(row[2])
                adj = float(row[6]) if len(row) > 6 and row[6] is not None else np.nan
                glist = row[5] if len(row) > 5 else ""
                hits.append(dict(
                    list=desc, library=lib, term=term, p=p, adj_p=adj,
                    genes=glist if isinstance(glist, str) else ",".join(glist),
                    aging_token=any(tok in str(term).lower() for tok in AGING_TOKENS),
                ))
        except Exception as e:
            log(f"   Enrichr enrich failed ({desc} / {lib}): {e}")
    return hits


def _oof_consensus(Y, obs, tech, types, sets, variant, k, clean, rng, log):
    log("   L3 OOF consensus vs own TARGET (nested train-fold spaces)")
    folds_all = make_folds(obs, np.random.default_rng(LD_SEED))
    y = obs.age.to_numpy(float)
    out_rows = []
    for scheme, folds in (("within_site", folds_all["within_site"]), ("loso", folds_all["loso"])):
        pred_own = np.full(len(obs), np.nan)
        pred_cons = np.full(len(obs), np.nan)
        for fi, fold in enumerate(folds):
            tr, te = fold[0], fold[1]
            pack = fit_pack(Y, obs, tech, tr, types, variant, sets)
            Xte = (Y[te] - pack["mu"]) / pack["sd"]
            ct_te = obs.celltype.astype(str).to_numpy()[te]
            basis, _, _ = identity_basis(pack["Xz"], pack["ct"], types)
            dirs, _ = _age_gene_dirs(pack, types, k, clean)
            W, _geom = _targets_from_dirs(dirs, types, basis)
            st = _pairwise_stats(W)
            w_cons = st["w_cons"]
            for j, t in enumerate(types):
                tr_m = pack["ct"] == t
                te_m = ct_te == t
                if int(tr_m.sum()) < MIN_N_TYPE or int(te_m.sum()) < 1:
                    continue
                sp = pack["spaces"].get(t)
                if sp is None:
                    continue
                sl = apply_slice(sp, pack["Xz"][tr_m], Xte[te_m], k, pack["tech"][tr_m], clean)
                if sl["k_keep"] < 1:
                    continue
                ytr = pack["y"][tr_m]
                w_t = W[:, j]
                # project test gene-z onto TARGET (gene weights) and onto consensus
                if np.isfinite(w_t).all() and float(np.linalg.norm(w_t)) > 1e-12:
                    s_tr = pack["Xz"][tr_m] @ w_t
                    s_te = Xte[te_m] @ w_t
                    pte, _, _ = calibrate_1d(s_tr, ytr, s_te)
                    pred_own[te[te_m]] = pte
                if w_cons is not None:
                    s_tr = pack["Xz"][tr_m] @ w_cons
                    s_te = Xte[te_m] @ w_cons
                    pte, _, _ = calibrate_1d(s_tr, ytr, s_te)
                    pred_cons[te[te_m]] = pte
            log(f"      {scheme} fold {fi+1}/{len(folds)}")
        sc_own = score_age_predictions(obs, pred_own)
        sc_cons = score_age_predictions(obs, pred_cons)
        prim_own = sc_own["site_strat_median_type_r2"]
        prim_cons = sc_cons["site_strat_median_type_r2"]
        pt_own = sc_own["per_type"]
        pt_cons = sc_cons["per_type"]
        own_m = pt_own[pt_own.site == "mean_of_sites"].set_index("celltype")["r2"]
        cons_m = pt_cons[pt_cons.site == "mean_of_sites"].set_index("celltype")["r2"]
        per = []
        for t in types:
            a = float(own_m.get(t, np.nan)) if t in own_m.index else np.nan
            b = float(cons_m.get(t, np.nan)) if t in cons_m.index else np.nan
            per.append(dict(scheme=scheme, celltype=t, r2_own=a, r2_consensus=b,
                            loss=(a - b) if np.isfinite(a) and np.isfinite(b) else np.nan,
                            ratio=(b / a) if np.isfinite(a) and abs(a) > 1e-12 else np.nan))
        per_df = pd.DataFrame(per)
        per_df.to_csv(LD_DIR / f"l3_v3c_{scheme}_per_type.csv", index=False)
        med_loss = float(per_df.loss.median()) if per_df.loss.notna().any() else np.nan
        med_ratio = float(per_df.ratio.median()) if per_df.ratio.notna().any() else np.nan
        log(f"      {scheme}: own R²={prim_own:+.3f}  consensus R²={prim_cons:+.3f}  "
            f"median loss={med_loss:+.3f}  median ratio={med_ratio:.3f}")
        out_rows.append(dict(
            scheme=scheme, r2_own=prim_own, r2_consensus=prim_cons,
            median_loss=med_loss, median_ratio=med_ratio,
        ))
    return pd.DataFrame(out_rows)


def _run(log, rng):
    ld_log_banner(log, "L3 — TARGET in the winning space")
    winner = _winner_or_none(log)
    if winner is None:
        dump_json(LD_DIR / "l3_summary.json", dict(skipped=True, reason="STOP L2"))
        (LD_DIR / "l3_STOP.txt").write_text("skipped=True\nreason=STOP L2\n", encoding="utf-8")
        return dict(skipped=True)

    variant = winner["variant"]
    k = int(winner["k"])
    clean = bool(winner["clean"])
    cfg = winner["config"]
    log(f"[L3] winner {cfg}  variant={variant} k={k} clean={clean}")

    data = load_bundle(log)
    Y, genes, obs, tech, types = data["Y"], data["genes"], data["obs"], data["tech"], data["types"]
    sets, _ = load_gene_sets(genes, log)
    pack = fit_pack(Y, obs, tech, np.arange(len(obs)), types, variant, sets)
    basis, S_id, used = identity_basis(pack["Xz"], pack["ct"], types)
    log(f"[L3] identity subspace: {basis.shape[1]} axes from {len(used)} type centroids")
    dirs, meta = _age_gene_dirs(pack, types, k, clean)
    W, geom = _targets_from_dirs(dirs, types, basis)
    geom = pd.concat([meta.reset_index(drop=True), geom.drop(columns=["celltype"], errors="ignore")], axis=1)
    geom.to_csv(LD_DIR / "l3_geom_full.csv", index=False)
    log(f"[L3] median norm retained after identity residualization: {geom.norm_retained.median():.3f}  "
        f"median angle to age dir: {geom.angle_to_age_deg.median():.2f}°")

    st = _pairwise_stats(W)
    log(f"[L3] pairwise TARGET median angle={st['median_pairwise_deg']:.1f}°  "
        f"PC1={st['shared_pc1_frac']:.3f}  n_ok={st['n_ok']}")

    log(f"[L3] pairwise permutation null ({N_PERM_PAIR}) — space frozen, ages shuffled")
    null_rows = []
    for i in range(N_PERM_PAIR):
        obs_p = permute_age_within_site_obs(obs, rng)
        yp = obs_p.age.to_numpy(float)
        dirs_p, _ = _age_gene_dirs(pack, types, k, clean, y_override=yp)
        Wp, _ = _targets_from_dirs(dirs_p, types, basis)
        stp = _pairwise_stats(Wp)
        null_rows.append(dict(perm=i, pairwise_median_deg=stp["median_pairwise_deg"],
                              shared_pc1_frac=stp["shared_pc1_frac"]))
        if (i + 1) % 10 == 0:
            log(f"      perm {i+1}/{N_PERM_PAIR}")
    null_df = pd.DataFrame(null_rows)
    null_df.to_csv(LD_DIR / "l3_pairwise_null.csv", index=False)
    p_ang = permutation_p(st["median_pairwise_deg"], null_df.pairwise_median_deg, greater=False)
    p_pc1 = permutation_p(st["shared_pc1_frac"], null_df.shared_pc1_frac, greater=True)
    ns_ang = summarize_null_col(null_df.pairwise_median_deg)
    ns_pc1 = summarize_null_col(null_df.shared_pc1_frac)
    log(f"      observed {st['median_pairwise_deg']:.1f}° vs null {ns_ang['mean']:.1f}°  p={p_ang:.3f}")
    log(f"      PC1 {st['shared_pc1_frac']:.3f} vs null {ns_pc1['mean']:.3f}  p={p_pc1:.3f}")

    type_specific = bool(np.isfinite(p_ang) and p_ang > 0.05)
    # below null (smaller angle) would mean more shared than chance
    shared_excess = bool(np.isfinite(p_ang) and p_ang <= 0.05)
    if shared_excess:
        verdict = "shared excess vs null — type-specificity of gene-space V3 was inflated by non-identifiability"
    else:
        verdict = "pairwise still at null — type-specificity survives the low-dim fix (or both are noise)"
    log(f"[L3] V3 retest: {verdict}")

    symbols = genes.symbol.astype(str).to_numpy()
    gids = genes.gene_id.astype(str).to_numpy()
    w_cons = st["w_cons"]
    frames = []
    if w_cons is not None:
        frames.append(_top_genes(w_cons, symbols, gids).assign(source="consensus"))
    acc = np.zeros(W.shape[0])
    n = 0
    for j in range(W.shape[1]):
        w = W[:, j]
        if np.isfinite(w).all() and float(np.linalg.norm(w)) > 1e-12:
            acc = acc + w / np.linalg.norm(w)
            n += 1
    w_mean = acc / max(n, 1)
    frames.append(_top_genes(w_mean, symbols, gids).assign(source="mean_per_type"))
    for want in ("oligodendrocyte", "astrocyte", "microglial cell",
                 "L2/3 intratelencephalic projecting glutamatergic neuron"):
        if want in types:
            j = types.index(want)
            frames.append(_top_genes(W[:, j], symbols, gids).assign(source=want))
    top = pd.concat(frames, ignore_index=True)
    top.to_csv(LD_DIR / "l3_top_genes.csv", index=False)
    top["junk"] = top.symbol.map(junk_symbol)
    cons = top[top.source == "consensus"]
    n_junk = int(cons.junk.sum()) if len(cons) else 0
    log(f"[L3] consensus top-50 ± junk (LINC/pseudogene/OR/ENSG): {n_junk}/{len(cons)}")

    def _tech(sub):
        sy = sub.symbol.astype(str)
        return dict(
            n=int(len(sy)),
            n_mito=int(sy.map(is_mito).sum()),
            n_ribo=int(sy.map(is_ribo).sum()),
            n_cell_cycle=int(sy.isin(CELL_CYCLE_ALL).sum()),
            n_junk=int(sy.map(junk_symbol).sum()),
        )
    tech_flags = dict(
        consensus_pos=_tech(cons[cons.sign == "+"]) if len(cons) else {},
        consensus_neg=_tech(cons[cons.sign == "-"]) if len(cons) else {},
    )
    dump_json(LD_DIR / "l3_tech_flags.json", tech_flags)

    hits = []
    for src in top.source.unique():
        for sgn, lab in (("+", "pos"), ("-", "neg")):
            sy = top[(top.source == src) & (top.sign == sgn)].symbol.astype(str).tolist()
            desc = f"{src}_{lab}"
            log(f"   Enrichr {desc} (n={len(sy)})")
            hits.extend(_enrichr(sy, desc, log))
    hit_df = pd.DataFrame(hits)
    n_aging = 0
    if len(hit_df):
        hit_df.to_csv(LD_DIR / "l3_enrichr.csv", index=False)
        aging = hit_df[hit_df.aging_token & (hit_df.p < 0.05)]
        n_aging = int(len(aging))
        log(f"      Enrichr terms with aging tokens & p<0.05: {n_aging}")
    else:
        log("      Enrichr returned no terms.")

    v3c = _oof_consensus(Y, obs, tech, types, sets, variant, k, clean, rng, log)

    np.savez_compressed(
        LD_DIR / "target_W.npz",
        W_tgt=W, consensus=w_cons if w_cons is not None else np.zeros(W.shape[0]),
        gene_id=gids, symbol=symbols, types=np.array(types, dtype=object),
        variant=np.array(variant), k=np.array(k), clean=np.array(clean),
        config=np.array(cfg),
        mu=pack["mu"], sd=pack["sd"], basis=basis, S_id=S_id,
    )
    log("[L3] wrote target_W.npz")

    # figures
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    if len(null_df):
        ax.hist(null_df.pairwise_median_deg, bins=15, color="0.75", label="null")
    ax.axvline(st["median_pairwise_deg"], c="tab:red", lw=2,
               label=f"observed {st['median_pairwise_deg']:.1f}°")
    ax.axvline(86.0, c="k", ls="--", lw=1, label="gene-space V3 86.0°")
    ax.set_xlabel("median pairwise TARGET angle (deg)")
    ax.set_title("L3  V3 retest in the winning space")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(LD_FIG / "l3_pairwise_vs_null.png", dpi=140)
    plt.close(fig)

    if len(cons):
        fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.8))
        for ax, sgn, title in zip(axes, ("+", "-"), ("positive (older)", "negative")):
            sub = cons[cons.sign == sgn].head(25)
            y_pos = np.arange(len(sub))
            cols = ["tab:red" if j else "tab:blue" for j in sub.junk]
            ax.barh(y_pos, sub.weight.abs(), color=cols)
            ax.set_yticks(y_pos)
            ax.set_yticklabels(sub.symbol, fontsize=6)
            ax.invert_yaxis()
            ax.set_title(title)
            ax.set_xlabel("|weight|")
        fig.suptitle("L3  consensus top weights (red = LINC/pseudogene/OR/ENSG)", fontsize=10)
        fig.tight_layout()
        fig.savefig(LD_FIG / "l3_top_weights.png", dpi=140)
        plt.close(fig)

    summary = dict(
        seed=LD_SEED, dataset_id=DATASET_ID, skipped=False,
        winner=winner, n_identity_axes=int(basis.shape[1]),
        median_norm_retained=float(geom.norm_retained.median()) if geom.norm_retained.notna().any() else None,
        median_angle_to_age=float(geom.angle_to_age_deg.median()) if geom.angle_to_age_deg.notna().any() else None,
        pairwise_median_deg=st["median_pairwise_deg"],
        pairwise_null_mean=ns_ang["mean"], pairwise_p=p_ang,
        shared_pc1_frac=st["shared_pc1_frac"], pc1_null_mean=ns_pc1["mean"], pc1_p=p_pc1,
        n_ok=st["n_ok"], n_perm=N_PERM_PAIR,
        type_specific_at_null=type_specific, shared_excess=shared_excess, verdict=verdict,
        n_junk_consensus_top100=n_junk, n_consensus_listed=int(len(cons)),
        n_enrichr=int(len(hit_df)), n_aging_token_p05=n_aging,
        tech_flags=tech_flags, v3c=v3c.to_dict(orient="records"),
        gene_space_v3=dict(pairwise=86.0, null=85.8, p=0.73),
    )
    dump_json(LD_DIR / "l3_summary.json", summary)
    log("[L3] wrote results/lowdim/l3_*")
    log("[L3] done.")
    return summary


def run():
    rng = np.random.default_rng(LD_SEED)
    log = Logger(LD_DIR / "l3_report.txt")
    try:
        return _run(log, rng)
    finally:
        log.close()


if __name__ == "__main__":
    run()

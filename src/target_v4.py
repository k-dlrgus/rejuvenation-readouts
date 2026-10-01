"""STAGE V4 — how trustworthy is the target?

V4a. Site transfer: TARGET fit on HBCC → MSSM and vice versa.
V4b. Age-range: fit in 20–50 vs 50–89; angle between the two directions.
V4c. Biological sanity: top 50 ± weights (interpretation only) + pathway enrichment.

Usage: python src/target_v4.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from target_common import (  # noqa: E402
    TARGET_DIR, TARGET_FIG, TARGET_SEED, DATASET_ID, N_PERM, STOP_NULL, MIN_N_TYPE,
    Logger, dump_json, tgt_log_banner, load_phase1_matrix,
    zscore_train, identity_basis, per_type_directions, per_type_targets,
    consensus_direction, calibrate_1d, permute_age_within_site_obs,
    score_age_predictions, primary_r2, permutation_p, summarize_null_col,
    unsigned_angle_deg, unit, r2_mae, fit_full,
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


def _fit_on_mask(Y, obs, mask, method="ridge"):
    idx = np.flatnonzero(mask)
    sub = obs.iloc[idx].reset_index(drop=True)
    Ys = Y[idx]
    ct = sub.celltype.astype(str).to_numpy()
    types_all = sorted(pd.unique(obs.celltype.astype(str)))
    y = sub.age.to_numpy(float)
    Xs, mu, sd = zscore_train(Ys)
    basis, S_id, used = identity_basis(Xs, ct, types_all)
    W_age, meta, built = per_type_directions(Xs, y, ct, types_all, method=method)
    W_tgt, geom = per_type_targets(W_age, basis)
    return dict(
        idx=idx, obs=sub, X=Xs, mu=mu, sd=sd, y=y, ct=ct, types=types_all,
        basis=basis, W_age=W_age, W_tgt=W_tgt, meta=meta, geom=geom, svd=built,
        used=used, S_id=S_id,
    )


def _apply_and_score(pack_tr, Y, obs, test_mask, types):
    """Z-score test with train mu/sd; 1D OLS calibrated on train projections."""
    te = np.flatnonzero(test_mask)
    Xte = (Y[te] - pack_tr["mu"]) / pack_tr["sd"]
    yte = obs.age.to_numpy(float)[te]
    ct_te = obs.celltype.astype(str).to_numpy()[te]
    ct_tr = pack_tr["ct"]
    ytr = pack_tr["y"]
    Xtr = pack_tr["X"]
    pred = np.full(len(te), np.nan)
    rows = []
    for j, t in enumerate(types):
        w = pack_tr["W_tgt"][:, j]
        tr_m = ct_tr == t
        te_m = ct_te == t
        if int(tr_m.sum()) < MIN_N_TYPE or int(te_m.sum()) < 3:
            rows.append(dict(celltype=t, r2=np.nan, mae=np.nan, n_te=int(te_m.sum()), n_tr=int(tr_m.sum())))
            continue
        if not np.isfinite(w).all() or float(np.linalg.norm(w)) < 1e-12:
            rows.append(dict(celltype=t, r2=np.nan, mae=np.nan, n_te=int(te_m.sum()), n_tr=int(tr_m.sum())))
            continue
        pte, _, _ = calibrate_1d(Xtr[tr_m] @ w, ytr[tr_m], Xte[te_m] @ w)
        pred[te_m] = pte
        met = r2_mae(yte[te_m], pte)
        rows.append(dict(celltype=t, r2=met["r2"], mae=met["mae"], n_te=met["n"], n_tr=int(tr_m.sum())))
    tab = pd.DataFrame(rows)
    med = float(tab.r2.median()) if tab.r2.notna().any() else np.nan
    return pred, yte, tab, med


def _site_transfer(Y, obs, rng, log):
    site = obs.Source.astype(str).to_numpy()
    sites = sorted(pd.unique(site))
    log(f"[V4a] sites={ {s: int((obs.loc[site==s, 'donor'].nunique())) for s in sites} }")
    types = sorted(pd.unique(obs.celltype.astype(str)))
    recs = []
    per_tabs = []
    for tr_s in sites:
        te_s = [s for s in sites if s != tr_s]
        if not te_s:
            continue
        te_s = te_s[0]
        log(f"   fit on {tr_s} → test {te_s}")
        pack = _fit_on_mask(Y, obs, site == tr_s)
        pred, yte, tab, med = _apply_and_score(pack, Y, obs, site == te_s, types)
        tab = tab.assign(train_site=tr_s, test_site=te_s)
        per_tabs.append(tab)
        log(f"      median-over-types R²={med:+.3f}")
        # null: permute train ages (one site → within-site perm = any perm of train donors)
        log(f"      null ({N_PERM} train-age shuffles)...")
        null_med = []
        obs_tr = pack["obs"]
        for i in range(N_PERM):
            obs_p = permute_age_within_site_obs(obs_tr, rng)
            yp = obs_p.age.to_numpy(float)
            W_age, _, _ = per_type_directions(
                pack["X"], yp, pack["ct"], types, method="ridge", svd_cache=pack["svd"])
            W_tgt, _ = per_type_targets(W_age, pack["basis"])
            pack_p = dict(pack)
            pack_p["W_tgt"] = W_tgt
            pack_p["y"] = yp
            _, _, tnull, mnull = _apply_and_score(pack_p, Y, obs, site == te_s, types)
            null_med.append(mnull)
        pval = permutation_p(med, null_med, greater=True)
        ns = summarize_null_col(null_med)
        log(f"      null mean={ns['mean']:+.3f}  p={pval:.3f}")
        recs.append(dict(
            train_site=tr_s, test_site=te_s, median_r2=med,
            null_mean=ns["mean"], null_p95=ns["p95"], p=pval, n_perm=N_PERM,
            n_train_donors=int(obs.loc[site == tr_s, "donor"].nunique()),
            n_test_donors=int(obs.loc[site == te_s, "donor"].nunique()),
        ))
        pd.DataFrame(dict(perm=np.arange(N_PERM), median_r2=null_med)).to_csv(
            TARGET_DIR / f"v4a_{tr_s}_to_{te_s}_null.csv", index=False)
    per = pd.concat(per_tabs, ignore_index=True)
    per.to_csv(TARGET_DIR / "v4a_per_type.csv", index=False)
    rec_df = pd.DataFrame(recs)
    rec_df.to_csv(TARGET_DIR / "v4a_summary.csv", index=False)
    return rec_df, per


def _age_split(Y, obs, log):
    age = obs.age.to_numpy(float)
    don = obs.groupby("donor").agg(age=("age", "first"), site=("Source", "first"))
    # young: 20 <= age < 50; old: 50 <= age <= 89. Age==50 goes to old (no overlap).
    young_don = set(don.index[don.age < 50].astype(str))
    old_don = set(don.index[don.age >= 50].astype(str))
    donors = obs.donor.astype(str).to_numpy()
    m_y = np.array([d in young_don for d in donors])
    m_o = np.array([d in old_don for d in donors])
    log(f"[V4b] young (age<50): n_donors={len(young_don)}  rows={int(m_y.sum())}  "
        f"age {don.loc[don.age<50, 'age'].min():.0f}–{don.loc[don.age<50, 'age'].max():.0f}")
    log(f"      old   (age>=50): n_donors={len(old_don)}  rows={int(m_o.sum())}  "
        f"age {don.loc[don.age>=50, 'age'].min():.0f}–{don.loc[don.age>=50, 'age'].max():.0f}")
    log(f"      young sites={don.loc[don.age<50].groupby('site').size().to_dict()}  "
        f"old sites={don.loc[don.age>=50].groupby('site').size().to_dict()}")

    pack_y = _fit_on_mask(Y, obs, m_y)
    pack_o = _fit_on_mask(Y, obs, m_o)
    types = pack_y["types"]
    rows = []
    for j, t in enumerate(types):
        wy, wo = pack_y["W_tgt"][:, j], pack_o["W_tgt"][:, j]
        n_y = int((pack_y["ct"] == t).sum())
        n_o = int((pack_o["ct"] == t).sum())
        ang = unsigned_angle_deg(wy, wo)
        rows.append(dict(celltype=t, n_young=n_y, n_old=n_o, angle_deg=ang,
                         young_ok=bool(n_y >= MIN_N_TYPE and np.isfinite(wy).all()),
                         old_ok=bool(n_o >= MIN_N_TYPE and np.isfinite(wo).all())))
    tab = pd.DataFrame(rows)
    tab.to_csv(TARGET_DIR / "v4b_per_type_angles.csv", index=False)
    med = float(tab.loc[tab.young_ok & tab.old_ok, "angle_deg"].median())
    log(f"      median angle young-TARGET vs old-TARGET = {med:.1f}°")

    # consensus-consensus
    def _cons(pack):
        W = pack["W_tgt"]
        ok = [j for j in range(W.shape[1])
              if np.isfinite(W[:, j]).all() and np.linalg.norm(W[:, j]) > 1e-12]
        if len(ok) < 2:
            return None
        Wc = np.column_stack([W[:, j] / np.linalg.norm(W[:, j]) for j in ok])
        w, pc1 = consensus_direction(Wc)
        return w, pc1, len(ok)
    cy, co = _cons(pack_y), _cons(pack_o)
    ang_c = unsigned_angle_deg(cy[0], co[0]) if cy and co else np.nan
    log(f"      consensus young vs old angle = {ang_c:.1f}°  "
        f"(young PC1={cy[1] if cy else np.nan:.3f}, old PC1={co[1] if co else np.nan:.3f})")

    # does young-fit TARGET predict old ages (and reverse)? same-site mix — descriptive, with null later if cheap
    _, _, tab_yo, med_yo = _apply_and_score(pack_y, Y, obs, m_o, types)
    _, _, tab_oy, med_oy = _apply_and_score(pack_o, Y, obs, m_y, types)
    tab_yo.assign(direction="young_to_old").to_csv(TARGET_DIR / "v4b_young_to_old_per_type.csv", index=False)
    tab_oy.assign(direction="old_to_young").to_csv(TARGET_DIR / "v4b_old_to_young_per_type.csv", index=False)
    log(f"      young→old median R²={med_yo:+.3f}  old→young median R²={med_oy:+.3f}")
    return dict(
        n_donors_young=len(young_don), n_donors_old=len(old_don),
        median_angle_deg=med, consensus_angle_deg=float(ang_c) if np.isfinite(ang_c) else None,
        young_to_old_median_r2=float(med_yo), old_to_young_median_r2=float(med_oy),
        young_pc1=float(cy[1]) if cy else None, old_pc1=float(co[1]) if co else None,
        per_type=tab.to_dict(orient="records"),
    )


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
            # {lib: [[rank, term, p, z, combined, genes, adj, old_p, old_adj], ...]}
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


def _v4c(genes, W_tgt, w_cons, types, log):
    symbols = genes.symbol.astype(str).to_numpy()
    gids = genes.gene_id.astype(str).to_numpy()
    frames = []
    if w_cons is not None and np.isfinite(w_cons).all():
        frames.append(_top_genes(w_cons, symbols, gids).assign(source="consensus"))
    # mean of per-type targets (age-aligned already)
    acc = np.zeros(W_tgt.shape[0])
    n = 0
    for j in range(W_tgt.shape[1]):
        w = W_tgt[:, j]
        if np.isfinite(w).all() and float(np.linalg.norm(w)) > 1e-12:
            acc = acc + w / np.linalg.norm(w)
            n += 1
    w_mean = acc / max(n, 1)
    frames.append(_top_genes(w_mean, symbols, gids).assign(source="mean_per_type"))
    # two high-n interpretable types if present
    for want in ("oligodendrocyte", "astrocyte", "microglial cell",
                 "L2/3 intratelencephalic projecting glutamatergic neuron"):
        if want in types:
            j = types.index(want)
            frames.append(_top_genes(W_tgt[:, j], symbols, gids).assign(source=want))

    top = pd.concat(frames, ignore_index=True)
    top.to_csv(TARGET_DIR / "v4c_top_genes.csv", index=False)

    # technical composition of consensus top 50 ±
    cons = top[top.source == "consensus"]
    def _tech(sub):
        sy = sub.symbol.astype(str)
        return dict(
            n=int(len(sy)),
            n_mito=int(sy.map(is_mito).sum()),
            n_ribo=int(sy.map(is_ribo).sum()),
            n_cell_cycle=int(sy.isin(CELL_CYCLE_ALL).sum()),
            n_mt_prefix=int(sy.str.startswith("MT-").sum()),
        )
    tech = dict(
        consensus_pos=_tech(cons[cons.sign == "+"]),
        consensus_neg=_tech(cons[cons.sign == "-"]),
        mean_pos=_tech(top[(top.source == "mean_per_type") & (top.sign == "+")]),
        mean_neg=_tech(top[(top.source == "mean_per_type") & (top.sign == "-")]),
    )
    log(f"[V4c] consensus top50+ mito/ribo/cellcycle = "
        f"{tech['consensus_pos']['n_mito']}/{tech['consensus_pos']['n_ribo']}/{tech['consensus_pos']['n_cell_cycle']}")
    log(f"      consensus top50- mito/ribo/cellcycle = "
        f"{tech['consensus_neg']['n_mito']}/{tech['consensus_neg']['n_ribo']}/{tech['consensus_neg']['n_cell_cycle']}")

    hits = []
    for src in top.source.unique():
        for sgn, lab in (("+", "pos"), ("-", "neg")):
            sy = top[(top.source == src) & (top.sign == sgn)].symbol.astype(str).tolist()
            desc = f"{src}_{lab}"
            log(f"   Enrichr {desc} (n={len(sy)})")
            hits.extend(_enrichr(sy, desc, log))
    hit_df = pd.DataFrame(hits)
    if len(hit_df):
        hit_df.to_csv(TARGET_DIR / "v4c_enrichr.csv", index=False)
        aging = hit_df[hit_df.aging_token & (hit_df.p < 0.05)]
        log(f"      Enrichr terms with aging tokens & p<0.05: {len(aging)}")
        show = hit_df.sort_values("p").groupby("list").head(3)
        for _, row in show.iterrows():
            log(f"         {row.list} | {row.library}: {row.term}  p={row.p:.2e}")
    else:
        log("      Enrichr returned no terms (network or empty lists).")
        aging = pd.DataFrame()
    dump_json(TARGET_DIR / "v4c_tech_flags.json", tech)
    return dict(tech=tech, n_enrichr=int(len(hit_df)),
                n_aging_token_p05=int(len(aging)) if len(hit_df) else 0)


def run():
    rng = np.random.default_rng(TARGET_SEED)
    log = Logger(TARGET_DIR / "v4_report.txt")
    try:
        return _run(log, rng)
    finally:
        log.close()


def _run(log, rng):
    tgt_log_banner(log, "V4 — trustworthiness")
    v2stop = TARGET_DIR / "v2_STOP.txt"
    if v2stop.exists() and "stop_V2=True" in v2stop.read_text(encoding="utf-8"):
        log("[V4] STOP V2 is set — skipping trustworthiness of a non-predictive target.")
        dump_json(TARGET_DIR / "v4_summary.json", dict(skipped=True, reason="STOP V2"))
        return dict(skipped=True)

    data = load_phase1_matrix(log)
    Y, genes, obs = data["Y"], data["genes"], data["obs"]

    v4a_sum, v4a_pt = _site_transfer(Y, obs, rng, log)
    v4b = _age_split(Y, obs, log)

    pack = fit_full(Y, obs, method="ridge")
    cons_path = TARGET_DIR / "target_consensus.npz"
    w_cons = None
    if cons_path.exists():
        z = np.load(cons_path, allow_pickle=True)
        w_cons = z["consensus"]
    else:
        W = pack["W_tgt"]
        ok = [j for j in range(W.shape[1]) if np.isfinite(W[:, j]).all() and np.linalg.norm(W[:, j]) > 1e-12]
        if len(ok) >= 2:
            Wc = np.column_stack([W[:, j] / np.linalg.norm(W[:, j]) for j in ok])
            w_cons, _ = consensus_direction(Wc)
    v4c = _v4c(genes, pack["W_tgt"], w_cons, pack["types"], log)

    # figures
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    labs, vals, nulls = [], [], []
    for _, row in v4a_sum.iterrows():
        labs.append(f"{row.train_site}→{row.test_site}")
        vals.append(row.median_r2)
        nulls.append(row.null_mean)
    x = np.arange(len(labs))
    ax.bar(x, vals, color="tab:purple", label="observed")
    ax.scatter(x, nulls, c="k", zorder=3, label="null mean")
    ax.axhline(0, c="0.7", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(labs)
    ax.set_ylabel("median-over-types R²")
    ax.set_title("V4a  site transfer of TARGET")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(TARGET_FIG / "v4_site_transfer.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    tb = pd.DataFrame(v4b["per_type"]).dropna(subset=["angle_deg"]).sort_values("angle_deg")
    y_pos = np.arange(len(tb))
    ax.barh(y_pos, tb.angle_deg, color="tab:orange", height=0.7)
    ax.axvline(60, c="k", ls="--", lw=1, label="60°")
    ax.axvline(v4b["median_angle_deg"], c="tab:red", lw=1.5,
               label=f"median {v4b['median_angle_deg']:.1f}°")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(tb.celltype, fontsize=6)
    ax.set_xlabel("angle young TARGET vs old TARGET (deg)")
    ax.set_title("V4b  is it one line across the lifespan?")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(TARGET_FIG / "v4_age_range_angles.png", dpi=140)
    plt.close(fig)

    enr_path = TARGET_DIR / "v4c_enrichr.csv"
    if enr_path.exists():
        en = pd.read_csv(enr_path)
        sub = en[(en.list.str.contains("consensus")) & (en.p < 0.05)].sort_values("p").head(20)
        if len(sub):
            fig, ax = plt.subplots(figsize=(8.4, 5.2))
            y_pos = np.arange(len(sub))
            ax.barh(y_pos, -np.log10(np.clip(sub.p, 1e-30, 1)), color="tab:green", height=0.7)
            ax.set_yticks(y_pos)
            ax.set_yticklabels([f"{a} | {b}" for a, b in zip(sub.list, sub.term)], fontsize=6)
            ax.set_xlabel("-log10 p (Enrichr)")
            ax.set_title("V4c  consensus top-gene enrichment (p<0.05, top 20)")
            fig.tight_layout()
            fig.savefig(TARGET_FIG / "v4_enrichr.png", dpi=140)
            plt.close(fig)

    summary = dict(
        seed=TARGET_SEED, dataset_id=DATASET_ID,
        v4a=v4a_sum.to_dict(orient="records"),
        v4b=v4b, v4c=v4c,
    )
    dump_json(TARGET_DIR / "v4_summary.json", summary)
    log("[V4] wrote results/target/v4_*")
    log("[V4] done.")
    return summary


if __name__ == "__main__":
    run()

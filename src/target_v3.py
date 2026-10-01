"""STAGE V3 — is there ONE target, or twenty?

V3a. Pairwise angles between per-type TARGET vectors, with a null from permuted ages.
V3b. CONSENSUS target = PC1 of the per-type targets (joint within-type fit as companion).
     Report the fraction of variance across types it explains.
V3c. Does the consensus predict age within each cell type, held out?
     Compare per-type R² using the consensus vs that type's own target.

(a) One shared target — consensus performs near per-type targets.
(b) Type-specific targets — consensus loses most of the signal.

Usage: python src/target_v3.py
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
from target_common import (  # noqa: E402
    TARGET_DIR, TARGET_FIG, TARGET_SEED, DATASET_ID,
    N_PERM, N_PERM_PAIR, V3_RATIO_NEAR, V3_LOSS_NEAR, V3_RATIO_LOST, STOP_NULL,
    Logger, dump_json, load_json, tgt_log_banner, load_phase1_matrix, make_folds,
    fit_full, per_type_directions, per_type_targets, permute_age_within_site_obs,
    pairwise_median_deg, pairwise_angles_deg, consensus_direction,
    joint_within_type_direction, target_from_age, oof_consensus_vs_own,
    per_type_primary_r2, score_age_predictions, primary_r2, primary_mae,
    permutation_p, summarize_null_col, type_loses_age, unit, unsigned_angle_deg,
    zscore_train, identity_basis, save_matrix_bundle,
)


def _v3a_null(pack, obs, rng, n_perm, log):
    """Identity basis is age-independent; permute ages, refit ridge TARGETs (SVD cached)."""
    X, ct, types, y, basis = pack["X"], pack["ct"], pack["types"], pack["y"], pack["basis"]
    cache = pack["svd"]
    rows = []
    log(f"   V3a permutation null ({n_perm} donor-level within-site shuffles, ridge TARGET)")
    for i in range(n_perm):
        obs_p = permute_age_within_site_obs(obs, rng)
        yp = obs_p.age.to_numpy(float)
        W_age, _, _ = per_type_directions(X, yp, ct, types, method="ridge", svd_cache=cache)
        W_tgt, _ = per_type_targets(W_age, basis)
        med, _ = pairwise_median_deg(W_tgt)
        ok = [j for j in range(len(types))
              if np.isfinite(W_tgt[:, j]).all() and np.linalg.norm(W_tgt[:, j]) > 1e-12]
        if len(ok) >= 2:
            Wc = np.column_stack([W_tgt[:, j] / np.linalg.norm(W_tgt[:, j]) for j in ok])
            _, pc1 = consensus_direction(Wc)
        else:
            pc1 = np.nan
        rows.append(dict(perm=i, pairwise_median_deg=med, shared_pc1_frac=pc1))
        if (i + 1) % 10 == 0:
            log(f"      perm {i+1}/{n_perm}")
    return pd.DataFrame(rows)


def _align_consensus(W, w_cons):
    dots = []
    for j in range(W.shape[1]):
        w = W[:, j]
        if np.isfinite(w).all() and float(np.linalg.norm(w)) > 1e-12:
            dots.append(float(np.dot(w / np.linalg.norm(w), w_cons)))
    if dots and float(np.mean(dots)) < 0:
        return -w_cons
    return w_cons


def _v3c_scheme(Y, obs, folds, scheme, rng, log):
    log(f"   V3c OOF consensus vs own TARGET  scheme={scheme}")
    out = oof_consensus_vs_own(Y, obs, folds, method="ridge")
    cache = out.get("svd_cache")
    own_map, sc_own = per_type_primary_r2(obs, out["pred_own"], scheme)
    cons_map, sc_cons = per_type_primary_r2(obs, out["pred_cons"], scheme)
    prim_own = primary_r2(scheme, sc_own)
    prim_cons = primary_r2(scheme, sc_cons)
    log(f"      own TARGET primary R²={prim_own:+.3f}  consensus R²={prim_cons:+.3f}  "
        f"fold-mean PC1={out['pc1_frac']:.3f}")
    log(f"      consensus permutation null ({N_PERM} shuffles)...")
    null_rows = []
    for i in range(N_PERM):
        obs_p = permute_age_within_site_obs(obs, rng)
        yp = obs_p.age.to_numpy(float)
        o = oof_consensus_vs_own(Y, obs, folds, method="ridge", y_override=yp, svd_cache=cache)
        _, scn = per_type_primary_r2(obs_p, o["pred_cons"], scheme, y=yp)
        null_rows.append(dict(perm=i, primary=primary_r2(scheme, scn),
                              pooled=scn["pooled_median_type_r2"],
                              site_strat=scn["site_strat_median_type_r2"]))
        if (i + 1) % 5 == 0:
            log(f"         perm {i+1}/{N_PERM}")
    nt = pd.DataFrame(null_rows)
    p_cons = permutation_p(prim_cons, nt.primary, greater=True)
    null_sum = summarize_null_col(nt.primary.to_numpy(float))
    log(f"      consensus null mean={null_sum['mean']:+.3f}  p={p_cons:.3f}")

    rows = []
    for t in own_map:
        o, c = own_map[t], cons_map.get(t, np.nan)
        rows.append(dict(
            celltype=t, r2_own=o, r2_consensus=c,
            loss=float(o - c) if np.isfinite(o) and np.isfinite(c) else np.nan,
            ratio=float(c / o) if np.isfinite(o) and np.isfinite(c) and abs(o) > 1e-12 else np.nan,
        ))
    tab = pd.DataFrame(rows)
    tab.to_csv(TARGET_DIR / f"v3c_{scheme}_per_type.csv", index=False)
    sc_own["per_type"].assign(which="own").to_csv(
        TARGET_DIR / f"v3c_{scheme}_own_per_type.csv", index=False)
    sc_cons["per_type"].assign(which="consensus").to_csv(
        TARGET_DIR / f"v3c_{scheme}_cons_per_type.csv", index=False)
    nt.to_csv(TARGET_DIR / f"v3c_{scheme}_cons_null.csv", index=False)
    med_loss = float(tab.loss.median())
    med_ratio = float(tab.ratio.median())
    return dict(
        scheme=scheme,
        own_primary_r2=float(prim_own),
        cons_primary_r2=float(prim_cons),
        cons_null_mean=null_sum["mean"], cons_p=float(p_cons),
        cons_null_p95=null_sum["p95"], n_perm=N_PERM,
        median_loss=med_loss, median_ratio=med_ratio,
        fold_mean_pc1=float(out["pc1_frac"]),
        per_type=tab.to_dict(orient="records"),
    )


def _verdict(v3c_ws, pc1, pc1_null, pair_med, pair_null):
    own = v3c_ws["own_primary_r2"]
    cons = v3c_ws["cons_primary_r2"]
    loss = v3c_ws["median_loss"]
    ratio = (cons / own) if abs(own) > 1e-12 else np.nan
    near = (np.isfinite(cons) and np.isfinite(own)
            and cons >= V3_RATIO_NEAR * own and loss <= V3_LOSS_NEAR)
    lost = (not np.isfinite(cons)) or (np.isfinite(own) and cons < V3_RATIO_LOST * max(own, 1e-6))
    # also: consensus not above its null, while own is, counts as lost
    cons_above_null = cons > v3c_ws["cons_null_mean"] + 0.02
    if near:
        code, text = "a", (
            "One shared target — consensus performs near per-type targets. "
            "A single intervention is conceivable (as a linear move in expression space)."
        )
    elif lost or not cons_above_null:
        code, text = "b", (
            "Type-specific targets — consensus loses most of the signal. "
            "There is no universal rejuvenation direction in cortex; interventions would "
            "have to be cell-type-specific."
        )
    else:
        code, text = "partial", (
            "Partial sharing: consensus keeps some but not most of the per-type signal. "
            "A single cortical direction is not a substitute for type-specific targets."
        )
    return code, text, dict(ratio=ratio, loss=loss, near=near, lost=lost,
                            cons_above_null=cons_above_null)


def run():
    rng = np.random.default_rng(TARGET_SEED)
    log = Logger(TARGET_DIR / "v3_report.txt")
    try:
        return _run(log, rng)
    finally:
        log.close()


def _run(log, rng):
    tgt_log_banner(log, "V3 — one target or twenty")
    v2stop = TARGET_DIR / "v2_STOP.txt"
    if v2stop.exists() and "stop_V2=True" in v2stop.read_text(encoding="utf-8"):
        log("[V3] STOP V2 is set — refusing to interpret a TARGET that does not predict age.")
        dump_json(TARGET_DIR / "v3_summary.json", dict(skipped=True, reason="STOP V2"))
        return dict(skipped=True)

    data = load_phase1_matrix(log)
    Y, genes, obs = data["Y"], data["genes"], data["obs"]
    pack = fit_full(Y, obs, method="ridge")
    types = pack["types"]
    W = pack["W_tgt"]
    basis = pack["basis"]

    pair_med, ang = pairwise_median_deg(W)
    ok = [j for j in range(len(types))
          if np.isfinite(W[:, j]).all() and float(np.linalg.norm(W[:, j])) > 1e-12]
    Wc = np.column_stack([W[:, j] / np.linalg.norm(W[:, j]) for j in ok])
    w_cons, pc1 = consensus_direction(Wc)
    w_cons = _align_consensus(Wc, w_cons)
    log(f"[V3a] pairwise median angle of per-type TARGETs = {pair_med:.1f}°  (n_types_ok={len(ok)})")
    log(f"[V3b] PC1 of per-type TARGETs explains {pc1:.3f} of variance across types")

    if ang is not None:
        # ang is over `ok` types only
        ok_types = [types[j] for j in ok]
        pd.DataFrame(ang, index=ok_types, columns=ok_types).to_csv(
            TARGET_DIR / "v3a_pairwise_angles.csv")

    null = _v3a_null(pack, obs, rng, N_PERM_PAIR, log)
    null.to_csv(TARGET_DIR / "v3a_perm_null.csv", index=False)
    p_pair = permutation_p(pair_med, null.pairwise_median_deg, greater=False)
    p_pc1 = permutation_p(pc1, null.shared_pc1_frac, greater=True)
    pair_null = float(null.pairwise_median_deg.mean())
    pc1_null = float(null.shared_pc1_frac.mean())
    log(f"     pairwise null mean={pair_null:.1f}°  p(angle≤obs)={p_pair:.3f}")
    log(f"     PC1 null mean={pc1_null:.3f}  p={p_pc1:.3f}")

    w_joint, jrec = joint_within_type_direction(pack["X"], pack["y"], pack["ct"], types, method="ridge")
    w_joint_tgt, retained_j, ang_j, pr_j = target_from_age(w_joint, basis)
    # angles of joint TARGET to per-type TARGETs
    j_angs = [unsigned_angle_deg(w_joint_tgt, W[:, j]) for j in ok]
    log(f"[V3b] jointly-fit within-type ridge TARGET: angle to identity after projection "
        f"retained={retained_j:.3f}; median angle to per-type TARGETs={np.nanmedian(j_angs):.1f}°")

    # save consensus / joint for later
    save_matrix_bundle(
        TARGET_DIR / "target_consensus.npz",
        genes, types,
        W_target_ridge=W,
        consensus=w_cons,
        joint_within_type=w_joint_tgt,
        joint_raw=w_joint,
        basis=basis, mu=pack["mu"], sd=pack["sd"],
        pc1_frac=np.array([pc1]),
    )

    rng_folds = np.random.default_rng(TARGET_SEED)
    schemes = make_folds(obs, rng_folds)
    v3c = {}
    for name, folds in schemes.items():
        log(f"\n[V3c] scheme={name}")
        v3c[name] = _v3c_scheme(Y, obs, folds, name, rng, log)

    code, text, flags = _verdict(v3c["within_site"], pc1, pc1_null, pair_med, pair_null)
    log(f"\n[V3 verdict] ({code}) {text}")
    log(f"   within-site own R²={v3c['within_site']['own_primary_r2']:+.3f}  "
        f"consensus R²={v3c['within_site']['cons_primary_r2']:+.3f}  "
        f"median loss={v3c['within_site']['median_loss']:+.3f}  "
        f"median ratio={v3c['within_site']['median_ratio']:+.3f}")
    (TARGET_DIR / "v3_verdict.txt").write_text(
        f"verdict={code}\n{text}\n"
        f"pairwise_median_deg={pair_med}\nnull_pairwise={pair_null}\n"
        f"pc1={pc1}\nnull_pc1={pc1_null}\n"
        f"own_r2={v3c['within_site']['own_primary_r2']}\n"
        f"cons_r2={v3c['within_site']['cons_primary_r2']}\n",
        encoding="utf-8",
    )

    # figures
    if ang is not None:
        fig, ax = plt.subplots(figsize=(7.6, 6.6))
        im = ax.imshow(ang, cmap="viridis", vmin=0, vmax=90)
        short = [t if len(t) < 28 else t[:25] + "…" for t in ok_types]
        ax.set_xticks(range(len(ok_types)))
        ax.set_yticks(range(len(ok_types)))
        ax.set_xticklabels(short, rotation=90, fontsize=6)
        ax.set_yticklabels(short, fontsize=6)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="unsigned angle (deg)")
        ax.set_title(f"V3a  pairwise TARGET angles\nmedian {pair_med:.1f}°  (null {pair_null:.1f}°)")
        fig.tight_layout()
        fig.savefig(TARGET_FIG / "v3_pairwise.png", dpi=140)
        plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.8))
    ax = axes[0]
    ax.hist(null.pairwise_median_deg, bins=15, color="0.75", label="perm null")
    ax.axvline(pair_med, c="tab:purple", lw=2, label=f"obs {pair_med:.1f}°")
    ax.axvline(pair_null, c="k", ls=":")
    ax.set_xlabel("median pairwise TARGET angle (deg)")
    ax.set_title("V3a")
    ax.legend(fontsize=7)
    ax = axes[1]
    ax.hist(null.shared_pc1_frac, bins=15, color="0.75", label="perm null")
    ax.axvline(pc1, c="tab:purple", lw=2, label=f"obs {pc1:.3f}")
    ax.axvline(pc1_null, c="k", ls=":")
    ax.set_xlabel("PC1 fraction of per-type TARGETs")
    ax.set_title("V3b  shared?")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(TARGET_FIG / "v3_pairwise_pc1_vs_null.png", dpi=140)
    plt.close(fig)

    tab = pd.DataFrame(v3c["within_site"]["per_type"]).sort_values("r2_own")
    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    y_pos = np.arange(len(tab))
    ax.barh(y_pos - 0.18, tab.r2_own, height=0.36, color="tab:purple", label="own TARGET")
    ax.barh(y_pos + 0.18, tab.r2_consensus, height=0.36, color="tab:orange", label="consensus")
    ax.axvline(0, c="0.7", lw=0.6)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(tab.celltype, fontsize=6)
    ax.set_xlabel("site-stratified R² (within-site CV)")
    ax.set_title("V3c  consensus vs per-type TARGET")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(TARGET_FIG / "v3_consensus_vs_own.png", dpi=140)
    plt.close(fig)

    summary = dict(
        seed=TARGET_SEED, dataset_id=DATASET_ID,
        n_types_ok=len(ok), n_genes=int(Y.shape[1]),
        pairwise_median_deg=float(pair_med),
        pairwise_null_mean=pair_null, p_pairwise=float(p_pair),
        shared_pc1_frac=float(pc1), pc1_null_mean=pc1_null, p_pc1=float(p_pc1),
        n_perm_pair=N_PERM_PAIR,
        joint_median_angle_to_types=float(np.nanmedian(j_angs)),
        joint_norm_retained=float(retained_j),
        v3c=v3c, verdict=code, verdict_text=text, verdict_flags=flags,
    )
    dump_json(TARGET_DIR / "v3_summary.json", summary)
    log("[V3] wrote results/target/v3_*")
    log("[V3] done.")
    return summary


if __name__ == "__main__":
    run()

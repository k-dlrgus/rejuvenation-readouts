"""STAGE V1 — recover the age direction and identity subspace; bootstrap stability.

Per cell type, on the DLPFC cohort:
  - age direction (ridge; PLS-1 companion)
  - identity (type) subspace (centroid SVD over the 20 types)

Stability: donor bootstrap; distribution of unsigned angles between bootstrap
replicates of the age direction. A direction that is not stable under resampling
cannot be a target.

STOP V1: median-over-types of the per-type median pairwise bootstrap angle
above ~60°. Halt — the target vector would be noise.

Usage: python src/target_v1.py
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
    N_BOOT_RIDGE, N_BOOT_PLS, STOP_BOOT_ANGLE,
    Logger, dump_json, tgt_log_banner, load_phase1_matrix,
    fit_full, donor_bootstrap_index, per_type_directions, identity_basis,
    zscore_train, unsigned_angle_deg, pairwise_angles_deg, fmt, fmt_u,
)


def _bootstrap_angles(Y, obs, types, method, n_boot, rng, log, W_ref=None):
    """For each bootstrap: resample donors, global z-score, per-type age dirs.

    Returns per-type table + a long pairwise-angle frame.
    """
    ct = obs.celltype.astype(str).to_numpy()
    y = obs.age.to_numpy(float)
    p = Y.shape[1]
    # stack: boot x type x p would be huge; store list of W (p x T) is 20*25526*20*8 ~ 80MB
    Ws = []
    n_ok = np.zeros((n_boot, len(types)), dtype=bool)
    for b in range(n_boot):
        idx = donor_bootstrap_index(obs, rng)
        Xs, _, _ = zscore_train(Y[idx])
        W, meta, _ = per_type_directions(Xs, y[idx], ct[idx], types, method=method)
        Ws.append(W)
        n_ok[b] = meta.ok.to_numpy() if "ok" in meta.columns else np.isfinite(W).all(0)
        if (b + 1) % 5 == 0 or b == 0:
            log(f"   {method} bootstrap {b+1}/{n_boot}")
    # pairwise angles among replicates, per type
    pair_rows = []
    to_ref_rows = []
    summary = []
    for j, t in enumerate(types):
        cols = []
        for b in range(n_boot):
            w = Ws[b][:, j]
            if np.isfinite(w).all() and float(np.linalg.norm(w)) > 1e-12:
                cols.append(w / np.linalg.norm(w))
        if len(cols) < 2:
            summary.append(dict(
                celltype=t, method=method, n_boot_ok=len(cols), n_boot=n_boot,
                median_pairwise_deg=np.nan, mean_pairwise_deg=np.nan,
                p05_pairwise_deg=np.nan, p95_pairwise_deg=np.nan,
                median_to_full_deg=np.nan, n_pairs=0,
            ))
            continue
        Wc = np.column_stack(cols)
        ang = pairwise_angles_deg(Wc)
        iu = np.triu_indices(ang.shape[0], 1)
        vals = ang[iu]
        for a in vals:
            pair_rows.append(dict(celltype=t, method=method, angle_deg=float(a)))
        med_ref = np.nan
        if W_ref is not None:
            wr = W_ref[:, j]
            refs = [unsigned_angle_deg(Wc[:, k], wr) for k in range(Wc.shape[1])]
            med_ref = float(np.nanmedian(refs))
            for a in refs:
                to_ref_rows.append(dict(celltype=t, method=method, angle_deg=float(a)))
        summary.append(dict(
            celltype=t, method=method, n_boot_ok=int(Wc.shape[1]), n_boot=n_boot,
            median_pairwise_deg=float(np.median(vals)),
            mean_pairwise_deg=float(np.mean(vals)),
            p05_pairwise_deg=float(np.percentile(vals, 5)),
            p95_pairwise_deg=float(np.percentile(vals, 95)),
            median_to_full_deg=med_ref, n_pairs=int(len(vals)),
        ))
    return pd.DataFrame(summary), pd.DataFrame(pair_rows), pd.DataFrame(to_ref_rows)


def _compare_frozen(W_ridge, types, log):
    path = Path(__file__).resolve().parents[1] / "results" / "trajectory" / "frozen_scores.npz"
    if not path.exists():
        log("[V1] no results/trajectory/frozen_scores.npz — skip G1 freeze comparison")
        return None
    z = np.load(path, allow_pickle=True)
    fr_types = [str(t) for t in z["types"].tolist()]
    Wf = z["W_ridge"]
    rows = []
    for j, t in enumerate(types):
        if t not in fr_types:
            rows.append(dict(celltype=t, angle_to_g1_frozen_deg=np.nan))
            continue
        k = fr_types.index(t)
        rows.append(dict(celltype=t, angle_to_g1_frozen_deg=unsigned_angle_deg(W_ridge[:, j], Wf[:, k])))
    tab = pd.DataFrame(rows)
    log(f"[V1] median angle to G1 frozen ridge dirs: {tab.angle_to_g1_frozen_deg.median():.2f}° "
        f"(expected small; spaces differ by a second per-type z-score in G1)")
    return tab


def run():
    rng = np.random.default_rng(TARGET_SEED)
    log = Logger(TARGET_DIR / "v1_report.txt")
    try:
        return _run(log, rng)
    finally:
        log.close()


def _run(log, rng):
    tgt_log_banner(log, "V1 — recover axes + bootstrap stability")
    data = load_phase1_matrix(log)
    Y, genes, obs = data["Y"], data["genes"], data["obs"]
    n, p = Y.shape
    log(f"[filter] V1 uses ALL {p:,} genes. No threshold, no A/I bucket.")
    log(f"[filter] cohort: {data['n_donors']} donors, {data['n_types']} cell types, "
        f"sites={data['site_donor_counts']}, rows={n}")
    log(f"[filter] donor bootstrap: resample {data['n_donors']} donors with replacement; "
        f"all pseudobulks of each draw are kept (multiplicity). seed={TARGET_SEED}")
    log(f"[filter] n_boot ridge={N_BOOT_RIDGE}  n_boot pls1={N_BOOT_PLS}")

    pack_r = fit_full(Y, obs, method="ridge")
    pack_p = fit_full(Y, obs, method="pls1")
    types = pack_r["types"]
    log(f"[V1] identity subspace: {pack_r['basis'].shape[1]} axes from "
        f"{len(pack_r['used'])} type centroids; S head={np.round(pack_r['S_id'][:5], 3).tolist()}")

    meta_r = pack_r["meta"].copy()
    meta_r.to_csv(TARGET_DIR / "v1_ridge_full_meta.csv", index=False)
    pack_p["meta"].to_csv(TARGET_DIR / "v1_pls_full_meta.csv", index=False)
    log("[V1] full-cohort ridge age directions:")
    for _, row in meta_r.iterrows():
        log(f"   {row.celltype}: n={int(row.n)} ok={bool(row['ok'])} alpha={row.get('alpha', np.nan)} "
            f"norm={row.get('norm', np.nan)}")

    np.savez_compressed(
        TARGET_DIR / "v1_axes.npz",
        mu=pack_r["mu"], sd=pack_r["sd"],
        basis=pack_r["basis"], S_id=pack_r["S_id"],
        W_age_ridge=pack_r["W_age"], W_age_pls=pack_p["W_age"],
        W_tgt_ridge=pack_r["W_tgt"], W_tgt_pls=pack_p["W_tgt"],
        gene_id=genes.gene_id.astype(str).to_numpy(),
        symbol=genes.symbol.astype(str).to_numpy(),
        types=np.array(types, dtype=object),
        centroid_types=np.array(pack_r["used"], dtype=object),
    )
    log("[V1] wrote v1_axes.npz (age dirs + TARGET preview + identity basis)")

    frozen = _compare_frozen(pack_r["W_age"], types, log)
    if frozen is not None:
        frozen.to_csv(TARGET_DIR / "v1_vs_g1_frozen.csv", index=False)

    log(f"\n[V1] donor bootstrap of the RIDGE age direction ({N_BOOT_RIDGE} replicates)")
    sum_r, pair_r, ref_r = _bootstrap_angles(
        Y, obs, types, "ridge", N_BOOT_RIDGE, rng, log, W_ref=pack_r["W_age"])
    log(f"[V1] donor bootstrap of the PLS-1 age direction ({N_BOOT_PLS} replicates)")
    sum_p, pair_p, ref_p = _bootstrap_angles(
        Y, obs, types, "pls1", N_BOOT_PLS, rng, log, W_ref=pack_p["W_age"])

    boot_sum = pd.concat([sum_r, sum_p], ignore_index=True)
    boot_sum.to_csv(TARGET_DIR / "v1_bootstrap_summary.csv", index=False)
    pd.concat([pair_r, pair_p], ignore_index=True).to_csv(
        TARGET_DIR / "v1_bootstrap_pairwise.csv", index=False)
    pd.concat([ref_r, ref_p], ignore_index=True).to_csv(
        TARGET_DIR / "v1_bootstrap_to_full.csv", index=False)

    med_ridge = float(sum_r.median_pairwise_deg.median())
    med_pls = float(sum_p.median_pairwise_deg.median())
    n_unstable_r = int((sum_r.median_pairwise_deg > STOP_BOOT_ANGLE).sum())
    n_ok_r = int(sum_r.median_pairwise_deg.notna().sum())
    log(f"\n[V1] ridge median-over-types of median pairwise bootstrap angle = {med_ridge:.1f}°")
    log(f"     pls1  median-over-types of median pairwise bootstrap angle = {med_pls:.1f}°")
    log(f"     ridge types with median pairwise > {STOP_BOOT_ANGLE:.0f}°: {n_unstable_r}/{n_ok_r}")
    log("     per-type ridge:")
    for _, row in sum_r.sort_values("median_pairwise_deg").iterrows():
        log(f"       {row.celltype}: pairwise {row.median_pairwise_deg:.1f}°  "
            f"to-full {row.median_to_full_deg:.1f}°  n_boot={int(row.n_boot_ok)}")

    stop = bool(np.isfinite(med_ridge) and med_ridge > STOP_BOOT_ANGLE)
    (TARGET_DIR / "v1_STOP.txt").write_text(
        f"stop_V1={stop}\n"
        f"ridge_median_pairwise_deg={med_ridge}\n"
        f"pls_median_pairwise_deg={med_pls}\n"
        f"threshold_deg={STOP_BOOT_ANGLE}\n"
        f"n_types_above_threshold_ridge={n_unstable_r}/{n_ok_r}\n"
        f"reason={'median bootstrap angle between ridge replicates > ~60° — target would be noise' if stop else 'age direction stable enough to define a target'}\n",
        encoding="utf-8",
    )
    log(f"[STOP V1] {stop}  (median pairwise ridge {med_ridge:.1f}° vs {STOP_BOOT_ANGLE:.0f}°)")

    # figures
    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    sr = sum_r.sort_values("median_pairwise_deg")
    y_pos = np.arange(len(sr))
    ax.barh(y_pos, sr.median_pairwise_deg, color="tab:red", height=0.7,
            xerr=[sr.median_pairwise_deg - sr.p05_pairwise_deg,
                  sr.p95_pairwise_deg - sr.median_pairwise_deg],
            error_kw=dict(ecolor="0.4", lw=0.8, capsize=2))
    ax.axvline(STOP_BOOT_ANGLE, c="k", ls="--", lw=1, label=f"STOP {STOP_BOOT_ANGLE:.0f}°")
    ax.axvline(90, c="0.75", lw=0.6)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(sr.celltype, fontsize=6)
    ax.set_xlabel("median pairwise angle between ridge bootstrap replicates (deg)")
    ax.set_title("V1  age-direction stability (donor bootstrap)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(TARGET_FIG / "v1_bootstrap_per_type.png", dpi=140)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.8))
    ax = axes[0]
    if len(pair_r):
        ax.hist(pair_r.angle_deg, bins=30, color="tab:red", alpha=0.75)
    ax.axvline(med_ridge, c="k", lw=2, label=f"median-over-types {med_ridge:.1f}°")
    ax.axvline(STOP_BOOT_ANGLE, c="k", ls="--")
    ax.set_xlabel("pairwise bootstrap angle (deg)")
    ax.set_ylabel("pairs")
    ax.set_title("ridge")
    ax.legend(fontsize=7)
    ax = axes[1]
    if len(pair_p):
        ax.hist(pair_p.angle_deg, bins=30, color="tab:blue", alpha=0.75)
    ax.axvline(med_pls, c="k", lw=2, label=f"median-over-types {med_pls:.1f}°")
    ax.axvline(STOP_BOOT_ANGLE, c="k", ls="--")
    ax.set_xlabel("pairwise bootstrap angle (deg)")
    ax.set_title("PLS-1")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(TARGET_FIG / "v1_bootstrap_hist.png", dpi=140)
    plt.close(fig)

    summary = dict(
        seed=TARGET_SEED, dataset_id=DATASET_ID,
        n_donors=int(data["n_donors"]), n_types=len(types), n_genes=int(p),
        n_identity_axes=int(pack_r["basis"].shape[1]),
        n_boot_ridge=N_BOOT_RIDGE, n_boot_pls=N_BOOT_PLS,
        ridge_median_pairwise_deg=med_ridge,
        pls_median_pairwise_deg=med_pls,
        ridge_median_to_full_deg=float(sum_r.median_to_full_deg.median()),
        pls_median_to_full_deg=float(sum_p.median_to_full_deg.median()),
        n_types_unstable_ridge=n_unstable_r, n_types_ok_ridge=n_ok_r,
        stop_V1=stop, stop_threshold_deg=STOP_BOOT_ANGLE,
        site_donor_counts=data["site_donor_counts"],
        frozen_median_angle_deg=(float(frozen.angle_to_g1_frozen_deg.median())
                                 if frozen is not None else None),
        per_type_ridge=sum_r.to_dict(orient="records"),
        per_type_pls=sum_p.to_dict(orient="records"),
    )
    dump_json(TARGET_DIR / "v1_summary.json", summary)
    log("[V1] wrote results/target/v1_*")
    log("[V1] done.")
    return summary


if __name__ == "__main__":
    run()

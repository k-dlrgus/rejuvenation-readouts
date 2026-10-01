"""PART B — gene-count control (is the construct doomed by arithmetic?).

B1. 50 random gene sets of size |A|, predict age, both CV schemes, permutation null.
B2. Size sweep 10…1500 for A-subsets, I-subsets, and random, with a plot.
B3. Verdict (a) biology vs (b) method.

Does not derive new A/I definitions. Uses frozen P2 sets as the pools to subsample.

Usage: python src/geometry_partB.py
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
from geometry_common import (  # noqa: E402
    GEO_DIR, GEO_FIG, GEO_SEED, N_RANDOM_DRAWS, SIZE_SWEEP,
    NULL_USABLE, DATASET_ID, Logger, dump_json, load_json,
    geo_log_banner, load_phase1_matrix, make_folds,
    age_oof_svd, score_age_predictions, permute_age_within_site_obs,
)
N_SWEEP_DRAWS = 10  # B1 keeps 50; sweep uses 10 observed draws + one null per cell


def _metric(scheme, sc):
    if scheme == "within_site":
        return sc["site_strat_median_type_r2"], sc["pooled_median_type_r2"]
    return sc["pooled_median_type_r2"], sc["site_strat_median_type_r2"]


def _eval_set(Y, obs, idx, folds, scheme, y=None):
    if y is None:
        out = age_oof_svd(Y, obs, idx, folds)
        pred, yv = out["pred"], out["y"]
        sc = score_age_predictions(obs, pred, y=yv)
    else:
        pred = age_oof_svd(Y, obs, idx, folds, y_override=y)["pred"]
        sc = score_age_predictions(obs, pred, y=y)
    primary, other = _metric(scheme, sc)
    return dict(primary=primary, pooled=sc["pooled_median_type_r2"],
                site_strat=sc["site_strat_median_type_r2"],
                joint_pooled=sc["joint_pooled_r2"],
                joint_site_strat=sc["joint_site_strat_r2"])


def _draw_subset(pool, k, rng):
    pool = np.asarray(pool, dtype=int)
    k = int(min(k, pool.size))
    if k <= 0:
        return np.array([], dtype=int)
    if k == pool.size:
        return pool.copy()
    return rng.choice(pool, size=k, replace=False)


def run():
    rng = np.random.default_rng(GEO_SEED)
    log = Logger(GEO_DIR / "partB_report.txt")
    try:
        return _run(log, rng)
    finally:
        log.close()


def _run(log, rng):
    geo_log_banner(log, "PART B — gene-count control")
    a_sum = GEO_DIR / "partA_summary.json"
    if a_sum.exists():
        a = load_json(a_sum)
        if a.get("stop_A"):
            log("[B] STOP A is set — refusing to run Part B on an unsound pipeline.")
            dump_json(GEO_DIR / "partB_summary.json", dict(skipped=True, reason="STOP A"))
            return dict(skipped=True)
        log(f"[B] Part A passed. stop_A={a.get('stop_A')}  change={a.get('p4c_change')}")
    else:
        log("[B] WARNING: partA_summary.json missing; proceeding only if you ran A in-process.")

    data = load_phase1_matrix(log)
    Y, obs = data["Y"], data["obs"]
    A_idx, I_idx = data["A_idx"], data["I_idx"]
    n_genes = Y.shape[1]
    all_idx = np.arange(n_genes, dtype=int)
    n_A, n_I = int(len(A_idx)), int(len(I_idx))
    log(f"[filter] universe for random draws: all {n_genes:,} genes (no extra filter)")
    log(f"[filter] A-pool={n_A}  I-pool={n_I}  (frozen P2; not re-derived)")

    rng_folds = np.random.default_rng(GEO_SEED)
    schemes = make_folds(obs, rng_folds)

    # ------------------------------------------------------------------ B1
    log("\n" + "-" * 80)
    log(f"[B1] {N_RANDOM_DRAWS} random sets of size n_A={n_A} from all {n_genes} genes")
    b1_rows = []
    for name, folds in schemes.items():
        log(f"   scheme={name}")
        for d in range(N_RANDOM_DRAWS):
            sub = _draw_subset(all_idx, n_A, rng)
            sc = _eval_set(Y, obs, sub, folds, name)
            obs_p = permute_age_within_site_obs(obs, rng)
            scn = _eval_set(Y, obs_p, sub, folds, name, y=obs_p.age.to_numpy(float))
            b1_rows.append(dict(scheme=name, family="random", draw=d, k=n_A,
                                r2_primary=sc["primary"], r2_pooled=sc["pooled"],
                                r2_site_strat=sc["site_strat"],
                                null_primary=scn["primary"], null_pooled=scn["pooled"],
                                null_site_strat=scn["site_strat"]))
            if (d + 1) % 10 == 0:
                log(f"      draw {d+1}/{N_RANDOM_DRAWS}  r2={sc['primary']:+.3f}  "
                    f"null={scn['primary']:+.3f}")
        # frozen A and I as single points at this size
        scA = _eval_set(Y, obs, A_idx, folds, name)
        scI = _eval_set(Y, obs, I_idx, folds, name)
        obs_p = permute_age_within_site_obs(obs, rng)
        nA = _eval_set(Y, obs_p, A_idx, folds, name, y=obs_p.age.to_numpy(float))
        nI = _eval_set(Y, obs_p, I_idx, folds, name, y=obs_p.age.to_numpy(float))
        b1_rows.append(dict(scheme=name, family="A_frozen", draw=0, k=n_A,
                            r2_primary=scA["primary"], r2_pooled=scA["pooled"],
                            r2_site_strat=scA["site_strat"],
                            null_primary=nA["primary"], null_pooled=nA["pooled"],
                            null_site_strat=nA["site_strat"]))
        b1_rows.append(dict(scheme=name, family="I_frozen_full", draw=0, k=n_I,
                            r2_primary=scI["primary"], r2_pooled=scI["pooled"],
                            r2_site_strat=scI["site_strat"],
                            null_primary=nI["primary"], null_pooled=nI["pooled"],
                            null_site_strat=nI["site_strat"]))
        rand = [r for r in b1_rows if r["scheme"] == name and r["family"] == "random"]
        med = float(np.median([r["r2_primary"] for r in rand]))
        nmed = float(np.median([r["null_primary"] for r in rand]))
        log(f"   [{name}] random median R²={med:+.3f}  null median={nmed:+.3f}  "
            f"A_frozen={scA['primary']:+.3f}  I_full={scI['primary']:+.3f}")
    b1 = pd.DataFrame(b1_rows)
    b1.to_csv(GEO_DIR / "partB_b1_random_nA.csv", index=False)

    # ------------------------------------------------------------------ B2 size sweep
    log("\n" + "-" * 80)
    log(f"[B2] size sweep {list(SIZE_SWEEP)}  draws={N_SWEEP_DRAWS} (A/I random subsets + random genes)")
    b2_rows = []
    families = {
        "A": A_idx,
        "I": I_idx,
        "random": all_idx,
    }
    for name, folds in schemes.items():
        for fam, pool in families.items():
            for k in SIZE_SWEEP:
                k_eff = int(min(k, len(pool)))
                n_draw = 1 if k_eff == len(pool) and fam != "random" else N_SWEEP_DRAWS
                log(f"   {name}  {fam}  k={k} (eff={k_eff}) draws={n_draw}")
                obs_r2, null_r2 = [], []
                # One paired shuffle per (scheme, family, k), not per draw: the null is
                # about labels, not which genes were drawn. Every observed draw is still scored.
                sub0 = _draw_subset(pool, k_eff, rng)
                obs_p = permute_age_within_site_obs(obs, rng)
                scn0 = _eval_set(Y, obs_p, sub0, folds, name, y=obs_p.age.to_numpy(float))
                for d in range(n_draw):
                    sub = sub0 if d == 0 else _draw_subset(pool, k_eff, rng)
                    sc = _eval_set(Y, obs, sub, folds, name)
                    b2_rows.append(dict(
                        scheme=name, family=fam, k_requested=k, k_effective=k_eff, draw=d,
                        r2_primary=sc["primary"], r2_pooled=sc["pooled"],
                        r2_site_strat=sc["site_strat"],
                        null_primary=scn0["primary"], null_pooled=scn0["pooled"],
                        null_site_strat=scn0["site_strat"],
                    ))
                    obs_r2.append(sc["primary"])
                    null_r2.append(scn0["primary"])
                log(f"      median R²={float(np.median(obs_r2)):+.3f}  "
                    f"null={float(np.median(null_r2)):+.3f}")
    b2 = pd.DataFrame(b2_rows)
    b2.to_csv(GEO_DIR / "partB_b2_size_sweep.csv", index=False)

    # summary table at each k
    summ_rows = []
    for name in ("within_site", "loso"):
        for k in SIZE_SWEEP:
            for fam in ("A", "I", "random"):
                sub = b2[(b2.scheme == name) & (b2.family == fam) & (b2.k_requested == k)]
                if len(sub) == 0:
                    continue
                summ_rows.append(dict(
                    scheme=name, family=fam, k=k, k_effective=int(sub.k_effective.iloc[0]),
                    n_draws=len(sub),
                    r2_median=float(sub.r2_primary.median()),
                    r2_mean=float(sub.r2_primary.mean()),
                    r2_p05=float(sub.r2_primary.quantile(0.05)),
                    r2_p95=float(sub.r2_primary.quantile(0.95)),
                    null_median=float(sub.null_primary.median()),
                    null_mean=float(sub.null_primary.mean()),
                    null_usable=bool(abs(float(sub.null_primary.median())) <= NULL_USABLE
                                     or float(sub.null_primary.median()) <= NULL_USABLE),
                ))
    sweep_summ = pd.DataFrame(summ_rows)
    sweep_summ.to_csv(GEO_DIR / "partB_b2_size_sweep_summary.csv", index=False)

    # B1 summary
    b1s = []
    for name in ("within_site", "loso"):
        rand = b1[(b1.scheme == name) & (b1.family == "random")]
        Af = b1[(b1.scheme == name) & (b1.family == "A_frozen")]
        If = b1[(b1.scheme == name) & (b1.family == "I_frozen_full")]
        b1s.append(dict(
            scheme=name, n_A=n_A, n_draws=int(len(rand)),
            random_median=float(rand.r2_primary.median()),
            random_p05=float(rand.r2_primary.quantile(0.05)),
            random_p95=float(rand.r2_primary.quantile(0.95)),
            random_null_median=float(rand.null_primary.median()),
            A_frozen=float(Af.r2_primary.iloc[0]),
            A_null=float(Af.null_primary.iloc[0]),
            I_full=float(If.r2_primary.iloc[0]),
            I_null=float(If.null_primary.iloc[0]),
            A_beats_random=bool(float(Af.r2_primary.iloc[0]) > float(rand.r2_primary.quantile(0.95))),
            I_beats_A=bool(float(If.r2_primary.iloc[0]) > float(Af.r2_primary.iloc[0])),
            random_matches_A=bool(abs(float(Af.r2_primary.iloc[0]) - float(rand.r2_primary.median())) < 0.05),
        ))
    b1_summ = pd.DataFrame(b1s)
    b1_summ.to_csv(GEO_DIR / "partB_b1_summary.csv", index=False)
    log("\n[B1] summary:")
    log(b1_summ.to_string(index=False, float_format=lambda x: f"{x:+.3f}"))

    # ------------------------------------------------------------------ B3 verdict
    # Use within-site (corrected primary) as the informative scheme; LOSO must agree on the comparison.
    ws = b1_summ[b1_summ.scheme == "within_site"].iloc[0]
    lo = b1_summ[b1_summ.scheme == "loso"].iloc[0]
    # At matched size in the sweep, compare medians
    def _at(scheme, fam, k):
        sub = sweep_summ[(sweep_summ.scheme == scheme) & (sweep_summ.family == fam) & (sweep_summ.k == k)]
        return None if len(sub) == 0 else sub.iloc[0]

    # A vs random at k closest to n_A (1000 and 1500)
    k_match = 1500 if n_A >= 1500 else (1000 if n_A >= 1000 else 500)
    notes = []
    a_beats_rand_ws = False
    i_beats_a_ws = False
    random_matches_both = False

    # B1 is the direct matched-size test
    a_beats_rand_ws = bool(ws.A_frozen > ws.random_p95)
    # weaker: A above random median
    a_above_rand_med = bool(ws.A_frozen > ws.random_median + 0.05)
    i_beats_a_ws = bool(ws.I_full > ws.A_frozen + 0.02)
    rand_near_A = bool(abs(ws.A_frozen - ws.random_median) <= 0.05)
    rand_near_I = bool(abs(ws.I_full - ws.random_median) <= 0.08)

    if a_above_rand_med and i_beats_a_ws:
        verdict = "a"
        verdict_text = (
            "(a) A-genes beat random at matched size, and I-genes still beat A-genes → "
            "the fusion finding is about biology (the partition is informative, and identity "
            "genes still carry as much or more age signal)."
        )
    elif rand_near_A and (rand_near_I or abs(ws.I_full - ws.random_median) <= 0.10):
        verdict = "b"
        verdict_text = (
            "(b) Random matches both → the per-gene partition carries no information and "
            "the finding is about the method, not the tissue."
        )
    elif (not a_above_rand_med) and i_beats_a_ws:
        verdict = "b_leaning"
        verdict_text = (
            "(b-leaning) A-genes do not beat random at matched size, so the A-bucket is not "
            "a privileged age axis; I-genes still predict age, which is the methods failure "
            "mode (aggregating many genes with tiny age slopes) more than a unique biology of A-genes."
        )
    elif a_above_rand_med and not i_beats_a_ws:
        verdict = "a_incomplete"
        verdict_text = (
            "A-genes beat random at matched size, but I-genes do not beat A-genes under the "
            "corrected metric — the partition is not empty, and the cross-test does not "
            "show I ≥ A after leakage correction."
        )
    else:
        verdict = "mixed"
        verdict_text = (
            "Mixed: inspect the size-sweep. Neither clean (a) nor clean (b)."
        )
    notes.append(f"within-site B1: A={ws.A_frozen:+.3f} random_med={ws.random_median:+.3f} "
                 f"random_p95={ws.random_p95:+.3f} I={ws.I_full:+.3f} null_random={ws.random_null_median:+.3f}")
    notes.append(f"LOSO B1: A={lo.A_frozen:+.3f} random_med={lo.random_median:+.3f} "
                 f"I={lo.I_full:+.3f} null_random={lo.random_null_median:+.3f}")
    log("\n[B3] " + verdict_text)
    for n in notes:
        log("   " + n)

    # ------------------------------------------------------------------ figures
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.2), sharey=True)
    colors = dict(A="tab:red", I="tab:blue", random="0.45")
    for ax, name in zip(axes, ("within_site", "loso")):
        for fam in ("A", "I", "random"):
            sub = sweep_summ[(sweep_summ.scheme == name) & (sweep_summ.family == fam)].sort_values("k")
            if len(sub) == 0:
                continue
            ax.plot(sub.k, sub.r2_median, color=colors[fam], marker="o", lw=1.6, label=fam)
            ax.fill_between(sub.k, sub.r2_p05, sub.r2_p95, color=colors[fam], alpha=0.15, lw=0)
            ax.plot(sub.k, sub.null_median, color=colors[fam], ls=":", lw=1.0, alpha=0.8)
        ax.axhline(0, c="k", lw=0.5)
        ax.set_xlabel("set size (genes)")
        ax.set_title(name)
        ax.set_ylabel("age R² (primary metric)" if name == "within_site" else "")
        ax.legend(fontsize=8, loc="lower right")
        ylab = ("site-stratified median-over-types R²" if name == "within_site"
                else "pooled median-over-types R² (LOSO)")
        ax.set_ylabel(ylab, fontsize=8)
    fig.suptitle("B2  age R² vs gene-set size (dotted = permutation null)", fontsize=10)
    fig.tight_layout()
    fig.savefig(GEO_FIG / "partB_size_sweep.png", dpi=140)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.8), sharey=True)
    for ax, name in zip(axes, ("within_site", "loso")):
        rand = b1[(b1.scheme == name) & (b1.family == "random")]
        Af = b1[(b1.scheme == name) & (b1.family == "A_frozen")].iloc[0]
        If = b1[(b1.scheme == name) & (b1.family == "I_frozen_full")].iloc[0]
        ax.hist(rand.r2_primary, bins=15, color="0.7", edgecolor="k", linewidth=0.4,
                label="random size-matched")
        ax.axvline(Af.r2_primary, c="tab:red", lw=2, label=f"A-genes {Af.r2_primary:+.3f}")
        ax.axvline(If.r2_primary, c="tab:blue", lw=2, label=f"I-genes {If.r2_primary:+.3f}")
        ax.axvline(float(rand.null_primary.median()), c="k", ls=":", lw=1.2,
                   label=f"random null {float(rand.null_primary.median()):+.3f}")
        ax.set_title(name)
        ax.set_xlabel("age R² (primary)")
        ax.legend(fontsize=7)
    fig.suptitle(f"B1  {N_RANDOM_DRAWS} random sets of {n_A} genes", fontsize=10)
    fig.tight_layout()
    fig.savefig(GEO_FIG / "partB_random_nA.png", dpi=140)
    plt.close(fig)

    summary = dict(
        seed=GEO_SEED, dataset_id=DATASET_ID, n_A=n_A, n_I=n_I, n_genes=int(n_genes),
        n_random_draws=N_RANDOM_DRAWS, n_sweep_draws=N_SWEEP_DRAWS, size_sweep=list(SIZE_SWEEP),
        b1=[dict(r) for _, r in b1_summ.iterrows()],
        b3_verdict=verdict, b3_text=verdict_text, b3_notes=notes,
        k_match=k_match,
    )
    dump_json(GEO_DIR / "partB_summary.json", summary)
    (GEO_DIR / "partB_B3.txt").write_text(verdict + "\n" + verdict_text + "\n", encoding="utf-8")
    log("\n[B] wrote results/geometry/partB_*")
    log("[B] done.")
    return summary


if __name__ == "__main__":
    run()

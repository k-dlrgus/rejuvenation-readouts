"""STAGE P4 — shuffles + cross-test. Verdict per FALSIFICATION.md (unrevised).

Usage: python src/brain_phase1_p4.py
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
from brain_phase1_common import (  # noqa: E402
    PHASE1_DIR, PHASE1_FIG, PHASE1_SEED, Logger, HASH_POOL_LIMITATION,
    N_SITE_FOLDS, N_MATCH_DRAWS, ensure_logcpm, merge_p0_into_obs,
    load_json, dump_json, loso_folds, within_site_donor_folds,
)
from brain_phase1_models import (  # noqa: E402
    gene_index, age_oof_for_genes, identity_oof_for_genes, size_matched_I_age,
    permute_age_within_site_obs, permute_celltype,
)
from brain_phase1_p3 import _scheme_folds, _assert_no_donor_leak  # noqa: E402
from brain_phase1_findings import write_findings  # noqa: E402


def _verdict(r2_A, r2_I, acc_I, acc_A, shuffle_r2, shuffle_acc, chance):
    """Apply FALSIFICATION.md bands. Returns (label, text)."""
    notes = []
    sound = True
    if shuffle_r2 is not None and np.isfinite(shuffle_r2) and shuffle_r2 > 0.10:
        sound = False
        notes.append(f"shuffle-age R²={shuffle_r2:.3f} > 0.10 (pipeline soundness fail)")
    if shuffle_acc is not None and np.isfinite(shuffle_acc) and shuffle_acc > 2.0 * chance:
        sound = False
        notes.append(f"shuffle-type acc={shuffle_acc:.3f} > 2×chance (pipeline soundness fail)")

    ratio_r2 = (r2_I / r2_A) if (r2_A not in (0, None) and np.isfinite(r2_A) and r2_A != 0) else np.nan
    delta = (r2_A - r2_I) if (np.isfinite(r2_A) and np.isfinite(r2_I)) else np.nan
    ratio_acc = (acc_A / acc_I) if (acc_I not in (0, None) and np.isfinite(acc_I) and acc_I != 0) else np.nan

    falsified = False
    supported = False
    if np.isfinite(ratio_r2) and ratio_r2 >= 0.8:
        falsified = True
        notes.append(f"R²_I >= 0.8 × R²_A ({ratio_r2:.3f})")
    if np.isfinite(delta) and delta < 0.10:
        falsified = True
        notes.append(f"R²_A − R²_I = {delta:.3f} < 0.10")
    if np.isfinite(ratio_acc) and ratio_acc >= 0.9:
        falsified = True
        notes.append(f"Acc_A >= 0.9 × Acc_I ({ratio_acc:.3f})")

    if (np.isfinite(ratio_r2) and ratio_r2 <= 0.5 and np.isfinite(delta) and delta >= 0.15
            and np.isfinite(ratio_acc) and ratio_acc <= 0.7 and sound):
        supported = True

    if not sound:
        label = "soundness-fail"
    elif falsified:
        label = "falsified (not separable)"
    elif supported:
        label = "supported (separable)"
    else:
        label = "inconclusive"
    text = "; ".join(notes) if notes else (
        "within the inconclusive band" if label == "inconclusive" else label
    )
    return label, text, dict(ratio_r2=None if not np.isfinite(ratio_r2) else float(ratio_r2),
                             delta_r2=None if not np.isfinite(delta) else float(delta),
                             ratio_acc=None if not np.isfinite(ratio_acc) else float(ratio_acc))


def run_scheme(name, folds, Y, obs, A_idx, I_idx, rng, log):
    log("\n" + "=" * 80)
    log(f"[P4] scheme={name}")
    _assert_no_donor_leak(obs, folds, log, name)

    log("   age model on A-genes...")
    age_A = age_oof_for_genes(Y, obs, A_idx, folds, rng, method="ridge")
    log("   age model on full I-genes...")
    age_I_full = age_oof_for_genes(Y, obs, I_idx, folds, rng, method="ridge")
    log(f"   age model on size-matched I-genes ({N_MATCH_DRAWS} draws)...")
    matched = size_matched_I_age(Y, obs, I_idx, n_A=len(A_idx), folds=folds,
                                 rng=rng, n_draws=N_MATCH_DRAWS, method="ridge")
    log("   identity on I-genes...")
    id_I = identity_oof_for_genes(Y, obs, I_idx, folds, rng)
    log("   identity on A-genes (cross)...")
    id_A = identity_oof_for_genes(Y, obs, A_idx, folds, rng)

    obs_sh_age = permute_age_within_site_obs(obs, rng)
    sh_age = age_oof_for_genes(Y, obs_sh_age, A_idx, folds, rng, method="ridge")
    obs_sh_ct = permute_celltype(obs, rng)
    sh_id = identity_oof_for_genes(Y, obs_sh_ct, I_idx, folds, rng)

    r2_A = age_A["median_r2"]
    r2_I = matched["median_r2"]
    label, text, extra = _verdict(
        r2_A, r2_I, id_I["acc"], id_A["acc"],
        sh_age["median_r2"], sh_id["acc"], id_I["chance"],
    )
    log(f"   R²_A={r2_A:+.3f}  R²_I_matched={r2_I:+.3f}  R²_I_full={age_I_full['median_r2']:+.3f}")
    log(f"   Acc_I={id_I['acc']:.3f}  Acc_A={id_A['acc']:.3f}  chance={id_I['chance']:.3f}")
    log(f"   shuffle age R²={sh_age['median_r2']:+.3f}  shuffle type acc={sh_id['acc']:.3f}")
    log(f"   verdict: {label}  ({text})")
    matched["draws"].to_csv(PHASE1_DIR / f"p4_matched_I_draws_{name}.csv", index=False)
    return dict(
        r2_A=r2_A, r2_I_matched=r2_I, r2_I_full=age_I_full["median_r2"],
        acc_I=id_I["acc"], acc_A=id_A["acc"],
        shuffle_age_r2=sh_age["median_r2"], shuffle_id_acc=sh_id["acc"],
        id_chance=id_I["chance"],
        delta_r2=extra["delta_r2"], ratio_r2=extra["ratio_r2"], ratio_acc=extra["ratio_acc"],
        verdict=label, verdict_detail=text,
        joint_r2_A=age_A["joint_r2"], joint_r2_I_full=age_I_full["joint_r2"],
        n_A=int(len(A_idx)), n_I=int(len(I_idx)), n_matched=int(min(len(A_idx), len(I_idx))),
    )


def run():
    rng = np.random.default_rng(PHASE1_SEED)
    log = Logger(PHASE1_DIR / "p4_report.txt")
    try:
        return _run(log, rng)
    finally:
        log.close()


def _run(log, rng):
    log("=" * 100)
    log("BRAIN PHASE 1  STAGE P4 — negatives + cross-test.  seed =", PHASE1_SEED)
    log("=" * 100)
    log(HASH_POOL_LIMITATION)
    log("FALSIFICATION.md is pre-registration and is not revised.")
    p3 = load_json(PHASE1_DIR / "p3_summary.json")
    if p3.get("stop"):
        log("[P4] P3 STOPPED — not running cross-test")
        dump_json(PHASE1_DIR / "p4_summary.json", dict(skipped=True, reason="P3 stop"))
        write_findings(log=log)
        return dict(stop=True)

    frozen = load_json(PHASE1_DIR / "p2_frozen_gene_sets.json")
    Y, genes, obs = ensure_logcpm(log)
    obs = merge_p0_into_obs(obs)
    A_idx = gene_index(genes, frozen["A_genes"])
    I_idx = gene_index(genes, frozen["I_genes"])
    log(f"[P4] frozen A={len(A_idx)} I={len(I_idx)} (full-data sets; models retrained in CV, genes not re-derived)")

    schemes = _scheme_folds(obs, rng)
    rec = {}
    for name, folds in schemes.items():
        rec[name] = run_scheme(name, folds, Y, obs, A_idx, I_idx, rng, log)

    labels = {k: rec[k]["verdict"] for k in rec}
    if len(set(labels.values())) == 1:
        overall = list(labels.values())[0]
        overall_text = (
            f"Both CV schemes agree: **{overall}**. "
            + rec["within_site"]["verdict_detail"]
        )
    else:
        overall = "DISAGREEMENT between CV schemes"
        overall_text = (
            f"LOSO verdict = {rec['loso']['verdict']}; "
            f"within-site verdict = {rec['within_site']['verdict']}. "
            "Hard rule: report the disagreement rather than picking the better number. "
            f"LOSO: {rec['loso']['verdict_detail']}. Within-site: {rec['within_site']['verdict_detail']}."
        )
    log("\n[P4] THE DECISION: " + overall)
    log(overall_text)

    # figure
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.8), sharey=True)
    for ax, name in zip(axes, ("within_site", "loso")):
        s = rec[name]
        labs = ["A-genes", "I matched", "I full", "shuffle age"]
        vals = [s["r2_A"], s["r2_I_matched"], s["r2_I_full"], s["shuffle_age_r2"]]
        cols = ["tab:red", "tab:blue", "tab:cyan", "0.6"]
        ax.bar(labs, vals, color=cols)
        ax.axhline(0, c="k", lw=0.5)
        ax.set_title(f"{name}: {s['verdict']}")
        ax.set_ylabel("median-over-types age R²")
        ax.tick_params(axis="x", rotation=20, labelsize=8)
    fig.tight_layout()
    fig.savefig(PHASE1_FIG / "p4_cross_test_age.png", dpi=140)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.6), sharey=True)
    for ax, name in zip(axes, ("within_site", "loso")):
        s = rec[name]
        labs = ["I-genes", "A-genes", "shuffle type", "chance"]
        vals = [s["acc_I"], s["acc_A"], s["shuffle_id_acc"], s["id_chance"]]
        ax.bar(labs, vals, color=["tab:blue", "tab:red", "0.6", "0.85"])
        ax.set_title(name)
        ax.set_ylabel("identity accuracy")
        ax.tick_params(axis="x", rotation=20, labelsize=8)
    fig.tight_layout()
    fig.savefig(PHASE1_FIG / "p4_cross_test_identity.png", dpi=140)
    plt.close(fig)

    dump_json(PHASE1_DIR / "p4_summary.json", dict(
        loso=rec["loso"], within_site=rec["within_site"],
        verdict=overall, verdict_text=overall_text, seed=PHASE1_SEED,
        n_match_draws=N_MATCH_DRAWS,
    ))
    (PHASE1_DIR / "p4_verdict.txt").write_text(overall + "\n" + overall_text + "\n", encoding="utf-8")
    write_findings(log=log)
    log("[P4] done.")
    return dict(verdict=overall)


if __name__ == "__main__":
    run()

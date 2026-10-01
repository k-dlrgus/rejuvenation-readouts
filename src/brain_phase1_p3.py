"""STAGE P3 — A-gene age score, I-gene identity score, unrestricted ceiling.

Both CV schemes. Nested gene selection inside each outer training fold.
STOP P3: identity near chance, or A-gene age R² < ~0.15 while unrestricted > 0.35.

Usage: python src/brain_phase1_p3.py
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
    N_SITE_FOLDS, STOP_P3_AGE_R2, STOP_P3_UNREST_R2, IDENTITY_NEAR_CHANCE_MULT,
    ensure_logcpm, merge_p0_into_obs, load_p0_extras, load_json, dump_json,
    loso_folds, within_site_donor_folds, bootstrap_r2_mae,
)
from brain_phase1_models import (  # noqa: E402
    select_genes_on_train, age_oof_for_genes, identity_oof_for_genes, load_chrom_series,
)
from brain_phase1_findings import write_findings  # noqa: E402


def _scheme_folds(obs, rng):
    return {
        "loso": [(tr, te, site) for tr, te, site in loso_folds(obs)],
        "within_site": within_site_donor_folds(obs, N_SITE_FOLDS, rng),
    }


def _assert_no_donor_leak(obs, folds, log, name):
    for fold in folds:
        tr, te = fold[0], fold[1]
        ov = set(obs.iloc[tr].donor) & set(obs.iloc[te].donor)
        assert not ov, f"{name} donor leak {sorted(ov)[:5]}"
    log(f"[P3] {name}: {len(folds)} folds, no donor in both folds: True")


def run_scheme(name, folds, Y, obs, genes, extra, chrom_s, rng, log):
    log("\n" + "=" * 80)
    log(f"[P3] scheme={name}  n_folds={len(folds)}")
    _assert_no_donor_leak(obs, folds, log, name)
    infos = []
    for i, fold in enumerate(folds):
        info = select_genes_on_train(Y, obs, genes, extra, fold[0], chrom_s=chrom_s)
        infos.append(info)
        extra_tag = fold[2] if len(fold) > 2 else i
        log(f"   fold {i} {extra_tag}: nested A={info['n_A']} I={info['n_I']}")
    A_list = [i["A_idx"] for i in infos]
    I_list = [i["I_idx"] for i in infos]
    all_idx = np.arange(Y.shape[1])

    log("   fitting A-gene ridge (within type)...")
    age_A = age_oof_for_genes(Y, obs, A_list, folds, rng, method="ridge")
    log("   fitting unrestricted kernel ridge (within type)...")
    age_U = age_oof_for_genes(Y, obs, all_idx, folds, rng, method="kernel")
    log("   fitting I-gene identity classifier...")
    id_I = identity_oof_for_genes(Y, obs, I_list, folds, rng)

    pt = age_A["per_type"].merge(age_U["per_type"], on="celltype", suffixes=("_A", "_unrest"))
    pt["scheme"] = name
    log(f"   median-over-types R² A-genes={age_A['median_r2']:+.3f}  unrest={age_U['median_r2']:+.3f}  "
        f"joint A={age_A['joint_r2']:+.3f} unrest={age_U['joint_r2']:+.3f}")
    log(f"   identity acc={id_I['acc']:.3f}  macro-F1={id_I['macro_f1']:.3f}  chance={id_I['chance']:.3f}  "
        f"mean P(true)={id_I['mean_proba_true']:.3f}")

    # CIs on donor-level joint A
    don = age_A["donor_pred"]
    ci = bootstrap_r2_mae(don.y.to_numpy(), don.pred.to_numpy(), rng)
    summ = dict(
        age_A_median_r2=age_A["median_r2"], age_A_median_mae=age_A["median_mae"],
        age_A_joint_r2=age_A["joint_r2"], age_A_joint_mae=age_A["joint_mae"],
        age_unrest_median_r2=age_U["median_r2"], age_unrest_median_mae=age_U["median_mae"],
        age_unrest_joint_r2=age_U["joint_r2"],
        age_A_joint_r2_lo=ci["r2_lo"], age_A_joint_r2_hi=ci["r2_hi"],
        id_I_acc=id_I["acc"], id_I_f1=id_I["macro_f1"], id_chance=id_I["chance"],
        id_mean_proba_true=id_I["mean_proba_true"], id_mean_dist_own=id_I["mean_dist_own"],
        median_n_A_nested=float(np.median([i["n_A"] for i in infos])),
        median_n_I_nested=float(np.median([i["n_I"] for i in infos])),
        n_folds=len(folds),
    )
    return dict(summary=summ, per_type=pt, age_A=age_A, age_U=age_U, id_I=id_I, infos=infos)


def run():
    rng = np.random.default_rng(PHASE1_SEED)
    log = Logger(PHASE1_DIR / "p3_report.txt")
    try:
        return _run(log, rng)
    finally:
        log.close()


def _run(log, rng):
    log("=" * 100)
    log("BRAIN PHASE 1  STAGE P3 — scores.  seed =", PHASE1_SEED)
    log("=" * 100)
    log(HASH_POOL_LIMITATION)
    if not (PHASE1_DIR / "p1_DECLARED_BEFORE_SCORES.flag").exists():
        raise RuntimeError("P1 thresholds not declared before scores")
    p2 = load_json(PHASE1_DIR / "p2_summary.json")
    if p2.get("stop"):
        log("[P3] P2 STOPPED — not building scores")
        return dict(stop=True, stop_reason="P2 stop")
    frozen = load_json(PHASE1_DIR / "p2_frozen_gene_sets.json")
    log(f"[P3] frozen full-data sets (for listing/P5, not nested CV features): "
        f"A={frozen['n_A']} I={frozen['n_I']}")

    Y, genes, obs = ensure_logcpm(log)
    obs = merge_p0_into_obs(obs)
    extra = load_p0_extras()
    chrom_s = load_chrom_series(genes.gene_id.to_numpy())

    schemes = _scheme_folds(obs, rng)
    out_summ = {}
    per_type_rows = []
    id_rows = []
    stop = False
    stop_reason = ""
    interpretation = ""

    for name, folds in schemes.items():
        r = run_scheme(name, folds, Y, obs, genes, extra, chrom_s, rng, log)
        out_summ[name] = r["summary"]
        per_type_rows.append(r["per_type"])
        id_rows.append(dict(scheme=name, **{k: r["id_I"][k] for k in
                                            ("acc", "macro_f1", "chance", "mean_proba_true", "mean_dist_own", "n")}))
        # persist OOF
        pd.DataFrame({"donor": obs.donor, "celltype": obs.celltype, "Source": obs.Source,
                      "age": obs.age, "pred_A": r["age_A"]["pred"], "pred_unrest": r["age_U"]["pred"]}
                     ).to_csv(PHASE1_DIR / f"p3_oof_age_{name}.csv", index=False)

    pt = pd.concat(per_type_rows, ignore_index=True)
    pt.to_csv(PHASE1_DIR / "p3_age_per_type.csv", index=False)
    pd.DataFrame(id_rows).to_csv(PHASE1_DIR / "p3_identity.csv", index=False)

    # STOP uses within-site (informative scheme); LOSO is still reported
    ws = out_summ["within_site"]
    unrest = ws["age_unrest_median_r2"]
    a_r2 = ws["age_A_median_r2"]
    acc = ws["id_I_acc"]
    chance = ws["id_chance"]
    if acc < IDENTITY_NEAR_CHANCE_MULT * chance:
        stop = True
        stop_reason = (
            f"STOP P3: identity accuracy {acc:.3f} is near chance ({chance:.3f}, "
            f"threshold {IDENTITY_NEAR_CHANCE_MULT}×chance)."
        )
    if (a_r2 < STOP_P3_AGE_R2) and (unrest > STOP_P3_UNREST_R2):
        stop = True
        stop_reason = (
            (stop_reason + " " if stop_reason else "") +
            f"STOP P3: A-gene median-over-types R²={a_r2:.3f} < {STOP_P3_AGE_R2} while "
            f"unrestricted={unrest:.3f} > {STOP_P3_UNREST_R2} — gene-set restriction discards real signal."
        )
    if not stop:
        if abs(a_r2 - unrest) < 0.10 or a_r2 >= 0.70 * unrest:
            interpretation = (
                f"A-gene construct captures the within-type age signal "
                f"(A {a_r2:.3f} vs unrestricted {unrest:.3f} median-over-types, within-site). Proceed."
            )
        else:
            interpretation = (
                f"A-gene R²={a_r2:.3f} vs unrestricted {unrest:.3f} (within-site median-over-types). "
                "Not a STOP, but the restriction is lossy. Proceed to P4 with that caveat."
            )
        log("[P3] " + interpretation)
    else:
        log("[P3] " + stop_reason)
        interpretation = stop_reason

    # figure
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.0))
    for ax, name in zip(axes, ("within_site", "loso")):
        sub = pt[pt.scheme == name].sort_values("r2_unrest")
        y = np.arange(len(sub))
        ax.barh(y - 0.18, sub.r2_unrest, height=0.35, color="0.65", label="unrestricted")
        ax.barh(y + 0.18, sub.r2_A, height=0.35, color="tab:red", label="A-genes")
        ax.set_yticks(y)
        ax.set_yticklabels(sub.celltype, fontsize=6)
        ax.axvline(0, c="k", lw=0.5)
        ax.set_xlabel("held-out R²")
        ax.set_title(name)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(PHASE1_FIG / "p3_age_A_vs_unrestricted.png", dpi=140)
    plt.close(fig)

    dump_json(PHASE1_DIR / "p3_summary.json", dict(
        loso=out_summ["loso"], within_site=out_summ["within_site"],
        stop=stop, stop_reason=stop_reason, interpretation=interpretation,
        per_type_table=True, seed=PHASE1_SEED,
        frozen_n_A=frozen["n_A"], frozen_n_I=frozen["n_I"],
    ))
    (PHASE1_DIR / "p3_STOP_verdict.txt").write_text(
        ("STOP\n" + stop_reason) if stop else ("PASS\n" + interpretation), encoding="utf-8"
    )
    write_findings(log=log)
    log("[P3] done.")
    return dict(stop=stop)


if __name__ == "__main__":
    run()

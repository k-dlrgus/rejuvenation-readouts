"""STAGE P2 — confound filters on the declared A-gene set.

Same lists and fixed thresholds as blood 1b (Tirosh + data-driven cell cycle, ribo, mito,
sex-linked + sex DE, depth), plus genes correlated with included P0 covariates.

STOP P2: fewer than 30 A-genes survive. Halt (do not build scores).

Usage: python src/brain_phase1_p2.py
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
    STOP_P2_MIN_A, ensure_logcpm, merge_p0_into_obs, load_p0_extras,
    load_json, dump_json,
)
from brain_phase1_findings import write_findings  # noqa: E402
from confounds import build_flags, covariate_age_report, residualize_on_celltype, corr_with_covariate  # noqa: E402
from ensembl_chrom import fetch_chromosomes  # noqa: E402
from gene_lists import CELL_CYCLE_ALL  # noqa: E402

FILTER_STEPS = [
    ("cell-cycle (Tirosh S/G2M + proliferation list)", "flag_cc_list"),
    ("cell-cycle (data: |r|>=0.3 with cycling covariates)", "flag_cc_data"),
    ("ribosomal (RPL/RPS/MRP)", "flag_ribo"),
    ("mitochondrial (MT-)", "flag_mito"),
    ("sex-linked (chrX/chrY)", "flag_sex_chrom"),
    ("sex-DE (|d|>=0.5, p<1e-6, within cell type)", "flag_sex_data_only"),
    ("depth-correlated (|r|>=0.3 with mean UMI/genes per nucleus)", "flag_depth"),
    ("P0-covariate-correlated (|r|>=0.3, within type)", "flag_p0"),
]


def survival_table(genes_A, F):
    fa = F.loc[genes_A].copy()
    fa["flag_sex_data_only"] = fa.flag_sex_data & ~fa.flag_sex_chrom
    if "flag_p0" not in fa:
        fa["flag_p0"] = False
    n0 = len(fa)
    rows, removed = [], pd.Series(False, index=fa.index)
    for name, col in FILTER_STEPS:
        if col not in fa.columns:
            flag = pd.Series(False, index=fa.index)
        else:
            flag = fa[col].astype(bool)
        rows.append(dict(
            filter=name,
            flagged_marginal=int(flag.sum()),
            removed_incremental=int((flag & ~removed).sum()),
            surviving=n0 - int((removed | flag).sum()),
        ))
        removed |= flag
    return pd.DataFrame(rows), removed


def add_p0_flags(F, Y, obs, extra_cols, log):
    F = F.copy()
    F["flag_p0"] = False
    if not extra_cols:
        log("[P2] no P0 extras — flag_p0 is empty")
        return F
    ct = obs.celltype.to_numpy()
    Yr = residualize_on_celltype(np.asarray(Y, float), ct)
    any_flag = np.zeros(Y.shape[1], dtype=bool)
    for c in extra_cols:
        if c not in obs.columns:
            log(f"[P2] P0 extra {c} not on obs")
            continue
        v = obs[c].to_numpy(float)
        v = residualize_on_celltype(v[:, None], ct)[:, 0]
        r = corr_with_covariate(Yr, v)
        F[f"r_p0_{c}"] = r
        hit = np.abs(r) >= 0.30
        n_hit = int(np.nansum(hit))
        log(f"[P2] |r|>=0.3 with {c}: {n_hit} genes")
        any_flag |= np.nan_to_num(hit, nan=False).astype(bool)
    F["flag_p0"] = any_flag
    return F


def run():
    log = Logger(PHASE1_DIR / "p2_report.txt")
    try:
        return _run(log)
    finally:
        log.close()


def _run(log):
    log("=" * 100)
    log("BRAIN PHASE 1  STAGE P2 — confound filters.  seed =", PHASE1_SEED)
    log("=" * 100)
    log(HASH_POOL_LIMITATION)
    if not (PHASE1_DIR / "p1_DECLARED_BEFORE_SCORES.flag").exists():
        raise RuntimeError("P1 thresholds were not declared before scores — run p1 first")
    thresh = load_json(PHASE1_DIR / "p1_thresholds.json")
    decomp = pd.read_csv(PHASE1_DIR / "p1_variance_decomposition.csv", index_col=0)
    A_in = decomp.index[decomp.gene_class == "A"].tolist()
    I_in = decomp.index[decomp.gene_class == "I"].tolist()
    log(f"[P2] declared A-genes in: {len(A_in)}; I-genes: {len(I_in)}; MIXED: {int((decomp.gene_class=='MIXED').sum())}")

    Y, genes, obs = ensure_logcpm(log)
    obs = merge_p0_into_obs(obs)
    extra = load_p0_extras()
    gid = genes.gene_id.to_numpy()
    sym = genes.symbol.to_numpy()
    # chromosomes
    gid0 = pd.Index([g.split(".")[0] for g in gid])
    try:
        chrom = fetch_chromosomes(gid, verbose=True)
        chrom_s = pd.Series(chrom.chrom.astype(str).to_numpy(), index=chrom.gene_id.astype(str))
        chrom_map = chrom_s.copy()
        chrom_map.index = [i.split(".")[0] for i in chrom_map.index]
        chrom_aligned = chrom_map.reindex(gid0).fillna("NA")
        chrom_aligned.index = gid
        log(f"[P2] chromosomes mapped for {int((chrom_aligned != 'NA').sum()):,}/{len(gid):,} genes")
    except Exception as e:
        log(f"[P2] Ensembl chromosome fetch failed ({e}); sex-chrom flag will use symbol list + cache only")
        chrom_aligned = pd.Series("NA", index=gid)

    F = build_flags(Y, obs, gid, sym, chrom_aligned, log=log)
    F = add_p0_flags(F, Y, obs, extra, log)
    F["gene_class"] = decomp.reindex(F.index).gene_class
    F.to_csv(PHASE1_DIR / "p2_gene_flags.csv")

    log("\nCovariate–age (within cell type) for depth/cycle (orientation):")
    try:
        covariate_age_report(obs, log=log)
    except Exception as e:
        log(f"[P2] covariate_age_report skipped: {e}")

    surv, removed = survival_table(A_in, F)
    log("\nSequential A-gene survival:")
    log(surv.to_string(index=False))
    surv.to_csv(PHASE1_DIR / "p2_Agene_survival.csv", index=False)

    surviving = [g for g in A_in if not bool(removed.loc[g])]
    log(f"\n[P2] A surviving: {len(surviving)} / {len(A_in)}")
    # cell-cycle investigation
    fa = F.loc[A_in]
    n_cc = int((fa.flag_cc_list | fa.flag_cc_data).sum())
    cc_note = (
        f"{n_cc}/{len(A_in)} declared A-genes are cell-cycle flagged "
        f"(list={int(fa.flag_cc_list.sum())}, data={int(fa.flag_cc_data.sum())}). "
        "DLPFC is largely post-mitotic; OPCs are the main cycling population. "
    )
    if n_cc >= max(10, 0.15 * max(len(A_in), 1)):
        cc_syms = fa.loc[fa.flag_cc_list | fa.flag_cc_data, "symbol"].tolist()
        cc_note += f"Flagged symbols (n={len(cc_syms)}): " + ", ".join(cc_syms[:40])
        cc_note += " — a large CC hit in post-mitotic tissue is informative (OPC/progenitor leak or list overlap), not ignored."
        log("[P2] LARGE cell-cycle removal: " + cc_note)
    else:
        cc_note += "Small CC hit, as expected for post-mitotic cortex."
        log("[P2] " + cc_note)

    A_df = decomp.loc[surviving].sort_values("unique_age", ascending=False)
    A_df.to_csv(PHASE1_DIR / "p2_A_genes_surviving.csv")
    decomp.loc[I_in].sort_values("unique_ct", ascending=False).to_csv(PHASE1_DIR / "p2_I_genes.csv")
    frozen = dict(
        A_genes=list(map(str, surviving)),
        I_genes=list(map(str, I_in)),
        MIXED_genes=list(map(str, decomp.index[decomp.gene_class == "MIXED"])),
        A_symbols=A_df.symbol.astype(str).tolist(),
        n_A=len(surviving), n_I=len(I_in),
        seed=PHASE1_SEED, thresholds=thresh,
    )
    dump_json(PHASE1_DIR / "p2_frozen_gene_sets.json", frozen)

    fig, ax = plt.subplots(figsize=(8.2, 3.6))
    ax.plot(np.arange(len(surv)), surv.surviving, marker="o", c="tab:red")
    ax.set_xticks(np.arange(len(surv)))
    ax.set_xticklabels([s[:28] for s in surv["filter"]], rotation=35, ha="right", fontsize=7)
    ax.axhline(STOP_P2_MIN_A, c="k", ls="--", lw=0.8, label=f"STOP line {STOP_P2_MIN_A}")
    ax.set_ylabel("A-genes surviving")
    ax.set_title("P2 sequential survival of declared A-genes")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(PHASE1_FIG / "p2_Agene_survival.png", dpi=140)
    plt.close(fig)

    stop = len(surviving) < STOP_P2_MIN_A
    stop_reason = ""
    if stop:
        stop_reason = (
            f"STOP P2: {len(surviving)} A-genes survive filters (need >={STOP_P2_MIN_A}). "
            "With 233 donors and R3 expression R²=0.44 this means the per-gene unique_age construct "
            "is not recovering the multivariate age signal. Halt; do not build scores."
        )
        log("\n" + stop_reason)
        (PHASE1_DIR / "p2_STOP_verdict.txt").write_text(stop_reason, encoding="utf-8")
    else:
        log(f"\n[P2] {len(surviving)} A-genes survive (>= {STOP_P2_MIN_A}). Proceed to P3.")
        (PHASE1_DIR / "p2_STOP_verdict.txt").write_text("PASS\n", encoding="utf-8")

    dump_json(PHASE1_DIR / "p2_summary.json", dict(
        n_A_in=len(A_in), n_A_surviving=len(surviving), n_I=len(I_in),
        stop=stop, stop_reason=stop_reason, cc_note=cc_note,
        n_cc_flagged=n_cc, extra_covariates=extra, seed=PHASE1_SEED,
    ))
    write_findings(log=log)
    log("[P2] done.")
    return dict(stop=stop, n_A=len(surviving))


if __name__ == "__main__":
    run()

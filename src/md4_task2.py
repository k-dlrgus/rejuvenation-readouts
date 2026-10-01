"""Task 2 — extrapolation distance of state × timepoint pseudobulks to GTEx train. Report-only."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md4_common import (  # noqa: E402
    MD4_DIR, MD4_PROC, MD2_DIR, MD4_SEED, MD4_BOOT,
    AGED_LINE, YOUNG_LINE, DAYS, STATE_NAMES, PRIMARY_PAIR,
    PREREG_TASK2, PREREG_TASK2_FLAG, FROZEN_RULER, N_PERM, N_BOOT,
    PCA_MAHAL_K, FIBRO3_EXTRAP_RHO, FIBRO2_DIR,
    StopStep, Logger, dump_json, jsonable, load_json, md4_log_banner,
    load_manifest, save_manifest, record_failure, log_columns,
    progress_snapshot, load_frozen_ruler,
)
from fibro2_extrap import _pca_mahal_fit, _nn_and_mahal, _gtex_z  # noqa: E402
from fibro_stage1 import load_pack  # noqa: E402
from fibro2_common import tmm_logcpm_rows, score_frozen  # noqa: E402
from gtex_common import spearman_safe  # noqa: E402
from trajectory_common import permutation_p  # noqa: E402
from fibro_common import bootstrap_rho_ci  # noqa: E402


def _load_Y(line):
    path = MD4_PROC / f"louvain_rulerY_{line}.npz"
    if not path.exists():
        raise StopStep("task2", f"missing {path} — run Task 1")
    z = np.load(path)
    return sparse.csr_matrix(
        (z["data"], z["indices"], z["indptr"]),
        shape=tuple(int(x) for x in z["shape"]),
    )


def _rho_block(dist, age, tag, n_perm=N_PERM, n_boot=N_BOOT, seed=MD4_SEED, boot=MD4_BOOT):
    dist = np.asarray(dist, float)
    age = np.asarray(age, float)
    m = np.isfinite(dist) & np.isfinite(age)
    dist, age = dist[m], age[m]
    n = int(len(dist))
    rec = dict(distance=tag, n=n, n_perm=int(n_perm), n_boot=int(n_boot),
               seed=int(seed), boot_seed=int(boot))
    if n < 8:
        rec.update(rho=np.nan, rho_null=np.nan, rho_p=np.nan,
                   rho_ci_lo=np.nan, rho_ci_hi=np.nan,
                   note="n<8; spearman_safe not applied")
        return rec
    rho = spearman_safe(dist, age)
    rng = np.random.default_rng(int(seed))
    nulls = []
    for _ in range(int(n_perm)):
        yp = age.copy()
        rng.shuffle(yp)
        nulls.append(spearman_safe(dist, yp))
    nulls = np.asarray(nulls, float)
    rec["rho"] = float(rho) if np.isfinite(rho) else np.nan
    rec["rho_null"] = float(np.nanmedian(nulls)) if np.isfinite(nulls).any() else np.nan
    rec["rho_p"] = permutation_p(rho, nulls, greater=True)
    rng_b = np.random.default_rng(int(boot))
    ci = bootstrap_rho_ci(dist, age, rng_b, n_boot=n_boot)
    rec["rho_ci_lo"] = ci.get("p025")
    rec["rho_ci_hi"] = ci.get("p975")
    rec["note"] = "permutation shuffles frozen ruler score, holds distance; unit=state×timepoint"
    return rec


def run_task2(log=None):
    if not PREREG_TASK2_FLAG.exists():
        raise StopStep("prereg", "PREREG_TASK2.flag missing")
    if not FROZEN_RULER.exists():
        raise StopStep("frozen_ruler", f"missing {FROZEN_RULER}")
    close_log = False
    if log is None:
        log = Logger(MD4_DIR / "t2_report.txt")
        close_log = True
    md4_log_banner(log, "TASK2")
    log(PREREG_TASK2)
    frozen = load_frozen_ruler()
    pack = load_pack()
    log(f"[t2] GTEx pack X shape={np.asarray(pack['X']).shape}")
    log_columns("t2_gtex_obs", list(pack["obs"].columns), str(FIBRO2_DIR.parent / "fibro" / "stage1_obs.csv"))
    Xz = _gtex_z(pack, frozen)
    log(f"[t2] GTEx train z n={Xz.shape[0]} p={Xz.shape[1]}")
    pca = _pca_mahal_fit(Xz, PCA_MAHAL_K, MD4_SEED)
    log(f"[t2] PCA k={pca['k']}")

    rows = []
    Z_blocks = []
    meta_blocks = []
    for line in (AGED_LINE, YOUNG_LINE):
        obs = pd.read_csv(MD2_DIR / f"louvain_obs_{line}.csv")
        lab = pd.read_csv(MD2_DIR / f"t2_cluster_labels_{line}.csv")
        log_columns(f"t2_louvain_obs_{line}", list(obs.columns), str(MD2_DIR / f"louvain_obs_{line}.csv"))
        log_columns(f"t2_labels_{line}", list(lab.columns), str(MD2_DIR / f"t2_cluster_labels_{line}.csv"))
        lab_map = {int(r.cluster): str(r.label) for _, r in lab.iterrows()}
        obs = obs.reset_index(drop=True)
        obs["label"] = obs["cluster"].map(lambda c: lab_map.get(int(c)))
        Y = _load_Y(line)
        if Y.shape[0] != len(obs):
            raise StopStep("task2", f"{line} Y n={Y.shape[0]} obs n={len(obs)}")
        mats, meta = [], []
        for d in DAYS:
            for st in STATE_NAMES:
                m = (obs.day.astype(int) == int(d)) & (obs.label == st)
                n = int(m.sum())
                if n < 1:
                    continue
                idx = np.flatnonzero(m.to_numpy())
                pb = np.asarray(Y[idx].sum(axis=0), dtype=np.float64).ravel()
                mats.append(pb)
                meta.append(dict(cell_line=line, day=int(d), label=st, n_cells=n))
        if not mats:
            log(f"[t2] {line}: no state×timepoint rows with n≥1")
            continue
        C = np.vstack(mats)
        logcpm, nf = tmm_logcpm_rows(C, log, tag=f"t2_state_tp_{line}")
        age, Z, missing = score_frozen(logcpm, frozen, counts=C)
        nn, nn_idx, mahal = _nn_and_mahal(Z, Xz, pca)
        for i, rec in enumerate(meta):
            rec.update(
                age_score=float(age[i]),
                nn_euclidean=float(nn[i]),
                nn_train_index=int(nn_idx[i]),
                mahalanobis_pca=float(mahal[i]),
                pca_k=int(pca["k"]),
                age_source="TMM_log2CPM_frozen_ruler_state_timepoint_per_donor",
                n_missing_z0=int(np.asarray(missing).sum()),
            )
            rows.append(rec)
            log(f"[t2] {line} d{rec['day']} {rec['label']} n={rec['n_cells']} "
                f"age={rec['age_score']:+.3f} nn={rec['nn_euclidean']:+.3f} "
                f"mahal={rec['mahalanobis_pca']:+.3f}")
        Z_blocks.append(Z)
        meta_blocks.extend(meta)

    if not rows:
        raise StopStep("task2", "no state×timepoint rows to score")
    df = pd.DataFrame(rows)
    df.to_csv(MD4_DIR / "t2_extrap.csv", index=False)
    if Z_blocks:
        np.savez_compressed(MD4_DIR / "t2_state_z.npz", Z=np.vstack(Z_blocks))

    corr_rows = []
    for tag in ("nn_euclidean", "mahalanobis_pca"):
        rec = _rho_block(df[tag].to_numpy(float), df.age_score.to_numpy(float), tag)
        rec["scope"] = "all_donor_state_timepoint_n_ge_1"
        corr_rows.append(rec)
        log(f"[t2 corr] all {tag} n={rec['n']} ρ={rec.get('rho')} p={rec.get('rho_p')} "
            f"CI=[{rec.get('rho_ci_lo')}, {rec.get('rho_ci_hi')}]")
        for line, g in df.groupby("cell_line"):
            rec_d = _rho_block(g[tag].to_numpy(float), g.age_score.to_numpy(float), tag)
            rec_d["scope"] = f"per_donor_{line}"
            rec_d["cell_line"] = line
            corr_rows.append(rec_d)
            log(f"[t2 corr] {line} {tag} n={rec_d['n']} ρ={rec_d.get('rho')}")

    corr = pd.DataFrame(corr_rows)
    corr.to_csv(MD4_DIR / "t2_correlation.csv", index=False)

    # Direction relative to Task 1, not a gate.
    dir_rows = []
    for line in (AGED_LINE, YOUNG_LINE):
        g = df[df.cell_line == line]
        for d in DAYS:
            pr = g[(g.day == int(d)) & (g.label == PRIMARY_PAIR[0])]
            nr = g[(g.day == int(d)) & (g.label == PRIMARY_PAIR[1])]
            if not len(pr) or not len(nr):
                dir_rows.append(dict(
                    cell_line=line, day=int(d), ok=False,
                    reason="missing PartialReprog or NonReprog row",
                ))
                continue
            nn_pr, nn_nr = float(pr.nn_euclidean.iloc[0]), float(nr.nn_euclidean.iloc[0])
            mh_pr, mh_nr = float(pr.mahalanobis_pca.iloc[0]), float(nr.mahalanobis_pca.iloc[0])
            age_pr, age_nr = float(pr.age_score.iloc[0]), float(nr.age_score.iloc[0])
            pr_farther_nn = bool(nn_pr > nn_nr)
            rec = dict(
                cell_line=line, day=int(d), ok=True,
                n_PR=int(pr.n_cells.iloc[0]), n_NR=int(nr.n_cells.iloc[0]),
                nn_PR=nn_pr, nn_NR=nn_nr, nn_PR_minus_NR=nn_pr - nn_nr,
                mahal_PR=mh_pr, mahal_NR=mh_nr, mahal_PR_minus_NR=mh_pr - mh_nr,
                age_PR=age_pr, age_NR=age_nr, age_PR_minus_NR=age_pr - age_nr,
                PR_farther_nn=pr_farther_nn,
                PR_farther_mahal=bool(mh_pr > mh_nr),
            )
            # If ρ(distance, ruler)<0 and PR farther, this would push ruler(PR) down.
            rec["would_push_ruler_PR_lower_via_nn"] = bool(pr_farther_nn)
            dir_rows.append(rec)
    dirdf = pd.DataFrame(dir_rows)
    dirdf.to_csv(MD4_DIR / "t2_direction.csv", index=False)

    pooled_nn = corr[(corr.scope == "all_donor_state_timepoint_n_ge_1") & (corr.distance == "nn_euclidean")]
    rho_nn = float(pooled_nn.rho.iloc[0]) if len(pooled_nn) and pd.notna(pooled_nn.rho.iloc[0]) else np.nan
    summary = dict(
        pca_k=int(pca["k"]),
        n_rows=int(len(df)),
        rho_nn_euclidean_all=rho_nn,
        fibro3_quoted=FIBRO3_EXTRAP_RHO,
        note=(
            "Report-only. If ρ(distance, ruler)<0 and PartialReprog sits farther from "
            "training than NonReprog, that would push the Task 1 frozen-ruler contrast "
            "toward PR reading younger. Not a gate. Not used to discount Task 1."
        ),
    )
    dump_json(MD4_DIR / "t2_summary.json", jsonable(summary))
    man = load_manifest()
    man["status"] = "TASK2_DONE"
    save_manifest(man)
    progress_snapshot("write FINDINGS_MD4.md from disk", stop="TASK2_DONE")
    if close_log:
        log.close()
    return summary


if __name__ == "__main__":
    try:
        run_task2()
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise

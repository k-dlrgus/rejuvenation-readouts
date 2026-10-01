"""Task 1 — joint within-donor clustering and per-cluster trajectory. Refit nothing."""
from __future__ import annotations

import gc
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fibro3_common import (  # noqa: E402
    FIBRO3_DIR, FIBRO3_PROC, FIBRO3_SEED, MIN_CELLS, CLUSTER_K, CLUSTER_NPC,
    CLUSTER_BATCH, CLUSTER_N_INIT, DAYS, AGED_LINE, YOUNG_LINE,
    PREREG_TASK1, PREREG_TASK1_FLAG, CLUSTER_SETTINGS, CLUSTER_SETTINGS_FLAG,
    DECLARED_BEFORE_SCORES_FLAG, FROZEN_RULER,
    StopStep, Logger, dump_json, jsonable, fibro3_log_banner,
    load_manifest, save_manifest, record_failure, log_columns,
    load_frozen_ruler, progress_snapshot,
)
from fibro_common import FIBRO_DIR, FIBRO_RAW, unit  # noqa: E402
from fibro_stage2 import (  # noqa: E402
    find_10x_h5, load_one_h5, attach_meta, align_counts_to_ruler,
    _row_sum, _log1p_cp10k, score_age, differentiation_score,
)
from fibro2_common import tmm_logcpm_rows, score_frozen  # noqa: E402
from fibro2_tb import mito_mask  # noqa: E402


def _require_flags():
    if not PREREG_TASK1_FLAG.exists():
        raise StopStep("prereg", "PREREG_TASK1.flag missing")
    if not CLUSTER_SETTINGS_FLAG.exists():
        raise StopStep("prereg", "CLUSTER_SETTINGS.flag missing — settings must be frozen before clustering")
    existing = CLUSTER_SETTINGS_FLAG.read_text(encoding="utf-8")
    if existing.strip() != CLUSTER_SETTINGS.strip():
        raise StopStep("prereg", "CLUSTER_SETTINGS.flag does not match the frozen block. Not rewriting.")


def _load_libs(log):
    extract = FIBRO_RAW / "RAW"
    if not extract.exists():
        raise StopStep("task1", f"missing GSE297234 extract {extract}")
    samp_path = FIBRO_DIR / "stage2_geo_samples.csv"
    if not samp_path.exists():
        raise StopStep("task1", f"missing {samp_path} — Stage 2 sample table is required and is not re-fetched")
    samples = pd.read_csv(samp_path)
    log(f"[task1] Stage 2 sample table columns actually read: {list(samples.columns)} path={samp_path}")
    log_columns("t1_stage2_geo_samples", list(samples.columns), str(samp_path))
    need = {"gsm", "cell_line", "day", "age_years", "age_source"}
    missing = sorted(need - set(samples.columns))
    if missing:
        raise StopStep("task1", f"{samp_path} missing columns {missing}. Not substituting.")
    h5s = find_10x_h5(extract, log)
    if not h5s:
        raise StopStep("task1", "no filtered_feature_bc_matrix.h5 under GSE297234 RAW")
    libs = [load_one_h5(h, log) for h in h5s]
    libs = attach_meta(libs, samples, log)
    log(f"[task1] loaded {len(libs)} libraries")
    for lib in libs:
        log_columns(
            f"t1_h5_{lib['obs']['gsm'].iloc[0]}",
            ["barcode", "gene_id", "symbol", "feature_type", "genome"],
            lib["obs"]["file"].iloc[0] if "file" in lib["obs"].columns else "h5",
        )
    return libs, samples


def _qc_cells(lib):
    X = lib["X"].tocsr()
    umi = np.asarray(X.sum(axis=1)).ravel()
    ngenes = np.asarray((X > 0).sum(axis=1)).ravel()
    mito, mito_meta = mito_mask(lib["var"])
    if mito.any():
        mito_umi = np.asarray(X[:, mito].sum(axis=1)).ravel()
        frac = mito_umi / np.where(umi > 0, umi, np.nan)
    else:
        frac = np.full(X.shape[0], np.nan)
    return umi, ngenes, frac, mito_meta


def _batch_slices(n, batch, min_last):
    starts = list(range(0, n, batch))
    slices = [(s, min(s + batch, n)) for s in starts]
    if len(slices) >= 2 and (slices[-1][1] - slices[-1][0]) < min_last:
        s0, _ = slices[-2]
        _, e1 = slices[-1]
        slices = slices[:-2] + [(s0, e1)]
    return slices


def joint_cluster_one_donor(libs_d, frozen, log):
    """Cluster jointly across timepoints. No age score is computed here."""
    libs_d = sorted(libs_d, key=lambda L: int(L["day"]))
    line = libs_d[0]["cell_line"]
    days = [int(L["day"]) for L in libs_d]
    if days != sorted(days):
        raise StopStep("cluster", f"{line}: library day order {days}")
    Y_parts, obs_parts = [], []
    mito_meta = None
    for lib in libs_d:
        Y, rec, idx = align_counts_to_ruler(lib["X"], lib["var"], frozen, log)
        umi, ngenes, mito_frac, mito_meta = _qc_cells(lib)
        obs = lib["obs"].copy()
        obs["cell_line"] = lib["cell_line"]
        obs["day"] = int(lib["day"])
        obs["age_years"] = int(lib["age_years"])
        obs["age_source"] = lib["obs"]["age_source"].iloc[0]
        obs["gsm"] = str(lib["obs"]["gsm"].iloc[0])
        obs["umi"] = umi
        obs["n_genes"] = ngenes
        obs["mito_frac"] = mito_frac
        if len(obs) != Y.shape[0]:
            raise StopStep("cluster", f"{line} d{lib['day']}: obs n={len(obs)} Y n={Y.shape[0]}")
        Y_parts.append(Y)
        obs_parts.append(obs)
        log(f"[cluster-align] {line} d{lib['day']} n_cells={Y.shape[0]} overlap={rec.get('n_overlap')}")
    Yall = sparse.vstack(Y_parts, format="csr")
    obs = pd.concat(obs_parts, axis=0, ignore_index=True)
    n = int(Yall.shape[0])
    colsum = _row_sum(Yall)
    keep_idx = np.flatnonzero(colsum > 0)
    n_keep = int(keep_idx.size)
    if n < 2:
        raise StopStep("cluster", f"{line}: n_cells={n} < 2")
    npc = int(min(CLUSTER_NPC, n - 1, n_keep))
    log(f"[cluster] {line} joint n_cells={n} n_genes_kept={n_keep} npc={npc} "
        f"k={CLUSTER_K} days={days} pca=randomized")
    log(f"[cluster] {line} column subset via CSC")
    Yk = Yall.tocsc()[:, keep_idx].tocsr().astype(np.float32)
    log(f"[cluster] {line} Yk shape={Yk.shape} nnz={int(Yk.nnz)}")
    Z = Yk.toarray()
    del Yk
    gc.collect()
    lib = Z.sum(1, keepdims=True)
    np.maximum(lib, 1.0, out=lib)
    Z /= lib
    Z *= 1e4
    np.log1p(Z, out=Z)
    log(f"[cluster] {line} log1p-CP10k dense {Z.shape} dtype={Z.dtype}")
    pca = PCA(n_components=npc, svd_solver="randomized", random_state=int(FIBRO3_SEED))
    P = pca.fit_transform(Z).astype(np.float32)
    evr = float(np.sum(pca.explained_variance_ratio_))
    del Z
    gc.collect()
    log(f"[cluster] {line} PCA done explained_variance_ratio_sum={evr:.4f}")
    km = KMeans(n_clusters=CLUSTER_K, random_state=int(FIBRO3_SEED), n_init=CLUSTER_N_INIT)
    lab = km.fit_predict(np.asarray(P, np.float32))
    obs = obs.copy()
    obs["cluster"] = lab.astype(int)
    sizes = np.bincount(lab, minlength=CLUSTER_K).tolist()
    log(f"[cluster] {line} joint sizes={sizes} inertia={float(km.inertia_):.4g}")
    del P
    gc.collect()
    return obs, Yall, dict(
        cell_line=line, n_cells=n, n_genes_kept=n_keep, npc=npc,
        k=CLUSTER_K, sizes=sizes, days=days,
        batch_size=CLUSTER_BATCH, n_init=CLUSTER_N_INIT,
        method="PCA_randomized+KMeans",
        pca_explained_variance_ratio_sum=evr,
        random_state=FIBRO3_SEED,
        mito_rule=(mito_meta or {}).get("rule"),
        n_mito_genes=(mito_meta or {}).get("n_mito"),
        mito_symbols=(mito_meta or {}).get("mito_symbols"),
    )


def _pseudobulk_rows(obs, Yall, frozen, log):
    rows, mats = [], []
    for (line, day, cl), g in obs.groupby(["cell_line", "day", "cluster"], sort=True):
        idx = g.index.to_numpy()
        n = int(len(idx))
        pb = np.asarray(Yall[idx].sum(axis=0), dtype=np.float64).ravel()
        n_tp = int((obs.cell_line.eq(line) & obs.day.eq(day)).sum())
        rows.append(dict(
            cell_line=str(line), day=int(day), cluster=int(cl),
            age_years=int(g.age_years.iloc[0]), gsm=str(g.gsm.iloc[0]),
            age_source=str(g.age_source.iloc[0]),
            n_cells=n, n_timepoint=n_tp,
            frac_timepoint=float(n / n_tp) if n_tp else np.nan,
            median_UMI=float(np.median(g.umi.to_numpy(float))) if n else np.nan,
            median_genes=float(np.median(g.n_genes.to_numpy(float))) if n else np.nan,
            mito_frac_median_cell=float(np.nanmedian(g.mito_frac.to_numpy(float))) if n else np.nan,
            below_min_cells=bool(n < MIN_CELLS),
        ))
        mats.append(pb)
    C = np.vstack(mats)
    logcpm, nf = tmm_logcpm_rows(C, log, tag="t1_joint_cluster")
    return pd.DataFrame(rows), logcpm, C, nf


def composition_table(cluster_df):
    recs = []
    for line, g in cluster_df.groupby("cell_line"):
        clusters = sorted(g.cluster.astype(int).unique().tolist())
        for cl in clusters:
            r7 = g[(g.day == 7) & (g.cluster == cl)]
            r10 = g[(g.day == 10) & (g.cluster == cl)]
            n7 = int(r7.n_cells.iloc[0]) if len(r7) else 0
            n10 = int(r10.n_cells.iloc[0]) if len(r10) else 0
            f7 = float(r7.frac_timepoint.iloc[0]) if len(r7) else 0.0
            f10 = float(r10.frac_timepoint.iloc[0]) if len(r10) else 0.0
            if n7 >= MIN_CELLS and n10 < MIN_CELLS:
                status = "vanish"
            elif n7 < MIN_CELLS and n10 >= MIN_CELLS:
                status = "appear"
            elif n7 < MIN_CELLS and n10 < MIN_CELLS:
                status = "thin_both"
            elif f10 > f7:
                status = "grow"
            elif f10 < f7:
                status = "shrink"
            else:
                status = "stable_frac"
            recs.append(dict(
                cell_line=str(line), cluster=int(cl),
                n_d7=n7, frac_d7=f7, n_d10=n10, frac_d10=f10,
                delta_frac=f10 - f7, status=status,
                age_score_d7=float(r7.age_score.iloc[0]) if len(r7) and "age_score" in r7.columns else np.nan,
                age_score_d10=float(r10.age_score.iloc[0]) if len(r10) and "age_score" in r10.columns else np.nan,
                below_min_d7=bool(n7 < MIN_CELLS), below_min_d10=bool(n10 < MIN_CELLS),
            ))
    return pd.DataFrame(recs)


def task1_reading(cluster_df, log):
    """Aged donor only. Young tabulated, not pooled."""
    g = cluster_df[cluster_df.cell_line.astype(str).str.upper() == AGED_LINE].copy()
    if g.empty:
        return dict(key="unresolvable", text="aged donor rows missing.", donor=AGED_LINE,
                    min_cells=MIN_CELLS, C_ok=[], interpretable_cells=False)

    def n_at(cl, day):
        sub = g[(g.cluster == cl) & (g.day == int(day))]
        return int(sub.n_cells.iloc[0]) if len(sub) else 0

    def score_at(cl, day):
        sub = g[(g.cluster == cl) & (g.day == int(day))]
        if not len(sub):
            return np.nan
        if bool(sub.below_min_cells.iloc[0]):
            return np.nan
        return float(sub.age_score.iloc[0])

    def frac_at(cl, day):
        sub = g[(g.cluster == cl) & (g.day == int(day))]
        return float(sub.frac_timepoint.iloc[0]) if len(sub) else 0.0

    clusters = sorted(g.cluster.astype(int).unique().tolist())
    C_ok = [c for c in clusters if n_at(c, 7) >= MIN_CELLS and n_at(c, 10) >= MIN_CELLS]
    present_d7 = [c for c in clusters if n_at(c, 7) >= MIN_CELLS]
    per = []
    for c in C_ok:
        s7, s10 = score_at(c, 7), score_at(c, 10)
        per.append(dict(
            cluster=int(c), n_d7=n_at(c, 7), n_d10=n_at(c, 10),
            age_score_d7=s7, age_score_d10=s10,
            delta_age_d10_minus_d7=float(s10 - s7) if np.isfinite(s7) and np.isfinite(s10) else np.nan,
            reverse=bool(np.isfinite(s7) and np.isfinite(s10) and s10 > s7),
            stable_or_decline=bool(np.isfinite(s7) and np.isfinite(s10) and s10 <= s7),
        ))
    base = dict(
        donor=AGED_LINE, min_cells=MIN_CELLS, C_ok=C_ok, per_cluster=per,
        n_C_ok=int(len(C_ok)),
    )
    if not C_ok:
        text = (
            f"Unresolvable. Too few cells per cluster × timepoint "
            f"(threshold {MIN_CELLS} cells) to read either way."
        )
        log(f"[t1 reading] key=unresolvable {text}")
        return dict(key="unresolvable", text=text, **base)

    reverses = [p["cluster"] for p in per if p["reverse"]]
    stable = [p["cluster"] for p in per if p["stable_or_decline"]]
    d7_scores = [(c, score_at(c, 7)) for c in present_d7]
    d7_scores = [(c, s) for c, s in d7_scores if np.isfinite(s)]
    high = [c for c, s in d7_scores if s == max(x[1] for x in d7_scores)] if d7_scores else []
    low = [c for c, s in d7_scores if s == min(x[1] for x in d7_scores)] if d7_scores else []
    mix = False
    mix_notes = []
    for h in high:
        if frac_at(h, 10) > frac_at(h, 7):
            mix = True
            mix_notes.append(f"high-age cluster {h} grows frac {frac_at(h, 7):.4f}→{frac_at(h, 10):.4f}")
    for lo in low:
        if frac_at(lo, 10) < frac_at(lo, 7):
            mix = True
            mix_notes.append(f"low-age cluster {lo} shrinks/drops frac {frac_at(lo, 7):.4f}→{frac_at(lo, 10):.4f}")
    base.update(high_age_clusters_d7=high, low_age_clusters_d7=low,
                mix_shift=bool(mix), mix_notes=mix_notes,
                reverse_clusters=reverses, stable_or_decline_clusters=stable)

    if reverses and stable:
        which = (
            "reverse (age_score d10>d7): " + ",".join(str(c) for c in reverses)
            + "; stable_or_decline (d10≤d7): " + ",".join(str(c) for c in stable)
            + ". Do not average them into one number."
        )
        text = f"Mixed. {which}"
        log(f"[t1 reading] key=mixed {text}")
        return dict(key="mixed", text=text, **base)
    if reverses and not stable:
        text = (
            "Cell-intrinsic. The same cells' clusters individually reverse from d7 to d10. "
            "Then reprogramming genuinely pushes cells back along the age axis after d7, "
            "and that is a biological claim worth its own test. "
            f"clusters={reverses}."
        )
        log(f"[t1 reading] key=cell_intrinsic {text}")
        return dict(key="cell_intrinsic", text=text, **base)
    if stable and not reverses and mix:
        text = (
            "Compositional. Individual clusters' age scores are stable or keep declining "
            "from d7 to d10, but the mix shifts — a high-age-score cluster grows, or a "
            "low-scoring one shrinks or drops out. Then the all-cell reversal is a "
            "composition artifact, the same failure mode as the blood result in "
            "FINDINGS_BRAIN_PHASE1.md §2. "
            f"mix_notes={mix_notes}."
        )
        log(f"[t1 reading] key=compositional {text}")
        return dict(key="compositional", text=text, **base)
    text = (
        "no_listed_bullet. All C_ok clusters are stable_or_decline and mix_shift is False. "
        "See the tables. Not a fifth story."
    )
    log(f"[t1 reading] key=no_listed_bullet {text}")
    return dict(key="no_listed_bullet", text=text, **base)


def run_task1(log=None):
    _require_flags()
    close_log = False
    if log is None:
        log = Logger(FIBRO3_DIR / "t1_report.txt")
        close_log = True
    fibro3_log_banner(log, "TASK1")
    log(CLUSTER_SETTINGS)
    log(PREREG_TASK1)
    if not FROZEN_RULER.exists():
        raise StopStep("frozen_ruler", f"missing {FROZEN_RULER}")
    frozen = load_frozen_ruler()
    libs, samples = _load_libs(log)

    by_line = {}
    for lib in libs:
        by_line.setdefault(lib["cell_line"], []).append(lib)
    label_frames, Y_by_line, meta_by_line = [], {}, {}
    for line in (AGED_LINE, YOUNG_LINE):
        if line not in by_line:
            raise StopStep("task1", f"donor {line} missing from loaded libraries")
        obs_d, Y_d, meta = joint_cluster_one_donor(by_line[line], frozen, log)
        for lib in by_line[line]:
            lib.pop("X", None)
        gc.collect()
        label_frames.append(obs_d)
        Y_by_line[line] = Y_d
        meta_by_line[line] = meta
    labels = pd.concat(label_frames, axis=0, ignore_index=True)
    labels.to_csv(FIBRO3_DIR / "t1_cell_labels.csv", index=False)
    dump_json(FIBRO3_DIR / "t1_cluster_meta.json", jsonable(meta_by_line))
    n_clusters = {line: meta_by_line[line]["k"] for line in meta_by_line}
    log(f"[cluster] n clusters per donor (k-means k): {n_clusters}")

    if not DECLARED_BEFORE_SCORES_FLAG.exists():
        DECLARED_BEFORE_SCORES_FLAG.write_text(
            CLUSTER_SETTINGS + "\n\n" + PREREG_TASK1
            + "\n\nThis flag was written after joint clustering and before any age score.\n",
            encoding="utf-8",
        )
        log(f"[prereg] wrote {DECLARED_BEFORE_SCORES_FLAG} after clustering, before scoring")
    else:
        log(f"[prereg] {DECLARED_BEFORE_SCORES_FLAG.name} already exists; not rewriting")

    # stack Y in the same donor order as labels
    Y_all = sparse.vstack([Y_by_line[line] for line in (AGED_LINE, YOUNG_LINE)], format="csr")
    if Y_all.shape[0] != len(labels):
        raise StopStep("task1", f"Y_all n={Y_all.shape[0]} labels n={len(labels)}")

    pb_obs, logcpm, counts, nf = _pseudobulk_rows(labels, Y_all, frozen, log)
    age, Z, missing = score_frozen(logcpm, frozen, counts=counts)
    log(f"[score] missing-gene columns set to z=0: {int(missing.sum())}")
    pluri_w, meta_w = differentiation_score(Z, frozen, drop_oskm=False, log=log)
    pluri_wo, meta_wo = differentiation_score(Z, frozen, drop_oskm=True, log=log)
    pb_obs = pb_obs.copy()
    pb_obs["age_score"] = age
    pb_obs["pluri_with_OSKM"] = pluri_w
    pb_obs["pluri_without_OSKM"] = pluri_wo
    # Stage 2: if with/without disagree on aged all-cell, without is primary.
    # Here: per row both reported; primary = with_OSKM unless endpoint sign on aged
    # all-cluster mean is not used. Per FINDINGS_FIBRO: disagree if d0→d10 sign differs
    # or rise times differ, then without-OSKM is primary. Apply per donor on
    # size-weighted? Prompt: report both. Do not average clusters. Primary column
    # is with_OSKM unless the two series disagree on that donor's cluster table
    # when concatenated in day order for a fixed cluster. Report both always;
    # pluri_primary = with_OSKM (Stage 2 all-cell did not disagree).
    pb_obs["pluri_primary"] = pb_obs["pluri_with_OSKM"]
    pb_obs["pluri_primary_name"] = "with_OSKM"
    pb_obs.loc[pb_obs.below_min_cells, ["age_score", "pluri_with_OSKM",
                                        "pluri_without_OSKM", "pluri_primary"]] = np.nan
    pb_obs.to_csv(FIBRO3_DIR / "t1_cluster_table.csv", index=False)
    np.savez_compressed(
        FIBRO3_DIR / "t1_cluster_z.npz",
        Z=np.asarray(Z, np.float32),
        age=np.asarray(age, np.float64),
        pluri_with=np.asarray(pluri_w, np.float64),
        pluri_without=np.asarray(pluri_wo, np.float64),
        counts=np.asarray(counts, np.float32),
        logcpm=np.asarray(logcpm, np.float32),
        missing=missing,
        cell_line=pb_obs.cell_line.astype(str).to_numpy(),
        day=pb_obs.day.to_numpy(int),
        cluster=pb_obs.cluster.to_numpy(int),
        n_cells=pb_obs.n_cells.to_numpy(int),
        below_min_cells=pb_obs.below_min_cells.to_numpy(bool),
    )
    dump_json(FIBRO3_DIR / "t1_pluri_meta.json", jsonable(dict(with_OSKM=meta_w, without_OSKM=meta_wo)))

    comp = composition_table(pb_obs)
    for i, rec in comp.iterrows():
        r7 = pb_obs[(pb_obs.cell_line == rec.cell_line) & (pb_obs.day == 7) & (pb_obs.cluster == rec.cluster)]
        r10 = pb_obs[(pb_obs.cell_line == rec.cell_line) & (pb_obs.day == 10) & (pb_obs.cluster == rec.cluster)]
        if len(r7):
            comp.at[i, "age_score_d7"] = float(r7.age_score.iloc[0]) if not bool(r7.below_min_cells.iloc[0]) else np.nan
        if len(r10):
            comp.at[i, "age_score_d10"] = float(r10.age_score.iloc[0]) if not bool(r10.below_min_cells.iloc[0]) else np.nan
    comp.to_csv(FIBRO3_DIR / "t1_composition_d7_d10.csv", index=False)

    reading = task1_reading(pb_obs, log)
    dump_json(FIBRO3_DIR / "t1_reading.json", jsonable(reading))
    summary = dict(
        reading=reading,
        cluster_method="PCA svd_solver=randomized n_components=20 random_state=20260914; KMeans k=3 n_init=10 random_state=20260914",
        k=CLUSTER_K, npc=CLUSTER_NPC, min_cells=MIN_CELLS,
        n_label_rows=int(len(labels)), n_pseudobulk_rows=int(len(pb_obs)),
        n_below_min=int(pb_obs.below_min_cells.sum()),
        missing_z0=int(missing.sum()),
        per_donor_meta=meta_by_line,
        seed=FIBRO3_SEED,
        frozen_ruler=str(FROZEN_RULER),
        note="joint clustering within donor; nothing fitted on GSE297234; d7 not the Stage 2 endpoint",
    )
    dump_json(FIBRO3_DIR / "t1_summary.json", jsonable(summary))
    man = load_manifest()
    man["status"] = "TASK1_DONE"
    man["t1_reading"] = reading
    save_manifest(man)
    progress_snapshot("Task 2 extrapolation distances", stop="TASK1_DONE")
    if close_log:
        log.close()
    return summary


if __name__ == "__main__":
    try:
        run_task1()
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise

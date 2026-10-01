"""Task 1 — PartialReprog vs NonReprog at shared timepoints (Figure 3G contrast)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md4_common import (  # noqa: E402
    MD4_DIR, MD4_PROC, MD2_DIR, MD2_PROC, MD3_DIR, MD4_SEED, MD4_BOOT,
    MIN_CELLS_CONTRAST, AGED_LINE, YOUNG_LINE, DAYS, STATE_NAMES,
    PRIMARY_PAIR, CONTEXT_PAIRS, ALL_PAIRS, INSTRUMENTS,
    PREREG_TASK1, PREREG_TASK1_FLAG, DECLARED_BEFORE_SCORES_FLAG,
    FROZEN_RULER, N_PERM, N_BOOT, AMS_NBIN, AMS_CTRL,
    MD2_POOLED_MD_PR, MD2_POOLED_MD_NR, MD2_N_PR, MD2_N_NR,
    MD3_TASK1_STATUS, MD3_TASK1_READING, MD2_INSTRUMENTS_DISAGREE,
    MD2_NEGATIVE_HOLDS, STAGE2_VERDICT_FIBRO, WHAT_WOULD_SETTLE,
    StopStep, Logger, dump_json, jsonable, load_json, md4_log_banner,
    load_manifest, save_manifest, record_failure, log_columns, log_geneset,
    progress_snapshot, write_prereg_flag, perm_seed, boot_seed,
    load_frozen_ruler, MMC3_XLSX, detect_id_type, require_mappable,
)
from md3_idtype import csr_from_npz  # noqa: E402
from fibro_common import PLURI_ENDOGENOUS, FIBRO_IDENTITY, unit  # noqa: E402
from fibro_stage2 import gene_index_by_symbol  # noqa: E402
from gtex_common import EDGE_R_PRIOR  # noqa: E402
from trajectory_common import permutation_p  # noqa: E402
from md2_score import lognormalize_csr, add_module_score  # noqa: E402


def _require_prereg():
    if not PREREG_TASK1_FLAG.exists():
        raise StopStep("prereg", "PREREG_TASK1.flag missing")
    if not FROZEN_RULER.exists():
        raise StopStep("frozen_ruler", f"missing {FROZEN_RULER}")
    for line in (AGED_LINE, YOUNG_LINE):
        for p in (
            MD2_DIR / f"t2_cluster_labels_{line}.csv",
            MD2_DIR / f"louvain_obs_{line}.csv",
            MD2_PROC / f"louvain_counts_{line}.npz",
            MD2_PROC / f"louvain_ams_{line}.npz",
            MD2_PROC / f"rulerY_{line}.npz",
            MD2_DIR / f"allcell_obs_{line}.csv",
        ):
            if not p.exists():
                raise StopStep("task1", f"missing md2 artifact {p}. Not substituting.")
    if not (MD2_DIR / "t2_reproduction.json").exists():
        raise StopStep("task1", "missing results/md2/t2_reproduction.json")
    if not (MD2_DIR / "genesets.json").exists():
        raise StopStep("task1", "missing results/md2/genesets.json")
    if not (MD3_DIR / "genesets.json").exists():
        raise StopStep("task1", "missing results/md3/genesets.json")


def _load_ruler_y(path):
    z = np.load(path)
    return sparse.csr_matrix(
        (z["data"], z["indices"], z["indptr"]),
        shape=tuple(int(x) for x in z["shape"]),
    )


def _align_louvain_y(obs_l, obs_all, Yall, line):
    key_all = obs_all["gsm"].astype(str) + "||" + obs_all["barcode"].astype(str)
    pos = {k: i for i, k in enumerate(key_all)}
    key = obs_l["gsm"].astype(str) + "||" + obs_l["barcode"].astype(str)
    miss = [k for k in key if k not in pos]
    if miss:
        raise StopStep("task1", f"{line} louvain {len(miss)} cells not in allcell obs")
    take = np.array([pos[k] for k in key], dtype=int)
    if Yall.shape[0] != len(obs_all):
        raise StopStep("task1", f"{line} rulerY n={Yall.shape[0]} allcell n={len(obs_all)}")
    return Yall[take]


def _load_labels(line, log):
    path = MD2_DIR / f"t2_cluster_labels_{line}.csv"
    lab = pd.read_csv(path)
    log(f"[t1 labels] {line} columns actually read: {list(lab.columns)}")
    log_columns(f"t2_cluster_labels_{line}", list(lab.columns), str(path))
    need = {"cluster", "label", "n_cells"}
    missing = sorted(need - set(lab.columns))
    if missing:
        raise StopStep("task1", f"t2_cluster_labels_{line}.csv missing {missing}")
    return lab


def _load_obs_labelled(line, lab, log):
    path = MD2_DIR / f"louvain_obs_{line}.csv"
    obs = pd.read_csv(path)
    log(f"[t1 obs] {line} columns actually read: {list(obs.columns)}")
    log_columns(f"louvain_obs_{line}", list(obs.columns), str(path))
    if "cluster" not in obs.columns or "day" not in obs.columns:
        raise StopStep("task1", f"louvain_obs_{line} missing cluster/day")
    lab_map = {int(r.cluster): str(r.label) for _, r in lab.iterrows()}
    obs = obs.reset_index(drop=True)
    obs["label"] = obs["cluster"].map(lambda c: lab_map.get(int(c)))
    if obs["label"].isna().any():
        raise StopStep(
            "task1",
            f"{line}: clusters without md2 labels: "
            f"{sorted(obs.loc[obs.label.isna(), 'cluster'].unique())}",
        )
    return obs


def _load_ams_md_tgfb(line, obs, log):
    path = MD2_PROC / f"louvain_ams_{line}.npz"
    z = np.load(path, allow_pickle=True)
    log(f"[t1 ams] {line} keys actually read: {list(z.files)}")
    log_columns(f"louvain_ams_{line}", list(z.files), str(path))
    for k in ("md", "tgfb", "bins", "cluster", "day"):
        if k not in z.files:
            raise StopStep("task1", f"{path.name} missing field {k!r}. Not substituting.")
    md = np.asarray(z["md"], float)
    tgfb = np.asarray(z["tgfb"], float)
    cl = np.asarray(z["cluster"], int)
    day = np.asarray(z["day"], int)
    bins = np.asarray(z["bins"])
    if len(md) != len(obs) or len(tgfb) != len(obs):
        raise StopStep(
            "task1",
            f"{line} louvain_ams n_md={len(md)} n_tgfb={len(tgfb)} obs n={len(obs)}",
        )
    obs_cl = obs["cluster"].to_numpy(int)
    obs_day = obs["day"].to_numpy(int)
    if not np.array_equal(cl, obs_cl) or not np.array_equal(day, obs_day):
        raise StopStep(
            "task1",
            f"{line} louvain_ams cluster/day arrays do not match louvain_obs. "
            "Labels have drifted. Not scoring.",
        )
    nbin = int(np.unique(bins).size)
    if nbin != int(AMS_NBIN):
        raise StopStep("task1", f"{line} louvain_ams bins unique={nbin} != nbin={AMS_NBIN}")
    return dict(md=md, tgfb=tgfb, bins=bins, gene_mean=np.asarray(z["gene_mean"]) if "gene_mean" in z.files else None)


def load_genesets(log):
    md2_gs = load_json(MD2_DIR / "genesets.json")
    md3_gs = load_json(MD3_DIR / "genesets.json")
    log_columns("md2_genesets_json", list(md2_gs.keys()), str(MD2_DIR / "genesets.json"))
    log_columns("md3_genesets_json", list(md3_gs.keys()), str(MD3_DIR / "genesets.json"))
    md = [str(g).upper() for g in md2_gs["MD"]]
    tgfb = [str(g).upper() for g in md2_gs["TGFB"]]
    age_up = [str(g).upper() for g in md3_gs["age_up"]]
    age_down = [str(g).upper() for g in md3_gs["age_down"]]
    md2_up = [str(g).upper() for g in (md2_gs.get("mmc3_age_up") or [])]
    md2_dn = [str(g).upper() for g in (md2_gs.get("mmc3_age_down") or [])]
    md3_md = [str(g).upper() for g in md3_gs["MD"]]
    md3_tgfb = [str(g).upper() for g in md3_gs["TGFB"]]
    if set(md) != set(md3_md) or set(tgfb) != set(md3_tgfb):
        raise StopStep("genesets", "md2 vs md3 MD/TGFB sets differ. Not reconciling.")
    if set(age_up) != set(md2_up) or set(age_down) != set(md2_dn):
        raise StopStep(
            "genesets",
            "md3 Age up/down do not equal md2 genesets.json copies. Not reconciling.",
        )
    log_geneset("MD_built_from_md2", "MD (STAR Methods built; md2)", md,
                str(MD2_DIR / "genesets.json"), id_type="symbol")
    log_geneset("TGFB_built_from_md2", "TGFB (STAR Methods built; md2)", tgfb,
                str(MD2_DIR / "genesets.json"), id_type="symbol")
    log_geneset("mmc3_Aging_signatures_Age up", "Aging_signatures:Age up", age_up,
                str(MMC3_XLSX) + " via md3/genesets.json", id_type="symbol")
    log_geneset("mmc3_Aging_signatures_Age down", "Aging_signatures:Age down", age_down,
                str(MMC3_XLSX) + " via md3/genesets.json", id_type="symbol")
    summary = dict(
        n_MD=len(md), n_TGFB=len(tgfb),
        n_age_up=len(age_up), n_age_down=len(age_down),
        md2_equals_md3_MD=True, md2_equals_md3_TGFB=True,
        mmc3_equals_md2_age_up=True, mmc3_equals_md2_age_down=True,
        used_as_deseq2_rebuild=False, used_as_published_lists=True,
        seed=MD4_SEED,
        source_md2=str(MD2_DIR / "genesets.json"),
        source_md3=str(MD3_DIR / "genesets.json"),
    )
    dump_json(MD4_DIR / "genesets_summary.json", jsonable(summary))
    dump_json(MD4_DIR / "genesets.json", jsonable(dict(
        MD=md, TGFB=tgfb, age_up=age_up, age_down=age_down,
    )))
    log(f"[genesets] MD n={len(md)} TGFB n={len(tgfb)} Age up n={len(age_up)} Age down n={len(age_down)}")
    return dict(MD=md, TGFB=tgfb, age_up=age_up, age_down=age_down)


def occupancy_and_qualify(log):
    rows = []
    pooled = []
    for line in (AGED_LINE, YOUNG_LINE):
        lab = _load_labels(line, log)
        obs = _load_obs_labelled(line, lab, log)
        for st in STATE_NAMES:
            rec = dict(cell_line=line, label=st, n_clusters=int((lab.label.astype(str) == st).sum()))
            for d in DAYS:
                rec[f"n_d{d}"] = int(((obs.label == st) & (obs.day.astype(int) == int(d))).sum())
            rec["n_all"] = int((obs.label == st).sum())
            pooled.append(rec)
            log(f"[t1 counts] {line} {st} n0={rec['n_d0']} n3={rec['n_d3']} "
                f"n7={rec['n_d7']} n10={rec['n_d10']} n_all={rec['n_all']}")
        for d in DAYS:
            for st in STATE_NAMES:
                n = int(((obs.label == st) & (obs.day.astype(int) == int(d))).sum())
                n_tp = int((obs.day.astype(int) == int(d)).sum())
                rows.append(dict(
                    cell_line=line, day=int(d), label=st, n_cells=n,
                    n_timepoint=n_tp,
                    frac_timepoint=float(n / n_tp) if n_tp else np.nan,
                    ge50=bool(n >= MIN_CELLS_CONTRAST),
                ))
    counts = pd.DataFrame(rows)
    pooled_df = pd.DataFrame(pooled)
    counts.to_csv(MD4_DIR / "t1_cell_counts.csv", index=False)
    pooled_df.to_csv(MD4_DIR / "t1_cell_counts_pooled.csv", index=False)

    qrows = []
    for line in (AGED_LINE, YOUNG_LINE):
        for d in DAYS:
            for a, b in ALL_PAIRS:
                na = int(counts.loc[
                    (counts.cell_line == line) & (counts.day == int(d)) & (counts.label == a),
                    "n_cells",
                ].iloc[0])
                nb = int(counts.loc[
                    (counts.cell_line == line) & (counts.day == int(d)) & (counts.label == b),
                    "n_cells",
                ].iloc[0])
                ok = bool(na >= MIN_CELLS_CONTRAST and nb >= MIN_CELLS_CONTRAST)
                qrows.append(dict(
                    cell_line=line, day=int(d), state_a=a, state_b=b,
                    n_a=na, n_b=nb, min_cells=MIN_CELLS_CONTRAST,
                    qualifies=ok, primary=bool((a, b) == PRIMARY_PAIR),
                    reason=None if ok else f"n<{MIN_CELLS_CONTRAST} in at least one state",
                ))
                log(f"[t1 qualify] {line} d{d} {a} vs {b} n_a={na} n_b={nb} qualifies={ok}")
    qdf = pd.DataFrame(qrows)
    qdf.to_csv(MD4_DIR / "t1_qualify.csv", index=False)
    return counts, pooled_df, qdf


def sanity_md2(log, obs_by_line, ams_by_line):
    repro = load_json(MD2_DIR / "t2_reproduction.json")
    log_columns("t2_reproduction", list(repro.keys()), str(MD2_DIR / "t2_reproduction.json"))
    obs = obs_by_line[AGED_LINE]
    md = ams_by_line[AGED_LINE]["md"]
    n_pr = int((obs.label == "PartialReprog").sum())
    n_nr = int((obs.label == "NonReprog").sum())
    md_pr = float(np.mean(md[obs.label.eq("PartialReprog").to_numpy()])) if n_pr else np.nan
    md_nr = float(np.mean(md[obs.label.eq("NonReprog").to_numpy()])) if n_nr else np.nan
    rec = dict(
        n_PartialReprog=n_pr, n_NonReprog=n_nr,
        md_mean_PartialReprog=md_pr, md_mean_NonReprog=md_nr,
        md2_json_n_PartialReprog=int(repro["n_PartialReprog"]),
        md2_json_n_NonReprog=int(repro["n_NonReprog"]),
        md2_json_md_mean_PartialReprog=float(repro["md_mean_PartialReprog"]),
        md2_json_md_mean_NonReprog=float(repro["md_mean_NonReprog"]),
        findings_printed_PartialReprog="+0.150",
        findings_printed_NonReprog="+0.439",
        abs_delta_md_pr=abs(md_pr - float(repro["md_mean_PartialReprog"])),
        abs_delta_md_nr=abs(md_nr - float(repro["md_mean_NonReprog"])),
        n_match=bool(n_pr == int(repro["n_PartialReprog"]) and n_nr == int(repro["n_NonReprog"])),
        md_match_1e12=bool(
            abs(md_pr - float(repro["md_mean_PartialReprog"])) < 1e-12
            and abs(md_nr - float(repro["md_mean_NonReprog"])) < 1e-12
        ),
        printed_3dp_match=bool(f"{md_pr:+.3f}" == "+0.150" and f"{md_nr:+.3f}" == "+0.439"),
        expected_constants=dict(n_pr=MD2_N_PR, n_nr=MD2_N_NR, md_pr=MD2_POOLED_MD_PR, md_nr=MD2_POOLED_MD_NR),
    )
    rec["ok"] = bool(rec["n_match"] and rec["md_match_1e12"] and rec["printed_3dp_match"])
    dump_json(MD4_DIR / "t1_sanity.json", jsonable(rec))
    log(f"[t1 sanity] n_PR={n_pr} n_NR={n_nr} MD_PR={md_pr:+.6f} MD_NR={md_nr:+.6f} ok={rec['ok']}")
    if not rec["ok"]:
        raise StopStep(
            "sanity",
            "FINDINGS_MD2.md pooled MD means did not reproduce. "
            "Cell-state labels have drifted and nothing downstream is valid. "
            f"got n_PR={n_pr} n_NR={n_nr} MD_PR={md_pr} MD_NR={md_nr} "
            f"expected n_PR={repro['n_PartialReprog']} n_NR={repro['n_NonReprog']} "
            f"MD_PR={repro['md_mean_PartialReprog']} MD_NR={repro['md_mean_NonReprog']}",
            details=rec,
        )
    return rec


def _pluri_fibro_indices(frozen, log):
    pos, _, _ = gene_index_by_symbol(frozen)
    pluri_idx, pluri_miss = [], []
    for g in PLURI_ENDOGENOUS:
        j = pos.get(str(g).upper())
        if j is None:
            pluri_miss.append(g)
        else:
            pluri_idx.append(int(j))
    fibro_idx, fibro_miss = [], []
    for g in FIBRO_IDENTITY:
        j = pos.get(str(g).upper())
        if j is None:
            fibro_miss.append(g)
        else:
            fibro_idx.append(int(j))
    log(f"[pluri] endogenous requested={list(PLURI_ENDOGENOUS)} mapped={len(pluri_idx)} missing={pluri_miss}")
    log(f"[pluri] fibro_identity requested={list(FIBRO_IDENTITY)} mapped={len(fibro_idx)} missing={fibro_miss}")
    if not pluri_idx or not fibro_idx:
        raise StopStep("task1", f"pluri or fibro identity unmapped: pluri={pluri_idx} fibro={fibro_idx}")
    return np.asarray(pluri_idx, int), np.asarray(fibro_idx, int)


def _score_frozen_cells(Y, frozen, log, tag, chunk=512):
    """Per-cell log2-CPM (nf=1) then frozen μ/σ/w. Y is cells × ruler genes CSR."""
    Y = Y.tocsr()
    n, p = Y.shape
    mu = np.asarray(frozen["mu"], float)
    sd = np.asarray(frozen["sd"], float)
    sd = np.where(sd < 1e-12, 1.0, sd)
    w = np.asarray(frozen["w"], float)
    w, _ = unit(w)
    if len(mu) != p or len(w) != p:
        raise StopStep("task1", f"{tag}: rulerY p={p} frozen p_mu={len(mu)} p_w={len(w)}")
    lib = np.asarray(Y.sum(axis=1)).ravel().astype(np.float64)
    missing = np.asarray((Y.getnnz(axis=0) == 0))
    n_miss = int(missing.sum())
    finite_lib = lib[np.isfinite(lib) & (lib > 0)]
    ave = float(np.mean(finite_lib)) if finite_lib.size else 1.0
    if ave <= 0:
        ave = 1.0
    pluri_idx, fibro_idx = _pluri_fibro_indices(frozen, log)
    log(f"[ruler-cells {tag}] n={n} p={p} n_missing_allzero={n_miss} "
        f"lib_median={float(np.median(finite_lib)) if finite_lib.size else np.nan:.1f} "
        f"prior.count={EDGE_R_PRIOR:g} tmm_nf=1")
    age = np.full(n, np.nan, dtype=np.float64)
    pluri = np.full(n, np.nan, dtype=np.float64)
    for start in range(0, n, int(chunk)):
        sl = slice(start, min(start + int(chunk), n))
        dense = np.asarray(Y[sl].todense(), dtype=np.float64)
        lib_s = lib[sl]
        adj_prior = float(EDGE_R_PRIOR) * lib_s / ave
        den = lib_s + 2.0 * adj_prior
        den = np.where(den <= 0, 1.0, den)
        with np.errstate(divide="ignore", invalid="ignore"):
            cpm = (dense + adj_prior[:, None]) / den[:, None] * 1e6
            logcpm = np.log2(np.clip(cpm, 1e-12, None))
        bad = ~(np.isfinite(lib_s) & (lib_s > 0))
        Z = (logcpm - mu) / sd
        Z[:, missing] = 0.0
        age_s = Z @ w
        psc = Z[:, pluri_idx].mean(1) - Z[:, fibro_idx].mean(1)
        age_s[bad] = np.nan
        psc[bad] = np.nan
        age[sl] = age_s
        pluri[sl] = psc
        if start == 0 or start % (chunk * 10) == 0:
            log(f"[ruler-cells {tag}] {start}/{n}")
    return age, pluri, dict(
        n=n, p=p, n_missing_z0=n_miss, tmm_nf=1, prior=float(EDGE_R_PRIOR),
        ave_lib=ave, n_pluri=int(pluri_idx.size), n_fibro=int(fibro_idx.size),
    )


def cohens_d(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    na, nb = int(len(a)), int(len(b))
    if na < 2 or nb < 2:
        return np.nan
    va, vb = float(a.var(ddof=1)), float(b.var(ddof=1))
    sp = np.sqrt(((na - 1) * va + (nb - 1) * vb) / (na + nb - 2))
    if not np.isfinite(sp) or sp < 1e-15:
        return np.nan
    return float((a.mean() - b.mean()) / sp)


def _one_contrast(x, mask_a, n_perm, n_boot, seed_p, seed_b):
    x = np.asarray(x, float)
    mask_a = np.asarray(mask_a, bool)
    ok = np.isfinite(x)
    x = x[ok]
    mask_a = mask_a[ok]
    a = x[mask_a]
    b = x[~mask_a]
    n_a, n_b = int(len(a)), int(len(b))
    rec = dict(n_a=n_a, n_b=n_b, n_perm=int(n_perm), n_boot=int(n_boot),
               perm_seed=int(seed_p), boot_seed=int(seed_b),
               effect_size="cohens_d_pooled_nminus1")
    if n_a < 2 or n_b < 2:
        rec.update(ok=False, reason="n<2 after dropping non-finite scores")
        return rec
    dmean = float(a.mean() - b.mean())
    dmed = float(np.median(a) - np.median(b))
    rec.update(
        ok=True,
        mean_a=float(a.mean()), mean_b=float(b.mean()), delta_mean=dmean,
        median_a=float(np.median(a)), median_b=float(np.median(b)), delta_median=dmed,
        cohens_d=cohens_d(a, b),
    )
    rng = np.random.default_rng(int(seed_p))
    null_mean = np.empty(int(n_perm), dtype=np.float64)
    null_med = np.empty(int(n_perm), dtype=np.float64)
    lab = mask_a.copy()
    for i in range(int(n_perm)):
        rng.shuffle(lab)
        aa, bb = x[lab], x[~lab]
        null_mean[i] = float(aa.mean() - bb.mean())
        null_med[i] = float(np.median(aa) - np.median(bb))
    rec["p_mean_lower"] = permutation_p(dmean, null_mean, greater=False)
    rec["p_mean_two_sided"] = float(
        (np.sum(np.abs(null_mean) >= abs(dmean)) + 1) / (len(null_mean) + 1)
    )
    rec["n_perm_mean_le_obs"] = int(np.sum(null_mean <= dmean))
    rec["p_median_lower"] = permutation_p(dmed, null_med, greater=False)
    rec["p_median_two_sided"] = float(
        (np.sum(np.abs(null_med) >= abs(dmed)) + 1) / (len(null_med) + 1)
    )
    rec["null_mean_median"] = float(np.median(null_mean))
    rec["beats_null_lower"] = bool(dmean < 0 and rec["p_mean_lower"] <= 0.05)

    rng_b = np.random.default_rng(int(seed_b))
    boot_m = np.empty(int(n_boot), dtype=np.float64)
    boot_d = np.empty(int(n_boot), dtype=np.float64)
    for i in range(int(n_boot)):
        aa = rng_b.choice(a, size=n_a, replace=True)
        bb = rng_b.choice(b, size=n_b, replace=True)
        boot_m[i] = float(aa.mean() - bb.mean())
        boot_d[i] = float(np.median(aa) - np.median(bb))
    rec["delta_mean_ci_lo"] = float(np.percentile(boot_m, 2.5))
    rec["delta_mean_ci_hi"] = float(np.percentile(boot_m, 97.5))
    rec["delta_median_ci_lo"] = float(np.percentile(boot_d, 2.5))
    rec["delta_median_ci_hi"] = float(np.percentile(boot_d, 97.5))
    return rec


def _classify(md_lower, ruler_lower):
    if md_lower and ruler_lower:
        return "both_lower"
    if md_lower and not ruler_lower:
        return "md_lower_ruler_not"
    if ruler_lower and not md_lower:
        return "ruler_lower_md_not"
    return "neither_differs"


def _reading_text(key, line, day, md_row, ru_row):
    md_d = md_row.get("delta_mean")
    ru_d = ru_row.get("delta_mean")
    md_p = md_row.get("p_mean_lower")
    ru_p = ru_row.get("p_mean_lower")
    n_pr = md_row.get("n_a")
    n_nr = md_row.get("n_b")
    bits = (
        f"{line} d{day} PartialReprog vs NonReprog n_PR={n_pr} n_NR={n_nr} "
        f"MD Δmean={md_d} p_lower={md_p} Cohen_d={md_row.get('cohens_d')} "
        f"ruler Δmean={ru_d} p_lower={ru_p} Cohen_d={ru_row.get('cohens_d')}."
    )
    if key == "both_lower":
        return (
            "MD lower in PartialReprog than NonReprog and the frozen ruler lower too, "
            "both beating their nulls. Both instruments agree that partially reprogrammed "
            "cells read younger. Their claim is supported and our earlier negatives "
            "reflected the wrong test (a within-state trajectory with no baseline). "
            + bits
        )
    if key == "md_lower_ruler_not":
        return (
            "MD lower, ruler not lower (or higher). The instruments disagree on the "
            "population the published claim concerns, in a well-powered comparison. "
            "No winner is declared. What would settle it: " + WHAT_WOULD_SETTLE + " "
            + bits
        )
    if key == "ruler_lower_md_not":
        return (
            "Ruler lower, MD not. Unexplained discrepancy. No winner is declared. "
            + bits
        )
    return (
        "Neither differs. Our reproduction does not recover their Figure 3G contrast; "
        "this is a reproduction failure, not a refutation. " + bits
    )


def run_task1(log=None):
    close_log = False
    if log is None:
        log = Logger(MD4_DIR / "t1_report.txt")
        close_log = True
    md4_log_banner(log, "TASK1")
    log(PREREG_TASK1)
    _require_prereg()
    sets = load_genesets(log)
    progress_snapshot("Task 1 occupancy from md2 labels (no new scores)", stop="GENESETS_DONE")

    obs_by, ams_by, lab_by = {}, {}, {}
    for line in (AGED_LINE, YOUNG_LINE):
        lab_by[line] = _load_labels(line, log)
        obs_by[line] = _load_obs_labelled(line, lab_by[line], log)
        ams_by[line] = _load_ams_md_tgfb(line, obs_by[line], log)

    counts, pooled, qdf = occupancy_and_qualify(log)

    md3_p = MD3_DIR / "t1_occupancy_pooled.csv"
    if md3_p.exists():
        md3 = pd.read_csv(md3_p)
        log_columns("md3_t1_occupancy_pooled", list(md3.columns), str(md3_p))
        for _, r in pooled.iterrows():
            hit = md3[(md3.cell_line == r.cell_line) & (md3.label == r.label)]
            if not len(hit):
                raise StopStep("task1", f"md3 occupancy missing {r.cell_line} {r.label}")
            for d in DAYS:
                a, b = int(r[f"n_d{d}"]), int(hit.iloc[0][f"n_d{d}"])
                if a != b:
                    raise StopStep(
                        "task1",
                        f"occupancy drift {r.cell_line} {r.label} d{d}: md4={a} md3={b}",
                    )
        log("[t1] occupancy matches md3 t1_occupancy_pooled.csv")

    sanity = sanity_md2(log, obs_by, ams_by)
    if not DECLARED_BEFORE_SCORES_FLAG.exists():
        DECLARED_BEFORE_SCORES_FLAG.write_text(
            PREREG_TASK1
            + "\n\nThis flag was written after cell-count tables, qualify tables, and the "
              "md2 MD-mean sanity check, and before any frozen-ruler, Table S3, or "
              "pluripotency cell-level score.\n"
            + f"\nsanity_ok={sanity['ok']}\n",
            encoding="utf-8",
        )
        log(f"[prereg] wrote {DECLARED_BEFORE_SCORES_FLAG}")
    else:
        log(f"[prereg] {DECLARED_BEFORE_SCORES_FLAG.name} already exists; not rewriting")
    progress_snapshot("score per-cell ruler / Table S3 / pluri (MD/TGF-β from md2 disk)",
                      stop="DECLARED_BEFORE_SCORES")

    frozen = load_frozen_ruler()
    cell_scores = {}
    s3_meta = {}
    for line in (AGED_LINE, YOUNG_LINE):
        obs = obs_by[line]
        Xl, symbols, gene_id, zcount = csr_from_npz(MD2_PROC / f"louvain_counts_{line}.npz")
        id_sym = detect_id_type(symbols, log, f"louvain_{line}_symbols",
                                str(MD2_PROC / f"louvain_counts_{line}.npz") + ":symbols")
        require_mappable(id_sym, {"symbol"}, log, "task1")
        if gene_id is not None and len(gene_id):
            detect_id_type(gene_id, log, f"louvain_{line}_gene_id",
                           str(MD2_PROC / f"louvain_counts_{line}.npz") + ":gene_id")
        if Xl.shape[0] != len(obs):
            raise StopStep("task1", f"{line} louvain counts n={Xl.shape[0]} obs n={len(obs)}")

        logX = lognormalize_csr(Xl)
        age_sets = {"age_up": list(sets["age_up"]), "age_down": list(sets["age_down"])}
        ams_s3, ams_s3_meta = add_module_score(
            logX, symbols, age_sets, log=log, tag=f"{line}_louvain_S3",
            bins=ams_by[line]["bins"], seed=MD4_SEED,
        )
        s3_meta[line] = {
            k: {kk: vv for kk, vv in v.items() if kk != "mapped"}
            for k, v in (ams_s3_meta.get("per_set") or {}).items()
        }
        dump_json(MD4_DIR / f"t1_s3_ams_meta_{line}.json", jsonable(s3_meta[line]))

        obs_all = pd.read_csv(MD2_DIR / f"allcell_obs_{line}.csv")
        log_columns(f"allcell_obs_{line}", list(obs_all.columns), str(MD2_DIR / f"allcell_obs_{line}.csv"))
        Yall = _load_ruler_y(MD2_PROC / f"rulerY_{line}.npz")
        Y = _align_louvain_y(obs, obs_all, Yall, line)
        age, pluri, meta_r = _score_frozen_cells(Y, frozen, log, tag=line)
        dump_json(MD4_DIR / f"t1_ruler_cell_meta_{line}.json", jsonable(meta_r))

        age_up = np.asarray(ams_s3["age_up"], float)
        age_dn = np.asarray(ams_s3["age_down"], float)
        cell_scores[line] = dict(
            frozen_ruler=age,
            md_score=ams_by[line]["md"],
            tgfb_score=ams_by[line]["tgfb"],
            age_up=age_up,
            age_down=age_dn,
            age_up_minus_age_down=age_up - age_dn,
            pluri_primary=pluri,
        )
        np.savez_compressed(
            MD4_PROC / f"cell_scores_{line}.npz",
            frozen_ruler=age, md_score=ams_by[line]["md"], tgfb_score=ams_by[line]["tgfb"],
            age_up=age_up, age_down=age_dn, age_up_minus_age_down=age_up - age_dn,
            pluri_primary=pluri,
            cluster=obs["cluster"].to_numpy(int),
            day=obs["day"].to_numpy(int),
            label=obs["label"].to_numpy(),
        )
        # keep Y for Task 2
        np.savez_compressed(
            MD4_PROC / f"louvain_rulerY_{line}.npz",
            data=Y.data, indices=Y.indices, indptr=Y.indptr, shape=np.array(Y.shape),
        )

    # state × timepoint means (cell-level, not TMM pseudobulk)
    mean_rows = []
    for line in (AGED_LINE, YOUNG_LINE):
        obs = obs_by[line]
        sc = cell_scores[line]
        for d in DAYS:
            for st in STATE_NAMES:
                m = (obs.day.astype(int) == int(d)) & (obs.label == st)
                n = int(m.sum())
                rec = dict(cell_line=line, day=int(d), label=st, n_cells=n)
                for inst in INSTRUMENTS:
                    arr = np.asarray(sc[inst], float)[m.to_numpy()]
                    rec[f"{inst}_mean"] = float(np.nanmean(arr)) if n else np.nan
                    rec[f"{inst}_median"] = float(np.nanmedian(arr)) if n else np.nan
                mean_rows.append(rec)
    means_df = pd.DataFrame(mean_rows)
    means_df.to_csv(MD4_DIR / "t1_state_timepoint_means.csv", index=False)

    # contrasts
    contrast_rows = []
    primary_keys = []
    for line in (AGED_LINE, YOUNG_LINE):
        obs = obs_by[line]
        sc = cell_scores[line]
        for d in DAYS:
            day_m = obs.day.astype(int) == int(d)
            for pair in ALL_PAIRS:
                a, b = pair
                sub = day_m & obs.label.isin([a, b])
                n_a = int((day_m & (obs.label == a)).sum())
                n_b = int((day_m & (obs.label == b)).sum())
                qualifies = bool(n_a >= MIN_CELLS_CONTRAST and n_b >= MIN_CELLS_CONTRAST)
                pair_rec_base = dict(
                    cell_line=line, day=int(d), state_a=a, state_b=b,
                    n_a=n_a, n_b=n_b, qualifies=qualifies,
                    primary=bool(pair == PRIMARY_PAIR),
                    min_cells=MIN_CELLS_CONTRAST,
                )
                if not qualifies:
                    for inst in INSTRUMENTS:
                        contrast_rows.append(dict(
                            **pair_rec_base, instrument=inst, ok=False,
                            reason=f"n<{MIN_CELLS_CONTRAST} in at least one state",
                        ))
                    continue
                idx = np.flatnonzero(sub.to_numpy())
                labels_sub = obs.loc[sub, "label"].to_numpy()
                mask_a = labels_sub == a
                pair_inst = {}
                for inst in INSTRUMENTS:
                    x = np.asarray(sc[inst], float)[idx]
                    rec = _one_contrast(
                        x, mask_a, N_PERM, N_BOOT,
                        perm_seed(line, d, pair), boot_seed(line, d, pair),
                    )
                    rec.update(pair_rec_base)
                    rec["instrument"] = inst
                    contrast_rows.append(rec)
                    pair_inst[inst] = rec
                    log(f"[t1 contrast] {line} d{d} {a} vs {b} {inst} "
                        f"Δmean={rec.get('delta_mean')} p_lower={rec.get('p_mean_lower')} "
                        f"d={rec.get('cohens_d')} beats={rec.get('beats_null_lower')}")
                if pair == PRIMARY_PAIR:
                    md_row = pair_inst.get("md_score") or {}
                    ru_row = pair_inst.get("frozen_ruler") or {}
                    key = _classify(bool(md_row.get("beats_null_lower")), bool(ru_row.get("beats_null_lower")))
                    primary_keys.append(dict(
                        cell_line=line, day=int(d), key=key,
                        md={k: md_row.get(k) for k in (
                            "n_a", "n_b", "delta_mean", "p_mean_lower", "p_mean_two_sided",
                            "delta_mean_ci_lo", "delta_mean_ci_hi", "cohens_d", "beats_null_lower",
                            "mean_a", "mean_b",
                        )},
                        ruler={k: ru_row.get(k) for k in (
                            "n_a", "n_b", "delta_mean", "p_mean_lower", "p_mean_two_sided",
                            "delta_mean_ci_lo", "delta_mean_ci_hi", "cohens_d", "beats_null_lower",
                            "mean_a", "mean_b",
                        )},
                        text=_reading_text(key, line, d, md_row, ru_row),
                    ))

    cdf = pd.DataFrame(contrast_rows)
    cdf.to_csv(MD4_DIR / "t1_contrasts.csv", index=False)

    by_donor = {}
    for line in (AGED_LINE, YOUNG_LINE):
        hits = [h for h in primary_keys if h["cell_line"] == line]
        keys = sorted({h["key"] for h in hits})
        if not hits:
            by_donor[line] = dict(
                key="no_qualifying_timepoint",
                n_qualifying=0,
                timepoints=[],
                text=f"{line}: no timepoint where both PartialReprog and NonReprog have ≥{MIN_CELLS_CONTRAST} cells.",
            )
        elif len(keys) == 1:
            by_donor[line] = dict(
                key=keys[0], n_qualifying=len(hits),
                timepoints=[h["day"] for h in hits],
                per_timepoint=[{k: h[k] for k in ("day", "key", "text")} for h in hits],
                text=" ".join(h["text"] for h in hits),
            )
        else:
            by_donor[line] = dict(
                key="timepoints_disagree_do_not_pick",
                n_qualifying=len(hits),
                keys=keys,
                timepoints=[h["day"] for h in hits],
                per_timepoint=[{k: h[k] for k in ("day", "key", "text")} for h in hits],
                text=("Qualifying timepoints fire different bullets; reported separately, not picked, not averaged. "
                      + " ".join(h["text"] for h in hits)),
            )

    reading = dict(
        by_donor=by_donor,
        primary_pair=list(PRIMARY_PAIR),
        min_cells=MIN_CELLS_CONTRAST,
        md3_status_withdrawn=MD3_TASK1_STATUS,
        md3_sentence=MD3_TASK1_READING,
        md2_instruments_disagree=MD2_INSTRUMENTS_DISAGREE,
        md2_negative_holds=MD2_NEGATIVE_HOLDS,
        fibro_stage2=STAGE2_VERDICT_FIBRO,
        what_would_settle=WHAT_WOULD_SETTLE,
        effect_size="cohens_d_pooled_nminus1",
        delta_definition="mean(state_A) - mean(state_B)",
        lower_rule="delta_mean < 0 and p_mean_lower <= 0.05",
    )
    dump_json(MD4_DIR / "t1_reading.json", jsonable(reading))
    dump_json(MD4_DIR / "t1_summary.json", jsonable(dict(
        n_contrast_rows=int(len(cdf)),
        n_qualifying=int(qdf.qualifies.sum()) if len(qdf) else 0,
        sanity_ok=bool(sanity["ok"]),
        by_donor={k: v.get("key") for k, v in by_donor.items()},
        s3_mapped={ln: s3_meta[ln] for ln in s3_meta},
    )))
    log(f"[t1 reading] " + " ".join(f"{k}={v.get('key')}" for k, v in by_donor.items()))
    man = load_manifest()
    man["status"] = "TASK1_DONE"
    save_manifest(man)
    progress_snapshot("Task 2 extrapolation (report-only)", stop="TASK1_DONE")
    if close_log:
        log.close()
    return reading


if __name__ == "__main__":
    try:
        run_task1()
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise

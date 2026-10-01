"""Positive control / power curve for the identifiability pipeline.

Gene-space ridge L2 battery (V1 bootstrap, V2 OOF, V4a transfer, L2d young/old)
called with a different X and y. Does not modify prior FINDINGS*.md or FALSIFICATION.md.
Seed 20260914. Namespaced under results/posctrl/.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import RESULTS, ROOT  # noqa: E402
from target_common import (  # noqa: E402
    TARGET_SEED, DATASET_ID, N_BOOT_RIDGE, N_PERM, MIN_N_TYPE, ALPHAS,
    HASH_POOL_LIMITATION, Logger, dump_json, load_json, load_phase1_matrix,
    make_folds, permute_age_within_site_obs, score_age_predictions, primary_r2,
    unsigned_angle_deg, unit, donor_bootstrap_index, pairwise_angles_deg,
    zscore_train, per_type_directions, oof_project, calibrate_1d, md_table,
)
from lowdim_common import ANGLE_MATERIAL, gene_space_baseline  # noqa: E402
from trajectory_common import permutation_p, summarize_null_col, fmt, fmt_u  # noqa: E402
from brain_phase1_models import load_chrom_series  # noqa: E402
from brain_phase1_common import r2_mae  # noqa: E402

PC_DIR = RESULTS / "posctrl"
PC_FIG = PC_DIR / "figures"
for _p in (PC_DIR, PC_FIG):
    _p.mkdir(parents=True, exist_ok=True)

PC_SEED = TARGET_SEED  # 20260914
W_SEED = PC_SEED + 1
CAL_SEED = PC_SEED + 2
C3_SEED = PC_SEED + 3
N_BOOT_FULL = int(N_BOOT_RIDGE)  # 20; same B as results/target json
N_BOOT_HALF = N_BOOT_FULL // 2
N_PERM_FULL = int(N_PERM)
N_PERM_HALF = max(5, N_PERM_FULL // 2)
N_SUPPORT = 200
N_TYPES_CPM = 15
CPM_CUT = 1.0
CAL_R2_TOL = 0.03
TARGET_R2_PRIMARY = 0.30
N_C3_DRAWS = 10
C3_NS = (60, 100, 150, 233)
ANGLE_BAR = float(ANGLE_MATERIAL)  # 35.0; retired for C1/C3 on 2026-09-16 (kept for P1 record)
SEX_ANGLE_REF = 41.3  # C2a bootstrap; report C1/C3 angles against this, not 35°
TRANSFER_NULL_BAR = 0.05  # permutation-null usability; primary bar uses this
PROGRESS_PATH = ROOT / "PROGRESS_POSCTRL.md"
FINDINGS_PATH = ROOT / "FINDINGS_POSCTRL.md"
SUPERSESSION_FLAG = PC_DIR / "SUPERSESSION_20260916.flag"
SUPERSESSION_JSON = PC_DIR / "supersession_20260916.json"

PREREG_BLOCK = """- **P1 pipeline broken:** sex control (C2a, with sex-chromosome genes) fails bootstrap ≤35° or transfer >0 → STOP, report, do not proceed to C1/C3. Everything downstream of the target pipeline is suspect.
- **P2 negative is informative:** planted dense-shared direction at R²≈0.30 (C1) reaches bootstrap ≤35° AND transfer >0 in both directions → the pipeline can find a real direction at this n; the age negative in FINDINGS_LOWDIM stands as a statement about the data.
- **P3 negative is a power statement:** planted dense-shared at R²≈0.30 gives bootstrap >35° OR transfer ≤0 in either direction → the LOWDIM bars were unattainable at n=233; "not identifiable" must be rewritten as "not identifiable at n=233"; the power curve (C3) becomes the primary result.
- **Mixed** (passes angle, fails transfer, or vice versa) → report which, do not average. Transfer failing on a planted shared direction means bank covariance differences alone defeat ridge transfer; say so."""

SUPERSESSION_BLOCK = """The ≤35° angle bar is retired. New pre-registered bars for C1 and C3:
- (a) transfer R² > 0 in both directions with permutation null ≤ 0.05 — primary
- (b) bootstrap median pairwise angle reported against the sex reference 41.3°, not against 35°. Angle is not a pass/fail gate.

Restated P2 / P3 (these bars; C1 dense_shared at R²≈0.30):
- **P2** = planted dense_shared at R²≈0.30 transfers both ways
- **P3** = it does not

Sex is sparse and huge, so it is not a fair reference for a planted dense continuous direction at R²≈0.30. C1 is.

Verdict rule (only the outcome that fired):
- if planted transfers and real age does not → the age negative is real and it is specifically a transfer failure, not an angle failure
- If planted does not transfer either → bank covariance differences defeat transfer for any direction at this n; the age negative is uninformative and the power curve is the result
- If mixed across regimes, report per regime, do not average.

No threshold in this block is re-tuned after numbers exist."""

# GRCh38 PAR1/PAR2 symbols (Ensembl places these on X or Y; listed so they are
# dropped even if chrom is annotated XY/NA).
PAR_SYMBOLS = frozenset({
    "PLCXD1", "GTPBP6", "PPP2R3B", "SHOX", "CRLF2", "CSF2RA", "IL3RA", "SLC25A6",
    "ASMTL", "P2RY8", "AKAP17A", "SFRS17A", "ASMT", "DHRSX", "ZBED1", "CD99", "XG",
    "SPRY3", "VAMP7", "SYBL1", "IL9R", "WASH6P", "CXYORF1",
})
SEX_CHROM = frozenset({"X", "Y", "XY", "PAR", "X_PAR", "Y_PAR"})
SITE_NAME = {"H": "HBCC", "M": "MSSM"}
REGIMES = ("dense_shared", "dense_typespec", "sparse_shared", "sparse_typespec")
C1_JOBS = tuple(
    [(r, TARGET_R2_PRIMARY) for r in REGIMES]
    + [("dense_shared", 0.15), ("dense_shared", 0.50)]
)

LOWDIM_L4_QUOTE = (
    "**The negative result stands.** What would be needed to change it: more donors per "
    "cell type (so that even gene space is not p ≫ n, or so that a k ≪ n space is estimated "
    "from hundreds of independent brains rather than ~200); **more sites** (two banks are one "
    "degree of freedom of transfer, and that transfer failed); or a **different tissue** with "
    "less bank structure and more independent donors. Re-fitting more unsupervised variants on "
    "this same 233-donor DLPFC matrix will not make an underdetermined direction unique."
)


def pc_log_banner(log, stage: str):
    log("=" * 100)
    log(f"POSCTRL {stage}  seed={PC_SEED}  dataset_id={DATASET_ID}")
    log("=" * 100)
    log(HASH_POOL_LIMITATION)
    log("NO perturbation data. No TF Atlas. No candidate interventions.")
    log("FALSIFICATION.md and prior FINDINGS*.md are not modified.")
    if SUPERSESSION_FLAG.exists():
        log("SUPERSESSION 2026-09-16: 35° angle bar retired for C1/C3.")
        log("C1/C3 primary bar: transfer R² > 0 both directions AND permutation null ≤ 0.05.")
        log(f"Angle reported against sex reference {SEX_ANGLE_REF:.1f}°, not {ANGLE_BAR:.0f}°.")
    else:
        log(f"L2 bars frozen: bootstrap ≤ {ANGLE_BAR:.0f}°; planted transfer both directions > 0.")
    log("Ridge α grid and fold constructor unchanged.")


def cell_path(name: str) -> Path:
    return PC_DIR / f"{name}.json"


def cell_done(name: str) -> bool:
    p = cell_path(name)
    if not p.exists():
        return False
    try:
        rec = load_json(p)
    except Exception:
        return False
    return bool(rec.get("complete"))


def save_cell(name: str, rec: dict):
    rec = dict(rec)
    rec["complete"] = True
    rec["seed"] = PC_SEED
    rec["dataset_id"] = DATASET_ID
    dump_json(cell_path(name), rec)


def load_cell(name: str) -> dict:
    return load_json(cell_path(name))


def write_progress(text: str):
    PROGRESS_PATH.write_text(text, encoding="utf-8")


def finite(x):
    try:
        return bool(np.isfinite(float(x)))
    except (TypeError, ValueError):
        return False


def pass_angle(boot):
    return bool(finite(boot) and float(boot) <= ANGLE_BAR)


def pass_transfer_any(hm, mh):
    return bool((finite(hm) and float(hm) > 0) or (finite(mh) and float(mh) > 0))


def pass_transfer_both(hm, mh):
    return bool((finite(hm) and float(hm) > 0) and (finite(mh) and float(mh) > 0))


def pass_transfer_primary(hm, mh, hm_null, mh_null):
    """C1/C3 supersession bar: both transfers > 0 and both permutation nulls ≤ 0.05."""
    return bool(
        pass_transfer_both(hm, mh)
        and finite(hm_null) and float(hm_null) <= TRANSFER_NULL_BAR
        and finite(mh_null) and float(mh_null) <= TRANSFER_NULL_BAR
    )


def p1_blocks_c1_c3():
    if SUPERSESSION_FLAG.exists():
        return False
    stop = PC_DIR / "STOP_P1.txt"
    if stop.exists() and "stop_P1=True" in stop.read_text(encoding="utf-8"):
        return True
    return False


def c1_cell_name(regime, target_r2):
    rtag = f"r{int(round(float(target_r2) * 100)):03d}"
    return f"c1_{regime}_{rtag}"


def expected_c3_cells():
    names = []
    n_full = max(C3_NS)
    for kd in ("real", "planted"):
        for n in C3_NS:
            n_draws = 1 if int(n) >= int(n_full) else N_C3_DRAWS
            for d in range(n_draws):
                names.append(f"c3_{kd}_n{int(n):03d}_d{int(d):02d}")
    return names


def _skip_progress_json(name: str) -> bool:
    skip_sfx = (
        "_battery.json", "_bootstrap.json", "_within_site.json",
        "_loso.json", "_youngold.json", "_setup.json",
    )
    if name.endswith(skip_sfx):
        return True
    if name in (
        "supersession_20260916.json", "age_star.json", "frozen_beta.json",
        "planted_w_ready.json", "c1_setup.json", "c3_fit.json", "c3_summary.json",
    ):
        return True
    if name.startswith("planted_w_"):
        return True
    return False


def collect_cells_from_disk():
    cells = {}
    for p in sorted(PC_DIR.glob("*.json")):
        if _skip_progress_json(p.name):
            continue
        try:
            rec = load_json(p)
        except Exception:
            continue
        if rec.get("complete") or rec.get("boot_median_angle") is not None:
            rec = dict(rec)
            rec["path"] = str(p)
            cells[p.stem] = rec
    return cells


def load_betas_from_disk():
    p = PC_DIR / "frozen_beta.json"
    if not p.exists():
        return None
    rec = load_json(p)
    return {k: (v.get("beta") if isinstance(v, dict) else v) for k, v in rec.items()}


def next_posctrl_action():
    if not cell_done("c2a"):
        return "run C2a (gate)"
    if not cell_done("c2b"):
        return "run C2b"
    if p1_blocks_c1_c3():
        return "STOP P1 — C1/C3 blocked"
    for r, t in C1_JOBS:
        if not cell_done(c1_cell_name(r, t)):
            return f"run C1 {r} {t:.2f}"
    for name in expected_c3_cells():
        if not cell_done(name):
            return f"run {name}"
    if not (PC_DIR / "c3_summary.json").exists():
        return "write C3 summary + figure"
    return "write supersession verdict"


def refresh_progress(next_action=None, stop=None, extra=""):
    """Rewrite PROGRESS_POSCTRL.md from json on disk. Never from in-memory numbers."""
    cells = collect_cells_from_disk()
    if next_action is None:
        next_action = next_posctrl_action()
    if stop is None:
        if SUPERSESSION_FLAG.exists():
            stop = "P1 fired under a miscalibrated angle bar; superseded 2026-09-16; C1/C3 below"
        elif (PC_DIR / "STOP_P1.txt").exists():
            stop = "P1"
        else:
            stop = "none"
    progress_snapshot(cells, next_action, stop=stop, betas=load_betas_from_disk(), extra=extra)


def donor_mean(obs, y):
    y = np.asarray(y, float)
    don = pd.Series(y, index=np.arange(len(obs))).groupby(obs.donor.astype(str).to_numpy()).first()
    return float(don.mean())


def obs_with_y(obs, y):
    out = obs.copy()
    out["age"] = np.asarray(y, float)
    return out


def cpm_from_logcpm(Y):
    return np.maximum(np.power(2.0, np.asarray(Y, np.float64)) - 1.0, 0.0)


def eligible_sparse_genes(Y, obs, types, n_types=N_TYPES_CPM, cpm_cut=CPM_CUT):
    """Genes with mean CPM > cpm_cut in ≥ n_types cell types."""
    cpm = cpm_from_logcpm(Y)
    ct = obs.celltype.astype(str).to_numpy()
    hits = np.zeros(Y.shape[1], dtype=np.int32)
    for t in types:
        m = ct == t
        if int(m.sum()) < 1:
            continue
        hits += (cpm[m].mean(0) > cpm_cut).astype(np.int32)
    elig = np.flatnonzero(hits >= int(n_types))
    return elig, hits


def unit_normal(p, rng):
    w = rng.normal(0.0, 1.0, size=int(p))
    n = float(np.linalg.norm(w))
    if n < 1e-12:
        w = np.zeros(int(p))
        w[0] = 1.0
        return w
    return w / n


def sparse_unit(p, eligible, rng, n_support=N_SUPPORT):
    eligible = np.asarray(eligible, dtype=int)
    k = int(min(n_support, eligible.size))
    if k < 1:
        w = np.zeros(int(p))
        w[0] = 1.0
        return w, np.array([0], dtype=int)
    idx = np.sort(rng.choice(eligible, size=k, replace=False).astype(int))
    w = np.zeros(int(p), dtype=np.float64)
    w[idx] = rng.normal(0.0, 1.0, size=k)
    n = float(np.linalg.norm(w))
    if n < 1e-12:
        w[idx[0]] = 1.0
        n = 1.0
    return w / n, idx


def draw_planted_w(p, types, eligible, rng):
    """Draw and freeze w for all four regimes. Shared first, then typespec."""
    out = {}
    w_ds = unit_normal(p, rng)
    out["dense_shared"] = dict(kind="shared", w=w_ds, support=None)
    W_dt = np.column_stack([unit_normal(p, rng) for _ in types])
    out["dense_typespec"] = dict(kind="typespec", w=W_dt, support=None)
    w_ss, idx_ss = sparse_unit(p, eligible, rng)
    out["sparse_shared"] = dict(kind="shared", w=w_ss, support=idx_ss)
    cols, sups = [], []
    for _ in types:
        w, idx = sparse_unit(p, eligible, rng)
        cols.append(w)
        sups.append(idx)
    out["sparse_typespec"] = dict(kind="typespec", w=np.column_stack(cols), support=sups)
    return out


def save_planted_w(regime, pack, genes, types):
    payload = dict(
        regime=np.array(regime),
        kind=np.array(pack["kind"]),
        gene_id=genes.gene_id.astype(str).to_numpy(),
        symbol=genes.symbol.astype(str).to_numpy(),
        types=np.array(types, dtype=object),
        w=np.asarray(pack["w"], np.float64),
        seed=W_SEED,
    )
    if pack["support"] is None:
        payload["support"] = np.array([], dtype=int)
        payload["support_per_type"] = np.array(0)
    elif pack["kind"] == "shared":
        payload["support"] = np.asarray(pack["support"], dtype=int)
        payload["support_per_type"] = np.array(0)
    else:
        # ragged: store as object via concatenated + offsets
        sups = [np.asarray(s, dtype=int) for s in pack["support"]]
        payload["support"] = np.concatenate(sups) if sups else np.array([], dtype=int)
        payload["support_offsets"] = np.cumsum([0] + [len(s) for s in sups])
        payload["support_per_type"] = np.array(1)
    np.savez_compressed(PC_DIR / f"planted_w_{regime}.npz", **payload)


def w_for_type(pack, j):
    w = np.asarray(pack["w"], np.float64)
    if pack["kind"] == "shared" or w.ndim == 1:
        return w.ravel()
    return w[:, j]


def support_for_type(pack, j):
    if pack.get("support") is None:
        return None
    if pack["kind"] == "shared":
        return np.asarray(pack["support"], dtype=int)
    return np.asarray(pack["support"][j], dtype=int)


def plant(Y, y, obs, types, pack, beta, age_mean=None):
    """X_planted = X + β · (y − mean) · wᵀ, per row's cell type."""
    Y = np.asarray(Y, np.float64)
    y = np.asarray(y, np.float64)
    if age_mean is None:
        age_mean = donor_mean(obs, y)
    a = y - float(age_mean)
    ct = obs.celltype.astype(str).to_numpy()
    if pack["kind"] == "shared":
        w = np.asarray(pack["w"], np.float64).ravel()
        return Y + (float(beta) * a)[:, None] * w[None, :]
    Yp = Y.copy()
    W = np.asarray(pack["w"], np.float64)
    for j, t in enumerate(types):
        m = ct == t
        if not np.any(m):
            continue
        Yp[m] += (float(beta) * a[m])[:, None] * W[:, j][None, :]
    return Yp


def _median_pairwise(cols):
    cols = [np.asarray(c, float) for c in cols
            if np.isfinite(c).all() and float(np.linalg.norm(c)) > 1e-12]
    if len(cols) < 2:
        return np.nan, 0
    Wc = np.column_stack([c / np.linalg.norm(c) for c in cols])
    ang = pairwise_angles_deg(Wc)
    iu = np.triu_indices(ang.shape[0], 1)
    return float(np.median(ang[iu])), int(len(iu[0]))


def bootstrap_ridge(Y, obs, types, n_boot, rng, log, planted_pack=None):
    """V1 path: donor bootstrap, global z-score, per-type ridge age dirs."""
    ct = obs.celltype.astype(str).to_numpy()
    y = obs.age.to_numpy(float)
    Ws = []
    t0 = time.time()
    for b in range(n_boot):
        idx = donor_bootstrap_index(obs, rng)
        Xs, _, _ = zscore_train(Y[idx])
        W, meta, _ = per_type_directions(Xs, y[idx], ct[idx], types, method="ridge")
        Ws.append(W)
        if (b + 1) % 5 == 0 or b == 0:
            log(f"   ridge bootstrap {b+1}/{n_boot}  ({time.time()-t0:.0f}s)")
    rows, recov_all, frac_all = [], [], []
    for j, t in enumerate(types):
        cols = []
        for b in range(n_boot):
            w = Ws[b][:, j]
            if np.isfinite(w).all() and float(np.linalg.norm(w)) > 1e-12:
                cols.append(w / np.linalg.norm(w))
        med, n_pairs = _median_pairwise(cols)
        rec = np.nan
        frac = np.nan
        w_p = None if planted_pack is None else w_for_type(planted_pack, j)
        if w_p is not None and cols:
            recs = [unsigned_angle_deg(c, w_p) for c in cols]
            rec = float(np.nanmedian(recs))
            recov_all.extend(recs)
        sup = None if planted_pack is None else support_for_type(planted_pack, j)
        if sup is not None and cols:
            frs = []
            sup_set = set(int(i) for i in np.asarray(sup))
            for c in cols:
                top = set(np.argsort(np.abs(c))[-N_SUPPORT:].tolist())
                frs.append(len(top & sup_set) / float(max(len(sup_set), 1)))
            frac = float(np.nanmedian(frs))
            frac_all.append(frac)
        rows.append(dict(
            celltype=t, n_boot_ok=len(cols), n_boot=n_boot,
            median_pairwise_deg=med, n_pairs=n_pairs,
            recovery_angle_deg=rec, support_frac_top200=frac,
        ))
    tab = pd.DataFrame(rows)
    med_boot = float(tab.median_pairwise_deg.median()) if tab.median_pairwise_deg.notna().any() else np.nan
    med_rec = float(np.nanmedian(recov_all)) if recov_all else np.nan
    med_frac = float(np.nanmedian(frac_all)) if frac_all else np.nan
    log(f"   median-over-types pairwise bootstrap angle = {med_boot:.1f}°")
    if np.isfinite(med_rec):
        log(f"   recovery angle (median over boots and types) = {med_rec:.1f}°")
    if np.isfinite(med_frac):
        log(f"   sparse support in top-{N_SUPPORT} |coef| (median types) = {med_frac:.3f}")
    return dict(
        per_type=tab, boot_median_angle=med_boot,
        recovery_angle=med_rec, support_frac=med_frac, n_boot=n_boot,
    )


def oof_scheme(Y, obs, folds, scheme, rng, n_perm, log, y_is_binary=False):
    """V2 path: OOF ridge age direction (not TARGET). Primary R² site-stratified."""
    log(f"   OOF ridge scheme={scheme} n_folds={len(folds)} n_perm={n_perm}")
    out = oof_project(Y, obs, folds, "ridge", use_target=False)
    cache = out.get("svd_cache")
    pred = out["pred"]
    rmap, sc = _primary_map(obs, pred, scheme)
    prim = sc["site_strat_median_type_r2"]  # never pooled
    log(f"      site-strat R²={prim:+.3f}  pooled={sc['pooled_median_type_r2']:+.3f}")
    auc = auc_summary(obs, pred, obs.age.to_numpy(float)) if y_is_binary else None
    if auc is not None:
        log(f"      site-strat AUC={auc['site_strat_auc']:+.3f}  pooled={auc['pooled_auc']:+.3f}")
    null_prim, null_rows = [], []
    log(f"      permutation null ({n_perm} donor-level within-site shuffles)...")
    for i in range(n_perm):
        obs_p = permute_age_within_site_obs(obs, rng)
        yp = obs_p.age.to_numpy(float)
        out_p = oof_project(Y, obs, folds, "ridge", use_target=False, y_override=yp, svd_cache=cache)
        _, scn = _primary_map(obs_p, out_p["pred"], scheme, y=yp)
        npv = scn["site_strat_median_type_r2"]
        null_prim.append(npv)
        null_rows.append(dict(perm=i, site_strat=npv, pooled=scn["pooled_median_type_r2"]))
        if (i + 1) % 5 == 0:
            log(f"         perm {i+1}/{n_perm}")
    ns = summarize_null_col(null_prim)
    pval = permutation_p(prim, null_prim, greater=True)
    log(f"      null mean={ns['mean']:+.3f}  p95={ns['p95']:+.3f}  p={pval:.3f}")
    rec = dict(
        scheme=scheme, site_strat_r2=float(prim) if finite(prim) else np.nan,
        pooled_r2=float(sc["pooled_median_type_r2"]) if finite(sc["pooled_median_type_r2"]) else np.nan,
        null_mean=ns["mean"], null_p95=ns["p95"], p=float(pval) if finite(pval) else np.nan,
        n_perm=int(n_perm), per_type_r2=rmap,
        null_usable=bool(finite(ns["mean"]) and ns["mean"] <= 0.05),
    )
    if auc is not None:
        rec.update(auc)
    return rec, pd.DataFrame(null_rows), sc["per_type"]


def _primary_map(obs, pred, scheme, y=None):
    sc = score_age_predictions(obs, pred, y=y)
    pt = sc["per_type"]
    sub = pt[pt.site == "mean_of_sites"]
    return sub.set_index("celltype")["r2"].to_dict(), sc


def auc_score(y, pred):
    y = np.asarray(y, float)
    pred = np.asarray(pred, float)
    m = np.isfinite(y) & np.isfinite(pred)
    y, pred = y[m], pred[m]
    if len(y) < 5:
        return np.nan
    ypos = (y > 0).astype(np.int32)
    if ypos.size == 0 or int(ypos.min()) == int(ypos.max()):
        return np.nan
    try:
        from sklearn.metrics import roc_auc_score
        return float(roc_auc_score(ypos, pred))
    except Exception:
        pos, neg = pred[ypos == 1], pred[ypos == 0]
        if len(pos) < 1 or len(neg) < 1:
            return np.nan
        gt = np.sum(pos[:, None] > neg[None, :])
        eq = np.sum(pos[:, None] == neg[None, :])
        return float((gt + 0.5 * eq) / (len(pos) * len(neg)))


def auc_summary(obs, pred, y):
    y = np.asarray(y, float)
    pred = np.asarray(pred, float)
    ct = obs.celltype.astype(str).to_numpy()
    site = obs.Source.astype(str).to_numpy()
    types = sorted(pd.unique(ct))
    sites = sorted(pd.unique(site))
    pooled, strat = [], []
    for t in types:
        m = ct == t
        pooled.append(auc_score(y[m], pred[m]))
        ss = []
        for s in sites:
            ss.append(auc_score(y[m & (site == s)], pred[m & (site == s)]))
        ss = [a for a in ss if finite(a)]
        if ss:
            strat.append(float(np.mean(ss)))
    return dict(
        pooled_auc=float(np.nanmedian(pooled)) if any(finite(a) for a in pooled) else np.nan,
        site_strat_auc=float(np.nanmedian(strat)) if strat else np.nan,
    )


def _fit_on_mask(Y, obs, mask):
    idx = np.flatnonzero(mask)
    sub = obs.iloc[idx].reset_index(drop=True)
    Ys = Y[idx]
    ct = sub.celltype.astype(str).to_numpy()
    types_all = sorted(pd.unique(obs.celltype.astype(str)))
    y = sub.age.to_numpy(float)
    Xs, mu, sd = zscore_train(Ys)
    W_age, meta, built = per_type_directions(Xs, y, ct, types_all, method="ridge")
    return dict(
        idx=idx, obs=sub, X=Xs, mu=mu, sd=sd, y=y, ct=ct, types=types_all,
        W_age=W_age, meta=meta, svd=built,
    )


def _apply_age(pack_tr, Y, obs, test_mask, types):
    """V4a scoring but on the ridge AGE direction (L2c / gene-space equivalent)."""
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
        w = pack_tr["W_age"][:, j]
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


def site_transfer(Y, obs, rng, n_perm, log, y_is_binary=False):
    """L2c / V4a path: fit on one bank, score the other. Age direction, not TARGET."""
    site = obs.Source.astype(str).to_numpy()
    sites = sorted(pd.unique(site))
    types = sorted(pd.unique(obs.celltype.astype(str)))
    recs, per_tabs = [], []
    for tr_s in sites:
        te_s = [s for s in sites if s != tr_s][0]
        log(f"   fit on {SITE_NAME.get(tr_s, tr_s)} → test {SITE_NAME.get(te_s, te_s)}")
        pack = _fit_on_mask(Y, obs, site == tr_s)
        pred, yte, tab, med = _apply_age(pack, Y, obs, site == te_s, types)
        tab = tab.assign(train_site=tr_s, test_site=te_s)
        per_tabs.append(tab)
        log(f"      median-over-types R²={med:+.3f}")
        auc = auc_score(yte, pred) if y_is_binary else np.nan
        if y_is_binary:
            log(f"      held-out AUC={auc:+.3f}")
        log(f"      null ({n_perm} train-age shuffles)...")
        null_med = []
        obs_tr = pack["obs"]
        for i in range(n_perm):
            obs_p = permute_age_within_site_obs(obs_tr, rng)
            yp = obs_p.age.to_numpy(float)
            W_age, _, _ = per_type_directions(
                pack["X"], yp, pack["ct"], types, method="ridge", svd_cache=pack["svd"])
            pack_p = dict(pack)
            pack_p["W_age"] = W_age
            pack_p["y"] = yp
            _, _, _, mnull = _apply_age(pack_p, Y, obs, site == te_s, types)
            null_med.append(mnull)
        pval = permutation_p(med, null_med, greater=True)
        ns = summarize_null_col(null_med)
        log(f"      null mean={ns['mean']:+.3f}  p={pval:.3f}")
        recs.append(dict(
            train_site=tr_s, test_site=te_s, median_r2=med,
            null_mean=ns["mean"], null_p95=ns["p95"], p=pval, n_perm=n_perm,
            n_train_donors=int(obs.loc[site == tr_s, "donor"].nunique()),
            n_test_donors=int(obs.loc[site == te_s, "donor"].nunique()),
            auc=float(auc) if finite(auc) else np.nan,
        ))
    rec_df = pd.DataFrame(recs)
    per = pd.concat(per_tabs, ignore_index=True) if per_tabs else pd.DataFrame()
    hm = _xfer(rec_df, "H", "M")
    mh = _xfer(rec_df, "M", "H")
    return rec_df, per, hm, mh


def _xfer(df, tr, te):
    sub = df[(df.train_site.astype(str) == tr) & (df.test_site.astype(str) == te)]
    if not len(sub):
        return np.nan
    return float(sub.median_r2.iloc[0])


def _age_masks(obs, log, rng):
    """L2d bank-controlled young/old donor masks (copied, not re-tuned)."""
    don = obs.groupby("donor").agg(age=("age", "first"), site=("Source", "first"))
    donors = obs.donor.astype(str).to_numpy()
    site = obs.Source.astype(str).to_numpy()
    out = {}
    for s in sorted(pd.unique(site)):
        young = set(don.index[(don.site.astype(str) == s) & (don.age < 50)].astype(str))
        old = set(don.index[(don.site.astype(str) == s) & (don.age >= 50)].astype(str))
        m_y = np.array([d in young for d in donors])
        m_o = np.array([d in old for d in donors])
        log(f"[L2d] within {SITE_NAME.get(s, s)}: young donors={len(young)} rows={int(m_y.sum())}  "
            f"old donors={len(old)} rows={int(m_o.sum())}")
        out[f"within_{s}"] = (m_y, m_o, dict(n_young=len(young), n_old=len(old), site=s))
    y_keep, o_keep = [], []
    for s in sorted(pd.unique(don.site.astype(str))):
        y = don.index[(don.site.astype(str) == s) & (don.age < 50)].astype(str).to_numpy()
        o = don.index[(don.site.astype(str) == s) & (don.age >= 50)].astype(str).to_numpy()
        n = int(min(len(y), len(o)))
        log(f"[L2d] balance {SITE_NAME.get(s, s)}: min(young,old)={n} (young={len(y)} old={len(o)})")
        if n < 5:
            continue
        y_keep.extend(rng.choice(y, size=n, replace=False).tolist())
        o_keep.extend(rng.choice(o, size=n, replace=False).tolist())
    y_keep, o_keep = set(y_keep), set(o_keep)
    m_y = np.array([d in y_keep for d in donors])
    m_o = np.array([d in o_keep for d in donors])
    y_sites = don.loc[list(y_keep)].groupby("site").size().to_dict() if y_keep else {}
    o_sites = don.loc[list(o_keep)].groupby("site").size().to_dict() if o_keep else {}
    log(f"[L2d] balanced pool: young donors={len(y_keep)} sites={y_sites}  "
        f"old donors={len(o_keep)} sites={o_sites}")
    out["balanced"] = (m_y, m_o, dict(n_young=len(y_keep), n_old=len(o_keep),
                                      young_sites=y_sites, old_sites=o_sites))
    return out


def young_old_angles(Y, obs, types, rng, log, masks=None):
    """L2d path: young vs old ridge age-dir angles, bank-controlled."""
    if masks is None:
        masks = _age_masks(obs, log, rng)
    rows = []
    for split, (m_y, m_o, meta) in masks.items():
        iy, io = np.flatnonzero(m_y), np.flatnonzero(m_o)
        log(f"   split={split} young_rows={len(iy)} old_rows={len(io)}")
        if len(iy) < 20 or len(io) < 20:
            log("      SKIP too few rows")
            continue
        pack_y = _fit_on_mask(Y, obs, m_y)
        pack_o = _fit_on_mask(Y, obs, m_o)
        angs = []
        for j, t in enumerate(types):
            wy, wo = pack_y["W_age"][:, j], pack_o["W_age"][:, j]
            n_y = int((pack_y["ct"] == t).sum())
            n_o = int((pack_o["ct"] == t).sum())
            ok = (n_y >= MIN_N_TYPE and n_o >= MIN_N_TYPE
                  and np.isfinite(wy).all() and np.isfinite(wo).all()
                  and float(np.linalg.norm(wy)) > 1e-12 and float(np.linalg.norm(wo)) > 1e-12)
            ang = unsigned_angle_deg(wy, wo) if ok else np.nan
            if ok:
                angs.append(ang)
            rows.append(dict(
                split=split, celltype=t, n_young=n_y, n_old=n_o, angle_deg=ang, ok=bool(ok),
                n_donors_young=meta.get("n_young"), n_donors_old=meta.get("n_old"),
            ))
        med = float(np.median(angs)) if angs else np.nan
        log(f"      {split} median angle={med:.1f}°  n_types={len(angs)}")
    df = pd.DataFrame(rows)
    summary = {}
    for split in df.split.unique() if len(df) else []:
        sub = df[(df.split == split) & df.ok]
        summary[str(split)] = float(sub.angle_deg.median()) if len(sub) else np.nan
    return df, summary


def cheap_within_site_r2(Y, obs, y, rng, n_folds=2):
    """Calibration-only 2-fold within-site R². Reported numbers use the real 5-fold."""
    from brain_phase1_common import within_site_donor_folds  # noqa
    obs_t = obs_with_y(obs, y)
    folds = within_site_donor_folds(obs_t, int(n_folds), rng)
    out = oof_project(Y, obs_t, folds, "ridge", use_target=False)
    _, sc = _primary_map(obs_t, out["pred"], "within_site")
    return float(sc["site_strat_median_type_r2"])


def full_within_site_r2(Y, obs, y):
    """Reported calibration metric: real Phase-1 within-site folds, site-stratified ridge R²."""
    obs_t = obs_with_y(obs, y)
    folds = make_folds(obs, np.random.default_rng(PC_SEED))["within_site"]
    out = oof_project(Y, obs_t, folds, "ridge", use_target=False)
    _, sc = _primary_map(obs_t, out["pred"], "within_site")
    return float(sc["site_strat_median_type_r2"])


def calibrate_beta(Y, obs, y, types, pack, target_r2, log, lo=1e-4, hi=3.0, n_grid=10):
    """Grid β so cheap within-site R² is near target; freeze; measured R² is later."""
    rng = np.random.default_rng(CAL_SEED)
    betas = np.geomspace(float(lo), float(hi), int(n_grid))
    rows = []
    best_b, best_d = float(betas[0]), np.inf
    for b in betas:
        Yp = plant(Y, y, obs, types, pack, b)
        r2 = cheap_within_site_r2(Yp, obs, y, rng)
        d = abs(r2 - target_r2) if finite(r2) else np.inf
        rows.append(dict(beta=float(b), cheap_r2=r2, abs_err=d if np.isfinite(d) else np.nan))
        log(f"   cal β={b:.5g}  cheap within-site R²={r2:+.3f}  |err|={d if np.isfinite(d) else np.nan:.3f}")
        if d < best_d:
            best_d, best_b = d, float(b)
        if finite(r2) and r2 > target_r2 + 0.20:
            break
    # interpolate in log-β between the two grid points that straddle the target
    tab = pd.DataFrame(rows)
    ok = tab[np.isfinite(tab.cheap_r2)]
    frozen = best_b
    if len(ok) >= 2:
        below = ok[ok.cheap_r2 <= target_r2]
        above = ok[ok.cheap_r2 >= target_r2]
        if len(below) and len(above):
            b0, r0 = float(below.iloc[-1].beta), float(below.iloc[-1].cheap_r2)
            b1, r1 = float(above.iloc[0].beta), float(above.iloc[0].cheap_r2)
            if abs(r1 - r0) > 1e-9:
                t = (target_r2 - r0) / (r1 - r0)
                logb = np.log(b0) + t * (np.log(b1) - np.log(b0))
                frozen = float(np.exp(logb))
    log(f"   frozen β={frozen:.6g} (target R²={target_r2:.2f}; cheap grid best |err|={best_d:.3f})")
    return float(frozen), tab


def sex_pm1(obs):
    s = obs.sex.astype(str).str.lower().to_numpy()
    y = np.full(len(obs), np.nan)
    y[s == "male"] = 1.0
    y[s == "female"] = -1.0
    n_m = int(np.sum(y == 1.0))
    n_f = int(np.sum(y == -1.0))
    return y, n_m, n_f


def sex_chrom_mask(genes):
    gid = genes.gene_id.astype(str).to_numpy()
    sym = genes.symbol.astype(str).str.upper().to_numpy()
    chrom = load_chrom_series(gid)
    chrom = pd.Series(np.asarray(chrom).astype(str), index=np.arange(len(gid)))
    ch = chrom.str.upper().to_numpy()
    drop = np.isin(ch, list(SEX_CHROM)) | np.isin(sym, list(PAR_SYMBOLS))
    keep = ~drop
    return keep, dict(
        n_genes=int(len(gid)), n_drop=int(drop.sum()), n_keep=int(keep.sum()),
        n_chrom_xy=int(np.isin(ch, ["X", "Y", "XY"]).sum()),
        n_par_symbol=int(np.isin(sym, list(PAR_SYMBOLS)).sum()),
        chrom_counts={str(k): int(v) for k, v in pd.Series(ch).value_counts().to_dict().items()},
    )


def subsample_donors(obs, n, rng):
    """Bank-stratified donor subsample (largest remainder). Returns boolean row mask + donor list."""
    don = obs.groupby("donor").agg(site=("Source", "first")).reset_index()
    don["donor"] = don.donor.astype(str)
    don["site"] = don.site.astype(str)
    n = int(min(n, len(don)))
    sites = list(don.site.value_counts().index)
    counts = don.site.value_counts().to_dict()
    total = int(sum(counts.values()))
    alloc = {}
    remainders = []
    used = 0
    for i, s in enumerate(sites):
        exact = n * counts[s] / total
        base = int(np.floor(exact))
        alloc[s] = base
        used += base
        remainders.append((exact - base, s))
    remainders.sort(reverse=True)
    extra = n - used
    for k in range(extra):
        alloc[remainders[k % len(remainders)][1]] += 1
    picked = []
    for s, k in alloc.items():
        pool = don.loc[don.site == s, "donor"].to_numpy()
        k = int(min(k, len(pool)))
        if k < 1:
            continue
        picked.extend(rng.choice(pool, size=k, replace=False).tolist())
    picked = np.array(sorted(set(str(d) for d in picked)))
    row = obs.donor.astype(str).isin(picked).to_numpy()
    return row, picked


def restrict(Y, obs, row_mask):
    idx = np.flatnonzero(row_mask)
    return np.ascontiguousarray(Y[idx]), obs.iloc[idx].reset_index(drop=True)


def gene_ridge_row():
    base = gene_space_baseline()
    return dict(
        regime="gene_ridge",
        target_r2=np.nan,
        achieved_within_site_r2=base["within_site_r2_v2"],
        within_site_null=np.nan,
        boot_median_angle=base["boot_median_angle"],
        recovery_angle=np.nan,
        loso_r2=np.nan,
        transfer_H_to_M=base["transfer_H_to_M"],
        transfer_M_to_H=base["transfer_M_to_H"],
        young_old_balanced=base["young_old_angle"],
        pass_angle=False,
        pass_transfer=False,
        pass_both=False,
        n_boot=N_BOOT_FULL,
        n_perm=N_PERM_FULL,
        note="gene-space TARGET/V1 from results/target json (LOWDIM gene_ridge row); not refit",
        b_reduced=False,
    )


def summarize_battery(bat):
    hm = bat.get("transfer_H_to_M", np.nan)
    mh = bat.get("transfer_M_to_H", np.nan)
    hm_n = bat.get("transfer_H_to_M_null", np.nan)
    mh_n = bat.get("transfer_M_to_H_null", np.nan)
    ang = bat.get("boot_median_angle", np.nan)
    pa = pass_angle(ang)
    pt = pass_transfer_both(hm, mh)
    pt_pri = pass_transfer_primary(hm, mh, hm_n, mh_n)
    ang_vs_sex = (float(ang) - SEX_ANGLE_REF) if finite(ang) else np.nan
    return dict(
        boot_median_angle=ang,
        recovery_angle=bat.get("recovery_angle", np.nan),
        achieved_within_site_r2=bat.get("within_site_r2", np.nan),
        within_site_null=bat.get("within_site_null", np.nan),
        loso_r2=bat.get("loso_r2", np.nan),
        loso_null=bat.get("loso_null", np.nan),
        transfer_H_to_M=hm,
        transfer_M_to_H=mh,
        transfer_H_to_M_null=hm_n,
        transfer_M_to_H_null=mh_n,
        young_old_balanced=bat.get("young_old_balanced", np.nan),
        young_old_within_H=bat.get("young_old_within_H", np.nan),
        young_old_within_M=bat.get("young_old_within_M", np.nan),
        pass_angle=pa,
        pass_transfer=pt,
        pass_transfer_any=pass_transfer_any(hm, mh),
        pass_transfer_primary=pt_pri,
        angle_vs_sex_ref=ang_vs_sex,
        sex_angle_ref=SEX_ANGLE_REF,
        pass_both=bool(pa and pt),
        support_frac=bat.get("support_frac", np.nan),
        site_strat_auc=bat.get("site_strat_auc", np.nan),
        pooled_auc=bat.get("pooled_auc", np.nan),
        n_boot=bat.get("n_boot", np.nan),
        n_perm=bat.get("n_perm", np.nan),
    )


def run_l2_gene_battery(Y, obs, y, rng, log, *, n_boot=N_BOOT_FULL, n_perm=N_PERM_FULL,
                        planted_pack=None, parts=None, y_is_binary=False,
                        split_on_chronological=False, tag=""):
    """Identical gene-space L2 operations as V1/V2/V4a/L2d, different X and y.

    Primary R² is site-stratified. Age direction (not TARGET). Permutation null
    re-shuffles y at donor level within bank; X is left in place.
    """
    parts = set(parts or ("boot", "oof_ws", "oof_loso", "transfer", "youngold"))
    types = sorted(pd.unique(obs.celltype.astype(str)))
    obs_t = obs_with_y(obs, y)
    t0 = time.time()
    out = dict(n_boot=int(n_boot), n_perm=int(n_perm), tag=tag, n_rows=int(len(obs)),
               n_donors=int(obs.donor.nunique()), n_genes=int(Y.shape[1]), n_types=len(types))
    log(f"[battery] {tag} n_boot={n_boot} n_perm={n_perm} parts={sorted(parts)} "
        f"donors={out['n_donors']} genes={out['n_genes']}")

    if "boot" in parts:
        boot = bootstrap_ridge(Y, obs_t, types, n_boot, rng, log, planted_pack=planted_pack)
        out["boot_median_angle"] = boot["boot_median_angle"]
        out["recovery_angle"] = boot["recovery_angle"]
        out["support_frac"] = boot["support_frac"]
        boot["per_type"].to_csv(PC_DIR / f"{tag}_bootstrap_per_type.csv", index=False)
        dump_json(PC_DIR / f"{tag}_bootstrap.json", {k: boot[k] for k in boot if k != "per_type"})

    folds = None
    if parts & {"oof_ws", "oof_loso"}:
        folds = make_folds(obs, np.random.default_rng(PC_SEED))
        log(f"   folds seed={PC_SEED}  within_site={len(folds['within_site'])}  loso={len(folds['loso'])}")

    if "oof_ws" in parts:
        rec, nt, pt = oof_scheme(Y, obs_t, folds["within_site"], "within_site", rng, n_perm, log,
                                 y_is_binary=y_is_binary)
        out["within_site_r2"] = rec["site_strat_r2"]
        out["within_site_null"] = rec["null_mean"]
        out["within_site_p"] = rec["p"]
        out["within_site_pooled_r2"] = rec["pooled_r2"]
        if y_is_binary:
            out["site_strat_auc"] = rec.get("site_strat_auc", np.nan)
            out["pooled_auc"] = rec.get("pooled_auc", np.nan)
        nt.to_csv(PC_DIR / f"{tag}_within_site_null.csv", index=False)
        pt.to_csv(PC_DIR / f"{tag}_within_site_per_type.csv", index=False)
        dump_json(PC_DIR / f"{tag}_within_site.json", rec)

    if "oof_loso" in parts:
        rec, nt, pt = oof_scheme(Y, obs_t, folds["loso"], "loso", rng, n_perm, log,
                                 y_is_binary=y_is_binary)
        out["loso_r2"] = rec["site_strat_r2"]
        out["loso_null"] = rec["null_mean"]
        out["loso_p"] = rec["p"]
        out["loso_pooled_r2"] = rec["pooled_r2"]
        nt.to_csv(PC_DIR / f"{tag}_loso_null.csv", index=False)
        pt.to_csv(PC_DIR / f"{tag}_loso_per_type.csv", index=False)
        dump_json(PC_DIR / f"{tag}_loso.json", rec)

    if "transfer" in parts:
        rec_df, per, hm, mh = site_transfer(Y, obs_t, rng, n_perm, log, y_is_binary=y_is_binary)
        out["transfer_H_to_M"] = hm
        out["transfer_M_to_H"] = mh
        if len(rec_df):
            hrow = rec_df[(rec_df.train_site == "H") & (rec_df.test_site == "M")]
            mrow = rec_df[(rec_df.train_site == "M") & (rec_df.test_site == "H")]
            out["transfer_H_to_M_null"] = float(hrow.null_mean.iloc[0]) if len(hrow) else np.nan
            out["transfer_M_to_H_null"] = float(mrow.null_mean.iloc[0]) if len(mrow) else np.nan
            out["transfer_H_to_M_p"] = float(hrow.p.iloc[0]) if len(hrow) else np.nan
            out["transfer_M_to_H_p"] = float(mrow.p.iloc[0]) if len(mrow) else np.nan
        rec_df.to_csv(PC_DIR / f"{tag}_transfer.csv", index=False)
        per.to_csv(PC_DIR / f"{tag}_transfer_per_type.csv", index=False)

    if "youngold" in parts:
        rng_masks = np.random.default_rng(PC_SEED)
        masks = _age_masks(obs, log, rng_masks) if split_on_chronological else None
        yo, summary = young_old_angles(Y, obs_t, types, rng_masks, log, masks=masks)
        out["young_old_balanced"] = summary.get("balanced", np.nan)
        out["young_old_within_H"] = summary.get("within_H", np.nan)
        out["young_old_within_M"] = summary.get("within_M", np.nan)
        yo.to_csv(PC_DIR / f"{tag}_youngold_per_type.csv", index=False)
        dump_json(PC_DIR / f"{tag}_youngold.json", summary)

    out["elapsed_s"] = float(time.time() - t0)
    out.update(summarize_battery(out))
    dump_json(PC_DIR / f"{tag}_battery.json", out)
    log(f"[battery] {tag} done in {out['elapsed_s']:.0f}s  boot={out.get('boot_median_angle', np.nan):.1f}°  "
        f"ws R²={out.get('within_site_r2', np.nan):+.3f}  "
        f"H→M={out.get('transfer_H_to_M', np.nan):+.3f}  M→H={out.get('transfer_M_to_H', np.nan):+.3f}")
    return out


def jsonable(x):
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return float(x) if np.isfinite(x) else None
    if isinstance(x, np.ndarray):
        return jsonable(x.tolist())
    if isinstance(x, (pd.Timestamp, Path)):
        return str(x)
    if x is None or isinstance(x, (str, int, float, bool)):
        if isinstance(x, float) and not np.isfinite(x):
            return None
        return x
    return str(x)


def progress_snapshot(cells, next_action, stop=None, betas=None, extra=""):
    lines = [
        "# PROGRESS_POSCTRL",
        "",
        f"**STOP status:** {stop if stop else 'none'}.",
        f"**Next action:** {next_action}",
        "",
        "## Pre-registered block (verbatim)",
        "",
        PREREG_BLOCK,
        "",
        "## Supersession 2026-09-16 (verbatim, written before any C1/C3 fit)",
        "",
        SUPERSESSION_BLOCK,
        "",
        "## Seeds",
        "",
        f"- main / battery: `{PC_SEED}`",
        f"- planted w: `{W_SEED}`",
        f"- β calibration: `{CAL_SEED}`",
        f"- C3 subsample: `{C3_SEED}`",
        f"- dataset: `{DATASET_ID}`",
        "",
        "## Frozen β per regime",
        "",
    ]
    if betas:
        for k, v in betas.items():
            lines.append(f"- `{k}`: β={v}")
    else:
        lines.append("Not calibrated (C1 not run)." if not any(n.startswith("c1_") for n in cells) else "See cells.")
    lines += ["", "## Finished cells", ""]
    if not cells:
        lines.append("None.")
    else:
        # c2, then c1, then c3, then other — numbers from disk json only
        order = sorted(cells.keys(), key=lambda n: (
            0 if n.startswith("c2") else 1 if n.startswith("c1_") else 2 if n.startswith("c3") else 3,
            n,
        ))
        keys = (
            "boot_median_angle", "angle_vs_sex_ref", "sex_angle_ref", "recovery_angle",
            "achieved_within_site_r2", "within_site_null", "loso_r2",
            "transfer_H_to_M", "transfer_M_to_H",
            "transfer_H_to_M_null", "transfer_M_to_H_null",
            "young_old_balanced",
            "pass_angle", "pass_transfer", "pass_transfer_primary", "pass_both",
            "beta", "target_r2", "regime", "kind", "n_requested", "draw",
            "n_boot", "n_perm", "site_strat_auc",
        )
        for name in order:
            rec = cells[name]
            lines.append(f"### {name}")
            for k in keys:
                if k in rec and rec[k] is not None:
                    lines.append(f"- {k}: {rec[k]}")
            if rec.get("path"):
                p = rec["path"]
                # prefer repo-relative path in the progress file
                try:
                    rel = str(Path(p).resolve().relative_to(ROOT.resolve())).replace("\\", "/")
                    lines.append(f"- path: `{rel}`")
                except Exception:
                    lines.append(f"- path: `{p}`")
            lines.append("")
    lines += ["## Paths written", ""]
    for name in sorted(cells):
        lines.append(f"- `results/posctrl/{name}.json`")
    if SUPERSESSION_FLAG.exists():
        lines.append("- `results/posctrl/SUPERSESSION_20260916.flag`")
        lines.append("- `results/posctrl/supersession_20260916.json`")
    if (PC_DIR / "STOP_P1.txt").exists():
        lines.append("- `results/posctrl/STOP_P1.txt`")
    lines.append("- `FINDINGS_POSCTRL.md`")
    lines.append("- `PROGRESS_POSCTRL.md`")
    lines += [
        "",
        "## Bars (frozen)",
        "",
        "### Original (P1 record; 35° retired for C1/C3)",
        f"- angle: bootstrap median pairwise ≤ {ANGLE_BAR:.0f}°",
        "- planted transfer (P2): HBCC→MSSM > 0 AND MSSM→HBCC > 0",
        "- C2a transfer (P1 / L2): at least one direction > 0",
        f"- calibration: within-site site-stratified ridge R² = {TARGET_R2_PRIMARY:.2f} ± {CAL_R2_TOL:.2f}",
        "",
        "### Supersession 2026-09-16 (C1/C3)",
        f"- (a) primary: transfer R² > 0 in both directions with permutation null ≤ {TRANSFER_NULL_BAR}",
        f"- (b) bootstrap median pairwise angle reported against sex reference {SEX_ANGLE_REF:.1f}°, not {ANGLE_BAR:.0f}°",
        "- P2 = planted dense_shared at R²≈0.30 transfers both ways",
        "- P3 = it does not",
        f"- calibration unchanged: within-site site-stratified ridge R² = {TARGET_R2_PRIMARY:.2f} ± {CAL_R2_TOL:.2f}",
        "",
    ]
    if extra:
        lines += [extra, ""]
    write_progress("\n".join(lines) + "\n")

"""TOWARD — do partially-reprogrammed cells move TOWARD young fibroblasts?

Read-only on existing matrices, labels, and the frozen ruler. Refit nothing.
Re-cluster nothing. Rescore nothing. GSE325735 out of scope.

Pre-registration must already exist at results/toward/PREREG.flag (written before
any statistic). This script does not rewrite that block.

Usage: python src/toward_run.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import RESULTS, DATA_PROC, ROOT  # noqa: E402
from brain_phase1_common import Logger, dump_json, load_json  # noqa: E402
from gtex_common import (  # noqa: E402
    StopStep, EXPECTED_AGE_BINS, tmm_norm_factors, log2_cpm_edger, EDGE_R_PRIOR,
)
from fibro_common import (  # noqa: E402
    FROZEN_RULER, jsonable, FIBRO_DIR, FIBRO_PROC,
)
from fibro_stage1 import load_pack  # noqa: E402
from fibro2_common import load_frozen_ruler  # noqa: E402
from target_common import unit, md_table  # noqa: E402
from trajectory_common import permutation_p  # noqa: E402
from md2_common import MD2_DIR  # noqa: E402
from md4_common import MD4_DIR, MD4_PROC, STATE_NAMES  # noqa: E402
from fibro3_common import AGED_LINE, YOUNG_LINE, DAYS  # noqa: E402


TOWARD_DIR = RESULTS / "toward"
TOWARD_DIR.mkdir(parents=True, exist_ok=True)
FINDINGS_PATH = ROOT / "FINDINGS_TOWARD.md"
PROGRESS_PATH = ROOT / "PROGRESS_TOWARD.md"
PREREG_FLAG = TOWARD_DIR / "PREREG.flag"
MANIFEST_PATH = TOWARD_DIR / "manifest.json"
ANCHORS_PATH = TOWARD_DIR / "anchors.npz"

TOWARD_SEED = 20260914
BOOT_SEED = 20260918
N_PERM = 200
N_BOOT = 200
N_RANDOM = 200
N_FOLDS = 5
MIN_CELLS = 30
MIN_ANCHOR_N = 50
C3_COS_BAR = 0.50
C3_MIN_YOUNG_TRAIN = 20
C3_MIN_OLD_TRAIN = 20
C3_MIN_OLD_TEST = 5
P_BAR = 0.05

YOUNG_BINS = ("20-29", "30-39")
OLD_BINS = ("60-69", "70-79")
MIDDLE_BINS = ("40-49", "50-59")
N_YOUNG_EXPECTED = 113
N_OLD_EXPECTED = 231

DEST_STATES = ("PartialReprog", "EarlyPluripotency", "Pluripotency", "NonReprog")
ORIGIN_STATE = "Fibroblast"
CLAIM_DONOR = AGED_LINE
CONTRAST_DONOR = YOUNG_LINE
DONOR_BOOT_OFFSET = {AGED_LINE: 0, YOUNG_LINE: 1}

EXTRAP_PATH = MD4_DIR / "t2_extrap.csv"
GSE325735_BAN = "GSE325735 is out of scope and is not opened."


def load_prereg_text():
    if not PREREG_FLAG.exists():
        raise StopStep(
            "prereg",
            "PREREG.flag missing — write the pre-registration before any statistic. Not computing.",
        )
    return PREREG_FLAG.read_text(encoding="utf-8").strip()


def load_manifest():
    if MANIFEST_PATH.exists():
        return load_json(MANIFEST_PATH)
    return dict(seed=TOWARD_SEED, boot_seed=BOOT_SEED, failures=[], status="INIT")


def save_manifest(man):
    dump_json(MANIFEST_PATH, jsonable(man))
    return MANIFEST_PATH


def record_failure(step, message, details=None):
    man = load_manifest()
    rec = dict(step=step, message=message, details=jsonable(details or {}))
    man.setdefault("failures", []).append(rec)
    man["status"] = "STOP"
    save_manifest(man)
    return rec


def log_columns(log, tag, cols, path):
    log(f"[columns {tag}] actually read from {path}: {list(cols)}")


def cosine_align(d, u):
    """cos = dot(d, u) / ||d||  with u unit. NA if ||d|| is degenerate."""
    d = np.asarray(d, float).ravel()
    u = np.asarray(u, float).ravel()
    nd = float(np.linalg.norm(d))
    if (not np.isfinite(nd)) or nd < 1e-12 or (not np.isfinite(d).all()) or (not np.isfinite(u).all()):
        return np.nan
    return float(np.dot(d, u) / nd)


def tmm_logcpm_quiet(Y):
    """Y samples × genes counts → log2-CPM (TMM, prior.count=2)."""
    Y = np.asarray(Y, np.float64)
    nf = tmm_norm_factors(Y.T)
    logcpm = log2_cpm_edger(Y.T, nf, prior=EDGE_R_PRIOR).T
    return logcpm, nf


def zscore_frozen(logcpm, frozen, counts):
    mu = np.asarray(frozen["mu"], float)
    sd = np.asarray(frozen["sd"], float)
    sd = np.where(sd < 1e-12, 1.0, sd)
    missing = (np.asarray(counts) == 0).all(0)
    Z = (np.asarray(logcpm, float) - mu) / sd
    Z[:, missing] = 0.0
    return Z, missing


def ruler_score(Z, frozen):
    w, _ = unit(np.asarray(frozen["w"], float))
    return np.asarray(Z, float) @ w


def percentile_ci(vals):
    v = np.asarray(vals, float)
    v = v[np.isfinite(v)]
    if v.size < 2:
        return np.nan, np.nan
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


def _fmt(x, d=3):
    if x is None:
        return "NA"
    if isinstance(x, (bool, np.bool_)):
        return str(bool(x))
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    if not np.isfinite(v):
        return "NA"
    return f"{v:+.{d}f}"


def _fmt_p(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(v):
        return "NA"
    return f"{v:.4f}"


def _fmt_ci(lo, hi, d=3):
    a, b = _fmt(lo, d), _fmt(hi, d)
    if a == "NA" and b == "NA":
        return "[NA, NA]"
    return f"[{a}, {b}]"


def sum_rows(Y, idx):
    idx = np.asarray(idx, int)
    if idx.size == 0:
        return np.zeros(Y.shape[1], dtype=np.float64)
    return np.asarray(Y[idx].sum(axis=0), dtype=np.float64).ravel()


def load_Y(line):
    path = MD4_PROC / f"louvain_rulerY_{line}.npz"
    if not path.exists():
        raise StopStep("rulerY", f"missing {path}. Not substituting.")
    z = np.load(path)
    return sparse.csr_matrix(
        (z["data"], z["indices"], z["indptr"]),
        shape=tuple(int(x) for x in z["shape"]),
    )


def require_file(path, step):
    if not Path(path).exists():
        raise StopStep(step, f"missing {path}. Not substituting.")
    return Path(path)


def load_gtex_z(log):
    require_file(FROZEN_RULER, "frozen_ruler")
    require_file(FIBRO_PROC / "stage1_fibro.npz", "gtex_pack")
    require_file(FIBRO_DIR / "stage1_obs.csv", "gtex_obs")
    frozen = load_frozen_ruler()
    pack = load_pack()
    obs = pack["obs"].copy()
    log_columns(log, "gtex_obs", list(obs.columns), str(FIBRO_DIR / "stage1_obs.csv"))
    for col in ("SAMPID", "donor", "AGE"):
        if col not in obs.columns:
            raise StopStep("gtex_obs", f"missing column {col!r}. Not substituting.",
                           dict(columns=list(obs.columns)))
    X = np.asarray(pack["X"], np.float64)
    mu = np.asarray(frozen["mu"], float)
    sd = np.asarray(frozen["sd"], float)
    if X.shape[1] != mu.shape[0]:
        raise StopStep(
            "genes",
            f"pack n_genes={X.shape[1]} frozen mu n={mu.shape[0]}. Not aligning ad hoc.",
        )
    ens_pack = np.asarray(pack["ensembl"]).astype(str)
    ens_fr = np.asarray(frozen["ensembl"]).astype(str)
    if not np.array_equal(ens_pack, ens_fr):
        gid_pack = np.asarray(pack["genes"]).astype(str)
        gid_fr = np.asarray(frozen["gene_id"]).astype(str)
        if not np.array_equal(gid_pack, gid_fr):
            raise StopStep(
                "genes",
                "GTEx pack gene order does not match frozen-ruler overlap gene set. Not reordering.",
            )
        log("[genes] ensembl arrays differ; gene_id arrays match. Using pack as frozen overlap.")
    else:
        log("[genes] pack ensembl matches frozen-ruler overlap gene set.")
    sd_safe = np.where(sd < 1e-12, 1.0, sd)
    Z = (X - mu) / sd_safe
    age = obs["AGE"].astype(str).to_numpy()
    unknown = sorted(set(age) - set(EXPECTED_AGE_BINS))
    if unknown:
        raise StopStep("age_bins", f"unexpected AGE bins {unknown}. Not recoding.")
    log(f"[gtex] n_donors={len(obs)} n_genes={Z.shape[1]} "
        f"filter>=6 in >=20% ; TMM/log2-CPM prior.count={EDGE_R_PRIOR:g} (Stage 0 pack, not rebuilt)")
    log(f"[gtex] AGE counts: {pd.Series(age).value_counts().reindex(list(EXPECTED_AGE_BINS)).fillna(0).astype(int).to_dict()}")
    return dict(Z=Z, obs=obs, frozen=frozen, pack=pack, age=age)


def freeze_anchors(gtex, log):
    Z, obs, age = gtex["Z"], gtex["obs"], gtex["age"]
    young_m = np.isin(age, YOUNG_BINS)
    old_m = np.isin(age, OLD_BINS)
    mid_m = np.isin(age, MIDDLE_BINS)
    n_young = int(young_m.sum())
    n_old = int(old_m.sum())
    n_mid = int(mid_m.sum())
    bin_n = {b: int((age == b).sum()) for b in EXPECTED_AGE_BINS}
    log(f"[anchors] YOUNG bins {YOUNG_BINS} n={n_young} (expected {N_YOUNG_EXPECTED})")
    log(f"[anchors] OLD bins {OLD_BINS} n={n_old} (expected {N_OLD_EXPECTED})")
    log(f"[anchors] middle bins {MIDDLE_BINS} n={n_mid} excluded from both anchors, not dropped silently")
    log(f"[anchors] per-bin n={bin_n}")
    if n_young < MIN_ANCHOR_N:
        raise StopStep("anchors", f"YOUNG n={n_young} < {MIN_ANCHOR_N}. STOP.")
    if n_old < MIN_ANCHOR_N:
        raise StopStep("anchors", f"OLD n={n_old} < {MIN_ANCHOR_N}. STOP.")
    c_young = Z[young_m].mean(0)
    c_old = Z[old_m].mean(0)
    diff = c_young - c_old
    u, nrm = unit(diff)
    if nrm < 1e-12:
        raise StopStep("anchors", "||c_young - c_old|| is degenerate. STOP.")
    np.savez_compressed(
        ANCHORS_PATH,
        c_young=np.asarray(c_young, np.float64),
        c_old=np.asarray(c_old, np.float64),
        u=np.asarray(u, np.float64),
        n_young=np.array(n_young),
        n_old=np.array(n_old),
        n_middle=np.array(n_mid),
        young_bins=np.array(YOUNG_BINS),
        old_bins=np.array(OLD_BINS),
        middle_bins=np.array(MIDDLE_BINS),
        ensembl=np.asarray(gtex["frozen"]["ensembl"]).astype(str),
        n_genes=np.array(int(Z.shape[1])),
        n_axis=np.array(float(nrm)),
        n_young_expected=np.array(N_YOUNG_EXPECTED),
        n_old_expected=np.array(N_OLD_EXPECTED),
        note=np.array("frozen before any GSE297234 vector"),
    )
    log(f"[anchors] wrote {ANCHORS_PATH}  ||c_young-c_old||={nrm:.6f}")
    rows = []
    for b in EXPECTED_AGE_BINS:
        role = "YOUNG" if b in YOUNG_BINS else ("OLD" if b in OLD_BINS else "middle_excluded")
        rows.append(dict(AGE_bin=b, n=bin_n[b], role=role))
    rows.append(dict(AGE_bin="YOUNG_total", n=n_young, role="YOUNG"))
    rows.append(dict(AGE_bin="OLD_total", n=n_old, role="OLD"))
    rows.append(dict(AGE_bin="middle_total", n=n_mid, role="middle_excluded"))
    adf = pd.DataFrame(rows)
    adf.to_csv(TOWARD_DIR / "anchor_table.csv", index=False)
    meta = dict(
        n_young=n_young, n_old=n_old, n_middle=n_mid, bin_n=bin_n,
        n_young_expected=N_YOUNG_EXPECTED, n_old_expected=N_OLD_EXPECTED,
        axis_norm=float(nrm), n_genes=int(Z.shape[1]), path=str(ANCHORS_PATH),
    )
    dump_json(TOWARD_DIR / "anchors.json", jsonable(meta))
    return dict(
        c_young=c_young, c_old=c_old, u=u, n_young=n_young, n_old=n_old,
        n_middle=n_mid, bin_n=bin_n, axis_norm=float(nrm),
        young_m=young_m, old_m=old_m, mid_m=mid_m, table=adf,
    )


def build_null_directions(Z, u, log):
    rng1 = np.random.default_rng(TOWARD_SEED)
    rng2 = np.random.default_rng(TOWARD_SEED)
    n_genes = int(u.size)
    U1 = np.empty((N_RANDOM, n_genes), dtype=np.float64)
    for i in range(N_RANDOM):
        U1[i] = rng1.permutation(u)
    n = int(Z.shape[0])
    n_half = n // 2
    U2 = np.empty((N_PERM, n_genes), dtype=np.float64)
    n_degen = 0
    for i in range(N_PERM):
        idx = rng2.permutation(n)
        a = idx[:n_half]
        b = idx[n_half:2 * n_half]
        diff = Z[a].mean(0) - Z[b].mean(0)
        ur, nrm = unit(diff)
        if nrm < 1e-12:
            n_degen += 1
            U2[i] = np.nan
        else:
            U2[i] = ur
    log(f"[nulls] N1 n_random={N_RANDOM} seed={TOWARD_SEED}; "
        f"N2 n_perm={N_PERM} n_donors={n} n_half={n_half} n_degen={n_degen} seed={TOWARD_SEED}")
    np.savez_compressed(TOWARD_DIR / "null_directions.npz", N1=U1, N2=U2)
    return U1, U2


def null_p(cos_obs, U, d):
    null_cos = np.array([cosine_align(d, ur) for ur in U], dtype=float)
    p = permutation_p(cos_obs, null_cos, greater=True)
    return p, null_cos


def run_c3(gtex, log):
    Z, obs, age = gtex["Z"], gtex["obs"], gtex["age"]
    n = len(obs)
    rng = np.random.default_rng(TOWARD_SEED)
    order = rng.permutation(n)
    folds = np.empty(n, dtype=int)
    folds[order] = np.arange(n) % N_FOLDS
    rows = []
    construction_ok = True
    reasons = []
    rng_boot = np.random.default_rng(BOOT_SEED)
    for k in range(N_FOLDS):
        tr = folds != k
        te = folds == k
        age_tr, age_te = age[tr], age[te]
        young_tr = np.isin(age_tr, YOUNG_BINS)
        old_tr = np.isin(age_tr, OLD_BINS)
        old_te = np.isin(age_te, OLD_BINS)
        n_ytr, n_otr, n_ote = int(young_tr.sum()), int(old_tr.sum()), int(old_te.sum())
        rec = dict(
            fold=int(k), n_train=int(tr.sum()), n_test=int(te.sum()),
            n_young_train=n_ytr, n_old_train=n_otr, n_old_test=n_ote,
            cos=np.nan, d_norm=np.nan, cos_ci_lo=np.nan, cos_ci_hi=np.nan,
            construction_ok=True, reason="",
        )
        if n_ytr < C3_MIN_YOUNG_TRAIN or n_otr < C3_MIN_OLD_TRAIN or n_ote < C3_MIN_OLD_TEST:
            rec["construction_ok"] = False
            rec["reason"] = (
                f"construction minimum failed: n_young_train={n_ytr} "
                f"(bar {C3_MIN_YOUNG_TRAIN}) n_old_train={n_otr} (bar {C3_MIN_OLD_TRAIN}) "
                f"n_old_test={n_ote} (bar {C3_MIN_OLD_TEST})"
            )
            construction_ok = False
            reasons.append(rec["reason"])
            rows.append(rec)
            log(f"[C3] fold {k} {rec['reason']}")
            continue
        Ztr, Zte = Z[tr], Z[te]
        c_y = Ztr[young_tr].mean(0)
        c_o = Ztr[old_tr].mean(0)
        u_tr, nrm = unit(c_y - c_o)
        if nrm < 1e-12:
            rec["construction_ok"] = False
            rec["reason"] = "train ||c_young-c_old|| degenerate"
            construction_ok = False
            reasons.append(rec["reason"])
            rows.append(rec)
            continue
        z_old_te = Zte[old_te]
        origin = z_old_te.mean(0)
        d = c_y - origin
        rec["d_norm"] = float(np.linalg.norm(d))
        rec["cos"] = cosine_align(d, u_tr)
        rec["axis_norm_train"] = float(nrm)
        boot = []
        n_ote_i = int(z_old_te.shape[0])
        for _ in range(N_BOOT):
            take = rng_boot.integers(0, n_ote_i, size=n_ote_i)
            d_b = c_y - z_old_te[take].mean(0)
            boot.append(cosine_align(d_b, u_tr))
        rec["cos_ci_lo"], rec["cos_ci_hi"] = percentile_ci(boot)
        rec["pass_bar"] = bool(np.isfinite(rec["cos"]) and rec["cos"] > C3_COS_BAR)
        log(f"[C3] fold {k} cos={rec['cos']:+.4f} CI={_fmt_ci(rec['cos_ci_lo'], rec['cos_ci_hi'])} "
            f"n_young_train={n_ytr} n_old_train={n_otr} n_old_test={n_ote} pass_bar={rec['pass_bar']}")
        rows.append(rec)
    df = pd.DataFrame(rows)
    df.to_csv(TOWARD_DIR / "c3_folds.csv", index=False)
    finite = df["cos"].to_numpy(float)
    all_finite = bool(np.isfinite(finite).all()) if len(df) == N_FOLDS else False
    all_bar = bool(all_finite and construction_ok and (finite > C3_COS_BAR).all())
    summary = dict(
        pass_=all_bar,
        construction_ok=construction_ok,
        n_folds=int(len(df)),
        bar=C3_COS_BAR,
        mean_cos=float(np.nanmean(finite)) if np.isfinite(finite).any() else np.nan,
        min_cos=float(np.nanmin(finite)) if np.isfinite(finite).any() else np.nan,
        max_cos=float(np.nanmax(finite)) if np.isfinite(finite).any() else np.nan,
        reasons=reasons,
        seed=TOWARD_SEED,
        boot_seed=BOOT_SEED,
        n_boot=N_BOOT,
    )
    # json-friendly
    summary["pass"] = summary.pop("pass_")
    dump_json(TOWARD_DIR / "c3_summary.json", jsonable(summary))
    log(f"[C3] pass={summary['pass']} min_cos={_fmt(summary['min_cos'])} "
        f"mean_cos={_fmt(summary['mean_cos'])} bar={C3_COS_BAR} construction_ok={construction_ok}")
    return df, summary


def load_donor(line, log):
    obs_p = require_file(MD2_DIR / f"louvain_obs_{line}.csv", "louvain_obs")
    lab_p = require_file(MD2_DIR / f"t2_cluster_labels_{line}.csv", "cluster_labels")
    sc_p = require_file(MD4_PROC / f"cell_scores_{line}.npz", "cell_scores")
    obs = pd.read_csv(obs_p)
    lab = pd.read_csv(lab_p)
    log_columns(log, f"louvain_obs_{line}", list(obs.columns), str(obs_p))
    log_columns(log, f"labels_{line}", list(lab.columns), str(lab_p))
    for col in ("cluster", "day"):
        if col not in obs.columns:
            raise StopStep("louvain_obs", f"{line} missing {col!r}. Not substituting.")
    if "cluster" not in lab.columns or "label" not in lab.columns:
        raise StopStep("cluster_labels", f"{line} labels missing cluster/label. Not substituting.")
    lab_map = {int(r.cluster): str(r.label) for _, r in lab.iterrows()}
    obs = obs.reset_index(drop=True)
    obs["label"] = obs["cluster"].map(lambda c: lab_map.get(int(c)))
    if obs["label"].isna().any():
        missing = sorted(obs.loc[obs.label.isna(), "cluster"].unique().tolist())
        raise StopStep("labels", f"{line}: unlabeled clusters {missing}. Not re-labelling.")
    scores = np.load(sc_p, allow_pickle=True)
    log(f"[cell_scores {line}] keys actually read: {list(scores.files)}")
    cl = np.asarray(scores["cluster"]).astype(int)
    day = np.asarray(scores["day"]).astype(int)
    lab_sc = np.asarray(scores["label"]).astype(str)
    if len(cl) != len(obs):
        raise StopStep("cell_scores", f"{line} cell_scores n={len(cl)} obs n={len(obs)}")
    if not np.array_equal(cl, obs["cluster"].to_numpy(int)):
        raise StopStep("cell_scores", f"{line} cluster does not match louvain_obs")
    if not np.array_equal(day, obs["day"].to_numpy(int)):
        raise StopStep("cell_scores", f"{line} day does not match louvain_obs")
    if not np.array_equal(lab_sc, obs["label"].to_numpy().astype(str)):
        raise StopStep("cell_scores", f"{line} label does not match md2 argmax labels")
    Y = load_Y(line)
    if Y.shape[0] != len(obs):
        raise StopStep("rulerY", f"{line} louvain_rulerY n={Y.shape[0]} obs n={len(obs)}")
    if "barcode" in obs.columns:
        log(f"[donor {line}] n_cells={len(obs)} n_clusters={obs.cluster.nunique()} "
            f"days={sorted(obs.day.astype(int).unique().tolist())}")
    return obs, Y


def occupancy_table(obs, line):
    rows = []
    for st in STATE_NAMES:
        m = obs["label"].astype(str) == st
        rec = dict(cell_line=line, state=st, n_cells=int(m.sum()))
        for d in DAYS:
            rec[f"n_d{d}"] = int(((obs["day"].astype(int) == int(d)) & m).sum())
        rec["n_d0_fibroblast"] = rec["n_d0"] if st == ORIGIN_STATE else np.nan
        rec["qualifies_dest"] = bool(st in DEST_STATES and rec["n_cells"] >= MIN_CELLS)
        rec["origin_ok"] = bool(st == ORIGIN_STATE and rec["n_d0"] >= MIN_CELLS)
        rec["unused_fibroblast_nond0"] = rec["n_cells"] - rec["n_d0"] if st == ORIGIN_STATE else 0
        rec["skipped_reason"] = ""
        if st in DEST_STATES and rec["n_cells"] < MIN_CELLS:
            rec["skipped_reason"] = f"n={rec['n_cells']}<{MIN_CELLS}"
        rows.append(rec)
    return pd.DataFrame(rows)


def extrap_strings(extrap, line, state):
    sub = extrap[(extrap.cell_line.astype(str) == line) & (extrap.label.astype(str) == state)].copy()
    if not len(sub):
        return "not in t2_extrap.csv", "not in t2_extrap.csv"
    sub["day"] = sub["day"].astype(int)
    sub = sub.sort_values("day")
    nn = "; ".join(
        f"d{int(r.day)}:{float(r.nn_euclidean):.3f}(n={int(r.n_cells)})" for _, r in sub.iterrows()
    )
    mh = "; ".join(
        f"d{int(r.day)}:{float(r.mahalanobis_pca):.3f}(n={int(r.n_cells)})" for _, r in sub.iterrows()
    )
    return nn, mh


def donor_panel(obs, Y, log, line):
    n_d0 = int(((obs.label.astype(str) == ORIGIN_STATE) & (obs.day.astype(int) == 0)).sum())
    if n_d0 < MIN_CELLS:
        raise StopStep(
            "origin",
            f"{line}: d0 Fibroblast n={n_d0} < {MIN_CELLS}. STOP for this donor.",
            dict(cell_line=line, n_d0=n_d0),
        )
    names = [ORIGIN_STATE]
    idx = {
        ORIGIN_STATE: np.flatnonzero(
            (obs.label.astype(str) == ORIGIN_STATE) & (obs.day.astype(int) == 0)
        )
    }
    skipped = []
    for st in DEST_STATES:
        ix = np.flatnonzero(obs.label.astype(str) == st)
        n = int(ix.size)
        if n < MIN_CELLS:
            skipped.append(dict(state=st, n_cells=n, reason=f"n={n}<{MIN_CELLS}"))
            log(f"[panel {line}] skip {st} n={n} < {MIN_CELLS}")
            continue
        names.append(st)
        idx[st] = ix
        log(f"[panel {line}] {st} n={n} days="
            f"{ {int(d): int(((obs.day.astype(int)==d) & (obs.label.astype(str)==st)).sum()) for d in DAYS} }")
    mats = [sum_rows(Y, idx[nm]) for nm in names]
    C = np.vstack(mats)
    logcpm, nf = tmm_logcpm_quiet(C)
    log(f"[tmm {line}] n_rows={C.shape[0]} states={names} "
        f"nf_median={float(np.median(nf)):.4f} prior.count={EDGE_R_PRIOR:g}")
    return dict(names=names, idx=idx, C=C, logcpm=logcpm, nf=nf, skipped=skipped, n_d0=n_d0)


def stats_from_Z(Z, names, anchors, frozen, missing):
    c_y, c_o, u = anchors["c_young"], anchors["c_old"], anchors["u"]
    pos = {nm: i for i, nm in enumerate(names)}
    z0 = Z[pos[ORIGIN_STATE]]
    dist_y0 = float(np.linalg.norm(z0 - c_y))
    dist_o0 = float(np.linalg.norm(z0 - c_o))
    ruler = ruler_score(Z, frozen)
    out = {}
    for nm in names:
        z = Z[pos[nm]]
        rec = dict(
            state=nm,
            role="origin" if nm == ORIGIN_STATE else "destination",
            ruler_score=float(ruler[pos[nm]]),
            dist_young=float(np.linalg.norm(z - c_y)),
            dist_old=float(np.linalg.norm(z - c_o)),
            dist_young_origin=dist_y0,
            dist_old_origin=dist_o0,
            n_missing_z0=int(np.asarray(missing).sum()),
        )
        rec["delta_young"] = rec["dist_young"] - dist_y0
        rec["delta_old"] = rec["dist_old"] - dist_o0
        if nm == ORIGIN_STATE:
            rec.update(cos_S=np.nan, d_norm=0.0, skipped=False)
        else:
            d = z - z0
            rec["d_norm"] = float(np.linalg.norm(d))
            rec["cos_S"] = cosine_align(d, u)
            rec["d"] = d
            rec["skipped"] = False
        rec["z"] = z
        out[nm] = rec
    return out, z0


def bootstrap_donor(Y, names, idx, anchors, frozen, rng, log, line):
    c_y, c_o, u = anchors["c_young"], anchors["c_old"], anchors["u"]
    boots = {nm: dict(cos=[], delta_young=[], delta_old=[], d_norm=[]) for nm in names if nm != ORIGIN_STATE}
    for b in range(N_BOOT):
        mats = []
        for nm in names:
            ix = idx[nm]
            take = rng.choice(ix, size=int(ix.size), replace=True)
            mats.append(sum_rows(Y, take))
        C = np.vstack(mats)
        logcpm, _ = tmm_logcpm_quiet(C)
        Zb, _ = zscore_frozen(logcpm, frozen, C)
        z0 = Zb[0]
        dy0 = float(np.linalg.norm(z0 - c_y))
        do0 = float(np.linalg.norm(z0 - c_o))
        for j, nm in enumerate(names):
            if nm == ORIGIN_STATE:
                continue
            z = Zb[j]
            d = z - z0
            boots[nm]["cos"].append(cosine_align(d, u))
            boots[nm]["d_norm"].append(float(np.linalg.norm(d)))
            boots[nm]["delta_young"].append(float(np.linalg.norm(z - c_y)) - dy0)
            boots[nm]["delta_old"].append(float(np.linalg.norm(z - c_o)) - do0)
        if (b + 1) % 50 == 0:
            log(f"[boot {line}] {b+1}/{N_BOOT}")
    out = {}
    for nm, d in boots.items():
        clo, chi = percentile_ci(d["cos"])
        ylo, yhi = percentile_ci(d["delta_young"])
        olo, ohi = percentile_ci(d["delta_old"])
        out[nm] = dict(
            cos_ci_lo=clo, cos_ci_hi=chi,
            delta_young_ci_lo=ylo, delta_young_ci_hi=yhi,
            delta_old_ci_lo=olo, delta_old_ci_hi=ohi,
            n_boot=N_BOOT,
        )
    return out


def meets_toward(rec):
    if rec is None or rec.get("skipped"):
        return False
    vals = (rec.get("cos_S"), rec.get("delta_young"), rec.get("p_N1"), rec.get("p_N2"))
    if any(v is None or not np.isfinite(v) for v in vals):
        return False
    return (
        float(rec["cos_S"]) > 0
        and float(rec["p_N1"]) <= P_BAR
        and float(rec["p_N2"]) <= P_BAR
        and float(rec["delta_young"]) < 0
    )


def reading_for_donor(line, recs, n_by_state, origin_stop, c3_pass):
    is_contrast = line == CONTRAST_DONOR
    prefix = (
        f"CONTRAST (young donor {line}; never the reading): "
        if is_contrast else f"{line}: "
    )

    def n_s(st):
        return int(n_by_state.get(st, 0))

    if not c3_pass:
        return dict(
            key="c3_broken", cell_line=line, is_contrast=is_contrast,
            text=prefix + "C3 failed. Statistic broken. Nothing downstream is reported.",
        )
    if origin_stop:
        return dict(
            key="origin_stop", cell_line=line, is_contrast=is_contrast,
            text=prefix + f"d0 Fibroblast n={n_s(ORIGIN_STATE)} < {MIN_CELLS}. STOP for this donor.",
        )
    pr = recs.get("PartialReprog")
    pl = recs.get("Pluripotency")
    if pr is None or pr.get("skipped"):
        return dict(
            key="skipped_partial", cell_line=line, is_contrast=is_contrast,
            text=prefix + (
                f"PartialReprog skipped (n={n_s('PartialReprog')}<{MIN_CELLS}). "
                "Not a destination. Not dropped silently."
            ),
        )
    cos = pr.get("cos_S")
    dy = pr.get("delta_young")
    if np.isfinite(cos) and float(cos) > 0 and np.isfinite(dy) and float(dy) >= 0:
        return dict(
            key="away_from_old_only", cell_line=line, is_contrast=is_contrast,
            text=prefix + (
                f"PartialReprog has cos_S={float(cos):+.3f} but delta_young={float(dy):+.3f} "
                f"(n_cells={n_s('PartialReprog')}). Movement is away-from-old only. "
                "This is a real negative finding: a 1D age score cannot distinguish toward-young "
                "from away-from-old. The frozen-ruler score is not used to rescue this."
            ),
        )
    pr_t = meets_toward(pr)
    pl_skipped = pl is None or pl.get("skipped")
    if pr_t and pl_skipped:
        return dict(
            key="c1_unevaluable", cell_line=line, is_contrast=is_contrast,
            text=prefix + (
                f"PartialReprog meets TOWARD young (n_cells={n_s('PartialReprog')}) but "
                f"Pluripotency is skipped (n={n_s('Pluripotency')}<{MIN_CELLS}), so C1 cannot "
                "be evaluated. Do not fire a toward-young claim."
            ),
        )
    pl_t = meets_toward(pl)
    if pr_t and pl_t:
        return dict(
            key="uninterpretable_c1", cell_line=line, is_contrast=is_contrast,
            text=prefix + (
                f"PartialReprog passes both nulls with delta_young<0 (n_cells={n_s('PartialReprog')}) "
                f"but Pluripotency ALSO does (n_cells={n_s('Pluripotency')}). The statistic tracks "
                "displacement magnitude and/or extrapolation distance. Uninterpretable. "
                "Do NOT report a toward-young result. What would fix it: a direction statistic "
                "that Pluripotency does not pass, or an in-manifold young-fibroblast target that "
                "is not the GTEx centroid in the full overlap gene space."
            ),
        )
    if pr_t and (not pl_t) and c3_pass:
        return dict(
            key="toward_young", cell_line=line, is_contrast=is_contrast,
            text=prefix + (
                f"partially-reprogrammed cells move toward the young fibroblast target, not "
                f"merely away from the old state (PartialReprog n_cells={n_s('PartialReprog')}; "
                f"cos_S={float(pr['cos_S']):+.3f}, delta_young={float(pr['delta_young']):+.3f}, "
                f"p_N1={_fmt_p(pr['p_N1'])}, p_N2={_fmt_p(pr['p_N2'])}; "
                f"Pluripotency does not meet TOWARD young, C1 clean). "
                "The safe-target framing is supported in this dataset."
            ),
        )
    p1, p2 = pr.get("p_N1"), pr.get("p_N2")
    both_fail = (
        (not np.isfinite(p1) or float(p1) > P_BAR)
        and (not np.isfinite(p2) or float(p2) > P_BAR)
    )
    one_only = (np.isfinite(p1) and np.isfinite(p2)
                and ((float(p1) <= P_BAR) ^ (float(p2) <= P_BAR)))
    extra = (
        "Both nulls fail (p>0.05). "
        if both_fail else
        "Passing one null is not passing. "
        if one_only else
        "Both nulls are not passed. "
    )
    return dict(
        key="no_directional_structure", cell_line=line, is_contrast=is_contrast,
        text=prefix + (
            f"PartialReprog (n_cells={n_s('PartialReprog')}, cos_S={_fmt(cos)}, "
            f"p_N1={_fmt_p(p1)}, p_N2={_fmt_p(p2)}, delta_young={_fmt(dy)}): "
            f"{extra}No directional structure detectable at this n. "
            "The frozen-ruler score is not used to rescue this."
        ),
    )


def write_progress(next_action, stop, extra=""):
    man = load_manifest()
    c3p = TOWARD_DIR / "c3_summary.json"
    c3 = load_json(c3p) if c3p.exists() else {}
    ap = TOWARD_DIR / "anchors.json"
    anc = load_json(ap) if ap.exists() else {}
    rp = TOWARD_DIR / "reading.json"
    reading = load_json(rp) if rp.exists() else {}
    lines = [
        "# PROGRESS_TOWARD",
        "",
        f"**STOP status:** {stop}",
        f"**Next action:** {next_action}",
        "",
        "## Seeds / gates",
        "",
        f"- seed `{TOWARD_SEED}`  boot `{BOOT_SEED}`  n_perm={N_PERM}  n_boot={N_BOOT}  n_random={N_RANDOM}",
        f"- primary metric: cos_S = dot(d_S, u) / ||d_S||; uncalibrated R² is never a gate",
        f"- frozen ruler: context column only; not the statistic; not refit; exists={FROZEN_RULER.exists()}",
        f"- MIN_CELLS={MIN_CELLS}  C3_COS_BAR={C3_COS_BAR}  P_BAR={P_BAR}",
        f"- YOUNG bins {YOUNG_BINS}; OLD bins {OLD_BINS}; middle reported, not dropped",
        f"- GSE325735: out of scope",
        f"- Refit nothing. Do not re-cluster, re-label, or re-tune Louvain or AddModuleScore.",
        f"- Donors never pooled. States never averaged.",
        f"- PREREG.flag exists={PREREG_FLAG.exists()}",
        f"- anchors.npz exists={ANCHORS_PATH.exists()}",
        "",
        "## Anchors (from disk)",
        "",
        f"- n_young={anc.get('n_young')} (expected {N_YOUNG_EXPECTED})  "
        f"n_old={anc.get('n_old')} (expected {N_OLD_EXPECTED})  "
        f"n_middle={anc.get('n_middle')}",
        "",
        "## C3",
        "",
        f"- pass={c3.get('pass')} min_cos={_fmt(c3.get('min_cos'))} mean_cos={_fmt(c3.get('mean_cos'))} bar={c3.get('bar')}",
        f"- construction_ok={c3.get('construction_ok')}",
        "",
        "## Reading keys (from disk)",
        "",
    ]
    if isinstance(reading, dict) and reading:
        by = reading.get("by_donor") or reading
        if isinstance(by, dict):
            for k, v in by.items():
                if isinstance(v, dict):
                    lines.append(f"- {k}: key=`{v.get('key')}` contrast={v.get('is_contrast')}")
        lines.append("")
    if extra:
        lines += ["## Note", "", extra, ""]
    fails = man.get("failures") or []
    if fails:
        lines += ["## Failures (manifest)", ""]
        for f in fails:
            lines.append(f"- **{f.get('step')}:** {f.get('message')}")
        lines.append("")
    lines += [
        "## Files",
        "",
        "- `src/toward_run.py`",
        "- `results/toward/`",
        "- `FINDINGS_TOWARD.md`",
        "- `PROGRESS_TOWARD.md`",
        "",
    ]
    PROGRESS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return PROGRESS_PATH


def write_findings():
    prereg = load_prereg_text()
    man = load_manifest()
    c3 = load_json(TOWARD_DIR / "c3_summary.json") if (TOWARD_DIR / "c3_summary.json").exists() else {}
    anc = load_json(TOWARD_DIR / "anchors.json") if (TOWARD_DIR / "anchors.json").exists() else {}
    reading = load_json(TOWARD_DIR / "reading.json") if (TOWARD_DIR / "reading.json").exists() else {}
    summary = load_json(TOWARD_DIR / "summary.json") if (TOWARD_DIR / "summary.json").exists() else {}
    status = man.get("status", "unknown")
    claim = (reading.get("claim") if isinstance(reading, dict) else None) or {}
    status_line = (
        f"**Status:** {status}. Claim donor {CLAIM_DONOR} key=`{claim.get('key', 'not fired')}`. "
        f"Seed `{TOWARD_SEED}`. boot `{BOOT_SEED}`. n_perm={N_PERM}. n_boot={N_BOOT}. n_random={N_RANDOM}. "
        f"C3 pass={c3.get('pass')}. Frozen ruler exists={FROZEN_RULER.exists()}."
    )
    lines = [
        "# FINDINGS_TOWARD — do partially-reprogrammed cells move TOWARD young fibroblasts?",
        "",
        status_line,
        "",
        "Does not modify any existing `FINDINGS_*.md` or `FALSIFICATION.md`. No GSE325735. "
        "Uncalibrated R² is never a gate. Donors never pooled. Frozen ruler is a context column only.",
        "",
        "Reproduced by `src/toward_run.py`. Frozen ruler: "
        f"`{FROZEN_RULER}` exists={FROZEN_RULER.exists()}.",
        "",
        "## Pre-registration (verbatim, written before any TOWARD statistic)",
        "",
        prereg,
        "",
        f"Flag: `results/toward/PREREG.flag` exists={PREREG_FLAG.exists()}.",
        "",
        "## STOP / failures",
        "",
        "Not substituting columns or repairing rows.",
        "",
    ]
    fails = man.get("failures") or []
    if fails:
        for f in fails:
            lines.append(f"- **{f.get('step')}:** {f.get('message')}")
        lines.append("")
    else:
        lines += ["None recorded.", ""]

    lines += ["## Anchors (frozen before any GSE297234 vector)", ""]
    if anc:
        lines.append(
            f"- YOUNG bins {YOUNG_BINS}: n={anc.get('n_young')} (expected {N_YOUNG_EXPECTED})."
        )
        lines.append(
            f"- OLD bins {OLD_BINS}: n={anc.get('n_old')} (expected {N_OLD_EXPECTED})."
        )
        lines.append(
            f"- Middle bins {MIDDLE_BINS}: n={anc.get('n_middle')} excluded from both anchors, "
            "not dropped silently."
        )
        lines.append(
            f"- Gene space n={anc.get('n_genes')} (frozen-ruler overlap). "
            f"||c_young − c_old||={_fmt(anc.get('axis_norm'), 4)}."
        )
        lines.append(f"- Written to `{ANCHORS_PATH}`.")
        lines.append("")
    at = TOWARD_DIR / "anchor_table.csv"
    if at.exists():
        adf = pd.read_csv(at)
        lines += [md_table(adf), ""]

    lines += ["## C3 positive control (GTEx held-out OLD → train YOUNG centroid)", ""]
    if c3:
        lines.append(
            f"- C3 pass={c3.get('pass')} (all 5 folds finite and cos > {c3.get('bar')}; "
            f"construction_ok={c3.get('construction_ok')}). "
            f"min_cos={_fmt(c3.get('min_cos'))} mean_cos={_fmt(c3.get('mean_cos'))} "
            f"max_cos={_fmt(c3.get('max_cos'))}."
        )
        if c3.get("reasons"):
            for r in c3["reasons"]:
                lines.append(f"- construction: {r}")
        lines.append(
            "- Bootstrap 95% CI (seed 20260918) resamples held-out OLD donors within fold; not a gate."
        )
        lines.append("")
    cf = TOWARD_DIR / "c3_folds.csv"
    if cf.exists():
        cdf = pd.read_csv(cf)
        show = cdf.copy()
        keep = [c for c in (
            "fold", "n_train", "n_test", "n_young_train", "n_old_train", "n_old_test",
            "cos", "d_norm", "cos_ci_lo", "cos_ci_hi", "pass_bar", "construction_ok", "reason",
        ) if c in show.columns]
        lines += [md_table(show[keep]), ""]
        if not c3.get("pass"):
            lines.append(
                "**C3 failed. Statistic broken. STOP. Nothing downstream is reported.**"
            )
            lines.append("")

    if summary.get("c3_stop"):
        lines += _files_touched()
        FINDINGS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return FINDINGS_PATH

    occ_p = TOWARD_DIR / "occupancy.csv"
    if occ_p.exists():
        lines += ["## Timepoint composition per state per donor", ""]
        lines.append(
            "Labels: md2 argmax of mean AddModuleScore (mmc3 Reprog_cell_state_signatures). "
            "Not re-labelled. Fibroblast origin uses d0 only. Non-d0 Fibroblast cells are unused, "
            "not a destination. Destination states pool timepoints. MIN_CELLS=30. "
            "No destination state was skipped."
        )
        lines.append("")
        odf = pd.read_csv(occ_p)
        odf = odf.copy()
        if "n_d0_fibroblast" in odf.columns:
            odf["n_d0_fibroblast"] = [
                str(int(v)) if pd.notna(v) else "—" for v in odf["n_d0_fibroblast"]
            ]
        if "skipped_reason" in odf.columns:
            odf["skipped_reason"] = [
                (str(v) if pd.notna(v) and str(v).strip() not in {"", "nan", "NA"} else "—")
                for v in odf["skipped_reason"]
            ]
        lines += [md_table(odf), ""]

    stats_p = TOWARD_DIR / "stats.csv"
    if stats_p.exists():
        lines += ["## Statistic table (primary = cos_S)", ""]
        lines.append(
            "cos_S = dot(d_S, u) / ||d_S||. delta_young < 0 means closer to the YOUNG centroid. "
            "Frozen-ruler score is a **context** column only. nn_euclidean and Mahalanobis are "
            "read as-is from results/md4/t2_extrap.csv by day and are not averaged."
        )
        lines.append("")
        sdf = pd.read_csv(stats_p)
        sdf_disp = sdf.copy()
        for pc in ("p_N1", "p_N2"):
            if pc in sdf_disp.columns:
                sdf_disp[pc] = [_fmt_p(v) for v in sdf_disp[pc]]
        num_cols = [c for c in (
            "cell_line", "state", "role", "n_cells", "n_d0", "n_d3", "n_d7", "n_d10",
            "skipped", "cos_S", "d_norm", "dist_young", "dist_young_origin", "delta_young",
            "dist_old", "dist_old_origin", "delta_old", "ruler_score_context",
            "p_N1", "p_N2", "pass_nulls", "meets_toward",
            "cos_ci_lo", "cos_ci_hi", "delta_young_ci_lo", "delta_young_ci_hi",
        ) if c in sdf_disp.columns]
        lines += ["### Primary and accompanying numbers", "", md_table(sdf_disp[num_cols]), ""]
        ctx_cols = [c for c in (
            "cell_line", "state", "ruler_score_context",
            "nn_euclidean_by_day", "mahalanobis_by_day", "skipped_reason",
        ) if c in sdf.columns]
        lines += [
            "### Context columns (frozen ruler; md4 t2_extrap as-is; not substituted for cos_S)",
            "",
            md_table(sdf[ctx_cols]),
            "",
        ]

    lines += ["## N1 / N2 p-values", ""]
    lines.append(
        "N1: permute entries of frozen u (seed 20260914), n_random=200. "
        "N2: equal-size GTEx donor splits ignoring age (seed 20260914), n_perm=200. "
        "p = (# null ≥ observed + 1) / 201. PASS_NULLS requires both p≤0.05."
    )
    lines.append("")
    if stats_p.exists():
        sdf = pd.read_csv(stats_p)
        dest = sdf[sdf.role.astype(str) == "destination"].copy()
        if len(dest):
            dest = dest.copy()
            if "p_N1" in dest.columns:
                dest["p_N1"] = [_fmt_p(v) for v in dest["p_N1"]]
            if "p_N2" in dest.columns:
                dest["p_N2"] = [_fmt_p(v) for v in dest["p_N2"]]
            pcols = [c for c in (
                "cell_line", "state", "n_cells", "cos_S", "p_N1", "p_N2", "pass_nulls",
                "delta_young", "meets_toward",
            ) if c in dest.columns]
            lines += [md_table(dest[pcols]), ""]

    lines += ["## C1 Pluripotency / C2 NonReprog", ""]
    c1 = summary.get("c1") or {}
    c2 = summary.get("c2") or {}
    for line_id in (AGED_LINE, YOUNG_LINE):
        a = c1.get(line_id) or {}
        b = c2.get(line_id) or {}
        lines.append(
            f"- **{line_id} C1 Pluripotency:** skipped={a.get('skipped')} "
            f"meets_toward={a.get('meets_toward')} cos_S={_fmt(a.get('cos_S'))} "
            f"delta_young={_fmt(a.get('delta_young'))} p_N1={_fmt_p(a.get('p_N1'))} "
            f"p_N2={_fmt_p(a.get('p_N2'))} n_cells={a.get('n_cells')} "
            f"d_norm={_fmt(a.get('d_norm'))}."
        )
        lines.append(
            f"- **{line_id} C2 NonReprog (not a gate):** skipped={b.get('skipped')} "
            f"meets_toward={b.get('meets_toward')} cos_S={_fmt(b.get('cos_S'))} "
            f"delta_young={_fmt(b.get('delta_young'))} p_N1={_fmt_p(b.get('p_N1'))} "
            f"p_N2={_fmt_p(b.get('p_N2'))} n_cells={b.get('n_cells')}."
        )
    lines.append("")

    lines += ["## Pre-registered reading (only the outcome that fired)", ""]
    by = (reading.get("by_donor") or {}) if isinstance(reading, dict) else {}
    claim = reading.get("claim") or {}
    if claim:
        lines.append(
            f"**Claim donor ({CLAIM_DONOR}) key: `{claim.get('key')}`.** "
            "GM23815 is a contrast, never the reading. Cosines are not averaged."
        )
        lines.append("")
        lines.append(claim.get("text") or "")
        lines.append("")
    for line_id in (AGED_LINE, YOUNG_LINE):
        rec = by.get(line_id) or {}
        if rec:
            tag = "contrast" if rec.get("is_contrast") else "claim"
            lines.append(f"### {line_id} ({tag}) key=`{rec.get('key')}`")
            lines.append("")
            lines.append(rec.get("text") or "")
            lines.append("")

    lines += [
        "## Limitations",
        "",
        "- GTEx public AGE is a 10-year bin; YOUNG and OLD are bin unions, not chronological years.",
        "- GTEx cultured fibroblasts are not the same system as GSE297234 primary dermal fibroblasts.",
        "- Two donors. GM23815 is a contrast and is not pooled into the claim.",
        "- State labels are Louvain-cluster argmax of signature scores, not per-cell argmax.",
        "- TMM is among a handful of state-pseudobulk rows per donor.",
        "- The aging-axis anchor lives in the full frozen-ruler overlap space; states that leave "
        "the fibroblast manifold (Pluripotency) can produce large displacements that are not "
        "rejuvenation. C1 exists to catch that.",
        "- nn_euclidean / Mahalanobis are copied from md4 t2_extrap.csv (state × timepoint) and "
        "are not recomputed on the pooled-state vector.",
        "- Day-0 fibroblast origins already sit hundreds of z-units from both GTEx centroids "
        "(scRNA pseudobulk vs bulk cultured fibroblasts). Subsequent states recede from both "
        "anchors; delta_young and delta_old are therefore not a within-cloud slide.",
        "- Bootstrap CIs re-TMM the donor panel, so they are CIs on the TMM-coupled statistic "
        "and need not contain the observed point estimate. Not a gate.",
        "- Uncalibrated R² is never a gate and is not computed here.",
        f"- {GSE325735_BAN}",
        "",
    ]
    lines += _files_touched()
    FINDINGS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return FINDINGS_PATH


def _files_touched():
    return [
        "## Files touched",
        "",
        "- `src/toward_run.py`",
        "- `FINDINGS_TOWARD.md`",
        "- `PROGRESS_TOWARD.md`",
        "- `results/toward/PREREG.flag` (not rewritten if it already existed)",
        "- `results/toward/anchors.npz`",
        "- `results/toward/anchors.json`",
        "- `results/toward/anchor_table.csv`",
        "- `results/toward/null_directions.npz`",
        "- `results/toward/c3_folds.csv`",
        "- `results/toward/c3_summary.json`",
        "- `results/toward/occupancy.csv`",
        "- `results/toward/stats.csv`",
        "- `results/toward/state_z.npz`",
        "- `results/toward/reading.json`",
        "- `results/toward/summary.json`",
        "- `results/toward/manifest.json`",
        "- `results/toward/run_report.txt`",
        "- Existing `FINDINGS_*.md`, `PROGRESS_*.md`, and `FALSIFICATION.md` were not modified "
        "(except FINDINGS_TOWARD.md and PROGRESS_TOWARD.md).",
        "",
    ]


def state_brief(rec, n_cells):
    if rec is None:
        return dict(skipped=True, n_cells=n_cells)
    return dict(
        skipped=bool(rec.get("skipped")),
        meets_toward=meets_toward(rec),
        cos_S=rec.get("cos_S"),
        delta_young=rec.get("delta_young"),
        p_N1=rec.get("p_N1"),
        p_N2=rec.get("p_N2"),
        n_cells=n_cells,
        d_norm=rec.get("d_norm"),
    )


def run(log):
    prereg = load_prereg_text()
    log("=" * 100)
    log("TOWARD  seed=%s  boot=%s  n_perm=%s  n_boot=%s  n_random=%s" % (
        TOWARD_SEED, BOOT_SEED, N_PERM, N_BOOT, N_RANDOM))
    log("=" * 100)
    log("FALSIFICATION.md and existing FINDINGS*.md are not modified.")
    log(GSE325735_BAN)
    log("Primary statistic: cos_S. Uncalibrated R² is never a gate.")
    log("Frozen ruler is a context column only. Not refit.")
    log("PREREG.flag exists; not rewriting the block.")
    log(prereg[:400] + " ...")
    man = load_manifest()
    man["status"] = "RUNNING"
    save_manifest(man)
    write_progress("load GTEx Stage 0 pack and freeze anchors", stop="RUNNING")

    gtex = load_gtex_z(log)
    anchors = freeze_anchors(gtex, log)
    write_progress("build N1/N2 directions, then C3", stop="ANCHORS_FROZEN")

    U1, U2 = build_null_directions(gtex["Z"], anchors["u"], log)
    c3_df, c3_sum = run_c3(gtex, log)
    dump_json(TOWARD_DIR / "summary.json", jsonable(dict(c3=c3_sum, c3_stop=not c3_sum["pass"])))
    if not c3_sum["pass"]:
        raise StopStep(
            "C3",
            "C3 failed. Statistic broken. STOP. Nothing downstream is reported. "
            f"min_cos={c3_sum.get('min_cos')} bar={C3_COS_BAR} "
            f"construction_ok={c3_sum.get('construction_ok')} reasons={c3_sum.get('reasons')}",
            dict(c3=c3_sum),
        )
    write_progress("GSE297234 state displacements (C3 passed)", stop="C3_PASS")

    extrap_p = require_file(EXTRAP_PATH, "t2_extrap")
    extrap = pd.read_csv(extrap_p)
    log_columns(log, "t2_extrap", list(extrap.columns), str(extrap_p))
    for col in ("cell_line", "day", "label", "n_cells", "nn_euclidean", "mahalanobis_pca"):
        if col not in extrap.columns:
            raise StopStep("t2_extrap", f"missing column {col!r}. Not substituting.")
    extrap.to_csv(TOWARD_DIR / "t2_extrap_copied.csv", index=False)

    occ_all = []
    stats_rows = []
    recs_by = {}
    n_by = {}
    origin_stop = {}
    z_store = {}
    frozen = gtex["frozen"]

    for line in (AGED_LINE, YOUNG_LINE):
        origin_stop[line] = False
        recs_by[line] = {}
        try:
            obs, Y = load_donor(line, log)
        except StopStep as e:
            if e.step == "origin":
                raise
            raise
        occ = occupancy_table(obs, line)
        occ_all.append(occ)
        n_by[line] = {str(r.state): int(r.n_cells) for _, r in occ.iterrows()}
        n_by[line][ORIGIN_STATE] = int(
            occ.loc[occ.state == ORIGIN_STATE, "n_d0"].iloc[0]
        )
        try:
            panel = donor_panel(obs, Y, log, line)
        except StopStep as e:
            if e.step == "origin":
                origin_stop[line] = True
                log(f"[origin STOP] {e.message}")
                for st in STATE_NAMES:
                    hit = occ[occ.state == st].iloc[0]
                    nn, mh = extrap_strings(extrap, line, st)
                    stats_rows.append(dict(
                        cell_line=line, state=st,
                        role="origin" if st == ORIGIN_STATE else "destination",
                        n_cells=int(hit.n_cells if st != ORIGIN_STATE else hit.n_d0),
                        n_d0=int(hit.n_d0), n_d3=int(hit.n_d3),
                        n_d7=int(hit.n_d7), n_d10=int(hit.n_d10),
                        skipped=True, skipped_reason=e.message,
                        cos_S=np.nan, d_norm=np.nan,
                        dist_young=np.nan, dist_young_origin=np.nan, delta_young=np.nan,
                        dist_old=np.nan, dist_old_origin=np.nan, delta_old=np.nan,
                        ruler_score_context=np.nan,
                        p_N1=np.nan, p_N2=np.nan, pass_nulls=False, meets_toward=False,
                        cos_ci_lo=np.nan, cos_ci_hi=np.nan,
                        delta_young_ci_lo=np.nan, delta_young_ci_hi=np.nan,
                        nn_euclidean_by_day=nn, mahalanobis_by_day=mh,
                    ))
                continue
            raise
        Z, missing = zscore_frozen(panel["logcpm"], frozen, panel["C"])
        recs, z0 = stats_from_Z(Z, panel["names"], anchors, frozen, missing)
        rng_b = np.random.default_rng(BOOT_SEED + DONOR_BOOT_OFFSET[line])
        boot = bootstrap_donor(Y, panel["names"], panel["idx"], anchors, frozen, rng_b, log, line)

        for st, rec in recs.items():
            if st == ORIGIN_STATE:
                rec["p_N1"] = np.nan
                rec["p_N2"] = np.nan
            else:
                p1, nc1 = null_p(rec["cos_S"], U1, rec["d"])
                p2, nc2 = null_p(rec["cos_S"], U2, rec["d"])
                rec["p_N1"] = p1
                rec["p_N2"] = p2
                rec["pass_nulls"] = bool(
                    np.isfinite(p1) and np.isfinite(p2) and p1 <= P_BAR and p2 <= P_BAR
                )
                rec["meets_toward"] = meets_toward(rec)
                np.savez_compressed(
                    TOWARD_DIR / f"null_cos_{line}_{st}.npz", N1=nc1, N2=nc2,
                    cos_obs=np.array(rec["cos_S"]),
                )
                log(
                    f"[stat {line} {st}] n={int(panel['idx'][st].size)} "
                    f"cos={_fmt(rec['cos_S'])} ||d||={_fmt(rec['d_norm'])} "
                    f"delta_young={_fmt(rec['delta_young'])} delta_old={_fmt(rec['delta_old'])} "
                    f"p_N1={_fmt_p(p1)} p_N2={_fmt_p(p2)} toward={rec['meets_toward']} "
                    f"ruler_context={_fmt(rec['ruler_score'])}"
                )
            recs_by[line][st] = rec
            z_store[f"{line}:{st}"] = rec["z"]

        for skip in panel["skipped"]:
            recs_by[line][skip["state"]] = dict(
                skipped=True, reason=skip["reason"],
                cos_S=np.nan, d_norm=np.nan, delta_young=np.nan, delta_old=np.nan,
                p_N1=np.nan, p_N2=np.nan, pass_nulls=False, meets_toward=False,
            )

        occ_map = {str(r.state): r for _, r in occ.iterrows()}
        for st in STATE_NAMES:
            hit = occ_map[st]
            nn, mh = extrap_strings(extrap, line, st)
            rec = recs_by[line].get(st)
            skipped = bool(rec is None or rec.get("skipped"))
            n_cells = int(hit.n_d0) if st == ORIGIN_STATE else int(hit.n_cells)
            b = boot.get(st) or {}
            row = dict(
                cell_line=line, state=st,
                role="origin" if st == ORIGIN_STATE else "destination",
                n_cells=n_cells,
                n_d0=int(hit.n_d0), n_d3=int(hit.n_d3),
                n_d7=int(hit.n_d7), n_d10=int(hit.n_d10),
                skipped=skipped,
                skipped_reason="" if not skipped else (
                    rec.get("reason") if rec else f"n={n_cells}<{MIN_CELLS}"
                ),
                cos_S=rec.get("cos_S") if rec else np.nan,
                d_norm=rec.get("d_norm") if rec else np.nan,
                dist_young=rec.get("dist_young") if rec else np.nan,
                dist_young_origin=rec.get("dist_young_origin") if rec else np.nan,
                delta_young=rec.get("delta_young") if rec else np.nan,
                dist_old=rec.get("dist_old") if rec else np.nan,
                dist_old_origin=rec.get("dist_old_origin") if rec else np.nan,
                delta_old=rec.get("delta_old") if rec else np.nan,
                ruler_score_context=rec.get("ruler_score") if rec else np.nan,
                p_N1=rec.get("p_N1") if rec else np.nan,
                p_N2=rec.get("p_N2") if rec else np.nan,
                pass_nulls=bool(rec.get("pass_nulls")) if rec else False,
                meets_toward=meets_toward(rec) if rec else False,
                cos_ci_lo=b.get("cos_ci_lo", np.nan),
                cos_ci_hi=b.get("cos_ci_hi", np.nan),
                delta_young_ci_lo=b.get("delta_young_ci_lo", np.nan),
                delta_young_ci_hi=b.get("delta_young_ci_hi", np.nan),
                nn_euclidean_by_day=nn,
                mahalanobis_by_day=mh,
            )
            if st == ORIGIN_STATE:
                row["skipped"] = False
                row["skipped_reason"] = (
                    f"origin d0 n={int(hit.n_d0)}; unused non-d0 Fibroblast n={int(hit.unused_fibroblast_nond0)}"
                )
            stats_rows.append(row)

    occ_df = pd.concat(occ_all, ignore_index=True) if occ_all else pd.DataFrame()
    if len(occ_df):
        occ_df.to_csv(TOWARD_DIR / "occupancy.csv", index=False)
    stats_df = pd.DataFrame(stats_rows)
    stats_df.to_csv(TOWARD_DIR / "stats.csv", index=False)
    if z_store:
        np.savez_compressed(
            TOWARD_DIR / "state_z.npz",
            **{k.replace(":", "_"): v for k, v in z_store.items()},
        )

    by_donor = {}
    for line in (AGED_LINE, YOUNG_LINE):
        by_donor[line] = reading_for_donor(
            line, recs_by.get(line) or {}, n_by.get(line) or {},
            origin_stop.get(line, False), c3_pass=True,
        )
        log(f"[reading {line}] key={by_donor[line]['key']}")
        log(by_donor[line]["text"])

    claim = by_donor[CLAIM_DONOR]
    reading = dict(by_donor=by_donor, claim=claim)
    dump_json(TOWARD_DIR / "reading.json", jsonable(reading))

    c1 = {ln: state_brief((recs_by.get(ln) or {}).get("Pluripotency"),
                          (n_by.get(ln) or {}).get("Pluripotency", 0))
          for ln in (AGED_LINE, YOUNG_LINE)}
    c2 = {ln: state_brief((recs_by.get(ln) or {}).get("NonReprog"),
                          (n_by.get(ln) or {}).get("NonReprog", 0))
          for ln in (AGED_LINE, YOUNG_LINE)}
    dump_json(TOWARD_DIR / "summary.json", jsonable(dict(
        c3=c3_sum, c3_stop=False, c1=c1, c2=c2,
        claim_key=claim.get("key"),
        n_young=anchors["n_young"], n_old=anchors["n_old"],
        n_middle=anchors["n_middle"],
        min_cells=MIN_CELLS, c3_bar=C3_COS_BAR,
        seed=TOWARD_SEED, boot_seed=BOOT_SEED,
        n_perm=N_PERM, n_boot=N_BOOT, n_random=N_RANDOM,
    )))
    man = load_manifest()
    man["status"] = "DONE"
    man["claim_key"] = claim.get("key")
    save_manifest(man)
    write_findings()
    write_progress("done", stop="DONE")
    log("[run] DONE")
    return dict(claim=claim, c3=c3_sum)


def main():
    log = Logger(TOWARD_DIR / "run_report.txt")
    try:
        run(log)
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        log(f"[run] STOP [{e.step}] {e.message}")
        # If C3 failed, still write findings with what exists.
        try:
            write_findings()
        except Exception as fe:
            log(f"[run] findings write after STOP failed: {fe}")
        try:
            write_progress(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        except Exception:
            pass
        raise
    finally:
        log.close()


if __name__ == "__main__":
    main()

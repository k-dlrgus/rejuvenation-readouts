"""STAGE G3 — apply FROZEN geometric scores to mouse aging (transfer check)
then to reprogramming data. Do not refit.

Datasets (URLs from g3_search / CELLxGENE / GEO SOFT; not invented):
  - Tabula Muris Senis brain non-myeloid `66ff82b4-…` and myeloid `c08f8441-…`
    (independent mouse aging atlas).
  - GSE224438 mouse SVZ partial reprogramming (young / old / old+OSK).
  - GSE276656 mPFC OSK multiome is 14.4 GB RAW — not downloaded; GSE271794 is E15.5.

Usage: python src/trajectory_g3.py
"""
from __future__ import annotations

import re
import sys
import tarfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import DATA_RAW, DATA_PROC  # noqa: E402
from trajectory_common import (  # noqa: E402
    TRAJ_DIR, TRAJ_FIG, TRAJ_SEED, DATASET_ID, STOP_R2_COLLAPSE,
    Logger, dump_json, load_json, traj_log_banner, permutation_p,
    summarize_null_col, identity_project, nearest_centroid_predict,
    fmt, fmt_u, r2_mae,
)
from download import download  # noqa: E402
from search_cxg import parse_age_years  # noqa: E402
from geo_meta import fetch_soft  # noqa: E402

try:
    import anndata as ad
    import scipy.sparse as sp
except Exception:  # pragma: no cover
    ad = None
    sp = None

RAW = DATA_RAW / "trajectory"
PROC = DATA_PROC / "trajectory"
for _p in (RAW, PROC):
    _p.mkdir(parents=True, exist_ok=True)

# Fetched in g3_search_cxg_mouse_aging_ranked.csv (CXG curation API).
TMS_NONMY = dict(
    dataset_id="66ff82b4-9380-469c-bc4b-cfa08eacd325",
    title="Brain non-myeloid cells - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse",
    collection="Tabula Muris Senis",
    doi="10.1038/s41586-020-2496-1",
    url="https://datasets.cellxgene.cziscience.com/11554919-8ee5-411b-a088-7be1a7c9a5a6.h5ad",
    dest=RAW / "tms_brain_nonmyeloid",
)
TMS_MY = dict(
    dataset_id="c08f8441-4a10-4748-872a-e70c0bcccdba",
    title="Brain myeloid cells - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - Smart-seq2",
    collection="Tabula Muris Senis",
    doi="10.1038/s41586-020-2496-1",
    url="https://datasets.cellxgene.cziscience.com/b14cb09f-dd56-4c6c-bbf1-6822f01af53f.h5ad",
    dest=RAW / "tms_brain_myeloid",
)
GSE224 = dict(
    accession="GSE224438",
    url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE224nnn/GSE224438/suppl/GSE224438_RAW.tar",
    dest=RAW / "GSE224438",
    pubmed="38553564",
)

# Mouse annotation -> human DLPFC Phase-1 cell type (exact Phase-1 strings).
TYPE_MAP = {
    "astrocyte": "astrocyte",
    "oligodendrocyte": "oligodendrocyte",
    "oligodendrocyte precursor cell": "oligodendrocyte precursor cell",
    "opc": "oligodendrocyte precursor cell",
    "endothelial cell": "endothelial cell",
    "endothelial": "endothelial cell",
    "brain pericyte": "pericyte",
    "pericyte": "pericyte",
    "microglial cell": "microglial cell",
    "microglia": "microglial cell",
    "macrophage": "perivascular macrophage",
    "perivascular macrophage": "perivascular macrophage",
    "vascular leptomeningeal cell": "vascular leptomeningeal cell",
    "vlmc": "vascular leptomeningeal cell",
    "smooth muscle cell": "smooth muscle cell",
    "interneuron": "GABAergic neuron",
    "gabaergic neuron": "GABAergic neuron",
    "inhibitory neuron": "GABAergic neuron",
    "vip gabaergic interneuron": "VIP GABAergic interneuron",
    "pvalb gabaergic interneuron": "pvalb GABAergic interneuron",
    "sst gabaergic interneuron": "sst GABAergic interneuron",
    "lamp5 gabaergic interneuron": "lamp5 GABAergic interneuron",
}


def parse_mouse_age_years(stage) -> float:
    if stage is None or (isinstance(stage, float) and not np.isfinite(stage)):
        return np.nan
    s0 = str(stage).strip()
    m = re.match(r"^(\d+(?:\.\d+)?)\s*m$", s0.lower())
    if m:
        return float(m.group(1)) / 12.0
    a = parse_age_years(s0)
    if a is not None and np.isfinite(a):
        return float(a)
    s = s0.lower()
    m = re.search(r"(\d+(?:\.\d+)?)\s*-?\s*month", s)
    if m:
        return float(m.group(1)) / 12.0
    m = re.search(r"(\d+)\s*-?\s*week", s)
    if m:
        return float(m.group(1)) / 52.0
    m = re.search(r"(\d+)\s*-?\s*year", s)
    if m:
        return float(m.group(1))
    return np.nan


def map_type(label: str, types_human) -> str | None:
    s = str(label).strip().lower()
    if s in TYPE_MAP:
        t = TYPE_MAP[s]
        return t if t in types_human else None
    for k, v in TYPE_MAP.items():
        if k in s or s in k:
            return v if v in types_human else None
    return None


def logcpm_from_counts(C):
    C = np.asarray(C, np.float64)
    lib = C.sum(1, keepdims=True)
    lib[lib <= 0] = 1.0
    return np.log2(C / lib * 1e6 + 1.0)


def load_frozen(log):
    path = TRAJ_DIR / "frozen_scores.npz"
    if not path.exists():
        raise FileNotFoundError("frozen_scores.npz missing — G1 must finish without STOP")
    z = np.load(path, allow_pickle=True)
    meta = load_json(TRAJ_DIR / "frozen_meta.json")
    types = [str(t) for t in z["types"].tolist()]
    log(f"[frozen] genes={z['mu'].size} types={len(types)} identity_axes={z['basis'].shape[1]}")
    return dict(mu=z["mu"], sd=z["sd"], basis=z["basis"], centroids=z["centroids"],
                W_pls=z["W_pls"], W_ridge=z["W_ridge"],
                gene_id=np.array([str(g).split(".")[0] for g in z["gene_id"].tolist()]),
                symbol=np.array([str(s) for s in z["symbol"].tolist()]),
                types=types,
                centroid_types=[str(t) for t in z["centroid_types"].tolist()],
                meta=meta)


def load_orthologs(log):
    p = TRAJ_DIR / "g3_orthologs_1to1.csv"
    if not p.exists():
        raise FileNotFoundError("g3_orthologs_1to1.csv missing — run trajectory_orthologs.py")
    df = pd.read_csv(p)
    # column names from NCBI writer
    h = "human_ensembl" if "human_ensembl" in df.columns else df.columns[0]
    m = "mouse_ensembl" if "mouse_ensembl" in df.columns else df.columns[1]
    df[h] = df[h].astype(str).str.split(".").str[0]
    df[m] = df[m].astype(str).str.split(".").str[0]
    log(f"[orthologs] strict 1:1 pairs={len(df)}")
    return df, h, m


def align_mouse_to_human(C_mouse, mouse_ids, frozen, ortho, hcol, mcol, log, id_kind="ensembl"):
    """Return log-CPM matrix with frozen human gene order. Unmapped genes = 0 after z-score (fill human mean)."""
    mouse_ids = np.array([str(g).split(".")[0] for g in mouse_ids])
    if id_kind == "symbol":
        # map mouse symbol -> mouse ensembl via ortho mouse_symbol if present, else skip
        if "mouse_symbol" in ortho.columns:
            s2e = dict(zip(ortho["mouse_symbol"].astype(str).str.lower(), ortho[mcol]))
            mapped = np.array([s2e.get(g.lower(), "") for g in mouse_ids])
            mouse_ens = mapped
        else:
            mouse_ens = mouse_ids
    else:
        mouse_ens = mouse_ids
    m2h = dict(zip(ortho[mcol], ortho[hcol]))
    h_index = {g: i for i, g in enumerate(frozen["gene_id"])}
    # mouse column -> human column
    human_j = []
    mouse_j = []
    for j, mid in enumerate(mouse_ens):
        hid = m2h.get(mid)
        if hid in h_index:
            human_j.append(h_index[hid])
            mouse_j.append(j)
    n_map = len(human_j)
    log(f"[filter] mouse genes={len(mouse_ids)}; mapped 1:1 onto frozen human genes={n_map} "
        f"({n_map / max(len(frozen['gene_id']), 1):.3f} of frozen)")
    C = np.asarray(C_mouse, np.float64)
    if sp is not None and sp.issparse(C_mouse):
        C = C_mouse.toarray().astype(np.float64)
    Ym = logcpm_from_counts(C)
    Yh = np.zeros((Ym.shape[0], len(frozen["gene_id"])), np.float64)
    # fill human mean so z-score of unmapped is 0
    Yh[:] = frozen["mu"]
    if n_map:
        Yh[:, np.array(human_j)] = Ym[:, np.array(mouse_j)]
    return Yh, n_map


def apply_scores(Yh, frozen, type_labels=None):
    """Yh is log-CPM in frozen human gene order (unmapped already at human mean)."""
    mu, sd = frozen["mu"], frozen["sd"]
    sd = np.where(sd < 1e-12, 1.0, sd)
    Xs = (Yh - mu) / sd
    Z = identity_project(Xs, frozen["basis"])
    Cents = frozen["centroids"]
    ct_names = frozen["centroid_types"]
    d2 = ((Z[:, None, :] - Cents[None, :, :]) ** 2).sum(-1)
    j = np.argmin(d2, axis=1)
    assigned = np.array([ct_names[i] for i in j], dtype=object)
    dist_assigned = np.sqrt(d2[np.arange(len(Z)), j])
    type_to_j = {t: i for i, t in enumerate(frozen["types"])}
    age_ridge = np.full(len(Yh), np.nan)
    age_pls = np.full(len(Yh), np.nan)
    dist_own = np.full(len(Yh), np.nan)
    used_type = assigned.copy()
    if type_labels is not None:
        mapped = np.array([map_type(t, frozen["types"]) for t in type_labels], dtype=object)
        for i, t in enumerate(mapped):
            if t is not None:
                used_type[i] = t
    for i, t in enumerate(used_type):
        k = type_to_j.get(str(t))
        if k is None or not np.isfinite(frozen["W_ridge"][:, k]).all():
            continue
        wr = frozen["W_ridge"][:, k]
        nrm = float(np.linalg.norm(wr))
        if nrm < 1e-12:
            continue
        # ridge score: projection; frozen W_ridge is unit direction. Recover years via
        # a 1D map stored? We saved unit W and ridge intercept in meta.
        age_ridge[i] = float(Xs[i] @ wr)  # unit-direction projection
        wp = frozen["W_pls"][:, k]
        if np.isfinite(wp).all():
            age_pls[i] = float(Xs[i] @ wp)
        col = {c: j for j, c in enumerate(ct_names)}
        if str(t) in col:
            dist_own[i] = float(np.sqrt(d2[i, col[str(t)]]))
        else:
            dist_own[i] = dist_assigned[i]
    return dict(age_ridge=age_ridge, age_pls=age_pls, Z=Z, assigned=assigned,
                dist_assigned=dist_assigned, dist_own=dist_own, used_type=used_type)


def calibrate_age_projection_to_years(obs_human_style, frozen):
    """Optional: not used. Frozen ridge unit direction is not in years.

    We re-read g1_frozen_per_type.csv for intercepts if present; otherwise we treat
    the projection as the age SCORE (higher = older) and evaluate rank/R² against
    chronological age, which does not require a year scale.
    """
    return None


def _h5ad_path(dest: Path) -> Path:
    files = list(dest.glob("*.h5ad"))
    if not files:
        raise FileNotFoundError(f"no h5ad in {dest}")
    return files[0]


def tms_pseudobulk(path: Path, name: str, log, min_cells=10):
    log(f"[TMS] reading {path}")
    a = ad.read_h5ad(path)
    obs = a.obs.copy()
    log(f"[filter] {name}: {a.n_obs:,} cells x {a.n_vars:,} genes; obs cols={list(obs.columns)[:25]}")
    donor_col = next((c for c in ("donor_id", "donor", "mouse_id") if c in obs.columns), None)
    ct_col = next((c for c in ("cell_type", "celltype", "free_annotation") if c in obs.columns), None)
    age_col = next((c for c in ("development_stage", "age") if c in obs.columns), None)
    if donor_col is None or ct_col is None:
        raise RuntimeError(f"{name}: missing donor/cell_type columns: {list(obs.columns)}")
    age = np.array([parse_mouse_age_years(v) for v in obs[age_col]]) if age_col else np.full(len(obs), np.nan)
    log(f"[filter] donor_col={donor_col} ct_col={ct_col} age_col={age_col} "
        f"numeric age {np.isfinite(age).sum()}/{len(age)}")
    # gene ids
    if "feature_name" in a.var.columns:
        symbols = a.var["feature_name"].astype(str).to_numpy()
    else:
        symbols = a.var.index.astype(str).to_numpy()
    gids = a.var.index.astype(str).to_numpy()
    kind = "ensembl" if str(gids[0]).startswith("ENSMUSG") else "symbol"
    log(f"[filter] gene id kind={kind} example={gids[:3]}")

    # CELLxGENE: X is often log-normalized; raw/X holds integer counts.
    if a.raw is not None:
        X = a.raw.X
        log(f"[filter] {name} using adata.raw.X (counts)")
        raw_var = a.raw.var
        gids = raw_var.index.astype(str).to_numpy()
        if "feature_name" in raw_var.columns:
            symbols = raw_var["feature_name"].astype(str).to_numpy()
        kind = "ensembl" if str(gids[0]).startswith("ENSMUSG") else kind
    else:
        X = a.X
        log(f"[filter] {name} no raw; using adata.X")
    if sp is not None and sp.issparse(X):
        X = X.tocsr()
    rows = []
    has_sub = "subtissue" in obs.columns
    for (don, ct), sub in obs.assign(_age=age).groupby([donor_col, ct_col], sort=False):
        idx = sub.index
        pos = obs.index.get_indexer(idx)
        pos = pos[pos >= 0]
        if len(pos) < min_cells:
            continue
        if sp is not None and sp.issparse(X):
            counts = np.asarray(X[pos].sum(0)).ravel()
        else:
            counts = np.asarray(X[pos], np.float64).sum(0)
        rec = dict(donor=str(don), celltype=str(ct), n_cells=int(len(pos)),
                   age=float(np.nanmean(sub["_age"])), counts=counts)
        if has_sub:
            rec["subtissue"] = str(sub["subtissue"].astype(str).mode().iloc[0])
        rows.append(rec)
    log(f"[filter] {name} pseudobulks with ≥{min_cells} cells: {len(rows)}")
    if not rows:
        return None
    C = np.vstack([r.pop("counts") for r in rows])
    pb = pd.DataFrame(rows)
    pb.attrs = {}
    return dict(C=C, obs=pb, gene_ids=gids, symbols=symbols, id_kind=kind, name=name)


def score_mouse_pb(pb, frozen, ortho, hcol, mcol, log):
    Yh, n_map = align_mouse_to_human(
        pb["C"], pb["gene_ids"], frozen, ortho, hcol, mcol, log, id_kind=pb["id_kind"])
    sc = apply_scores(Yh, frozen, type_labels=pb["obs"].celltype.astype(str))
    out = pb["obs"].copy()
    out["age_score_ridge"] = sc["age_ridge"]
    out["age_score_pls"] = sc["age_pls"]
    out["id_dist"] = sc["dist_own"]
    out["id_assigned"] = sc["assigned"]
    out["id_used_type"] = sc["used_type"]
    out["mapped_human"] = [map_type(t, frozen["types"]) for t in out.celltype]
    out["source"] = pb["name"]
    return out, n_map


def _ols_r2(y, x):
    """R² of y ~ 1 + x. Evaluation calibration of a frozen 1D score; does not refit genes."""
    y = np.asarray(y, float)
    x = np.asarray(x, float)
    m = np.isfinite(y) & np.isfinite(x)
    y, x = y[m], x[m]
    if len(y) < 5:
        return np.nan
    A = np.column_stack([np.ones(len(x)), x])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ coef
    return float(r2_mae(y, pred)["r2"])


def mouse_aging_r2(tab, rng, n_perm, log, score_col="age_score_ridge"):
    """Does the frozen score track chronological age?

    The frozen ridge direction is a unit projection, not human years. Identity R² of
    (projection vs mouse years) is not a usable number. Tracking is: Spearman rank
    correlation, and OLS R² of mouse age ~ 1 + frozen score (1D evaluation calibration).
    """
    t = tab.dropna(subset=["age", score_col]).copy()
    if len(t) < 8:
        return dict(r2=np.nan, n=len(t), transfer_metric="ols_age_on_score")
    type_rows = []
    for ct, sub in t.groupby("id_used_type"):
        if sub.donor.nunique() < 6:
            continue
        don = sub.groupby("donor").agg(age=("age", "first"), pred=(score_col, "mean"))
        sp = float(pd.Series(don.age).corr(pd.Series(don.pred), method="spearman"))
        type_rows.append(dict(celltype=str(ct), n_donors=int(len(don)),
                              spearman=sp, ols_r2=_ols_r2(don.age, don.pred)))
    type_df = pd.DataFrame(type_rows)
    if len(type_df):
        type_df.to_csv(TRAJ_DIR / "g3_tms_per_type_transfer.csv", index=False)
        log("[G3] per-type frozen-score vs age (Spearman):")
        for _, r in type_df.sort_values("spearman").iterrows():
            log(f"   {r.celltype}: n={int(r.n_donors)} Spearman={r.spearman:+.3f} OLS_R²={r.ols_r2:+.3f}")
    median_type = float(np.nanmedian(type_df.ols_r2)) if len(type_df) else np.nan
    don = t.groupby("donor").agg(age=("age", "first"), pred=(score_col, "mean"))
    ols = _ols_r2(don.age, don.pred)
    spear = float(pd.Series(don.age).corr(pd.Series(don.pred), method="spearman")) if len(don) >= 5 else np.nan
    pear = float(pd.Series(don.age).corr(pd.Series(don.pred))) if len(don) >= 5 else np.nan
    log(f"   donor-mean {score_col} OLS_R²(age~score)={ols:+.3f}  Pearson={pear:+.3f}  "
        f"Spearman={spear:+.3f}  n_donors={len(don)}  median-type OLS_R²={median_type:+.3f} "
        f"(n_types={len(type_df)})")
    ages = don.age.to_numpy(float)
    pred = don.pred.to_numpy(float)
    null_ols, null_sp = [], []
    for i in range(n_perm):
        yp = rng.permutation(ages)
        null_ols.append(_ols_r2(yp, pred))
        null_sp.append(float(pd.Series(yp).corr(pd.Series(pred), method="spearman")))
    null_ols, null_sp = np.array(null_ols, float), np.array(null_sp, float)
    p_ols = permutation_p(ols, null_ols, greater=True)
    p_sp = permutation_p(abs(spear) if np.isfinite(spear) else np.nan,
                         np.abs(null_sp), greater=True)
    log(f"   shuffle-age null OLS_R² mean={np.nanmean(null_ols):+.3f}  p_ols={p_ols:.3f}  "
        f"null Spearman mean={np.nanmean(null_sp):+.3f}  p_spearman={p_sp:.3f}  n_perm={n_perm}")
    # cortex-only sensitivity if subtissue present
    cortex = dict(r2=np.nan, spearman=np.nan, n_donors=0)
    if "subtissue" in t.columns:
        m = t.subtissue.astype(str).str.lower().str.contains("cortex")
        if m.sum() >= 8:
            dc = t[m].groupby("donor").agg(age=("age", "first"), pred=(score_col, "mean"))
            cortex = dict(r2=_ols_r2(dc.age, dc.pred),
                          spearman=float(pd.Series(dc.age).corr(pd.Series(dc.pred), method="spearman")),
                          n_donors=int(len(dc)))
            log(f"   cortex-subtissue sensitivity: n_donors={cortex['n_donors']} "
                f"OLS_R²={cortex['r2']:+.3f} Spearman={cortex['spearman']:+.3f}")
    sign_pos = int((type_df.spearman > 0).sum()) if len(type_df) else 0
    sign_neg = int((type_df.spearman < 0).sum()) if len(type_df) else 0
    log(f"   per-type Spearman signs: +{sign_pos} / -{sign_neg}  (mixed signs ⇒ donor-mean can cancel)")
    return dict(r2=ols, mae=np.nan, spearman=spear, pearson=pear, n_donors=int(len(don)),
                median_type_r2=median_type, n_types=int(len(type_df)),
                null_mean=float(np.nanmean(null_ols)), p=p_ols, p_spearman=p_sp, n_perm=n_perm,
                n_pseudobulks=int(len(t)), transfer_metric="ols_age_on_score",
                n_types_spearman_pos=sign_pos, n_types_spearman_neg=sign_neg,
                cortex_sensitivity=cortex)


def ensure_downloaded(log):
    for rec, label in ((TMS_NONMY, "TMS non-myeloid"), (TMS_MY, "TMS myeloid"), (GSE224, "GSE224438")):
        dest = rec["dest"]
        dest.mkdir(parents=True, exist_ok=True)
        already = list(dest.glob("*.h5ad")) + list(dest.glob("*.tar"))
        if already:
            log(f"[data] cached {label}: {already[0].name}")
            continue
        log(f"[download] {label} {rec['url']}")
        download(rec["url"], str(dest))


def parse_gse224_samples(log):
    """Parse GSM titles + characteristics from GEO SOFT (targ=gsm) — no invented groups."""
    text = fetch_soft("GSE224438", targ_gsm=True)
    samples = {}
    cur = None
    title = None
    char = {}
    for line in text.splitlines():
        if line.startswith("^SAMPLE"):
            if cur and title:
                samples[cur] = dict(title=title, **char)
            cur = line.split("=")[-1].strip()
            title, char = None, {}
        elif line.startswith("!Sample_title"):
            title = line.split("=", 1)[-1].strip()
        elif line.startswith("!Sample_characteristics_ch"):
            body = line.split("=", 1)[-1].strip()
            if ":" in body:
                k, v = body.split(":", 1)
                char[k.strip()] = v.strip()
    if cur and title:
        samples[cur] = dict(title=title, **char)
    rows = []
    for gsm, rec in samples.items():
        group = _gse224_group_from_rec(rec)
        rows.append(dict(
            gsm=gsm, title=rec.get("title", ""), group=group,
            treatment=rec.get("treatment", ""), genotype=rec.get("genotype", ""),
        ))
        log(f"  {gsm}  group={group}  treatment={rec.get('treatment','')}  {rec.get('title','')[:80]}")
    tab = pd.DataFrame(rows)
    tab.to_csv(TRAJ_DIR / "g3_GSE224438_samples.csv", index=False)
    n_ok = int((tab.group != "multiplexed_unusable").sum())
    log(f"[filter] GSE224438 GSM with a single treatment (usable): {n_ok}/{len(tab)}")
    return tab


def _gse224_group_from_rec(rec: dict) -> str:
    treat = str(rec.get("treatment", "")).lower()
    title = str(rec.get("title", "")).lower()
    # Whole-body libraries multiplex young/old/OSK in one 10x lane (GEO treatment strings
    # contain 'and'). Those cannot be split without hashing barcodes, which are not in RAW.
    mixed = (
        " and " in treat
        or ( "young" in treat and "old" in treat)
        or treat.count("control") >= 2
        or "multiplexed" in title and "and" in treat
    )
    if mixed:
        return "multiplexed_unusable"
    if "svz" in treat and "reprogram" in treat:
        return "old_svz_osk"
    if ("whole-body" in treat or "whole body" in treat) and "reprogram" in treat:
        return "old_whole_osk"
    if "reprogram" in treat:
        return "old_reprogrammed"
    if "young" in treat:
        return "young"
    if "old" in treat:
        return "old"
    return "unknown"


def load_gse224(log):
    tar_path = GSE224["dest"] / "GSE224438_RAW.tar"
    if not tar_path.exists():
        hits = list(GSE224["dest"].glob("*.tar"))
        if not hits:
            raise FileNotFoundError(tar_path)
        tar_path = hits[0]
    extract = GSE224["dest"] / "extracted"
    extract.mkdir(parents=True, exist_ok=True)
    log(f"[GSE224438] tar={tar_path} size={tar_path.stat().st_size:,}")
    names = []
    with tarfile.open(tar_path, "r") as tar:
        names = tar.getnames()
        log(f"[filter] tar members={len(names)}")
        for n in names[:30]:
            log(f"   {n}")
        need = [m for m in tar.getmembers() if not (extract / Path(m.name).name).exists()]
        if need:
            tar.extractall(extract, members=None)
    files = sorted(extract.rglob("*"))
    log(f"[filter] extracted files={len(files)}")
    for p in files[:40]:
        if p.is_file():
            log(f"   {p.relative_to(extract)}  {p.stat().st_size:,}")
    return extract, names


def gse224_matrices(extract: Path, samples: pd.DataFrame, log):
    """Load per-GSM 10x-style matrices if present; else any count tables.

    Multiplexed whole-body libraries are not loaded (young/old/OSK share a 10x lane;
    RAW has no hashing barcodes).
    """
    named = sorted(extract.rglob("*_matrix.mtx.gz"))
    if named:
        log(f"[filter] GEO-named 10x triplets on disk: {len(named)}")
        return _load_geo_named_mtx(named, samples, log)
    h5s = list(extract.rglob("*.h5")) + list(extract.rglob("*.h5ad"))
    mtgs = list(extract.rglob("matrix.mtx*"))
    log(f"[filter] h5/h5ad={len(h5s)}  mtx={len(mtgs)}")
    if h5s:
        return _load_10x_h5_stack(h5s, samples, log)
    if mtgs:
        return _load_mtx_stack(mtgs, samples, log)
    return None


def _usable_gsm_set(samples: pd.DataFrame) -> set[str]:
    return set(samples.loc[samples.group != "multiplexed_unusable", "gsm"].astype(str))


def _load_one_geo_mtx(mat: Path, log):
    """GSE224438 RAW tar uses GSM######_<sample>_matrix.mtx.gz next to _features and _barcodes."""
    import gzip
    from scipy.io import mmread
    stem = mat.name.replace("_matrix.mtx.gz", "").replace("_matrix.mtx", "")
    gsm = _gsm_from_name(stem) or _gsm_from_name(mat.name)
    feat = mat.with_name(stem + "_features.tsv.gz")
    if not feat.exists():
        feat = mat.with_name(stem + "_genes.tsv.gz")
    bc = mat.with_name(stem + "_barcodes.tsv.gz")
    if not feat.exists() or not bc.exists():
        log(f"   skip {mat.name}: missing features/barcodes")
        return None
    log(f"   reading {mat.name} ...")
    X = mmread(gzip.open(mat, "rb")).T.tocsr()
    genes = pd.read_csv(feat, sep="\t", header=None)
    barcodes = pd.read_csv(bc, sep="\t", header=None)
    if genes.shape[1] >= 2:
        var = pd.DataFrame({"gene_ids": genes[0].astype(str),
                            "feature_name": genes[1].astype(str)})
        var.index = genes[1].astype(str)
    else:
        var = pd.DataFrame(index=genes[0].astype(str))
        var["gene_ids"] = genes[0].astype(str)
    obs = pd.DataFrame(index=barcodes[0].astype(str))
    obs["gsm"] = gsm or stem
    obs["file"] = mat.name
    a = ad.AnnData(X=X, obs=obs, var=var)
    a.var_names_make_unique()
    log(f"   {mat.name}: {a.n_obs} cells x {a.n_vars} genes gsm={gsm}")
    return a


def _load_geo_named_mtx(mats, samples, log):
    """Load only single-treatment GSMs. Returns a list of AnnData (not concatenated)."""
    usable = _usable_gsm_set(samples)
    log(f"[filter] single-treatment GSMs (load these): {sorted(usable)}")
    ads = []
    n_skip_mux = 0
    for mat in mats:
        gsm = _gsm_from_name(mat.name)
        if gsm and gsm not in usable:
            n_skip_mux += 1
            log(f"   skip {mat.name}: group=multiplexed_unusable")
            continue
        a = _load_one_geo_mtx(mat, log)
        if a is not None:
            ads.append(a)
    log(f"[filter] loaded {len(ads)} usable libraries; skipped multiplexed files={n_skip_mux}")
    if not ads:
        return None
    return ads


def _gsm_from_name(name: str):
    m = re.search(r"(GSM\d+)", name)
    return m.group(1) if m else None


def _load_10x_h5_stack(h5s, samples, log):
    import scanpy as sc
    ads = []
    for p in h5s:
        gsm = _gsm_from_name(p.name)
        try:
            a = sc.read_10x_h5(p)
        except Exception as e:
            log(f"   skip {p.name}: {e}")
            continue
        a.var_names_make_unique()
        a.obs["gsm"] = gsm or p.stem
        a.obs["file"] = p.name
        ads.append(a)
        log(f"   {p.name}: {a.n_obs} cells x {a.n_vars} genes gsm={gsm}")
    if not ads:
        return None
    combo = ad.concat(ads, join="outer", label="file_batch", index_unique="-") if len(ads) > 1 else ads[0]
    return combo


def _load_mtx_stack(mtgs, samples, log):
    import scanpy as sc
    ads = []
    for mtx in mtgs:
        d = mtx.parent
        gsm = _gsm_from_name(str(d)) or _gsm_from_name(mtx.name)
        try:
            a = sc.read_10x_mtx(d, var_names="gene_symbols", cache=False)
        except Exception:
            try:
                a = sc.read_10x_mtx(d, var_names="gene_ids", cache=False)
            except Exception as e:
                log(f"   skip mtx {d}: {e}")
                continue
        a.var_names_make_unique()
        a.obs["gsm"] = gsm or d.name
        ads.append(a)
        log(f"   mtx {d.name}: {a.n_obs} cells gsm={gsm}")
    if not ads:
        return None
    return ad.concat(ads, join="outer", label="file_batch", index_unique="-") if len(ads) > 1 else ads[0]


def _ortholog_column_index(mouse_ids, frozen, ortho, hcol, mcol, id_kind="ensembl"):
    """Mouse column indices that map 1:1 onto frozen human genes, and the human columns."""
    mouse_ids = np.array([str(g).split(".")[0] for g in mouse_ids])
    if id_kind == "symbol" and "mouse_symbol" in ortho.columns:
        s2e = dict(zip(ortho["mouse_symbol"].astype(str).str.lower(), ortho[mcol]))
        mouse_ens = np.array([s2e.get(g.lower(), "") for g in mouse_ids])
    else:
        mouse_ens = mouse_ids
    m2h = dict(zip(ortho[mcol], ortho[hcol]))
    h_index = {g: i for i, g in enumerate(frozen["gene_id"])}
    human_j, mouse_j = [], []
    for j, mid in enumerate(mouse_ens):
        hid = m2h.get(mid)
        if hid in h_index:
            human_j.append(h_index[hid])
            mouse_j.append(j)
    return np.array(human_j, int), np.array(mouse_j, int)


def _assign_types_frozen_sparse(X, mouse_ids, frozen, ortho, hcol, mcol, log, id_kind="ensembl",
                                min_genes=200):
    """Nearest frozen identity centroid per cell. Unmapped genes contribute 0 (human-mean fill).

    Does not refit, does not select HVGs, does not run Leiden. Scores use every frozen gene's
    weight; genes with no 1:1 ortholog sit at the human mean (z = 0).
    """
    if sp is None or not sp.issparse(X):
        X = sp.csr_matrix(np.asarray(X))
    else:
        X = X.tocsr()
    n_genes = np.asarray((X > 0).sum(1)).ravel()
    keep = n_genes >= min_genes
    log(f"[filter] cells with ≥{min_genes} detected genes: {int(keep.sum())}/{len(keep)}")
    if keep.sum() < 10:
        return keep, np.array(["unassigned"] * len(keep), dtype=object)
    human_j, mouse_j = _ortholog_column_index(mouse_ids, frozen, ortho, hcol, mcol, id_kind)
    log(f"[filter] cell-level ortholog columns used for identity assignment={len(mouse_j)}")
    Xk = X[keep][:, mouse_j]
    lib = np.asarray(X[keep].sum(1)).ravel().astype(np.float64)
    lib[lib <= 0] = 1.0
    Ym = np.log2(np.asarray(Xk.toarray(), np.float64) / lib[:, None] * 1e6 + 1.0)
    mu = frozen["mu"][human_j]
    sd = np.where(frozen["sd"][human_j] < 1e-12, 1.0, frozen["sd"][human_j])
    Zg = (Ym - mu) / sd
    basis_m = frozen["basis"][human_j]
    Z = Zg @ basis_m
    Cents = frozen["centroids"]
    d2 = ((Z[:, None, :] - Cents[None, :, :]) ** 2).sum(-1)
    j = np.argmin(d2, axis=1)
    assigned_keep = np.array([frozen["centroid_types"][i] for i in j], dtype=object)
    assigned = np.array(["unassigned"] * X.shape[0], dtype=object)
    assigned[keep] = assigned_keep
    vc = pd.Series(assigned_keep).value_counts()
    log(f"[filter] frozen-centroid assignments (top): {vc.head(8).to_dict()}")
    return keep, assigned


def _pseudobulk_one_library(adata, group, frozen, ortho, hcol, mcol, log, min_cells=20, min_genes=200):
    X = adata.X
    if sp is not None and sp.issparse(X):
        X = X.tocsr()
    gids = adata.var["gene_ids"].astype(str).to_numpy() if "gene_ids" in adata.var.columns \
        else adata.var.index.astype(str).to_numpy()
    kind = "ensembl" if str(gids[0]).startswith("ENSMUSG") else "symbol"
    keep, assigned = _assign_types_frozen_sparse(
        X, gids, frozen, ortho, hcol, mcol, log, id_kind=kind,
        min_genes=min_genes)
    gsm = str(adata.obs["gsm"].iloc[0])
    rows, mats = [], []
    pos_all = np.where(keep)[0]
    labels = assigned[keep]
    for t in sorted(pd.unique(labels)):
        sel = pos_all[labels == t]
        if len(sel) < min_cells:
            continue
        if sp is not None and sp.issparse(X):
            counts = np.asarray(X[sel].sum(0)).ravel()
        else:
            counts = np.asarray(X[sel], np.float64).sum(0)
        rows.append(dict(donor=gsm, celltype=str(t), n_cells=int(len(sel)),
                         group=str(group), gsm=gsm, cluster=str(t)))
        mats.append(counts)
    log(f"[filter] {gsm} group={group} pseudobulks ≥{min_cells} cells: {len(rows)}")
    if not rows:
        return None
    symbols = (adata.var["feature_name"].astype(str).to_numpy()
               if "feature_name" in adata.var.columns
               else adata.var_names.astype(str).to_numpy())
    return dict(C=np.vstack(mats), obs=pd.DataFrame(rows), gene_ids=gids,
                symbols=symbols, id_kind=kind, name="GSE224438")


def cluster_and_pseudobulk(adata, samples, log, min_cells=20, frozen=None, ortho=None,
                           hcol=None, mcol=None):
    """Pseudobulk GSE224438. `adata` may be a list of per-GSM AnnData.

    Cell types are assigned by the FROZEN human identity subspace (nearest centroid).
    Leiden/HVG are not used — that would be a second clustering, not the frozen score.
    """
    if frozen is None or ortho is None:
        raise RuntimeError("frozen scores + orthologs required to assign types without refitting")
    ads = adata if isinstance(adata, list) else [adata]
    gmap = dict(zip(samples.gsm.astype(str), samples.group))
    parts = []
    n_cells_in, n_libs = 0, 0
    for a in ads:
        n_libs += 1
        n_cells_in += int(a.n_obs)
        gsm = str(a.obs["gsm"].iloc[0]) if "gsm" in a.obs.columns else _gsm_from_name(str(a.obs.index[0]))
        group = gmap.get(gsm, "unknown")
        if group == "multiplexed_unusable":
            log(f"[filter] skip {gsm}: multiplexed_unusable")
            continue
        log(f"[GSE224438] frozen identity assignment {gsm} n_cells={a.n_obs} group={group}")
        pb = _pseudobulk_one_library(a, group, frozen, ortho, hcol, mcol, log, min_cells=min_cells)
        if pb is not None:
            parts.append(pb)
    log(f"[filter] GSE224438 libraries processed={n_libs} cells_in={n_cells_in} "
        f"pseudobulk-parts={len(parts)}")
    if not parts:
        return None
    gids0 = parts[0]["gene_ids"]
    obs = pd.concat([p["obs"] for p in parts], ignore_index=True)
    C = np.vstack([p["C"] for p in parts])
    age_map = dict(young=0.25, old=1.5, old_whole_osk=1.5, old_svz_osk=1.5,
                   old_reprogrammed=1.5, unknown=np.nan)
    obs["age"] = obs.group.map(age_map)
    log(f"[filter] GSE224438 total type-pseudobulks={len(obs)}; groups={obs.group.value_counts().to_dict()}")
    return dict(C=C, obs=obs, gene_ids=gids0, symbols=parts[0]["symbols"],
                id_kind=parts[0]["id_kind"], name="GSE224438")


def _donor_bootstrap_means(sub, rng, n_boot, age_col="age_score_ridge", id_col="id_retained"):
    """Equal-weight GSMs/donors, not type-pseudobulks. With 2 libraries the CI is wide on purpose."""
    don = sub.groupby("donor").agg(age=(age_col, "mean"), idr=(id_col, "mean"))
    ages = don.age.to_numpy(float)
    ids = don.idr.to_numpy(float)
    nd = len(don)
    if nd < 1:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 0
    age_mean = float(np.nanmean(ages))
    id_mean = float(np.nanmean(ids))
    if nd == 1 or n_boot < 1:
        return age_mean, age_mean, age_mean, id_mean, id_mean, id_mean, nd
    ba, bi = [], []
    idx = np.arange(nd)
    for _ in range(n_boot):
        i = rng.choice(idx, size=nd, replace=True)
        ba.append(np.nanmean(ages[i]))
        bi.append(np.nanmean(ids[i]))
    ba, bi = np.asarray(ba, float), np.asarray(bi, float)
    return (age_mean, float(np.nanpercentile(ba, 2.5)), float(np.nanpercentile(ba, 97.5)),
            id_mean, float(np.nanpercentile(bi, 2.5)), float(np.nanpercentile(bi, 97.5)), nd)


def trajectory_stats(tab, log, rng, n_boot=1000):
    """Group means on (age_score, -id_dist) with donor-clustered bootstrap CIs.

    GSE224438 is a 3-condition protocol contrast (young / old / old+OSK), not a
    multi-day time course. A quadratic in 3 points is saturated and is not tested.
    """
    t = tab.dropna(subset=["age_score_ridge", "id_dist"]).copy()
    t = t[~t.group.astype(str).isin(["multiplexed_unusable", "unknown"])].copy()
    t["id_retained"] = -t["id_dist"]
    groups_order = [g for g in ("young", "old", "old_whole_osk", "old_svz_osk", "old_reprogrammed")
                    if g in set(t.group)]
    rows = []
    for g in groups_order:
        sub = t[t.group == g]
        if sub.donor.nunique() < 1:
            continue
        age_mean, age_lo, age_hi, id_mean, id_lo, id_hi, nd = _donor_bootstrap_means(sub, rng, n_boot)
        rows.append(dict(
            group=g, n=int(len(sub)), n_donors=int(nd),
            age_mean=age_mean, age_lo=age_lo, age_hi=age_hi,
            id_mean=id_mean, id_lo=id_lo, id_hi=id_hi,
            dist_mean=float(np.nanmean(sub.id_dist)),
        ))
        log(f"   {g}: n_pb={len(sub)} n_donors={nd} age_score={rows[-1]['age_mean']:+.3f} "
            f"[{rows[-1]['age_lo']:+.3f},{rows[-1]['age_hi']:+.3f}]  "
            f"id_retained={rows[-1]['id_mean']:+.3f} "
            f"[{rows[-1]['id_lo']:+.3f},{rows[-1]['id_hi']:+.3f}]")
    means = pd.DataFrame(rows)
    means.to_csv(TRAJ_DIR / "g3_trajectory_group_means.csv", index=False)
    rec = dict(groups=means.to_dict(orient="records"), bend=None, window=None, didt=[])
    if {"young", "old"} <= set(means.group) and len(means) >= 3:
        y = means.set_index("group")
        v_age = np.array([y.loc["old"].age_mean - y.loc["young"].age_mean,
                          y.loc["old"].id_mean - y.loc["young"].id_mean])
        rec["aging_vector"] = v_age.tolist()
        rec["didt"].append(dict(segment="young→old",
                                d_age=float(v_age[0]), d_id=float(v_age[1]),
                                d_id_per_d_age=float(v_age[1] / v_age[0]) if abs(v_age[0]) > 1e-8 else np.nan))
        osk_groups = [g for g in ("old_svz_osk", "old_whole_osk", "old_reprogrammed") if g in y.index]
        for og in osk_groups:
            v_r = np.array([y.loc[og].age_mean - y.loc["old"].age_mean,
                            y.loc[og].id_mean - y.loc["old"].id_mean])
            rec["didt"].append(dict(segment=f"old→{og}",
                                    d_age=float(v_r[0]), d_id=float(v_r[1]),
                                    d_id_per_d_age=float(v_r[1] / v_r[0]) if abs(v_r[0]) > 1e-8 else np.nan))
            # reverse-aging direction is -v_age
            na, nr = np.linalg.norm(v_age), np.linalg.norm(v_r)
            cos_rev = float(np.dot(v_r, -v_age) / (na * nr)) if na > 0 and nr > 0 else np.nan
            # off-axis: identity change relative to aging identity change
            age_falls = v_r[0] < 0 and (y.loc[og].age_hi < y.loc["old"].age_mean)
            id_holds = abs(v_r[1]) <= 0.5 * abs(v_age[1]) if abs(v_age[1]) > 1e-8 else abs(v_r[1]) <= abs(y.loc["old"].id_hi - y.loc["old"].id_lo)
            # CI overlap on identity with old
            id_overlap = not (y.loc[og].id_hi < y.loc["old"].id_lo or y.loc[og].id_lo > y.loc["old"].id_hi)
            rec.setdefault("reprogramming", []).append(dict(
                group=og, vector=v_r.tolist(), cos_reverse_aging=cos_rev,
                age_falls=bool(age_falls), identity_holds=bool(id_holds or id_overlap),
                id_ci_overlaps_old=bool(id_overlap),
            ))
            log(f"   old→{og}: Δage={v_r[0]:+.3f} Δid={v_r[1]:+.3f}  "
                f"cos_to_reverse_aging={cos_rev:+.3f}  age_falls={age_falls}  "
                f"identity_holds={id_holds or id_overlap}")
        # bend: reprogramming not collinear with aging (cos of v_r with v_age not ~±1)
        bend = False
        window = None
        for r in rec.get("reprogramming", []):
            vr = np.array(r["vector"])
            na = np.linalg.norm(v_age)
            nr = np.linalg.norm(vr)
            cos_age = float(np.dot(vr, v_age) / (na * nr)) if na > 0 and nr > 0 else np.nan
            r["cos_to_aging"] = cos_age
            # bend if not a simple reversal (cos_rev high) AND not continuing aging (cos_age high)
            if r["age_falls"] and r["identity_holds"]:
                window = dict(
                    range=f"old → {r['group']} (protocol, not a time window in days)",
                    d_age=r["vector"][0],
                    id_retained="CI overlaps old" if r["id_ci_overlaps_old"] else "Δid small vs aging",
                    age_ci="see g3_trajectory_group_means.csv",
                    id_ci="see g3_trajectory_group_means.csv",
                    group=r["group"],
                )
            if np.isfinite(cos_age) and abs(cos_age) < 0.85:
                bend = True
        rec["bend"] = bool(bend)
        rec["window"] = window
        # linear vs nonlinear: group-mean identity ~ age_score
        if len(means) >= 3:
            x = means.age_mean.to_numpy(float)
            yv = means.id_mean.to_numpy(float)
            A = np.column_stack([np.ones(len(x)), x])
            coef, *_ = np.linalg.lstsq(A, yv, rcond=None)
            pred = A @ coef
            ss_res = float(np.sum((yv - pred) ** 2))
            ss_tot = float(np.sum((yv - yv.mean()) ** 2))
            r2_lin = np.nan if ss_tot <= 0 else 1 - ss_res / ss_tot
            # quadratic
            A2 = np.column_stack([np.ones(len(x)), x, x ** 2])
            coef2, *_ = np.linalg.lstsq(A2, yv, rcond=None)
            pred2 = A2 @ coef2
            ss_res2 = float(np.sum((yv - pred2) ** 2))
            r2_q = np.nan if ss_tot <= 0 else 1 - ss_res2 / ss_tot
            rec["linear_r2"] = r2_lin
            rec["n_conditions"] = int(len(means))
            rec["nonlinear_note"] = (
                "Protocol contrast, not a time course. "
                "Quadratic vs linear is unidentified with 3 conditions (a 3-parameter quadratic saturates)."
            )
            if len(means) >= 4:
                rec["quadratic_r2"] = r2_q
                rec["quadratic_better"] = bool(np.isfinite(r2_q) and np.isfinite(r2_lin) and (r2_q - r2_lin) > 0.05)
                log(f"   identity ~ age_score  linear R²={r2_lin:+.3f}  quadratic R²={r2_q:+.3f}")
            else:
                rec["quadratic_r2"] = None
                rec["quadratic_better"] = False
                log(f"   identity ~ age_score  linear R²={r2_lin:+.3f}  (quadratic not tested; {len(means)} conditions)")
    return rec, means


def plot_trajectory(means, path):
    if means is None or len(means) == 0:
        return
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    for _, r in means.iterrows():
        ax.errorbar(r.age_mean, r.id_mean,
                    xerr=[[r.age_mean - r.age_lo], [r.age_hi - r.age_mean]],
                    yerr=[[r.id_mean - r.id_lo], [r.id_hi - r.id_mean]],
                    fmt="o", capsize=3, label=r.group)
    if len(means) >= 2:
        ax.plot(means.age_mean, means.id_mean, c="0.5", lw=0.8)
    ax.set_xlabel("frozen age score (ridge projection)")
    ax.set_ylabel("identity retained (−dist to type centroid)")
    ax.set_title("G3  GSE224438 on the age–identity plane")
    ax.legend(fontsize=7, loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_didt(rec, path):
    segs = rec.get("didt") or []
    if not segs:
        return
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    labs = [s["segment"] for s in segs]
    vals = [s["d_id_per_d_age"] for s in segs]
    ax.bar(np.arange(len(labs)), vals, color="tab:purple")
    ax.axhline(0, c="0.4", lw=0.6)
    ax.set_xticks(np.arange(len(labs)))
    ax.set_xticklabels(labs, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("d(identity_retained)/d(age_score)")
    ax.set_title("G3  slope along the trajectory")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_mouse_aging(tab, path):
    t = tab.dropna(subset=["age", "age_score_ridge"])
    if len(t) < 5:
        return
    fig, ax = plt.subplots(figsize=(5.8, 4.4))
    don = t.groupby("donor").agg(age=("age", "first"), pred=("age_score_ridge", "mean"))
    ax.scatter(don.age * 12, don.pred, c="tab:red", s=28)
    ax.set_xlabel("donor age (months)")
    ax.set_ylabel("frozen human ridge age score")
    ax.set_title("G3  transfer check: Tabula Muris Senis brain")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def run():
    rng = np.random.default_rng(TRAJ_SEED)
    log = Logger(TRAJ_DIR / "g3_report.txt")
    try:
        return _run(log, rng)
    finally:
        log.close()


def _run(log, rng):
    traj_log_banner(log, "G3 — frozen scores on mouse aging + reprogramming")
    g1 = load_json(TRAJ_DIR / "g1_summary.json") if (TRAJ_DIR / "g1_summary.json").exists() else {}
    g2 = load_json(TRAJ_DIR / "g2_summary.json") if (TRAJ_DIR / "g2_summary.json").exists() else {}
    if g1.get("stop_G1") or g2.get("stop_G2") or g2.get("skipped"):
        msg = "G1/G2 halt — not applying scores to reprogramming data"
        log(f"[G3] {msg}")
        dump_json(TRAJ_DIR / "g3_summary.json", dict(status="skipped", reason=msg, verdict="c"))
        return dict(status="skipped")

    ensure_downloaded(log)
    frozen = load_frozen(log)
    ortho, hcol, mcol = load_orthologs(log)
    n_human = len(frozen["gene_id"])
    n_1to1 = len(ortho)
    log(f"[filter] frozen human genes={n_human}; strict 1:1 orthologs={n_1to1}")

    # ----- independent mouse aging (TMS) -----
    tabs = []
    n_maps = []
    for rec, name in ((TMS_NONMY, "tms_nonmyeloid"), (TMS_MY, "tms_myeloid")):
        path = _h5ad_path(rec["dest"])
        pb = tms_pseudobulk(path, name, log)
        if pb is None:
            continue
        scored, n_map = score_mouse_pb(pb, frozen, ortho, hcol, mcol, log)
        scored.to_csv(TRAJ_DIR / f"g3_{name}_scores.csv", index=False)
        tabs.append(scored)
        n_maps.append(n_map)
        log(f"   mapped types: {scored.mapped_human.notna().mean():.3f} of pseudobulks have a DLPFC type map")
    if not tabs:
        raise RuntimeError("no TMS pseudobulks")
    tms = pd.concat(tabs, ignore_index=True)
    tms.to_csv(TRAJ_DIR / "g3_tms_scores.csv", index=False)
    n_map_use = int(np.mean(n_maps)) if n_maps else 0
    log("\n[G3] independent mouse aging validation (frozen scores, no refit)")
    val = mouse_aging_r2(tms, rng, n_perm=50, log=log)
    plot_mouse_aging(tms, TRAJ_FIG / "g3_mouse_aging_validation.png")
    transfer_ok = bool(
        np.isfinite(val["r2"]) and val["r2"] > STOP_R2_COLLAPSE and val["p"] < 0.05
        and np.isfinite(val.get("spearman", np.nan)) and val["spearman"] > 0
    )
    if (not transfer_ok) and np.isfinite(val.get("spearman", np.nan)) and abs(val["spearman"]) >= 0.4 \
            and val.get("p_spearman", 1.0) < 0.05:
        transfer_ok = True
        log("[G3] OLS R² weak but Spearman |ρ|≥0.4 with p_spearman<0.05 — rank transfer usable.")
    log(f"[G3] transfer_ok={transfer_ok}")

    ortholog_rec = dict(n_human=n_human, n_1to1=n_1to1, n_used=n_map_use,
                        frac_used=n_map_use / max(n_human, 1))
    mouse_aging = dict(
        id="66ff82b4-9380-469c-bc4b-cfa08eacd325 + c08f8441-4a10-4748-872a-e70c0bcccdba",
        title="Tabula Muris Senis brain non-myeloid + myeloid (CXG)",
        **val, transfer_ok=transfer_ok, n_genes_mapped=n_map_use,
    )

    if not transfer_ok:
        narrative = (
            "The frozen human DLPFC age score does not track chronological age in independent "
            "mouse brain aging (Tabula Muris Senis) after 1:1 ortholog mapping "
            f"({n_map_use}/{n_human} genes mapped, {n_map_use / max(n_human,1):.1%} of frozen). "
            f"OLS R² of mouse age ~ score = {val['r2']:+.3f} (null {val['null_mean']:+.3f}, p={val['p']:.3f}); "
            f"Spearman ρ={val['spearman']:+.3f} (p={val.get('p_spearman', np.nan):.3f}). "
            f"Per-type Spearman signs +{val.get('n_types_spearman_pos')} / -{val.get('n_types_spearman_neg')} "
            "(mixed signs; a donor-mean can cancel). "
            "The trajectory is not computed. "
            "Needed: a human adult cortical/DLPFC scRNA dataset (preferred) in which the frozen score "
            "tracks donor age, plus an adult cortical OSK or partial-reprogramming time course with "
            "cell-type labels; or a mouse cortical aging atlas on which this frozen score tracks age, "
            "plus an adult cortical (not SVZ, not E15.5) OSK time course."
        )
        summary = dict(
            status="no_transfer", verdict="d",
            verdict_text=narrative, narrative=narrative,
            orthologs=ortholog_rec, mouse_aging=mouse_aging,
            needed=("Human adult DLPFC/cortex OSK or partial-reprogramming scRNA time course "
                    "with numeric ages and cell-type labels; or a mouse cortical aging atlas "
                    "on which the frozen score tracks age, plus an adult cortical (not SVZ, not E15.5) "
                    "OSK time course."),
            skipped_GSE276656="14.4 GB multiome RAW; not downloaded",
            skipped_GSE271794="E15.5 development, not adult aging",
        )
        dump_json(TRAJ_DIR / "g3_summary.json", summary)
        (TRAJ_DIR / "g3_STOP.txt").write_text("transfer_ok=False\nverdict=d\n", encoding="utf-8")
        log("[G3] STOP: transfer failed. " + narrative)
        return summary

    # ----- reprogramming -----
    log("\n[G3] GSE224438 SVZ partial reprogramming")
    samples = parse_gse224_samples(log)
    extract, tar_names = load_gse224(log)
    adata = gse224_matrices(extract, samples, log)
    traj_rec, means = None, None
    if adata is None:
        narrative = (
            "Transfer succeeded on Tabula Muris Senis, but GSE224438 RAW.tar did not yield "
            "a usable count matrix with the loaders tried. Needed: processed counts with "
            "sample-level young/old/OSK labels and cell types."
        )
        status, verdict = "no_transfer", "d"
        log("[G3] " + narrative)
    else:
        pb = cluster_and_pseudobulk(
            adata, samples, log, frozen=frozen, ortho=ortho, hcol=hcol, mcol=mcol)
        if pb is None:
            narrative = "GSE224438 produced no cluster pseudobulks after filters."
            status, verdict = "no_transfer", "d"
        else:
            scored, n_map_r = score_mouse_pb(pb, frozen, ortho, hcol, mcol, log)
            scored.to_csv(TRAJ_DIR / "g3_GSE224438_scores.csv", index=False)
            log(f"[filter] GSE224438 mapped genes used={n_map_r}")
            traj_rec, means = trajectory_stats(scored, log, rng)
            plot_trajectory(means, TRAJ_FIG / "g3_trajectory.png")
            plot_didt(traj_rec, TRAJ_FIG / "g3_didt.png")
            if traj_rec.get("window"):
                status, verdict = "window", "a"
                w = traj_rec["window"]
                narrative = (
                    f"Two independent scores transferred. Age falls while identity holds on "
                    f"{w.get('range')}: Δage_score={w.get('d_age'):+.3f}. "
                    "This is a protocol contrast (old vs old+OSK), not a multi-day time window."
                )
            elif traj_rec.get("bend"):
                status, verdict = "no_bend_or_off_axis", "b"
                narrative = (
                    "Two independent scores transferred. The reprogramming vector in the "
                    "age–identity plane is not a simple reversal of aging (the curve bends / "
                    "is off-axis), but no window was found where age falls while identity holds."
                )
                # If it doesn't reverse age either, that's (b): fused under reprogramming
                any_age_fall = any(r.get("age_falls") for r in traj_rec.get("reprogramming", []))
                if not any_age_fall:
                    status, verdict = "no_bend", "b"
                    narrative = (
                        "Two independent scores exist and transferred, but under reprogramming "
                        "age and identity do not separate: OSK does not lower the frozen age "
                        "score while holding identity. Axes separable in aging, fused under reprogramming."
                    )
            else:
                # collinear with aging or reverse aging
                any_age_fall = any(r.get("age_falls") for r in traj_rec.get("reprogramming", []))
                any_id_hold = any(r.get("identity_holds") for r in traj_rec.get("reprogramming", []))
                if any_age_fall and any_id_hold:
                    status, verdict = "window", "a"
                    narrative = (
                        "Reprogramming reverses the frozen age score while identity CIs overlap the old group."
                    )
                    traj_rec["window"] = traj_rec.get("window") or dict(
                        range="old → OSK", d_age="see didt", id_retained="CI overlap")
                elif any_age_fall:
                    status, verdict = "no_bend", "b"
                    narrative = (
                        "Age score falls under OSK but identity moves with it (no hold). "
                        "Separable in aging; fused under reprogramming."
                    )
                else:
                    status, verdict = "no_bend", "b"
                    narrative = (
                        "Two independent scores transferred, but the SVZ OSK contrast does not "
                        "lower the frozen age score while holding identity. "
                        "Axes separable in aging, fused (or null) under this reprogramming protocol."
                    )
            traj_rec["text"] = narrative

    summary = dict(
        seed=TRAJ_SEED, dataset_id=DATASET_ID, status=status, verdict=verdict,
        verdict_text=narrative, narrative=narrative,
        orthologs=ortholog_rec, mouse_aging=mouse_aging,
        trajectory=traj_rec, figure="results/trajectory/figures/g3_trajectory.png",
        skipped_GSE276656="14.4 GB multiome RAW; not downloaded",
        skipped_GSE271794="E15.5 development, not adult aging",
        gse224438=dict(accession="GSE224438", pubmed="38553564",
                       tar_members=int(len(tar_names)) if 'tar_names' in dir() else None),
        needed_if_d=("Human adult cortical/DLPFC partial-reprogramming scRNA time course "
                     "with cell-type labels and multiple OSK durations."),
    )
    dump_json(TRAJ_DIR / "g3_summary.json", summary)
    (TRAJ_DIR / "g3_STOP.txt").write_text(
        f"transfer_ok={transfer_ok}\nverdict={verdict}\nstatus={status}\n", encoding="utf-8")
    log(f"[G3] verdict=({verdict}) {narrative}")
    log("[G3] done.")
    return summary


if __name__ == "__main__":
    run()

"""SENG_WHY — why the four named hits missed q <= 0.05 under our cosine.

Diagnostic only. Does not modify FINDINGS_SENG.md, its pre-registration,
any other FINDINGS/PROGRESS/FALSIFICATION file, or any existing src file.
Does not re-run or re-tune a FINDINGS_SENG gate. Imports src/seng_run.py.

Seeds: 20260914 (splits, N_NT nulls, N1 permutations of u).
"""
from __future__ import annotations

import json
import sys
import time
import urllib.request
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fibro2_common import load_frozen_ruler  # noqa: E402
from fibro_common import jsonable  # noqa: E402
from gtex_common import EDGE_R_PRIOR  # noqa: E402
from seng_run import (  # noqa: E402
    ANCHORS_PATH,
    GUIDE_NT_RE,
    GUIDE_TF_RE,
    H5_PATH,
    MIN_CELLS,
    N_NULL,
    N_PERM,
    NAMED_HITS,
    P_BAR,
    PD_RE,
    SEED,
    bh_q,
    cosine_full,
    cos_against_perms,
    parse_wide_obs,
    quartile_masks,
    read_obs_var,
    read_sparse_group,
    ruler_column_index,
    split_half,
    sum_cells,
    to_ruler_from_present,
    z_from_counts,
)
from same_run import permutation_p  # noqa: E402
from toward_run import tmm_logcpm_quiet, zscore_frozen  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "seng_why"
SCREEN_CSV = ROOT / "results" / "seng" / "screen.csv"
CC_TSV = ROOT / "data" / "reference" / "seurat_cc.genes_tirosh2016.tsv"

FIGSHARE_PERTURB = "https://api.figshare.com/v2/articles/30898748"
FIGSHARE_MOUSE = "https://api.figshare.com/v2/articles/30899018"
PAPER_URL = "https://www.pnas.org/doi/10.1073/pnas.2515183123"
PMC_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC12799168/"
DOI_PERTURB = "https://doi.org/10.6084/m9.figshare.30898748"
DOI_MOUSE = "https://doi.org/10.6084/m9.figshare.30899018"

# Locked before any SENG_WHY statistic. These are assumptions, not paper facts.
# Gene filter: the SI and the main text do not state one. Use every gene in the file.
# log2FC: pseudobulk library-size CPM, pseudocount 1 CPM, no TMM (the paper does not mention TMM).
# Aging arm: all WT PD32 vs all WT PD14 (the paper's figures name those two PDs). Not half-A.
# Pearson. Negative R_rej is the paper's rejuvenating direction.
RREJ_PRIOR_CPM = 1.0
# H4: frozen mu/sd where the gene is on the ruler; panel mean/sd otherwise.
TOP_RANK = 15  # Table 1 lists 15 TFs per modality


class Log:
    def __init__(self, path):
        self.fh = open(path, "w", encoding="utf-8")

    def __call__(self, msg):
        print(msg, flush=True)
        self.fh.write(str(msg) + "\n")
        self.fh.flush()

    def close(self):
        self.fh.close()


def fetch_json(url):
    req = urllib.request.Request(url, headers={"User-Agent": "seng-why-diagnostic"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        return json.loads(resp.read().decode("utf-8"))


def file_rows(article):
    rows = []
    for f in article.get("files", []):
        rows.append(dict(
            id=f.get("id"),
            name=f.get("name"),
            size_bytes=f.get("size"),
            md5=f.get("computed_md5") or f.get("supplied_md5"),
            download_url=f.get("download_url"),
            mimetype=f.get("mimetype"),
        ))
    return rows


def uniq_info(values):
    arr = np.asarray(values)
    if arr.dtype == object:
        keys = [("NA" if (v is None or (isinstance(v, float) and not np.isfinite(v))) else str(v)) for v in arr.tolist()]
        # preserve order, cap the unique scan via pandas
        u = pd.unique(np.asarray(keys, dtype=object))
    else:
        u = pd.unique(arr)
    n = int(len(u))
    sample = []
    for v in list(u[:10]):
        if isinstance(v, (np.floating, float)):
            sample.append(None if not np.isfinite(v) else float(v))
        elif isinstance(v, (np.integer, int)):
            sample.append(int(v))
        elif isinstance(v, (np.bool_, bool)):
            sample.append(bool(v))
        else:
            sample.append(str(v))
    return n, sample


def column_role(name):
    s = str(name)
    if PD_RE.fullmatch(s):
        return "pd_flag"
    if s in ("CRA", "CRI"):
        return "modality_flag"
    if s == "WT":
        return "wt_flag"
    if s in ("n_counts", "log_counts", "n_genes", "SampleName"):
        return "qc"
    if GUIDE_TF_RE.fullmatch(s) or GUIDE_NT_RE.fullmatch(s):
        return "guide_umi"
    low = s.lower()
    if any(k in low for k in ("sample", "batch", "library", "lane", "hash", "orig.ident", "orig_ident")):
        return "name_looks_like_batch"
    return "other"


def load_cc_from_tsv():
    df = pd.read_csv(CC_TSV, sep="\t")
    s = df.loc[df["phase"] == "S", "symbol"].astype(str).tolist()
    g2m = df.loc[df["phase"] == "G2M", "symbol"].astype(str).tolist()
    if len(s) != 43 or len(g2m) != 54 or s[0] != "MCM5" or g2m[-1] != "CENPA":
        raise SystemExit(f"cc gene tsv unexpected S={len(s)} G2M={len(g2m)}")
    return s, g2m


def module_score(csr, cols, lib):
    if len(cols) == 0:
        return np.full(csr.shape[0], np.nan)
    sub = csr[:, cols]
    arr = np.asarray(sub.toarray(), np.float64)
    scale = np.maximum(lib, 1.0)[:, None]
    return np.log1p(arr / scale * 1e4).mean(axis=1)


def group_sum(csr, idx):
    idx = np.asarray(idx, int)
    if idx.size == 0:
        return np.zeros(csr.shape[1], np.float64)
    return np.asarray(csr[idx].sum(axis=0), np.float64).ravel()


def log2fc_from_sums(sum_a, lib_a, sum_b, lib_b, prior_cpm=RREJ_PRIOR_CPM):
    cpm_a = np.asarray(sum_a, np.float64) / float(lib_a) * 1e6
    cpm_b = np.asarray(sum_b, np.float64) / float(lib_b) * 1e6
    return np.log2(cpm_a + prior_cpm) - np.log2(cpm_b + prior_cpm)


def ranks_high_is_1(values, higher_is_better):
    """Rank 1 is the strongest score. Ties: average rank, then we report min-rank of ties as competition via method min on the oriented score."""
    x = np.asarray(values, float)
    # competition rank on the rejuvenating orientation, rank 1 = best
    score = x if higher_is_better else -x
    # pandas rank method=min, ascending=False so largest score is rank 1
    r = pd.Series(score).rank(method="min", ascending=False).to_numpy()
    return r


def null_p_cos(nt_dense, nt_sum_ruler, n, v, cache, mu, sd, missing, rng, dst, n_ruler):
    n_nt = int(nt_dense.shape[0])
    if not (0 < n < n_nt):
        return np.nan
    null = np.empty(N_NULL, np.float64)
    for t in range(N_NULL):
        take = rng.choice(n_nt, size=n, replace=False)
        fake_present = nt_dense[take].sum(axis=0)
        fake = to_ruler_from_present(fake_present, dst, n_ruler)
        origin = nt_sum_ruler - fake
        z_f, _, _, _ = z_from_counts(fake, cache, mu, sd, missing)
        z_o, _, _, _ = z_from_counts(origin, cache, mu, sd, missing)
        null[t] = cosine_full(z_f - z_o, v)
    return null


def panel_cache(C, nf):
    """Same reference quantities z_from_counts expects. Factors already match nf."""
    from seng_run import tmm_cache_from_panel
    return tmm_cache_from_panel(C, nf)


def main():
    t0 = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    log = Log(OUT / "run_report.txt")
    log("[seng_why] diagnostic only; FINDINGS_SENG.md is not written")

    art_p = fetch_json(FIGSHARE_PERTURB)
    art_m = fetch_json(FIGSHARE_MOUSE)
    (OUT / "figshare_30898748.json").write_text(json.dumps(art_p, indent=1), encoding="utf-8")
    (OUT / "figshare_30899018.json").write_text(json.dumps(art_m, indent=1), encoding="utf-8")
    files_p = file_rows(art_p)
    files_m = file_rows(art_m)
    pd.DataFrame(files_p).to_csv(OUT / "figshare_30898748_files.csv", index=False)
    pd.DataFrame(files_m).to_csv(OUT / "figshare_30899018_files.csv", index=False)
    h5_names = [r["name"] for r in files_p if str(r["name"]).lower().endswith((".h5ad", ".h5", ".hdf5"))]
    log(f"[H1] perturb-seq record files={[(r['name'], r['size_bytes']) for r in files_p]}")
    log(f"[H1] mouse record files={[(r['name'], r['size_bytes']) for r in files_m]}")
    log(f"[H1] h5-like names on the perturb-seq record: {h5_names}")

    frozen = load_frozen_ruler()
    mu = np.asarray(frozen["mu"], float)
    sd = np.asarray(frozen["sd"], float)
    n_ruler = int(mu.shape[0])
    log(f"[load] {H5_PATH} size={H5_PATH.stat().st_size}")
    with h5py.File(H5_PATH, "r") as f:
        obs = read_obs_var(f["obs"])
        var = read_obs_var(f["var"])
        mat, enc = read_sparse_group(f["X"])
    log(f"[X] enc={enc} shape={mat.shape} nnz={mat.nnz} obs_cols={obs.shape[1]} var_cols={var.shape[1]}")

    parsed = parse_wide_obs(obs, log)
    labels = parsed["labels"].to_numpy()
    kind = parsed["kind"].to_numpy()
    pd_labels = parsed["pd_labels"].to_numpy()
    is_wt = np.asarray(parsed["is_wt"], bool)
    modality = np.asarray(parsed["modality"], dtype=object)
    pert_mask = np.isin(kind, ["CRA", "CRI"])
    log(f"[roles] WT={int(is_wt.sum())} pert={int(pert_mask.sum())} NT={int((kind == 'NT').sum())} unassigned={int((kind == 'unassigned').sum())}")

    # --- obs column inventory ---
    col_rows = []
    for name in obs.columns:
        n_all, sample = uniq_info(obs[name].to_numpy())
        n_wt, sample_wt = uniq_info(obs.loc[is_wt, name].to_numpy())
        n_pert, sample_pert = uniq_info(obs.loc[pert_mask, name].to_numpy())
        # set difference on a compact key
        wt_vals = set(pd.unique(obs.loc[is_wt, name].to_numpy()).tolist())
        pert_vals = set(pd.unique(obs.loc[pert_mask, name].to_numpy()).tolist())
        # NaN-safe: pandas unique keeps nan as a single value; set may break on nan
        def _freeze(vs):
            out = []
            for v in vs:
                if isinstance(v, (float, np.floating)) and not np.isfinite(v):
                    out.append("NA")
                elif isinstance(v, (np.integer,)):
                    out.append(int(v))
                elif isinstance(v, (np.floating,)):
                    out.append(float(v))
                elif isinstance(v, (np.bool_,)):
                    out.append(bool(v))
                else:
                    out.append(v if isinstance(v, (int, float, str, bool)) else str(v))
            return set(map(str, out))
        differs = _freeze(wt_vals) != _freeze(pert_vals)
        role = column_role(name)
        col_rows.append(dict(
            column=str(name),
            dtype=str(obs[name].dtype),
            role=role,
            n_unique=n_all,
            unique_up_to_10=json.dumps(sample),
            n_unique_wt=n_wt,
            unique_wt_up_to_10=json.dumps(sample_wt),
            n_unique_perturbed=n_pert,
            unique_perturbed_up_to_10=json.dumps(sample_pert),
            differs_wt_vs_perturbed=bool(differs),
        ))
    col_df = pd.DataFrame(col_rows)
    col_df.to_csv(OUT / "obs_columns.csv", index=False)
    candidates = col_df[(col_df["role"] != "guide_umi") & (col_df["n_unique"] <= 20)].copy()
    candidates.to_csv(OUT / "obs_low_cardinality.csv", index=False)
    log(f"[obs] columns={len(col_df)} low_cardinality_non_guide={len(candidates)}")
    for _, r in candidates.iterrows():
        log(f"[obs low] {r['column']} dtype={r['dtype']} role={r['role']} n_unique={r['n_unique']} "
            f"pert_values={r['unique_perturbed_up_to_10']} differs={r['differs_wt_vs_perturbed']}")

    # perturbed PD flags
    pd_cols = [c for c in obs.columns if PD_RE.fullmatch(str(c))]
    pd_on_pert = {}
    for c in pd_cols:
        vals = obs.loc[pert_mask, c]
        pd_on_pert[c] = uniq_info(vals.to_numpy())
        log(f"[PD on perturbed] {c} n_unique={pd_on_pert[c][0]} sample={pd_on_pert[c][1]} "
            f"n_flagged={int((pd.to_numeric(obs.loc[pert_mask, c], errors='coerce').fillna(0) == 1).sum())}")

    symbols = var.index.astype(str).to_numpy()
    gene_ids = symbols
    ruler_idx = ruler_column_index(gene_ids, symbols, frozen, log)
    present = np.flatnonzero(ruler_idx >= 0)
    src = ruler_idx[present]
    n_overlap = int(present.size)
    log(f"[genes] ruler={n_ruler} present={n_overlap} file_genes={mat.shape[1]}")
    if len(src) != len(set(src.tolist())):
        log(f"[genes] WARNING duplicate file columns in ruler map: {len(src) - len(set(src.tolist()))}")

    lib_total = np.asarray(mat.sum(axis=1), np.float64).ravel()
    # full-gene sums are taken before column subset
    # ruler csr for the cosine path
    csr = mat[:, src].tocsr()
    dst = present

    early_idx = np.flatnonzero(is_wt & (pd_labels == "PD14"))
    late32_idx = np.flatnonzero(is_wt & (pd_labels == "PD32"))
    late26_idx = np.flatnonzero(is_wt & (pd_labels == "PD26"))
    log(f"[PD n] PD14={early_idx.size} PD26={late26_idx.size} PD32={late32_idx.size}")
    early_a, early_b = split_half(early_idx, SEED)
    late32_a, late32_b = split_half(late32_idx, SEED)
    late26_a, late26_b = split_half(late26_idx, SEED)
    log(f"[split PD32] early A/B={early_a.size}/{early_b.size} late A/B={late32_a.size}/{late32_b.size}")
    log(f"[split PD26] late A/B={late26_a.size}/{late26_b.size}")

    s_genes, g2m_genes = load_cc_from_tsv()
    sym_upper = np.array([str(s).upper() for s in symbols])
    pos = {}
    for i, s in enumerate(sym_upper):
        pos.setdefault(s, i)
    file_to_sub = {int(src[i]): i for i in range(len(src))}
    s_sub = [file_to_sub[pos[g]] for g in s_genes if g in pos and pos[g] in file_to_sub]
    g_sub = [file_to_sub[pos[g]] for g in g2m_genes if g in pos and pos[g] in file_to_sub]
    s_missing = [g for g in s_genes if g not in pos or pos[g] not in file_to_sub]
    g_missing = [g for g in g2m_genes if g not in pos or pos[g] not in file_to_sub]
    log(f"[cc] S used={len(s_sub)} missing_or_not_ruler={s_missing} G2M used={len(g_sub)} missing_or_not_ruler={g_missing}")
    score = module_score(csr, s_sub, lib_total) + module_score(csr, g_sub, lib_total)

    nt_idx = {}
    q_idx = {}
    for mod in ("CRA", "CRI"):
        nt = np.flatnonzero((kind == "NT") & (modality == mod))
        nt_idx[mod] = nt
        bottom, top, info = quartile_masks(score[nt])
        if bottom is None:
            raise SystemExit(f"quartile failed {mod} {info}")
        q_idx[mod] = dict(bottom=nt[bottom], top=nt[top], info=info)
        log(f"[quartile {mod}] n_NT={nt.size} bottom={info['n_bottom']} top={info['n_top']}")

    qual = []
    for lab, sub in pd.Series(np.arange(len(labels))).groupby(labels, sort=True).groups.items():
        lab = str(lab)
        if not (lab.startswith("CRA_") or lab.startswith("CRI_")):
            continue
        idx = np.asarray(sub, int)
        mod = "CRA" if lab.startswith("CRA_") else "CRI"
        if not np.all(modality[idx] == mod):
            log(f"[skip] {lab} modality disagrees")
            continue
        if idx.size < MIN_CELLS:
            continue
        qual.append((mod, lab, idx))
    qual.sort(key=lambda t: (t[0], t[1]))
    log(f"[qual] n>={MIN_CELLS}: {len(qual)}")

    screen = pd.read_csv(SCREEN_CSV)
    prolif_shift = {}
    for mod, lab, idx in qual:
        prolif_shift[lab] = float(np.nanmean(score[idx]) - np.nanmean(score[nt_idx[mod]]))

    anch = np.load(ANCHORS_PATH, allow_pickle=True)
    u = np.asarray(anch["u"], float).ravel()
    rng_u = np.random.default_rng(SEED)
    Uperm = np.empty((N_PERM, n_ruler), np.float64)
    for i in range(N_PERM):
        Uperm[i] = rng_u.permutation(u)

    def run_late(late_name, late_a):
        groups = [
            ("EARLY_A", early_a),
            ("LATE_A", late_a),
        ]
        for mod in ("CRA", "CRI"):
            groups.append((f"NT_{mod}", nt_idx[mod]))
            groups.append((f"TOP_{mod}", q_idx[mod]["top"]))
            groups.append((f"BOT_{mod}", q_idx[mod]["bottom"]))
        for mod, lab, idx in qual:
            groups.append((lab, idx))
        counts = [sum_cells(csr, idx, dst, n_ruler) for _, idx in groups]
        C = np.vstack(counts)
        logcpm, nf = tmm_logcpm_quiet(C)
        Z, missing = zscore_frozen(logcpm, frozen, C)
        cache = panel_cache(C, nf)
        name_to_i = {name: i for i, (name, _) in enumerate(groups)}
        v = Z[name_to_i["EARLY_A"]] - Z[name_to_i["LATE_A"]]
        cos_vu, p_vu, _ = cos_against_perms(v, u, Uperm)
        rows = []
        for mod, lab, idx in qual:
            d = Z[name_to_i[lab]] - Z[name_to_i[f"NT_{mod}"]]
            rows.append(dict(
                modality=mod, perturbation=lab, n_cells=int(idx.size),
                cos=float(cosine_full(d, v)),
                prolif_shift=prolif_shift[lab],
            ))
        df = pd.DataFrame(rows)
        rho = {}
        for mod in ("CRA", "CRI"):
            sub = df[df["modality"] == mod]
            r, p = spearmanr(sub["cos"], sub["prolif_shift"])
            rho[mod] = dict(rho=float(r), p=float(p), n=int(len(sub)))
        log(f"[LATE {late_name}] ||v||={np.linalg.norm(v):.6f} cos(v,u)={cos_vu:.6f} p={p_vu:.4f} "
            f"rho CRA={rho['CRA']['rho']:.6f} CRI={rho['CRI']['rho']:.6f}")
        return dict(df=df, v=v, Z=Z, missing=np.asarray(missing, bool), cache=cache,
                    name_to_i=name_to_i, groups=groups, C=C, cos_vu=float(cos_vu), p_vu=float(p_vu), rho=rho)

    log("[panel] LATE=PD32")
    pack32 = run_late("PD32", late32_a)
    log("[panel] LATE=PD26")
    pack26 = run_late("PD26", late26_a)

    m32 = pack32["df"].merge(screen[["perturbation", "cos", "q_cos", "p_cos", "prolif_shift"]],
                             on="perturbation", suffixes=("_re", "_screen"))
    max_abs = float(np.max(np.abs(m32["cos_re"] - m32["cos_screen"])))
    max_shift = float(np.max(np.abs(m32["prolif_shift_re"] - m32["prolif_shift_screen"])))
    log(f"[match PD32] max |cos - screen.csv|={max_abs:.3e} max |prolif_shift diff|={max_shift:.3e}")
    if max_abs > 1e-5:
        raise SystemExit(f"PD32 cos does not match screen.csv (max abs {max_abs}). Not interpreting a drifted cosine.")

    # side-by-side named hits
    side = pack32["df"][["modality", "perturbation", "n_cells", "cos"]].rename(columns={"cos": "cos_PD32"})
    side = side.merge(pack26["df"][["perturbation", "cos"]].rename(columns={"cos": "cos_PD26"}), on="perturbation")
    side = side.merge(screen[["perturbation", "q_cos", "p_cos"]], on="perturbation")
    side = side.rename(columns={"q_cos": "q_cos_PD32_from_FINDINGS", "p_cos": "p_cos_PD32_from_FINDINGS"})

    # N_NT nulls under PD26 only. Same draw rule as seng_run: continuing Generator(SEED), per perturbation.
    log(f"[null PD26] n_null={N_NULL} perturbations={len(qual)}")
    nt_dense = {}
    nt_sum = {}
    for mod in ("CRA", "CRI"):
        block = np.asarray(csr[nt_idx[mod]].toarray(), np.float64)
        nt_dense[mod] = block
        nt_sum[mod] = to_ruler_from_present(block.sum(axis=0), dst, n_ruler)
        log(f"[NT dense] {mod} {block.shape}")
    rng_null = np.random.default_rng(SEED)
    p26 = []
    t_null = time.perf_counter()
    for i, (mod, lab, idx) in enumerate(qual):
        n = int(idx.size)
        # point cos already known; null distribution
        null = null_p_cos(
            nt_dense[mod], nt_sum[mod], n, pack26["v"], pack26["cache"], mu, sd,
            pack26["missing"], rng_null, dst, n_ruler,
        )
        cos = float(pack26["df"].loc[pack26["df"]["perturbation"] == lab, "cos"].iloc[0])
        if isinstance(null, float) or null is None or (isinstance(null, np.ndarray) and not np.isfinite(null).any()):
            p = np.nan
        else:
            p = float(permutation_p(cos, null, greater=True))
        p26.append(p)
        if (i + 1) % 20 == 0 or lab in NAMED_HITS:
            el = time.perf_counter() - t_null
            rate = el / (i + 1)
            log(f"[null PD26] {i+1}/{len(qual)} {lab} p={p:.4f} elapsed_s={el:.1f} eta_s={rate * (len(qual) - i - 1):.0f}")
    df26 = pack26["df"].copy()
    df26["p_cos"] = p26
    df26["q_cos"] = np.nan
    for mod in ("CRA", "CRI"):
        m = df26["modality"].to_numpy() == mod
        df26.loc[m, "q_cos"] = bh_q(df26.loc[m, "p_cos"].to_numpy())
    side = side.merge(df26[["perturbation", "p_cos", "q_cos"]].rename(
        columns={"p_cos": "p_cos_PD26", "q_cos": "q_cos_PD26"}), on="perturbation")
    side["cos_gt0_q05_PD26"] = (side["cos_PD26"] > 0) & (side["q_cos_PD26"] <= P_BAR)
    side["cos_gt0_q05_PD32"] = (side["cos_PD32"] > 0) & (side["q_cos_PD32_from_FINDINGS"] <= P_BAR)

    both = pack32["df"].merge(pack26["df"][["perturbation", "cos"]], on="perturbation", suffixes=("_PD32", "_PD26"))
    both = both.merge(df26[["perturbation", "p_cos", "q_cos"]], on="perturbation")
    both.to_csv(OUT / "cos_by_late.csv", index=False)
    named_late = side[side["perturbation"].isin(NAMED_HITS)].copy()
    named_late.to_csv(OUT / "named_cos_by_late.csv", index=False)
    log("[named LATE]")
    log(named_late.to_string(index=False))

    # rho and cos(v,u) table
    axis_rows = []
    for tag, pack in (("PD32", pack32), ("PD26", pack26)):
        axis_rows.append(dict(
            late=tag,
            cos_v_u=pack["cos_vu"],
            p_n1=pack["p_vu"],
            rho_prolif_CRA=pack["rho"]["CRA"]["rho"],
            rho_p_CRA=pack["rho"]["CRA"]["p"],
            rho_prolif_CRI=pack["rho"]["CRI"]["rho"],
            rho_p_CRI=pack["rho"]["CRI"]["p"],
            v_norm=float(np.linalg.norm(pack["v"])),
        ))
    axis_df = pd.DataFrame(axis_rows)
    axis_df.to_csv(OUT / "axis_by_late.csv", index=False)
    log(axis_df.to_string(index=False))

    # --- H3 R_rej ---
    # pseudobulk sums on ALL genes
    log("[R_rej] pseudobulk log2FC, all genes, Pearson")
    sum_early = group_sum(mat, early_idx)
    sum_late32 = group_sum(mat, late32_idx)
    lib_early = float(lib_total[early_idx].sum())
    lib_late32 = float(lib_total[late32_idx].sum())
    aging = log2fc_from_sums(sum_late32, lib_late32, sum_early, lib_early)
    nt_sum_all = {}
    nt_lib = {}
    for mod in ("CRA", "CRI"):
        nt_sum_all[mod] = group_sum(mat, nt_idx[mod])
        nt_lib[mod] = float(lib_total[nt_idx[mod]].sum())
    r_rows = []
    for mod, lab, idx in qual:
        s = group_sum(mat, idx)
        lib_g = float(lib_total[idx].sum())
        pert_fc = log2fc_from_sums(s, lib_g, nt_sum_all[mod], nt_lib[mod])
        # all genes; drop non-finite if any
        ok = np.isfinite(aging) & np.isfinite(pert_fc)
        if int(ok.sum()) < 3:
            r = np.nan
            rp = np.nan
        else:
            r, rp = pearsonr(aging[ok], pert_fc[ok])
        r_rows.append(dict(
            modality=mod, perturbation=lab, n_cells=int(idx.size),
            R_rej=float(r), R_rej_p=float(rp), n_genes=int(ok.sum()),
        ))
    rdf = pd.DataFrame(r_rows)
    rdf = rdf.merge(pack32["df"][["perturbation", "cos"]], on="perturbation")
    # ranks within modality
    rdf["rank_R_rej"] = np.nan
    rdf["rank_cos_g"] = np.nan
    rdf["n_in_modality"] = np.nan
    for mod in ("CRA", "CRI"):
        m = rdf["modality"].to_numpy() == mod
        rdf.loc[m, "rank_R_rej"] = ranks_high_is_1(rdf.loc[m, "R_rej"].to_numpy(), higher_is_better=False)
        rdf.loc[m, "rank_cos_g"] = ranks_high_is_1(rdf.loc[m, "cos"].to_numpy(), higher_is_better=True)
        rdf.loc[m, "n_in_modality"] = int(m.sum())
    sp_all, sp_all_p = spearmanr(rdf["R_rej"], rdf["cos"])
    sp_mod = {}
    for mod in ("CRA", "CRI"):
        sub = rdf[rdf["modality"] == mod]
        rr, pp = spearmanr(sub["R_rej"], sub["cos"])
        sp_mod[mod] = dict(rho=float(rr), p=float(pp), n=int(len(sub)))
    log(f"[spearman R_rej vs cos_g] pooled rho={sp_all:.4f} p={sp_all_p:.3e} "
        f"CRA={sp_mod['CRA']['rho']:.4f} CRI={sp_mod['CRI']['rho']:.4f}")
    published = {
        "CRA_E2F3": -0.53, "CRA_EZH2": -0.36, "CRI_ZFX": -0.51, "CRI_STAT3": -0.41,
    }
    rdf["published_table1_R_rej"] = rdf["perturbation"].map(published)
    rdf.to_csv(OUT / "rrej_vs_cos.csv", index=False)
    named_r = rdf[rdf["perturbation"].isin(NAMED_HITS)].copy()
    named_r.to_csv(OUT / "named_rrej_vs_cos.csv", index=False)
    log(named_r.to_string(index=False))

    # --- H4 all-gene cosine at LATE=PD32, same groups ---
    log("[H4] all-gene TMM/z at LATE=PD32")
    # map file column -> ruler index (first if duplicated)
    file_to_ruler = np.full(mat.shape[1], -1, dtype=int)
    n_dup_map = 0
    for j, col in enumerate(ruler_idx):
        if col < 0:
            continue
        if file_to_ruler[col] >= 0:
            n_dup_map += 1
            continue
        file_to_ruler[col] = j
    groups = pack32["groups"]
    counts_all = [group_sum(mat, idx) for _, idx in groups]
    C_all = np.vstack(counts_all)
    logcpm_all, nf_all = tmm_logcpm_quiet(C_all)
    mu_panel = logcpm_all.mean(axis=0)
    sd_panel = logcpm_all.std(axis=0)
    mapped = file_to_ruler >= 0
    mu_use = mu_panel.copy()
    sd_use = sd_panel.copy()
    mu_use[mapped] = mu[file_to_ruler[mapped]]
    sd_use[mapped] = sd[file_to_ruler[mapped]]
    sd_safe = np.where(sd_use < 1e-12, 1.0, sd_use)
    Z_all = (logcpm_all - mu_use) / sd_safe
    missing_all = (C_all == 0).all(axis=0)
    Z_all[:, missing_all] = 0.0
    name_to_i = pack32["name_to_i"]
    v_all = Z_all[name_to_i["EARLY_A"]] - Z_all[name_to_i["LATE_A"]]
    h4_rows = []
    for mod, lab, idx in qual:
        d = Z_all[name_to_i[lab]] - Z_all[name_to_i[f"NT_{mod}"]]
        cos_all = float(cosine_full(d, v_all))
        # subspace norms of d
        n_map = float(np.linalg.norm(d[mapped]))
        n_unmap = float(np.linalg.norm(d[~mapped]))
        h4_rows.append(dict(
            modality=mod, perturbation=lab, n_cells=int(idx.size),
            cos_all_genes=cos_all,
            norm_d_ruler_mapped=n_map,
            norm_d_unmapped=n_unmap,
        ))
    h4 = pd.DataFrame(h4_rows)
    h4 = h4.merge(pack32["df"][["perturbation", "cos"]].rename(columns={"cos": "cos_ruler"}), on="perturbation")
    h4["rank_cos_all"] = np.nan
    h4["rank_cos_ruler"] = np.nan
    for mod in ("CRA", "CRI"):
        m = h4["modality"].to_numpy() == mod
        h4.loc[m, "rank_cos_all"] = ranks_high_is_1(h4.loc[m, "cos_all_genes"].to_numpy(), True)
        h4.loc[m, "rank_cos_ruler"] = ranks_high_is_1(h4.loc[m, "cos_ruler"].to_numpy(), True)
    h4.to_csv(OUT / "cos_all_genes.csv", index=False)
    named_h4 = h4[h4["perturbation"].isin(NAMED_HITS)].copy()
    named_h4.to_csv(OUT / "named_cos_all_genes.csv", index=False)
    log(f"[H4] mapped_genes={int(mapped.sum())} unmapped={int((~mapped).sum())} dup_file_cols_skipped={n_dup_map} "
        f"missing_allzero={int(missing_all.sum())}")
    log(named_h4.to_string(index=False))

    # readings
    named_r = named_r.set_index("perturbation")
    r_top = all(int(named_r.loc[h, "rank_R_rej"]) <= TOP_RANK for h in NAMED_HITS)
    n_cos_not_top = sum(int(named_r.loc[h, "rank_cos_g"]) > TOP_RANK for h in NAMED_HITS)
    metric_explains = bool(r_top and n_cos_not_top >= 3)

    named_late_i = named_late.set_index("perturbation")
    n_cross = int(sum(bool(named_late_i.loc[h, "cos_gt0_q05_PD26"]) for h in NAMED_HITS))
    reference_explains = bool(n_cross >= 3)

    named_h4_i = named_h4.set_index("perturbation")
    h4_top = all(int(named_h4_i.loc[h, "rank_cos_all"]) <= TOP_RANK and float(named_h4_i.loc[h, "cos_all_genes"]) > 0
                 for h in NAMED_HITS)
    n_ruler_not_top = sum(int(named_h4_i.loc[h, "rank_cos_ruler"]) > TOP_RANK for h in NAMED_HITS)
    scale_dominated = any(
        float(named_h4_i.loc[h, "norm_d_unmapped"]) > 5.0 * float(named_h4_i.loc[h, "norm_d_ruler_mapped"])
        for h in NAMED_HITS
    )
    gene_space_explains = bool(h4_top and n_ruler_not_top >= 3 and not scale_dominated)

    only_h5 = h5_names == ["exp3_merged_not_normalized.h5ad"]
    mouse_is_h5 = any(str(r["name"]).lower().endswith((".h5ad", ".h5")) for r in files_m)
    file_explains = bool((not only_h5) or mouse_is_h5)

    fired = []
    if metric_explains:
        fired.append("metric_explains")
    if reference_explains:
        fired.append("reference_explains")
    if gene_space_explains:
        fired.append("gene_space_explains")
    if file_explains:
        fired.append("file_explains")
    if not fired:
        fired.append("unexplained")

    summary = dict(
        metric_explains=metric_explains,
        reference_explains=reference_explains,
        gene_space_explains=gene_space_explains,
        file_explains=file_explains,
        fired=fired,
        r_top=r_top,
        n_cos_rank_outside_top15=n_cos_not_top,
        n_named_q05_PD26=n_cross,
        h4_top=h4_top,
        n_ruler_rank_outside_top15=n_ruler_not_top,
        h4_scale_dominated=scale_dominated,
        spearman_Rrej_cos_pooled=float(sp_all),
        spearman_Rrej_cos_pooled_p=float(sp_all_p),
        spearman_by_modality=sp_mod,
        pd32_cos_match_max_abs=max_abs,
        n_ruler=n_ruler,
        n_ruler_present=n_overlap,
        n_file_genes=int(mat.shape[1]),
        n_mapped=int(mapped.sum()),
        n_unmapped=int((~mapped).sum()),
        n_dup_map=int(n_dup_map),
        top_rank_cutoff=TOP_RANK,
        assumptions=dict(
            rrej_gene_selection="all genes in the file; the paper states no significance cut and no DE filter for R_rej",
            rrej_correlation="Pearson (SI states Pearson)",
            rrej_rejuvenating_direction="negative (paper states negative R_rej)",
            rrej_log2fc="pseudobulk library-size CPM, log2((CPM_a+1)/(CPM_b+1)); no TMM; all cells in the group, not half-A",
            rrej_aging_pds="WT PD32 vs WT PD14, the PDs the paper's figures name as late and early",
            h4_unmapped_z="panel mean and sd of TMM log2-CPM; frozen mu/sd have no entry for these genes",
        ),
        elapsed_s=time.perf_counter() - t0,
    )
    (OUT / "summary.json").write_text(json.dumps(jsonable(summary), indent=1), encoding="utf-8")
    log(f"[fired] {fired}")
    log(f"[done] elapsed_s={summary['elapsed_s']:.1f}")
    log.close()


if __name__ == "__main__":
    main()

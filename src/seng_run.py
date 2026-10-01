"""SENG — do Sengstack TF perturbations move fibroblasts toward young cells?

Reads the figshare h5ad. Reuses src/toward_run.py and src/same_run.py.
Does not modify those files, any other src file, or any other FINDINGS/PROGRESS.

Pre-registration must already exist at results/seng/PREREG.flag.
Southard/Hs27 is not opened. GSE325735 is not opened.

Usage: python src/seng_run.py
"""
from __future__ import annotations

import bz2
import hashlib
import re
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import DATA_RAW, RESULTS, ROOT  # noqa: E402
from brain_phase1_common import Logger, dump_json, load_json, strip_ensembl  # noqa: E402
from fibro_common import FROZEN_RULER, jsonable  # noqa: E402
from fibro2_common import align_counts_to_ruler, load_frozen_ruler  # noqa: E402
from gtex_common import (  # noqa: E402
    EDGE_R_PRIOR,
    TMM_LOGRATIO_TRIM,
    TMM_SUM_TRIM,
    StopStep,
    tmm_norm_factors,
)
from target_common import md_table  # noqa: E402
from toward_run import (  # noqa: E402
    _fmt,
    _fmt_ci,
    _fmt_p,
    percentile_ci,
    ruler_score,
    tmm_logcpm_quiet,
    zscore_frozen,
)
from same_run import (  # noqa: E402
    axis_stats,
    ci_covers_point,
    ci_excludes_zero,
    cosine_full,
    frac_along,
    monotonic_increasing,
    permutation_p,
    resample_idx,
    split_half,
)


SENG_DIR = RESULTS / "seng"
SENG_DIR.mkdir(parents=True, exist_ok=True)
RAW_DIR = DATA_RAW / "seng"
RAW_DIR.mkdir(parents=True, exist_ok=True)
FINDINGS_PATH = ROOT / "FINDINGS_SENG.md"
PROGRESS_PATH = ROOT / "PROGRESS_SENG.md"
PREREG_FLAG = SENG_DIR / "PREREG.flag"
MANIFEST_PATH = SENG_DIR / "manifest.json"
H5_NAME = "exp3_merged_not_normalized.h5ad"
H5_PATH = RAW_DIR / H5_NAME
FIGSHARE_URL = "https://ndownloader.figshare.com/files/60428597"
FIGSHARE_API_SIZE = 1426871504
CC_RDA = ROOT / "data" / "reference" / "cc.genes.rda"
CC_TSV = ROOT / "data" / "reference" / "seurat_cc.genes_tirosh2016.tsv"
ANCHORS_PATH = RESULTS / "toward" / "anchors.npz"

SEED = 20260914
BOOT_SEED = 20260918
N_NULL = 200
N_PERM = 200
N_BOOT = 200
MIN_CELLS = 30
MIN_WT = 200
P_BAR = 0.05
F_VALUES = (0.0, 0.25, 0.50)
RHO_PROLIF_BAR = 0.5
IDENTITY_DROP = 0.5
NAMED_HITS = ("CRA_E2F3", "CRA_EZH2", "CRI_STAT3", "CRI_ZFX")
IDENTITY_GENES = (
    "COL1A1", "COL1A2", "FN1", "LUM", "PDGFRA", "PDGFRB",
    "POSTN", "PRRX1", "SERPINH1", "VIM",
)
GSE325735_BAN = "GSE325735 is out of scope and is not opened."
SOUTHARD_BAN = "Southard/Hs27 is out of scope and is not opened."

CRA_RE = re.compile(r"^CRA_.+")
CRI_RE = re.compile(r"^CRI_.+")
NT_RE = re.compile(r"^NT_.*")
PERT_RE = re.compile(r"^(CRA_.+|CRI_.+|NT_.*)$")
PD_RE = re.compile(r"^PD(\d+)$")
WT_TOKENS = {
    "WT", "WILDTYPE", "WILD-TYPE", "WILD_TYPE",
    "UNTRANSDUCED", "UNINFECTED", "NONTRANSDUCED", "NON-TRANSDUCED",
    "NOTRANSDUCTION", "NO_TRANSDUCTION", "UNTRANSFECTED",
}
MODALITY_TOKENS = {
    "CRA": "CRA", "CRISPRA": "CRA", "CRISPR-A": "CRA", "A": "CRA",
    "CRI": "CRI", "CRISPRI": "CRI", "CRISPR-I": "CRI", "I": "CRI",
}
GUIDE_NAMES = {
    "guide", "guides", "guide_id", "guide_identity", "sgrna", "sgrnas",
    "grna", "feature_call", "protospacer", "guide_name", "target_guide",
}
SYMBOL_COL_PRIORITY = (
    "feature_name", "gene_symbols", "gene_symbol", "gene_name", "symbol",
)
ID_COL_PRIORITY = ("gene_ids", "gene_id", "ensembl", "ensembl_id", "gene")


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
    return dict(seed=SEED, boot_seed=BOOT_SEED, failures=[], status="INIT")


def save_manifest(man):
    dump_json(MANIFEST_PATH, jsonable(man))
    return MANIFEST_PATH


def record_failure(step, message, details=None):
    man = load_manifest()
    rec = dict(step=step, message=str(message), details=jsonable(details or {}))
    man.setdefault("failures", []).append(rec)
    man["status"] = "STOP"
    save_manifest(man)
    return rec


def _s(x):
    if x is None:
        return ""
    if isinstance(x, (bytes, np.bytes_)):
        return x.decode("utf-8", "replace")
    if isinstance(x, str):
        return x
    if isinstance(x, (float, np.floating)) and not np.isfinite(x):
        return ""
    return str(x)


def _attr(obj, key, default=None):
    if not hasattr(obj, "attrs") or key not in obj.attrs:
        return default
    v = obj.attrs[key]
    if isinstance(v, np.ndarray):
        if v.shape == ():
            v = v.item()
        else:
            return [_s(x) for x in v.tolist()]
    if isinstance(v, (bytes, np.bytes_)):
        return v.decode("utf-8", "replace")
    return v


def _to_object_str(arr):
    arr = np.asarray(arr)
    flat = arr.ravel()
    out = np.empty(flat.shape[0], dtype=object)
    for i, x in enumerate(flat):
        if x is None:
            out[i] = ""
        elif isinstance(x, (bytes, np.bytes_)):
            out[i] = x.decode("utf-8", "replace")
        elif isinstance(x, str):
            out[i] = x
        elif isinstance(x, (float, np.floating)) and not np.isfinite(x):
            out[i] = ""
        else:
            out[i] = str(x)
    return out.reshape(arr.shape)


def _read_1d(node):
    if isinstance(node, h5py.Group):
        enc = str(_attr(node, "encoding-type", "") or "")
        if "categorical" in enc or ("categories" in node and "codes" in node):
            cats = _to_object_str(node["categories"][:]).ravel()
            codes = np.asarray(node["codes"][:]).ravel()
            out = np.array([""] * codes.shape[0], dtype=object)
            ok = codes >= 0
            if ok.any():
                out[ok] = cats[codes[ok]]
            return out
        raise StopStep("h5", f"cannot read group {node.name} encoding={enc!r}. Not substituting.")
    arr = node[:]
    if getattr(arr, "dtype", None) is not None and arr.dtype.kind in ("S", "O", "U"):
        return _to_object_str(arr).ravel()
    return np.asarray(arr).ravel()


def read_obs_var(group):
    if group is None or not isinstance(group, h5py.Group):
        raise StopStep("h5", "obs/var group missing. Not substituting.")
    order = _attr(group, "column-order", None)
    keys = [k for k in group.keys() if not str(k).startswith("_")]
    if order:
        names = [c for c in order if c in group]
        for k in keys:
            if k not in names and k not in ("index",):
                names.append(k)
    else:
        names = keys
    index = None
    for key in ("_index", "index"):
        if key in group:
            index = _read_1d(group[key])
            break
    data = {}
    for name in names:
        if name in ("_index", "index"):
            continue
        data[str(name)] = _read_1d(group[name])
    n = None
    for v in data.values():
        n = len(v)
        break
    if index is None:
        index = np.arange(n if n is not None else 0)
    if n is not None and len(index) != n:
        raise StopStep(
            "h5",
            f"obs/var index length {len(index)} != column length {n}. Not repairing rows.",
        )
    df = pd.DataFrame(data)
    df.index = pd.Index(_to_object_str(np.asarray(index)).ravel(), name="index")
    return df


def dump_tree(f, path, max_depth=3):
    lines = []

    def walk(g, depth):
        if depth > max_depth:
            return
        for k in g.keys():
            obj = g[k]
            kind = "group" if isinstance(obj, h5py.Group) else "dataset"
            extra = ""
            if isinstance(obj, h5py.Dataset):
                extra = f" shape={obj.shape} dtype={obj.dtype}"
            enc = _attr(obj, "encoding-type", "")
            lines.append(f"{'  ' * depth}{k}: {kind}{extra} encoding={enc}")
            if isinstance(obj, h5py.Group):
                walk(obj, depth + 1)

    walk(f, 0)
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return lines


def scan_numeric(ds, chunk=5_000_000):
    """Return raw-count audit of a 1-d values dataset or a dense matrix flattened in chunks."""
    dtype = ds.dtype
    integer_dtype = np.issubdtype(dtype, np.integer)
    n_frac = 0
    n_neg = 0
    n_nan = 0
    n_seen = 0
    vmin = None
    vmax = None
    if ds.ndim == 1:
        n = int(ds.shape[0])
        spans = [(i, min(i + chunk, n)) for i in range(0, n, chunk)]
        getter = lambda a, b: ds[a:b]
    else:
        # Row chunks. Avoid loading a dense multi-GB matrix at once.
        n_rows = int(ds.shape[0])
        width = int(np.prod(ds.shape[1:])) if ds.ndim > 1 else 1
        row_chunk = max(1, chunk // max(width, 1))
        spans = [(i, min(i + row_chunk, n_rows)) for i in range(0, n_rows, row_chunk)]
        getter = lambda a, b: np.asarray(ds[a:b]).ravel()
    for a, b in spans:
        sl = np.asarray(getter(a, b)).ravel()
        n_seen += int(sl.size)
        if sl.size == 0:
            continue
        if integer_dtype:
            mn = int(sl.min())
            mx = int(sl.max())
            n_neg += int((sl < 0).sum())
        else:
            finite = np.isfinite(sl)
            n_nan += int((~finite).sum())
            slf = sl[finite]
            if slf.size == 0:
                continue
            mn = float(slf.min())
            mx = float(slf.max())
            n_neg += int((slf < -1e-8).sum())
            n_frac += int((np.abs(slf - np.round(slf)) > 1e-6).sum())
        vmin = mn if vmin is None else min(vmin, mn)
        vmax = mx if vmax is None else max(vmax, mx)
    raw_ok = (n_frac == 0) and (n_neg == 0) and (n_nan == 0) and (n_seen > 0)
    return dict(
        dtype=str(dtype), integer_dtype=bool(integer_dtype), n_seen=int(n_seen),
        n_fractional=int(n_frac), n_negative=int(n_neg), n_nonfinite=int(n_nan),
        min=vmin, max=vmax, raw_integer_counts=bool(raw_ok),
    )


def read_sparse_group(g):
    enc = str(_attr(g, "encoding-type", "") or "")
    shape = _attr(g, "shape", None)
    if shape is None:
        raise StopStep("h5", f"{g.name} sparse group has no shape. Not substituting.")
    shape = tuple(int(x) for x in shape)
    data = g["data"][:]
    indices = g["indices"][:]
    indptr = g["indptr"][:]
    if "csr" not in enc and "csc" not in enc:
        if len(indptr) == shape[0] + 1:
            enc = "csr_matrix"
        elif len(indptr) == shape[1] + 1:
            enc = "csc_matrix"
        else:
            raise StopStep(
                "h5",
                f"{g.name} encoding {enc!r} indptr={len(indptr)} shape={shape}. Not substituting.",
            )
    if "csr" in enc:
        mat = sparse.csr_matrix((data, indices, indptr), shape=shape)
    else:
        mat = sparse.csc_matrix((data, indices, indptr), shape=shape).tocsr()
    return mat, enc


def locate_count_matrix(f):
    """Return (matrix csr or None if dense path, audit, source_name).

    X is used when it is raw integer counts. Otherwise a raw layer, if one exists.
    """
    audits = []
    if "X" not in f:
        raise StopStep("h5", "h5ad has no X. Not substituting.")
    X = f["X"]
    if isinstance(X, h5py.Group) and "data" in X:
        audit = scan_numeric(X["data"])
        audit["where"] = "X/data"
        audits.append(audit)
        if audit["raw_integer_counts"]:
            mat, enc = read_sparse_group(X)
            return mat, audit, f"X ({enc})", audits
    elif isinstance(X, h5py.Dataset):
        audit = scan_numeric(X)
        audit["where"] = "X"
        audits.append(audit)
        if audit["raw_integer_counts"]:
            # Dense raw matrix. Caller loads ruler columns; signal with a sentinel.
            return ("DENSE_X", X), audit, "X (dense)", audits
    layer_names = list(f["layers"].keys()) if "layers" in f else []
    raw_group = "raw" in f
    candidates = []
    if "layers" in f:
        for name in layer_names:
            node = f["layers"][name]
            if isinstance(node, h5py.Group) and "data" in node:
                audit = scan_numeric(node["data"])
            elif isinstance(node, h5py.Dataset):
                audit = scan_numeric(node)
            else:
                continue
            audit["where"] = f"layers/{name}"
            audits.append(audit)
            if audit["raw_integer_counts"]:
                candidates.append(name)
    preferred = [n for n in ("counts", "raw_counts", "raw") if n in candidates]
    chosen = preferred[0] if preferred else (candidates[0] if len(candidates) == 1 else None)
    if chosen is not None:
        node = f["layers"][chosen]
        if isinstance(node, h5py.Group):
            mat, enc = read_sparse_group(node)
            return mat, next(a for a in audits if a["where"] == f"layers/{chosen}"), f"layers/{chosen} ({enc})", audits
        return ("DENSE_LAYER", node), next(a for a in audits if a["where"] == f"layers/{chosen}"), f"layers/{chosen} (dense)", audits
    if raw_group and "X" in f["raw"]:
        node = f["raw"]["X"]
        if isinstance(node, h5py.Group) and "data" in node:
            audit = scan_numeric(node["data"])
            audit["where"] = "raw/X"
            audits.append(audit)
            if audit["raw_integer_counts"]:
                mat, enc = read_sparse_group(node)
                return mat, audit, f"raw/X ({enc})", audits
    msg = (
        "X is not raw integer counts and no raw layer exists. "
        f"audits={jsonable(audits)} layers={layer_names} raw_group={raw_group}."
    )
    raise StopStep("no_raw_counts", msg, dict(audits=audits, layers=layer_names))


def file_md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(8 * 1024 * 1024)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def load_cc_genes(log):
    if not CC_RDA.exists():
        raise StopStep("cc_genes", f"missing {CC_RDA}. Seurat cc.genes list was not substituted from memory.")
    raw = bz2.open(CC_RDA, "rb").read()
    strs = [s.decode() for s in re.findall(rb"[\x20-\x7e]{2,}", raw)]
    if "cc.genes" not in strs or "names" not in strs:
        raise StopStep("cc_genes", f"cc.genes.rda strings unexpected: {strs[:8]} ... {strs[-6:]}")
    genes = strs[strs.index("cc.genes") + 1: strs.index("names")]
    if not (
        len(genes) == 97 and genes[0] == "MCM5" and genes[42] == "E2F8"
        and genes[43] == "HMGB2" and genes[-1] == "CENPA"
    ):
        raise StopStep(
            "cc_genes",
            "Parsed cc.genes.rda does not match Seurat cc.genes boundaries "
            f"(n={len(genes)} first={genes[:1]} cut={genes[42:44]} last={genes[-1:]}). Not substituting a hand list.",
        )
    s_genes = genes[:43]
    g2m_genes = genes[43:]
    CC_TSV.parent.mkdir(parents=True, exist_ok=True)
    lines = ["phase\tsymbol\tsource"]
    for g in s_genes:
        lines.append(f"S\t{g}\tSeurat cc.genes / Tirosh 2016 (cc.genes.rda v4.4.0)")
    for g in g2m_genes:
        lines.append(f"G2M\t{g}\tSeurat cc.genes / Tirosh 2016 (cc.genes.rda v4.4.0)")
    text = "\n".join(lines) + "\n"
    CC_TSV.write_text(text, encoding="utf-8")
    (SENG_DIR / "seurat_cc.genes_tirosh2016.tsv").write_text(text, encoding="utf-8")
    log(f"[cc.genes] S={len(s_genes)} G2M={len(g2m_genes)} written {CC_TSV}")
    return s_genes, g2m_genes


def classify_token(val):
    s = _s(val).strip()
    if PERT_RE.fullmatch(s):
        if s.startswith("CRA_"):
            return "CRA", s
        if s.startswith("CRI_"):
            return "CRI", s
        return "NT", s
    m = PD_RE.fullmatch(s)
    if m:
        return "PD", s
    if s.upper() in WT_TOKENS:
        return "WT", s
    return "", s


def column_fractions(series):
    vals = series.map(lambda x: _s(x).strip()).to_numpy()
    n = len(vals)
    n_cra = n_cri = n_nt = n_pd = n_wt = 0
    for v in vals:
        kind, _ = classify_token(v)
        if kind == "CRA":
            n_cra += 1
        elif kind == "CRI":
            n_cri += 1
        elif kind == "NT":
            n_nt += 1
        elif kind == "PD":
            n_pd += 1
        elif kind == "WT":
            n_wt += 1
    return dict(
        n=n, n_cra=n_cra, n_cri=n_cri, n_nt=n_nt, n_pert=n_cra + n_cri + n_nt,
        n_pd=n_pd, n_wt=n_wt, nunique=int(pd.Series(vals).nunique()),
    )


def pick_column(df, key, log):
    rows = []
    for c in df.columns:
        fr = column_fractions(df[c])
        fr["column"] = c
        rows.append(fr)
        log(f"[inventory col] {c} nunique={fr['nunique']} CRA={fr['n_cra']} CRI={fr['n_cri']} "
            f"NT={fr['n_nt']} PD={fr['n_pd']} WT={fr['n_wt']}")
    tab = pd.DataFrame(rows)
    if key == "perturbation":
        hit = tab[tab["n_pert"] > 0].sort_values("n_pert", ascending=False)
    elif key == "pd":
        hit = tab[tab["n_pd"] > 0].sort_values("n_pd", ascending=False)
    elif key == "wt":
        hit = tab[tab["n_wt"] > 0].sort_values("n_wt", ascending=False)
    else:
        raise StopStep("inventory", f"unknown pick key {key}")
    if hit.empty:
        return None, tab
    top = int(hit.iloc[0][key if key != "perturbation" else "n_pert"] if key != "perturbation" else hit.iloc[0]["n_pert"])
    # tie: more than one column shares the max count
    count_col = "n_pert" if key == "perturbation" else ("n_pd" if key == "pd" else "n_wt")
    tied = hit[hit[count_col] == hit.iloc[0][count_col]]
    if len(tied) > 1:
        raise StopStep(
            "inventory",
            f"Ambiguous {key} column. Tied columns={list(tied['column'])} "
            f"counts={list(tied[count_col])}. Not substituting.",
            dict(columns=list(tied["column"])),
        )
    return str(hit.iloc[0]["column"]), tab


def assign_nt_modality(df, pert_col, log):
    """Return a modality series (CRA/CRI/'') aligned to df, plus a note. STOP if NT unassigned."""
    labels = df[pert_col].map(lambda x: _s(x).strip())
    is_nt = labels.map(lambda s: bool(NT_RE.fullmatch(s))).to_numpy()
    is_cra = labels.map(lambda s: bool(CRA_RE.fullmatch(s))).to_numpy()
    is_cri = labels.map(lambda s: bool(CRI_RE.fullmatch(s))).to_numpy()
    modality = np.array([""] * len(df), dtype=object)
    modality[is_cra] = "CRA"
    modality[is_cri] = "CRI"
    note = ""
    if not is_nt.any():
        raise StopStep("nt_unassignable", "No NT_* labels in the perturbation column. Not assigning another control.")

    # 1. A column on NT rows whose values are modality tokens, and both modalities appear.
    for c in df.columns:
        if c == pert_col:
            continue
        vals = df[c].map(lambda x: _s(x).strip().upper()).to_numpy()
        mapped = np.array([MODALITY_TOKENS.get(v, "") for v in vals], dtype=object)
        nt_mapped = mapped[is_nt]
        kinds = set(k for k in nt_mapped if k)
        if kinds == {"CRA", "CRI"} and (nt_mapped == "").sum() == 0:
            modality[is_nt] = mapped[is_nt]
            note = f"NT modality from column {c!r} (values map to CRA/CRI on every NT row)."
            log(f"[NT] {note}")
            return modality, note

    # 2. The NT label itself encodes modality.
    def label_mod(s):
        u = s.upper()
        if "CRISPRA" in u or u.startswith("NT_CRA") or "_CRA" in u or u.endswith("_A"):
            return "CRA"
        if "CRISPRI" in u or u.startswith("NT_CRI") or "_CRI" in u or u.endswith("_I"):
            return "CRI"
        return ""

    mapped_lab = labels.map(label_mod).to_numpy()
    if is_nt.any() and set(mapped_lab[is_nt]) <= {"CRA", "CRI"} and "" not in set(mapped_lab[is_nt]):
        if set(mapped_lab[is_nt]) == {"CRA", "CRI"}:
            modality[is_nt] = mapped_lab[is_nt]
            note = "NT modality encoded in the NT label itself."
            log(f"[NT] {note}")
            return modality, note

    # 3. A batch/sample column: each partition that contains NT contains perturbed cells of only one modality.
    best = None
    for c in df.columns:
        if c == pert_col:
            continue
        vals = df[c].map(lambda x: _s(x).strip()).to_numpy()
        if pd.Series(vals).nunique() < 2 or pd.Series(vals).nunique() > 500:
            continue
        assigned = modality.copy()
        ok = True
        n_assigned = 0
        for level in pd.unique(vals):
            m = vals == level
            if not (m & is_nt).any():
                continue
            pert_mods = set(modality[m & (is_cra | is_cri)].tolist())
            pert_mods.discard("")
            if len(pert_mods) != 1:
                ok = False
                break
            assigned[m & is_nt] = next(iter(pert_mods))
            n_assigned += int((m & is_nt).sum())
        if ok and n_assigned == int(is_nt.sum()) and set(assigned[is_nt]) == {"CRA", "CRI"}:
            best = (c, assigned)
            break
    if best is not None:
        c, assigned = best
        note = (
            f"NT modality from column {c!r}: each level that contains NT cells contains "
            "perturbed cells of exactly one of CRA or CRI."
        )
        log(f"[NT] {note}")
        return assigned, note

    nt_labels = sorted(set(labels[is_nt].tolist()))
    raise StopStep(
        "nt_unassignable",
        "NT cells cannot be assigned to CRISPRa vs CRISPRi. "
        f"NT labels={nt_labels}. No modality column, no modality-encoded NT name, "
        "and no batch column that is pure for one modality. Not guessing.",
        dict(nt_labels=nt_labels),
    )


PD_NAME_RE = re.compile(r"(passage|population.?doubling|^pd$|doubling|pdl)", re.I)


def pd_number(label):
    s = _s(label).strip()
    m = PD_RE.fullmatch(s)
    if m:
        return int(m.group(1))
    if re.fullmatch(r"\d{1,3}", s):
        v = int(s)
        if 0 <= v <= 200:
            return v
    return None


META_OBS = {"n_counts", "log_counts", "n_genes", "SampleName"}
GUIDE_TF_RE = re.compile(r"^(.+)_(CRA|CRI)(\d+)$")
GUIDE_NT_RE = re.compile(r"^NT_(\d+)$")


def is_wide_guide_obs(columns):
    cols = {str(c) for c in columns}
    if not {"CRA", "CRI", "WT"} <= cols:
        return False
    n_pd = sum(1 for c in cols if PD_RE.fullmatch(c))
    n_g = sum(1 for c in cols if GUIDE_TF_RE.fullmatch(c) or GUIDE_NT_RE.fullmatch(c))
    return n_pd >= 2 and n_g >= 10


def _as_flag(series):
    if series.dtype == bool:
        return series.to_numpy(dtype=bool)
    vals = pd.to_numeric(series, errors="coerce").to_numpy()
    return np.isfinite(vals) & (vals == 1)


def parse_wide_obs(obs, log):
    """Guide identity is a UMI column per guide, not a CRA_<TF> label column.

    Call, fixed before any cosine: a non-WT cell is assigned the guide column
    whose count is a unique strict maximum and > 0. Ties are not given a winner.
    File column {TF}_{CRA|CRI}{k} pools to perturbation {CRA|CRI}_{TF}.
    NT_* columns are not line-specific; the cell's modality is the CRA or CRI
    indicator. A TF-guide call that disagrees with that indicator is not assigned.
    """
    pd_cols = sorted([c for c in obs.columns if PD_RE.fullmatch(str(c))], key=lambda c: pd_number(c))
    known = set(META_OBS) | set(pd_cols) | {"CRA", "CRI", "WT"}
    guide_cols = [c for c in obs.columns if GUIDE_TF_RE.fullmatch(str(c)) or GUIDE_NT_RE.fullmatch(str(c))]
    unknown = [c for c in obs.columns if str(c) not in known and c not in guide_cols]
    if unknown:
        raise StopStep(
            "inventory",
            f"Unrecognized obs columns {unknown[:40]}. Not dropping them.",
        )
    pd_flags = np.column_stack([_as_flag(obs[c]) for c in pd_cols])
    n_pd_on = pd_flags.sum(1)
    if int((n_pd_on > 1).sum()):
        raise StopStep(
            "inventory",
            f"{int((n_pd_on > 1).sum())} cells are marked True in more than one PD column. Not picking one.",
        )
    is_wt = _as_flag(obs["WT"])
    any_pd = n_pd_on == 1
    if not np.array_equal(is_wt, any_pd):
        raise StopStep(
            "inventory",
            "WT column and the PD indicator columns do not mark the same cells. "
            f"n_WT={int(is_wt.sum())} n_any_PD={int(any_pd.sum())} "
            f"WT_only={int((is_wt & ~any_pd).sum())} PD_only={int((~is_wt & any_pd).sum())}. "
            "Not forcing them to agree.",
        )
    cra = _as_flag(obs["CRA"])
    cri = _as_flag(obs["CRI"])
    if int((cra & cri).sum()):
        raise StopStep("inventory", "CRA and CRI indicators overlap. Not assigning both.")
    if int((is_wt & (cra | cri)).sum()):
        raise StopStep("inventory", "WT cells also have CRA or CRI set. Not assigning both roles.")
    if int((~is_wt & ~cra & ~cri).sum()):
        raise StopStep(
            "nt_unassignable",
            "Some non-WT cells have neither CRA nor CRI set. Not guessing a modality.",
        )
    cell_mod = np.where(cra, "CRA", np.where(cri, "CRI", "")).astype(object)
    A = np.column_stack([
        pd.to_numeric(obs[c], errors="coerce").fillna(0.0).to_numpy(dtype=np.float64)
        for c in guide_cols
    ])
    mx = A.max(axis=1)
    second = np.partition(A, -2, axis=1)[:, -2]
    top = np.argmax(A, axis=1)
    called = (mx > second) & (mx > 0)
    n_tie = int(((mx == second) & (mx > 0)).sum())
    gmod = []
    pooled = []
    guide_id = []
    for name in guide_cols:
        m = GUIDE_TF_RE.fullmatch(str(name))
        if m:
            tf, mod, k = m.group(1), m.group(2), m.group(3)
            gmod.append(mod)
            pooled.append(f"{mod}_{tf}")
            guide_id.append(str(name))
        else:
            gmod.append("NT")
            pooled.append(str(name))
            guide_id.append(str(name))
    gmod = np.array(gmod, dtype=object)
    pooled = np.array(pooled, dtype=object)
    guide_id = np.array(guide_id, dtype=object)
    top_mod = gmod[top]
    top_pooled = pooled[top]
    top_guide = guide_id[top]
    disagree = called & np.isin(top_mod, ["CRA", "CRI"]) & (top_mod != cell_mod)
    nt_ok = called & (top_mod == "NT") & np.isin(cell_mod, ["CRA", "CRI"])
    tf_ok = called & np.isin(top_mod, ["CRA", "CRI"]) & (top_mod == cell_mod)
    if int((is_wt & called).sum()):
        raise StopStep(
            "inventory",
            f"{int((is_wt & called).sum())} WT cells have a positive guide count. Not assigning both roles.",
        )
    labels = np.array([""] * len(obs), dtype=object)
    guides_called = np.array([""] * len(obs), dtype=object)
    kind = np.array([""] * len(obs), dtype=object)
    modality = np.array([""] * len(obs), dtype=object)
    labels[tf_ok] = top_pooled[tf_ok]
    labels[nt_ok] = top_pooled[nt_ok]
    guides_called[tf_ok | nt_ok] = top_guide[tf_ok | nt_ok]
    kind[is_wt] = "WT"
    kind[tf_ok] = top_mod[tf_ok]
    kind[nt_ok] = "NT"
    kind[(~is_wt) & ~(tf_ok | nt_ok)] = "unassigned"
    modality[tf_ok | nt_ok] = cell_mod[tf_ok | nt_ok]
    pd_labels = np.array([""] * len(obs), dtype=object)
    for j, c in enumerate(pd_cols):
        pd_labels[pd_flags[:, j]] = str(c)
    # Perturbed cells carry no PD flag. That fact is reported by the caller.
    n_unassigned = int((kind == "unassigned").sum())
    log(
        f"[wide obs] PD={pd_cols} n_WT={int(is_wt.sum())} n_CRA={int(cra.sum())} n_CRI={int(cri.sum())} "
        f"guide_cols={len(guide_cols)} unique_argmax_TF={int(tf_ok.sum())} unique_argmax_NT={int(nt_ok.sum())} "
        f"ties={n_tie} modality_disagreements={int(disagree.sum())} unassigned={n_unassigned}"
    )
    nt_note = (
        "NT_1..NT_6 are shared guide-count columns: a unique-maximum call to each NT_* occurs in both "
        "CRA and CRI cells, so the NT label is not a line id. "
        "Modality of an NT cell is the CRA or CRI indicator (those two columns do not overlap and cover every non-WT cell). "
        f"Unique-argmax NT cells: CRA={int((nt_ok & cra).sum())} CRI={int((nt_ok & cri).sum())}."
    )
    wt_note = (
        "WT column is 0/1 and marks the same cells as exactly one of the boolean PD columns "
        + ",".join(pd_cols)
        + f". n_WT={int(is_wt.sum())}."
    )
    gnote = (
        "No single guide-label column. Each guide is a UMI count column named {TF}_CRA{k}, {TF}_CRI{k}, or NT_{k}. "
        "A cell is that guide when the column is a unique strict maximum and > 0. "
        f"Ties with a positive maximum: {n_tie} (not given a winner). "
        f"TF-guide calls that disagree with the CRA/CRI indicator: {int(disagree.sum())} (not assigned). "
        "Pooled perturbation id is {CRA|CRI}_{TF} from columns {TF}_{CRA|CRI}{k}."
    )
    pert_note = (
        "Perturbation identity is the pooled id {CRA|CRI}_{TF} over guide columns {TF}_{CRA|CRI}{k}. "
        "There is no obs column whose values are the strings CRA_<TF>."
    )
    return dict(
        labels=pd.Series(labels, index=obs.index),
        kind=pd.Series(kind, index=obs.index),
        pd_labels=pd.Series(pd_labels, index=obs.index),
        is_wt=is_wt,
        modality=modality,
        guide_series=guides_called,
        pert_col=pert_note,
        pd_col=",".join(pd_cols),
        wt_col="WT",
        wt_note=wt_note,
        gcol=None,
        gnote=gnote,
        nt_note=nt_note,
        n_tie=n_tie,
        n_disagree=int(disagree.sum()),
        n_unassigned=n_unassigned,
    )


def guide_column(df, pert_col, log):
    labels = df[pert_col].map(lambda x: _s(x).strip())
    name_hits = [c for c in df.columns if c.lower() in GUIDE_NAMES and c != pert_col]
    if len(name_hits) == 1:
        log(f"[guide] column {name_hits[0]!r} by name")
        return name_hits[0], "name"
    if len(name_hits) > 1:
        log(f"[guide] ambiguous name hits {name_hits}. Per-guide counts not resolved.")
        return None, f"ambiguous names {name_hits}"
    return None, "no guide-named column"


def unscaled_tmm_factor(obs, n_o, ref_c, n_r, logratio_trim=TMM_LOGRATIO_TRIM,
                        sum_trim=TMM_SUM_TRIM, acutoff=-1e10):
    with np.errstate(divide="ignore", invalid="ignore"):
        log_r = np.log2((obs / n_o) / (ref_c / n_r))
        abs_e = 0.5 * (np.log2(obs / n_o) + np.log2(ref_c / n_r))
        var = (n_o - obs) / n_o / obs + (n_r - ref_c) / n_r / ref_c
    keep = np.isfinite(log_r) & np.isfinite(abs_e) & (abs_e > acutoff) & np.isfinite(var) & (var > 0)
    n = int(keep.sum())
    if n < 10:
        return 1.0
    lr, ae, v = log_r[keep], abs_e[keep], var[keep]
    r_m = rankdata(lr, method="average")
    r_a = rankdata(ae, method="average")
    lo_l = n * logratio_trim + 1.0
    hi_l = n + 1.0 - lo_l
    lo_s = n * sum_trim + 1.0
    hi_s = n + 1.0 - lo_s
    trim = (r_m >= lo_l) & (r_m <= hi_l) & (r_a >= lo_s) & (r_a <= hi_s)
    if int(trim.sum()) < 2:
        return 1.0
    w = 1.0 / v[trim]
    return float(2 ** (np.sum(lr[trim] * w) / np.sum(w)))


def tmm_cache_from_panel(C, nf):
    """C is samples × genes. nf is edgeR geom-mean-scaled factors from tmm_norm_factors."""
    C = np.asarray(C, np.float64)
    lib = C.sum(1)
    Ct = C.T
    q75 = np.quantile(Ct, 0.75, axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        f75 = q75 / np.where(lib <= 0, np.nan, lib)
    ref = int(np.nanargmin(np.abs(f75 - np.nanmean(f75))))
    raw = np.ones(C.shape[0], dtype=np.float64)
    for j in range(C.shape[0]):
        if j == ref or not np.isfinite(lib[j]) or lib[j] <= 0:
            continue
        raw[j] = unscaled_tmm_factor(C[j], float(lib[j]), C[ref], float(lib[ref]))
    pos = raw[np.isfinite(raw) & (raw > 0)]
    geom = float(np.exp(np.mean(np.log(pos)))) if pos.size else 1.0
    scaled = raw / geom
    if np.max(np.abs(scaled - np.asarray(nf, float))) > 1e-6:
        raise StopStep(
            "tmm",
            "Fast TMM factors do not match tmm_norm_factors "
            f"(max abs {np.max(np.abs(scaled - np.asarray(nf, float))):.3e}). Not using a drifted null.",
        )
    ave = float(np.mean(lib * np.asarray(nf, float)))
    return dict(ref=ref, ref_counts=C[ref].copy(), ref_lib=float(lib[ref]), geom=geom, ave=ave, lib=lib, nf=np.asarray(nf, float))


def logcpm_one(counts, nf, ave, prior=EDGE_R_PRIOR):
    counts = np.asarray(counts, np.float64).ravel()
    lib = float(counts.sum())
    lib_eff = lib * float(nf)
    adj_prior = prior * lib_eff / ave if ave > 0 else prior
    den = lib_eff + 2.0 * adj_prior
    if den <= 0:
        den = 1.0
    return np.log2(np.clip((counts + adj_prior) / den * 1e6, 1e-12, None)), lib


def z_from_counts(counts, cache, mu, sd, missing):
    fac = unscaled_tmm_factor(counts, float(np.sum(counts)), cache["ref_counts"], cache["ref_lib"])
    nf = fac / cache["geom"]
    logcpm, lib = logcpm_one(counts, nf, cache["ave"])
    sd_safe = np.where(sd < 1e-12, 1.0, sd)
    z = (logcpm - mu) / sd_safe
    z = np.asarray(z, np.float64)
    z[missing] = 0.0
    return z, logcpm, nf, lib


def stats_from_d(d, v, p_unit):
    d = np.asarray(d, np.float64).ravel()
    v = np.asarray(v, np.float64).ravel()
    p_unit = np.asarray(p_unit, np.float64).ravel()
    d_r = d - np.dot(d, p_unit) * p_unit
    v_r = v - np.dot(v, p_unit) * p_unit
    return dict(
        cos=cosine_full(d, v),
        frac=frac_along(d, v),
        delta=float(np.linalg.norm(v - d) - np.linalg.norm(v)),
        cos_resid=cosine_full(d_r, v_r),
        delta_resid=float(np.linalg.norm(v_r - d_r) - np.linalg.norm(v_r)),
        d_norm=float(np.linalg.norm(d)),
        v_norm=float(np.linalg.norm(v)),
    )


def bh_q(p):
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if int(ok.sum()) == 0:
        return q
    pp = p[ok]
    m = int(pp.size)
    order = np.argsort(pp)
    ranked = pp[order] * m / np.arange(1, m + 1)
    rev = np.minimum.accumulate(ranked[::-1])[::-1]
    qq = np.empty(m)
    qq[order] = np.minimum(rev, 1.0)
    q[ok] = qq
    return q


def module_score(counts_csr, cols, lib):
    """Mean log1p(CP10k) over the listed columns. counts_csr is cells × genes."""
    if len(cols) == 0:
        return np.full(counts_csr.shape[0], np.nan)
    sub = counts_csr[:, cols]
    arr = np.asarray(sub.toarray() if sparse.issparse(sub) else sub, np.float64)
    scale = np.maximum(lib, 1.0)[:, None]
    return np.log1p(arr / scale * 1e4).mean(axis=1)


def quartile_masks(scores):
    scores = np.asarray(scores, float)
    finite = np.isfinite(scores)
    if int(finite.sum()) < 4:
        return None, None, dict(reason="fewer than 4 finite NT scores")
    q25, q75 = np.quantile(scores[finite], [0.25, 0.75])
    bottom = finite & (scores <= q25)
    top = finite & (scores >= q75)
    if bottom.sum() == 0 or top.sum() == 0 or np.any(bottom & top):
        return None, None, dict(
            reason="quartile sets empty or overlapping",
            q25=float(q25), q75=float(q75),
            n_bottom=int(bottom.sum()), n_top=int(top.sum()),
            n_overlap=int((bottom & top).sum()),
        )
    return bottom, top, dict(q25=float(q25), q75=float(q75), n_bottom=int(bottom.sum()), n_top=int(top.sum()))


def choose_symbol_column(var, log):
    cols = list(var.columns)
    for name in SYMBOL_COL_PRIORITY:
        if name in var.columns:
            log(f"[genes] symbol column {name!r} by pre-registered priority {SYMBOL_COL_PRIORITY}")
            return name
    index_vals = var.index.astype(str)
    n_ens = int(sum(strip_ensembl(x).startswith("ENS") for x in index_vals[:500]))
    if n_ens < 0.5 * min(500, len(index_vals)):
        log("[genes] no priority symbol column; using var index (not Ensembl-like)")
        return "__index__"
    raise StopStep(
        "genes",
        f"No gene-symbol column among {SYMBOL_COL_PRIORITY} and the var index looks like Ensembl. "
        f"var columns={cols}. Not substituting.",
    )


def choose_id_column(var, log):
    for name in ID_COL_PRIORITY:
        if name in var.columns:
            log(f"[genes] id column {name!r}")
            return name
    log("[genes] no id column; gene ids fall back to the symbol column for align")
    return None


def sum_cells(csr, idx, dst, n_ruler):
    idx = np.asarray(idx, int)
    out = np.zeros(n_ruler, np.float64)
    if idx.size == 0:
        return out
    s = np.asarray(csr[idx].sum(axis=0), np.float64).ravel()
    out[dst] = s
    return out


def run_posctrl(csr, dst, n_ruler, early_b, late_b, frozen, log):
    mixes = []
    late_b = np.asarray(late_b, int)
    early_b = np.asarray(early_b, int)
    n_tot = int(min(late_b.size, early_b.size))
    rng = np.random.default_rng(SEED)
    for f in F_VALUES:
        n_early = int(np.round(float(f) * n_tot))
        n_late = int(n_tot - n_early)
        take_l = rng.choice(late_b, size=n_late, replace=False) if n_late else np.array([], dtype=int)
        take_e = rng.choice(early_b, size=n_early, replace=False) if n_early else np.array([], dtype=int)
        mixes.append(dict(f=float(f), n_late=n_late, n_early=n_early, n_total=n_tot,
                          late_idx=np.asarray(take_l, int), early_idx=np.asarray(take_e, int)))
        log(f"[P] f={f:.2f} n_late={n_late} n_early={n_early} n_total={n_tot}")
    c_early = sum_cells(csr, early_b, dst, n_ruler)
    c_late = sum_cells(csr, late_b, dst, n_ruler)
    rows = [c_early, c_late]
    names = ["EARLY_B", "LATE_B"]
    for m in mixes:
        rows.append(sum_cells(csr, m["late_idx"], dst, n_ruler) + sum_cells(csr, m["early_idx"], dst, n_ruler))
        names.append(f"mix_f{m['f']:.2f}")
    C = np.vstack(rows)
    logcpm, nf = tmm_logcpm_quiet(C)
    Z, missing = zscore_frozen(logcpm, frozen, C)
    log(f"[P tmm] rows={names} nf={np.round(nf, 4).tolist()} prior.count={EDGE_R_PRIOR:g} "
        f"n_missing_z0={int(np.asarray(missing).sum())} half-A excluded")
    z_e, z_l = Z[0], Z[1]
    obs = []
    for i, m in enumerate(mixes):
        st = axis_stats(z_l, z_e, Z[2 + i])
        obs.append(st)
        log(f"[P stat] f={m['f']:.2f} frac={_fmt(st['frac'])} delta={_fmt(st['delta'])} cos={_fmt(st['cos'])}")
    ruler = ruler_score(Z, frozen)
    cache = tmm_cache_from_panel(C, nf)
    mu = np.asarray(frozen["mu"], float)
    sd = np.asarray(frozen["sd"], float)
    for i, m in enumerate(mixes):
        z_chk, _, _, _ = z_from_counts(rows[2 + i], cache, mu, sd, np.asarray(missing, bool))
        err = float(np.max(np.abs(z_chk - Z[2 + i])))
        log(f"[P z check] f={m['f']:.2f} max abs z err vs panel={err:.3e}")
    rng_b = np.random.default_rng(BOOT_SEED)
    boot = {m["f"]: dict(frac=[], delta=[], cos=[]) for m in mixes}
    # Freeze EARLY_B and LATE_B z. Resample only mixture cells. v stays the point-estimate direction.
    for b in range(N_BOOT):
        for i, m in enumerate(mixes):
            c_mix = (
                sum_cells(csr, resample_idx(m["late_idx"], rng_b), dst, n_ruler)
                + sum_cells(csr, resample_idx(m["early_idx"], rng_b), dst, n_ruler)
            )
            z_mix, _, _, _ = z_from_counts(c_mix, cache, mu, sd, np.asarray(missing, bool))
            st = axis_stats(z_l, z_e, z_mix)
            boot[m["f"]]["frac"].append(st["frac"])
            boot[m["f"]]["delta"].append(st["delta"])
            boot[m["f"]]["cos"].append(st["cos"])
        if (b + 1) % 50 == 0:
            log(f"[P boot] {b + 1}/{N_BOOT}")
    out_rows = []
    for i, m in enumerate(mixes):
        st = obs[i]
        flo, fhi = percentile_ci(boot[m["f"]]["frac"])
        dlo, dhi = percentile_ci(boot[m["f"]]["delta"])
        clo, chi = percentile_ci(boot[m["f"]]["cos"])
        out_rows.append(dict(
            f=m["f"], n_late=m["n_late"], n_early=m["n_early"], n_total=m["n_total"],
            frac=st["frac"], delta=st["delta"], cos=st["cos"],
            frac_ci_lo=flo, frac_ci_hi=fhi,
            delta_ci_lo=dlo, delta_ci_hi=dhi,
            cos_ci_lo=clo, cos_ci_hi=chi,
            frac_ci_valid=bool(ci_covers_point(st["frac"], flo, fhi)),
            delta_ci_valid=bool(ci_covers_point(st["delta"], dlo, dhi)),
            cos_ci_valid=bool(ci_covers_point(st["cos"], clo, chi)),
            ruler_score_context=float(ruler[2 + i]),
            noise_floor=bool(abs(m["f"]) < 1e-12),
        ))
    pdf = pd.DataFrame(out_rows)
    pdf.to_csv(SENG_DIR / "posctrl.csv", index=False)
    fracs = [float(r["frac"]) for r in out_rows]
    row50 = next(r for r in out_rows if abs(r["f"] - 0.50) < 1e-12)
    row0 = next(r for r in out_rows if abs(r["f"]) < 1e-12)
    mono = monotonic_increasing(fracs)
    frac_ci_ok = bool(row50["frac_ci_valid"]) and ci_excludes_zero(row50["frac_ci_lo"], row50["frac_ci_hi"])
    delta_neg = bool(np.isfinite(row50["delta"]) and row50["delta"] < 0)
    passed = bool(mono and frac_ci_ok and delta_neg)
    reasons = []
    if not mono:
        reasons.append(f"frac not monotonically increasing in f: {[_fmt(x) for x in fracs]}")
    if not row50["frac_ci_valid"]:
        reasons.append(
            "frac(f=0.50) CI INVALID (does not contain the point estimate) "
            f"point={_fmt(row50['frac'])} {_fmt_ci(row50['frac_ci_lo'], row50['frac_ci_hi'])}"
        )
    elif not ci_excludes_zero(row50["frac_ci_lo"], row50["frac_ci_hi"]):
        reasons.append(
            f"frac(f=0.50) CI does not exclude 0: {_fmt_ci(row50['frac_ci_lo'], row50['frac_ci_hi'])}"
        )
    if not delta_neg:
        reasons.append(f"delta(f=0.50)={_fmt(row50['delta'])} is not < 0")
    summary = dict(
        pass_=passed, monotonic=mono, frac50_ci_excludes_0=bool(frac_ci_ok),
        delta50_neg=delta_neg, frac50=row50["frac"], delta50=row50["delta"],
        frac50_ci_lo=row50["frac_ci_lo"], frac50_ci_hi=row50["frac_ci_hi"],
        frac50_ci_valid=row50["frac_ci_valid"],
        frac0=row0["frac"], delta0=row0["delta"], n_total=n_tot,
        fracs=fracs, reason="; ".join(reasons),
        n_boot=N_BOOT, seed=SEED, boot_seed=BOOT_SEED,
    )
    dump_json(SENG_DIR / "posctrl_summary.json", jsonable(summary))
    log(f"[P] pass={passed} monotonic={mono} frac50_ci_excludes_0={frac_ci_ok} "
        f"delta50_neg={delta_neg} reason={summary['reason'] or 'NA'}")
    return pdf, summary


def panel_z(rows, frozen):
    C = np.vstack([np.asarray(r, np.float64).ravel() for r in rows])
    logcpm, nf = tmm_logcpm_quiet(C)
    Z, missing = zscore_frozen(logcpm, frozen, C)
    return Z, missing, nf, C


def classify_hit(row):
    """Exactly one label. Order is the pre-registration. Returns the label."""
    cos_sig = bool(np.isfinite(row["q_cos"]) and row["q_cos"] <= P_BAR)
    delta_neg = bool(np.isfinite(row["delta"]) and row["delta"] < 0)
    delta_sig = bool(np.isfinite(row["q_delta"]) and row["q_delta"] <= P_BAR)
    dres_neg = bool(np.isfinite(row["delta_resid"]) and row["delta_resid"] < 0)
    dres_p = bool(np.isfinite(row["p_delta_resid"]) and row["p_delta_resid"] <= P_BAR)
    cres_null = bool(np.isfinite(row["p_cos_resid"]) and row["p_cos_resid"] > P_BAR)
    if delta_neg and delta_sig and dres_neg and dres_p:
        return "toward_young"
    if cos_sig and cres_null:
        return "proliferation_explained"
    if cos_sig and ((not delta_neg) or (not delta_sig)):
        return "away_from_old_only"
    return "no_effect"


def meets_toward(row):
    return classify_hit(row) == "toward_young" and (
        bool(np.isfinite(row["delta"]) and row["delta"] < 0)
        and bool(np.isfinite(row["q_delta"]) and row["q_delta"] <= P_BAR)
        and bool(np.isfinite(row["delta_resid"]) and row["delta_resid"] < 0)
        and bool(np.isfinite(row["p_delta_resid"]) and row["p_delta_resid"] <= P_BAR)
    )


def write_progress(stop, next_action, extra_lines=None):
    man = load_manifest()
    lines = [
        "# PROGRESS_SENG",
        "",
        f"**STOP status:** {stop}",
        f"**Next action:** {next_action}",
        "",
        "## Seeds / gates",
        "",
        f"- seed `{SEED}` (splits, nulls, permutations)  boot `{BOOT_SEED}`",
        f"- n_null={N_NULL}  n_perm={N_PERM}  n_boot={N_BOOT}",
        "- Uncalibrated R2 is never a gate. Frozen-ruler score is context only.",
        f"- MIN_CELLS={MIN_CELLS}  MIN_WT={MIN_WT}  P_BAR={P_BAR}",
        f"- {SOUTHARD_BAN}",
        f"- {GSE325735_BAN}",
        f"- PREREG.flag exists={PREREG_FLAG.exists()} (not rewritten)",
        f"- failures recorded: {len(man.get('failures') or [])}",
        "",
    ]
    if extra_lines:
        lines.extend(extra_lines)
    PROGRESS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _failures_md(man):
    fails = man.get("failures") or []
    if not fails:
        return ["None recorded.", ""]
    lines = []
    for f in fails:
        lines.append(f"- **{f.get('step')}:** {f.get('message')}")
    lines.append("")
    return lines


def write_findings(status, sections):
    prereg = load_prereg_text()
    lines = [
        "# FINDINGS_SENG — do Sengstack et al. (PNAS 2026) \"rejuvenating\" TF perturbations move fibroblasts TOWARD young cells?",
        "",
        f"**Status:** {status}",
        "",
        "## Pre-registration (verbatim, written before any SENG statistic)",
        "",
        "Written before any SENG statistic. No threshold in this block is re-tuned after numbers exist. The block below is the task text, verbatim.",
        "",
        "```",
        prereg,
        "```",
        "",
        f"Flag: `results/seng/PREREG.flag` exists={PREREG_FLAG.exists()}.",
        "",
    ]
    lines.extend(sections)
    FINDINGS_PATH.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def limitations_md():
    return [
        "## Limitations",
        "",
        "- Neonatal cells aged by passaging, not donor aging. Passage reversal is not a donor-age effect.",
        "- One cell line.",
        "- Additivity assumption in delta: delta applies the TF-minus-NT vector to a WT LATE cell. "
        "That assumes the TF effect is additive to the construct effect (the construct cancels in "
        "d_g = z(g) − z(NT)) and that the same additive shift applies on the WT background.",
        "- Pseudobulk per perturbation: each perturbation is one summed count vector, so cell-to-cell "
        "heterogeneity inside the guide pool is not the statistic.",
        "- TMM for positive control P uses half-B rows only. The analysis panel (half-A, NT, quartile "
        "pseudobulks, qualifying perturbations) is a separate TMM. A joint TMM would put half-A inside P "
        "and half-B inside the young direction.",
        "- Null and bootstrap z-scores recompute a TMM factor for the resampled pseudobulk against the "
        "analysis-panel reference. v and p stay the panel point estimates.",
        "- Proliferation score is the sum of mean log1p(CP10k) of Seurat cc.genes S genes and of G2M genes "
        "present in this file. Quartile cut is the 25th and 75th percentile of NT scores within modality.",
        "- Frozen-ruler score is a context column only.",
        "- This file cannot test Sengstack's in vivo mouse results.",
        f"- {SOUTHARD_BAN}",
        f"- {GSE325735_BAN}",
        "",
    ]


def inventory_sections(inv):
    lines = [
        "## Stage 0 inventory (descriptive; no statistic)",
        "",
        f"- File: `{inv.get('path')}`",
        f"- Size: {inv.get('size_bytes')} bytes. md5: `{inv.get('md5')}`.",
        f"- Figshare API size: {FIGSHARE_API_SIZE}. URL: {FIGSHARE_URL}.",
        f"- Matrix used: {inv.get('matrix_source')}.",
        f"- Raw-count audit: `{jsonable(inv.get('raw_audit'))}`.",
        "",
        "### Columns",
        "",
        f"- Perturbation labels: `{inv.get('pert_col')}`.",
        f"- Guide identity: `{inv.get('guide_col')}` ({inv.get('guide_note')}).",
        f"- Population doubling: `{inv.get('pd_col')}`.",
        f"- WT status: `{inv.get('wt_note')}`.",
        f"- NT assignment: {inv.get('nt_note')}",
        f"- Perturbed-cell PD: {inv.get('perturbed_pd_note')}",
        f"- EARLY: {inv.get('early_note')}",
        f"- LATE: {inv.get('late_note')}",
        "",
    ]
    for title, key in (
        ("### WT cells per PD", "wt_pd_csv"),
        ("### NT cells per modality", "nt_csv"),
        ("### Perturbation counts (pooled over guides)", "pert_csv"),
    ):
        lines.append(title)
        lines.append("")
        p = inv.get(key)
        if p and Path(p).exists():
            lines.append(md_table(pd.read_csv(p)))
            lines.append("")
        else:
            lines.append("Not written.")
            lines.append("")
    if inv.get("guide_csv") and Path(inv["guide_csv"]).exists():
        lines.append("### Per-guide counts")
        lines.append("")
        lines.append(f"Full table: `{inv['guide_csv']}`.")
        gdf = pd.read_csv(inv["guide_csv"])
        lines.append(f"Rows: {len(gdf)}.")
        lines.append("")
        if len(gdf) <= 80:
            lines.append(md_table(gdf))
            lines.append("")
    other_csv = inv.get("other_csv")
    if other_csv and Path(other_csv).exists():
        odf = pd.read_csv(other_csv)
        lines.append("### Other labels in the perturbation column (not CRA_/CRI_/NT_/WT)")
        lines.append("")
        if len(odf) == 0:
            lines.append("None.")
        elif len(odf) <= 40:
            lines.append(md_table(odf))
        else:
            lines.append(f"{len(odf)} labels. Full table: `{other_csv}`.")
        lines.append("")
    return lines


def main():
    t0 = time.perf_counter()
    log = Logger(SENG_DIR / "run_report.txt", mode="w")
    log(f"[start] PREREG.flag exists={PREREG_FLAG.exists()} {FROZEN_RULER} exists={FROZEN_RULER.exists()}")
    log(SOUTHARD_BAN)
    log(GSE325735_BAN)
    load_prereg_text()
    man = load_manifest()
    man["status"] = "RUNNING"
    man["seed"] = SEED
    man["boot_seed"] = BOOT_SEED
    save_manifest(man)

    if not H5_PATH.exists():
        record_failure("download", f"missing {H5_PATH}")
        write_progress("STOP", f"Place the figshare h5ad at {H5_PATH}.")
        write_findings("STOP. h5ad missing.", [
            "## STOP / failures", "",
            *_failures_md(load_manifest()),
        ])
        log.close()
        return

    size = H5_PATH.stat().st_size
    log(f"[download] path={H5_PATH} size={size} expected_api_size={FIGSHARE_API_SIZE}")
    if size != FIGSHARE_API_SIZE:
        msg = (
            f"File size {size} != figshare API size {FIGSHARE_API_SIZE}. "
            "Not analyzing a truncated or substituted file."
        )
        record_failure("download", msg)
        write_progress("STOP", "Re-download the figshare h5ad; size does not match the API record.")
        write_findings("STOP. Download size mismatch.", [
            "## STOP / failures", "", *_failures_md(load_manifest()),
        ])
        log.close()
        return
    md5 = file_md5(H5_PATH)
    log(f"[download] md5={md5}")
    dump_json(SENG_DIR / "download.json", dict(
        path=str(H5_PATH), size_bytes=int(size), md5=md5,
        url=FIGSHARE_URL, api_size=FIGSHARE_API_SIZE, filename=H5_NAME,
    ))

    frozen = load_frozen_ruler()
    mu = np.asarray(frozen["mu"], float)
    sd = np.asarray(frozen["sd"], float)
    n_ruler = int(mu.shape[0])

    with h5py.File(H5_PATH, "r") as f:
        tree = dump_tree(f, SENG_DIR / "h5_tree.txt")
        log(f"[h5] top keys={list(f.keys())} tree_lines={len(tree)}")
        obs = read_obs_var(f["obs"])
        var = read_obs_var(f["var"])
        log(f"[obs] n={len(obs)} columns={list(obs.columns)}")
        log(f"[var] n={len(var)} columns={list(var.columns)}")
        try:
            mat, raw_audit, matrix_source, audits = locate_count_matrix(f)
        except StopStep as e:
            record_failure(e.step, str(e), getattr(e, "details", None) if hasattr(e, "details") else None)
            # StopStep in this repo is (step, message) or (step, message, details)
            write_inventory_stop(obs, var, None, str(e), log, size, md5)
            log.close()
            return
        if isinstance(mat, tuple) and mat[0] in ("DENSE_X", "DENSE_LAYER"):
            record_failure(
                "h5",
                f"Count matrix {matrix_source} is dense. This runner reads sparse csr/csc raw counts. "
                "Not densifying a multi-GB matrix and not substituting a normalized copy.",
            )
            write_progress("STOP", "Sparse raw counts were not found; dense path was not improvised.")
            write_findings("STOP. Dense count matrix was not improvised.", [
                "## STOP / failures", "", *_failures_md(load_manifest()),
                f"Matrix source: {matrix_source}. Audit: `{jsonable(raw_audit)}`.",
                "",
            ])
            log.close()
            return
        log(f"[X] source={matrix_source} shape={mat.shape} nnz={mat.nnz} audit={raw_audit}")
        if mat.shape[0] != len(obs):
            record_failure(
                "h5",
                f"X rows {mat.shape[0]} != obs rows {len(obs)}. Not repairing rows.",
            )
            write_progress("STOP", "X/obs row mismatch. Do not repair.")
            write_findings("STOP. X/obs row mismatch.", [
                "## STOP / failures", "", *_failures_md(load_manifest()),
            ])
            log.close()
            return

        wide = is_wide_guide_obs(obs.columns)
        log(f"[obs layout] wide_guide_columns={wide}")
        guide_series = None
        if wide:
            parsed = parse_wide_obs(obs, log)
            labels = parsed["labels"]
            kind = parsed["kind"]
            pd_labels = parsed["pd_labels"]
            is_wt = parsed["is_wt"]
            modality = parsed["modality"]
            guide_series = parsed["guide_series"]
            pert_col = parsed["pert_col"]
            pd_col = parsed["pd_col"]
            wt_col = parsed["wt_col"]
            wt_note = parsed["wt_note"]
            gcol = parsed["gcol"]
            gnote = parsed["gnote"]
            nt_note = parsed["nt_note"]
            _write_count_tables(
                obs, labels, kind, pd_labels, is_wt, modality, gcol, log, guide_series=guide_series,
            )
        else:
            try:
                pert_col, frac_tab = pick_column(obs, "perturbation", log)
                pd_col, _ = pick_column(obs, "pd", log)
                wt_col, _ = pick_column(obs, "wt", log)
            except StopStep as e:
                record_failure(e.step, e.args[0] if e.args else str(e))
                write_progress("STOP", "Column assignment is ambiguous. Inventory only.")
                write_findings("STOP. Ambiguous columns.", [
                    "## STOP / failures", "", *_failures_md(load_manifest()),
                ])
                log.close()
                return
            frac_tab.to_csv(SENG_DIR / "inventory_column_fractions.csv", index=False)

            if pert_col is None:
                record_failure("inventory", "No column contains CRA_*/CRI_*/NT_* labels. Not substituting.")
                write_progress("STOP", "Perturbation column not found.")
                write_findings("STOP. No perturbation column.", [
                    "## STOP / failures", "", *_failures_md(load_manifest()),
                ])
                log.close()
                return
            if pd_col is None:
                name_hits = []
                for c in obs.columns:
                    if not PD_NAME_RE.search(str(c)):
                        continue
                    nums = [pd_number(v) for v in obs[c].map(lambda x: _s(x).strip()).unique()]
                    nums = [n for n in nums if n is not None]
                    if len(set(nums)) >= 2:
                        name_hits.append(c)
                if len(name_hits) == 1:
                    pd_col = name_hits[0]
                    log(f"[PD] no PD14-style strings; using numeric passage column {pd_col!r}")
                elif len(name_hits) > 1:
                    record_failure(
                        "no_wt_reference",
                        f"Multiple passage-like columns {name_hits}. Not substituting one.",
                    )
                    write_progress("STOP", "key=`no_wt_reference`. Ambiguous passage columns.")
                    write_findings("STOP. key=`no_wt_reference`.", [
                        "## Failures (verbatim)", "", *_failures_md(load_manifest()),
                    ])
                    log.close()
                    return
            if pd_col is None:
                record_failure(
                    "no_wt_reference",
                    "No column contains PD14/PD26/PD32-style labels, so WT cells at >= 2 PDs cannot be identified. "
                    "Not substituting a passage column.",
                )
                write_progress("STOP", "PD column not found. key=no_wt_reference.")
                write_findings("STOP. key=`no_wt_reference`.", [
                    "## STOP / failures", "", *_failures_md(load_manifest()),
                    "Stage 0 could not find a population-doubling column. Inventory fractions:",
                    "",
                    md_table(frac_tab),
                    "",
                ])
                log.close()
                return

            labels = obs[pert_col].map(lambda x: _s(x).strip())
            pd_labels = obs[pd_col].map(lambda x: _s(x).strip())
            kind = labels.map(lambda s: classify_token(s)[0])
            # WT: explicit WT column, else WT token inside the perturbation column.
            if wt_col is not None and wt_col != pert_col:
                wt_vals = obs[wt_col].map(lambda x: _s(x).strip())
                is_wt = wt_vals.map(lambda s: s.upper() in WT_TOKENS).to_numpy()
                wt_note = f"column {wt_col!r}; WT tokens {sorted(WT_TOKENS)}"
            elif (kind == "WT").any():
                is_wt = (kind == "WT").to_numpy()
                wt_note = f"WT token inside perturbation column {pert_col!r}"
                wt_col = pert_col
            else:
                is_wt = np.zeros(len(obs), dtype=bool)
                wt_note = "no WT token found in any column"
            log(f"[WT] n={int(is_wt.sum())} note={wt_note}")

            try:
                modality, nt_note = assign_nt_modality(obs, pert_col, log)
            except StopStep as e:
                record_failure("nt_unassignable", e.args[0] if e.args else str(e))
                _write_partial_inventory(
                    obs, labels, kind, pd_labels, is_wt, None, pert_col, pd_col, wt_col, wt_note,
                    "NT unassigned", size, md5, matrix_source, raw_audit, log,
                )
                write_progress("STOP", "key=`nt_unassignable`. Inventory only.")
                inv = load_json(SENG_DIR / "inventory.json")
                write_findings("STOP. key=`nt_unassignable`.", [
                    "## STOP / failures", "", *_failures_md(load_manifest()),
                    *inventory_sections(inv),
                ])
                log.close()
                return

            gcol, gnote = guide_column(obs, pert_col, log)
            _write_count_tables(obs, labels, kind, pd_labels, is_wt, modality, gcol, log)

        wt_pd = (
            pd.DataFrame({"pd": pd_labels[is_wt]})
            .groupby("pd", dropna=False)
            .size()
            .rename("n_cells")
            .reset_index()
        )
        wt_pd_levels = []
        for _, r in wt_pd.iterrows():
            num = pd_number(r["pd"])
            if num is not None and int(r["n_cells"]) > 0:
                wt_pd_levels.append((num, str(r["pd"]), int(r["n_cells"])))
        wt_pd_levels.sort()
        log(f"[WT PD] {wt_pd_levels}")
        if len(wt_pd_levels) < 2:
            record_failure(
                "no_wt_reference",
                f"WT (untransduced) cells do not exist at >= 2 distinct PDs. "
                f"WT PD counts={wt_pd_levels}. n_WT={int(is_wt.sum())}.",
            )
            inv = _finish_inv(
                size, md5, matrix_source, raw_audit, pert_col, pd_col, wt_col, wt_note,
                gcol, gnote, nt_note, "not set", "not set", "not set",
            )
            write_progress("STOP", "key=`no_wt_reference`. Inventory only.")
            write_findings("STOP. key=`no_wt_reference`.", [
                "## STOP / failures", "", *_failures_md(load_manifest()),
                *inventory_sections(inv),
            ])
            log.close()
            return

        pert_mask = ((kind == "CRA") | (kind == "CRI")).to_numpy()
        pert_pds = sorted(set(pd_labels[pert_mask].tolist()) - {""})
        pert_pd_nums = [(pd_number(p), p) for p in pert_pds]
        pert_pd_nums = [(n, p) for n, p in pert_pd_nums if n is not None]
        n_pert_blank_pd = int(((pd_labels == "") & pert_mask).sum())
        if len(pert_pd_nums) == 1:
            late_pd_num, late_pd_label = pert_pd_nums[0]
            late_note = (
                f"LATE = WT at {late_pd_label}, the only PD on perturbed cells. "
                f"Perturbed cells with blank PD: {n_pert_blank_pd}."
            )
            late_determined = True
        else:
            late_pd_num, late_pd_label = wt_pd_levels[-1][0], wt_pd_levels[-1][1]
            late_note = (
                "The PD of the perturbed cells cannot be determined "
                f"(perturbed PD values={pert_pds}, blank PD n={n_pert_blank_pd}). "
                f"LATE = WT at the highest PD, {late_pd_label}."
            )
            late_determined = False
        early_pd_num, early_pd_label, early_n = wt_pd_levels[0]
        early_note = f"EARLY = WT at the lowest PD {early_pd_label} (n={early_n})."
        log(f"[EARLY] {early_note}")
        log(f"[LATE] {late_note} determined={late_determined}")

        early_idx = np.flatnonzero(is_wt & (pd_labels.to_numpy() == early_pd_label))
        late_idx = np.flatnonzero(is_wt & (pd_labels.to_numpy() == late_pd_label))
        if early_pd_label == late_pd_label:
            record_failure(
                "no_wt_reference",
                f"EARLY and LATE resolved to the same PD {early_pd_label}. "
                "Not using one PD as both ends.",
            )
            inv = _finish_inv(
                size, md5, matrix_source, raw_audit, pert_col, pd_col, wt_col, wt_note,
                gcol, gnote, nt_note, early_note, late_note, late_note,
            )
            write_progress("STOP", "key=`no_wt_reference` (EARLY PD equals LATE PD).")
            write_findings("STOP. key=`no_wt_reference`.", [
                "## STOP / failures", "", *_failures_md(load_manifest()),
                *inventory_sections(inv),
            ])
            log.close()
            return
        if early_idx.size < MIN_WT or late_idx.size < MIN_WT:
            record_failure(
                "wt_too_small",
                f"EARLY n={early_idx.size} LATE n={late_idx.size}; minimum is {MIN_WT}.",
            )
            inv = _finish_inv(
                size, md5, matrix_source, raw_audit, pert_col, pd_col, wt_col, wt_note,
                gcol, gnote, nt_note, early_note, late_note,
                f"perturbed PD determined={late_determined}",
            )
            write_progress("STOP", "key=`wt_too_small`.")
            write_findings("STOP. key=`wt_too_small`.", [
                "## STOP / failures", "", *_failures_md(load_manifest()),
                *inventory_sections(inv),
            ])
            log.close()
            return

        # Overlap of roles. Do not assign a cell twice.
        nt_mask = (kind == "NT").to_numpy()
        overlap = int((is_wt & (nt_mask | pert_mask)).sum())
        if overlap:
            record_failure(
                "inventory",
                f"WT labels overlap NT or perturbed labels on {overlap} cells. Not assigning a cell to both roles.",
            )
            inv = _finish_inv(
                size, md5, matrix_source, raw_audit, pert_col, pd_col, wt_col, wt_note,
                gcol, gnote, nt_note, early_note, late_note, "role overlap",
            )
            write_progress("STOP", "WT overlaps NT/perturbed labels. Not repairing rows.")
            write_findings("STOP. WT overlaps another role.", [
                "## STOP / failures", "", *_failures_md(load_manifest()),
                *inventory_sections(inv),
            ])
            log.close()
            return

        sym_col = choose_symbol_column(var, log)
        id_col = choose_id_column(var, log)
        if sym_col == "__index__":
            symbols = var.index.astype(str).to_numpy()
        else:
            symbols = var[sym_col].astype(str).to_numpy()
        if id_col is None:
            gene_ids = symbols
        else:
            gene_ids = var[id_col].astype(str).to_numpy()
        if len(symbols) != mat.shape[1]:
            record_failure(
                "genes",
                f"var rows {len(symbols)} != X columns {mat.shape[1]}. Not repairing.",
            )
            write_progress("STOP", "var/X mismatch.")
            write_findings("STOP. var/X mismatch.", [
                "## STOP / failures", "", *_failures_md(load_manifest()),
            ])
            log.close()
            return

        ruler_idx = ruler_column_index(gene_ids, symbols, frozen, log)
        present = np.flatnonzero(ruler_idx >= 0)
        src = ruler_idx[present]
        n_overlap = int(present.size)
        log(f"[genes] ruler={n_ruler} present_in_file={n_overlap}")
        # Verify against align_counts_to_ruler on a one-column total.
        totals = np.asarray(mat.sum(axis=0), np.float64).ravel()
        Y_check, rec = align_counts_to_ruler(
            totals.reshape(-1, 1), gene_ids, symbols, frozen, log, tag="seng_check",
        )
        mapped = np.zeros(n_ruler, np.float64)
        mapped[present] = totals[src]
        if int(rec["n_overlap"]) != n_overlap or np.max(np.abs(mapped - Y_check[0])) > 1e-6:
            record_failure(
                "genes",
                "Ruler alignment does not match fibro2_common.align_counts_to_ruler. Not continuing.",
            )
            write_progress("STOP", "Gene alignment mismatch.")
            write_findings("STOP. Gene alignment mismatch.", [
                "## STOP / failures", "", *_failures_md(load_manifest()),
            ])
            log.close()
            return

        lib_total = np.asarray(mat.sum(axis=1), np.float64).ravel()
        csr = mat[:, src].tocsr()
        dst = present
        del mat

        inv = _finish_inv(
            size, md5, matrix_source, raw_audit, pert_col, pd_col, wt_col, wt_note,
            gcol, gnote, nt_note, early_note, late_note,
            f"perturbed PD values={pert_pds}; determined={late_determined}; blank={n_pert_blank_pd}",
        )
        inv["n_ruler"] = n_ruler
        inv["n_ruler_present"] = n_overlap
        inv["symbol_column"] = sym_col
        inv["id_column"] = id_col
        inv["align"] = rec
        dump_json(SENG_DIR / "inventory.json", jsonable(inv))

        # Split. Independent Generator(SEED) inside split_half, once per anchor.
        early_a, early_b = split_half(early_idx, SEED)
        late_a, late_b = split_half(late_idx, SEED)
        log(f"[split] EARLY A={early_a.size} B={early_b.size} LATE A={late_a.size} B={late_b.size}")
        if len(set(early_a).intersection(early_b)) or len(set(late_a).intersection(late_b)):
            raise StopStep("split", "Half A and half B overlap. Not continuing.")
        split_rec = dict(
            early_pd=early_pd_label, late_pd=late_pd_label,
            n_early=int(early_idx.size), n_early_a=int(early_a.size), n_early_b=int(early_b.size),
            n_late=int(late_idx.size), n_late_a=int(late_a.size), n_late_b=int(late_b.size),
            seed=SEED, late_determined=bool(late_determined),
        )
        dump_json(SENG_DIR / "split.json", split_rec)
        np.savez_compressed(
            SENG_DIR / "split_idx.npz",
            early_a=early_a, early_b=early_b, late_a=late_a, late_b=late_b,
        )
        pd.DataFrame([split_rec]).to_csv(SENG_DIR / "split_table.csv", index=False)

        s_genes, g2m_genes = load_cc_genes(log)
        sym_upper = np.array([str(s).upper() for s in symbols])
        # First-occurrence map, same rule as symbol fallback (no alias substitution).
        pos = {}
        for i, s in enumerate(sym_upper):
            pos.setdefault(s, i)
        s_cols_file = [pos[g] for g in s_genes if g in pos]
        g_cols_file = [pos[g] for g in g2m_genes if g in pos]
        s_found = [g for g in s_genes if g in pos]
        g_found = [g for g in g2m_genes if g in pos]
        log(f"[cc.genes present] S {len(s_found)}/{len(s_genes)} G2M {len(g_found)}/{len(g2m_genes)} "
            f"missing_S={[g for g in s_genes if g not in pos]} missing_G2M={[g for g in g2m_genes if g not in pos]}")
        # Score from the ruler-aligned matrix only if the cc genes sit in it; they may not all be
        # ruler genes. Score from the file columns directly. Rebuild a thin csr of those columns
        # from csr only when the file column is one of src. File columns not in src are absent
        # from csr. Use the original column positions via a second slice stored above — csr is
        # already restricted. Map file col -> position in csr.
        file_to_sub = {int(src[i]): i for i in range(len(src))}
        s_sub = [file_to_sub[c] for c in s_cols_file if c in file_to_sub]
        g_sub = [file_to_sub[c] for c in g_cols_file if c in file_to_sub]
        # Genes present in the file but dropped because they are not ruler genes still exist in
        # `pos`. If a cc gene is in the file and not in the ruler subset, say so; do not drop it
        # silently from the score. Re-open is avoided: csr only has ruler columns. Record missing.
        s_not_in_ruler = [g for g in s_found if pos[g] not in file_to_sub]
        g_not_in_ruler = [g for g in g_found if pos[g] not in file_to_sub]
        if s_not_in_ruler or g_not_in_ruler:
            log(f"[cc.genes] in file but not in ruler subset used for counts: S={s_not_in_ruler} G2M={g_not_in_ruler}")
        # CP10k uses each cell's total UMI (all genes in X), not the ruler subset.
        score = module_score(csr, s_sub, lib_total) + module_score(csr, g_sub, lib_total)
        nt_modality = {}
        quartile_info = dict(
            s_found=s_found, g2m_found=g_found,
            s_missing=[g for g in s_genes if g not in pos],
            g2m_missing=[g for g in g2m_genes if g not in pos],
            library="total UMI in the raw count matrix",
        )
        q_idx = {}
        for mod in ("CRA", "CRI"):
            nt_idx = np.flatnonzero(nt_mask & (modality == mod))
            nt_modality[mod] = nt_idx
            bottom, top, info = quartile_masks(score[nt_idx])
            quartile_info[mod] = info
            quartile_info[mod]["n_nt"] = int(nt_idx.size)
            if bottom is None:
                record_failure(
                    "proliferation",
                    f"NT {mod} quartile sets could not be formed. {info}",
                )
                write_progress("STOP", f"Proliferation quartiles failed for {mod}. Cut was not changed.")
                write_findings("STOP. Proliferation quartiles could not be formed.", [
                    "## STOP / failures", "", *_failures_md(load_manifest()),
                    *inventory_sections(inv),
                    "",
                    f"Ruler genes present in this file: {n_overlap} / {n_ruler}.",
                    "",
                ])
                log.close()
                return
            q_idx[mod] = dict(bottom=nt_idx[bottom], top=nt_idx[top])
            log(f"[quartile {mod}] n_NT={nt_idx.size} bottom={info['n_bottom']} top={info['n_top']} "
                f"q25={info['q25']:.4f} q75={info['q75']:.4f}")
        dump_json(SENG_DIR / "prolif_quartiles.json", jsonable(quartile_info))
        np.savez_compressed(
            SENG_DIR / "prolif_quartile_idx.npz",
            CRA_bottom=q_idx["CRA"]["bottom"], CRA_top=q_idx["CRA"]["top"],
            CRI_bottom=q_idx["CRI"]["bottom"], CRI_top=q_idx["CRI"]["top"],
            score=score,
        )
        log("[quartile] sets frozen before positive control P and before any perturbation vector")

        # Positive control P BEFORE any perturbation vector.
        pdf, psum = run_posctrl(csr, dst, n_ruler, early_b, late_b, frozen, log)
        if not psum["pass_"]:
            record_failure("pos_control_broken", psum["reason"] or "P failed")
            write_progress("STOP", "key=`pos_control_broken`. No perturbation vector was computed.")
            write_findings(
                "STOP. key=`pos_control_broken`.",
                _pre_perturb_sections(inv, split_rec, n_overlap, n_ruler, pdf, psum, load_manifest()),
            )
            log.close()
            return

        # Analysis panel. Half-B is not in this panel.
        groups = []
        def add_group(name, idx, kind_name):
            groups.append(dict(name=name, idx=np.asarray(idx, int), kind=kind_name))

        add_group("EARLY_A", early_a, "anchor")
        add_group("LATE_A", late_a, "anchor")
        for mod in ("CRA", "CRI"):
            add_group(f"NT_{mod}", nt_modality[mod], "nt")
            add_group(f"TOP_{mod}", q_idx[mod]["top"], "quartile")
            add_group(f"BOT_{mod}", q_idx[mod]["bottom"], "quartile")

        # Perturbations with n >= MIN_CELLS, pooled by the label as stored.
        pert_labels = labels.to_numpy()
        cell_mod = modality
        qual = []
        for lab, sub in pd.Series(np.arange(len(labels))).groupby(pert_labels, sort=True).groups.items():
            lab = str(lab)
            if not (CRA_RE.fullmatch(lab) or CRI_RE.fullmatch(lab)):
                continue
            idx = np.asarray(sub, int)
            mod = "CRA" if lab.startswith("CRA_") else "CRI"
            # Do not re-label. Modality of the cells must agree with the label prefix.
            if not np.all(cell_mod[idx] == mod):
                record_failure(
                    "inventory",
                    f"{lab}: cell modality disagrees with the label prefix. Not re-labelling. Reported and skipped.",
                )
                continue
            if idx.size < MIN_CELLS:
                continue
            qual.append((mod, lab, idx))
        qual.sort(key=lambda t: (t[0], t[1]))
        for mod, lab, idx in qual:
            add_group(lab, idx, "perturbation")
        log(f"[panel] qualifying perturbations n>={MIN_CELLS}: {len(qual)}")

        counts = [sum_cells(csr, g["idx"], dst, n_ruler) for g in groups]
        Z, missing, nf, C = panel_z(counts, frozen)
        cache = tmm_cache_from_panel(C, nf)
        log(f"[panel tmm] n_rows={C.shape[0]} ref={cache['ref']} ave={cache['ave']:.4f} "
            f"n_missing_z0={int(missing.sum())}")
        miss_chk = np.asarray(missing, bool)
        for i in range(min(5, C.shape[0])):
            z_chk, _, _, _ = z_from_counts(C[i], cache, mu, sd, miss_chk)
            err = float(np.max(np.abs(z_chk - Z[i])))
            log(f"[panel z check] row={groups[i]['name']} max abs z err={err:.3e}")
            if err > 1e-4:
                record_failure("tmm", f"z_from_counts disagrees with panel z on {groups[i]['name']} (max abs {err:.3e})")
                write_progress("STOP", "Null z-score path does not match the panel. Not running nulls.")
                write_findings("STOP. z-score path mismatch.", [
                    "## Failures (verbatim)", "", *_failures_md(load_manifest()),
                    *inventory_sections(inv),
                ])
                log.close()
                return
        # Confirm single-sample logcpm matches the panel for every row.
        max_err = 0.0
        for i in range(C.shape[0]):
            logcpm_i, _ = logcpm_one(C[i], cache["nf"][i], cache["ave"])
            # Rebuild from z path's official logcpm: invert is unnecessary; recompute official row.
            max_err = max(max_err, float(np.max(np.abs(logcpm_i - _official_logcpm_row(C, nf)[i]))))
        log(f"[panel tmm] single-sample logcpm max abs err vs panel={max_err:.3e}")
        if max_err > 1e-5:
            record_failure("tmm", f"single-sample logcpm disagrees with the panel (max abs {max_err:.3e})")
            write_progress("STOP", "TMM path mismatch. Nulls were not run.")
            write_findings("STOP. TMM path mismatch.", [
                "## STOP / failures", "", *_failures_md(load_manifest()),
                *inventory_sections(inv),
            ])
            log.close()
            return

        name_to_i = {g["name"]: i for i, g in enumerate(groups)}
        z_early = Z[name_to_i["EARLY_A"]]
        z_late = Z[name_to_i["LATE_A"]]
        v = z_early - z_late
        logcpm_panel = _official_logcpm_row(C, nf)
        ruler_ctx = ruler_score(Z, frozen)

        # p before any d_g is formed.
        p_unit = {}
        cos_vp = {}
        for mod in ("CRA", "CRI"):
            z_top = Z[name_to_i[f"TOP_{mod}"]]
            z_bot = Z[name_to_i[f"BOT_{mod}"]]
            pvec = z_top - z_bot
            nrm = float(np.linalg.norm(pvec))
            if nrm < 1e-12:
                record_failure("proliferation", f"p_{mod} has zero norm. Not renormalizing by another rule.")
                write_progress("STOP", "Proliferation axis has zero length.")
                write_findings("STOP. Proliferation axis degenerate.", [
                    "## STOP / failures", "", *_failures_md(load_manifest()),
                ])
                log.close()
                return
            p_unit[mod] = pvec / nrm
            cos_vp[mod] = cosine_full(v, p_unit[mod])
            log(f"[axis] cos(v, p_{mod})={_fmt(cos_vp[mod])} ||p|| before norm={nrm:.4f} ||v||={np.linalg.norm(v):.4f}")

        anch = np.load(ANCHORS_PATH, allow_pickle=True)
        u = np.asarray(anch["u"], float).ravel()
        if u.shape[0] != n_ruler:
            record_failure("anchors", f"u_TOWARD length {u.shape[0]} != ruler {n_ruler}. Not aligning ad hoc.")
            write_progress("STOP", "u_TOWARD gene length mismatch.")
            write_findings("STOP. anchors length mismatch.", [
                "## STOP / failures", "", *_failures_md(load_manifest()),
            ])
            log.close()
            return
        rng_u = np.random.default_rng(SEED)
        Uperm = np.empty((N_PERM, n_ruler), np.float64)
        for i in range(N_PERM):
            Uperm[i] = rng_u.permutation(u)
        cos_vu, p_vu, _ = cos_against_perms(v, u, Uperm)
        beats_n1 = bool(np.isfinite(p_vu) and p_vu <= P_BAR and np.isfinite(cos_vu) and cos_vu > 0)
        log(f"[axis] cos(v, u)={_fmt(cos_vu)} p_N1={_fmt_p(p_vu)} beats_N1={beats_n1}")
        dump_json(SENG_DIR / "axes.json", jsonable(dict(
            cos_v_p=cos_vp, cos_v_u=cos_vu, p_n1_v_u=p_vu, beats_n1=beats_n1,
            v_norm=float(np.linalg.norm(v)),
        )))
        np.savez_compressed(
            SENG_DIR / "axes.npz",
            v=v, u=u, p_CRA=p_unit["CRA"], p_CRI=p_unit["CRI"],
        )

        # Identity gene positions in the ruler logcpm (upper-case ruler symbols).
        ruler_sym = np.array([str(s).upper() for s in np.asarray(frozen["symbol"])], dtype=object)
        rpos = {}
        for i, s in enumerate(ruler_sym):
            rpos.setdefault(s, i)
        id_present = [g for g in IDENTITY_GENES if g in rpos and ruler_idx[rpos[g]] >= 0]
        id_missing = [g for g in IDENTITY_GENES if g not in id_present]
        id_cols = np.array([rpos[g] for g in id_present], dtype=int)
        log(f"[identity] present={id_present} missing={id_missing}")

        # NT dense blocks for null sums (cells × present genes), float64 sums later.
        nt_dense = {}
        nt_sum = {}
        for mod in ("CRA", "CRI"):
            block = np.asarray(csr[nt_modality[mod]].toarray(), np.float64)
            nt_dense[mod] = block
            nt_sum[mod] = to_ruler_from_present(block.sum(axis=0), dst, n_ruler)
            log(f"[NT dense] {mod} shape={block.shape}")

        # Point estimates. d_g is formed here, after P and after p.
        mu_f = mu
        sd_f = sd
        miss = np.asarray(missing, bool)
        rows_out = []
        rng_null = np.random.default_rng(SEED)
        rng_boot = np.random.default_rng(BOOT_SEED)
        # Precompute official z for NT and groups from the panel (point estimate).
        for mod, lab, idx in qual:
            i = name_to_i[lab]
            i_nt = name_to_i[f"NT_{mod}"]
            d = Z[i] - Z[i_nt]
            st = stats_from_d(d, v, p_unit[mod])
            # Nulls: subset of n NT cells; origin = remaining.
            n = int(idx.size)
            n_nt = int(nt_dense[mod].shape[0])
            nulls = {k: np.empty(N_NULL) for k in ("cos", "frac", "delta", "cos_resid", "delta_resid")}
            null_ok = n < n_nt and n > 0
            if not null_ok:
                record_failure(
                    "null",
                    f"{lab}: cannot draw n={n} from NT n={n_nt} with a non-empty remainder. "
                    "p-values left NA. Not drawing with replacement.",
                )
                for k in nulls:
                    nulls[k][:] = np.nan
            else:
                for t in range(N_NULL):
                    take = rng_null.choice(n_nt, size=n, replace=False)
                    fake_present = nt_dense[mod][take].sum(axis=0)
                    fake = to_ruler_from_present(fake_present, dst, n_ruler)
                    origin = nt_sum[mod] - fake
                    z_f, _, _, _ = z_from_counts(fake, cache, mu_f, sd_f, miss)
                    z_o, _, _, _ = z_from_counts(origin, cache, mu_f, sd_f, miss)
                    stn = stats_from_d(z_f - z_o, v, p_unit[mod])
                    for k in nulls:
                        nulls[k][t] = stn[k]
            # Bootstrap: resample g and NT.
            g_block = np.asarray(csr[idx].toarray(), np.float64)
            boots = {k: np.empty(N_BOOT) for k in ("cos", "frac", "delta", "cos_resid", "delta_resid")}
            for b in range(N_BOOT):
                g_take = rng_boot.integers(0, g_block.shape[0], size=g_block.shape[0])
                nt_take = rng_boot.integers(0, n_nt, size=n_nt)
                c_g = to_ruler_from_present(g_block[g_take].sum(axis=0), dst, n_ruler)
                c_nt = to_ruler_from_present(nt_dense[mod][nt_take].sum(axis=0), dst, n_ruler)
                z_g, _, _, _ = z_from_counts(c_g, cache, mu_f, sd_f, miss)
                z_n, _, _, _ = z_from_counts(c_nt, cache, mu_f, sd_f, miss)
                stb = stats_from_d(z_g - z_n, v, p_unit[mod])
                for k in boots:
                    boots[k][b] = stb[k]
            rec = dict(
                modality=mod, perturbation=lab, n_cells=n,
                cos=st["cos"], frac=st["frac"], delta=st["delta"],
                cos_resid=st["cos_resid"], delta_resid=st["delta_resid"],
                d_norm=st["d_norm"],
            )
            for k in ("cos", "frac", "delta", "cos_resid", "delta_resid"):
                lo, hi = percentile_ci(boots[k])
                rec[f"{k}_ci_lo"] = lo
                rec[f"{k}_ci_hi"] = hi
                rec[f"{k}_ci_valid"] = bool(ci_covers_point(rec[k], lo, hi))
                if k in ("delta", "delta_resid"):
                    rec[f"p_{k}"] = permutation_p(rec[k], nulls[k], greater=False)
                else:
                    rec[f"p_{k}"] = permutation_p(rec[k], nulls[k], greater=True)
            rec["prolif_shift"] = float(np.nanmean(score[idx]) - np.nanmean(score[nt_modality[mod]]))
            cos_u, p_u, _ = cos_against_perms(d, u, Uperm)
            rec["cos_d_u"] = cos_u
            rec["p_cos_d_u"] = p_u
            if id_cols.size:
                diff = logcpm_panel[i, id_cols] - logcpm_panel[i_nt, id_cols]
                rec["identity_diff"] = float(np.mean(diff))
                rec["identity_drop"] = float(-rec["identity_diff"])
                rec["identity_n_genes"] = int(id_cols.size)
                rec["identity_loss"] = bool(rec["identity_drop"] > IDENTITY_DROP)
            else:
                rec["identity_diff"] = np.nan
                rec["identity_drop"] = np.nan
                rec["identity_n_genes"] = 0
                rec["identity_loss"] = False
            rec["ruler_score_context"] = float(ruler_ctx[i])
            rec["ruler_score_nt_context"] = float(ruler_ctx[i_nt])
            rows_out.append(rec)
            if len(rows_out) % 10 == 0:
                log(f"[screen] {len(rows_out)}/{len(qual)} last={lab} cos={_fmt(rec['cos'])} "
                    f"delta={_fmt(rec['delta'])}")

        sdf = pd.DataFrame(rows_out)
        # BH within modality, per statistic.
        for mod in ("CRA", "CRI"):
            m = sdf["modality"].to_numpy() == mod
            for stat, pcol, qcol in (
                ("cos", "p_cos", "q_cos"),
                ("frac", "p_frac", "q_frac"),
                ("delta", "p_delta", "q_delta"),
                ("cos_resid", "p_cos_resid", "q_cos_resid"),
                ("delta_resid", "p_delta_resid", "q_delta_resid"),
            ):
                if qcol not in sdf.columns:
                    sdf[qcol] = np.nan
                sdf.loc[m, qcol] = bh_q(sdf.loc[m, pcol].to_numpy())
        sdf["meets_toward_young"] = [meets_toward(r) for _, r in sdf.iterrows()]
        sdf.to_csv(SENG_DIR / "screen.csv", index=False)

        # Reproduction, then labels only if it passed.
        named_rows = []
        n_repro = 0
        for hit in NAMED_HITS:
            sub = sdf[sdf["perturbation"] == hit]
            if sub.empty:
                named_rows.append(dict(perturbation=hit, present=False, n_cells=0, cos=np.nan, q_cos=np.nan))
                record_failure("reproduction", f"{hit} is not in the n>={MIN_CELLS} screen. Not renamed.")
                continue
            r = sub.iloc[0].to_dict()
            r["present"] = True
            ok = bool(np.isfinite(r["cos"]) and r["cos"] > 0 and np.isfinite(r["q_cos"]) and r["q_cos"] <= P_BAR)
            r["repro_cos"] = ok
            if ok:
                n_repro += 1
            named_rows.append(r)
        repro_pass = n_repro >= 3
        log(f"[reproduction] {n_repro}/4 cos>0 and q<=0.05; pass={repro_pass}")
        if not repro_pass:
            record_failure(
                "pipeline_disagrees",
                f"Fewer than 3 of 4 named hits have cos_g>0 with N_NT FDR q<=0.05 ({n_repro}/4). "
                "Our pipeline does not reproduce their direction. Named hits are not interpreted as toward/away.",
            )
        labels_on = repro_pass
        for r in named_rows:
            if not r.get("present"):
                r["label"] = "absent"
                continue
            if not labels_on:
                r["label"] = "not_interpreted"
            else:
                r["label"] = classify_hit(r)
        ndf = pd.DataFrame(named_rows)
        ndf.to_csv(SENG_DIR / "named_hits.csv", index=False)

        screen_sum = {}
        flags = []
        for mod in ("CRA", "CRI"):
            sub = sdf[sdf["modality"] == mod]
            rho, rho_p = spearmanr(sub["cos"], sub["prolif_shift"]) if len(sub) >= 3 else (np.nan, np.nan)
            n_hit = int((sub["q_cos"] <= P_BAR).sum())
            n_hit_pos = int(((sub["q_cos"] <= P_BAR) & (sub["cos"] > 0)).sum())
            n_toward = int(sub["meets_toward_young"].sum())
            n_toward_among_hits = int(((sub["q_cos"] <= P_BAR) & sub["meets_toward_young"]).sum())
            screen_sum[mod] = dict(
                n=int(len(sub)), rho_prolif=float(rho) if np.isfinite(rho) else None,
                rho_p=float(rho_p) if np.isfinite(rho_p) else None,
                n_cos_q=n_hit, n_cos_q_and_positive=n_hit_pos,
                n_toward_young=n_toward, n_toward_among_cos_q=n_toward_among_hits,
            )
            log(f"[screen {mod}] n={len(sub)} rho={rho} cos_q={n_hit} toward_among={n_toward_among_hits}")
        rhos = [screen_sum[m]["rho_prolif"] for m in ("CRA", "CRI")]
        if all(r is not None and r >= RHO_PROLIF_BAR for r in rhos):
            flags.append("screen_tracks_proliferation")
        if not beats_n1:
            flags.append("passage_axis_not_donor_age")
        if not repro_pass:
            flags.append("pipeline_disagrees")
        dump_json(SENG_DIR / "screen_summary.json", jsonable(dict(
            by_modality=screen_sum, flags=flags, reproduction_pass=repro_pass,
            n_repro=n_repro, cos_v_p=cos_vp, cos_v_u=cos_vu, p_n1=p_vu, beats_n1=beats_n1,
            identity_genes_present=id_present, identity_genes_missing=id_missing,
        )))
        dump_json(SENG_DIR / "reading.json", jsonable(dict(
            flags=flags, reproduction_pass=repro_pass, n_repro=n_repro,
            named_labels={r["perturbation"]: r.get("label") for r in named_rows},
        )))

        man = load_manifest()
        man["status"] = "DONE"
        save_manifest(man)
        sections = _pre_perturb_sections(inv, split_rec, n_overlap, n_ruler, pdf, psum, man)
        sections += _post_sections(
            sdf, ndf, screen_sum, flags, repro_pass, n_repro, cos_vp, cos_vu, p_vu, beats_n1,
            id_present, id_missing, quartile_info,
        )
        write_findings("DONE.", sections)
        write_progress(
            "DONE",
            "none",
            [
                "## Reading",
                "",
                f"- reproduction_pass={repro_pass} ({n_repro}/4)",
                f"- flags={flags}",
                f"- elapsed_s={time.perf_counter() - t0:.1f}",
                "",
            ],
        )
        log(f"[done] elapsed_s={time.perf_counter() - t0:.1f} flags={flags}")
        log.close()


def _official_logcpm_row(C, nf):
    from gtex_common import log2_cpm_edger
    return log2_cpm_edger(np.asarray(C, np.float64).T, np.asarray(nf, float)).T


def to_ruler_from_present(sum_present, dst, n_ruler):
    out = np.zeros(n_ruler, np.float64)
    out[dst] = np.asarray(sum_present, np.float64).ravel()
    return out


def cos_against_perms(d, u, Uperm):
    cos = cosine_full(d, u)
    dn = float(np.linalg.norm(d))
    if dn < 1e-12 or not np.isfinite(cos):
        return cos, np.nan, None
    un = np.linalg.norm(Uperm, axis=1)
    null = (Uperm @ np.asarray(d, float).ravel()) / (un * dn)
    p = permutation_p(cos, null, greater=True)
    return float(cos), float(p), null


def ruler_column_index(gene_ids, symbols, frozen, log):
    ruler_ens = np.asarray(frozen["ensembl"]).astype(str)
    ruler_sym = np.array([str(s).upper() for s in np.asarray(frozen["symbol"])], dtype=object)
    gid = np.asarray(gene_ids).astype(str)
    ens = np.array([strip_ensembl(g) for g in gid])
    sym_u = np.array([str(s).upper() for s in np.asarray(symbols)], dtype=object)
    is_ens = np.array([e.startswith("ENS") for e in ens])
    pos_e, pos_s = {}, {}
    for i, e in enumerate(ens):
        if is_ens[i]:
            pos_e.setdefault(e, i)
        pos_s.setdefault(sym_u[i], i)
    idx = np.full(len(ruler_ens), -1, dtype=int)
    n_ens = n_sym = 0
    for j, e in enumerate(ruler_ens):
        i = pos_e.get(e)
        if i is not None:
            idx[j] = i
            n_ens += 1
            continue
        i = pos_s.get(ruler_sym[j])
        if i is not None:
            idx[j] = i
            n_sym += 1
    log(f"[align index] overlap={int((idx >= 0).sum())} ensembl={n_ens} symbol={n_sym}")
    return idx


def _write_count_tables(obs, labels, kind, pd_labels, is_wt, modality, gcol, log, guide_series=None):
    wt = (
        pd.DataFrame({"pd": pd_labels[is_wt].to_numpy()})
        .groupby("pd", dropna=False).size().rename("n_cells").reset_index()
    )
    wt.to_csv(SENG_DIR / "inventory_wt_pd.csv", index=False)
    nt = pd.DataFrame({
        "label": labels[kind == "NT"].to_numpy(),
        "modality": modality[kind.to_numpy() == "NT"],
    })
    if len(nt):
        nt_sum = nt.groupby(["modality", "label"], dropna=False).size().rename("n_cells").reset_index()
    else:
        nt_sum = pd.DataFrame(columns=["modality", "label", "n_cells"])
    nt_sum.to_csv(SENG_DIR / "inventory_nt.csv", index=False)
    pert = pd.DataFrame({
        "perturbation": labels.to_numpy(),
        "kind": kind.to_numpy(),
        "modality": modality,
    })
    pert = pert[pert["kind"].isin(["CRA", "CRI"])]
    pooled = pert.groupby(["modality", "perturbation"], dropna=False).size().rename("n_cells").reset_index()
    pooled.to_csv(SENG_DIR / "inventory_perturbations.csv", index=False)
    if guide_series is not None:
        guide_vals = np.asarray(guide_series, dtype=object)
    elif gcol is not None:
        guide_vals = obs[gcol].map(lambda x: _s(x).strip()).to_numpy()
    else:
        guide_vals = None
    if guide_vals is not None:
        g = pd.DataFrame({
            "perturbation": labels.to_numpy(),
            "guide": guide_vals,
            "modality": modality,
            "kind": kind.to_numpy(),
        })
        g = g[g["kind"].isin(["CRA", "CRI", "NT"])]
        guide_tab = g.groupby(["modality", "perturbation", "guide"], dropna=False).size().rename("n_cells").reset_index()
    else:
        guide_tab = pd.DataFrame(columns=["modality", "perturbation", "guide", "n_cells"])
    guide_tab.to_csv(SENG_DIR / "inventory_guides.csv", index=False)
    other = pd.DataFrame({"label": labels.to_numpy(), "kind": kind.to_numpy()})
    other = other[~other["kind"].isin(["CRA", "CRI", "NT", "WT"])]
    if len(other):
        oc = other.groupby("label", dropna=False).size().rename("n_cells").reset_index()
    else:
        oc = pd.DataFrame(columns=["label", "n_cells"])
    oc.to_csv(SENG_DIR / "inventory_other.csv", index=False)
    log(f"[counts] perturbations={len(pooled)} guide_rows={len(guide_tab)} wt_rows={len(wt)} nt_rows={len(nt_sum)} other={len(oc)}")


def _finish_inv(size, md5, matrix_source, raw_audit, pert_col, pd_col, wt_col, wt_note,
                gcol, gnote, nt_note, early_note, late_note, perturbed_pd_note):
    other = []
    inv = dict(
        path=str(H5_PATH), size_bytes=int(size), md5=md5,
        matrix_source=matrix_source, raw_audit=raw_audit,
        pert_col=pert_col, pd_col=pd_col, wt_col=wt_col, wt_note=wt_note,
        guide_col=gcol, guide_note=gnote, nt_note=nt_note,
        early_note=early_note, late_note=late_note, perturbed_pd_note=perturbed_pd_note,
        other_labels=other,
        wt_pd_csv=str(SENG_DIR / "inventory_wt_pd.csv"),
        nt_csv=str(SENG_DIR / "inventory_nt.csv"),
        pert_csv=str(SENG_DIR / "inventory_perturbations.csv"),
        guide_csv=str(SENG_DIR / "inventory_guides.csv"),
        other_csv=str(SENG_DIR / "inventory_other.csv"),
    )
    dump_json(SENG_DIR / "inventory.json", jsonable(inv))
    return inv


def _write_partial_inventory(obs, labels, kind, pd_labels, is_wt, modality, pert_col, pd_col,
                             wt_col, wt_note, nt_note, size, md5, matrix_source, raw_audit, log):
    if modality is None:
        modality = np.array([""] * len(obs), dtype=object)
    _write_count_tables(obs, labels, kind, pd_labels, is_wt, modality, None, log)
    _finish_inv(
        size, md5, matrix_source, raw_audit, pert_col, pd_col, wt_col, wt_note,
        None, "not resolved", nt_note, "not set", "not set", "not set",
    )


def write_inventory_stop(obs, var, raw_audit, message, log, size, md5):
    write_progress("STOP", "key=`no_raw_counts`.")
    write_findings("STOP. key=`no_raw_counts`.", [
        "## STOP / failures", "",
        *_failures_md(load_manifest()),
        f"- obs columns: {list(obs.columns)}",
        f"- var columns: {list(var.columns)}",
        f"- size={size} md5=`{md5}`",
        "",
    ])


def _pre_perturb_sections(inv, split_rec, n_overlap, n_ruler, pdf, psum, man):
    lines = [
        "## Failures (verbatim)",
        "",
        *_failures_md(man),
        *inventory_sections(inv),
        "## Gene space",
        "",
        f"Frozen-ruler overlap genes present in this file: **{n_overlap}** / {n_ruler}.",
        "Missing ruler genes are z=0 after the frozen GTEx mu/sd. "
        f"prior.count={EDGE_R_PRIOR:g}.",
        "",
        "## Split sizes (frozen before any statistic)",
        "",
        f"Seed `{SEED}`. Each of EARLY and LATE uses an independent `Generator({SEED})` "
        "(same construction as `split_half`). Half A defines v. Half B is used only in P.",
        "",
    ]
    if split_rec:
        lines.append(md_table(pd.DataFrame([split_rec])))
        lines.append("")
    lines += [
        "## Positive control P (half-B only; run before any perturbation vector)",
        "",
        "v_B = z(EARLY_B) − z(LATE_B). d_f = z(mix_f) − z(LATE_B). "
        f"Young fraction f={list(F_VALUES)}. Total cells fixed at min(n_EARLY_B, n_LATE_B). "
        f"Seed `{SEED}` for the mixture draw. Bootstrap seed `{BOOT_SEED}`, B={N_BOOT}, "
        "resamples cells inside each mixture; EARLY_B and LATE_B membership stays the full half. "
        "The P panel is TMM'd on its own. An INVALID interval does not count as excluding 0.",
        "",
    ]
    if psum:
        lines.append(
            f"- P pass={psum.get('pass_')}. monotonic={psum.get('monotonic')}. "
            f"frac(0.50) CI excludes 0 = {psum.get('frac50_ci_excludes_0')} "
            f"{_fmt_ci(psum.get('frac50_ci_lo'), psum.get('frac50_ci_hi'))} "
            f"valid={psum.get('frac50_ci_valid')}. "
            f"delta(0.50)={_fmt(psum.get('delta50'))} < 0 = {psum.get('delta50_neg')}."
        )
        lines.append(
            f"- Noise floor f=0: frac={_fmt(psum.get('frac0'))} delta={_fmt(psum.get('delta0'))}."
        )
        if psum.get("reason"):
            lines.append(f"- fail reason: {psum.get('reason')}")
        lines.append("")
    if pdf is not None and len(pdf):
        show = pdf.copy()
        for c in ("frac_ci_valid", "delta_ci_valid", "cos_ci_valid", "noise_floor"):
            if c in show.columns:
                show[c] = show[c].map(lambda x: str(bool(x)))
        keep = [c for c in (
            "f", "n_late", "n_early", "n_total", "frac", "delta", "cos",
            "frac_ci_lo", "frac_ci_hi", "frac_ci_valid",
            "delta_ci_lo", "delta_ci_hi", "delta_ci_valid", "noise_floor",
        ) if c in show.columns]
        lines.append(md_table(show[keep]))
        lines.append("")
    if psum and not psum.get("pass_"):
        lines += [
            "**P failed. key=`pos_control_broken`. No perturbation vector is reported.**",
            "",
            *limitations_md(),
        ]
    return lines


def _fmt_bool(x):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    return str(bool(x))


def _post_sections(sdf, ndf, screen_sum, flags, repro_pass, n_repro, cos_vp, cos_vu, p_vu,
                   beats_n1, id_present, id_missing, quartile_info):
    lines = [
        "## Reproduction check",
        "",
        "Named hits: CRA_E2F3, CRA_EZH2, CRI_STAT3, CRI_ZFX. "
        "A hit counts if cos_g > 0 and N_NT BH q within its modality is <= 0.05.",
        "",
        f"- {n_repro} / 4 counted. reproduction_pass={repro_pass}.",
        "",
    ]
    if not repro_pass:
        lines.append(
            "Our pipeline does not reproduce their direction. "
            "Named hits are not interpreted as toward or away."
        )
        lines.append("")
    lines += [
        "## Named hits",
        "",
        "delta assumes the TF effect is additive to the construct effect: "
        "d_g = z(g) − z(NT same modality) is applied to a LATE cell. "
        "Negative delta means that sum is closer to EARLY than LATE was. "
        "cos_g alone is not evidence of rejuvenation.",
        "",
        f"Identity genes in the panel log2-CPM: {id_present}. Missing: {id_missing}. "
        f"`identity_loss` if the mean drop (NT − g) is > {IDENTITY_DROP}.",
        "",
    ]
    if len(ndf):
        show = ndf.copy()
        cols = [c for c in (
            "perturbation", "present", "n_cells", "label", "repro_cos",
            "cos", "frac", "delta", "cos_resid", "delta_resid",
            "p_cos", "q_cos", "p_frac", "q_frac", "p_delta", "q_delta",
            "p_cos_resid", "q_cos_resid", "p_delta_resid", "q_delta_resid",
            "cos_ci_lo", "cos_ci_hi", "cos_ci_valid",
            "frac_ci_lo", "frac_ci_hi", "frac_ci_valid",
            "delta_ci_lo", "delta_ci_hi", "delta_ci_valid",
            "cos_resid_ci_lo", "cos_resid_ci_hi", "cos_resid_ci_valid",
            "delta_resid_ci_lo", "delta_resid_ci_hi", "delta_resid_ci_valid",
            "prolif_shift", "cos_d_u", "p_cos_d_u",
            "identity_diff", "identity_drop", "identity_loss",
            "ruler_score_context",
        ) if c in show.columns]
        # p/q as strings so md_table does not force +.3f
        for c in cols:
            if c.startswith("p_") or c.startswith("q_"):
                show[c] = [_fmt_p(v) for v in show[c]]
            elif c.endswith("_valid") or c in ("present", "identity_loss", "repro_cos"):
                show[c] = [_fmt_bool(v) for v in show[c]]
        lines.append(md_table(show[cols]))
        lines.append("")
        n_invalid = 0
        for c in show.columns:
            if c.endswith("_ci_valid"):
                n_invalid += int((ndf[c] == False).sum()) if c in ndf.columns else 0
        lines.append(f"Named-hit CI flags marked INVALID (valid=False): {n_invalid} column-entries.")
        lines.append("")
    lines += [
        "## Screen-level reading",
        "",
        "Modalities are not averaged. Full per-perturbation table: `results/seng/screen.csv`.",
        "",
        f"cos(v, p_CRA)={_fmt(cos_vp.get('CRA'))}. cos(v, p_CRI)={_fmt(cos_vp.get('CRI'))}.",
        f"cos(v, u)={_fmt(cos_vu)} N1 p={_fmt_p(p_vu)} beats N1={beats_n1}.",
        "",
        "Quartile sizes (NT only, cut frozen before P): "
        + ", ".join(
            f"{mod} bottom={quartile_info.get(mod, {}).get('n_bottom')} "
            f"top={quartile_info.get(mod, {}).get('n_top')} n_NT={quartile_info.get(mod, {}).get('n_nt')}"
            for mod in ("CRA", "CRI")
        )
        + ".",
        "",
    ]
    sum_rows = []
    for mod in ("CRA", "CRI"):
        s = screen_sum[mod]
        sum_rows.append(dict(
            modality=mod, n_perturbations=s["n"],
            rho_prolif=s["rho_prolif"] if s["rho_prolif"] is not None else np.nan,
            n_cos_q=s["n_cos_q"],
            n_cos_q_and_positive=s["n_cos_q_and_positive"],
            n_toward_young=s["n_toward_young"],
            n_toward_among_their_style=s["n_toward_among_cos_q"],
        ))
    lines.append(md_table(pd.DataFrame(sum_rows)))
    lines.append("")
    lines.append(
        "Their-style hits in this table are perturbations with cos_g BH q <= 0.05 "
        "(the screen bullet). The column n_cos_q_and_positive also requires cos_g > 0, "
        "which is the reproduction definition."
    )
    lines.append("")
    lines += [
        "## Fired keys",
        "",
    ]
    if flags:
        for fl in flags:
            lines.append(f"- `{fl}`")
        lines.append("")
    else:
        lines.append("None of `screen_tracks_proliferation`, `passage_axis_not_donor_age`, `pipeline_disagrees` fired.")
        lines.append("")
    if not repro_pass:
        lines.append("key=`pipeline_disagrees`.")
        lines.append("")
    if "passage_axis_not_donor_age" in flags:
        lines.append(
            "`passage_axis_not_donor_age`: aging-in-a-dish does not resemble GTEx donor aging in this gene space."
        )
        lines.append("")
    lines += limitations_md()
    lines += [
        "## Files",
        "",
        "- `src/seng_run.py`",
        "- `src/toward_run.py` and `src/same_run.py` imported, not modified",
        "- `results/seng/`",
        "- `FINDINGS_SENG.md`",
        "- `PROGRESS_SENG.md`",
        "- `data/reference/seurat_cc.genes_tirosh2016.tsv`",
        "",
    ]
    return lines


if __name__ == "__main__":
    try:
        main()
    except StopStep as e:
        msg = e.args[0] if e.args else str(e)
        step = getattr(e, "step", None) or "stop"
        # gtex StopStep signature is StopStep(step, message, details=None) via Exception args
        if len(e.args) >= 2:
            step, msg = e.args[0], e.args[1]
        record_failure(str(step), str(msg))
        write_progress("STOP", str(msg)[:300])
        write_findings(f"STOP. {step}.", [
            "## STOP / failures", "",
            *_failures_md(load_manifest()),
            *limitations_md(),
        ])
        raise

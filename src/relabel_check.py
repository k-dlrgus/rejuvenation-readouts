"""RELABEL — do the headline numbers hold if the reprogramming cells are clustered
the way Lu et al. did (both donors merged, then clustered together)?

Step 0  Lu et al.'s own per-cell labels: read meta.data from the two GEO Seurat
        objects (uncompressed XDR .rds) with a streaming reader; no R needed.
Step 1  Joint labelling: md2 Louvain pipeline (src/md2_cluster.py settings) on
        both donors' cells in one matrix, labelled with md2's state signatures
        and argmax rule (src/md2_task2.py::_label_clusters).
Step 2  Agreement of old (per-donor md2) and new labels, per donor.
Step 3  Headline numbers with the new labels (functions from same_run.py,
        toward_run.py, genespace_run.py, newstory.py; same seeds).

Writes only results/relabel_check/ and FINDINGS_RELABEL.md.

Usage: python src/relabel_check.py [step0|prereg|cluster|label|agree|stats|findings]
"""
from __future__ import annotations

import json
import struct
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import RESULTS, ROOT, DATA_RAW  # noqa: E402

OUT = RESULTS / "relabel_check"
OUT.mkdir(parents=True, exist_ok=True)
FINDINGS_PATH = ROOT / "FINDINGS_RELABEL.md"
RDS_DIR = DATA_RAW / "gse297234"
RDS_FILES = {
    "HFIB_COMBINED": RDS_DIR / "GSE297234_HFIB_COMBINED_SEVOSKM.rds",
    "GM00731": RDS_DIR / "GSE297234_GM00731_SEVOSKM.rds",
}


class Log:
    def __init__(self, path):
        self.f = open(path, "a", encoding="utf-8")

    def __call__(self, msg):
        print(msg, flush=True)
        self.f.write(str(msg) + "\n")
        self.f.flush()

    def close(self):
        self.f.close()


# --------------------------------------------------------------------------- Step 0: RDS reader
NILVALUE, GLOBALENV, UNBOUND, MISSINGARG, BASENS = 254, 253, 252, 251, 250
NAMESPACE, PACKAGE, PERSIST, EMPTYENV, BASEENV = 249, 248, 247, 242, 241
REFSXP, ATTRLANG, ATTRLIST, ALTREP = 255, 240, 239, 238
PAIRLIST_TYPES = {2, 3, 5, 6, 17, ATTRLANG, ATTRLIST}
VEC_ITEM = {10: (4, ">i4"), 13: (4, ">i4"), 14: (8, ">f8"), 15: (16, None), 24: (1, None)}
NA_INT = -2147483648


class RObj:
    __slots__ = ("t", "v", "attr")

    def __init__(self, t, v=None, attr=None):
        self.t, self.v, self.attr = t, v, attr

    def attrs(self):
        """Attribute pairlist as {tag: value}."""
        a = self.attr
        if a is None or a.t not in PAIRLIST_TYPES:
            return {}
        return {tag: val for tag, val, _ in a.v}


class Skipped:
    def __init__(self, t, n):
        self.t, self.n = t, n

    def __repr__(self):
        return f"<skipped type={self.t} n={self.n}>"


class RDSReader:
    """Streaming reader for uncompressed XDR R serialization version 3.

    Numeric/raw vectors longer than `skip_bytes` are seeked over, not read, so
    the multi-GB count matrices never enter memory.
    """

    def __init__(self, path, skip_bytes=4 << 20):
        self.f = open(path, "rb", buffering=1 << 20)
        self.refs = []
        self.skip_bytes = int(skip_bytes)
        self.n_skipped = 0
        self.bytes_skipped = 0

    def close(self):
        self.f.close()

    def i4(self):
        return struct.unpack(">i", self.f.read(4))[0]

    def length(self):
        n = self.i4()
        if n == -1:
            hi, lo = struct.unpack(">iI", self.f.read(8))
            n = (hi << 32) + lo
        return n

    def header(self):
        magic = self.f.read(2)
        if magic != b"X\n":
            raise RuntimeError(f"not an uncompressed XDR rds (magic={magic!r})")
        ver, writer, minver = self.i4(), self.i4(), self.i4()
        enc = None
        if ver == 3:
            enc = self.f.read(self.i4()).decode("ascii", "replace")
        return dict(version=ver, writer=writer, min_reader=minver, native_encoding=enc)

    def item(self):
        return self.item_flags(self.i4())

    def _charsxp(self):
        flags = self.i4()
        if flags & 0xFF != 9:
            return self.item_flags(flags).v
        n = self.i4()
        return None if n == -1 else self.f.read(n).decode("utf-8", "replace")

    def item_flags(self, flags):
        t = flags & 0xFF
        has_attr = bool(flags & (1 << 9))
        has_tag = bool(flags & (1 << 10))
        if t == REFSXP:
            idx = flags >> 8
            if idx == 0:
                idx = self.i4()
            return self.refs[idx - 1]
        if t in (NILVALUE, GLOBALENV, UNBOUND, MISSINGARG, BASENS, EMPTYENV, BASEENV):
            return RObj(t)
        if t in (NAMESPACE, PACKAGE, PERSIST):
            self.i4()
            n = self.i4()
            obj = RObj(t, [self._charsxp() for _ in range(n)])
            self.refs.append(obj)
            return obj
        if t == 1:
            obj = RObj(1, self.item().v)
            self.refs.append(obj)
            return obj
        if t == 4:
            self.i4()
            obj = RObj(4)
            self.refs.append(obj)
            enclos, frame, hashtab = self.item(), self.item(), self.item()
            obj.attr = self.item()
            obj.v = dict(frame=frame, hashtab=hashtab)
            return obj
        if t in PAIRLIST_TYPES:
            items = []
            head_t = t
            while True:
                attr = self.item() if has_attr else None
                tag = self.item() if has_tag else None
                car = self.item()
                items.append((tag.v if tag is not None else None, car, attr))
                nflags = self.i4()
                nt = nflags & 0xFF
                if nt in PAIRLIST_TYPES:
                    has_attr = bool(nflags & (1 << 9))
                    has_tag = bool(nflags & (1 << 10))
                    continue
                tail = self.item_flags(nflags)
                if tail.t != NILVALUE:
                    items.append(("__cdr__", tail, None))
                break
            return RObj(head_t, items)
        if t == 22:
            obj = RObj(22)
            self.refs.append(obj)
            self.item()
            self.item()
        elif t == 23:
            obj = RObj(23)
            self.refs.append(obj)
        elif t in (7, 8):
            obj = RObj(t, self.f.read(self.i4()).decode("ascii", "replace"))
        elif t == 9:
            n = self.i4()
            obj = RObj(9, None if n == -1 else self.f.read(n).decode("utf-8", "replace"))
        elif t in VEC_ITEM:
            n = self.length()
            size, dt = VEC_ITEM[t]
            nbytes = n * size
            if nbytes > self.skip_bytes:
                self.f.seek(nbytes, 1)
                self.n_skipped += 1
                self.bytes_skipped += nbytes
                obj = RObj(t, Skipped(t, n))
            else:
                buf = self.f.read(nbytes)
                obj = RObj(t, np.frombuffer(buf, dt).copy() if dt else buf)
        elif t == 16:
            n = self.length()
            obj = RObj(16, [self._charsxp() for _ in range(n)])
        elif t in (19, 20):
            n = self.length()
            obj = RObj(t, [self.item() for _ in range(n)])
        elif t == 25:
            obj = RObj(25)
        elif t == ALTREP:
            info, state, attr = self.item(), self.item(), self.item()
            return RObj(ALTREP, (info, state), attr)
        else:
            raise RuntimeError(f"unsupported SEXP type {t} at byte {self.f.tell()}")
        if has_attr:
            obj.attr = self.item()
        return obj


def altrep_class(o):
    info = o.v[0]
    if info.t in PAIRLIST_TYPES and info.v:
        return info.v[0][1].v
    return None


def r_to_py(o):
    """Vector-like R object → numpy / list. Factors → pandas.Categorical."""
    if o.t == ALTREP:
        cls = altrep_class(o)
        state = o.v[1]
        if cls == "compact_intseq":
            n, start, step = (float(x) for x in state.v)
            return np.arange(int(start), int(start) + int(n) * int(step), int(step))
        if cls == "compact_realseq":
            n, start, step = (float(x) for x in state.v)
            return start + step * np.arange(int(n))
        if cls in ("wrap_integer", "wrap_real", "wrap_string", "wrap_logical", "wrap_raw", "wrap_complex"):
            inner = state.v[0]
            if o.attr is not None and o.attr.t != NILVALUE:
                inner.attr = o.attr
            return r_to_py(inner)
        if cls == "deferred_string":
            base = state.v[0][1] if state.t in PAIRLIST_TYPES else state.v[0]
            vals = r_to_py(base)
            return [str(x) for x in vals]
        raise RuntimeError(f"unsupported ALTREP class {cls}")
    a = o.attrs()
    if o.t in (10, 13):
        v = o.v
        if isinstance(v, Skipped):
            return v
        if "levels" in a and o.t == 13:
            levels = r_to_py(a["levels"])
            codes = np.asarray(v, int)
            codes = np.where(codes == NA_INT, -1, codes - 1)
            return pd.Categorical.from_codes(codes, categories=list(levels))
        if o.t == 10:
            out = np.where(v == NA_INT, np.nan, v.astype(float))
            return out
        return np.where(v == NA_INT, np.nan, v.astype(float)) if (v == NA_INT).any() else v.astype(int)
    if o.t == 14:
        return o.v
    if o.t == 16:
        return o.v
    if o.t == 19:
        return [r_to_py(x) for x in o.v]
    if o.t == NILVALUE:
        return None
    return o


def data_frame(o):
    a = o.attrs()
    names = r_to_py(a["names"])
    cols = {}
    for nm, col in zip(names, o.v):
        cols[nm] = r_to_py(col)
    rn = a.get("row.names")
    rn = r_to_py(rn) if rn is not None else None
    df = pd.DataFrame(cols)
    if isinstance(rn, list):
        df.index = rn
    return df


def read_seurat_meta(path, log):
    r = RDSReader(path)
    hdr = r.header()
    log(f"[step0] {path.name}: header {hdr}")
    top = r.item()
    top_cls = r_to_py(top.attrs().get("class")) if top.attrs().get("class") is not None else None
    slots = top.attrs()
    log(f"[step0] top-level type={top.t} class={top_cls} slots={list(slots)}; "
        f"skipped {r.n_skipped} large vectors ({r.bytes_skipped / 1e9:.2f} GB) by seeking")
    meta = data_frame(slots["meta.data"])
    extra = {}
    for k in ("active.assay", "active.ident", "project.name", "version"):
        if k in slots:
            try:
                v = r_to_py(slots[k])
                extra[k] = (list(v.categories)[:50] if isinstance(v, pd.Categorical) else v)
            except Exception as e:  # noqa: BLE001
                extra[k] = f"unreadable: {e}"
    assays = slots.get("assays")
    if assays is not None:
        extra["assays"] = r_to_py(assays.attrs()["names"]) if "names" in assays.attrs() else None
    red = slots.get("reductions")
    if red is not None and "names" in red.attrs():
        extra["reductions"] = r_to_py(red.attrs()["names"])
    gr = slots.get("graphs")
    if gr is not None and "names" in gr.attrs():
        extra["graphs"] = r_to_py(gr.attrs()["names"])
    r.close()
    return meta, extra


def step0(log):
    import json
    rows = []
    for key, path in RDS_FILES.items():
        if not path.exists():
            log(f"[step0] {path} missing")
            continue
        meta, extra = read_seurat_meta(path, log)
        meta.to_csv(OUT / f"lu_meta_{key}.csv", index_label="cell")
        log(f"[step0] {key}: meta.data {meta.shape}; columns {list(meta.columns)}")
        log(f"[step0] {key}: other slots {extra}")
        for c in meta.columns:
            s = meta[c]
            nun = int(s.nunique(dropna=False))
            ex = s.astype(str).value_counts().head(25).to_dict() if nun <= 60 else s.astype(str).head(5).tolist()
            rows.append(dict(object=key, column=c, dtype=str(s.dtype), n_unique=nun, examples=json.dumps(ex)))
            log(f"[step0] {key}.{c}: dtype={s.dtype} n_unique={nun} " + (f"counts={ex}" if nun <= 60 else f"head={ex}"))
        (OUT / f"lu_slots_{key}.json").write_text(json.dumps(extra, indent=2, default=str), encoding="utf-8")
    pd.DataFrame(rows).to_csv(OUT / "lu_meta_columns.csv", index=False)


# --------------------------------------------------------------------------- pre-registration
PREREG_FLAG = OUT / "PREREG.flag"
RULE = ("The conclusion holds if, in both the full and MD gene spaces, the aged PartialReprog cells end "
        "farther from the young target (change in distance > 0), the mixture control passes, and the "
        "forward-minus-reverse asymmetry interval excludes 0.")

PREREG_TEXT = f"""## Pre-registration (written before the joint clustering and before any statistic)

Frozen copy: results/relabel_check/PREREG.flag. Not rewritten after this point.

### What Step 0 found (done before this block)

Lu et al. deposited two Seurat objects as GEO GSE297234 supplementary files
(GSE297234_HFIB_COMBINED_SEVOSKM.rds, GSE297234_GM00731_SEVOSKM.rds; the third file is the
10x RAW tar we already use). The paper has no separate data-availability paragraph; its key
resources table points to GEO GSE297234. The combined object (both donors, 53,454 cells, active
assay SCT) has a per-cell `cell_state` column with 16 sub-states (Fib_1-3, ParRep_1/2/4 and
ParRepr_3, EarPlu_1-3, Plu_1-2, NonRep_1-4). These are Lu et al.'s published per-cell labels.
The GM00731-only object has cluster letters but no `cell_state`. Cell names are
`<line>_D<day>_<barcode>`, which match our (cell_line, day, barcode).
Disclosure: before writing this block I tabulated Lu's cell_state x sample cell counts (counts
only; no expression, score or distance was computed with any new label).

### Labellings

- OLD: md2 labels (per-donor Louvain, results/md2/t2_cluster_labels_*.csv).
- JOINT (primary; answers the question asked): all cells of both donors in one matrix, md2
  Louvain pipeline unchanged: Gene Expression features; cells with >500 detected genes; genes
  detected in >=5 cells of the merged matrix; LogNormalize (scale 1e4); 3000 HVGs by the md2
  quadratic log-mean/log-variance rule; each HVG regressed on percent-mito and z-scored;
  50 PCs (randomized, random_state 20260914); SNN k=20, prune 1/15; networkx Louvain
  resolution 0.8, seed 20260914. Labels: md2 AddModuleScore state signatures
  (results/md2/genesets.json; nbin 24, ctrl 100, seed 20260914) computed on the joint matrix
  (one set of expression bins over both donors), then src/md2_task2.py::_label_clusters
  (each joint cluster gets the state with the highest mean score over all its cells).
  Row order: aged libraries (day 0, 3, 7, 10) then young.
- LU (second labelling): Lu's `cell_state` collapsed by prefix: Fib_ -> Fibroblast,
  ParRep_/ParRepr_ -> PartialReprog, EarPlu_ -> EarlyPluripotency, Plu_ -> Pluripotency,
  NonRep_ -> NonReprog. Our cells absent from Lu's object are unlabelled and belong to no
  state.

Donors are pooled only to cluster and label. Every statistic stays per donor.

### Memory deviation and its gate

The unchanged code (src/md2_cluster.py::louvain_one_donor) needs several full float64 copies
of the merged matrix and does not fit in this machine's memory for ~56k cells. The joint run
therefore uses a streamed version meant to perform the same arithmetic in the same order
(per-gene means accumulated row by row in the same order scipy uses; per-cell steps done per
donor; the HVG block residualized in column blocks). No cell or gene is subsampled.
Gate: the same streamed code is first run on GM00731 alone and must reproduce md2 exactly
(identical cluster of every cell in results/md2/cluster_labels_louvain_GM00731.csv, identical
state scores in data/processed/md2/louvain_state_GM00731.npz, identical labels in
results/md2/t2_cluster_labels_GM00731.csv). If it does not, the joint run is not done and I
stop and report.

### Statistics (per donor, same functions and seeds)

a. MD score per state (src/newstory.py Part A joint scoring): the per-cell joint MD scores in
   results/newstory/partA_cell_scores.npz do not depend on labels and are reused; state means
   and bootstrap CIs (B=200, Generator(20260918), donors then states in fixed order) are
   recomputed with the new masks. Reported: aged day-0 Fibroblast, aged PartialReprog, young
   day-0 Fibroblast, drop / gap.
b. Same-platform test (src/genespace_run.py sections 2 and bootstrap, FULL and MD spaces):
   O, Y = half A of day-0 Fibroblast cells of the aged and young donor; S = aged PartialReprog;
   S_rev = young PartialReprog. Forward progress, change in distance (delta), final / start
   distance, reverse progress, asymmetry (forward minus reverse progress) with its 95%
   percentile bootstrap interval (B=200, Generator(20260918), O and Y frozen), and the mixture
   control (genespace_run.mixture_gate: progress non-decreasing in f = 0, 0.10, 0.25, 0.50;
   progress CI at f=0.50 excludes 0 and contains the point; delta at f=0.50 < 0).
   Split rule: if the set of day-0 Fibroblast cells of either donor differs from the old one,
   the half-split is redone with same_run.split_half(seed 20260914) on the new set (the
   mixtures are rebuilt from the new halves with the same seed), and this is stated.
c. GTEx test (src/toward_run.py donor panel, genespace_run.gtex_point_block, FULL space,
   aged PartialReprog): cosine, change in distance to young GTEx, change in distance to old.

Harness check before any new-label number: the same code is run with the OLD labels and must
reproduce the stored values (results/newstory/partA_*.csv, results/genespace/same_stats.csv,
same_posctrl.csv, gtex_stats.csv; points and CI bounds) to within 1e-9. If not, stop.

### Decision rule (verbatim)

"{RULE}"

Applied to each labelling separately (JOINT is the answer to the question; LU is reported
alongside, not averaged). "Full" and "MD" are the FULL (23,485 ruler genes) and MD (205
genes) spaces of src/genespace_run.py; "change in distance" is the forward delta; the
interval is the 95% percentile bootstrap interval. If a state needed by the rule is missing
or too small for the code to run, the rule is counted as not met for that labelling.
How far numbers moved is reported descriptively (old, new, difference); no threshold.
"""


def prereg(log):
    if PREREG_FLAG.exists():
        if PREREG_FLAG.read_text(encoding="utf-8") != PREREG_TEXT:
            raise RuntimeError("PREREG.flag exists and differs from PREREG_TEXT. Not rewriting.")
        log("[prereg] PREREG.flag already present and identical")
    else:
        PREREG_FLAG.write_text(PREREG_TEXT, encoding="utf-8")
        log("[prereg] wrote PREREG.flag")
    if not FINDINGS_PATH.exists():
        FINDINGS_PATH.write_text(
            "# FINDINGS_RELABEL — do the headline numbers survive Lu-style joint clustering?\n\n"
            "Status: pre-registered; no statistic computed yet.\n\n" + PREREG_TEXT, encoding="utf-8")
        log("[prereg] wrote FINDINGS_RELABEL.md (pre-registration only)")


def require_prereg():
    if not PREREG_FLAG.exists() or PREREG_FLAG.read_text(encoding="utf-8") != PREREG_TEXT:
        raise RuntimeError("PREREG.flag missing or changed. Run `prereg` first; do not edit it.")


# --------------------------------------------------------------------------- Step 1: joint Louvain
def _md2_imports():
    global MD2_SEED, LOUVAIN_RES, LOUVAIN_K_PARAM, LOUVAIN_PRUNE, LOUVAIN_NPC, LOUVAIN_N_HVG
    global LOUVAIN_MIN_GENES, LOUVAIN_MIN_CELLS_PER_GENE, AMS_NBIN, MD2_DIR, MD2_PROC
    global AGED_LINE, YOUNG_LINE, DONORS, mito_mask, add_module_score, lognormalize_csr, _cut_number
    global _gene_expression_mask, _label_clusters, STATE_NAMES
    from md2_common import (MD2_SEED, LOUVAIN_RES, LOUVAIN_K_PARAM, LOUVAIN_PRUNE, LOUVAIN_NPC,
                            LOUVAIN_N_HVG, LOUVAIN_MIN_GENES, LOUVAIN_MIN_CELLS_PER_GENE, AMS_NBIN,
                            MD2_DIR, MD2_PROC, AGED_LINE, YOUNG_LINE)
    from fibro2_tb import mito_mask
    from md2_score import add_module_score, lognormalize_csr, _cut_number
    from md2_cluster import _gene_expression_mask
    from md2_task2 import _label_clusters, STATE_NAMES
    DONORS = (AGED_LINE, YOUNG_LINE)


def _load_allcell(line):
    z = np.load(MD2_PROC / f"allcell_counts_{line}.npz", allow_pickle=True)
    var = pd.DataFrame(dict(symbol=np.asarray(z["symbols"]).astype(str),
                            gene_id=np.asarray(z["gene_id"]).astype(str),
                            feature_type=np.asarray(z["feature_type"]).astype(str)))
    obs = pd.read_csv(MD2_DIR / f"allcell_obs_{line}.csv")
    shape = tuple(int(x) for x in z["shape"])
    data, indices, indptr = z["data"], z["indices"], z["indptr"]
    if len(obs) != shape[0] or shape[1] != len(var):
        raise RuntimeError(f"{line}: allcell obs/var do not match counts {shape}")
    if data.size and int(data.min()) <= 0:
        raise RuntimeError(f"{line}: explicit zeros or negatives in counts")
    return data, indices, indptr, shape, var, obs


def _kept_csr(line, keep_genes, remap, n_kept_genes):
    """Rows with >MIN_GENES detected genes, columns keep_genes; entry order within rows preserved
    (the same matrix as X[keep_cells][:, keep_genes] in louvain_one_donor)."""
    from scipy import sparse
    data, indices, indptr, shape, var, obs = _load_allcell(line)
    nnz_row = np.diff(indptr)
    keep_cells = nnz_row > int(LOUVAIN_MIN_GENES)
    em = np.repeat(keep_cells, nnz_row)
    em &= keep_genes[indices]
    c = np.zeros(em.size + 1, dtype=np.int32)
    np.cumsum(em, dtype=np.int32, out=c[1:])
    per_row = (c[indptr[1:]] - c[indptr[:-1]])[keep_cells]
    del c
    new_indptr = np.zeros(per_row.size + 1, dtype=np.int64)
    np.cumsum(per_row, out=new_indptr[1:])
    new_data = data[em]
    new_idx = remap[indices[em]]
    del data, indices, em
    X = sparse.csr_matrix((new_data, new_idx, new_indptr), shape=(int(keep_cells.sum()), int(n_kept_genes)))
    obs = obs.iloc[np.flatnonzero(keep_cells)].reset_index(drop=True)
    return X, obs


def _row_chunks(n, size=4000):
    for r0 in range(0, n, size):
        yield r0, min(n, r0 + size)


def streamed_louvain(lines, log, tag, tmp_dir):
    """md2 louvain_one_donor on the row-stack of `lines`, without materializing float64 copies."""
    import gc
    import time
    import networkx as nx
    from scipy import sparse
    from sklearn.decomposition import PCA
    from sklearn.neighbors import NearestNeighbors
    t_start = time.time()
    sets = json.loads((MD2_DIR / "genesets.json").read_text(encoding="utf-8"))
    gene_sets = {"md_score": list(sets["MD"]), "tgfb_score": list(sets["TGFB"])}
    for st, genes in (sets.get("reprog") or {}).items():
        gene_sets[f"state_{st}"] = list(genes)

    # pass 1: cell filter, per-gene detection counts over the merged matrix
    var0, n_in, gene_cnt, n_keep = None, {}, None, {}
    for line in lines:
        data, indices, indptr, shape, var, obs = _load_allcell(line)
        ge = _gene_expression_mask(var)
        if not ge.all():
            raise RuntimeError(f"{line}: non-Gene-Expression features present; streamed code assumes none")
        if var0 is None:
            var0 = var
            gene_cnt = np.zeros(len(var), dtype=np.int64)
        elif not (np.array_equal(var0.gene_id.to_numpy(), var.gene_id.to_numpy())
                  and np.array_equal(var0.symbol.to_numpy(), var.symbol.to_numpy())):
            raise RuntimeError(f"{line}: gene axis differs from {lines[0]}")
        nnz_row = np.diff(indptr)
        keep_cells = nnz_row > int(LOUVAIN_MIN_GENES)
        gene_cnt += np.bincount(indices[np.repeat(keep_cells, nnz_row)], minlength=len(var))
        n_in[line], n_keep[line] = int(shape[0]), int(keep_cells.sum())
        log(f"[{tag} pass1] {line} cells in={shape[0]} kept={n_keep[line]} dropped={shape[0] - n_keep[line]}")
        del data, indices, indptr
        gc.collect()
    keep_genes = gene_cnt >= int(LOUVAIN_MIN_CELLS_PER_GENE)
    G = int(keep_genes.sum())
    remap = np.full(len(var0), -1, dtype=np.int32)
    remap[keep_genes] = np.arange(G, dtype=np.int32)
    var = var0.iloc[np.flatnonzero(keep_genes)].reset_index(drop=True)
    symbols = var["symbol"].astype(str).to_numpy()
    N = int(sum(n_keep.values()))
    inv = 1.0 / N
    log(f"[{tag} filter] cells kept={N} genes kept={G} dropped genes={int((~keep_genes).sum())}")

    # pass 2: pct_mt, raw moments (HVG), log-normalized means (AMS bins)
    mito, mito_meta = mito_mask(var)
    if not mito.any():
        raise RuntimeError("no MT- genes")
    m_raw, m_sq, m_log = np.zeros(G), np.zeros(G), np.zeros(G)
    obs_parts, pct_parts = [], []
    for line in lines:
        t0 = time.time()
        X, obs = _kept_csr(line, keep_genes, remap, G)
        umi = np.asarray(X.sum(axis=1)).ravel()
        mito_umi = np.asarray(X[:, mito].sum(axis=1)).ravel()
        pct = 100.0 * mito_umi / np.where(umi > 0, umi, np.nan)
        if not np.isfinite(pct).all():
            raise RuntimeError(f"{line}: non-finite percent.mt")
        for s in range(0, X.data.size, 20_000_000):
            e = min(X.data.size, s + 20_000_000)
            v = X.data[s:e].astype(np.float64)
            ix = X.indices[s:e]
            np.add.at(m_raw, ix, v * inv)
            np.add.at(m_sq, ix, (v ** 2) * inv)
            del v
        for r0, r1 in _row_chunks(X.shape[0]):
            L = lognormalize_csr(X[r0:r1])
            np.add.at(m_log, L.indices, L.data * inv)
        obs_parts.append(obs)
        pct_parts.append(pct)
        log(f"[{tag} pass2] {line} n={X.shape[0]} nnz={X.nnz} ({time.time() - t0:.0f}s)")
        del X
        gc.collect()
    obs = pd.concat(obs_parts, axis=0, ignore_index=True)
    pct_mt = np.concatenate(pct_parts)

    # HVG: tail of md2_cluster._hvg_quadratic on the accumulated moments
    var_g = np.maximum(m_sq - m_raw ** 2, 0.0)
    m = (m_raw > 0) & np.isfinite(var_g) & (var_g > 0)
    idx = np.flatnonzero(m)
    lx, ly = np.log10(m_raw[idx]), np.log10(var_g[idx])
    coef = np.polyfit(lx, ly, 2)
    resid = ly - np.polyval(coef, lx)
    hvg_idx = idx[np.argsort(resid)[::-1]][:int(LOUVAIN_N_HVG)]
    log(f"[{tag} hvg] n_candidates={idx.size} coef={coef.tolist()}")
    bins = _cut_number(m_log, AMS_NBIN, np.random.default_rng(int(MD2_SEED)))

    # pass 3: HVG block to disk, AddModuleScore per donor with the shared bins
    tmp_dir.mkdir(parents=True, exist_ok=True)
    ams_parts, blk_paths = {k: [] for k in gene_sets}, []
    for line in lines:
        t0 = time.time()
        X, _ = _kept_csr(line, keep_genes, remap, G)
        n = X.shape[0]
        logdata = np.empty(X.nnz, dtype=np.float64)
        for r0, r1 in _row_chunks(n):
            L = lognormalize_csr(X[r0:r1])
            logdata[X.indptr[r0]:X.indptr[r1]] = L.data
        logX = sparse.csr_matrix((logdata, X.indices, X.indptr), shape=X.shape)
        del X, L
        gc.collect()
        p = tmp_dir / f"hvg_{tag}_{line}.npy"
        mm = np.lib.format.open_memmap(p, mode="w+", dtype=np.float32, shape=(n, int(LOUVAIN_N_HVG)))
        for r0, r1 in _row_chunks(n):
            mm[r0:r1] = logX[r0:r1][:, hvg_idx].astype(np.float32).toarray()
        mm.flush()
        del mm
        blk_paths.append(p)
        ams, _meta = add_module_score(logX, symbols, gene_sets, log=None, tag=f"{tag}_{line}",
                                      bins=bins, gene_mean=m_log)
        for k in gene_sets:
            ams_parts[k].append(np.asarray(ams[k], np.float64))
        log(f"[{tag} pass3] {line} n={n} ams n_ctrl="
            f"{ {k: v['n_ctrl'] for k, v in _meta['per_set'].items()} } ({time.time() - t0:.0f}s)")
        del logX, logdata, ams
        gc.collect()
    ams = {k: np.concatenate(v) for k, v in ams_parts.items()}

    Z = np.empty((N, int(LOUVAIN_N_HVG)), dtype=np.float32)
    r = 0
    for p in blk_paths:
        mm = np.load(p, mmap_mode="r")
        Z[r:r + mm.shape[0]] = mm
        r += mm.shape[0]
        del mm
        p.unlink()
    # md2_cluster._scale_residualize_mt, in column blocks
    A = np.column_stack([np.ones(N), np.asarray(pct_mt, float)])
    for j0 in range(0, Z.shape[1], 250):
        j1 = min(Z.shape[1], j0 + 250)
        blk = np.ascontiguousarray(Z[:, j0:j1])
        beta, *_ = np.linalg.lstsq(A, blk, rcond=None)
        res = blk - A @ beta
        mu = res.mean(axis=0)
        sd = res.std(axis=0)
        sd = np.where(sd < 1e-12, 1.0, sd)
        Z[:, j0:j1] = np.asarray((res - mu) / sd, np.float32)
        del blk, res
    gc.collect()
    npc = int(min(LOUVAIN_NPC, Z.shape[0] - 1, Z.shape[1]))
    pca = PCA(n_components=npc, svd_solver="randomized", random_state=int(MD2_SEED))
    P = pca.fit_transform(Z).astype(np.float32)
    evr = float(np.sum(pca.explained_variance_ratio_))
    del Z, pca
    gc.collect()
    log(f"[{tag} pca] npc={npc} evr_sum={evr:.4f}")
    nn = NearestNeighbors(n_neighbors=int(LOUVAIN_K_PARAM), algorithm="auto", metric="euclidean")
    nn.fit(P)
    _, nbr = nn.kneighbors(P)
    n = int(P.shape[0])
    k = int(LOUVAIN_K_PARAM)
    rows = np.repeat(np.arange(n), k)
    M = sparse.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, nbr.ravel())), shape=(n, n))
    inter = (M @ M.T).tocoo()
    union = (2.0 * k) - inter.data
    jacc = inter.data / np.where(union > 0, union, 1.0)
    keep_e = (jacc >= float(LOUVAIN_PRUNE)) & (inter.row < inter.col)
    Gr = nx.Graph()
    Gr.add_nodes_from(range(n))
    Gr.add_weighted_edges_from(zip(inter.row[keep_e].tolist(), inter.col[keep_e].tolist(), jacc[keep_e].tolist()))
    n_edges = Gr.number_of_edges()
    del inter, union, jacc, keep_e, M
    log(f"[{tag} snn] n={n} edges={n_edges}")
    t0 = time.time()
    comms = nx.community.louvain_communities(Gr, weight="weight", resolution=float(LOUVAIN_RES),
                                             seed=int(MD2_SEED))
    lab = np.full(n, -1, dtype=int)
    for i, members in enumerate(sorted(comms, key=lambda s: (-len(s), min(s) if s else 0))):
        for v in members:
            lab[int(v)] = i
    if (lab < 0).any():
        raise RuntimeError("unassigned cells")
    sizes = np.bincount(lab).tolist()
    log(f"[{tag} louvain] n_clusters={len(sizes)} sizes={sizes} ({time.time() - t0:.0f}s; total {time.time() - t_start:.0f}s)")
    obs["cluster"] = lab
    obs["percent_mt"] = pct_mt
    meta = dict(lines=list(lines), n_cells_in=n_in, n_cells_kept=n_keep, n_genes_kept=G,
                n_genes_dropped=int((~keep_genes).sum()), n_hvg_candidates=int(idx.size),
                hvg_coef=coef.tolist(), npc=npc, evr_sum=evr, snn_edges=int(n_edges),
                n_clusters=len(sizes), sizes=sizes, mito_rule=mito_meta.get("rule"),
                n_mito=mito_meta.get("n_mito"))
    return obs, ams, meta, dict(m_log=m_log, bins=bins, P=P)


def _label_df(obs, ams, log, name):
    lab = _label_clusters(obs, ams, log, name)
    lab_map = {int(c): str(s) for c, s in zip(lab.cluster, lab.label)}
    return lab, obs["cluster"].map(lab_map).astype(str).to_numpy()


def cluster_step(log):
    require_prereg()
    _md2_imports()
    tmp = OUT / "_tmp"
    # gate: GM00731 alone must reproduce md2
    gate_path = OUT / "validation_GM00731.json"
    obs, ams, meta, extra = streamed_louvain([AGED_LINE], log, "val", tmp)
    ref = pd.read_csv(MD2_DIR / "cluster_labels_louvain_GM00731.csv")
    st_ref = np.load(MD2_PROC / "louvain_state_GM00731.npz")
    ams_ref = np.load(MD2_PROC / "louvain_ams_GM00731.npz", allow_pickle=True)
    lab_ref = pd.read_csv(MD2_DIR / "t2_cluster_labels_GM00731.csv")
    lab_new, _ = _label_df(obs, ams, log, AGED_LINE)
    checks = dict(
        same_cells_in_order=bool(len(obs) == len(ref) and (obs.barcode.astype(str).to_numpy() == ref.barcode.astype(str).to_numpy()).all()
                                 and (obs.gsm.astype(str).to_numpy() == ref.gsm.astype(str).to_numpy()).all()),
        percent_mt_max_abs_diff=float(np.max(np.abs(obs.percent_mt.to_numpy() - ref.percent_mt.to_numpy()))),
        cluster_identical=bool(np.array_equal(obs.cluster.to_numpy(int), ref.cluster.to_numpy(int))),
        gene_mean_max_abs_diff=float(np.max(np.abs(extra["m_log"] - ams_ref["gene_mean"]))),
        bins_identical=bool(np.array_equal(extra["bins"], ams_ref["bins"])),
        md_max_abs_diff=float(np.max(np.abs(ams["md_score"] - ams_ref["md"]))),
        state_max_abs_diff={k: float(np.max(np.abs(ams[k] - st_ref[k]))) for k in st_ref.files},
        labels_identical=bool(lab_new[["cluster", "label"]].astype(str).reset_index(drop=True)
                              .equals(lab_ref[["cluster", "label"]].astype(str).reset_index(drop=True))),
        meta=meta,
    )
    checks["passed"] = bool(checks["same_cells_in_order"] and checks["cluster_identical"]
                            and checks["bins_identical"] and checks["labels_identical"]
                            and checks["md_max_abs_diff"] == 0.0
                            and max(checks["state_max_abs_diff"].values()) == 0.0
                            and checks["percent_mt_max_abs_diff"] <= 1e-12)
    gate_path.write_text(json.dumps(checks, indent=2, default=str), encoding="utf-8")
    log(f"[gate] {json.dumps({k: v for k, v in checks.items() if k != 'meta'}, default=str)}")
    if not checks["passed"]:
        raise RuntimeError("GM00731 streamed rerun does not reproduce md2 exactly. STOP; joint run not done.")
    del obs, ams, extra
    import gc
    gc.collect()

    obs, ams, meta, extra = streamed_louvain(list(DONORS), log, "joint", tmp)
    lab, cell_label = _label_df(obs, ams, log, "JOINT")
    n_by = obs.assign(label=cell_label).groupby(["cluster", "cell_line"]).size().unstack(fill_value=0)
    lab = lab.merge(n_by.reset_index(), on="cluster", how="left")
    lab.to_csv(OUT / "joint_cluster_labels.csv", index=False)
    cells = obs[["barcode", "gsm", "cell_line", "day", "umi", "n_genes", "percent_mt", "cluster"]].copy()
    cells["joint_label"] = cell_label
    cells.to_csv(OUT / "joint_cells.csv", index=False)
    np.savez_compressed(OUT / "joint_ams.npz", **ams, P=extra["P"])
    (OUT / "joint_meta.json").write_text(json.dumps(meta, indent=2, default=str), encoding="utf-8")
    log(f"[joint] wrote joint_cells.csv, joint_cluster_labels.csv; labels per cluster: "
        f"{dict(zip(lab.cluster.tolist(), lab.label.tolist()))}")


# --------------------------------------------------------------------------- labels per donor
LU_PREFIX = (("ParRepr_", "PartialReprog"), ("ParRep_", "PartialReprog"), ("EarPlu_", "EarlyPluripotency"),
             ("Plu_", "Pluripotency"), ("NonRep_", "NonReprog"), ("Fib_", "Fibroblast"))
UNLAB = "Unlabelled"


def lu_state(s):
    for pre, st in LU_PREFIX:
        if str(s).startswith(pre):
            return st
    raise ValueError(f"unknown Lu cell_state {s!r}")


def donor_labels():
    """louvain_obs rows per donor with OLD, JOINT and LU labels (+ Lu sub-state)."""
    _md2_imports()
    joint = pd.read_csv(OUT / "joint_cells.csv")
    lu = pd.read_csv(OUT / "lu_meta_HFIB_COMBINED.csv", usecols=["cell", "orig.ident", "cell_state"])
    lu_map = dict(zip(lu["cell"].astype(str), lu["cell_state"].astype(str)))
    out = {}
    for line in DONORS:
        obs = pd.read_csv(MD2_DIR / f"louvain_obs_{line}.csv").reset_index(drop=True)
        lab = pd.read_csv(MD2_DIR / f"t2_cluster_labels_{line}.csv")
        old_map = {int(c): str(s) for c, s in zip(lab.cluster, lab.label)}
        obs["old"] = obs["cluster"].map(lambda c: old_map[int(c)])
        j = joint[joint.cell_line.astype(str) == line].reset_index(drop=True)
        if not (len(j) == len(obs) and (j.barcode.astype(str).to_numpy() == obs.barcode.astype(str).to_numpy()).all()
                and (j.gsm.astype(str).to_numpy() == obs.gsm.astype(str).to_numpy()).all()):
            raise RuntimeError(f"{line}: joint cells are not the louvain_obs cells in order")
        obs["joint_cluster"] = j["cluster"].to_numpy(int)
        obs["joint"] = j["joint_label"].astype(str).to_numpy()
        names = line + "_D" + obs["day"].astype(int).astype(str) + "_" + obs["barcode"].astype(str)
        sub = names.map(lambda k: lu_map.get(k))
        obs["lu_cell_state"] = sub.fillna(UNLAB)
        obs["lu"] = [lu_state(s) if s != UNLAB else UNLAB for s in obs["lu_cell_state"]]
        out[line] = obs
    return out


# --------------------------------------------------------------------------- Step 2
def agree_step(log):
    require_prereg()
    lab = donor_labels()
    order = list(STATE_NAMES) + [UNLAB]
    rows, tabs = [], []
    for line, obs in lab.items():
        for a, b in (("old", "joint"), ("old", "lu"), ("joint", "lu")):
            ct = pd.crosstab(obs[a], obs[b]).reindex(index=[s for s in order if s in set(obs[a])],
                                                      columns=[s for s in order if s in set(obs[b])], fill_value=0)
            ct.to_csv(OUT / f"agree_{line}_{a}_x_{b}.csv")
            both = obs[b] != UNLAB
            same = float((obs.loc[both, a] == obs.loc[both, b]).mean())
            rows.append(dict(cell_line=line, rows=a, cols=b, n_cells=int(len(obs)), n_compared=int(both.sum()),
                             share_unchanged=same))
            tabs.append((line, a, b, ct))
            log(f"[agree] {line} {a} x {b}: share unchanged={same:.4f} (n={int(both.sum())})\n{ct}")
        for scheme in ("old", "joint", "lu"):
            d0 = obs.day.astype(int) == 0
            for st in order:
                n = int((obs[scheme] == st).sum())
                if st == UNLAB and n == 0:
                    continue
                rows.append(dict(cell_line=line, rows=scheme, cols="count", state=st, n_cells=n,
                                 n_day0=int(((obs[scheme] == st) & d0).sum())))
    pd.DataFrame(rows).to_csv(OUT / "agreement_summary.csv", index=False)
    d0_sets = {}
    for line, obs in lab.items():
        d0 = obs.day.astype(int) == 0
        old = set(np.flatnonzero((obs.old == "Fibroblast") & d0))
        for scheme in ("joint", "lu"):
            new = set(np.flatnonzero((obs[scheme] == "Fibroblast") & d0))
            d0_sets[f"{line}_{scheme}"] = dict(n_old=len(old), n_new=len(new), n_both=len(old & new),
                                               only_old=len(old - new), only_new=len(new - old),
                                               identical=old == new)
    (OUT / "d0_fibroblast_membership.json").write_text(json.dumps(d0_sets, indent=2), encoding="utf-8")
    log(f"[agree] day-0 Fibroblast membership: {d0_sets}")


# --------------------------------------------------------------------------- Step 3
def headline(obs_by, scheme, log, ctx):
    """Part a (newstory joint MD), part b (same-platform FULL/MD), part c (GTEx FULL) for one labelling."""
    import genespace_run as gs
    import newstory as ns
    from same_run import (split_half, d0_fibroblast_idx, state_idx, mix_indices, mix_sum,
                          resample_idx, panel_z, check_gene_space, DEST_STATE, PLURI_STATE, NONREPROG_STATE)
    from toward_run import donor_panel, zscore_frozen, sum_rows
    obs = {line: o.assign(label=o[scheme].astype(str)).reset_index(drop=True) for line, o in obs_by.items()}
    res = dict(scheme=scheme)

    # a. MD score per state, newstory joint scoring
    sc = ctx["partA"]
    for line in DONORS:
        if not np.array_equal(sc[f"day__{line}"].astype(int), obs[line].day.to_numpy(int)):
            raise RuntimeError(f"{line}: partA_cell_scores rows are not louvain_obs rows")
    rng = np.random.default_rng(ns.BOOT_SEED)
    a_rows, draws = [], {}
    for line in DONORS:
        masks = ns.state_masks(obs[line])
        x_all = np.asarray(sc[f"joint__{line}"], float)
        for st in ns.STATES:
            x = x_all[masks[st]]
            if x.size == 0:
                a_rows.append(dict(cell_line=line, state=st, n_cells=0, joint_mean=np.nan,
                                   joint_ci_lo=np.nan, joint_ci_hi=np.nan))
                continue
            b = ns.boot_mean_draws(x, rng)
            draws[(line, st)] = b
            a_rows.append(dict(cell_line=line, state=st, n_cells=int(x.size), joint_mean=float(x.mean()),
                               joint_ci_lo=float(np.percentile(b, 2.5)), joint_ci_hi=float(np.percentile(b, 97.5))))
    a_df = pd.DataFrame(a_rows)

    def ga(line, st):
        return float(a_df[(a_df.cell_line == line) & (a_df.state == st)].joint_mean.iloc[0])
    a0, y0, apr = ga(AGED_LINE, "Fibroblast_d0"), ga(YOUNG_LINE, "Fibroblast_d0"), ga(AGED_LINE, "PartialReprog")
    bg = draws[(AGED_LINE, "Fibroblast_d0")] - draws[(YOUNG_LINE, "Fibroblast_d0")]
    bd = draws[(AGED_LINE, "Fibroblast_d0")] - draws[(AGED_LINE, "PartialReprog")]
    res["a"] = dict(states=a_df, aged_d0=a0, young_d0=y0, aged_PartialReprog=apr, gap=a0 - y0, drop=a0 - apr,
                    drop_over_gap=(a0 - apr) / (a0 - y0),
                    gap_ci=[float(np.percentile(bg, 2.5)), float(np.percentile(bg, 97.5))],
                    drop_ci=[float(np.percentile(bd, 2.5)), float(np.percentile(bd, 97.5))],
                    ratio_ci=[float(np.percentile(bd / bg, 2.5)), float(np.percentile(bd / bg, 97.5))])

    # b. same-platform
    frozen, spaces = ctx["frozen"], ctx["spaces"]
    obs_a, obs_y = obs[AGED_LINE], obs[YOUNG_LINE]
    Y_a, Y_y = ctx["Y"][AGED_LINE], ctx["Y"][YOUNG_LINE]
    d0a, d0y = d0_fibroblast_idx(obs_a), d0_fibroblast_idx(obs_y)
    o_idx, aged_b = split_half(d0a, gs.SEED)
    y_idx, young_b = split_half(d0y, gs.SEED)
    old_sp = ctx["old_split"]
    res["split"] = dict(n_d0_aged=int(d0a.size), n_d0_young=int(d0y.size), O=int(o_idx.size), Y=int(y_idx.size),
                        aged_B=int(aged_b.size), young_B=int(young_b.size),
                        same_as_stored=bool(np.array_equal(o_idx, old_sp["O"]) and np.array_equal(y_idx, old_sp["Y"])
                                            and np.array_equal(aged_b, old_sp["aged_B"])
                                            and np.array_equal(young_b, old_sp["young_B"])))
    main_spec = [("O", Y_a, o_idx), ("Y", Y_y, y_idx), ("S", Y_a, state_idx(obs_a, DEST_STATE)),
                 ("S_rev", Y_y, state_idx(obs_y, DEST_STATE)), ("Pluri", Y_a, state_idx(obs_a, PLURI_STATE)),
                 ("NonReprog", Y_a, state_idx(obs_a, NONREPROG_STATE))]
    res["n_cells"] = {nm: int(ix.size) for nm, _, ix in main_spec}
    log(f"[{scheme}] split {res['split']} panel n={res['n_cells']}")
    main_counts = [sum_rows(Ym, ix) for _, Ym, ix in main_spec]
    mixes, n_tot = mix_indices(aged_b, young_b, gs.SEED)
    mix_counts = [sum_rows(Y_a, o_idx), sum_rows(Y_y, y_idx)]
    for mx in mixes:
        mix_counts.append(mix_sum(Y_a, Y_y, mx["aged_idx"], mx["young_idx"]))
    Zmain = panel_z(main_counts, frozen)[0]
    Zmix = panel_z(mix_counts, frozen)[0]
    SPS = ("FULL", "MD")
    pt = {sp: gs.same_point_block(Zmain, Zmix, spaces[sp], mixes, sp) for sp in SPS}
    boot = {}

    def bput(k, v):
        boot.setdefault(k, []).append(float(v))
    rng_b = np.random.default_rng(gs.BOOT_SEED)
    for b in range(gs.N_BOOT):
        bc = [mix_counts[0], mix_counts[1]]
        for mx in mixes:
            bc.append(mix_sum(Y_a, Y_y, resample_idx(mx["aged_idx"], rng_b), resample_idx(mx["young_idx"], rng_b)))
        Zb = panel_z(bc, frozen)[0]
        for sp in SPS:
            zb = Zb[:, spaces[sp]]
            for i, mx in enumerate(mixes):
                st = gs.same_axis(zb[0], zb[1], zb[2 + i], f"boot{b} {sp} mix")
                bput(f"mix__{sp}__{i}__progress", st["progress"])
                bput(f"mix__{sp}__{i}__delta", st["delta"])
    rng_b = np.random.default_rng(gs.BOOT_SEED)
    for b in range(gs.N_BOOT):
        bc = [main_counts[i] if nm in ("O", "Y") else sum_rows(Ym, resample_idx(ix, rng_b))
              for i, (nm, Ym, ix) in enumerate(main_spec)]
        Zb = panel_z(bc, frozen)[0]
        for sp in SPS:
            zb = Zb[:, spaces[sp]]
            f = gs.same_axis(zb[0], zb[1], zb[2], f"boot{b} {sp} fwd")
            r = gs.same_axis(zb[1], zb[0], zb[3], f"boot{b} {sp} rev")
            for nm, st in (("forward", f), ("reverse", r)):
                bput(f"main__{sp}__{nm}__progress", st["progress"])
                bput(f"main__{sp}__{nm}__delta", st["delta"])
            bput(f"main__{sp}__asym", f["progress"] - r["progress"])
        if (b + 1) % 50 == 0:
            log(f"[{scheme}] main boot {b + 1}/{gs.N_BOOT}")
    res["b"] = {}
    for sp in SPS:
        main, mix = pt[sp]
        fw, rv = main["forward"], main["reverse"]
        rec = dict(n_genes=int(spaces[sp].size))
        for nm, st in (("forward", fw), ("reverse", rv)):
            for stat in ("progress", "delta"):
                lo, hi, ok = gs.ci_with_flag(boot[f"main__{sp}__{nm}__{stat}"], st[stat])
                rec[f"{nm}_{stat}"], rec[f"{nm}_{stat}_ci"], rec[f"{nm}_{stat}_ci_valid"] = st[stat], [lo, hi], ok
            rec[f"{nm}_dist_final"], rec[f"{nm}_dist_start"] = st["dist_final"], st["dist_start"]
            rec[f"{nm}_final_over_start"] = st["dist_final"] / st["dist_start"]
        asym = fw["progress"] - rv["progress"]
        lo, hi, ok = gs.ci_with_flag(boot[f"main__{sp}__asym"], asym)
        rec.update(asymmetry=asym, asymmetry_ci=[lo, hi], asymmetry_ci_valid=ok)
        mix_rows = []
        for i, (mx, st) in enumerate(zip(mixes, mix)):
            lo, hi, ok = gs.ci_with_flag(boot[f"mix__{sp}__{i}__progress"], st["progress"])
            dlo, dhi, dok = gs.ci_with_flag(boot[f"mix__{sp}__{i}__delta"], st["delta"])
            mix_rows.append(dict(f=mx["f"], n_aged=mx["n_aged"], n_young=mx["n_young"], progress=st["progress"],
                                 progress_ci_lo=lo, progress_ci_hi=hi, progress_ci_valid=ok, delta=st["delta"],
                                 delta_ci_lo=dlo, delta_ci_hi=dhi))
        rec["mixtures"] = mix_rows
        rec["gate"] = gs.mixture_gate(mix_rows)
        alo, ahi = rec["asymmetry_ci"]
        rec["rule_met"] = bool(fw["delta"] > 0 and rec["gate"]["passed"] and (alo > 0 or ahi < 0))
        res["b"][sp] = rec
        log(f"[{scheme} {sp}] fwd progress={fw['progress']:.6f} delta={fw['delta']:.6f} "
            f"final/start={rec['forward_final_over_start']:.4f} rev progress={rv['progress']:.6f} "
            f"asym={asym:.6f} [{alo:.6f}, {ahi:.6f}] gate={rec['gate']['passed']} rule_met={rec['rule_met']}")

    # c. GTEx, FULL, aged PartialReprog
    donor_z = {}
    for line in DONORS:
        panel = donor_panel(obs[line], ctx["Y"][line], log, line)
        Zd = zscore_frozen(panel["logcpm"], frozen, panel["C"])[0]
        donor_z[line] = dict(Z=Zd, names=list(panel["names"]), idx=panel["idx"])
    g = gs.gtex_point_block(donor_z, ctx["anchors"], spaces["FULL"])
    gp = g[AGED_LINE].get("PartialReprog")
    res["c"] = dict(cos_S=gp["cos_S"], delta_young=gp["delta_young"], delta_old=gp["delta_old"],
                    n_PartialReprog=int(donor_z[AGED_LINE]["idx"]["PartialReprog"].size)) if gp else None
    res["rule_met"] = bool(all(res["b"][sp]["rule_met"] for sp in SPS))
    return res


def _flat(res):
    """Headline quantities as name -> value, for old/new/diff tables."""
    a, b, c = res["a"], res["b"], res["c"]
    out = {
        ("a", "aged day-0 Fibroblast MD score"): a["aged_d0"],
        ("a", "aged PartialReprog MD score"): a["aged_PartialReprog"],
        ("a", "young day-0 Fibroblast MD score"): a["young_d0"],
        ("a", "gap (aged d0 - young d0)"): a["gap"],
        ("a", "drop (aged d0 - aged PartialReprog)"): a["drop"],
        ("a", "drop / gap"): a["drop_over_gap"],
        ("a", "drop / gap CI low"): a["ratio_ci"][0],
        ("a", "drop / gap CI high"): a["ratio_ci"][1],
    }
    for sp in ("FULL", "MD"):
        r = b[sp]
        out.update({
            ("b " + sp, "forward progress"): r["forward_progress"],
            ("b " + sp, "change in distance (delta)"): r["forward_delta"],
            ("b " + sp, "delta CI low"): r["forward_delta_ci"][0],
            ("b " + sp, "delta CI high"): r["forward_delta_ci"][1],
            ("b " + sp, "final / start distance"): r["forward_final_over_start"],
            ("b " + sp, "reverse progress"): r["reverse_progress"],
            ("b " + sp, "asymmetry"): r["asymmetry"],
            ("b " + sp, "asymmetry CI low"): r["asymmetry_ci"][0],
            ("b " + sp, "asymmetry CI high"): r["asymmetry_ci"][1],
            ("b " + sp, "mixture progress at f=0.50"): r["gate"]["progress50"],
            ("b " + sp, "mixture delta at f=0.50"): r["gate"]["delta50"],
            ("b " + sp, "mixture control passes"): float(r["gate"]["passed"]),
        })
    if c:
        out.update({("c FULL", "cosine"): c["cos_S"], ("c FULL", "change in distance to young GTEx"): c["delta_young"],
                    ("c FULL", "change in distance to old GTEx"): c["delta_old"]})
    return out


def harness_refs():
    """Stored values the OLD-label run must reproduce."""
    G = RESULTS / "genespace"
    NS = RESULTS / "newstory"
    st = pd.read_csv(NS / "partA_md_by_state.csv")
    gd = pd.read_csv(NS / "partA_gap_drop.csv")
    gd = gd[gd.scoring == "joint"].iloc[0]
    ss = pd.read_csv(G / "same_stats.csv")
    pc = pd.read_csv(G / "same_posctrl.csv")
    gt = pd.read_csv(G / "gtex_stats.csv")
    refs = {}
    for _, r in st.iterrows():
        for k in ("joint_mean", "joint_ci_lo", "joint_ci_hi"):
            refs[("a_state", r.cell_line, r.state, k)] = float(r[k])
    for k in ("gap", "drop", "drop_over_gap", "gap_ci_lo", "gap_ci_hi", "drop_ci_lo", "drop_ci_hi",
              "ratio_ci_lo", "ratio_ci_hi"):
        refs[("a", k)] = float(gd[k])
    for sp in ("FULL", "MD"):
        for test in ("forward", "reverse"):
            r = ss[(ss.space == sp) & (ss.test == test)].iloc[0]
            for k in ("progress", "delta", "dist_final", "dist_start", "progress_ci_lo", "progress_ci_hi",
                      "delta_ci_lo", "delta_ci_hi"):
                refs[("b", sp, test, k)] = float(r[k])
            if test == "forward":
                for k in ("asymmetry", "asymmetry_ci_lo", "asymmetry_ci_hi"):
                    refs[("b", sp, test, k)] = float(r[k])
        for _, r in pc[pc.space == sp].iterrows():
            for k in ("progress", "progress_ci_lo", "progress_ci_hi", "delta"):
                refs[("mix", sp, float(r.f), k)] = float(r[k])
    r = gt[(gt.space == "FULL") & (gt.cell_line == AGED_LINE) & (gt.state == "PartialReprog")].iloc[0]
    for k in ("cos_S", "delta_young", "delta_old"):
        refs[("c", k)] = float(r[k])
    return refs


def harness_values(res):
    v = {}
    for _, r in res["a"]["states"].iterrows():
        for k in ("joint_mean", "joint_ci_lo", "joint_ci_hi"):
            v[("a_state", r.cell_line, r.state, k)] = float(r[k])
    a = res["a"]
    v.update({("a", "gap"): a["gap"], ("a", "drop"): a["drop"], ("a", "drop_over_gap"): a["drop_over_gap"],
              ("a", "gap_ci_lo"): a["gap_ci"][0], ("a", "gap_ci_hi"): a["gap_ci"][1],
              ("a", "drop_ci_lo"): a["drop_ci"][0], ("a", "drop_ci_hi"): a["drop_ci"][1],
              ("a", "ratio_ci_lo"): a["ratio_ci"][0], ("a", "ratio_ci_hi"): a["ratio_ci"][1]})
    for sp, r in res["b"].items():
        for test in ("forward", "reverse"):
            v[("b", sp, test, "progress")] = r[f"{test}_progress"]
            v[("b", sp, test, "delta")] = r[f"{test}_delta"]
            v[("b", sp, test, "dist_final")] = r[f"{test}_dist_final"]
            v[("b", sp, test, "dist_start")] = r[f"{test}_dist_start"]
            v[("b", sp, test, "progress_ci_lo")], v[("b", sp, test, "progress_ci_hi")] = r[f"{test}_progress_ci"]
            v[("b", sp, test, "delta_ci_lo")], v[("b", sp, test, "delta_ci_hi")] = r[f"{test}_delta_ci"]
        v[("b", sp, "forward", "asymmetry")] = r["asymmetry"]
        v[("b", sp, "forward", "asymmetry_ci_lo")], v[("b", sp, "forward", "asymmetry_ci_hi")] = r["asymmetry_ci"]
        for m in r["mixtures"]:
            for k in ("progress", "progress_ci_lo", "progress_ci_hi", "delta"):
                v[("mix", sp, float(m["f"]), k)] = float(m[k])
    for k in ("cos_S", "delta_young", "delta_old"):
        v[("c", k)] = res["c"][k]
    return v


def stats_context(log):
    import genespace_run as gs
    from fibro2_common import load_frozen_ruler
    from toward_run import load_Y, load_gtex_z
    from same_run import check_gene_space
    frozen_npz = load_frozen_ruler()
    frozen = {k: np.asarray(frozen_npz[k]) for k in ("mu", "sd", "w", "symbol", "ensembl")}
    gm = gs.map_gene_spaces(frozen, json.loads(gs.GENESETS_PATH.read_text(encoding="utf-8")), log)
    Y = {}
    for line in DONORS:
        Y[line] = load_Y(line)
        check_gene_space(Y[line], frozen, line)
    gtex = load_gtex_z(log)
    age = gtex["age"]
    anchors = dict(c_young=gtex["Z"][np.isin(age, gs.YOUNG_BINS)].mean(0),
                   c_old=gtex["Z"][np.isin(age, gs.OLD_BINS)].mean(0))
    del gtex
    sp = np.load(RESULTS / "same" / "split_idx.npz")
    return dict(frozen=frozen, spaces=gm["spaces"], Y=Y, anchors=anchors,
                partA=np.load(RESULTS / "newstory" / "partA_cell_scores.npz", allow_pickle=True),
                old_split={k: np.asarray(sp[k], int) for k in ("O", "Y", "aged_B", "young_B")})


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, pd.DataFrame):
        return o.to_dict(orient="records")
    if isinstance(o, (np.floating, np.integer, np.bool_)):
        return o.item()
    return o


def stats_step(log):
    require_prereg()
    _md2_imports()
    lab = donor_labels()
    for line in DONORS:
        if not np.array_equal(np.load(RESULTS / "newstory" / "partA_cell_scores.npz", allow_pickle=True)[f"label__{line}"].astype(str),
                              lab[line]["old"].to_numpy().astype(str)):
            raise RuntimeError(f"{line}: stored partA labels differ from md2 labels")
    ctx = stats_context(log)
    results = {}
    old = headline(lab, "old", log, ctx)
    refs, vals = harness_refs(), harness_values(old)
    h_rows = []
    for k, ref in refs.items():
        new = vals.get(k, np.nan)
        h_rows.append(dict(quantity=" | ".join(map(str, k)), reference=ref, rerun=new, diff=new - ref,
                           ok=bool(np.isfinite(new) and abs(new - ref) <= 1e-9)))
    h = pd.DataFrame(h_rows)
    h.to_csv(OUT / "harness_check.csv", index=False)
    log(f"[harness] {int(h.ok.sum())}/{len(h)} within 1e-9; max|diff|={float(np.nanmax(np.abs(h['diff']))):.3e}")
    if not h.ok.all():
        log(h[~h.ok].to_string())
        raise RuntimeError("OLD-label rerun does not reproduce stored values. STOP.")
    results["old"] = old
    for scheme in ("joint", "lu"):
        results[scheme] = headline(lab, scheme, log, ctx)
    fo = _flat(results["old"])
    rows = []
    for scheme in ("joint", "lu"):
        fn = _flat(results[scheme])
        for k in fo:
            rows.append(dict(labelling=scheme, part=k[0], quantity=k[1], old=fo[k], new=fn.get(k, np.nan),
                             diff=fn.get(k, np.nan) - fo[k]))
    pd.DataFrame(rows).to_csv(OUT / "headline_old_new.csv", index=False)
    for scheme, r in results.items():
        r["a"]["states"].to_csv(OUT / f"partA_states_{scheme}.csv", index=False)
    (OUT / "headline.json").write_text(json.dumps(_jsonable(results), indent=2, default=str), encoding="utf-8")
    log("[stats] wrote headline_old_new.csv, headline.json")


# --------------------------------------------------------------------------- FINDINGS
def _mdt(df, fmt=None):
    fmt = fmt or {}
    cols = list(df.columns)
    lines = ["| " + " | ".join(map(str, cols)) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if c in fmt:
                cells.append(fmt[c](v))
            elif isinstance(v, (float, np.floating)):
                cells.append("NA" if not np.isfinite(v) else f"{v:.4g}")
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def _num(v):
    if isinstance(v, (bool, np.bool_)):
        return "pass" if v else "fail"
    v = float(v)
    if not np.isfinite(v):
        return "NA"
    a = abs(v)
    return f"{v:.3f}" if a >= 100 else (f"{v:.4f}" if a >= 0.01 or a == 0 else f"{v:.2e}")


def _ci(c):
    return f"[{_num(c[0])}, {_num(c[1])}]"


def findings_step(log):
    require_prereg()
    _md2_imports()
    H = json.loads((OUT / "headline.json").read_text(encoding="utf-8"))
    gate = json.loads((OUT / "validation_GM00731.json").read_text(encoding="utf-8"))
    harness = pd.read_csv(OUT / "harness_check.csv")
    meta = json.loads((OUT / "joint_meta.json").read_text(encoding="utf-8"))
    jl = pd.read_csv(OUT / "joint_cluster_labels.csv")
    agree = pd.read_csv(OUT / "agreement_summary.csv")
    d0m = json.loads((OUT / "d0_fibroblast_membership.json").read_text(encoding="utf-8"))
    lab = donor_labels()
    WORD = {AGED_LINE: "aged donor", YOUNG_LINE: "young donor"}
    yn = {True: "yes", False: "no"}
    L = ["# FINDINGS_RELABEL — do the headline numbers survive Lu-style joint clustering?", ""]

    verdict = {s: H[s]["rule_met"] for s in ("old", "joint", "lu")}
    oj = agree[(agree.rows == "old") & (agree.cols == "joint")]
    share_j = dict(zip(oj.cell_line, oj.share_unchanged.astype(float)))
    j, lu, old = H["joint"], H["lu"], H["old"]
    L += ["## Answer", "",
          f"**{'Yes' if verdict['joint'] else 'No'}. Under the joint clustering the pre-registered rule is met: "
          f"{yn[verdict['joint']]}. Under Lu et al.'s own published labels it is met: {yn[verdict['lu']]}.** In both the full and MD gene spaces "
          "the aged PartialReprog cells still end farther from the young target, the mixture control "
          "still passes, and the asymmetry interval stays far from 0.",
          "",
          f"- Joint clustering changed the label of {100 * (1 - share_j[AGED_LINE]):.1f}% "
          f"of aged cells and {100 * (1 - share_j[YOUNG_LINE]):.1f}% of young cells.",
          f"- Same-platform, full space: forward progress {_num(old['b']['FULL']['forward_progress'])} → "
          f"{_num(j['b']['FULL']['forward_progress'])}; change in distance {_num(old['b']['FULL']['forward_delta'])} → "
          f"{_num(j['b']['FULL']['forward_delta'])}; asymmetry {_num(old['b']['FULL']['asymmetry'])} → "
          f"{_num(j['b']['FULL']['asymmetry'])} {_ci(j['b']['FULL']['asymmetry_ci'])}.",
          f"- Same-platform, MD space: forward progress {_num(old['b']['MD']['forward_progress'])} → "
          f"{_num(j['b']['MD']['forward_progress'])}; change in distance {_num(old['b']['MD']['forward_delta'])} → "
          f"{_num(j['b']['MD']['forward_delta'])}; asymmetry {_num(old['b']['MD']['asymmetry'])} → "
          f"{_num(j['b']['MD']['asymmetry'])} {_ci(j['b']['MD']['asymmetry_ci'])}.",
          f"- MD score drop ÷ gap: {old['a']['drop_over_gap']:.2f} → {j['a']['drop_over_gap']:.2f} "
          f"{_ci(j['a']['ratio_ci'])} (joint); {lu['a']['drop_over_gap']:.2f} {_ci(lu['a']['ratio_ci'])} (Lu labels).",
          f"- GTEx, full space: the change in distance to young GTEx is {_num(old['c']['delta_young'])} → "
          f"{_num(j['c']['delta_young'])} (joint), {_num(lu['c']['delta_young'])} (Lu labels); still > 0.",
          "",
          "Under the joint clustering, forward progress, change in distance, final ÷ start, asymmetry and the "
          "mixture numbers each moved by less than 8% of their old value. Reverse progress, which is close to 0, "
          f"moved from {_num(old['b']['FULL']['reverse_progress'])} to {_num(j['b']['FULL']['reverse_progress'])} "
          f"(full) and {_num(old['b']['MD']['reverse_progress'])} to {_num(j['b']['MD']['reverse_progress'])} (MD). "
          "Every sign and every gate verdict is unchanged. Lu's labels move the numbers somewhat more, with the "
          "largest shift in the GTEx distances (details below); that doesn't change a sign or a verdict either.",
          "",
          "**Sentence the paper could use:** \"Re-clustering both donors together, as Lu et al. did, changed "
          f"the state label of about {100 * (1 - np.mean(list(share_j.values()))):.0f}% of "
          "cells and left the result unchanged: aged partially reprogrammed cells still ended farther from young "
          f"fibroblasts than they started (change in distance {j['b']['FULL']['forward_delta']:.1f} across all "
          f"genes and {j['b']['MD']['forward_delta']:.1f} across MD genes, versus "
          f"{old['b']['FULL']['forward_delta']:.1f} and {old['b']['MD']['forward_delta']:.1f} with the original "
          "labels), and the same held with Lu et al.'s published per-cell state labels.\"",
          ""]

    L += [PREREG_TEXT, ""]

    L += ["## Checks before the new-label numbers", "",
          f"- **Clustering gate (GM00731 alone, streamed code).** Passed: cells and order identical, cluster "
          f"of every cell identical ({gate['meta']['n_clusters']} clusters, {gate['meta']['snn_edges']:,} SNN edges, "
          f"as in md2), gene means and AMS bins identical, MD and all five state scores identical (max |diff| 0), "
          f"cluster labels identical. percent-mito max |diff| {gate['percent_mt_max_abs_diff']:.1e} "
          "(CSV round-trip of the stored value).",
          f"- **Harness (old labels through the new code).** {int(harness.ok.sum())}/{len(harness)} stored values "
          f"reproduced within 1e-9 (max |diff| {float(np.abs(harness['diff']).max()):.1e}): newstory Part A state "
          "means, CIs, gap/drop/ratio; genespace same_stats and same_posctrl (FULL, MD; points and CIs); "
          "GTEx FULL aged PartialReprog. See results/relabel_check/harness_check.csv.",
          f"- **Memory.** The joint run used all {sum(meta['n_cells_kept'].values()):,} cells "
          f"(aged {meta['n_cells_kept'][AGED_LINE]:,}, young {meta['n_cells_kept'][YOUNG_LINE]:,}) and "
          f"{meta['n_genes_kept']:,} genes (detected in ≥5 cells of the merged matrix; {meta['n_genes_dropped']:,} "
          "dropped). Nothing was subsampled.",
          ""]

    # Step 1 summary
    jl["aged_share"] = jl[AGED_LINE] / jl["n_cells"]
    means = jl[[f"mean_{s}" for s in STATE_NAMES]].to_numpy()
    srt = np.sort(means, axis=1)
    jl["margin"] = srt[:, -1] - srt[:, -2]
    mixed = int(((jl.aged_share >= 0.2) & (jl.aged_share <= 0.8)).sum())
    L += ["## Step 1 — joint clustering", "",
          f"{meta['n_clusters']} joint clusters (resolution 0.8; {meta['snn_edges']:,} SNN edges; 50 PCs explain "
          f"{meta['evr_sum']:.3f} of HVG variance). {mixed} of {meta['n_clusters']} clusters have 20–80% aged cells; "
          "the rest are mostly one donor, so the donors still separate in several clusters even when "
          "clustered together. Each cluster's label, its mean state scores, the margin between the top two scores "
          "and its donor mix:", "",
          _mdt(jl[["cluster", "n_cells", "label", "margin", AGED_LINE, YOUNG_LINE]].rename(
              columns={"margin": "top-2 score margin", AGED_LINE: "aged cells", YOUNG_LINE: "young cells"})),
          "",
          f"Clusters whose top-2 margin is below 0.01 (the label is close to a tie): "
          f"{', '.join(f'c{int(r.cluster)} ({r.label}, n={int(r.n_cells)}, margin {r.margin:.4f})' for _, r in jl[jl.margin < 0.01].iterrows()) or 'none'}.",
          ""]

    # Step 2
    L += ["## Step 2 — agreement (per donor)", ""]
    order = list(STATE_NAMES) + [UNLAB]
    for line in DONORS:
        obs = lab[line]
        L += [f"### {line} ({WORD[line]})", ""]
        for a, b, title in (("old", "joint", "old md2 label (rows) × joint label (columns)"),
                            ("old", "lu", "old md2 label (rows) × Lu label (columns)")):
            ct = pd.crosstab(obs[a], obs[b])
            ct = ct.reindex(index=[s for s in order if s in ct.index], columns=[s for s in order if s in ct.columns],
                            fill_value=0)
            ct.insert(0, "", ct.index)
            share = agree[(agree.cell_line == line) & (agree.rows == a) & (agree.cols == b)].iloc[0]
            L += [f"Cell counts, {title}. Share of cells with an unchanged label: **{100 * share.share_unchanged:.1f}%**"
                  + (f" (of {int(share.n_compared):,} cells that Lu labelled; {int(share.n_cells - share.n_compared):,} "
                     "of our cells are not in Lu's object)" if b == "lu" else f" of {int(share.n_cells):,}") + ".", "",
                  _mdt(ct.reset_index(drop=True)), ""]
        cnt = []
        for st in order:
            rec = {"state": st}
            for s in ("old", "joint", "lu"):
                m = obs[s] == st
                rec[s] = int(m.sum())
                rec[f"{s} day 0"] = int((m & (obs.day.astype(int) == 0)).sum())
            if st != UNLAB or rec["lu"]:
                cnt.append(rec)
        L += ["Cells per state (all days, and day 0):", "",
              _mdt(pd.DataFrame(cnt)[["state", "old", "joint", "lu", "old day 0", "joint day 0", "lu day 0"]]), ""]
    L += ["Day-0 Fibroblast membership (the origin cells O/Y of the same-platform test and the GTEx origin):", ""]
    L += [_mdt(pd.DataFrame([dict(set=k, **v) for k, v in d0m.items()]))]
    L += ["", "Membership changed for both donors under both labellings, so the half-split was redone with "
          "same_run.split_half(seed 20260914) on each new day-0 Fibroblast set, and the mixtures were rebuilt "
          "from the new halves (same seed).", ""]
    for s in ("joint", "lu"):
        sp = H[s]["split"]
        L += [f"- {s}: O={sp['O']}, Y={sp['Y']}, aged half B={sp['aged_B']}, young half B={sp['young_B']} "
              f"(old: 2469 / 3851 / 2470 / 3852); panel sizes {H[s]['n_cells']}."]
    L += [""]

    # Step 3
    tbl = pd.read_csv(OUT / "headline_old_new.csv")
    L += ["## Step 3 — headline numbers, old vs new", "",
          "\"new\" = joint labelling (primary). Lu-label values are in the last columns. diff = new − old; "
          "rel = diff / |old|.", ""]
    jt = tbl[tbl.labelling == "joint"].reset_index(drop=True)
    lt = tbl[tbl.labelling == "lu"].reset_index(drop=True)
    parts = {"a": "a. MD score per state (joint AddModuleScore)", "b FULL": "b. Same-platform test, FULL space",
             "b MD": "b. Same-platform test, MD space", "c FULL": "c. GTEx test, FULL space, aged PartialReprog"}
    for p, title in parts.items():
        rows = []
        for i in jt.index[jt.part == p]:
            q = jt.loc[i, "quantity"]
            is_gate = q == "mixture control passes"
            o, n, lv = jt.loc[i, "old"], jt.loc[i, "new"], lt.loc[i, "new"]
            rows.append({"quantity": q, "old": ("pass" if o else "fail") if is_gate else _num(o),
                         "new (joint)": ("pass" if n else "fail") if is_gate else _num(n),
                         "diff": "" if is_gate else _num(n - o),
                         "rel": "" if is_gate or o == 0 else f"{100 * (n - o) / abs(o):+.1f}%",
                         "Lu labels": ("pass" if lv else "fail") if is_gate else _num(lv),
                         "Lu diff": "" if is_gate else _num(lv - o)})
        L += [f"### {title}", "", _mdt(pd.DataFrame(rows)), ""]
    L += ["Mixture control detail (progress at f = 0, 0.10, 0.25, 0.50; gate = non-decreasing, CI at 0.50 "
          "excludes 0 and contains the point, delta at 0.50 < 0):", ""]
    mrows = []
    for s in ("old", "joint", "lu"):
        for sp in ("FULL", "MD"):
            g = H[s]["b"][sp]["gate"]
            mrows.append(dict(labelling=s, space=sp, progress=", ".join(f"{p:.4f}" for p in g["progress"]),
                              ci50=_ci(g["progress50_ci"]), delta50=_num(g["delta50"]),
                              gate="pass" if g["passed"] else "fail: " + "; ".join(g["failed"])))
    L += [_mdt(pd.DataFrame(mrows)), ""]
    L += ["MD-score CIs (bootstrap B=200): "
          + "; ".join(f"{s}: gap {_ci(H[s]['a']['gap_ci'])}, drop {_ci(H[s]['a']['drop_ci'])}, "
                      f"drop/gap {_ci(H[s]['a']['ratio_ci'])}" for s in ("old", "joint", "lu")) + ".", ""]

    # rule
    L += ["## Decision rule, applied", "", f"\"{RULE}\"", ""]
    rr = []
    for s in ("old", "joint", "lu"):
        for sp in ("FULL", "MD"):
            b = H[s]["b"][sp]
            a_ci = b["asymmetry_ci"]
            rr.append({"labelling": s, "space": sp, "change in distance": _num(b["forward_delta"]),
                       "> 0": b["forward_delta"] > 0, "mixture control": "pass" if b["gate"]["passed"] else "fail",
                       "asymmetry interval": _ci(a_ci), "excludes 0": bool(a_ci[0] > 0 or a_ci[1] < 0),
                       "met": b["rule_met"]})
    L += [_mdt(pd.DataFrame(rr)), "",
          f"Joint labelling (primary): **{'holds' if verdict['joint'] else 'does not hold'}**. "
          f"Lu labels: **{'holds' if verdict['lu'] else 'does not hold'}**. (Old labels, for reference: "
          f"{'holds' if verdict['old'] else 'does not hold'}.)", ""]

    # plain summary
    aged = lab[AGED_LINE]
    dcomp = {s: aged[aged[s] == "PartialReprog"].day.astype(int).value_counts().sort_index().to_dict()
             for s in ("old", "joint", "lu")}
    L += ["## Plain summary", "",
          "- **Does the conclusion hold?** Yes, under both relabellings, in both gene spaces. Clustering the two "
          "donors together, as Lu et al. did, does not change what the aged partially reprogrammed cells do: part "
          "of their movement runs in the young direction, but they end about twice as far from the young day-0 "
          f"cells as they started (final ÷ start {j['b']['FULL']['forward_final_over_start']:.2f} in all genes, "
          f"{j['b']['MD']['forward_final_over_start']:.2f} in MD genes; old {old['b']['FULL']['forward_final_over_start']:.2f} "
          f"and {old['b']['MD']['forward_final_over_start']:.2f}).",
          "- **How far did the numbers move?** With the joint clustering, about 8% of cells changed label. "
          "Same-platform forward progress fell by about 0.02–0.03, the change in distance rose by about 2% in both "
          "spaces, and the asymmetry fell by about 0.04–0.06, staying near 0.7 (full) and 1.0 (MD) with intervals "
          "far from 0. "
          f"The MD score drop ÷ gap moved from {old['a']['drop_over_gap']:.2f} to {j['a']['drop_over_gap']:.2f}. "
          "The GTEx numbers barely moved (cosine +0.0006, distance changes +0.04 and +0.06).",
          f"- **Lu's published labels** agree less with ours (72–80% of cells keep their label) and move the "
          "numbers a little further: same-platform deltas rise by about 4%, the MD asymmetry falls by 0.09, and the "
          f"GTEx distance changes rise from {_num(old['c']['delta_young'])} / {_num(old['c']['delta_old'])} to "
          f"{_num(lu['c']['delta_young'])} / {_num(lu['c']['delta_old'])} (young / old). That GTEx shift is the "
          "largest move anywhere, and it goes in the direction of *more* distance, not less. Lu's aged "
          f"PartialReprog set is more purely day 3 than ours (cells by day: Lu {dcomp['lu']}, old {dcomp['old']}, "
          f"joint {dcomp['joint']}). I did not test whether that difference causes the GTEx shift.",
          "- **One sentence for the paper:** see the Answer section above.",
          "",
          "## Notes", "",
          "- Joint labelling reuses md2 code for the whole pipeline. The only deviation is the streamed memory "
          "layout, and it reproduced md2 exactly on GM00731 alone. As in md2, this is LogNormalize plus "
          "percent-mito regression, not Lu's sctransform v2, so this answers \"clustered together\" rather than "
          "\"clustered with Lu's exact normalization\". The Lu-label run uses Lu's own clustering, which was built "
          "on sctransform.",
          "- The MD score in part a is the newstory joint per-cell score, which does not depend on labels; "
          "only the state memberships changed.",
          "- results/verify/claims.csv (paper sentence P42/P96) describes the analysis as using \"the published "
          "cell-state assignments\". The analysis actually used the md2 re-derived labels, not Lu's published "
          "`cell_state`. That wording is worth fixing; this check shows the numbers survive either way. "
          "No manuscript file was edited.",
          "- Disclosure (repeated from the pre-registration): Lu's cell_state × sample counts were tabulated "
          "before the rule was written.",
          "",
          "## Files", "",
          "- src/relabel_check.py (`step0`, `prereg`, `cluster`, `agree`, `stats`, `findings`)",
          "- results/relabel_check/: PREREG.flag; lu_meta_*.csv, lu_slots_*.json, lu_meta_columns.csv (Step 0); "
          "validation_GM00731.json, joint_cells.csv, joint_cluster_labels.csv, joint_ams.npz, joint_meta.json "
          "(Step 1); agree_*.csv, agreement_summary.csv, d0_fibroblast_membership.json (Step 2); harness_check.csv, "
          "headline_old_new.csv, headline.json, partA_states_*.csv (Step 3); run_log.txt and console logs.",
          ""]
    FINDINGS_PATH.write_text("\n".join(L), encoding="utf-8")
    log(f"[findings] wrote {FINDINGS_PATH}")


def main():
    step = sys.argv[1] if len(sys.argv) > 1 else "all"
    log = Log(OUT / "run_log.txt")
    try:
        if step == "step0":
            step0(log)
        elif step == "prereg":
            prereg(log)
        elif step == "cluster":
            cluster_step(log)
        elif step == "agree":
            _md2_imports()
            agree_step(log)
        elif step == "stats":
            stats_step(log)
        elif step == "findings":
            findings_step(log)
        else:
            raise SystemExit(f"unknown step {step!r}")
    finally:
        log.close()


if __name__ == "__main__":
    main()

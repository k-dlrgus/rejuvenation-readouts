"""Single-pass streaming conversion of the gzipped DENSE TF Atlas h5ad (172 GB decompressed; does not fit on disk)
into a sparse CSR h5ad, reading the HDF5 file directly from the gzip stream (contiguous layout => sequential row reads).

Output: data/processed/GSE217460_TFAtlas_raw_sparse.h5ad  (anndata-compatible; X csr float32, obs/var DataFrames).
Then prints shape and obs columns (Phase 0 deliverable). Nothing is analysed.
"""
import gzip
import os
import sys
import time
import numpy as np
import scipy.sparse as sp
import h5py
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from config import TF_ATLAS, DATA_PROC, RESULTS  # noqa
from peek_gz_h5ad import SeekableGzip  # noqa

GZ = TF_ATLAS["local_dir"] / "GSE217460_210322_TFAtlas_raw.h5ad.gz"
OUT = DATA_PROC / "GSE217460_TFAtlas_raw_sparse.h5ad"
BLOCK = 2000  # rows per read (2000 x 37528 x 4 B = 300 MB)


def log(*a):
    print(*a, flush=True)


def read_legacy_df(g: h5py.Group) -> pd.DataFrame:
    """Read an anndata obs/var group written in the legacy (<0.8) or current encoding."""
    try:
        from anndata.io import read_elem
        return read_elem(g)
    except Exception as e:  # legacy format with __categories
        log(f"   read_elem failed ({type(e).__name__}); falling back to legacy reader")
    idx_key = g.attrs.get("_index", "_index")
    idx_key = idx_key.decode() if isinstance(idx_key, bytes) else idx_key
    cats = g["__categories"] if "__categories" in g else {}
    cols = {}
    for k in g.keys():
        if k in ("__categories", idx_key):
            continue
        v = g[k][()]
        if k in cats:
            c = cats[k][()]
            c = np.array([x.decode() if isinstance(x, bytes) else x for x in c], dtype=object)
            v = pd.Categorical.from_codes(v.astype(int), c)
        elif v.dtype.kind in ("S", "O"):
            v = np.array([x.decode() if isinstance(x, bytes) else x for x in v], dtype=object)
        cols[k] = v
    idx = g[idx_key][()]
    idx = np.array([x.decode() if isinstance(x, bytes) else x for x in idx], dtype=object)
    return pd.DataFrame(cols, index=idx)


def main():
    if OUT.exists():
        log(f"[cached] {OUT}")
    else:
        t0 = time.time()
        tmp = str(OUT) + ".part"
        sg = SeekableGzip(str(GZ))
        h = h5py.File(sg, "r")
        X = h["X"]
        n, m = X.shape
        off = X.id.get_offset()
        x_end = off + n * m * X.dtype.itemsize
        sg.tail_from = x_end  # everything after the dense block (obs, var, ...) is kept in RAM -> no rewinds for metadata
        log(f"[layout] X shape={X.shape} dtype={X.dtype} chunks={X.chunks} contiguous offset={off} -> {n*m*X.dtype.itemsize/1e9:.1f} GB dense; "
            f"tail cached from byte {x_end:,}")
        assert X.chunks is None and off is not None, "expected contiguous layout for sequential streaming"

        # ---- step 1: stream X -> CSR (resumable: skipped if a complete X exists in the .part file) ----
        x_done = False
        if os.path.exists(tmp):
            try:
                with h5py.File(tmp, "r") as o:
                    x_done = "X/indptr" in o and o["X/indptr"].shape[0] == n + 1
            except OSError:
                x_done = False
            if not x_done:
                os.remove(tmp)
        if x_done:
            log("[X] already converted in .part file; skipping (one decompression pass is still needed to reach obs/var)")
        else:
            with h5py.File(tmp, "w") as o:
                gX = o.create_group("X")
                gX.attrs["encoding-type"] = "csr_matrix"; gX.attrs["encoding-version"] = "0.1.0"
                gX.attrs["shape"] = np.array([n, m], dtype=np.int64)
                d_data = gX.create_dataset("data", shape=(0,), maxshape=(None,), dtype=np.float32, chunks=(1 << 20,), compression="lzf")
                d_idx = gX.create_dataset("indices", shape=(0,), maxshape=(None,), dtype=np.int32, chunks=(1 << 20,), compression="lzf")
                indptr = [0]
                nnz = 0
                for r0 in range(0, n, BLOCK):
                    r1 = min(r0 + BLOCK, n)
                    blk = X[r0:r1]
                    csr = sp.csr_matrix(blk)
                    k = csr.nnz
                    d_data.resize((nnz + k,)); d_idx.resize((nnz + k,))
                    d_data[nnz:nnz + k] = csr.data.astype(np.float32); d_idx[nnz:nnz + k] = csr.indices.astype(np.int32)
                    # nnz exceeds 2^31 for this matrix: keep all offset arithmetic in int64 (C long is 32-bit on Windows)
                    indptr.extend((np.int64(nnz) + csr.indptr[1:].astype(np.int64)).tolist())
                    nnz = int(nnz + k)
                    if (r0 // BLOCK) % 25 == 0:
                        el = time.time() - t0
                        log(f"   rows {r1:,}/{n:,} ({100*r1/n:.1f}%)  nnz={nnz:,}  density={nnz/(r1*m):.4f}  streamed {sg.bytes_read/1e9:.1f} GB "
                            f"@ {sg.bytes_read/el/1e6:.0f} MB/s  rewinds={sg.rewinds}  elapsed {el:.0f}s")
                gX.create_dataset("indptr", data=np.array(indptr, dtype=np.int64))
                log(f"[X] done: nnz={nnz:,} density={nnz/(n*m):.4f}  ({time.time()-t0:.0f}s, rewinds={sg.rewinds})")

        # ---- step 2: obs / var (located after X). Stream the rest of the file once into the tail cache so that every
        # HDF5 metadata read (object headers, B-trees, heaps, column datasets) is a RAM hit and can never trigger a rewind.
        eof = sg.prefetch_to_eof()
        log(f"[tail] streamed to EOF at byte {eof:,}; tail cached {len(sg.tail)/1e6:.0f} MB  (rewinds={sg.rewinds}, {time.time()-t0:.0f}s)")
        obs = read_legacy_df(h["obs"]); var = read_legacy_df(h["var"])
        log(f"[obs/var] obs={obs.shape} var={var.shape}  tail cached {len(sg.tail)/1e6:.0f} MB  (rewinds={sg.rewinds}, {time.time()-t0:.0f}s)")
        h.close(); sg.close()
        from anndata.io import write_elem
        with h5py.File(tmp, "a") as o:
            for k in ("obs", "var", "uns", "obsm", "varm", "obsp", "varp", "layers"):
                if k in o:
                    del o[k]
            write_elem(o, "obs", obs); write_elem(o, "var", var)
            o.create_group("uns").attrs.update({"encoding-type": "dict", "encoding-version": "0.1.0"})
            o["uns"].create_dataset("source", data=str(GZ.name))
            for k in ("obsm", "varm", "obsp", "varp", "layers"):
                o.create_group(k).attrs.update({"encoding-type": "dict", "encoding-version": "0.1.0"})
        os.replace(tmp, OUT)
        log(f"[sparse] wrote {OUT} ({os.path.getsize(OUT)/1e9:.2f} GB) in {time.time()-t0:.0f}s")

    import anndata as ad
    a = ad.read_h5ad(OUT, backed="r")
    lines = ["=" * 90,
             f"GSE216481 / GSE217460 TF Atlas raw counts (sparse-converted): shape (cells, genes) = {a.shape}",
             f"obs columns: {list(a.obs.columns)}",
             f"var columns: {list(a.var.columns)}  (var index = gene symbols, e.g. {list(a.var.index[:5])})"]
    if "TF" in a.obs:
        lines.append(f"n distinct TF labels: {a.obs['TF'].nunique():,}")
    if "batch" in a.obs:
        lines.append(f"n batches: {a.obs['batch'].nunique()}")
    pd.set_option("display.width", 200)
    lines.append(a.obs.head(3).T.to_string())
    with h5py.File(OUT, "r") as f:
        d = f["X/data"][:1_000_000]
        lines.append(f"X sample: dtype={d.dtype} min={d.min()} max={d.max()} integer-fraction={np.mean(np.abs(d-np.round(d))<1e-6):.4f}")
    lines.append("=" * 90)
    txt = "\n".join(lines)
    log(txt)
    (RESULTS / "phase0_GSE217460_full_verify.txt").write_text(txt, encoding="utf-8")


if __name__ == "__main__":
    main()

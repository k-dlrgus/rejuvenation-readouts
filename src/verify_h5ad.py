"""Verify an h5ad (optionally .gz) loads; print shape, obs columns, var columns, X dtype/format.

Usage: python src/verify_h5ad.py path/to/file.h5ad[.gz] [--obs-head N] [--decompress-to DIR]
Uses backed='r' mode so it never loads X into RAM. If the file is gzipped it is
decompressed once next to the source (cached; skipped if the target already exists).
"""
import argparse
import gzip
import os
import shutil
import sys
import h5py
import anndata as ad
import pandas as pd


def gunzip_once(src: str) -> str:
    dst = src[:-3]
    if os.path.exists(dst) and os.path.getsize(dst) > 0:
        print(f"[cached] {dst} already decompressed", flush=True)
        return dst
    print(f"[gunzip] {src} -> {dst}", flush=True)
    tmp = dst + ".part"
    with gzip.open(src, "rb") as fi, open(tmp, "wb") as fo:
        shutil.copyfileobj(fi, fo, 1 << 24)
    os.replace(tmp, dst)
    return dst


def describe(path: str, obs_head: int = 5):
    print("=" * 90)
    print("FILE:", path, f"({os.path.getsize(path)/1e9:.2f} GB on disk)")
    with h5py.File(path, "r") as f:
        print("top-level keys:", list(f.keys()))
        X = f.get("X")
        if isinstance(X, h5py.Group):
            print("X: sparse group, encoding:", dict(X.attrs), "data dtype:", X["data"].dtype, "nnz:", X["data"].shape[0])
        elif X is not None:
            print("X: dense dataset", X.shape, X.dtype)
        if "raw" in f:
            rx = f["raw"].get("X")
            if isinstance(rx, h5py.Group):
                print("raw/X: sparse", dict(rx.attrs), "data dtype:", rx["data"].dtype)
            elif rx is not None:
                print("raw/X: dense", rx.shape, rx.dtype)
        if "layers" in f:
            print("layers:", list(f["layers"].keys()))
    a = ad.read_h5ad(path, backed="r")
    print("shape (cells, genes):", a.shape)
    print("obs columns:", list(a.obs.columns))
    print("var columns:", list(a.var.columns))
    print("uns keys:", list(a.uns.keys())[:30])
    print("obsm keys:", list(a.obsm.keys()))
    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 40)
    print(f"obs head({obs_head}):")
    print(a.obs.head(obs_head).T.to_string())
    print("=" * 90)
    return a


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--obs-head", type=int, default=3)
    a = ap.parse_args()
    p = a.path
    if p.endswith(".gz"):
        p = gunzip_once(p)
    describe(p, a.obs_head)

"""Inflate one Cell Ranger matrix member and write Gene Expression feature names."""
import zlib
from pathlib import Path

import h5py

src = Path("results/survey/opened/southard/L4A_matrix.h5.deflate")
dst = Path("results/survey/opened/southard/L4A_matrix.h5")
raw = zlib.decompress(src.read_bytes(), -15)
dst.write_bytes(raw)
print("inflated", len(raw))

with h5py.File(dst, "r") as f:
    def walk(g, p=""):
        for k in g.keys():
            item = g[k]
            path = f"{p}/{k}"
            if isinstance(item, h5py.Dataset):
                print(path, item.shape, item.dtype)
            else:
                walk(item, path)
    walk(f)
    names = [x.decode() if isinstance(x, bytes) else str(x) for x in f["matrix/features/name"][:]]
    ids = [x.decode() if isinstance(x, bytes) else str(x) for x in f["matrix/features/id"][:]]
    types = [x.decode() if isinstance(x, bytes) else str(x) for x in f["matrix/features/feature_type"][:]]
from collections import Counter
print(Counter(types))
genes = [n for n, t in zip(names, types) if t == "Gene Expression"]
print("gene_expression", len(genes), "unique", len(set(genes)))
out = Path("results/survey/gene_lists/southard_hs27_L4A_genes.txt")
out.write_text("\n".join(genes) + "\n", encoding="utf-8")
print("wrote", out)

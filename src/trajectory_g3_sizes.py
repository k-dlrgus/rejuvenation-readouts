"""Peek remote sizes for G3 downloads."""
from download import remote_size
urls = [
    "https://datasets.cellxgene.cziscience.com/11554919-8ee5-411b-a088-7be1a7c9a5a6.h5ad",
    "https://datasets.cellxgene.cziscience.com/b14cb09f-dd56-4c6c-bbf1-6822f01af53f.h5ad",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE224nnn/GSE224438/suppl/GSE224438_RAW.tar",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE276nnn/GSE276656/suppl/GSE276656_RAW.tar",
]
for u in urls:
    n = remote_size(u)
    print(n, None if n is None else round(n / 1e6, 1), "MB", u.rsplit("/", 1)[-1], flush=True)

"""Download G3 files whose URLs were fetched in g3_search (not invented)."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import DATA_RAW
from download import download

RAW = DATA_RAW / "trajectory"
JOBS = [
    ("tms_brain_nonmyeloid",
     "https://datasets.cellxgene.cziscience.com/11554919-8ee5-411b-a088-7be1a7c9a5a6.h5ad"),
    ("tms_brain_myeloid",
     "https://datasets.cellxgene.cziscience.com/b14cb09f-dd56-4c6c-bbf1-6822f01af53f.h5ad"),
    ("GSE224438",
     "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE224nnn/GSE224438/suppl/GSE224438_RAW.tar"),
]
if __name__ == "__main__":
    for name, url in JOBS:
        dest = RAW / name
        dest.mkdir(parents=True, exist_ok=True)
        print(f"=== {name} ===", flush=True)
        download(url, str(dest))

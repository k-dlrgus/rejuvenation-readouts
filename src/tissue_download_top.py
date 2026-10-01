"""Download named CELLxGENE h5ads into data/raw/tissue/<name>/. Accessions and URLs are read from
the fetched table results/tissue/t1_cxg_primary_candidates.csv — never invented.
Usage: python src/tissue_download_top.py retina_sn brain_aging muscle
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tissue_common import TISSUE_DIR, TISSUE_RAW  # noqa: E402
from download import download  # noqa: E402

# Short names -> CXG dataset_id (IDs taken from the fetched primary-candidate table).
REGISTRY = {
    "retina_sn": "d6505c89-c43d-4c28-8c4f-7351a5fd5528",       # 104 donors, 37 GB
    "retina_sc": "be41a86a-b606-4b1c-8055-32f334898775",       # 20 donors, 2.5 GB
    "brain_aging": "4442d412-91cb-4261-acca-8adf5fa04c11",     # 291 donors, 12 GB
    "brain_wm": "c05e6940-729c-47bd-a2a6-6ce3730c4919",        # 20 donors, 0.43 GB
    "muscle": "15d374d6-0dfd-4d8e-ade7-81f73dc921ee",          # 32 donors, 1.4 GB
    "heart_hca": "364bd0c7-f7fd-48ed-99c1-ae26872b1042",       # 50 donors, 6.4 GB
    "heart_myocarditis": "fe7aae33-6f7c-41a5-8d29-9996a9ddf1ab",
    "liver_macrophage": "e84f2780-51e8-4cfa-8aa0-13bbfef677c7",
}


def resolve(name: str) -> dict:
    tab = pd.read_csv(TISSUE_DIR / "t1_cxg_primary_candidates.csv")
    did = REGISTRY[name]
    hit = tab[tab.dataset_id == did]
    if hit.empty:
        raise SystemExit(f"{name} dataset_id {did} not in fetched primary table")
    r = hit.iloc[0]
    return dict(name=name, dataset_id=did, url=r.h5ad_url, n_donors=int(r.n_donors),
                cell_count=int(r.cell_count), h5ad_bytes=int(r.h5ad_bytes),
                title=r.title, tissue=r.tissue_group)


def local_path(name: str) -> Path:
    info = resolve(name)
    dest = TISSUE_RAW / name
    dest.mkdir(parents=True, exist_ok=True)
    fname = info["url"].rstrip("/").rsplit("/", 1)[-1]
    return dest / fname


def main(names=None):
    names = names or sys.argv[1:]
    if not names:
        names = ["muscle", "brain_wm", "retina_sc"]  # smallest first if called bare
    for name in names:
        info = resolve(name)
        print(f"=== {name} ===", flush=True)
        print(f"  dataset_id={info['dataset_id']}\n  title={info['title']}\n  "
              f"donors={info['n_donors']} cells={info['cell_count']:,} bytes={info['h5ad_bytes']:,}", flush=True)
        print(f"  url={info['url']}", flush=True)
        download(info["url"], str(TISSUE_RAW / name))


if __name__ == "__main__":
    main()

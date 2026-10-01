"""Human–mouse 1:1 orthologs from NCBI Gene (gene_orthologs + gene2ensembl).

BioMart was too slow/unreliable for 25k IDs. NCBI files are the public gene_orthologs
dump; IDs are not invented. Cached under data/processed/.

Usage: python src/trajectory_orthologs.py
"""
from __future__ import annotations

import gzip
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import DATA_PROC  # noqa: E402
from download import download  # noqa: E402
from trajectory_common import TRAJ_DIR, TRAJ_SEED, Logger, dump_json, load_phase1_matrix  # noqa: E402

NCBI_ORTH = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz"
NCBI_ENS = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz"
HUMAN_TAX = 9606
MOUSE_TAX = 10090
RAW = DATA_PROC.parent / "raw" / "ncbi_gene"
CACHE = DATA_PROC / "human_mouse_orthologs_ncbi.tsv"


def _read_orthologs(path: Path, log) -> pd.DataFrame:
    log(f"[read] {path}")
    rows = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline()
        log(f"  header: {header.strip()[:120]}")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            try:
                # NCBI gene_orthologs: tax_id, GeneID, relationship, Other_tax_id, Other_GeneID
                t1, g1, rel, t2, g2 = int(p[0]), p[1], p[2], int(p[3]), p[4]
            except ValueError:
                continue
            if t1 == HUMAN_TAX and t2 == MOUSE_TAX:
                rows.append((g1, g2, rel))
            elif t1 == MOUSE_TAX and t2 == HUMAN_TAX:
                rows.append((g2, g1, rel))
    df = pd.DataFrame(rows, columns=["human_geneid", "mouse_geneid", "relationship"])
    log(f"  human-mouse rows={len(df)}")
    return df


def _read_ensembl(path: Path, tax_ids, log) -> pd.DataFrame:
    log(f"[read] {path} tax={tax_ids}")
    want = set(tax_ids)
    rows = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline()
        log(f"  header: {header.strip()[:160]}")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            try:
                tax = int(p[0])
            except ValueError:
                continue
            if tax not in want:
                continue
            # tax_id, GeneID, Ensembl_gene_identifier, RNA_nucleotide_accession, ...
            rows.append((tax, p[1], p[2].split(".")[0]))
    df = pd.DataFrame(rows, columns=["tax_id", "geneid", "ensembl"])
    df = df.drop_duplicates(["tax_id", "geneid", "ensembl"])
    log(f"  rows={len(df)}")
    return df


def run():
    log = Logger(TRAJ_DIR / "g3_orthologs.txt")
    try:
        return _run(log)
    finally:
        log.close()


def _run(log):
    log(f"TRAJECTORY orthologs  seed={TRAJ_SEED}  source=NCBI gene_orthologs")
    data = load_phase1_matrix(log)
    genes = data["genes"]
    human = [str(g).split(".")[0] for g in genes.gene_id.astype(str)]
    human_set = set(human)
    log(f"[filter] Phase-1 human Ensembl IDs: {len(human)}")

    RAW.mkdir(parents=True, exist_ok=True)
    p_orth = Path(download(NCBI_ORTH, str(RAW)))
    p_ens = Path(download(NCBI_ENS, str(RAW)))

    orth = _read_orthologs(p_orth, log)
    ens = _read_ensembl(p_ens, {HUMAN_TAX, MOUSE_TAX}, log)
    h_ens = ens[ens.tax_id == HUMAN_TAX][["geneid", "ensembl"]].rename(
        columns={"geneid": "human_geneid", "ensembl": "human_ensembl"})
    m_ens = ens[ens.tax_id == MOUSE_TAX][["geneid", "ensembl"]].rename(
        columns={"geneid": "mouse_geneid", "ensembl": "mouse_ensembl"})
    # one NCBI gene can have multiple Ensembl IDs; keep protein-coding-looking ENSG/ENSMUSG
    h_ens = h_ens[h_ens.human_ensembl.str.startswith("ENSG")]
    m_ens = m_ens[m_ens.mouse_ensembl.str.startswith("ENSMUSG")]

    merged = orth.merge(h_ens, on="human_geneid", how="inner").merge(m_ens, on="mouse_geneid", how="inner")
    log(f"[join] orthologs with both Ensembl IDs: {len(merged)}")
    merged = merged[merged.human_ensembl.isin(human_set)]
    log(f"[filter] restricted to Phase-1 human genes: {len(merged)} rows; "
        f"unique human={merged.human_ensembl.nunique()}")

    # strict 1:1 on Ensembl IDs
    h_n = merged.groupby("human_ensembl").mouse_ensembl.nunique()
    m_n = merged.groupby("mouse_ensembl").human_ensembl.nunique()
    h_ok = set(h_n[h_n == 1].index)
    m_ok = set(m_n[m_n == 1].index)
    strict = merged[merged.human_ensembl.isin(h_ok) & merged.mouse_ensembl.isin(m_ok)]
    strict = strict.drop_duplicates(["human_ensembl", "mouse_ensembl"])
    frac = len(strict) / max(len(human), 1)
    log(f"[filter] strict 1:1 Ensembl pairs: {len(strict)}  fraction_of_phase1={frac:.3f}")

    # also map human symbol -> mouse, using Phase-1 symbols
    sym = {str(g).split(".")[0]: s for g, s in zip(genes.gene_id.astype(str), genes.symbol.astype(str))}
    strict = strict.copy()
    strict["human_symbol"] = strict.human_ensembl.map(sym)
    strict.to_csv(CACHE, sep="\t", index=False)
    strict.to_csv(TRAJ_DIR / "g3_orthologs_1to1.csv", index=False)
    dump_json(TRAJ_DIR / "g3_orthologs_summary.json", dict(
        seed=TRAJ_SEED, source="NCBI gene_orthologs + gene2ensembl",
        n_human_phase1=len(human), n_strict_1to1=int(len(strict)),
        frac_phase1_strict=float(frac),
        n_orth_rows=int(len(orth)), n_joined=int(len(merged)),
        urls=[NCBI_ORTH, NCBI_ENS], cache=str(CACHE),
    ))
    log("[orthologs] wrote results/trajectory/g3_orthologs_*")
    return dict(n_strict=int(len(strict)), frac=float(frac))


if __name__ == "__main__":
    run()

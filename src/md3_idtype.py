"""ID-type checks. Fail loudly if a matrix cannot be mapped. Never guess silently."""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md3_common import log_id_type  # noqa: E402
from gtex_common import StopStep  # noqa: E402


_ENS_RE = re.compile(r"^ENS[A-Z]*[GT]\d+", re.I)
_NM_RE = re.compile(r"^(NM_|NR_|XM_|XR_)\d+", re.I)
_SYM_RE = re.compile(r"^[A-Z][A-Z0-9\-.]{1,18}$")


def detect_id_type(ids, log, tag, source):
    ids = [str(x) for x in np.asarray(ids).astype(str)]
    n = int(len(ids))
    if n == 0:
        rec = dict(tag=tag, source=str(source), n=0, kind="empty",
                   n_ensembl=0, n_refseq=0, n_symbol_like=0, examples=[])
        log_id_type(tag, rec)
        if log is not None:
            log(f"[id_type {tag}] EMPTY source={source}")
        return rec
    n_ens = int(sum(1 for x in ids if _ENS_RE.match(x.split(".")[0])))
    n_ref = int(sum(1 for x in ids if _NM_RE.match(x)))
    n_sym = int(sum(1 for x in ids if _SYM_RE.match(x.upper()) and not _ENS_RE.match(x) and not _NM_RE.match(x)))
    frac_ens = n_ens / n
    frac_ref = n_ref / n
    frac_sym = n_sym / n
    if frac_ens >= 0.8:
        kind = "ensembl"
    elif frac_ref >= 0.5:
        kind = "refseq_transcript"
    elif frac_sym >= 0.5:
        kind = "symbol"
    else:
        kind = "unknown"
    rec = dict(
        tag=tag, source=str(source), n=n, kind=kind,
        n_ensembl=n_ens, n_refseq=n_ref, n_symbol_like=n_sym,
        frac_ensembl=frac_ens, frac_refseq=frac_ref, frac_symbol_like=frac_sym,
        examples=ids[:12],
    )
    log_id_type(tag, rec)
    if log is not None:
        log(f"[id_type {tag}] kind={kind} n={n} n_ensembl={n_ens} n_refseq={n_ref} "
            f"n_symbol_like={n_sym} source={source} examples={ids[:8]}")
    return rec


def require_mappable(rec, allowed, log, step):
    if rec.get("kind") not in allowed:
        msg = (
            f"{rec.get('tag')}: ID type {rec.get('kind')!r} not in {allowed}. "
            f"n={rec.get('n')} n_ensembl={rec.get('n_ensembl')} n_refseq={rec.get('n_refseq')} "
            f"n_symbol_like={rec.get('n_symbol_like')} examples={rec.get('examples')}. "
            "Not mapping. This step stops."
        )
        if log is not None:
            log(f"[id_type STOP] {msg}")
        raise StopStep(step, msg, details=rec)
    return rec


def csr_from_npz(path):
    z = np.load(path, allow_pickle=True)
    X = sparse.csr_matrix(
        (z["data"], z["indices"], z["indptr"]),
        shape=tuple(int(x) for x in z["shape"]),
    )
    symbols = np.asarray(z["symbols"]).astype(str)
    gene_id = np.asarray(z["gene_id"]).astype(str) if "gene_id" in z.files else None
    return X, symbols, gene_id, z

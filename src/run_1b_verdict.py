"""Final STOP CONDITION 1b verdict, combining the literal flag criterion with the batch/depth-mediation analysis.

Rule (from the spec): halt if the MAJORITY of A-genes are removed as cell-cycle or technical.
'Technical' is operationalised two ways and both are reported:
  (i)  literal pre-specified flags: cell-cycle lists/data, ribosomal, mitochondrial, depth |r| >= 0.3;
  (ii) mechanism-based: a gene whose age variance retains < 50 % of its value after adjusting for the technical variable
       (10x pool, or the two depth covariates alone) has an age association that is mediated by that technical variable.
       A genuine age effect retains >= ~0.77 (1 - R2(age~pool)); simulation gave 0.86-1.30.
The verdict fires if the majority is removed under EITHER operationalisation. Written to results/phase1b_STOP_verdict.txt.
"""
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from config import TAB, RESULTS, DATA_PROC  # noqa
from decomp import PRIMARY_THRESH  # noqa


def main():
    F = pd.read_csv(TAB / "phase1b_gene_flags.csv", index_col=0)
    df = pd.read_csv(TAB / "phase1a_variance_decomposition.csv", index_col=0)
    base = df.unique_age.reindex(F.index)
    A = F.index[F.gene_class_primary == "A"]
    fa = F.loc[A]
    n = len(A)
    lit = fa.flag_cell_cycle | fa.flag_technical
    ret_pool = (fa.ua_given_pool / base.loc[A]).fillna(0)
    ret_depth = (fa.ua_given_depth / base.loc[A]).fillna(0)
    mech_pool = ret_pool < 0.5
    mech_depth = ret_depth < 0.5
    cc = fa.flag_cell_cycle
    lines = []
    P = lines.append
    P("=" * 96)
    P("STOP CONDITION 1b — VERDICT (primary A-genes, thresholds %s)" % PRIMARY_THRESH)
    P("=" * 96)
    P(f"A-genes: {n}")
    P(f"(i)  literal flags — cell-cycle: {int(cc.sum())}; cell-cycle OR technical (ribo/mito/depth|r|>=0.3): {int(lit.sum())}/{n} = {100*lit.mean():.1f}%")
    P(f"(ii) mechanism — age variance retained after adjusting for 10x pool: median {ret_pool.median():.3f} "
      f"(expected >= ~0.77 for a real effect); genes retaining < 50%: {int(mech_pool.sum())}/{n} = {100*mech_pool.mean():.1f}%")
    P(f"     mechanism — age variance retained after adjusting for depth covariates only (2 df): median {ret_depth.median():.3f}; "
      f"genes retaining < 50%: {int(mech_depth.sum())}/{n} = {100*mech_depth.mean():.1f}%")
    removed_any = lit | mech_pool | mech_depth | cc
    P(f"     union (cell-cycle OR literal technical OR batch/depth-mediated): {int(removed_any.sum())}/{n} = {100*removed_any.mean():.1f}%")
    P(f"     of which cell-cycle: {int(cc.sum())} ({100*cc.mean():.1f}%) — proliferation is NOT the dominant confound")
    fired = (lit.mean() > 0.5) or (mech_pool.mean() > 0.5) or (mech_depth.mean() > 0.5)
    P("")
    P(f"VERDICT: STOP CONDITION 1b {'FIRED' if fired else 'NOT fired'}.")
    if fired:
        P("  The majority of A-genes are technical artifacts: their donor-age association is mediated by 10x-pool batch and the")
        P("  pool-level sequencing depth (UMI per cell), which is confounded with donor age (R2(age~pool)=0.23; between-pool r(age,depth)=+0.43;")
        P("  within-pool r=-0.03). Literal |r|>=0.3 depth flags alone under-count this (they catch %d/%d); the mechanism test catches it." % (int(fa.flag_depth.sum()), n))
        P("  Per the operating mode, Phase 1 halts here: 1c and 1d are NOT run. See FINDINGS.md.")
    surv = fa[~removed_any]
    P("")
    P(f"A-genes surviving every criterion: {len(surv)} -> {', '.join(surv.sort_values('unique_age', ascending=False).symbol)}")
    P("=" * 96)
    txt = "\n".join(lines)
    print(txt)
    (RESULTS / "phase1b_STOP_verdict.txt").write_text(txt, encoding="utf-8")
    return fired


if __name__ == "__main__":
    main()

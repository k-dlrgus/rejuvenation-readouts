# PROGRESS_TOWARD

**STOP status:** DONE
**Next action:** done

## Seeds / gates

- seed `20260914`  boot `20260918`  n_perm=200  n_boot=200  n_random=200
- primary metric: cos_S = dot(d_S, u) / ||d_S||; uncalibrated R² is never a gate
- frozen ruler: context column only; not the statistic; not refit; exists=True
- MIN_CELLS=30  C3_COS_BAR=0.5  P_BAR=0.05
- YOUNG bins ('20-29', '30-39'); OLD bins ('60-69', '70-79'); middle reported, not dropped
- GSE325735: out of scope
- Refit nothing. Do not re-cluster, re-label, or re-tune Louvain or AddModuleScore.
- Donors never pooled. States never averaged.
- PREREG.flag exists=True
- anchors.npz exists=True

## Anchors (from disk)

- n_young=113 (expected 113)  n_old=231 (expected 231)  n_middle=308

## C3

- pass=True min_cos=+0.615 mean_cos=+0.775 bar=0.5
- construction_ok=True

## Reading keys (from disk)

- GM00731: key=`away_from_old_only` contrast=False
- GM23815: key=`away_from_old_only` contrast=True

## Files

- `src/toward_run.py`
- `results/toward/`
- `FINDINGS_TOWARD.md`
- `PROGRESS_TOWARD.md`


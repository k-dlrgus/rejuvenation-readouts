# PROGRESS_SAME

**STOP status:** DONE
**Next action:** done

## Seeds / gates

- seed `20260914`  boot `20260918`  n_perm=200  n_boot=200
- primary: frac_S, delta_S, Asymmetry = frac_S − frac_rev; uncalibrated R² is never a gate
- frozen ruler: context column only; not the statistic; not refit; exists=True
- MIN_HALF=500  P_BAR=0.05  f=[0.0, 0.1, 0.25, 0.5]
- GSE325735: out of scope
- Refit nothing. Do not re-cluster, re-label, or re-tune Louvain or AddModuleScore.
- Donors never pooled. States never averaged. PartialReprog pooled across timepoints only.
- PREREG.flag exists=True (not rewritten)

## Split (from disk)

- GM00731 d0 Fibroblast half A n=2469  half B n=2470
- GM23815 d0 Fibroblast half A n=3851  half B n=3852
- split_ok=True

## Positive control P

- pass=True  monotonic=True  frac_f0.50_ci_excludes_0=True  delta_f0.50<0=True

## Reading (from disk)

- key=`no_approach_same_platform`  flags=[]

## Files

- `src/same_run.py`
- `src/toward_run.py` (functions reused, not modified)
- `results/same/`
- `FINDINGS_SAME.md`
- `PROGRESS_SAME.md`


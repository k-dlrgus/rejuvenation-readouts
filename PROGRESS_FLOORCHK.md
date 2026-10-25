# PROGRESS_FLOORCHK

**Status:** COMPLETE. No STOP fired. This matches the status in FINDINGS_FLOORCHK.md.

**Primary readings:** `CLAIM1-PARTIAL: 18 of 1,836` (S1) and `CLAIM2-HOLDS` (all genes). The full list
is in `results/floorchk/readings.txt`.

| step | result |
|---|---|
| Pre-registration | `results/floorchk/PREREG.flag`, written 2026-10-25T23:17:26+09:00, before any new statistic. Rulings U1–U3 are in it verbatim. |
| Reproduction (R1–R4) | 30 gating checks, 0 mismatches (`reproduction.csv`) |
| Claim-1 floor | 212 sizes × 200 draws × 4 spaces |
| Claim-1 gate | passed in S1, S2_k20, S2_k50 and S2_k100 |
| Claim-2 floor and fit | 5 ladder sizes × 200 draws. The 1/n fit passes in all genes and fails in MD and age, which therefore use the observed 2,470-cell floor only. |
| Claim-2 gate | passed in all, MD and age genes |
| Claim-2 gap spread | 200 re-splits, descriptive |
| Readings | written; `verify` gives 0 mismatches over 30 lines |
| Diagnostic (not pre-registered) | chance baseline for LENGTH-STANDS: 5.05 expected vs 18 observed (S1) |

Commands, from the repository root, with Anaconda Python and `-I`:

```bash
python -I src/floorchk_run.py run
```

```bash
python -I src/floorchk_run.py verify
```

```bash
python -I src/floorchk_run.py diagnostic
```

`run` took about 11 minutes. It needs the git-ignored, regenerable SOUTH6 intermediates
(`results/south6/stage1.npz` and `matched_*`) and the processed GSE297234 matrices (`data/processed/md4/`).
`verify` needs only the tracked CSVs.

**Next.** Nothing more within this task. Whether and how to amend paper 1's Supplementary Note 3 is a
decision for outside the session. The preprint was not edited.

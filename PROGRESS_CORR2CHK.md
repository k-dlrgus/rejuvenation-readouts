# PROGRESS_CORR2CHK

**Status:** COMPLETE. No STOP fired. This matches the status in FINDINGS_CORR2CHK.md.

**Primary readings:** `PAPER2-DEPTH` for A1b (paper 2 paragraph 31, the Southard clause) and
`PAPER2-DEPTH` for C5 (paragraph 30, "all 80 move away from them"), plus 17
`PAPER2-NOT-AFFECTED` verdicts. The full list is in `results/corr2chk/readings.txt`.

**The two failures are not alike.** A1b has 0 of 4 units beyond their own no-effect 95th
percentile and 3 below their own median, which is what no effect looks like. C5 has 56 of 80
beyond it against 4 expected by chance, so only the universal word "all" fails there. The
findings file says so in its first section; the single pre-registered key does not.

| step | result |
|---|---|
| Step 1, inventory | 21 paper-2 statements, 22 clause-level verdicts after the Q1 split. Reported to the user in chat in two passes, before the flag. Two of the session's own groupings were wrong in the first pass and corrected in the second (C5 moved to the checked list; C7 re-grounded as having no threshold at all). |
| Four questions to the user | Q1–Q4, answered in writing and reproduced verbatim in the flag. |
| Pre-registration | `results/corr2chk/PREREG.flag`, written 2026-10-26T11:21:38+09:00, before any statistic. Part 0 discloses that the inventory had already established both outcomes. |
| Reproduction (R1–R9) | 41 gating checks, 0 mismatches; 6 further descriptive rows (`reproduction.csv`). R5 reproduced every one of AXANG's stored 5th percentiles to 2.2e-16. |
| A1b | nothing drawn; the stored baseline read at each factor's own cell count (`a1b_units.csv`) |
| C5 | AXANG's floor loop re-run at its own seed with the median and 95th percentile recorded, 161 (modality, cell count) rows, plus a fresh-seed replicate that gates nothing (`c5_baseline.csv`, `c5_units.csv`) |
| Readings | written; `verify` gives 0 mismatches over 26 lines |
| Step 3 | corrected wording for both clauses, in FINDINGS_CORR2CHK.md. No paper-2 source edited, no correction notice written. |
| Step 4 | replacement wording for paper 1's paragraph-142 bullet, in FINDINGS_CORR2CHK.md. No file edited. |

Commands, from the repository root, with Anaconda Python and `-I`:

```bash
python -I src/corr2chk_run.py run
```

```bash
python -I src/corr2chk_run.py verify
```

```bash
python -I src/corr2chk_findings.py
```

`run` took about 23 minutes, almost all of it C5's two baseline passes. It needs the
git-ignored, regenerable intermediates `results/south6/stage1.npz`,
`results/south7/basis.npz`, `results/seng2/panel_counts.npz` and `results/axang/axes.npz`,
and the raw Sengstack h5ad at `data/raw/seng/exp3_merged_not_normalized.h5ad`. `verify` needs
only the tracked CSVs. `corr2chk_findings.py` rebuilds FINDINGS_CORR2CHK.md from the tracked
tables and the flag, and computes nothing.

**One ambiguity the session resolved rather than stopping on**, because it affects no
conclusion and both readings are reported: the ruling on Q4 asked for AXANG's loop
"unchanged", for "fresh seeds", and for the stored 5th percentiles to reproduce "from the same
re-run", which cannot all hold for a resampled quantile. Flag D7 fixes the resolution — the
primary re-run keeps AXANG's seed 20260914, a fresh-seed replicate at 20261030 is reported
beside it — and deviation 4 in the findings states it plainly.

**An arithmetic slip in the flag is disclosed, not edited:** its three Part 3 counts are each
labelled "of 21" and sum to 22, because paragraph 31's sentence is one statement among the 21
and two clause-level entries once Q1 splits it. The flag is the pre-registration and stands as
written.

**Next.** Nothing more within this task. Whether and how to amend paper 2 before it is made
public, and whether to fold the paragraph-142 wording into the paper-1 revision, are decisions
for outside the session. Nothing was pushed.

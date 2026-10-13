# NOT A RESULT — computed in violation of the flag's ordering

`results/transfer/PREREG.flag` section 4 says, of the same-context qualification
gate:

> If the gate fails: STOP `unmeasurable`. Reading = UNMEASURABLE. The primary
> statistic is not computed and no cross-screen number is reported.

**The gate failed** (Sengstack side: M_same below the 0.85 bar fixed in advance;
see `../gate_summary.json`). The first pass of `src/transfer_run.py read`
nevertheless computed the primary cross-screen statistic and the cross-screen
diagnostics in the same pass, because the implementation evaluated the gate and
then ran the readout unconditionally instead of stopping. That is a deviation
from the flag, it is disclosed in `FINDINGS_TRANSFER.md`, and the session that
ran it has seen those numbers.

The files are kept here rather than deleted, because the record of what was
computed must stay. They are **not a result**:

- no pre-registered reading is drawn from them;
- no number in them is quoted in `FINDINGS_TRANSFER.md`;
- the flag forbids reporting them, and that holds whatever they say.

Contents:

- `per_factor.csv` — the per-factor cross-screen table from the first pass.
- `diagnostics.json` — the cross-screen diagnostics D1-D6 from the first pass.
- `summary_first_pass.json` — the first pass's summary, including the primary
  block.
- `read_first_pass.log` — the first pass's log, including the printed primary
  line.

After this quarantine, `src/transfer_run.py` was corrected so that the `read`
stage stops at a failed gate and computes nothing past it. The corrected stage
was re-run; the gate is deterministic at seed 20260914 and its outcome is
unchanged. That re-run is not a resumption of a stopped run and did not change
any threshold: the stop still fires, with the same key.

The within-screen gate tables (`../gate_south.csv`, `../gate_seng.csv`) are
**not** quarantined. They are the gate itself, they contain no cross-screen
comparison, and they are this task's result.

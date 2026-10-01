# PROGRESS_SOUTH5

**STOP status:** STOP `control_ambiguous`

**Next action:** The `control` column does not induce a clean non-targeting / targeting partition. Not guessing the non-targeting set from guide_target strings. A corrected control annotation (or an explicit, pre-registered mapping) would be a different task.

**Exact terminal command to run the full job:**

```
python src/south5_run.py
```

**Expected runtime:** Stops in Stage 0 / guide assignment; wall time is dominated by the md5 checksum (~1-3 min on this file) plus reading obs (~30 s).

PREREG.flag exists=True and was not rewritten by the run.
Failures recorded: 1.
Flags recorded: none.

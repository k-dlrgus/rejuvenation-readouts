# PREREG_NOTE — the unfilled "[Paste the SOUTH4 prereg text …]" line in the SOUTH5 and SOUTH6 flags

Written 2026-10-02. The flag files are timestamped pre-registration records and were not edited. This note explains one line in each of them.

## The placeholder

Both flags contain this unfilled instruction, left over from the task text:

```
[Paste the SOUTH4 prereg text from "--- STAGE 1: VECTORS AND SPACES ---"
through "--- OUTPUT ---" here verbatim.]
```

| flag | written | placeholder | SOUTH4 text reproduced |
| --- | --- | --- | --- |
| `results/south5/PREREG.flag` | 2026-09-27 | lines 72–73 | banner at line 82, text from line 86 ("--- STAGE 1: VECTORS AND SPACES ---") to line 177, the end of the file |
| `results/south6/PREREG.flag` | 2026-09-27 16:08 | lines 90–91 | banner at line 100, text from line 104 to line 203, the end of the file |

The placeholder was never filled in place. The text it asks for is reproduced further down in the same flag, under the banner "SOUTH4 pre-registration, Stages 1–5, reproduced verbatim (governed by the two changes stated above: n_perm = 5000, and the named guide columns)". Read the placeholder as a pointer to that block, not as missing content.

The reproduced block runs from "--- STAGE 1: VECTORS AND SPACES ---" through the end of SOUTH4's "--- OUTPUT ---" section ("Record every failure verbatim. Do not substitute columns or repair rows."). As a text comparison, it is identical in both flags to `FINDINGS_SOUTH4.md` lines 57–148 (5,381 characters, no differences). In the SOUTH6 flag it sits inside the SOUTH5 pre-registration, which SOUTH6 adopts for Stages 0–5 ("Identical to the SOUTH5 pre-registration, reproduced verbatim below this line …").

The reproduced SOUTH4 text names `results/south4/` paths (for example `results/south4/top20.csv`). Those are SOUTH4's own words. The output locations that govern each run are set by that flag's own `--- OUTPUT ---` section above the banner, and by `src/south5_run.py` and `src/south6_run.py`. SOUTH6 wrote to `results/south6/`.

## SOUTH4 is absent and withdrawn

- `results/south4/` does not exist. `results/south5/PREREG.flag` lines 9–12 instructed: "Delete results/south4/ entirely; its run continued past its own recorded STOP and its outputs are not trustworthy. Add one line at the top of FINDINGS_SOUTH4.md marking it superseded and withdrawn".
- `FINDINGS_SOUTH4.md` line 1 records it: "**SUPERSEDED AND WITHDRAWN by SOUTH5 (FINDINGS_SOUTH5.md).** … `results/south4/` has been deleted as untrustworthy. Nothing below is changed." Its stop key was `guide_assignment`.
- SOUTH4's own flag, `results/south4/PREREG.flag`, was deleted with the directory.

## Where the SOUTH4 pre-registration text lives

1. `FINDINGS_SOUTH4.md`, the fenced block at lines 13–149: the full SOUTH4 pre-registration, written before any statistic. This is the only complete copy, including SOUTH4's own task header, DATA and STAGE 0 sections.
2. `results/south5/PREREG.flag` lines 86–177 and `results/south6/PREREG.flag` lines 104–203: Stages 1–5 through OUTPUT, verbatim, as described above.
3. `FINDINGS_SOUTH5.md` and `FINDINGS_SOUTH6.md` repeat their flags verbatim, including the same placeholder and the reproduced block.

`src/south4_run.py` remains as code. It does not contain the pre-registration text.

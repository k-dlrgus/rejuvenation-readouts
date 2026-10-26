# PROGRESS_CORR3

**Status:** COMPLETE. No STOP fired. This matches the status in FINDINGS_CORR3.md.

Writing only. No statistic was computed; every number was copied from `FINDINGS_CORR2CHK.md`,
`FINDINGS_FLOORCHK.md` or their `results/` tables.

| step | result |
|---|---|
| 0. Confirm the paper-2 source | **v0.7 is the source of the posted PDF.** 99 of 99 body paragraphs verbatim and in order, 52 of 52 Table 1 cells present, no paragraph differing. Hyphens had to be normalised out of both sides, because `pdftotext` drops a hyphen at a line break. |
| 1. Grammar fix in v0.8 | "1 move away less far" → "1 moves away less far", paragraph 30, nothing else. `results/corr2chk/readings.txt` and `FINDINGS_CORR2CHK.md` untouched. |
| 2. `CORRECTIONS_PAPER2.md` | written, in `CORRECTIONS_PAPER1.md`'s form, stating plainly that paper 2 is unpublished so this is a pre-publication revision |
| 3. Paper 1, paragraph 142 | folded into v1.3 beside CORR1's change; `CORRECTIONS_PAPER1.md` updated in seven places; both reversible judgement calls accepted |
| 4. Repository record | paper-2 entry added to `results/verify/WITHDRAWN_READINGS.md` section D, paper-1 entry extended to paragraph 142, README pointers for both papers |

Commands, from the repository root, with Anaconda Python and `-I`:

```bash
python -I src/corr3_run.py
```

It rebuilds both revised docx files from their originals, writes both paragraph-level diffs, and
asserts that the paragraph-142 wording is byte-identical to the blockquote in
`FINDINGS_CORR2CHK.md`. It needs `lxml`, which is in the machine's Anaconda environment but not in
`requirements.txt` — the same is true of `src/corr1_run.py` and `src/corr2_run.py`.

**Verification, independent of that script.** Both source files are untouched and their sha256
values match the ones CORR1 and CORR2 recorded in their own committed diffs (`v1.2 5dfeaf0e…` at
`d8d892a`, `v0.7 74f5c339…` at `46a05c1`). In each revision only `word/document.xml` differs; other
members are byte-identical 38 of 38 and 53 of 53. Paper 1's order-aware diff shows exactly three
changes: an inserted paragraph at position 134, v1.2 paragraph 142 replaced, v1.2 paragraph 144
replaced. Paper 2 changes paragraphs 30 and 31 and nothing else.

**Two things a reader should be told, both in FINDINGS_CORR3.md:**

1. **A standing rule was departed from, on instruction, without asking first.** CLAUDE.md forbids
   writing into a previous task's results directory; this task regenerated
   `results/corr1/paper1_diff.txt` and `results/corr2/paper2_diff.txt` because each describes a
   docx that changed. Each now carries a provenance header naming CORR3 and the commit where the
   previous version can be read. Deviation 1.
2. **`SUPPLEMENTARY_NOTE_COORDSPACE.md:52` still carries the superseded paragraph-142 bullet**, in
   a tracked file headed "text for the preprint". It was left alone, following CORR1's convention,
   and is listed as row 11 of `CORRECTIONS_PAPER1.md` Appendix A. Whether to revise it is a
   decision outside this task.

**Next.** Nothing more within this task. No PDF was produced for either paper, nothing was
uploaded, and nothing was pushed. Typesetting v1.3 and v0.8, and the timing of any release, are
decisions outside the session.

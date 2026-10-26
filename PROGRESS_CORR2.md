# PROGRESS_CORR2

**Status:** COMPLETE. No STOP fired. This matches the status in FINDINGS_CORR2.md.

**What this task is.** A writing task: fold CORR2CHK's corrected wording for paper-2 paragraphs
30 and 31 into the draft, on the CORR1 pattern. **No statistic was computed.** Every number in
the new text is copied from `FINDINGS_CORR2CHK.md` or `results/corr2chk/`.

**Paper 2 is not public**, so this is a pre-publication revision of a draft, not a correction
notice. There is no `CORRECTIONS_PAPER2.md`, and no notice paragraph was inserted into the draft.
That is the one structural difference from CORR1, which revised a published paper.

| step | result |
|---|---|
| Read first | `CLAUDE.md`, `FINDINGS_CORR2CHK.md`, `CORRECTIONS_PAPER1.md`, `src/corr1_run.py`. |
| Locate the clauses | v0.7 body paragraphs 30 (C5) and 31 (A1b), by docx body-paragraph index, text-box content excluded. Both old strings occur exactly once in `word/document.xml`. |
| Revision | 3 string edits: paragraph 30, one sentence; paragraph 31, its last sentence in two pieces, either side of the superscript Southard citation run, which is kept in place. CORR2CHK's wording used verbatim. |
| Output | `paper/paper2_passage_draft_v0.8.docx`, beside an untouched v0.7. 117 body paragraphs in both; no paragraph inserted. |
| Verification | each old string unique (3 of 3); other zip members byte-identical (38 of 38); undoing the three edits reproduces v0.7's XML exactly; paragraph count unchanged. v0.7's hash taken before and after the run and unchanged. |
| Tracked diff | `results/corr2/paper2_diff.txt`, paragraph-level, with both files' sha256. |
| Consistency sweep | 21 rows in `FINDINGS_CORR2.md`, across tracked and untracked text, both manuscripts, every paper-2 draft and PDF, `bundle/`, `paper_package/`, `notebooks/` and `results/paper2_figs*/`. One row edited (the v0.8 draft); everything else not edited, with the reason per row. Nothing became inconsistent. |
| Two readings resolved rather than assumed | the Summary's "all four move them away from young donors" and paragraph 34's Figure 2 caption. The Summary statement is the *Sengstack* four-factor test, CORR2CHK id C4, ruled `PAPER2-NOT-AFFECTED` with a size-matched threshold and power 0.945–0.990; it is **not** corrected. The caption's "all lie above the horizontal line" is a call-based statement — zero on that axis *is* the threshold — so it stays consistent with the revised paragraph 30. |
| Figure provenance | `results/paper2_figs_fix/` carries "all 80 move away from GTEx young (delta_gtex_young > 0)" in its plan and READMEs. **Not edited:** a previous task's results directory, which CLAUDE.md forbids writing into, and its wording is already the sign of `delta_gtex_young`, which is the part that survives. `paper2_numbers.md` carries neither clause. |

Commands, from the repository root, with Anaconda Python (which has `lxml`; the repository's
`.venv` does not, and `src/corr1_run.py` imports it too):

```bash
python -I src/corr2_run.py
```

The script refuses to overwrite an existing `paper/paper2_passage_draft_v0.8.docx`, so delete that
file before re-running. It reads one docx and writes one docx plus
`results/corr2/paper2_diff.txt`; it touches no dataset and no other file.

**Deviations are listed in full in `FINDINGS_CORR2.md`.** The ones worth naming here: no notice
paragraph (paper 2 is unpublished); the typographic apostrophe substituted for the straight one in
CORR2CHK's blockquote; a subject-verb slip in the pre-registered wording ("1 move away less far")
kept as written and disclosed rather than improved silently; and paragraph 31 replaced as two
string edits so its citation marker stays where it was.

**Next.** Nothing more within this task. Whether to typeset a v0.8 PDF, whether to publish paper 2
at all, and whether to fold CORR2CHK's step-4 wording into the paper-1 revision are decisions for
outside the session. One local commit was made. **Nothing was pushed.**

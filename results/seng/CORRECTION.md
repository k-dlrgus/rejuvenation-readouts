# SENG — corrections

Written 2026-10-07. Housekeeping pass, no new analysis and no new statistic.
`FINDINGS_SENG.md`, `DEBUG_SENG.md` and the files under `results/seng/` and
`results/seng_why/` are **not edited**. The one exception is a single line added
at the very top of `FINDINGS_SENG.md` and of `DEBUG_SENG.md` that points to this
file. Each item below quotes the original text with its file and line, and says
what is wrong with it or how it should be read. The original text stays as
written.

**Line numbers** refer to the two files as committed in `d79ecf3`, before the
notice line was added. In the current files, add 1.

---

## C1 — "does not reproduce their direction" means the significance gate failed, not that the effect reversed

`FINDINGS_SENG.md` L136:

> - **pipeline_disagrees:** Fewer than 3 of 4 named hits have cos_g>0 with N_NT FDR q<=0.05 (1/4). Our pipeline does not reproduce their direction. Named hits are not interpreted as toward/away.

`FINDINGS_SENG.md` L630:

> Our pipeline does not reproduce their direction. Named hits are not interpreted as toward or away.

This is the pre-registered wording (`FINDINGS_SENG.md` L93, inside the verbatim
pre-registration; `results/seng/PREREG.flag` L84). It stands, and
`pipeline_disagrees` stands as the fired reading.

It should **not** be read as the effect reversing. All four named factors moved
in the direction the paper reports: each had a positive `cos_g`, that is, a
displacement aligned with the late-to-early passage direction `v`
(`DEBUG_SENG.md` L9). What failed was significance: 1 of 4 reached q ≤ 0.05 under
this project's null (N_NT, non-targeting cells drawn as fake perturbations, BH
within modality).

| named factor | cos_g | BH q (N_NT) | q ≤ 0.05 |
|---|---|---|---|
| CRA_E2F3 | +0.248 | 0.092 | no |
| CRA_EZH2 | +0.197 | 0.092 | no |
| CRI_STAT3 | +0.155 | 0.100 | no |
| CRI_ZFX | +0.316 | 0.034 | yes |

Values as stored: `DEBUG_SENG.md` L9 and L68–71, `results/seng/screen.csv`.

---

## C2 — "That is their metric" is incorrect

`FINDINGS_SENG.md` L113:

> - cos_g alone as evidence of rejuvenation. That is their metric; ours is delta.

The line sits inside the verbatim pre-registration (also
`results/seng/PREREG.flag` L104), so it stays as written.

It is incorrect. The paper's metric is **R_rej**, a Pearson correlation between
two vectors of log2 fold changes: late versus early passage, and TF perturbation
versus non-targeting control (`DEBUG_SENG.md` L83–95, which quotes the paper's
main text and SI Appendix). `cos_g` is **this project's** statistic, defined at
`FINDINGS_SENG.md` L55. The two rank the screen similarly (Spearman −0.984,
`DEBUG_SENG.md` L111, a diagnostic; see C6), but they are different statistics,
and `cos_g` should not be described as the paper's metric.

---

## C3 — "their-style hits" is this project's label, and the count should not be cited

`FINDINGS_SENG.md` L107:

> - Report: number of perturbations with cos_g q <= 0.05 ("their-style hits"),

`FINDINGS_SENG.md` L661:

> Their-style hits in this table are perturbations with cos_g BH q <= 0.05 (the screen bullet). The column n_cos_q_and_positive also requires cos_g > 0, which is the reproduction definition.

L107 sits inside the verbatim pre-registration (also `results/seng/PREREG.flag`
L98), so it stays as written.

"Their-style hits" is this project's label for perturbations with `cos_g` BH
q ≤ 0.05 under N_NT. It is **not** the paper's hit rule, which selects TFs with
R_rej ≤ −0.3 (`DEBUG_SENG.md` L83). The count reported under this label
(`FINDINGS_SENG.md` L658–659, column `n_cos_q`: CRA 0 of 167, CRI 25 of 173) was
reported without a chance baseline and **should not be cited**.

---

## C4 — the 9-of-15 and 7-of-15 figures cannot be verified and should not be cited

`DEBUG_SENG.md` L115:

> The absolute values are weaker than Table 1. Nine of the 15 CRA names in Table 1, and seven of the 15 CRI names, fall in our top 15 by `R_rej`; those same names are also high by `cos_g`.

These figures cannot be verified from tracked files, because the Table 1 names
are not tracked. The only published Table 1 values in the repository are those
of the four named factors (`src/seng_why.py` L555–557; column
`published_table1_R_rej` of `results/seng_why/rrej_vs_cos.csv`, which is empty
for the other 336 rows). The 9-of-15 and 7-of-15 figures, and the clause that
"those same names are also high by `cos_g`", **are not to be cited**.

---

## C5 — the mouse-liver statement rests only on the figshare listing

`DEBUG_SENG.md` L160:

> The passage of the perturbed cells is not in this file and not in the mouse-liver file.

The mouse-liver file was never downloaded (`DEBUG_SENG.md` L13: "No additional
file was downloaded."). The statement about it rests only on the figshare
listing: `results/seng_why/figshare_30899018.json` and
`results/seng_why/figshare_30899018_files.csv`, which list two CSV files of 233
and 1,669,439 bytes, and the record description quoted at `DEBUG_SENG.md` L18.
No contents of that file were inspected.

---

## C6 — every verdict in DEBUG_SENG.md and results/seng_why/ is a diagnostic

`DEBUG_SENG.md` uses the language of gated results: "Verdict" at L13, L37, L77
and L121; `metric_explains`, `reference_explains`, `gene_space_explains` and
`file_explains` "does not fire"; and at L149–151:

> ## Fired reading
>
> **`unexplained`.**

`results/seng_why/summary.json` carries the same keys and `"fired":
["unexplained"]`, and `results/seng_why/run_report.txt` L82 prints
`[fired] ['unexplained']`.

Every one of these is a **diagnostic, not a pre-registered reading**. No flag
file exists for `seng_why`: there is no pre-registration or flag file in
`results/seng_why/`, and none was committed elsewhere. No committed file fixes
the hypotheses, the keys or the reading `unexplained` in advance; the header of
`src/seng_why.py` also calls the run "Diagnostic only". `DEBUG_SENG.md` L3 labels the file "Diagnostic only", and that label
applies to every verdict in it.

---

## C7 — the pre-registration's order cannot be shown from the record

`FINDINGS_SENG.md` L5 and L7:

> ## Pre-registration (verbatim, written before any SENG statistic)

> Written before any SENG statistic. No threshold in this block is re-tuned after numbers exist. The block below is the task text, verbatim.

`results/seng/PREREG.flag` has **no timestamp**. It entered the repository in
commit `d79ecf3` (2026-10-02), the same commit that added `FINDINGS_SENG.md` and
the SENG results in `results/seng/`. The record therefore cannot show that the
pre-registration came before the numbers. It does not show the opposite either.
The order is simply not established by the repository.

---

## What this file does not change

No gate in `FINDINGS_SENG.md` was re-run, and no number in it or in
`DEBUG_SENG.md` was recomputed. The fired keys listed at `FINDINGS_SENG.md`
L665–667 (`screen_tracks_proliferation`, `passage_axis_not_donor_age`,
`pipeline_disagrees`) stand as recorded. C7 applies to all three: the record
cannot show that the pre-registration that defines them came first.

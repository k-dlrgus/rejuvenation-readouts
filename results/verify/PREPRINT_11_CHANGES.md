# PREPRINT_11_CHANGES — preprint_10 to preprint_11

Written 2026-10-02. `paper/rejuvenation_readouts_preprint_11.docx` was made by copying `paper/rejuvenation_readouts_preprint_10.docx` and replacing text inside `word/document.xml` only. No analysis was run and no number was recomputed. preprint_10 was not modified.

## Files

| | bytes | MD5 |
| --- | --- | --- |
| preprint_10 | 1,297,833 | `4CBE1CF551E20B8BF168D46CC793C127` |
| preprint_11 | 1,298,789 | `65F57FAD709672ABEAFB28344E802029` |

## Result

- Zip members: 55 in each, same names in the same order. The only member whose bytes differ is `word/document.xml`. Styles, numbering, headers, footers, comments, relationships, properties and media are byte-identical.
- Images: all 8 files under `word/media/` have identical SHA-256 (table below).
- `word/document.xml`: the XML outside the two target paragraphs is byte-identical. Each target paragraph was one run (one `<w:r>`, one `<w:t>`) in preprint_10 and still is. Its paragraph and run properties are unchanged, so no runs were merged. Only the text inside `<w:t>` changed.
- Body: 159 top-level paragraphs and 4 tables (126 cells) in each. Compared paragraph by paragraph and cell by cell (XML of each element): 2 differ, P89 and P94. Paragraph numbers are 1-indexed top-level body paragraphs, as in `results/verify/claims_v9.csv`. No table cell differs.
- P89 is the body paragraph under the heading "Pre-registration and decision rules" (P88). P94 begins "Analysis unit in the fibroblast atlas."

## P89 (edits B1 and B2)

Exact replacements. Applying these two to the preprint_10 text gives the preprint_11 text character for character, and nothing else in the paragraph changed:

- B1: "Where a later analysis superseded an earlier reading, the earlier text was retained and marked withdrawn rather than edited; three readings were withdrawn in this way." → "Where a later analysis superseded an earlier reading, the earlier text was retained and marked withdrawn or superseded rather than edited. Nine readings across the project are marked this way, three of which concern analyses reported here; all are listed in the repository."
- B2: "were carried out in the same project under the same protocols;" → "were carried out in the same project and pre-registered in the same way, with any deviations recorded alongside them;"

preprint_10:

> Each gated analysis was specified in writing before any number was computed, including its statistics, nulls, positive controls, decision thresholds and the exact wording of each possible conclusion. A flag file was written before execution and the pre-registration text is reproduced verbatim in the results file for each analysis. Only the pre-registered reading whose conditions fired is reported as the conclusion of that analysis. Where a later analysis superseded an earlier reading, the earlier text was retained and marked withdrawn rather than edited; three readings were withdrawn in this way. These are self-administered protocols recorded in the project repository, not deposits with an external registry. Diagnostics and descriptive summaries are not treated as gated results. Analyses of brain tissue and of a published transcription-factor screen were carried out in the same project under the same protocols; they do not bear on the question tested here and are not reported in this paper. Their pre-registrations, code and results, including a transfer that failed its gate, are in the archive.

preprint_11:

> Each gated analysis was specified in writing before any number was computed, including its statistics, nulls, positive controls, decision thresholds and the exact wording of each possible conclusion. A flag file was written before execution and the pre-registration text is reproduced verbatim in the results file for each analysis. Only the pre-registered reading whose conditions fired is reported as the conclusion of that analysis. Where a later analysis superseded an earlier reading, the earlier text was retained and marked withdrawn or superseded rather than edited. Nine readings across the project are marked this way, three of which concern analyses reported here; all are listed in the repository. These are self-administered protocols recorded in the project repository, not deposits with an external registry. Diagnostics and descriptive summaries are not treated as gated results. Analyses of brain tissue and of a published transcription-factor screen were carried out in the same project and pre-registered in the same way, with any deviations recorded alongside them; they do not bear on the question tested here and are not reported in this paper. Their pre-registrations, code and results, including a transfer that failed its gate, are in the archive.

## P94 (edit B3)

Exact change: text appended after "Seeds and resampling settings were unchanged.", inside the same run. Nothing else in the paragraph changed:

- B3: appended "This grouping was not pre-registered. The pre-registered analysis treated each of the 93 pseudobulks as a unit and gave ρ = +0.455; we report the per-person analysis because repeated samples from one person are not independent, and it gives the lower value. The ruler-minus-MD difference computed on the 93 pseudobulks is superseded by the per-person value in Table 2."

preprint_10:

> Analysis unit in the fibroblast atlas. The atlas's donor field holds one identifier per sequencing sample in several constituent studies. Each identifier was traced to its source study and the original author's donor identifier, and pseudobulk scores were averaged within each person before any correlation, permutation or bootstrap, giving 65 donors (64 after dropping pseudobulks with fewer than 30 cells, which gave the same conclusion). Seeds and resampling settings were unchanged.

preprint_11:

> Analysis unit in the fibroblast atlas. The atlas's donor field holds one identifier per sequencing sample in several constituent studies. Each identifier was traced to its source study and the original author's donor identifier, and pseudobulk scores were averaged within each person before any correlation, permutation or bootstrap, giving 65 donors (64 after dropping pseudobulks with fewer than 30 cells, which gave the same conclusion). Seeds and resampling settings were unchanged. This grouping was not pre-registered. The pre-registered analysis treated each of the 93 pseudobulks as a unit and gave ρ = +0.455; we report the per-person analysis because repeated samples from one person are not independent, and it gives the lower value. The ruler-minus-MD difference computed on the 93 pseudobulks is superseded by the per-person value in Table 2.

## Values filled in or confirmed before writing

- B1 counts, from `results/verify/WITHDRAWN_READINGS.md` after the 2026-10-02 update (section "Counts used by preprint_11 P89"): 2 marked withdrawn (A1, A2) + 7 marked superseded (B1–B7) = **9**, written "Nine". Of these, **3** concern analyses reported here (A2, B5, B7), written "three".
- B3 `93`: `results/fibro2/ta_cxg_a19d1667_result.json` has `n_samples` = 93 and `n_donors` = 93 (one pseudobulk per donor_id). Confirmed, unchanged.
- B3 `+0.455`: the same file has `rho` = 0.45512666633316906. At three decimals with an explicit sign, as the manuscript prints ρ elsewhere (e.g. "ρ = +0.433" in P47, "+0.433" in Table 1), it is +0.455. Confirmed, unchanged.
- B3 "lower": the per-person ρ is 0.43338041757870716 (`results/paper_figs/cxg_true_donor_results.csv`, `analysis=true_donor`, `contrast=ruler_spearman`, 65 donors). 0.433 < 0.455, so the statement is true.

## Images under word/media/

| member | bytes | SHA-256 preprint_10 | SHA-256 preprint_11 | identical |
| --- | --- | --- | --- | --- |
| `word/media/fig_5.png` | 122,565 | `0a5beb49843a92811157f265a89e31605ae31fa5c61525cd1cfb4f552e0837b9` | `0a5beb49843a92811157f265a89e31605ae31fa5c61525cd1cfb4f552e0837b9` | True |
| `word/media/78ebe06391334a0ed3a4bf99452cdf53b2d947a7.png` | 300,064 | `478a83e1d1e7d7ed0b0a76214a7d7dd297bdb6d1d86bb975bf5b2b6efa26efcd` | `478a83e1d1e7d7ed0b0a76214a7d7dd297bdb6d1d86bb975bf5b2b6efa26efcd` | True |
| `word/media/fig_n1.png` | 118,005 | `d5977dfca083893f7b9a408060eafeb3194fd89e314adf167fb5b5ddc0521993` | `d5977dfca083893f7b9a408060eafeb3194fd89e314adf167fb5b5ddc0521993` | True |
| `word/media/495d3938b8fcc77cde51b9ec03c8ed7cb46f9e90.png` | 200,791 | `febcfa799415904c329656ce63dbe4dce9cfec89270324aaa74f5f4b0977f36a` | `febcfa799415904c329656ce63dbe4dce9cfec89270324aaa74f5f4b0977f36a` | True |
| `word/media/65d61ce1b305cc08eaf63cc3b5d06609592f3d03.png` | 226,392 | `3fbd45410adbdf96702f08d0d2cc40ed72219df43ac079ce78348f5e4aa067a5` | `3fbd45410adbdf96702f08d0d2cc40ed72219df43ac079ce78348f5e4aa067a5` | True |
| `word/media/f504d7eff6beddab9ebd68e32921dd51381d73ae.png` | 77,467 | `c414721bce3e00a0198d6deb627d2bd9075e3344b6f3db91d1f8111e4db32c01` | `c414721bce3e00a0198d6deb627d2bd9075e3344b6f3db91d1f8111e4db32c01` | True |
| `word/media/2dde3ab0536538ed3b9714aeb12c226bc6244ca5.png` | 119,576 | `f104321a5364e8db8e7343fbb667938582a11ec77a66c31a05c76e5d29a47d71` | `f104321a5364e8db8e7343fbb667938582a11ec77a66c31a05c76e5d29a47d71` | True |
| `word/media/fig_n2.png` | 207,518 | `ddec4b3883d75fb5a1e947519aec11aa734430bfb880404a3f43116056a480cd` | `ddec4b3883d75fb5a1e947519aec11aa734430bfb880404a3f43116056a480cd` | True |

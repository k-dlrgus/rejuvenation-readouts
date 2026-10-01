# PREPRINT_12_CHANGES — preprint_11 to preprint_12

Written 2026-10-02. **This record supersedes the earlier preprint_12 build of the same day** (SHA-256 `4201a17b…03d0e7`, 1,299,063 bytes). That file was overwritten. The rebuild applies the same thirteen changes with two differences. Change 3 now uses corrected replacement text, because the earlier text made a claim that the data do not support (see "Why change 3 was corrected"). Change 13 now also converts six more British spelling families.

`paper/rejuvenation_readouts_preprint_12.docx` was made by copying `paper/rejuvenation_readouts_preprint_11.docx` and replacing text inside `word/document.xml` only. No analysis was run. No figure or image was touched. preprint_11 was not modified: its SHA-256 was the same before and after the build, and it matches the hash recorded for the earlier build. The only numbers computed were the three change-3 checks below. They were read from an existing results table, and no analysis was rerun to get them.

## Files

| | bytes | MD5 | SHA-256 |
| --- | --- | --- | --- |
| preprint_11 | 1,298,789 | `65F57FAD709672ABEAFB28344E802029` | `96e11d391ff3233a2d1614bec852509456e3f9b201e1d0c65cfdabb7f4a0b475` |
| preprint_12 (this build) | 1,299,087 | `49E1589E5D13CDBAEB4B88EA4C70DFB5` | `64105e483a9aab07d44f3d9e09751193dd2814f3aee5f0af37188db6f251107a` |
| preprint_12 (earlier build, superseded) | 1,299,063 | `D625A5457E6F2A097B84C2A8E4B3F896` | `4201a17b53d65aa21b85d6b86ba3b1f30e6e5231d4a504d77656cd724503d0e7` |

## Why change 3 was corrected

The earlier change-3 text said: "In every one of the 400 sets the cells travel more than twice the length of the donor line … and those mixtures travel only a quarter to a sixth as far". `results/newstory/partB_random_sets_per_set.csv` contradicts the first claim. With ρ = √(fwd_progress² + fwd_perp²), 8 of the 400 sets travel less than twice the donor line:

| space | n_genes | set_index | ρ |
| --- | --- | --- | --- |
| MD | 205 | 121 | 1.952 |
| AGE | 3089 | 5 | 1.989 |
| AGE | 3089 | 15 | 1.998 |
| AGE | 3089 | 119 | 1.994 |
| AGE | 3089 | 137 | 1.994 |
| AGE | 3089 | 142 | 1.966 |
| AGE | 3089 | 186 | 1.994 |
| AGE | 3089 | 192 | 1.988 |

None of these 8 ends closer. Starting one donor-line length from the target, a move of length ρ at cosine cos ends closer only if ρ < 2·cos, and no set meets that condition. So "more than twice" was not the right criterion for "cannot end closer" in these sets. The corrected text states both the true minimum length (1.95) and the direction-aware condition. The "quarter to a sixth" figure for the mixtures was never traced to any results file, and the CSV has no sideways component for the mixtures to check it against. It is dropped from the corrected text.

### The three checks, run before writing the file

Source: `results/newstory/partB_random_sets_per_set.csv`, 400 rows (MD 205 genes × 200 sets, AGE 3,089 genes × 200 sets), no missing values in the columns used.

| check | expected | found | agrees |
| --- | --- | --- | --- |
| minimum ρ over 400 rows | 1.952 (rounds to 1.95) | 1.95224 | yes |
| rows with ρ < 2·cos, where cos = fwd_progress / ρ | 0 | 0 | yes |
| rows with mix50_delta < 0 | 400 | 400 | yes |

All three agreed, so the file was written. As a consistency check, √((1 − fwd_progress)² + fwd_perp²) equals the CSV's `fwd_final_over_start` to within 1.3e-15, and that column is above 1 in all 400 rows with minimum 1.814. That matches the "closest any comes is 1.81 times its starting distance" already printed in P31.

Changes 2, 5 and 9a still say the cells travel "2.3 times" or "more than twice" the donor gap. Those statements are about the actual MD-gene displacement in P22 (progress 0.790, sideways 2.16), not about the 400 random sets, so this correction does not affect them.

## Result

- Zip members: 55 in each, same names in the same order, same per-entry date, compression method and attributes. The only member whose bytes differ is `word/document.xml`. Styles, numbering, headers, footers, comments, footnotes, endnotes, relationships, content types, properties and media are byte-identical.
- Images: all 8 files under `word/media/` have identical SHA-256 (table at the end).
- `word/document.xml`: the XML before the first body element and after the last (namespaces, `w:body` wrapper) is byte-identical, and so is every body element not listed below. The file is well-formed XML (parsed with lxml). `python-docx` is not installed in this environment, so the file was not also opened through it.
- Body: preprint_11 has 159 top-level paragraphs and 4 tables (126 cells). preprint_12 has 158 paragraphs and the same 4 tables: 2 paragraphs deleted, 1 inserted. Paragraph numbers are 1-indexed top-level body paragraphs, as in `PREPRINT_11_CHANGES.md` and `results/verify/claims_v9.csv`.
- Alignment: P1 to P123 keep their numbers. preprint_11 P124 and P125 are deleted. preprint_11 P126 to P144 become preprint_12 P124 to P142. The new reference is preprint_12 P143. preprint_11 P145 to P159 become preprint_12 P144 to P158.
- Compared element by element on the raw XML of each aligned pair, exactly 26 paragraphs differ, and they are exactly the paragraphs targeted by changes 1 to 13 (list below). No other paragraph differs. All 4 tables are byte-identical as whole elements, and 0 of their 126 cells differ (T1 25, T2 18, T3 18, T4 65 cells).
- In every changed paragraph only the characters inside `<w:t>` changed. With the `<w:t>` contents blanked, the XML of each preprint_11 paragraph equals that of its preprint_12 counterpart: same paragraph properties, same run properties, same run count (1 run each; 2 runs for the bullet P86, of which only the text run changed). No runs were merged or split.
- For every changed paragraph, applying the exact replacements below and the spelling substitutions listed for it to the preprint_11 text gives the preprint_12 text character for character.
- The inserted reference paragraph has the same XML as the Cahan entry it follows (preprint_11 P144), except for the `<w:t>` text.
- Against the earlier build: the earlier build was regenerated from its original script into a scratch file, and that file's SHA-256 matched the recorded `4201a17b…03d0e7` exactly. Compared with it, this build differs only in `word/document.xml`, and within it only in P31 (change 3) and in the 7 new spelling substitutions in P27, P61, P65, P71, P124 and P149 (preprint_12 numbering). Every word-level difference in those six paragraphs is one of the substitutions in the table under change 13.

## Changes

Paragraph numbers are given as preprint_11 → preprint_12. Changes 1, 2 and 4 to 12 are the same as in the earlier build.

**(1) P89 → P89**, the body paragraph under "Pre-registration and decision rules" (P88).
"Analyses of brain tissue and of a published transcription-factor screen were carried out in the same project and pre-registered in the same way, with any deviations recorded alongside them; they do not bear on the question tested here and are not reported in this paper." → "A published transcription-factor screen was analysed in the same project and pre-registered in the same way, with its deviation recorded alongside it; it does not bear on the question tested here and is not reported in this paper."
Change 13 then turned the new "analysed" into "analyzed" (see the spelling table), so preprint_12 reads "was analyzed in the same project".

**(2) P22 → P22**, the Results paragraph ending "Complete reprogramming goes further in the same way (pluripotent cells: progress 1.33, final distance 3.9 times the starting distance)."
Inserted immediately before that final sentence, after "(Figure 1B).": "The cells also travel a long way: 2.3 times the line between the two donors. A move longer than twice that line cannot end closer to the target whichever way it points, so these cells overshoot on length alone."
The 2.3 was supplied with the change, not computed here. It is consistent with the progress (0.790) and sideways component (2.16) already printed in the same paragraph.

**(3) P31 → P31**, under "The result does not depend on which genes are read". *Replaced in this build.*
"Yet all 400 register the 50% young-cell mixture as moving closer. Every gene set can see a real approach, and in none of them do the reprogrammed cells make one. Moving sideways is not a property of the MD gene list, nor of any list we could draw at these sizes: it is where these cells go." → "In all 400 sets the cells travel at least 1.95 times the length of the donor line, and in none of them is the move short enough to end closer given the direction it takes, so the overshoot does not depend on which genes are read. All 400 nonetheless register the 50% young-cell mixture as moving closer, so a gene set of this size is not too small to detect an approach."
The three quantities in this text were checked against `results/newstory/partB_random_sets_per_set.csv` before the file was written (see "Why change 3 was corrected"). The new text contains none of the change-13 families.

**(4) P69 → P69**, Discussion, the paragraph beginning "This is not a claim that partial reprogramming fails to rejuvenate."
"and in the data examined here the positive evidence for approach is absent where it should be strongest." → "and in the data examined here the movement is too long to end closer to a young cell of the same type, whatever its direction."

**(5) P86 → P86**, the final Limitations bullet. The bullet glyph run ("•  ") was not changed. Only the text run was replaced:
"Our evidence is that approach is not detectable in this dataset at this sample size, which is weaker than evidence that approach does not occur." → "In the same-platform test the displacement exceeds twice the donor gap, so no approach could have been registered whatever its direction; that test shows overshoot rather than a failed search for approach. The GTEx-anchored test, where an approach was attainable, gives a near-orthogonal result instead."

**(6) P30 → P30**, the paragraph beginning "One reading of Figure 1 is that MD is simply the wrong list of genes".
"3,089 of their 3,538 genes are in the reference panel" → "3,089 of their 3,538 unique gene symbols are in the reference panel (the two published lists hold 1,533 and 2,007 entries, with one duplicated symbol each)".
These numbers already appear in the manuscript: 1,533 and 2,007 in P57 and P105, and 1,312 + 1,777 = 3,089 mapped genes in P105. 1,533 + 2,007 − 2 = 3,538.

**(7) P108 → P108**, under "Nulls, intervals and multiple testing" (P106). P108 is the second body paragraph of that section and the one containing the sentence. P107 is the first.
"this happens for distance intervals in the high-dimensional displacement tests" → "this happens for distance, cosine and sideways-component intervals in the high-dimensional displacement tests".

**(8) P21 → P21**, the Results paragraph reporting the MD score drop.
"it falls from 0.515 (95% CI 0.513 to 0.517)" → "it falls from 0.515 (95% cell-resampling interval 0.513 to 0.517; these intervals describe within-sample precision, not donor-level uncertainty)".
In the document the sentence reads "In the aged donor it falls from 0.515 …", with a lowercase "it" mid-sentence. The match was on that text and the lowercase "it" was kept. Only this first interval was relabelled. The later "(95% CI 7.3 to 8.2)" in the same paragraph is unchanged.

**(9a) P71 → P71**, Discussion, the paragraph beginning "The field already uses the right form of statistic".
Appended at the end, after "this cell type, and young.": "Distance to a young control has been used before: Roux et al. report that transiently reprogrammed aged mouse cells end closer to young controls. Those cells were profiled after the reprogramming factors were withdrawn, and their displacement is small relative to the age gap; the cells tested here are assayed during OSKM expression and travel more than twice that gap."
P71 also takes the change-13 substitution "programmes" → "programs" earlier in the paragraph.

**(9b) New paragraph preprint_12 P143**, inserted immediately after the Cahan et al. CellNet entry (preprint_11 P144 → preprint_12 P142):
"17. Roux, A. E. et al. Diverse partial reprogramming strategies restore youthful gene expression and transiently suppress cell identity. Cell Syst. 13, 574–587.e11 (2022)."
- The page range uses an en dash (574–587), matching every other entry in the list (e.g. "5895–5911.e17", "903–915"). The supplied text had an ASCII hyphen. This is the only character that departs from the supplied wording.
- The paragraph and run properties are copied from the Cahan entry, so font, size, spacing and hanging indent match the rest of the list. The number is followed by a single space, as for entries 10 to 16.
- Renumbering: Cahan was the last entry (16), so Roux is 17 and no existing entry changes number. The list now runs 1 to 17 consecutively. The body cites references by author name, never by number: there are no superscript runs (`vertAlign` count 0) and no bracketed numeric citations. So there were no in-text reference numbers to update.
- The section-break paragraph that followed the Cahan entry (preprint_11 P145, now P144) is unchanged and still ends the references section.

**(10) P119 → P119**, DATA AND CODE AVAILABILITY. "are available at [GITHUB REPOSITORY URL] and archived at [ZENODO DOI]" is unchanged. Each placeholder occurs once in preprint_12. P119 differs from preprint_11 only because of the spelling change "un-normalised" → "un-normalized" (change 13).

**(11) preprint_11 P124 and P125 deleted**: the heading "COMPETING INTERESTS" (Heading1) and its paragraph "The author declares no competing interests." Neither held a section break, bookmark or comment anchor. "competing interests" no longer occurs in `document.xml`. ACKNOWLEDGMENTS now follows USE OF AI TOOLS directly.

**(12) P123 → P123**, the USE OF AI TOOLS paragraph, replaced in full:
"Analysis code was written by AI coding agents operating under the author's written pre-registration protocols, and the manuscript text and figures were drafted with the assistance of an AI language model (Claude, Anthropic). The author reviewed the analyses and the manuscript and takes full responsibility for their content." → "Analysis code was written by AI coding agents under the author's written pre-registration protocols, and the manuscript was drafted with AI assistance. The author reviewed all analyses and takes full responsibility."
The apostrophe is ASCII ('), as in the rest of the document.

**(13) British → US spelling**: 29 substitutions in 19 paragraphs. That is 28 already in preprint_11, plus 1 in the sentence that change 1 introduced. The 22 substitutions from the earlier build are unchanged. Rows marked *new* are the 7 from the families added in this build (programme, neighbour, artefact, grey, travelled/travelling, acknowledgement).

| preprint_11 ¶ | preprint_12 ¶ | before | after | |
| --- | --- | --- | --- | --- |
| P17 | P17 | Operationalising | Operationalizing | |
| P17 | P17 | behaviour | behavior | |
| P17 | P17 | behaviour | behavior | |
| P20 | P20 | labelled | labeled | |
| P22 | P22 | normalised | normalized | |
| P22 | P22 | standardised | standardized | |
| P27 | P27 | Grey | Gray | *new* |
| P33 | P33 | Relabelling | Relabeling | |
| P61 | P61 | nearest-neighbour | nearest-neighbor | *new* |
| P61 | P61 | artefact | artifact | *new* |
| P65 | P65 | normalised | normalized | |
| P65 | P65 | nearest-neighbour | nearest-neighbor | *new* |
| P71 | P71 | programmes | programs | *new* |
| P89 | P89 | analysed (text from change 1) | analyzed | |
| P91 | P91 | regularisation | regularization | |
| P91 | P91 | minimising | minimizing | |
| P91 | P91 | unit-normalised | unit-normalized | |
| P92 | P92 | TMM-normalised | TMM-normalized | |
| P96 | P96 | TMM-normalised | TMM-normalized | |
| P102 | P102 | log-normalisation | log-normalization | |
| P102 | P102 | labelling | labeling | |
| P102 | P102 | labellings | labelings | |
| P104 | P104 | log-normalised | log-normalized | |
| P105 | P105 | TMM-normalised | TMM-normalized | |
| P105 | P105 | normalisation | normalization | |
| P119 | P119 | un-normalised | un-normalized | |
| P126 | P124 | ACKNOWLEDGEMENTS | ACKNOWLEDGMENTS | *new* |
| P127 | P125 | analysed | analyzed | |
| P150 | P149 | travelled | traveled | *new* |

P27 is the Figure 1 caption ("Gray points are young-cell mixtures"). P61 is the Discussion paragraph on the ruler's training distribution. P65 is the Figure 5 caption. P126 → P124 is the ACKNOWLEDGMENTS heading (Heading1, the only change to it is the text). P127 → P125 is the acknowledgments paragraph ("datasets analyzed here"). P150 → P149 is the paragraph beginning "Along the full trajectory" ("The distance traveled is not monotonic").

Scope and exclusions:
- The search covered all body paragraphs, captions and all 126 table cells. No table cell contained any of the listed families.
- "programme" was matched as the whole word "programme"/"programmes" only. A substring match would also hit "reprogrammed", "non-reprogrammed" and "reprogramming", which are US spelling already: 54 occurrences of "reprogramm…" in preprint_11 body text, including table cells T4C7, T4C22, T4C27, T4C42, T4C47 and T4C62. None of them was changed. The count is 56 in preprint_12. The change-3 replacement removes 1 ("the reprogrammed cells make one"), the change-9a sentence adds 2 ("reprogrammed aged mouse cells", "reprogramming factors"), and the Roux reference adds 1 ("partial reprogramming"). Only P71's "programmes" was converted.
- "grey" was matched as a whole word: the only occurrence is "Grey" in P27. There was no "travelling" and no lowercase or singular "acknowledgement" anywhere in the body.
- In `document.xml` as a whole, each of the new families occurs only inside `<w:t>` text, never in bookmark names (0 bookmarks), field codes (0 `instrText`) or attributes.
- No occurrence was found of modell*, colour*, centre*, favour*, summaris*, characteris*, recognis*, maximis*, utilis* or emphasise/emphasised/emphasising. "emphasis*" was read as the verb forms, because the noun "emphasis" is already US spelling. The noun does not occur either.
- The noun "analysis/analyses" was not touched. Body counts went from 24 to 23 only because change 1 removed the sentence beginning "Analyses of brain tissue".
- Reference list (preprint_11 P129–P144, preprint_12 P127–P143): left alone, and it contains none of the families anyway. "normalization" in the Robinson & Oshlack title is already US spelling.
- None of the substituted words is a gene symbol, accession, software package name or author name. "TMM-", "log-", "unit-" and "nearest-" are kept, and only the British part of the word changed.
- There is no text in quotation marks in the body. The only apostrophes are possessives.
- After the build, a re-scan of the preprint_12 body found 0 remaining occurrences of any listed family, old or new. Headers, footers, footnotes, endnotes, comments and every other XML part contain none of the new families either. They were not edited and are byte-identical to preprint_11.

## Paragraphs that differ (aligned preprint_11 → preprint_12)

P17, P20, P21, P22, P27, P30, P31, P33, P61, P65, P69, P71, P86, P89, P91, P92, P96, P102, P104, P105, P108, P119, P123, P126 → P124, P127 → P125, and P150 → P149. In addition, P124 and P125 were deleted and P143 was inserted. Content changes 1–12 landed in P21, P22, P30, P31, P69, P71, P86, P89, P108, P123 and the new P143. The remaining paragraphs differ only by the spelling substitutions above. Every other paragraph, all 4 tables and all 126 cells are byte-identical.

## Images under word/media/

| member | bytes | SHA-256 preprint_11 | SHA-256 preprint_12 | identical |
| --- | --- | --- | --- | --- |
| `word/media/fig_5.png` | 122,565 | `0a5beb49843a92811157f265a89e31605ae31fa5c61525cd1cfb4f552e0837b9` | `0a5beb49843a92811157f265a89e31605ae31fa5c61525cd1cfb4f552e0837b9` | True |
| `word/media/78ebe06391334a0ed3a4bf99452cdf53b2d947a7.png` | 300,064 | `478a83e1d1e7d7ed0b0a76214a7d7dd297bdb6d1d86bb975bf5b2b6efa26efcd` | `478a83e1d1e7d7ed0b0a76214a7d7dd297bdb6d1d86bb975bf5b2b6efa26efcd` | True |
| `word/media/fig_n1.png` | 118,005 | `d5977dfca083893f7b9a408060eafeb3194fd89e314adf167fb5b5ddc0521993` | `d5977dfca083893f7b9a408060eafeb3194fd89e314adf167fb5b5ddc0521993` | True |
| `word/media/495d3938b8fcc77cde51b9ec03c8ed7cb46f9e90.png` | 200,791 | `febcfa799415904c329656ce63dbe4dce9cfec89270324aaa74f5f4b0977f36a` | `febcfa799415904c329656ce63dbe4dce9cfec89270324aaa74f5f4b0977f36a` | True |
| `word/media/65d61ce1b305cc08eaf63cc3b5d06609592f3d03.png` | 226,392 | `3fbd45410adbdf96702f08d0d2cc40ed72219df43ac079ce78348f5e4aa067a5` | `3fbd45410adbdf96702f08d0d2cc40ed72219df43ac079ce78348f5e4aa067a5` | True |
| `word/media/f504d7eff6beddab9ebd68e32921dd51381d73ae.png` | 77,467 | `c414721bce3e00a0198d6deb627d2bd9075e3344b6f3db91d1f8111e4db32c01` | `c414721bce3e00a0198d6deb627d2bd9075e3344b6f3db91d1f8111e4db32c01` | True |
| `word/media/2dde3ab0536538ed3b9714aeb12c226bc6244ca5.png` | 119,576 | `f104321a5364e8db8e7343fbb667938582a11ec77a66c31a05c76e5d29a47d71` | `f104321a5364e8db8e7343fbb667938582a11ec77a66c31a05c76e5d29a47d71` | True |
| `word/media/fig_n2.png` | 207,518 | `ddec4b3883d75fb5a1e947519aec11aa734430bfb880404a3f43116056a480cd` | `ddec4b3883d75fb5a1e947519aec11aa734430bfb880404a3f43116056a480cd` | True |

## Amendment, 2026-10-02 05:36: P30 parentheses

After the build above (SHA-256 `64105e48…51107a`), one further edit was made to `paper/rejuvenation_readouts_preprint_12.docx`, in `word/document.xml` only. It removes the nested parentheses that change 6 created in P30.

"(3,089 of their 3,538 unique gene symbols are in the reference panel (the two published lists hold 1,533 and 2,007 entries, with one duplicated symbol each))" → "(3,089 of their 3,538 unique gene symbols are in the reference panel; the two published lists hold 1,533 and 2,007 entries, with one duplicated symbol each)"

So P30 now has one pair of parentheses, and the semicolon that follows it in the sentence ("…one duplicated symbol each); all 23,485 genes of the panel; …") is unchanged. The wording and numbers of change 6 are unchanged. Only the inner opening parenthesis became "; " and one of the two closing parentheses was removed.

| | bytes | MD5 | SHA-256 |
| --- | --- | --- | --- |
| preprint_12 before amendment | 1,299,087 | `49E1589E5D13CDBAEB4B88EA4C70DFB5` | `64105e483a9aab07d44f3d9e09751193dd2814f3aee5f0af37188db6f251107a` |
| preprint_12 after amendment (current) | 1,299,087 | `A5740861E87EE0EFBFFFC139375FDDA8` | `bfaf1b1f68bfa29c3a16cd2940f9f97f5f8919b9f4ea4f9d7ca50614b7a8e36f` |

The two files are the same size by coincidence of compression, but their contents differ.

Verification against the file as it was before the amendment:
- Zip members: 55 in each, same names, order, per-entry date, compression method and attributes. `word/document.xml` is the only member whose bytes differ.
- Images: all 8 files under `word/media/` have identical SHA-256, so the image table above still holds.
- `word/document.xml`: the XML before and after the body is byte-identical. Of the 163 top-level body elements (158 paragraphs, 4 tables, 1 section properties element), only P30 differs. With the `<w:t>` contents blanked, P30 is identical (1 run each), and its text equals the earlier text with exactly this replacement applied. The file is well-formed XML (lxml).
- preprint_11 is unchanged (SHA-256 `96e11d39…a0b475`).

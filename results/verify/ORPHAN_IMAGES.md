# ORPHAN_IMAGES — images inside preprint_10 that no paragraph displays

Written 2026-10-02. The docx was opened read-only as a zip archive and was not modified.

Source: `paper/rejuvenation_readouts_preprint_10.docx`. `word/media/` holds 8 images. `word/_rels/document.xml.rels` relates all 8 to the document. `word/document.xml` references (`r:embed`) only 6 of them: Figures 1–5 and Supplementary Figure S1. No other part of the package refers to media. The two relationships never referenced are `rId36` and `rId39`.

Both files were extracted byte-for-byte to `paper/orphan_images/`, and each extracted copy was checked equal to its zip entry.

| file in `paper/orphan_images/` | docx entry | relationship | bytes | SHA-256 | pixels (w × h) | mode | DPI | print size at that DPI |
| --- | --- | --- | ---: | --- | --- | --- | --- | --- |
| `f504d7eff6beddab9ebd68e32921dd51381d73ae.png` | `word/media/f504d7eff6beddab9ebd68e32921dd51381d73ae.png` | `rId36` | 77,467 | `c414721bce3e00a0198d6deb627d2bd9075e3344b6f3db91d1f8111e4db32c01` | 1080 × 870 | RGBA | 300 | 3.60 × 2.90 in |
| `78ebe06391334a0ed3a4bf99452cdf53b2d947a7.png` | `word/media/78ebe06391334a0ed3a4bf99452cdf53b2d947a7.png` | `rId39` | 300,064 | `478a83e1d1e7d7ed0b0a76214a7d7dd297bdb6d1d86bb975bf5b2b6efa26efcd` | 2100 × 1680 | RGBA | 300 | 7.00 × 5.60 in |

Both are PNG. The zip entries are dated 2026-09-29 18:35:34.

## What they are

The same bytes (same SHA-256) are in the two older drafts:

- `f504d7…png` = `paper/rejuvenation_readouts_preprint_3.docx` `word/media/image2.png`. In preprint_3 it is displayed as **Figure 2**, captioned "The frozen ruler reads cells outside its training distribution as young. Each point is one donor × cluster × timepoint pseudobulk from the reprogramming dataset with at least 20 cells (n = 22). …". The plot shows frozen-ruler score against nearest-neighbour distance to GTEx training donors (PCA, k = 50), annotated "ρ = −0.625, 95% CI [−0.834, −0.320], n = 22 pseudobulks". preprint_10 describes the same comparison as Figure 5B (caption P65) in the displayed `fig_5.png`.
- `78ebe0…png` = preprint_3 `word/media/image5.png`. In preprint_3 it is displayed as **Figure 5**, captioned "Same-platform test against a young target measured in the same experiment. …". It has four panels: (A) forward test, aged cells toward the young target; (B) reverse test; (C) positive control progress against the young fraction mixed in; (D) positive control Δ distance.
- `paper/rejuvenation_readouts_preprint_9.docx` carries both files under the same names and relationship IDs, and does not display them either.

So these are preprint_3 figures whose relationships were carried into preprint_9 and preprint_10 after the figures themselves were replaced. They are not displayed when the document is opened. They do travel inside the file, and anyone unzipping it can read them.

A caution for reuse: the plotted values in `78ebe0…png` are from the preprint_3 analysis and do not all match preprint_10's text. Read off the plot, the 50% mixture's progress is about 0.44 (panel C) and its Δ distance about −34.5 (panel D). preprint_10 P23 reports 0.388 and −5.143. Do not use this image as a figure for preprint_10.

Apart from the extracted copies in `paper/orphan_images/` and the older drafts named above, no file in the working tree is byte-identical to either image. Before this extraction, `bundle/MANIFEST.md` section 4 recorded "identical file in repo: none found".

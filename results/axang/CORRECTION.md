# AXANG — corrections

Written 2026-10-07. Housekeeping pass, no new analysis. The original
`d1_angles.csv`, `d1_summary.json`, `d2_summary.json` and `FINDINGS_AXANG.md`
are **not edited**. Every value this file replaces is reproduced below and
marked **superseded**, per CLAUDE.md: "A reading that is superseded is marked
superseded, never edited. The old text stays."

Corrected values live in `results/axang/d1_angles_corrected_C1lu.csv` and
`results/axang/correction_values.json`; the recomputation is
`src/axang_correction.py`.

---

## C1 — the Lu / GSE297234 reference-panel coverage row was wrong

**AXANG reported 16,781 present / 6,704 absent. That is wrong. The correct
figures are 20,452 present / 3,033 absent, and the preprint's 20,452 is right.**

### Why it was wrong

AXANG's D2 coverage test matched ruler genes to the Lu feature list **by gene
symbol only** (`src/axang_run.py`, `d2`: `luset = set(... zl["symbols"] ...)`).
The project's standard alignment is `ruler_column_index`
(`src/seng_run.py:2125`), which matches **Ensembl ID first and falls back to
symbol**. The Lu count pack carries Ensembl IDs for all 36,601 features, so
symbol-only matching discarded most of the available mapping:

| matching rule | ruler genes present | absent |
|---|---|---|
| symbol only (**what AXANG used**) | 16,781 | 6,704 |
| Ensembl first, symbol fallback (**project standard**) | **20,452** (20,443 by Ensembl + 9 by symbol) | **3,033** |

The error should have been caught when the Lu row came out **identical** to the
Sengstack row. That coincidence was not a coincidence: both datasets use the
same 36,601-feature 10x human reference, so a symbol-only match gives the same
overlap for both. The Sengstack figure of 16,781 is **correct** because the
Sengstack h5ad carries no Ensembl IDs at all — SENGINV established that
(0 of 36,601 entries match `^ENSG\d{11}`) — so for that dataset symbol matching
is the only rule available and Ensembl-first changes nothing. The Southard
figure of 20,264 is also correct: it was read from `results/south6/stage1.npz`
key `matched_ruler`, which SOUTH6 built with the standard alignment.

So **only the Lu row was affected**, and only because that one dataset had
Ensembl IDs that AXANG failed to use.

### Superseded and corrected values

**Coverage table, D2 (section 6 of FINDINGS_AXANG.md):**

| quantity | superseded | corrected |
|---|---|---|
| Lu panel genes present | 16,781 | **20,452** |
| Lu panel genes absent | 6,704 | **3,033** |
| % of the GTEx age gap's squared length in genes Lu does not measure | 21.97% | **8.97%** |

The corrected Lu figure changes the D2 conclusion's range: the share of the
donor-age direction that is unreachable by construction is **9.0% (Lu), 10.0%
(Southard), 22.0% (Sengstack)**, not "between a tenth and a fifth" as
FINDINGS_AXANG.md section 6 and section 7 state. The correct statement is
**between a tenth and a fifth for Southard and Sengstack, and under a tenth for
Lu** — Lu is the best-covered of the three, not the worst.

### `C1_lu` was built with the same faulty map, so every cosine involving it moves

`lu_c1` in `src/axang_run.py` used the same symbol-only mapping to scatter the
Lu quartile pseudobulks into the 23,485 ruler frame, so the axis itself was
built on 16,781 coordinates instead of 20,452. Rebuilt with the standard
alignment, the old and new axes are **cos 0.913871 apart (23.95°)** — the error
moved the axis materially, even though it leaves every qualitative reading
unchanged.

**S2A_k20** (the primary space):

| pair | cos superseded | cos corrected | angle from parallel, corrected | % removed, corrected |
|---|---|---|---|---|
| A1 (GTEx donor age) vs C1_lu | −0.3264 | **−0.3361** | 70.36° | 5.8% |
| A2 (frozen ruler) vs C1_lu | −0.1395 | **−0.1313** | 82.45° | 0.9% |
| A3 (Sengstack passage) vs C1_lu | +0.8421 | **+0.8486** | 31.94° | 47.1% |
| **A4 (Lu donor pair) vs C1_lu** | **+0.5355** | **+0.5304** | **57.97°** | **15.2%** |
| C2 (Tirosh gene list) vs C1_lu | +0.6907 | **+0.6904** | 46.34° | 27.7% |
| C1_south vs C1_lu | +0.9478 | **+0.9522** | 17.79° | 69.5% |
| C1_seng_CRA vs C1_lu | +0.8814 | **+0.8877** | 27.41° | 54.0% |
| C1_seng_CRI vs C1_lu | +0.8034 | **+0.8133** | 35.58° | 41.8% |
| C1_gtex vs C1_lu | +0.6346 | **+0.6469** | 49.69° | 23.7% |
| C1_gtex_anchors vs C1_lu | +0.6614 | **+0.6744** | 47.59° | 26.2% |

**Full-gene space** (23,485 dimensions; not comparable to the table above):

| pair | cos superseded | cos corrected | angle from parallel, corrected | % removed, corrected |
|---|---|---|---|---|
| A1 vs C1_lu | −0.1287 | **−0.1316** | 82.44° | 0.9% |
| A2 vs C1_lu | −0.0225 | **−0.0232** | 88.67° | 0.03% |
| A3 vs C1_lu | +0.3585 | **+0.3278** | 70.87° | 5.5% |
| A4 vs C1_lu | +0.2612 | **+0.2636** | 74.71° | 3.5% |
| C2 vs C1_lu | +0.2377 | **+0.2175** | 77.44° | 2.4% |
| C1_south vs C1_lu | +0.6259 | **+0.5961** | 53.41° | 19.7% |
| C1_seng_CRA vs C1_lu | +0.4484 | **+0.4093** | 65.84° | 8.8% |
| C1_seng_CRI vs C1_lu | +0.4612 | **+0.4214** | 65.08° | 9.3% |
| C1_gtex vs C1_lu | +0.3706 | **+0.3745** | 68.01° | 7.3% |
| C1_gtex_anchors vs C1_lu | +0.3841 | **+0.3882** | 67.16° | 7.8% |

Per-donor detail after correction: GM23815 7,782 day-0 cells, 1,946 per
quartile; GM00731 5,021 day-0 cells, 1,256 per quartile; both map 20,452 ruler
genes (20,443 by Ensembl, 9 by symbol fallback); 94 Tirosh genes found in both.

### Which AXANG results are affected, and which are not

**No AXANG gate is affected.**

- **D1's gate** is `cos(A1, C1_gtex)` = −0.8108, 35.83° from parallel,
  INTERMEDIATE. `C1_lu` does not enter it. Unchanged.
- **D1's C2 comparison** is `cos(A1, C2)` = −0.6588, SEPARABLE. Unchanged.
- **D1's triangle checks** use A1, A3, `C1_seng_CRA` and `C1_gtex`. Unchanged,
  both still consistent.
- **D3a** (QUALIFIED, p25 = 0.945) and **D3b** (NO APPROACH, 0 of 4, 5 of 340)
  do not involve `C1_lu` or the Lu dataset at all. Unchanged.
- **Stage 0** does not involve `C1_lu`. Unchanged, still PASS.

**Affected and corrected above:** the Lu row of the D2 coverage table and its
unmeasured-gap share; the D2 sentence about the range across the three datasets;
every cosine involving `C1_lu` in both spaces; and the A4 row of the one-scale
summary tables in the headline and in section 8, where the angle changes from
57.63° to **57.97°** and the share removed from 15.5% to **15.2%**.

The one-scale ordering that is AXANG's result is unchanged: passage aging
19.66°, GTEx donor age 35.83°, Lu donor pair **57.97°**, frozen ruler 85.51°.

---

## C2 — the S2A_k20 donor count: 308 is correct, and no file was wrong

**The true number is 308 GTEx donors aged 40–59.** Verified three independent
ways:

| source | value |
|---|---|
| `results/south7/basis.npz` key `n_mid` | 308 |
| `results/south7/basis.npz` key `mid_mask`, summed | 308 |
| `results/toward/anchors.npz` key `n_middle` | 308 |
| GTEx age bins in the frozen pack: 40-49 = 106, 50-59 = 202 | 308 |

The code agrees: `anchor_free_basis` (`src/south7_run.py:117`) selects
`mid = np.isin(age, list(MID_BINS))` with `MID_BINS = ("40-49", "50-59")`
(`src/south7_run.py:53`), which is exactly those two bins.

**Neither findings file is wrong.** `FINDINGS_SENG3.md` says 308 in all four
places it appears (lines 316, 326, 517, 548), `results/seng3/AMENDMENTS.md` says
308, and `FINDINGS_AXANG.md` says 308 in both places (lines 408, 510). There is
no 319 anywhere in the committed record.

The figure 319 appeared once, in a chat message during the SENG3 session, in the
preamble to the question that asked how the basis should be built. It was a
guess made before the stored `n_mid` had been read, and it was corrected to 308
in the same session before any file was written — the correction is itself
recorded in `results/seng3/AMENDMENTS.md` under B1, which notes the donor count
as "308 donors, `n_mid` in the stored basis". So the premise that SENG3
"describes it as 319" does not hold for the record; only the transient chat
message said so.

Logged here because the discrepancy was raised and because the check is cheap
to reproduce: nothing in the repository needs changing for it.

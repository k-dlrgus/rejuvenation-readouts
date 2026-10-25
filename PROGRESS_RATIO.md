# PROGRESS_RATIO

**STOP status:** complete for what was workable. No pre-registered STOP fired (`reproduction_failed`,
`flag_ambiguous` and `control_pool_empty` all did not). GSE314342 was declared **not workable** under the
flag's D13(b), which is not a STOP. Per the task ("complete what is workable, and stop"), the run ended
after GSE278572.

**Outcome (pre-registered wording; full measured-gene space only, because paper 2's compressed space
cannot be built for T cells):**

- `RATIO-OPEN: in GSE278572 Treg, 19 of 19 perturbations have ratio < 2 in full measured-gene space; arrival was geometrically attainable for them.` [within_floor = 9]
- `RATIO-OPEN: in GSE278572 Teff, 17 of 17 perturbations have ratio < 2 in full measured-gene space; arrival was geometrically attainable for them.` [within_floor = 11]
- RATIO-BLOCKED and RATIO-UNMEASURABLE did not fire in either sub-screen.

**The caveat that matters:** the size-matched control floor sits just under the perturbation ratios
(median ratio / own-size P95 ≈ 1.03 and 1.02). At this screen's depth (~900–1,150 cells per
perturbation), ratio < 2 is close to guaranteed by depth. At 30 cells, noise alone would read 1.5–1.8,
and a full-gap move 2.4–2.6. The task's premise holds for SOUTH6 in S1 (1,836 of 1,836 ≥ 2), but not
in its compressed spaces (1–2 of 1,836).

**Single next action:** none in this task. FINDINGS_RATIO.md names two follow-ups:
(1) a size-matched ratio floor for SOUTH6 S1;
(2) GSE314342 where 171 GB is workable.

## Record of what was decided when (2026-10-25, +09:00)

| time | event |
|---|---|
| before 19:20 | Read CLAUDE.md, SCOPING_DISEASE.md, FINDINGS_SOUTH6.md (lines 1–511, including stored Stage 4 medians), FINDINGS_SENG3.md B1, SCOPEDIS dataset rows, code. Fetched GEO FTP listings and sizes. |
| ~19:21 | Downloaded 3 small GSE278572 metadata files; read structure only. |
| ~19:23 | **Asked the user 4 questions** (GSE92872 arm contradiction, firing rules, floor test, full-space scaling) and stated that the compressed space cannot be built for T cells. Rulings U1–U4. |
| 19:24:21 | `results/ratio/PREREG.flag` written. Not edited afterwards. |
| 19:26–20:57 | GSE278572 matrix downloaded (3.86 GB, single stream ~0.7 MB/s; parallel range requests were throttled by NCBI). |
| 19:31 | Reproduction: 6 gating checks, 0 mismatches. R1b premise check recorded. |
| 19:32 | Synthetic smoke test of the analysis (scratchpad only). |
| 21:00–21:16 | Two streaming passes over 1,027,798,128 entries; stream output cached. |
| 21:17 | **First analysis attempt crashed** (Windows pickling, OSError 22) at the first floor launch. Only cell counts and the Treg ‖gap‖ had been logged. |
| 21:18 | Infrastructure fix (pools passed by file; 8 workers); rerun from cache; ‖gap‖ identical. |
| 21:33 | GSE278572 complete: both sub-screens RATIO-OPEN. `verify`: 0 mismatches. |
| 21:35 | GSE314342 workability, attempt 1: wrong URL path, error pages, invalid. Kept, not used. |
| 21:35 | Labelled non-pre-registered diagnostic computed from tracked CSVs. |
| 21:35–22:07 | GSE314342 workability, attempt 2: 1,917,820,392 of 1,917,820,392 bytes in 1,885 s = 1.02 MB/s → ~46.7 h > 24 h. **Not workable.** Files deleted unopened. Task stopped. |

## Commands

```bash
python -I src/ratio_run.py reproduce
```

```bash
python -I src/ratio_run.py gse278572
```

```bash
python -I src/ratio_run.py verify
```

```bash
python -I src/ratio_run.py diagnostic
```

`gse278572` needs the four GSE278572 files in `data/raw/ratio/gse278572/` (URLs in the script header).
Streaming takes ~17 min the first time, and the floor and positive control ~15 min on 8 workers. `verify`
and `diagnostic` need only tracked files in `results/ratio/`.

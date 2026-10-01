# Path scrub

Removes the absolute repo root, which contained the local Windows account name, from text files in the working tree. Text edit only. No analysis was rerun, no file was regenerated, nothing was committed.

Scope: every file on disk under the repo root, whether tracked, untracked, or git-ignored, except `.git/`, `.venv/` and `data/raw/`. Files were searched as raw bytes (UTF-8 and UTF-16LE), and zip containers (`.zip`, `.docx`) were also searched member by member.

In this file the account name is written `<USER>` and the drive as `{drive}`, so that the report does not reintroduce the string it removes.

## Replacement rule

Each occurrence of the absolute repo root is replaced by the literal token `<repo>`. Bytes before and after the match are unchanged, so `<root>\results\x.csv` becomes `<repo>\results\x.csv`, and the original separator style and escaping of the rest of the path are kept.

Root components, in order: `{drive}` = `C:`, `Users`, `<USER>`, `Desktop`, `Immortality Game`, `Age` U+2013 `Identity Separability Benchmark`.

UTF-8 and ASCII files. Byte regex, case-insensitive, with one separator style per match:

```
C(?::|%3A)(?P<s>\\+|/|%5C|%2F)Users(?P=s)<USER>(?P=s)Desktop(?P=s)Immortality(?: |%20)Game(?P=s)Age(?:\xe2\x80\x93|\\+u2013|%E2%80%93)Identity(?: |%20)Separability(?: |%20)Benchmark
```

The separator `s` covers a raw backslash, JSON-escaped `\\`, JSON-in-JSON `\\\\`, `/`, and URL-encoded `%5C` / `%2F`. The en-dash covers literal UTF-8, `\u2013` at any escape depth, and `%E2%80%93`. The replacement is `<repo>` in UTF-8.

UTF-16LE files (BOM `FF FE`, PowerShell `>` captures). Literal match at even byte offsets, replaced by `<repo>` in UTF-16LE:

- `{drive}\Users\<USER>\Desktop\Immortality Game\Age` U+2013 `Identity Separability Benchmark`
- `{drive}\Users\<USER>\Desktop\Immortality Game\Age?` U+BC12 `dentity Separability Benchmark`. This is the form actually present: a cp949 round trip turned the en-dash and the `I` into `?` U+BC12.
- Both of the above with `/` as the separator and with a lowercase `c:` drive. Neither variant occurs.

One special case: `bundle/MANIFEST.md` embeds a raw-byte preview of `results/south6/smoke_stdout.txt`. That preview is UTF-16LE, lossily decoded, and cut off at `...Separability Se` before an ellipsis, so the root inside it is truncated. The truncated UTF-16LE root bytes, up to the ellipsis, were replaced by `<repo>` in UTF-16LE (1 occurrence).

Not changed: any other string. That includes other paths under the user's home directory, the email address, the `.docx`, `.zip`, `.pyc` and `.h5ad` binaries, and line endings.

JSON check: every changed `.json` and `.ipynb` file parsed before and after. The parsed object after equals the parsed object before with the same rule applied to every string, keys included. All of them passed. In the whole tree, only two `.json` files fail to parse: the two copies of `results/survey/opened/encode_search_limit0.json`. Both are saved HTML pages, contain no path, and were not touched.

## Before

Files that contained the user path (`{drive}`, then `Users`, then `<USER>`) in any encoding: **744**, with **5655** occurrences. Separately, 15 more files contained `<USER>` only as a substring of the correspondence email or as a bare word (the 8 `.docx`, plus 7 text files). One `.bin` was a scanner false positive: it has a zip signature but is not a zip and contains no path.

| file type | files | occurrences |
|---|---:|---:|
| `.json` | 404 | 1028 |
| `.txt` | 102 | 640 |
| `.pyc` | 93 | 93 |
| `.md` | 91 | 1651 |
| `.log` | 20 | 106 |
| `.csv` | 16 | 258 |
| `.h5ad` | 5 | 5 |
| `.err` | 4 | 8 |
| `.ipynb` | 4 | 14 |
| `.py` | 4 | 4 |
| `.zip` | 1 | 1848 |
| `.docx` | 0 (8 with email only) | 0 |

By git status: 635 untracked, 98 ignored, 11 tracked.

Encodings of the repo root seen (before):

| form | occurrences replaced |
|---|---:|
| raw backslash, literal en-dash | 2384 |
| `\\` separators, `\u2013` (JSON) | 1026 |
| `\\` separators, literal en-dash | 218 |
| UTF-16LE, mangled en-dash | 64 |
| `\\\\` separators, `\u2013` (JSON inside a JSON string) | 2 |
| UTF-16LE truncated preview in `bundle/MANIFEST.md` | 1 |
| **total** | **3695** |

No forward-slash or URL-encoded form of the root occurs anywhere.

## After

Files changed: **641**. Occurrences replaced: **3695**. Repo root left in text files: **0**.

The full re-scan of the same scope searched for the username `<USER>`, for `{drive}` followed by `Users` (any separator or encoding), and, in text files, for any absolute drive path `X:\` or `X:/`. Remaining hits are below; none of them is the repo root in a text file.

Text files, not the repo root, left as is:

| what | occurrences | files |
|---|---:|---|
| PowerShell temp script path `{drive}\Users\<USER>\AppData\Local\Temp\ps-script-*.ps1` (UTF-16LE notebook-exec logs) | 8 | `results/nb_phase0_exec.err`, `results/nb_phase1_exec.err`, `results/rerun/nb_stageA_exec.log`, `results/rerun/nb_stageB_exec.log` (all tracked), and the 4 `bundle/` copies |
| Cursor agent-store path `{drive}\Users\<USER>\AppData\Local\Cursor\AgentStores\...` | 2 | `results/posctrl/s2_report.txt`, `bundle/results/posctrl/s2_report.txt` |
| Hard-coded reference PDF path `{drive}\Users\<USER>\Desktop\Prevalent mesenchymal drift ... .pdf` (line 37) | 4 | `src/md2_common.py`, `bundle/src/md2_common.py`, `paper_package/src/md2_common.py`, `bundle/paper_package/src/md2_common.py` |
| Bare account name, `user segment(s): <USER>` (output of an earlier privacy check) | 310 | `bundle/_build/check_c_privacy.csv` |
| Correspondence email containing `<USER>` as a substring | 16 | `README.md`, `FINDINGS_VERIFY.md`, `results/verify/claims.csv` (2), `bundle/README.md`, `bundle/FINDINGS_VERIFY.md`, `bundle/results/verify/claims.csv` (2), `bundle/_build/check_c_privacy.csv` (8) |
| Absolute drive path without the username: `{drive}\nrn\mingw\usr\bin\unzip.exe` (line 14) | 1 | `bundle/_build/zip_check.py` |

A drive-path false positive was checked and dismissed: `notebooks/rerun_stageA_age_axis.ipynb` line 349, where the prose `...not A:` is followed by an escaped `\n`.

Binary files, not edited (see the Binaries section):

| file | git status | content |
|---|---|---|
| `bundle.zip` | untracked | 1848 user-path occurrences in 321 members (stale copies of the pre-scrub `bundle/` files), plus the email |
| 92 `src/__pycache__/*.pyc` and 1 `bundle/_build/__pycache__/*.pyc` | ignored | `co_filename` holds the repo root (92 of them); `md2_common` holds the PDF path |
| `data/processed/{onek1k,tissue_brain_aging,tissue_brain_wm,tissue_muscle,tissue_retina_sc}_pseudobulk.h5ad` | ignored | 1 HDF5 attribute each: the source path under `data/raw/` |
| 8 `.docx` | see below | email only |

## .docx

Not opened for writing. None of them contains the repo root or any `{drive}` path. Each contains the account name only inside the correspondence email in `word/document.xml`:

- `paper/rejuvenation_readouts_preprint_3.docx`
- `paper/rejuvenation_readouts_preprint_9.docx`
- `paper/rejuvenation_readouts_preprint_10.docx`
- `paper/rejuvenation_readouts_preprint_11.docx`
- `paper/rejuvenation_readouts_preprint_12.docx`
- `bundle/paper/rejuvenation_readouts_preprint_3.docx`
- `bundle/paper/rejuvenation_readouts_preprint_10.docx`
- `bundle/paper/rejuvenation_readouts_preprint_11.docx`

## Binaries

A byte-length change breaks these formats, so they were left as they are:

- `.pyc` files are git-ignored and are rebuilt on the next import. Deleting `__pycache__/` removes the hits.
- `.h5ad` files are git-ignored and were regenerated from `data/raw/`. Removing the hits means rewriting one `uns`/attribute string, or rebuilding the files, which counts as rerunning.
- `bundle.zip` is an untracked archive of `bundle/`. Rebuilding it from the scrubbed `bundle/` removes all of its path hits.

## Files changed

Counts are occurrences replaced. `t`/`u`/`i` = tracked / untracked / git-ignored.

641 files: 11 tracked, 630 untracked, 0 ignored.

| file | git | replaced |
|---|---|---:|
| `FINDINGS_FIBRO.md` | u | 1 |
| `FINDINGS_FIBRO2.md` | u | 1 |
| `FINDINGS_FIBRO3.md` | u | 10 |
| `FINDINGS_MD.md` | u | 1 |
| `FINDINGS_MD2.md` | u | 31 |
| `FINDINGS_MD3.md` | u | 34 |
| `FINDINGS_MD4.md` | u | 26 |
| `FINDINGS_MD5.md` | u | 28 |
| `FINDINGS_PLANE.md` | u | 30 |
| `FINDINGS_SAME.md` | u | 1 |
| `FINDINGS_SENG.md` | u | 2 |
| `FINDINGS_SOUTH.md` | u | 1 |
| `FINDINGS_SOUTH2.md` | u | 1 |
| `FINDINGS_SOUTH3.md` | u | 1 |
| `FINDINGS_TOWARD.md` | u | 2 |
| `PROGRESS_FIBRO.md` | u | 3 |
| `PROGRESS_FIBRO2.md` | u | 3 |
| `PROGRESS_FIBRO3.md` | u | 15 |
| `PROGRESS_MD.md` | u | 32 |
| `PROGRESS_MD2.md` | u | 39 |
| `PROGRESS_MD3.md` | u | 37 |
| `PROGRESS_MD4.md` | u | 31 |
| `PROGRESS_MD5.md` | u | 48 |
| `PROGRESS_PLANE.md` | u | 36 |
| `bundle/FINDINGS_FIBRO.md` | u | 1 |
| `bundle/FINDINGS_FIBRO2.md` | u | 1 |
| `bundle/FINDINGS_FIBRO3.md` | u | 10 |
| `bundle/FINDINGS_MD.md` | u | 1 |
| `bundle/FINDINGS_MD2.md` | u | 31 |
| `bundle/FINDINGS_MD3.md` | u | 34 |
| `bundle/FINDINGS_MD4.md` | u | 26 |
| `bundle/FINDINGS_MD5.md` | u | 28 |
| `bundle/FINDINGS_PLANE.md` | u | 30 |
| `bundle/FINDINGS_SAME.md` | u | 1 |
| `bundle/FINDINGS_SENG.md` | u | 2 |
| `bundle/FINDINGS_SOUTH.md` | u | 1 |
| `bundle/FINDINGS_SOUTH2.md` | u | 1 |
| `bundle/FINDINGS_SOUTH3.md` | u | 1 |
| `bundle/FINDINGS_TOWARD.md` | u | 2 |
| `bundle/MANIFEST.md` | u | 1 |
| `bundle/PROGRESS_FIBRO.md` | u | 3 |
| `bundle/PROGRESS_FIBRO2.md` | u | 3 |
| `bundle/PROGRESS_FIBRO3.md` | u | 15 |
| `bundle/PROGRESS_MD.md` | u | 32 |
| `bundle/PROGRESS_MD2.md` | u | 39 |
| `bundle/PROGRESS_MD3.md` | u | 37 |
| `bundle/PROGRESS_MD4.md` | u | 31 |
| `bundle/PROGRESS_MD5.md` | u | 48 |
| `bundle/PROGRESS_PLANE.md` | u | 36 |
| `bundle/paper_package/FINDINGS_FIBRO.md` | u | 1 |
| `bundle/paper_package/FINDINGS_FIBRO2.md` | u | 1 |
| `bundle/paper_package/FINDINGS_FIBRO3.md` | u | 10 |
| `bundle/paper_package/FINDINGS_MD.md` | u | 1 |
| `bundle/paper_package/FINDINGS_MD2.md` | u | 31 |
| `bundle/paper_package/FINDINGS_MD3.md` | u | 34 |
| `bundle/paper_package/FINDINGS_MD4.md` | u | 26 |
| `bundle/paper_package/FINDINGS_MD5.md` | u | 28 |
| `bundle/paper_package/FINDINGS_PLANE.md` | u | 30 |
| `bundle/paper_package/FINDINGS_SAME.md` | u | 1 |
| `bundle/paper_package/FINDINGS_SENG.md` | u | 2 |
| `bundle/paper_package/FINDINGS_TOWARD.md` | u | 2 |
| `bundle/paper_package/PROGRESS_FIBRO.md` | u | 3 |
| `bundle/paper_package/PROGRESS_FIBRO2.md` | u | 3 |
| `bundle/paper_package/PROGRESS_FIBRO3.md` | u | 15 |
| `bundle/paper_package/PROGRESS_MD.md` | u | 32 |
| `bundle/paper_package/PROGRESS_MD2.md` | u | 39 |
| `bundle/paper_package/PROGRESS_MD3.md` | u | 37 |
| `bundle/paper_package/PROGRESS_MD4.md` | u | 31 |
| `bundle/paper_package/PROGRESS_MD5.md` | u | 48 |
| `bundle/paper_package/PROGRESS_PLANE.md` | u | 36 |
| `bundle/paper_package/results/fibro/freeze.json` | u | 1 |
| `bundle/paper_package/results/fibro/stage1_summary.json` | u | 1 |
| `bundle/paper_package/results/md3/t2_instruments.csv` | u | 4 |
| `bundle/paper_package/results/toward/anchors.json` | u | 1 |
| `bundle/results/brain/r0_peek.txt` | u | 1 |
| `bundle/results/brain/r1_r3_analyze.txt` | u | 2 |
| `bundle/results/brain_phase1/p0_obs_audit.json` | u | 1 |
| `bundle/results/brain_phase1/p0_report.txt` | u | 1 |
| `bundle/results/brain_phase1/p1_report.txt` | u | 4 |
| `bundle/results/brain_phase1/p2_report.txt` | u | 1 |
| `bundle/results/brain_phase1/p3_report.txt` | u | 1 |
| `bundle/results/brain_phase1/p4_report.txt` | u | 1 |
| `bundle/results/brain_phase1/p5_report.txt` | u | 1 |
| `bundle/results/external/freeze.json` | u | 1 |
| `bundle/results/external/freeze_report.txt` | u | 1 |
| `bundle/results/external/manifest.json` | u | 2 |
| `bundle/results/external/project_report.txt` | u | 2 |
| `bundle/results/external/stage0_report.txt` | u | 2 |
| `bundle/results/external/stage0_summary.json` | u | 1 |
| `bundle/results/fibro/freeze.json` | u | 1 |
| `bundle/results/fibro/manifest.json` | u | 5 |
| `bundle/results/fibro/stage1_gene_overlap.json` | u | 1 |
| `bundle/results/fibro/stage1_report.txt` | u | 1 |
| `bundle/results/fibro/stage1_summary.json` | u | 1 |
| `bundle/results/fibro/stage2_STOP.json` | u | 2 |
| `bundle/results/fibro/stage2_downloads.json` | u | 1 |
| `bundle/results/fibro/stage2_h5_files.json` | u | 8 |
| `bundle/results/fibro/stage2_report.txt` | u | 3 |
| `bundle/results/fibro/stage2_summary.json` | u | 1 |
| `bundle/results/fibro2/manifest.json` | u | 9 |
| `bundle/results/fibro2/run_report.txt` | u | 23 |
| `bundle/results/fibro2/ta_GSE113957_downloads.json` | u | 2 |
| `bundle/results/fibro2/ta_GSE226189_downloads.json` | u | 2 |
| `bundle/results/fibro2/ta_GSE307377_downloads.json` | u | 4 |
| `bundle/results/fibro2/ta_results.csv` | u | 1 |
| `bundle/results/fibro2/tb_summary.json` | u | 1 |
| `bundle/results/fibro3/manifest.json` | u | 10 |
| `bundle/results/fibro3/run_report.txt` | u | 29 |
| `bundle/results/fibro3/t1_summary.json` | u | 1 |
| `bundle/results/fibro3/t4_result.csv` | u | 3 |
| `bundle/results/fibro3/t4_result.json` | u | 3 |
| `bundle/results/fibro3/t4_summary.json` | u | 3 |
| `bundle/results/genespace/run_report.txt` | u | 5 |
| `bundle/results/genespace_diag/diag_log.txt` | u | 1 |
| `bundle/results/genespace_diag/loading_check.json` | u | 1 |
| `bundle/results/gtex/manifest.json` | u | 4 |
| `bundle/results/gtex/stage0_report.txt` | u | 2 |
| `bundle/results/gtex/stage1_report.txt` | u | 1 |
| `bundle/results/lowdim/l1_report.txt` | u | 1 |
| `bundle/results/md/manifest.json` | u | 2 |
| `bundle/results/md/run_report.txt` | u | 4 |
| `bundle/results/md/s1_gate_reval.txt` | u | 1 |
| `bundle/results/md/sources.csv` | u | 35 |
| `bundle/results/md/sources.json` | u | 35 |
| `bundle/results/md2/cluster_meta.json` | u | 1 |
| `bundle/results/md2/genesets_summary.json` | u | 2 |
| `bundle/results/md2/manifest.json` | u | 30 |
| `bundle/results/md2/run_report.txt` | u | 13 |
| `bundle/results/md2/t3_instruments.csv` | u | 3 |
| `bundle/results/md2/t3_report.txt` | u | 5 |
| `bundle/results/md3/genesets.json` | u | 1 |
| `bundle/results/md3/genesets_summary.json` | u | 1 |
| `bundle/results/md3/manifest.json` | u | 33 |
| `bundle/results/md3/run_report.txt` | u | 19 |
| `bundle/results/md3/t1_ams_meta_GM00731.json` | u | 1 |
| `bundle/results/md3/t1_ams_meta_GM23815.json` | u | 1 |
| `bundle/results/md3/t2_gse113957_repair.txt` | u | 3 |
| `bundle/results/md3/t2_instruments.csv` | u | 4 |
| `bundle/results/md4/genesets_summary.json` | u | 2 |
| `bundle/results/md4/manifest.json` | u | 25 |
| `bundle/results/md4/run_report.txt` | u | 8 |
| `bundle/results/md5/genesets_summary.json` | u | 2 |
| `bundle/results/md5/manifest.json` | u | 27 |
| `bundle/results/md5/run_report.txt` | u | 8 |
| `bundle/results/md5/t2_extrap.csv` | u | 17 |
| `bundle/results/nb_phase0_exec.err` | u | 1 |
| `bundle/results/nb_phase0_exec.log` | u | 8 |
| `bundle/results/nb_phase1_exec.err` | u | 1 |
| `bundle/results/newstory/run_log.txt` | u | 15 |
| `bundle/results/paper_figs/environment_versions.txt` | u | 1 |
| `bundle/results/phase0_GSE217460_stream_convert.log` | u | 1 |
| `bundle/results/phase1_pseudobulk_build.log` | u | 1 |
| `bundle/results/plane/genesets_summary.json` | u | 5 |
| `bundle/results/plane/manifest.json` | u | 27 |
| `bundle/results/plane/run_report.txt` | u | 11 |
| `bundle/results/plane/t3_figures.json` | u | 2 |
| `bundle/results/plane/t3_report.txt` | u | 3 |
| `bundle/results/posctrl/c1_dense_shared_r015.json` | u | 1 |
| `bundle/results/posctrl/c1_dense_shared_r030.json` | u | 1 |
| `bundle/results/posctrl/c1_dense_shared_r050.json` | u | 1 |
| `bundle/results/posctrl/c1_dense_typespec_r030.json` | u | 1 |
| `bundle/results/posctrl/c1_report.txt` | u | 7 |
| `bundle/results/posctrl/c1_sparse_shared_r030.json` | u | 1 |
| `bundle/results/posctrl/c1_sparse_typespec_r030.json` | u | 1 |
| `bundle/results/posctrl/c2_report.txt` | u | 2 |
| `bundle/results/posctrl/c2a.json` | u | 1 |
| `bundle/results/posctrl/c2b.json` | u | 1 |
| `bundle/results/posctrl/c3_draws.csv` | u | 62 |
| `bundle/results/posctrl/c3_planted_n060_d00.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n060_d01.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n060_d02.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n060_d03.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n060_d04.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n060_d05.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n060_d06.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n060_d07.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n060_d08.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n060_d09.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n100_d00.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n100_d01.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n100_d02.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n100_d03.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n100_d04.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n100_d05.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n100_d06.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n100_d07.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n100_d08.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n100_d09.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n150_d00.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n150_d01.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n150_d02.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n150_d03.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n150_d04.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n150_d05.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n150_d06.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n150_d07.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n150_d08.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n150_d09.json` | u | 1 |
| `bundle/results/posctrl/c3_planted_n233_d00.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n060_d00.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n060_d01.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n060_d02.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n060_d03.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n060_d04.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n060_d05.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n060_d06.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n060_d07.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n060_d08.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n060_d09.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n100_d00.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n100_d01.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n100_d02.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n100_d03.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n100_d04.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n100_d05.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n100_d06.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n100_d07.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n100_d08.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n100_d09.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n150_d00.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n150_d01.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n150_d02.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n150_d03.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n150_d04.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n150_d05.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n150_d06.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n150_d07.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n150_d08.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n150_d09.json` | u | 1 |
| `bundle/results/posctrl/c3_real_n233_d00.json` | u | 1 |
| `bundle/results/posctrl/s2_c1_dense_shared_r015.json` | u | 2 |
| `bundle/results/posctrl/s2_c1_dense_shared_r030.json` | u | 2 |
| `bundle/results/posctrl/s2_c1_dense_shared_r050.json` | u | 2 |
| `bundle/results/posctrl/s2_c1_dense_typespec_r030.json` | u | 2 |
| `bundle/results/posctrl/s2_c1_sparse_shared_r030.json` | u | 2 |
| `bundle/results/posctrl/s2_c1_sparse_typespec_r030.json` | u | 2 |
| `bundle/results/posctrl/s2_c2a.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n060_d00.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n060_d01.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n060_d02.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n060_d03.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n060_d04.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n060_d05.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n060_d06.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n060_d07.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n060_d08.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n060_d09.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n100_d00.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n100_d01.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n100_d02.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n100_d03.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n100_d04.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n100_d05.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n100_d06.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n100_d07.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n100_d08.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n100_d09.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n150_d00.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n150_d01.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n150_d02.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n150_d03.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n150_d04.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n150_d05.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n150_d06.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n150_d07.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n150_d08.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n150_d09.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_planted_n233_d00.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n060_d00.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n060_d01.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n060_d02.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n060_d03.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n060_d04.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n060_d05.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n060_d06.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n060_d07.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n060_d08.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n060_d09.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n100_d00.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n100_d01.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n100_d02.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n100_d03.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n100_d04.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n100_d05.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n100_d06.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n100_d07.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n100_d08.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n100_d09.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n150_d00.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n150_d01.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n150_d02.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n150_d03.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n150_d04.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n150_d05.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n150_d06.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n150_d07.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n150_d08.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n150_d09.json` | u | 2 |
| `bundle/results/posctrl/s2_c3_real_n233_d00.json` | u | 2 |
| `bundle/results/posctrl/s2_gene_ridge.json` | u | 2 |
| `bundle/results/posctrl/s2_report.txt` | u | 70 |
| `bundle/results/posctrl/s3_c1_dense_shared_r030.json` | u | 2 |
| `bundle/results/posctrl/s3_c2a.json` | u | 2 |
| `bundle/results/posctrl/s3_gene_ridge.json` | u | 2 |
| `bundle/results/posctrl/s3_idresid.json` | u | 2 |
| `bundle/results/posctrl/s3_report.txt` | u | 4 |
| `bundle/results/posctrl/s3_summary.json` | u | 8 |
| `bundle/results/relabel_check/run_log.txt` | u | 3 |
| `bundle/results/relabel_check/stats_console.txt` | u | 1 |
| `bundle/results/rerun/nb_stageA_exec.log` | u | 1 |
| `bundle/results/rerun/nb_stageB_exec.log` | u | 1 |
| `bundle/results/same/run_report.txt` | u | 6 |
| `bundle/results/seng/download.json` | u | 1 |
| `bundle/results/seng/inventory.json` | u | 6 |
| `bundle/results/seng/run_report.txt` | u | 3 |
| `bundle/results/seng_why/run_report.txt` | u | 1 |
| `bundle/results/south/inventory.json` | u | 1 |
| `bundle/results/south/manifest.json` | u | 1 |
| `bundle/results/south/south_run.log` | u | 4 |
| `bundle/results/south2/console.txt` | u | 12 |
| `bundle/results/south2/inventory.json` | u | 3 |
| `bundle/results/south2/manifest.json` | u | 1 |
| `bundle/results/south2/south2_run.log` | u | 12 |
| `bundle/results/south3/console.txt` | u | 13 |
| `bundle/results/south3/inventory.json` | u | 3 |
| `bundle/results/south3/manifest.json` | u | 1 |
| `bundle/results/south3/south3_run.log` | u | 13 |
| `bundle/results/south5/south5_run.log` | u | 4 |
| `bundle/results/south6/full_run_stdout.txt` | u | 4 |
| `bundle/results/south6/inventory.json` | u | 1 |
| `bundle/results/south6/smoke_stdout.txt` | u | 2 |
| `bundle/results/south6/south6_run.log` | u | 6 |
| `bundle/results/tissue/brain_wm_analyze.txt` | u | 2 |
| `bundle/results/tissue/brain_wm_summary.json` | u | 1 |
| `bundle/results/tissue/muscle_analyze.txt` | u | 2 |
| `bundle/results/tissue/muscle_summary.json` | u | 1 |
| `bundle/results/tissue/retina_sc_analyze.txt` | u | 2 |
| `bundle/results/tissue/retina_sc_summary.json` | u | 1 |
| `bundle/results/toward/anchors.json` | u | 1 |
| `bundle/results/toward/run_report.txt` | u | 7 |
| `bundle/results/trajectory/g3_orthologs.txt` | u | 2 |
| `bundle/results/trajectory/g3_orthologs_summary.json` | u | 1 |
| `bundle/results/trajectory/g3_report.txt` | u | 2 |
| `notebooks/phase0_data_acquisition.ipynb` | t | 9 |
| `notebooks/phase1_scores.ipynb` | t | 2 |
| `notebooks/rerun_stageA_age_axis.ipynb` | t | 2 |
| `notebooks/rerun_stageB_composition.ipynb` | t | 1 |
| `paper_package/FINDINGS_FIBRO.md` | u | 1 |
| `paper_package/FINDINGS_FIBRO2.md` | u | 1 |
| `paper_package/FINDINGS_FIBRO3.md` | u | 10 |
| `paper_package/FINDINGS_MD.md` | u | 1 |
| `paper_package/FINDINGS_MD2.md` | u | 31 |
| `paper_package/FINDINGS_MD3.md` | u | 34 |
| `paper_package/FINDINGS_MD4.md` | u | 26 |
| `paper_package/FINDINGS_MD5.md` | u | 28 |
| `paper_package/FINDINGS_PLANE.md` | u | 30 |
| `paper_package/FINDINGS_SAME.md` | u | 1 |
| `paper_package/FINDINGS_SENG.md` | u | 2 |
| `paper_package/FINDINGS_TOWARD.md` | u | 2 |
| `paper_package/PROGRESS_FIBRO.md` | u | 3 |
| `paper_package/PROGRESS_FIBRO2.md` | u | 3 |
| `paper_package/PROGRESS_FIBRO3.md` | u | 15 |
| `paper_package/PROGRESS_MD.md` | u | 32 |
| `paper_package/PROGRESS_MD2.md` | u | 39 |
| `paper_package/PROGRESS_MD3.md` | u | 37 |
| `paper_package/PROGRESS_MD4.md` | u | 31 |
| `paper_package/PROGRESS_MD5.md` | u | 48 |
| `paper_package/PROGRESS_PLANE.md` | u | 36 |
| `paper_package/results/fibro/freeze.json` | u | 1 |
| `paper_package/results/fibro/stage1_summary.json` | u | 1 |
| `paper_package/results/md3/t2_instruments.csv` | u | 4 |
| `paper_package/results/toward/anchors.json` | u | 1 |
| `results/brain/r0_peek.txt` | u | 1 |
| `results/brain/r1_r3_analyze.txt` | u | 2 |
| `results/brain_phase1/p0_obs_audit.json` | u | 1 |
| `results/brain_phase1/p0_report.txt` | u | 1 |
| `results/brain_phase1/p1_report.txt` | u | 4 |
| `results/brain_phase1/p2_report.txt` | u | 1 |
| `results/brain_phase1/p3_report.txt` | u | 1 |
| `results/brain_phase1/p4_report.txt` | u | 1 |
| `results/brain_phase1/p5_report.txt` | u | 1 |
| `results/external/freeze.json` | u | 1 |
| `results/external/freeze_report.txt` | u | 1 |
| `results/external/manifest.json` | u | 2 |
| `results/external/project_report.txt` | u | 2 |
| `results/external/stage0_report.txt` | u | 2 |
| `results/external/stage0_summary.json` | u | 1 |
| `results/fibro/freeze.json` | u | 1 |
| `results/fibro/manifest.json` | u | 5 |
| `results/fibro/stage1_gene_overlap.json` | u | 1 |
| `results/fibro/stage1_report.txt` | u | 1 |
| `results/fibro/stage1_summary.json` | u | 1 |
| `results/fibro/stage2_STOP.json` | u | 2 |
| `results/fibro/stage2_downloads.json` | u | 1 |
| `results/fibro/stage2_h5_files.json` | u | 8 |
| `results/fibro/stage2_report.txt` | u | 3 |
| `results/fibro/stage2_summary.json` | u | 1 |
| `results/fibro2/manifest.json` | u | 9 |
| `results/fibro2/run_report.txt` | u | 23 |
| `results/fibro2/ta_GSE113957_downloads.json` | u | 2 |
| `results/fibro2/ta_GSE226189_downloads.json` | u | 2 |
| `results/fibro2/ta_GSE307377_downloads.json` | u | 4 |
| `results/fibro2/ta_results.csv` | u | 1 |
| `results/fibro2/tb_summary.json` | u | 1 |
| `results/fibro3/manifest.json` | u | 10 |
| `results/fibro3/run_report.txt` | u | 29 |
| `results/fibro3/t1_summary.json` | u | 1 |
| `results/fibro3/t4_result.csv` | u | 3 |
| `results/fibro3/t4_result.json` | u | 3 |
| `results/fibro3/t4_summary.json` | u | 3 |
| `results/genespace/run_report.txt` | u | 5 |
| `results/genespace_diag/diag_log.txt` | u | 1 |
| `results/genespace_diag/loading_check.json` | u | 1 |
| `results/gtex/manifest.json` | u | 4 |
| `results/gtex/stage0_report.txt` | u | 2 |
| `results/gtex/stage1_report.txt` | u | 1 |
| `results/lowdim/l1_report.txt` | u | 1 |
| `results/md/manifest.json` | u | 2 |
| `results/md/run_report.txt` | u | 4 |
| `results/md/s1_gate_reval.txt` | u | 1 |
| `results/md/sources.csv` | u | 35 |
| `results/md/sources.json` | u | 35 |
| `results/md2/cluster_meta.json` | u | 1 |
| `results/md2/genesets_summary.json` | u | 2 |
| `results/md2/manifest.json` | u | 30 |
| `results/md2/run_report.txt` | u | 13 |
| `results/md2/t3_instruments.csv` | u | 3 |
| `results/md2/t3_report.txt` | u | 5 |
| `results/md3/genesets.json` | u | 1 |
| `results/md3/genesets_summary.json` | u | 1 |
| `results/md3/manifest.json` | u | 33 |
| `results/md3/run_report.txt` | u | 19 |
| `results/md3/t1_ams_meta_GM00731.json` | u | 1 |
| `results/md3/t1_ams_meta_GM23815.json` | u | 1 |
| `results/md3/t2_gse113957_repair.txt` | u | 3 |
| `results/md3/t2_instruments.csv` | u | 4 |
| `results/md4/genesets_summary.json` | u | 2 |
| `results/md4/manifest.json` | u | 25 |
| `results/md4/run_report.txt` | u | 8 |
| `results/md5/genesets_summary.json` | u | 2 |
| `results/md5/manifest.json` | u | 27 |
| `results/md5/run_report.txt` | u | 8 |
| `results/md5/t2_extrap.csv` | u | 17 |
| `results/nb_phase0_exec.err` | t | 1 |
| `results/nb_phase0_exec.log` | t | 8 |
| `results/nb_phase1_exec.err` | t | 1 |
| `results/newstory/run_log.txt` | u | 15 |
| `results/paper_figs/environment_versions.txt` | u | 1 |
| `results/phase0_GSE217460_stream_convert.log` | t | 1 |
| `results/phase1_pseudobulk_build.log` | t | 1 |
| `results/plane/genesets_summary.json` | u | 5 |
| `results/plane/manifest.json` | u | 27 |
| `results/plane/run_report.txt` | u | 11 |
| `results/plane/t3_figures.json` | u | 2 |
| `results/plane/t3_report.txt` | u | 3 |
| `results/posctrl/c1_dense_shared_r015.json` | u | 1 |
| `results/posctrl/c1_dense_shared_r030.json` | u | 1 |
| `results/posctrl/c1_dense_shared_r050.json` | u | 1 |
| `results/posctrl/c1_dense_typespec_r030.json` | u | 1 |
| `results/posctrl/c1_report.txt` | u | 7 |
| `results/posctrl/c1_sparse_shared_r030.json` | u | 1 |
| `results/posctrl/c1_sparse_typespec_r030.json` | u | 1 |
| `results/posctrl/c2_report.txt` | u | 2 |
| `results/posctrl/c2a.json` | u | 1 |
| `results/posctrl/c2b.json` | u | 1 |
| `results/posctrl/c3_draws.csv` | u | 62 |
| `results/posctrl/c3_planted_n060_d00.json` | u | 1 |
| `results/posctrl/c3_planted_n060_d01.json` | u | 1 |
| `results/posctrl/c3_planted_n060_d02.json` | u | 1 |
| `results/posctrl/c3_planted_n060_d03.json` | u | 1 |
| `results/posctrl/c3_planted_n060_d04.json` | u | 1 |
| `results/posctrl/c3_planted_n060_d05.json` | u | 1 |
| `results/posctrl/c3_planted_n060_d06.json` | u | 1 |
| `results/posctrl/c3_planted_n060_d07.json` | u | 1 |
| `results/posctrl/c3_planted_n060_d08.json` | u | 1 |
| `results/posctrl/c3_planted_n060_d09.json` | u | 1 |
| `results/posctrl/c3_planted_n100_d00.json` | u | 1 |
| `results/posctrl/c3_planted_n100_d01.json` | u | 1 |
| `results/posctrl/c3_planted_n100_d02.json` | u | 1 |
| `results/posctrl/c3_planted_n100_d03.json` | u | 1 |
| `results/posctrl/c3_planted_n100_d04.json` | u | 1 |
| `results/posctrl/c3_planted_n100_d05.json` | u | 1 |
| `results/posctrl/c3_planted_n100_d06.json` | u | 1 |
| `results/posctrl/c3_planted_n100_d07.json` | u | 1 |
| `results/posctrl/c3_planted_n100_d08.json` | u | 1 |
| `results/posctrl/c3_planted_n100_d09.json` | u | 1 |
| `results/posctrl/c3_planted_n150_d00.json` | u | 1 |
| `results/posctrl/c3_planted_n150_d01.json` | u | 1 |
| `results/posctrl/c3_planted_n150_d02.json` | u | 1 |
| `results/posctrl/c3_planted_n150_d03.json` | u | 1 |
| `results/posctrl/c3_planted_n150_d04.json` | u | 1 |
| `results/posctrl/c3_planted_n150_d05.json` | u | 1 |
| `results/posctrl/c3_planted_n150_d06.json` | u | 1 |
| `results/posctrl/c3_planted_n150_d07.json` | u | 1 |
| `results/posctrl/c3_planted_n150_d08.json` | u | 1 |
| `results/posctrl/c3_planted_n150_d09.json` | u | 1 |
| `results/posctrl/c3_planted_n233_d00.json` | u | 1 |
| `results/posctrl/c3_real_n060_d00.json` | u | 1 |
| `results/posctrl/c3_real_n060_d01.json` | u | 1 |
| `results/posctrl/c3_real_n060_d02.json` | u | 1 |
| `results/posctrl/c3_real_n060_d03.json` | u | 1 |
| `results/posctrl/c3_real_n060_d04.json` | u | 1 |
| `results/posctrl/c3_real_n060_d05.json` | u | 1 |
| `results/posctrl/c3_real_n060_d06.json` | u | 1 |
| `results/posctrl/c3_real_n060_d07.json` | u | 1 |
| `results/posctrl/c3_real_n060_d08.json` | u | 1 |
| `results/posctrl/c3_real_n060_d09.json` | u | 1 |
| `results/posctrl/c3_real_n100_d00.json` | u | 1 |
| `results/posctrl/c3_real_n100_d01.json` | u | 1 |
| `results/posctrl/c3_real_n100_d02.json` | u | 1 |
| `results/posctrl/c3_real_n100_d03.json` | u | 1 |
| `results/posctrl/c3_real_n100_d04.json` | u | 1 |
| `results/posctrl/c3_real_n100_d05.json` | u | 1 |
| `results/posctrl/c3_real_n100_d06.json` | u | 1 |
| `results/posctrl/c3_real_n100_d07.json` | u | 1 |
| `results/posctrl/c3_real_n100_d08.json` | u | 1 |
| `results/posctrl/c3_real_n100_d09.json` | u | 1 |
| `results/posctrl/c3_real_n150_d00.json` | u | 1 |
| `results/posctrl/c3_real_n150_d01.json` | u | 1 |
| `results/posctrl/c3_real_n150_d02.json` | u | 1 |
| `results/posctrl/c3_real_n150_d03.json` | u | 1 |
| `results/posctrl/c3_real_n150_d04.json` | u | 1 |
| `results/posctrl/c3_real_n150_d05.json` | u | 1 |
| `results/posctrl/c3_real_n150_d06.json` | u | 1 |
| `results/posctrl/c3_real_n150_d07.json` | u | 1 |
| `results/posctrl/c3_real_n150_d08.json` | u | 1 |
| `results/posctrl/c3_real_n150_d09.json` | u | 1 |
| `results/posctrl/c3_real_n233_d00.json` | u | 1 |
| `results/posctrl/s2_c1_dense_shared_r015.json` | u | 2 |
| `results/posctrl/s2_c1_dense_shared_r030.json` | u | 2 |
| `results/posctrl/s2_c1_dense_shared_r050.json` | u | 2 |
| `results/posctrl/s2_c1_dense_typespec_r030.json` | u | 2 |
| `results/posctrl/s2_c1_sparse_shared_r030.json` | u | 2 |
| `results/posctrl/s2_c1_sparse_typespec_r030.json` | u | 2 |
| `results/posctrl/s2_c2a.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n060_d00.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n060_d01.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n060_d02.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n060_d03.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n060_d04.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n060_d05.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n060_d06.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n060_d07.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n060_d08.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n060_d09.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n100_d00.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n100_d01.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n100_d02.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n100_d03.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n100_d04.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n100_d05.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n100_d06.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n100_d07.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n100_d08.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n100_d09.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n150_d00.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n150_d01.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n150_d02.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n150_d03.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n150_d04.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n150_d05.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n150_d06.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n150_d07.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n150_d08.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n150_d09.json` | u | 2 |
| `results/posctrl/s2_c3_planted_n233_d00.json` | u | 2 |
| `results/posctrl/s2_c3_real_n060_d00.json` | u | 2 |
| `results/posctrl/s2_c3_real_n060_d01.json` | u | 2 |
| `results/posctrl/s2_c3_real_n060_d02.json` | u | 2 |
| `results/posctrl/s2_c3_real_n060_d03.json` | u | 2 |
| `results/posctrl/s2_c3_real_n060_d04.json` | u | 2 |
| `results/posctrl/s2_c3_real_n060_d05.json` | u | 2 |
| `results/posctrl/s2_c3_real_n060_d06.json` | u | 2 |
| `results/posctrl/s2_c3_real_n060_d07.json` | u | 2 |
| `results/posctrl/s2_c3_real_n060_d08.json` | u | 2 |
| `results/posctrl/s2_c3_real_n060_d09.json` | u | 2 |
| `results/posctrl/s2_c3_real_n100_d00.json` | u | 2 |
| `results/posctrl/s2_c3_real_n100_d01.json` | u | 2 |
| `results/posctrl/s2_c3_real_n100_d02.json` | u | 2 |
| `results/posctrl/s2_c3_real_n100_d03.json` | u | 2 |
| `results/posctrl/s2_c3_real_n100_d04.json` | u | 2 |
| `results/posctrl/s2_c3_real_n100_d05.json` | u | 2 |
| `results/posctrl/s2_c3_real_n100_d06.json` | u | 2 |
| `results/posctrl/s2_c3_real_n100_d07.json` | u | 2 |
| `results/posctrl/s2_c3_real_n100_d08.json` | u | 2 |
| `results/posctrl/s2_c3_real_n100_d09.json` | u | 2 |
| `results/posctrl/s2_c3_real_n150_d00.json` | u | 2 |
| `results/posctrl/s2_c3_real_n150_d01.json` | u | 2 |
| `results/posctrl/s2_c3_real_n150_d02.json` | u | 2 |
| `results/posctrl/s2_c3_real_n150_d03.json` | u | 2 |
| `results/posctrl/s2_c3_real_n150_d04.json` | u | 2 |
| `results/posctrl/s2_c3_real_n150_d05.json` | u | 2 |
| `results/posctrl/s2_c3_real_n150_d06.json` | u | 2 |
| `results/posctrl/s2_c3_real_n150_d07.json` | u | 2 |
| `results/posctrl/s2_c3_real_n150_d08.json` | u | 2 |
| `results/posctrl/s2_c3_real_n150_d09.json` | u | 2 |
| `results/posctrl/s2_c3_real_n233_d00.json` | u | 2 |
| `results/posctrl/s2_gene_ridge.json` | u | 2 |
| `results/posctrl/s2_report.txt` | u | 70 |
| `results/posctrl/s3_c1_dense_shared_r030.json` | u | 2 |
| `results/posctrl/s3_c2a.json` | u | 2 |
| `results/posctrl/s3_gene_ridge.json` | u | 2 |
| `results/posctrl/s3_idresid.json` | u | 2 |
| `results/posctrl/s3_report.txt` | u | 4 |
| `results/posctrl/s3_summary.json` | u | 8 |
| `results/relabel_check/run_log.txt` | u | 3 |
| `results/relabel_check/stats_console.txt` | u | 1 |
| `results/rerun/nb_stageA_exec.log` | t | 1 |
| `results/rerun/nb_stageB_exec.log` | t | 1 |
| `results/same/run_report.txt` | u | 6 |
| `results/seng/download.json` | u | 1 |
| `results/seng/inventory.json` | u | 6 |
| `results/seng/run_report.txt` | u | 3 |
| `results/seng_why/run_report.txt` | u | 1 |
| `results/south/inventory.json` | u | 1 |
| `results/south/manifest.json` | u | 1 |
| `results/south/south_run.log` | u | 4 |
| `results/south2/console.txt` | u | 12 |
| `results/south2/inventory.json` | u | 3 |
| `results/south2/manifest.json` | u | 1 |
| `results/south2/south2_run.log` | u | 12 |
| `results/south3/console.txt` | u | 13 |
| `results/south3/inventory.json` | u | 3 |
| `results/south3/manifest.json` | u | 1 |
| `results/south3/south3_run.log` | u | 13 |
| `results/south5/south5_run.log` | u | 4 |
| `results/south6/full_run_stdout.txt` | u | 4 |
| `results/south6/inventory.json` | u | 1 |
| `results/south6/smoke_stdout.txt` | u | 2 |
| `results/south6/south6_run.log` | u | 6 |
| `results/tissue/brain_wm_analyze.txt` | u | 2 |
| `results/tissue/brain_wm_summary.json` | u | 1 |
| `results/tissue/muscle_analyze.txt` | u | 2 |
| `results/tissue/muscle_summary.json` | u | 1 |
| `results/tissue/retina_sc_analyze.txt` | u | 2 |
| `results/tissue/retina_sc_summary.json` | u | 1 |
| `results/toward/anchors.json` | u | 1 |
| `results/toward/run_report.txt` | u | 7 |
| `results/trajectory/g3_orthologs.txt` | u | 2 |
| `results/trajectory/g3_orthologs_summary.json` | u | 1 |
| `results/trajectory/g3_report.txt` | u | 2 |

## Public-push preparation (2026-10-02)

Second pass, before the first public push. Text edits and `.gitignore` only. No analysis was rerun, no file was regenerated, nothing was committed. Same notation as above: `<USER>` is the account name, `{drive}` is the drive.

### 1. `.gitignore`

Appended under the comment `# review packages, not public`:

```
bundle/
bundle.zip
paper_package/
paper/rejuvenation_readouts_preprint_3.docx
paper/rejuvenation_readouts_preprint_9.docx
paper/rejuvenation_readouts_preprint_10.docx
paper/rejuvenation_readouts_preprint_11.docx
paper/orphan_images/
```

None of these paths was tracked, so nothing needs `git rm --cached`. `git check-ignore -v` confirms all of them are now ignored. The file keeps LF line endings.

Already ignored before this pass: `__pycache__/` (line 5), `*.pyc` (line 6), `data/processed/` (line 16).

Still pushed from `paper/`: `paper/PREPRINT_DRAFT.md` and `paper/rejuvenation_readouts_preprint_12.docx`. The `.docx` contains the correspondence email only (see the `.docx` section above).

### 2. Home prefix in UTF-16 logs and `s2_report.txt`

The home prefix `{drive}\Users\<USER>` (uppercase drive, backslashes) was replaced by the literal token `<home>`, in the file's own encoding. Exactly one occurrence per file. BOM, line-ending counts and every other byte are unchanged; each file shrank by exactly the length difference of the two strings.

| file | git | encoding | line endings | what followed the prefix | bytes before → after |
|---|---|---|---|---|---|
| `results/nb_phase0_exec.err` | t | UTF-16LE, BOM `FF FE` | 16 CRLF | `\AppData\Local\Temp\ps-script-*.ps1` | 2510 → 2494 |
| `results/nb_phase1_exec.err` | t | UTF-16LE, BOM | 16 CRLF | same | 2468 → 2452 |
| `results/rerun/nb_stageA_exec.log` | t | UTF-16LE, BOM | 16 CRLF | same | 2496 → 2480 |
| `results/rerun/nb_stageB_exec.log` | t | UTF-16LE, BOM | 16 CRLF | same | 2506 → 2490 |
| `bundle/results/nb_phase0_exec.err` | i | UTF-16LE, BOM | 16 CRLF | same | 2510 → 2494 |
| `bundle/results/nb_phase1_exec.err` | i | UTF-16LE, BOM | 16 CRLF | same | 2468 → 2452 |
| `bundle/results/rerun/nb_stageA_exec.log` | i | UTF-16LE, BOM | 16 CRLF | same | 2496 → 2480 |
| `bundle/results/rerun/nb_stageB_exec.log` | i | UTF-16LE, BOM | 16 CRLF | same | 2506 → 2490 |
| `results/posctrl/s2_report.txt` | u | ASCII, no BOM | 1335 CRLF | `\AppData\Local\Cursor\AgentStores\...` | 68283 → 68275 |
| `bundle/results/posctrl/s2_report.txt` | i | ASCII, no BOM | 1335 CRLF | same | 68283 → 68275 |

`t`/`u`/`i` = tracked / untracked / git-ignored after step 1. Each `bundle/` copy is byte-identical to its live copy after the edit. UTF-16LE matches were all at even byte offsets.

### 3. `src/md2_common.py`

Copies on disk: `src/md2_common.py`, `paper_package/src/md2_common.py`, `bundle/src/md2_common.py`, `bundle/paper_package/src/md2_common.py`. After step 1 only `src/md2_common.py` is not ignored, so only that file was edited. The three ignored copies still hold the old line.

Edit: original line 37, the argument of `PAPER_PDF = Path(...)`, was the raw string `{drive}\Users\<USER>\Desktop\Prevalent mesenchymal drift in aging and disease is reversed by partial reprogramming.pdf`, written with a lowercase drive letter. It is now:

```python
    os.environ.get("MD2_PDF_PATH", str(Path.home() / r"Desktop\Prevalent mesenchymal drift in aging and disease is reversed by partial reprogramming.pdf"))
```

`import os` was missing and was added after `import json`, so the edited line is now line 38. `Path` was already imported. Nothing else in the file changed. CRLF line endings were kept (435 → 436 lines), with no BOM.

Resolved path on this machine, with `MD2_PDF_PATH` unset: identical in length and in every byte except the first. Before, the drive letter was lowercase (`0x63`). After, it is uppercase (`0x43`), because `Path.home()` returns the profile directory with an uppercase drive. Windows paths are case-insensitive, so the two name the same file and `pathlib` compares them as equal, but the strings are not byte-identical. Setting `MD2_PDF_PATH` overrides the default (checked).

What the PDF is: the published article Lu et al., "Prevalent mesenchymal drift in aging and disease is reversed by partial reprogramming", *Cell* 2025, 188:5895–5911, DOI 10.1016/j.cell.2025.07.031 (`TITLE` and `DOI` in `src/md_common.py`). It is a local reference copy and is not in the repo.

Which function uses it: none. `PAPER_PDF` is defined in `src/md2_common.py` and is never referenced anywhere else in the repo: no imports of the name, no `md2_common.PAPER_PDF`, no PDF-reading library. Removing it, or leaving the file absent, changes no output.

### 4. Frozen instruments

| promised in the manuscript | file | git | ignored | contents |
|---|---|---|---|---|
| (a) one frozen fibroblast ruler file | `results/fibro/frozen_ruler_ridge_raw.npz` (1,107,792 bytes) | untracked | no | `w` unit weight vector (norm 1.000000), `coef` un-normalised coefficients, `mu` training mean, `sd` training SD, each float64 over 23,485 genes; plus `intercept` 52.3528, `alpha` 0.01, `gene_id`, `symbol`, `ensembl` |
| (b) young and old GTEx anchors, separate file | `results/toward/anchors.npz` (651,118 bytes) | untracked | no | `c_young` (n = 113, ages 20–39) and `c_old` (n = 231, ages 60–79) over 23,485 genes; `u` unit young-minus-old axis; middle bins 40–59 (n = 308) excluded; note "frozen before any GSE297234 vector". Counts sidecar: `results/toward/anchors.json` |

Both exist outside `paper_package/` and `bundle/`, are not git-ignored (`git check-ignore` exits 1), and will be pushed. The `paper_package/` copies of all three are byte-identical (SHA-256). No copy was made and no `instruments/` folder was created.

Only in `paper_package/`, now ignored: `paper_package/instruments/fibroblast_ruler_weights.csv` (23,485 × 7: `symbol, ensembl, gene_id, weight_unit, coef, mu, sd`) and `paper_package/instruments/fibroblast_ruler_top50.csv`. These are CSV views of the same ruler. They are not needed to meet (a), and they will not be pushed.

README check: **not met.** The "Frozen instruments" table names only directories: `results/fibro/` freeze for the ruler and `results/toward/` for the anchors. Neither `frozen_ruler_ridge_raw.npz` nor `anchors.npz` is named. `results/toward/anchors.json` appears only in the Table 3 donor-count row. README was not edited.

### 5. Re-scan of the push set

Scope: `git ls-files --cached` (121) together with `git ls-files --others --exclude-standard` (2266), giving 2387 files, all present on disk. Raw bytes were searched as UTF-8 and as UTF-16LE, and 121 zip containers (`.npz`, `.docx`, `.xlsx`, …) were also searched member by member.

| pattern | hits |
|---|---:|
| `<USER>`, case-insensitive, excluding the correspondence email | **0** |
| `{drive}` + `Users`, any separator (`\`, `\\`, `/`, `%5C`, `%2F`), any case | **0** |
| any absolute drive path (letter, colon, then a backslash or slash), text files | 7, all false positives |
| any absolute drive path with a path-like tail, binary files | 6, all false positives |

False positives:

- `notebooks/lowdim_l1.ipynb`, `lowdim_l2.ipynb`, `lowdim_l3.ipynb`, `lowdim_l4.ipynb`: a `"k:"` print label in JSON source, where the colon is followed by an escaped quote.
- `notebooks/rerun_stageA_age_axis.ipynb`: prose `...not A:` followed by an escaped `\n`, the same case as above.
- `results/verify/PATH_SCRUB.md` line 74 (2): this report's own description of the drive-path pattern.
- `results/external/frozen_directions.npz`, `results/fibro3/t1_cluster_z.npz`, `results/relabel_check/joint_ams.npz` (2), `results/target/target_W.npz`, `results/target/v1_axes.npz`: random byte runs inside compressed float arrays, not strings.

Real remaining content hits in the push set: **zero**.

**Not covered by this scan: git history.** A push of `master` sends all 4 existing commits (`cb31762`, `aaa6264`, `fa285f1`, `d53662a`). Their committed blobs still contain `<USER>` inside the pre-scrub repo root, in `notebooks/phase0_data_acquisition.ipynb`, `phase1_scores.ipynb`, `rerun_stageA_age_axis.ipynb`, `rerun_stageB_composition.ipynb`, `results/nb_phase0_exec.log`, `results/nb_phase1_exec.err`, `results/phase0_GSE217460_stream_convert.log` and `results/phase1_pseudobulk_build.log`. That is what `git grep` finds in text blobs; UTF-16 logs are treated as binary and may undercount. Commit author and committer are `Age-Identity Benchmark <benchmark@local>` in all 4 commits. The working-tree scrub does not remove these hits. Publishing without them needs a history rewrite or a fresh root commit.

### 6. Size of the push set

2387 files, 244,408,991 bytes (233.09 MiB). These are working-tree sizes measured before this section was appended; the section adds about 10 KB. Git compresses blobs, and `core.autocrlf=true` (system gitconfig, no `.gitattributes`) stores text files with LF instead of CRLF, so the pushed pack will be smaller. No file reaches 50 MB.

Largest 10:

| bytes | MiB | git | file |
|---:|---:|---|---|
| 19,913,805 | 18.99 | u | `results/target/v1_axes.npz` |
| 19,913,223 | 18.99 | u | `results/target/target_W.npz` |
| 16,455,789 | 15.69 | u | `results/relabel_check/lu_meta_HFIB_COMBINED.csv` |
| 13,439,347 | 12.82 | u | `results/relabel_check/joint_ams.npz` |
| 12,103,096 | 11.54 | u | `results/external/frozen_directions.npz` |
| 12,060,910 | 11.50 | u | `results/trajectory/frozen_scores.npz` |
| 9,358,306 | 8.92 | u | `results/south6/stage3_full_table.csv` |
| 8,773,819 | 8.37 | u | `results/target/target_consensus.npz` |
| 6,338,011 | 6.04 | t | `results/tables/phase1b_gene_flags.csv` |
| 5,998,243 | 5.72 | u | `results/md4/t2_state_z.npz` |

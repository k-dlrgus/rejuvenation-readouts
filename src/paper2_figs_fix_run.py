"""PAPER2FIGS-fix: cosmetic re-render of the PAPER2FIGS figures. No plotted value may change.

Run 3 (default). Plan: results/paper2_figs_fix/run3/PLAN.flag. Writes only results/paper2_figs_fix/run3/.
Run 1 (results/paper2_figs_fix/) stopped `reproduction_failed` under a raw-byte comparison;
`--run1` repeats that comparison (prints only, writes nothing outside a deleted temp dir).
Run 2 (results/paper2_figs_fix/run2/) stopped `label_unplaceable` in Fig 4; `--run2` repeats it
into run2/ (Fig 4 labels by search, STOP on an exhausted search; G4 a STOP).
Run 3: Fig 4 labels by hand at the positions in its flag; label checks and G4 are reported,
not STOPs; the only STOPs are reproduction_failed and data_changed.

Comparison rule C (run 2 flag), for a written file F and its blob B at 7a1efb9:
  C1 F with CRLF -> LF equals B exactly; C2 same number of LF; C3 no lone CR in F, no CR in B.

  S1 reproduce: PAPER2FIGS as committed (paper2_figs.FIX = False), rule C on 9 files
  S2 re-render with FIX = True
  G1 rule C on the 7 *_data.csv and angle_check.csv
  G2 label overlap (label_check.csv); G3 Fig1 connectors (label_check.csv)
  new row 76 of paper2_numbers.md, counted from the re-rendered Fig4_data.csv
  G4 paper2_numbers.md differs from the blob only in the lines the flag allows

Usage: python src/paper2_figs_fix_run.py [--run1 | --run2]
"""
from __future__ import annotations

import hashlib
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paper2_figs as pf  # noqa: E402
from config import RESULTS, ROOT  # noqa: E402

BASE = "7a1efb9"
OLD = "results/paper2_figs"
RUN = 2 if "--run2" in sys.argv else 3
OUT = RESULTS / "paper2_figs_fix" / f"run{RUN}"
# Run 3 flag: Fig 4 text anchors (data coordinates), fixed before drawing
FIG4_TEXT = dict(CRI_ZFX=(-22.0, 21.0), CRA_EZH2=(-10.0, 21.0), CRA_E2F3=(-22.0, 3.2), CRI_STAT3=(-20.0, 0.0))
if RUN == 3:
    pf.FIG4_MANUAL = FIG4_TEXT
    pf.LABEL_FALLBACK = True
DATA = ["Fig1_data.csv", "Fig2_data.csv", "Fig3_data.csv", "Fig4_data.csv",
        "FigS1_data.csv", "FigS2_data.csv", "FigS3_data.csv"]


def committed(name):
    return subprocess.run(["git", "show", f"{BASE}:{OLD}/{name}"], cwd=ROOT, capture_output=True,
                          check=True).stdout


def sha(b):
    return hashlib.sha256(b).hexdigest()[:16]


def lf(b):
    return b.replace(b"\r\n", b"\n")


def draw_all(out, fix):
    pf.OUT, pf.FIX = out, fix
    out.mkdir(parents=True, exist_ok=True)
    for f in (pf.fig1, pf.fig2, pf.fig3, pf.fig4, pf.figS1, pf.figS2, pf.figS3):
        f()


def compare(out, names, stage):
    """Rule C of the run 2 flag."""
    rows = []
    for n in names:
        b, f = committed(n), (out / n).read_bytes()
        c1 = lf(f) == b
        c2 = f.count(b"\n") == b.count(b"\n")
        c3 = lf(f).count(b"\r") == 0 and b.count(b"\r") == 0
        rows.append(dict(stage=stage, file=n, bytes_blob=len(b), bytes_written=len(f),
                         lines_blob=b.count(b"\n"), lines_written=f.count(b"\n"),
                         crlf_written=f.count(b"\r\n"), sha256_blob=sha(b), sha256_written_raw=sha(f),
                         sha256_written_lf=sha(lf(f)), C1_lf_identical=c1, C2_line_count=c2,
                         C3_no_other_cr=c3, PASS=c1 and c2 and c3))
    return rows


def failed(rows):
    return [f"{r['file']} (" + ",".join(k[:2] for k in ("C1_lf_identical", "C2_line_count", "C3_no_other_cr")
                                         if not r[k]) + ")" for r in rows if not r["PASS"]]


def new_row():
    """Run 2 flag R1-R3, from the re-rendered Fig4_data.csv."""
    p = OUT / "Fig4_data.csv"
    d = pd.read_csv(p)
    d["_line"] = d.index + 2  # line 1 is the header
    d = d[d.perturbation != "all"]
    sub = d[d.called_early_passage_i.astype(str) == "True"]
    n_sub = len(sub)
    n_called = int((sub.called_gtex_young.astype(str) == "True").sum())
    n_pos = int((sub.delta_gtex_young > 0).sum())
    lines = sub._line.tolist()
    rng, s = [], lines[0]
    for a, b in zip(lines, lines[1:] + [None]):
        if b != a + 1:
            rng.append(f"{s}" if s == a else f"{s}-{a}")
            s = b
    ok = n_called == 0 and n_pos == 80
    row = dict(label="of the 80 called against early-passage cells (origin i): number called against "
                     "GTEx young, and number with delta_gtex_young > 0",
               stated="0 and 80", stored=f"{n_called} and {n_pos}",
               where=f"{pf.rel(p)}:{','.join(rng)}", result="MATCH" if ok else "MISMATCH",
               note=f"Count over plotted rows (arithmetic on stored values). Subset = "
                    f"called_early_passage_i True: {n_sub} rows (stated 80, "
                    f"{'MATCH' if n_sub == 80 else 'MISMATCH'}). Chance baseline for the GTEx count: "
                    f"5th-percentile floor, 4.0 of 80 expected with no effect. Each row's lines in "
                    f"results/seng3/d4_table.csv and results/axang/d3b_table.csv: "
                    f"{pf.rel(OUT / 'new_row_sources.csv')}.")
    pd.DataFrame(dict(fig4_data_line=sub._line, perturbation=sub.perturbation,
                      called_early_passage_i=sub.called_early_passage_i,
                      called_gtex_young=sub.called_gtex_young, delta_gtex_young=sub.delta_gtex_young,
                      source=sub.source)).to_csv(OUT / "new_row_sources.csv", index=False)
    return row, dict(R1=n_sub, R2=n_called, R3=n_pos)


def g4(new_bytes):
    """paper2_numbers.md: only title, Written-by, tally and the appended row 76 may differ."""
    old = committed("paper2_numbers.md").decode("utf-8").split("\n")
    new = lf(new_bytes).decode("utf-8").split("\n")
    i75 = next(i for i, x in enumerate(old) if x.startswith("| 75 |"))
    if len(new) != len(old) + 1 or not new[i75 + 1].startswith("| 76 |"):
        return [f"structure: {len(old)} -> {len(new)} lines, line after row 75 = {new[i75 + 1][:20]!r}"], []
    new_wo = new[:i75 + 1] + new[i75 + 2:]
    allowed = {0: "title", next(i for i, x in enumerate(old) if x.startswith("Written by")): "written-by",
               next(i for i, x in enumerate(old) if x.startswith("**") and "MATCH" in x): "tally"}
    diffs = [(i, allowed.get(i, "NOT ALLOWED")) for i, (a, b) in enumerate(zip(old, new_wo)) if a != b]
    return [f"line {i + 1}: {k}" for i, k in diffs if k == "NOT ALLOWED"], diffs


def main():
    run1 = "--run1" in sys.argv
    OUT.mkdir(parents=True, exist_ok=True)
    # S1 reproduce PAPER2FIGS as committed
    with tempfile.TemporaryDirectory(dir=OUT, prefix="_reproduce_") as tmp:  # inside ROOT; removed on exit
        tmp = Path(tmp)
        draw_all(tmp, fix=False)
        pf.write_numbers(pf.numbers(), pf.angle_check())
        rep = compare(tmp, DATA + ["paper2_numbers.md", "angle_check.csv"], "S1 reproduce (FIX=False)")
        if run1:
            for n in DATA + ["paper2_numbers.md", "angle_check.csv"]:
                print(n, "raw-identical:", committed(n) == (tmp / n).read_bytes())
            return
    pd.DataFrame(rep).to_csv(OUT / "gate_bytes.csv", index=False)
    if failed(rep):
        sys.exit(f"STOP reproduction_failed: {failed(rep)}")

    # S2 re-render with the fixes
    pf.LABEL_CHECK.clear()
    draw_all(OUT, fix=True)

    # G1
    g1 = compare(OUT, DATA, "G1 (FIX=True)")
    pd.DataFrame(rep + g1).to_csv(OUT / "gate_bytes.csv", index=False)
    if failed(g1):
        sys.exit(f"STOP data_changed: {failed(g1)}")

    # G2, G3
    lc = [r for k in ("Fig1", "Fig3C", "Fig4") for r in pf.LABEL_CHECK[k]]
    pd.DataFrame(lc).to_csv(OUT / "label_check.csv", index=False)

    # new row 76, then G4
    R = pf.numbers()
    row, counts = new_row()
    R.append(row)
    pf.write_numbers(R, pf.angle_check())
    p = OUT / "paper2_numbers.md"
    txt = p.read_text(encoding="utf-8")
    txt = txt.replace("# PAPER2FIGS: manuscript numbers against the stored results",
                      "# PAPER2FIGS-fix: manuscript numbers against the stored results", 1)
    txt = txt.replace("Written by `python src/paper2_figs.py`.",
                      "Written by `python src/paper2_figs_fix_run.py`. Rows 1-75 are PAPER2FIGS's; row 76 is new.", 1)
    p.write_text(txt, encoding="utf-8")
    ac = compare(OUT, ["angle_check.csv"], "G1 angle_check (FIX=True)")
    pd.DataFrame(rep + g1 + ac).to_csv(OUT / "gate_bytes.csv", index=False)
    if failed(ac):
        sys.exit(f"STOP data_changed: {failed(ac)}")
    bad, diffs = g4(p.read_bytes())
    print("G4 differing lines (0-based, kind):", diffs, "| not allowed:", bad)
    pd.DataFrame([dict(line=i + 1, kind=k) for i, k in diffs] + [dict(line="", kind=f"PROBLEM: {b}") for b in bad]
                 ).to_csv(OUT / "numbers_diff.csv", index=False)
    if bad and RUN == 2:
        sys.exit(f"STOP numbers_changed: {bad}")

    print("S1:", sum(r["PASS"] for r in rep), "of", len(rep), "PASS")
    print("G1:", sum(r["PASS"] for r in g1 + ac), "of", len(g1 + ac), "PASS")
    print(pd.DataFrame(lc).to_string(index=False))
    print("new row:", row["stated"], "->", row["stored"], row["result"], counts)


if __name__ == "__main__":
    main()

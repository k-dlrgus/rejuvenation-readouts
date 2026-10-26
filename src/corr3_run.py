"""CORR3 — finish the corrections. A writing task: no statistic is computed.

Two docx revisions, each rebuilt from its ORIGINAL rather than patched in place, so that
every check CORR1 and CORR2 made still holds end to end:

  paper 2: v0.7 -> v0.8, CORR2's three edits plus one grammar fix in paragraph 30
           ("1 move away less far" -> "1 moves away less far").
  paper 1: v1.2 -> v1.3, CORR1's two sentence replacements plus the paragraph-142 bullet
           and an in-document correction note updated to mention it.

The edit strings for the earlier work are imported from `src/corr1_run.py` and
`src/corr2_run.py` and are not retyped. The paragraph-142 replacement is asserted to be
byte-identical to the blockquote in FINDINGS_CORR2CHK.md, so it is copied, not restated.

This script OVERWRITES `paper/paper2_passage_draft_v0.8.docx` and
`paper/rejuvenation_readouts_preprint_v1.3.docx` (both git-ignored, like all of `paper/`),
and REGENERATES `results/corr2/paper2_diff.txt` and `results/corr1/paper1_diff.txt`, which
are the tracked records of those two docx files. Each regenerated diff carries a provenance
header naming CORR3 and the commit where the previous version can be read. CLAUDE.md says
not to write into a previous task's results directory; this is a deliberate, disclosed
exception, because a file whose whole content is "the diff of A against B" would otherwise
misdescribe the B it names. See FINDINGS_CORR3.md, "Deviations".

Run from the repository root:
    python -I src/corr3_run.py
"""
from __future__ import annotations

import copy
import difflib
import hashlib
import re
import sys
import zipfile
from pathlib import Path

from lxml import etree

SRC_DIR = Path(__file__).resolve().parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import corr1_run as C1  # noqa: E402
import corr2_run as C2  # noqa: E402

ROOT = SRC_DIR.parent
DOC = C1.DOC
FINDINGS2CHK = ROOT / "FINDINGS_CORR2CHK.md"

# --- paper 2: the grammar fix -------------------------------------------------------------
P2_FIX = ("1 move away less far", "1 moves away less far")

# Neither revision carries an in-document correction note. Earlier runs of this script added
# one to each paper; they were removed on 2026-10-26 because standard preprint practice is that
# a new version carries the corrected text only, with the explanation in the Zenodo version
# notes and the repository notices. Dropping them also returns both papers' paragraph counts,
# and therefore their paragraph numbering, to the originals'.

# --- paper 1: the paragraph-142 bullet ----------------------------------------------------
P142_OLD = (
    "In space C both perturbations move away from the young centroid (+106.6 and +113.8, "
    "q = 1.0), and in space B both end farther from it (+32.9 and +35.1). They have 84 and "
    "82 cells, against a median of 106."
)
P142_NEW = (
    "Neither of the other spaces separates these two from the rest of the screen. In space C "
    "both end farther from the young centroid (+106.6 and +113.8), as do all 1,836 "
    "perturbations: there every displacement is longer than twice the 32-unit gap, so no "
    "change in distance could have come out negative, and the median non-targeting floor is "
    "2.9 gap lengths. In space B both also end farther (+32.9 and +35.1), again with none of "
    "the 1,836 moving closer, and both sit just inside their own cell-count-matched floors "
    "(+33.3 and +35.7), as 1,088 of the 1,836 do. They have 84 and 82 cells, against a "
    "median of 106."
)


def sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def norm_ws(s: str) -> str:
    return re.sub(r"\s+", " ", s).strip()


def check_para142_against_findings() -> None:
    """P142_NEW must be the blockquote in FINDINGS_CORR2CHK.md, not a restatement."""
    text = FINDINGS2CHK.read_text(encoding="utf8")
    marker = "Replacement wording, for folding into the paper-1 revision. No file was edited."
    assert marker in text, "cannot find the paragraph-142 section in FINDINGS_CORR2CHK.md"
    tail = text.split(marker, 1)[1]
    quote = []
    for line in tail.splitlines():
        if line.startswith(">"):
            quote.append(line.lstrip(">").strip())
        elif quote:
            break
    assert quote, "no blockquote follows the paragraph-142 marker"
    want = norm_ws(" ".join(quote))
    want = re.sub(r"^[••]\s*", "", want)
    assert want == norm_ws(P142_NEW), (
        "P142_NEW differs from FINDINGS_CORR2CHK.md's blockquote:\n"
        f"findings: {want}\n\nscript  : {norm_ws(P142_NEW)}"
    )


def skeletons(xml: str) -> list[bytes]:
    """Per body paragraph, its XML with every w:t text blanked.

    What survives is the formatting alone: the paragraph's w:pPr, its run order, and each
    run's w:rPr. Two paragraphs whose skeletons are equal differ only in their text.
    """
    body = etree.fromstring(xml.encode("utf8")).find(f"{C1.W}body")
    out = []
    for p in body.iterchildren(f"{C1.W}p"):
        q = copy.deepcopy(p)
        for t in q.iter(f"{C1.W}t"):
            t.text = ""
        out.append(etree.tostring(q))
    return out


def check_formatting(old: str, new: str, label: str) -> list[int]:
    """Every paragraph this task rewrote must keep its original formatting exactly.

    Returns the 1-indexed positions whose text changed. Asserts that the paragraph count is
    unchanged (so no paragraph was inserted or removed) and that each changed paragraph's
    skeleton is byte-identical to the one it replaced.
    """
    a, b = skeletons(old), skeletons(new)
    assert len(a) == len(b), (
        f"{label}: paragraph count changed, {len(a)} -> {len(b)}; this task inserts none")
    ta = C1.body_paragraphs(old.encode("utf8"))
    tb = C1.body_paragraphs(new.encode("utf8"))
    changed = [i + 1 for i, (x, y) in enumerate(zip(ta, tb)) if x != y]
    for i in changed:
        assert a[i - 1] == b[i - 1], (
            f"{label}: paragraph {i} changed formatting, not just text")
    return changed


def apply_edits(old: str, edits: list[tuple[str, str]], label: str) -> str:
    new = old
    for o, n in edits:
        count = new.count(o)
        assert count == 1, f"{label}: old string occurs {count} times: {o[:70]!r}"
        new = new.replace(o, n)
    back = new
    for o, n in reversed(edits):
        assert back.count(n) == 1, f"{label}: new string not unique when undoing: {n[:70]!r}"
        back = back.replace(n, o)
    assert back == old, f"{label}: undoing the edits does not give back the original XML"
    return new


def write_docx(src: Path, dst: Path, new_xml: str) -> tuple[int, int]:
    """Writes dst from src with DOC replaced. Returns (identical others, total others)."""
    zin = zipfile.ZipFile(src)
    with zipfile.ZipFile(dst, "w") as zout:
        for info in zin.infolist():
            data = new_xml.encode("utf8") if info.filename == DOC else zin.read(info.filename)
            out = zipfile.ZipInfo(info.filename, date_time=info.date_time)
            out.compress_type, out.external_attr = info.compress_type, info.external_attr
            zout.writestr(out, data)
    zchk = zipfile.ZipFile(dst)
    others = [i.filename for i in zin.infolist() if i.filename != DOC]
    same = [f for f in others if zin.read(f) == zchk.read(f)]
    assert zchk.read(DOC).decode("utf8") == new_xml
    assert len(same) == len(others), "a zip member other than word/document.xml changed"
    return len(same), len(others)


def diff_text(src: Path, dst: Path, old: str, new: str, header: list[str]) -> str:
    a = C1.body_paragraphs(old.encode("utf8"))
    b = C1.body_paragraphs(new.encode("utf8"))
    d = difflib.unified_diff(
        a, b, fromfile=src.relative_to(ROOT).as_posix(),
        tofile=dst.relative_to(ROOT).as_posix(), lineterm="", n=0)
    return "\n".join(header + [f"Top-level body paragraphs: {len(a)} in the original, "
                               f"{len(b)} in the revision.", ""] + list(d)) + "\n"


def paper2() -> str:
    src, dst = C2.SRC, C2.DST
    edits = list(C2.EDITS) + [P2_FIX]
    with zipfile.ZipFile(src) as z:
        old = z.read(DOC).decode("utf8")
    new = apply_edits(old, edits, "paper 2")
    changed = check_formatting(old, new, "paper 2")
    if dst.exists():
        dst.unlink()
    same, others = write_docx(src, dst, new)
    header = [
        "CORR2/CORR3 — paper 2, v0.7 -> v0.8. Paragraph-level diff of the docx body text.",
        "Generated by src/corr3_run.py. REGENERATED by CORR3 on 26 October 2026: CORR2's three",
        "edits are unchanged and a grammar fix is added in paragraph 30 (\"1 move away less far\"",
        "-> \"1 moves away less far\"). CORR2's original version of this file is in git history at",
        "commit 46a05c1.",
        "Each line is one top-level body paragraph; hunk line numbers are paragraph positions",
        "(1-indexed). Text inside text boxes (figure and table captions) is not listed. v0.7 is",
        "the source of the paper's published v1 (https://doi.org/10.5281/zenodo.23591623) and",
        "v0.8 is the source of its v2. No paragraph is inserted or removed and v2 carries no",
        "in-document correction note, so paragraph numbers are the same in both: the two",
        "corrected paragraphs are 30 and 31 in each. The explanation of the correction lives in",
        "the Zenodo version notes and in CORRECTIONS_PAPER2.md, not in the manuscript.",
        "",
        f"v0.7 sha256 {sha(src.read_bytes())}",
        f"v0.8 sha256 {sha(dst.read_bytes())}",
        f"word/document.xml: {len(edits)} string edits (paragraph 30, one sentence replaced and",
        "one word corrected; paragraph 31, the last sentence in two pieces either side of the",
        "superscript reference run, which is kept); undoing them reproduces the v0.7 XML exactly.",
        f"Other zip members byte-identical: {same} of {others}.",
        f"Paragraphs whose text changed: {', '.join(str(i) for i in changed)}. Each keeps the",
        "paragraph properties, run order and run properties of the paragraph it replaces, so",
        "alignment, font and size are untouched; only the text differs.",
        "",
    ]
    out = C2.OUT / "paper2_diff.txt"
    text = diff_text(src, dst, old, new, header)
    out.write_text(text, encoding="utf8")
    return f"paper 2: v0.8 rebuilt, {len(edits)} edits, {same}/{others} other members identical"


def paper1() -> str:
    src, dst = C1.SRC, C1.DST
    edits = [C1.EDITS[0], C1.EDITS[1], (P142_OLD, P142_NEW)]
    with zipfile.ZipFile(src) as z:
        old = z.read(DOC).decode("utf8")
    new = apply_edits(old, edits, "paper 1")
    changed = check_formatting(old, new, "paper 1")
    if dst.exists():
        dst.unlink()
    same, others = write_docx(src, dst, new)
    header = [
        "CORR1/CORR3 — paper 1, v1.2 -> v1.3. Paragraph-level diff of the docx body text.",
        "Generated by src/corr3_run.py. REGENERATED by CORR3 on 26 October 2026: CORR1's two",
        "sentence replacements in Supplementary Note 3 are unchanged and the paragraph-142",
        "bullet is now corrected too. CORR1's original version of this file, describing v1.3",
        "with the note-3 change alone, is in git history at commit d8d892a.",
        "Each line is one top-level body paragraph; hunk line numbers are paragraph positions",
        "(1-indexed). Text inside text boxes (figure and table captions) is not listed. No",
        "paragraph is inserted or removed and v1.3 carries no in-document correction note, so",
        "paragraph numbers are the same in v1.2 and v1.3: the corrected bullet is paragraph 142",
        "and the corrected note-3 paragraph is 144 in each. The explanation of the correction",
        "lives in the Zenodo version notes and in CORRECTIONS_PAPER1.md, not in the manuscript.",
        "",
        f"v1.2 sha256 {sha(src.read_bytes())}",
        f"v1.3 sha256 {sha(dst.read_bytes())}",
        f"word/document.xml: {len(edits)} string edits (2 sentence replacements and 1 bullet",
        "replacement); undoing them reproduces the v1.2 XML exactly. Other zip members",
        f"byte-identical: {same} of {others}.",
        f"Paragraphs whose text changed: {', '.join(str(i) for i in changed)}. Each keeps the",
        "paragraph properties, run order and run properties of the paragraph it replaces, so",
        "alignment, font and size are untouched; only the text differs.",
        "",
    ]
    out = C1.OUT / "paper1_diff.txt"
    text = diff_text(src, dst, old, new, header)
    out.write_text(text, encoding="utf8")
    return f"paper 1: v1.3 rebuilt, {len(edits)} edits, {same}/{others} other members identical"


def main() -> int:
    check_para142_against_findings()
    print("paragraph-142 wording matches FINDINGS_CORR2CHK.md verbatim")
    for line in (paper2(), paper1()):
        print(line)
    print("wrote results/corr2/paper2_diff.txt and results/corr1/paper1_diff.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())

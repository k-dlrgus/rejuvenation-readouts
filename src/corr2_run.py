"""CORR2 — revised paper-2 source and its diff, folding in CORR2CHK's wording.

A writing task, not an analysis: no statistic is computed. Every number in the
replacement text is copied from FINDINGS_CORR2CHK.md (sections "Corrected
wording (step 3)", "A1b — the Southard clause of paragraph 31" and "C5 — the
young-donor clause of paragraph 30") or from results/corr2chk/.

Paper 2 is not public. This is a pre-publication revision of a draft, not a
correction notice, so — unlike CORR1 — no notice paragraph is inserted and no
public notice file is written.

Reads paper/paper2_passage_draft_v0.7.docx and writes a new file beside it,
paper/paper2_passage_draft_v0.8.docx. The original is not touched. paper/ is
git-ignored, so neither docx is tracked; the tracked record of the change is
results/corr2/paper2_diff.txt, which this script writes.

The change is three string edits to word/document.xml, in two body paragraphs:

  - paragraph 30, one sentence replaced (CORR2CHK's C5, the young-donor clause);
  - paragraph 31, the last sentence replaced in two pieces, because a
    superscript reference run ("5", the Southard citation) sits inside it and is
    kept where it is: the text before the citation is edit 2 and the text after
    it is edit 3.

The script checks that each old string occurs exactly once, that every other zip
member is byte-identical, and that undoing the three edits gives back the
original XML.

Run from the repository root:
    python -I src/corr2_run.py
"""
from __future__ import annotations

import difflib
import hashlib
import sys
import zipfile
from pathlib import Path

from lxml import etree

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "paper" / "paper2_passage_draft_v0.7.docx"
DST = ROOT / "paper" / "paper2_passage_draft_v0.8.docx"
OUT = ROOT / "results" / "corr2"
DOC = "word/document.xml"
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"

# Paragraph 30 — CORR2CHK's C5. The clause "and all 80 move away from them" is
# the part that fails; the sentence is replaced whole so the replacement reads
# as one sentence pair.
OLD1 = (
    "None of the 80 perturbations called closer to early-passage cells is called closer to young "
    "donors, and all 80 move away from them."
)
NEW1 = (
    "None of the 80 perturbations called closer to early-passage cells is called closer to young "
    "donors, and all 80 have a positive change in distance to them. Against no-effect samples of "
    "the same size, 56 of the 80 move away further than the 95th percentile, 23 lie inside the "
    "no-effect spread and 1 move away less far than the median no-effect sample, so the departure "
    "is established for 56 of them and not for the rest."
)
# Paragraph 31 — CORR2CHK's A1b. The failing clause is "and move away from young
# donors in the Southard screen"; CORR2CHK's replacement text also restates the
# activation-only clause that closes the sentence, so the whole sentence from
# "Separately," onward is replaced and nothing of it is left duplicated.
OLD2 = (
    "Separately, the four genes are also uncalled and move away from young donors in the Southard "
    "screen"
)
NEW2 = "Separately, the four genes are also uncalled in the Southard screen"
OLD3 = (
    ", a different neonatal line; that screen used activation only, so for STAT3 and ZFX it tests "
    "the opposite perturbation."
)
NEW3 = (
    ", a different neonatal line, where the test detects a quarter-gap approach toward young "
    "donors at their cell counts with power 0.925 to 1.000. Their changes in distance to young "
    "donors are positive, but so are those of non-targeting cells with no effect at the same cell "
    "counts: the no-effect median runs from +1.3 to +2.7 against the four factors’ +0.5 to "
    "+2.4. Three of the four move away less far than the median no-effect sample of their own "
    "size, and none exceeds its 95th percentile, so that screen shows no approach and does not "
    "show a departure either. That screen used activation only, so for STAT3 and ZFX it tests the "
    "opposite perturbation."
)
EDITS = [(OLD1, NEW1), (OLD2, NEW2), (OLD3, NEW3)]


def sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def body_paragraphs(xml: bytes) -> list[str]:
    """Text of each top-level body paragraph, 1-indexed by position; text-box content excluded."""
    body = etree.fromstring(xml).find(f"{W}body")
    out = []
    for p in body.iterchildren(f"{W}p"):
        ts = [t.text or "" for t in p.iter(f"{W}t")
              if not any(a.tag == f"{W}txbxContent" for a in t.iterancestors())]
        out.append("".join(ts))
    return out


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    if DST.exists():
        print(f"{DST} exists; not overwriting", file=sys.stderr)
        return 1
    zin = zipfile.ZipFile(SRC)
    old = zin.read(DOC).decode("utf8")
    new = old
    for o, n in EDITS:
        assert new.count(o) == 1, f"old string occurs {new.count(o)} times: {o[:60]}"
        new = new.replace(o, n)
    back = new
    for o, n in reversed(EDITS):
        assert back.count(n) == 1
        back = back.replace(n, o)
    assert back == old, "undoing the edits does not give back the original XML"

    with zipfile.ZipFile(DST, "w") as zout:
        for info in zin.infolist():
            data = new.encode("utf8") if info.filename == DOC else zin.read(info.filename)
            # a fresh ZipInfo: writestr rewrites header_offset on the one it is given
            out = zipfile.ZipInfo(info.filename, date_time=info.date_time)
            out.compress_type, out.external_attr = info.compress_type, info.external_attr
            zout.writestr(out, data)
    zchk = zipfile.ZipFile(DST)
    same = [i.filename for i in zin.infolist()
            if i.filename != DOC and zin.read(i.filename) == zchk.read(i.filename)]
    n_other = len(zin.infolist()) - 1
    assert len(same) == n_other and zchk.read(DOC).decode("utf8") == new

    a = body_paragraphs(old.encode("utf8"))
    b = body_paragraphs(new.encode("utf8"))
    assert len(a) == len(b), "paragraph count changed; this revision inserts no paragraph"
    diff = difflib.unified_diff(
        a, b,
        fromfile=SRC.relative_to(ROOT).as_posix(), tofile=DST.relative_to(ROOT).as_posix(),
        lineterm="", n=0)
    header = [
        "CORR2 — paper 2, v0.7 -> v0.8. Paragraph-level diff of the docx body text.",
        "Generated by src/corr2_run.py. Each line is one top-level body paragraph; hunk line",
        "numbers are paragraph positions (1-indexed). Text inside text boxes (figure and table",
        "captions) is not listed. No paragraph is inserted or removed, so positions do not shift.",
        "Paper 2 is not public: this is a pre-publication revision, not a correction notice, and",
        "no notice paragraph is added to the draft.",
        "",
        f"v0.7 sha256 {sha(SRC.read_bytes())}",
        f"v0.8 sha256 {sha(DST.read_bytes())}",
        "word/document.xml: 3 string edits (paragraph 30, one sentence; paragraph 31, the last",
        "sentence in two pieces, either side of the superscript reference run, which is kept);",
        "undoing them reproduces the v0.7 XML exactly. Other zip members byte-identical: "
        f"{len(same)} of {n_other}.",
        f"Top-level body paragraphs: {len(a)} in v0.7, {len(b)} in v0.8.",
        "",
    ]
    text = "\n".join(header + list(diff)) + "\n"
    (OUT / "paper2_diff.txt").write_text(text, encoding="utf8")
    sys.stdout.buffer.write(text.encode("utf8"))
    return 0


if __name__ == "__main__":
    sys.exit(main())

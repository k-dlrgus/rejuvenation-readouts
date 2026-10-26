"""CORR1 — revised paper-1 source and its diff, for the correction notice.

A writing task, not an analysis: no statistic is computed. Every number in the
replacement text is copied from FINDINGS_FLOORCHK.md (section "What no longer
holds in paper 1") or from results/floorchk/.

Reads paper/rejuvenation_readouts_preprint_v1.2.docx and writes a new file
beside it, paper/rejuvenation_readouts_preprint_v1.3.docx. The original is not
touched. paper/ is git-ignored, so neither docx is tracked; the tracked record
of the change is results/corr1/paper1_diff.txt, which this script writes.

The change is three string edits to word/document.xml: two sentences of the
Supplementary Note 3 paragraph that begins "Space C remains the screen's
pre-registered primary" are replaced, and one correction line is inserted as a
new paragraph directly under the Supplementary Note 3 heading. The script
checks that each old string occurs exactly once, that every other zip member is
byte-identical, and that undoing the three edits gives back the original XML.

Run from the repository root:
    python src/corr1_run.py
"""
from __future__ import annotations

import difflib
import hashlib
import sys
import zipfile
from pathlib import Path

from lxml import etree

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "paper" / "rejuvenation_readouts_preprint_v1.2.docx"
DST = ROOT / "paper" / "rejuvenation_readouts_preprint_v1.3.docx"
OUT = ROOT / "results" / "corr1"
DOC = "word/document.xml"
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"

OLD1 = (
    "In space C, for the GTEx old-to-young comparison, with the GTEx old centroid (60 to 79 years) "
    "as the origin and the GTEx young centroid (20 to 39 years) as the target, every one of the "
    "1,836 perturbations moved more than twice the distance between the two centroids (smallest "
    "ratio 2.10; per-perturbation values in results/south6/ in the repository), so from that origin "
    "no approach to that target could have been registered whatever its direction."
)
NEW1 = (
    "In space C, for the GTEx old-to-young comparison, with the GTEx old centroid (60 to 79 years) "
    "as the origin and the GTEx young centroid (20 to 39 years) as the target, every one of the "
    "1,836 measured displacements is longer than twice the distance between the two centroids "
    "(smallest ratio 2.10; per-perturbation values in results/south6/ in the repository), so from "
    "that origin the test could not register an approach to that target whatever the direction."
)
OLD2 = (
    "This holds for that origin and target, which sit only 32 units apart, and is not a statement "
    "about the screen in general; it is the same limitation the main text discloses for the "
    "same-platform test."
)
NEW2 = (
    "This holds for that origin and target, which sit only 32 units apart, and is not a statement "
    "about the screen in general. It reflects depth, not effect size: non-targeting cells with no "
    "effect, sampled at the same group sizes, also read more than twice that distance (2.1 to 12.9 "
    "at the 95th percentile). After the noise is subtracted, 18 perturbations are estimated to move "
    "more than twice the distance (about 5 expected by chance), and the median estimate is 0.6 "
    "times it (results/floorchk/ in the repository). The space C null for the screen is therefore a "
    "null of an instrument that could not see approach at this depth, not evidence that the "
    "perturbations overshoot, and it is not the same limitation as the one the main text discloses "
    "for the same-platform test, whose length statement survives the same check."
)
HEADING = (
    "<w:t>Supplementary Note 3. The same-platform test in three coordinate systems</w:t></w:r></w:p>"
)
NOTE = (
    "Correction (v1.3): two sentences on the transcription-factor screen's displacement lengths in "
    "space C, in the last paragraph of this note, were corrected on 26 October 2026; see "
    "CORRECTIONS_PAPER1.md in the repository."
)
NOTE_PARA = (
    '<w:p><w:pPr><w:spacing w:after="160"/><w:jc w:val="both"/></w:pPr><w:r><w:rPr>'
    '<w:rFonts w:eastAsia="Times New Roman"/><w:i/></w:rPr><w:t>' + NOTE + "</w:t></w:r></w:p>"
)
EDITS = [(OLD1, NEW1), (OLD2, NEW2), (HEADING, HEADING + NOTE_PARA)]


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
    diff = difflib.unified_diff(
        a, b,
        fromfile=SRC.relative_to(ROOT).as_posix(), tofile=DST.relative_to(ROOT).as_posix(),
        lineterm="", n=0)
    header = [
        "CORR1 — paper 1, v1.2 -> v1.3. Paragraph-level diff of the docx body text.",
        "Generated by src/corr1_run.py. Each line is one top-level body paragraph; hunk line",
        "numbers are paragraph positions (1-indexed). Text inside text boxes (figure and table",
        "captions) is not listed. v1.3 has one more paragraph, so later positions shift by one.",
        "",
        f"v1.2 sha256 {sha(SRC.read_bytes())}",
        f"v1.3 sha256 {sha(DST.read_bytes())}",
        f"word/document.xml: 3 string edits (2 sentence replacements, 1 inserted paragraph);",
        f"undoing them reproduces the v1.2 XML exactly. Other zip members byte-identical: "
        f"{len(same)} of {n_other}.",
        f"Top-level body paragraphs: {len(a)} in v1.2, {len(b)} in v1.3.",
        "",
    ]
    text = "\n".join(header + list(diff)) + "\n"
    (OUT / "paper1_diff.txt").write_text(text, encoding="utf8")
    sys.stdout.buffer.write(text.encode("utf8"))
    return 0


if __name__ == "__main__":
    sys.exit(main())

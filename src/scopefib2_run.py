"""SCOPEFIB2: metadata-only record fetcher for one unverified entry left by SCOPEFIB.

Scope is the ADAM12 study's human fibroblast Perturb-seq and sibling human accessions of
that study only (results/scopefib2/PREREG.flag). No matrix download, no statistic beyond the
cell counts bar C requires. Every request is logged with its URL and HTTP status to
results/scopefib2/url_records.json so each fact can be traced to the text it was read from.

Usage:
  python src/scopefib2_run.py url <label> <url>          -> url_records.json (<=50 MB, text)
  python src/scopefib2_run.py post <label> <url> <json>  -> url_records.json (POST body)
  python src/scopefib2_run.py file <label> <url> <out>   -> opened/<file> (<=50 MB, binary)
  python src/scopefib2_run.py series GSE... [...]        -> geo_records.json
  python src/scopefib2_run.py reproduce                  -> reproduction.csv (pre-reg check)
"""
from __future__ import annotations

import csv
import json
import sys
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

# The Windows console here is cp949; record text is UTF-8 (Chinese repository pages).
# Printing must not be able to abort a fetch that already succeeded and was logged.
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except AttributeError:  # pragma: no cover - non-reconfigurable stream
        pass

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "results" / "scopefib2"
UA = "AgeIdentitySeparabilityBenchmark/scopefib2 (metadata only)"
MAX_BYTES = 50 * 1024 * 1024


def append(name: str, record: dict) -> None:
    path = OUT / name
    log = json.loads(path.read_text(encoding="utf-8")) if path.exists() else []
    record["fetched_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    log.append(record)
    path.write_text(json.dumps(log, indent=1, ensure_ascii=False), encoding="utf-8")


def fetch(label: str, url: str, data: bytes | None = None, ctype: str | None = None,
          log: str = "url_records.json", store: bool = True) -> dict:
    """GET or POST, log URL + HTTP status (including errors), never exceed 50 MB."""
    headers = {"User-Agent": UA, "Accept": "*/*"}
    if ctype:
        headers["Content-Type"] = ctype
    req = urllib.request.Request(url, data=data, headers=headers)
    rec: dict = {"label": label, "url": url, "method": "POST" if data else "GET"}
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            rec["status"] = r.status
            n = r.headers.get("Content-Length")
            if n is not None and int(n) > MAX_BYTES:
                rec["error"] = f"refused: {n} bytes (>50 MB)"
                append(log, rec)
                return rec
            raw = r.read(MAX_BYTES + 1)
            if len(raw) > MAX_BYTES:
                rec["error"] = "refused: exceeded 50 MB while reading"
                append(log, rec)
                return rec
            rec["bytes"] = len(raw)
            rec["content_type"] = r.headers.get("Content-Type", "")
            rec["final_url"] = r.url
            text = raw.decode("utf-8", errors="replace")
    except urllib.error.HTTPError as e:
        rec["status"] = e.code
        body = e.read(200000).decode("utf-8", errors="replace")
        rec["error"] = f"HTTPError {e.code} {e.reason}"
        rec["text"] = body if store else ""
        append(log, rec)
        return rec
    except Exception as e:  # noqa: BLE001 - URLError, timeouts, TLS
        rec["status"] = None
        rec["error"] = f"{type(e).__name__}: {e}"
        append(log, rec)
        return rec
    if store:
        rec["text"] = text
    append(log, rec)
    rec["text"] = text
    return rec


def fetch_binary(label: str, url: str, out: Path) -> dict:
    """Save a binary record file (e.g. the study's public metadata xlsx) under
    results/scopefib2/opened/, logging URL + status like every other request. The
    50 MB cap applies. Downloaded files are untrusted input: they are written to
    their own directory and read by a separate script."""
    headers = {"User-Agent": UA, "Accept": "*/*"}
    rec: dict = {"label": label, "url": url, "method": "GET", "saved_to": str(out)}
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers=headers),
                                    timeout=120) as r:
            rec["status"] = r.status
            n = r.headers.get("Content-Length")
            if n is not None and int(n) > MAX_BYTES:
                rec["error"] = f"refused: {n} bytes (>50 MB)"
                append("url_records.json", rec)
                return rec
            raw = r.read(MAX_BYTES + 1)
            if len(raw) > MAX_BYTES:
                rec["error"] = "refused: exceeded 50 MB while reading"
                append("url_records.json", rec)
                return rec
            rec["bytes"] = len(raw)
            rec["content_type"] = r.headers.get("Content-Type", "")
            rec["content_disposition"] = r.headers.get("Content-Disposition", "")
    except urllib.error.HTTPError as e:
        rec["status"] = e.code
        rec["error"] = f"HTTPError {e.code} {e.reason}"
        append("url_records.json", rec)
        return rec
    except Exception as e:  # noqa: BLE001
        rec["status"] = None
        rec["error"] = f"{type(e).__name__}: {e}"
        append("url_records.json", rec)
        return rec
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(raw)
    append("url_records.json", rec)
    return rec


def geo_text(acc: str) -> dict:
    q = urllib.parse.urlencode({"acc": acc, "targ": "self", "form": "text", "view": "quick"})
    return fetch(acc, "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?" + q,
                 log="geo_records.json")


def reproduce() -> int:
    """Pre-registered reproduction of the stored SCOPEFIB record. See PREREG.flag."""
    rows: list[dict] = []

    def chk(key, expected, observed):
        rows.append({"check": key, "expected": expected, "observed": observed,
                     "match": "yes" if str(expected) == str(observed) else "NO"})

    with open(REPO / "results/scopefib/datasets.csv", newline="", encoding="utf-8") as fh:
        data = list(csv.DictReader(fh))
    chk("R1 datasets.csv data rows", 23, len(data))

    scoping = (REPO / "SCOPING_FIBROSIS.md").read_text(encoding="utf-8")
    progress = (REPO / "PROGRESS_SCOPEFIB.md").read_text(encoding="utf-8")
    readings = [ln.strip().strip("- ").strip("`") for ln in progress.split("\n")
                if ln.strip().startswith("- `")]
    chk("R2a n readings in PROGRESS_SCOPEFIB.md", 2, len(readings))
    for i, r in enumerate(readings, 1):
        chk(f"R2b reading {i} also in SCOPING_FIBROSIS.md verbatim", True, r in scoping)
        rows.append({"check": f"R2c reading {i} text", "expected": r, "observed": r,
                     "match": "yes"})

    unver = scoping.split("## unverified", 1)[-1]
    chk("R3 ADAM12 in unverified block", True,
        "ADAM12" in unver and "10.1016/j.ccell.2025.12.018" in unver)
    chk("R3b GSE278990 named mouse-only in unverified block", True,
        "GSE278990 contains mouse data only" in unver)

    src = (REPO / "results/scopefib/PREREG.flag").read_text(encoding="utf-8").split("\n")
    carried = "\n".join(src[9:87])
    new = (REPO / "results/scopefib2/PREREG.flag").read_text(encoding="utf-8")
    chk("R4 carried lines 10-87 byte-identical in SCOPEFIB2 flag", True, carried in new)
    chk("R4b carried block n lines", 78, len(carried.split("\n")))

    with open(OUT / "reproduction.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["check", "expected", "observed", "match"])
        w.writeheader()
        w.writerows(rows)
    bad = [r for r in rows if r["match"] == "NO"]
    for r in rows:
        print(f"{r['match']:>3}  {r['check']}")
    print(f"{len(rows)} checks, {len(bad)} mismatches")
    return 1 if bad else 0


def main(argv: list[str]) -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    mode, args = argv[0], argv[1:]
    if mode == "url":
        r = fetch(args[0], args[1])
        print(f"[{r.get('status')}] {r['url']}")
        print(r.get("text", r.get("error", ""))[:200000])
    elif mode == "post":
        r = fetch(args[0], args[1], data=args[2].encode(), ctype="application/json")
        print(f"[{r.get('status')}] {r['url']}")
        print(r.get("text", r.get("error", ""))[:200000])
    elif mode == "file":
        r = fetch_binary(args[0], args[1], Path(args[2]))
        print(f"[{r.get('status')}] {r.get('bytes')} bytes {r.get('content_type','')} "
              f"{r.get('error','')} -> {r.get('saved_to')}")
    elif mode == "series":
        for acc in args:
            r = geo_text(acc)
            print(f"===== {acc} [{r.get('status')}]\n{r.get('text', r.get('error', ''))}")
    elif mode == "reproduce":
        return reproduce()
    else:
        raise SystemExit(__doc__)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

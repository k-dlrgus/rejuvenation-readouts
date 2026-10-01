"""Parallel HTTP Range download of the SEA-AD MTG h5ad.

Stage 0 only needs obs (range-fetched). Projection needs raw/X (~50 GB).
A single-connection GET saturates around 2 MB/s; N parallel ranges are
the same published URL, not a substitute dataset.

Usage: python src/external_download.py
"""
from __future__ import annotations

import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from external_common import EXT_RAW, load_json, EXT_DIR, StopStep  # noqa: E402
from download import remote_size, CHUNK  # noqa: E402

N_PARTS = 8
RETRIES = 20


def _part_path(dest: Path, i: int) -> Path:
    return dest.with_name(dest.name + f".part{i:02d}")


def _fetch_range(url: str, dest_part: Path, start: int, end: int, label: str) -> int:
    """Inclusive start, exclusive end. Resume if dest_part already has a prefix of this range."""
    need = end - start
    have = dest_part.stat().st_size if dest_part.exists() else 0
    if have > need:
        dest_part.unlink()
        have = 0
    if have == need:
        print(f"[{label}] cached {need/1e9:.2f} GB", flush=True)
        return need
    attempt = 0
    while attempt < RETRIES:
        attempt += 1
        byte0 = start + have
        headers = {"Range": f"bytes={byte0}-{end - 1}"}
        try:
            with requests.get(url, headers=headers, stream=True, timeout=120) as r:
                if r.status_code != 206:
                    raise OSError(f"{label}: expected 206, got {r.status_code}")
                r.raise_for_status()
                t0 = time.time()
                last = t0
                mode = "ab" if have > 0 else "wb"
                with open(dest_part, mode) as f:
                    for chunk in r.iter_content(CHUNK):
                        if not chunk:
                            continue
                        f.write(chunk)
                        have += len(chunk)
                        now = time.time()
                        if now - last > 20:
                            rate = have / max(now - t0, 1e-9) / 1e6
                            print(f"  [{label}] {have/1e9:.2f}/{need/1e9:.2f} GB  ({rate:.1f} MB/s)", flush=True)
                            last = now
            if have == need:
                print(f"[{label}] done {need/1e9:.2f} GB", flush=True)
                return need
            print(f"[{label}] short ({have} vs {need}); retry", flush=True)
        except (requests.RequestException, OSError) as e:
            print(f"[{label} retry {attempt}/{RETRIES}] {e}", flush=True)
            time.sleep(min(60, 5 * attempt))
            have = dest_part.stat().st_size if dest_part.exists() else 0
    raise RuntimeError(f"{label} failed")


def assemble(dest: Path, parts: list[Path], total: int):
    """Concatenate parts without a second 50 GB temp copy (disk is tight)."""
    print(f"[assemble] writing {dest} from {len(parts)} parts", flush=True)
    first = parts[0]
    if dest.exists():
        dest.unlink()
    os.replace(first, dest)
    with open(dest, "ab") as out:
        for p in parts[1:]:
            print(f"[assemble] append {p.name} ({p.stat().st_size/1e9:.2f} GB)", flush=True)
            with open(p, "rb") as inp:
                while True:
                    buf = inp.read(CHUNK)
                    if not buf:
                        break
                    out.write(buf)
            p.unlink()
    got = dest.stat().st_size
    if got != total:
        raise RuntimeError(f"assembled size {got} != remote {total}")
    print(f"[done] {dest} ({got/1e9:.2f} GB)", flush=True)


def main():
    summ = load_json(EXT_DIR / "stage0_summary.json")
    url = summ["h5ad_url"]
    expected = int(summ["h5ad_bytes"])
    dest = EXT_RAW / url.rstrip("/").rsplit("/", 1)[-1]
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"[dl] {url}\n     -> {dest}\n     expected={expected}", flush=True)
    if dest.exists() and dest.stat().st_size == expected:
        print(f"[cached] complete {dest.stat().st_size/1e9:.2f} GB", flush=True)
        return str(dest)
    if dest.exists() and dest.stat().st_size != expected:
        print(f"[warn] incomplete {dest.stat().st_size} vs {expected}; removing", flush=True)
        dest.unlink()
    live = remote_size(url)
    if live is None:
        raise StopStep("download", "HEAD/Range did not return a size for the h5ad")
    if int(live) != expected:
        raise StopStep(
            "download",
            f"live size {live} != Stage 0 h5ad_bytes {expected}. Not substituting a different object.",
            dict(live=live, expected=expected),
        )
    spans = []
    for i in range(N_PARTS):
        a = (expected * i) // N_PARTS
        b = (expected * (i + 1)) // N_PARTS
        spans.append((i, a, b))
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=N_PARTS) as ex:
        futs = [
            ex.submit(_fetch_range, url, _part_path(dest, i), a, b, f"p{i:02d}")
            for i, a, b in spans
        ]
        for fut in as_completed(futs):
            fut.result()
    parts = [_part_path(dest, i) for i, _, _ in spans]
    assemble(dest, parts, expected)
    print(f"[dl] elapsed {time.time()-t0:.0f}s", flush=True)
    return str(dest)


if __name__ == "__main__":
    try:
        main()
    except StopStep as e:
        print(f"[STOP] {e.step}: {e.message}", flush=True)
        sys.exit(2)

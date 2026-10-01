"""Resumable, cache-aware downloader. Never re-downloads a complete file.

Usage: python src/download.py URL [URL ...] [--out data/raw]
A file is considered complete when its on-disk size equals the server's
Content-Length. Partial files are resumed with HTTP Range requests.
"""
import argparse
import os
import sys
import time
import requests

if hasattr(sys.stdout, "reconfigure"):
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

CHUNK = 1 << 22  # 4 MiB


def remote_size(url: str) -> int | None:
    try:
        r = requests.head(url, allow_redirects=True, timeout=60)
        if r.status_code == 200 and "Content-Length" in r.headers:
            return int(r.headers["Content-Length"])
    except requests.RequestException:
        pass
    # Fallback: GET with range 0-0
    try:
        r = requests.get(url, headers={"Range": "bytes=0-0"}, stream=True, timeout=60)
        cr = r.headers.get("Content-Range")
        if cr and "/" in cr:
            return int(cr.rsplit("/", 1)[1])
    except requests.RequestException:
        pass
    return None


def download(url: str, out_dir: str, retries: int = 20) -> str:
    os.makedirs(out_dir, exist_ok=True)
    fname = os.path.join(out_dir, url.rstrip("/").rsplit("/", 1)[-1])
    total = remote_size(url)
    have = os.path.getsize(fname) if os.path.exists(fname) else 0
    if total is not None and have == total:
        print(f"[cached] {fname} ({have/1e9:.2f} GB) — complete, skipping", flush=True)
        return fname
    if total is not None and have > total:
        print(f"[warn] local file larger than remote; removing and restarting: {fname}", flush=True)
        os.remove(fname)
        have = 0
    def _p(msg):
        try:
            print(msg, flush=True)
        except UnicodeEncodeError:
            print(str(msg).encode("ascii", "replace").decode("ascii"), flush=True)
    _p(f"[download] {url}\n    -> {fname}\n    remote={total} local={have}")
    attempt = 0
    while attempt < retries:
        attempt += 1
        headers = {"Range": f"bytes={have}-"} if have > 0 else {}
        try:
            with requests.get(url, headers=headers, stream=True, timeout=120) as r:
                if have > 0 and r.status_code != 206:
                    print("[warn] server ignored Range; restarting from 0", flush=True)
                    have = 0
                    mode = "wb"
                else:
                    mode = "ab" if have > 0 else "wb"
                r.raise_for_status()
                t0 = time.time()
                last = t0
                with open(fname, mode) as f:
                    for chunk in r.iter_content(CHUNK):
                        if not chunk:
                            continue
                        f.write(chunk)
                        have += len(chunk)
                        now = time.time()
                        if now - last > 15:
                            pct = f"{100*have/total:5.1f}%" if total else "?"
                            rate = have / max(now - t0, 1e-9) / 1e6
                            print(f"    {pct} {have/1e9:.2f} GB  ({rate:.1f} MB/s)", flush=True)
                            last = now
            if total is None or have == total:
                print(f"[done] {fname} ({have/1e9:.2f} GB)", flush=True)
                return fname
            print(f"[warn] size mismatch after stream ({have} vs {total}); retrying", flush=True)
        except (requests.RequestException, OSError) as e:
            print(f"[retry {attempt}/{retries}] {e}", flush=True)
            time.sleep(min(60, 5 * attempt))
            have = os.path.getsize(fname) if os.path.exists(fname) else 0
    raise RuntimeError(f"failed to download {url}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("urls", nargs="+")
    ap.add_argument("--out", default="data/raw")
    a = ap.parse_args()
    for u in a.urls:
        download(u, a.out)

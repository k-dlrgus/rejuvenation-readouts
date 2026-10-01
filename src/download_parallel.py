"""Multi-connection, resumable downloader for throttled hosts (e.g., NCBI GEO).

Splits the file into N segments, downloads each with HTTP Range in its own thread
into <name>.segK files (resumable), then concatenates. Skips if the final file is
already complete (size == Content-Length). Never re-downloads a completed file.

Usage: python src/download_parallel.py URL --out DIR [--n 16]
"""
import argparse
import os
import sys
import threading
import time
import requests

CHUNK = 1 << 20


def remote_size(url):
    r = requests.head(url, allow_redirects=True, timeout=60)
    if r.status_code == 200 and "Content-Length" in r.headers:
        return int(r.headers["Content-Length"])
    r = requests.get(url, headers={"Range": "bytes=0-0"}, stream=True, timeout=60)
    return int(r.headers["Content-Range"].rsplit("/", 1)[1])


def seg_worker(url, path, start, end, progress, idx, stop):
    """Download bytes [start, end] into path, resuming from existing size."""
    while not stop.is_set():
        have = os.path.getsize(path) if os.path.exists(path) else 0
        want = end - start + 1
        if have >= want:
            progress[idx] = want
            return
        try:
            with requests.get(url, headers={"Range": f"bytes={start+have}-{end}"}, stream=True, timeout=120) as r:
                if r.status_code != 206:
                    raise RuntimeError(f"expected 206, got {r.status_code}")
                with open(path, "ab") as f:
                    for chunk in r.iter_content(CHUNK):
                        if stop.is_set():
                            return
                        if chunk:
                            f.write(chunk)
                            have += len(chunk)
                            progress[idx] = have
            if have >= want:
                return
        except Exception as e:  # noqa
            time.sleep(5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("url")
    ap.add_argument("--out", default="data/raw")
    ap.add_argument("--n", type=int, default=16)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    name = a.url.rsplit("/", 1)[-1]
    final = os.path.join(a.out, name)
    total = remote_size(a.url)
    if os.path.exists(final) and os.path.getsize(final) == total:
        print(f"[cached] {final} complete ({total/1e9:.2f} GB); skipping", flush=True)
        return
    print(f"[parallel download] {a.url}\n  -> {final}  total={total/1e9:.2f} GB  connections={a.n}", flush=True)
    seg = total // a.n
    bounds = [(i * seg, (total - 1) if i == a.n - 1 else (i + 1) * seg - 1) for i in range(a.n)]
    paths = [f"{final}.seg{i:02d}" for i in range(a.n)]
    progress = [0] * a.n
    stop = threading.Event()
    threads = [threading.Thread(target=seg_worker, args=(a.url, paths[i], s, e, progress, i, stop), daemon=True)
               for i, (s, e) in enumerate(bounds)]
    for t in threads:
        t.start()
    t0 = time.time()
    last_bytes = sum(progress)
    try:
        while any(t.is_alive() for t in threads):
            time.sleep(20)
            done = sum(progress)
            rate = (done - last_bytes) / 20 / 1e6
            last_bytes = done
            print(f"  {100*done/total:5.1f}%  {done/1e9:.2f}/{total/1e9:.2f} GB  {rate:.2f} MB/s  elapsed {int(time.time()-t0)}s", flush=True)
    except KeyboardInterrupt:
        stop.set()
        raise
    # verify and concatenate
    for i, (s, e) in enumerate(bounds):
        sz = os.path.getsize(paths[i])
        if sz != e - s + 1:
            raise RuntimeError(f"segment {i} incomplete: {sz} vs {e-s+1}")
    print("[concat] assembling segments ...", flush=True)
    with open(final + ".part", "wb") as fo:
        for p in paths:
            with open(p, "rb") as fi:
                while True:
                    b = fi.read(1 << 24)
                    if not b:
                        break
                    fo.write(b)
    os.replace(final + ".part", final)
    assert os.path.getsize(final) == total
    for p in paths:
        os.remove(p)
    print(f"[done] {final} ({total/1e9:.2f} GB)", flush=True)


if __name__ == "__main__":
    main()

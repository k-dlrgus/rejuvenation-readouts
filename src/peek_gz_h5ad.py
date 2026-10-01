"""Peek at the layout of a gzipped h5ad WITHOUT decompressing to disk, by handing h5py a gzip stream.
Reports top-level keys, X layout (dense vs sparse), shape, dtype, and gzip ISIZE (uncompressed size mod 2^32).
Backward seeks force a rewind, so this can be slow for metadata stored at the end of the file; we set a time budget."""
import gzip
import os
import struct
import sys
import time
import h5py


class SeekableGzip:
    """Seekable read-only view of a gzip stream for h5py's file-object driver.

    - Forward reads stream through the decompressor.
    - The first `head_cache` bytes (HDF5 superblock, root group, dataset headers) are cached in memory.
    - Once `tail_from` is set (e.g. the byte where the big dense dataset ends), every byte streamed at or beyond
      it is kept in memory, so all obs/var metadata seeks are served from RAM.
    - Backward seeks outside both caches rewind and re-decompress (slow but correct); counted in `rewinds`.
    """
    def __init__(self, path, head_cache=64 << 20):
        self.path = path
        self.f = gzip.open(path, "rb")
        self.pos = 0          # decompressor position
        self.vpos = 0         # position h5py believes it is at
        self.rewinds = 0
        self.bytes_read = 0
        self.head = bytearray()
        self.head_cache = head_cache
        self.tail_from = None
        self.tail = bytearray()
        self.tail_start = None

    # ---- internal streaming with cache maintenance ----
    def _stream(self, n):
        b = self.f.read(n)
        if not b:
            return b
        # extend each cache by exactly the not-yet-cached suffix of this block (robust to re-streaming after a rewind)
        if self.pos < self.head_cache:
            take = min(len(b), self.head_cache - self.pos)
            if self.pos <= len(self.head) <= self.pos + take:
                self.head += b[len(self.head) - self.pos:take]
        if self.tail_from is not None and self.pos + len(b) > self.tail_from:
            s = max(0, self.tail_from - self.pos)
            if self.tail_start is None:
                self.tail_start = self.pos + s
            tail_end = self.tail_start + len(self.tail)
            if self.pos + s <= tail_end <= self.pos + len(b):
                self.tail += b[tail_end - self.pos:]
        self.pos += len(b)
        self.bytes_read += len(b)
        return b

    def _position_decompressor(self, target):
        if target < self.pos:
            # caches stay valid (file content is immutable); the contiguity checks in _stream prevent double-appends
            self.f.close(); self.f = gzip.open(self.path, "rb"); self.pos = 0; self.rewinds += 1
        while self.pos < target:
            if not self._stream(min(1 << 24, target - self.pos)):
                break

    def prefetch_to_eof(self, block=1 << 24):
        """Stream forward to the end of the gzip member(s), filling the tail cache (call after tail_from is set)."""
        while self._stream(block):
            pass
        return self.pos

    def _cached_prefix(self, p, n):
        """Bytes [p, p+k) available from a cache, k <= n (possibly 0)."""
        if p < len(self.head):
            k = min(n, len(self.head) - p)
            return self.head[p:p + k]
        if self.tail_start is not None and self.tail_start <= p < self.tail_start + len(self.tail):
            o = p - self.tail_start
            k = min(n, len(self.tail) - o)
            return self.tail[o:o + k]
        return b""

    def read(self, n=-1):
        if n is None or n < 0:
            raise ValueError("unbounded read not supported")
        p = self.vpos
        out = bytearray(self._cached_prefix(p, n))
        rem = n - len(out)
        if rem > 0:
            # a read that starts inside a cache and runs past its end must NOT rewind: the cache end is exactly the
            # decompressor position, so the remainder streams forward. Only a genuinely uncached backward seek rewinds.
            self._position_decompressor(p + len(out))
            out += self._stream(rem)
        self.vpos += len(out)
        return bytes(out)

    def readinto(self, b):
        data = self.read(len(b))
        b[:len(data)] = data
        return len(data)

    def seek(self, off, whence=0):
        if whence == 1:
            off = self.vpos + off
        elif whence == 2:
            # HDF5 asks for the file size once (seek(0, 2); tell()). Uncompressed size is unknown without a full
            # pass, so report a large fake EOF; HDF5 only uses it as an upper bound for addresses.
            self.vpos = (1 << 41) + off
            return self.vpos
        self.vpos = off
        return self.vpos

    def tell(self):
        return self.vpos

    def seekable(self):
        return True

    def readable(self):
        return True

    def writable(self):
        return False

    def flush(self):
        pass

    def close(self):
        self.f.close()


def main(path):
    with open(path, "rb") as f:
        f.seek(-4, 2)
        isize = struct.unpack("<I", f.read(4))[0]
    print(f"gzip ISIZE (uncompressed size mod 2^32): {isize:,} bytes; compressed size {os.path.getsize(path)/1e9:.2f} GB")
    t0 = time.time()
    sg = SeekableGzip(path)
    h = h5py.File(sg, "r")
    print("top-level keys:", list(h.keys()), f"({time.time()-t0:.0f}s, streamed {sg.bytes_read/1e9:.2f} GB, rewinds={sg.rewinds})")
    X = h["X"]
    if isinstance(X, h5py.Group):
        print("X: sparse group", dict(X.attrs), "data dtype", X["data"].dtype, "nnz", X["data"].shape)
    else:
        print("X: DENSE dataset shape", X.shape, "dtype", X.dtype, "chunks", X.chunks, "compression", X.compression,
              f"-> {X.shape[0]*X.shape[1]*X.dtype.itemsize/1e9:.1f} GB uncompressed")
    for k in ("obs", "var"):
        if k in h:
            g = h[k]
            print(f"{k} keys:", list(g.keys())[:40], f"({time.time()-t0:.0f}s, streamed {sg.bytes_read/1e9:.2f} GB, rewinds={sg.rewinds})")
            idx = g.attrs.get("_index")
            if idx is not None:
                print(f"{k} n =", g[idx.decode() if isinstance(idx, bytes) else idx].shape)
    h.close()


if __name__ == "__main__":
    main(sys.argv[1])

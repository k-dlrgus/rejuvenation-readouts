"""Peek an HDF5 object header stored at the start of a ranged download."""
import struct
from pathlib import Path

BASE = 592121592
data = Path("results/survey/opened/sengstack_var.bin").read_bytes()
print("ver", data[0], "nmsg", struct.unpack_from("<H", data, 2)[0], "hdr", struct.unpack_from("<I", data, 8)[0])
nmsg = struct.unpack_from("<H", data, 2)[0]
hdr = struct.unpack_from("<I", data, 8)[0]
pos = 16
end = 16 + hdr
i = 0
while pos + 8 <= end and i < nmsg:
    mtype, msize, flags = struct.unpack_from("<HHB", data, pos)
    body = data[pos + 8 : pos + 8 + msize]
    print("MSG", mtype, "size", msize, "flags", flags, "body_hex", body[:48].hex())
    if mtype == 17 and len(body) >= 16:
        btree, heap = struct.unpack_from("<QQ", body, 0)
        print(" btree", btree, "heap", heap, "in_chunk", 0 <= btree - BASE < len(data), 0 <= heap - BASE < len(data))
    pos += 8 + msize
    i += 1
print("parsed", i)

import struct
from pathlib import Path

BASE = 592120976
data = Path("results/survey/opened/sengstack_obs.bin").read_bytes()
OFF = 8


def u64(off):
    return struct.unpack_from("<Q", data, off)[0]


def heap_string(heap_abs, name_off):
    heap = heap_abs - BASE
    if data[heap : heap + 4] != b"HEAP":
        raise SystemExit("bad heap")
    data_addr = u64(heap + 8 + 8 + 8)
    start = (data_addr - BASE) + name_off
    end = data.index(b"\x00", start)
    return data[start:end].decode()


btree = 592121016 - BASE
heap = 592121560
nent = struct.unpack_from("<H", data, btree + 6)[0]
print("nent", nent, "type", data[btree + 4], "level", data[btree + 5])
pos = btree + 8 + OFF + OFF
keys = []
children = []
for i in range(nent + 1):
    keys.append(u64(pos))
    pos += OFF
    if i < nent:
        children.append(u64(pos))
        pos += OFF
print("children", children)
for child in children:
    rel = child - BASE
    print("child", child, "in", 0 <= rel < len(data))
    if not (0 <= rel < len(data)):
        continue
    print("sig", data[rel : rel + 4])
    nsym = struct.unpack_from("<H", data, rel + 6)[0]
    epos = rel + 8
    for _ in range(nsym):
        name_off = u64(epos)
        obj = u64(epos + OFF)
        name = heap_string(heap, name_off)
        print(name, obj)
        epos += OFF * 2 + 24

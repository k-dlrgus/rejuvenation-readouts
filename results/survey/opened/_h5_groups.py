"""List root-group link names from a local HDF5 prefix, following small range fetches.

Used only to locate gene-name datasets without downloading multi-GB matrices.
"""
import struct
import sys
from pathlib import Path

path = Path(sys.argv[1])
data = path.read_bytes()
OFF = 8


def u64(off):
    return struct.unpack_from("<Q", data, off)[0]


def need(n):
    if n > len(data):
        raise SystemExit(f"need {n} have {len(data)}")


def heap_string(heap_addr, name_off):
    # local heap: HEAP, version u8, reserved 3, data_seg_size u64, free_off u64, data_addr u64
    sig = data[heap_addr : heap_addr + 4]
    if sig != b"HEAP":
        raise SystemExit(f"bad heap {sig} at {heap_addr}")
    data_addr = u64(heap_addr + 8 + 8 + 8)
    start = data_addr + name_off
    end = data.index(b"\x00", start)
    return data[start:end].decode("utf-8", errors="replace")


def parse_group(obj_addr, prefix, depth=0):
    if depth > 4:
        return
    need(obj_addr + 16)
    ver = data[obj_addr]
    if ver != 1:
        print(f"skip obj ver {ver} at {obj_addr} {prefix}")
        return
    nmsg = struct.unpack_from("<H", data, obj_addr + 2)[0]
    hdr_size = struct.unpack_from("<I", data, obj_addr + 8)[0]
    pos = obj_addr + 16
    end = pos + hdr_size
    need(end)
    while pos + 8 <= end:
        mtype, msize, mflags = struct.unpack_from("<HHB", data, pos)
        pos += 8
        body = pos
        pos += msize
        if mtype == 17:  # symbol table
            btree = struct.unpack_from("<Q", data, body)[0]
            heap = struct.unpack_from("<Q", data, body + 8)[0]
            walk_btree(btree, heap, prefix, depth)
        elif mtype == 6:  # link info? skip
            pass


def walk_btree(btree, heap, prefix, depth):
    need(btree + 24)
    if data[btree : btree + 4] != b"TREE":
        print("bad tree", btree)
        return
    nent = struct.unpack_from("<H", data, btree + 6)[0]
    pos = btree + 8 + OFF + OFF  # skip sibling addrs
    # nent+1 keys, nent children
    keys = []
    children = []
    for i in range(nent + 1):
        keys.append(u64(pos))
        pos += OFF
        if i < nent:
            children.append(u64(pos))
            pos += OFF
    for child, key in zip(children, keys[1:]):
        # symbol table node SNOD
        need(child + 8)
        if data[child : child + 4] != b"SNOD":
            print("not snod", child, data[child : child + 4])
            continue
        nsym = struct.unpack_from("<H", data, child + 6)[0]
        epos = child + 8
        for _ in range(nsym):
            name_off = u64(epos)
            obj = u64(epos + OFF)
            epos += OFF * 2 + 4 + 4 + 16
            name = heap_string(heap, name_off)
            full = f"{prefix}/{name}" if prefix else name
            print("GROUP", full, "obj", obj)


root = u64(64)
print("root", root)
parse_group(root, "")

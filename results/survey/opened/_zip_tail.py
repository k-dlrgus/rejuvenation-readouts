"""Parse a ZIP end-of-central-directory from a downloaded tail."""
import struct
from pathlib import Path

data = Path("results/survey/opened/southard/zip_tail.bin").read_bytes()
# EOCD
sig = b"PK\x05\x06"
idx = data.rfind(sig)
print("EOCD_OFFSET_IN_TAIL", idx, "TAIL", len(data))
if idx < 0:
    z64 = data.rfind(b"PK\x06\x06")
    loc = data.rfind(b"PK\x06\x07")
    print("ZIP64_EOCD", z64, "ZIP64_LOC", loc)
    # print nearby ascii names
    text = data.decode("latin1", errors="ignore")
    hits = [line for line in text.split("PK") if "feature" in line.lower() or "gene" in line.lower()]
    print("HITS", len(hits))
    for h in hits[:20]:
        print(repr(h[:120]))
else:
    eocd = data[idx : idx + 22]
    (
        _sig,
        disk,
        cd_disk,
        n_disk,
        n_total,
        cd_size,
        cd_offset,
        comment_len,
    ) = struct.unpack("<IHHHHIIH", eocd)
    print("N", n_total, "CD_SIZE", cd_size, "CD_OFFSET", cd_offset, "COMMENT", comment_len)
    loc_off = idx - 20
    loc = data[loc_off : loc_off + 20]
    print("LOC_SIG", loc[:4])
    if loc[:4] == b"PK\x06\x07":
        _s, disk_z, z64_off, ndisk = struct.unpack("<IIQI", loc)
        print("ZIP64_EOCD_ABS", z64_off)
        file_size = 21586390757
        tail_start = file_size - len(data)
        rel = z64_off - tail_start
        print("ZIP64_EOCD_IN_TAIL", rel)
        if 0 <= rel < len(data):
            zsig = data[int(rel) : int(rel) + 4]
            print("ZSIG", zsig)
            # zip64 EOCD fixed part
            size, ver, ver_need, disk, cd_disk, n1, n2, cd_size64, cd_off64 = struct.unpack(
                "<QHHIIQQQQ", data[int(rel) + 4 : int(rel) + 4 + 52]
            )
            print("CD_SIZE64", cd_size64, "CD_OFF64", cd_off64, "N64", n2)
            cd_rel = cd_off64 - tail_start
            print("CD_IN_TAIL", cd_rel)
            if 0 <= cd_rel and cd_rel + cd_size64 <= len(data):
                cd = data[int(cd_rel) : int(cd_rel + cd_size64)]
                pos = 0
                names = []
                while pos + 46 <= len(cd):
                    if cd[pos : pos + 4] != b"PK\x01\x02":
                        print("BAD_CD", pos)
                        break
                    name_len, extra_len, cmt_len = struct.unpack("<HHH", cd[pos + 28 : pos + 34])
                    local_off = struct.unpack("<I", cd[pos + 42 : pos + 46])[0]
                    comp_size = struct.unpack("<I", cd[pos + 20 : pos + 24])[0]
                    name = cd[pos + 46 : pos + 46 + name_len].decode("utf-8", errors="replace")
                    extra = cd[pos + 46 + name_len : pos + 46 + name_len + extra_len]
                    # zip64 extra 0x0001
                    if local_off == 0xFFFFFFFF or comp_size == 0xFFFFFFFF:
                        # parse extra
                        epos = 0
                        while epos + 4 <= len(extra):
                            hid, hsz = struct.unpack("<HH", extra[epos : epos + 4])
                            blob = extra[epos + 4 : epos + 4 + hsz]
                            epos += 4 + hsz
                            if hid == 1:
                                vals = []
                                bpos = 0
                                if comp_size == 0xFFFFFFFF and bpos + 8 <= len(blob):
                                    comp_size = struct.unpack("<Q", blob[bpos : bpos + 8])[0]
                                    bpos += 8
                                if local_off == 0xFFFFFFFF and bpos + 8 <= len(blob):
                                    # uncompressed then offset; skip uncomp if present
                                    pass
                                # standard order: uncomp, comp, offset, disk
                                bpos = 0
                                if struct.unpack("<I", cd[pos + 24 : pos + 28])[0] == 0xFFFFFFFF and bpos + 8 <= len(blob):
                                    bpos += 8
                                if struct.unpack("<I", cd[pos + 20 : pos + 24])[0] == 0xFFFFFFFF and bpos + 8 <= len(blob):
                                    comp_size = struct.unpack("<Q", blob[bpos : bpos + 8])[0]
                                    bpos += 8
                                if struct.unpack("<I", cd[pos + 42 : pos + 46])[0] == 0xFFFFFFFF and bpos + 8 <= len(blob):
                                    local_off = struct.unpack("<Q", blob[bpos : bpos + 8])[0]
                    names.append((name, int(comp_size), int(local_off)))
                    pos += 46 + name_len + extra_len + cmt_len
                print("NFILES", len(names))
                for name, csz, off in names:
                    if "feature" in name.lower() or name.lower().endswith("genes.tsv.gz") or "gene" in name.lower():
                        print("HIT", csz, off, name)
                print("SAMPLE")
                for name, csz, off in names[:8]:
                    print(csz, off, name)

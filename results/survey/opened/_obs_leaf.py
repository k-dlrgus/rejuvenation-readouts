import struct
from pathlib import Path

obs = Path("results/survey/opened/sengstack_obs.bin").read_bytes()
OBS_BASE = 592120976
leaf = Path("results/survey/opened/sengstack_obs_leaf.bin").read_bytes()
OFF = 8


def heap_string(name_off):
    heap = 592121560 - OBS_BASE
    data_addr = struct.unpack_from("<Q", obs, heap + 24)[0]
    start = (data_addr - OBS_BASE) + name_off
    end = obs.index(b"\x00", start)
    return obs[start:end].decode()


nent = struct.unpack_from("<H", leaf, 6)[0]
pos = 8 + OFF + OFF
keys = []
children = []
for i in range(nent + 1):
    keys.append(struct.unpack_from("<Q", leaf, pos)[0])
    pos += OFF
    if i < nent:
        children.append(struct.unpack_from("<Q", leaf, pos)[0])
        pos += OFF
print("keys", [heap_string(k) if k != 0 else "" for k in keys])
print("snod_addrs")
for c in children:
    print(c)

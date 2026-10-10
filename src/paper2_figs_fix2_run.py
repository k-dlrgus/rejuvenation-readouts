"""PAPER2FIGS-fix2: redraw Fig 3 only into results/paper2_figs_fix/run4/ (plan: run4/PLAN.flag).
Hard rule: run4/Fig3_data.csv equals run3/Fig3_data.csv after CRLF->LF."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paper2_figs as pf  # noqa: E402
from config import RESULTS  # noqa: E402

OUT = RESULTS / "paper2_figs_fix" / "run4"
OLD = RESULTS / "paper2_figs_fix" / "run3" / "Fig3_data.csv"
pf.FIX, pf.LABEL_FALLBACK, pf.OUT = True, True, OUT
pf.FIG4_MANUAL = dict(CRI_ZFX=(-22.0, 21.0), CRA_EZH2=(-10.0, 21.0), CRA_E2F3=(-22.0, 3.2), CRI_STAT3=(-20.0, 0.0))
pf.FIG3_MANUAL = dict(CRA_E2F3=(0.22, -3.6))
OUT.mkdir(parents=True, exist_ok=True)
pf.fig3()
a = (OUT / "Fig3_data.csv").read_bytes()
b = OLD.read_bytes()
la, lb = a.replace(b"\r\n", b"\n"), b.replace(b"\r\n", b"\n")
ok = la == lb and b"\r" not in la and b"\r" not in lb
print("Fig3_data.csv LF-equal to run3:", ok)
for k, v in pf.LABEL_CHECK.items():
    print(k, v)
sys.exit(0 if ok else 1)

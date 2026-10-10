"""PAPER2FIGS-fix4: redraw Fig 3 and Fig S2 into results/paper2_figs_fix/run6/ (plan: run6/PLAN.flag).
Hard rule: *_data.csv equal to the previous version after CRLF->LF."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paper2_figs as pf  # noqa: E402
from config import RESULTS  # noqa: E402

B = RESULTS / "paper2_figs_fix"
OUT = B / "run6"
pf.FIX, pf.LABEL_FALLBACK, pf.OUT = True, True, OUT
pf.FIG4_MANUAL = dict(CRI_ZFX=(-22.0, 21.0), CRA_EZH2=(-10.0, 21.0), CRA_E2F3=(-22.0, 3.2), CRI_STAT3=(-20.0, 0.0))
pf.FIG3_MANUAL = dict(CRA_E2F3=(0.22, -3.6), CRA_EZH2=(0.16, 1.4), CRI_ZFX=(0.008, -20.4))
OUT.mkdir(parents=True, exist_ok=True)
pf.fig3()
pf.figS2()
ok = True
for name, old in (("Fig3_data.csv", B / "run5"), ("FigS2_data.csv", B / "run3")):
    a, b = (OUT / name).read_bytes(), (old / name).read_bytes()
    la, lb = a.replace(b"\r\n", b"\n"), b.replace(b"\r\n", b"\n")
    r = la == lb and b"\r" not in la and b"\r" not in lb
    print(name, "LF-equal to", old.name, r)
    ok &= r
sys.exit(0 if ok else 1)

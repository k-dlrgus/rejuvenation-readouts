"""Write results/paper_figs/environment_versions.txt: Python and every third-party package imported under src/."""
from __future__ import annotations

import ast
import importlib.metadata as md
import platform
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
OUT = ROOT / "results" / "paper_figs" / "environment_versions.txt"
REQUESTED = ("numpy", "scipy", "pandas", "scikit-learn", "anndata", "h5py", "matplotlib")


def imported_top_modules():
    found = {}
    for p in sorted(SRC.glob("*.py")):
        tree = ast.parse(p.read_text(encoding="utf-8-sig"), filename=str(p))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                names = [node.module]
            else:
                continue
            for n in names:
                found.setdefault(n.split(".")[0], set()).add(p.name)
    return found


def main():
    local = {p.stem for p in SRC.glob("*.py")}
    stdlib = set(sys.stdlib_module_names)
    mods = imported_top_modules()
    third = {m: f for m, f in mods.items() if m not in local and m not in stdlib and m != "__future__"}
    dist_map = md.packages_distributions()
    lines = [
        f"python {sys.version.split()[0]} ({platform.python_implementation()}, {sys.version})",
        f"executable {sys.executable}",
        f"platform {platform.platform()}",
        "",
        "# requested",
    ]
    for d in REQUESTED:
        try:
            lines.append(f"{d}=={md.version(d)}")
        except md.PackageNotFoundError:
            lines.append(f"{d}: NOT INSTALLED")
    lines += ["", "# every third-party top-level module imported under src/ (module -> distribution==version; files)"]
    for m in sorted(third):
        dists = dist_map.get(m, [])
        if dists:
            ver = ", ".join(f"{d}=={md.version(d)}" for d in sorted(set(dists)))
        else:
            ver = "NOT INSTALLED in this environment"
        files = sorted(third[m])
        shown = ", ".join(files[:6]) + (f", ... (+{len(files) - 6})" if len(files) > 6 else "")
        lines.append(f"{m} -> {ver}  [{shown}]")
    lines += ["", "# stdlib modules imported under src/ (versioned with Python)",
              ", ".join(sorted(m for m in mods if m in stdlib))]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()

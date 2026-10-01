"""Write gene symbols from the Hs27 mean-population AnnData that was opened."""
import h5py
from pathlib import Path

src = Path("results/survey/opened/southard/fibroblast_CRISPRa_mean_pop.h5ad")
out = Path("results/survey/gene_lists/southard_hs27_mean_pop_genes.txt")
out.parent.mkdir(parents=True, exist_ok=True)
with h5py.File(src, "r") as f:
    names = [x.decode() if isinstance(x, bytes) else str(x) for x in f["var/gene_name"][:]]
    targets = [x.decode() if isinstance(x, bytes) else str(x) for x in f["obs/target_gene/categories"][:]]
    guides = f["obs/guide_identity"].shape[0]
print("n_genes", len(names), "unique", len(set(names)))
print("n_target_categories", len(targets))
print("n_obs", guides)
print("targets_head", targets[:8])
out.write_text("\n".join(names) + "\n", encoding="utf-8")
Path("results/survey/opened/southard/mean_pop_targets.txt").write_text(
    "\n".join(targets) + "\n", encoding="utf-8"
)

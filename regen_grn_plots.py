"""
Regenerate GRN plots with the updated node filter (abs(r)>0.1 in both conditions).
"""
import sys, json
from pathlib import Path

sys.path.insert(0, "/vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbScope_v1")

H5AD = "/vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbScope_v1/results_replogle_weissman_k562_essential_mini20/final_adata.h5ad"
RESULTS = "/vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbScope_v1/results_replogle_weissman_k562_essential_mini20"

import anndata
print("Loading adata ...")
adata = anndata.read_h5ad(H5AD)
print(f"Loaded: {adata.shape}")
print("obs cols:", list(adata.obs.columns)[:10])
print("perturbation counts:\n", adata.obs["perturbation"].value_counts().head())

with open(f"{RESULTS}/checkpoint.json") as f:
    ckpt = json.load(f)
top_perts = ckpt.get("top_perturbations", None)
print("top_perturbations:", top_perts)

from perturbscope.genenet import run_gene_network
print("Re-running GRN plots with node filter abs(r)>0.1 in both conditions ...")
run_gene_network(
    adata,
    output_dir=RESULTS,
    perturbations=top_perts,
    n_top_genes=50,
    n_gene_clusters=5,
    corr_threshold=0.5,
)
print("Done. GRN plots regenerated.")

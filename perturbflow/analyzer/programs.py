# =============================================================================
# perturbflow/analyzer/programs.py
#
# Gene-program and pathway activity inference module.
#
# Scores each perturbation group against a set of curated gene programs
# (e.g. cell cycle, interferon response, stress).  The score is the mean
# log-normalised expression of the program genes that are present in the
# dataset, giving a quick read-out of which biological programs are activated
# or suppressed by each perturbation.
#
# Users can pass any custom gene_sets dict to infer_programs() to override
# the built-in defaults.
#
# Outputs written to adata.uns["program_scores"] and returned as a DataFrame
# with columns: perturbation, program, program_score, n_genes_used.
# =============================================================================

from __future__ import annotations

import numpy as np
import pandas as pd


def _gene_index(adata):
    # Build a name → integer-index lookup dict for fast gene access.
    return {g: i for i, g in enumerate(adata.var_names.astype(str))}


def infer_programs(adata, gene_sets=None):
    # Score every perturbation group against a dictionary of gene programs.
    #
    # Args:
    #   gene_sets -- dict mapping program name to list of gene symbols.
    #                If None, three built-in programs are used:
    #                  cell_cycle  : MKI67, TOP2A, PCNA, CDK1
    #                  interferon  : STAT1, ISG15, IFI6, IFIT1
    #                  stress      : JUN, FOS, ATF3, HSP90AA1
    #
    # Returns (adata, program_df).
    if gene_sets is None:
        gene_sets = {
            "cell_cycle": ["MKI67", "TOP2A", "PCNA", "CDK1"],
            "interferon": ["STAT1", "ISG15", "IFI6", "IFIT1"],
            "stress": ["JUN", "FOS", "ATF3", "HSP90AA1"],
        }

    gidx = _gene_index(adata)
    rows = []
    for perturb in sorted(adata.obs["perturbation"].astype(str).unique()):
        pm = adata.obs["perturbation"].values == perturb
        if pm.sum() == 0:
            continue
        for pname, genes in gene_sets.items():
            # Restrict to genes actually present in the dataset.
            idx = [gidx[g] for g in genes if g in gidx]
            if not idx:
                score = np.nan  # none of the program genes are in the data
            else:
                x = adata.X[pm][:, idx]
                score = float((x.A if hasattr(x, "A") else x).mean())
            rows.append(
                {
                    "perturbation": perturb,
                    "program": pname,
                    "program_score": score,
                    "n_genes_used": len(idx),
                }
            )

    program_df = pd.DataFrame(rows)
    adata.uns["program_scores"] = program_df.to_dict(orient="records")
    return adata, program_df

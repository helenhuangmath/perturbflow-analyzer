# =============================================================================
# perturbscope/interaction.py
#
# Combinatorial perturbation interaction analysis.
#
# For perturbations targeting multiple genes simultaneously (e.g. "GENE1+GENE2"),
# this module tests whether the combined effect deviates from the expected
# additive null model:
#
#   interaction_score = observed_effect - expected_additive_effect
#
# where expected_additive_effect is the mean transcriptional score of the
# individual single-gene perturbations (from state.py).
#
# Classification thresholds (±0.15 from zero):
#   synergistic  -- score > +0.15  (genes amplify each other)
#   antagonistic -- score < -0.15  (genes dampen each other)
#   additive     -- score within ±0.15 (independent effects)
#
# NOTE: This module requires that compute_effect_decomposition() has already
# been run so that adata.uns["effect_scores"] is populated.
#
# Returns (adata, interaction_df) where interaction_df has one row per
# combinatorial perturbation.
# =============================================================================

from __future__ import annotations

import numpy as np
import pandas as pd


def _parse_combo(name: str) -> list[str]:
    # Split a perturbation label into its component gene names.
    # Supports "+" and "," delimiters (e.g. "GENE1+GENE2" or "GENE1,GENE2").
    return [p.strip() for p in str(name).replace(",", "+").split("+") if p.strip()]


def analyze_interactions(adata):
    # Identify combinatorial perturbations and classify their interaction type.
    #
    # Iterates over all perturbation labels, skips single-gene and control
    # entries, then compares the observed transcriptional effect against the
    # mean of the individual single-gene effects.
    effects = adata.uns.get("effect_scores", {})
    rows = []

    for perturb in sorted(adata.obs["perturbation"].astype(str).unique()):
        genes = _parse_combo(perturb)
        if len(genes) <= 1:
            continue  # not a combinatorial perturbation
        if perturb not in effects:
            continue  # no effect score was computed (too few cells)

        observed = float(effects[perturb].get("transcriptional", 0.0))
        # Collect the transcriptional scores for each component single gene.
        singles = [float(effects[g].get("transcriptional", np.nan)) for g in genes if g in effects]
        if not singles:
            continue  # no single-gene reference data available

        # Additive null: average effect of the individual components.
        expected = float(np.nanmean(singles))
        score = observed - expected

        if score > 0.15:
            cls = "synergistic"
        elif score < -0.15:
            cls = "antagonistic"
        else:
            cls = "additive"

        rows.append(
            {
                "perturbation": perturb,
                "components": ";".join(genes),
                "observed_effect": observed,
                "expected_additive_effect": expected,
                "interaction_score": score,
                "interaction_class": cls,
            }
        )

    df = pd.DataFrame(rows)
    adata.uns["interaction_results"] = df.to_dict(orient="records")
    return adata, df

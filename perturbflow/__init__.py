"""PerturbFlow: open infrastructure for perturbation biology.

AnnData-native, scverse-interoperable tooling for standardized analysis,
mechanistic (rewiring-aware) interpretation, AI-ready outputs, and community
benchmarks across Perturb-seq, CRISPR screens, and single-cell multi-omics.
"""

from __future__ import annotations

from .workflow import prepare_h5ad, run_analysis
from .ai import write_agent_handoff

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "prepare_h5ad",
    "run_analysis",
    "write_agent_handoff",
]

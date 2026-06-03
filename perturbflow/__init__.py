"""PerturbFlow: user-friendly workflows for Perturb-seq analysis."""

from __future__ import annotations

from .ai import write_agent_handoff
from .workflow import prepare_h5ad, run_analysis

__version__ = "0.1.0"

__all__ = ["__version__", "prepare_h5ad", "run_analysis", "write_agent_handoff"]

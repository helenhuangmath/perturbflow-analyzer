"""Public workflow helpers for PerturbFlow."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable


def prepare_h5ad(
    input_path: str | Path,
    output_path: str | Path,
    perturbation_col: str,
    control_labels: str | None = None,
    cell_state_col: str | None = None,
) -> Path:
    """Standardize a user AnnData file for the PerturbFlow pipeline."""
    from perturbflow.analyzer.process_data import prepare_h5ad as _prepare_h5ad

    return _prepare_h5ad(
        input_path=input_path,
        output_path=output_path,
        perturbation_col=perturbation_col,
        control_labels=control_labels,
        cell_state_col=cell_state_col,
    )


def run_analysis(
    input_path: str | Path,
    output_dir: str | Path,
    config_path: str | Path | None = None,
    steps: Iterable[str] | None = None,
    perturbation_col: str | None = None,
    resume: bool = True,
    force_steps: Iterable[str] | None = None,
    clear_from: str | None = None,
):
    """Run the end-to-end PerturbFlow analysis.

    This is a stable public wrapper around the analyzer engine.
    """
    from perturbflow.analyzer.config import PipelineConfig
    from perturbflow.analyzer.pipeline import run_pipeline

    config = PipelineConfig.from_json(config_path) if config_path else PipelineConfig()
    return run_pipeline(
        input_path=str(input_path),
        output_dir=str(output_dir),
        config=config,
        steps=steps,
        perturbation_col=perturbation_col,
        resume=resume,
        force_steps=force_steps,
        clear_from=clear_from,
    )

"""Programmatic API for PerturbFlow workflows."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

from perturbflow.analyzer.config import PipelineConfig
from perturbflow.workflow import prepare_h5ad, run_analysis


def _coerce_config(config: PipelineConfig | Mapping[str, Any] | str | Path | None) -> PipelineConfig | None:
    if config is None:
        return None
    if isinstance(config, PipelineConfig):
        return config
    if isinstance(config, (str, Path)):
        return PipelineConfig.from_json(config)
    if isinstance(config, Mapping):
        valid = {
            key: value
            for key, value in config.items()
            if key in PipelineConfig.__dataclass_fields__
        }
        return PipelineConfig(**valid)
    raise TypeError(
        "config must be a PipelineConfig, mapping, path, or None; "
        f"got {type(config).__name__}"
    )


@dataclass
class PerturbFlowAPI:
    """Small Python interface for using PerturbFlow from other programs.

    Parameters set on the API instance are used as defaults for each call.
    Per-call arguments override those defaults.
    """

    config: PipelineConfig | Mapping[str, Any] | str | Path | None = None
    perturbation_col: str | None = None
    resume: bool = True

    def pipeline_config(self) -> PipelineConfig:
        """Return a concrete pipeline configuration."""
        return _coerce_config(self.config) or PipelineConfig()

    def list_steps(self) -> list[str]:
        """Return the default analysis steps for this API configuration."""
        return list(self.pipeline_config().default_steps)

    def prepare(
        self,
        input_path: str | Path,
        output_path: str | Path,
        perturbation_col: str | None = None,
        control_labels: str | None = None,
        cell_state_col: str | None = None,
    ) -> Path:
        """Standardize an AnnData file for PerturbFlow."""
        pert_col = perturbation_col or self.perturbation_col
        if not pert_col:
            raise ValueError("perturbation_col is required for prepare().")
        return prepare_h5ad(
            input_path=input_path,
            output_path=output_path,
            perturbation_col=pert_col,
            control_labels=control_labels,
            cell_state_col=cell_state_col,
        )

    def analyze(
        self,
        input_path: str | Path,
        output_dir: str | Path,
        steps: Iterable[str] | None = None,
        perturbation_col: str | None = None,
        resume: bool | None = None,
        force_steps: Iterable[str] | None = None,
        clear_from: str | None = None,
        config: PipelineConfig | Mapping[str, Any] | str | Path | None = None,
    ):
        """Run the analyzer workflow and return the final AnnData object."""
        cfg = _coerce_config(config) or self.pipeline_config()
        return run_analysis(
            input_path=input_path,
            output_dir=output_dir,
            config_path=None,
            steps=steps,
            perturbation_col=perturbation_col or self.perturbation_col,
            resume=self.resume if resume is None else resume,
            force_steps=force_steps,
            clear_from=clear_from,
            config=cfg,
        )

    def run(
        self,
        input_path: str | Path,
        output_dir: str | Path,
        prepared_path: str | Path | None = None,
        perturbation_col: str | None = None,
        control_labels: str | None = None,
        cell_state_col: str | None = None,
        steps: Iterable[str] | None = None,
        resume: bool | None = None,
        force_steps: Iterable[str] | None = None,
        clear_from: str | None = None,
    ):
        """Optionally prepare input data, then run the analyzer workflow.

        If ``prepared_path`` is provided, ``input_path`` is first standardized
        and the analysis runs on that prepared file. If it is omitted, analysis
        runs directly on ``input_path``.
        """
        analysis_input = Path(input_path)
        if prepared_path is not None:
            analysis_input = self.prepare(
                input_path=input_path,
                output_path=prepared_path,
                perturbation_col=perturbation_col,
                control_labels=control_labels,
                cell_state_col=cell_state_col,
            )
        return self.analyze(
            input_path=analysis_input,
            output_dir=output_dir,
            steps=steps,
            perturbation_col=perturbation_col,
            resume=resume,
            force_steps=force_steps,
            clear_from=clear_from,
        )


__all__ = ["PerturbFlowAPI"]

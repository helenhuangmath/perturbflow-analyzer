"""PerturbFlow Analyzer: state-aware Perturb-seq analysis engine."""

__version__ = "1.0.0"

__all__ = ["PipelineConfig", "run_pipeline"]


def __getattr__(name):
    if name == "PipelineConfig":
        from .config import PipelineConfig

        return PipelineConfig
    if name == "run_pipeline":
        from .pipeline import run_pipeline

        return run_pipeline
    raise AttributeError(name)

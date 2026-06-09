"""Command-line interface for the PerturbFlow ecosystem."""

from __future__ import annotations

import argparse

from .ai import write_agent_handoff
from .workflow import prepare_h5ad, run_analysis

_DEFAULT_STEPS = [
    "qc",
    "preprocess",
    "eda",
    "score",
    "effects",
    "trajectory",
    "programs",
    "interaction",
    "state_enrich",
    "deg",
    "genenet",
    "tf_genenet",
    "cscore",
    "regulatory",
    "report",
    "bundle",
]


def _parse_steps(text: str | None):
    if not text:
        return None
    return [step.strip() for step in text.split(",") if step.strip()]


def cmd_prepare(args) -> None:
    out = prepare_h5ad(
        input_path=args.input,
        output_path=args.output,
        perturbation_col=args.perturbation_col,
        control_labels=args.control_labels,
        cell_state_col=args.cell_state_col,
    )
    print(f"Wrote prepared AnnData: {out}")


def cmd_run(args) -> None:
    run_analysis(
        input_path=args.input,
        output_dir=args.output,
        config_path=args.config,
        steps=_parse_steps(args.steps),
        perturbation_col=args.perturbation_col,
        resume=args.resume,
        force_steps=_parse_steps(args.force_steps),
        clear_from=args.clear_from,
    )


def cmd_interpret(args) -> None:
    files = write_agent_handoff(
        results_dir=args.results,
        output_dir=args.output,
        project_name=args.project_name,
        max_rows=args.max_rows,
    )
    print("Wrote PerturbFlow agent interpretation files:")
    for label, path in files.items():
        print(f"  {label}: {path}")


def cmd_list_steps(_args) -> None:
    print("Available pipeline steps (default order):")
    for index, step in enumerate(_DEFAULT_STEPS, 1):
        print(f"  {index:2}. {step}")


def cmd_predict(_args) -> None:
    raise SystemExit(
        "PerturbFlow Predictor is reserved for future perturbation response "
        "prediction features and is not implemented yet."
    )


def cmd_benchmark(_args) -> None:
    raise SystemExit(
        "PerturbFlow Benchmark is reserved for community evaluation tooling "
        "(baseline-calibrated, distance-aware, and rewiring-aware metrics) and "
        "is not implemented yet."
    )


def _add_run_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--input", required=True, help="Input .h5ad file")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--config", default=None, help="Optional JSON config")
    parser.add_argument("--steps", default=None, help="Comma-separated steps, e.g. qc,preprocess,deg,report,bundle")
    parser.add_argument("--perturbation-col", default=None, help="obs column with perturbation labels")
    parser.add_argument("--resume", action="store_true", default=True, help="Resume from checkpoint.json when available")
    parser.add_argument("--no-resume", action="store_false", dest="resume", help="Ignore checkpoint.json and rerun requested steps")
    parser.add_argument("--force-steps", default=None, help="Comma-separated steps to rerun even if checkpointed")
    parser.add_argument("--clear-from", default=None, help="Invalidate checkpoint from this step onward")
    parser.set_defaults(func=cmd_run)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="perturbflow",
        description="Open infrastructure for perturbation biology: analysis, reporting, agent interpretation, and benchmarking.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    prep = sub.add_parser("prepare", help="Standardize a user .h5ad file")
    prep.add_argument("--input", required=True, help="Input .h5ad file")
    prep.add_argument("--output", required=True, help="Output standardized .h5ad file")
    prep.add_argument("--perturbation-col", required=True, help="obs column containing perturbation labels")
    prep.add_argument("--control-labels", default=None, help="Comma-separated labels to map to control")
    prep.add_argument("--cell-state-col", default=None, help="Optional obs column to copy to cell_state")
    prep.set_defaults(func=cmd_prepare)

    analyzer = sub.add_parser("analyzer", help="Run the analyzer workflow")
    _add_run_args(analyzer)

    analyze = sub.add_parser("analyze", help="Alias for analyzer")
    _add_run_args(analyze)

    run = sub.add_parser("run", help="Alias for analyzer")
    _add_run_args(run)

    predict = sub.add_parser("predict", help="Reserved for future perturbation response prediction")
    predict.set_defaults(func=cmd_predict)

    benchmark = sub.add_parser("benchmark", help="Reserved for community evaluation/benchmarking (Aim 3)")
    benchmark.set_defaults(func=cmd_benchmark)

    interpret = sub.add_parser("interpret", help="Create agent interpretation files")
    interpret.add_argument("--results", required=True, help="Completed PerturbFlow result directory")
    interpret.add_argument("--output", default=None, help="Output directory for agent interpretation files")
    interpret.add_argument("--project-name", default="PerturbFlow run", help="Human-readable analysis name")
    interpret.add_argument("--max-rows", type=int, default=10, help="Rows to include per summary table")
    interpret.set_defaults(func=cmd_interpret)

    list_steps = sub.add_parser("list-steps", help="Print available pipeline steps")
    list_steps.set_defaults(func=cmd_list_steps)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()

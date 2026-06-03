# =============================================================================
# perturbflow/analyzer/cli.py
#
# Command-line interface for PerturbScope.
#
# Entry point:  python -m perturbflow.analyzer.cli  (or the `perturbscope` script after
#               `pip install -e .`)
#
# Subcommands:
#   run         -- Run the full pipeline or a subset of steps.
#   list-steps  -- Print the default step order and exit.
#   qc / score / effects / trajectory / interaction / report
#               -- Placeholder stubs; redirect users to `run --steps <name>`.
#
# Checkpoint / resume flags on `run`:
#   --resume              Skip steps already recorded in checkpoint.json.
#   --force-steps a,b     Re-run these steps even if checkpointed.
#   --clear-from step     Invalidate checkpoint from `step` onward, then run.
#   --list-steps          Print available pipeline steps and exit.
# =============================================================================

from __future__ import annotations

import argparse
from pathlib import Path

from .config import PipelineConfig
from .pipeline import run_pipeline
from .process_data import prepare_h5ad

_DEFAULT_STEPS = [
    "qc", "preprocess", "eda", "score", "effects",
    "trajectory", "programs", "interaction", "deg",
    "genenet", "tf_genenet", "report", "bundle",
]


def _parse_steps(text: str | None):
    if not text:
        return None
    return [s.strip() for s in text.split(",") if s.strip()]


def cmd_run(args):
    cfg = PipelineConfig.from_json(args.config) if args.config else PipelineConfig()
    steps = _parse_steps(args.steps)
    force_steps = _parse_steps(args.force_steps)
    run_pipeline(
        input_path=args.input,
        output_dir=args.output,
        config=cfg,
        steps=steps,
        perturbation_col=args.perturbation_col,
        resume=args.resume,
        force_steps=force_steps,
        clear_from=args.clear_from,
    )


def cmd_prepare(args):
    out = prepare_h5ad(
        args.input,
        args.output,
        perturbation_col=args.perturbation_col,
        control_labels=args.control_labels,
        cell_state_col=args.cell_state_col,
    )
    print(f"Wrote prepared AnnData: {out}")


def cmd_list_steps(_args):
    print("Available pipeline steps (default order):")
    for i, s in enumerate(_DEFAULT_STEPS, 1):
        print(f"  {i:2}. {s}")


def build_parser():
    parser = argparse.ArgumentParser(prog="perturbscope", description="PerturbScope CLI")
    sub = parser.add_subparsers(dest="command", required=True)

    # --- `run` subcommand ---
    run_p = sub.add_parser("run", help="Run end-to-end pipeline")
    run_p.add_argument("--input", required=True, help="Input .h5ad file")
    run_p.add_argument("--output", required=True, help="Output directory")
    run_p.add_argument("--config", default=None, help="Optional JSON config")
    run_p.add_argument(
        "--steps",
        default=None,
        help="Comma-separated steps to run, e.g. qc,preprocess,deg,report",
    )
    run_p.add_argument("--perturbation-col", default=None, help="Column in .obs with perturbation labels")
    run_p.add_argument(
        "--resume",
        action="store_true",
        default=False,
        help="Skip steps already in checkpoint.json (resume interrupted run)",
    )
    run_p.add_argument(
        "--force-steps",
        default=None,
        metavar="STEPS",
        help="Comma-separated steps to re-run even if checkpointed, e.g. deg,report",
    )
    run_p.add_argument(
        "--clear-from",
        default=None,
        metavar="STEP",
        help="Invalidate checkpoint from STEP onward, then continue running",
    )
    run_p.set_defaults(func=cmd_run)

    prep_p = sub.add_parser("prepare", help="Standardize a user .h5ad for PerturbScope")
    prep_p.add_argument("--input", required=True, help="Input .h5ad file")
    prep_p.add_argument("--output", required=True, help="Output standardized .h5ad file")
    prep_p.add_argument("--perturbation-col", required=True, help="obs column containing perturbation labels")
    prep_p.add_argument("--control-labels", default=None, help="Comma-separated extra labels to map to control")
    prep_p.add_argument("--cell-state-col", default=None, help="Optional obs column to copy to cell_state")
    prep_p.set_defaults(func=cmd_prepare)

    # --- `list-steps` subcommand ---
    ls_p = sub.add_parser("list-steps", help="Print available pipeline steps and exit")
    ls_p.set_defaults(func=cmd_list_steps)

    # --- Stub subcommands ---
    for name in ["qc", "score", "effects", "trajectory", "interaction", "report"]:
        p = sub.add_parser(name, help=f"Placeholder — use 'run --steps {name}'")
        p.add_argument("--input", required=False)
        p.add_argument("--output", required=False)
        p.set_defaults(func=lambda _args, step=name: print(f"Use 'run --steps {step}' for step-only execution."))

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    if hasattr(args, "output") and args.output:
        Path(args.output).mkdir(parents=True, exist_ok=True)
    args.func(args)


if __name__ == "__main__":
    main()

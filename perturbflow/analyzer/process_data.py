from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad


CONTROL_ALIASES = {
    "control",
    "ctrl",
    "non-targeting",
    "nontargeting",
    "nt",
    "scramble",
    "safe-targeting",
    "safe_targeting",
}


def prepare_h5ad(
    input_path: str | Path,
    output_path: str | Path,
    perturbation_col: str,
    control_labels: str | None = None,
    cell_state_col: str | None = None,
) -> Path:
    """Prepare a user AnnData file for PerturbFlow Analyzer.

    The pipeline expects perturbation labels in ``adata.obs["perturbation"]``.
    This helper keeps the original column and writes a standardized copy.
    """
    input_path = Path(input_path)
    output_path = Path(output_path)
    adata = ad.read_h5ad(input_path)

    if perturbation_col not in adata.obs.columns:
        raise KeyError(
            f"Column {perturbation_col!r} not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)[:25]}"
        )

    ctrl = set(CONTROL_ALIASES)
    if control_labels:
        ctrl.update(x.strip().lower() for x in control_labels.split(",") if x.strip())

    pert = adata.obs[perturbation_col].astype(str).str.strip()
    adata.obs["perturbation_original"] = pert.values
    adata.obs["perturbation"] = [
        "control" if str(x).lower() in ctrl else str(x) for x in pert.values
    ]

    if cell_state_col:
        if cell_state_col not in adata.obs.columns:
            raise KeyError(f"Column {cell_state_col!r} not found in adata.obs.")
        adata.obs["cell_state"] = adata.obs[cell_state_col].astype(str).values

    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)
    return output_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare a .h5ad file for PerturbFlow Analyzer.")
    parser.add_argument("--input", required=True, help="Input .h5ad file")
    parser.add_argument("--output", required=True, help="Output standardized .h5ad file")
    parser.add_argument("--perturbation-col", required=True, help="obs column containing guide/gene perturbation labels")
    parser.add_argument("--control-labels", default=None, help="Comma-separated additional labels to map to control")
    parser.add_argument("--cell-state-col", default=None, help="Optional obs column to copy to cell_state")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    out = prepare_h5ad(
        args.input,
        args.output,
        perturbation_col=args.perturbation_col,
        control_labels=args.control_labels,
        cell_state_col=args.cell_state_col,
    )
    print(f"Wrote prepared AnnData: {out}")


if __name__ == "__main__":
    main()

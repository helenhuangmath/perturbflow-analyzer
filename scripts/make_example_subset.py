#!/usr/bin/env python3
"""Build a small real-data AnnData fixture for PerturbFlow quick starts.

The source is the scPerturb-hosted Adamson/Weissman 2016 pilot Perturb-seq
dataset. This script downloads the public source .h5ad when needed and writes a
balanced subset that is small enough to run through PerturbFlow quickly.
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import tempfile
import urllib.request
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


SOURCE_URL = (
    "https://zenodo.org/records/13350497/files/"
    "AdamsonWeissman2016_GSM2406675_10X001.h5ad?download=1"
)
SOURCE_MD5 = "232f7e3756d41602bbe434b50662a76f"
SOURCE_NAME = "AdamsonWeissman2016_GSM2406675_10X001.h5ad"
DEFAULT_OUTPUT = "examples/data/adamson_2016_upr_360x1000.h5ad"
CONTROL_LABEL = "62(mod)_pBA581"
PERTURBATIONS = [
    CONTROL_LABEL,
    "SPI1_pDS255",
    "EP300_pDS268",
    "SNAI1_pDS266",
    "ZNF326_pDS262",
    "BHLHE40_pDS258",
]


def _md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp = destination.with_suffix(destination.suffix + ".part")
    with urllib.request.urlopen(url) as response, tmp.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    tmp.replace(destination)


def _gene_names(perturbation: str) -> str:
    if perturbation == CONTROL_LABEL:
        return "non-targeting"
    return perturbation.split("_", 1)[0]


def _select_balanced_cells(adata, cells_per_group: int, seed: int) -> list[str]:
    rng = np.random.default_rng(seed)
    selected: list[str] = []
    labels = adata.obs["perturbation"].astype(str)
    for label in PERTURBATIONS:
        names = adata.obs_names[labels == label].to_numpy()
        if len(names) < cells_per_group:
            raise ValueError(f"Only {len(names)} cells available for {label!r}")
        selected.extend(rng.choice(names, size=cells_per_group, replace=False).tolist())
    return selected


def _select_genes(adata, n_genes: int) -> list[str]:
    required = {
        _gene_names(label)
        for label in PERTURBATIONS
        if _gene_names(label) != "non-targeting"
    }
    var = adata.var.copy()
    if "ncells" in var:
        ranked = var["ncells"].astype(float).sort_values(ascending=False).index.tolist()
    elif "ncounts" in var:
        ranked = var["ncounts"].astype(float).sort_values(ascending=False).index.tolist()
    else:
        x = adata.X
        detected = np.asarray((x > 0).sum(axis=0)).ravel() if sparse.issparse(x) else (x > 0).sum(axis=0)
        ranked = pd.Series(detected, index=adata.var_names).sort_values(ascending=False).index.tolist()

    genes = [gene for gene in ranked if gene in required]
    genes.extend(gene for gene in ranked if gene not in genes)
    return genes[: min(n_genes, len(genes))]


def build_subset(source: Path, output: Path, cells_per_group: int, n_genes: int, seed: int) -> Path:
    adata = ad.read_h5ad(source)
    missing = sorted(set(PERTURBATIONS) - set(adata.obs["perturbation"].astype(str)))
    if missing:
        raise ValueError(f"Missing expected perturbations in source data: {missing}")

    cells = _select_balanced_cells(adata, cells_per_group=cells_per_group, seed=seed)
    genes = _select_genes(adata, n_genes=n_genes)
    subset = adata[cells, genes].copy()

    source_labels = subset.obs["perturbation"].astype(str)
    subset.obs["source_perturbation"] = source_labels.values
    subset.obs["guide_gene"] = source_labels.map(_gene_names).astype(str).values
    subset.obs["cell_state_hint"] = subset.obs.get("celltype", pd.Series("K562", index=subset.obs_names)).astype(str).values
    subset.obs["example_split"] = "quickstart"

    subset.uns["perturbflow_example"] = {
        "source": "scPerturb Zenodo record 13350497",
        "source_file": SOURCE_NAME,
        "source_url": SOURCE_URL,
        "source_md5": SOURCE_MD5,
        "publication": "Adamson et al., Cell 2016",
        "subset_strategy": (
            f"{cells_per_group} cells per perturbation/control group and top "
            f"{len(genes)} detected genes, with target genes retained when present"
        ),
        "random_seed": seed,
    }

    subset.obs_names_make_unique()
    subset.var_names_make_unique()
    output.parent.mkdir(parents=True, exist_ok=True)
    subset.write_h5ad(output, compression="gzip")
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create the PerturbFlow quick-start real-data subset.")
    parser.add_argument("--output", default=DEFAULT_OUTPUT, help="Output subset .h5ad path")
    parser.add_argument("--source", default=None, help="Existing source .h5ad path; skips download when present")
    parser.add_argument("--cache-dir", default=None, help="Directory for the downloaded source file")
    parser.add_argument("--cells-per-group", type=int, default=60, help="Cells per selected perturbation/control group")
    parser.add_argument("--n-genes", type=int, default=1000, help="Number of genes to keep")
    parser.add_argument("--seed", type=int, default=0, help="Reproducible cell sampling seed")
    parser.add_argument("--no-md5", action="store_true", help="Skip source MD5 validation")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = Path(args.output)
    if args.source:
        source = Path(args.source)
    else:
        cache_root = Path(args.cache_dir) if args.cache_dir else Path(tempfile.gettempdir()) / "perturbflow-example-data"
        source = cache_root / SOURCE_NAME
        if not source.exists():
            print(f"Downloading {SOURCE_NAME} from Zenodo...")
            _download(SOURCE_URL, source)

    if not source.exists():
        raise FileNotFoundError(source)
    if not args.no_md5:
        observed = _md5(source)
        if observed != SOURCE_MD5:
            raise ValueError(f"MD5 mismatch for {source}: expected {SOURCE_MD5}, observed {observed}")

    result = build_subset(
        source=source,
        output=output,
        cells_per_group=args.cells_per_group,
        n_genes=args.n_genes,
        seed=args.seed,
    )
    print(f"Wrote {result}")


if __name__ == "__main__":
    main()

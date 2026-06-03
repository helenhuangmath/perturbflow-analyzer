# Data Input

PerturbFlow expects an AnnData `.h5ad` file with cells as observations and genes as variables.

## Required Metadata

Provide one `.obs` column containing perturbation labels and pass it with `--perturbation-col`.

```bash
perturbflow prepare \
  --input raw.h5ad \
  --output prepared.h5ad \
  --perturbation-col target_gene
```

## Recommended Metadata

- `cell_state`, cluster, lineage, or annotation labels.
- Guide IDs if target genes and guides are separate.
- Replicate or batch labels.
- Precomputed QC metrics when available.

Control-like labels are mapped to `control` when passed through `--control-labels`.

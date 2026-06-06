# PerturbFlow

**Open infrastructure for perturbation biology.**

PerturbFlow is an open-source, AnnData-native platform that provides a unified
infrastructure layer for perturbation experiments — Perturb-seq, pooled CRISPR
screens, and single-cell multi-omics. Rather than being one more standalone
analysis method, it standardizes data representation, reproducible workflows,
mechanistic interpretation, and AI-ready outputs, and complements the scverse
ecosystem rather than replacing it.

A distinctive capability is **rewiring-aware interpretation**: PerturbFlow
distinguishes perturbations that *amplify* existing transcriptional programs
from those that *reorganize* regulatory relationships — interpretable signal
that prediction-accuracy metrics alone miss.

The project is organized around three aims (see [`ROADMAP.md`](ROADMAP.md)):

1. **Community-standard infrastructure** — stable AnnData schemas, reproducible
   workflows, interactive reports, reference datasets, and interoperable APIs.
2. **AI-ready perturbation biology** — model-ready exports, scalable processing,
   and interfaces for AI agents and workflow managers.
3. **Community benchmarks** — baseline-calibrated, distance-aware, and
   rewiring-aware evaluation with reproducible benchmarking workflows.

The current release ships the `perturbflow.analyzer` subpackage for QC,
perturbation scoring, differential expression, trajectory effects, gene-network
rewiring, regulatory analysis, interactive reports, and AI/agent-ready
interpretation handoff. The `perturbflow.predictor` and `perturbflow.benchmark`
namespaces are reserved for the prediction (Aim 2) and benchmarking (Aim 3)
work on the roadmap.

## What PerturbFlow Produces

- Standardized `.h5ad` input with consistent perturbation and optional cell-state annotations.
- Reproducible pipeline outputs: QC plots, DEG tables, trajectory summaries, program scores, gene networks, C-scores, regulatory results, and final AnnData.
- `report.html` and `interactive_report.html` for browser-based review.
- A viewer-ready `bundle/` directory for downstream web apps.
- `agent_handoff/` files that summarize the run for LLMs or analysis agents without including raw count matrices.

## Install

```bash
git clone https://github.com/helenhuangmath/PerturbFlow.git
cd PerturbFlow
python -m pip install -e ".[bundle]"
```

<details>
<summary>Cluster install (existing conda environment, no dependency reinstall)</summary>

On an HPC cluster with a pre-built environment, install in place without
re-resolving dependencies:

```bash
source /path/to/anaconda3/etc/profile.d/conda.sh
conda activate /path/to/your/perturbflow_env
python -m pip install -e /path/to/PerturbFlow --no-deps
```

</details>

## Quick Start

Prepare an AnnData file:

```bash
perturbflow prepare \
  --input my_raw_data.h5ad \
  --output prepared/my_data.perturbflow.h5ad \
  --perturbation-col guide_gene \
  --control-labels control,non-targeting,NT \
  --cell-state-col leiden
```

Run the analysis:

```bash
perturbflow analyzer \
  --input prepared/my_data.perturbflow.h5ad \
  --output results/my_run \
  --config configs/cluster_default.json \
  --resume
```

Open the main report:

```text
results/my_run/interactive_report.html
```

Create the AI/agent interpretation handoff:

```bash
perturbflow interpret \
  --results results/my_run \
  --project-name "K562 essential gene Perturb-seq"
```

This writes:

```text
results/my_run/agent_handoff/
├── agent_manifest.json
├── agent_prompt.md
├── interpretation_context.md
└── machine_context.json
```

Review these files before sharing them with an external LLM provider.

## Expected Input

Minimum input is an AnnData `.h5ad` file with cells in rows, genes in columns, and one `.obs` column containing perturbation labels.

Recommended optional columns:

- Cell state, cluster, or lineage label for state-aware interpretation.
- Guide ID when target gene and guide are separate.
- Replicate or batch labels for downstream review.
- Existing QC metrics if already computed.

PerturbFlow standardizes common control labels such as `control`, `ctrl`, `NT`, `non-targeting`, and `scramble`.

## Main Commands

```bash
perturbflow prepare      # Standardize input .h5ad metadata
perturbflow analyzer     # Run the analyzer workflow
perturbflow analyze      # Alias for analyzer
perturbflow run          # Legacy alias for analyzer
perturbflow predict      # Reserved for future predictor features (Aim 2)
perturbflow benchmark    # Reserved for community evaluation tooling (Aim 3)
perturbflow interpret    # Export LLM/agent-ready interpretation context
perturbflow list-steps   # Show available pipeline steps
```

The historical `perturbscope` command remains available for compatibility.

## Pipeline Steps

Default analysis steps include:

```text
qc -> preprocess -> eda -> score -> effects -> trajectory -> programs
-> interaction -> state_enrich -> deg -> genenet -> tf_genenet
-> cscore -> regulatory -> report -> bundle
```

Step-only reruns are useful while tuning reports:

```bash
perturbflow analyzer --input prepared/my_data.perturbflow.h5ad --output results/my_run --steps deg,report,bundle
perturbflow analyzer --input prepared/my_data.perturbflow.h5ad --output results/my_run --force-steps report --resume
```

## Repository Layout

```text
perturbflow/
├── perturbflow/          # Public package namespace and CLI
├── perturbflow/analyzer/ # Current analysis engine
├── perturbflow/data/     # Data preparation namespace
├── perturbflow/predictor/ # Reserved prediction namespace (Aim 2)
├── perturbflow/benchmark/ # Reserved benchmarking namespace (Aim 3)
├── perturbflow/workflows/ # End-to-end workflow namespace
├── perturbflow/viz/      # Visualization/reporting namespace
├── configs/              # Default and test pipeline configs
├── examples/             # Notebook templates for common workflows
├── scripts/              # Companion scripts, including Seurat/Mixscape
├── docs/                 # MkDocs documentation site
├── ROADMAP.md            # Aims, milestones, and status
├── CONTRIBUTING.md
├── README.md
├── QUICK_START.md
├── METHOD.md
├── RESULT.md
└── pyproject.toml
```

## Project Documents

- [`ROADMAP.md`](ROADMAP.md) — aims mapped to concrete work and status.
- [`CONTRIBUTING.md`](CONTRIBUTING.md) — how to get involved.
- [`DESIGN.md`](DESIGN.md) / [`METHOD.md`](METHOD.md) — architecture and methods.

## Web Documentation

The docs are built with MkDocs Material:

```bash
python -m pip install -e ".[docs]"
mkdocs serve
```

Then open the local URL printed by MkDocs. The docs structure is inspired by practical package documentation such as Seurat: installation, quick start, data preparation, analysis workflow, interpretation, and examples.

## Example Notebooks

Notebook templates are available in [`examples/`](examples/):

- `01_prepare_and_run.ipynb`: prepare data and run the full pipeline.
- `02_step_rerun_and_config.ipynb`: customize config and rerun selected steps.
- `03_interpret_with_agents.ipynb`: create AI/agent-ready interpretation files.
- `04_explore_outputs.ipynb`: inspect result tables, reports, and bundles.

## AI And Agent Design

PerturbFlow does not send data to any LLM service automatically. Instead, `perturbflow interpret` creates a compact handoff package with:

- A human-readable interpretation context.
- A reusable agent prompt.
- A machine-readable JSON summary.
- A manifest describing suggested agent roles.

This makes it possible to connect outputs to local LLMs, OpenAI-compatible APIs, custom agents, or collaborative report-writing workflows while preserving analyst control over privacy and provenance.

## Development

```bash
python -m pip install -e ".[dev,bundle]"
pytest
```

Generated result folders, large `.h5ad` files, caches, logs, and local notebooks are ignored by git by default.

# PerturbFlow

PerturbFlow is a modular Perturb-seq framework for moving from user-provided AnnData files to analysis, reporting, AI-assisted interpretation, and future perturbation response prediction.

The current release ships the `perturbflow.analyzer` subpackage for QC, perturbation scoring, differential expression, trajectory effects, gene-network rewiring, regulatory analysis, interactive reports, and AI/agent-ready interpretation handoff. The `perturbflow.predictor` namespace is reserved for future model-based prediction features.

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

On the Wherry cluster, install into the existing environment without reinstalling dependencies:

```bash
source /vast/parcc/spack/sw/apps/linux-sapphirerapids/anaconda3-2023.09-0-ieilyrkph5mewqcum3ajc4odlt2vakri/etc/profile.d/conda.sh
conda activate /vast/projects/wherry/foundation-models-immuno/hhua/tools/perturbscope_env
python -m pip install -e /vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbFlow --no-deps
```

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
perturbflow predict      # Reserved for future predictor features
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
├── perturbflow/predictor/ # Reserved future prediction namespace
├── perturbflow/workflows/ # End-to-end workflow namespace
├── perturbflow/viz/      # Visualization/reporting namespace
├── configs/              # Default and test pipeline configs
├── examples/             # Notebook templates for common workflows
├── scripts/              # Companion scripts, including Seurat/Mixscape
├── docs/                 # MkDocs documentation site
├── README.md
├── QUICK_START.md
├── METHOD.md
├── RESULT.md
└── pyproject.toml
```

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

# PerturbFlow Roadmap

This roadmap connects the project's strategic aims to concrete engineering work
in this repository. It is a living document; dates are
targets, not commitments, and the issue tracker is the source of truth for
status.

PerturbFlow is **open infrastructure for perturbation biology**, organized
around three aims.

## Aim 1 — Community-standard infrastructure for perturbation analysis

Stable AnnData schemas, reproducible workflows, interactive reports, reference
datasets, and interoperable APIs.

**Milestone — `v1.0`:** validated workflows across at least three public
perturbation datasets.

| Status | Item |
| --- | --- |
| ✅ Done | AnnData-native analyzer pipeline (`perturbflow analyzer`): QC, preprocess, EDA, scoring, effects, trajectory, programs, interaction, DEG, gene/TF networks, C-scores, regulatory |
| ✅ Done | `perturbflow prepare` schema standardization (perturbation + control-label normalization, optional cell-state) |
| ✅ Done | Static and interactive HTML reports |
| ✅ Done | Versioned results `bundle/` schema for downstream viewers |
| 🔜 Planned | Documented, versioned AnnData schema contract (`.obs` / `.uns` keys) |
| 🔜 Planned | Curated reference dataset loaders (≥ 3 public Perturb-seq datasets) |
| 🔜 Planned | Stable public Python API surface and semantic-versioning policy |

## Aim 2 — Model-ready perturbation biology

Standardized exports, scalable processing, and interfaces for workflow systems.

**Milestone (month 18):** public benchmark-ready dataset collection and an
interoperability toolkit.

| Status | Item |
| --- | --- |
| ✅ Done | `perturbflow interpret` structured handoff exports (no automatic data upload) |
| 🔜 Planned | Standardized model-ready export format for perturbation datasets |
| 🔜 Planned | Scalable / chunked processing for large screens |
| 🔜 Planned | Workflow-manager interfaces (e.g. Nextflow / Snakemake adapters) |

## Aim 3 — Community benchmarks for perturbation biology

Baseline-calibrated evaluation, biological-distance-aware diagnostics,
rewiring-aware metrics, and reproducible benchmarking workflows.

**Milestone (month 24):** public benchmark suite and evaluation portal.

| Status | Item |
| --- | --- |
| 🧪 Reserved | `perturbflow.benchmark` namespace reserved for evaluation tooling |
| 🔜 Planned | Baseline-calibrated metrics (compare methods against simple baselines) |
| 🔜 Planned | Biological-distance-aware and rewiring-aware evaluation metrics |
| 🔜 Planned | Reproducible benchmarking workflows and a public evaluation portal |

## How to get involved

See [`CONTRIBUTING.md`](CONTRIBUTING.md). Proposals that advance any of the
three aims are especially welcome — open an issue describing the use case before
sending a large change.

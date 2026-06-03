# PerturbFlow-Analyzer — Synthesized Design

Combines the scientific modules from `tmp_plan_1` and the phased build from `tmp_plan_2` (already partially implemented in [perturbflow/analyzer/](perturbflow/analyzer/)) with the bundle + interactive-viewer architecture from `tmp_plan_3`.

## Plan evaluation

- **plan_1** defines the science (effective perturbation scoring, dual-level effect decomposition, trajectory, programs, interactions). Adopted as the module catalog.
- **plan_2** turns those modules into a 9-phase build with concrete code, tests, and CI. Adopted as the build sequence — the current [perturbflow/analyzer/](perturbflow/analyzer/) tree already covers Phases 0–7 and a static-HTML Phase 8.
- **plan_3** addresses the part the other two miss: how the user actually *explores* results — a results-bundle contract plus a static viewer with perturbation, gene, and network pages. Adopted as the delivery architecture.

The user-facing requirements (interactive exploration, per-gene / per-perturbation lookup and comparison, TF and complex networks) are squarely the plan_3 layer. So the right move is **not** to pick one plan but to extend the existing plan_1/2 implementation with a plan_3 viewer.

## Target architecture

```
h5ad
  │
  ▼
PerturbFlow-Analyzer pipeline (existing perturbflow/analyzer/, Python)
  - ingest → qc → preprocess → score → effects → trajectory → programs → interaction
  │
  ▼
RESULTS BUNDLE (versioned directory, schema-stable)
  manifest.json, metadata.json, embeddings/*.parquet,
  perturbations/<pert>.parquet, de/<pert>.parquet,
  genes/<gene>.parquet, modules.json, tf_network.json,
  perturbation_similarity.parquet, search_index.json
  │
  ▼
STATIC VIEWER (HTML+JS, no backend)
  landing · browse perturbations · perturbation page · browse genes ·
  gene page · modules · TF/co-reg network · benchmarks · about
```

The pipeline and viewer are decoupled by the bundle schema. Same viewer works for any dataset. No live compute server required (one of plan_3's main advantages — easy to host, easy to share, archivable per-paper).

## What's already implemented

From [perturbflow/analyzer/](perturbflow/analyzer/) and [pipeline.py](perturbflow/analyzer/pipeline.py):

- ingest/QC/preprocess/score/effects/trajectory/programs/interaction modules
- CLI entry (`perturbflow run`)
- Static-HTML report ([report.py](perturbflow/analyzer/report.py))
- Bundle emitter ([bundle.py](perturbflow/analyzer/bundle.py)) — parquet + JSON results bundle
- **EDA module** ([eda.py](perturbflow/analyzer/eda.py)) — cells-per-group bar chart, cluster
  proportion stacked bar, gene×cell heatmap, gene×perturbation pseudobulk
  heatmap, gene–gene Pearson correlation heatmap, UMAP by cell-state cluster
- **DEG module** ([deg.py](perturbflow/analyzer/deg.py)) — identifies top perturbations, runs
  Welch t-test DEG analysis per perturbation vs control, volcano plots, DEG
  tables (CSV), and a top-perturbations log₂FC summary heatmap
- SLURM submission scripts

## What needs to be added

Two remaining items from the original plan:

### 1. TF / complex network module (extend `perturbflow/analyzer/programs.py`)

A new pipeline step `bundle` that runs after `interaction` and writes a typed, versioned results bundle. Inputs: the AnnData and DataFrames already in memory at end of [pipeline.py](perturbflow/analyzer/pipeline.py). Outputs:

| File | Source in current code | Notes |
|---|---|---|
| `manifest.json` | new | version, schema_version, list of available artifacts, dataset id |
| `metadata.json` | adata, summary.json | n_cells, n_genes, n_perturbations, qc summary |
| `embeddings/umap.parquet` | adata.obsm['X_umap'] + obs | cell × (UMAP1, UMAP2, perturbation, perturb_class, cell_state) |
| `perturbations/<pert>.parquet` | per-pert rows | one row per perturbed cell with effect scores, escape prob, signature score |
| `perturbations_summary.parquet` | adata.uns['perturbation_stats'] + effect_df + traj_df | one row per perturbation: n_cells, transcriptional_score, state_shift_score, fate_bias, top DE gene, etc. — drives the browse table |
| `de/<pert>.parquet` | within-state DEG output | gene × (log2fc, pval, padj, cell_state) — lazy-loaded per page view |
| `genes/<gene>.parquet` | inverted index over DE tables | which perturbations affect this gene + stats — drives the gene page |
| `modules.json` | program_df | module id → genes, scores, perturbation associations |
| `tf_network.json` | regulator activity from programs | TF → targets (with weight, direction); plus complex/module memberships |
| `coreg.parquet` | optional gene-gene similarity | sparse, only if computed |
| `perturbation_similarity.parquet` | new (cosine over effect profiles) | drives "similar perturbations" lookup |
| `search_index.json` | precomputed | FlexSearch/Lunr-shaped index of all genes + perturbations |

Bundle schema is versioned (`schema_version: "1.0"`); viewer reads `manifest.json` first and hides any panel whose artifact is missing — that keeps the same viewer reusable for partial runs.

### 2. TF / complex network module (extend `perturbflow/analyzer/programs.py`)

Current [programs.py](perturbflow/analyzer/programs.py) (3 KB) is light. To support the user's "find TF or complex networks" ask, extend it to emit:

- **TF activity per cell** — DoRothEA/CollecTRI prior × scaled expression (pure Python via [decoupler-py](https://github.com/saezlab/decoupler-py)), stored in `adata.obsm['tf_activity']`.
- **TF–target edges** for the network panel — TF, target, weight (prior), perturbation_effect (delta TF activity perturbed vs control), serialized as `tf_network.json`.
- **Complexes / modules** — NMF on the residualized expression (perturbed minus control mean per state) → modules with member genes and per-perturbation activation scores → `modules.json`.

This is the only real addition to the science layer; everything else for the user's ask is delivery.

### 3. Static viewer (`viewer/`, separate from `perturbflow/analyzer/`)

**Recommendation: Quarto for v1.** Reasons matching plan_3: a single command builds the site, Plotly/Observable plots are first-class, GitHub Pages deploy is trivial, days-not-weeks to ship. Move to a Svelte/React frontend with Arrow.js + DuckDB-WASM only if/when interactivity outgrows it (large UMAPs, cross-page state, scale to 10k+ perturbations).

Pages, all reading the bundle, none requiring a server:

- **Landing** — counts, dataset summary, UMAP overview, citation, bundle download.
- **Browse perturbations** — searchable/sortable table from `perturbations_summary.parquet`. Columns: perturbation, n_cells, n_perturbed, transcriptional_score, state_shift_score, dominant_effect, top DE gene, fate_bias. Click row → perturbation page.
- **Perturbation page** (one per pert; lazy-loads `perturbations/<pert>.parquet` and `de/<pert>.parquet`) — volcano, DE table, UMAP highlighting that perturbation vs control, effect-decomposition bar (transcriptional vs state-shift), trajectory shift, top similar perturbations (click-through), pathway enrichment.
- **Browse genes** — search box; jumps to gene page.
- **Gene page** — which perturbations affect this gene (table with log2fc, padj, dominant cell state), expression distribution, network neighborhood subgraph from `tf_network.json` / `coreg.parquet`.
- **Modules** — module list from `modules.json`; each module shows member genes, top perturbations, GO/Reactome enrichment.
- **Networks** — TF-target view + co-regulation view, rendered with Cytoscape.js. Critical: never render the full graph — always subgraph around a search query (gene or TF) to a depth-2 neighborhood. Full networks crash browsers.
- **Compare** — pick 2–N perturbations from the browse table, get a side-by-side: shared/unique DEGs, cosine similarity, effect-score scatter, shared module activations. This is the "extract data and compare" ask.
- **Benchmarks** (optional) — model × perturbation × metric table if/when prediction models are added.
- **About** — methods, version, citation, bundle schema link.

What "interactive" gets you in this static setup (per plan_3): hover/zoom/click on every plot via Plotly, cross-filtering when you use linked Vega-Lite/Bokeh selections, lazy parquet loading per page, URL-state for bookmarkable perturbation/gene pages, in-browser search via the precomputed index. What you don't get: re-running DE with new thresholds, user accounts. If those become requirements later, that's the moment to add a Streamlit-on-HF-Spaces backend; not before.

## Phased build (delta on top of existing code)

Each phase is one PR-sized chunk. Phases 0–7 of `tmp_plan_2` are already done; EDA and DEG are now implemented; this tracks remaining work.

| Phase | Deliverable | New code | Status |
|---|---|---|---|
| **A** | Bundle schema + emitter | [bundle.py](perturbflow/analyzer/bundle.py) | ✅ done |
| **EDA** | Cells-per-group, 3 heatmaps, cluster proportions, UMAP by cluster | [eda.py](perturbflow/analyzer/eda.py) | ✅ done |
| **DEG** | Top-perturbation DEG, volcano plots, summary heatmap | [deg.py](perturbflow/analyzer/deg.py) | ✅ done |
| **B** | TF / module / network artifacts | extend [programs.py](perturbflow/analyzer/programs.py): decoupler TF activity, NMF modules, TF-network JSON serializer | pending |
| **C** | Quarto viewer skeleton | `viewer/` with landing, browse-perturbations, perturbation-page (Plotly), reading bundle | pending |
| **D** | Gene page + search | inverted index in bundle, gene page, FlexSearch wiring | ~2 days |
| **E** | Network page (Cytoscape.js, subgraph-on-query) | `viewer/network.qmd` + JS | ~2–3 days |
| **F** | Compare page | multi-select from browse table, side-by-side DE / similarity / module overlap | ~2 days |
| **G** | CI: pipeline runs on tagged release → emits bundle → builds + deploys viewer to GitHub Pages | one workflow file | ~1 day |

A through D get you the interactive perturbation/gene exploration the user asked for; E adds TF/network; F adds comparison; G makes it shareable.

## Bundle schema — concrete starting point

```json
// manifest.json
{
  "schema_version": "1.0",
  "perturbscope_version": "0.1.0",
  "dataset_id": "ReplogleWeissman2022_K562_essential",
  "created_utc": "2026-04-26T00:00:00Z",
  "artifacts": {
    "metadata": "metadata.json",
    "embeddings": ["embeddings/umap.parquet"],
    "perturbations_summary": "perturbations_summary.parquet",
    "perturbations": "perturbations/",
    "de": "de/",
    "genes": "genes/",
    "modules": "modules.json",
    "tf_network": "tf_network.json",
    "perturbation_similarity": "perturbation_similarity.parquet",
    "search_index": "search_index.json"
  }
}
```

```text
# perturbations_summary.parquet — one row per perturbation
perturbation: str
n_cells: int32
n_perturbed: int32
n_escaped: int32
transcriptional_score: float32
state_shift_score: float32
dominant_effect: enum{transcriptional, state_shift, both, neither}
fate_bias_score: float32
top_de_gene: str
top_de_log2fc: float32
n_significant_de: int32
```

```text
# de/<pert>.parquet — gene × stats, one file per pert (lazy-loaded)
gene: str
cell_state: str        # within-state DEG result
log2fc: float32
pval: float32
padj: float32
mean_expr_perturbed: float32
mean_expr_control: float32
```

```text
# genes/<gene>.parquet — inverted index for the gene page
perturbation: str
cell_state: str
log2fc: float32
padj: float32
direction: enum{up, down}
```

```text
# tf_network.json (subset)
{
  "tfs": [{"name": "STAT1", "n_targets": 124, "perturbation_delta": {"IFNG": 0.42, ...}}, ...],
  "edges": [{"tf": "STAT1", "target": "IRF1", "weight": 0.81, "source": "CollecTRI"}, ...],
  "complexes": [{"name": "NF-kB", "members": ["RELA","RELB","NFKB1"], ...}]
}
```

## Decisions to lock before Phase A

These are the small forks that change the design downstream — worth picking now rather than refactoring later:

1. **One bundle per dataset, or one bundle per (dataset, model_version)?** I'd recommend the latter — keeps reproducibility tight and lets the viewer show two runs side-by-side later.
2. **Quarto vs Streamlit-on-Stlite vs custom Svelte.** Default Quarto unless you already know you need cross-page reactive state.
3. **Bundle hosting.** Same repo as code, separate `bundles/` repo, or release artifact attached to a GitHub release. Release artifact is the cleanest if bundles get large (>100 MB).
4. **Network prior source.** CollecTRI for human (default), DoRothEA for legacy mouse work, or both with a switch. Affects [programs.py](perturbflow/analyzer/programs.py).
5. **TF activity scoring method.** decoupler `ulm`/`mlm`/`viper`. `ulm` is fast and the de-facto default; flag if you want viper.

## What I would not do

- Don't build a real database (Postgres/SQLite). Parquet + a static viewer covers the "extract per-gene/per-perturbation and compare" use case at zero ops cost. Add a DB only if you outgrow this — concretely, when you want server-side joins across hundreds of bundles, which you don't yet.
- Don't pursue Stlite/Pyodide unless you specifically need Python in the browser. The 5–15 s cold-start is felt by every visitor.
- Don't render full TF networks; subgraph-on-query is the only thing that scales.
- Don't change the existing [pipeline.py](perturbflow/analyzer/pipeline.py) module contracts in this work; just add `bundle` as an additional terminal step.

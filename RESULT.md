# PerturbFlow Analyzer — Results & Interpretation Guide

This document is a reader's manual for the analyses produced by
**PerturbFlow Analyzer** (the analysis pipeline in this repository, also
referred to internally as `PerturbFlow Analyzer`). For every analysis step it
describes:

* **What it shows** — what the underlying computation is and what the
  output represents.
* **How to read it** — concrete things to look for in the figure / table.
* **Caveats** — failure modes and common ways to be misled.

The companion document [`METHOD.md`](METHOD.md) describes the algorithms
and parameters; this document is purely about **interpretation**. All
file paths below are relative to a run's output directory (e.g.
`results_test/`).

The pipeline writes three top-level subfolders:

| folder        | contents                                                                                    |
|---------------|---------------------------------------------------------------------------------------------|
| `plots/`      | publication-ready PNG/PDF figures, one per analysis                                         |
| `csv/`        | underlying tables (every figure has a CSV behind it)                                        |
| `bundle/`     | versioned, viewer-ready Parquet/JSON dump for downstream apps                               |
| `report.html` + `interactive_report.html` | static and interactive summary reports that embed the figures and tables |

The `summary.json` at the top of the output directory lists every CSV
and PNG that was actually produced — handy for scripted comparisons.

---

## 1. Quality control — `qc_summary.png`

**What it shows.** Per-cell distributions of the standard scRNA-seq QC
metrics (`n_genes`, `n_counts`, `pct_mt`) before and after the QC filter
defined in `configs/*.json` (`min_genes`, `max_pct_mt`,
`min_cells_per_perturbation`).

**How to read it.**

* The "after" violins should be visibly trimmed compared to "before".
* `pct_mt` should be ≤ `max_pct_mt` for every retained cell — if not, the
  filter did not run or the `mt_pattern` regex did not match the gene
  symbols.
* Per-perturbation cell counts (`csv/eda_cells_per_perturbation.csv` and
  `eda_cells_per_perturbation.png`) should all be ≥
  `min_cells_per_perturbation`.

**Caveats.** Very low `n_counts` cells often survive the filter when the
threshold is set to the default 200; this can dominate downstream
correlations if the perturbation library is shallow.

## 2. Embedding & UMAPs — `umap_cell_state.png`, `umap_perturbation.png`

**What it shows.** Two views of the same UMAP embedding (computed in the
`preprocess` step from a 50-PC PCA on log-normalised counts):

* `umap_cell_state.png` — coloured by Leiden cluster (the unsupervised
  cell-state partition, resolution `leiden_resolution`).
* `umap_perturbation.png` — coloured by perturbation label.

**How to read it.**

* If perturbations co-locate in the embedding (e.g. all RPL/RPS knock-downs
  pile together), the perturbations have a coherent transcriptomic
  signature — a precondition for the rest of the pipeline.
* If perturbations spread *across* the cell-state UMAP, expect strong
  **state-conditional** effects (the perturbation's effect depends on
  which Leiden cluster the cell is in). The `state_enrich` step
  quantifies that.

**Caveats.** UMAP geometry is not a metric; do not over-interpret cluster
shapes or distances. Use the `eda_perturbation_similarity.png` heatmap
(cosine similarity in PCA space) for quantitative comparisons.

## 3. Cells per perturbation / cluster proportions

`eda_cells_per_perturbation.png` + `csv/eda_cells_per_perturbation.csv`
— bar chart of cell counts per perturbation, control highlighted in
neutral grey (`#9aa0a6`). The first sanity check: are any
perturbations under-sampled?

`eda_cluster_proportions.png` + `csv/eda_cluster_proportions.csv` —
stacked bar of cell-state cluster fractions for the top 30 perturbations.

* Bars that look very different from the control bar's stack indicate
  cell-state shifts — perturbation pushes cells into / out of a
  particular Leiden state.
* Use this together with §6 (`state_enrich`) to test whether the shift
  is statistically significant.

## 4. EDA heatmaps

The EDA step writes five heatmaps. All of them use **z-scored** values
(per gene across the visible columns), `vmin=−3`, `vmax=3`.

### 4a. `eda_gene_by_cell_heatmap.png`
Top `eda_n_top_genes` HVGs (rows) × cells (columns, sub-sampled to
`eda_max_cells_heatmap`). Column-side sidebar colours mark perturbation
membership; control cells are in neutral grey.

**Read it as:** "Do cells of the same perturbation produce visibly
similar expression columns?" — vertical bands of consistent colour =
strong, coherent perturbation signature.

### 4b. `eda_gene_by_perturbation_heatmap.png`
Genes × perturbation pseudobulk means. No row/column clustering.

### 4c. `eda_clustered_gene_by_pert_heatmap.png` ("All Genes × Perturbations")
Same as 4b but rows and columns are hierarchically clustered (Ward),
then flat-cut into `eda_n_gene_clusters` and `eda_n_pert_clusters`
groups; cluster IDs appear as coloured sidebars and are written to
`csv/eda_gene_clusters.csv` / `csv/eda_pert_clusters.csv`.

**Read it as:** the row sidebar groups co-regulated genes (the
"programs" of the system), the column sidebar groups perturbations with
similar phenotypic effects. The two cluster-CSVs are the formal output
that downstream analyses (`programs`, `interaction`) consume.

### 4d. `eda_gene_pert_cluster_summary_heatmap.png`
Reduced view of 4c: each cell of the heatmap is the **mean z-score of
all genes in gene-cluster GCi for all perturbations in pert-cluster
PCj**.

**Read it as:** a "what programme moves under what kind of
perturbation" matrix. Strong red cells (large positive mean z) mean
the gene module is up-regulated by that perturbation cluster as a
group; blue means down-regulated. Annotated values inside cells are
the actual mean z-scores.

### 4e. `eda_gene_correlation_heatmap.png` ("Gene–Gene Correlation (all cells)")
Square heatmap of the top-`eda_n_top_genes` HVG × HVG Pearson
correlation across all cells, sorted by hierarchical-clustering modules
(distance = `1 − |r|`, average linkage). The top + left sidebars show
module membership; a "Gene module" legend lives below the colour bar.

**Read it as:** the diagonal blocks ARE the gene modules. A clean
block-diagonal pattern with sharp boundaries means modules are
well-defined; a blurry plot means the HVG panel is dominated by one or
two large correlated programmes (often cell cycle).

`eda_gene_corr_vs_<pert>.png` — same plot but on cells of `<pert>`
only, for the top perturbations. Compare side-by-side with the all-cell
version: blocks that *appear* under perturbation are induced
co-regulation; blocks that *vanish* are loss of co-regulation.

### Common caveats
* All EDA heatmaps depend on the **HVG selection** — if cell cycle
  dominates the variance ranking (as it usually does in K562/RPE1
  essential screens), the plots will be cell-cycle heavy.
* Z-scoring is per-column (per gene) so absolute magnitudes are not
  comparable across the row of a heatmap.

## 5. Perturbation similarity — `eda_perturbation_similarity.png`

**What it shows.** Pairwise cosine similarity of perturbation pseudobulk
profiles in PCA space; perturbations on rows and columns; symmetric.

**How to read it.**

* High similarity blocks (red squares) = perturbations with the same
  transcriptomic phenotype. Confirms biology (e.g. all spliceosome
  knock-downs grouped together) or flags duplicates.
* The **control row/column** should be similar to itself only —
  perturbations strongly correlated with control are *non-effective*.
* Use it together with `cscore_vs_deg.png` and the C-score ranking
  (§11) to find the strongest perturbations to focus on.

## 6. Cell-state enrichment / depletion (`state_enrich`)

When the `state_enrich` step ran, the interactive report's "Cell State
Enrichment / Depletion" tab shows a heatmap of (perturbation × Leiden
state) odds-ratios with a chi-square FDR.

**How to read it.**

* Red cells = perturbation pushes cells *into* that state (enriched).
* Blue cells = perturbation pushes cells *out of* that state (depleted).
* Stars / outlines mark cells passing
  `state_enrich_fdr_threshold`; un-marked colour is descriptive only.

**Caveats.** Enrichment is computed against the control composition;
small per-perturbation cell counts give very wide odds-ratio confidence
intervals — focus on FDR-passing cells.

## 7. Effects: `perturbation_score_bar.png` + `csv/effect_decomposition.csv`

**What it shows.** The `effects` step decomposes the per-perturbation
effect into

```
total = bulk + state-conditional + interaction + residual
```

(see [`perturbflow/analyzer/state.py`](perturbflow/analyzer/state.py)). The bar plot
ranks perturbations by total effect size; the CSV is the data behind
it.

**How to read it.**

* High **bulk** but low **state-conditional** = the perturbation has the
  same effect everywhere → genuine transcriptional driver.
* High **state-conditional** = the effect of the perturbation is
  state-dependent → look at the corresponding row of the
  `state_enrich` heatmap to see *which* state.
* High **residual** = the perturbation produces heterogeneous effects
  not captured by the cell-state partition; consider re-running with a
  finer Leiden resolution or a different embedding.

## 8. Programs (`csv/program_scores.csv`)

The `programs` step fits gene-set programmes (PCA of the cluster-mean
matrix from §4c) and scores every cell. The CSV is `(cell × program)`;
the interactive report's "Programs" tab plots the per-perturbation
distribution of each program score.

**How to read it.**

* A perturbation with a sharply non-zero distribution for program N (vs
  the control's distribution) selectively activates / deactivates that
  programme.
* Cross-reference the loadings of program N with the gene-cluster
  membership in `csv/eda_gene_clusters.csv` to give the program a
  biological label.

## 9. Trajectory effects — `trajectory_shift.png`

**What it shows.** The `trajectory` step fits a 1-D principal-curve
pseudotime through the embedding and reports, for each perturbation,
the median pseudotime shift relative to control.

**How to read it.**

* Perturbations with positive shifts move cells "forward" along the
  trajectory, negative shifts move them "backward". In a
  proliferation-vs-differentiation axis (the typical K562 trajectory),
  positive shifts ≈ pushing cells toward the proliferation tail.
* Use this view sparingly — it depends on a unimodal pseudotime, which
  is rarely true beyond simple lineages.

## 10. DEG (`deg`) — `deg_volcano_<pert>.png`, `deg_top_perturbations_heatmap.png`

**Volcano plots.** x-axis = `log2FC` of perturbation vs control,
y-axis = `−log10(adj. p-value)`; significance threshold drawn at
`deg_logfc_threshold` and `deg_pval_threshold`. Significant up-regulated
genes are red, down-regulated are blue, non-significant are grey. Two
vertical lines mark the **mean** and **median** of `log2FC` so you can
see whether the perturbation is dominantly inducing (right of zero) or
dominantly repressing (left of zero).

**Read it as:**

* Many points in the upper-right and upper-left = strong, balanced
  bidirectional response.
* Asymmetric volcano (most points on one side) = transcriptional
  programme is monotonic — typical for ribosome / proteasome knock-downs
  which broadly down-regulate translation.
* The CSV [`csv/deg_<pert>.csv`](csv/) has the underlying gene list.

**Top-perturbations heatmap** (`deg_top_perturbations_heatmap.png`).
Top-`deg_n_top_heatmap` DEGs per perturbation, z-scored across
perturbations. Lets you see at a glance which DEGs are perturbation-
specific vs shared.

**Pathway enrichment** (`csv/deg_enrichment_<pert>.csv` and bar plot in
the report's "Perturbation" tab). MSigDB / KEGG enrichment per pert,
−log10(FDR) on the x-axis with a dashed line at FDR = 0.05.

**Caveats.** Volcano `padj` comes from a per-pert Wilcoxon test against
control with BH correction; with very few cells per perturbation, the
test is conservative and many true effects fail to reach significance.
Always inspect `log2FC` magnitudes alongside the p-value column.

## 11. Connectivity score (C-score) — `cscore_*`

The `cscore` step measures **how much the gene–gene co-expression
network rewires** under each perturbation. For each perturbation it
computes a scalar score and decomposes it into module-level
contributions.

* `cscore_ranked_bar.png` — perturbations ranked by total C-score.
  Top of the bar = most-rewiring perturbations.
* `cscore_decomposition.png` — stacked bars showing how much of the
  total C-score comes from gain-of-edge vs loss-of-edge events.
* `cscore_module_rewiring.png` — heatmap (modules × perturbations) of
  per-module rewiring; pinpoints *which* module is being remodelled.
* `cscore_gene_hub_heatmap.png` — top hub genes per perturbation; rows
  are genes whose connectivity changes the most under the perturbation.
* `cscore_vs_deg.png` — scatter plot of C-score vs total DEG count per
  perturbation. Useful to find perturbations that **rewire the network
  without changing many genes' bulk levels** (high C-score, low DEG) —
  these are typically regulatory perturbations rather than direct
  transcriptional drivers.

**Read it as:** C-score answers "*is the wiring different?*", the DEG
table answers "*are the levels different?*" — the two are complementary
and the scatter is the most useful single view.

**Caveats.** The C-score is sensitive to per-perturbation cell counts
(small samples produce noisy correlations and inflate apparent
rewiring). Always read `cscore_vs_deg.png` alongside the per-pert cell
counts plot.

## 12. Gene co-expression networks (`genenet`) — see §1 of `METHOD.md`

For each top perturbation the `genenet` step produces:

* `genenet_control_heatmap.png` — control gene-gene |r| matrix sorted
  by the 5 modules; module sidebar; "Module N" legend.
* `genenet_<pert>_ctrl_network.png`, `genenet_<pert>_pert_network.png`,
  `genenet_<pert>_diff_network.png` — three square spring-layout
  panels, sharing node positions across the row.
* `genenet_<pert>_heatmap_comparison.png` — side-by-side row-z-scored
  control / perturbation / log2FC heatmap, gene rows sorted by module.
* `csv/genenet_<pert>_network.json` — node positions + edges for the
  interactive viewer.

**How to read each panel.**

* **Control panel** — the baseline wiring. Modules should appear as
  visually-clumped groups of same-coloured nodes; densely-connected
  hubs are likely co-regulated by a shared upstream factor.
* **Perturbation panel** — same gene set, same layout, but recomputed
  edges. Read it side-by-side with the control panel and look for
  edges that **change**; the topology summary is in the differential
  panel.
* **Differential panel** — Okabe-Ito blue (`#0072B2`) edges = **gained**
  in perturbation, vermilion (`#D55E00`) = **lost**, neutral grey =
  **unchanged**. Node size is `|Δ mean expression|`. A perturbation
  that "tears apart" a module shows a cluster of orange edges around
  that module's nodes.

**Caveats** (also in `METHOD.md` §10).

* Edges are thresholded at `|r| ≥ corr_threshold`; an edge crossing the
  threshold under perturbation is *not* a significance test. Treat
  individual gained/lost edges as hypotheses.
* Small per-pert cell counts increase the noise floor of the
  correlation; a `cscore_vs_deg.png` outlier with a *low* DEG count and
  *high* C-score may simply be sampling noise.

## 13. TF-anchored gene network (`tf_genenet`) — new in this version

A sibling of `genenet` whose node panel is **the top 20 highly variable
TFs (HumanTFs v1.01) plus their top correlated partner genes in
control**. Outputs:

* `plots/tfnet_control_network.png` — shared control panel (all
  perturbations have the same control here, since the seed is built
  from control cells).
* `plots/tfnet_<pert>_pert_network.png` — same panel under
  perturbation.
* `plots/tfnet_<pert>_diff_network.png` — set-difference of the two
  graphs (gained / lost / shared, same colour code as `genenet` diff).
* `csv/tfnet_seed_tfs.csv` — TF, control variance, control mean.
* `csv/tfnet_partners.csv` — partner gene → top-TF + signed Pearson r.
* `csv/tfnet_<pert>_network.json` — same schema as `genenet_*_network.json`
  plus a per-node `is_tf` flag.

**How to read it.**

* TFs are drawn as **squares**, partner genes as **circles**; node
  colour reflects gene-module membership computed from the panel's
  control |r| matrix; node size = mean expression in the relevant
  condition (or |Δ expression| in the differential panel).
* Edges follow the same rules as `genenet` — sign of `r` (red /
  blue) in the static panels, status (blue / orange / grey) in the
  differential panel.
* The **TF-vs-partner shape distinction** is the value-add over
  `genenet`: a "lost" edge between a TF (square) and a partner
  (circle) under perturbation is a candidate **direct regulatory
  link** that has been disrupted; a "gained" edge between two TFs
  (two squares) is a candidate **rewired regulatory hub**.

**Caveats.**

* The seed list is purely "highly variable in control" — variance ≠
  regulatory importance. A constitutively-high TF will be missed.
* Partners are picked by simple `|r|` ranking; the same gene can be a
  partner of multiple TFs.
* This is still **co-expression**, not causality. To turn the
  TF-anchored panel into a regulatory call you still need an external
  TF→target prior (CollecTRI / DoRothEA). See §11 of `METHOD.md` for
  the proposed overlay.

## 14. Regulatory analysis (`regulatory`) — `tf_regulatory_heatmap.png`

The `regulatory` step builds a **TF–TF regulatory map from
perturbation effects** rather than from cell-level correlations. For
every perturbed TF, it reads the DEG table and assigns a signed
−log10(q) to every other gene; the heatmap is the resulting (perturbed
TF × gene) matrix.

**How to read it.**

* Red cell = the perturbed TF's KO **increased** the target's
  expression → the TF normally **inhibits** that target.
* Blue cell = the KO **decreased** the target → the TF normally
  **activates** the target.
* Rows that look identical = TFs whose KO has the same downstream
  effect = candidate redundant regulators.

`csv/tf_regulatory_matrix.csv` is the underlying signed −log10(q)
matrix; `csv/tf_regulatory_network.json` is the directed-edge list
consumed by the interactive report's "Regulatory" tab (which renders a
network where edges go from a perturbed TF to its strongest targets).

**Caveats.**

* This is a **functional** regulatory map — it tells you what changes
  when you knock the TF out, not what the TF binds. Unlike `genenet`
  (co-expression) or `tf_genenet` (TF-seeded co-expression), each cell
  here is testable causal evidence at the level of "perturbing X
  changes Y".
* Non-TFs that happened to be perturbed will not have a regulatory row
  — the matrix is restricted to perturbations whose target gene is in
  HumanTFs.

## 15. Interactive report (`interactive_report.html`)

A self-contained HTML viewer that bundles every figure and CSV from the
above steps into nine tabs:

| Tab | Sources |
|---|---|
| Home | summary metrics, top-perturbations chart, DEG-per-pert bar, effect scatter |
| QC | `qc_summary` + per-pert violin grid |
| UMAP | `umap_perturbation` / `umap_cell_state` colour toggles + free-text gene colouring |
| Heatmaps | the five EDA heatmaps + per-pert gene-correlation switcher |
| Perturbation | volcano, UMAP highlight, per-pert DEG bar, pathway bar, paginated DEG table |
| Gene Network | three square network panels (control / pert / diff) + control co-expression + per-pert expression heatmap |
| Gene | per-gene boxplot + per-pert log2FC bar |
| States | `state_enrich` heatmap |
| Regulatory | `tf_regulatory_heatmap` + clickable TF-TF network |

The HTML is static (no server needed) and ships with every PNG embedded
as base64. Open it directly in a browser; everything is keyboard-
navigable.

## 16. Bundle (`bundle/`)

Versioned, viewer-ready dump suitable for downstream apps. Files:

* `manifest.json` — schema version + list of every artefact + cell /
  gene / perturbation counts.
* `cell_level.parquet` — `(cell × {perturbation, cell_state, embedding,
  total_counts, ...})` long-form table.
* `pert_level.parquet` — `(perturbation × derived metrics)`.
* `deg.parquet` — top-`bundle_top_de_per_pert` DEGs per perturbation.
* `gene_index.parquet` — `(gene × n_perturbations_in_DEG)` to power
  the gene browse list.

The interactive report consumes the bundle for fast paginated table
rendering; you can also load it directly into any DataFrame-aware
viewer.

---

## Reading order for a new dataset

A practical sequence to interpret a fresh PerturbFlow Analyzer run:

1. **Sanity checks.** Open `qc_summary.png`,
   `eda_cells_per_perturbation.png`, and `eda_perturbation_similarity.png`
   — confirm cells passed QC, every perturbation has enough cells, and
   the similarity matrix is not dominated by a single block (which
   would indicate a confounder).
2. **Where do perturbations live?** Open `umap_perturbation.png` and
   `umap_cell_state.png` side-by-side. Are perturbations co-located in
   embedding space?
3. **What does each perturbation do?** Open the "Perturbation" tab of
   the interactive report (or `deg_volcano_<pert>.png` per pert) for the
   strongest perturbations from step 1; read off direction +
   magnitude + pathway enrichment.
4. **Why does it do that?** Open the "Gene Network" tab → look at the
   differential network for the same perturbation. Edges in vermilion
   (lost) cluster around the genes the volcano flagged.
5. **Is it a regulator?** Open `tf_regulatory_heatmap.png` and the
   TF-anchored `tfnet_<pert>_diff_network.png` — if the perturbation
   moves many targets and the TF-anchored network shows TF–partner
   edges flipping, this is a regulatory perturbation; if neither, it is
   a structural / metabolic perturbation.
6. **Is the wiring different?** Read `cscore_vs_deg.png` to see whether
   the perturbation's effect is "level changes" (high DEG / low
   C-score) or "wiring changes" (low DEG / high C-score).

Steps 1–6 cover most workflow questions; the deeper analyses
(`programs`, `trajectory`, `state_enrich`, `interaction`) are useful for
follow-up questions about specific perturbation classes.

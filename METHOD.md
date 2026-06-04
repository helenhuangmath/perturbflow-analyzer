# PerturbFlow Analyzer — Methods

This file documents the analysis methods that PerturbFlow Analyzer runs on a
Perturb-seq / CRISPR-screen `AnnData`. It is meant to be the
publication-style "Methods" companion of `DESIGN.md` (which describes the
software architecture). For each step we state the inputs, the algorithm,
the parameters that are exposed in the JSON config, and where the
artefacts land on disk. The pipeline is wired in
[`perturbflow/analyzer/pipeline.py`](perturbflow/analyzer/pipeline.py); each step lives in
its own module under `perturbflow/analyzer/`.

> **Scope of this document.** Below we cover the gene co-expression /
> gene regulatory network (GRN) step in detail, because it is the most
> involved single-step analysis in the package. The other steps (QC,
> preprocessing, EDA, DEG, etc.) are described at a high level in
> `DESIGN.md` and `README.md`.

---

## Gene co-expression / regulatory network (GRN) step

**Module:** [`perturbflow/analyzer/genenet.py`](perturbflow/analyzer/genenet.py)
**Public entry point:** `run_gene_network(adata, output_dir, perturbations, n_top_genes, n_gene_clusters, corr_threshold)`
**Pipeline call site:** [`perturbflow/analyzer/pipeline.py:204-217`](perturbflow/analyzer/pipeline.py#L204-L217)
**Step name (for `--steps` / `default_steps`):** `genenet`
**Test config:** `genenet_n_top_genes=30`, `genenet_n_gene_clusters=4`, `genenet_corr_threshold=0.4` (see [`configs/test.json`](configs/test.json)).

### 1. Inputs and assumptions

| Input | Source | Notes |
|---|---|---|
| `adata.X` | log1p-normalised counts after the `preprocess` step | Cell-level expression; the GRN step does *not* re-normalise. |
| `adata.obs["perturbation"]` | string column | One row per cell; `control` is reserved for the unperturbed pool. |
| `adata.obs["is_control"]` | boolean | Set by `qc` / `preprocess`. Identifies the cells used as the control pool. |
| `adata.var["highly_variable"]` | boolean (optional) | If present, the gene panel is taken from this flag; otherwise we fall back to top variance over control cells (see §2). |
| `top_perts` (parameter) | from [`identify_top_perturbations`](perturbflow/analyzer/perturbation.py) | The `n_top` perturbations ranked by the `effects` step's effect-size table. The test run uses `deg_n_top_perturbations = 3`. |

The GRN step requires **≥ 5 control cells and ≥ 5 cells per perturbation**; smaller groups are skipped with no error.

### 2. Node set — gene panel selection

The same node set is used in the control graph and in every per-perturbation
graph, so node positions are comparable across panels.

```python
if "highly_variable" in adata.var.columns:
    hv_idx = np.where(adata.var["highly_variable"].values)[0][:n_top_genes]
else:
    xc_full = _to_dense(adata.X[ctrl_mask, :])
    hv_idx = np.argsort(xc_full.var(axis=0))[::-1][:n_top_genes]
gene_names = adata.var_names[hv_idx].tolist()    # length = n_top_genes
```

* `n_top_genes` is the only knob; default `50`, test value `30`.
* The fallback branch ranks genes by variance computed **on control cells only**, not on the full pool, so the panel is not biased by any single perturbation.
* If `n_top_genes` exceeds `n_vars`, all genes are kept; if `< 4` genes survive, the step exits early.

### 3. Gene modules — hierarchical clustering on the control correlation

A single set of gene modules is fitted on the **control** correlation
matrix and reused in all panels (so the same gene has the same module
colour everywhere) ([`genenet.py:548-554`](perturbflow/analyzer/genenet.py#L548-L554)):

```python
ctrl_expr = _to_dense(adata.X[ctrl_mask, :][:, hv_idx])      # (n_ctrl_cells, n_genes)
ctrl_corr = np.corrcoef(ctrl_expr.T)                          # (n_genes, n_genes), NaNs → 0, diag = 1
dist      = np.clip(1.0 - np.abs(ctrl_corr), 0, None)         # signed correlations collapsed to magnitude
Z         = linkage(squareform(dist, checks=False),
                    method="average")
cluster_labels = fcluster(Z, t=n_gene_clusters,
                          criterion="maxclust")               # 1-indexed cluster IDs
```

* Distance is `1 − |r|`, so positively- and negatively-correlated genes are placed in the same module (we treat sign as orientation, not similarity).
* Linkage is **average** (UPGMA). `n_gene_clusters` is the only knob; default `5`, test `4`.
* Cluster IDs are written to `csv/genenet_gene_clusters.csv` together with mean expression columns per perturbation and per-perturbation log₂ fold-changes (added inside the for-pert loop).

### 4. Edges — thresholded Pearson correlation per condition

For each condition (control + each top perturbation) we compute a fresh
gene–gene Pearson correlation matrix on the cells of that condition,
then build an undirected `networkx.Graph` ([`genenet.py:104-121`](perturbflow/analyzer/genenet.py#L104-L121)):

```python
def _build_network(corr, gene_names, threshold):
    G = nx.Graph()
    G.add_nodes_from(gene_names)              # always all n_top_genes nodes
    for i in range(n):
        for j in range(i + 1, n):              # upper triangle only — no self-loops
            r = float(corr[i, j])
            if abs(r) >= threshold:            # default threshold = 0.4
                G.add_edge(gene_names[i], gene_names[j],
                           weight=r,            # signed Pearson r
                           positive=(r >= 0))
```

* The edge predicate is **a hard threshold on `|r|`**; we do *not* compute p-values or FDR. The default `0.4` is empirical and is the only knob.
* Edges are undirected; no self-loops.
* Edge attributes: `weight` = signed `r`, `positive` = sign of `r`.
* This is done once for the control cells (`G_ctrl`) and once per perturbation `p` on cells with `obs.perturbation == p` (`G_pert`). The same `gene_names` set guarantees node correspondence.

### 5. Layout — spring layout on the union graph

To keep gene positions comparable across the three panels of a
perturbation (control / perturbation / differential), the layout is
computed **once** on the union of edges and reused
([`genenet.py:325-329`](perturbflow/analyzer/genenet.py#L325-L329)):

```python
G_union = nx.Graph()
G_union.add_nodes_from(gene_names)
G_union.add_edges_from(G_ctrl.edges())
G_union.add_edges_from(G_pert.edges())
pos = nx.spring_layout(G_union,
                       weight="weight",
                       seed=0,
                       k=1.5)
```

* Algorithm: NetworkX's Fruchterman–Reingold force-directed layout.
* `weight="weight"` makes higher-`|r|` edges act as stiffer springs, pulling correlated genes closer together.
* `k=1.5` is the target inter-node distance (larger ⇒ more spread).
* `seed=0` makes the layout reproducible.

### 6. Differential network — edge-set difference

The differential graph is built directly from the two edge sets
([`genenet.py:341-355`](perturbflow/analyzer/genenet.py#L341-L355)):

```python
ctrl_edges = {tuple(sorted((u, v))) for u, v in G_ctrl.edges()}
pert_edges = {tuple(sorted((u, v))) for u, v in G_pert.edges()}

for u, v in ctrl_edges | pert_edges:
    if   (u, v) in ctrl_edges and (u, v) in pert_edges:  status = "shared"
    elif (u, v) in ctrl_edges:                           status = "lost"
    else:                                                status = "gained"
    G_diff.add_edge(u, v, status=status)
```

| status   | meaning                                                              | colour     |
|----------|----------------------------------------------------------------------|------------|
| `shared` | edge present in both control and perturbation                        | grey  `#adb5bd` |
| `lost`   | edge present in control, absent in perturbation                      | red   `#e63946` |
| `gained` | edge absent in control, present in perturbation                      | green `#2dc653` |

This is a set-difference on the *thresholded* graphs — it captures whether
a co-expression relationship crosses the `|r| ≥ corr_threshold` boundary
under perturbation. It is *not* a per-edge significance test; small per-
perturbation cell counts can flip individual edges across the threshold
just from sampling noise.

### 7. Visual encodings (per-panel)

| element        | encodes                                                            | source |
|----------------|--------------------------------------------------------------------|--------|
| node position  | spring layout on `G_union`                                         | shared across all three panels for one perturbation |
| node colour    | gene module (4 modules from §3), `tab10` palette                   | shared across all three panels |
| node outline   | thin dark blue (`#1d3557`) so nodes stay legible on light backgrounds | constant |
| node size      | control panel: `max(700, |mean_ctrl|·900 + 700)`; pert panel: same with `mean_pert`; differential: `max(700, |mean_pert − mean_ctrl|·1100 + 700)` | per-gene mean expression in the relevant condition |
| edge colour (ctrl, pert panels) | red `#e63946` if `r ≥ 0`, blue `#457b9d` if `r < 0` | sign of Pearson `r` |
| edge width     | `max(0.8, |r|·3.0)`                                                | magnitude of `r`     |
| edge colour (diff panel)        | grey / red / green per status (see §6)              | edge-set difference  |

Labels are drawn at font-size 14 bold whenever the gene panel has
≤ 80 genes (always true with the default `n_top_genes ≤ 50`).

### 8. Outputs

For every perturbation `p`, the GRN step writes (under `<output_dir>/`):

```
plots/
  genenet_control_heatmap.png             clustermap of the control correlation matrix,
                                          rows/cols sorted by gene module, with module
                                          sidebars and a "Module N" legend.
  genenet_<p>_ctrl_network.png            10 × 10 square spring-layout graph (control)
  genenet_<p>_pert_network.png            10 × 10 square spring-layout graph (perturbation)
  genenet_<p>_diff_network.png            10 × 10 square spring-layout graph (differential)
  genenet_<p>_heatmap_comparison.png      side-by-side row-z-scored heatmap of mean expression:
                                            Control (z) | <p> (z) | log₂FC, gene rows
                                            sorted by module, module sidebar.

csv/
  genenet_gene_clusters.csv               gene, cluster_id, mean_ctrl,
                                          mean_<p>, log2fc_<p>  (one column block per pert)
  genenet_<p>_network.json                nodes (gene_names + cluster_ids + ctrl/pert/log2fc means),
                                          positions (spring layout), and the lists of
                                          control and perturbation edges with signed `r` —
                                          consumed by the interactive HTML report.
```

The interactive report ([`perturbflow/analyzer/interactive.py`](perturbflow/analyzer/interactive.py)) base64-embeds the three per-perturbation networks side-by-side in row 1 of the *Gene Network* tab, and the `control_heatmap` + `<p>_heatmap_comparison` heatmaps in row 2.

### 9. Parameters and defaults

| Config key (JSON)             | Default | Test | Description |
|-------------------------------|--------:|------:|-------------|
| `genenet_n_top_genes`         |      50 |    30 | Number of HVGs (or top-variance-on-control fallback) included as nodes. |
| `genenet_n_gene_clusters`     |       5 |     4 | Number of gene modules (`fcluster` `maxclust` cut). |
| `genenet_corr_threshold`      |    0.40 |  0.40 | `|r|` threshold for keeping an edge in any condition. |
| `deg_n_top_perturbations`     |       — |     3 | How many perturbations the pipeline picks for GRN analysis (chosen by the `effects` step). |

Cell-count gates (hard-coded): control pool must have ≥ 5 cells; each
perturbation must have ≥ 5 cells, otherwise the perturbation is skipped.

### 10. What this network does and does not capture

* **Captures.** Gene–gene relationships that show up as cell-level
  Pearson co-expression `|r| ≥ corr_threshold` in the cells of one
  condition. Module structure recovered by hierarchical clustering on
  the *control* correlation matrix typically resolves clear biological
  blocks (translation, OxPhos, cell-cycle / E2F targets, etc.) in
  K562/RPE1 essential screens.
* **Does not capture.** Direction of regulation, causality, or
  whether two correlated genes are linked directly or via an
  intermediate. The differential edge classification is a set-difference
  on a thresholded graph, not a statistical test of "edge changed
  significantly", so individual gained/lost edges should be treated as
  hypotheses unless they are corroborated (consistent across replicates,
  enriched within a known module, supported by an external prior).
* **Sensitive to.** Cell counts (smaller per-perturbation samples ⇒
  noisier `r`); HVG selection (only the top `n_top_genes` are visible);
  cell-cycle heterogeneity (often dominates HVG variance and so
  dominates the network); library-size / dropout patterns (can induce
  spurious correlations if the upstream `preprocess` step's
  normalisation is loose).

### 11. Suggested extensions

* **TF–target overlay.** Mark each edge as "supported by a known
  TF→target relationship" using CollecTRI (or DoRothEA) and render
  those edges with a different style. This adds mechanistic
  interpretation without changing the topology. The biological-priors
  catalogue used by the sibling `PerturbVerse_v4` pipeline already
  contains a CollecTRI human edge list and would slot in directly here.
* **Background null.** Compute `r` on shuffled cell labels (or on a
  matched random subsample of cells) and report each edge's empirical
  rank, so the threshold becomes data-driven rather than `0.4` fixed.
* **TF–TF regulatory view.** The
  [`perturbflow/analyzer/regulatory.py`](perturbflow/analyzer/regulatory.py) module
  already produces a directed TF→target network from perturbation
  effect sizes; surfacing its output alongside the co-expression
  network would give a "phenotypic vs mechanistic" pair of views per
  perturbation.

### 12. Reproducibility

* Spring layout uses `seed=0` (fixed).
* Hierarchical clustering uses `linkage(method="average")` and
  `fcluster(criterion="maxclust", t=n_gene_clusters)`, both
  deterministic given the input matrix.
* `np.corrcoef` is deterministic; cells are masked by `obs.perturbation`
  (no random sampling inside the GRN step).
* Output filenames sanitise the perturbation label by replacing
  `"/"`, `" "`, and `"+"` with `"_"`.

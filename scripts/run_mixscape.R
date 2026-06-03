#!/usr/bin/env Rscript
# =============================================================================
# Seurat/Mixscape perturb-seq workflow.
#
# Implements the Seurat Mixscape analysis sequence:
#   1. Standard RNA normalization, variable genes, PCA, UMAP.
#   2. CalcPerturbSig to create a perturbation signature assay (PRTB).
#   3. PRTB PCA/UMAP to visualize local perturbation signatures.
#   4. RunMixscape to classify NT, NP, and KO cells.
#   5. KO/NP summaries, perturbation score plots, posterior plots, heatmaps.
#   6. MixscapeLDA and LDA UMAP for perturbation response visualization.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(reshape2)
  library(jsonlite)
  library(reticulate)
})

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs) || length(lhs) == 0 || (length(lhs) == 1 && is.na(lhs))) rhs else lhs
}

parse_args <- function(argv) {
  out <- list()
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    name <- sub("^--", "", key)
    if (i == length(argv) || startsWith(argv[[i + 1]], "--")) {
      out[[name]] <- TRUE
      i <- i + 1
    } else {
      out[[name]] <- argv[[i + 1]]
      i <- i + 2
    }
  }
  out
}

require_pkg <- function(pkg, reason = pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required R package is not installed: ", reason, call. = FALSE)
  }
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

save_plot <- function(plot, path, width = 7, height = 5) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot = plot, width = width, height = height, dpi = 220, bg = "white")
  path
}

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE)
  path
}

first_present <- function(candidates, choices) {
  hit <- candidates[!is.na(candidates) & candidates %in% choices]
  if (length(hit)) hit[[1]] else NULL
}

classic_theme <- function() {
  theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 12, color = "#203040"),
      axis.title = element_text(size = 13, color = "#203040"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11),
      strip.text = element_text(size = 12, face = "bold")
    )
}

load_object <- function(input) {
  if (dir.exists(input)) {
    counts <- Read10X(input)
    return(CreateSeuratObject(counts = counts, project = basename(normalizePath(input))))
  }
  ext <- tolower(tools::file_ext(input))
  if (ext == "rds") {
    return(readRDS(input))
  }
  if (ext %in% c("h5seurat", "h5seur")) {
    require_pkg("SeuratDisk", "SeuratDisk for h5Seurat input")
    return(SeuratDisk::LoadH5Seurat(input, assays = NULL, reductions = NULL, graphs = FALSE, verbose = FALSE))
  }
  if (ext == "h5ad") {
    require_pkg("reticulate", "reticulate plus Python anndata for h5ad input")
    require_pkg("Matrix", "Matrix for Matrix Market h5ad import")
    pybin <- Sys.getenv("RETICULATE_PYTHON", unset = "")
    if (nzchar(pybin)) reticulate::use_python(pybin, required = FALSE)
    prefix <- file.path(tempdir(), paste0("mixscape_h5ad_", Sys.getpid()))
    mtx_out <- paste0(prefix, ".mtx")
    obs_out <- paste0(prefix, "_obs.csv")
    features_out <- paste0(prefix, "_features.txt")
    barcodes_out <- paste0(prefix, "_barcodes.txt")
    code <- sprintf(
"import anndata
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

h5ad_input = %s
mtx_out = %s
obs_out = %s
features_out = %s
barcodes_out = %s

adata = anndata.read_h5ad(h5ad_input)
if 'counts' in adata.layers:
    X = adata.layers['counts']
    var_names = list(adata.var_names)
elif adata.raw is not None:
    X = adata.raw.X
    var_names = list(adata.raw.var_names)
else:
    X = adata.X
    var_names = list(adata.var_names)
if not sp.issparse(X):
    X = sp.csr_matrix(X)
sio.mmwrite(mtx_out, X.T.tocoo())
pd.Series(var_names).to_csv(features_out, index=False, header=False)
pd.Series(list(adata.obs_names)).to_csv(barcodes_out, index=False, header=False)
obs = adata.obs.copy()
obs.index = list(adata.obs_names)
obs.to_csv(obs_out)
",
      jsonlite::toJSON(input, auto_unbox = TRUE),
      jsonlite::toJSON(mtx_out, auto_unbox = TRUE),
      jsonlite::toJSON(obs_out, auto_unbox = TRUE),
      jsonlite::toJSON(features_out, auto_unbox = TRUE),
      jsonlite::toJSON(barcodes_out, auto_unbox = TRUE)
    )
    reticulate::py_run_string(code)
    counts <- Matrix::readMM(mtx_out)
    features <- readLines(features_out, warn = FALSE)
    barcodes <- readLines(barcodes_out, warn = FALSE)
    rownames(counts) <- make.unique(features)
    colnames(counts) <- make.unique(barcodes)
    obs <- utils::read.csv(obs_out, row.names = 1, check.names = FALSE)
    rownames(obs) <- make.unique(rownames(obs))
    obs <- obs[colnames(counts), , drop = FALSE]
    return(CreateSeuratObject(counts = counts, meta.data = obs, project = tools::file_path_sans_ext(basename(input))))
  }
  stop("Unsupported input type: ", input)
}

coerce_metadata <- function(obj, cfg, perturbation_override = NULL) {
  meta_names <- colnames(obj@meta.data)
  pert_col <- perturbation_override %||% cfg$perturbation_col %||%
    first_present(c("perturbation", "guide", "gRNA", "sgRNA", "crispr", "gene", "target"), meta_names)
  if (is.null(pert_col) || !pert_col %in% meta_names) {
    stop("Could not find perturbation metadata column. Set perturbation_col in config or pass the fourth sbatch argument.")
  }

  target_col <- cfg$target_col %||% first_present(c("target", "gene", "target_gene", "feature"), meta_names)
  guide_col <- cfg$guide_col %||% first_present(c("guide", "gRNA", "sgRNA", "crispr", "guide_id"), meta_names)
  replicate_col <- cfg$replicate_col %||% first_present(c("replicate", "batch", "sample", "orig.ident"), meta_names)
  nt_class <- cfg$nt_class %||% "NT"
  controls <- tolower(as.character(cfg$control_labels %||% c("control", "NT", "non-targeting", "nontargeting")))

  raw_pert <- as.character(obj@meta.data[[pert_col]])
  is_nt <- tolower(raw_pert) %in% controls

  if (!is.null(target_col) && target_col %in% meta_names) {
    target <- as.character(obj@meta.data[[target_col]])
  } else {
    target <- raw_pert
  }
  target[is_nt | is.na(target) | target == ""] <- nt_class

  if (!is.null(guide_col) && guide_col %in% meta_names) {
    guide <- as.character(obj@meta.data[[guide_col]])
  } else {
    guide <- raw_pert
  }
  guide[is_nt | is.na(guide) | guide == ""] <- nt_class

  obj$mixscape_target <- target
  obj$mixscape_guide <- guide
  obj$mixscape_is_nt <- obj$mixscape_target == nt_class

  list(
    object = obj,
    perturbation_col = pert_col,
    target_col = "mixscape_target",
    guide_col = "mixscape_guide",
    replicate_col = replicate_col,
    nt_class = nt_class
  )
}

plot_dim_if_present <- function(obj, group.by, reduction, title, path, point_size, cols = NULL) {
  if (!group.by %in% colnames(obj@meta.data) || !reduction %in% names(obj@reductions)) return(NULL)
  p <- DimPlot(obj, group.by = group.by, reduction = reduction, pt.size = point_size, label = FALSE, repel = TRUE, cols = cols) +
    ggtitle(title) + xlab("UMAP 1") + ylab("UMAP 2") + classic_theme()
  save_plot(p, path, width = 7, height = 5.5)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(args$input) || is.null(args$output)) {
  stop("Usage: run_mixscape.R --input <input> --output <output_dir> [--config config.json] [--perturbation-col col]")
}

cfg <- if (!is.null(args$config)) jsonlite::fromJSON(args$config, simplifyVector = FALSE) else list()
set.seed(as.integer(cfg$seed %||% 1))

out <- normalizePath(args$output, mustWork = FALSE)
plots <- file.path(out, "plots")
tables <- file.path(out, "csv")
dir.create(plots, recursive = TRUE, showWarnings = FALSE)
dir.create(tables, recursive = TRUE, showWarnings = FALSE)

manifest <- list(
  input = normalizePath(args$input, mustWork = FALSE),
  config = cfg,
  plots = list(),
  tables = list(),
  warnings = list()
)

obj <- load_object(args$input)
meta <- coerce_metadata(obj, cfg, args[["perturbation-col"]])
obj <- meta$object

assay <- cfg$assay %||% "RNA"
if (!assay %in% names(obj@assays)) {
  assay <- DefaultAssay(obj)
  manifest$warnings <- c(manifest$warnings, paste("Configured assay was not present; using", assay))
}
DefaultAssay(obj) <- assay

dims <- seq_len(as.integer(cfg$dims %||% 40))
dims <- dims[dims <= max(2, min(ncol(obj) - 1, 50))]
point_size <- as.numeric(cfg$point_size %||% 0.2)

obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = max(dims), verbose = FALSE)
dims <- dims[dims <= ncol(Embeddings(obj, "pca"))]
obj <- RunUMAP(obj, dims = dims, reduction = "pca", reduction.name = "umap", reduction.key = "UMAP_", verbose = FALSE)

manifest$plots <- c(manifest$plots, list(
  rna_umap_target = plot_dim_if_present(obj, meta$target_col, "umap", "RNA UMAP by target", file.path(plots, "seurat_rna_umap_target.png"), point_size),
  rna_umap_guide = plot_dim_if_present(obj, meta$guide_col, "umap", "RNA UMAP by guide", file.path(plots, "seurat_rna_umap_guide.png"), point_size),
  rna_umap_replicate = plot_dim_if_present(obj, meta$replicate_col %||% "", "umap", "RNA UMAP by replicate", file.path(plots, "seurat_rna_umap_replicate.png"), point_size)
))

finish_partial_report <- function(reason) {
  manifest$status <- "partial"
  manifest$warnings <- c(manifest$warnings, reason)
  capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
  plot_files <- sort(list.files(plots, pattern = "\\.png$", full.names = FALSE))
  table_files <- sort(list.files(tables, pattern = "\\.csv$", full.names = FALSE))
  manifest$plots <- plot_files
  manifest$tables <- table_files
  jsonlite::write_json(manifest, file.path(out, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
  plot_cards <- paste(
    sprintf('<figure><img src="plots/%s" alt="%s"><figcaption>%s</figcaption></figure>',
            plot_files, tools::file_path_sans_ext(plot_files), tools::file_path_sans_ext(plot_files)),
    collapse = "\n"
  )
  table_links <- paste(sprintf('<li><a href="csv/%s">%s</a></li>', table_files, table_files), collapse = "\n")
  html <- sprintf(
'<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Seurat Mixscape Report</title>
  <style>
    body{font-family:Arial,Helvetica,sans-serif;max-width:1180px;margin:2rem auto;color:#203040;background:#f6f9fc;}
    h1{border-bottom:2px solid #2c6fbb;padding-bottom:.45rem;color:#1e5f9e;}
    h2{color:#1e5f9e;margin-top:2rem;}
    .notice{background:#fff;border:1px solid #cdd9e5;border-left:4px solid #b58b1d;border-radius:6px;padding:1rem;margin:1rem 0;}
    figure{background:#fff;border:1px solid #cdd9e5;border-radius:6px;padding:1rem;margin:1rem 0;}
    img{max-width:100%%;}
    figcaption{color:#637487;font-size:.9rem;margin-top:.4rem;}
  </style>
</head>
<body>
<h1>Seurat Mixscape Perturb-seq Report</h1>
<div class="notice"><strong>Partial run:</strong> %s</div>
<p>RNA-level Seurat preprocessing and UMAP plots were written. Mixscape requires enough non-targeting and perturbed cells per class for local nearest-neighbor signatures.</p>
<h2>Tables</h2>
<ul>%s</ul>
<h2>Plots</h2>
%s
</body>
</html>',
    reason,
    table_links,
    plot_cards
  )
  writeLines(html, file.path(out, "report.html"))
  message("Partial Mixscape report written to: ", file.path(out, "report.html"))
  quit(save = "no", status = 0)
}

requested_neighbors <- as.integer(cfg$num_neighbors %||% 20)
target_sizes <- table(obj@meta.data[[meta$target_col]])
nt_n <- sum(obj@meta.data[[meta$target_col]] == meta$nt_class)
cell_limit <- min(min(target_sizes), nt_n)
if (cell_limit < 8) {
  finish_partial_report(
    paste("Skipped Mixscape because the smallest target/control class has only", cell_limit, "cells.")
  )
}
small_data_cap <- max(2, floor((cell_limit - 1) / 4))
safe_neighbors <- if (cell_limit <= (requested_neighbors * 2 + 1)) {
  min(requested_neighbors, small_data_cap)
} else {
  requested_neighbors
}
if (safe_neighbors < requested_neighbors) {
  manifest$warnings <- c(
    manifest$warnings,
    paste("Reduced CalcPerturbSig num.neighbors from", requested_neighbors, "to", safe_neighbors, "for available cell counts")
  )
}
message("CalcPerturbSig num.neighbors = ", safe_neighbors, " (requested ", requested_neighbors, ")")

calc_args <- list(
  object = obj,
  assay = assay,
  slot = "data",
  gd.class = meta$target_col,
  nt.cell.class = meta$nt_class,
  reduction = "pca",
  ndims = max(dims),
  num.neighbors = safe_neighbors,
  new.assay.name = "PRTB"
)
if (!is.null(meta$replicate_col) && meta$replicate_col %in% colnames(obj@meta.data)) {
  calc_args$split.by <- meta$replicate_col
}
obj <- tryCatch(
  do.call(CalcPerturbSig, calc_args),
  error = function(e) {
    finish_partial_report(paste("Skipped Mixscape because CalcPerturbSig failed:", conditionMessage(e)))
  }
)

DefaultAssay(obj) <- "PRTB"
VariableFeatures(obj) <- VariableFeatures(obj[[assay]])
obj <- ScaleData(obj, do.scale = FALSE, do.center = TRUE, verbose = FALSE)
obj <- RunPCA(obj, reduction.key = "prtbpca_", reduction.name = "prtbpca", npcs = max(dims), verbose = FALSE)
prtb_dims <- dims[dims <= ncol(Embeddings(obj, "prtbpca"))]
if (!length(prtb_dims)) prtb_dims <- seq_len(ncol(Embeddings(obj, "prtbpca")))
obj <- RunUMAP(obj, dims = prtb_dims, reduction = "prtbpca", reduction.key = "prtbumap_", reduction.name = "prtbumap", verbose = FALSE)

manifest$plots <- c(manifest$plots, list(
  prtb_umap_target = plot_dim_if_present(obj, meta$target_col, "prtbumap", "PRTB UMAP by target", file.path(plots, "mixscape_prtb_umap_target.png"), point_size),
  prtb_umap_guide = plot_dim_if_present(obj, meta$guide_col, "prtbumap", "PRTB UMAP by guide", file.path(plots, "mixscape_prtb_umap_guide.png"), point_size),
  prtb_umap_replicate = plot_dim_if_present(obj, meta$replicate_col %||% "", "prtbumap", "PRTB UMAP by replicate", file.path(plots, "mixscape_prtb_umap_replicate.png"), point_size)
))

obj <- RunMixscape(
  object = obj,
  assay = "PRTB",
  slot = "scale.data",
  labels = meta$target_col,
  nt.class.name = meta$nt_class,
  min.de.genes = as.integer(cfg$mixscape_min_de_genes %||% 5),
  iter.num = as.integer(cfg$mixscape_iter_num %||% 10),
  de.assay = assay,
  verbose = FALSE,
  prtb.type = cfg$mixscape_prtb_type %||% "KO"
)

manifest$tables$cell_metadata <- write_csv(obj@meta.data, file.path(tables, "mixscape_cell_metadata.csv"))

class_col <- "mixscape_class.global"
target_col <- meta$target_col
guide_col <- meta$guide_col
if (!class_col %in% colnames(obj@meta.data)) {
  stop("RunMixscape finished but did not create mixscape_class.global")
}

summary_df <- as.data.frame.matrix(table(obj@meta.data[[target_col]], obj@meta.data[[class_col]]))
summary_df$target <- rownames(summary_df)
summary_df$total_cells <- rowSums(summary_df[, setdiff(colnames(summary_df), "target"), drop = FALSE])
summary_df$ko_fraction <- if ("KO" %in% colnames(summary_df)) summary_df$KO / pmax(summary_df$total_cells, 1) else 0
summary_df$np_fraction <- if ("NP" %in% colnames(summary_df)) summary_df$NP / pmax(summary_df$total_cells, 1) else 0
summary_df$nt_fraction <- if ("NT" %in% colnames(summary_df)) summary_df$NT / pmax(summary_df$total_cells, 1) else 0
summary_df <- summary_df[order(summary_df$ko_fraction, decreasing = TRUE), ]
manifest$tables$summary_by_target <- write_csv(summary_df, file.path(tables, "mixscape_summary_by_target.csv"))

guide_tab <- as.data.frame(prop.table(table(obj@meta.data[[class_col]], obj@meta.data[[guide_col]]), 2))
colnames(guide_tab) <- c("mixscape_class", "guide", "fraction")
guide_tab$target <- obj@meta.data[[target_col]][match(guide_tab$guide, obj@meta.data[[guide_col]])]
manifest$tables$summary_by_guide <- write_csv(guide_tab, file.path(tables, "mixscape_summary_by_guide.csv"))

plot_guides <- guide_tab[guide_tab$guide != meta$nt_class, , drop = FALSE]
plot_guides$mixscape_class <- factor(plot_guides$mixscape_class, levels = c("NT", "NP", "KO"))
p_guide <- ggplot(plot_guides, aes(x = guide, y = fraction * 100, fill = mixscape_class)) +
  geom_col(width = 0.8) +
  facet_wrap(vars(target), scales = "free_x") +
  scale_fill_manual(values = c(NT = "#6b7280", NP = "#b8c3cf", KO = "#b04a5a"), drop = FALSE) +
  ylab("% of cells") + xlab("Guide") + labs(fill = "Mixscape class") +
  classic_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
manifest$plots$ko_np_bar <- save_plot(p_guide, file.path(plots, "mixscape_ko_np_nt_by_guide.png"), width = 12, height = 7)

manifest$plots$mixscape_umap_global <- plot_dim_if_present(
  obj, class_col, "prtbumap", "PRTB UMAP by Mixscape class",
  file.path(plots, "mixscape_prtb_umap_class_global.png"), point_size,
  cols = c(NT = "#6b7280", NP = "#b8c3cf", KO = "#b04a5a")
)

top_targets <- summary_df$target[summary_df$target != meta$nt_class]
top_targets <- head(top_targets, as.integer(cfg$top_targets %||% 12))
heatmap_targets <- head(top_targets, as.integer(cfg$max_heatmap_targets %||% 8))

for (target in top_targets) {
  safe <- safe_name(target)
  perturb_path <- file.path(plots, paste0("mixscape_perturb_score_", safe, ".png"))
  try({
    p <- PlotPerturbScore(
      object = obj,
      target.gene.ident = target,
      mixscape.class = "mixscape_class",
      col = "#b04a5a"
    ) + labs(fill = "Mixscape class") + classic_theme()
    manifest$plots[[paste0("perturb_score_", safe)]] <- save_plot(p, perturb_path, width = 7, height = 5)
  }, silent = TRUE)

  ids <- grep(paste0("^", gsub("([.|()\\^{}+$*?\\[\\]\\\\])", "\\\\\\1", target), " "), unique(obj$mixscape_class), value = TRUE)
  ids <- c(meta$nt_class, ids)
  ids <- ids[ids %in% unique(obj$mixscape_class)]
  if ("mixscape_class_p_ko" %in% colnames(obj@meta.data) && length(ids) >= 2) {
    try({
      p <- VlnPlot(obj, features = "mixscape_class_p_ko", idents = ids, group.by = "mixscape_class", pt.size = 0) +
        ggtitle(paste(target, "posterior probability")) + NoLegend() + classic_theme() +
        theme(axis.text.x = element_text(angle = 20, hjust = 1))
      manifest$plots[[paste0("posterior_", safe)]] <- save_plot(p, file.path(plots, paste0("mixscape_posterior_", safe, ".png")), width = 7, height = 5)
    }, silent = TRUE)
  }
}

Idents(obj) <- target_col
for (target in heatmap_targets) {
  safe <- safe_name(target)
  try({
    p <- MixscapeHeatmap(
      object = obj,
      ident.1 = meta$nt_class,
      ident.2 = target,
      balanced = FALSE,
      assay = assay,
      max.genes = as.integer(cfg$mixscape_heatmap_max_genes %||% 20),
      angle = 0,
      group.by = "mixscape_class",
      max.cells.group = as.integer(cfg$mixscape_heatmap_max_cells_group %||% 300),
      size = 6
    ) + NoLegend() + theme(axis.text.y = element_text(size = 12))
    manifest$plots[[paste0("heatmap_", safe)]] <- save_plot(p, file.path(plots, paste0("mixscape_heatmap_", safe, ".png")), width = 9, height = 7)
  }, silent = TRUE)
}

try({
  Idents(obj) <- class_col
  lda_obj <- subset(obj, idents = c("KO", "NT"))
  if (length(unique(lda_obj@meta.data[[target_col]])) >= 2) {
    lda_obj <- MixscapeLDA(
      object = lda_obj,
      assay = assay,
      pc.assay = "PRTB",
      labels = target_col,
      nt.label = meta$nt_class,
      npcs = as.integer(cfg$lda_npcs %||% 10),
      logfc.threshold = as.numeric(cfg$lda_logfc_threshold %||% 0.25),
      verbose = FALSE
    )
    lda_dims <- seq_len(min(ncol(Embeddings(lda_obj, "lda")), max(1, length(unique(lda_obj@meta.data[[target_col]])) - 1)))
    lda_obj <- RunUMAP(lda_obj, dims = lda_dims, reduction = "lda", reduction.key = "ldaumap_", reduction.name = "ldaumap", verbose = FALSE)
    Idents(lda_obj) <- "mixscape_class"
    p <- DimPlot(lda_obj, reduction = "ldaumap", label = TRUE, repel = TRUE, pt.size = point_size) +
      xlab("UMAP 1") + ylab("UMAP 2") + NoLegend() + classic_theme()
    manifest$plots$lda_umap <- save_plot(p, file.path(plots, "mixscape_lda_umap.png"), width = 8, height = 6)
    if (isTRUE(cfg$save_rds %||% TRUE)) saveRDS(lda_obj, file.path(out, "mixscape_lda_object.rds"))
  }
}, silent = TRUE)

if (isTRUE(cfg$save_rds %||% TRUE)) {
  saveRDS(obj, file.path(out, "mixscape_seurat_object.rds"))
}

capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))

plot_files <- sort(list.files(plots, pattern = "\\.png$", full.names = FALSE))
table_files <- sort(list.files(tables, pattern = "\\.csv$", full.names = FALSE))
manifest$plots <- plot_files
manifest$tables <- table_files
jsonlite::write_json(manifest, file.path(out, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)

plot_cards <- paste(
  sprintf('<figure><img src="plots/%s" alt="%s"><figcaption>%s</figcaption></figure>',
          plot_files, tools::file_path_sans_ext(plot_files), tools::file_path_sans_ext(plot_files)),
  collapse = "\n"
)
table_links <- paste(sprintf('<li><a href="csv/%s">%s</a></li>', table_files, table_files), collapse = "\n")

html <- sprintf(
'<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Seurat Mixscape Report</title>
  <style>
    body{font-family:Arial,Helvetica,sans-serif;max-width:1180px;margin:2rem auto;color:#203040;background:#f6f9fc;}
    h1{border-bottom:2px solid #2c6fbb;padding-bottom:.45rem;color:#1e5f9e;}
    h2{color:#1e5f9e;margin-top:2rem;}
    .summary{display:flex;gap:1rem;flex-wrap:wrap;margin:1rem 0;}
    .card{background:#fff;border:1px solid #cdd9e5;border-radius:6px;padding:1rem;min-width:150px;}
    .value{font-size:1.5rem;font-weight:700;color:#2c6fbb;}
    figure{background:#fff;border:1px solid #cdd9e5;border-radius:6px;padding:1rem;margin:1rem 0;}
    img{max-width:100%%;}
    figcaption{color:#637487;font-size:.9rem;margin-top:.4rem;}
    code{background:#eef4f8;padding:2px 5px;border-radius:3px;}
  </style>
</head>
<body>
<h1>Seurat Mixscape Perturb-seq Report</h1>
<p>Workflow follows Seurat Mixscape: RNA embedding, perturbation signatures, Mixscape KO/NP/NT classification, perturbation score inspection, heatmaps, and LDA.</p>
<div class="summary">
  <div class="card"><div>Cells</div><div class="value">%s</div></div>
  <div class="card"><div>Targets</div><div class="value">%s</div></div>
  <div class="card"><div>Guides</div><div class="value">%s</div></div>
  <div class="card"><div>Top KO target</div><div class="value">%s</div></div>
</div>
<h2>Tables</h2>
<ul>%s</ul>
<h2>Plots</h2>
%s
</body>
</html>',
  format(ncol(obj), big.mark = ","),
  length(unique(obj@meta.data[[target_col]])),
  length(unique(obj@meta.data[[guide_col]])),
  if (nrow(summary_df)) as.character(summary_df$target[[1]]) else "NA",
  table_links,
  plot_cards
)
writeLines(html, file.path(out, "report.html"))

message("Mixscape report written to: ", file.path(out, "report.html"))

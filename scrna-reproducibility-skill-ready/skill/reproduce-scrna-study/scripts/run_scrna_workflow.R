#!/usr/bin/env Rscript
bootstrap_args <- commandArgs(trailingOnly = FALSE)
bootstrap_file <- sub("^--file=", "", grep("^--file=", bootstrap_args, value = TRUE)[[1]])
source(file.path(dirname(normalizePath(bootstrap_file, winslash = "/")), "common.R"))
cli <- parse_cli("all")
if (isTRUE(cli$help)) {
  cat("Usage: Rscript run_scrna_workflow.R --config CONFIG [--stage all|validate|qc|analysis|figures] [--dry-run]\n")
  quit(status = 0)
}
cfg <- load_config(cli$config)
meta <- validate_config(cfg, require_data = !cli$dry_run)
if (cli$stage == "validate" || cli$dry_run) {
  cat(sprintf("Configuration valid: %d samples, %d cohorts.\n", nrow(meta), length(unique(meta$cohort))))
  quit(status = 0)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

out <- cfg$project$output_dir
ensure_dirs(out)
set.seed(cfg$analysis$seed %||% 1L)
options(stringsAsFactors = FALSE, future.globals.maxSize = (cfg$analysis$future_globals_gb %||% 8) * 1024^3)

read_matrix <- function(path) {
  if (dir.exists(path)) return(Read10X(path))
  if (grepl("\\.h5$", path, ignore.case = TRUE)) return(Read10X_h5(path, use.names = TRUE, unique.features = TRUE))
  stop("Unsupported input. Use a 10X H5 file or matrix directory: ", path)
}

run_qc <- function() {
  objects <- vector("list", nrow(meta)); rows <- vector("list", nrow(meta))
  q <- cfg$qc
  for (i in seq_len(nrow(meta))) {
    message(sprintf("[%d/%d] %s", i, nrow(meta), meta$sample_id[[i]]))
    counts <- read_matrix(sample_paths(cfg, meta)[[i]])
    obj <- CreateSeuratObject(counts, project = meta$sample_id[[i]], min.cells = 0, min.features = 0)
    obj <- RenameCells(obj, add.cell.id = meta$sample_id[[i]])
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = q$mitochondrial_pattern %||% "^MT-")
    for (column in setdiff(names(meta), "file_name")) obj[[column]] <- meta[[column]][[i]]
    keep <- obj$nFeature_RNA >= (q$min_features %||% 200) & obj$nFeature_RNA <= (q$max_features %||% Inf) &
      obj$nCount_RNA >= (q$min_counts %||% 0) & obj$nCount_RNA <= (q$max_counts %||% Inf) &
      obj$percent.mt <= (q$max_mito_percent %||% 20)
    rows[[i]] <- data.frame(sample_id = meta$sample_id[[i]], sample_group = meta$sample_group[[i]], cohort = meta$cohort[[i]],
      cells_before = ncol(obj), cells_after = sum(keep), retention_pct = round(100 * mean(keep), 2),
      median_features = median(obj$nFeature_RNA), median_counts = median(obj$nCount_RNA), median_mito = median(obj$percent.mt))
    objects[[i]] <- subset(obj, cells = colnames(obj)[keep])
  }
  names(objects) <- meta$sample_id
  qc <- do.call(rbind, rows)
  write.csv(qc, file.path(out, "tables", "qc_summary.csv"), row.names = FALSE)
  merged <- merge(objects[[1]], y = objects[-1], merge.data = FALSE)
  saveRDS(merged, file.path(out, "objects", "01_qc_merged.rds"), compress = FALSE)
  p <- ggplot(qc, aes(x = reorder(sample_group, retention_pct), y = retention_pct, fill = cohort)) +
    geom_col() + coord_flip() + scale_y_continuous(limits = c(0, 100)) + theme_classic(base_size = 11) +
    labs(x = NULL, y = "Cells retained (%)", title = "QC retention by sample")
  ggsave(file.path(out, "figures", "Figure_01_QC_retention.png"), p, width = 8, height = 6, dpi = 240)
  merged
}

run_analysis <- function(obj) {
  a <- cfg$analysis
  obj <- NormalizeData(obj, scale.factor = a$normalization_scale_factor %||% 10000, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = a$variable_features %||% 2000, verbose = FALSE)
  regress <- unlist(a$regress %||% character())
  regress <- regress[regress %in% colnames(obj[[]])]
  obj <- ScaleData(obj, features = VariableFeatures(obj), vars.to.regress = if (length(regress)) regress else NULL, verbose = FALSE)
  obj <- RunPCA(obj, npcs = a$pcs %||% 30, verbose = FALSE)
  dims <- seq_len(a$graph_dims %||% 20)
  obj <- FindNeighbors(obj, dims = dims, verbose = FALSE)
  obj <- FindClusters(obj, resolution = a$resolution %||% 0.5, random.seed = a$seed %||% 1L, verbose = FALSE)
  reductions <- tolower(unlist(a$reductions %||% c("tsne", "umap")))
  if ("tsne" %in% reductions) obj <- RunTSNE(obj, dims = dims, seed.use = a$seed %||% 1L,
    perplexity = a$tsne_perplexity %||% 30, check_duplicates = FALSE, verbose = FALSE)
  if ("umap" %in% reductions) obj <- RunUMAP(obj, dims = dims, seed.use = a$seed %||% 1L, verbose = FALSE)
  for (nm in names(cfg$gene_sets %||% list())) {
    genes <- intersect(unlist(cfg$gene_sets[[nm]]), rownames(obj))
    if (length(genes)) obj <- AddModuleScore(obj, features = list(genes), name = paste0(nm, "_"),
      ctrl = a$module_control_genes %||% 50, seed = a$seed %||% 1L)
  }
  saveRDS(obj, file.path(out, "objects", "02_analyzed.rds"), compress = FALSE)
  obj
}

make_figures <- function(obj) {
  reductions <- intersect(tolower(unlist(cfg$analysis$reductions)), names(obj@reductions))
  groups <- unlist(cfg$figures$group_by %||% c("sample_group", "seurat_clusters"))
  for (red in reductions) for (group in groups) {
    if (!group %in% colnames(obj[[]])) next
    p <- DimPlot(obj, reduction = red, group.by = group, label = identical(group, "seurat_clusters"), repel = TRUE,
      raster = FALSE, pt.size = cfg$figures$point_size %||% 0.25) + theme_classic(base_size = 11) +
      ggtitle(sprintf("%s by %s", toupper(red), group))
    ggsave(file.path(out, "figures", sprintf("Figure_%s_%s.png", toupper(red), group)), p, width = 9, height = 7, dpi = 260)
  }
  markers <- intersect(unlist(cfg$figures$markers %||% character()), rownames(obj))
  if (length(markers) && length(reductions)) {
    p <- FeaturePlot(obj, reduction = reductions[[1]], features = markers, ncol = min(4, length(markers)), order = TRUE,
      min.cutoff = "q05", max.cutoff = "q95", cols = c("#D9D9D9", "#2C7FB8", "#08306B"))
    ggsave(file.path(out, "figures", "Figure_marker_features.png"), p, width = 13, height = 3.2 * ceiling(length(markers) / 4), dpi = 240)
  }
  score_cols <- grep("_1$", colnames(obj[[]]), value = TRUE)
  if (length(score_cols)) {
    means <- aggregate(obj[[]][score_cols], list(sample_group = obj$sample_group), mean)
    write.csv(means, file.path(out, "tables", "module_score_summary.csv"), row.names = FALSE)
    mat <- scale(as.matrix(means[-1])); mat[!is.finite(mat)] <- 0; rownames(mat) <- means$sample_group
    if (requireNamespace("pheatmap", quietly = TRUE)) {
      png(file.path(out, "figures", "Figure_module_score_heatmap.png"), 1800, 1300, res = 200)
      pheatmap::pheatmap(mat, border_color = NA, main = "Relative module scores (column-wise z-score)")
      dev.off()
    }
  }
}

allowed <- c("all", "qc", "analysis", "figures")
if (!cli$stage %in% allowed) stop("Unknown stage: ", cli$stage)
obj_file_qc <- file.path(out, "objects", "01_qc_merged.rds")
obj_file_analysis <- file.path(out, "objects", "02_analyzed.rds")
obj <- if (cli$stage %in% c("all", "qc")) run_qc() else readRDS(obj_file_qc)
if (cli$stage %in% c("all", "analysis")) obj <- run_analysis(obj)
if (cli$stage %in% c("all", "figures")) {
  if (cli$stage == "figures") obj <- readRDS(obj_file_analysis)
  make_figures(obj)
}
writeLines(capture.output(sessionInfo()), file.path(out, "logs", "sessionInfo.txt"))
file.copy(cfg$._config_path, file.path(out, "logs", "config_used.yaml"), overwrite = TRUE)
message("Completed stage: ", cli$stage, ". Outputs: ", out)

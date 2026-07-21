#!/usr/bin/env Rscript
packages <- c("Seurat", "yaml", "ggplot2", "patchwork", "pheatmap")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
cat("Required packages are available.\n")

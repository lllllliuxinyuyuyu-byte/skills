if (!requireNamespace("Seurat", quietly = TRUE) || !requireNamespace("yaml", quietly = TRUE)) {
  stop("Install Seurat and yaml before running the end-to-end test.")
}
suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = FALSE)
this <- normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[[1]]), winslash = "/")
root <- normalizePath(file.path(dirname(this), ".."), winslash = "/")
fixture <- file.path(tempdir(), "scrna-skill-fixture")
unlink(fixture, recursive = TRUE)
dir.create(fixture, recursive = TRUE)

gzip_file <- function(source, destination) {
  input <- file(source, "rb"); output <- gzfile(destination, "wb")
  on.exit({ close(input); close(output) })
  repeat {
    chunk <- readBin(input, what = "raw", n = 1024 * 1024)
    if (!length(chunk)) break
    writeBin(chunk, output)
  }
}

genes <- c("PECAM1", "CDH5", "KDR", "ACTA2", "MKI67", paste0("GENE", 1:55))
set.seed(42)
for (sample in c("S1", "S2")) {
  d <- file.path(fixture, sample); dir.create(d)
  counts <- rsparsematrix(length(genes), 80, density = 0.35)
  counts@x <- as.numeric(pmax(1, round(abs(counts@x) * 5)))
  writeMM(counts, file.path(d, "matrix.mtx"))
  write.table(data.frame(genes, genes, "Gene Expression"), file.path(d, "features.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(paste0(sample, "_cell", seq_len(ncol(counts))), file.path(d, "barcodes.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  for (name in c("matrix.mtx", "features.tsv", "barcodes.tsv")) {
    gzip_file(file.path(d, name), file.path(d, paste0(name, ".gz")))
    unlink(file.path(d, name))
  }
}

meta <- data.frame(
  sample_id = c("S1", "S2"), file_name = c("S1", "S2"), sample_group = c("D0", "D1"),
  cohort = "test", day = c(0, 1), cell_line = "synthetic", protocol = "test",
  stage = c("early", "late"), sample_kind = "derived"
)
write.csv(meta, file.path(fixture, "metadata.csv"), row.names = FALSE)
cfg <- list(
  project = list(name = "synthetic smoke test", data_dir = fixture, metadata = "metadata.csv", output_dir = "results"),
  qc = list(min_features = 1, max_features = 1000, min_counts = 0, max_counts = 100000, max_mito_percent = 100),
  analysis = list(seed = 42, variable_features = 30, pcs = 10, graph_dims = 5, resolution = 0.3, reductions = c("tsne", "umap"), tsne_perplexity = 10, module_control_genes = 1),
  figures = list(group_by = c("sample_group", "seurat_clusters"), markers = c("PECAM1", "ACTA2", "MKI67")),
  gene_sets = list(),
  assertions = list(list(label = "synthetic cohort", sample_ids = c("S1", "S2"), expected_cells_before = 160)),
  narrative_claims = list(list(label = "synthetic trajectory", sample_ids = c("S1", "S2"), stated_days = c(0, 1)))
)
yaml::write_yaml(cfg, file.path(fixture, "config.yaml"))

runner <- file.path(root, "skill", "reproduce-scrna-study", "scripts", "run_scrna_workflow.R")
auditor <- file.path(root, "skill", "reproduce-scrna-study", "scripts", "audit_consistency.R")
status <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(runner), "--config", shQuote(file.path(fixture, "config.yaml")), "--stage", "all"))
stopifnot(status == 0)
status <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(auditor), "--config", shQuote(file.path(fixture, "config.yaml"))))
stopifnot(status == 0)
expected <- c("objects/02_analyzed.rds", "tables/qc_summary.csv", "tables/consistency_audit.csv", "figures/Figure_TSNE_sample_group.png", "figures/Figure_UMAP_sample_group.png")
stopifnot(all(file.exists(file.path(fixture, "results", expected))))
cat("End-to-end synthetic test passed.\n")

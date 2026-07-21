#!/usr/bin/env Rscript
bootstrap_args <- commandArgs(trailingOnly = FALSE)
bootstrap_file <- sub("^--file=", "", grep("^--file=", bootstrap_args, value = TRUE)[[1]])
source(file.path(dirname(normalizePath(bootstrap_file, winslash = "/")), "common.R"))
cli <- parse_cli(NULL)
if (isTRUE(cli$help)) {
  cat("Usage: Rscript audit_consistency.R --config CONFIG\n")
  quit(status = 0)
}
cfg <- load_config(cli$config)
meta <- validate_config(cfg, require_data = FALSE)
out <- cfg$project$output_dir
ensure_dirs(out)
qc_path <- file.path(out, "tables", "qc_summary.csv")
qc <- if (file.exists(qc_path)) read.csv(qc_path) else NULL
rows <- list(); k <- 0L
add <- function(category, label, status, detail) {
  k <<- k + 1L
  rows[[k]] <<- data.frame(category, label, status, detail, stringsAsFactors = FALSE)
}

for (a in cfg$assertions %||% list()) {
  ids <- unlist(a$sample_ids)
  unknown <- setdiff(ids, meta$sample_id)
  if (length(unknown)) {
    add("cohort", a$label, "FAIL", paste("Unknown sample IDs:", paste(unknown, collapse = ", ")))
    next
  }
  declared_days <- sort(unique(meta$day[meta$sample_id %in% ids & !is.na(meta$day)]))
  if (!is.null(qc) && !is.null(a$expected_cells_before)) {
    observed <- sum(qc$cells_before[qc$sample_id %in% ids])
    status <- if (observed == a$expected_cells_before) "PASS" else "FAIL"
    add("cell_total", a$label, status, sprintf("expected=%s; observed=%s; days=%s", a$expected_cells_before, observed, paste(declared_days, collapse = ",")))
  } else {
    add("cohort", a$label, "INFO", sprintf("samples=%d; days=%s", length(ids), paste(declared_days, collapse = ",")))
  }
}

for (claim in cfg$narrative_claims %||% list()) {
  ids <- unlist(claim$sample_ids)
  actual <- sort(unique(meta$day[meta$sample_id %in% ids & !is.na(meta$day)]))
  stated <- sort(as.numeric(unlist(claim$stated_days)))
  status <- if (length(actual) == length(stated) && all(actual == stated)) "PASS" else "FAIL"
  add("narrative", claim$label, status, sprintf("stated_days=%s; matrix_days=%s; samples=%s", paste(stated, collapse = ","), paste(actual, collapse = ","), paste(ids, collapse = ",")))
}

audit <- if (length(rows)) do.call(rbind, rows) else data.frame(category = character(), label = character(), status = character(), detail = character())
write.csv(audit, file.path(out, "tables", "consistency_audit.csv"), row.names = FALSE)
cat(paste(capture.output(print(audit, row.names = FALSE)), collapse = "\n"), "\n")
if (any(audit$status == "FAIL")) quit(status = 2)

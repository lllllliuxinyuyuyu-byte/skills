`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit)) normalizePath(sub("^--file=", "", hit[[1]]), winslash = "/") else normalizePath(".", winslash = "/")
}

parse_cli <- function(default_stage = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(config = NULL, stage = default_stage, dry_run = FALSE)
  i <- 1L
  while (i <= length(args)) {
    if (args[[i]] == "--config" && i < length(args)) { out$config <- args[[i + 1L]]; i <- i + 2L; next }
    if (args[[i]] == "--stage" && i < length(args)) { out$stage <- args[[i + 1L]]; i <- i + 2L; next }
    if (args[[i]] == "--dry-run") { out$dry_run <- TRUE; i <- i + 1L; next }
    if (args[[i]] %in% c("-h", "--help")) { out$help <- TRUE; i <- i + 1L; next }
    stop("Unknown or incomplete argument: ", args[[i]])
  }
  out
}

expand_env <- function(x) {
  if (!is.character(x)) return(x)
  vapply(x, function(value) {
    repeat {
      m <- regexpr("\\$\\{[A-Za-z_][A-Za-z0-9_]*\\}", value)
      if (m[[1]] < 0) break
      token <- regmatches(value, m)
      key <- substr(token, 3, nchar(token) - 1)
      replacement <- chartr("\\", "/", Sys.getenv(key, unset = token))
      if (identical(replacement, token)) break
      value <- sub("\\$\\{[A-Za-z_][A-Za-z0-9_]*\\}", replacement, value)
    }
    value
  }, character(1), USE.NAMES = FALSE)
}

resolve_path <- function(path, base_dir) {
  path <- expand_env(path)
  if (grepl("^[A-Za-z]:[/\\\\]|^/", path)) normalizePath(path, winslash = "/", mustWork = FALSE)
  else normalizePath(file.path(base_dir, path), winslash = "/", mustWork = FALSE)
}

load_config <- function(path) {
  if (is.null(path)) stop("Provide --config <config.yaml>.")
  if (!requireNamespace("yaml", quietly = TRUE)) stop("Install the R package 'yaml'.")
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  cfg <- yaml::read_yaml(path)
  cfg$._config_path <- path
  cfg$._base_dir <- dirname(path)
  cfg$project$metadata <- resolve_path(cfg$project$metadata, cfg$._base_dir)
  cfg$project$data_dir <- resolve_path(cfg$project$data_dir, cfg$._base_dir)
  cfg$project$output_dir <- resolve_path(cfg$project$output_dir %||% "results", cfg$._base_dir)
  cfg
}

read_metadata <- function(cfg) {
  meta <- read.csv(cfg$project$metadata, na.strings = c("", "NA"), check.names = FALSE)
  required <- c("sample_id", "file_name", "sample_group", "cohort", "day", "cell_line", "protocol", "stage", "sample_kind")
  missing <- setdiff(required, names(meta))
  if (length(missing)) stop("Metadata is missing columns: ", paste(missing, collapse = ", "))
  if (anyDuplicated(meta$sample_id)) stop("sample_id values must be unique.")
  meta
}

sample_paths <- function(cfg, meta) file.path(cfg$project$data_dir, meta$file_name)

ensure_dirs <- function(out) {
  for (d in c(out, file.path(out, "objects"), file.path(out, "figures"), file.path(out, "tables"), file.path(out, "logs"))) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}

write_tsv <- function(x, path) write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")

validate_config <- function(cfg, require_data = TRUE) {
  meta <- read_metadata(cfg)
  problems <- character()
  paths <- sample_paths(cfg, meta)
  if (require_data && any(!file.exists(paths))) problems <- c(problems, paste("Missing matrix:", paths[!file.exists(paths)]))
  reductions <- tolower(unlist(cfg$analysis$reductions %||% c("tsne", "umap")))
  bad_red <- setdiff(reductions, c("tsne", "umap"))
  if (length(bad_red)) problems <- c(problems, paste("Unsupported reduction:", bad_red))
  ids <- meta$sample_id
  for (a in cfg$assertions %||% list()) {
    unknown <- setdiff(unlist(a$sample_ids), ids)
    if (length(unknown)) problems <- c(problems, paste0("Assertion '", a$label, "' contains unknown IDs: ", paste(unknown, collapse = ", ")))
  }
  if (length(problems)) stop(paste(problems, collapse = "\n"))
  invisible(meta)
}

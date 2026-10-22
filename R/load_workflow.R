# Load script and report dependencies independently of the launch directory.

# Entries may sit at the root or in simulation/ and reports/.
slads_workflow_root <- function(input_path) {
  input_path <- normalizePath(input_path, winslash = "/", mustWork = TRUE)
  if (dir.exists(input_path)) stop("input_path must name an entry file.")
  root <- dirname(input_path)
  if (!file.exists(file.path(root, "R", "load_workflow.R")) &&
      basename(root) %in% c("simulation", "reports")) root <- dirname(root)
  if (!file.exists(file.path(root, "R", "load_workflow.R"))) {
    stop("Cannot locate R/load_workflow.R beside this entry or its parent directory.")
  }
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

# Resolve user paths from the repository root; outputs need not exist yet.
slads_resolve_path <- function(path, root = .slads_source_root) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("path must be one nonempty string.")
  }
  path <- path.expand(path)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) path <- file.path(root, path)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

slads_workflow_files <- function(workflow) {
  modules <- switch(workflow,
    tables = c("simulation_core.R", "table_reproduction.R"),
    comparison = c("simulation_core.R", "test_helpers.R", "lc_test.R", "clx_test.R",
      "czzw_test.R", "wzz_test.R", "comparison_runner.R"),
    appendix = "appendix_plots.R",
    stop("Unknown workflow: ", workflow)
  )
  file.path("R", c("load_workflow.R", "entrypoint.R", "runtime.R", modules,
                   "workflows.R"))
}

# Load into the caller's environment; pass a fresh one to isolate multiple workflows.
slads_load_workflow <- function(workflow, input_path, envir = parent.frame()) {
  input_path <- normalizePath(input_path, winslash = "/", mustWork = TRUE)
  root <- slads_workflow_root(input_path)
  files <- slads_workflow_files(workflow)
  paths <- file.path(root, files)
  missing <- files[!file.exists(paths) | dir.exists(paths)]
  if (length(missing)) {
    stop("Missing workflow dependencies: ", paste(missing, collapse = ", "),
      ". Keep R/, simulation/, reports/ and required data/ files in the repository.")
  }
  prefix <- paste0(sub("/+$", "", root), "/")
  if (any(!startsWith(normalizePath(paths, winslash = "/", mustWork = TRUE), prefix))) {
    stop("A workflow dependency resolves outside the repository root.")
  }
  assign(".slads_input_path", input_path, envir = envir)
  assign(".slads_source_root", root, envir = envir)
  assign(".slads_source_files", files, envir = envir)
  # Track the bootstrap without invoking it again; load execution helpers last.
  source_files <- files[grepl("\\.R$", files) & files != "R/entrypoint.R"]
  for (file in source_files) sys.source(file.path(root, file), envir = envir)
  invisible(envir)
}

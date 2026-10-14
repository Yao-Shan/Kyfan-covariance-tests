# Load one report's modules and data relative to its Rmd, independent of launch directory.

slads_workflow_files <- function(workflow) {
  modules <- switch(workflow,
    tables = c("simulation_core.R", "table_reproduction.R"),
    comparison = c("simulation_core.R", "entrywise_tests.R", "wzz_test.R",
      "comparison_runner.R"),
    mean_shift = c("simulation_core.R", "mean_shift_runner.R", "mean_shift_reference.R"),
    appendix = "appendix_plots.R",
    stop("Unknown workflow: ", workflow)
  )
  data <- switch(workflow,
    tables = c("table_archive.csv", "manuscript_reference.csv"),
    mean_shift = "manuscript_reference.csv",
    character()
  )
  c(file.path("R", c("load_workflow.R", "runtime.R", modules)),
    if (length(data)) file.path("data", data))
}

# Load into the caller's environment; pass a fresh one to isolate multiple workflows.
slads_load_workflow <- function(workflow, input_path, envir = parent.frame()) {
  input_path <- normalizePath(input_path, mustWork = TRUE)
  root <- dirname(input_path)
  files <- slads_workflow_files(workflow)
  missing <- files[!file.exists(file.path(root, files))]
  if (length(missing)) {
    stop("Missing workflow dependencies: ", paste(missing, collapse = ", "),
      ". Keep the Rmd, R/ and data/ folders together.")
  }
  assign(".slads_input_path", input_path, envir = envir)
  assign(".slads_source_root", root, envir = envir)
  assign(".slads_source_files", files, envir = envir)
  source_files <- files[grepl("\\.R$", files) & files != "R/load_workflow.R"]
  for (file in source_files) sys.source(file.path(root, file), envir = envir)
  invisible(envir)
}

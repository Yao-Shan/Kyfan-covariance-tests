# Shared orchestration for simulation scripts and optional reports.

slads_tables_workflow <- function(params) {
  stopifnot(params$mode %in% c("rerun", "archived"))
  params$output_root <- slads_resolve_path(params$output_root, .slads_source_root)
  grid <- slads_grid(params$profile)
  reference <- NULL
  reference_info <- NULL
  if (!is.character(params$reference_csv) || length(params$reference_csv) != 1L ||
      is.na(params$reference_csv)) stop("reference_csv must be one path or an empty string.")
  if (nzchar(params$reference_csv)) {
    if (params$mode != "rerun") stop("Reference comparison requires mode='rerun'.")
    reference_path <- normalizePath(slads_resolve_path(params$reference_csv, .slads_source_root), mustWork = TRUE)
    reference_raw <- read.csv(reference_path, stringsAsFactors = FALSE, check.names = FALSE)
    rate_column <- params$reference_rate_column
    keys <- c("p", "m", "covariance", "pair", "k")
    if (!is.character(rate_column) || length(rate_column) != 1L || is.na(rate_column) ||
        !all(c(keys, rate_column) %in% names(reference_raw))) {
      stop("Reference CSV needs p, m, covariance, pair, k and the selected rate column.")
    }
    reference_n <- params$reference_R
    if (!is.numeric(reference_n) || length(reference_n) != 1L || !is.finite(reference_n) ||
        reference_n < 1 || reference_n != floor(reference_n)) {
      stop("reference_R must be the positive integer repetition count of the fixed reference.")
    }
    reference <- reference_raw[keys]
    reference$reference_rate <- reference_raw[[rate_column]]
    if (!nrow(reference) || anyNA(reference) ||
        !all(vapply(reference[c("p", "m", "k")], function(x)
          is.numeric(x) && all(is.finite(x) & x >= 1 & x == floor(x)), logical(1))) ||
        !is.numeric(reference$reference_rate) ||
        any(!is.finite(reference$reference_rate) | reference$reference_rate < 0 |
            reference$reference_rate > 1) || anyDuplicated(archive_key(reference))) {
      stop("Reference keys must be unique and complete; rates must be finite and within [0,1].")
    }
    reference$reference_R <- reference_n
    reference$reference_mcse <- sqrt(reference$reference_rate *
                                      (1 - reference$reference_rate) / reference_n)
    reference_info <- list(path = reference_path, md5 = unname(tools::md5sum(reference_path)),
                           rate_column = rate_column, R = reference_n)
    cat("Fixed reference:", reference_path, "\nRate column:", rate_column,
        "\nReference repetitions:", reference_n, "\n")
  } else {
    cat("No reference supplied. Simulation and table generation require no saved results.\n")
  }
  if (params$mode == "rerun") {
    result <- slads_run(
      grid, R = params$R, B = params$B, alpha = params$alpha, seed = params$seed,
      workers = params$workers, resume = params$resume, output_root = params$output_root,
      protected_paths = if (is.null(reference_info)) character() else reference_info$path
    )
    table_data <- result$summary
    table_data$k <- as.integer(sub("KyFan_k", "", table_data$method, fixed = TRUE))
    table_data$R <- params$R
    table_data$B <- params$B
    result_dir <- result$directory
    provenance <- "Computed simulation results (including matching resumed repetitions, if any)"
  } else {
    if (!params$profile %in% c("manuscript", "archive")) {
      stop("Saved-results viewing supports profile='manuscript' or 'archive'.")
    }
    result <- NULL
    archived <- read_archived_results()
    table_data <- if (params$profile == "archive") archived else select_manuscript(archived)
    archive_path <- file.path(.slads_source_root, "data", "table_archive.csv")
    archive_id <- slads_hash(list(profile = params$profile, md5 = unname(tools::md5sum(archive_path))))
    result_dir <- file.path(params$output_root, paste0("archived_", substr(archive_id, 1, 12)))
    dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(archive_path, file.path(result_dir, "archive_input.csv"), overwrite = TRUE)) {
      stop("Could not preserve the viewed archive input.")
    }
    provenance <- "Saved results only: 360 later rerun main-table values and 660 historical extensions; no new simulation"
  }
  write.csv(table_data, file.path(result_dir, "table_values.csv"), row.names = FALSE)
  write_paper_tex(table_data, file.path(result_dir, "tables.tex"))
  cat("Source:", provenance, "\nOutput directory:", result_dir, "\n")
  diagnostics <- NULL
  comparison <- NULL
  comparison_dir <- NULL
  if (!is.null(result)) {
    diagnostics <- table_data[c("p", "m", "covariance", "pair", "k", "requested", "valid",
                                "failed", "rejections", "rejection_rate", "mcse", "ci_low", "ci_high")]
    diagnostics$alpha <- params$alpha
    diagnostics$rate_minus_alpha <- diagnostics$rejection_rate - params$alpha
    write.csv(diagnostics, file.path(result_dir, "calibration_summary.csv"), row.names = FALSE)
  } else {
    cat("Saved table view: archived estimates used R=B=1000 and alpha=0.05.\n")
  }
  if (!is.null(reference) && !is.null(result)) {
    at <- match(archive_key(table_data), archive_key(reference))
    comparison <- table_data[c("p", "m", "covariance", "pair", "k", "R", "B", "alpha",
                                "rejection_rate", "mcse")]
    names(comparison)[names(comparison) == "rejection_rate"] <- "new_rate"
    names(comparison)[names(comparison) == "mcse"] <- "new_mcse"
    comparison$reference_matched <- !is.na(at)
    comparison$reference_rate <- reference$reference_rate[at]
    comparison$reference_R <- reference$reference_R[at]
    comparison$reference_mcse <- reference$reference_mcse[at]
    comparison$difference <- comparison$new_rate - comparison$reference_rate
    comparison$rate_exact_match <- comparison$new_rate == comparison$reference_rate
    comparison$independent_difference_mcse <- sqrt(comparison$new_mcse^2 + comparison$reference_mcse^2)
    comparison_id <- slads_hash(reference_info)
    comparison_dir <- file.path(result_dir, "reference_comparisons", substr(comparison_id, 1, 12))
    dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)
    if (unname(tools::md5sum(reference_info$path)) != reference_info$md5) {
      stop("Reference changed during the run. Restore the fixed input before comparing.")
    }
    reference_copy <- file.path(comparison_dir, "reference_input.csv")
    if (!file.copy(reference_info$path, reference_copy, overwrite = TRUE) ||
        unname(tools::md5sum(reference_copy)) != reference_info$md5) {
      stop("Could not preserve the fixed reference input.")
    }
    saveRDS(reference_info, file.path(comparison_dir, "reference_info.rds"))
    writeLines(capture.output(str(reference_info)), file.path(comparison_dir, "reference_info.txt"))
    write.csv(comparison, file.path(comparison_dir, "new_vs_reference.csv"), row.names = FALSE)
    cat("Reference matches:", sum(comparison$reference_matched), "of", nrow(comparison),
        "new result cells.\nComparison directory:", comparison_dir, "\n")
  } else {
    cat("Reference comparison skipped; no fixed reference was supplied for a simulation run.\n")
  }
  list(grid = grid, result = result, table_data = table_data, result_dir = result_dir,
    provenance = provenance, diagnostics = diagnostics, reference = reference,
    reference_info = reference_info, comparison = comparison, comparison_dir = comparison_dir)
}

slads_comparison_workflow <- function(params) {
  params$output_root <- slads_resolve_path(params$output_root, .slads_source_root)
  grid <- slads_grid(params$profile)
  result <- slads_run(
    grid, R = params$R, B = params$B, alpha = params$alpha, seed = params$seed,
    workers = params$workers, output_root = params$output_root,
    resume = params$resume, checkpoint_every = params$checkpoint_every,
    input_path = .slads_input_path
  )
  cat("Output directory:", result$directory, "\n")
  list(grid = grid, result = result)
}

slads_appendix_workflow <- function(params) {
  plot_spec <- appendix_plot_spec()
  plot_paths <- appendix_plot_paths(params, .slads_source_root)
  summary_data <- read.csv(plot_paths$input_csv, stringsAsFactors = FALSE, check.names = FALSE)
  checked <- appendix_validate_summary(summary_data, plot_spec)
  output_files <- appendix_export_pdfs(checked$data, plot_paths$figure_dir, plot_spec)
  cat("Validated", nrow(checked$data), "cells; wrote eight vector PDFs to:",
      plot_paths$figure_dir, "\n")
  list(plot_spec = plot_spec, plot_paths = plot_paths, checked = checked,
    output_files = output_files)
}

slads_execute_workflow <- function(workflow, params) {
  switch(workflow,
    tables = slads_tables_workflow(params),
    comparison = slads_comparison_workflow(params),
    appendix = slads_appendix_workflow(params),
    stop("Unknown workflow: ", workflow)
  )
}

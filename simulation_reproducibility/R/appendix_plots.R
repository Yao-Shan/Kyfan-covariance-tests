# Appendix C: validate archived counts and draw Wilson-interval PDFs; no simulation.

# Fixed Appendix C design, method aliases, labels, and axes.
appendix_plot_spec <- function() {
  method_order <- c("KyFan_k1", "KyFan_k3", "LC", "CLX", "CZZW", "WZZ")
  method_titles <- c("Ky-Fan (k = 1)", "Ky-Fan (k = 3)", "LC", "CLX",
                     "CZZW", "WZZ")
  method_description <- c("Ky-Fan, order 1", "Ky-Fan, order 3",
                          "Li and Chen (2012)", "Cai, Liu and Xia (2013)",
                          "Chang, Zhou, Zhou and Wang (2017)",
                          "Wang, Zhu and Zhang (2024)")
  list(
    required_columns = c("m", "p", "covariance", "pair", "method", "alpha",
                         "rejections", "valid", "rejection_rate", "ci_low",
                         "ci_high", "confidence_level"),
    method_aliases = c(KyFan_k1 = "KyFan_k1", KyFan_k3 = "KyFan_k3",
                       LC = "LC", CLX = "CLX", Chang = "CZZW", CZZW = "CZZW",
                       Wang = "WZZ", WZZ = "WZZ"),
    method_order = method_order,
    method_titles = method_titles,
    method_table = data.frame(Method = method_titles,
                              Reference = method_description),
    pair_order = c("N-N", "T-T", "E-E", "G-G", "L-L", "N-T", "N-E", "N-G",
                   "N-L", "T-E", "T-G", "T-L", "E-G", "E-L", "G-L"),
    dimensions = c(10L, 40L),
    covariances = paste0("C", 1:4),
    axis_max = setNames(c(0.12, 0.12, 0.28, 0.12, 0.12, 0.12), method_order)
  )
}

# Resolve relative paths against base; expand absolute and tilde-prefixed paths.
appendix_absolute_path <- function(path, base, must_work = FALSE) {
  path <- path.expand(path)
  if (!grepl("^(/|[A-Za-z]:[/\\\\]|~)", path)) path <- file.path(base, path)
  normalizePath(path, winslash = "/", mustWork = must_work)
}

# Resolve input/output paths from the notebook location; require an existing input.
appendix_plot_paths <- function(params, source_root) {
  if (is.null(params$input_csv)) {
    params$input_csv <- "../manuscript/supplementary/comparison_summary.csv"
  }
  if (is.null(params$figure_dir)) params$figure_dir <- "../manuscript/figures"
  list(input_csv = appendix_absolute_path(params$input_csv, source_root,
                                         must_work = TRUE),
       figure_dir = appendix_absolute_path(params$figure_dir, source_root))
}

# Require finite nonnegative integer counts (strictly positive when requested).
appendix_assert_integer_count <- function(x, name, positive = FALSE) {
  if (any(!is.finite(x)) || any(x != floor(x)) ||
      any(x < if (positive) 1 else 0)) {
    stop(name, " must contain finite ", if (positive) "positive" else "nonnegative",
         " integer counts.")
  }
  invisible(NULL)
}

# Recompute rates, 95% Wilson intervals, and MCSE using valid counts as denominators.
appendix_wilson_from_counts <- function(rejections, valid) {
  z <- qnorm(0.975)
  rate_recomputed <- rejections / valid
  denominator <- 1 + z^2 / valid
  center <- (rate_recomputed + z^2 / (2 * valid)) / denominator
  half_width <- z * sqrt(rate_recomputed * (1 - rate_recomputed) /
                          valid + z^2 / (4 * valid^2)) / denominator
  data.frame(rate_recomputed = rate_recomputed,
             lower_recomputed = pmax(0, center - half_width),
             upper_recomputed = pmin(1, center + half_width),
             mcse_recomputed = sqrt(rate_recomputed * (1 - rate_recomputed) / valid))
}

# Check saved summaries against counts; skip absent optional columns.
# Stop on non-finite values or discrepancies above tolerance.
appendix_check_saved <- function(data, column, calculated, tolerance = 1e-10) {
  if (!column %in% names(data)) {
    return(data.frame(column = character(), max_absolute_difference = numeric()))
  }
  saved <- data[[column]]
  if (any(!is.finite(saved))) stop("Non-finite saved values in ", column, ".")
  discrepancy <- max(abs(saved - calculated))
  if (discrepancy > tolerance) {
    stop("Saved ", column, " disagrees with counts; max difference = ",
         format(discrepancy, digits = 6), ".")
  }
  data.frame(column = column, max_absolute_difference = discrepancy)
}

# Build row keys from the ordered design identifiers.
appendix_cell_keys <- function(data, key_columns) {
  do.call(paste, c(data[key_columns], sep = "|"))
}

# Validate all 720 design cells, counts, summaries, and axes before plotting.
# Also check optional MCSE and conditional summaries when present.
appendix_validate_summary <- function(summary_data, spec = appendix_plot_spec()) {
  missing_columns <- setdiff(spec$required_columns, names(summary_data))
  if (length(missing_columns)) {
    stop("Missing required CSV columns: ", paste(missing_columns, collapse = ", "))
  }
  # Normalize aliases before validation so mixed aliases cannot hide duplicates.
  selected <- with(summary_data,
                   m == 150 & p %in% spec$dimensions &
                   covariance %in% spec$covariances & pair %in% spec$pair_order &
                   method %in% names(spec$method_aliases))
  if (anyNA(selected)) stop("Missing design identifiers in the input CSV.")
  data <- summary_data[selected, , drop = FALSE]
  data$method <- unname(spec$method_aliases[data$method])

  key_columns <- c("p", "covariance", "pair", "method")
  expected <- expand.grid(p = spec$dimensions, covariance = spec$covariances,
                          pair = spec$pair_order, method = spec$method_order,
                          stringsAsFactors = FALSE)
  if (anyDuplicated(appendix_cell_keys(data, key_columns))) stop("Duplicate selected design cells.")
  if (nrow(data) != nrow(expected) ||
      !setequal(appendix_cell_keys(data, key_columns),
                appendix_cell_keys(expected, key_columns))) {
    stop("Selected data must contain exactly 720 design cells (90 per PDF).")
  }
  if (any(!is.finite(data$alpha)) || any(abs(data$alpha - 0.05) > 1e-12) ||
      any(!is.finite(data$confidence_level)) ||
      any(abs(data$confidence_level - 0.95) > 1e-12)) {
    stop("The Appendix C figures require alpha = 0.05 and confidence level = 0.95.")
  }
  appendix_assert_integer_count(data$valid, "valid", positive = TRUE)
  appendix_assert_integer_count(data$rejections, "rejections")
  if (any(data$rejections > data$valid)) stop("Rejections exceed valid replications.")
  if (any(data$valid != 1000L)) stop("Appendix C requires 1,000 valid replications.")
  if ("failed" %in% names(data) && any(is.na(data$failed) | data$failed != 0)) {
    stop("The archived Appendix C subset must have no failed replications.")
  }
  if ("complete" %in% names(data) && any(is.na(data$complete) | !data$complete)) {
    stop("The archived Appendix C subset must contain only complete cells.")
  }
  if ("B" %in% names(data) && any(is.na(data$B) | data$B != 1000L)) {
    stop("Appendix C requires the saved B = 1,000 design.")
  }

  calculated <- appendix_wilson_from_counts(data$rejections, data$valid)
  data$rate_recomputed <- calculated$rate_recomputed
  data$lower_recomputed <- calculated$lower_recomputed
  data$upper_recomputed <- calculated$upper_recomputed
  summary_columns <- list(
    rate_recomputed = c("rejection_rate", "type1_error_rate",
                        "conditional_rejection_rate"),
    lower_recomputed = c("ci_low", "conditional_ci_low"),
    upper_recomputed = c("ci_high", "conditional_ci_high"),
    mcse_recomputed = c("mcse", "conditional_mcse")
  )
  validation <- data.frame(column = character(), max_absolute_difference = numeric())
  for (quantity in names(summary_columns)) {
    for (column in summary_columns[[quantity]]) {
      validation <- rbind(validation,
                          appendix_check_saved(data, column, calculated[[quantity]]))
    }
  }
  # Refuse to silently clip an interval against the fixed Appendix C axes.
  if (any(data$upper_recomputed > spec$axis_max[data$method] + 1e-12)) {
    stop("An interval exceeds the fixed Appendix C axis limits; inspect the input.")
  }
  list(data = data, validation = validation)
}

# Draw one six-panel PDF from validated data; destination directory must exist.
appendix_draw_wilson_figure <- function(data, dimension, covariance, destination,
                                       spec = appendix_plot_spec()) {
  pair_order <- spec$pair_order
  method_order <- spec$method_order
  method_titles <- spec$method_titles
  axis_max <- spec$axis_max
  grDevices::pdf(destination, width = 160 / 25.4, height = 160 / 25.4,
                family = "Helvetica", pointsize = 9, useDingbats = FALSE,
                title = sprintf("Wilson intervals: p=%d, m=150, %s",
                                dimension, covariance))
  on.exit(grDevices::dev.off(), add = TRUE)
  colors <- c("#21618C", "#D35420")
  y <- 16 - seq_along(pair_order)
  row_colors <- colors[ifelse(seq_along(pair_order) <= 5, 1L, 2L)]
  # Physical panel locations reproduce the original square six-panel layout.
  panel_left <- c(0.061, 0.548)
  panel_right <- c(0.498, 0.986)
  panel_bottom <- c(0.747, 0.456, 0.165)
  panel_top <- c(0.967, 0.676, 0.385)
  for (i in seq_along(method_order)) {
    row <- (i - 1L) %/% 2L + 1L
    column <- (i - 1L) %% 2L + 1L
    method <- method_order[i]
    panel <- data[data$p == dimension & data$covariance == covariance &
                    data$method == method, , drop = FALSE]
    # Restore the prescribed row order independently of the CSV row order.
    panel <- panel[match(pair_order, panel$pair), , drop = FALSE]
    xmax <- axis_max[method]
    ticks <- seq(0, xmax, length.out = if (method == "LC") 5L else 4L)
    par(fig = c(panel_left[column], panel_right[column],
                panel_bottom[row], panel_top[row]),
        mar = c(0, 0, 0, 0), new = i > 1L, xaxs = "i", yaxs = "i",
        mgp = c(0, 0.15, 0), tcl = -0.15, lend = "butt")
    plot.new()
    plot.window(xlim = c(-0.0025 / 0.12 * xmax, xmax * 1.025),
                ylim = c(0.4, 15.6))
    abline(v = ticks, col = "#E0E0E0", lwd = 0.45)
    abline(h = y, col = "#E8E8E8", lwd = 0.4)
    abline(v = 0.05, col = "#333333", lty = "44", lwd = 0.8)
    abline(h = 10.5, col = "#999999", lty = "33", lwd = 0.55)
    segments(panel$lower_recomputed, y, panel$upper_recomputed, y,
             col = row_colors, lwd = 0.7)
    segments(panel$lower_recomputed, y - 0.18, panel$lower_recomputed, y + 0.18,
             col = row_colors, lwd = 0.45)
    segments(panel$upper_recomputed, y - 0.18, panel$upper_recomputed, y + 0.18,
             col = row_colors, lwd = 0.45)
    points(panel$rate_recomputed, y, pch = 16, col = row_colors, cex = 0.66)
    axis(1, at = ticks, labels = sprintf("%.2f", ticks), cex.axis = 0.78,
         col = "#666666", col.axis = "#222222", lwd = 0.45)
    if (column == 1L) {
      axis(2, at = y, labels = pair_order, las = 1, tick = FALSE,
           cex.axis = 0.78, col.axis = "#222222", hadj = 1, line = -0.2)
    }
    box(bty = "l", col = "#888888", lwd = 0.45)
    mtext(sprintf("(%s) %s", letters[i], method_titles[i]), side = 3,
          adj = 0, line = 0.12, cex = 0.96, font = 2)
  }
  # Separate page-level labels keep all panels the same size.
  par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE,
      xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(0.54, 0.112, "Empirical Type I error rate", cex = 0.90)
  legend_x <- c(0.064, 0.290, 0.539)
  legend_colors <- c(colors, "#333333")
  segments(legend_x, 0.066, legend_x + 0.030, 0.066,
           col = legend_colors, lwd = 0.8)
  points(legend_x + 0.015, rep(0.066, 3L), pch = 16,
         col = legend_colors, cex = 0.66)
  segments(c(0.539, 0.569), 0.062, c(0.539, 0.569), 0.070,
           col = "#333333", lwd = 0.5)
  text(legend_x + 0.042, rep(0.066, 3L),
       c("Same distribution", "Different distributions", "95% Wilson CI"),
       adj = 0, cex = 0.78)
  segments(0.786, 0.056, 0.786, 0.076, col = "#333333", lty = "44", lwd = 0.8)
  text(0.804, 0.066, "Nominal: 0.05", adj = 0, cex = 0.78)
  text(0.54, 0.029,
       "Axes differ by method; fixed across all eight figures.  m = 150, R = B = 1,000.",
       cex = 0.74)
  invisible(destination)
}

# Export all eight PDFs from validated data, replacing existing files with fixed names.
appendix_export_pdfs <- function(data, figure_dir, spec = appendix_plot_spec()) {
  if (!dir.exists(figure_dir)) dir.create(figure_dir, recursive = TRUE)
  if (!dir.exists(figure_dir)) stop("Could not create figure directory: ", figure_dir)
  output_files <- character()
  for (dimension in spec$dimensions) {
    for (covariance in spec$covariances) {
      destination <- file.path(figure_dir,
                               sprintf("wilson_p%d_%s.pdf", dimension, covariance))
      appendix_draw_wilson_figure(data, dimension, covariance, destination, spec)
      output_files <- c(output_files, destination)
    }
  }
  stopifnot(length(output_files) == 8L, all(file.exists(output_files)),
            all(file.info(output_files)$size > 0))
  output_files
}

# Mean-shift experiment: equal covariance is retained under every translation.
# Load through R/load_workflow.R together with R/simulation_core.R.
# ---- Scenario and mean design ------------------------------------------------

ms_pair_levels <- c("N-N", "T-T", "E-E", "G-G", "L-L", "N-T", "N-E", "N-G", "N-L", "T-E", "T-G",
  "T-L", "E-G", "E-L", "G-L")
# Build the manuscript/smoke grid with the original scenario identifiers and RNG keys.
ms_grid <- function(profile = c("manuscript", "smoke")) {
  profile <- match.arg(profile)
  sizes <- data.frame(p = c(10L, 10L, 40L), m = c(50L, 150L, 150L))
  g <- do.call(rbind, lapply(seq_len(nrow(sizes)), function(i) {
    z <- expand.grid(covariance = paste0("C", 1:4), pair = ms_pair_levels, stringsAsFactors = FALSE)
    z$p <- sizes$p[i]
    z$m <- sizes$m[i]
    z
  }))
  g$law_x <- substr(g$pair, 1, 1)
  g$law_y <- substr(g$pair, 3, 3)
  g$original_law_equal <- g$law_x == g$law_y
  g$scenario <- with(g, paste0("p", p, "_m", m, "_", covariance, "_", pair))
  if (profile == "smoke") {
    ids <- c("p10_m50_C1_N-N", "p10_m50_C2_N-T", "p10_m50_C3_E-E", "p10_m50_C4_N-L", "p10_m150_C2_G-L",
      "p40_m150_C4_N-N")
    g <- g[match(ids, g$scenario), , drop = FALSE]
  }
  rownames(g) <- NULL
  g[c("scenario", "p", "m", "covariance", "pair", "law_x", "law_y", "original_law_equal")]
}

# Validate a scalar integer against inclusive bounds.
ms_integer <- function(x, name, minimum = 1L, maximum = .Machine$integer.max) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x != floor(x) || x < minimum ||
    x > maximum)
    stop(name, " must be an integer in [", minimum, ", ", maximum, "].")
  as.integer(x)
}

# Put the paired baseline first, followed by common and unequal mean shifts.
# mean_deltas sets Euclidean lengths; common_shift applies to every coordinate.
ms_mean_design <- function(mean_deltas = c(0.5, 1, 2), mean_directions = c("sparse", "dense"),
  common_shift = 1) {
  if (!is.numeric(mean_deltas) || !length(mean_deltas) || any(!is.finite(mean_deltas)) || any(mean_deltas <=
    0) || anyDuplicated(mean_deltas))
    stop("mean_deltas must contain distinct, positive, finite numbers.")
  if (!is.character(mean_directions) || !length(mean_directions) || anyNA(mean_directions) ||
    any(!mean_directions %in% c("sparse", "dense")) || anyDuplicated(mean_directions))
    stop("mean_directions must be a unique subset of sparse and dense.")
  if (length(common_shift) != 1L || !is.numeric(common_shift) || !is.finite(common_shift) ||
    common_shift == 0)
    stop("common_shift must be a nonzero finite scalar.")
  fixed <- data.frame(mean_id = c("baseline", "common_shift"), mean_type = c("baseline", "common"),
    direction = c("none", "all_coordinates"), delta = 0, mean_equal = TRUE, common_shift = c(0,
      common_shift), mu_x_definition = c("0", "common_shift * 1_p"), mu_y_definition = c("0",
      "common_shift * 1_p"), stringsAsFactors = FALSE)
  extra <- do.call(rbind, lapply(mean_directions, function(direction) {
    data.frame(mean_id = paste0("unequal_", direction, "_d", seq_along(mean_deltas)), mean_type = "unequal",
      direction = direction, delta = mean_deltas, mean_equal = FALSE, common_shift = 0,
      mu_x_definition = "0", mu_y_definition = if (direction == "sparse") "delta * e_1" else "delta * 1_p / sqrt(p)", stringsAsFactors = FALSE)
  }))
  out <- rbind(fixed, extra)
  rownames(out) <- NULL
  out
}

# Expand one mean-design row into length-p mean vectors; covariance is unchanged.
ms_vectors <- function(p, design_row) {
  mu_x <- rep(design_row$common_shift, p)
  mu_y <- mu_x
  if (design_row$mean_type == "unequal") {
    direction <- if (design_row$direction == "sparse") c(1, rep(0, p - 1L)) else rep(1/sqrt(p), p)
    mu_y <- mu_x + design_row$delta * direction
  }
  list(mu_x = mu_x, mu_y = mu_y)
}

# Tabulate mean vectors and differences for each distinct dimension.
ms_vector_table <- function(grid, design) {
  do.call(rbind, lapply(sort(unique(grid$p)), function(p) do.call(rbind, lapply(seq_len(nrow(design)),
    function(j) {
      mu <- ms_vectors(p, design[j, , drop = FALSE])
      data.frame(p = p, mean_id = design$mean_id[j], coordinate = seq_len(p), mu_x = mu$mu_x,
        mu_y = mu$mu_y, mean_difference = mu$mu_y - mu$mu_x)
    }))))
}

# ---- Persistence and scalar summaries ----------------------------------------

# Hash scientific settings using RDS version 2; the caller excludes worker count.
ms_hash <- function(x) {
  f <- tempfile()
  on.exit(unlink(f))
  saveRDS(x, f, version = 2)
  unname(tools::md5sum(f))
}

# Commit an RDS checkpoint by renaming a temporary file beside its target.
ms_write_rds <- function(x, path) {
  tmp <- tempfile(pattern = ".checkpoint-", tmpdir = dirname(path))
  on.exit(unlink(tmp))
  saveRDS(x, tmp, version = 2)
  if (!file.rename(tmp, path))
    stop("Cannot commit checkpoint: ", path)
}

# Commit a UTF-8 CSV by renaming a temporary file beside its target.
ms_write_csv <- function(x, path) {
  tmp <- tempfile(pattern = ".csv-", tmpdir = dirname(path))
  on.exit(unlink(tmp))
  write.csv(x, tmp, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  if (!file.rename(tmp, path))
    stop("Cannot commit CSV: ", path)
}

# Wilson interval using valid repetitions; return NA bounds when none are valid.
ms_wilson <- function(s, n, level = 0.95) {
  if (n < 1L)
    return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - level)/2)
  ph <- s/n
  denominator <- 1 + z * z/n
  middle <- (ph + z * z/(2 * n))/denominator
  half <- z * sqrt(ph * (1 - ph)/n + z * z/(4 * n * n))/denominator
  c(max(0, middle - half), min(1, middle + half))
}

# Largest finite diagnostic discrepancy, or NA if none exists.
ms_max_finite <- function(x) if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_
# ---- Coupled simulation and statistical summaries ----------------------------

# Reuse zero-mean samples and reset the same permutation stream for every shift.
# Common shifts preserve all statistics up to roundoff; unequal shifts preserve
# only the observed covariance difference.
ms_compare_one <- function(X0, Y0, g, r, config) {
  rows <- vector("list", nrow(config$mean_design))
  baseline <- NULL
  for (j in seq_len(nrow(config$mean_design))) {
    d <- config$mean_design[j, , drop = FALSE]
    mu <- ms_vectors(g$p, d)
    # Shift the original samples; separate group-centering would remove the mean gap.
    X <- sweep(X0, 2L, mu$mu_x, "+")
    Y <- sweep(Y0, 2L, mu$mu_y, "+")
    slads_set_seed(slads_seed(config$seed, g$scenario, r, "permutation"))
    started <- proc.time()[["elapsed"]]
    ans <- tryCatch(kyfan_permutation(X, Y, ks = config$ks, B = config$B, keep_spectra = d$mean_type %in%
      c("baseline", "common")), error = function(e) e)
    seconds <- proc.time()[["elapsed"]] - started
    good <- !inherits(ans, "error") && length(ans$p_value) == length(config$ks) && all(is.finite(ans$p_value)) &&
      all(ans$p_value >= 0 & ans$p_value <= 1) && all(is.finite(ans$statistic)) && all(is.finite(ans$perm))
    status <- if (good) "ok" else "failed"
    detail <- if (good) "" else if (inherits(ans, "error"))
      conditionMessage(ans) else "Non-finite or malformed permutation result."
    if (d$mean_id == "baseline" && good)
      baseline <- ans
    z <- data.frame(
      g[rep(1L, length(config$ks)), , drop = FALSE],
      d[rep(1L, length(config$ks)), , drop = FALSE],
      k = config$ks, method = paste0("KyFan_k", config$ks), replicate = r,
      exchangeable = g$original_law_equal && d$mean_equal,
      mean_difference_norm = sqrt(sum((mu$mu_y - mu$mu_x)^2)),
      mu_x_norm = sqrt(sum(mu$mu_x^2)), mu_y_norm = sqrt(sum(mu$mu_y^2)),
      data_seed = slads_seed(config$seed, g$scenario, r, "data"),
      calibration_seed = slads_seed(config$seed, g$scenario, r, "permutation"),
      statistic = if (good) unname(ans$statistic) else NA_real_,
      p_value = if (good) unname(ans$p_value) else NA_real_,
      reject = if (good) unname(ans$p_value <= config$alpha) else NA,
      augmented_critical_value = NA_real_, observed_abs_diff_from_baseline = NA_real_,
      common_p_value_mismatch = NA, common_reject_mismatch = NA,
      common_perm_statistic_max_abs_diff = NA_real_,
      common_perm_kyfan_curve_max_abs_diff = NA_real_, status = status, detail = detail,
      seconds = seconds, stringsAsFactors = FALSE
    )
    if (good)
      for (h in seq_along(config$ks)) {
        z$augmented_critical_value[h] <- unname(quantile(c(ans$statistic[h], ans$perm[,
          h]), probs = 1 - config$alpha, type = 1L, names = FALSE))
      }
    if (good && !is.null(baseline)) {
      z$observed_abs_diff_from_baseline <- abs(unname(ans$statistic - baseline$statistic))
      if (d$mean_type %in% c("baseline", "common")) {
        z$common_p_value_mismatch <- unname(ans$p_value != baseline$p_value)
        z$common_reject_mismatch <- unname((ans$p_value <= config$alpha) != (baseline$p_value <=
          config$alpha))
        z$common_perm_statistic_max_abs_diff <- apply(abs(ans$perm - baseline$perm),
          2L, max)
        # These are cumulative Ky--Fan curves over ALL orders, scaled by sqrt(m).
        z$common_perm_kyfan_curve_max_abs_diff <- sqrt(g$m) * max(abs(ans$perm_spectra -
          baseline$perm_spectra))
      }
    }
    rows[[j]] <- z
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

ms_meta_columns <- c("scenario", "p", "m", "covariance", "pair", "law_x", "law_y", "original_law_equal",
  "mean_id", "mean_type", "direction", "delta", "mean_equal", "common_shift", "mu_x_definition",
  "mu_y_definition", "k", "method", "exchangeable", "mean_difference_norm", "mu_x_norm", "mu_y_norm")
# Main rates and intervals require all R trials to be valid; report valid-only estimates.
# Failed trials never count as non-rejections.
ms_summarize <- function(raw, config) {
  split_rows <- split(seq_len(nrow(raw)), interaction(raw$mean_id, raw$k, drop = TRUE))
  out <- do.call(rbind, lapply(split_rows, function(ii) {
    d <- raw[ii, , drop = FALSE]
    ok <- d$status == "ok" & !is.na(d$reject)
    n <- sum(ok)
    rejected <- sum(d$reject[ok])
    failed <- nrow(d) - n
    complete <- nrow(d) == config$R && failed == 0L
    rate <- if (n > 0L) rejected/n else NA_real_
    se <- if (n > 0L) sqrt(rate * (1 - rate)/n) else NA_real_
    ci <- ms_wilson(rejected, n)
    data.frame(
      d[1, ms_meta_columns, drop = FALSE], R = config$R, B = config$B, alpha = config$alpha,
      seed = config$seed, attempted = nrow(d), valid = n, failed = failed, rejections = rejected,
      complete = complete,
      type1_error_rate = if (complete) rate else NA_real_,
      mcse = if (complete) se else NA_real_,
      ci_low = if (complete) ci[1] else NA_real_,
      ci_high = if (complete) ci[2] else NA_real_, confidence_level = 0.95,
      conditional_type1_error_rate = rate, conditional_mcse = se,
      conditional_ci_low = ci[1], conditional_ci_high = ci[2],
      mean_augmented_critical_value = if (n > 0L) mean(d$augmented_critical_value[ok]) else NA_real_,
      mean_seconds = mean(d$seconds), stringsAsFactors = FALSE
    )
  }))
  out <- out[order(match(out$mean_id, config$mean_design$mean_id), out$k), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Match replicate indices within each k to compare shifts with the coupled baseline.
# Paired MCSE = sd(shift - baseline)/sqrt(number of valid pairs).
ms_paired <- function(raw, config) {
  base <- raw[raw$mean_id == "baseline", , drop = FALSE]
  groups <- split(seq_len(nrow(raw)), interaction(raw$mean_id, raw$k, drop = TRUE))
  out <- do.call(rbind, lapply(groups, function(ii) {
    d <- raw[ii, , drop = FALSE]
    b <- base[base$k == d$k[1], , drop = FALSE]
    b <- b[match(d$replicate, b$replicate), , drop = FALSE]
    ok <- d$status == "ok" & b$status == "ok" & !is.na(d$reject) & !is.na(b$reject)
    n <- sum(ok)
    diff <- as.integer(d$reject[ok]) - as.integer(b$reject[ok])
    estimate <- if (n > 0L) mean(diff) else NA_real_
    se <- if (n > 1L) sd(diff)/sqrt(n) else NA_real_
    complete <- nrow(d) == config$R && n == config$R
    data.frame(
      d[1, ms_meta_columns, drop = FALSE], R = config$R, B = config$B, alpha = config$alpha,
      matched_valid = n, unmatched_or_failed = nrow(d) - n, complete = complete,
      paired_rate_difference = if (complete) estimate else NA_real_,
      paired_mcse = if (complete) se else NA_real_,
      conditional_paired_rate_difference = estimate, conditional_paired_mcse = se,
      reject_only_shifted = sum(diff == 1L), reject_only_baseline = sum(diff == -1L),
      paired_baseline_rate = if (n > 0L) mean(b$reject[ok]) else NA_real_,
      paired_shifted_rate = if (n > 0L) mean(d$reject[ok]) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  out <- out[order(match(out$mean_id, config$mean_design$mean_id), out$k), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Summarize observed and common-shift permutation discrepancies by mean setting and k.
ms_diagnostics <- function(raw, config) {
  groups <- split(seq_len(nrow(raw)), interaction(raw$mean_id, raw$k, drop = TRUE))
  out <- do.call(rbind, lapply(groups, function(ii) {
    d <- raw[ii, , drop = FALSE]
    is_common <- d$mean_type[1] %in% c("baseline", "common")
    data.frame(
      d[1, ms_meta_columns, drop = FALSE], attempted = nrow(d),
      observed_comparisons = sum(is.finite(d$observed_abs_diff_from_baseline)),
      max_observed_statistic_abs_diff = ms_max_finite(d$observed_abs_diff_from_baseline),
      common_permutation_comparisons = if (is_common)
        sum(!is.na(d$common_p_value_mismatch)) else NA_integer_,
      common_p_value_mismatches = if (is_common)
        sum(d$common_p_value_mismatch, na.rm = TRUE) else NA_integer_,
      common_reject_mismatches = if (is_common)
        sum(d$common_reject_mismatch, na.rm = TRUE) else NA_integer_,
      max_common_permutation_statistic_abs_diff = ms_max_finite(d$common_perm_statistic_max_abs_diff),
      max_common_permutation_kyfan_curve_abs_diff = ms_max_finite(d$common_perm_kyfan_curve_max_abs_diff),
      stringsAsFactors = FALSE
    )
  }))
  out <- out[order(match(out$mean_id, config$mean_design$mean_id), out$k), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# ---- Checkpoints, scenario workers, and experiment orchestration -------------

# Validate checkpoint settings, mean/order rows, RNG seeds, and decisions before resuming.
ms_validate_checkpoint <- function(cp, g, config, signature, path) {
  bad <- function(message) stop("Invalid checkpoint ", path, ": ", message)
  if (!is.list(cp) || !identical(cp$signature, signature) || !identical(cp$scenario, g$scenario))
    bad("configuration/scenario mismatch")
  if (!is.list(cp$repetitions) || length(cp$repetitions) > config$R)
    bad("invalid repetition list")
  expected_mean <- rep(config$mean_design$mean_id, each = length(config$ks))
  expected_k <- rep(config$ks, times = nrow(config$mean_design))
  for (r in seq_along(cp$repetitions)) {
    d <- cp$repetitions[[r]]
    if (!is.data.frame(d) || !all(c(ms_meta_columns, "replicate", "data_seed", "calibration_seed",
      "p_value", "statistic", "reject", "status", "detail", "observed_abs_diff_from_baseline",
      "augmented_critical_value", "seconds", "common_p_value_mismatch", "common_reject_mismatch",
      "common_perm_statistic_max_abs_diff", "common_perm_kyfan_curve_max_abs_diff") %in%
      names(d)))
      bad(paste("missing fields at replicate", r))
    if (!identical(d$mean_id, expected_mean) || !identical(d$k, expected_k))
      bad("invalid mean/order rows")
    if (anyNA(d[c("replicate", "data_seed", "calibration_seed", "status")]) || any(d$replicate !=
      r) || any(d$data_seed != slads_seed(config$seed, g$scenario, r, "data")) || any(d$calibration_seed !=
      slads_seed(config$seed, g$scenario, r, "permutation")))
      bad("invalid replicate or random seeds")
    for (nm in names(g)) if (anyNA(d[[nm]]) || any(d[[nm]] != g[[nm]]))
      bad(paste("invalid field", nm))
    md <- config$mean_design[match(d$mean_id, config$mean_design$mean_id), , drop = FALSE]
    for (nm in names(md)) if (anyNA(d[[nm]]) || any(d[[nm]] != md[[nm]]))
      bad(paste("invalid mean field", nm))
    if (anyNA(d$exchangeable) || any(d$exchangeable != (d$original_law_equal & d$mean_equal)))
      bad("invalid exchangeability flag")
    if (any(!d$status %in% c("ok", "failed")))
      bad("invalid status")
    ok <- d$status == "ok"
    if (any(ok & (!is.finite(d$p_value) | d$p_value < 0 | d$p_value > 1 | !is.finite(d$statistic) |
      is.na(d$reject))) || any(!ok & !is.na(d$reject)) || any(d$reject[ok] != (d$p_value[ok] <=
      config$alpha)))
      bad("inconsistent p-values/rejections/status")
  }
  invisible(TRUE)
}

# Run/resume one scenario and save its checkpoint and summaries.
# Data seeds depend on master seed, scenario, repetition, and data tag.
ms_worker <- function(i, grid, config, signature, directory, resume, checkpoint_every) {
  g <- grid[i, , drop = FALSE]
  path <- file.path(directory, "checkpoints", paste0(g$scenario, ".rds"))
  cp <- list(signature = signature, scenario = g$scenario, repetitions = list())
  if (resume && file.exists(path)) {
    cp <- tryCatch(readRDS(path), error = function(e) stop("Cannot read ", path, ": ", conditionMessage(e)))
    ms_validate_checkpoint(cp, g, config, signature, path)
  }
  done <- cp$repetitions
  root <- slads_sqrt(slads_covariance(g$p, g$covariance))
  if (length(done) < config$R)
    for (r in seq.int(length(done) + 1L, config$R)) {
      slads_set_seed(slads_seed(config$seed, g$scenario, r, "data"))
      X0 <- slads_draw(g$m, g$law_x, root)
      Y0 <- slads_draw(g$m, g$law_y, root)
      done[[r]] <- ms_compare_one(X0, Y0, g, r, config)
      if (r%%checkpoint_every == 0L || r == config$R) {
        cp$repetitions <- done
        ms_write_rds(cp, path)
      }
    }
  raw <- do.call(rbind, done)
  rownames(raw) <- NULL
  summary <- ms_summarize(raw, config)
  paired <- ms_paired(raw, config)
  diagnostics <- ms_diagnostics(raw, config)
  ms_write_csv(summary, file.path(directory, "scenario_summaries", paste0(g$scenario, ".csv")))
  message(g$scenario, ": ", config$R, " repetitions complete; failed evaluations = ", sum(summary$failed))
  list(summary = summary, paired = paired, diagnostics = diagnostics, checkpoint = path)
}

# Export repetitions and failures one scenario at a time to bound memory.
ms_stream <- function(paths, directory) {
  files <- file.path(directory, c("replicates.csv", "failures.csv"))
  tmp <- vapply(seq_along(files), function(i) tempfile(pattern = ".replicates-", tmpdir = directory),
    character(1))
  on.exit(unlink(tmp))
  first <- TRUE
  for (path in paths) {
    d <- do.call(rbind, readRDS(path)$repetitions)
    for (j in seq_along(files)) {
      z <- if (j == 1L) d else d[d$status != "ok", , drop = FALSE]
      write.table(z, tmp[j], sep = ",", row.names = FALSE, col.names = first, append = !first,
        quote = TRUE, na = "", fileEncoding = "UTF-8")
    }
    first <- FALSE
  }
  for (j in seq_along(files)) if (!file.rename(tmp[j], files[j]))
    stop("Cannot commit ", files[j])
}

# Run the mean-shift experiment with checkpoints and source snapshots.
# Settings and loaded sources identify the run; worker count and execution order
# do not alter data or permutation seeds.
ms_run <- function(grid = ms_grid(), R = 1000L, B = 1000L, alpha = 0.05, seed = 20260925L, workers = 7L,
  mean_deltas = c(0.5, 1, 2), mean_directions = c("sparse", "dense"), common_shift = 1, output_root = "../simulation_results/mean_shift",
  resume = TRUE, checkpoint_every = 10L, input_path) {
  R <- ms_integer(R, "R")
  B <- ms_integer(B, "B", minimum = 19L)
  seed <- ms_integer(seed, "seed", minimum = 0L, maximum = 2147483646L)
  workers <- ms_integer(workers, "workers", maximum = 7L)
  checkpoint_every <- ms_integer(checkpoint_every, "checkpoint_every")
  design <- ms_mean_design(mean_deltas, mean_directions, common_shift)
  if (length(alpha) != 1L || !is.numeric(alpha) || !is.finite(alpha) || alpha <= 0 || alpha >=
    1)
    stop("alpha must be strictly between zero and one.")
  if (length(resume) != 1L || !is.logical(resume) || is.na(resume))
    stop("resume must be TRUE or FALSE.")
  if (length(output_root) != 1L || !is.character(output_root) || !nzchar(output_root))
    stop("output_root must be a path.")
  if (missing(input_path) || length(input_path) != 1L || !file.exists(input_path))
    stop("input_path must name this Rmd.")
  input_path <- normalizePath(input_path, mustWork = TRUE)
  expected_grid <- ms_grid("manuscript")
  if (!is.data.frame(grid) || !nrow(grid) || anyDuplicated(grid$scenario) || !identical(names(grid),
    names(expected_grid)) || anyNA(grid) || !all(grid$scenario %in% expected_grid$scenario))
    stop("grid must be a nonempty subset of ms_grid().")
  target <- expected_grid[match(grid$scenario, expected_grid$scenario), , drop = FALSE]
  rownames(target) <- NULL
  rownames(grid) <- NULL
  if (!isTRUE(all.equal(grid, target, check.attributes = FALSE)))
    stop("grid fields do not match manuscript scenarios.")
  if (1/(B + 1) > alpha)
    warning("The plus-one test cannot reject at this alpha with this B.")
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
    BLIS_NUM_THREADS = "1", RCPP_PARALLEL_NUM_THREADS = "1")
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
  }
  # Hash the report, loaded R modules, and external reference data together.
  config <- list(engine = "kyfan-mean-shift-v1", grid = grid, mean_design = design, R = R,
    B = B, alpha = as.numeric(alpha), seed = seed, ks = c(1L, 3L), code_hash = slads_source_hash(input_path),
    R_version = as.character(getRversion()), platform = R.version$platform, rng = "Mersenne-Twister/Inversion/Rejection")
  signature <- ms_hash(config)
  run_id <- paste0("means_", substr(signature, 1, 16))
  if (!resume)
    run_id <- paste0(run_id, "_fresh_", basename(tempfile(pattern = "run_")))
  directory <- file.path(output_root, run_id)
  for (path in c(directory, file.path(directory, "checkpoints"), file.path(directory, "scenario_summaries"))) if (!dir.exists(path) &&
    !dir.create(path, recursive = TRUE, showWarnings = FALSE))
    stop("Cannot create ", path)
  directory <- normalizePath(directory, mustWork = TRUE)
  cfgpath <- file.path(directory, "configuration.rds")
  if (file.exists(cfgpath)) {
    if (!identical(readRDS(cfgpath), config))
      stop("Existing output configuration mismatch: ", cfgpath)
  } else ms_write_rds(config, cfgpath)
  slads_snapshot_sources(directory, input_path)
  ms_write_csv(grid, file.path(directory, "scenario_grid.csv"))
  ms_write_csv(design, file.path(directory, "mean_design.csv"))
  ms_write_csv(ms_vector_table(grid, design), file.path(directory, "mean_vectors.csv"))
  writeLines(capture.output(str(config)), file.path(directory, "configuration.txt"))
  actual_workers <- min(workers, nrow(grid))
  message("Output: ", directory, "\nBase scenarios: ", nrow(grid), "; mean configurations: ",
    nrow(design), "; R=", R, "; B=", B, "; workers=", actual_workers, " (one BLAS thread per worker)")
  common <- list(grid = grid, config = config, signature = signature, directory = directory,
    resume = resume, checkpoint_every = checkpoint_every)
  if (actual_workers > 1L) {
    cl <- parallel::makePSOCKcluster(actual_workers, outfile = file.path(directory, "workers.log"))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    env <- environment(ms_run)
    candidates <- ls(envir = env, all.names = TRUE)
    exported <- candidates[grepl("^(ms_|slads_|kyfan_)", candidates)]
    parallel::clusterExport(cl, exported, envir = env)
    answer <- do.call(parallel::parLapplyLB, c(list(cl = cl, X = seq_len(nrow(grid)), fun = ms_worker),
      common))
  } else answer <- do.call(lapply, c(list(X = seq_len(nrow(grid)), FUN = ms_worker), common))
  collect <- function(name) {
    z <- do.call(rbind, lapply(answer, `[[`, name))
    rownames(z) <- NULL
    z$run_id <- run_id
    z
  }
  summary <- collect("summary")
  paired <- collect("paired")
  diagnostics <- collect("diagnostics")
  expected <- nrow(grid) * nrow(design) * length(config$ks)
  if (nrow(summary) != expected || anyDuplicated(paste(summary$scenario, summary$mean_id, summary$k)))
    stop("Unexpected number of mean/order result cells.")
  ms_write_csv(summary, file.path(directory, "summary.csv"))
  ms_write_csv(paired, file.path(directory, "paired_vs_baseline.csv"))
  ms_write_csv(diagnostics, file.path(directory, "translation_diagnostics.csv"))
  ms_write_csv(summary[summary$failed > 0L, , drop = FALSE], file.path(directory, "failed_cells.csv"))
  ms_stream(vapply(answer, `[[`, character(1), "checkpoint"), directory)
  writeLines(capture.output(sessionInfo()), file.path(directory, "sessionInfo.txt"))
  writeLines(c("Completed requested repetitions. Inspect failures before interpreting rates.",
    paste("run_id:", run_id), paste("failed_evaluations:", sum(summary$failed))), file.path(directory,
    "COMPLETED.txt"))
  if (any(summary$failed > 0L))
    warning("Failed evaluations: affected main rates and intervals are NA; inspect failures.csv.")
  list(summary = summary, paired = paired, diagnostics = diagnostics, config = config, directory = directory,
    workers = actual_workers)
}

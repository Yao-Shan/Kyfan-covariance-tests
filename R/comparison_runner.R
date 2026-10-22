# Deterministic scenario/method streams and checkpointed comparison output.
# Load dependencies together through R/load_workflow.R.

# Method metadata and stable scenario design ---------------------------------

slads_default_methods <- c("KyFan_k1", "KyFan_k3", "LC", "CLX", "CZZW", "WZZ")
# CZZW/WZZ retain historical Chang/Wang RNG tags to preserve random streams.
slads_method_seed_token <- function(method) {
  switch(method, CZZW = "Chang", WZZ = "Wang", method)
}

slads_method_labels <- c(
  KyFan_k1 = "Ky-Fan (k=1)", KyFan_k3 = "Ky-Fan (k=3)",
  LC = "Li-Chen (2012)", CLX = "Cai-Liu-Xia (2013)",
  CZZW = "CZZW (2017)", WZZ = "WZZ (2024), normal-reference"
)
slads_calibrations <- c(
  KyFan_k1 = "balanced permutation; plus-one upper-tail p",
  KyFan_k3 = "balanced permutation; plus-one upper-tail p",
  LC = "upper standard-normal tail", CLX = "asymptotic Gumbel tail",
  CZZW = "Gaussian multiplier bootstrap; empirical upper tail count/B",
  WZZ = "three-cumulant matched shifted/scaled chi-square; 2022 author manuscript"
)
pair_levels <- c(
  "N-N", "T-T", "E-E", "G-G", "L-L", "N-T", "N-E", "N-G", "N-L", "T-E", "T-G",
  "T-L", "E-G", "E-L", "G-L"
)

# Build the ordered H0 grid (or six-scenario smoke subset) with stable RNG keys.
slads_grid <- function(profile = c("manuscript", "smoke")) {
  profile <- match.arg(profile)
  sizes <- data.frame(p = c(10L, 10L, 40L), m = c(50L, 150L, 150L))
  g <- do.call(rbind, lapply(seq_len(nrow(sizes)), function(i) {
    z <- expand.grid(covariance = paste0("C", 1:4), pair = pair_levels, stringsAsFactors = FALSE)
    z$p <- sizes$p[i]
    z$m <- sizes$m[i]
    z
  }))
  g$law_x <- substr(g$pair, 1, 1)
  g$law_y <- substr(g$pair, 3, 3)
  g$exchangeable <- g$law_x == g$law_y
  g$scenario <- with(g, paste0("p", p, "_m", m, "_", covariance, "_", pair))
  if (profile == "smoke") {
    ids <- c(
      "p10_m50_C1_N-N", "p10_m50_C2_N-T", "p10_m50_C3_E-E",
      "p10_m50_C4_N-L", "p10_m150_C2_G-L", "p40_m150_C4_N-N"
    )
    g <- g[match(ids, g$scenario), , drop = FALSE]
  }
  rownames(g) <- NULL
  stopifnot(!anyNA(g), !anyDuplicated(g$scenario))
  g[c("scenario", "p", "m", "covariance", "pair", "law_x", "law_y", "exchangeable")]
}

# Record method/scenario interpretation limits without skipping the test.
slads_applicability <- function(method, g) {
  laws <- c(g$law_x, g$law_y)
  if (method %in% c("KyFan_k1", "KyFan_k3"))
    return(if (g$exchangeable) "equal-law exchangeability control" else
      "covariance-only null; finite-sample calibration/stress evaluation")
  if (method == "LC" && any(laws %in% c("T", "E", "G")))
    return("robustness: common-scale dependence is outside the independent-component model; t8 also lacks eighth moments")
  if (method == "CZZW" && "T" %in% laws)
    return("robustness: t8 lacks finite eighth moments required by the polynomial-tail sufficient condition")
  if (method == "WZZ" && "T" %in% laws)
    return("robustness: t8 lacks the raw eighth moments required by the quadratic-feature condition")
  if (method == "CLX")
    return("finite-sample calibration; tail, dependence, fourth-moment and dimension assumptions are not certified by a finite design")
  "finite-sample calibration; remaining moment, spectral and dimension assumptions are not automatically certified"
}

# Method execution and result records ---------------------------------------

# Normalize replicate results; failed or invalid p-values have missing decisions.
slads_result_row <- function(method, res, seconds, g, alpha) {
  status <- if (is.null(res$status)) "invalid_result" else as.character(res$status)[1]
  pv <- if (is.null(res$p_value)) NA_real_ else as.numeric(res$p_value)[1]
  valid_p <- length(pv) == 1L && is.finite(pv) && pv >= 0 && pv <= 1
  if (status == "ok" && !valid_p)
    status <- "invalid_result"
  decision <- if (status == "ok") pv <= alpha else NA
  detail <- if (is.character(res$detail)) paste(res$detail, collapse = "; ") else ""
  data.frame(
    method = method, p_value = pv,
    statistic = if (is.null(res$statistic)) NA_real_ else as.numeric(res$statistic)[1],
    reject = decision, status = status, detail = detail, seconds = as.numeric(seconds),
    calibration = unname(slads_calibrations[method]),
    applicability = slads_applicability(method, g), stringsAsFactors = FALSE
  )
}

# Apply all six methods to shared samples with deterministic calibration seeds.
slads_compare_one <- function(X, Y, g, r, seed, B, alpha) {
  # Both fixed orders use the same observations and permutation assignments.
  slads_set_seed(slads_seed(seed, g$scenario, r, "permutation"))
  started <- proc.time()[["elapsed"]]
  ky <- tryCatch(
    kyfan_permutation(X, Y, ks = c(1L, 3L), B = B),
    error = function(e) list(error = conditionMessage(e))
  )
  ky_seconds <- proc.time()[["elapsed"]] - started
  rows <- lapply(slads_default_methods, function(method) {
    slads_set_seed(slads_seed(seed, g$scenario, r, slads_method_seed_token(method)))
    started <- proc.time()[["elapsed"]]
    res <- tryCatch({
      if (method %in% c("KyFan_k1", "KyFan_k3")) {
        if (!is.null(ky$error))
          stop(ky$error)
        k <- if (method == "KyFan_k1") "k1" else "k3"
        list(p_value = unname(ky$p_value[k]), statistic = unname(ky$statistic[k]), status = "ok")
      } else switch(method,
        LC = method_li(X, Y), CLX = method_cai(X, Y),
        CZZW = method_chang(X, Y, B = B, mc_plus_one = FALSE), WZZ = method_wang(X, Y)
      )
    }, error = function(e) list(p_value = NA_real_, statistic = NA_real_, status = "error",
      detail = conditionMessage(e)))
    seconds <- proc.time()[["elapsed"]] - started
    if (method %in% c("KyFan_k1", "KyFan_k3"))
      seconds <- seconds + ky_seconds
    slads_result_row(method, res, seconds, g, alpha)
  })
  rows <- do.call(rbind, rows)
  for (nm in names(g)) rows[[nm]] <- g[[nm]]
  rows$replicate <- as.integer(r)
  rows$data_seed <- slads_seed(seed, g$scenario, r, "data")
  rows$calibration_seed <- vapply(
    rows$method,
    function(method) slads_seed(
      seed, g$scenario, r,
      if (method %in% c("KyFan_k1", "KyFan_k3")) "permutation" else
        slads_method_seed_token(method)
    ),
    integer(1)
  )
  rows
}

# Per-scenario summaries and Monte Carlo uncertainty ------------------------

# Two-sided Wilson interval; return missing bounds when the denominator is zero.
slads_wilson <- function(s, n, level = 0.95) {
  if (n < 1L)
    return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - level)/2)
  ph <- s/n
  denominator <- 1 + z * z/n
  middle <- (ph + z * z/(2 * n))/denominator
  half <- z * sqrt(ph * (1 - ph)/n + z * z/(4 * n * n))/denominator
  c(max(0, middle - half), min(1, middle + half))
}

# Main rates require all requested trials to be valid.
# Conditional rates and intervals use only successful trials.
slads_summary <- function(raw, config) {
  groups <- split(raw, interaction(raw$scenario, raw$method, drop = TRUE, lex.order = TRUE))
  out <- do.call(rbind, lapply(groups, function(d) {
    ok <- d$status == "ok" & !is.na(d$reject)
    valid <- sum(ok)
    rejected <- sum(d$reject[ok])
    failed <- nrow(d) - valid
    complete <- nrow(d) == config$R && failed == 0L
    rate <- if (valid > 0L) rejected/valid else NA_real_
    mcse <- if (valid > 0L) sqrt(rate * (1 - rate)/valid) else NA_real_
    ci <- slads_wilson(rejected, valid)
    data.frame(
      d[1, c(
        "scenario", "p", "m", "covariance", "pair", "law_x", "law_y",
        "exchangeable", "method"
      )],
      method_label = unname(slads_method_labels[d$method[1]]),
      k = if (d$method[1] == "KyFan_k1") 1L else
        if (d$method[1] == "KyFan_k3") 3L else NA_integer_,
      R = config$R, B = config$B, alpha = config$alpha, seed = config$seed,
      requested = config$R, attempted = nrow(d), valid = valid,
      failed = failed, rejections = rejected,
      complete = complete,
      type1_error_rate = if (complete) rate else NA_real_,
      rejection_rate = if (complete) rate else NA_real_,
      mcse = if (complete) mcse else NA_real_,
      ci_low = if (complete) ci[1] else NA_real_,
      ci_high = if (complete) ci[2] else NA_real_, confidence_level = 0.95,
      conditional_rejection_rate = rate, conditional_mcse = mcse,
      conditional_ci_low = ci[1], conditional_ci_high = ci[2], calibration = d$calibration[1],
      applicability = d$applicability[1], mean_seconds = mean(d$seconds), stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))
  out <- out[
    order(match(out$scenario, config$grid$scenario),
      match(out$method, slads_default_methods)),
    , drop = FALSE
  ]
  rownames(out) <- NULL
  out
}

# Configuration identity and atomic output ----------------------------------

# Validate a scalar integer against inclusive bounds.
slads_integer <- function(x, name, minimum = 1L, maximum = .Machine$integer.max) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x != floor(x) || x < minimum ||
    x > maximum)
    stop(name, " must be an integer in [", minimum, ", ", maximum, "].")
  as.integer(x)
}

# Hash an object using version-2 RDS serialization.
slads_hash <- function(x) {
  f <- tempfile()
  on.exit(unlink(f))
  saveRDS(x, f, version = 2)
  unname(tools::md5sum(f))
}

# Hash the notebook and loaded dependencies through the shared provenance helper.
slads_code_hash <- function(path) {
  slads_source_hash(path)
}

# Commit an RDS file by renaming a temporary file beside its target.
slads_write_rds <- function(x, path) {
  tmp <- tempfile(pattern = ".checkpoint-", tmpdir = dirname(path))
  on.exit(unlink(tmp))
  saveRDS(x, tmp, version = 2)
  if (!file.rename(tmp, path))
    stop("Cannot commit checkpoint: ", path)
}

# Commit a UTF-8 CSV by renaming a temporary file beside its target.
slads_write_csv <- function(x, path) {
  tmp <- tempfile(pattern = ".csv-", tmpdir = dirname(path))
  on.exit(unlink(tmp))
  write.csv(x, tmp, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  if (!file.rename(tmp, path))
    stop("Cannot commit CSV: ", path)
}

# Checkpoint validation and scenario workers --------------------------------

# Validate checkpoint identity, method rows, seeds, and decisions before resuming.
slads_validate_checkpoint <- function(cp, g, config, signature, path) {
  bad <- function(message) stop("Invalid checkpoint ", path, ": ", message)
  if (!is.list(cp) || !identical(cp$signature, signature) || !identical(cp$scenario, g$scenario))
    bad("configuration or scenario mismatch; do not mix output directories.")
  if (!is.list(cp$repetitions) || length(cp$repetitions) > config$R)
    bad("invalid repetition count.")
  for (r in seq_along(cp$repetitions)) {
    d <- cp$repetitions[[r]]
    required <- c("method", "replicate", "scenario", "data_seed", "calibration_seed", "p_value",
      "reject", "status", "statistic", "detail", "seconds", "calibration", "applicability",
      names(g))
    if (!is.data.frame(d) || !all(required %in% names(d)) ||
        nrow(d) != length(slads_default_methods) ||
        !identical(d$method, slads_default_methods))
      bad(paste("incomplete method rows at repetition", r))
    if (anyNA(d[c("replicate", "data_seed", "calibration_seed", "status")]) || any(d$replicate !=
      r) || any(d$data_seed != slads_seed(config$seed, g$scenario, r, "data")))
      bad(paste("repetition or seed mismatch at", r))
    for (nm in names(g)) if (anyNA(d[[nm]]) || any(d[[nm]] != g[[nm]]))
      bad(paste("scenario field", nm, "mismatch"))
    expected <- vapply(
      d$method,
      function(method) slads_seed(
        config$seed, g$scenario, r,
        if (method %in% c("KyFan_k1", "KyFan_k3")) "permutation" else
          slads_method_seed_token(method)
      ),
      integer(1)
    )
    if (!identical(unname(d$calibration_seed), unname(expected)))
      bad("calibration seed mismatch")
    ok <- d$status == "ok"
    if (any(ok & (!is.finite(d$p_value) | d$p_value < 0 | d$p_value > 1 | is.na(d$reject))) ||
      any(!ok & !is.na(d$reject)) || any(d$reject[ok] != (d$p_value[ok] <= config$alpha)))
      bad("inconsistent p-value/rejection/status")
  }
  invisible(TRUE)
}

# Complete one scenario, resuming validated repetitions and saving checkpoints.
slads_worker <- function(i, grid, config, signature, directory, resume, checkpoint_every) {
  g <- grid[i, , drop = FALSE]
  path <- file.path(directory, "checkpoints", paste0(g$scenario, ".rds"))
  cp <- list(signature = signature, scenario = g$scenario, repetitions = list())
  if (resume && file.exists(path)) {
    cp <- tryCatch(readRDS(path), error = function(e) stop("Cannot read checkpoint ", path,
      ": ", conditionMessage(e)))
    slads_validate_checkpoint(cp, g, config, signature, path)
  }
  done <- cp$repetitions
  root <- slads_sqrt(slads_covariance(g$p, g$covariance))
  if (length(done) < config$R)
    for (r in seq.int(length(done) + 1L, config$R)) {
      slads_set_seed(slads_seed(config$seed, g$scenario, r, "data"))
      X <- slads_draw(g$m, g$law_x, root)
      Y <- slads_draw(g$m, g$law_y, root)
      done[[r]] <- slads_compare_one(X, Y, g, r, config$seed, config$B, config$alpha)
      if (r%%checkpoint_every == 0L || r == config$R) {
        cp$repetitions <- done
        slads_write_rds(cp, path)
      }
    }
  raw <- do.call(rbind, done)
  rownames(raw) <- NULL
  summary <- slads_summary(raw, config)
  slads_write_csv(summary, file.path(directory, "scenario_summaries", paste0(g$scenario, ".csv")))
  message(
    g$scenario, ": ", config$R, " repetitions complete; failed evaluations = ",
    sum(summary$failed)
  )
  list(summary = summary, checkpoint = path)
}

# Final replicate output and public runner ----------------------------------

# Export ordered checkpoints to one CSV, reading one scenario at a time.
slads_stream_replicates <- function(paths, file) {
  tmp <- tempfile(pattern = ".replicates-", tmpdir = dirname(file))
  on.exit(unlink(tmp))
  first <- TRUE
  for (path in paths) {
    d <- do.call(rbind, readRDS(path)$repetitions)
    write.table(d, tmp, sep = ",", row.names = FALSE, col.names = first, append = !first,
      quote = TRUE, na = "", fileEncoding = "UTF-8")
    first <- FALSE
  }
  if (!file.rename(tmp, file))
    stop("Cannot commit replicate-level CSV: ", file)
}

# Run the comparison with deterministic streams, source snapshots, and checkpoints.
# R/B count repetitions/resamples; worker count does not change scientific settings.
slads_run <- function(grid, R = 1000L, B = 1000L, alpha = 0.05, seed = 20260925L, workers = 7L,
  output_root = "../simulation_results/comparison", resume = TRUE, checkpoint_every = 10L,
  input_path = .slads_input_path) {
  R <- slads_integer(R, "R")
  B <- slads_integer(B, "B", minimum = 19L)
  seed <- slads_integer(seed, "seed", minimum = 0L, maximum = 2147483646L)
  workers <- slads_integer(workers, "workers", maximum = 7L)
  checkpoint_every <- slads_integer(checkpoint_every, "checkpoint_every")
  if (length(alpha) != 1L || !is.numeric(alpha) || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1)
    stop("alpha must be strictly between zero and one.")
  if (length(resume) != 1L || !is.logical(resume) || is.na(resume))
    stop("resume must be TRUE or FALSE.")
  if (length(output_root) != 1L || !is.character(output_root) || !nzchar(output_root))
    stop("output_root must be a path.")
  if (1/(B + 1) > alpha)
    warning("The plus-one permutation test cannot reject at this alpha with this B.")
  # New R workers inherit these settings before loading their BLAS runtimes.
  Sys.setenv(
    OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1", BLIS_NUM_THREADS = "1", RCPP_PARALLEL_NUM_THREADS = "1"
  )
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
  }
  config <- list(
    engine = "selected-six-methods-v4", grid = grid, R = R, B = B,
    alpha = as.numeric(alpha), seed = seed, methods = slads_default_methods,
    calibrations = slads_calibrations, code_hash = slads_code_hash(input_path),
    R_version = as.character(getRversion()), platform = R.version$platform,
    rng = "Mersenne-Twister/Inversion/Rejection"
  )
  signature <- slads_hash(config)
  run_id <- paste0("h0_", substr(signature, 1, 16))
  if (!resume)
    run_id <- paste0(run_id, "_fresh_", basename(tempfile(pattern = "run_")))
  directory <- file.path(output_root, run_id)
  for (path in c(
    directory, file.path(directory, "checkpoints"),
    file.path(directory, "scenario_summaries")
  ))
    if (!dir.exists(path) &&
        !dir.create(path, recursive = TRUE, showWarnings = FALSE))
      stop("Cannot create ", path)
  directory <- normalizePath(directory, mustWork = TRUE)
  cfgpath <- file.path(directory, "configuration.rds")
  if (file.exists(cfgpath)) {
    if (!identical(readRDS(cfgpath), config))
      stop("Existing output configuration does not match: ", cfgpath)
  } else slads_write_rds(config, cfgpath)
  # Snapshot the notebook and every source dependency used for this workflow.
  slads_snapshot_sources(directory, input_path)
  slads_write_csv(grid, file.path(directory, "scenario_grid.csv"))
  writeLines(capture.output(str(config)), file.path(directory, "configuration.txt"))
  actual_workers <- min(workers, nrow(grid))
  message(
    "Output: ", directory, "\nScenarios: ", nrow(grid),
    "; methods: ", length(slads_default_methods), "; R=", R, "; B=", B,
    "; workers=", actual_workers, " (one BLAS thread per worker)"
  )
  common <- list(grid = grid, config = config, signature = signature, directory = directory,
    resume = resume, checkpoint_every = checkpoint_every)
  if (actual_workers > 1L) {
    # PSOCK also works in RStudio and on Windows; no forked GUI process.
    cl <- parallel::makePSOCKcluster(actual_workers, outfile = file.path(directory, "workers.log"))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    env <- environment(slads_run)
    candidates <- ls(envir = env, all.names = TRUE)
    # Export private helpers too; method seed aliases are required on every worker.
    exported <- candidates[grepl("^(slads_|\\.slads_|kyfan_|method_|\\.qt_|\\.wang_|pair_levels$)",
      candidates)]
    parallel::clusterExport(cl, exported, envir = env)
    answer <- do.call(
      parallel::parLapplyLB,
      c(list(cl = cl, X = seq_len(nrow(grid)), fun = slads_worker), common)
    )
  } else answer <- do.call(lapply, c(list(X = seq_len(nrow(grid)), FUN = slads_worker), common))
  summary <- do.call(rbind, lapply(answer, `[[`, "summary"))
  rownames(summary) <- NULL
  summary$run_id <- run_id
  expected <- nrow(grid) * length(slads_default_methods)
  if (nrow(summary) != expected || anyDuplicated(paste(summary$scenario, summary$method)))
    stop("Unexpected number of scenario/method results.")
  slads_write_csv(summary, file.path(directory, "summary.csv"))
  slads_stream_replicates(vapply(answer, `[[`, character(1), "checkpoint"), file.path(directory,
    "replicates.csv"))
  slads_write_csv(
    summary[summary$failed > 0L, , drop = FALSE],
    file.path(directory, "failed_cells.csv")
  )
  writeLines(capture.output(sessionInfo()), file.path(directory, "sessionInfo.txt"))
  writeLines(
    c(
      "Completed all requested repetitions. Inspect status/failed counts before using results.",
      paste("run_id:", run_id),
      paste("failed_method_evaluations:", sum(summary$failed))
    ),
    file.path(directory, "COMPLETED.txt")
  )
  if (any(summary$failed > 0L))
    warning("Some methods failed: main rates/intervals are NA for these cells. See failed_cells.csv and replicates.csv.")
  list(summary = summary, config = config, directory = directory, workers = actual_workers)
}

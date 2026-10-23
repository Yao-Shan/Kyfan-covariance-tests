# Archived table values and Ky-Fan reruns; load via slads_load_workflow("tables", input_path).

pair_levels <- c("N-N", "T-T", "E-E", "G-G", "L-L", "N-T", "N-E", "N-G", "N-L", "T-E", "T-G",
  "T-L", "E-G", "E-L", "G-L")

# Saved values and table formatting

# Match records by design and Ky-Fan order, independently of row order.
archive_key <- function(x) paste(x$p, x$m, x$covariance, x$pair, x$k, sep = ":")

# Read 1,020 archived cells and check their keys and rounded MCSEs.
# Archive mode uses these saved values; reruns generate their own results.
read_archived_results <- function() {
  archive_path <- file.path(.slads_source_root, "data", "table_archive.csv")
  if (!file.exists(archive_path)) {
    stop("Archived table values are unavailable: ", archive_path,
      ". Use mode = 'rerun' to generate new results, or restore data/table_archive.csv ",
      "for mode = 'archived'.")
  }
  x <- read.csv(archive_path, stringsAsFactors = FALSE)
  x$R <- 1000L
  x$B <- 1000L
  x$alpha <- 0.05
  # Distinguish cells from the completed rerun from historical TeX values.
  x$provenance <- ifelse(
    x$source == "rerun_main_20261006",
    "embedded_completed_rerun_manuscript", "embedded_historical_tex_values"
  )
  # Archived MCSEs were rounded to three decimals; allow half a rounding unit.
  stopifnot(
    nrow(x) == 1020L,
    !anyDuplicated(archive_key(x)),
    all(abs(x$mcse - sqrt(x$rejection_rate * (1 - x$rejection_rate) / 1000)) <= 0.00050001)
  )
  # Recover counts from saved rates; the archive does not contain replicate-level decisions.
  x$rejections <- as.integer(round(x$rejection_rate * x$R))
  x$exchangeable <- substr(x$pair, 1, 1) == substr(x$pair, 3, 3)
  x[order(x$p, x$m, x$k, match(x$pair, pair_levels), x$covariance), ]
}

# Select manuscript sample sizes and k = 1, 3; covariance rows pass through.
select_manuscript <- function(x) {
  x[
    ((x$p == 10 & x$m %in% c(50, 150)) | (x$p == 40 & x$m == 150)) & x$k %in% c(1, 3),
  ]
}

# Keep C1--C4 first and append any supplementary covariance designs.
format_paper_table <- function(x, p, m, k) {
  z <- x[x$p == p & x$m == m & x$k == k, ]
  out <- data.frame(Pair = pair_levels, check.names = FALSE)
  covariances <- c(paste0("C", 1:4), setdiff(unique(z$covariance), paste0("C", 1:4)))
  for (cv in covariances) {
    zz <- z[z$covariance == cv, ]
    ii <- match(pair_levels, zz$pair)
    # Mark missing scenarios explicitly so partial reruns cannot inherit archived values.
    out[[cv]] <- ifelse(
      is.na(ii), "—", sprintf("%.3f (%.3f)", zz$rejection_rate[ii], zz$mcse[ii])
    )
  }
  out
}

# Write a separate TeX table file from the supplied results.
write_paper_tex <- function(x, file) {
  grids <- unique(x[c("p", "m", "k")])
  out <- character()
  for (i in seq_len(nrow(grids))) {
    g <- grids[i, ]
    t <- format_paper_table(x, g$p, g$m, g$k)
    t$Pair <- gsub("-", "--", t$Pair, fixed = TRUE)
    q <- x[x$p == g$p & x$m == g$m & x$k == g$k, , drop = FALSE]
    counts <- if (all(c("R", "B") %in% names(q))) {
      sprintf(" R=%d, B=%d.", q$R[1], q$B[1])
    } else {
      ""
    }
    out <- c(
      out, "\\begin{table}[htbp]", "\\centering",
      sprintf(
        "\\caption{Null rejection proportions (Monte Carlo standard errors): p=%d, m=%d, k=%d.%s}",
        g$p, g$m, g$k, counts
      ),
      paste0("\\begin{tabular}{l", paste(rep("c", ncol(t) - 1L), collapse = ""), "}"),
      "\\toprule", paste0(paste(names(t), collapse = " & "), " \\\\"),
      "\\midrule", apply(t, 1, function(r) paste0(paste(r, collapse = " & "), " \\\\")),
      "\\bottomrule", "\\end{tabular}", "\\end{table}", ""
    )
  }
  writeLines(out, file, useBytes = TRUE)
}

# Scenario design, simulation and output

# Build manuscript, archive, or diagnostic grids independently of saved results.
slads_grid <- function(profile = c("manuscript", "archive", "smoke", "spectral")) {
  profile <- match.arg(profile)
  if (profile == "smoke") {
    g <- data.frame(p = 10L, m = 50L, covariance = c("C1", "C1", "C3", "C4"), pair = c("N-N",
      "N-T", "E-E", "N-L"), stringsAsFactors = FALSE)
    g$ks <- "1,3"
  } else if (profile == "spectral") {
    # Supplementary dimensions and spectra still use the same covariance for both samples.
    g <- data.frame(p = c(100L, 100L, 500L, 500L), m = 100L, covariance = c("SPIKE10",
      "SPIKE10", "C1", "C1"), pair = c("N-N", "L-L", "N-N", "L-L"), ks = "1,3", stringsAsFactors = FALSE)
  } else {
    # Preserve the historical sample sizes, design order, and requested Ky-Fan orders.
    sizes <- data.frame(
      p = c(10L, 10L, 20L, 20L, 40L, 40L),
      m = c(50L, 150L, 50L, 150L, 50L, 150L),
      ks = c("1,3", "1,3", "1,2,3", "1,2,5", "1,2,3", "1,2,3,5"),
      stringsAsFactors = FALSE
    )
    if (profile == "manuscript") {
      sizes <- sizes[c(1L, 2L, 6L), , drop = FALSE]
      sizes$ks <- "1,3"
    }
    # Store orders together so a scenario shares data and permutations across all k values.
    g <- do.call(rbind, lapply(seq_len(nrow(sizes)), function(i) {
      combinations <- expand.grid(covariance = paste0("C", 1:4), pair = pair_levels,
        stringsAsFactors = FALSE)
      data.frame(p = sizes$p[i], m = sizes$m[i], combinations, ks = sizes$ks[i],
        stringsAsFactors = FALSE)
    }))
  }
  g$law_x <- substr(g$pair, 1, 1)
  g$law_y <- substr(g$pair, 3, 3)
  g$exchangeable <- g$law_x == g$law_y
  g$scenario <- with(g, paste0("p", p, "_m", m, "_", covariance, "_", pair))
  rownames(g) <- NULL
  g
}

slads_hash <- function(x) {
  f <- tempfile()
  on.exit(unlink(f))
  saveRDS(x, f, version = 2)
  unname(tools::md5sum(f))
}

slads_code_hash <- function() {
  slads_source_hash(.slads_input_path)
}

slads_write_rds <- function(x, path) {
  # Rename only after serialization finishes so an interrupted write preserves the old checkpoint.
  temp <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(x, temp)
  if (!file.rename(temp, path))
    stop("Cannot commit checkpoint: ", path)
}

slads_wilson <- function(s, n, level = 0.95) {
  # Use valid repetitions as the binomial denominator; no valid trials gives missing bounds.
  if (n < 1L)
    return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - level) / 2)
  ph <- s / n
  d <- 1 + z * z / n
  mid <- (ph + z * z / (2 * n)) / d
  half <- z * sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / d
  c(max(0, mid - half), min(1, mid + half))
}

# Rates, MCSEs, and intervals condition on valid evaluations in the supplied rows.
# Unlike the comparison runner, this summary does not require a complete cell.
slads_summary <- function(raw, alpha = 0.05) {
  groups <- split(raw, interaction(raw$scenario, raw$method, drop = TRUE, lex.order = TRUE))
  result <- lapply(groups, function(d) {
    # Retain excluded evaluations in the status counts, not in the rate denominator.
    ok <- d$status == "ok" & !is.na(d$reject)
    n <- sum(ok)
    s <- sum(d$reject[ok])
    ci <- slads_wilson(s, n)
    rate <- if (n) s / n else NA_real_
    data.frame(
      d[1, c("scenario", "p", "m", "covariance", "pair", "exchangeable", "method")],
      requested = nrow(d), valid = n, not_applicable = sum(d$status == "not_applicable"),
      failed = sum(!ok & d$status != "not_applicable"), rejections = s, rejection_rate = rate,
      mcse = if (n) sqrt(rate * (1 - rate) / n) else NA_real_,
      ci_low = ci[1], ci_high = ci[2], alpha = alpha, mean_seconds = mean(d$seconds),
      row.names = NULL
    )
  })
  do.call(rbind, result)
}

# Run Ky-Fan tables with R sample pairs per scenario and B permutations per pair.
slads_run <- function(grid, R = 1000L, B = 1000L, alpha = 0.05, seed = 20260925L,
                      workers = 1L, output_root = "outputs", resume = TRUE,
                      checkpoint_every = 10L, protected_paths = character()) {
  stopifnot(R >= 1L, R == as.integer(R), length(B) == 1L, is.finite(B),
    B >= 19L, B == as.integer(B), alpha > 0, alpha < 1)
  # Avoid nested BLAS/OpenMP parallelism inside the scenario workers.
  Sys.setenv(
    OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1", BLIS_NUM_THREADS = "1", RCPP_PARALLEL_NUM_THREADS = "1"
  )
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
  }
  # Sources and scientific settings determine the run directory; worker count is excluded.
  config <- list(
    grid = grid, R = R, B = B, alpha = alpha, seed = seed, code_hash = slads_code_hash(),
    R_version = as.character(getRversion()),
    rng = "Mersenne-Twister/Inversion/Rejection", engine = "1.1-modular"
  )
  run_id <- paste0("reproduce_", substr(slads_hash(config), 1, 12))
  outdir <- file.path(output_root, run_id)
  # Keep the fixed comparison reference outside this run directory to prevent overwriting it.
  if (length(protected_paths)) {
    protected <- vapply(protected_paths, normalizePath, character(1),
      winslash = "/", mustWork = TRUE)
    destination <- normalizePath(outdir, winslash = "/", mustWork = FALSE)
    if (any(protected == destination | startsWith(protected, paste0(destination, "/")))) {
      stop("The fixed reference is inside this run's output directory. ",
        "Choose a fixed copy outside this run directory before starting the simulation.")
    }
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  slads_snapshot_sources(outdir, .slads_input_path)
  saveRDS(config, file.path(outdir, "configuration.rds"))
  writeLines(capture.output(str(config)), file.path(outdir, "configuration.txt"))
  # Resume each scenario from its saved rows in the configuration-specific run directory.
  one_scenario <- function(i) {
    g <- grid[i, , drop = FALSE]
    file <- file.path(outdir, paste0(g$scenario, ".rds"))
    done <- if (resume && file.exists(file)) readRDS(file) else list()
    if (length(done) > R)
      stop("Checkpoint larger than requested run.")
    root <- slads_sqrt(slads_covariance(g$p, g$covariance))
    ks <- unique(as.numeric(strsplit(g$ks, ",", fixed = TRUE)[[1]]))
    if (length(done) < R) {
      for (r in seq.int(length(done) + 1L, R)) {
        # Scenario/repetition seeds reproduce both streams after resuming or reordering work.
        slads_set_seed(slads_seed(seed, g$scenario, r, "data"))
        X <- slads_draw(g$m, g$law_x, root)
        Y <- slads_draw(g$m, g$law_y, root)
        slads_set_seed(slads_seed(seed, g$scenario, r, "permutation"))
        tm <- system.time(ky <- kyfan_permutation(X, Y, ks, B))["elapsed"]
        keys <- paste0("k", ks)
        rows <- data.frame(
          method = paste0("KyFan_k", ks), p_value = unname(ky$p_value[keys]),
          statistic = unname(ky$statistic[keys]), reject = unname(ky$p_value[keys]) <= alpha,
          status = "ok", detail = "", seconds = unname(tm), stringsAsFactors = FALSE
        )
        for (nm in c("scenario", "p", "m", "covariance", "pair", "exchangeable")) {
          rows[[nm]] <- g[[nm]]
        }
        rows$replicate <- r
        rows$data_seed <- slads_seed(seed, g$scenario, r, "data")
        done[[r]] <- rows
        if (r %% checkpoint_every == 0L || r == R) {
          slads_write_rds(done, file)
        }
      }
    }
    message(g$scenario, ": ", R, " repetitions complete")
    do.call(rbind, done)
  }
  workers <- as.integer(workers)
  if (is.na(workers) || workers < 1L)
    workers <- 1L
  workers <- min(workers, 7L)
  if (.Platform$OS.type == "windows" && workers > 1L) {
    warning("Windows uses serial execution; set workers=1.")
    workers <- 1L
  }
  # Explicit replicate seeds provide reproducibility without fork-assigned RNG streams.
  all <- if (workers > 1L) {
    parallel::mclapply(
      seq_len(nrow(grid)), one_scenario,
      mc.cores = min(workers, nrow(grid)), mc.set.seed = FALSE
    )
  } else {
    lapply(seq_len(nrow(grid)), one_scenario)
  }
  if (any(vapply(all, inherits, logical(1), "try-error")))
    stop("A worker failed; see checkpoints.")
  raw <- do.call(rbind, all)
  rownames(raw) <- NULL
  summary <- slads_summary(raw, alpha)
  write.csv(raw, file.path(outdir, "replicates.csv"), row.names = FALSE, na = "")
  write.csv(summary, file.path(outdir, "summary.csv"), row.names = FALSE, na = "")
  writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
  list(raw = raw, summary = summary, config = config, directory = outdir)
}

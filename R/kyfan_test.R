# Supplied-data Ky-Fan interface; requires R/simulation_core.R.

# Require equal sample dimensions and matching variable names in the same order.
kyfan_check_samples <- function(X, Y) {
  as_sample <- function(x, name) {
    if (is.data.frame(x)) {
      if (!all(vapply(x, is.numeric, logical(1)))) {
        stop(name, " must contain only numeric measurement columns.", call. = FALSE)
      }
      x <- as.matrix(x)
    }
    if (!is.matrix(x) || !is.numeric(x) || is.complex(x)) {
      stop(name, " must be a numeric matrix or data frame.", call. = FALSE)
    }
    if (nrow(x) < 2L || ncol(x) < 1L || any(!is.finite(x))) {
      stop(name, " needs at least two observations, one variable, and finite values.",
           call. = FALSE)
    }
    x
  }
  X <- as_sample(X, "X")
  Y <- as_sample(Y, "Y")
  if (!identical(dim(X), dim(Y))) {
    stop("X and Y must have equal sample sizes and the same number of variables.",
         call. = FALSE)
  }
  nx <- colnames(X)
  ny <- colnames(Y)
  if (!identical(nx, ny) ||
      (!is.null(nx) && (anyNA(nx) || any(!nzchar(nx)) || anyDuplicated(nx)))) {
    stop("Use identical, unique column names in the same order, or leave both matrices unnamed.",
         call. = FALSE)
  }
  list(X = X, Y = Y)
}

# Read either two supplied samples or two CSV files, preserving column names.
kyfan_read_samples <- function(X = NULL, Y = NULL, x_csv = "", y_csv = "",
                               root = getwd(), csv_header = TRUE, csv_sep = ",") {
  paths <- list(x_csv, y_csv)
  if (!all(vapply(paths, function(x) is.character(x) && length(x) == 1L &&
                  !is.na(x), logical(1)))) {
    stop("x_csv and y_csv must be single file paths.", call. = FALSE)
  }
  has_files <- nzchar(c(x_csv, y_csv))
  has_objects <- c(!is.null(X), !is.null(Y))
  if (any(has_files) && any(has_objects)) {
    stop("Supply either X and Y or two CSV paths, not both.", call. = FALSE)
  }
  if (any(has_objects)) {
    if (!all(has_objects)) stop("Supply both X and Y.", call. = FALSE)
    return(kyfan_check_samples(X, Y))
  }
  if (!all(has_files)) stop("Supply both x_csv and y_csv.", call. = FALSE)
  if (!is.logical(csv_header) || length(csv_header) != 1L || is.na(csv_header) ||
      !is.character(csv_sep) || length(csv_sep) != 1L || is.na(csv_sep) ||
      nchar(csv_sep) != 1L) {
    stop("csv_header must be TRUE/FALSE; csv_sep must be one character.", call. = FALSE)
  }
  read_sample <- function(path) {
    path <- path.expand(path)
    if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) path <- file.path(root, path)
    path <- normalizePath(path, mustWork = TRUE)
    read.table(path, header = csv_header, sep = csv_sep, row.names = NULL,
               check.names = FALSE, stringsAsFactors = FALSE, comment.char = "",
               fill = FALSE, blank.lines.skip = FALSE, quote = "\"")
  }
  kyfan_check_samples(read_sample(x_csv), read_sample(y_csv))
}

# Test prespecified orders with shared permutations and a relative tolerance for ties.
# Each order has its own unadjusted p-value.
method_kyfan <- function(X, Y, k = 1L, B = 1000L, seed = 20260925L,
                         alpha = 0.05,
                         relative_tolerance = 100 * .Machine$double.eps) {
  samples <- kyfan_check_samples(X, Y)
  X <- samples$X
  Y <- samples$Y
  whole_numbers <- function(x, lower, upper) {
    is.numeric(x) && !is.complex(x) && length(x) > 0L && all(is.finite(x)) &&
      all(x == trunc(x) & x >= lower & x <= upper)
  }
  if (!whole_numbers(k, 1L, ncol(X))) {
    stop("k must contain prespecified integers between 1 and p.", call. = FALSE)
  }
  if (length(B) != 1L || !whole_numbers(B, 1L, .Machine$integer.max - 1L)) {
    stop("B must be a positive integer.", call. = FALSE)
  }
  if (length(seed) != 1L || !whole_numbers(seed, 0L, .Machine$integer.max)) {
    stop("seed must be a nonnegative integer within R's integer range.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) stop("alpha must lie in (0, 1).", call. = FALSE)
  if (!is.numeric(relative_tolerance) || length(relative_tolerance) != 1L ||
      !is.finite(relative_tolerance) || relative_tolerance < 0 ||
      relative_tolerance >= 1) {
    stop("relative_tolerance must lie in [0, 1).", call. = FALSE)
  }
  k <- sort(unique(as.integer(k)))
  B <- as.integer(B)

  # Restore RNG settings and seed presence even if the calculation fails.
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  slads_set_seed(as.integer(seed))
  fit <- tryCatch(
    kyfan_permutation(X, Y, ks = k, B = B),
    error = function(e) stop(
      "Ky-Fan calculation failed: ", conditionMessage(e),
      ". For very large values, rescale both samples by the same factor.",
      call. = FALSE
    )
  )
  if (any(!is.finite(fit$statistic)) || any(!is.finite(fit$perm))) {
    stop("Non-finite statistics: check the data scale before testing.", call. = FALSE)
  }

  # Recount numerical ties with a relative threshold and the same plus-one correction.
  # No absolute tolerance is added, so common rescaling preserves the threshold rule.
  tolerance <- relative_tolerance * abs(fit$statistic)
  exceedances <- colSums(sweep(fit$perm, 2L, fit$statistic - tolerance, ">="))
  p_value <- setNames((1 + exceedances) / (B + 1), names(fit$statistic))
  summary <- data.frame(
    m = nrow(X), p = ncol(X), k = k,
    kyfan_norm = unname(fit$observed[k]),
    statistic = unname(fit$statistic), B = B,
    exceedances = unname(exceedances), p_value = unname(p_value),
    alpha = alpha, reject = unname(p_value <= alpha),
    tie_tolerance = unname(tolerance)
  )
  colnames(fit$perm) <- names(p_value)
  list(
    statistic = fit$statistic, p_value = p_value, status = "ok",
    summary = summary, perm = fit$perm,
    permutation_statistics = data.frame(permutation = seq_len(B), fit$perm),
    B = B, k = k, seed = seed, alpha = alpha,
    relative_tolerance = relative_tolerance
  )
}

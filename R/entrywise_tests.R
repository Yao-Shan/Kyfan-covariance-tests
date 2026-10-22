# Base-R LC, CLX, and CZZW tests; observations are rows and variables are columns.
# Sources and assumptions: 02_compare_methods_H0.Rmd.

# Input validation and common failure records --------------------------------

# Failed tests return missing p-values and statistics, never rejection decisions.
.slads_entry_failure <- function(status, detail) {
  list(p_value = NA_real_, statistic = NA_real_, status = status, detail = detail)
}

# Validate dimensions and finite observations; return NULL or a failure record.
.slads_entry_validate <- function(X, Y, min_n = 4L, min_p = 1L) {
  if (!is.matrix(X) || !is.matrix(Y) || !is.numeric(X) || !is.numeric(Y))
    return(.slads_entry_failure("invalid_input", "X and Y must be numeric matrices."))
  if (ncol(X) != ncol(Y) || ncol(X) < min_p || nrow(X) < min_n || nrow(Y) < min_n)
    return(.slads_entry_failure("invalid_dimensions", paste0("Equal column counts, p >= ",
      min_p, ", and each n >= ", min_n, " are required.")))
  if (any(!is.finite(X)) || any(!is.finite(Y)))
    return(.slads_entry_failure("invalid_input", "Non-finite observations are not allowed."))
  NULL
}

# Li--Chen: full U-statistic with a standard-normal upper tail ----------------

# Exact O(n^2 p) form of Li--Chen (2012), Eq. (2.1), retaining all mean corrections.
# Translation invariance permits centering. With zero-diagonal H and r = rowSums(H):
# ordered triple sum = sum(r^2) - sum(H^2);
# ordered quadruple sum = sum(H)^2 - 4 sum(r^2) + 2 sum(H^2).
.slads_trace_square_u <- function(X, already_centered = FALSE) {
  n <- nrow(X)
  if (n < 4L)
    stop("The exact trace-square U-statistic requires n >= 4.")
  if (!already_centered)
    X <- sweep(X, 2L, colMeans(X), "-")
  H <- tcrossprod(X)
  diag(H) <- 0
  s2 <- sum(H * H)
  r2 <- sum(rowSums(H)^2)
  s3 <- r2 - s2
  s4 <- sum(H)^2 - 4 * r2 + 2 * s2
  s2 / (n * (n - 1)) -
    2 * s3 / (n * (n - 1) * (n - 2)) +
    s4 / (n * (n - 1) * (n - 2) * (n - 3))
}

# Combine the within-group and cross-group trace estimates.
.slads_li_components <- function(X, Y) {
  n1 <- nrow(X)
  n2 <- nrow(Y)
  Xc <- sweep(X, 2L, colMeans(X), "-")
  Yc <- sweep(Y, 2L, colMeans(Y), "-")
  A1 <- .slads_trace_square_u(Xc, already_centered = TRUE)
  A2 <- .slads_trace_square_u(Yc, already_centered = TRUE)
  # Eq. (2.2), including all mean-correction terms, reduces exactly to
  # tr(cov(X) cov(Y)); use the Gram form to avoid forming p-by-p matrices.
  C12 <- sum(tcrossprod(Xc, Yc)^2)/((n1 - 1) * (n2 - 1))
  list(A1 = A1, A2 = A2, C12 = C12, T = A1 + A2 - 2 * C12, n1 = n1, n2 = n2)
}

# LC full U-statistic with a standard-normal upper tail; invalid null scales fail.
method_li <- function(X, Y) {
  invalid <- .slads_entry_validate(X, Y)
  if (!is.null(invalid))
    return(invalid)
  a <- .slads_li_components(X, Y)
  # Cross weights follow Li--Chen (2012), p. 914, and authors' PEtests code.
  sd0 <- 2 * a$A1/a$n2 + 2 * a$A2/a$n1
  if (!is.finite(sd0) || sd0 <= 0 || !is.finite(a$T))
    return(.slads_entry_failure(
      "invalid_null_scale",
      "Nonpositive or non-finite Li--Chen null standard-error estimate."
    ))
  z <- a$T/sd0
  list(
    p_value = pnorm(z, lower.tail = FALSE), statistic = z, status = "ok",
    detail = list(
      raw_statistic = a$T, null_sd = sd0, A1 = a$A1, A2 = a$A2, C12 = a$C12,
      log_p_value = pnorm(z, lower.tail = FALSE, log.p = TRUE),
      calibration = "upper standard-normal tail; exact U-statistic"
    )
  )
}

# Shared CLX/CZZW covariance-entry standardization ---------------------------

# Upper-triangle covariance products include the diagonal and use 1/n,
# as required by CLX/CZZW, rather than cov()'s 1/(n-1); invalid entry variances fail.
.slads_entry_products <- function(X, Y) {
  p <- ncol(X)
  n1 <- nrow(X)
  n2 <- nrow(Y)
  ij <- which(upper.tri(matrix(0, p, p), diag = TRUE), arr.ind = TRUE)
  Xc <- sweep(X, 2L, colMeans(X), "-")
  Yc <- sweep(Y, 2L, colMeans(Y), "-")
  U <- Xc[, ij[, 1L], drop = FALSE] * Xc[, ij[, 2L], drop = FALSE]
  V <- Yc[, ij[, 1L], drop = FALSE] * Yc[, ij[, 2L], drop = FALSE]
  sx <- colMeans(U)
  sy <- colMeans(V)
  U <- sweep(U, 2L, sx, "-")
  V <- sweep(V, 2L, sy, "-")
  denominator <- colMeans(U * U)/n1 + colMeans(V * V)/n2
  if (any(!is.finite(denominator)) || any(denominator <= 0))
    return(.slads_entry_failure(
      "invalid_entry_variance",
      "At least one tested covariance entry has zero or non-finite estimated variance; no entries were silently removed."
    ))
  z <- (sx - sy)/sqrt(denominator)
  if (any(!is.finite(z)))
    return(.slads_entry_failure(
      "invalid_statistic", "Non-finite standardized covariance difference."
    ))
  list(U = U, V = V, denominator = denominator, z = z, n1 = n1, n2 = n2, p = p)
}

# Cai--Liu--Xia: squared maximum with asymptotic Gumbel calibration ----------

# Stable Gumbel upper tail for squared maximum M in dimension p.
.slads_cai_gumbel <- function(M, p) {
  shifted <- M - 4 * log(p) + log(log(p))
  log_a <- -shifted/2 - log(8 * pi)/2
  # log(1-exp(-a)), evaluated stably even for extremely small tail areas.
  log_p <- if (log_a < -36) log_a else log(-expm1(-exp(log_a)))
  list(p_value = exp(log_p), log_p_value = log_p, shifted = shifted)
}

# CLX squared maximum with asymptotic Gumbel calibration, including the diagonal.
method_cai <- function(X, Y) {
  invalid <- .slads_entry_validate(X, Y, min_n = 2L, min_p = 2L)
  if (!is.null(invalid))
    return(invalid)
  a <- .slads_entry_products(X, Y)
  if (!is.null(a$status))
    return(a)
  M <- max(a$z^2)
  g <- .slads_cai_gumbel(M, a$p)
  list(
    p_value = g$p_value, statistic = M, status = "ok",
    detail = list(
      centered_statistic = g$shifted, log_p_value = g$log_p_value,
      calibration = "Cai--Liu--Xia asymptotic Gumbel tail"
    )
  )
}

# Chang et al.: absolute maximum with Gaussian multipliers ------------------

# CZZW Gaussian multiplier test; bootstrap_batch bounds resampling memory.
# Default tail: count/B; mc_plus_one selects (count+1)/(B+1).
method_chang <- function(X, Y, B = 999L, mc_plus_one = FALSE, bootstrap_batch = 200L) {
  invalid <- .slads_entry_validate(X, Y, min_n = 2L, min_p = 2L)
  if (!is.null(invalid))
    return(invalid)
  if (length(B) != 1L || !is.finite(B) || B < 1L ||
      B > .Machine$integer.max || B != floor(B) ||
      length(bootstrap_batch) != 1L || !is.finite(bootstrap_batch) ||
      bootstrap_batch < 1L || bootstrap_batch > .Machine$integer.max ||
      bootstrap_batch != floor(bootstrap_batch))
    return(.slads_entry_failure(
      "invalid_bootstrap_count",
      "B and bootstrap_batch must be positive finite integers."
    ))
  a <- .slads_entry_products(X, Y)
  if (!is.null(a$status))
    return(a)
  B <- as.integer(B)
  bootstrap_batch <- as.integer(bootstrap_batch)
  observed <- max(abs(a$z))
  # Share each independent observation-level Gaussian multiplier across its entries.
  W <- rbind(a$U/a$n1, -a$V/a$n2)
  W <- sweep(W, 2L, sqrt(a$denominator), "/")
  boot <- numeric(B)
  for (start in seq.int(1L, B, by = bootstrap_batch)) {
    idx <- start:min(B, start + bootstrap_batch - 1L)
    g <- matrix(rnorm((a$n1 + a$n2) * length(idx)), nrow = a$n1 + a$n2)
    values <- crossprod(g, W)
    boot[idx] <- apply(abs(values), 1L, max)
  }
  # Include ties in the empirical conditional upper tail.
  count <- sum(boot >= observed)
  pv <- if (isTRUE(mc_plus_one)) (count + 1)/(B + 1) else count/B
  list(
    p_value = pv, statistic = observed, status = "ok",
    detail = list(
      B = B, mc_plus_one = isTRUE(mc_plus_one), bootstrap_statistics = boot,
      calibration = "Chang et al. Gaussian multiplier bootstrap"
    )
  )
}

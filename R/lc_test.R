# LC (Li--Chen, 2012): full mean-corrected U-statistic with normal calibration.
# Requires R/test_helpers.R; observations are rows and variables are columns.

# Estimate tr(Sigma^2) using Li--Chen (2012), Eq. (2.1), in O(n^2 p) time.
# Translation invariance permits centering; zeroing H's diagonal lets row sums
# recover the ordered triple and quadruple terms with all indices distinct.
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

# Combine trace estimates into an unbiased estimate of ||Sigma1 - Sigma2||_F^2.
.slads_li_components <- function(X, Y) {
  n1 <- nrow(X)
  n2 <- nrow(Y)
  Xc <- sweep(X, 2L, colMeans(X), "-")
  Yc <- sweep(Y, 2L, colMeans(Y), "-")
  A1 <- .slads_trace_square_u(Xc, already_centered = TRUE)
  A2 <- .slads_trace_square_u(Yc, already_centered = TRUE)
  # Eq. (2.2) reduces to tr(cov(X) cov(Y)), with denominators n1 - 1 and n2 - 1.
  # The Gram form avoids constructing p-by-p sample covariance matrices.
  C12 <- sum(tcrossprod(Xc, Yc)^2)/((n1 - 1) * (n2 - 1))
  list(A1 = A1, A2 = A2, C12 = C12, T = A1 + A2 - 2 * C12, n1 = n1, n2 = n2)
}

# Calibrate the full U-statistic with a standard-normal upper tail.
# A nonpositive or non-finite null scale returns a failure record.
method_lc <- function(X, Y) {
  invalid <- .slads_entry_validate(X, Y)
  if (!is.null(invalid))
    return(invalid)
  a <- .slads_li_components(X, Y)
  # Cross-sample weights follow Li--Chen (2012), p. 914, and the authors' PEtests code.
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


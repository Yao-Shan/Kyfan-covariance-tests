# WZZ covariance test with three-cumulant chi-square calibration.
# Source: 2022 author manuscript arXiv:2212.12338 (journal article: 2024).
# The journal version is unaudited; see reports/h0_comparison.Rmd.
# Base R; observations are rows and variables are columns.

# Convert samples to matrices; invalid dimensions or values raise an error.
.qt_validate <- function(X, Y, min_n = 4L, min_p = 1L) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  if (!is.numeric(X) || !is.numeric(Y) || ncol(X) != ncol(Y) ||
      ncol(X) < min_p || nrow(X) < min_n || nrow(Y) < min_n ||
      any(!is.finite(X)) || any(!is.finite(Y))) {
    stop("X and Y must be finite numeric matrices with equal columns and adequate sample sizes.")
  }
  list(X = X, Y = Y)
}

# Return a missing p-value with the failure status and available diagnostics.
.qt_fail <- function(status, detail, statistic = NA_real_, ...) {
  c(list(p_value = NA_real_, statistic = statistic, status = status, detail = detail), list(...))
}

# Center both feature sets in a within- or cross-group Gram matrix.
.qt_center_gram <- function(K) {
  sweep(sweep(K, 1L, rowMeans(K), "-"), 2L, colMeans(K), "-") + mean(K)
}

# Compute the statistic and feature-covariance traces after sample-mean centering.
# The identity vec(xx')'vec(yy') = (x'y)^2 avoids forming p^2-by-p^2 matrices.
.wang_trace_components <- function(X, Y) {
  n1 <- nrow(X)
  n2 <- nrow(Y)
  X <- sweep(X, 2L, colMeans(X), "-")
  Y <- sweep(Y, 2L, colMeans(Y), "-")
  K11 <- tcrossprod(X)^2
  K22 <- tcrossprod(Y)^2
  K12 <- tcrossprod(X, Y)^2
  # Feature covariances use n - 1; cross-group traces need both groups' factors.
  A <- .qt_center_gram(K11)/(n1 - 1)
  D <- .qt_center_gram(K22)/(n2 - 1)
  B <- .qt_center_gram(K12)/sqrt((n1 - 1) * (n2 - 1))
  a1 <- sum(diag(A))
  d1 <- sum(diag(D))
  list(
    n1 = n1, n2 = n2,
    statistic = mean(K11) + mean(K22) - 2 * mean(K12) - a1/n1 - d1/n2,
    a1 = a1, d1 = d1, a2 = sum(A * A), d2 = sum(D * D), ad = sum(B * B),
    a3 = sum((A %*% A) * A), d3 = sum((D %*% D) * D),
    a2d = sum((A %*% B) * B), ad2 = sum((B %*% D) * B)
  )
}

# Match a shifted, scaled chi-square to zero mean and the estimated cumulants.
# Nonpositive cumulants or invalid matching parameters return a failure record.
method_wzz <- function(X, Y) {
  xy <- .qt_validate(X, Y, min_n = 4L)
  z <- .wang_trace_components(xy$X, xy$Y)
  n1 <- z$n1
  n2 <- z$n2
  # Bias-correct the second- and third-order feature-covariance trace estimates.
  est2 <- function(t1, t2, n) {
    (n - 1)^2/((n - 2) * (n + 1)) * (t2 - t1^2/(n - 1))
  }
  est3 <- function(t1, t2, t3, n) {
    (n - 1)^4 / ((n^2 + n - 6) * (n^2 - 2 * n - 3)) *
      (t3 - 3 * t1 * t2/(n - 1) + 2 * t1^3/(n - 1)^2)
  }
  aa <- est2(z$a1, z$a2, n1)
  dd <- est2(z$d1, z$d2, n2)
  aaa <- est3(z$a1, z$a2, z$a3, n1)
  ddd <- est3(z$d1, z$d2, z$d3, n2)
  aad <- (n1 - 1)/((n1 - 2) * (n1 + 1)) * ((n1 - 1) * z$a2d - z$ad * z$a1)
  add <- (n2 - 1)/((n2 - 2) * (n2 + 1)) * ((n2 - 1) * z$ad2 - z$ad * z$d1)
  K2 <- 2 * (aa/(n1 * (n1 - 1)) + 2 * z$ad/(n1 * n2) + dd/(n2 * (n2 - 1)))
  K3 <- 8 * ((n1 - 2) * aaa/(n1^2 * (n1 - 1)^2) +
    3 * aad/(n1^2 * n2) + 3 * add/(n1 * n2^2) +
    (n2 - 2) * ddd/(n2^2 * (n2 - 1)^2))
  if (any(!is.finite(c(K2, K3))) || K2 <= 0 || K3 <= 0) {
    return(.qt_fail(
      "invalid_cumulants",
      paste0(
        "The published 3-cumulant chi-square approximation requires K2 > 0 and K3 > 0; ",
        "the unmodified estimates do not satisfy this. No truncation or fallback was used."
      ),
      statistic = z$statistic, cumulant2 = K2, cumulant3 = K3
    ))
  }
  # Match the estimated cumulants without truncation or a fallback calibration.
  beta0 <- -2 * K2^2/K3
  beta1 <- K3/(4 * K2)
  df <- 8 * K2^3/K3^2
  if (any(!is.finite(c(beta0, beta1, df))) || beta1 <= 0 || df <= 0) {
    return(.qt_fail("invalid_approximation", "Non-finite 3-cumulant matching parameters.",
      statistic = z$statistic, cumulant2 = K2, cumulant3 = K3))
  }
  pval <- pchisq((z$statistic - beta0)/beta1, df = df, lower.tail = FALSE)
  list(
    p_value = pval, statistic = z$statistic, status = "ok",
    detail = "Published three-cumulant matched chi-square normal-reference calibration; sample-mean centering.",
    cumulant2 = K2, cumulant3 = K3, beta0 = beta0, beta1 = beta1, df = df
  )
}

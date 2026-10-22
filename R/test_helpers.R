# Shared input checks and CLX/CZZW covariance-entry products.
# Observations are rows and variables are columns.

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


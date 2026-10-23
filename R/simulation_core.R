# Base-R simulation generators and fixed-order Ky-Fan permutation tests.

# Derive seeds from the master seed and ordered scenario/repetition/purpose tags.
# Worker IDs are excluded; changing the key encoding changes the random streams.
slads_seed <- function(master, ...) {
  key <- paste(..., collapse = "|")
  h <- as.double(master %% 2147483646)
  for (v in utf8ToInt(key)) h <- (131 * h + v) %% 2147483646
  as.integer(h + 1)
}

# Fix all RNG algorithms, including sample.int(), before setting each stream.
slads_set_seed <- function(seed) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
}

# C1--C4 are the paper covariance designs; SPIKE10 supports
# slads_grid(profile = "spectral") in R/table_reproduction.R.
slads_covariance <- function(p, model = "C1") {
  stopifnot(p >= 5L, p == as.integer(p))
  if (model == "C1")
    return(diag(p))
  if (model == "C2")
    return(0.5^abs(outer(seq_len(p), seq_len(p), "-")))
  if (model == "C3")
    return(diag(c(rep(4, 5), rep(1, p - 5))))
  if (model == "C4") {
    stopifnot(p %% 2L == 0L)
    S <- diag(p)
    start <- 1L
    # One size-four block followed by size-two blocks, each with correlation 0.6.
    for (d in c(4L, rep(2L, (p - 4L) / 2L))) {
      ii <- seq.int(start, length.out = d)
      S[ii, ii] <- 0.4 * diag(d) + 0.6
      start <- start + d
    }
    return(S)
  }
  if (model == "SPIKE10")
    return(diag(c(10, rep(1, p - 1))))
  stop("Unknown covariance model: ", model)
}

# Use the symmetric positive-definite square root to preserve the generator.
slads_sqrt <- function(Sigma) {
  e <- eigen(Sigma, symmetric = TRUE)
  stopifnot(min(e$values) > 0)
  tcrossprod(sweep(e$vectors, 2L, sqrt(e$values), "*"), e$vectors)
}

# Draw independent rows with the target covariance from standardized laws.
# T/E/G share one random scale per row; draw order determines reproducible streams.
slads_draw <- function(n, law, sqrt_sigma) {
  p <- ncol(sqrt_sigma)
  if (law == "L") {
    Z <- matrix(rlogis(n * p), n, p) * sqrt(3) / pi
  } else {
    Z <- matrix(rnorm(n * p), n, p)
    scale <- switch(
      law,
      N = rep(1, n),
      T = sqrt(6 / rchisq(n, 8)),
      E = sqrt(rexp(n)),
      G = sqrt(rgamma(n, shape = 3, scale = 1 / 3)),
      stop("Unknown standardized law: ", law)
    )
    # Column-wise recycling applies the same scale to every coordinate in a row.
    Z <- Z * scale
  }
  Z %*% sqrt_sigma
}

# For symmetric D, position k is its Ky-Fan norm: the k largest singular values summed.
kyfan_curve <- function(D) cumsum(sort(
  abs(eigen(D, symmetric = TRUE, only.values = TRUE)$values), decreasing = TRUE
))

# Test prespecified orders for equal-size m x p samples using shared permutations.
# Rows of assignments select the first group's pooled indices; their count overrides B.
# Orders are not selected from the data; each p-value is unadjusted.
kyfan_permutation <- function(X, Y, ks = c(1L, 3L), B = 1000L,
                             keep_spectra = FALSE, assignments = NULL) {
  stopifnot(
    is.matrix(X), is.matrix(Y), nrow(X) == nrow(Y), ncol(X) == ncol(Y),
    all(is.finite(X)), all(is.finite(Y))
  )
  m <- nrow(X)
  p <- ncol(X)
  stopifnot(m > 1L, length(ks) > 0L, all(is.finite(ks)),
    all(ks >= 1L & ks <= p), all(ks == trunc(ks)))
  ks <- sort(unique(as.integer(ks)))
  if (is.null(assignments)) {
    stopifnot(length(B) == 1L, is.finite(B), B >= 1L, B == as.integer(B))
  } else {
    stopifnot(is.matrix(assignments), ncol(assignments) == m, nrow(assignments) >= 1L,
      all(is.finite(assignments)), all(assignments == trunc(assignments)),
      all(assignments >= 1L & assignments <= 2L * m))
    B <- nrow(assignments)
  }
  Z <- rbind(X, Y)
  # Pooled centering leaves covariance differences unchanged and improves stability.
  Z <- sweep(Z, 2L, colMeans(Z), "-")
  total_sum <- colSums(Z)
  total_cross <- crossprod(Z)
  spectrum <- function(ii) {
    A <- Z[ii, , drop = FALSE]
    a <- colSums(A)
    b <- total_sum - a
    AA <- crossprod(A)
    # Recover the covariance difference from pooled sums, using m - 1 in both groups.
    D <- (2 * AA - total_cross - (tcrossprod(a) - tcrossprod(b)) / m) / (m - 1)
    kyfan_curve(D)
  }
  observed <- spectrum(seq_len(m))
  out <- matrix(NA_real_, B, if (keep_spectra) p else length(ks))
  for (b in seq_len(B)) {
    ii <- if (is.null(assignments)) {
      sample.int(2L * m, m, replace = FALSE)
    } else {
      assignments[b, ]
    }
    if (!is.null(assignments)) {
      stopifnot(!anyDuplicated(ii))
    }
    v <- spectrum(ii)
    out[b, ] <- if (keep_spectra) v else v[ks]
  }
  chosen <- if (keep_spectra) out[, ks, drop = FALSE] else out
  # The common sqrt(m) factor cancels in upper-tail comparisons.
  # Include ties and apply the plus-one correction separately to each order.
  pv <- (1 + colSums(sweep(chosen, 2L, observed[ks], ">="))) / (B + 1)
  names(pv) <- paste0("k", ks)
  # statistic and perm are scaled by sqrt(m); observed and perm_spectra are not.
  list(p_value = pv, statistic = setNames(sqrt(m) * observed[ks], names(pv)), observed = observed,
    perm = sqrt(m) * chosen, perm_spectra = if (keep_spectra) out else NULL, B = B, ks = ks)
}

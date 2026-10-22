# Base-R sampling and fixed-order Ky-Fan permutation tests.

# Seeds use master seed and ordered scenario/repetition/purpose tags, never worker IDs.
# Keep key encoding and RNG kinds fixed to preserve streams.
slads_seed <- function(master, ...) {
  key <- paste(..., collapse = "|")
  h <- as.double(master %% 2147483646)
  for (v in utf8ToInt(key)) h <- (131 * h + v) %% 2147483646
  as.integer(h + 1)
}

# Restore all RNG algorithms before setting a stream, including sample.int().
slads_set_seed <- function(seed) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
}

# C1--C4 are the paper covariance designs; SPIKE10 serves the optional profile in 01.
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
    # One four-coordinate block, followed by two-coordinate blocks.
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

# Use the design's symmetric positive-definite root to preserve the generator.
slads_sqrt <- function(Sigma) {
  e <- eigen(Sigma, symmetric = TRUE)
  stopifnot(min(e$values) > 0)
  tcrossprod(sweep(e$vectors, 2L, sqrt(e$values), "*"), e$vectors)
}

# Draw independent rows with target covariance; T/E/G share one random scale per row.
# Preserve random-draw order to reproduce the original streams.
slads_draw <- function(n, law, sqrt_sigma) {
  # Each row is one observation; all laws have standardized covariance I_p.
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
    # R recycles the length-n vector down each column, sharing one scale per row.
    Z <- Z * scale
  }
  Z %*% sqrt_sigma
}

# Position k contains the sum of the k largest absolute eigenvalues.
kyfan_curve <- function(D) cumsum(sort(
  abs(eigen(D, symmetric = TRUE, only.values = TRUE)$values), decreasing = TRUE
))

# Fixed-order test for equal-size m x p samples; all orders share permutations.
# Supplied B x m assignments identify first-group pooled indices and override B.
# No data-dependent order selection or nested permutation.
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
  # Global centering is algebraically neutral and improves numerical stability.
  Z <- sweep(Z, 2L, colMeans(Z), "-")
  total_sum <- colSums(Z)
  total_cross <- crossprod(Z)
  spectrum <- function(ii) {
    A <- Z[ii, , drop = FALSE]
    a <- colSums(A)
    b <- total_sum - a
    AA <- crossprod(A)
    # Recover both centered covariances from pooled sums; both use denominator m - 1.
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
  # sqrt(m) is positive and common, so unscaled comparisons give identical p-values.
  # Count ties in the right tail and apply the same plus-one correction for every k.
  pv <- (1 + colSums(sweep(chosen, 2L, observed[ks], ">="))) / (B + 1)
  names(pv) <- paste0("k", ks)
  list(p_value = pv, statistic = setNames(sqrt(m) * observed[ks], names(pv)), observed = observed,
    perm = sqrt(m) * chosen, perm_spectra = if (keep_spectra) out else NULL, B = B, ks = ks)
}

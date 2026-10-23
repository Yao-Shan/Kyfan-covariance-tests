# CZZW (Chang--Zhou--Zhou--Wang, 2017): Gaussian multiplier test for covariance equality.
# Requires R/test_helpers.R; observations are rows and variables are columns.

# Calibrate the maximum absolute standardized difference, including diagonal entries.
# The default empirical tail is count/B; mc_plus_one selects (count + 1)/(B + 1).
method_czzw <- function(X, Y, B = 999L, mc_plus_one = FALSE, bootstrap_batch = 200L) {
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
  # Each observation's independent Gaussian multiplier is shared across its entries.
  W <- rbind(a$U/a$n1, -a$V/a$n2)
  W <- sweep(W, 2L, sqrt(a$denominator), "/")
  boot <- numeric(B)
  # Batching bounds multiplier memory while preserving the Gaussian draw order.
  for (start in seq.int(1L, B, by = bootstrap_batch)) {
    idx <- start:min(B, start + bootstrap_batch - 1L)
    g <- matrix(rnorm((a$n1 + a$n2) * length(idx)), nrow = a$n1 + a$n2)
    values <- crossprod(g, W)
    boot[idx] <- apply(abs(values), 1L, max)
  }
  # Include ties in the conditional bootstrap upper tail.
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

# CLX (Cai--Liu--Xia, 2013): squared entrywise maximum with Gumbel calibration.
# Requires R/test_helpers.R; observations are rows and variables are columns.

# Stable Gumbel upper tail for squared maximum M in dimension p.
.slads_cai_gumbel <- function(M, p) {
  shifted <- M - 4 * log(p) + log(log(p))
  log_a <- -shifted/2 - log(8 * pi)/2
  # log(1-exp(-a)), evaluated stably even for extremely small tail areas.
  log_p <- if (log_a < -36) log_a else log(-expm1(-exp(log_a)))
  list(p_value = exp(log_p), log_p_value = log_p, shifted = shifted)
}

# CLX squared maximum with asymptotic Gumbel calibration, including the diagonal.
method_clx <- function(X, Y) {
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


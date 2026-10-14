# Share the Table 1--3 reference with report 01.
.ms_paper_reference <- local({
  ref <- read.csv(
    file.path(.slads_source_root, "data", "manuscript_reference.csv"),
    stringsAsFactors = FALSE
  )
  # Reference runs use R=1000; recover counts and unrounded MCSEs from rates.
  data.frame(
    scenario = with(ref, paste0("p", p, "_m", m, "_", covariance, "_", pair)),
    method = paste0("KyFan_k", ref$k),
    paper_rejections = as.integer(round(1000 * ref$rejection_rate)),
    paper_rate = ref$rejection_rate,
    paper_mcse = sqrt(ref$rejection_rate * (1 - ref$rejection_rate) / 1000),
    stringsAsFactors = FALSE
  )
})

# Match simulated baseline rows to historical counts, rates, and MCSEs.
# Check simulation parameters and recorded R version separately.
ms_compare_paper_reference <- function(result) {
  base <- result$summary[result$summary$mean_id == "baseline", , drop = FALSE]
  key <- paste(base$scenario, base$method, sep = "|")
  ref_key <- paste(.ms_paper_reference$scenario, .ms_paper_reference$method, sep = "|")
  at <- match(key, ref_key)
  if (anyNA(at))
    stop("A baseline scenario is absent from the Table 1--3 reference.")
  ref <- .ms_paper_reference[at, , drop = FALSE]
  comparable <- result$config$R == 1000L && result$config$B == 1000L && result$config$alpha ==
    0.05 && result$config$seed == 20260925L
  data.frame(
    base[c("scenario", "p", "m", "covariance", "pair", "method", "k")],
    new_rejections = base$rejections,
    new_rate = base$type1_error_rate, new_mcse = base$mcse, paper_rejections = ref$paper_rejections,
    paper_rate = ref$paper_rate, paper_mcse = ref$paper_mcse,
    rate_difference = base$type1_error_rate - ref$paper_rate,
    same_simulation_parameters = comparable,
    same_R_version = as.character(getRversion()) == "4.5.1",
    valid_new_cell = base$complete,
    exact_rejection_count_match = if (comparable)
      ifelse(base$complete, base$rejections == ref$paper_rejections, NA) else NA,
    reference_source = "Table1-3 completed rerun 20261006, checked against current v8 on 2026-10-10",
    row.names = NULL, stringsAsFactors = FALSE
  )
}

# Draw Appendix C Wilson intervals from the saved comparison summary.
# Set parameters below; relative paths resolve from the repository root.
params <- list(
  input_csv = "data/comparison_summary.csv",
  figure_dir = "results/figures"
)

# Find this entry file for both Rscript and source() calls.
.slads_entry_file <- local({
  files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)[1L]
  }, character(1))
  files <- files[nzchar(files)]
  if (length(files)) tail(files, 1L) else {
    arg <- grep("^--file=", commandArgs(), value = TRUE)
    if (!length(arg)) stop("Run this script with source() or Rscript.")
    sub("^--file=", "", arg[1L])
  }
})
# Handle source(..., chdir = TRUE), which changes the working directory.
if (!file.exists(.slads_entry_file) && file.exists(basename(.slads_entry_file))) {
  .slads_entry_file <- basename(.slads_entry_file)
}
.slads_entry_file <- normalizePath(.slads_entry_file, winslash = "/", mustWork = TRUE)
source(file.path(dirname(.slads_entry_file), "..", "R", "entrypoint.R"), local = TRUE)
result <- slads_entry("appendix", params, .slads_entry_file)

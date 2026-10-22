# Compare Ky-Fan, LC, CLX, CZZW and WZZ under equal covariance.
# Edit these parameters; relative paths start at the repository root.
params <- list(
  profile = "manuscript",
  R = 1000,
  B = 1000,
  alpha = 0.05,
  seed = 20260925,
  workers = 7,
  resume = TRUE,
  checkpoint_every = 10,
  output_root = "results/comparison"
)

# Locate this script when run with Rscript or source() from any directory.
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
# source(..., chdir = TRUE) evaluates in the script's directory.
if (!file.exists(.slads_entry_file) && file.exists(basename(.slads_entry_file))) {
  .slads_entry_file <- basename(.slads_entry_file)
}
.slads_entry_file <- normalizePath(.slads_entry_file, winslash = "/", mustWork = TRUE)
source(file.path(dirname(.slads_entry_file), "..", "R", "entrypoint.R"), local = TRUE)
result <- slads_entry("comparison", params, .slads_entry_file)

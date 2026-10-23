# Compare two Ky-Fan orders with LC, CLX, CZZW, and WZZ under H0.
# Set parameters below; relative paths resolve from the repository root.
params <- list(
  profile = "manuscript",
  R = 1000, # Independent sample pairs per scenario.
  B = 1000, # Ky-Fan permutations and CZZW Gaussian multiplier draws.
  alpha = 0.05,
  seed = 20260925,
  workers = 7, # Maximum number of parallel scenario workers.
  resume = TRUE, # Reuse checkpoints with matching configuration and source.
  checkpoint_every = 10, # Save after this many repetitions per scenario.
  output_root = "results/comparison"
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
result <- slads_entry("comparison", params, .slads_entry_file)

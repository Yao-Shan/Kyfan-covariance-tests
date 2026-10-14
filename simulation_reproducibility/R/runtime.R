# Source provenance without RNG draws; hash relative paths and contents for portable IDs.

slads_source_manifest <- function(input_path) {
  input_path <- normalizePath(input_path, mustWork = TRUE)
  if (dirname(input_path) != .slads_source_root) {
    stop("Load the workflow for this input_path before starting a run.")
  }
  relative <- c(basename(input_path), .slads_source_files)
  paths <- file.path(.slads_source_root, relative)
  if (any(!file.exists(paths))) stop("A loaded source dependency is missing.")
  data.frame(path = relative, md5 = unname(tools::md5sum(paths)),
    stringsAsFactors = FALSE)
}

# Hash the Rmd, modules, and reference data; text edits create a new run ID.
slads_source_hash <- function(input_path) {
  manifest <- slads_source_manifest(input_path)
  tmp <- tempfile(pattern = "slads-sources-")
  on.exit(unlink(tmp))
  saveRDS(manifest, tmp, version = 2)
  unname(tools::md5sum(tmp))
}

# Snapshot runnable sources before simulation; existing copies must match byte-for-byte.
slads_snapshot_sources <- function(directory, input_path) {
  manifest <- slads_source_manifest(input_path)
  destination <- file.path(directory, "source_snapshot")
  for (i in seq_len(nrow(manifest))) {
    source <- file.path(.slads_source_root, manifest$path[i])
    target <- file.path(destination, manifest$path[i])
    if (!dir.exists(dirname(target)) &&
        !dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)) {
      stop("Cannot create source snapshot directory: ", dirname(target))
    }
    if (!file.exists(target) && !file.copy(source, target)) {
      stop("Cannot preserve source dependency: ", manifest$path[i])
    }
    if (!identical(unname(tools::md5sum(target)), manifest$md5[i])) {
      stop("Source snapshot mismatch: ", manifest$path[i])
    }
  }
  write.csv(manifest, file.path(destination, "source_manifest.csv"), row.names = FALSE)
  invisible(destination)
}

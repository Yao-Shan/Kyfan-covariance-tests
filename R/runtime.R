# Source provenance without RNG draws; hash relative paths and contents for portable IDs.

slads_source_manifest <- function(input_path) {
  input_path <- normalizePath(input_path, winslash = "/", mustWork = TRUE)
  root <- normalizePath(.slads_source_root, winslash = "/", mustWork = TRUE)
  prefix <- paste0(sub("/+$", "", root), "/")
  if (!identical(input_path, .slads_input_path) || !startsWith(input_path, prefix)) {
    stop("Load the workflow for this input_path before starting a run.")
  }
  relative <- c(substring(input_path, nchar(prefix) + 1L), .slads_source_files)
  if (anyNA(relative) || any(!nzchar(relative)) ||
      any(grepl("^(/|[A-Za-z]:|~)|(^|[/\\\\])\\.\\.([/\\\\]|$)", relative))) {
    stop("Source dependencies must use paths within the repository root.")
  }
  paths <- file.path(root, relative)
  if (any(!file.exists(paths) | dir.exists(paths))) {
    stop("A loaded source dependency is missing or is not a file.")
  }
  resolved <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  if (any(!startsWith(resolved, prefix))) {
    stop("A source dependency resolves outside the repository root.")
  }
  data.frame(path = relative, md5 = unname(tools::md5sum(paths)),
    stringsAsFactors = FALSE)
}

# Hash the entry and loaded modules; text edits create a new run ID.
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
  if (!dir.exists(destination) &&
      !dir.create(destination, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create source snapshot directory: ", destination)
  }
  destination <- normalizePath(destination, winslash = "/", mustWork = TRUE)
  prefix <- paste0(sub("/+$", "", destination), "/")
  for (i in seq_len(nrow(manifest))) {
    source <- file.path(.slads_source_root, manifest$path[i])
    target <- file.path(destination, manifest$path[i])
    if (!dir.exists(dirname(target)) &&
        !dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)) {
      stop("Cannot create source snapshot directory: ", dirname(target))
    }
    parent <- normalizePath(dirname(target), winslash = "/", mustWork = TRUE)
    if (!(identical(parent, destination) || startsWith(parent, prefix)) ||
        (file.exists(target) && !startsWith(
          normalizePath(target, winslash = "/", mustWork = TRUE), prefix))) {
      stop("Source snapshot path resolves outside its directory: ", manifest$path[i])
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

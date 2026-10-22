# Small command-line overrides; edit each script's params list for regular runs.
slads_parameters <- function(params, args = character()) {
  overrides <- getOption("slads.params", list())
  if (!is.list(overrides) || (length(overrides) &&
      (is.null(names(overrides)) || any(!nzchar(names(overrides))) ||
       any(!names(overrides) %in% names(params))))) {
    stop("options(slads.params=...) must be a named list of script parameters.")
  }
  params <- modifyList(params, overrides)
  if ("--smoke" %in% args) {
    if (!"profile" %in% names(params)) stop("This script does not run simulations.")
    params$profile <- "smoke"
    params$R <- 3
    params$B <- 39
    params$workers <- 1
    args <- args[args != "--smoke"]
  }
  for (arg in args) {
    if (!grepl("^--[A-Za-z_][A-Za-z_0-9]*=", arg)) {
      stop("Use --smoke or --parameter=value; received: ", arg)
    }
    name <- sub("^--([^=]+)=.*$", "\\1", arg)
    value <- sub("^--[^=]+=", "", arg)
    if (!name %in% names(params)) stop("Unknown parameter: ", name)
    original <- params[[name]]
    if (is.logical(original)) {
      if (!tolower(value) %in% c("true", "false")) {
        stop(name, " must be true or false.")
      }
      value <- tolower(value) == "true"
    } else if (is.numeric(original)) {
      value <- suppressWarnings(as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]]))
      if (!length(value) || anyNA(value) || any(!is.finite(value))) {
        stop(name, " must be numeric.")
      }
    } else if (length(original) > 1L) {
      value <- strsplit(value, ",", fixed = TRUE)[[1L]]
    }
    params[[name]] <- value
  }
  params
}

# Load each workflow in its own environment because runners share function names.
slads_entry <- function(workflow, params, input_path,
                        args = commandArgs(trailingOnly = TRUE)) {
  input_path <- normalizePath(input_path, winslash = "/", mustWork = TRUE)
  root <- dirname(dirname(input_path))
  env <- new.env(parent = globalenv())
  sys.source(file.path(root, "R", "load_workflow.R"), envir = env)
  env$slads_load_workflow(workflow, input_path, envir = env)
  params <- slads_parameters(params, args)
  invisible(env$slads_execute_workflow(workflow, params))
}

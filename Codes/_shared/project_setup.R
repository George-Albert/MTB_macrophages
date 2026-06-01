# Shared project utilities to improve reproducibility across scripts.

locate_project_setup <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = FALSE)

  repeat {
    setup_path <- file.path(current, "Codes", "_shared", "project_setup.R")
    if (file.exists(setup_path)) {
      return(setup_path)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate Codes/_shared/project_setup.R")
    }
    current <- parent
  }
}

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = FALSE)

  repeat {
    has_codes <- dir.exists(file.path(current, "Codes"))
    has_inputs <- dir.exists(file.path(current, "Inputs")) ||
      dir.exists(file.path(current, "Analyses", "Inputs"))
    has_outputs <- dir.exists(file.path(current, "Outputs")) ||
      dir.exists(file.path(current, "Analyses", "Outputs"))

    if (has_codes && has_inputs && has_outputs) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Project root not found. Run scripts from inside the project tree.")
    }
    current <- parent
  }
}

resolve_inputs_dir <- function(project_root) {
  modern <- file.path(project_root, "Inputs")
  legacy <- file.path(project_root, "Analyses", "Inputs")

  if (dir.exists(modern)) return(modern)
  if (dir.exists(legacy)) return(legacy)
  stop("Inputs directory not found. Expected 'Inputs/' or 'Analyses/Inputs/'.")
}

resolve_outputs_dir <- function(project_root) {
  modern <- file.path(project_root, "Outputs")
  legacy <- file.path(project_root, "Analyses", "Outputs")

  if (dir.exists(modern)) return(modern)
  if (dir.exists(legacy)) return(legacy)
  stop("Outputs directory not found. Expected 'Outputs/' or 'Analyses/Outputs/'.")
}

assert_required_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      paste0(
        "Missing R packages: ",
        paste(missing, collapse = ", "),
        ". Install them before running this script."
      )
    )
  }
}

load_required_packages <- function(pkgs) {
  assert_required_packages(pkgs)
  invisible(lapply(pkgs, function(pkg) {
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  }))
}

bootstrap_project <- function(start = getwd(), set_workdir = TRUE) {
  project_root <- find_project_root(start)
  if (isTRUE(set_workdir)) {
    setwd(project_root)
  }

  list(
    project_root = project_root,
    input_dir = resolve_inputs_dir(project_root),
    outputs_root = resolve_outputs_dir(project_root)
  )
}

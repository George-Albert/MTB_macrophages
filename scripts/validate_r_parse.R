bootstrap_root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
repeat {
  setup_candidate <- file.path(bootstrap_root, "Codes", "_shared", "project_setup.R")
  if (file.exists(setup_candidate)) {
    source(setup_candidate)
    break
  }
  parent <- dirname(bootstrap_root)
  if (identical(parent, bootstrap_root)) {
    stop("Could not locate Codes/_shared/project_setup.R")
  }
  bootstrap_root <- parent
}

paths <- bootstrap_project()
tracked_r_files <- system("git ls-files *.R", intern = TRUE)

parse_errors <- list()
for (file in tracked_r_files) {
  ok <- tryCatch(
    {
      parse(file)
      TRUE
    },
    error = function(err) {
      parse_errors[[file]] <<- conditionMessage(err)
      FALSE
    }
  )

  message(if (ok) "OK  " else "BAD ", file)
}

if (length(parse_errors) > 0) {
  print(parse_errors)
  stop("One or more tracked R files failed to parse.")
}

message(length(tracked_r_files), " tracked R files parse OK.")

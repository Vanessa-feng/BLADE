find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)
  markers <- c("code", "data", "demo", "reproduce")

  while (TRUE) {
    if (all(file.exists(file.path(current, markers)))) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not find the BLADE project root. Run scripts from inside the repository.")
    }
    current <- parent
  }
}

project_root <- find_project_root()
setwd(project_root)

project_library <- file.path(project_root, ".Rlib")
if (dir.exists(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}

load_required_packages <- function(packages) {
  missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(
      "Missing R package(s): ",
      paste(missing_packages, collapse = ", "),
      ". Install them before running this script.",
      call. = FALSE
    )
  }

  invisible(lapply(packages, library, character.only = TRUE))
}

# Bootstrap shared configuration and helper functions for scripts in model/.

find_model_root <- function(start_paths) {
  starts <- unique(normalizePath(start_paths[!is.na(start_paths)], mustWork = FALSE))
  for (start in starts) {
    current <- if (dir.exists(start)) start else dirname(start)
    repeat {
      if (file.exists(file.path(current, "src", "config.R")) &&
          dir.exists(file.path(current, "scripts")) &&
          dir.exists(file.path(current, "data"))) {
        return(normalizePath(current, mustWork = TRUE))
      }
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }
  stop("Could not locate the model/ project root.", call. = FALSE)
}

get_script_path <- function() {
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(script_arg) > 0) {
    return(sub("^--file=", "", script_arg[1]))
  }
  NULL
}

check_packages <- function(packages) {
  missing_packages <- setdiff(packages, rownames(installed.packages()))
  if (length(missing_packages) > 0) {
    stop("Install required packages before running this script: ",
         paste(missing_packages, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

write_notes <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, con = path)
}

script_path <- get_script_path()
model_root <- find_model_root(c(getwd(), script_path))
src_root <- file.path(model_root, "src")

source(file.path(src_root, "config.R"))
source(file.path(src_root, "feature_helpers.R"))
source(file.path(src_root, "data_helpers.R"))
source(file.path(src_root, "simulation_helpers.R"))
source(file.path(src_root, "stan_helpers.R"))
source(file.path(src_root, "diagnostics_helpers.R"))
source(file.path(src_root, "loading_visualization.R"))

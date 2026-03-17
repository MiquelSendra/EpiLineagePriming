script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[1]))
} else {
  normalizePath("figures_scripts/setup_behaviors_environment.R")
}

setwd(dirname(script_path))
source("R/behaviors_helpers.R")

install_missing_behavior_packages()
setup_behavior_project(project_root = "..")

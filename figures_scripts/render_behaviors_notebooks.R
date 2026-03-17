script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[1]))
} else {
  normalizePath("figures_scripts/render_behaviors_notebooks.R")
}

setwd(dirname(script_path))
source("R/behaviors_helpers.R")

install_missing_behavior_packages()

render_notebook <- function(input_file) {
  rmarkdown::render(
    input = input_file,
    output_dir = "output",
    envir = new.env(parent = globalenv())
  )
}

render_notebook("Behaviors_clean.Rmd")
render_notebook("Behaviors_explainer.Rmd")

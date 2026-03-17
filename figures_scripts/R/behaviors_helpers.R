required_behavior_packages <- c(
  celltrackR = "CRAN",
  circlize = "CRAN",
  ComplexHeatmap = "BIOC",
  DataExplorer = "CRAN",
  dplyr = "CRAN",
  fractaldim = "CRAN",
  ggplot2 = "CRAN",
  gridExtra = "CRAN",
  gtools = "CRAN",
  kableExtra = "CRAN",
  lme4 = "CRAN",
  lmerTest = "CRAN",
  pals = "CRAN",
  patchwork = "CRAN",
  plyr = "CRAN",
  plotly = "CRAN",
  pracma = "CRAN",
  RColorBrewer = "CRAN",
  readxl = "CRAN",
  scales = "CRAN",
  Seurat = "CRAN",
  stringr = "CRAN",
  tibble = "CRAN",
  tidyr = "CRAN",
  zoo = "CRAN"
)

install_missing_behavior_packages <- function(package_sources = required_behavior_packages) {
  repos <- getOption("repos")
  cran_repo <- unname(repos["CRAN"])
  if (is.null(repos) || length(cran_repo) == 0 || is.na(cran_repo) || cran_repo %in% c("@CRAN@", "")) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }

  installed <- rownames(installed.packages())
  missing <- setdiff(names(package_sources), installed)
  if (!length(missing)) {
    return(invisible(missing))
  }

  cran_packages <- missing[package_sources[missing] == "CRAN"]
  bioc_packages <- missing[package_sources[missing] == "BIOC"]

  if (length(cran_packages)) {
    install.packages(cran_packages)
  }

  if (length(bioc_packages)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(bioc_packages, ask = FALSE, update = FALSE)
  }

  invisible(missing)
}

load_behavior_packages <- function(package_sources = required_behavior_packages) {
  invisible(lapply(names(package_sources), function(package_name) {
    suppressPackageStartupMessages(
      library(package_name, character.only = TRUE)
    )
  }))
}

rgl_is_available <- function() {
  isTRUE(tryCatch({
    suppressPackageStartupMessages(library(rgl))
    TRUE
  }, error = function(...) FALSE))
}

plot_tracks_3d <- function(track_list, tracking_df, cell_type, rgl_available = FALSE) {
  if (rgl_available) {
    plot3d(track_list[[cell_type]], main = paste(cell_type, "tracks"), tick.marks = FALSE)
    return(invisible(NULL))
  }

  plot_df <- tracking_df |>
    dplyr::filter(cell_type == !!cell_type) |>
    dplyr::arrange(node_id, timepoint) |>
    dplyr::mutate(
      hover_label = paste0(
        "node_id: ", node_id,
        "<br>timepoint: ", timepoint,
        "<br>X: ", round(X, 2),
        "<br>Y: ", round(Y, 2),
        "<br>Z: ", round(Z, 2)
      )
    )

  plotly::plot_ly(
    data = plot_df,
    x = ~X,
    y = ~Y,
    z = ~Z,
    split = ~factor(node_id),
    type = "scatter3d",
    mode = "lines",
    text = ~hover_label,
    hoverinfo = "text",
    showlegend = FALSE
  ) |>
    plotly::layout(
      title = paste(cell_type, "tracks"),
      scene = list(
        xaxis = list(title = "X"),
        yaxis = list(title = "Y"),
        zaxis = list(title = "Z")
      )
    )
}

setup_behavior_project <- function(project_root = "..") {
  repos <- getOption("repos")
  cran_repo <- unname(repos["CRAN"])
  if (is.null(repos) || length(cran_repo) == 0 || is.na(cran_repo) || cran_repo %in% c("@CRAN@", "")) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }

  if (!requireNamespace("renv", quietly = TRUE)) {
    install.packages("renv")
  }

  project_root <- normalizePath(project_root, mustWork = TRUE)

  if (!file.exists(file.path(project_root, "renv.lock"))) {
    renv::init(project = project_root, bare = TRUE)
  } else {
    renv::activate(project = project_root)
  }

  renv::settings$snapshot.type("all", project = project_root)
  renv::snapshot(project = project_root, prompt = FALSE)

  invisible(project_root)
}

behavior_paths <- function(analysis_dir = ".", file_stem = "KM1_tracking_extended_withDists.csv") {
  analysis_dir <- normalizePath(analysis_dir, mustWork = TRUE)
  data_dir <- normalizePath(file.path(analysis_dir, "..", "data"), mustWork = TRUE)
  output_root <- file.path(analysis_dir, "output", "behaviors_v16")
  cache_dir <- file.path(output_root, "cache")
  figure_dir <- file.path(output_root, "plots")

  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  list(
    analysis_dir = analysis_dir,
    data_dir = data_dir,
    output_root = output_root,
    cache_dir = cache_dir,
    figure_dir = figure_dir,
    file_stem = file_stem,
    alltracks_csv = file.path(data_dir, paste0(file_stem, "_alltracks.csv")),
    track_details_xlsx = file.path(data_dir, "KM1_tracking.xlsx"),
    legacy_clean_wide = file.path(data_dir, paste0(file_stem, "_final_table_timepoints_TOT_wide_clean_nodups.csv")),
    legacy_stats = file.path(data_dir, "stats", "stat_lmer_ns5.csv")
  )
}

configure_knitr_for_behaviors <- function(paths) {
  knitr::opts_chunk$set(
    echo = TRUE,
    message = FALSE,
    warning = FALSE,
    fig.align = "center",
    fig.width = 8,
    fig.height = 5,
    dpi = 120,
    out.width = "100%",
    comment = NA,
    fig.path = paste0(normalizePath(paths$figure_dir, mustWork = TRUE), "/")
  )
}

read_tracking_table <- function(alltracks_csv) {
  tracking_df <- read.csv(alltracks_csv, check.names = FALSE)

  if (names(tracking_df)[1] %in% c("X", "X.1", "")) {
    tracking_df <- tracking_df[, -1, drop = FALSE]
  }

  expected_columns <- c(
    "node_id", "track_name", "X", "Y", "Z",
    "cell_type", "timepoint", "DistBorder"
  )
  missing_columns <- setdiff(expected_columns, names(tracking_df))
  if (length(missing_columns)) {
    stop("Missing expected columns in tracking table: ", paste(missing_columns, collapse = ", "))
  }

  tracking_df <- tracking_df[, expected_columns]
  tracking_df$node_id <- as.numeric(tracking_df$node_id)
  tracking_df$timepoint <- as.numeric(tracking_df$timepoint)
  tracking_df$DistBorder <- as.numeric(tracking_df$DistBorder)
  tracking_df$cell_type <- as.character(tracking_df$cell_type)

  tracking_df
}

build_track_list <- function(tracking_df) {
  cell_types <- unique(tracking_df$cell_type)
  track_list <- lapply(cell_types, function(cell_type) {
    celltrackR::as.tracks(
      tracking_df[tracking_df$cell_type == cell_type, ],
      id.column = 1,
      time.column = 7,
      pos.columns = 3:5
    )
  })
  names(track_list) <- cell_types
  track_list
}

resolve_duplicate_timepoints <- function(tracking_df) {
  duplicated_rows <- duplicated(tracking_df[, c("node_id", "timepoint")]) |
    duplicated(tracking_df[, c("node_id", "timepoint")], fromLast = TRUE)

  duplicate_table <- tracking_df[duplicated_rows, , drop = FALSE]
  corrected_table <- tracking_df[!duplicated(tracking_df[, c("node_id", "timepoint")]), , drop = FALSE]
  rownames(corrected_table) <- NULL

  list(
    duplicate_rows = duplicate_table,
    corrected_tracks = corrected_table
  )
}

measure_subtrack <- function(track) {
  c(
    trackLength = celltrackR::trackLength(track),
    maxDisplacement = celltrackR::maxDisplacement(track),
    speed = celltrackR::speed(track),
    displacement = celltrackR::displacement(track),
    squareDisplacement = celltrackR::squareDisplacement(track),
    displacementRatio = celltrackR::displacementRatio(track),
    outreachRatio = celltrackR::outreachRatio(track),
    straightness = celltrackR::straightness(track),
    asphericity = celltrackR::asphericity(track),
    overallAngle = pracma::rad2deg(celltrackR::overallAngle(track)),
    meanTurningAngle = pracma::rad2deg(celltrackR::meanTurningAngle(track)),
    overallDot = celltrackR::overallDot(track),
    overallNormDot = celltrackR::overallNormDot(track),
    fractalDimension = celltrackR::fractalDimension(track)
  )
}

compute_measurements_long <- function(track_list, corrected_tracks, windows) {
  measurement_tables <- list()
  distborder_tables <- list()

  for (window_size in windows) {
    message("Computing measurements for window ", window_size)

    window_tables <- lapply(names(track_list), function(cell_type) {
      subtrack_list <- celltrackR::subtracks(track_list[[cell_type]], window_size)
      measurement_matrix <- sapply(subtrack_list, measure_subtrack)

      if (is.null(dim(measurement_matrix))) {
        measurement_matrix <- matrix(
          measurement_matrix,
          ncol = 1,
          dimnames = list(names(measure_subtrack(track_list[[cell_type]][[1]])), names(subtrack_list)[1])
        )
      }

      measurement_table <- as.data.frame(t(measurement_matrix), stringsAsFactors = FALSE)
      measurement_table <- tibble::rownames_to_column(measurement_table, "track_time_key")
      measurement_table <- tidyr::separate(
        measurement_table,
        track_time_key,
        into = c("node_id", "t"),
        sep = "[^[:alnum:]]+",
        convert = TRUE
      )

      observed_timepoints <- unlist(lapply(track_list[[cell_type]], function(track) {
        if (nrow(track) > window_size) {
          track[seq_len(nrow(track) - window_size), "t"]
        } else {
          numeric(0)
        }
      }))

      measurement_table$timepoint <- observed_timepoints
      measurement_table$cell_type <- cell_type
      measurement_table$window <- window_size
      measurement_table <- dplyr::select(
        measurement_table,
        window, cell_type, node_id, timepoint, t, dplyr::everything()
      )

      measurement_table
    })

    measurement_tables[[as.character(window_size)]] <- dplyr::bind_rows(window_tables)

    distborder_window <- corrected_tracks |>
      dplyr::group_by(cell_type, node_id) |>
      dplyr::arrange(node_id, timepoint, .by_group = TRUE) |>
      dplyr::filter(dplyr::n() >= window_size) |>
      dplyr::mutate(
        window = window_size,
        DistBorderMean = c(zoo::rollmean(DistBorder, window_size), rep(NA_real_, window_size - 1))
      ) |>
      dplyr::ungroup() |>
      dplyr::filter(!is.na(DistBorderMean)) |>
      dplyr::select(window, node_id, cell_type, timepoint, DistBorderMean)

    distborder_tables[[as.character(window_size)]] <- distborder_window
  }

  measurement_table <- dplyr::bind_rows(measurement_tables)
  distborder_table <- dplyr::bind_rows(distborder_tables)

  dplyr::inner_join(
    measurement_table,
    distborder_table,
    by = c("window", "node_id", "cell_type", "timepoint")
  )
}

pivot_measurements_wide <- function(measurements_long) {
  value_columns <- setdiff(
    names(measurements_long),
    c("window", "cell_type", "node_id", "timepoint", "t")
  )

  measurements_long |>
    dplyr::arrange(window, cell_type, node_id, timepoint) |>
    tidyr::pivot_wider(
      names_from = window,
      values_from = dplyr::all_of(value_columns),
      names_glue = "{.value}_win{window}"
    )
}

clean_measurements_wide <- function(measurements_wide, max_missing_fraction = 0.16) {
  non_empty <- measurements_wide[, colSums(is.na(measurements_wide)) < nrow(measurements_wide), drop = FALSE]
  filtered <- non_empty[, colSums(is.na(non_empty)) < max_missing_fraction * nrow(non_empty), drop = FALSE]
  filtered[stats::complete.cases(filtered), , drop = FALSE]
}

collapse_duplicate_feature_rows <- function(clean_wide) {
  duplicated_rows <- duplicated(clean_wide[, 3:ncol(clean_wide)]) |
    duplicated(clean_wide[, 3:ncol(clean_wide)], fromLast = TRUE)

  duplicate_table <- clean_wide[duplicated_rows, , drop = FALSE]

  collapsed_duplicates <- duplicate_table |>
    dplyr::group_by(dplyr::across(-node_id)) |>
    dplyr::summarise(node_id = sum(node_id), .groups = "drop") |>
    dplyr::relocate(cell_type, node_id, timepoint, t)

  deduplicated <- dplyr::bind_rows(clean_wide[!duplicated_rows, , drop = FALSE], collapsed_duplicates)
  deduplicated <- deduplicated |>
    dplyr::mutate(
      cell_type = factor(cell_type),
      node_id = factor(node_id),
      timepoint = factor(timepoint),
      t = factor(t)
    )

  list(
    duplicate_rows = duplicate_table,
    clustering_table = deduplicated
  )
}

fit_longitudinal_models <- function(
  measurements_wide,
  time_var = "timepoint",
  effect_var = "cell_type",
  id_var = "node_id",
  response_vars = setdiff(names(measurements_wide), c("cell_type", "node_id", "timepoint", "t")),
  spline_df = 5
) {
  model_results <- lapply(seq_along(response_vars), function(index) {
    response_var <- response_vars[index]
    if (index %% 10 == 1 || index == length(response_vars)) {
      message("Fitting longitudinal model ", index, " / ", length(response_vars), ": ", response_var)
    }

    model_df <- measurements_wide |>
      dplyr::select(dplyr::all_of(c(effect_var, id_var, time_var, response_var))) |>
      dplyr::filter(is.finite(.data[[response_var]])) |>
      dplyr::mutate(
        "{effect_var}" := factor(.data[[effect_var]], levels = c("CM", "ET")),
        "{id_var}" := factor(.data[[id_var]])
      )

    if (nrow(model_df) < 10 || dplyr::n_distinct(model_df[[effect_var]]) < 2) {
      return(NULL)
    }

    model_df$response_value <- as.numeric(scale(model_df[[response_var]]))

    model_formula <- stats::as.formula(
      paste0(
        "response_value ~ ",
        effect_var,
        " * ns(",
        time_var,
        ", ",
        spline_df,
        ") + (1 + ",
        time_var,
        " || ",
        id_var,
        ")"
      )
    )

    fitted_model <- lmerTest::lmer(model_formula, data = model_df, REML = FALSE)
    coefficient_table <- as.data.frame(summary(fitted_model)$coefficients)
    coefficient_table$Effect <- rownames(coefficient_table)
    coefficient_table <- coefficient_table[grepl(paste0("^", effect_var), coefficient_table$Effect), , drop = FALSE]

    if (!nrow(coefficient_table)) {
      return(NULL)
    }

    coefficient_table |>
      dplyr::transmute(
        Param = response_var,
        Effect = Effect,
        Nobs = nrow(model_df),
        Nsamp = dplyr::n_distinct(.data[[id_var]]),
        Coef = Estimate,
        CI_low = Estimate - 1.96 * `Std. Error`,
        CI_upp = Estimate + 1.96 * `Std. Error`,
        pval = `Pr(>|t|)`,
        Model = "lmer"
      )
  })

  dplyr::bind_rows(model_results) |>
    dplyr::mutate(pval.adj = p.adjust(pval, method = "holm"))
}

derive_cell_type <- function(track_name) {
  track_name_lower <- tolower(track_name)
  ifelse(
    stringr::str_detect(track_name_lower, "ec|et|aorta|uc"),
    "ET",
    "CM"
  )
}

derive_cell_type_detailed <- function(track_name) {
  track_name_lower <- tolower(track_name)

  detailed_label <- ifelse(
    stringr::str_detect(track_name_lower, "ec|et|aorta|uc"),
    "EC_alone_or_with_EmE-EmEt",
    "CM_alone_or_with_EmE"
  )

  detailed_label[track_name_lower %in% c("partial_jcf1", "mes")] <- "EmE_without_ET_nor_CM"
  detailed_label[stringr::str_detect(
    track_name_lower,
    "full_et_aorta_noetmes1|partial_et_aorta|partial_uc2|partial_uc3"
  )] <- "EmET_noEC"

  detailed_label
}

load_track_details_metadata <- function(paths) {
  if (!file.exists(paths$track_details_xlsx)) {
    return(NULL)
  }

  detailed_tracks <- readxl::read_excel(paths$track_details_xlsx)
  detailed_tracks <- as.data.frame(detailed_tracks)

  if (!all(c("node_id", "track_name", "cell_type", "cell_type_detailed") %in% names(detailed_tracks))) {
    return(NULL)
  }

  detailed_tracks |>
    dplyr::transmute(
      node_id_orig = as.numeric(node_id),
      track_name,
      cell_type,
      cell_type_detailed
    )
}

build_node_id_correspondence <- function(alltracks_df, clean_wide_df, detailed_tracks_df = NULL) {
  unique_positions <- alltracks_df[!duplicated(alltracks_df[, c("node_id", "timepoint")]), , drop = FALSE]
  rownames(unique_positions) <- NULL

  collapsed_ids <- clean_wide_df |>
    dplyr::group_by(dplyr::across(-node_id)) |>
    dplyr::mutate(
      node_id_sum = sum(node_id),
      N_rep = dplyr::n_distinct(node_id)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(node_id, timepoint, node_id_sum, N_rep)

  mapping_df <- unique_positions |>
    dplyr::left_join(collapsed_ids, by = c("node_id", "timepoint")) |>
    dplyr::mutate(
      node_id_orig = node_id
    )

  if (!is.null(detailed_tracks_df)) {
    mapping_df <- mapping_df |>
      dplyr::left_join(
        detailed_tracks_df,
        by = c("node_id_orig", "track_name"),
        suffix = c("", ".detailed")
      ) |>
      dplyr::mutate(
        cell_type = dplyr::coalesce(cell_type.detailed, cell_type),
        cell_type_detailed = dplyr::coalesce(cell_type_detailed, derive_cell_type_detailed(track_name))
      ) |>
      dplyr::select(-dplyr::any_of("cell_type.detailed"))
  } else {
    mapping_df <- mapping_df |>
      dplyr::mutate(
        cell_type = derive_cell_type(track_name),
        cell_type_detailed = derive_cell_type_detailed(track_name)
      )
  }

  mapping_df <- mapping_df |>
    dplyr::mutate(
      node_id = as.factor(node_id_sum),
      timepoint = as.factor(timepoint)
    ) |>
    dplyr::select(
      node_id, node_id_orig, timepoint, track_name, cell_type,
      cell_type_detailed, dplyr::everything()
    ) |>
    dplyr::arrange(dplyr::desc(N_rep), dplyr::desc(node_id), timepoint, node_id_orig)

  mapping_df
}

add_annotations_to_seurat <- function(seurat_object, mapping_df, ntimepoints = 5) {
  metadata_df <- seurat_object@meta.data
  metadata_df$node_id <- as.numeric(as.character(metadata_df$node_id))
  metadata_df$timepoint <- as.numeric(as.character(metadata_df$timepoint))
  mapping_df$node_id <- as.numeric(as.character(mapping_df$node_id))
  mapping_df$timepoint <- as.numeric(as.character(mapping_df$timepoint))

  annotated_metadata <- dplyr::left_join(
    metadata_df,
    mapping_df,
    by = c("node_id", "timepoint", "track_name", "cell_type", "X", "Y", "Z", "DistBorder"),
    multiple = "first"
  )

  seurat_object@meta.data <- annotated_metadata
  seurat_object$Embryo <- "EMB1"
  seurat_object$node_id2 <- as.numeric(as.character(seurat_object$node_id))
  seurat_object$ntimepoint <- as.numeric(as.character(seurat_object$timepoint))
  seurat_object$intervalt <- cut(seurat_object$ntimepoint, ntimepoints)

  seurat_object
}

scrollable_kable <- function(data_frame, caption, height = "400px", digits = 3) {
  knitr::kable(
    data_frame,
    caption = caption,
    digits = digits,
    row.names = FALSE,
    booktabs = TRUE,
    longtable = TRUE
  ) |>
    kableExtra::kable_styling(
      bootstrap_options = c("striped", "hover", "condensed", "responsive"),
      fixed_thead = TRUE,
      position = "center"
    ) |>
    kableExtra::scroll_box(height = height)
}

load_or_create_rds <- function(cache_path, builder) {
  if (file.exists(cache_path)) {
    readRDS(cache_path)
  } else {
    object <- builder()
    saveRDS(object, cache_path)
    object
  }
}

load_or_create_csv <- function(cache_path, builder) {
  if (file.exists(cache_path)) {
    read.csv(cache_path, check.names = FALSE)
  } else {
    object <- builder()
    write.csv(object, cache_path, row.names = FALSE)
    object
  }
}

compare_with_reference <- function(current_df, reference_path) {
  if (!file.exists(reference_path)) {
    return(
      data.frame(
        reference_exists = FALSE,
        same_nrow = NA,
        same_ncol = NA,
        same_columns = NA
      )
    )
  }

  reference_df <- read.csv(reference_path, check.names = FALSE)
  if (names(reference_df)[1] %in% c("X", "X.1", "")) {
    reference_df <- reference_df[, -1, drop = FALSE]
  }

  data.frame(
    reference_exists = TRUE,
    same_nrow = identical(nrow(current_df), nrow(reference_df)),
    same_ncol = identical(ncol(current_df), ncol(reference_df)),
    same_columns = identical(names(current_df), names(reference_df))
  )
}

estimate_elbow_point <- function(values) {
  if (length(values) < 3) {
    return(seq_along(values)[1])
  }

  x_coords <- seq_along(values)
  first_point <- c(x_coords[1], values[1])
  last_point <- c(x_coords[length(values)], values[length(values)])
  line_vector <- last_point - first_point
  line_length <- sqrt(sum(line_vector^2))

  distances <- vapply(seq_along(values), function(index) {
    point <- c(x_coords[index], values[index])
    numerator <- abs(
      line_vector[2] * point[1] -
        line_vector[1] * point[2] +
        last_point[1] * first_point[2] -
        last_point[2] * first_point[1]
    )
    numerator / line_length
  }, numeric(1))

  which.max(distances)
}

build_behavior_analysis_state <- function(
  paths,
  windows = c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100),
  spline_df = 5,
  reuse_legacy_clean_wide = TRUE,
  reuse_legacy_stats = TRUE
) {
  alltracks_df <- read_tracking_table(paths$alltracks_csv)
  detailed_tracks_df <- load_track_details_metadata(paths)

  track_list <- load_or_create_rds(
    file.path(paths$cache_dir, "track_list_all.rds"),
    function() build_track_list(alltracks_df)
  )

  duplicate_resolution <- load_or_create_rds(
    file.path(paths$cache_dir, "duplicate_resolution.rds"),
    function() resolve_duplicate_timepoints(alltracks_df)
  )

  corrected_track_list <- load_or_create_rds(
    file.path(paths$cache_dir, "track_list_corrected.rds"),
    function() build_track_list(duplicate_resolution$corrected_tracks)
  )

  if (reuse_legacy_clean_wide && file.exists(paths$legacy_clean_wide)) {
    clean_wide <- read.csv(paths$legacy_clean_wide, check.names = FALSE)
    if (names(clean_wide)[1] %in% c("X", "X.1", "")) {
      clean_wide <- clean_wide[, -1, drop = FALSE]
    }
    measurements_long <- NULL
    measurements_wide <- clean_wide
    collapsed_duplicates <- list(
      duplicate_rows = data.frame(),
      clustering_table = clean_wide |>
        dplyr::mutate(
          cell_type = factor(cell_type),
          node_id = factor(node_id),
          timepoint = factor(timepoint),
          t = factor(t)
        )
    )
  } else {
    measurements_long <- load_or_create_rds(
      file.path(paths$cache_dir, "measurements_long.rds"),
      function() compute_measurements_long(
        track_list = corrected_track_list,
        corrected_tracks = duplicate_resolution$corrected_tracks,
        windows = windows
      )
    )

    measurements_wide <- load_or_create_rds(
      file.path(paths$cache_dir, "measurements_wide.rds"),
      function() pivot_measurements_wide(measurements_long)
    )

    clean_wide <- load_or_create_rds(
      file.path(paths$cache_dir, "measurements_wide_clean.rds"),
      function() clean_measurements_wide(measurements_wide)
    )

    collapsed_duplicates <- load_or_create_rds(
      file.path(paths$cache_dir, "measurements_wide_clustering.rds"),
      function() collapse_duplicate_feature_rows(clean_wide)
    )
  }

  if (reuse_legacy_stats && file.exists(paths$legacy_stats)) {
    model_results <- read.csv(paths$legacy_stats, check.names = FALSE)
  } else {
    model_results <- load_or_create_rds(
      file.path(paths$cache_dir, paste0("longitudinal_models_ns", spline_df, ".rds")),
      function() fit_longitudinal_models(
        measurements_wide = measurements_wide,
        spline_df = spline_df
      )
    )
  }

  node_mapping <- load_or_create_rds(
    file.path(paths$cache_dir, "node_id_mapping.rds"),
    function() build_node_id_correspondence(
      alltracks_df = duplicate_resolution$corrected_tracks,
      clean_wide_df = clean_wide,
      detailed_tracks_df = detailed_tracks_df
    )
  )

  write.csv(
    measurements_wide,
    file.path(paths$output_root, paste0(paths$file_stem, "_final_table_timepoints_TOT_wide.csv")),
    row.names = FALSE
  )
  write.csv(
    collapsed_duplicates$clustering_table,
    file.path(paths$output_root, paste0(paths$file_stem, "_final_table_timepoints_TOT_wide_clean_nodups.csv")),
    row.names = FALSE
  )
  write.csv(
    model_results,
    file.path(paths$output_root, paste0("stat_lmer_ns", spline_df, ".csv")),
    row.names = FALSE
  )
  write.csv(
    node_mapping,
    file.path(paths$output_root, paste0(paths$file_stem, "_alltracks_alloriginal_node_id.csv")),
    row.names = FALSE
  )

  if (!reuse_legacy_clean_wide) {
    write.csv(
      clean_wide,
      file.path(paths$output_root, paste0(paths$file_stem, "_final_table_timepoints_TOT_wide_clean.csv")),
      row.names = FALSE
    )
  }

  list(
    windows = windows,
    spline_df = spline_df,
    alltracks = alltracks_df,
    track_list = track_list,
    duplicate_resolution = duplicate_resolution,
    corrected_track_list = corrected_track_list,
    measurements_long = measurements_long,
    measurements_wide = measurements_wide,
    clean_wide = clean_wide,
    collapsed_duplicates = collapsed_duplicates,
    model_results = model_results,
    node_mapping = node_mapping
  )
}

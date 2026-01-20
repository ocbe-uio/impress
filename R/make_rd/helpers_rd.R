analysis_cohortcd <- function(cfg, default = 2) {
  cohortcd <- cfg$analysis$cohortcd %||% default
  as.character(cohortcd)
}

filter_cohort <- function(data, cfg, cohort_col = "cohortcd") {
  cohortcd <- analysis_cohortcd(cfg)
  data %>%
    dplyr::filter(.data[[cohort_col]] == .env$cohortcd)
}

get_plot_dodge_width <- function(cfg, default = 0.2) {
  width <- cfg$plot$dodge_width
  if (is.null(width) || !is.numeric(width) || length(width) != 1 || is.na(width)) {
    return(default)
  }
  width
}

format_estimate_ci <- function(estimate, conf.low, conf.high, digits = 2) {
  paste0(
    format(round(estimate, digits), nsmall = digits),
    " (",
    format(round(conf.low, digits), nsmall = digits),
    " to ",
    format(round(conf.high, digits), nsmall = digits),
    ")"
  )
}

binary_yesno <- function(x) {
  dplyr::case_when(
    toupper(as.character(x)) == "YES" ~ 1,
    toupper(as.character(x)) == "NO" ~ 0,
    TRUE ~ NA_real_
  )
}

steroid_baseline_map <- function(data, cfg, paramcd) {
  base_data <- data %>%
    filter_cohort(cfg) %>%
    dplyr::filter(.data$paramcd == .env$paramcd)

  base_rows <- base_data %>%
    dplyr::filter(.data$ablfl == "Y")

  if (!nrow(base_rows)) {
    base_rows <- base_data %>%
      dplyr::filter(.data$avisitn == -7)
  }

  if (!nrow(base_rows)) {
    base_rows <- base_data %>%
      dplyr::group_by(.data$usubjid) %>%
      dplyr::filter(.data$avisitn == min(.data$avisitn, na.rm = TRUE)) %>%
      dplyr::ungroup()
  }

  base_rows %>%
    dplyr::mutate(base_flag = binary_yesno(.data$avalc)) %>%
    dplyr::group_by(.data$usubjid) %>%
    dplyr::summarise(
      base_flag = if (any(!is.na(base_flag))) max(base_flag, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )
}

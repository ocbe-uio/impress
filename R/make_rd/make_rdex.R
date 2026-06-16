# Reporting for extent of exposure (SAP §7).
#
# Descriptive summaries on the configured cohort (cohortcd == 2), grouped by planned
# dose (FAS mock-up; _EX.csv dosing is empty, see make_adex.R).
#
# Outputs (all knitr::kable): N/% exposed per planned dose level, treatment duration
# overall and by long-term dose, short-term-only vs continued-to-long-term, switched
# dose / discontinued early, and duration categories.

make_rdex <- function(adex, cfg) {
  if (is.null(adex) || nrow(adex) == 0) {
    return(NULL)
  }

  dt <- adex %>% filter_cohort(cfg)
  if (nrow(dt) == 0) {
    return(NULL)
  }

  Ntot <- nrow(dt)

  # ---- N / % exposed per planned dose level --------------------------------
  exposed_tbl <- function(dose_var, lbl) {
    dt %>%
      dplyr::mutate(dose = factor(.data[[dose_var]], levels = c(0, 25, 50, 100))) %>%
      dplyr::filter(!is.na(dose)) %>%
      dplyr::group_by(dose, .drop = FALSE) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(txt = sprintf("%d (%.1f%%)", n, 100 * n / Ntot)) %>%
      dplyr::transmute(`Dose (mg)` = dose, !!lbl := txt)
  }

  exposed <- dplyr::full_join(
    exposed_tbl("dose01p", "Short-term (randomised)"),
    exposed_tbl("dose02p", "Long-term (maintenance)"),
    by = "Dose (mg)"
  ) %>%
    knitr::kable(caption = NULL)

  # ---- Duration overall and by long-term dose ------------------------------
  dur_summary <- function(df) {
    d <- df$trtdurd[!is.na(df$trtdurd)]
    if (length(d) == 0) {
      return(tibble::tibble(N = 0L, `Median (range), days` = "-", `Mean (SD), days` = "-"))
    }
    tibble::tibble(
      N = length(d),
      `Median (range), days` = sprintf("%.0f (%.0f-%.0f)", stats::median(d), min(d), max(d)),
      `Mean (SD), days`      = sprintf("%.1f (%.1f)", mean(d), stats::sd(d))
    )
  }

  duration_overall <- dur_summary(dt) %>%
    dplyr::mutate(Group = "Overall", .before = 1) %>%
    knitr::kable(caption = NULL)

  duration_by_dose <- dt %>%
    dplyr::mutate(dose = factor(dose02p, levels = c(0, 25, 50, 100))) %>%
    dplyr::filter(!is.na(dose)) %>%
    dplyr::group_by(dose) %>%
    dplyr::group_modify(~ dur_summary(.x)) %>%
    dplyr::ungroup() %>%
    dplyr::rename(`Long-term dose (mg)` = dose) %>%
    knitr::kable(caption = NULL)

  # ---- Short-term-only vs continued, switched, discontinued early ----------
  flag_pct <- function(flag) {
    n <- sum(dt[[flag]] == "Y", na.rm = TRUE)
    sprintf("%d (%.1f%%)", n, 100 * n / Ntot)
  }

  status <- tibble::tibble(
    Measure = c("Continued to long-term treatment",
                "Short-term assignment only",
                "Switched dose (short- vs long-term)",
                "Discontinued treatment early"),
    `N (%)` = c(
      flag_pct("ltfl"),
      sprintf("%d (%.1f%%)", sum(dt$ltfl != "Y", na.rm = TRUE),
              100 * sum(dt$ltfl != "Y", na.rm = TRUE) / Ntot),
      flag_pct("switchfl"),
      flag_pct("discfl")
    )
  ) %>%
    knitr::kable(caption = NULL)

  # ---- Duration categories --------------------------------------------------
  durcat <- dt %>%
    dplyr::filter(!is.na(durcat)) %>%
    dplyr::group_by(durcat, .drop = FALSE) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(`N (%)` = sprintf("%d (%.1f%%)", n, 100 * n / Ntot)) %>%
    dplyr::transmute(`Duration category` = durcat, `N (%)`) %>%
    knitr::kable(caption = NULL)

  list(
    nsubj            = Ntot,
    exposed          = exposed,
    duration_overall = duration_overall,
    duration_by_dose = duration_by_dose,
    status           = status,
    durcat           = durcat
  )
}

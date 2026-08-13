# Reporting for extent of exposure (SAP §7.7).
#
# Descriptive summaries on the configured cohort (cohortcd == 2), grouped by randomised 
# arm and planned dose 
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

  Ntot <- length(unique(dt$subjid))

  # ---- N / % exposed per planned dose level --------------------------------
  exposed_tbl <- function(dose_var, lbl) {
    dt %>%
      group_by(subjid) |> 
      filter(row_number() == 1) |> 
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



 

 exposure_table <- dt |>
   group_by(arm) |> 
  summarise(
    `N`                            = n(),
    `Duration (days), Mean (SD)`          = sprintf("%.1f (%.1f)", mean(totdur, na.rm = TRUE), sd(totdur, na.rm = TRUE)),
    `Duration (days), Median (Min - Max)`  = sprintf("%.1f (%.1f - %.1f)", median(totdur, na.rm = TRUE), min(totdur, na.rm = TRUE), max(totdur, na.rm = TRUE)),
    `Total dose (mg), Mean (SD)`        = sprintf("%.1f (%.1f)", mean(totdose, na.rm = TRUE), sd(totdose, na.rm = TRUE)),
    `Total dose(mg), Median (Min - Max)`= sprintf("%.1f (%.1f - %.1f)", median(totdose, na.rm = TRUE), min(totdose, na.rm = TRUE), max(totdose, na.rm = TRUE)),
    `Daily dose (mg), Mean (SD)`        = sprintf("%.1f (%.1f)", mean(dailydose, na.rm = TRUE), sd(dailydose, na.rm = TRUE)),
    `Daily dose(mg), Median (Min - Max)`= sprintf("%.1f (%.1f - %.1f)", median(dailydose, na.rm = TRUE), min(dailydose, na.rm = TRUE), max(dailydose, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  rename(Arm = arm) |> 
   knitr::kable(
  align = c("l", "c", "c", "c", "c", "c"),
  caption = "Exposure to randomised arm"
  )
exposure_table2 <- dt |>
  mutate(Dose = paste0(dose01p, " mg"),
         Dose = factor(Dose, levels = c("25 mg", "50 mg", "100 mg"))) |> 
   group_by(Dose) |> 
  summarise(
    `N`                            = n(),
    `Daily dose (mg), Mean (SD)`        = sprintf("%.1f (%.1f)", mean(dailydose, na.rm = TRUE), sd(dailydose, na.rm = TRUE)),
    `Daily dose(mg), Median (Min - Max)`= sprintf("%.1f (%.1f - %.1f)", median(dailydose, na.rm = TRUE), min(dailydose, na.rm = TRUE), max(dailydose, na.rm = TRUE)),
    .groups = "drop"
  ) |>
   knitr::kable(
  align = c("l", "c", "c", "c", "c", "c"),
  caption = "Exposure to treatment (mean daily dose) by dose"
  )
dose_mod_table <- dt |>
  distinct(arm, subjid, moddose) |>   # adjust if moddose is already one row per subject
  count(arm, moddose) |>
  group_by(arm) |>
  mutate(pct = n / sum(n) * 100, 
         N = sum(n)) |>
  ungroup() |>
  filter(moddose == TRUE) |>
  transmute(
    Arm = arm,
    N = N,
    `Dose modified, n (%)` = sprintf("%d (%.1f%%)", n, pct)
  )

dosemodtbl <- knitr::kable(
  dose_mod_table,
  align = c("l", "c"),
  caption = "Dose modification by treatment arm"
)

 
  list(
    nsubj            = Ntot,
    exposed          = exposed,
    exposure_table = exposure_table,
    exposure_table2 = exposure_table2,
    dosemodtbl           = dosemodtbl
  )
}

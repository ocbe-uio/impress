# ADVS dataset creation — vital signs (SAP §7)
#
# Source file (active Viedoc export):
#   _VS.csv  — one row per visit; columns VSSYS / VSDIA / VSPULSE (SBP / DBP / pulse)
#
# Output: BDS-style long dataset with one row per subject x parameter x visit,
# mirroring make_adeff.R (visit mapping, baseline on/before first treatment,
# change from baseline, stepped-wedge treatment assignment).
#
# Mock-up note: grouped by planned dose (FAS); see make_adae.R for context.
#
# PARAMCD: SBP / DBP / PULSE   PARAM: Systolic BP / Diastolic BP / Pulse rate

make_advs <- function(raw, adsl, cfg) {

  vs <- raw |> get_raw("vs")
  if (is.null(vs) || nrow(vs) == 0) {
    return(tibble::tibble())
  }

  visits <- visits_from_yaml(cfg)
  visit_levels_df <- visits$map %>%
    dplyr::arrange(avisitn, avisit) %>%
    dplyr::distinct(avisit, avisitn)
  avisit_levels <- visit_levels_df$avisit[!is.na(visit_levels_df$avisit)]
  visit_schedule <- visit_levels_df %>%
    dplyr::filter(!is.na(avisitn), !is.na(avisit)) %>%
    dplyr::arrange(avisitn) %>%
    dplyr::distinct(avisitn, .keep_all = TRUE)

  param_map <- tibble::tribble(
    ~col,       ~paramcd, ~param,          ~avalu,
    "vssys",    "SBP",    "Systolic BP",   "mmHg",
    "vsdia",    "DBP",    "Diastolic BP",  "mmHg",
    "vspulse",  "PULSE",  "Pulse rate",    "beats/min"
  )

  build_long <- function(col, paramcd, param, avalu) {
    vs %>%
      dplyr::mutate(
        subjid = subjectid,
        adt = suppressWarnings(as.Date(vsdat)),
        adt = dplyr::coalesce(adt, as.Date(eventdate))
      ) %>%
      dplyr::filter(!is.na(.data[[col]])) %>%
      dplyr::transmute(
        subjid, eventname, adt,
        paramcd = paramcd, param = param,
        aval = as.numeric(.data[[col]]), avalu = avalu
      )
  }

  advs_long <- purrr::pmap_dfr(param_map, build_long)
  if (!nrow(advs_long)) {
    return(tibble::tibble())
  }

  advs_base <- advs_long %>%
    dplyr::left_join(adsl, by = "subjid") %>%
    dplyr::mutate(
      studyid = cfg$studyid,
      domain  = "ADVS",
      avisit  = lookup_avisit(eventname, visits$map, fallback_title = TRUE),
      avisitn = lookup_avisitn(eventname, visits$map, unmapped = visits$defaults$avisitn_unmapped)
    )

  advs <- admiral::derive_vars_dy(
    dataset        = advs_base %>% dplyr::rename_with(toupper),
    reference_date = RFSTDT,
    source_vars    = admiral::exprs(ADT)
  ) %>%
    dplyr::rename_with(tolower) %>%
    dplyr::arrange(usubjid, paramcd, adt)

  # Baseline restricted to observations on/before first treatment (RFSTDT)
  advs <- advs %>%
    dplyr::group_by(usubjid, paramcd) %>%
    dplyr::mutate(
      baseline_date = {
        eligible <- !is.na(aval) & !is.na(rfstdt) & !is.na(adt) & adt <= rfstdt
        if (any(eligible)) max(adt[eligible], na.rm = TRUE) else as.Date(NA)
      },
      ablfl = dplyr::if_else(!is.na(adt) & !is.na(baseline_date) & adt == baseline_date,
                             "Y", "N", missing = "N"),
      base_raw = dplyr::if_else(ablfl == "Y", aval, NA_real_),
      base = if (any(!is.na(base_raw))) max(base_raw, na.rm = TRUE) else NA_real_,
      chg  = dplyr::if_else(!is.na(aval) & !is.na(base), aval - base, NA_real_),
      pchg = dplyr::if_else(!is.na(aval) & !is.na(base) & base != 0,
                            100 * (aval - base) / base, NA_real_)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-baseline_date, -base_raw)

  # Remap unmapped visits to the nearest scheduled visit by study day
  remap_needed <- (is.na(advs$avisitn) | advs$avisitn == visits$defaults$avisitn_unmapped) &
    !is.na(advs$ady)
  closest_visits <- nearest_visit_from_ady(advs$ady, visit_schedule)
  advs <- advs %>%
    dplyr::mutate(
      avisitn = dplyr::if_else(remap_needed, closest_visits$avisitn, avisitn),
      avisit  = dplyr::if_else(remap_needed, closest_visits$avisit, avisit)
    )

  # Stepped-wedge treatment assignment by visit day
  trt_map <- trt_map_from_yaml(cfg)
  advs <- advs %>%
    dplyr::mutate(
      trtcd = trt_from_map(armcd, avisitn, trt_map),
      trtcd = factor(trtcd, levels = c(0, 25, 50, 100)),
      trt   = factor(paste0(trtcd, " mg"), levels = c("0 mg", "25 mg", "50 mg", "100 mg"))
    )

  param_levels <- param_map$param
  paramcd_levels <- param_map$paramcd

  advs %>%
    dplyr::mutate(
      avisit  = factor(avisit, levels = add_observed_levels(avisit, avisit_levels)),
      paramcd = factor(paramcd, levels = add_observed_levels(paramcd, paramcd_levels)),
      param   = factor(param, levels = add_observed_levels(param, param_levels)),
      ablfl   = factor(ablfl, levels = c("N", "Y"))
    ) %>%
    dplyr::select(
      studyid, domain, usubjid, subjid,
      age, sex, cohort, cohortcd, arm, armcd,
      step, dose01p, dose02p, trt01p, trt02p, randdt, rfstdt,
      trt, trtcd,
      paramcd, param,
      aval, avalu, ablfl, base, chg, pchg,
      avisit, avisitn, adt, ady
    )
}

# Standalone execution -------------------------------------------------------
if (sys.nframe() == 0L) {
  suppressMessages({
    library(dplyr); library(tidyr); library(stringr); library(purrr)
    library(tibble); library(lubridate); library(yaml); library(admiral)
    library(targets)
  })
  source("R/external/functions.R")
  source("R/make_ad/helpers_ad.R")

  cfg <- yaml::read_yaml("config/cfg.yml")
  tar_load(dat); tar_load(adsl)
  advs <- make_advs(dat, adsl, cfg)
  message("advs: ", nrow(advs), " rows, ",
          dplyr::n_distinct(advs$usubjid), " subjects")
  print(advs %>% dplyr::count(paramcd))
}

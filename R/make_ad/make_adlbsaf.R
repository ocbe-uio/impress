# ADLBSAF dataset creation — conventional safety laboratory parameters (SAP §7)
#
# Source file (active Viedoc export):
#   _LBE.csv  — local safety-lab form. In this export only the creatinine block
#               (lbcr*) is populated; the potassium (lbk*) and urea (lbu*) blocks
#               are empty, so the mock-up implements creatinine only (hard-coded).
#               Reference ranges are carried in-data (lbcrr_lower / lbcrr_upper).
#
# This dataset is kept deliberately separate from the cytokine biomarker `adlb`
# (Bio-Plex panels); it is the SAP §7 "conventional clinical laboratory parameters"
# safety dataset (domain "ADLB", PARCAT1 = "SAFETY CHEM").
#
# BDS-style long, mirroring make_advs.R for visit/baseline/treatment derivation.
#
# PARAMCD CREAT (Creatinine). ANRLO/ANRHI reference limits; ANRIND LOW/NORMAL/HIGH.

make_adlbsaf <- function(raw, adsl, cfg) {

  lbe <- raw |> get_raw("lbe")
  if (is.null(lbe) || nrow(lbe) == 0) {
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

  # ---- Creatinine block (lbcr*) ----------------------------------------------
  lbsaf_long <- lbe %>%
    dplyr::mutate(
      subjid = subjectid,
      adt = suppressWarnings(as.Date(lbbdat)),
      adt = dplyr::coalesce(adt, as.Date(eventdate))
    ) %>%
    dplyr::filter(!is.na(suppressWarnings(as.numeric(lbcrres)))) %>%
    dplyr::transmute(
      subjid, eventname, adt,
      paramcd = "CREAT", param = "Creatinine", parcat1 = "SAFETY CHEM",
      aval  = as.numeric(lbcrres),
      avalu = as.character(lbcru),
      anrlo = suppressWarnings(as.numeric(lbcrr_lower)),
      anrhi = suppressWarnings(as.numeric(lbcrr_upper)),
      lbclsig = as.character(lbcrcs)   # investigator clinical-significance flag
    )

  if (!nrow(lbsaf_long)) {
    return(tibble::tibble())
  }

  lbsaf_base <- lbsaf_long %>%
    dplyr::left_join(adsl, by = "subjid") %>%
    dplyr::mutate(
      studyid = cfg$studyid,
      domain  = "ADLB",
      usubjid = dplyr::coalesce(usubjid, paste0(cfg$studyid, "-", subjid)),
      avisit  = lookup_avisit(eventname, visits$map, fallback_title = TRUE),
      avisitn = lookup_avisitn(eventname, visits$map, unmapped = visits$defaults$avisitn_unmapped),
      # Reference-range indicator
      anrind  = dplyr::case_when(
        is.na(aval)                          ~ NA_character_,
        !is.na(anrlo) & aval < anrlo         ~ "LOW",
        !is.na(anrhi) & aval > anrhi         ~ "HIGH",
        TRUE                                 ~ "NORMAL"
      ),
      anrind  = factor(anrind, levels = c("LOW", "NORMAL", "HIGH"))
    )

  lbsaf <- admiral::derive_vars_dy(
    dataset        = lbsaf_base %>% dplyr::rename_with(toupper),
    reference_date = RFSTDT,
    source_vars    = admiral::exprs(ADT)
  ) %>%
    dplyr::rename_with(tolower) %>%
    dplyr::arrange(usubjid, paramcd, adt)

  # Baseline on/before first treatment (RFSTDT)
  lbsaf <- lbsaf %>%
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

  # Remap unmapped visits to nearest scheduled visit by study day
  remap_needed <- (is.na(lbsaf$avisitn) | lbsaf$avisitn == visits$defaults$avisitn_unmapped) &
    !is.na(lbsaf$ady)
  closest_visits <- nearest_visit_from_ady(lbsaf$ady, visit_schedule)
  lbsaf <- lbsaf %>%
    dplyr::mutate(
      avisitn = dplyr::if_else(remap_needed, closest_visits$avisitn, avisitn),
      avisit  = dplyr::if_else(remap_needed, closest_visits$avisit, avisit)
    )

  # Stepped-wedge treatment assignment by visit day
  trt_map <- trt_map_from_yaml(cfg)
  lbsaf <- lbsaf %>%
    dplyr::mutate(
      trtcd = trt_from_map(armcd, avisitn, trt_map),
      trtcd = factor(trtcd, levels = c(0, 25, 50, 100)),
      trt   = factor(paste0(trtcd, " mg"), levels = c("0 mg", "25 mg", "50 mg", "100 mg"))
    )

  lbsaf %>%
    dplyr::mutate(
      avisit = factor(avisit, levels = add_observed_levels(avisit, avisit_levels)),
      ablfl  = factor(ablfl, levels = c("N", "Y"))
    ) %>%
    dplyr::select(
      studyid, domain, usubjid, subjid,
      age, sex, cohort, cohortcd, arm, armcd,
      step, dose01p, dose02p, trt01p, trt02p, randdt, rfstdt,
      trt, trtcd,
      paramcd, param, parcat1,
      aval, avalu, anrlo, anrhi, anrind, lbclsig,
      ablfl, base, chg, pchg,
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
  adlbsaf <- make_adlbsaf(dat, adsl, cfg)
  message("adlbsaf: ", nrow(adlbsaf), " rows, ",
          dplyr::n_distinct(adlbsaf$usubjid), " subjects")
  print(adlbsaf %>% dplyr::count(anrind))
}

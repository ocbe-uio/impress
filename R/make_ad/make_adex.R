# ADEX dataset creation — extent of exposure (SAP §7)
#
# IMPORTANT (mock-up): the _EX.csv dosing export is empty for this cohort, so actual
# administered dose/duration are not available. For the mock-up report we derive
# exposure from the *planned* (sham-randomised) dose in ADSL and the end-of-treatment
# / end-of-study dates, using Full Analysis Set definitions. Assessment by treatment
# actually received is deferred to database lock / unblinding, when _EX.csv is populated.
#
# One row per subject (ADSL-like), domain "ADEX".
#
# Key variables:
#   DOSE01P / DOSE02P  ← planned short-term (randomised) and long-term dose (mg)
#   STEP               ← randomisation step
#   TRTSDT             ← first treatment date (= RFSTDT)
#   TRTEDT             ← last treatment date (EOT, else EOS / death / last contact)
#   TRTDURD            ← treatment duration in days (TRTEDT - TRTSDT + 1)
#   DURCAT             ← duration category (<=14 / >14-<=42 / >42 days)
#   LTFL               ← continued to long-term treatment ("Y" if TRTDURD > 43)
#   SWITCHFL           ← planned dose switch between short- and long-term ("Y")
#   DISCFL             ← discontinued treatment early (EOT for a reason other than
#                        disease progression or death)

make_adex <- function(raw, adsl, cfg) {

  eot <- raw |> get_raw("eot")
  eos <- raw |> get_raw("eos")
  sa  <- raw |> get_raw("sa")

  # ---- End-of-treatment date + reason ----------------------------------------
  eot_dates <- if (is.null(eot)) {
    tibble::tibble(subjid = character(), eotdt = as.Date(character()), eotrea = character())
  } else {
    eot %>%
      dplyr::transmute(
        subjid = subjectid,
        eotdt  = suppressWarnings(as.Date(eotdat)),
        eotrea = as.character(eotrea)
      ) %>%
      dplyr::filter(!is.na(subjid)) %>%
      dplyr::group_by(subjid) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
  }

  # ---- End-of-study / death date ---------------------------------------------
  eos_dates <- if (is.null(eos)) {
    tibble::tibble(subjid = character(), eosdt = as.Date(character()), eosdthdt = as.Date(character()))
  } else {
    eos %>%
      dplyr::transmute(
        subjid   = subjectid,
        eosdt    = suppressWarnings(as.Date(eosdat)),
        eosdthdt = suppressWarnings(as.Date(eosdtdat))
      ) %>%
      dplyr::filter(!is.na(subjid)) %>%
      dplyr::group_by(subjid) %>%
      dplyr::summarise(
        eosdt    = if (any(!is.na(eosdt)))    max(eosdt, na.rm = TRUE)    else as.Date(NA),
        eosdthdt = if (any(!is.na(eosdthdt))) min(eosdthdt, na.rm = TRUE) else as.Date(NA),
        .groups  = "drop"
      )
  }

  sa_dates <- if (is.null(sa)) {
    tibble::tibble(subjid = character(), sadthdt = as.Date(character()), salast = as.Date(character()))
  } else {
    sa %>%
      dplyr::transmute(
        subjid  = subjectid,
        sadthdt = suppressWarnings(as.Date(sadtdat)),
        salast  = suppressWarnings(as.Date(saendat))
      ) %>%
      dplyr::filter(!is.na(subjid)) %>%
      dplyr::group_by(subjid) %>%
      dplyr::summarise(
        sadthdt = if (any(!is.na(sadthdt))) min(sadthdt, na.rm = TRUE) else as.Date(NA),
        salast  = if (any(!is.na(salast)))  max(salast,  na.rm = TRUE) else as.Date(NA),
        .groups = "drop"
      )
  }

  prog_terms <- c("Progressive disease according to RANO criteria",
                  "Extracranial progression")

  adex <- adsl %>%
    dplyr::select(
      studyid, usubjid, subjid, cohort, cohortcd, arm, armcd,
      step, dose01p, dose01pu, dose02p, dose02pu, trt01p, trt02p, randdt, rfstdt
    ) %>%
    dplyr::left_join(eot_dates, by = "subjid") %>%
    dplyr::left_join(eos_dates, by = "subjid") %>%
    dplyr::left_join(sa_dates,  by = "subjid") %>%
    dplyr::mutate(
      domain  = "ADEX",
      deathdt = dplyr::coalesce(eosdthdt, sadthdt),
      trtsdt  = dplyr::coalesce(rfstdt, randdt),
      # Last treatment date: EOT preferred; else death / study end / last contact
      trtedt  = dplyr::coalesce(eotdt, deathdt, eosdt, salast),
      trtdurd = dplyr::if_else(!is.na(trtsdt) & !is.na(trtedt),
                               as.integer(trtedt - trtsdt + 1L), NA_integer_),
      durcat  = dplyr::case_when(
        is.na(trtdurd)   ~ NA_character_,
        trtdurd <= 14    ~ "<=14 days",
        trtdurd <= 42    ~ ">14-<=42 days",
        TRUE             ~ ">42 days"
      ),
      durcat  = factor(durcat, levels = c("<=14 days", ">14-<=42 days", ">42 days")),
      ltfl    = dplyr::if_else(!is.na(trtdurd) & trtdurd > 43 & !is.na(dose02p),
                               "Y", "N", missing = "N"),
      switchfl = dplyr::if_else(!is.na(dose01p) & !is.na(dose02p) & dose01p != dose02p,
                                "Y", "N", missing = "N"),
      discfl  = dplyr::if_else(!is.na(eotrea) & !(eotrea %in% prog_terms) & eotrea != "Death",
                               "Y", "N", missing = "N"),
      # Planned-dose grouping factors
      trt01cd = factor(dose01p, levels = c(0, 25, 50, 100)),
      trt02cd = factor(dose02p, levels = c(0, 25, 50, 100))
    )

  adex %>%
    dplyr::select(
      studyid, domain, usubjid, subjid, cohort, cohortcd, arm, armcd,
      step, dose01p, dose01pu, dose02p, dose02pu, trt01p, trt02p,
      trt01cd, trt02cd,
      randdt, rfstdt, trtsdt, trtedt, deathdt, eotrea,
      trtdurd, durcat, ltfl, switchfl, discfl
    )
}

# Standalone execution -------------------------------------------------------
if (sys.nframe() == 0L) {
  suppressMessages({
    library(dplyr); library(tidyr); library(stringr); library(purrr)
    library(tibble); library(lubridate); library(yaml); library(targets)
  })
  source("R/external/functions.R")
  source("R/make_ad/helpers_ad.R")

  cfg <- yaml::read_yaml("config/cfg.yml")
  tar_load(dat); tar_load(adsl)
  adex <- make_adex(dat, adsl, cfg)
  message("adex: ", nrow(adex), " subjects")
  print(adex %>% dplyr::count(trt02cd))
  print(adex %>% dplyr::count(durcat))
}

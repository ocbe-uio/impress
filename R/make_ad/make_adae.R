# ADAE dataset creation — adverse events (SAP §7)
#
# Source files (active Viedoc export):
#   _AE.csv      — one row per adverse event (term, dates, severity, relationship,
#                  action taken, seriousness, outcome, fatal flag)
#   _MedDRA.csv  — MedDRA coding; the rows with itemname == "Adverse Event term"
#                  carry the SOC / PT coding for each AE form
#
# Join: AE <-> MedDRA on subjectid + subjectformseq (1:1 for AE-term rows).
#
# Mock-up note (sham data, EX export empty): the Safety Set is approximated by the
# Full Analysis Set and adverse events are grouped by *planned* dose. Assessment by
# treatment actually received is deferred to database lock / unblinding.
#
# Key ADAE variables:
#   STUDYID, DOMAIN ("ADAE"), USUBJID, SUBJID            ← identifiers (from adsl)
#   AGE, SEX, COHORT, COHORTCD, ARM, ARMCD               ← demographics / arm (adsl)
#   STEP, DOSE01P, DOSE02P, TRT01P, TRT02P, RANDDT, RFSTDT
#   TRTCD        ← stepped-wedge dose (mg) in effect at AE onset (trt_from_map)
#   AETERM       ← reported AE term
#   SOC_NAME, PT_NAME ← MedDRA system organ class / preferred term
#   PARAMCD/PARAM ← "AE" / "Adverse event"
#   ASTDT, AENDT ← AE start / end dates;  ASTDY ← study day of onset (rel. RFSTDT)
#   TRTEMFL      ← treatment-emergent flag ("Y" when ASTDT >= RFSTDT)
#   ASEV/ASEVN   ← severity (CTCAE grade, from AETOXGRI; ordered)
#   AREL/ARELFL  ← relationship to treatment / related flag
#   AESDFL       ← serious flag (AESER == "Yes")
#   ADTHFL       ← fatal flag (AESDTH == "Yes")
#   ADISCONFL    ← AE leading to treatment discontinuation (AEACN == "Drug withdrawn")
#   AEOUT        ← outcome

make_adae <- function(raw, adsl, cfg) {

  ae <- raw |> get_raw("ae")
  if (is.null(ae) || nrow(ae) == 0) {
    return(tibble::tibble())
  }

  meddra <- raw |> get_raw("meddra")

  # ---- MedDRA SOC/PT for AE-term rows ----------------------------------------
  ae_coding <- if (is.null(meddra)) {
    tibble::tibble(subjectid = character(), subjectformseq = double(),
                   soc_name = character(), pt_name = character())
  } else {
    meddra %>%
      dplyr::filter(stringr::str_squish(as.character(itemname)) == "Adverse Event term") %>%
      dplyr::transmute(
        subjectid,
        subjectformseq,
        soc_name = as.character(soc_name),
        pt_name  = as.character(pt_name)
      ) %>%
      dplyr::distinct(subjectid, subjectformseq, .keep_all = TRUE)
  }

  # ---- ADSL keys -------------------------------------------------------------
  adsl_keys <- adsl %>%
    dplyr::select(
      subjid, usubjid, studyid, age, sex,
      cohort, cohortcd, arm, armcd,
      step, dose01p, dose02p, trt01p, trt02p, randdt, rfstdt
    )

  sev_levels <- c("Mild", "Moderate", "Severe", "Life-threatening", "Death")
  related_vals <- c("Possibly related", "Probable related", "Definite")

  adae_base <- ae %>%
    dplyr::mutate(subjid = subjectid) %>%
    dplyr::left_join(ae_coding, by = c("subjectid", "subjectformseq")) %>%
    dplyr::left_join(adsl_keys, by = "subjid") %>%
    dplyr::mutate(
      studyid = cfg$studyid,
      domain  = "ADAE",
      usubjid = dplyr::coalesce(usubjid, paste0(cfg$studyid, "-", subjid)),
      astdt   = suppressWarnings(as.Date(aestdat)),
      aendt   = suppressWarnings(as.Date(aeendat)),
      aeterm  = as.character(aeterm),
      soc_name = dplyr::coalesce(soc_name, "Uncoded"),
      pt_name  = dplyr::coalesce(pt_name, aeterm),
      paramcd = "AE",
      param   = "Adverse event",
      # severity (CTCAE grade) — AESEV is empty in this export; use AETOXGRI
      asev    = factor(as.character(aetoxgri), levels = sev_levels, ordered = TRUE),
      asevn   = as.integer(asev),
      arel    = as.character(aerel),
      arelfl  = dplyr::if_else(arel %in% related_vals, "Y", "N", missing = "N"),
      aesdfl  = dplyr::if_else(toupper(as.character(aeser)) == "YES", "Y", "N", missing = "N"),
      adthfl  = dplyr::if_else(toupper(as.character(aesdth)) == "YES", "Y", "N", missing = "N"),
      adisconfl = dplyr::if_else(as.character(aeacn) == "Drug withdrawn", "Y", "N", missing = "N"),
      aeout   = as.character(aeout)
    )

  # ---- Study day of onset relative to first treatment ------------------------
  # admiral derives ASTDY from the ASTDT source variable (DT suffix -> DY)
  adae <- admiral::derive_vars_dy(
    dataset        = adae_base %>% dplyr::rename_with(toupper),
    reference_date = RFSTDT,
    source_vars    = admiral::exprs(ASTDT)
  ) %>%
    dplyr::rename_with(tolower)

  # ---- Treatment-emergent flag + stepped-wedge dose at onset -----------------
  trt_map <- trt_map_from_yaml(cfg)
  adae <- adae %>%
    dplyr::mutate(
      trtemfl = dplyr::if_else(!is.na(astdt) & !is.na(rfstdt) & astdt >= rfstdt,
                               "Y", "N", missing = "N"),
      trtcd   = trt_from_map(armcd, astdy, trt_map),
      trtcd   = factor(trtcd, levels = c(0, 25, 50, 100))
    )

  # ---- Final selection -------------------------------------------------------
  adae %>%
    dplyr::select(
      studyid, domain, usubjid, subjid,
      age, sex, cohort, cohortcd, arm, armcd,
      step, dose01p, dose02p, trt01p, trt02p, randdt, rfstdt,
      trtcd,
      paramcd, param,
      aeterm, soc_name, pt_name,
      astdt, aendt, astdy,
      trtemfl, asev, asevn, arel, arelfl, aesdfl, adthfl, adisconfl, aeout
    ) %>%
    dplyr::arrange(usubjid, astdt, soc_name, pt_name)
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
  adae <- make_adae(dat, adsl, cfg)
  message("adae: ", nrow(adae), " rows, ",
          dplyr::n_distinct(adae$usubjid), " subjects, ",
          sum(adae$trtemfl == "Y"), " TEAEs")
  print(adae %>% dplyr::count(soc_name, sort = TRUE))
}

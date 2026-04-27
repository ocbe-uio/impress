# ADSL dataset creation

make_adsl <- function(raw, cfg) {

  dm <- raw |> get_raw("dm")
  ran <- raw |> get_raw("ran")

  refdate <- raw %>%
    get_raw("event_dates") |>
    dplyr::filter(
      stringr::str_squish(eventname) == "Cycle 1",
      eventstatus == "Initiated"
    ) |>
    dplyr::arrange(eventinitiateddate) |>
    dplyr::group_by(subjectid) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      subjid = subjectid,
      rfstdt = lubridate::as_date(eventinitiateddate)
    )

  add_observed_levels <- function(values, base_levels = NULL) {
    base_levels <- if (is.null(base_levels)) character() else base_levels
    base_levels <- unique(as.character(base_levels))
    observed <- unique(as.character(values[!is.na(values)]))
    extra <- setdiff(observed, base_levels)
    unique(c(base_levels, extra))
  }

  adsl <- dm |>
    transmute(
      studyid = cfg$studyid,
      domain = "ADSL",
      usubjid = paste0(studyid, "-", subjectid),
      subjid = subjectid,
      cohort = cohort,
      cohortcd = cohortcd,
      country = cfg$country_code,
      age = dmage,
      ageu = "years",
      sex = sex
    ) %>%
    left_join(
      ran %>%
        transmute(
          subjid = subjectid,
          randdt = as_date(randat),
          step = as.character(ranst),
          dose01p = str_remove(randos, "mg") %>% as.numeric(),
          dose02p = str_remove(ranfudos, "mg") %>% as.numeric(),
          armcd_src = as.character(rantrt)
        ),
      by = "subjid"
    ) %>%
    left_join(
      refdate,
      by = "subjid"
    ) %>%
    mutate(
      rfstdt   = dplyr::coalesce(rfstdt, randdt),
      dose01pu = "mg",
      dose02pu = "mg",
      trt01p   = if_else(!is.na(dose01p), paste0("losartan ", dose01p, " mg; step ", step), NA_character_),
      trt02p   = if_else(!is.na(dose02p), paste0("losartan ", dose02p, " mg"), NA_character_),
      arm      = if_else(!is.na(dose01p) & !is.na(dose02p),
                        paste0("losartan ", dose01p, " mg; step ", step, " > losartan ", dose02p, " mg"),
                        NA_character_),
      armcd    = if_else(str_starts(armcd_src %||% "", "A"),
                        paste0(armcd_src, "_", dose02p),
                        armcd_src %||% NA_character_)
    )

  step_levels <- adsl %>%
    distinct(step) %>%
    mutate(step_num = suppressWarnings(as.numeric(step))) %>%
    arrange(step_num, step) %>%
    pull(step) %>%
    purrr::discard(is.na)

  trt01p_levels <- adsl %>%
    distinct(step, trt01p) %>%
    arrange(as.numeric(step), trt01p) %>%
    pull(trt01p) %>%
    purrr::discard(is.na)

  trt02p_levels <- adsl %>%
    distinct(dose02p, trt02p) %>%
    arrange(dose02p, trt02p) %>%
    pull(trt02p) %>%
    purrr::discard(is.na)

  arm_levels <- adsl %>%
    distinct(dose01p, dose02p, arm) %>%
    arrange(dose01p, dose02p, arm) %>%
    pull(arm) %>%
    purrr::discard(is.na)
  
  armcd_levels <- adsl %>%
    distinct(armcd) %>%
    arrange(armcd) %>%
    pull(armcd) %>%
    purrr::discard(is.na)

  adsl <- adsl %>%
    select(studyid, usubjid, subjid, cohort, cohortcd, country, age, sex,
          randdt, rfstdt, step, dose01p, dose01pu, trt01p,
          dose02p, dose02pu, trt02p, arm, armcd) %>%
    mutate(
      step   = factor(step,   levels = add_observed_levels(step,   step_levels)),
      trt01p = factor(trt01p, levels = add_observed_levels(trt01p, trt01p_levels)),
      trt02p = factor(trt02p, levels = add_observed_levels(trt02p, trt02p_levels)),
      arm    = factor(arm,    levels = add_observed_levels(arm,    arm_levels)),
      armcd  = factor(armcd,  levels = add_observed_levels(armcd,  armcd_levels))
    )

  return(adsl)
}

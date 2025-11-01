# ADSL dataset creation
make_adsl <- function(raw, cfg) {

  dm <- raw |> get_raw("dm")
  ran <- raw |> get_raw("ran")

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
        subjid   = subjectid,
        randdt   = as_date(randat),
        step     = as.character(ranst),
        dose01p  = str_remove(randos, 'mg') %>% as.numeric(),
        dose02p  = str_remove(ranfudos, 'mg') %>% as.numeric(),
        armcd_src = as.character(rantrt)
      ),
    by = "subjid"
  ) %>%
    mutate(
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
  adsl <- adsl %>%
    select(studyid, usubjid, subjid, cohort, cohortcd, country, age, sex,
          randdt, step, dose01p, dose01pu, trt01p,
          dose02p, dose02pu, trt02p, arm, armcd)

  return(adsl)
}

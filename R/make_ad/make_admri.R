# ============================================================
# 5) admri builder (BDS-style)
# ============================================================
make_admri <- function(raw, adsl, cfg) {

  visits <- visits_from_yaml(cfg)
  mri <- raw %>% get_raw("mri")
  mri_vars <- c(
    "mrt1cf", "mrt2cf", "mrt1ss", "mrt2ss", "mrt1ca", "mrt2ps", "mrpta",
    "mrrtt", "mrmef", "mrdsc", "mrrvsi", "mrrvc",  "mrrcva",
    "mrrdv", "mrdpta", "mrdrfa", "mrdrsi", "mrdiim", "mrdts", "mrdtv",
    "mrdtp", "mrrvt", "mrdtrv", "mrdte"
  )

  visit_levels_df <- visits$map %>%
    arrange(avisitn, avisit) %>%
    distinct(avisit, avisitn)
  avisit_levels <- visit_levels_df$avisit[!is.na(visit_levels_df$avisit)]
  avisitn_levels <- visit_levels_df$avisitn[!is.na(visit_levels_df$avisitn)]
  visit_schedule <- visit_levels_df %>%
    filter(!is.na(avisitn), !is.na(avisit)) %>%
    arrange(avisitn) %>%
    distinct(avisitn, .keep_all = TRUE)
  eventname_levels <- visits$map %>%
    arrange(avisitn, avisit, eventname) %>%
    distinct(eventname) %>%
    dplyr::pull(eventname)

  admri_base <-
    adsl %>%
    left_join(
      mri %>%
        transmute(
          subjid    = subjectid,
          eventname = eventname,
          adt       = as_date(mridat),
          avisit    = lookup_avisit(eventname, visits$map, fallback_title = TRUE),
          avisitn   = lookup_avisitn(eventname, visits$map, unmapped = visits$defaults$avisitn_unmapped),
          across(all_of(mri_vars))
        ),
      by = "subjid"
    ) |>
    mutate(
      randdt = as_date(randdt),
      ablfl = if_else(avisit == "Randomization And Baseline MRI", "Y", "N", missing = "N")
    ) 

  label_for <- function(var) {
    lab <- attr(admri_base[[var]], "label", exact = TRUE)
    if (is.null(lab) || !nzchar(lab)) var else as.character(lab)
  }

  paramcd_map <- purrr::set_names(purrr::map_chr(mri_vars, label_for), mri_vars)
  paramcd_levels <- mri_vars
  param_levels <- unname(paramcd_map[paramcd_levels])

  admri_long <-
    admri_base %>%
    tidyr::pivot_longer(
      cols = all_of(mri_vars),
      names_to = "paramcd",
      values_to = "aval"
    ) %>%
    mutate(
      param = unname(paramcd_map[paramcd])
    ) %>%
    rename_with(tolower, any_of(c("USUBJID", "RANDDT", "ADT", "AVISIT", "AVISITN")))

  admri <-
    derive_vars_dy(
      dataset = admri_long |> rename_with(toupper),
      reference_date = RFSTDT,
      source_vars = exprs(ADT)
    ) %>%
    arrange(USUBJID, ADT, PARAM) %>%
    rename_with(tolower) |>
    mutate(base_raw = if_else(ablfl == "Y", aval, NA_real_)) %>%
    group_by(usubjid, param) %>%
    mutate(
      base = if (any(!is.na(base_raw))) max(base_raw, na.rm = TRUE) else NA_real_,
      chg  = if_else(!is.na(aval) & !is.na(base), aval - base, NA_real_),
      pchg = if_else(!is.na(aval) & !is.na(base) & base != 0, 100 * (aval - base) / base, NA_real_)
    ) %>%
    ungroup() %>%
    select(-base_raw)

  # Map previously unmapped visits to the closest scheduled visit using ADY
  remap_needed <- (is.na(admri$avisitn) | admri$avisitn == visits$defaults$avisitn_unmapped) &
    !is.na(admri$ady)
  closest_visits <- nearest_visit_from_ady(admri$ady, visit_schedule)
  admri <- admri %>%
    mutate(
      avisitn = if_else(remap_needed, closest_visits$avisitn, avisitn),
      avisit  = if_else(remap_needed, closest_visits$avisit, avisit)
    )
  
  
  # Add treatment variable according to the randomisation and the dose levels.
  trt_map <- trt_map_from_yaml(cfg)
  admri <- admri %>%
    mutate(
      visit_day = avisitn,
      trt = trt_from_map(armcd, visit_day, trt_map)
    ) %>%
    select(-visit_day) |>
    mutate(
      trtcd = trt,
      trtcd = factor(trtcd, levels = c(0, 25, 50, 100)),
      trt = paste0(trt, " mg"),
      trt = factor(trt, levels = c("0 mg", "25 mg", "50 mg", "100 mg"))
    )
    



 admri <- admri %>%
   select(usubjid, eventname, avisit, avisitn, adt, ady, trt, trtcd, param, paramcd, aval, ablfl, randdt, everything()) |> 
   mutate(
     eventname = factor(eventname, levels = add_observed_levels(eventname, eventname_levels)),
     avisit  = factor(avisit,  levels = add_observed_levels(avisit,  avisit_levels)),
#     avisitn = factor(avisitn, levels = add_observed_levels(avisitn, avisitn_levels)),
     paramcd = factor(paramcd, levels = add_observed_levels(paramcd, paramcd_levels)),
     param   = factor(param,   levels = add_observed_levels(param,   param_levels)),
     ablfl   = factor(ablfl,   levels = add_observed_levels(ablfl,   c("N", "Y"))),
     studyid = cfg$studyid,
     domain  = "ADMRI"
   ) %>%
     select(studyid, domain, everything())
   
 
 


  return(admri)
}

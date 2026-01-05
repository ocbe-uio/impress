

visits_from_yaml <- function(ln) {
  v <- ln$visits
  if (is.null(v) || is.null(v$map)) stop("No 'visits: map:' section in cfg.yml")
  map <- tibble::tibble(
    eventname = stringr::str_squish(as.character(purrr::map_chr(v$map, "eventname"))),
    avisit    = stringr::str_squish(as.character(purrr::map_chr(v$map, "avisit", .null = NA_character_))),
    avisitn   = suppressWarnings(as.integer(purrr::map_int(v$map, "avisitn")))
  ) %>% distinct()
  defaults <- list(
    avisit_fallback_from_eventname = isTRUE(v$defaults$avisit_fallback_from_eventname),
    avisitn_unmapped               = as.integer(v$defaults$avisitn_unmapped %||% 999L)
  )
  list(map = map, defaults = defaults)
}


lookup_avisit  <- function(eventname, MAP, fallback_title = TRUE) {
  x <- tibble(eventname = str_squish(as.character(eventname))) %>%
    left_join(MAP, by = "eventname")
  out <- x$avisit
  if (fallback_title) out[is.na(out)] <- str_to_title(x$eventname[is.na(out)])
  out
}

lookup_avisitn <- function(eventname, MAP, unmapped = 999L) {
  x <- tibble(eventname = str_squish(as.character(eventname))) %>%
    left_join(MAP, by = "eventname")
  out <- x$avisitn
  out[is.na(out)] <- as.integer(unmapped)
  as.integer(out)
}



# ============================================================
# 5) adeff builder (BDS-style)
# ============================================================
make_adeff <- function(raw, adsl, cfg) {

  visits <- visits_from_yaml(cfg)
  mri <- raw %>% get_raw("mri")
  mri_vars <- c(
    "mrt1cf", "mrt2cf", "mrt1ss", "mrt2ss", "mrt1ca", "mrt2ps", "mrpta",
    "mrrtt", "mrmef", "mrdsc", "mrrvsi", "mrrvc", "mrrva", "mrrcva",
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

  add_observed_levels <- function(values, base_levels) {
    base_levels <- unique(as.character(base_levels))
    observed <- unique(as.character(values[!is.na(values)]))
    extra <- setdiff(observed, base_levels)
    unique(c(base_levels, extra))
  }

  nearest_visit_from_ady <- function(days) {
    if (!nrow(visit_schedule)) {
      return(list(
        avisitn = rep(NA_integer_, length(days)),
        avisit  = rep(NA_character_, length(days))
      ))
    }
    idx <- vapply(days, function(d) {
      if (is.na(d)) return(NA_integer_)
      which.min(abs(d - visit_schedule$avisitn))
    }, integer(1))
    list(
      avisitn = ifelse(is.na(idx), NA_integer_, visit_schedule$avisitn[idx]),
      avisit  = ifelse(is.na(idx), NA_character_, visit_schedule$avisit[idx])
    )
  }

  adeff_base <-
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
    lab <- attr(adeff_base[[var]], "label", exact = TRUE)
    if (is.null(lab) || !nzchar(lab)) var else as.character(lab)
  }

  paramcd_map <- purrr::set_names(purrr::map_chr(mri_vars, label_for), mri_vars)
  paramcd_levels <- mri_vars
  param_levels <- unname(paramcd_map[paramcd_levels])

  adeff_long <-
    adeff_base %>%
    tidyr::pivot_longer(
      cols = all_of(mri_vars),
      names_to = "paramcd",
      values_to = "aval"
    ) %>%
    mutate(
      param = unname(paramcd_map[paramcd])
    ) %>%
    rename_with(tolower, any_of(c("USUBJID", "RANDDT", "ADT", "AVISIT", "AVISITN")))

  adeff <-
    derive_vars_dy(
      dataset = adeff_long |> rename_with(toupper),
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
  remap_needed <- (is.na(adeff$avisitn) | adeff$avisitn == visits$defaults$avisitn_unmapped) &
    !is.na(adeff$ady)
  closest_visits <- nearest_visit_from_ady(adeff$ady)
  adeff <- adeff %>%
    mutate(
      avisitn = if_else(remap_needed, closest_visits$avisitn, avisitn),
      avisit  = if_else(remap_needed, closest_visits$avisit, avisit)
    )
  
  
  # Add treatment variable according to the randomisation and the dose levels.
  ar_trt <- function(day, start_day, dose, post_cycle4_dose) {
    case_when(
      is.na(day) ~ NA_real_,
      day <= start_day ~ 0,
      day > start_day & day <= 43 ~ dose,
      day > 43 & day <= 169 ~ post_cycle4_dose,
      day > 169 ~ 0
    )
  }

   an_trt <- function(day, start_day, dose, post_cycle4_dose) {
    case_when(
      is.na(day) ~ NA_real_,
      day <= start_day ~ 0,
      day > start_day & day <= 43 ~ dose,
      day > 43 & day <= 239 ~ post_cycle4_dose,
      day > 239 ~ 0
    )
  }

  bm_trt <- function(day, start_day, dose) {
    case_when(
      is.na(day) ~ NA_real_,
      day < start_day ~ 0,
      day >= start_day & day <= 1270 ~ dose,
      day > 270 ~ 0
    )
  }

  adeff <- adeff %>%
    mutate(
      visit_day = avisitn,
      armcd_trt = armcd,
      trt = case_when(
        armcd_trt == "AR1_1_0" ~ ar_trt(visit_day, 1, 25, 0),
        armcd_trt == "AR1_1_25" ~ ar_trt(visit_day, 1, 25, 25),
        armcd_trt == "AR1_2_0" ~ ar_trt(visit_day, 15, 25, 0),
        armcd_trt == "AR1_2_25" ~ ar_trt(visit_day, 15, 25, 25),
        armcd_trt == "AR1_3_0" ~ ar_trt(visit_day, 29, 25, 0),
        armcd_trt == "AR1_3_25" ~ ar_trt(visit_day, 29, 25, 25),
        armcd_trt == "AR2_4_0" ~ ar_trt(visit_day, 1, 50, 0),
        armcd_trt == "AR2_4_50" ~ ar_trt(visit_day, 1, 50, 50),
        armcd_trt == "AR2_5_0" ~ ar_trt(visit_day, 15, 50, 0),
        armcd_trt == "AR2_5_50" ~ ar_trt(visit_day, 15, 50, 50),
        armcd_trt == "AR2_6_0" ~ ar_trt(visit_day, 29, 50, 0),
        armcd_trt == "AR2_6_50" ~ ar_trt(visit_day, 29, 50, 50),
        armcd_trt == "AR3_1_100" ~ ar_trt(visit_day, 1, 100, 100),
        armcd_trt == "AR3_2_100" ~ ar_trt(visit_day, 15, 100, 100),
        armcd_trt == "AR3_3_100" ~ ar_trt(visit_day, 29, 100, 100),
        armcd_trt == "AN1_1_0" ~ an_trt(visit_day, 1, 25, 0),
        armcd_trt == "AN1_1_25" ~ an_trt(visit_day, 1, 25, 25),
        armcd_trt == "AN1_2_0" ~ an_trt(visit_day, 15, 25, 0),
        armcd_trt == "AN1_2_25" ~ an_trt(visit_day, 15, 25, 25),
        armcd_trt == "AN1_3_0" ~ an_trt(visit_day, 29, 25, 0),
        armcd_trt == "AN1_3_25" ~ an_trt(visit_day, 29, 25, 25),
        armcd_trt == "AN2_4_0" ~ an_trt(visit_day, 1, 50, 0),
        armcd_trt == "AN2_4_50" ~ an_trt(visit_day, 1, 50, 50),
        armcd_trt == "AN2_5_0" ~ an_trt(visit_day, 15, 50, 0),
        armcd_trt == "AN2_5_50" ~ an_trt(visit_day, 15, 50, 50),
        armcd_trt == "AN2_6_0" ~ an_trt(visit_day, 29, 50, 0),
        armcd_trt == "AN2_6_50" ~ an_trt(visit_day, 29, 50, 50),
        armcd_trt == "AN3_1_100" ~ an_trt(visit_day, 1, 100, 100),
        armcd_trt == "AN3_2_100" ~ an_trt(visit_day, 15, 100, 100),
        armcd_trt == "AN3_3_100" ~ an_trt(visit_day, 29, 100, 100),
        armcd_trt == "BM1" ~ bm_trt(visit_day, 1, 50),
        armcd_trt == "BM2" ~ bm_trt(visit_day, 1091, 50),
        armcd_trt == "BM3" ~ bm_trt(visit_day, 1181, 50),
        TRUE ~ NA_real_
      )
    ) %>%
    select(-visit_day, -armcd_trt) |>
    mutate(
      trtcd = trt,
      trtcd = factor(trtcd, levels = c(0, 25, 50, 100)),
      trt = paste0(trt, " mg"),
      trt = factor(trt, levels = c("0 mg", "25 mg", "50 mg", "100 mg"))
    )
    



 adeff <- adeff %>%
   select(usubjid, eventname, avisit, avisitn, adt, ady, trt, trtcd, param, paramcd, aval, ablfl, randdt, everything()) |> 
   mutate(
     eventname = factor(eventname, levels = add_observed_levels(eventname, eventname_levels)),
     avisit  = factor(avisit,  levels = add_observed_levels(avisit,  avisit_levels)),
#     avisitn = factor(avisitn, levels = add_observed_levels(avisitn, avisitn_levels)),
     paramcd = factor(paramcd, levels = add_observed_levels(paramcd, paramcd_levels)),
     param   = factor(param,   levels = add_observed_levels(param,   param_levels)),
     ablfl   = factor(ablfl,   levels = add_observed_levels(ablfl,   c("N", "Y"))),
     studyid = cfg$studyid,
     domain  = "ADEFF"
   ) %>%
     select(studyid, domain, everything())
   
 
 


  return(adeff)
}

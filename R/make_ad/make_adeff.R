



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
      reference_date = RANDDT,
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
  
  


 adeff <- adeff %>%
   select(usubjid, eventname, avisit, avisitn, adt, ady, param, paramcd, aval, ablfl, randdt, everything()) |> 
   mutate(
     studyid = cfg$studyid,
     domain  = "ADEFF"
   ) %>%
     select(studyid, domain, everything())
   
 
 


  return(adeff)
}

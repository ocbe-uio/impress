# Return the appropriate raw data object based on cfg$mode.
# When mode == "sham", shamraw carries sham-randomised treatment assignments;
# otherwise raw and shamraw are identical and either may be used.
# Always call this rather than choosing raw/shamraw manually.
effective_raw <- function(raw, shamraw, cfg) {
  if (isTRUE(cfg$mode == "sham")) shamraw else raw
}

visits_from_yaml <- function(ln) {
  v <- ln$visits
  if (is.null(v) || is.null(v$map)) stop("No 'visits: map:' section in cfg.yml")
  map <- tibble::tibble(
    eventname = stringr::str_squish(as.character(purrr::map_chr(v$map, "eventname"))),
    avisit    = stringr::str_squish(as.character(purrr::map_chr(v$map, "avisit", .null = NA_character_))),
    avisitn   = suppressWarnings(as.integer(purrr::map_int(v$map, "avisitn")))
  ) %>% dplyr::distinct()
  defaults <- list(
    avisit_fallback_from_eventname = isTRUE(v$defaults$avisit_fallback_from_eventname),
    avisitn_unmapped               = as.integer(v$defaults$avisitn_unmapped %||% 999L)
  )
  list(map = map, defaults = defaults)
}

lookup_avisit <- function(eventname, map, fallback_title = TRUE) {
  x <- tibble::tibble(eventname = stringr::str_squish(as.character(eventname))) %>%
    dplyr::left_join(map, by = "eventname")
  out <- x$avisit
  if (fallback_title) out[is.na(out)] <- stringr::str_to_title(x$eventname[is.na(out)])
  out
}

lookup_avisitn <- function(eventname, map, unmapped = 999L) {
  x <- tibble::tibble(eventname = stringr::str_squish(as.character(eventname))) %>%
    dplyr::left_join(map, by = "eventname")
  out <- x$avisitn
  out[is.na(out)] <- as.integer(unmapped)
  as.integer(out)
}

add_observed_levels <- function(values, base_levels) {
  base_levels <- unique(as.character(base_levels))
  observed <- unique(as.character(values[!is.na(values)]))
  extra <- setdiff(observed, base_levels)
  unique(c(base_levels, extra))
}

nearest_visit_from_ady <- function(days, visit_schedule) {
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

ar_trt <- function(day, start_day, dose, post_cycle4_dose) {
  dplyr::case_when(
    is.na(day) ~ NA_real_,
    day <= start_day ~ 0,
    day > start_day & day <= 43 ~ dose,
    day > 43 & day <= 169 ~ post_cycle4_dose,
    day > 169 ~ 0
  )
}

an_trt <- function(day, start_day, dose, post_cycle4_dose) {
  dplyr::case_when(
    is.na(day) ~ NA_real_,
    day <= start_day ~ 0,
    day > start_day & day <= 43 ~ dose,
    day > 43 & day <= 239 ~ post_cycle4_dose,
    day > 239 ~ 0
  )
}

bm_trt <- function(day, start_day, dose) {
  dplyr::case_when(
    is.na(day) ~ NA_real_,
    day < start_day ~ 0,
    day >= start_day & day <= 1270 ~ dose,
    day > 270 ~ 0
  )
}

trt_map_from_yaml <- function(cfg) {
  trt_map <- cfg$trt_map
  if (is.null(trt_map)) stop("No 'trt_map' section in cfg.yml")
  tibble::tibble(
    armcd = purrr::map_chr(trt_map, "armcd"),
    diagnosis = stringr::str_to_upper(purrr::map_chr(trt_map, "diagnosis")),
    start_day = as.integer(purrr::map_int(trt_map, "start_day")),
    dose = as.numeric(purrr::map_dbl(trt_map, "dose")),
    post_cycle4_dose = as.numeric(purrr::map_dbl(trt_map, "post_cycle4_dose", .null = NA_real_))
  ) %>%
    dplyr::distinct()
}

trt_from_map <- function(armcd, day, trt_map) {
  idx <- match(armcd, trt_map$armcd)
  diagnosis <- trt_map$diagnosis[idx]
  start_day <- trt_map$start_day[idx]
  dose <- trt_map$dose[idx]
  post_cycle4_dose <- trt_map$post_cycle4_dose[idx]
  dplyr::case_when(
    diagnosis == "AR" ~ ar_trt(day, start_day, dose, post_cycle4_dose),
    diagnosis == "AN" ~ an_trt(day, start_day, dose, post_cycle4_dose),
    diagnosis == "BM" ~ bm_trt(day, start_day, dose),
    TRUE ~ NA_real_
  )
}

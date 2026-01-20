#' Build ADEFF for neurologic status, steroids, and QoL endpoints
#'
#' Constructs a BDS-style ADaM dataset for continuous neurologic performance
#' scores (KPS, ECOG, NANO), QoL measures (QLQ-C30, QLQ-BN20), and steroid use.
#' Visits are mapped from `cfg.yml`, baseline is restricted to observations on or
#' before randomization, and stepped-wedge treatment assignment is derived from
#' `armcd` and visit day.
#'
#' @param raw Raw data list from `make_raw()`.
#' @param adsl Subject-level dataset containing randomization and arm metadata.
#' @param cfg Configuration list (including visit mapping).
#'
#' @return A tibble with ADEFF records including derived baseline and change
#'   variables, and time-varying treatment indicators.
#'
#' @details Manual verification: Inge Christoffer Olsen checked patient 2054,
#'   Cycle 4 against the eCRF on 2025-03-08.
#'
#' @export
# ============================================================
# 6) adeff builder (neurologic status, steroids, QoL; BDS-style)
# ============================================================
make_adeff <- function(raw, adsl, cfg) {

  visits <- visits_from_yaml(cfg)

  visit_levels_df <- visits$map %>%
    arrange(avisitn, avisit) %>%
    distinct(avisit, avisitn)
  avisit_levels <- visit_levels_df$avisit[!is.na(visit_levels_df$avisit)]
  visit_schedule <- visit_levels_df %>%
    filter(!is.na(avisitn), !is.na(avisit)) %>%
    arrange(avisitn) %>%
    distinct(avisitn, .keep_all = TRUE)
  eventname_levels <- visits$map %>%
    arrange(avisitn, avisit, eventname) %>%
    distinct(eventname) %>%
    dplyr::pull(eventname)

  safe_label <- function(data, var, default) {
    if (is.null(data)) {
      default
    } else {
      lab <- attr(data[[var]], "label", exact = TRUE)
      if (is.null(lab) || !nzchar(lab)) default else as.character(lab)
    }
  }

  build_long <- function(df, value_col, date_col, paramcd, param, avalc = NA_character_) {
    if (is.null(df)) {
      return(tibble::tibble())
    }
    df %>%
      mutate(
        subjid = subjectid,
        adt = suppressWarnings(as.Date(.data[[date_col]])),
        adt = dplyr::coalesce(adt, as.Date(eventdate))
      ) %>%
      filter(!is.na(.data[[value_col]])) %>%
      transmute(
        subjid,
        eventname,
        adt,
        paramcd = paramcd,
        param = param,
        aval = as.numeric(.data[[value_col]]),
        avalc = avalc
      )
  }

  build_long_cat <- function(df, value_col, date_col, paramcd, param) {
    if (is.null(df)) {
      return(tibble::tibble())
    }
    df %>%
      mutate(
        subjid = subjectid,
        adt = suppressWarnings(as.Date(.data[[date_col]])),
        adt = dplyr::coalesce(adt, as.Date(eventdate))
      ) %>%
      filter(!is.na(.data[[value_col]])) %>%
      transmute(
        subjid,
        eventname,
        adt,
        paramcd = paramcd,
        param = param,
        aval = NA_real_,
        avalc = as.character(.data[[value_col]])
      )
  }

  ecog <- raw |> get_raw("ecog")
  kps <- raw |> get_raw("kps")
  nano <- raw |> get_raw("nano")
  cs <- raw |> get_raw("cs")
  tr <- raw |> get_raw("tr")
  qlq <- raw |> get_raw("qlq")
  qlqbn <- raw |> get_raw("qlqbn")

  paramcd_map <- c(
    ecog = sub(" - Code$", "", safe_label(ecog, "ecogscd", "ECOG Score")),
    kps = sub(" - Code$", "", safe_label(kps, "kpsscd", "Karnofsky Performance Scale")),
    nano_tot = safe_label(nano, "nanotot", "Total NANO score"),
    cs_steroids_mri = safe_label(cs, "csmriyn", "Steroids at time of MRI"),
    qlq_c30 = safe_label(qlq, "qlqc3sc", "EORTC QLQ-C30 score"),
    qlq_bn20 = safe_label(qlqbn, "qlqbnsc", "EORTC QLQ-BN20 score"),
    trpspro = sub(" - Code$", "", safe_label(tr, "trpsprocd", "Occurrence of pseudoprogression")),
    trsdisea = sub(" - Code$", "", safe_label(tr, "trsdiseacd", "State of disease according to RANO")),
    trorresp = sub(" - Code$", "", safe_label(tr, "trorrespcd", "Best overall radiographic response"))
  )

  # Build a unified long-format dataset across neuro/QoL/steroid sources.
  cs_long <- if (!is.null(cs)) {
    cs %>%
      mutate(
        subjid = subjectid,
        adt = suppressWarnings(as.Date(eventdate)),
        steroid_flag = case_when(
          toupper(as.character(csmriyn)) == "YES" ~ 1,
          toupper(as.character(csmriyn)) == "NO" ~ 0,
          TRUE ~ NA_real_
        )
      ) %>%
      filter(!is.na(steroid_flag)) %>%
      transmute(
        subjid,
        eventname,
        adt,
        paramcd = "cs_steroids_mri",
        param = paramcd_map[["cs_steroids_mri"]],
        aval = NA_real_,
        avalc = as.character(csmriyn)
      )
  } else {
    tibble::tibble()
  }

  tr_long <- if (!is.null(tr)) {
    dplyr::bind_rows(
      tr %>%
        mutate(subjid = subjectid, adt = suppressWarnings(as.Date(eventdate))) %>%
        filter(!is.na(trpsprocd)) %>%
        transmute(
          subjid,
          eventname,
          adt,
          paramcd = "trpspro",
          param = paramcd_map[["trpspro"]],
          aval = NA_real_,
          avalc = as.character(trpspro)
        ),
      tr %>%
        mutate(subjid = subjectid, adt = suppressWarnings(as.Date(eventdate))) %>%
        filter(!is.na(trsdiseacd)) %>%
        transmute(
          subjid,
          eventname,
          adt,
          paramcd = "trsdisea",
          param = paramcd_map[["trsdisea"]],
          aval = NA_real_,
          avalc = as.character(trsdisea)
        ),
      tr %>%
        mutate(subjid = subjectid, adt = suppressWarnings(as.Date(eventdate))) %>%
        filter(!is.na(trorrespcd)) %>%
        transmute(
          subjid,
          eventname,
          adt,
          paramcd = "trorresp",
          param = paramcd_map[["trorresp"]],
          aval = NA_real_,
          avalc = as.character(trorresp)
        )
    )
  } else {
    tibble::tibble()
  }

  adeff_long <- dplyr::bind_rows(
    build_long(ecog, "ecogscd", "eventdate", "ecog", paramcd_map[["ecog"]]),
    build_long(kps, "kpsscd", "eventdate", "kps", paramcd_map[["kps"]]),
    build_long(nano, "nanotot", "eventdate", "nano_tot", paramcd_map[["nano_tot"]]),
    cs_long,
    build_long(qlq, "qlqc3sc", "qlqc3dat", "qlq_c30", paramcd_map[["qlq_c30"]]),
    build_long(qlqbn, "qlqbnsc", "qlqbndat", "qlq_bn20", paramcd_map[["qlq_bn20"]]),
    tr_long
  )

  if (!nrow(adeff_long)) {
    return(tibble::tibble())
  }

  adeff_base <- adeff_long %>%
    left_join(adsl, by = "subjid") %>%
    mutate(
      avisit = lookup_avisit(eventname, visits$map, fallback_title = TRUE),
      avisitn = lookup_avisitn(eventname, visits$map, unmapped = visits$defaults$avisitn_unmapped)
    )

  adeff <-
    derive_vars_dy(
      dataset = adeff_base |> rename_with(toupper),
      reference_date = RFSTDT,
      source_vars = exprs(ADT)
    ) %>%
    rename_with(tolower) %>%
    arrange(usubjid, adt, paramcd)

  # Baseline is restricted to observations on/before randomization.
  adeff <- adeff %>%
    group_by(usubjid, paramcd) %>%
    mutate(
      baseline_date = {
        eligible <- !is.na(aval) & !is.na(randdt) & !is.na(adt) & adt <= randdt
        if (any(eligible)) max(adt[eligible], na.rm = TRUE) else as.Date(NA)
      },
      ablfl = if_else(!is.na(adt) & !is.na(baseline_date) & adt == baseline_date, "Y", "N", missing = "N"),
      base_raw = if_else(ablfl == "Y", aval, NA_real_),
      base = if (any(!is.na(base_raw))) max(base_raw, na.rm = TRUE) else NA_real_,
      chg = if_else(!is.na(aval) & !is.na(base), aval - base, NA_real_),
      pchg = if_else(!is.na(aval) & !is.na(base) & base != 0, 100 * (aval - base) / base, NA_real_)
    ) %>%
    ungroup() %>%
    select(-baseline_date, -base_raw)

  remap_needed <- (is.na(adeff$avisitn) | adeff$avisitn == visits$defaults$avisitn_unmapped) &
    !is.na(adeff$ady)
  closest_visits <- nearest_visit_from_ady(adeff$ady, visit_schedule)
  adeff <- adeff %>%
    mutate(
      avisitn = if_else(remap_needed, closest_visits$avisitn, avisitn),
      avisit = if_else(remap_needed, closest_visits$avisit, avisit)
    )

  # Apply stepped-wedge treatment assignment based on armcd and visit day.
  trt_map <- trt_map_from_yaml(cfg)
  adeff <- adeff %>%
    mutate(
      visit_day = avisitn,
      trt = trt_from_map(armcd, visit_day, trt_map)
    ) %>%
    select(-visit_day) %>%
    mutate(
      trtcd = trt,
      trtcd = factor(trtcd, levels = c(0, 25, 50, 100)),
      trt = paste0(trt, " mg"),
      trt = factor(trt, levels = c("0 mg", "25 mg", "50 mg", "100 mg"))
    )

  paramcd_levels <- names(paramcd_map)
  param_levels <- unname(paramcd_map[paramcd_levels])

  # Manual verification note: Inge Christoffer Olsen checked patient 2054, Cycle 4 against the eCRF on 2025-03-08.
  adeff <- adeff %>%
    mutate(
      eventname = factor(eventname, levels = add_observed_levels(eventname, eventname_levels)),
      avisit = factor(avisit, levels = add_observed_levels(avisit, avisit_levels)),
      paramcd = factor(paramcd, levels = add_observed_levels(paramcd, paramcd_levels)),
      param = factor(param, levels = add_observed_levels(param, param_levels)),
      ablfl = factor(ablfl, levels = add_observed_levels(ablfl, c("N", "Y"))),
      studyid = cfg$studyid,
      domain = "ADEFF"
    ) %>%
    select(studyid, domain, everything())

  return(adeff)
}

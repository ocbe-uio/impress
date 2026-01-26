make_adtte <- function(raw, adsl, cfg) {
  dp <- raw |> get_raw("dp")
  sa <- raw |> get_raw("sa")

  prog_dates <- if (is.null(dp)) {
    tibble::tibble(subjid = character(), progdt = as.Date(character()))
  } else {
    dp %>%
      mutate(
        subjid = subjectid,
        dpdat = suppressWarnings(as.Date(dpdat)),
        eventdate = suppressWarnings(as.Date(eventdate)),
        progdt = dplyr::coalesce(dpdat, eventdate)
      ) %>%
      filter(!is.na(progdt)) %>%
      group_by(subjid) %>%
      summarise(progdt = min(progdt, na.rm = TRUE), .groups = "drop")
  }

  surv_dates <- if (is.null(sa)) {
    tibble::tibble(subjid = character(), deathdt = as.Date(character()), lastdt = as.Date(character()))
  } else {
    sa %>%
      mutate(
        subjid = subjectid,
        deathdt = suppressWarnings(as.Date(sadtdat)),
        lastdt = suppressWarnings(as.Date(saendat))
      ) %>%
      group_by(subjid) %>%
      summarise(
        deathdt = if (any(!is.na(deathdt))) min(deathdt, na.rm = TRUE) else as.Date(NA),
        lastdt = if (any(!is.na(lastdt))) max(lastdt, na.rm = TRUE) else as.Date(NA),
        .groups = "drop"
      )
  }

  base <- adsl %>%
    mutate(
      startdt = dplyr::coalesce(randdt, rfstdt)
    ) %>%
    select(
      studyid, usubjid, subjid, cohort, cohortcd, randdt, rfstdt, startdt,
      arm, armcd, trt01p, trt02p, dose01p, dose02p
    ) %>%
    left_join(prog_dates, by = "subjid") %>%
    left_join(surv_dates, by = "subjid")

  make_tte <- function(df, paramcd, param, eventdt, censor_dt, cutoff_days = NULL,
                       evntdesc = NULL, srcdom = NULL, srcvar = NULL) {
    startdt <- df$startdt
    cutoff_dt <- if (is.null(cutoff_days)) {
      rep(as.Date(NA), length(startdt))
    } else {
      startdt + as.integer(cutoff_days)
    }

    event_dt <- eventdt
    if (!is.null(cutoff_days)) {
      event_dt <- ifelse(!is.na(event_dt) & event_dt <= cutoff_dt, event_dt, as.Date(NA))
      event_dt <- as.Date(event_dt, origin = "1970-01-01")
    }

    censor_dt_use <- censor_dt
    if (!is.null(cutoff_days)) {
      censor_dt_use <- pmin(censor_dt_use, cutoff_dt, na.rm = TRUE)
    }

    adt <- ifelse(!is.na(event_dt), event_dt, censor_dt_use)
    adt <- as.Date(adt, origin = "1970-01-01")
    cnsr <- ifelse(!is.na(event_dt), 0L, 1L)
    aval <- ifelse(!is.na(adt) & !is.na(startdt), as.integer(adt - startdt + 1L), NA_integer_)

    df %>%
      mutate(
        paramcd = paramcd,
        param = param,
        astdt = startdt,
        astdy = ifelse(!is.na(startdt) & !is.na(randdt), as.integer(startdt - randdt + 1L), NA_integer_),
        adt = adt,
        ady = ifelse(!is.na(adt) & !is.na(randdt), as.integer(adt - randdt + 1L), NA_integer_),
        aval = aval,
        avalu = "DAYS",
        cnsr = cnsr,
        evntdesc = evntdesc,
        evntstat = ifelse(cnsr == 0L, "EVENT", "CENSORED"),
        srcdom = srcdom,
        srcvar = srcvar,
        cutoff_days = if (is.null(cutoff_days)) NA_integer_ else as.integer(cutoff_days)
      )
  }

  pfs_event <- pmin(base$progdt, base$deathdt, na.rm = TRUE)
  pfs_event[is.infinite(pfs_event)] <- as.Date(NA)
  pfs_srcdom <- ifelse(!is.na(base$progdt) & (is.na(base$deathdt) | base$progdt <= base$deathdt), "DP", "SA")
  pfs_srcvar <- ifelse(pfs_srcdom == "DP", "DPDAT", "SADTDAT")

  pfs <- make_tte(
    base, "PFS", "Progression-free survival", pfs_event, base$lastdt,
    evntdesc = "Progression or death", srcdom = pfs_srcdom, srcvar = pfs_srcvar
  )
  os <- make_tte(
    base, "OS", "Overall survival", base$deathdt, base$lastdt,
    evntdesc = "Death", srcdom = ifelse(!is.na(base$deathdt), "SA", NA_character_), srcvar = "SADTDAT"
  )

  pfs6 <- make_tte(
    base, "PFS6M", "6-month progression-free survival", pfs_event, base$lastdt, cutoff_days = 180L,
    evntdesc = "Progression or death by day 180", srcdom = pfs_srcdom, srcvar = pfs_srcvar
  )
  os12 <- make_tte(
    base, "OS12M", "Overall survival at 12 months", base$deathdt, base$lastdt, cutoff_days = 365L,
    evntdesc = "Death by day 365", srcdom = ifelse(!is.na(base$deathdt), "SA", NA_character_), srcvar = "SADTDAT"
  )
  os24 <- make_tte(
    base, "OS24M", "Overall survival at 24 months", base$deathdt, base$lastdt, cutoff_days = 730L,
    evntdesc = "Death by day 730", srcdom = ifelse(!is.na(base$deathdt), "SA", NA_character_), srcvar = "SADTDAT"
  )

  adtte <- bind_rows(pfs, pfs6, os, os12, os24) %>%
    mutate(
      domain = "ADTTE"
    ) %>%
    select(studyid, domain, usubjid, subjid, cohort, cohortcd, arm, armcd,
           trt01p, trt02p, dose01p, dose02p, astdt, astdy, adt, ady,
           aval, avalu, cnsr, evntdesc, evntstat, srcdom, srcvar,
           paramcd, param, cutoff_days, randdt, rfstdt, progdt, deathdt, lastdt)

  adtte
}

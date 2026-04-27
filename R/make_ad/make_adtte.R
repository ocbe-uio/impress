make_adtte <- function(raw, adsl, cfg) {
  dp  <- raw |> get_raw("dp")
  sa  <- raw |> get_raw("sa")
  eos <- raw |> get_raw("eos")

  # ---- DP form: earliest progression date ----------------------------------
  prog_dates <- if (is.null(dp)) {
    tibble::tibble(subjid = character(), progdt_dp = as.Date(character()))
  } else {
    dp %>%
      mutate(
        subjid    = subjectid,
        dpdat     = suppressWarnings(as.Date(dpdat)),
        eventdate = suppressWarnings(as.Date(eventdate)),
        progdt_dp = dplyr::coalesce(dpdat, eventdate)
      ) %>%
      filter(!is.na(progdt_dp)) %>%
      group_by(subjid) %>%
      summarise(progdt_dp = min(progdt_dp, na.rm = TRUE), .groups = "drop")
  }

  # ---- SA form: death date and last known alive date -----------------------
  surv_dates <- if (is.null(sa)) {
    tibble::tibble(
      subjid    = character(),
      deathdt_sa = as.Date(character()),
      lastdt_sa  = as.Date(character())
    )
  } else {
    sa %>%
      mutate(
        subjid     = subjectid,
        deathdt_sa = suppressWarnings(as.Date(sadtdat)),
        lastdt_sa  = suppressWarnings(as.Date(saendat))
      ) %>%
      group_by(subjid) %>%
      summarise(
        deathdt_sa = if (any(!is.na(deathdt_sa))) min(deathdt_sa, na.rm = TRUE) else as.Date(NA),
        lastdt_sa  = if (any(!is.na(lastdt_sa)))  max(lastdt_sa,  na.rm = TRUE) else as.Date(NA),
        .groups = "drop"
      )
  }

  # ---- EOS form: death date, progression date, last contact, cause of death
  eos_dates <- if (is.null(eos)) {
    tibble::tibble(
      subjid      = character(),
      deathdt_eos = as.Date(character()),
      progdt_eos  = as.Date(character()),
      lastdt_eos  = as.Date(character()),
      dthcaus     = character()
    )
  } else {
    eos %>%
      mutate(
        subjid      = subjectid,
        deathdt_eos = suppressWarnings(as.Date(eosdtdat)),
        progdt_eos  = suppressWarnings(as.Date(eospddat)),
        lastdt_eos  = suppressWarnings(as.Date(eosdat)),
        dthcaus     = as.character(eosdth)
      ) %>%
      filter(!is.na(subjid)) %>%
      group_by(subjid) %>%
      summarise(
        deathdt_eos = if (any(!is.na(deathdt_eos))) min(deathdt_eos, na.rm = TRUE) else as.Date(NA),
        progdt_eos  = if (any(!is.na(progdt_eos)))  min(progdt_eos,  na.rm = TRUE) else as.Date(NA),
        lastdt_eos  = if (any(!is.na(lastdt_eos)))  max(lastdt_eos,  na.rm = TRUE) else as.Date(NA),
        dthcaus     = dplyr::first(dthcaus[!is.na(dthcaus) & dthcaus != "NA"]),
        .groups     = "drop"
      )
  }

  # ---- Combine all sources -------------------------------------------------
  base <- adsl %>%
    mutate(
      startdt = dplyr::coalesce(rfstdt, randdt)
    ) %>%
    select(
      studyid, usubjid, subjid, cohort, cohortcd, randdt, rfstdt, startdt,
      arm, armcd, trt01p, trt02p, dose01p, dose02p
    ) %>%
    left_join(prog_dates, by = "subjid") %>%
    left_join(surv_dates, by = "subjid") %>%
    left_join(eos_dates,  by = "subjid") %>%
    mutate(
      # Earliest progression date across DP and EOS; track which source won
      progdt = {
        x <- pmin(progdt_dp, progdt_eos, na.rm = TRUE)
        x[is.infinite(x)] <- NA
        as.Date(x, origin = "1970-01-01")
      },
      progdt_src = dplyr::case_when(
        !is.na(progdt_dp) & (is.na(progdt_eos) | progdt_dp <= progdt_eos) ~ "DP",
        !is.na(progdt_eos)                                                 ~ "EOS",
        TRUE                                                               ~ NA_character_
      ),
      # Earliest death date across SA and EOS; track which source won
      deathdt = {
        x <- pmin(deathdt_sa, deathdt_eos, na.rm = TRUE)
        x[is.infinite(x)] <- NA
        as.Date(x, origin = "1970-01-01")
      },
      deathdt_src = dplyr::case_when(
        !is.na(deathdt_sa) & (is.na(deathdt_eos) | deathdt_sa <= deathdt_eos) ~ "SA",
        !is.na(deathdt_eos)                                                    ~ "EOS",
        TRUE                                                                   ~ NA_character_
      ),
      # Latest last-contact date across SA and EOS
      lastdt = {
        x <- pmax(lastdt_sa, lastdt_eos, na.rm = TRUE)
        x[is.infinite(x)] <- NA
        as.Date(x, origin = "1970-01-01")
      }
    )

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
  pfs_event <- as.Date(pfs_event, origin = "1970-01-01")

  # PFS source: determine whether progression or death triggered the event,
  # then report which form (DP, EOS, or SA) provided that date
  pfs_prog_wins <- !is.na(base$progdt) & (is.na(base$deathdt) | base$progdt <= base$deathdt)
  pfs_srcdom <- ifelse(pfs_prog_wins, base$progdt_src,
                       ifelse(!is.na(base$deathdt), base$deathdt_src, NA_character_))
  pfs_srcvar <- dplyr::case_when(
    pfs_prog_wins & base$progdt_src == "DP"  ~ "DPDAT",
    pfs_prog_wins & base$progdt_src == "EOS" ~ "EOSPDDAT",
    base$deathdt_src == "SA"                 ~ "SADTDAT",
    base$deathdt_src == "EOS"                ~ "EOSDTDAT",
    TRUE                                     ~ NA_character_
  )

  pfs <- make_tte(
    base, "PFS", "Progression-free survival", pfs_event, base$lastdt,
    evntdesc = "Progression or death", srcdom = pfs_srcdom, srcvar = pfs_srcvar
  )
  os <- make_tte(
    base, "OS", "Overall survival", base$deathdt, base$lastdt,
    evntdesc = "Death",
    srcdom = ifelse(!is.na(base$deathdt), base$deathdt_src, NA_character_),
    srcvar = ifelse(!is.na(base$deathdt) & base$deathdt_src == "EOS", "EOSDTDAT", "SADTDAT")
  )

  pfs6 <- make_tte(
    base, "PFS6M", "6-month progression-free survival", pfs_event, base$lastdt, cutoff_days = 180L,
    evntdesc = "Progression or death by day 180", srcdom = pfs_srcdom, srcvar = pfs_srcvar
  )
  os12 <- make_tte(
    base, "OS12M", "Overall survival at 12 months", base$deathdt, base$lastdt, cutoff_days = 365L,
    evntdesc = "Death by day 365",
    srcdom = ifelse(!is.na(base$deathdt), base$deathdt_src, NA_character_),
    srcvar = ifelse(!is.na(base$deathdt) & base$deathdt_src == "EOS", "EOSDTDAT", "SADTDAT")
  )
  os24 <- make_tte(
    base, "OS24M", "Overall survival at 24 months", base$deathdt, base$lastdt, cutoff_days = 730L,
    evntdesc = "Death by day 730",
    srcdom = ifelse(!is.na(base$deathdt), base$deathdt_src, NA_character_),
    srcvar = ifelse(!is.na(base$deathdt) & base$deathdt_src == "EOS", "EOSDTDAT", "SADTDAT")
  )

  adtte <- bind_rows(pfs, pfs6, os, os12, os24) %>%
    mutate(domain = "ADTTE") %>%
    select(
      studyid, domain, usubjid, subjid, cohort, cohortcd, arm, armcd,
      trt01p, trt02p, dose01p, dose02p, astdt, astdy, adt, ady,
      aval, avalu, cnsr, evntdesc, evntstat, srcdom, srcvar,
      paramcd, param, cutoff_days, randdt, rfstdt,
      progdt, progdt_src, deathdt, deathdt_src, dthcaus, lastdt
    )

  adtte
}

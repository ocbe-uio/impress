

`%||%` <- function(a, b) if (!is.null(a)) a else b

get_raw <- function(raw, name) {
  i <- which(tolower(raw$id) == tolower(name))
  if (length(i) == 0L) return(NULL)
  raw$data[[i[1]]]
}

make_adsl <- function(raw, cfg_adam, cfg, identifiers) {

  adsl_spec <- cfg_adam$datasets$ADSL
  if (is.null(adsl_spec)) stop("No ADSL section found in adam.yml")

  # ---------- collect sources from YAML ----------
  get_sources <- function(var) {
    srcs <- var$sources
    if (is.null(srcs)) return(NULL)
    map(srcs, function(s) {
      if (is.list(s)) tibble(dataset = tolower(s$dataset), id = s$id)
      else            tibble(dataset = NA_character_, id = as.character(s))
    }) %>% purrr::list_rbind()
  }

  all_sources <- map(adsl_spec$variables, get_sources) %>% compact() %>%
    purrr::list_rbind() %>% distinct()

  if (nrow(all_sources) && any(is.na(all_sources$dataset))) {
    missing_ds <- unique(all_sources$id[is.na(all_sources$dataset)])
    stop("Some sources in ADSL YAML lack a 'dataset': ",
         paste(missing_ds, collapse = ", "),
         "\nPlease specify both id and dataset (e.g., {id: \"sex\", dataset: \"DM\"}).")
  }

  needed <- all_sources %>%
    group_by(dataset) %>%
    summarise(cols = list(unique(id)), .groups = "drop")

  have_ids <- tolower(raw$id)

  # ---------- subject index from DM ----------
  dm <- get_raw(raw, "dm")
  if (is.null(dm)) stop("DM not found in 'raw'.")
  if (!all(identifiers %in% names(dm))) {
    stop("Identifiers ", paste(identifiers, collapse=", "),
         " not found in DM.")
  }

  subject_index <- dm %>% select(all_of(identifiers)) %>% distinct()
  wide <- subject_index

  # ---------- pull all needed source cols and join ----------
  for (i in seq_len(nrow(needed))) {
    ds_name <- needed$dataset[i]
    cols    <- needed$cols[[i]]
    df <- get_raw(raw, ds_name)
    if (is.null(df)) stop(sprintf("Dataset '%s' in YAML not found in 'raw'.", ds_name))
    df <- df %>% select(any_of(c(identifiers, cols))) %>% distinct()
    wide <- wide %>% left_join(df, by = identifiers)
  }

  # ---------- add RFSTDT_EVT from event_dates (Cycle 1) ----------
  evt <- get_raw(raw, "event_dates")
  if (!is.null(evt)) {
    evt_cycle1 <- evt %>%
      mutate(eventinitiateddate = as.Date(eventinitiateddate),
             eventname_norm = str_squish(str_to_lower(eventname))) %>%
      filter(eventname_norm == "cycle 1") %>%
      group_by(across(all_of(identifiers))) %>%
      summarise(
        RFSTDT_EVT = {
          d <- eventinitiateddate[!is.na(eventinitiateddate)]
          if (length(d)) min(d) else as.Date(NA)
        },
        .groups = "drop"
      )
    wide <- wide %>% left_join(evt_cycle1, by = identifiers)
  } else {
    wide <- wide %>% mutate(RFSTDT_EVT = as.Date(NA))
  }

  # ---------- output shell + labels ----------
  id_first <- identifiers[[1]]
  out <- subject_index %>% transmute(USUBJID = as.character(.data[[id_first]]))
  labels <- list(USUBJID = "Unique Subject Identifier")

  cast_to <- function(x, type) {
    type <- tolower(type %||% "char")
    if (type == "char")      return(as.character(x))
    if (type == "num")       return(suppressWarnings(as.numeric(x)))
    if (type == "date")      return(suppressWarnings(as.Date(x)))
    if (type == "datetime")  return(suppressWarnings(as.POSIXct(x, tz = "UTC")))
    as.character(x)
  }

  # ---------- derivation env ----------
  eval_env <- rlang::env(parent = baseenv())
  eval_env$cfg <- cfg
  for (nm in names(wide)) eval_env[[nm]] <- wide[[nm]]
  for (nm in names(out))  eval_env[[nm]] <- out[[nm]]

  # Optional: expose evt_date() if your YAML calls it directly
  # (here we assume you derive RFSTDT from RFSTDT_EVT in YAML, so not required)
  if (!is.null(evt)) {
    eval_env$EVT <- evt
    eval_env$evt_date <- function(USUBJID, EVENTNAME) {
      rows <- evt[evt[[id_first]] == USUBJID & tolower(trimws(evt$eventname)) == tolower(trimws(EVENTNAME)), , drop = FALSE]
      if (!nrow(rows)) return(as.Date(NA))
      d <- suppressWarnings(as.Date(rows$eventinitiateddate))
      d <- d[!is.na(d)]
      if (!length(d)) return(as.Date(NA))
      min(d)
    }
  }

  # ---------- build variables per YAML ----------
  for (v in adsl_spec$variables) {
    vname  <- v$name
    vlabel <- v$label  %||% vname
    vtype  <- v$type   %||% "char"
    vlen   <- v$length %||% NULL

    col <- rep(NA, nrow(wide))

    if (!is.null(v$sources) && is.null(v$derivation)) {
      vals <- map(v$sources, function(s) {
        stopifnot(is.list(s) && !is.null(s$id) && !is.null(s$dataset))
        pref_col <- s$id
        if (!(pref_col %in% names(wide))) {
          stop(sprintf("Source column %s not present. Check YAML source for %s.", pref_col, vname))
        }
        wide[[pref_col]]
      })
      col <- if (length(vals) == 1) vals[[1]] else do.call(dplyr::coalesce, vals)
    }

    if (!is.null(v$derivation)) {
      for (nm in names(out)) eval_env[[nm]] <- out[[nm]]  # allow dependency on earlier vars
      expr <- rlang::parse_expr(v$derivation)
      res  <- try(rlang::eval_tidy(expr, env = eval_env), silent = TRUE)
      if (inherits(res, "try-error")) {
        stop(sprintf("Failed to evaluate derivation for %s: %s", vname, as.character(res)))
      }
      col <- res
    }

    col <- cast_to(col, vtype)
    if (!is.null(vlen) && !is.na(vlen) && tolower(vtype) == "char") {
      col <- ifelse(is.na(col), col, stringr::str_sub(col, 1L, vlen))
    }

    out[[vname]] <- col
    labels[[vname]] <- vlabel
  }

  # If USUBJID not explicitly derived, keep identifier value
  if (!("USUBJID" %in% names(out))) {
    out <- out %>% mutate(USUBJID = as.character(subject_index[[id_first]]))
    labels[["USUBJID"]] <- "Unique Subject Identifier"
  }

  out %>% select(adsl_spec$keys, dplyr::everything()) %>%
    labelled::set_variable_labels(.labels = labels)
}

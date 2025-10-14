# This code makes the adsl.rds file for the Impress trial

# library(tidyverse)
# source("src/external/functions.R")

# library(purrr)
# library(dplyr)



`%||%` <- function(a, b) if (!is.null(a)) a else b

make_adsl <- function(raw, cfg_adam, cfg, identifiers ) {
  

    adsl_spec <- cfg_adam$datasets$ADSL
    if (is.null(adsl_spec)) stop("No ADSL section found in adam.yml")
  
    # ---- 2) Determine all source fields needed from RAWs
    get_sources <- function(var) {
      srcs <- var$sources
      if (is.null(srcs)) return(NULL)
      # each element is either a scalar id or a list with id+dataset
      map(srcs, function(s) {
        if (is.list(s)) {
          tibble(dataset = s$dataset, id = s$id)
        } else {
          # if dataset omitted, we can't fetch; encourage explicit dataset in YAML
          tibble(dataset = NA_character_, id = as.character(s))
        }
      }) %>% list_rbind()
    }
    
    all_sources <- map(adsl_spec$variables, get_sources) %>% compact() %>% list_rbind() %>% distinct()
    # sanity: require dataset for each source
    if (nrow(all_sources) && any(is.na(all_sources$dataset))) {
      missing_ds <- unique(all_sources$id[is.na(all_sources$dataset)])
      stop("Some sources in ADSL YAML lack a 'dataset': ", paste(missing_ds, collapse=", "),
           "\nPlease specify both id and dataset for each source (e.g., {id: \"sex\", dataset: \"DM\"}).")
    }
    
    # per raw dataset, the set of needed columns
    needed <- all_sources %>% group_by(dataset) %>% summarise(cols = list(unique(id)), .groups = "drop")
    
    # ---- 3) Build a wide subject-level frame with all required source columns
    # Start with union of subjects across all involved datasets
    subject_index <- pick(raw, "dm") %>% 
      select(all_of(identifiers)) %>%
      distinct()
    
    wide <- subject_index
  
    for (i in seq_len(nrow(needed))) {
      ds <- needed$dataset[i]
      cols <- needed$cols[[i]]
      if (!(tolower(ds) %in% raw$id)) stop(sprintf("Dataset '%s' referenced in YAML not found in 'raw' list.", ds))
      df <- pick(raw,tolower(ds)) %>%
        select(all_of(c(identifiers, cols))) %>%
        distinct()
      # prefix columns with ds. to avoid name collisions across datasets with same id
      #df <- df %>% rename_with(~ paste0(ds, ".", .x), -SUBJECT_KEY)
      wide <- wide %>% left_join(df, by = identifiers)
    }
    
 if ("event_dates" %in% raw$id) {
      evt_cycle1 <- pick(raw, "event_dates") %>%
        mutate(eventinitiateddate = as.Date(eventinitiateddate)) %>%
        filter(stringr::str_squish(stringr::str_to_lower(eventname)) == "cycle 1") %>%
        group_by(across(all_of(identifiers))) %>%
        summarise(
          RFSTDT_EVT = {
            dates <- eventinitiateddate[!is.na(eventinitiateddate)]
            if (length(dates)) min(dates) else as.Date(NA)
          },
          .groups = "drop"
        )
      wide <- wide %>%
        left_join(evt_cycle1, by = identifiers)
    } else {
      wide <- wide %>%
        mutate(RFSTDT_EVT = as.Date(NA))
    }
    # ---- 4) Compute each ADSL variable
    out <- tibble(USUBJID = subject_index) # temporary; will overwrite if spec contains USUBJID derivation
    labels <- list()
    
    cast_to <- function(x, type) {
      type <- tolower(type %||% "char")
      if (type == "char")      return(as.character(x))
      if (type == "num")       return(suppressWarnings(as.numeric(x)))
      if (type == "date")      return(suppressWarnings(as.Date(x)))
      if (type == "datetime")  return(suppressWarnings(as.POSIXct(x, tz = "UTC")))
      # default
      as.character(x)
    }
    
    # Bring cfg constants in environment for derivations
    eval_env <- rlang::env(parent = baseenv())
    # expose cfg as a list
    eval_env$cfg <- cfg
    
    for (v in adsl_spec$variables) {
      vname  <- v$name
      vlabel <- v$label %||% vname
      vtype  <- v$type  %||% "char"
      vlen   <- v$length %||% NULL
      vorig  <- v$origin %||% "Collected"
      vcode  <- v$codelist %||% NULL
      
      # default column
      col <- rep(NA, nrow(wide))
      
      # 4a) Collected: copy from first source (or combine if you want to extend)
      if (!is.null(v$sources) && is.null(v$derivation)) {
        # If only one source, use it; if multiple, take first non-missing by row
        vals <- map(v$sources, function(s) {
          stopifnot(is.list(s) && !is.null(s$id) && !is.null(s$dataset))
          pref_col <- s$id
          if (!(pref_col %in% names(wide))) {
            stop(sprintf("Source column %s not present in assembled data. Check YAML source for %s.", pref_col, vname))
          }
          wide[[pref_col]]
        })
        if (length(vals) == 1) {
          col <- vals[[1]]
        } else {
          # rowwise coalesce across multiple source columns
          col <- do.call(coalesce, vals)
        }
      }
      
      # 4b) Derived: evaluate expression in a safe environment where
      # current data columns are available as symbols
      if (!is.null(v$derivation)) {
  
        for (nm in names(wide)) eval_env[[nm]] <- wide[[nm]]
        # also expose already-computed ADSL columns so expressions can reference earlier vars
        for (nm in names(out))  eval_env[[nm]] <- out[[nm]]
        
        expr <- rlang::parse_expr(v$derivation)
        res <- try(rlang::eval_tidy(expr, env = eval_env),  silent = TRUE)
        if (inherits(res, "try-error")) {
          stop(sprintf("Failed to evaluate derivation for %s: %s", vname, as.character(res)))
        }
        col <- res
      }
      
      # type cast
      col <- cast_to(col, vtype)
      
      # apply length (soft: truncate chars)
      if (!is.null(vlen) && !is.na(vlen) && tolower(vtype) == "char") {
        col <- ifelse(is.na(col), col, str_sub(col, 1L, vlen))
      }
      
      out[[vname]] <- col
      labels[[vname]] <- vlabel
    }
    
    # If USUBJID wasn't explicitly built, try to keep the SUBJECT_KEY as USUBJID
    if (!("USUBJID" %in% names(out))) {
      out <- out %>% mutate(USUBJID = as.character(wide$SUBJECT_KEY))
      labels[["USUBJID"]] <- "Unique Subject Identifier"
    }
    
    # ---- 5) Keep only specified keys order (USUBJID first) + any others
    # Place USUBJID first if present
    out <- out %>% 
      select(adsl_spec$keys, everything())
    # ensure USUBJID is first column
    
    # ---- 6) Set variable labels
    out <- labelled::set_variable_labels(out, .labels = labels)
    
    # clean names to lower snake_case if you prefer; for ADaM we usually keep upper-case
    out
}
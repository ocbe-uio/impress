# install.packages(c("dplyr","purrr","tibble","stringr"))
library(dplyr)
library(purrr)
library(tibble)
library(stringr)

# Map your file ids to the dataset names used in YAML (ADEFF sources)
# Adjust as needed for your project.
.default_domain_map <- c(
  # Safety / events
  ae   = "AE",        # Adverse Events
  sa   = "SA",        # Serious AE listing (custom; AE subset) — keep separate if you like

  # Meds / exposure
  cm   = "CM",        # Concomitant Medication
  ex   = "EX",        # Exposure
  exbm = "EXBM",      # Exposure biomarkers (custom) — clarify if needed
  atc_without_ddd = "ATC",  # ATC ref table (aux)

  # Disposition / protocol
  ds   = "DS",        # Disposition
  dv   = "DV",        # Protocol Deviations
  eot  = "DS",        # End of Treatment (maps into DS)
  eos  = "DS",        # End of Study (maps into DS)
  ie   = "IE",        # Inclusion/Exclusion Criteria
  elig = "IE",        # Eligibility (alias to IE)

  # Demography / visits
  dm   = "DM",        # Demographics / randomization data
  sv   = "SV",        # Subject Visits
  vi   = "SV",        # Visit info (alias to SV)
  event_dates = "EVT",# Custom event date table (custom domain EVT)

  # Vitals / clinical scales / exams
  vs   = "VS",        # Vitals
  ecog = "VS",        # ECOG scores (treat as VS or a custom scale domain if you prefer)
  kps  = "VS",        # Karnofsky (same reasoning)
  pe   = "PE",        # Physical Exam

  # Imaging / tumor assessments
  mri  = "RS",        # Response/assessment (MRI-derived measures)
  ct   = "CM",        # Other Cancer therapy
  tr   = "TR",        # Tumor/Lesion Identification (SDTM TR equivalent)

  # Labs / biomarkers / ECG
  lbe  = "LB",        # Labs (check: if “LB Events” aggregate; still map to LB)
  lb   = "LB",        # (in case you ever have a plain lb)
  ecg  = "EG",        # ECG (SDTM domain EG)

  # Questionnaires / PROs
  qlq   = "QS",       # EORTC QLQ scales → QS
  qlqbn = "QS",       # Brain module → QS
  nano  = "QS",       # NANO neurological scale → QS
  sq    = "QS",       # Study questionnaire (generic) → QS

  # Medical history & others
  mh   = "MH",        # Medical History
  pre  = "PRE",       # Pre-study info (custom)
  pt   = "PT",        # Performance test? (custom; clarify later)

  # Trial management / reference
  ran  = "RAND",      # Randomization export (custom). You can also map to DM if preferred.
  ch   = "CH",        # Custom (clarify later)
  chb  = "CHB",       # Custom
  co   = "CO",        # Custom
  cs   = "CS",        # Custom
  dp   = "DP",        # Custom
  ul   = "UL",        # Custom
  meddra   = "MEDDRA"    # Reference dictionary
)

# raw_index: your tibble with columns id (like "dm","ex","mri",...) and list-col 'data'
as_raw_domains <- function(raw_index,
                           id_col   = "id",
                           data_col = "data",
                           domain_map = .default_domain_map,
                           subject_col = "subjectid",
                           drop_ids = c("codelist", "items")) {

  stopifnot(id_col %in% names(raw_index), data_col %in% names(raw_index))

  keep <- !tolower(raw_index[[id_col]]) %in% drop_ids
  raw_index <- raw_index[keep, , drop = FALSE]

  # build named list: names = mapped domain (e.g., "DM","EX","RS")
  out <- vector("list", nrow(raw_index))
  nm  <- character(nrow(raw_index))

  for (i in seq_len(nrow(raw_index))) {
    id_val <- tolower(as.character(raw_index[[id_col]][[i]]))
    dsname <- domain_map[[id_val]] %||% toupper(id_val)  # fallback to uppercased id

    df <- raw_index[[data_col]][[i]]
    # Ungroup if grouped_df, keep as tibble
    if (inherits(df, "grouped_df")) df <- dplyr::ungroup(df)
    df <- tibble::as_tibble(df)

    # sanity: ensure subjectid exists
    if (!subject_col %in% names(df)) {
      stop(sprintf("Domain '%s' (file id '%s') lacks '%s' column.", dsname, id_val, subject_col))
    }

    out[[i]] <- df
    nm[i]    <- dsname
  }

  # consolidate duplicates (e.g., multiple rows mapping to VS): bind rows
  domains <- split(out, nm) %>%
    purrr::imap(~ dplyr::bind_rows(.x))  # .x is a list of dfs for that name

  domains
}

# ----------------------------------------------------------------------------
# Helper functions for defining visits
#---------------------------------------------------------------------------

visits_from_yaml <- function(ln) {
  v <- ln$visits
  if (is.null(v) || is.null(v$map)) stop("No 'visits: map:' section in adam.yml")
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

lookup_avisit_yaml  <- function(eventname, MAP, fallback_title = TRUE) {
  x <- tibble(eventname = str_squish(as.character(eventname))) %>%
       left_join(MAP, by = "eventname")
  out <- x$avisit
  if (fallback_title) out[is.na(out)] <- str_to_title(x$eventname[is.na(out)])
  out
}

lookup_avisitn_yaml <- function(eventname, MAP, unmapped = 999L) {
  x <- tibble(eventname = str_squish(as.character(eventname))) %>%
       left_join(MAP, by = "eventname")
  out <- x$avisitn
  out[is.na(out)] <- as.integer(unmapped)
  as.integer(out)
}



#
`%||%` <- function(a,b) if (!is.null(a)) a else b

# ------------------------------------------------------------------------------
# make_adeff(): build ADEFF (BDS) from raw domains and lineage YAML
# ------------------------------------------------------------------------------

make_adeff <- function(
    lineage_yml,
    raw,
    adsl = NULL,
    cfg = list(studyid = "ImPRESS")
) {
  
  stopifnot(file.exists(lineage_yml))
  ln <- yaml::read_yaml(lineage_yml)
  adeff_spec <- ln$datasets$ADEFF
  if (is.null(adeff_spec)) stop("No ADEFF section found in YAML.")

  raw_list <- raw
  if (!is.null(adsl)) {
    derive_subjectid <- function(df) {
      if ("subjectid" %in% names(df)) return(df$subjectid)
      if ("SUBJECTID" %in% names(df)) return(df$SUBJECTID)
      if ("SUBJID" %in% names(df)) return(df$SUBJID)
      if ("USUBJID" %in% names(df)) {
        return(stringr::str_trim(sub("^[^-]+-\\s*", "", df$USUBJID)))
      }
      stop("Unable to derive subject identifier from supplied ADSL data.")
    }
    subj_vals <- derive_subjectid(adsl)
    adsl_prepped <- adsl %>%
      mutate(subjectid = subj_vals)
    raw_list[["ADSL"]] <- adsl_prepped
  }

  subject_keys <- setNames(rep("subjectid", length(raw_list)), names(raw_list))
  
  # --------------------------------------------------------------------------
  # 1) Detect subject key per raw dataset
  guess_subj_col <- function(df) {
    cand <- c("USUBJID","SUBJID","subject_id","SubjectID","SUBJECT_ID")
    hit <- cand[cand %in% names(df)]
    if (length(hit)) hit[1] else NA_character_
  }
  
  raws <- imap(raw_list, function(df, nm) {
    df <- as_tibble(df)
    subj_col <- if (!is.null(subject_keys) && nm %in% names(subject_keys)) subject_keys[[nm]] else guess_subj_col(df)
    if (is.na(subj_col)) {
      stop(sprintf("Can't find subject id column for dataset '%s'.", nm))
    }
    rename(df, SUBJECT_KEY = !!sym(subj_col))
  })
  
  # --------------------------------------------------------------------------
  # 2) Gather sources from YAML
  get_sources <- function(var) {
    srcs <- var$sources
    if (is.null(srcs)) return(NULL)
    map(srcs, ~tibble(dataset = .x$dataset, id = .x$id))
  }
  all_sources <- map(adeff_spec$variables, get_sources) %>% compact() %>% purrr::list_flatten() |>  purrr::list_rbind()
  
  needed <- all_sources %>%
    group_by(dataset) %>%
    summarise(cols = list(unique(id)), .groups="drop")
  
  # --------------------------------------------------------------------------
  # 3) Merge source datasets
  subject_index <- map(raws, ~ select(.x, SUBJECT_KEY) %>% distinct()) %>% purrr::list_rbind() %>% distinct()
  wide <- subject_index
  
  for (i in seq_len(nrow(needed))) {
    ds <- needed$dataset[i]
    cols <- needed$cols[[i]]
    if (!(ds %in% names(raws))) stop(sprintf("Dataset '%s' not found in raw list.", ds))
    df <- raws[[ds]] %>%
      select(SUBJECT_KEY, all_of(cols)) %>%
      distinct()
    df <- rename_with(df, ~ paste0(ds, ".", .x), -SUBJECT_KEY)
    wide <- left_join(wide, df, by="SUBJECT_KEY")
  }
  
  # alias single-source vars (optional)
  alias_map <- all_sources %>% count(id, dataset) %>% group_by(id) %>% filter(n()==1) %>% ungroup()
  for (j in seq_len(nrow(alias_map))) {
    id <- alias_map$id[j]; ds <- alias_map$dataset[j]
    src <- paste0(ds, ".", id)
    if (src %in% names(wide) && !(id %in% names(wide))) wide[[id]] <- wide[[src]]
  }
  
  # --------------------------------------------------------------------------
  # 4) Evaluate derivations
  vm <- visits_from_yaml(ln)
  eval_env <- rlang::env(parent = baseenv())
  eval_env$cfg <- cfg
  
  # Expose lookups for YAML derivations
  eval_env$VISIT_MAP  <- vm$map
  eval_env$lookup_avisit  <- function(eventname)  
    lookup_avisit_yaml(eventname, MAP = eval_env$VISIT_MAP,
                  fallback_title = vm$defaults$avisit_fallback_from_eventname)
  eval_env$lookup_avisitn <- function(eventname)  
    lookup_avisitn_yaml(eventname, MAP = eval_env$VISIT_MAP,
                        unmapped = vm$defaults$avisitn_unmapped)


  out <- tibble(USUBJID = wide$SUBJECT_KEY)
  labels <- list()
  
  cast_to <- function(x, type) {
    type <- tolower(type %||% "char")
    if (type == "char") return(as.character(x))
    if (type == "num")  return(as.numeric(x))
    if (type == "date") return(as.Date(x))
    if (type == "datetime") return(as.POSIXct(x, tz="UTC"))
    as.character(x)
  }
  
  intermediate_vars <- character()

  for (v in adeff_spec$variables) {
    vname <- v$name; vlabel <- v$label %||% vname; vtype <- v$type %||% "char"; vlen <- v$length %||% NULL
    col <- rep(NA, nrow(wide))
    
    # collected
    if (!is.null(v$sources) && is.null(v$derivation)) {
      vals <- map(v$sources, function(s) {
        pref <- paste0(s$dataset, ".", s$id)
        if (!(pref %in% names(wide))) stop(sprintf("Missing column %s", pref))
        wide[[pref]]
      })
      col <- if (length(vals)==1) vals[[1]] else do.call(coalesce, vals)
    }
    
    # derived
    if (!is.null(v$derivation)) {
      for (nm in names(wide)) eval_env[[nm]] <- wide[[nm]]
      for (nm in names(out))  eval_env[[nm]] <- out[[nm]]
      expr <- parse_expr(v$derivation)
      res <- try(eval(expr, envir = eval_env, enclos = baseenv()), silent=TRUE)
      if (inherits(res,"try-error")) stop(sprintf("Failed derivation for %s: %s", vname, as.character(res)))
      col <- res
    }
    
    col <- cast_to(col, vtype)
    if (!is.null(vlen) && !is.na(vlen) && tolower(vtype)=="char")
      col <- ifelse(is.na(col), col, str_sub(col, 1L, vlen))
    
    out[[vname]] <- col
    labels[[vname]] <- vlabel
    if (isTRUE(v$intermediate)) {
      intermediate_vars <- c(intermediate_vars, vname)
    }
  }

  if (length(intermediate_vars)) {
    keep_cols <- setdiff(names(out), intermediate_vars)
    out <- out[, keep_cols, drop = FALSE]
    labels <- labels[names(labels) %in% keep_cols]
  }
  
  # --------------------------------------------------------------------------
  # 5) Reshape: each PARAMCD x AVISIT
  # assumes PARAMCD and AVISIT (and maybe PARAM) exist
  # if PARAMCD not yet long, you can pivot_longer here.
  # For now, return "as is"; you'll normally have one row per PARAMCD × visit already.
  
  out <- labelled::set_variable_labels(out, .labels = labels)
  out
}

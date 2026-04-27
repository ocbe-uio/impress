# ADLB dataset creation — Bio-Plex multiplex cytokine/angiogenesis panels
#
# Source files (data/raw/biomarkers/):
#   August 2025 Summary PL1-4.xlsx  — measured concentrations (pg/mL), wide format
#   Boks 1-3 Impress Losartan.xlsx  — sample inventory: barcode → donor ID + collection date + visit
#   Prøve ID og plate fordeling.xlsx — plate mapping: barcode → patient ID (fallback)
#
# ADLB variable mapping:
#   --- Subject identifiers ---
#   STUDYID   ← cfg$studyid  ("ImPRESS")
#   DOMAIN    ← "LB"
#   USUBJID   ← cfg$studyid + "-" + subjid  (joined from adsl)
#   SUBJID    ← Donor Id (Boks) or Pasient ID (plate mapping), normalised to "ImPR-XXXX"
#   --- Demographics (from adsl) ---
#   AGE       ← age at enrolment
#   SEX       ← sex
#   --- Treatment / dosing (from adsl) ---
#   COHORT    ← cohort label
#   COHORTCD  ← cohort code
#   ARM       ← treatment arm description
#   ARMCD     ← treatment arm code
#   STEP      ← treatment step (1, 2, 3)
#   DOSE01P   ← planned dose in step 1 (mg)
#   DOSE01PU  ← unit for DOSE01P
#   TRT01P    ← planned treatment label in step 1
#   DOSE02P   ← planned follow-up dose (mg)
#   DOSE02PU  ← unit for DOSE02P
#   TRT02P    ← planned follow-up treatment label
#   RANDDT    ← randomisation date
#   RFSTDT    ← first treatment date (reference for ADY)
#   TRT       ← stepped-wedge actual treatment at visit ("0 mg", "25 mg", …)
#   TRTCD     ← numeric dose for TRT (factor: 0, 25, 50, 100)
#   --- ADaM BDS parameter (required) ---
#   PARAMCD   ← = LBTESTCD  (analysis parameter code; ≤ 8 chars)
#   PARAM     ← = LBTEST    (analysis parameter label)
#   PARCAT1   ← panel  ("19-PLEX" | "14-PLEX"; BDS name for LBCAT)
#   --- SDTM LB carry-forwards (for traceability) ---
#   LBTESTCD  ← short code from analyte column name  (e.g. "IL6", "VEGFA")
#   LBTEST    ← analyte name stripped of bead-lot number  (e.g. "IL-6", "VEGF-A")
#   LBORRES   ← raw character value from summary sheet  (preserves "*" and "OOR>" flags)
#   LBORRESU  ← "pg/mL"
#   LBSTRESC  ← cleaned character result: "*" stripped; "OOR>" → ">ULOQ"; "OOR<" → "<LLOQ"
#   LBSTRESN  ← numeric result (NA for OOR/outlier rows)
#   LBSTRESU  ← "pg/mL"
#   LBNRIND   ← NA  (no reference ranges available in source)
#   --- ADaM BDS analysis values (required) ---
#   AVAL      ← = LBSTRESN  (numeric analysis value; NA for OOR/outlier)
#   AVALC     ← = LBSTRESC  (character analysis value)
#   AVALU     ← = LBSTRESU  ("pg/mL")
#   --- Flags, derived analysis, and timing ---
#   ABLFL     ← "Y" when ADT <= RFSTDT (any pre-treatment sample; ADaM BDS baseline flag)
#   BASE      ← AVAL at baseline (Screening) for each subject × analyte
#   CHG       ← AVAL − BASE
#   PCHG      ← 100 × CHG / BASE (NA when BASE == 0)
#   AVISIT    ← canonical visit label from cfg visit map; falls back to title-case event name
#   AVISITN   ← derived from event_dates via visit schedule (e.g. -14 Screening, 1 Cycle 1, 15 Cycle 2)
#   ADT       ← Collection Date from Boks inventory
#   ADTM      ← NA  (sample-level time not available)
#   ADY       ← derived via admiral::derive_vars_dy, relative to rfstdt from adsl

make_adlb <- function(raw, adsl, cfg) {

  biomarker_dir <- "data/raw/biomarkers"

  visits <- visits_from_yaml(cfg)

  # ---- Analyte metadata: column name in summary file → ADLB codes ----------
  analyte_map <- tibble::tribble(
    ~col_name,                   ~lbtestcd,   ~lbtest,                ~lbcat,
    # 19-plex panel
    "Angiopoietin-1 (64)",       "ANGPT1",    "Angiopoietin-1",       "19-PLEX",
    "Angiopoietin-2 (26)",       "ANGPT2",    "Angiopoietin-2",       "19-PLEX",
    "ANGPTL4 (73)",              "ANGPTL4",   "ANGPTL4",              "19-PLEX",
    "CA9 (45)",                  "CA9",       "CA9",                  "19-PLEX",
    "CD163 (28)",                "CD163",     "CD163",                "19-PLEX",
    "Collagen IV alpha 1 (55)",  "COLIVA1",   "Collagen IV alpha 1",  "19-PLEX",
    "MIG (52)",                  "MIG",       "MIG (CXCL9)",          "19-PLEX",
    "IP-10 (21)",                "IP10",      "IP-10 (CXCL10)",       "19-PLEX",
    "I-TAC (63)",                "ITAC",      "I-TAC (CXCL11)",       "19-PLEX",
    "EGF (25)",                  "EGF",       "EGF",                  "19-PLEX",
    "ErbB2 (75)",                "ERBB2",     "ErbB2",                "19-PLEX",
    "TNFSF6 (39)",               "TNFSF6",    "TNFSF6 (FasL)",        "19-PLEX",
    "FGF Basic (47)",            "FGFB",      "FGF Basic",            "19-PLEX",
    "TNFRSF18 (61)",             "TNFRSF18",  "TNFRSF18 (GITR)",      "19-PLEX",
    "IL-2 (43)",                 "IL2",       "IL-2",                 "19-PLEX",
    "IL-12 p70 (56)",            "IL12P70",   "IL-12 p70",            "19-PLEX",
    "VEGFR2 (12)",               "VEGFR2",    "VEGFR2",               "19-PLEX",
    "Tenascin C (35)",           "TENASC",    "Tenascin C",           "19-PLEX",
    "PIGF (72)",                 "PIGF",      "PIGF",                 "19-PLEX",
    # 14-plex panel
    "MCP-1 (25)",                "MCP1",      "MCP-1",                "14-PLEX",
    "PECAM-1 (45)",              "PECAM1",    "PECAM-1",              "14-PLEX",
    "Fractalkine (46)",          "FRACTKN",   "Fractalkine",          "14-PLEX",
    "Granzyme B (57)",           "GZMB",      "Granzyme B",           "14-PLEX",
    "IL-1b (28)",                "IL1B",      "IL-1b",                "14-PLEX",
    "IL-6 (13)",                 "IL6",       "IL-6",                 "14-PLEX",
    "IL-8 (18)",                 "IL8",       "IL-8",                 "14-PLEX",
    "MMP-10 (66)",               "MMP10",     "MMP-10",               "14-PLEX",
    "PDGF-DD (54)",              "PDGFDD",    "PDGF-DD",              "14-PLEX",
    "TGF-alpha (63)",            "TGFA",      "TGF-alpha",            "14-PLEX",
    "Thrombospondin-2 (52)",     "THBS2",     "Thrombospondin-2",     "14-PLEX",
    "VEGF-A (26)",               "VEGFA",     "VEGF-A",               "14-PLEX",
    "VEGF-C (38)",               "VEGFC",     "VEGF-C",               "14-PLEX",
    "IFN-g (29)",                "IFNG",      "IFN-g",                "14-PLEX"
  )

  # ---- 1. Measurements (summary file, wide → long) -------------------------
  # col_types = "text" preserves raw flags: "*" (extrapolated), "OOR>" (above range),
  # "OOR<" (below range), "0" (below detection), "---" (outlier)
  summary_raw <- readxl::read_excel(
    file.path(biomarker_dir, "August 2025 Summary PL1-4.xlsx"),
    sheet     = "pgml",
    col_types = "text"
  )

  # The first column header is a merged cell artefact ("PL 1"); it holds barcodes
  names(summary_raw)[1] <- "barcode"

  meas_long <- summary_raw %>%
    dplyr::select(barcode, dplyr::any_of(analyte_map$col_name)) %>%
    dplyr::filter(
      !is.na(barcode),
      stringr::str_starts(barcode, "E0")   # patient barcodes only; drop QC/note rows
    ) %>%
    tidyr::pivot_longer(
      cols      = -barcode,
      names_to  = "col_name",
      values_to = "lborres"
    ) %>%
    dplyr::filter(!is.na(lborres)) %>%
    dplyr::left_join(analyte_map, by = "col_name") %>%
    dplyr::select(-col_name) %>%
    dplyr::mutate(
      lborresu = "pg/mL",
      # LBSTRESC: human-readable standardised character result
      lbstresc = dplyr::case_when(
        stringr::str_starts(lborres, "OOR>") ~ ">ULOQ",
        stringr::str_starts(lborres, "OOR<") ~ "<LLOQ",
        stringr::str_starts(lborres, "---")  ~ NA_character_,
        TRUE ~ stringr::str_remove(lborres, "^\\*")
      ),
      # LBSTRESN: numeric; OOR rows carry no numeric value
      lbstresn = suppressWarnings(as.numeric(lbstresc)),
      lbstresn = dplyr::if_else(
        lbstresc %in% c(">ULOQ", "<LLOQ"), NA_real_, lbstresn
      ),
      lbstresu = "pg/mL",
      lbnrind  = NA_character_
    )

  # ---- 2. Sample inventory (Boks files → barcode + collection date + visit) -
  read_boks <- function(path) {
    readxl::read_excel(path, sheet = 1) %>%
      dplyr::rename_with(
        ~ stringr::str_to_lower(.) %>%
          stringr::str_replace_all(" ", "_") %>%
          stringr::str_remove_all("[^a-z0-9_]")
      ) %>%
      dplyr::transmute(
        barcode         = as.character(barcode),
        donor_id        = as.character(donor_id),
        # readxl returns collection_date as POSIXct; coerce to Date
        collection_date = as.Date(collection_date),
        visit           = stringr::str_squish(as.character(visit))
      ) %>%
      dplyr::filter(!is.na(barcode), nzchar(barcode))
  }

  boks_files <- list.files(
    biomarker_dir, pattern = "^Boks.*\\.xlsx$", full.names = TRUE
  )

  boks <- purrr::map(boks_files, read_boks) %>%
    dplyr::bind_rows() %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      donor_id = stringr::str_replace(donor_id, "^[Ii][Mm][Pp][Rr]-", "ImPR-")
    )

  # ---- 2a. Event dates: canonical visit names per subject per date ----------
  # avisitn is mapped from the raw event name via the cfg visit schedule so that
  # the date fallback lookup can match by avisitn (numeric) rather than by event
  # name string — robust even when Viedoc event names differ from cfg names.
  event_lookup <- raw |>
    get_raw("event_dates") |>
    dplyr::filter(eventstatus == "Initiated") |>
    dplyr::transmute(
      subjid     = stringr::str_replace(as.character(subjectid), "^[Ii][Mm][Pp][Rr]-", "ImPR-"),
      event_date = as.Date(eventinitiateddate),
      eventname  = stringr::str_squish(as.character(eventname)),
      avisitn    = lookup_avisitn(stringr::str_squish(as.character(eventname)),
                                  visits$map, unmapped = NA_integer_)
    ) |>
    dplyr::filter(!is.na(event_date))

  # ---- 3. Plate mapping (barcode → patient ID fallback) --------------------
  # Sheet layout: four blocks of (DAG, Unik ID, Pasient ID, spacer) across cols 1–16
  plate_raw <- readxl::read_excel(
    file.path(biomarker_dir, "Prøve ID og plate fordeling.xlsx"),
    col_types = "text"
  )
  nms <- names(plate_raw)

  extract_id_block <- function(id_col, pat_col) {
    plate_raw %>%
      dplyr::select(
        unik_id    = dplyr::all_of(id_col),
        pasient_id = dplyr::all_of(pat_col)
      ) %>%
      dplyr::filter(
        !is.na(unik_id), !is.na(pasient_id),
        nzchar(unik_id), nzchar(pasient_id),
        stringr::str_starts(unik_id, "E0")
      )
  }

  plate_map <- dplyr::bind_rows(
    extract_id_block(nms[2],  nms[3]),
    extract_id_block(nms[6],  nms[7]),
    extract_id_block(nms[10], nms[11]),
    extract_id_block(nms[14], nms[15])
  ) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      pasient_id = stringr::str_replace(pasient_id, "^[Ii][Mm][Pp][Rr]-", "ImPR-")
    )

  # ---- 4. ADSL keys --------------------------------------------------------
  adsl_keys <- adsl %>%
    dplyr::select(
      subjid, usubjid, studyid,
      age, sex,
      cohort, cohortcd, arm, armcd,
      step, dose01p, dose01pu, trt01p, dose02p, dose02pu, trt02p,
      randdt, rfstdt
    )

  # ---- 5. Join all sources -------------------------------------------------
  adlb_base <- meas_long %>%
    # Primary subject info from Boks inventory
    dplyr::left_join(boks, by = "barcode") %>%
    # Fallback patient ID from plate mapping where Boks donor_id is missing
    dplyr::left_join(plate_map, by = c("barcode" = "unik_id")) %>%
    dplyr::mutate(
      subjid = dplyr::coalesce(donor_id, pasient_id),
      subjid = stringr::str_replace(subjid, "^[Ii][Mm][Pp][Rr]-", "ImPR-"),
      adt    = collection_date
    ) %>%
    dplyr::filter(!is.na(subjid)) %>%
    # Derive canonical event name from event_dates: nearest initiated event within ±7 days
    dplyr::mutate(
      # Map boks visit labels to canonical cfg event names where they differ.
      # ("Cycle 2", "Cycle 4", "Cycle 17" already match cfg names; only Screening differs)
      visit_mapped = dplyr::recode(visit,
        "Screening" = "Cycle 1 -4 days",
        .default    = visit
      ),
      # avisitn from boks visit label (via cfg map) — used for the date fallback below
      avisitn_boks = lookup_avisitn(visit_mapped, visits$map, unmapped = NA_integer_),
      # (a) When collection date is available: nearest event_dates match within ±7 days.
      #     Only consider events with a mapped avisitn to avoid unmapped events
      #     (e.g. "Study start", "End of Study") winning ties on the same date.
      eventname_ev = purrr::map2_chr(subjid, adt, function(s, d) {
        if (is.na(d)) return(NA_character_)
        ev <- event_lookup[event_lookup$subjid == s & !is.na(event_lookup$avisitn), ]
        if (nrow(ev) == 0L) return(NA_character_)
        diffs <- abs(as.integer(d - ev$event_date))
        i     <- which.min(diffs)
        if (diffs[[i]] > 7L) NA_character_ else ev$eventname[[i]]
      }),
      # (b) When collection date is missing: look up event date by avisitn (numeric match
      #     is robust to raw event name differences between Viedoc and the cfg map)
      adt_from_visit = as.Date(purrr::map2_dbl(subjid, avisitn_boks, function(s, n) {
        if (is.na(n) || is.na(s)) return(NA_real_)
        ev <- event_lookup[
          event_lookup$subjid == s &
          !is.na(event_lookup$avisitn) &
          event_lookup$avisitn == n, ]
        if (nrow(ev) == 0L) return(NA_real_)
        as.numeric(ev$event_date[[1]])
      }), origin = "1970-01-01"),
      # Fill missing collection date from event date (avisitn lookup)
      adt       = dplyr::coalesce(adt, adt_from_visit),
      # Prefer event name from date match; fall back to mapped boks visit label
      visit_key = dplyr::coalesce(eventname_ev, visit_mapped),
      avisit    = lookup_avisit(visit_key, visits$map, fallback_title = TRUE),
      avisitn   = lookup_avisitn(visit_key, visits$map, unmapped = NA_integer_)
    ) %>%
    dplyr::select(-visit_mapped, -avisitn_boks, -eventname_ev, -adt_from_visit, -visit_key, -visit) %>%
    dplyr::left_join(adsl_keys, by = "subjid") %>%
    dplyr::mutate(
      studyid  = cfg$studyid,
      domain   = "LB",
      usubjid  = dplyr::coalesce(usubjid, paste0(cfg$studyid, "-", subjid)),
      # ADaM BDS required parameter variables (= SDTM carry-forwards)
      paramcd  = lbtestcd,
      param    = lbtest,
      parcat1  = lbcat,
      # ADaM BDS analysis values (= SDTM standardised results)
      aval     = lbstresn,
      avalc    = lbstresc,
      avalu    = lbstresu,
      # Any sample collected on or before first treatment date is baseline
      # Use "N" (not NA) for post-treatment rows so filter(ablfl != "Y") works
      ablfl    = dplyr::if_else(!is.na(adt) & !is.na(rfstdt) & adt <= rfstdt, "Y", "N", missing = "N"),
      adtm     = NA_character_
    )

  # ---- 6. Derive ADY via admiral (days relative to first treatment date) ----
  adlb <- admiral::derive_vars_dy(
    dataset        = adlb_base %>% dplyr::rename_with(toupper),
    reference_date = RFSTDT,
    source_vars    = admiral::exprs(ADT)
  ) %>%
    dplyr::rename_with(tolower)

  # ---- 6a. Stepped-wedge treatment at each visit (mirrors adeff logic) -----
  trt_map <- trt_map_from_yaml(cfg)
  adlb <- adlb %>%
    dplyr::mutate(
      trtcd = trt_from_map(armcd, ady, trt_map),
      trtcd = factor(trtcd, levels = c(0, 25, 50, 100)),
      trt   = factor(paste0(trtcd, " mg"), levels = c("0 mg", "25 mg", "50 mg", "100 mg"))
    )

  # ---- 6b. Baseline, change, and percent change from baseline --------------
  # Base = aval from the most recent pre-treatment sample (latest adt where ablfl == "Y")
  adlb <- adlb %>%
    dplyr::group_by(usubjid, paramcd) %>%
    dplyr::mutate(
      baseline_adt = {
        bl <- ablfl == "Y" & !is.na(adt) & !is.na(aval)
        if (any(bl)) max(adt[bl], na.rm = TRUE) else as.Date(NA)
      },
      base = {
        bl_val <- aval[ablfl == "Y" & !is.na(adt) & adt == baseline_adt & !is.na(aval)]
        if (length(bl_val) > 0) bl_val[[length(bl_val)]] else NA_real_
      },
      chg  = dplyr::if_else(!is.na(aval) & !is.na(base), aval - base, NA_real_),
      pchg = dplyr::if_else(!is.na(aval) & !is.na(base) & base != 0,
                            100 * (aval - base) / base, NA_real_)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-baseline_adt)

  # ---- 7. Final column selection and sort ----------------------------------
  adlb <- adlb %>%
    dplyr::select(
      # Subject identifiers
      studyid, domain, usubjid, subjid,
      # Demographics
      age, sex,
      # Treatment / dosing (aligned with ADEFF)
      cohort, cohortcd, arm, armcd,
      step, dose01p, dose01pu, trt01p, dose02p, dose02pu, trt02p,
      randdt, rfstdt,
      trt, trtcd,
      # ADaM BDS parameter (required)
      paramcd, param, parcat1,
      # SDTM LB carry-forwards (traceability)
      lbtestcd, lbtest,
      lborres, lborresu,
      lbstresc, lbstresn, lbstresu,
      lbnrind,
      # ADaM BDS analysis values (required)
      aval, avalc, avalu,
      # Flags and derived
      ablfl, base, chg, pchg,
      # Timing
      avisit, avisitn,
      adt, adtm, ady,
      # Source
      barcode
    ) %>%
    dplyr::arrange(usubjid, parcat1, paramcd, adt)

  adlb
}

# Standalone execution -------------------------------------------------------
# Runs only when this script is called directly (e.g. Rscript make_adlb.R),
# not when sourced by _targets.R or another function.
if (sys.nframe() == 0L) {
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readxl)
  library(readr)
  library(tibble)
  library(lubridate)
  library(yaml)
  library(admiral)
  library(fs)
  library(haven)
  library(labelled)
  library(targets)

  source("R/external/functions.R")
  source("R/make_ad/helpers_ad.R")
  source("R/make_raw/make_raw.R")

  cfg  <- yaml::read_yaml("config/cfg.yml")
  tar_load(raw)
  tar_load(shamraw)
  tar_load(adsl)

  dat  <- effective_raw(raw, shamraw, cfg)
  adlb <- make_adlb(dat, adsl, cfg)

  out_dir <- "data/ad"
  fs::dir_create(out_dir)
  saveRDS(adlb, file.path(out_dir, "adlb.rds"))
  message(
    "Saved ", file.path(out_dir, "adlb.rds"),
    " — ", nrow(adlb), " rows, ",
    dplyr::n_distinct(adlb$usubjid), " subjects, ",
    dplyr::n_distinct(adlb$lbtestcd), " analytes"
  )

  # ---- Validation 1: cross-check against LBE (exploratory sample collection) ----
  # lbelabyn = "Yes" records in LBE confirm a biomarker sample was collected;
  # lbedat is the exploratory collection date. Each such record should align with
  # an ADLB row for the same subject and date, and visit names should agree.
  lbe <- dat |>
    get_raw("lbe") |>
    dplyr::filter(toupper(as.character(lbelabyn)) == "YES", !is.na(lbedat)) |>
    dplyr::transmute(
      subjid    = stringr::str_replace(as.character(subjectid), "^[Ii][Mm][Pp][Rr]-", "ImPR-"),
      lbe_date  = as.Date(lbedat),
      lbe_visit = stringr::str_squish(as.character(eventname))
    )

  adlb_dates <- adlb |>
    dplyr::distinct(subjid, adt, avisit) |>
    dplyr::filter(!is.na(adt))

  lbe_unmatched <- lbe |>
    dplyr::anti_join(adlb_dates, by = c("subjid", "lbe_date" = "adt"))

  if (nrow(lbe_unmatched) > 0) {
    message("WARNING: ", nrow(lbe_unmatched), " LBE exploratory-sample records have no matching ADLB date:")
    print(lbe_unmatched)
  } else {
    message("LBE validation OK: all exploratory-sample dates accounted for in ADLB")
  }

  lbe_visit_mismatch <- lbe |>
    dplyr::inner_join(adlb_dates, by = c("subjid", "lbe_date" = "adt")) |>
    dplyr::filter(
      !is.na(lbe_visit), !is.na(avisit),
      stringr::str_to_lower(lbe_visit) != stringr::str_to_lower(avisit)
    )

  if (nrow(lbe_visit_mismatch) > 0) {
    message("WARNING: ", nrow(lbe_visit_mismatch), " subject-date pairs have mismatched visit labels between LBE and ADLB:")
    print(lbe_visit_mismatch)
  } else {
    message("LBE visit-label validation OK")
  }

  # Validation results: There are 178 of 322 with unmatched dates. We use a ±7-day window to link ADLB dates to event_dates, 
  # so this suggests many exploratory samples were collected outside that window (e.g. unscheduled visits). Visit labels do not match, Visit 1 is used for start visit 2 etc. 
  

  # ---- Validation 2: cross-check VISITNUM against ADEFF ----------------------
  # ADEFF uses avisitn from the same cfg visit schedule; visitnums should agree
  # for any subject+date pair that appears in both datasets.
  adeff <- readRDS(file.path(cfg$paths$adam, "adeff.rds"))

  adeff_visits <- adeff |>
    dplyr::distinct(subjid, adt, avisitn) |>
    dplyr::filter(!is.na(adt), !is.na(avisitn))

  adlb_visits <- adlb |>
    dplyr::distinct(subjid, adt, avisitn) |>
    dplyr::filter(!is.na(adt), !is.na(avisitn))

  visit_conflicts <- adlb_visits |>
    dplyr::inner_join(adeff_visits, by = c("subjid", "adt")) |>
    dplyr::filter(avisitn.x != avisitn.y)

  if (nrow(visit_conflicts) > 0) {
    message("WARNING: ", nrow(visit_conflicts), " subject-date pairs differ in AVISITN between ADLB and ADEFF:")
    print(visit_conflicts)
  } else {
    message("ADEFF visit-number validation OK")
  }

  # ---- Validation 3: CDISC ADaM BDS compliance checks ----------------------
  cdisc_issues <- list()

  # Required BDS identifiers
  if (anyNA(adlb$studyid))  cdisc_issues <- c(cdisc_issues, "STUDYID has NA values")
  if (anyNA(adlb$usubjid))  cdisc_issues <- c(cdisc_issues, "USUBJID has NA values")
  if (any(adlb$domain != "LB", na.rm = TRUE)) cdisc_issues <- c(cdisc_issues, "DOMAIN != 'LB'")
  # Required BDS parameter variables
  if (anyNA(adlb$paramcd))  cdisc_issues <- c(cdisc_issues, "PARAMCD has NA values")
  if (any(nchar(adlb$paramcd) > 8, na.rm = TRUE))
    cdisc_issues <- c(cdisc_issues, paste("PARAMCD exceeds 8 chars:", paste(unique(adlb$paramcd[nchar(adlb$paramcd) > 8]), collapse = ", ")))
  if (anyNA(adlb$param))    cdisc_issues <- c(cdisc_issues, "PARAM has NA values")
  # Required BDS timing
  if (anyNA(adlb$adt))      cdisc_issues <- c(cdisc_issues, "ADT has NA values (missing collection dates)")
  # Required BDS analysis values (at least AVAL or AVALC must be non-NA per row)
  both_na <- is.na(adlb$aval) & is.na(adlb$avalc)
  if (any(both_na))
    cdisc_issues <- c(cdisc_issues, paste(sum(both_na), "rows have both AVAL and AVALC missing"))
  # SDTM carry-forward checks
  if (anyNA(adlb$lborres))  cdisc_issues <- c(cdisc_issues, "LBORRES has NA values")
  # Flag and visit checks
  if (!all(adlb$ablfl %in% c("Y", "N")))
    cdisc_issues <- c(cdisc_issues, "ABLFL contains values other than 'Y' or 'N'")
  if (!is.integer(adlb$avisitn) && !is.numeric(adlb$avisitn))
    cdisc_issues <- c(cdisc_issues, "AVISITN is not numeric")

  n_no_avisitn <- sum(is.na(adlb$avisitn))
  if (n_no_avisitn > 0)
    message("NOTE: ", n_no_avisitn, " rows (", round(100 * n_no_avisitn / nrow(adlb), 1),
            "%) have no AVISITN (unmatched visit label)")

  if (length(cdisc_issues) == 0) {
    message("CDISC ADaM BDS validation OK")
  } else {
    for (issue in cdisc_issues) message("CDISC issue: ", issue)
  }
}

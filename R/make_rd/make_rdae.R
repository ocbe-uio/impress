# Reporting for adverse events (SAP §7).
#
# Treatment-emergent AE incidence on the configured cohort (cohortcd == 2),
# descriptive only, split into two phases by AE onset study day (ASTDY):
#   * Early  (onset Day <= 44): grouped by the stepped-wedge dose in effect at
#            onset (TRTCD), so a "0 mg" column captures AEs occurring before a
#            subject's randomised step has started. Dose columns are NOT mutually
#            exclusive (a subject can contribute to 0 mg and to their active dose).
#            Denominator per column = subjects ever at that dose during Days 1-44.
#   * Late   (onset Day >= 45): grouped by the long-term/maintenance planned dose
#            (TRT02P). Denominator = subjects per long-term dose in the analysis set.
#
# Mock-up grouping is by planned dose (FAS; _EX.csv empty), see make_adae.R.
#
# Returns list(early = <tables>, late = <tables>); each <tables> is a list of
# knitr::kable objects (summary, severity, relationship, soc_pt) plus nsubj, or
# NULL when a phase has no events.

# Build the four AE tables for one phase, given a TEAE set carrying a `grp` factor
# column, the per-group denominators (named integer vector over grp_levels), the
# column levels, and the total-column denominator (distinct subjects in the analysis
# set; for the early phase this differs from sum(Nvec) because dose columns overlap).
build_ae_tables <- function(teae, Nvec, grp_levels, Ntot) {
  if (is.null(teae) || nrow(teae) == 0) {
    return(NULL)
  }

  Nvec <- Nvec[grp_levels]
  Nvec[is.na(Nvec)] <- 0L
  if (is.null(Ntot) || is.na(Ntot)) Ntot <- sum(Nvec)
  if (Ntot == 0) {
    return(NULL)
  }

  fmt_cell <- function(n, N) ifelse(N > 0, sprintf("%d (%.1f%%)", n, 100 * n / N), "0")

  # n (%) of distinct subjects per group + Total, for a filtered AE set
  subj_counts <- function(df) {
    cc <- df %>%
      dplyr::distinct(grp, usubjid) %>%
      dplyr::count(grp, .drop = FALSE, name = "n")
    n <- stats::setNames(cc$n, as.character(cc$grp))[grp_levels]
    n[is.na(n)] <- 0L
    tot <- dplyr::n_distinct(df$usubjid)
    cells <- mapply(fmt_cell, n, Nvec, USE.NAMES = FALSE)
    stats::setNames(c(cells, fmt_cell(tot, Ntot)), c(grp_levels, "Total"))
  }

  row_tbl <- function(label, df) {
    tibble::tibble(Category = label) %>%
      dplyr::bind_cols(tibble::as_tibble_row(subj_counts(df)))
  }

  header_row <- function() {
    cells <- stats::setNames(
      c(sprintf("N=%d", Nvec), sprintf("N=%d", Ntot)),
      c(grp_levels, "Total")
    )
    tibble::tibble(Category = "Subjects in analysis set") %>%
      dplyr::bind_cols(tibble::as_tibble_row(cells))
  }

  # ---- Overall AE summary ---------------------------------------------------
  summary_tbl <- dplyr::bind_rows(
    header_row(),
    row_tbl("Subjects with >=1 TEAE", teae),
    row_tbl("  with >=1 serious TEAE", dplyr::filter(teae, aesdfl == "Y")),
    row_tbl("  with >=1 treatment-related TEAE", dplyr::filter(teae, arelfl == "Y")),
    row_tbl("  with >=1 TEAE leading to discontinuation", dplyr::filter(teae, adisconfl == "Y")),
    row_tbl("  with >=1 fatal TEAE", dplyr::filter(teae, adthfl == "Y"))
  ) %>%
    knitr::kable()

  # ---- Maximum severity per subject -----------------------------------------
  sev_levels <- c("Mild", "Moderate", "Severe", "Life-threatening", "Death")
  max_sev <- teae %>%
    dplyr::filter(!is.na(asevn)) %>%
    dplyr::group_by(grp, usubjid) %>%
    dplyr::summarise(msevn = max(asevn, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(msev = factor(sev_levels[msevn], levels = sev_levels))

  severity_tbl <- if (nrow(max_sev) == 0) NULL else {
    purrr::map_dfr(sev_levels, function(s) {
      row_tbl(s, max_sev %>% dplyr::filter(msev == s))
    }) %>%
      knitr::kable()
  }

  # ---- Relationship to treatment (subjects with >=1 TEAE in category) -------
  rel_levels <- teae %>%
    dplyr::filter(!is.na(arel)) %>%
    dplyr::distinct(arel) %>%
    dplyr::pull(arel)
  relationship_tbl <- if (length(rel_levels) == 0) NULL else {
    purrr::map_dfr(rel_levels, function(r) {
      row_tbl(r, dplyr::filter(teae, arel == r))
    }) %>%
      knitr::kable()
  }

  # ---- SOC x PT incidence ---------------------------------------------------
  socs <- teae %>%
    dplyr::distinct(soc_name) %>%
    dplyr::arrange(soc_name) %>%
    dplyr::pull(soc_name)

  socpt_rows <- purrr::map_dfr(socs, function(soc) {
    soc_df <- teae %>% dplyr::filter(soc_name == soc)
    soc_row <- tibble::tibble(`SOC / Preferred term` = soc) %>%
      dplyr::bind_cols(tibble::as_tibble_row(subj_counts(soc_df)))
    pts <- soc_df %>% dplyr::distinct(pt_name) %>% dplyr::arrange(pt_name) %>% dplyr::pull(pt_name)
    pt_rows <- purrr::map_dfr(pts, function(pt) {
      tibble::tibble(`SOC / Preferred term` = paste0("  ", pt)) %>%
        dplyr::bind_cols(tibble::as_tibble_row(subj_counts(soc_df %>% dplyr::filter(pt_name == pt))))
    })
    dplyr::bind_rows(soc_row, pt_rows)
  })
  socpt_tbl <- socpt_rows %>% knitr::kable()

  list(
    nsubj        = dplyr::n_distinct(teae$usubjid),
    summary      = summary_tbl,
    severity     = severity_tbl,
    relationship = relationship_tbl,
    soc_pt       = socpt_tbl
  )
}

make_rdae <- function(adae, adsl, cfg) {
  if (is.null(adae) || nrow(adae) == 0) {
    return(NULL)
  }

  trt_map <- trt_map_from_yaml(cfg)
  grp_levels <- c("losartan 0 mg", "losartan 25 mg", "losartan 50 mg", "losartan 100 mg")

  teae <- adae %>%
    filter_cohort(cfg) %>%
    dplyr::filter(trtemfl == "Y")
  if (nrow(teae) == 0) {
    return(NULL)
  }

  # ---- Late denominators: subjects per long-term planned dose ---------------
  late_denom <- adsl %>%
    filter_cohort(cfg) %>%
    dplyr::mutate(grp = factor(as.character(trt02p), levels = grp_levels)) %>%
    dplyr::filter(!is.na(grp)) %>%
    dplyr::count(grp, .drop = FALSE, name = "N")
  late_Nvec <- stats::setNames(late_denom$N, as.character(late_denom$grp))

  # ---- Early denominators: subjects ever at each stepped-wedge dose, Days 1-44
  early_long <- adsl %>%
    filter_cohort(cfg) %>%
    dplyr::transmute(subjid, armcd = as.character(armcd)) %>%
    dplyr::filter(!is.na(armcd)) %>%
    tidyr::crossing(day = 1:44) %>%
    dplyr::mutate(d = trt_from_map(armcd, day, trt_map)) %>%
    dplyr::filter(!is.na(d)) %>%
    dplyr::distinct(subjid, d)
  early_denom <- early_long %>%
    dplyr::mutate(grp = factor(paste0("losartan ", d, " mg"), levels = grp_levels)) %>%
    dplyr::count(grp, .drop = FALSE, name = "N")
  early_Nvec <- stats::setNames(early_denom$N, as.character(early_denom$grp))
  early_Ntot <- dplyr::n_distinct(early_long$subjid)
  late_Ntot  <- sum(late_Nvec, na.rm = TRUE)

  # ---- Phase datasets with grouping factor ----------------------------------
  early <- teae %>%
    dplyr::filter(!is.na(astdy), astdy <= 44, !is.na(trtcd)) %>%
    dplyr::mutate(grp = factor(paste0("losartan ", as.character(trtcd), " mg"),
                               levels = grp_levels)) %>%
    dplyr::filter(!is.na(grp))

  late <- teae %>%
    dplyr::filter(!is.na(astdy), astdy >= 45) %>%
    dplyr::mutate(grp = factor(as.character(trt02p), levels = grp_levels)) %>%
    dplyr::filter(!is.na(grp))

  list(
    early = build_ae_tables(early, early_Nvec, grp_levels, early_Ntot),
    late  = build_ae_tables(late,  late_Nvec,  grp_levels, late_Ntot)
  )
}

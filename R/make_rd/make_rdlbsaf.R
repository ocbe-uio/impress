# Reporting for conventional safety laboratory parameters (SAP §7).
#
# Same dual approach as vital signs (descriptive + repeated-measures LMM, reusing
# the safety_cont_* helpers from make_rdvs.R), plus categorisation against reference
# ranges with counts/proportions outside the normal range.
#
# Mock-up: creatinine only (the only populated LBE block in this export).

# Counts / proportions LOW / NORMAL / HIGH vs reference range, by dose group.
safety_lab_refrange <- function(data, var, cfg) {
  dt <- data %>%
    filter_cohort(cfg) %>%
    dplyr::filter(paramcd == var, avisitn >= -7, avisitn < 999) %>%
    dplyr::filter(!is.na(anrind), !is.na(trtcd))

  if (nrow(dt) == 0) {
    return(NULL)
  }

  dt %>%
    dplyr::mutate(trt = factor(paste0(trtcd, " mg"),
                               levels = c("0 mg", "25 mg", "50 mg", "100 mg"))) %>%
    dplyr::group_by(trt, .drop = TRUE) %>%
    dplyr::mutate(N = dplyr::n()) %>%
    dplyr::group_by(trt, anrind, N, .drop = FALSE) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(N > 0) %>%
    dplyr::mutate(txt = sprintf("%d (%.1f%%)", n, 100 * n / N)) %>%
    dplyr::select(trt, anrind, txt) %>%
    tidyr::pivot_wider(names_from = anrind, values_from = txt, values_fill = "0 (0.0%)") %>%
    knitr::kable(caption = NULL)
}

make_rdlbsaf <- function(adlbsaf, cfg) {
  if (is.null(adlbsaf) || nrow(adlbsaf) == 0) {
    return(NULL)
  }

  vars <- adlbsaf %>%
    dplyr::distinct(paramcd) %>%
    dplyr::pull(paramcd) %>%
    as.character()

  purrr::set_names(vars) %>%
    purrr::map(function(v) {
      sec <- make_rdvs_section(adlbsaf, v, cfg)
      refrange <- safety_lab_refrange(adlbsaf, v, cfg)
      if (is.null(sec) && is.null(refrange)) {
        return(NULL)
      }
      if (is.null(sec)) {
        label <- adlbsaf %>%
          dplyr::filter(paramcd == v) %>%
          dplyr::pull(param) %>% unique() %>% as.character()
        sec <- list(descriptive = NULL, early = NULL, late = NULL,
                    label = if (length(label)) label[[1]] else v, nsubj = 0)
      }
      sec$refrange <- refrange
      sec
    })
}

# Reporting for vital signs (SAP §7): descriptive + repeated-measures LMM.
#
# Mirrors the efficacy continuous-endpoint pattern in make_rdeff.R
# (summarise_cont_early / summarise_cont_lmm). The shared helpers defined here
# (safety_cont_*) are also reused by make_rdlbsaf.R for conventional labs.
#
# All analyses are descriptive/exploratory on the configured cohort (cohortcd == 2),
# grouped by planned dose (FAS mock-up; see make_adae.R / make_adex.R).

# Descriptive: mean (SD) of observed value and change from baseline, by dose group
# (stepped-wedge dose in effect) and visit.
safety_cont_descriptive <- function(data, var, cfg) {
  dt <- data %>%
    filter_cohort(cfg) %>%
    dplyr::filter(paramcd == var, avisitn >= -7, avisitn < 999) %>%
    dplyr::filter(!is.na(aval), !is.na(trtcd))

  if (nrow(dt) == 0) {
    return(NULL)
  }

  tbl <- dt %>%
    dplyr::mutate(trt = factor(paste0(trtcd, " mg"),
                               levels = c("0 mg", "25 mg", "50 mg", "100 mg"))) %>%
    dplyr::group_by(avisitn, trt, .drop = TRUE) %>%
    dplyr::summarise(
      n        = dplyr::n(),
      obs_txt  = sprintf("%.1f (%.1f)", mean(aval, na.rm = TRUE), stats::sd(aval, na.rm = TRUE)),
      chg_txt  = {
        cv <- chg[!is.na(chg)]
        if (length(cv) == 0) "-" else sprintf("%.1f (%.1f)", mean(cv), stats::sd(cv))
      },
      .groups = "drop"
    ) %>%
    dplyr::arrange(avisitn, trt) %>%
    knitr::kable(
      col.names = c("Visit (day)", "Dose", "N", "Observed mean (SD)", "Change mean (SD)")
    )

  tbl
}

# Early window (Days 15/29/43): change from baseline, dose as the stepped-wedge
# treatment in effect.  Model: chg ~ base + trtcd + avisitn + (1 | usubjid).
safety_cont_early <- function(data, var, cfg) {
  dt <- data %>%
    filter_cohort(cfg) %>%
    dplyr::filter(paramcd == var, avisitn %in% c(15, 29, 43)) %>%
    dplyr::filter(ablfl != "Y", !is.na(chg), !is.na(base), !is.na(trtcd)) %>%
    dplyr::mutate(trtcd = droplevels(factor(trtcd, levels = c(0, 25, 50, 100))))

  if (nrow(dt) == 0 || dplyr::n_distinct(dt$trtcd) < 2 ||
      dplyr::n_distinct(dt$usubjid) >= nrow(dt)) {
    return(NULL)
  }

  m <- lme4::lmer(chg ~ base + trtcd + avisitn + (1 | usubjid), data = dt)
  emm <- emmeans::emmeans(m, specs = "trtcd")

  means <- broom::tidy(emm, conf.int = TRUE) %>%
    dplyr::mutate(
      trt = paste0(trtcd, " mg"),
      mean_txt = format_estimate_ci(estimate, conf.low, conf.high, digits = 2)
    )

  meanstbl <- means %>%
    dplyr::select(trt, mean_txt) %>%
    knitr::kable(col.names = c("Dose (mg)", "Estimated change from baseline (95% CI)"), digits = 2)

  diffstbl <- emmeans::contrast(emm, method = "trt.vs.ctrl") %>%
    broom::tidy(conf.int = TRUE) %>%
    dplyr::mutate(
      contrast = gsub("^trtcd", "losartan ", contrast),
      contrast = gsub(" - trtcd", " mg vs losartan 0 mg", contrast),
      estimate_txt = format_estimate_ci(estimate, conf.low, conf.high, digits = 2)
    ) %>%
    dplyr::select(contrast, estimate_txt) %>%
    knitr::kable(col.names = c("Comparison", "Difference (95% CI)"), digits = 2)

  plt <- means %>%
    dplyr::mutate(trtcd_num = suppressWarnings(as.numeric(as.character(trtcd)))) %>%
    ggplot2::ggplot(ggplot2::aes(x = trtcd_num, y = estimate)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = conf.low, ymax = conf.high), width = 1) +
    ggplot2::labs(x = "Dose (mg)", y = "Estimated change from baseline") +
    ggplot2::theme_minimal()

  list(means = meanstbl, diffs = diffstbl, plot = plt,
       nsubj = dplyr::n_distinct(dt$usubjid))
}

# Late window (Days 155/239): change from baseline by long-term planned dose.
# Model: chg ~ base + trt02p * avisitn + (1 | usubjid).
safety_cont_late <- function(data, var, cfg) {
  dt <- data %>%
    filter_cohort(cfg) %>%
    dplyr::filter(paramcd == var, avisitn %in% c(155, 239)) %>%
    dplyr::filter(ablfl != "Y", !is.na(chg), !is.na(base), !is.na(trt02p)) %>%
    dplyr::mutate(avisitn = factor(avisitn), trt02p = droplevels(factor(trt02p)))

  trt02p_levels <- dt %>%
    dplyr::distinct(dose02p, trt02p) %>%
    dplyr::arrange(dose02p, trt02p) %>%
    dplyr::pull(trt02p) %>%
    purrr::discard(is.na)

  if (length(trt02p_levels)) {
    dt <- dt %>% dplyr::mutate(trt02p = factor(trt02p, levels = trt02p_levels))
  }

  if (nrow(dt) == 0 || dplyr::n_distinct(dt$trt02p) < 2 ||
      dplyr::n_distinct(dt$usubjid) >= nrow(dt)) {
    return(NULL)
  }

  target_visits <- intersect(c("155", "239"), levels(dt$avisitn))
  if (length(target_visits) == 0) {
    return(NULL)
  }

  m <- lme4::lmer(chg ~ base + trt02p * avisitn + (1 | usubjid), data = dt)
  emm <- emmeans::emmeans(m, specs = c("trt02p", "avisitn"), at = list(avisitn = target_visits))

  means <- broom::tidy(emm, conf.int = TRUE) %>%
    dplyr::mutate(
      trt02p = factor(trt02p, levels = trt02p_levels),
      mean_txt = format_estimate_ci(estimate, conf.low, conf.high, digits = 2)
    )

  meanstbl <- means %>%
    dplyr::select(avisitn, trt02p, mean_txt) %>%
    knitr::kable(col.names = c("Visit", "Treatment", "Estimated change from baseline (95% CI)"), digits = 2)

  diffstbl <- emmeans::contrast(emm, method = "trt.vs.ctrl", by = "avisitn") %>%
    broom::tidy(conf.int = TRUE) %>%
    dplyr::mutate(estimate_txt = format_estimate_ci(estimate, conf.low, conf.high, digits = 2)) %>%
    dplyr::select(avisitn, contrast, estimate_txt) %>%
    knitr::kable(col.names = c("Visit", "Comparison", "Difference (95% CI)"), digits = 2)

  dodge <- ggplot2::position_dodge(width = get_plot_dodge_width(cfg))
  plt <- means %>%
    ggplot2::ggplot(ggplot2::aes(x = avisitn, y = estimate, color = trt02p, group = trt02p)) +
    ggplot2::geom_point(position = dodge) +
    ggplot2::geom_line(position = dodge) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = conf.low, ymax = conf.high), width = 0.15, position = dodge) +
    ggplot2::labs(x = "Visit", y = "Estimated change from baseline", color = "Long term treatment") +
    ggplot2::scale_color_discrete(limits = trt02p_levels) +
    ggplot2::theme_minimal()

  list(means = meanstbl, diffs = diffstbl, plot = plt,
       nsubj = dplyr::n_distinct(dt$usubjid))
}

make_rdvs_section <- function(data, var, cfg) {
  descriptive <- safety_cont_descriptive(data, var, cfg)
  early <- safety_cont_early(data, var, cfg)
  late  <- safety_cont_late(data, var, cfg)

  if (is.null(descriptive) && is.null(early) && is.null(late)) {
    return(NULL)
  }

  label <- data %>%
    dplyr::filter(paramcd == var) %>%
    dplyr::pull(param) %>%
    unique() %>%
    as.character()
  label <- if (length(label)) label[[1]] else as.character(var)

  nsubj <- max(early$nsubj %||% 0, late$nsubj %||% 0, na.rm = TRUE)

  list(descriptive = descriptive, early = early, late = late,
       label = label, nsubj = nsubj)
}

make_rdvs_batch <- function(advs, vs_vars, cfg) {
  purrr::set_names(vs_vars) %>%
    purrr::map(~ make_rdvs_section(advs, .x, cfg))
}

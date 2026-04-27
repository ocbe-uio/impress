# make_rdlb.R — Reporting / analysis functions for ADLB (biomarker) endpoints
#
# Analytes are serum cytokines and angiogenic factors measured in pg/mL from
# Bio-Plex 19-plex and 14-plex panels.  aval is the raw concentration; chg is
# change from the baseline sample (adt <= rfstdt); ablfl == "Y" flags baseline
# rows and ablfl == "N" flags all post-treatment rows.
#
# Two analysis windows:
#   Early  (avisitn %in% c(15, 43))   — stepped-wedge dose-response via LMM
#   Late   (avisitn >= 197, Cycle 15+) — long-term trt02p analysis via LM
#
# For each analyte a linear model and a log-transformed model are fitted; AIC
# selects the preferred scale in the early analysis.  The late analysis uses
# the same preferred scale to keep results comparable.

# ---- helper: per-analyte offset for log transform -------------------------
# offset = half the smallest positive observed value, so that zeros can be
# log-transformed without producing -Inf.  Returns 0 if all values are NA.
lb_log_offset <- function(x) {
  pos <- x[!is.na(x) & x > 0]
  if (length(pos) == 0L) return(0)
  min(pos) / 2
}


# ---- analysis: stepped-wedge dose-response at Days 15 and 43 -------------
summarise_lb_early <- function(data, var, cfg) {
  dt <- data %>%
    filter_cohort(cfg) %>%
    dplyr::filter(paramcd == var, avisitn %in% c(15L, 43L)) %>%
    dplyr::filter(ablfl != "Y") %>%
    dplyr::filter(!is.na(chg), !is.na(base), !is.na(trtcd))

  if (nrow(dt) == 0 || dplyr::n_distinct(dt$trtcd) < 2) {
    return(NULL)
  }

  label <- as.character(unique(dt$param))
  nsubj <- dplyr::n_distinct(dt$usubjid)

  # ---- linear model on change from baseline --------------------------------
  m_lin <- tryCatch(
    lme4::lmer(chg ~ base + trtcd + avisitn + (1 | usubjid), data = dt),
    error = function(e) NULL
  )
  if (is.null(m_lin)) return(NULL)

  emm_lin   <- emmeans::emmeans(m_lin, specs = "trtcd")
  means_lin <- broom::tidy(emm_lin, conf.int = TRUE) %>%
    dplyr::mutate(
      trt      = paste0(trtcd, " mg"),
      mean_txt = format_estimate_ci(estimate, conf.low, conf.high, digits = 2)
    )

  # ---- log-transformed model: log(aval) - log(base) ~ log(base) + trtcd + avisitn
  # offset handles occasional zero-concentration samples
  offset <- lb_log_offset(dt$aval)
  dt_log <- dt %>%
    dplyr::filter(!is.na(aval), !is.na(base)) %>%
    dplyr::mutate(
      log_chg  = log(aval  + offset) - log(base + offset),
      log_base = log(base  + offset)
    )

  m_log     <- NULL
  emm_log   <- NULL
  means_log <- NULL
  if (nrow(dt_log) >= 10 && dplyr::n_distinct(dt_log$trtcd) >= 2) {
    m_log <- tryCatch(
      lme4::lmer(log_chg ~ log_base + trtcd + avisitn + (1 | usubjid),
                 data = dt_log),
      error = function(e) NULL
    )
    if (!is.null(m_log)) {
      emm_log   <- emmeans::emmeans(m_log, specs = "trtcd")
      means_log <- broom::tidy(emm_log, conf.int = TRUE) %>%
        dplyr::mutate(
          trt      = paste0(trtcd, " mg"),
          mean_txt = format_estimate_ci(estimate, conf.low, conf.high, digits = 3)
        )
    }
  }

  # ---- preferred model by AIC ---------------------------------------------
  aic_lin <- AIC(m_lin)
  aic_log <- if (!is.null(m_log)) AIC(m_log) else NA_real_
  preferred_scale <- if (!is.na(aic_log) && aic_log < aic_lin) "log" else "linear"
  preferred_emm   <- if (preferred_scale == "log") emm_log else emm_lin
  preferred_means <- if (preferred_scale == "log") means_log else means_lin
  preferred_model <- if (preferred_scale == "log") m_log else m_lin

  y_label <- if (preferred_scale == "log") {
    "Change from baseline on the log scale (95% CI)"
  } else {
    "Estimated change from baseline, pg/mL (95% CI)"
  }

  col_label <- if (preferred_scale == "log") {
    "Change from baseline on the natural logarithmic scale (95% CI)"
  } else {
    "Change from baseline, pg/mL (95% CI)"
  }

  diff_col_label <- if (preferred_scale == "log") {
    "Difference of the log-transformed values (95% CI)"
  } else {
    "Difference (95% CI)"
  }

  # ---- tables -------------------------------------------------------------
  meanstbl <- preferred_means %>%
    dplyr::select(trt, mean_txt) %>%
    knitr::kable(col.names = c("Dose (mg)", col_label), digits = 2)

  diffstbl <- emmeans::contrast(preferred_emm, method = "trt.vs.ctrl") %>%
    broom::tidy(conf.int = TRUE) %>%
    dplyr::mutate(
      contrast     = gsub("^trtcd", "losartan ", contrast),
      contrast     = gsub(" - trtcd", " mg vs losartan 0 mg", contrast),
      estimate_txt = format_estimate_ci(estimate, conf.low, conf.high, digits = 2)
    ) %>%
    dplyr::select(contrast, estimate_txt) %>%
    knitr::kable(col.names = c("Comparison", diff_col_label), digits = 2)

  aic_tbl <- tibble::tibble(
    Model     = c("Linear (chg ~ base + trtcd + avisitn)",
                  "Log-transformed (log(value) - log(base) ~ log(base) + trtcd + avisitn)"),
    AIC       = c(round(aic_lin, 1), round(aic_log, 1)),
    Preferred = c(preferred_scale == "linear", preferred_scale == "log")
  ) %>%
    knitr::kable()

  # ---- plot ---------------------------------------------------------------
  plt <- preferred_means %>%
    dplyr::mutate(trtcd_num = suppressWarnings(as.numeric(as.character(trtcd)))) %>%
    ggplot2::ggplot(ggplot2::aes(x = trtcd_num, y = estimate)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = conf.low, ymax = conf.high), width = 1) +
    ggplot2::labs(x = "Dose (mg)", y = y_label,
                  title = label) +
    ggplot2::theme_minimal()

  list(
    means           = meanstbl,
    diffs           = diffstbl,
    aic             = aic_tbl,
    plot            = plt,
    label           = label,
    nsubj           = nsubj,
    preferred_scale = preferred_scale,
    model           = preferred_model
  )
}


# ---- batch over all paramcds in adlb: early analysis -----------------------
#' @param data adlb dataset (output of make_adlb)
#' @param vars named character vector of paramcds (e.g. from rlang::set_names)
#' @param cfg  config list
make_rdlb_batch <- function(data, vars, cfg) {
  data <- filter_cohort(data, cfg)
  purrr::imap(vars, ~ summarise_lb_early(data, var = .x, cfg = cfg))
}


# ---- analysis: long-term assessment at Cycle 15+ ---------------------------
# Uses the same scale (linear or log) as the corresponding early analysis so
# that results are comparable.  A simple LM is used instead of LMM because
# most participants have only one late sample (Cycle 17).
# Treatment variable is trt02p (long-term assigned dose after Cycle 4).
summarise_lb_late <- function(data, var, cfg, scale = "linear") {
  dt <- data %>%
    filter_cohort(cfg) %>%
    dplyr::filter(paramcd == var, avisitn >= 197L) %>%
    dplyr::filter(ablfl != "Y") %>%
    dplyr::filter(!is.na(chg), !is.na(base), !is.na(trt02p))

  if (nrow(dt) == 0 || dplyr::n_distinct(dt$trt02p) < 2) {
    return(NULL)
  }

  label <- as.character(unique(dt$param))
  nsubj <- dplyr::n_distinct(dt$usubjid)

  # ---- linear model on change from baseline ---------------------------------
  m_lin <- tryCatch(
    lm(chg ~ base + trt02p, data = dt),
    error = function(e) NULL
  )
  if (is.null(m_lin)) return(NULL)

  emm_lin   <- emmeans::emmeans(m_lin, specs = "trt02p")
  means_lin <- broom::tidy(emm_lin, conf.int = TRUE) %>%
    dplyr::mutate(mean_txt = format_estimate_ci(estimate, conf.low, conf.high, digits = 2))

  # ---- log-transformed model ------------------------------------------------
  offset <- lb_log_offset(dt$aval)
  dt_log <- dt %>%
    dplyr::filter(!is.na(aval)) %>%
    dplyr::mutate(
      log_chg  = log(aval + offset) - log(base + offset),
      log_base = log(base + offset)
    )

  m_log     <- NULL
  emm_log   <- NULL
  means_log <- NULL
  if (nrow(dt_log) >= 10 && dplyr::n_distinct(dt_log$trt02p) >= 2) {
    m_log <- tryCatch(
      lm(log_chg ~ log_base + trt02p, data = dt_log),
      error = function(e) NULL
    )
    if (!is.null(m_log)) {
      emm_log   <- emmeans::emmeans(m_log, specs = "trt02p")
      means_log <- broom::tidy(emm_log, conf.int = TRUE) %>%
        dplyr::mutate(mean_txt = format_estimate_ci(estimate, conf.low, conf.high, digits = 3))
    }
  }

  # ---- AIC table (for reference; scale is imposed from early analysis) ------
  aic_lin <- AIC(m_lin)
  aic_log <- if (!is.null(m_log)) AIC(m_log) else NA_real_
  aic_tbl <- tibble::tibble(
    Model = c(
      "Linear (chg ~ base + trt02p)",
      "Log-transformed (log(value) - log(base) ~ log(base) + trt02p)"
    ),
    AIC  = c(round(aic_lin, 1), round(aic_log, 1)),
    Used = c(scale == "linear", scale == "log")
  ) %>%
    knitr::kable()

  # Use scale imposed by early analysis; fall back to linear if log failed
  use_scale       <- if (scale == "log" && !is.null(emm_log)) "log" else "linear"
  preferred_emm   <- if (use_scale == "log") emm_log else emm_lin
  preferred_means <- if (use_scale == "log") means_log else means_lin

  col_label <- if (use_scale == "log") {
    "Change from baseline on the natural logarithmic scale (95% CI)"
  } else {
    "Change from baseline, pg/mL (95% CI)"
  }
  diff_col_label <- if (use_scale == "log") {
    "Difference of the log-transformed values (95% CI)"
  } else {
    "Difference (95% CI)"
  }
  y_label <- if (use_scale == "log") {
    "Change from baseline on the log scale (95% CI)"
  } else {
    "Estimated change from baseline, pg/mL (95% CI)"
  }

  # ---- tables ---------------------------------------------------------------
  meanstbl <- preferred_means %>%
    dplyr::select(trt02p, mean_txt) %>%
    knitr::kable(col.names = c("Long-term dose", col_label), digits = 2)

  diffstbl <- emmeans::contrast(preferred_emm, method = "trt.vs.ctrl") %>%
    broom::tidy(conf.int = TRUE) %>%
    dplyr::mutate(
      estimate_txt = format_estimate_ci(estimate, conf.low, conf.high, digits = 2)
    ) %>%
    dplyr::select(contrast, estimate_txt) %>%
    knitr::kable(col.names = c("Comparison", diff_col_label), digits = 2)

  # ---- plot -----------------------------------------------------------------
  plt <- preferred_means %>%
    dplyr::mutate(dose_num = readr::parse_number(as.character(trt02p))) %>%
    ggplot2::ggplot(ggplot2::aes(x = dose_num, y = estimate)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = conf.low, ymax = conf.high), width = 1) +
    ggplot2::labs(x = "Long-term dose (mg)", y = y_label, title = label) +
    ggplot2::theme_minimal()

  list(
    means           = meanstbl,
    diffs           = diffstbl,
    aic             = aic_tbl,
    plot            = plt,
    label           = label,
    nsubj           = nsubj,
    preferred_scale = use_scale
  )
}


# ---- batch over all paramcds: late analysis (uses early preferred scale) ----
#' @param data          adlb dataset
#' @param vars          named character vector of paramcds
#' @param cfg           config list
#' @param early_summaries output of make_rdlb_batch (to extract preferred scale)
make_rdlb_late_batch <- function(data, vars, cfg, early_summaries) {
  data <- filter_cohort(data, cfg)
  purrr::imap(vars, function(var, nm) {
    scale <- if (!is.null(early_summaries[[nm]])) early_summaries[[nm]]$preferred_scale else "linear"
    summarise_lb_late(data, var = var, cfg = cfg, scale = scale)
  })
}

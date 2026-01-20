# library(dplyr)
# library(lme4)
# library(emmeans)
# library(ggplot2)
# library(broom)

cont_paramcds <- function() {
  c("kps", "ecog", "nano_tot", "qlq_c30", "qlq_bn20")
}

neuro_paramcds <- function() {
  c("kps", "ecog", "nano_tot")
}

qol_paramcds <- function() {
  c("qlq_c30", "qlq_bn20")
}

steroid_paramcd <- function() {
  "cs_steroids_mri"
}

rano_paramcds <- function() {
  c("trpspro", "trsdisea", "trorresp")
}
summarise_cont_early <- function(data, var, cfg) {
  dt <- data %>%
    filter_cohort(cfg) %>%
    filter(paramcd == var, avisitn %in% c(15, 29, 43)) %>%
    filter(ablfl != "Y") %>%
    filter(!is.na(chg)) %>%
    filter(!is.na(trtcd))

  if (nrow(dt) == 0 || dplyr::n_distinct(dt$trtcd) < 2) {
    return(NULL)
  }

  m <- lme4::lmer(chg ~ base + trtcd + avisitn + (1 | usubjid), data = dt)
  emm <- emmeans::emmeans(m, specs = "trtcd")

  means <- broom::tidy(emm, conf.int = TRUE) %>%
    mutate(
      trt = if (is.numeric(trtcd)) paste0(trtcd, " mg") else as.character(trtcd),
      mean_txt = format_estimate_ci(estimate, conf.low, conf.high, digits = 2)
    )

  meanstbl <- means %>%
    select(trt, mean_txt) %>%
    knitr::kable(col.names = c("Dose (mg)", "Estimated change from baseline (95% CI)"), digits = 2)

  diffstbl <- emmeans::contrast(emm, method = "trt.vs.ctrl") %>%
    broom::tidy(conf.int = TRUE) %>%
    mutate(
      contrast = gsub("^trtcd", "losartan ", contrast),
      contrast = gsub(" - trtcd", " mg vs losartan 0 mg", contrast),
      estimate_txt = format_estimate_ci(estimate, conf.low, conf.high, digits = 2)
    ) %>%
    select(contrast, estimate_txt) %>%
    knitr::kable(col.names = c("Comparison", "Difference (95% CI)"), digits = 2)

  plt <- means %>%
    mutate(trtcd_num = suppressWarnings(as.numeric(as.character(trtcd)))) %>%
    ggplot(aes(x = trtcd_num, y = estimate)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 1) +
    labs(x = "Dose (mg)", y = "Estimated change from baseline") +
    theme_minimal()

  label <- as.character(unique(dt$param))
  nsubj <- dt %>% distinct(usubjid) %>% nrow()

  plot_step <- NULL
  if ("step" %in% names(dt)) {
    m_step <- lme4::lmer(chg ~ base + step * avisitn + (1 | usubjid), data = dt)
    step_emm <- emmeans::emmeans(m_step, specs = c("avisitn", "step"))
    plot_df0 <- as.data.frame(step_emm) %>%
      mutate(
        avisitn = as.character(avisitn),
        step = as.character(step)
      )

    baseline_rows <- tibble::tibble(
      avisitn = "Baseline",
      step = levels(droplevels(dt$step)),
      emmean = 0,
      lower.CL = NA_real_,
      upper.CL = NA_real_
    )
    plot_df <- bind_rows(plot_df0, baseline_rows) %>%
      mutate(
        avisitn = factor(avisitn, levels = c("Baseline", as.character(unique(dt$avisitn)))),
        step = factor(step, levels = levels(dt$step)),
        step = droplevels(step)
      )

    dodge <- position_dodge(width = get_plot_dodge_width(cfg))
    plot_step <- ggplot(plot_df, aes(x = avisitn, y = emmean, color = step, group = step)) +
      geom_point(position = dodge) +
      geom_line(position = dodge) +
      geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.15, position = dodge) +
      labs(x = "Visit", y = "Estimated change from baseline", color = "Treatment start") +
      theme_minimal()
  }
  #list(plot_step = plot_step, label = label, nsubj = nsubj)  
  list(means = meanstbl, diffs = diffstbl, plot = plt, plot_step = plot_step, label = label, nsubj = nsubj)
}

summarise_cont_lmm <- function(data, var, cfg) {
  dt <- data %>%
    filter_cohort(cfg) %>%
    filter(paramcd == var, avisitn > -7) %>%
    filter(ablfl != "Y") %>%
    filter(!is.na(chg), !is.na(base), !is.na(trt02p)) %>%
    mutate(
      avisitn = factor(avisitn),
      trt02p = droplevels(factor(trt02p))
    )

  trt02p_levels <- dt %>%
    distinct(dose02p, trt02p) %>%
    arrange(dose02p, trt02p) %>%
    pull(trt02p) %>%
    purrr::discard(is.na)

  if (length(trt02p_levels)) {
    dt <- dt %>%
      mutate(trt02p = factor(trt02p, levels = trt02p_levels))
  }

  if (nrow(dt) == 0 || dplyr::n_distinct(dt$trt02p) < 2) {
    return(NULL)
  }

  if (dplyr::n_distinct(dt$usubjid) >= nrow(dt)) {
    return(NULL)
  }

  target_visits <- intersect(c("155", "239"), levels(dt$avisitn))
  if (length(target_visits) == 0) {
    return(NULL)
  }

  m <- lme4::lmer(chg ~ base + trt02p * avisitn + (1 | usubjid), data = dt)
  emm <- emmeans::emmeans(m, specs = c("trt02p", "avisitn"), at = list(avisitn = target_visits))

  means <- broom::tidy(emm, conf.int = TRUE) %>%
    mutate(trt02p = factor(trt02p, levels = trt02p_levels)) %>%
    mutate(
      mean_txt = format_estimate_ci(estimate, conf.low, conf.high, digits = 2)
    )

  meanstbl <- means %>%
    select(avisitn, trt02p, mean_txt) %>%
    knitr::kable(col.names = c("Visit", "Treatment", "Estimated change from baseline (95% CI)"), digits = 2)

  diffstbl <- emmeans::contrast(emm, method = "trt.vs.ctrl", by = "avisitn") %>%
    broom::tidy(conf.int = TRUE) %>%
    mutate(
      estimate_txt = format_estimate_ci(estimate, conf.low, conf.high, digits = 2)
    ) %>%
    select(avisitn, contrast, estimate_txt) %>%
    knitr::kable(col.names = c("Visit", "Comparison", "Difference (95% CI)"), digits = 2)

  dodge <- position_dodge(width = get_plot_dodge_width(cfg))
  plt <- means %>%
    ggplot(aes(x = avisitn, y = estimate, color = trt02p, group = trt02p)) +
    geom_point(position = dodge) +
    geom_line(position = dodge) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15, position = dodge) +
    labs(x = "Visit", y = "Estimated change from baseline", color = "Long term treatment") +
    scale_color_discrete(limits = trt02p_levels) +
    theme_minimal()

  label <- as.character(unique(dt$param))
  nsubj <- dt %>% distinct(usubjid) %>% nrow()

  list(means = meanstbl, diffs = diffstbl, plot = plt, label = label, nsubj = nsubj)
}

make_cont_section <- function(data, var, cfg) {
  early <- if (var %in% qol_paramcds()) NULL else summarise_cont_early(data, var, cfg)
  late <- summarise_cont_lmm(data, var, cfg)

  if (is.null(early) && is.null(late)) {
    return(NULL)
  }

  label <- unique(na.omit(c(early$label, late$label)))
  label <- if (length(label)) as.character(label[[1]]) else as.character(var)
  nsubj <- max(
    early$nsubj %||% 0,
    late$nsubj %||% 0,
    na.rm = TRUE
  )

  list(
    early = early,
    late = late,
    label = label,
    nsubj = nsubj
  )
}

summarise_steroid_early <- function(data, cfg) {
  dt <- data %>%
    filter_cohort(cfg) %>%
    filter(paramcd == steroid_paramcd(), avisitn %in% c(15, 29, 43)) %>%
    filter(ablfl != "Y") %>%
    filter(!is.na(trtcd), !is.na(step)) %>%
    mutate(
      avisitn = factor(avisitn),
      trtcd = droplevels(factor(trtcd, levels = c(0, 25, 50, 100))),
      step = droplevels(factor(step)),
      aval_num = binary_yesno(avalc)
    )

  if (nrow(dt) == 0 || dplyr::n_distinct(dt$trtcd) < 2) {
    return(NULL)
  }

  base_map <- steroid_baseline_map(data, cfg, steroid_paramcd())

  dt <- dt %>%
    left_join(base_map, by = "usubjid") %>%
    filter(!is.na(aval_num), !is.na(base_flag))

  if (nrow(dt) == 0) {
    return(NULL)
  }

  m <- lme4::glmer(aval_num ~ base_flag + trtcd + avisitn + (1 | usubjid),
                   data = dt, family = binomial())

  preds <- marginaleffects::avg_predictions(
    m,
    variables = "trtcd",
    type = "response",
    re.form = NA
  )

  pred_tbl <- preds %>%
    mutate(
      trt = paste0("losartan ", trtcd, " mg"),
      est_txt = format_estimate_ci(estimate, conf.low, conf.high, digits = 3)
    ) %>%
    select(trt, est_txt) %>%
    knitr::kable(col.names = c("Treatment", "Risk (95% CI)"), digits = 3)

  diffs <- marginaleffects::avg_comparisons(
    m,
    variables = "trtcd",
    type = "response",
    re.form = NA
  )

  difftbl <- diffs %>%
    mutate(
      contrast = gsub("^trtcd", "losartan ", contrast),
      contrast = gsub(" - trtcd", " mg vs losartan ", contrast),
      estimate_txt = format_estimate_ci(estimate, conf.low, conf.high, digits = 3)
    ) %>%
    select( contrast, estimate_txt) %>%
    knitr::kable(col.names = c("Comparison", "Risk difference (95% CI)"), digits = 3)

  plot_df <- preds %>%
    mutate(trtcd_num = suppressWarnings(as.numeric(as.character(trtcd))))

  plt <- plot_df %>%
    ggplot(aes(x = trtcd_num, y = estimate)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 1) +
    labs(x = "Dose (mg)", y = "Risk of steroid use") +
    theme_minimal()

  nsubj <- dt %>% distinct(usubjid) %>% nrow()

  list(means = pred_tbl, diffs = difftbl, plot = plt, nsubj = nsubj)
}

summarise_steroid_lmm <- function(data, cfg) {
  dt <- data %>%
    filter_cohort(cfg) %>%
    filter(paramcd == steroid_paramcd(), avisitn %in% c(155, 239)) %>%
    filter(ablfl != "Y") %>%
    filter(!is.na(trt02p)) %>%
    mutate(
      avisitn = factor(avisitn),
      trt02p = droplevels(factor(trt02p)),
      aval_num = binary_yesno(avalc)
    )

  trt02p_levels <- dt %>%
    distinct(dose02p, trt02p) %>%
    arrange(dose02p, trt02p) %>%
    pull(trt02p) %>%
    purrr::discard(is.na)

  if (length(trt02p_levels)) {
    dt <- dt %>%
      mutate(trt02p = factor(trt02p, levels = trt02p_levels))
  }

  if (nrow(dt) == 0 || dplyr::n_distinct(dt$trt02p) < 2) {
    return(NULL)
  }

  base_map <- steroid_baseline_map(data, cfg, steroid_paramcd())

  dt <- dt %>%
    left_join(base_map, by = "usubjid") %>%
    filter(!is.na(aval_num), !is.na(base_flag))

  if (nrow(dt) == 0) {
    return(NULL)
  }

  if (dplyr::n_distinct(dt$usubjid) >= nrow(dt)) {
    return(NULL)
  }

  target_visits <- intersect(c("155", "239"), levels(dt$avisitn))
  if (length(target_visits) == 0) {
    return(NULL)
  }

  m <- stats::glm(aval_num ~ base_flag + trt02p * avisitn,
                  data = dt, family = binomial())

  preds <- marginaleffects::avg_predictions(
    m,
    by = c("trt02p", "avisitn"),
    type = "response",
    vcov = ~ usubjid,
    newdata = marginaleffects::datagrid(
      trt02p = trt02p_levels,
      avisitn = target_visits,
      base_flag = mean(dt$base_flag, na.rm = TRUE)
    )
  )

  pred_tbl <- preds %>%
    mutate(
      est_txt = format_estimate_ci(estimate, conf.low, conf.high, digits = 3)
    ) %>%
    select(avisitn, trt02p, est_txt) %>%
    knitr::kable(col.names = c("Visit", "Treatment", "Risk (95% CI)"), digits = 3)

  diffs <- marginaleffects::avg_comparisons(
    m,
    variables = "trt02p",
    by = "avisitn",
    type = "response",
    vcov = ~ usubjid
  )

  difftbl <- diffs %>%
    mutate(
      estimate_txt = format_estimate_ci(estimate, conf.low, conf.high, digits = 3)
    ) %>%
    select(avisitn, contrast, estimate_txt) %>%
    knitr::kable(col.names = c("Visit", "Comparison", "Risk difference (95% CI)"), digits = 3)

  baseline_mean <- base_map %>%
    summarise(mean = mean(base_flag, na.rm = TRUE)) %>%
    pull(mean)

  plot_df <- preds %>%
    mutate(
      avisitn = as.character(avisitn),
      trt02p = as.character(trt02p)
    )

  if (!is.na(baseline_mean)) {
    baseline_rows <- tibble::tibble(
      avisitn = "Baseline",
      trt02p = trt02p_levels,
      estimate = baseline_mean,
      conf.low = NA_real_,
      conf.high = NA_real_
    )
    plot_df <- bind_rows(plot_df, baseline_rows)
  }

  plot_df <- plot_df %>%
    mutate(
      avisitn = factor(avisitn, levels = c("Baseline", levels(dt$avisitn))),
      trt02p = factor(trt02p, levels = trt02p_levels)
    )

  dodge <- position_dodge(width = get_plot_dodge_width(cfg))
  plt <- plot_df %>%
    ggplot(aes(x = avisitn, y = estimate, color = trt02p, group = trt02p)) +
    geom_point(position = dodge) +
    geom_line(position = dodge) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15, position = dodge) +
    labs(x = "Visit", y = "Risk of steroid use", color = "Long term treatment") +
    scale_color_discrete(limits = trt02p_levels) +
    theme_minimal()

  nsubj <- dt %>% distinct(usubjid) %>% nrow()

  list(means = pred_tbl, diffs = difftbl, plot = plt, nsubj = nsubj)
}

make_steroid_section <- function(data, cfg) {
  early <- summarise_steroid_early(data, cfg)
  late <- summarise_steroid_lmm(data, cfg)

  if (is.null(early) && is.null(late)) {
    return(NULL)
  }

  nsubj <- max(
    early$nsubj %||% 0,
    late$nsubj %||% 0,
    na.rm = TRUE
  )

  list(
    early = early,
    late = late,
    label = "Steroid use",
    nsubj = nsubj
  )
}

summarise_rano_polr <- function(data, var, cfg) {
  dt <- data %>%
    filter_cohort(cfg) %>%
    filter(paramcd == var) %>%
    filter(!is.na(avalc), !is.na(trt02p))

  visits <- if (var == "trorresp") 225 else c(141, 225)
  dt <- dt %>%
    filter(avisitn %in% visits)

  if (nrow(dt) == 0) {
    return(NULL)
  }

  rano_levels <- switch(
    var,
    trpspro = c("No", "Unknown", "Yes"),
    trsdisea = c("Progressive disease (PD)", "Stable disease (SD)", "Partial response (PR)", "Complete response (CR)"),
    trorresp = c("Progressive disease (PD)", "Stable disease (SD)", "Partial response (PR)", "Complete response (CR)"),
    unique(as.character(dt$avalc))
  )

  dt <- dt %>%
    mutate(
      trt02p = droplevels(factor(trt02p, levels = c("losartan 0 mg", "losartan 25 mg", "losartan 50 mg", "losartan 100 mg"))),
      avisitn = factor(avisitn),
      avalc_norm = dplyr::case_when(
        var == "trpspro" ~ stringr::str_to_sentence(tolower(as.character(avalc))),
        TRUE ~ as.character(avalc)
      ),
      rano = factor(avalc_norm, levels = rano_levels, ordered = TRUE)
    ) %>%
    filter(!is.na(rano))

  if (nrow(dt) == 0) {
    return(NULL)
  }

  visit_levels <- sort(unique(as.numeric(as.character(dt$avisitn))))
  visit_levels <- visit_levels[!is.na(visit_levels)]

  tables <- purrr::map(
    visit_levels,
    function(v) {
      dti <- dt %>%
        filter(as.numeric(as.character(avisitn)) == v)

      if (nrow(dti) == 0 || dplyr::n_distinct(dti$trt02p) < 2) {
        return(NULL)
      }

      if (dplyr::n_distinct(dti$rano) < 2) {
        return(NULL)
      }

      mod <- MASS::polr(rano ~ trt02p, data = dti, Hess = TRUE)
      tidy <- tryCatch(
        broom::tidy(mod, conf.int = TRUE, exponentiate = TRUE) |>
          dplyr::filter(coef.type =="coefficient") |> 
          dplyr::mutate(
            term = gsub("^trt02p", "", term),
            term = gsub("mg", " vs 0 mg", term),
            estimate_txt = format_estimate_ci(estimate, conf.low, conf.high, digits = 3)
          ) %>%
          dplyr::select(term, estimate_txt) %>%
          knitr::kable(col.names = c("Comparison", "Propoprtional Odds Ratio (95% CI)"), digits = 3),
        error = function(e) {
          broom::tidy(mod, conf.int = FALSE, exponentiate = TRUE) |>
            dplyr::filter(coef.type =="coefficient") |> 
            dplyr::mutate(
              term = gsub("^trt02p", "" , term),
              term = gsub("mg", "mg vs 0 mg", term),
              estimate_txt = paste0(format(round(estimate, 3), nsmall = 3), " (CI not available)")
            ) %>%
            dplyr::select(term, estimate_txt) %>%
            knitr::kable(col.names = c("Comparison", "Propoprtional Odds Ratio (95% CI)"), digits = 3)
        }
      )

      list(avisitn = v, table = tidy)
    }
  )

  tables <- purrr::compact(tables)
  if (!length(tables)) {
    return(NULL)
  }

  list(tables = tables, label = as.character(unique(dt$param)))
}




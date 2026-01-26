
# library(dplyr)
# library(lme4)
# library(ggplot2)
# library(purrr)

make_lineplot1 <- function(data, var, cfg) {

  dt0 <- data %>%
    filter(paramcd == var & avisitn <= 43)

  dt <- dt0 %>%
    filter(avisitn != -7) %>%
    mutate(
      vis = factor(
        avisitn,
        levels = c(15, 29, 43),
        labels = c("Day 15", "Day 29", "Day 43")
      ),
      step = droplevels(factor(step))
    )

  base_mean <- dt0 %>%
    filter(avisitn == -7) %>%
    summarise(mean = mean(base, na.rm = TRUE)) %>%
    pull(mean)

  m <- lme4::lmer(aval ~ base + step*vis + (1 | usubjid), data = dt)
  raw_est <- emmeans::emmeans(
    m,
    specs = c("vis", "step"),
    at = list(
      vis = levels(dt$vis),
      base = base_mean
    )
  )

  plot_df <- as.data.frame(raw_est) %>%
    mutate(
      vis = factor(vis, levels = levels(dt$vis))
    )

  base_rows <- tibble::tibble(
    vis = factor("Baseline", levels = c("Baseline", levels(dt$vis))),
    step = levels(dt$step),
    emmean = base_mean,
    lower.CL = NA_real_,
    upper.CL = NA_real_
  )

  plot_df <- plot_df %>%
    mutate(vis = factor(as.character(vis), levels = c("Baseline", levels(dt$vis)))) %>%
    bind_rows(base_rows)

  dodge <- position_dodge(width = get_plot_dodge_width(cfg))
  plot <- ggplot(plot_df, aes(x = vis, y = emmean, color = factor(step), group = factor(step))) +
    geom_point(position = dodge) +
    geom_line(position = dodge) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, position = dodge) +
    labs(
      x = "Visit",
      y = "Mean Response",
      color = "Treatment start"
    ) +
    theme_minimal() +
    ggtitle(paste("Mean Response over Time for", unique(dt$param)))
  
}

make_cmodels <- function(doses, cmodels) {
  fmodels <- Mods(
    linear = NULL, emax = cmodels$emax, sigEmax = cmodels$sigEmax, linInt = cmodels$linInt,
    doses = doses
  )
  fmodels_names <- names(fmodels)
  names(fmodels_names) <- names(fmodels)
  plot <- plot(fmodels, ylab = "Dose response", main = "Figure 1: Pre-specified candidate dose-response models")
  return(list(fmodels = fmodels, plot = plot))
}


# Make estimates for continuous longitudinal endpoints of the first three timepoints
make_contest <- function(data, var) {
  dt <- data %>%
    filter(paramcd == var & avisitn <= 43) |>
    filter(ablfl != "Y") 

  m <- lme4::lmer(chg ~ base + trtcd + avisitn + (1 | usubjid), data = dt)

  est_raw <- emmeans::emmeans(m, specs = "trtcd")
  S <- vcov(est_raw)

  est <- broom::tidy(est_raw, conf.int = TRUE) %>%
    mutate(
      trt = paste0(trtcd, " mg"),
      mean_txt = format_estimate_ci(estimate, conf.low, conf.high, digits = 2)
    )

table <- est %>%
  select(trt, mean_txt) |>
  knitr::kable(col.names = c("Treatment", "Estimated Change from Baseline (95% CI)"), digits = 2)

  list(model = m, estimates = est, S = S, emmeans = est_raw, table = table)
}

# Differences vs 0 mg reference
diff_vs_0mg <- function(est_obj) {
  if (is.null(est_obj$emmeans)) stop("est_obj must contain emmeans component")

  contr <- emmeans::contrast(est_obj$emmeans, method = "trt.vs.ctrl")
  broom::tidy(contr, conf.int = TRUE) %>%
    mutate(
      contrast = paste0(c(25, 50, 100), " mg vs 0 mg"),
      estimate_txt = paste0(
        round(estimate, 2),
        " (",
        round(conf.low, 2),
        " to ",
        round(conf.high, 2),
        ")"
      )
    )
}

any_dose_vs_0mg <- function(data, var) {
  dt <- data %>%
    filter(paramcd == var & avisitn <= 43) |>
    filter(!is.na(chg), !is.na(trtcd)) |>
    filter(ablfl != "Y") |>
    mutate(
      trt_num = suppressWarnings(readr::parse_number(as.character(trtcd))),
      trt_any = if_else(trt_num > 0, "Any dose", "0 mg"),
      trt_any = factor(trt_any, levels = c("0 mg", "Any dose"))
    )

  if (nrow(dt) == 0 || dplyr::n_distinct(dt$trt_any) < 2) {
    return(NULL)
  }

  m <- lme4::lmer(chg ~ base + trt_any + avisitn + (1 | usubjid), data = dt)
  emm <- emmeans::emmeans(m, specs = "trt_any")
  contr <- emmeans::contrast(emm, method = list("Any dose vs 0 mg" = c(-1, 1)))

  estimate <- broom::tidy(contr, conf.int = TRUE) %>%
    mutate(
      estimate_txt = format_estimate_ci(estimate, conf.low, conf.high, digits = 2)
    )

  table <- estimate %>%
    select(contrast, estimate_txt) %>%
    knitr::kable(col.names = c("Comparison", "Difference (95% CI)"), digits = 2)

  list(model = m, estimate = estimate, table = table)
}

f_fitmod <- function(estobj, doses, fmodels) {
  est_df <- estobj$estimates %>%
    mutate(trtcd = suppressWarnings(as.numeric(trtcd))) %>%
    arrange(match(trtcd, doses))

  cov_ord <- estobj$S[as.character(est_df$trtcd), as.character(est_df$trtcd), drop = FALSE]
  if (any(is.na(cov_ord))) stop("Covariance matrix could not be aligned to treatment levels.")

  model_names <- names(fmodels)
  fit_one <- function(model_name) {
    fitmod <- fitMod(doses, est_df$estimate, S = cov_ord, model = model_name, type = "general")
    plot <- plot(
      fitmod,
      CI = TRUE,
      plotData = "meansCI",
      ylab = "Estimated Change from baseline (95% CI)",
      xlab = "Dose (mg)"
    )

    fitparms <- c(
      AIC = gAIC(fitmod),
      ED50 = ED(fitmod, 0.5),
      ED70 = ED(fitmod, 0.7),
      ED90 = ED(fitmod, 0.9)
    )
    list(fitmod = fitmod, plot = plot, fitparms = fitparms)
  }

  fitmod <- purrr::set_names(lapply(model_names, fit_one), model_names)

  table <- purrr::imap_dfr(
    fitmod,
    ~ tibble::enframe(.x$fitparms, name = "name", value = "value") %>%
      mutate(Model = .y)
  ) %>%
    select(Model, name, value) %>%
    tidyr::pivot_wider(names_from = name, values_from = value) %>%
    arrange(AIC) %>%
    knitr::kable(col.names = c("Model", "AIC-value", "ED~50~", "ED~70~", "ED~90~"), digits = 2)

  list(fitmod = fitmod, table = table)
}

f_contmat <- function(estobj, fmodels) {

  contMat <- optContr(fmodels, S = estobj$S)
  table <- contMat$contMat %>%
    as_tibble(rownames = NA) %>%
    rownames_to_column(var="Treatment") %>%
    mutate(Treatment = paste0(Treatment, " mg")) %>%
    knitr::kable(digits = 3)

  list(contmat = contMat, table = table)
}


f_mcttest <- function(cfg, estobj, fmodels, contmat_obj){

  MCTtest <- MCTtest(cfg$doses, estobj$estimates$estimate,
                     S = estobj$S,
                     models = fmodels$fmodels,
                     type = "general",
                     contMat = contmat_obj$contMat,
                     pVal = TRUE, alternative = "one.sided",
                     alpha = 0.025,
                     critV = TRUE)
  p_table <- MCTtest$tStat %>%
    as_tibble(rownames=NA) %>%
    rownames_to_column(var="cm") %>%
    mutate(pval_ = attr(value,"pVal")) %>%
    mutate(pval = if_else(pval_<0.001, "<0.001",
                  if_else(pval_ <0.1, as.character(round(pval_,digits=3), digits=3),
                         as.character(round(pval_,digits=2),digits=2)))) %>%
    arrange(pval_) %>%
    select(-pval_) %>%
    knitr::kable(digits=2, align=c("lrr"),col.names = c("Candidate model", "t-statistic", "Adjusted p-value"))

  list(mcttest=MCTtest, p_table = p_table)

}

#' Bundle MCP-Mod outputs for R Markdown
#'
#' Given a dataset already filtered to the desired cohort/paramcd, compute
#' candidate models, estimates, contrasts, MCT tests, and fitted models. The
#' returned object contains plots/tables you can drop straight into an Rmd.
make_rdmri_section <- function(data, var, cfg) {
  data <- filter_cohort(data, cfg)
  lineplot <- make_lineplot1(data, var, cfg)
  cmodels  <- make_cmodels(cfg$doses, cfg$cmodels)
  est      <- make_contest(data, var)
  anydose  <- any_dose_vs_0mg(data, var)
  contmat  <- f_contmat(est, cmodels$fmodels)
  mct      <- f_mcttest(cfg, est, cmodels$fmodels, contmat$contmat)
  fitmods <- f_fitmod(est, cfg$doses, cmodels$fmodels)
  nsubj <- data %>%
    filter(paramcd == var) %>%
    distinct(usubjid) %>%
    nrow()

  list(
    text = list(
      intro = "The primary analysis follows the MCP-Mod methodology as detailed in the SAP, using pre-specified candidate dose-response models.",
      nsubj = paste0("A total of ", nsubj, " subjects were included in the analysis for this endpoint."),
      models = "The four candidate models (linear, piecewise linear, Emax, sigmoid Emax) are shown in Figure 1."
    ),
    plots = list(
      candidate_models = cmodels$plot,
      lineplot = lineplot,
      fit_models = purrr::map(fitmods$fitmod, "plot")
    ),
    tables = list(
      estimates = est$table,
      anydose = if (is.null(anydose)) NULL else anydose$table,
      contmat = contmat$table,
      mct = mct$p_table,
      fit = fitmods$table
    ),
    objects = list(
      est = est,
      anydose = anydose,
      cmodels = cmodels,
      contmat = contmat,
      mct = mct,
      fitmods = fitmods
    )
  )
}

# Batch MCP-Mod sections for multiple paramcds in one call
make_rdmri_batch <- function(data, vars, cfg) {
  data <- filter_cohort(data, cfg)
  purrr::imap(vars, ~ make_rdmri_section(data, var = .x, cfg = cfg))
}

# Summaries for endpoints where only mean by dose, diff vs 0 mg, and a plot are needed
summarize_mri_endpoint <- function(data, var, cfg) {
  data <- filter_cohort(data, cfg)
  dt <- data %>%
    filter(paramcd == var, avisitn <= 43) %>%
    filter(!is.na(chg), !is.na(trtcd)) |>
    filter(ablfl != "Y")

  if (nrow(dt) == 0 || dplyr::n_distinct(dt$trtcd) < 2) {
    return(NULL)
  }
  
  lineplot <- make_lineplot1(data, var, cfg)

  m <- lme4::lmer(chg ~ base + trtcd + avisitn + (1 | usubjid), data = dt)
  emm <- emmeans::emmeans(m, specs = "trtcd")

  means <- broom::tidy(emm, conf.int = TRUE) %>%
    mutate(
      trt = paste0(trtcd, " mg"),
      mean_txt = format_estimate_ci(estimate, conf.low, conf.high, digits = 2)
    )
  meanstbl <- means |>
    select(trt, mean_txt) |>
    knitr::kable(col.names = c("Dose", "Estimated change from baseline (95% CI)"), digits = 2)

  diffstbl <- emmeans::contrast(emm, method = "trt.vs.ctrl") %>%
    broom::tidy(conf.int = TRUE) %>%
    mutate(
      contrast = paste0(c(25, 50, 100), " mg vs 0 mg"),
      estimate_txt = format_estimate_ci(estimate, conf.low, conf.high, digits = 2)
    ) |>
    select(contrast, estimate_txt) |>
    knitr::kable(col.names = c("Comparison", "Difference (95% CI)"), digits = 2)

  plt <- means %>%
    mutate(trtcd = as.numeric(trtcd)) %>%
    ggplot(aes(x = trtcd, y = estimate)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 1) +
    labs(x = "Dose (mg)", y = "Estimated mean") +
    theme_minimal()
  
  label <- as.character(unique(dt$param))
  
  nsubj <- dt %>%
    distinct(usubjid) %>%
    nrow()

  anydose <- any_dose_vs_0mg(data, var)

  list(
    means = meanstbl,
    lineplot = lineplot,
    diffs = diffstbl,
    anydose = if (is.null(anydose)) NULL else anydose$table,
    plot = plt,
    label = label,
    nsubj = nsubj
  )
}

summarize_mri_cycle11 <- function(data, var, cfg) {
  data <- filter_cohort(data, cfg)
  dt <- data %>%
    filter(paramcd == var, avisitn == 141) %>%
    filter(!is.na(chg), !is.na(trt02p)) |>
    filter(ablfl != "Y")

  if (nrow(dt) == 0 || dplyr::n_distinct(dt$trt02p) < 2) {
    return(NULL)
  }

  m <- lm(chg ~ base + trt02p , data = dt)
  emm <- emmeans::emmeans(m, specs = "trt02p")

  means <- broom::tidy(emm, conf.int = TRUE) %>%
    mutate(
      trt = trt02p,
      trtn = parse_number(trt),
      mean_txt = paste0(
        format(round(estimate, 2), nsmall = 2),
        " (",
        format(round(conf.low, 2), nsmall = 2),
        " to ",
        format(round(conf.high, 2), nsmall = 2),
        ")"
      )
    )
  meanstbl <- means |>
    select(trt, mean_txt) |>
    knitr::kable(col.names = c("Treatment", "Estimated change from baseline (95% CI)"), digits = 2)

  diffstbl <- emmeans::contrast(emm, method = "trt.vs.ctrl", ref = 1) %>%
    broom::tidy(conf.int = TRUE) %>%
    mutate(
      estimate_txt = paste0(
        format(round(estimate, 2), nsmall = 2),
        " (",
        format(round(conf.low, 2), nsmall = 2),
        " to ",
        format(round(conf.high, 2), nsmall = 2),
        ")"
      )
    ) |>
    select(contrast, estimate_txt) |>
    knitr::kable(col.names = c("Comparison", "Difference (95% CI)"), digits = 2)

  plt <- means %>%
    ggplot(aes(x = trtn, y = estimate, group = 1)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
    labs(x = "Dose (mg)", y = "Estimated mean") +
    theme_minimal()

  label <- as.character(unique(dt$param))

  nsubj <- dt %>%
    distinct(usubjid) %>%
    nrow()

  list(means = meanstbl, diffs = diffstbl, plot = plt, label = label, nsubj = nsubj)
}
  

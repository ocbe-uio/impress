
# library(dplyr)
# library(lme4)
# library(ggplot2)
# library(purrr)



#' Build a linear spline basis from analysis days and knot locations.
#'
#' @param ady Numeric vector of analysis days.
#' @param knots Numeric vector of knot locations; will be sorted and deduplicated.
#'
#' @return Tibble with one column per spline segment (s1, s2, ...).
#' @examples
#' make_splines(c(0, 10, 20), c(0, 15, 30))
#'
# Linear spline helper that adapts to knots supplied in cfg$linsp$knots
make_splines <- function(ady, knots) {
  k <- sort(unique(as.numeric(knots)))
  if (length(k) < 2) stop("At least two knots are required.")

  t0 <- pmax(0, ady) # no post-baseline effect before day 0

  segs <- purrr::map2(
    head(k, -1),
    tail(k, -1),
    ~ pmax(0, pmin(t0 - .x, .y - .x))
  )
  tail_seg <- pmax(0, t0 - max(k))

  out <- dplyr::bind_cols(segs)
  names(out) <- paste0("s", seq_along(segs))
  out[[paste0("s", length(segs) + 1)]] <- tail_seg
  tibble::as_tibble(out)
}


#' Fit linear spline mixed model and return predictions and diagnostics.
#'
#' @param dat Data frame already filtered to the desired cohort/parameter, with
#'   columns `aval` (response), `ady` (analysis day), `arm` (treatment factor),
#'   and `subjid` (subject id).
#' @param knots Numeric vector of knot locations for the linear splines.
#'
#' @return List with elements: `fit` (lmer model), `knots` (sorted knots),
#'   `preds` (grid with predictions), `plot` (mean trajectories),
#'   `res_plot` (fitted vs Pearson residuals), `qq_plot` (QQ plot).
#'
# Main entry point for targets
# dat: already filtered dataset (e.g., for cohort/paramcd)
# cfg: list from cfg.yml; expects cfg$linsp$knots
make_linsp <- function(dat, knots) {
  if (is.null(knots)) stop("Knots must be provided")
  knots <- sort(unique(as.numeric(knots)))

  # Build spline basis
  basis <- make_splines(dat$ady, knots)
  basis_names <- names(basis)

  dat_spl <- dat %>%
    mutate(
      arm = factor(arm),
      subjid = factor(subjid)
    ) %>%
    bind_cols(basis)

  rhs <- paste(
    c(
      paste(basis_names, collapse = " + "),
      sprintf("arm:(%s)", paste(basis_names, collapse = " + "))
    ),
    collapse = " + "
  )
  form <- as.formula(sprintf("aval ~ %s + (1 | subjid)", rhs))

  fit <- lmer(form, data = dat_spl)

  newdat <- tidyr::expand_grid(
    ady = seq(min(knots), max(knots), by = 1),
    arm = levels(dat_spl$arm)
  ) %>%
    mutate(arm = factor(arm, levels = levels(dat_spl$arm))) %>%
    bind_cols(make_splines(.$ady, knots))

  newdat <- newdat %>%
    mutate(pred = predict(fit, newdata = ., re.form = NA))

  plt <- ggplot(newdat, aes(x = ady, y = pred, colour = arm)) +
    geom_line() +
    labs(x = "Day", y = "Predicted mean", colour = "ARM") +
    theme_minimal() +
    theme(legend.position = "bottom")
  

diag_df <- tibble(
  fitted  = fitted(fit),
  pearson = residuals(fit, type = "pearson")
)

res_plot <- ggplot(diag_df, aes(x = fitted, y = pearson)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(x = "Fitted values", y = "Pearson residuals") +
  theme_minimal()



qq_plot <- ggplot(diag_df, aes(sample = pearson)) +
  stat_qq(alpha = 0.5) +
  stat_qq_line(color = "red") +
  labs(x = "Theoretical quantiles", y = "Sample quantiles") +
  theme_minimal()

  list(
    fit = fit,
    knots = knots,
    preds = newdat,
    plot = plt, 
    res_plot = res_plot,
    qq_plot = qq_plot
  )
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
    filter(paramcd == var & avisitn <= 43)

  m <- lme4::lmer(aval ~ trtcd + avisitn + (1 | usubjid), data = dt)

  est_raw <- emmeans::emmeans(m, specs = "trtcd")
  S <- vcov(est_raw)

  est <- broom::tidy(est_raw, conf.int = TRUE) %>%
    mutate(
      trt = paste0(trtcd, " mg"),
      mean_txt = paste0(
        round(estimate, digits = 2),
        " (",
        round(conf.low, digits = 2),
        " to ",
        round(conf.high, digits = 2),
        ")"
      )
    )

table <- est %>%
  select(trt, mean_txt) |>
  knitr::kable(col.names = c("Treatment", "Estimated mean (95% CI)"), digits = 2)

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
      ylab = "Mean response (95% CI)",
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
make_rdeff_section <- function(data, var, cfg) {
  cmodels  <- make_cmodels(cfg$doses, cfg$cmodels)
  est      <- make_contest(data, var)
  contmat  <- f_contmat(est, cmodels$fmodels)
  mct      <- f_mcttest(cfg, est, cmodels$fmodels, contmat$contmat)
  fitmods  <- f_fitmod(est, cfg$doses, cmodels$fmodels)

  list(
    text = list(
      intro = "The primary analysis follows the MCP-Mod methodology as detailed in the SAP, using pre-specified candidate dose-response models.",
      models = "The four candidate models (linear, Emax, sigmoid Emax, logistic) are shown in Figure 1."
    ),
    plots = list(
      candidate_models = cmodels$plot,
      fit_models = purrr::map(fitmods$fitmod, "plot")
    ),
    tables = list(
      estimates = est$table,
      contmat = contmat$table,
      mct = mct$p_table,
      fit = fitmods$table
    ),
    objects = list(
      est = est,
      cmodels = cmodels,
      contmat = contmat,
      mct = mct,
      fitmods = fitmods
    )
  )
}

# Batch MCP-Mod sections for multiple paramcds in one call
make_rdeff_batch <- function(data, vars, cfg) {
  purrr::imap(vars, ~ make_rdeff_section(data, var = .x, cfg = cfg))
}

# Summaries for endpoints where only mean by dose, diff vs 0 mg, and a plot are needed
summarize_mri_endpoint <- function(data, var) {
  dt <- data %>%
    filter(paramcd == var, avisitn <= 43) %>%
    filter(!is.na(aval), !is.na(trtcd))
  if (nrow(dt) == 0 || dplyr::n_distinct(dt$trtcd) < 2) return(NULL)

  m <- lme4::lmer(aval ~ trtcd + avisitn + (1 | usubjid), data = dt)
  emm <- emmeans::emmeans(m, specs = "trtcd")

  means <- broom::tidy(emm, conf.int = TRUE) %>%
    mutate(
      trt = paste0(trtcd, " mg"),
      mean_txt = paste0(
        round(estimate, 2),
        " (",
        round(conf.low, 2),
        " to ",
        round(conf.high, 2),
        ")"
      )
    )

  diffs <- emmeans::contrast(emm, method = "trt.vs.ctrl") %>%
    broom::tidy(conf.int = TRUE) %>%
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

  plt <- means %>%
    mutate(trtcd = as.numeric(trtcd)) %>%
    ggplot(aes(x = trtcd, y = estimate)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 1) +
    labs(x = "Dose (mg)", y = "Estimated mean") +
    theme_minimal()

  list(means = means, diffs = diffs, plot = plt)
}
  

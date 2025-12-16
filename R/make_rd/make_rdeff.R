
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



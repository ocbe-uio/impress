make_tte_section <- function(data, paramcd, cfg) {
  dt <- data %>%
    filter_cohort(cfg) %>%
    filter(paramcd == .env$paramcd) %>%
    filter(!is.na(aval), !is.na(cnsr)) %>%
    mutate(
      trtcd = factor(dose02p, levels = c(0, 25, 50, 100)),
      trt = factor(paste0(trtcd, " mg"), levels = c("0 mg", "25 mg", "50 mg", "100 mg"))
    )

  if (nrow(dt) == 0 || dplyr::n_distinct(dt$trtcd) < 2) {
    return(NULL)
  }

  fit <- survival::survfit(survival::Surv(aval, 1L - cnsr) ~ trt, data = dt)

  ggs <- survminer::ggsurvplot(
    fit,
    data = dt,
    risk.table = TRUE,
    pval = FALSE,
    conf.int = FALSE,
    risk.table.height = 0.22,
    risk.table.fontsize = 4,
    legend.title = "Dose",
    legend.labs = levels(dt$trt),
    xlab = "Days since randomisation",
    ylab = "Survival probability"
  )

  sanitize_theme <- function(p) {
    if (is.null(p) || is.null(p$theme)) return(p)
    th <- p$theme
    for (nm in names(th)) {
      if (inherits(th[[nm]], "element_markdown")) {
        th[[nm]] <- ggplot2::element_text()
      }
    }
    p$theme <- th
    p
  }

  ggs$plot <- sanitize_theme(ggs$plot)
  ggs$table <- sanitize_theme(ggs$table)

  logrank <- survival::survdiff(survival::Surv(aval, 1L - cnsr) ~ trt, data = dt)
  pval <- stats::pchisq(logrank$chisq, df = length(logrank$n) - 1, lower.tail = FALSE)
  pval_txt <- ifelse(pval < 0.001, "<0.001", format(round(pval, 3), nsmall = 3))

  cox <- survival::coxph(survival::Surv(aval, 1L - cnsr) ~ trtcd, data = dt)
  cox_tbl <- broom::tidy(cox, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      term = gsub("^trtcd", "", term),
      term = paste0(term, " mg vs 0 mg"),
      estimate_txt = format_estimate_ci(estimate, conf.low, conf.high, digits = 2)
    ) %>%
    select(term, estimate_txt, p.value) %>%
    mutate(p.value = ifelse(p.value < 0.001, "<0.001", format(round(p.value, 3), nsmall = 3))) %>%
    knitr::kable(col.names = c("Comparison", "HR (95% CI)", "P-value"))

  label <- unique(na.omit(dt$param))
  label <- if (length(label)) as.character(label[[1]]) else paramcd

  list(
    label = label,
    km_plot = ggs,
    logrank_p = pval_txt,
    cox_table = cox_tbl
  )
}

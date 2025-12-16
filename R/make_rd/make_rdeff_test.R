
library(targets)
library(tidyverse)
library(lme4)
library(nlme)
library(broom)
library(mmrm)
library(marginaleffects)

adeff <- tar_read(adeff, store = "_targets")




tmp <- adeff |>
    filter(paramcd == "mrt1cf" & cohortcd == 2) |>
    filter(avisitn <= 43) |>
    mutate(trt = case_when(
        avisitn == 15 & str_starts(step, "Cycle 1") ~ paste0(dose01p, " ", dose01pu),
        avisitn == 29 & str_starts(step, "Cycle 1") ~ paste0(dose01p, " ", dose01pu),
        avisitn == 43 & str_starts(step, "Cycle 1") ~ paste0(dose01p, " ", dose01pu),
        avisitn == 29 & str_starts(step, "Cycle 2") ~ paste0(dose01p, " ", dose01pu),
        avisitn == 43 & str_starts(step, "Cycle 2") ~ paste0(dose01p, " ", dose01pu),
        avisitn == 43 & str_starts(step, "Cycle 3") ~ paste0(dose01p, " ", dose01pu),
        TRUE ~ "0mg"
    )) |>
    mutate(
        trt = factor(trt, levels = c("0mg", "25 mg", "50 mg", "100 mg")),
        subjid = factor(subjid),
        avisitn = factor(avisitn, levels = c(-7, 15, 29, 43)
        )
    )
   


m <- lme4::lmer(aval ~ trt + avisitn + (1 | subjid),
    data = tmp
)

summary(m4)

coef(m)
vcov(m)
summary(m)
residuals(m)
 
adeff |>
    group_by(avisit) |>
    summarise(mean = mean(ady, na.rm = TRUE), n= n())





tmp2 <- adeff |>
    filter(paramcd == "mrt1cf" & cohortcd == 2) |>
    filter(avisit %in% c("Randomization And Baseline MRI", "Cycle 2", "Cycle 3", "Cycle 4", "Cycle 11", "Cycle 17"))


m4 <- lme4::lmer(aval ~ trt + avisitn + (1 | subjid),
    data = tmp
)


#m4 %>%
#    avg_predictions(variables = list(armcd = c("AN2_4_50", "AN_4_25"), avisit = c("Cycle 2", "Cycle 3", "Cycle 4", "Cycle 7")))


library(dplyr)
library(splines2) # or splines

knots <- c(0, 15, 29, 45, 146, 238)
bnds  <- range(adeff$ady, na.rm = TRUE)

tmp2 <- adeff %>%
  filter( cohortcd == 2) |>
  mutate(
    # degree = 1 -> linear; use degree = 3 for cubic splines
    bs1 = bSpline(ady, knots = knots, degree = 1, Boundary.knots = bnds)[, 1],
    bs2 = bSpline(ady, knots = knots, degree = 1, Boundary.knots = bnds)[, 2],
    bs3 = bSpline(ady, knots = knots, degree = 1, Boundary.knots = bnds)[, 3],
    bs4 = bSpline(ady, knots = knots, degree = 1, Boundary.knots = bnds)[, 4], 
    bs5 = bSpline(ady, knots = knots, degree = 1, Boundary.knots = bnds)[, 5],
    bs6 = bSpline(ady, knots = knots, degree = 1, Boundary.knots = bnds)[, 6], 
    bs7 = bSpline(ady, knots = knots, degree = 1, Boundary.knots = bnds)[, 7]
  )


m_spline <- lmer(
  aval ~  armcd + bs1 + bs2 + bs3 + bs4 + bs5 + bs6 + (1 | subjid),
  data = tmp2
)

# Counterfactual / standardized marginal predictions by armcd over ADY
ady_grid <- seq(min(tmp2$ady, na.rm = TRUE), max(tmp2$ady, na.rm = TRUE), length.out = 200)

pred_cf_data <- expand_grid(
  armcd = unique(tmp2$armcd),
  ady = ady_grid
) |>
  mutate(
    armcd = factor(armcd, levels = levels(tmp2$armcd)),
    subjid = factor("std", levels = levels(tmp2$subjid)),
    bs1 = bSpline(ady, knots = knots, degree = 1, Boundary.knots = bnds)[, 1],
    bs2 = bSpline(ady, knots = knots, degree = 1, Boundary.knots = bnds)[, 2],
    bs3 = bSpline(ady, knots = knots, degree = 1, Boundary.knots = bnds)[, 3],
    bs4 = bSpline(ady, knots = knots, degree = 1, Boundary.knots = bnds)[, 4],
    bs5 = bSpline(ady, knots = knots, degree = 1, Boundary.knots = bnds)[, 5],
    bs6 = bSpline(ady, knots = knots, degree = 1, Boundary.knots = bnds)[, 6]
  )

pred_cf <- predictions(
  m_spline,
  newdata = pred_cf_data,
  re.form = NA,
  by = c("armcd", "ady")
)

ggplot(pred_cf, aes(x = ady, y = estimate, color = armcd, fill = armcd)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  labs(
    x = "Analysis Day (ADY)",
    y = "Standardized predicted aval",
    color = "ARMCD",
    fill = "ARMCD",
    title = "m_spline standardized trajectories by treatment arm"
  ) +
  theme_minimal()



library(dplyr)
library(lme4)

knots <- c(0, 15, 29, 43, 146, 233)

make_splines <- function(df) {
  df %>%
    mutate(
      s1 = pmin(ady, 15),
      s2 = pmax(pmin(ady, 29)  - 15,  0),
      s3 = pmax(pmin(ady, 43)  - 29,  0),
      s4 = pmax(pmin(ady, 146) - 43,  0),
      s5 = pmax(pmin(ady, 233) - 146, 0),
      s6 = pmax(ady - 233, 0)
    )
}

# Your original data frame is called `dat`
dat2 <- make_splines(tmp)
dat2$armcd <- factor(dat2$armcd)

# Mixed model with random intercept for subjid
fit <- lmer(
  aval ~ (s1 + s2 + s3 + s4 + s5 + s6) * armcd + (1 | subjid),
  data = dat2
)
summary(fit)


library(ggplot2)

# Prediction grid over time and groups
newdat <- expand.grid(
  ady = seq(0, 233, by = 1),
  armcd   = levels(dat2$armcd)
)

newdat <- make_splines(newdat)
newdat$armcd <- factor(newdat$armcd, levels = levels(dat2$armcd))
# Predictions using fixed effects only
newdat$y_hat <- predict(fit, newdata = newdat, re.form = NA)

ggplot(newdat, aes(x = ady, y = y_hat, colour = armcd)) +
  geom_line() +
  labs(x = "Day (ADY)", y = "Predicted mean of y") +
  theme_minimal() +
  #make y axis vary between 0 and 10 
  ylim(0, 10)


library(dplyr)

knots <- c(0, 15, 29, 43, 146, 233)

dat <- tmp2 %>%
  mutate(
    # time, but not allowed to be negative (no post-baseline time < 0)
    t0 = pmax(0, ady),

    # linear spline basis (each term is 0 at t0 = 0)
    s1 = pmin(t0, 15),
    s2 = pmax(0, pmin(t0 - 15, 29 - 15)),
    s3 = pmax(0, pmin(t0 - 29, 43 - 29)),
    s4 = pmax(0, pmin(t0 - 43, 146 - 43)),
    s5 = pmax(0, pmin(t0 - 146, 233 - 146)),
    s6 = pmax(0, t0 - 233)  # last segment beyond 233
  )


library(lme4)

fit <- lmer(
  aval ~ s1 + s2 + s3 + s4 + s5 + s6 + 
      armcd:(s1 + s2 + s3 + s4 + s5 + s6) + 
      (1 | subjid),
  data = dat
)


newdat <- expand.grid(
  ady = seq(-1, 234, by = 1),
  armcd   = levels(dat$armcd)
) %>%
  mutate(
    t0 = pmax(0, ady),
    s1 = pmin(t0, 15),
    s2 = pmax(0, pmin(t0 - 15, 29 - 15)),
    s3 = pmax(0, pmin(t0 - 29, 43 - 29)),
    s4 = pmax(0, pmin(t0 - 43, 146 - 43)),
    s5 = pmax(0, pmin(t0 - 146, 233 - 146)),
    s6 = pmax(0, t0 - 233)
  )

newdat$y_hat <- predict(fit, newdata = newdat, re.form = NA)



library(dplyr)
library(lme4)

# your knots
knots <- c(0, 15, 29, 43, 146, 233)

make_splines <- function(ady) {
  t0 <- pmax(0, ady)  # truncation at 0 ⇒ no post-baseline effect before 0
  
  tibble(
    s1 = pmin(t0, 15),
    s2 = pmax(0, pmin(t0 - 15, 29 - 15)),
    s3 = pmax(0, pmin(t0 - 29, 43 - 29)),
    s4 = pmax(0, pmin(t0 - 43, 146 - 43)),
    s5 = pmax(0, pmin(t0 - 146, 233 - 146)),
    s6 = pmax(0, t0 - 233)
  )
}


dat_spl <- adeff %>%
  filter(paramcd == "mrt1cf" & cohortcd == 2) |>
  mutate(armcd = factor(armcd)) %>%
  bind_cols(make_splines(.$ady))

# model: one common baseline curve, treatment differences only via splines
fit <- lmer(
  aval ~ s1 + s2 + s3 + s4 + s5 + s6 + 
      armcd:(s1 + s2 + s3 + s4 + s5 + s6) +  # differences only after 0
      (1 | subjid),
  data = dat_spl
)

newdat <- expand.grid(
  ady   = seq(0, 233, by = 1),
  armcd = levels(dat_spl$armcd)
) %>%
  as_tibble() %>%
  mutate(
    armcd = factor(armcd, levels = levels(dat_spl$armcd))
  ) %>%
  bind_cols(make_splines(.$ady))

newdat$y_hat <- predict(fit, newdata = newdat, re.form = NA)

ggplot(newdat, aes(x = ady, y = y_hat, colour = armcd)) +
  geom_line() +
  labs(x = "Day (ADY)", y = "Predicted mean of y") +
  theme_minimal() 

plot(fit)
qqnorm(residuals(fit))
qqline(residuals(fit))


diag_df <- tibble(
  fitted  = fitted(fit),
  pearson = residuals(fit, type = "pearson")
)

ggplot(diag_df, aes(x = fitted, y = pearson)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(x = "Fitted values", y = "Pearson residuals") +
  theme_minimal()




ggplot(diag_df, aes(sample = pearson)) +
  stat_qq(alpha = 0.5) +
  stat_qq_line(color = "red") +
  labs(x = "Theoretical quantiles", y = "Sample quantiles") +
  theme_minimal()

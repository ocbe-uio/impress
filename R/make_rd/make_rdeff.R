
library(targets)
library(tidyverse)
library(lme4)
library(nlme)
library(broom)
library(mmrm)
library(marginaleffects)

tar_load(adeff)




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
   
m <- mmrm(aval ~ trt + avisitn + us(avisitn | subjid),
    data = tmp
)



m2 <- mmrm(aval ~ trt + avisitn + ar1(avisitn | subjid),
    data = tmp
)


m3 <- mmrm(aval ~ trt + avisitn + cs(avisitn | subjid),
    data = tmp
)

m4 <- lme4::lmer(aval ~ trt + avisitn + (1 | subjid),
    data = tmp
)

summary(m4)

coef(m)
vcov(m)
summary(m)
residuals(m)
 

tmp2 <- adeff |>
    filter(paramcd == "mrt1cf" & cohortcd == 2) |>
    filter(avisit %in% c("Randomization And Baseline MRI", "Cycle 2", "Cycle 3", "Cycle 4", "Cycle 11", "Cycle 17"))


m21 <- mmrm(aval ~ avisit + us(avisit | subjid),
    data = tmp2
)

m22 <- mmrm(aval ~ avisit + ar1(avisit | subjid),
    data = tmp2
)

m23 <- mmrm(aval ~ avisit + cs(avisit | subjid),
    data = tmp2
)

summary(m21)
summary(m22)
summary(m23)

tmp2 |>
    group_by(avisit) |>
    summarise(n = n())



m2 %>%
    avg_predictions(variables = list(armcd = c("AN2_4_50", "AN_4_25"), avisit = c("Cycle 2", "Cycle 3", "Cycle 4", "Cycle 7")))


library(dplyr)
library(splines2) # or splines

knots <- c(-7, 15, 29, 43)
bnds  <- range(tmp$ady, na.rm = TRUE)

tmp2 <- tmp %>%
  mutate(
    # degree = 1 -> linear; use degree = 3 for cubic splines
    bs1 = bSpline(ady, knots = knots, degree = 1, Boundary.knots = bnds)[, 1],
    bs2 = bSpline(ady, knots = knots, degree = 1, Boundary.knots = bnds)[, 2],
    bs3 = bSpline(ady, knots = knots, degree = 1, Boundary.knots = bnds)[, 3],
    bs4 = bSpline(ady, knots = knots, degree = 1, Boundary.knots = bnds)[, 4]
  )


m_spline <- mmrm(
  aval ~ trt + bs1 + bs2 + bs3 + bs4,
  data = tmp2,
  covariance = cs(ady | subjid),  # or cs(avisitn | subjid) if you prefer
  reml = FALSE
)

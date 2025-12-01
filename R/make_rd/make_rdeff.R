library(targets)
library(tidyverse)
library(lme4)
library(nlme)
library(broom)
library(mmrm)

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
    
    

m <- mmrm(aval ~trt + avisitn + us(avisitn | subjid),
        data = tmp
)
summary(m3)

coef(m3)
vcov(m3)
summary(m)
residuals(m)
 |>
    lmer


library(mmrm)
fit <- mmrm(
    formula = FEV1 ~ RACE + SEX + ARMCD * AVISIT + us(AVISIT | USUBJID),
    data = fev_data
)

fev_data
View(fev_data)


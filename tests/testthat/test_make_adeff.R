# Unit tests for make_adeff()

library(testthat)
library(dplyr)

source(file.path("R", "make_ad", "make_adeff.R"))

test_that("make_adeff produces at most one baseline per usubjid/paramcd", {
  skip_if_not_installed("targets")
  skip_if_not_installed("admiral")
  suppressPackageStartupMessages(library(admiral))

  cfg <- targets::tar_read(cfg, store = "_targets")
  shamraw <- targets::tar_read(shamraw, store = "_targets")
  adsl <- targets::tar_read(adsl, store = "_targets")

  adeff <- make_adeff(shamraw, adsl, cfg)

  baseline_counts <- adeff %>%
    filter(ablfl == "Y") %>%
    count(usubjid, paramcd, name = "n_baseline") %>%
    filter(n_baseline > 1)

  expect_equal(nrow(baseline_counts), 0)
})

test_that("make_adeff baseline occurs on/before randomization date", {
  skip_if_not_installed("targets")
  skip_if_not_installed("admiral")
  suppressPackageStartupMessages(library(admiral))

  cfg <- targets::tar_read(cfg, store = "_targets")
  shamraw <- targets::tar_read(shamraw, store = "_targets")
  adsl <- targets::tar_read(adsl, store = "_targets")

  adeff <- make_adeff(shamraw, adsl, cfg)

  baseline_after_rand <- adeff %>%
    filter(ablfl == "Y") %>%
    filter(!is.na(adt), !is.na(randdt)) %>%
    filter(adt > randdt)

  expect_equal(nrow(baseline_after_rand), 0)
})

test_that("make_adeff trtcd matches stepped-wedge armcd schedule", {
  skip_if_not_installed("targets")
  skip_if_not_installed("admiral")
  suppressPackageStartupMessages(library(admiral))

  cfg <- targets::tar_read(cfg, store = "_targets")
  shamraw <- targets::tar_read(shamraw, store = "_targets")
  adsl <- targets::tar_read(adsl, store = "_targets")

  adeff <- make_adeff(shamraw, adsl, cfg)

  ar_trt <- function(day, start_day, dose, post_cycle4_dose) {
    case_when(
      is.na(day) ~ NA_real_,
      day <= start_day ~ 0,
      day > start_day & day <= 43 ~ dose,
      day > 43 & day <= 169 ~ post_cycle4_dose,
      day > 169 ~ 0
    )
  }

  an_trt <- function(day, start_day, dose, post_cycle4_dose) {
    case_when(
      is.na(day) ~ NA_real_,
      day <= start_day ~ 0,
      day > start_day & day <= 43 ~ dose,
      day > 43 & day <= 239 ~ post_cycle4_dose,
      day > 239 ~ 0
    )
  }

  bm_trt <- function(day, start_day, dose) {
    case_when(
      is.na(day) ~ NA_real_,
      day < start_day ~ 0,
      day >= start_day & day <= 1270 ~ dose,
      day > 270 ~ 0
    )
  }

  comp <- adeff %>%
    mutate(
      visit_day = suppressWarnings(as.numeric(as.character(avisitn))),
      trtcd_num = suppressWarnings(as.numeric(as.character(trtcd))),
      expected_trtcd = case_when(
        armcd == "AR1_1_0" ~ ar_trt(visit_day, 1, 25, 0),
        armcd == "AR1_1_25" ~ ar_trt(visit_day, 1, 25, 25),
        armcd == "AR1_2_0" ~ ar_trt(visit_day, 15, 25, 0),
        armcd == "AR1_2_25" ~ ar_trt(visit_day, 15, 25, 25),
        armcd == "AR1_3_0" ~ ar_trt(visit_day, 29, 25, 0),
        armcd == "AR1_3_25" ~ ar_trt(visit_day, 29, 25, 25),
        armcd == "AR2_4_0" ~ ar_trt(visit_day, 1, 50, 0),
        armcd == "AR2_4_50" ~ ar_trt(visit_day, 1, 50, 50),
        armcd == "AR2_5_0" ~ ar_trt(visit_day, 15, 50, 0),
        armcd == "AR2_5_50" ~ ar_trt(visit_day, 15, 50, 50),
        armcd == "AR2_6_0" ~ ar_trt(visit_day, 29, 50, 0),
        armcd == "AR2_6_50" ~ ar_trt(visit_day, 29, 50, 50),
        armcd == "AR3_1_100" ~ ar_trt(visit_day, 1, 100, 100),
        armcd == "AR3_2_100" ~ ar_trt(visit_day, 15, 100, 100),
        armcd == "AR3_3_100" ~ ar_trt(visit_day, 29, 100, 100),
        armcd == "AN1_1_0" ~ an_trt(visit_day, 1, 25, 0),
        armcd == "AN1_1_25" ~ an_trt(visit_day, 1, 25, 25),
        armcd == "AN1_2_0" ~ an_trt(visit_day, 15, 25, 0),
        armcd == "AN1_2_25" ~ an_trt(visit_day, 15, 25, 25),
        armcd == "AN1_3_0" ~ an_trt(visit_day, 29, 25, 0),
        armcd == "AN1_3_25" ~ an_trt(visit_day, 29, 25, 25),
        armcd == "AN2_4_0" ~ an_trt(visit_day, 1, 50, 0),
        armcd == "AN2_4_50" ~ an_trt(visit_day, 1, 50, 50),
        armcd == "AN2_5_0" ~ an_trt(visit_day, 15, 50, 0),
        armcd == "AN2_5_50" ~ an_trt(visit_day, 15, 50, 50),
        armcd == "AN2_6_0" ~ an_trt(visit_day, 29, 50, 0),
        armcd == "AN2_6_50" ~ an_trt(visit_day, 29, 50, 50),
        armcd == "AN3_1_100" ~ an_trt(visit_day, 1, 100, 100),
        armcd == "AN3_2_100" ~ an_trt(visit_day, 15, 100, 100),
        armcd == "AN3_3_100" ~ an_trt(visit_day, 29, 100, 100),
        armcd == "BM1" ~ bm_trt(visit_day, 1, 50),
        armcd == "BM2" ~ bm_trt(visit_day, 1091, 50),
        armcd == "BM3" ~ bm_trt(visit_day, 1181, 50),
        TRUE ~ NA_real_
      )
    ) %>%
    filter(!is.na(expected_trtcd), !is.na(trtcd_num))

  if (nrow(comp) == 0) {
    skip("No records with mapped armcd and trtcd to compare.")
  }

  mismatch <- comp %>%
    filter(trtcd_num != expected_trtcd)

  expect_equal(nrow(mismatch), 0)
})

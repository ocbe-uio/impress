# Unit tests for make_rdmri helpers

library(testthat)
library(dplyr)

source(file.path("R", "make_rd", "make_rdmri.R"))

test_that("make_lineplot1 returns a ggplot object", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("emmeans")
  skip_if_not_installed("ggplot2")

  set.seed(1)
  toy <- tidyr::expand_grid(
    usubjid = factor(1:8),
    avisitn = c(-7, 15, 29, 43),
    step = factor(c("Cycle 1", "Cycle 2")),
    paramcd = "mrt1cf",
    param = "Test"
  ) %>%
    mutate(
      base_raw = if_else(avisitn == -7, rnorm(n(), 10, 0.5), NA_real_)
    )

  base_map <- toy %>%
    filter(avisitn == -7) %>%
    group_by(usubjid) %>%
    summarise(base = mean(base_raw, na.rm = TRUE), .groups = "drop")

  toy <- toy %>%
    left_join(base_map, by = "usubjid") %>%
    mutate(
      aval = if_else(avisitn == -7, base, base + rnorm(n(), 1, 0.3))
    ) %>%
    select(-base_raw)

  plt <- make_lineplot1(toy, "mrt1cf")
  expect_s3_class(plt, "ggplot")
})

test_that("make_contest returns expected components", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("emmeans")

  set.seed(2)
  toy <- tidyr::expand_grid(
    usubjid = factor(1:10),
    avisitn = c(15, 29, 43),
    trtcd = factor(c(0, 25, 50)),
    paramcd = "mrt1cf"
  ) %>%
    mutate(
      base = rnorm(n(), 10, 0.5),
      chg = rnorm(n(), 0, 1),
      aval = base + chg,
      ablfl = "N"
    )

  res <- make_contest(toy, "mrt1cf")
  expect_named(res, c("model", "estimates", "S", "emmeans", "table"))
  expect_true(inherits(res$model, "merMod"))
  expect_true(all(c("estimate", "conf.low", "conf.high") %in% names(res$estimates)))
})

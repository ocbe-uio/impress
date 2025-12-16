# Unit tests for make_splines() and make_linsp()

library(testthat)
library(dplyr)
library(lme4)

# Source the functions under test
source(file.path("R", "make_rd", "make_rdeff.R"))

test_that("make_splines builds expected number of segments and values", {
  knots <- c(0, 15, 30)
  ady <- c(-5, 0, 5, 15, 25, 40)
  spl <- make_splines(ady, knots)

  expect_s3_class(spl, "tbl_df")
  expect_equal(ncol(spl), length(knots)) # s1..s3
  expect_true(all(spl >= 0))

  # Spot-check: before 0 should be zero
  expect_equal(spl$s1[1], 0)
  # Between 0 and 15 lives entirely in s1
  expect_gt(spl$s1[3], 0)
  expect_equal(spl$s2[3], 0)
  # After last knot flows into tail segment
  expect_gt(spl$s3[length(ady)], 0)
})

test_that("make_linsp returns model, predictions, and diagnostics", {
  set.seed(123)
  knots <- c(0, 15, 30)

  toy <- tidyr::expand_grid(
    subjid = factor(1:6),
    arm = factor(c("A", "B")),
    ady = c(0, 10, 20, 30)
  ) %>%
    mutate(
      aval = 10 + 0.5 * ady + if_else(arm == "B", 2, 0) + rnorm(n(), sd = 0.1)
    )

  res <- make_linsp(toy, knots)

  expect_named(res, c("fit", "knots", "preds", "plot", "res_plot", "qq_plot"))
  expect_true(inherits(res$fit, "merMod"))
  expect_equal(res$knots, sort(unique(knots)))

  # Pred grid should span knots range and both arms
  expect_true(all(c(min(res$preds$ady), max(res$preds$ady)) == range(knots)))
  expect_setequal(unique(res$preds$arm), levels(toy$arm))
  expect_false(any(is.na(res$preds$pred)))

  # Plots are ggplot objects
  expect_s3_class(res$plot, "ggplot")
  expect_s3_class(res$res_plot, "ggplot")
  expect_s3_class(res$qq_plot, "ggplot")
})

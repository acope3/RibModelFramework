library(testthat)
library(ribModel)

context("Utility")

test_that("general utility functions", {
  expect_equal(testUtility(), 0)
})

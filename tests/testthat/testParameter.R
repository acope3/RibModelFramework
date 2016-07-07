library(testthat)
library(ribModel)

context("Parameter")

# TODO: Fix the R vs C side bug in initParameterSet
test_that("general parameter functions", {
  expect_equal(testParameter(), 0)
})

# TODO: Implement the following
p <-new(Parameter)

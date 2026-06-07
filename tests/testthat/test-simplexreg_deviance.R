
library(testthat)
library(SimplexRegression)

# ==============================================================================
# Tests for dev.unit.simplex function
# ==============================================================================

test_that("dev.unit.simplex validates inputs correctly", {
  # y must be in (0, 1)
  expect_error(dev.unit.simplex(y = 0, mu = 0.5), "'y' must always be in \\(0, 1\\)")
  expect_error(dev.unit.simplex(y = 1, mu = 0.5), "'y' must always be in \\(0, 1\\)")
  expect_error(dev.unit.simplex(y = -0.5, mu = 0.5), "'y' must always be in \\(0, 1\\)")
  expect_error(dev.unit.simplex(y = 1.5, mu = 0.5), "'y' must always be in \\(0, 1\\)")

  # mu must be in (0, 1)
  expect_error(dev.unit.simplex(y = 0.5, mu = 0), "parameter 'mu' must always be in \\(0, 1\\)")
  expect_error(dev.unit.simplex(y = 0.5, mu = 1), "parameter 'mu' must always be in \\(0, 1\\)")
  expect_error(dev.unit.simplex(y = 0.5, mu = -0.5), "parameter 'mu' must always be in \\(0, 1\\)")
  expect_error(dev.unit.simplex(y = 0.5, mu = 1.5), "parameter 'mu' must always be in \\(0, 1\\)")
})

test_that("dev.unit.simplex accepts valid inputs without error", {
  expect_no_error(dev.unit.simplex(y = 0.5, mu = 0.5))
  expect_no_error(dev.unit.simplex(y = 0.3, mu = 0.7))
  expect_no_error(dev.unit.simplex(y = 0.8, mu = 0.2))
  expect_no_error(dev.unit.simplex(y = 0.1, mu = 0.9))
})

test_that("dev.unit.simplex returns correct type and dimensions", {
  # scalar input
  result <- dev.unit.simplex(y = 0.5, mu = 0.5)
  expect_type(result, "double")
  expect_length(result, 1)

  # vector input same length
  result <- dev.unit.simplex(y = c(0.2, 0.5, 0.8), mu = c(0.3, 0.5, 0.7))
  expect_type(result, "double")
  expect_length(result, 3)

  # recycling (y vector, mu scalar)
  result <- dev.unit.simplex(y = c(0.2, 0.5, 0.8), mu = 0.5)
  expect_length(result, 3)

  # recycling (y scalar, mu vector)
  result <- dev.unit.simplex(y = 0.5, mu = c(0.2, 0.5, 0.8))
  expect_length(result, 3)
})

test_that("dev.unit.simplex returns zero when y equals mu", {
  # Perfect fit should return zero deviance
  expect_equal(dev.unit.simplex(y = 0.5, mu = 0.5), 0)
  expect_equal(dev.unit.simplex(y = 0.3, mu = 0.3), 0)
  expect_equal(dev.unit.simplex(y = 0.7, mu = 0.7), 0)
  expect_equal(dev.unit.simplex(y = 0.1, mu = 0.1), 0)
  expect_equal(dev.unit.simplex(y = 0.9, mu = 0.9), 0)
})

test_that("dev.unit.simplex returns non-negative values", {
  # All deviances should be >= 0
  expect_true(dev.unit.simplex(y = 0.6, mu = 0.5) >= 0)
  expect_true(dev.unit.simplex(y = 0.4, mu = 0.5) >= 0)
  expect_true(all(dev.unit.simplex(y = c(0.2, 0.5, 0.8), mu = c(0.3, 0.5, 0.7)) >= 0))
})

test_that("dev.unit.simplex handles numeric precision near boundaries", {
  # Values close to boundaries (but still valid)
  expect_true(is.finite(dev.unit.simplex(y = 1e-8, mu = 0.5)))
  expect_true(is.finite(dev.unit.simplex(y = 1 - 1e-8, mu = 0.5)))
  expect_true(is.finite(dev.unit.simplex(y = 0.5, mu = 1e-8)))
  expect_true(is.finite(dev.unit.simplex(y = 0.5, mu = 1 - 1e-8)))

  # Should still be positive
  expect_true(dev.unit.simplex(y = 1e-8, mu = 0.5) > 0)
  expect_true(dev.unit.simplex(y = 0.5, mu = 1e-8) > 0)
})

test_that("dev.unit.simplex is consistent with known values", {
  # Test some specific values for correctness
  # For y = 0.6, mu = 0.5
  # diff = 0.1
  # muonemu = 0.5 * 0.5 = 0.25
  # yoneminy = 0.6 * 0.4 = 0.24
  # deviance = (0.1 / 0.25)^2 / 0.24 = (0.4)^2 / 0.24 = 0.16 / 0.24 = 0.6666667

  result <- dev.unit.simplex(y = 0.6, mu = 0.5)
  expect_equal(result, 0.16 / 0.24, tolerance = 1e-12)
  expect_equal(result, 2/3, tolerance = 1e-12)

  # For y = 0.7, mu = 0.5
  # diff = 0.2
  # muonemu = 0.25
  # yoneminy = 0.7 * 0.3 = 0.21
  # deviance = (0.2 / 0.25)^2 / 0.21 = (0.8)^2 / 0.21 = 0.64 / 0.21 ≈ 3.047619

  result2 <- dev.unit.simplex(y = 0.7, mu = 0.5)
  expect_equal(result2, 0.64 / 0.21, tolerance = 1e-12)
})

test_that("dev.unit.simplex handles NA values appropriately", {
  # Note: This behavior depends on your stopifnot implementation
  # If you want to allow NA, adjust accordingly
  expect_error(dev.unit.simplex(y = NA, mu = 0.5))
  expect_error(dev.unit.simplex(y = 0.5, mu = NA))
})

test_that("dev.unit.simplex formula matches expected mathematical properties", {
  # The deviance should be scale invariant in a specific way
  y <- 0.6
  mu <- 0.4

  # Direct calculation
  direct <- dev.unit.simplex(y, mu)

  # Alternative calculation using formula components
  diff <- y - mu
  mu_term <- mu * (1 - mu)
  y_term <- y * (1 - y)
  alt <- (diff / mu_term)^2 / y_term

  expect_equal(direct, alt, tolerance = 1e-12)
})

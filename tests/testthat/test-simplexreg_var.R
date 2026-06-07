library(testthat)
library(SimplexRegression)

test_that("variance.simplex validates inputs correctly", {
  # mu must be in (0, 1)
  expect_error(variance.simplex(mu = 0,    sigma2 = 0.5), "parameter 'mu' must always be in \\(0, 1\\)")
  expect_error(variance.simplex(mu = 1,    sigma2 = 0.5), "parameter 'mu' must always be in \\(0, 1\\)")
  expect_error(variance.simplex(mu = -0.5, sigma2 = 0.5), "parameter 'mu' must always be in \\(0, 1\\)")
  expect_error(variance.simplex(mu = 1.5,  sigma2 = 0.5), "parameter 'mu' must always be in \\(0, 1\\)")

  # sigma2 must be positive
  expect_error(variance.simplex(mu = 0.5, sigma2 = 0),    "parameter 'sigma2' must always be positive")
  expect_error(variance.simplex(mu = 0.5, sigma2 = -0.1), "parameter 'sigma2' must always be positive")

  # NA values
  expect_error(variance.simplex(mu = NA,  sigma2 = 0.5), "parameter 'mu' must always be in \\(0, 1\\)")
  expect_error(variance.simplex(mu = 0.5, sigma2 = NA),  "parameter 'sigma2' must always be positive")
})

test_that("variance.simplex accepts valid inputs without error", {
  expect_no_error(variance.simplex(mu = 0.5, sigma2 = 0.5))
  expect_no_error(variance.simplex(mu = 0.1, sigma2 = 0.1))
  expect_no_error(variance.simplex(mu = 0.9, sigma2 = 1.0))
})

test_that("variance.simplex returns correct type and dimensions", {
  # scalar input
  result <- variance.simplex(mu = 0.5, sigma2 = 0.1)
  expect_type(result, "double")
  expect_length(result, 1)

  # vector input
  result <- variance.simplex(mu = c(0.3, 0.5, 0.7), sigma2 = c(0.1, 0.15, 0.2))
  expect_type(result, "double")
  expect_length(result, 3)

  # recycling
  result <- variance.simplex(mu = c(0.3, 0.5, 0.7), sigma2 = 0.1)
  expect_length(result, 3)
})

test_that("variance.simplex returns values within theoretical bounds", {
  mu <- 0.5
  bernoulli_var <- mu * (1 - mu)  # theoretical upper bound

  result <- variance.simplex(mu = mu, sigma2 = 0.5)
  expect_true(result > 0)
  expect_true(result < bernoulli_var)

  # all vector values must respect the bounds
  mu_vec <- c(0.2, 0.5, 0.8)
  result_vec <- variance.simplex(mu = mu_vec, sigma2 = 0.3)
  expect_true(all(result_vec > 0))
  expect_true(all(result_vec < mu_vec * (1 - mu_vec)))
})

test_that("variance.simplex is symmetric around mu = 0.5", {
  expect_equal(
    variance.simplex(mu = 0.3, sigma2 = 0.5),
    variance.simplex(mu = 0.7, sigma2 = 0.5),
    tolerance = 1e-10
  )
  expect_equal(
    variance.simplex(mu = 0.2, sigma2 = 0.3),
    variance.simplex(mu = 0.8, sigma2 = 0.3),
    tolerance = 1e-10
  )
})

test_that("variance.simplex is monotonically increasing in sigma2", {
  sigma2_vals <- c(0.1, 0.5, 1, 10)
  vars <- sapply(sigma2_vals, function(s) variance.simplex(mu = 0.5, sigma2 = s))

  for (i in seq_len(length(vars) - 1)) {
    expect_true(vars[i] < vars[i + 1])
  }
})

test_that("variance.simplex converges to mu*(1-mu) as sigma2 grows large", {
  mu <- 0.5
  bernoulli_var <- mu * (1 - mu)

  # Variance should increase toward bernoulli_var as sigma2 grows
  result_small  <- variance.simplex(mu = mu, sigma2 = 1)
  result_medium <- variance.simplex(mu = mu, sigma2 = 100)
  result_large  <- variance.simplex(mu = mu, sigma2 = 1e6)

  expect_true(result_small < result_medium)
  expect_true(result_medium < result_large)
  expect_true(result_large < bernoulli_var)
  expect_true(abs(result_large - bernoulli_var) < 0.01)
})

test_that("variance.simplex handles extreme values of mu and sigma2", {
  expect_true(is.finite(variance.simplex(mu = 0.001, sigma2 = 0.5)))
  expect_true(is.finite(variance.simplex(mu = 0.999, sigma2 = 0.5)))
  expect_true(is.finite(variance.simplex(mu = 0.5,   sigma2 = 1e-6)))
  expect_true(is.finite(variance.simplex(mu = 0.5,   sigma2 = 1e6)))
  expect_true(variance.simplex(mu = 0.5, sigma2 = 1e-6) > 0)
})

test_that("variance.simplex uses asymptotic approximation when a > 700", {
  sigma2_small <- 1e-8
  mu <- 0.5
  # a = 1 / (2 * sigma2 * (mu*(1-mu))^2) must be > 700
  a <- 1 / (2 * sigma2_small * (mu * (1 - mu))^2)
  expect_true(a > 700)

  result <- variance.simplex(mu = mu, sigma2 = sigma2_small)
  expect_true(is.finite(result))
  expect_true(result > 0)
})

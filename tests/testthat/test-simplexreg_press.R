library(testthat)
library(SimplexRegression)

# ---------------------------------------------------------------------------
# Fixture: simulate data and fit models once for all tests
# ---------------------------------------------------------------------------
set.seed(2026)
n <- 100
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
z1 <- runif(n, 0, 1)
mu <- parametric_mean_link_inv(0.6 - 2*x1 - 1.5*x2, 0.5, "plogit1")
sigma2 <- dispersion_link_inv(-2 - 2.5*z1, "log")
y <- rsimplex(n, mu, sigma2)
data <- data.frame(y = y, x1 = x1, x2 = x2, z1 = z1)

fit_plogit1 <- simplexreg(y ~ x1 + x2 | z1, data = data, link.mu = "plogit1")
fit_logit   <- simplexreg(y ~ x1 + x2 | z1, data = data, link.mu = "logit")

# ---------------------------------------------------------------------------

test_that("press requires objects of class simplexregression", {
  expect_error(press(list()),           "class 'simplexregression'")
  expect_error(press(fit_plogit1, "x"), "class 'simplexregression'")
})

test_that("press rejects invalid type", {
  expect_error(press(fit_plogit1, type = "invalid"))
})

test_that("press returns a named numeric vector for a single model", {
  result <- press(fit_plogit1)
  expect_true(is.numeric(result))
  expect_named(result, c("P2", "P2_c", "PRESS"))
})

test_that("press returns a data.frame for multiple models", {
  result <- press(fit_plogit1, fit_logit)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("P2", "P2_c", "PRESS"))
  expect_equal(nrow(result), 2)
})

test_that("press uses named arguments as row names", {
  result <- press(m1 = fit_plogit1, m2 = fit_logit)
  expect_equal(rownames(result), c("m1", "m2"))
})

test_that("press returns finite values for a single model", {
  result <- press(fit_plogit1)
  expect_true(all(is.finite(result)))
})

test_that("press returns finite values for multiple models", {
  result <- press(fit_plogit1, fit_logit)
  # data.frame: check each column
  expect_true(all(is.finite(result$P2)))
  expect_true(all(is.finite(result$P2_c)))
  expect_true(all(is.finite(result$PRESS)))
})

test_that("press PRESS statistic is positive", {
  result <- press(fit_plogit1)
  expect_true(result["PRESS"] > 0)
})

test_that("press works with type = 'biasvariance'", {
  result <- press(fit_plogit1, type = "biasvariance")
  expect_named(result, c("P2", "P2_c", "PRESS"))
  expect_true(all(is.finite(result)))
})

test_that("press works with fixed link model", {
  result <- press(fit_logit)
  expect_named(result, c("P2", "P2_c", "PRESS"))
  expect_true(all(is.finite(result)))
  expect_true(result["PRESS"] > 0)
})

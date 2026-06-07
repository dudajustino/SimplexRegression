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
fit_plogit2 <- simplexreg(y ~ x1 + x2 | z1, data = data, link.mu = "plogit2")
fit_logit   <- simplexreg(y ~ x1 + x2 | z1, data = data, link.mu = "logit")
fit_probit  <- simplexreg(y ~ x1 + x2 | z1, data = data, link.mu = "probit")

# ---------------------------------------------------------------------------

test_that("penalized.ss requires at least 2 models", {
  expect_error(penalized.ss(fit_plogit1), "At least 2 models are required")
})

test_that("penalized.ss requires objects of class simplexregression", {
  expect_error(penalized.ss(fit_plogit1, list()),    "class 'simplexregression'")
  expect_error(penalized.ss(fit_plogit1, "model"),   "class 'simplexregression'")
})

test_that("penalized.ss returns a data.frame with correct columns", {
  result <- penalized.ss(fit_plogit1, fit_plogit2, kappa = 0.1, verbose = FALSE)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("df", "SS"))
  expect_equal(nrow(result), 2)
})

test_that("penalized.ss uses named arguments as row names", {
  result <- penalized.ss(m1 = fit_plogit1, m2 = fit_plogit2,
                         kappa = 0.1, verbose = FALSE)
  expect_equal(rownames(result), c("m1", "m2"))
})

test_that("penalized.ss returns finite numeric SS values", {
  result <- penalized.ss(fit_plogit1, fit_plogit2, kappa = 0.1, verbose = FALSE)
  expect_true(all(is.finite(result$SS)))
  expect_true(all(is.numeric(result$SS)))
})

test_that("penalized.ss works with kappa = 0 (standard Scout Score)", {
  result <- penalized.ss(fit_logit, fit_probit, kappa = 0, verbose = FALSE)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_true(all(is.finite(result$SS)))
})

test_that("penalized.ss works with more than 2 models", {
  result <- penalized.ss(fit_plogit1, fit_plogit2, fit_logit, fit_probit,
                         kappa = 0, verbose = FALSE)
  expect_equal(nrow(result), 4)
  expect_true(all(is.finite(result$SS)))
})

test_that("penalized.ss warns and sets kappa = 0 when models lack parametric link", {
  expect_warning(
    penalized.ss(fit_plogit1, fit_logit, kappa = 0.1, verbose = FALSE),
    "kappa > 0 is only valid"
  )
})

test_that("penalized.ss returns invisibly when verbose = TRUE", {
  result <- withVisible(
    penalized.ss(fit_plogit1, fit_plogit2, kappa = 0.1, verbose = TRUE)
  )
  expect_false(result$visible)
})

test_that("penalized.ss returns visibly when verbose = FALSE", {
  result <- withVisible(
    penalized.ss(fit_plogit1, fit_plogit2, kappa = 0.1, verbose = FALSE)
  )
  expect_true(result$visible)
})

test_that("penalized.ss df column matches number of estimated parameters", {
  result <- penalized.ss(fit_plogit1, fit_plogit2, kappa = 0.1, verbose = FALSE)
  expect_equal(result$df[1], n - fit_plogit1$df.residual)
  expect_equal(result$df[2], n - fit_plogit2$df.residual)
})

test_that("penalized.ss silently sets kappa = 0 when models lack parametric link and kappa not specified", {
  result <- expect_no_warning(
    penalized.ss(fit_logit, fit_probit, verbose = TRUE)
  )
  expect_true(all(is.finite(result$SS)))
})

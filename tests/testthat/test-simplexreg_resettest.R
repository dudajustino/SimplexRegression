library(testthat)
library(SimplexRegression)

# ---------------------------------------------------------------------------
# Fixture
# ---------------------------------------------------------------------------
set.seed(2026)
n <- 100
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
z1 <- runif(n, 0, 1)
mu     <- parametric_mean_link_inv(0.6 - 2*x1 - 1.5*x2, 0.5, "plogit1")
sigma2 <- dispersion_link_inv(-2 - 2.5*z1, "log")
y      <- rsimplex(n, mu, sigma2)
data   <- data.frame(y = y, x1 = x1, x2 = x2, z1 = z1)

fit_plogit1 <- simplexreg(y ~ x1 + x2 | z1, data = data, link.mu = "plogit1")
fit_logit   <- simplexreg(y ~ x1 + x2 | z1, data = data, link.mu = "logit")

# Pre-compute results reused across tests
res_default  <- resettest(fit_plogit1)
res_nodisp   <- resettest(fit_plogit1, dispersion = FALSE)
res_power23  <- resettest(fit_plogit1, power = 2:3)
res_fitted   <- resettest(fit_plogit1, type = "fitted")
res_scalar   <- resettest(fit_plogit1, power = 2)   # scalar power → triggers is.vector branch

# ---------------------------------------------------------------------------

test_that("resettest requires an object of class simplexregression", {
  expect_error(resettest(list()),   "class 'simplexregression'")
  expect_error(resettest("model"),  "class 'simplexregression'")
})

test_that("resettest rejects power < 2", {
  expect_error(resettest(fit_plogit1, power = 1),    "integers >= 2")
  expect_error(resettest(fit_plogit1, power = 0),    "integers >= 2")
  expect_error(resettest(fit_plogit1, power = c(2, 1)), "integers >= 2")
})

test_that("resettest rejects non-integer power", {
  expect_error(resettest(fit_plogit1, power = 2.5),  "integers >= 2")
  expect_error(resettest(fit_plogit1, power = "2"),  "integers >= 2")
})

test_that("resettest returns an object of class htest", {
  expect_s3_class(res_default, "htest")
})

test_that("resettest returns all required htest components", {
  expect_named(res_default,
               c("statistic", "parameter", "p.value", "method", "data.name"),
               ignore.order = TRUE)
})

test_that("resettest method description is correct", {
  expect_equal(res_default$method, "RESET test")
})

test_that("resettest returns a finite non-negative test statistic", {
  expect_true(is.finite(res_default$statistic))
  expect_true(res_default$statistic >= 0)
})

test_that("resettest returns a p-value in [0, 1]", {
  expect_true(res_default$p.value >= 0 && res_default$p.value <= 1)
})

test_that("resettest degrees of freedom is positive", {
  expect_true(res_default$parameter > 0)
})

test_that("resettest works with dispersion = FALSE", {
  expect_s3_class(res_nodisp, "htest")
  expect_true(is.finite(res_nodisp$statistic))
  expect_true(res_nodisp$p.value >= 0 && res_nodisp$p.value <= 1)
})

test_that("resettest dispersion = TRUE has >= df than dispersion = FALSE", {
  expect_true(res_default$parameter >= res_nodisp$parameter)
})

test_that("resettest power = 2:3 returns larger df than power = 2", {
  expect_true(res_power23$parameter > res_scalar$parameter)
})

test_that("resettest power = 2:3 returns finite statistic and valid p-value", {
  expect_true(is.finite(res_power23$statistic))
  expect_true(res_power23$p.value >= 0 && res_power23$p.value <= 1)
})

test_that("resettest type = 'fitted' returns valid htest", {
  expect_s3_class(res_fitted, "htest")
  expect_true(is.finite(res_fitted$statistic))
  expect_true(res_fitted$p.value >= 0 && res_fitted$p.value <= 1)
})

test_that("resettest type = 'lp' and type = 'fitted' give different statistics", {
  expect_false(identical(res_default$statistic, res_fitted$statistic))
})

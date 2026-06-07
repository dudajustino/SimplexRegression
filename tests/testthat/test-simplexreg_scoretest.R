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

fit_logit   <- simplexreg(y ~ x1 + x2 | z1, data = data,
                           link.mu = "logit")
fit_plogit1 <- simplexreg(y ~ x1 + x2 | z1, data = data, link.mu = "plogit1")
fit_plogit2 <- simplexreg(y ~ x1 + x2 | z1, data = data, link.mu = "plogit2")

# ---------------------------------------------------------------------------

test_that("scoretest requires an object of class simplexregression", {
  expect_error(scoretest(list()),   "class 'simplexregression'")
  expect_error(scoretest("model"),  "class 'simplexregression'")
})

test_that("scoretest requires logit link in the model under H0", {
  expect_error(
    scoretest(fit_plogit1, link.mu = "plogit2"),
    "only valid for models fitted with the logit mean link"
  )
  expect_error(
    scoretest(fit_plogit2, link.mu = "plogit1"),
    "only valid for models fitted with the logit mean link"
  )
})

test_that("scoretest returns an object of class htest", {
  result <- scoretest(fit_logit, link.mu = "plogit1")
  expect_s3_class(result, "htest")
})

test_that("scoretest returns all required htest components", {
  result <- scoretest(fit_logit, link.mu = "plogit1")
  expect_named(result, c("statistic", "parameter", "p.value", "method", "data.name"),
               ignore.order = TRUE)
})

test_that("scoretest returns a finite and non-negative test statistic", {
  result1 <- scoretest(fit_logit, link.mu = "plogit1")
  result2 <- scoretest(fit_logit, link.mu = "plogit2")

  expect_true(is.finite(result1$statistic))
  expect_true(result1$statistic >= 0)

  expect_true(is.finite(result2$statistic))
  expect_true(result2$statistic >= 0)
})

test_that("scoretest returns a p-value between 0 and 1", {
  result1 <- scoretest(fit_logit, link.mu = "plogit1")
  result2 <- scoretest(fit_logit, link.mu = "plogit2")

  expect_true(result1$p.value >= 0 && result1$p.value <= 1)
  expect_true(result2$p.value >= 0 && result2$p.value <= 1)
})

test_that("scoretest has 1 degree of freedom", {
  result <- scoretest(fit_logit, link.mu = "plogit1")
  expect_equal(unname(result$parameter), 1L)
})

test_that("scoretest data.name contains the alternative link function", {
  result1 <- scoretest(fit_logit, link.mu = "plogit1")
  result2 <- scoretest(fit_logit, link.mu = "plogit2")

  expect_true(grepl("plogit1", result1$data.name))
  expect_true(grepl("plogit2", result2$data.name))
})

test_that("scoretest produces a finite p-value when data come from a non-logit link", {
  # Data generated with plogit1 — the test should detect departure from logit
  result <- scoretest(fit_logit, link.mu = "plogit1")
  expect_true(is.finite(result$p.value))
})

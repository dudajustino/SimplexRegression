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
fit_logit <- simplexreg(y ~ x1 + x2 | z1, data = data, link.mu = "logit")

all_types <- c("quantile", "pearson", "deviance", "standardized",
               "weighted", "variance", "biasvariance", "score",
               "dualscore", "response")

# ---------------------------------------------------------------------------

test_that("residuals.simplexregression requires an object of class simplexregression", {
  # After adding inherits() check to the function, both lines work correctly
  expect_error(residuals.simplexregression(list()),  "class 'simplexregression'")
  expect_error(residuals.simplexregression("model"), "class 'simplexregression'")
})

test_that("residuals.simplexregression rejects invalid type", {
  expect_error(residuals(fit_logit, type = "invalid"))
})

test_that("residuals.simplexregression returns a numeric vector of length n for all types", {
  for (type in all_types) {
    res <- residuals(fit_plogit1, type = type)
    expect_true(is.numeric(res),   label = paste(type, "is numeric"))
    expect_equal(length(res), n,   label = paste(type, "has length n"))
  }
})

test_that("residuals.simplexregression returns finite values for all types", {
  for (type in all_types) {
    res <- residuals(fit_plogit1, type = type)
    expect_true(all(is.finite(res)), label = paste(type, "all finite"))
  }
})

test_that("residuals.simplexregression default type is quantile", {
  res_default  <- residuals(fit_plogit1)
  res_quantile <- residuals(fit_plogit1, type = "quantile")
  expect_equal(res_default, res_quantile)
})

test_that("residuals.simplexregression quantile residuals are approximately standard normal", {
  res <- residuals(fit_plogit1, type = "quantile")
  expect_true(abs(mean(res)) < 0.5)
  expect_true(abs(sd(res) - 1) < 0.5)
})

test_that("residuals.simplexregression response residuals equal y - fitted", {
  res  <- residuals(fit_plogit1, type = "response")
  diff <- as.vector(fit_plogit1$y) - as.vector(fit_plogit1$mu.fv)
  expect_equal(res, diff)
})

test_that("residuals.simplexregression pearson residuals have mean close to 0", {
  res <- residuals(fit_plogit1, type = "pearson")
  expect_true(abs(mean(res)) < 1)
})

test_that("residuals.simplexregression works with fixed link model", {
  for (type in all_types) {
    res <- residuals(fit_logit, type = type)
    expect_true(is.numeric(res),     label = paste("logit -", type, "is numeric"))
    expect_equal(length(res), n,     label = paste("logit -", type, "has length n"))
    expect_true(all(is.finite(res)), label = paste("logit -", type, "all finite"))
  }
})

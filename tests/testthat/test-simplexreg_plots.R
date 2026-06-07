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
# plot.simplexregression
# ---------------------------------------------------------------------------

test_that("plot.simplexregression rejects invalid which values", {
  expect_error(plot(fit_plogit1, which = 0))
  expect_error(plot(fit_plogit1, which = 8))
})

test_that("plot.simplexregression runs without error for all plots (fixed link)", {
  expect_no_error({
    pdf(NULL)
    plot(fit_logit, which = 1:7, ask = FALSE)
    dev.off()
  })
})

test_that("plot.simplexregression runs without error for all plots (parametric link)", {
  expect_no_error({
    pdf(NULL)
    plot(fit_plogit1, which = 1:7, ask = FALSE)
    dev.off()
  })
})

test_that("plot.simplexregression works for each individual plot index", {
  for (i in 1:7) {
    expect_error(
      { pdf(NULL); plot(fit_plogit1, which = i, ask = FALSE); dev.off() },
      NA,
      info = paste("which =", i)
    )
  }
})

test_that("plot.simplexregression works with different residual types", {
  for (type in c("quantile", "pearson", "deviance", "response")) {
    expect_error(
      { pdf(NULL); plot(fit_plogit1, which = 1, type = type, ask = FALSE); dev.off() },
      NA,
      info = paste("type =", type)
    )
  }
})

test_that("plot.simplexregression returns invisibly", {
  pdf(NULL)
  result <- withVisible(plot(fit_plogit1, which = 1, ask = FALSE))
  dev.off()
  expect_false(result$visible)
})

test_that("plot.simplexregression threshold flags observations for Cook's distance (which = 6)", {
  # threshold = 0 garantees all observations exceed it, triggering the text() branch
  expect_no_error({
    pdf(NULL)
    plot(fit_plogit1, which = 6, ask = FALSE, threshold = 0.01)
    dev.off()
  })
})

test_that("plot.simplexregression threshold flags observations for gleverage (which = 7)", {
  expect_no_error({
    pdf(NULL)
    plot(fit_plogit1, which = 7, ask = FALSE, threshold = 0.1)
    dev.off()
  })
})

# ---------------------------------------------------------------------------
# halfnormal.plot
# ---------------------------------------------------------------------------

test_that("halfnormal.plot requires an object of class simplexregression", {
  expect_error(halfnormal.plot(list()),  "class 'simplexregression'")
  expect_error(halfnormal.plot("model"), "class 'simplexregression'")
})

test_that("halfnormal.plot runs without error with default arguments", {
  expect_no_error({
    pdf(NULL)
    halfnormal.plot(fit_plogit1, nsim = 10, seed = 42)
    dev.off()
  })
})

test_that("halfnormal.plot runs without error with parametric link", {
  expect_no_error({
    pdf(NULL)
    halfnormal.plot(fit_logit, nsim = 10, seed = 42)
    dev.off()
  })
})

test_that("halfnormal.plot works with different residual types", {
  for (type in c("weighted", "quantile", "pearson")) {
    expect_error(
      { pdf(NULL); halfnormal.plot(fit_plogit1, type = type, nsim = 10, seed = 42); dev.off() },
      NA,
      info = paste("type =", type)
    )
  }
})

library(testthat)
library(SimplexRegression)

# ==============================================================================
# test-simplexreg_influence.R
# Tests for: local.influence, gleverage, diag.im, diag.distances
# ==============================================================================

# ---------------------------------------------------------------------------
# Fixture: simulate data and fit models once for all tests
# ---------------------------------------------------------------------------
set.seed(2026)
n <- 50
x1 <- runif(n, 0, 1)
z1 <- runif(n, 0, 1)
mu <- fixed_mean_link_inv(1 - 1.5*x1, "loglog")
sigma2 <- dispersion_link_inv(-2 - 2.5*z1, "log")
y <- rsimplex(n, mu, sigma2)
data <- data.frame(y = y, x1 = x1, z1 = z1)

fit_loglog <- simplexreg(y ~ x1 | z1, data = data, link.mu = "loglog")
fit_logit  <- simplexreg(y ~ x1 | z1, data = data, link.mu = "logit")

# ------------------------------------------------------------------------------
# diag.im fixtures: computed once, reused across all diag.im tests
# Each unique (type x parameter x interval) combination = one LOO run
# ------------------------------------------------------------------------------
im_s3        <- diag.im(fit_loglog, data = data, type = "s3",
                        parameter = "theta", interval = "I1", verbose = FALSE)
im_s5        <- diag.im(fit_loglog, data = data, type = "s5",
                        parameter = "theta", interval = "I1", verbose = FALSE)
im_both      <- diag.im(fit_loglog, data = data, type = c("s3", "s5"),
                        parameter = "theta", interval = "I1", verbose = FALSE)
im_s3_beta   <- diag.im(fit_loglog, data = data, type = "s3",
                        parameter = "beta",  interval = "I1", verbose = FALSE)
im_s3_gamma  <- diag.im(fit_loglog, data = data, type = "s3",
                        parameter = "gamma", interval = "I1", verbose = FALSE)

# ==============================================================================
# 1. local.influence
# ==============================================================================

test_that("local.influence requires an object of class simplexregression", {
  expect_error(local.influence(list()),   "class 'simplexregression'")
  expect_error(local.influence("model"),  "class 'simplexregression'")
})

test_that("local.influence returns a list with all named components", {
  result <- local.influence(fit_loglog, scheme = "case.weight")
  expect_type(result, "list")
  expect_named(result,
               c("dmax.beta", "dmax.gamma", "dmax.theta",
                 "Ci.beta",   "Ci.gamma",   "Ci.theta"),
               ignore.order = TRUE)
})

test_that("local.influence all components have length n", {
  result <- local.influence(fit_loglog, scheme = "case.weight")
  for (nm in names(result)) {
    expect_length(result[[nm]], n)
  }
})

test_that("local.influence all components are finite numeric vectors", {
  result <- local.influence(fit_loglog, scheme = "case.weight")
  for (nm in names(result)) {
    expect_true(is.numeric(result[[nm]]),
                label = paste(nm, "is numeric"))
    expect_true(all(is.finite(result[[nm]])),
                label = paste(nm, "all finite"))
  }
})

test_that("local.influence Ci components are non-negative", {
  result <- local.influence(fit_loglog, scheme = "case.weight")
  expect_true(all(result$Ci.beta  >= 0))
  expect_true(all(result$Ci.gamma >= 0))
  expect_true(all(result$Ci.theta >= 0))
})

test_that("local.influence dmax components are non-negative", {
  result <- local.influence(fit_loglog, scheme = "case.weight")
  expect_true(all(result$dmax.beta  >= 0))
  expect_true(all(result$dmax.gamma >= 0))
  expect_true(all(result$dmax.theta >= 0))
})

test_that("local.influence works with response perturbation scheme", {
  result <- local.influence(fit_loglog, scheme = "response")
  expect_type(result, "list")
  expect_true(all(is.finite(result$Ci.theta)))
})

test_that("local.influence case.weight and response give different results", {
  cw <- local.influence(fit_loglog, scheme = "case.weight")
  rp <- local.influence(fit_loglog, scheme = "response")
  expect_false(identical(cw$Ci.theta, rp$Ci.theta))
})

test_that("local.influence plot = TRUE returns result invisibly", {
  pdf(NULL)
  result <- withVisible(
    local.influence(fit_loglog, scheme = "case.weight",
                    parameter = "beta", type = "Ci", plot = TRUE)
  )
  dev.off()
  expect_false(result$visible)
  expect_type(result$value, "list")
})

test_that("local.influence plot = TRUE runs without error for all parameter/type combos", {
  params <- c("theta", "beta", "gamma")
  types  <- c("Ci", "dmax")
  for (par in params) {
    for (tp in types) {
      expect_no_error({
        pdf(NULL)
        local.influence(fit_loglog, scheme = "case.weight",
                        parameter = par, type = tp, plot = TRUE)
        dev.off()
      })
    }
  }
})

test_that("local.influence plot with threshold runs without error", {
  expect_no_error({
    pdf(NULL)
    local.influence(fit_logit, scheme = "case.weight",
                    parameter = "beta", type = "Ci",
                    plot = TRUE, threshold = 0.5)
    dev.off()
  })
})

test_that("local.influence both schemes give consistent structure", {
  cw <- local.influence(fit_logit, scheme = "case.weight")
  rp <- local.influence(fit_logit, scheme = "response")
  expect_equal(names(cw), names(rp))
})

# ==============================================================================
# 2. gleverage
# ==============================================================================

test_that("gleverage.simplexregression returns numeric vector of length n", {
  glev <- gleverage(fit_loglog)
  expect_true(is.numeric(glev))
  expect_length(glev, n)
})

test_that("gleverage returns finite values", {
  glev <- gleverage(fit_loglog)
  expect_true(all(is.finite(glev)))
})

test_that("gleverage values are non-negative", {
  glev <- gleverage(fit_loglog)
  expect_true(all(glev >= 0))
})

test_that("gleverage gives same result on repeated calls (deterministic)", {
  expect_equal(gleverage(fit_loglog), gleverage(fit_loglog))
})

# ==============================================================================
# 3. diag.im
# Fixtures im_s3, im_s5, im_both, im_s3_I2, im_s3_beta, im_s3_gamma
# are pre-computed above — no LOO runs happen inside these tests.
# ==============================================================================

test_that("diag.im requires an object of class simplexregression", {
  expect_error(diag.im(list(),  data = data), "class 'simplexregression'")
  expect_error(diag.im("model", data = data), "class 'simplexregression'")
})

test_that("diag.im rejects invalid ncores", {
  expect_error(diag.im(fit_logit, data = data, ncores = 0),   "positive integer")
  expect_error(diag.im(fit_logit, data = data, ncores = -1),  "positive integer")
  expect_error(diag.im(fit_logit, data = data, ncores = 1.5), "positive integer")
})

test_that("diag.im parameter = 'gamma' errors for intercept-only dispersion", {
  fit_nodispers <- simplexreg(y ~ x1, data = data, link.mu = "logit")
  expect_error(
    diag.im(fit_nodispers, data = data, type = "s3", parameter = "gamma",
            verbose = FALSE),
    "gamma"
  )
})

test_that("diag.im returns list with s3 components when type = 's3'", {
  expect_type(im_s3, "list")
  expect_true("s3_i"        %in% names(im_s3))
  expect_true("limits_s3"   %in% names(im_s3))
  expect_true("outliers_s3" %in% names(im_s3))
  expect_false("s5_i" %in% names(im_s3))
})

test_that("diag.im returns list with s5 components when type = 's5'", {
  expect_type(im_s5, "list")
  expect_true("s5_i"        %in% names(im_s5))
  expect_true("limits_s5"   %in% names(im_s5))
  expect_true("outliers_s5" %in% names(im_s5))
  expect_false("s3_i" %in% names(im_s5))
})

test_that("diag.im returns both s3 and s5 when both requested", {
  expect_true("s3_i" %in% names(im_both))
  expect_true("s5_i" %in% names(im_both))
})

test_that("diag.im s3_i has length n and is finite numeric", {
  expect_length(im_s3$s3_i, n)
  expect_true(is.numeric(im_s3$s3_i))
  expect_true(all(is.finite(im_s3$s3_i)))
})

test_that("diag.im s5_i has length n and is finite numeric", {
  expect_length(im_s5$s5_i, n)
  expect_true(is.numeric(im_s5$s5_i))
  expect_true(all(is.finite(im_s5$s5_i)))
})

test_that("diag.im s3_i values are positive", {
  expect_true(all(im_s3$s3_i > 0))
})

test_that("diag.im limits_s3 has lower and upper named elements", {
  expect_named(im_s3$limits_s3, c("lower", "upper"), ignore.order = TRUE)
  expect_true(im_s3$limits_s3["lower"] < im_s3$limits_s3["upper"])
})

test_that("diag.im outliers_s3 is a data.frame", {
  expect_s3_class(im_s3$outliers_s3, "data.frame")
})

test_that("diag.im outliers_s3 obs indices are within 1:n", {
  if (nrow(im_s3$outliers_s3) > 0)
    expect_true(all(im_s3$outliers_s3$Obs %in% seq_len(n)))
})

test_that("diag.im stores interval, parameter, and n in result", {
  expect_equal(im_s3$interval,  "I1")
  expect_equal(im_s3$parameter, "theta")
  expect_equal(im_s3$n,         n)
})

test_that("diag.im parameter = 'beta' returns finite s3_i", {
  expect_true(all(is.finite(im_s3_beta$s3_i)))
})

test_that("diag.im parameter = 'gamma' returns finite s3_i", {
  expect_true(all(is.finite(im_s3_gamma$s3_i)))
})

test_that("diag.im plot = TRUE on pre-computed result returns invisibly and plots without error", {
  pdf(NULL)
  # Re-run with plot = TRUE reusing same config as im_both fixture
  result <- withVisible(
    diag.im(fit_logit, data = data, type = c("s3", "s5"),
            parameter = "theta", verbose = FALSE, plot = TRUE)
  )
  dev.off()
  expect_false(result$visible)
  expect_type(result$value, "list")
})

# ==============================================================================
# 4. diag.distances
# H is the fastest type — used as main fixture for all structural tests.
# W1 and W2 get minimal smoke tests (type label only) to avoid LOO + quadgk
# overhead from slower numerical integration.
# ==============================================================================

dd_H <- diag.distances(fit_logit, data = data, type = "H", verbose = FALSE)
dd_W1 <- diag.distances(fit_logit, data = data, type = "W1", verbose = FALSE)

test_that("diag.distances requires an object of class simplexregression", {
  expect_error(diag.distances(list(),  data = data), "class 'simplexregression'")
  expect_error(diag.distances("model", data = data), "class 'simplexregression'")
})

test_that("diag.distances rejects invalid ncores", {
  expect_error(diag.distances(fit_logit, data = data, ncores = 0),   "positive integer")
  expect_error(diag.distances(fit_logit, data = data, ncores = -1),  "positive integer")
  expect_error(diag.distances(fit_logit, data = data, ncores = 1.5), "positive integer")
})

test_that("diag.distances rejects invalid type string", {
  expect_error(diag.distances(fit_logit, data = data, type = "X", verbose = FALSE))
})

test_that("diag.distances returns a list with required components", {
  expect_type(dd_H, "list")
  expect_named(dd_H, c("distances", "threshold", "outliers", "type", "n"),
               ignore.order = TRUE)
})

test_that("diag.distances distances vector has length n", {
  expect_length(dd_H$distances, n)
})

test_that("diag.distances distances are non-negative and finite", {
  expect_true(all(dd_H$distances >= 0))
  expect_true(all(is.finite(dd_H$distances)))
})

test_that("diag.distances threshold is a positive named scalar", {
  expect_length(dd_H$threshold, 1L)
  expect_named(dd_H$threshold, "upper")
  expect_true(dd_H$threshold > 0)
})

test_that("diag.distances outliers is a data.frame", {
  expect_s3_class(dd_H$outliers, "data.frame")
})

test_that("diag.distances flagged observations exceed threshold", {
  if (nrow(dd_H$outliers) > 0) {
    expect_true(all(dd_H$outliers$Obs %in% seq_len(n)))
    expect_true(all(dd_H$outliers$distance > dd_H$threshold["upper"]))
  }
})

test_that("diag.distances n field equals number of observations", {
  expect_equal(dd_H$n, n)
})

test_that("diag.distances type label is correct for H", {
  expect_equal(dd_H$type, "Hellinger")
})

test_that("diag.distances plot = TRUE returns list invisibly", {
  pdf(NULL)
  result <- withVisible(
    diag.distances(fit_logit, data = data, type = "H",
                   verbose = FALSE, plot = TRUE)
  )
  dev.off()
  expect_false(result$visible)
  expect_type(result$value, "list")
})

test_that("diag.distances W1 returns correct type label", {
  expect_equal(dd_W1$type, "Wasserstein-1")
})

test_that("diag.distances W1 distances are non-negative and finite", {
  expect_true(all(dd_W1$distances >= 0))
  expect_true(all(is.finite(dd_W1$distances)))
})

test_that("diag.distances H and W1 give different distances", {
  expect_false(identical(dd_H$distances, dd_W1$distances))
})

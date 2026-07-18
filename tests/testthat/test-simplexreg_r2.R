library(testthat)
library(SimplexRegression)

# ==============================================================================
# test-simplexreg_r2.R
# Tests for: r2
# ==============================================================================

# ---------------------------------------------------------------------------
# Fixture
# ---------------------------------------------------------------------------
set.seed(42)
n   <- 100
x1  <- runif(n, 0, 1)
x2  <- runif(n, 0, 1)
mu  <- parametric_mean_link_inv(0.8 - 1.2*x1 - 1.5*x2, 0.25, "plogit2")
y   <- rsimplex(n, mu, 0.5)
dat <- data.frame(y = y, x1 = x1, x2 = x2)

fit_plogit2 <- simplexreg(y ~ x1 + x2 | 1, data = dat, link.mu = "plogit2")
fit_logit   <- simplexreg(y ~ x1 + x2 | 1, data = dat, link.mu = "logit")
fit_vardisp <- simplexreg(y ~ x1 + x2 | x1, data = dat, link.mu = "logit")

# ==============================================================================
# 1. Input validation
# ==============================================================================

test_that("r2 rejects non-simplexregression objects", {
  expect_error(r2(list()),    "class 'simplexregression'")
  expect_error(r2(fit_logit, list()), "class 'simplexregression'")
})

test_that("r2 rejects alpha1 outside [0, 1]", {
  expect_error(r2(fit_logit, alpha1 = -0.1), "alpha1")
  expect_error(r2(fit_logit, alpha1 =  1.1), "alpha1")
})

test_that("r2 rejects non-positive alpha2", {
  expect_error(r2(fit_logit, alpha2 = 0),    "alpha2")
  expect_error(r2(fit_logit, alpha2 = -1),   "alpha2")
})

# ==============================================================================
# 2. Single model — return structure
# ==============================================================================

test_that("r2 returns a named numeric vector for a single model", {
  res <- r2(fit_logit)
  expect_true(is.numeric(res))
  expect_named(res, c("R2_FC", "R2_FC_c", "R2_N", "R2_N_c",
                      "R2_Nw_c", "R2", "R2_HS"))
})

test_that("r2 returns a vector of length 7 for a single model", {
  expect_length(r2(fit_logit), 7L)
})

test_that("r2 all components are finite for parametric link", {
  res <- r2(fit_plogit2)
  expect_true(all(is.finite(res)))
})

test_that("r2 all components are finite for fixed link", {
  res <- r2(fit_logit)
  expect_true(all(is.finite(res)))
})

# ==============================================================================
# 3. Single model — value ranges
# ==============================================================================

test_that("r2 R2_FC is in [0, 1]", {
  res <- r2(fit_logit)
  expect_true(res["R2_FC"] >= 0 & res["R2_FC"] <= 1)
})

test_that("r2 R2_N is in [0, 1]", {
  res <- r2(fit_logit)
  expect_true(res["R2_N"] >= 0 & res["R2_N"] <= 1)
})

test_that("r2 R2 (conventional) is at most 1", {
  res <- r2(fit_logit)
  expect_true(res["R2"] <= 1)
})

test_that("r2 R2_HS is at most 1", {
  res <- r2(fit_logit)
  expect_true(res["R2_HS"] <= 1)
})

# ==============================================================================
# 4. Single model — consistency with fitted object
# ==============================================================================

test_that("r2 R2_FC matches value stored in fitted object", {
  res <- r2(fit_logit)
  expect_equal(unname(res["R2_FC"]), fit_logit$R2_FC, tolerance = 1e-10)
})

test_that("r2 R2_N matches value stored in fitted object", {
  res <- r2(fit_logit)
  expect_equal(unname(res["R2_N"]), fit_logit$R2_N, tolerance = 1e-10)
})

test_that("r2 R2 equals 1 minus ratio of SSres to SStot", {
  res  <- r2(fit_logit)
  y    <- fit_logit$y
  mu   <- fit_logit$mu.fv
  ref  <- 1 - sum((y - mu)^2) / sum((y - mean(y))^2)
  expect_equal(unname(res["R2"]), ref, tolerance = 1e-10)
})

# ==============================================================================
# 5. Corrections — direction and alpha sensibility
# ==============================================================================

test_that("r2 R2_FC_c is a finite-sample correction of R2_FC", {
  res <- r2(fit_logit)
  # correction moves R2_FC toward 1 when n >> k (both should be close)
  expect_true(is.finite(res["R2_FC_c"]))
})

test_that("r2 R2_N_c is close to R2_Nw_c when alpha1 = 0 and alpha2 = 1", {
  res_default <- r2(fit_logit, alpha1 = 0, alpha2 = 1)
  # With alpha1 = 0 and alpha2 = 1, R2_Nw_c reduces to R2_N_c
  expect_equal(unname(res_default["R2_N_c"]),
               unname(res_default["R2_Nw_c"]),
               tolerance = 1e-10)
})

test_that("r2 R2_Nw_c varies with alpha1", {
  res1 <- r2(fit_vardisp, alpha1 = 0.0)
  res2 <- r2(fit_vardisp, alpha1 = 0.8)
  # Different alpha1 should produce different corrections when p != q
  expect_false(isTRUE(all.equal(res1["R2_Nw_c"], res2["R2_Nw_c"])))
})

test_that("r2 R2_Nw_c varies with alpha2", {
  res1 <- r2(fit_logit, alpha2 = 1)
  res2 <- r2(fit_logit, alpha2 = 2)
  expect_false(isTRUE(all.equal(res1["R2_Nw_c"], res2["R2_Nw_c"])))
})

test_that("r2 alpha3 = 'log' gives different R2_HS than alpha3 = '1'", {
  res1 <- r2(fit_logit, alpha3 = "1")
  res2 <- r2(fit_logit, alpha3 = "log")
  expect_false(isTRUE(all.equal(res1["R2_HS"], res2["R2_HS"])))
})

# ==============================================================================
# 6. Multiple models — return structure
# ==============================================================================

test_that("r2 returns a data.frame for multiple models", {
  res <- r2(fit_plogit2, fit_logit)
  expect_s3_class(res, "data.frame")
})

test_that("r2 data.frame has correct dimensions for two models", {
  res <- r2(fit_plogit2, fit_logit)
  expect_equal(nrow(res), 2L)
  expect_equal(ncol(res), 7L)
})

test_that("r2 data.frame has correct column names", {
  res <- r2(fit_plogit2, fit_logit)
  expect_named(res, c("R2_FC", "R2_FC_c", "R2_N", "R2_N_c",
                      "R2_Nw_c", "R2", "R2_HS"))
})

test_that("r2 data.frame row values match single-model calls", {
  res_multi  <- r2(fit_plogit2, fit_logit)
  res_single <- r2(fit_logit)
  expect_equal(as.numeric(res_multi[2, ]), as.numeric(res_single),
               tolerance = 1e-10)
})

test_that("r2 data.frame all values are finite", {
  res <- r2(fit_plogit2, fit_logit)
  expect_true(all(is.finite(unlist(res))))
})

test_that("r2 works for three models", {
  res <- r2(fit_plogit2, fit_logit, fit_vardisp)
  expect_equal(nrow(res), 3L)
})

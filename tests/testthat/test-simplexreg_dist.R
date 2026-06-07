library(testthat)
library(SimplexRegression)

# ==============================================================================
# test-simplexreg_dist.R
# Tests for: dsimplex, psimplex, qsimplex, rsimplex
# ==============================================================================

# ==============================================================================
# 1. DENSITY FUNCTION — dsimplex
# ==============================================================================

test_that("dsimplex rejects invalid mu", {
  expect_error(dsimplex(0.5, mu = 0,    sigma2 = 0.5), "'mu' must be in \\(0, 1\\)")
  expect_error(dsimplex(0.5, mu = 1,    sigma2 = 0.5), "'mu' must be in \\(0, 1\\)")
  expect_error(dsimplex(0.5, mu = -0.1, sigma2 = 0.5), "'mu' must be in \\(0, 1\\)")
  expect_error(dsimplex(0.5, mu = 1.5,  sigma2 = 0.5), "'mu' must be in \\(0, 1\\)")
})

test_that("dsimplex rejects non-positive sigma2", {
  expect_error(dsimplex(0.5, mu = 0.3, sigma2 = 0),    "'sigma2' must be positive")
  expect_error(dsimplex(0.5, mu = 0.3, sigma2 = -0.1), "'sigma2' must be positive")
})

test_that("dsimplex returns non-negative values for x in (0, 1)", {
  x <- seq(0.1, 0.9, by = 0.1)
  res <- dsimplex(x, mu = 0.4, sigma2 = 0.5)
  expect_true(all(res >= 0))
})

test_that("dsimplex returns a numeric vector of the correct length", {
  x   <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  res <- dsimplex(x, mu = 0.5, sigma2 = 0.5)
  expect_true(is.numeric(res))
  expect_length(res, 5)
})

test_that("dsimplex with log = TRUE equals log of density", {
  x  <- c(0.2, 0.5, 0.8)
  d  <- dsimplex(x, mu = 0.4, sigma2 = 0.3)
  ld <- dsimplex(x, mu = 0.4, sigma2 = 0.3, log = TRUE)
  expect_equal(ld, log(d), tolerance = 1e-10)
})

test_that("dsimplex integrates to approximately 1", {
  result <- integrate(dsimplex, lower = 1e-8, upper = 1 - 1e-8,
                      mu = 0.4, sigma2 = 0.5)$value
  expect_equal(result, 1, tolerance = 1e-4)
})

test_that("dsimplex integrates to approximately 1 for various parameters", {
  params <- list(
    list(mu = 0.2, sigma2 = 0.3),
    list(mu = 0.5, sigma2 = 1.0),
    list(mu = 0.8, sigma2 = 0.5)
  )
  for (p in params) {
    result <- integrate(dsimplex, lower = 1e-8, upper = 1 - 1e-8,
                        mu = p$mu, sigma2 = p$sigma2)$value
    expect_equal(result, 1, tolerance = 5e-4,
                 label = paste("mu =", p$mu, "sigma2 =", p$sigma2))
  }
})

test_that("dsimplex is maximised near mu", {
  mu    <- 0.4
  x_seq <- seq(0.05, 0.95, by = 0.05)
  dens  <- dsimplex(x_seq, mu = mu, sigma2 = 0.3)
  peak  <- x_seq[which.max(dens)]
  expect_true(abs(peak - mu) < 0.15)
})

test_that("dsimplex is symmetric about mu = 0.5 when mu = 0.5", {
  x   <- c(0.2, 0.3, 0.4)
  d1  <- dsimplex(x, mu = 0.5, sigma2 = 0.5)
  d2  <- dsimplex(1 - x, mu = 0.5, sigma2 = 0.5)
  expect_equal(d1, d2, tolerance = 1e-10)
})

# ==============================================================================
# 2. CUMULATIVE DISTRIBUTION FUNCTION — psimplex
# ==============================================================================

test_that("psimplex rejects invalid mu", {
  expect_error(psimplex(0.5, mu = 0,    sigma2 = 0.5), "'mu' must be in \\(0, 1\\)")
  expect_error(psimplex(0.5, mu = 1,    sigma2 = 0.5), "'mu' must be in \\(0, 1\\)")
  expect_error(psimplex(0.5, mu = -0.2, sigma2 = 0.5), "'mu' must be in \\(0, 1\\)")
})

test_that("psimplex rejects non-positive sigma2", {
  expect_error(psimplex(0.5, mu = 0.3, sigma2 = 0),    "'sigma2' must be positive")
  expect_error(psimplex(0.5, mu = 0.3, sigma2 = -1.0), "'sigma2' must be positive")
})

test_that("psimplex returns values in [0, 1]", {
  q   <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  res <- psimplex(q, mu = 0.4, sigma2 = 0.5)
  expect_true(all(res >= 0 & res <= 1))
})

test_that("psimplex is non-decreasing", {
  q   <- seq(0.1, 0.9, by = 0.1)
  res <- psimplex(q, mu = 0.4, sigma2 = 0.5)
  expect_true(all(diff(res) >= 0))
})

test_that("psimplex(mu) is approximately 0.5", {
  mu  <- 0.5
  res <- psimplex(mu, mu = mu, sigma2 = 0.5)
  expect_equal(res, 0.5, tolerance = 0.05)
})

test_that("psimplex lower.tail = FALSE equals 1 - lower.tail = TRUE", {
  q <- c(0.2, 0.5, 0.8)
  p_lo <- psimplex(q, mu = 0.4, sigma2 = 0.5, lower.tail = TRUE)
  p_up <- psimplex(q, mu = 0.4, sigma2 = 0.5, lower.tail = FALSE)
  expect_equal(p_lo + p_up, rep(1, 3), tolerance = 1e-8)
})

test_that("psimplex log.p = TRUE equals log of probability", {
  q  <- c(0.3, 0.5, 0.7)
  p  <- psimplex(q, mu = 0.4, sigma2 = 0.5)
  lp <- psimplex(q, mu = 0.4, sigma2 = 0.5, log.p = TRUE)
  expect_equal(lp, log(p), tolerance = 1e-8)
})

test_that("psimplex uses normal approximation for small sigma2 without error", {
  expect_no_error(psimplex(0.5, mu = 0.4, sigma2 = 1e-5))
})

# ==============================================================================
# 3. QUANTILE FUNCTION — qsimplex
# ==============================================================================

test_that("qsimplex rejects invalid mu", {
  expect_error(qsimplex(0.5, mu = 0,    sigma2 = 0.5), "'mu' must be in \\(0, 1\\)")
  expect_error(qsimplex(0.5, mu = 1,    sigma2 = 0.5), "'mu' must be in \\(0, 1\\)")
  expect_error(qsimplex(0.5, mu = -0.2, sigma2 = 0.5), "'mu' must be in \\(0, 1\\)")
})

test_that("qsimplex rejects non-positive sigma2", {
  expect_error(qsimplex(0.5, mu = 0.3, sigma2 = 0),    "'sigma2' must be positive")
  expect_error(qsimplex(0.5, mu = 0.3, sigma2 = -0.5), "'sigma2' must be positive")
})

test_that("qsimplex returns values in (0, 1)", {
  p   <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  res <- qsimplex(p, mu = 0.4, sigma2 = 0.5)
  expect_true(all(res > 0 & res < 1))
})

test_that("qsimplex returns non-decreasing values", {
  p   <- seq(0.1, 0.9, by = 0.1)
  res <- qsimplex(p, mu = 0.4, sigma2 = 0.5)
  expect_true(all(diff(res) >= 0))
})

test_that("qsimplex(0.5) is close to mu (median ≈ mean)", {
  mu  <- 0.5
  q50 <- qsimplex(0.5, mu = mu, sigma2 = 0.5)
  expect_equal(q50, mu, tolerance = 0.1)
})

test_that("qsimplex is the inverse of psimplex", {
  p    <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  q    <- qsimplex(p, mu = 0.4, sigma2 = 0.5)
  p_back <- psimplex(q, mu = 0.4, sigma2 = 0.5)
  expect_equal(p_back, p, tolerance = 1e-4)
})

test_that("qsimplex log.p = TRUE accepts log-scale probabilities", {
  p    <- c(0.25, 0.5, 0.75)
  q1   <- qsimplex(p,      mu = 0.4, sigma2 = 0.5, log.p = FALSE)
  q2   <- qsimplex(log(p), mu = 0.4, sigma2 = 0.5, log.p = TRUE)
  expect_equal(q1, q2, tolerance = 1e-6)
})

test_that("qsimplex lower.tail = FALSE is consistent with lower.tail = TRUE", {
  p  <- 0.3
  q1 <- qsimplex(p,       mu = 0.4, sigma2 = 0.5, lower.tail = TRUE)
  q2 <- qsimplex(1 - p,   mu = 0.4, sigma2 = 0.5, lower.tail = FALSE)
  expect_equal(q1, q2, tolerance = 1e-6)
})

test_that("qsimplex uses normal approximation for small sigma2 without error", {
  expect_no_error(qsimplex(0.5, mu = 0.4, sigma2 = 0.05))
})

test_that("qsimplex warns and caps sigma2 > 200", {
  expect_warning(qsimplex(0.5, mu = 0.4, sigma2 = 300), "sigma2 > 200")
})

test_that("qsimplex.norm recycles mu and sigma2 to match length of p", {
  p      <- c(0.1, 0.3, 0.5, 0.7, 0.9)   # length 5
  result <- qsimplex.norm(p, mu = 0.4, sigma2 = 0.5)  # mu e sigma2 escalares → reciclados
  expect_length(result, 5L)
  expect_true(all(is.finite(result)))
})

test_that("qsimplex.norm recycles sigma2 when mu is vector", {
  p      <- c(0.25, 0.5, 0.75)
  mu     <- c(0.3, 0.4, 0.5)              # length 3
  result <- qsimplex.norm(p, mu = mu, sigma2 = 0.5)   # sigma2 escalar → reciclado
  expect_length(result, 3L)
  expect_true(all(is.finite(result)))
})

# ==============================================================================
# 4. RANDOM GENERATION — rsimplex
# ==============================================================================

test_that("rsimplex rejects invalid mu", {
  expect_error(rsimplex(5, mu = 0,    sigma2 = 0.5), "'mu' must be in \\(0, 1\\)")
  expect_error(rsimplex(5, mu = 1,    sigma2 = 0.5), "'mu' must be in \\(0, 1\\)")
  expect_error(rsimplex(5, mu = -0.5, sigma2 = 0.5), "'mu' must be in \\(0, 1\\)")
  expect_error(rsimplex(5, mu = 1.5,  sigma2 = 0.5), "'mu' must be in \\(0, 1\\)")
})

test_that("rsimplex rejects non-positive sigma2", {
  expect_error(rsimplex(5, mu = 0.5, sigma2 = 0),    "'sigma2' must be positive")
  expect_error(rsimplex(5, mu = 0.5, sigma2 = -1.0), "'sigma2' must be positive")
})

test_that("rsimplex returns a numeric vector of length n", {
  res <- rsimplex(10, mu = 0.5, sigma2 = 0.5)
  expect_true(is.numeric(res))
  expect_length(res, 10)
})

test_that("rsimplex returns values in (0, 1)", {
  res <- rsimplex(200, mu = 0.5, sigma2 = 0.5)
  expect_true(all(res > 0 & res < 1))
})

test_that("rsimplex sample mean is close to mu for large n", {
  set.seed(1)
  mu  <- 0.4
  res <- rsimplex(5000, mu = mu, sigma2 = 0.5)
  expect_equal(mean(res), mu, tolerance = 0.03)
})

test_that("rsimplex accepts vector mu and sigma2 of length n", {
  n      <- 20
  mu_vec <- rep(c(0.3, 0.7), n / 2)
  s2_vec <- rep(c(0.3, 0.6), n / 2)
  res    <- rsimplex(n, mu = mu_vec, sigma2 = s2_vec)
  expect_length(res, n)
  expect_true(all(res > 0 & res < 1))
})

test_that("rsimplex is reproducible with the same seed", {
  set.seed(7)
  r1 <- rsimplex(10, mu = 0.5, sigma2 = 0.3)
  set.seed(7)
  r2 <- rsimplex(10, mu = 0.5, sigma2 = 0.3)
  expect_equal(r1, r2)
})

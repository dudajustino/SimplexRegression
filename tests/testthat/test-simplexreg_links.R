library(testthat)
library(SimplexRegression)

# ==============================================================================
# test-simplexreg_links.R
# Tests for: parametric_mean_links, fixed_mean_links, dispersion_links
# ==============================================================================

# ==============================================================================
# 1. PARAMETRIC MEAN LINK FUNCTIONS
# ==============================================================================

# --- parametric_mean_link -----------------------------------------------------

test_that("parametric_mean_link rejects mu outside (0, 1)", {
  expect_error(parametric_mean_link(0,    lambda = 1, type = "plogit2"), "interval \\(0, 1\\)")
  expect_error(parametric_mean_link(1,    lambda = 1, type = "plogit2"), "interval \\(0, 1\\)")
  expect_error(parametric_mean_link(-0.1, lambda = 1, type = "plogit1"), "interval \\(0, 1\\)")
  expect_error(parametric_mean_link(1.5,  lambda = 1, type = "plogit1"), "interval \\(0, 1\\)")
})

test_that("parametric_mean_link returns a finite numeric vector", {
  result <- parametric_mean_link(c(0.2, 0.5, 0.8), lambda = 1.5, type = "plogit2")
  expect_true(is.numeric(result))
  expect_length(result, 3)
  expect_true(all(is.finite(result)))
})

test_that("parametric_mean_link returns finite values for plogit1", {
  result <- parametric_mean_link(c(0.2, 0.5, 0.8), lambda = 1.5, type = "plogit1")
  expect_true(is.numeric(result))
  expect_length(result, 3)
  expect_true(all(is.finite(result)))
})

test_that("parametric_mean_link plogit2 equals logit when lambda = 1", {
  mu  <- c(0.2, 0.5, 0.8)
  res <- parametric_mean_link(mu, lambda = 1, type = "plogit2")
  ref <- log(mu / (1 - mu))
  expect_equal(res, ref, tolerance = 1e-8)
})

test_that("parametric_mean_link plogit1 equals logit when lambda = 1", {
  mu  <- c(0.2, 0.5, 0.8)
  res <- parametric_mean_link(mu, lambda = 1, type = "plogit1")
  ref <- log(mu / (1 - mu))
  expect_equal(res, ref, tolerance = 1e-8)
})

test_that("parametric_mean_link is monotone increasing in mu", {
  mu  <- seq(0.1, 0.9, by = 0.1)
  for (type in c("plogit1", "plogit2")) {
    res <- parametric_mean_link(mu, lambda = 1.5, type = type)
    expect_true(all(diff(res) > 0), label = paste(type, "is monotone increasing"))
  }
})

test_that("parametric_mean_link default type argument is matched correctly", {
  expect_no_error(parametric_mean_link(0.5, lambda = 1))
})

# --- parametric_mean_link_inv -------------------------------------------------

test_that("parametric_mean_link_inv returns values in (0, 1)", {
  eta <- c(-2, -1, 0, 1, 2)
  for (type in c("plogit1", "plogit2")) {
    res <- parametric_mean_link_inv(eta, lambda = 1.5, type = type)
    expect_true(all(res > 0 & res < 1), label = paste(type, "inverse in (0,1)"))
  }
})

test_that("parametric_mean_link_inv is the inverse of parametric_mean_link", {
  mu  <- c(0.2, 0.4, 0.6, 0.8)
  for (type in c("plogit1", "plogit2")) {
    eta <- parametric_mean_link(mu, lambda = 1.3, type = type)
    mu_back <- parametric_mean_link_inv(eta, lambda = 1.3, type = type)
    expect_equal(mu_back, mu, tolerance = 1e-8,
                 label = paste(type, "round-trip"))
  }
})

test_that("parametric_mean_link_inv plogit2 equals expit(eta) when lambda = 1", {
  eta <- c(-1, 0, 1)
  res <- parametric_mean_link_inv(eta, lambda = 1, type = "plogit2")
  ref <- 1 / (1 + exp(-eta))
  expect_equal(res, ref, tolerance = 1e-8)
})

# --- parametric_mean_link_deriv1 ----------------------------------------------

test_that("parametric_mean_link_deriv1 rejects mu outside (0, 1)", {
  expect_error(parametric_mean_link_deriv1(0,   lambda = 1, type = "plogit2"), "interval \\(0, 1\\)")
  expect_error(parametric_mean_link_deriv1(1.2, lambda = 1, type = "plogit1"), "interval \\(0, 1\\)")
})

test_that("parametric_mean_link_deriv1 returns positive finite values", {
  mu  <- c(0.2, 0.5, 0.8)
  for (type in c("plogit1", "plogit2")) {
    res <- parametric_mean_link_deriv1(mu, lambda = 1.5, type = type)
    expect_true(all(is.finite(res)), label = paste(type, "finite"))
    expect_true(all(res > 0),        label = paste(type, "positive"))
  }
})

test_that("parametric_mean_link_deriv1 plogit2 equals 1/(mu*(1-mu)) when lambda = 1", {
  mu  <- c(0.2, 0.5, 0.8)
  res <- parametric_mean_link_deriv1(mu, lambda = 1, type = "plogit2")
  ref <- 1 / (mu * (1 - mu))
  expect_equal(res, ref, tolerance = 1e-8)
})

test_that("parametric_mean_link_deriv1 approximates numerical derivative", {
  mu  <- 0.4
  lam <- 1.5
  h   <- 1e-6
  for (type in c("plogit1", "plogit2")) {
    num_deriv <- (parametric_mean_link(mu + h, lam, type) -
                    parametric_mean_link(mu - h, lam, type)) / (2 * h)
    ana_deriv <- parametric_mean_link_deriv1(mu, lam, type)
    expect_equal(ana_deriv, num_deriv, tolerance = 1e-5,
                 label = paste(type, "numerical vs analytical deriv1"))
  }
})

# --- parametric_mean_link_inv_deriv1 ------------------------------------------

test_that("parametric_mean_link_inv_deriv1 returns positive finite values", {
  eta <- c(-2, 0, 2)
  for (type in c("plogit1", "plogit2")) {
    res <- parametric_mean_link_inv_deriv1(eta, lambda = 1.5, type = type)
    expect_true(all(is.finite(res)), label = paste(type, "finite"))
    expect_true(all(res > 0),        label = paste(type, "positive"))
  }
})

test_that("parametric_mean_link_inv_deriv1 approximates numerical derivative", {
  eta <- 0.3
  lam <- 1.5
  h   <- 1e-6
  for (type in c("plogit1", "plogit2")) {
    num_deriv <- (parametric_mean_link_inv(eta + h, lam, type) -
                    parametric_mean_link_inv(eta - h, lam, type)) / (2 * h)
    ana_deriv <- parametric_mean_link_inv_deriv1(eta, lam, type)
    expect_equal(ana_deriv, num_deriv, tolerance = 1e-5,
                 label = paste(type, "numerical vs analytical inv_deriv1"))
  }
})

# --- parametric_mean_link_deriv2 ----------------------------------------------

test_that("parametric_mean_link_deriv2 rejects mu outside (0, 1)", {
  expect_error(parametric_mean_link_deriv2(0,   lambda = 1, type = "plogit2"), "interval \\(0, 1\\)")
  expect_error(parametric_mean_link_deriv2(1.1, lambda = 1, type = "plogit1"), "interval \\(0, 1\\)")
})

test_that("parametric_mean_link_deriv2 returns finite values", {
  mu  <- c(0.2, 0.5, 0.8)
  for (type in c("plogit1", "plogit2")) {
    res <- parametric_mean_link_deriv2(mu, lambda = 1.5, type = type)
    expect_true(all(is.finite(res)), label = paste(type, "finite"))
  }
})

test_that("parametric_mean_link_deriv2 approximates numerical second derivative", {
  mu  <- 0.4
  lam <- 1.5
  h   <- 1e-5
  for (type in c("plogit1", "plogit2")) {
    num_deriv2 <- (parametric_mean_link(mu + h, lam, type) -
                     2 * parametric_mean_link(mu, lam, type) +
                     parametric_mean_link(mu - h, lam, type)) / h^2
    ana_deriv2 <- parametric_mean_link_deriv2(mu, lam, type)
    expect_equal(ana_deriv2, num_deriv2, tolerance = 1e-4,
                 label = paste(type, "numerical vs analytical deriv2"))
  }
})

# ==============================================================================
# 2. FIXED MEAN LINK FUNCTIONS
# ==============================================================================

fixed_types <- c("logit", "probit", "loglog", "cloglog", "cauchit")

# --- fixed_mean_link ----------------------------------------------------------

test_that("fixed_mean_link rejects mu outside (0, 1)", {
  for (type in fixed_types) {
    expect_error(fixed_mean_link(0,    type = type), "interval \\(0, 1\\)",
                 label = paste(type, "mu = 0"))
    expect_error(fixed_mean_link(1,    type = type), "interval \\(0, 1\\)",
                 label = paste(type, "mu = 1"))
    expect_error(fixed_mean_link(-0.1, type = type), "interval \\(0, 1\\)",
                 label = paste(type, "mu < 0"))
  }
})

test_that("fixed_mean_link returns finite numeric values for all types", {
  mu <- c(0.2, 0.5, 0.8)
  for (type in fixed_types) {
    res <- fixed_mean_link(mu, type = type)
    expect_true(is.numeric(res),       label = paste(type, "numeric"))
    expect_length(res, 3)
    expect_true(all(is.finite(res)),   label = paste(type, "finite"))
  }
})

test_that("fixed_mean_link logit equals log(mu / (1 - mu))", {
  mu  <- c(0.2, 0.5, 0.8)
  res <- fixed_mean_link(mu, type = "logit")
  expect_equal(res, log(mu / (1 - mu)), tolerance = 1e-12)
})

test_that("fixed_mean_link is monotone increasing in mu for all types", {
  mu <- seq(0.1, 0.9, by = 0.1)
  for (type in fixed_types) {
    res <- fixed_mean_link(mu, type = type)
    expect_true(all(diff(res) > 0),
                label = paste(type, "monotone increasing"))
  }
})

# --- fixed_mean_link_inv ------------------------------------------------------

test_that("fixed_mean_link_inv returns values in (0, 1) for all types", {
  eta <- c(-3, -1, 0, 1, 3)
  for (type in fixed_types) {
    res <- fixed_mean_link_inv(eta, type = type)
    expect_true(all(res > 0 & res < 1), label = paste(type, "in (0,1)"))
  }
})

test_that("fixed_mean_link_inv is the inverse of fixed_mean_link", {
  mu <- c(0.2, 0.4, 0.6, 0.8)
  for (type in fixed_types) {
    eta    <- fixed_mean_link(mu, type = type)
    mu_back <- fixed_mean_link_inv(eta, type = type)
    expect_equal(mu_back, mu, tolerance = 1e-10,
                 label = paste(type, "round-trip"))
  }
})

# --- fixed_mean_link_deriv1 ---------------------------------------------------

test_that("fixed_mean_link_deriv1 rejects mu outside (0, 1)", {
  for (type in fixed_types) {
    expect_error(fixed_mean_link_deriv1(0,   type = type), "interval \\(0, 1\\)",
                 label = paste(type, "mu = 0"))
    expect_error(fixed_mean_link_deriv1(1.5, type = type), "interval \\(0, 1\\)",
                 label = paste(type, "mu > 1"))
  }
})

test_that("fixed_mean_link_deriv1 returns positive finite values for all types", {
  mu <- c(0.2, 0.5, 0.8)
  for (type in fixed_types) {
    res <- fixed_mean_link_deriv1(mu, type = type)
    expect_true(all(is.finite(res)), label = paste(type, "finite"))
    expect_true(all(res > 0),        label = paste(type, "positive"))
  }
})

test_that("fixed_mean_link_deriv1 logit equals 1 / (mu * (1 - mu))", {
  mu  <- c(0.2, 0.5, 0.8)
  res <- fixed_mean_link_deriv1(mu, type = "logit")
  expect_equal(res, 1 / (mu * (1 - mu)), tolerance = 1e-12)
})

test_that("fixed_mean_link_deriv1 approximates numerical derivative for all types", {
  mu <- 0.4
  h  <- 1e-6
  for (type in fixed_types) {
    num_deriv <- (fixed_mean_link(mu + h, type) -
                    fixed_mean_link(mu - h, type)) / (2 * h)
    ana_deriv <- fixed_mean_link_deriv1(mu, type)
    expect_equal(ana_deriv, num_deriv, tolerance = 1e-5,
                 label = paste(type, "numerical vs analytical"))
  }
})

# --- fixed_mean_link_deriv2 ---------------------------------------------------

test_that("fixed_mean_link_deriv2 rejects mu outside (0, 1)", {
  for (type in fixed_types) {
    expect_error(fixed_mean_link_deriv2(0,   type = type), "interval \\(0, 1\\)",
                 label = paste(type, "mu = 0"))
    expect_error(fixed_mean_link_deriv2(1.2, type = type), "interval \\(0, 1\\)",
                 label = paste(type, "mu > 1"))
  }
})

test_that("fixed_mean_link_deriv2 returns finite values for all types", {
  mu <- c(0.2, 0.5, 0.8)
  for (type in fixed_types) {
    res <- fixed_mean_link_deriv2(mu, type = type)
    expect_true(all(is.finite(res)), label = paste(type, "finite"))
  }
})

test_that("fixed_mean_link_deriv2 approximates numerical second derivative", {
  mu <- 0.4
  h  <- 1e-5
  for (type in fixed_types) {
    num_d2 <- (fixed_mean_link(mu + h, type) -
                 2 * fixed_mean_link(mu, type) +
                 fixed_mean_link(mu - h, type)) / h^2
    ana_d2 <- fixed_mean_link_deriv2(mu, type)
    expect_equal(ana_d2, num_d2, tolerance = 1e-4,
                 label = paste(type, "numerical vs analytical"))
  }
})

# --- fixed_mean_link_inv_deriv1 -----------------------------------------------

test_that("fixed_mean_link_inv_deriv1 returns positive finite values for all types", {
  eta <- c(-2, -0.5, 0, 0.5, 2)
  for (type in fixed_types) {
    res <- fixed_mean_link_inv_deriv1(eta, type = type)
    expect_true(all(is.finite(res)), label = paste(type, "finite"))
    expect_true(all(res > 0),        label = paste(type, "positive"))
  }
})

test_that("fixed_mean_link_inv_deriv1 approximates numerical derivative of inverse", {
  eta <- 0.3
  h   <- 1e-6
  for (type in fixed_types) {
    num_deriv <- (fixed_mean_link_inv(eta + h, type) -
                    fixed_mean_link_inv(eta - h, type)) / (2 * h)
    ana_deriv <- fixed_mean_link_inv_deriv1(eta, type)
    expect_equal(ana_deriv, num_deriv, tolerance = 1e-5,
                 label = paste(type, "numerical vs analytical"))
  }
})

# ==============================================================================
# 3. DISPERSION LINK FUNCTIONS
# ==============================================================================

disp_types <- c("log", "sqrt", "identity")

# --- dispersion_link ----------------------------------------------------------

test_that("dispersion_link rejects non-positive sigma2", {
  for (type in disp_types) {
    expect_error(dispersion_link(0,    type = type), "positive",
                 label = paste(type, "sigma2 = 0"))
    expect_error(dispersion_link(-0.5, type = type), "positive",
                 label = paste(type, "sigma2 < 0"))
  }
})

test_that("dispersion_link returns finite numeric values for all types", {
  sigma2 <- c(0.5, 1.0, 2.0)
  for (type in disp_types) {
    res <- dispersion_link(sigma2, type = type)
    expect_true(is.numeric(res),       label = paste(type, "numeric"))
    expect_length(res, 3)
    expect_true(all(is.finite(res)),   label = paste(type, "finite"))
  }
})

test_that("dispersion_link log equals log(sigma2)", {
  sigma2 <- c(0.5, 1.0, 2.0)
  expect_equal(dispersion_link(sigma2, type = "log"), log(sigma2), tolerance = 1e-12)
})

test_that("dispersion_link sqrt equals sqrt(sigma2)", {
  sigma2 <- c(0.5, 1.0, 4.0)
  expect_equal(dispersion_link(sigma2, type = "sqrt"), sqrt(sigma2), tolerance = 1e-12)
})

test_that("dispersion_link identity returns sigma2 unchanged", {
  sigma2 <- c(0.3, 1.0, 2.5)
  expect_equal(dispersion_link(sigma2, type = "identity"), sigma2)
})

# --- dispersion_link_inv ------------------------------------------------------

test_that("dispersion_link_inv returns appropriate values based on link type", {
  for (type in disp_types) {
    if (type == "log") {
      eta <- c(-2, -1, 0, 1, 2)
      res <- dispersion_link_inv(eta, type = type)
      expect_true(all(res > 0), label = paste(type, "returns positive values"))

    } else if (type == "sqrt") {
      eta <- c(-2, -1, 1, 2)
      res <- dispersion_link_inv(eta, type = type)
      expect_true(all(res > 0), label = paste(type, "returns positive values when eta != 0"))

      expect_equal(dispersion_link_inv(0, type = "sqrt"), 0)

    } else if (type == "identity") {
      eta <- c(0.1, 0.5, 1, 2)
      res <- dispersion_link_inv(eta, type = type)
      expect_true(all(res > 0), label = paste(type, "returns positive values when eta > 0"))

    }
  }
})

test_that("dispersion_link_inv is the inverse of dispersion_link for valid inputs", {
  sigma2 <- c(0.3, 0.8, 1.5, 3.0)
  for (type in disp_types) {
    eta <- dispersion_link(sigma2, type = type)
    sigma2_back <- dispersion_link_inv(eta, type = type)
    expect_equal(sigma2_back, sigma2, tolerance = 1e-10,
                 label = paste(type, "round-trip"))

    if (type %in% c("identity", "sqrt")) {
      expect_true(all(sigma2 > 0), label = paste(type, "sigma2 must be positive"))
    }
  }
})

test_that("dispersion_link_inv requires appropriate eta for identity and sqrt links", {
  expect_true(all(dispersion_link_inv(c(0.1, 0.5, 2), type = "identity") > 0))

  expect_true(all(dispersion_link_inv(c(-2, -1, 1, 2), type = "sqrt") > 0))
  expect_equal(dispersion_link_inv(0, type = "sqrt"), 0)
})

# --- dispersion_link_deriv1 ---------------------------------------------------

test_that("dispersion_link_deriv1 rejects non-positive sigma2", {
  for (type in disp_types) {
    expect_error(dispersion_link_deriv1(0,    type = type), "positive",
                 label = paste(type, "sigma2 = 0"))
    expect_error(dispersion_link_deriv1(-0.1, type = type), "positive",
                 label = paste(type, "sigma2 < 0"))
  }
})

test_that("dispersion_link_deriv1 returns finite values for all types", {
  sigma2 <- c(0.5, 1.0, 2.0)
  for (type in disp_types) {
    res <- dispersion_link_deriv1(sigma2, type = type)
    expect_true(all(is.finite(res)), label = paste(type, "finite"))
  }
})

test_that("dispersion_link_deriv1 log equals 1 / sigma2", {
  sigma2 <- c(0.5, 1.0, 2.0)
  expect_equal(dispersion_link_deriv1(sigma2, type = "log"), 1 / sigma2,
               tolerance = 1e-12)
})

test_that("dispersion_link_deriv1 identity returns a vector of ones", {
  sigma2 <- c(0.3, 1.0, 2.5)
  res    <- dispersion_link_deriv1(sigma2, type = "identity")
  expect_equal(res, rep(1, 3))
})

test_that("dispersion_link_deriv1 approximates numerical derivative for all types", {
  sigma2 <- 1.0
  h      <- 1e-6
  for (type in disp_types) {
    num_deriv <- (dispersion_link(sigma2 + h, type) -
                    dispersion_link(sigma2 - h, type)) / (2 * h)
    ana_deriv <- dispersion_link_deriv1(sigma2, type)
    expect_equal(ana_deriv, num_deriv, tolerance = 1e-5,
                 label = paste(type, "numerical vs analytical"))
  }
})

# --- dispersion_link_inv_deriv1 -----------------------------------------------

test_that("dispersion_link_inv_deriv1 returns appropriate values based on link type", {
  for (type in disp_types) {
    if (type == "log") {
      eta <- c(-2, -1, 0, 1, 2)
      res <- dispersion_link_inv_deriv1(eta, type = type)
      expect_true(all(is.finite(res)), label = paste(type, "finite"))
      expect_true(all(res > 0),        label = paste(type, "positive"))

    } else if (type == "sqrt") {
      eta <- c(-2, -1, 1, 2)
      res <- dispersion_link_inv_deriv1(eta, type = type)
      expect_true(all(is.finite(res)), label = paste(type, "finite"))

      expect_true(all(res != 0), label = paste(type, "non-zero derivative"))

      expect_equal(dispersion_link_inv_deriv1(0, type = "sqrt"), 0)

    } else if (type == "identity") {
      eta <- c(-2, -1, 0, 1, 2)
      res <- dispersion_link_inv_deriv1(eta, type = type)
      expect_true(all(is.finite(res)), label = paste(type, "finite"))
      expect_true(all(res == 1),       label = paste(type, "equals 1"))

    }
  }
})

test_that("dispersion_link_inv_deriv1 log equals exp(eta)", {
  eta <- c(-1, 0, 1, 2)
  expect_equal(dispersion_link_inv_deriv1(eta, type = "log"), exp(eta),
               tolerance = 1e-12)
})

test_that("dispersion_link_inv_deriv1 identity returns a vector of ones", {
  eta <- c(-1, 0, 1, 2)
  res <- dispersion_link_inv_deriv1(eta, type = "identity")
  expect_equal(res, rep(1, 4))
})

test_that("dispersion_link_inv_deriv1 approximates numerical derivative of inverse", {
  test_cases <- list(
    log = 0.5,
    sqrt = 2.0,
    identity = 0.5
  )

  h <- 1e-6
  for (type in disp_types) {
    eta <- test_cases[[type]]
    num_deriv <- (dispersion_link_inv(eta + h, type) -
                    dispersion_link_inv(eta - h, type)) / (2 * h)
    ana_deriv <- dispersion_link_inv_deriv1(eta, type)
    expect_equal(ana_deriv, num_deriv, tolerance = 1e-5,
                 label = paste(type, "numerical vs analytical"))
  }
})

library(testthat)
library(SimplexRegression)

# ==============================================================================
# test-simplexreg_fit.R
# Tests for: simplexreg, simplexreg.fit, simplexreg.control
# ==============================================================================

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
dat <- data.frame(y = y, x1 = x1, x2 = x2, z1 = z1)

fit_plogit1 <- simplexreg(y ~ x1 + x2 | z1, data = dat, link.mu = "plogit1")
fit_plogit2 <- simplexreg(y ~ x1 + x2 | z1, data = dat, link.mu = "plogit2")
fit_logit   <- simplexreg(y ~ x1 + x2 | z1, data = dat, link.mu = "logit")
fit_probit  <- simplexreg(y ~ x1 + x2 | z1, data = dat, link.mu = "probit")
fit_fixdisp <- simplexreg(y ~ x1 + x2,  data = dat, link.mu = "plogit1")
fit_null    <- simplexreg(y ~ 1,  data = dat, link.mu = "logit")

# ==============================================================================
# 1. simplexreg.control
# ==============================================================================

test_that("simplexreg.control returns a list with expected fields", {
  ctrl <- simplexreg.control()
  expect_type(ctrl, "list")
  expect_true(all(c("method", "maxit", "gradient", "hessian", "start",
                    "fsmaxit", "fstol", "reltol", "fnscale") %in% names(ctrl)))
})

test_that("simplexreg.control defaults are correct", {
  ctrl <- simplexreg.control()
  expect_equal(ctrl$method,   "BFGS")
  expect_equal(ctrl$maxit,    5000L)
  expect_true(ctrl$gradient)
  expect_false(ctrl$hessian)
  expect_equal(ctrl$fsmaxit,  500L)
  expect_equal(ctrl$fstol,    1e-8)
  expect_equal(ctrl$fnscale, -1)   # must be -1 (maximization)
})

test_that("simplexreg.control rejects invalid method", {
  expect_error(simplexreg.control(method = "Newton"))
})

test_that("simplexreg.control accepts valid methods", {
  for (m in c("BFGS", "Nelder-Mead", "CG", "L-BFGS-B", "SANN")) {
    expect_no_error(simplexreg.control(method = m))
  }
})

test_that("simplexreg.control warns if fnscale is modified", {
  expect_warning(simplexreg.control(fnscale = 1), "fnscale must not be modified")
})

test_that("simplexreg.control L-BFGS-B drops reltol and sets factr/pgtol", {
  ctrl <- simplexreg.control(method = "L-BFGS-B")
  expect_null(ctrl$reltol)
  expect_false(is.null(ctrl$factr))
  expect_false(is.null(ctrl$pgtol))
})

test_that("simplexreg.control accepts user-supplied start values", {
  ctrl <- simplexreg.control(start = c(0, 1, 2))
  expect_equal(ctrl$start, c(0, 1, 2))
})

test_that("simplexreg.control sets default reltol when reltol = NULL and method != L-BFGS-B", {
  ctrl <- simplexreg.control(method = "BFGS", reltol = NULL)
  expect_equal(ctrl$reltol, .Machine$double.eps^(0.5))
})

# ==============================================================================
# 2. simplexreg — object class and structure
# ==============================================================================

test_that("simplexreg returns an object of class simplexregression", {
  expect_s3_class(fit_plogit1, "simplexregression")
  expect_s3_class(fit_logit,   "simplexregression")
})

test_that("simplexreg object contains all required components", {
  required <- c("coefficients", "fitted.values", "optim", "scoring",
                "mu.fv", "mu.lp", "mu.x", "mu.link", "mu.df",
                "sigma2.fv", "sigma2.lp", "sigma2.x", "sigma2.link", "sigma2.df",
                "lambda.fv", "df.residual", "nobs", "loglik", "vcov",
                "residuals", "AIC", "BIC", "HQIC", "R2_FC", "R2_N",
                "zstat", "pvalues", "y", "converged")
  for (comp in required) {
    expect_true(comp %in% names(fit_plogit1),
                label = paste("component", comp, "present"))
  }
})

test_that("simplexreg coefficients list has correct sub-elements", {
  cf <- fit_plogit1$coefficients
  expect_true(!is.null(cf$mean))
  expect_true(!is.null(cf$dispersion))
  expect_true(!is.null(cf$lambda))    # parametric: lambda present
})

test_that("simplexreg coefficients lambda is NULL or NA for fixed link", {
  lam <- fit_logit$coefficients$lambda
  expect_true(is.null(lam) || (length(lam) == 1L && is.na(lam)))
})

test_that("simplexreg coefficients lambda is numeric for parametric link", {
  expect_true(is.numeric(fit_plogit1$coefficients$lambda))
  expect_true(is.finite(fit_plogit1$coefficients$lambda))
  expect_true(fit_plogit1$coefficients$lambda > 0)
})

# ==============================================================================
# 3. simplexreg — fitted values and linear predictors
# ==============================================================================

test_that("simplexreg mu.fv is in (0, 1) for all link types", {
  for (fit in list(fit_plogit1, fit_plogit2, fit_logit, fit_probit)) {
    expect_true(all(fit$mu.fv > 0 & fit$mu.fv < 1),
                label = paste(fit$mu.link, "mu.fv in (0,1)"))
  }
})

test_that("simplexreg sigma2.fv is positive for all link types", {
  for (fit in list(fit_plogit1, fit_logit)) {
    expect_true(all(fit$sigma2.fv > 0),
                label = paste(fit$mu.link, "sigma2.fv > 0"))
  }
})

test_that("simplexreg mu.fv and fitted.values are identical", {
  expect_equal(fit_plogit1$mu.fv, fit_plogit1$fitted.values)
  expect_equal(fit_logit$mu.fv,   fit_logit$fitted.values)
})

test_that("simplexreg mu.fv has length n", {
  expect_length(fit_plogit1$mu.fv, n)
  expect_length(fit_logit$mu.fv,   n)
})

test_that("simplexreg sigma2.fv has length n", {
  expect_length(fit_plogit1$sigma2.fv, n)
})

test_that("simplexreg linear predictor mu.lp has length n", {
  expect_length(fit_plogit1$mu.lp, n)
  expect_true(all(is.finite(fit_plogit1$mu.lp)))
})

# ==============================================================================
# 4. simplexreg — design matrices
# ==============================================================================

test_that("simplexreg mu.x has correct dimensions", {
  expect_equal(nrow(fit_plogit1$mu.x), n)
  expect_equal(ncol(fit_plogit1$mu.x), 3L)  # intercept + x1 + x2
})

test_that("simplexreg sigma2.x has correct dimensions (constant dispersion)", {
  expect_equal(nrow(fit_fixdisp$sigma2.x), n)
  expect_equal(ncol(fit_fixdisp$sigma2.x), 1L)  # intercept only
})

test_that("simplexreg sigma2.x has correct dimensions (variable dispersion)", {
  expect_equal(ncol(fit_plogit1$sigma2.x), 2L)  # intercept + z1
})

test_that("simplexreg mu.x first column is all ones (intercept)", {
  expect_true(all(fit_plogit1$mu.x[, 1] == 1))
})

# ==============================================================================
# 5. simplexreg — degrees of freedom and model dimensions
# ==============================================================================

test_that("simplexreg nobs equals n", {
  expect_equal(fit_plogit2$nobs, n)
  expect_equal(fit_logit$nobs,   n)
})

test_that("simplexreg df.residual is correct for parametric link", {
  r <- ncol(fit_plogit2$mu.x) + ncol(fit_plogit2$sigma2.x) + 1L  # +1 for lambda
  expect_equal(fit_plogit2$df.residual, n - r)
})

test_that("simplexreg df.residual is correct for fixed link", {
  r <- ncol(fit_logit$mu.x) + ncol(fit_logit$sigma2.x)
  expect_equal(fit_logit$df.residual, n - r)
})

test_that("simplexreg mu.df equals number of mean parameters", {
  expect_equal(fit_plogit2$mu.df, ncol(fit_plogit2$mu.x))
})

test_that("simplexreg sigma2.df equals number of dispersion parameters", {
  expect_equal(fit_plogit2$sigma2.df, ncol(fit_plogit2$sigma2.x))
})

# ==============================================================================
# 6. simplexreg — log-likelihood and information criteria
# ==============================================================================

test_that("simplexreg loglik is finite", {
  expect_true(is.finite(fit_plogit1$loglik))
  expect_true(is.finite(fit_logit$loglik))
})

test_that("simplexreg AIC, BIC, HQIC are finite", {
  for (ic in c("AIC", "BIC", "HQIC")) {
    expect_true(is.finite(fit_plogit1[[ic]]),
                label = paste(ic, "is finite"))
  }
})

# ==============================================================================
# 7. simplexreg — variance-covariance matrix
# ==============================================================================

test_that("simplexreg vcov is a square numeric matrix", {
  vc <- fit_plogit1$vcov
  expect_true(is.matrix(vc))
  expect_equal(nrow(vc), ncol(vc))
})

test_that("simplexreg vcov dimension matches number of parameters", {
  r_par <- ncol(fit_plogit1$mu.x) + ncol(fit_plogit1$sigma2.x) + 1L  # +1 lambda
  expect_equal(dim(fit_plogit1$vcov), c(r_par, r_par))
})

test_that("simplexreg vcov is positive definite (all eigenvalues > 0)", {
  ev <- eigen(fit_plogit1$vcov, only.values = TRUE)$values
  expect_true(all(ev > 0))
})

test_that("simplexreg vcov is symmetric", {
  vc <- fit_plogit1$vcov
  expect_equal(vc, t(vc), tolerance = 1e-12)
})

test_that("simplexreg vcov has named rows and columns", {
  vc <- fit_plogit1$vcov
  expect_false(is.null(rownames(vc)))
  expect_false(is.null(colnames(vc)))
})

# ==============================================================================
# 8. simplexreg — residuals, R-squared, zstat
# ==============================================================================

test_that("simplexreg residuals have length n", {
  expect_length(fit_plogit1$residuals, n)
  expect_length(fit_logit$residuals,   n)
})

test_that("simplexreg quantile residuals are approximately standard normal", {
  res <- fit_plogit1$residuals
  expect_true(abs(mean(res)) < 0.5)
  expect_true(abs(sd(res) - 1) < 0.5)
})

test_that("simplexreg R2_N is in (0, 1)", {
  expect_true(fit_plogit1$R2_N > 0 & fit_plogit1$R2_N < 1)
  expect_true(fit_logit$R2_N   > 0 & fit_logit$R2_N   < 1)
})

test_that("simplexreg R2_FC is in (0, 1)", {
  expect_true(fit_plogit1$R2_FC > 0 & fit_plogit1$R2_FC < 1)
})

test_that("simplexreg zstat and pvalues have correct length", {
  r_par <- ncol(fit_plogit2$mu.x) + ncol(fit_plogit2$sigma2.x) + 1L
  expect_length(fit_plogit2$zstat,   r_par)
  expect_length(fit_plogit2$pvalues, r_par)
})

test_that("simplexreg zstat values are non-negative", {
  expect_true(all(fit_plogit2$zstat >= 0))
})

test_that("simplexreg pvalues are in [0, 1]", {
  expect_true(all(fit_plogit2$pvalues >= 0 & fit_plogit2$pvalues <= 1))
})

# ==============================================================================
# 9. simplexreg — formula and call
# ==============================================================================

test_that("simplexreg stores the call", {
  expect_true(!is.null(fit_plogit2$call))
  expect_true(inherits(fit_plogit2$call, "call"))
})

test_that("simplexreg stores the formula", {
  expect_true(!is.null(fit_plogit2$formula))
  expect_true(inherits(fit_plogit2$formula, "formula"))
})

test_that("simplexreg mu.link is stored correctly", {
  expect_equal(fit_plogit2$mu.link, "plogit2")
  expect_equal(fit_logit$mu.link,   "logit")
})

test_that("simplexreg sigma2.link is stored", {
  expect_true(!is.null(fit_plogit2$sigma2.link))
  expect_true(fit_plogit2$sigma2.link %in%
                c("log", "sqrt", "identity"))
})

# ==============================================================================
# 10. simplexreg — convergence and optim info
# ==============================================================================

test_that("simplexreg converged flag is logical", {
  expect_true(is.logical(fit_plogit1$converged))
  expect_true(is.logical(fit_logit$converged))
})

test_that("simplexreg optim list contains required fields", {
  opt <- fit_plogit1$optim
  expect_true(all(c("start", "convergence", "counts", "method") %in% names(opt)))
})

test_that("simplexreg optim method is a valid string", {
  expect_true(fit_plogit1$optim$method %in%
                c("BFGS", "Nelder-Mead", "CG", "L-BFGS-B", "SANN"))
})

test_that("simplexreg scoring field is a non-negative integer", {
  expect_true(is.numeric(fit_plogit1$scoring))
  expect_true(fit_plogit1$scoring >= 0)
})

# ==============================================================================
# 11. simplexreg — input validation
# ==============================================================================

test_that("simplexreg rejects response outside (0, 1)", {
  dat_bad      <- dat
  dat_bad$y[1] <- 0
  expect_error(simplexreg(y ~ x1 + x2 | z1, data = dat_bad, link.mu = "logit"),
               "open interval \\(0, 1\\)")

  dat_bad$y[1] <- 1
  expect_error(simplexreg(y ~ x1 + x2 | z1, data = dat_bad, link.mu = "logit"),
               "open interval \\(0, 1\\)")
})

test_that("simplexreg requires a formula object", {
  expect_error(simplexreg("y ~ x1", data = dat), "formula")
})

# ==============================================================================
# 12. simplexreg — all mean link functions fit without error
# ==============================================================================

all_links <- c("logit", "probit", "loglog", "cloglog", "cauchit",
               "plogit1", "plogit2")

test_that("simplexreg fits without error for all mean link functions", {
  for (lnk in all_links) {
    expect_no_error(
      simplexreg(y ~ x1 + x2 | z1, data = dat, link.mu = lnk)
    )
  }
})

test_that("simplexreg produces finite log-likelihood for all mean link functions", {
  for (lnk in all_links) {
    fit <- simplexreg(y ~ x1 + x2 | z1, data = dat, link.mu = lnk)
    expect_true(is.finite(fit$loglik),
                label = paste("link.mu =", lnk, "loglik finite"))
  }
})

# ==============================================================================
# 13. simplexreg — all dispersion link functions fit without error
# ==============================================================================

disp_links <- c("log", "sqrt", "identity")

test_that("simplexreg fits without error for all dispersion link functions", {
  for (lnk in disp_links) {
    expect_no_error(
      simplexreg(y ~ x1 + x2 | z1, data = dat,
                 link.mu = "logit", link.sigma2 = lnk)
    )
  }
})

test_that("simplexreg produces finite log-likelihood for all dispersion links", {
  for (lnk in disp_links) {
    fit <- simplexreg(y ~ x1 + x2 | z1, data = dat,
                      link.mu = "logit", link.sigma2 = lnk)
    expect_true(is.finite(fit$loglik),
                label = paste("link.sigma2 =", lnk, "loglik finite"))
  }
})

# ==============================================================================
# 14. simplexreg — null (intercept-only) model
# ==============================================================================

test_that("simplexreg fits null model without error", {
  expect_no_error(simplexreg(y ~ 1, data = dat, link.mu = "logit"))
})

test_that("simplexreg null model has 1 mean coefficient", {
  expect_equal(length(fit_null$coefficients$mean), 1L)
})

test_that("simplexreg null model mu.x has 1 column (intercept only)", {
  expect_equal(ncol(fit_null$mu.x), 1L)
})

test_that("simplexreg null model loglik is less than full model loglik", {
  # Adding predictors should improve the log-likelihood
  expect_true(fit_null$loglik <= fit_logit$loglik)
})

# ==============================================================================
# 15. simplexreg — variable dispersion model
# ==============================================================================

test_that("simplexreg variable-dispersion model has correct sigma2.df", {
  expect_equal(fit_plogit1$sigma2.df, 2L)  # intercept + z1
})

test_that("simplexreg variable-dispersion sigma2.fv values are not all equal", {
  expect_true(length(unique(round(fit_plogit1$sigma2.fv, 6))) > 1L)
})

test_that("simplexreg variable-dispersion loglik >= constant-dispersion loglik", {
  # More flexible model should have at least as good a log-likelihood
  expect_true(fit_plogit1$loglik >= fit_fixdisp$loglik - 0.5)
})

# ==============================================================================
# 16. simplexreg — formula without | separator (simple formula)
# ==============================================================================

test_that("simplexreg accepts simple formula without | separator", {
  expect_no_error(simplexreg(y ~ x1 + x2, data = dat, link.mu = "logit"))
})

test_that("simplexreg simple formula produces constant dispersion model", {
  fit_simple <- simplexreg(y ~ x1 + x2, data = dat, link.mu = "logit")
  expect_equal(ncol(fit_simple$sigma2.x), 1L)
  expect_equal(length(unique(round(fit_simple$sigma2.fv, 8))), 1L)
})

# ==============================================================================
# 17. simplexreg — starting values via control
# ==============================================================================

test_that("simplexreg accepts user-supplied starting values via control", {
  # Use estimated parameters from a converged fit as starting values
  beta0  <- fit_logit$coefficients$mean
  delta0 <- fit_logit$coefficients$dispersion
  start  <- c(beta0, delta0)

  expect_no_error(
    simplexreg(y ~ x1 + x2 | z1, data = dat, link.mu = "logit",
               control = simplexreg.control(start = start))
  )
})

test_that("simplexreg with user start produces consistent log-likelihood", {
  beta0  <- fit_logit$coefficients$mean
  delta0 <- fit_logit$coefficients$dispersion
  start  <- c(beta0, delta0)

  fit_warm <- simplexreg(y ~ x1 + x2 | z1, data = dat, link.mu = "logit",
                         control = simplexreg.control(start = start))
  expect_equal(fit_warm$loglik, fit_logit$loglik, tolerance = 1e-4)
})

# ==============================================================================
# 18. simplexreg — Fisher scoring can be disabled
# ==============================================================================

test_that("simplexreg with fsmaxit = 0 runs without error", {
  expect_no_error(
    simplexreg(y ~ x1 + x2 | z1, data = dat, link.mu = "logit",
               control = simplexreg.control(fsmaxit = 0L))
  )
})

test_that("simplexreg with fsmaxit = 0 has scoring = 0", {
  fit_no_fs <- simplexreg(y ~ x1 + x2 | z1, data = dat, link.mu = "logit",
                          control = simplexreg.control(fsmaxit = 0L))
  expect_equal(fit_no_fs$scoring, 0L)
})

# ==============================================================================
# 19. simplexreg — parametric vs fixed link: lambda presence
# ==============================================================================

test_that("simplexreg lambda.fv is numeric and positive for plogit1", {
  expect_true(is.numeric(fit_plogit1$lambda.fv))
  expect_true(fit_plogit1$lambda.fv > 0)
})

test_that("simplexreg lambda.fv is NA or NULL for fixed link", {
  expect_true(is.null(fit_logit$lambda.fv) || is.na(fit_logit$lambda.fv))
  expect_true(is.null(fit_probit$lambda.fv) || is.na(fit_probit$lambda.fv))
})

test_that("simplexreg works when data is not supplied (uses environment)", {
  y_env  <- dat$y
  x1_env <- dat$x1
  x2_env <- dat$x2
  z1_env <- dat$z1
  expect_no_error(simplexreg(y_env ~ x1_env + x2_env | z1_env, link.mu = "logit"))
})

test_that("simplexreg rejects response outside (0,1) for simple formula", {
  dat_bad      <- dat
  dat_bad$y[1] <- 0
  expect_error(
    simplexreg(y ~ x1 + x2, data = dat_bad, link.mu = "logit"),
    "open interval \\(0, 1\\)"
  )
})

test_that("simplexreg accepts offset argument", {
  off <- rep(0, nrow(dat))
  expect_no_error(
    simplexreg(y ~ x1 + x2 | z1, data = dat, link.mu = "logit", offset = off)
  )
})

test_that("simplexreg stores design matrices when x = TRUE", {
  fit_x <- simplexreg(y ~ x1 + x2 | z1, data = dat, link.mu = "logit", x = TRUE)
  expect_true(!is.null(fit_x$x))
  expect_true(!is.null(fit_x$x$mean))
  expect_true(!is.null(fit_x$x$dispersion))
})

# test_that("simplexreg fits null model with parametric link", {
#   expect_no_error(
#     simplexreg(y ~ 1, data = dat, link.mu = "plogit1")
#   )
# })
#
# test_that("simplexreg null model with parametric link has finite loglik", {
#   fit_null_par <- simplexreg(y ~ 1, data = dat, link.mu = "plogit1")
#   expect_true(is.finite(fit_null_par$loglik))
#   expect_true(!is.null(fit_null_par$coefficients$lambda))
# })

# ==============================================================================
# 20. simplexreg.fit — low-level interface
# ==============================================================================

test_that("simplexreg.fit returns an object of class simplexregression", {
  X_mat <- cbind(x1, x2)
  Z_mat <- matrix(0, nrow = n, ncol = 0)
  fit_low <- simplexreg.fit(y, X_mat, Z_mat, link.mu = "logit",
                            link.sigma2 = "log")
  expect_s3_class(fit_low, "simplexregression")
})

test_that("simplexreg.fit produces finite log-likelihood", {
  X_mat <- cbind(x1, x2)
  Z_mat <- matrix(0, nrow = n, ncol = 0)
  fit_low <- simplexreg.fit(y, X_mat, Z_mat, link.mu = "logit",
                            link.sigma2 = "log")
  expect_true(is.finite(fit_low$loglik))
})

test_that("simplexreg.fit results are consistent with simplexreg", {
  X_mat <- cbind(x1, x2)
  Z_mat <- cbind(z1)
  fit_low <- simplexreg.fit(y, X_mat, Z_mat, link.mu = "logit",
                            link.sigma2 = "log",
                            x_names = c("(Intercept)", "x1", "x2"),
                            z_names = c("(Intercept)", "z1"))
  # Log-likelihoods should be very close
  expect_equal(fit_low$loglik, fit_logit$loglik, tolerance = 1e-3)
})

test_that("simplexreg.fit handles NULL x (intercept-only mean model)", {
  expect_no_error(
    simplexreg.fit(y, x = NULL, z = NULL,
                   link.mu = "logit", link.sigma2 = "log")
  )
})

test_that("simplexreg.fit accepts scalar weights", {
  X_mat <- cbind(x1, x2)
  Z_mat <- cbind(z1)
  expect_no_error(
    simplexreg.fit(y, X_mat, Z_mat, weights = 1,
                   link.mu = "logit", link.sigma2 = "log")
  )
})

test_that("simplexreg.fit accepts data.frame as x", {
  X_df <- as.data.frame(cbind(x1, x2))
  Z_mat <- cbind(z1)
  expect_no_error(
    simplexreg.fit(y, X_df, Z_mat, link.mu = "logit", link.sigma2 = "log")
  )
})

test_that("simplexreg.fit accepts data.frame as z", {
  X_mat <- cbind(x1, x2)
  Z_df  <- as.data.frame(cbind(z1))
  expect_no_error(
    simplexreg.fit(y, X_mat, Z_df, link.mu = "logit", link.sigma2 = "log")
  )
})

test_that("simplexreg.fit accepts scalar offset", {
  X_mat <- cbind(x1, x2)
  Z_mat <- cbind(z1)
  expect_no_error(
    simplexreg.fit(y, X_mat, Z_mat, offset = 1,
                   link.mu = "logit", link.sigma2 = "log")
  )
})

test_that("simplexreg.fit rejects offset with wrong length", {
  X_mat <- cbind(x1, x2)
  Z_mat <- cbind(z1)
  expect_error(
    simplexreg.fit(y, X_mat, Z_mat, offset = rep(0, n - 1),
                   link.mu = "logit", link.sigma2 = "log"),
    "offset must have length"
  )
})

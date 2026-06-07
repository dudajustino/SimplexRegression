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

# --- print ---

test_that("print.simplexregression returns the object invisibly", {
  result <- withVisible(print(fit_plogit1))
  expect_false(result$visible)
  expect_s3_class(result$value, "simplexregression")
})

# --- summary ---

test_that("summary.simplexregression returns a summary object", {
  s <- summary(fit_plogit1)
  expect_s3_class(s, "summary.simplexregression")
})

test_that("summary.simplexregression contains required components", {
  s <- summary(fit_plogit1)
  expect_true(!is.null(s$coefficients$mean))
  expect_true(!is.null(s$coefficients$dispersion))
  expect_true(is.numeric(s$loglik))
  expect_true(is.numeric(s$nobs))
})

test_that("summary.simplexregression includes lambda for parametric model", {
  s <- summary(fit_plogit1)
  expect_true(s$parametric)
  expect_true(!is.null(s$coefficients$lambda))
})

test_that("summary.simplexregression lambda is NULL for fixed link model", {
  s <- summary(fit_logit)
  expect_false(s$parametric)
  expect_null(s$coefficients$lambda)
})

# --- coef ---

test_that("coef.simplexregression returns named numeric vector", {
  cf <- coef(fit_plogit1)
  expect_true(is.numeric(cf))
  expect_true(!is.null(names(cf)))
})

test_that("coef.simplexregression model = 'mean' returns mean coefficients only", {
  cf <- coef(fit_plogit1, model = "mean")
  expect_equal(length(cf), ncol(fit_plogit1$mu.x))
})

test_that("coef.simplexregression model = 'dispersion' returns dispersion coefficients only", {
  cf <- coef(fit_plogit1, model = "dispersion")
  expect_equal(length(cf), ncol(fit_plogit1$sigma2.x))
})

test_that("coef.simplexregression full model includes lambda for parametric", {
  cf <- coef(fit_plogit1, model = "full")
  expect_true("lambda" %in% names(cf))
})

test_that("coef.simplexregression full model for fixed link has no lambda", {
  cf <- coef(fit_logit, model = "full")
  expect_false("lambda" %in% names(cf))
})

# --- vcov ---

test_that("vcov.simplexregression returns a square matrix", {
  vc <- vcov(fit_plogit1)
  expect_true(is.matrix(vc))
  expect_equal(nrow(vc), ncol(vc))
})

test_that("vcov.simplexregression is positive definite (all eigenvalues > 0)", {
  vc <- vcov(fit_plogit1)
  expect_true(all(eigen(vc)$values > 0))
})

test_that("vcov.simplexregression model = 'mean' returns mean submatrix", {
  vc <- vcov(fit_plogit1, model = "mean")
  p  <- ncol(fit_plogit1$mu.x)
  expect_equal(dim(vc), c(p, p))
})

test_that("vcov.simplexregression dispersion model has correct names", {
  vc <- vcov(fit_plogit1, model = "dispersion")
  expect_equal(colnames(vc), names(fit_plogit1$coefficients$dispersion))
})

# --- logLik ---

test_that("logLik.simplexregression returns a logLik object", {
  ll <- logLik(fit_plogit1)
  expect_s3_class(ll, "logLik")
  expect_true(is.finite(as.numeric(ll)))
})

test_that("logLik.simplexregression df attribute equals number of parameters", {
  ll <- logLik(fit_plogit1)
  p  <- length(coef(fit_plogit1, model = "mean"))
  q  <- length(coef(fit_plogit1, model = "dispersion"))
  expect_equal(attr(ll, "df"), p + q + 1L)  # +1 for lambda
})

# --- fitted ---

test_that("fitted.simplexregression returns numeric vector of length n", {
  fv <- fitted(fit_plogit1)
  expect_true(is.numeric(fv))
  expect_length(fv, n)
  expect_true(all(fv > 0 & fv < 1))
})

# --- predict ---

test_that("predict.simplexregression with no newdata equals fitted values", {
  expect_equal(predict(fit_plogit1, type = "response"), fitted(fit_plogit1))
})

test_that("predict.simplexregression with no newdata type = 'dispersion'", {
  expect_equal(predict(fit_plogit1, type = "dispersion"), fit_plogit1$sigma2.fv)
})

test_that("predict.simplexregression with no newdata type = 'link'", {
  pred <- predict(fit_plogit1, type = "link")
  expect_type(pred, "list")
  expect_named(pred, c("mean", "dispersion"))
})

test_that("predict.simplexregression with newdata returns correct length", {
  newdat <- data.frame(x1 = c(0.2, 0.5, 0.8), x2 = c(0.3, 0.7, 0.9),
                       z1 = c(0.1, 0.4, 0.7))
  pred <- predict(fit_logit, newdata = newdat, type = "response")
  expect_length(pred, 3)
  expect_true(all(pred > 0 & pred < 1))
})

test_that("predict.simplexregression with newdata type = 'dispersion'", {
  newdat <- data.frame(x1 = c(0.2, 0.5), x2 = c(0.3, 0.7), z1 = c(0.1, 0.4))
  pred <- predict(fit_logit, newdata = newdat, type = "dispersion")
  expect_length(pred, 2)
  expect_true(all(pred > 0))
})

test_that("predict.simplexregression with newdata type = 'link'", {
  newdat <- data.frame(x1 = c(0.2, 0.5), x2 = c(0.3, 0.7), z1 = c(0.1, 0.4))
  pred <- predict(fit_logit, newdata = newdat, type = "link")
  expect_type(pred, "list")
  expect_length(pred$mean, 2)
  expect_length(pred$dispersion, 2)
})

test_that("predict.simplexregression with newdata returns parametric response", {
  # covers parametric_mean_link_inv branch with newdata
  newdat <- data.frame(x1 = c(0.2, 0.5), x2 = c(0.3, 0.7), z1 = c(0.1, 0.4))
  pred <- predict(fit_plogit1, newdata = newdat, type = "response")
  expect_length(pred, 2)
  expect_true(all(pred > 0 & pred < 1))
})

test_that("predict.simplexregression works for model without dispersion covariates", {
  # covers else branch: no pipe → z_new = matrix(1, ...)
  fit_nodispers <- simplexreg(y ~ x1 + x2, data = data, link.mu = "logit")
  newdat <- data.frame(x1 = c(0.2, 0.5), x2 = c(0.3, 0.7))
  pred <- predict(fit_nodispers, newdata = newdat, type = "response")
  expect_length(pred, 2)
  expect_true(all(pred > 0 & pred < 1))
})

test_that("predict.simplexregression works with intercept-only dispersion '| 1'", {
  # covers the trimws(parts[2]) == '1' branch → z_new = matrix(1, ...)
  fit_disp1 <- simplexreg(y ~ x1 + x2 | 1, data = data, link.mu = "logit")
  newdat <- data.frame(x1 = c(0.2, 0.5), x2 = c(0.3, 0.7))
  pred <- predict(fit_disp1, newdata = newdat, type = "response")
  expect_length(pred, 2)
  expect_true(all(pred > 0 & pred < 1))
})

# --- nobs ---

test_that("nobs.simplexregression returns correct number of observations", {
  expect_equal(nobs(fit_plogit1), n)
})

# --- deviance ---

test_that("deviance.simplexregression returns a positive finite scalar", {
  d <- deviance(fit_plogit1)
  expect_length(d, 1)
  expect_true(is.finite(d))
  expect_true(d > 0)
})

test_that("deviance.simplexregression works when y is NULL", {
  fit_no_y <- fit_plogit1
  fit_no_y$y <- NULL
  d <- deviance(fit_no_y)
  expect_true(is.finite(d))
  expect_true(d > 0)
})

# --- AIC / BIC / HQIC ---

test_that("AIC.simplexregression returns a finite scalar for a single model", {
  a <- AIC(fit_plogit1)
  expect_length(a, 1)
  expect_true(is.finite(a))
})

test_that("AIC.simplexregression returns a data.frame for multiple models", {
  result <- AIC(fit_plogit1, fit_logit)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("df", "AIC"))
  expect_equal(nrow(result), 2)
})

test_that("BIC.simplexregression returns a finite scalar for a single model", {
  b <- BIC(fit_plogit1)
  expect_length(b, 1)
  expect_true(is.finite(b))
})

test_that("HQIC.simplexregression returns a finite scalar for a single model", {
  h <- HQIC(fit_plogit1)
  expect_length(h, 1)
  expect_true(is.finite(h))
})

test_that("BIC >= AIC for the same model (penalty is larger)", {
  expect_true(BIC(fit_plogit1) >= AIC(fit_plogit1))
})

test_that("BIC.simplexregression returns data.frame for multiple models", {
  result <- BIC(fit_plogit1, fit_logit)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("df", "BIC"))
  expect_equal(nrow(result), 2)
})

test_that("HQIC.simplexregression returns data.frame for multiple models", {
  result <- HQIC(fit_plogit1, fit_logit)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("df", "HQIC"))
  expect_equal(nrow(result), 2)
})

# --- hatvalues ---

test_that("hatvalues.simplexregression returns numeric vector of length n", {
  hv <- hatvalues(fit_plogit1)
  expect_true(is.numeric(hv))
  expect_length(hv, n)
})

test_that("hatvalues.simplexregression values are in (0, 1)", {
  hv <- hatvalues(fit_plogit1)
  expect_true(all(hv > 0 & hv < 1))
})

# --- cooks.distance ---

test_that("cooks.distance.simplexregression returns non-negative numeric vector", {
  cd <- cooks.distance(fit_plogit1)
  expect_true(is.numeric(cd))
  expect_length(cd, n)
  expect_true(all(cd >= 0))
})

test_that("cooks.distance.simplexregression works with type = 'weighted'", {
  cd <- cooks.distance(fit_plogit1, type = "weighted")
  expect_true(is.numeric(cd))
  expect_length(cd, n)
})

test_that("cooks.distance type pearson returns finite values", {
  cd <- cooks.distance(fit_plogit1, type = "pearson")
  expect_true(all(is.finite(cd)))
  expect_true(all(cd >= 0))
})

# --- simulate ---

test_that("simulate.simplexregression returns a data.frame with correct dimensions", {
  sim <- simulate(fit_plogit1, nsim = 3, seed = 1)
  expect_s3_class(sim, "data.frame")
  expect_equal(nrow(sim), n)
  expect_equal(ncol(sim), 3)
})

test_that("simulate.simplexregression simulated values are in (0, 1)", {
  sim <- simulate(fit_plogit1, nsim = 1, seed = 1)
  expect_true(all(sim[[1]] > 0 & sim[[1]] < 1))
})

test_that("simulate.simplexregression works with seed = NULL", {
  sim <- simulate(fit_plogit1, nsim = 2, seed = NULL)
  expect_s3_class(sim, "data.frame")
  expect_equal(ncol(sim), 2)
})

test_that("simulate.simplexregression initializes RNG when .Random.seed does not exist", {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old_seed <- .Random.seed
    rm(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv))
  }
  sim <- simulate(fit_plogit1, nsim = 1, seed = NULL)
  expect_s3_class(sim, "data.frame")
  expect_true(all(sim[[1]] > 0 & sim[[1]] < 1))
})

# --- formula / terms / model.matrix ---

test_that("formula.simplexregression returns the model formula", {
  f <- formula(fit_plogit1)
  expect_true(inherits(f, "formula") || inherits(f, "Formula"))
})

test_that("model.matrix.simplexregression returns the mean design matrix", {
  mm <- model.matrix(fit_plogit1, model = "mean")
  expect_equal(nrow(mm), n)
  expect_equal(ncol(mm), ncol(fit_plogit1$mu.x))
})

test_that("terms.simplexregression defaults to mean model", {
  trms <- terms(fit_plogit1)
  expect_true(!is.null(trms))
})

# --- update ---
test_that("update.simplexregression refits the model with new arguments", {
  fit_updated <- update(fit_logit, link.mu = "probit")
  expect_s3_class(fit_updated, "simplexregression")
  expect_equal(fit_updated$mu.link, "probit")
})

test_that("update.simplexregression errors when object has no call component", {
  obj_no_call <- fit_logit
  obj_no_call$call <- NULL
  expect_error(
    update.simplexregression(obj_no_call),
    "need an object with call component"
  )
})

test_that("update.simplexregression works with evaluate = FALSE", {
  result <- update(fit_logit, evaluate = FALSE)
  expect_true(is.call(result))
})

test_that("update.simplexregression updates formula", {
  fit_updated <- update(fit_logit, formula. = . ~ x1)
  expect_s3_class(fit_updated, "simplexregression")
})

test_that("update.simplexregression adds new argument not in original call", {
  fit_updated <- update(fit_logit, link.sigma2 = "sqrt")
  expect_s3_class(fit_updated, "simplexregression")
  expect_equal(fit_updated$sigma2.link, "sqrt")
})

# --- terms ---

test_that("terms.simplexregression works for mean model", {
  trms <- terms(fit_plogit1, model = "mean")
  expect_true(!is.null(trms))
})

test_that("terms.simplexregression works for dispersion model", {
  trms <- terms(fit_plogit1, model = "dispersion")
  expect_true(!is.null(trms))
})

# --- model.frame ---

test_that("model.frame.simplexregression returns stored model frame", {
  fit_with_model <- simplexreg(y ~ x1 + x2 | z1, data = data,
                               link.mu = "logit", model = TRUE)
  mf <- model.frame(fit_with_model)
  expect_s3_class(mf, "data.frame")
  expect_equal(nrow(mf), n)
})

test_that("model.frame.simplexregression reconstructs frame when not stored", {
  fit_no_model <- simplexreg(y ~ x1 + x2 | z1, data = data,
                             link.mu = "logit", model = FALSE)
  expect_null(fit_no_model$model)
  mf <- model.frame(fit_no_model)
  expect_s3_class(mf, "data.frame")
  expect_equal(nrow(mf), n)
})

# --- model.matrix for dispersion ---

test_that("model.matrix.simplexregression returns dispersion design matrix", {
  mm <- model.matrix(fit_plogit1, model = "dispersion")
  expect_equal(nrow(mm), n)
  expect_equal(ncol(mm), ncol(fit_plogit1$sigma2.x))
})

# --- df.residual ---

test_that("df.residual.simplexregression returns correct value", {
  dfr <- df.residual(fit_plogit1)
  expect_true(is.numeric(dfr))
  expect_equal(dfr, fit_plogit1$df.residual)
})

# --- print.summary ---

test_that("print.summary.simplexregression works without errors", {
  s <- summary(fit_plogit1)
  expect_output(print(s), "Call:")
  expect_output(print(s), "Coefficients")
})

# --- bread (sandwich) ---

test_that("bread.simplexregression returns a matrix", {
  skip_if_not_installed("sandwich")
  br <- bread(fit_plogit1)
  expect_true(is.matrix(br))
  expect_equal(dim(br), dim(vcov(fit_plogit1)))
})

test_that("bread.simplexregression returns correct dimension", {
  skip_if_not_installed("sandwich")
  br <- bread(fit_plogit1)
  npar <- length(coef(fit_plogit1, model = "full"))
  expect_equal(dim(br), c(npar, npar))
})

# --- estfun (sandwich) ---

test_that("estfun.simplexregression returns a matrix", {
  skip_if_not_installed("sandwich")
  ef <- estfun(fit_logit)
  expect_true(is.matrix(ef))
  expect_equal(nrow(ef), n)
})

test_that("estfun.simplexregression has correct column names", {
  skip_if_not_installed("sandwich")
  ef <- estfun(fit_plogit1)
  expected_names <- names(coef(fit_plogit1, model = "full"))
  expect_equal(colnames(ef), expected_names)
})

test_that("estfun.simplexregression works with plogit2 parametric link", {
  skip_if_not_installed("sandwich")
  fit_plogit2 <- simplexreg(y ~ x1 + x2 | z1, data = data, link.mu = "plogit2")
  ef <- estfun(fit_plogit2)
  expect_true(is.matrix(ef))
  expect_equal(nrow(ef), n)
  expect_equal(ncol(ef), length(coef(fit_plogit2, model = "full")))
})

# --- coeftest (lmtest) ---

test_that("coeftest.simplexregression works", {
  skip_if_not_installed("lmtest")
  ct <- coeftest(fit_plogit1)
  expect_s3_class(ct, "coeftest")
})

# --- lrtest (lmtest) ---

test_that("lrtest.simplexregression works", {
  skip_if_not_installed("lmtest")
  lt <- lrtest(fit_plogit1, fit_logit)
  expect_s3_class(lt, "anova")
})

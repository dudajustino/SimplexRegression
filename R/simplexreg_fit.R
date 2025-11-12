################################################################################
#                 SIMPLEX REGRESSION - MAIN FITTING FUNCTIONS                  #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2025-11-08                                                             #
# Description: Main functions for fitting simplex regression models with       #
#              parametric or fixed mean link functions    .                    #
################################################################################

# ==============================================================================
# 1. MAIN USER-FACING FUNCTION
# ==============================================================================

#' @title Simplex Regression with Parametric or Fixed Mean Link
#' @description Fits a simplex regression model with parametric or fixed mean link
#' functions.
#'
#' @param formula A two-part formula: y ~ mean_predictors | dispersion_predictors
#' @param data A data frame containing the variables in the model
#' @param diag Diagnostic level (default = 1)
#' @param link.mu Mean link function (parametric: "plogit1", "plogit2"; fixed:
#' "logit", "probit", "loglog", "cloglog", "cauchit")
#' @param link.sigma2 Dispersion link function ("log", "sqrt", "identity")
#'
#' @importFrom stats model.frame model.response model.matrix terms plogis qnorm pnorm cor
#'
#' @return An object of class \code{"simplexregression"} containing:
#' \itemize{
#'   \item \code{coefficients}: List with mean, dispersion, and lambda estimates
#'   \item \code{fitted.values}: Fitted mean values
#'   \item \code{sigma2.fv}: Fitted dispersion values
#'   \item \code{mu.lp}: Predicted linear mean
#'   \item \code{sigma2.lp}: Predicted linear dispersion
#'   \item \code{optim}: List with initial values and convergence of optim
#'   \item \code{residuals}: Quantile residuals
#'   \item \code{vcov}: Variance-covariance matrix
#'   \item \code{loglik}: Log-likelihood value
#'   \item \code{AIC}, \code{BIC}, \code{HQ}: Information criteria
#'   \item \code{R2_N}, \code{R2_FC}: Pseudo R-squared measures
#'
#' }
#'
#' @details
#' The model is specified using a two-part formula:
#' \itemize{
#'   \item Mean submodel: \eqn{g(\mu_i, \lambda) = x_i'\beta}
#'   or \eqn{g(\mu_i) = x_i'\beta}
#'   \item Dispersion submodel: \eqn{h(\sigma^2_i) = z_i'\gamma}
#' }
#'
#' The parametric mean link functions include a parameter \eqn{\lambda} that
#' is estimated along with other model parameters.
#'
#' @examples
#' # Simulate data
#' n <- 100
#' x1 <- runif(n, 0, 1)
#' x2 <- runif(n, 0, 1)
#' mu <- parametric_mean_link_inv(0.8 - 1.2*x1 - 1.5*x2 , 0.25, "plogit2")
#' sigma2 <- 0.5
#' summary(mu)
#' y <- rsimplex_opt(n, mu, sigma2)
#' data <- data.frame(y = y, x = cbind(x1, x2))
#'
#' # Fit model
#' fit <- simplexreg(y ~ x1 + x2 | 1, data = data,
#'                      link.mu = "plogit2", link.sigma2 = "identity")
#' summary(fit)
#'
#' @references
#' Justino, M. E. C. and Cribari-Neto, F. (2025).
#' Simplex regression with a flexible logit link: Inference and application
#' to cross-country impunity data.
#' \emph{Working Paper}, Departamento de Estatística,
#' Universidade Federal de Pernambuco.
#' \url{https://github.com/dudajustino/SimplexRegressionAMM}
#'
#' @export
simplexreg <- function(formula, data, diag = 1,
                          link.mu = c("plogit1", "plogit2", "logit", "probit",
                                      "loglog", "cloglog", "cauchit"),
                          link.sigma2 = c("log", "sqrt", "identity")) {

  # Input validation
  if (missing(data)) {
    stop("The 'data' argument is required")
  }

  if (!inherits(formula, "formula")) {
    stop("The 'formula' argument must be a formula object")
  }

  # Check for separator |
  formula_str <- as.character(formula)
  if (length(formula_str) < 3 || !grepl("\\|", formula_str[3])) {
    stop("Formula must have the format 'response ~ mean_predictors | dispersion_predictors'")
  }

  # Extract formula components
  parts <- strsplit(formula_str[3], "\\|")[[1]]
  formula_mean <- as.formula(paste(formula_str[2], formula_str[1], parts[1]))
  formula_disp <- as.formula(paste("~", parts[2]))

  # Check if null model (intercept only in both parts)
  is_null_model <- (trimws(parts[1]) == "1" && trimws(parts[2]) == "1")

  # Extract design matrices
  mf_mean <- model.frame(formula_mean, data = data)
  y <- model.response(mf_mean)
  x_matrix <- model.matrix(formula_mean, data = data)

  # Create x matrix without intercept
  if (ncol(x_matrix) > 1) {
    x <- x_matrix[, -1, drop = FALSE]
  } else {
    x <- matrix(0, nrow = length(y), ncol = 0)
  }

  # Save column names (includes "(Intercept)")
  x_names <- colnames(x_matrix)

  mf_disp <- model.frame(formula_disp, data = data)
  z_matrix <- model.matrix(formula_disp, data = data)

  # Create z matrix without intercept
  if (ncol(z_matrix) > 1) {
    z <- z_matrix[, -1, drop = FALSE]
  } else {
    z <- matrix(0, nrow = length(y), ncol = 0)
  }

  # Save column names (includes "(Intercept)")
  z_names <- colnames(z_matrix)

  link.mu <- match.arg(link.mu)
  link.sigma2 <- match.arg(link.sigma2)
  parametric <- link.mu %in% c("plogit1","plogit2")

  # Handle null model with special initialization
  if (is_null_model) {

    if(parametric) {
      ystar <- log(y / (1 - y))
      betaols_null <- mean(ystar)
      muols_null <- exp(betaols_null) / (1 + exp(betaols_null))
      devtrans_null <- (y - muols_null)^2 / (y * (1 - y) * muols_null^2 * (1 - muols_null)^2)
      deltaols_null <- sum(devtrans_null)/length(y)
      ini_nul <- c(betaols_null, deltaols_null, 1)
    } else {
      ystar <- fixed_mean_link(y, link.mu)
      betaols_null <- mean(ystar)
      muols_null <- fixed_mean_link_inv(betaols_null, link.mu)
      devtrans_null <- (y - muols_null)^2 / (y * (1 - y) * muols_null^2 * (1 - muols_null)^2)
      deltaols_null <- sum(devtrans_null)/length(y)
      ini_nul <- c(betaols_null, deltaols_null)
    }

    result <- simplexreg_fit(y, x, z, diag = 0,
                                link.mu = link.mu, link.sigma2 = link.sigma2,
                                x_names = x_names, z_names = z_names, ini_values = ini_nul)

  } else {
    # Call standard fitting function
    result <- simplexreg_fit(y, x, z, diag = 0,
                                link.mu = link.mu, link.sigma2 = link.sigma2,
                                x_names = x_names, z_names = z_names)
  }

  # Add formula information to result
  result$formula <- formula
  result$formula_mean <- formula_mean
  result$formula_disp <- formula_disp
  result$terms <- list(mean = terms(formula_mean), dispersion = terms(formula_disp))
  result$model <- mf_mean
  result$call <- match.call()

  return(result)
}

# ==============================================================================
# 2. CORE FITTING FUNCTION
# ==============================================================================

#' @title Simplex Regression with Parametric or Fixed Mean Link
#' @description Internal function that performs the actual model fitting.
#'
#' @param y Response variable (numeric vector, 0 < y < 1)
#' @param x Design matrix for mean model (without intercept)
#' @param z Design matrix for dispersion model (without intercept)
#' @param diag Diagnostic level
#' @param link.mu Mean link function (parametric: "plogit1", "plogit2"; fixed:
#' "logit", "probit", "loglog", "cloglog", "cauchit")
#' @param link.sigma2 Dispersion link function ("log", "sqrt", "identity")
#' @param x_names Column names for mean design matrix (includes intercept)
#' @param z_names Column names for dispersion design matrix (includes intercept)
#' @param ini_values Optional initial values for optimization
#'
#' @importFrom stats plogis qnorm pnorm cor optim lm.fit
#' @return A list with model fitting results
#' @keywords internal
simplexreg_fit <- function(y, x, z, diag = 1,
                              link.mu = c("plogit1", "plogit2", "logit", "probit",
                                           "loglog", "cloglog", "cauchit"),
                              link.sigma2 = c("log", "sqrt", "identity"),
                              x_names = NULL, z_names = NULL, ini_values = NULL){

  # Ensure x is a matrix
  if (is.null(x)) {
    x <- matrix(0, nrow = length(y), ncol = 0)
  } else if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  # Create x1 (with intercept)
  if (ncol(x) > 0) {
    x1 <- cbind(rep(1, length(y)), x)
  } else {
    x <- matrix(0, nrow = length(y), ncol = 0)
    x1 <- matrix(1, nrow = length(y), ncol = 1)
  }

  # Ensure z is a matrix
  if (is.null(z)) {
    z <- matrix(0, nrow = length(y), ncol = 0)
  } else if (!is.matrix(z)) {
    z <- as.matrix(z)
  }

  # Create z1 (with intercept)
  if (ncol(z) > 0) {
    z1 <- cbind(rep(1, length(y)), z)
  } else {
    z <- matrix(0, nrow = length(y), ncol = 0)
    z1 <- matrix(1, nrow = length(y), ncol = 1)
  }

  link.mu <- match.arg(link.mu)
  link.sigma2 <- match.arg(link.sigma2)
  parametric <- link.mu %in% c("plogit1","plogit2")

  y <- as.vector(y)

  # Model dimensions
  p <- ncol(x1)  # Number of mean parameters (including intercept)
  q <- ncol(z1)  # Number of dispersion parameters (including intercept)
  if(parametric) r <- (p + q + 1) else r <- (p + q)
  n <- length(y)

  # ============================
  # LOG-LIKELIHOOD FUNCTION
  # ============================
  loglik <- function(theta){
    beta <- theta[1:p]
    delta <- theta[(p+1):(p+q)]

    eta1 <- as.vector(x1%*%beta)
    eta2 <- as.vector(z1%*%delta)

    # Soft safeguards for eta2
    if (link.sigma2 == "log") {
      eta2 <- pmin(pmax(eta2, -20), 20)
    } else if (link.sigma2 == "sqrt") {
      eta2 <- pmax(eta2, 0.01)
    } else if (link.sigma2 == "identity") {
      eta2 <- pmax(eta2, 1e-6)
    }

    if(parametric){
      lambda <- theta[r]
      mu <- as.vector(parametric_mean_link_inv(eta1, lambda, link.mu))
      lambda <- max(lambda, 0.001)
    } else {
      mu <- as.vector(fixed_mean_link_inv(eta1, link.mu))
    }

    sigma2 <- as.vector(dispersion_link_inv(eta2, link.sigma2))

    # Standard safeguards
    y <- pmin(pmax(y, 1e-6), 1 - 1e-6)
    mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
    sigma2 <- pmax(sigma2, 1e-6)

    dev <- deviance_simplexreg(y, mu)

    # Check for non-finite values
    if (any(!is.finite(c(mu, sigma2, dev)))) {
      warning("Non-finite values in log-likelihood estimates.")
      return(-1e10)
    }

    # Log-likelihood function
    adFunc <- -0.5 * sum(log(2 * pi) + log(sigma2) + 3 * log(y * (1 - y)) + dev / sigma2)

    if (!is.finite(adFunc)) {
      warning("Non-finite log-likelihood value.")
      return(-1e10)
    }

    return(adFunc)
  }

  # ============================
  # SCORE FUNCTION
  # ============================
  escore <- function(theta){
    beta <- theta[1:p]
    delta <- theta[(p+1):(p+q)]

    eta1 <- as.vector(x1%*%beta)
    eta2 <- as.vector(z1%*%delta)

    #  Soft safeguards for eta2
    if (link.sigma2 == "log") {
      eta2 <- pmin(pmax(eta2, -20), 20)
    } else if (link.sigma2 == "sqrt") {
      eta2 <- pmax(eta2, 0.01)
    } else if (link.sigma2 == "identity") {
      eta2 <- pmax(eta2, 1e-6)
    }

    if(parametric){
      lambda <- theta[r]
      mu <- parametric_mean_link_inv(eta1, lambda, link.mu)
      lambda <- max(lambda, 0.001)
    } else {
      mu <- fixed_mean_link_inv(eta1, link.mu)
    }
    sigma2 <- as.vector(dispersion_link_inv(eta2, link.sigma2))

    y <- pmin(pmax(y, 1e-6), 1 - 1e-6)
    mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
    sigma2 <- pmax(sigma2, 1e-6)

    diff <- as.vector(y - mu)
    muonemu <- as.vector(mu * (1 - mu))

    dev <- deviance_simplexreg(y, mu)

    # Check for non-finite values
    if (any(!is.finite(c(mu, sigma2, diff, muonemu, dev)))) {
      warning("Non-finite values in score vector estimates.")
    }

    U <- diag(as.vector(dev / (sigma2 * muonemu) + 1 / (sigma2 * (muonemu)^3)))
    a <- as.vector((-1 / (2 * sigma2)) + (dev / (2 * sigma2^2)))

    T2 <- diag(as.vector(dispersion_link_inv_deriv1(eta2, link.sigma2)))

    if(parametric){
      T1 <- diag(as.vector(parametric_mean_link_inv_deriv1(eta1, lambda, link.mu)))

      # Compute rho based on link function
      if(link.mu == "plogit2"){
        exp_aval_frac <- plogis(eta1)^(1/lambda)
        rho <- as.vector(- exp_aval_frac * (log(plogis(eta1))) / (lambda^2))
      } else{
        rho <- as.vector((-1/(lambda^2)) * ((exp(eta1) + 1) ^ (-1/lambda)) * (log(exp(eta1) + 1)))
      }

      if (any(!is.finite(rho))) {
        warning("Non-finite values in rho.")
      }

      Ubeta <- t(x1) %*% T1 %*% U %*% diff
      Udelta <- t(z1) %*% T2 %*% a
      Ulambda <- t(rho) %*% U %*% diff
      avScore <- c(Ubeta, Udelta, Ulambda)

    } else {
      T1 <- diag(as.vector(fixed_mean_link_inv_deriv1(eta1, link.mu)))
      Ubeta <- t(x1) %*% T1 %*% U %*% diff
      Udelta <- t(z1) %*% T2 %*% a
      avScore <- c(Ubeta, Udelta)
    }

    return(as.vector(avScore))
  }

  # ============================
  # INITIAL VALUES
  # ============================
  if (!is.null(ini_values)) {
    ini <- ini_values
  } else {
    # Standard initialization with OLS
    if(parametric){
      ystar <- log(y / (1 - y))
      betaols <- (lm.fit(x1, ystar))$coefficients
      muols <- exp(x1 %*% betaols) / (1 + exp(x1 %*% betaols))
      devtrans <- (y - muols)^2 / (y * (1 - y) * muols^2 * (1 - muols)^2)
      deltaols <- (lm.fit(z1, dispersion_link(devtrans, link.sigma2)))$coefficients

      lambdaini <- 1
      ini <- c(as.vector(betaols), as.vector(deltaols), as.numeric(lambdaini))
    } else {
      ystar <- fixed_mean_link(y, link.mu)
      betaols <- (lm.fit(x1, ystar))$coefficients
      muols <- fixed_mean_link_inv(x1 %*% betaols, link.mu)
      devtrans <- (y - muols)^2 / (y * (1 - y) * muols^2 * (1 - muols)^2)
      deltaols <- (lm.fit(z1, dispersion_link(devtrans, link.sigma2)))$coefficients

      ini <- c(as.vector(betaols), as.vector(deltaols))
    }

    # ============================
    # OPTIMIZATION
    # ============================
    opt <- optim(ini, loglik, escore, method = "BFGS",
                 control=list(fnscale = -1, maxit = 500, reltol=1e-8))

    if (opt$convergence != 0) warning("FUNCTION DID NOT CONVERGE!")
  }

  k <- list()

  # Store results
  coef <- opt$par
  k$coef <- coef
  k$conv <- opt$convergence
  k$loglik <- opt$value
  k$counts <- as.numeric(opt$counts[1])
  k$betas <- coef[1:p]
  k$deltas <- coef[(p+1):(p+q)]
  if(parametric) k$lambda <- coef[r] else k$lambda <- NA
  k$y <- y
  k$mu.x <- x1
  k$sigma2.x <- z1

  eta1_hat <- as.vector(x1%*%(k$betas))
  eta2_hat <- as.vector(z1%*%(k$deltas))

  k$mu.lp <- eta1_hat
  k$sigma2.lp <- eta2_hat

  if(parametric) {
    mu_hat <- as.vector(parametric_mean_link_inv(eta1_hat, k$lambda, link.mu))
  } else {
    mu_hat <- as.vector(fixed_mean_link_inv(eta1_hat, link.mu))
  }

  sigma2_hat <- as.vector(dispersion_link_inv(eta2_hat, link.sigma2))
  lambda_hat <- k$lambda

  k$mu.fv <- mu_hat
  k$sigma2.fv <- sigma2_hat

  # ============================
  # VARIANCE-COVARIANCE
  # ============================
  # Compute variance-covariance matrix (if not null model)
  if (is.null(ini_values)) {
    if(parametric) {
      l1id1 <- as.vector(parametric_mean_link_inv_deriv1(eta1_hat, lambda_hat, link.mu))

      # Compute rho
      if(link.mu == "plogit2"){
        exp_aval_frac <- plogis(eta1_hat)
        rho <- (-1/(lambda_hat^2)) * (exp_aval_frac ^ (1/lambda_hat)) * (log(exp_aval_frac))
      } else{
        rho <- (-1/(lambda_hat^2)) * ((exp(eta1_hat) + 1) ^ (-1/lambda_hat)) * (log(exp(eta1_hat) + 1))
      }

      if (any(!is.finite(rho))) {
        warning("Non-finite values in rho.")
      }

    } else {
      l1id1 <- as.vector(fixed_mean_link_inv_deriv1(eta1_hat, link.mu))
    }

    l2id1 <- as.vector(dispersion_link_inv_deriv1(eta2_hat, link.sigma2))

    muonemu_hat <- as.vector(mu_hat * (1 - mu_hat))
    wi <- as.vector((3.0 * sigma2_hat) / muonemu_hat + (1 / (muonemu_hat^3)))
    vi <- as.vector(1 / (2 * sigma2_hat^2))

    if (any(!is.finite(mu_hat)) || any(!is.finite(sigma2_hat)) || any(!is.finite(muonemu_hat))) {
      warning("Non-finite values in final estimates.")
    }

    # Information matrix components
    Wbetabeta <- diag(as.vector((1.0 / sigma2_hat) * wi * (l1id1^2)))
    Wdeltadelta <- diag(as.vector(vi * (l2id1^2)))

    Kbetabeta <- t(x1) %*% Wbetabeta %*% x1
    Kdeltadelta <- t(z1) %*% Wdeltadelta %*% z1

    if(parametric){
      Wbetalambda <- diag(as.vector((1.0 / sigma2_hat) * wi * l1id1))
      Wlambdalambda <- diag(as.vector((1.0 / sigma2_hat) * wi))

      Kbetalambda <- t(x1) %*% Wbetalambda %*% rho
      Klambdabeta <- t(Kbetalambda)
      Klambdalambda <- t(rho) %*% Wlambdalambda %*% rho

      # Assemble information matrix
      K <- rbind(
        cbind(Kbetabeta, matrix(0, p, q), Kbetalambda),
        cbind(matrix(0, q, p), Kdeltadelta, matrix(0, q, 1)),
        cbind(Klambdabeta, matrix(0, 1, q), Klambdalambda)
      )
    } else {
      K <- rbind(
        cbind(Kbetabeta, matrix(0, p, q)),
        cbind(matrix(0, q, p), Kdeltadelta)
      )
    }

    k$K <- K

    # Variance-covariance matrix
    vcov <- chol2inv(chol(K))
    k$vcov <- vcov
    stderror <- sqrt(diag(vcov))
    k$stderror <- stderror
  }

  # Test statistics
  k$zstat <- abs(coef[1:(p+q)]/stderror[1:(p+q)])
  k$pvalues <- 2*(1 - pnorm(k$zstat) )

  # Model diagnostics
  k$quantile.res <- qnorm(psimplex_opt(y, mu_hat, sigma2_hat))
  k$R2_N <- 1 - exp((-2/n)*(opt$value - simplexreg.nul(y, link.mu)))
  if(parametric)  gy <- parametric_mean_link(y, k$lambda, link.mu) else gy <- fixed_mean_link(y, link.mu)
  k$R2_FC <- cor(eta1_hat, gy)^2
  k$aic <- -2*k$loglik+2*r
  k$bic <- -2*k$loglik+log(n)*r
  k$hq <- -2*k$loglik+2*r*log(log(n))

  # Construct result object
  result = list(
    y = structure(y, .Names = seq(1:n)),
    coefficients = list(
      mean = structure(k$betas, .Names = seq(1:p)),
      dispersion = structure(k$deltas, .Names = seq(1:q)),
      lambda = k$lambda),
    fitted.values = structure(mu_hat, .Names = seq(1:n)),
    optim = list(start = ini, convergence = k$conv, counts = k$counts),
    mu.fv = structure(k$mu.fv, .Names = seq(1:n)),
    mu.lp = structure(k$mu.lp, .Names = seq(1:n)),
    mu.x = k$mu.x,
    mu.link = link.mu,
    mu.df = length(k$betas),
    sigma2.fv = structure(k$sigma2.fv, .Names = seq(1:n)),
    sigma2.lp = structure(k$sigma2.lp, .Names = seq(1:n)),
    sigma2.x = k$sigma2.x,
    sigma2.link = link.sigma2,
    sigma2.df = length(k$deltas),
    lambda.fv = k$lambda,
    df.residual = n - length(coef),
    nobs = n,
    loglik = k$loglik,
    vcov = k$vcov,
    residuals = structure(k$quantile.res, .Names = seq(1:n)),
    AIC = k$aic,
    BIC = k$bic,
    HQ = k$hq,
    R2_FC = k$R2_FC,
    R2_N = k$R2_N,
    zstat = k$zstat,
    pvalues = k$pvalues,
    call = match.call()
  )

  result$x_names <- x_names
  result$z_names <- z_names

  class(result) <- c("simplexregression", "glm", "lm")

  return(result)
}

# ==============================================================================
# 3. NULL MODEL LOG-LIKELIHOOD
# ==============================================================================

#' @title Null Model Log-Likelihood for Simplex Regression
#' @description Computes the log-likelihood for the null (intercept-only) for a
#' simplex regression with a parametric or fixed mean link.
#'
#' @param y Numeric response vector (0 < y < 1)
#' @param link.mu Mean link function: parametric ("plogit1", "plogit2") or fixed
#'   ("logit", "probit", "loglog", "cloglog", "cauchit")
#'
#' @return Numeric value of the null model log-likelihood
#' @keywords internal
simplexreg.nul <- function(y, link.mu){
  y <- as.vector(y)
  n <- length(y)

  parametric <- link.mu %in% c("plogit1","plogit2")

  # Null model log-likelihood function
  fLogLiknull <- function(theta){
    sigma2 <- theta[2]
    eta1 <- cbind(rep(1, n)) %*% theta[1]

    if(parametric){
      lambda = theta[3]
      mu <- parametric_mean_link_inv(eta1, lambda, link.mu)
    } else {
      mu <- fixed_mean_link_inv(eta1, link.mu)
    }

    y <- pmin(pmax(y, 1e-6), 1 - 1e-6)
    mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
    sigma2 <- pmax(sigma2, 1e-6)

    diff <- as.vector(y - mu)
    yoneminy <- as.vector(y * (1 - y))
    muonemu <- as.vector(mu * (1 - mu))
    dev <- (diff / muonemu)^2 / yoneminy

    # Log-likelihood function
    adFunc <- -0.5 * sum(log(2*pi) + log(sigma2) + 3*log(yoneminy) + dev/sigma2)
    adFunc
  }

  # Initial values
  if(parametric){
    ystar <- log(y / (1 - y))
    betaols_null <- mean(ystar)
    muols_null <- exp(betaols_null) / (1 + exp(betaols_null))
    devtrans_null <- (y - muols_null)^2 / (y * (1 - y) * muols_null^2 * (1 - muols_null)^2)
    deltaols_null <- sum(devtrans_null)/length(y)

    ini_nul <- c(betaols_null, deltaols_null, 1)
  } else {
    ystar <- fixed_mean_link(y, link.mu)
    betaols_null <- mean(ystar)
    muols_null <- fixed_mean_link_inv(betaols_null, link.mu)
    devtrans_null <- (y - muols_null)^2 / (y * (1 - y) * muols_null^2 * (1 - muols_null)^2)
    deltaols_null <- sum(devtrans_null)/length(y)

    ini_nul <- c(betaols_null, deltaols_null)
  }

  # Optimize
  opt_nul <- optim(ini_nul, fLogLiknull, method = "BFGS",
                   control=list(fnscale = -1, maxit = 500, reltol = 1e-8))

  k <- opt_nul$value
  return(k)
}

################################################################################
#                 SIMPLEX REGRESSION - MAIN FITTING FUNCTIONS                  #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2025-11-08                                                             #
# Description: Main functions for fitting simplex regression models with       #
#              parametric or fixed mean link functions.                        #
################################################################################

# ==============================================================================
# 0. CONTROL FUNCTION
# ==============================================================================

#' @title Control Parameters for Simplex Regression
#' @description Auxiliary function for controlling simplex regression fitting.
#'
#' @param method Character string specifying the optimization method
#' passed to \code{optim} (default: \code{"BFGS"}).
#' @param maxit Integer specifying the maximum number of iterations for \code{optim}
#' (default: \code{5000}).
#' @param gradient Logical; use analytical gradient? (default: \code{TRUE}).
#' @param hessian Logical; compute Hessian via \code{optim}? (default: \code{FALSE}).
#' @param trace Logical; trace optimization? (default: \code{FALSE}).
#' @param start An optional vector with starting values for all parameters.
#' @param fsmaxit Integer specifying maximal number of additional Fisher scoring
#' iterations (default: \code{500}).
#' @param fstol Numeric tolerance for convergence in Fisher scoring (default: \code{1e-8}).
#' @param reltol Relative convergence tolerance (default: \code{sqrt(.Machine$double.eps)}).
#' @param ... Additional parameters passed to \code{optim}.
#'
#' @details
#' All parameters in \code{\link{simplexreg}} are estimated by maximum likelihood using
#' \code{\link{optim}} with control options set in \code{\link{simplexreg.control}}.
#' Most arguments are passed on directly to \code{optim}, and \code{start} controls
#' how \code{optim} is called.
#'
#' After the \code{optim} maximization, an additional Fisher scoring iteration can
#' be performed to further enhance the result by moving the gradient even closer to zero.
#' If \code{fsmaxit} is greater than zero, this additional optimization is performed and
#' it converges if the threshold \code{fstol} is attained for the absolute value of the
#' step size.
#'
#' Starting values can be supplied via \code{start} or estimated by \code{\link{lm.wfit}},
#' using the link-transformed response. For parametric mean link functions (\code{"plogit1"},
#' \code{"plogit2"}), the link parameter \eqn{\lambda} is jointly estimated with the regression
#' coefficients. Covariances are derived analytically using the expected Fisher information
#' matrix. The Fisher scoring uses analytical gradients and the expected information matrix to
#' refine the maximum likelihood estimates obtained from \code{optim}.
#'
#' The main parameters of interest are the coefficient vector \eqn{\boldsymbol{\beta}}
#' in the linear predictor of the mean submodel and the coefficient vector
#' \eqn{\boldsymbol{\gamma}} in the linear predictor of the dispersion submodel.
#' For parametric links, the additional link parameter \eqn{\lambda} is also
#' estimated and reported. The dispersion parameter \eqn{\sigma^2} can be modeled either
#' as constant (when the dispersion formula contains only an intercept) or as varying
#' across observations through a linear predictor.
#'
#' @seealso \code{\link{simplexreg}}
#'
#' @return A list of control parameters.
#' @export
simplexreg.control <- function(method = "BFGS",
                               maxit = 5000,
                               gradient = TRUE,
                               hessian = FALSE,
                               trace = FALSE,
                               start = NULL,
                               fsmaxit = 500,
                               fstol = 1e-8,
                               reltol = .Machine$double.eps^(0.5),
                               ...) {

  method <- match.arg(method, c("BFGS", "Nelder-Mead", "CG", "L-BFGS-B", "SANN"))

  rval <- list(
    method = method,
    maxit = maxit,
    gradient = gradient,
    hessian = hessian,
    trace = trace,
    start = start,
    fsmaxit = fsmaxit,
    fstol = fstol,
    reltol = reltol
  )

  rval <- c(rval, list(...))

  # Configuration specific for optim (maximization)
  if(!is.null(rval$fnscale)) warning("fnscale must not be modified")
  rval$fnscale <- -1  # Maximize

  if(method == "L-BFGS-B") {
    rval$reltol <- NULL
    if(is.null(rval$factr)) rval$factr <- 1e7
    if(is.null(rval$pgtol)) rval$pgtol <- 0
  } else {
    if(is.null(rval$reltol)) rval$reltol <- .Machine$double.eps^(0.5)
  }

  return(rval)
}

# ==============================================================================
# 1. MAIN USER-FACING FUNCTION
# ==============================================================================

#' @title Simplex Regression with Parametric or Fixed Mean Link
#' @description Fit simplex regression models for rates and proportions via
#' maximum likelihood estimation, modeling both the mean (via parametric
#' or fixed link function) and the dispersion parameter.
#'
#' @param formula A two-part formula: \code{y ~ x} or \code{y ~ x | 1}
#' (mean submodel, constant dispersion), or \code{y ~ x | z} (submodels for both
#' mean and dispersion).
#' @param data A data frame containing the variables in formula.
#' @param subset A specification of the rows/observations to be used: defaults to all.
#' @param na.action An optional (name of a) function for treating missing values (NAs).
#' @param weights An optional numeric vector of case weights.
#' @param offset Optional numeric vector specifying a known component to be
#' included in the linear predictor of the mean submodel during fitting.
#' One or more \code{\link{offset}} terms can be included in the formula instead
#' or as well, and if more than one is specified their sum is used. See
#' \code{\link{model.offset}}.
#' @param link.mu Character specification of the link function in the mean submodel
#' (parametric functions: \code{"plogit1"}, \code{"plogit2"}; or fixed functions:
#' \code{"logit"}, \code{"probit"}, \code{"loglog"}, \code{"cloglog"}, \code{"cauchit"}).
#' @param link.sigma2 Character specification of the link function in the dispersion
#' submodel (\code{"log"}, \code{"sqrt"}, \code{"identity"}).
#' @param control A list of control arguments specified via \code{\link{simplexreg.control}}.
#' @param contrasts An optional list. See the \code{contrasts.arg} argument of
#' \code{\link[stats]{model.matrix.default}}.
#' @param model,y,x  Logicals. If \code{TRUE} the corresponding components of the fit
#' (model frame, response, model matrix) are returned. For \code{\link{simplexreg.fit}},
#' \code{x} should be a numeric regressor matrix and \code{y} should be the numeric response
#' vector (with values in (0, 1)).
#' @param x_names Column names for mean design matrix (includes intercept).
#' @param z_names Column names for dispersion design matrix (includes intercept).
#' @param z Design matrix for dispersion model (without intercept).
#' @param ... Additional arguments passed to \code{\link{simplexreg.control}}.
#'
#' @name simplexreg.fit
#'
#' @aliases simplexreg simplexreg.fit
#'
#' @details
#' Simplex regression, introduced by Song and Tan (2000) and extended by
#' Song, Qiu and Tan (2004), is useful for modeling continuous response variables
#' restricted to the unit interval (0, 1), such as rates and proportions. The model
#' assumes that the dependent variable follows the simplex distribution, originally
#' proposed by Barndorff-Nielsen and Jørgensen (1991), which is indexed by the mean
#' \eqn{\mu} and the dispersion parameter \eqn{\sigma^2}.
#'
#' Similar to generalized linear models (GLMs), simplex regression relates
#' the mean of the response variable and the dispersion parameter to their respective
#' linear predictors through link functions. This package implements five fixed link
#' functions (\code{"logit"}, \code{"probit"}, \code{"loglog"}, \code{"cloglog"},
#' \code{"cauchit"}) and two parametric link functions (\code{"plogit1"},
#' \code{"plogit2"}) for the mean submodel. For the dispersion submodel, the
#' links \code{"log"}, \code{"sqrt"} and \code{"identity"} are supported.
#'
#' The model is specified through a two-part formula separated by \code{|}. The
#' left side contains the predictors for the mean submodel and the right side contains
#' the predictors for the dispersion submodel:
#' \itemize{
#'   \item Mean submodel: \eqn{g(\mu_i, \lambda) = \boldsymbol{x}_i'\boldsymbol{\beta}}
#'   (parametric link) or \eqn{g(\mu_i) = \boldsymbol{x}_i'\boldsymbol{\beta}} (fixed link),
#'   where \eqn{g(\cdot)} is the mean link function, \eqn{\lambda} is the
#'   extra shape parameter of the parametric link,
#'   \eqn{\boldsymbol{x}_i} is the vector of covariates for the \eqn{i}-th observation in
#'   the mean submodel, and \eqn{\boldsymbol{\beta}} is the corresponding
#'   vector of regression coefficients.
#'   \item Dispersion submodel: \eqn{h(\sigma^2_i) = \boldsymbol{z}_i'\boldsymbol{\gamma}},
#'   where \eqn{h(\cdot)} is the dispersion link function, \eqn{\boldsymbol{z}_i} is the
#'   vector of covariates for the \eqn{i}-th observation in the dispersion
#'   submodel, and \eqn{\boldsymbol{\gamma}} is the corresponding vector of
#'   regression coefficients.
#' }
#'
#' Formula examples: \code{y ~ x1 + x2 | z1 + z2} (variable dispersion) or
#' \code{y ~ x1 + x2} (constant dispersion). The link functions for both
#' submodels are specified using \code{link.mu} and \code{link.sigma2}.
#'
#' The parametric mean link functions include a parameter \eqn{\lambda} that
#' is estimated along with other model parameters. Parameter estimation is performed
#' via maximum likelihood using the \code{optim} function with analytical gradient
#' and initial values obtained from an auxiliary linear regression of the transformed
#' response. Subsequently, the \code{optim} result may be enhanced by an additional Fisher
#' scoring iteration using analytical gradients and expected information. The Fisher
#' scoring is just a refinement to move the gradients even closer to zero and can be
#' disabled by setting \code{fsmaxit = 0} in the control arguments.
#'
#' Methods for extracting and analyzing results are implemented for objects
#' of class \code{"simplexregression"}, allowing the use of generic functions such as
#' \code{\link{summary}}, \code{\link{print}}, \code{\link{fitted}}, \code{\link{coef}},
#' \code{\link{formula}}, \code{\link{logLik}}, \code{\link{vcov}}, \code{\link{predict}},
#' \code{\link{terms}}, \code{\link{model.frame}}, \code{\link{model.matrix}},
#' \code{\link{plot}}, \code{\link{residuals}}, \code{\link{cooks.distance}},
#' \code{\link{gleverage}}, \code{\link{hatvalues}}, \code{\link{update}},
#' \code{\link{simulate}}, \code{\link{AIC}}, \code{\link[lmtest]{coeftest}} (from the
#' lmtest package), \code{\link[sandwich]{bread}} and \code{\link[sandwich]{estfun}} (from the sandwich package).
#'
#' @seealso \code{\link{summary.simplexregression}},
#' \code{\link{predict.simplexregression}}, \code{\link{residuals.simplexregression}},
#' \code{\link{penalized.ic}}, \code{\link{penalized.ss}},
#' \code{\link{scoretest}}
#'
#' @return An object of class \code{"simplexregression"}, i.e.,
#' a list with components as follows:
#' \itemize{
#'   \item \code{coefficients}: A list with elements \code{mean} and \code{dispersion},
#'   containing the estimated regression coefficients \eqn{\hat{\boldsymbol{\beta}}} and
#'   \eqn{\hat{\boldsymbol{\gamma}}} of the mean and dispersion submodels, respectively.
#'   For parametric mean link functions, the list also includes an additional element,
#'   \code{lambda}, containing the estimated link parameter \eqn{\hat{\lambda}}.
#'   \item \code{fitted.values}: a vector of fitted mean values,
#'   \item \code{optim}: a list containing \code{start} (initial values),
#'   \code{convergence} (convergence code), \code{counts} (number of iterations) and
#'   \code{method} (optimization method) from the optimization procedure,
#'   \item \code{scoring}: number of iterations from the optimization procedure via
#'   Fisher scoring,
#'   \item \code{mu.fv}: a vector of fitted mean values,
#'   \item \code{mu.lp}: a vector of fitted mean linear predictor,
#'   \item \code{mu.x}: design matrix for the mean model (with intercept),
#'   \item \code{mu.link}: character string specifying the mean link function,
#'   \item \code{mu.df}: degrees of freedom for the mean model,
#'   \item \code{sigma2.fv}: a vector of fitted dispersion values,
#'   \item \code{sigma2.lp}: a vector of fitted dispersion linear predictor,
#'   \item \code{sigma2.x}: design matrix for the dispersion model (with intercept),
#'   \item \code{sigma2.link}: character string specifying the dispersion link function,
#'   \item \code{sigma2.df}: degrees of freedom for the dispersion model,
#'   \item \code{lambda.fv}: estimated value of the parametric link function parameter
#'   (\code{NA} for mean fixed links),
#'   \item \code{df.residual}: residual degrees of freedom,
#'   \item \code{nobs}: number of observations,
#'   \item \code{loglik}: maximized log-likelihood value,
#'   \item \code{vcov}: variance-covariance matrix of the parameter estimates,
#'   \item \code{residuals}: a vector of quantile residuals,
#'   \item \code{AIC}, \code{BIC}, \code{HQIC}: Akaike, Schwarz, and Hannan-Quinn
#'   information criteria,
#'   \item \code{R2_N}, \code{R2_FC}: Nagelkerke, and Ferrari and Cribari-Neto
#'   pseudo R-squared measures,
#'   \item \code{zstat}: \eqn{z}-statistics for the coefficient tests,
#'   \item \code{pvalues}: \eqn{p}-values for the coefficient tests,
#'   \item \code{y}: the response vector,
#'   \item \code{x_names}: column names of the mean design matrix,
#'   \item \code{z_names}: column names of the dispersion design matrix,
#'   \item \code{control}: the control arguments passed to the \code{optim} call,
#'   \item \code{converged}: logical indicating successful convergence of \code{optim},
#'   \item \code{call}: the original function call,
#'   \item \code{formula}: the original two-part formula,
#'   \item \code{formula_mean}: formula for the mean submodel,
#'   \item \code{formula_disp}: formula for the dispersion submodel,
#'   \item \code{terms}: a list with \code{mean} and \code{dispersion} terms
#'   objects,
#'   \item \code{weights}: the weights used in fitting (if any),
#'   \item \code{offset}: the offset used in fitting (if any),
#'   \item \code{na.action}: the na.action attribute from the model frame,
#'   \item \code{subset}: the subset used in fitting (if any),
#'   \item \code{model}: the full model frame.
#'
#' }
#'
#' @examples
#' data(ReadingSkills, package = "SimplexRegression")
#' fit <- simplexreg(accuracy ~ dyslexia * iq | dyslexia + iq + I(iq^2),
#'                  data = ReadingSkills)
#' summary(fit)
#'
#' @references
#' Justino, M. E. C. and Cribari-Neto, F. (2026).
#' Simplex regression with a flexible logit link: Inference and application
#' to cross-country impunity data.
#' \emph{Applied Mathematical Modelling}, \bold{154}, 116713. \doi{10.1016/j.apm.2025.116713}
#'
#' Barndorff-Nielsen, O. E. and Jørgensen, B. (1991).
#' Some parametric models on the simplex.
#' \emph{Journal of Multivariate Analysis}, \bold{39}(1), 106--116.
#' \doi{10.1016/0047-259X(91)90008-P}
#'
#' Jørgensen, B. (1997). \emph{The Theory of Dispersion Models}.
#' Chapman and Hall, London.
#'
#' Song, P. X.-K. and Tan, M. (2000).
#' Marginal models for longitudinal continuous proportional data.
#' \emph{Biometrics}, \bold{56}(2), 496--502.
#' \doi{10.1111/j.0006-341X.2000.00496.x}
#'
#' Song, P. X.-K.; Qiu, Z. and Tan, M. (2004).
#' Modelling heterogeneous dispersion in marginal models for longitudinal
#' proportional data.
#' \emph{Biometrical Journal}, \bold{46}(5), 540--553.
#' \doi{10.1002/bimj.200110052}
#'
#' Song, P. X.-K. (2009).
#' Dispersion models in regression analysis.
#' \emph{Pakistan Journal of Statistics}, \bold{25}(4), 529--551.
#'
#' Zhang, P. and Qiu, Z. G. (2014).
#' Regression analysis of proportional data using simplex distribution.
#' \emph{SCIENTIA SINICA Mathematica}, \bold{44}(1), 89--104.
#' \doi{10.1360/012013-200}
NULL

#' @rdname simplexreg.fit
#' @importFrom stats model.frame model.response model.matrix terms plogis qnorm
#' @importFrom stats pnorm cor optim lm.fit weighted.mean delete.response model.weights model.offset
#' @export
simplexreg <- function(formula, data, subset, na.action, weights, offset,
                       link.mu = c("logit", "probit", "loglog", "cloglog", "cauchit",
                                   "plogit1", "plogit2"),
                       link.sigma2 = NULL, contrasts = NULL,
                       control = simplexreg.control(...),
                       model = TRUE, y = TRUE, x = FALSE, ...) {

  # Call
  cl <- match.call()

  # Set up model.frame call similar to betareg
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights",
               "offset", "contrasts"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  # Store original formula
  oformula <- as.formula(formula)

  # Handle formula - convert to formula object if needed
  if (!inherits(formula, "formula")) {
    stop("The 'formula' argument must be a formula object")
  }

  # Check for separator |
  formula_str <- as.character(formula)
  has_separator <- length(formula_str) >= 3 && grepl("\\|", formula_str[3])

  if (has_separator) {
    # Create Formula object for two-part formula
    formula <- Formula::as.Formula(formula)
    mf$formula <- formula

    # Evaluate model.frame
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    # Extract components
    Y <- model.response(mf, "numeric")

    if(any(Y <= 0 | Y >= 1, na.rm = TRUE)) {
      stop("Response variable must be in the open interval (0, 1)")
    }

    # Extract design matrices using Formula package
    X <- model.matrix(formula, data = mf, rhs = 1, contrasts.arg = contrasts)
    Z <- model.matrix(formula, data = mf, rhs = 2, contrasts.arg = contrasts)

    # Extract terms
    mtX <- terms(formula, data = mf, rhs = 1)
    mtZ <- delete.response(terms(formula, data = mf, rhs = 2))

    # Create formula objects for mean and dispersion
    formula_mean <- formula(formula, rhs = 1)
    formula_disp <- formula(formula, rhs = 2)

  } else {
    # Simple formula without separator
    mf$formula <- formula

    # Evaluate model.frame
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    # Extract components
    Y <- model.response(mf, "numeric")

    if(any(Y <= 0 | Y >= 1, na.rm = TRUE)) {
      stop("Response variable must be in the open interval (0, 1)")
    }

    # Extract design matrices
    X <- model.matrix(formula, data = mf)
    Z <- matrix(1, nrow = length(Y), ncol = 1)
    colnames(Z) <- "(Intercept)"
    rownames(Z) <- rownames(mf)

    # Extract terms
    mtX <- terms(formula, data = mf)
    mtZ <- terms(~ 1, data = mf)

    # Create formula objects
    formula_mean <- formula
    formula_disp <- as.formula("~ 1")
  }

  # Extract weights
  weights <- model.weights(mf)
  if(is.null(weights)) weights <- 1
  if(length(weights) == 1) weights <- rep.int(weights, length(Y))
  weights <- as.vector(weights)
  names(weights) <- rownames(mf)

  # Extract offset
  offset <- model.offset(mf)
  if(is.null(offset)) {
    offset <- rep(0, length(Y))
  } else {
    offset <- as.vector(offset)
  }
  names(offset) <- rownames(mf)

  # Save column names
  x_names <- colnames(X)
  z_names <- colnames(Z)

  is_constant_disp <- (ncol(Z) == 1 && all(Z[,1] == 1))
  # Set default link for sigma2 if not specified
  if (is.null(link.sigma2)) {
    if (!has_separator && is_constant_disp) {
      link.sigma2 <- "identity"
    } else {
      link.sigma2 <- "log"
    }
  } else {
    link.sigma2 <- match.arg(link.sigma2, c("log", "sqrt", "identity"))
  }

  parametric <- link.mu %in% c("plogit1","plogit2")

  # Handle null model with special initialization
  is_null_model <- (ncol(X) == 1 && all(X[,1] == 1)) && is_constant_disp
  if (is_null_model) {
    if(parametric) {
      ystar <- log(Y / (1 - Y))
      betaols_null <- weighted.mean(ystar, weights)
      muols_null <- exp(betaols_null) / (1 + exp(betaols_null))
      devtrans_null <- (Y - muols_null)^2 / (Y * (1 - Y) * muols_null^2 *
                                               (1 - muols_null)^2)
      deltaols_null <- weighted.mean(devtrans_null, weights)
      ini_nul <- c(betaols_null, deltaols_null, 1)
    } else {
      ystar <- fixed_mean_link(Y, link.mu)
      betaols_null <- weighted.mean(ystar, weights)
      muols_null <- fixed_mean_link_inv(betaols_null, link.mu)
      devtrans_null <- (Y - muols_null)^2 / (Y * (1 - Y) * muols_null^2 *
                                               (1 - muols_null)^2)
      deltaols_null <- weighted.mean(devtrans_null, weights)
      ini_nul <- c(betaols_null, deltaols_null)
    }

    rval <- simplexreg.fit(Y, X[, -1, drop = FALSE], Z[, -1, drop = FALSE],
                           weights = weights, offset = offset,
                           link.mu = link.mu, link.sigma2 = link.sigma2,
                           x_names = x_names, z_names = z_names,
                           start = ini_nul, control = control)

  } else {
    # Call standard fitting function
    rval <- simplexreg.fit(Y, X[, -1, drop = FALSE], Z[, -1, drop = FALSE],
                           weights = weights, offset = offset,
                           link.mu = link.mu, link.sigma2 = link.sigma2,
                           x_names = x_names, z_names = z_names, control = control)
  }

  # Further model information
  rval$call <- cl
  rval$formula <- oformula
  rval$formula_mean <- formula_mean
  rval$formula_disp <- formula_disp
  rval$terms <- list(mean = mtX, dispersion = mtZ)
  rval$weights <- if(identical(as.vector(weights), rep.int(1, length(Y)))) NULL else weights
  rval$offset <- if(identical(as.vector(offset), rep.int(0, length(Y)))) NULL else offset
  rval$na.action <- attr(mf, "na.action")
  rval$subset <- if(missing(subset)) NULL else subset
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- list(mean = X, dispersion = Z)

  return(rval)
}

#' @rdname simplexreg.fit
#' @importFrom stats lm.wfit sd
#' @export
simplexreg.fit <- function(y, x, z, weights = NULL, offset = NULL,
                           link.mu = c("logit", "probit", "loglog", "cloglog", "cauchit",
                                       "plogit1", "plogit2"),
                           link.sigma2 = c("log", "sqrt", "identity"),
                           x_names = NULL, z_names = NULL,
                           control = simplexreg.control(...), ...){

  y <- as.vector(y)
  n <- length(y)

  # Handle weights
  if(is.null(weights)) weights <- rep.int(1, n)
  if(length(weights) == 1) weights <- rep.int(weights, n)
  weights <- as.vector(weights)
  nobs <- sum(weights > 0)

  # Handle offset
  if(is.null(offset)) {
    offset <- rep(0, n)
  } else {
    if(length(offset) == 1) offset <- rep.int(offset, n)
    offset <- as.vector(offset)
    if(length(offset) != n) stop("offset must have length equal to number of observations")
  }

  # Ensure x is a matrix
  if (is.null(x)) {
    x <- matrix(0, nrow = n, ncol = 0)
  } else if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  # Create x1 (with intercept)
  if (ncol(x) > 0) {
    x1 <- cbind(rep(1, n), x)
    colnames(x1) <- if(!is.null(x_names)) x_names else c("(Intercept)", paste0("X", seq_len(ncol(x))))
  } else {
    x1 <- matrix(1, nrow = n, ncol = 1)
    colnames(x1) <- if(!is.null(x_names)) x_names else "(Intercept)"
  }

  # Ensure z is a matrix
  if (is.null(z)) {
    z <- matrix(0, nrow = n, ncol = 0)
  } else if (!is.matrix(z)) {
    z <- as.matrix(z)
  }

  # Create z1 (with intercept)
  if (ncol(z) > 0) {
    z1 <- cbind(1, z)
    colnames(z1) <- if(!is.null(z_names)) z_names else c("(Intercept)", paste0("Z", seq_len(ncol(z))))
  } else {
    z1 <- matrix(1, nrow = n, ncol = 1)
    colnames(z1) <- if(!is.null(z_names)) z_names else "(Intercept)"
  }

  link.mu <- match.arg(link.mu)
  link.sigma2 <- match.arg(link.sigma2)
  parametric <- link.mu %in% c("plogit1","plogit2")

  # Model dimensions
  p <- ncol(x1)  # Number of mean parameters (including intercept)
  q <- ncol(z1)  # Number of dispersion parameters (including intercept)
  r <- if(parametric) (p + q + 1) else (p + q)

  # Control parameters
  ocontrol <- control
  method <- control$method
  gradient <- control$gradient
  hessian <- control$hessian
  start <- control$start
  fsmaxit <- control$fsmaxit
  fstol <- control$fstol
  trace <- control$trace

  # Remove control parameters not needed by optim
  control_optim <- control
  control_optim$method <- control_optim$gradient <- control_optim$hessian <- NULL
  control_optim$start <- control_optim$fsmaxit <- control_optim$fstol <- NULL
  control_optim$trace <- NULL

  # Set x_names and z_names if not provided
  if(is.null(x_names)) x_names <- colnames(x1)
  if(is.null(z_names)) z_names <- colnames(z1)

  # ============================
  # FITTED VALUES FUNCTION
  # ============================
  fitfun <- function(par, deriv = 0L) {
    beta <- par[1:p]
    delta <- par[(p+1):(p+q)]

    eta1 <- as.vector(x1 %*% beta) + offset
    eta2 <- as.vector(z1 %*% delta)

    # Soft safeguards for eta2
    if (link.sigma2 == "log") {
      eta2 <- pmin(pmax(eta2, -20), 20)
    } else if (link.sigma2 == "sqrt") {
      eta2 <- pmax(eta2, 0.01)
    } else if (link.sigma2 == "identity") {
      eta2 <- pmax(eta2, 1e-6)
    }

    if(parametric) {
      lambda <- pmax(par[r], 0.001)
      mu <- as.vector(parametric_mean_link_inv(eta1, lambda, link.mu))
    } else {
      lambda <- NULL
      mu <- as.vector(fixed_mean_link_inv(eta1, link.mu))
    }

    sigma2 <- as.vector(dispersion_link_inv(eta2, link.sigma2))

    # Standard safeguards
    mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
    sigma2 <- pmax(sigma2, 1e-6)

    result <- list(
      beta = beta,
      delta = delta,
      lambda = lambda,
      eta1 = eta1,
      eta2 = eta2,
      mu = mu,
      sigma2 = sigma2
    )

    # First derivatives (for score function)
    if(deriv >= 1L) {
      if(parametric) {
        result$dmu_deta <- as.vector(parametric_mean_link_inv_deriv1(eta1, lambda, link.mu))

        # Compute rho (derivative w.r.t. lambda)
        if(link.mu == "plogit2") {
          exp_aval_frac <- plogis(eta1)^(1/lambda)
          result$rho <- as.vector(-exp_aval_frac * log(plogis(eta1)) / (lambda^2))
        } else {
          result$rho <- as.vector((-1/(lambda^2)) * ((exp(eta1) + 1)^(-1/lambda)) *
                                    log(exp(eta1) + 1))
        }
      } else {
        result$dmu_deta <- as.vector(fixed_mean_link_inv_deriv1(eta1, link.mu))
      }
      result$dsigma2_deta <- as.vector(dispersion_link_inv_deriv1(eta2, link.sigma2))
    }

    # Second derivatives (for Fisher information)
    if(deriv >= 2L) {
      result$muonemu <- as.vector(mu * (1 - mu))
      result$diff <- as.vector(y - mu)
    }

    return(result)
  }

  # ============================
  # LOG-LIKELIHOOD FUNCTION
  # ============================
  loglikfun <- function(par, fit = NULL){
    # Extract fitted parameters
    if(is.null(fit)) fit <- fitfun(par)

    with(fit, {
      # Check for problematic values
      if(!all(is.finite(mu)) || !all(is.finite(sigma2))) {
        return(-1e10)
      }

      # Safeguard y
      y_safe <- pmin(pmax(y, 1e-6), 1 - 1e-6)

      # Compute deviance
      diff <- y_safe - mu
      yoneminy <- y_safe * (1 - y_safe)
      muonemu <- mu * (1 - mu)
      dev <- (diff / muonemu)^2 / yoneminy

      if(!all(is.finite(dev))) {
        return(-1e10)
      }

      # Log-likelihood
      ll <- -0.5 * sum( weights *  (log(2 * pi) + log(sigma2) + 3 * log(y_safe * (1 - y_safe)) +
                         dev / sigma2))

      if(!is.finite(ll)) return(-1e10)

      return(ll)
    })
  }

  # ============================
  # SCORE FUNCTION (GRADIENT)
  # ============================
  gradfun <- function(par, sum = FALSE, fit = NULL){
    # Extract fitted means/precisions
    if(is.null(fit)) fit <- fitfun(par, deriv = 1L)

    with(fit, {
      # Safeguard y
      y_safe <- pmin(pmax(y, 1e-6), 1 - 1e-6)

      diff <- y_safe - mu
      yoneminy <- y_safe * (1 - y_safe)
      muonemu <- mu * (1 - mu)
      dev <- (diff / muonemu)^2 / yoneminy

      # Check for non-finite values
      if(!all(is.finite(c(mu, sigma2, diff, muonemu, dev)))) {
        warning("Non-finite values in score vector estimates.")
        if(parametric) {
          return(rep(0, p + q + 1))
        } else {
          return(rep(0, p + q))
        }
      }

      U <- weights * (dev / (sigma2 * muonemu) + 1 / (sigma2 * (muonemu)^3))
      a <- weights * ((-1 / (2 * sigma2)) + (dev / (2 * sigma2^2)))

      Ubeta <- crossprod(x1, dmu_deta * U * diff)
      Udelta <- crossprod(z1, dsigma2_deta * a)

      if(parametric) {
        if(!all(is.finite(rho))) {
          warning("Non-finite values in rho.")
          return(rep(0, p + q + 1))
        }
        Ulambda <- sum(U * rho * diff)
        rval <- c(Ubeta, Udelta, Ulambda)
      } else {
        rval <- c(Ubeta, Udelta)
      }

      return(as.vector(rval))
    })
  }

  # ============================
  # FISHER INFORMATION MATRIX
  # ============================
  hessfun <- function(par, fit) {
    with(fit, {
      # Safeguard y
      y_safe <- pmin(pmax(y, 1e-6), 1 - 1e-6)

      diff <- y_safe - mu
      yoneminy <- y_safe * (1 - y_safe)
      muonemu <- mu * (1 - mu)
      dev <- (diff / muonemu)^2 / yoneminy

      wi <- weights * ((3.0 * sigma2) / muonemu + (1 / (muonemu^3)))
      vi <- weights * (1 / (2 * sigma2^2))

      # Information matrix components
      Wbetabeta <- (1.0 / sigma2) * wi * (dmu_deta^2)
      Wdeltadelta <- vi * (dsigma2_deta^2)

      Kbetabeta <- crossprod(x1, Wbetabeta * x1)
      Kdeltadelta <- crossprod(z1, Wdeltadelta * z1)

      if(parametric) {
        Wbetalambda <- (1.0 / sigma2) * wi * dmu_deta
        Wlambdalambda <- (1.0 / sigma2) * wi

        Kbetalambda <- crossprod(x1, Wbetalambda * rho)
        Klambdabeta <- t(Kbetalambda)
        Klambdalambda <- sum(Wlambdalambda * rho^2)

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
      chol2inv(chol(K))
    })
  }

  # ============================
  # STARTING VALUES
  # ============================
  if (!is.null(start)) {
    ini <- start
  } else {
    # Standard initialization with OLS
    if(parametric){
      ystar <- log(y / (1 - y))
      auxreg <- lm.wfit(x1, ystar, weights, offset)
      betaols <- auxreg$coefficients
      muols <- exp(x1 %*% betaols) / (1 + exp(x1 %*% betaols))
      devtrans <- (y - muols)^2 / (y * (1 - y) * muols^2 * (1 - muols)^2)
      deltaols <- lm.wfit(z1, dispersion_link(devtrans, link.sigma2),
                          weights)$coefficients

      lambdaini <- 1
      ini <- c(as.vector(betaols), as.vector(deltaols), as.numeric(lambdaini))
    } else {
      ystar <- fixed_mean_link(y, link.mu)
      auxreg <- lm.wfit(x1, ystar, weights, offset)
      betaols <- auxreg$coefficients
      muols <- fixed_mean_link_inv(x1 %*% betaols, link.mu)
      devtrans <- (y - muols)^2 / (y * (1 - y) * muols^2 * (1 - muols)^2)
      deltaols <- lm.wfit(z1, dispersion_link(devtrans, link.sigma2),
                          weights)$coefficients

      ini <- c(as.vector(betaols), as.vector(deltaols))
    }
  }

  # ============================
  # OPTIMIZATION (STEP 1: BFGS)
  # ============================

  opt <- optim(par = ini,
               fn = loglikfun,
               gr = if(gradient) gradfun else NULL,
               method = method,
               hessian = hessian,
               control = control_optim)

  par <- opt$par
  ll_opt <- opt$value

  # ============================
  # FISHER SCORING (STEP 2)
  # ============================

  iter <- 0
  converged_fs <- TRUE

  if(fsmaxit > 0) {
    if(trace) cat("Starting Fisher scoring iterations...\n")

    for(iter in 1:fsmaxit) {
      # Compute step
      fit <- fitfun(par, deriv = 2L)
      scores <- gradfun(par, fit = fit)
      InfoInv <- try(hessfun(par, fit = fit), silent = TRUE)

      if(inherits(InfoInv, "try-error")) {
        warning("Failed to invert information matrix")
        converged_fs <- FALSE
        break
      }

      step <- InfoInv %*% scores
      par_new <- par + step

      # Check improvement with backtracking
      step_factor <- 0
      improved <- FALSE

      while(step_factor <= 10 && !improved) {
        ll_new <- loglikfun(par_new)

        if(ll_new > ll_opt) {
          par <- par_new
          ll_opt <- ll_new
          improved <- TRUE
        } else {
          # Reduce step size
          par_new <- par + (0.5^step_factor) * step
          step_factor <- step_factor + 1
        }
      }

      # Convergence check
      if(all(abs(step) < fstol)) {
        if(trace) cat(sprintf("Converged after %d Fisher scoring iterations\n", iter))
        break
      }

      if(!improved && step_factor > 10) {
        if(trace) cat("No improvement after backtracking\n")
        break
      }
    }

    if(iter >= fsmaxit) {
      converged_fs <- FALSE
      warning(sprintf("Fisher scoring did not converge in %d iterations", fsmaxit))
    }
  }

  # ============================
  # EXTRACT FINAL ESTIMATES
  # ============================
  fit <- fitfun(par, deriv = 2L)
  beta <- fit$beta
  delta <- fit$delta
  lambda <- fit$lambda
  mu <- fit$mu
  sigma2 <- fit$sigma2
  eta1 <- fit$eta1
  eta2 <- fit$eta2

  # Final log-likelihood
  ll <- loglikfun(par, fit = fit)

  # Variance-covariance matrix
  vcov <- hessfun(par, fit = fit)

  # ============================
  # NAMING AND DIAGNOSTICS
  # ============================
  names(beta) <- x_names
  names(delta) <- z_names

  if(parametric) {
    vcov_names <- c(
      x_names,
      if(q == 1L && z_names[1] == "(Intercept)") {
        "(dispersion)"
      } else {
        paste("(dispersion)", z_names, sep = "_")
      },
      "lambda"
    )
  } else {
    vcov_names <- c(
      x_names,
      if(q == 1L && z_names[1] == "(Intercept)") {
        "(dispersion)"
      } else {
        paste("(dispersion)", z_names, sep = "_")
      }
    )
  }

  rownames(vcov) <- colnames(vcov) <- vcov_names
  stderror <- sqrt(diag(vcov))

  # Test statistics
  zstat <- abs(par / stderror)
  pvalues <- 2 * (1 - pnorm(zstat))

  # Residuals
  quantile_res <- qnorm(psimplex(y, mu, sigma2))

  # Pseudo R-squared
  R2_N <- 1 - exp((-2/sum(weights)) * (ll - simplexreg.nul(y, link.mu, weights)))

  gy <- if(parametric) {
    parametric_mean_link(y, lambda, link.mu)
  } else {
    fixed_mean_link(y, link.mu)
  }

  if(sd(eta1) > 0 && sd(gy) > 0) {
    R2_FC <- cor(eta1, gy)^2
  } else {
    R2_FC <- NA
  }

  # Information criteria
  aic <- -2 * ll + 2 * r
  bic <- -2 * ll + log(n) * r
  hqic <- -2 * ll + 2 * r * log(log(n))

  # ============================
  # RETURN OBJECT
  # ============================
  mean_coefs <- beta
  names(mean_coefs) <- x_names

  disp_coefs <- delta
  names(disp_coefs) <- z_names

  # Construct result object
  result <- list(
    coefficients = list(
      mean = mean_coefs,
      dispersion = disp_coefs,
      lambda = if(parametric) lambda else NA
    ),
    fitted.values = structure(mu, .Names = seq_len(n)),
    optim = list(
      start = ini,
      convergence = opt$convergence,
      counts = opt$counts, # mudei aqui
      method = method
    ),
    scoring = iter,
    mu.fv = structure(mu, .Names = seq_len(n)),
    mu.lp = structure(eta1, .Names = seq_len(n)),
    mu.x = x1,
    mu.link = link.mu,
    mu.df = p,
    sigma2.fv = structure(sigma2, .Names = seq_len(n)),
    sigma2.lp = structure(eta2, .Names = seq_len(n)),
    sigma2.x = z1,
    sigma2.link = link.sigma2,
    sigma2.df = q,
    lambda.fv = lambda,
    df.residual = n - r,
    nobs = n,
    loglik = ll,
    vcov = vcov,
    residuals = structure(quantile_res, .Names = seq_len(n)),
    AIC = aic,
    BIC = bic,
    HQIC = hqic,
    R2_FC = R2_FC,
    R2_N = R2_N,
    zstat = zstat,
    pvalues = pvalues,
    y = structure(y, .Names = seq_len(n)),
    x_names = x_names,
    z_names = z_names,
    control = ocontrol,
    converged = (opt$convergence == 0) && converged_fs
  )

  class(result) <- "simplexregression"

  return(result)
}

# ==============================================================================
# 3. NULL MODEL LOG-LIKELIHOOD
# ==============================================================================

#' @title Null Model Log-Likelihood for Simplex Regression
#' @description Computes the log-likelihood for the null (intercept-only) model in
#' simplex regression with a parametric or fixed mean link.
#'
#' @param y Numeric response vector \eqn{(0 < y < 1)}.
#' @param link.mu Mean link function: parametric ("plogit1", "plogit2") or fixed
#'   ("logit", "probit", "loglog", "cloglog", "cauchit").
#' @param weights Optional vector of weights (default: \code{NULL}).
#'
#' @return Numeric value of the null model log-likelihood
#' @keywords internal
simplexreg.nul <- function(y, link.mu, weights = NULL){
  y <- as.vector(y)
  n <- length(y)

  parametric <- link.mu %in% c("plogit1","plogit2")

  # Null model log-likelihood function
  fLogLiknull <- function(theta){
    sigma2 <- theta[2]
    eta1 <- cbind(rep(1, n)) %*% theta[1]

    if(parametric){
      lambda = pmax(theta[3], 0.001)
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
    adFunc <- -0.5 * sum( weights * (log(2*pi) + log(sigma2) + 3*log(yoneminy) + dev/sigma2))
    adFunc
  }

  # Initial values
  if(parametric){
    ystar <- log(y / (1 - y))
    betaols_null <- weighted.mean(ystar, weights)
    muols_null <- exp(betaols_null) / (1 + exp(betaols_null))
    devtrans_null <- (y - muols_null)^2 / (y * (1 - y) * muols_null^2 *
                                             (1 - muols_null)^2)
    deltaols_null <- weighted.mean(devtrans_null, weights)

    ini_nul <- c(betaols_null, deltaols_null, 1)
  } else {
    ystar <- fixed_mean_link(y, link.mu)
    betaols_null <- weighted.mean(ystar, weights)
    muols_null <- fixed_mean_link_inv(betaols_null, link.mu)
    devtrans_null <- (y - muols_null)^2 / (y * (1 - y) * muols_null^2 *
                                             (1 - muols_null)^2)
    deltaols_null <- weighted.mean(devtrans_null, weights)

    ini_nul <- c(betaols_null, deltaols_null)
  }

  # Optimize
  opt_nul <- optim(ini_nul, fLogLiknull, method = "BFGS",
                   control=list(fnscale = -1, maxit = 5000, reltol = .Machine$double.eps^(0.5)))

  k <- opt_nul$value
  return(k)
}

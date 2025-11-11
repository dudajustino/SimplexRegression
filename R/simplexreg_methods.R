################################################################################
#                     SIMPLEX REGRESSION - S3 METHODS                          #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2025-11-08                                                             #
# Description: S3 methods for simplexregression objects including summary,     #
#              print, predict, coef, vcov, logLik, AIC, BIC, and others        #
################################################################################

# ==============================================================================
# HELPER: Check if model uses parametric mean link function
# ==============================================================================
is_parametric <- function(object) {
  !is.na(object$coefficients$lambda)
}

# ==============================================================================
# 1. SUMMARY METHOD
# ==============================================================================

#' @title Summary Method for Simplex Regression
#' @description Produces a summary of a fitted simplex regression model.
#'
#' @param object An object of class \code{"simplexregression"}
#' @param ... Additional arguments (currently not used)
#'
#' @return An object of class \code{"summary.simplexregression"} containing coefficient tables and diagnostics
#'
#' @export
summary.simplexregression <- function(object, ...) {
  #  Extract coefficients
  coef_mean <- object$coefficients$mean
  coef_disp <- object$coefficients$dispersion
  coef_lambda <- object$coefficients$lambda

  parametric <- is_parametric(object)

  # Get standard errors from variance-covariance matrix
  vcov_matrix <- object$vcov
  se <- sqrt(diag(vcov_matrix))

  # Dimensions
  p <- length(coef_mean)
  q <- length(coef_disp)

  # Standard errors
  se_mean <- se[1:p]
  se_disp <- se[(p+1):(p+q)]

  # lambda statistics (only for models with parametric mean link)
  if (parametric) {
    se_lambda <- se[p+q+1]
    zstat_lambda <- coef_lambda / se_lambda
    pval_lambda <- 2 * pnorm(-abs(zstat_lambda))
  } else {
    se_lambda <- NA
    zstat_lambda <- NA
    pval_lambda <- NA
  }

  # Z-statistics and p-values
  zstat_mean <- coef_mean / se_mean
  pval_mean <- 2 * pnorm(-abs(zstat_mean))

  zstat_disp <- coef_disp / se_disp
  pval_disp <- 2 * pnorm(-abs(zstat_disp))

  # Use regressor names if available
  if (!is.null(object$x_names)) {
    names(coef_mean) <- object$x_names
  }

  if (!is.null(object$z_names)) {
    names(coef_disp) <- object$z_names
  }

  # Create coefficient tables
  coef_table_mean <- cbind(
    Estimate = coef_mean,
    `Std. Error` = se_mean,
    `z value` = zstat_mean,
    `Pr(>|z|)` = pval_mean
  )

  coef_table_disp <- cbind(
    Estimate = coef_disp,
    `Std. Error` = se_disp,
    `z value` = zstat_disp,
    `Pr(>|z|)` = pval_disp
  )

  # Lambda parameter (if parametric)
  if (parametric) {
    lambda_table <- cbind(
      Estimate = coef_lambda,
      `Std. Error` = se_lambda,
      `z value` = zstat_lambda,
      `Pr(>|z|)` = pval_lambda
    )
    rownames(lambda_table) <- "lambda"
  } else {
    lambda_table <- NULL
  }

  result <- list(
    call = object$call,
    coef_mean = coef_table_mean,
    coef_disp = coef_table_disp,
    lambda = lambda_table,
    parametric = parametric,
    mu.link = object$mu.link,
    sigma2.link = object$sigma2.link,
    loglik = object$loglik,
    aic = object$AIC,
    bic = object$BIC,
    hqic = object$HQ,
    nobs = object$nobs,
    df.residual = object$df.residual,
    counts = object$optim$counts,
    R2_RV = object$R2_RV,
    R2_FC = object$R2_FC
  )

  class(result) <- "summary.simplexregression"
  return(result)
}

#' @title Print Method for Summary of Simplex Regression
#' @description Prints a summary of a fitted simplex regression model.
#'
#' @param x An object of class \code{"summary.simplexregression"}
#' @param digits Number of digits to display (default: 3)
#' @param ... Additional arguments passed to \code{printCoefmat}
#' @importFrom stats printCoefmat
#' @export
print.summary.simplexregression <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:\n")
  print(x$call)

  parametric <- x$parametric

  cat("\n")
  if (parametric) {
    cat("Simplex Regression with Parametric Mean Link:", x$mu.link, "\n")
  } else {
    cat("Simplex Regression with Fixed Mean Link:", x$mu.link, "\n")
  }
  cat("Dispersion Link:", x$sigma2.link, "\n")

  cat("\nMean model coefficients:\n")
  printCoefmat(round(x$coef_mean,4), P.values = TRUE, has.Pvalue = TRUE)

  cat("\nDispersion model coefficients:\n")
  printCoefmat(round(x$coef_disp,4), P.values = TRUE, has.Pvalue = TRUE)

  if (parametric) {
    cat("\nLink function parameter:\n")
    printCoefmat(x$lambda, P.values = TRUE, has.Pvalue = TRUE, digits = digits)
  }

  cat("\nLog-likelihood:", round(x$loglik, digits), "\n")
  cat("AIC:", round(x$aic, digits), "\n")
  cat("BIC:", round(x$bic, digits), "\n")
  cat("HQIC:", round(x$hqic, digits), "\n")
  cat("Number of observations:", x$nobs, "\n")
  cat("Degrees of freedom:", x$df.residual, "\n")
  cat("Number of iterations in BFGS optim:", x$counts, "\n")
  cat("Pseudo-R2 Nagelkerke:", round(x$R2_RV, digits), "\n")
  cat("Pseudo-R2 Ferrari e Cribari-Neto:", round(x$R2_FC, digits), "\n")

  invisible(x)
}

# ==============================================================================
# 2.AUXILIARY METHODS: model.matrix, terms, print
# ==============================================================================

#' @title Extract Fitted Values
#' @description Extracts fitted values from a simplex regression model.
#'
#' @param object An object of class \code{"simplexregression"}
#' @param ... Additional arguments (currently not used)
#'
#' @return Numeric vector of fitted mean values
#' @export
fitted.simplexregression <- function(object, ...) {
  object$fitted.values
}

#' @title Extract Model Formula
#' @description Extracts the model formula from a simplex regression model.
#'
#' @param x An object of class \code{"simplexregression"}
#' @param ... Additional arguments (currently not used)
#' @importFrom stats as.formula
#' @return Formula object
#' @export
formula.simplexregression <- function(x, ...) {
  x$formula
}

#' @title Extract Coefficients
#' @description Extracts all coefficients from a simplex regression model.
#'
#' @param object An object of class \code{"simplexregression"}
#' @param ... Additional arguments (currently not used)
#' @importFrom stats coef
#' @return Named numeric vector with all coefficients (mean, dispersion, lambda)
#' or (mean, dispersion)
#' @export
coef.simplexregression <- function(object, ...) {
  if (is_parametric(object)) {
    c(
      mean = object$coefficients$mean,
      dispersion = object$coefficients$dispersion,
      lambda = object$coefficients$lambda
    )
  } else {
    c(
      mean = object$coefficients$mean,
      dispersion = object$coefficients$dispersion
    )
  }
}

#' @title Extract Log-Likelihood
#' @description Extracts the log-likelihood from a fitted model.
#'
#' @param object An object of class \code{"simplexregression"}
#' @param ... Additional arguments (currently not used)
#'
#' @return Object of class \code{"logLik"} with attributes df and nobs
#' @export
logLik.simplexregression <- function(object, ...) {
  val <- object$loglik
  attr(val, "df") <- length(object$coefficients$mean) +
    length(object$coefficients$dispersion) +
    ifelse(is_parametric(object), 1, 0)
  attr(val, "nobs") <- object$nobs
  class(val) <- "logLik"
  val
}

#' @title Extract Variance-Covariance Matrix
#' @description Extracts the variance-covariance matrix of parameter estimates.
#'
#' @param object An object of class \code{"simplexregression"}
#' @param ... Additional arguments (currently not used)
#'
#' @return Variance-covariance matrix
#' @export
vcov.simplexregression <- function(object, ...) {
  object$vcov
}

#' @title Extract Number of Observations
#' @description Extracts the number of observations used in model fitting.
#'
#' @param object An object of class \code{"simplexregression"}
#' @param ... Additional arguments (currently not used)
#' @importFrom stats nobs
#' @return Integer number of observations
#' @export
nobs.simplexregression <- function(object, ...) {
  object$nobs
}

#' @title Extract Residual Degrees of Freedom
#' @description Extracts the residual degrees of freedom.
#'
#' @param object An object of class \code{"simplexregression"}
#' @param ... Additional arguments (currently not used)
#'
#' @return Integer residual degrees of freedom
#' @export
df.residual.simplexregression <- function(object, ...) {
  object$nobs - (length(object$coefficients$mean) +
                   length(object$coefficients$dispersion) +
                   ifelse(is_parametric(object), 1, 0))
}

#' @title Model Matrix Extraction
#' @description Extracts the design matrix used in the mean or dispersion submodel.
#'
#' @param object An object of class \code{"simplexregression"}
#' @param type Character string indicating which design matrix to extract:
#'   \code{"mean"} (default) or \code{"dispersion"}
#' @param ... Additional arguments (currently not used)
#'
#' @return Model matrix (numeric matrix)
#' @export
model.matrix.simplexregression <- function(object, type = c("mean", "dispersion"), ...) {
  type <- match.arg(type)
  if (type == "mean") return(object$mu.x)
  else return(object$sigma2.x)
}

#' @title Terms Extraction
#' @description Extracts the terms object for the mean or dispersion submodel.
#'
#' @param x An object of class \code{"simplexregression"}
#' @param type Character string indicating which terms to extract:
#'   \code{"mean"} (default) or \code{"dispersion"}
#' @param ... Additional arguments (currently not used)
#'
#' @return An object of class \code{"terms"}
#' @export
terms.simplexregression <- function(x, type = c("mean", "dispersion"), ...) {
  type <- match.arg(type)

  if (type == "mean") return(x$terms$mean)
  else return(x$terms$dispersion)
}

#' @title Print Method for Simplex Regression
#' @description Prints a concise summary of a fitted simplex regression model.
#'
#' @param x An object of class \code{"simplexregression"}
#' @param digits Number of digits to print
#' @param ... Additional arguments (currently not used)
#' @export
print.simplexregression <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:\n")
  print(x$call)

  if (is_parametric(x)) {
    cat("\nSimplex Regression with Parametric Mean Link:", x$mu.link, "\n")
  } else {
    cat("\nSimplex Regression with Fixed Mean Link:", x$mu.link, "\n")
  }
  cat("Dispersion Link:", x$sigma2.link, "\n")

  cat("\nCoefficients (Mean model):\n")
  print(round(x$coefficients$mean, digits))

  cat("\nCoefficients (Dispersion model):\n")
  print(round(x$coefficients$dispersion, digits))

  if (is_parametric(x)) {
    cat("\nLambda:", round(x$coefficients$lambda, digits), "\n")
  }

  invisible(x)
}

# ==============================================================================
# 3. INFORMATION CRITERIA
# ==============================================================================

#' @title Akaike Information Criterion
#' @description Computes the AIC for a simplex regression model.
#'
#' @param object An object of class \code{"simplexregression"}
#' @param ... Additional arguments (currently not used)
#' @param k Penalty parameter (default: 2)
#' @importFrom stats AIC
#' @return Numeric AIC value
#' @export
AIC.simplexregression <- function(object, ..., k = 2) {
  -2 * object$loglik + k * (length(object$coefficients$mean) +
                              length(object$coefficients$dispersion) +
                              ifelse(is_parametric(object), 1, 0))
}

#' @title Bayesian Information Criterion
#' @description Computes the BIC for a simplex regression model.
#'
#' @param object An object of class \code{"simplexregression"}
#' @param ... Additional arguments (currently not used)
#' @importFrom stats BIC
#' @return Numeric BIC value
#' @export
BIC.simplexregression <- function(object, ...) {
  -2 * object$loglik + log(object$nobs) * (length(object$coefficients$mean) +
                                             length(object$coefficients$dispersion) +
                                             ifelse(is_parametric(object), 1, 0))
}

# ==============================================================================
# 4. PREDICTION METHOD
# ==============================================================================

#' @title Predict Method for Simplex Regression
#' @description Computes predictions from a fitted simplex regression model.
#'
#' @param object An object of class \code{"simplexregression"}
#' @param newdata Optional data frame with new predictor values. If \code{NULL},
#'   predictions are made for the original data
#' @param type Character string specifying the type of prediction:
#' \itemize{
#'   \item \code{"response"}: Predicted mean values (default)
#'   \item \code{"link"}: Linear predictors for mean and dispersion
#'   \item \code{"dispersion"}: Predicted dispersion values
#' }
#' @param ... Additional arguments (currently not used)
#'
#' @return Numeric vector or list of predictions depending on \code{type}
#'
#' @export
predict.simplexregression <- function(object, newdata = NULL, type = c("link", "response", "dispersion"), ...) {
  type <- match.arg(type)

  parametric <- is_parametric(object)

  # If no new data, use fitted values
  if (is.null(newdata)) {
    if (type == "link") {
      return(list(mean = object$mu.lp, dispersion = object$sigma2.lp))
    } else if (type == "response") {
      return(object$fitted.values)
    } else if (type == "dispersion") {
      return(object$sigma2.fv)
    }
  } else {
    # Predictions for new data
    if (type == "response") {
      # Extract design matrix for mean
      x_matrix_new <- model.matrix(object$formula_mean, data = newdata)
      beta <- object$coefficients$mean

      # Check dimensions
      if (ncol(x_matrix_new) != length(beta)) {
        stop("Incompatible dimensions between newdata and fitted model")
      }

      eta1 <- as.vector(x_matrix_new %*% beta)

      # Use appropriate link function
      if (parametric) {
        mu_pred <- parametric_mean_link_inv(eta1, object$lambda.fv, object$mu.link)
      } else {
        mu_pred <- fixed_mean_link_inv(eta1, object$mu.link)
      }

      return(mu_pred)

    } else if (type == "dispersion") {
      # Extract design matrix for dispersion
      z_matrix_new <- model.matrix(object$formula_disp, data = newdata)
      delta <- object$coefficients$dispersion

      # Check dimensions
      if (ncol(z_matrix_new) != length(delta)) {
        stop("Incompatible dimensions between newdata and fitted model")
      }

      eta2 <- as.vector(z_matrix_new %*% delta)
      sigma2_pred <- dispersion_link_inv(eta2, object$sigma2.link)
      return(sigma2_pred)

    } else if (type == "link") {
      # Return linear predictors for both
      x_matrix_new <- model.matrix(object$formula_mean, data = newdata)
      z_matrix_new <- model.matrix(object$formula_disp, data = newdata)

      beta <- object$coefficients$mean
      delta <- object$coefficients$dispersion

      eta1 <- as.vector(x_matrix_new %*% beta)
      eta2 <- as.vector(z_matrix_new %*% delta)

      return(list(mean = eta1, dispersion = eta2))
    }
  }
}

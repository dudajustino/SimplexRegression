################################################################################
#               SIMPLEX REGRESSION - PENALIZED INFORMATION CRITERIA            #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2026-05-26                                                             #
# Description: Penalized information criterion for simplex regression model    #
#              selection                                                       #
################################################################################

#' @title Penalized Information Criteria for Simplex Regression Model Selection
#' @description Implements the Akaike, Schwarz, and Hannan–Quinn information
#' criteria with a penalty term for selecting among competing simplex regression
#' models with parametric mean link functions.
#'
#' @param ... One or more objects of class \code{simplexregression} fitted with a
#' parametric mean link functions ("plogit1" or "plogit2").
#' @param kappa A numeric value controlling the additional penalty for the link
#' mean parameter. Default is 0.1.
#' @param verbose Logical. If \code{TRUE} (default), prints the criteria values.
#' If \code{FALSE}, returns results silently.
#' @param digits Integer specifying the number of decimal places for output.
#' Default is \code{max(3, getOption("digits") - 3)}.
#'
#' @details
#' The penalized information criteria, as proposed by Justino and
#' Cribari-Neto (2026), extend the classical Akaike, Schwarz and
#' Hannan--Quinn criteria with an additional penalty term for the link
#' parameter \eqn{\lambda}:
#' \deqn{AIC^{(\lambda)} = -2 \ell + (2 + \kappa \, |\log(\lambda)|)r}
#' \deqn{BIC^{(\lambda)} = -2 \ell + (\log(n) + \kappa \, |\log(\lambda)|)r}
#' \deqn{HQIC^{(\lambda)}  = -2 \ell + (2 \log(\log(n)) + \kappa \, |\log(\lambda)|)r}
#' where:
#' \itemize{
#'   \item \eqn{\ell} denotes the maximized log-likelihood function;
#'   \item \eqn{\kappa \geq 0} controls the additional penalty associated with
#'   the link parameter;
#'   \item \eqn{\lambda} is the parameter of the parametric mean link function;
#'   \item \eqn{r} indicate the dimension of their parameter vector;
#'   \item \eqn{n} is the number of observations.
#' }
#'
#' \strong{Important}: These penalized versions of the criteria should only be used
#' when the candidate model employ a parametric link function in the mean submodel
#' (use \code{kappa = 0.1}). When candidate models includes specifications
#' with fixed link functions, the standard unpenalized versions of these criteria
#' should be applied instead (use \code{kappa = 0}).
#'
#' @return A data frame with rows named after the candidate models and four columns:
#' \describe{
#'   \item{\code{df}}{Number of estimated parameters.}
#'   \item{\code{AICc}}{Penalized AIC value.}
#'   \item{\code{BICc}}{Penalized BIC value.}
#'   \item{\code{HQICc}}{Penalized HQIC value.}
#' }
#' When \code{verbose = TRUE}, the results are also printed to the console and
#' the data frame is returned invisibly. When \code{verbose = FALSE}, the data
#' frame is returned visibly without printing.
#'
#' @examples
#' # Simulate data
#' set.seed(2026)
#' n <- 100
#' x1 <- runif(n, 0, 1)
#' x2 <- runif(n, 0, 1)
#' z1 <- runif(n, 0, 1)
#' mu <- parametric_mean_link_inv(0.6 - 2*x1 - 1.5*x2, 0.5, "plogit1")
#' sigma2 <- dispersion_link_inv(-2 - 2.5*z1, "log")
#' y <- rsimplex(n, mu, sigma2)
#' data <- data.frame(y = y, x1 = x1, x2 = x2, z1 = z1)
#'
#' # Fit two models with parametric mean link functions
#' fit1 <- simplexreg(y ~ x1 + x2 | z1, data = data, link.mu = "plogit1")
#' fit2 <- simplexreg(y ~ x1 + x2 | z1, data = data, link.mu = "plogit2")
#'
#' # Compute penalized criteria
#' penalized.ic(fit1, fit2)
#'
#' @references
#' Justino, M. E. C. and Cribari-Neto, F. (2026).
#' Simplex regression with a flexible logit link: Inference and application
#' to cross-country impunity data.
#' \emph{Applied Mathematical Modelling}, \bold{154}, 116713. \doi{10.1016/j.apm.2025.116713}
#'
#' Akaike, H. (1973). Information theory and an extension of the maximum
#' likelihood principle. \emph{Akadémiai Kiadó}, 267--281.
#'
#' Schwarz, G. E. (1978). Estimating the dimension of a model. \emph{Annals of
#' Statistics}, \bold{6}(2), 461-–464. \doi{10.1214/aos/1176344136}
#'
#' Hannan, E. J. and Quinn, B. G. (1979). The Determination of the Order of an
#' Autoregression. \emph{Journal of the Royal Statistical Society Series B:
#' Statistical Methodology}, \bold{41}(2), 190--195.
#' \doi{10.1111/j.2517-6161.1979.tb01072.x}
#'
#' @seealso \code{\link{penalized.ss}}
#' @export
penalized.ic <- function(..., kappa = 0.1, verbose = TRUE,
                         digits = max(3, getOption("digits") - 3)) {

  models <- list(...)
  M <- length(models)

  # Verify all objects are simplexregression models
  if (!all(sapply(models, function(x) inherits(x, "simplexregression")))) {
    stop("All arguments must be objects of class 'simplexregression'.")
  }

  # Check parametric mean link
  has_lambda <- sapply(models, function(x) !is.null(x$lambda.fv))

  call <- match.call()
  user_specified_kappa <- "kappa" %in% names(call)

  # Prevent misuse of kappa
  if (!all(has_lambda)) {
    if (user_specified_kappa && kappa > 0) {
      warning("kappa > 0 is only valid when all models have parametric mean links. ",
              "Setting kappa = 0 automatically for this comparison.")
      kappa <- 0
    } else if (!user_specified_kappa && kappa > 0) {
      # User didn't specify kappa, using default but models don't support it
      kappa <- 0
      # No warning here - it's expected behavior
    }
  }

  # Get model names
  model_names <- names(models)
  if (is.null(model_names)) {
    model_names <- vapply(as.list(substitute(list(...)))[-1L], deparse, character(1))
  }

  # Initialize storage
  df_values <- aic_values <- bic_values <- hqic_values <- numeric(M)

  for (i in seq_len(M)) {
    model <- models[[i]]

    # Calculate degrees of freedom (number of parameters)
    df_values[i] <- length(model$coefficients$mean) +
      length(model$coefficients$dispersion) +
      as.integer(!is.null(model$lambda.fv))

    n <- length(model$fitted.values)
    ll <- model$loglik

    if (has_lambda[i]) {
      lambda <- model$lambda.fv
      penalty <- kappa * abs(log(lambda))
    } else {
      penalty <- 0
    }

    aic_values[i] <- -2 * ll + (2 + penalty) * df_values[i]
    bic_values[i] <- -2 * ll + (log(n) + penalty) * df_values[i]
    hqic_values[i] <- -2 * ll + (2 * log(log(n)) + penalty) * df_values[i]
  }

  result <- data.frame(
    df = df_values,
    AICc = aic_values,
    BICc = bic_values,
    HQICc = hqic_values,
    row.names = model_names
  )

  if (kappa == 0) {
    names(result)[-1] <- c("AIC", "BIC", "HQIC")
  }

  # Print message
  if (verbose) {
    cat("\n")
    if (kappa == 0) {
      cat("Information criteria values:\n")
    } else {
      cat(sprintf("Penalized information criteria values (kappa = %s):\n",
                  formatC(kappa, digits = 3, format = "f")))
    }

    result_print <- data.frame(
      df = formatC(result[, 1], digits = 0, format = "f"),
      AICc = formatC(result[, 2], digits = digits, format = "f"),
      BICc = formatC(result[, 3], digits = digits, format = "f"),
      HQICc = formatC(result[, 4], digits = digits, format = "f"),
      row.names = rownames(result)
    )

    names(result_print) <- names(result)

    print(result_print, quote = FALSE, right = TRUE)
    return(invisible(result))
  } else {
    return(result)
  }
}

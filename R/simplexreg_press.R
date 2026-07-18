################################################################################
#               SIMPLEX REGRESSION - PRESS AND P2                              #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2026-05-26                                                             #
# Description: PRESS-Based P2 Statistics for Simplex Regression                #
################################################################################

#' @title PRESS-Based \eqn{P2} Statistics for Simplex Regression
#' @description Computes the PRESS (Predicted Residual Error Sum of Squares)
#' statistic and the associated \eqn{P^2} and adjusted \eqn{P^2} measures for simplex
#' regression models with parametric or fixed mean link function, as proposed by
#' Espinheira and Silva (2026).
#'
#' @param ... One or more objects of class \code{"simplexregression"}.
#' @param type Character string specifying the type of residual to use.
#' Options are \code{"standardized"} (default) or \code{"biasvariance"}
#' (see \code{\link{residuals.simplexregression}}).
#'
#' @details
#' The PRESS statistic for the simplex regression model is given by:
#' \deqn{\text{PRESS} = \sum_{i=1}^{n} \left(\frac{r_i}{1 - h_{ii}}\right)^2,}
#' where \eqn{r_i} denotes the residual for observation \eqn{i} and
#' \eqn{\hat{h}_{ii}} is the \eqn{i}-th diagonal element of the hat matrix
#' (see \code{\link{hatvalues.simplexregression}}).
#'
#' The \eqn{P^2} statistic is a cross-validation analog of \eqn{R^2}, defined
#' for the simplex regression mode as:
#' \deqn{P^2 = 1 - \frac{\text{PRESS}}{\left(\frac{n}{n-r}\right)^2 \text{SST}},}
#' where \eqn{r} is the number of estimated parameters,
#' \eqn{\text{SST} = \sum_{i=1}^{n}(\check{y}_i - \bar{\check{\boldsymbol{y}}})^2},
#' with \eqn{\bar{\check{y}}} is the mean of the transformed fitted values
#' \eqn{\check{y}_i} defined by Espinheira and Silva (2026).
#'
#' The adjusted \eqn{P^2} is given by:
#' \deqn{P^2_c = 1 - (1 - P^2)\frac{n-1}{n-r}.}
#'
#' The \code{type} argument controls which residuals \eqn{r_i} are used in
#' the PRESS computation. Only \code{"standardized"} and \code{"biasvariance"}
#' residuals are supported, as these are the residual types for which the PRESS-based
#' cross-validation analog is defined in Espinheira and Silva (2026)
#'
#' Values of \eqn{P^2} and \eqn{P^2_c} closer to 1 indicate better predictive
#' performance of the model. Since PRESS is nonnegative, both measures are bounded
#' above by 1. Hence,
#' \deqn{P^2, P^2_c \in (-\infty, 1].}
#'
#' @return When a single model is provided, a named numeric vector with components
#' \code{P2}, \code{P2_c}, and \code{PRESS}. When multiple models are provided,
#' a data frame with one row per model and columns \code{P2}, \code{P2_c}, and
#' \code{PRESS}.
#'
#' @references
#' Espinheira, P. L. and Silva, A. O. (2020). Residual and influence analysis to a
#' general class of simplex regression. \emph{TEST}, \bold{29}, 523--552.
#' \doi{10.1007/s11749-019-00665-3}
#'
#' Espinheira, P. L. and Silva, A. O. (2026). Prediction in the nonlinear simplex model.
#' \emph{International Journal of Data Science and Analytics}, \bold{22}, 161.
#' \doi{10.1007/s41060-026-01114-9}
#'
#' @examples
#' data(ReadingSkills, package = "SimplexRegression")
#' fit <- simplexreg(accuracy ~ dyslexia * iq | dyslexia + iq + I(iq^2),
#'                  data = ReadingSkills)
#'
#' fit1 <- simplexreg(accuracy ~ dyslexia * iq | dyslexia + iq + I(iq^2),
#'                  data = ReadingSkills, link.mu = "loglog")
#' # Single model
#' press(fit)
#'
#' # Comparing multiple models
#' press(fit, fit1)
#'
#' # Using bias-variance residuals
#' press(fit, fit1, type = "biasvariance")
#'
#' @seealso \code{\link{simplexreg}}, \code{\link{residuals.simplexregression}}.
#' @export
press <- function(..., type = c("standardized", "biasvariance")) {

  models <- list(...)
  M <- length(models)

  # Verify all objects are simplexregression models
  if (!all(sapply(models, function(x) inherits(x, "simplexregression")))) {
    stop("All arguments must be objects of class 'simplexregression'.")
  }

  type <- match.arg(type)

  # Get model names (supports named arguments)
  model_names <- names(models)
  if (is.null(model_names)) {
    model_names <- vapply(as.list(substitute(list(...)))[-1L], deparse, character(1))
  }

  p2_values <- p2_c_values <- press_values <- numeric(M)

  for (i in seq_len(M)) {
    model <- models[[i]]
    parametric <- !is.na(model$coefficients$lambda)

    res    <- residuals(model, type)
    n      <- model$nobs
    p      <- ncol(model$mu.x)
    q      <- ncol(model$sigma2.x)
    r      <- if (parametric) p + q + 1 else p + q

    eta1   <- model$mu.lp
    y      <- as.vector(model$y)
    mu     <- as.vector(model$mu.fv)
    sigma2 <- as.vector(model$sigma2.fv)

    diff     <- y - mu
    muonemu  <- mu * (1 - mu)
    yoneminy <- y * (1 - y)
    dev      <- (diff / muonemu)^2 / yoneminy

    if (parametric) {
      Dlink <- parametric_mean_link_inv_deriv1(eta1, model$lambda.fv, model$mu.link)
    } else {
      Dlink <- fixed_mean_link_inv_deriv1(eta1, model$mu.link)
    }

    Wi    <- (3 * sigma2 / muonemu) + (1 / (muonemu^3))
    Ui    <- (dev / muonemu) + (1 / (muonemu^3))
    ychap <- sqrt((1 / sigma2) * Wi * Dlink^2) * (eta1 + (Ui * diff) / (Dlink * Wi))

    hii <- hatvalues(model)
    press_values[i] <- sum((res / (1 - hii))^2)

    ychap_mean      <- mean(ychap)
    sst             <- sum((ychap - ychap_mean)^2)
    n_factor        <- (n / (n - r))^2
    p2_values[i]    <- 1 - (press_values[i] / (n_factor * sst))
    p2_c_values[i]  <- 1 - (1 - p2_values[i]) * ((n - 1) / (n - r))
  }

  if (M == 1) {
    return(c(P2 = p2_values, P2_c = p2_c_values, PRESS = press_values))
  } else {
    return(data.frame(P2 = p2_values,
                      P2_c =  p2_c_values,
                      PRESS = press_values,
                      row.names = model_names))
  }
}

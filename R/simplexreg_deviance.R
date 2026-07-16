################################################################################
#                SIMPLEX REGRESSION - DEVIANCE FUNCTION                        #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2026-05-26                                                             #
# Description: Computes the unit deviance for the simplex distribution         #
################################################################################

#' @title Unit Deviance Function of the Simplex Distribution
#' @description Computes the unit deviance (scaled deviance component) for the
#' simplex distribution.
#'
#' @param y Numeric scalar or vector of observed values (\eqn{0 < y < 1}).
#' If a vector, must be the same length as \code{mu} or recyclable.
#' @param mu Numeric scalar or vector of mean values (\eqn{0 < \mu < 1}).
#' If a vector, must be the same length as \code{y} or recyclable.
#'
#' @details
#' The unit deviance for the simplex distribution is defined as:
#' \deqn{d(y, \mu) = \frac{(y - \mu)^2}{y(1-y)[\mu(1-\mu)]^2}.}
#'
#' This function is used internally in maximum likelihood estimation and
#' model diagnostics for simplex regression.
#'
#' @return A numeric scalar or vector of unit deviance values.
#'
#' @examples
#' # Single value
#' dev.unit.simplex(y = 0.6, mu = 0.5)
#'
#' # Vector of values
#' y_vec <- c(0.2, 0.5, 0.8)
#' mu_vec <- c(0.3, 0.5, 0.7)
#' dev.unit.simplex(y = y_vec, mu = mu_vec)
#'
#' # Perfect fit returns zero deviance
#' dev.unit.simplex(y = 0.5, mu = 0.5)
#'
#' @references
#' Jørgensen, B. (1997). \emph{The Theory of Dispersion Models}.
#' Chapman and Hall, London.
#'
#' Song, P. X.-K. and Tan, M. (2000).
#' Marginal models for longitudinal continuous proportional data.
#' \emph{Biometrics}, \bold{56}(2), 496--502.
#' \doi{10.1111/j.0006-341X.2000.00496.x}
#'
#' @seealso \code{\link{variance.simplex}}, \code{\link{dsimplex}}.
#' @export
dev.unit.simplex <- function(y, mu){
  # Input validation
  stopifnot(
    "'y' must always be in (0, 1)" = all(y > 0 & y < 1),
    "parameter 'mu' must always be in (0, 1)" = all(mu > 0 & mu < 1)
  )

  # Compute unit deviance
  diff <- y - mu
  yoneminy <- y * (1 - y)
  muonemu <- mu * (1 - mu)
  deviance <- (diff / muonemu)^2 / yoneminy

  return(deviance)
}

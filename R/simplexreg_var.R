################################################################################
#                 SIMPLEX REGRESSION - VARIANCE FUNCTION                       #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2026-05-26                                                             #
# Description: Computes the variance function for the simplex distribution     #
################################################################################

#' @title Variance Function of the Simplex Distribution
#' @description Computes the variance of the simplex distribution as a function
#' of the mean parameter \eqn{\mu} and dispersion parameter \eqn{\sigma^2}.
#'
#' @param mu Numeric scalar or vector of mean parameters (\eqn{0 < \mu < 1}).
#' If a vector, must be the same length as \code{sigma2} or recyclable.
#' @param sigma2 Numeric scalar or vector of dispersion parameters
#' (\eqn{\sigma^2 > 0}).
#' If a vector, must be the same length as \code{mu} or recyclable.
#'
#' @details
#' The variance function for the simplex distribution is given by:
#' \deqn{Var(Y) = \mu(1-\mu) - \frac{1}{\sqrt{2\sigma^2}} \exp(a) \Gamma(0.5, a),}
#' where \eqn{a = \frac{1}{2\sigma^2[\mu(1-\mu)]^2}} and \eqn{\Gamma(0.5, a)}
#' is the upper incomplete gamma function.
#'
#' For large values of \eqn{a} (> 700), an asymptotic approximation is used
#' to avoid numerical overflow:
#' \deqn{Var(Y) \approx \mu(1-\mu) - \frac{1}{\sqrt{2\sigma^2}} \sqrt{\frac{1}{a}}.}
#'
#' @return A numeric scalar or vector of variance values.
#'
#' @examples
#' # Single value
#' variance.simplex(mu = 0.5, sigma2 = 0.1)
#'
#' # Vector of values
#' mu_vec <- c(0.3, 0.5, 0.7)
#' sigma2_vec <- c(0.1, 0.15, 0.2)
#' variance.simplex(mu = mu_vec, sigma2 = sigma2_vec)
#'
#' @references
#' Jørgensen, B. (1997).
#' \emph{The Theory of Dispersion Models}.
#' Chapman and Hall, London.
#'
#' Song, P. X.-K. and Tan, M. (2000).
#' Marginal models for longitudinal continuous proportional data.
#' \emph{Biometrics}, \bold{56}(2), 496--502.
#' \doi{10.1111/j.0006-341X.2000.00496.x}
#'
#' @importFrom expint gammainc
#'
#' @export
variance.simplex <- function(mu, sigma2) {
  # Input validation
  stopifnot(
    "parameter 'mu' must always be in (0, 1)" = all(mu > 0 & mu < 1),
    "parameter 'sigma2' must always be positive" = all(sigma2 > 0)
  )

  term1 <- mu * (1 - mu)
  a <- 1/(2*sigma2*term1^2)
  
  # Compute adjustment term with numerical stability
  # For large a (> 700), use asymptotic approximation to avoid overflow
  term2 <- ifelse(a <= 700, 
                  1/sqrt(2*sigma2) * exp(a) * gammainc(0.5, a),
                  1/sqrt(2*sigma2) * sqrt(1 / a))
  
  tol <- 1e-10
  simplex_variance <- ifelse(term2 < term1 - tol,
                            term1 - term2,
                            term1 - .Machine$double.eps)
  return(simplex_variance)
}

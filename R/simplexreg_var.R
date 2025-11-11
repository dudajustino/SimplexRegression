################################################################################
#                 SIMPLEX REGRESSION - VARIANCE FUNCTION                       #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2025-11-08                                                             #
# Description: Computes the variance function for the simplex distribution     #
################################################################################

# ==============================================================================
# SIMPLEX VARIANCE FUNCTION
# ==============================================================================

#' @title Simplex Variance Function
#' @description Computes the variance of the simplex distribution as a function
#' of the mean parameter \eqn{\mu} and dispersion parameter \eqn{\sigma^2}.
#'
#' @param mu Numeric vector of mean parameters (\eqn{0 < \mu < 1})
#' @param sigma2 Numeric vector of dispersion parameters (\eqn{\sigma^2 > 0})
#'
#' @details
#' The variance function for the simplex distribution is given by:
#' \deqn{Var(Y) = \mu(1-\mu) - \frac{1}{\sqrt{2\sigma^2}} \exp(a) \Gamma(0.5, a)}
#' where \eqn{a = \frac{1}{2\sigma^2[\mu(1-\mu)]^2}} and \eqn{\Gamma(0.5, a)}
#' is the upper incomplete gamma function.
#'
#' For large values of \eqn{a} (> 700), an asymptotic approximation is used
#' to avoid numerical overflow:
#' \deqn{Var(Y) \approx \mu(1-\mu) - \frac{1}{\sqrt{2\sigma^2}} \sqrt{\frac{1}{a}}}
#'
#' @return Numeric vector of variance values
#'
#' @importFrom expint gammainc
#'
#'@references
#' O. E. Barndorff-Nielsen, B. Jørgensen, Some parametric models on the simplex, Journal of
#' Multivariate Analysis 39 (1991) 106–116. https://doi.org/10.1016/0047-259X(91)90008-P.
#'
#' @examples
#' # Single value
#' variance_simplexreg(mu = 0.5, sigma2 = 0.1)
#'
#' # Vector of values
#' mu_vec <- c(0.3, 0.5, 0.7)
#' sigma2_vec <- c(0.1, 0.15, 0.2)
#' variance_simplexreg(mu = mu_vec, sigma2 = sigma2_vec)
#'
#' @export
variance_simplexreg <- function(mu, sigma2) {
  # Input validation
  if(any(mu <= 0 | mu >= 1)) {
    stop("mu must be in the interval (0, 1)")
  }
  if(any(sigma2 <= 0)) {
    stop("sigma2 must be positive")
  }

  term1 <- mu * (1 - mu)
  a = 1/(2*sigma2*term1^2)
  # Compute adjustment term with numerical stability
  # For large a (> 700), use asymptotic approximation to avoid overflow
  term2 = ifelse(a <= 700, 1/sqrt(2*sigma2) * exp(a) * expint::gammainc(0.5, a),
                 1/sqrt(2*sigma2) * sqrt(1 / a))
  simplex.variance = ifelse(term2 < term1,
                            term1 - term2,
                            term1 - .Machine$double.eps)
  return(simplex.variance)
}

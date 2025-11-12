################################################################################
#                SIMPLEX REGRESSION - LOG-LIKELIHOOD FUNCTION                  #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2025-11-08                                                             #
# Description: Computes the log-likelihood function for the simplex regression #
#              models with parametric or fixed mean link functions            #
################################################################################

#' @title Log-likelihood for the Simplex Regression
#' @description Computes the log-likelihood function for given parameters in the
#'   simplex regression models with parametric or fixed mean link functions.
#' @param theta Parameter vector
#' @param y Response vector (0 < y < 1)
#' @param x1 Design matrix for mean model (with intercept)
#' @param z1 Design matrix for dispersion model (with intercept)
#' @param link.mu Mean link function ("logit", "probit", "cloglog", "loglog",
#' "cauchit", "plogit1" or "plogit2")
#' @param link.sigma2 Dispersion link function ("log", "sqrt", or "identity")
#'
#' @details
#' The log-likelihood function for the simplex regression model with
#' parametric or fixed mean link function is given by
#' \deqn{
#'   \ell(\boldsymbol{\theta}) =
#'   -\frac{1}{2} \sum_{i=1}^{n}
#'   \left[
#'     \log(2\pi) +
#'     \log(\sigma_i^2) +
#'     3\log(y_i(1 - y_i)) +
#'     \frac{d(y_i; \mu_i)}{\sigma_i^2}
#'   \right],
#' }
#' where \eqn{d(y_i; \mu_i) = \frac{(y_i - \mu_i)^2}{y_i(1 - y_i)\mu_i^2(1 - \mu_i)^2}} is the
#' deviance component of the simplex model, \eqn{\mu_i} is obtained from the mean link function,
#' and \eqn{\sigma_i^2} from the dispersion link function.
#'
#' The parameter vector \eqn{\boldsymbol{\theta}} has different structures
#' depending on the type of mean link:
#' \itemize{
#'   \item For parametric mean link functions ("plogit1" or "plogit2"):
#'     \eqn{\boldsymbol{\theta} = (\boldsymbol{\beta}^\top, \boldsymbol{\gamma}^\top, \lambda)^\top},
#'     with \eqn{p + q + 1} elements.
#'   \item For fixed mean link functions ("logit", "probit", "cloglog", "loglog", "cauchit"):
#'     \eqn{\boldsymbol{\theta} = (\boldsymbol{\beta}^\top, \boldsymbol{\gamma}^\top)^\top},
#'     with \eqn{p + q} elements.
#' }
#'
#' Here, \eqn{\boldsymbol{\beta}} and \eqn{\boldsymbol{\gamma}} are regression
#' coefficients for the mean and dispersion submodels, respectively.
#' @return Numeric scalar (log-likelihood)
#' @export
simplex_loglik <- function(theta, y, x1, z1, link.mu, link.sigma2) {
  p <- ncol(x1)
  q <- ncol(z1)

  beta <- as.vector(theta[1:p])
  delta <- as.vector(theta[(p+1):(p+q)])

  eta1 <- as.vector(x1%*%beta)
  eta2 <- as.vector(z1%*%delta)

  if (link.mu %in% c("plogit1", "plogit2")) {
    lambda <- as.numeric(theta[p + q + 1])
    mu <- as.vector(parametric_mean_link_inv(eta1, lambda, link.mu))
  } else {
    mu <- as.vector(fixed_mean_link_inv(eta1, link.mu))
  }

  sigma2 <- as.vector(dispersion_link_inv(eta2, link.sigma2))

  dev <- deviance_simplexreg(y, mu)

  # Log-likelihood function
  ll <- -0.5 * sum(log(2 * pi) + log(sigma2) + 3 * log(y * (1 - y)) + dev / sigma2)

  return(ll)
}

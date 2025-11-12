################################################################################
#                    SIMPLEX REGRESSION - SCORE FUNCTION                       #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2025-11-08                                                             #
# Description: Computes the log-likelihood function for the simplex regression #
#              models with parametric or fixed mean link functions             #
################################################################################

#' @title Score Vector for the Simplex Regression
#' @description Computes the score vector for given parameters in the
#'   simplex regression models with parametric or fixed mean link functins.
#' @param theta Parameter vector
#' @param y Response vector (0 < y < 1)
#' @param x1 Design matrix for mean model (with intercept)
#' @param z1 Design matrix for dispersion model (with intercept)
#' @param link.mu Mean link function ("logit", "probit", "cloglog", "loglog",
#' "cauchit", "plogit1" or "plogit2")
#' @param link.sigma2 Dispersion link function ("log", "sqrt", or "identity")
#'
#' @details
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
#'
#' @return Numeric vector (gradient)
#' @export
simplex_score <- function(theta, y, x1, z1, link.mu, link.sigma2) {
  p <- ncol(x1)
  q <- ncol(z1)

  beta <- as.vector(theta[1:p])
  delta <- as.vector(theta[(p+1):(p+q)])

  eta1 <- as.vector(x1%*%beta)
  eta2 <- as.vector(z1%*%delta)

  if (link.mu %in% c("plogit1", "plogit2")) {
    lambda <- as.numeric(theta[p + q + 1])
    mu <- as.vector(parametric_mean_link_inv(eta1, lambda, link.mu))
    T1 <- diag(as.vector(parametric_mean_link_inv_deriv1(eta1, lambda, link.mu)))
  } else {
    mu <- as.vector(fixed_mean_link_inv(eta1, link.mu))
    T1 <- diag(as.vector(fixed_mean_link_inv_deriv1(eta1, link.mu)))
  }

  sigma2 <- as.vector(dispersion_link_inv(eta2, link.sigma2))
  T2 <- diag(as.vector(dispersion_link_inv_deriv1(eta2, link.sigma2)))

  diff <- as.vector(y - mu)
  muonemu <- as.vector(mu * (1 - mu))

  dev <- deviance_simplexreg(y, mu)

  U <- diag(as.vector(dev / (sigma2 * muonemu) + 1 / (sigma2 * (muonemu)^3)))
  a <- as.vector((-1 / (2 * sigma2)) + (dev / (2 * sigma2^2)))

  Ubeta <- t(x1) %*% T1 %*% U %*% diff
  Udelta <- t(z1) %*% T2 %*% a

  if (link.mu %in% c("plogit1", "plogit2")) {
    # Compute rho based on link function
    if(link.mu == "plogit2"){
      exp_aval_frac <- plogis(eta1)^(1/lambda)
      rho <- as.vector(- exp_aval_frac * (log(plogis(eta1))) / (lambda^2))
    } else{
      rho <- as.vector((-1/(lambda^2)) * ((exp(eta1) + 1) ^ (-1/lambda)) * (log(exp(eta1) + 1)))
    }
    Ulambda <- t(rho) %*% U %*% diff
    return(as.vector(c(Ubeta, Udelta, Ulambda)))
  } else{
    return(as.vector(c(Ubeta, Udelta)))
  }
}

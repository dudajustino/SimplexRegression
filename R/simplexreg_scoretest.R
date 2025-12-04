################################################################################
#                     SIMPLEX REGRESSION - SCORE TEST                          #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2025-11-08                                                             #
# Description: Rao score test for testing the link parameter lambda = 1        #
################################################################################

# ==============================================================================
# SCORE TEST
# ==============================================================================

#' @title Rao Score Test for Simplex Regression with Parametric Mean Link Function
#' @description Performs a Rao score test to test whether the link parameter
#' \eqn{\lambda} equals 1, which corresponds to testing between different
#' parametric mean link functions.
#'
#' @param model An object of class \code{"simplexregression"}
#' @param link.mu Character string specifying the link function under the
#' alternative hypothesis. Options are \code{"plogit1"} or \code{"plogit2"}
#'
#' @details
#' The score test (also known as the Rao test or Lagrange multiplier test)
#' tests the null hypothesis \eqn{H_0: \lambda = 1} against the alternative
#' that \eqn{\lambda \neq 1}.
#'
#' The test statistic is based on the score function and the information matrix
#' evaluated at the restricted maximum likelihood estimate under the null
#' hypothesis. The advantage of the score test is that it only requires fitting
#' the model under the null hypothesis.
#'
#' The test statistic is:
#' \deqn{S_R = \boldsymbol{U}_\lambda(\boldsymbol{\theta})^\top \boldsymbol{K}^{\lambda \lambda}
#' \boldsymbol{U}_\lambda(\boldsymbol{\theta})}
#'
#' where \eqn{\boldsymbol{U}_\lambda} is the score function with respect to \eqn{\lambda}
#' and \eqn{\boldsymbol{K}^{\lambda \lambda}} is the corresponding element of the
#' inverse Fisher information matrix.
#'
#' Under the null hypothesis, the test statistic follows a chi-squared
#' distribution with 1 degree of freedom.
#'
#' @return An object of class \code{"htest"} containing:
#' \itemize{
#'   \item \code{statistic}: The score test statistic
#'   \item \code{parameter}: Degrees of freedom (always 1)
#'   \item \code{p.value}: The p-value of the test
#'   \item \code{method}: Description of the test
#'   \item \code{data.name}: Model name and link function being tested
#' }
#'
#' @references
#' Rao, C. R. (1948). Large sample tests of statistical hypotheses concerning
#' several parameters with applications to problems of estimation.
#' \emph{Mathematical Proceedings of the Cambridge Philosophical Society},
#' 44(1), 50-57.
#'
#' @examples
#' \dontrun{
#' # Fit a simplex regression model with plogit1 link
#' model <- simplexregpar(y ~ x1 + x2, link.mu = "plogit1", data = mydata)
#'
#' # Test if lambda = 1 (i.e., test against plogit2 link)
#' scoretest.simplexregpar(model, link.mu = "plogit2")
#' }
#'
#' @importFrom stats pchisq plogis
#'
#' @export
scoretest.simplexreg <- function(model, link.mu = c("plogit1", "plogit2")) {

  METHOD = "Rao score test"
  DNAME  = deparse(substitute(model))
  DNAME  = paste(DNAME, "vs", link.mu)

  y <- as.vector(model$y)
  mu <- as.vector(model$mu.fv)
  sigma2 <- as.vector(model$sigma2.fv)
  lambda <- 1
  eta1 <- as.vector(model$mu.lp)
  mu_x <- model$mu.x

  diff <- as.vector(y - mu)

  yoneminy <- as.vector(y * (1 - y))
  muonemu <- as.vector(mu * (1 - mu))
  dev <- (diff / muonemu)^2 / yoneminy

  if(link.mu == "plogit2"){
    exp_aval_frac <- plogis(eta1)^(1/lambda)
    rho <- as.vector(- exp_aval_frac * (log(plogis(eta1))) / (lambda^2))
  } else{
    rho <- as.vector((-1/(lambda^2)) * ((exp(eta1) + 1) ^ (-1/lambda)) * (log(exp(eta1) + 1)))
  }

  Ui <- (dev / muonemu) + (1 / (muonemu^3))
  Ulambda <- sum(diff * Ui * rho / sigma2)

  l1id1 <- as.vector(fixed_mean_link_inv_deriv1(eta1, model$mu.link))
  wi <- (3 * sigma2) / muonemu + (1 / (muonemu^3))

  Klambdalambda <- sum(wi * rho^2 / sigma2)
  Klambdabeta <- colSums(wi * rho * l1id1 / sigma2 * mu_x)
  Kbetabeta <- crossprod(mu_x, (wi * l1id1^2 / sigma2) * mu_x)

  vcovlambda_inv <- Klambdalambda - sum(Klambdabeta * solve(Kbetabeta, Klambdabeta))

  S <- Ulambda^2 / vcovlambda_inv
  gl = 1
  PVAL = pchisq(S, gl, lower.tail= F)

  names(gl) = "df"
  names(S) = "S"

  RVAL <- list(statistic = S, parameter = gl, p.value = PVAL,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"

  return(RVAL)
}

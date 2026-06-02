################################################################################
#                     SIMPLEX REGRESSION - RAO SCORE TEST                      #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2026-05-26                                                             #
# Description: Rao score test for testing the link parameter lambda = 1        #
################################################################################

#' @title Rao's Score Test for Simplex Regression with Parametric Mean Link Function
#' @description Performs a Rao's score test to test whether the link parameter
#' \eqn{\lambda} equals 1, which corresponds to evaluating the null hypothesis
#' that the mean link function is the standard logit.
#'
#' @param model An object of class \code{simplexregression}.
#' @param link.mu Character string specifying the link function under the
#' alternative hypothesis. Options are "plogit1" or "plogit2".
#'
#' @details
#' Given that the fixed logit link function is a particular case of the parametric
#' logit link functions when \eqn{\lambda = 1}, it is possible to test whether the
#' mean link function is logit by testing \eqn{H_0: \lambda = 1} against
#' \eqn{H_1: \lambda \neq 1}.
#'
#' The score test statistic for this hypothesis test is given by:
#'
#' \deqn{S_R = \boldsymbol{U}_\lambda(\boldsymbol{\tilde{\theta}})^\top
#' \boldsymbol{\tilde{K}}^{\lambda \lambda}
#' \boldsymbol{U}_\lambda(\boldsymbol{\tilde{\theta}}),}
#'
#' where \eqn{\boldsymbol{U}_\lambda} and \eqn{\boldsymbol{\tilde{K}}^{\lambda \lambda}}
#' denote, respectively, the component of the score vector and the corresponding
#' element of the inverse Fisher information matrix associated with \eqn{\lambda},
#' both evaluated at \eqn{\boldsymbol{\tilde{\theta}} = (\boldsymbol{\tilde{\beta}}^\top,
#' \boldsymbol{\tilde{\gamma}}^\top, 1)^\top}, the maximum likelihood estimator under
#' the null hypothesis. For more details see Justino and Cribari-Neto (2025).
#'
#' Under regularity conditions, the null hypothesis, and when \eqn{n} is large,
#' the test statistic follows a chi-squared distribution with 1 degree of freedom.
#'
#' @return An object of class \code{"htest"} containing:
#' \itemize{
#'   \item \code{statistic}: The score test statistic,
#'   \item \code{parameter}: Degrees of freedom (always 1),
#'   \item \code{p.value}: The p-value of the test,
#'   \item \code{method}: Description of the test,
#'   \item \code{data.name}: Model name and link function being tested.
#' }
#'
#' @references
#' Justino, M. E. C. and Cribari-Neto, F. (2026).
#' Simplex regression with a flexible logit link: Inference and application
#' to cross-country impunity data.
#' \emph{Applied Mathematical Modelling}, \bold{154}, 116713. \doi{10.1016/j.apm.2025.116713}
#'
#' Rao, C. R. (1948). Large sample tests of statistical hypotheses concerning
#' several parameters with applications to problems of estimation.
#' \emph{Mathematical Proceedings of the Cambridge Philosophical Society},
#'  \bold{44}(1), 50--57. \doi{10.1017/S0305004100023987}
#'
#' @examples
#' # Simulate data with plogit2
#' n <- 100
#' x1 <- runif(n, 0, 1)
#' x2 <- runif(n, 0, 1)
#' mu <- parametric_mean_link_inv(0.8 - 1.2*x1 - 1.5*x2 , 0.25, "plogit2")
#' sigma2 <- 0.5
#' y <- rsimplex(n, mu, sigma2)
#' data <- data.frame(y = y, x1 = x1, x2 = x2)
#'
#' # Fit model with logit
#' model <- simplexreg(y ~ x1 + x2 | 1, data = data,
#'                      link.mu = "logit", link.sigma2 = "identity")
#'
#' # Test if lambda = 1
#' scoretest(model, link.mu = "plogit2")
#'
#' @importFrom stats pchisq plogis
#'
#' @export
scoretest <- function(model, link.mu = c("plogit1", "plogit2")) {
  
  if (!inherits(model, "simplexregression")) {
    stop("'model' must be an object of class 'simplexregression'")
  }

  link.mu <- match.arg(link.mu)

  # Validate that the fitted model uses logit link
  if (model$mu.link != "logit") {
    stop("The score test is only valid for models fitted with the logit mean link function. ",
         "The null hypothesis H0: lambda = 1 corresponds to the logit link, which must ",
         "be the link function of the model under H0. ",
         "Model provided uses link: '", model$mu.link, "'.")
  }

  METHOD <- "Rao score test"
  DNAME <- paste("Logit vs", link.mu)

  y <- as.vector(model$y)
  mu <- as.vector(model$mu.fv)
  sigma2 <- as.vector(model$sigma2.fv)
  lambda <- 1 # Under H0
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
    rho <- as.vector((-1/(lambda^2)) * ((exp(eta1) + 1) ^ (-1/lambda)) *
                       (log(exp(eta1) + 1)))
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
  df <- 1L
  PVAL <- pchisq(S, df, lower.tail = FALSE)

  names(df) <- "df"
  names(S) <- "S"

  RVAL <- list(statistic = S, parameter = df, p.value = PVAL,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"

  return(RVAL)
}

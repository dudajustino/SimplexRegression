################################################################################
#                      SIMPLEX REGRESSION - RESIDUALS                          #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2026-05-26                                                             #
# Description: Diagnostic residuals                                            #
################################################################################

#' @title Residuals for Simplex Regression Models
#' @description Extracts various types of residuals for diagnostic analysis in
#' simplex regression models with parametric or fixed mean link functions.
#'
#' @param object An object of class \code{simplexregression}.
#' @param type Character string specifying the type of residual to extract.
#' Options:
#' \code{"quantile"} (default), \code{"pearson"}, \code{"deviance"},
#' \code{"standardized"}, \code{"weighted"}, \code{"variance"},
#' \code{"biasvariance"}, \code{"score"}, \code{"dualscore"}, \code{"response"}.
#' @param ... Additional arguments (currently not used).
#'
#' @return A numeric vector of residuals.
#'
#' @details
#' Several types of residuals are available for model diagnostics:
#'
#' \strong{Quantile residuals} (\code{"quantile"}): Proposed by Dunn and Smyth
#' (1996) as \eqn{r_i^Q = \Phi^{-1}(F(y_i; \hat{\mu}_i, \hat{\sigma}^2_i))},
#' where \eqn{\Phi(\cdot)} is the standard normal CDF and
#' \eqn{F(\cdot; \cdot)} is the simplex CDF (see \code{\link{psimplex}}).
#' Under correct model specification, these residuals are approximately standard
#' normal and are therefore recommended for general diagnostic use.
#'
#' \strong{Pearson residuals} (\code{"pearson"}): Defined in McCullagh and
#' Nelder (1989) as \eqn{r_i^P = (y_i - \hat{\mu}_i) /
#' \sqrt{\widehat{\text{Var}}(y_i)}}, where \eqn{\widehat{\text{Var}}(y_i)} is
#' the estimated variance of the response
#' (see \code{\link{variance.simplex}}).
#'
#' \strong{Deviance residuals} (\code{"deviance"}): Defined in Jørgensen (1997,
#' p. 115) as \eqn{r_i^D = (y_i - \hat{\mu}_i) /
#' (\hat{\mu}_i(1-\hat{\mu}_i)\sqrt{y_i(1-y_i)})}.
#'
#' \strong{Standardized residuals} (\code{"standardized"}): Proposed by Espinheira
#' and Silva (2020, Eq. 15) as \eqn{r_i^\beta = \hat{u}_i(y_i - \hat{\mu}_i) /
#' \sqrt{\hat{\sigma}^2_i \hat{w}_i}}, where
#' \eqn{\hat{w}_i} and \eqn{\hat{u}_i} are weight functions.
#'
#' \strong{Weighted residuals} (\code{"weighted"}): Proposed by Espinheira and
#' Silva (2020, Eq. 16) as \eqn{r_i^{\beta*} = \hat{u}_i(y_i - \hat{\mu}_i) /
#' \sqrt{\hat{\sigma}^2_i \hat{w}_i(1-\hat{h}_{ii})}}, where \eqn{\hat{h}_{ii}}
#' are the diagonal elements of the hat matrix (see
#' \code{\link{hatvalues.simplexregression}}).
#' These residuals are recommended for simulated envelope plots.
#'
#' \strong{Variance residuals} (\code{"variance"}): Proposed by Espinheira et al.
#' (2021, Eq. 7) as \eqn{r_i^\gamma = (\hat{d}_i - \hat{\sigma}^2_i) /
#' (\hat{\sigma}^2_i\sqrt{2})}, where \eqn{\hat{d}_i} is the estimated unit
#' deviance.
#'
#' \strong{Bias–variance residuals} (\code{"biasvariance"}): Proposed by Espinheira
#' et al. (2021, Eq. 8) as \eqn{r_i^{\beta \gamma} =
#' (\hat{u_i}(y_i - \hat{\mu}_i) + \hat{a}_i) / \sqrt{\hat{\sigma}^2_i
#' \hat{w}_i + 1/(2\hat{\sigma}^4_i)}}, where \eqn{\hat{a}_i} is a correction term.
#'
#' \strong{Score residuals} (\code{"score"}): Defined in Jørgensen (1997, p. 115) as
#' \eqn{r_i^{S} = (y_i - \hat{\mu}_i)(\hat{\mu}_i^2 + y_i - 2y_i\hat{\mu}_i) /
#' (y_i(1-y_i)\hat{\mu}_i^{1.5}(1-\hat{\mu}_i)^{1.5})}.
#'
#' \strong{Dual score residuals} (\code{"dualscore"}): Defined in Jørgensen
#' (1997, p. 115) as \eqn{r_i^{DS} = (y_i - \hat{\mu}_i)(y_i + \hat{\mu}_i -
#' 2y_i\hat{\mu}_i) / (2\sqrt{y_i(1-y_i)}\hat{\mu}_i^2(1-\hat{\mu}_i)^2)}.
#'
#' \strong{Response residuals} (\code{"response"}): Simple difference between
#' observed and fitted values, \eqn{r_i^R = y_i - \hat{\mu}_i}.
#'
#' \strong{Recommendations}: Quantile residuals are recommended for general model
#' diagnostics due to their theoretical properties. Weighted residuals are
#' particularly useful for constructing simulated envelope plots, as they account
#' for both the variance structure and leverage effects.
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
#' # Fit model with parametric mean link functions
#' fit <- simplexreg(y ~ x1 + x2 | z1, data = data, link.mu = "plogit1")
#'
#' # Compute different types of residuals
#' res_quantile <- residuals(fit, type = "quantile")
#' res_pearson <- residuals(fit, type = "pearson")
#' res_weighted <- residuals(fit, type = "weighted")
#'
#' @importFrom stats pnorm qnorm
#'
#' @references
#'
#' McCullagh, P. and Nelder, J. A. (1989). \emph{Generalized Linear Models}.
#' 2nd ed.  Chapman and Hall, London.

#' Jørgensen, B. (1997).
#' \emph{The Theory of Dispersion Models}.
#' Chapman and Hall, London.
#'
#' Dunn, P. K. and Smyth, G. K. (1996). Randomized Quantile Residuals.
#' \emph{Journal of Computational and Graphical Statistics}, \bold{5}(3),
#' 236-–244. \doi{10.2307/1390802}
#'
#' Espinheira, P. L., Silva, L. C. M. and Cribari-Neto, F. (2021). Bias and
#' variance residuals for machine learning nonlinear simplex regressions.
#' \emph{Expert Systems With Applications}, \bold{185}, 115656.
#' \doi{10.1016/j.eswa.2021.115656}
#'
#' Espinheira, P. L., Silva, A. O. (2020). Residual and influence analysis to a
#' general class of simplex regression. \emph{TEST}, \bold{29}, 523–-552.
#' \doi{10.1007/s11749-019-00665-3}
#'
#' @seealso \code{\link{plot.simplexregression}}, \code{\link{halfnormal.plot}},
#' \code{\link{hatvalues.simplexregression}}.
#'
#' @export
residuals.simplexregression <- function(object, type = c("quantile", "pearson",
                                                     "deviance", "standardized",
                                                     "weighted", "variance",
                                                     "biasvariance", "score",
                                                     "dualscore", "response"), ...) {

  if (!inherits(object, "simplexregression")) {
    stop("'object' must be an object of class 'simplexregression'")
  }

  type <- match.arg(type)

  # Quantile residuals are already stored in the object
  if (type == "quantile") {
    return(object$residuals)
  }

  y <- as.vector(object$y)
  mu <- as.vector(object$mu.fv)
  sigma2 <- as.vector(object$sigma2.fv)

  diff <- y - mu
  muonemu <- mu * (1 - mu)
  yoneminy <- y * (1 - y)
  dev <- (diff / muonemu)^2 / yoneminy

  # Calculate weights for weighted residuals
  if(type %in% c("standardized", "weighted", "biasvariance")) {
    wi <- (3*sigma2 / muonemu) + (1 / (muonemu^3))
    ui <- (dev/muonemu) + (1/(muonemu^3))
  }

  res <- switch(type,

                "pearson" = {
                  var_hat <- variance.simplex(mu, sigma2)
                  diff / sqrt(var_hat)
                },

                "deviance" = {
                  diff / (muonemu*sqrt(yoneminy))
                },

                "standardized" = {
                  (ui*diff) / sqrt(sigma2*wi)
                },

                "weighted" = {
                  (ui*diff) / sqrt(sigma2*wi*(1-hatvalues(object)))
                },

                "variance" = {
                  (dev - sigma2) / (sigma2*sqrt(2))
                },

                "biasvariance" = {
                  ai <- - 1/(2*sigma2) + dev/(2*sigma2^2)
                  (ui*diff + ai) / sqrt(sigma2*wi + 1/(2*sigma2^2))
                },

                "response" = {
                  diff
                },

                "score" = {
                  (diff * (mu^2 + y - 2*y*mu)) / (yoneminy * muonemu^1.5)
                },

                "dualscore" = {
                  (diff * (y + mu - 2*y*mu)) / (2 * sqrt(yoneminy) * muonemu^2)
                })

  return(res)
}

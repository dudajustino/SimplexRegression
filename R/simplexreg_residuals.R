################################################################################
#                      SIMPLEX REGRESSION - RESIDUALS                          #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2025-11-08                                                             #
# Description: Diagnostic residuals                                            #
################################################################################

#' @title Residuals for Simplex Regression
#' @description Computes various types of residuals for simplex regression models
#' with Parametric or Fixed Mean Link Functions.
#'
#' @param object An object of class \code{"simplexregression"}
#' @param type Character string specifying the type of residual. Options:
#' \itemize{
#'   \item \code{"quantile"}: Quantile residuals (default)
#'   \item \code{"pearson"}: Pearson residuals
#'   \item \code{"pearson P"}: Standardized Pearson residuals
#'   \item \code{"deviance"}: Deviance residuals
#'   \item \code{"deviance P"}: Standardized deviance residuals
#'   \item \code{"standardized"}: Standardized residuals
#'   \item \code{"weighted"}: Weighted residuals
#'   \item \code{"variance"}: Variance residuals
#'   \item \code{"variance P"}: Standardized variance residuals
#'   \item \code{"biasvariance"}: Bias-variance residuals
#'   \item \code{"anscombe"}: Anscombe residuals
#'   \item \code{"williams"}: Williams residuals
#'   \item \code{"response"}: Response (ordinary) residuals
#'   \item \code{"score"}: Score residuals
#'   \item \code{"dualscore"}: Dual score residuals
#' }
#' @param ... Additional arguments (currently not used)
#'
#' @return Numeric vector of residuals
#'
#' @details
#' Different residual types are suitable for different diagnostic purposes:
#' \itemize{
#'   \item Quantile residuals follow standard normal distribution under correct
#'   model
#'   \item Pearson and deviance residuals are classical GLM residuals
#'   \item Standardized versions (with "P") account for leverage
#'   \item Standardized and Weighted residuals are useful for detecting mean model
#'   misspecification
#'   \item Variance residuals detect dispersion model misspecification
#' }
#'
#' @examples
#' \dontrun{
#' fit <- simplexreg(y ~ x | 1, data = mydata)
#' res_quantile <- residuals(fit, type = "quantile")
#' res_pearson <- residuals(fit, type = "pearson P")
#' }
#'
#' @importFrom stats pnorm qnorm
#' @importFrom expint gammainc
#' @export
residuals.simplexregression <- function(object, type = c("quantile", "pearson", "pearson P",
                                                     "deviance", "deviance P", "standardized",
                                                     "weighted", "variance", "variance P",
                                                     "biasvariance", "anscombe", "williams",
                                                     "response", "score", "dualscore"), ...) {
  type <- match.arg(type)

  y <- as.vector(object$y)
  mu <- as.vector(object$mu.fv)
  sigma2 <- as.vector(object$sigma2.fv)

  diff <- y - mu
  muonemu <- mu * (1 - mu)
  yoneminy <- y * (1 - y)
  dev <- (diff / muonemu)^2 / yoneminy

  # Response residuals
  if(type == "response") return(diff)

  # Calculate hat values only when needed
  hat_beta_needed <- type %in% c("pearson P", "deviance P", "weighted", "williams")
  hat_beta <- if(hat_beta_needed) {
    hatvalues.simplexregression(object)
  } else NULL

  # Calculate weights for weighted residuals
  if(type %in% c("standardized", "weighted", "biasvariance")) {
    wi <- (3*sigma2 / muonemu) + (1 / (muonemu^3))
    ui <- (dev/muonemu) + (1/(muonemu^3))
  }

  res <- switch(type,

                "quantile" = {
                  object$residuals
                },

                "pearson" = {
                  var_hat = muonemu - (1/sqrt(2*sigma2) * exp(1/(2*sigma2*muonemu^2)) *
                                         expint::gammainc(0.5, 1/(2*sigma2*muonemu^2)))
                  diff / sqrt(var_hat)
                },

                "pearson P" = {
                  diff / sqrt(sigma2*muonemu^3*(1-hat_beta))
                },

                "deviance" = {
                  diff / (muonemu*sqrt(yoneminy))
                },

                "deviance P" = {
                  diff / (muonemu*sqrt(sigma2*yoneminy*(1-hat_beta)))
                },

                "anscombe" = {
                  log((y*(1-mu)) / (mu*(1-y))) / (sqrt(muonemu))

                },

                "standardized" = {
                  (ui*diff) / sqrt(sigma2*wi)
                },

                "weighted" = {
                  (ui*diff) / sqrt(sigma2*wi*(1-hat_beta))
                },

                "variance" = {
                  (dev - sigma2) / (sigma2*sqrt(2))
                },

                "variance P" = {
                  Z <- object$sigma2.x
                  eta2 <- object$sigma2.lp

                  Vi <- 1/(2*sigma2^2)
                  Dlink <- dispersion_link_inv_deriv1(object$sigma2.lp,  object$sigma2.link)
                  ai <- -1/(2*sigma2) + dev*Vi

                  weights <- Vi * Dlink^2
                  Zw <- Z * sqrt(weights)
                  Inv <- chol2inv(chol(t(Zw) %*% Zw))

                  hat_gamma <- diag(Zw %*% Inv %*% t(Zw))

                  (dev - sigma2) / (sigma2*sqrt(2*(1-hat_gamma)))
                },

                "biasvariance" = {
                  ai <- - 1/(2*sigma2) + dev/(2*sigma2^2)
                  (ui*diff + ai) / sqrt(sigma2*wi + 1/(2*sigma2^2))
                },

                "score" = {
                  (diff * (mu^2 + y - 2*y*mu)) / (yoneminy * muonemu^1.5)
                },

                "dualscore" = {
                  (diff * (y + mu - 2*y*mu)) / (2 * sqrt(yoneminy) * muonemu^2)
                },

                "williams" = {
                  resP <- residuals.simplexregression(object, type="pearson P")
                  resD <- residuals.simplexregression(object, type="deviance P")
                  sign(diff) * sqrt((1-hat_beta)*(resD^2) + hat_beta*(resP^2))
                })

  return(res)
}

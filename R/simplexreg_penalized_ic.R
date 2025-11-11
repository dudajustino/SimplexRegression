################################################################################
#               PENALIZED INFORMATION CRITERIA FOR SIMPLEX REGRESSION          #
################################################################################

#' @title Penalized Information Criteria for Simplex Regression
#' @description Computes AIC, BIC, and HQ information criteria with a penalty term
#' depending on the estimated parameter \eqn{\lambda} of the parametric mean link
#' functions ("plogit1" or "plogit2").
#'
#' @param model An object of class \code{"simplexregression"} fitted with a
#' parametric mean link functions ("plogit1" or "plogit2").
#' @param kappa A numeric constant controlling the penalization strength
#'   (default is 0.1).
#'
#' @details
#' The penalized criteria are computed as:
#' \deqn{AIC^{(\lambda)} = -2 \ell + (2 + c \, |\log(\lambda)|)r}
#' \deqn{BIC^{(\lambda)} = -2 \ell + (\log(n) + c \, |\log(\lambda)|)r}
#' \deqn{HQIC^{(\lambda)}  = -2 \ell + (2 \log(\log(n)) + c \, |\log(\lambda)|)r}
#' where:
#' \itemize{
#'   \item \eqn{\ell} denotes the maximized log-likelihood function;
#'   \item \eqn{\lambda} is the parameter of the parametric mean link function;
#'   \item \eqn{r} is the number of parameters in the model.
#' }
#'
#' These penalized versions add a smooth penalty on the magnitude of
#' \eqn{\lambda}, encouraging simpler link structures.
#'
#' @return A named vector with components \eqn{AIC^{(\lambda)}}, \eqn{BIC^{(\lambda)}},
#' and \eqn{HQIC^{(\lambda)}}.
#'
#' @examples
#' \dontrun{
#' # Fit two models with parametric mean link functions
#' fit1 <- simplexreg(y ~ x1 + x2, link = "plogit1")
#' fit2 <- simplexreg(y ~ x1 + x2, link = "plogit2")
#'
#' # Compute penalized criteria
#' penalized_ic(fit1)
#' penalized_ic(fit2, kappa = 0.2)
#' }
#'
#' @export
penalized_ic <- function(model, kappa = 0.1) {

  if (is.null(model$lambda.fv)) {
    stop("Model does not have a parametric mean link (lambda.fv missing).")
  }

  n <- length(model$fitted.values)
  r <- n - model$df.residual
  lambda <- model$lambda.fv
  ll <- model$loglik
  c <- kappa

  # penalized terms
  aic_c <- -2 * ll + (2 + c * abs(log(lambda))) * r
  bic_c <- -2 * ll + (log(n) + c * abs(log(lambda))) * r
  hq_c  <- -2 * ll + (2 * log(log(n)) + c * abs(log(lambda))) * r

  out <- c(AICc = aic_c, BICc = bic_c, HQc = hq_c)
  return(out)
}

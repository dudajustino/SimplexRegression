################################################################################
#                     SIMPLEX REGRESSION - RESET TEST                          #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2026-05-26                                                             #
# Description: RESET (Regression Equation Specification Error Test) for        #
#              testing functional form misspecification                        #
################################################################################

#' @title RESET Test for Simplex Regression
#' @description Performs the Ramsey's RESET test for functional form in simplex
#' regression models.
#'
#' @param model An object of class \code{simplexregression}.
#' @param dispersion Logical. If \code{TRUE}, includes the squared linear predictor
#' in the dispersion submodel as well. Default is \code{TRUE}.
#'
#' @details
#' The RESET test augments the original model by adding the squared linear
#' predictor as an additional covariate. Under the null hypothesis of correct
#' functional form, this additional term should not be significant.
#'
#' The test statistic follows a chi-squared distribution with degrees of freedom
#' equal to the difference in the number of parameters between the augmented
#' and original models.
#'
#' If \code{dispersion = TRUE}, the squared linear predictor is added to both
#' the mean and dispersion submodels. If \code{FALSE}, it is only added to the
#' mean submodel.
#'
#' @return An object of class \code{"htest"} containing:
#' \itemize{
#'   \item \code{statistic}: The score test statistic,
#'   \item \code{parameter}: Degrees of freedom,
#'   \item \code{p.value}: The p-value of the test,
#'   \item \code{method}: Description of the test,
#'   \item \code{data.name}: Model name and link function being tested.
#' }
#'
#' @references
#' Ramsey, J. B. (1969). Tests for specification errors in classical linear
#' least-squares regression analysis. \emph{Journal of the Royal Statistical
#' Society: Series B}, 31(2), 350-371.
#'
#' @examples
#' \dontrun{
#' # Fit a simplex regression model
#' model <- simplexreg(y ~ x1 + x2, data = mydata)
#'
#' # Perform RESET test
#' resettest(model)
#'
#' # RESET test only for mean submodel
#' resettest(model, dispersion = FALSE)
#' }
#'
#' @importFrom stats pchisq
#'
#' @seealso \code{\link{simplexreg}}.
#' @export
resettest <- function(model, dispersion = TRUE){

  if (!inherits(model, "simplexregression")) {
    stop("'model' must be an object of class 'simplexregression'")
  }
  
  METHOD = "RESET test"
  DNAME <- paste(deparse(model$formula), collapse = "")

  y <- as.vector(model$y)
  mu_lp_squared <- model$mu.lp^2

  x <- cbind(model$mu.x, mu_lp_squared)
  
  if(dispersion == TRUE){
    z <- cbind(model$sigma2.x, mu_lp_squared)
  } else {
    z <- model$sigma2.x
  }

  ctrl <- model$control
  ctrl$hessian <- FALSE

  modelh1 <- simplexreg.fit(y, x[, colnames(x) != "(Intercept)", drop = FALSE], 
                            z[, colnames(z) != "(Intercept)", drop = FALSE],
                            link.mu = model$mu.link, link.sigma2 = model$sigma2.link,
                            x_names = colnames(x), z_names = colnames(z), control = ctrl,
                            weights = model$weights)

  lH0 <- model$loglik
  lH1 <- modelh1$loglik
  LR <- 2*(lH1 - lH0)
  df <- model$df.residual - modelh1$df.residual
  PVAL <- pchisq(LR, df, lower.tail = FALSE)

  names(df) <- "df"
  names(LR) <- "RESET"

  RVAL <- list(statistic = LR, parameter = df, p.value = PVAL,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"

  return(RVAL)
}

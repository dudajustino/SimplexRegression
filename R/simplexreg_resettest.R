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
#' @param dispersion Logical. If \code{TRUE}, includes the augmented terms
#' in the dispersion submodel as well. Default is \code{TRUE}.
#' @param power Integer vector specifying which powers of the linear predictor
#' (or fitted values) to include as additional regressors. Default is \code{2}
#' (squared term only). Use \code{power = 2:3} to include both squared and
#' cubic terms, following the convention of \code{lmtest::resettest}.
#' @param type Character string specifying the base for the augmented terms.
#' \code{"lp"} (default) uses the mean linear predictor \eqn{\hat{\eta}_1};
#' \code{"fitted"} uses the fitted mean values \eqn{\hat{\mu}} on the (0, 1)
#' scale. For models with parametric link functions, \code{"lp"} is generally
#' preferred as it operates on an unrestricted scale.
#'
#' @details
#' The RESET test augments the original model by adding powers of the linear
#' predictor (or fitted values) as additional covariates. Under the null
#' hypothesis of correct functional form, these additional terms should not
#' be significant.
#'
#' The likelihood ratio statistic follows a chi-squared distribution with
#' degrees of freedom equal to the number of augmented terms added (i.e.,
#' \code{length(power)} if \code{dispersion = FALSE}, or
#' \code{2 * length(power)} if \code{dispersion = TRUE}).
#'
#' If \code{dispersion = TRUE}, the augmented terms are added to both the
#' mean and dispersion submodels. If \code{FALSE}, they are only added to
#' the mean submodel.
#'
#' @return An object of class \code{"htest"} containing:
#' \itemize{
#'   \item \code{statistic}: The likelihood ratio test statistic,
#'   \item \code{parameter}: Degrees of freedom,
#'   \item \code{p.value}: The p-value of the test,
#'   \item \code{method}: Description of the test,
#'   \item \code{data.name}: Model formula.
#' }
#'
#' @references
#' Ramsey, J. B. (1969). Tests for specification errors in classical linear
#' least-squares regression analysis. \emph{Journal of the Royal Statistical
#' Society: Series B}, 31(2), 350-371.
#'
#' @examples
#' data(ReadingSkills, package = "SimplexRegression")
#' fit <- simplexreg(accuracy ~ dyslexia * iq | dyslexia + iq + I(iq^2),
#'                  data = ReadingSkills)
#'
#' # Default: squared linear predictor in both submodels
#' resettest(fit)
#'
#' # RESET test only for mean submodel
#' resettest(fit, dispersion = FALSE)
#'
#' # Include squared and cubic terms
#' resettest(fit, power = 2:3)
#'
#' # Use fitted values instead of linear predictor
#' resettest(fit, type = "fitted")
#'
#' @importFrom stats pchisq
#'
#' @seealso \code{\link{simplexreg}}.
#' @export
resettest <- function(model, dispersion = TRUE, power = 2,
                      type = c("lp", "fitted")) {

  if (!inherits(model, "simplexregression"))
    stop("'model' must be an object of class 'simplexregression'")

  if (!is.numeric(power) || any(power < 2) || any(power != floor(power)))
    stop("'power' must be a vector of integers >= 2 (e.g., 2 or 2:3).")

  type <- match.arg(type)

  METHOD <- "RESET test"
  DNAME  <- paste(deparse(model$formula), collapse = "")

  y <- as.vector(model$y)

  # Base for augmented terms: linear predictor or fitted values
  base <- if (type == "lp") model$mu.lp else model$mu.fv

  # Build matrix of augmented terms: base^2, base^3, ...
  aug_cols <- matrix(
    sapply(power, function(pw) base^pw),
    ncol = length(power)
  )
  colnames(aug_cols) <- paste0(if (type == "lp") "lp^" else "mu^", power)

  x <- cbind(model$mu.x, aug_cols)
  z <- if (dispersion) cbind(model$sigma2.x, aug_cols) else model$sigma2.x

  ctrl          <- model$control
  ctrl$hessian  <- FALSE

  modelh1 <- tryCatch(
    simplexreg.fit(y,
                   x[, colnames(x) != "(Intercept)", drop = FALSE],
                   z[, colnames(z) != "(Intercept)", drop = FALSE],
                   link.mu     = model$mu.link,
                   link.sigma2 = model$sigma2.link,
                   x_names     = colnames(x),
                   z_names     = colnames(z),
                   control     = ctrl,
                   weights     = model$weights),
    error = function(e) {
      stop("Failed to fit the augmented model for the RESET test. ",
           "The information matrix is singular. ",
           "This may occur with parametric link functions or near-collinear predictors. ",
           "Try fitting the model with a fixed link function (e.g., link.mu = 'logit').")
    }
  )

  lH0  <- model$loglik
  lH1  <- modelh1$loglik
  LR   <- 2 * (lH1 - lH0)
  df   <- model$df.residual - modelh1$df.residual
  PVAL <- pchisq(LR, df, lower.tail = FALSE)

  names(df) <- "df"
  names(LR) <- "RESET"

  RVAL <- list(statistic = LR, parameter = df, p.value = PVAL,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"

  return(RVAL)
}

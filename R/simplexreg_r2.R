################################################################################
#               SIMPLEX REGRESSION - R2 MEASURES                              #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto            #
# Date: 2026-07-14                                                            #
# Description: Pseudo-R2 measures for simplex regression models               #
################################################################################

#' @title Pseudo-R2 Measures for Simplex Regression
#' @description Extracts the Ferrari and Cribari-Neto pseudo-\eqn{R^2}
#' (\eqn{R^2_{FC}}), the likelihood-ratio pseudo-\eqn{R^2}
#' (\eqn{R^2_N}), the conventional coefficient of determination
#' (\eqn{R^2}), and computes their finite-sample corrections. The
#' corrections for \eqn{R^2_{FC}} and \eqn{R^2_N} follow Bayer and
#' Cribari-Neto (2017), whereas the correction for the conventional
#' \eqn{R^2} follows Hu and Shao (2008).
#'
#' @param ... One or more objects of class \code{"simplexregression"} to be
#' evaluated.
#' @param alpha1 Numeric, \eqn{\alpha_1 \in [0, 1]} controlling the relative penalty
#' weight given to the mean and dispersion submodels in the weighted
#' correction \eqn{R^2_{Nw_c}}. Default is \code{0.4}, as
#' recommended by Bayer and Cribari-Neto (2017).
#' @param alpha2 Numeric, \eqn{\alpha_2 > 0}, controlling the penalization
#' intensity in \eqn{R^2_{Nw_c}}. Default is \code{1}, as recommended by
#' Bayer and Cribari-Neto (2017).
#' @param alpha3 Penalization constant for the correction of
#' \eqn{R^2_{HS}}. One of \code{"1"} (default, \eqn{\alpha_3 = 1}) or
#' \code{"log"} (\eqn{\alpha_3 = \log(n)}).
#'
#' @details
#' \eqn{R^2_{FC}} is the Ferrari and Cribari-Neto (2004) pseudo-R2, defined as
#' the squared correlation between the fitted mean linear predictor and the
#' link-transformed response, i.e.,
#' \deqn{R^2_{FC} = \mathrm{corr}^2\!\left(
#' \boldsymbol{\hat\eta_1}, g(\boldsymbol{y})\right),}
#' where \eqn{\boldsymbol{\hat\eta_1}} denotes the vector of fitted mean linear
#' predictors, \eqn{g(\cdot)} is the mean link function, and \eqn{\boldsymbol{y}}
#' is the vector of observed values of the response variable.
#'
#' \eqn{R^2_N} is a likelihood-ratio-based pseudo-R2 (Nagelkerke, 1991), defined as
#' \deqn{R^2_N = 1 - \left(\frac{L_{null}}{L_{fit}}\right)^{2/n},}
#' where \eqn{L_{null}} and \eqn{L_{fit}} are the maximized likelihoods of the
#' null (intercept-only) and fitted models, respectively.
#'
#' The conventional coefficient of determination is
#' \deqn{R^2 = 1 - \frac{\sum_{i=1}^n (y_i - \hat\mu_i)^2}{\sum_{i=1}^n (y_i - \bar y)^2},}
#' where \eqn{\bar y} is the sample mean of the responses.
#'
#' \strong{Corrections:}
#' Both finite-sample corrections implemented here for
#' \eqn{R^2_{FC}} and \eqn{R^2_N} were proposed by Bayer and Cribari-Neto
#' (2017) for beta regression with varying precision and are used here in
#' their direct simplex regression analogue.
#'
#' The simple correction, a function of the total number of estimated
#' parameters \eqn{r} only,
#' is applied to both \eqn{R^2_{FC}} and \eqn{R^2_N}:
#' \deqn{R^2_{FC_c} = 1 - (1 - R^2_{FC})\frac{n-1}{n-r}, \qquad
#' R^2_{N_c} = 1 - (1 - R^2_{N})\frac{n-1}{n-r}.}
#'
#' A second, weighted correction is defined for \eqn{R^2_N} only, penalizing
#' the mean and dispersion submodels asymmetrically through \code{alpha1} and
#' controlling penalization intensity through \code{alpha2}:
#' \deqn{R^2_{Nw_c} = 1 - (1 - R^2_N)\left(\frac{n-1}{n-(1+\alpha_1)p-
#' (1-\alpha_1)q}\right)^{\alpha2},}
#' where \eqn{p} and \eqn{q} being the number of parameters in the mean and
#' dispersion submodels, respectively.
#' Setting \code{alpha1 = 0} and \code{alpha2 = 1} reduces \eqn{R^2_{Nw_c}} to
#' the simple \eqn{R^2_{N_c}} above.
#' Bayer and Cribari-Neto (2017) recommend
#' \code{alpha1 = 0.4} and \code{alpha2 = 1} as sensible defaults; both arguments
#' are exposed here for users who wish to specify different values.
#'
#' Finally, the conventional coefficient of determination admits the
#' finite-sample correction proposed by Hu and Shao (2008):
#' \deqn{R^2_{HS} = 1 - \frac{n-1}{n - \alpha_3 r}
#' \frac{\sum_{i=1}^n (y_i - \hat\mu_i)^2}{\sum_{i=1}^n (y_i - \bar y)^2}.}
#' where \eqn{\alpha_3} is a penalization constant.
#' When \code{alpha_3 = "1"} (default), \eqn{\alpha_3 = 1}, \eqn{R^2_{HS}}
#' reduces to the modified R2 of Mittlböck and Schemper (2002). When
#' \code{alpha_3 = "log"}, \eqn{\alpha_3 = \log(n)}, as evaluated by Hu and Shao
#' (2008).
#'
#' The measures \eqn{R^2_{FC}} and \eqn{R^2_N} take values in \eqn{[0,1]}.
#' The conventional coefficient of determination \eqn{R^2}, as well as the
#' corrected measures \eqn{R^2_{FC_c}}, \eqn{R^2_{N_c}}, \eqn{R^2_{Nw_c}}, and
#' \eqn{R^2_{HS}}, take values in \eqn{(-\infty,1]}.
#' Larger values indicate better model fit.
#'
#' @return When a single model is provided, a named numeric vector with
#' components \code{R2_FC}, \code{R2_FC_c}, \code{R2_N}, \code{R2_N_c},
#' \code{R2_Nw_c}, \code{R2}, and \code{R2_HS}. When multiple models are
#' provided, a data frame with one row per model.
#'
#' @references
#' Ferrari, S. L. P. and Cribari-Neto, F. (2004). Beta regression for
#' modelling rates and proportions. \emph{Journal of Applied Statistics},
#' \bold{31}(7), 799--815. \doi{10.1080/0266476042000214501}
#'
#' Nagelkerke, N. J. D. (1991). A note on a general definition of the
#' coefficient of determination. \emph{Biometrika}, \bold{78}(3), 691--692.
#' \doi{10.1093/biomet/78.3.691}
#'
#' Hu, B. and Shao, J. (2008). Generalized linear model selection using R2.
#' \emph{Journal of Statistical Planning and Inference}, \bold{138}(12),
#' 3705--3712. \doi{10.1016/j.jspi.2007.12.009}
#'
#' Mittlböck, M. and Schemper, M. (2002). Explained variation for logistic
#' regression -- small sample adjustments, confidence intervals and other
#' issues. \emph{Statistics in Medicine}, \bold{21}(23), 3547--3562.
#' \doi{10.1002/1521-4036(200204)44:3<263::AID-BIMJ263>3.0.CO;2-7}
#'
#' Bayer, F. M. and Cribari-Neto, F. (2017). Model selection criteria in beta
#' regression with varying dispersion. \emph{Communications in Statistics -
#' Simulation and Computation}, \bold{46}(1),
#' 729--746. \doi{10.1080/03610918.2014.977918}
#'
#' @examples
#' data(ReadingSkills, package = "SimplexRegression")
#' fit1 <- simplexreg(accuracy ~ dyslexia * iq | dyslexia + iq + I(iq^2),
#'                    data = ReadingSkills)
#'
#' # Single model, default alpha1, alpha2 and alpha3
#' r2(fit1)
#'
#' # Comparing multiple models
#' fit2 <- simplexreg(accuracy ~ dyslexia * iq | dyslexia + iq + I(iq^2),
#'                    data = ReadingSkills, link.mu = "loglog")
#' r2(fit1, fit2)
#'
#' # Custom alpha1/alpha2 for the weighted correction, and alpha3 = log(n)
#' r2(fit1, fit2, alpha1 = 0.6, alpha2 = 1.5, alpha3 = "log")
#'
#' @seealso \code{\link{simplexreg}}, \code{\link{press}}
#' @export
r2 <- function(..., alpha1 = 0.4, alpha2 = 1, alpha3 = c("1", "log")) {

  models <- list(...)
  M <- length(models)

  if (!all(sapply(models, function(x) inherits(x, "simplexregression")))) {
    stop("All arguments must be objects of class 'simplexregression'.")
  }
  if (alpha1 < 0 || alpha1 > 1) {
    stop("'alpha1' must be between 0 and 1.")
  }
  if (alpha2 <= 0) {
    stop("'alpha2' must be positive.")
  }

  alpha3 <- match.arg(alpha3)

  # Get model names (supports named arguments)
  model_names <- names(models)
  if (is.null(model_names)) {
    model_names <- vapply(as.list(substitute(list(...)))[-1L], deparse, character(1))
  }

  r2_fc <- r2_fc_c <- r2_n <- r2_n_c <- r2_nw_c <- r2 <- r2_hs <- numeric(M)

  for (i in seq_len(M)) {
    model <- models[[i]]

    n <- model$nobs
    p <- ncol(model$mu.x)
    q <- ncol(model$sigma2.x)
    parametric <- !is.na(model$coefficients$lambda)
    k <- if (parametric) p + q + 1 else p + q

    # --- R2 Ferrari-Cribari-Neto: pulled from the fitted object ---
    r2_fc[i] <- model$R2_FC
    r2_fc_c[i] <- 1 - (1 - r2_fc[i]) * ((n - 1) / (n - k))

    # --- R2 likelihood ratio: pulled from the fitted object ---
    r2_n[i] <- model$R2_N

    # Simple correction, function of k only (Bayer and Cribari-Neto, 2017)
    r2_n_c[i] <- 1 - (1 - r2_n[i]) * ((n - 1) / (n - k))

    # Weighted correction, separately penalizing mean (p) and dispersion (q)
    # submodels via alpha1/alpha2 (Bayer and Cribari-Neto, 2017)
    r2_nw_c[i] <- 1 - (1 - r2_n[i]) *
      ((n - 1) / (n - (1 + alpha1) * p - (1 - alpha1) * q))^alpha2

    # --- Coefficient of determination ---
    y    <- as.vector(model$y)
    mu   <- as.vector(model$mu.fv)
    ybar <- mean(y)

    r2[i] <- 1 - (sum((y - mu)^2) / sum((y - ybar)^2))

    # --- R2 Hu and Shao (2008) ---
    alpha3_value <- switch(alpha3, "1" = 1, "log" = log(n))

    r2_hs[i] <- 1 - ((n - 1) / (n - alpha3_value * k)) *
      (sum((y - mu)^2) / sum((y - ybar)^2))
  }

  if (M == 1) {
    return(c(R2_FC = r2_fc, R2_FC_c = r2_fc_c,
             R2_N = r2_n, R2_N_c = r2_n_c, R2_Nw_c = r2_nw_c,
             R2 = r2, R2_HS = r2_hs))
  } else {
    data.frame(R2_FC = r2_fc, R2_FC_c = r2_fc_c,
               R2_N = r2_n, R2_N_c = r2_n_c, R2_Nw_c = r2_nw_c,
               R2 = r2, R2_HS = r2_hs,
               row.names = model_names)
  }
}

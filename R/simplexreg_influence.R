################################################################################
#               SIMPLEX REGRESSION - LOCAL AND GLOBAL INFLUENCE                #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2025-11-08                                                             #
# Description: Influence measures: Local influence via case-weight and response#
#                               perturbation scheme, and generalized leverage  #
################################################################################

# ==============================================================================
# HELPER: Check if model uses parametric mean link function
# ==============================================================================
is_parametric <- function(object) {
  !is.na(object$coefficients$lambda)
}

# ==============================================================================
# 1. LOCAL INFLUENCE
# ==============================================================================

#' @title Local Influence for Simplex Regression Models
#' @description Computes local influence measures under case-weight and
#' response perturbation schemes for simplex regression models with
#' parametric or fixed mean link function.
#'
#' @param model An object of class \code{simplexregression}.
#' @param scheme Character string specifying the perturbation scheme:
#' \code{"case.weight"} or \code{"response"}.
#' @param parameter Character string indicating the parameter block:
#' \code{"theta"} (default, all parameters), \code{"beta"} (mean submodel),
#' or \code{"gamma"} (dispersion submodel).
#' @param type Character string specifying the influence measure to compute:
#' \code{"Ci"} (default, total local influence / normal curvature) or
#' \code{"dmax"} (maximum influence direction).
#' @param plot Logical; if \code{TRUE}, produces an index plot of the selected
#' measure. Default is \code{FALSE}.
#' @param threshold Numeric threshold for identifying influential observations in.
#' If \code{NULL} (default), no observations are highlighted.
#' @param label.pos Position(s) for outlier labels in plot. Can be a single value
#' (applied to all labels) or a vector. Values: 1=below, 2=left, 3=above, 4=right.
#' @param plot.type Character string controlling the plot style when
#' \code{plot = TRUE}. If \code{NULL} (default), uses \code{"h"} (vertical
#' lines). Passed to the \code{type} argument of \code{plot()}.
#' @param ... Additional graphical parameters passed to \code{plot()}.
#'
#' @return If \code{plot = FALSE} (default), a list containing:
#' \itemize{
#'   \item \code{dmax.beta}: Maximum influence direction for mean parameters;
#'   \item \code{dmax.gamma}: Maximum influence direction for dispersion parameters;
#'   \item \code{dmax.theta}: Maximum influence direction for all parameters;
#'   \item \code{Ci.beta}: Total local influence for mean parameters;
#'   \item \code{Ci.gamma}: Total local influence for dispersion parameters;
#'   \item \code{Ci.theta}: Total local influence for all parameters.
#' }
#' If \code{plot = TRUE}, the same list is returned invisibly.
#'
#' @details
#' Measures local influence based on the curvature of the log-likelihood surface
#' under small perturbations. Two perturbation schemes are implemented:
#' \itemize{
#'   \item \strong{Case-weight}: Perturbs observation weights;
#'   \item \strong{Response perturbation}: Perturbs response values.
#' }
#'
#' The index plot of \code{dmax} can be used to detect observations that are
#' jointly influential for parameters. The index plot of the normal curvature
#' \code{Ci} can be used to detect observations that are individually influential
#' for parameters.
#'
#' @examples
#' # Simulate and fit data
#' n <- 100
#' x1 <- runif(n)
#' x2 <- runif(n)
#' mu <- parametric_mean_link_inv(0.8 - 1.2*x1 - 1.5*x2, 0.25, "plogit2")
#' y <- rsimplex(n, mu, 0.5)
#' data <- data.frame(y = y, x1 = x1, x2 = x2)
#' fit <- simplexreg(y ~ x1 + x2 | 1, data = data, link.mu = "plogit2")
#'
#' # Local influence under case-weight perturbation — return results
#' infl_cw <- local.influence(fit, scheme = "case.weight")
#'
#' # Plot Ci for beta directly from the function call
#' local.influence(fit, scheme = "case.weight",
#'                 parameter = "beta", type = "Ci", plot = TRUE)
#'
#' # Plot dmax for all parameters under response perturbation
#' local.influence(fit, scheme = "response",
#'                 parameter = "theta", type = "dmax", plot = TRUE)
#'
#' @references
#' Espinheira, P. L., Silva, A. O. (2020). Residual and influence analysis to a
#' general class of simplex regression. \emph{TEST}, \bold{29}, 523–-552.
#' \doi{10.1007/s11749-019-00665-3}
#'
#' @seealso \code{\link{hatvalues.simplexregression}},
#' \code{\link{cooks.distance.simplexregression}},
#' \code{\link{gleverage.simplexregression}}.
#'
#' @importFrom stats plogis quantile
#' @export
local.influence <- function(model, scheme = c("case.weight", "response"),
                            parameter = c("theta", "beta", "gamma"),
                            type      = c("Ci", "dmax"),
                            plot      = FALSE, threshold = NULL,
                            label.pos = 3, plot.type = NULL, ...) {

  if (!inherits(model, "simplexregression")) {
    stop("'model' must be an object of class 'simplexregression'")
  }

  scheme <- match.arg(scheme)
  parameter <- match.arg(parameter)
  type      <- match.arg(type)

  parametric <- is_parametric(model)

  n <- model$nobs
  y <- model$y
  mu <- model$mu.fv
  sigma2 <- model$sigma2.fv
  link_mu <- model$mu.link
  link_sigma2 <- model$sigma2.link
  eta1 <- model$mu.lp
  eta2 <- model$sigma2.lp
  X <- model$mu.x
  Z <- model$sigma2.x

  p <- ncol(X)
  q <- ncol(Z)
  r <- if (parametric) p + q + 1 else p + q

  diff <- y - mu
  yoneminy <- y * (1 - y)
  muonemu <- mu * (1 - mu)
  dev <- (diff / muonemu)^2 / yoneminy

  Ui <- (dev/muonemu) + (1/(muonemu^3))
  UI <- diag(Ui, ncol=n, nrow = n)
  Vi <- 1/(2*sigma2^2)
  VI <- diag(Vi, ncol=n, nrow = n)
  ai <- -1/(2*sigma2) + dev*Vi
  AI <- diag(ai, ncol=n, nrow = n)
  Mi <- 2 * ai / sigma2
  YDAG <- diag(diff, ncol=n, nrow = n)
  Sy <- diag(sqrt(variance.simplex(mu, sigma2)), n, n)

  sigma2i <- 1/sigma2
  SIG <- diag(sigma2i, ncol=n, nrow = n)
  Dlink.sigma2 <- dispersion_link_inv_deriv1(eta2,  link_sigma2)
  HI <- diag(Dlink.sigma2, ncol=n, nrow = n)

  deriv2_linksigma2 <- switch(link_sigma2,
                              log = -1 / sigma2^2,
                              sqrt = -1 / (4 * sigma2 * sqrt(sigma2)),
                              identity = 0,
                              softplus.inv = -exp(sigma2) / (expm1(sigma2))^2)

  sdag <- ai * deriv2_linksigma2

  deriv2_dev <- - (2 / ((1 - y) * muonemu^3)) *
    (10*y - (3*(y^2 + mu^4) / (y*muonemu)) + (2*(1 + 6*mu^2 - 4*mu) / (1-mu)))
  dev2 <- deriv2_dev / 2

  b_1 <- ( 2 / (y * (1-mu)^3) + (1 - 3*mu) / (mu^2 * (1-mu)^3) + (diff * Ui)) / yoneminy
  B_1 <- diag(b_1, ncol=n, nrow = n)
  b_2 <- (dev + (2*diff / (y * mu * (1-mu)^2))) / yoneminy
  B_2 <- diag(b_2, ncol=n, nrow = n)

  B <- function(Delta,I,M) {
    t(Delta)%*%(I-M)%*%Delta
  }

  if(parametric){
    lambda <- model$lambda.fv
    Dlink.mu <- parametric_mean_link_inv_deriv1(eta1, lambda, link_mu)

    if(link_mu == "plogit2") {
      exp_aval_frac <- plogis(eta1)^(1/lambda)
      log_aval <- log(plogis(eta1))
      rho <- as.vector(- exp_aval_frac * log_aval / (lambda^2))
      deriv2_linkmu <- as.vector((lambda * mu^lambda * (1+lambda) - lambda) /
                                   (mu^2 * (1-mu^lambda)^2))
      VARTHETA <- as.vector(- (1 / (lambda^3 * (1 + exp(eta1)))) * exp_aval_frac *
                              (lambda + log_aval))
      VARSIG <- as.vector(1/lambda^4 * exp_aval_frac * log_aval * (log_aval + 2*lambda))
    } else {
      log_aval <- log(1+exp(eta1))
      rho <- as.vector((-1/(lambda^2)) * ((1 + exp(eta1)) ^ (-1/lambda)) * log_aval)
      deriv2_linkmu <- as.vector(lambda * (1 - (1+lambda) * (1-mu)^lambda) /
                                   ((1-mu)^2 * (1 - (1-mu)^lambda)^2))
      VARTHETA <- as.vector(exp(eta1) * (1 + exp(eta1))^(-1 - 1/lambda) *
                              (log_aval - lambda) / (lambda^3))
      VARSIG <- as.vector((1 + exp(eta1))^(-1/lambda) * log_aval *
                            (2*lambda - log_aval) / lambda^4)
    }

    L1 <- crossprod(X, sigma2i^2 * diff * Ui * Dlink.mu * Dlink.sigma2 * Z)
    L2 <- crossprod(X, sigma2i * (dev2 * Dlink.mu * rho - diff * Ui * VARTHETA))
    L3 <- crossprod(Z, sigma2i^2 * diff * Ui * Dlink.sigma2 * rho)
    L <- rbind(
      cbind(crossprod(X, sigma2i * (dev2 + diff * Ui * Dlink.mu * deriv2_linkmu) *
                        Dlink.mu^2 * X), L1, L2),
      cbind(t(L1), crossprod(Z, (Vi + Mi + sdag * Dlink.sigma2) * Dlink.sigma2^2 * Z), L3),
      cbind(t(L2), t(L3), sum(sigma2i * (dev2 * rho^2 - diff * Ui * VARSIG)))
    )

    Lbetabeta <- L[1:p,1:p]
    Lbetalambda <- L[1:p, r]
    Llambdabeta <- t(Lbetalambda)
    Llambdalambda <- L[r, r]
    Ldeltadelta <- L[(p+1):(p+q),(p+1):(p+q)]
    Ldeltalambda <- L[(p+1):(p+q),r]
    Llambdadelta <- t(Ldeltalambda)

    # For influence on BETA
    sub_deltalambda <- rbind(
      cbind(Ldeltadelta, Ldeltalambda),
      cbind(Llambdadelta, Llambdalambda)
    )

    sub_deltalambda_inv <- solve(sub_deltalambda)
    B1 <- matrix(0, nrow = r, ncol = r)
    B1[(p+1):(r), (p+1):(r)] <- sub_deltalambda_inv

    # For influence on DELTA
    sub_betalambda <- rbind(
      cbind(Lbetabeta, Lbetalambda),
      cbind(Llambdabeta, Llambdalambda)
    )

    sub_betalambda_inv <- solve(sub_betalambda)
    B2 <- matrix(0, nrow = r, ncol = r)
    B2[c(1:p, r), c(1:p, r)] <- sub_betalambda_inv

    # For influence on THETA
    B3 <- matrix(0, nrow = r, ncol = r)

  } else {
    Dlink.mu <- fixed_mean_link_inv_deriv1(eta1, link_mu)
    deriv2_linkmu <- fixed_mean_link_deriv2(mu, link_mu)

    L1 <- crossprod(X, sigma2i^2 * diff * Ui * Dlink.mu * Dlink.sigma2 * Z)
    L <- rbind(
      cbind(crossprod(X, sigma2i * (dev2 + diff * Ui * Dlink.mu * deriv2_linkmu) *
                        Dlink.mu^2 * X), L1),
      cbind(t(L1), crossprod(Z, (Vi + Mi + sdag * Dlink.sigma2) * Dlink.sigma2^2 * Z)))

    Lbetabeta <- L[1:p,1:p]
    Ldeltadelta <- L[(p+1):r,(p+1):r]

    # For influence on BETA
    sub_delta_inv <- solve(Ldeltadelta)
    B1 <- matrix(0, nrow = r, ncol = r)
    B1[(p+1):r, (p+1):r] <- sub_delta_inv

    # For influence on DELTA
    sub_beta_inv <- solve(Lbetabeta)
    B2 <- matrix(0, nrow = r, ncol = r)
    B2[1:p, 1:p] <- sub_beta_inv

    # For influence on THETA
    B3 <- matrix(0, nrow = r, ncol = r)
  }

  TI <- diag(Dlink.mu, ncol=n, nrow = n)

  if(scheme=="case.weight")
  {
    # Case Weight Perturbation Scheme
    Deltamu <- t(X) %*% SIG %*% TI %*% UI %*% YDAG
    Deltasigma2 <- t(Z) %*% HI %*% AI
    if(parametric){
      Deltarho <- t(rho) %*% SIG %*% UI %*% YDAG
      Delta <- rbind(Deltamu, Deltasigma2, Deltarho)
    } else {
      Delta <- rbind(Deltamu, Deltasigma2)
    }
  } else {
    # Response Perturbation Scheme
    Deltamu <- t(X) %*% SIG %*% TI %*% B_1 %*% Sy
    Deltasigma2 <- t(Z) %*% VI %*% HI %*% B_2 %*% Sy
    if(parametric){
      Deltarho <- t(rho) %*% SIG %*% B_1 %*% Sy
      Delta <- rbind(Deltamu, Deltasigma2, Deltarho)
    } else{
      Delta <- rbind(Deltamu, Deltasigma2)
    }
  }

  L_inv <- solve(L)

  compute_block <- function(Bmat) {
    BX <- B(Delta, L_inv, Bmat)
    eig <- eigen(BX, symmetric = TRUE)
    list(
      dmax = abs(eig$vec[, 1]),
      Ci   = 2 * abs(diag(BX))
    )
  }

  res_theta <- compute_block(B3)
  res_beta  <- compute_block(B1)
  res_gamma <- compute_block(B2)

  result <- list(
    dmax.beta  = res_beta$dmax,
    dmax.gamma = res_gamma$dmax,
    dmax.theta = res_theta$dmax,
    Ci.beta    = res_beta$Ci,
    Ci.gamma   = res_gamma$Ci,
    Ci.theta   = res_theta$Ci
  )

  op <- par(no.readonly = TRUE)
  on.exit(par(op))

  if (plot) {
    measure_name <- paste0(type, ".", parameter)
    values <- result[[measure_name]]

    ylab_expr <- if (type == "Ci") expression(C[i]) else expression(d[max])

    par(mar = c(3, 3, 2, 3), oma = c(0.5, 0.5, 0.5, 0.5), mgp = c(1.92, 0.6, 0))

    pt <- if (is.null(plot.type)) "h" else plot.type

    # Use string ylab to avoid expression serialization issues in do.call
    plot_args <- modifyList(
      list(type     = pt,
           xlab     = "Observation index",
           ylab     = ylab_expr,
           cex.lab  = 1.2,
           cex      = 1,
           cex.axis = 0.8,
           ylim     = c(0, max(values, na.rm = TRUE) * 1.05)),
      list(...)
    )

    do.call(graphics::plot, c(list(values), plot_args))

    if (!is.null(threshold)) {
      outliers <- which(values > threshold)
      if (length(outliers) > 0)
        text(outliers, values[outliers], labels = outliers,
             pos = label.pos, cex = 0.8, col = "red")
    }

    invisible(result)
  } else {
    result
  }
}


# ==============================================================================
# 2. GENERALIZED LEVERAGE
# ==============================================================================

#' @title Generalized Leverage Values for Simplex Regression Models
#' @description Compute the generalized leverage values for simplex regression
#' models with parametric or fixed mean link function.
#'
#' @param model An object of class \code{simplexregression}.
#'
#' @return A numeric vector of generalized leverage values.
#'
#' @details
#' \code{gleverage} computing generalized leverage values as suggested by Wei, Hu,
#' and Fung (1998). Generalized leverage extends the concept of hat values to account for both
#' mean and dispersion parameters. High leverage values indicate observations
#' that have potentially large influence on parameter estimates.
#'
#' @examples
#' # Simulate and fit data
#' n <- 100
#' x1 <- runif(n)
#' x2 <- runif(n)
#' mu <- parametric_mean_link_inv(0.8 - 1.2*x1 - 1.5*x2, 0.25, "plogit2")
#' y <- rsimplex(n, mu, 0.5)
#' data <- data.frame(y = y, x1 = x1, x2 = x2)
#' fit <- simplexreg(y ~ x1 + x2 | 1, data = data, link.mu = "plogit2")
#'
#' # Compute generalized leverage
#' glev <- gleverage(fit)
#'
#' # Plot leverage values
#' plot(glev, type = "h", ylab = "Generalized leverage",
#'      xlab = "Observation index")
#' abline(h = 2 * mean(glev), lty = 2, col = "red")
#'
#' @seealso \code{\link{hatvalues.simplexregression}},
#' \code{\link{cooks.distance.simplexregression}},
#' \code{\link{local.influence}}.
#'
#' @references
#' Justino, M. E. C. and Cribari-Neto, F. (2026).
#' Simplex regression with a flexible logit link: Inference and application
#' to cross-country impunity data.
#' \emph{Applied Mathematical Modelling}, \bold{154}, 116713. \doi{10.1016/j.apm.2025.116713}
#'
#' Wei, B. C., Hu, Y. Q., and Fung, W. K. (1998).
#' Generalized Leverage and Its Applications.
#' \emph{Scandinavian Journal of Statistics}, \bold{25}, 25–37.
#'
#' @importFrom stats plogis
#' @export
gleverage <- function(model) {
  UseMethod("gleverage")
}

#' @rdname gleverage
#' @method gleverage simplexregression
#' @export
gleverage.simplexregression <- function(model){

  parametric <- is_parametric(model)

  n <- model$nobs
  y <- model$y
  mu <- model$mu.fv
  sigma2 <- model$sigma2.fv
  link_mu <- model$mu.link
  link_sigma2 <- model$sigma2.link
  eta1 <- model$mu.lp
  eta2 <- model$sigma2.lp
  X <- model$mu.x
  Z <- model$sigma2.x

  p <- ncol(X)
  q <- ncol(Z)
  r <- if (parametric) p + q + 1 else p + q

  diff <- y - mu
  yoneminy <- y * (1 - y)
  muonemu <- mu * (1 - mu)
  dev <- (diff / muonemu)^2 / yoneminy

  Ui <- (dev/muonemu) + (1/(muonemu^3))
  Vi <- 1/(2*sigma2^2)
  ai <- -1/(2*sigma2) + dev*Vi
  Mi <- 2 * ai / sigma2

  ci <- Ui + (diff / (yoneminy * muonemu)) * (dev + 2*diff / (y * mu * (1 - mu)^2))
  cidag <- (1 / yoneminy) * (dev + 2*diff / (y * mu * (1 - mu)^2))

  sigma2i <- 1/sigma2
  Dlink.sigma2 <- dispersion_link_inv_deriv1(eta2,  link_sigma2)

  deriv2_linksigma2 <- switch(link_sigma2,
                              log = -1 / sigma2^2,
                              sqrt = -1 / (4 * sigma2 * sqrt(sigma2)),
                              identity = 0,
                              softplus.inv = -exp(sigma2) / (expm1(sigma2))^2)

  sdag <- ai * deriv2_linksigma2

  deriv2_dev <- - (2 / ((1 - y) * muonemu^3)) *
    (10*y - (3*(y^2 + mu^4) / (y*muonemu)) + (2*(1 + 6*mu^2 - 4*mu) / (1-mu)))
  dev2 <- deriv2_dev / 2

  if(parametric){
    lambda <- model$lambda.fv
    Dlink.mu <- parametric_mean_link_inv_deriv1(eta1, lambda, link_mu)

    if(link_mu == "plogit2") {
      exp_aval_frac <- plogis(eta1)^(1/lambda)
      log_aval <- log(plogis(eta1))
      rho <- as.vector(- exp_aval_frac * log_aval / (lambda^2))
      deriv2_linkmu <- as.vector((lambda * mu^lambda * (1+lambda) - lambda) /
                                   (mu^2 * (1-mu^lambda)^2))
      VARTHETA <- as.vector(- (1 / (lambda^3 * (1 + exp(eta1)))) * exp_aval_frac *
                              (lambda + log_aval))
      VARSIG <- as.vector(1/lambda^4 * exp_aval_frac * log_aval * (log_aval + 2*lambda))
    } else {
      log_aval <- log(1+exp(eta1))
      rho <- as.vector((-1/(lambda^2)) * ((1 + exp(eta1)) ^ (-1/lambda)) * log_aval)
      deriv2_linkmu <- as.vector(lambda * (1 - (1+lambda) * (1-mu)^lambda) /
                                   ((1-mu)^2 * (1 - (1-mu)^lambda)^2))
      VARTHETA <- as.vector(exp(eta1) * (1 + exp(eta1))^(-1 - 1/lambda) *
                              (log_aval - lambda) / (lambda^3))
      VARSIG <- as.vector((1 + exp(eta1))^(-1/lambda) * log_aval *
                            (2*lambda - log_aval) / lambda^4)
    }

    L1 <- crossprod(X, sigma2i^2 * diff * Ui * Dlink.mu * Dlink.sigma2 * Z)
    L2 <- crossprod(X, sigma2i * (dev2 * Dlink.mu * rho - diff * Ui * VARTHETA))
    L3 <- crossprod(Z, sigma2i^2 * diff * Ui * Dlink.sigma2 * rho)
    L <- rbind(
      cbind(crossprod(X, sigma2i * (dev2 + diff * Ui * Dlink.mu * deriv2_linkmu) *
                        Dlink.mu^2 * X), L1, L2),
      cbind(t(L1), crossprod(Z, (Vi + Mi + sdag * Dlink.sigma2) * Dlink.sigma2^2 * Z), L3),
      cbind(t(L2), t(L3), sum(sigma2i * (dev2 * rho^2 - diff * Ui * VARSIG)))
    )

    D <- cbind(Dlink.mu * X, matrix(0, nrow = n, ncol = q), rho)
    Lty <- t(cbind(sigma2i * Dlink.mu * ci * X, Vi * Dlink.sigma2 * cidag * Z,
                   sigma2i * ci * rho))

  } else {
    Dlink.mu <- fixed_mean_link_inv_deriv1(eta1, link_mu)
    deriv2_linkmu <- fixed_mean_link_deriv2(mu, link_mu)

    L1 <- crossprod(X, sigma2i^2 * diff * Ui * Dlink.mu * Dlink.sigma2 * Z)
    L <- rbind(
      cbind(crossprod(X, sigma2i * (dev2 + diff * Ui * Dlink.mu * deriv2_linkmu) *
                        Dlink.mu^2 * X), L1),
      cbind(t(L1), crossprod(Z, (Vi + Mi + sdag * Dlink.sigma2) * Dlink.sigma2^2 * Z))
    )

    D <- cbind(Dlink.mu * X, matrix(0, nrow = n, ncol = q))
    Lty <- t(cbind(sigma2i * Dlink.mu * ci * X, Vi * Dlink.sigma2 * cidag * Z))
  }

  leverage <- D %*% solve(L) %*% Lty
  diag(leverage)
}

# ==============================================================================
# 3. s_{j,i} MEASURES
# ==============================================================================

#' @title Sample Influence Measures for Simplex Regression Models
#' @description Computes leave-one-out sample influence measures \eqn{s_{3,i}}
#' and \eqn{s_{5,i}} for simplex regression models, based on the information-matrix-based
#' criteria measures proposed by Cribari-Neto, Vasconcellos and Santana
#' e Silva (2025).
#'
#' @param model An object of class \code{simplexregression}.
#' @param data The data frame used to fit \code{model}.
#' @param type Character vector specifying measure(s): \code{"s3"}, \code{"s5"},
#' or both (default).
#' @param interval Character string specifying the outlier detection threshold:
#'   \code{"I1"} (default, moderate) or \code{"I2"} (strict).
#' @param parameter Character string indicating the parameter block: \code{"theta"}
#' (default, all parameters), \code{"beta"} (mean submodel), or \code{"gamma"}
#' (dispersion submodel).
#' @param plot Logical; if \code{TRUE}, produces index plots of \eqn{s_{3,i}}
#'   and \eqn{s_{5,i}} with threshold lines and flagged-observation labels.
#'   Default is \code{FALSE}.
#' @param verbose Logical; if \code{TRUE} (default), prints progress during
#'   leave-one-out refitting.
#' @param label.pos Position(s) for outlier labels in plot. Can be a single value
#' (applied to all labels) or a vector. Values: 1=below, 2=left, 3=above, 4=right.
#' @param plot.type Character string controlling the plot style when
#' \code{plot = TRUE}. If \code{NULL} (default), uses \code{"h"} for
#' \eqn{n \le 150} and \code{"p"} for \eqn{n > 150} (automatic). Passed
#' to the \code{type} argument of \code{plot()}.
#' @param ... Additional graphical parameters passed to \code{plot()}.
#'
#' @return If \code{plot = FALSE} (default), a list containing only the requested
#' measures:
#' \describe{
#'   \item{\code{s3_i}}{(if requested) Numeric vector of \eqn{s_{3,i}} values.}
#'   \item{\code{s5_i}}{(if requested) Numeric vector of \eqn{s_{5,i}} values.}
#'   \item{\code{outliers_s3}}{(if requested) Data frame of flagged observations for \eqn{s_{3,i}}.}
#'   \item{\code{outliers_s5}}{(if requested) Data frame of flagged observations for \eqn{s_{5,i}}.}
#'   \item{\code{limits_s3}}{(if requested) Named vector for \eqn{s_{3,i}} thresholds.}
#'   \item{\code{limits_s5}}{(if requested) Named vector for \eqn{s_{5,i}} thresholds.}
#'   \item{\code{interval}}{Interval type used.}
#'   \item{\code{parameter}}{Parameter block used.}
#'   \item{\code{n}}{Number of observations.}
#' }
#'
#' @details
#' For each observation \eqn{i}, the model is refit on the dataset with
#' observation \eqn{i} removed. Let \eqn{\hat\theta} and \eqn{\hat\theta_{(i)}}
#' denote the full and leave-one-out MLEs, and let \eqn{A_{n,(i)}} and
#' \eqn{B_{n,(i)}} be the corresponding information matrices.
#'
#' \strong{Measures computed:}
#' \deqn{s_{3,i} = m_{3,(i)} / m_3}
#' \deqn{s_{5,i} = D_i^{\mathrm{mod}} - D_i^{\mathrm{gen}}}
#'
#' where
#' \deqn{D_i^{\mathrm{gen}} =
#'   (\hat\theta - \hat\theta_{(i)})^\top (-A_{n,(i)})
#'   (\hat\theta - \hat\theta_{(i)})}
#' \deqn{D_i^{\mathrm{mod}} =
#'   (\hat\theta - \hat\theta_{(i)})^\top
#'   \tfrac{1}{2}(-A_{n,(i)} + B_{n,(i)})
#'   (\hat\theta - \hat\theta_{(i)})}
#'
#' and \eqn{m_3 = \|\mathrm{vech}(P_n^{-1} B_n P_n^{-\top} - I)\|_2},
#' with \eqn{-A_n = P_n P_n^\top} (Cholesky factorisation).
#'
#' If \eqn{-A_{n,(i)}} is not positive definite, \code{nearPD} is used to
#' find the nearest positive-definite matrix and a message is printed.
#'
#' \strong{Threshold intervals} use two asymmetric IQR spreads:
#' \eqn{IQR_1 = Q(0.50) - Q(0.125)} (left) and
#' \eqn{IQR_2 = Q(0.875) - Q(0.50)} (right).
#' Limits are \eqn{v - z \cdot IQR_1} (lower) and \eqn{v + z \cdot IQR_2}
#' (upper), with reference value \eqn{v = 1} for \eqn{s_3} and \eqn{v = 0}
#' for \eqn{s_5}:
#' \itemize{
#'   \item \code{I1} (strict): \eqn{z = 2.5} for \eqn{s_3}; \eqn{z = 4.0} for \eqn{s_5}.
#'   \item \code{I2} (moderate): \eqn{z = 5.0} for \eqn{s_3}; \eqn{z = 8.0} for \eqn{s_5}.
#' }
#'
#' @examples
#' n  <- 100
#' x1 <- runif(n); x2 <- runif(n)
#' mu <- parametric_mean_link_inv(0.8 - 1.2*x1 - 1.5*x2, 0.25, "plogit2")
#' y  <- rsimplex(n, mu, 0.5)
#' dat <- data.frame(y = y, x1 = x1, x2 = x2)
#' fit <- simplexreg(y ~ x1 + x2 | 1, data = dat, link.mu = "plogit2")
#'
#' \donttest{
#' # Return numeric results
#' im <- diag.im(fit, data = dat, type = "s3", interval = "I2", parameter = "theta")
#'
#' # Produce index plots directly
#' diag.im(fit, data = dat, type = "s3", interval = "I2", parameter = "theta", plot = TRUE)
#' }
#'
#' @references
#' Cribari-Neto, F.; Vasconcellos, K. L. P.; Santana e Silva, J. J. (2025).
#' New strategies for detecting atypical observations based on the information
#' matrix equality. \emph{Journal of Applied Statistics}, \bold{52}, 2873--2893.
#' \doi{10.1080/02664763.2025.2487914}
#'
#' @seealso \code{\link{local.influence}}, \code{\link{gleverage}},
#' \code{\link{cooks.distance.simplexregression}}.
#'
#' @importFrom stats quantile
#' @importFrom Matrix nearPD
#' @export
diag.im <- function(model, data, type = c("s3", "s5"), interval  = c("I1", "I2"),
                    parameter = c("theta", "beta", "gamma"),
                    plot = FALSE, verbose = TRUE,
                    label.pos = 3, plot.type = NULL, ...) {

  if (!inherits(model, "simplexregression")) {
    stop("'model' must be an object of class 'simplexregression'")
  }

  type      <- match.arg(type, several.ok = TRUE)
  interval  <- match.arg(interval)
  parameter <- match.arg(parameter)

  parametric <- is_parametric(model)

  # Dispersion submodel has a single intercept-only term when sigma2.x has one
  # column whose only unique fitted value is constant — i.e. fixed dispersion.
  # Influence on gamma is undefined in that case.
  dispersion_is_fixed <- ncol(model$sigma2.x) == 1L &&
    length(unique(round(model$sigma2.fv, 10))) == 1L

  if (parameter == "gamma" && dispersion_is_fixed)
    stop(
      "parameter = \"gamma\" is not applicable when the dispersion is fixed ",
      "(intercept-only dispersion submodel with a single estimated value). ",
      "Use parameter = \"beta\" or parameter = \"theta\" instead.",
      call. = FALSE
    )

  # Z multipliers for thresholds
  z_vals <- list(
    I1 = list(s3 = 2.5, s5 = 4.0),
    I2 = list(s3 = 5.0, s5 = 8.0)
  )

  # 1. Global measures on full data
  n   <- model$nobs
  p <- ncol(model$mu.x)
  q <- ncol(model$sigma2.x)
  full_sj  <- compute_m3(model, parameter)
  m3_full  <- full_sj$v_m3
  An_neg <- full_sj$An_neg
  B_n <- full_sj$Bn

  # ------------------------------------------------------------------
  # 2. Leave-one-out loop
  # ------------------------------------------------------------------
  s3_i <- s5_i <- cook_gen <- cook_mod <- numeric(n)
  An_neg_i <- Bn_i <- array(NA, c(ncol(An_neg), ncol(An_neg), n))

  for (i in seq_len(n)) {

    if (verbose && (i %% 10 == 0 || i == 1))
      message(sprintf("  Leave-one-out refit: observation %d / %d", i, n))

    data_i  <- data[-i, , drop = FALSE]

    # Refit model without observation i
    fit_i <- tryCatch(
      simplexreg(
        formula  = model$formula,
        link.mu  = model$mu.link,
        link.sigma2 = model$sigma2.link,
        data = data_i
      ),
      error = function(e) {
        warning(sprintf("Refit failed for observation %d: %s", i, conditionMessage(e)))
        NULL
      }
    )

    if (is.null(fit_i)) {
      s3_i[i] <- s5_i[i] <- NA
      next
    }

    # --- s3_i: ratio of m3 measures ---
    sj_i <- compute_m3(fit_i, parameter = parameter)
    s3_i[i] <- sj_i$v_m3 / m3_full

    An_neg_i[,,i] <- sj_i$An_neg
    Bn_i[,,i] <- sj_i$Bn

    if(parametric){
      coef_full <- unlist(model$coefficients)
      coef_loo <- unlist(fit_i$coefficients)
    } else {
      coef_full <- unlist(model$coefficients[c("mean","dispersion")])
      coef_loo  <- unlist(fit_i$coefficients[c("mean","dispersion")])
    }

    if(parameter == "beta") {
      theta_hat_des <- coef_full[1:p] - coef_loo[1:p]
      theta_hat_des_lin <- t(theta_hat_des)
    } else if(parameter == "gamma") {
      theta_hat_des <- coef_full[(p+1):(p+q)] - coef_loo[(p+1):(p+q)]
      theta_hat_des_lin <- t(theta_hat_des)
    } else {  # theta
      theta_hat_des <- coef_full - coef_loo
      theta_hat_des_lin <- t(theta_hat_des)
    }

    cook_gen[i] <- theta_hat_des_lin %*% An_neg_i[,,i] %*% theta_hat_des
    cook_mod[i] <- theta_hat_des_lin %*% ( 0.5*(An_neg_i[,,i] +  Bn_i[,,i]) ) %*% theta_hat_des

    s5_i[i] <- cook_mod[i] - cook_gen[i]
  }

  # ------------------------------------------------------------------
  # 3. Robust threshold computation
  # ------------------------------------------------------------------
  iqr_limits <- function(x, ref, z) {
    q <- quantile(x, c(0.125, 0.5, 0.875), na.rm = TRUE)
    iqr_left <- q[2] - q[1]
    iqr_right <- q[3] - q[2]
    c(lower = ref - z * iqr_left, upper = ref + z * iqr_right)
  }

  limits <- list(
    s3 = if ("s3" %in% type) iqr_limits(s3_i, 1, z_vals[[interval]]$s3) else NULL,
    s5 = if ("s5" %in% type) iqr_limits(s5_i, 0, z_vals[[interval]]$s5) else NULL
  )

  # Flag outliers
  outliers <- list()
  if ("s3" %in% type) {
    idx <- which(s3_i < limits$s3[1] | s3_i > limits$s3[2])
    outliers$s3 <- if (length(idx)) data.frame(Obs = idx, s3_i = s3_i[idx])
  }
  if ("s5" %in% type) {
    idx <- which(s5_i < limits$s5[1] | s5_i > limits$s5[2])
    outliers$s5 <- if (length(idx)) data.frame(Obs = idx, s5_i = s5_i[idx])
  }

  # ------------------------------------------------------------------
  # 5. Return
  # ------------------------------------------------------------------
  result <- list(
    interval = interval,
    parameter = parameter,
    n = n
  )

  if ("s3" %in% type) {
    result$s3_i <- s3_i
    result$outliers_s3 <- outliers$s3
    result$limits_s3 <- limits$s3
  }

  if ("s5" %in% type) {
    result$s5_i <- s5_i
    result$outliers_s5 <- outliers$s5
    result$limits_s5 <- limits$s5
  }

  # ------------------------------------------------------------------
  # 6. Plot (if requested) or return
  # ------------------------------------------------------------------
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  if (plot) {
    n_plots <- length(type)
    if (n_plots == 2L) par(mfrow = c(1L, 2L))

    pt <- if (is.null(plot.type)) (if (n > 150L) "p" else "h") else plot.type

    make_ylim <- function(vals, lims) {
      rng <- diff(range(vals, lims, na.rm = TRUE))
      pad <- rng * 0.08
      c(min(vals, lims, na.rm = TRUE) - pad * 0.4,
        max(vals, lims, na.rm = TRUE) + pad)
    }

    label_obs <- function(vals, lims, lpos) {
      idx <- which(vals < lims[1] | vals > lims[2])
      for (k in seq_along(idx)) {
        ii  <- idx[k]
        pos <- lpos[((k - 1L) %% length(lpos)) + 1L]
        text(ii, vals[ii], labels = ii, pos = pos,
             cex = 0.8, col = "red", offset = 0.3)
      }
    }

    make_plot <- function(vals, lims, ylab_expr) {
      par(mar = c(3, 3, 2, 3), oma = c(0.5, 0.5, 0.5, 0.5), mgp = c(2, 0.6, 0))
      plot_args <- modifyList(
        list(type = pt, pch = 1, xlab = "Observation index", ylab = ylab_expr,
             cex.lab = 1.2, cex = 1, cex.axis = 0.8,
             ylim = make_ylim(vals, lims)),
        list(...)
      )
      do.call(graphics::plot, c(list(vals), plot_args))
      graphics::abline(h = lims, lty = 2, col = "gray60")
      label_obs(vals, lims, label.pos)
    }

    if ("s3" %in% type)
      make_plot(result$s3_i, result$limits_s3, expression(s[3 * i]))
    if ("s5" %in% type)
      make_plot(result$s5_i, result$limits_s5, expression(s[5 * i]))

    invisible(result)
  } else {
    result
  }
}

#' @keywords internal
compute_m3 <- function(model, parameter = c("theta", "beta", "gamma")) {

  parametric <- is_parametric(model)

  n <- model$nobs
  y <- model$y
  mu <- model$mu.fv
  sigma2 <- model$sigma2.fv

  X <- model$mu.x
  Z <- model$sigma2.x

  p <- ncol(X)
  q <- ncol(Z)

  eta1 <- model$mu.lp
  eta2 <- model$sigma2.lp

  link_mu <- model$mu.link
  link_sigma2 <- model$sigma2.link

  diff <- y - mu
  yoneminy <- y * (1 - y)
  muonemu <- mu * (1 - mu)
  dev <- (diff / muonemu)^2 / yoneminy

  Ui <- (dev/muonemu) + (1/(muonemu^3))
  Vi <- 1/(2*sigma2^2)
  ai <- -1/(2*sigma2) + dev*Vi
  Mi <- 2 * ai / sigma2

  sigma2i <- 1/sigma2
  Dlink.sigma2 <- dispersion_link_inv_deriv1(eta2,  link_sigma2)

  deriv2_linksigma2 <- switch(link_sigma2,
                              log = -1 / sigma2^2,
                              sqrt = -1 / (4 * sigma2 * sqrt(sigma2)),
                              identity = 0,
                              softplus.inv = -exp(sigma2) / (expm1(sigma2))^2)

  sdag <- ai * deriv2_linksigma2

  deriv2_dev <- - (2 / ((1 - y) * muonemu^3)) *
    (10*y - (3*(y^2 + mu^4) / (y*muonemu)) + (2*(1 + 6*mu^2 - 4*mu) / (1-mu)))
  dev2 <- deriv2_dev / 2

  # Parametric vs fixed link mean
  if(parametric){
    lambda <- model$lambda.fv
    Dlink.mu <- parametric_mean_link_inv_deriv1(eta1, lambda, link_mu)

    if(link_mu == "plogit2") {
      exp_aval_frac <- plogis(eta1)^(1/lambda)
      log_aval <- log(plogis(eta1))
      rho <- as.vector(- exp_aval_frac * log_aval / (lambda^2))
      deriv2_linkmu <- as.vector((lambda * mu^lambda * (1+lambda) - lambda) /
                                   (mu^2 * (1-mu^lambda)^2))
      VARTHETA <- as.vector(- (1 / (lambda^3 * (1 + exp(eta1)))) * exp_aval_frac *
                              (lambda + log_aval))
      VARSIG <- as.vector(1/lambda^4 * exp_aval_frac * log_aval * (log_aval + 2*lambda))
    } else {
      log_aval <- log(1+exp(eta1))
      rho <- as.vector((-1/(lambda^2)) * ((1 + exp(eta1)) ^ (-1/lambda)) * log_aval)
      deriv2_linkmu <- as.vector(lambda * (1 - (1+lambda) * (1-mu)^lambda) /
                                   ((1-mu)^2 * (1 - (1-mu)^lambda)^2))
      VARTHETA <- as.vector(exp(eta1) * (1 + exp(eta1))^(-1 - 1/lambda) *
                              (log_aval - lambda) / (lambda^3))
      VARSIG <- as.vector((1 + exp(eta1))^(-1/lambda) * log_aval *
                            (2*lambda - log_aval) / lambda^4)
    }

    L1 <- crossprod(X, sigma2i^2 * diff * Ui * Dlink.mu * Dlink.sigma2 * Z)
    L2 <- crossprod(X, sigma2i * (dev2 * Dlink.mu * rho - diff * Ui * VARTHETA))
    L3 <- crossprod(Z, sigma2i^2 * diff * Ui * Dlink.sigma2 * rho)
    L <- rbind(
      cbind(crossprod(X, sigma2i * (dev2 + diff * Ui * Dlink.mu * deriv2_linkmu) *
                        Dlink.mu^2 * X), L1, L2),
      cbind(t(L1), crossprod(Z, (Vi + Mi + sdag * Dlink.sigma2) * Dlink.sigma2^2 * Z), L3),
      cbind(t(L2), t(L3), sum(sigma2i * (dev2 * rho^2 - diff * Ui * VARSIG)))
    )

    M <- cbind(
      sigma2i * diff * Ui * Dlink.mu * X,
      ai * Dlink.sigma2 * Z,
      sigma2i * Ui * diff * rho
    )

  } else {
    Dlink.mu <- fixed_mean_link_inv_deriv1(eta1, link_mu)
    deriv2_linkmu <- fixed_mean_link_deriv2(mu, link_mu)

    L1 <- crossprod(X, sigma2i^2 * diff * Ui * Dlink.mu * Dlink.sigma2 * Z)
    L <- rbind(
      cbind(crossprod(X, sigma2i * (dev2 + diff * Ui * Dlink.mu * deriv2_linkmu) *
                        Dlink.mu^2 * X), L1),
      cbind(t(L1), crossprod(Z, (Vi + Mi + sdag * Dlink.sigma2) * Dlink.sigma2^2 * Z)))

    M <- cbind(
      sigma2i * diff * Ui * Dlink.mu * X,
      ai * Dlink.sigma2 * Z
    )
  }

  Bn <- crossprod(M) / n
  An <- - L/n
  An_neg <- -An

  L_complete <- L
  Bn_complete <- Bn * n

  # Select sub-block according to 'parameter'
  if(parameter == "theta"){
    Bn <- Bn
    An <- An
    An_neg <- An_neg
    r <- if (parametric) p + q + 1 else p + q
  } else if (parameter == "beta") {
    Bn <- Bn[1:p, 1:p]
    An <- An[1:p, 1:p]
    An_neg <- An_neg[1:p, 1:p]
    r <- p
  } else {
    Bn <- Bn[(p+1):(p+q), (p+1):(p+q)]
    An <- An[(p+1):(p+q), (p+1):(p+q)]
    An_neg <- An_neg[(p+1):(p+q), (p+1):(p+q)]
    r <- q
  }

  eigenvalues <- eigen(An_neg, only.values = TRUE)$values

  if (min(eigenvalues) <= 1e-10) {
    message("The matrix -A_n is not positive definite. Using nearPD approximation.")
    An_neg = as.matrix(Matrix::nearPD(An_neg)$mat)
  }

  Pn <- t(chol(An_neg))
  Pn_inv <- solve(Pn)

  C3n <- Pn_inv %*% Bn %*% t(Pn_inv) - diag(1, nrow = r, ncol = r)
  v_m3 <- norm(C3n[lower.tri(C3n, diag = TRUE)], type = "2")

  return(list(v_m3 = v_m3, An_neg =  An_neg, Bn =  Bn))
}

# ==============================================================================
# 4. HELLINGER AND WASSERSTEIN DISTANCES
# ==============================================================================

#' @title Distance-Based Influence Diagnostics for Simplex Regression
#' @description Computes leave-one-out influence measures based on distributional
#' distances (Wasserstein W1, W2, or Hellinger) for simplex regression models.
#'
#' @param model An object of class \code{simplexregression}.
#' @param data The data frame used to fit \code{model}.
#' @param type Character string or integer specifying the distance measure:
#'   \code{"W1"} or \code{1} (default, Wasserstein with p_W = 1), \code{"W2"} or \code{2}
#'   (Wasserstein with p_W = 2), \code{"H"} or \code{3} (Hellinger).
#' @param plot Logical; if \code{FALSE} (default), returns the numeric vector of
#'   distances. If \code{TRUE}, produces an index plot with the ad hoc threshold.
#' @param verbose Logical; if \code{TRUE} (default), prints progress during
#'   leave-one-out refitting.
#' @param label.pos Position(s) for outlier labels in plot. Can be a single value
#' (applied to all labels) or a vector. Values: 1=below, 2=left, 3=above, 4=right.
#' @param plot.type Character string controlling the plot style when
#' \code{plot = TRUE}. If \code{NULL} (default), uses \code{"h"} for
#' \eqn{n \le 150} and \code{"p"} for \eqn{n > 150} (automatic). Passed
#' to the \code{type} argument of \code{plot()}.
#' @param ... Additional graphical parameters passed to \code{plot()}.
#'
#' @return If \code{plot = FALSE}, a list containing:
#' \describe{
#'   \item{\code{distances}}{Numeric vector of length \eqn{n} with the
#'     leave-one-out distances.}
#'   \item{\code{threshold}}{Named numeric vector with the ad hoc upper threshold.}
#'   \item{\code{outliers}}{Data frame of flagged observations (index and distance).}
#'   \item{\code{type}}{Distance type used.}
#'   \item{\code{n}}{Number of observations.}
#' }
#' If \code{plot = TRUE}, the same list is returned invisibly.
#'
#' @details
#' For each observation \eqn{i}, the model is refit on the dataset with
#' observation \eqn{i} removed. The distance between the full and leave-one-out
#' fitted distributions is then computed pointwise and summed.
#'
#' The ad hoc threshold uses an asymmetric IQR spread:
#' \deqn{\text{threshold} = Q(0.75) + (1 + a) \cdot (Q(0.75) - Q(0.25))},
#' where a is its sample skewness.
#'
#' \strong{Warning:} For \eqn{n > 500}, numerical integration used internally
#' by the distance functions may be slow. Consider using a subset or a
#' faster integration method.
#'
#' @examples
#' n  <- 100
#' x1 <- runif(n); x2 <- runif(n)
#' mu <- parametric_mean_link_inv(0.8 + 1.8*x1 - 1.5*x2, 0.25, "plogit2")
#' y  <- rsimplex(n, mu, 0.5)
#' dat <- data.frame(y = y, x1 = x1, x2 = x2)
#' fit <- simplexreg(y ~ x1 + x2 | 1, data = dat, link.mu = "plogit2")
#'
#' # Produce index plot
#' # diag.distances(fit, data = dat, type = "W1", plot = TRUE)
#'
#' @seealso \code{\link{diag.im}}
#'
#' @importFrom pracma quadgk
#'
#' @export
diag.distances <- function(model, data, type = c("W1", "W2", "H"),
                           plot = FALSE, verbose = TRUE,
                           label.pos = 3, plot.type = NULL, ...) {

  if (!inherits(model, "simplexregression"))
    stop("'model' must be an object of class 'simplexregression'")

  # ---- Normalise 'type' (accept integer shortcuts) -------------------------
  if (is.numeric(type)) {
    type <- switch(as.character(as.integer(type)),
                   "1" = "W1",
                   "2" = "W2",
                   "3" = "H",
                   stop("'type' must be 1 (W1), 2 (W2), or 3 (H)."))
  } else {
    type <- match.arg(type)
  }

  n <- model$nobs

  # ---- Large-sample warning ------------------------------------------------
  if (n > 500)
    warning(
      sprintf(paste0("n = %d > 500: numerical integration used by the distance ",
                     "functions may be slow for large samples."), n),
      call. = FALSE
    )

  # ---- Pre-compute full-model predictions ----------------------------------
  mu_full    <- predict.simplexregression(model, newdata = data, type = "response")
  sigma2_full <- predict.simplexregression(model, newdata = data, type = "dispersion")

  # ---- Distance function selector ------------------------------------------
  dist_fun <- switch(type,
                     W1 = function(mf, sf, mr, sr)
                       wasserstein_simplex(mf, sf, mr, sr, p_W = 1, region = "complete")$sum,
                     W2 = function(mf, sf, mr, sr)
                       wasserstein_simplex(mf, sf, mr, sr, p_W = 2)$sum,
                     H  = function(mf, sf, mr, sr)
                       hellinger_simplex(mf, sf, mr, sr)$sum
  )

  # ---- Leave-one-out loop --------------------------------------------------
  distances <- numeric(n)

  for (i in seq_len(n)) {

    if (verbose && (i == 1 || i %% 10 == 0))
      message(sprintf("  Leave-one-out refit: observation %d / %d", i, n))

    fit_i <- tryCatch(
      simplexreg(
        formula     = model$formula,
        link.mu     = model$mu.link,
        link.sigma2 = model$sigma2.link,
        data        = data[-i, , drop = FALSE]
      ),
      error = function(e) {
        warning(sprintf("Refit failed for observation %d: %s",
                        i, conditionMessage(e)), call. = FALSE)
        NULL
      }
    )

    if (is.null(fit_i)) {
      distances[i] <- NA_real_
      next
    }

    mu_red     <- predict.simplexregression(fit_i, newdata = data, type = "response")
    sigma2_red <- predict.simplexregression(fit_i, newdata = data, type = "dispersion")

    distances[i] <- dist_fun(mu_full, sigma2_full, mu_red, sigma2_red)
  }

  # ---- Ad hoc threshold (asymmetric IQR) -----------------------------------
  ad_hoc_threshold <- function(x) {
    q <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
    iqr_right <- q[2] - q[1]
    unname(q[2] + (1 + moments::skewness(x)) * iqr_right)
  }

  thr <- ad_hoc_threshold(distances)

  idx      <- which(distances > thr)
  outliers <- if (length(idx))
    data.frame(Obs = idx, distance = distances[idx])
  else
    data.frame(Obs = integer(0), distance = numeric(0))

  # ---- Assemble result -----------------------------------------------------
  type_label <- c(W1 = "Wasserstein-1", W2 = "Wasserstein-2", H = "Hellinger")[type]

  result <- list(
    distances = distances,
    threshold = c(upper = thr),
    outliers  = outliers,
    type      = type_label,
    n         = n
  )

  # ---- Plot ----------------------------------------------------------------
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  if (plot) {
    ylab_expr <- switch(type,
                        W1 = "Wasserstein (p_W = 1) distance",
                        W2 = "Wasserstein (p_W = 2) distance",
                        H  = "Hellinger distance"
    )

    pt <- if (is.null(plot.type)) {
      if (n > 150) "p" else "h"
    } else {
      plot.type
    }

    y_max <- if (nrow(outliers) > 0) {
      max(distances, na.rm = TRUE) * 1.12
    } else {
      max(distances, na.rm = TRUE) * 1.05
    }

    par(mar = c(3, 3, 2, 3), oma = c(0.5, 0.5, 0.5, 0.5), mgp = c(2, 0.6, 0))

    plot_args <- modifyList(
      list(type = pt, pch = 1, xlab = "Observation index", ylab = ylab_expr,
           cex = 1, cex.lab = 1.2, cex.axis = 0.8,
           ylim = c(0, y_max)),
      list(...)
    )
    do.call(graphics::plot, c(list(distances), plot_args))
    graphics::abline(h = thr, lty = 2, col = "gray60", lwd = 1.5)

    if (nrow(outliers) > 0L)
      text(outliers$Obs, outliers$distance,
           labels = outliers$Obs, pos = label.pos,
           cex = 0.8, col = "red", offset = 0.3)

    invisible(result)
  } else {
    result
  }
}

#' @keywords internal
hellinger_simplex <- function(mu1, sig1, mu2, sig2) {

  rho_simplex <- function(mu1, sig1, mu2, sig2) {
    integrand <- function(y) {
      sqrt(dsimplex(y, mu1, sig1) * dsimplex(y, mu2, sig2))
    }
    pracma::quadgk(integrand, 0, 1, tol = 1e-10)
  }

  n <- length(mu1)

  rho <- sapply(seq_len(n), function(i) {
    rho_simplex(mu1[i], sig1[i], mu2[i], sig2[i])
  })

  distances <- sqrt(1 - rho)

  return(list(distances = distances, sum = sum(distances)))
}

#' @keywords internal
wasserstein_simplex_p_W1 <- function(mu1, sig1, mu2, sig2, lower = 0, upper = 1) {
  integrand <- function(y) {
    abs(psimplex(y, mu1, sig1) - psimplex(y, mu2, sig2))
  }

  pracma::quadgk(integrand, lower, upper, tol = 1e-10)
}

#' @keywords internal
wasserstein_simplex_general <- function(mu1, sig1, mu2, sig2, p_W = 2) {
  if (abs(p_W - 1) < 1e-10) {
    return(wasserstein_simplex_p_W1(mu1, sig1, mu2, sig2))
  }

  integrand <- function(q) {
    abs(qsimplex(q, mu1, sig1) - qsimplex(q, mu2, sig2))^p_W
  }

  integral_result <- pracma::quadgk(integrand, 0, 1, tol = 1e-10)
  return(integral_result^(1/p_W))
}

#' @keywords internal
wasserstein_simplex <- function(mu1, sig1, mu2, sig2, p_W = 1,
                                region = c("complete", "center", "tails")) {

  if (p_W < 1) stop("p_W must be >= 1 to be a valid distance")

  region <- match.arg(region)
  n <- length(mu1)

  if (p_W != 1 && region != "complete") {
    stop("Region calculation only available for p_W = 1")
  }

  if (p_W == 1 && region == "complete") {
    distances <- sapply(seq_len(n), function(i) {
      wasserstein_simplex_p_W1(mu1[i], sig1[i], mu2[i], sig2[i], lower = 0, upper = 1)
    })
  } else if (p_W == 1 && region != "complete") {
    distances <- sapply(seq_len(n), function(i) {

      q1 <- qsimplex(0.25, mu1[i], sig1[i])
      q3 <- qsimplex(0.75, mu1[i], sig1[i])

      if (region == "center") {
        wasserstein_simplex_p_W1(mu1[i], sig1[i], mu2[i], sig2[i], q1, q3)
      } else {
        d_inf <- wasserstein_simplex_p_W1(mu1[i], sig1[i], mu2[i], sig2[i], 0, q1)
        d_sup <- wasserstein_simplex_p_W1(mu1[i], sig1[i], mu2[i], sig2[i], q3, 1)
        d_inf + d_sup
      }
    })
  } else {
    distances <- sapply(seq_len(n), function(i) {
      wasserstein_simplex_general(mu1[i], sig1[i], mu2[i], sig2[i], p_W = p_W)
    })
  }

  return(list(
    distances = distances,
    sum = sum(distances),
    region = region,
    p_W = p_W
  ))
}

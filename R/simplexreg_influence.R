################################################################################
#               SIMPLEX REGRESSION - LOCAL AND GLOBAL INFLUENCE                #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2025-11-08                                                             #
# Description: Influence measures: Hat values, Cook's distances, and Local     #
#              infuence case-weight and response perturbation scheme           #                                         #
################################################################################

# ==============================================================================
# HELPER: Check if model uses parametric mean link function
# ==============================================================================
is_parametric <- function(object) {
  !is.na(object$coefficients$lambda)
}

# ==============================================================================
# 1. HAT VALUES
# ==============================================================================

#' @title Hat Values for Simplex Regression
#' @description Computes hat (leverage) values for the mean submodel.
#'
#' @param model An object of class \code{"simplexregression"}
#' @param ... Additional arguments (currently not used)
#'
#' @importFrom stats hatvalues
#' @return Numeric vector of hat values
#' @export
hatvalues.simplexregression <- function(model, ...) {
  parametric <- is_parametric(model)

  X <- model$mu.x
  mu <- as.vector(model$mu.fv)
  sigma2 <- as.vector(model$sigma2.fv)
  muonemu <- mu * (1 - mu)

  if(parametric){
    weights <- (1/sigma2) * ((3*sigma2 / muonemu) + ( 1 / (muonemu^3))) *
      (parametric_mean_link_inv_deriv1(model$mu.lp, model$lambda.fv, model$mu.link)^2)
  } else {
    weights <- (1/sigma2) * ((3*sigma2 / muonemu) + ( 1 / (muonemu^3))) *
      (fixed_mean_link_inv_deriv1(model$mu.lp, model$mu.link)^2)
  }

  Xw <- X * sqrt(weights)
  Inv <- chol2inv(chol(t(Xw) %*% Xw))

  Hat <- Xw %*% Inv %*% t(Xw)
  hat <- diag(Hat)

  return(hat)
}

# ==============================================================================
# 2. COOK'S DISTANCE
# ==============================================================================

#' @title Approximate Cook's Distance for Simplex Regression
#' @description Computes approximate Cook's distance based on the likelihood
#' displacement approximation.
#'
#' @param model An object of class \code{"simplexregression"}
#' @param ... Additional arguments (currently not used)
#'
#' @return Numeric vector of Cook's distances
#' @importFrom stats cooks.distance
#'
#' @export
cooks.distance.simplexregression <- function(model, ...) {
  h <- hatvalues.simplexregression(model)
  res <- residuals.simplexregression(model, type = "sweighted2")
  cook <- (res^2) * (h / (1 - h))

  return(cook)
}

# ==============================================================================
# 3. LOCAL INFLUENCE
# ==============================================================================

#' @title Local Influence for Simplex Regression
#' @description Computes local influence measures under case-weight and
#' response perturbation schemes.
#'
#' @param model An object of class \code{"simplexregression"}
#' @param scheme Character string specifying perturbation scheme:
#' "case.weight" or "response"
#'
#' @return List containing:
#' \itemize{
#'   \item \code{dmax.beta}: Maximum influence direction for beta parameters
#'   \item \code{dmax.alpha}: Maximum influence direction for delta parameters
#'   \item \code{dmax.theta}: Maximum influence direction for all parameters
#'   \item \code{Ci.beta}: Total local influence for beta parameters
#'   \item \code{Ci.alpha}: Total local influence for delta parameters
#'   \item \code{Ci.theta}: Total local influence for all parameters
#' }
#'
#' @importFrom stats plogis
#' @export
local.influence_simplexreg <- function(model, scheme = c("case.weight", "response")) {

  scheme <- match.arg(scheme)

  parametric <- is_parametric(model)

  n <- model$nobs
  y <- model$y
  mu <- model$mu.fv
  sigma2 <- model$sigma2.fv
  link_mu <- model$mu.link
  link_sigma2 <- model$sigma2.link
  eta1 <- model$mu.lp
  eta2 <- model$sigma2.lp

  diff <- y - mu
  yoneminy <- y * (1 - y)
  muonemu <- mu * (1 - mu)
  dev <- deviance_simplexreg(y, mu)

  X <- model$mu.x
  Z <- model$sigma2.x

  p <- ncol(X)
  q <- ncol(Z)

  Ui <- (dev/muonemu) + (1/(muonemu^3))
  UI <- diag(Ui, ncol=n, nrow = n)
  Vi <- 1/(2*sigma2^2)
  VI <- diag(Vi, ncol=n, nrow = n)
  ai <- -1/(2*sigma2) + dev*Vi
  AI <- diag(ai, ncol=n, nrow = n)
  Mi <- 2 * ai / sigma2
  YDAG <- diag(diff, ncol=n, nrow = n)
  Sy <- diag(sqrt(variance_simplexreg(mu, sigma2)), n, n)

  sigma2i <- 1/sigma2
  SIG <- diag(sigma2i, ncol=n, nrow = n)
  Dlink.sigma2 <- dispersion_link_inv_deriv1(eta2,  link_sigma2)
  HI <- diag(Dlink.sigma2, ncol=n, nrow = n)

  if(link_sigma2 == "log") {
    deriv2_linksigma2 <- as.vector(- 1 / (sigma2^2))
  } else {
    deriv2_linksigma2 <- as.vector(- 1 / (4 * sigma2 * sqrt(sigma2)))
  }

  sdag <- ai * deriv2_linksigma2

  deriv2_dev <- - (2 / ((1 - y) * muonemu^3)) *
    (10*y - (3*(y^2 + mu^4) / (y*muonemu)) + (2*(1 + 6*mu^2 - 4*mu) / (1-mu)))
  dev2 <- deriv2_dev / 2

  b_1 <- ( 2 / (y * (1-mu)^3) + (1 - 3*mu) / (mu^2 * (1-mu)^3) + (diff * Ui)) / yoneminy
  B_1 <- diag(b_1, ncol=n, nrow = n)
  b_2 <- (dev + (2*diff / (y * mu * (1-mu)^2))) / yoneminy
  B_2 <- diag(b_2, ncol=n, nrow = n)

  B <- function(Delta,I,M) {
    B = (t(Delta)%*%(I-M)%*%Delta)
    return(B)
  }

  if(parametric){
    r <- p + q + 1
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
    r <- p + q
    Dlink.mu <- fixed_mean_link_inv_deriv1(eta1, link_mu)
    deriv2_linkmu <- fixed_mean_link_deriv2(mu, link_mu)

    L1 <- crossprod(X, sigma2i^2 * diff * Ui * Dlink.mu * Dlink.sigma2 * Z)
    L <- rbind(
      cbind(crossprod(X, sigma2i * (dev2 + diff * Ui * Dlink.mu * deriv2_linkmu) * Dlink.mu^2 * X), L1),
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

    # Theta
    BT <- B(Delta, solve(L), B3)
    autovmaxthetaPC <- eigen(BT,symmetric=TRUE)$val[1]
    vetorpcthetaPC <- eigen(BT,symmetric=TRUE)$vec[,1]
    dmaxG.theta <- abs(vetorpcthetaPC)
    vCithetaPC <- 2*abs(diag(BT))
    Cb0 <- vCithetaPC
    Cb.theta <- Cb0

    # Beta
    BM <- B(Delta,solve(L),B1)
    autovmaxbetaPC <- eigen(BM,symmetric=TRUE)$val[1]
    vetorpcbetaPC <- eigen(BM,symmetric=TRUE)$vec[,1]
    dmaxG.beta <- abs(vetorpcbetaPC)
    vCibetaPC <- 2*abs(diag(BM))
    Cb1 <- vCibetaPC
    Cb.beta <- Cb1

    # Delta
    BD <- B(Delta,solve(L),B2)
    autovmaxdeltaPC <- eigen(BD,symmetric=TRUE)$val[1]
    vetordeltaPC <- eigen(BD,symmetric=TRUE)$vec[,1]
    dmaxG.alpha <- abs(vetordeltaPC)
    vCideltaPC <- 2*abs(diag(BD))
    Cb2 <- vCideltaPC
    Cb.alpha <- Cb2

    result <- list(dmax.beta = dmaxG.beta,
                   dmax.alpha = dmaxG.alpha,
                   dmax.theta = dmaxG.theta,
                   Ci.beta = Cb.beta,
                   Ci.alpha = Cb.alpha,
                   Ci.theta = Cb.theta)
    return(result)
  }

  if(scheme=="response") {
    # Response Perturbation Scheme
    Deltamu <- t(X) %*% SIG %*% TI %*% B_1 %*% Sy
    Deltasigma2 <- t(Z) %*% VI %*% HI %*% B_2 %*% Sy
    if(parametric){
      Deltarho <- t(rho) %*% SIG %*% B_1 %*% Sy
      Delta <- rbind(Deltamu, Deltasigma2, Deltarho)
    } else{
      Delta <- rbind(Deltamu, Deltasigma2)
    }

    # Theta
    BT <- B(Delta,solve(L),B3)
    autovmaxthetaPC <- eigen(BT,symmetric=TRUE)$val[1]
    vetorthetaRP <- eigen(BT,symmetric=TRUE)$vec[,1]
    dmaxG.theta <- abs(vetorthetaRP)
    vCithetaRP <- 2*abs(diag(BT))
    Cb0 <- vCithetaRP
    Cb.theta <- Cb0

    # Beta
    BM <- B(Delta,solve(L),B1)
    autovmaxbetaRP <- eigen(BM,symmetric=TRUE)$val[1]
    vetorbetaRP <- eigen(BM,symmetric=TRUE)$vec[,1]
    dmaxG.beta <- abs(vetorbetaRP)
    vCibetaRP <- 2*abs(diag(BM))
    Cb1 <- vCibetaRP
    Cb.beta <- Cb1

    # Delta
    BD <- B(Delta,solve(L),B2)
    autovmaxdeltaRP <- eigen(BD,symmetric=TRUE)$val[1]
    vetordeltaRP <- eigen(BD,symmetric=TRUE)$vec[,1]
    dmaxG.alpha <- abs(vetordeltaRP)
    vCideltaRP <- 2*abs(diag(BD))
    Cb2 <- vCideltaRP
    Cb.alpha <- Cb2


    result <- list(dmax.beta = dmaxG.beta,
                   dmax.alpha = dmaxG.alpha,
                   dmax.theta = dmaxG.theta,
                   Ci.beta = Cb.beta,
                   Ci.alpha = Cb.alpha,
                   Ci.theta = Cb.theta)
    return(result)
  }
}

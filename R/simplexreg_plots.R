################################################################################
#          SIMPLEX REGRESSION - PLOTS RESIDUALS AND SIMULATED ENVELOPES        #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2025-11-08                                                             #
# Description: Visual diagnostic tools: plots residuals, half normal plots and #
#              simulated envelopes plot.
################################################################################

# ==============================================================================
# 1. PLOTs RESIDUALS
# ==============================================================================

#' @title Diagnostic Plots for Simplex Regression
#' @description Produces diagnostic plots for the simplex regression model with
#' parametric or fixed mean link function.
#'
#' @param x An object of class \code{"simplexregression"}
#' @param which Numeric vector indicating which plots to produce (1:7)
#' @param type Character string specifying residual type
#' (see \code{residuals.simplexreg}). Default is \code{"quantile"}.
#' @param ask Logical. If \code{TRUE}, the user is asked before each plot.
#' Default is \code{TRUE} when multiple plots are requested.
#' @param ... Additional graphical parameters
#'
#' @details
#' Available plots:
#' \itemize{
#'   \item 1: Residuals vs observation index
#'   \item 2: Residuals vs fitted values
#'   \item 3: Residuals vs linear predictor
#'   \item 4: Observed vs fitted values
#'   \item 5: Normal Q-Q plot
#'   \item 6: Worm plot
#'   \item 7: Cook's distance
#' }
#'
#' @importFrom stats qqnorm qqline
#' @importFrom graphics par abline legend text
#' @importFrom gamlss wp
#' @importFrom grDevices dev.interactive
#' @export
plot.simplexregression <- function(x, which = 1:7,
                               type = c("quantile", "pearson", "pearson P",
                                        "deviance", "deviance P", "standardized",
                                        "weighted", "variance", "variance P",
                                        "biasvariance", "anscombe", "williams",
                                        "response", "score", "dualscore"),
                               ask = prod(par("mfcol")) < length(which) && dev.interactive(),
                               ...) {
  if(!is.numeric(which) || any(which < 1) || any(which > 7))
    stop("`which' must be in 1:7")

  type <- match.arg(type)

  resid <- residuals.simplexregression(x, type = type)

  n <- length(resid)
  show <- rep(FALSE, 7)
  show[which] <- TRUE
  one.fig <- prod(par("mfcol")) == 1

  # Configure ask mode
  if(ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  }

  par(mar = c(3, 3, 2, 3), oma = c(0.5, 0.5, 0.5, 0.5), mgp = c(2, 0.6, 0))

  # 1. Residuals vs indices
  if(show[1]) {
    plot(1:n, resid,
         xlab = "Observation index", ylab = "Residuals",
         pch = 1, cex = 1, cex.axis = 0.8, cex.lab = 1.2,
         ylim = c(min(resid), max(resid)),...)
    abline(h = c(-3, -2, 0, 2, 3), lty = 2, col = "gray60")
  }

  # 2. Residuals vs fitted values
  if(show[2]) {
    plot(x$fitted.values, resid,
         xlab = "Fitted values", ylab = "Residuals",
         pch = 1, cex = 1, cex.axis = 0.8, cex.lab = 1.2,
         ylim = c(min(resid), max(resid)),...)
    abline(h = c(-3, -2, 0, 2, 3), lty = 2, col = "gray60")
  }

  # 3. Residuals vs linear predictor
  if(show[3]) {
    plot(x$mu.lp, resid,
         xlab = "Estimated mean linear predictor", ylab = "Residuals",
         pch = 1, cex = 1, cex.axis = 0.8, cex.lab = 1.2,
         ylim = c(min(resid), max(resid)),...)
    abline(h = c(-3, -2, 0, 2, 3), lty = 2, col = "gray60")
  }

  # 4. Observed vs fitted
  if(show[4]) {
    y <- x$y
    plot(y, x$fitted.values,
         xlab = "Observed values", ylab = "Fitted values",
         pch = 1, cex = 1, cex.axis = 0.8, cex.lab = 1.2, ...)
    abline(a = 0, b = 1, lty = 1, col = "gray60")
  }

  # 5. Q-Q plot
  if(show[5]) {
    qqnorm(resid, xlab = "Normal quantiles", ylab = "Empirical quantiles",
           main = NULL, cex = 1, cex.axis = 0.8, cex.lab = 1.2)
    qqline(resid, col = "gray60")
  }

  # 6. Worm plot
  if(show[6]) {
    gamlss::wp(resid = resid, main = "Worm plot")
  }

  # 7. Cook's distance
  if(show[7]) {
    cook <- cooks.distance.simplexregression(x)
    plot(cook, type = "h", ylim = c(min(cook), max(cook)),
         xlab = "Observation index", ylab = "Cook's distance",
         pch = 1, cex = 1, cex.axis = 0.8, cex.lab = 1.2,...)

    cook_maior_1 <- which(cook > 1)
    if(length(cook_maior_1) > 0) {
      text(cook_maior_1, cook[cook_maior_1], labels = cook_maior_1,
           pos = c(2,3,3,3,3), cex = 0.8, col = "black")
    }
  }

  invisible()
}

# ==============================================================================
# 2. SIMULATED ENVELOPES
# ==============================================================================

#' @title Simulated Envelope Plot for Simplex Regression
#' @description Produces simulated envelope plots for residual diagnostics for the
#' simplex regression model with parametric or fixed mean link function.
#'
#' @param model An object of class \code{"simplexregression"}
#' @param type Character string specifying residual type
#' (see \code{residuals.simplexregression}). Default is \code{"weighted"}.
#' @param nsim Number of simulations (default: 100)
#' @param seed Random seed (default: 2025)
#' @param level Confidence level (default: 0.95)
#' @param ... Additional graphical parameters
#'
#' @importFrom stats qqnorm quantile median
#' @importFrom graphics par legend
#' @export
envelope.simplexreg <- function(model, type = c("weighted", "quantile", "pearson", "pearson P",
                                                "deviance", "deviance P", "standardized",
                                                "variance", "variance P", "biasvariance",
                                                "anscombe", "williams", "response", "score",
                                                "dualscore"),
                                   nsim = 100, seed = 2025, level = 0.95, ...) {
  set.seed(seed)

  type <- match.arg(type)

  n <- model$nobs
  alpha <- (1 - level)/2
  td <- residuals.simplexregression(model, type)
  X <- model$mu.x
  Z <- model$sigma2.x
  mu <- as.vector(model$mu.fv)
  sigma2 <- as.vector(model$sigma2.fv)

  re <- matrix(0, nrow = n, ncol = nsim)

  i = 1
  while (i <= nsim) {
    y1 <- rsimplex_opt(n, mu, sigma2)

    model1 <- suppressWarnings(
      simplexreg_fit(y1, X[,-1, drop=FALSE], Z[,-1, drop=FALSE],
                        diag = 0, link.mu = model$mu.link, link.sigma2 = model$sigma2.link,
                        x_names = model$x_names, z_names = model$z_names))

    if(model1$optim$convergence == 0) {
      re[,i] <- sort(residuals.simplexregression(model1, type))
      i = i + 1
    }
  }

  e1 <- apply(re, 1, quantile, probs = alpha)
  e2 <- apply(re, 1, quantile, probs = 1 - alpha)
  med <- apply(re, 1, median)
  faixa <- range(td,e1,e2)

  # Points outside envelope
  td_sorted <- sort(td)
  td_order <- order(td)
  outside_bands <- (td_sorted < e1) | (td_sorted > e2)
  c95 <- sum(outside_bands)
  posi <- if(c95 > 0) td_order[outside_bands] else NULL
  prop95 <- round(c95/n * 100, 2)

  # Plot
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mar = c(3, 3, 2, 3), oma = c(0.5, 0.5, 0.5, 0.5), mgp = c(2, 0.6, 0))

  r <- qqnorm(td, xlab = "Normal quantiles", ylab = "Empirical quantiles", ylim = faixa,
              pch = 1, main = "", cex = 1, cex.axis = 0.8, cex.lab = 1.2)$x
  par(new = TRUE)
  qqnorm(e1, axes = FALSE, xlab = "", ylab = "", type = "l", ylim = faixa, lty = 1, main = "")
  par(new = TRUE)
  qqnorm(e2, axes = FALSE, xlab = "", ylab = "", type = "l", ylim = faixa, lty = 1, main = "")
  par(new = TRUE)
  qqnorm(med, axes = FALSE, xlab = "", ylab = "", type = "l", ylim = faixa, lty = 2, main = "")
  legend("topleft",
         legend = c(paste("Points outside the envelope:", c95, "(", prop95, "%)"),
                    paste("Total of points:", n)),
         bty="n", cex = 0.8)
}

#' @title Half-Normal Plot with Simulated Envelope for Simplex Regression
#' @description Produces half-normal plots with simulated envelopes for the simplex
#' regression model with parametric or fixed mean link function.
#'
#' @param model An object of class \code{"simplexregression"}
#' @param type Character string specifying residual type
#' (see \code{residuals.simplexregression}). Default is \code{"weighted"}.
#' @param nsim Number of simulations (default: 100)
#' @param level Confidence level (default: 0.95)
#' @param seed Random seed (default: 2025)
#' @param ... Additional graphical parameters
#'
#' @importFrom stats qnorm quantile median
#' @importFrom graphics matplot points legend
#' @export
hnp.simplexreg = function (model, type = c("weighted", "quantile", "pearson", "pearson P",
                                           "deviance", "deviance P", "standardized",
                                           "variance", "variance P", "biasvariance",
                                           "anscombe", "williams", "response", "score",
                                           "dualscore"),
                              nsim = 100, level = 0.95, seed = 2025, ...) {
  set.seed(seed)

  type <- match.arg(type)

  n <- model$nobs
  alpha <- (1 - level)/2
  td <- residuals.simplexregression(model, type)
  X <- model$mu.x
  Z <- model$sigma2.x
  mu <- as.vector(model$mu.fv)
  sigma2 <- as.vector(model$sigma2.fv)

  re <- matrix(0, nrow = n, ncol = nsim)
  e1 <- numeric(n)
  e2 <- numeric(n)
  cOut = 0

  i = 1
  while (i <= nsim) {
    ysim = rsimplex_opt(n, mu, sigma2)

    fit <- suppressWarnings(
      simplexreg_fit(ysim, X[,-1, drop=FALSE], Z[,-1, drop=FALSE], diag = 0,
                        link.mu = model$mu.link, link.sigma2 = model$sigma2.link,
                        x_names = model$x_names, z_names = model$z_names))

    if(fit$optim$convergence == 0){
      re[, i] <- sort(abs(residuals.simplexregression(fit, type)))
      i = i + 1}
  }

  for (i in 1:n) {
    eo = sort(re[i, ])
    e1[i] = quantile(eo, alpha)
    e2[i] = quantile(eo, 1 - alpha)
  }

  e0 = apply(re, 1, median)

  qq = qnorm((n + 1:n + 0.5)/(2 * n + 1.125))
  xx = cbind(qq, qq)
  yy = cbind(e1, e2)

  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mar = c(3, 3, 2, 3), oma = c(0.5, 0.5, 0.5, 0.5), mgp = c(2, 0.6, 0))

  matplot(xx, yy, type = "l", lty = c(1,1), col = c("black", "black"),
          xlab = "Normal quantiles", ylab = "Empirical quantiles",
          cex = 1, cex.axis = 0.8, cex.lab = 1.2, ...)
  res.sorted.abs = sort(abs(td))
  points(qq, res.sorted.abs)

  outside_bands <- (res.sorted.abs < e1) | (res.sorted.abs > e2)
  cOut <- sum(outside_bands)

  prop95 <- round(cOut /n*100, 2)
  legend("topleft",
         legend = c(paste("Points outside the envelope:", cOut, "(", prop95, "%)"),
                    paste("Total of points:", n)), bty="n", cex = 0.8)
}

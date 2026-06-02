################################################################################
#          SIMPLEX REGRESSION - DIAGNOSTIC PLOTS AND SIMULATED ENVELOPES       #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2026-05-26                                                             #
# Description: Visual diagnostic tools: residuals plots, half-normal plots and #
#              simulated envelopes.                                            #
################################################################################

# ==============================================================================
# 1. PLOTs RESIDUALS
# ==============================================================================

#' @title Diagnostic Plots for Simplex Regression Models
#' @description Produces diagnostic plots for fitted simplex regression models with
#' parametric or fixed mean link function.
#'
#' @param x An object of class \code{simplexregression}.
#' @param which Numeric vector indicating which plots to display (\code{1:8}).
#' @param type Character string specifying the residual type (default: \code{"quantile"}).
#' See \code{\link{residuals.simplexregression}} for available options.
#' @param ask Logical; if \code{TRUE}, the user is asked before each plot.
#' Default is \code{TRUE} when multiple plots are requested.
#' @param reset.par Logical; if \code{TRUE}, resets graphical parameters before plotting.
#' Set to \code{FALSE} to preserve user-defined \code{par()} settings such as \code{mfrow}.
#' Default is \code{TRUE}.
#' @param threshold Numeric threshold for identifying influential observations in.
#' If \code{NULL} (default), no observations are highlighted.
#' @param label.pos Position(s) for outlier labels in plots 7 and 8. Can be a single value
#' (applied to all labels) or a vector. Values: 1=below, 2=left, 3=above, 4=right.
#' Default is \code{3} (above). See \code{\link[graphics]{text}} for details.
#' @param plot.type Controls the plot symbol/type for scatter plots and index plots.
#' If \code{NULL} (default), uses \code{pch = 1} (open circles) for residual plots
#' (which = 1--6) and \code{type = "h"} (vertical lines) for Cook's distance and
#' generalized leverage plots (which = 7--8). Otherwise, the value is passed
#' directly to \code{pch} (for scatter plots) or \code{type} (for index plots).
#' @param ... Additional graphical parameters.
#'
#' @details
#' Eight diagnostic plots are available:
#' \itemize{
#'   \item Residuals vs observation index (\code{which = 1}): Identifies outliers
#'   and temporal patterns;
#'
#'   \item Residuals vs fitted values (\code{which = 2}): Checks for heteroscedasticity
#'   and patterns;
#'
#'   \item Residuals vs linear predictor (\code{which = 3}): Evaluates link function
#'   adequacy;
#'
#'   \item Observed vs fitted values (\code{which = 4}): Assesses overall model fit;
#'
#'   \item Normal Q–Q plot (\code{which = 5}): Evaluates residual normality (especially
#'   useful for quantile residuals);
#'
#'   \item Worm plot (\code{which = 6}): Implemented via \code{gamlss::wp};
#'
#'   \item Cook’s distance vs indices of observations (\code{which = 7}): Identifies
#'   influential observations;
#'
#'   \item Generalized leverage vs indices of observations (\code{which = 8}): Identifies
#'   influential observations.
#' }
#'
#' @examples
#' # Simulate data
#' n <- 100
#' x1 <- runif(n, 0, 1)
#' x2 <- runif(n, 0, 1)
#' mu <- parametric_mean_link_inv(0.8 - 1.2*x1 - 1.5*x2, 0.25, "plogit2")
#' y <- rsimplex(n, mu, 0.5)
#' data <- data.frame(y = y, x1 = x1, x2 = x2)
#'
#' # Fit model with parametric mean link functions
#' fit <- simplexreg(y ~ x1 + x2 | 1, data = data, link.mu = "plogit2")
#'
#' # Display all diagnostic plots
#' par(mfrow = c(3, 2))
#' plot(fit, which = 1:8)
#'
#' @importFrom stats qqnorm qqline residuals cooks.distance
#' @importFrom graphics par abline text
#' @importFrom gamlss wp
#' @importFrom grDevices dev.interactive
#'
#' @seealso \code{\link{residuals.simplexregression}},
#' \code{\link{cooks.distance.simplexregression}},
#' \code{\link{halfnormal.plot}}.
#'
#' @importFrom graphics lines
#' @importFrom utils modifyList
#' @export
plot.simplexregression <- function(x, which = 1:8,
                                   type = c("quantile", "pearson",
                                            "deviance", "standardized",
                                            "weighted", "variance",
                                            "biasvariance", "score",
                                            "dualscore", "response"),
                                   ask = prod(par("mfcol")) < length(which) && dev.interactive(),
                                   reset.par = TRUE, threshold = NULL,
                                   label.pos = 3, plot.type = NULL, ...) {

  if (!is.numeric(which) || any(which < 1) || any(which > 8))
    stop("`which' must be in 1:8")

  type <- match.arg(type)

  resid <- residuals(x, type = type)

  n <- length(resid)
  show <- rep(FALSE, 8)
  show[which] <- TRUE
  one.fig <- prod(par("mfcol")) == 1

  # Automatic plot.type defaults
  pt_scatter <- if (is.null(plot.type)) 1    else plot.type   # pch for plots 1-6
  pt_index   <- if (is.null(plot.type)) "h"  else plot.type   # type for plots 7-8

  # Configure ask mode
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  if (ask) par(ask = TRUE)
  if (reset.par) par(mar = c(3,3,2,3), oma = c(0.5,0.5,0.5,0.5), mgp = c(2,0.6,0))

  # Helper to merge defaults with user ...
  plot_args <- function(..., defaults) modifyList(defaults, list(...))
  
  scatter_defaults <- list(pch = pt_scatter, cex = 1, cex.axis = 0.8, cex.lab = 1.2)
  index_defaults   <- list(type = pt_index,  cex = 1, cex.axis = 0.8, cex.lab = 1.2)
  user_args <- list(...)
  
  # 1. Residuals vs indices
  if (show[1]) {
    args <- modifyList(c(list(x = 1:n, y = resid,
                              xlab = "Observation index", ylab = "Residuals"),
                         scatter_defaults), user_args)
    do.call(plot, args)
    abline(h = c(-3, -2, 0, 2, 3), lty = 2, col = "gray60")
  }
  
  # 2. Residuals vs fitted values
  if (show[2]) {
    args <- modifyList(c(list(x = x$fitted.values, y = resid,
                              xlab = "Fitted values", ylab = "Residuals"),
                         scatter_defaults), user_args)
    do.call(plot, args)
    abline(h = c(-3, -2, 0, 2, 3), lty = 2, col = "gray60")
  }

  # 3. Residuals vs linear predictor
  if (show[3]) {
    args <- modifyList(c(list(x = x$mu.lp, y = resid,
                              xlab = "Linear predictor", ylab = "Residuals"),
                         scatter_defaults), user_args)
    do.call(plot, args)
    abline(h = c(-3, -2, 0, 2, 3), lty = 2, col = "gray60")
  }
  
  # 4. Observed vs fitted
  if (show[4]) {
    args <- modifyList(c(list(x = x$y, y = x$fitted.values,
                              xlab = "Observed values", ylab = "Fitted values"),
                         scatter_defaults), user_args)
    do.call(plot, args)
    abline(a = 0, b = 1, lty = 1, col = "gray60")
  }

  # 5. Q-Q plot (qqnorm does not accept x/y, filter valid args)
  if (show[5]) {
    qq_valid <- c("pch", "cex", "cex.axis", "cex.lab", "col", "lwd")
    args <- modifyList(list(cex = 1, cex.axis = 0.8, cex.lab = 1.2,
                            xlab = "Normal quantiles", ylab = "Empirical quantiles",
                            main = NULL),
                       user_args[names(user_args) %in% qq_valid])
    do.call(qqnorm, c(list(y = resid), args))
    qqline(resid, col = "gray60")
  }

  # 6. Worm plot
  if(show[6]) {
    old_par <- par(cex.axis = 0.8)
    gamlss::wp(resid = resid, main = "Worm plot", cex = 1, cex.lab = 1.2, pch = 1)
    par(old_par)
  }

  # 7. Cook's distance
  if (show[7]) {
    cook <- cooks.distance(x, type = "pearson")
    args <- modifyList(c(list(x = cook,
                              xlab = "Observation index", ylab = "Cook's distance",
                              ylim = c(min(cook), max(cook, na.rm = TRUE) * 1.05)),
                         index_defaults), user_args)
    do.call(plot, args)
    if (!is.null(threshold)) {
      idx <- which(cook > threshold)
      if (length(idx) > 0)
        text(idx, cook[idx], labels = idx, pos = label.pos, cex = 0.8, col = "red")
    }
  }

  # 8. Generalized leverage
  if (show[8]) {
    glev <- gleverage(x)
    args <- modifyList(c(list(x = glev,
                              xlab = "Observation index", ylab = "Generalized leverage",
                              ylim = c(min(glev), max(glev, na.rm = TRUE) * 1.05)),
                         index_defaults), user_args)
    do.call(plot, args)
    if (!is.null(threshold)) {
      idx <- which(glev > threshold)
      if (length(idx) > 0)
        text(idx, glev[idx], labels = idx, pos = label.pos, cex = 0.8, col = "red")
    }
  }

  invisible()
}

# ==============================================================================
# 2. SIMULATED ENVELOPE - HALF-NORMAL PLOT
# ==============================================================================

#' @title Half-Normal Plots with Simulated Envelopes for Simplex Regression
#' @description Produces half-normal plots with simulated envelopes for simplex
#' regression model with parametric or fixed mean link function.
#'
#' @param model An object of class \code{simplexregression}.
#' @param type Character string specifying the residual type (default: \code{"weighted"}).
#' See \code{\link{residuals.simplexregression}} for available options.
#' @param nsim Number of simulations for envelope construction (default: 100).
#' @param seed Integer setting the random seed for reproducibility (default: 1987).
#' @param level Confidence level for envelope bounds (default: 0.95).
#' @param ... Additional graphical parameters.
#'
#' @details
#' The envelope is based on the following steps:
#' \enumerate{
#'   \item Simulate \code{nsim} response vectors from the fitted model
#'         using its estimated mean and dispersion parameter;
#'   \item Refitting the model to each simulated dataset;
#'   \item Computing absolute residuals and their order statistics;
#'   \item Obtaining envelope bounds from empirical quantiles.
#' }
#'
#' Points outside the envelope may indicate model inadequacy.
#'
#' @examples
#' # Simulate data (quick)
#' n <- 50  
#' x1 <- runif(n, 0, 1)
#' x2 <- runif(n, 0, 1)
#' mu <- parametric_mean_link_inv(0.8 - 1.2*x1 - 1.5*x2, 0.25, "plogit2")
#' y <- rsimplex(n, mu, 0.5)
#' data <- data.frame(y = y, x1 = x1, x2 = x2)
#'
#' # Fit model
#' fit <- simplexreg(y ~ x1 + x2 | 1, data = data, link.mu = "plogit2")
#'
#' \donttest{
#' # Half-normal plot with fewer simulations (faster)
#' halfnormal.plot(fit, type = "weighted", nsim = 20)
#' }
#'
#' @importFrom stats qnorm quantile median residuals
#' @importFrom graphics matplot points legend
#'
#' @seealso \code{\link{residuals.simplexregression}},
#' \code{\link{plot.simplexregression}}.
#'
#' @importFrom graphics lines
#' @importFrom utils modifyList
#' @export
halfnormal.plot <- function (model, type = c("weighted", "quantile",
                                             "pearson", "deviance",
                                             "standardized", "variance",
                                             "biasvariance", "score",
                                             "dualscore", "response"),
                             nsim = 100, level = 0.95, seed = 1987, ...) {

  if (!inherits(model, "simplexregression")) {
    stop("'model' must be an object of class 'simplexregression'")
  }

  set.seed(seed)
  type <- match.arg(type)

  n <- model$nobs
  alpha <- (1 - level)/2
  
  # Original residuals
  td <- residuals(model, type)
  
  X <- model$mu.x
  Z <- model$sigma2.x
  mu <- as.vector(model$mu.fv)
  sigma2 <- as.vector(model$sigma2.fv)

  re <- matrix(0, nrow = n, ncol = nsim)

  e1 <- numeric(n)
  e2 <- numeric(n)

  ctrl <- model$control
  ctrl$hessian <- FALSE
  ctrl$fsmaxit <- 0  # skip Fisher scoring in simulations for speed
  
  max_tries <- nsim * 5
  tries <- 0
  
  i <- 1
  while (i <= nsim) {
    tries <- tries + 1
    if (tries > max_tries)
      stop("Could not obtain ", nsim, " converged fits after ", max_tries, " attempts.")
    
    ysim <- rsimplex(n, mu, sigma2)

    fit <- suppressWarnings(
        simplexreg.fit(
          y = ysim,
          x = X[, colnames(X) != "(Intercept)", drop = FALSE],
          z = Z[, colnames(Z) != "(Intercept)", drop = FALSE],
          link.mu = model$mu.link,
          link.sigma2 = model$sigma2.link,
          x_names = model$x_names,
          z_names = model$z_names,
          control = ctrl,
          weights = model$weights
        )
      )

    if(fit$optim$convergence == 0){
      re[, i] <- sort(abs(residuals(fit, type)))
      i <- i + 1
      }
  }

  for (j in 1:n) {
    eo <- sort(re[j, ])
    e1[j] <- quantile(eo, alpha)
    e2[j] <- quantile(eo, 1 - alpha)
  }

  e0 <- apply(re, 1, median, na.rm = TRUE)

  qq <- qnorm((n + 1:n + 0.5)/(2 * n + 1.125))
  xx <- cbind(qq, qq)
  yy <- cbind(e1, e2)

  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mar = c(3, 3, 2, 3), oma = c(0.5, 0.5, 0.5, 0.5), mgp = c(2, 0.6, 0))

  # Plot defaults (can be overridden via ...)
  plot_args <- modifyList(
    list(type = "l", lty = c(1, 1), col = c("black", "black"),
         xlab = "Normal quantiles", ylab = "Empirical quantiles",
         cex = 1, cex.axis = 0.8, cex.lab = 1.2),
    list(...)
  )
  
  do.call(matplot, c(list(xx, yy), plot_args))
  
  res.sorted.abs <- sort(abs(td))

  # arguments valid for points() only
  points_valid <- c("pch", "cex", "col", "bg", "lwd")
  point_args <- modifyList(
    list(pch = 1, cex = plot_args$cex),
    list(...)[names(list(...)) %in% points_valid]
  )
  do.call(points, c(list(qq, res.sorted.abs), point_args))
  
  lines(qq, e0, lty = 2, col = "black")
  
  outside_bands <- (res.sorted.abs < e1) | (res.sorted.abs > e2)
  cOut <- sum(outside_bands)
  prop95 <- round(cOut /n*100, 2)
  
  legend("topleft",
         legend = c(paste("Points outside:", cOut, "(", prop95, "%)"),
                    paste("Total points:", n)), 
         bty="n", cex = 0.8)
}

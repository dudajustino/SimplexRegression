################################################################################
#                     SIMPLEX REGRESSION - S3 METHODS                          #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2026-05-26                                                             #
# Description: S3 methods for simplexregression objects including summary,     #
#              print, predict, coef, vcov, logLik, AIC, BIC, and others        #
################################################################################

# ==============================================================================
# HELPERS
# ==============================================================================
is_parametric <- function(object) {
  !is.na(object$coefficients$lambda)
}

.npar <- function(object) {
  length(object$coefficients$mean) +
    length(object$coefficients$dispersion) +
    if(is_parametric(object)) 1L else 0L
}

#' @title Methods for simplexregression Objects
#' @description
#' Methods for extracting information from fitted simplex regression model objects
#' of class \code{simplexregression}.
#'
#' @param object,x An object of class \code{simplexregression}.
#' @param model Character specifying for which component of the model
#' coefficients/covariance should be extracted.
#' @param digits Number of digits to printing.
#' @param newdata Optional data frame for prediction.
#' @param formula A model formula or terms object.
#' @param type Character indicating type of predictions: fitted means of the
#' response (default, "response"), corresponding linear
#' predictor ("link") or fitted dispersion parameter ("dispersion").
#' @param formula. Changes to the formula.
#' @param evaluate If true evaluate the new call else return the call.
#' @param nsim number of response vectors to simulate. Defaults to 1.
#' @param seed an object specifying if and how the random number generator
#' should be initialized.
#' @param k weight of the penalty term in AIC. Default is 2.
#' @param vcov. a specification of the covariance matrix of the estimated coefficients.
#' @param df the degrees of freedom to be used.
#' @param ... Additional arguments.
#'
#' @name simplexreg.methods
#'
#' @return The return value depends on the method called:
#' \itemize{
#'   \item \code{print}: returns \code{x} invisibly;
#'   \item \code{summary}: returns an object of class
#'     \code{"summary.simplexregression"};
#'   \item \code{print.summary}: returns \code{x} invisibly;
#'   \item \code{coef}: named numeric vector of coefficients;
#'   \item \code{vcov}: numeric matrix (variance-covariance);
#'   \item \code{logLik}: object of class \code{"logLik"};
#'   \item \code{fitted}: numeric vector of fitted mean values;
#'   \item \code{predict}: numeric vector or list depending on \code{type};
#'   \item \code{nobs}: integer scalar (number of observations);
#'   \item \code{df.residual}: integer scalar (residual degrees of freedom);
#'   \item \code{deviance}: numeric scalar (total deviance);
#'   \item \code{formula}: the model formula;
#'   \item \code{terms}: the model terms for the selected submodel;
#'   \item \code{model.frame}: a data frame;
#'   \item \code{model.matrix}: numeric matrix of regressors;
#'   \item \code{update}: fitted model or call depending on \code{evaluate};
#'   \item \code{simulate}: a data frame with \code{nsim} columns;
#'   \item \code{AIC}, \code{BIC}, \code{HQIC}: numeric scalar (single model)
#'     or data frame (multiple models);
#'   \item \code{hatvalues}: numeric vector of leverage values;
#'   \item \code{cooks.distance}: numeric vector of Cook's distances;
#'   \item \code{bread}: numeric matrix;
#'   \item \code{estfun}: numeric matrix of score contributions;
#'   \item \code{coeftest}: object of class \code{"coeftest"};
#'   \item \code{lrtest}: object of class \code{"anova"}.
#' }
#'
#' @seealso \code{\link{simplexreg}}.
#'
#' @examples
#' data(ReadingSkills, package = "SimplexRegression")
#' fit <- simplexreg(accuracy ~ dyslexia * iq | dyslexia + iq + I(iq^2),
#'                  data = ReadingSkills)
#'
#' # Extract information
#' summary(fit)
#' coef(fit)
#' vcov(fit)
#' logLik(fit)
#' fitted(fit)
#' AIC(fit)
#' BIC(fit)
#' HQIC(fit)
#' hatvalues(fit)
#' cooks.distance(fit)
#'
#' @aliases print.simplexregression summary.simplexregression coef.simplexregression
#' @aliases vcov.simplexregression logLik.simplexregression fitted.simplexregression
#' @aliases predict.simplexregression nobs.simplexregression df.residual.simplexregression
#' @aliases deviance.simplexregression formula.simplexregression terms.simplexregression
#' @aliases model.frame.simplexregression model.matrix.simplexregression update.simplexregression
#' @aliases simulate.simplexregression AIC.simplexregression BIC.simplexregression
#' @aliases hatvalues.simplexregression cooks.distance.simplexregression bread.simplexregression
#' @aliases estfun.simplexregression coeftest.simplexregression lrtest.simplexregression

NULL

#' @rdname simplexreg.methods
#' @export
print.simplexregression <- function(x, digits = max(3, getOption("digits") - 3),
                                    ...) {

  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)),
      "", sep = "\n")

  converged <- is.null(x$optim$convergence) || x$optim$convergence == 0

  if(!converged) {
    cat("model did not converge\n")
  } else {
    # Coefficients (Mean model)
    if(length(x$coefficients$mean)) {
      cat(sprintf("Coefficients (mean model with %s link):\n", x$mu.link))
      coef_mean <- x$coefficients$mean
      if (!is.null(x$x_names)) names(coef_mean) <- x$x_names
      print.default(format(coef_mean, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else {
      cat("No coefficients (in mean model)\n\n")
    }

    # Coefficients (Dispersion model)
    if(length(x$coefficients$dispersion)) {
      cat(sprintf("Coefficients (dispersion model with %s link):\n",
                  x$sigma2.link))
      coef_disp <- x$coefficients$dispersion
      if (!is.null(x$z_names)) names(coef_disp) <- x$z_names
      print.default(format(coef_disp, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else {
      cat("No coefficients (in dispersion model)\n\n")
    }

    if (is_parametric(x)) {
      cat(sprintf("Link function parameter (parametric %s model)\nlambda: %s\n\n",
                  x$mu.link,
                  format(x$coefficients$lambda, digits = digits)))
    }
  }

  invisible(x)
}

#' @rdname simplexreg.methods
#' @importFrom stats pnorm na.omit
#' @export
summary.simplexregression <- function(object, ...) {

  parametric <- is_parametric(object)

  P2 <- press(object)[1]

  # Extract coefficients
  coef_mean <- object$coefficients$mean
  coef_disp <- object$coefficients$dispersion
  coef_lambda <- object$coefficients$lambda

  # Get standard errors from variance-covariance matrix
  vcov_matrix <- object$vcov
  se <- sqrt(diag(vcov_matrix))

  # Dimensions
  p <- length(coef_mean)
  q <- length(coef_disp)

  # Standard errors
  se_mean <- se[1:p]
  se_disp <- se[(p+1):(p+q)]

  cf <- c(coef_mean, coef_disp)
  if(parametric) cf <- c(cf, coef_lambda)

  se_all <- se[seq_len(length(cf))]
  cf_table <- cbind(cf, se_all, cf/se_all, 2 * pnorm(-abs(cf/se_all)))
  colnames(cf_table) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

  coef_table_mean <- cf_table[1:p, , drop = FALSE]
  coef_table_disp <- cf_table[(p+1):(p+q), , drop = FALSE]

  # Use regressor names if available
  if (!is.null(object$x_names)) {
    rownames(coef_table_mean) <- object$x_names
  } else {
    rownames(coef_table_mean) <- names(coef_mean)
  }

  if (!is.null(object$z_names)) {
    rownames(coef_table_disp) <- object$z_names
  } else {
    rownames(coef_table_disp) <- names(coef_disp)
  }

  # Lambda parameter (if parametric)
  if (parametric) {
    lambda_table <- cf_table[p+q+1, , drop = FALSE]
    rownames(lambda_table) <- "lambda"
  } else {
    lambda_table <- NULL
  }

  mytail <- function(x) x[length(x)]
  if (!is.null(object$optim$counts)) {
    object$iterations <- c("optim" = as.vector(mytail(na.omit(object$optim$counts))))
  } else {
    object$iterations <- NA
  }

  converged <- is.null(object$optim$convergence) || (object$optim$convergence == 0)

  object$fitted.values <- object$terms <- object$model <- object$y <-
    object$x <- object$z <- object$levels <- object$contrasts <-
    object$start <- NULL

  result <- list(
    call = object$call,
    coefficients = list(
      mean = coef_table_mean,
      dispersion = coef_table_disp,
      lambda = lambda_table
    ),
    parametric = parametric,
    mu.link = object$mu.link,
    sigma2.link = object$sigma2.link,
    loglik = object$loglik,
    aic = object$AIC,
    bic = object$BIC,
    hqic = object$HQ,
    nobs = object$nobs,
    df.residual = object$df.residual,
    iterations = c("optim" = object$iterations, "scoring" = object$scoring),
    method = ifelse(is.null(object$method), "BFGS", object$method),
    R2_N = object$R2_N,
    R2_FC = object$R2_FC,
    P2 = P2,
    residuals = if(!is.null(object$residuals)) object$residuals else NULL,
    residuals.type = if(!is.null(object$residuals)) "quantile" else NULL,
    converged = converged
  )

  class(result) <- "summary.simplexregression"
  return(result)
}

#' @rdname simplexreg.methods
#' @importFrom stats printCoefmat
#' @export
print.summary.simplexregression <- function(x, digits = max(3, getOption("digits") - 3),
                                            ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)),
      "", sep = "\n")

  if(x$converged) {

    if(!is.null(x$residuals)) {
      cat(sprintf("%s:\n", "Quantile residuals"))
      qres <- quantile(x$residuals)
      names(qres) <- c("Min", "1Q", "Median", "3Q", "Max")
      print(formatC(qres, digits = digits, format = "f"), quote = FALSE)
    }

    # Coefficients (Mean model)
    if(NROW(x$coefficients$mean)) {
      cat(sprintf("\nCoefficients (mean model with %s link):\n", x$mu.link))
      coef_mean_print <- x$coefficients$mean
      printCoefmat(coef_mean_print, digits = digits, signif.legend = FALSE)
    } else {
      cat("\nNo coefficients (in mean model)\n")
    }

    # Coefficients (Dispersion model)
    if(NROW(x$coefficients$dispersion)) {
      cat(sprintf("\nDispersion coefficients (dispersion model with %s link):\n",
                  x$sigma2.link))
      coef_disp_print <- x$coefficients$dispersion
      printCoefmat(coef_disp_print, digits = digits, signif.legend = FALSE)
    } else {
      cat("\nNo coefficients (in dispersion model)\n")
    }

    # Lambda parameter (if parametric)
    if(x$parametric && !is.null(x$coefficients$lambda)) {
      cat(sprintf("\nLink function parameter (parametric %s model):\n", x$mu.link))
      lambda_print <- x$coefficients$lambda[, c("Estimate", "Std. Error"), drop = FALSE]
      print(formatC(lambda_print, digits = digits, format = "f"), quote = FALSE)
    }

    if(getOption("show.signif.stars") &&
       any(do.call("rbind", x$coefficients)[, 4L] < 0.1, na.rm = TRUE)) {
      cat("---\nSignif. codes: ",
          "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
    }

    # Model fit statistics
    cat("\nLog-likelihood:", formatC(x$loglik, digits = digits, format = "f"),
        "on", sum(sapply(x$coefficients, NROW)), "Df")
    cat("\nAIC:", formatC(x$aic, digits = digits, format = "f"))
    cat("  BIC:", formatC(x$bic, digits = digits, format = "f"))
    cat("  HQIC:", formatC(x$hqic, digits = digits, format = "f"))
    cat("\nPseudo R-squared (Nagelkerke):", formatC(x$R2_N, digits = digits, format = "f"))
    cat("\nPseudo R-squared (Ferrari and Cribari-Neto):", formatC(x$R2_FC, digits = digits, format = "f"))
    cat("\nP-squared (Espinheira, Silva and Lima):", formatC(x$P2, digits = digits, format = "f"))
    cat("\nNumber of observations:", x$nobs)
    cat(paste("\nNumber of iterations:", x$iterations[1L],
              sprintf("(%s) +", x$method), x$iterations[2L], paste("(Fisher scoring)", "\n")))
  } else{
    cat("model did not converge\n")
   }

  invisible(x)
}

#' @rdname simplexreg.methods
#' @importFrom stats coef
#' @export
coef.simplexregression <- function(object,
                                   model = c("full", "mean", "dispersion"),
                                   ...) {

  model <- match.arg(model)
  cf <- object$coefficients

  switch(model,
         "mean" = {
           cf$mean
         },
         "dispersion" = {
           cf$dispersion
         },
         "full" = {
           nam1 <- names(cf$mean)
           nam2 <- names(cf$dispersion)

           coefs <- c(cf$mean, cf$dispersion)
           names(coefs) <- c(
             nam1,
             if(identical(nam2, "(dispersion)")) "(dispersion)"
             else paste("(dispersion)", nam2, sep = "_")
           )

           if (is_parametric(object)) {
             coefs <- c(coefs, lambda = cf$lambda)
           }

           coefs
         }
  )
}

#' @rdname simplexreg.methods
#' @importFrom stats vcov
#' @export
vcov.simplexregression <- function(object,
                                   model = c("full", "mean", "dispersion"),
                                   ...) {

  model <- match.arg(model)

  vc <- object$vcov
  p <- length(object$coefficients$mean)
  q <- length(object$coefficients$dispersion)

  switch(model,
         "mean" = {
           vc[seq.int(length.out = p), seq.int(length.out = p), drop = FALSE]
         },
         "dispersion" = {
           vc_disp <- vc[seq.int(length.out = q) + p,
                         seq.int(length.out = q) + p,
                         drop = FALSE]
           colnames(vc_disp) <- rownames(vc_disp) <- names(object$coefficients$dispersion)
           vc_disp
         },
         "full" = {
           vc
         }
  )
}

#' @rdname simplexreg.methods
#' @importFrom stats logLik
#' @export
logLik.simplexregression <- function(object, ...) {
  structure(object$loglik,
            df = .npar(object),
            nobs = object$nobs,
            class = "logLik")
}

#' @rdname simplexreg.methods
#' @method fitted simplexregression
#' @export
fitted.simplexregression <- function(object, ...) {
  object$fitted.values
}

#' @rdname simplexreg.methods
#' @importFrom stats na.pass
#' @importFrom Formula as.Formula
#' @method predict simplexregression
#' @export
predict.simplexregression <- function(object, newdata = NULL,
                                      type = c("response", "link", "dispersion"),
                                      ...) {
  type <- match.arg(type)
  parametric <- is_parametric(object)

  # Use fitted values if no new data
  if (is.null(newdata)) {
    return(switch(type,
                  response   = object$fitted.values,
                  link       = list(mean = object$mu.lp, dispersion = object$sigma2.lp),
                  dispersion = object$sigma2.fv
    ))
  }

  formula_full <- object$formula

  deparse1 <- function(x) paste(deparse(x), collapse = " ")
  has_pipe <- function(f) grepl("\\|", deparse1(f), fixed = FALSE)

  # ---- Extract X (mean submodel) -------------------------------------------
  if (has_pipe(formula_full)) {
    f_str       <- deparse1(formula_full)
    rhs_str     <- sub("^.*~\\s*", "", f_str)
    lhs_part    <- trimws(strsplit(rhs_str, "\\|")[[1]][1])
    rhs_formula <- as.formula(paste("~", lhs_part))
  } else {
    rhs_formula <- formula_full
    if (length(rhs_formula) == 3L) rhs_formula <- rhs_formula[-2]
  }

  Terms <- terms(rhs_formula, data = newdata)
  mf    <- model.frame(Terms, data = newdata, na.action = na.pass)
  x_new <- model.matrix(Terms, data = mf)

  # ---- Extract Z (dispersion submodel) -------------------------------------
  if (has_pipe(formula_full)) {
    f_str   <- deparse1(formula_full)
    rhs_str <- sub("^.*~\\s*", "", f_str)
    parts   <- strsplit(rhs_str, "\\|")[[1]]
    if (length(parts) > 1L && trimws(parts[2]) != "1") {
      rhs_disp   <- as.formula(paste("~", trimws(parts[2])))
      Terms_disp <- terms(rhs_disp, data = newdata)
      mf_disp    <- model.frame(Terms_disp, data = newdata, na.action = na.pass)
      z_new      <- model.matrix(Terms_disp, data = mf_disp)
    } else {
      z_new <- matrix(1, nrow = nrow(newdata), ncol = 1L)
    }
  } else {
    z_new <- matrix(1, nrow = nrow(newdata), ncol = 1L)
  }

  # ---- Coefficients --------------------------------------------------------
  beta  <- object$coefficients$mean
  delta <- object$coefficients$dispersion

  # ---- Prediction ----------------------------------------------------------
  if (type %in% c("response", "link")) eta1 <- as.vector(x_new %*% beta)
  if (type %in% c("dispersion", "link")) eta2 <- as.vector(z_new %*% delta)

  switch(type,
         response   = {
           if (parametric)
             parametric_mean_link_inv(eta1, object$lambda.fv, object$mu.link)
           else
             fixed_mean_link_inv(eta1, object$mu.link)
         },
         dispersion = dispersion_link_inv(eta2, object$sigma2.link),
         link       = list(mean = eta1, dispersion = eta2)
  )
}

#' @rdname simplexreg.methods
#' @importFrom stats nobs
#' @export
nobs.simplexregression <- function(object, ...) {
  object$nobs
}

#' @rdname simplexreg.methods
#' @export
df.residual.simplexregression <- function(object, ...) {
  object$nobs - .npar(object)
}

#' @rdname simplexreg.methods
#' @importFrom stats deviance
#' @export
deviance.simplexregression <- function(object, ...) {
  y <- if(is.null(object$y)) model.response(model.frame(object)) else object$y
  mu <- object$mu.fv
  yoneminy <- y * (1 - y)
  muonemu <- mu * (1 - mu)
  sum(((y - mu) / muonemu)^2 / yoneminy)
}

#' @rdname simplexreg.methods
#' @importFrom stats as.formula
#' @export
formula.simplexregression <- function(x, ...) {
  x$formula
}

#' @rdname simplexreg.methods
#' @export
terms.simplexregression <- function(x, model = c("mean", "dispersion"), ...) {
  model <- match.arg(model)
  x$terms[[model]]
}

#' @rdname simplexreg.methods
#' @export
model.frame.simplexregression <- function(formula, ...) {
  if (!is.null(formula$model)) return(formula$model)

  # Reconstruct model frame from the original call when model was not stored
  call          <- formula$call
  call[[1L]]    <- quote(stats::model.frame)
  call$formula  <- formula$formula
  call$link.mu     <- NULL
  call$link.sigma2 <- NULL
  call$model       <- NULL
  eval(call, parent.frame())
}

#' @rdname simplexreg.methods
#' @export
model.matrix.simplexregression <- function(object, model = c("mean", "dispersion"),
                                           ...) {
  model <- match.arg(model)

  switch(model,
         mean = object$mu.x,
         dispersion = object$sigma2.x)
}

#' @rdname simplexreg.methods
#' @importFrom stats formula update
#' @export
update.simplexregression <- function(object, formula., ..., evaluate = TRUE) {
  call <- object$call

  if(is.null(call)) {
    stop("need an object with call component")
  }

  extras <- match.call(expand.dots = FALSE)$...

  if(!missing(formula.)) {
    call$formula <- formula(update(Formula::Formula(formula(object)), formula.))
  }

  if(length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if(!all(existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }

  if(evaluate) {
    eval(call, parent.frame())
  } else {
    call
  }
}

#' @rdname simplexreg.methods
#' @importFrom stats simulate
#' @export
simulate.simplexregression <- function(object, nsim = 1, seed = NULL, ...) {

  if (!is.null(seed)) set.seed(seed)
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
  RNGstate <- get(".Random.seed", envir = .GlobalEnv)

  mu <- object$mu.fv
  sigma2 <- object$sigma2.fv
  n <- length(mu)
  nm <- names(mu)

  s <- replicate(nsim, rsimplex(n, mu = mu, sigma2 = sigma2))

  s <- as.data.frame(s)
  names(s) <- paste("sim", seq_len(nsim), sep = "_")

  if (!is.null(nm)) {
    row.names(s) <- nm
  }

  attr(s, "seed") <- RNGstate
  return(s)
}

#' @rdname simplexreg.methods
#' @importFrom stats AIC
#' @export
AIC.simplexregression <- function(object, ..., k = 2) {
  objects <- list(object, ...)

  aic_fn <- function(obj) -2 * obj$loglik + k * .npar(obj)

  if (length(objects) == 1) return(aic_fn(object))

  model_names <- c(deparse(substitute(object)),
                   sapply(substitute(list(...)), deparse)[-1])
  data.frame(
    df  = sapply(objects, .npar),
    AIC = sapply(objects, aic_fn),
    row.names = model_names
  )
}

#' @rdname simplexreg.methods
#' @importFrom stats BIC
#' @export
BIC.simplexregression <- function(object, ...) {
  objects <- list(object, ...)

  bic_fn <- function(obj) -2 * obj$loglik + log(obj$nobs) * .npar(obj)

  if (length(objects) == 1) return(bic_fn(object))

  model_names <- c(deparse(substitute(object)),
                   sapply(substitute(list(...)), deparse)[-1])
  data.frame(
    df  = sapply(objects, .npar),
    BIC = sapply(objects, bic_fn),
    row.names = model_names
  )
}

#' @rdname simplexreg.methods
#' @export
HQIC <- function(object, ...) UseMethod("HQIC")

#' @rdname simplexreg.methods
#' @export
HQIC.simplexregression <- function(object, ...) {
  objects <- list(object, ...)

  hqic_fn <- function(obj) -2 * obj$loglik + 2 * .npar(obj) * log(log(obj$nobs))

  if (length(objects) == 1) return(hqic_fn(object))

  model_names <- c(deparse(substitute(object)),
                   sapply(substitute(list(...)), deparse)[-1])
  data.frame(
    df   = sapply(objects, .npar),
    HQIC = sapply(objects, hqic_fn),
    row.names = model_names
  )
}

#' @rdname simplexreg.methods
#' @importFrom stats hatvalues
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
  XtX <- crossprod(Xw)
  Inv <- tryCatch(
    chol2inv(chol(XtX)),
    error = function(e) stop("Hat matrix computation failed: design matrix may be singular.")
  )

  rowSums((Xw %*% Inv) * Xw)
}

#' @rdname simplexreg.methods
#' @importFrom stats residuals hatvalues cooks.distance
#' @export
cooks.distance.simplexregression <- function(model, type = c("pearson", "weighted"), ...) {
  type <- match.arg(type)
  h <- hatvalues(model)

  switch(type,
         "pearson" = {
           h * (residuals(model, type = "pearson") ^ 2) / (2*((1-h)^2))
         },
         "weighted" = {
           (residuals(model, type = "weighted") ^ 2) * (h / (1 - h))
         }
  )
}

#' @rdname simplexreg.methods
#' @importFrom sandwich bread
#' @export
bread.simplexregression <- function(x, ...) {
  x$nobs * vcov(x)
}

#' @rdname simplexreg.methods
#' @importFrom sandwich estfun
#' @importFrom stats plogis
#' @export
estfun.simplexregression <- function(x, ...) {

  y <- if(is.null(x$y)) model.response(model.frame(x)) else x$y
  mu <- x$mu.fv
  sigma2 <- x$sigma2.fv
  eta1 <- x$mu.lp
  eta2 <- x$sigma2.lp
  X <- x$mu.x
  Z <- x$sigma2.x
  link.mu <- x$mu.link
  n <- length(y)
  wts <- if(!is.null(x$weights)) x$weights else rep(1, n)
  parametric <- is_parametric(x)

  # Derivative of mu with respect to eta1
  if(parametric) {
    lambda <- x$lambda.fv
    dmu_deta <- as.vector(parametric_mean_link_inv_deriv1(eta1, lambda, link.mu))

    # Compute rho (derivative w.r.t. lambda)
    if(link.mu == "plogit2") {
      exp_aval_frac <- plogis(eta1)^(1/lambda)
      rho <- as.vector(-exp_aval_frac * log(plogis(eta1)) / (lambda^2))
    } else {
      rho <- as.vector((-1/(lambda^2)) * ((exp(eta1) + 1)^(-1/lambda)) *
                         log(exp(eta1) + 1))
    }
  } else {
    dmu_deta <- fixed_mean_link_inv_deriv1(eta1, link.mu)
  }

  # Derivative of sigma2 with respect to eta2
  dsigma2_deta <- dispersion_link_inv_deriv1(eta2, x$sigma2.link)

  diff <- as.vector(y - mu)
  muonemu <- as.vector(mu * (1 - mu))

  yoneminy <- y * (1 - y)
  dev <- ((y - mu) / muonemu)^2 / yoneminy

  U <- wts * (dev / (sigma2 * muonemu) + 1 / (sigma2 * (muonemu)^3))
  a <- wts * ((-1 / (2 * sigma2)) + (dev / (2 * sigma2^2)))

  score_beta <- (dmu_deta * U * diff) * X
  score_delta <- (dsigma2_deta * a) * Z

  if(parametric) {
    score_lambda <- U * rho * diff
    rval <- cbind(score_beta, score_delta, score_lambda)
  } else {
    rval <- cbind(score_beta, score_delta)
  }

  colnames(rval) <- names(coef(x, model = "full"))
  rval[wts <= 0, ] <- 0

  return(rval)
}

#' @rdname simplexreg.methods
#' @importFrom lmtest coeftest
#' @export
coeftest.simplexregression <- function(x, vcov. = NULL, df = Inf, ...) {
  if (is.null(vcov.)) vcov. <- x$vcov
  lmtest::coeftest.default(x, vcov. = vcov., df = df, ...)
}

#' @rdname simplexreg.methods
#' @importFrom lmtest lrtest
#' @export
lrtest.simplexregression <- function(object, ...) {
  lmtest::lrtest.default(object, ...)
}

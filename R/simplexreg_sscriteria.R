################################################################################
#                  SIMPLEX REGRESSION - SCOUT SCORE CRITERION                  #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2026-05-26                                                             #
# Description: Scout Score criterion for simplex regression model selection    #
################################################################################

#' @title Scout Score Criterion for Simplex Regression Model Selection
#' @description Implements the Scout Score (SS) criterion for selecting among
#' competing simplex regression models with parametric and fixed mean link functions.
#'
#' @param ... Two or more objects of class \code{"simplexregression"} to be compared.
#' @param kappa A numeric value controlling the additional penalty for the link mean
#' parameter, \eqn{\kappa \geq 0}. Default is \code{0.1}. Use \code{kappa = 0} for
#' standard Scout Score.
#' @param verbose Logical. If \code{TRUE} (default), prints the SS values for
#' all models and the selected model. If \code{FALSE}, returns results silently.
#' @param digits Integer specifying the number of decimal places for output.
#' Default is \code{max(3, getOption("digits") - 3)}.
#'
#' @details
#' The Scout Score criterion, originally proposed by Costa et al. (2024) for
#' selecting link functions in the \eqn{\beta\text{ARMA}} (beta autoregressive
#' moving average) model, extends Vuong's test statistic to compare \eqn{M \geq 2}
#' competing non-nested models using their individual log-likelihood contributions.
#'
#' For each candidate model \eqn{j \in {1, \ldots, M}}, the Scout Score is defined as:
#' \deqn{SS_j = 1 - M + \sum_{k=1, k \neq j}^M (1 + \dot{\Delta}_{jk})^2,}
#'
#' where \eqn{\dot{\Delta}_{jk} = \max\{0, \Delta_{jk}\}} and \eqn{\Delta_{jk}}
#' is Vuong's (1989) test statistic comparing models \eqn{j} and \eqn{k}, penalized
#' by a term \eqn{\delta_{jk}} that combines the difference in parameter-vector
#' dimensions with, when applicable, a link-complexity penalty controlled by
#' \code{kappa}. The model with the highest Scout Score is selected as the most
#' adequate. For the full derivation of \eqn{\Delta_{jk}}, \eqn{\delta_{jk}}, and
#' the rationale behind the link-complexity penalty, see Justino and
#' Cribari-Neto (2026) and Costa et al. (2024).
#'
#' When at least one of the candidate models does not use a parametric mean
#' link function, \eqn{\kappa} is internally set to \code{0} (see \strong{Important}
#' below), so that \eqn{\delta_{jk}} reduces to the classical dimension penalty
#' with no link-complexity term.
#'
#' The model with the highest Scout Score is selected as the most adequate.
#'
#' \strong{Important}: The penalty term \eqn{\delta_{jk}} is only applied when
#' \emph{all} candidate models employ a parametric mean link function, in which
#' case \code{kappa} (default \code{0.1}) controls the additional penalty. In
#' any other case --- whether all candidate models use fixed mean links, or the
#' set of candidate models mixes parametric and fixed mean links --- the penalty
#' is disabled and the standard, unpenalized Scout Score is computed for all
#' models (\code{kappa} is internally reset to \code{0}). If the user explicitly
#' requested \code{kappa > 0} in either of these situations, a warning is issued;
#' if \code{kappa} was left at its default value, no warning is issued, since
#' falling back to the standard Scout Score is the expected behavior.
#'
#' \strong{Note}: all candidate models must be fitted to the same response
#' vector \code{y}; the function verifies this and stops with an error if the
#' response vectors differ.
#'
#' @return A data frame with rows named after the candidate models and two columns:
#' \describe{
#'   \item{\code{df}}{Number of estimated parameters in each model.}
#'   \item{\code{SS}}{Scout Score value. The model with the highest value is
#'   the selected one.}
#' }
#' When \code{verbose = TRUE}, the selected model is also printed to the console.
#' The data frame is returned invisibly in this case, and visibly when
#' \code{verbose = FALSE}.
#'
#' @references
#' Justino, M. E. C. and Cribari-Neto, F. (2026).
#' Simplex regression with a flexible logit link: Inference and application
#' to cross-country impunity data.
#' \emph{Applied Mathematical Modelling}, \bold{154}, 116713. \doi{10.1016/j.apm.2025.116713}
#'
#' Costa, E., Cribari-Neto, F. and Scher, V. T. (2024).
#' Test inferences and link function selection in dynamic beta modeling of seasonal
#' hydro-environmental time series with temporary abnormal regimes.
#' \emph{Journal of Hydrology}, \bold{638}, 131489.
#' \doi{10.1016/j.jhydrol.2024.131489}
#'
#' Vuong, Q. H. (1989). Likelihood ratio tests for model selection and non-nested
#' hypotheses. \emph{Econometrica}, \bold{57}(2), 307--333. \doi{10.2307/1912557}
#'
#' @examples
#' # Simulate data
#' set.seed(2026)
#' n <- 100
#' x1 <- runif(n, 0, 1)
#' x2 <- runif(n, 0, 1)
#' z1 <- runif(n, 0, 1)
#' mu <- parametric_mean_link_inv(0.6 - 2*x1 - 1.5*x2, 0.5, "plogit1")
#' sigma2 <- dispersion_link_inv(-2 - 2.5*z1, "log")
#' y <- rsimplex(n, mu, sigma2)
#' data <- data.frame(y = y, x1 = x1, x2 = x2, z1 = z1)
#'
#' # Fit models
#' fit1 <- simplexreg(y ~ x1 + x2 | z1, data = data, link.mu = "plogit1")
#' fit2 <- simplexreg(y ~ x1 + x2 | z1, data = data, link.mu = "plogit2")
#' fit3 <- simplexreg(y ~ x1 + x2 | z1, data = data, link.mu = "logit")
#' fit4 <- simplexreg(y ~ x1 + x2 | z1, data = data, link.mu = "probit")
#'
#' # Compare models with verbose output
#' result <- penalized.ss(fit1, fit2, kappa = 0.1)
#'
#' # Compare models silently
#' result <- penalized.ss(fit1, fit2, kappa = 0.1, verbose = FALSE)
#'
#' # Use standard Scout Score (no parametric link penalty)
#' result <- penalized.ss(fit1, fit2, fit3, fit4, kappa = 0)
#'
#' @seealso \code{\link{penalized.ic}}
#' @export
penalized.ss <- function(..., kappa = 0.1, verbose = TRUE,
                         digits = max(3, getOption("digits") - 3)) {

  models <- list(...)
  M <- length(models)

  if (M < 2) stop("At least 2 models are required for comparison")

  # Verify all objects are simplexregression models
  if (!all(sapply(models, function(x) inherits(x, "simplexregression")))) {
    stop("All arguments must be objects of class 'simplexregression'")
  }

  if (!is.numeric(kappa) || length(kappa) != 1L || is.na(kappa) || kappa < 0) {
    stop("'kappa' must be a single non-negative numeric value.")
  }

  # Verify all models were fitted to the same response vector
  y_list <- lapply(models, function(m) as.vector(m$y))
  if (!all(vapply(y_list[-1], function(yy) isTRUE(all.equal(yy, y_list[[1]])),
                  logical(1L)))) {
    stop("All models must be fitted to the same response vector 'y' ",
         "for the Scout Score comparison to be meaningful.")
  }

  # Get model names (supports named arguments)
  model_names <- names(models)
  if (is.null(model_names)) {
    model_names <- sapply(as.list(substitute(list(...)))[-1], deparse)
  }

  y <- models[[1]]$y
  n <- models[[1]]$nobs

  # Initialize storage lists
  mu_list <- sigma2_list <- lambda_list <- f_list <- vector("list", M)
  has_lambda <- logical(M)
  r_list <- numeric(M)

  # Extract values from all models
  for (i in 1:M) {
    mu_list[[i]] <- models[[i]]$mu.fv
    sigma2_list[[i]] <- models[[i]]$sigma2.fv
    f_list[[i]] <- dsimplex(y, mu_list[[i]], sigma2_list[[i]])
    r_list[i] <- n - models[[i]]$df.residual

    # check if model has lambda
    if (!is.null(models[[i]]$lambda.fv) && !all(is.na(models[[i]]$lambda.fv))) {
      lambda_list[[i]] <- models[[i]]$lambda.fv
      has_lambda[i] <- TRUE
    } else {
      lambda_list[[i]] <- NA
      has_lambda[i] <- FALSE
    }
  }

  all_parametric <- all(has_lambda)

  # Detect if user explicitly provided kappa
  call <- match.call()
  user_kappa <- "kappa" %in% names(call)

  # Adjust kappa intelligently
  if (!all_parametric) {
    if (user_kappa && kappa > 0) {
      warning("kappa > 0 is only valid when all models have parametric mean links. ",
              "Setting kappa = 0 automatically for this comparison.")
      kappa <- 0
    } else if (!user_kappa && kappa > 0) {
      # User didn't specify kappa, using default but models don't support it
      kappa <- 0
      # No warning - it's expected behavior
    }
  }

  # Calculate SS for each model
  SS_values <- numeric(M)

  for (i in 1:M) {
    delta_squared_sum <- 0

    for (j in 1:M) {
      if (i != j) {
        # Calculate log-likelihood ratio
        log_ratio <- log(f_list[[i]]/f_list[[j]])

        # Calculate omega_ij^2
        omega2_ij <- 1/n * sum((log_ratio)^2) - (mean(log_ratio)^2)

        # Protect against numerical issues (near-zero variance of log-ratio)
        if (omega2_ij < 1e-10) omega2_ij <- 1e-10

        # Penalization delta_jk
        if (all_parametric) {
          penalty <- kappa * (abs(log(lambda_list[[i]])) -
                                abs(log(lambda_list[[j]])))
        } else {
          penalty <- 0
        }

        # Calculate delta_ij
        delta_ij <- (sum(log_ratio) - 0.5 * ((r_list[i] - r_list[j]) + penalty)
                     * log(n)) / (sqrt(n) * sqrt(omega2_ij))

        # Apply max(0, delta_jk)
        delta_ij_max <- max(0, delta_ij)

        # Add to sum
        delta_squared_sum <- delta_squared_sum + (1 + delta_ij_max)^2
      }
    }

    # Calculate final SS value
    SS_values[i] <- 1 - M + delta_squared_sum
  }

  result <- data.frame(
    df = r_list,
    SS = SS_values,
    row.names = model_names
  )

  # Verbose output
  if (verbose) {
    cat("\n")
    if (kappa == 0) {
      cat("Scout Score values:\n")
    } else {
      cat(sprintf("Penalized Scout Score values (kappa = %s):\n",
                  formatC(kappa, digits = 3, format = "f")))
    }

    result_print <- data.frame(
      df = formatC(result$df, digits = 0, format = "f"),
      SS = formatC(result$SS, digits = digits, format = "f"),
      row.names = rownames(result)
    )
    print(result_print, quote = FALSE, right = TRUE)

    best_model_index <- which.max(SS_values)
    cat(sprintf("\nSelected model: %s\n", model_names[best_model_index]))

    return(invisible(result))
  } else {
    return(result)
  }
}

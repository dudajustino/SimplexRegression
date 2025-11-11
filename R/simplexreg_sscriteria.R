################################################################################
#                  SIMPLEX REGRESSION - SCOUT SCORE CRITERION                  #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2025-11-08                                                             #
# Description: Scout Score criterion for simplex regression model selection    #
################################################################################

# ==============================================================================
# SCOUT SCORE CRITERION
# ==============================================================================

#' @title Scout Score Criterion for Simplex Regression Model Selection
#' @description Implements the Scout Score (SS) criterion for selecting among
#' competing simplex regression models with parametric and fixed mean link functions.
#'
#' @param ... Two or more objects of class \code{"simplexregression"} to be compared
#' @param kappa Numeric value controlling the additional penalty for the link mean
#' parameter. Default is 0.1. Use \code{kappa = 0} for standard Scout Score
#' @param verbose Logical. If \code{TRUE} (default), prints the SS values for
#' all models and the selected model. If \code{FALSE}, returns results silently
#'
#' @details
#' The Scout Score criterion, originally proposed by Costa et al. (2022) for
#' \eqn{\beta\text{ARMA}} A (beta autoregressive moving average) models,
#' extends Vuong's test statistic to compare M ≥ 2 competing
#' non-nested models using their individual log-likelihood contributions.
#'
#' For each candidate model \eqn{j \in {1, ..., M}}, the Scout Score is defined as:
#' \deqn{SS_j = 1 - M + \sum_{k=1, k \neq j}^M (1 + \dot{\Delta}_{jk})^2}
#'
#' where \eqn{\dot{\Delta}_{jk} = \max\{0, \Delta_{jk}\}} and \eqn{\Delta_{jk}}
#' is Vuong's test statistic comparing models j and k:
#' \deqn{\Delta_{jk} = n^{-1/2} \hat{\omega}_{jk}^{-1} \left[\sum_{i=1}^n
#' \log\frac{f_{ji}(\hat{\theta}_j)}{f_{ki}(\hat{\theta}_k)} - \delta_{jk}\right]}
#'
#' For simplex regression models with parametric mean link function, the penalization
#' term is:
#' \deqn{\delta_{jk} = \frac{1}{2}[(r_j - r_k) + \kappa(|\log(\hat{\lambda}_j)| -
#' |\log(\hat{\lambda}_k)|)]\log(n)}
#'
#' where \eqn{\kappa \geq 0} controls the additional penalty associated with the link
#' parameter. The term \eqn{\kappa(|\log(\hat{\lambda}_j)| - |\log(\hat{\lambda}_k)|)}
#' measures link complexity on a logarithmic scale, ensuring symmetry around \eqn{\lambda = 1}.
#'
#' The model with the highest Scout Score is selected as most adequate.
#'
#' **Important**: This penalized version \eqn{SS^{(\lambda)}} should only be used
#' when all candidate models employ a parametric mean link function. For models with fixed
#' mean links, use \code{kappa = 0}.
#'
#' @return A list of class \code{"ss.simplexregression"} containing:
#' \itemize{
#'   \item \code{ss_values}: Numeric vector of Scout Score values for all models
#'   \item \code{best_model}: Index of the selected model
#'   \item \code{kappa}: The penalty parameter used
#'   \item \code{n_models}: Number of models compared
#' }
#'
#' @references
#' E. Costa, F. Cribari-Neto, V. T. Scher, Test inferences and link function selection
#'  in dynamic beta modeling of seasonal hydro-environmental time series with temporary
#'  abnormal regimes, Journal of Hydrology 638 (2024) 131489.
#'  https://doi.org/10.1016/j.jhydrol.2024.131489.
#'
#'Q. H. Vuong, Likelihood ratio tests for model selection and non-nested hypotheses,
#'Econometrica 57 (1989) 307–333. https://doi.org/10.2307/1912557.
#'
#' @examples
#' \dontrun{
#' # Fit multiple models
#' model1 <- simplexreg(y ~ x1, link.mu = "plogit1", data = mydata)
#' model2 <- simplexreg(y ~ x1 + x2, link.mu = "plogit2", data = mydata)
#' model3 <- simplexreg(y ~ x1, link.mu = "logit", data = mydata)
#' model4 <- simplexreg(y ~ x1 + x2, link.mu = "probit", data = mydata)
#'
#' # Compare models with verbose output
#' result <- ss.simplexreg(model1, model2, kappa = 0.1)
#'
#' # Compare models silently
#' result <- ss.simplexreg(model1, model2, kappa = 0.1, verbose = FALSE)
#'
#' # Use standard Scout Score (no parametric link penalty)
#' result <- ss.simplexreg(model3, model4, kappa = 0)
#' }
#'
#' @export
ss.simplexreg <- function(..., kappa = 0.1, verbose = TRUE) {
  # Collect models
  models <- list(...)
  M <- length(models)

  if (M < 2) {
    stop("At least 2 models are required for comparison")
  }

  # Verify all objects are simplexregpar models
  if (!all(sapply(models, function(x) inherits(x, "simplexregression")))) {
    stop("All arguments must be objects of class 'simplexregression'")
  }

  y <- models[[1]]$y
  n <- models[[1]]$nobs

  # Initialize storage lists
  mu_list <- vector("list", M)
  sigma2_list <- vector("list", M)
  lambda_list <- vector("list", M)
  has_lambda <- logical(M)
  f_list <- vector("list", M)
  r_list <- numeric(M)

  # Extract values from all models
  for (i in 1:M) {
    mu_list[[i]] <- models[[i]]$mu.fv
    sigma2_list[[i]] <- models[[i]]$sigma2.fv
    f_list[[i]] <- dsimplex_opt(y, mu_list[[i]], sigma2_list[[i]])
    r_list[i] <- models[[i]]$df.residual

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
  # Prevent misuse of kappa
  if (!all_parametric && kappa > 0) {
    stop("kappa > 0 is only allowed when all models have parametric mean links. Set kappa = 0 for fixed-link models.")
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
        omega2_ij <- 1/n * sum((log_ratio)^2) - (1/n * sum(log_ratio))^2

        # Penalization delta_jk
        if (all_parametric) {
          penalty <- kappa * (abs(log(lambda_list[[i]])) - abs(log(lambda_list[[j]])))
        } else {
          penalty <- 0
        }

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

  # Find index of model with highest SS value
  best_model_index <- which.max(SS_values)

  # Verbose output
  if (verbose) {
    cat("\n")
    cat("============================================================\n")
    cat("           Scout Score Model Selection Results\n")
    cat("============================================================\n")
    cat("\nPenalty parameter (kappa):", kappa, "\n")
    cat("Number of models compared:", M, "\n\n")
    cat("Scout Score values:\n")
    cat("------------------------------------------------------------\n")
    for (i in 1:M) {
      marker <- if (i == best_model_index) " *" else "  "
      cat(sprintf("  Model %d: %8.4f%s\n", i, SS_values[i], marker))
    }
    cat("------------------------------------------------------------\n")
    cat("\nSelected model: Model", best_model_index, "\n")
    cat("(* indicates selected model)\n")
    cat("============================================================\n\n")
  }

  return(invisible(list(
    ss_values = SS_values,
    best_model = best_model_index,
    kappa = kappa,
    n_models = M
  )))
}


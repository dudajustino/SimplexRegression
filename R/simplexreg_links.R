################################################################################
#                 LINK FUNCTIONS - SIMPLEX REGRESSION MODELS                   #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2025-11-08                                                             #
# Description: Link functions, their inverses, and derivatives for simplex     #
#              regression models with parametric and fixed mean link functions #
################################################################################

# ==============================================================================
# 1. PARAMETRIC MEAN LINK FUNCTIONS (plogit1 and plogit2)
# ==============================================================================

#' @title Parametric Mean Link Functions and Derivatives
#'
#' @description
#' Provides the parametric mean link functions, their inverses, and derivatives
#' for the simplex regression model. Two parametric link types are supported:
#' \code{"plogit1"} and \code{"plogit2"}.
#'
#' @details
#' Two parametric mean link functions are available:
#' \itemize{
#'   \item \strong{plogit2}: \eqn{g(\mu; \lambda) = \log(\mu^\lambda / (1 - \mu^\lambda))}
#'   \item \strong{plogit1}: \eqn{g(\mu; \lambda) = \log((1-\mu)^{-\lambda} - 1)}
#' }
#' Their inverses and derivatives with respect to \eqn{\mu} are also implemented.
#'
#' @param mu Mean parameter (numeric vector, 0 < mu < 1)
#' @param eta Linear predictor of mean (numeric vector)
#' @param lambda Power parameter (numeric scalar, positive)
#' @param link_type Type of link function: \code{"plogit1"} or \code{"plogit2"}
#'
#' @importFrom stats plogis qlogis
#' @return Numeric vector with transformed values
#' @examples
#' parametric_mean_link(0.2, lambda = 1.2, link_type = "plogit2")
#' parametric_mean_link(c(0.2, 0.5, 0.8), lambda = 1.5, link_type = "plogit1")
#' parametric_mean_link_inv(0, lambda = 1, link_type = "plogit2")
#' parametric_mean_link_deriv1(0.5, lambda = 1, link_type = "plogit2")
#' parametric_mean_link_inv_deriv1(0, lambda = 1, link_type = "plogit2")
#' parametric_mean_link_deriv2(0.5, lambda = 1, link_type = "plogit2")
#'
#' @name parametric_mean_links
NULL

#' @rdname parametric_mean_links
#' @export
parametric_mean_link <- function(mu, lambda, link_type = c("plogit1", "plogit2")) {
  link_type <- match.arg(link_type)

  # Input validation
  if (any(mu <= 0 | mu >= 1)) {
    stop("All values of 'mu' must be in the interval (0, 1)")
  }

  if(link_type == "plogit2") {
    mu_lambda <- pmin(pmax(mu^lambda, .Machine$double.eps), 1-.Machine$double.eps)
    result <- qlogis(mu_lambda)
  } else {
    one_minus_mu_lambda <- pmin(pmax(1-(1-mu)^lambda, .Machine$double.eps),
                                1-.Machine$double.eps)
    result <- qlogis(one_minus_mu_lambda)
  }

  if (any(is.nan(result))) {
    stop("NaN values detected in parametric_mean_link. Check input values.")
  }

  return(result)
}

#' @rdname parametric_mean_links
#' @export
parametric_mean_link_inv <- function(eta, lambda, link_type = c("plogit1", "plogit2")) {
  link_type <- match.arg(link_type)

  if(link_type == "plogit2") {
    result <- pmin(pmax(plogis(eta)^(1/lambda), .Machine$double.eps), 1-.Machine$double.eps)
  } else {
    result <- pmin(pmax(1-(exp(eta)+1)^(-1/lambda), .Machine$double.eps), 1-.Machine$double.eps)
  }

  return(result)
}

#' @rdname parametric_mean_links
#' @export
parametric_mean_link_deriv1 <- function(mu, lambda, link_type = c("plogit1", "plogit2")) {
  link_type <- match.arg(link_type)

  if (any(mu <= 0 | mu >= 1)) {
    stop("All values of 'mu' must be in the interval (0, 1)")
  }

  if(link_type== "plogit2") {
    mu_lambda <- pmin(pmax(mu^lambda, .Machine$double.eps), 1-.Machine$double.eps)
    result <- lambda / (mu * (1 - mu_lambda))
  } else {
    one_minus_mu_lambda <- pmin(pmax(1-(1-mu)^lambda, .Machine$double.eps),
                                1-.Machine$double.eps)
    result <- lambda / ((1-mu) * one_minus_mu_lambda)
  }

  return(result)
}

#' @rdname parametric_mean_links
#' @export
parametric_mean_link_inv_deriv1 <- function(eta, lambda, link_type = c("plogit1", "plogit2")) {
  link_type <- match.arg(link_type)

  if(link_type == "plogit2") {
    result <- parametric_mean_link_inv(eta, lambda, link_type) / (lambda * (1 + exp(eta)))
  } else {
    result <- exp(eta) / (lambda * (1 + exp(eta))^(1/lambda + 1))
  }

  return(result)
}

#' @rdname parametric_mean_links
#' @export
parametric_mean_link_deriv2 <- function(mu, lambda, link_type = c("plogit1", "plogit2")) {
  link_type <- match.arg(link_type)

  if (any(mu <= 0 | mu >= 1)) {
    stop("All values of 'mu' must be in the interval (0, 1)")
  }

  if(link_type == "plogit2") {
    mu_lambda <- pmin(pmax(mu^lambda, .Machine$double.eps), 1-.Machine$double.eps)
    result <- (lambda * mu_lambda * (1+lambda) - lambda) / ((mu*(1-mu_lambda))^2)
  } else {
    one_minus_mu_lambda <- pmin(pmax(1-(1-mu)^lambda, .Machine$double.eps),
                                1-.Machine$double.eps)
    result <- (lambda * (one_minus_mu_lambda - lambda*(1-mu)^lambda )) /
      (((1-mu) * one_minus_mu_lambda)^2)
  }

  return(result)
}

# ==============================================================================
# 2. FIXED MEAN LINK FUNCTIONS
# ==============================================================================

#' @title Fixed Mean Link Functions and Derivatives
#'
#' @description
#' Provides the fixed mean link functions, their inverses, and derivatives
#' for the simplex regression model. Supported link types are:
#' \code{"logit"}, \code{"probit"}, \code{"cloglog"}, \code{"loglog"}, and \code{"cauchit"}.
#'
#' @param mu Mean parameter (numeric vector, 0 < mu < 1)
#' @param eta Linear predictor of mean (numeric vector)
#' @param link_type Type of link function: \code{"logit"}, \code{"probit"},
#'   \code{"cloglog"}, \code{"loglog"}, or \code{"cauchit"}.
#'
#' @details
#' Available link functions:
#' \itemize{
#'   \item \strong{Logit}: \eqn{g(\mu) = \log(\mu/(1-\mu))}
#'   \item \strong{Probit}: \eqn{g(\mu) = \Phi^{-1}(\mu)}
#'   \item \strong{Complementary loglog}: \eqn{g(\mu) = \log(-\log(1-\mu))}
#'   \item \strong{Log-log}: \eqn{g(\mu) = -\log(-\log(\mu))}
#'   \item \strong{Cauchit}: \eqn{g(\mu) = \tan(\pi(\mu - 0.5))}
#' }
#'
#' @importFrom stats pnorm qnorm dnorm pcauchy qcauchy dcauchy plogis dlogis
#' @return A numeric vector corresponding to the evaluated link,
#' its inverse, or derivative depending on the function.
#'
#' @examples
#' fixed_mean_link(0.5, link_type = "logit")
#' fixed_mean_link(c(0.2, 0.5, 0.8), link_type = "probit")
#' fixed_mean_link_inv(eta = 0.2, link_type = "logit")
#' fixed_mean_link_deriv1(mu = 0.5, link_type = "logit")
#'
#' @name fixed_mean_links
NULL

#' @rdname fixed_mean_links
#' @export
fixed_mean_link <- function(mu, link_type = c("logit", "probit", "cloglog",
                                              "loglog", "cauchit")) {
  link_type <- match.arg(link_type)

  if (any(mu <= 0 | mu >= 1)) {
    stop("All values of 'mu' must be in the interval (0, 1)")
  }

  result <- switch(link_type,
                   loglog = -log(-log(mu)),
                   cloglog = log(-log(1.0 - mu)),
                   probit = qnorm(mu),
                   cauchit = qcauchy(mu),
                   logit = log(mu/(1.0 - mu))
  )

  return(result)
}

#' @rdname fixed_mean_links
#' @export
fixed_mean_link_inv <- function(eta, link_type = c("logit", "probit", "cloglog",
                                                   "loglog", "cauchit")) {
  link_type <- match.arg(link_type)

  result <- switch(link_type,
                   loglog = exp(-exp(-eta)),
                   cloglog = 1 - exp(-exp(eta)),
                   probit = pnorm(eta),
                   cauchit = pcauchy(eta),
                   logit = plogis(eta)
  )

  return(result)
}

#' @rdname fixed_mean_links
#' @export
fixed_mean_link_deriv1 <- function(mu, link_type = c("logit", "probit", "cloglog",
                                                     "loglog", "cauchit")) {
  link_type <- match.arg(link_type)

  if (any(mu <= 0 | mu >= 1)) {
    stop("All values of 'mu' must be in the interval (0, 1)")
  }

  result <- switch(link_type,
                   loglog = -1 / (mu * log(mu)),
                   cloglog = -1 / ((1 - mu) * log(1 - mu)),
                   probit = 1 / dnorm(qnorm(mu)),
                   cauchit = 1 / dcauchy(qcauchy(mu)),
                   logit = 1 / (mu * (1 - mu))
  )

  return(result)
}

#' @rdname fixed_mean_links
#' @export
fixed_mean_link_deriv2 <- function(mu, link_type = c("logit", "probit", "cloglog",
                                                     "loglog", "cauchit")) {
  link_type <- match.arg(link_type)

  if (any(mu <= 0 | mu >= 1)) {
    stop("All values of 'mu' must be in the interval (0, 1)")
  }

  result <- switch(link_type,
                   loglog = (log(mu) + 1) / (mu^2 * log(mu)^2),
                   cloglog = - (log(1 - mu) + 1) / ((1 - mu)^2 * log(1 - mu)^2),
                   probit = qnorm(mu) / (dnorm(qnorm(mu))^2),
                   cauchit = 2*pi^2 * qcauchy(mu) * (1+qcauchy(mu)^2),
                   logit = (2*mu - 1) / (mu^2 * (1 - mu)^2)
  )

  return(result)
}

#' @rdname fixed_mean_links
#' @export
fixed_mean_link_inv_deriv1 <- function(eta, link_type = c("logit", "probit", "cloglog",
                                                          "loglog", "cauchit")) {
  link_type <- match.arg(link_type)

  result <- switch(link_type,
                   loglog = exp(-eta - exp(-eta)),
                   cloglog = exp(eta - exp(eta)),
                   probit = dnorm(eta),
                   cauchit = dcauchy(eta),
                   logit = dlogis(eta)
  )
  return(result)
}

# ==============================================================================
# 3. DISPERSION LINK FUNCTIONS
# ==============================================================================

#' @title Dispersion Link Functions and Their Derivatives
#' @description
#' Provides the link function, its inverse, and derivative for the dispersion
#' submodel in the simplex regression.
#' Supported link types are: \code{"log"}, \code{"sqrt"}, and \code{"identity"}.
#'
#' @param sigma2 Dispersion parameter (numeric vector, sigma2 > 0)
#' @param eta Linear predictor of dispersion (numeric vector)
#' @param link_type Type of link: \code{"log"}, \code{"sqrt"}, or \code{"identity"}.
#'
#' @details
#' Available link functions:
#' \itemize{
#'   \item \strong{Log}: \eqn{h(\sigma^2) = \log(\sigma^2)} (ensures positivity)
#'   \item \strong{Sqrt}: \eqn{h(\sigma^2) = \sqrt{\sigma^2}}
#'   \item \strong{Identity}: \eqn{h(\sigma^2) = \sigma^2} (no transformation)
#' }
#'
#' @return Numeric vector with transformed values
#' @examples
#' dispersion_link(1.5, link_type = "log")
#' dispersion_link(c(0.5, 1, 2), link_type = "sqrt")
#' dispersion_link_inv(0, link_type = "log")
#' dispersion_link_deriv1(1, link_type = "log")
#' dispersion_link_inv_deriv1(0, link_type = "log")
#'
#' @name dispersion_links
NULL

#' @rdname dispersion_links
#' @export
dispersion_link <- function(sigma2, link_type = c("log", "sqrt", "identity")) {
  link_type <- match.arg(link_type)

  if (any(sigma2 <= 0)) {
    stop("All values of 'sigma2' must be positive")
  }

  result <- switch(link_type,
                   log = log(sigma2),
                   sqrt = sqrt(sigma2),
                   identity = sigma2
  )

  return(result)
}

#' @rdname dispersion_links
#' @export
dispersion_link_inv <- function(eta, link_type = c("log", "sqrt", "identity")) {
  link_type <- match.arg(link_type)

  result <- switch(link_type,
                   log = exp(eta),
                   sqrt = eta^2,
                   identity = eta
  )

  return(result)
}

#' @rdname dispersion_links
#' @export
dispersion_link_deriv1 <- function(sigma2, link_type = c("log", "sqrt", "identity")) {
  link_type <- match.arg(link_type)

  if (any(sigma2 <= 0)) {
    stop("All values of 'sigma2' must be positive")
  }

  result <- switch(link_type,
                   log = 1 / sigma2,
                   sqrt = 0.5 * sigma2^(-0.5),
                   identity = rep.int(1, length(sigma2))
  )

  return(result)
}

#' @rdname dispersion_links
#' @export
dispersion_link_inv_deriv1 <- function(eta, link_type = c("log", "sqrt", "identity")) {
  link_type <- match.arg(link_type)

  result <- switch(link_type,
                   log = exp(eta),
                   sqrt = 2 * eta,
                   identity = rep.int(1, length(eta))
  )

  return(result)
}

################################################################################
#       SIMPLEX DISTRIBUTION FUNCTIONS - DENSITY, CDF, QUANTILE, RANDOM        #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2025-11-08                                                             #
# Description: Probability density, cumulative distribution, quantile, and     #
#              random generation functions for the simplex distribution        #
################################################################################

# ==============================================================================
# 1. DENSITY FUNCTION
# ==============================================================================

#' @title Simplex Distribution Functions
#' @name simplex_opt
#' @description
#' Provides the main distributional functions for the simplex distribution:
#' \itemize{
#'   \item \code{dsimplex_opt()}: density function
#'   \item \code{psimplex_opt()}: cumulative distribution function (CDF)
#'   \item \code{qsimplex_opt()}: quantile function (inverse CDF)
#'   \item \code{rsimplex_opt()}: random number generation
#' }
#'
#' The simplex distribution is defined on the unit interval \eqn{(0, 1)},
#' with parameters mean \eqn{\mu} and dispersion \eqn{\sigma^2 > 0}
#'
#' @param x Numeric vector of quantiles (for density).
#' @param q Numeric vector of quantiles (for CDF).
#' @param p Numeric vector of probabilities (for quantile function).
#' @param mu Mean parameter (\eqn{0 < \mu < 1})
#' @param sigma2 Dispersion parameter (\eqn{\sigma^2 > 0})
#' @param n Number of observations (for random generation).
#'
#' @details
#' The probability density function of the simplex distribution is:
#' \deqn{f(y; \mu, \sigma^2) = \frac{1}{\sqrt{2\pi\sigma^2[y(1-y)]^3}}
#'       \exp\left(-\frac{1}{2\sigma^2} d(y; \mu)\right),}
#' where \eqn{y \in (0,1)}, and \eqn{d(x; \mu) = \frac{(y - \mu)^2}{y(1 - y)\mu^2(1 - \mu)^2}} is the
#' deviance component of the simplex model.
#'
#' @importFrom stats pnorm qnorm integrate rchisq runif uniroot
#' @return Numeric vector with density, probability, quantile, or simulated values.
#'
#' @references
#' O. E. Barndorff-Nielsen, B. Jørgensen, Some parametric models on the simplex, Journal of
#' Multivariate Analysis 39 (1991) 106–116. https://doi.org/10.1016/0047-259X(91)90008-P.
#'
#' @examples
#' dsimplex_opt(0.5, mu = 0.3, sigma2 = 0.5)
#' psimplex_opt(0.5, mu = 0.3, sigma2 = 0.5)
#' qsimplex_opt(0.5, mu = 0.3, sigma2 = 0.5)
#' rsimplex_opt(5, mu = 0.5, sigma2 = 0.5)
#'
#' @export
dsimplex_opt <- function(x, mu, sigma2) {
  # Truncate extreme values to avoid numerical issues
  x <- pmin(pmax(x, 1e-8), 1 - 1e-8)

  # The distribution is only defined for x between 0 and 1.
  # Values outside this range should have a density of 0.

  # The vectorized `ifelse` function handles this condition for each element of 'x'.
  valid_x <- x > 0 & x < 1

  # The expression below is the simplex PDF, written in a vectorized way.
  # All calculations (numerator, denominator, exponential) are applied element-wise.
  # If mu and sig are scalars, R recycles them to match the length of x.

  numerator <- exp(-((x - mu)^2) / (2 * sigma2 * x * (1 - x) * mu^2 * (1 - mu)^2))
  denominator <- sqrt(2 * pi * sigma2 * (x * (1 - x))^3)

  result <- ifelse(valid_x, numerator / denominator, 0)

  return(result)
}

# ==============================================================================
# 2. CUMULATIVE DISTRIBUTION FUNCTION (CDF)
# ==============================================================================

#' @keywords internal
psimplex.norm_opt <-  function (q, mu, sigma2) {
  return(pnorm((q-mu)/sqrt(sigma2*mu^3*(1-mu)^3)))
}

#' @rdname simplex_opt
#' @export
psimplex_opt <- function(q, mu, sigma2) {
  sig <- sqrt(sigma2)
  # Ensures that all vectors have the same length
  n <- length(q)
  if (length(mu) != n) mu <- rep(mu, length.out = n)
  if (length(sig) != n) sig <- rep(sig, length.out = n)

  # Defines the internal density function
  dsimp <- function(x, mu_val, sig_val) {
    1 / sqrt(2 * pi * sig_val^2 * (x * (1 - x))^3) *
      exp(-0.5 / sig_val^2 * (x - mu_val)^2 / (x * (1 - x) * mu_val^2 * (1 - mu_val)^2))
  }

  # Uses `sapply` to iterate over indices, not directly over `q`
  pp <- sapply(1:n, function(i) {
    qi <- q[i]
    mui <- mu[i]
    sigi <- sig[i]

    # Applies the normal approximation if the condition is met
    if (sigi < 0.001 | (1 - mui) * sigi < 0.01) {
      return(psimplex.norm_opt(qi, mui, sigi^2))
    } else {
      # Calls integration for each set of parameters
      integrate(dsimp, lower = 1e-8, upper = qi, mu_val = mui, sig_val = sigi)$value
    }
  })

  return(pp)
}

# ==============================================================================
# 3. QUANTILE FUNCTION
# ==============================================================================

#' @keywords internal
qsimplex.norm_opt <- function(p, mu, sigma2) {
  # Ensures all input vectors have the same length.
  # This is crucial when 'mu' and 'sigma2' are vectors.
  n <- length(p)
  if (length(mu) != n) mu <- rep(mu, length.out = n)
  if (length(sigma2) != n) sigma2 <- rep(sigma2, length.out = n)

  # Calculates the standard deviation of the normal approximation in a vectorized way.
  # This operation works correctly even if 'mu' and 'sigma2' are vectors.
  sd_approx <- sqrt(sigma2 * mu^3 * (1 - mu)^3)

  # Returns the normal distribution quantile for each probability in 'p',
  # using the corresponding values of mean ('mu') and standard deviation ('sd_approx').
  return(qnorm(p, mean = mu, sd = sd_approx))
}

#' @rdname simplex_opt
#' @export
qsimplex_opt <- function(p, mu, sigma2) {
  # Ensures that all vectors have the same length
  n <- length(p)
  if (length(mu) != n) mu <- rep(mu, length.out = n)
  if (length(sigma2) != n) sigma2 <- rep(sigma2, length.out = n)

  # Uses `sapply` to iterate over indices
  qq <- sapply(1:n, function(i) {
    pi <- p[i]
    mui <- mu[i]
    sigma2i <- sigma2[i]

    # Handles the case of high dispersion
    if (sigma2i > 200) {
      sigma2i <- 200
    }

    if (sigma2i < 0.1) {
      return(qsimplex.norm_opt(pi, mui, sigma2i))
    } else {
      # Calls `uniroot` with the correct scalar parameters for each iteration
      tryCatch({
        uniroot(
        f = function(x) psimplex_opt(q = x, mu = mui, sigma2 = sigma2i) - pi,
        interval = c(1e-8, 1 - 1e-8),
        tol = 1e-6,
        extendInt = "no"
      )$root
      }, error = function(e) NA)
    }
  })

  return(qq)
}

# ==============================================================================
# 4. RANDOM GENERATION
# ==============================================================================

#' @keywords internal
rIG_opt <- function(n, epsilon, Tau) {
  # Generates vectors of length 'n' for z and u1.
  z <- rchisq(n, 1)
  u1 <- runif(n, 0, 1)

  # The vectorized mathematical operations `sqrt`, `*`, `/` are applied element-wise.
  # If 'epsilon' and 'Tau' are scalars, R recycles them to match the length of 'z'.
  ss <- sqrt(4 * epsilon * z / Tau + (epsilon * z)^2)
  z1 <- epsilon + (epsilon^2) * Tau * z / 2 - (epsilon * Tau / 2) * ss

  xxx <- z1

  # The logical indexing handles recycling automatically.
  idx <- (u1 > (epsilon / (epsilon + z1)))

  # Correction: apply the `[idx]` indexing to both the numerator and the denominator.
  # This ensures that the replacement operation is performed with vectors of the same length.
  # The operations on the right side of the assignment are applied only to the elements selected by `idx`.
  # This line is correct as is: R's recycling behavior makes 'epsilon^2' work correctly here.
  # A simple check (if-else) can be added for clarity, but it's not strictly necessary.
  xxx[idx] <- ((epsilon^2) / z1)[idx]

  return(as.numeric(xxx))
}

#' @keywords internal
rMIG_opt <- function(n, epsilon, Tau, mu) {
  # Passes the vectors (or recycled scalars) to rIG_opt.
  x1 <- rIG_opt(n, epsilon, Tau)

  # The other operations also work for vectors.
  x2 <- rchisq(n, 1)
  x3 <- x2 * Tau * (epsilon^2)
  u2 <- runif(n, 0, 1)

  xx <- x1

  # Logical indexing handles recycling automatically.
  idx <- which(u2 < mu)
  xx[idx] <- x1[idx] + x3[idx]

  return(as.numeric(xx))
}

#' @rdname simplex_opt
#' @export
rsimplex_opt <- function(n, mu, sigma2) {
  # The operations below are vectorized.
  # If 'mu' and 'sigma2' are scalars, R recycles them to the length 'n'.
  # If 'mu' and 'sigma2' are vectors of length 'n', the operations are element-wise.
  epsilon <- mu / (1 - mu)
  Tau <- sigma2 * ((1 - mu)^2)

  # The rMIG_opt function is now able to receive vectors for 'epsilon', 'Tau', and 'mu'.
  # This is handled by the next step.
  x <- rMIG_opt(n, epsilon, Tau, mu)

  # The final transformation is also vectorized.
  yy <- x / (1 + x)

  return(as.vector(yy))
}


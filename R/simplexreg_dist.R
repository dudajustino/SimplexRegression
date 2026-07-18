################################################################################
#       SIMPLEX DISTRIBUTION FUNCTIONS - DENSITY, CDF, QUANTILE, RANDOM        #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2026-05-26                                                             #
# Description: Probability density, cumulative distribution, quantile, and     #
#              random generation functions for the simplex distribution        #
################################################################################

# ==============================================================================
# 1. DENSITY FUNCTION
# ==============================================================================

#' @title Simplex Distribution Functions
#' @name simplex_opt
#' @description
#' Density, distribution function, quantile function and random generation for the
#' simplex distribution with parameters mean \eqn{\mu} and dispersion \eqn{\sigma^2}.
#'
#' @param x,q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param mu Mean parameter (\eqn{0 < \mu < 1}).
#' @param sigma2 Dispersion parameter (\eqn{\sigma^2 > 0}).
#' @param n Number of observations.
#' @param log,log.p Logical; if \code{TRUE}, probabilities/densities \eqn{p} are given
#' as \eqn{log(p)}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P[X <= x]},
#'        otherwise, \eqn{P[X > x]}.
#'
#' @details
#' The probability density function of the simplex distribution is given by:
#' \deqn{f(y; \mu, \sigma^2) = \frac{1}{\sqrt{2\pi\sigma^2[y(1-y)]^3}}
#'       \exp\left(-\frac{1}{2\sigma^2} d(y; \mu)\right),}
#' where \eqn{y \in (0, 1)}, and \eqn{d(y; \mu) = \frac{(y - \mu)^2}{y(1 - y)
#' \mu^2(1 - \mu)^2}} is the unit deviance.
#'
#' The cumulative distribution function and the quantile function of the simplex
#' distribution do not admit closed-form expressions. For small values of
#' \eqn{\sigma^2}, \code{psimplex()} and \code{qsimplex()} use the normal
#' approximation implied by the small-dispersion asymptotic theory (Jørgensen, 1997).
#' Otherwise, \code{psimplex()} is computed by numerical integration of the density and
#' \code{qsimplex()} is obtained by numerical root finding.
#'
#' Random generation in \code{rsimplex()} is based on the inverse Gaussian
#' mixture (M-IG) representation of the simplex distribution,
#' followed by the transformation \eqn{Y = X/(1 + X)}.
#'
#' @importFrom stats pnorm qnorm integrate rchisq runif uniroot
#'
#' @return \code{dsimplex} gives the density, \code{psimplex} the
#' distribution function, \code{qsimplex} the quantile function, and
#' \code{rsimplex} generates random deviates.
#'
#' For \code{sigma2} values requiring numerical root-finding (i.e., not small
#' enough for the normal approximation), \code{qsimplex} may return \code{NA}
#' if \code{uniroot} fails to find a root in the given interval.
#'
#' Invalid arguments (\code{mu} outside (0, 1) or non-positive \code{sigma2}) will
#' trigger an error.
#'
#' @references
#' Barndorff-Nielsen, O. E. and Jørgensen, B. (1991).
#' Some parametric models on the simplex.
#' \emph{Journal of Multivariate Analysis}, \bold{39}(1), 106--116.
#' \doi{10.1016/0047-259X(91)90008-P}
#'
#' Jørgensen B (1997). \emph{The Theory of Dispersion Models}.
#' Chapman and Hall, London.
#'
#' @examples
#' dsimplex(0.5, mu = 0.3, sigma2 = 0.5)
#' dsimplex(0.5, mu = 0.3, sigma2 = 0.5, log = TRUE)
#' psimplex(0.5, mu = 0.3, sigma2 = 0.5)
#' psimplex(0.5, mu = 0.3, sigma2 = 0.5, lower.tail = FALSE)
#' qsimplex(0.5, mu = 0.3, sigma2 = 0.5)
#' qsimplex(log(0.5), mu = 0.3, sigma2 = 0.5, log.p = TRUE)
#' rsimplex(5, mu = 0.5, sigma2 = 0.5)
#'
#' @export
dsimplex <- function(x, mu, sigma2, log = FALSE) {
  if (any(mu <= 0 | mu >= 1)) stop("'mu' must be in (0, 1)")
  if (any(sigma2 <= 0))       stop("'sigma2' must be positive")

  # Truncate extreme values to avoid numerical issues
  x <- pmin(pmax(x, 1e-8), 1 - 1e-8)

  # The distribution is only defined for x between 0 and 1.
  # Values outside this range should have a density of 0.

  # The vectorized `ifelse` function handles this condition for each element of 'x'.
  valid_x <- x > 0 & x < 1

  # The expression below is the simplex PDF, written in a vectorized way.
  # All calculations (numerator, denominator, exponential) are applied element-wise.
  # If mu and sig are scalars, R recycles them to match the length of x.

  # Calculate log-density to avoid numerical overflow/underflow
  log_numerator <- -((x - mu)^2) / (2 * sigma2 * x * (1 - x) * mu^2 * (1 - mu)^2)
  log_denominator <- 0.5 * log(2 * pi * sigma2) + 1.5 * log(x * (1 - x))

  log_density <- log_numerator - log_denominator

  # Set invalid values to -Inf (which becomes 0 on exp scale)
  log_density <- ifelse(valid_x, log_density, -Inf)

  # Return log-density or density depending on 'log' argument
  if (log) return(log_density)
  return(exp(log_density))
}

# ==============================================================================
# 2. CUMULATIVE DISTRIBUTION FUNCTION (CDF)
# ==============================================================================

#' @keywords internal
psimplex.norm <-  function (q, mu, sigma2, lower.tail = TRUE, log.p = FALSE) {
  prob <- pnorm(q, mean = mu, sd = sqrt(sigma2 * mu^3 * (1 - mu)^3),
                lower.tail = lower.tail, log.p = log.p)
  return(prob)
}

#' @rdname simplex_opt
#' @export
psimplex <- function(q, mu, sigma2, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu <= 0 | mu >= 1)) stop("'mu' must be in (0, 1)")
  if (any(sigma2 <= 0))       stop("'sigma2' must be positive")

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
      return(psimplex.norm(qi, mui, sigi^2, lower.tail = lower.tail, log.p = log.p))
    } else {
      # Calls integration for each set of parameters
      result <- integrate(dsimp, lower = 1e-8, upper = qi,
                          mu_val = mui, sig_val = sigi)$value

      if (!lower.tail) {
        result <- 1 - result
      }
      if (log.p) {
        result <- log(result)
      }

      return(result)
    }
  })

  return(pp)
}

# ==============================================================================
# 3. QUANTILE FUNCTION
# ==============================================================================

#' @keywords internal
qsimplex.norm <- function(p, mu, sigma2, lower.tail = TRUE, log.p = FALSE) {
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
  return(qnorm(p, mean = mu, sd = sd_approx, lower.tail = lower.tail, log.p = log.p))
}

#' @rdname simplex_opt
#' @export
qsimplex <- function(p, mu, sigma2, lower.tail = TRUE, log.p = FALSE) {

  if (any(mu <= 0 | mu >= 1)) stop("'mu' must be in (0, 1)")
  if (any(sigma2 <= 0))       stop("'sigma2' must be positive")

  # Ensures that all vectors have the same length
  n <- length(p)
  if (length(mu) != n) mu <- rep(mu, length.out = n)
  if (length(sigma2) != n) sigma2 <- rep(sigma2, length.out = n)

  # Uses `sapply` to iterate over indices
  qq <- sapply(1:n, function(i) {
    pi <- p[i]
    mui <- mu[i]
    sigma2i <- sigma2[i]

    pi_transformed <- pi
    if (log.p) {
      pi_transformed <- exp(pi_transformed)
    }
    if (!lower.tail) {
      pi_transformed <- 1 - pi_transformed
    }

    # Handles the case of high dispersion
    if (sigma2i > 200) {
      warning("sigma2 > 200 capped at 200 for numerical stability in qsimplex.")
      sigma2i <- 200
    }

    if (sigma2i < 0.1) {
      return(qsimplex.norm(pi_transformed, mui, sigma2i, lower.tail = TRUE, log.p = FALSE))
    } else {
      # Calls `uniroot` with the correct scalar parameters for each iteration
      tryCatch({
        uniroot(
        f = function(x) psimplex(q = x, mu = mui, sigma2 = sigma2i,
                                 lower.tail = TRUE, log.p = FALSE) - pi_transformed,
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
rIG <- function(n, epsilon, Tau) {
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
  # This ensures that the replacement operation is performed with vectors of the
  # same length.
  # The operations on the right side of the assignment are applied only to the elements
  # selected by `idx`.
  # This line is correct as is: R's recycling behavior makes 'epsilon^2' work correctly here.
  # A simple check (if-else) can be added for clarity, but it's not strictly necessary.
  xxx[idx] <- ((epsilon^2) / z1)[idx]

  return(as.numeric(xxx))
}

#' @keywords internal
rMIG <- function(n, epsilon, Tau, mu) {
  # Passes the vectors (or recycled scalars) to rIG.
  x1 <- rIG(n, epsilon, Tau)

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
rsimplex <- function(n, mu, sigma2) {
  if (any(mu <= 0 | mu >= 1)) stop("'mu' must be in (0, 1)")
  if (any(sigma2 <= 0))       stop("'sigma2' must be positive")

  # The operations below are vectorized.
  # If 'mu' and 'sigma2' are scalars, R recycles them to the length 'n'.
  # If 'mu' and 'sigma2' are vectors of length 'n', the operations are element-wise.
  epsilon <- mu / (1 - mu)
  Tau <- sigma2 * ((1 - mu)^2)

  # The rMIG function is now able to receive vectors for 'epsilon', 'Tau', and 'mu'.
  # This is handled by the next step.
  x <- rMIG(n, epsilon, Tau, mu)

  # The final transformation is also vectorized.
  yy <- x / (1 + x)

  return(as.vector(yy))
}


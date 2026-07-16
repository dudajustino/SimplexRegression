#' Biomass Allocation in Two Grass Species Under Different Nitrate Supply
#'
#' @description
#' This dataset examines biomass allocation patterns in plants, specifically
#' the proportional distribution of biomass to different plant organs (stems,
#' leaves, and roots). The data come from an experiment manipulating nitrate
#' supply in fast-growing and slow-growing grass species.
#'
#' The response variables (stem, leaves, roots) are proportions bounded in the
#' (0, 1) interval, making them suitable for simplex regression analysis.
#'
#' @docType data
#' @usage data(Biomass)
#'
#' @format A data frame with 500 observations and 13 variables:
#' \describe{
#'   \item{group}{Factor. Combined species-by-nitrate-treatment code
#'         (\code{DfH}, \code{DfL}, \code{HlH}, \code{HlL}), corresponding to the combination of
#'         \code{species} and \code{trt} where:
#'     \itemize{
#'       \item \code{DfH} = D. flexuosa (slow-growing), high nitrate
#'       \item \code{DfL} = D. flexuosa (slow-growing), low nitrate
#'       \item \code{HlH} = H. lanatus (fast-growing), high nitrate
#'       \item \code{HlL} = H. lanatus (fast-growing), low nitrate
#'     }
#'   }
#'   \item{species}{Factor. Species scientific name (\code{D. flexuosa} or
#'         \code{H. lanatus}).}
#'   \item{trt}{Factor. Nitrate treatment level (\code{high} or \code{low}).}
#'   \item{day}{Numeric. Experimental day (0 to 49).}
#'   \item{pl_num}{Numeric. Plant number (individual plant identifier,
#'         6-8 replicates per treatment).}
#'   \item{ldm_mg}{Numeric. Leaf dry mass (in mg).}
#'   \item{sdm_mg}{Numeric. Stem dry mass (in mg).}
#'   \item{rdm_mg}{Numeric. Root dry mass (in mg).}
#'   \item{tdm_mg}{Numeric. Total dry mass (in mg).}
#'   \item{lmf}{Numeric. Leaf mass fraction, proportion of biomass allocated
#'         to leaves, bounded in the open interval (0, 1).}
#'   \item{smf}{Numeric. Stem mass fraction, proportion of biomass allocated
#'         to stems, bounded in the open interval (0, 1).}
#'   \item{rmf}{Numeric. Root mass fraction, proportion of biomass allocated
#'         to roots, bounded in the open interval (0, 1).}
#'   \item{ln_tdm}{Numeric. Natural log of total dry mass (log-transformed
#'         for allometric analysis).}
#' }
#'
#' @source
#' bobdouma (2019). \emph{bobdouma/proportions_beta_Dirichlet: v.01}.
#' \doi{10.5281/zenodo.3234670}.
#'
#' @references
#' bobdouma (2019). \emph{bobdouma/proportions_beta_Dirichlet: v.01}.
#' \doi{10.5281/zenodo.3234670}.
#'
#' Poorter, H.; van de Vijver, C. A. D. M.; Boot, R. G. A. and Lambers, H. (1995).
#' Growth and carbon economy of a fast-growing and a slow-growing grass
#' species as dependent on nitrate supply. \emph{Plant and Soil}, \bold{171},
#' 217--227. \doi{10.1007/BF00010275}
#'
#' Poorter, H. and Sack, L. (2012). Pitfalls and possibilities in the analysis
#' of biomass allocation patterns in plants. \emph{Frontiers in Plant
#' Science}, \bold{3}, 259. \doi{10.3389/fpls.2012.00259}
#'
#' @keywords datasets
#'
#' @examples
#' # Load the data
#' data(Biomass)
#'
#' # Quick overview
#' head(Biomass)
#' str(Biomass)
#'
#' # Check that proportions sum to 1 (within rounding error)
#' summary(rowSums(Biomass[, c("lmf", "smf", "rmf")]))
#'
#' # Simple plot of root mass fraction by treatment
#' boxplot(rmf ~ trt * species, data = Biomass,
#'         main = "Root Mass Fraction by Species and Nitrate Treatment",
#'         las = 2)
"Biomass"

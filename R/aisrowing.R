#' Body Composition Data for Australian Rowers
#'
#' @description
#' Hematological, body composition, and anthropometric measurements on
#' 37 elite rowers (22 female, 15 male) at the Australian Institute of
#' Sport (AIS), a subset of a larger dataset collected across multiple
#' sports. The data are useful for investigating sex-based differences
#' in blood and body composition among highly trained athletes.
#'
#' @docType data
#' @usage data(AISRowing)
#'
#' @format A data frame with 37 observations and 12 variables:
#' \describe{
#'   \item{sex}{Factor. Sex of the athlete (\code{female} or \code{male}).}
#'   \item{rcc}{Numeric. Red blood cell count (in \eqn{10^{12}} per litre).}
#'   \item{wcc}{Numeric. White blood cell count (in \eqn{10^{12}} per litre).}
#'   \item{hc}{Numeric. Hematocrit (in percent).}
#'   \item{hb}{Numeric. Hemoglobin concentration (in g per decilitre).}
#'   \item{ferr}{Numeric. Plasma ferritin concentration (in ng per millilitre).}
#'   \item{bmi}{Numeric. Body mass index (in kg per metre-squared).}
#'   \item{ssf}{Numeric. Sum of skin folds (in mm).}
#'   \item{bfat}{Numeric. Body fat proportion, originally measured
#'         in percent and rescaled to the unit interval (0, 1).}
#'   \item{lbm}{Numeric. Lean body mass (in kg).}
#'   \item{ht}{Numeric. Height (in cm).}
#'   \item{wt}{Numeric. Weight (in kg).}
#' }
#'
#' @details
#' The original measurements were collected in the late 1980s by Richard
#' Telford and Ross Cunningham at the Australian Institute of Sport (AIS),
#' and were later compiled and popularized by Cook and Weisberg (1994).
#' \code{AISRowing} is the 37-athlete subset corresponding to rowing,
#' out of 202 athletes across all sports in the original dataset. The
#' full dataset (covering all sports, not just rowing) is also discussed
#' in Weisberg (2005, Section 6.4).
#'
#' @source
#' Telford, R. D. and Cunningham, R. B. (1991). Sex, sport, and body-size
#' dependency of hematology in highly trained athletes. \emph{Medicine &
#' Science in Sports & Exercise}, \bold{23}(7), 788--794.
#'
#' @references
#' Cook, R. D. and Weisberg, S. (1994). \emph{An Introduction to Regression
#' Graphics}. John Wiley & Sons, New York.
#'
#' Weisberg, S. (2005). \emph{Applied Linear Regression}, 3rd edition.
#' New York: Wiley, Section 6.4.
#'
#' @keywords datasets
#'
#' @examples
#' # Load the data
#' data(AISRowing)
#'
#' # Quick overview
#' head(AISRowing)
#' str(AISRowing)
#'
#' # Summary statistics
#' summary(AISRowing)
#'
#' # Sex-based comparison
#' aggregate(bfat ~ sex, data = AISRowing, summary)
#' boxplot(bfat ~ sex, data = AISRowing,
#'         main = "Body Fat Proportion by Sex",
#'         ylab = "Body fat (proportion)")
"AISRowing"

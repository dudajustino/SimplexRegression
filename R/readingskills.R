#' Reading Accuracy in Dyslexic and Non-Dyslexic Children
#'
#' @description
#' This dataset examines the relationship between non-verbal IQ and reading
#' accuracy in children diagnosed with dyslexia and in typical readers.
#' The reading scores are proportions bounded in the (0,1) interval, making
#' them suitable for beta regression and simplex regression analysis.
#'
#' @docType data
#' @usage data(ReadingSkills)
#'
#' @format A data frame with 44 observations and 4 variables:
#' \describe{
#'   \item{accuracy}{Numeric. Reading score transformed to the open (0,1) interval.
#'         Values originally equal to 1 were replaced with 0.99. Suitable for
#'         standard beta regression.}
#'   \item{dyslexia}{Factor. Indicates whether the child has dyslexia
#'         (levels: "no", "yes"). Note that sum contrasts are typically used
#'         instead of treatment contrasts in beta regression analyses of this data.}
#'   \item{iq}{Numeric. Non-verbal intelligence quotient, transformed to
#'         z-scores (mean = 0, standard deviation = 1).}
#'   \item{accuracy1}{Numeric. Unrestricted reading score in the [0,1] interval.
#'         This version preserves the original maximum value of 1 and can be used
#'         with extended-support beta mixture regression models.}
#' }
#'
#' @details
#' The data were originally collected by Pammer and Kevan (2004) and later
#' analyzed by Smithson and Verkuilen (2006) to demonstrate beta regression.
#'
#' The transformation procedure for \code{accuracy} was as follows:
#' \enumerate{
#'   \item The original test scores were scaled to the [0,1] interval using
#'         the minimum and maximum possible scores in the reading test,
#'         resulting in \code{accuracy1}.
#'   \item To avoid boundary values (0 and 1) that are problematic for standard
#'         beta regression, all observations with value 1 were replaced with 0.99,
#'         creating the \code{accuracy} variable.
#' }
#'
#' The unrestricted \code{accuracy1} variable can be analyzed using extended-support
#' beta regression methods (Kosmidis & Zeileis, 2025), which naturally accommodate
#' boundary observations.
#'
#' @source
#' Pammer, K. & Kevan, A. (2004). The Contribution of Visual Sensitivity,
#' Phonological Processing and Non-Verbal IQ to Children's Reading.
#' \emph{Unpublished manuscript}, The Australian National University, Canberra.
#'
#' Smithson, M. & Verkuilen, J. (2006). A Better Lemon Squeezer?
#' Maximum-Likelihood Regression with Beta-Distributed Dependent Variables.
#' \emph{Psychological Methods}, 11(1), 54-71.
#'
#' @references
#' Cribari-Neto, F. & Zeileis, A. (2010). Beta Regression in R.
#' \emph{Journal of Statistical Software}, 34(2), 1-24.
#' \doi{10.18637/jss.v034.i02}
#'
#' Kosmidis, I. & Zeileis, A. (2025). Extended-Support Beta Regression for
#' [0, 1] Responses. \emph{Journal of the Royal Statistical Society C},
#' forthcoming. \doi{10.1093/jrsssc/qlaf039}
#'
#' @keywords datasets
#'
#' @examples
#' # Load the data
#' data(ReadingSkills)
#'
#' # Quick overview
#' head(ReadingSkills)
#' str(ReadingSkills)
#'
#' # Summary statistics by dyslexia status
#' aggregate(accuracy ~ dyslexia, data = ReadingSkills, summary)
#' aggregate(iq ~ dyslexia, data = ReadingSkills, summary)
#'
#' # Visualize the relationship between IQ and reading accuracy
#' plot(accuracy ~ iq, data = ReadingSkills,
#'      col = c(4, 2)[dyslexia], pch = 19,
#'      main = "Reading Accuracy vs. IQ by Dyslexia Status",
#'      xlab = "IQ (z-scored)", ylab = "Reading Accuracy")
#' legend("topleft", legend = c("Non-dyslexic", "Dyslexic"),
#'        col = c(4, 2), pch = 19, bty = "n")
#'
#' # Check for boundary values
#' table(ReadingSkills$accuracy == 0.99)  # Values replaced from 1
#' table(ReadingSkills$accuracy1 == 1)    # Original boundary values
"ReadingSkills"

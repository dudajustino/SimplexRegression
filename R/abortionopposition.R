#' Public Opinion on Abortion Across U.S. States
#'
#' @description
#' This dataset examines the relationship between public opposition to abortion
#' and demographic, religious, and socioeconomic characteristics across the
#' 50 U.S. states and the District of Columbia. The response variable is the
#' proportion of adults who believe abortion should be illegal in all or most
#' cases, bounded in the open interval (0, 1), making the data suitable for simplex
#' regression and beta regression analyses.
#'
#' @docType data
#' @usage data(AbortionOpposition)
#'
#' @format A data frame with 51 observations and 9 variables:
#' \describe{
#'   \item{state}{Factor. U.S. state or the District of Columbia.}
#'   \item{abortion_opp}{Numeric. Proportion of adults who believe
#'         abortion should be illegal in all or most cases, based on the Pew
#'         Research Center 2023--24 Religious Landscape Study.}
#'   \item{relig_attend}{Numeric. Percentage of adults who report
#'         attending religious services at least once a week, based on the
#'         Pew Research Center 2023--24 Religious Landscape Study.}
#'   \item{income}{Numeric. Mean annual household income (USD), 2024.}
#'   \item{sex_ratio}{Numeric. Number of males per 100 females, 2024.}
#'   \item{female}{Numeric. Percentage of the population that is female, 2024.}
#'   \item{pop_18_24}{Numeric. Percentage of the population aged 18--24 years,
#'         obtained directly from the 2024 American Community Survey (ACS)
#'         1-Year Estimates.}
#'   \item{bachelors}{Numeric. Percentage of adults aged 25 and older with a
#'         bachelor's degree or higher, 2024.}
#'   \item{urban}{Numeric. Percentage of the population living in urban areas,
#'         based on 2020 U.S. Census Bureau urban area criteria (a densely
#'         settled core of census blocks meeting minimum housing unit and/or
#'         population density thresholds, requiring at least 2,000 housing
#'         units or 5,000 people).}
#' }
#'
#' @details
#' The dataset was assembled by the package authors from multiple publicly
#' available sources. The response variable, \code{abortion_opp},
#' was obtained from the Pew Research Center 2023--24 Religious Landscape
#' Study, originally published as a percentage (0--100 scale), and rescaled
#' to the (0, 1) interval by dividing by 100.
#'
#' The explanatory variables were obtained as follows:
#' \enumerate{
#'   \item \code{relig_attend} was obtained from the Pew Research Center
#'         2023--24 Religious Landscape Study.
#'   \item \code{income}, \code{sex_ratio}, and \code{bachelors} were obtained from
#'         World Population Review, which itself compiles estimates from the
#'         U.S. Census Bureau/American Community Survey. World Population
#'         Review is cited here as the immediate data source, following
#'         standard practice for reproducibility.
#'   \item \code{female} was obtained from the U.S. Census Bureau's American
#'         Community Survey (ACS) 1-Year Estimates, Table S0101 (Age and
#'         Sex), the same table used for \code{pop_18_24}. Note that
#'         \code{sex_ratio} and \code{female} are closely related (both
#'         describe the sex composition of each state's population, from
#'         different sources) and should generally not be included
#'         together as covariates in the same regression model, as doing
#'         so may introduce near-perfect collinearity.
#'   \item \code{pop_18_24} corresponds to the percent estimate of the
#'         population aged 18--24 years reported directly in Table S0101
#'         (Age and Sex) of the 2024 American Community Survey (ACS)
#'         1-Year Estimates.
#'   \item \code{urban} was obtained from World Population Review and is
#'         based on 2020 U.S. Census Bureau data, the most recent Census
#'         estimate of urbanization by state available at the time this
#'         dataset was compiled. Unlike the other explanatory variables,
#'         which reflect 2023--24 estimates, \code{urban} reflects the 2020
#'         Census urban area delineation (see \code{@format} above for
#'         definition details).
#' }
#'
#' @source
#' Pew Research Center (2025). 2023-24 Religious Landscape Study (RLS)
#' Dataset. \doi{10.58094/3kwb-bf52}. Data accessed in 2026.
#'
#' World Population Review (2024). Per Capita Income by State.
#' \url{https://worldpopulationreview.com/state-rankings/per-capita-income-by-state}.
#' Data accessed in 2026.
#'
#' World Population Review (2024). Sex Ratio by State.
#' \url{https://worldpopulationreview.com/state-rankings/sex-ratio-by-state}.
#' Data accessed in 2026.
#'
#' World Population Review (2024). Educational Attainment by State.
#' \url{https://worldpopulationreview.com/state-rankings/educational-attainment-by-state}.
#' Data accessed in 2026.
#'
#' World Population Review (2024). Most Urban States.
#' \url{https://worldpopulationreview.com/state-rankings/most-urban-states}.
#' Data accessed in 2026.
#'
#' U.S. Census Bureau (2024). Age and Sex. American Community Survey,
#' ACS 1-Year Estimates Subject Tables, Table S0101.
#' \url{https://data.census.gov/table/ACSST1Y2024.S0101?q=age+and+sex+by+state&moe=false}.
#' Data accessed in 2026.
#'
#' @keywords datasets
#'
#' @examples
#' # Load the data
#' data(AbortionOpposition)
#'
#' # Quick overview
#' head(AbortionOpposition)
#' str(AbortionOpposition)
#'
#' # Summary statistics
#' summary(AbortionOpposition)
#'
#' # Relationship between religious attendance and abortion opposition
#' plot(abortion_opp ~ relig_attend,
#'      data = AbortionOpposition,
#'      pch = 19,
#'      xlab = "Religious attendance (%)",
#'      ylab = "Abortion opposition (proportion)")
#'
#' # Correlation among numeric variables
#' cor(AbortionOpposition[, c("abortion_opp",
#'                         "relig_attend",
#'                         "income",
#'                         "sex_ratio",
#'                         "female",
#'                         "pop_18_24",
#'                         "bachelors",
#'                         "urban")])
"AbortionOpposition"

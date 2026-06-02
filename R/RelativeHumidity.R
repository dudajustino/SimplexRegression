#' Monthly Relative Humidity Data
#'
#' @description
#' Monthly meteorological data including relative humidity (RH),
#' temperature, insolation, cloudiness, precipitation, atmospheric pressure,
#' and wind speed for a Brazilian station (2000-2025).
#'
#' The dataset contains two missing observations in the variables \code{Ins}
#' and \code{Pre}, corresponding to September 2020 and November 2025. To address
#' these missing values, imputation was performed via seasonal interpolation using
#' the \pkg{imputeTS} package, with the imputed versions stored as \code{Ins2}
#' and \code{Pre2}.
#'
#' @format A data frame with 312 observations and 11 variables:
#' \describe{
#'   \item{Date}{Observation date (YYYY-MM-DD)}
#'   \item{RH}{Monthly mean relative humidity, originally measured in percent and
#'   rescaled to the unit interval (0,1)}
#'   \item{Ins}{Monthly total insolation (in hours), contains 2 missing values}
#'   \item{Ins2}{Monthly total insolation (in hours) with missing values imputed
#'   via seasonal interpolation}
#'   \item{Pre}{Monthly total precipitation (in mm), contains 2 missing values}
#'   \item{Pre2}{Monthly total precipitation (in mm) with missing values imputed
#'   via seasonal interpolation}
#'   \item{Neb}{Monthly mean cloudiness (in tenths)}
#'   \item{AP}{Monthly mean atmospheric pressure (in hPa)}
#'   \item{MT}{Monthly mean temperature (in Degrees Celsius)}
#'   \item{WS}{Monthly mean wind speed (in m/s)}
#'   \item{Dir}{Monthly predominant wind direction (in degrees)}
#' }
#' @source Brazilian meteorological station
"RelativeHumidity"

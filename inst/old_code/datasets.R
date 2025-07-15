#' Sounding data over Boston
#'
#' A dataset containing measurements of solar-induced cholophyl flourescence (SIF)
#' and albedo within polygons over Boston.
#'
#'
#' \itemize{
#'   \item Sounding. An integer denoting the order of soundings during retrieval (NEEDS DOUBLE CHECKED)
#'   \item Timestamp. A POSIXcd datetime.
#'   \item MeasurementMode
#'   \item ViewZenith
#'   \item ViewAzimuth
#'   \item SIF_757nm
#'   \item SIF_Uncertainty_757nm
#'   \item SIF_771nm
#'   \item SIF_Uncertainty_771nm
#'   \item Albedo
#'   \item cloud_flag_abp
#'   \item Quality_Flag
#'   \item sounding_land_fraction
#'   \item sounding_qual_flag
#'   \item geometry
#' }
#'
#' @docType data
#' @keywords datasets
#' @name soundings
#' @usage data(soundings)
#' @format An sf object with 6854 rows and 15 variables.
NULL

#' A regular grid overlapping with the sounding data over boston
#' (see `soundings`).
#'
#' A collection of polygons alongside corresponding measurments of albedo data and
#'
#'
#' \itemize{
#'   \item row. An integer denoting the row in the grid.
#'   \item col. An integer denoting the column in the grid.
#'   \item MeasurementMode
#'   \item ViewZenith
#'   \item ViewAzimuth
#'   \item SIF_757nm
#'   \item SIF_Uncertainty_757nm
#'   \item SIF_771nm
#'   \item SIF_Uncertainty_771nm
#'   \item Albedo
#'   \item cloud_flag_abp
#'   \item Quality_Flag
#'   \item sounding_land_fraction
#'   \item sounding_qual_flag
#'   \item geometry
#' }
#'
#' @docType data
#' @keywords datasets
#' @name target_grid
#' @usage data(target_grid)
#' @format An sf object with 8290 rows and 20 columns
NULL

#' A leaflet color palette that is useful for visualizing SIF values.
#'
#'
#'
#'
#' @docType data
#' @keywords color, palette
#' @name color_palette
#' @usage data(target_grid)
#' @format A numeric color palette from the leaflet package
NULL


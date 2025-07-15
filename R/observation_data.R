#' Create Observation Data Object for Downscaling
#'
#' Extracts the response vector and builds region representations
#' from an sf object, returning a structured object for modeling.
#'
#' @param sf An sf object with geometries and the response variable.
#' @param response Name of the response column (character).
#' @param region_method One of "quads", "triangles", or "auto". Determines how to build the region list.
#' @param ... Additional arguments passed to `build_region_list()`.
#'
#' @return An object of class "observation_data" with slots:
#'   - y: numeric vector of responses
#'   - region_list: list of region objects
#'   - sf_metadata: list with CRS and row IDs
#'   - sf_original: the input sf object (optional, can remove if you prefer)
#'
#' @export
make_observation_data <- function(
    sf,
    response,
    region_method = "raw",
    ...
) {
  stopifnot(inherits(sf, "sf"))
  stopifnot(response %in% names(sf))

  # Extract y
  y <- sf[[response]]
  if (!is.numeric(y)) {
    stop("The response variable must be numeric.")
  }

  # Build regions
  region_list <- build_region_list(sf, method = region_method, ...)

  # Create object
  out <- list(
    y = y,
    region_list = region_list,
    sf_metadata = list(
      crs = sf::st_crs(sf),
      n_obs = nrow(sf)
    ),
    sf_original = sf # optional, could omit to save memory
  )
  class(out) <- "observation_data"

  return(out)
}

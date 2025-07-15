#' Ensure an sf object is in a projected CRS using UTM zone if possible
#'
#' Converts an sf object to a projected coordinate reference system (CRS) if it is
#' currently in a geographic CRS (e.g., EPSG:4326). Automatically selects an appropriate
#' UTM zone based on the centroid of the geometry. If the CRS is already projected,
#' the object is returned unchanged.
#'
#' @param x An `sf` or `sfc` object.
#' @param fallback_epsg Integer EPSG code to use if UTM selection fails. Defaults to EPSG:6933 (Equal Earth).
#' @param verbose Logical; if `TRUE`, print messages when conversion occurs.
#'
#' @return An `sf` or `sfc` object transformed to a projected CRS if needed.
#'
#' @examples
#' x <- sf::st_as_sf(data.frame(lon = -120, lat = 35), coords = c("lon", "lat"), crs = 4326)
#' x_proj <- ensure_projected(x)
#'
#' @export
ensure_projected <- function(x, fallback_epsg = 6933, verbose = TRUE) {
  if (!inherits(x, c("sf", "sfc"))) {
    stop("`x` must be an object of class `sf` or `sfc`.")
  }
  crs <- sf::st_crs(x)
  if (is.na(crs)) {
    stop("Object has no CRS. Please assign one before using `ensure_projected()`.")
  }
  if (!crs$IsGeographic) {
    # Already projected
    return(x)
  }

  # Try to compute centroid coordinates
  centroid <- tryCatch(
    {
      cxy <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(x)))[1, ]
      list(lon = cxy[1], lat = cxy[2])
    },
    error = function(e) NULL
  )

  if (is.null(centroid)) {
    if (verbose) {
      message("Could not compute centroid; falling back to EPSG:", fallback_epsg)
    }
    return(sf::st_transform(x, fallback_epsg))
  }

  # Compute UTM zone
  lon <- centroid$lon
  lat <- centroid$lat
  utm_zone <- floor((lon + 180) / 6) + 1
  if (lat >= 0) {
    epsg <- 32600 + utm_zone  # Northern hemisphere
  } else {
    epsg <- 32700 + utm_zone  # Southern hemisphere
  }

  if (verbose) {
    message(
      sprintf(
        "Transforming object from geographic CRS (EPSG:%s) to UTM zone %d (EPSG:%d).",
        crs$epsg,
        utm_zone,
        epsg
      )
    )
  }

  # Transform
  x_proj <- tryCatch(
    sf::st_transform(x, epsg),
    error = function(e) {
      if (verbose) {
        message("UTM transformation failed; falling back to EPSG:", fallback_epsg)
      }
      sf::st_transform(x, fallback_epsg)
    }
  )

  return(x_proj)
}


#' Ensure Two CRS Objects Are Identical
#'
#' Checks whether two CRS definitions (already computed via `sf::st_crs()`)
#' are identical. If not, throws an error with contextual messaging.
#'
#' @param crs_a A coordinate reference system object created by `sf::st_crs()`.
#' @param crs_b Another CRS object created by `sf::st_crs()`.
#' @param context Optional string describing the context of the check (e.g., "integration", "evaluation").
#'
#' @return Invisibly returns `TRUE` if the CRSs match. Otherwise, throws an error.
#' @export
#'
#' @examples
#' crs1 <- sf::st_crs(4326)
#' crs2 <- sf::st_crs("EPSG:4326")
#' ensure_same_crs(crs1, crs2, context = "prediction")  # Passes
#'
#' \dontrun{
#' crs3 <- sf::st_crs(3857)
#' ensure_same_crs(crs1, crs3, context = "integration")  # Errors
#' }
  ensure_same_crs <- function(crs_a, crs_b, context = "operation") {
  if (!inherits(crs_a, "crs") || !inherits(crs_b, "crs")) {
    stop("Inputs to `ensure_same_crs()` must be `sf::st_crs()` objects.")
  }

  if (!isTRUE(crs_a == crs_b)) {
    stop(sprintf(
      "CRS mismatch during %s:\n  - CRS A: %s\n  - CRS B: %s\nPlease transform inputs to a common CRS before proceeding.",
      context,
      crs_a$input %||% "<unknown>",
      crs_b$input %||% "<unknown>"
    ))
  }

  invisible(TRUE)
}

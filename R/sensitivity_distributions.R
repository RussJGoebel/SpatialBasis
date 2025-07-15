#' Gaussian Sensitivity Weight Dispatcher
#'
#' Computes Gaussian weights for a set of coordinates given precomputed parameters
#' specific to each observation. This function is designed to replace a list of closures
#' and instead use a single shared function plus numeric parameter data for efficiency
#' and better parallel performance.
#'
#' @param coords A numeric matrix with two columns (x, y), representing the mapped
#'   integration points in projected coordinates.
#' @param index An integer index indicating which observation's parameters to use.
#' @param params A list of parameter lists, one per observation. Each parameter list
#'   must contain:
#'   - `x0`: numeric center x-coordinate
#'   - `y0`: numeric center y-coordinate
#'   - `angle`: numeric rotation angle in radians
#'   - `sig_x`: numeric variance parameter along the major axis
#'   - `sig_y`: numeric variance parameter along the minor axis
#'
#' @return A numeric vector of length equal to nrow(coords) containing normalized weights.
#'
#' @export
gaussian_sensitivity_dispatcher <- function(coords, index, params) {
  p <- params[[index]]
  delta_x <- coords[, 1] - p$x0
  delta_y <- coords[, 2] - p$y0
  d <- sqrt(delta_x^2 + delta_y^2)
  ang <- atan2(delta_y, delta_x)
  dx <- d * cos(ang - p$angle)
  dy <- d * sin(ang - p$angle)
  w <- exp(-dx^2 / p$sig_x) * exp(-dy^2 / p$sig_y)
  w / sum(w)
}

#' Generate Gaussian Sensitivity Parameters for Each Sounding
#'
#' Computes the center coordinates, rotation angle, and variance parameters
#' for each observation polygon. These parameters are used in the
#' `gaussian_sensitivity_dispatcher()` for efficient weight computation.
#'
#' @param soundings_proj An sf object with geometries in projected coordinates (e.g., meters).
#'   Must be a collection of polygons (e.g., footprints).
#' @param major_axis_sigma_scaler Numeric scalar multiplier for the variance along the major axis.
#' @param minor_axis_sigma_scaler Numeric scalar multiplier for the variance along the minor axis.
#'
#' @return A list of length equal to the number of soundings, each containing:
#'   - `x0`: numeric center x-coordinate
#'   - `y0`: numeric center y-coordinate
#'   - `angle`: numeric rotation angle (radians)
#'   - `sig_x`: numeric variance parameter along the major axis
#'   - `sig_y`: numeric variance parameter along the minor axis
#'
#' @export
generate_gaussian_sensitivity_params <- function(
    soundings_proj,
    major_axis_sigma_scaler = 1,
    minor_axis_sigma_scaler = 1
) {
  n_obs <- length(soundings_proj$geometry)

  params <- vector("list", n_obs)

  for (i in seq_len(n_obs)) {
    # Extract the three defining points
    coords <- sf::st_coordinates(soundings_proj$geometry[[i]])[1:3, ]

    x1 <- coords[1, 1]
    y1 <- coords[1, 2]
    x2 <- coords[2, 1]
    y2 <- coords[2, 2]
    x3 <- coords[3, 1]
    y3 <- coords[3, 2]

    # Compute axes
    minor_axis_length <- sqrt((x2 - x1)^2 + (y2 - y1)^2) / 2
    major_axis_length <- sqrt((x3 - x2)^2 + (y3 - y2)^2) / 2
    major_axis_angle  <- atan2(y3 - y2, x3 - x2)

    # Compute variances
    sig_x <- 2 * (major_axis_length * major_axis_sigma_scaler)^2
    sig_y <- 2 * (minor_axis_length * minor_axis_sigma_scaler)^2

    # Center
    x0 <- mean(c(x1, x3))
    y0 <- mean(c(y1, y3))

    # Store
    params[[i]] <- list(
      x0 = x0,
      y0 = y0,
      angle = major_axis_angle,
      sig_x = sig_x,
      sig_y = sig_y
    )
  }

  return(params)
}

#' Scale every geometry in an sf object around its centroid
#'
#' @param sf_obj An sf object with one or more rows
#' @param scale_factor Numeric scaling factor
#'
#' @return An sf object with all geometries scaled
#'
#' @export
#'
#' @examples
#' library(sf)
#' # Create example sf object with 2 rectangles
#' rects <- st_sf(
#'   id = 1:2,
#'   geometry = st_sfc(
#'     st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0)))),
#'     st_polygon(list(rbind(c(2,2), c(3,2), c(3,3), c(2,3), c(2,2))))
#'   ),
#'   crs = 4326
#' )
#' # Scale all geometries by 2
#' rects_scaled <- scale_sf_geometries(rects, 2)
#' # Plot
#' plot(st_geometry(rects), border = "blue")
#' plot(st_geometry(rects_scaled), border = "red", add = TRUE)
scale_sf_geometries <- function(sf_obj, scale_factor) {
  stopifnot(inherits(sf_obj, "sf"))

  # Extract geometry
  geom <- sf::st_geometry(sf_obj)

  # Scale each feature
  scaled_list <- lapply(seq_along(geom), function(i) {
    g <- geom[[i]]            # sfg geometry
    centroid <- sf::st_centroid(sf::st_sfc(g))
    coords <- sf::st_coordinates(centroid)

    # Center the geometry relative to the centroid
    g_centered <- g
    g_centered[[1]][,1] <- g_centered[[1]][,1] - coords[1]
    g_centered[[1]][,2] <- g_centered[[1]][,2] - coords[2]

    # Scale
    g_scaled <- g_centered
    g_scaled[[1]][,1] <- g_scaled[[1]][,1] * scale_factor
    g_scaled[[1]][,2] <- g_scaled[[1]][,2] * scale_factor

    # Re-translate back to centroid
    g_scaled[[1]][,1] <- g_scaled[[1]][,1] + coords[1]
    g_scaled[[1]][,2] <- g_scaled[[1]][,2] + coords[2]

    g_scaled
  })

  # Reassemble
  scaled_sfc <- sf::st_sfc(scaled_list, crs = sf::st_crs(geom))

  # Replace geometry
  sf_scaled <- sf_obj
  sf::st_geometry(sf_scaled) <- scaled_sfc

  return(sf_scaled)
}

#' Scale a single sfc geometry around its centroid
#'
#' @param geom An sfc geometry of length 1
#' @param scale_factor Numeric scaling factor
#'
#' @return A scaled sfc geometry
scale_single_geometry <- function(geom, scale_factor) {
  stopifnot(inherits(geom, "sfc"), length(geom) == 1)
  centroid <- st_centroid(geom)
  scaled_geom <- (geom - centroid) * scale_factor + centroid
  scaled_geom
}

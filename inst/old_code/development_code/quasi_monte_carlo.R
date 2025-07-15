#' Generate QMC Points in [0,1]^2
#'
#' Generates low-discrepancy quasi-Monte Carlo points using the Sobol sequence.
#'
#' @param n Number of points to generate.
#'
#'
#' @importFrom randtoolbox sobol
#' @return A matrix of size [n × 2] with points in [0,1]^2.
#' @export
generate_qmc_unit_square <- function(n) {

  randtoolbox::sobol(
    n = n,
    dim = 2,
    scrambling = 0,
    normal = FALSE
  )

}




#' Map points from unit square to a quadrilateral
#'
#' Transforms points from the unit square \eqn{[0,1]^2} to a general quadrilateral using bilinear interpolation.
#' The quadrilateral is defined by four corner points, ordered clockwise or counter-clockwise as:
#' bottom-left, bottom-right, top-right, top-left.
#'
#' This function is typically used to warp a quasi-Monte Carlo sample from the unit square
#' to an arbitrary 2D quadralateral (especially for integration or sampling tasks).
#'
#' @param qmc_points An \eqn{n × 2} matrix of points in \eqn{[0,1]^2}.
#' @param quad_coords A \eqn{4 × 2} matrix giving the coordinates of the quadrilateral's four corners.
#'        The corners must be ordered as: bottom-left, bottom-right, top-right, top-left.
#'
#' @return An \eqn{n × 2} matrix of transformed coordinates lying inside the quadrilateral.
#'
#' @examples
#' qmc <- generate_qmc_unit_square(100)
#' quad <- matrix(c(0,0, 1,0, 1,1, 0,1), ncol = 2, byrow = TRUE)
#' mapped <- map_unit_square_to_quad(qmc, quad)
#'
#' @export
map_unit_square_to_quad <- function(qmc_points, quad_coords) {
  # quad_coords: 4x2 matrix, assumed ordered as (bottom-left, bottom-right, top-right, top-left)
  stopifnot(nrow(quad_coords) == 4)

  # Extract (u,v) from unit square
  u <- qmc_points[, 1]
  v <- qmc_points[, 2]

  # Bilinear interpolation
  mapped <- (1 - u) * (1 - v) %*% quad_coords[1, , drop = FALSE] +
    u     * (1 - v) %*% quad_coords[2, , drop = FALSE] +
    u     *     v   %*% quad_coords[3, , drop = FALSE] +
    (1 - u) *     v   %*% quad_coords[4, , drop = FALSE]

  return(mapped)
}

#' Map unit square points to a unit triangle
#'
#' Transforms points from the unit square \eqn{[0,1]^2} to a triangle using the inverse transform method
#' that produces a uniform distribution over the triangle.
#'
#' This is useful when using quasi-Monte Carlo points (e.g., Sobol) to sample over a triangle.
#' The output triangle is defined by three vertices (2D coordinates), ordered arbitrarily.
#'
#' @param qmc_points An \eqn{n × 2} matrix of points in \eqn{[0,1]^2}.
#' @param tri_coords A \eqn{3 × 2} matrix giving the triangle's corner coordinates.
#'
#' @return An \eqn{n × 2} matrix of points uniformly mapped inside the triangle.
#'
#' @examples
#' qmc <- generate_qmc_unit_square(100)
#' tri <- matrix(c(0,0, 1,0, 0,1), ncol = 2, byrow = TRUE)
#' mapped <- map_unit_square_to_triangle(qmc, tri)
#'
#' @export
map_unit_square_to_triangle <- function(qmc_points, tri_coords) {
  stopifnot(nrow(tri_coords) == 3)

  u <- qmc_points[, 1]
  v <- qmc_points[, 2]

  # Apply the standard "flip" transform to ensure uniformity over triangle
  flip <- u + v > 1
  u[flip] <- 1 - u[flip]
  v[flip] <- 1 - v[flip]

  A <- tri_coords[1, ]
  B <- tri_coords[2, ]
  C <- tri_coords[3, ]

  # Affine combination: (1 - u - v) * A + u * B + v * C
  mapped <- (1 - u - v) %*% matrix(A, nrow = 1) +
    u %*% matrix(B, nrow = 1) +
    v %*% matrix(C, nrow = 1)

  return(mapped)
}




#' Extract coordinate matrices from an sf POLYGON/MULTIPOLYGON object
#'
#' Converts each geometry in the sf object into a matrix of polygon coordinates.
#' Assumes each geometry has a single outer ring. For quadrilateral-specific
#' routines, the first 4 unique coordinates are returned.
#'
#' @param sf_polygons An sf object with POLYGON or MULTIPOLYGON geometries.
#' @param expected_vertices Optional: expected number of unique vertices per polygon (default = 4).
#'
#' @return A list of matrices, each with dimensions [k × 2], where k is the number of unique vertices.
#' @export
extract_polygon_coords_list <- function(sf_polygons, expected_vertices = 4) {
  if (!inherits(sf_polygons, "sf")) {
    stop("Input must be an sf object.")
  }

  geom_type <- unique(sf::st_geometry_type(sf_polygons))
  if (!all(geom_type %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("Geometry must be POLYGON or MULTIPOLYGON.")
  }

  geom_list <- sf::st_geometry(sf_polygons)
  coords_list <- vector("list", length(geom_list))

  for (i in seq_along(geom_list)) {
    coords <- sf::st_coordinates(geom_list[[i]])
    unique_coords <- coords[!duplicated(coords[, c("X", "Y")]), c("X", "Y")]

    if (!is.null(expected_vertices) && nrow(unique_coords) != expected_vertices) {
      warning(sprintf("Polygon %d has %d unique vertices (expected %d).",
                      i, nrow(unique_coords), expected_vertices))
    }

    coords_list[[i]] <- as.matrix(unique_coords)
  }


  return(coords_list)

}




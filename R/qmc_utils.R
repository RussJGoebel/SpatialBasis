#' Generate Sobol QMC Points in the Unit Square
#'
#' Returns a matrix of quasi-Monte Carlo (QMC) points in the unit square \code{[0,1]^2}
#' using a Sobol sequence, suitable for mapping into triangle or quad regions.
#'
#' @param n_points Number of QMC points to generate
#' @return An n_points × 2 matrix of QMC points
#'
#' @examples
#' qmc <- generate_qmc_unit_square(128)
#' plot(qmc, asp = 1, main = "Sobol QMC Points in [0,1]^2")
#'
#' @export
generate_qmc_unit_square <- function(n_points) {
  qrng::sobol(n = n_points, d = 2)
}

#' Map unit square points to a triangle using barycentric coordinates
#'
#' Maps an n × 2 matrix of points from the unit square \code{[0,1]^2} into a triangle
#' defined by 3 coordinates using a uniform barycentric transformation.
#'
#' @param qmc_points An n × 2 matrix of points in \code{[0,1]^2}
#' @param triangle_coords A 3 × 2 matrix of triangle vertex coordinates (X, Y)
#' @return An n × 2 matrix of points mapped into the triangle
#' @export
map_unit_square_to_triangle <- function(qmc_points, triangle_coords) {
  stopifnot(is.matrix(qmc_points), ncol(qmc_points) == 2)
  stopifnot(is.matrix(triangle_coords), dim(triangle_coords) == c(3, 2))

  u <- qmc_points[, 1]
  v <- qmc_points[, 2]

  # Reflect (u, v) across diagonal if outside standard triangle
  reflect <- (u + v > 1)
  u[reflect] <- 1 - u[reflect]
  v[reflect] <- 1 - v[reflect]

  A <- triangle_coords[1, ]
  B <- triangle_coords[2, ]
  C <- triangle_coords[3, ]

  # Apply barycentric transformation
  mapped_x <- (1 - u - v) * A[1] + u * B[1] + v * C[1]
  mapped_y <- (1 - u - v) * A[2] + u * B[2] + v * C[2]

  cbind(mapped_x, mapped_y)
}

#' Map unit square points to a quadrilateral using bilinear interpolation
#'
#' Maps an n × 2 matrix of points from the unit square \code{[0,1]^2} into a quadrilateral
#' defined by 4 coordinates, ordered counterclockwise: bottom-left, bottom-right,
#' top-right, top-left.
#'
#' @param qmc_points An n × 2 matrix of points in \code{[0,1]^2}
#' @param quad_coords A 4 × 2 matrix of quad vertex coordinates (X, Y)
#' @return An n × 2 matrix of points mapped into the quadrilateral
#' @export
map_unit_square_to_quad <- function(qmc_points, quad_coords) {
  stopifnot(is.matrix(qmc_points), ncol(qmc_points) == 2)
  stopifnot(is.matrix(quad_coords), dim(quad_coords) == c(4, 2))

  u <- qmc_points[, 1]
  v <- qmc_points[, 2]

  # Vertices: P1 (bottom-left), P2 (bottom-right), P3 (top-right), P4 (top-left)
  P1 <- quad_coords[1, ]
  P2 <- quad_coords[2, ]
  P3 <- quad_coords[3, ]
  P4 <- quad_coords[4, ]

  # Bilinear interpolation formula:
  mapped_x <- (1 - u) * (1 - v) * P1[1] +
    u * (1 - v) * P2[1] +
    u * v * P3[1] +
    (1 - u) * v * P4[1]

  mapped_y <- (1 - u) * (1 - v) * P1[2] +
    u * (1 - v) * P2[2] +
    u * v * P3[2] +
    (1 - u) * v * P4[2]

  cbind(mapped_x, mapped_y)
}


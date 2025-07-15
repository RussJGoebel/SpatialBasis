#' Map unit square points to an arbitrary convex quadrilateral
#'
#' Applies a bilinear or affine transformation from the unit square [0,1]^2
#' to a user-supplied quadrilateral with 4 vertices.
#'
#' @param u Matrix of QMC points in the unit square. Must be [n × 2] with values in [0,1].
#' @param quad_coords A [4 × 2] matrix of x/y coordinates for the quadrilateral.
#'        Vertices should be ordered (e.g., clockwise), and define a convex quad.
#'
#' @return A matrix of transformed points, [n × 2], mapped into the quadrilateral.
#' @export
map_unit_square_to_quad <- function(u, quad_coords) {
  stopifnot(is.matrix(u), ncol(u) == 2)
  stopifnot(is.matrix(quad_coords), dim(quad_coords)[1] == 4, dim(quad_coords)[2] == 2)

  # Vertices: A, B, C, D (clockwise or counterclockwise)
  A <- quad_coords[1, ]
  B <- quad_coords[2, ]
  C <- quad_coords[3, ]
  D <- quad_coords[4, ]

  # Parametric coordinates
  u1 <- u[, 1]
  u2 <- u[, 2]

  # Bilinear interpolation formula:
  # P(u,v) = (1 - u)(1 - v)A + u(1 - v)B + uvC + (1 - u)vD
  P <- (1 - u1) * (1 - u2)
  Q <- u1       * (1 - u2)
  R <- u1       * u2
  S <- (1 - u1) * u2

  pts <- cbind(
    P * A[1] + Q * B[1] + R * C[1] + S * D[1],
    P * A[2] + Q * B[2] + R * C[2] + S * D[2]
  )

  return(pts)
}

#' Map [0,1]^2 points to an arbitrary triangle using affine transformation
#'
#' This uses the "triangle reflection" method to map QMC points in [0,1]^2
#' to a triangle with vertices (a, b, c).
#'
#' @param u A matrix of [n × 2] QMC points in [0,1]^2
#' @param tri A [3 × 2] matrix of triangle vertices (a, b, c)
#'
#' @return A matrix [n × 2] of transformed coordinates in the triangle
#' @export
map_unit_square_to_triangle <- function(u, tri) {
  stopifnot(ncol(u) == 2, nrow(tri) == 3, ncol(tri) == 2)

  # Reflect QMC points to fit inside the unit triangle
  u1 <- u[, 1]
  u2 <- u[, 2]
  reflect <- (u1 + u2) > 1
  u1[reflect] <- 1 - u1[reflect]
  u2[reflect] <- 1 - u2[reflect]

  a <- tri[1, ]
  b <- tri[2, ]
  c <- tri[3, ]

  # Affine map: T(u) = a + u1*(b - a) + u2*(c - a)
  pts <- matrix(NA_real_, nrow = nrow(u), ncol = 2)
  pts <- sweep(outer(u1, b - a, "*"), 1, a, "+") + outer(u2, c - a, "*")
  dim(pts) <- c(nrow(u), 2)
  pts
}

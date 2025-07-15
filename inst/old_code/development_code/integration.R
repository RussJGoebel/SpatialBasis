#' Fast integration for pixel-style basis functions
#'
#' Computes the matrix of integration weights between basis geometries and target geometries
#' using optimized spatial operations and sparse matrix construction.
#'
#' @param basis A `basis` object with sf POLYGON metadata in `metadata$geometry`
#' @param targets An sf POLYGON object representing regions to integrate over
#' @param mc.cores Number of cores to use
#'
#' @return A sparse matrix [n_targets × n_basis] of integration weights
#' @export
integrate_pixel_basis <- function(basis,
                                  targets,
                                  mc.cores = parallel::detectCores(logical = FALSE)) {
  stopifnot(inherits(basis, "basis"))
  stopifnot(!is.null(basis$metadata$geometry))

  basis_geom <- basis$metadata$geometry
  intersection <- st_parallel_intersection(targets, basis_geom, mc.cores)

  target_areas <- as.vector(sf::st_area(targets))
  intersect_areas <- sf::st_area(intersection)

  intersects <- sf::st_intersects(targets, basis_geom, sparse = FALSE)

  # Encode overlap areas into binary matrix and normalize
  intersects[c(intersects)] <- intersect_areas
  A <- sweep(intersects, 1, target_areas, FUN = "/")

  return(Matrix::Matrix(A, sparse = TRUE))
}

#' Extract coordinate matrices from an sf POLYGON/MULTIPOLYGON object
#'
#' Assumes each geometry has a single outer ring (L1 = 1) with no holes.
#' Used to retrieve coordinates for quadrilateral or triangle-based integration.
#'
#' @param sf_polygons An sf object with POLYGON or MULTIPOLYGON geometries.
#' @param expected_vertices Optional: expected number of unique vertices per polygon.
#'
#' @return A list of matrices, each with dimensions [k × 2], where k is the number of unique vertices.
#' @export
extract_polygon_coords_list <- function(sf_polygons, expected_vertices = NULL) {
  stopifnot(inherits(sf_polygons, "sf"))

  geom_list <- sf::st_geometry(sf_polygons)
  coords_list <- vector("list", length(geom_list))

  for (i in seq_along(geom_list)) {
    coords <- sf::st_coordinates(geom_list[i])
    if (!all(c("X", "Y") %in% colnames(coords))) {
      stop(sprintf("Geometry %d is malformed or missing X/Y coordinates.", i))
    }

    # L1 == 1 always assumed, so extract and deduplicate (first ring only)
    unique_coords <- coords[!duplicated(coords[, c("X", "Y")]), c("X", "Y"), drop = FALSE]

    if (!is.null(expected_vertices) && nrow(unique_coords) != expected_vertices) {
      warning(sprintf("Polygon %d has %d unique vertices (expected %d).",
                      i, nrow(unique_coords), expected_vertices))
    }

    coords_list[[i]] <- as.matrix(unique_coords)
  }

  coords_list
}

#' Integrate a function over many triangles, aggregating by polygon ID
#'
#' @param f Function taking [n × 2] coordinates and returning [n × k] values
#' @param triangle_list A list of triangle coordinate matrices (each [3 × 2])
#' @param qmc QMC points in [0,1]^2
#' @param id_map Optional vector mapping triangles to parent polygons
#' @param average Logical: if TRUE, return area-weighted averages; otherwise integrals
#'
#' @return A matrix [n_polygon × k] of integrals or averages
#' @export
integrate_triangles_qmc <- function(f, triangle_list, qmc, id_map = NULL, average = TRUE) {
  k <- ncol(f(qmc))
  n_tri <- length(triangle_list)

  if (is.null(id_map)) {
    id_map <- rep(1L, n_tri)
  }
  n_polys <- max(id_map)

  out <- matrix(0, nrow = n_polys, ncol = k)
  area_total <- numeric(n_polys)

  for (i in seq_len(n_tri)) {
    tri <- triangle_list[[i]]
    pts <- map_unit_square_to_triangle(qmc, tri)
    fx <- f(pts)

    area <- abs(det(cbind(tri[2, ] - tri[1, ], tri[3, ] - tri[1, ]))) / 2
    id <- id_map[i]

    out[id, ] <- out[id, ] + colMeans(fx) * area
    area_total[id] <- area_total[id] + area
  }

  if (average) {
    out <- sweep(out, 1, area_total, FUN = "/")
  }

  out
}

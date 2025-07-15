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

#' Integrate a function over a single polygon using QMC points
#'
#' Approximates the integral of a function over a polygon via quasi-Monte Carlo (QMC)
#' integration. The polygon is supplied as a matrix of coordinates, and the QMC
#' points are assumed to be defined on the unit square [0,1]^2. The function
#' maps these unit square points to the polygon via bilinear transformation.
#'
#' @param f A function taking an [n × 2] matrix of coordinates and returning a length-n vector of values.
#' @param polygon_coords A matrix of coordinates defining the polygon: must be [4 × 2] for now (convex quad).
#' @param qmc_points A [n × 2] matrix of points in the unit square [0,1]^2 (e.g., from `randtoolbox::sobol()`).
#' @param average A boolean denoting whether to return the average value of f over the region instead of the integral
#' @param polygon_area Optional polygon area (required if `average = FALSE`).
#'
#' @return Approximate integral of `f` over the polygon (a scalar).
#' @export
integrate_single_polygon <- function(f,
                                     qmc_points,
                                     polygon_coords,
                                     polygon_area = NULL,
                                     average = TRUE,
                                     checks = TRUE) {

  if(checks){
    if (!is.matrix(polygon_coords) || nrow(polygon_coords) != 4 || ncol(polygon_coords) != 2) {
      stop("`polygon_coords` must be a 4 × 2 matrix representing a quadrilateral.")
    }
    if (!is.matrix(qmc_points) || ncol(qmc_points) != 2) {
      stop("`qmc_points` must be a matrix with 2 columns (points in [0,1]^2).")
    }
    if (!average && is.null(polygon_area)){
      stop("Integral requires polygon area unless average = TRUE.")
    }
  }



  # Map unit square to polygon
  transformed_pts <- map_unit_square_to_quad(qmc_points, polygon_coords)

  # Evaluate function
  fx <- f(transformed_pts)


  # Return depends on whether only average is needed
  return(if (average) Matrix::colMeans(fx) else Matrix::colMeans(fx) * polygon_area)

}

### //

#' Integrate over a set of triangles with optional aggregation
#'
#' Aggregates triangle integrals or averages back to parent polygons.
#'
#' @param f Function accepting [n x 2] matrix of coordinates.
#' @param triangle_coords_list List of [3 x 2] triangle coordinate matrices.
#' @param qmc_points Matrix of QMC points in unit square [0,1]^2.
#' @param id_map Optional vector mapping triangles to parent polygons.
#' @param average Logical; return average (TRUE) or integral (FALSE).
#'
#' @return If `id_map` is NULL: matrix [1 × k]. Else: matrix [n_parent × k].
#' @export
integrate_triangles_to_polygon <- function(f,
                                        triangle_coords_list,
                                        qmc_points,
                                        id_map = NULL,
                                        average = TRUE) {
  k <- length(f(qmc_points)[1, ])  # infer output size
  n_tri <- length(triangle_coords_list)

  if (is.null(id_map)) {
    id_map <- rep(1L, n_tri)
    n_parent <- 1L
  } else {
    n_parent <- max(id_map)
  }

  result_matrix <- matrix(0, nrow = n_parent, ncol = k)
  area_total <- numeric(n_parent)

  for (i in seq_len(n_tri)) {
    coords <- triangle_coords_list[[i]]
    area <- abs(shoelace_area(coords))
    pts <- map_unit_square_to_triangle(qmc_points, coords)
    fx <- f(pts)
    mean_fx <- colMeans(fx)

    result_matrix[id_map[i], ] <- result_matrix[id_map[i], ] + area * mean_fx
    area_total[id_map[i]] <- area_total[id_map[i]] + area
  }

  if (average) {
    result_matrix <- sweep(result_matrix, 1, area_total, FUN = "/")
  }

  return(result_matrix)
}

#' Compute area of a polygon using the shoelace formula
#'
#' Assumes vertices are ordered (clockwise or counterclockwise).
#'
#' @param coords A matrix with 2 columns and at least 3 rows (x, y coordinates).
#'
#' @return A scalar giving the signed area of the polygon.
#' @export
shoelace_area <- function(coords) {
  x <- coords[, 1]
  y <- coords[, 2]
  n <- nrow(coords)
  sum1 <- sum(x * y[c(2:n, 1)])
  sum2 <- sum(y * x[c(2:n, 1)])
  0.5 * (sum1 - sum2)
}

#' Integrate a function over a triangle using QMC
#'
#' Places QMC points inside a triangle and evaluates a function at those points.
#' Used to compute area-weighted or average integrals over arbitrary polygons.
#'
#' @param f Function accepting [n x 2] matrix of coordinates.
#' @param triangle_coords A [3 x 2] matrix defining triangle vertices.
#' @param qmc_points A [n x 2] matrix of QMC points in unit square [0,1]^2.
#' @param triangle_area Optional precomputed area (otherwise computed internally).
#' @param average Logical; return average (TRUE) or integral (FALSE).
#'
#' @return A length-k numeric vector (or scalar) giving the estimated integral or average.
#' @export
integrate_single_triangle <- function(f,
                                      triangle_coords,
                                      qmc_points = 256,
                                      triangle_area = NULL,
                                      average = TRUE) {
  if (is.null(triangle_area)) {
    triangle_area <- abs(shoelace_area(triangle_coords))
  }

  mapped_pts <- map_unit_square_to_triangle(qmc_points, triangle_coords)
  fx <- f(mapped_pts)
  mean_fx <- colMeans(as.matrix(fx))

  if (average) return(mean_fx) else return(mean_fx * triangle_area)
}


#' Integrate a function over a polygon or list of triangles
#'
#' Supports either triangle-based QMC integration or quadrilateral mapping.
#'
#' @param f Function accepting [n x 2] matrix of coordinates.
#' @param polygon_coords Either a [4 x 2] matrix (for quad), or a list of [3 x 2] triangle matrices.
#' @param qmc_points Matrix of QMC points in unit square [0,1]^2.
#' @param polygon_area Optional area for the polygon.
#' @param average Logical; whether to return average or integral.
#' @param method Either \"triangle\" or \"quad\".
#' @param checks Logical; whether to perform sanity checks.
#' @param id_map Optional ID vector if `polygon_coords` is a list of triangles.
#'
#' @return Vector of length k (basis dimension)
#' @export
integrate_single_polygon_general <- function(f,
                                             polygon_coords,
                                             qmc_points,
                                             polygon_area = NULL,
                                             average = TRUE,
                                             method = c("triangle", "quad"),
                                             checks = TRUE,
                                             id_map = NULL) {
  method <- match.arg(method)

  if (checks) {
    if (!is.matrix(qmc_points) || ncol(qmc_points) != 2) {
      stop("`qmc_points` must be a [n x 2] matrix in the unit square.")
    }
    if (!average && is.null(polygon_area) && method == "quad") {
      stop("Integral requires polygon_area unless average = TRUE.")
    }
  }

  if (method == "quad") {
    if (!is.matrix(polygon_coords) || nrow(polygon_coords) != 4) {
      stop("Quadrilateral method requires a 4 × 2 matrix.")
    }
    transformed_pts <- map_unit_square_to_quad(qmc_points, polygon_coords)
    fx <- f(transformed_pts)
    return(if (average) colMeans(fx) else colMeans(fx) * polygon_area)
  } else {
    # polygon_coords should now be a list of triangle matrices
    return(integrate_triangles_to_polygon(
      f = f,
      triangle_coords_list = polygon_coords,
      qmc_points = qmc_points,
      id_map = id_map,
      average = average
    ))
  }
}


### //


#' Integrate a basis over polygons using QMC (single core)
#'
#' Approximates the integral or average of each basis function over each polygon
#' using quasi-Monte Carlo (QMC) integration. The basis must be a `basis` object.
#'
#' @param basis A basis object created via `make_basis()`, with `evaluate_fn`.
#' @param sf_polygons An sf object with convex quadrilateral POLYGONs or MULTIPOLYGONs.
#' @param n_qmc Number of QMC points per polygon.
#' @param average If TRUE, return average value over each polygon; otherwise return integral.
#' @param checks Logical; whether to perform input checks in each polygon integration.
#'
#' @return A [n_polygons × k] numeric matrix of integrated basis evaluations.
#' @export
integrate_sf_polygons_qmc_single_core <- function(basis,
                                                  sf_polygons,
                                                  n_qmc = 256,
                                                  average = TRUE,
                                                  checks = TRUE) {
  stopifnot(inherits(basis, "basis"))
  f <- basis$evaluate_fn

  # Generate QMC points once
  qmc <- generate_qmc_unit_square(n_qmc)

  # Get polygon coordinate matrices
  coords_list <- extract_polygon_coords_list(sf_polygons, expected_vertices = 4)

  # Compute areas if needed
  if (!average) {
    areas <- as.numeric(sf::st_area(sf_polygons))
  } else {
    areas <- rep(NA_real_, length(coords_list))  # dummy placeholder
  }

  results <- do.call(rbind, lapply(seq_along(coords_list), function(i) {
    integrate_single_polygon(
      f = f,
      qmc_points = qmc,
      polygon_coords = coords_list[[i]],
      polygon_area = areas[i],
      average = average,
      checks = checks
    )
  }))

  return(results)
}



#' Integrate a basis over all polygons using QMC (parallelized)
#'
#' Approximates the integral or average of each basis function over each polygon
#' using quasi-Monte Carlo (QMC) integration. Assumes each polygon is a convex quadrilateral.
#' Uses `future` and `future.apply` for parallelization.
#'
#' @param basis A basis object created via `make_basis()`, with an `evaluate_fn`.
#' @param sf_polygons An sf object with convex quadrilateral POLYGON or MULTIPOLYGON geometries.
#' @param n_qmc Number of QMC points to use per polygon.
#' @param average Logical; if TRUE, return average value over each polygon instead of integral.
#' @param checks Logical; whether to perform input checks inside each polygon integration call.
#' @param workers Integer; number of parallel workers to use (default = available cores - 1).
#'
#' @return A numeric matrix of dimension [n_polygons × k], where k is the number of basis functions.
#' @export
integrate_sf_polygons_qmc <- function(basis,
                                      sf_polygons,
                                      n_qmc = 256,
                                      average = TRUE,
                                      checks = TRUE,
                                      workers = future::availableCores() - 1,
                                      method = c("triangle","quad")) {
  stopifnot(inherits(basis, "basis"))
  method <- match.arg(method)


  sf_polygons <- ensure_projected(sf_polygons)

  if(method == "triangle"){

    message("Converting polygons into triangles...")

    sf_polygons <-  sf::st_triangulate(sf_polygons)


  message("There are ", dim(sf_polygons)[1], " triangular regions.")


  }



  f <- basis$evaluate_fn

  # Generate QMC points once
  qmc <- generate_qmc_unit_square(n_qmc)

  # Get polygon coordinate matrices

  if(method == "quad"){expected_vertices = 4}
  if(method == "triangle"){expected_vertices = 3}

  print("HI")

  coords_list <- extract_polygon_coords_list(sf_polygons, expected_vertices = expected_vertices)

  print("OK")

  # Compute areas if needed
  if (!average) {
    areas <- as.numeric(sf::st_area(sf_polygons))
  } else {
    areas <- rep(NA_real_, length(coords_list))  # dummy placeholder
  }

  # Set up temporary parallel plan if current plan is sequential
old_plan <- future::plan()
if (inherits(old_plan, "sequential")) {
  future::plan(future::multisession, workers = workers)
  on.exit(future::plan(old_plan), add = TRUE)
}

  # Parallel QMC integration
  results <- future.apply::future_lapply(seq_along(coords_list), function(i) {
    integrate_single_polygon_general(
      f = f,
      qmc_points = qmc,
      polygon_coords = coords_list[[i]],
      polygon_area = areas[i],
      average = average,
      checks = checks,
      method = method
    )
  }, future.seed = NULL)

  return(do.call(rbind, results))
}

#' Integrate a basis over all polygons using stratified QMC (parallelized)
#'
#' Approximates the integral or average of each basis function over each polygon
#' using quasi-Monte Carlo (QMC) integration stratified by auxiliary subregions.
#' Assumes each polygon is a convex quadrilateral. Uses `future` and `future.apply` for parallelization.
#'
#' @param basis A basis object created via `make_basis()`.
#' @param sf_targets An sf object with convex quadrilateral POLYGON or MULTIPOLYGON geometries.
#' @param stratify_by An sf object with subregions to stratify sampling (e.g., landcover).
#' @param n_qmc Number of QMC points to use per subregion.
#' @param average Logical; if TRUE, return average value over each polygon instead of integral.
#' @param checks Logical; whether to perform input checks inside each polygon integration call.
#' @param workers Integer; number of parallel workers to use (default = available cores - 1).
#'
#' @return A numeric matrix of dimension [n_targets × k], where k is the number of basis functions.
#' @export
integrate_sf_polygons_qmc_stratify <- function(basis,
                                      sf_targets,
                                      stratify_by,
                                      n_qmc = 256,
                                      average = TRUE,
                                      checks = TRUE,
                                      workers = future::availableCores() - 1) {
  stopifnot(inherits(basis, "basis"))
  f <- basis$evaluate_fn

  message("Intersecting targets with strata...")
  # Intersect targets with strata
  intersected <- st_parallel_intersection(sf_targets, stratify_by)
  intersected$sub_id <- seq_len(nrow(intersected))

  message(dim(intersected)[1]," polygons in intersection.")

  # Record mapping from subregions back to targets
  id_map <- sf::st_intersects(intersected, sf_targets)
  target_ids <- vapply(id_map, function(x) if (length(x)) x[1] else NA_integer_, integer(1))

  # Generate QMC points once
  qmc <- generate_qmc_unit_square(n_qmc)

  # Get polygon coordinate matrices
  coords_list <- extract_polygon_coords_list(intersected, expected_vertices = 4)
  areas <- as.numeric(sf::st_area(intersected))

  message("Applying qmc integration...")

  # Set up temporary parallel plan
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = workers)

  # Parallel QMC integration on subregions
  sub_results <- future.apply::future_lapply(seq_along(coords_list), function(i) {
    downscaling::integrate_single_polygon(
      f = f,
      qmc_points = qmc,
      polygon_coords = coords_list[[i]],
      polygon_area = areas[i],
      average = FALSE,
      checks = checks
    )
  }, future.seed = NULL)

  # Re-aggregate to original targets
  k <- basis$k
  n_targets <- nrow(sf_targets)
  result_matrix <- matrix(0, nrow = n_targets, ncol = k)
  total_area <- numeric(n_targets)

  for (i in seq_along(sub_results)) {
    tid <- target_ids[i]
    if (is.na(tid)) next
    result_matrix[tid, ] <- result_matrix[tid, ] + sub_results[[i]]
    total_area[tid] <- total_area[tid] + areas[i]
  }

  if (average) {
    result_matrix <- sweep(result_matrix, 1, total_area, FUN = "/")
  }

  return(result_matrix)
}


#' Integrate a basis over polygons, with optional triangle subdivision
#'
#' This function integrates a basis function over an sf object, with optional
#' triangulation of input polygons for higher accuracy in irregular shapes.
#'
#' @param basis A basis object with an `evaluate_fn()` function.
#' @param sf_polygons An sf object of polygons to integrate over.
#' @param n_qmc Number of QMC points per shape (or per triangle if triangulated).
#' @param average Logical; if TRUE, return average rather than integral.
#' @param checks Logical; perform polygon checks.
#' @param triangles Logical; if TRUE, triangulate polygons first.
#' @param workers Number of parallel workers to use.
#'
#' @return A matrix [n_polygons x k] of integrated basis evaluations.
#' @export
integrate_sf_polygons_qmc_general <- function(basis,
                                              sf_polygons,
                                              n_qmc = 256,
                                              average = TRUE,
                                              checks = TRUE,
                                              triangles = FALSE,
                                              workers = future::availableCores() - 1) {
  stopifnot(inherits(basis, "basis"))
  f <- basis$evaluate_fn

  # Optionally triangulate and record mapping back to parent polygons
  if (triangles) {

    messsage("Converting sf_polygon to triangles...")

    triangle_geom <- sf::st_triangulate(sf_polygons) |>
      sf::st_collection_extract("POLYGON")

    id_map <- rep(seq_len(nrow(sf_polygons)), lengths(sf::st_geometry(sf::st_triangulate(sf_polygons))))
    sf_polygons <- sf::st_sf(geometry = triangle_geom)

    message("There are ", dim(sf_polygons)[1]," triangular geometries.")

  }

  # Generate QMC points once
  qmc <- generate_qmc_unit_square(n_qmc)

  # Extract coordinates
  coords_list <- extract_polygon_coords_list(sf_polygons, expected_vertices = if (triangles) 3 else 4)

  # Compute areas
  areas <- if (!average) as.numeric(sf::st_area(sf_polygons)) else rep(NA_real_, length(coords_list))

  # Set up parallel plan
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = workers)

  # Integrate over all triangles or polygons
  results <- future.apply::future_lapply(seq_along(coords_list), function(i) {
    if (triangles) {
      integrate_single_triangle(
        f = f,
        triangle_coords = coords_list[[i]],
        qmc_points = qmc,
        triangle_area = areas[i],
        average = average
      )
    } else {
      integrate_single_polygon(
        f = f,
        qmc_points = qmc,
        polygon_coords = coords_list[[i]],
        polygon_area = areas[i],
        average = average,
        checks = checks
      )
    }
  }, future.seed = NULL)

  # Combine and optionally aggregate
  results_mat <- do.call(rbind, results)

  if (triangles) {
    # Aggregate back to original polygon structure
    k <- basis$k
    n_parent <- nrow(basis$geometry)  # or use nrow from original input if stored elsewhere
    final_mat <- matrix(0, nrow = n_parent, ncol = k)

    if (!average) {
      for (i in seq_along(id_map)) {
        final_mat[id_map[i], ] <- final_mat[id_map[i], ] + results_mat[i, ]
      }
    } else {
      triangle_areas <- as.numeric(sf::st_area(sf_polygons))
      total_area <- tapply(triangle_areas, id_map, sum)
      for (i in seq_along(id_map)) {
        final_mat[id_map[i], ] <- final_mat[id_map[i], ] + results_mat[i, ] * triangle_areas[i]
      }
      final_mat <- sweep(final_mat, 1, total_area, FUN = "/")
    }
    return(final_mat)
  }

  return(results_mat)
}


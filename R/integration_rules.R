###############################################################################
###################### High-level evaluation functions ########################
###############################################################################

#' Generic integration dispatcher for region × rule combinations
#' @export
evaluate_integral <- function(rule, region, basis, index = NULL) {
  if (inherits(region, "region_composite")) {
    return(evaluate_integral.region_composite(rule, region, basis, index = index))
  }
  UseMethod("evaluate_integral", rule)
}

#' Generic integration dispatcher for region × rule combinations
#' @export
evaluate_integral.default <- function(rule, region, basis, index = NULL) {
  stop("No integration method for rule class: ", class(rule)[1])
}



#' Streaming Composite Region Evaluator (Memory-Safe)
#'
#' Evaluates a composite region by streaming over its atomic parts,
#' accumulating a weighted average of basis evaluations.
#' Avoids memory blow-up by never materializing full 3D arrays.
#'
#' @param region A region_composite object (with $parts and $areas)
#' @param rule An integration rule (e.g., qmc_average)
#' @param basis A basis object (with $evaluate_fn and $k)
#'
#' @return A numeric vector of length basis$k
#'
#' @export
evaluate_integral.region_composite <- function(rule, region, basis,index = NULL) {
  parts <- region$parts
  weights <- region$areas / sum(region$areas)

  # Preallocate output accumulator
  out <- numeric(basis$k)

  # Stream through parts, accumulating weighted results
  for (i in seq_along(parts)) {
    val_i <- evaluate_integral(rule, parts[[i]], basis, index = index)
    out <- out + weights[i] * val_i
  }

  return(out)
}

###############################################################################
################ Specific Integration Rule / Evaluators #######################
###############################################################################

## ----------------------------------------------------------------------------


#' QMC-Based Integration Rule (Average of Basis Values)
#'
#' Constructs an integration rule that computes the average value of a basis function
#' over a polygonal region using quasi-Monte Carlo (QMC) integration.
#'
#' @param qmc_points An n×2 matrix of QMC points in the unit square.
#' @param sensitivity_params_list An optional list of parameter vectors, one per observation.
#' @param sensitivity_dispatcher An optional function(coords, params) returning weights.
#'
#' @return An integration rule object.
#' @export
integration_rule_qmc_average <- function(
    qmc_points,
    sensitivity_params_list = NULL,
    sensitivity_dispatcher = NULL
) {
  if (!is.null(sensitivity_params_list) && is.null(sensitivity_dispatcher)) {
    stop("If sensitivity_params_list is provided, you must also provide sensitivity_dispatcher.")
  }
  structure(
    list(
      type = "qmc_average",
      qmc_points = qmc_points,
      sensitivity_params = sensitivity_params_list,
      sensitivity_dispatcher = sensitivity_dispatcher
    ),
    class = c("qmc_average", "integration_rule")
  )
}

#' QMC Integration for an Atomic Region (with Optional Sensitivity Weights)
#'
#' Evaluates the (weighted) average value of basis functions over an atomic region
#' using quasi-Monte Carlo (QMC) integration.
#'
#' @param rule An integration rule object of class "qmc_average".
#' @param region A region_atomic object (must have $shape and $coords).
#' @param basis A basis object (with $evaluate_fn).
#' @param index An integer index of the observation (used to pick sensitivity function).
#' @param tolerance A tolerance for sensitivity to set overall weight to 0.
#'
#' @return A numeric vector of length equal to basis$k.
#'
#' @export
evaluate_integral.qmc_average <- function(rule, region, basis, index = NULL, tolerance = 1e-6) {
  stopifnot(inherits(region, "region_atomic"))
  stopifnot(!is.null(region$shape), region$shape %in% c("triangle", "quad"))

  qmc <- rule$qmc_points

  # Map QMC points
  coords_mapped <- switch(
    region$shape,
    "triangle" = map_unit_square_to_triangle(qmc, region$coords),
    "quad"     = map_unit_square_to_quad(qmc, region$coords),
    stop("Unsupported shape in region: must be 'triangle' or 'quad'")
  )

  phi_vals <- basis$evaluate_fn(coords_mapped)

  if (!is.null(rule$sensitivity_params)) {
    if (is.null(index)) stop("index must be provided when sensitivity_params are used.")
    sens_weights <- rule$sensitivity_dispatcher(coords_mapped, index, rule$sensitivity_params)
    stopifnot(length(sens_weights) == nrow(coords_mapped))
    if(sum(sens_weights) < 1e-6) {weighted_means <- rep(0,length(sens_weights))} else {
      weighted_means <- colSums(phi_vals * sens_weights) / sum(sens_weights)}
  } else {
    weighted_means <- colMeans(phi_vals)
  }

  return(weighted_means)
}

## ----------------------------------------------------------------------------


#' Area Overlap Integration Rule
#'
#' Constructs an integration rule that computes the proportion of overlap area
#' between an atomic region and each cell in the basis grid. This is typically used
#' when the basis is a geometry basis with indicator functions over polygons.
#'
#' @return An object of class \code{"area_overlap"} and \code{"integration_rule"}.
#'
#' @examples
#' rule_area <- integration_rule_area_overlap()
#'
#' @export
integration_rule_area_overlap <- function() {
  structure(
    list(
      type = "area_overlap"
    ),
    class = c("area_overlap", "integration_rule")
  )
}

#' Area Overlap Integration for an Atomic Region (Supports Weighted Piecewise Bases)
#'
#' Computes the area-weighted integral of a basis function over a region by computing
#' the proportion of overlap area between the region and each polygon in the basis.
#'
#' - If the basis is a standard indicator geometry basis, it returns the overlap vector.
#' - If the basis is a piecewise constant spatial covariate with known values,
#'   it returns the weighted sum over those values.
#'
#' @param rule An integration rule object of class "area_overlap"
#' @param region A region_atomic object (with a valid geometry)
#' @param basis A basis object (e.g., geometry or piecewise_constant basis)
#' @param index Ignored for area_overlap (included for compatibility)
#'
#' @return A numeric vector:
#'   - of length `basis$k` if using a piecewise-constant basis
#'   - of length equal to number of basis polygons if using indicator geometry basis
#'
#' @export
evaluate_integral.area_overlap <- function(rule, region, basis, index = NULL) {
  stopifnot(inherits(region, "region_atomic"))

  old_s2 <- sf::sf_use_s2()
  suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(old_s2)), add = TRUE)

  # Normalize region geometry
  geom_region <- region$geometry
  if (inherits(geom_region, "sfg")) geom_region <- sf::st_sfc(geom_region)
  if (inherits(geom_region, "sf")) geom_region <- sf::st_geometry(geom_region)
  stopifnot(inherits(geom_region, "sfc"))

  # Normalize basis geometry
  geom_basis <- basis$geometry
  if (inherits(geom_basis, "sfg")) geom_basis <- sf::st_sfc(geom_basis)
  if (inherits(geom_basis, "sf")) geom_basis <- sf::st_geometry(geom_basis)
  stopifnot(inherits(geom_basis, "sfc"))

  # Set CRS
  sf::st_crs(geom_region) <- sf::st_crs(geom_basis)

  # Identify which grid cells intersect the region
  idx_grid <- suppressMessages(
    which(sf::st_intersects(
      geom_basis,
      geom_region,
      sparse = TRUE,
      prepared = TRUE
    ) %>% lengths() > 0)
  )

  # If no intersections, return zeros
  if (length(idx_grid) == 0) {
    if (!is.null(basis$values)) return(rep(0, basis$k))
    return(rep(0, length(geom_basis)))
  }

  # Subset basis to only relevant cells and add ID
  basis_subset <- sf::st_sf(id = idx_grid, geometry = geom_basis[idx_grid])

  # Compute intersections
  intersection <- sf::st_intersection(basis_subset, geom_region)
  intersection <- sf::st_make_valid(intersection)

  if (nrow(intersection) == 0) {
    if (!is.null(basis$values)) return(rep(0, basis$k))
    return(rep(0, length(geom_basis)))
  }

  overlap_areas <- sf::st_area(intersection)

  # Aggregate area per cell index
  agg <- tapply(
    as.numeric(overlap_areas),
    intersection$id,
    sum,
    default = 0
  )

  # Full overlap vector (length = number of polygons)
  overlap <- numeric(length(geom_basis))
  overlap[as.integer(names(agg))] <- agg

  # Normalize so sum = 1
  total_area <- sum(overlap)
  if (total_area > 0) {
    overlap <- overlap / total_area
  }

  # Return either:
  # (a) vector of overlaps (standard geometry basis)
  # (b) weighted sum of values (piecewise constant)
  if (!is.null(basis$values)) {
    # basis$values is (n_polygons × k), overlap is length-n_polygons
    # return a length-k numeric vector
    return(as.vector(drop(overlap %*% basis$values)))
  } else {
    return(overlap)
  }
}

## ----------------------------------------------------------------------------

#' Centroid Gaussian Integration Rule
#'
#' Creates an integration rule that computes Gaussian weights over grid centroids.
#'
#' @param grid_centroids A matrix of grid cell centroids (n_cells × 2).
#' @param params_list A list of parameter vectors per observation (mu_x, mu_y, sig_x, sig_y, angle).
#'
#' @return An integration rule object.
#' @export
integration_rule_centroid_gaussian <- function(grid_centroids, params_list) {
  structure(
    list(
      type = "centroid_gaussian",
      centroids = grid_centroids,
      params_list = params_list
    ),
    class = c("centroid_gaussian", "integration_rule")
  )
}

#' Centroid Gaussian integration for a region (ignored)
#'
#' @export
evaluate_integral.centroid_gaussian <- function(rule, region, basis, index) {

  params <- rule$params_list[[index]]
  mu_x <- params[["x0"]]
  mu_y <- params[["y0"]]
  sig_x <- params[["sig_x"]]
  sig_y <- params[["sig_y"]]
  angle <- params[["angle"]]

  centroids <- sf::st_coordinates(rule$centroids)

  delta_x <- centroids[,1] - mu_x
  delta_y <- centroids[,2] - mu_y
  d <- sqrt(delta_x^2 + delta_y^2)
  ang <- atan2(delta_y, delta_x)
  dx <- d * cos(ang - angle)
  dy <- d * sin(ang - angle)

  w <- exp(-dx^2 / sig_x) * exp(-dy^2 / sig_y)
  w <- w / sum(w)

  return(w)
}






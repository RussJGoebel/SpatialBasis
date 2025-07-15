#' Pixel Basis Constructor (Compatible with Integration Framework)
#'
#' Constructs a pixel-level indicator basis from an sf object of polygons.
#' The returned object has an evaluate() method and stores the geometry metadata
#' required for integration via QMC or pixel overlap.
#'
#' @param pixel_sf An sf object with POLYGON/MULTIPOLYGON geometry. Must contain a unique `pixel_id` column.
#'
#' @return An object of class "basis" created with make_basis().
#' @export
make_pixel_basis <- function(pixel_sf) {
  if (!inherits(pixel_sf, "sf")) stop("Input must be an sf object.")
  if (!all(sf::st_geometry_type(pixel_sf) %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("All geometries must be POLYGON or MULTIPOLYGON.")
  }

  pixel_sf <- ensure_projected(pixel_sf)

  # Ensure pixel_id is present and valid
  pixel_sf$pixel_id <- seq_len(nrow(pixel_sf))

  k <- nrow(pixel_sf)
  crs <- sf::st_crs(pixel_sf)

  evaluate_fn <- function(locations) {
    if (inherits(locations, "sf")) {
      point_sf <- locations
    } else if (is.matrix(locations) || is.data.frame(locations)) {
      df <- as.data.frame(locations)
      names(df) <- c("x", "y")
      point_sf <- sf::st_as_sf(df, coords = c("x", "y"), crs = NA)
    } else {
      stop("Unsupported input type: must be sf, matrix, or data.frame.")
    }

    # Ensure CRS is assigned and matches
    if (is.na(sf::st_crs(point_sf))) {
      sf::st_crs(point_sf) <- crs
    } else if (sf::st_crs(point_sf) != crs) {
      message("Transforming evaluation points to match pixel basis CRS.")
      point_sf <- sf::st_transform(point_sf, crs)
    }

    joined <- sf::st_join(point_sf, pixel_sf[, "pixel_id"], join = sf::st_within, left = TRUE)
    match_idx <- joined$pixel_id

    out <- matrix(0, nrow = nrow(point_sf), ncol = k)
    valid <- !is.na(match_idx)
    out[cbind(which(valid), match_idx[valid])] <- 1
    return(out)
  }

  make_basis(
    evaluate_fn = evaluate_fn,
    k = k,
    crs = crs,
    metadata = list(
      type = "pixel",
      geometry = pixel_sf
    )
  )
}


#' Create a Matérn Basis Object
#'
#' Constructs a spatial basis using eigenfunctions of the Laplacian with Matérn-type decay.
#' Evaluates the top `k` eigenfunctions on a rectangular domain using separable cosine basis.
#'
#' @param k Number of basis functions.
#' @param L Square domain side length
#' @param nu Smoothness parameter of the Matérn (must be > 0).
#' @param rho Range parameter of the Matérn.
#' @param kappa Optional direct specification of the SPDE scale parameter (overrides rho if given).
#'
#' @return A basis object of class "basis" with an evaluate() method.
#' @export
make_matern_basis <- function(k, L = 1, nu = 1, rho = NULL, kappa = NULL) {
  stopifnot(k > 0, nu > 0)

  Lx <- Ly <- L

  # Determine kappa from rho if not supplied directly
  if (is.null(kappa)) {
    if (is.null(rho)) stop("Either `rho` or `kappa` must be provided.")
    stopifnot(rho > 0)
    kappa <- sqrt(8 * nu) / rho
  } else {
    stopifnot(kappa > 0)
  }

  # Get top k (m,n) pairs based on m^2 + n^2
  mn_pairs <- top_matern_indices(k)

  # Compute decay exponent
  alpha <- nu + 1

  # Eigenvalue decay function
  eval_lambda <- function(m, n) {
    val <- kappa^2 + (pi * m / Lx)^2 + (pi * n / Ly)^2
    val^(-alpha)
  }

  weights <- mapply(
    eval_lambda,
    m = mn_pairs$m,
    n = mn_pairs$n
  )

  evaluate_fn <- function(coords) {
    coords <- as.matrix(coords)
    x <- coords[, 1]
    y <- coords[, 2]

    basis_matrix <- matrix(0, nrow = length(x), ncol = k)
    for (j in seq_len(k)) {
      m <- mn_pairs$m[j]
      n <- mn_pairs$n[j]
      basis_matrix[, j] <- sqrt(2) * cos(pi * m * x / Lx) * sqrt(2) * cos(pi * n * y / Ly)
    }
    sweep(basis_matrix, 2, sqrt(weights), `*`)
  }

  metadata <- list(type = "matern", nu = nu, rho = rho, kappa = kappa, domain = c(Lx, Ly))
  make_basis(evaluate_fn = evaluate_fn, k = k, metadata = metadata)
}


#' Constant basis function for integration testing
#'
#' Always returns a matrix of 1s with one column.
#'
#' @return A `basis` object with a constant function.
#' @export
make_constant_basis <- function() {
  evaluate_fn <- function(coords) {
    stopifnot(is.matrix(coords), ncol(coords) == 2)
    matrix(1, nrow = nrow(coords), ncol = 1)
  }

  make_basis(evaluate_fn = evaluate_fn, k = 1, metadata = list(name = "constant"))
}

#' Efficiently compute the top-k (m, n) indices for Laplacian eigenvalues
#'
#' Uses a frontier-expanding method to avoid full grid enumeration.
#' Returns smallest k values of m^2 + n^2.
#'
#' @param k Number of indices to return
#' @return A data.frame with columns m, n, lambda_sq
#' @export
top_matern_indices <- function(k) {
  visited <- new.env(hash = TRUE, parent = emptyenv())
  encode <- function(m, n) paste0(m, ",", n)

  frontier <- list(c(1, 1))
  visited[[encode(1,1)]] <- TRUE

  result <- data.frame(m = integer(0), n = integer(0), lambda_sq = numeric(0))

  while (nrow(result) < k) {
    # Find min lambda_sq among frontier
    lambda_vals <- sapply(frontier, function(p) sum(p^2))
    i_min <- which.min(lambda_vals)
    p_min <- frontier[[i_min]]
    frontier[[i_min]] <- NULL

    m <- p_min[1]
    n <- p_min[2]
    result <- rbind(result, data.frame(m = m, n = n, lambda_sq = m^2 + n^2))

    # Add neighbors if not visited
    for (neighbor in list(c(m + 1, n), c(m, n + 1))) {
      key <- encode(neighbor[1], neighbor[2])
      if (is.null(visited[[key]])) {
        frontier[[length(frontier) + 1]] <- neighbor
        visited[[key]] <- TRUE
      }
    }
  }

  return(result)
}




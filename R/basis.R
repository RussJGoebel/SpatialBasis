#' Create a Spatial Basis Object
#'
#' Constructs a basis object used to evaluate spatial basis functions for downscaling.
#' A basis is defined by a function that maps coordinates to a matrix of basis evaluations,
#' the number of basis functions used (`k`), an optional coordinate reference system (`crs`),
#' and additional optional metadata.
#'
#' @param evaluate_fn A function of the form `f(coords)` that returns a matrix of basis
#'   evaluations. `coords` should be a matrix or data frame with two columns (x and y).
#'   The function must return a matrix of dimension `[nrow(coords), k]`.
#' @param k An integer specifying the number of basis functions (i.e., the number of columns
#'   in the matrix returned by `evaluate_fn`).
#' @param crs Optional. A coordinate reference system, either as an integer EPSG code,
#'   a character string (e.g., `"EPSG:4326"`), or an `sf::st_crs()` object. This defines
#'   the coordinate system in which the basis is constructed and evaluated. If `NULL`,
#'   the basis is assumed to be unitless (e.g., defined on \[0,1\]^2).
#' @param metadata An optional list of additional metadata, such as kernel parameters,
#'   construction logs, or user-defined tags. Ignored by default methods.
#'
#' @return A list of class `"basis"` with fields:
#' \describe{
#'   \item{evaluate_fn}{The user-defined function for evaluating basis functions at input coordinates.}
#'   \item{k}{The number of basis functions.}
#'   \item{crs}{The coordinate reference system (or `NULL` for unitless bases).}
#'   \item{metadata}{A list of optional metadata elements.}
#' }
#'
#' @examples
#' # Example: basis defined on the unit square with two basis functions
#' my_eval <- function(coords) {
#'   cbind(sin(pi * coords[, 1]), cos(pi * coords[, 2]))
#' }
#' b <- make_basis(my_eval, k = 2, crs = NA)
#' b$evaluate_fn(matrix(c(0.5, 0.5), ncol = 2))
#'
#' @export
make_basis <- function(evaluate_fn, k, crs = NULL, metadata = list()) {
  if (!is.function(evaluate_fn)) {
    stop("`evaluate_fn` must be a function that returns a matrix of basis evaluations.")
  }
  if (!is.numeric(k) || length(k) != 1 || k < 1 || k != as.integer(k)) {
    stop("`k` must be a positive integer.")
  }
  if (!is.null(crs)) {
    crs <- sf::st_crs(crs)  # Validates or converts to st_crs object
    if (is.na(crs)) {
      warning("Provided CRS could not be interpreted; using NULL.")
      crs <- NULL
    }
  }

  structure(
    list(
      evaluate_fn = evaluate_fn,
      k = k,
      crs = crs,
      metadata = metadata
    ),
    class = "basis"
  )
}



#' Construct a Geometry-Based Spatial Basis
#'
#' Creates a spatial basis where each function is an indicator over an `sf` polygon.
#' Evaluates to 1 inside the polygon, 0 otherwise. Integration is typically area-based.
#'
#' @param geometry An `sf` object of polygons (one per basis function).
#' @param crs Optional: coordinate reference system. Inferred from geometry if omitted.
#' @param metadata Optional: auxiliary metadata (e.g., resolution, tags).
#'
#' @return An object of class `geometry_basis`, `spatial_basis`, and `basis`.
#' @export
make_geometry_basis <- function(geometry, crs = NULL, metadata = list()) {
  stopifnot(inherits(geometry, "sf"))
  if (is.null(crs)) crs <- sf::st_crs(geometry)

  evaluate_fn <- function(coords) {
    pts <- sf::st_as_sf(data.frame(x = coords[, 1], y = coords[, 2]),
                        coords = c("x", "y"), crs = crs)

    idx_list <- sf::st_intersects(pts, geometry, sparse = TRUE, prepared = TRUE)
    idx_list <- lapply(idx_list, function(x) if (length(x) > 0) x[1] else integer(0)) # there's some goofy boundary behavior and this is really just for testing qmc
    n <- length(idx_list)
    k <- nrow(geometry)

    # Build sparse matrix of indicators
    i <- rep(seq_len(n), lengths(idx_list))
    j <- unlist(idx_list)
    out_sparse <- Matrix::sparseMatrix(i = i, j = j, x = 1, dims = c(n, k))

    as.matrix(out_sparse)
  }

  structure(
    list(
      evaluate_fn = evaluate_fn,
      k = nrow(geometry),
      crs = sf::st_crs(crs),
      geometry = geometry,
      type = "indicator",
      metadata = metadata
    ),
    class = c("geometry_basis", "basis")
  )
}

#' Construct a Piecewise Constant Basis from Polygonal Geometry and Values
#'
#' Creates a spatial basis object where a spatial function is defined piecewise constant
#' over a set of polygons (e.g., grid cells). This is typically used for representing
#' covariates or low-rank models with known spatial structure.
#'
#' The resulting basis has \code{k = 1} by default (scalar field). You can supply
#' a matrix of values to construct a multivariate basis with \code{k > 1}.
#'
#' @param geometry An \code{sf} object of polygons defining the support of the basis.
#'                 Each row defines one spatial tile (e.g., grid cell).
#' @param values A numeric vector of length \code{nrow(geometry)} (for scalar field),
#'               or a numeric matrix with \code{nrow(values) == nrow(geometry)}
#'               and \code{k} columns for multivariate fields.
#' @param crs Optional CRS; inferred from geometry if missing.
#' @param metadata Optional metadata list.
#'
#' @return A \code{"basis"} object (subclass \code{"piecewise_constant_basis"}) with:
#'   \itemize{
#'     \item \code{geometry}: the input sf polygons
#'     \item \code{values}: the piecewise constant values per polygon (matrix with k columns)
#'     \item \code{evaluate_fn}: returns values matched to points (slow; uses st_intersects)
#'     \item \code{k}: number of basis functions (columns of values)
#'     \item \code{crs}: spatial reference
#'     \item \code{metadata}: user-defined metadata
#'   }
#'
#' @export
make_piecewise_constant_basis <- function(geometry, values, crs = NULL, metadata = list()) {
  stopifnot(inherits(geometry, "sf"))

  n <- nrow(geometry)

  # Normalize values to matrix
  if (is.vector(values)) {
    values <- matrix(values, ncol = 1)
  }
  if (!is.matrix(values) || nrow(values) != n) {
    stop("`values` must be a vector of length nrow(geometry) or a matrix with that many rows.")
  }

  k <- ncol(values)

  # Determine CRS
  if (is.null(crs)) crs <- sf::st_crs(geometry)

  # Evaluation function (slow fallback)
  evaluate_fn <- function(coords) {
    pts <- sf::st_as_sf(data.frame(x = coords[, 1], y = coords[, 2]),
                        coords = c("x", "y"), crs = crs)
    idx_list <- sf::st_intersects(pts, geometry, sparse = TRUE, prepared = TRUE)
    idx_vec <- sapply(idx_list, function(x) if (length(x) > 0) x[1] else NA_integer_)

    out <- matrix(NA_real_, nrow = nrow(coords), ncol = k)
    valid <- !is.na(idx_vec)
    out[valid, ] <- values[idx_vec[valid], , drop = FALSE]
    out
  }

  structure(
    list(
      geometry = geometry,
      values = values,
      evaluate_fn = evaluate_fn,
      k = k,
      crs = crs,
      type = "piecewise_constant",
      metadata = metadata
    ),
    class = c("piecewise_constant_basis", "basis")
  )
}


#' Construct a Function-Based Spatial Basis
#'
#' Creates a spatial basis where each function is evaluated numerically at coordinates
#' (e.g., KL modes, splines, Fourier basis). Integration should be done using QMC or quadrature.
#'
#' @param evaluate_fn A function that takes an n × 2 coordinate matrix and returns an n × k matrix.
#' @param k Optional. Number of basis functions. Inferred if not provided.
#' @param crs A coordinate reference system (passed to `sf::st_crs()`).
#' @param type A label for the basis type (e.g., "kl", "fourier", "spline").
#' @param metadata Optional metadata list.
#'
#' @return An object of class `function_basis`, `spatial_basis`, and `basis`.
#' @export
make_function_basis <- function(evaluate_fn, k = NULL, crs = NULL,
                                type = "custom", metadata = list()) {
  stopifnot(is.function(evaluate_fn))

  if (is.null(k)) {
    test_coords <- matrix(c(0, 0), ncol = 2)
    test_output <- evaluate_fn(test_coords)
    stopifnot(is.matrix(test_output))
    k <- ncol(test_output)
  }

  structure(
    list(
      evaluate_fn = evaluate_fn,
      k = k,
      crs = sf::st_crs(crs),
      geometry = NULL,
      type = type,
      metadata = metadata
    ),
    class = c("function_basis", "basis")
  )
}

#' Constant Basis Function
#'
#' Returns a spatial basis that evaluates to 1 for all locations.
#' Used for testing performance of integration routines.
#'
#' @param k Number of basis functions. Defaults to 1.
#' @param crs Optional coordinate reference system.
#' @param metadata Optional list of metadata.
#'
#' @return A `basis` object with class `constant_basis` and `basis`.
#' @export
make_constant_basis <- function(k = 1, crs = NULL, metadata = list()) {
  evaluate_fn <- function(coords) {
    matrix(1, nrow = nrow(coords), ncol = k)
  }

  structure(
    list(
      evaluate_fn = evaluate_fn,
      k = k,
      crs = crs,
      type = "constant",
      metadata = metadata
    ),
    class = c("constant_basis", "basis")
  )
}

#' Construct a Sine/Cosine Spatial Basis Matching Matérn SPDE Eigenstructure
#' with Automatic Square Embedding
#'
#' Creates a cosine basis over a square embedding of the input rectangle,
#' sorted by decreasing eigenvalue.
#'
#' @param max_frequency Integer: maximum frequency index in each dimension.
#' @param kappa Positive numeric: kappa parameter controlling range.
#' @param alpha Positive numeric: alpha parameter controlling smoothness.
#' @param K Integer: number of basis functions to keep.
#' @param domain Numeric vector length 4: c(xmin, xmax, ymin, ymax).
#' @param crs Optional CRS.
#' @param metadata Optional metadata list.
#'
#' @return A basis object of class "function_basis".
#' @export
#' Construct a Sine/Cosine Spatial Basis Matching Matérn SPDE Eigenstructure
#' with Automatic Square Embedding
#'
#' Creates a cosine basis over a square embedding of the input rectangle,
#' sorted by decreasing Matérn covariance eigenvalue.
#'
#' @param max_frequency Integer: maximum frequency index in each dimension.
#' @param kappa Positive numeric: kappa parameter controlling range (used only for sorting).
#' @param alpha Positive numeric: alpha parameter controlling smoothness (used only for sorting).
#' @param K Integer: number of basis functions to keep.
#' @param domain Numeric vector length 4: c(xmin, xmax, ymin, ymax).
#' @param crs Optional CRS.
#' @param metadata Optional metadata list.
#'
#' @return A basis object of class "function_basis".
#' @export
make_matern_fourier_basis <- function(
    max_frequency = 30,
    kappa = 1,
    alpha = 1,
    K = 200,
    domain = c(0,1,0,1),
    crs = NULL,
    metadata = list()
) {
  stopifnot(length(domain)==4)
  xmin0 <- domain[1]
  xmax0 <- domain[2]
  ymin0 <- domain[3]
  ymax0 <- domain[4]

  dx <- xmax0 - xmin0
  dy <- ymax0 - ymin0
  L <- max(dx, dy)

  xmin <- xmin0
  xmax <- xmin0 + L
  ymin <- ymin0
  ymax <- ymin0 + L

  scale_x <- xmax - xmin
  scale_y <- ymax - ymin

  # All frequency combinations
  kx_all <- rep(0:max_frequency, each = max_frequency + 1)
  ky_all <- rep(0:max_frequency, times = max_frequency + 1)
  #keep <- !(kx_all==0 & ky_all==0)
  kx <- kx_all#[keep]
  ky <- ky_all#[keep]

  # Laplacian eigenvalues
  laplacian_lambda <- pi^2*(kx^2 + ky^2)/L^2

  # For sorting only
  matern_lambda <- (kappa^2 + laplacian_lambda)^(-alpha)

  # Sort by decreasing covariance eigenvalue
  o <- order(matern_lambda, decreasing=TRUE)
  kx <- kx[o]
  ky <- ky[o]
  laplacian_lambda <- laplacian_lambda[o]
  matern_lambda <- matern_lambda[o]

  # Truncate
  kx <- kx[seq_len(K)]
  ky <- ky[seq_len(K)]
  laplacian_lambda <- laplacian_lambda[seq_len(K)]
  matern_lambda <- matern_lambda[seq_len(K)]

  evaluate_fn <- function(coords) {
    x_scaled <- (coords[,1] - xmin)/scale_x
    y_scaled <- (coords[,2] - ymin)/scale_y
    n <- nrow(coords)
    mat <- matrix(NA_real_, nrow=n, ncol=K)
    for(i in seq_len(K)){
      fx <- if(kx[i]==0) rep(1,n) else cos(kx[i]*pi*x_scaled)
      fy <- if(ky[i]==0) rep(1,n) else cos(ky[i]*pi*y_scaled)
      mat[,i] <- fx * fy
    }
    mat
  }

  make_function_basis(
    evaluate_fn = evaluate_fn,
    k = K,
    crs = crs,
    type = "matern_fourier",
    metadata = c(
      list(
        kappa_for_sorting = kappa,
        alpha_for_sorting = alpha,
        kx = kx,
        ky = ky,
        laplacian_lambda = laplacian_lambda,
        matern_lambda_for_sorting = matern_lambda,
        domain_rectangle = domain,
        domain_square = c(xmin, xmax, ymin, ymax),
        side_length = L
      ),
      metadata
    )
  )
}

#' Project a Covariate onto a Spatial Field Basis
#'
#' Projects a covariate observed over geometry onto the basis defined by a
#' `spatial_field()` object, returning a basis-like object with fitted scores
#' and evaluation support. Intended for fixed-effect covariates.
#'
#' @param spatial_field A `spatial_field` object with basis and prior.
#' @param covariate_data An `sf` object with covariate values and geometry.
#' @param covariate_name Name of the covariate column in `covariate_data`.
#' @param ... Additional arguments passed to `prepare_spatial_data()`.
#'
#' @return A `basis`-like object of class `covariate_basis`, `basis`, and the original basis class,
#'         with fields:
#' \describe{
#'   \item{evaluate_fn}{A function to evaluate the smoothed covariate at new coordinates.}
#'   \item{k}{The number of basis functions.}
#'   \item{crs}{The coordinate reference system.}
#'   \item{scores}{The fitted coefficients (posterior mean from ridge or flat prior).}
#'   \item{basis}{The original basis object used for fitting.}
#'   \item{label}{The covariate name.}
#'   \item{field_type}{Always `"fixed"` for projected covariates.}
#'   \item{metadata}{Optional metadata, taken from the original basis.}
#' }
#'
#' @export
make_covariate_basis <- function(spatial_field,
                                 covariate_data,
                                 covariate_name,
                                 ...) {
  stopifnot(inherits(spatial_field, "spatial_field"))
  stopifnot(inherits(covariate_data, "sf"))
  stopifnot(covariate_name %in% names(covariate_data))

  basis <- spatial_field$basis
  prior <- spatial_field$prior

  if (is.null(basis)) stop("`spatial_field` must include a basis.")
  if (is.null(prior)) stop("`spatial_field` must include a prior (flat for fixed effects).")

  prepared <- prepare_spatial_data(
    spatial_field = spatial_field,
    observation_data = covariate_data,
    response = covariate_name,
    ...
  )

  fit <- fit_spatial_field(prepared)

  evaluate_fn <- function(coords) {
    Phi <- basis$evaluate_fn(coords)
    drop(Phi %*% fit$posterior_mean)
  }

  structure(
    list(
      evaluate_fn = evaluate_fn,
      k = basis$k,
      crs = basis$crs,
      geometry = basis$geometry,
      scores = fit$posterior_mean,
      basis = basis,
      label = covariate_name,
      field_type = "fixed",
      metadata = basis$metadata
    ),
    class = c("covariate_basis", class(basis))
  )
}





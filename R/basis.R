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
#' Creates an orthonormal cosine basis over a square embedding of the input rectangle,
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
    metadata = list(),
    include_constant = TRUE
) {
  stopifnot(length(domain) == 4)
  xmin0 <- domain[1]
  xmax0 <- domain[2]
  ymin0 <- domain[3]
  ymax0 <- domain[4]

  dx <- xmax0 - xmin0
  dy <- ymax0 - ymin0
  L <- max(dx, dy)
  cx <- (xmin0 + xmax0) / 2
  cy <- (ymin0 + ymax0) / 2

  xmin <- cx - L/2
  xmax <- cx + L/2
  ymin <- cy - L/2
  ymax <- cy + L/2

  scale_x <- xmax - xmin
  scale_y <- ymax - ymin

  # Frequency grid
  kx_all <- rep(0:max_frequency, each = max_frequency + 1)
  ky_all <- rep(0:max_frequency, times = max_frequency + 1)
  n_all <- length(kx_all)

  # Laplacian eigenvalues and Matérn eigenvalues
  laplacian_lambda <- pi^2 * (kx_all^2 + ky_all^2)
  matern_lambda <- (kappa^2 + laplacian_lambda)^(-alpha)

  # Sort by decreasing covariance eigenvalue
  o <- order(matern_lambda, decreasing = TRUE)
  kx <- kx_all[o]
  ky <- ky_all[o]
  laplacian_lambda <- laplacian_lambda[o]
  matern_lambda <- matern_lambda[o]

  # Truncate
  kx <- kx[seq_len(K)]
  ky <- ky[seq_len(K)]
  laplacian_lambda <- laplacian_lambda[seq_len(K)]
  matern_lambda <- matern_lambda[seq_len(K)]

  # Evaluation function with orthonormal scaling
  evaluate_fn <- function(coords) {
    x_scaled <- (coords[,1] - xmin) / scale_x
    y_scaled <- (coords[,2] - ymin) / scale_y
    n <- nrow(coords)
    mat <- matrix(NA_real_, nrow = n, ncol = K)
    for(i in seq_len(K)) {
      fx <- if (kx[i] == 0) rep(1, n) else sqrt(2) * cos(kx[i] * pi * x_scaled)
      fy <- if (ky[i] == 0) rep(1, n) else sqrt(2) * cos(ky[i] * pi * y_scaled)
      mat[, i] <- fx * fy
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
        side_length = L,
        constant_mode_included = include_constant,
        constant_mode_precision = kappa^(2 * alpha)  # could override to (1 - rho)^2 if desired
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


## NYSTROM ADDON ###############################################################

# Isotropic Matérn covariance (nu=1.5 example)
make_matern_cov_fn <- function(range = 1, sigma2 = 1, nu = 1.5) {
  function(x1, x2) {
    stopifnot(ncol(x1) == ncol(x2))  # allow any dimension
    dists <- fields::rdist(x1, x2)   # returns [nrow(x1), nrow(x2)] matrix

    if (nu == 0.5) {
      sigma2 * exp(-dists / range)
    } else if (nu == 1.5) {
      scaled <- sqrt(3) * dists / range
      sigma2 * (1 + scaled) * exp(-scaled)
    } else if (nu == 2.5) {
      scaled <- sqrt(5) * dists / range
      sigma2 * (1 + scaled + scaled^2 / 3) * exp(-scaled)
    } else {
      stop("Currently only ν = 0.5, 1.5, 2.5 supported")
    }
  }
}


#' Construct a Nyström Basis Using a Covariance Function and Landmark Points
#'
#' Builds a set of approximate eigenfunctions of a Gaussian process covariance operator
#' using the Nyström method. The basis functions are smooth and globally supported,
#' formed as linear combinations of kernel evaluations at landmark locations.
#'
#' You may optionally normalize the basis functions using L2 normalization to improve
#' numerical stability and interpretability. This preserves the diagonal structure
#' of the prior covariance and avoids basis rotations.
#'
#' If a `domain` is provided, landmark and evaluation coordinates are rescaled to
#' \[0, 1\]^2 to match the reference domain used when constructing the basis.
#'
#' @param cov_fn A function taking two matrices `x1` and `x2` (each of shape `[n, d]`)
#'   and returning a `[nrow(x1), nrow(x2)]` covariance matrix. Must support vectorized input.
#' @param landmarks A numeric matrix `[m, d]` of landmark locations.
#' @param epsilon Small scalar added to the diagonal of `K_UU` for numerical stability.
#' @param k Optional. Number of basis functions to retain (defaults to all available).
#' @param crs Optional coordinate reference system (passed to output `function_basis`).
#' @param metadata Optional list of metadata for downstream use.
#' @param domain Optional length-4 numeric vector `c(xmin, xmax, ymin, ymax)` for rescaling
#'   both landmarks and evaluation coordinates to the unit square.
#' @param normalize Character string specifying normalization strategy.
#'   Must be one of `"none"` (default) or `"l2"`.
#'   - `"none"`: return raw Nyström basis (scaled as-is).
#'   - `"l2"`: empirically normalize basis functions to have unit L2 norm over the unit square.
#' @param normalization_points Number of sample points used for computing L2 norms (only used if `normalize = "l2"`).
#'
#' @return A `function_basis` object with class `"nystrom"`, containing:
#'   - An `evaluate_fn` for computing the basis matrix at new locations,
#'   - The number of basis functions `k`,
#'   - The original domain, landmarks, eigenvectors/values, and scaling metadata.
#'
#' @export

make_nystrom_basis <- function(cov_fn,
                               landmarks,
                               epsilon = 1e-6,
                               k = NULL,
                               crs = NULL,
                               metadata = list(),
                               domain = NULL,
                               normalize = c("none","l2"),
                               normalization_points = 1e4) {
  normalize <- match.arg(normalize)
  stopifnot(is.function(cov_fn))
  stopifnot(is.matrix(landmarks))
  m <- nrow(landmarks)

  if (!is.null(domain)) {
    stopifnot(length(domain) == 4)
    xmin <- domain[1]
    xmax <- domain[2]
    ymin <- domain[3]
    ymax <- domain[4]
    scale_x <- xmax - xmin
    scale_y <- ymax - ymin

    rescale <- function(coords) {
      coords <- as.matrix(coords)
      cbind((coords[, 1] - xmin) / scale_x,
            (coords[, 2] - ymin) / scale_y)
    }

    landmarks <- rescale(landmarks)
  } else {
    rescale <- identity
  }

  K_UU <- cov_fn(landmarks, landmarks)
  K_UU <- K_UU + epsilon * diag(m)

  message("Computing eigendecomposition for Nyström basis...")
  eig <- eigen(K_UU, symmetric = TRUE)
  Lambda <- eig$values
  V <- eig$vectors

  keep <- if (is.null(k)) length(Lambda) else min(k, length(Lambda))
  Lambda_k <- Lambda[1:keep]
  V_k <- V[, 1:keep, drop = FALSE]

  if (any(Lambda_k <= 0)) stop("Non-positive eigenvalues encountered. Try increasing epsilon.")

  evaluate_fn <- function(coords) {
    coords_scaled <- rescale(coords)
    K_SU <- cov_fn(coords_scaled, landmarks)
    K_SU %*% V_k
  }

  if (normalize == "l2") {
    message("Applying L2 normalization")
    set.seed(1)
    grid_pts <- matrix(runif(2 * normalization_points), ncol = 2)

    K_SU <- cov_fn(grid_pts, landmarks)  # skip rescaling; grid_pts already in [0,1]^2
    Phi <- K_SU %*% V_k

    norms <- sqrt(colMeans(Phi^2))
    norms[norms == 0] <- 1
    evaluate_fn <- function(coords) {
      coords_scaled <- rescale(coords)
      K_SU <- cov_fn(coords_scaled, landmarks)
      sweep(K_SU %*% V_k, 2, norms, "/")
    }
  }

  make_function_basis(
    evaluate_fn = evaluate_fn,
    k = keep,
    crs = crs,
    type = "nystrom",
    metadata = c(
      list(
        landmarks = landmarks,
        epsilon = epsilon,
        eigenvalues = Lambda_k,
        eigenvectors = V_k,
        K_UU = K_UU,
        original_domain = domain,
        normalization = normalize
      ),
      metadata
    )
  )
}

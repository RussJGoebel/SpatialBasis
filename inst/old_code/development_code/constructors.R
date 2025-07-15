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

#' Create a Gaussian prior over basis coefficients
#'
#' Defines a Gaussian prior on the coefficients of a spatial basis expansion.
#' Users must supply either:
#' - A `variances` vector (to construct a diagonal sparse precision matrix), or
#' - A full `precision` matrix (can be dense or sparse; not necessarily invertible).
#'
#' Internally, precision matrices are coerced to `Matrix` class (sparse by default)
#' for consistency and downstream compatibility with sparse solvers.
#'
#' @param basis A `basis` object representing the basis functions
#' @param precision Optional `k × k` precision matrix (base or Matrix class)
#' @param variances Optional numeric vector of length `k` (strictly positive)
#' @param hyperparams Optional list of hyperparameters (e.g., range, nu)
#' @param metadata Optional list of metadata (e.g., kernel, notes)
#'
#' @return A `prior_on_basis` object with:
#' \describe{
#'   \item{basis}{The associated `basis` object}
#'   \item{precision}{A sparse Matrix-class precision matrix}
#'   \item{variances}{Vector of prior marginal variances (or `NULL` if not applicable)}
#'   \item{hyperparams}{List of optional hyperparameters}
#'   \item{metadata}{List of user-defined metadata}
#' }
#'
#' @importFrom Matrix Diagonal Matrix isSymmetric
#' @export
make_prior_on_basis <- function(basis,
                                precision = NULL,
                                variances = NULL,
                                hyperparams = list(),
                                metadata = list()) {
  if (!inherits(basis, "basis")) {
    stop("`basis` must be a `basis` object.")
  }

  k <- basis$k

  if (is.null(variances) && is.null(precision)) {
    stop("You must provide either `variances` or `precision`.")
  }

  # Case 1: Variances → build diagonal precision
  if (!is.null(variances)) {
    if (!is.numeric(variances) || length(variances) != k || any(variances <= 0)) {
      stop("`variances` must be a numeric vector of length k with strictly positive values.")
    }
    precision <- Matrix::Diagonal(x = 1 / variances)
  }

  # Case 2: Precision → validate and coerce
  if (!is.null(precision)) {
    if (!is.matrix(precision) && !inherits(precision, "Matrix")) {
      stop("`precision` must be a base matrix or a Matrix object.")
    }
    if (any(dim(precision) != k)) {
      stop("`precision` must be a square matrix of size k × k.")
    }

    # Allow for minor numerical asymmetries from construction
    if (!Matrix::isSymmetric(precision, tol = 1e-8)) {
      stop("`precision` must be symmetric.")
    }
    if (!inherits(precision, "Matrix")) {
      precision <- Matrix::Matrix(precision, sparse = TRUE)
    }


  }

  structure(
    list(
      basis = basis,
      precision = precision,
      variances = variances,
      hyperparams = hyperparams,
      metadata = metadata
    ),
    class = "prior_on_basis"
  )
}


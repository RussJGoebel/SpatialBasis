#' Construct a generic Basis object
#'
#' A Basis object provides a function to evaluate a set of spatial basis functions
#' at arbitrary coordinates. The basis can be defined however the user chooses
#' (analytic, pixel-based, derived), as long as the evaluate_fn returns an n × k matrix.
#'
#' @param evaluate_fn A function taking coordinates (matrix, data.frame, or sf POINT)
#'        and returning an [n × k] matrix of basis evaluations.
#' @param k Number of basis functions
#' @param metadata Optional list of descriptive info (e.g., type, kernel, support)
#'
#' @return An object of class "Basis"
#' @export
make_basis <- function(evaluate_fn, k, metadata = list()) {
  if (!is.function(evaluate_fn)) {
    stop("`evaluate_fn` must be a function that returns a matrix of basis evaluations.")
  }
  if (!is.numeric(k) || length(k) != 1 || k < 1) {
    stop("`k` must be a positive integer.")
  }

  structure(
    list(
      evaluate = evaluate_fn,
      k = k,
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
    if (!Matrix::isSymmetric(precision, tol = 1e-8)) {
      stop("`precision` must be symmetric.")
    }
    if (!inherits(precision, "Matrix")) {
      precision <- Matrix::Matrix(precision, sparse = TRUE)
    }

    # Drop variances unless they were user-supplied
    if (is.null(variances)) {
      variances <- NULL
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


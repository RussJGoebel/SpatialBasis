#' Construct a prior over pixel-level basis functions
#'
#' This function creates a prior object assuming each fine-resolution grid cell
#' has a distinct weight (i.e., a pixel-level basis). The prior precision can
#' be diagonal (independent) or structured via a user-supplied precision matrix (e.g., SAR).
#'
#' @param locations An sf object with one row per pixel; determines the spatial support.
#' @param precision Optional precision matrix (e.g., SAR precision). Default: identity.
#' @param tau Optional scalar multiplier for the precision matrix (default = 1).
#' @param integration_method Method for computing inner products (default = "area").
#'
#' @return An object of class `downscaling_prior`
#' @export
make_pixel_prior <- function(locations,
                             precision = NULL,
                             tau = 1,
                             integration_method = "area") {
  if (!inherits(locations, "sf")) stop("`locations` must be an sf object.")
  n <- nrow(locations)

  # Default to diagonal precision if none given
  if (is.null(precision)) {
    precision <- Matrix::Diagonal(n, tau)
  } else {
    if (!inherits(precision, "Matrix")) stop("`precision` must be a Matrix object.")
    precision <- tau * precision
  }

  # Define a pixel identity basis: predict() returns identity matrix
  identity_basis <- make_basis(
    predict_fn = function(locs) {
      if (nrow(locs) != n) stop("Pixel basis expects locations with same number of rows as original.")
      diag(n)
    },
    k = n,
    metadata = list(type = "pixel", original_locations = locations)
  )

  make_downscaling_prior(
    type = "pixel",
    basis_constructor = NULL,
    basis_args = list(),
    k = n,
    precision = precision,
    covariance = NULL,
    hyperparams = list(tau = tau),
    metadata = list(locations = locations, integration_method = integration_method)
  )
}

#' Create a functional prior for basis weight regularization
#'
#' This constructor defines a prior for a functional basis model, such as one using
#' eigenfunctions or KL expansions. It supports a flexible number of basis functions,
#' and allows specification of a prior precision or covariance structure.
#'
#' @param basis An object representing the basis; must support evaluation (e.g., via a `predict()` method).
#' @param k Integer truncation level (number of basis functions to retain).
#' @param precision A k x k prior precision matrix.
#' @param integration_method Character string; e.g., `"qmc"` or `"analytic"`.
#' @param hyperparameters Optional named list of prior hyperparameters (e.g., `list(tau = 1)`).
#'
#' @return An object of class `downscaling_prior` with `type = "functional"`.
#' @export
make_functional_prior <- function(basis,
                                  k,
                                  precision,
                                  integration_method = "qmc",
                                  hyperparameters = list()) {
  if (!is.numeric(k) || length(k) != 1 || k <= 0) stop("`k` must be a positive integer.")
  if (!(is.matrix(precision) || inherits(precision, "Matrix"))) stop("`precision` must be a matrix.")
  if (any(dim(precision) != k)) stop("`precision` must be a k x k matrix.")

  make_downscaling_prior(
    type = "functional",
    k = k,
    basis = basis,
    precision = precision,
    integration_method = integration_method,
    hyperparameters = hyperparameters,
    metadata = list(basis_class = class(basis))
  )
}

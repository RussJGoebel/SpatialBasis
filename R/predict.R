#' Predict from a fitted spatial field (single or stacked)
#'
#' Estimates the posterior mean at new coordinates for one or more spatial fields.
#'
#' @param object A fitted_spatial_field object (from either `fit_spatial_field()` or `fit_stacked_spatial_model()`).
#' @param newdata Matrix or data.frame of coordinates (n x 2), columns (x, y).
#' @param ... Unused.
#'
#' @return A numeric vector of predicted values.
#'
#' @export
predict.fitted_spatial_field <- function(object, newdata, ...) {
  stopifnot(inherits(object, "fitted_spatial_field"))

  if (missing(newdata) || is.null(newdata)) {
    stop("You must provide `newdata` as a matrix or data.frame of coordinates.")
  }

  coords <- as.matrix(newdata)
  if (ncol(coords) != 2) {
    stop("`newdata` must have 2 columns (x and y coordinates).")
  }

  # Handle single-field models
  if (!is.null(object$spatial_data)) {
    basis <- object$spatial_data$basis
    phi <- basis$evaluate_fn(coords)
    preds <- as.numeric(phi %*% object$posterior_mean)
    return(preds)
  }

  # Handle stacked models
  if (!is.null(object$fields) && !is.null(object$column_map)) {
    preds <- numeric(nrow(coords))

    for (i in seq_along(object$fields)) {
      field <- object$fields[[i]]
      label <- field$label %||% paste0("V", i)

      cols <- object$column_map[[i]]
      beta_i <- object$posterior_mean[cols]

      phi <- field$basis$evaluate_fn(coords)
      preds <- preds + as.numeric(phi %*% beta_i)
    }

    return(preds)
  }

  stop("Invalid fitted_spatial_field object: missing `spatial_data` or `fields`.")
}

#' Predict region-wise averages from a fitted spatial field
#'
#' Computes expected values over regions by integrating the spatial field basis functions.
#'
#' @param object A `fitted_spatial_field` object.
#' @param newregions An `sf` object (e.g. polygons) for which to predict region-averaged values.
#' @param integration_rule Optional integration rule to override the default.
#' @param parallel Logical or integer: whether/how to parallelize the integration.
#' @param ... Unused.
#'
#' @return A numeric vector of predicted region means.
#' @export
predict_regionwise.fitted_spatial_field <- function(
    object,
    newregions,
    integration_rule = NULL,
    parallel = FALSE,
    ...
) {
  stopifnot(inherits(object, "fitted_spatial_field"))
  stopifnot(inherits(newregions, "sf"))

  region_list <- build_region_list(newregions)

  if (!is.null(object$spatial_data)) {
    basis <- object$spatial_data$basis
    rule <- integration_rule %||% object$spatial_data$integration_rule
    X_new <- integrate_region_list_parallel(region_list, rule, basis,parallel = parallel,...)
    preds <- as.numeric(X_new %*% object$posterior_mean)
    return(preds)
  }

  if (!is.null(object$fields) && !is.null(object$column_map)) {
    preds <- numeric(length(region_list))

    for (i in seq_along(object$fields)) {
      field <- object$fields[[i]]
      label <- field$label %||% paste0("V", i)
      cols <- object$column_map[[label]]
      beta_i <- object$posterior_mean[cols]
      rule <- integration_rule %||% field$integration_rule
      X_i <- integrate_region_list_parallel(region_list, rule, field$basis,...)
      preds <- preds + as.numeric(X_i %*% beta_i)
    }

    return(preds)
  }

  stop("Invalid fitted_spatial_field object: missing `spatial_data` or `fields`.")
}

#' Predict from a fitted spatial field (coordinates or regions)
#'
#' Estimates the posterior mean at new coordinates (pointwise) or over spatial regions (averaged).
#'
#' @param object A `fitted_spatial_field` object.
#' @param newdata Either a matrix/data.frame of coordinates (n x 2), or an `sf` object of prediction regions.
#' @param integration_rule Optional: override the integration rule (only for regionwise).
#' @param parallel Logical or integer: whether/how to parallelize integration (only for regionwise).
#' @param ... Unused.
#'
#' @return A numeric vector of predicted values.
#'
#' @export
predict.fitted_spatial_field <- function(object, newdata,
                                         integration_rule = NULL,
                                         parallel = FALSE,
                                         n_workers = future::availableCores(),
                                         ...) {
  stopifnot(inherits(object, "fitted_spatial_field"))

  # Case 1: regionwise prediction (sf polygons)
  if (inherits(newdata, "sf")) {
    return(predict_regionwise.fitted_spatial_field(
      object,
      newregions = newdata,
      integration_rule = integration_rule,
      parallel = parallel,
      n_workers = n_workers,
      ...
    ))
  }

  # Case 2: pointwise prediction (x/y coordinates)
  coords <- as.matrix(newdata)
  if (ncol(coords) != 2) {
    stop("`newdata` must have 2 columns (x and y coordinates).")
  }

  # Single spatial field
  if (!is.null(object$spatial_data)) {
    phi <- object$spatial_data$basis$evaluate_fn(coords)
    return(as.numeric(phi %*% object$posterior_mean))
  }

  # Stacked fields
  if (!is.null(object$fields) && !is.null(object$column_map)) {
    preds <- numeric(nrow(coords))
    for (i in seq_along(object$fields)) {
      field <- object$fields[[i]]
      label <- field$label %||% paste0("V", i)
      cols <- object$column_map[[label]]
      beta_i <- object$posterior_mean[cols]
      phi <- field$basis$evaluate_fn(coords)
      preds <- preds + as.numeric(phi %*% beta_i)
    }
    return(preds)
  }

  stop("Invalid fitted_spatial_field object: missing `spatial_data` or `fields`.")
}

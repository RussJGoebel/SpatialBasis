#' Predict region-wise averages from a fitted spatial field
#'
#' Computes expected values over regions by integrating the spatial field basis functions.
#'
#' @param object A `fitted_spatial_field` object.
#' @param newregions An `sf` object (e.g. polygons) for which to predict region-averaged values.
#' @param integration_rule Optional integration rule to override the default.
#' @param parallel Logical or integer: whether/how to parallelize the integration.
#' @param return_variance Logical: if TRUE, also return prediction variance.
#' @param diagonal_only Logical: if TRUE, only compute the diagonal entries of the posterior.
#' @param ... Passed to integration or variance functions.
#'
#' @return list(preds = <vector>, posterior_variance = <vector or matrix, optional>)
#' @export
predict_regionwise.fitted_spatial_field <- function(
    object,
    newregions,
    integration_rule = NULL,
    parallel = FALSE,
    return_variance = FALSE,
    diagonal_only = TRUE,
    ...
) {
  stopifnot(inherits(object, "fitted_spatial_field"))
  stopifnot(inherits(object$spatial_data, "prepared_spatial_data") ||
              inherits(object$spatial_data, "stacked_spatial_data"))
  stopifnot(inherits(newregions, "sf"))

  region_list <- build_region_list(newregions)

  ## ------------------------ Single-field case ------------------------
  if (inherits(object$spatial_data, "prepared_spatial_data")) {
    basis <- object$spatial_data$basis
    rule  <- if (is.null(integration_rule)) object$spatial_data$integration_rule else integration_rule
    X_new <- integrate_region_list_parallel(region_list, rule, basis)
    preds <- as.numeric(X_new %*% object$posterior_mean)

    if (!return_variance) return(list(preds = preds))

    # Coefficient posterior covariance
    Vbeta <- object$posterior_variance
    if (is.null(Vbeta) || is.null(dim(Vbeta))) {
      Vbeta <- compute_posterior_variance(
        object$XtX + object$spatial_data$Lambda,
        cholesky_factor = object$cholesky_factor,
        diagonal_only   = FALSE
      )
      obs_var_mode <- compute_posterior_measurement_error(
        object$posterior_mean, object$X, object$Lambda, object$y
      )$sigma2_mode
      Vbeta <- Vbeta * obs_var_mode
    }

    XV  <- X_new %*% Vbeta
    Var <- tcrossprod(XV, X_new)  # X V X'
    if (diagonal_only) Var <- rowSums(X_new * (XV))  # diag(X V X') efficiently

    return(list(preds = preds, posterior_variance = Var))
  }

  ## ------------------------ Stacked (multi-field) ------------------------
  if (inherits(object$spatial_data, "stacked_spatial_data")) {
    stopifnot(!is.null(object$spatial_data$fields), !is.null(object$column_map))

    # 1) Build per-field region-averaged design
    X_block_list <- vector("list", length(object$spatial_data$fields))

    # Try to get labels from fields; fallback to column_map names
    field_names <- vapply(object$spatial_data$fields, function(f) {
      if (!is.null(f$label)) f$label else ""
    }, character(1))
    if (all(field_names == "")) field_names <- names(object$column_map)

    for (i in seq_along(object$spatial_data$fields)) {
      field <- object$spatial_data$fields[[i]]
      rule_i <- if (is.null(integration_rule)) field$integration_rule else integration_rule
      X_block_list[[i]] <- integrate_region_list_parallel(region_list, rule_i, field$basis)
    }
    names(X_block_list) <- field_names

    # 2) Reorder blocks to match training column_map order
    desired_order <- names(object$column_map)
    X_block_list  <- X_block_list[desired_order]

    # 3) Residualize target block using training A (reuse-A)
    A <- object$spatial_data$A
    if (!is.null(A)) {
      tgt_lab <- object$target_label
      stopifnot(tgt_lab %in% names(object$column_map))

      X_target_star <- X_block_list[[tgt_lab]]
      X_other_star  <- do.call(cbind, X_block_list[names(X_block_list) != tgt_lab])

      # Sanity: A dims (#other_cols x #target_cols)
      other_cols  <- unlist(object$column_map[names(object$column_map) != tgt_lab])
      target_cols <- object$column_map[[tgt_lab]]
      stopifnot(nrow(A) == length(other_cols), ncol(A) == length(target_cols))

      X_target_star_perp <- X_target_star - X_other_star %*% A
      X_block_list[[tgt_lab]] <- X_target_star_perp
    }

    # 4) Reassemble design in training column order
    X_new <- do.call(cbind, X_block_list)

    # 5) Predict mean
    preds <- as.numeric(X_new %*% object$posterior_mean)

    if (!return_variance) return(list(preds = preds))

    # Coefficient posterior covariance
    Vbeta <- object$posterior_variance
    if (is.null(Vbeta) || is.null(dim(Vbeta))) {
      Vbeta <- compute_posterior_variance(
        object$XtX + object$spatial_data$Lambda,
        cholesky_factor = object$cholesky_factor,
        diagonal_only   = FALSE
      )
      obs_var_mode <- compute_posterior_measurement_error(
        object$posterior_mean, object$X, object$Lambda, object$y
      )$sigma2_mode
      Vbeta <- Vbeta * obs_var_mode
    }

    XV  <- X_new %*% Vbeta
    Var <- tcrossprod(XV, X_new)  # X V X'
    if (diagonal_only) Var <- rowSums(X_new * (XV))

    return(list(preds = preds, posterior_variance = Var))
  }
}

#' Predict from a fitted spatial field (coordinates or regions)
#'
#' If `newdata` is an `sf` object, computes regional predictions via integration.
#' If `newdata` is a coordinate matrix, computes pointwise predictions.
#'
#' @param object A `fitted_spatial_field` object.
#' @param newdata Either a matrix/data.frame of coords (n x 2), or an `sf` object of regions.
#' @param return_variance Logical: if TRUE, return posterior prediction variance.
#' @param diagonal_only Logical: if TRUE, only compute diagonal of the variance matrix.
#' @param integration_rule Optional integration rule (regionwise only).
#' @param parallel Logical or integer: whether/how to parallelize integration (regionwise only).
#' @param n_workers Parallel workers (regionwise only).
#' @param ... Additional args.
#'
#' @return list(preds = <vector>, posterior_variance = <vector or matrix, optional>)
#' @export
predict.fitted_spatial_field <- function(object, newdata,
                                         return_variance = FALSE,
                                         diagonal_only = TRUE,
                                         integration_rule = NULL,
                                         parallel = FALSE,
                                         n_workers = future::availableCores(),
                                         ...) {
  stopifnot(inherits(object, "fitted_spatial_field"))

  # Regionwise path delegates
  if (inherits(newdata, "sf")) {
    return(predict_regionwise.fitted_spatial_field(
      object          = object,
      newregions      = newdata,
      return_variance = return_variance,
      diagonal_only   = diagonal_only,
      integration_rule = integration_rule,
      parallel        = parallel,
      n_workers       = n_workers,
      ...
    ))
  }

  # Pointwise predictions
  coords <- as.matrix(newdata)
  if (ncol(coords) != 2 || !is.numeric(coords)) {
    stop("`newdata` must be a numeric matrix/data.frame with 2 columns (x, y).")
  }

  ## ------------------------ Single-field case ------------------------
  if (inherits(object$spatial_data, "prepared_spatial_data")) {
    phi   <- object$spatial_data$basis$evaluate_fn(coords)
    preds <- as.numeric(phi %*% object$posterior_mean)

    if (!return_variance) return(list(preds = preds))

    Vbeta <- object$posterior_variance
    if (is.null(Vbeta) || is.null(dim(Vbeta))) {
      Vbeta <- compute_posterior_variance(
        object$XtX + object$Lambda,
        cholesky_factor = object$cholesky_factor,
        diagonal_only   = FALSE
      )
      obs_var_mode <- compute_posterior_measurement_error(
        object$posterior_mean, object$X, object$Lambda, object$y
      )$sigma2_mode
      Vbeta <- Vbeta * obs_var_mode
    }

    PV  <- phi %*% Vbeta
    Var <- tcrossprod(PV, phi)  # phi V phi'
    if (diagonal_only) Var <- rowSums(phi * (PV))

    return(list(preds = preds, posterior_variance = Var))
  }

  ## ------------------------ Stacked (multi-field) ------------------------
  if (inherits(object$spatial_data, "stacked_spatial_data")) {
    # Build per-field design at new coords
    phi_block_list <- vector("list", length(object$spatial_data$fields))

    field_names <- vapply(object$spatial_data$fields, function(f) {
      if (!is.null(f$label)) f$label else ""
    }, character(1))
    if (all(field_names == "")) field_names <- names(object$column_map)

    for (i in seq_along(object$spatial_data$fields)) {
      field <- object$spatial_data$fields[[i]]
      phi_block_list[[i]] <- field$basis$evaluate_fn(coords)
    }
    names(phi_block_list) <- field_names

    # Reorder to column_map order
    desired_order  <- names(object$column_map)
    phi_block_list <- phi_block_list[desired_order]

    # Residualize target block using reuse-A
    A <- object$spatial_data$A
    if (!is.null(A)) {
      tgt_lab <- object$target_label
      stopifnot(tgt_lab %in% names(object$column_map))

      Phi_target_star <- phi_block_list[[tgt_lab]]
      Phi_other_star  <- do.call(cbind, phi_block_list[names(phi_block_list) != tgt_lab])

      other_cols  <- unlist(object$column_map[names(object$column_map) != tgt_lab])
      target_cols <- object$column_map[[tgt_lab]]
      stopifnot(nrow(A) == length(other_cols), ncol(A) == length(target_cols))

      Phi_target_star_perp <- Phi_target_star - Phi_other_star %*% A
      phi_block_list[[tgt_lab]] <- Phi_target_star_perp
    }

    phi_all <- do.call(cbind, phi_block_list)

    preds <- as.numeric(phi_all %*% object$posterior_mean)

    if (!return_variance) return(list(preds = preds))

    Vbeta <- object$posterior_variance
    if (is.null(Vbeta) || is.null(dim(Vbeta))) {
      Vbeta <- compute_posterior_variance(
        object$XtX + object$Lambda,
        cholesky_factor = object$cholesky_factor,
        diagonal_only   = FALSE
      )
      obs_var_mode <- compute_posterior_measurement_error(
        object$posterior_mean, object$X, object$Lambda, object$y
      )$sigma2_mode
      Vbeta <- Vbeta * obs_var_mode
    }

    PV  <- phi_all %*% Vbeta
    Var <- tcrossprod(PV, phi_all)  # phi_all V phi_all'
    if (diagonal_only) Var <- rowSums(phi_all * (PV))

    return(list(preds = preds, posterior_variance = Var))
  }
}

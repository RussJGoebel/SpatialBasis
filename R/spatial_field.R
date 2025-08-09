#' Specify a Spatial Random Field Term
#'
#' Creates an object representing a latent spatial random field for downscaling models.
#' You can optionally supply a basis, prior, and integration rule.
#' If any are left NULL, defaults will be inferred inside `downscale()`.
#'
#' @param basis Optional. A basis object created by [make_basis()], [make_geometry_basis()], etc.
#' @param prior Optional. A prior object created by e.g. [make_sar_prior()].
#' @param integration_rule Optional. An integration rule object (e.g., from [integration_rule_qmc_average()]).
#' @param label Optional. A descriptive label for this spatial field (e.g., "Random effect").
#'
#' @return An object of class `"spatial_field"`.
#'
#' @examples
#' # Empty (defaults filled later)
#' sf1 <- make_spatial_field()
#'
#' # With a geometry basis and SAR prior
#' grid <- sf::st_make_grid(sf::st_bbox(c(xmin=0, ymin=0, xmax=10, ymax=10)), n = c(5,5))
#' b <- make_geometry_basis(sf::st_sf(geometry=grid))
#' p <- make_sar_prior(b)
#' sf2 <- make_spatial_field(basis=b, prior=p, label="SAR spatial field")
#'
#' @export
make_spatial_field <- function(basis = NULL,
                               prior = NULL,
                               integration_rule = NULL,
                               label = "Spatial Field") {
  # Validate basis
  if (!is.null(basis) && !inherits(basis, "basis")) {
    stop("`basis` must be a basis object (from make_basis(), make_geometry_basis(), etc.).")
  }

  # Validate prior
  if (!is.null(prior) && !inherits(prior, "prior")) {
    stop("`prior` must be a prior object (e.g., from make_sar_prior()).")
  }

  # Use default integration rule if needed
  if (is.null(integration_rule) && !is.null(basis)) {
    integration_rule <- default_integration_rule(basis)
  }

  # Validate integration rule
  if (!is.null(integration_rule) && !inherits(integration_rule, "integration_rule")) {
    stop("`integration_rule` must be an integration rule object.")
  }

  # Validate label
  if (!is.character(label) || length(label) != 1) {
    stop("`label` must be a character string.")
  }

  structure(
    list(
      basis = basis,
      prior = prior,
      integration_rule = integration_rule,
      label = label
    ),
    class = "spatial_field"
  )
}


#' Prepare Spatial Field for Downscaling
#'
#' Integrates the basis functions of a spatial field over observation regions,
#' and optionally extracts the response variable.
#'
#' @param spatial_field A spatial_field() object (basis + prior + integration rule).
#' @param observation_data An sf object containing the observations, including the response column.
#' @param response Character string naming the response column in `observation_data`.
#' @param region_list Optional: precomputed list of regions (overrides auto-building).
#' @param integration_rule Optional: override the integration rule in the spatial_field.
#' @param parallel Logical or integer: whether/how to parallelize integration.
#' @param compute_cholesky Logical: whether to compute and return the Cholesky of XtX + Lambda.
#' @param ... Additional arguments passed to integration functions.
#'
#' @return A list with:
#'   - y: the response vector.
#'   - X_spatial: the integrated design matrix.
#'   - XtX: XᵗX matrix
#'   - Lambda: prior precision matrix
#'   - posterior_precision: XtX + Lambda
#'   - cholesky_factor: Cholesky of XtX + Lambda (if requested)
#'   - basis, prior, integration_rule, region_list, etc.
#'
#' @export
prepare_spatial_data <- function(
    spatial_field,
    observation_data,
    response,
    region_list = NULL,
    integration_rule = NULL,
    parallel = FALSE,
    compute_cholesky = FALSE,
    ...
) {
  stopifnot(inherits(spatial_field, "spatial_field"))
  stopifnot(inherits(observation_data, "sf"))
  stopifnot(response %in% names(observation_data))

  # Build region list if needed
  if (is.null(region_list)) {
    message("Building region list from observation_data geometries...")
    region_list <- build_region_list(observation_data)
  }

  # Decide integration rule
  rule <- if (!is.null(integration_rule)) integration_rule else spatial_field$integration_rule
  stopifnot(inherits(rule, "integration_rule"))

  message("Integrating basis functions over region...")
  X_spatial <- if (parallel) {
    integrate_region_list_parallel(region_list, rule, spatial_field$basis, ...)
  } else {
    integrate_region_list(region_list, rule, spatial_field$basis)
  }

  y <- observation_data[[response]]
  X <- Matrix::Matrix(X_spatial)

  message("Computing X'X...")
  XtX <- Matrix::crossprod(X)
  Lambda <- compute_precision(spatial_field$prior)

  if (is.null(Lambda)) {
    stop("Spatial field must have a prior to compute posterior precision.")
  }

  posterior_precision <- XtX + Lambda

  # Optional Cholesky computation
  cholesky_factor <- NULL
  if (compute_cholesky) {
    message("Computing Cholesky of posterior precision...")
    cholesky_factor <- try(Matrix::chol(posterior_precision), silent = TRUE)
    if (inherits(cholesky_factor, "try-error")) {
      warning("Cholesky failed: posterior precision is not positive definite.")
      cholesky_factor <- NULL
    }
  }

  out <- list(
    y = y,
    X_spatial = X,
    XtX = XtX,
    Lambda = Lambda,
    posterior_precision = posterior_precision,
    cholesky_factor = cholesky_factor,
    basis = spatial_field$basis,
    prior = spatial_field$prior,
    integration_rule = rule,
    region_list = region_list,
    observation_data = observation_data,
    label = spatial_field$label,
    response = response
  )
  class(out) <- "prepared_spatial_data"
  return(out)
}




#' Convert a mixture of sf objects and spatial_fields into a unified list of spatial_field objects
#'
#' @param data_list A list of `sf` objects and/or `spatial_field` objects
#' @param covariate_names A character vector of covariate names to extract from `sf` objects
#' @param default_prior A function that takes `k` and returns a prior object (e.g., ridge prior)
#'
#' @return A list of `spatial_field` objects, one per covariate or explicitly passed field
#' @export
as_spatial_field_list <- function(data_list,
                                  covariate_names,
                                  default_prior = function(basis) make_ridge_prior(basis, tau = 1e-4)
) {
  stopifnot(is.list(data_list))
  stopifnot(is.character(covariate_names))

  spatial_fields <- list()
  sf_objects <- list()

  # Separate out spatial_field() and sf objects
  for (item in data_list) {
    if (inherits(item, "spatial_field")) {
      spatial_fields[[length(spatial_fields) + 1]] <- item
    } else if (inherits(item, "sf")) {
      sf_objects[[length(sf_objects) + 1]] <- item
    } else {
      stop("Each element of `data_list` must be an sf object or a spatial_field.")
    }
  }

  # Extract already used labels
  used_labels <- vapply(spatial_fields, function(f) f$label, character(1))

  # Process each covariate not already supplied
  for (cov_name in covariate_names) {
    if (cov_name %in% used_labels) {
      message(sprintf("Covariate '%s' already supplied as spatial_field; skipping.", cov_name))
      next
    }

    # Search for the covariate in the sf_objects
    found <- FALSE
    for (sf_obj in sf_objects) {
      if (cov_name %in% names(sf_obj)) {
        values <- sf_obj[[cov_name]]

        # Sanity check: must be numeric or integer
        if (!is.numeric(values)) {
          stop(sprintf("Covariate '%s' must be numeric to create piecewise basis.", cov_name))
        }

        basis <- make_piecewise_constant_basis(sf_obj, values = values)
        prior <- default_prior(basis)

        spatial_fields[[length(spatial_fields) + 1]] <-
          make_spatial_field(label = cov_name, basis = basis, prior = prior)

        found <- TRUE
        break
      }
    }

    if (!found) {
      stop(sprintf("Covariate '%s' not found in any provided sf object.", cov_name))
    }
  }

  return(spatial_fields)
}

#' Prepare a List of Spatial Fields for Downscaling
#'
#' Applies `prepare_spatial_data()` to each spatial field in a list,
#' extracting the design matrix and response variable for each term.
#'
#' This is useful for additive models where multiple spatial fields (e.g., covariates,
#' structured priors, or random effects) contribute to the overall linear predictor.
#'
#' @param field_list A list of `spatial_field()` objects (e.g., from [as_spatial_field_list()]).
#' @param observation_data An `sf` object containing the response and covariates.
#' @param response Character name of the response variable in `observation_data`.
#' @param region_list Optional: shared precomputed region list (applied to all fields).
#' @param integration_rule Optional: override the integration rule for all fields.
#' @param parallel Logical or integer: whether/how to parallelize integration.
#' @param ... Additional arguments passed to `prepare_spatial_data()`.
#'
#' @return A list of `prepared_spatial_data` objects (one per field).
#' @export
prepare_spatial_field_list <- function(
    field_list,
    observation_data,
    response,
    region_list = NULL,
    integration_rule = NULL,
    parallel = FALSE,
    ...
) {
  stopifnot(is.list(field_list))
  stopifnot(all(vapply(field_list, inherits, logical(1), "spatial_field")))
  stopifnot(inherits(observation_data, "sf"))
  stopifnot(response %in% names(observation_data))

  lapply(field_list, function(field) {
    prepare_spatial_data(
      spatial_field = field,
      observation_data = observation_data,
      response = response,
      region_list = region_list,
      integration_rule = integration_rule,
      parallel = parallel,
      ...
    )
  })
}


#' Stack Spatial Design Matrices and Priors
#'
#' Takes a list of prepared spatial field objects (from `prepare_spatial_field_list`)
#' and stacks their design matrices and precision priors into a joint model matrix and block diagonal prior.
#'
#' @param spatial_data_list A list of `prepared_spatial_data` objects (one per covariate).
#' @param compute_cholesky A boolean; if TRUE, the cholesky factor of the precision is saved using Matrix::chol
#'
#' @return A list with components:
#'   - `y`: the shared response vector
#'   - `X`: stacked design matrix
#'   - `Lambda`: block-diagonal precision matrix
#'   - `fields`: original list of spatial_data objects
#'   - `column_map`: a named list mapping covariate labels to column indices in `X`
#'
#' @export
stack_spatial_design_matrices <- function(spatial_data_list, compute_cholesky = FALSE) {
  stopifnot(is.list(spatial_data_list))
  if (length(spatial_data_list) == 0) stop("spatial_data_list must not be empty.")

  # Check all elements
  for (d in spatial_data_list) {
    if (!inherits(d, "prepared_spatial_data")) {
      stop("Each element of spatial_data_list must be a prepared_spatial_data object.")
    }
  }

  # Validate y consistency
  y_lengths <- vapply(spatial_data_list, function(d) length(d$y), integer(1))
  if (length(unique(y_lengths)) != 1) {
    stop("All spatial fields must share the same response length (i.e., same regions).")
  }

  y <- spatial_data_list[[1]]$y
  X_blocks <- lapply(spatial_data_list, function(d) d$X_spatial)
  Lambda_blocks <- lapply(spatial_data_list, function(d) compute_precision(d$prior))

  # Track which columns belong to which variable
  column_map <- list()
  col_start <- 1
  for (d in spatial_data_list) {
    k <- ncol(d$X_spatial)
    label <- d$label %||% "unnamed"
    column_map[[label]] <- col_start:(col_start + k - 1)
    col_start <- col_start + k
  }

  # Stack columns and construct block diagonal precision
  X_stacked <- do.call(cbind, X_blocks)
  XtX <- Matrix::crossprod(X_stacked)
  Lambda_block_diag <- Matrix::bdiag(Lambda_blocks)

  posterior_precision <- XtX+Lambda_block_diag

  cholesky_factor <- NULL
  if (compute_cholesky) {
    message("Computing Cholesky of posterior precision...")
    cholesky_factor <- try(Matrix::chol(posterior_precision), silent = TRUE)
    if (inherits(cholesky_factor, "try-error")) {
      warning("Cholesky failed: posterior precision is not positive definite.")
      cholesky_factor <- NULL
    }
  }

  stacked = list(
    y = y,
    X = X_stacked,
    XtX = XtX,
    cholesky_factor = cholesky_factor,
    posterior_precision = posterior_precision,
    Lambda = Lambda_block_diag,
    fields = spatial_data_list,
    column_map = column_map
  )

  class(stacked) <- "stacked_spatial_data"

  return(stacked)
}


#' Orthogonalize a design block against others (sparse/scalable version)
#'
#' Projects the target block onto the orthogonal complement of the other blocks
#' using Gram-matrix projection. Preserves sparsity and scales to large n.
#'
#' @param stacked A list from `stack_spatial_design_matrices()`.
#' @param target_label A string for the field to orthogonalize (e.g., "Spatial field").
#' @return A modified `stacked` object with orthogonalized design block and adjusted prior.
#' @export
orthogonalize_design_block <- function(stacked, target_label) {
  message("Orthogonalizing design block against all others (sparse): ", target_label)
  stopifnot("X" %in% names(stacked), "column_map" %in% names(stacked))

  O <- orthogonalize_matrices(stacked$X,stacked$column_map,target_label)

  stacked$X <- O$X_new
  stacked$A <- O$A
  message("Computing new XtX...")
  stacked$XtX <- Matrix::crossprod(Matrix::Matrix(stacked$X))
  #stacked$Lambda <- Lambda_new
  return(stacked)
}


#' Fit a Spatial Field Model Without Covariates
#'
#' Estimates the posterior mean of the spatial random field coefficients, and optionally
#' returns the Cholesky factor or full posterior variance matrix.
#'
#' @param spatial_data Either a `prepared_spatial_data` object (from `prepare_spatial_data()`)
#'   or a `spatial_field` object.
#' @param observation_data If `spatial_data` is a `spatial_field`, the observation data (sf).
#' @param response If `spatial_data` is a `spatial_field`, the name of the response variable.
#' @param region_list Optional: passed to `prepare_spatial_data()`.
#' @param integration_rule Optional: override the integration rule.
#' @param parallel Logical or integer: whether/how to parallelize integration.
#' @param return_cholesky Logical: if `TRUE`, return the Cholesky factor of the posterior precision matrix.
#'   If `spatial_data` already includes a `cholesky_factor` entry, it will be reused.
#' @param return_variance Logical: if `TRUE`, return the full posterior variance matrix (inverse of posterior precision).
#'   If `return_cholesky` is also `TRUE`, the Cholesky is used to compute the inverse efficiently.
#' @param ... Additional arguments passed to `prepare_spatial_data()` if needed.
#'
#' @return A list of class `"fitted_spatial_field"` containing:
#' \describe{
#'   \item{posterior_mean}{Estimated coefficients of the spatial field.}
#'   \item{fitted}{Fitted values at observed locations.}
#'   \item{residuals}{Observed minus fitted values.}
#'   \item{X}{Design matrix used for estimation.}
#'   \item{spatial_data}{Original input or prepared spatial data object.}
#'   \item{cholesky_factor}{(Optional) Cholesky factor of the posterior precision matrix.}
#'   \item{posterior_variance}{(Optional) Posterior variance matrix of the coefficients.}
#' }
#'
#' @export
fit_spatial_field <- function(
    spatial_data,
    observation_data = NULL,
    response = NULL,
    region_list = NULL,
    integration_rule = NULL,
    parallel = FALSE,
    return_cholesky = FALSE,
    return_variance = FALSE,
    diagonal_only = FALSE,
    ...
) {
  if (inherits(spatial_data, "prepared_spatial_data")) {
    design <- spatial_data
  } else if (inherits(spatial_data, "spatial_field")) {
    message("Preparing spatial design (via integration or areas of intersections)...")
    design <- prepare_spatial_data(
      spatial_field = spatial_data,
      observation_data = observation_data,
      response = response,
      region_list = region_list,
      integration_rule = integration_rule,
      parallel = parallel,
      ...
    )
  } else {
    stop("`spatial_data` must be either a `prepared_spatial_data` or `spatial_field` object.")
  }

  y <- design$y
  X <- design$X_spatial
  Lambda <- compute_precision(design$prior)
  if (is.null(Lambda)) stop("The spatial_field does not have a prior specified.")

  message("Computing posterior precision...")
  XtX <- Matrix::crossprod(Matrix::Matrix(X))
  posterior_precision <- XtX + Lambda

  # Initialize optional outputs
  cholesky_factor <- NULL
  posterior_variance <- NULL

  # <<< UPDATED: Check for cached Cholesky in design object
  if (!is.null(design$cholesky_factor)) {
    message("Using cached Cholesky factor from design.")
    cholesky_factor <- design$cholesky_factor
  } else if (return_cholesky || return_variance) {
    message("Computing Cholesky factor of posterior precision...")
    cholesky_factor <- try(Matrix::chol(posterior_precision), silent = TRUE)

    if (inherits(cholesky_factor, "try-error")) {
      warning("Cholesky failed: posterior precision not positive definite.")
      cholesky_factor <- NULL
    }
  }

  message("Computing posterior mean...")
  rhs <- Matrix::crossprod(X, y)
  beta_hat <- as.numeric(solve_system(posterior_precision, rhs))

  message("Computing fitted values...")
  fitted <- as.numeric(X %*% beta_hat)
  residuals <- y - fitted

  out <- list(
    posterior_mean = beta_hat,
    fitted = fitted,
    residuals = residuals,
    X = X,
    XtX = XtX,
    Lambda = Lambda,
    spatial_data = spatial_data,
    cholesky_factor = cholesky_factor
  )

  if (return_variance) {
    out$posterior_variance <- compute_posterior_variance(
      precision_matrix = posterior_precision,
      cholesky_factor = cholesky_factor,
      diagonal_only = diagonal_only
    )

  out$sigma2 <- compute_posterior_measurement_error(
      posterior_mean = out$posterior_mean,
      X = X,
      Lambda = Lambda,
      observations = y
    )

  out$posterior_variance <- out$sigma2$sigma2_mode * out$posterior_variance

  }




  class(out) <- "fitted_spatial_field"
  return(out)
}



#' Fit a Stacked Spatial Model
#'
#' Estimates the posterior mean and fitted values from a joint spatial model
#' with one or more spatial fields, using a stacked design matrix and prior.
#'
#' @param stacked A list returned by `stack_spatial_design_matrices()`.
#' @param orthogonalize Optional. A string indicating which label (from `column_map`) to orthogonalize
#' @param return_variance Optional. A boolean; if TRUE returns the entire posterior covariance, otherwise returns only the diagonal.
#' @param diagonal_only Optional.
#'   against all other fields. If NULL, no orthogonalization is done.
#'
#' @return An object of class `"fitted_spatial_field"` with components:
#'   - `posterior_mean`: estimated coefficients
#'   - `fitted`: fitted values
#'   - `residuals`: residual vector
#'   - `X`: design matrix
#'   - `Lambda`: precision matrix
#'   - `y`: response vector
#'   - `column_map`: index map to design columns
#'   - `fields`: original field metadata
#'
#' @export
fit_stacked_spatial_field <- function(stacked, orthogonalize = NULL,return_variance = FALSE, diagonal_only = FALSE) {
  stopifnot(is.list(stacked))
  stopifnot(all(c("y", "X", "Lambda", "column_map", "fields") %in% names(stacked)))

  if (!is.null(orthogonalize)) {

    if(is.character(orthogonalize)){
      target_label = orthogonalize
      stacked <- orthogonalize_design_block(stacked, target_label =target_label)}
    if(orthogonalize == TRUE){
      target_label <- names(stacked$column_map)[1]
      message("Orthogonalizing against first spatial field in list...")
      stacked <- orthogonalize_design_block(stacked, target_label = target_label)

    }



  }

  y <- stacked$y
  X <- stacked$X
  Lambda <- stacked$Lambda

  message("Computing posterior precision...")
  posterior_precision <- stacked$XtX + Lambda

  message("Solving for posterior mean...")
  rhs <- Matrix::crossprod(X, y)
  beta_hat <- as.numeric(solve_system(posterior_precision, rhs))

  message("Computing fitted values...")
  fitted <- as.numeric(X %*% beta_hat)
  residuals <- y - fitted

  out <- list(
    posterior_mean = beta_hat,
    fitted = fitted,
    residuals = residuals,
    X = X,
    XtX = stacked$XtX,
    Lambda = Lambda,
    y = y,
    column_map = stacked$column_map,
    fields = stacked$fields,
    spatial_data = stacked
  )

  if (return_variance) {
    out$posterior_variance <- compute_posterior_variance(
      precision_matrix = posterior_precision,
      cholesky_factor = stacked$cholesky_factor,
      diagonal_only = diagonal_only)
     measurement_error <- compute_posterior_measurement_error(out$posterior_mean,out$X,out$Lambda,out$y)
     out$posterior_variance <- out$posterior_variance*measurement_error$sigma2_mode
     out$sigma2_a_n <- measurement_error$a_n
     out$sigma2_b_n <- measurement_error$b_n

  }

  if (!is.null(orthogonalize)) {
    out$orthogonalize <- orthogonalize
    if(!identical(orthogonalize, FALSE)){
      out$target_label <- target_label
    }
  }




  class(out) <- "fitted_spatial_field"
  return(out)
}



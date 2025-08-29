#' Cross-Validation for Tau in a Single Spatial Field
#'
#'
#' Performs K-fold cross-validation to select the regularization parameter \eqn{\tau}
#' (precision multiplier) for a single spatial field model.
#'
#' This function uses the precomputed design matrix and prior structure from
#' a `prepared_spatial_data` object, evaluating a sequence of \eqn{\tau} values
#' and computing mean squared prediction error on held-out folds.
#'
#' @param prepared_spatial_data A `prepared_spatial_data` object from [prepare_spatial_data()].
#' @param tau_grid A numeric vector of \eqn{\tau} values to evaluate (controls prior strength).
#' @param linear_solver Solver to use for posterior mean: one of `"cg"`, `"pcg"`, `"cholesky"`, or `"matrix_free_cg"`.
#' @param k Number of cross-validation folds.
#' @param mc.cores Number of parallel workers to use (via `future::multisession`).
#' @param show_progress Logical; if `TRUE`, shows a progress bar (requires `progressr`).
#'
#' @return A data frame with one row per `tau`, containing:
#' \describe{
#'   \item{tau}{The value of the regularization parameter \eqn{\tau}.}
#'   \item{mean_cv_error}{Mean squared error across the CV folds.}
#'   \item{se_cv_error}{Standard error of the mean CV error (SD / sqrt(k)).}
#' }
#'
#' @examples
#' \dontrun{
#' grid <- st_make_grid(st_bbox(c(xmin=0, ymin=0, xmax=1, ymax=1)), n = c(10,10))
#' basis <- make_geometry_basis(st_sf(geometry = grid))
#' prior <- make_sar_prior(basis)
#' spatial_field <- make_spatial_field(basis = basis, prior = prior)
#'
#' obs <- simulate_observations(grid_sf, sf::st_sf(geometry = grid), noise_sd = 0.05)
#' prepared <- prepare_spatial_data(spatial_field, obs, response = "y")
#'
#' tau_grid <- 10^seq(-4, 2, length.out = 20)
#' results <- cv_tau_spatial_field(prepared, tau_grid = tau_grid)
#' }
#'
#' @export
cv_tau_spatial_field <- function(prepared_spatial_data, tau_grid,
                                 linear_solver = c("cg", "pcg", "cholesky", "matrix_free_cg"),
                                 k = 5, mc.cores = parallel::detectCores(), show_progress = TRUE,
                                 folds = sample(rep(seq_len(k), length.out = length(prepared_spatial_data$y)))) {

  linear_solver <- match.arg(linear_solver)
  if (!requireNamespace("future.apply", quietly = TRUE)) stop("Please install 'future.apply'.")

  if (inherits(future::plan("list")[[1]], "sequential")) {
    future::plan(future::multisession, workers = mc.cores)
  }

  message("Using cached design matrix and response...")
  X <- prepared_spatial_data$X_spatial
  y <- prepared_spatial_data$y
  prior_template <- prepared_spatial_data$prior

  eval_fn <- function(tau, X, y, folds, linear_solver, prior_template, k) {
    fold_errors <- numeric(k)

    prior_tau <- prior_template
    prior_tau$hyperparameters$tau <- tau
    Q <- compute_precision(prior_tau)

    for (i in seq_len(k)) {
      test_idx <- which(folds == i)
      train_idx <- setdiff(seq_len(nrow(X)), test_idx)

      X_train <- X[train_idx, , drop = FALSE]
      y_train <- y[train_idx]
      X_test <- X[test_idx, , drop = FALSE]
      y_test <- y[test_idx]

      XtX_train <- Matrix::crossprod(X_train)
      Xty_train <- Matrix::crossprod(X_train, y_train)
      precision <- XtX_train + Q

      mu_hat <- estimate_target_grid(precision_matrix = precision,
                                     Atb = Xty_train,
                                     linear_solver = linear_solver)

      preds <- as.vector(X_test %*% mu_hat)
      fold_errors[i] <- mean((preds - y_test)^2)
    }

    data.frame(
      tau = tau,
      mean_cv_error = mean(fold_errors),
      se_cv_error = sd(fold_errors) / sqrt(k)
    )
  }

  if (show_progress && requireNamespace("progressr", quietly = TRUE)) {
    progressr::handlers("txtprogressbar")
    progressr::with_progress({
      p <- progressr::progressor(steps = length(tau_grid))
      results <- future.apply::future_lapply(
        tau_grid,
        function(tau) {
          p(sprintf("tau = %.4f", tau))
          eval_fn(tau, X, y, folds, linear_solver, prior_template, k)
        },
        future.globals = list(X = X, y = y, folds = folds,
                              linear_solver = linear_solver,
                              prior_template = prior_template,
                              eval_fn = eval_fn, k = k),
        future.packages = c("Matrix", "downscaling", "SpatialBasis"),
        future.seed = TRUE
      )
      do.call(rbind, results)
    })
  } else {
    results <- future.apply::future_lapply(
      tau_grid,
      function(tau) {
        eval_fn(tau, X, y, folds, linear_solver, prior_template, k)
      },
      future.globals = list(X = X, y = y, folds = folds,
                            linear_solver = linear_solver,
                            prior_template = prior_template,
                            eval_fn = eval_fn, k = k),
      future.packages = c("Matrix", "downscaling", "SpatialBasis"),
      future.seed = TRUE
    )
    do.call(rbind, results)
  }
}

#' Cross-Validation for Tau in a Stacked Spatial Field
#'
#' Performs K-fold cross-validation to select the regularization parameter \eqn{\tau}
#' (precision multiplier) for a single spatial field model.
#'
#' This function uses the precomputed design matrix and prior structure from
#' a `prepared_spatial_data` object, evaluating a sequence of \eqn{\tau} values
#' and computing mean squared prediction error on held-out folds.
#'
#' @param stacked_spatial_data A `stacked_spatial_data` object from [stack_spatial_design_matrices()].
#' @param tau_grid A numeric vector of \eqn{\tau} values to evaluate (controls prior strength).
#' @param linear_solver Solver to use for posterior mean: one of `"cg"`, `"pcg"`, `"cholesky"`, or `"matrix_free_cg"`.
#' @param k Number of cross-validation folds.
#' @param mc.cores Number of parallel workers to use (via `future::multisession`).
#' @param show_progress Logical; if `TRUE`, shows a progress bar (requires `progressr`).
#'
#' @return A data frame with one row per `tau`, containing:
#' \describe{
#'   \item{tau}{The value of the regularization parameter \eqn{\tau}.}
#'   \item{mean_cv_error}{Mean squared error across the CV folds.}
#'   \item{se_cv_error}{Standard error of the mean CV error (SD / sqrt(k)).}
#' }
#'
#' @examples
#' \dontrun{
#' grid <- st_make_grid(st_bbox(c(xmin=0, ymin=0, xmax=1, ymax=1)), n = c(10,10))
#' basis <- make_geometry_basis(st_sf(geometry = grid))
#' prior <- make_sar_prior(basis)
#' spatial_field <- make_spatial_field(basis = basis, prior = prior)
#'
#' obs <- simulate_observations(grid_sf, sf::st_sf(geometry = grid), noise_sd = 0.05)
#' prepared <- prepare_spatial_data(spatial_field, obs, response = "y")
#'
#' tau_grid <- 10^seq(-4, 2, length.out = 20)
#' results <- cv_tau_spatial_field(prepared, tau_grid = tau_grid)
#' }
#'
#' @export
cv_tau_stacked_spatial_field <- function(prepared_spatial_field_list, tau_grid,
                                 linear_solver = c("cg", "pcg", "cholesky", "matrix_free_cg"),
                                 k = 5, mc.cores = parallel::detectCores(), show_progress = TRUE,
                                 folds = sample(rep(seq_len(k), length.out = length(prepared_spatial_field_list[[1]]$y)))) {

  linear_solver <- match.arg(linear_solver)
  if (!requireNamespace("future.apply", quietly = TRUE)) stop("Please install 'future.apply'.")

  if (inherits(future::plan("list")[[1]], "sequential")) {
    future::plan(future::multisession, workers = mc.cores)
  }

  message("Using cached design matrix and response...")

  stacked <- stack_spatial_design_matrices(prepared_spatial_field_list)

  X <- stacked$X
  y <- stacked$y
  Lambda <- stacked$Lambda
  column_map <- stacked$column_map[[1]]

  prior_template <- prepared_spatial_field_list[[1]]$prior

  eval_fn <- function(tau, X, y, Lambda, column_map, folds, linear_solver, prior_template, k) {
    fold_errors <- numeric(k)

    prior_tau <- prior_template
    prior_tau$hyperparameters$tau <- tau

    Q <- Lambda
    Q[column_map,column_map] <- compute_precision(prior_tau)

    for (i in seq_len(k)) {
      test_idx <- which(folds == i)
      train_idx <- setdiff(seq_len(nrow(X)), test_idx)

      X_train <- X[train_idx, , drop = FALSE]
      y_train <- y[train_idx]
      X_test <- X[test_idx, , drop = FALSE]
      y_test <- y[test_idx]

      XtX_train <- Matrix::crossprod(X_train)
      Xty_train <- Matrix::crossprod(X_train, y_train)
      precision <- XtX_train + Q

      mu_hat <- solve_system(precision_matrix = precision,
                                     Atb = Xty_train,
                                     linear_solver = linear_solver)

      preds <- as.vector(X_test %*% mu_hat)
      fold_errors[i] <- mean((preds - y_test)^2)
    }

    data.frame(
      tau = tau,
      mean_cv_error = mean(fold_errors),
      se_cv_error = sd(fold_errors) / sqrt(k)
    )
  }

  if (show_progress && requireNamespace("progressr", quietly = TRUE)) {
    progressr::handlers("txtprogressbar")
    progressr::with_progress({
      p <- progressr::progressor(steps = length(tau_grid))
      results <- future.apply::future_lapply(
        tau_grid,
        function(tau) {
          p(sprintf("tau = %.4f", tau))
          eval_fn(tau, X, y, Lambda, column_map, folds, linear_solver, prior_template, k)
        },
        future.globals = list(X = X, y = y, Lambda = Lambda, column_map = column_map, folds = folds,
                              linear_solver = linear_solver,
                              prior_template = prior_template,
                              eval_fn = eval_fn, k = k),
        future.packages = c("Matrix", "SpatialBasis"),
        future.seed = TRUE
      )
      do.call(rbind, results)
    })
  } else {
    results <- future.apply::future_lapply(
      tau_grid,
      function(tau) {
        eval_fn(tau, X, y, Lambda, column_map, folds, linear_solver, prior_template, k)
      },
      future.globals = list(X = X, y = y, Lambda = Lambda, column_map = column_map, folds = folds,
                            linear_solver = linear_solver,
                            prior_template = prior_template,
                            eval_fn = eval_fn, k = k),
      future.packages = c("Matrix", "SpatialBasis"),
      future.seed = TRUE
    )
    do.call(rbind, results)
  }
}


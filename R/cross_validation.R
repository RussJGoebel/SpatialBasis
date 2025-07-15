#' Parallel Grid Search Cross-Validation (per hyperparameter x fold)
#'
#' @param prepared_data A prepared_spatial_data object.
#' @param prior_template A prior object (e.g., from make_sar_prior()).
#' @param hyperparameter_grid A data frame of hyperparameter combinations.
#' @param folds Number of folds (integer) or list of test indices.
#' @param metric Function for error metric.
#' @param n_workers Number of parallel workers (default: all available cores).
#' @param use_progress Whether to show progress bar.
#' @param verbose Whether to print messages.
#'
#' @return Data frame of hyperparameters and mean errors.
#' @export
grid_search_cv_parallel <- function(
    prepared_data,
    prior_template,
    hyperparameter_grid,
    folds = 5,
    metric = function(truth, pred) sqrt(mean((truth - pred)^2)),
    n_workers = future::availableCores(),
    use_progress = TRUE,
    verbose = TRUE
) {
  stopifnot(inherits(prepared_data, "prepared_spatial_data"))
  stopifnot(is.data.frame(hyperparameter_grid))
  stopifnot(inherits(prior_template, "prior"))

  y <- prepared_data$y
  X_spatial <- prepared_data$X_spatial
  n <- length(y)

  # Detect X_spatial type and prepare storage
  if (inherits(X_spatial, "dgCMatrix")) {
    X_format <- "sparse"
    X_slots <- list(
      i = X_spatial@i,
      p = X_spatial@p,
      x = X_spatial@x,
      Dim = X_spatial@Dim
    )
  } else if (inherits(X_spatial, "dgeMatrix") || is.matrix(X_spatial)) {
    X_format <- "dense"
    X_dense <- as.matrix(X_spatial)
  } else {
    stop("X_spatial must be either a dgCMatrix, dgeMatrix, or base R matrix.")
  }

  # Make folds
  if (is.numeric(folds) && length(folds) == 1) {
    set.seed(123)
    fold_indices <- split(sample(seq_len(n)), rep(1:folds, length.out = n))
  } else if (is.list(folds)) {
    fold_indices <- folds
  } else {
    stop("folds must be integer or list of indices.")
  }

  # Create all combinations (hyperparameter x fold)
  param_fold_grid <- expand.grid(
    param_i = seq_len(nrow(hyperparameter_grid)),
    fold_k = seq_along(fold_indices)
  )

  total_tasks <- nrow(param_fold_grid)

  # Setup future plan
  if (inherits(future::plan("next"), "sequential")) {
    future::plan(future::multisession, workers = n_workers)
  }

  # Progress handling
  if (use_progress) {
    progressr::handlers("txtprogressbar")
    progressr::with_progress({
      p <- progressr::progressor(steps = total_tasks)

      results <- future.apply::future_lapply(
        seq_len(total_tasks),
        function(idx) {
          row <- param_fold_grid[idx, ]
          i <- row$param_i
          k <- row$fold_k

          # Reconstruct X inside the worker
          if (X_format == "sparse") {
            X <- new("dgCMatrix",
                     i = X_slots$i,
                     p = X_slots$p,
                     Dim = X_slots$Dim,
                     x = X_slots$x,
                     Dimnames = list(NULL, NULL),
                     factors = list())
          } else {
            X <- X_dense
          }

          hyperparams <- hyperparameter_grid[i, , drop = FALSE]

          if (verbose) {
            message("Evaluating grid row ", i, " fold ", k, "...")
          }

          prior <- prior_template
          prior$hyperparameters <- as.list(hyperparams)

          test_idx <- fold_indices[[k]]
          train_idx <- setdiff(seq_len(n), test_idx)

          y_train <- y[train_idx]
          y_test <- y[test_idx]
          X_train <- X[train_idx, , drop = FALSE]
          X_test <- X[test_idx, , drop = FALSE]

          Q <- compute_precision(prior)

          XtX <- Matrix::crossprod(X_train)
          rhs <- Matrix::crossprod(X_train, y_train)
          posterior_precision <- XtX + Q

          beta_hat <- as.numeric(solve_system(posterior_precision, rhs))
          y_hat <- as.numeric(X_test %*% beta_hat)

          fold_error <- metric(y_test, y_hat)

          p(sprintf("Row %d fold %d complete", i, k))

          # Remove big objects and force garbage collection
          rm(X, X_train, X_test, XtX, rhs, posterior_precision, Q, beta_hat, y_hat)
          gc()

          data.frame(
            param_i = i,
            fold_k = k,
            error = fold_error
          )
        },
        future.seed = TRUE,
        future.packages = "Matrix"
      )
    })
  } else {
    results <- future.apply::future_lapply(
      seq_len(total_tasks),
      function(idx) {
        row <- param_fold_grid[idx, ]
        i <- row$param_i
        k <- row$fold_k

        # Reconstruct X inside the worker
        if (X_format == "sparse") {
          X <- new("dgCMatrix",
                   i = X_slots$i,
                   p = X_slots$p,
                   Dim = X_slots$Dim,
                   x = X_slots$x,
                   Dimnames = list(NULL, NULL),
                   factors = list())
        } else {
          X <- X_dense
        }

        hyperparams <- hyperparameter_grid[i, , drop = FALSE]

        if (verbose) {
          message("Evaluating grid row ", i, " fold ", k, "...")
        }

        prior <- prior_template
        prior$hyperparameters <- as.list(hyperparams)

        test_idx <- fold_indices[[k]]
        train_idx <- setdiff(seq_len(n), test_idx)

        y_train <- y[train_idx]
        y_test <- y[test_idx]
        X_train <- X[train_idx, , drop = FALSE]
        X_test <- X[test_idx, , drop = FALSE]

        Q <- compute_precision(prior)

        XtX <- Matrix::crossprod(X_train)
        rhs <- Matrix::crossprod(X_train, y_train)
        posterior_precision <- XtX + Q

        beta_hat <- as.numeric(solve_system(posterior_precision, rhs))
        y_hat <- as.numeric(X_test %*% beta_hat)

        fold_error <- metric(y_test, y_hat)

         # Remove big objects and force garbage collection
         rm(X, X_train, X_test, XtX, rhs, posterior_precision, Q, beta_hat, y_hat)
         gc()

        data.frame(
          param_i = i,
          fold_k = k,
          error = fold_error
        )
      },
      future.seed = TRUE,
      future.packages = "Matrix"
    )
  }

  # Combine all results
  results_df <- do.call(rbind, results)

  # Aggregate mean error per hyperparameter combination
  summary_df <- aggregate(
    error ~ param_i,
    data = results_df,
    FUN = mean
  )

  # Add hyperparameters back
  hyperparams_df <- hyperparameter_grid
  hyperparams_df$param_i <- seq_len(nrow(hyperparams_df))
  final_df <- merge(hyperparams_df, summary_df, by = "param_i")
  final_df <- final_df[order(final_df$error), ]
  rownames(final_df) <- NULL
  names(final_df)[names(final_df) == "error"] <- "mean_error"

  if (verbose) {
    message("Grid search complete. Best config:")
    print(final_df[1, ])
  }

  return(final_df)
}

#' Combine per-field hyperparameter grids into a joint (row-wise) grid
#'
#' @param grids A named list of data.frames, one per field. Each data.frame must have 1 or N rows,
#'              where N is the maximum across all fields.
#'
#' @return A joint data.frame with column names like `param_field`.
#' @export
combine_stacked_hyperparameters <- function(grids) {
  stopifnot(is.list(grids))
  field_names <- names(grids)
  if (is.null(field_names)) stop("List must be named by field.")

  nrows_each <- vapply(grids, nrow, integer(1))
  max_rows <- max(nrows_each)

  # Expand/repeat 1-row entries or error
  grids_expanded <- mapply(function(df, field) {
    if (!is.data.frame(df)) stop(sprintf("Grid for '%s' is not a data.frame", field))
    if (nrow(df) == 0) return(NULL)  # treat as fixed (no params)
    if (nrow(df) == 1 && max_rows > 1) {
      df <- df[rep(1, max_rows), , drop = FALSE]
    } else if (nrow(df) != max_rows) {
      stop(sprintf("Grid for '%s' must have 1 or %d rows (has %d)", field, max_rows, nrow(df)))
    }

    # Rename columns to param_field
   # names(df) <- paste0(names(df), "_", field)
    df
  }, grids, field_names, SIMPLIFY = FALSE)

  result <- do.call(cbind, Filter(Negate(is.null), grids_expanded))
  rownames(result) <- NULL
  return(result)
}


#' Parallel Grid Search CV for Stacked Spatial Models (joint grid)
#'
#' @param prepared_list A list of prepared_spatial_data objects.
#' @param hyperparameter_grid A data.frame of joint hyperparameter combinations. Each column should be named as `fieldname.param`, e.g., `water.rho`, `land.tau`.
#' @param folds Number of folds (integer) or list of test indices.
#' @param metric Function for error metric.
#' @param n_workers Number of parallel workers.
#' @param use_progress Show progress bar.
#' @param verbose Print messages.
#' @param orthogonalize Optional. A string indicating which field to orthogonalize.
#'
#' @return Data frame of grid results (mean error per combination).
#' @export
grid_search_cv_stacked <- function(
    prepared_list,
    hyperparameter_grid,
    folds = 5,
    metric = function(truth, pred) sqrt(mean((truth - pred)^2)),
    n_workers = future::availableCores(),
    use_progress = TRUE,
    verbose = TRUE,
    orthogonalize = NULL
) {
  stopifnot(is.list(prepared_list))
  stopifnot(is.data.frame(hyperparameter_grid))

  field_labels <- vapply(prepared_list, function(p) p$label, character(1))

  y <- prepared_list[[1]]$y
  n <- length(y)

  if (is.numeric(folds)) {
    set.seed(123)
    fold_indices <- split(sample(seq_len(n)), rep(1:folds, length.out = n))
  } else {
    fold_indices <- folds
  }

  total_tasks <- nrow(hyperparameter_grid) * length(fold_indices)

  if (inherits(future::plan("next"), "sequential")) {
    future::plan(future::multisession, workers = n_workers)
  }

  if (use_progress) {
    progressr::handlers("txtprogressbar")
    progressr::with_progress({
      p <- progressr::progressor(steps = total_tasks)
      results <- future.apply::future_lapply(seq_len(total_tasks), function(idx) {
        param_i <- (idx - 1) %% nrow(hyperparameter_grid) + 1
        fold_k <- ((idx - 1) %/% nrow(hyperparameter_grid)) + 1

        test_idx <- fold_indices[[fold_k]]
        train_idx <- setdiff(seq_len(n), test_idx)

        for (j in seq_along(prepared_list)) {
          label <- prepared_list[[j]]$label

          matching_cols <- grep(paste0("^", label, "\\."), names(hyperparameter_grid), value = TRUE)

          if (length(matching_cols) > 0) {
            param_values <- hyperparameter_grid[param_i, matching_cols, drop = FALSE]
            names(param_values) <- sub(paste0("^", label, "\\."), "", matching_cols)

            prepared_list[[j]]$prior$hyperparameters <- as.list(param_values)
          }
        }

        stacked <- stack_spatial_design_matrices(prepared_list)

        if (!is.null(orthogonalize)) {
          stacked <- orthogonalize_design_block(stacked, target_label = orthogonalize)
        }

        X_train <- stacked$X[train_idx, , drop = FALSE]
        X_test  <- stacked$X[test_idx, , drop = FALSE]
        y_train <- stacked$y[train_idx]
        y_test  <- stacked$y[test_idx]
        Lambda  <- stacked$Lambda

        XtX <- Matrix::crossprod(X_train)
        rhs <- Matrix::crossprod(X_train, y_train)
        posterior_precision <- XtX + Lambda
        beta_hat <- as.numeric(solve_system(posterior_precision, rhs))
        y_hat <- as.numeric(X_test %*% beta_hat)

        fold_error <- metric(y_test, y_hat)
        p(sprintf("Row %d Fold %d done", param_i, fold_k))

        result <- data.frame(param_i = param_i, fold = fold_k, error = fold_error)
        result <- cbind(result, hyperparameter_grid[param_i, , drop = FALSE])
        result
      }, future.seed = TRUE)
    })
  } else {
    results <- future.apply::future_lapply(...) # fallback
  }

  results <- Filter(function(x) !is.null(x) && "error" %in% names(x), results)
  if (length(results) == 0) stop("All CV tasks failed. Check for upstream errors.")

  results_df <- do.call(rbind, results)
  summary_df <- aggregate(error ~ param_i, data = results_df, FUN = mean)
  final_grid <- merge(hyperparameter_grid, summary_df, by.x = 0, by.y = "param_i")
  final_grid <- final_grid[order(final_grid$error), ]
  rownames(final_grid) <- NULL
  names(final_grid)[names(final_grid) == "error"] <- "mean_error"

  return(final_grid)
}

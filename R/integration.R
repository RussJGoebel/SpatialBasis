#' Integrate a List of Regions Using a Spatial Basis
#'
#' Applies a specified integration rule and basis to each region in a list of
#' `region_atomic` or `region_composite` objects. Handles `NULL` entries by returning
#' `NA` vectors of the appropriate dimension.
#'
#' @param region_list A list of region objects as returned by [build_region_list()].
#' @param rule An integration rule object (e.g., from [integration_rule_qmc_average()]).
#' @param basis A basis object (e.g., from [make_geometry_basis()] or [make_function_basis()]).
#'
#' @return A matrix of dimension `length(region_list) × basis$k`, where each row
#'         contains the integrated basis values for one region.
#'
#' @export
integrate_region_list <- function(region_list, rule, basis) {
  stopifnot(inherits(basis, "basis"), !is.null(basis$k))

  results <- mapply(
    function(region, idx) {
      if (is.null(region)) {
        return(rep(NA_real_, basis$k))
      }
      evaluate_integral(rule, region, basis, index = idx)
    },
    region_list,
    seq_along(region_list),
    SIMPLIFY = FALSE
  )

  do.call(rbind, results)
}

#' Integrate a List of Regions Using a Spatial Basis (Parallel with Progress)
#'
#' Applies a specified integration rule and basis to each region in a list of
#' `region_atomic` or `region_composite` objects using parallel computation
#' via `future.apply::future_mapply()`, with optional progress reporting.
#'
#' @param region_list A list of region objects as returned by [build_region_list()].
#' @param rule An integration rule object (e.g., from [integration_rule_qmc_average()]).
#' @param basis A basis object (e.g., from [make_geometry_basis()] or [make_function_basis()]).
#' @param use_progress Logical; if TRUE, shows a progress bar using [furrr].
#' @param n_workers Number of parallel workers (default: all available cores).
#'
#' @return A matrix of dimension `length(region_list) × basis$k`.
#'
#' @export
integrate_region_list_parallel <- function(
    region_list,
    rule,
    basis,
    use_progress = TRUE,
    n_workers = future::availableCores()
) {
  stopifnot(inherits(basis, "basis"), !is.null(basis$k))
  stopifnot(inherits(rule, "integration_rule"))

  # Plan workers if sequential
  if (inherits(future::plan("list")[[1]], "sequential")) {
    future::plan(future::multisession, workers = n_workers)
  }

  indices <- seq_along(region_list)

  if (use_progress) {
    progressr::handlers("txtprogressbar")  # CLI bar
    progressr::with_progress({
      p <- progressr::progressor(steps = length(region_list))

      res <- future.apply::future_mapply(
        function(region, idx) {
          p()
          if (is.null(region)) {
            return(rep(NA_real_, basis$k))
          }
          downscaling::evaluate_integral(rule, region, basis, index = idx)
        },
        region_list,
        indices,
        SIMPLIFY = FALSE,
        future.globals = list(rule = rule, basis = basis),
        future.packages = c("downscaling", "sf"),
        future.seed = TRUE
      )

      do.call(rbind, res)
    })
  } else {
    res <- future.apply::future_mapply(
      function(region, idx) {
        if (is.null(region)) {
          return(rep(NA_real_, basis$k))
        }
        downscaling::evaluate_integral(rule, region, basis, index = idx)
      },
      region_list,
      indices,
      SIMPLIFY = FALSE,
      future.globals = list(rule = rule, basis = basis),
      future.packages = c("downscaling", "sf"),
      future.seed = TRUE
    )

    do.call(rbind, res)
  }
}


#' Batched Parallel Integrate Region List
#'
#' Processes regions in batches to reduce per-task overhead
#' while retaining dynamic load balancing.
#'
#' @param region_list List of regions to integrate.
#' @param rule Integration rule object.
#' @param basis Basis object.
#' @param n_batches Number of batches (default: 3000).
#' @param n_workers Number of workers (default: all available cores).
#' @param use_progress Show progress bar (default: TRUE).
#'
#' @return Matrix of results (length(region_list) × basis$k).
#'
#' @export
integrate_region_list_batched_parallel <- function(
    region_list,
    rule,
    basis,
    n_batches = 3000,
    n_workers = future::availableCores(),
    use_progress = TRUE
) {
  stopifnot(inherits(basis, "basis"), !is.null(basis$k))
  stopifnot(inherits(rule, "integration_rule"))

  if (identical(future::plan("list")[[1]], future::sequential)) {
    future::plan(future::multisession, workers = n_workers)
  }

  # Split regions and indices into batches
  indices <- seq_along(region_list)
  batch_size <- ceiling(length(region_list) / n_batches)
  batch_ids <- ceiling(indices / batch_size)
  batches <- split(seq_along(region_list), batch_ids)

  message("Processing ", length(batches), " batches...")

  # Main logic
  if (use_progress) {
    progressr::handlers("txtprogressbar")
    progressr::with_progress({
      p <- progressr::progressor(steps = length(region_list))

      res_batches <- future.apply::future_lapply(
        batches,
        function(batch_idx) {
          lapply(batch_idx, function(idx) {
            p()
            region <- region_list[[idx]]
            if (is.null(region)) {
              return(rep(NA_real_, basis$k))
            }
            downscaling::evaluate_integral(rule, region, basis, index = idx)
          })
        },
        future.seed = TRUE
      )
    })
  } else {
    res_batches <- future.apply::future_lapply(
      batches,
      function(batch_idx) {
        lapply(batch_idx, function(idx) {
          region <- region_list[[idx]]
          if (is.null(region)) {
            return(rep(NA_real_, basis$k))
          }
          downscaling::evaluate_integral(rule, region, basis, index = idx)
        })
      },
      future.seed = TRUE
    )
  }

  # Flatten nested list of lists
  res_flat <- unlist(res_batches, recursive = FALSE)

  # Row-bind into matrix
  res_mat <- do.call(rbind, res_flat)

  rm(res_batches)
  gc()

  return(res_mat)
}

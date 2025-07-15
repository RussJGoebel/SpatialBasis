#' st_fast_intersection (stack overflow)
#'
#' @param x An sf object, first parameter in sf:st_intersects
#' @param y An sf object, second parameter in sf:st_intersects
#' @param ... other arguments for sf::st_intersects
#'
#' @return Geometry if the the intersection between x and subset of y
#' @export
#'
#' @examples
st_fast_intersection <- function(x, y, ...) {
  y_idx <- sf::st_intersects(x, y) %>%
    unlist() %>%
    unique() %>%
    sort()

  if (length(y_idx) == 0) {
    # Return empty sf with same schema as y
    return(y[0, , drop = FALSE])
  }

  y_subset <- y[y_idx, , drop = FALSE]

  # Convert x to sf if needed
  if (inherits(x, "sfc")) {
    x <- sf::st_sf(geometry = x)
  }

  sf::st_intersection(x, y_subset, ...)
}

#' st_parallel_intersection
#'
#' Similar to st_intersection, but runs in parallel to make use of more processors.
#' Uses the future pacakge (for compatibility with windows systems)
#'
#'
#' @param soundings An sf object, The observed soundings
#' @param fine_grid An sf object, a fine grid on which we want to estimate SIF
#' @param mc.cores number of cores to using in mclapply
#' @param n_groups number core groups. defaults to mc.cores parameter
#'
#' @return sf dataframe of intersections.
#' @export
#'
#' @examples if(interactive()) intersections <- st_parallel_intersection(soundings, fine_grid,mc.cores = 2)
st_parallel_intersection <- function(soundings, fine_grid,
                                     mc.cores = parallel::detectCores(logical = FALSE),n_groups = mc.cores){

  n <- dim(fine_grid)[1]
  raster_pixels <- 1:n

  intersects <- sf::st_intersects(fine_grid,soundings)

  intersections <- parallel::mclapply(raster_pixels,function(x){sf::st_intersection(fine_grid[x,],
                                                                                    soundings[intersects[[x]],])},
                                      mc.cores = mc.cores)

  intersections <- dplyr::bind_rows(intersections)

  return(intersections)

}

#' st_parallel_intersection
#'
#' Parallel version of st_intersection: splits the fine grid into rows and intersects
#' with soundings in parallel.
#'
#' @param soundings An sf object: observed polygons.
#' @param fine_grid An sf object: grid polygons.
#' @param n_workers Integer: number of parallel workers (default: NULL to use the current plan).
#' @param use_progress Logical: show progress bar (default TRUE).
#'
#' @return sf dataframe of intersections.
#' @export
#'
#' @examples
#' \dontrun{
#' library(future)
#' plan(multisession)
#' intersections <- st_parallel_intersection(soundings, fine_grid, n_workers = 4)
#' }
st_parallel_intersection_future <- function(
    soundings,
    fine_grid,
    n_workers = NULL,
    use_progress = TRUE
) {
  n <- nrow(fine_grid)
  grid_indices <- seq_len(n)

  # Precompute which soundings intersect each grid cell
  intersects <- sf::st_intersects(fine_grid, soundings)

  # Temporarily set plan if specified
  if (!is.null(n_workers)) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = n_workers)
  }

  # Progress bar
  if (use_progress) {
    progressr::handlers("txtprogressbar")
    progressr::with_progress({
      p <- progressr::progressor(steps = n)
      intersection_list <- future.apply::future_lapply(
        grid_indices,
        function(i) {
          p()
          suppressWarnings(
            sf::st_intersection(
              fine_grid[i, , drop = FALSE],
              soundings[intersects[[i]], , drop = FALSE]
            ))
        }
      )
    })
  } else {
    intersection_list <- future.apply::future_lapply(
      grid_indices,
      function(i) {
        suppressWarnings(
          sf::st_intersection(
            fine_grid[i, , drop = FALSE],
            soundings[intersects[[i]], , drop = FALSE]
          ))
      }
    )
  }

  # Combine
  intersection_sf <- dplyr::bind_rows(intersection_list)
  return(intersection_sf)
}



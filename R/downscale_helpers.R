#' Parse downscale formula
#'
#' @param formula A model formula
#' @return A list:
#'   - fixed_terms: character vector of labels
#'   - random_bars: list of (1 | group) calls
#'   - spatial_terms: list of spatial_field objects
#' @export
parse_downscale_formula <- function(formula) {
  stopifnot(inherits(formula, "formula"))

  # Get random effects
  random_bars <- lme4::findbars(formula)

  # Remove random effects
  fixed_formula <- lme4::nobars(formula)

  # Get term labels
  tf <- terms(fixed_formula)
  term_labels <- attr(tf, "term.labels")
  intercept <- attr(tf, "intercept")

  fixed_terms <- character()
  spatial_terms <- list()

  formula_env <- environment(formula)

  for (label in term_labels) {
    obj <- tryCatch(
      eval(parse(text = label), envir = formula_env),
      error = function(e) NULL
    )

    if (inherits(obj, "spatial_field")) {
      spatial_terms[[length(spatial_terms) + 1]] <- obj
    } else {
      fixed_terms <- c(fixed_terms, label)
    }
  }

  if (intercept == 1 && length(fixed_terms) == 0) {
    fixed_terms <- "1"
  }

  list(
    fixed_terms = fixed_terms,
    random_bars = random_bars,
    spatial_terms = spatial_terms
  )
}

#' Validate downscaling inputs
#'
#' Checks that observation and grid data are well-formed.
#'
#' @param observation_data An sf object with the observed responses.
#' @param covariate_data A named list of sf objects with grids and covariates,
#'   or NULL if no grid covariates are provided.
#'
#' @return Invisibly TRUE if all checks pass; otherwise throws an error.
#' @export
validate_downscale_inputs <- function(observation_data, covariate_data) {
  # Check observation_data
  if (!inherits(observation_data, "sf")) {
    stop("`observation_data` must be an sf object.")
  }

  # If covariate_data is NULL, skip grid checks
  if (is.null(covariate_data)) {
    return(invisible(TRUE))
  }

  # Check covariate_data
  if (!is.list(covariate_data)) {
    stop("`covariate_data` must be a list of sf objects, or NULL.")
  }
  if (length(covariate_data) == 0) {
    # This is allowed if no grid covariates are used
    return(invisible(TRUE))
  }
  if (is.null(names(covariate_data)) || any(names(covariate_data) == "")) {
    stop("`covariate_data` must be a *named* list. Each grid needs a name.")
  }

  # Check that each entry is an sf object
  non_sf <- vapply(covariate_data, function(x) !inherits(x, "sf"), logical(1))
  if (any(non_sf)) {
    bad_names <- names(covariate_data)[non_sf]
    stop(sprintf("The following entries in `covariate_data` are not sf objects: %s",
                 paste(bad_names, collapse = ", ")))
  }

  # Check for duplicate names
  if (any(duplicated(names(covariate_data)))) {
    stop("`covariate_data` has duplicated names. Each grid must have a unique name.")
  }

  # Check CRS compatibility
  obs_crs <- sf::st_crs(observation_data)
  for (nm in names(covariate_data)) {
    grid_crs <- sf::st_crs(covariate_data[[nm]])
    if (obs_crs != grid_crs) {
      stop(sprintf("CRS mismatch: observation_data CRS (%s) does not match grid '%s' CRS (%s).",
                   obs_crs$input %||% "<unknown>",
                   nm,
                   grid_crs$input %||% "<unknown>"))
    }
  }

  invisible(TRUE)
}

#' Downscale coarse observations to a fine grid
#'
#' @param formula A model formula, e.g., y ~ 1 + spatial_field() + covariate.
#' @param data An sf object with observed polygons.
#' @param target_grid An sf object defining the grid to predict on.
#' @param region_method One of "raw", "triangulate", "intersect-then-triangulate".
#' @param strata Optional stratification polygons.
#' @param region_list Optional precomputed list of regions.
#' @param sampling "posterior_mode", "gibbs", or "elliptical_slice".
#' @param priors List of optional hyperparameters.
#' @param sigma2 Optional known observation variance.
#' @param ... Reserved for future.
#'
#' @return An object of class "downscale_fit".
fit_downscale <- function(
    formula,
    observation_data,
    covariate_data= NULL,
    target_grid = NULL,
    regionalization_method = c("raw", "triangulate", "intersect-then-triangulate"),
    integration_strata = NULL,
    region_list = NULL,
    sampling = c("posterior_mode", "gibbs", "elliptical_slice"),
    priors = list(),
    sigma2 = NULL,
    prepared_spatial_data = NULL,
    ...
) {

  regionalization_method <- match.arg(regionalization_method)

  if(is.null(integration_strata) & regionalization_method != "intersect-then-triangulate"){
    message("When regionalization_method is not `intersect-then-triangulate`, integration_strata are not used in integration.")
  }

  # 1. Parse formula
  parsed <- parse_downscale_formula(formula)
  # parsed$fixed_terms, parsed$spatial_terms

  # 2. Build or validate regions
  message("Building region list for integration...")

  if (!is.null(region_list)) {
    regions <- region_list
  } else {
    regions <- build_region_list(
      obs_polygons = observation_data,
      stratum_polygons = integration_strata,
      method = regionalization_method
    )
  }

  if(is.null(prepared_spatial_field)){
  message("Integrating regions...")
  prepared_spatial_field <- prepare_spatial_field(spatial_field = parsed$spatial_terms[[1]],
                                                  region_list = regions,
                                                  parallel = TRUE,n_workers = 16)
  }

  # 3. Fixed-effect design matrix
  #if (length(parsed$fixed_terms) > 0) {
  #  X <- model.matrix(
   #   reformulate(parsed$fixed_terms),
  #   #  data = data
  #   )
  # } else {
  #   X <- matrix(1, nrow=nrow(data), ncol=1)
  #   colnames(X) <- "(Intercept)"
  # }

  # 4. Random-effect design matrices and priors
  # B_list <- list()
  # Q_list <- list()
  #
  # for (i in seq_along(parsed$spatial_terms)) {
  #   field <- parsed$spatial_terms[[i]]
  #   B_i <- integrate_region_list(
  #     region_list = regions,
  #     rule = field$integration_rule,
  #     basis = field$basis
  #   )
  #   B_list[[i]] <- B_i
  #   Q_list[[i]] <- field$prior$precision
  # }
  #
  # # 5. Stack random effects
  # if (length(B_list) > 0) {
  #   B <- do.call(cbind, B_list)
  #   Q_block <- Matrix::bdiag(Q_list)
  # } else {
  #   B <- NULL
  #   Q_block <- NULL
  # }

  # 6. Prepare response
  lhs <- formula[[2]]
  y <- eval(lhs, envir = observation_data, enclos = environment(formula))
  n <- length(y)


  # 8. Package output
  result <- list(
    formula = formula,
    sampling = sampling,
    # X = X,
    # B = B,
    # Q = Q_block,
    # beta_hat = beta_hat,
    # w_hat = w_hat,
    # sigma2_hat = sigma2_hat,
    # samples = samples,
    region_list = regions,
    spatial_fields = parsed$spatial_terms,
    observation_data = observation_data,
    covariate_data = covariate_data,
    prepared_spatial_field = prepared_spatial_field
  )
  class(result) <- "downscale_fit"

  return(result)
}

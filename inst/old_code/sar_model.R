#' Fit Bayesian downscaling model using SAR prior
#'
#' This function wraps the complete fitting and inference pipeline using the SAR prior model.
#' You can provide either a `downscaling_model` object or the raw inputs (y, X, geometry).
#'
#' @param model A downscaling_model object from `make_downscaling_model()` (optional)
#' @param y Vector of coarse-resolution observations (optional if `model` is supplied)
#' @param X Design matrix at fine resolution (optional if `model` is supplied)
#' @param geometry Coarse-resolution sf geometry (optional if `model` is supplied)
#' @param target_grid An sf object containing the fine-resolution prediction grid
#' @param column The name of the column in the coarse data representing the observed variable
#' @param grid_shape Grid assumption: "rectangle" or "other"
#' @param include_variance Logical; whether to compute posterior variance
#' @param linear_solver Solver to use: "cg", "pcg", or "cholesky"
#' @param rho SAR spatial autocorrelation parameter (typically close to 1)
#' @param lambda Scalar or vector prior precision weight
#' @param only_intersects Logical; restrict grid to regions overlapping soundings
#' @param ... Additional arguments to pass into parameter preparation
#'
#' @return target_grid with `posterior_mean` column, and optionally `posterior_variance`
#' @export
fit_sar_downscale <- function(model = NULL,
                                   y = NULL,
                                   X = NULL,
                                   geometry = NULL,
                                   target_grid,
                                   column,
                                   grid_shape = c("other", "rectangle"),
                                   include_variance = FALSE,
                                   linear_solver = c("cg", "pcg", "cholesky"),
                                   rho = 0.99,
                                   lambda = 4,
                                   only_intersects = FALSE,
                                   ...) {

  grid_shape <- match.arg(grid_shape)
  linear_solver <- match.arg(linear_solver)

  if (!is.null(model)) {
    if (is.null(model$geometry)) stop("Model object lacks coarse geometry")
    soundings <- model$geometry
    soundings[[column]] <- model$y
    X <- model$X
    y <- model$y
  } else {
    if (is.null(y) || is.null(X) || is.null(geometry)) {
      stop("If no model is provided, y, X, and geometry must be specified")
    }
    soundings <- geometry
    soundings[[column]] <- y
  }

  if (only_intersects) {
    message("Restricting target grid to cells that intersect soundings")
    target_grid <- target_grid[sapply(sf::st_intersects(target_grid, soundings), function(x) length(x) > 0), ]
    grid_shape <- "other"
  }

  parameters <- complete_parameters(
    soundings = soundings,
    target_grid = target_grid,
    column = column,
    grid_shape = grid_shape,
    linear_solver = linear_solver,
    rho = rho,
    lambda = lambda,
    ...
  )

  message("Solving for posterior mean...")
  target_grid$posterior_mean <- do.call(estimate_target_grid, parameters)

  if (include_variance) {
    message("Computing posterior variance (conditional on posterior mode of sigma^2)...")
    a_n <- nrow(soundings) / 2 + 1
    quad <- crossprod(y) - crossprod(target_grid$posterior_mean, parameters$precision_matrix %*% target_grid$posterior_mean)
    b_n <- 1 + as.numeric(quad) / 2
    sigma2_mode <- b_n / (a_n + 1)

    target_grid$sigma2_mode <- sigma2_mode
    target_grid$posterior_variance <- sigma2_mode * compute_posterior_variance(parameters$precision_matrix)
  }

  return(target_grid)
}


#' Prepare parameters for Bayesian downscaling
#'
#' This function prepares the necessary matrices and hyperparameters to run the SAR-model-based
#' downscaling, including the aggregation matrix A, spatial weights W, and posterior precision.
#'
#' @param soundings An sf object with coarse-resolution observations
#' @param target_grid An sf object for the fine-resolution prediction grid
#' @param column The name of the column in `soundings` that holds the response
#' @param grid_shape One of "rectangle" or "other"
#' @param precision_matrix Optional; full precision matrix (overrides default construction)
#' @param W Optional; sparse neighbor matrix
#' @param A Optional; aggregation matrix from fine to coarse grid
#' @param AtA Optional; precomputed crossprod(A)
#' @param Atb Optional; precomputed A'b where b is the observed response
#' @param rho SAR autocorrelation parameter
#' @param tau SAR precision weight (scalar or vector)
#' @param linear_solver Solver for downstream posterior solve
#'
#' @return A list with `precision_matrix`, `Atb`, and `linear_solver`
#' @export
complete_parameters <- function(soundings,
                                target_grid,
                                column,
                                grid_shape = c("rectangle", "other"),
                                precision_matrix = NULL,
                                W = NULL,
                                A = NULL,
                                AtA = NULL,
                                Atb = NULL,
                                rho = 0.99,
                                tau = 4,
                                linear_solver = c("cg", "pcg", "cholesky")) {

  grid_shape <- match.arg(grid_shape)
  linear_solver <- match.arg(linear_solver)

  if (is.null(AtA) || is.null(Atb)) {
    message("At least one of AtA or Atb is missing. Computing...")
    if (is.null(A)) {
      message("Computing aggregation matrix A...")
      A <- compute_A_matrix(target_grid, soundings)
    }
    if (is.null(AtA)) {
      message("Computing AtA...")
      AtA <- Matrix::crossprod(A)
    }
    if (is.null(Atb)) {
      message("Computing Atb...")
      Atb <- Matrix::crossprod(A, soundings[[column]])
    }
  }

  if (is.null(precision_matrix)) {
    if (is.null(W)) {
      message("Computing spatial weights matrix W...")
      W <- compute_W_matrix(target_grid, grid_shape)
    }
    message("Computing posterior precision matrix...")
    I <- Matrix::Diagonal(nrow(W))
    if (length(tau) == 1) {
      tau <- tau * Matrix::crossprod(I - rho * W)
    } else {
      tau <- Matrix::crossprod(Matrix::Diagonal(x = sqrt(tau)) %*% (I - rho * W))
    }
    precision_matrix <- AtA + tau
  } else {
    precision_matrix <- AtA + precision_matrix
  }

  list(
    precision_matrix = precision_matrix,
    Atb = Atb,
    linear_solver = linear_solver
  )
}

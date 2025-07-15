### METHODS TO TRY:
# 1. Full-Bayesian
# 2. REML
# 3. Penalized REML
# 4. Cross-validation

#' Compute SAR Precision Matrix
#'
#' Computes the SAR prior precision matrix Q = (I - rho * W)'(I - rho * W)
#'
#' @param W A sparse row-normalized neighbor weight matrix (e.g., from compute_W_matrix)
#' @param rho Spatial dependence parameter, typically in [0, 1)
#'
#' @return A sparse precision matrix Q
#' @export
compute_SAR_precision <- function(W, rho) {
  stopifnot(methods::is(W, "Matrix"))
  stopifnot(methods::is(W, "dgCMatrix"))
  stopifnot(rho >= 0 && rho <= 1)

  I <- Matrix::Diagonal(n = nrow(W))
  I_minus_rhoW <- I - rho * W
  Q <- Matrix::crossprod(I_minus_rhoW)
  if (!inherits(Q, "dgCMatrix")) {
    Q <- Matrix::Matrix(Q, sparse = TRUE)
  }
  return(Q)
}

#' Complete Parameters for Downscaling
#'
#' Checks inputs and constructs needed matrices for downscaling estimation.
#'
#' @return A list of all model components for downstream estimation
#' @export
complete_parameters2 <- function(soundings,
                                 target_grid,
                                 column,
                                 method = c("bayes", "reml", "penalized_reml", "cv"),
                                 grid_shape = c("rectangle", "other"),
                                 A = NULL,
                                 AtA = NULL,
                                 Atb = NULL,
                                 W = NULL,
                                 Q = NULL,
                                 gamma = NULL,
                                 rho = 0.99,
                                 lambda = 1,
                                 linear_solver = c("cg", "pcg", "cholesky", "rcpp_cg"),
                                 tol = 1e-8,
                                 maxit = 1000) {
  method <- match.arg(method)
  grid_shape <- match.arg(grid_shape)
  linear_solver <- match.arg(linear_solver)

  if (is.null(A)) {
    message("Computing A matrix...")
    A <- compute_A_matrix(target_grid, soundings)
  }

  if (is.null(AtA)) {
    message("Computing AtA...")
    AtA <- Matrix::crossprod(A)
    AtA <- methods::as(Matrix::Matrix(AtA), "generalMatrix")
  }

  if (is.null(Atb)) {
    message("Computing Atb...")
    Atb <- Matrix::crossprod(A, soundings[[column]])
  }

  if (is.null(Q)) {
    if (is.null(W)) {
      message("Computing neighbor matrix W...")
      W <- compute_W_matrix(target_grid, grid_shape)

    }
    Q <- compute_SAR_precision(W, rho)
  }

  # Only construct precision_matrix if gamma is scalar
  precision_matrix <- if (!is.null(gamma) && length(gamma) == 1) {
    AtA + gamma * Q
  } else {
    NULL
  }

  return(list(
    A = A,
    AtA = AtA,
    Atb = Atb,
    W = W,
    Q = Q,
    gamma = gamma,
    rho = rho,
    lambda = lambda,
    precision_matrix = precision_matrix,
    linear_solver = linear_solver,
    method = method,
    tol = tol,
    maxit = maxit
  ))
}

call_estimator <- function(params) {
  linear_solver <- params$linear_solver

  args <- list(
    precision_matrix = params$precision_matrix,
    Atb = as.numeric(params$Atb),
    linear_solver = linear_solver,
    AtA = params$AtA
      )

  if (linear_solver == "rcpp_cg") {
    if (is.null(params$Q) || is.null(params$gamma)) stop("Q and gamma must be provided for rcpp_cg solver.")

    args$Q <- as(params$Q,"generalMatrix") #  No conversion needed
    args$gamma <- params$gamma  # ✅upports vector or scalar
    args$AtA <- Matrix::Matrix(params$AtA, sparse = TRUE)
    args$tol <- params$tol
    args$maxit <- params$maxit
  }

  do.call(estimate_target_grid, args)
}


#' Downscale using Bayesian Linear Model with SAR Prior
#'
#' @return An sf object containing posterior means (and optionally variances)
#' @export
downscale2 <- function(soundings,
                       target_grid,
                       column = "SIF_757nm",
                       include_variance = FALSE,
                       linear_solver = c("cg", "pcg", "cholesky", "rcpp_cg"),
                       grid_shape = c("other", "rectangle"),
                       only_intersects = FALSE,
                       gamma = NULL,
                       rho = 0.99,
                       tol = 1e-8,
                       maxit = 1000,
                       ...) {
  linear_solver <- match.arg(linear_solver)
  grid_shape <- match.arg(grid_shape)

  if (only_intersects) {
    message("Only grid cells intersecting the soundings will be used.")
    target_grid <- target_grid[sapply(sf::st_intersects(target_grid, soundings), function(x) length(x) > 0), ]
    message("Non-rectangular grid-shape assumed after removing non-intersecting cells.")
    grid_shape <- "other"
  }

  parameters <- complete_parameters2(soundings = soundings,
                                     target_grid = target_grid,
                                     grid_shape = grid_shape,
                                     linear_solver = linear_solver,
                                     column = column,
                                     gamma = gamma,
                                     rho = rho,
                                     tol = tol,
                                     maxit = maxit,
                                     ...)

  message("Computing posterior mean...")
  posterior <- call_estimator(parameters)
  stopifnot(is.numeric(posterior), length(posterior) == nrow(target_grid))
  target_grid$posterior_mean <- as.numeric(posterior)

  if (include_variance) {
    if (is.null(parameters$precision_matrix)) {
      warning("Variance estimation skipped because gamma is not scalar (precision matrix not available).")
    } else {
      message("Computing posterior variance conditional on posterior mode of sigma^2...")
      a_n <- dim(soundings)[1] / 2 + 1
      b_n <- 1 + (Matrix::crossprod(soundings[[column]]) -
                    Matrix::crossprod(target_grid$posterior_mean,
                                      parameters$precision_matrix %*% target_grid$posterior_mean)) / 2
      sigma2_mode <- as.vector(b_n / (a_n + 1))

      print(sigma2_mode)

      target_grid$sigma2_mode <- sigma2_mode
      target_grid$posterior_variance <- sigma2_mode * compute_posterior_variance(parameters$precision_matrix)
    }
  }

  return(target_grid)
}


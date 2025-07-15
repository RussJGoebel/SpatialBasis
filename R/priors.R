###################### Priors on Bases #########################################
## NOTE: each prior stores its parameters and useful cached objects.
## Then, the prior has a compute_precision method.
## If a gradient is possible,

#' Compute Precision Matrix for a Prior Object
#'
#' @param prior A prior object.
#' @export
compute_precision <- function(prior) {
  UseMethod("compute_precision")
}

#' Compute Gradient of the Prior Precision Matrix
#' @param prior A prior object.
#' @return A named list of sparse matrices: the derivative w.r.t. each hyperparameter.
compute_gradient <- function(prior) {
  UseMethod("compute_gradient")
}

#' Flat Prior (No Penalty)
#' @param basis A basis object (used only to get `k`)
#' @export
make_flat_prior <- function(basis) {
  stopifnot(inherits(basis, "basis"))
  k <- basis$k
  structure(list(k = k), class = c("flat", "prior"))
}

#' @export
compute_precision.flat <- function(prior) {
  Matrix::Diagonal(prior$k, 0)
}

#' @export
compute_gradient.flat <- function(prior) {
  list()  # no hyperparameters, no gradients
}


#' Ridge Prior (λ * I Penalty)
#' @param basis A basis object
#' @param tau Ridge penalty parameter
#' @export
make_ridge_prior <- function(basis, tau = 1) {
  stopifnot(inherits(basis, "basis"))
  if (tau <= 0) stop("tau must be positive.")
  k <- basis$k
  structure(
    list(k = k, tau = tau),
    class = c("ridge", "prior")
  )
}

#' @export
compute_precision.ridge <- function(prior) {
  Matrix::Diagonal(prior$k, prior$lambda)
}

#' @export
compute_gradient.ridge <- function(prior) {
  list(lambda = Matrix::Diagonal(prior$k, 1))
}


#' Make a SAR Prior
#'
#' Constructs a spatial autoregressive (SAR) prior specification.
#'
#' @param basis A `geometry_basis` object.
#' @param rho Spatial dependence parameter (|rho| <= 1).
#' @param tau Precision scaling parameter (>0).
#' @param W Optional neighborhood matrix. If NULL, computed from basis.
#' @param scaling_weights Optional variance scaling vector (length matching the number of geometries).
#'
#' @return An object of class c("sar", "prior").
#'
#' @export
make_sar_prior <- function(
    basis,
    rho = 0.99,
    tau = 1,
    W = NULL,
    scaling_weights = NULL
) {
  if (!inherits(basis, "geometry_basis")) {
    stop("make_sar_prior() requires a geometry_basis object.")
  }
  if (abs(rho) > 1) {
    stop("rho must be between -1 and 1.")
  }
  if (tau <= 0) {
    stop("tau must be positive.")
  }
  if (!is.null(scaling_weights)) {
    if (length(scaling_weights) != basis$k) {
      stop("scaling_weights must be a numeric vector with length matching the number of geometries.")
    }
  }
  n <- basis$k

  if (is.null(W)) {
    W <- compute_W_matrix(sf::st_as_sf(basis$geometry), grid_shape = "other")
  }

  # Always convert W to sparse Matrix
  W <- Matrix::Matrix(W)

  structure(
    list(
      type = "sar",
      structure = list(
        W = W,
        scaling_weights = scaling_weights,
        n = n
      ),
      hyperparameters = list(
        rho = rho,
        tau = tau
      ),
      precomputed = list(), # used for caching
      basis = basis
    ),
    class = c("sar", "prior")
  )
}

#' Compute Precision Matrix for SAR Prior
#'
#' Given a SAR prior object, returns the precision matrix Q.
#' Reuses cached components when possible.
#'
#' @param prior An object of class "sar".
#'
#' @return A sparse symmetric precision matrix (Matrix).
#'
#' @export
compute_precision.sar <- function(prior) {
  # Validate
  if (!inherits(prior, "sar")) {
    stop("compute_precision.sar() requires an object of class 'sar'.")
  }

  # Extract components
  rho <- prior$hyperparameters$rho
  tau <- prior$hyperparameters$tau
  W <- prior$structure$W
  scaling_weights <- prior$structure$scaling_weights
  n <- prior$structure$n

  if(rho > 1 | rho < -1){message("To be interpretable, rho should be in [-1,1].")}

  # Identity matrix
  I <- Matrix::Diagonal(n)

  # Check if we can reuse cached crossprod(M)
  reuse_cache <- (
    !is.null(prior$precomputed$rho) &&
      identical(prior$precomputed$rho, rho) &&
      !is.null(prior$precomputed$C)
  )

  if (reuse_cache) {
    # Use cached C = crossprod(M)
    C <- prior$precomputed$C
  } else {
    # Compute M
    M <- if (!is.null(scaling_weights)) {
      (I - rho * W) %*% Matrix::Diagonal(x = sqrt(scaling_weights))
    } else {
      (I - rho * W)
    }

    # Compute C = crossprod(M)
    C <- Matrix::crossprod(M)

    # Update cache (in-place assignment)
    prior$precomputed <- list(
      rho = rho,
      C = C
    )
  }

  # Final precision matrix
  Q <- tau * C

  Q
}

#' Compute Gradient of SAR Prior Precision Matrix
#'
#' Computes the derivative of the precision matrix w.r.t. each hyperparameter.
#'
#' @param prior An object of class "sar".
#'
#' @return A named list of sparse matrices.
#'
#' @export
compute_gradient.sar <- function(prior) {
  # Validate
  if (!inherits(prior, "sar")) {
    stop("compute_gradient.sar() requires an object of class 'sar'.")
  }

  # Extract components
  rho <- prior$hyperparameters$rho
  tau <- prior$hyperparameters$tau
  W <- prior$structure$W
  scaling_weights <- prior$structure$scaling_weights
  n <- prior$structure$n

  I <- Matrix::Diagonal(n)

  # Compute M
  M <- if (!is.null(scaling_weights)) {
    (I - rho * W) %*% Matrix::Diagonal(x = sqrt(scaling_weights))
  } else {
    (I - rho * W)
  }

  # dQ/dtau = crossprod(M)
  dQ_dtau <- Matrix::crossprod(M)

  # dM/drho
  dM_drho <- if (!is.null(scaling_weights)) {
    -W %*% Matrix::Diagonal(x = sqrt(scaling_weights))
  } else {
    -W
  }

  # Compute dQ/drho = tau * [ (dM_drho)^T M + M^T dM_drho ]
  term1 <- Matrix::crossprod(dM_drho, M)
  term2 <- Matrix::crossprod(M, dM_drho)
  dQ_drho <- tau * (term1 + term2)

  # Return as list
  list(
    tau = dQ_dtau,
    rho = dQ_drho
  )
}


#' Make a Matérn-Fourier Prior
#'
#' Constructs a prior over a Matérn-Fourier basis, parameterized by:
#' \itemize{
#'   \item \code{tau}: Precision scaling parameter (>0).
#'   \item \code{kappa}: Range parameter (>0), controls the decay of correlation.
#'   \item \code{alpha}: Smoothness parameter (>0).
#' }
#' Unlike simpler versions, this function does **not** precompute eigenvalues.
#' Instead, it stores the fixed Laplacian eigenvalues from the basis and computes
#' the eigenvalues dynamically during precision and gradient evaluation.
#'
#' @param basis A \code{function_basis} object created by \code{make_matern_fourier_basis()}.
#' @param tau Positive numeric, precision scaling parameter.
#' @param kappa Positive numeric, range parameter.
#' @param alpha Positive numeric, smoothness parameter.
#'
#' @return An object of class \code{c("matern_fourier", "prior")}.
#' @export
make_matern_fourier_prior <- function(
    basis,
    tau = 1,
    kappa = 1,
    alpha = 1
) {
  if (!inherits(basis, "function_basis") || basis$type != "matern_fourier") {
    stop("make_matern_fourier_prior() requires a basis of type 'matern_fourier'.")
  }
  if (tau <= 0) stop("tau must be positive.")
  if (kappa <= 0) stop("kappa must be positive.")
  if (alpha <= 0) stop("alpha must be positive.")

  structure(
    list(
      type = "matern_fourier",
      structure = list(
        laplacian_lambda = basis$metadata$laplacian_lambda,
        n = basis$k
      ),
      hyperparameters = list(
        tau = tau,
        kappa = kappa,
        alpha = alpha
      ),
      precomputed = list()  # Reserved for any future caching
    ),
    class = c("matern_fourier", "prior")
  )
}


#' Compute Precision Matrix for Matérn-Fourier Prior
#'
#' Given a Matérn-Fourier prior object, computes the precision matrix:
#' \deqn{
#'   Q = \tau \, \mathrm{diag}\Bigl[\bigl(\kappa^2 + \lambda_j^{Laplacian}\bigr)^{\alpha}\Bigr]
#' }
#' where \eqn{\lambda_j^{Laplacian}} are the stored Laplacian eigenvalues of the basis.
#'
#' @param prior An object of class "matern_fourier".
#'
#' @return A sparse diagonal precision matrix (class "Diagonal").
#' @export
#' Compute Precision Matrix for Matérn-Fourier Prior
#'
#' Computes the diagonal precision matrix:
#'   Q_j = tau * (kappa^2 + Laplacian_lambda_j)^alpha
#'
#' @param prior An object of class "matern_fourier".
#' @return A sparse diagonal precision matrix.
#' @export
compute_precision.matern_fourier <- function(prior) {
  if (!inherits(prior, "matern_fourier")) {
    stop("compute_precision.matern_fourier() requires an object of class 'matern_fourier'.")
  }

  tau <- prior$hyperparameters$tau
  kappa <- prior$hyperparameters$kappa
  alpha <- prior$hyperparameters$alpha
  lambda_lap <- prior$structure$laplacian_lambda
  n <- prior$structure$n

  # Compute diagonal entries: Q_j = tau * (kappa^2 + lambda_j)^alpha
  q_diag <- tau * (kappa^2 + lambda_lap)^alpha

  Matrix::Diagonal(n, x = q_diag)
}

#' Compute Gradient of Matérn-Fourier Prior Precision Matrix
#'
#' Computes the derivative of the precision matrix with respect to each hyperparameter:
#' \itemize{
#'   \item \code{tau}: scaling precision.
#'   \item \code{kappa}: range parameter.
#'   \item \code{alpha}: smoothness parameter.
#' }
#' All derivatives are returned as sparse diagonal matrices.
#'
#' The gradients are:
#' \deqn{
#'   \frac{\partial Q}{\partial \tau} = \mathrm{diag}\Bigl[\bigl(\kappa^2 + \lambda_j^{Laplacian}\bigr)^{\alpha}\Bigr]
#' }
#' \deqn{
#'   \frac{\partial Q}{\partial \kappa} = \tau \, \mathrm{diag}\Bigl[ \alpha \cdot 2\kappa \cdot \bigl(\kappa^2 + \lambda_j^{Laplacian}\bigr)^{\alpha - 1}\Bigr]
#' }
#' \deqn{
#'   \frac{\partial Q}{\partial \alpha} = \tau \, \mathrm{diag}\Bigl[\bigl(\kappa^2 + \lambda_j^{Laplacian}\bigr)^{\alpha} \cdot \log\bigl(\kappa^2 + \lambda_j^{Laplacian}\bigr)\Bigr]
#' }
#'
#' @param prior An object of class "matern_fourier".
#'
#' @return A named list with components:
#' \describe{
#'   \item{tau}{Derivative w.r.t. \code{tau}.}
#'   \item{kappa}{Derivative w.r.t. \code{kappa}.}
#'   \item{alpha}{Derivative w.r.t. \code{alpha}.}
#' }
#' @export
compute_gradient.matern_fourier <- function(prior) {
  if (!inherits(prior, "matern_fourier")) {
    stop("compute_gradient.matern_fourier() requires an object of class 'matern_fourier'.")
  }

  tau <- prior$hyperparameters$tau
  kappa <- prior$hyperparameters$kappa
  alpha <- prior$hyperparameters$alpha
  lambda_lap <- prior$structure$laplacian_lambda
  n <- prior$structure$n

  kappa_sq_plus_lap <- kappa^2 + lambda_lap

  # dQ/dtau
  dQ_dtau_diag <- kappa_sq_plus_lap^alpha
  dQ_dtau <- Matrix::Diagonal(n, x = dQ_dtau_diag)

  # dQ/dkappa
  dQ_dkappa_diag <- tau * alpha * 2 * kappa * kappa_sq_plus_lap^(alpha - 1)
  dQ_dkappa <- Matrix::Diagonal(n, x = dQ_dkappa_diag)

  # dQ/dalpha
  dQ_dalpha_diag <- tau * kappa_sq_plus_lap^alpha * log(kappa_sq_plus_lap)
  dQ_dalpha <- Matrix::Diagonal(n, x = dQ_dalpha_diag)

  list(
    tau = dQ_dtau,
    kappa = dQ_dkappa,
    alpha = dQ_dalpha
  )
}



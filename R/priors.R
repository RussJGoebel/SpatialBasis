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
    list(k = k,hyperparameters = list(
      rho = rho,
      tau = tau
    )),
    class = c("ridge", "prior")
  )
}

#' @export
compute_precision.ridge <- function(prior) {
  Matrix::Diagonal(prior$k, prior$hyperparameters$tau)
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
    scaling_weights = NULL,
    adjacency = "queen"
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
    W <- SpatialBasis::compute_W_matrix(sf::st_as_sf(basis$geometry), grid_shape = "other", adjacency = adjacency)
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
#'   tau: Precision scaling parameter (>0)
#'   kappa: Range parameter (>0)
#'   alpha: Smoothness parameter (>0)
#'
#' The constant mode precision (first entry) is recomputed as (kappa^2)^alpha if included.
#'
#' @param basis A function_basis created by make_matern_fourier_basis()
#' @param tau Positive numeric, precision scaling parameter
#' @param kappa2 Positive numeric, range parameter
#' @param alpha Positive numeric, smoothness parameter
#' @return An object of class c("matern_fourier", "prior")
#' @export
make_matern_fourier_prior <- function(
    basis,
    tau = 1,
    kappa2 = 1,
    alpha = 1
) {
  if (!inherits(basis, "function_basis") || basis$type != "matern_fourier") {
    stop("make_matern_fourier_prior() requires a basis of type 'matern_fourier'.")
  }
  if (tau <= 0) stop("tau must be positive.")
  if (kappa2 <= 0) stop("kappa2 must be positive.")
  if (alpha <= 0) stop("alpha must be positive.")

  constant_mode_included <- basis$metadata$constant_mode_included %||% FALSE
  constant_mode_precision <- if (constant_mode_included) tau*(kappa2)^alpha else NA_real_

  structure(
    list(
      type = "matern_fourier",
      structure = list(
        laplacian_lambda = basis$metadata$laplacian_lambda,
        n = basis$k,
        constant_mode_included = constant_mode_included
      ),
      hyperparameters = list(
        tau = tau,
        kappa2 = kappa2,
        alpha = alpha,
        constant_mode_precision = constant_mode_precision
      ),
      precomputed = list()
    ),
    class = c("matern_fourier", "prior")
  )
}



#' Compute Precision Matrix for Matérn-Fourier Prior
#'
#' Computes the diagonal precision matrix:
#'   Q_j = tau * (kappa^2 + Laplacian_lambda_j)^alpha
#'   The first entry (if constant_mode_included) may be overridden
#'
#' @param prior An object of class "matern_fourier".
#' @return A sparse diagonal precision matrix.
#' @export
compute_precision.matern_fourier <- function(prior) {
  if (!inherits(prior, "matern_fourier")) {
    stop("compute_precision.matern_fourier() requires an object of class 'matern_fourier'.")
  }

  tau <- prior$hyperparameters$tau
  kappa2 <- prior$hyperparameters$kappa2
  alpha <- prior$hyperparameters$alpha
  lambda_lap <- prior$structure$laplacian_lambda
  n <- prior$structure$n

  constant_mode_included <- prior$structure$constant_mode_included %||% FALSE
  constant_mode_precision <- prior$hyperparameters$constant_mode_precision %||% NA_real_

  q_diag <- tau * (kappa2 + lambda_lap)^alpha

  if (constant_mode_included && !is.na(constant_mode_precision)) {
    q_diag[1] <- constant_mode_precision
  }

  Matrix::Diagonal(n, x = q_diag)
}


#' Compute Gradient of Matérn-Fourier Prior Precision Matrix
#'
#' Computes the derivative of the precision matrix with respect to each hyperparameter:
#' \itemize{
#'   \item \code{tau}: scaling precision.
#'   \item \code{kappa2}: range parameter.
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
#'   \item{kappa2}{Derivative w.r.t. \code{kappa2}.}
#'   \item{alpha}{Derivative w.r.t. \code{alpha}.}
#' }
#' @export
compute_gradient.matern_fourier <- function(prior) {
  if (!inherits(prior, "matern_fourier")) {
    stop("compute_gradient.matern_fourier() requires an object of class 'matern_fourier'.")
  }

  tau <- prior$hyperparameters$tau
  kappa2 <- prior$hyperparameters$kappa2
  alpha <- prior$hyperparameters$alpha
  lambda_lap <- prior$structure$laplacian_lambda
  n <- prior$structure$n

  kappa_sq_plus_lap <- kappa2 + lambda_lap

  # dQ/dtau
  dQ_dtau_diag <- kappa_sq_plus_lap^alpha
  dQ_dtau <- Matrix::Diagonal(n, x = dQ_dtau_diag)

  # dQ/dkappa
  dQ_dkappa2_diag <- tau * alpha * 2 * kappa2 * kappa_sq_plus_lap^(alpha - 1)
  dQ_dkappa2 <- Matrix::Diagonal(n, x = dQ_dkappa2_diag)

  # dQ/dalpha
  dQ_dalpha_diag <- tau * kappa_sq_plus_lap^alpha * log(kappa_sq_plus_lap)
  dQ_dalpha <- Matrix::Diagonal(n, x = dQ_dalpha_diag)

  list(
    tau = dQ_dtau,
    kappa2 = dQ_dkappa2,
    alpha = dQ_dalpha
  )
}

### NYSTROM ====================================================================

#' Make a Nyström Prior
#'
#' Uses eigenvalues from the Nyström basis as diagonal precisions
#'
#' @param basis A basis object created by `make_nystrom_basis()`
#' @param tau A global precision multiplier (default = 1)
#' @export
make_nystrom_prior <- function(basis, tau = 1) {
  stopifnot(inherits(basis, "function_basis"))
  stopifnot("nystrom" %in% basis$type)
  if (tau <= 0) stop("tau must be positive.")

  lambda <- basis$metadata$eigenvalues
  if (is.null(lambda)) stop("No eigenvalues found in basis metadata.")

  structure(
    list(
      type = "nystrom",
      k = length(lambda),
      tau = tau,
      eigenvalues = lambda
    ),
    class = c("nystrom", "prior")
  )
}

#' @export
compute_precision.nystrom <- function(prior) {
  Matrix::Diagonal(prior$k, x = prior$tau * prior$eigenvalues)
}

#' @export
compute_gradient.nystrom <- function(prior) {
  list(tau = Matrix::Diagonal(prior$k, x = prior$eigenvalues))
}



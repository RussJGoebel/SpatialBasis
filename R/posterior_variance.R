#' Compute Posterior Variance of Spatial Coefficients
#'
#' Computes the posterior variance (or covariance matrix) of spatial coefficients
#' given a precision matrix or its Cholesky factor.
#'
#' @param precision_matrix The full posterior precision matrix (e.g., XtX + Lambda).
#' @param cholesky_factor Optional Cholesky factor of the precision matrix.
#' @param diagonal_only Logical: whether to return only the diagonal (marginal variances).
#' @param safe Logical: if TRUE, gracefully return NULL on failure; otherwise, stop.
#'
#' @return A vector of marginal variances if `diagonal_only = TRUE`, otherwise a full matrix.
#' @export
compute_posterior_variance <- function(precision_matrix,
                                       cholesky_factor = NULL,
                                       diagonal_only = TRUE,
                                       safe = TRUE) {
  if (!is.null(cholesky_factor)) {
    # Try chol2inv if Cholesky is available
    if (diagonal_only) {
      message("Computing marginal variances from Cholesky factor...")
      # Use solve(L, diag(n)) to get diag of inverse
      n <- ncol(cholesky_factor)
      I_n <- Matrix::Diagonal(n)
      L_inv <- try(Matrix::solve(cholesky_factor, I_n), silent = TRUE)
      if (inherits(L_inv, "try-error")) {
        if (safe) return(NULL) else stop("Cholesky inversion failed")
      }
      variances <- Matrix::colSums(L_inv^2)
      return(as.numeric(variances))
    } else {
      message("Computing full posterior covariance from Cholesky factor...")
      Sigma <- try(Matrix::chol2inv(cholesky_factor), silent = TRUE)
      if (inherits(Sigma, "try-error")) {
        if (safe) return(NULL) else stop("chol2inv failed")
      }
      return(Sigma)
    }
  } else {
    # Fallback to direct solve
    message("Computing posterior variance via matrix inverse...")
    if (diagonal_only) {
      # Diagonal only version via solve + column extraction
      n <- ncol(precision_matrix)
      I_n <- Matrix::Diagonal(n)
      V <- try(Matrix::solve(precision_matrix, I_n), silent = TRUE)
      if (inherits(V, "try-error")) {
        if (safe) return(NULL) else stop("Matrix inversion failed")
      }
      return(as.numeric(Matrix::diag(V)))
    } else {
      Sigma <- try(Matrix::solve(precision_matrix), silent = TRUE)
      if (inherits(Sigma, "try-error")) {
        if (safe) return(NULL) else stop("Matrix inversion failed")
      }
      return(Sigma)
    }
  }
}

#' compute_posterior_measurement_error
#'
#' @param posterior_mean The posterior mean calculated by fit_spatial_field() or fit_stacked_spatial_field()
#' @param X The design matrix computed by a prepare_ type function
#' @param Lambda The prior precision matrix before scaling
#' @param tau The scale factor on the prior precision matrix
#' @param y A vector of response values
#'
#' @return
#' @export
#'
#' @examples
compute_posterior_measurement_error <- function(posterior_mean, X, Lambda, y) {
  message("Computing posterior variance conditional on posterior mode of sigma^2...")

  a_0 <- 0 #1
  b_0 <- 0 # 1

  y <- y
  mu <- posterior_mean

  # Residual sum of squares
  residual <- y - X %*% mu
  residual_ss <- as.numeric(crossprod(residual))

  # Prior penalty term: muᵗ Λ mu
  penalty <- as.numeric(crossprod(mu, Lambda %*% mu))

  # Posterior shape and scale parameters
  a_n <- length(y) / 2 + a_0
  b_n <- b_0 + 0.5 * residual_ss + 0.5 * penalty

  # Posterior mode of sigma² (Inverse-Gamma)
  sigma2_mode <- b_n / (a_n + 1)

  cat(residual_ss, "\n")
  cat(b_n,"\n")


  message("sigma2_mode = ", sigma2_mode)

  cat("Penalty:", penalty, "\n")
  cat("RSS:", residual_ss, "\n")
  cat("Penalty/RSS ratio:", penalty / residual_ss, "\n")

  return(list(
    sigma2_mode = sigma2_mode,
    a_n = a_n,
    b_n = b_n,
    a_0 = a_0,
    b_0 = b_0
  ))
}


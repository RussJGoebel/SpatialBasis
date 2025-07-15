library(Matrix)
library(DownscalingPaper)
library(downscaling)


A <- compute_A_matrix(target_grid_paper, soundings_paper)
A <- as(A, "dgCMatrix")

W <- compute_W_matrix(target_grid_paper, "other")
Q <- Matrix::crossprod(Diagonal(nrow(W)) - W)
Q <- as(Q, "dgCMatrix")

y <- as.numeric(soundings$SIF_757nm)

cg_solve <- function(A, b, tol = 1e-6, maxiter = 1000) {
  # Ensure b is a matrix (required for Matrix::solve dispatch to work)
  if (is.vector(b)) {
    b <- matrix(b, ncol = 1)
  } else if (inherits(b, "Matrix") && is.null(dim(b))) {
    b <- Matrix::Matrix(b, ncol = 1)
  }
  # Sparse-aware solve
  res <- suppressWarnings(Matrix::solve(A, b, sparse = TRUE, tol = tol))
  return(res)
}


reml_objective_sparse <- function(log_gamma, A, y, Q, sigma2 = 1.0, log_prior_fn = NULL) {
  gamma <- exp(log_gamma)
  n <- ncol(A)
  g <- nrow(A)

  # Compute M = A^T A + gamma * Q
  AtA <- Matrix::crossprod(A)               # A^T A
  M <- AtA + gamma * Q              # M is positive definite even if Q is not

  # Solve M^{-1} A^T y
  At_y <- Matrix::crossprod(A, y)
  beta_hat <-  Matrix::solve(M, At_y, sparse = TRUE)

  # Compute residual r = y - A %*% beta_hat
  y_hat <- A %*% beta_hat
  r <- y - y_hat
  quad_form <- sum(r^2) / sigma2

  # Compute A %*% M^{-1} %*% A^T via solving each column
  # Let Z = solve(M, t(A)) --> A %*% Z = A %*% solve(M, t(A))
  Z <- Matrix::solve(M, t(A), sparse = TRUE)
  Sigma_y <- A %*% Z + Diagonal(g, sigma2)

  # Compute Cholesky and log determinant of Sigma_y
  Sigma_y_chol <- Cholesky(Sigma_y, LDL = FALSE)
  logdet_Sigma_y <- 2 * sum(log(diag(chol(Sigma_y))))

  # Total REML negative log-likelihood
  nll <- logdet_Sigma_y + quad_form

  # Optional log prior penalty on log(gamma)
  if (!is.null(log_prior_fn)) {
    nll <- nll - log_prior_fn(log_gamma)
  }

  return(nll)
}



reml_obj_grad_sparse <- function(log_gamma, A, y, Q, sigma2 = 1.0, log_prior_fn = NULL) {
  gamma <- exp(log_gamma)
  n <- ncol(A)
  g <- nrow(A)

  AtA <- crossprod(A)
  M <- Matrix::Matrix(AtA + gamma * Q)

  # Solve for beta_hat = M^{-1} A^T y
  At_y <- crossprod(A, y)
  beta_hat <- cg_solve(M, At_y)

  y_hat <- A %*% beta_hat
  r <- y - y_hat
  quad_form <- sum(r^2) / sigma2

  # Precompute Z = M^{-1} A^T
  Z <- cg_solve(M, t(A))  # This is n x g
  AZ <- A %*% Z           # This is g x g
  Sigma_y <- AZ + Diagonal(g, sigma2)

  # Cholesky + log det
  Sigma_y_chol <- Cholesky(Sigma_y, LDL = FALSE)
  logdet_Sigma_y <- 2 * sum(log(diag(chol(Sigma_y))))

  # REML objective
  nll <- logdet_Sigma_y + quad_form

  # Gradient:
  # dSigma_y/dgamma = -A M^{-1} Q M^{-1} A^T = -A (W) A^T
  W <- cg_solve(M, Q %*% Z)   # W = M^{-1} Q M^{-1} A^T
  dSigma_dgamma <- A %*% W    # g x g

  S_y_inv <- solve(Sigma_y_chol, solve(t(Sigma_y_chol), Diagonal(g)))

  trace_term <- sum(S_y_inv * dSigma_dgamma)
  alpha <- solve(Sigma_y_chol, solve(t(Sigma_y_chol), y))
  quad_term <- as.numeric(crossprod(alpha, dSigma_dgamma %*% alpha))

  grad_gamma <- -(trace_term - quad_term)  # d/d gamma
  grad_log_gamma <- grad_gamma * gamma     # chain rule for log gamma

  # Add optional log prior
  if (!is.null(log_prior_fn)) {
    nll <- nll - log_prior_fn(log_gamma)
  }

  return(list("value" = nll, "gradient" = grad_log_gamma))
}


##



res <- optim(
  par = log(1),
  fn = reml_objective_sparse,
  A = A,
  y = y,
  Q = Q,
  sigma2 = 1.0,
  method = "L-BFGS-B",
  lower = log(1e-3),
  upper = log(1e3),
  control = list(trace = 1)
)

gamma_hat <- exp(res$par)
cat("Estimated gamma (REML):", gamma_hat, "\n")

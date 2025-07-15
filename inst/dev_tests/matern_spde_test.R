#' Match SPDE Matérn Parameters to SAR Model (Normalized, Precision = (I - ρW)'(I - ρW))
#'
#' For a regular rectangular grid with normalized 4‑neighbor weights, the SAR precision
#' is Q = (I - ρ W)' (I - ρ W) / σ².  The eigenvalues of Q are (1 - ρ λ_W)² / σ², where λ_W
#' are the eigenvalues of W.  This function matches those covariance eigenvalues to the
#' SPDE Matérn covariance spectrum by least‑squares fitting on the log‑eigenvalues of the
#' lowest‑frequency modes.
#'
#' @param nx,ny Grid dimensions.
#' @param rho   SAR parameter (0 < rho < 1 with normalized weights).
#' @param sigma2 Initial variance guess for σ².
#' @param k_fit Number of low‑frequency eigenpairs used in the fit.
#' @return List(kappa, sigma2, loss, sar_eigen, matern_eigen) and a diagnostic plot.
#' @export
match_spde_to_sar <- function(nx, ny, rho, sigma2 = 1, k_fit = 500) {
  if (rho >= 1 || rho <= 0) stop("rho must be in (0,1) for normalized 4‑neighbor weights.")

  # Laplacian eigenvalues with Dirichlet BCs
  i <- 1:nx
  j <- 1:ny
  lambda_x <- 2 * (1 - cos(pi * i / (nx + 1)))
  lambda_y <- 2 * (1 - cos(pi * j / (ny + 1)))
  lambda_grid <- as.vector(outer(lambda_x, lambda_y, `+`))

  # Eigenvalues of normalized adjacency matrix W (largest = 1)
  W_eig <- lambda_grid / 4

  # SAR covariance eigenvalues:  σ² / (1 - ρ λ_W)²
  denom <- 1 - rho * W_eig
  if (any(abs(denom) < .Machine$double.eps)) stop("rho too close to 1/max(W_eig); precision becomes singular.")
  sar_eigen <- sigma2 / (denom^2)

  # Sort by increasing frequency (in lambda_grid)
  freq_order <- order(lambda_grid)
  lambda_sorted <- lambda_grid[freq_order]
  sar_sorted <- sar_eigen[freq_order]

  k_fit <- min(k_fit, length(sar_sorted))

  # Joint optimization over log(kappa) and log(sigma2)
  match_loss <- function(par) {
    log_kappa <- par[1]
    log_sigma2 <- par[2]
    kappa2 <- exp(2 * log_kappa)
    sigma2_m <- exp(log_sigma2)
    matern_eig <- sigma2_m / (kappa2 + lambda_sorted)^2
    sum((log(sar_sorted[1:k_fit]) - log(matern_eig[1:k_fit]))^2)
  }

  opt <- optim(par = c(0, log(sigma2)), fn = match_loss)
  kappa_est <- exp(opt$par[1])
  sigma2_est <- exp(opt$par[2])
  matern_eig_fit <- sigma2_est / (kappa_est^2 + lambda_sorted)^2

  # Plot
  plot(log(sar_sorted[1:k_fit]), type = "l", col = "blue", lwd = 2,
       ylab = "log eigenvalue", xlab = "Index (low → high freq)",
       main = "SPDE Matérn vs Normalized SAR Spectrum")
  lines(log(matern_eig_fit[1:k_fit]), col = "red", lwd = 2, lty = 2)
  legend("topright", legend = c("SAR (normalized)", "Fitted Matérn"),
         col = c("blue", "red"), lty = c(1, 2), lwd = 2)

  list(
    kappa = kappa_est,
    sigma2 = sigma2_est,
    loss = opt$value,
    sar_eigen = sar_sorted[1:k_fit],
    matern_eigen = matern_eig_fit[1:k_fit]
  )
}

# Example:
# result <- match_spde_to_sar(40, 40, rho = 0.95)
# result$kappa

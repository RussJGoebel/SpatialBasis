#' Create Matérn KL basis functions on a rectangular domain using sine eigenfunctions
#'
#' @param Lx Length of domain in x-direction
#' @param Ly Length of domain in y-direction
#' @param nx Number of basis functions in x-direction (i index)
#' @param ny Number of basis functions in y-direction (j index)
#' @param kappa Range-like parameter for Matérn (higher = smoother)
#' @param nu Smoothness parameter for Matérn
#'
#' @return A `downscaling_basis` object with evaluate and metadata
#' @export
make_matern_kl_basis <- function(Lx, Ly, nx, ny, kappa = 1, nu = 1.5) {
  stopifnot(Lx > 0, Ly > 0, nx >= 1, ny >= 1)

  # List of basis indices
  basis_ids <- expand.grid(i = 1:nx, j = 1:ny)
  k <- nrow(basis_ids)

  # Precompute eigenvalues (scaled variances)
  lambda <- with(basis_ids, {
    ((kappa^2 + (pi * i / Lx)^2 + (pi * j / Ly)^2))^(-nu - 1)
  })

  # Evaluation function: accepts a matrix of locations (s, t)
  evaluate <- function(coords) {
    if (!is.matrix(coords) || ncol(coords) != 2) stop("coords must be n x 2 matrix")
    s <- coords[, 1]
    t <- coords[, 2]

    phi <- matrix(NA_real_, nrow = nrow(coords), ncol = k)

    for (idx in seq_len(k)) {
      i <- basis_ids$i[idx]
      j <- basis_ids$j[idx]
      phi[, idx] <- sin(pi * i * s / Lx) * sin(pi * j * t / Ly)
    }

    phi
  }

  # Return as basis object
  structure(
    list(
      type = "matern_fourier_sin",
      evaluate = evaluate,
      k = k,
      lambda = lambda,
      domain = list(Lx = Lx, Ly = Ly),
      basis_ids = basis_ids,
      metadata = list(
        resolution = paste0(round(Lx / nx, 2), " x ", round(Ly / ny, 2), " spacing"),
        kappa = kappa,
        nu = nu
      )
    ),
    class = "downscaling_basis"
  )
}

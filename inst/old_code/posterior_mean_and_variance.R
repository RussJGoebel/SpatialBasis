#' estimate_target_grid
#'
#' @param precision_matrix matrix calculated using compute_W_matrix()
#' @param Atb Matrix multiplcation of A matrix and sounding vlaues (b)
#' @param linear_solver specific linear solver for estimation
#'
#' @return posterior_mean of downscaled estimate
#' @export
#'
#' @examples posterior_mean <- do.call(estimate_target_grid,parameters)
estimate_target_grid <- function(precision_matrix,Atb,linear_solver = c("cholesky","cg","pcg")){

  linear_solver <- match.arg(linear_solver)

  message("Using ",linear_solver," method...")

  if(linear_solver == "cholesky") value <- solve_system_cholesky(precision_matrix,Atb)
  if(linear_solver == "cg") value  <- cPCG::cgsolve(Matrix::as.matrix(precision_matrix),Matrix::as.matrix(Atb))
  if(linear_solver == "pcg") value  <- cPCG::pcgsolve(Matrix::as.matrix(precision_matrix),Matrix::as.matrix(Atb))

  return(value)


}

#' compute_posterior_variance
#'
#' @param precision_matrix matrix calculated using compute_W_matrix()
#'
#' @return diagonalized posterior variance matrix
#' @export
#'
#' @examples posterior_variance <- computer_posterior_variance(precision_matrix)
compute_posterior_variance <- function(precision_matrix){
  diag(solve(precision_matrix))
}

#' estimate_target_grid
#'
#' @param precision_matrix
#' @param Atb
#' @param linear_solver
#'
#' @return
#' @export
#'
#' @examples
solve_system <- function(precision_matrix,Atb,linear_solver = c("cg","pcg","cholesky")){

  linear_solver <- match.arg(linear_solver)

  if(linear_solver == "cholesky") value <- solve_system_cholesky(precision_matrix,Atb)
  if(linear_solver == "cg") value  <- cPCG::cgsolve(Matrix::as.matrix(precision_matrix),Matrix::as.matrix(Atb))
  if(linear_solver == "pcg") value  <- cPCG::pcgsolve(Matrix::as.matrix(precision_matrix),Matrix::as.matrix(Atb),preconditioner = "ICC")

  return(value)


}

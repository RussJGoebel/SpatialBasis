#' solve_system_cholesky
#' Solves a system of equations using the cholesky decomposition
#'
#' @param X A positive definite matrix
#' @param y A vector
#'
#' @return The solution $b$ to $Xb = y$.
#' @export
#'
#' @examples
#'
#' # Using solve_system_cholesky to perform least squares:
#'
#' # Simulate data
#' set.seed(1)
#'
#' n <- 1000
#' p <- 10#'
#' x <- matrix(rnorm(n*p),ncol = p)
#' b <- 1:10
#' y <- x %*% b +rnorm(n)
#'
#' # Use solve_system_cholesky to solve normal equations
#' bhat <- solve_system_cholesky(crossprod(x),crossprod(x,y))
#'
#' # Compare to linear model result
#' lm(y~x-1)
#'
#'
solve_system_cholesky <- function(X,y){

  R <- chol(X)
  b <- backsolve(R,forwardsolve(t(R),y))

  return(b)

}


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
estimate_target_grid <- function(precision_matrix,Atb,linear_solver = c("cg","pcg","cholesky")){

  linear_solver <- match.arg(linear_solver)

  if(linear_solver == "cholesky") value <- solve_system_cholesky(precision_matrix,Atb)
  if(linear_solver == "cg") value  <- cPCG::cgsolve(Matrix::as.matrix(precision_matrix),Matrix::as.matrix(Atb))
  if(linear_solver == "pcg") value  <- cPCG::pcgsolve(Matrix::as.matrix(precision_matrix),Matrix::as.matrix(Atb),preconditioner = "ICC")

  return(value)


}

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
estimate_target_grid_with_prior <- function(precision_matrix,Atb,linear_solver = c("cg","pcg","cholesky")){

  prior_mean <- rep(1,length(Atb))
  prior_covariance <-

  linear_solver <- match.arg(linear_solver)

  if(linear_solver == "cholesky") value <- solve_system_cholesky(precision_matrix,Atb)
  if(linear_solver == "cg") value  <- cPCG::cgsolve(Matrix::as.matrix(precision_matrix),Matrix::as.matrix(Atb))
  if(linear_solver == "pcg") value  <- cPCG::pcgsolve(Matrix::as.matrix(precision_matrix),Matrix::as.matrix(Atb),preconditioner = "ICC")

  return(value)


}



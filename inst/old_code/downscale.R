#' complete parameters
#' checks the parameters of the downscaling function and computes missing variables
#'
#' @param soundings An sf object containing sounding observations.
#' @param target_grid An sf object containing the underlying fine grid.
#' @param grid_shape One of "rectangle" or "other". If "rectangle" is specified, then a rectangular grid within the bounding box of `soundings` is created. If "other" is specified, only cells overlapping with the soundings are used.
#' @param precision_matrix The posterior precision matrix given by A'A+Lambda, where Lambda is the prior matrix.
#' @param W A sparse matrix of neighbor weights. The SAR model has covariance sigma^2(I-rho*W)(I-rho*W)^T.
#' @param A A sparse matrix of weights relating the target grid to the soundings.
#' @param AtA A'A with A above. (Matrix::crosprod(A))
#' @param Atb Where b is the values over the soundings, A'b (Matrix::crossprod(A,b))
#' @param rho A number between 0 and 1.
#' @param lambda A weight on the covariance. Closer to 0 means the fit is closer to a diagonal covariance.
#' @param linear_solver One of "cg" (conjugate gradient), "pcg" (preconditioned conjugate gradient) or "cholesky."
#'
#' @return list containing all parameteres needed for downscale()
#' @export
#'
#' @examples params <- complete_parameters(soundings, target_grid, "SIF_757nm")
complete_parameters <- function(soundings,
                                target_grid,
                                column,
                                grid_shape = c("rectangle","other"),
                                precision_matrix = NULL,
                                W = NULL,
                                A = NULL,
                                AtA = NULL,
                                Atb = NULL,
                                rho = 0.99,
                                lambda = 4,
                                linear_solver = c("cg","pcg","cholesky")){

  grid_shape <- match.arg(grid_shape)

  linear_solver <- match.arg(linear_solver)

  if(is.null(AtA) | is.null(Atb)){
    message("At least one of AtA or Atb is not supplied.")
    if(is.null(A)){
      message("Computing A...")
      A <- compute_A_matrix(target_grid,soundings)
    }
    if(is.null(AtA)){
      message("Computing AtA...")
      AtA <- Matrix::crossprod(A)
    }
    if(is.null(Atb)){
      message("Computing Atb...")
      Atb <- Matrix::crossprod(A,soundings[[column]])
    }
  }

  if(is.null(precision_matrix)){
    if(is.null(W)){
      message("Computing weight matrix...")
      W <- compute_W_matrix(target_grid,grid_shape)
    }
    if(length(lambda) == 1){message("Computing posterior precision matrix using lambda = ",lambda); precision_matrix <- AtA+lambda*Matrix::crossprod(Matrix::diag(dim(W)[1])-rho*W)
    }else{precision_matrix <- AtA+crossprod(diag(sqrt(lambda)) %*% (Matrix::diag(dim(W)[1])-rho*W))}


  }else{
    precision_matrix <- AtA+precision_matrix
  }

  parameters <- list(
    precision_matrix = precision_matrix,
    Atb = Atb,
    linear_solver = linear_solver
  )

}

# downscale <- function(target_grid,
#                       precision_matrix,
#                       Atb,
#                       linear_solver = c("pcg","cholesky")){
#
#   linear_solver <- match.arg(linear_solver)
#
#   target_grid$posterior_mean <- estimate_target_grid(precision_matrix,Atb,linear_solver)
#
# }

#' Title
#'
#' @param soundings An sf object containing sounding observations.
#' @param target_grid An sf object containing the underlying fine grid.
#' @param include_variance A boolean specifying whether or not to include the posterior variance.
#' @param linear_solver  One of "cg" (conjugate gradient), "pcg" (preconditioned conjugate gradient) or "cholesky."
#' @param grid_shape One of "rectangle" or "other". If "rectangle" is specified, then a rectangular grid within the bounding box of `soundings` is created. If "other" is specified, only cells overlapping with the soundings are used.
#' @param ... Additional arguments to be passed to complete the parameters as needed.
#'
#' @return sf_object with shape of target_grid containing posterior_means from bayesian downscaling technique
#' @export
#'
#' @examples d <- downscale(soundings, target_grid, column="XCO2_ppm", A=APG, lambda=10)
downscale <- function(soundings,
                      target_grid,
                      column = "SIF_757nm",
                      include_variance = FALSE,
                      linear_solver = c("cg","pcg","cholesky"),
                      grid_shape = c("other","rectangle"),
                      only_intersects = FALSE,
                      ...){

  linear_solver <- match.arg(linear_solver)
  grid_shape <- match.arg(grid_shape)

  if(only_intersects){
    message("Only grid cells intersecting the soundings will be used. (Set only_intersects to false to include more.)")
    target_grid <- target_grid[sapply(sf::st_intersects(target_grid,soundings),function(x){length(x)>0}), ]

    message("Non-rectangular grid-shape assumed after removing non-intersecting cells.")
    grid_shape <- "other"
  }

  parameters <- complete_parameters(soundings = soundings,
                                    target_grid = target_grid,
                                    grid_shape = grid_shape,
                                    linear_solver = linear_solver,
                                    column = column,
                                    ...)

  message("Computing posterior mean...")
  target_grid$posterior_mean <- do.call(estimate_target_grid,parameters)

  if(include_variance){
    message("Computing posterior variance conditional on posterior mode of sigma^2...");
    # Sigma^2 follows a gamma distribution with parameters an, bn.

    a_n <- dim(soundings)[1]/2+1
    b_n <- 1+(Matrix::crossprod(soundings[[column]])-Matrix::crossprod(target_grid$posterior_mean,
                                                                     parameters$precision_matrix %*% target_grid$posterior_mean))/2
    sigma2_mode <- as.vector((b_n)/(a_n+1))

    print(sigma2_mode)

    target_grid$sigma2_mode <- sigma2_mode
    target_grid$posterior_variance <-  sigma2_mode*compute_posterior_variance(parameters$precision_matrix)
    #target_grid$cov_matrix <- sigma2_mode*


  }

  ## CREATING DATASTRUCTURE
  #downscaled_result <- list(
  #  target_grid = target_grid,
  #  soundings = soundings
  #)

  return(target_grid)

}



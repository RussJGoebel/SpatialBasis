#' k_folds_cv_lambda_two_class
#'
#' @param soundings
#' @param target_grid
#' @param candidate_lambda
#' @param k
#' @param reps
#' @param column
#' @param grid_shape
#' @param precision_matrix
#' @param W
#' @param A
#' @param AtA
#' @param Atb
#' @param rho
#' @param lambda
#' @param linear_solver
#' @param mc.cores
#'
#' @return
#' @export
#'
#' @examples
k_folds_cv_lambda_two_class <- function(soundings,
                                        target_grid,
                                        candidate_lambda,
                                        k,
                                        reps,
                                        column = "SIF_757nm",
                                        grid_shape = c("rectangle","other"),
                                        precision_matrix = NULL,
                                        W = NULL,
                                        A = NULL,
                                        AtA = NULL,
                                        Atb = NULL,
                                        rho = 0.99,
                                        lambda = 4,
                                        linear_solver = c("cg","pcg","cholesky"),
                                        mc.cores = 16){

  ### Precompute parameters ####################################################
  grid_shape <- match.arg(grid_shape)
  linear_solver <- match.arg(linear_solver)

  if(any(is.na(soundings[[column]]))){
    message("Removed ",sum(is.na(soundings[[column]]))," soundings with missing values.")
    soundings <- soundings[!is.na(soundings[[column]]),]
  }

  if(is.null(A)){
    message("Computing A...")
    A <- compute_A_matrix(target_grid,soundings)
  }
  if(is.null(W)){
    message("Computing weight matrix...")
    W <- compute_W_matrix(target_grid,grid_shape = "other")
  }

  ### /Precompute Parameters ###################################################

  n <- dim(soundings)[1]
  fold_indices <- rep(1:k,
                      length.out = n)[sample(1:n)]

  compute_mse_function <- function(x){compute_mean_squared_error_over_folds(fold_indices,
                                                                            A,
                                                                            W,
                                                                            soundings,
                                                                            target_grid,
                                                                            column = column,
                                                                            x)}
  message("Using k-folds CV in parallel to pick overall lambda...")
  mean_square_error <- parallel::mclapply(candidate_lambda,compute_mse_function,
                                          mc.cores = mc.cores)

  best_lambda <- candidate_lambda[which.min(mean_square_error)]
  message("Best lambda found to be ",best_lambda)

  compute_mse_function <- function(x){compute_mean_squared_error_over_folds(fold_indices,
                                                                            A,
                                                                            W,
                                                                            soundings,
                                                                            target_grid,
                                                                            column = column,
                                                                            ifelse(target_grid$dominant_class == 6,x,best_lambda))}

  message("Using k-folds CV in parallel to pick water lambda...")
  mean_square_error_water <- parallel::mclapply(candidate_lambda,compute_mse_function,
                                          mc.cores = mc.cores)

  best_water_lambda <- candidate_lambda[which.min(mean_square_error_water)]
  best_lambda <- candidate_lambda[which.min(mean_square_error_water)]
  message("Best lambda found to be ",best_lambda)



  return(list(best_lambda = best_lambda,
              best_water_lambda = best_water_lambda,
              mean_square_error = mean_square_error,
              mean_square_error_water = mean_square_error_water))


}

#' k_folds_cv_lambda
#'
#' @param soundings
#' @param target_grid
#' @param candidate_lambda
#' @param k
#' @param reps
#' @param column
#' @param grid_shape
#' @param precision_matrix
#' @param W
#' @param A
#' @param AtA
#' @param Atb
#' @param rho
#' @param lambda
#' @param linear_solver
#' @param mc.cores
#'
#' @return
#' @export
#'
#' @examples
k_folds_cv_lambda <- function(soundings,
                              target_grid,
                              candidate_lambda,
                              k,
                              reps,
                              column = "SIF_757nm",
                              grid_shape = c("rectangle","other"),
                              precision_matrix = NULL,
                              W = NULL,
                              A = NULL,
                              AtA = NULL,
                              Atb = NULL,
                              rho = 0.99,
                              lambda = 4,
                              linear_solver = c("cg","pcg","cholesky"),
                              mc.cores = 16){

  ### Precompute parameters ####################################################
  grid_shape <- match.arg(grid_shape)
  linear_solver <- match.arg(linear_solver)

  if(any(is.na(soundings[[column]]))){
    message("Removed ",sum(is.na(soundings[[column]]))," soundings with missing values.")
    soundings <- soundings[!is.na(soundings[[column]]),]
    }

  if(is.null(A)){
    message("Computing A...")
    A <- compute_A_matrix(target_grid,soundings)
  }
  if(is.null(W)){
    message("Computing weight matrix...")
    W <- compute_W_matrix(target_grid,grid_shape = "other")
  }

  ### /Precompute Parameters ###################################################

  n <- dim(soundings)[1]
  fold_indices <- rep(1:k,
                      length.out = n)[sample(1:n)]

  compute_mse_function <- function(x){compute_mean_squared_error_over_folds(fold_indices,
                                                                            A,
                                                                            W,
                                                                            soundings,
                                                                            target_grid,
                                                                            column = column,
                                                                            x)}
  message("Using k-folds CV in parallel...")
  mean_square_error <- parallel::mclapply(candidate_lambda,compute_mse_function,
                     mc.cores = mc.cores)


  return(mean_square_error)


}


###############################################################################


#' Title compute_mean_squared_error_over_folds
#'
#' @param fold_indices
#' @param A
#' @param W
#' @param soundings
#' @param target_grid
#' @param column
#' @param lambda
#'
#' @return
#' @export
#'
#' @examples
compute_mean_squared_error_over_folds <- function(fold_indices,
                                                  A,
                                                  W,
                                                  soundings,
                                                  target_grid,
                                                  column = "SIF_757nm",
                                                  lambda
){



  mean_square_errors <- numeric(length(unique(fold_indices)))

  for(jj in unique(fold_indices)){

    message("Computing AtA...")
    AtA <- Matrix::crossprod(A[fold_indices != jj,])

    message("Computing Atb...")
    Atb <- Matrix::crossprod(A[fold_indices != jj,],soundings[[column]][fold_indices != jj])

    d <- downscale(soundings[fold_indices != jj,],
                   target_grid,
                   column = column,
                   AtA = AtA,
                   Atb = Atb,
                   W = W,
                   lambda = lambda)

    pred <- A %*% d$posterior_mean
    mean_square_errors[jj] <- mean((pred[fold_indices == jj]-soundings[[column]][fold_indices == jj])^2)

    message("Finished with fold ",jj,". Mean square error: ",mean_square_errors[jj])


  }

  return(mean(mean_square_errors))

}

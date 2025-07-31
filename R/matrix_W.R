#' Title
#'
#' @param nrow
#' @param ncol
#'
#' @return
#' @export
#'
#' @examples
generate_neighbor_matrix <- function(nrow,ncol){

  n <- nrow*ncol

  get_index <- function(i,j){
    # from coordinates i,j on the lattice, obtain index on lattice
    ncol*(i-1)+j
  }

  coords <- expand.grid(1:nrow,1:ncol)
  indices_of_ones <- apply(coords,1,function(x){
    indices <- outer(c((x[1]-1):(x[1]+1)),c((x[2]-1):(x[2]+1)),FUN = get_index)

    remove_indices <- c(-5,
                        -c(1,4,7)*(x[1]==1), #remove pixels above if x coordinate is on the top row of the lattice
                        -c(3,6,9)*(x[1]==nrow), #remove pixels below if x coordinate is on the bottom row of the lattice
                        -c(1:3)*(x[2]==1), # remove pixels to the left if x coordinate is on the left column of the lattice
                        -c(7:9)*(x[2]==ncol)) #remove pixels to the right if x coordinate is on the right column of the lattice
    indices <- indices[remove_indices]

  })

  indices_of_ones <- sapply(seq_along(indices_of_ones),function(x){cbind(get_index(coords[x,1],coords[x,2]),indices_of_ones[[x]])})
  indices_of_ones <- do.call(rbind,indices_of_ones)

  M <- Matrix::Matrix(0,nrow = n,ncol = n, sparse = TRUE)
  M[indices_of_ones] <- 1

  #M <- Matrix::Matrix(0,nrow = n, ncol = n)

  return(M)


}

#' Compute Neighborhood (W) Matrix for Spatial Grid
#'
#' @param target_grid An `sf` object representing the spatial grid.
#' @param grid_shape Either `"rectangle"` or `"other"`. Determines method of neighbor calculation.
#' @param adjacency `"queen"` (default) or `"rook"`; only relevant if `grid_shape = "other"`.
#'
#' @return A sparse row-normalized neighbor matrix `W`.
#' @export
compute_W_matrix <- function(target_grid,
                             grid_shape = c("rectangle", "other"),
                             adjacency = c("queen", "rook")) {

  grid_shape <- match.arg(grid_shape)
  adjacency <- match.arg(adjacency)

  if (grid_shape == "rectangle") {
    n_row <- max(target_grid$row)
    n_col <- max(target_grid$col)

    W <- generate_neighbor_matrix(n_row, n_col)  # assumes this returns binary matrix
    W <- row_normalize_matrix(W)

  } else if (grid_shape == "other") {
    message(sprintf("Computing neighbor matrix with %s adjacency.", adjacency))

    queen <- (adjacency == "queen")

    nb <- spdep::poly2nb(target_grid, queen = queen)
    mat <- spdep::nb2mat(nb, style = "W", zero.policy = TRUE)
    W <- Matrix::Matrix(mat, sparse = TRUE)
  }

  return(W)
}

row_normalize_matrix <- function(M){

  M_rowsums <- Matrix::rowSums(M)
  adjusted_M_rowsums <- ifelse(M_rowsums == 0,1,M_rowsums)
  sweep(M,1,adjusted_M_rowsums,"/")

}

################################################################################
## W matrices with different covariance structures
################################################################################

#' compute_dominant_class_W_matrix
#'
#' @param target_grid An sf object
#' @param class An integer from 1-8. 6 represents water.
#' @param epsilon A small nonzero value representing the weight of correlation between pixels with different dominant classes.
#'
#' @return
#' @export
#'
#' @examples
compute_dominant_class_W_matrix <- function(target_grid,
                                            class = 6,
                                            epsilon = 0.01){

  if(!("dominant_class" %in% names(target_grid))){stop("A 'dominant_class' variable is required.")}

  W <- compute_W_matrix(target_grid, grid_shape = "other")

  dominant <- target_grid$dominant

  scale_function <- function(x,y){
    return(epsilon*(x == class & y != class)+epsilon*(x !=class & y == class)+1*(x !=class& y !=class)+1*(x == class & y ==class ))
  }

  O <- outer(dominant,dominant, FUN = scale_function)

  W <- row_normalize_matrix(W*O)

  return(W)

}

compute_proportion_class_W_matrix <- function(target_grid,
                                              epsilon = 0.01){


  if(!("proportion_6" %in% names(target_grid))){stop("A 'proportion_6' variable is required.")}

  W <- compute_W_matrix(target_grid, grid_shape = "other")

  prop_water <- target_grid$proportion_6
  prop_water <- ifelse(is.na(prop_water),0,prop_water)


  scale_function_proportion <- function(x,y){
    return(pmax(1-abs(x-y),epsilon))
  }

  O <- outer(prop_water,prop_water, FUN = scale_function_proportion)

  W <- row_normalize_matrix(W*O)

  return(W)


}

compute_dominant_class_lambda <- function(dominant, water_scaling){

  lambda <- (dominant == 6)*water_scaling+(dominant != 6)
  return(lambda)

}

compute_proportion_class_lambda <- function(prop_water,water_lambda,land_lambda){

  water_lambda+(1-prop_water)*(land_lambda-water_lambda)

}


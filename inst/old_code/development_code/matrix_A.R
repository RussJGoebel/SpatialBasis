#' st_fast_intersection (stack overflow)
#'
#' @param x An sf object, first parameter in sf:st_intersects
#' @param y An sf object, second parameter in sf:st_intersects
#' @param ... other arguments for sf::st_intersects
#'
#' @return Geometry if the the intersection between x and subset of y
#' @export
#'
#' @examples
st_fast_intersection <- function(x,y,...){

  y_subset <-
    sf::st_intersects(x, y) %>%
    unlist() %>%
    unique() %>%
    sort() %>%
    {y[.,]}

  sf::st_intersection(x, y_subset, ...)
}

#' st_parallel_intersection
#'
#' Similar to st_intersection, but runs in parallel to make use of more processors.
#'
#'
#' @param soundings An sf object, The observed soundings
#' @param fine_grid An sf object, a fine grid on which we want to estimate SIF
#' @param mc.cores number of cores to using in mclapply
#' @param n_groups number core groups. defaults to mc.cores parameter
#'
#' @return sf dataframe of intersections.
#' @export
#'
#' @examples if(interactive()) intersections <- st_parallel_intersection(soundings, fine_grid,mc.cores = 2)
st_parallel_intersection <- function(soundings, fine_grid,
                                     mc.cores = parallel::detectCores(logical = FALSE),n_groups = mc.cores){

  n <- dim(fine_grid)[1]
  raster_pixels <- 1:n

  intersects <- sf::st_intersects(fine_grid,soundings)

  intersections <- parallel::mclapply(raster_pixels,function(x){sf::st_intersection(fine_grid[x,],
                                                                             soundings[intersects[[x]],])},
                                      mc.cores = mc.cores)

  intersections <- dplyr::bind_rows(intersections)

  return(intersections)

}
#' compute_A_matrix
#'
#' Computes the matrix of weights, denoted $A$, determining overlapping areas.
#'
#' @param fine_grid An sf object, a fine grid on which we want to estimate SIF
#' @param soundings An sf object, The observed soundings
#'
#' @return Matrix of wieghted area overlap between grid and soundings
#' @export
#'
#' @examples if(interactive()){APG <- compute_A_matrix(target_grid, soundings, mc.cores = 2)}
compute_A_matrix <- function(fine_grid,
                             soundings,
                             intersection = NULL,
                             mc.cores = parallel::detectCores(logical = FALSE))  {

  if(is.null(intersection)){
    intersection <- st_parallel_intersection(soundings,fine_grid, mc.cores)
  }

  intersects <- sf::st_intersects(soundings,fine_grid,sparse = FALSE)

  sounding_areas <- as.vector(sf::st_area(soundings))
  intersection_areas <- sf::st_area(intersection)

  intersects[c(intersects)] <- intersection_areas;
  A <- sweep(intersects,1,sounding_areas,FUN = "/")

  return(Matrix::Matrix(A))

}


#' compute_A_matrix_2d_gauss
#'
#' @param fine_grid An sf object, a fine grid on which we want to estimate SIF or XCO2
#' @param soundings An sf object, The observed soundings
#' @param major_axis_sigma_scaler a factor to scale the guassian blur in the major axis direction
#' @param minor_axis_sigma_scaler a factor to scale the gaussian blue in the minor axis direction
#'
#' @return Matrix of weights corresponding to 2d guassian centered on each sounding
#' @export
#'
#' @examples
#' # Note: if(interactive()){} does not run in checks and can temporarily wrap long-running examples.
#' if(interactive()){A2D <- compute_A_matrix_2d_gauss(target_grid, soundings)}
compute_A_matrix_2d_gauss <- function(fine_grid,soundings,
                                      major_axis_sigma_scaler=1,
                                      minor_axis_sigma_scaler=1){
  n_obs   <- length(soundings$geometry)
  n_state <- length(fine_grid$col)
  A <- matrix(0,nrow = n_obs, ncol = n_state)

  #Since we need to work in meters, find the best UTM zone and use that:
  x1 <- soundings$geometry[[1]][1][[1]][1,1]
  y1 <- soundings$geometry[[1]][1][[1]][1,2]
  utm_zone <- floor( (x1 + 180)*(60/360) )
  if( y1 > 0){
    epsg_crs <- utm_zone + 32601
  }else{
    epsg_crs <- utm_zome + 32701
  }
  proj_soundings <- sf::st_transform(soundings,epsg_crs)
  proj_grid <- sf::st_transform(fine_grid,epsg_crs)

  proj_grid$x0 <- 0
  proj_grid$y0 <- 0

  for (ii in 1:n_state){
    proj_grid$x0[ii] <-  ( ( proj_grid$x[[ii]][1][[1]][1,1] + proj_grid$x[[ii]][1][[1]][3,1] ) / 2 )
    proj_grid$y0[ii] <-  ( ( proj_grid$x[[ii]][1][[1]][1,2] + proj_grid$x[[ii]][1][[1]][3,2] ) / 2 )
  }

  for( x in 1:n_obs){
    x1 <- proj_soundings$geometry[[x]][1][[1]][1,1]
    x2 <- proj_soundings$geometry[[x]][1][[1]][2,1]
    x3 <- proj_soundings$geometry[[x]][1][[1]][3,1]
    y1 <- proj_soundings$geometry[[x]][1][[1]][1,2]
    y2 <- proj_soundings$geometry[[x]][1][[1]][2,2]
    y3 <- proj_soundings$geometry[[x]][1][[1]][3,2]

    minor_axis_length <- sqrt( (x2 - x1)^2 + (y2 - y1 )^2 ) / 2
    major_axis_length <- sqrt( (x3 - x2)^2 + (y3 - y2 )^2 ) / 2
    minor_axis_angle  <- atan2( (x2-x1), (y2-y1) )
    major_axis_angle  <- atan2( (y3-y2), (x3-x2) )

    sig_x <- 2*((major_axis_length*major_axis_sigma_scaler)^2)
    sig_y <- 2*((minor_axis_length*minor_axis_sigma_scaler)^2)

    x0 <- mean( c(x1,x3))
    y0 <- mean( c(y1,y3))
    proj_grid$delta_x <- proj_grid$x0 - x0
    proj_grid$delta_y <- proj_grid$y0 - y0

    proj_grid$d <- sqrt( proj_grid$delta_x^2 + proj_grid$delta_y^2 )
    proj_grid$angle <- atan2( proj_grid$delta_y, proj_grid$delta_x )

    proj_grid$dy <- proj_grid$d*sin(proj_grid$angle - major_axis_angle )
    proj_grid$dx <- proj_grid$d*cos(proj_grid$angle - major_axis_angle ) #Note: I'm setting it orthoginal here...probably shouldn't.

    proj_grid$gauss <- exp(-1*(proj_grid$dx^2)/(sig_x)) * exp(-1*(proj_grid$dy^2)/(sig_y))
    proj_grid$gauss <- proj_grid$gauss/sum(proj_grid$gauss)
    A[x,] <- proj_grid$gauss
  }
  A[ A<1e-6 ] <- 0
  return(Matrix::Matrix(A))
}

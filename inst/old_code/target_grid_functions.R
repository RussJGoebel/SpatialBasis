#' convert lat and lon to an Alber's Equal Area crs string
#'
#' @param lon_0 Longitude of projection center in degrees.
#' @param lat_0 Latitude of projection center in degrees
#' @param lat_1 First standard parallel.
#' @param lat_2 Second standard parallel.
#'
#' @return String appropriate for crs.
#' @export
#'
#' @examples crs_string <- convert_lat_and_lon_to_aea_crs_string(53.2, 19.5, 19.7, 53.4)
convert_lat_and_lon_to_aea_crs_string <- function(lon_0,lat_0,lat_1,lat_2){

  crs_string <- paste("+proj=aea ",
                      "+lon_0=",lon_0,
                      " +lat_0=",lat_0,
                      " +lat_1=",lat_1,
                      " +lat_2=",lat_2,
                      " +units=m",
                      sep = "")

}

#' create a crs string corresponding to Alber's Equal Area projection centered
#' at the centroid of soundings
#'
#' @param soundings An sf object (thought to be the original soundings)
#'
#' @return String appropriate for crs caluclated from soundings centroid
#' @export
#'
#' @examples crs_string <- create_aea_projection_crs_string_at_centroid(soundings, lat_1 = 29.5,lat_2 = 45.5)
create_aea_projection_crs_string_at_centroid <- function(soundings,
                                                         lat_1 = 29.5,
                                                         lat_2 = 45.5){

  sounding_centroid <- sf::st_centroid(sf::st_union(soundings))
  sounding_centroid_coords <- sf::st_coordinates(sounding_centroid)
  lon_0 <-  sounding_centroid_coords[1]
  lat_0 <- sounding_centroid_coords[2]

  crs_string <- convert_lat_and_lon_to_aea_crs_string(lon_0,lat_0,lat_1,lat_2)

  return(crs_string)

}

#' Title
#'
#' @param soundings An sf object (thought to be the original soundings)
#' @param grid An sf object containing the underlying fine grid.
#' @param buffer A numeric value specifying a buffer to put around the soundings (meters)
#'
#' @return original grid with added cells to increase bounding box beyond soundings bounding box
#' @export
#'
#' @examples buffered_grid <- buffered_convex_hull(soundings, target_grid, buffer=500)
buffered_convex_hull <- function(soundings, grid, buffer) {
    convex_hull <- sf::st_convex_hull(sf::st_union(soundings))

    within_distance_of_soundings <- sf::st_is_within_distance(grid, convex_hull, buffer)
    within_distance_of_soundings <- sapply(within_distance_of_soundings,
                                           FUN=function(x){length(x) > 0})

    grid <- grid[within_distance_of_soundings,]
    return (sf::st_as_sf(grid))
}

#' generate grid
#'
#' Uses the centroid of the union of the soundings as the projection center.
#'
#' @param soundings An sf object (thought to be the original soundings)
#' @param res A numeric value specifying the resolution (in meters)
#' @param buffer A numeric value specifying a buffer to put around the soundings
#' @param grid_shape One of "buffered_convex_hull" or "rectangle". If "buffered convex hull", then the the shape of the grid is a convex hull of the soundings, buffered by a distance of 'buffer'.
#' @param lon_0 Longitude of projection center.
#' @param lat_0 Latitude of projection center?
#' @param lat_1 First standard parallel
#' @param lat_2 Second standard parallel
#'
#' @return target_grid with cells of size res used for downscaling soundings
#' @export
#'
#' @examples grid <- generate_grid(soundings, res=330, buffer=2*330)
generate_grid <- function(soundings,
                          res,
                          buffer = 2*res,
                          grid_shape = c("buffered_convex_hull","rectangle"),
                          lon_0 = NULL,
                          lat_0 = NULL,
                          lat_1 = 29.5,
                          lat_2 = 45.5){

  grid_shape <- match.arg(grid_shape)

  # define crs -----------------------------------------------------------------

  if(is.null(lon_0) | is.null(lat_0)){
    message("Projection centered using soundings.")
    crs_string <- create_aea_projection_crs_string_at_centroid(soundings,
                                                               lat_1, lat_2)
  } else {
    message("Projection centered using (lon_0, lat_0)")
    crs_string <- convert_lat_and_lon_to_aea_crs_string(lon_0,lat_0,lat_1,lat_2)
  }

  message("Crs:", crs_string)

  soundings_proj <- sf::st_transform(soundings,crs = crs_string)
  crs <- sf::st_crs(soundings)

  # apply the buffer -----------------------------------------------------------

  buffered_soundings <- sf::st_buffer(sf::st_union(soundings_proj), dist = buffer)

  # Generate the grid within the bounding box with the calculated cell size ----
  grid <- sf::st_make_grid(buffered_soundings, cellsize = res, square = TRUE)

  n_rows <- diff(sf::st_bbox(grid)[c(2,4)])/res
  n_cols <- diff(sf::st_bbox(grid)[c(1,3)])/res

  cat(paste("Dimensions of grid before is",n_rows, "by ",n_cols,"\n"))

  grid <- sf::st_as_sf(grid)
  grid <- sf::st_transform(grid,crs)

  grid$row <- rep(1:n_rows,each = n_cols)
  grid$col <- rep(1:n_cols,times = n_rows)
  grid$geometry <- grid$x

  # Subset only a convex hull --------------------------------------------------

  if(grid_shape == "buffered_convex_hull"){
    message("Grid shape is a convex hull of the soundings with buffer applied.")
    grid <- buffered_convex_hull(soundings = soundings, grid=grid, buffer=buffer)
  }

  return(grid)
}

#' aggregate_SpatRaster_over_polygon
#'
#' @param SpatRaster_object SpatialRaster of landcover, albedo, or other dataset
#' @param sf_polygon sf object, usually a fine grid to resample SpatRaster_object
#'
#' @returns
#' @export
#'
#' @examples
aggregate_SpatRaster_over_polygon <- function(SpatRaster_object,sf_polygon){

  message("Aggregating using mean.")
  rasterized_sf <- sf_to_SpatRaster(sf_polygon,SpatRaster_object)
  sf_coordinates <- which(!is.na(terra::values(rasterized_sf)))

  sf_polygon$mean_value <- mean(terra::values(SpatRaster_object)[sf_coordinates])

  return(sf_polygon)

}

#' summarize_SpatRaster_class_representation_over_polygon
#' If a SpatRaster object's values uniquely correspond to classes (e.g, landcover class),
#' Aggregatete
#'
#' @param sf_grid sf object, usually a fine grid to resample SpatRaster_object
#' @param SpatRaster_object SpatialRaster of landcover, albedo, or other dataset
#'
#' @return A tibble whose first column is an index
#' @export
#'
#' @examples
summarize_SpatRaster_class_representation_over_grid <- function(SpatRaster_object,sf_grid){

  sf_grid$pixel <- 1:dim(sf_grid)[1] # Indexing each pixel
  message("Converting sf object into SpatRaster...")
  grid_rast <- terra::rasterize(sf_grid,SpatRaster_object,field = "pixel")

  message("Creating tibble to summarize classes...")
  t <- dplyr::tibble(pixel = terra::values(grid_rast),label = terra::values(SpatRaster_object))

  missing_pixels <- !(sf_grid$pixel %in% t$pixel)
  number_of_missing_pixels <- sum(missing_pixels)

  if(number_of_missing_pixels > 0){
    message(number_of_missing_pixels," pixels in the sf object don't overlap with the SpatRaster object, and will be discarded.")
    sf_grid <- sf_grid[!missing_pixels,]
  }

  message("Summarizing classes by count...")
  t <- dplyr::count(t,pixel,label)
  message("Summarizing classes by proportion...")
  t <- dplyr::group_by(t,pixel)
  t <- dplyr::mutate(t,proportion = n/sum(n))

  message("Computing dominant class...")
  dominant_class <- dplyr::summarize(t,label[which.max(proportion)])

  message("Adding everything to original sf object.")
  t <- tidyr::pivot_wider(t,names_from = label,values_from = c(n,proportion))

  sf_grid <- dplyr::full_join(sf_grid,t,by = "pixel")
  sf_grid <- dplyr::mutate(sf_grid,dominant_class = dominant_class$`label[which.max(proportion)]`)

  return(sf_grid)

}

#' summarize_SpatRaster_class_representation_over_polygon
#' If a SpatRaster object's values uniquely correspond to classes (e.g, landcover class),
#' Aggregatete
#'
#' @param sf_grid sf object, usually a fine grid to resample SpatRaster_object
#' @param SpatRaster_object SpatialRaster of landcover, albedo, or other dataset
#'
#' @return A tibble whose first column is an index
#' @export
#'
#' @examples
summarize_SpatRaster_mean_over_grid <- function(SpatRaster_object,sf_grid){

  sf_grid$pixel <- 1:dim(sf_grid)[1] # Indexing each pixel
  message("Converting sf object into SpatRaster...")
  grid_rast <- terra::rasterize(sf_grid,SpatRaster_object,field = "pixel")

  message("Creating tibble to summarize values...")
  t <- dplyr::tibble(pixel = terra::values(grid_rast),value = terra::values(SpatRaster_object))

  missing_pixels <- !(sf_grid$pixel %in% t$pixel)
  number_of_missing_pixels <- sum(missing_pixels)

  if(number_of_missing_pixels > 0){
    message(number_of_missing_pixels," pixels in the sf object don't overlap with the SpatRaster object, and will be discarded.")
    sf_grid <- sf_grid[!missing_pixels,]
  }

  message("Summarizing values by mean...")
  t <- dplyr::group_by(t,pixel)
  t <- dplyr::summarize(t,mean_value = mean(value, na.rm = T))

  sf_grid <- dplyr::full_join(sf_grid,t,by = "pixel")
 # sf_grid <- dplyr::mutate(sf_grid,dominant_class = dominant_class$`label[which.max(proportion)]`)

  return(sf_grid)

}

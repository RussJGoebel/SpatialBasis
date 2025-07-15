#' SpatRaster_to_sf
#'
#' Converts SpatRaster objects, such as those obtained from reading .tif files using
#' terra::rast(), into sf objects.
#'
#' @param SpatRaster SpatRaster to contert to polygons with terra
#' @param ... optional parameters for terra::as.polygons
#'
#' @return
#' @export
#'
#' @examples
SpatRaster_to_sf <- function(SpatRaster,...){

  SpatRaster_as_polygon <- terra::as.polygons(SpatRaster,...)
  sf_object <- sf::st_as_sf(SpatRaster_as_polygon)

  return(sf_object)

}

#' sf_to_SpatRaster
#'
#' Converts an sf object from the sf package into a SpatRaster object from the terra package
#'
#' @param sf_object
#' @param SpatRaster_object
#' @param field
#' @param res
#' @param ext
#'
#' @return
#' @export
#'
#' @examples
sf_to_SpatRaster <- function(sf_object,
                             SpatRaster_object = NULL,
                             field = "",
                             res = NULL,
                             ext = NULL){

  if(!is.null(SpatRaster_object)){

    ext <- terra::ext(SpatRaster_object)
    res <- terra::res(SpatRaster_object)
    r <- terra::rast(ext,resolution = res)

    vectorized_sf <- terra::vect(sf_object)
    rasterized_sf <- terra::rasterize(vectorized_sf,r,field = field,)

  }else{

    if(is.null(res) | is.null(ext)) stop("Either a SpatRaster object or both desired ext and res are required.")

    r <- terra::rast(ext,resolution = res)

    vectorized_sf <- terra::vect(sf_object)
    rasterized_sf <- terra::rasterize(vectorized_sf,r,field = field)

  }

  return(rasterized_sf)

}

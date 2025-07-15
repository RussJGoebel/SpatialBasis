#' Return sf object with both original and complement geometries
#'
#' @param x An sf object (e.g., water polygons)
#' @param bbox Optional bounding box (defaults to bbox of x)
#' @param type_names Character vector of length 2, e.g. c("water", "nonwater")
#' @return An sf object with a 'type' column distinguishing original and complement
#' @export
sf_with_complement <- function(x, bbox = NULL, type_names = c("original", "complement")) {
  stopifnot(inherits(x, "sf"))
  stopifnot(length(type_names) == 2)

  # Compute bounding box
  if (is.null(bbox)) {
    bbox <- sf::st_bbox(x)
  }

  bbox_poly <- sf::st_as_sfc(bbox)
  bbox_poly <- sf::st_set_crs(bbox_poly, sf::st_crs(x))

  # Union and validate original
  x_union <- sf::st_union(x)
  x_union <- sf::st_make_valid(x_union)

  # Compute complement
  complement_geom <- sf::st_difference(bbox_poly, x_union)
  complement_geom <- sf::st_make_valid(complement_geom)

  # Construct sf objects with 'type' column
  original_sf <- sf::st_sf(type = type_names[1], geometry = x_union)
  complement_sf <- sf::st_sf(type = type_names[2], geometry = complement_geom)

  # Combine into single sf
  combined <- rbind(original_sf, complement_sf)
  return(combined)
}

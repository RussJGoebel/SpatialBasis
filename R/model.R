
#' add_guassian_to_sf
#'
#' @param sf_obj An sf object with at least one numeric column.
#' @param lat_0 gaussian plume source latitude in degrees
#' @param lon_0 guassian plume source longitude in degrees
#' @param ws_mps local surface windspeed in meters per second
#' @param wd_deg local surface wind direction in degrees
#' @param x_0_m guassian plume length scale in meters (default 1 km)
#'
#' @return original sf_obj with guassian plume model eastimate at each geometry
#' @export
#'
#' @examples grid <- add_gaussian_to_sf(target_grid, lat_0=19.5, lon_0=53.3, ws_mps=2.5, wd_deg =270, x_0_m=1000)
add_gaussian_to_sf <- function(
  sf_obj,
  lat_0, lon_0, ws_mps, wd_deg,
  x_0_m = 1000) {

  source_lonlat <- c(lon_0, lat_0) # lon, lat
  wd_rad <- wd_deg * ( pi / 180 )  # rad

  gauss_cells <- c()
  r_cells <- c()
  x_cells <- c()
  y_cells <- c()
  for (pgon in sf_obj$geometry) {
    pgon_point <- sf::st_centroid(pgon)
    pgon_lonlat <- c(pgon_point[1], pgon_point[2])     # lon, lat

    theta_z <- geosphere::bearing(source_lonlat, pgon_lonlat) * ( pi / 180 ) # rad
    r       <- geosphere::distHaversine(source_lonlat, pgon_lonlat) # m

    dyg <- r * sin(theta_z - wd_rad)  #m
    dxg <- r * cos(theta_z - wd_rad)  #m

    a <- 213 # moenini etal 2024 varies with distance (a = 104, 213)
    X <- (dxg / x_0_m)
    sigma_y <- abs(a * (X^0.894))

    gauss <- (1/ ( sqrt(2 * pi) * sigma_y * ws_mps)) * exp(-(1/2)*((dyg/sigma_y)^2))

    gauss_cells <- c(gauss_cells, gauss)
    r_cells <- c(r_cells, r)
    x_cells <- c(x_cells, dxg)
    y_cells <- c(y_cells, dyg)
  }
  sf_obj$gaussian_gpm2 <- gauss_cells
  sf_obj$src_dist_m <- r_cells
  sf_obj$x_dist_m <- x_cells
  sf_obj$y_dist_m <- y_cells

  return (sf_obj)
}


#' custom_gaussian_mask
#'
#' @param sf_obj An sf object containing guassian_gpm2 estimate from add_gaussain_to_sf().
#' @param source_radius_m radius around source to remove from in-plume mask
#' @param inner_radius_m radius around source to remove from background mask
#' @param outer_radius_m radius around source containing valid samples
#' @param gaussian_cutoff_percent percent of plume maximum value used as lower limit boundary for in-plume mask
#' @param gaussian_buffer_percent percent of plume maximum value used for buffer mask. Must be lower than gaussian_cutoff_percent.
#' @param nucleous_cutoff_percent percent of plume maximum value used as upper limit bouardary for in-plume mask (Not Recommended)
#'
#' @return input sf_obj with boolean columns for in-plume, background, and buffer regions of guassian model.
#' @export
#'
#' @examples # grid <- custom_gaussian_mask(grid)
custom_gaussian_mask <- function(sf_obj,
                                 inner_radius_m =  3000,
                                 outer_radius_m = 30000,
                                 source_radius_m =1000,
                                 gaussian_cutoff_percent = 1,
                                 gaussian_buffer_percent = 0.1,
                                 nucleous_cutoff_percent = NULL) {

  if (is.null(source_radius_m)) { source_radius_m <- 0 }

  gaussian_cutoff_gpm2 <- max(sf_obj$gaussian_gpm2, na.rm = T) * gaussian_cutoff_percent/100
  buffer_cutoff_gpm2   <- max(sf_obj$gaussian_gpm2, na.rm = T) * gaussian_buffer_percent/100

  above_gauss_cutoff_mask  <- (sf_obj$gaussian_gpm2 > gaussian_cutoff_gpm2)
  above_gauss_buffer_mask  <- (sf_obj$gaussian_gpm2 > buffer_cutoff_gpm2)
  below_gauss_buffer_mask  <- !above_gauss_buffer_mask | is.na(sf_obj$gaussian_gpm2)
  #inside_gauss_buffer_mask <- above_gauss_buffer_mask & !above_gauss_cutoff_mask

  valid_outer_distance_mask <- (sf_obj$src_dist_m < outer_radius_m)
  if ("buffer" %in% colnames(sf_obj)) {
   valid_outer_distance_mask <- valid_outer_distance_mask & sf_obj$buffer
  }

  valid_plume_mask      <- above_gauss_cutoff_mask & (valid_outer_distance_mask & (sf_obj$src_dist_m > source_radius_m ))
  valid_plume_mask[(sapply(valid_plume_mask, is.na))] <- F # impotant!

  valid_background_mask <- below_gauss_buffer_mask & (valid_outer_distance_mask & (sf_obj$src_dist_m > inner_radius_m ))
  valid_other_mask      <- !valid_plume_mask & !valid_background_mask


  sf_obj$plume <- valid_plume_mask
  sf_obj$background <- valid_background_mask
  sf_obj$other <- valid_other_mask

  return (sf_obj)
}

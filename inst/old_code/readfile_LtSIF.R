#' read_oco2_v11_sif_lite_nowarn
#'
#' @param filename oco2 SIF Lite file (*.nc4)
#' @param target_name string name for target
#'
#' @return sf dataframe of SIF lite file. output is redirected to nullfile() to hide encoding warning
#' @export
#'
#' @examples data <- read_oco2_v11_sif_lite_nowarn("oco2_LtSIF_210424_B11012Ar_220804024105s.nc4", "bostonMA")
read_oco2_v11_sif_lite_nowarn <- function(filename, target_name) {
  # does this ignore errors too?
  capture.output(
    df <- read_oco2_v11_sif_lite(filename, target_name),
    file = nullfile())
  return (df)
}

#' read_oco2_v11_sif_lite
#'
#' @param filename oco2 SIF Lite file (*.nc4)
#' @param target_name string name for target
#'
#' @return sf dataframe of SIF lite file.
#' @export
#'
#' @examples data <- read_oco2_v11_sif_lite("oco2_LtSIF_210424_B11012Ar_220804024105s.nc4", "bostonMA")
read_oco2_v11_sif_lite <- function(filename, target_name) {

  nc <- ncdf4::nc_open(filename)

  tstamps    <- ncdf4::ncvar_get(nc,"Delta_Time")
  lat_cnrs   <- ncdf4::ncvar_get(nc,"Latitude_Corners")
  lng_cnrs   <- ncdf4::ncvar_get(nc,"Longitude_Corners")
  meas_mode  <- ncdf4::ncvar_get(nc,"Metadata/MeasurementMode")
  cloud_flag <- ncdf4::ncvar_get(nc,"Cloud/cloud_flag_abp")
  liteq_flag  <- ncdf4::ncvar_get(nc,"Quality_Flag")
  sif757        <- ncdf4::ncvar_get(nc,"Science/SIF_757nm")
  sif757_sigma  <- ncdf4::ncvar_get(nc,"Science/SIF_Uncertainty_757nm")
  sif771        <- ncdf4::ncvar_get(nc,"Science/SIF_771nm")
  sif771_sigma  <- ncdf4::ncvar_get(nc,"Science/SIF_Uncertainty_771nm")
  oco_landf  <- ncdf4::ncvar_get(nc,"Science/sounding_land_fraction")
  #sound_flag <- ncvar_get(nc,"Science/sounding_qual_flag")
  view_zenith <- ncdf4::ncvar_get(nc,"Geolocation/sensor_zenith_angle")
  view_azimuth <- ncdf4::ncvar_get(nc,"Geolocation/sensor_azimuth_angle")
  surf_albedo <- ncdf4::ncvar_get(nc,"Cloud/surface_albedo_abp")

  ncdf4::nc_close(nc)

  df <- data.frame(Timestamp= as.POSIXct(tstamps, origin = "1990-01-01" ),
                   MeasurementMode=meas_mode,
                   ViewZenith = view_zenith,
                   ViewAzimuth = view_azimuth,
                   SIF_757nm=sif757,
                   SIF_Uncertainty_757nm  = sif757_sigma,
                   SIF_771nm=sif771,
                   SIF_Uncertainty_771nm  = sif771_sigma,
                   Albedo = surf_albedo,
                   cloud_flag_abp = cloud_flag,
                   Quality_Flag = liteq_flag,
                   sounding_land_fraction  = oco_landf)
  #sounding_qual_flag = sound_flag)

  downscaling::get_target_bounds(target_name)

  lat_cnrs[is.na(lat_cnrs)] = 0
  lng_cnrs[is.na(lng_cnrs)] = 0
  n_soundings <- dim(lat_cnrs)[2]

  sounding_mask <- c()
  poly_list_soundings <- list()
  for(i in 1:n_soundings){
    if( (lat_cnrs[1,i] > lat_min ) & ( lat_cnrs[1,i] < lat_max ) & (lng_cnrs[1,i] > lon_min) &(lng_cnrs[1,i] < lon_max)){
      if(meas_mode[i] %in% c(0,2,4) ){
        pgon <- cbind( lng_cnrs[,i] , lat_cnrs[,i] )
        pgon <- sf::st_polygon(list( rbind( pgon, pgon[1,]) ) )
        if(sf::st_area(pgon) > 1e-6){
          poly_list_soundings <- append( poly_list_soundings, list( pgon) )
          sounding_mask <- rbind(sounding_mask, i)
        }
      }
    }
  }

  df <- df[sounding_mask, ]

  n <- length(df$SIF_757nm)

  df <- cbind(data.frame(Soundings=sounding_mask), df)
  sf_polys <- sf::st_as_sf(sf::st_sfc(poly_list_soundings))
  df$geometry <- sf_polys$x
  rownames(df) <- 1:nrow(df)

  sounding_df <- df
  #sounding_df <- df[df$MeasurementMode == 2,]
  #sounding_df <- sounding_df[(sounding_df$Quality_Flag != 0) | (sounding_df$Quality_Flag == 1),]

  soundings <- sf::st_as_sf(sounding_df, crs =  sf::st_crs(4326))
  sum(!sf::st_is_valid(soundings))
  soundings <- soundings[sf::st_is_valid(soundings),]

  return (soundings)
}

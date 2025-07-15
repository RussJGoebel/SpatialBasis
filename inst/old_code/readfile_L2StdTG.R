#' read_oco2_v11_co2_L2StdTG_nowarn
#'
#' @param filename oco2 L2 Standard file (*.h5)
#' @param target_name string name for target
#'
#' @return sf dataframe of L2 standard file. output is redirected to nullfile() to hide encoding warning
#' @export
#'
#' @examples oco2_L2StdTG_210513 <- read_oco2_v11_co2_L2StdTG_nowarn("oco2_L2StdTG_36511a_210513_B11006r_220522083538.h5", "belchatowPl")
read_oco2_v11_co2_L2StdTG_nowarn <- function(filename, target_name) {
  # does this ignore errors too?
  capture.output(
    df <- read_oco2_v11_co2_L2StdTG(filename, target_name),
    file = nullfile())
  return (df)
}

#' read_oco2_v11_co2_L2StdTG
#'
#' @param filename oco2 L2 Standard file (*.h5)
#' @param target_name string name for target
#'
#' @return sf dataframe of L2 standard file.
#' @export
#'
#' @examples oco2_L2StdTG_210513 <- read_oco2_v11_co2_L2StdTG("oco2_L2StdTG_36511a_210513_B11006r_220522083538.h5", "belchatowPl")
read_oco2_v11_co2_L2StdTG <- function(filename, target_name) {

  nc <- ncdf4::nc_open(filename)

  tstamps      <- ncdf4::ncvar_get(nc,"RetrievalHeader/retrieval_time_tai93")
  lat_cnrs     <- ncdf4::ncvar_get(nc,"RetrievalGeometry/retrieval_vertex_latitude")[,1,]
  lng_cnrs     <- ncdf4::ncvar_get(nc,"RetrievalGeometry/retrieval_vertex_longitude")[,1,]
  meas_mode    <- ncdf4::ncvar_get(nc,"Metadata/OperationMode")
  cloud_flag   <- ncdf4::ncvar_get(nc,"PreprocessingResults/cloud_flag_abp")
  oco_landf    <- ncdf4::ncvar_get(nc,"RetrievalGeometry/retrieval_land_fraction")
  view_zenith  <- ncdf4::ncvar_get(nc,"RetrievalGeometry/retrieval_solar_zenith")
  view_azimuth <- ncdf4::ncvar_get(nc,"RetrievalGeometry/retrieval_solar_azimuth")
  surf_albedo  <- ncdf4::ncvar_get(nc,"PreprocessingResults/albedo_o2_abp")
  xco2         <- ncdf4::ncvar_get(nc,"RetrievalResults/xco2") #moles/mole
  xco2_sigma   <- ncdf4::ncvar_get(nc,"RetrievalResults/xco2_uncert") # moles/mole
  pressure       <- ncdf4::ncvar_get(nc,"RetrievalResults/surface_pressure_fph") # Pa
  pressure_sigma <- ncdf4::ncvar_get(nc,"RetrievalResults/surface_pressure_uncert_fph") # Pa
  co2_col      <- ncdf4::ncvar_get(nc,"RetrievalResults/retrieved_co2_column") # molecules/m^2

  u_met        <- ncdf4::ncvar_get(nc,"RetrievalResults/wind_speed_u_met") # m/s
  v_met        <- ncdf4::ncvar_get(nc,"RetrievalResults/wind_speed_v_met") # m/s

  ncdf4::nc_close(nc)

  df <- data.frame(Timestamp=as.POSIXct(tstamps, origin="1990-01-01"),
                   MeasurementMode=meas_mode,
                   ViewZenith = view_zenith,
                   ViewAzimuth = view_azimuth,
                   XCO2=xco2,                               # moles/mole
                   XCO2_Uncertainty = xco2_sigma,           # moles/mole
                   XCO2_ppm=xco2 * 1E6,                     # ppm
                   XCO2_Uncertainty_ppm = xco2_sigma * 1E6, # ppm
                   Pressure_pascals = pressure,             # Pa
                   Pressure_Uncertainty_pascals = pressure, # Pa
                   CO2_Column_moleculespm2= co2_col,        # Molecules(co2)/m^2
                   Wind_U_mps = u_met,                      # m/s
                   Wind_V_mps = v_met,                      # m/s
                   #Albedo = surf_albedo,
                   Cloud_Flag = cloud_flag,
                   Land_Fraction  = oco_landf
                   )


  downscaling::get_target_bounds(target_name)

  lat_cnrs[is.na(lat_cnrs)] = 0
  lng_cnrs[is.na(lng_cnrs)] = 0
  n_soundings <- dim(lat_cnrs)[2]

  sounding_mask <- c()
  poly_list_soundings <- list()
  for(i in 1:n_soundings){
    #if( (lat_cnrs[1,i] > lat_min ) & ( lat_cnrs[1,i] < lat_max ) & (lng_cnrs[1,i] > lon_min) &(lng_cnrs[1,i] < lon_max)){
    #if(meas_mode[i] %in% c("TG") ){
    pgon <- cbind( lng_cnrs[,i] , lat_cnrs[,i] )
    pgon <- sf::st_polygon(list( rbind( pgon, pgon[1,]) ) )
    if(sf::st_area(pgon) > 1e-6){
      poly_list_soundings <- append( poly_list_soundings, list( pgon) )
      sounding_mask <- rbind(sounding_mask, i)
    }
    #}
    #}
  }

  df <- df[sounding_mask, ]

  n <- length(df$XCO2)

  df <- cbind(data.frame(Soundings=sounding_mask), df)
  sf_polys <- sf::st_as_sf(sf::st_sfc(poly_list_soundings))
  df$geometry <- sf_polys$x
  rownames(df) <- 1:nrow(df)

  sounding_df <- df

  soundings <- sf::st_as_sf(sounding_df, crs =  sf::st_crs(4326))
  sum(!sf::st_is_valid(soundings))
  soundings <- soundings[sf::st_is_valid(soundings),]

  return (soundings)
}

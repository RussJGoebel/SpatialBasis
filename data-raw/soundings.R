## Code to prepare `soundings` dataset #########################################
# This code is written by Leeza Moldavchuk
# Modified by Russell Goebel to use the larger bounding box & produce valid sf objects
# Most recently updated 09/26/24 as an example dataset for downscaling package
################################################################################

require(ncdf4)     # OCO data files are .nc
require(raster)    # needed for rasters
library(leaflet)   # Maps
library(sf)        # Simple Features library for geometries
library(ggplot2)   # plots
library(geosphere)
library(exactextractr)
library(terra)
library(mapview)
library(htmltools)
library(RColorBrewer)
library(Matrix)
library(viridis)
#library(gdalUtils)
library(tidyr)
source("/projectnb/buultra/tjones/legend_swap.r")


target <- "Boston"
#target <- "Toronto"
#target <- "Lamont"

# Boston Summer
filename <- "/projectnb/buultra/oco_data/oco2_sif/v11/2022/oco2_LtSIF_220614_B11012Ar_220910020814s.nc4"

# Lamont Winter 12/12/2018
#filename <- "/projectnb/buultra/oco_data/oco2_sif/v11/2018/oco2_LtSIF_181212_B11012Ar_221017173252s.nc4"

# Lamont Summer 07/31/2019
#filename <- "/projectnb/buultra/oco_data/oco2_sif/v11/2019/oco2_LtSIF_190731_B11012Ar_220928175944s.nc4"

# Toronto Winter 10/05/2022
#filename <- '/projectnb/buultra/oco_data/oco2_sif/v11/2022/oco2_LtSIF_221005_B11012Ar_221205192438s.nc4'

# Toronto Summer 06/24/2022
#filename <- "/projectnb/buultra/oco_data/oco2_sif/v11/2022/oco2_LtSIF_220624_B11012Ar_220910021248s.nc4"

df <- data.frame()
poly_list_soundings <- list()

if(target == "Boston"){
  # This is the *larger* bounding box.
  lat_min <- 42.2207
  lat_max <- 42.53601
  lon_min <- -71.26929
  lon_max <- -70.95319
} else if(target == "Lamont"){
  lat_min <- 36.24495
  lat_max <- 36.8359
  lon_min <- -97.8228
  lon_max <- -97.15583
} else if(target == "Toronto"){
  lat_min <- 43.52777
  lat_max <- 43.98577
  lon_min <- -79.69979
  lon_max <- -79.19504
}


#import the OCO-2 data
nc <- nc_open(filename)

tstamps    <- ncvar_get(nc,"Delta_Time")
lat_cnrs   <- ncvar_get(nc,"Latitude_Corners")
lng_cnrs   <- ncvar_get(nc,"Longitude_Corners")
meas_mode  <- ncvar_get(nc,"Metadata/MeasurementMode")
cloud_flag <- ncvar_get(nc,"Cloud/cloud_flag_abp")
liteq_flag  <- ncvar_get(nc,"Quality_Flag")
sif757        <- ncvar_get(nc,"Science/SIF_757nm")
sif757_sigma  <- ncvar_get(nc,"Science/SIF_Uncertainty_757nm")
sif771        <- ncvar_get(nc,"Science/SIF_771nm")
sif771_sigma  <- ncvar_get(nc,"Science/SIF_Uncertainty_771nm")
oco_landf  <- ncvar_get(nc,"Science/sounding_land_fraction")
sound_flag <- ncvar_get(nc,"Science/sounding_qual_flag")
view_zenith <- ncvar_get(nc,"Geolocation/sensor_zenith_angle")
view_azimuth <- ncvar_get(nc,"Geolocation/sensor_azimuth_angle")
surf_albedo <- ncvar_get(nc,"Cloud/surface_albedo_abp")
lat_cnrs[is.na(lat_cnrs)] = 0
lng_cnrs[is.na(lng_cnrs)] = 0
n_soundings <- dim(lat_cnrs)[2]


for(i in 1:n_soundings){
  if( (lat_cnrs[1,i] > lat_min ) & ( lat_cnrs[1,i] < lat_max ) & (lng_cnrs[1,i] > lon_min) &(lng_cnrs[1,i] < lon_max)){
    if(meas_mode[i] %in% c(0,2,4) ){
      pgon <- cbind( lng_cnrs[,i] , lat_cnrs[,i] )
      pgon <- st_polygon(list( rbind( pgon, pgon[1,]) ) )
      if(st_area(pgon) > 1e-6){
        poly_list_soundings <- append( poly_list_soundings, list( pgon) )

        df <- rbind( df, data.frame(Sounding=i,
                                    Timestamp= as.POSIXct(tstamps[i], origin = "1990-01-01" ),
                                    MeasurementMode=meas_mode[i],
                                    ViewZenith = view_zenith[i],
                                    ViewAzimuth = view_azimuth[i],
                                    SIF_757nm=sif757[i],
                                    SIF_Uncertainty_757nm  = sif757_sigma[i],
                                    SIF_771nm=sif771[i],
                                    SIF_Uncertainty_771nm  = sif771_sigma[i],
                                    Albedo = surf_albedo[i],
                                    cloud_flag_abp = cloud_flag[i],
                                    Quality_Flag = liteq_flag[i],
                                    sounding_land_fraction  = oco_landf[i],
                                    sounding_qual_flag = sound_flag[i]) )
      }
    }
  }
}

doi <- df$Timestamp %>%
  mean %>%
  as.POSIXct( origin = "1990-01-01" ) %>%
  strftime( format = "%Y%m%d")

n <- length(df$SIF_757nm)


sf_polys <- st_as_sf(st_sfc(poly_list_soundings))
df$geometry <- sf_polys$x

date <- paste0("20",substring(filename,nchar(filename) - 32,nchar(filename) - 27))
#st_write(df, paste0("/projectnb/buultra/SIF_downscaling/Boston_sample/Boston_sample_",date,"_v2.geojson"))
#st_write(df, paste0("/projectnb/buultra/SIF_downscaling/target_data_samples/",target,"_sample_",date,".geojson"))


sounding_df <- data.frame()
sounding_df <- df

pal2 <- colorNumeric(c( '#FEFBE9', '#FCF7D5', '#F5F3C1', '#EAF0B5', '#DDECBF', '#D0E7CA', '#C2E3D2', '#B5DDD8', '#A8D8DC', '#9BD2E1', '#8DCBE4', '#81C4E7', '#7BBCE7', '#7EB2E4', '#88A5DD', '#9398D2', '#9B8AC4', '#9D7DB2', '#9A709E', '#906388', '#805770', '#684957', '#46353A'),
                     c(-.5,2),
                     na.color="transparent")


tag.map.title <- tags$style(HTML("
  .leaflet-control.map-title {
    transform: translate(-50%,20%);
    left: 50%;
    text-align: center;
    padding-left: 10px;
    padding-right: 10px;
    background: rgba(255,255,255,0.9);
    font-weight: bold;
    font-size: 20px;
  }
"))

title <- tags$div(
  tag.map.title, HTML(as.character(paste0(target, "  ",as.Date(date, format = "%Y%m%d"))))
)


sounding_bbox <- st_bbox(c(xmin=lon_min,xmax=lon_max,ymin=lat_min,ymax=lat_max), crs = st_crs(4326))


#plot the polygons:
m <- leaflet() %>% addProviderTiles(providers$CartoDB.Positron)
m <- m %>%
  addPolygons(data=st_sfc( poly_list_soundings, crs =  st_crs(4326)),
              weight = .04,
              color = "grey",
              fillOpacity = .3,
              fillColor = pal2(sounding_df$SIF_757nm),
              group = "Target") %>%
  #addRectangles(sounding_bbox[1],sounding_bbox[2],sounding_bbox[3],sounding_bbox[4], fill = F, color = "red") %>%
  addLegend_decreasing("topright",
                       pal =pal2,
                       values = c(-.5,2),
                       opacity =  0.8,
                       decreasing = T,
                       title = "SIF<sub>757</sub> [mW m<sup>-2</sup>nm<sup>-1</sup>sr<sup>-1</sup>]",
                       group="All") %>%
  addScaleBar(position="bottomright")%>%
  addControl(title, position = "bottomleft", className="map-title")

print(m)

### Removing invalid soundings #################################################
soundings <- sf::st_as_sf(sounding_df, crs =  st_crs(4326))
# some of the soundings are invalid:
sum(!sf::st_is_valid(soundings))

# these are basically just soundings that don't form complete polygons so we remove them:
ggplot(data = soundings[!sf::st_is_valid(soundings),])+geom_sf()

soundings <- soundings[sf::st_is_valid(soundings),]

#################################################################################

soundings <- downscaling::new_soundings(soundings,"SIF_757nm")

usethis::use_data(soundings, overwrite = TRUE)

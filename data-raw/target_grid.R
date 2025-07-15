## code to prepare `target_grid` dataset goes here

### target_grid consists of a regular lattice located over Boston.
# this code prepares target_grid by:
# 1) generating a grid of cells using Alber's Equal Area projection centered at
# the centroid of the union of the soundings of interest in Boston. (the 'soundings')
# dataset.)
# 2) augmenting the grid using median albedo values gathered from Google Earth Engine
# 3) augmenting the grid using landcover values gathered from Google Earth engine.
# specificaly, the I95 worldcover data.
###############################################################################

library(downscaling)

### Parameters #################################################################
################################################################################

resolution_of_grid <- 330 # grid cell size (in meters)
buffer <- sqrt(2)*resolution_of_grid # buffer to add to edges of grid

#LANDCOVER INFORMATION:
# -var palette =
#[
#  '#c4281b', // 1.urban
#  '#e49634', // 2.croplands
#  '#397e48', // 3. forest
#  '#DFC35A', // 4. shrub & scrub
#  '#88af52', // 5. grass
#  '#429ae4', // 6. water
#  '#7c85c9', // 7. wetlands
#  '#a6a6a6' // 8. bare surface
#  //'#B39FE1' // 9. ice and snow
#];


landcover_key <- c("1" = "urban",
                    "2" = "croplands",
                    "3" = "forest",
                    "4" = "shrub & scrub",
                    "5" = "grass",
                    "6" = "water",
                    "7" = "wetlands",
                    "8" = "bare surface")


landcover_rast <- terra::rast("/projectnb/buultra/SIF_downscaling/russell/data/WorldCover_I95_large_bbox.tif")
#terra::values(landcover_rast) <- landcover_key[terra::values(landcover_rast)]

median_albedo_rast <- terra::rast("/projectnb/buultra/SIF_downscaling/russell/data/median_boston_albedo_june_2022_largebbox.tif")

################################################################################
################################################################################

### 1) create grid -------------------------------------------------------------

target_grid <- generate_grid(downscaling::soundings,resolution_of_grid,buffer = buffer)
usethis::use_data(target_grid, overwrite = TRUE)

### 2) augment grid using median albedo values ---------------------------------

target_grid <- summarize_SpatRaster_mean_over_grid(median_albedo_rast,target_grid)

### 3) augment grid using landcover values -------------------------------------

target_grid <- summarize_SpatRaster_class_representation_over_grid(landcover_rast,target_grid)

# DELETE THIS AFTER WE GET MORE ALBEDO DATA

target_grid <- target_grid[1:(dim(target_grid)[1]-1),]

################################################################################
########## SAVE DATA ###########################################################
usethis::use_data(target_grid, overwrite = TRUE)

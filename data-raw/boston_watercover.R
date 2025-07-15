library(terra)
library(sf)
library(downscaling)

# Load the raster
landcover_rast <- terra::rast("/projectnb/buultra/SIF_downscaling/russell/data/WorldCover_I95_large_bbox.tif")

# Identify water: class 6
is_water <- landcover_rast == 6

# Convert only water regions to polygons
water_poly <- terra::as.polygons(is_water, dissolve = TRUE)  # dissolve adjacent water

# Drop any non-water values (in case dissolve didn't remove them)
water_poly <- water_poly[water_poly[[1]] == 1, ]  # keep only cells where is_water == TRUE

# Convert to sf
water_sf <- sf::st_as_sf(water_poly)

# Merge to a single multipolygon
water_union <- sf::st_union(water_sf)

# (Optional) Drop Z-dimension if present
water_union <- sf::st_zm(water_union, drop = TRUE)

# Result: a single MULTIPOLYGON sf object
print(water_union)

boston_watercover <- sf::st_as_sf(water_union)

boston_watercover <- sf_with_complement(boston_watercover)

usethis::use_data(boston_watercover,overwrite = TRUE)

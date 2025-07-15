#' Plot a fitted spatial field with toggles
#'
#' Interactive leaflet visualization of the posterior mean and observations.
#'
#' @param x A fitted_spatial_field object (from fit_spatial_field()).
#' @param basemap Character. One of "OpenStreetMap", "Esri.WorldImagery", or "CartoDB.Positron".
#' @param opacity Numeric between 0 and 1. Opacity of the grid polygons.
#' @param ... Unused.
#'
#' @export
plot.fitted_spatial_field <- function(
    x,
    basemap = "OpenStreetMap",
    opacity = 0.7,
    ...
) {
  stopifnot(inherits(x, "fitted_spatial_field"))
  if (!requireNamespace("leaflet", quietly = TRUE)) {
    stop("The `leaflet` package is required. Please install it.")
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("The `sf` package is required. Please install it.")
  }

  valid_basemaps <- c("OpenStreetMap", "Esri.WorldImagery", "CartoDB.Positron")
  if (!basemap %in% valid_basemaps) {
    stop("`basemap` must be one of: ", paste(valid_basemaps, collapse = ", "))
  }

  # Reconstruct grid polygons (the basis geometry)
  grid_geom <- NULL
  if (inherits(x$spatial_data$basis, "geometry_basis")) {
    grid_geom <- x$spatial_data$basis$geometry
  } else {
    stop("Plotting is only supported for geometry_basis at the moment.")
  }


  # Posterior mean layer
  grid_sf <- sf::st_sf(
    posterior_mean = x$posterior_mean,
    geometry = sf::st_geometry(grid_geom),
    crs = sf::st_crs(grid_geom)
  )

  m <- leaflet::leaflet() %>%
    leaflet::addProviderTiles(provider = basemap) %>%
    leaflet::addPolygons(
      data = grid_sf,
      group = "Predicted",
      fillColor = ~leaflet::colorNumeric("viridis", domain = grid_sf$posterior_mean)(posterior_mean),
      fillOpacity = opacity,
      color = NA,   # <--- No border color
      weight = 0,   # <--- No border thickness
      label = ~sprintf("Posterior mean: %.3f", posterior_mean)
    ) %>%
    leaflet::addLegend(
      pal = leaflet::colorNumeric("viridis", domain = grid_sf$posterior_mean),
      values = grid_sf$posterior_mean,
      title = "Posterior Mean",
      group = "Predicted"
    )

  # Optional observation layer
  if (!is.null(x$spatial_data$observation_data)) {
    obs <- x$spatial_data$observation_data
    response_vals <- as.vector(x$spatial_data$y)

    obs_sf <- sf::st_sf(
      response = response_vals,
      geometry = sf::st_geometry(obs),
      crs = sf::st_crs(obs)
    )

    if (any(sf::st_geometry_type(obs_sf) %in% c("POINT", "MULTIPOINT"))) {
      m <- m %>%
        leaflet::addCircleMarkers(
          data = obs_sf,
          group = "Observed",
          radius = 4,
          fillOpacity = 0.9,
          stroke = FALSE,
          color = ~leaflet::colorNumeric("plasma", domain = obs_sf$response)(response),
          label = ~sprintf("Observed: %.3f", response)
        ) %>%
        leaflet::addLegend(
          pal = leaflet::colorNumeric("plasma", domain = obs_sf$response),
          values = obs_sf$response,
          title = "Observed",
          group = "Observed"
        )
    } else {
      m <- m %>%
        leaflet::addPolygons(
          data = obs_sf,
          group = "Observed",
          fillOpacity = 0.9,
          color = ~leaflet::colorNumeric("plasma", domain = obs_sf$response)(response),
          label = ~sprintf("Observed: %.3f", response)
        ) %>%
        leaflet::addLegend(
          pal = leaflet::colorNumeric("plasma", domain = obs_sf$response),
          values = obs_sf$response,
          title = "Observed",
          group = "Observed"
        )
    }
  }

  # Add layer control toggle
  m <- m %>%
    leaflet::addLayersControl(
      overlayGroups = c("Predicted", "Observed"),
      options = leaflet::layersControlOptions(collapsed = FALSE)
    )

  m
}

#' Plot a smooth fitted spatial field from function_basis
#'
#' Visualizes the posterior mean over a dense grid using Leaflet.
#'
#' @param x A fitted_spatial_field object with function_basis.
#' @param grid_resolution Integer. Number of grid cells per axis (default: 200).
#' @param basemap Character. One of "OpenStreetMap", "Esri.WorldImagery", or "CartoDB.Positron".
#' @param opacity Numeric between 0 and 1. Raster opacity.
#' @param ... Unused.
#'
#' @return A leaflet map.
#'
#' @export
plot_function_basis_field <- function(
    x,
    grid_resolution = 200,
    basemap = "OpenStreetMap",
    opacity = 0.8,
    ...
) {
  # Validation
  stopifnot(inherits(x, "fitted_spatial_field"))

  basis <- x$spatial_data$basis

  if (!inherits(basis, "function_basis")) {
    stop("This plotting function requires a fitted_spatial_field with a function_basis.")
  }

  valid_basemaps <- c("OpenStreetMap", "Esri.WorldImagery", "CartoDB.Positron")
  if (!basemap %in% valid_basemaps) {
    stop("`basemap` must be one of: ", paste(valid_basemaps, collapse = ", "))
  }

  # Determine domain (basis domain)
  domain <- basis$metadata$domain_square
  if (is.null(domain)) {
    stop("The basis metadata must include `domain_square`.")
  }
  xmin_basis <- domain[1]
  xmax_basis <- domain[2]
  ymin_basis <- domain[3]
  ymax_basis <- domain[4]

  # Default to basis domain
  obs_bbox <- sf::st_bbox(x$spatial_data$observation_data)


  # TEMPORARY - WORKS WELL WITH SINS/COS
  xmin <- obs_bbox$xmin
  xmax <- obs_bbox$xmax
  ymin <- obs_bbox$ymin
  ymax <- obs_bbox$ymax

  # If observation data exists, restrict to intersection of bbox
  # if (!is.null(x$spatial_data$observation_data)) {
  #   obs_bbox <- sf::st_bbox(x$spatial_data$observation_data)
  #   # Intersect bounding boxes
  #   xmin <- max(xmin_basis, obs_bbox$xmin)
  #   xmax <- min(xmax_basis, obs_bbox$xmax)
  #   ymin <- max(ymin_basis, obs_bbox$ymin)
  #   ymax <- min(ymax_basis, obs_bbox$ymax)
  #   # Validate
  #   if (xmin >= xmax || ymin >= ymax) {
  #     stop("Observation data bounding box does not overlap the basis domain.")
  #   }
  # }

  # Build grid
  xs <- seq(xmin, xmax, length.out = grid_resolution)
  ys <- seq(ymin, ymax, length.out = grid_resolution)
  grid <- expand.grid(x = xs, y = ys)

  # Evaluate field
  phi_mat <- basis$evaluate_fn(as.matrix(grid)) # [n_grid, k]
  preds <- as.vector(phi_mat %*% x$posterior_mean)

  # Create raster matrix (matrix expected by leaflet raster plotting)
  raster_matrix <- matrix(preds, nrow = grid_resolution, ncol = grid_resolution, byrow = TRUE)
  raster_matrix <- raster_matrix[grid_resolution:1, ]

  # Create raster object
  if (!requireNamespace("raster", quietly = TRUE)) {
    stop("The `raster` package is required. Please install it.")
  }

  # Robust CRS logic
  crs_obj <- NULL

  # Try basis CRS first
  if (!is.null(basis$crs)) {
    tmp <- sf::st_crs(basis$crs)
    if (!is.na(tmp)) {
      crs_obj <- tmp
    }
  }

  # Fallback to observation_data CRS
  if (is.null(crs_obj) && !is.null(x$spatial_data$observation_data)) {
    obs_geom <- sf::st_geometry(x$spatial_data$observation_data)
    tmp <- sf::st_crs(obs_geom)
    if (!is.na(tmp)) {
      crs_obj <- tmp
    }
  }

  # Extract PROJ4 string or NA
  if (!is.null(crs_obj) && !is.na(crs_obj)) {
    crs_string <- crs_obj$proj4string
    if (is.null(crs_string) || is.na(crs_string)) {
      crs_string <- NA
    }
  } else {
    crs_string <- NA
  }

  # Create raster
  r <- raster::raster(
    raster_matrix,
    xmn = xmin,
    xmx = xmax,
    ymn = ymin,
    ymx = ymax,
    crs = crs_string
  )

  # Color palette
  pal <- leaflet::colorNumeric("viridis", domain = raster::values(r), na.color = "transparent")

  # Build leaflet map
  m <- leaflet::leaflet() %>%
    leaflet::addProviderTiles(provider = basemap) %>%
    leaflet::addRasterImage(
      r,
      colors = pal,
      opacity = opacity,
      project = FALSE
    ) %>%
    leaflet::addLegend(
      pal = pal,
      values = raster::values(r),
      title = "Posterior Mean"
    )

  m
}

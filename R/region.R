#' Create an Atomic Spatial Region (Triangle or Quadrilateral)
#'
#' Constructs a region object representing a simple polygon (typically a triangle or quadrilateral)
#' used for spatial integration. Automatically infers the shape based on the number of coordinates.
#'
#' @param coords A numeric matrix with 3 or 4 rows and 2 columns representing the polygon vertices.
#'               Rows should be ordered counterclockwise (if possible) and represent (X, Y) pairs.
#' @param geometry Optional. An `sf` geometry object associated with the region (e.g., for area-based integration).
#'
#' @return An object of class `"region_atomic"` and `"region"` with fields:
#' \describe{
#'   \item{type}{The string `"atomic"`.}
#'   \item{coords}{The 3 × 2 or 4 × 2 matrix of coordinates.}
#'   \item{geometry}{Optional associated `sf` geometry.}
#'   \item{shape}{Either `"triangle"` or `"quad"` depending on `nrow(coords)`.}
#' }
#'
#' @examples
#' coords <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
#' r <- region_atomic(coords)
#' str(r)
#'
#' @export
region_atomic <- function(coords, geometry = NULL) {
  stopifnot(
    is.matrix(coords),
    ncol(coords) == 2,
    nrow(coords) >= 3
  )

  n <- nrow(coords)
  shape <- if (n == 3) {
    "triangle"
  } else if (n == 4) {
    "quad"
  } else if (n > 4) {
    "polygon"
  } else {
    stop("Invalid number of coordinates: must be at least 3.")
  }

  structure(
    list(
      type = "atomic",
      coords = coords,
      geometry = geometry,
      shape = shape
    ),
    class = c("region_atomic", "region")
  )
}

#' Create a composite region (e.g., triangulated polygon)
#'
#' @param parts A list of region_atomic() objects
#' @param areas A numeric vector of the same length as parts
#' @return An S3 object of class "region_composite"
#' @export
region_composite <- function(parts, areas) {
  # Validate types
  if (!is.list(parts) || length(parts) == 0) {
    stop("`parts` must be a non-empty list of region_atomic objects.")
  }

  # Validate length consistency
  if (length(parts) != length(areas)) {
    stop("`parts` and `areas` must have the same length.")
  }

  # Validate each part is region_atomic
  if (!all(sapply(parts, inherits, "region_atomic"))) {
    stop("All elements of `parts` must be of class 'region_atomic'.")
  }

  # Validate areas
  if (!is.numeric(areas) || any(is.na(areas)) || any(areas <= 0)) {
    stop("`areas` must be a numeric vector of positive values (no NAs).")
  }

  structure(
    list(
      type = "composite",
      parts = parts,
      areas = areas
    ),
    class = c("region_composite", "region")
  )
}

#' Triangulate an sf object into individual triangle geometries
#'
#' Given an sf geometry or sf object, returns a list of triangle geometries.
#' Polygons that are already triangles are returned as-is.
#'
#' @param geom An sf geometry object (POLYGON or MULTIPOLYGON)
#' @return An sfc collection of triangle geometries
#' @export
triangulate_sf <- function(geom) {
  if (inherits(geom, "sf")) {
    geom <- sf::st_geometry(geom)
  }

  # Ensure geom is always an sfc
  if (inherits(geom, "sfg")) {
    geom <- sf::st_sfc(geom)
  }

  geom_type <- unique(sf::st_geometry_type(geom))
  # If multiple types, fallback to triangulate
  if (length(geom_type) != 1 || !geom_type %in% c("POLYGON", "MULTIPOLYGON")) {
    warning("Geometry type not POLYGON/MULTIPOLYGON; forcing triangulation.")
    tris <- sf::st_triangulate(geom)
    return(sf::st_collection_extract(tris, "POLYGON"))
  }

  # If the geometry is a collection, flatten to separate polygons
  polys <- if (geom_type == "MULTIPOLYGON") {
    sf::st_cast(geom, "POLYGON")
  } else {
    geom
  }

  out <- list()
  for (i in seq_along(polys)) {
    p <- polys[i]
    coords <- sf::st_coordinates(p)
   # print(i)

    # Extract unique (X,Y) of the exterior ring (L2==1)
    ext_coords <- coords[coords[, "L2"] == 1, c("X", "Y"), drop = FALSE]

    # Remove the closing point duplication
    if (nrow(ext_coords) >= 2 && all(ext_coords[1, ] == ext_coords[nrow(ext_coords), ])) {
      ext_coords <- ext_coords[-nrow(ext_coords), , drop = FALSE]
    }

    if (nrow(ext_coords) == 3) {
      # Already a triangle
      out[[length(out) + 1]] <- p
    } else {
      # Triangulate this polygon
      tris <- sf::st_triangulate(p)
      tris <- sf::st_collection_extract(tris, "POLYGON")
      out <- c(out, sf::st_cast(tris, "POLYGON"))
    }



  }

  # Combine results
  out_flat <- unlist(lapply(out, function(g) {
    if (inherits(g, "sfc")) {
      as.list(g)  # extract sfgs from an sfc
    } else {
      list(g)     # wrap single sfg in list
    }
  }), recursive = FALSE)


  return(sf::st_sfc(out_flat, crs = sf::st_crs(geom)))
}

#' Wrap a triangle geometry into a region_atomic object
#' @param tri A triangle sf geometry (must be exactly 3 unique exterior points)
#' @return A region_atomic object
wrap_atomic <- function(tri) {
  stopifnot(inherits(tri, "sfg") || inherits(tri, "sfc"))

  coords <- sf::st_coordinates(tri)

  # Extract only the exterior ring
  ext_coords <- coords[coords[, "L2"] == 1, c("X", "Y"), drop = FALSE]

  # Remove the closing point duplication
  if (nrow(ext_coords) >= 2 && all(ext_coords[1, ] == ext_coords[nrow(ext_coords), ])) {
    ext_coords <- ext_coords[-nrow(ext_coords), , drop = FALSE]
  }

  # Validate that we have exactly 3 unique exterior coordinates
  if (nrow(ext_coords) != 3) {
    stop(sprintf(
      "wrap_atomic: Expected triangle geometry with 3 unique exterior coordinates; got %d.",
      nrow(ext_coords)
    ))
  }

  region_atomic(ext_coords, geometry = tri)
}

#' Wrap a list of triangle geometries into a composite region
#' @param triangles A list of triangle sf geometries
#' @return A region_composite object
#'
#' @export
make_region_composite <- function(triangles) {
  if (!inherits(triangles, "sfc")) {
    # Try coercion, assuming list of sfgs or sfcs
    triangles <- sf::st_sfc(unlist(lapply(triangles, function(g) {
      if (inherits(g, "sfc")) {
        as.list(g)
      } else {
        list(g)
      }
    }), recursive = FALSE))
  }

  parts <- lapply(triangles, wrap_atomic)
  areas <- as.numeric(sf::st_area(triangles))
  region_composite(parts[areas > 0], areas[areas > 0])
}

#' Parallel intersection and triangulation with progress reporting
#'
#' Intersects each observation polygon with stratum polygons, triangulates the result,
#' and wraps into `region_composite` objects. Uses `future.apply` and `progressr` for
#' parallelism and progress monitoring.
#'
#' @param obs_polygons An `sf` object of observation polygons.
#' @param stratum_polygons An `sf` object of stratification polygons.
#'
#' @return A list of `region_composite` objects (or `NULL` for no intersections).
#'
#' @details Requires a parallel backend set with `future::plan()` and activates a progress
#' handler with `progressr::handlers()`. Wrap calls in `with_progress()` to display updates.
#'
#' @examples
#' \dontrun{
#' library(future)
#' library(progressr)
#' plan(multisession)
#' handlers(global = TRUE)
#' with_progress({
#'   result <- intersect_then_triangulate_parallel(obs_sf, strata_sf)
#' })
#' }
#'
#' @export
intersect_then_triangulate_parallel <- function(obs_polygons, stratum_polygons) {
  stopifnot(!is.null(stratum_polygons))
  message("Parallel triangulation and region construction by observation...")

  p <- progressr::progressor(steps = nrow(obs_polygons))

  future.apply::future_lapply(seq_len(nrow(obs_polygons)), function(i) {
    print(i)
    p()
    obs <- obs_polygons[i, , drop = FALSE]

    # Step 1: Fast intersection
    sub_stratum <- st_fast_intersection(obs, stratum_polygons)
    if (nrow(sub_stratum) == 0) return(NULL)

    # Step 2: Triangulate
    tris <- triangulate_sf(sub_stratum)
    if (length(tris) == 0) return(NULL)

    # Step 3: Wrap
    if (any(sf::st_is_empty(tris))) stop("Some empty triangle exists i guess")
    make_region_composite(tris)
  })
}


#' Build region list from observation and (optional) stratum polygons
#'
#' @param obs_polygons An sf object of observation polygons
#' @param stratum_polygons Optional sf object of stratifying polygons
#' @param method One of "raw", "triangulate", or "intersect-then-triangulate"
#' @return A list of region_composite or region_atomic objects (some may be NULL)
#'
#' @export
build_region_list <- function(obs_polygons, stratum_polygons = NULL,
                              method = c("triangulate", "intersect-then-triangulate", "raw")) {
  method <- match.arg(method)
  n <- nrow(obs_polygons)

  if(!is.null(stratum_polygons)){stopifnot(sf::st_crs(obs_polygons) == sf::st_crs(stratum_polygons))}

  if (method == "raw") {
    return(lapply(seq_len(n), function(i) {
      geom <- sf::st_geometry(obs_polygons)[[i]]
      coords <- sf::st_coordinates(geom)[, c("X", "Y")]
      coords <- coords[!duplicated(coords), , drop = FALSE]
      if (nrow(coords) < 3) return(NULL)
      region_atomic(coords, geometry = geom)
    }))
  }

  if (method == "triangulate") {
    return(lapply(seq_len(n), function(i) {
      geom <- sf::st_geometry(obs_polygons)[[i]]
      tris <- triangulate_sf(geom)
      if (length(tris) == 0) return(NULL)
      make_region_composite(tris)
    }))
  }

  if (method == "intersect-then-triangulate") {
    stopifnot(!is.null(stratum_polygons))

    message("Using parallel intersection and triangulation...")
    return(intersect_then_triangulate_parallel(obs_polygons, stratum_polygons))
  }
}

#' Convert a composite region into an sf object for visualization
#'
#' @param composite A region_composite object
#'
#' @return An sf object containing all atomic parts of the region
#' @export


composite_region_to_sf <- function(composite) {
  stopifnot(inherits(composite, "region_composite"))

  geoms <- lapply(composite$parts, function(part) {
    stopifnot(inherits(part, "region_atomic"))
    part$geometry
  })

  # Remove NULLs (just in case)
  geoms <- Filter(Negate(is.null), geoms)

  sf::st_sf(geometry = sf::st_sfc(geoms))
}


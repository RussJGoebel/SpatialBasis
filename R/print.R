#' @export
print.region <- function(x, ...) {
  if (x$type == "atomic") {
    cat("Region (atomic) with", nrow(x$coords), "vertices\n")
  } else {
    cat("Region (composite) with", length(x$parts), "parts\n")
  }
  invisible(x)
}

#' Print Method for Spatial Basis Objects
#'
#' Displays a summary of a basis object, including its class, type, number of basis functions,
#' CRS, and any available geometry.
#'
#' @param x A basis object (typically of class `"geometry_basis"` or `"function_basis"`).
#' @param ... Ignored.
#' @export
print.basis <- function(x, ...) {
  cat("Spatial Basis Object\n")
  cat("  Class:    ", paste(class(x), collapse = " > "), "\n")

  # Print type if available
  if (!is.null(x$type)) {
    cat("  Type:     ", x$type, "\n")
  } else {
    cat("  Type:     [unspecified]\n")
  }

  cat("  k:        ", x$k, "basis functions\n")

  # Print CRS
  if (!is.null(x$crs)) {
    cat("  CRS:      ", sf::st_as_text(x$crs, pretty = TRUE), "\n")
  } else {
    cat("  CRS:      [none specified]\n")
  }

  # Print geometry information if applicable
  if (inherits(x, "geometry_basis") && !is.null(x$geometry)) {
    cat("  Geometry: ", nrow(x$geometry), "polygon(s)\n")
  } else if (inherits(x, "function_basis")) {
    cat("  Geometry: [not used]\n")
  }

  invisible(x)
}

#' Print Method for Prior Objects
#'
#' Nicely formats a prior object, showing its type, parameters, and dimension.
#'
#' @param x A prior object created by e.g. make_sar_prior().
#' @param ... Ignored.
#'
#' @export
print.prior <- function(x, ...) {
  cat("Spatial Prior\n")
  cat("  Type:        ", x$type, "\n")
  if (!is.null(x$params)) {
    cat("  Parameters:\n")
    for (nm in names(x$params)) {
      val <- x$params[[nm]]
      if (is.numeric(val) && length(val) == 1) {
        cat("    -", nm, "=", signif(val, 4), "\n")
      } else if (is.character(val)) {
        cat("    -", nm, "=", val, "\n")
      } else if (is.numeric(val)) {
        cat("    -", nm, "= Numeric vector of length", length(val), "\n")
      } else {
        cat("    -", nm, "= (complex object)\n")
      }
    }
  }
  cat("  Dimension:   ", nrow(x$precision), "×", ncol(x$precision), "\n")
  cat("  Sparsity:    ", round(100 * (1 - mean(as.vector(x$precision) != 0)), 2), "% zeros\n")
  invisible(x)
}

#' Print method for spatial_field objects
#'
#' @param x A spatial_field object
#' @param ... Ignored
#' @export
print.spatial_field <- function(x, ...) {
  cat("Spatial Field:\n")
  cat("  Label:          ", x$label, "\n")
  if (!is.null(x$basis)) {
    cat("  Basis:          ", class(x$basis)[1], sprintf("(%d functions)", x$basis$k), "\n")
  } else {
    cat("  Basis:           [None]\n")
  }
  if (!is.null(x$prior)) {
    cat("  Prior:          ", x$prior$type, "\n")
  } else {
    cat("  Prior:           [None]\n")
  }
  if (!is.null(x$integration_rule)) {
    cat("  Integration:    ", x$integration_rule$type, "\n")
  } else {
    cat("  Integration:     [None]\n")
  }
  invisible(x)
}

#' @export
print.prepared_spatial_data <- function(x, ...) {
  cat("Prepared Spatial Data\n")
  cat("  Basis:\n")
  cat("    - Class:      ", paste(class(x$basis), collapse = " > "), "\n")
  cat("    - k:          ", x$basis$k, "basis functions\n")
  if (!is.null(x$basis$crs)) {
    cat("    - CRS:        ", sf::st_as_text(x$basis$crs, pretty = TRUE), "\n")
  } else {
    cat("    - CRS:        [none specified]\n")
  }

  cat("  Prior:\n")
  if (inherits(x$prior, "Matrix") || is.matrix(x$prior)) {
    cat("    - Dimension:  ", paste(dim(x$prior), collapse = " × "), "\n")
    cat("    - Type:       Precision matrix\n")
  } else {
    cat("    - Type:       [unknown]\n")
  }

  cat("  Integration:\n")
  cat("    - Rule:       ", paste(class(x$integration_rule), collapse = " > "), "\n")
  cat("    - Observations:", nrow(x$X_spatial), "\n")
  cat("    - k columns:  ", ncol(x$X_spatial), "\n")
  cat("    - Region list:", length(x$region_list), "regions\n")

  invisible(x)
}

#' @export
print.observation_data <- function(x, ...) {
  cat("Observation Data Object\n")
  cat("  Observations:", length(x$y), "\n")
  cat("  CRS:", sf::st_as_text(x$sf_metadata$crs), "\n")
  cat("  Region list:", length(x$region_list), "regions\n")
  invisible(x)
}

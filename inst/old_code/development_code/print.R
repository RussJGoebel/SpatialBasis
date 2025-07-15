#' Print method for basis objects
#'
#' @param x A basis object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.basis <- function(x, ...) {
  cat("<basis object>\n")
  cat("  Number of basis functions:", x$k, "\n")
  if (!is.null(x$metadata) && length(x$metadata) > 0) {
    cat("  Metadata keys:\n")
    cat("   -", paste(names(x$metadata), collapse = ", "), "\n")
  }
  invisible(x)
}

#' Print method for prior_on_basis objects
#'
#' Provides a concise summary of the prior defined on basis coefficients,
#' including the number of basis functions, prior type (diagonal or full precision),
#' and any associated hyperparameters or metadata.
#'
#' @param x A `prior_on_basis` object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input `x`.
#' @export
print.prior_on_basis <- function(x, ...) {
  if (!inherits(x, "prior_on_basis")) {
    stop("Object is not a prior_on_basis.")
  }

  k <- x$basis$k
  cat("<prior_on_basis>\n")
  cat("  Basis dimension:", k, "\n")

  if (!is.null(x$variances)) {
    cat("  Prior: diagonal precision (from variances)\n")
  } else {
    cat("  Prior: full precision matrix\n")
  }

  if (!is.null(x$hyperparams) && length(x$hyperparams) > 0) {
    cat("  Hyperparameters:\n")
    for (nm in names(x$hyperparams)) {
      val <- x$hyperparams[[nm]]
      summary <- if (length(val) == 1 && is.atomic(val)) format(val) else paste0(class(val), " [", length(val), "]")
      cat("    -", nm, "=", summary, "\n")
    }
  }

  if (!is.null(x$metadata) && length(x$metadata) > 0) {
    cat("  Metadata keys:\n")
    cat("   -", paste(names(x$metadata), collapse = ", "), "\n")
  }

  invisible(x)
}

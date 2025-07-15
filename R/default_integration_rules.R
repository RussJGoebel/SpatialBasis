#' Determine Default Integration Rule for a Basis
#'
#' Returns the appropriate integration rule for a given basis type.
#'
#' @param basis A basis object (e.g., from make_geometry_basis(), make_matern_fourier_basis(), etc.)
#'
#' @return An integration_rule object
#' @export
default_integration_rule <- function(basis) {
  UseMethod("default_integration_rule")
}

#' @export
default_integration_rule.geometry_basis <- function(basis) {
  integration_rule_area_overlap()
}

#' @export
default_integration_rule.piecewise_constant_basis <- function(basis) {
  integration_rule_area_overlap()
}

#' @export
default_integration_rule.matern_fourier_basis <- function(basis) {
  integration_rule_qmc_average(generate_qmc_unit_square(256))
}

#' @export
default_integration_rule.function_basis <- function(basis) {
  integration_rule_qmc_average(generate_qmc_unit_square(256))
}

#' @export
default_integration_rule.default <- function(basis) {
  stop("No default integration rule defined for basis of class: ", class(basis)[1])
}

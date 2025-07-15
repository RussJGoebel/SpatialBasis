#'
#' @param object A basis object.
#' @param newdata A matrix, data.frame, or `sf` object of coordinates.
#' @param ... Ignored.
#'
#' @return A matrix of basis evaluations.
#' @export
predict.basis <- function(object, newdata, ...) {
  object$evaluate(newdata)
}

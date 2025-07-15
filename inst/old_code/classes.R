################################################################################
############################# CLASSES ##########################################
################################################################################
# There are two classes relevant to the downscaling in this package:
#
# 1) sounding // the satellite measurements consisting of potentially irregular
#    sf polygons and associated measurements
# 2) grid // the target grid, which is assumed to consist of square polygons
#    of equal area
#
################################################################################

#' new_soundings
#'
#' A constructor function for an S3 class to be used/interpreted as satellite soundings.
#' This is essentially an sf_object that knows which column is the measurement values.
#' This is useful for readable syntax both for downscaling and visualizing the object.
#' We check that the sf_object is a multipolygon.
#'
#' @param sf_object An 'sf' object of type 'MULTIPOLYGON'.
#' @param column Either an integer denoting the place of a column or the name of a column
#' @param column An integer or string denoting which column of sf_object is
#' of interest. This is extracted with `sf_object[[column]]` or `sf_object[,column]`.
#'
#' @return An object of class 'soundings'
#' @export
#'
#' @examples
#'
#' soundings
#'
new_soundings <- function(sf_object,
                          column){

  if (!inherits(sf_object, "sf") ||
      !(sf::st_geometry_type(sf_object, by_geometry = FALSE) %in% c("MULTIPOLYGON","POLYGON"))) {
    stop("The sf_object input must be an 'sf' object of type 'MULTIPOLYGON' or 'POLYGON'.")
  }

  if(is.null(sf_object[[column]])){
    stop("Column specified by 'column' in sf_object does not exist.
          Check that 'column' is either the name of a column in sf_object or
          an integer denoting which column to use.")
  }
  if(!is.numeric(sf_object[[column]])){
    stop("'column' corresponds to a column that is not numeric. Make sure that
           the corresponding column in sf_object is numeric.")
  }

  class(sf_object) <- c("soundings",class(sf_object))
  attr(sf_object, "column") <- column

  return(sf_object)

}

#' new_downscaling results
#'
#' A constructor function for an S3 class to be used/interpreted as the results from
#' the downscaling procedure.
#' .
#' This is essentially an sf_object that contains:
#' -geometries
#' -associated posterior means
#' -anything in the original sf object
#' -possibly associated posterior variance, and/or variance parameters
#' -
#'
#'
#' @param sf_object An 'sf' object of type 'MULTIPOLYGON'
#'
#' @return An object of class downscaling_results.
#' @export
#'
#' @examples
new_downscaling_results <- function(sf_object){

  if (!inherits(sf_object, "sf") || sf::st_geometry_type(sf_object, by_geometry = FALSE) != "MULTIPOLYGON") {
    stop("The sf_object input must be an 'sf' object of type 'MULTIPOLYGON'.")
  }

  if (!("posterior_mean" %in% colnames(sf_object))){
    stop("The sf_object input must have a column named 'posterior_mean'.")
  }

  class(sf_object) <- c("downscaling_results",class(sf_object))


  return(downscaling_results)


}



#' Downscaling Model Constructor
#'
#' @param soundings
#' @param target_grid
#' @param column
#'
#' @return
#' @export
#'
#' @examples
new_downscale_model <- function(soundings,
                                target_grid,
                                column = "SIF_757nm",
                                A = NULL,
                                AtA = NULL,
                                Atb = NULL,
                                W = NULL,
                                lambda = NULL,
                                rho = NULL,
                                linear_solver = "cg") {

  structure(
    list(
      data = list(soundings = soundings, target_grid = target_grid, column = column),
      matrices = list(A = A, AtA = AtA, Atb = Atb, W = W),
      hyperparams = list(lambda = lambda, rho = rho, method = "fixed"),
      solver = list(linear_solver = linear_solver),
      results = list(),
      diagnostics = list()
    ),
    class = "downscale_model"
  )
}


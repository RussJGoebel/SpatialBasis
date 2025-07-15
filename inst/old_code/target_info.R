

#' target_datanames
#'
#' @param product name of data product (L2StdTG, LtSIF, ...)
#' @param targetname target string name (belchatowPl, bostonMA, ...)
#' @param targetdate 6-character data string (YYMMDD)
#' @param lambda numeric or list of numerics used for sensitivity
#' @param resolution_m downscaled cell resolution in meters (330, 200, 500, ...)
#' @param Amat_type specify polygon (PG) or 2D Gaussian (2D) intersection matrix
#'
#' @return dataframe containing string names and parameters for specified downscaled data and intermediate products
#' @export
#'
#' @examples tg_belchatowPl_xco2 <- target_datanames(product = "L2StdTG", targetname = "belchatowPl", column = "XCO2_ppm", targetdate = c("220324"), resolution_m = c(250), lambda = c(5))
target_datanames <- function(product, targetname, targetdate, column, lambda, resolution_m)
{
  product <- match.arg(product, c("L2StdTG", "LtSIF", "LtCO2"))

  df <- data.frame()
  for (date in targetdate) {
    for (res in resolution_m) {
      for (lam in lambda) {
        str_resolution <- match.arg(paste(res), c("250", "330", "500", "750", "1000"))
        str_resolution <- paste0(str_resolution,"m")

        str_lambda <- downscaling::number2string(lam, multiplier = 10, pad=T, digits=3)

        strN <- paste(product, targetname, date, sep='_')

        df <-  dplyr::bind_rows(df, data.frame(
          strN   = strN,
          strN1  = paste("border", strN, sep='_'),
          strG0  = paste("grid", str_resolution, product, targetname, sep='_'),
          strG   = paste("grid", str_resolution, strN, sep='_'),
          strAPG = paste("A", "PG", str_resolution, strN, sep='_'),
          strA2D = paste("A", "2D", str_resolution, strN, sep='_'),
          strDPG = paste("D", "PG", str_resolution, str_lambda, strN, sep='_'),
          strD2D = paste("D", "2D", str_resolution, str_lambda, strN, sep='_'),
          lambda = lam,
          resolution_m = res, # m
          date   = date, # YYMMDD
          targetname   = targetname,
          column = column,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  return (df)
}

#' Title
#'
#' @param target_name string name for target
#'
#' @return NULL. target bounding box coordinates in degrees automatically added to current environment (lat_min, lat_max, lon_min, lon_max)
#' @export
#'
#' @examples bounds <- get_target_bounds("bostonMA")
get_target_bounds <- function(targetname) {
  str <- tolower(targetname)
  if(str %in% c("boston", "bostonma")) {
    # This is the *larger* bounding box.
    lat_min <<-  42.2207
    lat_max <<-  42.53601
    lon_min <<- -71.26929
    lon_max <<- -70.95319
  } else if(str %in% c("lamont", "lamontok")){
    lat_min <<-  36.24495
    lat_max <<-  36.8359
    lon_min <<- -97.8228
    lon_max <<- -97.15583
  } else if(str %in% c("toronto", 'torontoca')) {
    lat_min <<-  43.52777
    lat_max <<-  43.98577
    lon_min <<- -79.69979
    lon_max <<- -79.19504
  } else if (str %in% c("belchatow", "belchatowpl")) {
    lat_min <<-  51.0000
    lat_max <<-  51.6900
    lon_min <<-  19.0900
    lon_max <<-  19.5500
  } else {
    message(paste("unregistered name:", targetname))
    lat_min <<- NULL
    lat_max <<- NULL
    lon_min <<- NULL
    lon_max <<- NULL
  }
  return (c(lat_min, lat_max, lon_min, lon_max))
}

#' upscale
#'
#' @param sf_object_small
#' @param sf_object_large
#' @param column
#' @param gaussian_noise_variance
#'
#' @return Given a small sf_object with a values field, upscale with a small amount of gaussian nose.
#' @export
#'
#' @examples
upscale <- function(sf_object_small,
                    sf_object_large,
                    column,
                    gaussian_noise_variance){

n <- dim(sf_object_large)[1]

message("Computing A matrix...")
A <- compute_A_matrix(sf_object_small,sf_object_large)

message("Upscaling with gaussian noise...")
sf_object_large$upscaled_values <- as.vector(A %*% sf_object_small[[column]]+rnorm(n,0,sd = gaussian_noise_variance^2))

message("Upscaled with variance ",gaussian_noise_variance)

return(sf_object_large)

}

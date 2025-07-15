#' get_basepath
#'
#' @return path to local or BU cluster datasets. Eather /home/data/oco_data/ or /projectnb/buultra/oco_data
#' @export
#'
#' @examples dir <- get_basepath()
get_basepath <- function() {
  for (basepath in c("/home/data/oco_data",
                     "/projectnb/buultra/oco_data" )) {
    if (dir.exists(basepath)) { break }
  }
  return(normalizePath(basepath, winslash = "/", mustWork = T))
}

#' get_directory_files
#'
#' @param dirpath path to directory
#' @param pattern file search pattern
#' @param recursive search subdirectories in dirpath
#'
#' @return list of files matching with matching pattern (wildcard). recursive=TRUE (default) include all subdirectories in search. 
#' @export
#'
#' @examples fnames <- get_directory_files("oco_data/", pattern='*.h5') 
get_directory_files <- function(dirpath, pattern='*.nc4', recursive = T) {
  dirpath = normalizePath(dirpath, winslash='/', mustWork=T)
  fnames <- list.files(path=dirpath, pattern = pattern,
                       recursive = recursive, full.names = T)
  sapply(fnames, unique)
  return (fnames)
}

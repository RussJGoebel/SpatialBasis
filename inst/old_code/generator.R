
#' generator
#'
#' @param tg populated target_dataname
#' @param save_dir directory for temporary data storage (must exist)
#'
#' @return NULL. required values in target_dataframe are automatically saved to global environment in function
#' @export
#'
#' @examples generator(tg_xco2_belchatownPl, save_dir="../tmp/")
generator <- function(tg, save_dir = "../tmp/") {
  save_dir <- normalizePath(save_dir, winslash = '/', mustWork = TRUE)

  if (!verify_string_name(tg$strN)) {
    warning(paste("Not Found:", tg$strN)); return;
  }
  #overwrites original strN if found in save_dir
  if (verify_file_name(tg$strN, save_dir = save_dir)) {
    message(paste("-> Reading", tg$strN))
    assign(tg$strN, readRDS(file=paste0(save_dir, "/", tg$strN, ".rda")), envir = .GlobalEnv)
  }

  if (!verify_string_name(tg$strG0)) {
    warning(paste("Not Found:", tg$strG0)); return;
  }

  if (!verify_file_name(tg$strN1, save_dir = save_dir)) {
    message(paste("# Generating", tg$strN1))
    assign(tg$strN1, sf::st_convex_hull(sf::st_union(get(tg$strN))), .GlobalEnv)
    saveRDS(get(tg$strN1), file=paste0(save_dir,"/",tg$strN1,".rda"))
  }
  if (!verify_string_name(tg$strN1)) {
    message(paste("-> Reading", tg$strN1))
    assign(tg$strN1, readRDS(file=paste0(save_dir, "/", tg$strN1, ".rda")), envir = .GlobalEnv)
  }

  if (!verify_file_name(tg$strG, save_dir = save_dir)) {
    message(paste("# Generating", tg$strG))
    timer0$start(tg$strG)

    grid_tmp <- downscaling::buffered_convex_hull(soundings = get(tg$strN),
                                                  grid=get(tg$strG0),
                                                  buffer = 2*tg$resolution_m)

    overs <- sf::st_intersects(grid_tmp, get(tg$strN1), sparse=F)
    overs[(sapply(overs, length))] <- F
    overs[!(sapply(overs, length))] <- T
    grid_tmp$buffer <- unlist(overs)

    assign(tg$strG, grid_tmp, envir = .GlobalEnv)

    timer0$stop(tg$strG)
    saveRDS(get(tg$strG), file=paste0(save_dir,"/",tg$strG,".rda"))
  }
  if (!verify_string_name(tg$strG)) {
    message(paste("-> Reading", tg$strG))
    assign(tg$strG, readRDS(file=paste0(save_dir, "/", tg$strG, ".rda")), envir = .GlobalEnv)
  }

  # polynomial
  if (!verify_file_name(tg$strAPG, save_dir = save_dir)) {
    message(paste("# Generating", tg$strAPG))
    timer0$start(tg$strAPG)
    assign(tg$strAPG,
           compute_A_matrix(get(tg$strG), get(tg$strN)), envir = .GlobalEnv)
    timer0$stop(tg$strAPG)
    saveRDS(get(tg$strAPG), file=paste0(save_dir,"/",tg$strAPG,".rda"))
  }
  if (!verify_string_name(tg$strAPG)) {
    message(paste("-> Reading", tg$strAPG))
    assign(tg$strAPG, readRDS(file=paste0(save_dir, "/", tg$strAPG, ".rda")), envir = .GlobalEnv)
  }

  # 2d gaussian
  if (!verify_file_name(tg$strA2D, save_dir = save_dir)) {
    message(paste("# Generating", tg$strA2D))
    timer0$start(tg$strA2D)
    assign(tg$strA2D,
           compute_A_matrix_2d_gauss(get(tg$strG), get(tg$strN)), envir = .GlobalEnv)
    timer0$stop(tg$strA2D)
    saveRDS(get(tg$strA2D), file=paste0(save_dir,"/",tg$strA2D,".rda"))
  }
  if (!verify_string_name(tg$strA2D)) {
    message(paste("-> Reading", tg$strA2D))
    assign(tg$strA2D, readRDS(file=paste0(save_dir, "/", tg$strA2D, ".rda")), envir = .GlobalEnv)
  }

  # polynomial
  if (!verify_file_name(tg$strDPG, save_dir = save_dir)) {
    message(paste("# Generating", tg$strAPG))
    timer0$start(tg$strDPG)
    D <- downscaling::downscale(soundings = get(tg$strN),
                   target_grid = get(tg$strG), A = get(tg$strAPG),
                   column = tg$column, lambda = tg$lambda, rho = 1)
    D[[tg$column]] <- with(D, ifelse(buffer, posterior_mean, NA))
    timer0$stop(tg$strDPG)
    assign(tg$strDPG, D, envir = .GlobalEnv)
    saveRDS(get(tg$strDPG), file=paste0(save_dir,"/",tg$strDPG,".rda"))
  }
  message(paste("-> Reading", tg$strDPG))
  assign(tg$strDPG, readRDS(file=paste0(save_dir, "/", tg$strDPG, ".rda")), envir = .GlobalEnv)

  # 2d gaussian
  if (!verify_file_name(tg$strD2D, save_dir = save_dir)) {
    message(paste("# Generating", tg$strD2D))
    timer0$start(tg$strD2D)
    D <- downscaling::downscale(soundings = get(tg$strN),
                   target_grid = get(tg$strG), A = get(tg$strA2D),
                   column = tg$column, lambda = tg$lambda, rho = 1)
    D[[tg$column]] <- with(D, ifelse(buffer, posterior_mean, NA))
    timer0$stop(tg$strD2D)
    assign(tg$strD2D, D, envir = .GlobalEnv)
    saveRDS(get(tg$strD2D), file=paste0(save_dir,"/",tg$strD2D,".rda"))
  }
  message(paste("-> Reading", tg$strD2D))
  assign(tg$strD2D, readRDS(file=paste0(save_dir, "/", tg$strD2D, ".rda")), envir = .GlobalEnv)
}


#' cleaner
#'
#' @param tg populeted target_dataname
#' @param save_dir directory for temporary data storage (must exist)
#' @param delN boolean: remove native sounding file from save_dir
#' @param delG boolean: remove grid file from save_dir
#' @param delPG boolean: remove Polygon Method files from save_dir
#' @param del2D boolean: remove 2D-Gauss Method files from save_dir
#' @param dryrun boolean: TRUE (default) shows files that would be removed; FALSE permenentally removes same file
#'
#' @return NULL. values in target_dataframe are prompted and removed in function
#' @export
#'
#' @examples generator(tg_xco2_belchatownPl, save_dir="../tmp/")
cleaner <- function(tg, save_dir,
                      delN = T, delG=T,
                      delPG=F, del2D=F,
                      dryrun=T
                      ) {

  fnames <- list()
  for (idx in 1:nrow(tg) ) {

    if (delN) {
      fstr <- get_file_name(tg$strN[[idx]], save_dir = save_dir)
      if (!is.null(fstr)) { if (!(fstr %in% fnames)) { fnames <- append(fnames, fstr) }}
    }

    if (delG) {
      fstr <- get_file_name(tg$strG[[idx]], save_dir = save_dir)
      if (!is.null(fstr)) { if (!(fstr %in% fnames)) { fnames <- append(fnames, fstr) }}
    }

    if (delPG) {
      fstr <- get_file_name(tg$strAPG[[idx]], save_dir = save_dir)
      if (!is.null(fstr)) { if (!(fstr %in% fnames)) { fnames <- append(fnames, fstr) }}

      fstr <- get_file_name(tg$strDPG[[idx]], save_dir = save_dir)
      if (!is.null(fstr)) { if (!(fstr %in% fnames)) { fnames <- append(fnames, fstr) }}
    }

    if (del2D) {
      fstr <- get_file_name(tg$strA2D[[idx]], save_dir = save_dir)
      if (!is.null(fstr)) { if (!(fstr %in% fnames)) { fnames <- append(fnames, fstr) }}

      fstr <- get_file_name(tg$strD2D[[idx]], save_dir = save_dir)
      if (!is.null(fstr)) { if (!(fstr %in% fnames)) { fnames <- append(fnames, fstr) }}
    }
  }


  for (f in fnames) {
    if (dryrun) {
      message(paste("Would Remove", f))
      next
    }
    file.remove(f)
    message(paste("Removed", f))
  }

}



#' verify_string_name
#'
#' @param string_name string name of variable
#'
#' @return boolean value: TRUE if string_name is assigned in current environment
#' @export
#'
#' @examples verify_string_name("L2StdTG_belchatowPl_220324")
verify_string_name <- function(string_name) {
  if (!exists(string_name)) { return (FALSE) }
  return (TRUE)
}

#' verify_file_name
#'
#' @param string_name string name of variable
#' @param save_dir directory for temporary data storage (must exist)
#' @param ext file extension (default .rda)
#'
#' @return boolean value: true if <string_name>.<ext> is a file in <save_dir>
#' @export
#'
#' @examples verify_file_name("L2StdTG_belchatowPl_220324", save_dir="../tmp/", ext="rda")
verify_file_name <- function(string_name, save_dir="../tmp/", ext="rda") {
  if (!is.null(get_file_name(string_name, save_dir = save_dir, ext=ext))) {
    return (TRUE)
  }
  return (FALSE)
}

#' get_file_name
#'
#' @param string_name string name of variable
#' @param save_dir directory for temporary data storage (must exist)
#' @param ext file extension (default rda)
#'
#' @return string or null value: path to <string_name>.<ext> file in <save_dir>; null if does not exist.
#' @export
#'
#' @examples get_file_name("L2StdTG_belchatowPl_220324", save_dir="../tmp/", ext="rda")
get_file_name <- function(string_name, save_dir="../tmp/", ext="rda") {
  dirpath <- normalizePath(save_dir, winslash = "/", mustWork = T)
  fnames  <- list.files(path=dirpath, pattern = paste0("*.", ext), recursive = T, full.names = T)
  for (f in fnames) {
    fstr <- sub(paste0('\\.',ext,'$'), '', basename(f))
    if (grepl(string_name, fstr, fixed = TRUE)) {
      if (fstr == string_name) { return (f) }
    }
  }
  return (NULL)
}


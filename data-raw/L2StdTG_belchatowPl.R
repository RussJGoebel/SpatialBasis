library(devtools)
library(downscaling)

downscaling::create_timers()

targ_basedir <- paste(downscaling::get_basepath(),
                      "oco2_co2/targets/belchatowPl/", sep = '/')

overwrite     <- T
product       <- "L2StdTG"
pattern       <- "*.h5"
targetname    <- "belchatowPl"
resolution_m  <- c(250, 330, 500, 750)
buf_multipler <- 10

all_soundings = list()
fnames <- downscaling::get_directory_files(targ_basedir, pattern=pattern)
for (fname in fnames) {
  if (!grepl(product, fname)) next

  s <- strsplit(basename(fname), split="_")[[1]]
  prefix <- paste(s[2], sep="_")
  date   <- s[4]
  strN <- paste(prefix, targetname, date, sep="_")


  if (!downscaling::verify_string_name(strN)) {
    sounds      <- downscaling::read_oco2_v11_co2_L2StdTG(fname, targetname)
    sounds$date <- date

    assign(strN, sounds)
    do.call("use_data", list(as.name(strN), overwrite=overwrite))
  }
  all_soundings <- dplyr::bind_rows(all_soundings,  get(strN))
}

for (res in resolution_m) {
  strG0 <- paste("grid", paste0(res, "m"), product, targetname, sep="_")
  if (!exists(strG0)) {
    timer0$start(strG0)
    buffer <- buf_multipler * res
    grid   <- downscaling::generate_grid(soundings=all_soundings,
                                         res = res, buffer = buffer,
                                         grid_shape = "buffered_convex_hull")

    assign(strG0, grid)
    do.call("use_data", list(as.name(strG0), overwrite=overwrite))
    timer0$stop(strG0)
  }
}


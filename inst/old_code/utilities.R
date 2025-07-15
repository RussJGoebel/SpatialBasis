library(timeR)

#' Title
#'
#' @param timer_str string value for timeR assignment in global environment
#'
#' @return NULL. timer_str automatically assigned to timeR structure in global environment
#' @export
#'
#' @examples create_timers(timer_str = "my_timer")
create_timers <- function(timer_str = "timer0") {
  for (tstr in timer_str) {
    if (!exists(tstr)) {
      assign(tstr, timeR::createTimer(F), envir = .GlobalEnv)
      do.call("use_data", list(as.name(tstr), overwrite=T))
    }
  }
}


#' Title
#'
#' @param num number to be converter to string
#' @param multiplier premultiplier on num for string format
#' @param pad option to pad values with zeros for same length
#' @param digits number of digits in padded output. requires pad=TRUE>
#'
#' @return string value representing input number
#' @export
#'
#' @examples padded_num <- number2string(0.1, multipler=10, pad=T)
number2string <- function(num, multiplier=1, pad = F, digits=3) {

  if (pad & (num*multiplier >= 10^(digits))) {
    warning(paste(num, "*", multiplier, "has more than", digits, "digits"))
  }

  str <- paste(num*multiplier)

  if (pad & (digits > 0)) {
    for (d in (digits-1):1) {
      str <- ifelse(num < (10^(d-1)), paste0("0", str), str)
    }
  }
  return (str)
}


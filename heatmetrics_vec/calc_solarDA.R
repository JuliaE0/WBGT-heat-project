#' @title Calculate solar declination angle
#'
#' @description This function calculates the solar declination angle ("d") in degrees and time correction ("tc").
#'
#' @param jd Julian day of year (1-366, e.g., Feb. 1 = 32)
#' @param hour Hour (0-23 UTC)
#'
#' @return Returns list of outputs: solar declination angle ("d") and time offset ("tc").
#'
#' @examples calc_solarDA(40, 12)
calc_solarDA <- function(jd, hour) {

  # Calculate angular fraction of the year in radians
  #
  g <- (360 / 365.25) * (jd + (hour / 24))  # fractional year g in degrees
  g <- ifelse(g > 360, (g - 360), g)
  g_rad <-  g * (pi / 180) # convert to radians

  # Calculate the solar declination angle, lowercase delta, in degrees:
  #
  d <- 0.396372 - 22.91327 * cos(g_rad) + 4.025430 * sin(g_rad) - 0.387205 *
    cos(2 * g_rad) + 0.051967 * sin(2 * g_rad) - 0.154527 * cos(3 * g_rad) +
    0.084798 * sin(3 * g_rad)

  tc <- (0.004297 + 0.107029 * cos(g_rad) - 1.837877 * sin(g_rad) -
           0.837378 * cos(2 * g_rad) - 2.340475 * sin(2 * g_rad))

  outputs <- list("d" = d, "tc" = tc)
  return(outputs)
}

#' @title Calculate the cosine solar zenith angle integrated over the hour
#'
#' @description To calculate the integrated cosine of solar zenith angle (cza) for an hour.
#' The calc_cza() function can be used by itself instead, but it is less accurate for sunrise and sunset
#' hours. See Hogan and Hirahara (2016) [https://doi.org/10.1002/2015GL066868] and Di Napoli
#' (2020) [https://doi.org/10.1007/s00484-020-01900-5] for more details. This function is
#' only for use with hourly time steps.
#'
#' @param lat Degrees north latitude (-90 to 90)
#' @param lon Degrees east longitude (-180 to 180)
#' @param y Year (four digits, e.g., 2020)
#' @param mon Month (1-12)
#' @param d Day of month (whole number)
#' @param hr Hour (0-24 UTC)
#'
#' @return Returns the integrated cosine of the solar zenith angle (cza).
#'
#' @examples v_calc_cza_int(30, -100, 2020, 1, 1, 12)
v_calc_cza_int <- function(lat, lon, y, mon, d, hr) {

  E <- c(-sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0))
  W <- c((5.0 / 9.0), (8.0 / 9.0), (5.0 / 9.0))

  # Initialize output (same shape as grid)
  integral <- lat * 0
  
  # Single interval [-1, 1]
  ti <- -1
  tf <- 1
  deltat <- tf - ti
  jacob <- deltat / 2
  
  w <- jacob * W
  t <- jacob * E + ((tf + ti) / 2)

  for (n in seq_along(w)) {
    
    cza <- v_calc_cza(lat = lat,
                    lon = lon,
                    y = y,
                    mon = mon,
                    d = d,
                    hr = hr + t[n])
    integral <- integral + w[n] * cza
    
  }

  return(integral / 2)
}

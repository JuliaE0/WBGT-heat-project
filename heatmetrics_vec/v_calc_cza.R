#' @title Calculate the cosine of the solar zenith angle
#'
#' @description This function can be used by itself, but a more-accurate cza will
#' be obtained with calc_cza_int(), which calls this function and integrates over the hour.
#'
#' @param lat Degrees north latitude (-90 to 90)
#' @param lon Degrees east longitude (-180 to 180)
#' @param y Year (four digits, e.g., 2020)
#' @param mon Month (1-12)
#' @param d Day of month (whole number)
#' @param hr Hour (0-24 UTC)
#'
#' @return Returns cosine of the solar zenith angle (cza)
#' @export
#'
#' @examples v_calc_cza(30, -100, 2020, 1, 1, 12)
v_calc_cza <- function(lat, lon, y, mon, d, hr) {

  date0 <- as.Date(sprintf("%04d-%02d-%02d", y, mon, d))
  datetime <- as.POSIXct(date0, tz = "UTC") + hr * 3600

  # Calculate Julian Day
  jd <- as.integer(strftime(datetime, "%j"))

  hr <- ifelse(hr < 0, 24 + hr, hr)

  # declination angle + time correction for solar angle
  d_tc <- calc_solarDA(jd, hr)
  d <- d_tc$d
  tc <- d_tc$tc

  d_rad <- d * (pi / 180)

  lat_rad <- lat * (pi / 180)

  sindec_sinlat <- sin(d_rad) * sin(lat_rad)
  cosdec_coslat <- cos(d_rad) * cos(lat_rad)

  # solar hour angle [h.deg]
  sha_rad <- ((hr - 12) * 15 + lon + tc) * (pi / 180)
  csza <- sindec_sinlat + cosdec_coslat * cos(sha_rad)

  csza <- pmax(csza, 0)
  return(csza)
}

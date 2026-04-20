#' @title Wet-Bulb Globe Temperature (WBGT)
#'
#' @description Calculates the outdoor wet bulb-globe temperature (WBGT), which is the
#' weighted sum of the dry-bulb air temperature (Ta), the globe temperature (Tg), and
#' the natural wet bulb temperature (Tw):
#'
#' WBGT = (0.1 ⋅ Ta) + (0.7 ⋅ Tw) + (0.2 ⋅ Tg)
#'
#' The program predicts Tw and Tg using meteorological input data, and then combines
#' the results to produce WBGT.
#'
#' Reference: Liljegren, et al. Modeling the Wet Bulb Globe Temperature Using
#' Standard Meteorological Measurements. J. Occup. Environ. Hyg. 5, 645-655 (2008).
#' https://doi.org/10.1080/15459620802310770
#'
#' @param year 4-digit integer, e.g., 2007
#' @param month Month (1-12) or month = 0 if reporting day as day of year
#' @param dday Decimal day of month (1-31.96) -or- day of year (1-366.96), in UTC day-fractions
#' @param lat Degrees north latitude (-90 to 90)
#' @param lon Degrees east longitude (-180 to 180)
#' @param solar Solar irradiance (W/m2)
#' @param cza Cosine solar zenith angle (0-1); use calc_cza_int() or calc_solar_parameters()$cza if cza is not known
#' @param fdir Fraction of surface solar radiation that is direct (0-1)
#' @param pres Barometric pressure in millibars (equivalent to hPa)
#' @param Tair Dry-bulb air temperature (deg. C)
#' @param relhum Relative humidity (\%)
#' @param speed Wind speed (m/s)
#' @param zspeed Height of wind-speed measurement, meters (typically 10m)
#' @param dT Vertical temperature difference (upper minus lower) in degrees Celsius
#' @param urban 1 for urban locations or 0 for non-urban locations
#'
#' @return Returns the wet-bulb globe temperature in degrees C.
#' @export
#'
v_wbgt <- function(year, month, dday, lat, lon, solar, cza, fdir, pres, Tair,
                   relhum, speed, zspeed, dT, urban) {
  
  # --- solar adjustment ---
  solar_out <- v_calc_solar_parameters(
    year, month, dday, lat, lon, solar, cza, fdir
  )
  solar <- solar_out$solarRet

  # --- NA mask ---
  valid <- !is.na(solar) &
    !is.na(cza) &
    !is.na(fdir) &
    !is.na(pres) &
    !is.na(Tair) &
    !is.na(relhum) &
    !is.na(speed) &
    !is.na(urban)
  
  if (!any(valid)) {
    return(array(NA_real_, dim = dim(solar)))
  }
  
  # --- stability ---
  daytime <- cza > 0
  
  stability_class <- v_stab_srdt(
    daytime,
    speed,
    solar,
    dT
  )
  
  # --- wind adjustment ---
  speed_adj <- v_est_wind_speed(
    speed,
    zspeed,
    stability_class,
    urban
  )
  
  speed_adj <- pmax(speed_adj, 0.5)
  
  # --- conversions ---
  tk <- Tair + 273.15
  rh <- 0.01 * relhum
  
  # --- WBGT components ---
  Tg   <- v_Tglobe(tk, rh, pres, speed_adj, solar, fdir, cza)
  Tnwb <- v_Twb(tk, rh, pres, speed_adj, solar, fdir, cza)
  
  Twbg <- (0.1 * Tair) + (0.2 * Tg) + (0.7 * Tnwb)
  
  Twbg[!valid] <- NA
  
  return(Twbg)
}

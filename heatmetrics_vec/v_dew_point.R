#' @title Calculate dew-point temperature from pressure
#'
#' @description Calculates the dew point or frost point temperature
#' in units of Kelvin (K) from barometric pressure and vapor pressure. To calculate
#' dew-point temperature from ambient temperature and relative humidity, use
#' the td() function.
#'
#' @param e Vapor pressure in millibars (equivalent to hPa)
#' @param phase Indicator - 0 for dew point or 1 for frost point
#' @param Pair Barometric pressure in millibars (equivalent to hPa)
#'
#' @return Returns the dew-point temperature in units of Kelvin (K)
#' @export
#'
v_dew_point <- function(e, phase, Pair){

  tdk <- array(NA_real_, dim = dim(e))
  
  idx_water <- phase == 0
  idx_ice   <- !idx_water
  
  # --- water (dew point) ---
  if (any(idx_water)) {
    EF <- 1.0007 + (3.46e-6 * Pair[idx_water])
    z  <- log(e[idx_water] / (6.1121 * EF))
    tdk[idx_water] <- 273.15 + 240.97 * z / (17.502 - z)
  }
  
  # --- ice (frost point) ---
  if (any(idx_ice)) {
    EF <- 1.0003 + (4.18e-6 * Pair[idx_ice])
    z  <- log(e[idx_ice] / (6.1115 * EF))
    tdk[idx_ice] <- 273.15 + 272.55 * z / (22.452 - z)
  }
  
  return(tdk)
}

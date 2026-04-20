#' @title Saturation vapor pressure
#'
#' @description Calculates the saturation vapor pressure (mb) over liquid water
#'  (phase = 0) or ice (phase = 1).
#'
#' @param tk Air temperature in Kelvin (K)
#' @param phase Over liquid water (0) or ice (1)
#' @param Pair Barometric pressure in millibars (equivalent to hPa)
#'
#' @return Returns the saturation vapor pressure in millibars (equivalent
#' to hPa).
#'
v_esat <- function(tk, phase, Pair){
  
  es <- matrix(NA_real_, nrow = nrow(tk), ncol = ncol(tk))
  
  idx_water <- phase == 0
  idx_ice   <- !idx_water
  
  if (any(idx_water)) {
    y <- (tk[idx_water] - 273.15) / (tk[idx_water] - 32.18)
    es[idx_water] <- 6.1121 * exp(17.502 * y)
    es[idx_water] <- (1.0007 + 3.46e-6 * Pair[idx_water]) * es[idx_water]
  }
  
  if (any(idx_ice)) {
    y <- (tk[idx_ice] - 273.15) / (tk[idx_ice] - 0.6)
    es[idx_ice] <- 6.1115 * exp(22.452 * y)
    es[idx_ice] <- (1.0003 + 4.18e-6 * Pair[idx_ice]) * es[idx_ice]
  }
  
  return(es)
}

#' @title Atmospheric emissivity
#'
#' @description Calculates the atmospheric emissivity, a necessary input to the
#' calculation of globe temperature.
#'
#' @param Tair Air temperature in Kelvin (K)
#' @param rh Relative humidity as a proportion between 0 and 1
#' @param pres Barometric pressure in millibars (equivalent to hPa)
#'
#' @return Returns the atmospheric emissivity
#' @export
#'
v_emis_atm <- function(Tair, rh,pres){

  eee = rh * v_esat(Tair, 0, pres)
  return(0.575 * (eee^0.143))

}

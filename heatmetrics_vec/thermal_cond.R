#' @title Thermal conductivity
#'
#' @description Calculates the thermal conductivity of air in units of W/(m⋅K).
#' This value is used as an input to the calculation of the convective heat
#' transfer coefficient.
#'
#' @param Tair Air temperature (K)
#'
#' @return Returns the thermal conductivity in units of W/(m⋅K)
#' @export
#'
#' @examples thermal_cond(290)
thermal_cond <- function(Tair){

  # CONSTANTS ________________________________________________________________
     Cp <- 1003.5
  R_GAS <- 8314.34
  M_AIR <- 28.97
  R_AIR <- (R_GAS / M_AIR)

  return((Cp + 1.25 * R_AIR) * viscosity(Tair))
}

#' @title Heat of evaporation
#'
#' @description Calculates the heat of evaporation in units of J/kg.
#' Also known as the heat of vaporization or enthalpy of vaporization.
#'
#' @param Tair Air temperature in Kelvin (K)
#'
#' @return Returns the heat of evaporation in J/kg.
#' @export
#'
#' @examples evap(293)
evap <- function(Tair){

  # Use an updated algorithm for calculating heat of vaporization
  # Reference: Meyra et al. (2004), https://doi.org/10.1016/j.fluid.2003.12.011

     Zc <- 0.292   # Universal critical ratio
     Tc <- 647.3   # Critical temperature of H2O (K)
     Tt <- 273.16  # Triple temperature of H2O (K)
  dH_tp <- 2500900 # enthalpy of vaporization of H2O at its triple point (J/kg)

  H <- dH_tp * ((Tc - Tair) / (Tc - Tt)) ^ ((Zc * Zc) * ((Tair - Tt) / (Tc - Tt)) + Zc)
  return(H)

  # Original equation used, based on Van Wylen and Sonntag, Table A.1.1
  # return((313.15 - Tair) / 30 * (-71100) + 2.4073e6)

}

#' @title Convective heat transfer coefficient (sphere)
#'
#' @description To calculate the convective heat transfer coefficient
#' in units of W/(m2⋅K) for flow around a sphere.
#'
#' @param diameter Sphere diameter (m)
#' @param Tair Air temperature (K)
#' @param Pair Barometric pressure in millibars (equivalent to hPa)
#' @param speed Wind speed (m/s)
#'
#' @return Returns the convective heat transfer coefficient in units of W/(m2⋅K)
#' @export
#'
v_h_sphere_in_air <- function(diameter, Tair, Pair, speed){

  # CONSTANTS ___________________________________________________________________
      R_GAS <- 8314.34
      M_AIR <- 28.97
      R_AIR <- (R_GAS / M_AIR)
  MIN_SPEED <- 0.5   # originally was 0.13 m/s
         Cp <- 1003.5
         Pr <- (Cp / (Cp + 1.25 * R_AIR))

  density <- Pair * 100 / ( R_AIR * Tair )

  # Calculate Reynolds Number (Re)
  speed_eff <- pmax(speed, MIN_SPEED)
  Re <- speed_eff * density * diameter / viscosity(Tair)

  # Calculate Nusselt Number (Nu)
  Nu <- 2.0 + 0.6 * sqrt(Re) * (Pr^0.3333)

  return (Nu * thermal_cond(Tair) / diameter)
}

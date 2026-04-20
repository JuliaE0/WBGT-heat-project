#' @title Convective heat transfer coefficient (cylinder)
#'
#' @description Calculates the convective heat transfer coefficient in
#' units of W/(m2⋅K) for a long cylinder in cross flow.
#'
#' @param diameter Cylinder diameter (m)
#' @param length Cylinder length (m)
#' @param Tair Air temperature (K)
#' @param Pair Barometric pressure in millibars (equivalent to hPa)
#' @param speed Wind speed (m/s)
#'
#' @return Returns the convective heat transfer coefficient in units of W/(m2⋅K)
#' @export
#'
v_h_cylinder_in_air <- function(diameter, length, Tair, Pair, speed) {

  # CONSTANTS ___________________________________________________________________
          a <- 0.56
          b <- 0.281
          c <- 0.4
      R_GAS <- 8314.34
      M_AIR <- 28.97
      R_AIR <- (R_GAS / M_AIR)
  MIN_SPEED <- 0.5       # Originally was 0.13 m/s
         Cp <- 1003.5
         Pr <- (Cp / (Cp + 1.25 * R_AIR))

  density <- Pair * 100 / (R_AIR * Tair)
  Re <- pmax(speed, MIN_SPEED) * density * diameter / viscosity(Tair)
  Nu <- b * (Re^(1 - c)) * (Pr^(1 - a))
  return (Nu * thermal_cond(Tair) / diameter)
}

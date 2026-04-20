#' @title Diffusivity of water vapor in air
#'
#' @description To calculate the diffusivity of water vapor in air, m2/s.
#'
#' @param Tair Air temperature in Kelvin (K)
#' @param Pair Barometric pressure in millibars (equivalent to hPa)
#'
#' @return Returns the diffusivity in units of m2/s
#' @export
#'
#' @examples diffusivity(290, 1014)
diffusivity <- function(Tair,Pair) {

  # CONSTANTS ________________________________________________________________
      M_AIR <- 28.97
      M_H2O <- 18.015
  Pcrit_air <- 36.4
  Pcrit_h2o <- 218.0
  Tcrit_air <- 132.0
  Tcrit_h2o <- 647.3
          a <- 3.640e-4
          b <- 2.334

  Pcrit13  <- (Pcrit_air * Pcrit_h2o)^(1/3)
  Tcrit512 <- (Tcrit_air * Tcrit_h2o)^(5/12)
  Tcrit12  <- sqrt(Tcrit_air * Tcrit_h2o)
  Mmix <- sqrt(1 / M_AIR + 1 / M_H2O)
  Patm <- Pair / 1013.25  # convert pressure from mb (or hPa) to atmospheres (atm)

  return(a * ((Tair / Tcrit12)^b) * Pcrit13 * Tcrit512 * Mmix / Patm * 1e-4)
}

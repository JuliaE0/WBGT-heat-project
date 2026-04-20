#' @title Viscosity
#'
#' @description Calculates the viscosity of air in units of kg/(m⋅s). This is
#' an input into the calculation of thermal conductivity, wet-bulb temperature,
#' and the convective heat transfer coefficient.
#'
#' @param Tair Air temperature (K)
#'
#' @return Returns the air viscosity in units of kg/(m⋅s)
#' @export
#'
#' @examples viscosity(290)
viscosity <- function(Tair){

  # CONSTANTS ________________________________________________________________
      M_AIR <- 28.97
      sigma <- 3.617
  eps_kappa <- 97.0

  Tr <- Tair / eps_kappa
  omega <- (Tr - 2.9) / 0.4 * (-0.034) + 1.048
  return (2.6693e-6 * sqrt(M_AIR * Tair) / (sigma * sigma * omega))
}

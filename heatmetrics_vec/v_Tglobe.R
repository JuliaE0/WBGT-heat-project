#' @title Globe temperature
#'
#' @description Calculates the globe temperature as an input to the wet-bulb
#' globe temperature (WBGT).
#'
#' @param Tair Dry-bulb air temperature (Kelvin)
#' @param rh Relative humidity as proportion (0-1)
#' @param Pair Barometric pressure in millibars (equivalent to hPa)
#' @param speed Wind speed (m/s)
#' @param solar Solar irradiance (W/m2)
#' @param fdir Fraction of solar irradiance due to direct beam (0-1)
#' @param cza Cosine of solar zenith angle (0-1)
#'
#' @return Returns the globe temperature in degrees Celsius
#' @export
#'
v_Tglobe <- function(Tair, rh, Pair, speed, solar, fdir, cza) {

  # The equation for Tglobe_new has cza in the denominator, so it will result in
  # NaN for cza = 0. This should only be 0 at nighttime, in which case both fdir
  # and cza should both be zero. When fdir and cza are both 0, the value of cza
  # has no bearing on Tglobe, even when solar > 0. To avoid an
  # unnecessary NA value, replace cza with 0.01 when cza < 0.01

  cza_safe <- pmax(cza, 0.01)

  # CONSTANTS ______________________________________________________________________
   EMIS_GLOBE <- 0.95
    ALB_GLOBE <- 0.05
      ALB_SFC <- 0.45
      D_GLOBE <- 0.0508
     EMIS_SFC <- 0.999
      STEFANB <- 5.6696e-8
  CONVERGENCE <- 0.02
     MAX_ITER <- 50

  # VARIABLES ______________________________________________________________________
          Tsfc <- Tair
   Tglobe_prev <- Tair # first guess is the air temperature
   
  # --- active mask ---
  active <- !is.na(Tair) & !is.na(speed) & !is.na(Pair) &
   !is.na(rh) & !is.na(solar) & !is.na(fdir)
  
  if (!any(active)) {
   return(matrix(NA_real_, nrow = nrow(Tair), ncol = ncol(Tair)))
  }
  
  emis_air <- v_emis_atm(Tair, rh, Pair)
  
  converged <- matrix(FALSE, nrow = nrow(Tair), ncol = ncol(Tair))
  
  # --- iteration ---
  for (iter in 1:MAX_ITER) {
    
    # only calculate where needed
    idx <- active & !converged
    
    if (!any(idx)) break
    
    Tref <- 0.5 * (Tglobe_prev + Tair)
    
    h <- v_h_sphere_in_air(D_GLOBE, Tref, Pair, speed)
    
    Tglobe_new <- (
      0.5 * (emis_air * (Tair^4) + EMIS_SFC * (Tsfc^4)) -
        h / (STEFANB * EMIS_GLOBE) * (Tglobe_prev - Tair) +
        solar/(2 * STEFANB * EMIS_GLOBE) * (1 - ALB_GLOBE) *
        (fdir * (1 / (2 * cza_safe) - 1) + 1 + ALB_SFC)
    )^0.25
    
    # --- convergence ---
    diff <- abs(Tglobe_new - Tglobe_prev)
    
    newly_converged <- diff < CONVERGENCE
    newly_converged[is.na(newly_converged)] <- FALSE
    
    converged <- converged | newly_converged
    
    # update only unconverged valid cells
    update_idx <- idx & !newly_converged & !is.na(Tglobe_new)
    
    Tglobe_prev[update_idx] <- 
      0.9 * Tglobe_prev[update_idx] + 0.1 * Tglobe_new[update_idx]
  }
  
  # --- output ---
  result <- Tglobe_new - 273.15
  
  # cells that never converged → NA
  result[!converged] <- NA
  
  return(result)
}

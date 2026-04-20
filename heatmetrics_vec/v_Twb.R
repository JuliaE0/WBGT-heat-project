#' @title Natural wet-bulb temperature
#'
#' @description Calculates the natural wet-bulb temperature.
#'
#' @param Tair Air temperature (dry bulb) in Kelvin (K)
#' @param rh Relative humidity as a proportion (0-1)
#' @param Pair Barometric pressure in millibars	(equivalent to hPa)
#' @param speed Wind speed (m/s)
#' @param solar Solar irradiance (W/m2)
#' @param fdir Fraction of solar irradiance due to direct beam
#' @param cza Cosine of solar zenith angle
#'
#' @return Returns the natural wet-bulb temperature in degrees Celsius (°C)
#' @export
#'
v_Twb <- function(Tair, rh, Pair, speed, solar, fdir, cza){
  
  rad <- 1
  
  # --- safe cza ---
  cza_safe <- pmax(cza, 0.01)
  
  # --- constants ---
  CONVERGENCE <- 0.02
  MAX_ITER <- 30
  D_WICK <- 0.007
  L_WICK <- 0.0254
  PI <- 3.1415926535897932
  Cp <- 1003.5
  R_GAS <- 8314.34
  M_AIR <- 28.97
  M_H2O <- 18.015
  R_AIR <- (R_GAS / M_AIR)
  Pr <- (Cp / (Cp + 1.25 * R_AIR))
  STEFANB <- 5.6696e-8
  EMIS_SFC <- 0.999
  EMIS_WICK <- 0.95
  ALB_WICK <- 0.4
  ALB_SFC <- 0.45
  RATIO <- (Cp * M_AIR / M_H2O)
  a <- 0.56
  
  # --- variables ---
  Tsfc <- Tair
  sza <- acos(cza_safe)
  
  eair <- rh * v_esat(Tair, 0, Pair)
  Tdew <- v_dew_point(eair, 0, Pair)
  
  Twb_prev <- Tdew
  
  # valid cells onl
  active <- !is.na(Tair) & !is.na(speed) & !is.na(Pair) &
    !is.na(rh) & !is.na(solar) & !is.na(fdir)
  
  if (!any(active)) {
    return(matrix(NA_real_, nrow = nrow(Tair), ncol = ncol(Tair)))
  }
  
  emis_air <- v_emis_atm(Tair, rh, Pair)
  
  converged <- array(FALSE, dim = dim(Tair))
  
  # --- track convergence ---
  converged <- array(FALSE, dim = dim(Tair))
  
  for (iter in 1:MAX_ITER) {
    
    idx <- active & !converged
    
    if (!any(idx)) break
    
    Tref <- 0.5 * (Twb_prev + Tair)
    
    h <- v_h_cylinder_in_air(D_WICK, L_WICK, Tref, Pair, speed)
    
    Fatm <- STEFANB * EMIS_WICK *
      (0.5 * (emis_air * (Tair^4) + EMIS_SFC * (Tsfc^4)) -
         (Twb_prev^4)) +
      (1 - ALB_WICK) * solar *
      ((1 - fdir) * (1 + 0.25 * D_WICK / L_WICK) +
         fdir * ((tan(sza) / PI) + 0.25 * D_WICK / L_WICK) +
         ALB_SFC)
    
    ewick <- v_esat(Twb_prev, 0, Pair)
    
    density <- Pair * 100 / (R_AIR * Tref)
    
    Sc <- viscosity(Tref) / (density * diffusivity(Tref, Pair))
    
    Twb_new <- Tair -
      evap(Tref) / RATIO *
      (ewick - eair) / (Pair - ewick) *
      ((Pr / Sc)^a) +
      (Fatm / h * rad)
    
    # --- convergence ---
    diff <- abs(Twb_new - Twb_prev)
    
    newly_converged <- diff < CONVERGENCE
    newly_converged[is.na(newly_converged)] <- FALSE
    
    converged <- converged | newly_converged
    
    # update only unconverged + valid
    update_idx <- idx & !newly_converged & !is.na(Twb_new)
    
    Twb_prev[update_idx] <- 
      0.9 * Twb_prev[update_idx] + 0.1 * Twb_new[update_idx]
  }
  
  Twb_out <- Twb_prev - 273.15
  
  Twb_out[!converged] <- NA
  
  return(Twb_out)
}

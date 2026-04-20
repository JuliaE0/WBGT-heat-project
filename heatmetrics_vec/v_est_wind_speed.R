#' @title Estimate 2-meter wind speeds
#'
#' @description Estimates 2-m wind speeds for all stability conditions when wind speeds are known
#' at a different altitude. The minimum wind speed is set to 0.5 m/s.
#'
#' @param speed Wind speed in meters per second (m/s)
#' @param zspeed Height of wind-speed measurement in meters (typically 10m)
#' @param stability_class Stability class (0-6), defined by the stab_srdt() function
#' @param urban Designation for whether area is urban (1) or not urban (0)
#'
#' @return Returns the estimated 2-meter wind speed in meters per second (m/s)
#' @export
#'
v_est_wind_speed <- function(speed, zspeed, stability_class, urban){

  # Define constants:
  MIN_SPEED <- 0.5   # 0.13 m/s in the original code
  REF_HEIGHT <- 2.0

  urban_exp <- c(0.15, 0.15, 0.20, 0.25, 0.30, 0.30)
  rural_exp <- c(0.07, 0.07, 0.10, 0.15, 0.35, 0.55)

  # --- exponent lookup ---
  exponent <- matrix(NA, nrow = nrow(speed), ncol = ncol(speed))
  
  # --- valid mask (CRITICAL FIX) ---
  valid <- !is.na(stability_class) & stability_class > 0
  
  # urban
  idx_urban <- valid & urban == 1
  exponent[idx_urban] <- urban_exp[stability_class[idx_urban]]
  
  # rural
  idx_rural <- valid & urban != 1
  exponent[idx_rural] <- rural_exp[stability_class[idx_rural]]
  
  # --- wind scaling ---
  est_speed <- speed * ((REF_HEIGHT / zspeed)^exponent)
  
  # enforce minimum
  est_speed <- pmax(est_speed, MIN_SPEED)
  
  return(est_speed)
}

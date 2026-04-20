#' @title Stability class
#'
#' @description Estimates the stability class for calculating 2-meter wind
#' speed from higher-altitude wind speeds.
#'
#' @param daytime "1" for daytime, "0" for nighttime
#' @param speed Wind speed (m/s)
#' @param solar Irradiance (W/m2)
#' @param dT Temperature differential between wind-speed heights (deg C)
#'
#' @return Returns the stability class (0-6) used for adjusting wind speeds from
#' reference height to 2-meter height.
#' @export
#'
v_stab_srdt <- function(daytime, speed, solar, dT){

  lsrdt <- matrix(nrow = 6, ncol = 8)
  lsrdt[1,] <- c(1, 1, 2, 4, 0, 5, 6, 0)
  lsrdt[2,] <- c(1, 2, 3, 4, 0, 4, 5, 0)   # CORRECTED columns 6 & 7 from "5, 6" to "4, 5"
  lsrdt[3,] <- c(2, 2, 3, 4, 0, 4, 4, 0)
  lsrdt[4,] <- c(3, 3, 4, 4, 0, 0, 0, 0)
  lsrdt[5,] <- c(3, 4, 4, 4, 0, 0, 0, 0)
  lsrdt[6,] <- 0
  
  # initialize matrices
  i <- matrix(NA_integer_, nrow = nrow(speed), ncol = ncol(speed))
  j <- matrix(NA_integer_, nrow = nrow(speed), ncol = ncol(speed))
  
  # DAYTIME
  idx_day <- daytime == 1
  
  # j (solar class)
  j[idx_day & solar >= 925] <- 1
  j[idx_day & solar >= 675 & solar < 925] <- 2
  j[idx_day & solar >= 175 & solar < 675] <- 3
  j[idx_day & solar < 175] <- 4
  
  # i (wind class)
  i[idx_day & speed >= 6] <- 5
  i[idx_day & speed >= 5 & speed < 6] <- 4
  i[idx_day & speed >= 3 & speed < 5] <- 3
  i[idx_day & speed >= 2 & speed < 3] <- 2
  i[idx_day & speed < 2] <- 1
  
  # NIGHTTIME
  idx_night <- !idx_day
  
  j[idx_night & dT >= 0] <- 7
  j[idx_night & dT < 0]  <- 6
  
  i[idx_night & speed >= 2.5] <- 3
  i[idx_night & speed >= 2 & speed < 2.5] <- 2
  i[idx_night & speed < 2] <- 1
  
  # --- lookup ---
  out <- matrix(NA, nrow = nrow(i), ncol = ncol(i))
  
  valid <- !is.na(i) & !is.na(j)
  out[valid] <- lsrdt[cbind(i[valid], j[valid])]
  
  return(out)
}

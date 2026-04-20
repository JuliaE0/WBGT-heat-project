#' @title Calculate solar parameters
#'
#' @description To calculate the adjusted surface solar irradiance
#' and fraction of the solar irradiance due to the direct beam. Note that
#' calc_cza_int() provides more-accurate cza on hourly data and should be used when possible.
#'
#' @param year Year (4 digits)
#' @param month Month (1-12)
#' @param day Day-fraction of month based on UTC time. Day number must include
#' fractional day based on time, e.g., 4.5 = noon UTC on the 4th of the month.
#' @param lat Degrees north latitude (-90 to 90)
#' @param lon Degrees east longitude (-180 to 180)
#' @param solar Total surface solar irradiance (W/m2)
#' @param cza Cosine solar zenith angle (0-1); optional (supply "NA" if unknown)
#' @param fdir Fraction of the surface solar radiation from direct (0-1); optional (supply "NA" if unknown)
#'
#' @return Returns adjusted solar radiation ("solarRet") and the fraction of irradiance due to
#' direct beam ("fdir", unchanged if user-supplied).
#'
v_calc_solar_parameters <- function(year, month, day, lat, lon, solar, cza, fdir) {
  
  SOLAR_CONST <- 1367.0
  CZA_MIN <- 0.00873
  NORMSOLAR_MAX <- 0.85
  
  # --- Earth–Sun distance ---
  day_int  <- floor(day)
  frac_day <- day - day_int
  hour_full <- frac_day * 24     
  
  jd <- as.integer(strftime(
    as.POSIXct(sprintf("%04d-%02d-%02d", year, month, day_int), tz="UTC") + hour_full*3600,
    "%j"
  ))
  
  g <- 2 * pi * (jd - 1 + frac_day) / 365.25
  soldist <- 1.00014 - 0.01671 * cos(g) - 0.00014 * cos(2*g)
  
  dims <- dim(solar)
  
  # --- initialize ---
  solarRet <- solar
  fdir_out <- fdir
  
  soldist <- array(soldist, dim = dims)
  
  # --- clamp cza ---
  cza_clamped <- ifelse(!is.na(cza), pmax(cza, 0), NA)
  
  # --- TOA solar ---
  toasolar <- SOLAR_CONST * cza_clamped / (soldist^2)
  
  # horizon cutoff
  toasolar[cza_clamped < CZA_MIN] <- 0
  
  idx_day <- !is.na(toasolar) & toasolar > 0
  
  # --- normalized solar ---
  normsolar <- array(NA_real_, dim = dims)
  normsolar[idx_day] <- pmin(solar[idx_day] / toasolar[idx_day], NORMSOLAR_MAX)
  
  # --- solarRet update ONLY where valid ---
  solarRet[idx_day] <- normsolar[idx_day] * toasolar[idx_day]
  
  # --- preserve solar NA ---
  solarRet[is.na(solar)] <- NA
  
  # --- fdir ---
  idx_norm_pos <- idx_day & !is.na(normsolar) & normsolar > 0
  
  # estimate only where NA
  idx_est <- is.na(fdir_out) & idx_norm_pos
  fdir_out[idx_est] <- exp(3 - 1.34 * normsolar[idx_est] - 1.65 / normsolar[idx_est])
  
  # clamp
  valid_fdir <- !is.na(fdir_out)
  fdir_out[valid_fdir] <- pmax(pmin(fdir_out[valid_fdir], 0.9), 0)
  
  # --- night / invalid conditions ---
  idx_zero <- !idx_norm_pos
  fdir_out[idx_zero] <- 0
  
  # --- preserve original NA in fdir ONLY where truly missing ---
  fdir_out[is.na(fdir)] <- NA
  
  return(list(solarRet = solarRet, fdir = fdir_out))
}

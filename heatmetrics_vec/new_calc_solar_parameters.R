new_calc_solar_parameters <- function(year, month, day, lat, lon, solar, cza, fdir) {
  
  # --- constants ---
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
  
  # --- ensure matrix shape ---
  dims <- dim(solar)
  
  solarRet  <- array(NA_real_, dim = dims)
  normsolar <- array(NA_real_, dim = dims)
  toasolar  <- array(NA_real_, dim = dims)
  
  # broadcast soldist
  soldist <- array(soldist, dim = dims)
  
  # --- masks ---
  valid <- !is.na(cza) & !is.na(solar)
  
  # --- clamp cza ---
  cza_clamped <- ifelse(!is.na(cza), pmax(cza, 0), NA)
  
  # --- TOA solar ---
  toasolar[valid] <- SOLAR_CONST * cza_clamped[valid] / (soldist[valid]^2)
  
  # horizon cutoff
  idx_horizon <- valid & !is.na(cza_clamped) & cza_clamped < CZA_MIN
  toasolar[idx_horizon] <- 0
  
  # --- normalized solar ---
  idx <- valid & !is.na(toasolar) & toasolar > 0
  
  normsolar[idx] <- pmin(solar[idx] / toasolar[idx], NORMSOLAR_MAX)
  
  # --- adjusted solar ---
  solarRet[idx] <- normsolar[idx] * toasolar[idx]
  
  # --- fdir ---
  fdir_out <- fdir
  
  idx_fdir <- is.na(fdir) & !is.na(normsolar) & normsolar > 0
  
  fdir_out[idx_fdir] <- exp(3 - 1.34 * normsolar[idx_fdir] - 1.65 / normsolar[idx_fdir])
  
  # clamp
  not_na_fdir <- !is.na(fdir_out)
  fdir_out[not_na_fdir] <- pmax(pmin(fdir_out[not_na_fdir], 0.9), 0)
  
  # zero conditions
  idx_zero <- !is.na(toasolar) & (
    toasolar <= 0 | (!is.na(normsolar) & normsolar <= 0)
  )
  fdir_out[idx_zero] <- 0
  
  return(list(solarRet = solarRet, fdir = fdir_out))
}
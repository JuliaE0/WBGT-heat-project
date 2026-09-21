# Calculate WBGT for 2000 to 2025


library(ncdf4)
library(terra)

source("/mnt/HDD7/juliae/WBGT-heat-project/heatmetrics_vec/v_calc_solar_parameters.R")
source("/mnt/HDD7/juliae/WBGT-heat-project/heatmetrics_vec/v_stab_srdt.R")
source("/mnt/HDD7/juliae/WBGT-heat-project/heatmetrics_vec/v_est_wind_speed.R")
source("/mnt/HDD7/juliae/WBGT-heat-project/heatmetrics_vec/v_esat.R")
source("/mnt/HDD7/juliae/WBGT-heat-project/heatmetrics_vec/v_emis_atm.R")
source("/mnt/HDD7/juliae/WBGT-heat-project/heatmetrics_vec/v_h_sphere_in_air.R")
source("/mnt/HDD7/juliae/WBGT-heat-project/heatmetrics_vec/v_Tglobe.R")
source("/mnt/HDD7/juliae/WBGT-heat-project/heatmetrics_vec/v_dew_point.R")
source("/mnt/HDD7/juliae/WBGT-heat-project/heatmetrics_vec/v_h_cylinder_in_air.R")
source("/mnt/HDD7/juliae/WBGT-heat-project/heatmetrics_vec/v_Twb.R")
source("/mnt/HDD7/juliae/WBGT-heat-project/heatmetrics_vec/v_wbgt.R")
source("/mnt/HDD7/juliae/WBGT-heat-project/heatmetrics_vec/viscosity.R")
source("/mnt/HDD7/juliae/WBGT-heat-project/heatmetrics_vec/thermal_cond.R")
source("/mnt/HDD7/juliae/WBGT-heat-project/heatmetrics_vec/diffusivity.R")
source("/mnt/HDD7/juliae/WBGT-heat-project/heatmetrics_vec/evap.R")

############## define input and output folders
input_dir <- "wbgt_inputs_complete"
output_dir <- "wbgt"
dir.create(output_dir, showWarnings = FALSE)

############## get all files
files <- list.files(
  input_dir,
  pattern = "wbgt_inputs_\\d{4}_\\d{2}\\.nc$",
  full.names = TRUE)


for (file in files) {
  
  file_name <- basename(file)
  out_name <- sub("wbgt_inputs_", "wbgt_", file_name)
  out_file <- file.path(output_dir, out_name)
  
  # skip existing outputs
  if (file.exists(out_file)) {
    cat("Skipping:", out_name, "\n")
    next}
  
  inputs <- nc_open(file)
  
  file_year <- substr(file_name, 13, 16)
  file_month <- substr(file_name, 18, 19)
  cat("Calculating:", file_year, "-", file_month, "\n")
  
  ############## get input variables
  year  <- ncvar_get(inputs, "year")
  month <- ncvar_get(inputs, "month")
  dday  <- ncvar_get(inputs, "dday")
  lat <- ncvar_get(inputs, "latitude")
  lon <- ncvar_get(inputs, "longitude")
  solar  <- aperm(ncvar_get(inputs, "v_solarRet"), c(3,2,1))
  cza    <- aperm(ncvar_get(inputs, "v_cza"), c(3,2,1))
  fdir   <- aperm(ncvar_get(inputs, "v_fdir"), c(3,2,1))
  pres   <- aperm(ncvar_get(inputs, "pres"), c(3,2,1))
  Tair   <- aperm(ncvar_get(inputs, "Tair"), c(3,2,1))
  relhum <- aperm(ncvar_get(inputs, "relhum"), c(3,2,1))
  speed  <- aperm(ncvar_get(inputs, "speed"), c(3,2,1))
  urban  <- aperm(ncvar_get(inputs, "urban"), c(3,2,1))
  zspeed <- 10
  dT <- -0.052
  
  nlat <- length(lat)
  nlon <- length(lon)
  ntime <- dim(solar)[1]
  
  lat_mat <- matrix(lat, nrow = nlat, ncol = nlon)
  lon_mat <- matrix(lon, nrow = nlat, ncol = nlon, byrow = TRUE)
  
  wbgt_out <- array(NA_real_, dim = c(ntime, nlat, nlon))
  
  ############## WBGT calculation
  system.time({
  for (t in 1:ntime) {
    year_mat  <- matrix(year[t],  nlat, nlon)
    month_mat <- matrix(month[t], nlat, nlon)
    dday_mat  <- matrix(dday[t],  nlat, nlon)
    
    wbgt_out[t,,] <- v_wbgt(
      year_mat, month_mat, dday_mat,
      lat_mat, lon_mat,
      solar[t,,], cza[t,,], fdir[t,,],
      pres[t,,], Tair[t,,], relhum[t,,],
      speed[t,,], zspeed, dT, urban[t,,]
    )}
  })
  
  ########## save output
  wbgt_to_write <- aperm(wbgt_out, c(3,2,1))
  
  # define dimensions
  dim_lon  <- ncdim_def("longitude","degrees_east", lon)
  dim_lat  <- ncdim_def("latitude","degrees_north", lat)
  dim_time <- ncdim_def("valid_time","seconds since 1970-01-01", ncvar_get(inputs, "valid_time"))
  
  # define variable
  wbgt <- ncvar_def(
    name = "wbgt",
    units = "degrees_Celsius",
    dim = list(dim_lon, dim_lat, dim_time),
    missval = -9999,
    longname = "Wet-bulb globe temperature",
    prec = "float")

  # create NetCDF with wbgt variable
  wbgt_output <- nc_create(out_file, list(wbgt))
  ncvar_put(wbgt_output, "wbgt", wbgt_to_write)
  
  nc_close(inputs)
  nc_close(wbgt_output)
  
}

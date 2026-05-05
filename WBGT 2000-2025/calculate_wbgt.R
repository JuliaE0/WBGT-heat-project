# Calculate WBGT for 2000 to 2025


library(ncdf4)

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
  
  cat("Calculating:", basename(file), "\n")
  
  inputs <- nc_open(file)
  
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
    )
  }
  
  ########## crop to California shape
  times <- ncvar_get(inputs, "valid_time")
  
  wbgt_for_terra <- aperm(wbgt_out, c(2,3,1))
  
  r_wbgt <- rast(wbgt_for_terra)
  ext(r_wbgt) <- c(min(lon), max(lon), min(lat), max(lat))
  crs(r_wbgt) <- "EPSG:4326"
  
  ca_shape <- vect("ca_grid/ERA5_Land_grid_CAfull.shp")
  ca_shape <- project(ca_shape, crs(r_wbgt))
  
  r_wbgt <- mask(crop(r_wbgt, ca_shape), ca_shape)
  
  vals <- terra::values(r_wbgt, mat = TRUE)
  
  nrow_ca <- nrow(r_wbgt)
  ncol_ca <- ncol(r_wbgt)
  ntime <- nlyr(r_wbgt)
  
  wbgt_ca <- array(vals, dim = c(ncol_ca, nrow_ca, ntime))
  wbgt_ca <- aperm(wbgt_ca, c(3,2,1))
  
  ########## save output
  out_name <- sub("wbgt_inputs_", "wbgt_", basename(file))
  out_file <- file.path(output_dir, out_name)
  
  # define dimensions
  dim_lon  <- ncdim_def("longitude","degrees_east",inputs$dim$longitude$vals)
  dim_lat  <- ncdim_def("latitude","degrees_north",inputs$dim$latitude$vals)
  dim_time <- ncdim_def("valid_time","seconds since 1970-01-01", inputs$dim$valid_time$vals)
  
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
  ncvar_put(wbgt_output, "wbgt", wbgt_ca)
  
  nc_close(inputs)
  nc_close(wbgt_output)
  
}

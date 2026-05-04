# Inputs preprocessing, part 2

library(ncdf4)

############ load vectorized functions
source("C:/Users/new account/Desktop/WBGT-heat-project/heatmetrics_vec/v_calc_cza.R")
source("C:/Users/new account/Desktop/WBGT-heat-project/heatmetrics_vec/v_calc_cza_int.R")
source("C:/Users/new account/Desktop/WBGT-heat-project/heatmetrics_vec/v_calc_solar_parameters.R")

############ open wbgt_inputs_part1 files
files <- list.files("wbgt_inputs_part1",
                    pattern = "wbgt_inputs_\\d{4}_\\d{2}\\.nc$",
                    full.names = TRUE)

############ create output directory
dir.create("wbgt_inputs_complete", showWarnings = FALSE)

############ loop through files
for (file in files) {
  wbgt_inputs <- nc_open(file)  # open file
  
  ############ ----------- calculate cza -----------
  # get inputs for v_calc_cza_int
  lat <- ncvar_get(wbgt_inputs, "latitude")
  lon <- ncvar_get(wbgt_inputs, "longitude")
  year  <- ncvar_get(wbgt_inputs, "year")
  month <- ncvar_get(wbgt_inputs, "month")
  day   <- ncvar_get(wbgt_inputs, "day")
  hour  <- ncvar_get(wbgt_inputs, "hour")
  
  # make grid
  lon_grid <- matrix(rep(lon, each = length(lat)), nrow = length(lat))
  lat_grid <- matrix(rep(lat, times = length(lon)), nrow = length(lat))
  
  # preallocate output array
  ntime <- length(year)
  nlat  <- length(lat)
  nlon  <- length(lon)
  
  v_cza_array <- array(NA, dim = c(nlat, nlon, ntime))
  
  # apply function over time (vectorized function takes grid input and outputs a grid of values)
  cza_list <- lapply(seq_len(ntime), function(t) {
    v_calc_cza_int(
      lat_grid,
      lon_grid,
      year[t],
      month[t],
      day[t],
      hour[t]
    )
  })
  
  v_cza_array <- simplify2array(cza_list)
  v_cza_array <- aperm(v_cza_array, c(2, 1, 3))
  
  ############ ----------- calculate fdir and solarRet -----------
  # get inputs for v_calc_solar_parameters
  lat <- ncvar_get(wbgt_inputs, "latitude")
  lon <- ncvar_get(wbgt_inputs, "longitude")
  year  <- ncvar_get(wbgt_inputs, "year")
  month <- ncvar_get(wbgt_inputs, "month")
  day   <- ncvar_get(wbgt_inputs, "dday")
  solar <- ncvar_get(wbgt_inputs, "ssrd_W")  # ERA5 ssrd in W/m2
  cza <- v_cza_array  # calculated from calc_cza_int()
  fdir <- ncvar_get(wbgt_inputs, "fdir_frac") # (NN-inter ERA5 fdir)/(NN-inter ERA5 ssrd)
  
  # make grid
  lon_grid <- matrix(rep(lon, each = length(lat)), nrow = length(lat))
  lat_grid <- matrix(rep(lat, times = length(lon)), nrow = length(lat))
  
  # preallocate output array
  ntime <- length(year)
  nlat  <- length(lat)
  nlon  <- length(lon)
  
  v_solarRet_array <- array(NA, dim = c(nlat, nlon, ntime))
  v_fdir_array <- array(NA, dim = c(nlat, nlon, ntime))
  
  # apply function over time
  solar_params_list <- lapply(seq_len(ntime), function(t) {
    
    v_calc_solar_parameters(
      year  = year[t],
      month = month[t],
      day   = day[t],
      lat   = lat_grid,
      lon   = lon_grid,
      solar = solar[,,t],
      cza   = cza[,,t],
      fdir  = fdir[,,t]
    )
  })
  
  v_solarRet_array <- simplify2array(lapply(solar_params_list, `[[`, "solarRet"))
  v_fdir_array     <- simplify2array(lapply(solar_params_list, `[[`, "fdir"))
  
  
  ############ ----------- create v_cza, v_fdir, and v_solarRet variables -----------
  # define dimensions
  dim_lon  <- ncdim_def("longitude","degrees_east",wbgt_inputs$dim$longitude$vals)
  dim_lat  <- ncdim_def("latitude","degrees_north",wbgt_inputs$dim$latitude$vals)
  dim_time <- ncdim_def("valid_time","seconds since 1970-01-01", wbgt_inputs$dim$valid_time$vals)
  
  # define variables
  v_cza <- ncvar_def(
    name = "v_cza",
    units = "",
    dim = list(dim_lon, dim_lat, dim_time),
    missval = -9999,
    longname = "cosine solar zenith angle (0-1)",
    prec = "float")
  v_fdir <- ncvar_def(
    name = "v_fdir",
    units = "",   # no units because it's a ratio
    dim = list(dim_lon, dim_lat, dim_time),
    missval = -9999,
    longname = "Fraction of irradiance due to direct beam (0-1)",
    prec = "float")
  v_solarRet <- ncvar_def(
    name = "v_solarRet",
    units = "W m**-2",
    dim = list(dim_lon, dim_lat, dim_time),
    missval = -9999,
    longname = "adjusted solar radiation",
    prec = "float")
  
  ############ ----------- save v_cza, v_fdir, and v_solarRet to input files -----------
  output_file <- file.path("wbgt_inputs_complete", basename(file))
  file.copy(file, output_file, overwrite = TRUE)
  
  nc_out <- nc_open(output_file, write = TRUE)
  ncvar_add(nc_out, v_cza)
  ncvar_add(nc_out, v_fdir)
  ncvar_add(nc_out, v_solarRet)
  
  ncvar_put(nc_out, "v_cza", v_cza_array)
  ncvar_put(nc_out, "v_fdir", v_fdir_array)
  ncvar_put(nc_out, "v_solarRet", v_solarRet_array)
  
  nc_close(nc_out)
  
}
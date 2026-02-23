

# install required libraries
# !pip install xarray netCDF4 h5netcdf numpy rioxarray regionmask geopandas

import xarray as xr
import numpy as np
import geopandas as gpd
import regionmask

# loop over montly files in era5_land_data folder

# save preprocessed files into a new folder

out_dir = Path("era5_land_preprocessed")
out_dir.mkdir(parents=True, exist_ok=True)


# convert t2m from K to Celsius -- used as "Tair" input
t2m_c = era5land["t2m"] - 273.15
t2m_c.attrs["units"] = "C"

# convert d2m from K to Celsius -- used in "relhum" calculation
d2m_c = era5land["d2m"] - 273.15
d2m_c.attrs["units"] = "C"

# convert ssrd from J/m2 to W/m2 -- used as "solar" input
era5land["ssrd"].valid_time.diff("time") # check that time resolution is hourly
ssrd_W = era5land["ssrd"] / 3600         # dividy by 3600 seconds since hourly accumulated
ssrd_W.attrs["units"] = "W m**-2"

# convert sp from Pa to hPa -- used as "pres" input
sp_hpa = era5land["sp"] / 100
sp_hpa.attrs["units"] = "hPa"


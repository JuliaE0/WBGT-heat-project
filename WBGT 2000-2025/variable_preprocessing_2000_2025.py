# ---- Variable preprocessing 2000 to 2025 ----

# install required libraries
# !pip install xarray netCDF4 h5netcdf numpy rioxarray regionmask geopandas

import xarray as xr
import numpy as np
import geopandas as gpd
import regionmask

# create folder to hold new preprocessed files
out_dir = Path("era5_land_preprocessed")
out_dir.mkdir(parents=True, exist_ok=True)

# loop over monthly files in era5_land_data folder
IN_DIR = "/Users/new account/Desktop/WBGT-heat-project/WBGT 2000-2025/era5_land_data"

# CONTINUE HERE

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

# calculate relative humidity (%) -- used as "relhum" input
vapor_pres = 610.94*np.exp(17.625*d2m_c / (243.04+d2m_c))
sat_vapor_pres = 610.94*np.exp(17.625*t2m_c / (243.04+t2m_c))
rh = 100*(vapor_pres/sat_vapor_pres)
rh.attrs["long_name"] = "relative humidity (%)"
rh.attrs["units"] = "%"

# calculate wind speed (m/s) -- used as "speed" input
u = era5land["u10"]
v = era5land["v10"]
ws = np.sqrt(u**2 + v**2)
ws.attrs["long_name"] = "wind speed (m/s)"
ws.attrs["units"] = "m s**-1"



# save preprocessed files into the new folder


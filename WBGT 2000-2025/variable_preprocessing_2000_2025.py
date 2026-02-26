# ---- Variable preprocessing 2000 to 2025 ----

# Preprocess these WBGT input variables: solar, fdir, pres, Tair, relhum, speed, urban
# other WBGT input: cza will be derived in R

# install required libraries
# !pip install xarray netCDF4 h5netcdf numpy rioxarray regionmask geopandas

import xarray as xr
import numpy as np
import geopandas as gpd
import regionmask
from pathlib import Path


INPUT_DIR_ERA5_LAND = Path("era5_land_data")
INPUT_DIR_ERA5 = Path("era5_interpolated")
OUTPUT_DIR = Path("wbgt_inputs_preprocessed")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def preprocess_variables(era5_land, era5):
    """
    preprocess input variables
    """

    # convert t2m from K to Celsius -- used as "Tair" input
    Tair = era5_land["t2m"] - 273.15
    Tair.attrs["units"] = "C"

    # convert d2m from K to Celsius -- used in "relhum" calculation
    d2m_c = era5_land["d2m"] - 273.15
    d2m_c.attrs["units"] = "C"

    # convert ssrd from J/m2 to W/m2 -- used as "solar" input
    solar = era5_land["ssrd"] / 3600    # dividy by 3600 seconds since hourly accumulated
    solar.attrs["units"] = "W m**-2"

    # convert sp from Pa to hPa -- used as "pres" input
    pres = era5_land["sp"] / 100
    pres.attrs["units"] = "hPa"
    
    # calculate relative humidity (%) -- used as "relhum" input
    vapor_pres = 610.94*np.exp(17.625*d2m_c / (243.04+d2m_c))
    sat_vapor_pres = 610.94*np.exp(17.625*Tair / (243.04+Tair))
    relhum = 100*(vapor_pres/sat_vapor_pres)
    relhum.attrs["long_name"] = "relative humidity (%)"
    relhum.attrs["units"] = "%"

    # calculate wind speed (m/s) -- used as "speed" input
    u = era5_land["u10"]
    v = era5_land["v10"]
    speed = np.sqrt(u**2 + v**2)
    speed.attrs["long_name"] = "wind speed (m/s)"
    speed.attrs["units"] = "m s**-1"

    # calculate fraction of surface solar radiation that is direct (0-1) -- used as "fdir" input
    fdir = xr.where(era5_land["ssrd"].isnull(), np.nan, 
                    xr.where(era5_land["ssrd"] > 0, era5["fdir"].values / era5_land["ssrd"].values, 0.0))
    fdir.attrs["long_name"] = "fraction of surface solar radiation that is direct (0-1)"
    fdir.attrs["units"] = ""  # no units since it's a fraction

    # load urban variable that was created in QGIS -- used as "urban" input
    # TO DO
    # first need to create the yearly urban files in QGIS

    
    # build working dataset: save preprocessed variables into output dataset
    preprocessed = xr.Dataset(
        {"solar": solar,
         "fdir": fdir,
         "pres": pres,
         "Tair": Tair,
         "relhum": relhum,
         "speed": speed,
         "urban": TO DO
        }
        )
    preprocessed = preprocessed.sortby("latitude", ascending=False)
    
    # Create needed coordinates
    time = preprocessed.valid_time
    preprocessed = preprocessed.assign_coords(
        year = time.dt.year,
        month = time.dt.month,
        day=time.dt.day,
        hour=time.dt.hour,
        # decimal day of month (UTC)
        dday = (time.dt.day
                + time.dt.hour / 24
                + time.dt.minute / 1440
                + time.dt.second / 86400))
    
    return preprocessed


def clip_to_ca_boundary(preprocessed):
    
    shp_path = Path("ca_state/CA_state.shp")  # USE NEW SHAPE FILE. DOWNLOAD FROM GOOGLE DRIVE FIRST. FIX
    ca = gpd.read_file(shp_path)
    ca = ca.to_crs("EPSG:4326")
    mask = regionmask.mask_geopandas(ca, preprocessed.longitude, preprocessed.latitude)
    output = preprocessed.where(~mask.isnull()) # keep whole grid cells on the edges of CA
    
    return output


def main():
    files = sorted(INPUT_DIR_ERA5_LAND.glob("era5_land_*.nc"))

    for file in files:
        year = file.stem.split("_")[2]
        month = file.stem.split("_")[3]
        
        output_file = OUTPUT_DIR/f"wbgt_inputs_{year}_{month}.nc"
        
        if output_file.exists():
            print(f"Skipping existing: {output_file.name}")
            continue
        
        idw_file = INPUT_DIR_ERA5 / f"era5_{year}_{month}_idw.nc"

        print(f"Processing: {year}-{month}")

        with xr.open_dataset(file) as era5_land, xr.open_dataset(idw_file) as era5:
            ds = preprocess_variables(era5_land, era5)
            ds = clip_to_ca_boundary(ds)
            ds.to_netcdf(output_file)

    print("Done.")


if __name__ == "__main__":
    main()




# ---- Variable preprocessing 2000 to 2025 ----

# install required libraries
# !pip install xarray netCDF4 h5netcdf numpy rioxarray regionmask geopandas

import xarray as xr
import numpy as np
import geopandas as gpd
import regionmask
from pathlib import Path

# set directory and create folder to hold new preprocessed files
INPUT_DIR = Path("era5_land_data")
OUTPUT_DIR = Path("era5_land_preprocessed")
OUTPUT_DIR.mkdir(parents=TRUE, exist_ok=True)

def preprocess_dataset(ds):
    """
    preprocess ERA5-Land variables
    """

    # convert t2m from K to Celsius -- used as "Tair" input
    ds["Tair"] = ds["t2m"] - 273.15
    ds["Tair"].attrs["units"] = "C"

    # convert d2m from K to Celsius -- used in "relhum" calculation
    ds["d2m_c"] = ds["d2m"] - 273.15
    ds["d2m_c"].attrs["units"] = "C"

    # convert ssrd from J/m2 to W/m2 -- used as "solar" input
    ds["solar"] = ds["ssrd"] / 3600    # dividy by 3600 seconds since hourly accumulated
    ds["solar"].attrs["units"] = "W m**-2"

    # convert sp from Pa to hPa -- used as "pres" input
    ds["pres"] = ds["sp"] / 100
    ds["pres"].attrs["units"] = "hPa"
    
    # calculate relative humidity (%) -- used as "relhum" input
    vapor_pres = 610.94*np.exp(17.625*ds["d2m_c"] / (243.04+ds["d2m_c"]))
    sat_vapor_pres = 610.94*np.exp(17.625*ds["Tair"] / (243.04+ds["Tair"]))
    ds["relhum"] = 100*(vapor_pres/sat_vapor_pres)
    ds["relhum"].attrs["long_name"] = "relative humidity (%)"
    ds["relhum"].attrs["units"] = "%"

    # calculate wind speed (m/s) -- used as "speed" input
    u = ds["u10"]
    v = ds["v10"]
    ds["speed"] = np.sqrt(u**2 + v**2)
    ds["speed"].attrs["long_name"] = "wind speed (m/s)"
    ds["speed"].attrs["units"] = "m s**-1"

    # calculate fraction of surface solar radiation that is direct (0-1) -- used as "fdir" input
    # TO DO

    # load urban_variable.nc that was created in QGIS
    # TO DO

    # build working dataset: save preprocessed variables into output file
    preprocessed = xr.Dataset(
        {"solar": ds["solar"],
         "fdir": TODO,
         "pres": ds["pres"],
         "Tair": ds["Tair"],
         "relhum": ds["relhum"],
         "speed": ds["speed"],
         "urban": TODO
        }
        )
    preprocessed = preprocessed.sortby("latitude", ascending=False)
    
    # Create needed coordinates
    time = preprocessed.valid_time
    output = preprocessed.assign_coords(
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
    
    shp_path = Path("ca_state/CA_state.shp")
    ca = gpd.read_file(shp_path)
    ca = ca.to_crs("EPSG:4326")
    mask = regionmask.mask_geopandas(ca, preprocessed.longitude, preprocessed.latitude)
    output = preprocessed.where(~mask.isnull()) # keep whole grid cells on the edges of CA
    
    return output


def main():
    files = sorted(INPUT_DIR.glob("era5_land_*.nc"))

    for file in files:
        year = file.stem.split("_")[2]
        month = file.stem.split("_")[3]
        output_file = OUTPUT_DIR/f"wbgt_inputs_{year}_{month}.nc"

        if output_file.exists():
            print(f"Skipping existing: {output_file.name}")
            continue

        print(f"Processing: {file.name}")

        with xr.open_dataset(file) as ds:
            ds_processed = preprocess_dataset(ds)
            ds_processed = clip_to_ca_boundary(ds_processed)
            ds_processed.to_netcdf(output_file)

    print("Done.")


if __name__ == "__main__":
    main()




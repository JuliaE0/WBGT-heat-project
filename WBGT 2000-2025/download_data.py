# ---- Download data 2000 to 2025 ----

# install required libraries
# !pip install xarray netCDF4 h5netcdf numpy rioxarray regionmask geopandas "cdsapi>=0.7.7"

import cdsapi
from pathlib import Path

# Download ERA5-Land data using CDS API; loop over all months from 2000 to 2025

out_dir = Path("era5_land_data")
out_dir.mkdir(parents=True, exist_ok=True)
client = cdsapi.Client()
dataset = "reanalysis-era5-land"
request = {
    "variable": [
        "2m_dewpoint_temperature",
        "2m_temperature",
        "surface_solar_radiation_downwards",
        "10m_u_component_of_wind",
        "10m_v_component_of_wind",
        "surface_pressure"],
    "day": [f"{d:02d}" for d in range(1, 32)],
    "time": [f"{h:02d}:00" for h in range(0, 24)],
    "data_format": "netcdf",
    "download_format": "unarchived",
    "area": [42, -124.5, 32.5, -114.1]}

for y in range(2000, 2026):
    for m in range(1, 13):
        year = f"{y}"
        month = f"{m:02d}"
        target = out_dir/f"era5_land_{year}_{month}.nc"
        if target.exists():
            print("Skip existing:", target)
            continue
        req = request.copy()
        req["year"] = year
        req["month"] = month
        print("Downloading", target)
        client.retrieve(dataset, req, target)
print("Finished downloading ERA5-Land data")


# Download ERA5 data using CDS API; loop over all months from 2000 to 2025

out_dir = Path("era5_data")
out_dir.mkdir(parents=True, exist_ok=True)
client = cdsapi.Client()
dataset = "reanalysis-era5-single-levels"
request = {
    "product_type": ["reanalysis"],
    "variable": ["total_sky_direct_solar_radiation_at_surface"],
    "day": [f"{d:02d}" for d in range(1, 32)],
    "time": [f"{h:02d}:00" for h in range(0, 24)],
    "data_format": "netcdf",
    "download_format": "unarchived",
    "area": [42, -124.5, 32.5, -114.1]
}

for y in range(2000, 2026):
    for m in range(1, 13):
        year = f"{y}"
        month = f"{m:02d}"
        target = out_dir/f"era5_{year}_{month}.nc"
        if target.exists():
            print("Skip existing:", target)
            continue
        req = request.copy()
        req["year"] = year
        req["month"] = month
        print("Downloading", target)
        client.retrieve(dataset, req, target)
print("Finished downloading ERA5 data")



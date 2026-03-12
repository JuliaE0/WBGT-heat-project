#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Edited for nearest-neighbor interpolation (Spangler-style)

ERA5 variables -> subset -> nearest-neighbor interpolate to target grid defined by GeoTIFF
-> export one NetCDF per variable.

This is a minimal revision of the original IDW script: only the interpolation
method has been changed from IDW to nearest-neighbor, matching the method
described in Spangler et al. for ERA5 FDIR.

Dependencies:
  pip install rasterio pyproj xarray netcdf4 numpy
"""

import os
import glob
import numpy as np
import xarray as xr
import rasterio
from rasterio.transform import xy
from pyproj import Transformer

# =========================
# USER SETTINGS (EDIT THESE)
# =========================
IN_DIR = "/Users/new account/Desktop/WBGT-heat-project/Aug 2019 test"
OUT_DIR = "/Users/new account/Desktop/WBGT-heat-project/Aug 2019 test/nn_interpolated"
TARGET_TIF = "/Users/new account/Desktop/WBGT-heat-project/era5land_d2m_2019-08-01T00.tif"

# Subset buffer in degrees (applied in WGS84 lon/lat after converting bounds)
BUFFER_DEG = 0.20

# Variable names in ERA5
ERA5_VAR = ["fdir"]
# =========================


def read_target_grid_from_tif(tif_path):
    """
    Returns target grid definition derived from GeoTIFF:
      - tgt_lon2d, tgt_lat2d: 2D arrays of cell-center lon/lat in EPSG:4326
      - tgt_lons, tgt_lats: 1D vectors if grid is regular in lon/lat
      - bounds_wgs84: (minlon, minlat, maxlon, maxlat)
      - shape: (nlat, nlon)
    """
    with rasterio.open(tif_path) as src:
        transform = src.transform
        width = src.width
        height = src.height
        crs = src.crs
        bounds = src.bounds
        nodata = src.nodata

    rows = np.arange(height)
    cols = np.arange(width)
    col2d, row2d = np.meshgrid(cols, rows)

    xs, ys = xy(transform, row2d, col2d, offset="center")
    x2d = np.array(xs, dtype=np.float64)
    y2d = np.array(ys, dtype=np.float64)

    if crs is None:
        raise ValueError("TARGET_TIF has no CRS. Please define CRS for TARGET_TIF.")

    if crs.to_epsg() == 4326:
        lon2d, lat2d = x2d, y2d
        bounds_wgs84 = (bounds.left, bounds.bottom, bounds.right, bounds.top)
    else:
        transformer = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
        lon2d, lat2d = transformer.transform(x2d, y2d)

        bx = np.array([bounds.left, bounds.right, bounds.right, bounds.left], dtype=np.float64)
        by = np.array([bounds.bottom, bounds.bottom, bounds.top, bounds.top], dtype=np.float64)
        blon, blat = transformer.transform(bx, by)
        bounds_wgs84 = (
            float(np.min(blon)), float(np.min(blat)),
            float(np.max(blon)), float(np.max(blat))
        )

    if lon2d.ndim == 1:
        lon2d = lon2d.reshape((height, width))
        lat2d = lat2d.reshape((height, width))

    tgt_lons = lon2d[0, :]
    tgt_lats = lat2d[:, 0]

    return tgt_lons, tgt_lats, lon2d, lat2d, bounds_wgs84, (height, width), nodata


def subset_source(ds, varname, bounds_wgs84, buffer_deg):
    """Subset ERA5 data to bounds + buffer."""
    min_lon, min_lat, max_lon, max_lat = bounds_wgs84
    lon_min = min_lon - buffer_deg
    lon_max = max_lon + buffer_deg
    lat_min = min_lat - buffer_deg
    lat_max = max_lat + buffer_deg

    lat_vals = ds["latitude"].values
    if lat_vals[0] > lat_vals[-1]:
        lat_slice = slice(lat_max, lat_min)  # descending
    else:
        lat_slice = slice(lat_min, lat_max)  # ascending

    sub = ds[[varname]].sel(longitude=slice(lon_min, lon_max), latitude=lat_slice)
    return sub


def nearest_interp_to_target(sub, varname, tgt_lons, tgt_lats):
    """
    Nearest-neighbor interpolation to the target lon/lat grid using xarray.
    This matches the Spangler paper's stated interpolation method for ERA5 FDIR.
    """
    da = sub[varname]

    # Ensure source coords are monotonic for xarray.interp
    if da["latitude"].values[0] > da["latitude"].values[-1]:
        da = da.sortby("latitude")
    if da["longitude"].values[0] > da["longitude"].values[-1]:
        da = da.sortby("longitude")

    # Keep target latitude ordering consistent with xarray.interp, then flip back if needed
    tgt_lats_interp = np.array(tgt_lats)
    flip_lat_back = False
    if tgt_lats_interp[0] > tgt_lats_interp[-1]:
        tgt_lats_interp = tgt_lats_interp[::-1]
        flip_lat_back = True

    out = da.interp(
        longitude=xr.DataArray(tgt_lons, dims="longitude"),
        latitude=xr.DataArray(tgt_lats_interp, dims="latitude"),
        method="nearest"
    )

    if flip_lat_back:
        out = out.sel(latitude=out.latitude[::-1])

    return out


def export_timeseries_to_nc(out_path, data_3d, times, tgt_lons, tgt_lats, varname, units, fill_value=-9999.0):
    """
    data_3d shape: (time, lat, lon)
    """
    da = xr.DataArray(
        data_3d.astype(np.float32),
        dims=("time", "latitude", "longitude"),
        coords={
            "time": times,
            "latitude": tgt_lats,
            "longitude": tgt_lons,
        },
        name=varname,
        attrs={
            "long_name": f"ERA5 {varname} (nearest-neighbor interpolated to target grid)",
            "units": units,
        },
    )

    ds_out = xr.Dataset(
        {varname: da},
        attrs={
            "source": "ERA5",
            "target_grid": os.path.basename(TARGET_TIF),
            "interpolation": "Nearest-neighbor in lon/lat space (Spangler-style)",
            "crs": "EPSG:4326",
        },
    )

    encoding = {
        varname: {
            "zlib": True,
            "complevel": 4,
            "dtype": "float32",
            "_FillValue": np.float32(fill_value),
        }
    }

    ds_out.to_netcdf(out_path, encoding=encoding)


def process_one_year(nc_path, tgt_lons, tgt_lats, bounds_wgs84):
    print(f"\nOpening: {nc_path}")
    ds = xr.open_dataset(nc_path)

    for varname in ERA5_VAR:
        if varname not in ds.data_vars:
            continue

        units = ds[varname].attrs.get("units", "")
        sub = subset_source(ds, varname, bounds_wgs84, BUFFER_DEG)

        time_dim = None
        for cand in ("valid_time", "time"):
            if cand in sub.coords:
                time_dim = cand
                break
        if time_dim is None:
            raise KeyError("No time coordinate found.")

        times = sub[time_dim].values
        year_str = str(times[0])[:4]

        year_out_dir = os.path.join(OUT_DIR, year_str)
        os.makedirs(year_out_dir, exist_ok=True)

        print(f"Nearest-neighbor interpolating {len(times)} timesteps for {varname} ({year_str})...")
        out_da = nearest_interp_to_target(sub, varname, tgt_lons, tgt_lats)
        out_da = out_da.transpose(time_dim, "latitude", "longitude")
        out_3d = out_da.values

        out_path = os.path.join(year_out_dir, f"era5_{varname}_nn_{year_str}.nc")
        export_timeseries_to_nc(out_path, out_3d, times, tgt_lons, tgt_lats, varname, units)

        print(f"Finished {varname} {year_str}")

    ds.close()


#=========== main ===============
os.makedirs(OUT_DIR, exist_ok=True)

print(f"Reading target grid from GeoTIFF:\n  {TARGET_TIF}")
tgt_lons, tgt_lats, lon2d, lat2d, bounds_wgs84, shape, nodata = read_target_grid_from_tif(TARGET_TIF)

print("Target grid info (from TARGET_TIF):")
print(f"  shape (rows, cols): {shape}")
print(f"  bounds WGS84: minlon={bounds_wgs84[0]:.6f}, minlat={bounds_wgs84[1]:.6f}, "
      f"maxlon={bounds_wgs84[2]:.6f}, maxlat={bounds_wgs84[3]:.6f}")
print(f"  approx lon step: {float(np.nanmedian(np.diff(tgt_lons))):.6f} deg")
print(f"  approx lat step: {float(np.nanmedian(np.diff(tgt_lats))):.6f} deg")

nc_files = sorted(glob.glob(os.path.join(IN_DIR, "*.nc")))
if not nc_files:
    raise FileNotFoundError(f"No .nc files found in {IN_DIR}")

for nc_path in nc_files:
    process_one_year(nc_path, tgt_lons, tgt_lats, bounds_wgs84)

print("\nAll years processed.")

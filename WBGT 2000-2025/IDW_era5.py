"""
ERA5 monthly files -> IDW to target grid defined by GeoTIFF -> export one NetCDF per month.

Reads target grid geoinfo (bounds, resolution, CRS, shape) from: 
GeoTIFF file derived from ERA5-Land

Outputs a folder containing individual .nc file for each month.

Dependencies:
  pip install rasterio pyproj xarray netcdf4 numpy scipy
"""

import os
import glob
import numpy as np
import xarray as xr
import rasterio
from rasterio.transform import xy
from scipy.spatial import cKDTree
from pyproj import Transformer
from pathlib import Path

# =========================
# USER SETTINGS (EDIT THESE)
# =========================

IN_DIR = Path("era5_data")
OUT_DIR = Path("era5_interpolated")

TARGET_TIF = Path("target_tif.tif")

# IDW settings
K_NEIGHBORS = 8
POWER = 2.0

# Subset buffer in degrees (applied in WGS84 lon/lat after converting bounds)
BUFFER_DEG = 0.20

# Variable name in ERA5 to interpolate
ERA5_VAR = ["fdir"]

# =========================
# END USER SETTINGS
# =========================



def read_target_grid_from_tif(tif_path):
    """
    Returns target grid definition derived from GeoTIFF:
      - tgt_lon2d, tgt_lat2d: 2D arrays of cell-center lon/lat in EPSG:4326
      - tgt_lons, tgt_lats: 1D vectors if grid is regular in lon/lat (we still return them for convenience)
      - target_points: Nx2 lon/lat pairs (for KDTree query)
      - bounds_wgs84: (minlon, minlat, maxlon, maxlat)
      - shape: (nlat, nlon)
    """
    with rasterio.open(tif_path) as src:
        transform = src.transform
        width = src.width
        height = src.height
        crs = src.crs
        bounds = src.bounds  # in src CRS
        nodata = src.nodata

    # Build arrays of row/col indices
    rows = np.arange(height)
    cols = np.arange(width)

    # Get X/Y for cell centers:
    # rasterio.transform.xy can generate coordinates; use meshgrid for full grid.
    col2d, row2d = np.meshgrid(cols, rows)

    xs, ys = xy(transform, row2d, col2d, offset="center")
    x2d = np.array(xs, dtype=np.float64)
    y2d = np.array(ys, dtype=np.float64)

    # Convert target grid coordinates to lon/lat if needed
    if crs is None:
        raise ValueError("TARGET_TIF has no CRS. Please define CRS for Standard.tif.")

    if crs.to_epsg() == 4326:
        lon2d, lat2d = x2d, y2d
        bounds_wgs84 = (bounds.left, bounds.bottom, bounds.right, bounds.top)
    else:
        transformer = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
        lon2d, lat2d = transformer.transform(x2d, y2d)

        # Convert bounds too (do 4 corners and take min/max)
        bx = np.array([bounds.left, bounds.right, bounds.right, bounds.left], dtype=np.float64)
        by = np.array([bounds.bottom, bounds.bottom, bounds.top, bounds.top], dtype=np.float64)
        blon, blat = transformer.transform(bx, by)
        bounds_wgs84 = (float(np.min(blon)), float(np.min(blat)), float(np.max(blon)), float(np.max(blat)))

    # Derive 1D lon/lat vectors if the grid is rectilinear in lon/lat
    # (Even if it's not perfectly rectilinear, we can still use lon2d/lat2d for points.)
    # Many WGS84 rasters are regular, so this works:
    if lon2d.ndim == 1:
        lon2d = lon2d.reshape((height, width))
        lat2d = lat2d.reshape((height, width))

    tgt_lons = lon2d[0, :]
    tgt_lats = lat2d[:, 0]

    target_points = np.column_stack([lon2d.ravel(), lat2d.ravel()])

    return tgt_lons, tgt_lats, lon2d, lat2d, target_points, bounds_wgs84, (height, width), nodata


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


def build_kdtree_and_weights(src_lons, src_lats, target_points, k, power, eps=1e-12):
    """KDTree in lon/lat degree space; compute neighbor indices + IDW weights."""
    lon2d, lat2d = np.meshgrid(src_lons, src_lats)
    src_points = np.column_stack([lon2d.ravel(), lat2d.ravel()])

    tree = cKDTree(src_points)
    dists, idx = tree.query(target_points, k=k, workers=-1)

    if k == 1:
        dists = dists[:, None]
        idx = idx[:, None]

    dists = np.maximum(dists, eps)
    weights = 1.0 / (dists ** power)
    weights /= weights.sum(axis=1, keepdims=True)

    # Exact match handling
    zero_mask = dists <= (eps * 10.0)
    if np.any(zero_mask):
        rows = np.where(zero_mask.any(axis=1))[0]
        for r in rows:
            j = np.where(zero_mask[r])[0][0]
            weights[r, :] = 0.0
            weights[r, j] = 1.0

    return idx, weights


def idw_interpolate_day(src_day_2d, idx, weights, out_shape):
    """IDW interpolate one day 2D (lat, lon) to target points."""
    src_flat = src_day_2d.ravel()
    neigh_vals = src_flat[idx]

    nan_mask = np.isnan(neigh_vals)
    if np.any(nan_mask):
        w = weights.copy()
        w[nan_mask] = 0.0
        wsum = w.sum(axis=1, keepdims=True)
        all_nan = (wsum[:, 0] == 0.0)
        wsum[all_nan, :] = 1.0
        w /= wsum
        out = np.sum(neigh_vals * w, axis=1)
        out[all_nan] = np.nan
    else:
        out = np.sum(neigh_vals * weights, axis=1)

    return out.reshape(out_shape)


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
            "long_name": f"ERA5 {varname} (IDW interpolated to target grid)",
            "units": units,
        },
    )

    ds_out = xr.Dataset(
        {varname: da},
        attrs={
            "source": "ERA5",
            "target_grid": os.path.basename(TARGET_TIF),
            "interpolation": f"IDW (k={K_NEIGHBORS}, power={POWER}) in lon/lat degree space",
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


def process_one_file(nc_path, tgt_lons, tgt_lats, target_points, bounds_wgs84, idx, weights):
    
    print(f"\nOpening: {nc_path}")
    ds = xr.open_dataset(nc_path)
    
    filename = os.path.basename(nc_path)
    parts = filename.replace(".nc", "").split("_")
    year = parts[1]
    month = parts[2]
    
    for varname in ERA5_VAR:
        if varname not in ds.data_vars:
            continue

        units = ds[varname].attrs.get("units", "")
        sub = subset_source(ds, varname, bounds_wgs84, BUFFER_DEG)

        src_lons = sub["longitude"].values
        src_lats = sub["latitude"].values
        
        # time coordinate
        time_dim = None
        for cand in ("valid_time", "time"):
            if cand in sub.coords:
                time_dim = cand
                break
        if time_dim is None:
            raise KeyError("No time coordinate found.")
    
        times = sub[time_dim].values

        ntime = len(times)
        out_shape = (ntime, len(tgt_lats), len(tgt_lons))
        out_3d = np.full(out_shape, np.nan, dtype=np.float32)
    
        print(f"Interpolating {ntime} timesteps for {varname} ({year}-{month})...")
    
        for i, t in enumerate(times):
            src_2d = sub[varname].sel({time_dim: t}).values
            out_2d = idw_interpolate_day(src_2d, idx, weights, out_shape[1:])
            out_3d[i, :, :] = out_2d
    
            if (i + 1) % 24 == 0 or (i + 1) == ntime:
                print(f"  done {i+1}/{ntime}")
                
        out_path = os.path.join(OUT_DIR, f"era5_{year}_{month}_idw.nc")
    
        export_timeseries_to_nc(out_path, out_3d, times, tgt_lons, tgt_lats, varname, units,)

    ds.close()
    print(f"Finished {year}-{month}")




#===========main===============
os.makedirs(OUT_DIR, exist_ok=True)

print(f"Reading target grid from GeoTIFF:\n  {TARGET_TIF}")
tgt_lons, tgt_lats, lon2d, lat2d, target_points, bounds_wgs84, shape, nodata = read_target_grid_from_tif(TARGET_TIF)

print("Target grid info (from Standard.tif):")
print(f"  shape (rows, cols): {shape}")
print(f"  bounds WGS84: minlon={bounds_wgs84[0]:.6f}, minlat={bounds_wgs84[1]:.6f}, "
      f"maxlon={bounds_wgs84[2]:.6f}, maxlat={bounds_wgs84[3]:.6f}")
print(f"  approx lon step: {float(np.nanmedian(np.diff(tgt_lons))):.6f} deg")
print(f"  approx lat step: {float(np.nanmedian(np.diff(tgt_lats))):.6f} deg")

nc_files = sorted(glob.glob(os.path.join(IN_DIR, "*.nc")))
if not nc_files:
    raise FileNotFoundError(f"No .nc files found in {IN_DIR}")

# --------------------------------------------------
# Build KDTree + weights ONCE using first file
# --------------------------------------------------
print("\nBuilding KDTree + weights (once)...")

sample_ds = xr.open_dataset(nc_files[0])

# assume first variable in ERA5_VAR exists
sample_var = ERA5_VAR[0]
sub_sample = subset_source(sample_ds, sample_var, bounds_wgs84, BUFFER_DEG)

src_lons = sub_sample["longitude"].values
src_lats = sub_sample["latitude"].values

idx, weights = build_kdtree_and_weights(
    src_lons=src_lons,
    src_lats=src_lats,
    target_points=target_points,
    k=K_NEIGHBORS,
    power=POWER
)

sample_ds.close()

print("KDTree built successfully.")

# --------------------------------------------------
# Process each file
# --------------------------------------------------
for nc_path in nc_files:
    filename = os.path.basename(nc_path)
    parts = filename.replace(".nc", "").split("_")
    year = parts[1]
    month = parts[2]

    out_path = os.path.join(OUT_DIR, f"era5_{year}_{month}_idw.nc")

    if os.path.exists(out_path):
        print(f"Skipping existing: {out_path}")
        continue

    process_one_file(nc_path, tgt_lons, tgt_lats, target_points, bounds_wgs84, idx, weights)

print("\nAll months processed.")



"""
Regrid WBGT monthly data to a master GeoTIFF grid and write one NetCDF per month.

Inputs
------
- One WGS84 master grid GeoTIFF (e.g., 0.01 x 0.01 degree).
- WBGT monthly files in one folder, named like:
    wbgt_2000_01.nc
    wbgt_2000_02.nc
    ...

Outputs
-------
- One file per month, containing: wbgt
- The output grid exactly matches the master GeoTIFF.
- Optional 2-D lat/lon coordinate variables controlled by EXPORT_LAT_LON.

Required packages:
    pip install numpy xarray netCDF4 rasterio pyproj
"""

from __future__ import annotations
from typing import Tuple

from pathlib import Path
import numpy as np
import rasterio
from pyproj import CRS, Transformer
from rasterio.crs import CRS as RioCRS
from rasterio.enums import Resampling
from rasterio.transform import Affine
from rasterio.warp import reproject
import xarray as xr


# =============================================================================
# USER SETTINGS
# =============================================================================

WBGT_FOLDER = Path(r"/mnt/HDD7/juliae/WBGT-heat-project/WBGT 2000-2025/wbgt/")
MASTER_GRID_TIF = Path(r"/mnt/HDD7/juliae/WBGT-heat-project/WBGT 2000-2025/ERA5_Land_grid0pnt01_CAfull.tif")
OUTPUT_FOLDER = Path(r"/mnt/HDD7/juliae/WBGT-heat-project/WBGT 2000-2025/wbgt_regridded")

START_YEAR = 2000
END_YEAR = 2025

# Add 2-D variables named "lat" and "lon" to every NetCDF.
EXPORT_LAT_LON = False

# Apply the valid-data mask of the master GeoTIFF.
# Pixels invalid/NoData in the master grid become NoData in all output variables.
APPLY_MASTER_MASK = True

# use nearest resampling
# Other choices supported here: "nearest", "cubic".
RESAMPLING_METHOD = "nearest"

# Store monthly files in OUTPUT_FOLDER
YEAR_SUBFOLDERS = False

# Skip existing monthly files unless OVERWRITE is True.
OVERWRITE = False

# NetCDF compression level (1-9).
COMPRESSION_LEVEL = 4

# Output NoData / _FillValue
NODATA = np.float32(-9999.0)


# =============================================================================
# HELPERS
# =============================================================================

def grid_center_coordinates(transform: Affine, width: int, height: int):
    """Return 1-D x/y pixel-center coordinates for a north-up master grid."""
    if not np.isclose(transform.b, 0.0) or not np.isclose(transform.d, 0.0):
        raise ValueError("Rotated master grids are not supported by this script.")

    x = transform.c + (np.arange(width, dtype=np.float64) + 0.5) * transform.a
    y = transform.f + (np.arange(height, dtype=np.float64) + 0.5) * transform.e
    return x, y

def get_source_crs(ds: xr.Dataset) -> RioCRS:
    """
    Determine the CRS of the WBGT NetCDF.

    This expects either:
      1. a CF grid-mapping variable, or
      2. a CRS/WKT attribute.

    For a WGS84 WBGT file, EPSG:4326 is used when no CRS variable
    is present but the dataset contains longitude/latitude coordinates.
    """

    # -------------------------------------------------------------------------
    # Look for a CF grid mapping variable.
    # -------------------------------------------------------------------------
    for name in ds.variables:
        var = ds[name]

        grid_mapping_name = var.attrs.get("grid_mapping_name")

        if grid_mapping_name:
            try:
                crs = CRS.from_cf(dict(var.attrs))
                return RioCRS.from_wkt(crs.to_wkt())
            except Exception:
                pass

    # -------------------------------------------------------------------------
    # Look for common CRS attributes.
    # -------------------------------------------------------------------------
    for name in ("crs", "spatial_ref", "crs_wkt"):
        if name in ds.attrs:
            try:
                return RioCRS.from_user_input(ds.attrs[name])
            except Exception:
                pass

    # -------------------------------------------------------------------------
    # If the dataset has lat/lon coordinates, assume WGS84.
    # -------------------------------------------------------------------------
    coord_names = set(ds.coords) | set(ds.variables)

    has_lat = "lat" in coord_names or "latitude" in coord_names
    has_lon = "lon" in coord_names or "longitude" in coord_names

    if has_lat and has_lon:
        return RioCRS.from_epsg(4326)

    raise ValueError(
        "Could not determine the CRS of the WBGT NetCDF. "
        "The file needs CF CRS metadata or latitude/longitude coordinates."
    )

def get_xy_coordinates(
    ds: xr.Dataset,
):
    """
    Find the horizontal coordinates in the WBGT NetCDF.
    Supports either: x/y
    or: lon/lat
        longitude/latitude
    """
    if "x" in ds.coords and "y" in ds.coords:
        return "x", "y", ds["x"].values, ds["y"].values
    if "lon" in ds.coords and "lat" in ds.coords:
        return "lon", "lat", ds["lon"].values, ds["lat"].values
    if "longitude" in ds.coords and "latitude" in ds.coords:
        return ("longitude",
                "latitude",
                ds["longitude"].values,
                ds["latitude"].values,)
    raise ValueError("Could not find horizontal coordinates. "
                     "Expected x/y, lon/lat, or longitude/latitude.")

def prepare_source_geometry(ds: xr.Dataset,):
    """
    Build a raster transform for the source WBGT grid.
    Returns: src_transform, src_crs, flip_y, flip_x
    """
    x_name, y_name, x, y = get_xy_coordinates(ds)
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)

    if x.ndim != 1 or y.ndim != 1:
        raise ValueError("This script expects 1-D horizontal coordinates.")
    if len(x) < 2 or len(y) < 2:
        raise ValueError("Source x/y coordinates are too short to determine resolution.")

    dx = float(np.nanmedian(np.abs(np.diff(x))))
    dy = float(np.nanmedian(np.abs(np.diff(y))))

    if dx <= 0 or dy <= 0:
        raise ValueError("Invalid source coordinate spacing.")

    # Rasterio expects:
    #   columns = west -> east
    #   rows    = north -> south
    flip_x = x[1] < x[0]
    flip_y = y[1] > y[0]

    x_ordered = x[::-1] if flip_x else x
    y_ordered = y[::-1] if flip_y else y

    west = float(x_ordered[0] - dx / 2.0)
    north = float(y_ordered[0] + dy / 2.0)

    src_transform = Affine(dx, 0.0, west, 0.0, -dy, north,)

    src_crs = get_source_crs(ds)

    print(f"  Source coordinates: {x_name}, {y_name}")
    print(f"  Source CRS: {src_crs}")
    print(f"  Source resolution: {dx}, {dy}")
    print(f"  Source transform: {src_transform}")

    return src_transform, src_crs, flip_y, flip_x

def read_wbgt(ds: xr.Dataset, flip_y: bool, flip_x: bool,) -> np.ndarray:
    """
    Read WBGT as a 3-D float32 array with dimensions: (valid_time, y, x)
    Supports:
        (valid_time, y, x)
        (valid_time, lat, lon)
        (valid_time, latitude, longitude)
    Also allows a single 2-D field, which is treated as one time step.
    """
    if "wbgt" not in ds: raise KeyError("The NetCDF does not contain a variable named 'wbgt'.")
    da = ds["wbgt"]

    # -------------------------------------------------------------------------
    # Identify spatial dimensions.
    # -------------------------------------------------------------------------
    horizontal_pairs = [("y", "x"), ("lat", "lon"), ("latitude", "longitude"),]
    selected_pair = None

    for y_name, x_name in horizontal_pairs:
        if y_name in da.dims and x_name in da.dims:
            selected_pair = (y_name, x_name)
            break

    if selected_pair is None:
        raise ValueError(f"Could not identify spatial dimensions in wbgt: {da.dims}")

    y_name, x_name = selected_pair

    # -------------------------------------------------------------------------
    # Identify the time dimension.
    # -------------------------------------------------------------------------
    time_names = ["valid_time", "time", "datetime",]
    time_name = None
    for name in time_names:
        if name in da.dims:
            time_name = name
            break

    # -------------------------------------------------------------------------
    # Put dimensions into [time, y, x] order.
    # -------------------------------------------------------------------------
    if time_name is not None:
        da = da.transpose(time_name, y_name, x_name,)
        arr = da.values.astype(np.float32, copy=False,)
    else:
        # ---------------------------------------------------------------------
        # If there is no time dimension, treat the field as one timestep.
        # ---------------------------------------------------------------------
        da = da.transpose(y_name, x_name,)
        arr2d = da.values.astype(np.float32, copy=False,)
        arr = arr2d[np.newaxis, :, :]

    # -------------------------------------------------------------------------
    # Convert NaN to our internal NoData value.
    # -------------------------------------------------------------------------
    arr = np.where(np.isfinite(arr), arr, NODATA,).astype(np.float32, copy=False,)

    # -------------------------------------------------------------------------
    # Ensure raster orientation is north -> south and west -> east.
    # -------------------------------------------------------------------------
    if flip_y:
        arr = arr[:, ::-1, :]
    if flip_x:
        arr = arr[:, :, ::-1]
    return arr


def regrid_wbgt(src_array: np.ndarray, src_transform: Affine, src_crs: RioCRS,
    dst_shape, dst_transform: Affine, dst_crs: RioCRS, resampling: Resampling,) -> np.ndarray:
    """
    Regrid a 3-D WBGT array: (time, source_y, source_x)
    to: (time, destination_y, destination_x)
    """
    if src_array.ndim != 3:
        raise ValueError(f"Expected 3-D source array (time, y, x), "
                         f"got shape {src_array.shape}")
    ntime = src_array.shape[0]
    dst = np.full((ntime, dst_shape[0], dst_shape[1]), np.nan, dtype=np.float32,)

    for t in range(ntime):
        src = src_array[t]
        dst_t = np.full(dst_shape, NODATA, dtype=np.float32,)

        reproject(
            source=src,
            destination=dst_t,
            src_transform=src_transform,
            src_crs=src_crs,
            src_nodata=float(NODATA),
            dst_transform=dst_transform,
            dst_crs=dst_crs,
            dst_nodata=float(NODATA),
            resampling=resampling,
            num_threads=2,
            init_dest_nodata=True,
        )

        # Convert NoData to NaN.
        dst_t[dst_t == NODATA] = np.nan
        dst[t] = dst_t

    return dst



def make_lat_lon(x: np.ndarray, y: np.ndarray, dst_crs: RioCRS,) -> Tuple[np.ndarray, np.ndarray]:
    """Create 2-D latitude/longitude arrays at target pixel centers."""
    from pyproj import Transformer
    xx, yy = np.meshgrid(x, y)

    pyproj_dst = CRS.from_user_input(dst_crs)
    wgs84 = CRS.from_epsg(4326)

    if pyproj_dst == wgs84:
        lon = xx
        lat = yy
    else:
        transformer = Transformer.from_crs(pyproj_dst, wgs84, always_xy=True)
        lon, lat = transformer.transform(xx, yy)

    return lat.astype(np.float32), lon.astype(np.float32)


def build_output_dataset(wbgt: np.ndarray, valid_time, year: int, month: int,
                         x: np.ndarray, y: np.ndarray, dst_transform: Affine,
                         dst_crs: RioCRS, lat2d=None, lon2d=None,):
    """
    Create output Dataset with dimensions: valid_time, y, x
    """
    ds_out = xr.Dataset(
        data_vars={
            "wbgt": (
                ("valid_time", "y", "x"),
                wbgt,
            ),
        },
        coords={
            "valid_time": valid_time,
            "x": ("x", x),
            "y": ("y", y),
        },
        attrs={
            "title": "Hourly WBGT regridded to master GeoTIFF grid",
            "source": "Monthly WBGT NetCDF",
            "year": year,
            "month": month,
            "regridding": (
                f"Reprojected/resampled to master GeoTIFF "
                f"using {RESAMPLING_METHOD}"
            ),
        },
    )

    ds_out["wbgt"].attrs.update(
        long_name="Wet-bulb globe temperature",
    )

    # -------------------------------------------------------------------------
    # CRS
    # -------------------------------------------------------------------------
    pyproj_dst = CRS.from_user_input(dst_crs)

    crs_attrs = pyproj_dst.to_cf()

    crs_attrs.update(
        spatial_ref=pyproj_dst.to_wkt(),
        crs_wkt=pyproj_dst.to_wkt(),
        GeoTransform=" ".join(
            map(str, dst_transform.to_gdal())
        ),
    )

    ds_out["crs"] = xr.DataArray(
        np.int32(0),
        attrs=crs_attrs,
    )

    ds_out["wbgt"].attrs["grid_mapping"] = "crs"

    # -------------------------------------------------------------------------
    # Coordinate metadata
    # -------------------------------------------------------------------------
    if pyproj_dst.is_geographic:

        ds_out["x"].attrs.update(
            standard_name="longitude",
            long_name="longitude of pixel center",
            units="degrees_east",
            axis="X",
        )

        ds_out["y"].attrs.update(
            standard_name="latitude",
            long_name="latitude of pixel center",
            units="degrees_north",
            axis="Y",
        )

    else:

        ds_out["x"].attrs.update(
            standard_name="projection_x_coordinate",
            long_name="x coordinate of projection",
            units="m",
            axis="X",
        )

        ds_out["y"].attrs.update(
            standard_name="projection_y_coordinate",
            long_name="y coordinate of projection",
            units="m",
            axis="Y",
        )

    # -------------------------------------------------------------------------
    # Optional 2-D latitude/longitude
    # -------------------------------------------------------------------------
    if EXPORT_LAT_LON:

        ds_out = ds_out.assign_coords(
            lat=(("y", "x"), lat2d),
            lon=(("y", "x"), lon2d),
        )

        ds_out["lat"].attrs.update(
            standard_name="latitude",
            long_name="latitude",
            units="degrees_north",
        )

        ds_out["lon"].attrs.update(
            standard_name="longitude",
            long_name="longitude",
            units="degrees_east",
        )

    return ds_out


def write_netcdf(ds_out: xr.Dataset, output_path: Path,):
    """
    Write compressed NetCDF-4 output.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True,)

    encoding = {
        "wbgt": {
            "dtype": "float32",
            "zlib": True,
            "complevel": COMPRESSION_LEVEL,
            "_FillValue": float(NODATA),

            # Chunk by time so individual hours can be accessed efficiently.
            "chunksizes": (
                1,
                min(512, ds_out.sizes["y"]),
                min(512, ds_out.sizes["x"]),
            ),
        },
    }

    if EXPORT_LAT_LON:
        encoding["lat"] = {
            "dtype": "float32",
            "zlib": True,
            "complevel": COMPRESSION_LEVEL,
        }

        encoding["lon"] = {
            "dtype": "float32",
            "zlib": True,
            "complevel": COMPRESSION_LEVEL,
        }

    ds_out.to_netcdf(
        output_path,
        mode="w",
        format="NETCDF4",
        engine="netcdf4",
        encoding=encoding,
    )


# =============================================================================
# MAIN
# =============================================================================

def main():
    OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)

    resampling_lookup = {
        "nearest": Resampling.nearest,
        "bilinear": Resampling.bilinear,
        "cubic": Resampling.cubic,
    }
    method = RESAMPLING_METHOD.lower()
    if method not in resampling_lookup:
        raise ValueError(
            f"Unsupported RESAMPLING_METHOD={RESAMPLING_METHOD!r}. "
            f"Choose from {list(resampling_lookup)}."
        )
    resampling = resampling_lookup[method]

    # -------------------------------------------------------------------------
    # Read master grid once
    # -------------------------------------------------------------------------
    with rasterio.open(MASTER_GRID_TIF) as master:
        if master.crs is None:
            raise ValueError("Master GeoTIFF has no CRS.")

        dst_crs = master.crs
        dst_transform = master.transform
        dst_height = master.height
        dst_width = master.width
        dst_shape = (dst_height, dst_width)

        x_out, y_out = grid_center_coordinates(
            dst_transform, dst_width, dst_height
        )

        if APPLY_MASTER_MASK:
            master_valid_mask = master.dataset_mask() > 0
        else:
            master_valid_mask = np.ones(dst_shape, dtype=bool)

    print("Master grid:")
    print(f"  CRS: {dst_crs}")
    print(f"  Shape: {dst_shape}")
    print(f"  Resolution: ({dst_transform.a}, {abs(dst_transform.e)})")

    # Optional static 2-D lat/lon arrays, computed once.
    if EXPORT_LAT_LON:
        lat2d, lon2d = make_lat_lon(x_out, y_out, dst_crs)
    else:
        lat2d = lon2d = None

    # -------------------------------------------------------------------------
    # Process each month
    # -------------------------------------------------------------------------
    for year in range(START_YEAR, END_YEAR + 1):
        for month in range(1, 13):
            input_path = (WBGT_FOLDER / f"wbgt_{year}_{month:02d}.nc")
            # -------------------------------------------------------------
            # Respect YEAR_SUBFOLDERS
            # -------------------------------------------------------------
            if YEAR_SUBFOLDERS:
                output_path = (OUTPUT_FOLDER / str(year) / f"wbgt_{year}_{month:02d}_regridded.nc")
            else:
                output_path = (OUTPUT_FOLDER / f"wbgt_{year}_{month:02d}_regridded.nc")
            print(f"\n=== {year}-{month:02d} ===")
    
            # -------------------------------------------------------------
            # Check input.
            # -------------------------------------------------------------
            if not input_path.exists():
                print(f"  Missing input: {input_path}")    
                continue
    
            # -------------------------------------------------------------
            # Check existing output.
            # -------------------------------------------------------------
            if output_path.exists() and not OVERWRITE:
                print(f"  Exists, skip: {output_path}")
                continue
    
            # -------------------------------------------------------------
            # Open monthly WBGT file.
            # -------------------------------------------------------------
            with xr.open_dataset(input_path, mask_and_scale=True,) as ds:
                if "wbgt" not in ds:
                    raise KeyError(f"{input_path} does not contain "
                                   "a variable named 'wbgt'.")
                print(f"  Input variables: "
                      f"{list(ds.data_vars)}")
    
                # ---------------------------------------------------------
                # Check time dimension.
                # ---------------------------------------------------------
                if "valid_time" in ds["wbgt"].dims:
                    valid_time = ds["valid_time"].values
                elif "time" in ds["wbgt"].dims:
                    valid_time = ds["time"].values
                elif "datetime" in ds["wbgt"].dims:
                    valid_time = ds["datetime"].values
                else:
                    valid_time = np.array([
                        np.datetime64(f"{year}-{month:02d}-15T00:00:00")
                    ])
                    
                print(f"  Number of time steps: "
                      f"{len(valid_time)}")
    
                # ---------------------------------------------------------
                # Determine source grid.
                # ---------------------------------------------------------
                (src_transform, src_crs, flip_y, flip_x,) = prepare_source_geometry(ds)
    
                # ---------------------------------------------------------
                # Read WBGT.
                # ---------------------------------------------------------
                wbgt_src = read_wbgt(ds, flip_y=flip_y, flip_x=flip_x,)
                print(f"  Source shape: "
                      f"{wbgt_src.shape}")
    
                # ---------------------------------------------------------
                # Regrid every time step.
                # ---------------------------------------------------------
                wbgt_dst = regrid_wbgt(
                    src_array=wbgt_src,
                    src_transform=src_transform,
                    src_crs=src_crs,
                    dst_shape=dst_shape,
                    dst_transform=dst_transform,
                    dst_crs=dst_crs,
                    resampling=resampling,
                )
    
                print(f"  Output shape before mask: "
                      f"{wbgt_dst.shape}")
    
                # ---------------------------------------------------------
                # Apply master GeoTIFF mask.
                # ---------------------------------------------------------
                if APPLY_MASTER_MASK:
                    wbgt_dst[:, ~master_valid_mask] = np.nan
    
                # ---------------------------------------------------------
                # Create output Dataset.
                # ---------------------------------------------------------
                ds_out = build_output_dataset(
                    wbgt=wbgt_dst,
                    valid_time=valid_time,
                    year=year,
                    month=month,
                    x=x_out,
                    y=y_out,
                    dst_transform=dst_transform,
                    dst_crs=dst_crs,
                    lat2d=lat2d,
                    lon2d=lon2d,
                )
    
                # ---------------------------------------------------------
                # Copy WBGT attributes from input.
                # ---------------------------------------------------------
                for key, value in ds["wbgt"].attrs.items():
                    if key not in ds_out["wbgt"].attrs:
                        ds_out["wbgt"].attrs[key] = value
    
                # ---------------------------------------------------------
                # Write.
                # ---------------------------------------------------------
                write_netcdf(ds_out, output_path,)
                ds_out.close()
                print(f"  Wrote: {output_path}")

    print("\nFinished.")



if __name__ == "__main__":
    main()

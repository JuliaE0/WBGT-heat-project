# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 15:25:18 2025

@author: dell


Batch clip GeoTIFFs by polygon mask (no ArcPy)

Dependencies:
    pip install geopandas rasterio shapely fiona

Tested for Spyder / standard Python interpreter.
"""

import os
import rasterio
from rasterio.mask import mask
import geopandas as gpd
from pathlib import Path


def clip_raster_by_shape(in_raster_path, out_raster_path, mask_geom, mask_crs, force_nodata=None):
    """
    Open a raster, reproject the mask if needed, clip, and write result.

    Parameters
    ----------
    in_raster_path : str
        Path to source .tif
    out_raster_path : str
        Where clipped .tif should be written
    mask_geom : shapely geometry (multi or single)
        Union of all polygons from the shapefile
    mask_crs : CRS object compatible with GeoPandas (e.g. clip_gdf.crs)
        CRS of the mask_geom
    force_nodata : number or None
        If not None, pixels outside mask will be set to this value and written
        as the output nodata.
    """
    with rasterio.open(in_raster_path) as src:
        raster_crs = src.crs

        if mask_crs is None:
            raise ValueError("Mask shapefile has no CRS defined.")

        # Reproject mask if needed
        if mask_crs != raster_crs:
            mask_gdf = gpd.GeoDataFrame({"geometry": [mask_geom]}, crs=mask_crs)
            mask_gdf_proj = mask_gdf.to_crs(raster_crs)
            clip_geom_proj = mask_gdf_proj.geometry.iloc[0]
        else:
            clip_geom_proj = mask_geom

        # Handle nodata
        nodata_val = force_nodata if force_nodata is not None else src.nodata

        # safely read colormap (palette)
        colormap = None
        if src.count == 1:
            try:
                colormap = src.colormap(1)
            except Exception:
                colormap = None

        # Clip raster
        clipped_data, clipped_transform = mask(
            dataset=src,
            shapes=[clip_geom_proj.__geo_interface__],
            crop=True,
            nodata=nodata_val
        )

        # Update metadata
        out_meta = src.meta.copy()
        out_meta.update({
            "height": clipped_data.shape[1],
            "width": clipped_data.shape[2],
            "transform": clipped_transform,
            "nodata": nodata_val,
            "dtype": src.dtypes[0]   # for palette rasters
        })

    # Write output
    with rasterio.open(out_raster_path, "w", **out_meta) as dest:
        dest.write(clipped_data)

        # write colormap back
        if colormap:
            dest.write_colormap(1, colormap)


def main():
    # --- USER PATHS (edit these for your machine) ---
    mypath = Path("NLCD")              # root of input folders
    outpath = Path("NLCD_clipped")     # root of output folders
    shpname = Path("NLCD_clip_shape/NLCD_clip_extent.shp")    # clip polygon shapefile

    # Optional: set a nodata value for output. Use None to keep source nodata.
    forced_nodata_value = None
    # Example if you WANT to force nodata:
    # forced_nodata_value = -9999.0

    # 1. Read shapefile once
    clip_gdf = gpd.read_file(shpname)

    # dissolve all polygons into one geometry (acts like ExtractByMask_sa)
    union_geom = clip_gdf.unary_union   # shapely geometry (multi or single)
    union_crs = clip_gdf.crs            # save crs for reprojection later

    # 2. Make sure the top-level output directory exists
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    # 3. Walk through all subfolders and process matching rasters
    for root, _, files in os.walk(mypath):
        for tiffile in files:
            # mimic your arcpy condition: filename ends with "test.tif"
            if tiffile.lower().endswith(".tif"):
                in_raster = os.path.join(root, tiffile)

                # mirror the input folder structure under outpath
                rel_path = os.path.relpath(root, mypath)
                out_subdir = os.path.join(outpath, rel_path)
                if not os.path.exists(out_subdir):
                    os.makedirs(out_subdir)

                out_raster = os.path.join(out_subdir, tiffile)

                print("Clipping:")
                print("  IN :", in_raster)
                print("  OUT:", out_raster)

                clip_raster_by_shape(
                    in_raster_path=in_raster,
                    out_raster_path=out_raster,
                    mask_geom=union_geom,
                    mask_crs=union_crs,
                    force_nodata=forced_nodata_value
                )

    print("All done.")


# standard Python entry point so Spyder (or CLI) can run it cleanly
if __name__ == "__main__":
    main()


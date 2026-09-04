# -*- coding: utf-8 -*-
"""
Calculate distance to coast
"""

import geopandas as gpd
from pathlib import Path

# -----------------------------
# 1. Input / output paths
# -----------------------------
grid_path = Path(r"C:\Users\new account\Desktop\WBGT-heat-project\WBGT 0.01\distance to coast 0.01\ca_grid\ERA5_Land_grid0pnt01_CAfull.shp")
coast_path = Path(r"C:\Users\new account\Desktop\WBGT-heat-project\WBGT 0.01\distance to coast 0.01\tl_2025_us_coastline\tl_2025_us_coastline.shp")

out_path = Path(r"C:\Users\new account\Desktop\WBGT-heat-project\WBGT 0.01\distance to coast 0.01\CA_coast_dist_0.01.csv")

# -----------------------------
# 2. Read data
# -----------------------------
grid = gpd.read_file(grid_path)
coast = gpd.read_file(coast_path)

# -----------------------------
# 3. check CRS
# -----------------------------
print("Grid CRS:", grid.crs)
print("Coastline CRS:", coast.crs)

# -----------------------------
# 5. Project data only for distance calculation
# -----------------------------
target_crs = "EPSG:3310"

grid_proj = grid.to_crs(target_crs)
coast_proj = coast.to_crs(target_crs)

# -----------------------------
# 6. Create centroids in projected CRS
# -----------------------------
centroids_proj = grid_proj.copy()
centroids_proj["geometry"] = grid_proj.geometry.centroid

# -----------------------------
# 7. Merge coastline geometry
# -----------------------------
try:
    coast_union = coast_proj.geometry.union_all()
except AttributeError:
    coast_union = coast_proj.geometry.unary_union

# -----------------------------
# 8. Calculate distance in meters
# -----------------------------
centroids_proj["dist_coast_m"] = centroids_proj.geometry.distance(coast_union)
centroids_proj["dist_coast_km"] = centroids_proj["dist_coast_m"] / 1000

# ============================================================
# 9. Calculate centroid longitude / latitude
# ============================================================

# Convert centroids back to EPSG:4269 so longitude/latitude
# are stored in the original NAD83 geographic coordinates.

centroids_wgs84 = gpd.GeoSeries(
    centroids_proj.geometry,
    crs=target_crs
).to_crs("EPSG:4269")

centroids_proj["centroid_lon"] = centroids_wgs84.x
centroids_proj["centroid_lat"] = centroids_wgs84.y

# -----------------------------
# 10. Save csv
# -----------------------------
# Remove geometry because CSV cannot store shapefile geometry
output = centroids_proj.drop(columns="geometry")
output.to_csv(out_path, index=False)
print(f"Saved shapefile: {out_path}")

print("Done!")
print("Final output CRS:", grid.crs)
print(f"Output saved to: {out_path}")
print(centroids_proj[["dist_coast_m", "dist_coast_km"]].describe())
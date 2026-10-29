# -*- coding: utf-8 -*-
"""
review_05_uhi.py  -  surface urban heat island (SUHI) sensitivity for the V1 review
  SUHI = pixel LST - mean LST of rural pixels within a radius (1 km pixels)
  - urban / rural masks from MCD12Q1 (original), ESA WorldCover 2021 and GHSL 2025
  - radius 10, 15 and 20 km
  - share of urban pixels with no rural pixel inside the radius
  - Dhaka and Chattogram SUHI at 1 km and after averaging to the 0.05 degree grid
LST: mean of the monthly means of 2025 (months with >= 50 % coverage), Terra, day and night
Needs: BD_AQ_2025/LST_*_monthly_2025.tif, BD_AQ_review/controls_1km.tif
Output: Review_V1/05_uhi/
"""
import os
import numpy as np
import pandas as pd
from scipy import ndimage
import review_common as rc

RADII_KM = [10, 15, 20]
URBAN_MIN = 0.3          # built-up fraction of a 1 km pixel to call it urban (WorldCover, GHSL)
RURAL_MAX = 0.05         # built-up fraction below which a pixel is rural
WATER_MAX = 0.5          # pixels with more water than this are excluded
COVER_MIN = 50
CITIES = {"Dhaka": (23.8103, 90.4125), "Chattogram": (22.3569, 91.7832)}
CITY_RADIUS_KM = 15
OUT = rc.ensure_out("05_uhi")

grid = rc.Grid()
ctrl, trc, _, bc, names = rc.read_all(os.path.join(rc.REV_DIR, "controls_1km.tif"))
C = {nm: ctrl[i] for i, nm in enumerate(names)}


def lst_mean(var):
    a, tr, _, b, _ = rc.read_all(os.path.join(rc.DATA_DIR, "{}_monthly_2025.tif".format(var)))
    land = grid.land_on(a[0].shape, tr) > 0
    covs = np.array([100.0 * np.isfinite(a[m][land]).sum() / land.sum() for m in range(12)])
    use = [m for m in range(12) if covs[m] >= COVER_MIN]
    st = a[use]
    cnt = np.isfinite(st).sum(0)
    v = np.where(cnt > 0, np.nansum(st, 0) / np.maximum(cnt, 1), np.nan)
    return v, tr, b, use


def disk(radius_km, dy_km, dx_km):
    ry, rx = int(np.ceil(radius_km / dy_km)), int(np.ceil(radius_km / dx_km))
    yy, xx = np.mgrid[-ry:ry + 1, -rx:rx + 1]
    return ((yy * dy_km) ** 2 + (xx * dx_km) ** 2 <= radius_km ** 2).astype(float)


rows = []
for var in ["LST_Day", "LST_Night"]:
    lst, tr, b, use = lst_mean(var)
    if lst.shape != C["built_wc"].shape:
        raise SystemExit("controls_1km.tif and {} are on different grids".format(var))
    lat_c = (b.top + b.bottom) / 2.0
    dy, dx = abs(tr.e) * 111.32, tr.a * 111.32 * np.cos(np.radians(lat_c))
    r_idx, c_idx = np.indices(lst.shape)
    lon = tr.c + (c_idx + 0.5) * tr.a
    lat = tr.f + (r_idx + 0.5) * tr.e
    valid = np.isfinite(lst)

    masks = {
        "MCD12Q1 (original)": (np.nan_to_num(C["urban_mcd"]) >= 0.5,
                               (np.nan_to_num(C["urban_mcd"]) < 0.5) & (np.nan_to_num(C["water_mcd"], nan=1) < 0.5)),
        "WorldCover 2021": (np.nan_to_num(C["built_wc"]) >= URBAN_MIN,
                            (np.nan_to_num(C["built_wc"], nan=1) <= RURAL_MAX) & (np.nan_to_num(C["water_wc"], nan=1) < WATER_MAX)),
        "GHSL 2025": (np.nan_to_num(C["built_ghsl"]) >= URBAN_MIN,
                      (np.nan_to_num(C["built_ghsl"], nan=1) <= RURAL_MAX) & (np.nan_to_num(C["water_wc"], nan=1) < WATER_MAX)),
    }
    water = np.nan_to_num(C["water_wc"], nan=1) >= WATER_MAX
    for mname, (urban, rural) in masks.items():
        urban = urban & valid & ~water
        rural = rural & valid
        for R in RADII_KM:
            k = disk(R, dy, dx)
            s = ndimage.convolve(np.where(rural, lst, 0.0), k, mode="constant")
            c = ndimage.convolve(rural.astype(float), k, mode="constant")
            ref = np.where(c > 0, s / np.maximum(c, 1), np.nan)
            uhi = np.where(valid & ~water, lst - ref, np.nan)
            row = {"LST": var, "months_used": len(use), "mask": mname, "radius_km": R,
                   "urban_pixels": int(urban.sum()),
                   "urban_no_rural_ref_%": 100.0 * (urban & (c == 0)).sum() / max(urban.sum(), 1),
                   "SUHI_urban_mean": np.nanmean(uhi[urban]), "SUHI_urban_p95": np.nanpercentile(uhi[urban], 95),
                   "SUHI_max_1km": np.nanmax(uhi)}
            g = grid.grid_band(uhi, tr)
            row["SUHI_max_0.05deg"] = np.nanmax(g)
            for city, (la0, lo0) in CITIES.items():
                d = 6371.0 * np.hypot(np.radians(lon - lo0) * np.cos(np.radians(la0)), np.radians(lat - la0))
                cu = urban & (d <= CITY_RADIUS_KM)
                row[city + "_urban_pixels"] = int(cu.sum())
                row[city + "_SUHI_urban_mean_1km"] = np.nanmean(uhi[cu]) if cu.any() else np.nan
                row[city + "_no_rural_ref_%"] = 100.0 * (cu & (c == 0)).sum() / max(cu.sum(), 1)
                dg = 6371.0 * np.hypot(np.radians(grid.lon - lo0) * np.cos(np.radians(la0)),
                                       np.radians(grid.lat - la0))
                row[city + "_SUHI_mean_0.05deg_within_15km"] = np.nanmean(g[dg <= CITY_RADIUS_KM])
            rows.append(row)
            print("{} {:20s} R={:2d} km: urban px {:6d}, no rural ref {:5.1f}%, urban SUHI {:.2f}, "
                  "Dhaka {:.2f} (0.05 deg {:.2f}), max 1 km {:.2f}".format(
                      var, mname, R, row["urban_pixels"], row["urban_no_rural_ref_%"], row["SUHI_urban_mean"],
                      row["Dhaka_SUHI_urban_mean_1km"], row["Dhaka_SUHI_mean_0.05deg_within_15km"],
                      row["SUHI_max_1km"]))

# the original annual UHI raster for comparison
a, tr, _, _, _ = rc.read_all(os.path.join(rc.DATA_DIR, "UHI_annual_2025.tif"))
u0 = a[0]
g0 = grid.grid_band(u0, tr)
print("\nOriginal UHI_annual_2025.tif: max at 1 km {:.2f}, max on 0.05 deg grid {:.2f}".format(
    np.nanmax(u0), np.nanmax(g0)))
t = pd.DataFrame(rows)
t.to_csv(os.path.join(OUT, "05_suhi_sensitivity.csv"), index=False)
print("\nDone. Outputs in", OUT)

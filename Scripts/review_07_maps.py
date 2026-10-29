# -*- coding: utf-8 -*-
"""
review_07_maps.py  -  map checks for the V1 review
  (a) the straight diagonal line in the annual AOD map: is it a MODIS tile edge?
      MAIAC (MCD19A2) is produced on the MODIS sinusoidal grid, so a processing
      edge appears as a line of constant sinusoidal x. The script finds the
      x position with the largest step in AOD and draws it with the tile edges.
  (b) NO2 and SO2 near large coal power plants in and around Bangladesh
      (Farakka and Sagardighi in West Bengal, Rampal, Payra, Barapukuria, ...)
Output: Review_V1/07_maps/
"""
import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import ndimage
from review_common import DATA_DIR, BOUNDARY, PLANTS, ensure_out, read_all

OUT = ensure_out("07_maps")
bnd = gpd.read_file(BOUNDARY)

R_SIN = 6371007.181                  # MODIS sinusoidal sphere radius (m)
TILE = 1111950.5197665554            # MODIS tile width (m)
X_MIN = -20015109.354                # left edge of tile h00


def sin_x(lon, lat):
    return R_SIN * np.cos(np.radians(lat)) * np.radians(lon)


def lon_of_x(x, lat):
    return np.degrees(x / (R_SIN * np.cos(np.radians(lat))))


# ======================= (a) AOD line =======================
a, tr, crs, b, _ = read_all(os.path.join(DATA_DIR, "AOD_annual_2025.tif"))
aod = a[0]
rows, cols = np.indices(aod.shape)
lon = tr.c + (cols + 0.5) * tr.a
lat = tr.f + (rows + 0.5) * tr.e
xs = sin_x(lon, lat)

# search the x position (sinusoidal metres) with the largest AOD step across it
ok = np.isfinite(aod)
cands = np.arange(np.nanmin(xs[ok]) + 20000, np.nanmax(xs[ok]) - 20000, 1000.0)
steps = []
for x0 in cands:
    left = ok & (xs > x0 - 6000) & (xs <= x0)
    right = ok & (xs > x0) & (xs < x0 + 6000)
    if left.sum() > 200 and right.sum() > 200:
        steps.append((x0, np.nanmean(aod[right]) - np.nanmean(aod[left]), left.sum() + right.sum()))
steps = pd.DataFrame(steps, columns=["x_sin_m", "step_aod", "n_pixels"])
best = steps.iloc[steps["step_aod"].abs().idxmax()]
x_line = best["x_sin_m"]
h = int(np.floor((x_line - X_MIN) / TILE))
inside_km = (x_line - (X_MIN + h * TILE)) / 1000.0
# typical step size elsewhere, for comparison
typical = steps["step_aod"].abs().median()
txt = ("Largest AOD step across a line of constant sinusoidal x: x = {:.0f} m, step = {:+.3f}\n"
       "(median absolute step at other positions = {:.3f}).\n"
       "This line lies {:.0f} km east of the western edge of MODIS tile column h{:02d}; "
       "the nearest tile edge (h{:02d}/h{:02d}) runs through north-west Bangladesh.".format(
           x_line, best["step_aod"], typical, inside_km, h, h - 1, h))
print(txt)
steps.to_csv(os.path.join(OUT, "07a_aod_step_by_sinusoidal_x.csv"), index=False)

f = np.where(ok, aod, np.nanmean(aod))
hp = f - ndimage.gaussian_filter(f, 15)
hp[~ok] = np.nan
ext = [b.left, b.right, b.bottom, b.top]
fig, axes = plt.subplots(1, 2, figsize=(11, 6.6))
for ax, img, ttl, kw in [
        (axes[0], aod, "Annual AOD 2025", dict(cmap="viridis", vmin=np.nanpercentile(aod, 2),
                                                vmax=np.nanpercentile(aod, 98))),
        (axes[1], hp, "AOD minus its 15-pixel smooth (edges)", dict(cmap="RdBu_r", vmin=-0.04, vmax=0.04))]:
    im = ax.imshow(img, extent=ext, **kw)
    bnd.boundary.plot(ax=ax, color="black", lw=0.5)
    la = np.linspace(20.5, 26.8, 50)
    ax.plot(lon_of_x(x_line, la), la, "m--", lw=1.2, label="line of largest step")
    for hh in [h, h + 1]:
        ax.plot(lon_of_x(X_MIN + hh * TILE, la), la, "k-", lw=1.0, label="MODIS tile edge" if hh == h else None)
    ax.set_xlim(b.left, b.right)
    ax.set_ylim(b.bottom, b.top)
    ax.set_title(ttl, fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.04)
axes[0].legend(fontsize=7, loc="lower left")
plt.tight_layout()
plt.savefig(os.path.join(OUT, "07a_aod_line.png"), dpi=250)
plt.close(fig)
with open(os.path.join(OUT, "07a_aod_line.txt"), "w") as fh:
    fh.write(txt + "\nx_line_m = {:.0f}\n".format(x_line))

# ======================= (b) power plants =======================
pl = pd.DataFrame(PLANTS, columns=["plant", "lat", "lon", "MW", "country"])
summary = []
fig, axes = plt.subplots(1, 2, figsize=(11, 6.6))
for ax, var in zip(axes, ["NO2", "SO2"]):
    a, tr, crs, b, _ = read_all(os.path.join(DATA_DIR, "{}_annual_2025.tif".format(var)))
    v = a[0] * 1e6                                   # umol/m2
    rows, cols = np.indices(v.shape)
    lo = tr.c + (cols + 0.5) * tr.a
    la = tr.f + (rows + 0.5) * tr.e
    nat_med = np.nanmedian(v)
    nat_p90 = np.nanpercentile(v, 90)
    for _, p in pl.iterrows():
        dkm = 6371.0 * np.hypot(np.radians(lo - p["lon"]) * np.cos(np.radians(p["lat"])),
                                np.radians(la - p["lat"]))
        for rad in [25, 50]:
            m = (dkm <= rad) & np.isfinite(v)
            summary.append({"variable": var, "plant": p["plant"], "country": p["country"],
                            "MW": p["MW"], "radius_km": rad, "n_pixels_in_BD": int(m.sum()),
                            "mean_umol_m2": np.nanmean(v[m]) if m.any() else np.nan,
                            "national_median": nat_med, "national_p90": nat_p90})
    im = ax.imshow(v, extent=[b.left, b.right, b.bottom, b.top], cmap="viridis",
                   vmin=np.nanpercentile(v, 2), vmax=np.nanpercentile(v, 98))
    bnd.boundary.plot(ax=ax, color="black", lw=0.5)
    ax.scatter(pl["lon"], pl["lat"], s=pl["MW"] / 25.0, facecolor="none", edgecolor="red", lw=1.2)
    for _, p in pl.iterrows():
        ax.annotate(p["plant"].split(" (")[0], (p["lon"], p["lat"]), fontsize=6, color="red",
                    xytext=(3, 3), textcoords="offset points")
    ax.set_xlim(b.left - 0.5, b.right + 0.1)
    ax.set_ylim(b.bottom - 0.1, b.top + 0.1)
    ax.set_title(u"Annual {} 2025 (µmol m⁻²) and coal power plants".format(var), fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.04)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "07b_no2_so2_power_plants.png"), dpi=250)
plt.close(fig)
s = pd.DataFrame(summary)
s.to_csv(os.path.join(OUT, "07b_plants_no2_so2.csv"), index=False)
print("\nNO2 / SO2 (umol/m2) within 25 and 50 km of each plant (pixels inside Bangladesh only)")
print(s.round(1).to_string(index=False))
print("\nDone. Outputs in", OUT)

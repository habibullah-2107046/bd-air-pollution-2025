import matplotlib
matplotlib.use("Agg")

import os
import glob
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
from matplotlib.ticker import FuncFormatter

# ---------------- settings ----------------
ROOT      = r"C:\Users\User\Desktop\AirPollutionResearchWithRakibVai"
BOUNDARY  = os.path.join(ROOT, "boundary", "BGD_adm0.shp")
ADM1_DIR  = os.path.join(ROOT, "bgd_adm_bbs_20201113")     # divisions (optional)
WORLD_URL = "https://naciscdn.org/naturalearth/50m/cultural/ne_50m_admin_0_countries.zip"
WORLD_LOCAL = os.path.join(ROOT, "boundary", "ne_50m_admin_0_countries.gpkg")
OUT_DIR   = os.path.join(ROOT, "Figures_v2", "StudyArea")
# ------------------------------------------

# label offset (dx, dy) in degrees for cities whose label needs moving
LABEL_SHIFT = {"Barishal": (-0.55, 0.05), "Rangpur": (-0.62, 0.05)}

# divisional headquarters (lon, lat)
CITIES = {
    "Dhaka":      (90.41, 23.81), "Chattogram": (91.83, 22.36),
    "Rajshahi":   (88.60, 24.37), "Khulna":     (89.56, 22.85),
    "Sylhet":     (91.87, 24.90), "Barishal":   (90.35, 22.70),
    "Rangpur":    (89.25, 25.75), "Mymensingh": (90.41, 24.75),
}
DIV_COLOURS = ["#e5f0d8", "#f6e8c3", "#dcebf5", "#f2dede",
               "#e8e0f0", "#fbe6d0", "#d9efe9", "#eeeeee"]


def dms(v, axis):
    hemi = ("N" if v >= 0 else "S") if axis == "y" else ("E" if v >= 0 else "W")
    d = int(abs(v)); m = int(round((abs(v) - d) * 60))
    if m == 60:
        d, m = d + 1, 0
    return u"{}\u00b0{}'0\"{}".format(d, m, hemi)


def north_arrow(ax, x=0.90, y=0.86, size=0.05):
    t = ax.transAxes; h = size * 0.30
    ax.add_patch(Polygon([[x, y + size], [x + h, y - size * 0.6], [x, y - size * 0.2]],
                         fc="black", ec="black", lw=0.5, transform=t, zorder=6))
    ax.add_patch(Polygon([[x, y + size], [x - h, y - size * 0.6], [x, y - size * 0.2]],
                         fc="white", ec="black", lw=0.5, transform=t, zorder=6))
    ax.text(x, y + size + 0.01, "N", transform=t, ha="center", fontsize=10, fontweight="bold")


def scale_bar(ax, total_km=100, x=0.60, y=0.05):
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    km_per_deg = 111.32 * np.cos(np.radians((y0 + y1) / 2))
    w = total_km / km_per_deg / (x1 - x0) / 4        # one segment, axes fraction
    t = ax.transAxes
    for k in range(4):
        ax.add_patch(Rectangle((x + k * w, y), w, 0.012, transform=t,
                               fc="black" if k % 2 == 0 else "white", ec="black", lw=0.5, zorder=6))
    for k in (0, 2, 4):
        ax.text(x + k * w, y + 0.02, "{:g}".format(total_km * k / 4), transform=t,
                ha="center", fontsize=7)
    ax.text(x + 4 * w + 0.015, y + 0.006, "km", transform=t, va="center", fontsize=7)


# ---- load data ----
bd = gpd.read_file(BOUNDARY).to_crs(4326)

adm1 = None
shp = glob.glob(os.path.join(ADM1_DIR, "**", "*adm1*.shp"), recursive=True)
if shp:
    adm1 = gpd.read_file(shp[0]).to_crs(4326)
    name_col = next((c for c in ["ADM1_EN", "NAME_1", "ADM1_NAME", "name"] if c in adm1.columns), None)
    print("divisions: {} ({} polygons)".format(os.path.basename(shp[0]), len(adm1)))
else:
    print("no division shapefile found - map will show the national boundary only")

if os.path.exists(WORLD_LOCAL):
    world = gpd.read_file(WORLD_LOCAL)
else:
    print("downloading world countries for the inset (one time only)...")
    world = gpd.read_file(WORLD_URL)
    world.to_file(WORLD_LOCAL, driver="GPKG")
for c in ["ADMIN", "NAME", "SOVEREIGNT"]:              # draw Bangladesh only from your own boundary
    if c in world.columns:
        world = world[world[c] != "Bangladesh"]
        break

# ---- main map ----
fig = plt.figure(figsize=(7.5, 8.5))
ax = fig.add_axes([0.10, 0.07, 0.85, 0.88])

world.plot(ax=ax, color="#f4f4f4", edgecolor="#9a9a9a", lw=0.5, zorder=1)
if adm1 is not None:
    adm1.plot(ax=ax, color=[DIV_COLOURS[i % len(DIV_COLOURS)] for i in range(len(adm1))],
              edgecolor="#555555", lw=0.6, zorder=2)
    if name_col:
        for _, r in adm1.iterrows():
            p = r.geometry.representative_point()
            ax.text(p.x, p.y - 0.12, r[name_col], ha="center", fontsize=7.5,
                    color="#444444", style="italic", zorder=5)
else:
    bd.plot(ax=ax, color="#e5f0d8", zorder=2)
bd.boundary.plot(ax=ax, color="black", lw=1.0, zorder=3)

for name, (lon, lat) in CITIES.items():
    ax.plot(lon, lat, marker="*" if name == "Dhaka" else "o",
            ms=11 if name == "Dhaka" else 5, color="#a32020", mec="black", mew=0.5, zorder=6)
    dx, dy = LABEL_SHIFT.get(name, (0.07, 0.05))
    ax.text(lon + dx, lat + dy, name, fontsize=8, fontweight="bold", zorder=6,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.6, pad=0.5))

ax.text(90.9, 21.0, "Bay of Bengal", fontsize=10, style="italic", color="#2c5f9e", ha="center")
ax.text(88.3, 23.0, "INDIA", fontsize=9, color="#777777")
ax.text(92.45, 21.6, "MYANMAR", fontsize=9, color="#777777", rotation=90)

minx, miny, maxx, maxy = bd.total_bounds
ax.set_xlim(minx - 0.4, maxx + 0.4)
ax.set_ylim(miny - 0.4, maxy + 0.3)
ax.set_facecolor("#dcebf5")                             # sea colour
ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: dms(v, "x")))
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: dms(v, "y")))
ax.tick_params(labelsize=7)
ax.tick_params(axis="y", labelrotation=90)
north_arrow(ax)
scale_bar(ax)

# ---- inset: location in South Asia ----
ins = ax.inset_axes([0.02, 0.02, 0.30, 0.20])   # inside the map frame, bottom-left   # bottom-left corner (sea / India)
world.plot(ax=ins, color="#e6e6e6", edgecolor="#8a8a8a", lw=0.4)
bd.plot(ax=ins, color="#a32020")
ins.add_patch(Rectangle((minx, miny), maxx - minx, maxy - miny, fill=False,
                        ec="black", lw=0.8))
ins.set_xlim(60, 105)
ins.set_ylim(5, 38)
ins.set_facecolor("#dcebf5")
ins.set_xticks([]); ins.set_yticks([])
ins.set_title("South Asia", fontsize=7, pad=2)

if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR)
out = os.path.join(OUT_DIR, "study_area_map.png")
plt.savefig(out, dpi=400)
plt.close(fig)
print("saved: " + out)
print("Done.")

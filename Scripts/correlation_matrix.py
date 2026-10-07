import matplotlib
matplotlib.use("Agg")

import os
import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
import matplotlib.pyplot as plt
from scipy import stats

# ---------------- settings ----------------
YEAR     = 2025
NODATA   = -9999
CELL_DEG = 0.05            # common grid ~5.5 km (close to TROPOMI pixel size)
MIN_LAND = 0.5             # a cell needs >= 50% land pixels to be used

ROOT     = r"C:\Users\User\Desktop\AirPollutionResearchWithRakibVai"
DATA_DIR = os.path.join(ROOT, "BD_AQ_2025")
OUT_DIR  = os.path.join(ROOT, "Figures_v2", "Statistics")
# ------------------------------------------

O3_TO_DU = 1.0 / 4.4615e-4

# file prefix -> (label for the plot, unit factor)
VARIABLES = {
    "LST_Day":   ("LST Day",       1.0),
    "LST_Night": ("LST Night",     1.0),
    "UHI":       ("UHI",           1.0),
    "AOD":       ("AOD",           1.0),
    "UVAI":      ("UVAI",          1.0),
    "CH4":       (u"CH\u2084",     1.0),
    "HCHO":      ("HCHO",          1.0),
    "O3":        (u"O\u2083",      O3_TO_DU),
    "NO2":       (u"NO\u2082",     1.0),
    "CO":        ("CO",            1.0),
    "SO2":       (u"SO\u2082",     1.0),
}
MODIS_1KM = ["LST_Day", "LST_Night", "UHI"]   # masked to land before averaging


def read(var):
    path = os.path.join(DATA_DIR, "{}_annual_{}.tif".format(var, YEAR))
    with rasterio.open(path) as src:
        a = src.read(1).astype("float64")
        a[(~np.isfinite(a)) | (a <= NODATA + 1)] = np.nan
        return a, src.transform, src.crs, src.bounds


# ---- 1. land mask from AOD (AOD is already water-masked in GEE) ----
aod, aod_tr, crs, b = read("AOD")
land = np.isfinite(aod).astype("float64")

# ---- 2. common grid ----
W = int(np.ceil((b.right - b.left) / CELL_DEG))
H = int(np.ceil((b.top - b.bottom) / CELL_DEG))
dst_tr = from_origin(b.left, b.top, CELL_DEG, CELL_DEG)


def to_grid(arr, tr):
    out = np.full((H, W), np.nan)
    reproject(arr, out, src_transform=tr, src_crs=crs, dst_transform=dst_tr,
              dst_crs=crs, src_nodata=np.nan, dst_nodata=np.nan,
              resampling=Resampling.average)
    return out


land_frac = to_grid(land, aod_tr)
good_cell = land_frac >= MIN_LAND

cols = {}
for var, (label, factor) in VARIABLES.items():
    a, tr, _, _ = read(var)
    a = a * factor
    if var in MODIS_1KM:
        if a.shape == land.shape:
            a[land == 0] = np.nan
        else:   # different grid: bring the land mask onto this grid first
            lm = np.zeros(a.shape)
            reproject(land, lm, src_transform=aod_tr, src_crs=crs,
                      dst_transform=tr, dst_crs=crs, resampling=Resampling.nearest)
            a[lm == 0] = np.nan
    g = to_grid(a, tr)
    g[~good_cell] = np.nan
    cols[label] = g.ravel()
    print("{:10s} valid cells: {}".format(var, int(np.isfinite(g).sum())))

rows, cols_idx = np.indices((H, W))
lon = (dst_tr.c + (cols_idx + 0.5) * CELL_DEG).ravel()
lat = (dst_tr.f - (rows + 0.5) * CELL_DEG).ravel()

df = pd.DataFrame(cols)
df.insert(0, "lat", lat)
df.insert(0, "lon", lon)
df = df.dropna().reset_index(drop=True)   # keep cells valid in ALL variables
n = len(df)
print("cells used (valid in all variables): {}".format(n))

if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR)
df.to_csv(os.path.join(OUT_DIR, "annual_pixels_{}.csv".format(YEAR)), index=False)

# ---- 3. correlation (Pearson and Spearman) ----
labels = [v[0] for v in VARIABLES.values()]
X = df[labels]

for method in ["pearson", "spearman"]:
    k = len(labels)
    r = np.zeros((k, k))
    p = np.zeros((k, k))
    for i in range(k):
        for j in range(k):
            if method == "pearson":
                r[i, j], p[i, j] = stats.pearsonr(X.iloc[:, i], X.iloc[:, j])
            else:
                r[i, j], p[i, j] = stats.spearmanr(X.iloc[:, i], X.iloc[:, j])

    pd.DataFrame(r, index=labels, columns=labels).round(3).to_csv(
        os.path.join(OUT_DIR, "corr_{}_r_{}.csv".format(method, YEAR)))
    pd.DataFrame(p, index=labels, columns=labels).to_csv(
        os.path.join(OUT_DIR, "corr_{}_p_{}.csv".format(method, YEAR)))

    # heatmap: lower triangle only
    show = np.where(np.tril(np.ones((k, k), bool)), r, np.nan)
    fig, ax = plt.subplots(figsize=(0.62 * k + 1.4, 0.55 * k + 1.0))
    im = ax.imshow(show, cmap="RdBu_r", vmin=-1, vmax=1)
    for i in range(k):
        for j in range(i + 1):
            star = "***" if p[i, j] < 0.001 else "**" if p[i, j] < 0.01 else "*" if p[i, j] < 0.05 else ""
            txt = "1" if i == j else "{:.2f}{}".format(r[i, j], star)
            ax.text(j, i, txt, ha="center", va="center", fontsize=7,
                    color="white" if abs(r[i, j]) > 0.6 else "black")
    ax.set_xticks(range(k))
    ax.set_yticks(range(k))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(labels, fontsize=8)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("{} r".format(method.capitalize()), fontsize=8)
    cb.ax.tick_params(labelsize=7)
    ax.set_title("{} correlation, annual {} (n = {} cells, {:.2f}\u00b0 grid)\n"
                 "* p<0.05  ** p<0.01  *** p<0.001".format(
                     method.capitalize(), YEAR, n, CELL_DEG), fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "corr_{}_{}.png".format(method, YEAR)), dpi=400)
    plt.close(fig)
    print("{}: saved".format(method))

print("Done.")

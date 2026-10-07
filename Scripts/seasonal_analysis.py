import matplotlib
matplotlib.use("Agg")

import os
import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon, Rectangle
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import stats

# ---------------- settings ----------------
YEAR      = 2025
NODATA    = -9999          # fill value written by the GEE export
COVER_MIN = 50             # seasons below this % coverage are labelled
PCT_LO    = 1              # shared colour range = 1st to 99th percentile
PCT_HI    = 99             #   of all valid pixels across the 4 seasons
MARGIN    = 0.08
CELL_DEG  = 0.05           # grid for the seasonal correlation (same as correlation_matrix.py)
MIN_LAND  = 0.5
MIN_CELLS = 30             # a correlation needs at least this many cells

ROOT     = r"C:\Users\User\Desktop\AirPollutionResearchWithRakibVai"
DATA_DIR = os.path.join(ROOT, "BD_AQ_2025")
BOUNDARY = os.path.join(ROOT, "boundary", "BGD_adm0.shp")
OUT_DIR  = os.path.join(ROOT, "Figures_v2", "Seasonal")
# ------------------------------------------

# Seasons of Bangladesh (month numbers). Winter uses Jan, Feb and Dec of the same year.
SEASONS = [
    ("Winter",       "Dec\u2013Feb", [12, 1, 2]),
    ("Pre-monsoon",  "Mar\u2013May", [3, 4, 5]),
    ("Monsoon",      "Jun\u2013Sep", [6, 7, 8, 9]),
    ("Post-monsoon", "Oct\u2013Nov", [10, 11]),
]
SEASON_COLOURS = ["#2c5f9e", "#e8955a", "#35978f", "#bf812d"]

O3_TO_DU = 1.0 / 4.4615e-4   # GEE O3 is mol/m2; 1 DU = 4.4615e-4 mol/m2
UMOL     = 1e6               # mol/m2 -> micromol/m2 (NO2, SO2)

CMAP_POLL = LinearSegmentedColormap.from_list(
    "greenbrown", ["#1a6b2e", "#7fb03a", "#dcd97a", "#c68b3c", "#7d3b1f"])
CMAP_TEMP = LinearSegmentedColormap.from_list(
    "bluered", ["#2c5f9e", "#7fb8dd", "#f2efc7", "#e8955a", "#a32020"])

# file prefix -> (label, unit, colormap, multiply factor)
VARIABLES = {
    "CH4":       (u"CH\u2084",   "ppb",                 CMAP_POLL, 1.0),
    "HCHO":      ("HCHO",        u"mol/m\u00b2",        CMAP_POLL, 1.0),
    "O3":        (u"O\u2083",    "DU",                  CMAP_POLL, O3_TO_DU),
    "UVAI":      ("UVAI",        "index",               CMAP_POLL, 1.0),
    "NO2":       (u"NO\u2082",   u"\u00b5mol/m\u00b2",  CMAP_POLL, UMOL),
    "CO":        ("CO",          u"mol/m\u00b2",        CMAP_POLL, 1.0),
    "SO2":       (u"SO\u2082",   u"\u00b5mol/m\u00b2",  CMAP_POLL, UMOL),
    "AOD":       ("AOD",         "550 nm",              CMAP_POLL, 1.0),
    "LST_Day":   ("LST Day",     u"\u00b0C",            CMAP_TEMP, 1.0),
    "LST_Night": ("LST Night",   u"\u00b0C",            CMAP_TEMP, 1.0),
}
LST_VARS  = ["LST_Day", "LST_Night"]
MODIS_1KM = ["LST_Day", "LST_Night"]     # masked to land before gridding


def fmt_val(v):
    a = abs(v)
    if a >= 100:
        return "{:.0f}".format(v)
    if a >= 1:
        return "{:.1f}".format(v)
    if a >= 0.01:
        return "{:.3f}".format(v)
    if a == 0:
        return "0"
    return "{:.2e}".format(v)


def dms(value, axis):
    hemi = ("N" if value >= 0 else "S") if axis == "y" else ("E" if value >= 0 else "W")
    v = abs(value)
    d = int(v)
    m = int(round((v - d) * 60))
    if m == 60:
        d, m = d + 1, 0
    return u"{}\u00b0{}'0\"{}".format(d, m, hemi)


def north_arrow(ax, x=0.90, y=0.83, size=0.055):
    t = ax.transAxes
    half = size * 0.30
    ax.add_patch(Polygon([[x, y + size], [x + half, y - size * 0.6], [x, y - size * 0.2]],
                         closed=True, facecolor="black", edgecolor="black",
                         lw=0.5, transform=t, zorder=6))
    ax.add_patch(Polygon([[x, y + size], [x - half, y - size * 0.6], [x, y - size * 0.2]],
                         closed=True, facecolor="white", edgecolor="black",
                         lw=0.5, transform=t, zorder=6))
    ax.text(x, y + size + 0.008, "N", transform=t, ha="center", va="bottom",
            fontsize=8, fontweight="bold", zorder=6)


def scale_bar(ax, ext, y=0.035, frac=0.34, n_seg=4):
    lat_mid = (ext[2] + ext[3]) / 2.0
    km_per_deg = 111.320 * np.cos(np.radians(lat_mid))
    span_km = (ext[1] - ext[0]) * (1 + 2 * MARGIN) * km_per_deg
    nice = np.array([25, 50, 100, 150, 200, 250])
    total_km = float(nice[np.argmin(np.abs(nice - span_km * frac))])
    bar_w = total_km / span_km
    seg_w = bar_w / n_seg
    bar_h = 0.010
    x = 0.5 - bar_w / 2.0
    t = ax.transAxes
    for k in range(n_seg):
        ax.add_patch(Rectangle((x + k * seg_w, y), seg_w, bar_h,
                               facecolor="black" if k % 2 == 0 else "white",
                               edgecolor="black", lw=0.5, transform=t, zorder=6))
    for k in range(0, n_seg + 1, 2):
        ax.text(x + k * seg_w, y + bar_h + 0.004, "{:g}".format(total_km * k / n_seg),
                transform=t, ha="center", va="bottom", fontsize=6, zorder=6)
    ax.text(x + bar_w + 0.010, y + bar_h * 0.4, "km", transform=t,
            ha="left", va="center", fontsize=6, zorder=6)


def read_months(var):
    """Return 12 monthly arrays (NaN = no data, raw units), transform, crs, bounds."""
    path = os.path.join(DATA_DIR, "{}_monthly_{}.tif".format(var, YEAR))
    if not os.path.exists(path):
        return None
    with rasterio.open(path) as src:
        a = src.read().astype("float64")
        a[(~np.isfinite(a)) | (a <= NODATA + 1)] = np.nan
        if src.nodata is not None:
            a[a == src.nodata] = np.nan
        return a, src.transform, src.crs, src.bounds


def season_mean(months, nums):
    """Pixel mean of the months of one season (months with no data are ignored)."""
    stack = months[[n - 1 for n in nums]]
    count = np.isfinite(stack).sum(axis=0)
    total = np.nansum(stack, axis=0)
    out = np.full(total.shape, np.nan)
    ok = count > 0
    out[ok] = total[ok] / count[ok]
    return out


# ---------------------------------------------------------------- main
if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR)

bnd = gpd.read_file(BOUNDARY)
season_names = [s[0] for s in SEASONS]

seasonal = {}      # var -> (list of 4 arrays in display units, transform, crs, bounds)
rows = []

# ============ 1. seasonal composites, table and maps ============
for var, (label, unit, cmap, factor) in VARIABLES.items():
    got = read_months(var)
    if got is None:
        print("{}: monthly raster not found - skipped".format(var))
        continue
    months, tr, crs, b = got
    if months.shape[0] != 12:
        print("{}: expected 12 bands, found {} - skipped".format(var, months.shape[0]))
        continue

    comps = [season_mean(months, nums) * factor for _, _, nums in SEASONS]
    seasonal[var] = (comps, tr, crs, b)
    ext = [b.left, b.right, b.bottom, b.top]

    valid_any = np.any(np.isfinite(months), axis=0)
    n_any = max(int(valid_any.sum()), 1)
    cover = []
    for (name, span, nums), a in zip(SEASONS, comps):
        c = 100.0 * np.isfinite(a).sum() / n_any
        cover.append(c)
        ok = np.isfinite(a)
        rows.append({"variable": var, "unit": unit, "season": name, "months": span,
                     "mean": np.nanmean(a) if ok.any() else np.nan,
                     "std": np.nanstd(a) if ok.any() else np.nan,
                     "min": np.nanmin(a) if ok.any() else np.nan,
                     "max": np.nanmax(a) if ok.any() else np.nan,
                     "coverage_%": c})

    pooled = np.concatenate([a[np.isfinite(a)] for a in comps])
    vmin, vmax = [float(v) for v in np.percentile(pooled, [PCT_LO, PCT_HI])]

    fig, axes = plt.subplots(2, 2, figsize=(8.27, 9.9))
    for i, (a, ax) in enumerate(zip(comps, axes.ravel())):
        im = ax.imshow(np.ma.masked_invalid(a), extent=ext, cmap=cmap, vmin=vmin, vmax=vmax, zorder=2)
        bnd.boundary.plot(ax=ax, edgecolor="black", linewidth=0.5, zorder=3)
        dx = (ext[1] - ext[0]) * MARGIN
        dy = (ext[3] - ext[2]) * MARGIN
        ax.set_xlim(ext[0] - dx, ext[1] + dx)
        ax.set_ylim(ext[2] - dy, ext[3] + dy)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: dms(v, "x")))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: dms(v, "y")))
        ax.tick_params(labelsize=6, length=2, pad=1.5)
        ax.tick_params(axis="y", labelrotation=90)
        ax.tick_params(axis="y", labelright=True, right=True)
        for s in ax.spines.values():
            s.set_linewidth(0.6)

        ax.text(0.5, 0.975, "{} ({}) {}".format(SEASONS[i][0], SEASONS[i][1], YEAR),
                transform=ax.transAxes, ha="center", va="top",
                fontsize=8.5, fontweight="bold", zorder=6,
                bbox=dict(boxstyle="square,pad=0.25", facecolor="white",
                          edgecolor="black", linewidth=0.5))
        if cover[i] < COVER_MIN:
            ax.text(0.5, 0.62, "Insufficient data\n(cloud, {:.0f}% coverage)".format(cover[i]),
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=7.5, fontweight="bold", zorder=7,
                    bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                              edgecolor="#a32020", linewidth=0.6, alpha=0.9))
        north_arrow(ax)
        scale_bar(ax, ext)

        cax = inset_axes(ax, width="4%", height="26%", loc="lower left",
                         bbox_to_anchor=(0.045, 0.15, 1, 1),
                         bbox_transform=ax.transAxes, borderpad=0)
        cb = fig.colorbar(im, cax=cax, extend="both")
        cb.set_ticks([vmin, vmax])
        cb.set_ticklabels([fmt_val(vmin), fmt_val(vmax)])
        cb.ax.tick_params(labelsize=6.5, pad=1.5, length=1.5)
        cb.ax.yaxis.set_ticks_position("right")
        cb.outline.set_linewidth(0.5)
        cb.ax.set_title("{}\n{}".format(label, unit), fontsize=7, fontweight="bold", pad=3)

    plt.subplots_adjust(left=0.07, right=0.975, top=0.985, bottom=0.03, wspace=0.20, hspace=0.08)
    plt.savefig(os.path.join(OUT_DIR, "{}_{}_seasonal.png".format(var, YEAR)), dpi=400)
    plt.close(fig)
    low = [n for n, c in zip(season_names, cover) if c < COVER_MIN]
    print("{:10s} map saved  (scale {} to {}{})".format(
        var, fmt_val(vmin), fmt_val(vmax), ", labelled: " + ", ".join(low) if low else ""))

if not rows:
    raise SystemExit("No monthly rasters found in " + DATA_DIR)

table = pd.DataFrame(rows)
table.to_csv(os.path.join(OUT_DIR, "seasonal_means_{}.csv".format(YEAR)), index=False)

# ============ 2. seasonal means figure ============
vars_done = list(seasonal.keys())
n_col = 5 if len(vars_done) > 5 else len(vars_done)
n_row = int(np.ceil(len(vars_done) / float(n_col)))
fig, axes = plt.subplots(n_row, n_col, figsize=(2.3 * n_col, 2.6 * n_row), squeeze=False)
x = np.arange(len(SEASONS))
for ax, var in zip(axes.ravel(), vars_done):
    label, unit, _, _ = VARIABLES[var]
    t = table[table["variable"] == var].set_index("season").reindex(season_names)
    f = 1e4 if var == "HCHO" else 1.0          # HCHO shown as x10^-4 mol/m2
    u = u"\u00d710\u207b\u2074 mol/m\u00b2" if var == "HCHO" else unit
    low = (t["coverage_%"] < COVER_MIN).values
    m, sd = t["mean"].values * f, t["std"].values * f
    ax.plot(x, m, "-", color="grey", lw=0.8, zorder=1)
    for k in range(len(x)):                     # hollow marker = coverage below the threshold
        ax.errorbar(x[k], m[k], yerr=sd[k], fmt="o", ms=6.5, color=SEASON_COLOURS[k],
                    mfc="white" if low[k] else SEASON_COLOURS[k], mec=SEASON_COLOURS[k],
                    ecolor=SEASON_COLOURS[k], elinewidth=1.0, capsize=3, zorder=2)
    ax.set_xlim(-0.5, len(x) - 0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(["Win", "Pre", "Mon", "Post"], fontsize=6.5)
    ax.set_title(u"{} ({})".format(label, u), fontsize=7.5, fontweight="bold")
    ax.tick_params(labelsize=6.5)
    ax.grid(axis="y", alpha=0.25, lw=0.4)
for ax in axes.ravel()[len(vars_done):]:
    ax.axis("off")
fig.suptitle("Seasonal national mean, Bangladesh {}  (dot = mean, bar = \u00b11 SD, "
             "hollow = coverage < {}%)".format(YEAR, COVER_MIN), fontsize=8)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(os.path.join(OUT_DIR, "seasonal_means_{}.png".format(YEAR)), dpi=400)
plt.close(fig)

# ============ 3. seasonal correlation: LST vs each pollutant ============
corr_rows = []
if "AOD" in seasonal and any(v in seasonal for v in LST_VARS):
    # land mask from the annual AOD raster (AOD is water-masked in GEE)
    with rasterio.open(os.path.join(DATA_DIR, "AOD_annual_{}.tif".format(YEAR))) as src:
        aod = src.read(1).astype("float64")
        aod[(~np.isfinite(aod)) | (aod <= NODATA + 1)] = np.nan
        aod_tr, crs0, b0 = src.transform, src.crs, src.bounds
    land = np.isfinite(aod).astype("float64")

    W = int(np.ceil((b0.right - b0.left) / CELL_DEG))
    H = int(np.ceil((b0.top - b0.bottom) / CELL_DEG))
    dst_tr = from_origin(b0.left, b0.top, CELL_DEG, CELL_DEG)

    def to_grid(arr, tr):
        out = np.full((H, W), np.nan)
        reproject(arr, out, src_transform=tr, src_crs=crs0, dst_transform=dst_tr,
                  dst_crs=crs0, src_nodata=np.nan, dst_nodata=np.nan,
                  resampling=Resampling.average)
        return out

    good_cell = to_grid(land, aod_tr) >= MIN_LAND

    grid = {}
    for var, (comps, tr, _, _) in seasonal.items():
        for k, a in enumerate(comps):
            a = a.copy()
            if var in MODIS_1KM:
                if a.shape == land.shape:
                    a[land == 0] = np.nan
                else:
                    lm = np.zeros(a.shape)
                    reproject(land, lm, src_transform=aod_tr, src_crs=crs0,
                              dst_transform=tr, dst_crs=crs0, resampling=Resampling.nearest)
                    a[lm == 0] = np.nan
            g = to_grid(a, tr)
            g[~good_cell] = np.nan
            grid[(var, k)] = g.ravel()

    pollutants = [v for v in seasonal if v not in LST_VARS]
    for lst in [v for v in LST_VARS if v in seasonal]:
        for pol in pollutants:
            for k, name in enumerate(season_names):
                xg, yg = grid[(pol, k)], grid[(lst, k)]
                ok = np.isfinite(xg) & np.isfinite(yg)
                n = int(ok.sum())
                if n >= MIN_CELLS and np.std(xg[ok]) > 0 and np.std(yg[ok]) > 0:
                    r, p = stats.pearsonr(xg[ok], yg[ok])
                    rho, p_s = stats.spearmanr(xg[ok], yg[ok])
                else:
                    r = p = rho = p_s = np.nan
                corr_rows.append({"LST": VARIABLES[lst][0], "pollutant": VARIABLES[pol][0],
                                  "season": name, "n_cells": n, "pearson_r": r, "pearson_p": p,
                                  "spearman_rho": rho, "spearman_p": p_s})

    corr = pd.DataFrame(corr_rows)
    corr.to_csv(os.path.join(OUT_DIR, "seasonal_correlation_{}.csv".format(YEAR)), index=False)

    # heatmap: rows = pollutants, columns = seasons, one panel per LST
    lsts = list(corr["LST"].unique())
    pols = [VARIABLES[v][0] for v in pollutants]
    fig, axes = plt.subplots(1, len(lsts), figsize=(4.3 * len(lsts) + 0.8, 0.42 * len(pols) + 1.6),
                             squeeze=False)
    for ax, lst in zip(axes.ravel(), lsts):
        c = corr[corr["LST"] == lst]
        R = c.pivot(index="pollutant", columns="season", values="pearson_r").reindex(
            index=pols, columns=season_names)
        P = c.pivot(index="pollutant", columns="season", values="pearson_p").reindex(
            index=pols, columns=season_names)
        im = ax.imshow(np.ma.masked_invalid(R.values), cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
        for i in range(len(pols)):
            for j in range(len(season_names)):
                r, p = R.values[i, j], P.values[i, j]
                if not np.isfinite(r):
                    ax.text(j, i, "n/a", ha="center", va="center", fontsize=7, color="grey")
                    continue
                star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
                ax.text(j, i, "{:.2f}{}".format(r, star), ha="center", va="center", fontsize=7,
                        color="white" if abs(r) > 0.6 else "black")
        ax.set_xticks(range(len(season_names)))
        ax.set_xticklabels(season_names, fontsize=7.5)
        ax.set_yticks(range(len(pols)))
        ax.set_yticklabels(pols, fontsize=8)
        ax.set_title(lst, fontsize=9, fontweight="bold")
        for s in ax.spines.values():
            s.set_visible(False)
    cb = fig.colorbar(im, ax=axes.ravel().tolist(), fraction=0.03, pad=0.02)
    cb.set_label("Pearson r", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    fig.suptitle("Seasonal correlation between LST and pollutants, {} ({:.2f}\u00b0 grid)   "
                 "* p<0.05  ** p<0.01  *** p<0.001".format(YEAR, CELL_DEG), fontsize=8)
    plt.savefig(os.path.join(OUT_DIR, "seasonal_correlation_{}.png".format(YEAR)),
                dpi=400, bbox_inches="tight")
    plt.close(fig)

# ---- print ----
print("\nSeasonal means:")
print(table.pivot(index="variable", columns="season", values="mean")
      .reindex(index=vars_done, columns=season_names).to_string(float_format=lambda v: fmt_val(v)))
print("\nCoverage (%):")
print(table.pivot(index="variable", columns="season", values="coverage_%")
      .reindex(index=vars_done, columns=season_names).round(0).to_string())
if corr_rows:
    print("\nSeasonal Pearson r (LST Day):")
    c = corr[corr["LST"] == "LST Day"]
    print(c.pivot(index="pollutant", columns="season", values="pearson_r")
          .reindex(columns=season_names).round(2).to_string())
print("\nsaved to: " + OUT_DIR)
print("Done.")

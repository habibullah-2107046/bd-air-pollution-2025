import matplotlib
matplotlib.use("Agg")

import os
import numpy as np
import pandas as pd
import rasterio
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon, Rectangle
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# ---------------- settings ----------------
YEAR      = 2025
NODATA    = -9999          # fill value written by the GEE export
COVER_MIN = 50             # months below this % coverage get the monsoon label
PCT_LO    = 1              # shared colour range = 1st to 99th percentile
PCT_HI    = 99             #   of all valid pixels across the 12 months
MARGIN    = 0.08           # padding around the map (fraction of extent)

A4_W, A4_H = 8.27, 11.69   # inches, portrait

ROOT     = r"C:\Users\User\Desktop\AirPollutionResearchWithRakibVai"
DATA_DIR = os.path.join(ROOT, "BD_AQ_2025")
BOUNDARY = os.path.join(ROOT, "boundary", "BGD_adm0.shp")
COVERAGE = os.path.join(DATA_DIR, "coverage_monthly_{}.csv".format(YEAR))
FIG_DIR  = os.path.join(ROOT, "Figures_v2", "Monthly")
# ------------------------------------------

O3_TO_DU = 1.0 / 4.4615e-4   # GEE O3 is mol/m2; 1 DU = 4.4615e-4 mol/m2
UMOL     = 1e6                # mol/m2 -> micromol/m2 (NO2, SO2)

MONTHS = ["January", "February", "March", "April", "May", "June",
          "July", "August", "September", "October", "November", "December"]

CMAP_POLL = LinearSegmentedColormap.from_list(
    "greenbrown", ["#1a6b2e", "#7fb03a", "#dcd97a", "#c68b3c", "#7d3b1f"])
CMAP_TEMP = LinearSegmentedColormap.from_list(
    "bluered", ["#2c5f9e", "#7fb8dd", "#f2efc7", "#e8955a", "#a32020"])

# file prefix -> (legend label, unit, colormap, multiply factor)
VARIABLES = {
    "CH4":       (u"CH\u2084",   "ppb",           CMAP_POLL, 1.0),
    "HCHO":      ("HCHO",        u"mol/m\u00b2",  CMAP_POLL, 1.0),
    "O3":        (u"O\u2083",    "DU",            CMAP_POLL, O3_TO_DU),
    "UVAI":      ("UVAI",        "index",         CMAP_POLL, 1.0),
    "NO2":       (u"NO\u2082",   u"\u00b5mol/m\u00b2", CMAP_POLL, UMOL),
    "CO":        ("CO",          u"mol/m\u00b2",  CMAP_POLL, 1.0),
    "SO2":       (u"SO\u2082",   u"\u00b5mol/m\u00b2", CMAP_POLL, UMOL),
    "AOD":       ("AOD",         "550 nm",        CMAP_POLL, 1.0),
    "LST_Day":   ("LST Day",     u"\u00b0C",      CMAP_TEMP, 1.0),
    "LST_Night": ("LST Night",   u"\u00b0C",      CMAP_TEMP, 1.0),
}


def read_band(src, band, factor):
    """Read one band, mask the -9999 fill and NaN, apply unit factor."""
    a = src.read(band).astype("float64")
    bad = ~np.isfinite(a) | (a <= NODATA + 1)
    if src.nodata is not None:
        bad |= (a == src.nodata)
    return np.ma.array(a * factor, mask=bad)


def load_coverage(path):
    """Return {(variable, month): pct_covered}; empty dict if file missing."""
    if not os.path.exists(path):
        print("coverage file not found - no months will be labelled")
        return {}
    df = pd.read_csv(path)
    return {(str(r["variable"]), int(r["month"])): float(r["pct_covered"])
            for _, r in df.iterrows()}


def fmt_val(v):
    """Colourbar tick text suited to both large and very small values."""
    a = abs(v)
    if a >= 100:
        return "{:.0f}".format(v)
    if a >= 1:
        return "{:.1f}".format(v)
    if a >= 0.01:
        return "{:.3f}".format(v)
    return "{:.2e}".format(v)


def dms(value, axis):
    """Format a decimal degree as e.g. 88 deg 0'0\"E."""
    hemi = ("N" if value >= 0 else "S") if axis == "y" else ("E" if value >= 0 else "W")
    v = abs(value)
    d = int(v)
    m = int(round((v - d) * 60))
    if m == 60:
        d, m = d + 1, 0
    return u"{}\u00b0{}'0\"{}".format(d, m, hemi)


def north_arrow(ax, x=0.90, y=0.83, size=0.055):
    """Small split-needle compass, drawn in axes coordinates."""
    t = ax.transAxes
    half = size * 0.30
    ax.add_patch(Polygon([[x, y + size], [x + half, y - size * 0.6], [x, y - size * 0.2]],
                         closed=True, facecolor="black", edgecolor="black",
                         lw=0.5, transform=t, zorder=6))
    ax.add_patch(Polygon([[x, y + size], [x - half, y - size * 0.6], [x, y - size * 0.2]],
                         closed=True, facecolor="white", edgecolor="black",
                         lw=0.5, transform=t, zorder=6))
    ax.text(x, y + size + 0.008, "N", transform=t, ha="center", va="bottom",
            fontsize=6, fontweight="bold", zorder=6)


def scale_bar(ax, ext, y=0.035, frac=0.34, n_seg=4):
    """Segmented cartographic scale bar, horizontally centred."""
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
                transform=t, ha="center", va="bottom", fontsize=4.8, zorder=6)
    ax.text(x + bar_w + 0.010, y + bar_h * 0.4, "km", transform=t,
            ha="left", va="center", fontsize=4.8, zorder=6)


# ---------------------------------------------------------------- main
if not os.path.exists(FIG_DIR):
    os.makedirs(FIG_DIR)

bnd = gpd.read_file(BOUNDARY)
cover = {}   # filled per variable from the raster pixels

for var, (label, unit, cmap, factor) in VARIABLES.items():

    path = os.path.join(DATA_DIR, "{}_monthly_{}.tif".format(var, YEAR))
    if not os.path.exists(path):
        print("{}: {} not found - skipping".format(var, os.path.basename(path)))
        continue

    with rasterio.open(path) as src:
        if src.count != 12:
            print("{}: expected 12 bands, found {} - skipping".format(var, src.count))
            continue
        months = [read_band(src, b, factor) for b in range(1, 13)]
        ext = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]

    # coverage from the raster itself: % of pixels valid this month,
    # out of all pixels valid in at least one month (so water/no-data areas don't count)
    valid_any = np.zeros(months[0].shape, dtype=bool)
    for m in months:
        valid_any |= ~np.ma.getmaskarray(m)
    n_any = max(int(valid_any.sum()), 1)
    for k, m in enumerate(months):
        cover[(var, k + 1)] = 100.0 * (~np.ma.getmaskarray(m) & valid_any).sum() / n_any

    # one shared colour range for all 12 months
    pooled = np.concatenate([m.compressed() for m in months])
    vmin, vmax = np.percentile(pooled, [PCT_LO, PCT_HI])
    vmin, vmax = float(vmin), float(vmax)

    fig, axes = plt.subplots(4, 3, figsize=(A4_W, A4_H))

    for i, (arr, ax) in enumerate(zip(months, axes.ravel())):
        im = ax.imshow(arr, extent=ext, cmap=cmap, vmin=vmin, vmax=vmax, zorder=2)
        bnd.boundary.plot(ax=ax, edgecolor="black", linewidth=0.4, zorder=3)

        dx = (ext[1] - ext[0]) * MARGIN
        dy = (ext[3] - ext[2]) * MARGIN
        ax.set_xlim(ext[0] - dx, ext[1] + dx)
        ax.set_ylim(ext[2] - dy, ext[3] + dy)

        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: dms(v, "x")))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: dms(v, "y")))
        ax.tick_params(labelsize=4.5, length=2, pad=1.5)
        ax.tick_params(axis="y", labelrotation=90)
        ax.tick_params(axis="y", labelright=True, right=True)
        for s in ax.spines.values():
            s.set_linewidth(0.6)

        ax.text(0.5, 0.975, "{} {}".format(MONTHS[i], YEAR),
                transform=ax.transAxes, ha="center", va="top",
                fontsize=7, fontweight="bold", zorder=6,
                bbox=dict(boxstyle="square,pad=0.25", facecolor="white",
                          edgecolor="black", linewidth=0.5))

        # label months below the coverage threshold (still drawn, as agreed)
        pct = cover.get((var, i + 1))
        if pct is not None and pct < COVER_MIN:
            ax.text(0.5, 0.62,
                    "Insufficient data\n(monsoon cloud, {:.0f}% coverage)".format(pct),
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=5.5, fontweight="bold", zorder=7,
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
        cb.ax.tick_params(labelsize=5, pad=1.5, length=1.5)
        cb.ax.yaxis.set_ticks_position("right")
        cb.outline.set_linewidth(0.5)
        cb.ax.set_title("{}\n{}".format(label, unit), fontsize=5.5,
                        fontweight="bold", pad=3)

    plt.subplots_adjust(left=0.085, right=0.985, top=0.985, bottom=0.035,
                        wspace=0.28, hspace=0.10)
    out_png = os.path.join(FIG_DIR, "{}_{}_panel.png".format(var, YEAR))
    plt.savefig(out_png, dpi=400)
    plt.close(fig)

    low = [MONTHS[m - 1][:3] for (v, m), p in sorted(cover.items())
           if v == var and p < COVER_MIN]
    print("{}: saved  (scale {} to {}{})".format(
        var, fmt_val(vmin), fmt_val(vmax),
        ", labelled: " + ", ".join(low) if low else ""))

print("Done.")

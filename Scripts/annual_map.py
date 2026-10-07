import matplotlib
matplotlib.use("Agg")

import os
import numpy as np
import rasterio
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon, Rectangle
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# ---------------- settings ----------------
YEAR    = 2025
NODATA  = -9999            # fill value written by the GEE export
PCT_LO  = 1                # colour range = 1st to 99th percentile of valid pixels
PCT_HI  = 99
MARGIN  = 0.06             # padding around the map (fraction of extent)

AX_W = 6.2                 # width of the map frame, inches
PAD  = 0.85                # margin on every side, inches

ROOT     = r"C:\Users\User\Desktop\AirPollutionResearchWithRakibVai"
DATA_DIR = os.path.join(ROOT, "BD_AQ_2025")
BOUNDARY = os.path.join(ROOT, "boundary", "BGD_adm0.shp")
FIG_DIR  = os.path.join(ROOT, "Figures_v2", "Annual")
# ------------------------------------------

O3_TO_DU = 1.0 / 4.4615e-4   # GEE O3 is mol/m2; 1 DU = 4.4615e-4 mol/m2
UMOL     = 1e6                # mol/m2 -> micromol/m2 (NO2, SO2)

CMAP_POLL = LinearSegmentedColormap.from_list(
    "greenbrown", ["#1a6b2e", "#7fb03a", "#dcd97a", "#c68b3c", "#7d3b1f"])
CMAP_TEMP = LinearSegmentedColormap.from_list(
    "bluered", ["#2c5f9e", "#7fb8dd", "#f2efc7", "#e8955a", "#a32020"])

# file prefix -> (legend label, unit, colormap, title, multiply factor)
VARIABLES = {
    "CH4":       (u"CH\u2084", "ppb",          CMAP_POLL, u"Mean annual CH\u2084",   1.0),
    "HCHO":      ("HCHO",      u"mol/m\u00b2", CMAP_POLL, "Mean annual HCHO",        1.0),
    "O3":        (u"O\u2083",  "DU",           CMAP_POLL, u"Mean annual O\u2083",    O3_TO_DU),
    "UVAI":      ("UVAI",      "index",        CMAP_POLL, "Mean annual UVAI",        1.0),
    "NO2":       (u"NO\u2082", u"\u00b5mol/m\u00b2", CMAP_POLL, u"Mean annual NO\u2082", UMOL),
    "CO":        ("CO",        u"mol/m\u00b2", CMAP_POLL, "Mean annual CO",          1.0),
    "SO2":       (u"SO\u2082", u"\u00b5mol/m\u00b2", CMAP_POLL, u"Mean annual SO\u2082", UMOL),
    "AOD":       ("AOD",       "550 nm",       CMAP_POLL, "Annual median AOD",       1.0),
    "LST_Day":   ("LST",       u"\u00b0C",     CMAP_TEMP, "Mean annual LST (Day)",   1.0),
    "LST_Night": ("LST",       u"\u00b0C",     CMAP_TEMP, "Mean annual LST (Night)", 1.0),
    "UHI":       ("UHI",       u"\u00b0C",     CMAP_TEMP, "Annual UHI intensity",    1.0),
}


def load(path, factor):
    """Read band 1, mask the -9999 fill and NaN, apply unit factor."""
    with rasterio.open(path) as src:
        a = src.read(1).astype("float64")
        bad = ~np.isfinite(a) | (a <= NODATA + 1)
        if src.nodata is not None:
            bad |= (a == src.nodata)
        ext = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]
    return np.ma.array(a * factor, mask=bad), ext


def fmt_val(v):
    """Colourbar tick text suited to both large and very small values."""
    a = abs(v)
    if a >= 100:
        return "{:.0f}".format(v)
    if a >= 1:
        return "{:.2f}".format(v)
    if a >= 0.01:
        return "{:.3f}".format(v)
    if a == 0:
        return "0"
    return "{:.2e}".format(v)


def tick_labels(ticks):
    """Same number of decimals for every tick of one colourbar."""
    top = float(np.max(np.abs(ticks)))
    step = abs(float(ticks[1] - ticks[0]))
    if top < 0.01:
        return [fmt_val(t) for t in ticks]          # very small values: scientific
    if top >= 100:
        dec = 0 if step >= 5 else 1
    elif top >= 10:
        dec = 1 if step >= 1 else 2
    elif top >= 1:
        dec = 2
    else:
        dec = 3 if step >= 0.005 else 4
    out = ["{:.{d}f}".format(t, d=dec) for t in ticks]
    return ["0" if float(o) == 0 else o for o in out]


def dms(value, axis):
    """Format a decimal degree as e.g. 88 deg 0'0\"E."""
    hemi = ("N" if value >= 0 else "S") if axis == "y" else ("E" if value >= 0 else "W")
    v = abs(value)
    d = int(v)
    m = int(round((v - d) * 60))
    if m == 60:
        d, m = d + 1, 0
    return u"{}\u00b0{}'0\"{}".format(d, m, hemi)


def north_arrow(ax, x=0.92, y=0.78, size=0.060):
    """Split-needle compass, drawn in axes coordinates."""
    t = ax.transAxes
    half = size * 0.30
    ax.add_patch(Polygon([[x, y + size], [x + half, y - size * 0.6], [x, y - size * 0.2]],
                         closed=True, facecolor="black", edgecolor="black",
                         lw=0.7, transform=t, zorder=6))
    ax.add_patch(Polygon([[x, y + size], [x - half, y - size * 0.6], [x, y - size * 0.2]],
                         closed=True, facecolor="white", edgecolor="black",
                         lw=0.7, transform=t, zorder=6))
    ax.text(x, y + size + 0.010, "N", transform=t, ha="center", va="bottom",
            fontsize=10, fontweight="bold", zorder=6)


def scale_bar(ax, ext, y=0.045, frac=0.30, n_seg=4):
    """Segmented cartographic scale bar, horizontally centred."""
    lat_mid = (ext[2] + ext[3]) / 2.0
    km_per_deg = 111.320 * np.cos(np.radians(lat_mid))
    span_km = (ext[1] - ext[0]) * (1 + 2 * MARGIN) * km_per_deg

    nice = np.array([25, 50, 100, 150, 200, 250])
    total_km = float(nice[np.argmin(np.abs(nice - span_km * frac))])

    bar_w = total_km / span_km
    seg_w = bar_w / n_seg
    bar_h = 0.013
    x = 0.5 - bar_w / 2.0
    t = ax.transAxes

    for k in range(n_seg):
        ax.add_patch(Rectangle((x + k * seg_w, y), seg_w, bar_h,
                               facecolor="black" if k % 2 == 0 else "white",
                               edgecolor="black", lw=0.7, transform=t, zorder=6))
    for k in range(n_seg + 1):
        ax.text(x + k * seg_w, y + bar_h + 0.007, "{:g}".format(total_km * k / n_seg),
                transform=t, ha="center", va="bottom", fontsize=7, zorder=6)
    ax.text(x + bar_w + 0.014, y + bar_h * 0.4, "km", transform=t,
            ha="left", va="center", fontsize=7, zorder=6)


# ---------------------------------------------------------------- main
if not os.path.exists(FIG_DIR):
    os.makedirs(FIG_DIR)

bnd = gpd.read_file(BOUNDARY)

for var, (label, unit, cmap, title, factor) in VARIABLES.items():

    path = os.path.join(DATA_DIR, "{}_annual_{}.tif".format(var, YEAR))
    if not os.path.exists(path):
        print("{}: {} not found - skipping".format(var, os.path.basename(path)))
        continue

    arr, ext = load(path, factor)
    vmin, vmax = [float(v) for v in np.percentile(arr.compressed(), [PCT_LO, PCT_HI])]
    if var == "UHI":                      # symmetric around 0 so white = no UHI
        lim = max(abs(vmin), abs(vmax))
        vmin, vmax = -lim, lim

    dx = (ext[1] - ext[0]) * MARGIN
    dy = (ext[3] - ext[2]) * MARGIN
    x0, x1 = ext[0] - dx, ext[1] + dx
    y0, y1 = ext[2] - dy, ext[3] + dy

    ax_h = AX_W * (y1 - y0) / (x1 - x0)
    fig_w = AX_W + 2 * PAD
    fig_h = ax_h + 2 * PAD

    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([PAD / fig_w, PAD / fig_h, AX_W / fig_w, ax_h / fig_h])

    im = ax.imshow(arr, extent=ext, cmap=cmap, vmin=vmin, vmax=vmax, zorder=2)
    bnd.boundary.plot(ax=ax, edgecolor="black", linewidth=0.9, zorder=3)

    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("auto")

    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: dms(v, "x")))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: dms(v, "y")))
    ax.tick_params(labelsize=8, length=3, pad=2)
    ax.tick_params(axis="y", labelrotation=90)
    ax.tick_params(axis="y", labelright=True, right=True)
    ax.tick_params(axis="x", labeltop=True, top=True)
    for s in ax.spines.values():
        s.set_linewidth(0.9)

    ax.text(0.5, 0.975, "{} {}".format(title, YEAR),
            transform=ax.transAxes, ha="center", va="top",
            fontsize=10, fontweight="bold", zorder=6,
            bbox=dict(boxstyle="square,pad=0.35", facecolor="white",
                      edgecolor="black", linewidth=0.7))

    north_arrow(ax)
    scale_bar(ax, ext)

    cax = inset_axes(ax, width="4.5%", height="30%", loc="lower left",
                     bbox_to_anchor=(0.055, 0.15, 1, 1),
                     bbox_transform=ax.transAxes, borderpad=0)
    cb = fig.colorbar(im, cax=cax, extend="both")
    ticks = np.linspace(vmin, vmax, 5)
    cb.set_ticks(ticks)
    cb.set_ticklabels(tick_labels(ticks))
    cb.ax.tick_params(labelsize=7.5, pad=2, length=2)
    cb.ax.yaxis.set_ticks_position("right")
    cb.outline.set_linewidth(0.7)
    cb.ax.set_title("{}\n{}".format(label, unit), fontsize=8.5, fontweight="bold", pad=5)

    out_png = os.path.join(FIG_DIR, "{}_{}_annual.png".format(var, YEAR))
    plt.savefig(out_png, dpi=400)
    plt.close(fig)
    print("{}: saved  (scale {} to {})".format(var, fmt_val(vmin), fmt_val(vmax)))

print("Done.")

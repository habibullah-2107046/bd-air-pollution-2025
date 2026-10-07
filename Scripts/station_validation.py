"""
station_validation.py
Compare satellite monthly rasters with DoE CAMS station data (2025).

Inputs  (all in BD_AQ_2025):
    doe_cams_monthly_2025.csv        station monthly values (16 stations)
    <VAR>_monthly_2025.tif           12-band monthly rasters (AOD, UVAI, NO2, SO2, CO, O3)

Outputs (in Station_validation):
    station_validation_pairs.csv     every station-month pair (ground value + satellite value)
    station_validation_summary.csv   Pearson / Spearman for each pair
    station_who_comparison.csv       station PM2.5 / PM10 means vs WHO and Bangladesh standards
    station_validation_scatter.png   scatter plots

A pair is skipped automatically if its raster is not in the folder yet,
so you can run this now and again after NO2 / CO / SO2 are downloaded.
"""
import os
import numpy as np
import pandas as pd
import rasterio
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = r"C:\Users\User\Desktop\AirPollutionResearchWithRakibVai"
IN_DIR = os.path.join(ROOT, "BD_AQ_2025")
OUT_DIR = os.path.join(ROOT, "Station_validation")
STATION_CSV = os.path.join(IN_DIR, "doe_cams_monthly_2025.csv")

MIN_CAPTURE = 50      # use station-months with at least 50% data capture
WINDOW = 1            # 1 = 3x3 pixels around the station (about 3 km)
NODATA = -9999

# (ground parameter, satellite variable)
PAIRS = [
    ("PM2.5", "AOD"),
    ("PM10",  "AOD"),
    ("PM2.5", "UVAI"),
    ("NO2",   "NO2"),
    ("SO2",   "SO2"),
    ("CO",    "CO"),
    ("O3",    "O3"),
]

# annual limits in ug/m3: WHO 2021 guideline, Bangladesh standard
LIMITS = {"PM2.5": (5, 35), "PM10": (15, 50)}

MONTHS = ["January", "February", "March", "April", "May", "June", "July",
          "August", "September", "October", "November", "December"]


def month_number(x):
    """Accept 1..12, '1', 'Jan', 'January', '2025-01' ..."""
    try:
        return int(float(x))
    except (TypeError, ValueError):
        s = str(x).strip()
        for i, name in enumerate(MONTHS, 1):
            if s.lower().startswith(name[:3].lower()):
                return i
        return pd.to_datetime(s).month


def sample_raster(path, lats, lons):
    """Mean of a (2*WINDOW+1)^2 window around each point, for all 12 bands.
    Returns array [n_points, 12] with NaN where there is no data."""
    out = np.full((len(lats), 12), np.nan)
    with rasterio.open(path) as src:
        data = src.read().astype("float64")
        nod = src.nodata
        data[data == NODATA] = np.nan
        if nod is not None:
            data[data == nod] = np.nan
        nb = min(12, data.shape[0])
        for i, (lat, lon) in enumerate(zip(lats, lons)):
            r, c = src.index(lon, lat)
            if r < 0 or c < 0 or r >= src.height or c >= src.width:
                continue
            r0, r1 = max(r - WINDOW, 0), min(r + WINDOW + 1, src.height)
            c0, c1 = max(c - WINDOW, 0), min(c + WINDOW + 1, src.width)
            win = data[:nb, r0:r1, c0:c1].reshape(nb, -1)
            ok = np.isfinite(win).any(axis=1)
            out[i, :nb][ok] = np.nanmean(win[ok], axis=1)
    return out


def corr(x, y):
    if len(x) < 4 or np.std(x) == 0 or np.std(y) == 0:
        return dict(pearson_r=np.nan, pearson_p=np.nan, spearman_rho=np.nan, spearman_p=np.nan)
    pr, pp = stats.pearsonr(x, y)
    sr, sp = stats.spearmanr(x, y)
    return dict(pearson_r=pr, pearson_p=pp, spearman_rho=sr, spearman_p=sp)


os.makedirs(OUT_DIR, exist_ok=True)

# ---------- station data
d = pd.read_csv(STATION_CSV)
d["month"] = d["month"].map(month_number)
d["mean"] = pd.to_numeric(d["mean"], errors="coerce")
d = d[(d["capture_pct"] >= MIN_CAPTURE) & d["mean"].notna()].copy()
stations = d[["station", "lat", "lon"]].drop_duplicates("station").reset_index(drop=True)
print("stations: {}   usable station-months: {}".format(len(stations), len(d)))

# ---------- sample each satellite variable at the stations
sat = {}
for var in sorted(set(v for _, v in PAIRS)):
    path = os.path.join(IN_DIR, "{}_monthly_2025.tif".format(var))
    if not os.path.exists(path):
        print("  {:5s} raster not found -> skipped".format(var))
        continue
    vals = sample_raster(path, stations["lat"].values, stations["lon"].values)
    t = pd.DataFrame(vals, columns=range(1, 13))
    t["station"] = stations["station"]
    sat[var] = t.melt(id_vars="station", var_name="month", value_name="sat_value")
    print("  {:5s} sampled".format(var))

# ---------- build pairs
rows, summary = [], []
for par, var in PAIRS:
    if var not in sat:
        continue
    g = d[d["parameter"] == par][["station", "lat", "lon", "month", "mean"]]
    g = g.rename(columns={"mean": "ground_value"})
    j = g.merge(sat[var], on=["station", "month"]).dropna(subset=["ground_value", "sat_value"])
    if j.empty:
        continue
    j.insert(0, "pair", "{} vs {}".format(par, var))
    j["ground_parameter"], j["sat_variable"] = par, var
    rows.append(j)

    # monthly pairs (all station-months)
    s = dict(pair="{} vs {}".format(par, var), level="station-month", n=len(j))
    s.update(corr(j["ground_value"].values, j["sat_value"].values))
    summary.append(s)

    # station means (spatial agreement)
    a = j.groupby("station")[["ground_value", "sat_value"]].mean()
    s = dict(pair="{} vs {}".format(par, var), level="station mean", n=len(a))
    s.update(corr(a["ground_value"].values, a["sat_value"].values))
    summary.append(s)

if not rows:
    raise SystemExit("No pairs found. Check the rasters and the station CSV are in " + IN_DIR)

pairs = pd.concat(rows, ignore_index=True)
summary = pd.DataFrame(summary)
pairs.to_csv(os.path.join(OUT_DIR, "station_validation_pairs.csv"), index=False)
summary.round(4).to_csv(os.path.join(OUT_DIR, "station_validation_summary.csv"), index=False)
print("\n" + summary.round(3).to_string(index=False))

# ---------- WHO / Bangladesh standard comparison (ground data)
who = []
for par, (who_lim, bd_lim) in LIMITS.items():
    g = d[d["parameter"] == par].groupby("station")["mean"].agg(["mean", "count"])
    for st, r in g.iterrows():
        who.append(dict(station=st, parameter=par, mean_ug_m3=round(r["mean"], 1),
                        months_used=int(r["count"]),
                        WHO_guideline=who_lim, times_WHO=round(r["mean"] / who_lim, 1),
                        BD_standard=bd_lim, times_BD=round(r["mean"] / bd_lim, 1)))
who = pd.DataFrame(who)
who.to_csv(os.path.join(OUT_DIR, "station_who_comparison.csv"), index=False)
for par, (who_lim, bd_lim) in LIMITS.items():
    w = who[who["parameter"] == par]
    print("\n{}: {} of {} stations above WHO ({}), {} above Bangladesh standard ({}); range {}-{} ug/m3".format(
        par, (w["mean_ug_m3"] > who_lim).sum(), len(w), who_lim,
        (w["mean_ug_m3"] > bd_lim).sum(), bd_lim, w["mean_ug_m3"].min(), w["mean_ug_m3"].max()))

# ---------- scatter figure
names = list(pairs["pair"].unique())
ncol = min(3, len(names))
nrow = int(np.ceil(len(names) / ncol))
plt.rcParams["font.family"] = ["Times New Roman", "DejaVu Serif"]
fig, axes = plt.subplots(nrow, ncol, figsize=(4.4 * ncol, 3.9 * nrow), squeeze=False,
                         constrained_layout=True)
for ax, name in zip(axes.ravel(), names):
    j = pairs[pairs["pair"] == name]
    x, y = j["sat_value"].values, j["ground_value"].values
    sc = ax.scatter(x, y, c=j["month"], cmap="twilight", vmin=1, vmax=12, s=22,
                    edgecolor="k", linewidth=0.3)
    c = corr(x, y)
    if np.isfinite(c["pearson_r"]):
        k, b = np.polyfit(x, y, 1)
        xx = np.linspace(x.min(), x.max(), 50)
        ax.plot(xx, k * xx + b, "r-", lw=1)
        ax.text(0.04, 0.96, "r = {:.2f}\nn = {}".format(c["pearson_r"], len(j)),
                transform=ax.transAxes, va="top", fontsize=10)
    ax.set_xlabel("Satellite " + j["sat_variable"].iloc[0])
    ax.set_ylabel("Station " + j["ground_parameter"].iloc[0])
    ax.set_title(name, fontsize=11)
    ax.grid(alpha=0.3)
for ax in axes.ravel()[len(names):]:
    ax.axis("off")
cb = fig.colorbar(sc, ax=axes.ravel().tolist(), shrink=0.6, ticks=range(1, 13))
cb.set_label("Month")
fig.savefig(os.path.join(OUT_DIR, "station_validation_scatter.png"), dpi=300, bbox_inches="tight")

print("\nsaved to: " + OUT_DIR)
print("Done.")

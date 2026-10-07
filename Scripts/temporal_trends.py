import matplotlib
matplotlib.use("Agg")

import os
import numpy as np
import pandas as pd
import rasterio
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# ---------------- settings ----------------
YEAR      = 2025
NODATA    = -9999
COVER_MIN = 50           # AOD days / months below this % coverage are treated as unreliable
ROLL_DAYS = 15           # smoothing window for the daily line
OUTLIER_K = 3            # drop days more than K robust SDs from the local median
ROOT      = r"C:\Users\User\Desktop\AirPollutionResearchWithRakibVai"
DATA_DIR  = os.path.join(ROOT, "BD_AQ_2025")
OUT_DIR   = os.path.join(ROOT, "Figures_v2", "Temporal")
# ------------------------------------------

O3_TO_DU = 1.0 / 4.4615e-4

# file prefix -> (axis label, factor, colour)
VARIABLES = {
    "AOD":       ("AOD (550 nm)",                              1.0,      "#8c510a"),
    "UVAI":      ("UVAI (index)",                              1.0,      "#bf812d"),
    "CH4":       (u"CH\u2084 (ppb)",                           1.0,      "#35978f"),
    "HCHO":      (u"HCHO (\u00d710\u207b\u2074 mol/m\u00b2)",  1e4,      "#01665e"),
    "O3":        (u"O\u2083 (DU)",                             O3_TO_DU, "#762a83"),
    "LST_Day":   (u"LST Day (\u00b0C)",                        1.0,      "#a32020"),
    "LST_Night": (u"LST Night (\u00b0C)",                      1.0,      "#2c5f9e"),
}
MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
          "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR)


def monsoon(ax):
    ax.axvspan(pd.Timestamp(YEAR, 6, 1), pd.Timestamp(YEAR, 9, 30),
               color="#dbe9f6", zorder=0, lw=0)


# ============ 1. daily national means ============
fig, axes = plt.subplots(len(VARIABLES), 1, figsize=(7.5, 11.0), sharex=True)
daily_stats = []

for ax, (var, (label, f, colour)) in zip(axes, VARIABLES.items()):
    path = os.path.join(DATA_DIR, "{}_daily_{}.csv".format(var, YEAR))
    if not os.path.exists(path):
        ax.text(0.5, 0.5, "{} not found".format(os.path.basename(path)),
                transform=ax.transAxes, ha="center")
        print("{}: daily file not found - skipped".format(var))
        continue
    d = pd.read_csv(path)
    d["date"] = pd.to_datetime(d["date"])
    d["value"] = pd.to_numeric(d["mean_value"], errors="coerce") * f
    if "pct_valid" in d.columns:                      # AOD: keep well-covered days only
        d.loc[d["pct_valid"] < COVER_MIN, "value"] = np.nan
    d = d.sort_values("date").set_index("date")

    # remove outlier days (usually cloudy days where only a few pixels were seen)
    win = "{}D".format(ROLL_DAYS)
    med = d["value"].rolling(win, center=True, min_periods=3).median()
    mad = (d["value"] - med).abs().rolling(win, center=True, min_periods=3).median()
    out = (d["value"] - med).abs() > OUTLIER_K * 1.4826 * mad
    n_out = int(out.sum())
    d.loc[out, "value"] = np.nan
    smooth = d["value"].rolling(win, center=True, min_periods=3).median()

    monsoon(ax)
    ax.plot(d.index, d["value"], ".", ms=2.2, color=colour, alpha=0.35)
    ax.plot(smooth.index, smooth, "-", lw=1.3, color=colour)
    ax.set_ylabel(label, fontsize=7)
    ax.tick_params(labelsize=6.5)
    ax.grid(alpha=0.25, lw=0.4)

    n_ok = int(d["value"].notna().sum())
    daily_stats.append({"variable": var, "days_used": n_ok,
                        "mean": d["value"].mean(), "min": d["value"].min(),
                        "max": d["value"].max()})
    print("{:10s} days used: {}   outliers removed: {}".format(var, n_ok, n_out))

axes[-1].xaxis.set_major_locator(mdates.MonthLocator())
axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%b"))
axes[0].set_title("Daily national mean, Bangladesh {}  (dots = daily, line = {}-day median, "
                  "blue band = monsoon)".format(YEAR, ROLL_DAYS), fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "daily_trends_{}.png".format(YEAR)), dpi=400)
plt.close(fig)
pd.DataFrame(daily_stats).round(4).to_csv(
    os.path.join(OUT_DIR, "daily_summary_{}.csv".format(YEAR)), index=False)


# ============ 2. monthly means from the monthly rasters ============
rows = []
for var, (label, f, colour) in VARIABLES.items():
    path = os.path.join(DATA_DIR, "{}_monthly_{}.tif".format(var, YEAR))
    if not os.path.exists(path):
        continue
    with rasterio.open(path) as src:
        bands = []
        for b in range(1, 13):
            a = src.read(b).astype("float64")
            a[(~np.isfinite(a)) | (a <= NODATA + 1)] = np.nan
            bands.append(a * f)
    valid_any = np.any([np.isfinite(a) for a in bands], axis=0)
    n_any = max(int(valid_any.sum()), 1)
    for m, a in enumerate(bands, start=1):
        cover = 100.0 * (np.isfinite(a) & valid_any).sum() / n_any
        rows.append({"variable": var, "month": m, "mean": np.nanmean(a),
                     "std": np.nanstd(a), "coverage_%": cover})

monthly = pd.DataFrame(rows)
monthly.round(4).to_csv(os.path.join(OUT_DIR, "monthly_means_{}.csv".format(YEAR)), index=False)

# ---- LST Day vs Night figure ----
fig, ax = plt.subplots(figsize=(6.5, 3.8))
x = np.arange(1, 13)
for var, name in [("LST_Day", "LST Day"), ("LST_Night", "LST Night")]:
    s = monthly[monthly["variable"] == var].set_index("month").reindex(x)
    colour = VARIABLES[var][2]
    low = s["coverage_%"] < COVER_MIN
    ax.fill_between(x, s["mean"] - s["std"], s["mean"] + s["std"],
                    color=colour, alpha=0.15, lw=0)
    ax.plot(x, s["mean"], "-", color=colour, lw=1.4, label=name)
    ax.plot(x[~low], s["mean"][~low], "o", color=colour, ms=5)
    ax.plot(x[low], s["mean"][low], "o", mfc="white", mec=colour, ms=5)

dn = monthly[monthly["variable"] == "LST_Day"].set_index("month")["mean"] - \
     monthly[monthly["variable"] == "LST_Night"].set_index("month")["mean"]
ax2 = ax.twinx()
ax2.bar(x, dn.reindex(x), width=0.5, color="grey", alpha=0.25, label="Day \u2212 Night")
ax2.set_ylim(0, np.nanmax(dn.values) * 2.5)   # keep bars in the lower part
ax2.set_ylabel(u"Day \u2212 Night difference (\u00b0C)", fontsize=8)
ax2.tick_params(labelsize=7)
ax.set_zorder(ax2.get_zorder() + 1)
ax.patch.set_visible(False)

ax.set_xticks(x)
ax.set_xticklabels(MONTHS, fontsize=7.5)
ax.set_ylabel(u"LST (\u00b0C)", fontsize=8)
ax.tick_params(labelsize=7)
ax.grid(alpha=0.25, lw=0.4)
h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, fontsize=7, loc="upper left")
ax.set_title("Monthly mean LST, Bangladesh {}  (shading = \u00b11 SD, "
             "hollow = coverage < {}%)".format(YEAR, COVER_MIN), fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "lst_day_night_monthly_{}.png".format(YEAR)), dpi=400)
plt.close(fig)

print("\nMonthly means (LST):")
print(monthly[monthly["variable"].str.startswith("LST")]
      .pivot(index="month", columns="variable", values="mean").round(2).to_string())
print("Done.")

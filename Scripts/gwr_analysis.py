import matplotlib
matplotlib.use("Agg")

import os
import time
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from mgwr.gwr import GWR
from mgwr.sel_bw import Sel_BW

# ---------------- settings ----------------
YEAR     = 2025
CELL_DEG = 0.05            # must match correlation_matrix.py
ROOT     = r"C:\Users\User\Desktop\AirPollutionResearchWithRakibVai"
OUT_DIR  = os.path.join(ROOT, "Figures_v2", "Statistics")
CSV      = os.path.join(OUT_DIR, "annual_pixels_{}.csv".format(YEAR))
BOUNDARY = os.path.join(ROOT, "boundary", "BGD_adm0.shp")

TARGETS    = ["LST Day", "LST Night"]            # one GWR model for each
PREDICTORS = ["AOD", u"CH\u2084", "HCHO"]        # UVAI overlaps AOD; O3 is mostly stratospheric
# ------------------------------------------

df = pd.read_csv(CSV, encoding="utf-8")
bnd = gpd.read_file(BOUNDARY)
coords = list(zip(df["lon"], df["lat"]))


def zscore(a):
    return (a - a.mean(axis=0)) / a.std(axis=0)


# standardised, so coefficients are comparable between pollutants
X = zscore(df[PREDICTORS].values)

# grid for turning cell values back into a map
lon0, lat1 = df["lon"].min(), df["lat"].max()
col = np.round((df["lon"] - lon0) / CELL_DEG).astype(int).values
row = np.round((lat1 - df["lat"]) / CELL_DEG).astype(int).values
H, W = row.max() + 1, col.max() + 1
ext = [lon0 - CELL_DEG / 2, lon0 + (W - 0.5) * CELL_DEG,
       lat1 - (H - 0.5) * CELL_DEG, lat1 + CELL_DEG / 2]


def to_map(values):
    m = np.full((H, W), np.nan)
    m[row, col] = values
    return m


summary = []

for target in TARGETS:
    y = zscore(df[target].values).reshape(-1, 1)
    t0 = time.time()
    print("\n{}: searching bandwidth (this can take several minutes)...".format(target))

    bw = Sel_BW(coords, y, X, spherical=True, kernel="bisquare", fixed=False).search()
    res = GWR(coords, y, X, bw, spherical=True, kernel="bisquare", fixed=False).fit()

    # global OLS for comparison
    Xc = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(Xc, y, rcond=None)
    ols_r2 = 1 - ((y - Xc @ beta) ** 2).sum() / ((y - y.mean()) ** 2).sum()

    print("  bandwidth: {} neighbours   time: {:.0f} s".format(int(bw), time.time() - t0))
    print("  OLS R2: {:.3f}   GWR R2: {:.3f}   GWR adj R2: {:.3f}   AICc: {:.1f}".format(
        ols_r2, res.R2, res.adj_R2, res.aicc))
    summary.append({"target": target, "bandwidth": int(bw), "OLS_R2": ols_r2,
                    "GWR_R2": res.R2, "GWR_adjR2": res.adj_R2, "AICc": res.aicc})

    # significance (corrected for multiple testing)
    sig = res.filter_tvals() != 0

    out = pd.DataFrame({"lon": df["lon"], "lat": df["lat"], "local_R2": res.localR2.ravel()})
    for k, p in enumerate(PREDICTORS):
        out["coef_" + p] = res.params[:, k + 1]
        out["sig_" + p] = sig[:, k + 1]
        pct_sig = 100 * sig[:, k + 1].mean()
        pos = 100 * ((res.params[:, k + 1] > 0) & sig[:, k + 1]).mean()
        print("  {:5s} significant in {:.0f}% of cells ({:.0f}% positive)".format(p, pct_sig, pos))
    tag = target.replace(" ", "_")
    out.to_csv(os.path.join(OUT_DIR, "gwr_{}_{}.csv".format(tag, YEAR)), index=False)

    # ---- maps: local R2 + one coefficient map per pollutant ----
    n_pan = 1 + len(PREDICTORS)
    fig, axes = plt.subplots(2, 2, figsize=(8.27, 9.5))
    axes = axes.ravel()

    panels = [("Local R\u00b2", res.localR2.ravel(), None)] + \
             [(u"{} coefficient".format(p), res.params[:, k + 1], sig[:, k + 1])
              for k, p in enumerate(PREDICTORS)]
    lim = np.nanpercentile(np.abs(res.params[:, 1:]), 98)

    for ax, (title, vals, s) in zip(axes, panels):
        if s is None:
            im = ax.imshow(to_map(vals), extent=ext, cmap="viridis", vmin=0, vmax=1)
        else:
            ax.imshow(to_map(np.where(s, np.nan, 1)), extent=ext, cmap="Greys",
                      vmin=0, vmax=4)                     # non-significant cells in light grey
            im = ax.imshow(to_map(np.where(s, vals, np.nan)), extent=ext, cmap="RdBu_r",
                           norm=TwoSlopeNorm(0, -lim, lim))
        bnd.boundary.plot(ax=ax, color="black", lw=0.4)
        ax.set_title(title, fontsize=9, fontweight="bold")
        ax.tick_params(labelsize=6)
        cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
        cb.ax.tick_params(labelsize=6)
    for ax in axes[n_pan:]:
        ax.axis("off")

    fig.suptitle("GWR: {} ~ {}   (grey = not significant)".format(
        target, " + ".join(PREDICTORS)), fontsize=10)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "gwr_{}_{}.png".format(tag, YEAR)), dpi=400)
    plt.close(fig)

pd.DataFrame(summary).round(3).to_csv(os.path.join(OUT_DIR, "gwr_summary_{}.csv".format(YEAR)),
                                      index=False)
print("\nDone.")

import matplotlib
matplotlib.use("Agg")

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

# ---------------- settings ----------------
YEAR    = 2025
ROOT    = r"C:\Users\User\Desktop\AirPollutionResearchWithRakibVai"
OUT_DIR = os.path.join(ROOT, "Figures_v2", "Statistics")
CSV     = os.path.join(OUT_DIR, "annual_pixels_{}.csv".format(YEAR))
# ------------------------------------------

# column in CSV -> (axis label, multiply factor for display)
POLLUTANTS = {
    "AOD":       ("AOD (550 nm)",                       1.0),
    "UVAI":      ("UVAI (index)",                       1.0),
    u"CH\u2084": (u"CH\u2084 (ppb)",                    1.0),
    "HCHO":      (u"HCHO (\u00d710\u207b\u2074 mol/m\u00b2)", 1e4),
    u"O\u2083":  (u"O\u2083 (DU)",                      1.0),
    u"NO\u2082": (u"NO\u2082 (\u00b5mol/m\u00b2)",       1e6),
    "CO":        (u"CO (mol/m\u00b2)",                  1.0),
    u"SO\u2082": (u"SO\u2082 (\u00b5mol/m\u00b2)",       1e6),
}
LST = {"LST Day": "#a32020", "LST Night": "#2c5f9e"}

df = pd.read_csv(CSV, encoding="utf-8")
rows = []

# each pollutant takes two panels (LST Day, LST Night);
# up to 5 pollutants: one pair per row, above that: two pairs per row
n_pol = len(POLLUTANTS)
NPAIR = 1 if n_pol <= 5 else 2
NROW = int(np.ceil(n_pol / float(NPAIR)))
fig, axes = plt.subplots(NROW, 2 * NPAIR, figsize=(7.5 if NPAIR == 1 else 11.5, 2.3 * NROW),
                         squeeze=False)

for i, (pcol, (plabel, f)) in enumerate(POLLUTANTS.items()):
    x = df[pcol].values * f
    for j, (lcol, colour) in enumerate(LST.items()):
        y = df[lcol].values
        res = stats.linregress(x, y)
        r2 = res.rvalue ** 2
        rows.append({"LST": lcol, "pollutant": pcol, "n": len(x),
                     "slope": res.slope, "intercept": res.intercept,
                     "r": res.rvalue, "R2": r2, "p": res.pvalue})

        ax = axes[i // NPAIR, 2 * (i % NPAIR) + j]
        ax.scatter(x, y, s=2, alpha=0.25, color=colour, edgecolors="none", rasterized=True)
        xs = np.linspace(np.percentile(x, 0.5), np.percentile(x, 99.5), 50)
        ax.plot(xs, res.intercept + res.slope * xs, color="black", lw=1.2)

        sign = "+" if res.intercept >= 0 else "\u2212"
        ax.text(0.03, 0.96,
                "y = {:.3g}x {} {:.3g}\nR\u00b2 = {:.2f}   r = {:.2f}".format(
                    res.slope, sign, abs(res.intercept), r2, res.rvalue),
                transform=ax.transAxes, ha="left", va="top", fontsize=6.5,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1.5))

        ax.set_xlim(np.percentile(x, 0.5), np.percentile(x, 99.5))
        ax.set_xlabel(plabel, fontsize=7.5)
        ax.set_ylabel(u"{} (\u00b0C)".format(lcol), fontsize=7.5)
        ax.tick_params(labelsize=6.5)
        if i < NPAIR:
            ax.set_title(lcol, fontsize=9, fontweight="bold")

for k in range(2 * n_pol, axes.size):             # empty panels
    axes.ravel()[k].axis("off")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "regression_LST_pollutants_{}.png".format(YEAR)), dpi=400)
plt.close(fig)

table = pd.DataFrame(rows)
table.to_csv(os.path.join(OUT_DIR, "regression_LST_pollutants_{}.csv".format(YEAR)), index=False)
print(table[["LST", "pollutant", "slope", "r", "R2"]].round(3).to_string(index=False))
print("Done.")

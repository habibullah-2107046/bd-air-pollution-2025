# -*- coding: utf-8 -*-
"""
review_03_seasonal.py  -  seasonal analysis redone for the V1 review
  - winter = ONE continuous winter, December 2024 - February 2025
  - a month enters a season only if its coverage >= threshold (50 % and 75 %)
  - seasonal coverage = mean of the monthly coverages of the months used
  - correlations on a COMMON set of cells (cells valid in every variable that
    has data in that season), with n reported, and the original pairwise way
    for comparison; p values corrected for spatial autocorrelation
  - Aqua LST version
Needs: review_01_data.py run first
Output: Review_V1/03_seasonal/
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import review_common as rc

THRESHOLDS = [50, 75]
MIN_CELLS = 30
OUT = rc.ensure_out("03_seasonal")

grid = rc.Grid()
cube = rc.load_cube(grid)
cov = {n: rc.coverage(cube, n) for n in rc.variables(cube)}


def season_value(n, idx, thr):
    """Cell mean over the season months that pass the coverage threshold."""
    use = [k for k in idx if cov[n][k] >= thr]
    if not use:
        return None, use
    st = cube[n][use]
    cnt = np.isfinite(st).sum(0)
    v = np.where(cnt > 0, np.nansum(st, 0) / np.maximum(cnt, 1), np.nan)
    v[~grid.good] = np.nan
    return v, use


VERSIONS = [("terra", {}), ("aqua", {"LST_Day": "LSTAqua_Day", "LST_Night": "LSTAqua_Night"})]
mean_rows, corr_rows = [], []
for thr in THRESHOLDS:
    for ver, sub in VERSIONS:
        if any(a not in cube for a in sub.values()):
            continue
        for sname, span, idx in rc.SEASONS:
            sv, used = {}, {}
            for n in rc.POLLUTANTS + rc.LSTS:
                nn = sub.get(n, n)
                v, u = season_value(nn, idx, thr)
                used[n] = u
                if v is not None:
                    sv[n] = v
            if ver == "terra":
                for n in rc.POLLUTANTS + rc.LSTS:
                    u = used[n]
                    row = {"threshold_%": thr, "season": sname, "variable": rc.label(n), "unit": rc.unit(n),
                           "months_used": ", ".join(rc.PERIOD_LABEL[k] for k in u),
                           "coverage_%": np.mean([cov[n][k] for k in idx]) if idx else np.nan,
                           "coverage_used_months_%": np.mean([cov[n][k] for k in u]) if u else np.nan}
                    if n in sv:
                        a = sv[n][grid.good]
                        row.update({"mean_all_cells": np.nanmean(a), "sd_all_cells": np.nanstd(a)})
                    mean_rows.append(row)
            have = [n for n in rc.POLLUTANTS + rc.LSTS if n in sv]
            common = np.all([np.isfinite(sv[n]) for n in have], axis=0) & grid.good
            n_common = int(common.sum())
            if ver == "terra":
                for r in mean_rows:
                    if r["threshold_%"] == thr and r["season"] == sname and \
                            r["variable"] in [rc.label(n) for n in have]:
                        n0 = [n for n in have if rc.label(n) == r["variable"]][0]
                        r["mean_common_cells"] = np.nanmean(sv[n0][common]) if n_common else np.nan
                        r["n_common_cells"] = n_common
            sc = rc.SpatialCorr(grid.lon[common], grid.lat[common]) if n_common >= MIN_CELLS else None
            for lst in rc.LSTS:
                if lst not in sv:
                    continue
                for pol in rc.POLLUTANTS:
                    row = {"threshold_%": thr, "LST_sensor": ver, "season": sname, "LST": lst,
                           "pollutant": pol, "variables_in_common_set": " ".join(have)}
                    if pol not in sv:
                        row["note"] = "no month passes the threshold"
                        corr_rows.append(row)
                        continue
                    pw = np.isfinite(sv[pol]) & np.isfinite(sv[lst]) & grid.good
                    row["n_pairwise"] = int(pw.sum())
                    if pw.sum() >= MIN_CELLS:
                        row["r_pairwise"] = np.corrcoef(sv[pol][pw], sv[lst][pw])[0, 1]
                    row["n_common"] = n_common
                    if sc is not None:
                        r, p, ne, pc = sc.test(sv[pol][common], sv[lst][common], "pearson",
                                               kx=pol + sname, ky=lst + sname)
                        rho, _, ne_s, pc_s = sc.test(sv[pol][common], sv[lst][common], "spearman",
                                                     kx=pol + sname, ky=lst + sname)
                        row.update({"r_common": r, "p_naive": p, "n_eff": ne, "p_corrected": pc,
                                    "rho_common": rho, "p_rho_corrected": pc_s})
                    corr_rows.append(row)
            print("thr {} {} {}: common cells {} ({} variables)".format(thr, ver, sname, n_common, len(have)))

means = pd.DataFrame(mean_rows)
means.to_csv(os.path.join(OUT, "03_seasonal_means.csv"), index=False)
corr = pd.DataFrame(corr_rows)
corr.to_csv(os.path.join(OUT, "03_seasonal_correlations.csv"), index=False)

order = [s[0] for s in rc.SEASONS]
m50 = means[means["threshold_%"] == 50]
print("\nSeasonal national means, all cells (threshold 50 %)")
print(m50.pivot(index="variable", columns="season", values="mean_all_cells")[order].round(3).to_string())
print("\nSeasonal coverage (mean of monthly coverage, %)")
print(m50.pivot(index="variable", columns="season", values="coverage_%")[order].round(0).to_string())
print("\nMonths used (threshold 50 %)")
print(m50.pivot(index="variable", columns="season", values="months_used")[order].to_string())

for thr in THRESHOLDS:
    for ver in ["terra", "aqua"]:
        c = corr[(corr["threshold_%"] == thr) & (corr["LST_sensor"] == ver)]
        if c.empty or "r_common" not in c:
            continue
        print("\nLST-pollutant r on common cells, threshold {} %, {} LST (corrected p in brackets)".format(thr, ver))
        c = c.assign(txt=[("{:+.2f} ({:.2g})".format(r, p) if np.isfinite(r) else "n/a")
                          for r, p in zip(c["r_common"].fillna(np.nan), c["p_corrected"].fillna(np.nan))])
        print(c.pivot_table(index=["LST", "pollutant"], columns="season", values="txt",
                            aggfunc="first").reindex(columns=order).to_string())
        print("common cells per season:",
              c.groupby("season")["n_common"].first().reindex(order).to_dict())

# heatmap of the main version (threshold 50 %, Terra, common cells)
c = corr[(corr["threshold_%"] == 50) & (corr["LST_sensor"] == "terra")]
if "r_common" in c:
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.6))
    for ax, lst in zip(axes, rc.LSTS):
        t = c[c["LST"] == lst]
        R = t.pivot(index="pollutant", columns="season", values="r_common").reindex(
            index=rc.POLLUTANTS, columns=order)
        P = t.pivot(index="pollutant", columns="season", values="p_corrected").reindex(
            index=rc.POLLUTANTS, columns=order)
        im = ax.imshow(np.ma.masked_invalid(R.values.astype(float)), cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
        for i in range(R.shape[0]):
            for j in range(R.shape[1]):
                r, p = R.values[i, j], P.values[i, j]
                if not np.isfinite(r):
                    ax.text(j, i, "n/a", ha="center", va="center", fontsize=7, color="grey")
                    continue
                star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
                ax.text(j, i, "{:.2f}{}".format(r, star), ha="center", va="center", fontsize=7,
                        color="white" if abs(r) > 0.6 else "black")
        ns = t.groupby("season")["n_common"].first().reindex(order)
        ax.set_xticks(range(4))
        ax.set_xticklabels(["{}\n(n = {})".format(s, int(n) if np.isfinite(n) else 0) for s, n in ns.items()], fontsize=7)
        ax.set_yticks(range(len(rc.POLLUTANTS)))
        ax.set_yticklabels([rc.label(n) for n in rc.POLLUTANTS], fontsize=8)
        ax.set_title(rc.label(lst), fontsize=9)
    fig.colorbar(im, ax=list(axes), fraction=0.03, label="Pearson r")
    fig.suptitle("Seasonal LST-pollutant correlation on a common set of cells; winter = Dec 2024-Feb 2025\n"
                 "p corrected for spatial autocorrelation: * < 0.05, ** < 0.01, *** < 0.001", fontsize=8)
    plt.savefig(os.path.join(OUT, "03_seasonal_correlation_common_cells.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

print("\nDone. Outputs in", OUT)

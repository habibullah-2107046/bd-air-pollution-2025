# -*- coding: utf-8 -*-
"""
review_01_data.py  -  data checks for the V1 review (run this first)
  - builds the 0.05 degree monthly cube (Dec 2024 + Jan-Dec 2025) used by
    review_02 / 03 / 04, in one unit convention
  - monthly coverage (% of land cells with data) and national monthly means
  - Terra and Aqua overpass (view) time            (Terra orbit drift)
  - Terra vs Aqua LST, monthly                      (Terra orbit drift)
  - day-night LST difference on pixels with BOTH    (day-night difference)
  - AOD: original vs strict QA, daily vs monthly    (June AOD)
  - SO2: original vs cloud fraction < 0.3           (SO2 swing)
  - CH4: original vs bias-corrected                 (albedo artefact)
Needs: BD_AQ_2025/ and BD_AQ_review/ (exports of review_exports_GEE.js)
Output: Review_V1/01_data/  and the cache Review_V1/cube.npz
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import review_common as rc

OUT = rc.ensure_out("01_data")
grid = rc.Grid()
print("Gridding (first run takes a few minutes)...")
cube = rc.load_cube(grid)

# ---------------- coverage and national means ----------------
rows = []
for n in rc.variables(cube):
    c = cube[n]
    cov = rc.coverage(cube, n)
    for k, per in enumerate(rc.PERIODS):
        b = c[k][grid.good]
        rows.append({"variable": n, "label": rc.label(n), "unit": rc.unit(n), "period": per,
                     "period_label": rc.PERIOD_LABEL[k], "coverage_%": cov[k],
                     "mean": np.nanmean(b) if np.isfinite(b).any() else np.nan,
                     "std": np.nanstd(b) if np.isfinite(b).any() else np.nan})
tab = pd.DataFrame(rows)
tab.to_csv(os.path.join(OUT, "01_monthly_coverage_and_means.csv"), index=False)
cov_w = tab.pivot(index="period_label", columns="variable", values="coverage_%").reindex(rc.PERIOD_LABEL)
mean_w = tab.pivot(index="period_label", columns="variable", values="mean").reindex(rc.PERIOD_LABEL)
print("\nCoverage (% of land pixels with data, native resolution)")
print(cov_w.round(0).to_string())
print("\nNational monthly means (unified units)")
print(mean_w.round(3).to_string())

# ---------------- overpass time ----------------
ov_path = os.path.join(rc.REV_DIR, "overpass_time.csv")
if os.path.exists(ov_path):
    ov = pd.read_csv(ov_path)
    keep = [c for c in ov.columns if c.startswith(("Day_view", "Night_view")) or c in ("sensor", "period")]
    ov = ov[keep]
    ov.to_csv(os.path.join(OUT, "01_overpass_time.csv"), index=False)
    print("\nOverpass (view) time, local solar time in hours, national monthly mean")
    print(ov.round(2).to_string(index=False))
    fig, ax = plt.subplots(figsize=(7, 3.6))
    for s, col in [("Terra", "#2c5f9e"), ("Aqua", "#bf812d")]:
        t = ov[ov["sensor"] == s].set_index("period").reindex(rc.PERIODS)
        ax.plot(range(13), t["Day_view_time_mean"], "o-", color=col, label=s + " day")
        ax.plot(range(13), t["Night_view_time_mean"] - 12, "s--", color=col, label=s + " night (minus 12 h)")
    ax.set_xticks(range(13))
    ax.set_xticklabels(rc.PERIOD_LABEL, rotation=60, fontsize=7)
    ax.set_ylabel("Local solar time (h)")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "01_overpass_time.png"), dpi=300)
    plt.close(fig)
else:
    print("\noverpass_time.csv not found in BD_AQ_review - skipped")

# ---------------- Terra vs Aqua ----------------
comp = []
for t_name, a_name in [("LST_Day", "LSTAqua_Day"), ("LST_Night", "LSTAqua_Night")]:
    if a_name not in cube:
        continue
    for k, per in enumerate(rc.PERIODS):
        t, a = cube[t_name][k], cube[a_name][k]
        both = np.isfinite(t) & np.isfinite(a)
        if both.sum() > 30:
            r = np.corrcoef(t[both], a[both])[0, 1]
            comp.append({"pair": t_name + " vs " + a_name, "period": rc.PERIOD_LABEL[k],
                         "n_cells_both": int(both.sum()), "terra_mean": t[both].mean(),
                         "aqua_mean": a[both].mean(), "aqua_minus_terra": (a[both] - t[both]).mean(),
                         "spatial_r": r})
if comp:
    comp = pd.DataFrame(comp)
    comp.to_csv(os.path.join(OUT, "01_terra_vs_aqua_lst.csv"), index=False)
    print("\nTerra vs Aqua LST on cells with both (monthly)")
    print(comp.round(2).to_string(index=False))

# ---------------- day-night difference on pixels with both (native 1 km) ----------------
dn = []
dayM, trD, _, _, _ = rc.read_all(os.path.join(rc.DATA_DIR, "LST_Day_monthly_2025.tif"))
nigM, trN, _, _, _ = rc.read_all(os.path.join(rc.DATA_DIR, "LST_Night_monthly_2025.tif"))
landD = grid.land_on(dayM[0].shape, trD) > 0
for m in range(12):
    d, n_ = dayM[m].copy(), nigM[m].copy()
    d[~landD] = np.nan
    n_[~landD] = np.nan
    both = np.isfinite(d) & np.isfinite(n_)
    dn.append({"month": m + 1,
               "day_mean_all": np.nanmean(d), "night_mean_all": np.nanmean(n_),
               "diff_of_means_all": np.nanmean(d) - np.nanmean(n_),
               "day_cover_%": 100.0 * np.isfinite(d).sum() / landD.sum(),
               "night_cover_%": 100.0 * np.isfinite(n_).sum() / landD.sum(),
               "both_cover_%": 100.0 * both.sum() / landD.sum(),
               "diff_common_pixels": (d[both] - n_[both]).mean() if both.any() else np.nan})
dn = pd.DataFrame(dn)
dn.to_csv(os.path.join(OUT, "01_day_night_difference.csv"), index=False)
print("\nDay-night LST difference: all pixels vs pixels with both (1 km, 2025)")
print(dn.round(2).to_string(index=False))

# ---------------- AOD: strict QA and daily vs monthly ----------------
if "AOD_strict" in cube:
    a1 = tab[tab["variable"] == "AOD"].set_index("period_label")
    a2 = tab[tab["variable"] == "AOD_strict"].set_index("period_label")
    aod = pd.DataFrame({"AOD_mean": a1["mean"], "AOD_cover_%": a1["coverage_%"],
                        "strictQA_mean": a2["mean"], "strictQA_cover_%": a2["coverage_%"]}).reindex(rc.PERIOD_LABEL)
    # same cells only
    same = []
    for k in range(13):
        o, s = cube["AOD"][k], cube["AOD_strict"][k]
        b = np.isfinite(o) & np.isfinite(s)
        same.append((o[b].mean() if b.any() else np.nan, s[b].mean() if b.any() else np.nan))
    aod["AOD_mean_same_cells"] = [x[0] for x in same]
    aod["strictQA_mean_same_cells"] = [x[1] for x in same]
else:
    aod = tab[tab["variable"] == "AOD"].set_index("period_label")[["mean", "coverage_%"]]
daily = pd.read_csv(os.path.join(rc.DATA_DIR, "AOD_daily_2025.csv"), parse_dates=["date"])
daily["month"] = daily["date"].dt.month
dm = daily.groupby("month").agg(
    daily_mean_all_days=("mean_value", "mean"),
    daily_mean_days_50pct=("mean_value", lambda s: s[daily.loc[s.index, "pct_valid"] >= 50].mean()),
    n_days_50pct=("pct_valid", lambda s: int((s >= 50).sum())),
    mean_pct_valid=("pct_valid", "mean"))
dm.index = ["{} 2025".format(m) for m in
            ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]]
aod = aod.join(dm)
aod.to_csv(os.path.join(OUT, "01_aod_checks.csv"))
print("\nAOD: monthly composite (original / strict QA) vs daily national means")
print(aod.round(3).to_string())

# ---------------- SO2 and CH4 versions ----------------
for base, alt in [("SO2", "SO2_cf03"), ("CH4", "CH4_bc")]:
    if alt not in cube:
        continue
    t = pd.DataFrame({
        base + "_mean": tab[tab["variable"] == base].set_index("period_label")["mean"],
        base + "_cover_%": tab[tab["variable"] == base].set_index("period_label")["coverage_%"],
        alt + "_mean": tab[tab["variable"] == alt].set_index("period_label")["mean"],
        alt + "_cover_%": tab[tab["variable"] == alt].set_index("period_label")["coverage_%"],
    }).reindex(rc.PERIOD_LABEL)
    r = []
    for k in range(13):
        x, y = cube[base][k], cube[alt][k]
        b = np.isfinite(x) & np.isfinite(y)
        r.append(np.corrcoef(x[b], y[b])[0, 1] if b.sum() > 30 else np.nan)
    t["spatial_r"] = r
    t.to_csv(os.path.join(OUT, "01_{}_versions.csv".format(base)))
    print("\n{}: original vs {} (national monthly means)".format(base, alt))
    print(t.round(3).to_string())

print("\nDone. Outputs in", OUT)

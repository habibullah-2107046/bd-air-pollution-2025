# -*- coding: utf-8 -*-
"""
review_06_station.py  -  station checks for the V1 review
  (a) PM2.5 vs AOD and PM2.5 vs UVAI on the SAME station-months, plus the
      partial correlation controlling for month (shared seasonal cycle)
  (b) station ranking: only stations with >= 9 valid months, and a
      seasonally adjusted mean
  (c) station-months where PM2.5 > PM10 (physically inconsistent), and
      PM2.5 / PM10 means on common months only
Needs: BD_AQ_2025/doe_cams_monthly_2025.csv and
       Station_validation/station_validation_pairs.csv (from station_validation.py)
Output: Review_V1/06_station/
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
from review_common import ROOT, DATA_DIR, ensure_out

MIN_CAPTURE = 50       # same rule as station_validation.py
MIN_MONTHS = 9         # for the ranking
OUT = ensure_out("06_station")

d = pd.read_csv(os.path.join(DATA_DIR, "doe_cams_monthly_2025.csv"))
d["mean"] = pd.to_numeric(d["mean"], errors="coerce")
d = d[(d["capture_pct"] >= MIN_CAPTURE) & d["mean"].notna()].copy()
pairs = pd.read_csv(os.path.join(ROOT, "Station_validation", "station_validation_pairs.csv"))


def partial_r(x, y, month):
    """Correlation of x and y after removing the month means (month dummies)."""
    df = pd.DataFrame({"x": x, "y": y, "m": month})
    rx = df["x"] - df.groupby("m")["x"].transform("mean")
    ry = df["y"] - df.groupby("m")["y"].transform("mean")
    n_m = df["m"].nunique()
    r = np.corrcoef(rx, ry)[0, 1]
    dof = len(df) - 2 - (n_m - 1)
    t = r * np.sqrt(dof / max(1e-12, 1 - r * r))
    return r, 2 * stats.t.sf(abs(t), dof), dof


# ======================= (a) same station-months =======================
rows = []
for gp in ["PM2.5", "PM10"]:
    a = pairs[(pairs["ground_parameter"] == gp) & (pairs["sat_variable"] == "AOD")]
    u = pairs[(pairs["ground_parameter"] == gp) & (pairs["sat_variable"] == "UVAI")]
    for nm, t in [("AOD", a), ("UVAI", u)]:
        if len(t) > 3:
            r, p = stats.pearsonr(t["ground_value"], t["sat_value"])
            rows.append({"ground": gp, "satellite": nm, "sample": "all available station-months",
                         "n": len(t), "r": r, "p": p})
    j = a.merge(u, on=["station", "month"], suffixes=("_aod", "_uvai"))
    if len(j) > 3:
        for nm in ["aod", "uvai"]:
            r, p = stats.pearsonr(j["ground_value_aod"], j["sat_value_" + nm])
            rows.append({"ground": gp, "satellite": nm.upper(), "sample": "common station-months",
                         "n": len(j), "r": r, "p": p})
            pr, pp, dof = partial_r(j["ground_value_aod"], j["sat_value_" + nm], j["month"])
            rows.append({"ground": gp, "satellite": nm.upper(),
                         "sample": "common station-months, partial r (month removed)",
                         "n": len(j), "r": pr, "p": pp})
        # correlation of the monthly mean cycles (12 points at most): how much is seasonal
        mm = j.groupby("month")[["ground_value_aod", "sat_value_aod", "sat_value_uvai"]].mean()
        for nm in ["aod", "uvai"]:
            if len(mm) > 3:
                r, p = stats.pearsonr(mm["ground_value_aod"], mm["sat_value_" + nm])
                rows.append({"ground": gp, "satellite": nm.upper(),
                             "sample": "monthly network means (seasonal cycle only)",
                             "n": len(mm), "r": r, "p": p})
        # months present for each satellite variable
        print("{}: AOD station-months {}, UVAI station-months {}, common {}".format(
            gp, len(a), len(u), len(j)))
        miss = u.merge(a[["station", "month"]], on=["station", "month"], how="left", indicator=True)
        miss = miss[miss["_merge"] == "left_only"]["month"].value_counts().sort_index()
        print("   station-months with UVAI but no AOD, by month:", miss.to_dict())

res_a = pd.DataFrame(rows)
res_a.to_csv(os.path.join(OUT, "06a_pm_aod_uvai_same_months.csv"), index=False)
print("\n(a) PM vs AOD / UVAI")
print(res_a.round(3).to_string(index=False))

# ======================= (b) ranking =======================
rank_rows = []
for gp in ["PM2.5", "PM10"]:
    g = d[d["parameter"] == gp][["station", "month", "mean"]]
    month_mean = g.groupby("month")["mean"].mean()
    grand = g["mean"].mean()
    g = g.assign(anom=g["mean"] - g["month"].map(month_mean))
    s = g.groupby("station").agg(n_months=("mean", "size"), mean_raw=("mean", "mean"),
                                 months=("month", lambda m: ",".join(str(int(x)) for x in sorted(m))),
                                 anom=("anom", "mean"))
    s["mean_season_adjusted"] = grand + s["anom"]
    s["rank_raw_all"] = s["mean_raw"].rank(ascending=False).astype(int)
    ok = s["n_months"] >= MIN_MONTHS
    s["rank_raw_9plus"] = np.nan
    s.loc[ok, "rank_raw_9plus"] = s.loc[ok, "mean_raw"].rank(ascending=False)
    s["rank_adjusted"] = s["mean_season_adjusted"].rank(ascending=False).astype(int)
    s.insert(0, "parameter", gp)
    rank_rows.append(s.drop(columns="anom").sort_values("mean_season_adjusted", ascending=False))
rank = pd.concat(rank_rows)
rank.to_csv(os.path.join(OUT, "06b_station_ranking.csv"))
print("\n(b) station ranking (raw mean, >= {} months, seasonally adjusted)".format(MIN_MONTHS))
print(rank.round(1).to_string())

# ======================= (c) PM2.5 > PM10 =======================
pm = d[d["parameter"].isin(["PM2.5", "PM10"])].pivot_table(
    index=["station", "month"], columns="parameter", values="mean")
both = pm.dropna()
both = both.assign(ratio=both["PM2.5"] / both["PM10"])
bad = both[both["ratio"] > 1].reset_index()
bad.to_csv(os.path.join(OUT, "06c_pm25_gt_pm10_station_months.csv"), index=False)
print("\n(c) station-months with PM2.5 > PM10 (both with >= {}% capture): {} of {}".format(
    MIN_CAPTURE, len(bad), len(both)))
print(bad.round(2).to_string(index=False))

only = pm.reset_index()
summ = []
for st, g in only.groupby("station"):
    c = g.dropna(subset=["PM2.5", "PM10"])
    summ.append({"station": st,
                 "PM2.5_months": int(g["PM2.5"].notna().sum()), "PM2.5_mean_all": g["PM2.5"].mean(),
                 "PM10_months": int(g["PM10"].notna().sum()), "PM10_mean_all": g["PM10"].mean(),
                 "common_months": len(c), "PM2.5_mean_common": c["PM2.5"].mean(),
                 "PM10_mean_common": c["PM10"].mean(),
                 "ratio_mean_common": (c["PM2.5"] / c["PM10"]).mean() if len(c) else np.nan,
                 "months_ratio_gt_1": int((c["PM2.5"] > c["PM10"]).sum())})
summ = pd.DataFrame(summ)
summ.to_csv(os.path.join(OUT, "06c_pm25_pm10_common_months.csv"), index=False)
print("\nPM2.5 and PM10 on all months vs common months")
print(summ.round(2).to_string(index=False))

# months reported per station (missing May and low data capture)
avail = d.pivot_table(index="station", columns="parameter", values="month", aggfunc="count")
avail.to_csv(os.path.join(OUT, "06d_months_per_station.csv"))
print("\nValid months per station and parameter (capture >= {}%)".format(MIN_CAPTURE))
print(avail.to_string())
print("\nDone. Outputs in", OUT)

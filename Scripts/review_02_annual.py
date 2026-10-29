# -*- coding: utf-8 -*-
"""
review_02_annual.py  -  annual analysis redone for the V1 review
  - ONE definition of the annual value for every variable: the mean of the
    monthly values over the months that all variables share (months with
    coverage >= threshold, and the cell valid in all variables in that month).
    This removes the different time weighting of the original annual composites.
  - coverage threshold 50 % and 75 %
  - correlations with p values corrected for spatial autocorrelation
    (effective sample size, Clifford et al. 1989 / Dutilleul 1993)
  - sensitivity versions: Aqua LST, bias-corrected CH4, strict-QA AOD,
    cells near the straight AOD line removed
  - CH4 vs SWIR surface albedo (albedo artefact check)
  - PCA of the pollutants only, then LST regressed on the PC scores
  - writes the annual cell tables (with control variables) for review_04 / 05
Needs: review_01_data.py run first (it builds Review_V1/cube.npz)
Output: Review_V1/02_annual/
"""
import os
import numpy as np
import pandas as pd
from scipy import stats, ndimage
import review_common as rc

THRESHOLDS = [50, 75]
MIN_COMMON_MONTHS = 6       # a cell needs at least this many shared months
AOD_LINE_BUFFER_KM = 10     # cells this close to the AOD line are removed in one test
MAIN = rc.POLLUTANTS + rc.LSTS
OUT = rc.ensure_out("02_annual")

grid = rc.Grid()
cube = rc.load_cube(grid)
cov = {n: rc.coverage(cube, n) for n in rc.variables(cube)}


# ---------------- control variables on the 0.05 degree grid ----------------
def load_controls():
    path = os.path.join(rc.REV_DIR, "controls_1km.tif")
    if not os.path.exists(path):
        print("controls_1km.tif not found - control variables skipped")
        return {}
    a, tr, crs, b, names = rc.read_all(path)
    out = {}
    for i, nm in enumerate(names):
        out[nm] = grid.grid_band(a[i], tr)
    # distance to permanent water (WorldCover class 80 fraction >= 0.5), km
    if "water_wc" in names:
        w = a[names.index("water_wc")]
        water = np.nan_to_num(w) >= 0.5
        lat_c = (b.top + b.bottom) / 2.0
        dy = abs(tr.e) * 111.32
        dx = tr.a * 111.32 * np.cos(np.radians(lat_c))
        dist = ndimage.distance_transform_edt(~water, sampling=(dy, dx))
        dist[~np.isfinite(w)] = np.nan
        out["dist_water_km"] = grid.grid_band(dist, tr)
    return out


controls = load_controls()


# ---------------- annual value on shared months ----------------
def annual(names, thr):
    usable = [k for k in rc.IDX_2025 if all(cov[n][k] >= thr for n in names)]
    if not usable:
        return None, usable, None
    st = np.stack([cube[n][usable] for n in names])            # (V, M, H, W)
    allfin = np.all(np.isfinite(st), axis=0)                    # (M, H, W)
    nmon = allfin.sum(axis=0)
    ok = (nmon >= MIN_COMMON_MONTHS) & grid.good
    vals = {}
    for i, n in enumerate(names):
        s = np.where(allfin, st[i], 0.0).sum(axis=0) / np.maximum(nmon, 1)
        s[~ok] = np.nan
        vals[n] = s
    return vals, usable, np.where(ok, nmon, 0)


def sin_x(lon, lat):
    return 6371007.181 * np.cos(np.radians(lat)) * np.radians(lon)


x_line = 9173740.0       # from review_07_maps.py, read from its output if present
p_line = os.path.join(rc.OUT_DIR, "07_maps", "07a_aod_line.txt")
if os.path.exists(p_line):
    for ln in open(p_line):
        if ln.startswith("x_line_m"):
            x_line = float(ln.split("=")[1])

VARIANTS = [
    # name, threshold, substitutions {main variable: alternative}, extra cell filter
    ("main_thr50", 50, {}, None),
    ("main_thr75", 75, {}, None),
    ("aqua_thr50", 50, {"LST_Day": "LSTAqua_Day", "LST_Night": "LSTAqua_Night"}, None),
    ("ch4bc_thr50", 50, {"CH4": "CH4_bc"}, None),
    ("aodstrict_thr50", 50, {"AOD": "AOD_strict"}, None),
    ("no_aod_line_thr50", 50, {}, "aod_line"),
]

corr_rows, month_rows, desc_rows = [], [], []
tables = {}
for vname, thr, sub, filt in VARIANTS:
    names = [sub.get(n, n) for n in MAIN]
    if any(n not in cube for n in names):
        print("{}: variable missing - skipped".format(vname))
        continue
    vals, usable, nmon = annual(names, thr)
    if vals is None:
        print("{}: no month passes the threshold - skipped".format(vname))
        continue
    ok = np.isfinite(vals[names[0]])
    if filt == "aod_line":
        ok &= np.abs(sin_x(grid.lon, grid.lat) - x_line) > AOD_LINE_BUFFER_KM * 1000.0
    df = pd.DataFrame({"lon": grid.lon[ok], "lat": grid.lat[ok], "n_months": nmon[ok]})
    for main_n, n in zip(MAIN, names):
        df[main_n] = vals[n][ok]
    for c, g in controls.items():
        df[c] = g[ok]
    tables[vname] = df
    df.to_csv(os.path.join(OUT, "annual_cells_{}.csv".format(vname)), index=False)
    used = [rc.PERIOD_LABEL[k].split()[0] for k in usable]
    month_rows.append({"variant": vname, "threshold_%": thr, "months_shared": ", ".join(used),
                       "n_months_shared": len(usable), "n_cells": len(df),
                       "median_months_per_cell": float(np.median(df["n_months"]))})
    print("\n=== {} === months: {}  cells: {}".format(vname, ", ".join(used), len(df)))

    sc = rc.SpatialCorr(df["lon"].values, df["lat"].values)
    for lst in rc.LSTS:
        for pol in rc.POLLUTANTS:
            row = {"variant": vname, "LST": lst, "pollutant": pol, "n": len(df)}
            for meth, tag in [("pearson", "r"), ("spearman", "rho")]:
                r, p, ne, pc = sc.test(df[pol].values, df[lst].values, meth, kx=pol, ky=lst)
                row.update({tag: r, "p_" + tag: p, "n_eff_" + tag: ne, "p_" + tag + "_corrected": pc})
            sl, ic, _, _, _ = stats.linregress(df[pol].values, df[lst].values)
            row.update({"slope": sl, "intercept": ic, "R2": row["r"] ** 2})
            corr_rows.append(row)
    if vname == "main_thr50":
        # full correlation matrix with corrected p (revised Figure 11)
        k = len(MAIN)
        R = np.zeros((k, k)); P = np.zeros((k, k)); NE = np.zeros((k, k))
        for i in range(k):
            for j in range(k):
                if i == j:
                    R[i, j], P[i, j], NE[i, j] = 1, 0, len(df)
                    continue
                R[i, j], _, NE[i, j], P[i, j] = sc.test(df[MAIN[i]].values, df[MAIN[j]].values,
                                                       "pearson", kx=MAIN[i], ky=MAIN[j])
        labs = [rc.label(n) for n in MAIN]
        pd.DataFrame(R, index=labs, columns=labs).round(3).to_csv(os.path.join(OUT, "02_corr_matrix_r.csv"))
        pd.DataFrame(P, index=labs, columns=labs).to_csv(os.path.join(OUT, "02_corr_matrix_p_corrected.csv"))
        pd.DataFrame(NE, index=labs, columns=labs).round(0).to_csv(os.path.join(OUT, "02_corr_matrix_n_eff.csv"))
        # descriptive statistics (revised Table 3)
        for n in MAIN:
            desc_rows.append({"variable": rc.label(n), "unit": rc.unit(n), "mean": df[n].mean(),
                              "sd": df[n].std(), "min": df[n].min(), "max": df[n].max(),
                              "n_cells": len(df)})

pd.DataFrame(month_rows).to_csv(os.path.join(OUT, "02_months_used.csv"), index=False)
corr = pd.DataFrame(corr_rows)
corr.to_csv(os.path.join(OUT, "02_lst_pollutant_correlations_all_variants.csv"), index=False)
pd.DataFrame(desc_rows).to_csv(os.path.join(OUT, "02_descriptive_main_thr50.csv"), index=False)
print("\nMonths shared by all variables")
print(pd.DataFrame(month_rows).to_string(index=False))

# ---------------- compare with the original annual correlations ----------------
old_p = os.path.join(rc.ROOT, "Figures_v2", "Statistics", "corr_pearson_r_2025.csv")
if os.path.exists(old_p):
    old = pd.read_csv(old_p, index_col=0, encoding="utf-8")
    cmp_rows = []
    for _, r in corr[corr["variant"] == "main_thr50"].iterrows():
        try:
            o = old.loc[rc.label(r["LST"]), rc.label(r["pollutant"])]
        except KeyError:
            o = np.nan
        cmp_rows.append({"LST": r["LST"], "pollutant": r["pollutant"], "r_original": o,
                         "r_new_thr50": r["r"], "p_new_corrected": r["p_r_corrected"],
                         "n_eff": r["n_eff_r"]})
    cmp = pd.DataFrame(cmp_rows)
    cmp.to_csv(os.path.join(OUT, "02_compare_original_vs_new.csv"), index=False)
    print("\nOriginal annual composite vs new shared-month annual (Pearson r)")
    print(cmp.round(3).to_string(index=False))

wide = corr.pivot_table(index=["LST", "pollutant"], columns="variant", values="r")
wide.to_csv(os.path.join(OUT, "02_r_by_variant.csv"))
print("\nPearson r of LST with each pollutant, all sensitivity versions")
print(wide.round(2).to_string())
wp = corr.pivot_table(index=["LST", "pollutant"], columns="variant", values="p_r_corrected")
print("\nCorrected p values")
print(wp.map(lambda v: "{:.3g}".format(v)).to_string())

# ---------------- CH4 and albedo ----------------
alb_rows = []
for vname in ["main_thr50", "ch4bc_thr50"]:
    if vname not in tables or "albedo_b7" not in tables[vname]:
        continue
    df = tables[vname].dropna(subset=["albedo_b7"])
    sc = rc.SpatialCorr(df["lon"].values, df["lat"].values)
    r, p, ne, pc = sc.test(df["CH4"].values, df["albedo_b7"].values, "pearson")
    row = {"variant": vname, "n": len(df), "r_CH4_albedo": r, "n_eff": ne, "p_corrected": pc}
    for lst in rc.LSTS:
        # partial correlation of LST and CH4 controlling for albedo
        def resid(v):
            b = np.polyfit(df["albedo_b7"], df[v], 1)
            return df[v] - np.polyval(b, df["albedo_b7"])
        rr, _, ne2, pc2 = sc.test(resid("CH4").values, resid(lst).values, "pearson")
        row["partial_r_{}_CH4_given_albedo".format(lst)] = rr
        row["p_corrected_partial_{}".format(lst)] = pc2
        row["r_{}_CH4".format(lst)] = np.corrcoef(df["CH4"], df[lst])[0, 1]
    alb_rows.append(row)
if alb_rows:
    alb = pd.DataFrame(alb_rows)
    alb.to_csv(os.path.join(OUT, "02_ch4_albedo.csv"), index=False)
    print("\nCH4 vs SWIR albedo (MODIS band 7) and LST-CH4 partial correlation given albedo")
    print(alb.round(3).T.to_string())

# ---------------- PCA of the pollutants only ----------------
if "main_thr50" in tables:
    df = tables["main_thr50"]
    P = rc.POLLUTANTS
    X = df[P].values
    Z = (X - X.mean(0)) / X.std(0, ddof=1)
    ev, evec = np.linalg.eigh(np.corrcoef(Z, rowvar=False))
    o = np.argsort(ev)[::-1]
    ev, evec = ev[o], evec[:, o]
    for k in range(evec.shape[1]):
        if evec[np.argmax(np.abs(evec[:, k])), k] < 0:
            evec[:, k] *= -1
    pcs = ["PC{}".format(i + 1) for i in range(len(P))]
    var = pd.DataFrame({"eigenvalue": ev, "variance_%": 100 * ev / ev.sum(),
                        "cumulative_%": np.cumsum(100 * ev / ev.sum())}, index=pcs)
    load = pd.DataFrame(evec * np.sqrt(ev), index=[rc.label(n) for n in P], columns=pcs)
    var.round(3).to_csv(os.path.join(OUT, "02_pca_pollutants_variance.csv"))
    load.round(3).to_csv(os.path.join(OUT, "02_pca_pollutants_loadings.csv"))
    n_keep = int((ev > 1).sum())
    S = Z @ evec[:, :n_keep]
    S = S / S.std(0, ddof=1)
    sc = rc.SpatialCorr(df["lon"].values, df["lat"].values)
    reg_rows = []
    for lst in rc.LSTS:
        y = df[lst].values
        yz = (y - y.mean()) / y.std(ddof=1)
        A = np.column_stack([np.ones(len(y)), S])
        beta, *_ = np.linalg.lstsq(A, yz, rcond=None)
        res = yz - A @ beta
        r2 = 1 - (res ** 2).sum() / ((yz - yz.mean()) ** 2).sum()
        for j in range(n_keep):
            r, p, ne, pc = sc.test(S[:, j], y, "pearson")
            reg_rows.append({"LST": lst, "PC": pcs[j], "std_coef": beta[j + 1], "r": r,
                             "n_eff": ne, "p_corrected": pc, "R2_all_PCs": r2})
    reg = pd.DataFrame(reg_rows)
    reg.to_csv(os.path.join(OUT, "02_lst_on_pc_scores.csv"), index=False)
    print("\nPCA of the pollutants only (eigenvalue > 1 kept: {})".format(n_keep))
    print(var.round(2).to_string())
    print(load.round(2).to_string())
    print("\nLST regressed on the PC scores")
    print(reg.round(3).to_string(index=False))

print("\nDone. Outputs in", OUT)

# -*- coding: utf-8 -*-
"""
review_04_gwr.py  -  regression models redone for the V1 review
Annual data: Review_V1/02_annual/annual_cells_main_thr50.csv (from review_02)
For LST Day and LST Night, standardised variables:
  OLS and GWR with (1) pollutants, (2) controls only, (3) pollutants + controls
  intercept-only GWR (own bandwidth and the bandwidth of the pollutant model)
  -> AICc, R2, adjusted R2 for OLS and GWR, Moran's I of residuals,
     share of cells with a significant coefficient, local condition number
  Monte Carlo test of spatial variability (fixed bandwidth, locations permuted)
  MGWR of the pollutants + controls model
Controls: NDVI, elevation, built-up fraction (GHSL 2025), distance to permanent water
Output: Review_V1/04_gwr/
Run time: GWR searches take minutes; Monte Carlo and MGWR can take hours.
Switch them off below if needed.
"""
import os
import time
import numpy as np
import pandas as pd
from mgwr.gwr import GWR, MGWR
from mgwr.sel_bw import Sel_BW
import review_common as rc
import functools
print = functools.partial(print, flush=True)   # show progress immediately

POLL = ["AOD", "CH4", "HCHO", "NO2", "CO"]
CTRL = ["NDVI", "elev", "built_ghsl", "dist_water_km"]
RUN_MONTE_CARLO = True
N_PERM = 99
RUN_MGWR = False          # MGWR is run separately on a 0.1 degree grid: review_04b_mgwr.py
OUT = rc.ensure_out("04_gwr")


def z(a):
    a = np.asarray(a, float)
    return (a - a.mean(axis=0)) / a.std(axis=0)


def aicc(rss, n, k):
    """AICc exactly as mgwr computes it for Gaussian models (k = trace of the hat matrix)."""
    llf = -np.log(rss) * n / 2.0 - (1 + np.log(np.pi / (n / 2.0))) * n / 2.0
    return -2.0 * llf + 2.0 * n * (k + 1.0) / (n - k - 2.0)


def vif(X):
    R = np.corrcoef(X, rowvar=False)
    return np.diag(np.linalg.inv(R))


def ols(y, X):
    n = len(y)
    A = np.column_stack([np.ones(n), X]) if X.shape[1] else np.ones((n, 1))
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    res = y - A @ beta
    k = A.shape[1]
    rss = (res ** 2).sum()
    tss = ((y - y.mean()) ** 2).sum()
    r2 = 1 - rss / tss
    return {"R2": r2, "adjR2": 1 - (1 - r2) * (n - 1) / (n - k), "AICc": aicc(rss, n, k),
            "k": k}, res.ravel()


def gwr_fit(coords, y, X, bw=None):
    if X.shape[1] == 0:          # intercept only
        Xg, const = np.ones((len(y), 1)), False
    else:
        Xg, const = X, True
    if bw is None:
        bw = Sel_BW(coords, y, Xg, spherical=True, kernel="bisquare", fixed=False,
                    constant=const).search()
    res = GWR(coords, y, Xg, bw, spherical=True, kernel="bisquare", fixed=False,
              constant=const).fit()
    return bw, res


def save(summ, coef_rows, mc_rows):
    """Write the tables after every model, so nothing is lost if the run is stopped."""
    pd.DataFrame(summ).to_csv(os.path.join(OUT, "04_model_summary.csv"), index=False)
    pd.DataFrame(coef_rows).to_csv(os.path.join(OUT, "04_coefficients.csv"), index=False)
    if mc_rows:
        pd.DataFrame(mc_rows).to_csv(os.path.join(OUT, "04_monte_carlo.csv"), index=False)


def main():
    path = os.path.join(rc.OUT_DIR, "02_annual", "annual_cells_main_thr50.csv")
    df = pd.read_csv(path)
    missing = [c for c in CTRL if c not in df]
    if missing:
        raise SystemExit("control variables missing in {}: {}".format(path, missing))
    df = df.dropna(subset=POLL + CTRL + rc.LSTS).reset_index(drop=True)
    coords = list(zip(df["lon"], df["lat"]))
    lon, lat = df["lon"].values, df["lat"].values
    n = len(df)
    print("cells:", n)

    models = [("pollutants", POLL), ("controls", CTRL), ("pollutants+controls", POLL + CTRL)]
    summ, coef_rows, mc_rows = [], [], []
    print("VIF pollutants+controls:", dict(zip(POLL + CTRL, np.round(vif(df[POLL + CTRL].values), 2))))

    for lst in rc.LSTS:
        y = z(df[lst].values).reshape(-1, 1)
        bw_poll = None
        for mname, cols in models:
            X = z(df[cols].values)
            t0 = time.time()
            print("{}: {} - searching bandwidth...".format(lst, mname))
            o, ores = ols(y, X)
            bw, res = gwr_fit(coords, y, X)
            if mname == "pollutants":
                bw_poll = bw
            chk = aicc(res.resid_ss, n, res.tr_S)
            sig = res.filter_tvals() != 0
            cn = res.local_collinearity()[2].ravel()   # local condition number
            mi_o, p_o = rc.morans_i(ores, lon, lat, n_perm=199)
            mi_g, p_g = rc.morans_i(res.resid_response.ravel(), lon, lat, n_perm=199)
            row = {"LST": lst, "model": mname, "n": n, "k": len(cols),
                   "max_VIF": float(vif(df[cols].values).max()) if len(cols) > 1 else 1.0,
                   "OLS_R2": o["R2"], "OLS_adjR2": o["adjR2"], "OLS_AICc": o["AICc"],
                   "OLS_resid_MoransI": mi_o, "OLS_resid_MoransI_p": p_o,
                   "GWR_bandwidth": int(bw), "GWR_R2": res.R2, "GWR_adjR2": res.adj_R2,
                   "GWR_AICc": res.aicc, "GWR_AICc_check": chk,
                   "GWR_resid_MoransI": mi_g, "GWR_resid_MoransI_p": p_g,
                   "local_CN_median": float(np.median(cn)), "local_CN_pct_gt30": float(100 * (cn > 30).mean()),
                   "seconds": time.time() - t0}
            summ.append(row)
            for j, c in enumerate(cols):
                b = res.params[:, j + 1]
                s = sig[:, j + 1]
                coef_rows.append({"LST": lst, "model": mname, "variable": c,
                                  "OLS_coef": float(np.linalg.lstsq(np.column_stack([np.ones(n), X]), y,
                                                                    rcond=None)[0][j + 1, 0]),
                                  "GWR_median": float(np.median(b)), "GWR_q25": float(np.percentile(b, 25)),
                                  "GWR_q75": float(np.percentile(b, 75)),
                                  "pct_cells_significant": float(100 * s.mean()),
                                  "pct_sig_positive": float(100 * (s & (b > 0)).mean()),
                                  "pct_sig_negative": float(100 * (s & (b < 0)).mean())})
            out = pd.DataFrame({"lon": lon, "lat": lat, "local_R2": res.localR2.ravel(), "local_CN": cn})
            for j, c in enumerate(cols):
                out["coef_" + c] = res.params[:, j + 1]
                out["sig_" + c] = sig[:, j + 1]
            out["intercept"] = res.params[:, 0]
            out.to_csv(os.path.join(OUT, "gwr_cells_{}_{}.csv".format(lst, mname.replace("+", "_"))), index=False)
            print("{} {:20s} bw {:4d}  OLS R2 {:.3f} AICc {:.0f} | GWR R2 {:.3f} adjR2 {:.3f} AICc {:.0f} "
                  "(check {:.0f})  Moran I resid OLS {:.2f} GWR {:.2f}  [{:.0f} s]".format(
                      lst, mname, int(bw), o["R2"], o["AICc"], res.R2, res.adj_R2, res.aicc, chk,
                      mi_o, mi_g, time.time() - t0))
            save(summ, coef_rows, mc_rows)

            # Monte Carlo test of spatial variability, fixed bandwidth
            if RUN_MONTE_CARLO and mname == "pollutants+controls":
                t1 = time.time()
                obs_sd = res.params.std(axis=0)
                rng = np.random.default_rng(42)
                exceed = np.zeros(res.params.shape[1])
                for i in range(N_PERM):
                    perm = rng.permutation(n)
                    pc = [coords[k] for k in perm]
                    r2 = GWR(pc, y, X, bw, spherical=True, kernel="bisquare", fixed=False).fit()
                    exceed += r2.params.std(axis=0) >= obs_sd
                    if i == 0:
                        print("   Monte Carlo: first permutation took {:.0f} s, {} permutations".format(
                            time.time() - t1, N_PERM))
                for j, c in enumerate(["intercept"] + cols):
                    mc_rows.append({"LST": lst, "model": mname, "variable": c, "observed_sd": obs_sd[j],
                                    "p_value": (exceed[j] + 1.0) / (N_PERM + 1.0), "n_perm": N_PERM,
                                    "bandwidth_fixed": int(bw)})
                save(summ, coef_rows, mc_rows)

        # intercept-only GWR: own bandwidth and the pollutant model's bandwidth
        X0 = np.zeros((n, 0))
        for tag, bwi in [("intercept-only (own bandwidth)", None),
                         ("intercept-only (bandwidth of pollutant model)", bw_poll)]:
            t0 = time.time()
            bw, res = gwr_fit(coords, y, X0, bw=bwi)
            mi_g, p_g = rc.morans_i(res.resid_response.ravel(), lon, lat, n_perm=199)
            summ.append({"LST": lst, "model": tag, "n": n, "k": 0, "GWR_bandwidth": int(bw),
                         "GWR_R2": res.R2, "GWR_adjR2": res.adj_R2, "GWR_AICc": res.aicc,
                         "GWR_resid_MoransI": mi_g, "GWR_resid_MoransI_p": p_g,
                         "seconds": time.time() - t0})
            print("{} {:45s} bw {:4d} R2 {:.3f} adjR2 {:.3f} AICc {:.0f}".format(
                lst, tag, int(bw), res.R2, res.adj_R2, res.aicc))
            save(summ, coef_rows, mc_rows)

        # MGWR
        if RUN_MGWR:
            t0 = time.time()
            cols = POLL + CTRL
            X = z(df[cols].values)
            print("{}: MGWR running (can take a long time)...".format(lst))
            sel = Sel_BW(coords, y, X, multi=True, spherical=True, kernel="bisquare", fixed=False)
            sel.search(verbose=False)
            mres = MGWR(coords, y, X, sel, spherical=True, kernel="bisquare", fixed=False, n_jobs=-1).fit()
            msig = mres.filter_tvals() != 0
            summ.append({"LST": lst, "model": "MGWR pollutants+controls", "n": n, "k": len(cols),
                         "GWR_bandwidth": -1, "GWR_R2": mres.R2, "GWR_adjR2": mres.adj_R2,
                         "GWR_AICc": mres.aicc, "seconds": time.time() - t0})
            for j, c in enumerate(["intercept"] + cols):
                b = mres.params[:, j]
                s = msig[:, j]
                coef_rows.append({"LST": lst, "model": "MGWR pollutants+controls", "variable": c,
                                  "MGWR_bandwidth": int(sel.bw[0][j]),
                                  "GWR_median": float(np.median(b)),
                                  "GWR_q25": float(np.percentile(b, 25)), "GWR_q75": float(np.percentile(b, 75)),
                                  "pct_cells_significant": float(100 * s.mean()),
                                  "pct_sig_positive": float(100 * (s & (b > 0)).mean()),
                                  "pct_sig_negative": float(100 * (s & (b < 0)).mean())})
            print("{} MGWR R2 {:.3f} adjR2 {:.3f} AICc {:.0f}, bandwidths {} [{:.0f} s]".format(
                lst, mres.R2, mres.adj_R2, mres.aicc,
                dict(zip(["intercept"] + cols, [int(b) for b in sel.bw[0]])), time.time() - t0))

        save(summ, coef_rows, mc_rows)

    print("\nModel summary")
    print(pd.DataFrame(summ).drop(columns=["seconds"]).round(3).to_string(index=False))
    print("\nCoefficients")
    print(pd.DataFrame(coef_rows).round(3).to_string(index=False))
    if mc_rows:
        print("\nMonte Carlo test of spatial variability (fixed bandwidth)")
        print(pd.DataFrame(mc_rows).round(3).to_string(index=False))
    print("\nDone. Outputs in", OUT)


if __name__ == "__main__":
    main()

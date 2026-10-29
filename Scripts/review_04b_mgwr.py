# -*- coding: utf-8 -*-
"""
review_04b_mgwr.py  -  multiscale GWR (MGWR) for the V1 review
MGWR on the full 0.05 degree grid (about 4,700 cells, 9 predictors) takes far too
long on a PC, so the annual table is averaged to a 0.1 degree grid (2 x 2 cells,
about 1,200 cells). On the same 0.1 degree data the script also fits OLS, GWR and
an intercept-only GWR, so the four models are compared on identical data.
Model: LST ~ AOD + CH4 + HCHO + NO2 + CO + NDVI + elevation + built-up + distance to water
Needs: Review_V1/02_annual/annual_cells_main_thr50.csv (from review_02)
Output: Review_V1/04b_mgwr/
"""
import os
import time
import functools
import numpy as np
import pandas as pd
from mgwr.gwr import GWR, MGWR
from mgwr.sel_bw import Sel_BW
import review_common as rc

print = functools.partial(print, flush=True)

POLL = ["AOD", "CH4", "HCHO", "NO2", "CO"]
CTRL = ["NDVI", "elev", "built_ghsl", "dist_water_km"]
BLOCK = 2            # 2 x 2 cells of 0.05 degree = 0.1 degree
MIN_CELLS = 2        # a 0.1 degree cell needs at least 2 of its 4 cells
OUT = rc.ensure_out("04b_mgwr")


def z(a):
    a = np.asarray(a, float)
    return (a - a.mean(axis=0)) / a.std(axis=0)


def aicc(rss, n, k):
    """AICc as mgwr computes it for Gaussian models (k = trace of the hat matrix)."""
    llf = -np.log(rss) * n / 2.0 - (1 + np.log(np.pi / (n / 2.0))) * n / 2.0
    return -2.0 * llf + 2.0 * n * (k + 1.0) / (n - k - 2.0)


def main():
    df = pd.read_csv(os.path.join(rc.OUT_DIR, "02_annual", "annual_cells_main_thr50.csv"))
    cols = POLL + CTRL + rc.LSTS
    df = df.dropna(subset=cols)
    c = np.round((df["lon"] - df["lon"].min()) / rc.CELL_DEG).astype(int) // BLOCK
    r = np.round((df["lat"].max() - df["lat"]) / rc.CELL_DEG).astype(int) // BLOCK
    g = df.assign(bc=c, br=r).groupby(["br", "bc"])
    agg = g[["lon", "lat"] + cols].mean()
    agg["n_cells"] = g.size()
    agg = agg[agg["n_cells"] >= MIN_CELLS].reset_index(drop=True)
    n = len(agg)
    coords = list(zip(agg["lon"], agg["lat"]))
    print("0.1 degree cells:", n)
    agg.to_csv(os.path.join(OUT, "cells_0.1deg.csv"), index=False)

    names = POLL + CTRL
    X = z(agg[names].values)
    summ, coef = [], []
    for lst in rc.LSTS:
        y = z(agg[lst].values).reshape(-1, 1)

        # OLS
        A = np.column_stack([np.ones(n), X])
        beta, *_ = np.linalg.lstsq(A, y, rcond=None)
        rss = float(((y - A @ beta) ** 2).sum())
        r2 = 1 - rss / float(((y - y.mean()) ** 2).sum())
        k = A.shape[1]
        summ.append({"LST": lst, "model": "OLS", "n": n, "R2": r2,
                     "adjR2": 1 - (1 - r2) * (n - 1) / (n - k), "AICc": aicc(rss, n, k)})

        # GWR and intercept-only GWR
        t0 = time.time()
        bw = Sel_BW(coords, y, X, spherical=True, kernel="bisquare", fixed=False).search()
        res = GWR(coords, y, X, bw, spherical=True, kernel="bisquare", fixed=False).fit()
        summ.append({"LST": lst, "model": "GWR", "n": n, "bandwidth": int(bw), "R2": res.R2,
                     "adjR2": res.adj_R2, "AICc": res.aicc})
        one = np.ones((n, 1))
        bw0 = Sel_BW(coords, y, one, spherical=True, kernel="bisquare", fixed=False, constant=False).search()
        r0 = GWR(coords, y, one, bw0, spherical=True, kernel="bisquare", fixed=False, constant=False).fit()
        summ.append({"LST": lst, "model": "GWR intercept-only", "n": n, "bandwidth": int(bw0),
                     "R2": r0.R2, "adjR2": r0.adj_R2, "AICc": r0.aicc})
        print("{}: OLS R2 {:.3f} | GWR bw {} R2 {:.3f} AICc {:.0f} | intercept-only bw {} R2 {:.3f} AICc {:.0f} [{:.0f} s]".format(
            lst, r2, int(bw), res.R2, res.aicc, int(bw0), r0.R2, r0.aicc, time.time() - t0))
        pd.DataFrame(summ).to_csv(os.path.join(OUT, "04b_model_summary.csv"), index=False)

        # MGWR
        t0 = time.time()
        print("{}: MGWR bandwidth search (progress is printed per iteration)...".format(lst))
        sel = Sel_BW(coords, y, X, multi=True, spherical=True, kernel="bisquare", fixed=False)
        sel.search(verbose=True)
        m = MGWR(coords, y, X, sel, spherical=True, kernel="bisquare", fixed=False, n_jobs=-1).fit()
        sig = m.filter_tvals() != 0
        bws = [int(b) for b in sel.bw[0]]
        summ.append({"LST": lst, "model": "MGWR", "n": n, "bandwidth": -1, "R2": m.R2,
                     "adjR2": m.adj_R2, "AICc": m.aicc})
        for j, nm in enumerate(["intercept"] + names):
            b = m.params[:, j]
            s = sig[:, j]
            coef.append({"LST": lst, "variable": nm, "MGWR_bandwidth": bws[j],
                         "median": float(np.median(b)), "q25": float(np.percentile(b, 25)),
                         "q75": float(np.percentile(b, 75)),
                         "pct_significant": float(100 * s.mean()),
                         "pct_sig_positive": float(100 * (s & (b > 0)).mean()),
                         "pct_sig_negative": float(100 * (s & (b < 0)).mean())})
        out = pd.DataFrame({"lon": agg["lon"], "lat": agg["lat"]})
        for j, nm in enumerate(["intercept"] + names):
            out["coef_" + nm] = m.params[:, j]
            out["sig_" + nm] = sig[:, j]
        out.to_csv(os.path.join(OUT, "mgwr_cells_{}.csv".format(lst)), index=False)
        pd.DataFrame(summ).to_csv(os.path.join(OUT, "04b_model_summary.csv"), index=False)
        pd.DataFrame(coef).to_csv(os.path.join(OUT, "04b_mgwr_coefficients.csv"), index=False)
        print("{}: MGWR R2 {:.3f} adjR2 {:.3f} AICc {:.0f}  bandwidths {}  [{:.0f} s]".format(
            lst, m.R2, m.adj_R2, m.aicc, dict(zip(["intercept"] + names, bws)), time.time() - t0))

    print("\nModel summary (0.1 degree grid)")
    print(pd.DataFrame(summ).round(3).to_string(index=False))
    print("\nMGWR coefficients")
    print(pd.DataFrame(coef).round(3).to_string(index=False))
    print("\nDone. Outputs in", OUT)


if __name__ == "__main__":
    main()

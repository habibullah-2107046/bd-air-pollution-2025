# -*- coding: utf-8 -*-
"""
review_common.py
Shared settings and helpers for the review_0X scripts (response to the V1 review).
Not run on its own: the other review scripts import it, so keep it in the Scripts folder.
"""
import os
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin

# ---------------- paths ----------------
ROOT     = os.environ.get("AQ_ROOT", r"C:\Users\User\Desktop\AirPollutionResearchWithRakibVai")
DATA_DIR = os.path.join(ROOT, "BD_AQ_2025")       # original exports
REV_DIR  = os.path.join(ROOT, "BD_AQ_review")     # exports of review_exports_GEE.js
OUT_DIR  = os.path.join(ROOT, "Review_V1")        # all outputs of the review scripts
BOUNDARY = os.path.join(ROOT, "boundary", "BGD_adm0.shp")

# ---------------- grid ----------------
NODATA   = -9999
CELL_DEG = 0.05      # same 0.05 degree grid as the original analysis
MIN_LAND = 0.5       # a cell needs >= 50 % land pixels

# ---------------- periods ----------------
# P0 = December 2024, P1..P12 = January..December 2025
PERIODS = ["D2412"] + ["M{:02d}".format(m) for m in range(1, 13)]
PERIOD_LABEL = ["Dec 2024"] + ["{} 2025".format(m) for m in
                ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]]
IDX_2025 = list(range(1, 13))                     # the 12 months of 2025

# Seasons: winter now uses ONE continuous winter, Dec 2024 - Feb 2025
SEASONS = [
    ("Winter",       "Dec 2024\u2013Feb 2025", [0, 1, 2]),
    ("Pre-monsoon",  "Mar\u2013May 2025",      [3, 4, 5]),
    ("Monsoon",      "Jun\u2013Sep 2025",      [6, 7, 8, 9]),
    ("Post-monsoon", "Oct\u2013Nov 2025",      [10, 11]),
]

# ---------------- variables ----------------
O3_DU = 1.0 / 4.4615e-4   # mol/m2 -> Dobson units

# One unit convention for every table and figure (review comment on units):
#   NO2, SO2, HCHO: umol/m2    CO: mmol/m2    O3: DU    CH4: ppb    LST: deg C
# name: (label, unit, factor, file of Jan-Dec 2025, file of Dec 2024, band in Dec 2024 file)
VARS = {
    "AOD":       ("AOD",       "\u2013",                     1.0,   "AOD_monthly_2025.tif",       "MODIS_2024_12.tif", "AOD"),
    "UVAI":      ("UVAI",      "\u2013",                     1.0,   "UVAI_monthly_2025.tif",      "S5P_2024_12.tif",   "UVAI"),
    "NO2":       ("NO\u2082",  "\u00b5mol m\u207b\u00b2",    1e6,   "NO2_monthly_2025.tif",       "S5P_2024_12.tif",   "NO2"),
    "SO2":       ("SO\u2082",  "\u00b5mol m\u207b\u00b2",    1e6,   "SO2_monthly_2025.tif",       "S5P_2024_12.tif",   "SO2"),
    "CO":        ("CO",        "mmol m\u207b\u00b2",         1e3,   "CO_monthly_2025.tif",        "S5P_2024_12.tif",   "CO"),
    "CH4":       ("CH\u2084",  "ppb",                        1.0,   "CH4_monthly_2025.tif",       "S5P_2024_12.tif",   "CH4"),
    "HCHO":      ("HCHO",      "\u00b5mol m\u207b\u00b2",    1e6,   "HCHO_monthly_2025.tif",      "S5P_2024_12.tif",   "HCHO"),
    "O3":        ("O\u2083",   "DU",                         O3_DU, "O3_monthly_2025.tif",        "S5P_2024_12.tif",   "O3"),
    "LST_Day":   ("LST Day",   "\u00b0C",                    1.0,   "LST_Day_monthly_2025.tif",   "MODIS_2024_12.tif", "LST_Day"),
    "LST_Night": ("LST Night", "\u00b0C",                    1.0,   "LST_Night_monthly_2025.tif", "MODIS_2024_12.tif", "LST_Night"),
}
# alternative versions used in sensitivity tests (13-band files in BD_AQ_review)
ALT_VARS = {
    "LSTAqua_Day":   ("LST Day (Aqua)",   "\u00b0C",                 1.0, "LSTAqua_Day_monthly.tif"),
    "LSTAqua_Night": ("LST Night (Aqua)", "\u00b0C",                 1.0, "LSTAqua_Night_monthly.tif"),
    "AOD_strict":    ("AOD (strict QA)",  "\u2013",                  1.0, "AOD_strictQA_monthly.tif"),
    "CH4_bc":        ("CH\u2084 (bias-corrected)", "ppb",           1.0, "CH4_bc_monthly.tif"),
    "SO2_cf03":      ("SO\u2082 (cloud fraction < 0.3)", "\u00b5mol m\u207b\u00b2", 1e6, "SO2_cf03_monthly.tif"),
}
POLLUTANTS = ["AOD", "UVAI", "NO2", "SO2", "CO", "CH4", "HCHO", "O3"]
LSTS       = ["LST_Day", "LST_Night"]
MODIS_1KM  = ["LST_Day", "LST_Night", "LSTAqua_Day", "LSTAqua_Night"]   # land-masked before gridding

# coal power plants (Global Energy Monitor / Wikipedia coordinates), used in review_07
PLANTS = [
    # name, lat, lon, MW, country
    ("Farakka (NTPC)",         24.7725, 87.8942, 2100, "India"),
    ("Sagardighi",             24.3789, 88.0956, 1600, "India"),
    ("Barapukuria",            25.5538, 88.9486,  525, "Bangladesh"),
    ("Rampal (Maitree)",       22.5925, 89.5564, 1320, "Bangladesh"),
    ("Payra (BCPCL)",          22.0000, 90.3057, 1320, "Bangladesh"),
    ("Patuakhali (RPCL-NORINCO)", 22.0251, 90.3123, 1320, "Bangladesh"),
    ("Banshkhali (S Alam)",    21.9708, 91.8889, 1320, "Bangladesh"),
    ("Matarbari",              21.7018, 91.8834, 1200, "Bangladesh"),
]


def ensure_out(*sub):
    p = os.path.join(OUT_DIR, *sub)
    if not os.path.exists(p):
        os.makedirs(p)
    return p


def read_all(path):
    """All bands of a GeoTIFF as float64 with no-data -> NaN."""
    with rasterio.open(path) as src:
        a = src.read().astype("float64")
        names = list(src.descriptions)
        a[(~np.isfinite(a)) | (a <= NODATA + 1)] = np.nan
        return a, src.transform, src.crs, src.bounds, names


def band_by_name(path, name):
    a, tr, crs, b, names = read_all(path)
    if name in names:
        return a[names.index(name)], tr
    raise KeyError("band {} not found in {} (bands: {})".format(name, path, names))


class Grid(object):
    """0.05 degree grid and land mask, built exactly as in correlation_matrix.py."""

    def __init__(self):
        a, tr, crs, b, _ = read_all(os.path.join(DATA_DIR, "AOD_annual_2025.tif"))
        self.crs = crs
        self.land = np.isfinite(a[0]).astype("float64")
        self.land_tr = tr
        self.W = int(np.ceil((b.right - b.left) / CELL_DEG))
        self.H = int(np.ceil((b.top - b.bottom) / CELL_DEG))
        self.tr = from_origin(b.left, b.top, CELL_DEG, CELL_DEG)
        self.good = self.to_grid(self.land, tr) >= MIN_LAND
        r, c = np.indices((self.H, self.W))
        self.lon = self.tr.c + (c + 0.5) * CELL_DEG
        self.lat = self.tr.f - (r + 0.5) * CELL_DEG
        self.ext = [b.left, b.left + self.W * CELL_DEG, b.top - self.H * CELL_DEG, b.top]

    def to_grid(self, arr, tr):
        out = np.full((self.H, self.W), np.nan)
        reproject(arr, out, src_transform=tr, src_crs=self.crs, dst_transform=self.tr,
                  dst_crs=self.crs, src_nodata=np.nan, dst_nodata=np.nan,
                  resampling=Resampling.average)
        return out

    def land_on(self, shape, tr):
        """Land mask brought onto another raster grid."""
        if shape == self.land.shape:
            return self.land
        lm = np.zeros(shape)
        reproject(self.land, lm, src_transform=self.land_tr, src_crs=self.crs,
                  dst_transform=tr, dst_crs=self.crs, resampling=Resampling.nearest)
        return lm

    def grid_band(self, a, tr, mask_land=False, min_valid=None):
        """
        Average a raster onto the 0.05 degree grid.
        min_valid: a cell is kept only if at least this share of its land pixels has data
        (None = keep any cell with at least one valid pixel, as in the original scripts).
        """
        a = a.copy()
        land = self.land_on(a.shape, tr) > 0
        if mask_land:
            a[~land] = np.nan
        g = self.to_grid(a, tr)
        if min_valid is not None:
            vf = self.to_grid(np.where(land, np.isfinite(a), np.nan).astype("float64"), tr)
            g[~(vf >= min_valid)] = np.nan
        g[~self.good] = np.nan
        return g

    def native_coverage(self, a, tr):
        """% of land pixels (native resolution) with data, as in the GEE coverage table."""
        land = self.land_on(a.shape, tr) > 0
        return 100.0 * np.isfinite(a[land]).sum() / max(land.sum(), 1)


MIN_VALID_FRAC = 0.5   # a monthly 0.05 degree cell needs data in >= 50 % of its land pixels


def load_cube(grid, names=None, verbose=True):
    """
    13-period cube on the 0.05 degree grid for each variable:
    cube[name] = array (13, H, W) in the unified units, and
    cube["cov__" + name] = native-resolution coverage (%) of the 13 periods.
    Main variables: Dec 2024 from BD_AQ_review, Jan-Dec 2025 from BD_AQ_2025.
    Alternative variables: all 13 bands from BD_AQ_review.
    Cached in Review_V1/cube.npz, delete that file to rebuild.
    """
    cache = os.path.join(ensure_out(), "cube.npz")
    cube = {}
    if os.path.exists(cache):
        z = np.load(cache)
        cube = {k: z[k] for k in z.files}
    wanted = names or (list(VARS) + list(ALT_VARS))
    todo = [n for n in wanted if n not in cube]
    for n in todo:
        mask_land = n in MODIS_1KM
        if n in VARS:
            label_, unit_, f, fm, fd, bd = VARS[n]
            a, tr, _, _, _ = read_all(os.path.join(DATA_DIR, fm))
            if a.shape[0] != 12:
                raise ValueError("{} has {} bands, expected 12".format(fm, a.shape[0]))
            dec, trd = band_by_name(os.path.join(REV_DIR, fd), bd)
            raw = [(dec, trd)] + [(a[i], tr) for i in range(12)]
        else:
            label_, unit_, f, fm = ALT_VARS[n]
            path = os.path.join(REV_DIR, fm)
            if not os.path.exists(path):
                if verbose:
                    print("  {}: {} not found - skipped".format(n, fm))
                continue
            a, tr, _, _, _ = read_all(path)
            if a.shape[0] != 13:
                raise ValueError("{} has {} bands, expected 13".format(fm, a.shape[0]))
            raw = [(a[i], tr) for i in range(13)]
        cube[n] = np.stack([grid.grid_band(x * f, t, mask_land, MIN_VALID_FRAC) for x, t in raw])
        cube["cov__" + n] = np.array([grid.native_coverage(x, t) for x, t in raw])
        if verbose:
            print("  gridded {}".format(n))
    if todo:
        np.savez_compressed(cache, **cube)
    return {k: v for k, v in cube.items() if k in wanted or k[5:] in wanted}


def variables(cube):
    return [k for k in cube if not k.startswith("cov__")]


def coverage(cube, name):
    """Native-resolution coverage (% of land pixels with data) of the 13 periods."""
    return cube["cov__" + name]


def label(n):
    if n in VARS:
        return VARS[n][0]
    if n in ALT_VARS:
        return ALT_VARS[n][0]
    return n


def unit(n):
    if n in VARS:
        return VARS[n][1]
    if n in ALT_VARS:
        return ALT_VARS[n][1]
    return ""


# ---------------- statistics ----------------
def km_xy(lon, lat):
    """Equirectangular projection to km (fine for distances inside Bangladesh)."""
    lat0 = np.radians(23.7)
    x = 6371.0 * np.radians(lon) * np.cos(lat0)
    y = 6371.0 * np.radians(lat)
    return x, y


class SpatialCorr(object):
    """
    Distance classes for one set of cells, reused for every variable.
    Effective sample size of a correlation between two spatially autocorrelated
    variables (Clifford, Richardson & Hemon 1989; Dutilleul 1993):
        var(r) ~ 1/n + sum_k N_k I_x(k) I_y(k) / n^2,    n_eff = 1 + 1/var(r)
    N_k = number of ordered pairs in distance class k, I = Moran's I in class k.
    Equal-width distance classes, number of classes from Sturges' rule.
    """

    def __init__(self, lon, lat, chunk=500):
        X, Y = km_xy(np.asarray(lon, float), np.asarray(lat, float))
        self.n = n = len(X)
        dmax = np.hypot(X.max() - X.min(), Y.max() - Y.min()) + 1e-6
        self.n_class = int(np.ceil(1 + np.log2(n * (n - 1) / 2.0)))
        edges = np.linspace(0, dmax, self.n_class + 1)
        self.blocks = []
        self.cnt = np.zeros(self.n_class)
        for i0 in range(0, n, chunk):
            i1 = min(n, i0 + chunk)
            d = np.hypot(X[i0:i1, None] - X[None, :], Y[i0:i1, None] - Y[None, :])
            k = np.clip(np.searchsorted(edges, d, side="right") - 1, 0, self.n_class - 1).astype(np.int16)
            idx = np.arange(i0, i1)
            k[idx - i0, idx] = -1                      # drop self pairs
            self.blocks.append((i0, i1, k))
            kk = k[k >= 0]
            self.cnt += np.bincount(kk, minlength=self.n_class)
        self._cache = {}

    def moran_by_class(self, v, key=None):
        if key is not None and key in self._cache:
            return self._cache[key]
        v = np.asarray(v, float)
        z = (v - v.mean()) / v.std()
        s = np.zeros(self.n_class)
        for i0, i1, k in self.blocks:
            prod = z[i0:i1, None] * z[None, :]
            m = k >= 0
            s += np.bincount(k[m], weights=prod[m], minlength=self.n_class)
        I = np.where(self.cnt > 0, s / np.maximum(self.cnt, 1), 0.0)
        if key is not None:
            self._cache[key] = I
        return I

    def n_eff(self, x, y, kx=None, ky=None):
        Ix = self.moran_by_class(x, kx)
        Iy = self.moran_by_class(y, ky)
        var_r = 1.0 / self.n + np.sum(self.cnt * Ix * Iy) / float(self.n) ** 2
        ne = 1.0 + 1.0 / var_r if var_r > 0 else self.n
        return float(min(max(ne, 3.0), self.n))

    def test(self, x, y, method="pearson", kx=None, ky=None):
        """r (or rho), naive p, effective n, corrected p."""
        from scipy import stats
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        if method == "pearson":
            r, p = stats.pearsonr(x, y)
            xr, yr = x, y
        else:
            r, p = stats.spearmanr(x, y)
            xr, yr = stats.rankdata(x), stats.rankdata(y)
            kx = None if kx is None else kx + "_rank"
            ky = None if ky is None else ky + "_rank"
        ne = self.n_eff(xr, yr, kx, ky)
        df = max(ne - 2.0, 1.0)
        t = r * np.sqrt(df / max(1e-12, 1 - r * r))
        return float(r), float(p), ne, float(2 * stats.t.sf(abs(t), df))


def knn_weights(lon, lat, k=8):
    """Row-standardised k-nearest-neighbour weights as (index array, weight)."""
    from scipy.spatial import cKDTree
    X, Y = km_xy(np.asarray(lon), np.asarray(lat))
    _, idx = cKDTree(np.column_stack([X, Y])).query(np.column_stack([X, Y]), k=k + 1)
    return idx[:, 1:], 1.0 / k


def morans_i(v, lon, lat, k=8, n_perm=999, seed=1):
    """Moran's I with k-nearest-neighbour weights and a permutation p value."""
    v = np.asarray(v, float).ravel()
    nb, w = knn_weights(lon, lat, k)
    z = v - v.mean()
    I = (z * (w * z[nb].sum(axis=1))).sum() / (z * z).sum()
    rng = np.random.default_rng(seed)
    perm = np.empty(n_perm)
    for i in range(n_perm):
        zp = rng.permutation(z)
        perm[i] = (zp * (w * zp[nb].sum(axis=1))).sum() / (zp * zp).sum()
    p = (np.sum(np.abs(perm) >= abs(I)) + 1.0) / (n_perm + 1.0)
    return float(I), float(p)


def fmt(v, d=3):
    if v is None or not np.isfinite(v):
        return "n/a"
    return "{:.{}f}".format(v, d)

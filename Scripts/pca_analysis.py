import matplotlib
matplotlib.use("Agg")

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------------- settings ----------------
YEAR    = 2025
ROOT    = r"C:\Users\User\Desktop\AirPollutionResearchWithRakibVai"
OUT_DIR = os.path.join(ROOT, "Figures_v2", "Statistics")
CSV     = os.path.join(OUT_DIR, "annual_pixels_{}.csv".format(YEAR))
# ------------------------------------------

# UHI is left out: it is calculated from LST Day, so it would double-count
VARS = ["LST Day", "LST Night", "AOD", "UVAI", u"CH\u2084", "HCHO", u"O\u2083"]

df = pd.read_csv(CSV, encoding="utf-8")
X = df[VARS].values
Z = (X - X.mean(axis=0)) / X.std(axis=0, ddof=1)        # standardise

# PCA from the correlation matrix
R = np.corrcoef(Z, rowvar=False)
eigval, eigvec = np.linalg.eigh(R)
order = np.argsort(eigval)[::-1]
eigval, eigvec = eigval[order], eigvec[:, order]

# make the largest loading of each PC positive (sign is arbitrary)
for k in range(eigvec.shape[1]):
    if eigvec[np.argmax(np.abs(eigvec[:, k])), k] < 0:
        eigvec[:, k] *= -1

k_all   = len(VARS)
pcs     = ["PC{}".format(i + 1) for i in range(k_all)]
var_pct = 100 * eigval / eigval.sum()
cum_pct = np.cumsum(var_pct)
loadings = eigvec * np.sqrt(eigval)                      # correlation of variable with PC
scores   = Z @ eigvec

# ---- tables ----
if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR)

summary = pd.DataFrame({"eigenvalue": eigval, "variance_%": var_pct,
                        "cumulative_%": cum_pct}, index=pcs)
summary.round(3).to_csv(os.path.join(OUT_DIR, "pca_variance_{}.csv".format(YEAR)))

load_df = pd.DataFrame(loadings, index=VARS, columns=pcs)
load_df.round(3).to_csv(os.path.join(OUT_DIR, "pca_loadings_{}.csv".format(YEAR)))

score_df = pd.DataFrame(scores[:, :3], columns=pcs[:3])
score_df.insert(0, "lat", df["lat"])
score_df.insert(0, "lon", df["lon"])
score_df.to_csv(os.path.join(OUT_DIR, "pca_scores_{}.csv".format(YEAR)), index=False)

# ---- figure 1: scree plot ----
fig, ax = plt.subplots(figsize=(5.5, 3.8))
ax.bar(range(1, k_all + 1), var_pct, color="#7fb03a", edgecolor="black", lw=0.5,
       label="Explained variance")
ax.plot(range(1, k_all + 1), cum_pct, "o-", color="#a32020", ms=4, lw=1.2,
        label="Cumulative")
for i, v in enumerate(var_pct):
    ax.text(i + 1, v + 1.5, "{:.1f}".format(v), ha="center", fontsize=7)
ax.set_xticks(range(1, k_all + 1))
ax.set_xticklabels(pcs, fontsize=8)
ax.set_ylabel("Variance explained (%)", fontsize=8)
ax.set_ylim(0, 105)
ax.tick_params(labelsize=7)
ax.legend(fontsize=7, loc="center right")
ax.set_title("PCA scree plot, annual {} (n = {} cells)".format(YEAR, len(df)), fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "pca_scree_{}.png".format(YEAR)), dpi=400)
plt.close(fig)

# ---- figure 2: biplot PC1 vs PC2 ----
fig, ax = plt.subplots(figsize=(5.8, 5.4))
s = scores[:, :2] / np.abs(scores[:, :2]).max(axis=0)    # scale scores to [-1, 1]
ax.scatter(s[:, 0], s[:, 1], s=2, color="grey", alpha=0.25, edgecolors="none",
           rasterized=True)
for i, v in enumerate(VARS):
    x, y = loadings[i, 0], loadings[i, 1]
    colour = "#2c5f9e" if "LST" in v else "#a32020"
    ax.arrow(0, 0, x, y, color=colour, width=0.004, head_width=0.03,
             length_includes_head=True, zorder=3)
    ax.text(x * 1.12, y * 1.12, v, color=colour, fontsize=8, fontweight="bold",
            ha="center", va="center", zorder=4)
circle = plt.Circle((0, 0), 1, fill=False, ls="--", lw=0.6, color="black")
ax.add_patch(circle)
ax.axhline(0, color="black", lw=0.4)
ax.axvline(0, color="black", lw=0.4)
ax.set_xlim(-1.2, 1.2)
ax.set_ylim(-1.2, 1.2)
ax.set_aspect("equal")
ax.set_xlabel("PC1 ({:.1f}%)".format(var_pct[0]), fontsize=8)
ax.set_ylabel("PC2 ({:.1f}%)".format(var_pct[1]), fontsize=8)
ax.tick_params(labelsize=7)
ax.set_title("PCA biplot, annual {}".format(YEAR), fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "pca_biplot_{}.png".format(YEAR)), dpi=400)
plt.close(fig)

# ---- print ----
print("Explained variance:")
print(summary.round(2).to_string())
print("\nLoadings (PC1-PC3):")
print(load_df.iloc[:, :3].round(2).to_string())
print("\nPCs with eigenvalue > 1 (Kaiser rule): {}".format(int((eigval > 1).sum())))
print("Done.")

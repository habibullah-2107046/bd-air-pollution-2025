import os, glob
import pandas as pd

ROOT = r"C:\Users\User\Desktop\AirPollutionResearchWithRakibVai"
IN_DIR = os.path.join(ROOT, "BD_AQ_2025")          # where the CSVs were downloaded
OUT = os.path.join(IN_DIR, "AOD_daily_2025.csv")

files = sorted(glob.glob(os.path.join(IN_DIR, "AOD_daily_2025_M*.csv")))
print("found {} monthly files".format(len(files)))

df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
df["date"] = pd.to_datetime(df["date"])
df = df.sort_values("date").reset_index(drop=True)

print("days: {}  (expected 365)".format(len(df)))
print("days with any AOD:      {}".format(df["mean_value"].notna().sum()))
print("days with >=50% cover:  {}".format((df["pct_valid"] >= 50).sum()))

df.to_csv(OUT, index=False)
print("saved: " + OUT)

# Satellite-based air pollution and land surface temperature in Bangladesh (2025)

Scripts, tables and figures for a national-scale study of air pollutants and land surface temperature (LST) over Bangladesh for 1 January – 31 December 2025. Satellite data are processed in Google Earth Engine (GEE); maps and statistics are produced in Python. The workflow adapts the approach of Sameh et al. (2026, *The Egyptian Journal of Remote Sensing and Space Sciences*, 29(1), 158–178, https://doi.org/10.1016/j.ejrs.2026.01.008) to Bangladesh.

**Status:** work in progress. The study has not yet been peer-reviewed or published, so the results are preliminary and may change. Please contact the maintainer before using or citing them.

## Key findings

- **Two spatial patterns.** AOD and CO follow a regional west–east gradient; NO2, HCHO and CH4 are concentrated over the Dhaka conurbation, where the annual NO2 column reaches 189 µmol/m² against a national mean of 36.6 µmol/m². Surface urban heat island intensity reaches 3.7 °C.
- **Strong seasonality.** NO2, CO, SO2 and AOD are highest in winter and the pre-monsoon season and lowest in the monsoon. Monsoon cloud leaves too few CH4, AOD and LST observations from June to September.
- **LST and pollutants (annual).** Daytime LST is positively correlated with CH4 (r = 0.48), NO2 (0.36), SO2 (0.31) and HCHO (0.30), but no single pollutant explains more than 23% of its spatial variance.
- **Relationships change with season.** The AOD–daytime LST correlation is −0.52 in winter and +0.29 in the monsoon.
- **PCA.** Three components explain 74.9% of the variance of ten variables.
- **GWR.** Explained variance of LST rises from 0.24 (global OLS) to 0.91 (day) and 0.93 (night), so the relationships vary strongly across the country.
- **Ground comparison.** Station PM2.5 correlates with AOD (r = 0.50) and UVAI (r = 0.72). Every DoE station exceeds the WHO guideline and the national standard for PM2.5 (station means 39.8–141.0 µg/m³).

## Workflow diagram

![Workflow of the study](Figures_v2/Flowchart/methodology_flowchart_v2.png)

## Variables and data sources

| Variable | Product | Export scale |
|---|---|---|
| CH4, HCHO, O3, UVAI | Sentinel-5P TROPOMI OFFL L3 | 1,113 m |
| NO2, CO, SO2 | Sentinel-5P TROPOMI OFFL L3 | 1,113 m |
| AOD (550 nm) | MODIS MCD19A2 (MAIAC, C6.1) | 1,000 m |
| LST Day / Night | MODIS MOD11A1 (C6.1) | 1,000 m |
| Land cover (urban, water) | MODIS MCD12Q1 | 500 m |
| Ground PM2.5, PM10, NO2, SO2, CO, O3 | DoE CAMS monthly reports, 16 stations | – |
| National boundary | `BGD_adm0.shp` (GADM, level 0) | – |
| Neighbouring countries (study area map) | Natural Earth, 1:50 m and 1:110 m | – |

Ground data come from the monthly air quality reports of the Department of Environment (DoE): <https://doe.gov.bd/pages/static-pages/6922db96933eb65569e0b11a>. Reports were available for 11 months of 2025; the May report could not be obtained. The values were compiled into `BD_AQ_2025/doe_cams_monthly_2025.csv`.

UHI is computed from annual daytime LST as pixel LST minus the mean rural LST within 10 km.

## Folder structure

```
bd-air-pollution-2025/
├── README.md
├── .gitignore
├── environment.yml            Python environment (package versions)
├── Script for GEE/            GEE scripts (JavaScript): raster and daily CSV exports
├── Scripts/                   Python scripts for maps and statistics
├── BD_AQ_2025/                Tables exported from GEE and station data (CSV only)
├── boundary/                  Natural Earth country boundaries (see "Data not in this repository")
├── Figures_v2/                Final figures
│   ├── Annual/                Annual maps
│   ├── Monthly/               Monthly 4 × 3 panel maps
│   ├── Seasonal/              Seasonal maps, means and correlations
│   ├── Temporal/              Daily trends, monthly LST
│   ├── Statistics/            Correlation, regression, PCA, GWR
│   ├── StudyArea/             Study area map
│   └── Flowchart/             Workflow diagram
├── Station_validation/        Satellite vs DoE station comparison (tables and figure)
└── Writtings/                 AI-generated reference draft of the manuscript and the progress report (not the final manuscript)
```

## Data not in this repository

| Data | Why | Where to get it | Where to put it |
|---|---|---|---|
| Raster files (`.tif`): monthly (12-band) and annual composites of every variable | Too large for GitHub | [Google Drive folder](https://drive.google.com/drive/folders/1fJYs1cqs_V31OdvbfbCfgoshujVFL4cM?usp=sharing), or re-export with the GEE scripts | `BD_AQ_2025/` |
| National boundary `BGD_adm0.shp` (with `.dbf`, `.shx`, `.prj`, `.cpg`) | GADM data may not be redistributed | <https://gadm.org> (Bangladesh, level 0); save it under this name | `boundary/` |

The GEE scripts clip to the same national boundary, uploaded to GEE as an asset. Upload the boundary to your own GEE account and change the asset path at the top of each script.

## Workflow

### 1. Google Earth Engine

| Script | Purpose |
|---|---|
| `test-gee.js` | Exploratory tests: coverage, quality masks, compositing, UHI method |
| `for_all_Pollutants_GEE_script.js` | Main export: monthly (12-band) and annual rasters, UHI, coverage table, daily national means |
| `aod_daily_fix.js` | AOD daily means as 12 monthly tasks |
| `aod_check.js` | Visual check of AOD observation count and water mask (no exports) |
| `NO2_CO_SO2_GEE_script.js` | Monthly, annual and daily exports for NO2, CO and SO2 |

All rasters: EPSG:4326, no-data = −9999, Drive folder `BD_AQ_2025`.

### 2. Python

Run in this order:

| Step | Script | Output |
|---|---|---|
| 1 | `merge_aod_daily.py` | `AOD_daily_2025.csv` |
| 2 | `study_area_map.py` | Study area map |
| 3 | `panel_map.py` | Monthly 4 × 3 panel maps |
| 4 | `annual_map.py` | Annual maps |
| 5 | `temporal_trends.py` | Daily trends, monthly LST day/night |
| 6 | `correlation_matrix.py` | 0.05° pixel table, Pearson and Spearman matrices |
| 7 | `regression_plots.py` | LST vs pollutant regressions |
| 8 | `pca_analysis.py` | PCA scree plot, biplot, loadings |
| 9 | `gwr_analysis.py` | GWR coefficient maps and summary (main model and UVAI check model) |
| 10 | `seasonal_analysis.py` | Seasonal maps, seasonal means, seasonal LST–pollutant correlation |
| 11 | `station_validation.py` | Satellite vs station comparison, WHO / Bangladesh standard table |

Steps 7–9 read the pixel table made in step 6, so run step 6 first.

## How to run

1. Download the rasters and the boundary shapefile (see "Data not in this repository").
2. Create the Python environment. On Windows:

```
conda env create -f environment.yml
```

On other systems, install the main packages yourself: Python 3.11, `rasterio` 1.4, `geopandas` 1.1, `pandas` 3.0, `numpy` 2.4, `scipy` 1.17, `matplotlib` 3.11, `scikit-learn` 1.9, `mgwr` 2.2, `matplotlib-scalebar` 0.9.

3. Open each script and set `ROOT` (near the top) to your own project folder.
4. Run the scripts in the order given above, for example:

```
conda activate gis
python panel_map.py
```

Each script prints `Done.` when it finishes.

## Output tables

| File | Content |
|---|---|
| `BD_AQ_2025/<VAR>_daily_2025.csv` | Daily national mean of each variable: `date`, `mean_value`, `pct_valid` (% of the country with valid pixels), `n_images` |
| `BD_AQ_2025/doe_cams_monthly_2025.csv` | DoE station data: `station`, `lat`, `lon`, `month`, `parameter`, `mean` (monthly mean), `capture_pct` (% data capture) |
| `Figures_v2/Temporal/monthly_means_2025.csv` | Monthly national mean, standard deviation and spatial coverage (%) of each variable |
| `Figures_v2/Temporal/daily_summary_2025.csv` | Number of days used and mean, minimum and maximum of the daily series |
| `Figures_v2/Statistics/annual_pixels_2025.csv` | Analysis table: annual value of every variable on the 0.05° grid (4,901 cells; `lon`, `lat` = cell centre) |
| `Figures_v2/Statistics/corr_pearson_r/p_2025.csv`, `corr_spearman_r/p_2025.csv` | Correlation coefficients and p-values between all variables |
| `Figures_v2/Statistics/regression_LST_pollutants_2025.csv` | Simple linear regression of LST Day and LST Night on each variable: slope, intercept, r, R², p |
| `Figures_v2/Statistics/pca_variance_2025.csv`, `pca_loadings_2025.csv`, `pca_scores_2025.csv` | PCA eigenvalues and explained variance, loadings, and cell scores (PC1–PC3) |
| `Figures_v2/Statistics/gwr_summary_2025.csv` | GWR model summary: bandwidth, OLS R², GWR R², adjusted R², AICc, maximum VIF |
| `Figures_v2/Statistics/gwr_LST_Day_2025.csv`, `gwr_LST_Night_2025.csv` | GWR results per cell: `local_R2`, standardised coefficient (`coef_*`) and significance after multiple-testing correction (`sig_*`) |
| `Figures_v2/Statistics/gwr_UVAIcheck_*.csv` | The same for the check model with UVAI in place of AOD |
| `Figures_v2/Seasonal/seasonal_means_2025.csv` | Seasonal national mean, standard deviation, range and coverage of each variable |
| `Figures_v2/Seasonal/seasonal_correlation_2025.csv` | Pearson and Spearman correlation between LST and each variable, by season |
| `Station_validation/station_validation_pairs.csv` | Paired station and satellite values for each station-month |
| `Station_validation/station_validation_summary.csv` | Correlation between station and satellite values (station-months and station means) |
| `Station_validation/station_who_comparison.csv` | Station mean PM2.5 and PM10 against the WHO guideline and the Bangladesh standard |

Units: AOD and UVAI unitless; CH4 in ppb; CO in mol/m²; O3 in Dobson Units; LST and UHI in °C; station PM2.5 and PM10 in µg/m³. HCHO, NO2 and SO2 differ between tables:

| Table | HCHO | NO2, SO2 |
|---|---|---|
| `annual_pixels_2025.csv` | mol/m² | mol/m² |
| `seasonal_means_2025.csv` (has a `unit` column) | mol/m² | µmol/m² |
| `monthly_means_2025.csv`, `daily_summary_2025.csv` | ×10⁻⁴ mol/m² | µmol/m² |

## Processing notes

- AOD: values ≤ 0 removed, at least 2 observations per pixel per month, median composite, inland water masked.
- LST: only good-quality pixels with LST error ≤ 2 K.
- Months with < 50% valid pixels are labelled as insufficient data (monsoon cloud).
- One shared colour scale per variable across all months.
- O3 is converted from mol/m² to Dobson Units for mapping.
- NO2 and SO2 are shown in µmol/m² (mol/m² × 10⁶); SO2 values below −0.001 mol/m² are removed as noise.
- PM2.5 and PM10 are not estimated from satellite data; AOD is analysed directly.
- Seasons: winter (Dec–Feb, using Jan, Feb and Dec 2025), pre-monsoon (Mar–May), monsoon (Jun–Sep), post-monsoon (Oct–Nov).
- GWR: LST Day and LST Night ~ AOD + CH4 + HCHO + NO2 + CO (adaptive bisquare kernel, AICc bandwidth). UVAI, O3 and SO2 are left out of the main model; a second model with UVAI in place of AOD is run as a check.
- Station comparison: 3 × 3 pixel mean at each DoE station, station-months with at least 50% data capture.

## Data sources and acknowledgement

- Sentinel-5P TROPOMI products: European Space Agency / Copernicus programme.
- MODIS MCD19A2, MOD11A1 and MCD12Q1 products: NASA.
- Processing platform: Google Earth Engine.
- Ground measurements: Department of Environment (DoE), Government of Bangladesh.
- Country boundaries: GADM and Natural Earth.

Each dataset remains under the terms of its provider.

## Use and citation

No licence has been set yet. Until the paper is published, please ask the maintainer before reusing the code, tables or figures. A citation will be added here after publication.

## Research team

| Name | Role in this project | Affiliation |
|---|---|---|
| [Abdullah Al Rakib](https://www.linkedin.com/in/abdullahal-rakib/) | Study lead and supervision | Research Associate, BRAC James P Grant School of Public Health, BRAC University, Dhaka, Bangladesh |
| [Rubayet Arafin Rimon](https://www.linkedin.com/in/rubayet-arafin-rimon-255203361/) | Supervision and review | Graduate Teaching Assistant, Department of Geography and Environmental Studies, Texas Tech University, USA |
| [Md. Habibullah Masbah](https://www.linkedin.com/in/habibullah046/) | Data processing, analysis, code and repository maintenance | Undergraduate student (BURP), Department of Urban and Regional Planning, Rajshahi University of Engineering & Technology (RUET), Bangladesh |
| [Md Khadem Ali](https://www.linkedin.com/in/md-khadem-ali-10849b384/) | Literature review and manuscript writing | Undergraduate student (B.Sc. in Geography and Environment), National University, Bangladesh |

## Contact

**[Md. Habibullah Masbah](https://www.linkedin.com/in/habibullah046/)** (repository maintainer)\
Department of Urban and Regional Planning, Rajshahi University of Engineering & Technology (RUET), Rajshahi, Bangladesh\
Email: habibullah.ruet.urp@gmail.com

For questions about the code or data, you can also open an issue in this repository.

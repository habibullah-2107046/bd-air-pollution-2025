// ============================================================
// review_exports_GEE.js
// Extra exports needed to answer Rimon vai's review (V1).
// Same settings as the original scripts: EPSG:4326, no-data = -9999.
// All files go to a NEW Drive folder: BD_AQ_review
//
// 13-band stacks use the band names D2412 (December 2024) and M01..M12 (2025).
//
// Tasks (9):
//   1  S5P_2024_12            Dec 2024 composites: CH4, HCHO, O3, UVAI, NO2, CO, SO2   (winter fix)
//   2  MODIS_2024_12          Dec 2024 composites: AOD, LST_Day, LST_Night (Terra)     (winter fix)
//   3  LSTAqua_Day_monthly    Aqua MYD11A1 daytime LST, 13 bands                       (Terra drift)
//   4  LSTAqua_Night_monthly  Aqua MYD11A1 night-time LST, 13 bands                    (Terra drift)
//   5  AOD_strictQA_monthly   MAIAC AOD, clear + not cloud-adjacent + best QA          (June AOD)
//   6  CH4_bc_monthly         TROPOMI CH4, bias-corrected band                         (albedo artefact)
//   7  SO2_cf03_monthly       TROPOMI SO2, only cloud fraction < 0.3                    (SO2 swing)
//   8  controls_1km           NDVI, elevation, built-up and water fractions,
//                             MCD12Q1 urban/water fractions, SWIR albedo                (controls, UHI)
//   9  overpass_time          monthly mean view time of Terra and Aqua (table)        (Terra drift)
// ============================================================

var bd = ee.FeatureCollection('projects/ee-2107046habibullah/assets/BangladeshBoundary').geometry();
var FOLDER = 'BD_AQ_review';
var NODATA = -9999;

// ---------- periods: Dec 2024 + Jan..Dec 2025 ----------
var periods = [{y: 2024, m: 12, name: 'D2412'}];
for (var mm = 1; mm <= 12; mm++) {
  periods.push({y: 2025, m: mm, name: 'M' + (mm < 10 ? '0' + mm : mm)});
}
var allStart = ee.Date.fromYMD(2024, 12, 1);
var allEnd = ee.Date.fromYMD(2026, 1, 1);

// ---------- land cover used in the original scripts (latest MCD12Q1 year) ----------
var lcImg = ee.ImageCollection('MODIS/061/MCD12Q1').sort('system:time_start', false).first();
var lc = lcImg.select('LC_Type1');
var notWater = lc.neq(17);
print('MCD12Q1 year used in the original scripts:',
      ee.Date(lcImg.get('system:time_start')).format('YYYY'));

// ---------- helpers ----------
function composite(col, s, e, stat, minObs) {
  var c = col.filterDate(s, e);
  var img = (stat === 'median') ? c.median() : c.mean();
  return img.updateMask(c.count().gte(minObs)).rename('v');
}

function stack13(col, stat, minObs) {
  var imgs = periods.map(function (p) {
    var s = ee.Date.fromYMD(p.y, p.m, 1);
    return composite(col, s, s.advance(1, 'month'), stat, minObs).rename(p.name);
  });
  return ee.Image.cat(imgs).toFloat().clip(bd).unmask(NODATA);
}

function exportImg(img, name, scale) {
  Export.image.toDrive({
    image: img, description: name, folder: FOLDER, fileNamePrefix: name,
    region: bd.bounds(), scale: scale, crs: 'EPSG:4326', maxPixels: 1e13
  });
}

function s5p(id, band, minValid) {
  return ee.ImageCollection(id).filterDate(allStart, allEnd).filterBounds(bd)
    .map(function (img) {
      var x = img.select([band], ['v']);
      if (minValid !== null) { x = x.updateMask(x.gte(minValid)); }
      return x.set('system:time_start', img.get('system:time_start'));
    });
}

// AOD exactly as in the original script (values > 0, inland water masked)
function aodPrep(img) {
  var a = img.select('Optical_Depth_055').multiply(0.001);
  return a.updateMask(a.gt(0)).updateMask(notWater).rename('v')
    .set('system:time_start', img.get('system:time_start'));
}

// AOD with strict QA: cloud mask = clear (1), adjacency = normal (0), AOD QA = best (0)
function aodStrict(img) {
  var qa = img.select('AOD_QA');
  var clear = qa.bitwiseAnd(7).eq(1);
  var noAdj = qa.rightShift(5).bitwiseAnd(7).eq(0);
  var best = qa.rightShift(8).bitwiseAnd(15).eq(0);
  return aodPrep(img).updateMask(clear.and(noAdj).and(best));
}

// LST exactly as in the original script (mandatory QA <= 1 and LST error <= 2 K)
function lstPrep(band, qcBand) {
  return function (img) {
    var t = img.select(band).multiply(0.02).subtract(273.15);
    var qc = img.select(qcBand).unmask(0);
    var ok = qc.bitwiseAnd(3).lte(1).and(qc.rightShift(6).bitwiseAnd(3).lte(1));
    return t.updateMask(ok).rename('v')
      .set('system:time_start', img.get('system:time_start'));
  };
}

var aodCol = ee.ImageCollection('MODIS/061/MCD19A2_GRANULES')
  .filterDate(allStart, allEnd).filterBounds(bd);
var terra = ee.ImageCollection('MODIS/061/MOD11A1').filterDate(allStart, allEnd).filterBounds(bd);
var aqua = ee.ImageCollection('MODIS/061/MYD11A1').filterDate(allStart, allEnd).filterBounds(bd);

// ============================================================
// 1. Dec 2024 composites of the TROPOMI variables (same settings as the originals)
// ============================================================
var dS = ee.Date.fromYMD(2024, 12, 1);
var dE = ee.Date.fromYMD(2025, 1, 1);
var S5P = [
  ['CH4',  'COPERNICUS/S5P/OFFL/L3_CH4',    'CH4_column_volume_mixing_ratio_dry_air', null],
  ['HCHO', 'COPERNICUS/S5P/OFFL/L3_HCHO',   'tropospheric_HCHO_column_number_density', null],
  ['O3',   'COPERNICUS/S5P/OFFL/L3_O3',     'O3_column_number_density', null],
  ['UVAI', 'COPERNICUS/S5P/OFFL/L3_AER_AI', 'absorbing_aerosol_index', null],
  ['NO2',  'COPERNICUS/S5P/OFFL/L3_NO2',    'tropospheric_NO2_column_number_density', null],
  ['CO',   'COPERNICUS/S5P/OFFL/L3_CO',     'CO_column_number_density', null],
  ['SO2',  'COPERNICUS/S5P/OFFL/L3_SO2',    'SO2_column_number_density', -0.001]
];
var s5pDec = ee.Image.cat(S5P.map(function (r) {
  return composite(s5p(r[1], r[2], r[3]), dS, dE, 'mean', 1).rename(r[0]);
}));
exportImg(s5pDec.toFloat().clip(bd).unmask(NODATA), 'S5P_2024_12', 1113.2);

// ============================================================
// 2. Dec 2024 composites of AOD and Terra LST (same settings as the originals)
// ============================================================
var modisDec = ee.Image.cat([
  composite(aodCol.map(aodPrep), dS, dE, 'median', 2).rename('AOD'),
  composite(terra.map(lstPrep('LST_Day_1km', 'QC_Day')), dS, dE, 'mean', 1).rename('LST_Day'),
  composite(terra.map(lstPrep('LST_Night_1km', 'QC_Night')), dS, dE, 'mean', 1).rename('LST_Night')
]);
exportImg(modisDec.toFloat().clip(bd).unmask(NODATA), 'MODIS_2024_12', 1000);

// ============================================================
// 3-4. Aqua LST (same QA as Terra)
// ============================================================
exportImg(stack13(aqua.map(lstPrep('LST_Day_1km', 'QC_Day')), 'mean', 1), 'LSTAqua_Day_monthly', 1000);
exportImg(stack13(aqua.map(lstPrep('LST_Night_1km', 'QC_Night')), 'mean', 1), 'LSTAqua_Night_monthly', 1000);

// ============================================================
// 5. AOD with strict QA
// ============================================================
exportImg(stack13(aodCol.map(aodStrict), 'median', 2), 'AOD_strictQA_monthly', 1000);

// ============================================================
// 6. CH4 bias-corrected band
// ============================================================
exportImg(stack13(s5p('COPERNICUS/S5P/OFFL/L3_CH4',
  'CH4_column_volume_mixing_ratio_dry_air_bias_corrected', null), 'mean', 1), 'CH4_bc_monthly', 1113.2);

// ============================================================
// 7. SO2 with cloud fraction < 0.3
// ============================================================
var so2cf = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_SO2').filterDate(allStart, allEnd).filterBounds(bd)
  .map(function (img) {
    var x = img.select(['SO2_column_number_density'], ['v']);
    x = x.updateMask(x.gte(-0.001)).updateMask(img.select('cloud_fraction').lt(0.3));
    return x.set('system:time_start', img.get('system:time_start'));
  });
exportImg(stack13(so2cf, 'mean', 1), 'SO2_cf03_monthly', 1113.2);

// ============================================================
// 8. Control variables and land-cover fractions at 1 km
// ============================================================
var p1km = ee.Projection('EPSG:4326').atScale(1000);
function toKm(img) {
  return img.reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 65535}).reproject(p1km);
}
var y25s = ee.Date.fromYMD(2025, 1, 1);
var y25e = ee.Date.fromYMD(2026, 1, 1);

var ndvi = ee.ImageCollection('MODIS/061/MOD13A2').filterDate(y25s, y25e)
  .select('NDVI').mean().multiply(0.0001).rename('NDVI');
var elev = toKm(ee.Image('USGS/SRTMGL1_003').select('elevation').toFloat()).rename('elev');
var wc = ee.ImageCollection('ESA/WorldCover/v200').first().select('Map');
var builtWC = toKm(wc.eq(50).toFloat()).rename('built_wc');
var waterWC = toKm(wc.eq(80).toFloat()).rename('water_wc');
var ghsl = ee.Image('JRC/GHSL/P2023A/GHS_BUILT_S/2025').select('built_surface');
var builtGHSL = toKm(ghsl.divide(10000).toFloat()).rename('built_ghsl');
var urbanMCD = toKm(lc.eq(13).toFloat()).rename('urban_mcd');
var waterMCD = toKm(lc.eq(17).toFloat()).rename('water_mcd');
var alb = ee.ImageCollection('MODIS/061/MCD43A3').filterDate(y25s, y25e)
  .map(function (img) {
    var ok = img.select('BRDF_Albedo_Band_Mandatory_Quality_Band7').eq(0);
    return img.select('Albedo_WSA_Band7').multiply(0.001).updateMask(ok);
  }).mean().rename('albedo_b7');

var controls = ee.Image.cat([ndvi, elev, builtWC, waterWC, builtGHSL, urbanMCD, waterMCD, alb]);
exportImg(controls.toFloat().clip(bd).unmask(NODATA), 'controls_1km', 1000);

// ============================================================
// 9. Overpass (view) time of Terra and Aqua, national monthly mean (local solar time, hours)
// ============================================================
function viewTime(col, sensor) {
  return periods.map(function (p) {
    var s = ee.Date.fromYMD(p.y, p.m, 1);
    var img = col.filterDate(s, s.advance(1, 'month'))
      .select(['Day_view_time', 'Night_view_time']).mean().multiply(0.1);
    var st = img.reduceRegion({
      reducer: ee.Reducer.mean().combine({reducer2: ee.Reducer.stdDev(), sharedInputs: true}),
      geometry: bd, scale: 5000, maxPixels: 1e13, tileScale: 16
    });
    return ee.Feature(null, st).set({sensor: sensor, period: p.name});
  });
}
Export.table.toDrive({
  collection: ee.FeatureCollection(viewTime(terra, 'Terra').concat(viewTime(aqua, 'Aqua'))),
  description: 'overpass_time', folder: FOLDER, fileNamePrefix: 'overpass_time', fileFormat: 'CSV'
});

// ---------- quick check: smallest SO2 value in 2025 over Bangladesh ----------
// If this is >= -0.001, the -0.001 cut in the original script removed nothing
// (the L3 product already drops values below -0.001).
var so2raw = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_SO2')
  .filterDate(ee.Date.fromYMD(2025, 1, 1), y25e).filterBounds(bd).select('SO2_column_number_density');
print('Smallest SO2 value in 2025 (mol/m2):', so2raw.min().reduceRegion({
  reducer: ee.Reducer.min(), geometry: bd, scale: 5000, maxPixels: 1e13, tileScale: 16}));

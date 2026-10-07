// ============================================================
// Bangladesh Air Quality 2025 - FINAL EXPORT (Step 6)
// ============================================================
var bd = ee.FeatureCollection('projects/ee-2107046habibullah/assets/BangladeshBoundary').geometry();
var YEAR = 2025;
var start = ee.Date.fromYMD(YEAR, 1, 1);
var end = start.advance(1, 'year');
var FOLDER = 'BD_AQ_2025';
var NODATA = -9999;
var monthNums = [1,2,3,4,5,6,7,8,9,10,11,12];
var notWater = ee.ImageCollection('MODIS/061/MCD12Q1').sort('system:time_start', false).first().select('LC_Type1').neq(17);


// ---------- 1. Preparation functions (each returns one band named 'v') ----------
function prepS5P(id, band) {
  return ee.ImageCollection(id).filterDate(start, end).filterBounds(bd)
    .select([band], ['v']);
}

function aodPrep(img) {
  var a = img.select('Optical_Depth_055').multiply(0.001);
  return a.updateMask(a.gt(0)).updateMask(notWater).rename('v')
    .set('system:time_start', img.get('system:time_start'));
}

function lstPrep(band, qcBand) {
  return function(img) {
    var t = img.select(band).multiply(0.02).subtract(273.15);
    var qc = img.select(qcBand).unmask(0);
    var ok = qc.bitwiseAnd(3).lte(1).and(qc.rightShift(6).bitwiseAnd(3).lte(1));
    return t.updateMask(ok).rename('v')
      .set('system:time_start', img.get('system:time_start'));
  };
}

var modLST = ee.ImageCollection('MODIS/061/MOD11A1').filterDate(start, end).filterBounds(bd);

var vars = {
  CH4:  {col: prepS5P('COPERNICUS/S5P/OFFL/L3_CH4', 'CH4_column_volume_mixing_ratio_dry_air'), stat: 'mean', minObs: 1, scale: 1113.2},
  HCHO: {col: prepS5P('COPERNICUS/S5P/OFFL/L3_HCHO', 'tropospheric_HCHO_column_number_density'), stat: 'mean', minObs: 1, scale: 1113.2},
  O3:   {col: prepS5P('COPERNICUS/S5P/OFFL/L3_O3', 'O3_column_number_density'), stat: 'mean', minObs: 1, scale: 1113.2},
  UVAI: {col: prepS5P('COPERNICUS/S5P/OFFL/L3_AER_AI', 'absorbing_aerosol_index'), stat: 'mean', minObs: 1, scale: 1113.2},
  AOD:  {col: ee.ImageCollection('MODIS/061/MCD19A2_GRANULES').filterDate(start, end).filterBounds(bd).map(aodPrep), stat: 'median', minObs: 2, scale: 1000},
  LST_Day:   {col: modLST.map(lstPrep('LST_Day_1km', 'QC_Day')), stat: 'mean', minObs: 1, scale: 1000},
  LST_Night: {col: modLST.map(lstPrep('LST_Night_1km', 'QC_Night')), stat: 'mean', minObs: 1, scale: 1000}
};

// ---------- 2. Helpers ----------
function composite(v, s, e) {
  var col = v.col.filterDate(s, e);
  var img = (v.stat === 'median') ? col.median() : col.mean();
  return img.updateMask(col.count().gte(v.minObs)).rename('v').clip(bd);
}

function exportImg(img, name, scale) {
  Export.image.toDrive({
    image: img, description: name, folder: FOLDER, fileNamePrefix: name,
    region: bd.bounds(), scale: scale, crs: 'EPSG:4326', maxPixels: 1e13
  });
}

// ---------- 3. Monthly (12-band) and annual rasters ----------
Object.keys(vars).forEach(function(k) {
  var v = vars[k];
  var bands = monthNums.map(function(m) {
    var s = ee.Date.fromYMD(YEAR, m, 1);
    return composite(v, s, s.advance(1, 'month')).rename('M' + (m < 10 ? '0' + m : m));
  });
  exportImg(ee.Image.cat(bands).toFloat().unmask(NODATA), k + '_monthly_' + YEAR, v.scale);
  exportImg(composite(v, start, end).rename('annual').toFloat().unmask(NODATA), k + '_annual_' + YEAR, v.scale);
});

// ---------- 4. UHI (annual, day LST vs rural LST within 10 km, water masked) ----------
var projLST = modLST.first().select('LST_Day_1km').projection();
var lstDay = vars.LST_Day.col.mean().setDefaultProjection(projLST);
var lc = ee.ImageCollection('MODIS/061/MCD12Q1').sort('system:time_start', false).first().select('LC_Type1');
var rural = lc.neq(13).and(lc.neq(17));
var ruralLST = lstDay.updateMask(rural);
var kernel = ee.Kernel.circle(10000, 'meters');
var ruralSum = ruralLST.unmask(0).reduceNeighborhood({reducer: ee.Reducer.sum(), kernel: kernel}).reproject(projLST);
var ruralCount = ruralLST.mask().unmask(0).reduceNeighborhood({reducer: ee.Reducer.sum(), kernel: kernel}).reproject(projLST);
var uhi = lstDay.subtract(ruralSum.divide(ruralCount))
  .updateMask(lc.neq(17)).clip(bd).rename('UHI');
exportImg(uhi.toFloat().unmask(NODATA), 'UHI_annual_' + YEAR, 1000);

// ---------- 5. Monthly coverage table (for the 50% rule) ----------
var cov = [];
Object.keys(vars).forEach(function(k) {
  monthNums.forEach(function(m) {
    var s = ee.Date.fromYMD(YEAR, m, 1);
    var pct = composite(vars[k], s, s.advance(1, 'month')).mask().gt(0).rename('c').unmask(0)
      .reduceRegion({reducer: ee.Reducer.mean(), geometry: bd, scale: 5000, maxPixels: 1e13, tileScale: 16})
      .getNumber('c').multiply(100);
    cov.push(ee.Feature(null, {variable: k, month: m, pct_covered: pct}));
  });
});
Export.table.toDrive({collection: ee.FeatureCollection(cov), description: 'coverage_monthly_' + YEAR, folder: FOLDER, fileFormat: 'CSV'});

// ---------- 6. Daily national-mean CSVs (for time series graphs) ----------
var nDays = end.difference(start, 'day');
Object.keys(vars).forEach(function(k) {
  var v = vars[k];
  var fc = ee.FeatureCollection(ee.List.sequence(0, nDays.subtract(1)).map(function(d) {
    var s = start.advance(d, 'day');
    var col = v.col.filterDate(s, s.advance(1, 'day'));
    var val = ee.Algorithms.If(col.size().gt(0),
      col.mean().reduceRegion({reducer: ee.Reducer.mean(), geometry: bd, scale: 5000, maxPixels: 1e13, tileScale: 16}).get('v'),
      null);
    return ee.Feature(null, {variable: k, date: s.format('YYYY-MM-dd'), mean_value: val});
  }));
  Export.table.toDrive({collection: fc, description: k + '_daily_' + YEAR, folder: FOLDER, fileFormat: 'CSV'});
});
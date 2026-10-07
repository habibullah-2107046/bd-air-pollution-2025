// ============================================================
// NO2_CO_SO2_GEE_script.js
// Adds NO2, CO and SO2 (Sentinel-5P TROPOMI OFFL L3) for Bangladesh, 2025
// Same settings as for_all_Pollutants_GEE_script.js:
//   Drive folder BD_AQ_2025, EPSG:4326, 1113 m, no-data = -9999
//
// Outputs (9 tasks):
//   NO2_monthly_2025.tif  (12 bands, M01..M12)   NO2_annual_2025.tif   NO2_daily_2025.csv
//   CO_monthly_2025.tif                          CO_annual_2025.tif    CO_daily_2025.csv
//   SO2_monthly_2025.tif                         SO2_annual_2025.tif   SO2_daily_2025.csv
//
// Units of all three: mol/m2 (converted later in Python for the maps)
// ============================================================

var bdFull = ee.FeatureCollection('projects/ee-2107046habibullah/assets/BangladeshBoundary').geometry();
var bd = bdFull.simplify(500);          // lighter geometry for the daily means only

var YEAR = 2025;
var FOLDER = 'BD_AQ_2025';
var SCALE = 1113;                       // raster export scale (same as CH4, HCHO, O3, UVAI)
var DAILY_SCALE = 5000;                 // same as the other daily CSVs
var NODATA = -9999;

var start = ee.Date.fromYMD(YEAR, 1, 1);
var end = start.advance(1, 'year');

// minValid: values below this are dropped (only SO2 needs it; the product
// documentation recommends removing values below -0.001 mol/m2 as noise)
var VARS = [
  {name: 'NO2', id: 'COPERNICUS/S5P/OFFL/L3_NO2', band: 'tropospheric_NO2_column_number_density', minValid: null,
   vis: {min: 0, max: 0.00015, palette: ['blue', 'cyan', 'yellow', 'red']}},
  {name: 'CO',  id: 'COPERNICUS/S5P/OFFL/L3_CO',  band: 'CO_column_number_density', minValid: null,
   vis: {min: 0.03, max: 0.05, palette: ['blue', 'cyan', 'yellow', 'red']}},
  {name: 'SO2', id: 'COPERNICUS/S5P/OFFL/L3_SO2', band: 'SO2_column_number_density', minValid: -0.001,
   vis: {min: 0, max: 0.0005, palette: ['blue', 'cyan', 'yellow', 'red']}}
];

Map.centerObject(bdFull, 7);

VARS.forEach(function (v) {

  // ---- collection: one band called 'v'
  var col = ee.ImageCollection(v.id)
    .filterDate(start, end)
    .filterBounds(bdFull)
    .select(v.band)
    .map(function (img) {
      var x = img.rename('v');
      if (v.minValid !== null) { x = x.updateMask(x.gte(v.minValid)); }
      return x.set('system:time_start', img.get('system:time_start'));
    });

  print(v.name + ' images in ' + YEAR + ':', col.size());

  // ---- annual mean
  var annual = col.mean().rename(v.name).clip(bdFull);
  Map.addLayer(annual, v.vis, v.name + ' annual mean', v.name === 'NO2');

  Export.image.toDrive({
    image: annual.unmask(NODATA).toFloat(),
    description: v.name + '_annual_' + YEAR,
    folder: FOLDER,
    fileNamePrefix: v.name + '_annual_' + YEAR,
    region: bdFull.bounds(),
    scale: SCALE,
    crs: 'EPSG:4326',
    maxPixels: 1e13,
    formatOptions: {noData: NODATA}
  });

  // ---- monthly means, 12 bands (M01..M12)
  var stack = null;
  for (var m = 1; m <= 12; m++) {
    var s = ee.Date.fromYMD(YEAR, m, 1);
    var band = 'M' + (m < 10 ? '0' + m : m);
    var im = col.filterDate(s, s.advance(1, 'month')).mean().rename(band);
    stack = (stack === null) ? im : stack.addBands(im);
  }

  Export.image.toDrive({
    image: stack.clip(bdFull).unmask(NODATA).toFloat(),
    description: v.name + '_monthly_' + YEAR,
    folder: FOLDER,
    fileNamePrefix: v.name + '_monthly_' + YEAR,
    region: bdFull.bounds(),
    scale: SCALE,
    crs: 'EPSG:4326',
    maxPixels: 1e13,
    formatOptions: {noData: NODATA}
  });

  // ---- daily national mean (CSV)
  var nDays = end.difference(start, 'day');
  var days = ee.List.sequence(0, nDays.subtract(1));

  var daily = ee.FeatureCollection(days.map(function (d) {
    var ds = start.advance(d, 'day');
    var c = col.filterDate(ds, ds.advance(1, 'day'));
    var n = c.size();
    var props = ee.Dictionary(ee.Algorithms.If(n.gt(0),
      (function () {
        var img = c.mean();
        // 'v' = value (masked), 'valid' = 1 where data exists, 0 elsewhere
        var both = img.addBands(img.mask().gt(0).unmask(0).rename('valid'));
        return both.reduceRegion({
          reducer: ee.Reducer.mean(), geometry: bd,
          scale: DAILY_SCALE, maxPixels: 1e13, tileScale: 4
        });
      })(),
      ee.Dictionary({v: null, valid: 0})));
    return ee.Feature(null, {
      variable: v.name,
      date: ds.format('YYYY-MM-dd'),
      mean_value: props.get('v'),
      pct_valid: ee.Number(props.get('valid', 0)).multiply(100),
      n_images: n
    });
  }));

  Export.table.toDrive({
    collection: daily,
    description: v.name + '_daily_' + YEAR,
    folder: FOLDER,
    fileNamePrefix: v.name + '_daily_' + YEAR,
    fileFormat: 'CSV',
    selectors: ['variable', 'date', 'mean_value', 'pct_valid', 'n_images']
  });
});

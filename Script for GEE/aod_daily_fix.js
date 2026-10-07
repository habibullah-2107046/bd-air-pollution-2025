// ============================================================
// Bangladesh AOD 2025 - DAILY national-mean CSV (timeout fix)
// Splits the year into 12 monthly export tasks.
// Output: AOD_daily_2025_M01.csv ... AOD_daily_2025_M12.csv
// ============================================================
var bdFull = ee.FeatureCollection('projects/ee-2107046habibullah/assets/BangladeshBoundary').geometry();
var bd = bdFull.simplify(500);          // lighter geometry for reduceRegion only
var YEAR = 2025;
var FOLDER = 'BD_AQ_2025';
var SCALE = 5000;                       // same as the other daily CSVs

// Same prep as the main script (scale factor + drop non-positive values)
function aodPrep(img) {
  var a = img.select('Optical_Depth_055').multiply(0.001);
  return a.updateMask(a.gt(0)).rename('v')
    .set('system:time_start', img.get('system:time_start'));
}

var aodAll = ee.ImageCollection('MODIS/061/MCD19A2_GRANULES')
  .filterDate(ee.Date.fromYMD(YEAR, 1, 1), ee.Date.fromYMD(YEAR + 1, 1, 1))
  .filterBounds(bd)
  .select('Optical_Depth_055');

function dailyFeature(s) {
  s = ee.Date(s);
  var col = aodAll.filterDate(s, s.advance(1, 'day')).map(aodPrep);
  var n = col.size();
  var props = ee.Algorithms.If(n.gt(0),
    (function () {
      var img = col.mean();
      // band 'v' = AOD (masked), band 'valid' = 1 where AOD exists, 0 elsewhere
      var both = img.addBands(img.mask().gt(0).unmask(0).rename('valid'));
      return both.reduceRegion({
        reducer: ee.Reducer.mean(), geometry: bd,
        scale: SCALE, maxPixels: 1e13, tileScale: 16
      });
    })(),
    ee.Dictionary({v: null, valid: 0}));
  props = ee.Dictionary(props);
  return ee.Feature(null, {
    variable: 'AOD',
    date: s.format('YYYY-MM-dd'),
    mean_value: props.get('v'),
    pct_valid: ee.Number(props.get('valid', 0)).multiply(100),
    n_granules: n
  });
}

for (var m = 1; m <= 12; m++) {
  var s = ee.Date.fromYMD(YEAR, m, 1);
  var e = s.advance(1, 'month');
  var nDays = e.difference(s, 'day');
  var dates = ee.List.sequence(0, nDays.subtract(1)).map(function (d) {
    return s.advance(d, 'day');
  });
  var tag = 'AOD_daily_' + YEAR + '_M' + (m < 10 ? '0' + m : m);
  Export.table.toDrive({
    collection: ee.FeatureCollection(dates.map(dailyFeature)),
    description: tag, folder: FOLDER, fileNamePrefix: tag, fileFormat: 'CSV',
    selectors: ['variable', 'date', 'mean_value', 'pct_valid', 'n_granules']
  });
}

// Optional quick test before running all 12 tasks
// (prints one week of January to the console):
// print(ee.FeatureCollection(ee.List.sequence(0, 6).map(function (d) {
//   return dailyFeature(ee.Date.fromYMD(YEAR, 1, 1).advance(d, 'day'));
// })));
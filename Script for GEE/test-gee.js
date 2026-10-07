// var bd = ee.FeatureCollection('projects/ee-2107046habibullah/assets/BangladeshBoundary').geometry();
// var ch4 = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_CH4')
//   .select('CH4_column_volume_mixing_ratio_dry_air');

// var cov = ee.FeatureCollection(ee.List.sequence(1, 12).map(function(m) {
//   var s = ee.Date.fromYMD(2025, m, 1);
//   var col = ch4.filterDate(s, s.advance(1, 'month')).filterBounds(bd);
//   var pct = ee.Algorithms.If(col.size().gt(0),
//     col.count().rename('n').unmask(0).gt(0).reduceRegion({
//       reducer: ee.Reducer.mean(), geometry: bd, scale: 5000, maxPixels: 1e13
//     }).getNumber('n').multiply(100),
//     0);
//   return ee.Feature(null, {month: m, scenes: col.size(), pct_covered: pct});
// }));
// print('CH4 coverage 2025', cov);



// -------------------------------------------------------------------------------
// Step 1b: Check coverage for all variables

// var bd = ee.FeatureCollection('projects/ee-2107046habibullah/assets/BangladeshBoundary').geometry();

// var sources = {
//   CH4:       ['COPERNICUS/S5P/OFFL/L3_CH4', 'CH4_column_volume_mixing_ratio_dry_air'],
//   HCHO:      ['COPERNICUS/S5P/OFFL/L3_HCHO', 'tropospheric_HCHO_column_number_density'],
//   O3:        ['COPERNICUS/S5P/OFFL/L3_O3', 'O3_column_number_density'],
//   UVAI:      ['COPERNICUS/S5P/OFFL/L3_AER_AI', 'absorbing_aerosol_index'],
//   AOD:       ['MODIS/061/MCD19A2_GRANULES', 'Optical_Depth_055'],
//   LST_Day:   ['MODIS/061/MOD11A1', 'LST_Day_1km'],
//   LST_Night: ['MODIS/061/MOD11A1', 'LST_Night_1km']
// };

// function coverage(id, band) {
//   var col0 = ee.ImageCollection(id).select(band);
//   return ee.List.sequence(1, 12).map(function(m) {
//     var s = ee.Date.fromYMD(2025, m, 1);
//     var col = col0.filterDate(s, s.advance(1, 'month')).filterBounds(bd);
//     return ee.Algorithms.If(col.size().gt(0),
//       col.count().rename('n').unmask(0).gt(0).reduceRegion({
//         reducer: ee.Reducer.mean(), geometry: bd, scale: 5000, maxPixels: 1e13
//       }).getNumber('n').multiply(100).round(),
//       0);
//   });
// }

// Object.keys(sources).forEach(function(k) {
//   print(k + ' % covered (Jan–Dec)', coverage(sources[k][0], sources[k][1]));
  
// });



// --------------------------------------------------------------------------
// Step 1c: Check AOD coverage one month at a time
  
// var bd = ee.FeatureCollection('projects/ee-2107046habibullah/assets/BangladeshBoundary').geometry();
// var aod = ee.ImageCollection('MODIS/061/MCD19A2_GRANULES')
//   .select('Optical_Depth_055');

// for (var m = 1; m <= 12; m++) {
//   var s = ee.Date.fromYMD(2025, m, 1);
//   var col = aod.filterDate(s, s.advance(1, 'month')).filterBounds(bd);
//   var pct = col.map(function(img) { return img.mask().gt(0); })
//     .max().rename('n').unmask(0)
//     .reduceRegion({
//       reducer: ee.Reducer.mean(), geometry: bd,
//       scale: 10000, maxPixels: 1e13, tileScale: 16
//     }).getNumber('n').multiply(100).round();
//   print('AOD month ' + m + ' % covered', pct);
// }

// --------------------------------------------------------------------------------------------------
// Step 3a: Test the AOD quality mask
// var bd = ee.FeatureCollection('projects/ee-2107046habibullah/assets/BangladeshBoundary').geometry();
// var col = ee.ImageCollection('MODIS/061/MCD19A2_GRANULES');

// function rawAOD(img) {
//   return img.select('Optical_Depth_055').multiply(0.001);
// }
// function qaAOD(img) {
//   var qa = img.select('AOD_QA');
//   var clear = qa.bitwiseAnd(7).eq(1);                 // bits 0-2: clear sky
//   var best  = qa.rightShift(8).bitwiseAnd(15).eq(0);  // bits 8-11: best quality
//   return img.select('Optical_Depth_055').multiply(0.001)
//     .updateMask(clear.and(best));
// }

// function report(m, fn, label) {
//   var s = ee.Date.fromYMD(2025, m, 1);
//   var img = col.filterDate(s, s.advance(1, 'month')).filterBounds(bd)
//     .map(fn).mean().rename('aod');
//   var cover = img.mask().gt(0).unmask(0).reduceRegion({
//     reducer: ee.Reducer.mean(), geometry: bd,
//     scale: 10000, maxPixels: 1e13, tileScale: 16
//   }).getNumber('aod').multiply(100).round();
//   var st = img.reduceRegion({
//     reducer: ee.Reducer.minMax().combine({
//       reducer2: ee.Reducer.percentile([50, 99]), sharedInputs: true
//     }),
//     geometry: bd, scale: 5000, maxPixels: 1e13, tileScale: 16
//   });
//   print('Month ' + m + ' ' + label + ' | % covered:', cover, 'AOD stats:', st);
// }

// [1, 7].forEach(function(m) {
//   report(m, rawAOD, 'RAW');
//   report(m, qaAOD, 'QA');
// });



// -------------------------------------------------------------------------------------------------
// Step 3a-2: Test a relaxed mask on May and July
// var bd = ee.FeatureCollection('projects/ee-2107046habibullah/assets/BangladeshBoundary').geometry();
// var col = ee.ImageCollection('MODIS/061/MCD19A2_GRANULES');

// function rawAOD(img) {
//   return img.select('Optical_Depth_055').multiply(0.001);
// }
// function clearAOD(img) {
//   var qa = img.select('AOD_QA');
//   var clear = qa.bitwiseAnd(7).eq(1);          // bits 0-2: clear sky only
//   var aod = img.select('Optical_Depth_055').multiply(0.001);
//   return aod.updateMask(clear.and(aod.gt(0)));
// }

// function report(m, fn, label) {
//   var s = ee.Date.fromYMD(2025, m, 1);
//   var img = col.filterDate(s, s.advance(1, 'month')).filterBounds(bd)
//     .map(fn).mean().rename('aod');
//   var cover = img.mask().gt(0).unmask(0).reduceRegion({
//     reducer: ee.Reducer.mean(), geometry: bd,
//     scale: 10000, maxPixels: 1e13, tileScale: 16
//   }).getNumber('aod').multiply(100).round();
//   var st = img.reduceRegion({
//     reducer: ee.Reducer.minMax().combine({
//       reducer2: ee.Reducer.percentile([50, 99]), sharedInputs: true
//     }),
//     geometry: bd, scale: 1000, maxPixels: 1e13, tileScale: 16
//   });
//   print('Month ' + m + ' ' + label + ' | % covered:', cover, 'AOD stats:', st);
// }

// [5, 7].forEach(function(m) {
//   report(m, rawAOD, 'RAW');
//   report(m, clearAOD, 'CLEAR');
// });



// -------------------------------------------------------------------------------------------------
// Step 3a-3: Test mean vs median for May and July
// var bd = ee.FeatureCollection('projects/ee-2107046habibullah/assets/BangladeshBoundary').geometry();
// var col = ee.ImageCollection('MODIS/061/MCD19A2_GRANULES');

// function validAOD(img) {
//   var aod = img.select('Optical_Depth_055').multiply(0.001);
//   return aod.updateMask(aod.gt(0));
// }

// function report(m, label) {
//   var s = ee.Date.fromYMD(2025, m, 1);
//   var monthly = col.filterDate(s, s.advance(1, 'month')).filterBounds(bd)
//     .map(validAOD);
//   var img = (label === 'MEAN' ? monthly.mean() : monthly.median()).rename('aod');
//   var cover = img.mask().gt(0).unmask(0).reduceRegion({
//     reducer: ee.Reducer.mean(), geometry: bd,
//     scale: 10000, maxPixels: 1e13, tileScale: 16
//   }).getNumber('aod').multiply(100).round();
//   var st = img.reduceRegion({
//     reducer: ee.Reducer.minMax().combine({
//       reducer2: ee.Reducer.percentile([50, 99, 99.9]), sharedInputs: true
//     }),
//     geometry: bd, scale: 1000, maxPixels: 1e13, tileScale: 16
//   });
//   print('Month ' + m + ' ' + label + ' | % covered:', cover, 'AOD stats:', st);
// }

// [5, 7].forEach(function(m) {
//   report(m, 'MEAN');
//   report(m, 'MEDIAN');
// });



// -------------------------------------------------------------------------------------------------
// Step 3a-4: Test a minimum-observation rule
// var bd = ee.FeatureCollection('projects/ee-2107046habibullah/assets/BangladeshBoundary').geometry();
// var col = ee.ImageCollection('MODIS/061/MCD19A2_GRANULES');

// function validAOD(img) {
//   var aod = img.select('Optical_Depth_055').multiply(0.001);
//   return aod.updateMask(aod.gt(0));
// }

// function report(m, minObs) {
//   var s = ee.Date.fromYMD(2025, m, 1);
//   var monthly = col.filterDate(s, s.advance(1, 'month')).filterBounds(bd)
//     .map(validAOD);
//   var n = monthly.count();
//   var img = monthly.median().updateMask(n.gte(minObs)).rename('aod');
//   var cover = img.mask().gt(0).unmask(0).reduceRegion({
//     reducer: ee.Reducer.mean(), geometry: bd,
//     scale: 10000, maxPixels: 1e13, tileScale: 16
//   }).getNumber('aod').multiply(100).round();
//   var st = img.reduceRegion({
//     reducer: ee.Reducer.minMax().combine({
//       reducer2: ee.Reducer.percentile([50, 99]), sharedInputs: true
//     }),
//     geometry: bd, scale: 1000, maxPixels: 1e13, tileScale: 16
//   });
//   print('Month ' + m + ' | min obs ' + minObs + ' | % covered:', cover, 'AOD stats:', st);
// }

// [5, 7].forEach(function(m) {
//   report(m, 3);
//   report(m, 5);
// });




// -------------------------------------------------------------------------------------------------
// Step 3b: Test the LST quality mask
// var bd = ee.FeatureCollection('projects/ee-2107046habibullah/assets/BangladeshBoundary').geometry();
// var col = ee.ImageCollection('MODIS/061/MOD11A1');

// function makeFn(band, qcBand, mode) {
//   return function(img) {
//     var lst = img.select(band).multiply(0.02).subtract(273.15);
//     var qc = img.select(qcBand).unmask(0);
//     var quality = qc.bitwiseAnd(3);            // bits 0-1
//     var err = qc.rightShift(6).bitwiseAnd(3);  // bits 6-7
//     var mask = ee.Image(1);
//     if (mode === 'GOOD')  mask = quality.eq(0);
//     if (mode === 'ERR2K') mask = quality.lte(1).and(err.lte(1));
//     return lst.updateMask(mask).rename('lst');
//   };
// }

// function report(m, name, band, qcBand, mode) {
//   var s = ee.Date.fromYMD(2025, m, 1);
//   var img = col.filterDate(s, s.advance(1, 'month')).filterBounds(bd)
//     .map(makeFn(band, qcBand, mode)).mean().rename('lst');
//   var cover = img.mask().gt(0).unmask(0).reduceRegion({
//     reducer: ee.Reducer.mean(), geometry: bd,
//     scale: 10000, maxPixels: 1e13, tileScale: 16
//   }).getNumber('lst').multiply(100).round();
//   var st = img.reduceRegion({
//     reducer: ee.Reducer.percentile([1, 50]).combine({
//       reducer2: ee.Reducer.max(), sharedInputs: true
//     }),
//     geometry: bd, scale: 1000, maxPixels: 1e13, tileScale: 16
//   });
//   print('M' + m + ' ' + name + ' ' + mode + ' | % covered:', cover, st);
// }

// [1, 7].forEach(function(m) {
//   ['RAW', 'GOOD', 'ERR2K'].forEach(function(mode) {
//     report(m, 'Day',   'LST_Day_1km',   'QC_Day',   mode);
//     report(m, 'Night', 'LST_Night_1km', 'QC_Night', mode);
//   });
// });



// -------------------------------------------------------------------------------------------------
// Step 4c: Check CAMS PM data in GEE
// var bd = ee.FeatureCollection('projects/ee-2107046habibullah/assets/BangladeshBoundary').geometry();
// var cams = ee.ImageCollection('ECMWF/CAMS/NRT');

// print('CAMS band names', cams.first().bandNames());

// var bands = ['particulate_matter_d_less_than_25_um_surface',
//             'particulate_matter_d_less_than_10_um_surface'];

// [1, 7].forEach(function(m) {
//   var s = ee.Date.fromYMD(2025, m, 1);
//   var col = cams.filterDate(s, s.advance(1, 'month')).select(bands);
//   var img = col.mean().multiply(1e9);   // kg/m3 -> µg/m3
//   var st = img.reduceRegion({
//     reducer: ee.Reducer.mean().combine({
//       reducer2: ee.Reducer.minMax(), sharedInputs: true
//     }),
//     geometry: bd, scale: 10000, maxPixels: 1e13
//   });
//   print('Month ' + m + ' | images:', col.size(), 'PM µg/m3:', st);
// });




// -------------------------------------------------------------------------------------------------
// Step 5a: Test the new UHI method
var bd = ee.FeatureCollection('projects/ee-2107046habibullah/assets/BangladeshBoundary').geometry();

var mod = ee.ImageCollection('MODIS/061/MOD11A1')
  .filterDate('2025-01-01', '2026-01-01').filterBounds(bd);
var proj = mod.first().select('LST_Day_1km').projection();

function qcDay(img) {
  var lst = img.select('LST_Day_1km').multiply(0.02).subtract(273.15);
  var qc = img.select('QC_Day').unmask(0);
  var ok = qc.bitwiseAnd(3).lte(1).and(qc.rightShift(6).bitwiseAnd(3).lte(1));
  return lst.updateMask(ok);
}
var lst = mod.map(qcDay).mean().setDefaultProjection(proj).rename('lst');

var lcCol = ee.ImageCollection('MODIS/061/MCD12Q1').sort('system:time_start', false);
var lc = lcCol.first().select('LC_Type1');
print('Land cover year:', ee.Date(lcCol.first().get('system:time_start')).format('YYYY'));
var urban = lc.eq(13);
var rural = lc.neq(13).and(lc.neq(17));   // not urban, not water

var ruralLST = lst.updateMask(rural);
var kernel = ee.Kernel.circle(10000, 'meters');
var ruralSum = ruralLST.unmask(0).reduceNeighborhood({
  reducer: ee.Reducer.sum(), kernel: kernel
}).reproject(proj);
var ruralCount = ruralLST.mask().unmask(0).reduceNeighborhood({
  reducer: ee.Reducer.sum(), kernel: kernel
}).reproject(proj);
var ruralMean = ruralSum.divide(ruralCount);

var uhi = lst.subtract(ruralMean).rename('uhi').clip(bd);

var stats = function(mask, label) {
  print(label, uhi.updateMask(mask).reduceRegion({
    reducer: ee.Reducer.mean().combine({reducer2: ee.Reducer.minMax(), sharedInputs: true}),
    geometry: bd, scale: 1000, maxPixels: 1e13, tileScale: 16
  }));
};
stats(urban, 'UHI over URBAN pixels (°C)');
stats(rural, 'UHI over RURAL pixels (°C)');

Map.centerObject(bd, 7);
Map.addLayer(uhi, {min: -3, max: 3, palette: ['blue', 'white', 'red']}, 'UHI 2025');


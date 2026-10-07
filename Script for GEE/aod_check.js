// AOD check: observation count + water-masked annual
var bd = ee.FeatureCollection('projects/ee-2107046habibullah/assets/BangladeshBoundary').geometry();
var start = ee.Date.fromYMD(2025, 1, 1), end = start.advance(1, 'year');

function aodPrep(img) {
  var a = img.select('Optical_Depth_055').multiply(0.001);
  return a.updateMask(a.gt(0)).rename('v');
}
var aod = ee.ImageCollection('MODIS/061/MCD19A2_GRANULES')
  .filterDate(start, end).filterBounds(bd).map(aodPrep);

// 1. How many observations each pixel has
var count = aod.count().clip(bd);

// 2. Annual median AOD with water removed
var lc = ee.ImageCollection('MODIS/061/MCD12Q1').sort('system:time_start', false).first().select('LC_Type1');
var aodAnnual = aod.median().updateMask(aod.count().gte(2)).updateMask(lc.neq(17)).clip(bd);

Map.centerObject(bd, 7);
Map.addLayer(count, {min: 50, max: 400, palette: ['red', 'yellow', 'green']}, 'AOD observation count');
Map.addLayer(aodAnnual, {min: 0.16, max: 0.68, palette: ['green', 'yellow', 'brown']}, 'AOD annual (water masked)', false);
// =====================================================================
// PAPER 1 · STEP 7 — DAILY RAIN FOR THE LOCAL CATCHMENT OF EVERY CROSSING  (GEE run C)
//
// Purpose: independent rain data, used in Python to
//   (1) confirm each satellite flood event (was there rain before it?),
//   (2) find dry periods at every crossing (Check F),
//   (3) build the rain-based flood episodes for small rivers (timing observability, Step 9).
// Big rivers get river flow from your GloFAS files (Python); rain is exported for all crossings anyway.
//
// Catchments: the HydroBASINS level-12 basin (local, about 100-300 km2) and level-10 basin (larger, about
// 1,000-3,000 km2) that contain each crossing. Rain products:
//   chirps_mm  CHIRPS daily (0.05 deg, satellite + rain gauges)              -> main rain product
//   imerg_mm   GPM IMERG V07 (0.1 deg); daily total = mean half-hourly rate x 24 h (days in UTC) -> second product
//              (in Earth Engine IMERG V07 is final up to Sep 2025; Oct-Dec 2025 are provisional values)
//
// Outputs: Google Drive folder P1_Step07_Rain (13 tasks)
//   P1_S07_Rain_<year>.csv (2014-2025; 2014 gives the rain before the first 2015 passes)
//       basin_id, level (12 or 10), date, chirps_mm, imerg_mm   (mean over the basin, mm per day)
//   P1_S07_Rain_basin_lookup.csv
//       cross_id, hybas12, area12_km2, up12_km2, hybas10, area10_km2, up10_km2
// How to run: NEW Code Editor script -> paste -> Run -> press RUN on the 13 tasks. Nothing is drawn.
// =====================================================================

// ---------- 1. Settings ----------
var CROSSINGS_ASSET = 'projects/ee-zabih/assets/CPEC_Crossings_FINAL';
var DRIVE_FOLDER    = 'P1_Step07_Rain';
var PREFIX          = 'P1_S07_Rain';
var YEARS = [2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];
var SCALE_M = 1000;                       // fine sampling grid, so small basins get an area-weighted mean

// ---------- 2. Crossings and their basins ----------
var crossings = ee.FeatureCollection(CROSSINGS_ASSET);
var crossPts  = crossings.geometry();
var B12 = ee.FeatureCollection('WWF/HydroSHEDS/v1/Basins/hybas_12').filterBounds(crossPts)
            .map(function (f) { return f.set('level', 12); });
var B10 = ee.FeatureCollection('WWF/HydroSHEDS/v1/Basins/hybas_10').filterBounds(crossPts)
            .map(function (f) { return f.set('level', 10); });
var BASINS = B12.merge(B10);
print('Crossings:', crossings.size(), '| level-12 basins:', B12.size(), '| level-10 basins:', B10.size());

function basinValue(fc, g, prop) {
  var hit = fc.filterBounds(g);
  return ee.Algorithms.If(hit.size().gt(0), ee.Feature(hit.first()).get(prop), -1);
}
var lookup = crossings.map(function (f) {
  var g = f.geometry();
  return ee.Feature(null, {
    cross_id: f.get('cross_id'),
    hybas12: basinValue(B12, g, 'HYBAS_ID'), area12_km2: basinValue(B12, g, 'SUB_AREA'), up12_km2: basinValue(B12, g, 'UP_AREA'),
    hybas10: basinValue(B10, g, 'HYBAS_ID'), area10_km2: basinValue(B10, g, 'SUB_AREA'), up10_km2: basinValue(B10, g, 'UP_AREA')
  });
});
Export.table.toDrive({
  collection: lookup, description: PREFIX + '_basin_lookup', folder: DRIVE_FOLDER, fileFormat: 'CSV',
  selectors: ['cross_id', 'hybas12', 'area12_km2', 'up12_km2', 'hybas10', 'area10_km2', 'up10_km2']
});

// ---------- 3. One image per day: CHIRPS and IMERG (mm/day) ----------
var CHIRPS = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY').select('precipitation');
var IMERG  = ee.ImageCollection('NASA/GPM_L3/IMERG_V07').select('precipitation');      // mm/hr, every 30 min
var EMPTY  = ee.Image.constant(0).rename('precipitation').toFloat().updateMask(ee.Image(0));

function dayImage(d) {
  d = ee.Date(d);
  var next = d.advance(1, 'day');
  var c = CHIRPS.filterDate(d, next).merge(ee.ImageCollection([EMPTY])).max();          // masked if the day is missing
  var i = IMERG.filterDate(d, next).merge(ee.ImageCollection([EMPTY])).mean().multiply(24);
  return c.rename('chirps').addBands(i.rename('imerg')).set('date', d.format('YYYY-MM-dd'));
}

function round2(v) {                                                    // 2 decimals; empty stays empty
  return ee.Algorithms.If(ee.Algorithms.IsEqual(v, null), null, ee.Number(v).multiply(100).round().divide(100));
}

// ---------- 4. Exports: one task per year ----------
YEARS.forEach(function (y) {
  var start = ee.Date.fromYMD(y, 1, 1);
  var nDays = ee.Date.fromYMD(y + 1, 1, 1).difference(start, 'day');
  var days  = ee.List.sequence(0, nDays.subtract(1)).map(function (k) { return start.advance(k, 'day'); });
  var rows = ee.FeatureCollection(days.map(function (d) {
    var img = dayImage(d);
    return img.reduceRegions({collection: BASINS, reducer: ee.Reducer.mean(), scale: SCALE_M, tileScale: 2})
      .map(function (f) {
        return ee.Feature(null, {basin_id: f.get('HYBAS_ID'), level: f.get('level'), date: img.get('date'),
                                 chirps_mm: round2(f.get('chirps')), imerg_mm: round2(f.get('imerg'))});
      });
  })).flatten();
  Export.table.toDrive({
    collection: rows, description: PREFIX + '_' + y, folder: DRIVE_FOLDER, fileFormat: 'CSV',
    selectors: ['basin_id', 'level', 'date', 'chirps_mm', 'imerg_mm']
  });
});

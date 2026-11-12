// =====================================================================
// PAPER 1 · STEP 2a — Describe the FINAL units and crossings (run once, after the fixed notebook)
// Google Earth Engine Code Editor version of kit script gee_00_study_area.py
//
// What it makes (Google Drive folder DRIVE_FOLDER):
//   P1_2a_units_attributes.csv       one row per HydroBASINS L6 unit
//   P1_2a_crossings_attributes.csv   one row per road–river crossing
//   P1_2a_aridity_index_1km.tif      aridity index raster (P / PET, 1991–2020)
// and one new asset with a fixed crossing ID used by every later step:
//   CROSSINGS_OUT (below)
//
// How to run: paste into the Code Editor → Run → Tasks tab → press RUN on the 4 tasks.
// Nothing heavy is drawn on the map, so "capacity exceeded" should not appear;
// all heavy work happens inside the export tasks.
// =====================================================================

// ---------- 1. Settings (edit only this block) ----------
var UNITS_ASSET     = 'projects/ee-zabih/assets/CPEC_HydroBasins_L6';
var CROSSINGS_ASSET = 'projects/ee-zabih/assets/CPEC_Road_River_Crossings';
var BUFFER_ASSET    = 'projects/ee-zabih/assets/CPEC_10km_Corridor_Buffer';
var CROSSINGS_OUT   = 'projects/ee-zabih/assets/CPEC_Crossings_FINAL';  // permanent ID list
var UNIT_ID_FIELD   = 'HYBAS_ID';    // ID field in your L6 file
var DRIVE_FOLDER    = 'CPEC_Paper1';
var UPA_SEARCH_M    = 1000;          // search radius for the MERIT pixel whose upstream area best matches HydroRIVERS
var SHOW_MAP        = false;         // true = draw the three vector layers only

// ---------- 2. Inputs ----------
var units        = ee.FeatureCollection(UNITS_ASSET);
var crossingsRaw = ee.FeatureCollection(CROSSINGS_ASSET);
var buffer       = ee.FeatureCollection(BUFFER_ASSET);

print('Number of L6 units:', units.size());
print('Number of crossings:', crossingsRaw.size());
print('Fields in units file:', units.first().propertyNames());
print('Fields in crossings file:', crossingsRaw.first().propertyNames());

// ---------- 3. Aridity index and dryness class (UNEP thresholds) ----------
// AI = mean annual precipitation / mean annual PET, TerraClimate 1991–2020.
// TerraClimate 'pet' has scale factor 0.1; 'pr' has none.
var tc  = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE').filterDate('1991-01-01', '2021-01-01');
var nYears = 30;
var P   = tc.select('pr').sum().divide(nYears);
var PET = tc.select('pet').sum().multiply(0.1).divide(nYears);
// Capped at 5: on glaciers TerraClimate PET is ~0, which gives meaningless values up to ~1000.
// Every value above 0.65 is 'humid' anyway, so the cap does not change any class.
var AI  = P.divide(PET.max(1)).min(5).rename('AI');

// 1 hyper-arid <0.05 | 2 arid 0.05–0.20 | 3 semi-arid 0.20–0.50 | 4 dry sub-humid 0.50–0.65 | 5 humid >=0.65
var aridClass = ee.Image(0)
  .where(AI.lt(0.05), 1)
  .where(AI.gte(0.05).and(AI.lt(0.20)), 2)
  .where(AI.gte(0.20).and(AI.lt(0.50)), 3)
  .where(AI.gte(0.50).and(AI.lt(0.65)), 4)
  .where(AI.gte(0.65), 5)
  .updateMask(AI.mask())
  .rename('arid_class').toByte();

// ---------- 4. Terrain, hydrography, land cover ----------
var demCol = ee.ImageCollection('COPERNICUS/DEM/GLO30').select('DEM');
var dem    = demCol.mosaic().setDefaultProjection(demCol.first().projection()).rename('elev');
var slope  = ee.Terrain.slope(dem).rename('slope');
var merit  = ee.Image('MERIT/Hydro/v1_0_1');
var hand   = merit.select('hnd').rename('hand');     // height above nearest drainage (m)
var upa    = merit.select('upa').rename('upa');      // upstream drainage area (km2)
var wc     = ee.ImageCollection('ESA/WorldCover/v200').first().select('Map');

// ---------- 5. Unit attributes ----------
var unitStack = AI
  .addBands(dem).addBands(slope).addBands(hand)
  .addBands(wc.eq(40).rename('crop_share'))
  .addBands(wc.eq(60).rename('bare_share'))
  .addBands(wc.eq(50).rename('built_share'))
  .addBands(wc.eq(10).rename('tree_share'));

var aridShares = ee.Image.cat([1, 2, 3, 4, 5].map(function (k) {
  return aridClass.eq(k).rename('arid_sh' + k);
}));

var unitsMean = unitStack.reduceRegions({
  collection: units, reducer: ee.Reducer.mean(), scale: 250, tileScale: 4
});
var unitsAll = aridShares.reduceRegions({
  collection: unitsMean, reducer: ee.Reducer.mean(), scale: 1000, tileScale: 4
});

var unitsFinal = unitsAll.map(function (f) {
  var sh = ee.List([1, 2, 3, 4, 5].map(function (k) {
    var v = f.get('arid_sh' + k);
    return ee.Number(ee.Algorithms.If(v, v, 0));
  }));
  var maxSh = ee.Number(sh.reduce(ee.Reducer.max()));
  var cls = ee.Number(ee.Algorithms.If(maxSh.gt(0), ee.Number(sh.indexOf(maxSh)).add(1), 0));
  var c = f.geometry().centroid(100).coordinates();
  return ee.Feature(null, f.toDictionary()).set({
    unit_id_str: ee.Number(f.get(UNIT_ID_FIELD)).toLong().format('%d'),
    arid_class: cls,                 // majority dryness class (0 = no data)
    dryland: cls.gte(1).and(cls.lte(3)),
    area_km2: f.geometry().area(100).divide(1e6),
    lon: c.get(0),
    lat: c.get(1)
  });
});

// ---------- 6. Crossings: fixed ID, unit, point attributes ----------
// Keep an existing 'cross_id'; otherwise number the crossings 1..N.
var hasId = ee.List(crossingsRaw.first().propertyNames()).contains('cross_id');
var nC = crossingsRaw.size();
var cList = crossingsRaw.toList(nC);
var withNewId = ee.FeatureCollection(ee.List.sequence(0, nC.subtract(1)).map(function (i) {
  return ee.Feature(cList.get(i)).set('cross_id', ee.Number(i).add(1));
}));
var crossings = ee.FeatureCollection(ee.Algorithms.If(hasId, crossingsRaw, withNewId));

// Which L6 unit contains each crossing
var joined = ee.Join.saveFirst({matchKey: 'unit', outer: true}).apply({
  primary: crossings,
  secondary: units,
  condition: ee.Filter.intersects({leftField: '.geo', rightField: '.geo', maxError: 10})
});
print('Crossings inside an L6 unit:',
      ee.FeatureCollection(joined).filter(ee.Filter.notNull(['unit'])).size());

var pointStack = AI.addBands(aridClass).addBands(dem).addBands(slope).addBands(hand)
  .addBands(wc.rename('landcover'));

var crossAttr = ee.FeatureCollection(joined).map(function (f) {
  var u = f.get('unit');
  var pt = f.geometry().centroid(1);
  var v = pointStack.reduceRegion({reducer: ee.Reducer.first(), geometry: pt, scale: 30});
  // MERIT pixel within UPA_SEARCH_M whose upstream area is closest to HydroRIVERS UPLAND_SKM
  var target = ee.Number(f.get('UPLAND_SKM'));
  var score = upa.divide(target).log().abs().rename('score');
  var best = score.addBands(upa).reduceRegion({
    reducer: ee.Reducer.min(2).setOutputs(['score_min', 'upa_best']),
    geometry: pt.buffer(UPA_SEARCH_M), scale: 90
  });
  var up = best.get('upa_best');
  var c = pt.coordinates();
  var props = f.toDictionary().remove(['unit'], true);
  return ee.Feature(pt, props).set(v).set({
    unit_id: ee.Algorithms.If(u, ee.Feature(u).get(UNIT_ID_FIELD), -1),
    unit_id_str: ee.Algorithms.If(u, ee.Number(ee.Feature(u).get(UNIT_ID_FIELD)).toLong().format('%d'), 'none'),
    upa_merit_km2: up,               // best-matching MERIT upstream area (km2), for GloFAS matching
    upa_ratio: ee.Algorithms.If(up, ee.Number(up).divide(target), null),   // ~1 = crossing sits on the right river
    lon: c.get(0),
    lat: c.get(1)
  });
});

// Preview one record of each (small, safe to compute interactively)
print('Example unit record:', unitsFinal.first());
print('Example crossing record:', crossAttr.first());

// ---------- 7. Exports (run them in the Tasks tab) ----------
Export.table.toAsset({
  collection: crossAttr,
  description: 'P1_2a_crossings_with_ID_asset',
  assetId: CROSSINGS_OUT
});

Export.table.toDrive({
  collection: crossAttr.map(function (f) { return ee.Feature(null, f.toDictionary()); }),
  description: 'P1_2a_crossings_attributes',
  folder: DRIVE_FOLDER,
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: unitsFinal,
  description: 'P1_2a_units_attributes',
  folder: DRIVE_FOLDER,
  fileFormat: 'CSV'
});

Export.image.toDrive({
  image: AI.toFloat(),
  description: 'P1_2a_aridity_index_1km',
  folder: DRIVE_FOLDER,
  region: units.geometry().bounds(1000),
  scale: 1000,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// ---------- 8. Optional map (vectors only) ----------
if (SHOW_MAP) {
  Map.centerObject(buffer, 5);
  Map.addLayer(units.style({color: '6b7780', fillColor: '00000000', width: 1}), {}, 'L6 units');
  Map.addLayer(buffer.style({color: '0f5e66', fillColor: '0f5e6622', width: 1}), {}, '10 km buffer');
  Map.addLayer(crossingsRaw.style({color: 'b9853f', pointSize: 2}), {}, 'Crossings');

// Add these to sanity-check the computed layers:
  Map.addLayer(AI.clip(buffer), {min: 0, max: 1, palette: ['brown', 'yellow', 'green', 'blue']}, 'Aridity Index', false);
  Map.addLayer(aridClass.clip(buffer), {min: 1, max: 5, palette: ['8B4513','D2B48C','FFD700','9ACD32','006400']}, 'Dryness Class', false);
}

// =====================================================================
// PAPER 1 · STEP 2 — Dryland-aware Sentinel-1 flood classifier, event E2022 (15 Jun – 15 Oct 2022)
// Code Editor version of kit script gee_02_flood_mapping.py (same rules, same classes)
//
// Class of every pixel on every satellite pass:
//   0 no observation | 1 not flooded | 2 flooded, high confidence | 3 flooded, low confidence
//   4 permanent water | 5 unobservable (cannot be judged; NEVER counted as dry)
//
// Run it three times; change ONLY the MODE line:
//   MODE = 'test'     → 5 tasks: (a) map of a test area (Dadu–Larkana–Jacobabad–Naseerabad, arid, 35 crossings)
//                                (b) every 2022 satellite scene at ALL crossings, one task per month (4 tasks)
//   MODE = 'corridor' → 1 task : event map for the whole 10 km corridor, saved as a GEE asset
//   MODE = 'tables'   → 2 tasks: flooded area per L6 unit and per dryness class (reads that asset)
//   MODE = 'passes'   → 2 tasks: classifier AND naive result for each pass 24 Aug – 3 Sep 2022 in the test
//                                area (to compare with UNOSAT maps of 27, 29 and 31 Aug), plus the pass list
// Nothing is drawn on the map, so "capacity exceeded" cannot appear; all work runs in the Tasks tab.
// =====================================================================

// ---------- 1. Settings ----------
var MODE = 'corridor';                 // 'test' → 'corridor' → 'tables'   ('passes' = UNOSAT comparison)

var CROSSINGS_ASSET = 'projects/ee-zabih/assets/CPEC_Crossings_FINAL';
var UNITS_ASSET     = 'projects/ee-zabih/assets/CPEC_HydroBasins_L6';
var BUFFER_ASSET    = 'projects/ee-zabih/assets/CPEC_10km_Corridor_Buffer';
var SUMMARY_ASSET   = 'projects/ee-zabih/assets/P1_flood_E2022_summary';   // written in 'corridor', read in 'tables'
var DRIVE_FOLDER    = 'CPEC_Paper1';
var TEST_BOX = ee.Geometry.Rectangle([67.3, 26.3, 68.5, 28.6]);

var EVENT_ID = 'E2022';
var EVENT    = ['2022-06-15', '2022-10-15'];                               // end date not included
var REF      = [['2020-06-15', '2020-10-15'], ['2021-06-15', '2021-10-15']]; // same season, earlier years
var LONG_REF = ['2021-01-01', '2022-01-01'];                               // 12 months, for the look-alike layer
var MONTHS   = [['2022-06-15', '2022-07-15'], ['2022-07-15', '2022-08-15'],
                ['2022-08-15', '2022-09-15'], ['2022-09-15', '2022-10-15']];  // crossing table: one task each
var CROSS_BUFFER_M = 300;                                                  // circle sampled around each crossing
var PASS_WINDOW = ['2022-08-24', '2022-09-04'];                            // MODE 'passes' (end date not included)

// Rules (same values as the kit's config.FLOOD; do not change without logging it)
var F = {
  speckle_radius_m: 30,       // focal median speckle filter
  border_noise_db: -30,       // VV below this = border noise / no data
  z_k: 2.5,                   // flood if VV dropped by >= 2.5 reference standard deviations ...
  z_k_vh: 2.0,                // ... and VH by >= 2.0 ...
  abs_vv_db: -15,             // ... and VV itself is <= -15 dB
  z_k_strict: 3.5,            // stricter drop needed on look-alike (always-dark) surfaces such as sand
  min_std_db: 0.5,            // floor on the reference standard deviation
  min_ref_count: 3,           // fewer reference passes -> unobservable
  lookalike_vv_db: -18,       // 'dark' threshold for the look-alike layer ...
  lookalike_freq: 0.5,        // ... dark in >= 50% of the 12-month reference = look-alike surface
  hand_max_m: 15,             // higher above the nearest drainage -> implausible (low confidence, not erased)
  slope_plaus_deg: 5,         // steeper -> implausible (low confidence)
  slope_terrain_deg: 15,      // steeper -> radar layover/shadow risk -> unobservable unless a clear drop
  perm_water_occ: 80,         // JRC water occurrence >= 80% -> permanent water
  baseline_static_vv_db: -18  // naive single-threshold method, kept only for comparison
};

// ---------- 2. Inputs ----------
var crossings  = ee.FeatureCollection(CROSSINGS_ASSET);
var units      = ee.FeatureCollection(UNITS_ASSET);
var bufferGeom = ee.FeatureCollection(BUFFER_ASSET).geometry();
print('Mode:', MODE, '| crossings:', crossings.size(), '| L6 units:', units.size());

// ---------- 3. Static layers ----------
var demCol = ee.ImageCollection('COPERNICUS/DEM/GLO30').select('DEM');
var dem    = demCol.mosaic().setDefaultProjection(demCol.first().projection());
var slope  = ee.Terrain.slope(dem);
var hand   = ee.Image('MERIT/Hydro/v1_0_1').select('hnd');
var wc     = ee.ImageCollection('ESA/WorldCover/v200').first().select('Map');
var perm        = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence').unmask(0).gte(F.perm_water_occ);
var implausible = hand.gt(F.hand_max_m).or(slope.gt(F.slope_plaus_deg)).unmask(0);
var terrainBad  = slope.gt(F.slope_terrain_deg).unmask(0);
var hardSurface = wc.eq(50).or(wc.eq(10)).unmask(0);           // built-up, tree cover

// Aridity class (identical to Step 2a), used in MODE 'tables'
var tc  = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE').filterDate('1991-01-01', '2021-01-01');
var P   = tc.select('pr').sum().divide(30);
var PET = tc.select('pet').sum().multiply(0.1).divide(30);
var AI  = P.divide(PET.max(1)).min(5);
var aridClass = ee.Image(0)
  .where(AI.lt(0.05), 1)
  .where(AI.gte(0.05).and(AI.lt(0.20)), 2)
  .where(AI.gte(0.20).and(AI.lt(0.50)), 3)
  .where(AI.gte(0.50).and(AI.lt(0.65)), 4)
  .where(AI.gte(0.65), 5)
  .rename('arid_class').toByte();

// ---------- 4. Sentinel-1 preparation ----------
var PROPS = ['system:time_start', 'relativeOrbitNumber_start', 'orbitProperties_pass', 'platform_number'];

function s1Base(region, start, end) {
  return ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(region)
    .filterDate(start, end)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
    .map(function (img) {
      var sig = img.select(['VV', 'VH']);
      sig = sig.updateMask(sig.select('VV').gt(F.border_noise_db));
      sig = sig.focalMedian({radius: F.speckle_radius_m, kernelType: 'circle', units: 'meters'});
      var day = img.date().format('YYYY-MM-dd');
      var ro  = ee.Number(img.get('relativeOrbitNumber_start')).format('%d');
      return ee.Image(sig.copyProperties(img, PROPS)).set({day: day, key: ee.String(day).cat('_').cat(ro)});
    });
}

// scene count only (metadata, no pixels): safe to print
function rawScenes(region) {
  return ee.ImageCollection('COPERNICUS/S1_GRD').filterBounds(region).filterDate(EVENT[0], EVENT[1])
    .filter(ee.Filter.eq('instrumentMode', 'IW'));
}

function s1Windows(region, windows) {
  var out = s1Base(region, windows[0][0], windows[0][1]);
  for (var w = 1; w < windows.length; w++) {
    out = out.merge(s1Base(region, windows[w][0], windows[w][1]));
  }
  return out;
}

// one image per (day, relative orbit): frames along the same pass are mosaicked (used for the maps)
function mosaicByPass(col) {
  var keys = col.aggregate_array('key').distinct();
  return ee.ImageCollection.fromImages(keys.map(function (k) {
    var sub = col.filter(ee.Filter.eq('key', k));
    var first = ee.Image(sub.first());
    return ee.Image(sub.mosaic().copyProperties(first, PROPS.concat(['day', 'key'])));
  }));
}

// reference passes from the same relative orbit and direction; a masked placeholder makes an empty match
// give count 0 (-> unobservable) instead of an image without bands
var EMPTY = ee.Image.constant([0, 0]).rename(['VV', 'VH']).toFloat().updateMask(ee.Image(0));
function sameGeometry(col, img) {
  return col.filter(ee.Filter.eq('relativeOrbitNumber_start', img.get('relativeOrbitNumber_start')))
            .filter(ee.Filter.eq('orbitProperties_pass', img.get('orbitProperties_pass')))
            .merge(ee.ImageCollection([EMPTY]));
}

// ---------- 5. Classifier (one pass) ----------
function classify(img, refCol, longCol) {
  var ref = sameGeometry(refCol, img);
  var lng = sameGeometry(longCol, img);
  var stats = ref.select(['VV', 'VH']).reduce(
    ee.Reducer.median().combine({reducer2: ee.Reducer.stdDev(), sharedInputs: true})
                       .combine({reducer2: ee.Reducer.count(), sharedInputs: true}));
  var zVV = img.select('VV').subtract(stats.select('VV_median')).divide(stats.select('VV_stdDev').max(F.min_std_db));
  var zVH = img.select('VH').subtract(stats.select('VH_median')).divide(stats.select('VH_stdDev').max(F.min_std_db));

  var lookalike = lng.map(function (i) { return i.select('VV').lt(F.lookalike_vv_db); })
                     .mean().unmask(0).gte(F.lookalike_freq);
  var refOk  = stats.select('VV_count').unmask(0).gte(F.min_ref_count);
  var dark   = img.select('VV').lte(F.abs_vv_db);
  var change = zVV.lte(-F.z_k).and(dark).and(zVH.lte(-F.z_k_vh));
  var strict = zVV.lte(-F.z_k_strict).and(dark);
  var valid  = img.select('VV').mask().gt(0);

  var cls = ee.Image(0)
    .where(valid, 1)
    .where(valid.and(refOk.not()), 5)                                     // no usable reference
    .where(valid.and(refOk).and(terrainBad).and(change.not()), 5)         // steep terrain without a signal
    .where(valid.and(refOk).and(lookalike).and(strict.not()), 5);         // always-dark surface: not dry, unknown
  var lowConf = lookalike.or(implausible).or(hardSurface).or(terrainBad);
  var flood   = valid.and(refOk).and(change).and(lookalike.not().or(strict));
  cls = cls.where(flood.and(lowConf.not()), 2)
           .where(flood.and(lowConf), 3)
           .where(valid.and(perm), 4);

  var baseline = valid.and(img.select('VV').lt(F.baseline_static_vv_db)).and(perm.not());
  return ee.Image(cls.rename('cls').toByte().addBands(baseline.rename('baseline').toByte())
    .copyProperties(img, PROPS.concat(['day', 'key'])));
}

function classified(region) {
  var ev  = mosaicByPass(s1Base(region, EVENT[0], EVENT[1]));
  var ref = s1Windows(region, REF);
  var lng = s1Base(region, LONG_REF[0], LONG_REF[1]);
  return ev.map(function (img) { return classify(img, ref, lng); });
}

// ---------- 6. Event summary (all passes together) ----------
function eventSummary(clsCol) {
  var c = clsCol.select('cls');
  function count(test) { return c.map(test).sum(); }
  var nValid = count(function (i) { return i.gte(1).and(i.lte(4)); });
  var nHigh  = count(function (i) { return i.eq(2); });
  var nAny   = count(function (i) { return i.eq(2).or(i.eq(3)); });
  var nNot   = count(function (i) { return i.eq(1); });
  var nUnobs = count(function (i) { return i.eq(5); });
  var nPerm  = count(function (i) { return i.eq(4); });
  var t0 = ee.Date(EVENT[0]).millis();
  var days = c.map(function (i) {
    var d = ee.Number(i.get('system:time_start')).subtract(t0).divide(86400000);
    return ee.Image.constant(d).toFloat().updateMask(i.eq(2).or(i.eq(3)));
  });
  var ev = ee.Image(0)
    .where(nUnobs.gt(0), 5)
    .where(nNot.gt(0), 1)
    .where(nPerm.gt(0), 4)
    .where(nAny.gt(0), 3)
    .where(nHigh.gt(0), 2);
  return ev.rename('event_class').toByte()
    .addBands(nValid.rename('n_valid').toByte())
    .addBands(nHigh.rename('n_high').toByte())
    .addBands(nAny.rename('n_any').toByte())
    .addBands(nUnobs.rename('n_unobs').toByte())
    .addBands(days.min().floor().unmask(255).rename('first_day').toByte())   // days after 15 Jun 2022; 255 = never
    .addBands(days.max().floor().unmask(255).rename('last_day').toByte())
    .addBands(clsCol.select('baseline').max().unmask(0).rename('baseline_flood').toByte());
}

// ---------- 7. Crossing table: every scene, class shares inside 300 m of each crossing ----------
// Scenes are classified one by one (no mosaics), so the memory stays small. A crossing that lies in the
// small overlap of two frames of the same pass gets two rows with the same day and orbit; Python keeps one.
var crossBuf = crossings.map(function (f) {
  return ee.Feature(f.geometry().buffer(CROSS_BUFFER_M, 5), {cross_id: f.get('cross_id')});
});
var ROW_COLS = ['cross_id', 'scene', 'day', 'time_utc', 'rel_orbit', 'pass', 'platform',
                'f_c0', 'f_c1', 'f_c2', 'f_c3', 'f_c4', 'f_c5', 'f_base'];
var crossPts = crossings.geometry();
var crossRef = s1Windows(crossPts, REF);
var crossLng = s1Base(crossPts, LONG_REF[0], LONG_REF[1]);

function crossingRows(scene) {
  var img = classify(scene, crossRef, crossLng);
  var c = img.select('cls');
  var shares = ee.Image.cat([0, 1, 2, 3, 4, 5].map(function (k) { return c.eq(k).rename('f_c' + k); }))
                 .addBands(img.select('baseline').rename('f_base'));
  var hit = crossBuf.filterBounds(scene.geometry());               // crossings inside this scene only
  return shares.reduceRegions({collection: hit, reducer: ee.Reducer.mean(), scale: 20, tileScale: 4})
    .map(function (f) {
      return ee.Feature(null, {
        cross_id: f.get('cross_id'), scene: scene.get('system:index'), day: img.get('day'),
        time_utc: ee.Date(img.get('system:time_start')).format("YYYY-MM-dd'T'HH:mm:ss"),
        rel_orbit: img.get('relativeOrbitNumber_start'), pass: img.get('orbitProperties_pass'),
        platform: img.get('platform_number'),
        f_c0: f.get('f_c0'), f_c1: f.get('f_c1'), f_c2: f.get('f_c2'), f_c3: f.get('f_c3'),
        f_c4: f.get('f_c4'), f_c5: f.get('f_c5'), f_base: f.get('f_base')
      });
    });
}

// ---------- 8. Exports ----------
if (MODE === 'test') {
  var clsTest = classified(TEST_BOX);
  print('Satellite passes over the test area, 2022 event:', clsTest.size());
  print('Sentinel-1 scenes over the crossings, 2022 event:', rawScenes(crossPts).size());
  Export.image.toDrive({
    image: eventSummary(clsTest).clip(TEST_BOX),
    description: 'P1_2_' + EVENT_ID + '_test_map',
    folder: DRIVE_FOLDER, region: TEST_BOX, scale: 20, crs: 'EPSG:4326', maxPixels: 1e13
  });

  for (var m = 0; m < MONTHS.length; m++) {                       // one small task per month
    var scenes = s1Base(crossPts, MONTHS[m][0], MONTHS[m][1]);
    Export.table.toDrive({
      collection: ee.FeatureCollection(scenes.map(crossingRows)).flatten(),
      description: 'P1_2_' + EVENT_ID + '_crossing_passes_m' + (m + 1),
      folder: DRIVE_FOLDER, fileFormat: 'CSV', selectors: ROW_COLS
    });
  }
}

if (MODE === 'corridor') {
  var clsCorr = classified(bufferGeom);
  print('Sentinel-1 scenes over the corridor, 2022 event:', rawScenes(bufferGeom).size());
  Export.image.toAsset({
    image: eventSummary(clsCorr).clip(bufferGeom),
    description: 'P1_2_' + EVENT_ID + '_summary_asset',
    assetId: SUMMARY_ASSET, region: bufferGeom.bounds(), scale: 20, crs: 'EPSG:4326', maxPixels: 1e13,
    pyramidingPolicy: {'.default': 'mode'}
  });
}

if (MODE === 'tables') {
  var s   = ee.Image(SUMMARY_ASSET);
  var ev  = s.select('event_class');
  var px  = ee.Image.pixelArea().divide(1e6);                       // km2
  var bands = [px.updateMask(ev.mask()).rename('area_buffer')];
  for (var k = 0; k <= 5; k++) { bands.push(px.updateMask(ev.eq(k)).rename('area_c' + k)); }
  bands.push(px.updateMask(s.select('baseline_flood').eq(1)).rename('area_baseline'));
  var perUnit = ee.Image.cat(bands).reduceRegions({collection: units, reducer: ee.Reducer.sum(), scale: 20, tileScale: 8});
  Export.table.toDrive({
    collection: perUnit, description: 'P1_2_' + EVENT_ID + '_units', folder: DRIVE_FOLDER, fileFormat: 'CSV',
    selectors: ['HYBAS_ID', 'area_buffer', 'area_c0', 'area_c1', 'area_c2', 'area_c3', 'area_c4', 'area_c5', 'area_baseline']
  });

  var code = ev.toInt().multiply(10).add(aridClass.toInt()).rename('code');
  var grouped = px.addBands(code).reduceRegion({
    reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'code'}),
    geometry: bufferGeom, scale: 60, maxPixels: 1e13, tileScale: 8
  });
  var byArid = ee.FeatureCollection(ee.List(grouped.get('groups')).map(function (d) {
    d = ee.Dictionary(d);
    var cd = ee.Number(d.get('code'));
    return ee.Feature(null, {event: EVENT_ID, event_class: cd.divide(10).floor(), arid_class: cd.mod(10),
                             area_km2: d.get('sum')});
  }));
  Export.table.toDrive({
    collection: byArid, description: 'P1_2_' + EVENT_ID + '_aridity', folder: DRIVE_FOLDER, fileFormat: 'CSV',
    selectors: ['event', 'event_class', 'arid_class', 'area_km2']
  });
}

if (MODE === 'passes') {
  var pv  = mosaicByPass(s1Base(TEST_BOX, PASS_WINDOW[0], PASS_WINDOW[1]));
  var pr  = s1Windows(TEST_BOX, REF);
  var pl  = s1Base(TEST_BOX, LONG_REF[0], LONG_REF[1]);
  var pc  = pv.map(function (img) { return classify(img, pr, pl); });
  print('Satellite passes over the test area, ' + PASS_WINDOW[0] + ' to ' + PASS_WINDOW[1] + ':', pc.size());
  var tags = pc.aggregate_array('key').map(function (k) { return ee.String(k).replace('-', '', 'g'); });
  var cBands = tags.map(function (t) { return ee.String('c_').cat(t); });    // classifier class 0-5
  var bBands = tags.map(function (t) { return ee.String('b_').cat(t); });    // naive threshold 0/1
  var stack = pc.select('cls').toBands().rename(cBands)
                .addBands(pc.select('baseline').toBands().rename(bBands)).toByte().clip(TEST_BOX);
  Export.image.toDrive({
    image: stack, description: 'P1_2_' + EVENT_ID + '_test_passes',
    folder: DRIVE_FOLDER, region: TEST_BOX, scale: 20, crs: 'EPSG:4326', maxPixels: 1e13
  });
  Export.table.toDrive({
    collection: pc.map(function (i) {
      return ee.Feature(null, {
        tag: ee.String(i.get('key')).replace('-', '', 'g'), day: i.get('day'),
        time_utc: ee.Date(i.get('system:time_start')).format("YYYY-MM-dd'T'HH:mm:ss"),
        rel_orbit: i.get('relativeOrbitNumber_start'), pass: i.get('orbitProperties_pass')
      });
    }),
    description: 'P1_2_' + EVENT_ID + '_test_passes_list', folder: DRIVE_FOLDER, fileFormat: 'CSV',
    selectors: ['tag', 'day', 'time_utc', 'rel_orbit', 'pass']
  });
}
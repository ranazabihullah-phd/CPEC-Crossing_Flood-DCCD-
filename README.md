// =====================================================================
// PAPER 1 · STEP 5 — FLOOD RECORD AT ALL 1,127 CROSSINGS, 2015-2025
// VERSION 2.2 — COMMON DENOMINATOR (FINAL)
//
// Frozen Rule (Step 4, Rule 3652):
// - Flooded: VV < -13 dB AND drop >= 6 dB AND HAND <= 15 m
// - Unobservable: Look-alike <= -24 dB (10%) OR reference count < 3
// - Permanent water: JRC occurrence >= 80%
// =====================================================================

// ---------- 1. SETTINGS ----------

var RUN = 'rest'; // 'check' = 2022 test first; 'rest' = remaining 10 years

var CROSSINGS_ASSET = 'projects/ee-zabih/assets/CPEC_Crossings_FINAL';
var DRIVE_FOLDER    = 'P1_Step05_Crossings';
var PREFIX          = 'P1_S05_Cross';
var CIRCLES_M       = [300, 150];

var ALL_YEARS = [2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];
var YEARS = RUN === 'check' ? [2022] : ALL_YEARS.filter(function(y) { return y !== 2022; });

// ---------- 2. FROZEN RULE ----------

var RULE = {
  water_vv_db: -13,
  drop_db: 6,
  hand_max_m: 15,
  season_half_days: 60,
  look_pct: 10,
  look_db: -24,
  look_db_strict: -20,
  min_ref: 3,
  perm_water_occ: 80,
  speckle_radius_m: 30,
  border_noise_db: -30
};

var STH_DB = -18;
var SNOW_DAYS = 3;
var SNOW_MODIS_MIN = 40;
var SNOW_ERA5_MIN = 50;

// ---------- 3. INPUTS AND STATIC LAYERS ----------

var crossings = ee.FeatureCollection(CROSSINGS_ASSET);
var crossPts  = crossings.geometry();

function circles(r) {
  return crossings.map(function(f) {
    return ee.Feature(f.geometry().buffer(r, 5), {cross_id: f.get('cross_id'), circle_m: r});
  });
}
var CIRCLES = circles(CIRCLES_M[0]).merge(circles(CIRCLES_M[1]));

var demCol   = ee.ImageCollection('COPERNICUS/DEM/GLO30').select('DEM');
var dem      = demCol.mosaic().setDefaultProjection(demCol.first().projection());
var slopeRad = ee.Terrain.slope(dem).multiply(Math.PI / 180);
var aspectRad= ee.Terrain.aspect(dem).multiply(Math.PI / 180);

var hand   = ee.Image('MERIT/Hydro/v1_0_1').select('hnd');
var handOk = hand.gt(RULE.hand_max_m).unmask(0, false).not(); // Unknown HAND counts as not too high

var perm   = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence').unmask(0, false).gte(RULE.perm_water_occ);

var MODIS_SNOW = ee.ImageCollection('MODIS/061/MOD10A1').merge(ee.ImageCollection('MODIS/061/MYD10A1')).select('NDSI_Snow_Cover');
var ERA5       = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR').select('snow_cover');

function toFloat(i) { return i.toFloat(); }
var EMPTY_SNOW = ee.Image.constant(0).rename('NDSI_Snow_Cover').toFloat().updateMask(ee.Image(0));
var EMPTY_ERA5 = ee.Image.constant(0).rename('snow_cover').toFloat().updateMask(ee.Image(0));

print('Run:', RUN, '| years:', YEARS.join(', '), '| crossings:', crossings.size());

// ---------- 4. SENTINEL-1 PREPARATION ----------

var PROPS = ['system:time_start', 'relativeOrbitNumber_start', 'orbitProperties_pass', 'platform_number'];

var RAW = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(crossPts)
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'));

function prep(img) {
  var sig = img.select(['VV', 'VH']);
  sig = sig.updateMask(sig.select('VV').gt(RULE.border_noise_db));
  sig = sig.focalMedian({radius: RULE.speckle_radius_m, kernelType: 'circle', units: 'meters'});
  return ee.Image(sig.copyProperties(img, PROPS)).set({day: img.date().format('YYYY-MM-dd')});
}

var ALL      = RAW.filterDate('2014-01-01', '2028-01-01').map(prep);
var EMPTY    = ee.Image.constant([0, 0]).rename(['VV', 'VH']).toFloat().updateMask(ee.Image(0));
var EMPTYCOL = ee.ImageCollection([EMPTY]);

// ---------- 5. REFERENCE IMAGE FUNCTIONS ----------

function sameOrbit(col, raw) {
  return col.filter(ee.Filter.eq('relativeOrbitNumber_start', raw.get('relativeOrbitNumber_start')))
            .filter(ee.Filter.eq('orbitProperties_pass', raw.get('orbitProperties_pass')))
            .filterBounds(raw.geometry());
}

function seasonRef(raw, yearA, yearB) {
  var t = raw.date(), h = RULE.season_half_days;
  var a = t.advance(yearA, 'year'), b = t.advance(yearB, 'year');
  return sameOrbit(ALL.filterDate(a.advance(-h, 'day'), a.advance(h, 'day'))
                      .merge(ALL.filterDate(b.advance(-h, 'day'), b.advance(h, 'day'))), raw);
}

function yearRef(raw, offset) {
  var y = raw.date().get('year').add(offset);
  return sameOrbit(ALL.filterDate(ee.Date.fromYMD(y, 1, 1), ee.Date.fromYMD(y.add(1), 1, 1)), raw);
}

// ---------- 6. REFERENCE REDUCERS ----------

var MED_COUNT = ee.Reducer.median().combine({reducer2: ee.Reducer.count(), sharedInputs: true});
var P_COUNT   = ee.Reducer.percentile([RULE.look_pct]).combine({reducer2: ee.Reducer.count(), sharedInputs: true});

function perPixel(prevCol, nextCol, reducer, band) {
  var p = prevCol.merge(EMPTYCOL).select('VV').reduce(reducer);
  var n = nextCol.merge(EMPTYCOL).select('VV').reduce(reducer);
  var prevOk = p.select('VV_count').unmask(0, false).gte(RULE.min_ref);
  return {
    value: n.select(band).blend(p.select(band).updateMask(prevOk)),
    count: n.select('VV_count').unmask(0, false).where(prevOk, p.select('VV_count').unmask(0, false)),
    useNext: prevOk.not()
  };
}

// ---------- 7. RADAR SHADOW AND LAYOVER ----------

function geometryFlags(raw) {
  var lookAz = ee.Terrain.aspect(raw.select('angle')).reduceRegion({
    reducer: ee.Reducer.median(), geometry: raw.geometry(), scale: 1000, maxPixels: 1e9, bestEffort: true
  }).get('aspect');
  var thetaI = raw.select('angle').multiply(Math.PI / 180);
  var phiR   = ee.Image.constant(ee.Number(lookAz)).multiply(Math.PI / 180).subtract(aspectRad);
  var alphaR = slopeRad.tan().multiply(phiR.cos()).atan();
  return {
    shadow: alphaR.lt(thetaI.subtract(Math.PI / 2)).unmask(0, false),
    layover: alphaR.gt(thetaI).unmask(0, false)
  };
}

// ---------- 8. SNOW FLAG ----------

function snowFlag(raw) {
  var day0 = ee.Date(raw.date().format('YYYY-MM-dd'));
  var m = MODIS_SNOW.filterDate(day0.advance(-SNOW_DAYS, 'day'), day0.advance(SNOW_DAYS + 1, 'day'))
                    .map(toFloat).merge(ee.ImageCollection([EMPTY_SNOW])).max();
  var e = ERA5.filterDate(day0, day0.advance(1, 'day')).map(toFloat).merge(ee.ImageCollection([EMPTY_ERA5])).max();
  return m.gte(SNOW_MODIS_MIN).unmask(0, false).or(e.multiply(100).gte(SNOW_ERA5_MIN).unmask(0, false));
}

// ---------- 9. EXPORT COLUMNS ----------

var COLS = [
  'cross_id', 'circle_m', 'scene', 'day', 'time_utc', 'rel_orbit', 'pass', 'platform',
  'f_valid', 'f_perm', 'f_flood', 'f_dry', 'f_unobs', 'f_look', 'f_noref', 'f_look20',
  'f_flood_look20', 'f_sth', 'f_shadow', 'f_layover', 'f_snow', 'vv_mean', 'f_ref_next', 'f_look_next'
];

// ---------- 10. CROSSING RECORD CALCULATION ----------

function crossingRows(raw) {
  var img = prep(raw);
  var ref = perPixel(seasonRef(raw, -1, -2), seasonRef(raw, 1, 2), MED_COUNT, 'VV_median');
  var lkQ = perPixel(yearRef(raw, -1), yearRef(raw, 1), P_COUNT, 'VV_p' + RULE.look_pct);
  var geo = geometryFlags(raw);

  var vv = img.select('VV');

  function circleShare(binaryImage, bandName) {
    return ee.Image(binaryImage).toFloat().unmask(0, false).rename(bandName);
  }

  var valid = vv.mask().gt(0).unmask(0, false);
  var land  = valid.and(perm.not());

  var look   = lkQ.value.lt(RULE.look_db).unmask(0, false);
  var look20 = lkQ.value.lt(RULE.look_db_strict).unmask(0, false);

  // Dual Reference Availability Check
  var seasonRefAvailable = ref.value.mask().gt(0).unmask(0, false);
  var lookRefAvailable   = lkQ.value.mask().gt(0).unmask(0, false);
  var seasonCountOk      = ref.count.unmask(0, false).gte(RULE.min_ref);
  var lookCountOk        = lkQ.count.unmask(0, false).gte(RULE.min_ref);

  var noref = seasonCountOk.not()
      .or(lookCountOk.not())
      .or(seasonRefAvailable.not())
      .or(lookRefAvailable.not())
      .unmask(1, false);

  var unobs  = land.and(look.or(noref));
  var judged = land.and(unobs.not());

  var drop   = vv.subtract(ref.value);
  var dropOk = drop.lte(-RULE.drop_db).unmask(0, false);

  var flood  = judged.and(vv.lt(RULE.water_vv_db)).and(dropOk).and(handOk).unmask(0, false);

  var bands = ee.Image.cat([
    circleShare(valid, 'f_valid'),
    circleShare(valid.and(perm), 'f_perm'),
    circleShare(flood, 'f_flood'),
    circleShare(judged.and(flood.not()), 'f_dry'),
    circleShare(unobs, 'f_unobs'),
    circleShare(land.and(look), 'f_look'),
    circleShare(land.and(noref), 'f_noref'),
    circleShare(land.and(look20), 'f_look20'),
    circleShare(flood.and(look20.not()), 'f_flood_look20'),
    circleShare(land.and(vv.lt(STH_DB)), 'f_sth'),
    circleShare(land.and(geo.shadow), 'f_shadow'),
    circleShare(land.and(geo.layover), 'f_layover'),
    circleShare(land.and(snowFlag(raw)), 'f_snow'),
    vv.rename('vv_mean'),
    circleShare(land.and(ref.useNext), 'f_ref_next'),
    circleShare(land.and(lkQ.useNext), 'f_look_next')
  ]);

  var info = {
    scene: raw.get('system:index'), day: img.get('day'),
    time_utc: raw.date().format("YYYY-MM-dd'T'HH:mm:ss"),
    rel_orbit: raw.get('relativeOrbitNumber_start'), pass: raw.get('orbitProperties_pass'),
    platform: raw.get('platform_number')
  };

  return bands.reduceRegions({
      collection: CIRCLES.filterBounds(raw.geometry()),
      reducer: ee.Reducer.mean(),
      scale: 20,
      tileScale: 4
    })
    .map(function(f) { return f.set(info); });
}

// ---------- 11. EXPORTS ----------

YEARS.forEach(function(y) {
  var scenes = RAW.filterDate(y + '-01-01', (y + 1) + '-01-01');
  print('Sentinel-1 scenes over the crossings, ' + y + ':', scenes.size());
  Export.table.toDrive({
    collection: ee.FeatureCollection(scenes.map(crossingRows)).flatten(),
    description: PREFIX + '_' + y, folder: DRIVE_FOLDER, fileFormat: 'CSV', selectors: COLS
  });
});

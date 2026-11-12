// =====================================================================
// PAPER 1 · STEP 6 (v2) — THE FROZEN FLOOD RULE AT ALL CROSSINGS, 2015-2025
// v2 FIX: the choice between the previous and the following reference years is now made PIXEL BY PIXEL.
//   (v1 chose for the whole scene; in 2016-2017 a scene could pass the 'at least 3 passes' test while the
//    crossing itself had fewer, which made 55-79% of those crossing-passes 'unobservable'.)
//   For 2019-2025 v1 and v2 give the same result (the previous years always have enough passes there),
//   so only 2015-2018 are run again.
// The rule chosen in Step 5 (rule 582) is applied to EVERY Sentinel-1 scene at EVERY crossing.
// For each crossing and scene: the share of a 300 m circle around the crossing that is
//   flooded / not flooded / unobservable / permanent water / not imaged.
// No maps are made: only small tables (this is what Paper 1 counts floods and observability from).
//
// THE RULE (do not change; it was chosen and tested in Step 5):
//   flooded       VV darker than -15 dB AND at least 6 dB darker than the same-season normal
//                 (median of the same orbit, +-60 days around the same date, in the 2 previous years)
//   unobservable  darker than -20 dB in at least 25% of the previous calendar year (look-alike surface)
//                 OR fewer than 3 reference passes
//   permanent water: JRC water occurrence >= 80%
//   Where the 2 previous years have fewer than 3 passes of the same orbit AT THAT PIXEL (early years), the
//   2 following years are used instead. Columns f_ref_next / f_look_next = share of the land circle where
//   the following years were used (0 = previous years everywhere).
//
// Extra columns for the sensitivity test and the comparison in the paper:
//   f_flood_d4, f_flood_d8   same rule with a 4 dB or 8 dB drop instead of 6 dB
//   f_simple18               the simple threshold (VV < -18 dB, no reference, no mask)
//
// HOW TO RUN:
//   RUN = 'fix'   -> 4 tasks: P1_6_Cross_2015_v2 ... P1_6_Cross_2018_v2  (this is the only run needed now)
//   ('check' = 2022 and 'rest' = the other years were run with v1; they are kept for the record)
// Outputs in Google Drive folder CPEC_Paper1_Crossings_2015_2025. Nothing is drawn on the map.
// =====================================================================

// ---------- 1. Settings ----------
var RUN = 'fix';                                         // 'fix' = 2015-2018 again with the v2 reference choice

var CROSSINGS_ASSET = 'projects/ee-zabih/assets/CPEC_Crossings_FINAL';
var DRIVE_FOLDER    = 'CPEC_Paper1_Crossings_2015_2025';
var PREFIX          = 'P1_6_Cross';
var CROSS_BUFFER_M  = 300;                               // same circle as Steps 2-3
var ALL_YEARS = [2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];
var YEARS = RUN === 'fix' ? [2015, 2016, 2017, 2018]
          : RUN === 'check' ? [2022] : ALL_YEARS.filter(function (y) { return y !== 2022; });
var SUFFIX = RUN === 'fix' ? '_v2' : '';                 // new file names: nothing earlier is overwritten

// Frozen rule (Step 5, rule 582)
var RULE = {
  water_vv_db: -15,       // flooded only if VV is darker than this ...
  drop_db: 6,             // ... and at least this many dB below the same-season normal
  season_half_days: 60,   // same season = +-60 days around the same date in the reference years
  look_pct: 25,           // look-alike: darker than look_db in at least 25% of the year
  look_db: -20,
  min_ref: 3,             // fewer reference passes -> unobservable
  perm_water_occ: 80,     // JRC occurrence >= 80% -> permanent water
  speckle_radius_m: 30,   // same radar preparation as Steps 2-5
  border_noise_db: -30
};
var ALT_DROPS = [4, 8];                                  // sensitivity columns
var SIMPLE_DB = -18;                                     // simple threshold, for comparison only

// ---------- 2. Inputs ----------
var crossings = ee.FeatureCollection(CROSSINGS_ASSET);
var crossPts  = crossings.geometry();
var crossBuf  = crossings.map(function (f) {
  return ee.Feature(f.geometry().buffer(CROSS_BUFFER_M, 5), {cross_id: f.get('cross_id')});
});
var perm = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence').unmask(0).gte(RULE.perm_water_occ);
print('Run:', RUN, '| years:', YEARS.join(', '), '| crossings:', crossings.size());

// ---------- 3. Sentinel-1 (same preparation as Steps 2-5) ----------
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
var ALL = RAW.filterDate('2014-01-01', '2027-01-01').map(prep);     // reference images (prepared on demand)

var EMPTY = ee.Image.constant([0, 0]).rename(['VV', 'VH']).toFloat().updateMask(ee.Image(0));

function sameOrbit(col, raw) {
  return col.filter(ee.Filter.eq('relativeOrbitNumber_start', raw.get('relativeOrbitNumber_start')))
            .filter(ee.Filter.eq('orbitProperties_pass', raw.get('orbitProperties_pass')))
            .filterBounds(raw.geometry());
}

// same season (+-60 days around the same date) in two other years, same orbit
function seasonRef(raw, yearA, yearB) {
  var t = raw.date();
  var h = RULE.season_half_days;
  var a = t.advance(yearA, 'year');
  var b = t.advance(yearB, 'year');
  var col = ALL.filterDate(a.advance(-h, 'day'), a.advance(h, 'day'))
               .merge(ALL.filterDate(b.advance(-h, 'day'), b.advance(h, 'day')));
  return sameOrbit(col, raw);
}

// one calendar year of the same orbit (look-alike layer)
function yearRef(raw, offset) {
  var y = raw.date().get('year').add(offset);
  return sameOrbit(ALL.filterDate(ee.Date.fromYMD(y, 1, 1), ee.Date.fromYMD(y.add(1), 1, 1)), raw);
}

// ---------- 4. One scene -> one row per crossing inside it ----------
var COLS = ['cross_id', 'scene', 'day', 'time_utc', 'rel_orbit', 'pass', 'platform', 'ref_years', 'look_years',
            'f_valid', 'f_perm', 'f_flood', 'f_dry', 'f_unobs', 'f_look', 'f_noref',
            'f_flood_d4', 'f_flood_d8', 'f_simple18', 'vv_mean', 'f_ref_next', 'f_look_next'];

var EMPTYCOL  = ee.ImageCollection([EMPTY]);
var MED_COUNT = ee.Reducer.median().combine({reducer2: ee.Reducer.count(), sharedInputs: true});
var P_COUNT   = ee.Reducer.percentile([RULE.look_pct]).combine({reducer2: ee.Reducer.count(), sharedInputs: true});

// per pixel: the previous years where they have >= 3 passes, otherwise the following years
function perPixel(prevCol, nextCol, reducer, band) {
  var p = prevCol.merge(EMPTYCOL).select('VV').reduce(reducer);
  var n = nextCol.merge(EMPTYCOL).select('VV').reduce(reducer);
  var prevOk = p.select('VV_count').unmask(0).gte(RULE.min_ref);
  return {
    value: n.select(band).blend(p.select(band).updateMask(prevOk)),
    count: n.select('VV_count').unmask(0).where(prevOk, p.select('VV_count').unmask(0)),
    useNext: prevOk.not()
  };
}

function crossingRows(raw) {
  var img  = prep(raw);
  var ref  = perPixel(seasonRef(raw, -1, -2), seasonRef(raw, 1, 2), MED_COUNT, 'VV_median');
  var lkQ  = perPixel(yearRef(raw, -1), yearRef(raw, 1), P_COUNT, 'VV_p' + RULE.look_pct);

  var vv     = img.select('VV');
  var valid  = vv.mask().gt(0).unmask(0);                                   // imaged by this scene
  var land   = valid.and(perm.not());
  var look   = lkQ.value.lt(RULE.look_db).unmask(0);
  var noref  = ref.count.lt(RULE.min_ref);
  var unobs  = land.and(look.or(noref));
  var judged = land.and(unobs.not());
  var drop   = vv.subtract(ref.value);
  var dark   = vv.lt(RULE.water_vv_db);
  function floodWith(d) { return judged.and(dark).and(drop.lte(-d)).unmask(0); }
  var flood  = floodWith(RULE.drop_db);

  var bands = ee.Image.cat([
    valid.rename('f_valid'),
    valid.and(perm).rename('f_perm'),
    flood.rename('f_flood'),
    judged.and(flood.not()).rename('f_dry'),
    unobs.rename('f_unobs'),
    land.and(look).rename('f_look'),
    land.and(noref).rename('f_noref'),
    floodWith(ALT_DROPS[0]).rename('f_flood_d' + ALT_DROPS[0]),
    floodWith(ALT_DROPS[1]).rename('f_flood_d' + ALT_DROPS[1]),
    land.and(vv.lt(SIMPLE_DB)).unmask(0).rename('f_simple18'),
    vv.rename('vv_mean'),                                                 // mean VV of the imaged part (dB)
    land.and(ref.useNext).rename('f_ref_next'),
    land.and(lkQ.useNext).rename('f_look_next')
  ]);

  var info = {
    scene: raw.get('system:index'), day: img.get('day'),
    time_utc: raw.date().format("YYYY-MM-dd'T'HH:mm:ss"),
    rel_orbit: raw.get('relativeOrbitNumber_start'), pass: raw.get('orbitProperties_pass'),
    platform: raw.get('platform_number'),
    ref_years: 'per-pixel',
    look_years: 'per-pixel'
  };
  return bands.reduceRegions({collection: crossBuf.filterBounds(raw.geometry()), reducer: ee.Reducer.mean(),
                              scale: 20, tileScale: 4})
    .map(function (f) { return f.set(info); });
}

// ---------- 5. Exports: one task per year ----------
YEARS.forEach(function (y) {
  var scenes = RAW.filterDate(y + '-01-01', (y + 1) + '-01-01');
  print('Sentinel-1 scenes over the crossings, ' + y + ':', scenes.size());
  Export.table.toDrive({
    collection: ee.FeatureCollection(scenes.map(crossingRows)).flatten(),
    description: PREFIX + '_' + y + SUFFIX, folder: DRIVE_FOLDER, fileFormat: 'CSV', selectors: COLS
  });
});
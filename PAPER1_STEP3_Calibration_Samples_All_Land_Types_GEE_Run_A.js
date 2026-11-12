// =====================================================================
// PAPER 1 · STEP 3 — CALIBRATION SAMPLES FOR ALL LAND TYPES  (GEE run A)
//
// Purpose: one data set that lets Step 4 (Python) choose the final flood rule (DCCD) with EVERY land type
// of the study area tested, and compare it with STH, STH-M and ZCD. Nothing is decided here.
//
// Random points in 7 test cases. For every point and every Sentinel-1 pass in the case window, one row with
// the radar values, the reference values and the extra layers:
//   case             area                                   window                question (check)
//   flood2022        Sindh test area                        24 Aug - 3 Sep 2022   finds real flood? (A, vs UNOSAT)
//   normal2021       Sindh test area (same points)          24 Aug - 3 Sep 2021   seasonal water, normal year (D)
//   chagai           Chagai desert, N-40                    Oct - Nov 2021        dry sand (C)
//   makran           Makran coast, N-10                     Oct - Nov 2021        coast, mudflats (C)
//   mountain_dry     Kaghan - Babusar valley, N-15          Oct - Nov 2021        mountains, dry (C)
//   mountain_spring  Kaghan - Babusar valley (same points)  Mar - Apr 2021        mountains, wet snow (C)
//   riverbed_dry     Kurram river near Bannu, N-55 + M-14   Oct - Nov 2021        moving river channels (C)
//
// Columns (radar values in dB x 100, e.g. -1834 = -18.34 dB; -9999 = no value):
//   vv, vh                         this pass (30 m speckle filter, as in all earlier steps)
//   s_vv_med ... s_n               same-season reference: same orbit, the 2 previous years
//   d_vv_med ... d_n               dry-season reference: same orbit (Mar-May of the same year; Oct-Nov 2020 for spring)
//   l_p10 ... l_p75, l_n, l_dark18 12-month reference: VV percentiles, number of passes, % of passes darker than -18 dB
//   angle                          radar incidence angle (degrees x 10)
//   shadow, layover                1 = the radar cannot see this pixel on this pass (from the terrain and the pass geometry)
//   snow_modis                     highest MODIS snow cover (NDSI x 100, Terra + Aqua) within +-3 days; -1 = always cloudy
//   snow_era5                      ERA5-Land snow cover on the pass day (%); t2m = ERA5-Land air temperature (deg C x 10)
//   slope (deg x 10), hand (m x 10), jrc (% of time water 1984-2021), jrc_max (1 = water at least once), wc (land cover)
//   cls                            class of the earlier ZCD classifier (for the comparison table)
//
// v2: snow and weather images are cast to one pixel type (Float) before they are combined
//     (v1 failed: MODIS snow is stored as 'Short' and cannot be mixed with Float images).
// Outputs: Google Drive folder P1_Step03_Calibration_Samples (8 tasks)
//   P1_S03_Samples_<case>.csv (7 files) and P1_S03_Rain_daily.csv (CHIRPS rain per case, to prove the dry windows)
// How to run: NEW Code Editor script -> paste -> Run -> press RUN on the 8 tasks in the Tasks tab. Nothing is drawn.
// =====================================================================

// ---------- 1. Settings ----------
var DRIVE_FOLDER = 'P1_Step03_Calibration_Samples';
var PREFIX       = 'P1_S03';
var NODATA       = -9999;
var SNOW_DAYS    = 3;                                               // MODIS snow is searched +-3 days around the pass

var SINDH    = ee.Geometry.Rectangle([67.3, 26.3, 68.5, 28.6]);     // same test area as before (same points)
var MOUNTAIN = ee.Geometry.Rectangle([73.4, 34.6, 74.1, 35.2]);     // Kaghan - Babusar, N-15 (1,336 - 3,419 m)
var RIVERBED = ee.Geometry.Rectangle([70.6, 32.45, 71.45, 32.95]);  // Kurram river near Bannu, N-55 and M-14

var DRY21  = {window: ['2021-10-01', '2021-12-01'],
              sameRef: [['2020-10-01', '2020-12-01'], ['2019-10-01', '2019-12-01']],
              dryRef: ['2021-03-01', '2021-06-01'], longRef: ['2020-10-01', '2021-10-01']};
var CASES = [
  {name: 'flood2022', geom: SINDH, points: 40000, seed: 51, window: ['2022-08-24', '2022-09-04'],
   sameRef: [['2020-06-15', '2020-10-15'], ['2021-06-15', '2021-10-15']],
   dryRef: ['2022-03-01', '2022-06-01'], longRef: ['2021-01-01', '2022-01-01']},
  {name: 'normal2021', geom: SINDH, points: 40000, seed: 51, window: ['2021-08-24', '2021-09-04'],
   sameRef: [['2019-06-15', '2019-10-15'], ['2020-06-15', '2020-10-15']],
   dryRef: ['2021-03-01', '2021-06-01'], longRef: ['2020-01-01', '2021-01-01']},
  {name: 'chagai', geom: ee.Geometry.Rectangle([63.8, 28.5, 64.8, 29.1]), points: 8000, seed: 52, ref: DRY21},
  {name: 'makran', geom: ee.Geometry.Rectangle([63.4, 25.15, 64.8, 25.5]), points: 12000, seed: 53, ref: DRY21},
  {name: 'mountain_dry', geom: MOUNTAIN, points: 12000, seed: 54, ref: DRY21},
  {name: 'mountain_spring', geom: MOUNTAIN, points: 12000, seed: 54, window: ['2021-03-01', '2021-05-01'],
   sameRef: [['2020-03-01', '2020-05-01'], ['2019-03-01', '2019-05-01']],
   dryRef: ['2020-10-01', '2020-12-01'], longRef: ['2020-03-01', '2021-03-01']},
  {name: 'riverbed_dry', geom: RIVERBED, points: 20000, seed: 55, ref: DRY21}
];
CASES.forEach(function (c) {                                        // cases that share the dry Oct-Nov 2021 windows
  if (c.ref) { c.window = c.ref.window; c.sameRef = c.ref.sameRef; c.dryRef = c.ref.dryRef; c.longRef = c.ref.longRef; }
});

// Same rules as the earlier steps (only used for the 'cls' column = ZCD)
var F = {
  speckle_radius_m: 30, border_noise_db: -30, z_k: 2.5, z_k_vh: 2.0, abs_vv_db: -15, z_k_strict: 3.5,
  min_std_db: 0.5, min_ref_count: 3, lookalike_vv_db: -18, lookalike_freq: 0.5, hand_max_m: 15,
  slope_plaus_deg: 5, slope_terrain_deg: 15, perm_water_occ: 80
};

// ---------- 2. Static layers ----------
var demCol = ee.ImageCollection('COPERNICUS/DEM/GLO30').select('DEM');
var dem    = demCol.mosaic().setDefaultProjection(demCol.first().projection());
var slope  = ee.Terrain.slope(dem);
var aspect = ee.Terrain.aspect(dem);
var hand   = ee.Image('MERIT/Hydro/v1_0_1').select('hnd');
var wc     = ee.ImageCollection('ESA/WorldCover/v200').first().select('Map');
var gsw    = ee.Image('JRC/GSW1_4/GlobalSurfaceWater');
var occ    = gsw.select('occurrence').unmask(0);
var perm        = occ.gte(F.perm_water_occ);
var implausible = hand.gt(F.hand_max_m).or(slope.gt(F.slope_plaus_deg)).unmask(0);
var terrainBad  = slope.gt(F.slope_terrain_deg).unmask(0);
var hardSurface = wc.eq(50).or(wc.eq(10)).unmask(0);
var slopeRad  = slope.multiply(Math.PI / 180);
var aspectRad = aspect.multiply(Math.PI / 180);

var STATIC = ee.Image.cat([
  slope.multiply(10).round().toInt().unmask(NODATA).rename('slope'),
  hand.multiply(10).round().toInt().unmask(NODATA).rename('hand'),
  occ.round().toInt().rename('jrc'),
  gsw.select('max_extent').unmask(0).toInt().rename('jrc_max'),
  wc.toInt().unmask(0).rename('wc')
]);

// snow and air temperature
var MODIS_SNOW = ee.ImageCollection('MODIS/061/MOD10A1').merge(ee.ImageCollection('MODIS/061/MYD10A1'))
                   .select('NDSI_Snow_Cover');                      // stored as 'Short': cast below before combining
var ERA5 = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR').select(['snow_cover', 'temperature_2m']);
var EMPTY_SNOW = ee.Image.constant(0).rename('NDSI_Snow_Cover').toFloat().updateMask(ee.Image(0));
function toFloat(i) { return i.toFloat(); }                       // one pixel type for every image that is combined
var EMPTY_ERA5 = ee.Image.constant([0, 0]).rename(['snow_cover', 'temperature_2m']).toFloat().updateMask(ee.Image(0));

// ---------- 3. Sentinel-1 (same preparation as all earlier steps; the incidence angle is kept) ----------
var PROPS = ['system:time_start', 'relativeOrbitNumber_start', 'orbitProperties_pass', 'platform_number'];

function s1Base(region, start, end) {
  return ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(region).filterDate(start, end)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
    .map(function (img) {
      var sig = img.select(['VV', 'VH']);
      sig = sig.updateMask(sig.select('VV').gt(F.border_noise_db));
      sig = sig.focalMedian({radius: F.speckle_radius_m, kernelType: 'circle', units: 'meters'});
      var day = img.date().format('YYYY-MM-dd');
      var ro  = ee.Number(img.get('relativeOrbitNumber_start')).format('%d');
      return ee.Image(sig.addBands(img.select('angle')).copyProperties(img, PROPS))
               .set({day: day, key: ee.String(day).cat('_').cat(ro)});
    });
}

function s1Windows(region, windows) {
  var out = s1Base(region, windows[0][0], windows[0][1]);
  for (var w = 1; w < windows.length; w++) { out = out.merge(s1Base(region, windows[w][0], windows[w][1])); }
  return out;
}

// one image per pass; look_az = direction from the ground towards the satellite (the downhill direction of the
// incidence-angle field), computed once per pass from its first scene (native projection)
function mosaicByPass(col, region) {
  var keys = col.aggregate_array('key').distinct();
  return ee.ImageCollection.fromImages(keys.map(function (k) {
    var sub = col.filter(ee.Filter.eq('key', k));
    var first = ee.Image(sub.first());
    var lookAz = ee.Terrain.aspect(first.select('angle')).reduceRegion({
      reducer: ee.Reducer.median(), geometry: region, scale: 1000, maxPixels: 1e9, bestEffort: true
    }).get('aspect');
    return ee.Image(sub.mosaic().copyProperties(first, PROPS.concat(['day', 'key']))).set('look_az', lookAz);
  }));
}

var EMPTY = ee.Image.constant([0, 0]).rename(['VV', 'VH']).toFloat().updateMask(ee.Image(0));
function sameGeometry(col, img) {
  return col.filter(ee.Filter.eq('relativeOrbitNumber_start', img.get('relativeOrbitNumber_start')))
            .filter(ee.Filter.eq('orbitProperties_pass', img.get('orbitProperties_pass')))
            .select(['VV', 'VH'])
            .merge(ee.ImageCollection([EMPTY]));
}

// ---------- 4. All values for one pass ----------
var STATS = ee.Reducer.median().combine({reducer2: ee.Reducer.stdDev(), sharedInputs: true})
                               .combine({reducer2: ee.Reducer.count(), sharedInputs: true});
var QUANT = ee.Reducer.percentile([10, 25, 50, 75]).combine({reducer2: ee.Reducer.count(), sharedInputs: true});

function cdb(x) { return x.multiply(100).round().toInt().unmask(NODATA); }    // dB -> dB x 100, gaps -> -9999

function refBands(stats, p) {
  return ee.Image.cat([
    cdb(stats.select('VV_median')).rename(p + 'vv_med'),
    cdb(stats.select('VV_stdDev')).rename(p + 'vv_sd'),
    cdb(stats.select('VH_median')).rename(p + 'vh_med'),
    cdb(stats.select('VH_stdDev')).rename(p + 'vh_sd'),
    stats.select('VV_count').toInt().unmask(0).rename(p + 'n')
  ]);
}

// radar shadow and layover from the terrain and this pass's geometry (Vollrath et al. 2020, Remote Sensing 12, 1867)
function geometryFlags(img) {
  var thetaI = img.select('angle').multiply(Math.PI / 180);
  var phiR   = ee.Image.constant(ee.Number(img.get('look_az'))).multiply(Math.PI / 180).subtract(aspectRad);
  var alphaR = slopeRad.tan().multiply(phiR.cos()).atan();              // terrain slope in the radar's range direction
  return ee.Image.cat([
    alphaR.lt(thetaI.subtract(Math.PI / 2)).unmask(0).toInt().rename('shadow'),   // slope facing away, steeper than 90 deg - incidence
    alphaR.gt(thetaI).unmask(0).toInt().rename('layover')                          // slope facing the radar, steeper than the incidence
  ]);
}

function snowBands(img) {
  var day0 = ee.Date(ee.String(img.get('day')));                    // calendar day of the pass (00:00 UTC)
  var m = MODIS_SNOW.filterDate(day0.advance(-SNOW_DAYS, 'day'), day0.advance(SNOW_DAYS + 1, 'day'))
                    .map(toFloat).merge(ee.ImageCollection([EMPTY_SNOW])).max();
  var e = ERA5.filterDate(day0, day0.advance(1, 'day')).map(toFloat).merge(ee.ImageCollection([EMPTY_ERA5])).max();
  return ee.Image.cat([
    m.round().toInt().unmask(-1).rename('snow_modis'),
    e.select('snow_cover').multiply(100).round().toInt().unmask(-1).rename('snow_era5'),
    e.select('temperature_2m').subtract(273.15).multiply(10).round().toInt().unmask(NODATA).rename('t2m')
  ]);
}

function passStack(img, sCol, dCol, lCol) {
  var sRef   = sameGeometry(sCol, img);
  var lRef   = sameGeometry(lCol, img);
  var sStats = sRef.reduce(STATS);
  var dStats = sameGeometry(dCol, img).reduce(STATS);
  var lQuant = lRef.select('VV').reduce(QUANT);
  var lkFreq = lRef.map(function (i) { return i.select('VV').lt(F.lookalike_vv_db); }).mean();

  // earlier ZCD classifier (identical rules to the earlier steps) -> 'cls' (0-5)
  var zVV = img.select('VV').subtract(sStats.select('VV_median')).divide(sStats.select('VV_stdDev').max(F.min_std_db));
  var zVH = img.select('VH').subtract(sStats.select('VH_median')).divide(sStats.select('VH_stdDev').max(F.min_std_db));
  var lookalike = lkFreq.unmask(0).gte(F.lookalike_freq);
  var refOk  = sStats.select('VV_count').unmask(0).gte(F.min_ref_count);
  var dark   = img.select('VV').lte(F.abs_vv_db);
  var change = zVV.lte(-F.z_k).and(dark).and(zVH.lte(-F.z_k_vh));
  var strict = zVV.lte(-F.z_k_strict).and(dark);
  var valid  = img.select('VV').mask().gt(0);
  var cls = ee.Image(0).where(valid, 1)
    .where(valid.and(refOk.not()), 5)
    .where(valid.and(refOk).and(terrainBad).and(change.not()), 5)
    .where(valid.and(refOk).and(lookalike).and(strict.not()), 5);
  var lowConf = lookalike.or(implausible).or(hardSurface).or(terrainBad);
  var flood   = valid.and(refOk).and(change).and(lookalike.not().or(strict));
  cls = cls.where(flood.and(lowConf.not()), 2).where(flood.and(lowConf), 3).where(valid.and(perm), 4);

  return ee.Image.cat([
    img.select('VV').multiply(100).round().toInt().rename('vv'),           // keeps the pass mask:
    img.select('VH').multiply(100).round().toInt().rename('vh'),           // points outside the pass are dropped
    refBands(sStats, 's_'),
    refBands(dStats, 'd_'),
    cdb(lQuant.select('VV_p10')).rename('l_p10'),
    cdb(lQuant.select('VV_p25')).rename('l_p25'),
    cdb(lQuant.select('VV_p50')).rename('l_p50'),
    cdb(lQuant.select('VV_p75')).rename('l_p75'),
    lQuant.select('VV_count').toInt().unmask(0).rename('l_n'),
    lkFreq.multiply(100).round().toInt().unmask(NODATA).rename('l_dark18'),
    img.select('angle').multiply(10).round().toInt().unmask(NODATA).rename('angle'),
    geometryFlags(img),
    snowBands(img),
    STATIC,
    cls.toInt().rename('cls')
  ]);
}

function passRows(img, caseName, pts, sCol, dCol, lCol) {
  var info = {case_name: caseName, day: img.get('day'),
              rel_orbit: img.get('relativeOrbitNumber_start'), pass: img.get('orbitProperties_pass')};
  return passStack(img, sCol, dCol, lCol)
    .sampleRegions({collection: pts, scale: 20, tileScale: 4, geometries: false})
    .map(function (f) { return f.set(info); });
}

function randomPoints(c) {
  return ee.FeatureCollection.randomPoints(c.geom, c.points, c.seed, 1).map(function (f) {
    var xy = f.geometry().coordinates();
    return f.set({pid: f.get('system:index'), lon: xy.get(0), lat: xy.get(1)});
  });
}

// ---------- 5. Exports ----------
var COLS = ['case_name', 'pid', 'lon', 'lat', 'day', 'rel_orbit', 'pass', 'vv', 'vh',
            's_vv_med', 's_vv_sd', 's_vh_med', 's_vh_sd', 's_n',
            'd_vv_med', 'd_vv_sd', 'd_vh_med', 'd_vh_sd', 'd_n',
            'l_p10', 'l_p25', 'l_p50', 'l_p75', 'l_n', 'l_dark18',
            'angle', 'shadow', 'layover', 'snow_modis', 'snow_era5', 't2m',
            'slope', 'hand', 'jrc', 'jrc_max', 'wc', 'cls'];

var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY').select('precipitation');
var rain = ee.FeatureCollection([]);

CASES.forEach(function (c) {
  var pts    = randomPoints(c);
  var passes = mosaicByPass(s1Base(c.geom, c.window[0], c.window[1]), c.geom);
  var sCol   = s1Windows(c.geom, c.sameRef);
  var dCol   = s1Base(c.geom, c.dryRef[0], c.dryRef[1]);
  var lCol   = s1Base(c.geom, c.longRef[0], c.longRef[1]);
  print('Satellite passes, ' + c.name + ':', passes.size());

  var rows = ee.FeatureCollection(passes.map(function (img) {
    return passRows(img, c.name, pts, sCol, dCol, lCol);
  })).flatten();
  Export.table.toDrive({
    collection: rows, description: PREFIX + '_Samples_' + c.name, folder: DRIVE_FOLDER, fileFormat: 'CSV', selectors: COLS
  });

  var days = chirps.filterDate(ee.Date(c.window[0]).advance(-15, 'day'), c.window[1]);
  rain = rain.merge(ee.FeatureCollection(days.map(function (d) {
    var mm = d.reduceRegion({reducer: ee.Reducer.mean(), geometry: c.geom, scale: 5566, maxPixels: 1e9})
              .get('precipitation');
    return ee.Feature(null, {case_name: c.name, date: d.date().format('YYYY-MM-dd'), rain_mm: mm});
  })));
});

Export.table.toDrive({
  collection: rain, description: PREFIX + '_Rain_daily', folder: DRIVE_FOLDER, fileFormat: 'CSV',
  selectors: ['case_name', 'date', 'rain_mm']
});
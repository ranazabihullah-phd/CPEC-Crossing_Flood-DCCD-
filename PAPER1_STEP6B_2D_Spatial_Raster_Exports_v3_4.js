// ============================================================================
// PAPER 1 · STEP 6B — 2D SPATIAL RASTER EXPORTS (v3.4 — one complete pass per area)
// Fix vs v3.3: rawCol.mosaic() joined ALL images of the date window, so it could mix
//   different passes (orbits) while the reference was taken from the first image only.
//   Now the mosaic holds ONE pass: all frames of the SAME day AND the SAME relative orbit,
//   so every pixel is compared with the reference of its own orbit.
// For each area the pass is chosen among the Step 3 dry-pass days (Sindh: 27 Aug 2022 only)
//   as the one that covers the LARGEST part of the box (coverage does not depend on the
//   flood result, so this choice is neutral). The console prints the chosen day:
//   copy it into AREAS (Cell 2) of the Step 6B notebook.
// RUN lists the areas exported now (default: Chagai and Makran). The other four keep their
//   current rasters (their date windows held one pass only).
// ============================================================================

var DRIVE_FOLDER = 'P1_Step06_Rasters';
var PREFIX       = 'P1_S06_Raster';
var RUN          = ['Chagai', 'Makran', 'Sindh', 'Kaghan_Autumn', 'Kaghan_Spring', 'Kurram'];   // add 'Sindh', 'Kaghan_Autumn', 'Kaghan_Spring', 'Kurram' to redo them

// ---------- 1. FROZEN RULE 3652 PARAMETERS ----------
var RULE = {
  water_vv_db: -13, drop_db: 6, hand_max_m: 15, season_half_days: 60,
  look_pct: 10, look_db: -24, min_ref: 3,
  perm_water_occ: 80, speckle_radius_m: 30, border_noise_db: -30
};
var STH_DB = -18, SNOW_DAYS = 3, SNOW_MODIS_MIN = 40, SNOW_ERA5_MIN = 50;

// ---------- 2. TEST BOXES (Step 3) AND CANDIDATE DAYS (Step 3 dry passes) ----------
var BOX = {
  Sindh:         ee.Geometry.Rectangle([67.3, 26.3, 68.5, 28.6]),
  Chagai:        ee.Geometry.Rectangle([63.8, 28.5, 64.8, 29.1]),
  Makran:        ee.Geometry.Rectangle([63.4, 25.15, 64.8, 25.5]),
  Kaghan_Autumn: ee.Geometry.Rectangle([73.4, 34.6, 74.1, 35.2]),
  Kaghan_Spring: ee.Geometry.Rectangle([73.4, 34.6, 74.1, 35.2]),
  Kurram:        ee.Geometry.Rectangle([70.6, 32.45, 71.45, 32.95])
};
var CANDIDATES = {
  Sindh:         ['2022-08-27'],                                   // same day as the UNOSAT map
  Chagai:        ['2021-10-05', '2021-10-08', '2021-10-12', '2021-10-17', '2021-10-20', '2021-10-24', '2021-10-29',
                  '2021-11-01', '2021-11-05', '2021-11-10', '2021-11-13', '2021-11-17', '2021-11-22', '2021-11-25', '2021-11-29'],
  Makran:        ['2021-10-03', '2021-10-05', '2021-10-08', '2021-10-12', '2021-10-15', '2021-10-17', '2021-10-20',
                  '2021-10-24', '2021-10-27', '2021-10-29', '2021-11-01', '2021-11-08', '2021-11-10', '2021-11-13',
                  '2021-11-17', '2021-11-20', '2021-11-22', '2021-11-25', '2021-11-29'],
  Kaghan_Autumn: ['2021-11-17', '2021-11-28', '2021-11-29'],
  Kaghan_Spring: ['2021-04-02', '2021-04-03'],
  Kurram:        ['2021-10-04', '2021-10-05', '2021-10-09', '2021-10-16', '2021-10-17', '2021-10-21', '2021-10-28',
                  '2021-10-29', '2021-11-02', '2021-11-09', '2021-11-10', '2021-11-14', '2021-11-21', '2021-11-22', '2021-11-26']
};

// ---------- 3. STATIC ANCILLARY LAYERS ----------
var demCol    = ee.ImageCollection('COPERNICUS/DEM/GLO30').select('DEM');
var dem       = demCol.mosaic().setDefaultProjection(demCol.first().projection());
var slopeRad  = ee.Terrain.slope(dem).multiply(Math.PI / 180);
var aspectRad = ee.Terrain.aspect(dem).multiply(Math.PI / 180);
var hand      = ee.Image('MERIT/Hydro/v1_0_1').select('hnd');
var handOk    = hand.gt(RULE.hand_max_m).unmask(0, false).not();
var perm      = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence').unmask(0, false).gte(RULE.perm_water_occ);
var MODIS_SNOW = ee.ImageCollection('MODIS/061/MOD10A1').merge(ee.ImageCollection('MODIS/061/MYD10A1')).select('NDSI_Snow_Cover');
var ERA5       = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR').select('snow_cover');
function toFloat(i) { return i.toFloat(); }
var EMPTY_SNOW = ee.Image.constant(0).rename('NDSI_Snow_Cover').toFloat().updateMask(ee.Image(0));
var EMPTY_ERA5 = ee.Image.constant(0).rename('snow_cover').toFloat().updateMask(ee.Image(0));

// ---------- 4. SENTINEL-1 PREPARATION ----------
var PROPS = ['system:time_start', 'relativeOrbitNumber_start', 'orbitProperties_pass', 'platform_number'];
var RAW = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'));
function prep(img) {
  var sig = img.select(['VV', 'VH']);
  sig = sig.updateMask(sig.select('VV').gt(RULE.border_noise_db));
  sig = sig.focalMedian({radius: RULE.speckle_radius_m, kernelType: 'circle', units: 'meters'});
  return ee.Image(sig.addBands(img.select('angle')).copyProperties(img, PROPS)).set({day: img.date().format('YYYY-MM-dd')});
}
var ALL      = RAW.filterDate('2014-01-01', '2028-01-01').map(prep);
var EMPTY    = ee.Image.constant([0, 0]).rename(['VV', 'VH']).toFloat().updateMask(ee.Image(0));
var EMPTYCOL = ee.ImageCollection([EMPTY]);

// ---------- 5. REFERENCE FUNCTIONS (unchanged from v3.3) ----------
function sameOrbit(col, refImg, geom) {
  return col.filter(ee.Filter.eq('relativeOrbitNumber_start', refImg.get('relativeOrbitNumber_start')))
            .filter(ee.Filter.eq('orbitProperties_pass', refImg.get('orbitProperties_pass')))
            .filterBounds(geom);
}
function seasonRef(refImg, yearA, yearB, geom) {
  var t = refImg.date(), h = RULE.season_half_days;
  var a = t.advance(yearA, 'year'), b = t.advance(yearB, 'year');
  return sameOrbit(ALL.filterDate(a.advance(-h, 'day'), a.advance(h, 'day'))
                      .merge(ALL.filterDate(b.advance(-h, 'day'), b.advance(h, 'day'))), refImg, geom);
}
function yearRef(refImg, offset, geom) {
  var y = refImg.date().get('year').add(offset);
  return sameOrbit(ALL.filterDate(ee.Date.fromYMD(y, 1, 1), ee.Date.fromYMD(y.add(1), 1, 1)), refImg, geom);
}
var MED_COUNT = ee.Reducer.median().combine({reducer2: ee.Reducer.count(), sharedInputs: true});
var P_COUNT   = ee.Reducer.percentile([RULE.look_pct]).combine({reducer2: ee.Reducer.count(), sharedInputs: true});
function perPixel(prevCol, nextCol, reducer, band) {
  var p = prevCol.merge(EMPTYCOL).select('VV').reduce(reducer);
  var n = nextCol.merge(EMPTYCOL).select('VV').reduce(reducer);
  var prevOk = p.select('VV_count').unmask(0, false).gte(RULE.min_ref);
  return {
    value: n.select(band).blend(p.select(band).updateMask(prevOk)),
    count: n.select('VV_count').unmask(0, false).where(prevOk, p.select('VV_count').unmask(0, false))
  };
}

// ---------- 6. GEOMETRY AND SNOW FLAGS ----------
function geometryFlags(mosaicImg, geom) {
  var lookDict = ee.Terrain.aspect(mosaicImg.select('angle')).reduceRegion({
    reducer: ee.Reducer.median(), geometry: geom, scale: 1000, maxPixels: 1e9, bestEffort: true
  });
  var lookVal = ee.Number(ee.Algorithms.If(lookDict.get('aspect'), lookDict.get('aspect'), 0));
  var thetaI  = mosaicImg.select('angle').multiply(Math.PI / 180);
  var phiR    = ee.Image.constant(lookVal).multiply(Math.PI / 180).subtract(aspectRad);
  var alphaR  = slopeRad.tan().multiply(phiR.cos()).atan();
  return {
    shadow:  alphaR.lt(thetaI.subtract(Math.PI / 2)).unmask(0, false),
    layover: alphaR.gt(thetaI).unmask(0, false)
  };
}
function snowFlag(refImg) {
  var day0 = ee.Date(refImg.date().format('YYYY-MM-dd'));
  var m = MODIS_SNOW.filterDate(day0.advance(-SNOW_DAYS, 'day'), day0.advance(SNOW_DAYS + 1, 'day'))
                    .map(toFloat).merge(ee.ImageCollection([EMPTY_SNOW])).max();
  var e = ERA5.filterDate(day0, day0.advance(1, 'day')).map(toFloat).merge(ee.ImageCollection([EMPTY_ERA5])).max();
  return m.gte(SNOW_MODIS_MIN).unmask(0, false).or(e.multiply(100).gte(SNOW_ERA5_MIN).unmask(0, false));
}

// ---------- 7. CHOOSE THE PASS: candidate day + orbit with the largest share of the box ----------
function bestPass(days, geom) {
  var boxArea = geom.area(100);
  var feats = ee.List(days).map(function (d) {
    var t0  = ee.Date(d);
    var col = RAW.filterBounds(geom).filterDate(t0, t0.advance(1, 'day'));   // raw scenes: exact footprints
    var orbits = ee.List(col.aggregate_array('relativeOrbitNumber_start')).distinct();
    return orbits.map(function (o) {
      var fp  = col.filter(ee.Filter.eq('relativeOrbitNumber_start', o)).geometry().dissolve(100);
      var cov = fp.intersection(geom, 100).area(100).divide(boxArea);
      return ee.Feature(null, {day: d, orbit: o, coverage_pct: cov.multiply(100)});
    });
  }).flatten();
  return ee.FeatureCollection(feats).sort('coverage_pct', false);
}

// ---------- 8. RASTER GENERATION FUNCTION ----------
function exportTestRaster(areaName) {
  var geom  = BOX[areaName];
  var ranks = bestPass(CANDIDATES[areaName], geom);
  var best  = ee.Feature(ranks.first());
  print(areaName + ': all candidate passes (coverage % of the box)', ranks);
  print(areaName + ': CHOSEN pass -> copy this day into Cell 2 of Step 6B', best);

  var t0 = ee.Date(best.get('day'));
  var passCol = ALL.filterBounds(geom).filterDate(t0, t0.advance(1, 'day'))
                   .filter(ee.Filter.eq('relativeOrbitNumber_start', best.get('orbit')))
                   .sort('system:time_start');
  var refImg    = ee.Image(passCol.first());   // date and orbit of the pass
  var mosaicImg = passCol.mosaic();            // all frames of this one pass

  var ref = perPixel(seasonRef(refImg, -1, -2, geom), seasonRef(refImg, 1, 2, geom), MED_COUNT, 'VV_median');
  var lkQ = perPixel(yearRef(refImg, -1, geom), yearRef(refImg, 1, geom), P_COUNT, 'VV_p' + RULE.look_pct);
  var geo = geometryFlags(mosaicImg, geom);
  var vv  = mosaicImg.select('VV');

  var valid = vv.mask().gt(0).unmask(0, false);
  var land  = valid.and(perm.not());
  var look  = lkQ.value.lt(RULE.look_db).unmask(0, false);
  var seasonRefAvailable = ref.value.mask().gt(0).unmask(0, false);
  var lookRefAvailable   = lkQ.value.mask().gt(0).unmask(0, false);
  var seasonCountOk      = ref.count.unmask(0, false).gte(RULE.min_ref);
  var lookCountOk        = lkQ.count.unmask(0, false).gte(RULE.min_ref);
  var noref  = seasonCountOk.not().or(lookCountOk.not()).or(seasonRefAvailable.not()).or(lookRefAvailable.not()).unmask(1, false);
  var unobs  = land.and(look.or(noref));
  var judged = land.and(unobs.not());
  var dropOk = vv.subtract(ref.value).lte(-RULE.drop_db).unmask(0, false);

  var flood_sth  = land.and(vv.lt(STH_DB)).unmask(0, false);
  var flood_dccd = judged.and(vv.lt(RULE.water_vv_db)).and(dropOk).and(handOk).unmask(0, false);

  var rasterStack = ee.Image.cat([
    vv.rename('vv_db').toFloat(),
    flood_sth.rename('flood_sth').toFloat(),
    flood_dccd.rename('flood_dccd').toFloat(),
    perm.rename('perm_water').toFloat(),
    unobs.rename('unobs').toFloat(),
    look.rename('look_alike').toFloat(),
    geo.shadow.rename('shadow').toFloat(),
    geo.layover.rename('layover').toFloat(),
    snowFlag(refImg).rename('snow').toFloat()
  ]).clip(geom);

  Map.addLayer(vv.clip(geom), {min: -25, max: 0}, areaName + ' VV (dB)');
  Map.addLayer(flood_dccd.updateMask(flood_dccd).clip(geom), {palette: ['red']}, areaName + ' DCCD flood', false);

  Export.image.toDrive({
    image: rasterStack,
    description: PREFIX + '_' + areaName,
    folder: DRIVE_FOLDER,
    scale: 20,
    crs: 'EPSG:4326',
    region: geom,
    maxPixels: 1e9
  });
}

// ---------- 9. RUN ----------
Map.centerObject(BOX[RUN[0]], 8);
RUN.forEach(exportTestRaster);
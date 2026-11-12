// =====================================================================================
// Paper 1 · Step 10 · GEE script A — Sentinel-2 image near each radar pass (candidates of 10.1)
// Input : table asset uploaded from P1_S10_1_candidates_for_GEE.csv  (X column = lon, Y column = lat)
// Output: Google Drive -> P1_S10_S2_availability.csv
//         put it in  D:\PhD Work\Data\Data required for Python analysis\Step-10\   then run 10.2
// For every candidate: all Sentinel-2 L1C images within +-2 days that touch the 300 m circle;
// keeps the best one = covers >= 95 % of the circle, then fewest days apart, then least cloud.
// No flood rule is applied here (rule 3652 results come from Step 8; this only finds optical images).
// =====================================================================================
var CAND_ASSET = 'projects/ee-ranazabihullah/assets/P1_S10_1_candidates_for_GEE';   // <-- change to your asset path
var RADIUS_M   = 300;   // same circle as Step 8
var MAX_DAYS   = 2;     // radar day +- 2 days
var CLOUD_PROB = 40;    // a pixel is cloudy if the Sentinel-2 cloud probability is >= 40

var cand = ee.FeatureCollection(CAND_ASSET);
var S2   = ee.ImageCollection('COPERNICUS/S2_HARMONIZED');
var S2CP = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY');
var sameIndex = ee.Filter.equals({leftField: 'system:index', rightField: 'system:index'});

function check(f) {
  var day    = ee.Date.fromYMD(ee.Number(f.get('y')), ee.Number(f.get('m')), ee.Number(f.get('d')));
  var circle = f.geometry().buffer(RADIUS_M);
  var t0 = day.advance(-MAX_DAYS, 'day');
  var t1 = day.advance(MAX_DAYS + 1, 'day');
  var imgs = ee.ImageCollection(ee.Join.saveFirst('cp').apply({
    primary:   S2.filterBounds(circle).filterDate(t0, t1),
    secondary: S2CP.filterBounds(circle).filterDate(t0, t1),
    condition: sameIndex
  })).filter(ee.Filter.notNull(['cp']));

  var rows = imgs.toList(100).map(function (x) {
    var img   = ee.Image(x);
    var b2    = img.select('B2');
    var valid = ee.Image(0).where(b2.mask().gt(0), 1).rename('v');                       // 1 = image has data
    var cloud = ee.Image(0).where(ee.Image(img.get('cp')).select('probability').gte(CLOUD_PROB), 1)
                  .multiply(valid).rename('c');                                            // 1 = cloudy pixel with data
    var r = valid.addBands(cloud).reduceRegion({
      reducer: ee.Reducer.mean(), geometry: circle, crs: b2.projection(), scale: 10, maxPixels: 1e8});
    var v = ee.Number(r.get('v'));
    var c = ee.Number(r.get('c')).divide(v.max(1e-6));
    var s2day = img.date().format('YYYY-MM-dd');
    var apart = ee.Date(s2day).difference(day, 'day').abs().round();
    return ee.Feature(null, {
      s2_id: img.get('system:index'), s2_date: s2day, days_apart: apart,
      cloud_pct: c.multiply(100), valid_pct: v.multiply(100),
      score: ee.Number(v.lt(0.95)).multiply(1000).add(apart.multiply(100)).add(c.multiply(100))
    });
  });

  var fc   = ee.FeatureCollection(rows);
  var none = ee.Feature(null, {s2_id: 'none', s2_date: '', days_apart: 99, cloud_pct: 100, valid_pct: 0});
  var best = ee.Feature(ee.Algorithms.If(fc.size().gt(0), fc.sort('score').first(), none));
  return ee.Feature(null, {
    cand_id: f.get('cand_id'), s2_id: best.get('s2_id'), s2_date: best.get('s2_date'),
    days_apart: best.get('days_apart'), cloud_pct: best.get('cloud_pct'), valid_pct: best.get('valid_pct')
  });
}

var out = cand.map(check);
print('Candidates in the asset:', cand.size());
print('Check of the first 3 results:', out.limit(3));

Export.table.toDrive({
  collection: out,
  description: 'P1_S10_S2_availability',
  fileNamePrefix: 'P1_S10_S2_availability',
  fileFormat: 'CSV',
  selectors: ['cand_id', 's2_id', 's2_date', 'days_apart', 'cloud_pct', 'valid_pct']
});

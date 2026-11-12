// =====================================================================================
// Paper 1 · Step 10 · GEE script B — blind labelling viewer
// Input : table asset uploaded from P1_S10_2_sample_for_GEE.csv (X column = lon, Y column = lat).
//         Upload the CSV exactly as Python wrote it (do not open and save it in Excel first).
// LEFT map : the Sentinel-2 image chosen in 10.2 (date shown).
// RIGHT map: the NORMAL view for this time of year = cloud-free median of Sentinel-2,
//            same +-20 days of the year, all other years 2017-2025.
// The radar result is NOT shown (blind). Write your label for each sample_id in
// P1_S10_2_label_sheet.csv:  water / no water / unclear
// =====================================================================================
var SAMPLE_ASSET = 'projects/ee-ranazabihullah/assets/P1_S10_2_sample_for_GEE';   // <-- change to your asset path
var RADIUS_M    = 300;   // same circle as Step 8
var CLOUD_PROB  = 40;    // cloudy pixel if cloud probability >= 40 (normal view)
var HALF_WINDOW = 20;    // normal view: +-20 days of the year

var sample = ee.FeatureCollection(SAMPLE_ASSET).sort('sample_id');
var list   = sample.toList(5000);
var S2     = ee.ImageCollection('COPERNICUS/S2_HARMONIZED');
var S2CP   = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY');
var sameIndex = ee.Filter.equals({leftField: 'system:index', rightField: 'system:index'});
var perm     = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence').gte(80).unmask(0);
var permLine = perm.focalMax(1).neq(perm.focalMin(1)).selfMask();          // outline of permanent water (JRC >= 80 %)

var TC   = {bands: ['B4', 'B3', 'B2'], min: 0, max: 3000, gamma: 1.2};
var SW   = {bands: ['B12', 'B8', 'B4'], min: 0, max: 4000};
var NDWI = {min: -0.3, max: 0.5, palette: ['ffffff', '0000ff']};

var RULES =
  'HOW TO LABEL (compare LEFT with RIGHT)\n\n' +
  'water: NEW standing water (open water, flooded fields or ground) in the LEFT image that is not water in the ' +
  'RIGHT normal view, covering about one tenth or more of the land inside the magenta circle.\n\n' +
  'no water: no new standing water. Water that is also in the normal view (river channel, ponds, rice fields, ' +
  'water inside the yellow permanent-water line) does not count.\n\n' +
  'unclear: cloud, cloud shadow or haze over the circle, or you cannot decide.\n\n' +
  'Tips: water is dark blue to black in true colour and black in the SWIR view. Snow and ice are white in ' +
  'true colour and bright cyan in the SWIR view (not water). Mountain shadow is dark but follows the slopes. ' +
  'Wet soil without standing water = no water. Use the layer list (top right of each map) to switch views.';

// ---- layout: control panel + two linked maps ----
var leftMap  = ui.Map();
var rightMap = ui.Map();
ui.Map.Linker([leftMap, rightMap]);
[leftMap, rightMap].forEach(function (m) {
  m.setControlVisibility({layerList: false, zoomControl: false, scaleControl: false, mapTypeControl: false, fullscreenControl: false});
});
var leftTitle  = ui.Label('', {position: 'bottom-left', fontWeight: 'bold'});
var rightTitle = ui.Label('', {position: 'bottom-left', fontWeight: 'bold'});
leftMap.add(leftTitle);
rightMap.add(rightTitle);
var split = ui.SplitPanel({firstPanel: leftMap, secondPanel: rightMap, orientation: 'horizontal', style: {stretch: 'both'}});

var i = 0;
var N = 0;
var info = ui.Label('Loading ...', {fontWeight: 'bold', fontSize: '18px'});
var prevBtn = ui.Button('< Previous', function () { if (i > 0) { show(i - 1); } });
var nextBtn = ui.Button('Next >', function () { if (i < N - 1) { show(i + 1); } });
var jump = ui.Textbox({placeholder: 'go to sample_id, press Enter', onChange: function (t) {
  var k = parseInt(t, 10);
  if (k >= 1 && k <= N) { show(k - 1); }
}});
var panel = ui.Panel([
  ui.Label('P1 Step 10 - blind labelling', {fontWeight: 'bold', fontSize: '20px'}),
  info,
  ui.Panel([prevBtn, nextBtn], ui.Panel.Layout.flow('horizontal')),
  jump,
  ui.Label(RULES, {whiteSpace: 'pre-wrap', fontSize: '13px'})
], ui.Panel.Layout.flow('vertical'), {width: '360px'});

ui.root.clear();
ui.root.setLayout(ui.Panel.Layout.flow('horizontal'));
ui.root.add(panel);
ui.root.add(split);

// ---- show one point ----
function show(k) {
  i = k;
  info.setValue('Loading point ' + (k + 1) + ' of ' + N + ' ...');
  ee.Feature(list.get(k)).evaluate(function (f, err) {
    if (err) { info.setValue('Error: ' + err); return; }
    var p   = f.properties;
    var pt  = ee.Geometry.Point([p.lon, p.lat]);
    var ring = ee.Image().byte().paint(ee.FeatureCollection([ee.Feature(pt.buffer(RADIUS_M))]), 1, 3);
    var img = ee.Image('COPERNICUS/S2_HARMONIZED/' + p.s2_id);

    var dt  = new Date(String(p.s2_date) + 'T00:00:00Z');
    var yr  = dt.getUTCFullYear();
    var doy = Math.round((dt.getTime() - Date.UTC(yr, 0, 1)) / 86400000) + 1;
    var a = doy - HALF_WINDOW;
    var b = doy + HALF_WINDOW;
    if (a < 1) { a += 365; }
    if (b > 365) { b -= 365; }
    var others = S2.filterBounds(pt)
      .filter(ee.Filter.calendarRange(2017, 2025, 'year'))
      .filter(ee.Filter.calendarRange(yr, yr, 'year').not())
      .filter(ee.Filter.calendarRange(a, b, 'day_of_year'));
    var normal = ee.ImageCollection(ee.Join.saveFirst('cp').apply({
      primary: others, secondary: S2CP.filterBounds(pt), condition: sameIndex
    })).filter(ee.Filter.notNull(['cp'])).map(function (im) {
      return im.updateMask(ee.Image(im.get('cp')).select('probability').lt(CLOUD_PROB));
    }).median();

    leftMap.layers().reset([
      ui.Map.Layer(img, TC, 'True colour', false),
      ui.Map.Layer(img, SW, 'SWIR false colour (water black, snow cyan)', false),
      ui.Map.Layer(img.normalizedDifference(['B3', 'B8']), NDWI, 'NDWI (water blue; snow also blue)', false),
      ui.Map.Layer(permLine, {palette: ['ffff00']}, 'Permanent water line (JRC >= 80 %)', false),
      ui.Map.Layer(ring, {palette: ['ff00ff']}, '300 m circle', false)
    ]);
    rightMap.layers().reset([
      ui.Map.Layer(normal, TC, 'Normal view - true colour', false),
      ui.Map.Layer(normal, SW, 'Normal view - SWIR false colour', false),
      ui.Map.Layer(
        normal.normalizedDifference(['B3', 'B8']),false
        NDWI,
        'Normal view - NDWI (water blue; snow also blue)',
        false
        ),
      ui.Map.Layer(permLine, {palette: ['ffff00']}, 'Permanent water line (JRC >= 80 %)', false),
      ui.Map.Layer(ring, {palette: ['ff00ff']}, '300 m circle', false)
    ]);
    leftMap.centerObject(pt, 15);
    leftTitle.setValue('Sentinel-2 on ' + p.s2_date);
    rightTitle.setValue('Normal view: day ' + doy + ' +-' + HALF_WINDOW + ', other years');
    info.setValue('sample_id ' + p.sample_id + '   (' + (k + 1) + ' of ' + N + ')');
  });
}

sample.size().evaluate(function (n) {
  N = n;
  if (N > 0) { show(0); } else { info.setValue('The sample asset is empty - check SAMPLE_ASSET.'); }
});

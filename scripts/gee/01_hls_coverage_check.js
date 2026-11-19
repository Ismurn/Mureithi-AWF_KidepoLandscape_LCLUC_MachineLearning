/* 
-  KIDEPO LANDSCAPE — LAND COVER / LAND USE CHANGE (LCLUC)
-  GEE Script 01 — HLS coverage check
- -----------------------------------------------------------------------------
-  PURPOSE
-  -------
-  Decide the TWO comparison dates with evidence instead of guesswork. For every
-  candidate year this script evaluates clear HLS observation coverage inside the
-  PA union and separately for the Uganda and South Sudan portions.
 *
-  ASSETS USED (in users/imnduta1/KidepoAOI-LULC-ML/):
-      kidepo_pa_union_shp
-      kidepo_pa_parts_shp
 *
-  Author: Isaac Mureithi
-  */


/* 
-  PHASE 0 — CONFIGURATION
-  */

//  0.1 Assets 
var ASSET_ROOT = 'users/imnduta1/KidepoAOI-LULC-ML/';
var aoi_union  = ee.FeatureCollection(ASSET_ROOT + 'kidepo_pa_union_shp');
var aoi_parts  = ee.FeatureCollection(ASSET_ROOT + 'kidepo_pa_parts_shp');

var aoi = aoi_union.geometry();

//  0.2 Candidate years 
var YEAR_START = 2016;
var YEAR_END   = 2025;

//  0.3 Seasonal windows 
var WINDOWS = {
  dry_early: {startMonth: 12, startDay: 1, lengthDays: 62},  // 01 Dec -> 31 Jan
  dry_late:  {startMonth:  2, startDay: 1, lengthDays: 59}   // 01 Feb -> 31 Mar
};

var WINDOW_NAME = 'dry_late';          // <<< switch to 'dry_early' and re-run to compare
var WINDOW      = WINDOWS[WINDOW_NAME];

//  0.4 Quality thresholds 
var MAX_SCENE_CLOUD = 60;   // %; drop whole scenes cloudier than this before masking
var MIN_CLEAR_OBS   = 4;    // clear observations per pixel for a usable median
var SCALE_STATS     = 270;  // m; coarse screening scale to avoid capacity/concurrency limits


/* 
-  PHASE 1 — HELPERS
-  */

function windowFor(year) {
  var start = ee.Date.fromYMD(year, WINDOW.startMonth, WINDOW.startDay);
  var end   = start.advance(WINDOW.lengthDays, 'day');
  return {start: start, end: end};
}

/**
- Fmask decoding:
-   bit 0: cirrus | bit 1: cloud | bit 2: adjacent-to-cloud | bit 3: cloud shadow
- Water (bit 5) is kept unmasked.
 */
function clearMask(img) {
  var fmask    = img.select('Fmask');
  var bitClear = fmask.bitwiseAnd(1 << 0).eq(0)
    .and(fmask.bitwiseAnd(1 << 1).eq(0))
    .and(fmask.bitwiseAnd(1 << 2).eq(0))
    .and(fmask.bitwiseAnd(1 << 3).eq(0));
  return bitClear.selfMask().rename('clear');
}

/**
- Build a lightweight date- and AOI-filtered Fmask collection for one product.
- Filtering by date first avoids mapping across the entire multi-year archive.
 */
function prepWindowProduct(collectionId, productTag, start, end) {
  return ee.ImageCollection(collectionId)
    .filterBounds(aoi)
    .filterDate(start, end)
    .filterMetadata('CLOUD_COVERAGE', 'less_than', MAX_SCENE_CLOUD)
    .map(function (img) {
      return clearMask(img)
        .set('hls_product', productTag);
    });
}

function hlsWindowCol(start, end) {
  var l30 = prepWindowProduct('NASA/HLS/HLSL30/v002', 'L30', start, end);
  var s30 = prepWindowProduct('NASA/HLS/HLSS30/v002', 'S30', start, end);
  return l30.merge(s30);
}

function clearObsFromCol(col) {
  return ee.Image(
      ee.Algorithms.If(col.size().gt(0),
                       col.select('clear').count(),
                       ee.Image.constant(0))
    )
    .unmask(0)
    .rename('clear_obs')
    .clip(aoi);
}


/* 
-  PHASE 2 — PER-YEAR COVERAGE STATISTICS
-  */

var SUBAREAS = {
  all:         aoi,
  uganda:      aoi_parts.filter(ee.Filter.eq('country', 'Uganda')).geometry(),
  south_sudan: aoi_parts.filter(ee.Filter.eq('country', 'South Sudan')).geometry()
};

var years = ee.List.sequence(YEAR_START, YEAR_END);

function statsForYear(year) {
  year    = ee.Number(year);
  var win = windowFor(year);
  var col = hlsWindowCol(win.start, win.end);

  var clearObs = clearObsFromCol(col);
  var enough   = clearObs.gte(MIN_CLEAR_OBS).rename('enough');
  var metrics  = clearObs.addBands(enough);

  function meanDict(geom) {
    return ee.Dictionary(
      metrics.reduceRegion({
        reducer:    ee.Reducer.mean(),
        geometry:   geom,
        scale:      SCALE_STATS,
        maxPixels:  1e10,
        bestEffort: true
      })
    );
  }

  var sAll = meanDict(SUBAREAS.all);
  var sUga = meanDict(SUBAREAS.uganda);
  var sSsd = meanDict(SUBAREAS.south_sudan);

  return ee.Feature(null, {
    year:         year,
    window:       WINDOW_NAME,
    date_start:   win.start.format('YYYY-MM-dd'),
    date_end:     win.end.format('YYYY-MM-dd'),
    n_scenes:     col.size(),
    n_scenes_L30: col.filterMetadata('hls_product', 'equals', 'L30').size(),
    n_scenes_S30: col.filterMetadata('hls_product', 'equals', 'S30').size(),
    mean_obs_all: sAll.get('clear_obs'),
    mean_obs_uga: sUga.get('clear_obs'),
    mean_obs_ssd: sSsd.get('clear_obs'),
    pct_ok_all:   ee.Number(sAll.get('enough')).multiply(100),
    pct_ok_uga:   ee.Number(sUga.get('enough')).multiply(100),
    pct_ok_ssd:   ee.Number(sSsd.get('enough')).multiply(100)
  });
}

var coverage = ee.FeatureCollection(years.map(statsForYear));


/* 
-  PHASE 3 — CHARTS (UNCOMMENT ONE AT A TIME IF RUNNING INTERACTIVELY)
-  */
// Running all 3 charts at once often triggers "Too many concurrent aggregations"
// or "Earth Engine capacity exceeded". Prefer the Drive export in Phase 5, or
// uncomment ONE chart block below per run.

// print(ui.Chart.feature.byFeature(coverage, 'year',
//         ['pct_ok_all', 'pct_ok_uga', 'pct_ok_ssd'])
//   .setChartType('LineChart')
//   .setOptions({
//     title:  '% of AOI with >= ' + MIN_CLEAR_OBS + ' clear HLS observations (' + WINDOW_NAME + ')',
//     hAxis:  {title: 'Year', format: '####'},
//     vAxis:  {title: '% of area', viewWindow: {min: 0, max: 100}},
//     series: {0: {color: '#333333'}, 1: {color: '#4C9F70'}, 2: {color: '#6A8EAE'}},
//     pointSize: 4, lineWidth: 2
//   }));


/* 
-  PHASE 4 — MAPS: WHERE the gaps are (TWO CANDIDATE YEARS ONLY)
-  */

var YEAR_A = 2019;   // <<< edit to candidate year 1
var YEAR_B = 2024;   // <<< edit to candidate year 2

function clearObsImage(year) {
  var win = windowFor(year);
  return clearObsFromCol(hlsWindowCol(win.start, win.end));
}

var visObs = {min: 0, max: 20, palette: ['#67001f', '#f7f7f7', '#1b7837']};

Map.centerObject(aoi, 8);
Map.addLayer(clearObsImage(YEAR_A), visObs, 'Clear obs ' + YEAR_A + ' (' + WINDOW_NAME + ')');
Map.addLayer(clearObsImage(YEAR_B), visObs, 'Clear obs ' + YEAR_B + ' (' + WINDOW_NAME + ')');

Map.addLayer(clearObsImage(YEAR_A).lt(MIN_CLEAR_OBS).selfMask(),
             {palette: ['#d7301f']}, 'Under-sampled ' + YEAR_A, false);
Map.addLayer(clearObsImage(YEAR_B).lt(MIN_CLEAR_OBS).selfMask(),
             {palette: ['#d7301f']}, 'Under-sampled ' + YEAR_B, false);

Map.addLayer(aoi_parts.style({color: '#000000', fillColor: '00000000', width: 1}),
             {}, 'PA polygons');
Map.addLayer(ee.Image().paint(aoi_union, 0, 2), {palette: ['#E07B39']}, 'PA union (buffered)');


/* 
-  PHASE 5 — EXPORT COVERAGE TABLE TO DRIVE (RECOMMENDED OVER CONSOLE CHARTS)
-  */

Export.table.toDrive({
  collection:  coverage,
  description: 'kidepo_hls_coverage_' + WINDOW_NAME + '_' + YEAR_START + '_' + YEAR_END,
  folder:      'kidepo_lcluc',
  fileFormat:  'CSV'
});
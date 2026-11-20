/*
 *  KIDEPO LANDSCAPE — LAND COVER / LAND USE CHANGE (LCLUC)
 *  GEE Script 01 — HLS coverage check
 *
 *  PURPOSE
 *  Choose the two comparison dates on evidence. For every candidate year AND both
 *  candidate seasonal windows, count the per-pixel CLEAR HLS observations inside
 *  the PA union, reported separately for the Uganda and South Sudan portions.
 *
 *  WHAT CHANGED IN THIS VERSION
 *   - Both windows (dry_early and dry_late) are evaluated in ONE run and land in
 *     ONE CSV, so the cloud-vs-fire trade-off can be read off a single table
 *     instead of two runs that have to be reconciled by hand.
 *   - Window ends are pinned to explicit calendar dates rather than a day count,
 *     which removes the one-day leap-year drift (Feb 1 + 59 days ends 31 Mar in a
 *     leap year and 1 Apr otherwise). Two dates must share an identical window.
 *   - Candidate years run to 2026 for dry_late (Feb-Mar 2026 has closed) and to
 *     2025 for dry_early (Dec 2025 - Jan 2026 has closed; Dec 2026 has not).
 *
 *  The date-filter-before-map structure, the Fmask-only collection, the single
 *  two-band reduceRegion and the coarse screening scale are kept as-is: they are
 *  what makes this run inside Earth Engine's interactive limits.
 *
 *  ASSETS USED (in users/imnduta1/KidepoAOI-LULC-ML/):
 *      kidepo_pa_union_shp    dissolved PA union
 *      kidepo_pa_parts_shp    individual PAs, carrying the `country` attribute
 *
 *  OUTPUT: one CSV in Drive/kidepo_lcluc, with a `window` column.
 *
 *  Author: Isaac Mureithi - IMureithi@awf.org
 */


// =============================================================================
// PHASE 0 — CONFIGURATION
// =============================================================================

// 0.1 Assets
var ASSET_ROOT = 'users/imnduta1/KidepoAOI-LULC-ML/';
var aoi_union  = ee.FeatureCollection(ASSET_ROOT + 'kidepo_pa_union_shp');
var aoi_parts  = ee.FeatureCollection(ASSET_ROOT + 'kidepo_pa_parts_shp');

var aoi = aoi_union.geometry();

// 0.2 Candidate years
var YEAR_START = 2016;   // first year both HLS sensors contribute

// 0.3 Seasonal windows
// Karamoja is unimodal: long rains roughly Apr-Oct, sharp dry season Dec-Mar.
// The two windows trade cloud against fire:
//
//   dry_early  1 Dec -> 31 Jan   Crosses the new year, so a row labelled 2019
//                                means Dec 2019 -> Jan 2020.
//   dry_late   1 Feb -> 31 Mar   Same calendar year.
//
// MEASURED, not assumed (runs of this script and 01b, 2016-2026). Both of the
// expectations this script was originally written around turned out to be wrong,
// which is the reason for measuring:
//
//   - dry_early is NOT the cloudier window. It is the CLEARER one: it reaches
//     >= MIN_CLEAR_OBS over ~100% of the AOI in 9 of 10 years (the exception is
//     Dec 2019 - Jan 2020 at 62%), against 95-100% in only 6 of 11 for dry_late,
//     and it carries more clear looks per pixel. Smoke from the Dec-Jan fires
//     plausibly costs dry_late some observations by tripping Fmask.
//   - dry_late is NOT the peak burning window. December is, followed by January
//     and November. By February the burning is largely over: 1-7% of the AOI
//     burns DURING dry_late against 12-37% during dry_early.
//
// The scars persist either way, so both windows image a landscape that is
// 18-44% burnt. See 01b and the window note in 02 for which window that argues
// for, and why.
//
// `endYearOffset` is 1 where the window crosses into the following year.
// `lastStartYear` stops the sequence before a window that has not closed yet.
var WINDOW_DEFS = {
  dry_early: {startMonth: 12, startDay: 1, endMonth: 1, endDay: 31,
              endYearOffset: 1, lastStartYear: 2025},
  dry_late:  {startMonth:  2, startDay: 1, endMonth: 3, endDay: 31,
              endYearOffset: 0, lastStartYear: 2026}
};

var WINDOW_NAMES = ['dry_early', 'dry_late'];

// 0.4 Quality thresholds
var MAX_SCENE_CLOUD = 60;   // %; drop whole scenes cloudier than this before masking
var MIN_CLEAR_OBS   = 4;    // clear observations per pixel for a usable median
var SCALE_STATS     = 270;  // m; coarse screening scale to stay inside EE limits


// =============================================================================
// PHASE 1 — HELPERS
// =============================================================================

/**
 * 1.1 Date range for one window in one year.
 * The end is built from an explicit calendar date and then advanced one day,
 * because filterDate's end is EXCLUSIVE. This makes the window the same length
 * in leap and non-leap years, which a fixed day-count does not.
 */
function windowFor(windowName, year) {
  var w     = WINDOW_DEFS[windowName];
  var start = ee.Date.fromYMD(year, w.startMonth, w.startDay);
  var end   = ee.Date
    .fromYMD(ee.Number(year).add(w.endYearOffset), w.endMonth, w.endDay)
    .advance(1, 'day');
  return {start: start, end: end};
}

/**
 * 1.2 Fmask decoding.
 * HLS v2 carries one 8-bit Fmask band. Bits, LSB first:
 *   0 cirrus | 1 cloud | 2 adjacent-to-cloud | 3 cloud shadow
 *   4 snow/ice | 5 water | 6-7 aerosol
 * Water (bit 5) is deliberately NOT masked: permanent water is a land-cover
 * class we want to classify, not a defect to remove.
 *
 * selfMask() matters: `clear` is a 0/1 band, and an unmasked version would make
 * ImageCollection.count() count every observation in the footprint, cloud
 * included. Masking the zeros away makes count() mean "clear looks".
 */
function clearMask(img) {
  var fmask = img.select('Fmask');
  return fmask.bitwiseAnd(1 << 0).eq(0)
    .and(fmask.bitwiseAnd(1 << 1).eq(0))
    .and(fmask.bitwiseAnd(1 << 2).eq(0))
    .and(fmask.bitwiseAnd(1 << 3).eq(0))
    .selfMask()
    .rename('clear');
}

/**
 * 1.3 A date- and AOI-filtered mask-only collection for one HLS product.
 * Filtering by date BEFORE map() is what keeps this cheap: there is no reason to
 * touch the whole multi-year archive to describe one 60-day window.
 */
function prepWindowProduct(collectionId, productTag, start, end) {
  return ee.ImageCollection(collectionId)
    .filterBounds(aoi)
    .filterDate(start, end)
    .filterMetadata('CLOUD_COVERAGE', 'less_than', MAX_SCENE_CLOUD)
    .map(function (img) {
      return clearMask(img).set('hls_product', productTag);
    });
}

function hlsWindowCol(start, end) {
  return prepWindowProduct('NASA/HLS/HLSL30/v002', 'L30', start, end)
    .merge(prepWindowProduct('NASA/HLS/HLSS30/v002', 'S30', start, end));
}

/**
 * 1.4 Per-pixel clear-observation count.
 * The ee.Algorithms.If guard matters: count() on an EMPTY collection returns an
 * image with no bands, and rename() then throws. A window outside HLS coverage
 * must report zero, not break the run.
 */
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


// =============================================================================
// PHASE 2 — PER-YEAR, PER-WINDOW COVERAGE STATISTICS
// =============================================================================

// 2.1 Sub-areas. The transboundary pieces are the ones at risk of being
// under-sampled, so they get their own columns.
// Note: the ~26 km2 where Kidepo GR (SS) overlaps Kidepo Valley NP (UG) falls in
// BOTH sub-areas, so uganda + south_sudan will not sum exactly to all.
var SUBAREAS = {
  all:         aoi,
  uganda:      aoi_parts.filter(ee.Filter.eq('country', 'Uganda')).geometry(),
  south_sudan: aoi_parts.filter(ee.Filter.eq('country', 'South Sudan')).geometry()
};

/**
 * 2.2 One row: one window, one year.
 * `windowName` stays a CLIENT-side string so WINDOW_DEFS can be indexed with it;
 * only the year varies server-side inside ee.List.map.
 */
function statsFor(windowName, year) {
  year     = ee.Number(year);
  var win  = windowFor(windowName, year);
  var col  = hlsWindowCol(win.start, win.end);

  var clearObs = clearObsFromCol(col);
  var enough   = clearObs.gte(MIN_CLEAR_OBS).rename('enough');
  var metrics  = clearObs.addBands(enough);

  // One reduction per sub-area returning BOTH statistics, rather than two.
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
    window:       windowName,
    year:         year,
    date_start:   win.start.format('YYYY-MM-dd'),
    date_end:     win.end.advance(-1, 'day').format('YYYY-MM-dd'),  // inclusive, for reading
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

// 2.3 Build every row. The window loop is a plain JS loop (two iterations) so
// that windowName stays client-side; only the year sequence is server-side.
var blocks = [];
WINDOW_NAMES.forEach(function (nm) {
  var yrs = ee.List.sequence(YEAR_START, WINDOW_DEFS[nm].lastStartYear);
  blocks.push(ee.FeatureCollection(
    yrs.map(function (y) { return statsFor(nm, y); })
  ));
});

var coverage = ee.FeatureCollection(blocks).flatten();


// =============================================================================
// PHASE 3 — CHARTS (OPTIONAL; UNCOMMENT ONE AT A TIME)
// =============================================================================
// Running several charts at once tends to trigger "Too many concurrent
// aggregations" or "Earth Engine capacity exceeded". The Drive export in Phase 5
// is the reliable route; these are only for a quick look.

// print(ui.Chart.feature.groups(coverage, 'year', 'pct_ok_all', 'window')
//   .setChartType('LineChart')
//   .setOptions({
//     title: '% of AOI with >= ' + MIN_CLEAR_OBS + ' clear HLS observations',
//     hAxis: {title: 'Year', format: '####'},
//     vAxis: {title: '% of area', viewWindow: {min: 0, max: 100}},
//     pointSize: 4, lineWidth: 2
//   }));


// =============================================================================
// PHASE 4 — MAPS: WHERE the gaps are, for the two candidate years
// =============================================================================
// An AOI-wide percentage hides spatial structure: persistent cloud sits on the
// Imatong and Agoro-Agu highlands, so a year can look fine overall while one
// massif is unusable. Check the red layers before committing to a pair.

var MAP_WINDOW = 'dry_late';
var YEAR_A = 2019;
var YEAR_B = 2024;

function clearObsImage(windowName, year) {
  var win = windowFor(windowName, year);
  return clearObsFromCol(hlsWindowCol(win.start, win.end));
}

var obsA = clearObsImage(MAP_WINDOW, YEAR_A);
var obsB = clearObsImage(MAP_WINDOW, YEAR_B);

var visObs = {min: 0, max: 20, palette: ['#67001f', '#f7f7f7', '#1b7837']};

Map.centerObject(aoi, 8);
Map.addLayer(obsA, visObs, 'Clear obs ' + YEAR_A + ' (' + MAP_WINDOW + ')');
Map.addLayer(obsB, visObs, 'Clear obs ' + YEAR_B + ' (' + MAP_WINDOW + ')');

Map.addLayer(obsA.lt(MIN_CLEAR_OBS).selfMask(), {palette: ['#d7301f']},
             'Under-sampled ' + YEAR_A, false);
Map.addLayer(obsB.lt(MIN_CLEAR_OBS).selfMask(), {palette: ['#d7301f']},
             'Under-sampled ' + YEAR_B, false);

// Pixels short of the threshold in EITHER year — these are the pixels where the
// change map would be weakest, and they are a reporting caveat either way.
Map.addLayer(obsA.lt(MIN_CLEAR_OBS).or(obsB.lt(MIN_CLEAR_OBS)).selfMask(),
             {palette: ['#6a3d9a']}, 'Under-sampled in either year', false);

Map.addLayer(aoi_parts.style({color: '#000000', fillColor: '00000000', width: 1}),
             {}, 'PA polygons');
Map.addLayer(ee.Image().paint(aoi_union, 0, 2), {palette: ['#E07B39']}, 'PA union');


// =============================================================================
// PHASE 5 — EXPORT THE COVERAGE TABLE
// =============================================================================

Export.table.toDrive({
  collection:  coverage,
  description: 'kidepo_hls_coverage_bothwindows_' + YEAR_START + '_2026',
  folder:      'kidepo_lcluc',
  fileFormat:  'CSV',
  selectors:   ['window', 'year', 'date_start', 'date_end',
                'n_scenes', 'n_scenes_L30', 'n_scenes_S30',
                'mean_obs_all', 'mean_obs_uga', 'mean_obs_ssd',
                'pct_ok_all', 'pct_ok_uga', 'pct_ok_ssd']
});


/*
 *  HOW TO READ THE RESULT
 *  1. Shortlist years where pct_ok_all is high AND pct_ok_ssd is not far behind.
 *     A year carried by the Uganda side alone is not a usable date.
 *  2. Prefer a comfortable margin over MIN_CLEAR_OBS, not a year that just clears it.
 *  3. Match the two years on mean_obs as well as pct_ok. A median built from 8
 *     looks against one built from 14 can manufacture apparent change.
 *  4. Separate the dates by 5+ years for a land-cover signal in semi-arid rangeland.
 *  5. Compare the two windows in the `window` column, then read
 *     01b_burn_season_profile.js, which quantifies how much of the AOI is burnt
 *     inside each window. That is the evidence for the cloud-vs-fire decision.
 */

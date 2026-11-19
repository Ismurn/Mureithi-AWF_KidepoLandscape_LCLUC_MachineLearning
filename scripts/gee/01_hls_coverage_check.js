/* 
 *  KIDEPO LANDSCAPE — LAND COVER / LAND USE CHANGE (LCLUC)
 *  GEE Script 01 — HLS coverage check
 * 
 *  PURPOSE
 *  
 *  Decide the TWO comparison dates with evidence instead of guesswork. For every
 *  candidate year this script counts, inside the PA union:
 *
 *    - how many HLS scenes intersect the AOI in the chosen seasonal window
 *    - how many of those are usable (scene cloud cover below a threshold)
 *    - the per-pixel count of CLEAR observations after Fmask masking, which is the
 *      number that actually governs whether a median composite is trustworthy
 *    - the % of the AOI reaching a minimum clear-observation count
 *    - the same, reported SEPARATELY for the Uganda and South Sudan portions,
 *      because the transboundary sub-areas are the ones at risk of under-sampling
 *
 *  It then charts clear-observation counts by year so two well-sampled years that
 *  share a phenological window can be picked.
 *
 *  WHY HLS: NASA/HLS/HLSL30 (Landsat 8/9, 30 m) + HLSS30 (Sentinel-2, 30 m) are
 *  harmonised to a common grid and surface reflectance, so they can be stacked.
 *  That density matters in cloud-prone Karamoja. HLSL30 runs from ~2013, HLSS30
 *  from ~2015, so 2016 onward is where BOTH sensors contribute.
 *
 *  RUN THIS FIRST, before any compositing or export.
 *
 *  PREREQUISITES — upload these from outputs/aoi/ (built by scripts/R/01_build_aoi.R):
 *      users/imnduta1/Ismurn/kidepo_pa_union   <- kidepo.pa.union.buffer500m.geojson
 *      users/imnduta1/Ismurn/kidepo_pa_parts   <- kidepo.pa.individual.geojson
 *      users/imnduta1/Ismurn/kidepo_bbox       <- kidepo.bbox.geojson        (optional here)
 *
 *  Save as: users/imnduta1/Ismurn/01_hls_coverage_check
 *  Author: Isaac Mureithi
 *  */


/* 
 *  PHASE 0 — CONFIGURATION
 *  */

//  0.1 Assets 
var ASSET_ROOT = 'users/imnduta1/KidepoAOI-LULC-ML/';
var aoi_union  = ee.FeatureCollection(ASSET_ROOT + 'kidepo_pa_union_shp');
var aoi_parts  = ee.FeatureCollection(ASSET_ROOT + 'kidepo_pa_parts_shp');

var aoi = aoi_union.geometry();

//  0.2 Candidate years 
// 2016 is the first year both HLS sensors contribute. The current year is excluded
// because its dry season may be incomplete at run time.
var YEAR_START = 2016;
var YEAR_END   = 2025;

//  0.3 Seasonal windows 
// Karamoja is unimodal: a long rains roughly Apr-Oct, a sharp dry season Dec-Mar.
// Two windows are evaluated because they trade cloud against fire:
//
//   dry_early  Dec-Jan  vegetation still senescing, FEWER burn scars, more cloud
//   dry_late   Feb-Mar  clearest skies, but PEAK burning — burn scars will read as
//                       bare/barren and can swamp a change signal
//
// CAVEAT (CLAUDE.md #1): whichever window is chosen must be the SAME for both
// comparison dates, or the change map measures phenology and fire, not land cover.
// `dry_early` crosses the new year, which the date arithmetic below handles.
var WINDOWS = {
  dry_early: {startMonth: 12, startDay: 1, lengthDays: 62},  // 01 Dec -> 31 Jan
  dry_late:  {startMonth:  2, startDay: 1, lengthDays: 59}   // 01 Feb -> 31 Mar
};

var WINDOW_NAME = 'dry_late';          // <<< switch to 'dry_early' and re-run to compare
var WINDOW      = WINDOWS[WINDOW_NAME];

//  0.4 Quality thresholds 
var MAX_SCENE_CLOUD = 60;   // %; drop whole scenes cloudier than this before masking
var MIN_CLEAR_OBS   = 4;    // clear observations per pixel for a usable median
var SCALE_STATS     = 90;   // m; stats resolution. 3x the 30 m pixel keeps the
                            // interactive reductions fast. Only the EXPORT needs 30 m.

//  0.5 Bands 
// HLSL30 (Landsat) and HLSS30 (Sentinel-2) use DIFFERENT band names for the same
// wavelengths. Renaming to a common scheme is what makes the two stackable — this
// mapping is the single most error-prone part of using HLS, so it is explicit.
//
//   common :  BLUE   GREEN  RED    NIR    SWIR1  SWIR2
//   HLSL30 :  B2     B3     B4     B5     B6     B7
//   HLSS30 :  B2     B3     B4     B8A    B11    B12     (B8A = narrow NIR, 20 m native)
var BANDS_COMMON = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2'];
var BANDS_L30    = ['B2', 'B3', 'B4', 'B5',  'B6',  'B7'];
var BANDS_S30    = ['B2', 'B3', 'B4', 'B8A', 'B11', 'B12'];


/* 
 *  PHASE 1 — HELPERS
 *  */

/**
 * 1.1 Build the seasonal window for one year.
 * Returns {start, end, label}. For a window that crosses the new year, `year` is
 * the year the window STARTS in, so dry_early 2020 means Dec 2020 -> Jan 2021.
 */
function windowFor(year) {
  var start = ee.Date.fromYMD(year, WINDOW.startMonth, WINDOW.startDay);
  var end   = start.advance(WINDOW.lengthDays, 'day');
  return {start: start, end: end};
}

/**
 * 1.2 Fmask decoding.
 * HLS v2 carries a single 8-bit `Fmask` band. Bits (LSB first):
 *   0 cirrus | 1 cloud | 2 adjacent-to-cloud | 3 cloud shadow
 *   4 snow/ice | 5 water | 6-7 aerosol level
 * A pixel is kept when cloud, adjacent-cloud, shadow and cirrus are all 0.
 * Water is NOT masked here: permanent water is a legitimate land-cover class, and
 * masking it would silently delete it from the classification.
 */
function clearMask(img) {
  var fmask    = img.select('Fmask');
  var bitClear = fmask.bitwiseAnd(1 << 0).eq(0)   // not cirrus
    .and(fmask.bitwiseAnd(1 << 1).eq(0))          // not cloud
    .and(fmask.bitwiseAnd(1 << 2).eq(0))          // not cloud-adjacent
    .and(fmask.bitwiseAnd(1 << 3).eq(0));         // not cloud shadow
  return bitClear.rename('clear');
}

/**
 * 1.3 Harmonise one collection: rename bands, scale reflectance, attach the mask.
 * HLS v2 surface reflectance is stored as int16 with a 0.0001 scale factor.
 */
function prepCollection(collectionId, bandsIn) {
  return ee.ImageCollection(collectionId)
    .filterBounds(aoi)
    .filterMetadata('CLOUD_COVERAGE', 'less_than', MAX_SCENE_CLOUD)
    .map(function (img) {
      var clear = clearMask(img);
      var sr    = img.select(bandsIn, BANDS_COMMON)
                     .multiply(0.0001)
                     .updateMask(clear);
      // selfMask() is essential: `clear` is a 0/1 band, so an UNMASKED version
      // would make ImageCollection.count() count every observation in the footprint,
      // cloud included. Masking the zeros away makes count() mean "clear looks".
      return sr
        .addBands(clear.selfMask())
        .copyProperties(img, ['system:time_start', 'CLOUD_COVERAGE', 'SENSOR'])
        .set('hls_product', collectionId.indexOf('HLSL30') > -1 ? 'L30' : 'S30');
    });
}

// 1.4 The merged, harmonised archive (filtered to the AOI, not yet to a date).
var hls = prepCollection('NASA/HLS/HLSL30/v002', BANDS_L30)
  .merge(prepCollection('NASA/HLS/HLSS30/v002', BANDS_S30));

print('HLS scenes intersecting the AOI, all dates, scene cloud < ' +
      MAX_SCENE_CLOUD + '%:', hls.size());


/**
 * 1.5 Per-pixel count of clear observations in a date range.
 * The ee.Algorithms.If guard matters: ImageCollection.count() on an EMPTY
 * collection returns an image with NO bands, and rename() then throws. A candidate
 * year outside HLS coverage should report zero, not break the whole script.
 */
function clearObsFor(start, end) {
  var col = hls.filterDate(start, end);
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
 *  PHASE 2 — PER-YEAR COVERAGE STATISTICS
 *  */

// 2.1 Sub-areas to report on separately. The transboundary pieces are the ones
// likely to be under-sampled (CLAUDE.md caveat #4), so they get their own rows.
var SUBAREAS = {
  all:         aoi,
  uganda:      aoi_parts.filter(ee.Filter.eq('country', 'Uganda')).geometry(),
  south_sudan: aoi_parts.filter(ee.Filter.eq('country', 'South Sudan')).geometry()
};

var years = ee.List.sequence(YEAR_START, YEAR_END);

/**
 * 2.2 Statistics for one year.
 * clearObs is the per-pixel count of unmasked observations — the quantity that
 * decides whether a median composite is defensible. A high scene count means
 * nothing if every scene is clouded over the same ridge.
 */
function statsForYear(year) {
  year    = ee.Number(year);
  var win = windowFor(year);
  var col = hls.filterDate(win.start, win.end);

  var clearObs = clearObsFor(win.start, win.end);

  // Fraction of each sub-area with at least MIN_CLEAR_OBS clear looks.
  var enough = clearObs.gte(MIN_CLEAR_OBS).rename('enough');

  function fracEnough(geom) {
    return ee.Number(
      enough.reduceRegion({
        reducer:   ee.Reducer.mean(),
        geometry:  geom,
        scale:     SCALE_STATS,
        maxPixels: 1e10,
        bestEffort: true
      }).get('enough')
    );
  }

  function meanObs(geom) {
    return ee.Number(
      clearObs.reduceRegion({
        reducer:   ee.Reducer.mean(),
        geometry:  geom,
        scale:     SCALE_STATS,
        maxPixels: 1e10,
        bestEffort: true
      }).get('clear_obs')
    );
  }

  return ee.Feature(null, {
    year:            year,
    window:          WINDOW_NAME,
    date_start:      win.start.format('YYYY-MM-dd'),
    date_end:        win.end.format('YYYY-MM-dd'),
    n_scenes:        col.size(),
    n_scenes_L30:    col.filterMetadata('hls_product', 'equals', 'L30').size(),
    n_scenes_S30:    col.filterMetadata('hls_product', 'equals', 'S30').size(),
    mean_obs_all:    meanObs(SUBAREAS.all),
    mean_obs_uga:    meanObs(SUBAREAS.uganda),
    mean_obs_ssd:    meanObs(SUBAREAS.south_sudan),
    pct_ok_all:      fracEnough(SUBAREAS.all).multiply(100),
    pct_ok_uga:      fracEnough(SUBAREAS.uganda).multiply(100),
    pct_ok_ssd:      fracEnough(SUBAREAS.south_sudan).multiply(100)
  });
}

var coverage = ee.FeatureCollection(years.map(statsForYear));

// 2.3 The table to read. Sort by pct_ok_all to see the best-sampled years at a glance.
print('Coverage by year — window "' + WINDOW_NAME + '" (' +
      WINDOW.startDay + '/' + WINDOW.startMonth + ' + ' + WINDOW.lengthDays + ' days), ' +
      'pct_ok = % of area with >= ' + MIN_CLEAR_OBS + ' clear observations:',
      coverage);


/* 
 *  PHASE 3 — CHARTS
 *  */

// 3.1 Usable-area fraction by year, split by country portion.
// Read this chart first: pick two years where ALL THREE series are high. A year
// where pct_ok_ssd lags is a year whose South Sudan result would be weakly supported.
print(ui.Chart.feature.byFeature(coverage, 'year',
        ['pct_ok_all', 'pct_ok_uga', 'pct_ok_ssd'])
  .setChartType('LineChart')
  .setOptions({
    title:  '% of AOI with >= ' + MIN_CLEAR_OBS + ' clear HLS observations (' + WINDOW_NAME + ')',
    hAxis:  {title: 'Year', format: '####'},
    vAxis:  {title: '% of area', viewWindow: {min: 0, max: 100}},
    series: {0: {color: '#333333'}, 1: {color: '#4C9F70'}, 2: {color: '#6A8EAE'}},
    pointSize: 4, lineWidth: 2
  }));

// 3.2 Mean clear observations per pixel — the margin above the threshold.
print(ui.Chart.feature.byFeature(coverage, 'year',
        ['mean_obs_all', 'mean_obs_uga', 'mean_obs_ssd'])
  .setChartType('ColumnChart')
  .setOptions({
    title:  'Mean clear observations per pixel (' + WINDOW_NAME + ')',
    hAxis:  {title: 'Year', format: '####'},
    vAxis:  {title: 'Clear observations'},
    series: {0: {color: '#333333'}, 1: {color: '#4C9F70'}, 2: {color: '#6A8EAE'}}
  }));

// 3.3 Scene counts by sensor — shows when Sentinel-2 density arrives.
print(ui.Chart.feature.byFeature(coverage, 'year', ['n_scenes_L30', 'n_scenes_S30'])
  .setChartType('ColumnChart')
  .setOptions({
    title:   'HLS scenes over the AOI by sensor (' + WINDOW_NAME + ')',
    hAxis:   {title: 'Year', format: '####'},
    vAxis:   {title: 'Scenes'},
    isStacked: true,
    series:  {0: {color: '#B07AA1'}, 1: {color: '#E07B39'}}
  }));


/* 
 *  PHASE 4 — MAPS: WHERE the gaps are
 *  */
//  A single AOI-wide percentage hides spatial structure: persistent cloud sits on
//  the Imatong and Agoro-Agu highlands, so a year can look fine overall while one
//  massif is unusable. Inspect two years side by side before committing.

var YEAR_A = 2019;   // <<< edit after reading the charts
var YEAR_B = 2024;   // <<< edit after reading the charts

function clearObsImage(year) {
  var win = windowFor(year);
  return clearObsFor(win.start, win.end);
}

var visObs = {min: 0, max: 20, palette: ['#67001f', '#f7f7f7', '#1b7837']};

Map.centerObject(aoi, 8);
Map.addLayer(clearObsImage(YEAR_A), visObs, 'Clear obs ' + YEAR_A + ' (' + WINDOW_NAME + ')');
Map.addLayer(clearObsImage(YEAR_B), visObs, 'Clear obs ' + YEAR_B + ' (' + WINDOW_NAME + ')');

// Pixels falling short of the threshold — the actual holes in a median composite.
Map.addLayer(clearObsImage(YEAR_A).lt(MIN_CLEAR_OBS).selfMask(),
             {palette: ['#d7301f']}, 'Under-sampled ' + YEAR_A, false);
Map.addLayer(clearObsImage(YEAR_B).lt(MIN_CLEAR_OBS).selfMask(),
             {palette: ['#d7301f']}, 'Under-sampled ' + YEAR_B, false);

// AOI outlines for context.
Map.addLayer(aoi_parts.style({color: '#000000', fillColor: '00000000', width: 1}),
             {}, 'PA polygons');
Map.addLayer(ee.Image().paint(aoi_union, 0, 2), {palette: ['#E07B39']}, 'PA union (buffered)');


/* 
 *  PHASE 5 — OPTIONAL: export the coverage table for the write-up
 *  */
//  Uncomment to get a CSV in Drive, then read it in R alongside the AOI audit.

// Export.table.toDrive({
//   collection:  coverage,
//   description: 'kidepo_hls_coverage_' + WINDOW_NAME + '_' + YEAR_START + '_' + YEAR_END,
//   folder:      'kidepo_lcluc',
//   fileFormat:  'CSV'
// });


/* 
 *  HOW TO READ THE RESULT — and what to decide
 * 
 *  1. From chart 3.1, shortlist years where pct_ok_all is high AND pct_ok_ssd is
 *     not far behind. A year carried only by the Uganda side is not a usable date.
 *  2. From chart 3.2, prefer a comfortable margin over MIN_CLEAR_OBS rather than a
 *     year that only just clears it.
 *  3. Separate the two chosen years by enough time for real change to accumulate —
 *     5+ years for a land-cover signal in a semi-arid rangeland.
 *  4. Re-run with WINDOW_NAME = 'dry_early' and compare. If dry_late wins on cloud
 *     but the burn-scar risk is judged too high, dry_early is the safer window even
 *     at some cost in clear observations. That trade-off is a decision to record,
 *     not one to leave implicit.
 *  5. Check the red "under-sampled" layers for BOTH candidate years. Persistent
 *     holes in the same place in both years are a reporting caveat; holes in only
 *     one year bias the change map.
 *
 *  Then set the two dates in scripts/gee/02_build_composites.js (next script).
 *  */

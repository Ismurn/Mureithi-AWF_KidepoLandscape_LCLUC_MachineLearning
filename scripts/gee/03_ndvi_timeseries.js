/*
 *  KIDEPO LANDSCAPE — LAND COVER / LAND USE CHANGE (LCLUC)
 *  GEE Script 03 — Annual dry_late NDVI, 2016-2026, with the two confounds
 *                  that would otherwise ruin it
 *
 *  WHY THIS EXISTS
 *  The central limitation of everything done so far is that it rests on TWO
 *  DATES. A +0.013 NDVI difference between 2019 and 2024 cannot be told apart
 *  from one drier year against one wetter one, and no amount of care with those
 *  two composites will fix it. An annual series across the same window can.
 *
 *  It asks two things the two-date design cannot:
 *
 *    1. Is the 2019->2024 greening part of a TREND, or are those two years
 *       simply placed either side of ordinary inter-annual swings?
 *    2. Do the candidate patches from script 09 move MONOTONICALLY, which is
 *       what thickening or land-use change looks like, or do they bounce, which
 *       is what rainfall looks like?
 *
 *  READ THIS BEFORE TRUSTING ANY TREND FROM THIS EXPORT
 *
 *  A naive NDVI trend over this record is confounded twice over, and both
 *  confounds are LARGER here than the signal being measured. The numbers below
 *  are from this project's own script 01 and 01b outputs, not from assumption:
 *
 *  CONFOUND 1 — OBSERVATION DEPTH RISES STEEPLY ACROSS THE RECORD.
 *  Mean clear looks per pixel in the dry_late window:
 *
 *      2016  4.05     2020  4.75     2024 10.41
 *      2017  3.15     2021  7.00     2025 13.85
 *      2018  6.10     2022 10.96     2026  7.87  (incomplete, see below)
 *      2019  8.13     2023 11.25
 *
 *  That is not weather, it is the satellite constellation growing — Sentinel-2B
 *  from 2017, Landsat 9 from 2021. Script 08 established that on THIS data more
 *  clear looks produce higher NDVI, and attributed 28% of the 2019->2024
 *  greening to a gap of only 2.3 looks. Across this record the gap is nearer
 *  TEN looks, so the depth artefact alone could manufacture a greening larger
 *  than the entire measured change. The annual clear-observation count is
 *  therefore exported as its own stack and MUST be carried as a covariate. A
 *  trend taken from the NDVI stack alone would be mostly an artefact of how
 *  many times each pixel was looked at.
 *
 *  CONFOUND 2 — FIRE VARIES ENORMOUSLY BETWEEN YEARS.
 *  Cumulative burnt fraction of the AOI by the window's end:
 *
 *      2016 23.4%    2020 31.1%    2024 22.8%
 *      2017 27.4%    2021 39.7%    2025 36.3%
 *      2018 43.9%    2022 36.8%    2026 38.3%
 *      2019 19.7%    2023 26.8%
 *
 *  2019 and 2024 are the two LEAST burnt years on record — which is exactly why
 *  script 02 chose them, and it means the chosen pair is systematically
 *  unrepresentative of the series. A raw annual series will show those two years
 *  as local maxima for a reason that has nothing to do with vegetation change.
 *  So the annual burnt mask is exported too, and must be carried as a covariate
 *  or used to drop burnt pixels year by year.
 *
 *  2016, 2017 AND 2026 ARE NOT FULLY USABLE.
 *  Script 01 found adequate coverage over only 65.7% of the AOI in 2016 and
 *  30.1% in 2017 — in 2017 the mean pixel got 3.15 clear looks, which is below
 *  the support floor this project uses. They are exported anyway, because
 *  dropping data at export time costs another GEE run and another download to
 *  reverse while dropping it in R is free, but the R stage builds its balanced
 *  panel over 2018-2026 and treats the first two years as context only. 2026 is
 *  the current year's window and its 7.87 mean looks against 2025's 13.85
 *  suggests HLS processing has not caught up; treat it as provisional.
 *
 *  WHAT IS EXPORTED — three rasters and one table.
 *
 *    1  kidepo_ndvi_annual_...   11 bands, NDVI x 10000, Int16
 *    2  kidepo_obs_annual_...    11 bands, clear-look COUNT, unscaled
 *    3  kidepo_burnt_annual_...  11 bands, 1 where burnt this fire season
 *    4  kidepo_ndvi_annual_aoi.csv   AOI-wide annual means, a cross-check
 *
 *  Band order in all three rasters is YEARS as listed below, ascending, and the
 *  R loader names them positionally from the same list rather than trusting the
 *  GeoTIFF's band descriptions.
 *
 *  METHOD MATCHES SCRIPT 02 EXACTLY, WHICH MATTERS MORE THAN IT LOOKS.
 *  Script 02 takes the MEDIAN OF THE REFLECTANCE BANDS and computes NDVI from
 *  that median — not the median of per-scene NDVI. Those two are NOT equal, and
 *  a series built the other way would not be comparable with the 2019 and 2024
 *  composites or with any number already reported. The compositing, the Fmask
 *  decoding, the cloud threshold and the window arithmetic are all reproduced
 *  from 02 verbatim for the same reason.
 *
 *  SIZE — this is the biggest download in the project, around 300-350 MB in
 *  total, nearly all of it the NDVI stack. The other two compress very well:
 *  counts are small integers and the burn product is 500 m replicated onto a
 *  30 m grid.
 *
 *  ENCODING
 *     dtype   Int16
 *     NDVI    value x 10000  (divide by 10000 in R)
 *     obs     raw count, no scaling
 *     burnt   0 or 1, no scaling
 *     nodata  -9999  (outside the AOI; and in NDVI where a pixel had no clear
 *             look at all that year)
 *     CRS     EPSG:32636, 30 m, same region as the composites
 *
 *  ASSETS USED (in users/imnduta1/KidepoAOI-LULC-ML/):
 *      kidepo_pa_union_shp
 *
 *  Author: Isaac Mureithi - IMureithi@awf.org
 */


//
// PHASE 0 — CONFIGURATION
//
// Every value here MUST match scripts/gee/02_build_composites.js.

var ASSET_ROOT = 'users/imnduta1/KidepoAOI-LULC-ML/';
var aoi_union  = ee.FeatureCollection(ASSET_ROOT + 'kidepo_pa_union_shp');

var WINDOW_NAME = 'dry_late';
var YEARS = [2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026];

var BUFFER_M   = 500;
var aoiExport  = aoi_union.geometry().buffer(BUFFER_M);

var MAX_SCENE_CLOUD = 60;
var EXPORT_CRS   = 'EPSG:32636';
var EXPORT_SCALE = 30;
var NODATA       = -9999;
var SCALE_INT    = 10000;

var WINDOW_DEFS = {
  dry_early: {startMonth: 12, startDay: 1, endMonth: 1, endDay: 31, endYearOffset: 1},
  dry_late:  {startMonth:  2, startDay: 1, endMonth: 3, endDay: 31, endYearOffset: 0}
};

var SEASON_START_MONTH = 11;    // fire season opens 1 Nov, as in 01b and 01c

var BANDS_COMMON = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2'];
var BANDS_L30    = ['B2', 'B3', 'B4', 'B5',  'B6',  'B7'];
var BANDS_S30    = ['B2', 'B3', 'B4', 'B8A', 'B11', 'B12'];

var BURN_COLL = 'MODIS/061/MCD64A1';

// NO SCALE FACTOR ON THE REFLECTANCE. HLS v002 in Earth Engine is already
// served in reflectance units; applying the documented 0.0001 again truncates
// every band to zero, and NDVI survives it looking perfect because a
// normalised difference is scale-invariant. See 02 for the full account.


//
// PHASE 1 — HELPERS, REPRODUCED FROM 02
//

/** 1.1 Leap-safe window dates. Identical construction to 01, 01b, 01c and 02. */
function windowFor(windowName, year) {
  var w     = WINDOW_DEFS[windowName];
  var start = ee.Date.fromYMD(year, w.startMonth, w.startDay);
  var end   = ee.Date
    .fromYMD(ee.Number(year).add(w.endYearOffset), w.endMonth, w.endDay)
    .advance(1, 'day');
  return {start: start, end: end};
}

/** 1.2 The 1 November immediately before the window opens. As 01c. */
function seasonStart(windowName, year) {
  var w      = WINDOW_DEFS[windowName];
  var offset = (w.startMonth >= SEASON_START_MONTH) ? 0 : -1;
  return ee.Date.fromYMD(year + offset, SEASON_START_MONTH, 1);
}

/**
 * 1.3 Fmask decoding — identical to 01 and 02.
 * Bits, LSB first: 0 cirrus | 1 cloud | 2 adjacent-to-cloud | 3 cloud shadow
 *                  4 snow/ice | 5 water | 6-7 aerosol
 * Water is NOT masked: it is a land-cover class to be classified, not noise.
 */
function clearBits(img) {
  var fmask = img.select('Fmask');
  return fmask.bitwiseAnd(1 << 0).eq(0)
    .and(fmask.bitwiseAnd(1 << 1).eq(0))
    .and(fmask.bitwiseAnd(1 << 2).eq(0))
    .and(fmask.bitwiseAnd(1 << 3).eq(0));
}

/** 1.4 One harmonised, masked, date-filtered product collection. */
function prepProduct(collectionId, bandsIn, start, end) {
  return ee.ImageCollection(collectionId)
    .filterBounds(aoiExport)
    .filterDate(start, end)
    .filterMetadata('CLOUD_COVERAGE', 'less_than', MAX_SCENE_CLOUD)
    .map(function (img) {
      var clear = clearBits(img);
      return img.select(bandsIn, BANDS_COMMON)
        .updateMask(clear)
        .copyProperties(img, ['system:time_start']);
    });
}

/** 1.5 Both products merged for one window. */
function hlsFor(windowName, year) {
  var win = windowFor(windowName, year);
  return prepProduct('NASA/HLS/HLSL30/v002', BANDS_L30, win.start, win.end)
    .merge(prepProduct('NASA/HLS/HLSS30/v002', BANDS_S30, win.start, win.end));
}


//
// PHASE 2 — ONE YEAR'S NDVI, OBSERVATION COUNT AND BURNT MASK
//

/**
 * 2.1 NDVI for one year, built the way 02 builds it: median the reflectance
 * first, then take the normalised difference of the medians. Doing it the
 * other way round — median of per-scene NDVI — gives a different number and
 * would silently break comparability with every figure already reported.
 */
function annualNdvi(year) {
  var col = hlsFor(WINDOW_NAME, year);
  var reflectance = col.select(BANDS_COMMON).median();
  return reflectance.normalizedDifference(['NIR', 'RED'])
    .multiply(SCALE_INT)
    .rename('NDVI_' + year);
}

/**
 * 2.2 Clear looks per pixel. The collection is already masked to clear by
 * prepProduct, so count() counts exactly the observations the median was taken
 * over. unmask(0) makes "never seen" an explicit zero rather than a hole.
 */
function annualObs(year) {
  var col = hlsFor(WINDOW_NAME, year);
  return col.select('RED').count().unmask(0).rename('obs_' + year);
}

/**
 * 2.3 Burnt at any point in the fire season up to the window's end.
 *
 * No ee.Algorithms.If and no ee.Image.constant anywhere in here, deliberately.
 * Script 01c lost two entire bands to that combination: a projectionless
 * constant sampled at 30 m evaluated to the constant and every fire vanished
 * while a sibling band from the same collection came out correct. Here
 * img.gt(0) inherits the MODIS projection, max() reduces within it, and
 * unmask(0) fills afterwards — so the projection is carried the whole way and
 * there is no branch to pick the wrong one.
 *
 * Phase 3 checks the result against the independently known burnt fractions
 * from 01b rather than assuming this reasoning is right.
 */
function annualBurnt(year) {
  var win    = windowFor(WINDOW_NAME, year);
  var season = seasonStart(WINDOW_NAME, year);

  return ee.ImageCollection(BURN_COLL)
    .filterBounds(aoiExport)
    .filterDate(season, win.end)
    .select('BurnDate')
    .map(function (img) { return img.gt(0); })
    .max()
    .unmask(0)
    .rename('burnt_' + year);
}

/** 2.4 Stack one band per year, in ascending year order. */
function stackYears(fn) {
  var imgs = YEARS.map(fn);
  return ee.Image(imgs.slice(1).reduce(
    function (acc, img) { return acc.addBands(img); }, imgs[0]));
}

var ndviStack  = stackYears(annualNdvi);
var obsStack   = stackYears(annualObs);
var burntStack = stackYears(annualBurnt);

/**
 * 2.5 Finish a stack for export. The valid-data mask comes from the OBSERVATION
 * stack, never from NDVI: a pixel with no clear look has no NDVI at all, and
 * using NDVI to define its own validity would be circular.
 */
function forExport(img, label) {
  return img
    .clip(aoiExport)       // clip BEFORE unmask — the lesson from 02
    .unmask(NODATA)
    .toInt16()
    .set({
      window:    WINDOW_NAME,
      years:     ee.List(YEARS).join(','),
      content:   label,
      scale_int: (label === 'ndvi') ? SCALE_INT : 1,
      nodata:    NODATA
    });
}

var ndviOut  = forExport(ndviStack,  'ndvi');
var obsOut   = forExport(obsStack,   'clear_obs');
var burntOut = forExport(burntStack, 'burnt_season');


//
// PHASE 3 — CHECK AGAINST WHAT IS ALREADY KNOWN, BEFORE SPENDING AN EXPORT
//
// These are not decorative. Script 01 and 01b already measured the mean clear
// looks and the burnt fraction for every one of these years by a different
// code path, so the two must agree. If they do not, something in the window
// arithmetic or the masking has drifted and the export is not worth making.
//
// EXPECTED mean clear looks, from kidepo_hls_coverage_bothwindows_2016_2026.csv:
//    2016 4.05  2017 3.15  2018 6.10  2019 8.13  2020 4.75  2021 7.00
//    2022 10.96 2023 11.25 2024 10.41 2025 13.85 2026 7.87
//
// EXPECTED burnt fraction, from kidepo_burn_by_window_2016_2026.csv:
//    2016 0.234  2017 0.274  2018 0.439  2019 0.197  2020 0.311  2021 0.397
//    2022 0.368  2023 0.268  2024 0.228  2025 0.363  2026 0.383
//
// A burnt fraction of ZERO for every year is the 01c failure repeating itself.
// An NDVI mean near zero means the reflectance scale factor is being applied
// twice. An NDVI mean in the tens of thousands means it is not being applied
// at all — stored units are NDVI x 10000, so expect roughly 1000-5000.

// No mask is applied in these checks, and that is deliberate rather than an
// oversight. reduceRegion's `geometry` already confines the statistic to the
// AOI, and the obvious-looking obsStack.gte(0) would NOT have confined
// anything: annualObs() ends in unmask(0), which clears the mask entirely and
// makes the band globally valid, so gte(0) is true over the whole planet. NDVI
// carries its own mask where a pixel had no clear look, so masked pixels drop
// out of its mean on their own.

function aoiMeans(img, label) {
  return img.reduceRegion({
    reducer:    ee.Reducer.mean(),
    geometry:   aoiExport,
    scale:      240,          // coarse on purpose: a smoke test, not a statistic
    maxPixels:  1e10,
    bestEffort: true
  }).set('content', label);
}

print('SANITY — mean clear looks per year (compare with the table above):',
      aoiMeans(obsStack, 'clear_obs'));
print('SANITY — burnt fraction per year (ZERO everywhere = the 01c bug again):',
      aoiMeans(burntStack, 'burnt_season'));
print('SANITY — mean NDVI per year, stored units, expect ~1000-5000:',
      aoiMeans(ndviStack, 'ndvi_x10000'));

print('Band names, NDVI stack (must be ascending years):', ndviStack.bandNames());

// Scene counts per year, as a second independent handle on coverage.
var sceneCounts = ee.FeatureCollection(YEARS.map(function (y) {
  return ee.Feature(null, {year: y, n_scenes: hlsFor(WINDOW_NAME, y).size()});
}));
print('Scenes contributing per year:', sceneCounts);

var ndviPalette = ['#a6611a', '#dfc27d', '#f5f5f5', '#80cdc1', '#018571'];

Map.centerObject(aoiExport, 8);
Map.addLayer(ndviStack.select('NDVI_2019').clip(aoiExport),
             {min: 1000, max: 5000, palette: ndviPalette}, 'NDVI 2019');
Map.addLayer(ndviStack.select('NDVI_2024').clip(aoiExport),
             {min: 1000, max: 5000, palette: ndviPalette}, 'NDVI 2024', false);
Map.addLayer(obsStack.select('obs_2017').clip(aoiExport),
             {min: 0, max: 12, palette: ['#000000', '#ffffff']},
             'clear looks 2017 (the thinnest year)', false);
Map.addLayer(burntStack.select('burnt_2018').clip(aoiExport).selfMask(),
             {palette: ['#d7301f']}, 'burnt 2018 (the worst fire year)', false);
Map.addLayer(ee.Image().paint(aoi_union, 0, 2), {palette: ['#1b7837']}, 'PA union');


//
// PHASE 4 — EXPORT
//

var span = YEARS[0] + '_' + YEARS[YEARS.length - 1];

function exportStack(img, kind) {
  var name = 'kidepo_' + kind + '_annual_' + WINDOW_NAME + '_' + span +
             '_utm36n_30m';
  Export.image.toDrive({
    image:          img,
    description:    name,
    fileNamePrefix: name,
    folder:         'kidepo_lcluc',
    region:         aoiExport,
    crs:            EXPORT_CRS,
    scale:          EXPORT_SCALE,
    maxPixels:      1e13,
    fileFormat:     'GeoTIFF',
    formatOptions:  {cloudOptimized: true}
  });
}

exportStack(ndviOut,  'ndvi');
exportStack(obsOut,   'obs');
exportStack(burntOut, 'burnt');

// The AOI-wide annual table. Independent of the rasters and far cheaper, so it
// is worth having even though R will recompute these numbers properly from the
// stacks: if the two disagree, something is wrong with the download rather than
// with the analysis.
var annualTable = ee.FeatureCollection(YEARS.map(function (y) {
  var ndvi  = annualNdvi(y);
  var obs   = annualObs(y);
  var burnt = annualBurnt(y);
  var stats = ndvi.addBands(obs).addBands(burnt)
    .reduceRegion({
      reducer:    ee.Reducer.mean(),
      geometry:   aoiExport,
      scale:      240,
      maxPixels:  1e10,
      bestEffort: true
    });
  return ee.Feature(null, {
    year:          y,
    n_scenes:      hlsFor(WINDOW_NAME, y).size(),
    mean_ndvi:     ee.Number(stats.get('NDVI_' + y)).divide(SCALE_INT),
    mean_obs:      stats.get('obs_' + y),
    burnt_frac:    stats.get('burnt_' + y)
  });
}));

Export.table.toDrive({
  collection:     annualTable,
  description:    'kidepo_ndvi_annual_aoi',
  fileNamePrefix: 'kidepo_ndvi_annual_aoi',
  folder:         'kidepo_lcluc',
  fileFormat:     'CSV',
  selectors:      ['year', 'n_scenes', 'mean_ndvi', 'mean_obs', 'burnt_frac']
});


/*
 *  AFTER THE EXPORTS FINISH
 *
 *  Download all four into the same folder as the composites, then run
 *      scripts/R/11_ndvi_timeseries.R
 *
 *  which does the following, in this order, because the order is the argument:
 *
 *   1. Rebuild the annual series and check it against the AOI table and against
 *      the 01/01b figures printed above. Three independent code paths agreeing
 *      is the only reason to believe any of them.
 *
 *   2. ADJUST FOR THE TWO CONFOUNDS BEFORE LOOKING AT THE TREND. Regress NDVI
 *      on year WITH clear looks and burnt status held constant. Report the raw
 *      and adjusted trends side by side. The gap between them is the headline
 *      result of this script, whichever way it falls — if the raw trend is
 *      mostly depth, then so was part of the two-date estimate.
 *
 *   3. Place 2019 and 2024 in the distribution of years. If the 2019->2024 step
 *      is ordinary against the year-to-year spread, the +0.013 is a sample of
 *      two from a noisy series and should be reported with that spread attached
 *      rather than as a five-year change.
 *
 *   4. Compare candidate trajectories with background ones. Monotone movement
 *      is thickening or land-use change; bouncing is rainfall. This is the one
 *      test that can distinguish the two explanations script 09 left open, and
 *      it is the reason this export is worth 350 MB.
 *
 *   5. A depth-adjusted per-pixel trend map, for the report.
 */

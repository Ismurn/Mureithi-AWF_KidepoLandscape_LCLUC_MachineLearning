/*
 *  KIDEPO LANDSCAPE — LAND COVER / LAND USE CHANGE (LCLUC)
 *  GEE Script 01c — Export MODIS burned-area rasters for the two comparison dates
 *
 *  WHY THIS EXISTS
 *  Script 04 tried to separate fire from land-cover change without a burnt-area
 *  raster, by taking the NBR threshold that flags the fraction of the AOI MODIS
 *  says burnt. It did not work: burn switching predicted change only 1.24x
 *  better than chance, and the residual still drifted out of the most burn-like
 *  cluster. The likely reason is scar AGE. By February a scar laid down in
 *  December has greened partway back and is no longer in the lowest NBR
 *  quantile, while rock, shadow and water margins sit there permanently — so
 *  the threshold flagged the wrong pixels.
 *
 *  The fix is to stop inferring where the fires were and export where MODIS
 *  actually recorded them.
 *
 *  WHAT IS EXPORTED, per comparison date — 3 bands:
 *     1 burnt_season      1 if the pixel burnt between the fire season opening
 *                         (1 Nov preceding the window) and the window's end
 *     2 burnt_in_window   1 if it burnt DURING the window itself, i.e. while
 *                         the composite was being collected — these are the
 *                         pixels whose own median straddles their fire
 *     3 days_since_burn   days between the MOST RECENT burn and the window's
 *                         end; -1 where the pixel did not burn in the season
 *
 *  Band 3 is the one script 04 lacked. Scar visibility decays with age, so a
 *  7-day-old scar and a 120-day-old scar are different objects spectrally even
 *  though a binary mask calls them both "burnt". With this band the R stage can
 *  weight, stratify or threshold by age instead of guessing.
 *
 *  RESOLUTION — READ THIS BEFORE USING IT AT FACE VALUE
 *  MCD64A1 is a 500 m product. Exporting on the 30 m composite grid is for
 *  stacking convenience ONLY; it does not create 30 m information. Each 500 m
 *  cell is replicated across roughly 278 thirty-metre pixels, so scar edges are
 *  accurate to about half a kilometre and small burns may be missed entirely.
 *  Use it to characterise and stratify, never to delineate.
 *
 *  ENCODING — must match what the R stage expects:
 *     dtype   Int16
 *     nodata  -9999  (outside the AOI)
 *     -1      in days_since_burn means "did not burn this season" (not nodata)
 *     CRS     EPSG:32636, 30 m, same region as the composites, so the rasters
 *             stack directly onto them with no resampling
 *
 *  ASSETS USED (in users/imnduta1/KidepoAOI-LULC-ML/):
 *      kidepo_pa_union_shp
 *
 *  Author: Isaac Mureithi - IMureithi@awf.org
 */


// 
// PHASE 0 — CONFIGURATION
// 
// These MUST match scripts/gee/02_build_composites.js, or the exported rasters
// will not align with the composites they are meant to explain.

var ASSET_ROOT = 'users/imnduta1/KidepoAOI-LULC-ML/';
var aoi_union  = ee.FeatureCollection(ASSET_ROOT + 'kidepo_pa_union_shp');

var WINDOW_NAME = 'dry_late';
var YEARS       = [2019, 2024];

var BUFFER_M   = 500;
var aoiExport  = aoi_union.geometry().buffer(BUFFER_M);

var EXPORT_CRS   = 'EPSG:32636';
var EXPORT_SCALE = 30;
var NODATA       = -9999;
var NO_BURN      = -1;      // days_since_burn where the pixel did not burn

// Seasonal windows — same definitions as 01, 01b and 02.
var WINDOW_DEFS = {
  dry_early: {startMonth: 12, startDay: 1, endMonth: 1, endDay: 31, endYearOffset: 1},
  dry_late:  {startMonth:  2, startDay: 1, endMonth: 3, endDay: 31, endYearOffset: 0}
};

// Fire season opens on 1 November, before either window. Same assumption as 01b.
var SEASON_START_MONTH = 11;

var BURN_COLL = 'MODIS/061/MCD64A1';
var MS_PER_DAY = 86400000;


// 
// PHASE 1 — HELPERS
// 

var burn = ee.ImageCollection(BURN_COLL)
  .filterBounds(aoiExport)
  .select('BurnDate');

/** 1.1 Leap-safe window dates, identical in construction to the other scripts. */
function windowFor(windowName, year) {
  var w     = WINDOW_DEFS[windowName];
  var start = ee.Date.fromYMD(year, w.startMonth, w.startDay);
  var end   = ee.Date
    .fromYMD(ee.Number(year).add(w.endYearOffset), w.endMonth, w.endDay)
    .advance(1, 'day');
  return {start: start, end: end};
}

/** 1.2 The 1 November immediately before the window opens. */
function seasonStart(windowName, year) {
  var w      = WINDOW_DEFS[windowName];
  var offset = (w.startMonth >= SEASON_START_MONTH) ? 0 : -1;
  return ee.Date.fromYMD(year + offset, SEASON_START_MONTH, 1);
}

/*
 * 1.3 There is deliberately NO separate "burnt at any point" function.
 *
 * The first version of this script computed the binary bands independently, as
 * col.map(i => i.gt(0).unmask(0)).max() wrapped in an ee.Algorithms.If guard.
 * Both binary bands exported as ZERO everywhere while days_since_burn came back
 * correct from the same collection and the same dates — so the fires were
 * found, and something between the reducer and the export discarded them. The
 * most likely culprit is the If: its two branches had different projections
 * (MODIS 500 m against a projectionless ee.Image.constant), and an image with
 * no fixed projection sampled at 30 m evaluates to the constant.
 *
 * Rather than diagnose that further, both binary bands are now DERIVED from
 * days_since_burn below. One computation path, no duplicated logic, and the
 * bands cannot disagree with the age band because they are made from it.
 *
 * The emptiness guards have gone with it. They were protecting against a
 * collection with no products, which cannot occur here: MCD64A1 runs from 2000
 * and every window this project uses falls inside that. A guard that is never
 * exercised is a guard that can only misbehave.
 */

/**
 * 1.4 Days between the most recent burn in the period and `end`.
 *
 * MCD64A1 stores BurnDate as a day-of-year WITHIN each monthly product, so the
 * DOY has to be anchored to that product's own year before it means anything —
 * a January product and a December product both carry small and large DOYs
 * respectively, and mixing them up would scramble the ordering.
 *
 * Taking the MINIMUM days-since across the period keeps the MOST RECENT fire,
 * which is the one that governs how the ground looks now. An older scar
 * underneath is irrelevant once it has burnt again.
 */
function daysSinceBurn(start, end) {

  var col = burn.filterDate(start, end);
  var endMillis = ee.Number(end.millis());

  var stamped = col.map(function (img) {
    var doy  = img.select('BurnDate');
    var hit  = doy.gt(0);
    var year = ee.Date(img.get('system:time_start')).get('year');

    // Absolute burn time = 1 Jan of the product's year + (DOY - 1) days.
    var janFirst   = ee.Number(ee.Date.fromYMD(year, 1, 1).millis());
    var burnMillis = ee.Image.constant(janFirst)
      .add(doy.subtract(1).multiply(MS_PER_DAY));

    return ee.Image.constant(endMillis)
      .subtract(burnMillis)
      .divide(MS_PER_DAY)
      .updateMask(hit)
      .rename('days_since_burn');
  });

  // min() keeps the most recent fire. Every pixel that never burnt stays masked
  // through the reduction and is then filled with the NO_BURN sentinel, so the
  // band is complete and carries the MODIS projection throughout.
  return stamped.min()
    .unmask(NO_BURN)
    .rename('days_since_burn');
}


// 
// PHASE 2 — BUILD ONE BURN RASTER PER COMPARISON DATE
// 

function burnRasterFor(year) {

  var win    = windowFor(WINDOW_NAME, year);
  var season = seasonStart(WINDOW_NAME, year);

  var daysSince = daysSinceBurn(season, win.end);

  // Both binary bands are derived from the age band — see note 1.3.
  // days_since_burn is measured back from win.end, so a scar that burnt ON the
  // window's opening day has an age equal to the window length, and anything
  // younger burnt inside the window. NO_BURN is negative, so gte(0) separates
  // burnt from unburnt without a second collection pass.
  var winDays = ee.Number(win.end.difference(win.start, 'day'));

  var burntSeason = daysSince.gte(0).rename('burnt_season');
  var burntWindow = daysSince.gte(0)
    .and(daysSince.lte(ee.Image.constant(winDays)))
    .rename('burnt_in_window');

  return burntSeason
    .addBands(burntWindow)
    .addBands(daysSince)
    .clip(aoiExport)        // clip BEFORE unmask, so the sentinel reaches
    .unmask(NODATA)         // every pixel outside the AOI (lesson from 02)
    .toInt16()
    .set({
      window:       WINDOW_NAME,
      year:         year,
      season_start: season.format('YYYY-MM-dd'),
      window_end:   win.end.advance(-1, 'day').format('YYYY-MM-dd')
    });
}

var burn2019 = burnRasterFor(YEARS[0]);
var burn2024 = burnRasterFor(YEARS[1]);


// 
// PHASE 3 — SANITY CHECK BEFORE EXPORTING
// 
// Expected, from the 01b burn profile: burnt_season should average close to
// 0.197 for 2019 and 0.228 for 2024. A mean far from those means the season
// window or the date arithmetic is wrong — check that before exporting.
//
// days_since_burn reports a MINIMUM of -1, not 0: unburnt pixels carry the
// NO_BURN sentinel and they are the majority, so its mean is meaningless here.
// What matters is the MAXIMUM, which should be at or just under 151 — the days
// from 1 Nov to 31 Mar. A maximum well above that means the DOY is being
// anchored to the wrong year somewhere in the date arithmetic.

function checkStats(img, year) {
  var valid = img.select('burnt_season').gte(0);   // -9999 outside the AOI
  return img.updateMask(valid).reduceRegion({
    reducer:    ee.Reducer.mean().combine({reducer2: ee.Reducer.minMax(), sharedInputs: true}),
    geometry:   aoiExport,
    scale:      480,          // near MODIS native; this is a check, not a statistic
    maxPixels:  1e10,
    bestEffort: true
  }).set('year', year);
}

print('SANITY CHECK ' + YEARS[0] + ':', checkStats(burn2019, YEARS[0]));
print('SANITY CHECK ' + YEARS[1] + ':', checkStats(burn2024, YEARS[1]));

Map.centerObject(aoiExport, 8);
Map.addLayer(burn2019.select('burnt_season').selfMask(),
             {palette: ['#d7301f']}, YEARS[0] + ' burnt this season');
Map.addLayer(burn2024.select('burnt_season').selfMask(),
             {palette: ['#6a3d9a']}, YEARS[1] + ' burnt this season', false);
Map.addLayer(burn2019.select('days_since_burn').updateMask(
               burn2019.select('burnt_season').eq(1)),
             {min: 0, max: 151, palette: ['#000000', '#d7301f', '#fdae61', '#f7f7f7']},
             YEARS[0] + ' days since burn (dark = recent)', false);
Map.addLayer(ee.Image().paint(aoi_union, 0, 2), {palette: ['#1b7837']}, 'PA union');


// 
// PHASE 4 — EXPORT
// 
// Same CRS, scale and region as the composites, so the two stack without
// resampling. About 9.4 million pixels, 3 Int16 bands — a small file, since a
// 500 m product replicated onto a 30 m grid compresses very well.

function exportBurn(img, year) {
  var name = 'kidepo_burn_' + WINDOW_NAME + '_' + year + '_utm36n_30m';
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

exportBurn(burn2019, YEARS[0]);
exportBurn(burn2024, YEARS[1]);


/*
 *  AFTER THE EXPORTS FINISH
 *
 *  Download both into the same folder as the composites, then the R stage can
 *  redo the attribution against recorded fires rather than an NBR quantile.
 *
 *  The questions worth putting to the data once it is there:
 *
 *   1. Does the calibrated proxy in 04 agree with MODIS at all? Cross-tabulate
 *      fire_state.tif against burnt_season. If agreement is poor, that confirms
 *      the proxy was flagging the wrong pixels and the 851 km2 residual was
 *      never meaningful.
 *
 *   2. Does scar AGE explain the failure? Compare NBR against days_since_burn.
 *      If NBR recovers towards unburnt values beyond roughly 60-90 days, that
 *      is the mechanism, and it also tells us which scars a February composite
 *      can still see at all.
 *
 *   3. Redo the attribution with the real mask, and re-run the association
 *      check in 04 Phase 4.5. If enrichment stays near 1.0 even against
 *      recorded fires, then fire genuinely is not driving the change and the
 *      residual deserves to be taken seriously as land cover.
 *
 *   4. Only then: the supervised classification, with burn handled explicitly —
 *      either as its own class or by excluding recently burnt pixels from
 *      training.
 */
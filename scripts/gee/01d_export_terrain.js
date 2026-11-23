/*
 *  KIDEPO LANDSCAPE — LAND COVER / LAND USE CHANGE (LCLUC)
 *  GEE Script 01d — Export terrain and observation timing, to test the
 *                   terrain-illumination explanation for the candidate change
 *
 *  WHY THIS EXISTS
 *  Script 09 examined the 277.5 km2 of candidate change and found that it
 *  concentrates in the small HILL reserves — Timu 26.5% of its area, Rom 23.5%,
 *  Imatong 7.2% — while Kidepo Valley, the large lowland park, has 0.5%. Three
 *  explanations survived that check, and one of them is an artefact:
 *
 *    (a) genuine woody thickening or regrowth
 *    (b) DIFFERENTIAL TERRAIN ILLUMINATION — the candidates favour rough,
 *        hilly ground, and reflectance on a slope depends on how the sun
 *        strikes it, so two composites could differ on slopes without the
 *        ground having changed at all
 *    (c) real land-use change
 *
 *  (b) is the one that can be settled with data already in orbit, and it must
 *  be settled first: if it holds, part of the 277.5 km2 is not change at all.
 *
 *  WHAT THE ILLUMINATION HYPOTHESIS ACTUALLY PREDICTS — read this, because the
 *  test has a predicted NULL and that is not an excuse made after the fact.
 *
 *  Both composites use the SAME calendar window, 1 Feb to 31 Mar. The sun is
 *  therefore in nearly the same place in both, so the local illumination
 *  geometry cos(i) is nearly the same in both, so it very largely CANCELS in
 *  the difference. Illumination sets the LEVEL of reflectance on a slope; it
 *  cannot by itself produce a CHANGE between two same-season composites. For it
 *  to produce one, something about the observations must differ between years,
 *  and there are two routes:
 *
 *    ROUTE 1  TIMING. If the clear looks in one year fall earlier or later in
 *             the window than in the other, the mean solar declination differs,
 *             so cos(i) differs, and the difference is largest on steep slopes.
 *             Band obs_doy_<year> below measures this per pixel.
 *
 *    ROUTE 2  SHADOW MASKING. Fmask flags cloud shadow. Terrain shadow on a
 *             steep slope looks much the same, so shaded slopes lose
 *             observations, and which observations survive depends on how many
 *             there were. Script 08 already found 28% of the greening is
 *             attributable to observation depth; if that artefact is
 *             concentrated on steep ground then terrain is its mechanism.
 *             Testing this needs slope crossed with the clear_obs bands the
 *             composites already carry.
 *
 *  So this export carries the two things the R stage lacks: TERRAIN, and WHEN
 *  within the window each pixel was actually seen.
 *
 *  WHAT IS EXPORTED — one image, 5 bands, in this order:
 *     1 elevation      metres above sea level, NO scaling
 *     2 slope          degrees x 10      (0-900)
 *     3 aspect         degrees x 10      (0-3600, clockwise from north)
 *     4 obs_doy_2019   mean day-of-year of the CLEAR looks that built the 2019
 *                      composite, NO scaling (expect ~32-91)
 *     5 obs_doy_2024   the same for 2024
 *
 *  ENCODING — the scale factors are NOT uniform and that is a hazard. This
 *  project has already lost one export to a scale factor applied twice (see
 *  02_build_composites.js). The R loader therefore does not trust these
 *  comments: it applies a per-band scale from one lookup table and then ASSERTS
 *  that each band lands in its physically possible range — slope 0-90, aspect
 *  0-360, elevation 500-4000 m for this landscape, DOY inside the window. A
 *  wrong scale cannot pass all four.
 *
 *     dtype   Int16
 *     nodata  -9999  (outside the AOI, and where a pixel had no clear look)
 *     CRS     EPSG:32636, 30 m, same region as the composites, so the rasters
 *             stack directly onto them with no resampling
 *
 *  ASPECT IS MEANINGLESS ON FLAT GROUND. Where slope is 0 the aspect band
 *  carries 0, which reads as "due north" and is not. The R stage excludes
 *  pixels below 2 degrees of slope from every aspect statistic.
 *
 *  ALSO EXPORTED — one small CSV, kidepo_hls_sun_geometry.csv: the per-scene
 *  solar and view angles from the HLS metadata, one row per contributing
 *  scene. The R stage uses it to compute each year's mean illumination geometry
 *  from the real recorded angles rather than from an astronomical estimate. If
 *  Earth Engine's HLS ingestion does not carry these property names the columns
 *  arrive EMPTY — the CSV still exports and the R stage falls back to computing
 *  the angles from date and latitude. Phase 5 prints the property names of a
 *  real scene so the truth is visible either way.
 *
 *  DEM CHOICE: NASADEM rather than raw SRTM. It is the same 1-arc-second
 *  mission reprocessed with the voids filled, and voids in SRTM sit exactly
 *  where this test is most sensitive — steep slopes and deep shadow. A void
 *  would silently drop the steepest ground out of the regression.
 *
 *  ASSETS USED (in users/imnduta1/KidepoAOI-LULC-ML/):
 *      kidepo_pa_union_shp
 *
 *  Author: Isaac Mureithi - IMureithi@awf.org
 */


//
// PHASE 0 — CONFIGURATION
//
// Every value here MUST match scripts/gee/02_build_composites.js. The whole
// point of this export is that it stacks onto the composites.

var ASSET_ROOT = 'users/imnduta1/KidepoAOI-LULC-ML/';
var aoi_union  = ee.FeatureCollection(ASSET_ROOT + 'kidepo_pa_union_shp');

var WINDOW_NAME = 'dry_late';
var YEARS       = [2019, 2024];

var BUFFER_M   = 500;
var aoiExport  = aoi_union.geometry().buffer(BUFFER_M);

var EXPORT_CRS   = 'EPSG:32636';
var EXPORT_SCALE = 30;
var NODATA       = -9999;

var MAX_SCENE_CLOUD = 60;    // same threshold as 01 and 02
var ANGLE_SCALE     = 10;    // slope and aspect are stored x10

var WINDOW_DEFS = {
  dry_early: {startMonth: 12, startDay: 1, endMonth: 1, endDay: 31, endYearOffset: 1},
  dry_late:  {startMonth:  2, startDay: 1, endMonth: 3, endDay: 31, endYearOffset: 0}
};

var BANDS_COMMON = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2'];
var BANDS_L30    = ['B2', 'B3', 'B4', 'B5',  'B6',  'B7'];
var BANDS_S30    = ['B2', 'B3', 'B4', 'B8A', 'B11', 'B12'];

var DEM_ID = 'NASA/NASADEM_HGT/001';

// The metadata fields wanted from each HLS scene. Names per the HLS v002
// product specification; see the caveat in the header.
var SUN_FIELDS = ['MEAN_SUN_ZENITH_ANGLE', 'MEAN_SUN_AZIMUTH_ANGLE',
                  'MEAN_VIEW_ZENITH_ANGLE', 'MEAN_VIEW_AZIMUTH_ANGLE'];


//
// PHASE 1 — HELPERS, COPIED FROM 02 SO THE OBSERVATIONS MATCH
//
// obs_doy must describe the SAME observations that formed the composite. Any
// difference in the date filter, the cloud threshold or the Fmask decoding
// would make it describe a different set of looks, and the regression would be
// controlling for the wrong thing. So these functions are reproduced verbatim
// rather than approximated.

/** 1.1 Leap-safe window dates. Identical construction to 01, 01b, 01c and 02. */
function windowFor(windowName, year) {
  var w     = WINDOW_DEFS[windowName];
  var start = ee.Date.fromYMD(year, w.startMonth, w.startDay);
  var end   = ee.Date
    .fromYMD(ee.Number(year).add(w.endYearOffset), w.endMonth, w.endDay)
    .advance(1, 'day');
  return {start: start, end: end};
}

/**
 * 1.2 Fmask decoding — identical to 01 and 02.
 * Bits, LSB first: 0 cirrus | 1 cloud | 2 adjacent-to-cloud | 3 cloud shadow
 *                  4 snow/ice | 5 water | 6-7 aerosol
 * Water is NOT masked, for the reasons given in 02.
 */
function clearBits(img) {
  var fmask = img.select('Fmask');
  return fmask.bitwiseAnd(1 << 0).eq(0)
    .and(fmask.bitwiseAnd(1 << 1).eq(0))
    .and(fmask.bitwiseAnd(1 << 2).eq(0))
    .and(fmask.bitwiseAnd(1 << 3).eq(0));
}

/** 1.3 One harmonised, masked, date-filtered product collection. */
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

/** 1.4 Both products merged for one window. */
function hlsFor(windowName, year) {
  var win = windowFor(windowName, year);
  return prepProduct('NASA/HLS/HLSL30/v002', BANDS_L30, win.start, win.end)
    .merge(prepProduct('NASA/HLS/HLSS30/v002', BANDS_S30, win.start, win.end));
}

/** 1.5 The same scenes WITHOUT the band rename, so the metadata survives. */
function hlsRawFor(windowName, year) {
  var win = windowFor(windowName, year);
  var one = function (id) {
    return ee.ImageCollection(id)
      .filterBounds(aoiExport)
      .filterDate(win.start, win.end)
      .filterMetadata('CLOUD_COVERAGE', 'less_than', MAX_SCENE_CLOUD);
  };
  return one('NASA/HLS/HLSL30/v002').merge(one('NASA/HLS/HLSS30/v002'));
}


//
// PHASE 2 — MEAN OBSERVATION DAY-OF-YEAR, PER PIXEL, PER YEAR
//
// This is Route 1 of the illumination hypothesis made measurable. A pixel whose
// clear looks all fall in early February was seen under a different sun from
// one seen in late March, and on a steep slope that matters.
//
// The DOY image is built by multiplying a real HLS band by zero and adding the
// scalar, NOT by ee.Image.constant. The difference is not cosmetic: a constant
// has no fixed projection, and 01c lost two whole bands to exactly that — a
// projectionless image sampled at 30 m evaluated to the constant and the fires
// vanished. Deriving it from the band inherits both the projection and the
// clear-look mask, so the mean is taken over precisely the observations that
// built the composite.
//
// dry_late runs 1 Feb to 31 Mar with endYearOffset 0, so the window never
// crosses a new year and plain day-of-year needs no unwrapping. If the window
// is ever changed to dry_early, which does cross, this arithmetic breaks.

function meanObsDoy(windowName, year) {

  var col = hlsFor(windowName, year);

  var stamped = col.map(function (img) {
    var doy = ee.Date(img.get('system:time_start')).getRelative('day', 'year').add(1);
    return img.select('RED')
      .multiply(0)
      .add(ee.Number(doy))
      .rename('obs_doy')
      .toFloat();
  });

  return stamped.mean().rename('obs_doy_' + year);
}

var doy2019 = meanObsDoy(WINDOW_NAME, YEARS[0]);
var doy2024 = meanObsDoy(WINDOW_NAME, YEARS[1]);


//
// PHASE 3 — TERRAIN
//
// Slope and aspect are computed from the DEM as served and exported at 30 m in
// UTM. NASADEM is 1 arc-second, about 30.9 m at this latitude, so the export
// resamples by a factor near 1 and the derivatives are effectively native. No
// .reproject() is used: it would force the whole DEM through a resampling at
// the requested scale for no gain in information and a large cost in memory.

var dem  = ee.Image(DEM_ID).select('elevation');
var terr = ee.Terrain.products(dem);

var elevation = terr.select('elevation').rename('elevation');
var slope     = terr.select('slope').multiply(ANGLE_SCALE).rename('slope');
var aspect    = terr.select('aspect').multiply(ANGLE_SCALE).rename('aspect');

var stack = elevation
  .addBands(slope)
  .addBands(aspect)
  .addBands(doy2019)
  .addBands(doy2024)
  .clip(aoiExport)        // clip BEFORE unmask, so the sentinel reaches every
  .unmask(NODATA)         // pixel outside the AOI — the lesson from 02
  .toInt16()
  .set({
    window:      WINDOW_NAME,
    year_a:      YEARS[0],
    year_b:      YEARS[1],
    dem:         DEM_ID,
    angle_scale: ANGLE_SCALE
  });


//
// PHASE 4 — SOLAR AND VIEW GEOMETRY, ONE ROW PER SCENE
//
// Exported as a table rather than printed, so the R stage reads the recorded
// angles instead of an astronomical approximation of them. One row per scene
// that passed the cloud filter, both years together.
//
// If a property name is wrong for this ingestion the column arrives empty and
// nothing else fails — the raster export is independent and the R stage has a
// fallback. The print in Phase 5 shows what the properties are really called.

function sunRows(year) {
  return hlsRawFor(WINDOW_NAME, year).map(function (img) {
    var d = ee.Date(img.get('system:time_start'));
    return ee.Feature(null, ee.Dictionary({
      scene_id:  img.get('system:index'),
      year:      year,
      date:      d.format('YYYY-MM-dd'),
      doy:       d.getRelative('day', 'year').add(1),
      cloud_pct: img.get('CLOUD_COVERAGE')
    }).combine(img.toDictionary(SUN_FIELDS), false));
  });
}

var sunFc = ee.FeatureCollection(sunRows(YEARS[0]))
  .merge(ee.FeatureCollection(sunRows(YEARS[1])));


//
// PHASE 5 — CHECKS BEFORE EXPORTING
//
// What to expect, and what a wrong answer means:
//
//   elevation   roughly 500-4000 m. Imatong carries the highest ground in the
//               landscape, Kidepo Valley the lowest. A minimum at or below 0
//               means NASADEM voids or a clipping failure.
//   slope       mean a few tens (degrees x10), maximum well under 900. A
//               maximum AT 900 means the band saturated.
//   obs_doy     these are UNSCALED, so expect 32-91. Anything outside that
//               means the window definition and the DOY arithmetic disagree.
//
// The difference between the two years' mean obs_doy is the headline number
// here. If it is small — a few days — then Route 1 of the illumination
// hypothesis is dead on arrival, and the R stage should say so plainly rather
// than hunting for an effect that has no mechanism behind it.

print('HLS scene property names (check the sun-angle fields exist):',
      hlsRawFor(WINDOW_NAME, YEARS[0]).first().propertyNames());

print('Scenes contributing, ' + YEARS[0] + ' then ' + YEARS[1] + ':',
      hlsRawFor(WINDOW_NAME, YEARS[0]).size(),
      hlsRawFor(WINDOW_NAME, YEARS[1]).size());

print('Sun geometry rows (first 10) — empty angle columns mean the property ' +
      'names differ in this ingestion:', sunFc.limit(10));

var valid = stack.select('elevation').gt(NODATA);

print('RASTER SANITY CHECK (elevation m; slope and aspect x10; obs_doy raw):',
      stack.updateMask(valid).reduceRegion({
        reducer: ee.Reducer.mean().combine({
          reducer2: ee.Reducer.minMax(), sharedInputs: true}),
        geometry:   aoiExport,
        scale:      240,        // a check, not a statistic
        maxPixels:  1e10,
        bestEffort: true
      }));

Map.centerObject(aoiExport, 8);
Map.addLayer(stack.select('elevation').updateMask(valid),
             {min: 600, max: 2800,
              palette: ['#2b83ba', '#abdda4', '#ffffbf', '#fdae61', '#d7191c']},
             'elevation');
Map.addLayer(stack.select('slope').updateMask(valid),
             {min: 0, max: 250, palette: ['#ffffff', '#000000']},
             'slope x10', false);
Map.addLayer(ee.Terrain.hillshade(dem.multiply(3)), {min: 0, max: 255},
             'hillshade (x3 relief, for looking only)', false);
Map.addLayer(stack.select('obs_doy_' + YEARS[1]).updateMask(valid)
               .subtract(stack.select('obs_doy_' + YEARS[0]).updateMask(valid)),
             {min: -20, max: 20, palette: ['#4575b4', '#ffffbf', '#d73027']},
             'obs DOY shift, ' + YEARS[1] + ' minus ' + YEARS[0], false);
Map.addLayer(ee.Image().paint(aoi_union, 0, 2), {palette: ['#1b7837']}, 'PA union');


//
// PHASE 6 — EXPORT
//
// About 32 million pixels, 5 Int16 bands. Terrain does not compress as well as
// the 500 m burn product did — expect a file in the tens of megabytes rather
// than the single digits.

var rasterName = 'kidepo_terrain_' + WINDOW_NAME + '_utm36n_30m';

Export.image.toDrive({
  image:          stack,
  description:    rasterName,
  fileNamePrefix: rasterName,
  folder:         'kidepo_lcluc',
  region:         aoiExport,
  crs:            EXPORT_CRS,
  scale:          EXPORT_SCALE,
  maxPixels:      1e13,
  fileFormat:     'GeoTIFF',
  formatOptions:  {cloudOptimized: true}
});

Export.table.toDrive({
  collection:     sunFc,
  description:    'kidepo_hls_sun_geometry',
  fileNamePrefix: 'kidepo_hls_sun_geometry',
  folder:         'kidepo_lcluc',
  fileFormat:     'CSV',
  selectors:      ['scene_id', 'year', 'date', 'doy', 'cloud_pct'].concat(SUN_FIELDS)
});


/*
 *  AFTER THE EXPORTS FINISH
 *
 *  Download both into the same folder as the composites, then run
 *      scripts/R/10_terrain_test.R
 *
 *  which asks, in order:
 *
 *   1. Do the candidates occupy a particular TERRAIN NICHE — steeper, higher,
 *      or facing a particular way — or are they spread over terrain the way the
 *      landscape is? An aspect preference is the fingerprint of illumination;
 *      a slope or elevation preference on its own is equally consistent with
 *      ecology, because woody vegetation genuinely prefers hills here.
 *
 *   2. Does the illumination geometry actually DIFFER between the two dates?
 *      cos(i) is computed per pixel for each year from slope, aspect and that
 *      year's mean solar position. If the difference is negligible, as the same
 *      calendar window predicts, illumination cannot produce the change and the
 *      hypothesis is closed on mechanism rather than on a weak coefficient.
 *
 *   3. Net of that, does slope predict dNDVI once protected area and
 *      observation depth are held constant — and does the observation-depth
 *      artefact from script 08 concentrate on steep ground? That is Route 2,
 *      and it is the version of the terrain hypothesis most likely to survive.
 *
 *   4. How much of the 277.5 km2 and the +0.013 dNDVI is left standing.
 */

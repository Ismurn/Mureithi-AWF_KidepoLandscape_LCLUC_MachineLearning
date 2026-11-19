/*
 *  KIDEPO LANDSCAPE — LAND COVER / LAND USE CHANGE (LCLUC)
 *  GEE Script 02 — Build the two matched annual composites
 *
 *  PURPOSE
 *  Stage 1 of the hybrid GEE -> R pipeline. Produce two cloud-masked median
 *  composites, one per comparison date, identical in bands, CRS, resolution and
 *  extent, and export them as GeoTIFFs for classification in R.
 *
 *  RUN 01 AND 01b FIRST. The two years and the window below are the result of
 *  those scripts, not a guess. Change them in one place, here.
 *
 *  WHAT IS EXPORTED, per date — 10 bands, in this order:
 *      1 BLUE   2 GREEN  3 RED  4 NIR  5 SWIR1  6 SWIR2
 *      7 NDVI   8 NDWI   9 NBR
 *     10 clear_obs
 *
 *  ENCODING — read this before loading in R:
 *     dtype              Int16
 *     reflectance        surface reflectance x 10000  (divide by 10000 in R)
 *     indices            index value x 10000          (divide by 10000 in R)
 *     clear_obs          raw count, no scaling
 *     nodata             -9999  (set to NA in R)
 *     CRS                EPSG:32636 (UTM 36N)
 *     resolution         30 m
 *
 *  WHY Int16 rather than Float32: it halves the download for an identical result
 *  once rescaled, and 10000 x reflectance fits comfortably inside Int16's range.
 *
 *  GRID ALIGNMENT (CLAUDE.md caveat 3): both exports use the SAME crs, the SAME
 *  scale and the SAME region. With `crs` + `scale` and no crsTransform, Earth
 *  Engine snaps the grid to the projection origin, so the two rasters land on an
 *  identical grid and terra can stack them without resampling. The R script
 *  still verifies this with compareGeom() rather than trusting it.
 *
 *  ASSETS USED (in users/imnduta1/KidepoAOI-LULC-ML/):
 *      kidepo_pa_union_shp    dissolved PA union (UNBUFFERED — buffered below)
 *
 *  Author: Isaac Mureithi - IMureithi@awf.org
 */


// 
// PHASE 0 — CONFIGURATION: the decisions from scripts 01 and 01b
// 

// 0.1 THE TWO DATES AND THE WINDOW — the only lines that normally change.
var WINDOW_NAME = 'dry_late';
var YEAR_1      = 2019;
var YEAR_2      = 2024;

// 0.2 Assets
var ASSET_ROOT = 'users/imnduta1/KidepoAOI-LULC-ML/';
var aoi_union  = ee.FeatureCollection(ASSET_ROOT + 'kidepo_pa_union_shp');

// 0.3 Extent.
// The classification and reporting extent is the PA union; the composite is cut
// 500 m wider so that focal operations, sieve filters and edge pixels in R have
// valid neighbours instead of running off the edge of the data. Buffering here
// means no second asset upload is needed.
var BUFFER_M   = 500;
var aoiExport  = aoi_union.geometry().buffer(BUFFER_M);

// 0.4 Compositing parameters
var MAX_SCENE_CLOUD = 60;    // %; same threshold the coverage check was run with.
                             // Could be relaxed to ~80 since per-pixel masking
                             // does the real work, but keeping it matches the
                             // evidence in 01.
var EXPORT_CRS   = 'EPSG:32636';
var EXPORT_SCALE = 30;
var NODATA       = -9999;
var SCALE_INT    = 10000;

// 0.5 Seasonal windows — MUST match scripts 01 and 01b.
var WINDOW_DEFS = {
  dry_early: {startMonth: 12, startDay: 1, endMonth: 1, endDay: 31, endYearOffset: 1},
  dry_late:  {startMonth:  2, startDay: 1, endMonth: 3, endDay: 31, endYearOffset: 0}
};

// 0.6 Band mapping.
// HLSL30 (Landsat 8/9) and HLSS30 (Sentinel-2) name the same wavelengths
// DIFFERENTLY. Renaming both to a common scheme is what makes them stackable,
// and it is the single most error-prone part of working with HLS — so it is
// written out explicitly rather than inferred.
//
//   common :  BLUE   GREEN  RED    NIR    SWIR1  SWIR2
//   HLSL30 :  B2     B3     B4     B5     B6     B7
//   HLSS30 :  B2     B3     B4     B8A    B11    B12    (B8A = narrow NIR)
//
// B8A is used rather than B8 because it is the narrow NIR band that HLS
// harmonises against Landsat's B5; B8 is a wider band and is not equivalent.
var BANDS_COMMON = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2'];
var BANDS_L30    = ['B2', 'B3', 'B4', 'B5',  'B6',  'B7'];
var BANDS_S30    = ['B2', 'B3', 'B4', 'B8A', 'B11', 'B12'];

// HLS v2 stores surface reflectance as Int16 with a 0.0001 scale factor.
var HLS_SR_SCALE = 0.0001;


// 
// PHASE 1 — HELPERS
// 

/** 1.1 Leap-safe window dates. Identical construction to scripts 01 and 01b. */
function windowFor(windowName, year) {
  var w     = WINDOW_DEFS[windowName];
  var start = ee.Date.fromYMD(year, w.startMonth, w.startDay);
  var end   = ee.Date
    .fromYMD(ee.Number(year).add(w.endYearOffset), w.endMonth, w.endDay)
    .advance(1, 'day');
  return {start: start, end: end};
}

/**
 * 1.2 Fmask decoding — identical to script 01, deliberately.
 * Bits, LSB first: 0 cirrus | 1 cloud | 2 adjacent-to-cloud | 3 cloud shadow
 *                  4 snow/ice | 5 water | 6-7 aerosol
 * Water (bit 5) is NOT masked: permanent water is a land-cover class to be
 * classified, and masking it would delete it from the map.
 */
function clearBits(img) {
  var fmask = img.select('Fmask');
  return fmask.bitwiseAnd(1 << 0).eq(0)
    .and(fmask.bitwiseAnd(1 << 1).eq(0))
    .and(fmask.bitwiseAnd(1 << 2).eq(0))
    .and(fmask.bitwiseAnd(1 << 3).eq(0));
}

/**
 * 1.3 One harmonised, masked, date-filtered product collection.
 * Date filtering precedes map() so the whole archive is never touched.
 */
function prepProduct(collectionId, bandsIn, start, end) {
  return ee.ImageCollection(collectionId)
    .filterBounds(aoiExport)
    .filterDate(start, end)
    .filterMetadata('CLOUD_COVERAGE', 'less_than', MAX_SCENE_CLOUD)
    .map(function (img) {
      var clear = clearBits(img);
      return img.select(bandsIn, BANDS_COMMON)
        .multiply(HLS_SR_SCALE)
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

/**
 * 1.5 Spectral indices.
 * NDVI  vegetation vigour — the main land-cover discriminator.
 * NDWI  McFeeters (GREEN, NIR) — open water, which Fmask leaves unmasked.
 * NBR   (NIR, SWIR2) — burn sensitivity. Included because Karamoja burns inside
 *       the dry-season windows; it gives the R classifier a direct handle on
 *       burn scars instead of leaving them to masquerade as bare ground.
 */
function addIndices(img) {
  return img
    .addBands(img.normalizedDifference(['NIR',   'RED'  ]).rename('NDVI'))
    .addBands(img.normalizedDifference(['GREEN', 'NIR'  ]).rename('NDWI'))
    .addBands(img.normalizedDifference(['NIR',   'SWIR2']).rename('NBR'));
}

/**
 * 1.6 The composite for one date.
 * A MEDIAN is used rather than a mean: it is robust to the residual cloud and
 * shadow that Fmask misses, which matters more than the slight loss of
 * efficiency. The clear-observation count rides along as a band so that R can
 * flag or mask thinly-sampled pixels rather than trusting every pixel equally.
 */
function buildComposite(windowName, year) {
  var col = hlsFor(windowName, year);

  var reflectance = col.select(BANDS_COMMON).median();
  var withIndices = addIndices(reflectance);

  var clearObs = col.select('RED').count().unmask(0).rename('clear_obs');

  return withIndices
    .multiply(SCALE_INT)          // reflectance and indices -> Int16 domain
    .addBands(clearObs)           // count stays unscaled
    .unmask(NODATA)               // explicit sentinel, not a silent 0
    .toInt16()
    .clip(aoiExport)
    .set({
      window:      windowName,
      year:        year,
      scale_int:   SCALE_INT,
      nodata:      NODATA,
      n_scenes:    col.size()
    });
}


// 
// PHASE 2 — BUILD BOTH COMPOSITES
// 

var comp1 = buildComposite(WINDOW_NAME, YEAR_1);
var comp2 = buildComposite(WINDOW_NAME, YEAR_2);

print('Composite ' + YEAR_1 + ' (' + WINDOW_NAME + '): bands', comp1.bandNames());
print('Composite ' + YEAR_2 + ' (' + WINDOW_NAME + '): bands', comp2.bandNames());

// A band-count mismatch would silently break the stack in R. Both are printed so
// the two can be compared at a glance before the exports are started.


// 
// PHASE 3 — VISUAL CHECK BEFORE EXPORTING
// 
// Export tasks take time and quota. Look at the composites first: a false-colour
// pair makes gross failures obvious (a black stripe, a cloud bank that survived
// masking, a seam between sensors).

var visFalse = {bands: ['SWIR2', 'NIR', 'GREEN'], min: 0, max: 3500};
var visNbr   = {bands: ['NBR'], min: -3000, max: 6000,
                palette: ['#67001f', '#d7301f', '#fdae61', '#f7f7f7', '#1b7837']};

Map.centerObject(aoiExport, 8);
Map.addLayer(comp1, visFalse, YEAR_1 + ' false colour (SWIR2/NIR/GREEN)');
Map.addLayer(comp2, visFalse, YEAR_2 + ' false colour (SWIR2/NIR/GREEN)');
Map.addLayer(comp1, visNbr, YEAR_1 + ' NBR (low = burnt)', false);
Map.addLayer(comp2, visNbr, YEAR_2 + ' NBR (low = burnt)', false);
Map.addLayer(ee.Image().paint(aoi_union, 0, 2), {palette: ['#E07B39']}, 'PA union');


// 
// PHASE 4 — EXPORT
// 
// Roughly 8,500 km2 at 30 m is about 9.4 million pixels; 10 Int16 bands is
// roughly 190 MB per date. If Drive returns multiple tiles, raise maxPixels or
// set fileDimensions, and mosaic them in R with terra::vrt().

function exportComposite(img, year) {
  Export.image.toDrive({
    image:          img,
    description:    'kidepo_hls_' + WINDOW_NAME + '_' + year + '_utm36n_30m',
    fileNamePrefix: 'kidepo_hls_' + WINDOW_NAME + '_' + year + '_utm36n_30m',
    folder:         'kidepo_lcluc',
    region:         aoiExport,     // SAME region for both — see grid alignment above
    crs:            EXPORT_CRS,
    scale:          EXPORT_SCALE,
    maxPixels:      1e13,
    fileFormat:     'GeoTIFF',
    formatOptions:  {cloudOptimized: true}
  });
}

exportComposite(comp1, YEAR_1);
exportComposite(comp2, YEAR_2);


/*
 *  AFTER THE EXPORTS FINISH
 *
 *  1. Download both GeoTIFFs into outputs/composites/ (gitignored).
 *  2. In R: rescale by 10000, set -9999 to NA, and confirm with compareGeom()
 *     that the two rasters share a grid before stacking them.
 *  3. Check the clear_obs band in both. Pixels below the threshold used in
 *     script 01 are weakly supported, and any change reported there needs the
 *     caveat — or masking out.
 *  4. Then classify: k-means first for a fast look at the class structure,
 *     Random Forest once the scheme is worth labelling for.
 *
 *  A caution carried from CLAUDE.md (caveat 2): the ARSET workflow trains one RF
 *  on date 1 and applies it to both dates, assuming spectral signatures are
 *  stable across dates. In a fire-driven semi-arid landscape that assumption is
 *  shaky. Check it with the RF proximity / training-similarity diagnostic before
 *  trusting the cross-date prediction — and if it fails, classify each date
 *  separately with its own training data.
 */
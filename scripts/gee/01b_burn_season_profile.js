/*
 *  KIDEPO LANDSCAPE — LAND COVER / LAND USE CHANGE (LCLUC)
 *  GEE Script 01b — Burn-season profile
 *
 *  PURPOSE
 *  Settle the seasonal-window decision with data instead of assertion.
 *
 *  Script 01 shows that the late dry season (Feb-Mar) has the clearest skies.
 *  But Karamoja burns in exactly that season, and a fresh burn scar classifies as
 *  bare/barren. If the two comparison dates carry different amounts of burning,
 *  the change map measures fire, not land-cover change (CLAUDE.md caveat 1).
 *
 *  This script quantifies that risk from MODIS burned area (MCD64A1, 500 m,
 *  monthly since 2000) and answers three questions:
 *
 *    1. WHEN does the AOI burn?            -> monthly table, every year x month
 *    2. HOW MUCH burns inside each window? -> per window-year, two measures:
 *         burn_in_window_pct   burned during the window itself
 *         burn_cumulative_pct  burned from the start of the fire season up to
 *                              the window's end. THIS is the number that matters,
 *                              because a scar is still visible for weeks after
 *                              the fire, so a composite sees the accumulated
 *                              burning, not just what burnt while it was being
 *                              collected.
 *    3. DOES IT DIFFER between the two candidate years? A pair that differs
 *       sharply in burnt fraction is a pair that will produce fire-driven change.
 *
 *  MCD64A1 is 500 m against HLS's 30 m, so treat these as landscape proportions
 *  for choosing a window, not as a per-pixel burn mask.
 *
 *  ASSETS USED (in users/imnduta1/KidepoAOI-LULC-ML/):
 *      kidepo_pa_union_shp    dissolved PA union
 *      kidepo_pa_parts_shp    individual PAs, carrying the `country` attribute
 *
 *  OUTPUT: two CSVs in Drive/kidepo_lcluc.
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

// 0.2 Period
var YEAR_START = 2016;
var YEAR_END   = 2026;

// 0.3 Seasonal windows — MUST match scripts/gee/01_hls_coverage_check.js.
// Kept as a literal copy so each script runs standalone in the Code Editor; if
// one is edited, edit the other.
var WINDOW_DEFS = {
  dry_early: {startMonth: 12, startDay: 1, endMonth: 1, endDay: 31,
              endYearOffset: 1, lastStartYear: 2025},
  dry_late:  {startMonth:  2, startDay: 1, endMonth: 3, endDay: 31,
              endYearOffset: 0, lastStartYear: 2026}
};
var WINDOW_NAMES = ['dry_early', 'dry_late'];

// 0.4 Fire-season start.
// Assumed to open on 1 November, before either window. This is an ASSUMPTION:
// check it against the monthly table this script exports, and move it if the
// monthly profile shows burning starting earlier.
var SEASON_START_MONTH = 11;

// 0.5 Burned-area source. MCD64A1 BurnDate: 0 = unburned, 1-366 = day of burn,
// negative = unmapped. So `.gt(0)` isolates burnt pixels under every encoding.
//
// LATENCY WARNING — do not misread the tail of the table. MCD64A1 is published
// a couple of months in arrears, and a month with no product yet returns 0%
// burnt here, which is indistinguishable from a month that genuinely did not
// burn. Treat the most recent months of the monthly table as missing, not as
// zero, and cross-check against the fire seasons of earlier years before using
// them for anything.
var BURN_COLL  = 'MODIS/061/MCD64A1';
var BURN_SCALE = 500;   // m; native MCD64A1 resolution


// 
// PHASE 1 — HELPERS
// 

var burn = ee.ImageCollection(BURN_COLL).select('BurnDate');

/**
 * 1.1 Binary "burnt at any point in this period" image.
 * max() over the monthly 0/1 layers, so a pixel that burns in any month of the
 * period counts once — this is burnt AREA, not a count of burn events.
 * The ee.Algorithms.If guard handles a period with no MCD64A1 products.
 */
function burntImage(start, end) {
  var col = burn.filterDate(start, end);
  return ee.Image(
      ee.Algorithms.If(col.size().gt(0),
                       col.map(function (i) { return i.gt(0).unmask(0); }).max(),
                       ee.Image.constant(0))
    )
    .unmask(0)
    .rename('burnt')
    .clip(aoi);
}

/** 1.2 Percentage of a sub-area burnt during a period. */
function burntPct(start, end, geom) {
  return ee.Number(
    burntImage(start, end).reduceRegion({
      reducer:    ee.Reducer.mean(),
      geometry:   geom,
      scale:      BURN_SCALE,
      maxPixels:  1e10,
      bestEffort: true
    }).get('burnt')
  ).multiply(100);
}

/** 1.3 Window dates, identical in construction to script 01 (leap-safe). */
function windowFor(windowName, year) {
  var w     = WINDOW_DEFS[windowName];
  var start = ee.Date.fromYMD(year, w.startMonth, w.startDay);
  var end   = ee.Date
    .fromYMD(ee.Number(year).add(w.endYearOffset), w.endMonth, w.endDay)
    .advance(1, 'day');
  return {start: start, end: end};
}

/**
 * 1.4 The 1 November immediately BEFORE a window starts.
 * dry_early starts in December, so its season opens in the same year;
 * dry_late starts in February, so its season opened the previous November.
 */
function seasonStart(windowName, year) {
  var w      = WINDOW_DEFS[windowName];
  var offset = (w.startMonth >= SEASON_START_MONTH) ? 0 : -1;
  return ee.Date.fromYMD(ee.Number(year).add(offset), SEASON_START_MONTH, 1);
}

var SUBAREAS = {
  all:         aoi,
  uganda:      aoi_parts.filter(ee.Filter.eq('country', 'Uganda')).geometry(),
  south_sudan: aoi_parts.filter(ee.Filter.eq('country', 'South Sudan')).geometry()
};


// 
// PHASE 2 — QUESTION 1: WHEN DOES THE AOI BURN?
// 
// One row per year x calendar month. Aggregate in R to get the climatology; kept
// disaggregated here so an unusual fire year stays visible rather than averaged
// away.

var monthlyRows = [];
for (var y = YEAR_START; y <= YEAR_END; y++) {
  for (var m = 1; m <= 12; m++) {
    monthlyRows.push(ee.Feature(null, {
      year:          y,
      month:         m,
      burnt_pct_all: burntPct(ee.Date.fromYMD(y, m, 1),
                              ee.Date.fromYMD(y, m, 1).advance(1, 'month'),
                              SUBAREAS.all)
    }));
  }
}
var monthlyBurn = ee.FeatureCollection(monthlyRows);


// 
// PHASE 3 — QUESTION 2: HOW MUCH BURNS INSIDE EACH WINDOW?
// 

function burnStatsFor(windowName, year) {
  var win    = windowFor(windowName, year);
  var season = seasonStart(windowName, year);

  return ee.Feature(null, {
    window:              windowName,
    year:                year,
    window_start:        win.start.format('YYYY-MM-dd'),
    window_end:          win.end.advance(-1, 'day').format('YYYY-MM-dd'),
    season_start:        season.format('YYYY-MM-dd'),
    // Burnt during the window itself.
    burn_in_window_pct:  burntPct(win.start, win.end, SUBAREAS.all),
    // Burnt from the season opening to the window's end — what a composite sees.
    burn_cumulative_pct: burntPct(season, win.end, SUBAREAS.all),
    burn_cum_uga_pct:    burntPct(season, win.end, SUBAREAS.uganda),
    burn_cum_ssd_pct:    burntPct(season, win.end, SUBAREAS.south_sudan)
  });
}

var windowRows = [];
WINDOW_NAMES.forEach(function (nm) {
  for (var yy = YEAR_START; yy <= WINDOW_DEFS[nm].lastStartYear; yy++) {
    windowRows.push(burnStatsFor(nm, yy));
  }
});
var windowBurn = ee.FeatureCollection(windowRows);


// 
// PHASE 4 — MAPS: WHERE AND HOW OFTEN IT BURNS
// 
// Burn frequency over the study period, for each window. Areas that burn almost
// every year are a persistent classification hazard, not a one-off: they will
// read as bare in BOTH dates, which is at least consistent. Areas that burn
// sporadically are the ones that create false change.

var MAP_WINDOW = 'dry_late';

function burnFrequency(windowName, yFrom, yTo) {
  var imgs = [];
  for (var yr = yFrom; yr <= yTo; yr++) {
    var w = windowFor(windowName, yr);
    imgs.push(burntImage(seasonStart(windowName, yr), w.end));
  }
  return ee.ImageCollection(imgs).sum().rename('burn_years').clip(aoi);
}

var freq = burnFrequency(MAP_WINDOW, YEAR_START, WINDOW_DEFS[MAP_WINDOW].lastStartYear);
var nYears = WINDOW_DEFS[MAP_WINDOW].lastStartYear - YEAR_START + 1;

Map.centerObject(aoi, 8);
Map.addLayer(freq,
             {min: 0, max: nYears, palette: ['#f7f7f7', '#fdae61', '#d7301f', '#67001f']},
             'Burn frequency, years burnt by end of ' + MAP_WINDOW);
Map.addLayer(aoi_parts.style({color: '#000000', fillColor: '00000000', width: 1}),
             {}, 'PA polygons');


// 
// PHASE 5 — EXPORTS
// 

Export.table.toDrive({
  collection:  monthlyBurn,
  description: 'kidepo_burn_monthly_' + YEAR_START + '_' + YEAR_END,
  folder:      'kidepo_lcluc',
  fileFormat:  'CSV',
  selectors:   ['year', 'month', 'burnt_pct_all']
});

Export.table.toDrive({
  collection:  windowBurn,
  description: 'kidepo_burn_by_window_' + YEAR_START + '_' + YEAR_END,
  folder:      'kidepo_lcluc',
  fileFormat:  'CSV',
  selectors:   ['window', 'year', 'window_start', 'window_end', 'season_start',
                'burn_in_window_pct', 'burn_cumulative_pct',
                'burn_cum_uga_pct', 'burn_cum_ssd_pct']
});


/*
 *  HOW TO READ THE RESULT — and what it decides
 *
 *  Read burn_cumulative_pct, not burn_in_window_pct: the composite sees every
 *  scar laid down since the season opened, not only the fires that happened to
 *  burn while the imagery was being collected.
 *
 *  Then:
 *   - If dry_early carries substantially less cumulative burning than dry_late,
 *     that is the argument for accepting dry_early's extra cloud. Check in
 *     script 01 that both candidate years still clear MIN_CLEAR_OBS in that
 *     window; if they do, dry_early is the better window.
 *   - If both windows burn similarly, the fire argument is neutral and dry_late
 *     wins on clear observations.
 *   - Whatever the level, what matters most is that the TWO CHOSEN YEARS are
 *     similar to each other. A pair at 12% and 13% burnt is far safer than a
 *     pair at 5% and 25%, even though the second pair has a lower minimum.
 *   - If no pair is well matched, the fallback is to carry a burn class through
 *     the classification rather than pretend it is not there, and to say so in
 *     the write-up.
 */
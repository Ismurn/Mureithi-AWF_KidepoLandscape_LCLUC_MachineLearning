# =============================================================================
#  KIDEPO LANDSCAPE — LCLUC
#  Script 10 — Is the candidate change an artefact of terrain illumination?
#
#  WHY
#  Script 09 left three explanations for the 277.5 km2 of candidate change:
#  genuine woody thickening, differential terrain illumination, or real
#  land-use change. The candidates concentrate in the small HILL reserves —
#  Timu 26.5% of its area, Rom 23.5%, Imatong 7.2% — against 0.5% for lowland
#  Kidepo Valley. Rough ground is exactly where an illumination artefact would
#  live, so this has to be settled before any of that area is reported.
#
#  THE MECHANISM MATTERS MORE THAN THE COEFFICIENT
#  Both composites use the SAME calendar window, 1 Feb to 31 Mar. The sun is in
#  nearly the same place in both. Illumination sets the LEVEL of reflectance on
#  a slope, so it very largely CANCELS in the difference between two
#  same-season composites. The hypothesis therefore needs something to differ
#  between the years, and there are only two candidates:
#
#    ROUTE 1  TIMING — if one year's clear looks fall earlier or later in the
#             window, the sun sat at a different declination, cos(i) differs,
#             and the difference is largest on steep slopes.
#
#    ROUTE 2  SHADOW MASKING — Fmask cannot tell terrain shadow from cloud
#             shadow, so shaded slopes lose observations, and WHICH looks
#             survive depends on how many there were. Script 08 already found
#             28% of the greening attributable to observation depth. If that
#             artefact concentrates on steep ground, terrain is its mechanism.
#
#  Route 1 has a PREDICTED NULL, and saying so in advance is the point: if the
#  two years were observed at nearly the same time of season then illumination
#  cannot produce a change, and the hypothesis dies on mechanism rather than on
#  a weak coefficient fished out afterwards. Route 2 is the version most likely
#  to survive, and it is not really about illumination at all — it is the
#  script 08 confound with a spatial pattern.
#
#  WHAT THIS SCRIPT WILL NOT DO
#  It cannot distinguish genuine woody thickening from land-use change. That
#  needs high-resolution imagery or a multi-year series. It can only remove, or
#  fail to remove, the artefact.
#
#  NEEDS  scripts/gee/01d_export_terrain.js run and downloaded:
#           kidepo_terrain_dry_late_utm36n_30m.tif   (required)
#           kidepo_hls_sun_geometry.csv              (optional, see Phase 2)
#
#  HOW TO RUN — from the repo root, in PowerShell. NOT in Git Bash: the MSYS
#  mingw64 DLLs shadow GDAL and sf segfaults.
#      & "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" scripts\R\10_terrain_test.R
#
#  OUTPUTS -> outputs/terrain/
#
#  Author: Isaac Mureithi - IMureithi@awf.org
# =============================================================================

source("scripts/R/00_functions.R")
pacman::p_load(tidyterra, ggspatial, ggnewscale)

terra::terraOptions(progress = 0)
set.seed(11)

path.out.terr <- file.path("outputs", "terrain")
dir.create(path.out.terr, recursive = TRUE, showWarnings = FALSE)

n.samp      <- 400000   # sample size for every regression below
slope.flat  <- 2        # degrees; below this, aspect is meaningless
slope.steep <- 10       # degrees; the "steep ground" cut for Route 2

## The HLS overpass time is `overpass.solar.hour` in 00_functions.R, beside the
## illumination helpers it belongs to.

cat("\n================================================================\n")
cat(" KIDEPO LCLUC — 10 Terrain illumination test\n")
cat(" run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n")


# =============================================================================
#  PHASE 0 — THE TERRAIN RASTER, AND A REFUSAL TO TRUST ITS SCALE FACTORS
# =============================================================================
#  The export stores slope and aspect x10 and elevation and DOY unscaled. Mixed
#  scale factors in one file are how this project already lost a whole set of
#  reflectance bands (see 02_build_composites.js). So the scales live in ONE
#  lookup below, and every band is then checked against the range it is
#  physically allowed to occupy. A scale applied twice, or not at all, cannot
#  pass all four checks.

terr.bands <- c("elevation", "slope", "aspect", "obs.doy.a", "obs.doy.b")

terr.scale <- c(elevation = 1, slope = 10, aspect = 10,
                obs.doy.a = 1, obs.doy.b = 1)

## Physically possible ranges, in TRUE units, for this landscape and window.
terr.limit <- list(
  elevation = c(300,  4500),   # Kidepo Valley floor to the Imatong massif
  slope     = c(0,      90),
  aspect    = c(0,     360),
  obs.doy.a = c(32,     91),   # 1 Feb to 31 Mar
  obs.doy.b = c(32,     91)
)

f.terr <- file.path(path.gee.out, "kidepo_terrain_dry_late_utm36n_30m.tif")
f.sun  <- file.path(path.gee.out, "kidepo_hls_sun_geometry.csv")

if (!file.exists(f.terr)) {
  stop(glue(
    "Terrain raster not found:\n  {f.terr}\n\n",
    "Run scripts/gee/01d_export_terrain.js in the Earth Engine Code Editor, ",
    "let both exports finish in Tasks, then download them into\n  {path.gee.out}\n",
    "The CSV is optional; the GeoTIFF is not."))
}

terr <- terra::rast(f.terr)

if (terra::nlyr(terr) != length(terr.bands)) {
  stop(glue("{basename(f.terr)} has {terra::nlyr(terr)} bands; expected ",
            "{length(terr.bands)} ({paste(terr.bands, collapse = ', ')}). ",
            "The export band order may have changed."))
}
names(terr) <- terr.bands
terra::NAflag(terr) <- comp.nodata

## The raster stays INTEGER. Dividing 5 bands of ~32 million cells would
## materialise half a gigabyte of float temp files for no gain: the scale is
## applied to the SAMPLE in Phase 1, exactly as the composites are handled.
##
## compute = TRUE forces a full read of all five bands, which takes a minute or
## so. That is the price of the check and it is worth paying once.
mm <- terra::minmax(terr, compute = TRUE)

cat("\n[0] BAND RANGES IN TRUE UNITS (after the lookup scale)\n\n")
range.tbl <- tibble(
  band    = terr.bands,
  scale   = as.numeric(terr.scale[terr.bands]),
  obs.min = as.numeric(mm[1, terr.bands]) / as.numeric(terr.scale[terr.bands]),
  obs.max = as.numeric(mm[2, terr.bands]) / as.numeric(terr.scale[terr.bands]),
  ok.min  = map_dbl(terr.bands, ~ terr.limit[[.x]][1]),
  ok.max  = map_dbl(terr.bands, ~ terr.limit[[.x]][2])
) %>%
  mutate(pass = obs.min >= ok.min & obs.max <= ok.max)

range.tbl %>%
  mutate(across(where(is.numeric), ~round(.x, 2))) %>%
  as.data.frame() %>% print(row.names = FALSE)

if (any(!range.tbl$pass)) {
  bad <- range.tbl %>% filter(!pass) %>% pull(band) %>% paste(collapse = ", ")
  stop(glue(
    "\nBand(s) outside their physically possible range: {bad}.\n",
    "Either the export scale factors differ from `terr.scale` in this script, ",
    "or the band order changed. Fix the mismatch before reading anything ",
    "below as a result — a wrong scale here would silently rescale every ",
    "illumination statistic in this script."))
}
say("\n    All bands land inside their possible ranges. Scales confirmed.")


# =============================================================================
#  PHASE 1 — ASSEMBLE THE SAMPLE
# =============================================================================
#  One sample serves every test, so the same pixels answer every question and
#  the numbers are mutually comparable. Latitude comes along because solar
#  zenith depends on it and the AOI spans about a degree.

comp <- load.composite.pair()
yr.a <- comp$years[1]
yr.b <- comp$years[2]

aoi <- load.aoi()
aoi.union.v <- terra::vect(aoi$union)
aoi.parts.v <- terra::vect(aoi$parts)
mask.to.aoi <- function(r) terra::mask(terra::crop(r, aoi.union.v), aoi.union.v)

f.cand <- file.path(path.out.cls, "change_strong_clusterfree.tif")
if (!file.exists(f.cand)) stop("Candidate raster missing — run 03c_k_sensitivity.R first.")

cand <- mask.to.aoi(terra::rast(f.cand))
names(cand) <- "candidate"

## compareGeom rather than trust: the terrain export used the same crs, scale
## and region as the composites, which SHOULD align it, but a half-pixel offset
## would quietly attach the wrong slope to every pixel.
if (!isTRUE(terra::compareGeom(comp$a, terr, stopOnError = FALSE, messages = FALSE))) {
  stop(glue("The terrain raster is NOT on the composite grid. Re-export it ",
            "from 01d with an identical `region`, `crs` and `scale`."))
}

terr.m <- mask.to.aoi(terr)

ndvi.a <- mask.to.aoi(comp$a[["NDVI"]])
ndvi.b <- mask.to.aoi(comp$b[["NDVI"]])
obs.a  <- mask.to.aoi(comp$a[["clear_obs"]])
obs.b  <- mask.to.aoi(comp$b[["clear_obs"]])

pa.r <- terra::rasterize(aoi.parts.v, ndvi.a, field = "pa.label")
names(pa.r) <- "pa"

px.km2 <- prod(terra::res(cand)) / 1e6

## Read the headline candidate area off the raster rather than hard-coding the
## 277.5 km2 figure from 03c, which would go stale the moment 03c is re-run.
cand.km2.total <- terra::global(cand, "sum", na.rm = TRUE)[1, 1] * px.km2

## NDVI arrives named "NDVI" from both dates and clear_obs likewise, so this
## stack carries duplicate layer names until the rename on the next line. That
## is harmless for terra but it broke as_tibble() in script 09 — rename first,
## sample second, and the problem never arises.
all.r <- c(cand, ndvi.a, ndvi.b, obs.a, obs.b, terr.m, pa.r)
names(all.r) <- c("candidate", "ndvi.a", "ndvi.b", "obs.a", "obs.b",
                  terr.bands, "pa")

samp <- terra::spatSample(all.r, n.samp, method = "random",
                          na.rm = TRUE, xy = TRUE, warn = FALSE) %>%
  as_tibble() %>%
  drop_na()

## Latitude, for the solar position. UTM metres -> geographic degrees.
lat <- samp %>%
  select(x, y) %>%
  sf::st_as_sf(coords = c("x", "y"), crs = crs.utm) %>%
  sf::st_transform(crs.geo) %>%
  sf::st_coordinates()

samp <- samp %>%
  mutate(
    lat        = lat[, "Y"],
    ## Scale factors applied HERE and nowhere else.
    elevation  = elevation / terr.scale[["elevation"]],
    slope.deg  = slope     / terr.scale[["slope"]],
    aspect.deg = aspect    / terr.scale[["aspect"]],
    obs.doy.a  = obs.doy.a / terr.scale[["obs.doy.a"]],
    obs.doy.b  = obs.doy.b / terr.scale[["obs.doy.b"]],
    ndvi.a     = to.reflectance(ndvi.a),
    ndvi.b     = to.reflectance(ndvi.b),
    dndvi      = ndvi.b - ndvi.a,
    obs.gap.px = obs.b - obs.a,
    doy.gap    = obs.doy.b - obs.doy.a
  ) %>%
  ## Same support filter as 08, so the two scripts' numbers are comparable.
  filter(obs.a >= comp.min.obs, obs.b >= comp.min.obs)

say(glue("
[1] SAMPLE
    pixels sampled, well supported   : {format(nrow(samp), big.mark = ',')}
    candidates in the sample         : {format(sum(samp$candidate), big.mark = ',')} ({round(100 * mean(samp$candidate), 2)}%)
    mean dNDVI                       : {sprintf('%+.5f', mean(samp$dndvi))}
    mean slope                       : {round(mean(samp$slope.deg), 2)} degrees
    slope above {slope.steep} degrees           : {round(100 * mean(samp$slope.deg > slope.steep), 1)}%
    elevation range                  : {round(min(samp$elevation))}-{round(max(samp$elevation))} m"))


# =============================================================================
#  PHASE 2 — WHERE WAS THE SUN?
# =============================================================================
#  Two independent routes to the solar position, because neither is beyond
#  doubt and they should agree:
#
#   RECORDED   the per-scene angles in the HLS metadata, exported to CSV by
#              01d. Authoritative where present, but it is a per-SCENE mean, so
#              it cannot vary with which pixels were clear.
#
#   COMPUTED   Cooper's declination and the standard hour-angle formula,
#              evaluated at each pixel's OWN mean observation day and latitude.
#              Approximate in absolute terms — it assumes a fixed overpass at
#              local solar time 10.25 and ignores the equation of time, worth
#              up to about 4 degrees of hour angle.
#
#  The absolute error in COMPUTED does not matter much here, and that is worth
#  being explicit about: a wrong overpass time is wrong by the SAME amount in
#  both years, so it cancels almost exactly in the difference, which is the only
#  quantity this script draws a conclusion from. Absolute cos(i) is indicative;
#  delta cos(i) is the measurement.

## solar.declination(), solar.position(), cos.illum(), deg2rad() and rad2deg()
## all live in 00_functions.R section 3.9, so that
## scripts/R/tests/test_solar_geometry.R validates the same code this script
## runs rather than a copy of it. The overpass assumption is the constant
## `overpass.solar.hour` there.

sun.a <- solar.position(samp$obs.doy.a, samp$lat)
sun.b <- solar.position(samp$obs.doy.b, samp$lat)

samp <- samp %>%
  mutate(
    zen.a  = sun.a$zenith,  azi.a = sun.a$azimuth,
    zen.b  = sun.b$zenith,  azi.b = sun.b$azimuth,
    cosi.a = cos.illum(slope.deg, aspect.deg, zen.a, azi.a),
    cosi.b = cos.illum(slope.deg, aspect.deg, zen.b, azi.b),
    dcosi  = cosi.b - cosi.a
  )

## The recorded angles, if the CSV made it and the property names were right.
sun.recorded <- NULL
if (file.exists(f.sun)) {
  raw <- suppressWarnings(readr::read_csv(f.sun, show_col_types = FALSE))
  zcol <- intersect("MEAN_SUN_ZENITH_ANGLE",  names(raw))
  acol <- intersect("MEAN_SUN_AZIMUTH_ANGLE", names(raw))
  have.keys <- all(c("year", "doy") %in% names(raw))
  if (have.keys && length(zcol) == 1 && length(acol) == 1 &&
      any(!is.na(raw[[zcol]])) && any(!is.na(raw[[acol]]))) {
    sun.recorded <- raw %>%
      filter(!is.na(.data[[zcol]])) %>%
      group_by(year) %>%
      summarise(scenes       = n(),
                mean.doy     = mean(doy, na.rm = TRUE),
                mean.zenith  = mean(.data[[zcol]], na.rm = TRUE),
                mean.azimuth = mean(.data[[acol]], na.rm = TRUE),
                .groups = "drop")
  }
}

cat("\n[2] SOLAR GEOMETRY\n\n")

computed <- tibble(
  year         = c(yr.a, yr.b),
  mean.doy     = c(mean(samp$obs.doy.a), mean(samp$obs.doy.b)),
  mean.zenith  = c(mean(samp$zen.a),     mean(samp$zen.b)),
  mean.azimuth = c(mean(samp$azi.a),     mean(samp$azi.b)),
  source       = "computed"
)

if (is.null(sun.recorded)) {
  say("    No usable recorded angles. Either the CSV was not downloaded, or
    Earth Engine's HLS ingestion does not expose MEAN_SUN_ZENITH_ANGLE under
    that name. Falling back to the computed angles alone — which is adequate,
    because the conclusion rests on the DIFFERENCE between years and the
    assumptions are common to both. Check the property names printed by 01d
    if you want the recorded angles instead.")
  sun.cmp <- computed
} else {
  sun.cmp <- bind_rows(computed, sun.recorded %>% mutate(source = "recorded"))
  dz <- abs(computed$mean.zenith - sun.recorded$mean.zenith[match(computed$year, sun.recorded$year)])
  say(glue("
    Recorded and computed zenith differ by {round(max(dz, na.rm = TRUE), 1)} degrees at most.
    Under about 3 degrees means the computed angles are sound and the
    per-pixel version below can be used with confidence. A larger gap means
    the overpass-time assumption is off — it still largely cancels in the
    difference, but say so when reporting."))
}

sun.cmp %>% mutate(across(where(is.numeric), ~round(.x, 2))) %>%
  as.data.frame() %>% print(row.names = FALSE)

doy.shift <- mean(samp$obs.doy.b) - mean(samp$obs.doy.a)

say(glue("
    ROUTE 1 PRECONDITION — did the two years get observed at the same time of
    season?

    mean observation day, {yr.a}      : {round(mean(samp$obs.doy.a), 1)} (day of year)
    mean observation day, {yr.b}      : {round(mean(samp$obs.doy.b), 1)}
    shift                        : {sprintf('%+.1f', doy.shift)} days
    per-pixel shift, 5th-95th    : {sprintf('%+.1f', quantile(samp$doy.gap, 0.05))} to {sprintf('%+.1f', quantile(samp$doy.gap, 0.95))} days
    mean solar zenith shift      : {sprintf('%+.2f', mean(samp$zen.b) - mean(samp$zen.a))} degrees

    mean delta cos(i)            : {sprintf('%+.5f', mean(samp$dcosi))}
    sd of delta cos(i)           : {round(sd(samp$dcosi), 5)}
    on slopes above {slope.steep} degrees   : mean {sprintf('%+.5f', mean(samp$dcosi[samp$slope.deg > slope.steep]))}, sd {round(sd(samp$dcosi[samp$slope.deg > slope.steep]), 5)}

    cos(i) itself runs 0-1, so a mean shift of {sprintf('%.4f', abs(mean(samp$dcosi)))} is a
    {round(100 * abs(mean(samp$dcosi)), 2)}% change in how hard the sun strikes the average slope."))


# =============================================================================
#  PHASE 3 — DO THE CANDIDATES OCCUPY A TERRAIN NICHE?
# =============================================================================
#  ASPECT is the diagnostic one. Slope and elevation preferences are ambiguous:
#  woody vegetation genuinely prefers hills in Karamoja, so "candidates are
#  steeper" is equally consistent with ecology and with artefact. But there is
#  no ecological reason for regrowth to prefer, say, west-facing ground over
#  east-facing at the same slope and elevation — whereas illumination has
#  exactly that signature, because the sun is in the east at a mid-morning
#  overpass and east-facing slopes are lit while west-facing ones are not.
#
#  Aspect is undefined on flat ground, where the export writes 0 and it would
#  read as due north. Everything below excludes slopes under `slope.flat`.

terrain.cmp <- samp %>%
  group_by(candidate) %>%
  summarise(n = n(),
            median.slope  = median(slope.deg),
            mean.slope    = mean(slope.deg),
            pct.steep     = 100 * mean(slope.deg > slope.steep),
            median.elev   = median(elevation),
            .groups = "drop")

cat("\n[3] TERRAIN OF CANDIDATES vs EVERYWHERE ELSE\n\n")
terrain.cmp %>% mutate(across(where(is.numeric), ~round(.x, 2))) %>%
  as.data.frame() %>% print(row.names = FALSE)

oct.levels <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")

asp <- samp %>%
  filter(slope.deg >= slope.flat) %>%
  mutate(octant = oct.levels[floor(((aspect.deg + 22.5) %% 360) / 45) + 1],
         octant = factor(octant, levels = oct.levels))

asp.tbl <- asp %>%
  count(octant, candidate) %>%
  pivot_wider(names_from = candidate, values_from = n,
              names_prefix = "n.", values_fill = 0) %>%
  mutate(pct.cand       = 100 * n.1 / sum(n.1),
         pct.background = 100 * n.0 / sum(n.0),
         enrichment     = pct.cand / pct.background)

cat(glue("\n    Aspect, slopes above {slope.flat} degrees only ",
         "({format(nrow(asp), big.mark = ',')} pixels)\n\n"))
asp.tbl %>% mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

asp.spread <- max(asp.tbl$enrichment) / min(asp.tbl$enrichment)
asp.worst  <- asp.tbl$octant[which.max(asp.tbl$enrichment)]

## Is the aspect preference the one illumination would produce? At a
## mid-morning overpass the sun sits in the east-southeast, so the LIT slopes
## face east and the shaded ones face west. An artefact driven by illumination
## should prefer one of those two groups, not a random octant.
lit.oct   <- c("NE", "E", "SE")
shade.oct <- c("SW", "W", "NW")

lit.enrich <- with(asp.tbl, sum(pct.cand[octant %in% lit.oct]) /
                            sum(pct.background[octant %in% lit.oct]))
shd.enrich <- with(asp.tbl, sum(pct.cand[octant %in% shade.oct]) /
                            sum(pct.background[octant %in% shade.oct]))

say(glue("
    strongest octant             : {asp.worst} ({round(max(asp.tbl$enrichment), 2)}x)
    spread, strongest / weakest  : {round(asp.spread, 2)}x
    sun-facing (NE,E,SE)         : {round(lit.enrich, 2)}x
    shaded (SW,W,NW)             : {round(shd.enrich, 2)}x

    A spread under about 1.3x is no meaningful aspect preference, and
    illumination then has no fingerprint here. A strong preference for the
    sun-facing or the shaded group is the fingerprint; a strong preference for
    north or south alone is NOT illumination — at 3-4 degrees north the sun
    passes nearly overhead through the year, so north and south slopes are lit
    about equally and a north-south asymmetry points at something else, most
    likely rainfall or aspect-linked land use."))


# =============================================================================
#  PHASE 4 — ROUTE 1: HOW MUCH GREENING CAN ILLUMINATION ACTUALLY EXPLAIN?
# =============================================================================
#  Two estimates, deliberately built to bracket the answer.
#
#  DIFFERENCED (the honest one). Regress dNDVI on delta cos(i) with protected
#  area fixed effects. Differencing the same pixel across dates removes
#  everything time-invariant about it — its vegetation, its soil, its slope —
#  so the coefficient is identified off the CHANGE in illumination alone. This
#  is the same identification strategy as script 08, and for the same reason:
#  within-date slopes measure terrain and ecology, not artefact.
#
#  CROSS-SECTIONAL (the generous one). Regress 2019 NDVI on 2019 cos(i) across
#  pixels, and multiply that sensitivity by the observed mean delta cos(i).
#  This is an UPPER BOUND and badly confounded: shaded slopes are genuinely
#  moister and more wooded here, so the cross-sectional slope absorbs real
#  ecology and overstates the instrument's sensitivity to light. It is included
#  precisely because it is generous — if even the generous estimate is small,
#  the hypothesis is finished.
#
#  WHAT THIS TEST CAN AND CANNOT SEE — established by simulation BEFORE the
#  terrain raster was exported, so it is a prior expectation and not a story
#  told around whatever came out. The differenced estimator above was checked
#  over 40 synthetic replications with a coefficient of 0.30 planted in: it
#  recovered 0.3002 with no detectable bias, and its reported standard error
#  matched the spread of the estimates. The estimator is sound. What it is up
#  against is the ARITHMETIC OF SMALL NUMBERS:
#
#    At n = 400,000 and a delta cos(i) spread of about 0.012, the standard
#    error on the coefficient is roughly 0.007, so anything above about 0.013
#    is detectable. That is ample power. But the greening illumination can
#    explain is the coefficient TIMES the mean shift in cos(i), and for a
#    timing shift of a few days that mean shift is of order 0.004:
#
#       coefficient 0.1  ->  +0.0004 NDVI,  3% of a +0.0127 greening
#       coefficient 0.3  ->  +0.0012 NDVI,  9%
#       coefficient 0.5  ->  +0.0020 NDVI, 16%
#       coefficient 1.0  ->  +0.0040 NDVI, 32%
#
#    NDVI's cross-sectional spread across aspects in rugged terrain is rarely
#    more than about 0.1 per unit cos(i), so coefficients near 1.0 are not
#    physically available. The ceiling on what illumination can explain here is
#    therefore SINGLE-DIGIT PERCENT unless the two years were observed many
#    days apart — which Phase 2 measures directly.
#
#    For reference, the worst case over aspect that a given timing shift can do
#    to cos(i) on a slope, computed at this latitude in early March:
#
#       shift     flat    5 deg   10 deg   20 deg   30 deg
#        1 day   0.0013   0.0018   0.0024   0.0034   0.0043
#        3 days  0.0037   0.0054   0.0071   0.0102   0.0130
#        5 days  0.0061   0.0089   0.0117   0.0169   0.0217
#       10 days  0.0113   0.0170   0.0227   0.0334   0.0431
#       20 days  0.0187   0.0305   0.0420   0.0641   0.0843
#
#  So the honest reading of this phase is: a large coefficient would be a real
#  finding about the sensor, but only a large TIMING SHIFT could make it matter
#  for the change estimate. Phase 2 is where the hypothesis lives or dies.

fit.diff <- lm(dndvi ~ dcosi + factor(pa), data = samp)
b.diff   <- coef(fit.diff)[["dcosi"]]
illum.diff <- b.diff * mean(samp$dcosi)

fit.cross <- lm(ndvi.a ~ cosi.a + factor(pa), data = samp)
b.cross   <- coef(fit.cross)[["cosi.a"]]
illum.cross <- b.cross * mean(samp$dcosi)

mean.dndvi <- mean(samp$dndvi)

cat("\n[4] ROUTE 1 — ILLUMINATION AS AN EXPLANATION FOR THE GREENING\n")
say(glue("
    observed mean dNDVI                        : {sprintf('%+.5f', mean.dndvi)}

    DIFFERENCED  dNDVI ~ delta cos(i) + PA
      coefficient on delta cos(i)              : {sprintf('%+.4f', b.diff)} NDVI per unit cos(i)
      mean delta cos(i)                        : {sprintf('%+.5f', mean(samp$dcosi))}
      greening attributable to illumination    : {sprintf('%+.5f', illum.diff)}
      share of the observed greening            : {sprintf('%+.1f', 100 * illum.diff / mean.dndvi)}%

    CROSS-SECTIONAL (upper bound, confounded by ecology)
      NDVI sensitivity to cos(i)               : {sprintf('%+.4f', b.cross)}
      implied greening at the observed shift   : {sprintf('%+.5f', illum.cross)}
      share of the observed greening            : {sprintf('%+.1f', 100 * illum.cross / mean.dndvi)}%"))

illum.share <- illum.diff / mean.dndvi
illum.share.max <- illum.cross / mean.dndvi


# =============================================================================
#  PHASE 5 — ROUTE 2: IS THE OBSERVATION-DEPTH ARTEFACT A TERRAIN ARTEFACT?
# =============================================================================
#  Script 08 found 28% of the greening attributable to the difference in clear
#  looks between the two dates. It did not ask WHERE. If terrain shadow is
#  being masked as cloud shadow then the depth artefact should be stronger on
#  steep ground, and the interaction below will show it.
#
#  Slope enters as a main effect too. In a DIFFERENCED regression a
#  time-invariant variable like slope has no business predicting the change
#  unless something time-varying is attached to it — so a significant slope
#  main effect is itself evidence of a terrain-linked artefact, not of ecology.

fit.r2 <- lm(dndvi ~ obs.gap.px * slope.deg + factor(pa), data = samp)
co.r2  <- coef(fit.r2)
se.r2  <- summary(fit.r2)$coefficients

b.gap   <- co.r2[["obs.gap.px"]]
b.slope <- co.r2[["slope.deg"]]
b.inter <- co.r2[["obs.gap.px:slope.deg"]]
p.inter <- se.r2["obs.gap.px:slope.deg", "Pr(>|t|)"]

## The depth effect evaluated on flat and on steep ground.
gap.flat  <- b.gap + b.inter * median(samp$slope.deg[samp$slope.deg <= slope.steep])
gap.steep <- b.gap + b.inter * median(samp$slope.deg[samp$slope.deg >  slope.steep])

## And what the whole terrain block explains: predict with the terrain-linked
## terms zeroed, exactly as 08 predicted at a zero observation gap.
pred.flat.zero <- mean(predict(fit.r2, newdata = samp %>%
                                 mutate(obs.gap.px = 0, slope.deg = 0)))

cat("\n[5] ROUTE 2 — DOES THE DEPTH ARTEFACT LIVE ON STEEP GROUND?\n")
say(glue("
    dNDVI ~ observation gap * slope + PA

      observation gap, main effect   : {sprintf('%+.6f', b.gap)} per clear look
      slope, main effect             : {sprintf('%+.6f', b.slope)} per degree
      gap x slope interaction        : {sprintf('%+.7f', b.inter)}  (p = {format.pval(p.inter, digits = 2)})

      depth effect on flat ground    : {sprintf('%+.6f', gap.flat)} per look
      depth effect on steep ground   : {sprintf('%+.6f', gap.steep)} per look
      steep / flat                   : {round(gap.steep / gap.flat, 2)}x

      mean observed depth gap        : {round(mean(samp$obs.gap.px), 2)} looks
      greening at zero gap AND zero slope : {sprintf('%+.5f', pred.flat.zero)}
      of an observed                      : {sprintf('%+.5f', mean.dndvi)}"))

if (!is.na(p.inter) && p.inter < 0.01 && abs(gap.steep) > 1.5 * abs(gap.flat)) {
  say("
    The depth artefact IS concentrated on steep ground. That is Route 2
    confirmed: Fmask is most likely masking terrain shadow as cloud shadow, so
    shaded slopes lose looks and the surviving looks are biased. This is the
    script 08 confound with a terrain mechanism behind it, and the candidate
    patches on steep ground should be treated as the least reliable.")
} else {
  say("
    No material concentration of the depth artefact on steep ground. The
    script 08 observation-depth confound is real but it is not terrain-driven,
    so correcting for depth (as 08 already does) does not need a separate
    terrain correction on top.")
}


# =============================================================================
#  PHASE 6 — HOW MUCH OF THE CANDIDATE AREA IS ILLUMINATION-EXPLICABLE?
# =============================================================================
#  A per-pixel version of Phase 4: what would illumination alone predict for
#  this pixel's dNDVI, and is that enough to account for the whole of it?
#
#  "Explicable" is a generous test on purpose — it uses the CROSS-SECTIONAL
#  sensitivity, the upper bound, and asks only whether the predicted
#  illumination effect has the same sign and at least the same magnitude as the
#  observed change. A pixel that fails this test cannot be an illumination
#  artefact; one that passes it merely might be.

samp <- samp %>%
  mutate(dndvi.illum = b.cross * dcosi,
         illum.explicable = sign(dndvi.illum) == sign(dndvi) &
                            abs(dndvi.illum) >= abs(dndvi))

cand.samp <- samp %>% filter(candidate == 1)

if (nrow(cand.samp) < 500) {
  stop(glue("Only {nrow(cand.samp)} candidate pixels landed in the sample — too ",
            "few to say anything about them. Raise `n.samp` above {n.samp}, or ",
            "check that {basename(f.cand)} really carries the candidate mask."))
}

expl.share <- mean(cand.samp$illum.explicable)

cat("\n[6] CANDIDATE AREA THAT ILLUMINATION COULD ACCOUNT FOR\n")
say(glue("
    candidate pixels in the sample            : {format(nrow(cand.samp), big.mark = ',')}
    illumination-explicable (generous test)   : {round(100 * expl.share, 2)}%
    of the {round(cand.km2.total, 1)} km2 candidate area, that is   : {round(cand.km2.total * expl.share, 1)} km2
    leaving                                   : {round(cand.km2.total * (1 - expl.share), 1)} km2

    candidate median slope                    : {round(median(cand.samp$slope.deg), 2)} degrees
    candidate median delta cos(i)             : {sprintf('%+.5f', median(cand.samp$dcosi))}
    candidate share on slopes above {slope.steep} deg     : {round(100 * mean(cand.samp$slope.deg > slope.steep), 1)}%

    The km2 figures scale the sample share onto the candidate area measured off
    the raster. They assume the sample is representative of the candidate
    population, which it is by construction, but they are NOT a re-delineation
    — no pixel has been removed from any map by this script."))


# =============================================================================
#  PHASE 7 — VERDICT
# =============================================================================

cat("\n================================================================\n")
cat(" VERDICT — can terrain illumination explain the candidate change?\n")
cat("================================================================\n\n")

route1.dead <- abs(illum.share.max) < 0.10 && asp.spread < 1.3
route1.live <- abs(illum.share) > 0.15 || asp.spread > 1.8
route2.live <- !is.na(p.inter) && p.inter < 0.01 &&
               abs(gap.steep) > 1.5 * abs(gap.flat)

cat(sprintf("  Observation timing shifted by %+.1f days between %d and %d,\n",
            doy.shift, yr.a, yr.b))
cat(sprintf("  moving the mean solar zenith %+.2f degrees and cos(i) by %+.5f.\n\n",
            mean(samp$zen.b) - mean(samp$zen.a), mean(samp$dcosi)))

cat(sprintf("  ROUTE 1  illumination explains %+.1f%% of the greening (differenced),\n",
            100 * illum.share))
cat(sprintf("           at most %+.1f%% on the generous cross-sectional bound.\n",
            100 * illum.share.max))
cat(sprintf("           Aspect preference: %.2fx spread across octants.\n\n", asp.spread))
cat(sprintf("  ROUTE 2  depth artefact is %.2fx stronger on steep ground than flat\n",
            gap.steep / gap.flat))
cat(sprintf("           (interaction p = %s).\n\n", format.pval(p.inter, digits = 2)))

if (route1.dead && !route2.live) {
  cat("  TERRAIN ILLUMINATION IS NOT THE EXPLANATION.\n\n")
  cat("  The mechanism was never available: both composites come from the same\n")
  cat("  calendar window, the two years were observed at nearly the same point\n")
  cat("  in that window, so the sun sat in nearly the same place and cos(i)\n")
  cat("  barely moved. There is also no aspect fingerprint — the candidates do\n")
  cat("  not favour lit or shaded slopes. The candidates ARE on steeper ground\n")
  cat("  than average, but in a differenced regression that is a statement\n")
  cat("  about where the vegetation is, not about the light.\n\n")
  cat("  Two explanations are left and this pipeline cannot separate them:\n")
  cat("  genuine woody thickening or regrowth, and real land-use change.\n")
  cat("  Separating those needs high-resolution imagery for a sample of\n")
  cat("  patches, or a multi-year NDVI series to tell a trend from one wet\n")
  cat("  year against one dry one. The two-date design cannot do it.\n")
} else if (route1.live) {
  cat("  ILLUMINATION IS MATERIAL. Part of the candidate area is an artefact of\n")
  cat("  how the sun struck these slopes, not of the ground changing. Report the\n")
  cat("  Phase 6 residual area rather than the 277.5 km2 headline, and apply a\n")
  cat("  topographic correction (C-correction or the SCS+C method) to the\n")
  cat("  composites before any supervised classification — the Random Forest\n")
  cat("  would otherwise learn illumination as if it were land cover.\n")
} else {
  cat("  ILLUMINATION PROPER IS NOT THE EXPLANATION, BUT TERRAIN IS INVOLVED.\n\n")
  cat("  cos(i) barely moved between the dates, so the light itself cannot have\n")
  cat("  produced the change. What terrain does do here is bias the OBSERVATIONS:\n")
  cat("  the depth artefact already identified in script 08 is concentrated on\n")
  cat("  steep ground, which points at terrain shadow being masked as cloud\n")
  cat("  shadow. The 28% depth correction in 08 already absorbs this; what it\n")
  cat("  adds is that the correction is not uniform across the landscape, and\n")
  cat("  candidates on steep ground deserve less confidence than those on flat.\n")
}


# =============================================================================
#  PHASE 8 — FIGURES
# =============================================================================

## 8.1 Aspect rose — the diagnostic figure.
p.rose <- asp.tbl %>%
  ggplot(aes(x = octant, y = enrichment, fill = enrichment)) +
  geom_col(width = 1, colour = "white", linewidth = 0.3) +
  geom_hline(yintercept = 1, colour = "grey20", linewidth = 0.4) +
  coord_polar(start = -pi / 8) +
  scale_fill_gradient2(low = "#2166ac", mid = "grey92", high = "#b2182b",
                       midpoint = 1, guide = "none") +
  labs(title = "Candidate change by slope aspect",
       subtitle = glue("Enrichment relative to the landscape; 1.0 = no preference. ",
                       "Slopes above {slope.flat} degrees only."),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA))

## 8.2 dNDVI against the illumination shift, binned. If illumination drove the
## greening this would slope; a flat line is the null made visible.
## unique() on the breaks: on near-flat ground delta cos(i) varies only with
## latitude and observation day, so the lower twentieths can tie and cut()
## refuses duplicate breaks.
dcosi.brk <- unique(quantile(samp$dcosi, seq(0, 1, 0.05), na.rm = TRUE))

bin.dcosi <- samp %>%
  mutate(bin = cut(dcosi, breaks = dcosi.brk,
                   include.lowest = TRUE, labels = FALSE)) %>%
  group_by(bin) %>%
  summarise(mid.dcosi = median(dcosi), mean.dndvi = mean(dndvi), n = n(),
            .groups = "drop")

p.dcosi <- bin.dcosi %>%
  ggplot(aes(mid.dcosi, mean.dndvi)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_point(aes(size = n), colour = "#1b7837", alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = "#762a83", linewidth = 0.6) +
  scale_size_continuous(range = c(1.5, 4), guide = "none") +
  labs(title = "Greening against the change in illumination",
       subtitle = glue("Twentieths of the delta cos(i) distribution. ",
                       "Differenced slope {sprintf('%+.4f', b.diff)} NDVI per unit cos(i)."),
       x = "delta cos(i), 2024 minus 2019", y = "mean dNDVI") +
  theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

## 8.3 Candidate frequency over the slope-elevation plane.
p.plane <- samp %>%
  mutate(s.bin = cut(slope.deg, breaks = c(-Inf, 2, 5, 10, 15, 20, Inf),
                     labels = c("0-2", "2-5", "5-10", "10-15", "15-20", "20+")),
         e.bin = cut(elevation, breaks = 6, dig.lab = 4)) %>%
  group_by(s.bin, e.bin) %>%
  summarise(pct.cand = 100 * mean(candidate), n = n(), .groups = "drop") %>%
  filter(n >= 200) %>%
  ggplot(aes(s.bin, e.bin, fill = pct.cand)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_viridis_c(option = "magma", direction = -1,
                       name = "% flagged") +
  labs(title = "Where candidates sit in slope and elevation",
       subtitle = "Cells with fewer than 200 sampled pixels are left blank",
       x = "slope (degrees)", y = "elevation (m)") +
  theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave(file.path(path.out.terr, "terrain_aspect_rose.png"), p.rose,
       width = 5.5, height = 5.5, dpi = 200, bg = "white")
ggsave(file.path(path.out.terr, "terrain_dndvi_vs_dcosi.png"), p.dcosi,
       width = 7, height = 4.6, dpi = 200, bg = "white")
ggsave(file.path(path.out.terr, "terrain_slope_elevation_plane.png"), p.plane,
       width = 7, height = 4.6, dpi = 200, bg = "white")

## 8.4 The visual version of the aspect test: candidates over a hillshade, in
## whichever protected area carries the highest candidate fraction. If the
## candidates sat consistently on one side of every ridge, it would be visible
## here and nowhere in a table.
pa.rank <- samp %>%
  group_by(pa) %>%
  summarise(pct = 100 * mean(candidate), n = n(), .groups = "drop") %>%
  filter(n >= 2000) %>%
  arrange(desc(pct))

pa.pick <- if (nrow(pa.rank) > 0) as.character(pa.rank$pa[1]) else NA_character_
one.pa  <- if (!is.na(pa.pick)) aoi.parts.v[aoi.parts.v$pa.label == pa.pick, ] else NULL

## If rasterize() handed back integer codes rather than labels the match fails
## and one.pa is empty. Say so and skip the map rather than dying here — every
## number above is already computed and worth printing.
if (is.null(one.pa) || nrow(one.pa) == 0) {

  say(glue("
    Skipping the hillshade map: could not match the protected area '{pa.pick}'
    back to a polygon in the AOI. Every statistic above is unaffected."))

} else {

  crop.to.pa <- function(r) terra::mask(terra::crop(r, one.pa), one.pa)

  slope.r  <- crop.to.pa(terr.m[["slope"]])  / terr.scale[["slope"]]
  aspect.r <- crop.to.pa(terr.m[["aspect"]]) / terr.scale[["aspect"]]

  hs <- terra::shade(deg2rad(slope.r), deg2rad(aspect.r),
                     angle = 45, direction = 110)   # sun in the ESE, as at overpass
  names(hs) <- "hillshade"

  cand.pa <- crop.to.pa(cand)
  cand.pa <- terra::ifel(cand.pa == 1, 1, NA)

  p.map <- ggplot() +
    geom_spatraster(data = hs, maxcell = 6e5) +
    scale_fill_gradient(low = "grey15", high = "white", na.value = NA,
                        guide = "none") +
    ggnewscale::new_scale_fill() +
    geom_spatraster(data = as.factor(cand.pa), maxcell = 6e5) +
    scale_fill_manual(values = c("1" = "#d7191c"), na.value = NA,
                      labels = c("candidate change"), name = NULL,
                      na.translate = FALSE) +
    geom_sf(data = sf::st_as_sf(one.pa), fill = NA, colour = "#1b7837",
            linewidth = 0.5) +
    annotation_scale(location = "bl", width_hint = 0.25, height = unit(0.15, "cm")) +
    coord_sf(expand = FALSE) +
    labs(title = glue("Candidate change over terrain — {pa.pick}"),
         subtitle = glue("{round(pa.rank$pct[1], 1)}% of this area flagged. ",
                         "Hillshade lit from the east-southeast, as at the HLS overpass. ",
                         "Candidates clinging to one side of every ridge would mean illumination.")) +
    theme_void(base_size = 9) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 8),
          plot.background = element_rect(fill = "white", colour = NA))

  f.map <- file.path(path.out.terr, glue("terrain_candidates_over_hillshade.png"))
  ggsave(f.map, p.map, width = 7.5, height = 7.5, dpi = 220, bg = "white")
  say(glue("    wrote {f.map}"))
}


# =============================================================================
#  PHASE 9 — OUTPUTS
# =============================================================================

readr::write_csv(range.tbl,   file.path(path.out.terr, "terrain_band_ranges.csv"))
readr::write_csv(sun.cmp,     file.path(path.out.terr, "terrain_solar_geometry.csv"))
readr::write_csv(terrain.cmp, file.path(path.out.terr, "terrain_candidates_vs_background.csv"))
readr::write_csv(asp.tbl,     file.path(path.out.terr, "terrain_aspect_enrichment.csv"))
readr::write_csv(bin.dcosi,   file.path(path.out.terr, "terrain_dndvi_by_dcosi.csv"))

readr::write_csv(
  tibble(
    quantity = c("mean.dndvi", "doy.shift.days", "mean.dcosi",
                 "b.differenced", "illum.attributable", "illum.share",
                 "b.crosssectional", "illum.upper.bound", "illum.share.max",
                 "aspect.spread", "lit.enrichment", "shaded.enrichment",
                 "gap.effect.flat", "gap.effect.steep", "gap.interaction.p",
                 "cand.illum.explicable.share"),
    value = c(mean.dndvi, doy.shift, mean(samp$dcosi),
              b.diff, illum.diff, illum.share,
              b.cross, illum.cross, illum.share.max,
              asp.spread, lit.enrich, shd.enrich,
              gap.flat, gap.steep, p.inter,
              expl.share)
  ),
  file.path(path.out.terr, "terrain_headline_numbers.csv")
)

cat("\n================================================================\n")
cat(" DONE — tables and figures in", path.out.terr, "\n")
cat("================================================================\n\n")

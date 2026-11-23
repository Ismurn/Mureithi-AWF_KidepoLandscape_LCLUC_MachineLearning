#  KIDEPO LANDSCAPE - LCLUC
#  Shared constants and helpers for the R stage
#
#  Sourced by the numbered scripts; it defines nothing that runs on its own.
#      source("scripts/R/00_functions.R")
#
#  Author: Isaac Mureithi - IMureithi@awf.org


#  PHASE 0 - PACKAGES

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(sf, terra, tidyverse, glue, scico, viridis, patchwork)


#  PHASE 1 - PATHS

## Raw and GEE-derived data live on OneDrive, outside the repo.
path.data.root <- "C:/Users/IsaacM/OneDrive - African Wildlife Foundation/R WORK/Uganda's Kidepo Region"
path.gee.out   <- file.path(path.data.root, "kidepo_lcluc")

## Derived artefacts stay in the repo but gitignored.
path.out.aoi <- file.path("outputs", "aoi")
path.out.qc  <- file.path("outputs", "qc")
path.out.cls <- file.path("outputs", "classification")

for (p in c(path.out.qc, path.out.cls)) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}


#  PHASE 2 - CONSTANTS
#  These MUST agree with scripts/gee/02_build_composites.js. If the GEE encoding
#  changes, change it here in the same commit.

crs.geo <- 4326    # WGS84 geographic
crs.utm <- 32636   # UTM Zone 36N - the working CRS for everything local

comp.window <- "dry_late"
comp.years  <- c(2019, 2024)

comp.nodata    <- -9999   # sentinel written by the GEE export
comp.scale.int <- 10000   # reflectance and indices were multiplied by this
comp.min.obs   <- 4       # clear observations below which a pixel is weakly supported

## Band order as exported. Bands 1-9 are scaled; clear_obs is a raw count.
comp.bands.sr    <- c("BLUE", "GREEN", "RED", "NIR", "SWIR1", "SWIR2")
comp.bands.index <- c("NDVI", "NDWI", "NBR")
comp.bands.all   <- c(comp.bands.sr, comp.bands.index, "clear_obs")

## Local solar time of the HLS overpass, in decimal hours. Landsat 8/9 cross
## the descending node near 10:00 and Sentinel-2 near 10:30, and HLS merges
## both, so 10.25 is the midpoint. Used by the illumination helpers in 3.9.
##
## This is an APPROXIMATION and the equation of time is ignored, together worth
## up to about 4 degrees of hour angle. That matters less than it looks: the
## error is the same in both comparison years, so it very nearly cancels in any
## difference between them, which is the only quantity drawn on. Absolute
## cos(i) from these helpers is indicative; a difference in cos(i) is a
## measurement.
overpass.solar.hour <- 10.25


#  PHASE 3 - HELPERS

#' 3.1 Path to one exported composite.
composite.path <- function(year, window = comp.window, dir = path.gee.out) {
  file.path(dir, glue("kidepo_hls_{window}_{year}_utm36n_30m.tif"))
}

#' 3.2 Load one composite.
#'
#' The raster is left as Int16 ON PURPOSE. Dividing 9 bands of ~32 million cells
#' by 10000 would materialise about a gigabyte of float temp files for no
#' analytical gain: every band shares the same scale factor, so band ratios,
#' distances and cluster structure are identical either way. Rescale only when a
#' number has to be read by a human - use to.reflectance() for that.
#'
#' NAflag() is used rather than a reclassify: it tells GDAL to treat the sentinel
#' as NA on read, which costs nothing, where `r[r == -9999] <- NA` would rewrite
#' every cell.
load.composite <- function(year, window = comp.window, dir = path.gee.out) {

  f <- composite.path(year, window, dir)
  if (!file.exists(f)) {
    stop(glue("Composite not found:\n  {f}\n",
              "Export it from scripts/gee/02_build_composites.js and download it, ",
              "or point `dir` at wherever it was saved."))
  }

  r <- terra::rast(f)

  if (terra::nlyr(r) != length(comp.bands.all)) {
    stop(glue("{basename(f)} has {terra::nlyr(r)} bands; expected ",
              "{length(comp.bands.all)} ({paste(comp.bands.all, collapse = ', ')})."))
  }
  names(r) <- comp.bands.all

  ## Treat the export sentinel as missing data.
  terra::NAflag(r) <- comp.nodata

  r
}

#' 3.3 Convert scaled bands to true units (reflectance 0-1, indices -1..1).
#' Applied to a small object - a summary table or a cropped tile - not to a
#' full-size raster.
to.reflectance <- function(x) x / comp.scale.int

#' 3.4 Load both dates and refuse to continue unless they share a grid.
#'
#' CLAUDE.md caveat 3. The GEE script exports both dates with the same CRS,
#' scale and region, which SHOULD align them - but "should" is not a check, and a
#' silent half-pixel offset would quietly corrupt every change statistic
#' downstream.
load.composite.pair <- function(years = comp.years, window = comp.window,
                                dir = path.gee.out) {

  r.1 <- load.composite(years[1], window, dir)
  r.2 <- load.composite(years[2], window, dir)

  ok <- terra::compareGeom(r.1, r.2, stopOnError = FALSE, messages = FALSE)
  if (!isTRUE(ok)) {
    stop(glue("The {years[1]} and {years[2]} composites are NOT on the same grid. ",
              "Re-export both from GEE with an identical `region`, `crs` and `scale`."))
  }

  list(a = r.1, b = r.2, years = years, window = window)
}

#' 3.5 The AOI layers written by 01_build_aoi.R.
load.aoi <- function() {
  f.union <- file.path(path.out.aoi, "kidepo.pa.union.geojson")
  f.parts <- file.path(path.out.aoi, "kidepo.pa.individual.geojson")

  if (!file.exists(f.union) || !file.exists(f.parts)) {
    stop("AOI files missing - run scripts/R/01_build_aoi.R first.")
  }

  list(
    union = sf::st_read(f.union, quiet = TRUE) %>% sf::st_transform(crs.utm),
    parts = sf::st_read(f.parts, quiet = TRUE) %>% sf::st_transform(crs.utm)
  )
}

#' 3.6 Console helper - glue() eats one trailing newline.
say <- function(x) cat(x, "\n", sep = "")

#' 3.7 Remove patches below a minimum mapping unit.
#'
#' This replaces the patches() + freq() + ifel() idiom the earlier scripts used,
#' which was the slowest step in the whole pipeline by a wide margin: labelling
#' every connected component across a 32-million-cell raster that contains
#' 150,000 scattered patches. terra::sieve() does the same job in C++ and is
#' 10-50x faster — 0.2 s against 7-21 s per protected area when measured, and
#' minutes against seconds over the full AOI.
#'
#' sieve() ALONE IS NOT EQUIVALENT. It replaces every sub-threshold patch with
#' the value of its largest neighbour, zeros included, so it also fills HOLES
#' smaller than the threshold inside large patches and quietly inflates the
#' area. Masking the result back against the input keeps only cells that were
#' flagged to begin with. That combination reproduces the patches() result
#' exactly — verified cell-for-cell over Kidepo Valley, Karenga, Imatong and
#' Zulia, with zero cells differing in each.
#'
#' @param mask01 SpatRaster of 0/1, NA outside the area of interest
#' @param mmu.ha minimum patch size to keep, in hectares
apply.mmu <- function(mask01, mmu.ha = 1) {

  cell.m2 <- prod(terra::res(mask01))
  mmu.px  <- ceiling(mmu.ha * 10000 / cell.m2)

  sieved <- terra::sieve(mask01, threshold = mmu.px, directions = 8)

  terra::mask(terra::as.int(sieved == 1 & mask01 == 1), mask01)
}

#' 3.8 Is a binary change layer spatially clustered, or randomly scattered?
#'
#' Real land-cover change has edges - a cleared field, an expanding settlement.
#' Pixel noise does not. Two layers can report the same "% changed" while one is
#' a signal and the other is nothing.
#'
#' The test: among flagged pixels, count how many of the 8 neighbours are also
#' flagged, and compare with what spatial randomness would give. If a fraction p
#' of the area is flagged and flagging were random, a flagged pixel would average
#' 8p flagged neighbours. Materially more than that means real patches.
#'
#' A ratio near 1 means noise. It does NOT follow that a high ratio means land
#' COVER change: burn scars are coherent patches too. This separates signal from
#' noise, not fire from clearing.
#'
#' @param mask01 SpatRaster of 0/1 (NA outside the area of interest)
#' @return a one-row tibble
change.coherence <- function(mask01) {

  mask01 <- terra::as.int(mask01)

  nb.sum <- terra::focal(mask01, w = matrix(1, 3, 3), fun = "sum",
                         na.policy = "omit", na.rm = TRUE)

  ## Neighbours flagged, excluding the focal pixel, over flagged pixels only.
  nb.of.flagged <- terra::mask(nb.sum - mask01, mask01, maskvalues = c(0, NA))

  p.flagged <- terra::global(mask01, "mean", na.rm = TRUE)[1, 1]
  obs.nb    <- terra::global(nb.of.flagged, "mean", na.rm = TRUE)[1, 1]
  exp.nb    <- 8 * p.flagged
  pct.iso   <- 100 * terra::global(nb.of.flagged <= 1, "mean", na.rm = TRUE)[1, 1]

  tibble::tibble(
    pct.flagged     = 100 * p.flagged,
    expected.nb     = exp.nb,
    observed.nb     = obs.nb,
    coherence.ratio = obs.nb / exp.nb,
    pct.isolated    = pct.iso
  )
}


#  PHASE 3.9 - SOLAR GEOMETRY AND TERRAIN ILLUMINATION
#
#  Used by 10_terrain_test.R to ask whether the candidate change is an artefact
#  of how the sun struck these slopes. They live here rather than in that
#  script so that scripts/R/tests/test_solar_geometry.R can validate the SAME
#  code the analysis runs - a test holding its own copy of the arithmetic would
#  drift away from it silently and then certify the wrong thing.
#
#  All angles are in DEGREES at the interface. Radians appear only inside.

deg2rad <- function(d) d * pi / 180
rad2deg <- function(r) r * 180 / pi

#' 3.9.1 Solar declination in degrees, from day of year (Cooper's equation).
#' Accurate to a few tenths of a degree, which is well inside what this is for.
solar.declination <- function(doy) 23.44 * sin(deg2rad(360 * (doy + 284) / 365))

#' 3.9.2 Solar zenith and azimuth at a place and time.
#'
#' Azimuth is measured CLOCKWISE FROM NORTH, the same convention GEE's
#' ee.Terrain.aspect() uses, so the two can be differenced directly. Getting
#' that convention wrong would put the sun in the west and invert every
#' conclusion about which slopes are lit, so it is stated rather than assumed.
#'
#' @param doy       day of year, may be a vector
#' @param lat.deg   latitude in degrees north, may be a vector
#' @param solar.hour local solar time in decimal hours
#' @return list of `zenith` and `azimuth`, both in degrees
solar.position <- function(doy, lat.deg, solar.hour = overpass.solar.hour) {

  dec <- deg2rad(solar.declination(doy))
  phi <- deg2rad(lat.deg)
  h   <- deg2rad(15 * (solar.hour - 12))   # negative before local noon

  cos.z <- pmin(pmax(sin(phi) * sin(dec) + cos(phi) * cos(dec) * cos(h), -1), 1)
  z     <- acos(cos.z)

  sin.a <- -cos(dec) * sin(h) / sin(z)
  cos.a <- (sin(dec) - sin(phi) * cos.z) / (cos(phi) * sin(z))

  list(zenith  = rad2deg(z),
       azimuth = (rad2deg(atan2(sin.a, cos.a)) + 360) %% 360)
}

#' 3.9.3 Cosine of the local solar incidence angle on an inclined surface.
#'
#' 1 means the sun is normal to the slope, 0 means grazing, below 0 means the
#' slope is turned away from the sun and lit only by sky and neighbours. On
#' flat ground it reduces to cos(zenith), which is the cheapest available check
#' that the formula has been typed correctly.
cos.illum <- function(slope.deg, aspect.deg, zenith.deg, azimuth.deg) {
  s <- deg2rad(slope.deg);  o <- deg2rad(aspect.deg)
  z <- deg2rad(zenith.deg); a <- deg2rad(azimuth.deg)
  cos(z) * cos(s) + sin(z) * sin(s) * cos(a - o)
}


#  PHASE 3.10 - LINEAR TREND IN CLOSED FORM
#
#  Used by 11_ndvi_timeseries.R, which needs the slope of an annual series both
#  for ~200,000 sampled pixels and for every cell of an 9-band raster. Fitting
#  lm() that many times is not an option, and neither is trusting hand-expanded
#  algebra that nothing checks - so the closed forms live here and
#  scripts/R/tests/test_trend_estimators.R asserts them equal to lm().

#' 3.10.1 OLS slope weights for a FIXED set of x values.
#'
#' The key fact this exploits: for a fixed x, the OLS slope is a weighted SUM of
#' the y values, with weights that do not depend on y at all. So the slope of a
#' raster stack is one linear combination of its bands - a handful of passes
#' rather than a regression per cell.
#'
#' It also means the slope operator is LINEAR, which is what lets a depth or
#' fire adjustment be applied to the slopes directly:
#'     slope(y - b*z)  ==  slope(y) - b*slope(z)
#' That identity is relied on in script 11 Phase 6 and is tested rather than
#' assumed.
ols.slope.weights <- function(x) {
  xd <- x - mean(x)
  xd / sum(xd^2)
}

#' 3.10.2 Slope of y on x. Equal to coef(lm(y ~ x))[2].
ols.slope <- function(x, y) sum(ols.slope.weights(x) * y)

#' 3.10.3 R-squared of that straight line. Equal to summary(lm(y ~ x))$r.squared.
#'
#' This is the second half of the question script 11 asks of a pixel. A steep
#' slope says a pixel ended higher than it started; a high R2 says it got there
#' STEADILY. Vegetation thickening gives both. A pixel swinging with the rains
#' that happens to be caught low in the first year and high in the last gives
#' the first without the second, and telling those apart is the whole point of
#' having a series rather than two dates.
#'
#' Returns NA for a flat series, where no line explains anything and the ratio
#' would be 0/0.
ols.r2 <- function(x, y) {
  xd  <- x - mean(x)
  sst <- sum((y - mean(y))^2)
  if (sst <= 0) return(NA_real_)
  b <- sum(xd * y) / sum(xd^2)
  min(1, (b^2 * sum(xd^2)) / sst)
}

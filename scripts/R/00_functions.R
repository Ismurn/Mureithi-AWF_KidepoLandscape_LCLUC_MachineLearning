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

#' 3.7 Is a binary change layer spatially clustered, or randomly scattered?
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

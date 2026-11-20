# =============================================================================
#  KIDEPO LANDSCAPE — LCLUC
#  Script 02 — Quality control of the two GEE composites
#
#  PURPOSE
#  Interrogate the composites BEFORE classifying anything. A classification will
#  happily run on a broken stack and return a confident, wrong answer, so the
#  point here is to find the problems while they are still cheap:
#
#    1. Do the two rasters share a grid, bands and encoding?   (hard stop if not)
#    2. Where is each date weakly sampled, and where do the two disagree about
#       how well sampled they are? Those pixels are where "change" is least
#       trustworthy.
#    3. Are the reflectance values physically plausible?
#    4. How much of each date is burnt, by NBR? The GEE burn profile said ~20% of
#       the AOI carries scars in both years; this checks that at 30 m, in the
#       actual composites, rather than taking MODIS's 500 m word for it.
#    5. Do the two dates differ in a way that looks like land cover, or in a way
#       that looks like fire and illumination?
#
#  HOW TO RUN — from the repo root, in PowerShell (NOT Git Bash: MSYS DLLs on
#  PATH shadow GDAL and R segfaults):
#      & "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" scripts\R\02_qc_composites.R
#
#  OUTPUTS -> outputs/qc/   (gitignored)
#
#  Author: Isaac Mureithi - IMureithi@awf.org
# =============================================================================

source("scripts/R/00_functions.R")

terra::terraOptions(progress = 0)

cat("\n================================================================\n")
cat(" KIDEPO LCLUC — 02 QC of composites\n")
cat(" run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n")


# =============================================================================
#  PHASE 1 — LOAD AND VERIFY
# =============================================================================

comp <- load.composite.pair()
r.a  <- comp$a      # earlier date
r.b  <- comp$b      # later date
yr.a <- comp$years[1]
yr.b <- comp$years[2]

aoi <- load.aoi()
aoi.union.v <- terra::vect(aoi$union)
aoi.parts.v <- terra::vect(aoi$parts)

say(glue("
[1] GEOMETRY  (compareGeom passed, or we would not be here)
    window        : {comp$window}
    dates         : {yr.a} and {yr.b}
    CRS           : {terra::crs(r.a, describe = TRUE)$name}
    resolution    : {paste(terra::res(r.a), collapse = ' x ')} m
    dimensions    : {terra::nrow(r.a)} rows x {terra::ncol(r.a)} cols x {terra::nlyr(r.a)} bands
    cells per band: {format(terra::ncell(r.a), big.mark = ',')}
    bands         : {paste(names(r.a), collapse = ', ')}
"))


# =============================================================================
#  PHASE 2 — CLEAR-OBSERVATION SUPPORT
# =============================================================================
#  The honest confidence map. A median composite built from 2 looks is a
#  different object from one built from 15, and the change between two such
#  pixels is mostly noise.

obs.a <- r.a[["clear_obs"]]
obs.b <- r.b[["clear_obs"]]

## Restrict every statistic to the reporting extent — the PA union — rather than
## the 500 m buffered export, so the numbers describe what will be reported on.
obs.a.m <- terra::mask(terra::crop(obs.a, aoi.union.v), aoi.union.v)
obs.b.m <- terra::mask(terra::crop(obs.b, aoi.union.v), aoi.union.v)

obs.stats <- function(x, label) {
  v <- terra::global(x, fun = c("mean", "min", "max"), na.rm = TRUE)
  frac.thin <- terra::global(x < comp.min.obs, fun = "mean", na.rm = TRUE)[1, 1]
  frac.zero <- terra::global(x == 0,           fun = "mean", na.rm = TRUE)[1, 1]
  tibble(
    date          = label,
    mean.obs      = round(v$mean, 2),
    min.obs       = v$min,
    max.obs       = v$max,
    pct.below.min = round(100 * frac.thin, 2),
    pct.zero.obs  = round(100 * frac.zero, 3)
  )
}

obs.table <- bind_rows(obs.stats(obs.a.m, as.character(yr.a)),
                       obs.stats(obs.b.m, as.character(yr.b)))

cat("\n[2] CLEAR-OBSERVATION SUPPORT inside the PA union\n")
cat("    pct.below.min = % of PA area with fewer than", comp.min.obs, "clear looks\n\n")
print(as.data.frame(obs.table), row.names = FALSE)

## Pixels thin in EITHER date — the union of weakness, which is what actually
## limits a change map.
thin.either <- (obs.a.m < comp.min.obs) | (obs.b.m < comp.min.obs)
pct.thin.either <- 100 * terra::global(thin.either, "mean", na.rm = TRUE)[1, 1]

say(glue("\n    Thin in EITHER date: {round(pct.thin.either, 2)}% of the PA area.",
         "\n    This is the ceiling on how much of the change map is well supported."))


# =============================================================================
#  PHASE 3 — ARE THE VALUES PHYSICALLY PLAUSIBLE?
# =============================================================================
#  Reflectance outside 0-1 and indices outside -1..1 mean the scale factor or the
#  nodata handling is wrong. Cheap to check, expensive to discover later.

band.stats <- function(r, label) {
  s <- terra::global(r[[c(comp.bands.sr, comp.bands.index)]],
                     fun = c("mean", "min", "max"), na.rm = TRUE)
  tibble(
    date = label,
    band = rownames(s),
    mean = round(to.reflectance(s$mean), 4),
    min  = round(to.reflectance(s$min),  4),
    max  = round(to.reflectance(s$max),  4)
  )
}

cat("\n[3] BAND VALUES (rescaled to true units; reflectance 0-1, indices -1..1)\n\n")
band.table <- bind_rows(band.stats(r.a, as.character(yr.a)),
                        band.stats(r.b, as.character(yr.b)))
print(as.data.frame(band.table), row.names = FALSE)

## Problems worth stopping for are collected and reported together at the end,
## so one run surfaces everything rather than failing at the first thing found.
qc.problems <- character(0)

## 3.1 Outside the physical range — scale factor or nodata handling is wrong.
impossible <- band.table %>%
  filter((band %in% comp.bands.sr    & (min < -0.01 | max > 1.6)) |
         (band %in% comp.bands.index & (min < -1.01 | max > 1.01)))

if (nrow(impossible) > 0) {
  cat("\n    !! VALUES OUTSIDE THE PHYSICAL RANGE — check the scale factor / nodata:\n")
  print(as.data.frame(impossible), row.names = FALSE)
  qc.problems <- c(qc.problems,
                   glue("{nrow(impossible)} band(s) outside the physical range"))
}

## 3.2 DEGENERATE BANDS — a band with no variation carries no information.
##
## This check exists because the first export passed 3.1 and was still broken:
## every reflectance band was uniformly ZERO, which is inside 0-1 and so looked
## "physically plausible". The cause was the 0.0001 HLS scale factor being
## applied twice; the indices hid it, because a normalised difference is
## scale-invariant and came out perfectly correct. A range test cannot catch
## that. A constancy test can.
degenerate <- band.table %>%
  filter(min == max) %>%
  mutate(issue = "constant — no variation")

## A reflectance band whose maximum rounds to zero is the specific signature of
## the double-scaling bug, worth naming explicitly.
collapsed <- band.table %>%
  filter(band %in% comp.bands.sr, abs(max) < 1e-4) %>%
  mutate(issue = "reflectance collapsed to zero — scale factor applied twice?")

degenerate.all <- bind_rows(degenerate, collapsed) %>% distinct(date, band, .keep_all = TRUE)

if (nrow(degenerate.all) > 0) {
  cat("\n    !! DEGENERATE BAND(S) — the stack is NOT usable for classification:\n")
  degenerate.all %>%
    transmute(date, band, mean, min, max, issue) %>%
    as.data.frame() %>%
    print(row.names = FALSE)
  qc.problems <- c(qc.problems,
                   glue("{nrow(degenerate.all)} degenerate band(s) carrying no information"))
} else if (nrow(impossible) == 0) {
  cat("\n    OK: all bands lie inside their physical range and all vary.\n")
}


# =============================================================================
#  PHASE 4 — BURN EXTENT AT 30 m
# =============================================================================
#  The GEE burn profile used MODIS at 500 m. Here the composites are asked
#  directly, via NBR. The threshold is a convention, not a truth: NBR < 0.1 is a
#  commonly used "burnt or very sparse" cut in savanna. It is reported as a
#  CONSISTENT YARDSTICK across the two dates, not as a burnt-area product.

nbr.a <- terra::mask(terra::crop(r.a[["NBR"]], aoi.union.v), aoi.union.v)
nbr.b <- terra::mask(terra::crop(r.b[["NBR"]], aoi.union.v), aoi.union.v)

nbr.cut <- 0.1 * comp.scale.int

pct.low.a <- 100 * terra::global(nbr.a < nbr.cut, "mean", na.rm = TRUE)[1, 1]
pct.low.b <- 100 * terra::global(nbr.b < nbr.cut, "mean", na.rm = TRUE)[1, 1]

## dNBR: positive means the later date is LESS burnt / greener than the earlier.
dnbr <- nbr.b - nbr.a
dnbr.mean <- terra::global(dnbr, "mean", na.rm = TRUE)[1, 1] / comp.scale.int

say(glue("
[4] BURN / BARE SIGNAL  (NBR < 0.1 inside the PA union)
    {yr.a}: {round(pct.low.a, 2)}% of PA area
    {yr.b}: {round(pct.low.b, 2)}% of PA area
    difference  : {round(pct.low.b - pct.low.a, 2)} percentage points
    mean dNBR   : {round(dnbr.mean, 4)}  ({yr.b} minus {yr.a}; positive = greener/less burnt)

    MODIS put cumulative burnt area at 19.7% ({yr.a}) and 22.8% ({yr.b}).
    A broadly similar gap here means the two dates are matched on fire about as
    well as the record allows. A much larger gap means the fire signal is
    stronger at 30 m than MODIS suggested, and a burnt class becomes essential
    rather than merely advisable.
"))


# =============================================================================
#  PHASE 5 — QC FIGURES
# =============================================================================
#  terra's base plotting is used deliberately: it downsamples a 32-million-cell
#  raster on the fly, where ggplot/tidyterra would try to build a data frame of
#  it. Publication maps come later, from the classified output, where tidyterra
#  is the better tool.

qc.png <- function(filename, width = 1800, height = 1100, res = 150) {
  grDevices::png(file.path(path.out.qc, filename),
                 width = width, height = height, res = res)
}

## 5.1 False colour, both dates. SWIR2/NIR/GREEN: burnt ground reads dark
## magenta-brown, vegetation green, water black.
qc.png("qc_false_colour_pair.png")
par(mfrow = c(1, 2), mar = c(2, 2, 3, 1))
terra::plotRGB(r.a, r = "SWIR2", g = "NIR", b = "GREEN",
               stretch = "lin", mar = c(2, 2, 3, 1))
title(main = glue("{yr.a} — SWIR2/NIR/GREEN"))
terra::plot(aoi.union.v, add = TRUE, border = "white", lwd = 1)
terra::plotRGB(r.b, r = "SWIR2", g = "NIR", b = "GREEN",
               stretch = "lin", mar = c(2, 2, 3, 1))
title(main = glue("{yr.b} — SWIR2/NIR/GREEN"))
terra::plot(aoi.union.v, add = TRUE, border = "white", lwd = 1)
dev.off()

## 5.2 Clear observations, both dates, on a shared scale so they are comparable.
obs.max <- max(terra::global(obs.a.m, "max", na.rm = TRUE)[1, 1],
               terra::global(obs.b.m, "max", na.rm = TRUE)[1, 1])

qc.png("qc_clear_obs_pair.png")
par(mfrow = c(1, 2), mar = c(2, 2, 3, 4))
terra::plot(obs.a.m, range = c(0, obs.max), col = viridis::viridis(50),
            main = glue("{yr.a} — clear observations"))
terra::plot(aoi.union.v, add = TRUE, border = "grey20", lwd = 0.6)
terra::plot(obs.b.m, range = c(0, obs.max), col = viridis::viridis(50),
            main = glue("{yr.b} — clear observations"))
terra::plot(aoi.union.v, add = TRUE, border = "grey20", lwd = 0.6)
dev.off()

## 5.3 Where the change map will be weak.
qc.png("qc_thin_support.png", width = 1200, height = 1100)
par(mar = c(2, 2, 3, 4))
terra::plot(thin.either, col = c("grey92", "#d7301f"),
            main = glue("Fewer than {comp.min.obs} clear looks in {yr.a} or {yr.b}"),
            plg = list(legend = c("well sampled", "thin")))
terra::plot(aoi.parts.v, add = TRUE, border = "grey20", lwd = 0.6)
dev.off()

## 5.4 dNBR — the first look at where the two dates actually differ.
qc.png("qc_dnbr.png", width = 1200, height = 1100)
par(mar = c(2, 2, 3, 4))
terra::plot(dnbr / comp.scale.int, range = c(-0.4, 0.4),
            col = scico::scico(50, palette = "vik", direction = -1),
            main = glue("dNBR: {yr.b} minus {yr.a}  (blue = greener/less burnt)"))
terra::plot(aoi.parts.v, add = TRUE, border = "grey10", lwd = 0.6)
dev.off()

say(glue("\n[5] wrote 4 QC figures to {path.out.qc}/"))


# =============================================================================
#  PHASE 6 — QC TABLES FOR THE WRITE-UP
# =============================================================================

readr::write_csv(obs.table,  file.path(path.out.qc, "qc_clear_obs_summary.csv"))
readr::write_csv(band.table, file.path(path.out.qc, "qc_band_summary.csv"))

summary.tbl <- tibble(
  metric = c("window", "date_a", "date_b", "pct_thin_either",
             "pct_nbr_low_a", "pct_nbr_low_b", "mean_dnbr"),
  value  = c(comp$window, yr.a, yr.b, round(pct.thin.either, 3),
             round(pct.low.a, 3), round(pct.low.b, 3), round(dnbr.mean, 5))
)
readr::write_csv(summary.tbl, file.path(path.out.qc, "qc_summary.csv"))

say(glue("[6] wrote 3 QC tables to {path.out.qc}/"))


# =============================================================================
#  PHASE 7 — VERDICT
# =============================================================================
#  The figures and tables above are written either way: a broken stack is easier
#  to diagnose with its pictures in front of you. But the script must not end
#  quietly as though everything passed, because the next step would then build a
#  classification on data that cannot support one.

cat("\n================================================================\n")

if (length(qc.problems) > 0) {
  cat(" QC FAILED\n")
  for (p in qc.problems) cat("   -", p, "\n")
  cat("\n Figures and tables were still written to", path.out.qc, "for diagnosis.\n")
  cat(" Fix the cause in scripts/gee/02_build_composites.js, re-export, re-run.\n")
  cat("================================================================\n\n")
  stop("QC failed — see the problems listed above.", call. = FALSE)
}

cat(" QC PASSED. Next: scripts/R/03_kmeans_firstlook.R\n")
cat("================================================================\n\n")

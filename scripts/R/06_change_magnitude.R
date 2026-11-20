# =============================================================================
#  KIDEPO LANDSCAPE — LCLUC
#  Script 06 — How big is the change, in spectral terms?
#
#  WHY THIS EXISTS
#  Scripts 03 to 05 measured change by counting how many clusters a pixel moved
#  along the NDVI ordering. That metric has a flaw which only became visible
#  once the cluster centroids were in front of us:
#
#      cluster   1      2      3      4      5      6
#      NDVI      0.252  0.259  0.307  0.387  0.521  0.792
#
#  The gaps are wildly uneven. A 1 -> 3 move counts as TWO steps but spans only
#  0.055 NDVI, while 5 -> 6 counts as ONE and spans 0.271. So "two or more
#  steps" does not mean "moved a lot"; for the crowded low clusters it can mean
#  almost nothing. The 1,600 km2 of "substantive change" and the 877 km2
#  residual both inherit that flaw.
#
#  This script re-measures the same change in the units the sensor actually
#  reports, and asks how much of it is large enough to be land-cover conversion
#  rather than a landscape that is simply a little greener than it was.
#
#  HOW TO RUN — from the repo root, in PowerShell:
#      & "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" scripts\R\06_change_magnitude.R
#
#  Requires 03 and 05 to have been run.
#  OUTPUTS -> outputs/classification/
#
#  Author: Isaac Mureithi - IMureithi@awf.org
# =============================================================================

source("scripts/R/00_functions.R")

terra::terraOptions(progress = 0)
set.seed(42)

n.sample <- 400000
mmu.ha   <- 1.0

## Thresholds on |dNDVI| at which a change starts to look like conversion
## rather than a greenness anomaly. Reported as a sensitivity range, because
## there is no single defensible cut.
mag.thresholds <- c(0.10, 0.15, 0.20)

cat("\n================================================================\n")
cat(" KIDEPO LCLUC — 06 Spectral magnitude of the change\n")
cat(" run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n")


# =============================================================================
#  PHASE 1 — LOAD
# =============================================================================

comp <- load.composite.pair()
yr.a <- comp$years[1]
yr.b <- comp$years[2]

aoi <- load.aoi()
aoi.union.v <- terra::vect(aoi$union)
aoi.parts.v <- terra::vect(aoi$parts)
mask.to.aoi <- function(r) terra::mask(terra::crop(r, aoi.union.v), aoi.union.v)

ndvi.a <- mask.to.aoi(comp$a[["NDVI"]])
ndvi.b <- mask.to.aoi(comp$b[["NDVI"]])
nbr.a  <- mask.to.aoi(comp$a[["NBR"]])
nbr.b  <- mask.to.aoi(comp$b[["NBR"]])

steps <- terra::rast(file.path(path.out.cls, "kmeans_change_steps.tif"))
resid <- terra::rast(file.path(path.out.cls, "change_residual_modis_mmu.tif"))

d.ndvi <- ndvi.b - ndvi.a      # still in stored units (x10000)
d.nbr  <- nbr.b  - nbr.a

px.km2 <- prod(terra::res(steps)) / 1e6


# =============================================================================
#  PHASE 2 — SPECTRAL MOVEMENT BY CLUSTER STEPS
# =============================================================================

samp <- terra::spatSample(c(d.ndvi, d.nbr, steps, resid), n.sample,
                          method = "random", na.rm = TRUE, warn = FALSE)
names(samp) <- c("dndvi", "dnbr", "steps", "resid")

samp <- as_tibble(samp) %>%
  mutate(dndvi = to.reflectance(dndvi),
         dnbr  = to.reflectance(dnbr))

by.steps <- samp %>%
  group_by(steps) %>%
  summarise(n              = n(),
            pct.of.area    = 100 * n() / nrow(samp),
            median.dNDVI   = median(dndvi),
            median.abs     = median(abs(dndvi)),
            q90.abs        = quantile(abs(dndvi), 0.90),
            .groups = "drop")

cat("\n[2] SPECTRAL MOVEMENT BY NUMBER OF CLUSTER STEPS\n\n")
by.steps %>% mutate(across(where(is.numeric), ~round(.x, 4))) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\n    Note how little 2 and 3 steps actually move. The step count is a rank,\n")
cat("    not a distance, and the three low clusters sit within 0.055 NDVI of\n")
cat("    each other - so a 2-step move among them is a very small spectral change.\n")


# =============================================================================
#  PHASE 3 — THE RESIDUAL AGAINST THE BACKGROUND
# =============================================================================
#  Two questions. Do residual pixels move further than everything else? And do
#  they move in BOTH directions? Land-cover conversion produces gains and
#  losses; a climate or rainfall difference pushes the whole landscape one way.

landscape.shift <- median(samp$dndvi, na.rm = TRUE)

by.resid <- samp %>%
  group_by(resid) %>%
  summarise(n            = n(),
            median.dNDVI = median(dndvi),
            median.abs   = median(abs(dndvi)),
            q90.abs      = quantile(abs(dndvi), 0.90),
            pct.greening = 100 * mean(dndvi > 0),
            .groups = "drop")

cat("\n[3] RESIDUAL (fire removed, size filtered) vs BACKGROUND\n\n")
by.resid %>% mutate(across(where(is.numeric), ~round(.x, 4))) %>%
  as.data.frame() %>% print(row.names = FALSE)

say(glue("
    landscape-wide median dNDVI : {round(landscape.shift, 4)}  ({yr.b} minus {yr.a})

    If the residual were land-cover conversion it would contain clearing as well
    as regrowth, so pct.greening would sit near 50. A value far above that says
    the residual is the landscape moving one way together - which is what a
    wetter or later-senescing season does, and what no burn mask can remove."))


# =============================================================================
#  PHASE 4 — HOW MUCH IS BIG ENOUGH TO BE CONVERSION?
# =============================================================================
#  Computed over the full rasters, not the sample, because these are the numbers
#  that would be quoted.

mmu.pixels <- ceiling(mmu.ha * 10000 / prod(terra::res(steps)))

magnitude.area <- map_dfr(mag.thresholds, function(thr) {

  big <- terra::as.int(abs(d.ndvi) > thr * comp.scale.int & resid == 1)

  ## Size filter, so scattered pixels over the threshold are not counted.
  ## apply.mmu() (00_functions.R) uses terra::sieve() — identical result to the
  ## patches() idiom this used to run, at a fraction of the cost.
  big.mmu <- terra::subst(apply.mmu(big, mmu.ha = mmu.ha), from = NA, to = 0)

  km2 <- terra::global(big.mmu, "sum", na.rm = TRUE)[1, 1] * px.km2

  ## Patch count, now that patch IDs are no longer produced as a side effect.
  n.patch <- terra::freq(terra::patches(terra::subst(big.mmu, from = 0, to = NA),
                                        directions = 8, zeroAsNA = TRUE)) %>%
    nrow()

  tibble(threshold.dNDVI = thr, km2 = km2, n.patches = n.patch)
})

area.total <- terra::global(!is.na(steps), "sum", na.rm = TRUE)[1, 1] * px.km2

magnitude.area <- magnitude.area %>%
  mutate(pct.of.area = 100 * km2 / area.total)

cat("\n[4] CANDIDATE CONVERSION: residual AND a large spectral move AND >= 1 ha\n\n")
magnitude.area %>% mutate(across(where(is.numeric), ~round(.x, 2))) %>%
  as.data.frame() %>% print(row.names = FALSE)

say(glue("
    For comparison:
      whole PA union analysed        : {format(round(area.total, 1), big.mark = ',')} km2
      'substantive change' (2+ steps): 1,600 km2
      residual after fire + MMU      : 877 km2

    The drop from those headline figures to the table above is the point of this
    script. Most of what the step metric called change is a small, one-directional
    greenness shift, not conversion."))


# =============================================================================
#  PHASE 5 — OUTPUTS
# =============================================================================

readr::write_csv(by.steps,       file.path(path.out.cls, "magnitude_by_steps.csv"))
readr::write_csv(by.resid,       file.path(path.out.cls, "magnitude_residual_vs_background.csv"))
readr::write_csv(magnitude.area, file.path(path.out.cls, "magnitude_candidate_conversion.csv"))

best <- terra::as.int(abs(d.ndvi) > 0.15 * comp.scale.int & resid == 1)
terra::writeRaster(best, file.path(path.out.cls, "change_candidate_strong.tif"),
                   overwrite = TRUE, datatype = "INT1U")

png(file.path(path.out.cls, "change_magnitude.png"), width = 1900, height = 800, res = 150)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
hist(samp$dndvi[samp$resid == 0], breaks = 80, col = "grey80", border = NA,
     xlim = c(-0.4, 0.4), freq = FALSE,
     main = "dNDVI: background vs residual", xlab = glue("NDVI {yr.b} - {yr.a}"))
hist(samp$dndvi[samp$resid == 1], breaks = 80, col = adjustcolor("#1b7837", 0.55),
     border = NA, freq = FALSE, add = TRUE)
abline(v = 0, lty = 2); abline(v = landscape.shift, col = "#d7301f", lty = 2)
legend("topright", bty = "n", cex = 0.8,
       fill = c("grey80", adjustcolor("#1b7837", 0.55)),
       legend = c("background", "residual"))
par(mar = c(2, 2, 3, 4))
terra::plot(best, col = c("grey92", "#762a83"), type = "classes",
            levels = c("no", "strong candidate"),
            main = "Residual with |dNDVI| > 0.15")
terra::plot(aoi.parts.v, add = TRUE, border = "grey20", lwd = 0.5)
dev.off()

say(glue("\n[5] wrote 3 tables, 1 raster and 1 figure to {path.out.cls}/"))

cat("\n================================================================\n")
cat(" DONE.\n")
cat("================================================================\n\n")

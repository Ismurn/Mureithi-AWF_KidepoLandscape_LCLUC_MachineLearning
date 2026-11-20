# =============================================================================
#  KIDEPO LANDSCAPE — LCLUC
#  Script 05 — Redo the fire attribution against recorded MODIS fires
#
#  PURPOSE
#  Script 04 had no burnt-area raster and inferred one from an NBR quantile.
#  The inference failed its own self-check: burn switching predicted change only
#  1.24x better than chance. This script replaces the guess with the MODIS
#  record exported by 01c, and asks three questions in order:
#
#    1. Was the 04 proxy wrong, and how wrong?
#    2. Does SCAR AGE explain the failure? The hypothesis was that a December
#       scar has greened past the lowest NBR quantile by February, so the
#       threshold flags permanently dark ground instead of actual fires.
#    3. With real fires, does fire explain the change or not?
#
#  Question 3 is the one that matters, and the answer is allowed to be "no".
#
#  A NOTE ON WHICH BANDS ARE USED
#  The 01c export has three bands, but only `days_since_burn` is trustworthy:
#  `burnt_season` and `burnt_in_window` came back as zero everywhere, while
#  days_since_burn carries the correct pattern (-1 unburnt, 80-130 days for
#  burnt pixels, totalling 19.7% of the AOI exactly as MODIS reported). Since
#  the age band strictly contains the binary ones, both are DERIVED from it here
#  rather than re-exported. 01c has been fixed to derive them the same way.
#
#  Outside the AOI the export wrote 0 rather than the -9999 sentinel, so every
#  layer is masked to the PA union before use and 0 is never read as "burnt
#  today".
#
#  HOW TO RUN — from the repo root, in PowerShell:
#      & "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" scripts\R\05_fire_attribution_modis.R
#
#  Requires 03 and 04 to have been run, and the 01c rasters downloaded.
#  OUTPUTS -> outputs/classification/
#
#  Author: Isaac Mureithi - IMureithi@awf.org
# =============================================================================

source("scripts/R/00_functions.R")

terra::terraOptions(progress = 0)
set.seed(42)


# =============================================================================
#  PHASE 0 — PARAMETERS
# =============================================================================

k.chosen  <- 6
min.steps <- 2
mmu.ha    <- 1.0

## The dry_late window runs 1 Feb to 31 Mar = 59 days. A scar with
## days_since_burn <= 58 therefore burnt DURING the window; anything older
## burnt before it opened.
window.days <- 59

n.sample <- 300000   # pixels for the age / NBR relationship

cat("\n================================================================\n")
cat(" KIDEPO LCLUC — 05 Fire attribution against MODIS\n")
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

burn.path <- function(year) {
  file.path(path.gee.out,
            glue("kidepo_burn_{comp.window}_{year}_utm36n_30m.tif"))
}

for (y in comp.years) {
  if (!file.exists(burn.path(y))) {
    stop(glue("Burn raster missing for {y}:\n  {burn.path(y)}\n",
              "Run scripts/gee/01c_export_burn_masks.js and download the output."))
  }
}

b.a <- terra::rast(burn.path(yr.a))
b.b <- terra::rast(burn.path(yr.b))

## Hard stop rather than a silent resample.
if (!isTRUE(terra::compareGeom(b.a, comp$a, stopOnError = FALSE, messages = FALSE))) {
  stop("Burn rasters are not on the composite grid — re-export with matching crs/scale/region.")
}

days.a <- mask.to.aoi(b.a[["days_since_burn"]])
days.b <- mask.to.aoi(b.b[["days_since_burn"]])

## DERIVED, not read: see the header. > 0 rather than >= 0, because 0 is what
## the export wrote outside the AOI and must never be read as a same-day burn.
burn.a <- terra::as.int(days.a > 0)
burn.b <- terra::as.int(days.b > 0)

px.km2 <- prod(terra::res(days.a)) / 1e6

say(glue("
[1] MODIS BURNED AREA, derived from days_since_burn
    {yr.a}: {round(100 * terra::global(burn.a, 'mean', na.rm = TRUE)[1,1], 2)}% of the PA union burnt this season
    {yr.b}: {round(100 * terra::global(burn.b, 'mean', na.rm = TRUE)[1,1], 2)}%
    (01b reported 19.7% and 22.8% from the same source)
"))


# =============================================================================
#  PHASE 2 — WAS THE 04 PROXY WRONG?
# =============================================================================
#  fire_state.tif encodes the NBR-quantile proxy: 0 neither, 1 = A only,
#  2 = B only, 3 = both. Compare it per date against the MODIS record.

f.fire.state <- file.path(path.out.cls, "fire_state.tif")
if (!file.exists(f.fire.state)) stop("fire_state.tif missing — run 04_separate_fire.R first.")

fire.state <- terra::rast(f.fire.state)
proxy.a <- terra::as.int(fire.state == 1 | fire.state == 3)   # burn-like in A
proxy.b <- terra::as.int(fire.state == 2 | fire.state == 3)   # burn-like in B

agreement <- function(proxy, truth, label) {

  tp <- terra::global(proxy == 1 & truth == 1, "sum", na.rm = TRUE)[1, 1]
  fp <- terra::global(proxy == 1 & truth == 0, "sum", na.rm = TRUE)[1, 1]
  fn <- terra::global(proxy == 0 & truth == 1, "sum", na.rm = TRUE)[1, 1]
  tn <- terra::global(proxy == 0 & truth == 0, "sum", na.rm = TRUE)[1, 1]

  n <- tp + fp + fn + tn
  po <- (tp + tn) / n
  ## Cohen's kappa: agreement corrected for what chance alone would produce.
  ## Raw agreement flatters a rare class badly — calling nothing burnt would
  ## score 80% here.
  pe <- ((tp + fp) * (tp + fn) + (fn + tn) * (fp + tn)) / n^2

  tibble(
    date      = label,
    precision = tp / (tp + fp),   # of pixels the proxy called burnt, how many were
    recall    = tp / (tp + fn),   # of pixels MODIS called burnt, how many were found
    agreement = po,
    kappa     = (po - pe) / (1 - pe)
  )
}

val <- bind_rows(agreement(proxy.a, burn.a, as.character(yr.a)),
                 agreement(proxy.b, burn.b, as.character(yr.b)))

cat("\n[2] THE 04 NBR PROXY vs THE MODIS RECORD\n\n")
val %>% mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\n    Kappa near 0 means the proxy was no better than chance at finding\n")
cat("    the fires, which would explain why filtering on it achieved nothing.\n")


# =============================================================================
#  PHASE 3 — DOES SCAR AGE EXPLAIN THE FAILURE?
# =============================================================================
#  The hypothesis from 04: by February a December scar has recovered enough to
#  leave the lowest NBR quantile, so the proxy cannot see it. If true, NBR
#  should climb back towards unburnt values as days_since_burn grows, and the
#  proxy's recall should collapse for older scars.

nbr.a <- mask.to.aoi(comp$a[["NBR"]])

samp <- terra::spatSample(c(days.a, nbr.a, proxy.a),
                          n.sample, method = "random", na.rm = TRUE, warn = FALSE)
names(samp) <- c("days", "nbr", "proxy")

samp <- samp %>%
  as_tibble() %>%
  filter(days != 0) %>%                       # drop the outside-AOI sentinel
  mutate(nbr = to.reflectance(nbr),
         burnt = days > 0)

nbr.unburnt <- median(samp$nbr[!samp$burnt], na.rm = TRUE)

age.profile <- samp %>%
  filter(burnt) %>%
  mutate(age.bin = cut(days, breaks = seq(0, 160, by = 20), include.lowest = TRUE)) %>%
  group_by(age.bin) %>%
  summarise(
    n            = n(),
    median.nbr   = median(nbr, na.rm = TRUE),
    pct.found    = 100 * mean(proxy == 1, na.rm = TRUE),
    .groups      = "drop"
  ) %>%
  mutate(nbr.vs.unburnt = median.nbr - nbr.unburnt)

cat(glue("\n[3] SCAR AGE vs SPECTRAL SIGNAL  (median NBR of UNBURNT ground = {round(nbr.unburnt, 4)})\n\n"))
age.profile %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\n    nbr.vs.unburnt is how far below unburnt ground a scar of that age sits.\n")
cat("    pct.found is the share of those scars the 04 proxy flagged.\n")

png(file.path(path.out.cls, "burn_age_vs_nbr.png"), width = 1500, height = 700, res = 150)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(age.profile$age.bin %>% as.integer() * 20 - 10, age.profile$median.nbr,
     type = "b", pch = 19, xlab = "days since burn", ylab = "median NBR",
     main = "Scar signal vs age")
abline(h = nbr.unburnt, col = "#1b7837", lty = 2)
text(x = 40, y = nbr.unburnt, labels = "unburnt ground", pos = 3, cex = 0.7, col = "#1b7837")
plot(age.profile$age.bin %>% as.integer() * 20 - 10, age.profile$pct.found,
     type = "b", pch = 19, xlab = "days since burn", ylab = "% flagged by the 04 proxy",
     main = "Proxy detection vs age", ylim = c(0, 100))
dev.off()


# =============================================================================
#  PHASE 4 — ATTRIBUTION WITH REAL FIRES
# =============================================================================

f.steps <- file.path(path.out.cls, "kmeans_change_steps.tif")
if (!file.exists(f.steps)) stop("kmeans_change_steps.tif missing — run 03 first.")
steps <- terra::rast(f.steps)

substantive <- terra::as.int(steps >= min.steps)

## Same logic as 04, now on recorded fires: burnt in both dates explains
## nothing, because fire is then a constant for that pixel.
switched <- terra::as.int(burn.a != burn.b)

area.total   <- terra::global(!is.na(substantive), "sum", na.rm = TRUE)[1, 1] * px.km2
sub.km2      <- terra::global(substantive, "sum", na.rm = TRUE)[1, 1] * px.km2
switched.km2 <- terra::global(switched,    "sum", na.rm = TRUE)[1, 1] * px.km2
fire.km2     <- terra::global(substantive == 1 & switched == 1, "sum", na.rm = TRUE)[1, 1] * px.km2

p.sub            <- sub.km2 / area.total
expected.overlap <- switched.km2 * p.sub
enrichment       <- fire.km2 / expected.overlap
risk.switched    <- fire.km2 / switched.km2
risk.not         <- (sub.km2 - fire.km2) / (area.total - switched.km2)
risk.ratio       <- risk.switched / risk.not

say(glue("
[4] ATTRIBUTION AGAINST RECORDED FIRES
    analysed area          : {format(round(area.total, 1), big.mark = ',')} km2
    substantive change     : {format(round(sub.km2, 1), big.mark = ',')} km2
    MODIS burn switched    : {format(round(switched.km2, 1), big.mark = ',')} km2
    both                   : {format(round(fire.km2, 1), big.mark = ',')} km2

    overlap if independent : {round(expected.overlap, 1)} km2
    enrichment             : {round(enrichment, 2)}x   (1.0 = no association)
    change rate where burn switched : {round(100 * risk.switched, 1)}%
    change rate where it did not    : {round(100 * risk.not, 1)}%
    risk ratio             : {round(risk.ratio, 2)}x
"))

if (enrichment < 1.5) {
  cat("    VERDICT: even against RECORDED fires, burning barely predicts where\n")
  cat("    the clusters changed. Fire is therefore not the main driver of the\n")
  cat("    change signal, and the earlier failure was not only the proxy's\n")
  cat("    fault. The change has some other cause - most plausibly differing\n")
  cat("    rainfall and greenness between the two dry seasons - which a burn\n")
  cat("    mask cannot remove and a supervised classifier must handle.\n")
} else {
  cat("    VERDICT: recorded fire does predict change, so the 04 failure was the\n")
  cat("    proxy's fault. Removing it is now a meaningful correction.\n")
}


# =============================================================================
#  PHASE 5 — RESIDUAL, SIZE FILTER, COHERENCE
# =============================================================================

residual <- terra::as.int(substantive == 1 & switched == 0)

mmu.pixels <- ceiling(mmu.ha * 10000 / prod(terra::res(steps)))
res.masked <- terra::subst(residual, from = 0, to = NA)
patch.id   <- terra::patches(res.masked, directions = 8, zeroAsNA = TRUE)

keep <- terra::freq(patch.id) %>% as_tibble() %>% filter(count >= mmu.pixels)

residual.mmu <- terra::ifel(patch.id %in% keep$value, 1, 0)
residual.mmu <- terra::subst(residual.mmu, from = NA, to = 0)
residual.mmu <- terra::mask(residual.mmu, substantive)

res.mmu.km2 <- terra::global(residual.mmu, "sum", na.rm = TRUE)[1, 1] * px.km2

coh.res <- change.coherence(residual.mmu)

say(glue("
[5] RESIDUAL AFTER REMOVING RECORDED FIRE
    before size filter : {format(round(terra::global(residual, 'sum', na.rm = TRUE)[1,1] * px.km2, 1), big.mark = ',')} km2
    after {mmu.ha} ha MMU   : {format(round(res.mmu.km2, 1), big.mark = ',')} km2
                         {round(100 * res.mmu.km2 / area.total, 2)}% of the analysed area
    coherence ratio    : {round(coh.res$coherence.ratio, 2)}
"))


# =============================================================================
#  PHASE 6 — OUTPUTS
# =============================================================================

terra::writeRaster(terra::as.int(burn.a * 1 + burn.b * 2),
                   file.path(path.out.cls, "modis_fire_state.tif"),
                   overwrite = TRUE, datatype = "INT1U")
terra::writeRaster(residual.mmu,
                   file.path(path.out.cls, "change_residual_modis_mmu.tif"),
                   overwrite = TRUE, datatype = "INT1U")

readr::write_csv(val,         file.path(path.out.cls, "modis_proxy_validation.csv"))
readr::write_csv(age.profile, file.path(path.out.cls, "burn_age_profile.csv"))
readr::write_csv(
  tibble(metric = c("analysed_km2", "substantive_km2", "switched_km2", "overlap_km2",
                    "expected_overlap_km2", "enrichment", "risk_ratio",
                    "residual_mmu_km2", "residual_coherence"),
         value  = round(c(area.total, sub.km2, switched.km2, fire.km2,
                          expected.overlap, enrichment, risk.ratio,
                          res.mmu.km2, coh.res$coherence.ratio), 3)),
  file.path(path.out.cls, "modis_attribution_summary.csv"))

png(file.path(path.out.cls, "modis_fire_vs_residual.png"), width = 1900, height = 1000, res = 150)
par(mfrow = c(1, 2), mar = c(2, 2, 3, 6))
terra::plot(terra::as.int(burn.a * 1 + burn.b * 2),
            col = c("grey92", "#f1a340", "#998ec3", "#4d4d4d"), type = "classes",
            levels = c("neither", glue("{yr.a} only"), glue("{yr.b} only"), "both"),
            main = "MODIS burnt, by date")
terra::plot(aoi.parts.v, add = TRUE, border = "grey20", lwd = 0.5)
terra::plot(residual.mmu, col = c("grey92", "#1b7837"), type = "classes",
            levels = c("no", "candidate"),
            main = glue("Change not explained by recorded fire (>= {mmu.ha} ha)"))
terra::plot(aoi.parts.v, add = TRUE, border = "grey20", lwd = 0.5)
dev.off()

cat("\n================================================================\n")
cat(" DONE.\n")
cat("================================================================\n\n")

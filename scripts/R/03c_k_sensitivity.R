#  KIDEPO LANDSCAPE - LCLUC
#  Script 03c - How much of the change was an artefact of choosing k=6?
#
#  THE QUESTION
#  03b showed that the landscape supports about three robust spectral strata,
#  not six, and that at k=6 a quarter of pixels sit ambiguously between
#  clusters. Those boundary pixels can flip between dates for no physical
#  reason. So: how much of the reported change survives if the clustering is
#  done at the k the data actually supports?
#
#  HOW THE COMPARISON IS MADE FAIR
#  "Cluster steps" cannot be compared across k - two steps out of six is not
#  two steps out of three. So the comparison rests on two k-independent
#  measures instead:
#
#    1. Did the pixel change cluster AT ALL, with fire removed and a size
#       filter applied? This still depends on k, and is expected to.
#    2. Did it change cluster AND move a large distance in NDVI? The magnitude
#       criterion is defined in sensor units and knows nothing about k. If the
#       method is sound, THIS number should be broadly stable across k. If it
#       swings wildly, then the estimate is an artefact of the clustering and
#       not a property of the landscape.
#
#  Measure 2 is the real test, and it is the one that would go in a write-up.
#
#  HOW TO RUN - from the repo root, in PowerShell:
#      & "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" scripts\R\03c_k_sensitivity.R
#
#  Requires the composites and the 01c burn rasters.
#  OUTPUTS -> outputs/classification/
#
#  Author: Isaac Mureithi - IMureithi@awf.org

source("scripts/R/00_functions.R")

terra::terraOptions(progress = 0)
set.seed(42)   # same draw as 03, so the k=6 column reproduces it

k.values      <- c(3, 6)
n.per.date    <- 60000
kmeans.nstart <- 25
mmu.ha        <- 1.0
mag.threshold <- 0.15      # |dNDVI| above which a change looks like conversion
feat.bands    <- c(comp.bands.sr, comp.bands.index)

cat("\n================================================================\n")
cat(" KIDEPO LCLUC - 03c sensitivity of the change signal to k\n")
cat(" run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n")


#  PHASE 1 - DATA

comp <- load.composite.pair()
yr.a <- comp$years[1]
yr.b <- comp$years[2]

aoi <- load.aoi()
aoi.union.v <- terra::vect(aoi$union)
mask.to.aoi <- function(r) terra::mask(terra::crop(r, aoi.union.v), aoi.union.v)

r.a <- mask.to.aoi(comp$a)
r.b <- mask.to.aoi(comp$b)

well.supported <- (r.a[["clear_obs"]] >= comp.min.obs) &
                  (r.b[["clear_obs"]] >= comp.min.obs)

r.a.f <- terra::mask(r.a[[feat.bands]], well.supported, maskvalues = c(FALSE, NA))
r.b.f <- terra::mask(r.b[[feat.bands]], well.supported, maskvalues = c(FALSE, NA))

d.ndvi <- r.b[["NDVI"]] - r.a[["NDVI"]]
px.km2 <- prod(terra::res(r.a)) / 1e6

## MODIS fire, derived from the age band (see 05 for why the binary bands are
## not used). > 0 rather than >= 0, because 0 is the outside-AOI fill.
burn.path <- function(year) {
  file.path(path.gee.out, glue("kidepo_burn_{comp.window}_{year}_utm36n_30m.tif"))
}
days.a <- mask.to.aoi(terra::rast(burn.path(yr.a))[["days_since_burn"]])
days.b <- mask.to.aoi(terra::rast(burn.path(yr.b))[["days_since_burn"]])
switched <- terra::as.int((days.a > 0) != (days.b > 0))

area.total <- terra::global(!is.na(well.supported) & well.supported,
                            "sum", na.rm = TRUE)[1, 1] * px.km2

say(glue("[1] analysed area: {format(round(area.total, 1), big.mark = ',')} km2 (both dates well supported)"))


#  PHASE 2 - ONE POOLED, STANDARDISED SAMPLE, SHARED BY BOTH k
#  Both fits see exactly the same data, so any difference is down to k alone.

samp <- bind_rows(
  as_tibble(terra::spatSample(r.a.f, n.per.date, "random", na.rm = TRUE, warn = FALSE)),
  as_tibble(terra::spatSample(r.b.f, n.per.date, "random", na.rm = TRUE, warn = FALSE))
) %>% drop_na()

X <- as.matrix(samp[, feat.bands])
feat.mu <- colMeans(X)
feat.sd <- apply(X, 2, sd)
X.z <- scale(X, center = feat.mu, scale = feat.sd)

say(glue("[2] pooled sample: {format(nrow(X.z), big.mark = ',')} pixels"))


#  PHASE 3 - CLASSIFY AND MEASURE CHANGE AT EACH k

mmu.pixels <- ceiling(mmu.ha * 10000 / prod(terra::res(r.a)))

## Assignment by raster algebra: one squared-distance surface per cluster, then
## which.min. Keeps the work in terra's C++ instead of an R callback per pixel.
classify <- function(r, centers.z) {
  r.z <- (r[[feat.bands]] - feat.mu) / feat.sd
  d <- lapply(seq_len(nrow(centers.z)), function(j) {
    sum((r.z - as.numeric(centers.z[j, ]))^2)
  })
  terra::which.min(terra::rast(d))
}

## The size filter now lives in 00_functions.R and uses terra::sieve(). The
## local patches()-based copy that used to sit here was the reason this script
## had to be run in the background; it is gone. Wrapped only to keep the
## subst(NA -> 0) that the rest of this script expects.
apply.mmu.local <- function(mask01) {
  terra::subst(apply.mmu(mask01, mmu.ha = mmu.ha), from = NA, to = 0)
}

results <- list()

for (k in k.values) {

  warned <- character(0)
  km <- withCallingHandlers(
    kmeans(X.z, centers = k, nstart = kmeans.nstart, iter.max = 100),
    warning = function(w) { warned <<- c(warned, conditionMessage(w))
                            invokeRestart("muffleWarning") })

  ## Order clusters by NDVI so cluster 1 is always the least vegetated.
  centers.true <- sweep(sweep(km$centers, 2, feat.sd, "*"), 2, feat.mu, "+")
  ord <- order(centers.true[, "NDVI"])
  centers.z <- km$centers[ord, , drop = FALSE]

  cls.a <- classify(r.a.f, centers.z)
  cls.b <- classify(r.b.f, centers.z)

  changed  <- terra::as.int(cls.a != cls.b)
  residual <- terra::as.int(changed == 1 & switched == 0)
  strong   <- terra::as.int(residual == 1 &
                            abs(d.ndvi) > mag.threshold * comp.scale.int)
  strong.mmu <- apply.mmu.local(strong)

  g <- function(x) terra::global(x, "sum", na.rm = TRUE)[1, 1] * px.km2

  results[[as.character(k)]] <- tibble(
    k                 = k,
    converged         = length(warned) == 0,
    ndvi.centroids    = paste(round(sort(to.reflectance(centers.true[, "NDVI"])), 3),
                              collapse = " / "),
    changed.km2       = g(changed),
    residual.km2      = g(residual),
    strong.km2        = g(strong),
    strong.mmu.km2    = g(strong.mmu)
  )

  terra::writeRaster(strong.mmu,
                     file.path(path.out.cls, glue("change_strong_k{k}.tif")),
                     overwrite = TRUE, datatype = "INT1U")

  say(glue("    k = {k}: classified, changed {round(100 * g(changed) / area.total, 1)}% of the area"))
}

res <- bind_rows(results) %>%
  mutate(changed.pct    = 100 * changed.km2 / area.total,
         residual.pct   = 100 * residual.km2 / area.total,
         strong.mmu.pct = 100 * strong.mmu.km2 / area.total)


#  PHASE 4 - THE COMPARISON

cat("\n[3] CHANGE AT EACH k\n\n")
res %>%
  transmute(k, converged,
            `changed km2`      = round(changed.km2, 1),
            `changed %`        = round(changed.pct, 1),
            `not fire km2`     = round(residual.km2, 1),
            `+ |dNDVI|>0.15`   = round(strong.km2, 1),
            `+ 1ha MMU km2`    = round(strong.mmu.km2, 1),
            `final %`          = round(strong.mmu.pct, 2)) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\n    NDVI centroids:\n")
for (i in seq_len(nrow(res))) {
  cat(sprintf("      k=%d: %s\n", res$k[i], res$ndvi.centroids[i]))
}

k3 <- res %>% filter(k == 3)
k6 <- res %>% filter(k == 6)

ratio.raw   <- k6$changed.km2    / k3$changed.km2
ratio.final <- k6$strong.mmu.km2 / k3$strong.mmu.km2

say(glue("
[4] VERDICT

    Raw 'changed cluster'   : {round(k3$changed.km2)} km2 at k=3 against {round(k6$changed.km2)} km2 at k=6
                              -> k=6 reports {round(ratio.raw, 2)}x as much
    After fire, magnitude
    and size filtering      : {round(k3$strong.mmu.km2, 1)} km2 at k=3 against {round(k6$strong.mmu.km2, 1)} km2 at k=6
                              -> {round(ratio.final, 2)}x
"))

if (abs(ratio.final - 1) < 0.35) {
  cat("    The filtered estimate is broadly STABLE across k. The raw cluster-change\n")
  cat("    figure is an artefact of how finely the continuum was sliced, but the\n")
  cat("    magnitude-filtered estimate is a property of the imagery and survives\n")
  cat("    the choice. That is the number to carry forward.\n")
} else {
  cat("    The filtered estimate moves substantially with k, so it is NOT yet a\n")
  cat("    property of the landscape. Change detection here should not depend on\n")
  cat("    clustering at all - threshold the spectral difference directly, or wait\n")
  cat("    for the supervised classifier.\n")
}

readr::write_csv(res, file.path(path.out.cls, "k_sensitivity.csv"))


#  PHASE 5 - THE ESTIMATE WITH NO CLUSTERING AT ALL
#  If the answer depends on k, stop using k. Nothing about detecting a large
#  spectral change requires clusters: the criteria are a big move in NDVI, no
#  fire switch to explain it, and a patch large enough to be credible. All
#  three are defined in sensor units and are reproducible by anyone.
#
#  This is the number that should go in a write-up, with the clustering kept
#  for what it is genuinely good at - suggesting what the classes ARE.
#
#  A note on comparing it with 06: that script started from "moved 2+ cluster
#  steps", whereas Phase 3 above starts from "changed cluster at all". Different
#  base sets, so their k=6 figures differ (123 km2 against 226 km2) without
#  either being wrong. The estimate below depends on neither.

strong.free <- terra::as.int(abs(d.ndvi) > mag.threshold * comp.scale.int &
                             switched == 0 &
                             well.supported)
strong.free.mmu <- apply.mmu.local(strong.free)

## Restrict to the analysed area before any RATE is computed.
##
## apply.mmu() ends with subst(NA -> 0), which fills the whole export rectangle
## with zeros - and that rectangle is about 3.7x the PA union. Areas are
## unaffected, because they are a sum of ones. But change.coherence() divides by
## the flagged FRACTION, so leaving the fill in place deflates p and inflates the
## coherence ratio absurdly: it reported 76.6 before this mask was added, where
## the true value is a few units.
strong.free.mmu <- terra::mask(strong.free.mmu, well.supported,
                               maskvalues = c(FALSE, NA))

free.km2 <- terra::global(strong.free.mmu, "sum", na.rm = TRUE)[1, 1] * px.km2

terra::writeRaster(strong.free.mmu,
                   file.path(path.out.cls, "change_strong_clusterfree.tif"),
                   overwrite = TRUE, datatype = "INT1U")

coh.free <- change.coherence(strong.free.mmu)

say(glue("
[5] CLUSTER-FREE ESTIMATE
    |dNDVI| > {mag.threshold}, no fire switch, >= {mmu.ha} ha patches
      area            : {format(round(free.km2, 1), big.mark = ',')} km2
      share of AOI    : {round(100 * free.km2 / area.total, 2)}%
      coherence ratio : {round(coh.free$coherence.ratio, 2)}

    For comparison, the cluster-dependent versions of the same idea gave
    {round(k3$strong.mmu.km2, 1)} km2 at k=3 and {round(k6$strong.mmu.km2, 1)} km2 at k=6. This figure replaces both.

    It remains a CANDIDATE set: a large greenness change with no fire to
    explain it is not yet land-cover conversion, and most of this landscape
    greened between the two dates for reasons that are probably climatic.
    Check these patches against high-resolution imagery before describing any
    of them as land-use change."))

readr::write_csv(
  tibble(metric = c("cluster_free_km2", "pct_of_aoi", "coherence_ratio",
                    "k3_km2", "k6_km2"),
         value  = round(c(free.km2, 100 * free.km2 / area.total,
                          coh.free$coherence.ratio,
                          k3$strong.mmu.km2, k6$strong.mmu.km2), 3)),
  file.path(path.out.cls, "change_clusterfree_summary.csv"))

cat("\n================================================================\n")
cat(" DONE.\n")
cat("================================================================\n\n")

sessionInfo()
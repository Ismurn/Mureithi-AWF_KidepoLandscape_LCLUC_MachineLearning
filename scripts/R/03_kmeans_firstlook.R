#  KIDEPO LANDSCAPE - LCLUC
#  Script 03 - Unsupervised (k-means) first look and change matrix
#
#  PURPOSE
#  A fast, training-data-free read on the class structure of the landscape and
#  on how much of it changed between the two dates. Following the ARSET
#  template: cluster both dates, cross-tabulate, derive a binary change mask,
#  and map the pair side by side.
#
#  This is a RECONNAISSANCE step, not the deliverable. k-means finds spectral
#  clusters, not land-cover classes, and it has no idea what a "cluster" means.
#  Its job here is to answer: how many separable classes does this landscape
#  actually support, which of them are burn scars, and is the change signal
#  large enough to be worth a supervised classification?
#
#  ONE DESIGN DECISION WORTH KNOWING ABOUT
#  k-means is fitted ONCE, on a pooled sample drawn from BOTH dates, and the
#  resulting centroids are then applied to each date. Clustering the dates
#  independently would be a mistake that is easy to miss: cluster 3 in 2019
#  would have no relationship to cluster 3 in 2024, so the cross-tabulation
#  would be counting coincidences. One shared feature space makes the change
#  matrix mean what it appears to mean.
#
#  HOW TO RUN - from the repo root, in PowerShell:
#      & "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" scripts\R\03_kmeans_firstlook.R
#
#  OUTPUTS -> outputs/classification/   (gitignored)
#
#  Author: Isaac Mureithi - IMureithi@awf.org

source("scripts/R/00_functions.R")

terra::terraOptions(progress = 0)
set.seed(42)   # reproducible sampling and cluster initialisation


#  PHASE 0 - PARAMETERS

k.candidates <- 3:9      # range scanned for the elbow plot

## k.chosen was originally set to 6 BEFORE the scan ran and never revisited,
## which 03b_kmeans_diagnostics.R has now corrected the record on. What that
## script found, so the choice is at least an informed one:
##
##     k   converged   % var   mean silhouette   % pixels ambiguous
##     3      yes       71.1        0.523              16.1
##     4      yes       78.7        0.490              18.9
##     5      NO        81.7        0.442              24.1
##     6      yes       84.1        0.440              24.8
##     7+     NO        86+         0.43 and falling   26+
##
## Silhouette falls monotonically as k rises, so the structure the data actually
## supports is about THREE strata, not six. At k=6 a QUARTER of pixels sit
## ambiguously between clusters - and those are precisely the pixels that flip
## between dates for no physical reason, which is a large part of why the raw
## change figure came out at 44.6%.
##
## k=6 is kept here deliberately, as reconnaissance: it resolves burnt ground,
## bare soil and the woodland gradient separately, which k=3 merges, and that
## detail informs the class scheme for the supervised stage. But NOTHING
## quantitative should be built on these six clusters. Use 06_change_magnitude.R
## for change in spectral units, and the supervised classifier for classes.
k.chosen     <- 6
n.sample     <- 60000    # pixels sampled PER DATE for fitting
kmeans.nstart <- 25      # restarts; k-means is sensitive to initialisation

## Features used for clustering. Reflectance plus the three indices.
feat.bands <- c(comp.bands.sr, comp.bands.index)

cat("\n================================================================\n")
cat(" KIDEPO LCLUC - 03 k-means first look\n")
cat(" run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n")


#  PHASE 1 - LOAD, MASK, AND RESTRICT TO WELL-SUPPORTED PIXELS

comp <- load.composite.pair()
yr.a <- comp$years[1]
yr.b <- comp$years[2]

aoi <- load.aoi()
aoi.union.v <- terra::vect(aoi$union)
aoi.parts.v <- terra::vect(aoi$parts)

## Mask to the PA union. Beyond restricting the analysis to the reporting
## extent, this removes the zero-fill that Earth Engine writes outside the clip
## geometry - zeros which would otherwise form a large, meaningless cluster of
## their own and distort every centroid.
mask.to.aoi <- function(r) terra::mask(terra::crop(r, aoi.union.v), aoi.union.v)
r.a <- mask.to.aoi(comp$a)
r.b <- mask.to.aoi(comp$b)

## Only cluster pixels that BOTH dates can support. A pixel thinly sampled in
## either date cannot produce a trustworthy change statement, so it is excluded
## from the fit and from the change matrix rather than quietly carried along.
well.supported <- (r.a[["clear_obs"]] >= comp.min.obs) &
                  (r.b[["clear_obs"]] >= comp.min.obs)

pct.supported <- 100 * terra::global(well.supported, "mean", na.rm = TRUE)[1, 1]

r.a.f <- terra::mask(r.a[[feat.bands]], well.supported, maskvalues = c(FALSE, NA))
r.b.f <- terra::mask(r.b[[feat.bands]], well.supported, maskvalues = c(FALSE, NA))

say(glue("
[1] DATA
    dates             : {yr.a} and {yr.b} ({comp$window})
    features          : {paste(feat.bands, collapse = ', ')}
    well supported    : {round(pct.supported, 2)}% of the PA union
                        (>= {comp.min.obs} clear observations in BOTH dates)
"))


#  PHASE 2 - SAMPLE AND STANDARDISE

samp.a <- terra::spatSample(r.a.f, n.sample, method = "random",
                            na.rm = TRUE, warn = FALSE)
samp.b <- terra::spatSample(r.b.f, n.sample, method = "random",
                            na.rm = TRUE, warn = FALSE)

samp.pooled <- bind_rows(
  as_tibble(samp.a) %>% mutate(date = yr.a),
  as_tibble(samp.b) %>% mutate(date = yr.b)
) %>%
  drop_na()

X <- as.matrix(samp.pooled[, feat.bands])

## STANDARDISATION IS NOT OPTIONAL HERE.
## k-means minimises Euclidean distance, so a feature with a wider numeric range
## dominates the result. Reflectance is stored on roughly 0-5000 while the
## indices span -10000 to 10000, so unstandardised clustering would be driven
## almost entirely by NDVI, NDWI and NBR, and the reflectance bands would barely
## participate. Centre and scale on the POOLED sample so both dates land in one
## shared feature space.
feat.mu <- colMeans(X)
feat.sd <- apply(X, 2, sd)

X.z <- scale(X, center = feat.mu, scale = feat.sd)

say(glue("[2] SAMPLE: {format(nrow(X), big.mark = ',')} pixels pooled across both dates, standardised."))


#  PHASE 3 - HOW MANY CLUSTERS DOES THE LANDSCAPE SUPPORT?
#  The elbow is a guide, not an oracle. Within-cluster sum of squares always
#  falls as k rises; what matters is where it stops falling steeply.

elbow <- map_dfr(k.candidates, function(k) {
  km <- kmeans(X.z, centers = k, nstart = 10, iter.max = 50)
  tibble(k = k,
         wss = km$tot.withinss,
         pct.explained = 100 * km$betweenss / km$totss)
})

cat("\n[3] CLUSTER SCAN\n\n")
elbow %>%
  mutate(wss = round(wss), pct.explained = round(pct.explained, 1)) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

png(file.path(path.out.cls, "kmeans_elbow.png"), width = 1400, height = 700, res = 150)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(elbow$k, elbow$wss, type = "b", pch = 19,
     xlab = "k", ylab = "within-cluster sum of squares", main = "Elbow")
abline(v = k.chosen, col = "#d7301f", lty = 2)
plot(elbow$k, elbow$pct.explained, type = "b", pch = 19,
     xlab = "k", ylab = "% variance explained", main = "Variance explained")
abline(v = k.chosen, col = "#d7301f", lty = 2)
dev.off()


#  PHASE 4 - FIT THE SHARED MODEL

km.fit <- kmeans(X.z, centers = k.chosen, nstart = kmeans.nstart, iter.max = 100)

## Centroids back in true units, so the clusters can actually be interpreted.
centers.true <- sweep(sweep(km.fit$centers, 2, feat.sd, "*"), 2, feat.mu, "+") %>%
  as_tibble() %>%
  mutate(across(everything(), to.reflectance)) %>%
  mutate(cluster = row_number(), .before = 1)

## Order clusters by NDVI so the numbering is interpretable rather than arbitrary:
## cluster 1 is the least vegetated, cluster k the most.
cluster.order <- order(centers.true$NDVI)
remap <- match(seq_len(k.chosen), cluster.order)   # old label -> new label

centers.ordered <- centers.true[cluster.order, ] %>%
  mutate(cluster.new = row_number(), .before = 1)

## A tentative reading of each cluster. These are SUGGESTIONS from the spectral
## signature, to be checked against imagery before they appear in any write-up -
## k-means knows nothing about land cover.
## NDVI alone does not separate these clusters: the three least-vegetated ones
## sit within 0.06 NDVI of each other and are distinguished by BRIGHTNESS and
## BURN sensitivity instead. A labelling rule that only looked at NDVI returned
## the same name three times, which told us nothing. NBR separates dark burn
## scars from bright bare soil at similar greenness, and SWIR1 separates bright
## mineral soil from dull dry grass.
suggest.label <- function(ndvi, ndwi, nbr, swir1) {
  case_when(
    ndwi >  0.00               ~ "water (check)",
    nbr  < -0.06               ~ "burnt / dark bare (check)",
    ndvi <  0.28 & swir1 > 0.28 ~ "bright bare / soil",
    ndvi <  0.28               ~ "sparse dry grass",
    ndvi <  0.35               ~ "dry grassland",
    ndvi <  0.45               ~ "grassland / shrub",
    ndvi <  0.62               ~ "wooded grassland",
    TRUE                       ~ "woodland / forest"
  )
}

centers.ordered <- centers.ordered %>%
  mutate(tentative = suggest.label(NDVI, NDWI, NBR, SWIR1))

cat("\n[4] CLUSTER CENTROIDS (true units, ordered by NDVI)\n\n")
centers.ordered %>%
  transmute(cluster = cluster.new,
            BLUE = round(BLUE, 3), GREEN = round(GREEN, 3), RED = round(RED, 3),
            NIR = round(NIR, 3), SWIR1 = round(SWIR1, 3), SWIR2 = round(SWIR2, 3),
            NDVI = round(NDVI, 3), NDWI = round(NDWI, 3), NBR = round(NBR, 3),
            tentative) %>%
  as.data.frame() %>%
  print(row.names = FALSE)


#  PHASE 5 - APPLY THE CENTROIDS TO BOTH DATES
#  Assignment is done with raster algebra rather than terra::predict() with an R
#  callback: for each cluster it is one squared-distance surface, then which.min
#  across them. That keeps the work inside terra's C++ rather than calling back
#  into R 32 million times.

centers.z <- km.fit$centers[cluster.order, , drop = FALSE]   # standardised, reordered

classify.kmeans <- function(r) {

  ## Standardise the raster with the SAME constants used for the fit.
  r.z <- (r[[feat.bands]] - feat.mu) / feat.sd

  ## One squared-distance layer per cluster.
  d <- lapply(seq_len(nrow(centers.z)), function(j) {
    sum((r.z - as.numeric(centers.z[j, ]))^2)
  })

  terra::which.min(terra::rast(d))
}

cls.a <- classify.kmeans(r.a.f)
cls.b <- classify.kmeans(r.b.f)

names(cls.a) <- glue("class_{yr.a}")
names(cls.b) <- glue("class_{yr.b}")

terra::writeRaster(cls.a, file.path(path.out.cls, glue("kmeans_k{k.chosen}_{yr.a}.tif")),
                   overwrite = TRUE, datatype = "INT1U")
terra::writeRaster(cls.b, file.path(path.out.cls, glue("kmeans_k{k.chosen}_{yr.b}.tif")),
                   overwrite = TRUE, datatype = "INT1U")

say(glue("\n[5] classified both dates; rasters written to {path.out.cls}/"))


#  PHASE 6 - CLASS AREAS AND THE CHANGE MATRIX

px.km2 <- prod(terra::res(cls.a)) / 1e6   # 30 x 30 m = 0.0009 km2

area.by.class <- function(cls, label) {
  terra::freq(cls) %>%
    as_tibble() %>%
    transmute(cluster = value, date = label, km2 = count * px.km2)
}

areas <- bind_rows(area.by.class(cls.a, as.character(yr.a)),
                   area.by.class(cls.b, as.character(yr.b))) %>%
  left_join(centers.ordered %>% select(cluster = cluster.new, tentative), by = "cluster") %>%
  pivot_wider(names_from = date, values_from = km2, names_prefix = "km2_") %>%
  mutate(
    change.km2 = .data[[glue("km2_{yr.b}")]] - .data[[glue("km2_{yr.a}")]],
    change.pct = 100 * change.km2 / .data[[glue("km2_{yr.a}")]]
  )

cat("\n[6] CLASS AREAS (km2, well-supported pixels only)\n\n")
areas %>%
  mutate(across(where(is.numeric), ~round(.x, 1))) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

## 6.2 The change matrix. Rows = date A class, columns = date B class.
xt <- terra::crosstab(c(cls.a, cls.b), useNA = FALSE)
xt.km2 <- as.matrix(xt) * px.km2

cat("\n[6.2] CHANGE MATRIX (km2). Rows =", yr.a, " Columns =", yr.b, "\n")
cat("      Diagonal = unchanged cluster membership.\n\n")
print(round(xt.km2, 1))

total.km2     <- sum(xt.km2)
unchanged.km2 <- sum(diag(xt.km2))
changed.km2   <- total.km2 - unchanged.km2

## 6.3 HOW BIG IS THE CHANGE, NOT JUST HOW MUCH?
##
## A raw "% changed" from clustered data flatters itself badly. Clusters are
## ordered by NDVI, so a pixel moving between NEIGHBOURING clusters has shifted
## one step along a greenness gradient - which in a semi-arid rangeland is what
## a slightly different rainfall year, or a burn scar at a different age, looks
## like. It is not land-cover change.
##
## Splitting the off-diagonal by how far a pixel moved separates the plausible
## signal from the churn. A jump of two or more ranks is the part worth
## investigating; single-step movement should be treated as noise until
## something demonstrates otherwise.
rank.gap <- abs(row(xt.km2) - col(xt.km2))

adjacent.km2     <- sum(xt.km2[rank.gap == 1])
nonadjacent.km2  <- sum(xt.km2[rank.gap >= 2])
big.jump.km2     <- sum(xt.km2[rank.gap >= 3])

say(glue("
      analysed area      : {format(round(total.km2, 1), big.mark = ',')} km2
      unchanged          : {format(round(unchanged.km2, 1), big.mark = ',')} km2  ({round(100*unchanged.km2/total.km2, 1)}%)
      changed            : {format(round(changed.km2, 1), big.mark = ',')} km2  ({round(100*changed.km2/total.km2, 1)}%)

      of which:
        one cluster step : {format(round(adjacent.km2, 1), big.mark = ',')} km2  ({round(100*adjacent.km2/total.km2, 1)}%)  <- treat as noise
        two or more      : {format(round(nonadjacent.km2, 1), big.mark = ',')} km2  ({round(100*nonadjacent.km2/total.km2, 1)}%)  <- worth investigating
        three or more    : {format(round(big.jump.km2, 1), big.mark = ',')} km2  ({round(100*big.jump.km2/total.km2, 1)}%)  <- strongest candidates
"))

change.summary <- tibble(
  metric = c("total_km2", "unchanged_km2", "changed_km2",
             "adjacent_step_km2", "two_or_more_steps_km2", "three_or_more_steps_km2"),
  km2    = round(c(total.km2, unchanged.km2, changed.km2,
                   adjacent.km2, nonadjacent.km2, big.jump.km2), 2),
  pct    = round(100 * c(total.km2, unchanged.km2, changed.km2,
                         adjacent.km2, nonadjacent.km2, big.jump.km2) / total.km2, 2)
)
readr::write_csv(change.summary, file.path(path.out.cls, "kmeans_change_summary.csv"))

readr::write_csv(as_tibble(xt.km2, rownames = glue("class_{yr.a}")),
                 file.path(path.out.cls, "kmeans_change_matrix_km2.csv"))
readr::write_csv(areas,           file.path(path.out.cls, "kmeans_class_areas.csv"))
readr::write_csv(centers.ordered, file.path(path.out.cls, "kmeans_centroids.csv"))
readr::write_csv(elbow,           file.path(path.out.cls, "kmeans_elbow.csv"))


#  PHASE 7 - BINARY CHANGE MASK AND MAPS

## as.int() + explicit levels: plotting a LOGICAL SpatRaster as classes makes
## terra coerce the labels and warn ("NAs introduced by coercion"). Giving it a
## proper categorical raster is both quieter and correctly labelled.
changed <- terra::as.int(cls.a != cls.b)
names(changed) <- "changed"
levels(changed) <- data.frame(value = c(0, 1), label = c("same", "changed"))
terra::writeRaster(changed, file.path(path.out.cls, "kmeans_change_mask.tif"),
                   overwrite = TRUE, datatype = "INT1U")

## How far each pixel moved along the greenness ordering - far more useful than
## the binary mask, because it shows WHERE the big transitions are rather than
## lumping a one-step drift in with a woodland-to-bare jump.
steps <- terra::as.int(abs(cls.a - cls.b))
names(steps) <- "cluster_steps"
terra::writeRaster(steps, file.path(path.out.cls, "kmeans_change_steps.tif"),
                   overwrite = TRUE, datatype = "INT1U")

cls.pal <- scico::scico(k.chosen, palette = "batlow")

png(file.path(path.out.cls, "kmeans_pair.png"), width = 1900, height = 1000, res = 150)
par(mfrow = c(1, 2), mar = c(2, 2, 3, 5))
terra::plot(cls.a, col = cls.pal, type = "classes",
            main = glue("k-means k={k.chosen} - {yr.a}"))
terra::plot(aoi.parts.v, add = TRUE, border = "grey20", lwd = 0.5)
terra::plot(cls.b, col = cls.pal, type = "classes",
            main = glue("k-means k={k.chosen} - {yr.b}"))
terra::plot(aoi.parts.v, add = TRUE, border = "grey20", lwd = 0.5)
dev.off()

png(file.path(path.out.cls, "kmeans_change_mask.png"), width = 1900, height = 1000, res = 150)
par(mfrow = c(1, 2), mar = c(2, 2, 3, 5))
terra::plot(changed, col = c("grey92", "#d7301f"),
            main = glue("Cluster membership changed, {yr.a} to {yr.b}"))
terra::plot(aoi.parts.v, add = TRUE, border = "grey20", lwd = 0.5)
terra::plot(steps, col = c("grey92", scico::scico(k.chosen - 1, palette = "lajolla",
                                                  direction = -1)),
            main = "Cluster steps moved (0 = none)")
terra::plot(aoi.parts.v, add = TRUE, border = "grey20", lwd = 0.5)
dev.off()

say(glue("[7] wrote change mask and 2 figures to {path.out.cls}/"))


#  PHASE 8 - IS THE CHANGE SPATIALLY COHERENT, OR IS IT SPECKLE?
#  The decisive test, and the cheapest one to skip by accident.
#
#  Real land-cover change has EDGES: a cleared field, an expanding settlement, a
#  thinned woodland are contiguous patches. Pixel noise, by contrast, is
#  scattered at random. Two landscapes can report the same "% changed" while one
#  is a genuine signal and the other is nothing at all.
#
#  The test: among changed pixels, count how many of the 8 neighbours also
#  changed, and compare that with what pure spatial randomness would give. If a
#  fraction p of the landscape changed and change were randomly scattered, a
#  changed pixel would average 8p changed neighbours. Substantially more than
#  that means clustering (real patches); close to it means speckle.

chg01 <- terra::as.int(cls.a != cls.b)

nb.sum <- terra::focal(chg01, w = matrix(1, 3, 3), fun = "sum",
                       na.policy = "omit", na.rm = TRUE)

## Neighbours that changed, excluding the focal pixel itself, over changed pixels.
nb.of.changed <- terra::mask(nb.sum - chg01, chg01, maskvalues = c(0, NA))

obs.neighbours <- terra::global(nb.of.changed, "mean", na.rm = TRUE)[1, 1]
p.changed      <- changed.km2 / total.km2
exp.neighbours <- 8 * p.changed

pct.isolated <- 100 * terra::global(nb.of.changed <= 1, "mean", na.rm = TRUE)[1, 1]

coherence.ratio <- obs.neighbours / exp.neighbours

say(glue("
[8] SPATIAL COHERENCE OF THE CHANGE
    changed fraction           : {round(100 * p.changed, 1)}%
    expected changed neighbours
      if change were random    : {round(exp.neighbours, 2)} of 8
    observed changed neighbours: {round(obs.neighbours, 2)} of 8
    coherence ratio            : {round(coherence.ratio, 2)}  (1.0 = indistinguishable from random)
    changed pixels with <= 1
      changed neighbour        : {round(pct.isolated, 1)}%  (isolated specks)
"))

if (coherence.ratio < 1.25) {
  cat("    VERDICT: the change is close to randomly scattered. It is dominated by\n")
  cat("    pixel-level noise - residual cloud, burn-scar age, seasonal greenness -\n")
  cat("    NOT by land-cover change. Do not report the headline % as change.\n")
} else if (coherence.ratio < 1.6) {
  cat("    VERDICT: mixed. Some genuine patches exist inside a great deal of\n")
  cat("    speckle. A minimum mapping unit is needed before anything is quoted.\n")
} else {
  cat("    VERDICT: the change is strongly clustered, not random noise.\n")
  cat("\n")
  cat("    But clustered is NOT the same as land-cover change. A burn scar is a\n")
  cat("    coherent patch with hard edges, and so is a flush of green after a\n")
  cat("    localised storm. This test rules out pixel noise; it cannot tell fire\n")
  cat("    from clearing. The next step is to separate those two, by overlaying\n")
  cat("    the change patches on the MODIS burn layers from 01b and on the NBR\n")
  cat("    difference, and keeping only what neither explains.\n")
}

readr::write_csv(
  tibble(metric = c("pct_changed", "expected_neighbours_if_random",
                    "observed_neighbours", "coherence_ratio", "pct_isolated"),
         value  = round(c(100 * p.changed, exp.neighbours, obs.neighbours,
                          coherence.ratio, pct.isolated), 3)),
  file.path(path.out.cls, "kmeans_change_coherence.csv"))


#  PHASE 8 - WHAT THIS DOES AND DOES NOT ESTABLISH

cat("\n================================================================\n")
cat(" DONE - but read this before quoting any of the numbers above.\n")
cat("\n")
cat(" - Clusters are SPECTRAL, not thematic. The `tentative` labels are read\n")
cat("   off centroid signatures and must be checked against imagery before use.\n")
cat(" - 'Changed' here means a pixel moved between clusters. In a fire-driven\n")
cat("   rangeland much of that is burn scars and seasonal greenness, not land\n")
cat("   use change. Compare the change mask against the burn layers before\n")
cat("   drawing any conclusion about land cover.\n")
cat(" - Thin pixels were excluded; the areas describe well-supported pixels only.\n")
cat(" - Next: a supervised Random Forest with real training data, and the\n")
cat("   training-similarity check on cross-date transfer (CLAUDE.md caveat 2).\n")
cat("================================================================\n\n")

sessionInfo()
# =============================================================================
#  KIDEPO LANDSCAPE — LCLUC
#  Script 03b — Did the k-means clustering actually work?
#
#  WHY THIS EXISTS
#  Script 03 ran, produced spatially coherent clusters and a change matrix, and
#  everything downstream was built on it. But "it ran" is not "it worked", and
#  three things were never checked:
#
#    1. CONVERGENCE. The elbow scan in 03 emitted five "Quick-TRANSfer stage
#       steps exceeded maximum" warnings. That message means Hartigan-Wong hit
#       its iteration ceiling without settling, so those fits are not at a
#       local optimum. Whether the FINAL k=6 fit converged was never checked.
#
#    2. k WAS NEVER ACTUALLY CHOSEN. k.chosen was set to 6 before the scan ran,
#       and the comment telling the next reader to revisit it after seeing the
#       elbow was never acted on. The elbow curve is smooth - 71, 79, 82, 84,
#       86, 87, 88% variance explained for k = 3..9 - with no obvious break, so
#       the choice needs a criterion beyond eyeballing.
#
#    3. SEPARABILITY. Three of the six clusters sit within 0.055 NDVI of each
#       other, and two within 0.007. That is what splitting a CONTINUUM looks
#       like rather than finding distinct classes, and it matters: pixels near
#       an arbitrary boundary flip between dates for no physical reason, which
#       is a strong candidate for the inflated change figures.
#
#  Silhouette is computed in its simplified, centroid-based form: for each point
#  the distance to its own centroid against the distance to the nearest other
#  centroid. True silhouette needs an n-by-n distance matrix, which at 60,000
#  points is 28 GB. The simplified version is the standard substitute and
#  answers the same question.
#
#  HOW TO RUN — from the repo root, in PowerShell:
#      & "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" scripts\R\03b_kmeans_diagnostics.R
#
#  OUTPUTS -> outputs/classification/
#
#  Author: Isaac Mureithi - IMureithi@awf.org
# =============================================================================

source("scripts/R/00_functions.R")

terra::terraOptions(progress = 0)
set.seed(7)   # deliberately NOT 42: if the conclusions depend on the draw,
              # they are not conclusions

k.candidates  <- 3:10
n.per.date    <- 30000
kmeans.nstart <- 25
feat.bands    <- c(comp.bands.sr, comp.bands.index)

cat("\n================================================================\n")
cat(" KIDEPO LCLUC — 03b k-means diagnostics\n")
cat(" run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n")


# =============================================================================
#  PHASE 1 — REBUILD THE SAME FEATURE SPACE AS 03
# =============================================================================

comp <- load.composite.pair()
aoi  <- load.aoi()
aoi.union.v <- terra::vect(aoi$union)
mask.to.aoi <- function(r) terra::mask(terra::crop(r, aoi.union.v), aoi.union.v)

r.a <- mask.to.aoi(comp$a)
r.b <- mask.to.aoi(comp$b)

well.supported <- (r.a[["clear_obs"]] >= comp.min.obs) &
                  (r.b[["clear_obs"]] >= comp.min.obs)

r.a.f <- terra::mask(r.a[[feat.bands]], well.supported, maskvalues = c(FALSE, NA))
r.b.f <- terra::mask(r.b[[feat.bands]], well.supported, maskvalues = c(FALSE, NA))

X <- bind_rows(
  as_tibble(terra::spatSample(r.a.f, n.per.date, "random", na.rm = TRUE, warn = FALSE)),
  as_tibble(terra::spatSample(r.b.f, n.per.date, "random", na.rm = TRUE, warn = FALSE))
) %>% drop_na() %>% as.matrix()

X.z <- scale(X)

say(glue("[1] pooled sample: {format(nrow(X.z), big.mark = ',')} pixels, {ncol(X.z)} standardised features"))


# =============================================================================
#  PHASE 2 — FIT EACH k, AND RECORD WHETHER IT CONVERGED
# =============================================================================
#  R signals the Quick-TRANSfer failure as a WARNING, not an error, so a script
#  that ignores warnings will happily use a fit that never settled. Each fit is
#  wrapped so the warning is captured rather than printed and forgotten.

fit.one <- function(k) {

  warned <- character(0)
  km <- withCallingHandlers(
    kmeans(X.z, centers = k, nstart = kmeans.nstart, iter.max = 100),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  ## Simplified silhouette, centroid-based.
  d2 <- apply(km$centers, 1, function(cen) rowSums(sweep(X.z, 2, cen)^2))
  d  <- sqrt(d2)                                   # n x k distances
  own <- d[cbind(seq_len(nrow(d)), km$cluster)]
  d.other <- d
  d.other[cbind(seq_len(nrow(d)), km$cluster)] <- Inf
  nearest <- apply(d.other, 1, min)
  sil <- (nearest - own) / pmax(own, nearest)

  tibble(
    k             = k,
    converged     = length(warned) == 0,
    pct.explained = 100 * km$betweenss / km$totss,
    mean.sil      = mean(sil),
    min.cluster.n = min(km$size),
    pct.sil.below.025 = 100 * mean(sil < 0.25)
  )
}

diag.tbl <- map_dfr(k.candidates, fit.one)

cat("\n[2] FIT QUALITY BY k\n\n")
diag.tbl %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\n    converged         FALSE means Hartigan-Wong hit its iteration ceiling.\n")
cat("    mean.sil          above ~0.5 is strong structure, 0.25-0.5 weak but real,\n")
cat("                      below 0.25 means no substantial structure was found.\n")
cat("    pct.sil.below.025 share of pixels sitting ambiguously between clusters.\n")

best.k <- diag.tbl %>% filter(converged) %>% slice_max(mean.sil, n = 1) %>% pull(k)
say(glue("\n    Highest mean silhouette among converged fits: k = {best.k}"))


# =============================================================================
#  PHASE 3 — WHERE k=6 IS WEAK
# =============================================================================

km6 <- withCallingHandlers(
  kmeans(X.z, centers = 6, nstart = kmeans.nstart, iter.max = 100),
  warning = function(w) invokeRestart("muffleWarning")
)

centers.true <- sweep(km6$centers, 2, attr(X.z, "scaled:scale"), "*") %>%
  sweep(2, attr(X.z, "scaled:center"), "+") %>%
  as_tibble() %>%
  mutate(across(everything(), to.reflectance))

ord <- order(centers.true$NDVI)

## Distance between every pair of centroids, in standardised feature space.
cen.d <- as.matrix(dist(km6$centers[ord, ]))
dimnames(cen.d) <- list(seq_len(6), seq_len(6))

cat("\n[3] DISTANCE BETWEEN CLUSTER CENTROIDS (standardised units, NDVI order)\n\n")
print(round(cen.d, 2))

nearest.pair <- which(cen.d == min(cen.d[cen.d > 0]), arr.ind = TRUE)[1, ]

say(glue("
    Closest pair: clusters {nearest.pair[1]} and {nearest.pair[2]}, {round(min(cen.d[cen.d > 0]), 2)} apart.
    For comparison the largest separation is {round(max(cen.d), 2)}.

    A ratio this wide means the six clusters are not six comparable classes: a
    few are well separated and the rest are slices of one continuum."))

cat("\n    Cluster NDVI in order:",
    paste(round(sort(centers.true$NDVI), 3), collapse = ", "), "\n")


# =============================================================================
#  PHASE 4 — OUTPUT
# =============================================================================

readr::write_csv(diag.tbl, file.path(path.out.cls, "kmeans_diagnostics.csv"))

png(file.path(path.out.cls, "kmeans_diagnostics.png"), width = 1500, height = 700, res = 150)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(diag.tbl$k, diag.tbl$mean.sil, type = "b", pch = 19,
     xlab = "k", ylab = "mean simplified silhouette",
     main = "Cluster structure", ylim = c(0, max(0.6, max(diag.tbl$mean.sil))))
abline(h = c(0.25, 0.5), lty = 2, col = c("#d7301f", "#1b7837"))
plot(diag.tbl$k, diag.tbl$pct.explained, type = "b", pch = 19,
     xlab = "k", ylab = "% variance explained", main = "Elbow (for reference)")
abline(v = 6, lty = 2, col = "grey50")
dev.off()

cat("\n================================================================\n")
cat(" DONE — read [2] before trusting anything built on k=6.\n")
cat("================================================================\n\n")

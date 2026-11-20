#  KIDEPO LANDSCAPE - LCLUC
#  Script 04 - Separate fire from candidate land-cover change
#
#  THE PROBLEM THIS SOLVES
#  Script 03 established that the change between the two dates is spatially
#  clustered rather than pixel noise (coherence ratio 1.62). That rules out one
#  explanation and leaves two standing, because a burn scar is every bit as
#  coherent as a cleared field. Roughly 20% of the landscape carries scars in
#  each date, and crucially they are not the SAME 20%: fire moves around the
#  landscape from year to year, so a pixel burnt in one date and not the other
#  will register as change while nothing about its land cover has altered.
#
#  THE APPROACH, AND ITS HONEST LIMIT
#  We have no burnt-area raster locally - MODIS burned area lives in Earth
#  Engine, and 01b exported only landscape totals. What those totals give us is
#  a CALIBRATION: MODIS measured 19.7% of the AOI burnt by the end of the 2019
#  window and 22.8% by the end of 2024. So for each date we take the NBR
#  threshold that flags exactly that fraction of the AOI as burn-like.
#
#  This is calibration, not detection. It assumes the lowest-NBR fraction of the
#  landscape is the burnt fraction, which will be wrong at the margins: deep
#  shadow, bare rock and the driest grass also sit at low NBR. It is defensible
#  because the fraction comes from an independent sensor rather than from a
#  threshold chosen to produce a pleasing answer, and because what matters here
#  is the DIFFERENCE between two dates treated identically.
#
#  To validate it properly, export the MODIS masks as rasters from GEE and
#  cross-tabulate. That is the natural next step, not a substitute for this one.
#
#  HOW TO RUN - from the repo root, in PowerShell:
#      & "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" scripts\R\04_separate_fire.R
#
#  Requires scripts/R/03_kmeans_firstlook.R to have been run.
#  OUTPUTS -> outputs/classification/
#
#  Author: Isaac Mureithi - IMureithi@awf.org

source("scripts/R/00_functions.R")

terra::terraOptions(progress = 0)
set.seed(42)


#  PHASE 0 - PARAMETERS

k.chosen <- 6          # must match what 03 was run with

## MODIS cumulative burnt fraction by the end of each window, from
## outputs of scripts/gee/01b_burn_season_profile.js. These calibrate the
## per-date NBR thresholds; they are measurements, not guesses.
burn.frac <- c("2019" = 0.197, "2024" = 0.228)

## A pixel counts as substantively changed if it moved two or more clusters
## along the greenness ordering. One step is within the range a different
## rainfall year produces (see 03).
min.steps <- 2

## Minimum mapping unit. Land-cover change smaller than this is not credible at
## 30 m from a two-date comparison, whatever the pixels say.
mmu.ha     <- 1.0
mmu.pixels <- ceiling(mmu.ha * 10000 / prod(terra::res(terra::rast(composite.path(comp.years[1])))))

n.quantile.sample <- 500000   # pixels sampled to estimate NBR thresholds

cat("\n================================================================\n")
cat(" KIDEPO LCLUC - 04 Separating fire from land-cover change\n")
cat(" run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n")


#  PHASE 1 - LOAD

comp <- load.composite.pair()
yr.a <- comp$years[1]
yr.b <- comp$years[2]

aoi <- load.aoi()
aoi.union.v <- terra::vect(aoi$union)
aoi.parts.v <- terra::vect(aoi$parts)

mask.to.aoi <- function(r) terra::mask(terra::crop(r, aoi.union.v), aoi.union.v)

## Classified layers from 03.
f.steps <- file.path(path.out.cls, "kmeans_change_steps.tif")
f.cls.a <- file.path(path.out.cls, glue("kmeans_k{k.chosen}_{yr.a}.tif"))
f.cls.b <- file.path(path.out.cls, glue("kmeans_k{k.chosen}_{yr.b}.tif"))

if (!all(file.exists(f.steps, f.cls.a, f.cls.b))) {
  stop("Outputs from 03_kmeans_firstlook.R not found - run that script first.")
}

steps <- terra::rast(f.steps)
cls.a <- terra::rast(f.cls.a)
cls.b <- terra::rast(f.cls.b)

nbr.a <- mask.to.aoi(comp$a[["NBR"]])
nbr.b <- mask.to.aoi(comp$b[["NBR"]])

## Restrict everything to the pixels 03 actually classified.
valid <- !is.na(steps)
nbr.a <- terra::mask(nbr.a, valid, maskvalues = c(FALSE, NA))
nbr.b <- terra::mask(nbr.b, valid, maskvalues = c(FALSE, NA))

say(glue("[1] loaded classified layers and NBR for {yr.a} and {yr.b}"))


#  PHASE 2 - CALIBRATE A BURN-LIKE THRESHOLD PER DATE
#  Quantiles are estimated from a large random sample rather than every cell:
#  an exact quantile over 32 million values would pull the whole band into
#  memory for a number that a half-million-pixel sample pins down to three
#  decimal places.

samp.nbr <- terra::spatSample(c(nbr.a, nbr.b), n.quantile.sample,
                              method = "random", na.rm = TRUE, warn = FALSE)
names(samp.nbr) <- c("a", "b")

thr.a <- quantile(samp.nbr$a, probs = burn.frac[[as.character(yr.a)]], na.rm = TRUE)
thr.b <- quantile(samp.nbr$b, probs = burn.frac[[as.character(yr.b)]], na.rm = TRUE)

burn.a <- terra::as.int(nbr.a < thr.a)
burn.b <- terra::as.int(nbr.b < thr.b)

say(glue("
[2] BURN-LIKE CALIBRATION
    {yr.a}: MODIS burnt {100 * burn.frac[[as.character(yr.a)]]}% -> NBR threshold {round(to.reflectance(thr.a), 4)}
    {yr.b}: MODIS burnt {100 * burn.frac[[as.character(yr.b)]]}% -> NBR threshold {round(to.reflectance(thr.b), 4)}

    The thresholds differ between dates because the whole NBR distribution
    shifts between years. Flagging a FRACTION rather than a fixed NBR value is
    what keeps the two dates comparable.
"))


#  PHASE 3 - THE FIRE TRANSITION
#  Four states. Only the two middle ones can generate spurious change.
#    0  burn-like in neither date        - fire cannot explain any change here
#    1  burn-like in A only              - scar present in A, gone by B
#    2  burn-like in B only              - scar present in B, absent in A
#    3  burn-like in both                - persistently low NBR; fire is a
#                                          constant here, so it explains little

fire.state <- burn.a * 1 + burn.b * 2
names(fire.state) <- "fire_state"

fire.labels <- c("neither", glue("{yr.a} only"), glue("{yr.b} only"), "both")

px.km2 <- prod(terra::res(steps)) / 1e6

fire.tab <- terra::freq(fire.state) %>%
  as_tibble() %>%
  transmute(state = value,
            label = fire.labels[value + 1],
            km2   = count * px.km2) %>%
  mutate(pct = 100 * km2 / sum(km2))

cat("\n[3] BURN-LIKE STATE ACROSS THE TWO DATES\n\n")
fire.tab %>%
  mutate(across(where(is.numeric), ~round(.x, 1))) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

fire.switched <- terra::as.int(fire.state == 1 | fire.state == 2)
pct.switched  <- 100 * terra::global(fire.switched, "mean", na.rm = TRUE)[1, 1]

say(glue("
    Burn state SWITCHED between dates over {round(pct.switched, 1)}% of the area.
    That is the share of the landscape where fire alone can manufacture an
    apparent change, with no land-cover change whatsoever."))


#  PHASE 4 - ATTRIBUTE THE SUBSTANTIVE CHANGE

substantive <- terra::as.int(steps >= min.steps)
sub.km2 <- terra::global(substantive, "sum", na.rm = TRUE)[1, 1] * px.km2

## Split substantive change by whether the burn state also switched.
fire.explained <- terra::as.int(substantive == 1 & fire.switched == 1)
residual       <- terra::as.int(substantive == 1 & fire.switched == 0)

fire.km2 <- terra::global(fire.explained, "sum", na.rm = TRUE)[1, 1] * px.km2
res.km2  <- terra::global(residual,       "sum", na.rm = TRUE)[1, 1] * px.km2

say(glue("
[4] ATTRIBUTION OF SUBSTANTIVE CHANGE ({min.steps}+ cluster steps)
    substantive change    : {format(round(sub.km2, 1), big.mark = ',')} km2
    burn state switched   : {format(round(fire.km2, 1), big.mark = ',')} km2  ({round(100*fire.km2/sub.km2, 1)}% of it)
    burn state unchanged  : {format(round(res.km2, 1), big.mark = ',')} km2  ({round(100*res.km2/sub.km2, 1)}% of it)

    The second line is change that fire can account for. The third is the
    residual: the candidate land-cover change, before any size filtering."))


#  PHASE 4.5 - IS THE FIRE FILTER ACTUALLY DOING ANYTHING?
#  A filter that removes a chunk of the change looks productive whether or not
#  it is removing the RIGHT chunk. If burn-state switching genuinely drives
#  change, then pixels that switched should change far more often than pixels
#  that did not. If the two rates are similar, the filter is deleting area
#  roughly at random and the "residual" is not a fire-corrected anything.
#
#  Two readings of the same question:
#    enrichment  observed overlap vs the overlap expected if the two layers
#                were independent. 1.0 means no association at all.
#    risk ratio  how much more likely a switched pixel is to have changed.

area.total   <- sum(fire.tab$km2)
switched.km2 <- fire.tab %>% filter(state %in% c(1, 2)) %>% pull(km2) %>% sum()

p.sub            <- sub.km2 / area.total
expected.overlap <- switched.km2 * p.sub
enrichment       <- fire.km2 / expected.overlap

risk.switched <- fire.km2 / switched.km2
risk.not      <- (sub.km2 - fire.km2) / (area.total - switched.km2)
risk.ratio    <- risk.switched / risk.not

say(glue("
[4.5] ASSOCIATION BETWEEN BURN SWITCHING AND CHANGE
    overlap if independent : {round(expected.overlap, 1)} km2
    overlap observed       : {round(fire.km2, 1)} km2
    enrichment             : {round(enrichment, 2)}x   (1.0 = no association)

    change rate where burn switched   : {round(100 * risk.switched, 1)}%
    change rate where it did not      : {round(100 * risk.not, 1)}%
    risk ratio                        : {round(risk.ratio, 2)}x
"))

if (enrichment < 1.5) {
  cat("    VERDICT: WEAK. Burn switching barely predicts change, so this filter\n")
  cat("    is removing area close to indiscriminately. Either the NBR-quantile\n")
  cat("    proxy is not finding the actual scars - likely, since a scar one to\n")
  cat("    four months old may have greened past the lowest NBR quantile while\n")
  cat("    rock, shadow and water margins sit permanently inside it - or fire is\n")
  cat("    not the main driver of the change. Do NOT present the residual as\n")
  cat("    land-cover change on this evidence. Export the MODIS burned-area\n")
  cat("    rasters from GEE and redo the attribution against real scars.\n")
} else {
  cat("    VERDICT: the burn layer carries real information about where change\n")
  cat("    occurs, so removing it is a meaningful correction.\n")
}


#  PHASE 5 - MINIMUM MAPPING UNIT
#  Scattered single pixels surviving the fire filter are not land-cover change;
#  they are the ragged edge of two imperfect classifications. patches() groups
#  contiguous residual pixels, and anything below the MMU is dropped.

say(glue("\n[5] applying a {mmu.ha} ha minimum mapping unit ({mmu.pixels} pixels, 8-connected)"))

res.masked <- terra::subst(residual, from = 0, to = NA)
patch.id   <- terra::patches(res.masked, directions = 8, zeroAsNA = TRUE)

patch.sizes <- terra::freq(patch.id) %>%
  as_tibble() %>%
  transmute(patch = value, px = count)

keep.patches <- patch.sizes %>% filter(px >= mmu.pixels)

say(glue("    residual patches        : {format(nrow(patch.sizes), big.mark = ',')}"))
say(glue("    patches >= {mmu.ha} ha       : {format(nrow(keep.patches), big.mark = ',')}"))

residual.mmu <- terra::ifel(patch.id %in% keep.patches$patch, 1, 0)
residual.mmu <- terra::mask(residual.mmu, valid, maskvalues = c(FALSE, NA))
residual.mmu <- terra::subst(residual.mmu, from = NA, to = 0)
residual.mmu <- terra::mask(residual.mmu, valid, maskvalues = c(FALSE, NA))

res.mmu.km2 <- terra::global(residual.mmu, "sum", na.rm = TRUE)[1, 1] * px.km2

say(glue("
    candidate land-cover change after filtering:
      {format(round(res.mmu.km2, 1), big.mark = ',')} km2
      {round(100 * res.mmu.km2 / sub.km2, 1)}% of the substantive change
      {round(100 * res.mmu.km2 / (sum(fire.tab$km2)), 2)}% of the analysed area"))


#  PHASE 6 - IS THE RESIDUAL COHERENT?
#  The same test script 03 applied to raw change, now applied to what survives.
#  If the residual is coherent, it is worth taking to a supervised
#  classification. If it has collapsed towards randomness, then the clustering
#  seen in 03 was the fire, and this exercise has said so.

coh.raw <- change.coherence(substantive)
coh.res <- change.coherence(residual.mmu)

cat("\n[6] SPATIAL COHERENCE, BEFORE AND AFTER REMOVING FIRE\n\n")
bind_rows(
  coh.raw %>% mutate(layer = glue("substantive change ({min.steps}+ steps)"), .before = 1),
  coh.res %>% mutate(layer = "residual after fire + MMU",                     .before = 1)
) %>%
  mutate(across(where(is.numeric), ~round(.x, 2))) %>%
  as.data.frame() %>%
  print(row.names = FALSE)


#  PHASE 7 - WHAT CHANGED INTO WHAT, IN THE RESIDUAL

cls.a.res <- terra::mask(cls.a, residual.mmu, maskvalues = c(0, NA))
cls.b.res <- terra::mask(cls.b, residual.mmu, maskvalues = c(0, NA))

xt.res <- terra::crosstab(c(cls.a.res, cls.b.res), useNA = FALSE)
xt.res.km2 <- as.matrix(xt.res) * px.km2

cat("\n[7] RESIDUAL CHANGE MATRIX (km2). Rows =", yr.a, " Columns =", yr.b, "\n")
cat("    Cluster 1 is least vegetated, cluster", k.chosen, "most.\n\n")
print(round(xt.res.km2, 1))

gain.loss <- tibble(
  cluster = seq_len(k.chosen),
  lost.km2   = rowSums(xt.res.km2),
  gained.km2 = colSums(xt.res.km2)
) %>%
  mutate(net.km2 = gained.km2 - lost.km2)

cat("\n    Net movement per cluster in the residual (km2):\n\n")
gain.loss %>%
  mutate(across(where(is.numeric), ~round(.x, 1))) %>%
  as.data.frame() %>%
  print(row.names = FALSE)


#  PHASE 8 - OUTPUTS

terra::writeRaster(fire.state, file.path(path.out.cls, "fire_state.tif"),
                   overwrite = TRUE, datatype = "INT1U")
terra::writeRaster(residual.mmu, file.path(path.out.cls, "change_residual_mmu.tif"),
                   overwrite = TRUE, datatype = "INT1U")

readr::write_csv(fire.tab,  file.path(path.out.cls, "fire_state_areas.csv"))
readr::write_csv(gain.loss, file.path(path.out.cls, "residual_gain_loss.csv"))
readr::write_csv(
  tibble(metric = c("substantive_km2", "fire_explained_km2", "residual_km2",
                    "residual_after_mmu_km2", "pct_burn_state_switched",
                    "coherence_substantive", "coherence_residual"),
         value  = round(c(sub.km2, fire.km2, res.km2, res.mmu.km2, pct.switched,
                          coh.raw$coherence.ratio, coh.res$coherence.ratio), 3)),
  file.path(path.out.cls, "fire_attribution_summary.csv"))

png(file.path(path.out.cls, "fire_vs_residual.png"), width = 1900, height = 1000, res = 150)
par(mfrow = c(1, 2), mar = c(2, 2, 3, 6))
terra::plot(fire.state, col = c("grey92", "#f1a340", "#998ec3", "#4d4d4d"),
            type = "classes", levels = fire.labels,
            main = glue("Burn-like state, {yr.a} vs {yr.b}"))
terra::plot(aoi.parts.v, add = TRUE, border = "grey20", lwd = 0.5)
terra::plot(residual.mmu, col = c("grey92", "#1b7837"),
            type = "classes", levels = c("no", "candidate"),
            main = glue("Candidate land-cover change (>= {mmu.ha} ha)"))
terra::plot(aoi.parts.v, add = TRUE, border = "grey20", lwd = 0.5)
dev.off()

say(glue("\n[8] wrote rasters, 3 tables and 1 figure to {path.out.cls}/"))

cat("\n================================================================\n")
cat(" DONE. What this does and does not settle:\n")
cat("\n")
cat("  - The burn-like layer is CALIBRATED to MODIS totals, not detected.\n")
cat("    Validate it by exporting the MODIS masks as rasters from GEE and\n")
cat("    cross-tabulating against fire_state.tif.\n")
cat("  - The residual is a CANDIDATE set, not a change map. It still contains\n")
cat("    rainfall-driven greenness differences that fire filtering cannot touch.\n")
cat("  - Check the residual patches against high-resolution imagery before any\n")
cat("    of it is described as land-use change.\n")
cat("  - Then: supervised Random Forest, with a burnt class in the scheme.\n")
cat("================================================================\n\n")

sessionInfo()
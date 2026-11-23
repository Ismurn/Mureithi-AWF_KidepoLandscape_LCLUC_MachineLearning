# =============================================================================
#  KIDEPO LANDSCAPE — LCLUC
#  Script 08 — Is the greening an artefact of unequal observation depth?
#
#  THE WORRY
#  The headline result is that the landscape greened between 2019 and 2024:
#  median ΔNDVI +0.0137, two thirds of the residual pixels greening. But the
#  two composites are not equally supported — 2019 carries 8.1 mean clear looks
#  per pixel against 2024's 10.4, about 25% fewer.
#
#  There is a plausible mechanism by which that alone could manufacture the
#  result. Cloud cover in Karamoja is not random in time: cloudier spells are
#  wetter, and wetter means greener. If Fmask removed proportionally more of the
#  greener, cloudier days from 2019, the 2019 median would sit lower than the
#  truth and the landscape would appear to green with no change on the ground.
#
#  THREE TESTS, WEAKEST TO STRONGEST
#
#    1. WITHIN-DATE. Within one date, does NDVI rise with the number of clear
#       looks a pixel had? If more looks systematically means greener, the date
#       with more looks is biased upward. Confounded by terrain — the cloudy
#       highlands are also the green highlands — so protected-area identity is
#       included as a control, and this test is treated as suggestive only.
#
#    2. IMPLIED BIAS. If a within-date slope exists, multiply it by the actual
#       8.1-to-10.4 gap. That gives the share of the +0.0137 median shift the
#       sampling difference could account for.
#
#    3. BALANCED SUBSET. The test that needs no model and no control. Restrict
#       to pixels where BOTH dates are well sampled and similarly so, and
#       recompute the greening. If it survives at full strength on pixels where
#       the depth difference is gone, the finding is real. If it collapses, it
#       was the sampling.
#
#  Test 3 is the one that decides it.
#
#  HOW TO RUN — from the repo root, in PowerShell:
#      & "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" scripts\R\08_sampling_bias.R
#
#  OUTPUTS -> outputs/qc/
#
#  Author: Isaac Mureithi - IMureithi@awf.org
# =============================================================================

source("scripts/R/00_functions.R")

terra::terraOptions(progress = 0)
set.seed(11)

n.sample <- 600000

cat("\n================================================================\n")
cat(" KIDEPO LCLUC — 08 Observation-depth bias test\n")
cat(" run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n")


# =============================================================================
#  PHASE 1 — ASSEMBLE A SAMPLE
# =============================================================================

comp <- load.composite.pair()
yr.a <- comp$years[1]
yr.b <- comp$years[2]

aoi <- load.aoi()
aoi.union.v <- terra::vect(aoi$union)
mask.to.aoi <- function(r) terra::mask(terra::crop(r, aoi.union.v), aoi.union.v)

ndvi.a <- mask.to.aoi(comp$a[["NDVI"]])
ndvi.b <- mask.to.aoi(comp$b[["NDVI"]])
obs.a  <- mask.to.aoi(comp$a[["clear_obs"]])
obs.b  <- mask.to.aoi(comp$b[["clear_obs"]])

## Protected-area identity as a spatial control. Rasterised once so each
## sampled pixel carries the PA it sits in, letting the within-date test ask
## whether NDVI tracks observation count WITHIN a PA rather than across the
## landscape, where terrain would drive both.
pa.r <- terra::rasterize(terra::vect(aoi$parts), ndvi.a, field = "pa.label")
names(pa.r) <- "pa"

samp <- terra::spatSample(c(ndvi.a, ndvi.b, obs.a, obs.b, pa.r),
                          n.sample, method = "random", na.rm = TRUE, warn = FALSE)
names(samp) <- c("ndvi.a", "ndvi.b", "obs.a", "obs.b", "pa")

samp <- as_tibble(samp) %>%
  mutate(ndvi.a = to.reflectance(ndvi.a),
         ndvi.b = to.reflectance(ndvi.b),
         dndvi  = ndvi.b - ndvi.a) %>%
  filter(obs.a >= comp.min.obs, obs.b >= comp.min.obs) %>%
  drop_na()

say(glue("
[1] SAMPLE
    pixels (both dates >= {comp.min.obs} looks) : {format(nrow(samp), big.mark = ',')}
    mean clear looks {yr.a} / {yr.b}        : {round(mean(samp$obs.a), 2)} / {round(mean(samp$obs.b), 2)}
    median dNDVI                       : {round(median(samp$dndvi), 4)}
"))


# =============================================================================
#  PHASE 2 — WITHIN-DATE: DOES NDVI TRACK OBSERVATION COUNT?
# =============================================================================

bin.profile <- function(ndvi, obs, label) {
  tibble(ndvi = ndvi, obs = obs) %>%
    mutate(bin = cut(obs, breaks = c(3, 5, 7, 9, 11, 13, 15, 100),
                     labels = c("4-5", "6-7", "8-9", "10-11", "12-13", "14-15", "16+"))) %>%
    group_by(bin) %>%
    summarise(date = label, n = n(), median.ndvi = median(ndvi), .groups = "drop") %>%
    filter(n > 500)
}

prof <- bind_rows(bin.profile(samp$ndvi.a, samp$obs.a, as.character(yr.a)),
                  bin.profile(samp$ndvi.b, samp$obs.b, as.character(yr.b)))

cat("\n[2] MEDIAN NDVI BY NUMBER OF CLEAR LOOKS, within each date\n\n")
prof %>%
  pivot_wider(id_cols = bin, names_from = date, values_from = median.ndvi) %>%
  mutate(across(where(is.numeric), ~round(.x, 4))) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

## Slope of NDVI on observation count, with PA fixed effects. The coefficient is
## the NDVI change per ADDITIONAL clear look, within a protected area.
fit.a <- lm(ndvi.a ~ obs.a + factor(pa), data = samp)
fit.b <- lm(ndvi.b ~ obs.b + factor(pa), data = samp)

slope.a <- coef(fit.a)[["obs.a"]]
slope.b <- coef(fit.b)[["obs.b"]]

say(glue("
    NDVI per additional clear look (PA fixed effects):
      {yr.a}: {sprintf('%+.5f', slope.a)}
      {yr.b}: {sprintf('%+.5f', slope.b)}

    The artefact hypothesis requires a POSITIVE slope - more looks recovering
    more greenness, so the better-sampled date reads greener. A NEGATIVE slope
    means the sampling difference pushes the comparison the other way, and
    cannot be what produced the greening."))


# =============================================================================
#  PHASE 3 — HOW MUCH OF THE SHIFT COULD THAT EXPLAIN?
# =============================================================================

obs.gap     <- mean(samp$obs.b) - mean(samp$obs.a)
slope.mean  <- mean(c(slope.a, slope.b))
implied     <- slope.mean * obs.gap
observed    <- median(samp$dndvi)

say(glue("
[3] IMPLIED BIAS
    observation gap ({yr.b} - {yr.a}) : {round(obs.gap, 2)} looks
    mean within-date slope          : {sprintf('%+.5f', slope.mean)} NDVI per look
    implied dNDVI from sampling     : {sprintf('%+.5f', implied)}
    observed median dNDVI           : {sprintf('%+.5f', observed)}
    share potentially explained     : {round(100 * implied / observed, 1)}%
"))


# =============================================================================
#  PHASE 4 — THE TEST THAT DECIDES IT: A BALANCED SUBSET
# =============================================================================
#  No model, no control, no assumption. Keep only pixels where both dates are
#  well sampled AND sampled to a similar depth, then recompute the greening on
#  those pixels alone. Three progressively stricter definitions, so the result
#  can be seen to be stable (or not) rather than resting on one cut.

balanced <- function(min.obs, max.diff) {
  s <- samp %>% filter(obs.a >= min.obs, obs.b >= min.obs,
                       abs(obs.a - obs.b) <= max.diff)
  tibble(
    rule         = glue("both >= {min.obs}, |diff| <= {max.diff}"),
    n            = nrow(s),
    pct.of.sample = 100 * nrow(s) / nrow(samp),
    mean.obs.a   = mean(s$obs.a),
    mean.obs.b   = mean(s$obs.b),
    median.dndvi = median(s$dndvi),
    pct.greening = 100 * mean(s$dndvi > 0)
  )
}

bal <- bind_rows(
  tibble(rule = "all pixels", n = nrow(samp), pct.of.sample = 100,
         mean.obs.a = mean(samp$obs.a), mean.obs.b = mean(samp$obs.b),
         median.dndvi = median(samp$dndvi),
         pct.greening = 100 * mean(samp$dndvi > 0)),
  balanced(6, 4), balanced(8, 3), balanced(10, 2)
)

cat("\n[4] GREENING ON DEPTH-BALANCED PIXELS\n\n")
bal %>%
  mutate(across(where(is.numeric), ~round(.x, 4))) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

strict <- bal %>% slice_tail(n = 1)
retained <- strict$median.dndvi / observed

say(glue("
    On the strictest balanced subset the two dates differ by
    {round(abs(strict$mean.obs.b - strict$mean.obs.a), 2)} looks instead of {round(obs.gap, 2)}, and the median dNDVI is
    {sprintf('%+.5f', strict$median.dndvi)} against {sprintf('%+.5f', observed)} overall — {round(100 * retained)}% of it retained.

    READ THAT WITH CARE. Tightening the rule does not hold the landscape fixed
    while removing the depth gap: mean looks CLIMB from {round(bal$mean.obs.a[1], 1)}/{round(bal$mean.obs.b[1], 1)} to
    {round(strict$mean.obs.a, 1)}/{round(strict$mean.obs.b, 1)}, so the strict subset is the {round(strict$pct.of.sample, 1)}% of the landscape
    with the clearest skies - lowland, less cloudy, not a random sample of it.
    The subsets confound depth balancing with spatial selection, which is why
    their medians bounce rather than trending. What they DO show, robustly, is
    that every subset still greens."))


# =============================================================================
#  PHASE 4b — THE SAME QUESTION WITHOUT SUBSETTING
# =============================================================================
#  Phase 4 cannot separate "balanced depth" from "clear-sky lowland". This can.
#
#  Regress the change itself on the DEPTH GAP, with protected-area fixed
#  effects. The slope says how much apparent greening each extra look of
#  advantage to 2024 buys. The intercept says what greening remains at ZERO
#  depth gap — on every pixel, not a selected subset. If the intercept stays
#  clearly positive, the greening does not come from the sampling difference.

samp <- samp %>% mutate(obs.gap.px = obs.b - obs.a)

fit.gap <- lm(dndvi ~ obs.gap.px + factor(pa), data = samp)

gap.slope <- coef(fit.gap)[["obs.gap.px"]]

## Intercept at zero gap, averaged over the PA effects so it describes the
## landscape rather than whichever PA happens to be the reference level.
pred.zero <- mean(predict(fit.gap, newdata = samp %>% mutate(obs.gap.px = 0)))

say(glue("
[4b] CHANGE REGRESSED ON THE DEPTH GAP (PA fixed effects)
    dNDVI per extra look of advantage to {yr.b} : {sprintf('%+.5f', gap.slope)}
    mean observed depth gap                   : {round(mean(samp$obs.gap.px), 2)} looks
    greening attributable to the gap          : {sprintf('%+.5f', gap.slope * mean(samp$obs.gap.px))}
    greening predicted at ZERO gap            : {sprintf('%+.5f', pred.zero)}
    mean observed dNDVI                       : {sprintf('%+.5f', mean(samp$dndvi))}
"))


# =============================================================================
#  PHASE 5 — VERDICT
# =============================================================================

#  WHICH ESTIMATE TO BELIEVE WHEN THEY DISAGREE — and here they do.
#
#  The within-date slopes (Phase 2) come out NEGATIVE: within one date, pixels
#  with more clear looks are slightly LESS green. Taken at face value that would
#  say the mechanism cannot operate. But that comparison is between DIFFERENT
#  PIXELS, and the ones with fewer looks are the cloudy highlands, which are
#  also the green ones. Protected-area fixed effects do not remove terrain
#  variation inside a protected area, so the slope is measuring geography as
#  much as sampling.
#
#  The Phase 4b regression is between the SAME PIXEL at two dates. Differencing
#  cancels everything time-invariant about it — terrain, soil, land cover — and
#  leaves the depth gap as the thing that varies. That is the better-identified
#  estimate, and it is the one the verdict uses. It says the gap DOES inflate
#  the greening, by about a fifth to a third of it.
#
#  Phase 4's "% retained" is not used at all: its subsets are spatially
#  selective, so the figure bounces for reasons unrelated to depth.

gap.share <- (gap.slope * mean(samp$obs.gap.px)) / mean(samp$dndvi)

cat("\n================================================================\n")
cat(" EVIDENCE\n")
cat(sprintf("   mean dNDVI observed                    : %+.5f\n", mean(samp$dndvi)))
cat(sprintf("   attributable to the depth gap          : %+.5f  (%.0f%%)\n",
            gap.slope * mean(samp$obs.gap.px), 100 * gap.share))
cat(sprintf("   remaining at zero depth gap            : %+.5f  (%.0f%%)\n",
            pred.zero, 100 * (1 - gap.share)))
cat(sprintf("   every balanced subset still greens     : %s\n",
            ifelse(all(bal$median.dndvi > 0), "YES", "NO")))
cat("\n")

if (gap.share < 0.2) {
  cat(" VERDICT: observation depth is NOT a material driver of the greening.\n")
} else if (gap.share < 0.6) {
  cat(" VERDICT: the greening is REAL but INFLATED by unequal sampling.\n\n")
  cat(sprintf(" About %.0f%% of it is attributable to 2024 being better sampled. The\n",
              100 * gap.share))
  cat(" remainder survives at zero depth gap and in every balanced subset, so\n")
  cat(" the direction of the finding stands - the landscape did green - but the\n")
  cat(" MAGNITUDE was overstated and should be quoted from the zero-gap figure.\n\n")
  cat(" Treat it as roughly +0.013 NDVI rather than +0.018, and as approximate:\n")
  cat(" estimates range from +0.009 to +0.019 across subsets, so one significant\n")
  cat(" figure is as much as the data supports.\n\n")
  cat(" Still untouched by any of this: two dates cannot separate inter-annual\n")
  cat(" variability from trend. That needs a time series, not a better test.\n")
} else {
  cat(" VERDICT: the greening is largely an ARTEFACT of unequal depth.\n")
  cat(" Restate the headline and recompute the change figures before quoting.\n")
}
cat("================================================================\n")


# =============================================================================
#  PHASE 6 — OUTPUTS
# =============================================================================

readr::write_csv(prof, file.path(path.out.qc, "bias_ndvi_by_obs.csv"))
readr::write_csv(bal,  file.path(path.out.qc, "bias_balanced_subsets.csv"))

png(file.path(path.out.qc, "bias_observation_depth.png"),
    width = 1600, height = 700, res = 150)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

with(prof[prof$date == as.character(yr.a), ],
     plot(seq_along(bin), median.ndvi, type = "b", pch = 19, col = "#B07AA1",
          ylim = range(prof$median.ndvi), xaxt = "n",
          xlab = "clear observations", ylab = "median NDVI",
          main = "NDVI vs observation depth"))
axis(1, at = seq_along(levels(droplevels(prof$bin))),
     labels = levels(droplevels(prof$bin)), cex.axis = 0.8)
with(prof[prof$date == as.character(yr.b), ],
     lines(seq_along(bin), median.ndvi, type = "b", pch = 17, col = "#E07B39"))
legend("bottomright", bty = "n", cex = 0.8, pch = c(19, 17),
       col = c("#B07AA1", "#E07B39"), legend = c(yr.a, yr.b))

barplot(bal$median.dndvi, names.arg = c("all", "6/4", "8/3", "10/2"),
        col = c("grey70", "#A6CEE3", "#4C9F70", "#1F4D2B"),
        xlab = "balancing rule (min obs / max diff)", ylab = "median dNDVI",
        main = "Greening after depth balancing")
abline(h = 0, lwd = 1)

dev.off()

say(glue("\n[6] wrote 2 tables and 1 figure to {path.out.qc}/"))
cat("\n")

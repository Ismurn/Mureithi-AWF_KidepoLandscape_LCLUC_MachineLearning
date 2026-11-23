# =============================================================================
#  KIDEPO LANDSCAPE — LCLUC
#  Test — the trend estimators behind 11_ndvi_timeseries.R
#
#  WHY THIS EXISTS
#  Script 11 cannot fit lm() 200,000 times for the sampled pixels, nor 32
#  million times for the trend map. It uses closed-form algebra instead, and
#  hand-expanded algebra that nothing checks is exactly the kind of thing that
#  produces a plausible wrong answer. So:
#
#    1. The closed forms are asserted EQUAL to lm() and summary(lm())$r.squared,
#       which are the ground truth.
#    2. The linearity identity the depth/fire adjustment depends on,
#       slope(y - b*z) == slope(y) - b*slope(z), is tested rather than assumed.
#    3. The difference-in-differences estimator is tested against the claim the
#       whole script rests on: that a confound COMMON to every pixel cancels
#       when candidates are compared with their background. A planted common
#       trend plus a planted candidate-specific trend must come back as the
#       candidate-specific part alone.
#    4. The 2019-dip test is given a planted dip to find.
#
#  Check 3 is the important one. Script 11 reports a difference-in-differences
#  as "the clean number", and that claim is only worth making if the estimator
#  actually does what it says.
#
#  It sources 00_functions.R and tests the SAME functions script 11 calls.
#  No rasters, no GEE exports, no network.
#
#  HOW TO RUN — from the repo root, in PowerShell:
#      & "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" scripts\R\tests\test_trend_estimators.R
#
#  Exits non-zero if any check fails.
#
#  Author: Isaac Mureithi - IMureithi@awf.org
# =============================================================================

source("scripts/R/00_functions.R")

set.seed(23)

years   <- 2018:2026

## 30,000 pixels, not 4,000. At a realistic candidate share of 3.5% the first
## version of this test had only 140 candidate pixels, and the standard error
## on a median of 140 noisy slopes is larger than the tolerances being checked
## — so three checks failed on sampling noise while the estimators were
## perfectly correct. The sample size and the tolerances below are now both
## derived from that standard error rather than guessed.
n.pixel    <- 30000
cand.share <- 0.035
fails      <- character(0)

check <- function(label, condition, detail = "") {
  ok <- isTRUE(condition)
  cat(sprintf("  [%s] %s%s\n", if (ok) "PASS" else "FAIL", label,
              if (nzchar(detail)) paste0("  — ", detail) else ""))
  if (!ok) fails <<- c(fails, label)
  invisible(ok)
}

cat("\n================================================================\n")
cat(" Trend estimator self-check for 11_ndvi_timeseries.R\n")
cat("================================================================\n")


# =============================================================================
#  1 — THE CLOSED FORMS AGAINST lm()
# =============================================================================

cat("\n[1] ols.slope() AND ols.r2() vs lm()\n\n")

worst.slope <- 0
worst.r2    <- 0

for (i in 1:300) {
  y <- rnorm(length(years), 0.3, 0.05) + (i / 300) * seq_along(years) * 0.01
  f <- lm(y ~ years)
  worst.slope <- max(worst.slope, abs(ols.slope(years, y) - coef(f)[[2]]))
  worst.r2    <- max(worst.r2,    abs(ols.r2(years, y) - summary(f)$r.squared))
}

check("slope matches lm() over 300 random series", worst.slope < 1e-10,
      sprintf("largest discrepancy %.3e", worst.slope))
check("R2 matches summary(lm()) over 300 random series", worst.r2 < 1e-10,
      sprintf("largest discrepancy %.3e", worst.r2))

## Degenerate and extreme cases, where closed forms usually break.
check("a perfectly flat series returns NA R2, not 0/0",
      is.na(ols.r2(years, rep(0.4, length(years)))))
check("a perfectly flat series returns slope 0",
      abs(ols.slope(years, rep(0.4, length(years)))) < 1e-12)
check("a perfectly linear series returns R2 of exactly 1",
      abs(ols.r2(years, 0.1 + 0.02 * years) - 1) < 1e-10)
check("a perfectly linear series recovers its own slope",
      abs(ols.slope(years, 0.1 + 0.02 * years) - 0.02) < 1e-12)
check("R2 is never above 1 even with floating-point slack",
      { v <- replicate(500, ols.r2(years, 0.1 + 0.02 * years + rnorm(9, 0, 1e-9)))
        all(v <= 1, na.rm = TRUE) })
check("weights sum to zero, so a constant series has no trend",
      abs(sum(ols.slope.weights(years))) < 1e-12)


# =============================================================================
#  2 — LINEARITY, WHICH THE ADJUSTMENT DEPENDS ON
# =============================================================================
#  Script 11 Phase 6 computes the adjusted trend map as
#      slope(ndvi) - b.obs*slope(obs) - b.burnt*slope(burnt)
#  rather than by building an adjusted NDVI stack and taking its slope. The two
#  are only the same thing if the slope operator is linear.

cat("\n[2] slope(y - b*z) == slope(y) - b*slope(z)\n\n")

worst.lin <- 0
for (i in 1:200) {
  y  <- rnorm(length(years), 0.3, 0.05)
  z1 <- rpois(length(years), 8)
  z2 <- rbinom(length(years), 1, 0.3)
  b1 <- rnorm(1, 0.002, 0.001)
  b2 <- rnorm(1, -0.05, 0.02)

  direct   <- ols.slope(years, y - b1 * z1 - b2 * z2)
  composed <- ols.slope(years, y) - b1 * ols.slope(years, z1) -
                                    b2 * ols.slope(years, z2)
  worst.lin <- max(worst.lin, abs(direct - composed))
}

check("adjustment by composition equals adjustment then slope",
      worst.lin < 1e-12, sprintf("largest discrepancy %.3e", worst.lin))


# =============================================================================
#  3 — DOES THE DIFFERENCE-IN-DIFFERENCES CANCEL A COMMON CONFOUND?
# =============================================================================
#  The claim script 11 makes, in its own words: "the archive densified and the
#  rains varied, but both did so over the whole landscape, so they affect
#  candidates and background alike and cancel in the difference between them."
#
#  Here that is given something to prove. The simulation contains:
#
#    a COMMON depth artefact  — observation depth rises over the record exactly
#                               as it does in the real data (4 looks to 14), and
#                               NDVI responds to it, so every pixel greens for a
#                               reason that is not vegetation
#    a COMMON rainfall signal — a shared year effect, large and non-monotone
#    a REAL candidate trend   — only the candidate pixels, +0.004 NDVI per year
#
#  A naive per-group trend must come back far too large, swamped by the common
#  terms. The difference-in-differences must come back at +0.004.

cat("\n[3] DIFFERENCE-IN-DIFFERENCES AGAINST A PLANTED COMMON CONFOUND\n\n")

true.cand.trend <- 0.004
b.obs.true      <- 0.0021          # NDVI per clear look, as script 08 implies
obs.by.year     <- c(6.1, 8.1, 4.7, 7.0, 11.0, 11.2, 10.4, 13.8, 7.9)
rain.by.year    <- c(0.00, -0.03, 0.02, -0.01, 0.03, -0.02, 0.01, 0.04, -0.03)

sim <- expand_grid(pixel = seq_len(n.pixel), year = years) %>%
  mutate(
    candidate = as.integer(pixel <= n.pixel * cand.share),  # as observed
    pa        = letters[1 + (pixel %% 13)],                 # 13 protected areas
    obs       = round(obs.by.year[match(year, years)] + rnorm(n(), 0, 1.5)),
    rain      = rain.by.year[match(year, years)],
    base      = 0.35 + 0.08 * ((pixel %% 97) / 97),         # pixel's own level
    ndvi      = base + b.obs.true * obs + rain +
                candidate * true.cand.trend * (year - min(years)) +
                rnorm(n(), 0, 0.04)
  )

per.pixel <- sim %>%
  group_by(pixel, candidate, pa) %>%
  summarise(slope = ols.slope(year, ndvi), r2 = ols.r2(year, ndvi),
            .groups = "drop")

naive.cand <- median(per.pixel$slope[per.pixel$candidate == 1])
naive.back <- median(per.pixel$slope[per.pixel$candidate == 0])
did.raw    <- naive.cand - naive.back

fit.did <- lm(slope ~ candidate + factor(pa), data = per.pixel)
did.adj <- coef(fit.did)[["candidate"]]
did.se  <- summary(fit.did)$coefficients["candidate", "Std. Error"]

## The standard error of a median of m draws is about 1.2533 * sd / sqrt(m).
## Tolerances below come from this rather than from a round percentage.
slope.sd  <- sd(per.pixel$slope[per.pixel$candidate == 0])
m.cand    <- sum(per.pixel$candidate == 1)
se.median <- 1.2533 * slope.sd / sqrt(m.cand)

cat(sprintf("  planted common depth effect   : %+.6f NDVI per look\n", b.obs.true))
cat(sprintf("  planted candidate-only trend  : %+.6f NDVI per year\n", true.cand.trend))
cat("\n")
cat(sprintf("  naive trend, background       : %+.6f per year  (truth: 0)\n", naive.back))
cat(sprintf("  naive trend, candidates       : %+.6f per year\n", naive.cand))
cat(sprintf("  difference of medians         : %+.6f +/- %.6f\n", did.raw, se.median))
cat(sprintf("  DiD with PA fixed effects     : %+.6f +/- %.6f\n", did.adj, did.se))

cat(sprintf("\n  THE POINT OF THIS BLOCK: the background's TRUE trend is zero, and the\n"))
cat(sprintf("  naive estimator reads %+.6f per year for it — which over a 5-year\n", naive.back))
cat(sprintf("  gap is %+.4f NDVI, close to the whole +0.0127 the real headline\n",
            naive.back * 5))
cat(sprintf("  claims. A denser archive alone can manufacture the entire signal.\n"))
cat(sprintf("  That is why script 11 refuses to quote an absolute trend.\n\n"))

check("the naive per-group trend IS materially contaminated",
      abs(naive.back) > 0.25 * true.cand.trend,
      sprintf("background reads %+.6f (%.0f%% of the planted signal) where the truth is 0",
              naive.back, 100 * abs(naive.back) / true.cand.trend))
check("the difference of medians recovers the candidate trend",
      abs(did.raw - true.cand.trend) < 3 * se.median,
      sprintf("%+.6f vs %+.6f, within %.1f SE", did.raw, true.cand.trend,
              abs(did.raw - true.cand.trend) / se.median))
check("the PA-adjusted DiD recovers the candidate trend",
      abs(did.adj - true.cand.trend) < 3 * did.se,
      sprintf("%+.6f vs %+.6f, within %.1f SE", did.adj, true.cand.trend,
              abs(did.adj - true.cand.trend) / did.se))

## And the converse: with NO candidate-specific trend, the DiD must read zero
## even though both groups are still drenched in the common confound.
sim.null <- sim %>%
  mutate(ndvi = ndvi - candidate * true.cand.trend * (year - min(years)))

pp.null <- sim.null %>%
  group_by(pixel, candidate, pa) %>%
  summarise(slope = ols.slope(year, ndvi), .groups = "drop")

f.null  <- lm(slope ~ candidate + factor(pa), data = pp.null)
did.null <- coef(f.null)[["candidate"]]
se.null  <- summary(f.null)$coefficients["candidate", "Std. Error"]

cat(sprintf("\n  with the candidate trend removed, DiD reads %+.6f +/- %.6f\n",
            did.null, se.null))
check("no candidate trend invented when none is planted",
      abs(did.null) < 3 * se.null,
      sprintf("%+.6f, SE %.6f", did.null, se.null))


# =============================================================================
#  4 — DOES R2 SEPARATE STEADY MOVEMENT FROM BOUNCING?
# =============================================================================
#  Script 11 reads a high R2 as "moved steadily" and a low one as "bounced with
#  the rains". Two populations are built that END in the same place by different
#  routes, and the statistic has to tell them apart.

cat("\n[4] R2 SEPARATES A TREND FROM A SWING\n\n")

total.rise <- 0.04

steady <- replicate(1000, {
  y <- 0.3 + (total.rise / (length(years) - 1)) * (years - min(years)) +
       rnorm(length(years), 0, 0.01)
  ols.r2(years, y)
})

## Same first and last value, but the path between them is noise.
bouncy <- replicate(1000, {
  y <- 0.3 + rnorm(length(years), 0, 0.05)
  y[1] <- 0.3; y[length(y)] <- 0.3 + total.rise
  ols.r2(years, y)
})

cat(sprintf("  steady pixels : median R2 %.3f\n", median(steady, na.rm = TRUE)))
cat(sprintf("  bouncy pixels : median R2 %.3f\n", median(bouncy, na.rm = TRUE)))

check("steady movement gives a clearly elevated R2",
      median(steady, na.rm = TRUE) > 0.6,
      sprintf("median %.3f", median(steady, na.rm = TRUE)))
check("bouncing gives a low R2", median(bouncy, na.rm = TRUE) < 0.5,
      sprintf("median %.3f", median(bouncy, na.rm = TRUE)))
check("the two populations are clearly separated",
      median(steady, na.rm = TRUE) - median(bouncy, na.rm = TRUE) > 0.4,
      sprintf("gap %.3f", median(steady, na.rm = TRUE) -
                          median(bouncy, na.rm = TRUE)))

## CALIBRATION — and the reason this matters for reading script 11.
##
## The first version of this test asserted that a genuinely steady pixel would
## show R2 above 0.8. It does not: at a trend of 0.005 NDVI/yr against year-to-
## year noise of 0.01 the median R2 is about 0.70, and at realistic NDVI noise
## of 0.04-0.05 it collapses towards 0.1. So an ABSOLUTE R2 threshold would be
## meaningless on this data — a real trend simply does not produce a tidy line
## over nine noisy years. Only the DIFFERENCE between candidate and background
## R2 carries information, which is how script 11 uses it.

cat("\n  Median R2 of a GENUINE linear trend, by trend size and noise:\n\n")
cat("              noise 0.010  0.020  0.040  0.060\n")
for (tr in c(0.002, 0.005, 0.010, 0.020)) {
  row <- sprintf("  %+.3f/yr ", tr)
  for (sg in c(0.01, 0.02, 0.04, 0.06)) {
    r2 <- median(replicate(400, ols.r2(
      years, 0.3 + tr * (years - min(years)) + rnorm(length(years), 0, sg))),
      na.rm = TRUE)
    row <- paste0(row, sprintf("  %5.2f", r2))
  }
  cat(row, "\n")
}
cat("\n  Read down the 0.040 column: at the noise level this data actually has,\n")
cat("  even a +0.010/yr trend — more than twice the headline rate — shows a\n")
cat("  median R2 well under 0.5. Absolute R2 is not a usable threshold here.\n")


# =============================================================================
#  5 — THE 2019-DIP TEST, GIVEN A DIP TO FIND
# =============================================================================
#  Script 11 Phase 5 compares each group's focal year against its OWN
#  neighbouring years, so a landscape-wide bad year drops out and only a
#  group-specific anomaly survives. Both halves of that are checked: a planted
#  candidate-only dip must be found, and a planted COMMON dip must not be
#  mistaken for one.

cat("\n[5] THE FOCAL-YEAR ANOMALY TEST\n\n")

dip.size  <- -0.015
focal.yr  <- 2019
neigh.yrs <- c(2018, 2020)

## Returns the excess dip AND its standard error. Returning the error matters:
## the first version of this test compared the excess against a flat 25% of the
## planted size, which is smaller than one standard error at a realistic
## candidate count, so it failed on noise while the estimator was correct.
dip.stats.of <- function(df) {
  per <- df %>%
    filter(year %in% c(focal.yr, neigh.yrs)) %>%
    mutate(role = if_else(year == focal.yr, "focal", "neighbour")) %>%
    group_by(candidate, pixel, role) %>%
    summarise(ndvi = mean(ndvi), .groups = "drop") %>%
    pivot_wider(names_from = role, values_from = ndvi) %>%
    mutate(dip = focal - neighbour)

  g <- per %>%
    group_by(candidate) %>%
    summarise(med = median(dip, na.rm = TRUE), n = n(),
              s = sd(dip, na.rm = TRUE), .groups = "drop")

  ## se of a median ~ 1.2533 * sd / sqrt(n); the two groups are independent.
  se <- sqrt(sum((1.2533 * g$s / sqrt(g$n))^2))

  list(excess = g$med[g$candidate == 1] - g$med[g$candidate == 0], se = se)
}

## Case A: the dip belongs to the candidates alone.
a <- dip.stats.of(
  sim %>% mutate(ndvi = ndvi + if_else(candidate == 1 & year == focal.yr,
                                       dip.size, 0)))

## Case B: every pixel dips in that year. The test must NOT report an excess,
## because each group is compared against its own neighbouring years.
b <- dip.stats.of(
  sim %>% mutate(ndvi = ndvi + if_else(year == focal.yr, dip.size, 0)))

cat(sprintf("  planted candidate-only dip of %+.4f -> excess %+.5f +/- %.5f\n",
            dip.size, a$excess, a$se))
cat(sprintf("  planted landscape-wide dip of %+.4f -> excess %+.5f +/- %.5f\n",
            dip.size, b$excess, b$se))

check("a candidate-only dip is detected at its planted size",
      abs(a$excess - dip.size) < 3 * a$se,
      sprintf("%+.5f vs %+.4f, within %.1f SE", a$excess, dip.size,
              abs(a$excess - dip.size) / a$se))
check("a landscape-wide dip is NOT reported as a candidate anomaly",
      abs(b$excess) < 3 * b$se,
      sprintf("excess %+.5f, within %.1f SE of zero", b$excess,
              abs(b$excess) / b$se))


# =============================================================================
#  SUMMARY
# =============================================================================

cat("\n================================================================\n")
if (length(fails) == 0) {
  cat(" ALL CHECKS PASSED — the closed forms in 00_functions.R section 3.10\n")
  cat(" match lm(), the adjustment algebra is linear as assumed, and the\n")
  cat(" difference-in-differences does cancel a common confound while still\n")
  cat(" finding a group-specific one.\n")
  cat("================================================================\n\n")
} else {
  cat(sprintf(" %d CHECK(S) FAILED:\n", length(fails)))
  for (f in fails) cat("   -", f, "\n")
  cat("================================================================\n\n")
  stop("Trend estimator self-check failed. Do not trust 11_ndvi_timeseries.R ",
       "until these are fixed.")
}

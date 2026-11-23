# =============================================================================
#  KIDEPO LANDSCAPE — LCLUC
#  Test — the solar geometry and illumination arithmetic behind 10_terrain_test.R
#
#  WHY THIS EXISTS
#  Script 10 decides whether 277.5 km2 of candidate change is an artefact of
#  terrain illumination. That decision rests on a chain of trigonometry — solar
#  declination, solar position, the local incidence angle on an inclined
#  surface — where a sign error or a degrees-for-radians slip would produce
#  numbers that look entirely reasonable and are entirely wrong. Putting the
#  sun in the west instead of the east, for instance, would invert every
#  statement about which slopes are lit without changing the magnitude of
#  anything.
#
#  So the arithmetic is checked against values that are known independently,
#  and the estimator is checked against data where the answer was planted.
#
#  It sources 00_functions.R and tests the SAME functions script 10 calls. A
#  test carrying its own copy of the formulae would drift away from the
#  analysis silently and then certify the wrong code.
#
#  It needs no rasters, no GEE exports and no network. It runs in a couple of
#  minutes and should be run after any change to section 3.9 of 00_functions.R.
#
#  HOW TO RUN — from the repo root, in PowerShell:
#      & "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" scripts\R\tests\test_solar_geometry.R
#
#  Exits non-zero if any check fails, so it can go in a hook or CI later.
#
#  Author: Isaac Mureithi - IMureithi@awf.org
# =============================================================================

source("scripts/R/00_functions.R")

set.seed(99)

lat.kidepo <- 3.9     # roughly the centre of the PA network
n.rep      <- 30      # replications for the estimator check
n.sim      <- 40000   # pixels per replication

fails <- character(0)

check <- function(label, condition, detail = "") {
  ok <- isTRUE(condition)
  cat(sprintf("  [%s] %s%s\n", if (ok) "PASS" else "FAIL", label,
              if (nzchar(detail)) paste0("  — ", detail) else ""))
  if (!ok) fails <<- c(fails, label)
  invisible(ok)
}

cat("\n================================================================\n")
cat(" Solar geometry self-check for 10_terrain_test.R\n")
cat("================================================================\n")


# =============================================================================
#  1 — DECLINATION AGAINST THE ASTRONOMICAL CALENDAR
# =============================================================================
#  Independently known: the declination is 0 at the equinoxes and +/-23.44 at
#  the solstices. Nothing about the project is involved in those values.

cat("\n[1] SOLAR DECLINATION\n\n")

dec <- sapply(c(mar = 80, jun = 172, sep = 266, dec = 355), solar.declination)
print(round(dec, 2))

check("equinox declinations near zero",
      max(abs(dec[c("mar", "sep")])) < 1.5,
      sprintf("max |dec| = %.2f deg", max(abs(dec[c("mar", "sep")]))))
check("June solstice near +23.44", abs(dec[["jun"]] - 23.44) < 0.5,
      sprintf("%.2f deg", dec[["jun"]]))
check("December solstice near -23.44", abs(dec[["dec"]] + 23.44) < 0.5,
      sprintf("%.2f deg", dec[["dec"]]))


# =============================================================================
#  2 — SOLAR POSITION OVER THE COMPOSITE WINDOW
# =============================================================================
#  The dry_late window is 1 Feb to 31 Mar. A mid-morning overpass at 4 degrees
#  north should put the sun 25-35 degrees from vertical and somewhere east of
#  south — if it comes out in the west, the azimuth convention is wrong.

cat("\n[2] SOLAR POSITION, dry_late window at Kidepo\n\n")

doys <- c(32, 46, 60, 75, 91)
pos  <- solar.position(doys, lat.kidepo)
print(data.frame(doy = doys,
                 date = format(as.Date(doys - 1, origin = "2019-01-01"), "%d %b"),
                 zenith = round(pos$zenith, 2),
                 azimuth = round(pos$azimuth, 2)))

check("zenith plausible for this latitude and season",
      all(pos$zenith > 20 & pos$zenith < 40),
      sprintf("%.1f-%.1f deg", min(pos$zenith), max(pos$zenith)))
check("sun is EAST of south at a morning overpass",
      all(pos$azimuth > 60 & pos$azimuth < 180),
      sprintf("azimuth %.1f-%.1f deg", min(pos$azimuth), max(pos$azimuth)))
check("zenith falls as the window advances towards the equinox",
      all(diff(pos$zenith) < 0))

## An afternoon overpass must mirror the morning one about due south. If this
## fails, the hour angle sign is wrong and the lit slopes are the wrong ones.
am <- solar.position(60, lat.kidepo, solar.hour = 10)
pm <- solar.position(60, lat.kidepo, solar.hour = 14)
check("morning and afternoon azimuths mirror about south",
      abs((180 - am$azimuth) - (pm$azimuth - 180)) < 1,
      sprintf("am %.1f, pm %.1f", am$azimuth, pm$azimuth))
check("solar noon puts the sun nearly due south at this latitude",
      abs(solar.position(60, lat.kidepo, solar.hour = 12)$azimuth - 180) < 1)


# =============================================================================
#  3 — THE INCIDENCE ANGLE ON A SLOPE
# =============================================================================

cat("\n[3] cos(i) ON AN INCLINED SURFACE\n\n")

s60 <- solar.position(60, lat.kidepo)

ci <- sapply(c(N = 0, E = 90, S = 180, W = 270),
             function(a) cos.illum(20, a, s60$zenith, s60$azimuth))
print(round(ci, 4))

flat <- cos.illum(0, 0, s60$zenith, s60$azimuth)

## The cheapest possible check that the formula is typed right.
check("flat ground reduces to cos(zenith)",
      abs(flat - cos(deg2rad(s60$zenith))) < 1e-12,
      sprintf("%.6f vs %.6f", flat, cos(deg2rad(s60$zenith))))
check("east-facing slope is the brightest at a morning overpass",
      which.max(ci) == which(names(ci) == "E"))
check("west-facing slope is the dimmest", which.min(ci) == which(names(ci) == "W"))
check("flat ground sits between the lit and shaded extremes",
      flat > min(ci) && flat < max(ci))
check("a slope normal to the sun returns 1",
      abs(cos.illum(s60$zenith, s60$azimuth, s60$zenith, s60$azimuth) - 1) < 1e-12)
check("cos(i) never leaves [-1, 1]",
      { v <- cos.illum(rep(seq(0, 80, 5), each = 72), rep(seq(0, 355, 5), 17),
                       s60$zenith, s60$azimuth)
        all(v >= -1 - 1e-12 & v <= 1 + 1e-12) })


# =============================================================================
#  4 — OCTANT BINNING
# =============================================================================
#  Script 10 bins aspect into eight compass octants. North wraps across 360,
#  which is where this kind of code usually breaks.

cat("\n[4] ASPECT OCTANT BINNING\n\n")

oct.levels <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
to.octant  <- function(a) oct.levels[floor(((a + 22.5) %% 360) / 45) + 1]

probe <- c(0, 22, 23, 45, 90, 135, 180, 225, 270, 315, 337, 338, 359, 360)
want  <- c("N", "N", "NE", "NE", "E", "SE", "S", "SW", "W", "NW", "NW",
           "N", "N", "N")
print(data.frame(aspect = probe, octant = to.octant(probe), expected = want))

check("octants correct, north wrap included", identical(to.octant(probe), want))
check("every octant gets an eighth of a uniform circle",
      { tb <- table(to.octant(seq(0, 359.9, 0.1)))
        length(tb) == 8 && max(abs(tb - mean(tb))) <= 1 })


# =============================================================================
#  5 — HOW FAR cos(i) CAN MOVE FOR A GIVEN SHIFT IN OBSERVATION DAY
# =============================================================================
#  This is the quantity the whole terrain test rests on, so it is tabulated
#  rather than asserted: the worst case over aspect that a timing shift can do
#  to the illumination of a slope. The numbers printed here are the ones quoted
#  in the header of 10_terrain_test.R and in section 7.3 of the write-up.

cat("\n[5] CEILING ON delta cos(i) BY TIMING SHIFT AND SLOPE\n\n")

worst.dcosi <- function(shift.days, slope.deg, doy0 = 60, lat = lat.kidepo) {
  a <- solar.position(doy0, lat)
  b <- solar.position(doy0 + shift.days, lat)
  d <- sapply(seq(0, 350, 10), function(asp)
    cos.illum(slope.deg, asp, b$zenith, b$azimuth) -
    cos.illum(slope.deg, asp, a$zenith, a$azimuth))
  d[which.max(abs(d))]
}

shifts <- c(1, 3, 5, 10, 20)
slopes <- c(0, 5, 10, 20, 30)
tabl <- outer(shifts, slopes, Vectorize(worst.dcosi))
dimnames(tabl) <- list(paste0(shifts, " d"), paste0(slopes, " deg"))
print(round(tabl, 5))

check("a one-day shift cannot move cos(i) by more than 0.01 on any slope",
      max(abs(tabl["1 d", ])) < 0.01,
      sprintf("max %.5f", max(abs(tabl["1 d", ]))))
check("the ceiling rises with slope", all(diff(abs(tabl["5 d", ])) > 0))
check("the ceiling rises with the size of the shift",
      all(diff(abs(tabl[, "20 deg"])) > 0))

cat("\n  Read these as a CEILING, not an expectation: they are the worst aspect\n")
cat("  at each slope, and the landscape mean is far smaller.\n")


# =============================================================================
#  6 — IS THE ESTIMATOR UNBIASED?
# =============================================================================
#  Script 10 identifies the illumination effect by regressing dNDVI on delta
#  cos(i) with protected-area fixed effects. Here a known coefficient is
#  planted in synthetic data built from the real helpers, and recovered.
#
#  Two cases, because a test that only checks one is half a test:
#    NULL     greening with no terrain in it at all — the estimator must not
#             manufacture an effect
#    PLANTED  a genuine illumination response of 0.30 — the estimator must
#             find it

cat("\n[6] ESTIMATOR,", n.rep, "replications of", format(n.sim, big.mark = ","),
    "pixels\n\n")

make.sim <- function(n) {
  d <- data.frame(
    lat   = runif(n, 3.3, 4.3),           # the AOI's latitude span
    slope = pmin(pmax(rgamma(n, 2, 0.3), 0), 60),
    asp   = runif(n, 0, 360),
    doy.a = rnorm(n, 58, 6),              # observation timing, 2019
    doy.b = rnorm(n, 61, 6),              # and 2024, three days later
    pa    = sample(letters[1:5], n, TRUE)
  )
  sa <- solar.position(d$doy.a, d$lat)
  sb <- solar.position(d$doy.b, d$lat)
  d$dcosi <- cos.illum(d$slope, d$asp, sb$zenith, sb$azimuth) -
             cos.illum(d$slope, d$asp, sa$zenith, sa$azimuth)
  d
}

b.null <- numeric(n.rep); b.plant <- numeric(n.rep)
se.plant <- numeric(n.rep); mean.dcosi <- numeric(n.rep)

for (i in seq_len(n.rep)) {
  d <- make.sim(n.sim)
  mean.dcosi[i] <- mean(d$dcosi)

  d$y.null  <- 0.0127 + rnorm(n.sim, 0, 0.05)
  d$y.plant <- 0.0127 + 0.30 * d$dcosi + rnorm(n.sim, 0, 0.05)

  b.null[i]  <- coef(lm(y.null  ~ dcosi + factor(pa), data = d))[["dcosi"]]
  f2 <- lm(y.plant ~ dcosi + factor(pa), data = d)
  b.plant[i]  <- coef(f2)[["dcosi"]]
  se.plant[i] <- summary(f2)$coefficients["dcosi", "Std. Error"]
}

t.null  <- mean(b.null)  / (sd(b.null)  / sqrt(n.rep))
t.plant <- (mean(b.plant) - 0.30) / (sd(b.plant) / sqrt(n.rep))

cat(sprintf("  mean delta cos(i) in the simulation : %+.5f\n", mean(mean.dcosi)))
cat(sprintf("  NULL    mean estimate %+.4f   (planted 0.0000)  t = %+.2f\n",
            mean(b.null), t.null))
cat(sprintf("  PLANTED mean estimate %+.4f   (planted 0.3000)  t = %+.2f\n",
            mean(b.plant), t.plant))
cat(sprintf("  spread of the planted estimates     : sd %.4f\n", sd(b.plant)))
cat(sprintf("  mean standard error it reported     : %.4f\n", mean(se.plant)))

check("no effect invented where none exists", abs(t.null) < 3,
      sprintf("t = %+.2f", t.null))
check("a planted coefficient of 0.30 is recovered", abs(t.plant) < 3,
      sprintf("recovered %.4f, t = %+.2f", mean(b.plant), t.plant))
check("reported standard errors match the real spread",
      abs(mean(se.plant) - sd(b.plant)) < 0.5 * sd(b.plant),
      sprintf("reported %.4f vs actual %.4f", mean(se.plant), sd(b.plant)))


# =============================================================================
#  7 — WHAT A RECOVERED COEFFICIENT WOULD MEAN
# =============================================================================
#  Not an assertion — the arithmetic that keeps the eventual result in
#  proportion. The greening illumination can explain is the coefficient times
#  the MEAN shift in cos(i), and that mean is small.

cat("\n[7] GREENING EXPLAINED, BY COEFFICIENT\n\n")

obs.greening <- 0.0127            # the depth-adjusted figure from script 08
mean.shift   <- mean(mean.dcosi)  # as simulated, for a three-day timing shift

for (b in c(0.1, 0.3, 0.5, 1.0)) {
  cat(sprintf("  coefficient %.1f  ->  %+.5f NDVI  =  %4.1f%% of a %+.4f greening\n",
              b, b * mean.shift, 100 * b * mean.shift / obs.greening, obs.greening))
}

cat("\n  NDVI's cross-sectional spread across aspects in rugged terrain is\n")
cat("  rarely more than about 0.1 per unit cos(i), so the top rows here are\n")
cat("  not physically available. The ceiling on what illumination can explain\n")
cat("  is single-digit percent unless the two years were observed many days\n")
cat("  apart — which is why the obs_doy bands matter more than the DEM does.\n")


# =============================================================================
#  SUMMARY
# =============================================================================

cat("\n================================================================\n")
if (length(fails) == 0) {
  cat(" ALL CHECKS PASSED — the illumination arithmetic in 00_functions.R\n")
  cat(" section 3.9 behaves, and the estimator in 10_terrain_test.R is\n")
  cat(" unbiased in both directions.\n")
  cat("================================================================\n\n")
} else {
  cat(sprintf(" %d CHECK(S) FAILED:\n", length(fails)))
  for (f in fails) cat("   -", f, "\n")
  cat("================================================================\n\n")
  stop("Solar geometry self-check failed. Do not trust 10_terrain_test.R ",
       "until these are fixed.")
}

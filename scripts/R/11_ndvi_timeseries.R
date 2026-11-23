# =============================================================================
#  KIDEPO LANDSCAPE — LCLUC
#  Script 11 — The annual NDVI series: trend, variability, and candidate
#              trajectories
#
#  WHY
#  Everything measured so far rests on two dates. A +0.013 NDVI difference
#  between 2019 and 2024 cannot be told apart from one drier year against one
#  wetter one. This script uses the annual dry_late series from
#  gee/03_ndvi_timeseries.js to ask what two dates cannot.
#
#  WHAT THIS SCRIPT CAN AND CANNOT DELIVER — read before the numbers.
#
#  It CANNOT give a trustworthy absolute trend, and pretending otherwise would
#  be the single easiest way to mislead with this data. Mean clear looks per
#  pixel go 4.05, 3.15, 6.10, 8.13, 4.75, 7.00, 10.96, 11.25, 10.41, 13.85 over
#  2016-2025. That is the satellite constellation growing, not the weather, and
#  script 08 established that on this data more looks produce higher NDVI. Year
#  and observation depth are therefore almost collinear across the record: a
#  regression asked to separate "greening" from "we looked more often" has
#  nothing in the aggregate data to separate them WITH. Any script that reports
#  a single depth-adjusted trend here is making a claim its data cannot carry.
#
#  It CAN answer the questions that are RELATIVE, and those are the ones that
#  matter. The depth trend and the rainfall history are COMMON to every pixel in
#  this small landscape, so they cancel when candidates are compared with their
#  own background:
#
#    Q1  Do the candidate patches move MONOTONICALLY, as thickening or land-use
#        change would, or do they bounce with the rains?
#    Q2  Is the candidate trend steeper than the background trend? That
#        difference-in-differences is clean — the common depth and rainfall
#        trends subtract out of it.
#    Q3  Was 2019 an anomalous DIP for the candidates specifically? If so, part
#        of the 2019->2024 "change" is a 2019 artefact rather than a 2024 gain.
#    Q4  How large is the 2019->2024 step against ordinary year-to-year spread?
#
#  Q2 and Q3 are the ones that can actually move the conclusion from script 09.
#
#  A SECOND CONFOUND: fire. Cumulative burnt fraction runs from 19.7% (2019) to
#  43.9% (2018). 2019 and 2024 are the two LEAST burnt years on the record,
#  which is why script 02 chose them, and it means the chosen pair sits at one
#  end of the fire distribution. The burnt stack is carried so this can be
#  handled rather than hoped away.
#
#  NEEDS  scripts/gee/03_ndvi_timeseries.js run and downloaded:
#           kidepo_ndvi_annual_dry_late_2016_2026_utm36n_30m.tif   (required)
#           kidepo_obs_annual_dry_late_2016_2026_utm36n_30m.tif    (required)
#           kidepo_burnt_annual_dry_late_2016_2026_utm36n_30m.tif  (required)
#           kidepo_ndvi_annual_aoi.csv                             (optional)
#
#  HOW TO RUN — from the repo root, in PowerShell. NOT in Git Bash: the MSYS
#  mingw64 DLLs shadow GDAL and sf segfaults.
#      & "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" scripts\R\11_ndvi_timeseries.R
#
#  Expect several minutes: the trend map touches 11 bands over ~32 million cells.
#
#  OUTPUTS -> outputs/timeseries/
#
#  Author: Isaac Mureithi - IMureithi@awf.org
# =============================================================================

source("scripts/R/00_functions.R")
pacman::p_load(tidyterra, ggspatial)

terra::terraOptions(progress = 0)
set.seed(13)

path.out.ts <- file.path("outputs", "timeseries")
dir.create(path.out.ts, recursive = TRUE, showWarnings = FALSE)

ts.years    <- 2016:2026    # as exported, ascending; band order depends on this
ts.panel    <- 2018:2026    # the balanced panel — 2016/2017 coverage is too thin
ts.provisional <- 2026      # HLS processing may not have caught up
n.samp      <- 200000       # pixels for the panel regressions

cat("\n================================================================\n")
cat(" KIDEPO LCLUC — 11 Annual NDVI series\n")
cat(" run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n")


# =============================================================================
#  PHASE 0 — LOAD THE THREE STACKS, AND CHECK THEM
# =============================================================================
#  Each stack carries one band per year in ascending order. Band names are
#  assigned POSITIONALLY from ts.years rather than read from the GeoTIFF: a
#  silently reordered export would otherwise attach the wrong year to every
#  value, and nothing downstream would complain.

span <- paste0(min(ts.years), "_", max(ts.years))

ts.files <- c(
  ndvi  = glue("kidepo_ndvi_annual_{comp.window}_{span}_utm36n_30m.tif"),
  obs   = glue("kidepo_obs_annual_{comp.window}_{span}_utm36n_30m.tif"),
  burnt = glue("kidepo_burnt_annual_{comp.window}_{span}_utm36n_30m.tif")
)

missing <- ts.files[!file.exists(file.path(path.gee.out, ts.files))]
if (length(missing)) {
  stop(glue(
    "Missing {length(missing)} of 3 annual stacks in\n  {path.gee.out}\n\n",
    "  {paste(missing, collapse = '\n  ')}\n\n",
    "Run scripts/gee/03_ndvi_timeseries.js, let all four exports finish in ",
    "Tasks, then download them. The CSV is optional; the three GeoTIFFs are ",
    "not."))
}

load.annual <- function(kind, prefix) {
  r <- terra::rast(file.path(path.gee.out, ts.files[[kind]]))
  if (terra::nlyr(r) != length(ts.years)) {
    stop(glue("{ts.files[[kind]]} has {terra::nlyr(r)} bands; expected ",
              "{length(ts.years)}, one per year {span}."))
  }
  names(r) <- paste0(prefix, ts.years)
  terra::NAflag(r) <- comp.nodata
  r
}

ndvi.r  <- load.annual("ndvi",  "ndvi.")
obs.r   <- load.annual("obs",   "obs.")
burnt.r <- load.annual("burnt", "burnt.")

comp <- load.composite.pair()

for (nm in c("ndvi.r", "obs.r", "burnt.r")) {
  if (!isTRUE(terra::compareGeom(comp$a, get(nm), stopOnError = FALSE,
                                 messages = FALSE))) {
    stop(glue("{nm} is NOT on the composite grid. Re-export it from 03 with ",
              "an identical `region`, `crs` and `scale`."))
  }
}

## Crop and mask to the AOI BEFORE the range checks rather than after. These
## are 33 bands over ~32 million cells; computing minmax on the full rectangle
## and then reading it all again for the annual means would be two complete
## passes over about a gigabyte for no gain.
aoi <- load.aoi()
aoi.union.v <- terra::vect(aoi$union)
aoi.parts.v <- terra::vect(aoi$parts)
mask.to.aoi <- function(r) terra::mask(terra::crop(r, aoi.union.v), aoi.union.v)

ndvi.m  <- mask.to.aoi(ndvi.r)
obs.m   <- mask.to.aoi(obs.r)
burnt.m <- mask.to.aoi(burnt.r)

## Range assertions, same discipline as script 10: the encoding comments are
## not evidence. A scale factor applied twice cannot pass these.
mm.ndvi  <- terra::minmax(ndvi.m,  compute = TRUE)
mm.obs   <- terra::minmax(obs.m,   compute = TRUE)
mm.burnt <- terra::minmax(burnt.m, compute = TRUE)

say(glue("
[0] ENCODING CHECKS
    NDVI  stored range : {min(mm.ndvi)} to {max(mm.ndvi)}  (expect roughly -10000..10000)
    NDVI  true range   : {round(min(mm.ndvi) / comp.scale.int, 3)} to {round(max(mm.ndvi) / comp.scale.int, 3)}
    obs   range        : {min(mm.obs)} to {max(mm.obs)} clear looks
    burnt range        : {min(mm.burnt)} to {max(mm.burnt)}  (must be 0 to 1)"))

if (min(mm.ndvi) < -comp.scale.int || max(mm.ndvi) > comp.scale.int) {
  stop("NDVI is outside +/-1 after scaling — the stored scale factor is not ",
       comp.scale.int, ". Check the export encoding before reading anything below.")
}
if (max(mm.burnt) > 1 || min(mm.burnt) < 0) {
  stop("The burnt stack is not 0/1. It should be a binary mask.")
}
if (max(mm.burnt) == 0) {
  stop("The burnt stack is ZERO everywhere. This is the script 01c failure ",
       "repeating: MODIS burn data lost between the reducer and the export. ",
       "Do not proceed — fire is the second confound and it cannot be handled ",
       "without this stack.")
}
say("    All three stacks pass. Encoding confirmed.")


# =============================================================================
#  PHASE 1 — THE ANNUAL SERIES, AND THE TWO CONFOUNDS BESIDE IT
# =============================================================================
#  Printed together on purpose. The NDVI column cannot be read honestly without
#  the two columns next to it.

annual <- tibble(
  year       = ts.years,
  mean.ndvi  = as.numeric(terra::global(ndvi.m,  "mean", na.rm = TRUE)[, 1]) /
                 comp.scale.int,
  mean.obs   = as.numeric(terra::global(obs.m,   "mean", na.rm = TRUE)[, 1]),
  burnt.frac = as.numeric(terra::global(burnt.m, "mean", na.rm = TRUE)[, 1])
) %>%
  mutate(in.panel    = year %in% ts.panel,
         provisional = year %in% ts.provisional)

cat("\n[1] THE ANNUAL SERIES\n\n")
annual %>%
  mutate(mean.ndvi = round(mean.ndvi, 4),
         mean.obs = round(mean.obs, 2),
         burnt.frac = round(burnt.frac, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)

## The cross-check against the GEE table, if it was downloaded. Three code
## paths agreeing is the only reason to believe any of them.
f.csv <- file.path(path.gee.out, "kidepo_ndvi_annual_aoi.csv")
if (file.exists(f.csv)) {
  gee.tbl <- readr::read_csv(f.csv, show_col_types = FALSE)
  cmp.csv <- annual %>%
    select(year, r.ndvi = mean.ndvi, r.obs = mean.obs, r.burnt = burnt.frac) %>%
    left_join(gee.tbl %>% select(year, gee.ndvi = mean_ndvi,
                                 gee.obs = mean_obs, gee.burnt = burnt_frac),
              by = "year") %>%
    mutate(d.ndvi = r.ndvi - gee.ndvi)
  say(glue("
    Cross-check against the GEE table (it was reduced at 240 m against R's
    30 m, so small differences are expected and large ones are not):
      largest NDVI discrepancy : {sprintf('%+.4f', cmp.csv$d.ndvi[which.max(abs(cmp.csv$d.ndvi))])}"))
  readr::write_csv(cmp.csv, file.path(path.out.ts, "annual_gee_vs_r.csv"))
} else {
  say("    No GEE cross-check table downloaded — optional, but it is the only
    independent check on the stacks. Everything below works without it.")
}

pan <- annual %>% filter(in.panel)

r.year.obs   <- cor(pan$year, pan$mean.obs)
r.year.burnt <- cor(pan$year, pan$burnt.frac)
r.obs.ndvi   <- cor(pan$mean.obs, pan$mean.ndvi)

say(glue("
    WHY AN ABSOLUTE TREND IS NOT AVAILABLE, over the {min(ts.panel)}-{max(ts.panel)} panel:

      correlation of year with observation depth : {sprintf('%+.3f', r.year.obs)}
      correlation of year with burnt fraction    : {sprintf('%+.3f', r.year.burnt)}
      correlation of depth with mean NDVI        : {sprintf('%+.3f', r.obs.ndvi)}

    If the first of those is large, year and depth are close to the same
    variable in this record and no regression on {length(ts.panel)} annual points can
    separate a real greening from a denser archive. That is a property of the
    data, not a shortcoming of the method, and it is why everything that
    follows is framed as a comparison rather than as a trend."))


# =============================================================================
#  PHASE 2 — HOW BIG IS THE 2019->2024 STEP AGAINST ORDINARY VARIATION?
# =============================================================================
#  The headline rests on one pair of years. If a step that size is common in
#  this series then it is a draw from a noisy distribution, and it should be
#  reported with that spread attached.

yr.a <- comp$years[1]
yr.b <- comp$years[2]

step.ab <- annual$mean.ndvi[annual$year == yr.b] -
           annual$mean.ndvi[annual$year == yr.a]

## Every pair of panel years the same distance apart, as a reference set.
gap.yrs <- yr.b - yr.a
pairs.same.gap <- tibble(from = ts.panel, to = ts.panel + gap.yrs) %>%
  filter(to <= max(ts.panel)) %>%
  left_join(annual %>% select(from = year, ndvi.from = mean.ndvi), by = "from") %>%
  left_join(annual %>% select(to   = year, ndvi.to   = mean.ndvi), by = "to") %>%
  mutate(step = ndvi.to - ndvi.from)

yoy <- diff(pan$mean.ndvi)

cat("\n[2] THE 2019->2024 STEP IN CONTEXT\n\n")
pairs.same.gap %>%
  transmute(from, to, step = round(step, 4),
            is.headline = (from == yr.a & to == yr.b)) %>%
  as.data.frame() %>% print(row.names = FALSE)

say(glue("
    the headline step, {yr.a} to {yr.b}         : {sprintf('%+.4f', step.ab)}
    mean of the {nrow(pairs.same.gap)} steps {gap.yrs} years apart : {sprintf('%+.4f', mean(pairs.same.gap$step))}
    sd of those steps                    : {round(sd(pairs.same.gap$step), 4)}
    steps at least as large              : {sum(abs(pairs.same.gap$step) >= abs(step.ab))} of {nrow(pairs.same.gap)}

    year-to-year change, sd              : {round(sd(yoy), 4)}
    largest single-year swing            : {sprintf('%+.4f', yoy[which.max(abs(yoy))])}

    Only {nrow(pairs.same.gap)} comparable steps exist in a {length(ts.panel)}-year panel, so that sd is
    itself a weak estimate and the year-to-year sd over {length(yoy)} differences is the
    sturdier of the two. Neither is a clean null: the reference steps carry the
    same depth confound as the headline one, which makes this a comparison of
    like with like rather than a significance test.

    A headline step that is ordinary against this spread does NOT mean the
    greening is false. It means two dates cannot establish a trend, and the
    figure should be quoted as a difference between two particular years with
    the spread beside it."))


# =============================================================================
#  PHASE 3 — THE PANEL SAMPLE
# =============================================================================
#  One sample, long format, pixel by year. Everything from here on is computed
#  on it so that the numbers are mutually comparable.

f.cand <- file.path(path.out.cls, "change_strong_clusterfree.tif")
if (!file.exists(f.cand)) stop("Candidate raster missing — run 03c_k_sensitivity.R first.")

cand <- mask.to.aoi(terra::rast(f.cand))
names(cand) <- "candidate"

pa.r <- terra::rasterize(aoi.parts.v, cand, field = "pa.label")
names(pa.r) <- "pa"

px.km2 <- prod(terra::res(cand)) / 1e6
cand.km2.total <- terra::global(cand, "sum", na.rm = TRUE)[1, 1] * px.km2

keep <- paste0(rep(c("ndvi.", "obs.", "burnt."), each = length(ts.panel)),
               rep(ts.panel, 3))

wide <- terra::spatSample(
  c(cand, pa.r, ndvi.m[[paste0("ndvi.",  ts.panel)]],
    obs.m[[paste0("obs.",   ts.panel)]],
    burnt.m[[paste0("burnt.", ts.panel)]]),
  n.samp, method = "random", na.rm = TRUE, warn = FALSE) %>%
  as_tibble() %>%
  mutate(pixel = row_number())

long <- wide %>%
  pivot_longer(all_of(keep),
               names_to = c(".value", "year"),
               names_pattern = "^(ndvi|obs|burnt)\\.(\\d{4})$") %>%
  mutate(year = as.integer(year),
         ndvi = ndvi / comp.scale.int)

say(glue("
[3] PANEL SAMPLE
    pixels                 : {format(nrow(wide), big.mark = ',')}
    years                  : {min(ts.panel)}-{max(ts.panel)} ({length(ts.panel)})
    pixel-years            : {format(nrow(long), big.mark = ',')}
    candidate pixels       : {format(sum(wide$candidate), big.mark = ',')} ({round(100 * mean(wide$candidate), 2)}%)
    candidate area on the map : {round(cand.km2.total, 1)} km2"))

if (sum(wide$candidate) < 500) {
  stop(glue("Only {sum(wide$candidate)} candidate pixels in the sample — too few. ",
            "Raise n.samp above {n.samp}."))
}


# =============================================================================
#  PHASE 4 — Q1 AND Q2: DO THE CANDIDATES MOVE DIFFERENTLY?
# =============================================================================
#  The centrepiece. Per pixel, fit NDVI against year and keep two things: the
#  SLOPE, and the share of that pixel's variance the straight line explains.
#
#    high slope, high R2  a pixel moving steadily in one direction — thickening,
#                         regrowth or land-use change
#    high slope, low R2   a pixel that bounces and happens to end higher than it
#                         started — rainfall
#
#  Both are computed in closed form rather than by fitting 200,000 models.

## ols.slope() and ols.r2() are in 00_functions.R section 3.10, where
## scripts/R/tests/test_trend_estimators.R asserts them equal to lm().
fit.pixel <- long %>%
  group_by(pixel, candidate, pa) %>%
  summarise(
    n         = n(),
    slope     = ols.slope(year, ndvi),
    r2        = ols.r2(year, ndvi),
    mean.ndvi = mean(ndvi),
    .groups   = "drop"
  )

traj <- fit.pixel %>%
  group_by(candidate) %>%
  summarise(n = n(),
            median.slope = median(slope),
            mean.slope   = mean(slope),
            pct.rising   = 100 * mean(slope > 0),
            median.r2    = median(r2, na.rm = TRUE),
            pct.r2.over.5 = 100 * mean(r2 > 0.5, na.rm = TRUE),
            .groups = "drop")

cat("\n[4] PER-PIXEL TRAJECTORIES, CANDIDATE vs BACKGROUND\n\n")
traj %>% mutate(across(where(is.numeric), ~round(.x, 5))) %>%
  as.data.frame() %>% print(row.names = FALSE)

slope.cand <- traj$median.slope[traj$candidate == 1]
slope.back <- traj$median.slope[traj$candidate == 0]
did        <- slope.cand - slope.back

## The difference-in-differences, with PA fixed effects so it is not just a
## statement about which reserves the candidates happen to be in. The common
## depth trend and the common rainfall history enter both groups equally and
## subtract out of this contrast — which is why it is the one number here that
## does not inherit the Phase 1 problem.
fit.did <- lm(slope ~ candidate + factor(pa), data = fit.pixel)
did.adj <- coef(fit.did)[["candidate"]]
did.se  <- summary(fit.did)$coefficients["candidate", "Std. Error"]

say(glue("
    median slope, candidates   : {sprintf('%+.6f', slope.cand)} NDVI per year
    median slope, background   : {sprintf('%+.6f', slope.back)} NDVI per year
    difference                 : {sprintf('%+.6f', did)} per year
                                 = {sprintf('%+.5f', did * gap.yrs)} over {gap.yrs} years

    with PA fixed effects      : {sprintf('%+.6f', did.adj)} +/- {sprintf('%.6f', did.se)} per year
                                 = {sprintf('%+.5f', did.adj * gap.yrs)} over {gap.yrs} years

    THIS IS THE CLEAN NUMBER in this script. The archive densified and the rains
    varied, but both did so over the whole landscape, so they affect candidates
    and background alike and cancel in the difference between them.

    linear fit quality is the second half of the question: candidates at
    median R2 {round(traj$median.r2[traj$candidate == 1], 3)} against background {round(traj$median.r2[traj$candidate == 0], 3)}. A
    candidate population that moves steadily should sit clearly higher — a
    straight line genuinely describes its pixels. Similar R2 in both groups
    means the candidates are not following a cleaner trajectory than anywhere
    else, only a steeper one, and steepness alone is what a wet-year/dry-year
    pair would also produce.

    READ THESE R2 FIGURES RELATIVELY, NEVER AGAINST AN ABSOLUTE THRESHOLD.
    Calibrated by simulation in tests/test_trend_estimators.R: over nine years
    at the year-to-year noise this data carries (0.04 NDVI), the median R2 of a
    GENUINE linear trend is 0.16 at +0.005/yr and only 0.39 at +0.010/yr, while
    pure noise sits near 0.09. A real trend does not produce a tidy line over
    nine noisy years, so an R2 of 0.15 is not evidence of no trend — only the
    gap between the two groups carries information."))


# =============================================================================
#  PHASE 5 — Q3: WAS 2019 AN ANOMALOUS DIP FOR THE CANDIDATES?
# =============================================================================
#  A sharp and cheap test, and the one most likely to revise the headline. If
#  the candidates were unusually BROWN in 2019 relative to their own
#  neighbouring years, while the background was not, then part of the
#  2019->2024 difference is a 2019 anomaly in those pixels rather than a gain by
#  2024. The comparison is each group against ITS OWN 2018 and 2020, so the
#  landscape-wide conditions of 2019 drop out.

neigh <- intersect(c(yr.a - 1, yr.a + 1), ts.panel)

if (length(neigh) == 0) {
  stop(glue("Neither {yr.a - 1} nor {yr.a + 1} is in the panel ({min(ts.panel)}-",
            "{max(ts.panel)}), so {yr.a} has no neighbouring years to compare ",
            "against. Widen ts.panel or skip this phase."))
}
if (length(neigh) == 1) {
  say(glue("    Only {neigh} is available as a neighbour of {yr.a}; this test is
    one-sided and weaker than intended."))
}

anom <- long %>%
  filter(year %in% c(yr.a, neigh)) %>%
  mutate(role = if_else(year == yr.a, "focal", "neighbour")) %>%
  group_by(candidate, pixel, role) %>%
  summarise(ndvi = mean(ndvi), .groups = "drop") %>%
  pivot_wider(names_from = role, values_from = ndvi) %>%
  mutate(dip = focal - neighbour) %>%
  group_by(candidate) %>%
  summarise(n = n(), median.dip = median(dip, na.rm = TRUE),
            pct.below = 100 * mean(dip < 0, na.rm = TRUE), .groups = "drop")

cat(glue("\n[5] WAS {yr.a} A DIP? each group against its own ",
         "{paste(neigh, collapse = ' and ')}\n\n"))
anom %>% mutate(across(where(is.numeric), ~round(.x, 5))) %>%
  as.data.frame() %>% print(row.names = FALSE)

dip.cand <- anom$median.dip[anom$candidate == 1]
dip.back <- anom$median.dip[anom$candidate == 0]
dip.excess <- dip.cand - dip.back

say(glue("
    candidate dip in {yr.a}     : {sprintf('%+.5f', dip.cand)}
    background dip in {yr.a}    : {sprintf('%+.5f', dip.back)}
    excess dip for candidates : {sprintf('%+.5f', dip.excess)}

    A negative excess means the candidates were specifically depressed in {yr.a}.
    Against a headline step of {sprintf('%+.4f', step.ab)}, an excess dip of
    {sprintf('%+.5f', dip.excess)} would account for {round(100 * abs(dip.excess) / abs(step.ab), 1)}% of it on its own."))


# =============================================================================
#  PHASE 6 — THE TREND MAP, DEPTH- AND FIRE-ADJUSTED
# =============================================================================
#  A per-pixel slope over the balanced panel, for the report.
#
#  The adjustment is algebra rather than another model. For a fixed set of years
#  the OLS slope is a WEIGHTED SUM of the annual values, with weights
#  (x - xbar) / sum((x - xbar)^2) that do not depend on the data. So the slope
#  of any stack is one linear combination of its bands — a handful of passes
#  instead of 32 million regressions.
#
#  Then, since slope is linear, the slope of NDVI net of depth and fire is
#
#      slope.adj = slope(ndvi) - b.obs * slope(obs) - b.burnt * slope(burnt)
#
#  where b.obs and b.burnt come from the panel below. Those two coefficients are
#  estimated WITHIN YEAR across pixels, with year fixed effects, so the common
#  time trend cannot leak into them. Per script 08's warning that within-date
#  slopes also absorb terrain and ecology, they are an UPPER BOUND on the
#  artefact, which makes the adjusted map a conservative one.

## Fitted on a subsample of pixel-years. lm() builds a dense model matrix, and
## the full panel at 1.8 million rows by a dozen columns would cost several
## hundred megabytes to estimate two coefficients that 400,000 rows pin down
## just as well.
wy.samp <- long %>% slice_sample(n = min(4e5, nrow(long)))

fit.wy <- lm(ndvi ~ obs + burnt + factor(year), data = wy.samp)
b.obs   <- coef(fit.wy)[["obs"]]
b.burnt <- coef(fit.wy)[["burnt"]]

say(glue("
[6] WITHIN-YEAR COEFFICIENTS (upper bounds on the two artefacts)
    NDVI per extra clear look  : {sprintf('%+.6f', b.obs)}
    NDVI if burnt this season  : {sprintf('%+.6f', b.burnt)}"))

w <- ols.slope.weights(ts.panel)

lin.comb <- function(stk, prefix) {
  layers <- paste0(prefix, ts.panel)
  Reduce(`+`, lapply(seq_along(layers),
                     function(i) stk[[layers[i]]] * w[i]))
}

## The balanced panel: every year present and adequately supported. Composition
## changing from year to year would put a trend in the trend map all by itself.
ok <- Reduce(`&`, lapply(ts.panel, function(y)
  !is.na(ndvi.m[[paste0("ndvi.", y)]]) &
   obs.m[[paste0("obs.", y)]] >= comp.min.obs))

slope.ndvi  <- lin.comb(ndvi.m,  "ndvi.") / comp.scale.int
slope.obs   <- lin.comb(obs.m,   "obs.")
slope.burnt <- lin.comb(burnt.m, "burnt.")

slope.adj <- slope.ndvi - b.obs * slope.obs - b.burnt * slope.burnt

slope.ndvi <- terra::mask(slope.ndvi, ok, maskvalues = c(0, NA))
slope.adj  <- terra::mask(slope.adj,  ok, maskvalues = c(0, NA))
names(slope.ndvi) <- "slope.raw"
names(slope.adj)  <- "slope.adj"

pct.balanced <- 100 * terra::global(ok, "mean", na.rm = TRUE)[1, 1]

say(glue("
    pixels in the balanced panel : {round(pct.balanced, 1)}% of the AOI
    mean raw slope               : {sprintf('%+.6f', terra::global(slope.ndvi, 'mean', na.rm = TRUE)[1,1])} NDVI per year
    mean adjusted slope          : {sprintf('%+.6f', terra::global(slope.adj, 'mean', na.rm = TRUE)[1,1])} NDVI per year

    The gap between those two is the share of the apparent greening trend that
    a denser archive and a quieter fire year can account for. Both are reported
    because the adjusted one is conservative, not because either is the truth."))

terra::writeRaster(c(slope.ndvi, slope.adj),
                   file.path(path.out.ts, "ndvi_trend_slope.tif"),
                   overwrite = TRUE, datatype = "FLT4S")


# =============================================================================
#  PHASE 7 — VERDICT
# =============================================================================

cat("\n================================================================\n")
cat(" VERDICT — what the series adds to the two-date result\n")
cat("================================================================\n\n")

step.ordinary  <- abs(step.ab) <= sd(pairs.same.gap$step)
cand.steeper   <- did.adj > 2 * did.se
cand.monotone  <- traj$median.r2[traj$candidate == 1] >
                  traj$median.r2[traj$candidate == 0] + 0.05
dip.matters    <- dip.excess < 0 && abs(dip.excess) > 0.2 * abs(step.ab)

cat(sprintf("  The %d->%d step is %+.4f, against an sd of %.4f across\n",
            yr.a, yr.b, step.ab, sd(pairs.same.gap$step)))
cat(sprintf("  comparable %d-year steps in this record.\n\n", gap.yrs))

cat(sprintf("  Candidates trend %+.6f NDVI/yr faster than background (PA-adjusted,\n",
            did.adj))
cat(sprintf("  SE %.6f), which is %+.5f over %d years.\n\n",
            did.se, did.adj * gap.yrs, gap.yrs))

cat(sprintf("  Candidate linear-fit R2 %.3f vs background %.3f.\n",
            traj$median.r2[traj$candidate == 1],
            traj$median.r2[traj$candidate == 0]))
cat(sprintf("  Excess %d dip for candidates: %+.5f.\n\n", yr.a, dip.excess))

if (cand.steeper && cand.monotone) {
  cat("  THE CANDIDATES ARE ON A REAL TRAJECTORY.\n\n")
  cat("  They trend faster than their background and a straight line describes\n")
  cat("  them better, which is what thickening, regrowth or land-use change\n")
  cat("  looks like and is not what rainfall looks like. The common confounds\n")
  cat("  cancel in that contrast. What this still cannot do is say WHICH of\n")
  cat("  those three it is — that needs high-resolution imagery.\n")
} else if (cand.steeper && !cand.monotone) {
  cat("  STEEPER, BUT NOT STEADIER.\n\n")
  cat("  The candidates do move faster than their background, so something\n")
  cat("  distinguishes them. But a straight line describes them no better than\n")
  cat("  it describes anywhere else, which is the signature of pixels that\n")
  cat("  swing with the rains and happened to be caught low in one year and\n")
  cat("  high in the other. Treat the candidate area as an upper bound and\n")
  cat("  quote the headline as a two-year difference, not a rate.\n")
} else {
  cat("  THE CANDIDATES ARE NOT DISTINGUISHABLE OVER THE SERIES.\n\n")
  cat("  Their trend is not separable from the background's once protected\n")
  cat("  area is held constant. The two-date comparison picked out pixels that\n")
  cat("  differed between 2019 and 2024 specifically, and over nine years they\n")
  cat("  do not behave differently from anywhere else. That points at\n")
  cat("  inter-annual variability rather than land-cover change, and the\n")
  cat("  candidate area should not be reported as change without independent\n")
  cat("  evidence.\n")
}

if (dip.matters) {
  cat(sprintf("\n  NOTE: %d was specifically depressed in the candidate pixels, by\n", yr.a))
  cat(sprintf("  %+.5f beyond the background's own dip — %.0f%% of the headline step.\n",
              dip.excess, 100 * abs(dip.excess) / abs(step.ab)))
  cat("  Part of what was measured as change is that anomaly.\n")
}

if (step.ordinary) {
  cat("\n  AND: the headline step is within one sd of comparable steps in this\n")
  cat("  record, so it should be quoted as a difference between two particular\n")
  cat("  years with that spread attached, never as an annual rate.\n")
}


# =============================================================================
#  PHASE 8 — FIGURES
# =============================================================================

hl <- annual %>% filter(year %in% c(yr.a, yr.b))

p.series <- annual %>%
  ggplot(aes(year, mean.ndvi)) +
  geom_line(colour = "grey60") +
  geom_point(aes(shape = in.panel, colour = burnt.frac, size = mean.obs)) +
  ## inherit.aes = FALSE on both: the top-level aes carries size and colour
  ## mappings, and letting them reach these two layers would collide with the
  ## fixed size and colour set here.
  geom_point(data = hl, aes(year, mean.ndvi), inherit.aes = FALSE,
             shape = 21, size = 6, stroke = 0.9, colour = "#762a83",
             fill = NA) +
  geom_text(data = hl, aes(year, mean.ndvi, label = year), inherit.aes = FALSE,
            vjust = -1.8, size = 3, colour = "#762a83") +
  scale_colour_viridis_c(option = "inferno", direction = -1, end = 0.9,
                         name = "burnt fraction") +
  scale_size_continuous(range = c(1.5, 6), name = "clear looks") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     name = "in panel") +
  scale_x_continuous(breaks = ts.years) +
  labs(title = "Mean dry-season NDVI by year",
       subtitle = glue("Circled: the two dates the headline rests on. ",
                       "Point size is observation depth, colour is fire — ",
                       "both confound the apparent trend."),
       x = NULL, y = "mean NDVI") +
  theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_blank())

traj.series <- long %>%
  group_by(candidate, year) %>%
  summarise(mean.ndvi = mean(ndvi), .groups = "drop") %>%
  mutate(group = if_else(candidate == 1, "candidate", "background"))

p.traj <- traj.series %>%
  ggplot(aes(year, mean.ndvi, colour = group)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  geom_vline(xintercept = c(yr.a, yr.b), linetype = 3, colour = "#762a83") +
  scale_colour_manual(values = c(candidate = "#d7191c",
                                 background = "#2166ac"), name = NULL) +
  scale_x_continuous(breaks = ts.panel) +
  labs(title = "Candidate patches against their background",
       subtitle = glue("Dotted: the two headline years. Parallel lines mean ",
                       "the candidates are not on a different trajectory."),
       x = NULL, y = "mean NDVI") +
  theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        legend.position = "bottom")

p.slope <- fit.pixel %>%
  mutate(group = if_else(candidate == 1, "candidate", "background")) %>%
  ggplot(aes(slope, fill = group)) +
  geom_vline(xintercept = 0, colour = "grey40") +
  geom_density(alpha = 0.5, colour = NA) +
  scale_fill_manual(values = c(candidate = "#d7191c",
                               background = "#2166ac"), name = NULL) +
  coord_cartesian(xlim = quantile(fit.pixel$slope, c(0.005, 0.995))) +
  labs(title = "Per-pixel NDVI trend, candidate vs background",
       subtitle = glue("Slopes over {min(ts.panel)}-{max(ts.panel)}. ",
                       "Overlapping distributions mean the two are not ",
                       "behaving differently."),
       x = "NDVI change per year", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        axis.text.y = element_blank(), legend.position = "bottom")

lim <- quantile(terra::spatSample(slope.adj, 2e5, na.rm = TRUE)[, 1],
                c(0.02, 0.98), na.rm = TRUE)
lim <- max(abs(lim))

p.map <- ggplot() +
  geom_spatraster(data = slope.adj, maxcell = 8e5) +
  scale_fill_gradient2(low = "#8c510a", mid = "#f6f6f6", high = "#01665e",
                       midpoint = 0, limits = c(-lim, lim), oob = scales::squish,
                       na.value = NA, name = "NDVI / yr") +
  geom_sf(data = aoi$parts, fill = NA, colour = "grey25", linewidth = 0.3) +
  annotation_scale(location = "bl", width_hint = 0.25,
                   height = unit(0.15, "cm")) +
  coord_sf(expand = FALSE) +
  labs(title = glue("NDVI trend {min(ts.panel)}-{max(ts.panel)}, adjusted for ",
                    "observation depth and fire"),
       subtitle = glue("Balanced panel only ({round(pct.balanced, 1)}% of the ",
                       "AOI). Adjustment is conservative — see Phase 6.")) +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8),
        plot.background = element_rect(fill = "white", colour = NA))

ggsave(file.path(path.out.ts, "ts_01_annual_series.png"), p.series,
       width = 8, height = 4.8, dpi = 200, bg = "white")
ggsave(file.path(path.out.ts, "ts_02_candidate_trajectories.png"), p.traj,
       width = 7, height = 4.8, dpi = 200, bg = "white")
ggsave(file.path(path.out.ts, "ts_03_slope_distributions.png"), p.slope,
       width = 7, height = 4.4, dpi = 200, bg = "white")
ggsave(file.path(path.out.ts, "ts_04_trend_map.png"), p.map,
       width = 7.5, height = 8, dpi = 220, bg = "white")


# =============================================================================
#  PHASE 9 — OUTPUTS
# =============================================================================

readr::write_csv(annual,         file.path(path.out.ts, "annual_series.csv"))
readr::write_csv(pairs.same.gap, file.path(path.out.ts, "comparable_steps.csv"))
readr::write_csv(traj,           file.path(path.out.ts, "trajectory_summary.csv"))
readr::write_csv(anom,           file.path(path.out.ts, "year_2019_anomaly.csv"))
readr::write_csv(traj.series,    file.path(path.out.ts, "candidate_vs_background_series.csv"))

readr::write_csv(
  tibble(
    quantity = c("headline.step", "sd.comparable.steps", "sd.year.on.year",
                 "cor.year.obs", "cor.year.burnt", "cor.obs.ndvi",
                 "median.slope.candidate", "median.slope.background",
                 "did.raw", "did.pa.adjusted", "did.se",
                 "median.r2.candidate", "median.r2.background",
                 "dip.candidate", "dip.background", "dip.excess",
                 "b.obs.within.year", "b.burnt.within.year",
                 "pct.balanced.panel", "cand.km2"),
    value = c(step.ab, sd(pairs.same.gap$step), sd(yoy),
              r.year.obs, r.year.burnt, r.obs.ndvi,
              slope.cand, slope.back,
              did, did.adj, did.se,
              traj$median.r2[traj$candidate == 1],
              traj$median.r2[traj$candidate == 0],
              dip.cand, dip.back, dip.excess,
              b.obs, b.burnt,
              pct.balanced, cand.km2.total)
  ),
  file.path(path.out.ts, "timeseries_headline_numbers.csv")
)

cat("\n================================================================\n")
cat(" DONE — tables, figures and the trend raster in", path.out.ts, "\n")
cat("================================================================\n\n")

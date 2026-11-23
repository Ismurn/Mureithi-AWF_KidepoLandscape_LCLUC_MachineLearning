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

#  THE SAMPLE COMES FIRST, and it carries ALL eleven years rather than just the
#  panel, because the annual series below needs the thin years in order to show
#  why they cannot be used.
#
#  Everything tabular in this script is computed from this one sample. At
#  400,000 pixels an annual mean is pinned to about +/-0.0001 NDVI, far finer
#  than any difference being discussed, and it costs one pass over the stacks
#  instead of the several that repeated terra::global() calls over 33 bands of
#  32 million cells would take. The full rasters are touched again only in
#  Phase 6, where a per-pixel map is genuinely the output.
#
#  na.rm = FALSE is deliberate. Dropping rows with any NA would silently impose
#  the intersection of all eleven years' footprints, which is exactly the
#  composition effect this phase exists to measure.

f.cand <- file.path(path.out.cls, "change_strong_clusterfree.tif")
if (!file.exists(f.cand)) stop("Candidate raster missing — run 03c_k_sensitivity.R first.")

cand <- mask.to.aoi(terra::rast(f.cand))
names(cand) <- "candidate"

pa.r <- terra::rasterize(aoi.parts.v, cand, field = "pa.label")
names(pa.r) <- "pa"

px.km2 <- prod(terra::res(cand)) / 1e6
cand.km2.total <- terra::global(cand, "sum", na.rm = TRUE)[1, 1] * px.km2

all.r <- c(cand, pa.r, ndvi.m, obs.m, burnt.m)
names(all.r) <- c("candidate", "pa",
                  paste0("ndvi.",  ts.years),
                  paste0("obs.",   ts.years),
                  paste0("burnt.", ts.years))

wide <- terra::spatSample(all.r, n.samp, method = "random",
                          na.rm = FALSE, warn = FALSE) %>%
  as_tibble() %>%
  filter(!is.na(candidate)) %>%        # inside the PA union
  mutate(pixel = row_number())

long.all <- wide %>%
  pivot_longer(matches("^(ndvi|obs|burnt)\\.\\d{4}$"),
               names_to = c(".value", "year"),
               names_pattern = "^(ndvi|obs|burnt)\\.(\\d{4})$") %>%
  mutate(year = as.integer(year),
         ndvi = ndvi / comp.scale.int,
         supported = !is.na(ndvi) & !is.na(obs) & obs >= comp.min.obs)

say(glue("
[1] SAMPLE
    pixels inside the PA union : {format(nrow(wide), big.mark = ',')}
    pixel-years, all {length(ts.years)} years  : {format(nrow(long.all), big.mark = ',')}
    candidate pixels           : {format(sum(wide$candidate, na.rm = TRUE), big.mark = ',')} ({round(100 * mean(wide$candidate, na.rm = TRUE), 2)}%)
    candidate area on the map  : {round(cand.km2.total, 1)} km2"))


# =============================================================================
#  PHASE 1b — THE ANNUAL SERIES, ON TWO DIFFERENT SETS OF PIXELS
# =============================================================================
#  This phase exists because of something the GEE cross-check table showed
#  before any of it was run: 2020 came back at mean NDVI 0.491 and 2026 at
#  0.528, against roughly 0.38 for every other year. Both are LOW-COVERAGE
#  years — 4.77 and 7.90 mean clear looks, and only 75.9% of the AOI
#  adequately covered in 2020.
#
#  That is not a greener landscape, it is a DIFFERENT landscape. A year's mean
#  taken over whatever pixels happened to be clear is a mean over a different
#  set of pixels each year, and in a year when a quarter of the AOI is missing
#  the surviving quarter is not a random quarter — script 08 already found that
#  clear-sky selection favours particular ground. Comparing such means across
#  years measures the weather's effect on WHICH PIXELS WERE SEEN as much as on
#  the vegetation.
#
#  So both versions are computed and printed side by side:
#
#    VARYING  each year averaged over its own supported pixels — what the GEE
#             table reports, and what any naive series would give
#    FIXED    every year averaged over the SAME pixels, those supported in all
#             panel years. Composition is then constant by construction and the
#             years are comparable.
#
#  The FIXED series is the one every later phase uses.

fixed.pixels <- long.all %>%
  filter(year %in% ts.panel) %>%
  group_by(pixel) %>%
  summarise(all.ok = all(supported), .groups = "drop") %>%
  filter(all.ok) %>%
  pull(pixel)

annual.varying <- long.all %>%
  filter(supported) %>%
  group_by(year) %>%
  summarise(n = n(), mean.ndvi = mean(ndvi), mean.obs = mean(obs),
            burnt.frac = mean(burnt, na.rm = TRUE), .groups = "drop")

annual.fixed <- long.all %>%
  filter(pixel %in% fixed.pixels) %>%
  group_by(year) %>%
  summarise(n = sum(!is.na(ndvi)), mean.ndvi = mean(ndvi, na.rm = TRUE),
            mean.obs = mean(obs, na.rm = TRUE),
            burnt.frac = mean(burnt, na.rm = TRUE), .groups = "drop")

annual <- annual.fixed %>%
  select(year, mean.ndvi, mean.obs, burnt.frac) %>%
  left_join(annual.varying %>% select(year, ndvi.varying = mean.ndvi,
                                      n.varying = n), by = "year") %>%
  mutate(in.panel    = year %in% ts.panel,
         provisional = year %in% ts.provisional,
         composition = ndvi.varying - mean.ndvi)

cat("\n[1b] THE ANNUAL SERIES — FIXED PIXELS vs VARYING\n\n")
annual %>%
  transmute(year,
            ndvi.fixed = round(mean.ndvi, 4),
            ndvi.varying = round(ndvi.varying, 4),
            composition = round(composition, 4),
            mean.obs = round(mean.obs, 2),
            burnt.frac = round(burnt.frac, 3),
            in.panel, provisional) %>%
  as.data.frame() %>% print(row.names = FALSE)

say(glue("
    pixels held fixed across {min(ts.panel)}-{max(ts.panel)} : {format(length(fixed.pixels), big.mark = ',')} of {format(nrow(wide), big.mark = ',')} ({round(100 * length(fixed.pixels) / nrow(wide), 1)}%)
    largest composition effect       : {sprintf('%+.4f', annual$composition[which.max(abs(annual$composition))])} NDVI, in {annual$year[which.max(abs(annual$composition))]}

    The `composition` column is what you get for free by averaging over
    whichever pixels were visible. Where it is large, the varying series is
    telling you about cloud cover, not about vegetation. Every phase below
    uses the fixed series."))

## The cross-check against the GEE table, if it was downloaded. Three code
## paths agreeing is the only reason to believe any of them.
#  The comparison is against the VARYING series, not the fixed one. GEE
#  averaged over whatever pixels it could see in each year, so that is the
#  quantity to hold it to; comparing it against the fixed series would show a
#  large difference that is this script's own doing.
f.csv <- file.path(path.gee.out, "kidepo_ndvi_annual_aoi.csv")
if (file.exists(f.csv)) {
  gee.tbl <- readr::read_csv(f.csv, show_col_types = FALSE)
  cmp.csv <- annual %>%
    select(year, r.ndvi.varying = ndvi.varying, r.ndvi.fixed = mean.ndvi,
           r.obs = mean.obs, r.burnt = burnt.frac) %>%
    left_join(gee.tbl %>% select(year, gee.ndvi = mean_ndvi,
                                 gee.obs = mean_obs, gee.burnt = burnt_frac),
              by = "year") %>%
    mutate(d.ndvi = r.ndvi.varying - gee.ndvi,
           d.burnt = r.burnt - gee.burnt)
  say(glue("
    Cross-check against the GEE table. GEE reduced at 240 m over the BUFFERED
    region while this is 30 m over the PA union, so small differences are
    expected; a large one would mean the download or the band order is wrong.
      largest NDVI discrepancy  : {sprintf('%+.4f', cmp.csv$d.ndvi[which.max(abs(cmp.csv$d.ndvi))])} (in {cmp.csv$year[which.max(abs(cmp.csv$d.ndvi))]})
      largest burnt discrepancy : {sprintf('%+.4f', cmp.csv$d.burnt[which.max(abs(cmp.csv$d.burnt))])} (in {cmp.csv$year[which.max(abs(cmp.csv$d.burnt))]})"))
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
#  PHASE 3 — THE BALANCED PANEL
# =============================================================================
#  `long` is the working table for every phase from here: the panel years only,
#  restricted to the pixels held fixed in Phase 1b, so each pixel contributes a
#  complete and adequately supported series. A pixel appearing in some years and
#  not others would put a trend into a trend estimate all by itself.

long <- long.all %>%
  filter(year %in% ts.panel, pixel %in% fixed.pixels)

n.cand.px <- long %>% filter(year == min(ts.panel)) %>% pull(candidate) %>% sum()

say(glue("
[3] BALANCED PANEL
    pixels                 : {format(length(fixed.pixels), big.mark = ',')}
    years                  : {min(ts.panel)}-{max(ts.panel)} ({length(ts.panel)})
    pixel-years            : {format(nrow(long), big.mark = ',')}
    candidate pixels       : {format(n.cand.px, big.mark = ',')} ({round(100 * n.cand.px / length(fixed.pixels), 2)}%)"))

if (n.cand.px < 500) {
  stop(glue("Only {n.cand.px} candidate pixels survive into the balanced panel ",
            "— too few to compare against background. Raise n.samp above ",
            "{n.samp}, or widen ts.panel if whole years are being lost."))
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
#  PHASE 5b — HOW MUCH OF ALL THAT IS THE SELECTION RULE ITSELF?
# =============================================================================
#  This phase exists because Phase 5 produced an impossible-looking number: an
#  excess 2019 dip of -0.067, which is THREE TIMES the headline step it is
#  supposed to help explain. An effect several times larger than the thing it
#  explains is usually a sign that the estimator is measuring the selection
#  rather than the world, and that is exactly what is happening.
#
#  The candidates were DEFINED as pixels with |dNDVI| > 0.15 between 2019 and
#  2024. Selecting pixels on a difference between two years guarantees two
#  things about them, whether or not anything real is going on:
#
#    they will look anomalously LOW in 2019 and HIGH in 2024, because that is
#    the criterion they were chosen by — ordinary regression to the mean
#
#    they will show a POSITIVE fitted trend across any window containing those
#    two years, because the fit is partly a restatement of the difference it
#    was selected on
#
#  So Phases 4 and 5 are both contaminated, and no amount of fixed effects
#  repairs it: the bias is in the sample, not the model.
#
#  TWO CORRECTIONS, both reported.
#
#  LEAVE-OUT   refit the trend with 2019 and 2024 DROPPED. The selection no
#              longer mechanically determines the fit. Not perfectly clean —
#              a pixel unusually brown in 2019 may be brown in 2018 too — but
#              it removes the dominant term.
#
#  PLACEBO     apply the SAME rule to other year pairs five years apart, and
#              see whether those placebo "candidates" behave the same way. If
#              they show the same excess dip and the same steeper trend, the
#              pattern belongs to the rule and not to the real patches. This is
#              the more decisive of the two.

nd.wide <- long %>%
  select(pixel, pa, candidate, year, ndvi) %>%
  pivot_wider(names_from = year, values_from = ndvi, names_prefix = "y")

mag.threshold <- 0.15   # as in 03c, which defined the candidate set

#' Behaviour of a set selected on the (from -> to) difference.
#' Returns the leave-out trend gap and the excess dip at `from`, each measured
#' against the pixels the same rule did NOT select.
selection.probe <- function(from, to) {

  yf <- paste0("y", from); yt <- paste0("y", to)
  flagged <- (nd.wide[[yt]] - nd.wide[[yf]]) > mag.threshold

  lo.years <- setdiff(ts.panel, c(from, to))
  lo.cols  <- paste0("y", lo.years)
  m <- as.matrix(nd.wide[, lo.cols])

  w  <- ols.slope.weights(lo.years)
  sl <- as.numeric(m %*% w)          # exact OLS slope, one matrix product

  ## Dip at `from`, against its own neighbouring years, excluding `to`.
  nb <- setdiff(intersect(c(from - 1, from + 1), ts.panel), to)
  dip <- if (length(nb) > 0) {
    nd.wide[[yf]] - rowMeans(as.matrix(nd.wide[, paste0("y", nb)]))
  } else rep(NA_real_, nrow(nd.wide))

  d <- tibble(flagged = flagged, sl = sl, dip = dip, pa = nd.wide$pa) %>%
    drop_na(sl)

  fit <- lm(sl ~ flagged + factor(pa), data = d)

  tibble(
    from = from, to = to,
    n.flagged = sum(d$flagged),
    pct.flagged = 100 * mean(d$flagged),
    slope.flagged = median(d$sl[d$flagged]),
    slope.rest = median(d$sl[!d$flagged]),
    slope.gap.pa.adj = coef(fit)[["flaggedTRUE"]],
    slope.gap.se = summary(fit)$coefficients["flaggedTRUE", "Std. Error"],
    excess.dip = median(d$dip[d$flagged], na.rm = TRUE) -
                 median(d$dip[!d$flagged], na.rm = TRUE)
  )
}

## The real pair, then every other pair the same distance apart. The earliest
## pair gets a one-sided dip test, since its `from` year has no panel year
## before it — weaker, but a third placebo is worth more than a tidy design
## with only two.
probe.pairs <- tibble(from = ts.panel, to = ts.panel + gap.yrs) %>%
  filter(to <= max(ts.panel))

probes <- purrr::pmap_dfr(list(probe.pairs$from, probe.pairs$to),
                          selection.probe) %>%
  mutate(is.real = (from == yr.a & to == yr.b))

cat("\n[5b] THE SAME RULE APPLIED TO OTHER YEAR PAIRS\n\n")
probes %>%
  transmute(pair = paste0(from, "->", to), is.real,
            pct.flagged = round(pct.flagged, 2),
            leaveout.slope.gap = round(slope.gap.pa.adj, 6),
            se = round(slope.gap.se, 6),
            excess.dip = round(excess.dip, 4)) %>%
  as.data.frame() %>% print(row.names = FALSE)

real.probe <- probes %>% filter(is.real)
fake.probe <- probes %>% filter(!is.real)

did.lo    <- real.probe$slope.gap.pa.adj
did.lo.se <- real.probe$slope.gap.se
dip.lo    <- real.probe$excess.dip

say(glue("
    REAL pair {yr.a}->{yr.b}
      leave-out trend gap    : {sprintf('%+.6f', did.lo)} +/- {sprintf('%.6f', did.lo.se)} per year
      against the in-sample  : {sprintf('%+.6f', did.adj)} per year from Phase 4
      excess dip at {yr.a}      : {sprintf('%+.5f', dip.lo)}

    PLACEBO pairs ({nrow(fake.probe)} of them), which select on years the real candidates
    were not chosen by:
      leave-out trend gap, mean : {sprintf('%+.6f', mean(fake.probe$slope.gap.pa.adj))} (range {sprintf('%+.6f', min(fake.probe$slope.gap.pa.adj))} to {sprintf('%+.6f', max(fake.probe$slope.gap.pa.adj))})
      excess dip, mean          : {sprintf('%+.5f', mean(fake.probe$excess.dip))} (range {sprintf('%+.5f', min(fake.probe$excess.dip))} to {sprintf('%+.5f', max(fake.probe$excess.dip))})

    HOW TO READ THIS. The placebo rows are what the selection rule produces
    from a landscape where nothing has been singled out — they are the null.
    If the real pair's numbers sit inside the placebo range, then everything
    Phases 4 and 5 found about the candidates is an artefact of how they were
    chosen, and the series says nothing about them. If the real pair stands
    outside that range, the excess is genuine."))

sel.explains.slope <- did.lo <= max(fake.probe$slope.gap.pa.adj)
sel.explains.dip   <- dip.lo >= min(fake.probe$excess.dip)


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
#  NO DEPTH CORRECTION IS APPLIED TO THIS MAP, and the reason is a result in
#  itself rather than an omission.
#
#  The plan was slope.adj = slope(ndvi) - b.obs*slope(obs) - b.burnt*slope(burnt),
#  with b.obs and b.burnt estimated within year across pixels so the common time
#  trend could not leak into them. Estimated that way, b.obs came out NEGATIVE:
#  about -0.0005 NDVI per extra clear look, against +0.0033 from the differenced
#  regression in script 08 that identified the artefact in the first place.
#
#  Opposite signs mean these two numbers are not measuring the same thing, and
#  script 08 said in advance which one is which: a within-date, cross-pixel
#  slope measures the LANDSCAPE — pixels that get seen more often are drier,
#  flatter, less cloudy ground — while differencing the same pixel across dates
#  measures the artefact. Subtracting a wrong-signed coefficient would not
#  correct this map, it would inflate it, which is exactly what the first run
#  did: it moved the mean slope from +0.0059 to +0.0064 per year.
#
#  And there is no better coefficient available. The within-PIXEL-over-time
#  estimator, which would be the right analogue of script 08, is collinear with
#  the year term by construction here — that is the whole thesis of Phase 1.
#  So the raw slope is the only honest map, and the coefficients below are
#  printed as a diagnostic, not used.

## Fitted on a subsample of pixel-years. lm() builds a dense model matrix, and
## the full panel at 1.8 million rows by a dozen columns would cost several
## hundred megabytes to estimate two coefficients that 400,000 rows pin down
## just as well.
wy.samp <- long %>% slice_sample(n = min(4e5, nrow(long)))

fit.wy <- lm(ndvi ~ obs + burnt + factor(year), data = wy.samp)
b.obs   <- coef(fit.wy)[["obs"]]
b.burnt <- coef(fit.wy)[["burnt"]]

say(glue("
[6] WITHIN-YEAR COEFFICIENTS — DIAGNOSTIC ONLY, NOT APPLIED
    NDVI per extra clear look  : {sprintf('%+.6f', b.obs)}
    NDVI if burnt this season  : {sprintf('%+.6f', b.burnt)}

    Script 08's differenced estimate of the depth effect is about +0.0033 per
    look. This one has the OPPOSITE SIGN, which is the signature script 08
    warned of: across pixels within a year, the ones seen more often are
    simply different ground. Using it to correct the trend map would inflate
    the trend rather than correct it, so the map below is the RAW slope and
    carries no depth correction. The fire coefficient is left unapplied for
    the same reason — it is identified the same way."))

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

slope.ndvi <- lin.comb(ndvi.m, "ndvi.") / comp.scale.int
slope.ndvi <- terra::mask(slope.ndvi, ok, maskvalues = c(0, NA))
names(slope.ndvi) <- "slope.raw"

pct.balanced <- 100 * terra::global(ok, "mean", na.rm = TRUE)[1, 1]
mean.slope   <- terra::global(slope.ndvi, "mean", na.rm = TRUE)[1, 1]

say(glue("
    pixels in the balanced panel : {round(pct.balanced, 1)}% of the AOI
    mean slope over those pixels : {sprintf('%+.6f', mean.slope)} NDVI per year
      implied over {gap.yrs} years       : {sprintf('%+.4f', mean.slope * gap.yrs)}

    Requiring every one of the {length(ts.panel)} years to be adequately supported keeps only
    {round(pct.balanced, 1)}% of the AOI, so this map is not a map of the landscape — it is a
    map of the part of it that is reliably visible every year, which is flatter
    and less cloudy than the whole. Read the pattern, not the coverage.

    And read the LEVEL with the Phase 1 warning attached: observation depth
    roughly triples across this record, so a positive mean slope here is partly
    the archive growing and cannot be separated from real greening."))

terra::writeRaster(slope.ndvi,
                   file.path(path.out.ts, "ndvi_trend_slope.tif"),
                   overwrite = TRUE, datatype = "FLT4S")


# =============================================================================
#  PHASE 7 — VERDICT
# =============================================================================

cat("\n================================================================\n")
cat(" VERDICT — what the series adds to the two-date result\n")
cat("================================================================\n\n")

step.ordinary <- abs(step.ab) <= sd(yoy)

cat(sprintf("  The %d->%d step on a FIXED set of pixels is %+.4f, against a\n",
            yr.a, yr.b, step.ab))
cat(sprintf("  year-to-year sd of %.4f in this record. %d of %d comparable\n",
            sd(yoy), sum(abs(pairs.same.gap$step) >= abs(step.ab)),
            nrow(pairs.same.gap)))
cat(sprintf("  %d-year steps are at least as large.\n\n", gap.yrs))

cat("  ON THE CANDIDATES, the in-sample contrasts and the selection-free ones\n")
cat("  disagree, and only the second kind can be believed:\n\n")
cat(sprintf("    trend gap, window includes %d and %d : %+.6f /yr\n",
            yr.a, yr.b, did.adj))
cat(sprintf("    trend gap, those two years dropped    : %+.6f /yr (SE %.6f)\n",
            did.lo, did.lo.se))
cat(sprintf("    same, from %d placebo year pairs       : %+.6f to %+.6f /yr\n",
            nrow(fake.probe), min(fake.probe$slope.gap.pa.adj),
            max(fake.probe$slope.gap.pa.adj)))
cat(sprintf("    excess %d dip                        : %+.5f\n", yr.a, dip.lo))
cat(sprintf("    same, from the placebo pairs           : %+.5f to %+.5f\n\n",
            min(fake.probe$excess.dip), max(fake.probe$excess.dip)))

if (sel.explains.slope && sel.explains.dip) {
  cat("  THE SERIES DOES NOT SINGLE THE CANDIDATES OUT.\n\n")
  cat("  Both things that looked like evidence for them — a steeper trend and a\n")
  cat("  2019 dip — are reproduced by applying the same threshold to year pairs\n")
  cat("  the candidates were not selected by. They are properties of the\n")
  cat("  SELECTION RULE, not of the ground: choose pixels on a difference\n")
  cat("  between two years and they will always look low in the first, high in\n")
  cat("  the second, and rising in between.\n\n")
  cat("  The candidate area therefore still rests on the two-date comparison\n")
  cat("  alone, and the series adds no independent support for it. It does add\n")
  cat("  something about the landscape: the step is ordinary against year-to-\n")
  cat("  year variation, so the figure belongs in a report as a difference\n")
  cat("  between two particular years with that spread beside it.\n")
} else if (!sel.explains.slope && !sel.explains.dip) {
  cat("  THE CANDIDATES SURVIVE THE PLACEBO.\n\n")
  cat("  Their trend is steeper than the placebo selections produce, and their\n")
  cat("  dip is deeper, so neither is explained by the selection rule. That is\n")
  cat("  real evidence of a different trajectory, and with the common depth and\n")
  cat("  rainfall trends cancelling in the contrast it is the strongest result\n")
  cat("  this pipeline has produced. It still cannot say whether the cause is\n")
  cat("  thickening, regrowth or land-use change.\n")
} else if (!sel.explains.slope) {
  cat("  A STEEPER TREND SURVIVES, THE DIP DOES NOT.\n\n")
  cat("  The candidates trend faster than the placebo selections manage, so\n")
  cat("  something distinguishes them. The 2019 dip, though, is within what the\n")
  cat("  selection rule manufactures on its own and should not be reported as a\n")
  cat("  finding about 2019.\n")
} else {
  cat("  A DEEPER DIP SURVIVES, THE TREND DOES NOT.\n\n")
  cat("  The candidates were specifically depressed in 2019 beyond what the\n")
  cat("  selection rule explains, while their trend is not distinguishable from\n")
  cat("  what that rule manufactures. That points at a 2019 anomaly in those\n")
  cat("  pixels rather than a trajectory, and the candidate area is an upper\n")
  cat("  bound.\n")
}

if (step.ordinary) {
  cat("\n  AND: the headline step is within one year-to-year sd, so it must never\n")
  cat("  be quoted as an annual rate.\n")
}


# =============================================================================
#  PHASE 8 — FIGURES
# =============================================================================

hl <- annual %>% filter(year %in% c(yr.a, yr.b))

p.series <- annual %>%
  ggplot(aes(year, mean.ndvi)) +
  ## The varying-footprint series, for contrast: where it departs from the
  ## fixed one it is reporting cloud cover rather than vegetation.
  geom_line(aes(y = ndvi.varying), colour = "#d7191c", linetype = 2,
            linewidth = 0.4) +
  geom_point(aes(y = ndvi.varying), colour = "#d7191c", size = 1,
             shape = 4) +
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
       subtitle = glue("Solid: the same {format(length(fixed.pixels), big.mark = ',')} pixels every year. ",
                       "Red dashed: each year over whatever pixels were clear ",
                       "— cloud cover, not vegetation.\nCircled: the two dates ",
                       "the headline rests on. Point size is observation depth, ",
                       "colour is fire."),
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

lim <- quantile(terra::spatSample(slope.ndvi, 2e5, na.rm = TRUE)[, 1],
                c(0.02, 0.98), na.rm = TRUE)
lim <- max(abs(lim))

p.map <- ggplot() +
  geom_spatraster(data = slope.ndvi, maxcell = 8e5) +
  scale_fill_gradient2(low = "#8c510a", mid = "#f6f6f6", high = "#01665e",
                       midpoint = 0, limits = c(-lim, lim), oob = scales::squish,
                       na.value = NA, name = "NDVI / yr") +
  geom_sf(data = aoi$parts, fill = NA, colour = "grey25", linewidth = 0.3) +
  annotation_scale(location = "bl", width_hint = 0.25,
                   height = unit(0.15, "cm")) +
  coord_sf(expand = FALSE) +
  labs(title = glue("NDVI trend {min(ts.panel)}-{max(ts.panel)} (uncorrected)"),
       subtitle = glue("Balanced panel only — {round(pct.balanced, 1)}% of the AOI, ",
                       "the part reliably visible every year.\nNo depth ",
                       "correction is available: see Phase 6. The LEVEL is ",
                       "partly the growing archive; read the pattern.")) +
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
readr::write_csv(probes,         file.path(path.out.ts, "selection_placebo.csv"))

readr::write_csv(
  tibble(
    quantity = c("headline.step.fixed.pixels", "sd.comparable.steps",
                 "sd.year.on.year", "largest.single.year.swing",
                 "cor.year.obs", "cor.year.burnt", "cor.obs.ndvi",
                 "median.slope.candidate", "median.slope.background",
                 "did.raw.in.sample", "did.pa.adjusted.in.sample", "did.se",
                 "did.leaveout", "did.leaveout.se",
                 "did.placebo.min", "did.placebo.max",
                 "median.r2.candidate", "median.r2.background",
                 "dip.candidate", "dip.background", "dip.excess.in.sample",
                 "dip.excess.leaveout",
                 "dip.placebo.min", "dip.placebo.max",
                 "selection.explains.trend", "selection.explains.dip",
                 "b.obs.within.year.unused", "b.burnt.within.year.unused",
                 "pct.balanced.panel", "mean.slope.balanced.panel",
                 "cand.km2"),
    value = c(step.ab, sd(pairs.same.gap$step),
              sd(yoy), yoy[which.max(abs(yoy))],
              r.year.obs, r.year.burnt, r.obs.ndvi,
              slope.cand, slope.back,
              did, did.adj, did.se,
              did.lo, did.lo.se,
              min(fake.probe$slope.gap.pa.adj), max(fake.probe$slope.gap.pa.adj),
              traj$median.r2[traj$candidate == 1],
              traj$median.r2[traj$candidate == 0],
              dip.cand, dip.back, dip.excess,
              dip.lo,
              min(fake.probe$excess.dip), max(fake.probe$excess.dip),
              as.numeric(sel.explains.slope), as.numeric(sel.explains.dip),
              b.obs, b.burnt,
              pct.balanced, mean.slope,
              cand.km2.total)
  ),
  file.path(path.out.ts, "timeseries_headline_numbers.csv")
)

cat("\n================================================================\n")
cat(" DONE — tables, figures and the trend raster in", path.out.ts, "\n")
cat("================================================================\n\n")

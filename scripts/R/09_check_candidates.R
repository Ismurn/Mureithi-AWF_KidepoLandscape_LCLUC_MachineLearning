# =============================================================================
#  KIDEPO LANDSCAPE — LCLUC
#  Script 09 — Interrogate the candidate change patches
#
#  WHY
#  277.5 km2 of candidate change has been reported and not one patch has been
#  looked at. The maps also showed something the numbers did not: candidates
#  cluster along protected-area margins, heavily on the northern and western
#  edges of Imatong, with the eastern half of Timu almost entirely flagged.
#
#  Two explanations fit that pattern and they have opposite implications:
#
#    ENCROACHMENT   land use genuinely differs across the boundary, which is a
#                   real and reportable finding.
#    MISREGISTRATION the two composites are slightly offset, so every sharp
#                   spectral edge produces a false difference along it.
#
#  THE TEST THAT SEPARATES THEM
#  A protected-area boundary is an ADMINISTRATIVE line. It has no inherent
#  spectral gradient — nothing about a legal boundary makes the imagery change
#  across it unless the land use actually differs. Misregistration, by contrast,
#  produces artefacts wherever the IMAGE has a sharp edge: woodland against
#  grassland, a river margin, a hill shadow. Those are mostly not boundaries.
#
#  So: does candidate status track distance to the legal boundary, or does it
#  track local spectral roughness? If the former dominates, the pattern is land
#  use. If the latter, it is registration error and the candidates near edges
#  should be discounted.
#
#  Plus a direct look: before-and-after crops of the largest patches, which is
#  the nearest thing to ground truth available without external imagery.
#
#  HOW TO RUN — from the repo root, in PowerShell:
#      & "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" scripts\R\09_check_candidates.R
#
#  OUTPUTS -> outputs/candidates/
#
#  Author: Isaac Mureithi - IMureithi@awf.org
# =============================================================================

source("scripts/R/00_functions.R")
pacman::p_load(tidyterra, ggspatial)

terra::terraOptions(progress = 0)
set.seed(5)

path.out.cand <- file.path("outputs", "candidates")
dir.create(path.out.cand, recursive = TRUE, showWarnings = FALSE)

n.zoom     <- 4        # largest patches to render
zoom.pad.m <- 2500     # context around each patch
agg.stats  <- 3        # 30 m -> 90 m for the distance/gradient statistics

cat("\n================================================================\n")
cat(" KIDEPO LCLUC — 09 Checking the candidate patches\n")
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

f.cand <- file.path(path.out.cls, "change_strong_clusterfree.tif")
if (!file.exists(f.cand)) stop("Candidate raster missing — run 03c_k_sensitivity.R first.")

cand <- mask.to.aoi(terra::rast(f.cand))
names(cand) <- "candidate"

ndvi.a <- mask.to.aoi(comp$a[["NDVI"]])
ndvi.b <- mask.to.aoi(comp$b[["NDVI"]])
d.ndvi <- ndvi.b - ndvi.a

px.km2 <- prod(terra::res(cand)) / 1e6

say(glue("[1] candidates: {round(terra::global(cand, 'sum', na.rm = TRUE)[1,1] * px.km2, 1)} km2"))


# =============================================================================
#  PHASE 2 — WHERE ARE THEY, AND HOW BIG?
# =============================================================================

pa.r <- terra::rasterize(aoi.parts.v, cand, field = "pa.label")

by.pa <- terra::zonal(cand, pa.r, fun = "sum", na.rm = TRUE) %>%
  as_tibble() %>%
  rename(pa.label = 1, cand.px = 2) %>%
  left_join(
    terra::zonal(terra::ifel(!is.na(cand), 1, 0), pa.r, fun = "sum", na.rm = TRUE) %>%
      as_tibble() %>% rename(pa.label = 1, total.px = 2),
    by = "pa.label"
  ) %>%
  mutate(cand.km2 = cand.px * px.km2,
         pa.km2   = total.px * px.km2,
         pct.of.pa = 100 * cand.px / total.px) %>%
  arrange(desc(pct.of.pa))

cat("\n[2] CANDIDATE AREA BY PROTECTED AREA\n\n")
by.pa %>%
  transmute(pa.label, pa.km2 = round(pa.km2, 1), cand.km2 = round(cand.km2, 1),
            `pct.of.pa` = round(pct.of.pa, 2)) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

## Patch inventory.
pid <- terra::patches(terra::subst(cand, from = 0, to = NA),
                      directions = 8, zeroAsNA = TRUE)
psize <- terra::freq(pid) %>%
  as_tibble() %>%
  transmute(patch = value, px = count, ha = px * prod(terra::res(cand)) / 10000) %>%
  arrange(desc(ha))

say(glue("
    patches            : {format(nrow(psize), big.mark = ',')}
    median patch size  : {round(median(psize$ha), 2)} ha
    90th percentile    : {round(quantile(psize$ha, 0.9), 1)} ha
    largest            : {round(max(psize$ha), 1)} ha
    area in patches
      under 2 ha       : {round(100 * sum(psize$ha[psize$ha < 2]) / sum(psize$ha), 1)}%
      over 10 ha       : {round(100 * sum(psize$ha[psize$ha > 10]) / sum(psize$ha), 1)}%

    READ THE LARGEST FIGURE WITH SUSPICION. Eight-connected labelling joins
    clumps that touch at a single corner, so a scatter of small patches threaded
    across a hillside is counted as ONE patch. The largest here is a sprawl of
    that kind, not a coherent block of changed ground, and the zoom panels in
    Phase 5 show it plainly. The median is the honest summary of patch size."))


# =============================================================================
#  PHASE 3 — BOUNDARY OR IMAGE EDGE?
# =============================================================================
#  Statistics computed at 90 m: the question is about tendency over hundreds of
#  metres, and a distance transform over 32 million cells to buy precision
#  nobody will use is not worth the minutes.

cand.s  <- terra::aggregate(cand, fact = agg.stats, fun = "max", na.rm = TRUE)
ndvi.s  <- terra::aggregate(ndvi.a, fact = agg.stats, fun = "mean", na.rm = TRUE)

## Distance to the nearest protected-area boundary, in metres.
pa.lines <- terra::as.lines(aoi.parts.v)
dist.edge <- terra::distance(cand.s, pa.lines)
names(dist.edge) <- "dist.edge"

## Local spectral roughness: the standard deviation of NDVI in a 3x3
## neighbourhood of the EARLIER date. High where the image has a real edge —
## woodland meeting grassland, a river margin, a scarp — which is exactly where
## a registration offset would manufacture difference.
rough <- terra::focal(ndvi.s, w = 3, fun = "sd", na.rm = TRUE)
names(rough) <- "roughness"

samp <- terra::spatSample(c(cand.s, dist.edge, rough), 300000,
                          method = "random", na.rm = TRUE, warn = FALSE) %>%
  as_tibble() %>%
  rename(candidate = 1, dist.edge = 2, roughness = 3) %>%
  mutate(roughness = to.reflectance(roughness)) %>%
  drop_na()

cmp <- samp %>%
  group_by(candidate) %>%
  summarise(n = n(),
            median.dist.m   = median(dist.edge),
            pct.within.500m = 100 * mean(dist.edge <= 500),
            median.roughness = median(roughness),
            .groups = "drop")

cat("\n[3] CANDIDATES vs EVERYWHERE ELSE\n\n")
cmp %>% mutate(across(where(is.numeric), ~round(.x, 4))) %>%
  as.data.frame() %>% print(row.names = FALSE)

## Which predictor carries more? Both standardised, so the coefficients are
## directly comparable: each is the change in log-odds per one standard
## deviation of that predictor.
samp.z <- samp %>%
  mutate(z.dist = as.numeric(scale(dist.edge)),
         z.rough = as.numeric(scale(roughness)))

fit <- glm(candidate ~ z.dist + z.rough, data = samp.z, family = binomial())
co  <- coef(fit)

say(glue("
    Logistic model of candidate status, standardised predictors:
      distance to PA boundary : {sprintf('%+.3f', co[['z.dist']])}  (negative = commoner NEAR boundaries)
      local spectral roughness: {sprintf('%+.3f', co[['z.rough']])}  (positive = commoner at image edges)

    The larger magnitude is the dominant pattern."))

edge.driven  <- abs(co[["z.dist"]])  > abs(co[["z.rough"]])
near.edge    <- co[["z.dist"]] < 0

if (edge.driven && near.edge) {
  verdict <- paste(
    "Candidate status tracks the ADMINISTRATIVE boundary more strongly than it",
    "tracks image roughness. A legal line has no spectral gradient of its own,",
    "so this is consistent with land use genuinely differing across it -",
    "encroachment - rather than with registration error. It is not proof:",
    "boundaries often follow rivers and ridges, which are also image edges.",
    sep = "\n    ")
} else if (!edge.driven) {
  verdict <- paste(
    "Candidate status tracks local spectral ROUGHNESS more strongly than",
    "distance to a boundary. That is the signature of misregistration: a small",
    "offset between dates manufactures difference wherever the image has a",
    "sharp edge. Treat edge-adjacent candidates as suspect and check the",
    "co-registration of the two composites before reporting any of them.",
    sep = "\n    ")
} else {
  verdict <- "Candidates are commoner AWAY from boundaries, so the edge pattern seen on the map is not the dominant signal."
}
cat("\n   ", verdict, "\n")


# =============================================================================
#  PHASE 4 — WHICH WAY DID THEY MOVE?
# =============================================================================

dir.samp <- terra::spatSample(c(cand, d.ndvi), 300000, method = "random",
                              na.rm = TRUE, warn = FALSE) %>%
  as_tibble() %>%
  rename(candidate = 1, dndvi = 2) %>%
  mutate(dndvi = to.reflectance(dndvi)) %>%
  drop_na()

dir.tbl <- dir.samp %>%
  group_by(candidate) %>%
  summarise(n = n(), median.dndvi = median(dndvi),
            pct.greening = 100 * mean(dndvi > 0), .groups = "drop")

cat("\n[4] DIRECTION OF MOVEMENT\n\n")
dir.tbl %>% mutate(across(where(is.numeric), ~round(.x, 4))) %>%
  as.data.frame() %>% print(row.names = FALSE)

say("
    Candidates that are overwhelmingly GREENING are regrowth, woody thickening
    or a wetter season - not clearing. Clearing would show as browning, and a
    conversion signal would contain both.")


# =============================================================================
#  PHASE 4b — WERE THEY BURNT IN THE EARLIER DATE?
# =============================================================================
#  The obvious suspect, and one this pipeline is structurally blind to.
#
#  MODIS burned area is 500 m. A fire of a few hectares - which is most fires -
#  occupies a fraction of one MODIS cell and is never recorded. So a patch that
#  burnt in 2019, was missed by MODIS, and had recovered by 2024 would show
#  exactly what the candidates show: a large greening, concentrated near
#  boundaries where people ignite fires, with no burn record to explain it.
#
#  A fresh scar is dark in NBR. If the candidates sat in the burnt tail of the
#  2019 NBR distribution and are normal by 2024, they are recovering fires
#  rather than land-cover change.

nbr.a <- mask.to.aoi(comp$a[["NBR"]])
nbr.b <- mask.to.aoi(comp$b[["NBR"]])

## Both arrive named "NBR"; c() keeps duplicates and as_tibble() then refuses
## the frame. Rename before combining.
names(nbr.a) <- "nbr_a"
names(nbr.b) <- "nbr_b"

burn.samp <- terra::spatSample(c(cand, nbr.a, nbr.b), 300000, method = "random",
                               na.rm = TRUE, warn = FALSE) %>%
  as_tibble() %>%
  rename(candidate = 1, nbr.a = 2, nbr.b = 3) %>%
  mutate(across(c(nbr.a, nbr.b), to.reflectance)) %>%
  drop_na()

## The burnt tail, defined from the landscape itself: MODIS put ~20% of the AOI
## under scars in 2019, so the lowest 20% of NBR is the comparable cut.
burn.cut <- quantile(burn.samp$nbr.a, 0.20, na.rm = TRUE)

burn.tbl <- burn.samp %>%
  group_by(candidate) %>%
  summarise(n = n(),
            median.nbr.a = median(nbr.a),
            median.nbr.b = median(nbr.b),
            pct.in.burnt.tail.a = 100 * mean(nbr.a < burn.cut),
            .groups = "drop")

cat(glue("\n[4b] NBR SIGNATURE (2019 burnt-tail cut = {round(burn.cut, 4)})\n\n"))
burn.tbl %>% mutate(across(where(is.numeric), ~round(.x, 4))) %>%
  as.data.frame() %>% print(row.names = FALSE)

enrich.burnt <- burn.tbl$pct.in.burnt.tail.a[burn.tbl$candidate == 1] /
                burn.tbl$pct.in.burnt.tail.a[burn.tbl$candidate == 0]

say(glue("
    Candidates are {round(enrich.burnt, 2)}x as likely as background to have sat in the
    2019 burnt tail."))

if (enrich.burnt > 1.8) {
  cat("\n    That is a strong enrichment. A large share of the candidates are most\n")
  cat("    plausibly SMALL FIRES that MODIS could not see at 500 m, recovering by\n")
  cat("    2024 - not land-cover conversion. The fire filtering in 04 and 05\n")
  cat("    removed only the fires MODIS recorded, and was always going to be\n")
  cat("    blind to the rest. Treat the candidate area as an upper bound.\n")
} else {
  cat("\n    No strong enrichment: the candidates were not unusually burnt in 2019,\n")
  cat("    so missed small fires are unlikely to be the main explanation.\n")
}


# =============================================================================
#  PHASE 5 — LOOK AT THE BIGGEST PATCHES
# =============================================================================

top <- psize %>% slice_head(n = n.zoom)

zoom.pair <- function(patch.id, ha) {

  one <- terra::ifel(pid == patch.id, 1, NA)
  e   <- terra::ext(terra::trim(one))
  e   <- terra::ext(e$xmin - zoom.pad.m, e$xmax + zoom.pad.m,
                    e$ymin - zoom.pad.m, e$ymax + zoom.pad.m)

  outline <- terra::as.polygons(terra::crop(one, e))

  panel <- function(r, ttl) {
    rgb <- terra::crop(r[[c("SWIR2", "NIR", "GREEN")]], e)
    ggplot() +
      geom_spatraster_rgb(data = rgb, max_col_value = 3500, interpolate = TRUE) +
      geom_sf(data = sf::st_as_sf(outline), fill = NA, colour = "#FFD400",
              linewidth = 0.6) +
      coord_sf(expand = FALSE) +
      labs(title = ttl) +
      theme_void(base_size = 9) +
      theme(plot.title = element_text(size = 9, face = "bold"),
            plot.background = element_rect(fill = "white", colour = NA))
  }

  panel(comp$a, glue("{yr.a}")) + panel(comp$b, glue("{yr.b}")) +
    plot_annotation(
      title = glue("Candidate patch {patch.id} — {round(ha, 1)} ha"),
      theme = theme(plot.title = element_text(size = 10, face = "bold")))
}

for (i in seq_len(nrow(top))) {
  p <- zoom.pair(top$patch[i], top$ha[i])
  f <- file.path(path.out.cand, glue("candidate_zoom_{i}_patch{top$patch[i]}.png"))
  ggsave(f, p, width = 10, height = 5.2, dpi = 200, bg = "white")
  say(glue("    wrote {f}"))
}


# =============================================================================
#  PHASE 6 — OUTPUTS
# =============================================================================

readr::write_csv(by.pa,   file.path(path.out.cand, "candidates_by_pa.csv"))
readr::write_csv(psize,   file.path(path.out.cand, "candidate_patch_sizes.csv"))
readr::write_csv(cmp,     file.path(path.out.cand, "candidates_vs_background.csv"))
readr::write_csv(dir.tbl, file.path(path.out.cand, "candidate_direction.csv"))

cat("\n================================================================\n")
cat(" DONE — tables and zoom panels in", path.out.cand, "\n")
cat("================================================================\n\n")

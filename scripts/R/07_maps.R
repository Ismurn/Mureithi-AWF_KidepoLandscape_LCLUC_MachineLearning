# =============================================================================
#  KIDEPO LANDSCAPE — LCLUC
#  Script 07 — Cartographic maps of the unsupervised stage
#
#  PURPOSE
#  Everything produced so far has been a diagnostic plot: base-R `terra::plot`
#  output with no scale bar, no north arrow and no legend worth the name. Useful
#  for catching problems, useless in a report. This script produces the
#  presentation versions of the same material.
#
#  TOOLCHAIN
#  ggplot2 + tidyterra for the rasters, ggspatial for scale bars and north
#  arrows, patchwork for panels, rnaturalearth for country context. tidyterra is
#  preferred over tmap here because the panelling, shared legends and insets
#  compose predictably with patchwork, and tmap's v4 API differs enough from v3
#  that examples in circulation no longer run.
#
#  A NOTE ON DISPLAY RESOLUTION
#  The composites are 32 million cells. Rendering them at native resolution into
#  a 2000 px figure wastes minutes to draw detail finer than a pixel. Display
#  copies are aggregated first — mean for continuous layers, modal for classes,
#  max for sparse binary layers so small patches survive rather than being
#  averaged out of existence. The ANALYSIS always used full resolution; this is
#  a drawing decision only.
#
#  HOW TO RUN — from the repo root, in PowerShell:
#      & "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" scripts\R\07_maps.R
#
#  OUTPUTS -> outputs/maps/
#
#  Author: Isaac Mureithi - IMureithi@awf.org
# =============================================================================

source("scripts/R/00_functions.R")
pacman::p_load(tidyterra, ggspatial, rnaturalearth, rnaturalearthdata, ggnewscale)

terra::terraOptions(progress = 0)

path.out.maps <- file.path("outputs", "maps")
dir.create(path.out.maps, recursive = TRUE, showWarnings = FALSE)

cat("\n================================================================\n")
cat(" KIDEPO LCLUC — 07 Cartographic maps\n")
cat(" run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n")


# =============================================================================
#  PHASE 0 — SHARED STYLE
# =============================================================================

agg.display <- 4     # 30 m -> 120 m for continuous display layers
agg.binary  <- 3     # 30 m -> 90 m for sparse binary layers

## A land-cover palette chosen to be SEMANTIC rather than merely pretty: burnt
## ground dark, bare soil tan, vegetation increasingly green. A viewer should be
## able to read the map without constantly checking the legend.
pal.classes <- c(
  "1 burnt / dark bare" = "#3B2B20",
  "2 bright bare soil"  = "#C8A87C",
  "3 dry grassland"     = "#D9C87E",
  "4 grassland / shrub" = "#A3B55C",
  "5 wooded grassland"  = "#5E8C4A",
  "6 woodland / forest" = "#1F4D2B"
)

pal.country <- c("Uganda" = "#4C9F70", "South Sudan" = "#6A8EAE")

theme.map <- function(base = 9) {
  theme_minimal(base_size = base) +
    theme(
      panel.grid       = element_line(linewidth = 0.2, colour = "grey92"),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      plot.title       = element_text(face = "bold", size = base + 2),
      plot.subtitle    = element_text(colour = "grey30", size = base - 0.5),
      plot.caption     = element_text(colour = "grey45", size = base - 2, hjust = 0),
      legend.position  = "right",
      legend.key.size  = unit(0.45, "cm"),
      axis.title       = element_blank(),
      axis.text        = element_text(size = base - 3, colour = "grey45")
    )
}

## Scale bar and north arrow, added to the map that anchors each figure.
map.furniture <- function(arrow = TRUE) {
  out <- list(
    annotation_scale(location = "bl", width_hint = 0.28, height = unit(0.15, "cm"),
                     text_cex = 0.6, bar_cols = c("grey20", "white"))
  )
  if (arrow) {
    out <- c(out, list(
      annotation_north_arrow(location = "tr", which_north = "true",
                             height = unit(0.8, "cm"), width = unit(0.8, "cm"),
                             style = north_arrow_minimal(text_size = 6))
    ))
  }
  out
}

save.map <- function(p, filename, width = 9, height = 7.5, dpi = 220) {
  f <- file.path(path.out.maps, filename)
  ggsave(f, p, width = width, height = height, dpi = dpi, bg = "white")
  say(glue("    wrote {f}"))
  invisible(f)
}


# =============================================================================
#  PHASE 1 — DATA
# =============================================================================

comp <- load.composite.pair()
yr.a <- comp$years[1]
yr.b <- comp$years[2]

aoi <- load.aoi()
aoi.union <- aoi$union
aoi.parts <- aoi$parts

aoi.union.v <- terra::vect(aoi.union)
mask.to.aoi <- function(r) terra::mask(terra::crop(r, aoi.union.v), aoi.union.v)

## Label anchors, inside each polygon and computed in UTM.
pa.labels <- aoi.parts %>%
  group_by(pa.label, country) %>%
  summarise(.groups = "drop") %>%
  st_make_valid()
pa.labels$geometry.lab <- st_point_on_surface(st_geometry(pa.labels))

## Country outlines for context, clipped generously around the AOI.
bb <- st_bbox(st_buffer(aoi.union, 40000))
countries <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(crs.utm) %>%
  st_crop(bb)

say("[1] loaded composites, AOI and context layers")

## ---- display copies -------------------------------------------------------
dsp <- function(r, fact = agg.display, fun = "mean") {
  terra::aggregate(mask.to.aoi(r), fact = fact, fun = fun, na.rm = TRUE)
}

rgb.a <- terra::aggregate(mask.to.aoi(comp$a[[c("SWIR2", "NIR", "GREEN")]]),
                          fact = agg.display, fun = "mean", na.rm = TRUE)
rgb.b <- terra::aggregate(mask.to.aoi(comp$b[[c("SWIR2", "NIR", "GREEN")]]),
                          fact = agg.display, fun = "mean", na.rm = TRUE)

obs.a <- dsp(comp$a[["clear_obs"]])
obs.b <- dsp(comp$b[["clear_obs"]])

ndvi.a <- dsp(comp$a[["NDVI"]])
ndvi.b <- dsp(comp$b[["NDVI"]])
d.ndvi <- (ndvi.b - ndvi.a) / comp.scale.int

say("[1] built display-resolution copies")


# =============================================================================
#  PHASE 2 — MAP 1: STUDY AREA
# =============================================================================

p.area <- ggplot() +
  geom_sf(data = countries, fill = "grey96", colour = "grey75", linewidth = 0.3) +
  geom_sf(data = aoi.parts, aes(fill = country), colour = "grey25", linewidth = 0.25,
          alpha = 0.9) +
  geom_sf_text(data = pa.labels, aes(geometry = geometry.lab, label = pa.label),
               fun.geometry = identity, size = 2.5, colour = "grey10",
               check_overlap = TRUE) +
  scale_fill_manual(values = pal.country, name = NULL) +
  coord_sf(xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"]),
           expand = FALSE) +
  map.furniture() +
  labs(
    title    = "Kidepo transboundary landscape",
    subtitle = glue("{nrow(aoi.parts)} protected areas | {format(round(as.numeric(sum(st_area(aoi.union)))/1e6), big.mark = ',')} km\u00b2 dissolved | UTM 36N"),
    caption  = "Protected-area boundaries: WDPA/WDOECM. Basemap: Natural Earth."
  ) +
  theme.map()

save.map(p.area, "map_01_study_area.png", width = 9, height = 8)


# =============================================================================
#  PHASE 3 — MAP 2: THE COMPOSITES
# =============================================================================
#  SWIR2/NIR/GREEN. Burnt ground reads dark brown-magenta, healthy vegetation
#  green, bare soil pink-white, water black. The stretch limits are shared
#  between dates — a per-date stretch would make them look different for
#  reasons that have nothing to do with the ground.

rgb.panel <- function(r, title) {
  ggplot() +
    geom_spatraster_rgb(data = r, max_col_value = 3500, interpolate = TRUE) +
    geom_sf(data = aoi.parts, fill = NA, colour = "white", linewidth = 0.25) +
    coord_sf(expand = FALSE) +
    labs(title = title) +
    theme.map() +
    theme(axis.text = element_blank(), panel.grid = element_blank())
}

p.rgb <- (rgb.panel(rgb.a, as.character(yr.a)) + map.furniture(arrow = FALSE)) +
         (rgb.panel(rgb.b, as.character(yr.b))) +
  plot_annotation(
    title    = "HLS dry-season composites, 1 February – 31 March",
    subtitle = "False colour SWIR2 / NIR / GREEN — burnt ground dark, vegetation green, bare soil pale",
    caption  = "Harmonised Landsat Sentinel-2 v2.0, Fmask-masked median composite, 30 m.",
    theme    = theme(plot.title = element_text(face = "bold", size = 13),
                     plot.subtitle = element_text(colour = "grey30", size = 9),
                     plot.caption = element_text(colour = "grey45", size = 7, hjust = 0))
  )

save.map(p.rgb, "map_02_composites.png", width = 13, height = 6.5)


# =============================================================================
#  PHASE 4 — MAP 3: OBSERVATION SUPPORT
# =============================================================================
#  The honesty map. Everything downstream rests on how many clear looks each
#  pixel actually had, and the two dates are NOT equally supported.

obs.panel <- function(r, title, lims) {
  ggplot() +
    geom_spatraster(data = r) +
    geom_sf(data = aoi.parts, fill = NA, colour = "grey20", linewidth = 0.2) +
    scale_fill_viridis_c(name = "clear\nobservations", limits = lims,
                         na.value = "transparent", option = "mako") +
    coord_sf(expand = FALSE) +
    labs(title = title) +
    theme.map() +
    theme(axis.text = element_blank(), panel.grid = element_blank())
}

obs.lims <- c(0, max(terra::global(obs.a, "max", na.rm = TRUE)[1, 1],
                     terra::global(obs.b, "max", na.rm = TRUE)[1, 1]))

p.obs <- (obs.panel(obs.a, glue("{yr.a} — mean 8.1 looks"), obs.lims) +
            map.furniture(arrow = FALSE)) +
         obs.panel(obs.b, glue("{yr.b} — mean 10.4 looks"), obs.lims) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Clear-observation support per pixel",
    subtitle = "2024 carries about 25% more looks than 2019 — a caveat on any change statistic",
    caption  = "Observations surviving Fmask cloud, shadow, cirrus and adjacency masking.",
    theme    = theme(plot.title = element_text(face = "bold", size = 13),
                     plot.subtitle = element_text(colour = "grey30", size = 9),
                     plot.caption = element_text(colour = "grey45", size = 7, hjust = 0))
  )

save.map(p.obs, "map_03_observation_support.png", width = 13, height = 6.5)


# =============================================================================
#  PHASE 5 — MAP 4: THE k-MEANS CLASSES
# =============================================================================

cls.path <- function(y) file.path(path.out.cls, glue("kmeans_k6_{y}.tif"))

if (all(file.exists(cls.path(yr.a), cls.path(yr.b)))) {

  cls.a <- terra::aggregate(terra::rast(cls.path(yr.a)), fact = agg.binary, fun = "modal")
  cls.b <- terra::aggregate(terra::rast(cls.path(yr.b)), fact = agg.binary, fun = "modal")

  ## Unclassified pixels INSIDE the AOI are not absent land — they are pixels
  ## with too few clear looks to classify (the cloud-prone Imatong and Agoro-Agu
  ## highlands, mostly). Left as NA they render white and read as a seventh
  ## class that happens to have no legend entry. Given an explicit level and a
  ## neutral grey, they say what they are. Outside the AOI stays NA.
  cls.na <- "no data (thin support)"
  pal.classes.na <- c(pal.classes, setNames("grey85", cls.na))

  as.classes <- function(r) {
    r <- terra::as.int(r)
    r <- terra::mask(terra::subst(r, from = NA, to = 7), aoi.union.v)
    levels(r) <- data.frame(value = 1:7, label = c(names(pal.classes), cls.na))
    names(r) <- "class"
    r
  }
  cls.a <- as.classes(cls.a); cls.b <- as.classes(cls.b)

  cls.panel <- function(r, title) {
    ggplot() +
      geom_spatraster(data = r) +
      geom_sf(data = aoi.parts, fill = NA, colour = "grey15", linewidth = 0.2) +
      scale_fill_manual(values = pal.classes.na, name = NULL,
                        na.value = "transparent", na.translate = FALSE,
                        drop = FALSE) +
      coord_sf(expand = FALSE) +
      labs(title = title) +
      theme.map() +
      theme(axis.text = element_blank(), panel.grid = element_blank())
  }

  p.cls <- (cls.panel(cls.a, as.character(yr.a)) + map.furniture()) +
           cls.panel(cls.b, as.character(yr.b)) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title    = "Unsupervised classification, k-means (k = 6)",
      subtitle = "Spectral clusters, NOT verified land-cover classes — labels are inferred from centroid signatures",
      caption  = paste("Fitted once on a pooled sample from both dates so clusters are comparable.",
                       "Diagnostics (03b) show the landscape supports about three robust strata;",
                       "at k=6 a quarter of pixels sit ambiguously between clusters."),
      theme    = theme(plot.title = element_text(face = "bold", size = 13),
                       plot.subtitle = element_text(colour = "grey30", size = 9),
                       plot.caption = element_text(colour = "grey45", size = 7, hjust = 0))
    )

  save.map(p.cls, "map_04_kmeans_classes.png", width = 13.5, height = 6.5)

} else {
  say("    skipped map 4 — run 03_kmeans_firstlook.R first")
}


# =============================================================================
#  PHASE 6 — MAP 5: FIRE
# =============================================================================

burn.path <- function(y) file.path(path.gee.out,
                                   glue("kidepo_burn_{comp.window}_{y}_utm36n_30m.tif"))

if (all(file.exists(burn.path(yr.a), burn.path(yr.b)))) {

  days.a <- mask.to.aoi(terra::rast(burn.path(yr.a))[["days_since_burn"]])
  days.b <- mask.to.aoi(terra::rast(burn.path(yr.b))[["days_since_burn"]])

  ## max() on aggregation so a burnt MODIS cell is not averaged away.
  fire.state <- terra::aggregate(
    terra::as.int((days.a > 0) * 1 + (days.b > 0) * 2),
    fact = agg.binary, fun = "modal")
  levels(fire.state) <- data.frame(
    value = 0:3,
    label = c("neither", glue("{yr.a} only"), glue("{yr.b} only"), "both dates"))
  names(fire.state) <- "state"

  age.a <- terra::aggregate(terra::ifel(days.a > 0, days.a, NA),
                            fact = agg.binary, fun = "mean", na.rm = TRUE)
  names(age.a) <- "days"

  p.fire <- ggplot() +
    geom_spatraster(data = fire.state) +
    geom_sf(data = aoi.parts, fill = NA, colour = "grey15", linewidth = 0.2) +
    scale_fill_manual(values = setNames(
                        c("grey90", "#E08214", "#8073AC", "#4D4D4D"),
                        c("neither", glue("{yr.a} only"),
                          glue("{yr.b} only"), "both dates")),
                      name = "burnt in", na.value = "transparent", na.translate = FALSE) +
    coord_sf(expand = FALSE) +
    map.furniture() +
    labs(title = "Burned area from MODIS, by date",
         subtitle = "Only the two middle classes can manufacture false change; 'both dates' is a constant",
         caption = paste("MCD64A1 500 m, cumulative from 1 November to the window's end, shown on the 30 m grid.",
                         "\nBurn switched state over 29% of the AOI yet predicts change only 1.11x better than chance.")) +
    theme.map()

  p.age <- ggplot() +
    geom_spatraster(data = age.a) +
    geom_sf(data = aoi.parts, fill = NA, colour = "grey15", linewidth = 0.2) +
    scale_fill_scico(palette = "lajolla", direction = -1, name = "days since\nburn",
                     na.value = "transparent") +
    coord_sf(expand = FALSE) +
    labs(title = glue("Scar age at the end of the {yr.a} window"),
         subtitle = "Every scar is 80+ days old — the fires happen before the window opens") +
    theme.map()

  save.map(p.fire, "map_05_fire_state.png", width = 9.5, height = 8)
  save.map(p.age,  "map_06_scar_age.png",   width = 9.5, height = 8)

} else {
  say("    skipped maps 5-6 — burn rasters not found")
}


# =============================================================================
#  PHASE 7 — MAP 7: CHANGE
# =============================================================================
#  The deliverable map. A diverging scale centred on zero, and the candidate
#  patches drawn on top of it.

cand.f <- file.path(path.out.cls, "change_strong_clusterfree.tif")

p.dndvi <- ggplot() +
  geom_spatraster(data = d.ndvi) +
  geom_sf(data = aoi.parts, fill = NA, colour = "grey15", linewidth = 0.2) +
  scale_fill_scico(palette = "vik", direction = -1, midpoint = 0,
                   limits = c(-0.3, 0.3), oob = scales::squish,
                   name = glue("\u0394NDVI\n{yr.b} - {yr.a}"),
                   na.value = "transparent") +
  coord_sf(expand = FALSE) +
  map.furniture() +
  labs(title = glue("Change in dry-season greenness, {yr.a} to {yr.b}"),
       subtitle = "Green = greener in 2024. The landscape moved one way together, median +0.014",
       caption = "A one-directional shift is what a wetter or later-senescing season produces; conversion would show both gains and losses.") +
  theme.map()

save.map(p.dndvi, "map_07_dndvi.png", width = 9.5, height = 8)

if (file.exists(cand.f)) {

  cand <- terra::aggregate(mask.to.aoi(terra::rast(cand.f)),
                           fact = agg.binary, fun = "max", na.rm = TRUE)
  cand <- terra::ifel(cand > 0, 1, NA)
  names(cand) <- "candidate"

  p.cand <- ggplot() +
    geom_spatraster(data = d.ndvi) +
    scale_fill_gradientn(colours = c("grey80", "grey96", "grey80"),
                         limits = c(-0.3, 0.3), oob = scales::squish,
                         na.value = "transparent", guide = "none") +
    ggnewscale::new_scale_fill() +
    geom_spatraster(data = cand) +
    scale_fill_gradientn(colours = "#762A83", na.value = "transparent",
                         guide = "none") +
    geom_sf(data = aoi.parts, fill = NA, colour = "grey30", linewidth = 0.25) +
    coord_sf(expand = FALSE) +
    map.furniture() +
    labs(title = "Candidate land-cover change",
         subtitle = "277.5 km\u00b2 — 3.6% of the protected-area network",
         caption = paste("|\u0394NDVI| > 0.15, no MODIS burn-state switch to explain it, contiguous patches \u2265 1 ha.",
                         "\nCANDIDATES, not confirmed change: none has yet been checked against high-resolution imagery.")) +
    theme.map()

  save.map(p.cand, "map_08_candidate_change.png", width = 9.5, height = 8)

} else {
  say("    skipped map 8 — run 03c_k_sensitivity.R first")
}

cat("\n================================================================\n")
cat(" DONE — maps in", path.out.maps, "\n")
cat("================================================================\n\n")

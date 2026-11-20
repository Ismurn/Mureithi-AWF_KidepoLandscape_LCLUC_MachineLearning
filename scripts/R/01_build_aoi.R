#  KIDEPO LANDSCAPE - LAND COVER / LAND USE CHANGE (LCLUC)

#  Script 01 - Build the Area Of Interest (AOI) from WDPA protected-area polygons

#  PURPOSE
#  Stage 0 of the hybrid GEE -> R pipeline. This script reads the per-protected-area
#  WDPA/WDOECM shapefiles held on OneDrive, audits them against the protected-area
#  register recorded in CLAUDE.md, dissolves them into a single AOI, and writes the
#  geometries that the Google Earth Engine stage needs:
#
#    (a) the individual PA polygons          -> per-PA reporting later in R
#    (b) the dissolved union of all PAs      -> the CLASSIFICATION / REPORTING extent
#    (c) a 500 m buffered union              -> composite slightly wider than (b) so
#                                               focal/edge operations have valid pixels
#    (d) the bounding box over all PAs       -> the COMPOSITING extent in GEE, and the
#                                               later "full landscape matrix" extension
#
#  Everything written for GEE is in EPSG:4326 (geographic), because Earth Engine
#  ingests and stores vector assets in WGS84. All AREA arithmetic is done in
#  EPSG:32636 (UTM Zone 36N), the project's working projected CRS.
#
#  OUTPUTS  ->  outputs/aoi/   (gitignored; fully regenerable by re-running this script)
#
#  HOW TO RUN
#  From the repo root, in PowerShell:
#      & "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" scripts\R\01_build_aoi.R
#  (Do NOT run this via Git Bash/MSYS - the MSYS `mingw64` DLLs on PATH shadow the
#   GDAL libraries that `sf` links against and R segfaults on st_read().)
#
#  Author: Isaac Mureithi - IMureithi@awf.org

#  PHASE 0 - SETUP: packages, paths, constants

## 0.1 Packages 
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(sf, tidyverse, glue, jsonlite)

sf::sf_use_s2(TRUE)  # spherical geometry for the WGS84 steps; areas are done in UTM

## 0.1b Console helper 
## glue() strips one trailing newline, which otherwise runs consecutive progress
## messages together on one line. say() puts it back.
say <- function(x) cat(x, "\n", sep = "")

## 0.2 Paths 
## The raw data live on OneDrive, OUTSIDE the repo. Note the apostrophe in
## "Uganda's" and the spaces - always build these with file.path(), never paste
## unquoted strings into a shell.
path.data.root <- "C:/Users/IsaacM/OneDrive - African Wildlife Foundation/R WORK/Uganda's Kidepo Region"
path.wdpa      <- file.path(path.data.root, "WDPA")

## Code runs from the repo root; derived artefacts stay inside the repo but gitignored.
path.out.aoi <- file.path("outputs", "aoi")
dir.create(path.out.aoi, recursive = TRUE, showWarnings = FALSE)

## 0.3 Coordinate reference systems 
crs.geo <- 4326   # WGS84 geographic - what GEE wants for vector assets
crs.utm <- 32636  # UTM Zone 36N     - project working CRS for all area/distance maths

## 0.4 Analysis constants 
buffer.m <- 500   # metres; buffer applied to the dissolved union for edge safety

## 0.5 Fail loudly if the data folder is not reachable 
## OneDrive files can be cloud-only placeholders. If a file is a placeholder,
## st_read() fails as though it were missing, so check up front.
if (!dir.exists(path.wdpa)) {
  stop(glue("WDPA folder not found:\n  {path.wdpa}\n",
            "Check the OneDrive path, and make sure the folder is downloaded locally ",
            "(right-click > 'Always keep on this device')."))
}

cat("\n================================================================\n")
cat(" KIDEPO LCLUC - 01 Build AOI\n")
cat(" run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n")


#  PHASE 1 - THE PROTECTED-AREA REGISTER (expected inventory)

#  Reported sizes transcribed from CLAUDE.md. This is the EXPECTED inventory; it is
#  deliberately kept separate from the shapefiles so that a PA missing from disk is
#  reported rather than silently dropped. `folder` is the subfolder under WDPA/.
#  `size.km2.reported` is the size as reported in the register (NA where not reported).

pa.register <- tribble(
  ~folder,                             ~pa.label,          ~country,      ~size.km2.reported,
  "KidepoValleyNationalPark_Uganda",   "Kidepo Valley",    "Uganda",      1430.00,
  "KidepoNationalReserve_SouthSudan",  "Kidepo",           "South Sudan", 1200.00,
  "Imatong_SouthSudan",                "Imatong",          "South Sudan", 1100.00,
  "Karenga_Uganda",                    "Karenga",          "Uganda",       956.00,
  "Nyangea-Napore_Uganda",             "Nyangea - Napore", "Uganda",       423.07,
  "Agoro-Agu_Uganda",                  "Agoro - Agu",      "Uganda",       263.52,
  "Lopeichubei_Uganda",                "Lopeichubei",      "Uganda",       166.64,
  "Timu_Uganda",                       "Timu",             "Uganda",       122.11,
  "Rom_Uganda",                        "Rom",              "Uganda",       108.98,
  "Lwala_Uganda",                      "Lwala",            "Uganda",        58.92,
  "Lotim-Puta_Uganda",                 "Lotim - Puta",     "Uganda",        18.97,
  "Lomej_Uganda",                      "Lomej",            "Uganda",         7.66,
  "Zulia_Uganda",                      "Zulia",            "Uganda",          NA
)

say(glue("\n[1] Register lists {nrow(pa.register)} protected areas ",
         "({sum(!is.na(pa.register$size.km2.reported))} with a reported size).\n"))


#  PHASE 2 - DISCOVER AND READ THE SHAPEFILES

## 2.1 Find every polygon shapefile under WDPA/ 
pa.files <- list.files(path.wdpa, pattern = "\\.shp$", recursive = TRUE,
                       full.names = TRUE, ignore.case = TRUE)

if (length(pa.files) == 0) {
  stop(glue("No .shp files found under {path.wdpa}"))
}

## The folder name (not the cryptic WDPA_WDOECM_* filename) is our join key.
pa.files.tbl <- tibble(
  path   = pa.files,
  folder = basename(dirname(pa.files)),
  file   = basename(pa.files)
)

say(glue("[2] Found {nrow(pa.files.tbl)} shapefile(s) in {n_distinct(pa.files.tbl$folder)} folder(s).\n"))

## 2.2 Read each one, keeping a harmonised attribute subset 
## These data are in the newer WDPA/WDOECM schema: SITE_ID / SITE_PID rather than
## the classic WDPAID / WDPA_PID. Select defensively so a schema change surfaces as
## an explicit warning instead of a crash.
fields.wanted <- c("SITE_ID", "SITE_PID", "NAME", "NAME_ENG", "DESIG_ENG",
                   "DESIG_TYPE", "IUCN_CAT", "ISO3", "STATUS", "STATUS_YR",
                   "REP_AREA", "GIS_AREA", "GOV_TYPE", "MANG_AUTH")

read.one.pa <- function(path, folder) {

  x <- st_read(path, quiet = TRUE)

  ## Flag (don't fail on) any expected field that this file lacks.
  fields.missing <- setdiff(fields.wanted, names(x))
  if (length(fields.missing) > 0) {
    warning(glue("{folder}: missing attribute(s): {paste(fields.missing, collapse = ', ')}"),
            call. = FALSE)
  }

  x %>%
    select(any_of(fields.wanted)) %>%
    mutate(folder = folder, .before = 1) %>%
    st_transform(crs.geo)
}

pa.list <- pmap(list(pa.files.tbl$path, pa.files.tbl$folder), read.one.pa)
pa.raw  <- bind_rows(pa.list)

## Guarantee every expected column exists (as NA if truly absent from the source),
## so the downstream summarise() can name them unconditionally.
fields.absent <- setdiff(fields.wanted, names(pa.raw))
for (f in fields.absent) pa.raw[[f]] <- NA

say(glue("    Read {nrow(pa.raw)} polygon feature(s); CRS = EPSG:{st_crs(pa.raw)$epsg}.\n"))

## 2.3 Repair geometry 
## WDPA polygons quite often carry self-intersections / ring-order problems that
## break st_union() later. Fix now and report whether anything actually needed it.
n.invalid <- sum(!st_is_valid(pa.raw))
pa.clean  <- st_make_valid(pa.raw)
say(glue("    Geometry validity: {n.invalid} of {nrow(pa.raw)} feature(s) needed repair.\n"))

## Drop anything that is not a (multi)polygon after repair. (Base-style subsetting:
## `st_geometry_type()` needs the object itself, which dplyr::filter() does not expose.)
keep.poly <- st_geometry_type(pa.clean) %in% c("POLYGON", "MULTIPOLYGON")
if (any(!keep.poly)) {
  say(glue("    dropped {sum(!keep.poly)} non-polygon feature(s) after repair.\n"))
}
pa.clean <- pa.clean[keep.poly, ]


#  PHASE 3 - AUDIT: what is on disk vs. what the register expects

## 3.1 Per-PA computed area, in UTM 36N 
pa.utm <- st_transform(pa.clean, crs.utm)

pa.areas <- pa.utm %>%
  mutate(size.km2.computed = as.numeric(st_area(pa.utm)) / 1e6) %>%
  st_drop_geometry() %>%
  group_by(folder) %>%
  summarise(
    n.features        = n(),
    name.wdpa         = first(NAME),
    desig             = first(DESIG_ENG),
    iucn              = first(IUCN_CAT),
    iso3              = first(ISO3),
    status.yr         = first(STATUS_YR),
    size.km2.repfield = sum(REP_AREA),
    size.km2.gisfield = sum(GIS_AREA),
    size.km2.computed = sum(size.km2.computed),
    .groups = "drop"
  )

## 3.2 Join register <-> disk, and compute the discrepancy 
## Two different comparisons, which must not be conflated:
##
##   diff.gis.pct  computed area  vs  WDPA's own GIS_AREA attribute.
##                 This is OUR check: it tests that the read, the repair and the
##                 UTM 36N area computation are right. It should be ~0 everywhere.
##
##   diff.rep.pct  WDPA's GIS_AREA  vs  the register's reported size (= REP_AREA,
##                 the gazetted figure from the founding legislation). Divergence
##                 here is a property of the SOURCE DATA, not a bug: gazetted and
##                 digitised extents often differ, sometimes by a lot.
pa.audit <- pa.register %>%
  full_join(pa.areas, by = "folder") %>%
  mutate(
    on.disk      = !is.na(n.features),
    diff.gis.pct = 100 * (size.km2.computed - size.km2.gisfield) / size.km2.gisfield,
    diff.rep.pct = 100 * (size.km2.gisfield - size.km2.reported)  / size.km2.reported
  ) %>%
  arrange(desc(coalesce(size.km2.computed, size.km2.reported)))

cat("\n[3] PROTECTED-AREA AUDIT (areas in km2, computed in UTM 36N)\n")
cat("    gazetted = register/REP_AREA | wdpa.gis = GIS_AREA attribute | computed = this script\n\n")
pa.audit %>%
  transmute(
    pa.label,
    country,
    on.disk,
    gazetted  = round(size.km2.reported, 2),
    wdpa.gis  = round(size.km2.gisfield, 2),
    computed  = round(size.km2.computed, 2),
    `vs.gis_%` = round(diff.gis.pct, 2),
    `vs.gaz_%` = round(diff.rep.pct, 1),
    iucn,
    status.yr
  ) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

## 3.3 Report register entries with no shapefile 
pa.absent <- pa.audit %>% filter(!on.disk)
if (nrow(pa.absent) > 0) {
  cat("\n    !! IN THE REGISTER BUT NOT ON DISK:\n")
  pa.absent %>%
    transmute(pa.label, country, reported.km2 = size.km2.reported) %>%
    as.data.frame() %>%
    print(row.names = FALSE)
  cat("       These PAs are EXCLUDED from the AOI. Download their WDPA polygons into\n",
      "      a new subfolder under WDPA/ and re-run to include them.\n", sep = "")
}

## 3.4 Report shapefiles not in the register 
pa.unexpected <- pa.audit %>% filter(on.disk, is.na(pa.label))
if (nrow(pa.unexpected) > 0) {
  cat("\n    !! ON DISK BUT NOT IN THE REGISTER (still included in the AOI):\n")
  print(as.data.frame(pa.unexpected[, c("folder", "name.wdpa", "size.km2.computed")]),
        row.names = FALSE)
}

## 3.5 OUR OWN CHECK: computed area must agree with WDPA's GIS_AREA 
## A real failure mode. >1% means the geometry read, the validity repair or the
## projection is wrong, and the AOI cannot be trusted.
pa.geom.bad <- pa.audit %>% filter(on.disk, !is.na(diff.gis.pct), abs(diff.gis.pct) > 1)
if (nrow(pa.geom.bad) > 0) {
  cat("\n    !! COMPUTED AREA DISAGREES WITH WDPA's GIS_AREA BY >1% - investigate:\n")
  pa.geom.bad %>%
    transmute(pa.label, wdpa.gis = round(size.km2.gisfield, 2),
              computed = round(size.km2.computed, 2), `diff_%` = round(diff.gis.pct, 2)) %>%
    as.data.frame() %>%
    print(row.names = FALSE)
} else {
  say(glue("\n    OK: computed areas agree with WDPA GIS_AREA to within ",
           "{round(max(abs(pa.audit$diff.gis.pct), na.rm = TRUE), 2)}% for all PAs.\n\n"))
}

## 3.6 SOURCE-DATA NOTE: gazetted vs digitised extent 
## Not a bug - but it matters for how PA areas are quoted in the write-up, so it
## is reported explicitly rather than left for a reader to trip over.
pa.gaz.gap <- pa.audit %>% filter(on.disk, !is.na(diff.rep.pct), abs(diff.rep.pct) > 10)
if (nrow(pa.gaz.gap) > 0) {
  cat("\n    NOTE: digitised boundary differs from the gazetted size by >10%:\n")
  pa.gaz.gap %>%
    transmute(pa.label, gazetted = round(size.km2.reported, 2),
              wdpa.gis = round(size.km2.gisfield, 2), `diff_%` = round(diff.rep.pct, 1)) %>%
    as.data.frame() %>%
    print(row.names = FALSE)
  cat("          The MAPPED extent is what the classification will actually cover, so\n",
      "         quote the computed areas in results, and cite the gazetted figure only\n",
      "         as the legal designation.\n", sep = "")
}

## 3.7 OVERLAP BETWEEN PA POLYGONS 
## WDPA designations are frequently stacked: a wildlife reserve can sit partly
## inside a national park, and both are legitimate separate records. This does NOT
## affect the dissolved AOI (st_union collapses overlaps), but it DOES matter for
## per-PA zonal statistics later - a pixel in an overlap would be counted once for
## each PA, so PA areas will not sum to the landscape total. Quantify it now.
pa.utm.lab <- pa.utm %>%
  left_join(pa.register %>% select(folder, pa.label), by = "folder") %>%
  mutate(pa.label = coalesce(pa.label, folder))

## Upper triangle only - each unordered pair tested once.
hits <- st_intersects(pa.utm.lab, sparse = FALSE) &
        upper.tri(matrix(NA, nrow(pa.utm.lab), nrow(pa.utm.lab)))
pairs.ij <- which(hits, arr.ind = TRUE)

pa.overlaps <- map_dfr(seq_len(nrow(pairs.ij)), function(k) {
  i <- pairs.ij[k, 1]; j <- pairs.ij[k, 2]
  ## Polygons that merely share a boundary intersect as a LINESTRING, whose area
  ## is 0 - so the area filter below separates true overlap from simple adjacency.
  g <- suppressWarnings(
    st_intersection(st_geometry(pa.utm.lab)[i], st_geometry(pa.utm.lab)[j])
  )
  tibble(
    pa.a        = pa.utm.lab$pa.label[i],
    pa.b        = pa.utm.lab$pa.label[j],
    overlap.km2 = if (length(g) == 0) 0 else as.numeric(sum(st_area(g))) / 1e6
  )
}) %>%
  filter(overlap.km2 > 0.01) %>%
  arrange(desc(overlap.km2))

area.sum.individual <- sum(pa.areas$size.km2.computed)

if (nrow(pa.overlaps) > 0) {
  cat("\n    OVERLAPPING PA POLYGONS (areas double-counted in per-PA reporting):\n")
  pa.overlaps %>%
    mutate(overlap.km2 = round(overlap.km2, 1)) %>%
    as.data.frame() %>%
    print(row.names = FALSE)
} else {
  cat("\n    No PA polygons overlap.\n")
}

## 4.1 Attach tidy labels to the individual PAs 
## This is the layer used later for per-PA zonal statistics (exactextractr).
aoi.pa <- pa.clean %>%
  left_join(pa.register %>% select(folder, pa.label, country), by = "folder") %>%
  mutate(pa.label = coalesce(pa.label, folder)) %>%
  relocate(pa.label, country, .after = folder)

## 4.2 Dissolved union - the classification / reporting extent 
## Union in UTM (planar) then transform back: avoids s2 edge cases where polygons
## merely touch, and keeps the buffer in Phase 4.3 in true metres.
aoi.union.utm <- pa.utm %>%
  st_geometry() %>%
  st_union() %>%
  st_make_valid()

aoi.union <- aoi.union.utm %>%
  st_transform(crs.geo) %>%
  st_sf(geometry = ., aoi = "kidepo.pa.union")

## 4.3 Buffered union - the safe compositing/clip extent 
aoi.union.buf <- aoi.union.utm %>%
  st_buffer(buffer.m) %>%
  st_union() %>%
  st_transform(crs.geo) %>%
  st_sf(geometry = ., aoi = glue("kidepo.pa.union.buffer{buffer.m}m"))

## 4.4 Bounding box over everything - the GEE compositing extent 
aoi.bbox <- aoi.union %>%
  st_bbox() %>%
  st_as_sfc() %>%
  st_sf(geometry = ., aoi = "kidepo.bbox")

## 4.5 Extent arithmetic 
area.pa.km2   <- as.numeric(st_area(aoi.union.utm)) / 1e6
area.buf.km2  <- as.numeric(st_area(st_transform(aoi.union.buf, crs.utm))) / 1e6
area.bbox.km2 <- as.numeric(st_area(st_transform(aoi.bbox, crs.utm))) / 1e6

bb <- st_bbox(aoi.union)

## Pixel counts at 30 m, to sanity-check export size before touching GEE.
px.pa   <- area.pa.km2   * 1e6 / (30 * 30)
px.bbox <- area.bbox.km2 * 1e6 / (30 * 30)

cat("\n[4] AOI GEOMETRY\n")
say(glue("
    PAs included            : {nrow(pa.areas)} of {nrow(pa.register)} in the register
    Sum of individual PAs   : {format(round(area.sum.individual, 1), big.mark = ',')} km2
    Overlap between PAs     : {format(round(area.sum.individual - area.pa.km2, 1), big.mark = ',')} km2 (counted once in the union below)
    Dissolved PA area       : {format(round(area.pa.km2, 1), big.mark = ',')} km2
    Buffered ({buffer.m} m) area  : {format(round(area.buf.km2, 1), big.mark = ',')} km2
    Bounding-box area       : {format(round(area.bbox.km2, 1), big.mark = ',')} km2
    PA share of bbox        : {round(100 * area.pa.km2 / area.bbox.km2, 1)}%
    Bbox (EPSG:4326)        : lon {round(bb['xmin'], 4)} to {round(bb['xmax'], 4)}, lat {round(bb['ymin'], 4)} to {round(bb['ymax'], 4)}
    Approx 30 m pixels      : {format(round(px.pa), big.mark = ',')} inside PAs | {format(round(px.bbox), big.mark = ',')} across the bbox

"))


#  PHASE 5 - EXPORT FOR THE GEE STAGE AND FOR LATER R SCRIPTS

#  GeoJSON in EPSG:4326 is the interchange format: Earth Engine ingests it directly,
#  and st_read() reads it back without the 10-character field-name truncation that
#  shapefiles impose. A zipped shapefile is also written because the Code Editor's
#  asset uploader accepts that form too.

write.geo <- function(x, filename) {
  f <- file.path(path.out.aoi, filename)
  st_write(x, f, delete_dsn = TRUE, quiet = TRUE)
  say(glue("    wrote  {f}  ({format(file.size(f) / 1024, digits = 3)} KB)\n\n"))
  invisible(f)
}

cat("[5] EXPORTS\n")

## 5.1 GeoJSONs 
write.geo(aoi.pa,        "kidepo.pa.individual.geojson")
write.geo(aoi.union,     "kidepo.pa.union.geojson")
write.geo(aoi.union.buf, glue("kidepo.pa.union.buffer{buffer.m}m.geojson"))
write.geo(aoi.bbox,      "kidepo.bbox.geojson")

## 5.1b Zipped shapefile of individual PAs, for the GEE asset uploader
path.pa.shp.dir <- file.path(path.out.aoi, "shp_pa_individual")
dir.create(path.pa.shp.dir, showWarnings = FALSE)
path.pa.shp <- file.path(path.pa.shp.dir, "kidepo_pa_parts.shp")
st_write(aoi.pa, path.pa.shp, delete_dsn = file.exists(path.pa.shp), quiet = TRUE)

path.pa.zip <- file.path(path.out.aoi, "kidepo_pa_parts_shp.zip")

old.wd <- getwd()
zip.ok.pa <- tryCatch({
  setwd(path.pa.shp.dir)
  utils::zip(zipfile = file.path("..", basename(path.pa.zip)),
             files = list.files(".", pattern = "^kidepo_pa_parts\\."),
             flags = "-q")
  TRUE
}, error   = function(e) FALSE,
warning = function(w) FALSE,
finally = setwd(old.wd))

if (isTRUE(zip.ok.pa) && file.exists(path.pa.zip)) {
  say(glue("    wrote  {path.pa.zip}  (ready for the GEE asset uploader)\n\n"))
} else {
  say(glue("    NOTE: could not create {path.pa.zip} automatically - zip the files in\n",
           "          {path.pa.shp.dir} manually.\n\n"))
}

## 5.2 Zipped shapefile of the dissolved union, for the GEE asset uploader -
path.shp.dir <- file.path(path.out.aoi, "shp")
dir.create(path.shp.dir, showWarnings = FALSE)
path.shp <- file.path(path.shp.dir, "kidepo_pa_union.shp")
st_write(aoi.union, path.shp, delete_dsn = file.exists(path.shp), quiet = TRUE)

path.zip <- file.path(path.out.aoi, "kidepo_pa_union_shp.zip")

## zip() needs a zip binary on PATH; Rtools supplies one. Zip from inside the folder
## so the archive has no directory prefix (the GEE uploader wants the files at the
## archive root). `finally` restores the working directory on every exit path -
## on.exit() would not fire here, because tryCatch's body is not a function.
old.wd <- getwd()
zip.ok <- tryCatch({
  setwd(path.shp.dir)
  utils::zip(zipfile = file.path("..", basename(path.zip)),
             files = list.files(".", pattern = "^kidepo_pa_union\\."),
             flags = "-q")
  TRUE
}, error   = function(e) FALSE,
   warning = function(w) FALSE,
   finally = setwd(old.wd))

if (isTRUE(zip.ok) && file.exists(path.zip)) {
  say(glue("    wrote  {path.zip}  (ready for the GEE asset uploader)\n\n"))
} else {
  say(glue("    NOTE: could not create {path.zip} automatically - no zip binary found.\n",
           "          Either upload the GeoJSON instead, or zip the files in\n",
           "          {path.shp.dir} by hand.\n\n"))
}

## 5.3 The audit table, as a CSV for the write-up 
path.csv <- file.path(path.out.aoi, "kidepo_pa_audit.csv")
write_csv(pa.audit %>% select(-matches("^geometry$")), path.csv)
say(glue("    wrote  {path.csv}\n\n"))

## 5.4 A machine-readable handle on the AOI extent, for the GEE script 
## The GEE script needs the bbox numerically. Writing it out means the two stages
## can never silently disagree about the extent.
aoi.meta <- list(
  generated       = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
  crs.geographic  = glue("EPSG:{crs.geo}"),
  crs.projected   = glue("EPSG:{crs.utm}"),
  buffer.m        = buffer.m,
  n.pa.on.disk    = nrow(pa.areas),
  n.pa.register   = nrow(pa.register),
  pa.absent       = pa.absent$pa.label,
  area.pa.km2     = round(area.pa.km2, 2),
  area.buf.km2    = round(area.buf.km2, 2),
  area.bbox.km2   = round(area.bbox.km2, 2),
  bbox.4326       = list(xmin = as.numeric(bb["xmin"]), ymin = as.numeric(bb["ymin"]),
                         xmax = as.numeric(bb["xmax"]), ymax = as.numeric(bb["ymax"]))
)

path.json <- file.path(path.out.aoi, "kidepo_aoi_meta.json")
writeLines(jsonlite::toJSON(aoi.meta, auto_unbox = TRUE, pretty = TRUE), path.json)
say(glue("    wrote  {path.json}\n\n"))


#  PHASE 6 - QC MAP

#  A deliberately plain overview: are the polygons where we expect, and is the
#  bounding box as oversized as the scattered geometry implies? Publication-quality
#  cartography comes later, with the classified rasters.

pa.labels.sf <- aoi.pa %>%
  group_by(pa.label, country) %>%
  summarise(.groups = "drop") %>%
  st_make_valid()

## Label anchors. point_on_surface() (not centroid) guarantees the point falls
## INSIDE the polygon, which matters for the crescent-shaped reserves here, and it
## is computed in UTM so sf does not warn about planar operations on lon/lat.
pa.labels.sf$centroid <- pa.labels.sf %>%
  st_transform(crs.utm) %>%
  st_geometry() %>%
  st_point_on_surface() %>%
  st_transform(crs.geo)

gg.aoi <- ggplot() +
  geom_sf(data = aoi.bbox, fill = NA, colour = "grey40", linetype = "22", linewidth = 0.4) +
  geom_sf(data = aoi.union.buf, fill = NA, colour = "#E07B39", linewidth = 0.3) +
  geom_sf(data = aoi.pa, aes(fill = country), colour = "grey20", linewidth = 0.2, alpha = 0.85) +
  # fun.geometry = identity: the anchors are ALREADY points (computed in UTM above),
  # so ggplot must not re-run st_point_on_surface() on lon/lat and warn about it.
  geom_sf_text(data = pa.labels.sf, aes(geometry = centroid, label = pa.label),
               fun.geometry = identity,
               size = 2.4, colour = "grey10", check_overlap = TRUE) +
  scale_fill_manual(values = c("Uganda" = "#4C9F70", "South Sudan" = "#6A8EAE"),
                    na.value = "grey70", name = NULL) +
  labs(
    title    = "Kidepo landscape AOI - WDPA protected areas",
    subtitle = glue("{nrow(pa.areas)} PAs on disk | dissolved area {round(area.pa.km2)} km2 | ",
                    "bounding box {format(round(area.bbox.km2), big.mark = ',')} km2"),
    caption  = glue("Dashed = GEE compositing bbox; orange = {buffer.m} m buffered union. ",
                    "Areas computed in EPSG:{crs.utm}.")
  ) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom",
        panel.grid = element_line(linewidth = 0.2, colour = "grey92"))

path.png <- file.path(path.out.aoi, "kidepo_aoi_overview.png")
ggsave(path.png, gg.aoi, width = 7.5, height = 8, dpi = 200)
say(glue("[6] wrote  {path.png}\n\n"))


#  PHASE 7 - NEXT STEP

cat("================================================================\n")
cat(" DONE. Next:\n")
cat("  1. Upload outputs/aoi/kidepo.bbox.geojson (or the zipped shapefile) as a GEE\n")
cat("     asset under users/imnduta1/Ismurn, e.g. .../kidepo_bbox\n")
cat("  2. Also upload kidepo.pa.union.buffer500m.geojson as .../kidepo_pa_union\n")
cat("  3. Run scripts/gee/01_hls_coverage_check.js to count clean scenes per\n")
cat("     candidate year, then choose the two comparison dates.\n")
cat("================================================================\n\n")

sessionInfo()
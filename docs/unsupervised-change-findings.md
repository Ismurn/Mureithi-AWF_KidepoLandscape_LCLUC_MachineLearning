# Kidepo landscape LCLUC — unsupervised stage findings

**Status:** unsupervised stage complete; supervised classification not started.
**Dates compared:** 1 Feb – 31 Mar, 2019 and 2024.
**Analysed area:** 7,771.7 km² (97.97% of the dissolved protected-area union).
**Author:** Isaac Mureithi — IMureithi@awf.org

---

## 1. Headline

Across five years, **the protected-area network got greener, and there is no
evidence of widespread land-cover loss.** For an audience expecting a
degradation narrative, that is the finding.

The change signal is dominated by a landscape-wide greenness shift between the
two dry seasons, not by conversion and not by fire. Once a change is required to
be large in spectral terms, not explained by fire, and big enough to be a real
patch, the candidate area falls to a small fraction of the landscape.

Three numbers that should **not** be quoted as change, and why:

| Figure | What it actually is |
|---|---|
| 44.6% of area "changed" | k-means cluster membership flipping, mostly at arbitrary cluster boundaries |
| 1,600 km² "substantive change" | pixels moving 2+ cluster ranks — a rank, not a distance |
| 877 km² "residual" | the above, minus fire; inherits the same flaw |

The defensible figure is **277.5 km² (3.57%)** of candidate change — large
spectral movement, no fire to explain it, in patches of 1 ha or more. Candidate,
not confirmed.

---

## 2. Processing chain

| Stage | Script | Output |
|---|---|---|
| AOI from WDPA polygons | `R/01_build_aoi.R` | 13 PAs, 7,910.5 km² dissolved |
| HLS coverage screening | `gee/01_hls_coverage_check.js` | clear-observation counts per year and window |
| Burn-season profiling | `gee/01b_burn_season_profile.js` | monthly and per-window burnt fractions |
| MODIS burn rasters | `gee/01c_export_burn_masks.js` | burnt extent and **scar age** on the composite grid |
| Composites | `gee/02_build_composites.js` | 2 × 10-band Int16 GeoTIFF, EPSG:32636, 30 m |
| Composite QC | `R/02_qc_composites.R` | pass/fail verdict, support maps |
| Clustering | `R/03_kmeans_firstlook.R` | k-means classes, change matrix |
| Cluster diagnostics | `R/03b_kmeans_diagnostics.R` | convergence, silhouette, separability |
| k sensitivity | `R/03c_k_sensitivity.R` | k=3 vs k=6, cluster-free estimate |
| Fire attribution | `R/04`, `R/05` | proxy vs MODIS record |
| Change magnitude | `R/06_change_magnitude.R` | change in NDVI rather than cluster ranks |
| Maps | `R/07_maps.R` | eight cartographic figures |
| Sampling-bias test | `R/08_sampling_bias.R` | how much of the greening is observation depth |
| Candidate examination | `R/09_check_candidates.R` | by-PA areas, patch sizes, boundary vs roughness, burn test |
| Terrain and observation timing | `gee/01d_export_terrain.js` | elevation, slope, aspect, mean observation day per year |
| Terrain-illumination test | `R/10_terrain_test.R` | whether the candidates are a light artefact |
| Solar-geometry self-check | `R/tests/test_solar_geometry.R` | validates the illumination arithmetic in `R/10` |

Composites are HLS (Landsat + Sentinel-2 harmonised), Fmask-masked, median
over the window, with NDVI/NDWI/NBR and a clear-observation count carried as
bands. Reflectance and indices are stored ×10000 as Int16 with −9999 nodata.

---

## 3. Why these dates

Both assumptions the project began with turned out to be wrong, and measuring
them changed the choice of window — though not in the direction expected.

**`dry_early` (Dec–Jan) is the clearer window, not the cloudier one.** It
reaches ≥4 clear observations over ~100% of the AOI in 9 of 10 years, against
95–100% in only 6 of 11 for `dry_late`.

**`dry_late` (Feb–Mar) is not peak burning. December is**, then January and
November. Only 1–7% of the AOI burns *during* `dry_late`, against 12–37% during
`dry_early`.

`dry_late` was chosen anyway, for a reason neither assumption anticipated: a
median composite assumes a **stationary surface**. In `dry_early` a pixel's
observations straddle its own fire; in `dry_late` the scars are already down and
ageing uniformly. A stable surface with fewer looks beats a shifting surface
with more.

**2019 and 2024** are the two least-burnt `dry_late` years on record here (19.8%
and 22.9% of the AOI carrying scars), five years apart, both above 96% clear
coverage in each country, with comparable observation density (8.1 and 10.4 mean
clear looks per pixel).

---

## 4. Composite quality

| | 2019 | 2024 |
|---|---|---|
| Mean clear observations | 8.11 | 10.38 |
| Area with <4 observations | 1.92% | 0.40% |
| Median RED / NIR / SWIR1 | 0.095 / 0.198 / 0.240 | 0.087 / 0.191 / 0.230 |
| Median NDVI | 0.324 | 0.333 |

2.03% of the PA union is thinly sampled in one date or the other, and is
excluded from all change statistics. Grids are identical (4559 × 7009, 30 m,
EPSG:32636), verified with `compareGeom` rather than assumed.

**This was tested, not just flagged** — see §7.1. About 28% of the greening is
attributable to 2024's better sampling; the remaining 72% is not.

---

## 5. What the clustering found, and what it cannot support

k-means was fitted once on a pooled sample from both dates and applied to each,
so clusters mean the same thing in both — clustering the dates independently
would have made the change matrix count coincidences.

The clusters are spatially coherent and interpretable: Imatong resolves as
woodland, the Karamoja reserves grade from burnt/bare through dry grassland.

But **the landscape supports about three robust strata, not six**:

| k | Converged | % variance | Mean silhouette | % pixels ambiguous |
|---|---|---|---|---|
| 3 | yes | 71.1 | **0.523** | **16.1** |
| 4 | yes | 78.7 | 0.490 | 18.9 |
| 5 | no | 81.7 | 0.442 | 24.1 |
| 6 | yes | 84.1 | 0.440 | 24.8 |
| 7–10 | no | 86+ | 0.43 falling | 26+ |

Silhouette falls monotonically as k rises while variance explained climbs —
which is why an elbow plot looked reassuring and the structure underneath was
not. At k=6 a quarter of pixels sit ambiguously between clusters; the closest
centroid pair is 1.88 apart in standardised units against 9.04 for the farthest.
Clusters 1–4 span NDVI 0.253–0.390 and are slices of one continuum.

Those boundary pixels flip between dates for no physical reason, which is a
large part of why the raw change figure came out as high as it did.

---

## 6. Fire does not explain the change

Attribution was done against the **MODIS record**, not a proxy. An earlier
NBR-quantile proxy scored kappa 0.086 (2019) and 0.262 (2024) against that
record — at best weak, at worst chance.

| | |
|---|---|
| MODIS burn switched between dates | 2,289 km² |
| Overlap with substantive change | 524 km² |
| Expected under independence | 471 km² |
| **Enrichment** | **1.11×** |
| Risk ratio | 1.17× |

Burning barely predicts where the clusters changed. Burnt ground covers ~20% of
the AOI in both dates, but it is largely the *same* ground in both, and fire is
therefore a constant rather than a driver for those pixels.

### Why burn scars are nearly invisible here

| Days since burn | Median NBR |
|---|---|
| 0–20 | **+0.294** |
| 40–60 | +0.088 |
| 80–100 | −0.020 |
| 120–160 | −0.066 |
| *(unburnt ground)* | *−0.001* |

The **youngest** scars look greenest. This is not scar recovery — it is
compositing. A pixel that burnt in mid-March has most of its February–March
observations *predating its own fire*, so the median shows unburnt ground. A
pixel that burnt in November has every observation after the fire. A Feb–Mar
median therefore hides recent fires and shows only old, faded ones.

This is a methodological finding worth carrying into any future dry-season
composite work in fire-driven systems: **a seasonal median composite is a poor
instrument for detecting fire that occurred within or just before its window.**

---

## 7. The change, measured properly

Cluster "steps" are a rank, not a distance. The k=6 centroids sit at NDVI 0.252,
0.259, 0.307, 0.387, 0.521, 0.792 — so a 1→3 move counts as two steps while
spanning 0.055 NDVI, and 5→6 counts as one while spanning 0.271.

| Cluster steps moved | % of area | Median \|ΔNDVI\| |
|---|---|---|
| 0 | 55.6 | 0.033 |
| 1 | 23.9 | 0.045 |
| 2 | 14.2 | 0.060 |
| 3 | 6.0 | 0.072 |
| 4+ | 0.37 | **0.220** |

Only the 0.37% moving four or more steps shift substantially.

### Sensitivity to k

| | k=3 | k=6 | ratio |
|---|---|---|---|
| Changed cluster | 1,568 km² (20.2%) | 3,467 km² (44.6%) | 2.21× |
| …not explained by fire | 1,076 km² | 2,289 km² | 2.13× |
| …and \|ΔNDVI\| > 0.15, ≥1 ha | 146.8 km² | 226.2 km² | **1.54×** |

Filtering removes most of the k-dependence but not all of it. A 50% swing from a
parameter the data cannot justify is not a measurement, which is why the
headline estimate below uses **no clustering at all**.

### The cluster-free estimate

Criteria, all defined in sensor units and reproducible without tuning:
|ΔNDVI| > 0.15, no MODIS burn-state switch to explain it, and contiguous
patches ≥ 1 ha.

> ### 277.5 km² — **3.57%** of the analysed area
>
> in contiguous patches of 1 ha or more.

This replaces both cluster-dependent versions (146.8 km² at k=3, 226.2 km² at
k=6). It is **larger** than either, because it drops the requirement that a
pixel also change cluster — the criterion is now purely spectral, so it catches
large spectral moves that happened to stay inside one cluster. That is the
intended behaviour: the clustering was never the thing doing the measuring.

A note against over-reading it: the layer is size-filtered by construction, so
its spatial coherence ratio (20.9, with 1.3% isolated pixels) is **not**
independent evidence of anything. A layer built from ≥1 ha patches is coherent
because it was made that way.

### Direction of change

| | Median ΔNDVI | % greening |
|---|---|---|
| Background | 0.012 | 59.2 |
| Residual (fire removed) | 0.037 | 65.8 |

Landscape-wide median ΔNDVI is **+0.0137** — 2024 greener than 2019.

Land-cover conversion produces clearing as well as regrowth, so a conversion
signal would sit near 50% greening. A landscape moving one way together is a
wetter or later-senescing season, and no burn mask can remove it.

### 7.1 Is the greening just unequal sampling?

2019 carries 8.2 mean clear looks against 2024's 10.5. If cloud screening
removed proportionally more of the greener, cloudier days from 2019, that alone
would manufacture the result. `R/08_sampling_bias.R` tests it three ways.

The decisive test regresses ΔNDVI on the per-pixel **depth gap** with
protected-area fixed effects. Differencing the same pixel across two dates
cancels everything time-invariant about it — terrain, soil, land cover — leaving
the depth gap as what varies:

| | ΔNDVI |
|---|---|
| Mean observed | +0.0176 |
| Attributable to the depth gap | **+0.0049 (28%)** |
| Remaining at zero depth gap | **+0.0127 (72%)** |

**The greening is real but was overstated.** Quote roughly **+0.013**, not
+0.018, and to one significant figure — estimates range from +0.009 to +0.019
across balanced subsets.

Two methodological notes, because the tests disagreed and the disagreement is
instructive:

- A *within-date* regression of NDVI on observation count comes out **negative**
  in both years, which would suggest the mechanism cannot operate at all. That
  comparison is between different pixels, and the ones with fewer looks are the
  cloudy highlands, which are also the green ones. It measures geography as much
  as sampling, and is the weaker estimate.
- Restricting to "balanced" subsets does **not** hold the landscape fixed: mean
  looks climb from 8.2/10.5 to 14.3/14.5 as the rule tightens, so the strictest
  subset is the 7.4% of the landscape with the clearest skies. Those subsets
  confound depth balancing with spatial selection, which is why their medians
  bounce rather than trend. They establish only that every subset still greens.

---

### 7.2 The candidate patches, examined

`R/09_check_candidates.R`. Four findings, one of which overturned my own
prediction.

**They concentrate in the small hill reserves, not the flagship parks:**

| Protected area | Area km² | Candidate km² | % of PA |
|---|---|---|---|
| Timu | 122.1 | 32.3 | **26.5** |
| Rom | 108.9 | 25.6 | **23.5** |
| Lwala | 58.7 | 6.1 | 10.3 |
| Imatong | 1001.2 | 72.1 | 7.2 |
| Kidepo (SS) | 2828.0 | 59.2 | 2.1 |
| Kidepo Valley | 889.6 | 4.2 | **0.5** |

**They sit nearer protected-area boundaries** — median 1,700 m against 3,026 m
for background. In a logistic model with standardised predictors, distance to
boundary (−0.85) outweighs local spectral roughness (+0.49). That matters
because a legal boundary has no inherent spectral gradient, whereas
misregistration produces artefacts at *image* edges. So boundary proximity is
the stronger pattern — though roughness is far from negligible, and boundaries
often follow rivers and ridges, which are image edges too.

**They green, strongly and almost exclusively** — median ΔNDVI **+0.18** with
**82%** greening, against +0.012 and 59% for background. Clearing would brown;
conversion would do both.

**They are not recovering fires.** This was my leading hypothesis — MODIS at
500 m cannot see a few-hectare fire, so small 2019 burns recovering by 2024
would look exactly like this. It is wrong: candidates are only **0.73×** as
likely as background to have been in the 2019 burnt tail, and their 2019 NBR
(+0.081) was already *above* the landscape median (−0.009). They were greener
than average before and greener still after.

**What remains, none of it established:** genuine woody thickening or regrowth;
differential terrain illumination, since candidates favour rough, hilly ground
in reserves that are mostly forested massifs; or real land-use change. The
illumination explanation is the one that can be settled from orbit, and §7.3
sets out the test. Separating the other two needs high-resolution imagery or a
multi-year series. One caution on the patch statistics: eight-connected
labelling threads scattered clumps into one component, so the "largest patch,
3,064 ha" is a sprawl across a hillside, not a block. The median patch is
2.34 ha.

---

### 7.3 The terrain-illumination test — designed, not yet run

**Status: awaiting the GEE export.** `gee/01d_export_terrain.js` and
`R/10_terrain_test.R` are written and syntax-checked; the test cannot run until
`kidepo_terrain_dry_late_utm36n_30m.tif` is exported and downloaded. What
follows is the design and the expectations set in advance of seeing the data —
recorded here so the eventual result cannot be read as a story told around
whatever came out.

**The hypothesis.** Reflectance on a slope depends on how the sun strikes it.
The candidates favour rough, hilly ground, so some of the "change" could be
light rather than land.

**Why the hypothesis is already in trouble, on mechanism.** Both composites use
the *same calendar window*, 1 Feb – 31 Mar. The sun is therefore in nearly the
same place in both, so the local illumination geometry cos(i) is nearly the same
in both, so it very largely **cancels in the difference**. Illumination sets the
*level* of reflectance on a slope; it cannot by itself produce a *change*
between two same-season composites. For it to do so, something about the
observations must differ between the years, which leaves two routes:

- **Route 1 — timing.** If one year's clear looks fall earlier or later in the
  window, the sun sat at a different declination and cos(i) shifts, most on
  steep slopes. The export's `obs_doy_<year>` bands measure this per pixel.
- **Route 2 — shadow masking.** Fmask cannot tell terrain shadow from cloud
  shadow, so shaded slopes lose observations. §7.1 already attributed 28% of the
  greening to observation depth; if that artefact concentrates on steep ground,
  terrain is its mechanism. This is the version most likely to survive — and it
  is not really about illumination at all, but about where the §7.1 confound
  lives.

**What the test can see.** Validated by simulation before the export, in
`R/tests/test_solar_geometry.R`: the differenced estimator recovered a planted
coefficient of 0.30 as **0.3002** over 40 replications with no detectable bias
(t = +0.09), and its reported standard error matched the empirical spread. At
n = 400,000 the standard error on the coefficient is about 0.007, so the test
has ample power.

**What it cannot see is a large effect, because there isn't one available.** The
greening illumination can explain is the coefficient times the mean shift in
cos(i), and for a timing shift of a few days that shift is of order 0.004:

| Coefficient | Greening explained | Share of +0.0127 |
|---|---|---|
| 0.1 | +0.0004 | 3% |
| 0.3 | +0.0012 | 9% |
| 0.5 | +0.0020 | 16% |
| 1.0 | +0.0040 | 32% |

NDVI's cross-sectional spread across aspects in rugged terrain is rarely more
than ~0.1 per unit cos(i), so coefficients near 1.0 are not physically
available. **The ceiling is single-digit percent unless the two years were
observed many days apart** — which is exactly what Phase 2 of the script
measures, and why the obs_doy bands matter more here than the DEM does.

**The diagnostic to watch is aspect, not slope.** "Candidates are steeper" is
ambiguous: woody vegetation genuinely prefers hills in Karamoja. But there is no
ecological reason for regrowth to prefer west-facing over east-facing ground at
the same slope and elevation, whereas illumination has exactly that signature —
at a mid-morning overpass the sun sits east-southeast (computed: zenith 26–34°,
azimuth 89–130° across the window), so east-facing slopes are lit and
west-facing ones are not. A north–south asymmetry would *not* be illumination:
at 3–4°N the sun passes nearly overhead and both are lit about equally.

**Reading the result.** If cos(i) barely moved and there is no aspect
fingerprint, the hypothesis is closed on mechanism and two explanations remain
that this two-date design cannot separate. If illumination does prove material,
the composites need a topographic correction (C-correction or SCS+C) *before*
the supervised stage — otherwise the Random Forest would learn illumination as
if it were land cover.

---

## 8. Maps and figures

All generated by `R/07_maps.R` into `outputs/maps/` (gitignored — re-run the
script to regenerate). Built with ggplot2 + tidyterra, ggspatial for scale bars
and north arrows, patchwork for panels, rnaturalearth for country context.

| File | Shows |
|---|---|
| `map_01_study_area.png` | The 13 protected areas by country, with context and labels |
| `map_02_composites.png` | False-colour SWIR2/NIR/GREEN composites, 2019 vs 2024 |
| `map_03_observation_support.png` | Clear looks per pixel — the honesty map |
| `map_04_kmeans_classes.png` | k-means classes both dates, semantic palette |
| `map_05_fire_state.png` | MODIS burnt area by date |
| `map_06_scar_age.png` | Days since burn at the window's end |
| `map_07_dndvi.png` | Greenness difference, diverging scale |
| `map_08_candidate_change.png` | The 277.5 km² of candidate change |

Zoom panels of the largest candidate patches, before and after, are written by
`R/09_check_candidates.R` into `outputs/candidates/`.

Diagnostic plots from earlier scripts remain in `outputs/qc/` and
`outputs/classification/`. Those are working figures, not presentation ones —
no scale bars, no north arrows — and should not go in a report.

---

## 9. What this does and does not establish

**Established:**

- The two composites are sound, matched and comparable, verified rather than assumed.
- Fire does not drive the change signal, tested against recorded fires.
- The raw unsupervised change figures are inflated by cluster-boundary effects.
- The dominant signal is a modest, one-directional greening.

**Not established:**

- That the remaining candidate patches *are* land-cover change. They are
  candidates: a large greenness change with no *recorded* fire to explain it.
  §7.2 rules out recovering small fires and shows the patches green from an
  already-vegetated baseline, concentrated in the small hill reserves. What
  that greening *is* remains open.
- Any thematic class labels. k-means finds spectral clusters; the tentative
  names in the scripts are read off centroid signatures and are unverified.
- Anything about causes. Greening is consistent with rainfall differences, woody
  encroachment, grazing changes or recovery, and this analysis separates none
  of them.

**Carry these caveats:**

1. 2019 has ~25% fewer clear observations than 2024. Tested (§7.1): it inflates
   the greening by about 28%, so the corrected figure is ~+0.013 ΔNDVI. The
   candidate-change areas in §7 have **not** been recomputed on a depth-adjusted
   basis and are therefore likely a little generous.
2. MODIS burned area is 500 m information — good for landscape attribution, not
   for per-patch verdicts.
3. Two dates cannot separate inter-annual variability from trend. A multi-year
   time series would.
4. Protected-area polygons overlap by 584 km², so per-PA areas double-count
   unless overlaps are assigned. The §7.2 per-PA table inherits this.
5. Terrain illumination has not been excluded as a driver of the candidate
   patches, which favour rough ground in hill reserves. Testing it needs a DEM.

---

## 10. Data quality problems found and fixed

Recorded because each was invisible until specifically tested for.

| Problem | How it surfaced | Fix |
|---|---|---|
| All six reflectance bands exported as zero | QC printed values; a range check had passed because zero is "plausible" reflectance | HLS v002 in GEE is already reflectance-scaled; the 0.0001 factor was applied twice. Band-mean smoke test added before export, and a degeneracy check (`min == max`) in QC |
| Band means 5× too low | Median over valid pixels disagreed with the mean | GEE writes pixels outside the clip as 0, not the sentinel; ~70% of the raster was zero-fill. All statistics now computed AOI-masked; GEE clips before unmasking |
| MODIS binary burn bands exported as zero | Age band was correct from the same collection and dates | Likely an `ee.Algorithms.If` branch-projection mismatch. Binary bands now derived from the age band |
| k never actually chosen | The script's own comment said to revisit it after the elbow; it was never acted on | `03b` supplies convergence, silhouette and separability; k=6 retained as reconnaissance only |

---

## 11. Next

1. **Run `gee/01d_export_terrain.js`** and download both outputs, then
   `R/10_terrain_test.R`. This is the blocking step: it settles §7.3, and it also
   decides whether the composites need a topographic correction *before* the
   supervised stage — a Random Forest trained on uncorrected slopes would learn
   illumination as if it were land cover. `R/tests/test_solar_geometry.R` already
   passes, so the arithmetic is not what is waiting.
2. **Agree a class scheme** before labelling. The k=3 strata are a sensible
   skeleton; burnt vs bare and grassland vs wooded grassland need training data.
3. **Training data** for the supervised stage — the open decision.
4. **Check candidate patches** against high-resolution imagery before any is
   described as land-use change. Neither §7.2 nor §7.3 can do this; it needs eyes
   on imagery this pipeline does not carry.
5. **Consider a time series.** The central limitation is two dates. An annual
   NDVI series across the same window would separate trend from the inter-annual
   variability that dominates here, and is cheap in GEE. No amount of care with
   the two dates we have will substitute for it.

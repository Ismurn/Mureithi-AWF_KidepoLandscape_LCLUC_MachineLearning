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

**Caveat:** 2019 has ~25% fewer clear looks than 2024. Some of the greening gap
below could be composite depth rather than ground.

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

---

## 8. What this does and does not establish

**Established:**

- The two composites are sound, matched and comparable, verified rather than assumed.
- Fire does not drive the change signal, tested against recorded fires.
- The raw unsupervised change figures are inflated by cluster-boundary effects.
- The dominant signal is a modest, one-directional greening.

**Not established:**

- That the remaining candidate patches *are* land-cover change. They are
  candidates: a large greenness change with no fire to explain it.
- Any thematic class labels. k-means finds spectral clusters; the tentative
  names in the scripts are read off centroid signatures and are unverified.
- Anything about causes. Greening is consistent with rainfall differences, woody
  encroachment, grazing changes or recovery, and this analysis separates none
  of them.

**Carry these caveats:**

1. 2019 has ~25% fewer clear observations than 2024.
2. MODIS burned area is 500 m information — good for landscape attribution, not
   for per-patch verdicts.
3. Two dates cannot separate inter-annual variability from trend. A multi-year
   time series would.
4. Protected-area polygons overlap by 584 km², so per-PA areas double-count
   unless overlaps are assigned.

---

## 9. Data quality problems found and fixed

Recorded because each was invisible until specifically tested for.

| Problem | How it surfaced | Fix |
|---|---|---|
| All six reflectance bands exported as zero | QC printed values; a range check had passed because zero is "plausible" reflectance | HLS v002 in GEE is already reflectance-scaled; the 0.0001 factor was applied twice. Band-mean smoke test added before export, and a degeneracy check (`min == max`) in QC |
| Band means 5× too low | Median over valid pixels disagreed with the mean | GEE writes pixels outside the clip as 0, not the sentinel; ~70% of the raster was zero-fill. All statistics now computed AOI-masked; GEE clips before unmasking |
| MODIS binary burn bands exported as zero | Age band was correct from the same collection and dates | Likely an `ee.Algorithms.If` branch-projection mismatch. Binary bands now derived from the age band |
| k never actually chosen | The script's own comment said to revisit it after the elbow; it was never acted on | `03b` supplies convergence, silhouette and separability; k=6 retained as reconnaissance only |

---

## 10. Next

1. **Agree a class scheme** before labelling. The k=3 strata are a sensible
   skeleton; burnt vs bare and grassland vs wooded grassland need training data.
2. **Training data** for the supervised stage — the open decision.
3. **Check candidate patches** against high-resolution imagery before any is
   described as land-use change.
4. **Consider a time series.** The central limitation is two dates. An annual
   NDVI series across the same window would separate trend from the inter-annual
   variability that dominates here, and is cheap in GEE.

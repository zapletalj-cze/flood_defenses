<!-- Author: Jakub Zapletal -->
<!-- Date: 2026-11-04 -->
<!-- Stage E, step E.2: outputs check, narrow corridor read benchmark, plan and measured runtime estimates (for approval). -->

# Stage E, step E.2: outputs check, read benchmark and plan

Status: **proposal, waiting for approval.** No Stage E code has been written yet. The prototypes behind the timings
below live only in a scratch folder.

---

## 1. Outputs of stages A–D (E.2 item 1)

**Configured country.** `config/defenses.toml` holds only the template country `XX`, and no data of yours is
reachable from this session. Each Stage E script will check its inputs when it starts:
- **Missing required input:** the script stops and lists it.
- **Missing optional input:** the script logs which features are off and continues.

| Input | Kind | Feature switched off if missing |
|---|---|---|
| `units_core.gpkg`, `units_window.gpkg`, per-unit `hand_undef.tif` | required | — |
| `rivers_sources.gpkg` (`is_source`) or the rivers input | required | — |
| building footprints or GHSL built-up surface (per `built_source`) | required | — |
| `zones.gpkg` | optional | context label; score weights fall back to the `rural` parameters |
| `barriers_all.gpkg`, `crest_lines.gpkg` | optional | crest and closure support. **Proposed fallback:** `barrier_lines.gpkg` of stage A for `crest_support`, logged |
| `compartments_rp.gpkg` | optional | precedence (`explicit_overlap_share`); RP integration |
| `gfm_evidence.csv` and GFM scenes | optional | GFM corroboration; the negative metric of the evaluation |
| GloFAS event RPs (`rp_bounds.py`) | optional | event-RP caps; errors are then counted as "potential" |
| `reference_protected_areas` | optional | the positive metric of the evaluation |

**Synthetic test country (T4/T4D).**
- **Required inputs:** all present, except the GHSL built-up surface. No script downloads that product yet; see
  section 3.1. Building footprints are available from OSM.
- **Optional inputs present:** `compartments_rp.gpkg`, `gfm_evidence.csv`, and `zones.gpkg` (T4D).
- **Optional inputs missing:**
  - `barriers_all.gpkg` and `crest_lines.gpkg` (Stage D step D2);
  - GloFAS bounds, because the test restores the Stage B result.

---

## 2. Narrow corridor reads along rivers (E.2 item 2)

**Setup.**
- **DTM:** synthetic benchmark world of Stage A, 24 × 24 km at 1 m.
- **Rivers:** 3 lines, 47 km, cut into 12 chunks of 5 km, each buffered by 30 m (`max_offset_m`). The corridors
  hold 11.4 MB of native pixels.
- **Cold reads:** files evicted from the page cache before each method.

| Layout | `gdal.Warp` average to 2 m over the chunk bounding box | Block-aligned reads, native | Block-aligned reads, decimated to 2 m |
|---|---|---|---|
| strip (12000 × 1 blocks, no overviews) | 5.57 s; 154.7 MB native-equivalent, **28 MB/s** | 1438 MB of full-width rows in 8.71 s, **165 MB/s** (126 × the corridor) | 11.02 s, **131 MB/s** |
| tiled 512², ZSTD, overviews (option A) | 0.71 s, **217 MB/s** | 122 MB of touched tiles in 0.61 s, **200 MB/s** (11 × the corridor) | 0.86 s, **143 MB/s** |

**Consequences.**
- **Cost per kilometre of river.** On a tiled DTM, snapping reads cost about 0.015 s per km of river (0.7 s per
  47 km). On a strip DTM they cost about 0.12 s per km, eight times more.
- **Block-aligned reads** pay off only on tiled files. On strip files they decode whole rows.
- **Proposal.** Use the chunk-bounding-box warp, as the specification says. Recommend conversion option A wherever
  the DTM is still in strip layout.

---

## 3. Plan

### 3.1 Files

| Step | File | Content |
|---|---|---|
| E1 | `processing/snap_rivers.py` (new) | thalweg snapping by chunks; `rivers_snapped.gpkg`, `river_points.gpkg`; plots `snap_profile`, `snap_map` |
| E1 | `plots.py` (new, next to `gis.py`) | shared plotting helpers; see 3.2 |
| E2 | `processing/implied_segments.py` (new) | envelope, nearest river point, bank side, subsegments per unit (checkpointed, parallel); `subsegments.tif` and lookup CSV per unit, merged `subsegments.gpkg` |
| E3 | `processing/implied_protection.py` (new), parts A–B | built-up fraction, curves (`curves.npz`), quantile and onset levels, sharpness, eligibility; `implied_subsegments.gpkg`; plots `curve`, `summary` |
| E3 | `download/download_ghsl.py` | optional GHSL built-up surface (`built_s`), per epoch for `ghsl_epochs` |
| E3 | `download/download_osm.py` | `man_made` in `osm_buildings.gpkg`; counts of `exclude_building_tags` |
| E4 | `processing/implied_protection.py`, parts C–F | implied cells, merging, precedence, temporal signal, corroboration (crest, closure, GFM), RP stratum `implied`, confidence; `implied_compartments.gpkg`, `run_config.toml`; plots `longitudinal`, `temporal`, `compartment_map` |
| E5 | `processing/implied_evaluate.py` (new) | sensitivity grid, positive and negative metrics; `sensitivity.csv`; `sensitivity` plots |
| E5 | integration | `defenses_clip_with_implied.gpkg` (`include_implied`), review queue entries (`include_in_review`), `use_snapped_rivers` switch in `compartments.py` (default off; Stage A results unchanged) |
| E1–E5 | `tests/test_stage_e.py` (new) | the checks of E.7, step by step |

### 3.2 Where the plotting helpers live: `plots.py` next to `gis.py`

**Why not in `gis.py`.** Every script imports `gis.py`. Putting matplotlib there would make it a dependency of the
whole package and would slow every script start. `plots.py` is imported only by the Stage E scripts.

**What it holds:**
- **Backend:** Agg, unless plots are to be shown.
- **Figure budget:**
  - flagged and anomalous items first, then a random sample with `sample_seed`, up to `max_per_type`;
  - all items in `plot_all_in_units`.
- **Saving:** save, then close every figure.
- **Showing:** non-blocking display with `plt.pause`, up to `show_max` figures, skipped without a display.
- **Look:** a consistent colour per level method.
- **Index:** `index.html` (thumbnails per type, links, metric table), written with the standard library.

**Hillshade crop.** `gdal.DEMProcessing` to MEM goes to `gis.Raster`, because it has no matplotlib dependency.

### 3.3 gis.py additions

All additions are extensions. Two moves of existing code are needed, and both are guarded by the existing tests.

| Function | Purpose |
|---|---|
| `Hydro.least_cost_path(cost, start, end)` | numba Dijkstra, 8-connected, diagonal cost × √2; NaN = impassable (reuses the heap of the priority flood) |
| `Vector.chaikin(line, iterations, max_shift)` | Chaikin smoothing that keeps vertices which would move more than `max_shift` |
| `Hydro.hand_above_streams(dtm_undef, stream, burn_m, burn_buffer_cells)` | **move** of the HAND code path of `compartments.py` (breaching, river reference, fallback), so that `recompute_hand` uses the same code; `compartments.py` calls it; stage A results unchanged (pipeline tests) |
| `Stats.binned_quantiles(area_per_bin, edges, quantiles)` | area-weighted quantiles interpolated within bins |
| `Stats.onset(share, smooth_bins, fraction, persist_bins)`, `Stats.sharpness(...)` | onset level and sharpness of a curve; shared by `implied_protection.py` and `implied_evaluate.py` |
| `Evidence.test_regions(...)` | **move** of the comparison logic of `gfm_evidence.py`, so that implied compartments are tested by the same code; Stage C tests guard it |
| `Raster.hillshade(array, geotransform, wkt)` | hillshade with `gdal.DEMProcessing` to MEM |

### 3.4 Configuration
- **Spec tables.** `[implied]`, `[implied.snap]` and `[implied.plots]` as in E.3, with all values marked as
  placeholders.
- **Additions:**
  - `[countries.XX.compartments] use_snapped_rivers = false`;
  - `sources.ghsl_url_templates.built_s` (not verified).
- **`--show-plots` becomes a run variable.** `SHOW_PLOTS = False` (`DEFENSES_SHOW_PLOTS`), because the scripts take
  no command-line parameters, as you asked earlier.
- **pixi.toml:** `matplotlib-base` (approved), plus tasks `snap-rivers`, `implied-segments`, `implied-protection`,
  `implied-evaluate`, `test-stage-e`.
- **Naming in all new code:**
  - `gdf_…` for GeoDataFrames, `df_…` for DataFrames, `list_…` for lists;
  - `gdf_sjoin_<a>_<b>` for spatial joins;
  - short names.

  Code of stages A–D stays as it is unless you want a separate mechanical renaming pass; the tests would guard it.

---

## 4. Runtime estimates, from measurements

**Prototypes, measured:**

| Item | Measurement | Estimate |
|---|---|---|
| Snapping, one 5 km chunk at 2 m (tiled source) | read 0.09 s, cost (minimum filter + corridor) 0.08 s, least-cost path 0.03 s (path of about 2000 cells) | about 40 s per 1000 km of river (tiled); about 2 min on a strip DTM (read 0.46 s per chunk, section 2) |
| Nearest river point (EDT with indices), 30.25 M cells | 3.0 s; indices int32, 242 MB | about 3 s per unit; with `recompute_hand` +13.9 s HAND (Stage A measurement) |
| Curves, 2000 subsegments × 40 bins on 30.25 M cells (`bincount`) | 2.8 s | about 3 s per unit |
| Buildings at a quarter of 5 m (1.25 m) | 100 000 buildings on 10 × 10 km: 2.5 s, block average 0.4 s | dense urban unit (27.5 km side): about 20 s. **Memory:** done in tiles of 4096², because a whole window at 1.25 m is 484 MB of uint8 |
| GHSL warps, totals preserved | 1.6–4.6 s per unit for two products (Stage B) | about 1–2 s per epoch per unit |
| Figures (Agg) | curve 0.16 s; hillshade map 2000² 0.53 s | budget of 200 per type: at most 32 s (curves) or 106 s (maps) per type and run |

**To be measured in the steps.**
- **E4:** merging and polygonizing implied compartments.
- **E5:** the sensitivity grid (4 methods × 3 shares × 3 lengths = 36 runs on the pilot units).
- **Pilot units first** (`implied.units`); a country estimate follows from their measured times.

---

## 5. Points to decide

1. **Dependency on Stage D.** `barriers_all.gpkg` and `crest_lines.gpkg` come from D2; `defenses_clip.gpkg` and
   `review_queue.py` come from D4, and these steps are not done yet.
   - **Proposal:** E1–E3 now; they do not depend on D2–D4.
   - **E4:** uses the D2 outputs when present, otherwise falls back to stage A barrier lines (logged).
   - **E5:** integration after D4.
   - **Alternative:** finish D2–D4 first.
2. **GHSL built-up surface.**
   - **Unverified:** the product and epochs (10 m for one year, 100 m multi-temporal) and the download templates.
     Local files would be used until they are confirmed.
   - **Consequence:** with footprints only, `built_source = "footprints"` works, but the temporal signal (part D)
     needs the multi-temporal surface.
3. **Two code moves to `gis.py`** (HAND path, GFM comparison logic), guarded by the existing tests.
4. **`--show-plots`** becomes the run variable `SHOW_PLOTS`.
5. **WhiteboxTools** is not part of the environment. The numba least-cost path measured 0.03 s per chunk, so it is
   used.
6. **`exclude_building_tags`.** taginfo is not reachable, so the counts come from the country extract (logged in
   E3).
7. **Renaming stages A–D to the new naming convention:** yes or no.

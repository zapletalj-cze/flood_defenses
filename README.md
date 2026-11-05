# flood_defenses — Defenses layer for the fluvial flood hazard pipeline

The package derives **compartments**: areas protected by flood levees. Each compartment
gets a return period (RP) of protection. The layer clips fluvial hazard maps in
post-processing, so it is built **only from terrain, levee lines, the river network and
free open data** (OpenStreetMap, GHSL, Copernicus GFM and GloFAS). No hazard map, hydraulic
model output, local hazard dataset or FABDEM is used anywhere. Satellite flood extents
(GFM) and reanalysed discharge (GloFAS) enter only as *observations* that bound the
modelled return period; they never draw a compartment.

**Czech quick guide:** [docs/navod_cs.md](docs/navod_cs.md): a short guide to running the
scripts.

**Error priority: recall.** When the evidence is ambiguous, a compartment is kept rather
than dropped. Every design choice below that trades precision for recall says so
explicitly.

| Stage | Content | Status |
|---|---|---|
| A | `gis.py` additions, `prepare_dtm.py`, `build_units.py`, `compartments.py` | done, tested, benchmarked |
| A+ | faster kernels; fixes B1–B4; gap filling for water surfaces | done, tested, benchmarked |
| B | `download_osm.py`, `osm_hydraulics.py`, `download_ghsl.py`, `rp_index.py` | done, tested on synthetic data; download URLs not verified |
| C | `download_gfm.py`, `gfm_evidence.py`, `download_glofas.py`, `rp_bounds.py`, `validate.py` | done, tested on synthetic data (26 checks); GFM endpoint, asset names and pixel codes, and the GloFAS day convention not verified |
| HAND_YZ | `download_worldcover.py`, `HAND_YZ.py`: land cover exposure supplementing GHSL | done, tested on synthetic data (7 checks); WorldCover source verified, Urban Atlas field and codes not verified |
| D | industrial sites, dense urban areas, large systems ([plan](docs/stage_d_plan.md)) | **D1 done** (DSM inspection, OSM and GHSL extensions, `download_clms.py`, `context_zones.py`; 10 checks); D2–D5 pending |
| E | implied protection in dense built-up reaches, experimental ([plan](docs/stage_e_plan.md)) | **E1–E3 done** (`snap_rivers.py`, `plots.py`, `implied_segments.py`, `implied_protection.py` A–B; 29 checks); E4–E5 in progress |

---

## Contents

1. [Running the scripts](#1-running-the-scripts)
2. [Repository layout](#2-repository-layout)
3. [Configuration](#3-configuration)
4. [Run order and outputs](#4-run-order-and-outputs)
5. [Method: from levee lines to protected compartments](#5-method-from-levee-lines-to-protected-compartments)
6. [Large DTMs: corridors, storage layout, units](#6-large-dtms-corridors-storage-layout-units)
7. [Parallel processing, checkpoints, merging](#7-parallel-processing-checkpoints-merging)
8. [Water sources and pumping stations from OpenStreetMap](#8-water-sources-and-pumping-stations-from-openstreetmap)
9. [Return period index](#9-return-period-index)
10. [Observational evidence: Copernicus GFM](#10-observational-evidence-copernicus-gfm)
11. [Observed bounds of the return period: GloFAS](#11-observed-bounds-of-the-return-period-glofas)
12. [Validation](#12-validation)
13. [Measured performance](#13-measured-performance)
14. [Tests](#14-tests)
15. [Bugs, fixes, limitations, open decisions](#15-bugs-fixes-limitations-open-decisions)
16. [Data sources and licences](#16-data-sources-and-licences)
17. [Stage D (in progress): industrial sites, dense urban areas, large systems](#17-stage-d-in-progress-industrial-sites-dense-urban-areas-large-systems)
18. [Stage E (experimental): implied protection in dense built-up reaches](#18-stage-e-experimental-implied-protection-in-dense-built-up-reaches)

---

## 1. Running the scripts

**No command-line parameters.** Each script starts with a block of *run variables*:

```python
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"
FORCE = False   # recompute finished units / footprint / conversion / extraction
UNITS = []      # compartments.py, rp_index.py, download_gfm.py, gfm_evidence.py: selected units, empty = all
```

`gfm_evidence.py` has three more: `HINTS` (missing-barrier hints, §10.3) and `DATE_FROM`,
`DATE_TO` (restrict the run to the scenes of one flood event).

Everything else comes from the TOML configuration: data paths, thresholds, resolutions,
URLs and tag lists. Code holds no paths, CRS, URLs or thresholds.

Batch runs and the tests override the run variables through environment variables named
`DEFENSES_<VARIABLE>`, for example `DEFENSES_COUNTRY=CZ` or
`DEFENSES_UNITS=012_034,012_035`. Without them the values in the script apply.

**Scripts run as modules from the repository root.** `python -m package.script` puts the
working directory on `sys.path`, so `import gis` works without `sys.path` manipulation or
packaging:

```
pixi run prepare-dtm        pixi run download-osm       pixi run download-gfm       pixi run test-kernels
pixi run build-units        pixi run osm-hydraulics     pixi run gfm-evidence       pixi run test-pipeline
pixi run compartments       pixi run download-ghsl      pixi run download-glofas    pixi run test-stage-b
                            pixi run rp-index           pixi run rp-bounds          pixi run test-stage-c
pixi run download-worldcover                            pixi run validate           pixi run test-hand-yz
pixi run hand-yz
```

`pixi run python -m processing.compartments` works the same way. Running a script by its
file path does **not** work, because Python then puts the script folder on the path
instead of the root.

**Dependencies** are listed in `pixi.toml`:
- `python ≥ 3.11` (for `tomllib` and `ProcessPoolExecutor(max_tasks_per_child)`),
- `gdal ≥ 3.11` (`MEM` vector datasets; the OSM driver; SpatiaLite SQL for the OSM
  extraction),
- `libgdal-netcdf`, `libgdal-grib`: GDAL's NetCDF and GRIB drivers, which conda-forge
  ships as separate plugin packages (GloFAS files),
- `numpy`, `scipy`, `numba`, `geopandas ≥ 1.0`, `shapely ≥ 2.0`, `pyproj`,
- `pandas ≥ 2.0` (already used by `gis.py`; 2.0 for ISO 8601 date parsing),
- `pyflwdir`, only because the existing `gis.py` imports it (it depends on affine,
  numba, numpy and scipy, not on rasterio).

Downloads and the STAC search use the standard library (`urllib`). No rasterio-based
library is used. `cdsapi` is **not** part of the environment: `download_glofas.py` uses it
only if you install it (§11.1); otherwise it reads files downloaded through the EWDS form.

## 2. Repository layout

```
gis.py                     shared GIS library (existing; new functions appended, nothing changed)
plots.py                   Stage E plotting helpers (matplotlib): budget, saving, colours, index.html
config/defenses.toml       configuration template
config/osmconf.ini         OSM driver configuration exposing the tag keys used
docs/navod_cs.md           short guide to running the scripts, in Czech
docs/stage_d_plan.md       Stage D, step D.1: inventory, inputs, plan and measured runtime impact
docs/stage_e_plan.md       Stage E, step E.2: outputs check, corridor read benchmark, plan (approved)
download/
  download_osm.py          Geofabrik extract -> hydraulic points, drains (single-layer GeoPackages)
  download_ghsl.py         GHSL building volume tiles -> VRT mosaics
  download_gfm.py          GFM scene search (STAC), screening, clipped layers per unit
  download_glofas.py       GloFAS daily discharge and upstream area -> compact cell series
  download_worldcover.py   ESA WorldCover tiles over the corridor -> worldcover.vrt
  download_clms.py         CORINE (legend mapping verified on known pixels) and Urban Atlas from local files
processing/
  prepare_dtm.py           DTM inspection, footprint, benchmark, estimates, optional conversion
  build_units.py           river corridor, processing units, memory estimate
  osm_hydraulics.py        water sources (closed streams), pumping stations
  compartments.py          per-unit compartments in parallel, checkpoints, merging
  rp_index.py              exposure, score, RP per compartment
  HAND_YZ.py               land cover exposure (WorldCover, Urban Atlas) supplementing GHSL
  context_zones.py         industrial / dense urban / urban / rural zones of the corridor (Stage D)
  snap_rivers.py           river lines snapped to the thalweg, river points with chainage (Stage E1)
  implied_segments.py      bank subsegments: envelope, nearest river point, side, chainage bins (Stage E2)
  implied_protection.py    built-up fraction, HAND curves and levels per subsegment (E3), implied compartments (E4)
  gfm_evidence.py          held / flooded tests per compartment and scene; missing-barrier hints
  rp_bounds.py             event return periods, observed bounds, combination with the model
  validate.py              area recall and precision; RP confusion matrix (was validate_defenses.py)
tests/
  test_kernels.py          numba kernels on small grids (16 checks)
  test_pipeline.py         end-to-end synthetic countries for Stage A (40 checks)
  test_stage_b.py          end-to-end synthetic country for Stage B (20 checks)
  test_stage_c.py          end-to-end Stage C on the Stage B country (26 checks)
  test_hand_yz.py          land cover exposure on the Stage B country (7 checks)
  test_stage_d.py          Stage D on the Stage B country, country code T4D (D1: 10 checks)
  test_stage_e.py          Stage E on its own synthetic country T5 (E1: 16, E2: 6, E3: 7 checks)
  make_benchmark_world.py  24 x 24 km synthetic DTM for performance runs
download_nld.py            US National Levee Database download (optional)
```

### New `gis.py` functions

Functions are appended to the existing classes or added as new classes; no existing
function was modified. The only change to existing lines is the import block.

| Class | Function | Purpose |
|---|---|---|
| `Raster` | `inspect`, `build_vrt` | Storage properties of a raster; VRT over several DTMs (stops on a CRS mismatch) |
| `Raster` | `read_window(src, bounds, res, resample, warp_mb, overview, multithread, dst_wkt, src_nodata)` | `gdal.Warp` of a map window to MEM, optionally reprojected; `src_nodata` overrides the file's nodata (GFM layers) |
| `Raster` | `read_native(src, bounds)` | Direct read of native pixels |
| `Raster` | `array_to_mem`, `write_array`, `best_compression` | Array ↔ GDAL; tiled compressed GeoTIFF; ZSTD → DEFLATE fallback |
| `Raster` | `rasterize_geometries`, `polygonize_labels`, `fill_nodata_cells`, `footprint` | Burn geometries; label raster → polygons; `gdal.FillNodata` on selected cells; valid-data footprint |
| `Raster` | `benchmark_reads(info, polygon, …)`, `map_codes(values, mapping)` | Read benchmark of one raster, shared by DTM and DSM; raster values → codes through a lookup (CORINE legend) |
| `Vector` | `check_single_layer`, `read_merged`, `write_gpkg`, `length_near_boundary` | Single-layer inputs; merged reads with `source`; single-layer GeoPackage output; levee length near compartment boundaries |
| `Vector` | `layer_names`, `single_layer_copy(path, layer, out_path)` | Layers of a vector source; one layer of a multi-layer source copied to its own GeoPackage with `gdal.VectorTranslate` |
| `Vector` | `zones_to_grid(zones, order, shape, geotransform, wkt)` | Context of every cell of a unit grid (position in `context.order`, 0 = none) |
| `Hydro` | `breach_to_streams(dem, stream, burn, burn_buffer_cells)` | River cell reached downstream (numba) |
| `Hydro` | `priority_flood(surface, seeds)` | Spill level, seed index and reach order per cell (numba) |
| `Hydro` | `distance_cells(mask, limit)` | Chamfer distance in cells (numba) |
| `Config` | `load`, `get`, `as_list`, `warn_placeholders` | TOML reading and validation |
| `Run` | `env_overrides`, `logger`, `peak_rss_mb`, `configure_gdal`, `run_units`, `format_seconds` | Run variables, run logs, memory measurement (Linux/macOS/Windows), process pool |
| `Http` | `download(url, path, retries, pause_s, headers)` | urllib download with size check, skip and backoff; `file://` supported |
| `Http` | `get_json(url, params, body, retries, pause_s, headers)` | JSON GET with query parameters or POST of a JSON body, with backoff (STAC search and its `next` links) |
| `Stats` | `percentile_rank`, `quantile_classes` | Ranks and quantile mapping for the RP score |
| `Stats` | `gumbel_lmoments(sample)`, `gumbel_rp(q, xi, alpha)` | Gumbel fit by L-moments; return period of a discharge |
| `Hydro` | `least_cost_path(cost, start, end)` | Dijkstra on an 8-connected grid (numba), NaN impassable; thalweg snapping |
| `Hydro` | `hand_above_streams(dem, stream, burn, burn_buffer_cells, valid)` | HAND code path of `compartments.py`, moved here so that Stage E recomputes HAND the same way |
| `Vector` | `chaikin(line, iterations, max_shift)`, `interpolate_xy(line, distances)` | Chaikin smoothing that keeps vertices which would move more than `max_shift`; points along a line in one numpy pass |
| `Raster` | `hillshade(array, geotransform, wkt)` | Hillshade with `gdal.DEMProcessing` to MEM (plots) |
| `Stats` | `binned_quantiles(area, edges, quantiles)`, `onset(share, smooth_bins, fraction, persist_bins)`, `sharpness(smoothed, onset_bin, window_bins)` | Area-weighted quantiles within bins; onset and sharpness of a built-share curve (Stage E) |

### Naming in the code

All scripts and tests name their variables by one scheme:
- **Prefixes:** `gdf_…` for GeoDataFrames, `df_…` for DataFrames, `list_…` for lists.
- **Spatial joins:** `gdf_sjoin_<a>_<b>`, for example `gdf_sjoin_src_rn` (river sources joined to the rivers).
- **Short abbreviations:**
  - `rn` rivers, `lv` levees, `bl` barrier lines, `br` bridges, `bar` barriers with crest heights;
  - `comp` compartments, `fp` footprints, `poly` polygons, `win` and `core` unit windows and cores;
  - `ua` Urban Atlas, `ev` GFM evidence, `ovr` overrides.

In `gis.py` the scheme applies to the local variables of the functions added for this package. The original
functions and all function parameters are unchanged.

## 3. Configuration

`config/defenses.toml` is the template:
- `[general]`: directories, workers, memory limits, GDAL cache and threads.
- `[benchmark]`: benchmark windows.
- `[countries.<code>]`: inputs, resolutions, units, conversion and OSM region, with
  sub-tables `compartments` and `rp`.
- `[sources]`: download templates, tag lists, retries.
- `[evidence]`: Stage C.

All path keys accept a string or a list. Every script validates exactly the keys it uses,
before doing any work:
- **presence and type;**
- **file existence;**
- **allowed values;**
- **consistency between keys:**
  - share sums and list lengths;
  - `distribution_classes`, `cap_rp` and `max_rp` must be values of `classes`;
  - the closing and pumping tags must be part of the hydraulic tags;
  - `unit_halo_m ≥ corridor_buffer_m` (reason in §6.4).

All problems are listed in one message and the script stops.

**Placeholders.** The top-level array `placeholders` lists keys whose values are not
confirmed yet. Every script logs a warning for each placeholder it uses.

Keys beyond the original schema, with the reason for each:

| Key | Reason |
|---|---|
| `general.gdal_threads_per_worker` | Avoid oversubscription: n workers × GDAL threads ≤ cores |
| `benchmark.*` | Window count, size and seed of the benchmark (reproducible positions) |
| `country.dtm_conversion`, `dtm_compress`, `footprint_resolution_m` | DTM conversion (§6.3) and footprint |
| `country.bytes_per_cell` | Memory model for unit sizing; replace with the value logged by `compartments.py` |
| `country.crest_box_size_m` | Size of the boxes read at native resolution around levees |
| `country.merge_overlap_ratio`, `drop_contradicted_truncated`, `contradiction_ratio`, `id_rounding_m` | Merging of units and ids (§7.3) |
| `country.keep_intermediate` | Write DTM, spill, margin, barrier and crest rasters per unit |
| `country.geofabrik_region`, `osm_pbf`, `endpoint_elevation_radius_m` | OSM extract and downstream-end detection (§8) |
| `compartments.touch_tolerance_m`, `self_close_min_length_factor` | Former code constants of the bridging |
| `compartments.gated_crossings`, `fill_nodata_holes` | Former command-line switches |
| `compartments.stream_burn_m`, `stream_burn_buffer_m` | B1 fix (§15.1) |
| `compartments.preserve_river_cells` | B2 fix |
| `compartments.cap_river_holes`, `max_fill_distance_m` | Gap filling for water surfaces (§5.5) |
| `rp.value_total_tolerance` | Tolerance of the value-warp totals check (spec: 1 %) |
| `sources.request_pause_s`, `max_retries` | Polite downloads with backoff |
| `sources.osm_tags_closing`, `osm_tags_pumping` | Which hydraulic structures close a river or pump water |
| `sources.ghsl_url_templates` (table `total`, `non_residential`), `ghsl_tile_schema`, `ghsl_tile_field`, `ghsl_nodata_as_zero` | GHSL tiles; a table instead of the list in the original schema, so that products cannot be swapped |
| `country.validation_extent` | Optional area where the reference is complete; recall and precision outside it are meaningless |
| `sources.gfm_collection`, `gfm_assets` (table `flood`, `exclusion`, `reference_water`, `likelihood`), `gfm_page_limit`, `gfm_token_env` | STAC collection, asset keys per role, page size, name of the environment variable holding an access token (never the token itself) |
| `sources.gfm_screen_resolution_m`, `gfm_download_threads` | Resolution of the screening read; scenes screened at the same time (§10.1) |
| `sources.gfm_flood_values`, `gfm_exclusion_values`, `gfm_reference_water_values`, `gfm_nodata_values` | Pixel codes of the GFM layers; placeholders until checked against the product description |
| `sources.glofas_request` (table), `glofas_time_keys`, `glofas_area_key`, `glofas_first_year`, `glofas_last_year`, `glofas_file_extension` | The request pasted from the EWDS form; the script only adds time and area |
| `sources.glofas_files`, `glofas_discharge_variable`, `glofas_upstream_area_file`, `glofas_upstream_area_url`, `glofas_upstream_area_variable`, `glofas_upstream_area_to_km2` | Local GRIB or NetCDF files; the upstream area map and its unit |
| `sources.glofas_min_days_per_year` | A year enters the annual maxima only if it is nearly complete |
| `sources.glofas_valid_time_ends_day` | Day attribution of daily values stamped at the end of their period (§11.1) |
| `evidence.river_level_percentile`, `hint_margin_m`, `hint_dry_share`, `hint_min_scenes`, `hint_exclusion_m` | Missing-barrier hints (§10.3); the specification fixed 90 % and "most scenes" in words |
| `country.lulc.*` (table) | Land cover exposure of `HAND_YZ.py` (§9.3): switch, WorldCover and Urban Atlas inputs, class-to-category mappings, values and curves per category |
| `sources.worldcover_url_template`, `worldcover_grid_url`, `worldcover_tile_field`, `worldcover_version`, `worldcover_year` | ESA WorldCover tiles (§9.3); verified |
| Stage D tables (`[countries.XX.context]`, `[countries.XX.context.params.<context>]`, `[crest_lines]`, `[barriers]`, `[coarse]`, `[review]`) and inputs (`dsm`, `building_footprints`, `manual_barriers`, `industrial_registry`, `ems_flood_extents`, `water_polygons`) | As specified for Stage D (§17) |
| `country.corine`, `corine_legend`, `corine_legend_index_field`, `corine_legend_code_field`, `corine_check` | CORINE from local files: the legend maps class indices to codes, verified on known pixels |
| `sources.osm_tags_water`, `barriers.osm_wall_candidates` | OSM water polygons; wall tags counted in the extract to choose `osm_tags_walls` (taginfo was not reachable) |
| `sources.ghsl_url_templates.smod`, `.pop` | Optional GHSL settlement model and population (Stage D) |

## 4. Run order and outputs

| Step | Script | Writes (`output_dir/<country>/` unless noted) |
|---|---|---|
| 1 | `prepare_dtm.py` | `dtm_inspection.csv`, `dtm_footprint.gpkg`, `dtm_benchmark.csv`, `dtm_estimates.json`, `work/<c>/dtm/dtm.vrt` |
| 2 | `build_units.py` | `corridor.gpkg`, `units_core.gpkg`, `units_window.gpkg`, `units.csv`, `work/<c>/units_domain.gpkg` |
| 3 (optional) | `prepare_dtm.py` with `convert_dtm = true` | tiled copies in `work/<c>/dtm/`; `dtm.vrt` repointed to them |
| 4 | `download_osm.py` | `osm_hydraulic_points.gpkg`, `osm_drains.gpkg` |
| 5 | `osm_hydraulics.py` | `rivers_sources.gpkg`, `osm_pumping.gpkg` |
| 6 | `compartments.py` (first two units, then all) | per unit in `work/<c>/units/<id>/`: `hand_undef.tif`, `labels.tif`, `compartments.gpkg`, `barrier_crest.csv`, `done.json`, `unit.log`; merged: `compartments.gpkg`, `barriers.gpkg`, `enclosed_streams.gpkg`, `performance.json` |
| 7 | `download_ghsl.py` | `ghsl_total.vrt`, `ghsl_non_residential.vrt`, `ghsl_sources.json` |
| D1 (Stage D) | `prepare_dtm.py` with `dsm`; `download_osm.py`; `download_ghsl.py` with `smod`, `pop` | `dsm_inspection.csv`, `dsm_benchmark.csv`, `dsm_footprint.gpkg`; `osm_buildings.gpkg`, `osm_water.gpkg`, `osm_crossings.gpkg`, `osm_walls.gpkg`, `osm_industrial.gpkg`; `ghsl_smod.vrt`, `ghsl_pop.vrt` |
| D1 (Stage D) | `download_clms.py`, then `context_zones.py` | `corine_codes.tif`, `urban_atlas.gpkg`, `clms_sources.json`; `zones.gpkg` |
| 7a (optional) | `download_worldcover.py` | `worldcover.vrt`, `worldcover_sources.json`; tiles in `work/<c>/worldcover/` |
| 7b (optional) | `HAND_YZ.py` | per unit `hand_yz.csv`, `hand_yz.json`; merged `exposure_lulc.csv` |
| 8 | `rp_index.py` | per unit `rp_index.csv`, `rp_index.json`; merged `compartments_rp.gpkg` (with `value_ghsl`, `value_lulc`) |
| 9 | `download_gfm.py` | `work/<c>/gfm/screening.csv` (checkpoint), clipped layers in `work/<c>/gfm/<unit>/`, `gfm_index.csv` |
| 10 | `gfm_evidence.py` | per unit `gfm_evidence.csv` (and `missing_barrier_hints.gpkg`); merged `gfm_evidence.csv`, `gfm_evidence.gpkg`, `missing_barrier_hints.gpkg` |
| 11 | `download_glofas.py` | `glofas_discharge.npz`; yearly downloads in `work/<c>/glofas/` |
| 12 | `rp_bounds.py` | `compartments_rp.gpkg` updated: `glofas_*` match, observed bounds, combined `rp`, `*_prior` model values |
| 13 | `validate.py` | `validation.json`, `validation_coverage.gpkg`, `validation_rp_confusion.csv` |

Run logs go to `output_dir/<country>/logs/<script>_<timestamp>.log`.

**Notes on the order:**
- `compartments.py` uses `rivers_sources.gpkg` (`is_source = True` only) and
  `osm_pumping.gpkg` when they exist; without them every river feature is a source and
  `has_pumping` is False.
- **`compartments.gpkg` columns:** `compartment_id`, `area_m2`, `protection_level_m`,
  `margin_median_m`, `margin_max_m`, `river_ref_median_m`, `entry_x`, `entry_y`,
  `source_x`, `source_y`, `truncated`, `touches_corridor_edge`, `touches_dtm_edge`,
  `has_pumping`, `no_levee`, `unit_id`, `unit_label`, geometry.
- **Added columns:** `touches_dtm_edge` (B3) and `unit_label`, which points to the
  compartment's cells in the unit's `labels.tif`.

## 5. Method: from levee lines to protected compartments

### 5.1 Principle

A levee protects an area if water from the river has to rise **higher** to reach that
area with the levee than without it. Both cases are computed on the same
height-above-river surface:

* `hand_undef` (Height Above Nearest Drainage on the *undefended* DTM): the elevation of
  a cell minus the elevation of the river cell its drainage path reaches (`river_ref`).
* **Spill level** of a cell: the lowest level above the river at which water starting at
  any river cell reaches the cell, moving between 4-neighbours. It is computed by a
  priority flood (Dijkstra with the max-operator):
  `spill(n) = min over paths of max(surface along the path)`, from 0 at river cells.
  Cells leave the queue in increasing spill order, so each cell is final at first
  removal.
* `spill_undef`: the flood over `hand_undef`. `spill_def`: the same flood over the
  *defended surface*, which equals `hand_undef` except on barrier cells, where it is
  `max(crest − river_ref, hand_undef)`.
* `margin = spill_def − spill_undef ≥ 0`. A cell is **protected** if
  `margin ≥ min_margin_m` and `spill_undef ≤ max_level_m`.
* 4-connected protected regions are compartments; regions smaller than `min_area_m2` are
  dropped. **`protection_level_m`** is the minimum `spill_def` in the compartment: the
  level above the river at which water first enters.

Consequences that look odd at first but are correct:
- **Upslope parts are not protected.** Parts of a ring that lie higher than its lowest
  crest are not protected, because water would reach them at their own height anyway.
- **Nested rings.** The inner compartment of nested rings has the higher of the two
  crests as its protection level (tested in both orders).

### 5.2 Barriers

1. **Merging.** All levee files are read, reprojected and split into single lines;
   `source` holds the file stem. Bridging and line ids are computed **once per country**
   (`work/<c>/barrier_lines.gpkg`), so all units see identical barriers.
2. **Bridging** (vectorised with a spatial index):
   - **Rule.** Each endpoint with another line within `gap_tolerance_m` is joined to the
     nearest such line: to its endpoint if one is within the tolerance, otherwise to its
     nearest point (T-junctions).
   - **Already connected.** Endpoints within `touch_tolerance_m` of another line count as
     connected.
   - **Self-closing.** A single line whose own ends are within the tolerance closes on
     itself if it is longer than `self_close_min_length_factor × gap_tolerance_m`.
3. **Rasterization** on the working grid uses `ALL_TOUCHED`. Diagonal-only contacts are
   then closed by adding the higher corner cell. A 4-connected flood cannot pass between
   two cells touching only at a corner anyway; the closing also guarantees that a river
   crossing a staircase levee cuts a barrier cell.
4. **Crest elevation** is the maximum of the *native* DTM near the line, because line
   geometry can be offset from the real crest:
   - **Where native pixels are read.** Only in boxes of `crest_box_size_m` around the
     levees; boxes in one row separated by at most one empty box are read together.
   - **How each box is read.** `gdal.Warp` with `resampleAlg = "max"` and **overviews
     disabled**, because `AVERAGE` overviews would lower the crest.
   - **Search radius.** A circular maximum filter with `crest_search_radius_m` follows.
     The effective search radius is up to one working cell larger than configured, which
     favours recall.
5. **Bridged cells** get crests interpolated linearly along the bridge between the crests
   at its ends.
6. **Rivers crossing barriers.** By default the river wins: crossing cells are not
   barrier. With `gated_crossings = true` crossings stay closed. Since Stage B, closures
   are decided per river feature from OSM (§8), which replaces this blanket switch.

### 5.3 Undefended DTM

All barrier cells, dilated by `barrier_removal_width_m`, are set to nodata and
interpolated with `gdal.FillNodata`. **River cells** inside that zone keep their
elevation, because they are the river bed (B2, `preserve_river_cells`). Only river cells
that cross a levee are interpolated, because there the DTM shows the levee body.

### 5.4 HAND on the undefended DTM

1. **Stream cells** are the water-source rivers rasterized along their centre line (not
   all-touched), so that rivers cut as few barrier cells as possible.
2. **Flow routing** runs on a *conditioned copy* of the undefended DTM:
   - **Tapered stream burning** (B1): river cells are lowered by `stream_burn_m`,
     decreasing linearly to 0 at `stream_burn_buffer_m` (chamfer distance). This gives
     every channel cell a gradient towards the rasterized river; without it, D8 paths
     run along the channel and reach the river only at its lowest point (§15.1).
   - **Complete breaching:** a priority flood from the river cells. A pit is carved back
     towards the outlet and then processed through a plain queue, because it lies below
     every queued cell. Flats get an ε gradient (1e-5 m).
   - **D8 steepest descent** on the conditioned surface. Receivers are always lower, so
     paths are acyclic; the river cell reached is resolved in O(n) by walking and
     assigning paths.
3. **`river_ref`** is the elevation of the river cell reached, from the *unconditioned*
   undefended DTM. `hand_undef = dtm_undefended − river_ref`.
4. **Fallback.** Cells not connected to any river cell take the nearest river cell; their
   number is logged.
5. **Only modelled water-source rivers are outlets**, so side catchments are routed
   (carved) to a modelled river. Every cell gets the reference of a river that the hazard
   maps contain.

### 5.5 Nodata in the DTM: coverage edge versus water surfaces

Lidar DTMs often have nodata on water surfaces. The two kinds of nodata must be treated
differently:

* **Missing coverage** (beyond the survey area) is unknown terrain. It stays nodata and
  acts as a wall. Compartments next to it are kept (recall) and flagged
  **`touches_dtm_edge`** (B3).
* **Gaps inside the coverage** (water surfaces, data gaps) are interpolated. Otherwise a
  river surveyed as nodata would have no river cells, no water sources and therefore no
  compartments at all.

**Telling them apart.** The coverage is the **morphological closing** of the valid DTM
cells with radius `max_fill_distance_m`, computed with two chamfer distance passes, so
the cost does not depend on the radius:
- **What it does:** it closes gaps narrower than twice the radius, including a gap where
  a river crosses the coverage edge, but does not move a straight coverage edge.
- **Window border:** the data mask is padded with "no data" before the closing, so a
  coverage edge near the window border is not closed against the border.
- **Search distance:** gaps are filled by `gdal.FillNodata` with a search distance equal
  to the largest gap half-width.

**Water surfaces** (`cap_river_holes`), i.e. gaps crossed by a river:
- **River cells** are interpolated **along each river line** between the nearest valid
  river cells (the river profile).
- **Where no river cell is valid in the window** (the whole river is nodata), and for the
  other cells of such gaps, the fill is capped at the lowest valid cell within reach:
  the waterline. Interpolation from the banks alone would lie above the water level;
  in a test this raised `river_ref` by about 1 m and produced a spurious compartment.

Test T3 makes the whole river nodata and checks that the result equals the one with a
complete DTM.

### 5.6 Seeds, entry and source

The priority flood also returns, for every cell:
- **the seed** (river cell) where its minimal spill path starts,
- **the order** in which cells were reached.

For each compartment:
- **`entry_x`, `entry_y`:** the compartment cell reached first by the defended flood,
  i.e. where water enters. Taking the cell with the minimal `spill_def` would be
  ambiguous, because every cell below the crest has the same spill level; the reach
  order resolves this.
- **`source_x`, `source_y`:** the seed of the entry cell, i.e. the river cell whose water
  enters first. Stage C matches it with GloFAS upstream areas.

### 5.7 Flags

| Flag | Meaning |
|---|---|
| `truncated` | Touches the edge of the unit window, where it may continue (§7.3) |
| `touches_corridor_edge` | Adjacent to cells outside the processing domain (corridor + halo) |
| `touches_dtm_edge` | Adjacent to missing DTM coverage, which closes it as a wall |
| `has_pumping` | A pumping station from `osm_pumping.gpkg` lies inside or within `levee_adjacent_distance_m` |
| `no_levee` | No barrier line within `levee_adjacent_distance_m` of its boundary |

## 6. Large DTMs: corridors, storage layout, units

National DTMs of 150–600 GB are never loaded whole. All reads go through `gis.Raster`
window functions on processing units.

### 6.1 Footprint

`Raster.footprint` reads the DTM decimated to `footprint_resolution_m` in a single call:
- **with overviews**, GDAL serves the read from the closest overview;
- **without overviews**, it decodes only the strips or tiles that contain the sampled
  rows and columns.

The valid mask is dilated by one coarse cell (recall) and polygonized. On the benchmark
DTM in strip layout this took 0.1 s per file.

### 6.2 Corridor

Corridor = rivers buffered by `corridor_buffer_m`, intersected with the footprint.
Compartments far from the modelled rivers cannot be flooded by them.

### 6.3 Storage layout and conversion

`prepare_dtm.py` reports for every DTM:
- size, pixel size, data type, nodata, compression, block size, overviews, file size;
- a **benchmark**:
  - 20 random 4096² native windows inside the footprint;
  - the same ground windows with `gdal.Warp` to the working resolution;
  - one sequential full-width read;
  - the compression ratio and speed of the target layout, measured on one window.

**Overview use is verified, not assumed.** An extra single-threaded warp with
`CPL_DEBUG` reports the source window the warp kernel reads, and a source window smaller
than the native one proves that an overview was used. A multithreaded warp cannot be
inspected this way, because its debug messages come from worker threads, whose error
handlers are separate per thread.

**Why layout matters.** In a strip layout (block height below 64, or block width equal
to the raster width) without overviews, every warp decodes every native pixel. Measured
per 30 M-cell unit (§13): terrain read **8–9 s** in strip layout, **1.7–2.3 s** after
conversion option A.

Conversion only runs with `convert_dtm = true`. Both options read the source once,
sequentially, and build `AVERAGE` overviews at factors 2–32 with the same compression and
predictor:
- **Option A:** the whole file as a tiled GeoTIFF (512², ZSTD or DEFLATE, predictor 3
  or 2, BigTIFF, all CPUs).
- **Option B:** only the corridor **plus the halo**, written in row blocks sized by
  `warp_memory_limit_mb`.

Before any conversion the script prints, from the measurements:
- **disk space:** raw size × measured compression ratio × 4/3 for overviews; for B, also
  × the domain share;
- **time:** sequential read time plus compression time.

It also checks the free disk space and stops if it is insufficient.

### 6.4 Units

- **Cores.** Square, `unit_core_size_m`, aligned to the working grid anchored at the DTM
  origin. Only cores that intersect the corridor are kept. Ids are `<column>_<row>`.
- **Windows.** Each core plus `unit_halo_m`, clipped to the DTM extent.
- **Domain.** Inside the window, cells outside the corridor plus halo are nodata and
  skipped by every algorithm.

**Halo rule: `unit_halo_m ≥ corridor_buffer_m`** (checked in the configuration).
1. A core cell can lie up to `corridor_buffer_m` from its river.
2. The window extends `unit_halo_m` beyond the core in every direction.
3. So the river point nearest to every core cell is inside the window only if
   halo ≥ buffer.
4. Without it, HAND has no valid reference; units whose window contains no river are
   listed.

**Memory model.** RAM per unit ≈ window cells × `bytes_per_cell`.
- **Check:** `build_units.py` stops if the largest window would exceed
  `max_ram_gb_per_worker` and prints the largest core size that fits.
- **Measured value:** `compartments.py` logs the measured bytes per window cell, now
  **93–95 B**.

## 7. Parallel processing, checkpoints, merging

### 7.1 Pool

`gis.Run.run_units` uses a `ProcessPoolExecutor` with the `spawn` start method and
`max_tasks_per_child = 1`:
- **A fresh process per unit.** Peak RSS is exact per unit, memory is returned to the
  system after each unit, and Linux and Windows behave the same.
- **Per-process GDAL handles.** Each worker opens its own GDAL datasets, because GDAL
  handles cannot be pickled, and sets `GDAL_CACHEMAX` and `GDAL_NUM_THREADS` from the
  configuration.
- **Compilation once.** The numba kernels are compiled once in the main process, so
  workers load them from the cache.
- **Submission order.** Units are submitted largest first; the ETA is based on the domain
  cells finished so far.
- **Failures.** A crashed worker, e.g. one killed for lack of memory, is reported; a rerun
  retries only units without a done marker.

### 7.2 Checkpoints

- **`compartments.py`:** each unit writes to `work_dir/<country>/units/<unit_id>/` and
  finishes with `done.json` (timings, peak RSS, counts), written atomically through a
  temporary file and a rename.
- **`rp_index.py`:** writes `rp_index.json` per unit and recomputes it when the unit's
  `done.json` is newer.
- **Reruns** skip finished units unless `FORCE = True`; `UNITS` restricts a run to the
  listed units.
- **Merging** always uses all finished units.

### 7.3 Merging

The same compartment can be found by several units, because windows overlap by twice the
halo.

1. **Grouping.** Candidates from different units that overlap by at least
   `merge_overlap_ratio` of the smaller one form a group (connected components).
2. **Groups with a complete candidate.** Keep the complete (not truncated) candidate
   whose representative point lies in its own unit's core; if several qualify, keep the
   one farthest from its window edge. A ring across a core boundary is kept exactly once.
3. **Groups truncated in every window** (compartments longer than a window):
   - the pieces are **united** (recall);
   - level, entry and source come from the piece with the **lowest** protection level. A
     window edge can only overestimate protection, because it acts as a wall.
4. **Contradiction rule** (`drop_contradicted_truncated`):
   - **What it catches:** window edges can create artificial protection, for example an
     area behind a long levee that looks enclosed because the window edges close it.
   - **Rule:** a compartment truncated in every window is dropped when finished units
     whose windows cover at least `contradiction_ratio` of it found nothing there. Those
     units see more of the levee.
5. **Identifiers.** `compartment_id = <country>_<x>_<y>`, from the representative point
   rounded to `id_rounding_m`, stable between runs. Duplicates get a suffix.

Barrier crest statistics in `barriers.gpkg` are sampled inside each unit's **core** only
and combined as the minimum.

## 8. Water sources and pumping stations from OpenStreetMap

### 8.1 Extraction (`download_osm.py`)

1. **Download.** The country `.osm.pbf` is downloaded from the URL template
   (`{region}` = `geofabrik_region`). It is skipped when the local file has the size the
   server reports. Alternatively a local extract can be set in `osm_pbf`.
2. **Custom `osmconf.ini`.** `config/osmconf.ini` (GDAL's default plus the tag keys used)
   exposes `waterway`, `man_made` and related keys as columns in the points, lines and
   multipolygons layers. The script checks that every key of the configured tag lists is
   exposed.
3. **Extraction.** One `gdal.VectorTranslate` call per output, each writing a single-layer
   GeoPackage:
   - **`osm_hydraulic_points.gpkg`:** SpatiaLite SQL over the three OSM layers
     (`UNION ALL`). Ways and areas become representative points (`ST_PointOnSurface`);
     `osm_id` and the tag columns are kept, and `osm_element` says node, way or area.
   - **`osm_drains.gpkg`:** the drain lines.
   - **Spatial filter:** the DTM footprint's bounding box in WGS 84. GDAL does not accept
     a filter CRS together with SQL, so the box is transformed first.
4. **Counts per tag value** are logged, to tune the tag lists.

### 8.2 Water sources (`osm_hydraulics.py`)

A river feature is **not a water source** (`is_source = False`, with a `reason`) when:
- **a closing structure** (`osm_tags_closing`, by default sluice and lock gates) lies
  within `gap_tolerance_m` of one of its crossings with a levee line, or
- **a pumping station** (`osm_tags_pumping`) lies within `gap_tolerance_m` of its
  **downstream end**.

**Downstream end.** It is taken from **DTM elevations**: the median of the native DTM
within `endpoint_elevation_radius_m` at both line ends, and the lower end is downstream.
The digitising direction of river lines is not reliable. The line direction is used only
when the elevations are equal or missing, and `downstream_method` records which source
decided.

**Closed streams are excluded** both as water sources and from the drainage network for
HAND. Otherwise HAND inside a compartment would refer to an internal drain instead of the
main river. A tributary that leaves a levee ring through a gate therefore no longer
floods the ring (test R1); without the OSM step it does, because the river wins.

**`osm_pumping.gpkg`** holds the pumping stations within `levee_adjacent_distance_m` of a
levee or river, with the flags `near_levee` and `near_river`.

**Granularity.** The decision applies to the whole river feature as delivered. A feature
that runs far beyond the levee is closed along its whole length, so river networks
should be split at confluences.

## 9. Return period index

### 9.1 Exposure value (`rp_index.py`, per unit in parallel)

1. **GHSL building volume.** `download_ghsl.py` downloads the total and non-residential
   volume tiles (100 m, Mollweide) that intersect the DTM footprint, found with the
   tiling schema, and mosaics each product as a VRT.
   - `ghsl_epoch = 0` takes the latest epoch the server provides that is not after the
     current year; it is found by probing.
   - Epochs after the present minus five years are logged as modelled.
2. **Value per cell:**
   - residential volume = total − non-residential, clipped at 0;
   - floor area = volume / `storey_height_m`;
   - value = floor area × maximum damage per m² (`max_damage_csv`: columns `country`,
     `residential`, `non_residential`; one row whose `country` equals the country code).
3. **Totals-preserving warp.** The source average (volume per source cell) times the
   ratio of working to source cell area is the volume per working cell; Mollweide is
   equal-area.
   - **Totals check, every unit:** the warped total is compared with an independent total,
     the source cells weighted by the **exact** share of their area inside the window
     (reprojected window ∩ cell boxes).
   - **Measured:** within 0.03–0.05 %, which is about the scale distortion of the national
     grid. The configured tolerance is 1 %.
4. **Only the bounding box of the unit's kept compartments** is valued: 2–5 s per unit
   instead of 15 s.
5. **Damage:**
   - depth per cell = max(0, `protection_level_m` − `hand_undef`);
   - damage fraction from `depth_damage_csv` (`depth_m`, `res`, `nonres`), interpolated
     linearly;
   - `value_exposed` V = Σ value × fraction over the compartment.

### 9.2 Score, classes, cap, interval (whole country)

1. **L**, `levee_length_m`: length of barrier lines within `levee_adjacent_distance_m` of
   the boundary. L = 0 sets `no_levee`; such compartments get no score.
2. **A**, `upstream_area_km2`: from `river_upstream_area_field` of the river feature
   nearest to the source cell, or from `upstream_area_raster` at the source cell.
   Compartments without A get a score from the value term only, and they are counted in
   the log.
3. **Score** = `w_v` × percentile rank of log(V/L) + `w_a` × percentile rank of log(A),
   ranked within the country. `v_rel` = V/L. V = 0 ranks lowest.
4. **Quantile mapping** (`gis.Stats.quantile_classes`): ascending scores get ascending
   `distribution_classes` with the configured shares; tied scores get the same class.
5. **Cap:**
   - **Rule.** `rp_cap` is a step function of `protection_level_m`: the first
     `cap_levels_m` value **above** the level selects the matching `cap_rp`. Above the
     last level there is no cap.
   - **Example.** With levels 0.5, 1, 2 and caps 10, 50, 100: a level of 0.3 m gives
     cap 10, a level of 1.4 m gives cap 100.
6. **RP:** `rp = min(rp_model, rp_cap, max_rp)`.
   - `rp_source = "cap"` when the cap or `max_rp` lowered the model value, or when there
     is no model value (`no_levee`); otherwise `"model"`.
   - A `no_levee` compartment without an applicable cap therefore gets `max_rp`, as the
     specification says ("from the cap only"). See the open decisions (§15.4).
7. **Interval:**
   - `rp_low` = the previous value in `classes` (the same value if `rp` is the lowest);
   - `rp_high` = the next value in `classes`, bounded by `rp_cap` and `max_rp`, never
     below `rp`.
8. **Overrides.** A compartment containing a point of the overrides file gets its `rp`,
   with `rp_low = rp_high = rp` and `rp_source = "override"`. If points disagree, the
   first is used and a warning is logged.
9. **Stage C columns.** `obs_rp_low`, `obs_rp_high`, `n_held`, `n_flooded`,
   `obs_conflict` and `model_obs_conflict` are present and empty; `rp_bounds.py` fills
   them (§11). `value_exposed` is an added column.

### 9.3 Land cover exposure (`HAND_YZ.py`)

GHSL measures buildings. Assets outside buildings carry no value in §9.1:
- crops;
- roads, railways and port areas;
- open-air equipment of industrial sites.

`HAND_YZ.py` adds them from land use and land cover data. The resulting value is added to
the building value: it **supplements** GHSL and never replaces it.

1. **Data.**
   - **ESA WorldCover 10 m** (2021, v200) everywhere. `download_worldcover.py` downloads
     the tiles of the official 3° tile grid that touch the corridor and builds
     `worldcover.vrt`. Local tiles can be listed in `lulc.worldcover` instead.
   - **Urban Atlas** where you provide it (manual download). It takes precedence wherever
     its polygons lie, also with classes that carry no value, because it is finer and
     distinguishes land use.
   - **Multi-layer files.** The configured layer (`urban_atlas_layer`) of a multi-layer
     file is copied to its own GeoPackage with `gdal.VectorTranslate`.
2. **Categories.** `worldcover_categories` and `urban_atlas_categories` map class codes
   to exposure categories such as `agriculture`, `built_up` and `transport`. Each category
   has two inputs:
   - a value per m² of land (`lulc.max_damage_csv`, columns `country` + categories);
   - a depth-damage curve (`lulc.depth_damage_csv`, `depth_m` + categories).
3. **Damage.** Depth per cell is `max(0, protection_level_m − hand_undef)`, as in §9.1.
   Value = cell area × value per m² of the cell's category; damage = value × the
   category's damage fraction at that depth. Sums per compartment and per category go to
   `exposure_lulc.csv`: `value_lulc`, `lulc_exposed_area_m2`, `lulc_<category>`.
4. **Grid.** WorldCover classes are read with nearest neighbour onto the working grid.
   Its pixels (about 6–9 m in Europe) are as large as or larger than the working cells
   (2–5 m). Urban Atlas polygons are rasterized by cell centre.
5. **Combination in `rp_index.py`.** With `lulc.enabled = true`:
   `value_exposed = value_ghsl + value_lulc`. The score then uses this sum, and both
   parts stay in `compartments_rp.gpkg`.
6. **No double counting.** The values per m² must cover only what GHSL does not hold.
   - **WorldCover `built_up` class:** it contains buildings, so its value should cover
     roads, utilities and open-air assets only.
   - **Unknown values:** set the value to 0 if unsure.
7. **Checkpoints.** Each unit stores a hash of its inputs in `hand_yz.json`: values,
   curves, mappings, file dates of the WorldCover tiles and Urban Atlas files. A rerun
   recomputes exactly the units whose inputs changed (tested).

**Verified and unverified sources.**
- **WorldCover, verified:**
  - URL pattern and tile grid from the official ESA-WorldCover download script;
  - a request to a real tile (HTTP 200);
  - class codes and licence (CC BY 4.0) from the tile's own metadata: 10 tree cover,
    20 shrubland, 30 grassland, 40 cropland, 50 built-up, 60 bare/sparse vegetation,
    70 snow and ice, 80 permanent water, 90 herbaceous wetland, 95 mangroves, 100 moss
    and lichen.
- **Urban Atlas, not verified:** the class field and code list could not be checked,
  because the EEA hosts were blocked in this session. A catalogue search shows the field
  `code_2018` and codes such as 11100 (continuous urban fabric), 12100 (industrial,
  commercial, public, military and private units) and 12210 (fast transit roads). The
  configuration entries are therefore empty.

## 10. Observational evidence: Copernicus GFM

The model RP of §9 is an *a priori* estimate built from exposure and topography. Stage C
adds what was actually observed: whether a compartment stayed dry while its surroundings
were flooded (Global Flood Monitoring, GFM, Sentinel-1 SAR flood extents since 2015), and
how rare the discharge was on that day (GloFAS, §11). Observations bound the RP; they never
create, move or delete a compartment.

### 10.1 Scene search and screening (`download_gfm.py`)

1. **One search for all units.** The densified unit cores are transformed to WGS84. Their
   common box goes to the STAC item search `GET {gfm_stac_url}/search` with `collections`,
   `bbox`, `datetime` (`gfm_date_from` to `gfm_date_to`; today when empty) and
   `limit = gfm_page_limit`.
   - The script follows `rel="next"` links to the last page.
   - A link that declares `method: POST` is sent as POST with its JSON `body`, because STAC
     APIs page either way.
   - Requests go through `gis.Http.get_json` (urllib, retries, exponential backoff); no STAC
     client library is needed.
2. **Token.** If `gfm_token_env` names an environment variable, its value is sent as
   `Authorization: Bearer …` with the search and with GDAL's HTTP reads
   (`GDAL_HTTP_HEADERS`). The token never appears in the configuration or the logs.
3. **Scenes to units.** Each item is assigned to the units whose core its footprint
   (`geometry`, else `bbox`) touches.
4. **Screening: one read per scene, at `gfm_screen_resolution_m` (100 m).** The flood layer
   is read once through `/vsicurl/`, with `average` resampling, over the union of the
   windows of all units the scene covers.
   - **Shared grid.** All unit windows are snapped to one grid of screening cells, so each
     unit's part is a plain slice of the scene's array.
   - **Cost.** A cloud-optimised GeoTIFF serves the read from its overviews: kilobytes
     instead of the 20 m data.
   - **Flooded share.** The mean of the 0/1 layer per screening cell is the flooded share;
     cells with a `gfm_nodata_values` code are left out of the mean.
   - **Flooded area.** The shares are summed over each unit's domain and converted to
     working cells.
   - **Skip rule.** A unit window with at most `min_valid_pixels` flooded working cells
     cannot produce a valid test, so it is skipped. Most acquisitions over a river are dry,
     so only the few scenes that matter get copied.
   - **Why one read per scene, not per unit and scene.** A country has hundreds of units
     and each unit about a thousand Sentinel-1 acquisitions since 2015. Reading per unit
     would mean hundreds of thousands of latency-bound requests; reading per scene cuts that
     by the number of units a scene covers.
5. **Clipped copies.** For each kept pair of scene and unit, the flood, exclusion and
   reference water layers are cut to the unit window with
   `gdal.Translate(projWin, projWinSRS)`. The optional likelihood layer is cut too when the
   item has it.
   - **Format.** Tiled, compressed GeoTIFFs in `work/<c>/gfm/<unit>/`.
   - **Grid.** The copies stay in the GFM grid and CRS. Reprojection happens once, in
     `gfm_evidence.py`, directly onto the unit grid; resampling twice would shift the 20 m
     pixel edges.
6. **Threads.** `gfm_download_threads` scenes are screened and copied at the same time.
   The work waits on the network, and GDAL releases the GIL while it waits, so threads are
   enough; no processes are needed.
7. **Checkpoint and incremental updates.**
   - Every screened pair of scene and unit, kept or not, is appended to
     `work/<c>/gfm/screening.csv`.
   - A rerun screens only new pairs, so running the script again a year later fetches only
     the new acquisitions. `FORCE` screens the selected units again.
   - A scene that cannot be read is logged and retried in the next run.
8. **Index.** `gfm_index.csv` lists the kept pairs: `unit_id`, `item_id`, `datetime`,
   `flooded_pixels` and the layer paths.

### 10.2 Tests per compartment and scene (`gfm_evidence.py`)

A dry compartment on a dry day proves nothing. Each scene therefore answers a comparative
question: **was the compartment dry while the undefended floodplain at the same height next
to it was wet?** Only that comparison, made within the same scene, turns a dry compartment
into evidence.

1. **Reading.** The three layers are warped with nearest neighbour onto the unit's working
   grid, which is the exact grid of `labels.tif` and `hand_undef.tif`.
   - **Valid cell:** it has data (`gfm_nodata_values` count as nodata even if the file
     declares none), is not excluded, is not reference (permanent) water, and has a HAND
     value.
   - **Wet cell:** valid and carrying one of the `gfm_flood_values`.
2. **Inside:** the valid cells of the compartment (its `unit_label` in `labels.tif`), giving
   `inside_valid_px` and `inside_wet_frac`.
3. **Comparison zone,** giving `comp_valid_px` and `comp_wet_frac`. It holds valid cells that
   meet all of these:
   - outside all compartments;
   - within `comparison_distance_m` (chamfer distance from the compartment);
   - `hand_undef` at most the median `hand_undef` inside the compartment;
   - `hand_undef` at most `max_level_m`.

   **Why the HAND condition.** It restricts the comparison to terrain that an undefended
   flood up to the compartment's median level would also reach. If those cells are wet,
   water stood at least that high next to the levee.
4. **Outcome:**
   - `no_test` unless both pixel counts reach `min_valid_pixels` and
     `comp_wet_frac ≥ f_wet_comparison`;
   - `held` if `inside_wet_frac ≤ f_dry_inside`;
   - `flooded` if `inside_wet_frac ≥ f_wet_inside`;
   - otherwise `inconclusive`.
5. **Efficiency.**
   - **Read extent.** Without hints, each scene is read only over the bounding box of the
     unit's kept compartments plus the comparison distance.
   - **Zones.** Each compartment's zone is computed once, in its own bounding box, and
     reused for every scene.
   - **Parallelism.** Units run in parallel with the same pool as `compartments.py`.
6. **Outputs:**
   - **`gfm_evidence.csv`:** one row per compartment and scene.
   - **`gfm_evidence.gpkg`:** the same rows with the compartment polygon, for QGIS. Filter
     by `item_id` or `datetime` and symbolise by `outcome`.
7. **Checking a named flood event.** Restrict the run to the event's dates:

   ```python
   DATE_FROM = "2024-09-12"   # first and last day of the event
   DATE_TO = "2024-09-25"
   ```

   The result is the event's evidence table plus a GeoPackage for visual inspection.
   Without the filter, all scenes in the index are tested.
8. **Consistency.**
   - The unit rasters must come from the same `compartments.py` run as `compartments.gpkg`,
     because `unit_label` points into them. After rerunning `compartments.py`, rerun this
     script.
   - Thresholds `f_*` are placeholders. With 0.05 and 0.3, "dry" tolerates 5 % wet cells
     (SAR noise, local ponding), and "flooded" requires 30 % of the valid area to be wet.

### 10.3 Missing-barrier hints (`HINTS = True`)

1. **Scene water level in HAND.** It is the `river_level_percentile` (90 %) of `hand_undef`
   among wet valid cells **connected to the river**.
   - **Connectivity** runs through wet cells, reference water and the rasterized river
     lines.
   - **Why reference water counts.** The river is permanent water and never "wet" in the
     flood layer. Labelling wet cells alone would cut the flood off from its own river.
2. **Tested cells.** A cell is tested if it meets all of these:
   - valid and outside compartments;
   - at least `hint_exclusion_m` from barrier lines and compartments;
   - `hand_undef < level − hint_margin_m`.

   Each scene adds to the cell's `test_count`, and also to its `dry_count` where the cell
   stays dry.
3. **Hints.**
   - **Selection:** cells with `test_count ≥ hint_min_scenes` and
     `dry_count / test_count ≥ hint_dry_share`.
   - **Clusters:** the selected cells are clustered, clusters of at least `min_area_m2`
     are kept, and they are polygonised into `missing_barrier_hints.gpkg` with `area_m2`,
     `scenes_tested_mean` and `dry_share_mean`.
   - **Units:** as with compartments, each hint is taken from the unit whose core holds its
     representative point.
4. **`hint_exclusion_m` (an addition to the specification).** A levee that holds is itself
   dry, and so are its toes and the 20 m GFM pixels that mix the levee with the floodplain.
   - **The problem.** These cells lie below the water level in `hand_undef`, because the
     levee was removed from that surface. Without the exclusion, every levee that held
     would show up as a "missing barrier" along its own body.
   - **The test shows it.** Without the exclusion, the levee bands of two rings become
     hints (3 areas instead of 1).
   - **Default.** 40 m = two GFM pixels.
5. **Limitation: HAND is relative to the drainage cell.**
   - **Backwater.** When a large river backs up into a tributary valley, cells along the
     tributary have low HAND relative to the tributary, while the main river sets the
     flood level. Such cells can stay dry legitimately and still appear as hints.
   - **Consequence.** Hints are candidates for inspection, not results, and they change no
     compartment.
   - **Better level, not implemented.** A level per drainage target would need the
     `breach_to_streams` target index stored per unit, which `compartments.py` does not keep.

## 11. Observed bounds of the return period: GloFAS

### 11.1 Discharge data (`download_glofas.py`)

1. **Source.** GloFAS historical river discharge (reanalysis), dataset
   `cems-glofas-historical` of the CEMS Early Warning Data Store. The identifier was
   checked on the dataset page.
   - **Data.** Daily `dis24` on a 0.05° grid from 1979.
   - **Upstream area.** The static upstream area map on the same grid comes from the GloFAS
     auxiliary data; check its units attribute against `glofas_upstream_area_to_km2`
     (1e-6 for m²). Its URL is not verified, so give a local file.
2. **Two ways in:**
   - **Local files** (`glofas_files`): GRIB or NetCDF files downloaded through the EWDS form.
     This needs nothing beyond the environment.
   - **Yearly download through `cdsapi`.** Paste the request shown by the EWDS form
     ("Show API request") into `glofas_request`, without the time and area keys.
     - **What the script adds:** the year, all months and days (`glofas_time_keys`) and the
       area (`glofas_area_key`: north, west, south, east, rounded outward to 0.01°).
     - **Approval needed.** `cdsapi` is not installed. Without it, the script stops and
       explains both options.
3. **Area.** The box of the DTM footprint in WGS84, extended by `glofas_search_radius_km`,
   so that river cells next to the footprint stay available for matching.
4. **Cells kept.** Every cell of the upstream area map with at least
   `glofas_min_upstream_area_km2` (after `glofas_upstream_area_to_km2`). Their daily series
   are read with GDAL alone (no xarray):
   - **GRIB:** each band's time is `GRIB_VALID_TIME`.
   - **NetCDF:** the band's `NETCDF_DIM_<time>` value with the `<time>#units` attribute
     (days, hours or seconds since a reference).
   - **Day attribution.** A daily mean stamped at the end of its 24 h period belongs to the
     day that ends. With `glofas_valid_time_ends_day = true`, one second is subtracted
     before the date is taken.
     - **Status:** a placeholder. Check one real file against the EWDS documentation.
     - **Why it matters:** one day of shift barely changes event discharges on large rivers,
       but it does matter on flashy small ones.
5. **Output `glofas_discharge.npz`:** cell coordinates, upstream areas, dates, and discharge
   (days × cells, float32). With the default minimum area of 0 km², every cell of the box is
   kept. A box of about 80 000 km² holds about 4000 cells: 46 years × 4000 cells ≈ 270 MB in
   memory, compressed on disk.

### 11.2 Matching compartments to GloFAS cells (`rp_bounds.py`)

1. **Candidates.** For each compartment with a source cell and an upstream area A (§9.2):
   GloFAS cells within `glofas_search_radius_km` of the source cell (KD-tree in the
   country's CRS).
2. **Ratio filter.** The candidate's upstream area divided by A must lie within
   [`upstream_ratio_min`, `upstream_ratio_max`].
3. **Choice.** Among the remaining cells, the one with the ratio closest to 1 (log scale)
   wins.
4. **Stored:** `glofas_cell`, `glofas_distance_km`, `glofas_area_ratio`.
5. **No match:** no event RP.

**Why the ratio and not the distance.** On a 0.05° grid, the nearest cell is often on
another river. The ratio test selects the cell that represents the same reach. The test
puts a second cell with a ratio of 0.92 within the radius; the cell with ratio 1 wins.

### 11.3 Frequency analysis

Annual maxima of each matched cell are taken over the full series, using years with at
least `glofas_min_days_per_year` days. A Gumbel distribution is fitted by L-moments:

- b0 = mean, b1 = Σ (i−1)/(n−1) · x(i) / n over the ascending sample, l1 = b0,
  l2 = 2 b1 − b0;
- α = l2 / ln 2, ξ = l1 − 0.5772 α;
- RP(Q) = 1 / (1 − exp(−exp(−(Q − ξ)/α))).

**Why Gumbel with L-moments.** It is stable for about 40 annual maxima and has no shape
parameter to over-fit.

**Known weakness: heavy-tailed rivers.** The Gumbel tail is lighter than theirs, so a large
discharge receives a higher RP than it really has.
- **Held tests** then give lower bounds that are too high, the unsafe direction for a lower
  bound.
- **Flooded tests** give upper bounds that are too high, which only loosens them.

A three-parameter GEV would capture the tail but is unstable for single cells with about 40
maxima; regional pooling of cells would be the next step. Cells with fewer than 10 complete
years are logged.

### 11.4 Event return periods and bounds

- **held:** discharge on the acquisition day. Its RP is a **lower bound** of the protection,
  because the levee withstood at least that.
- **flooded:** the maximum discharge over `flooded_window_days` up to and including the
  acquisition day. Its RP is an **upper bound**: the compartment flooded during an event of
  at most that RP, and the window covers a peak that may have passed before the overpass.
- **Per compartment:**
  - `n_held`, `n_flooded`;
  - `obs_rp_low` = maximum over held tests;
  - `obs_rp_high` = minimum over flooded tests;
  - `obs_conflict` if `obs_rp_low > obs_rp_high`.

The SAR image says *whether* the compartment was wet. The matched discharge series says
*how rare* the day was at the same river reach.

### 11.5 Combination with the model interval

| Case | Condition | Result |
|---|---|---|
| `override_kept` | the compartment has an override | observations counted; RP unchanged |
| `obs_conflict` | `obs_rp_low > obs_rp_high` | flag only; RP unchanged |
| `combined` | low = max(`rp_low`, `obs_rp_low`) ≤ high = min(`rp_high`, `obs_rp_high`) | `rp` clamped into [low, high], then snapped to the class inside [low, high] closest to it (log scale) |
| `model_obs_conflict` | low > high: the observations contradict the model interval | the observed interval alone is used; `rp` = the class inside it closest to the model value |

- **Open side.** An open side of the interval becomes the neighbouring class of the result,
  mirroring the ±1 class interval of the model. Two cases:
  - only held tests: rp is the smallest class at or above the observation, and the upper
    side is the next class;
  - only flooded tests: the mirror image.
- **No class in the interval.** This happens with a narrow interval between two classes.
  The closest class is taken and the reported interval is widened to include it, so
  `rp_low ≤ rp ≤ rp_high` always holds. The observed values stay in `obs_rp_low` and
  `obs_rp_high`.
- **Source tag.** `rp_source` gets `+obs` (e.g. `cap+obs`).
- **Rerunnable.** The model values are stored once as `rp_prior`, `rp_low_prior`,
  `rp_high_prior` and `rp_source_prior`, and every run starts from them. A rerun with new
  evidence therefore gives the same result as a first run (tested). After rerunning
  `rp_index.py`, rerun `rp_bounds.py`.
- **Log.** Counts per case are logged.

## 12. Validation

`validate.py` compares the compartments with reference protected areas.

1. **Areas,** as in the original script:
   - recall = |result ∩ reference| / |reference| and precision = |result ∩ reference| /
     |result|, both on the unions;
   - per reference polygon: covered area and covered share (`validation_coverage.gpkg`),
     plus the counts of fully missed polygons and polygons covered above 50 %;
   - `validation_extent` optionally restricts both layers to an area where the reference
     is complete.
2. **Return periods,** when `reference_rp_field` is set:
   - **Pairing.** Each compartment is compared with the reference polygon it overlaps most.
     Reference RPs are snapped to the classes on a log scale.
   - **Confusion matrix:** `validation_rp_confusion.csv`, reference class × assigned class.
   - **Interval share:** compartments whose reference RP lies within `rp_low`…`rp_high`.
   - **Class share:** compartments with the same class.
3. **Inputs and report.** `compartments_rp.gpkg` is used when it exists, otherwise
   `compartments.gpkg`. All numbers go to `validation.json`.

## 13. Measured performance

### 13.1 Where the numbers come from

All numbers are from this cloud container: 4 vCPUs, 15 GB RAM, virtual disk. The OS page
cache was dropped before timed runs. The synthetic benchmark DTM comes from
`tests/make_benchmark_world.py`:
- two GeoTIFFs of 12 000 × 24 000 px at 1 m;
- strip layout, no overviews.

Working resolution 2 m, cores 7 km, halo 2 km. That gives windows of
5500² = 30.25 M cells, close to a real unit of 20 + 2 × 3 km at 5 m (27 M).

**Disk and decompression speeds on your machines will differ.**

### 13.2 DTM reads (20 windows per file)

| Layout | Native 4096² window | `gdal.Warp` to 2 m (same ground) | Sequential | Overview used |
|---|---|---|---|---|
| Strip, no overviews | 0.31–0.45 s (148–214 MB/s) | 0.74–0.81 s (82–90 MB/s native eq.) | 530–640 MB/s | no |
| Option A (tiled, ZSTD, overviews) | 0.14–0.20 s (342–468 MB/s) | 0.34–0.38 s (178–197 MB/s) | 368–501 MB/s | yes (factor 2) |

Conversion results:
- **Size.** ZSTD with predictor 3 compressed the synthetic DTM to 0.47 of its raw size.
- **Time.** Option A took 8–11 s per 1.15 GB file plus 5–6 s for overviews.

### 13.3 Kernels on one unit (30.25 M cells)

| Kernel | Stage A | Now | Change |
|---|---|---|---|
| HAND (burn taper, breaching, D8) | 29.3 s (exact distance transform 8.8 s of it) | 13.9 s | chamfer taper, 4-ary array heap, plain queue for pits, O(n) path resolution |
| One priority flood | 14.6 s | 10.9–11.1 s | 4-ary array heap, plain queue for cells at or below the current level, 32-bit arrays |

The spill levels are bit-identical to Stage A.

Remaining potential:
- **Floor of the flood.** The flood's floor without any heap work is 5.0 s (memory
  access).
- **Radix heap.** A radix heap could save part of the remaining ~6 s per flood; it is
  not implemented.

### 13.4 `compartments.py` per unit (30.25 M window cells, 24.4 M domain cells, cold cache, option A source)

| Step | Stage A | Now |
|---|---|---|
| Terrain read | 1.7–2.3 s | 1.9–2.2 s |
| Domain, gap filling | 2.5–3.5 s | 2.1–2.2 s |
| Crest boxes (6–13 boxes) | 0.6–1.1 s | 1.3–1.5 s |
| HAND | 31–33 s | 15.9–16.6 s |
| Two floods | 28–33 s | 21.9–26.6 s |
| Compartment statistics | 2–3 s | 2.1–4.2 s |
| **Total** | **76–81 s** | **56–61 s** |
| Peak RSS | 3.0 GB (106–108 B/cell) | **2.6–2.7 GB (93–95 B/cell)** |

Worker time is **2.4 µs per domain cell**. At 5 m, 1000 km² of domain (40 M cells) takes
about 1.6 worker-minutes. `compartments.py` logs the extrapolation from your own test
units (`performance.json`).

### 13.5 `rp_index.py` per unit

With synthetic GHSL over the benchmark world it takes **1.6–4.6 s** with 0.6–0.9 GB peak.
Warped and source volume agree within 0.05 %.

### 13.6 Stage C

**`gfm_evidence.py`** was measured on a synthetic unit of production size:
- 6000 × 5000 = 30 M cells at 5 m;
- 32 compartments;
- GFM scenes at 20 m in EPSG:3035, warped onto the unit grid.

| Run | Time | Peak RSS |
|---|---|---|
| 1 scene, no hints | 5.4 s | 1.45 GB |
| 3 scenes, no hints | 6.4 s | 1.32 GB (45 B/cell) |
| 3 scenes, hints | 22.8 s | 2.85 GB (95 B/cell) |
| 6 scenes, hints | 27.7 s | 2.87 GB |

- **Without hints:** about 5 s per unit (labels, HAND, comparison zones) plus 0.5 s per
  scene. With 20 kept scenes per unit that is 15 s; 300 units take 75 worker-minutes, about
  20 minutes with 4 workers.
- **With hints:** about 18 s plus 1.6 s per scene, since each scene is read over the whole
  window. 300 units with 20 scenes take about 4 worker-hours.
- **Memory.** With hints, a worker needs about as much memory per cell as
  `compartments.py`; without hints, half. The worker count chosen for `compartments.py`
  therefore fits.

**`download_gfm.py`** depends on the service, which was not reachable from this environment.
Its cost is one screening read per scene plus the copies of the kept unit windows, with
`gfm_download_threads` requests at a time. On local files the test needs 5 s for 5 scenes,
4 units and 16 copied windows.

**`download_glofas.py`:**
- **Test.** 41 years × 63 cells from NetCDF take 3 s; 3 years from GRIB take 1.3 s.
- **Real data.** 46 yearly files of 365 bands with one small window per band take minutes.
  The EWDS queue, not the reading, sets the pace of a download.

**`rp_bounds.py` and `validate.py`** take about 1 s on the test country. Their loops run
over compartments and tests (KD-tree queries, one Gumbel fit per matched cell), so 10 000
compartments take seconds.

## 14. Tests

- **`pixi run test-kernels`, 16 checks:**
  - closed ring, gap, diagonal barrier;
  - nested rings in both crest orders;
  - seed and order propagation, entry cell;
  - breaching and nodata;
  - B1 reproduced without burning and fixed with it.
- **`pixi run test-pipeline`, 40 checks:** synthetic countries in two strip-layout files.
  - **Pipeline mechanics:** inspection, footprint, benchmark; conversions A and B;
    units; unit selection and checkpoints.
  - **Compartments and merging:** a seam ring kept once; bridged and wide gaps;
    staircase levee; nested rings to 0.01 m; a compartment longer than any window,
    merged, and the dropped window-edge artifact; unique ids; entry and source cells.
  - **Bug fixes:** B2, a levee on the river bank; B3, a U-shaped levee closed by missing
    coverage, kept and flagged; B4 in T2, a gated tributary entering from missing
    coverage.
  - **Gaps:** a river gap (river profile, no artifact); T3, the whole river without DTM
    data, giving the same compartments.
  - **T0, empty levees:** margin exactly 0.
- **`pixi run test-stage-b`, 20 checks:**
  - **OSM extraction:** from a `file://` template; nodes and ways as points; counts per
    tag.
  - **Water sources:** a gate at a levee crossing and a pumping station at a downstream
    end found from the DTM, so the rings are protected only after the OSM step;
    `has_pumping`; selection of pumping stations.
  - **GHSL:** zipped tiles and schema; latest-epoch probing.
  - **Value warp** totals within 1 %.
  - **RP:** exposure inside and outside a town; cap, `min` rule and class values;
    interval; override; `v_rel`.
  - **Quantile mapping:** shares within 1/N; ties.
- **`pixi run test-hand-yz`, 7 checks,** on the Stage B country:
  - **`download_worldcover`:** tiles chosen from the grid (the corridor crosses 15° E,
    so two tiles), VRT, licence;
  - **land cover value** per compartment equal to an independent computation from the
    unit rasters;
  - **Urban Atlas:** from a two-layer file; takes precedence (tree cover valued as
    transport, an unvalued class over half of ring 2);
  - **checkpoint hash:** an unchanged rerun skips every unit, doubled values recompute
    them;
  - **`rp_index`:** `value_exposed = value_ghsl + value_lulc`.
- **`pixi run test-stage-d`** (step D1: 10 checks), on the Stage B country under the code T4D:
  - **DSM:** inspected and benchmarked like the DTM;
  - **OSM:** buildings, water, crossings (a footway is excluded), walls, industrial areas (a
    closed `man_made=works` way counts as an area), candidate wall tag counts;
  - **GHSL:** SMOD and population of the probed epoch;
  - **CORINE:** the legend mapping is verified on a known pixel, and a wrong code stops the
    script;
  - **Urban Atlas:** taken from a two-layer file;
  - **context zones:** precedence of the order, small zones falling through to later
    contexts, rural as the rest, no overlaps, the corridor tiled exactly;
  - **`zones_to_grid`** on a unit grid.
- **`pixi run test-stage-e`** (E1: 16, E2: 6, E3: 7 checks), on its own synthetic country T5 (8 × 3.2 km, 1 m DTM, two levee rings):
  - **kernels:** least-cost path through the only gap of a wall; Chaikin keeps its end points, moves no
    vertex more than `max_shift` and keeps a sharper corner; hillshade with nodata 0;
  - **snapping:** a valley minimum 15 m off the line is found within one cell (2 m); a ditch 1 m deeper
    than the channel but beyond `max_offset_m` is not chosen; a line drawn upstream is flipped; a wide
    river with a nodata water surface follows the water polygon centre; a zigzag channel twice the length
    of its line keeps the original line (quality check); only `is_source` rivers are snapped;
  - **outputs:** two chunks stitched into one line; river points every `step_m` with chainage 0 at the
    downstream end and tangents pointing downstream;
  - **pilot units and plots:** rivers clipped to the pilot window; `plot_all_in_units` overriding
    `max_per_type`; `SHOW_PLOTS` without a display only warns; `index.html` lists every figure with its
    metrics.
  - **E2, after stage A on T5:**
    - every subsegment of the meandering main river on the correct bank (left = north, looking downstream);
    - core cells counted once across units;
    - 1000 m subsegments aggregated from 500 m equal a direct computation at 1000 m;
    - a second run skips all units through their checkpoints;
    - `recompute_hand` gives HAND 0 on the snapped rivers through the stage A code path.
  - **E3, with a town laid out from the subsegments of E2** (left bank, from 200 m off the thalweg upward, one
    park subsegment without buildings, isolated low mills, a boathouse, synthetic GHSL epochs, zones):
    - the onset of every town subsegment within 0.5 m of the HAND at the town edge, sampled from the rasters;
    - mills of about 6 % of the built area: the onset ignores them, q0.05 and q0.02 follow them;
    - sharpness of the town edge above `sharp_high`;
    - the park and the unbuilt right bank not eligible; the town labelled `dense_urban`;
    - the boathouse excluded through `exclude_building_tags` and logged; `curves.npz` with two GHSL epochs;
    - `curve` and `summary` figures listed in `index.html`.
- **`pixi run test-stage-c`, 26 checks,** on the Stage B country. The Stage B test runs
  first in the same folder if its outputs are missing.
  - **Synthetic inputs:**
    - **GFM:** five scenes at 20 m in EPSG:3035, with flood extents drawn from the
      undefended HAND of the units;
    - **STAC:** a local mock API with token check and GET and POST paging;
    - **GloFAS:** discharge as NetCDF (1981–2021) and GRIB (2019–2021), plus an upstream
      area map;
    - **references:** reference polygons with RPs.
  - **`download_gfm`:**
    - one search;
    - a scene with a few wet pixels screened out;
    - copies in the GFM CRS, the optional likelihood layer included;
    - a rerun screens nothing twice.
  - **Evidence:**
    - held, flooded and `no_test` (all excluded);
    - GFM nodata from the configuration (the file declares none);
    - the date filter;
    - the QGIS GeoPackage.
  - **Hints:**
    - an unmapped dry pocket found;
    - without `hint_exclusion_m`, the levee bands appear as well.
  - **GloFAS:**
    - day attribution;
    - the minimum upstream area;
    - GRIB equal to NetCDF;
    - Gumbel sanity: a fitted sample of 20 000 values gives RP(median) = 2.003, and the
      GloFAS series gives RP(median annual maximum) = 2.00.
  - **Bounds:**
    - matching by upstream area ratio rather than distance;
    - a combined interval;
    - `obs_conflict`;
    - an override kept;
    - an identical rerun;
    - `model_obs_conflict` with rp = 200.
  - **Validation:** recall and precision against an independent computation; the
    confusion matrix and both shares.

`DEFENSES_TEST_DIR` keeps a test folder for inspection.

## 15. Bugs, fixes, limitations, open decisions

### 15.1 Bugs in the original method: reported in Stage A, fixed now (switchable)

**B1. HAND drainage paths ran parallel to the rasterized river.**
- **Effect.** The whole window took the most downstream river elevation as
  `river_ref`, so HAND was too high by up to the river slope × window length.
- **Fix.** Tapered stream burning in the routing copy (`stream_burn_m = 5`,
  `stream_burn_buffer_m = 100`; `stream_burn_m = 0` restores the old behaviour).

**B2. The levee-removal zone overwrote river cells next to levees.**
- **Fix.** River cells keep their elevation, except where they cross a levee
  (`preserve_river_cells`).

**B3. Missing DTM coverage acts as a wall.**
- **Fix.** Compartments closed by it are kept (recall) and flagged `touches_dtm_edge`.

**B4. In gated mode, rivers entering from missing coverage were not treated as entering
from outside.**
- **Fix.** River cells next to nodata now count like river cells on the window edge.

### 15.2 Problems found and fixed in the new code

- **Option B** clipped the copy to the corridor without the halo; it now clips to
  corridor + halo.
- **Overview check** reported an overview *scan* as overview *use*; it now parses the
  source window size.
- **Halo smaller than the corridor buffer** silently produced units without rivers; this
  is now a validation error.
- **`scipy.ndimage` label statistics** scanned whole windows: 18 s → 3 s per unit.
- **Gap filling** first measured gap widths inside a crop without cells outside the gap,
  which capped the river at the lowest cell of the whole window. That is fixed. The rule
  "fill only nodata that does not touch the window edge" also failed for rivers that are
  nodata along their whole length; it was replaced by the coverage closing (§5.5).
- **`rp_index` totals check** first estimated cell coverage by rasterization (up to 3 %
  off on small windows); it now uses exact areas.

**Stage C, found while testing:**
- **Hint connectivity.** The river is reference water, never "wet", so wet cells alone did
  not touch the river lines. Connectivity now runs through wet cells, reference water and
  river cells.
- **Levee bodies as hints.** Levee bodies and mixed pixels turned into hints along every
  levee that held. Fixed with `hint_exclusion_m`.
- **`model_obs_conflict` with one-sided observations.**
  - **Bug.** The open side fell back to the model bound, so a held test above the model
    interval could end with an RP *below* the observation.
  - **Fix.** The open side now stays open while the class is chosen, and then becomes the
    neighbouring class.
- **`download_gfm` request count.** It read every scene once per unit, i.e. scenes × units
  requests. It now reads each scene once for all its units, in threads, with a checkpoint.
- **Matching columns.** The GloFAS matching wrote integers into a float frame row by row;
  it now fills arrays.
- **Dates.** Dates are parsed as ISO 8601 explicitly, which is why `pandas ≥ 2.0` is now
  required.

### 15.3 Limitations

- **L2. Mixed river references.** HAND levels of cells draining to different rivers are
  compared in one flood. Near confluences, a tributary's drainage area can look
  protected against the main river; T2 shows one such area next to a coverage edge. This
  is inherent to HAND. A per-origin reference was tested in an earlier session and
  created false protection downstream of floodplain obstacles.
- **L3. Truncated compartments** not contradicted by a neighbouring window may still be
  artificially closed by window edges. They are flagged; their share is reported.
- **L4. Crests from the DTM.** Floodwalls narrower than the DTM resolution are
  underestimated. The crest search can extend up to one cell beyond the configured
  radius.
- **L5. Platform coverage.** Everything was tested on Linux. On Windows, the spawn pool,
  `peak_rss_mb` (Windows API) and the OSM SQL extraction (needs SpatiaLite in the GDAL
  build) are written but untested.
- **L6. River features are closed as a whole** (§8.2).
- **L7. GHSL resolution.** At 100 m, buildings are spread over the cell, so exposure
  inside small compartments is approximate.
- **L8. What SAR sees.**
  - **Hidden water.** Water under dense vegetation and in built-up areas is hard to detect.
    GFM's exclusion mask removes such cells from the tests, so compartments in towns are
    tested less often.
  - **Wet but not failed.** Water inside a compartment can come from rain, seepage or
    ponding behind a levee that held. `flooded` then overstates failure. The thresholds
    `f_*` are placeholders.
- **L9. Hints use one HAND level per scene** (§10.3). Backwater into tributary valleys can
  create false hints.
- **L10. GloFAS resolution.** At 0.05°, rivers with catchments below a few hundred km² are
  poorly resolved. Their compartments often stay unmatched (`glofas_cell = −1`) and keep
  the model RP.
- **L11. Discharge is not water level.**
  - **Stage.** Overtopping depends on stage. Backwater from confluences, ice jams, changes
    of the channel and hysteresis break the link between discharge RP and stage RP.
  - **Tail.** The Gumbel tail is light (§11.3).
- **L12. Acquisition time.** A held test uses the mean discharge of the acquisition day.
  If the peak arrived *after* the overpass, that mean exceeds the flow the levee had faced
  so far, and the lower bound is too high. The specification accepts this convention.
- **L13. Day attribution of GloFAS values** is configurable but unverified (§11.1).
- **L14. Land cover classes are coarse proxies for asset value.** A WorldCover
  "built-up" pixel can be a road, a yard or a building. Category values are national
  averages per m² and need your calibration (§9.3).

### 15.4 Open decisions (assumptions taken, please confirm)

1. **Download URLs.** The Geofabrik template, GHSL tile templates and GHSL tiling schema
   are **not verified**: this session's network blocked the official hosts and the web
   search found only third-party descriptions. The config entries are empty; the
   scripts stop with a clear message and accept local files.
2. **GHSL epoch.** `0` = the latest epoch not after the current year (e.g. 2025, which is
   modelled), or set 2020 (the last epoch built from observations).
3. **`no_levee` compartments** without an applicable cap get `max_rp`, following "RP from
   the cap only". An alternative is the lowest class.
4. **`rp_source`.** `"cap"` covers both `rp_cap` and `max_rp` lowering the model value.
5. **Maximum damage table format:** columns `country`, `residential`,
   `non_residential`, matched on the country code.
6. **GFM source.** The STAC endpoint, collection and asset keys are **not verified**, and
   the host was not reachable from this environment.
   - **What a search showed.** EODC documentation names `https://stac.eodc.eu/api/v1`,
     collection `GFM` and the assets `ensemble_flood_extent`, `reference_water_mask` and
     `ensemble_likelihood`.
   - **Missing.** No asset name for the exclusion mask was found.
   - The configuration entries are empty; please confirm them.
7. **GFM pixel codes.** Flood 1, exclusion 1, reference water 1 and nodata 255 are
   placeholders.
8. **GloFAS:**
   - the request to paste from the EWDS form;
   - the day attribution (`glofas_valid_time_ends_day`);
   - the URL of the upstream area map;
   - approval of `cdsapi` if downloads should run from the script.
9. **Combination rules** (§11.5). The specification leaves these open; the following
   readings were taken:
   - rp snaps to a class *inside* the combined interval;
   - with `model_obs_conflict`, only the observed interval is used and rp is the class
     inside it closest to the model value; an open side becomes the neighbouring class;
   - overrides are never changed;
   - `obs_conflict` leaves the model values unchanged.
10. **`hint_exclusion_m`:** a new key (default 40 m, two GFM pixels), §10.3.
11. **Flood event for the QGIS check.** Please name one. The steps are:
    1. run `download_gfm.py`;
    2. run `gfm_evidence.py` with `DATE_FROM` and `DATE_TO` around the event;
    3. open `gfm_evidence.gpkg` in QGIS.
12. **Evidence thresholds.** `f_wet_comparison` 0.5, `f_dry_inside` 0.05, `f_wet_inside`
    0.3 and the upstream ratio range 0.7–1.3 are placeholders, as in the specification.
13. **Land cover (HAND_YZ):**
    - the class-to-category mapping and the values and curves per category are yours
      to set;
    - the Urban Atlas class field and codes need checking;
    - the switch `lulc.enabled` is off by default.

### 15.5 Issues in the existing `gis.py` (reported, not changed)

| Function | Issue |
|---|---|
| (module import) | Imports `pyflwdir`, added to `pixi.toml` |
| (module import) | Sets `gdal.SetCacheMax(4 GB)`; workers override it from the config |
| `Raster.update_nodata_vrt` | Calls `FlushCache()` on `None` |
| `Vector.load_ogr_ds_to_geodataframe` | Returns a tuple because of a trailing comma |
| `Raster.raster_to_vect`, `raster_to_mem_vector` | Fail with `filter_values=None` and with a path argument |
| `Files.get_root_path` | Loops forever without an `IF-HydroSim` parent |
| `Vector.extract_intersection` | Creates an exception without `raise` |
| Several functions | Assume EPSG codes; the new functions use WKT |

## 16. Data sources and licences

| Data | Use | Source (to be verified) | Licence |
|---|---|---|---|
| OpenStreetMap | Gates, weirs, pumping stations, drains | Geofabrik country extracts (`.osm.pbf`) | ODbL 1.0; © OpenStreetMap contributors |
| GHSL built-up volume, total and non-residential, R2023A, 100 m | Exposure value | JRC GHSL data package 2023 | to be confirmed from the data package |
| JRC global depth-damage functions | Maximum damage and damage curves | provided by the user as CSV | per the JRC publication |
| Copernicus Emergency Management Service, Global Flood Monitoring (GFM) | held / flooded tests, missing-barrier hints | STAC API, endpoint to be confirmed (§15.4) | Copernicus data policy (free and open); confirm the attribution text on the product page |
| GloFAS historical river discharge, `cems-glofas-historical` | event return periods | CEMS Early Warning Data Store (dataset page checked) | licence shown on the EWDS dataset page, accepted with the download |
| ESA WorldCover 10 m 2021 v200 | land cover exposure (§9.3) | public bucket `esa-worldcover` (eu-central-1), verified | CC BY 4.0; attribution from the file: "ESA WorldCover project 2021 / Contains modified Copernicus Sentinel data (2021) processed by ESA WorldCover consortium" |
| Copernicus Land Monitoring Service, Urban Atlas | land cover exposure where available (§9.3); context zones (§17) | manual download from the CLMS portal | Copernicus data policy; confirm the attribution text with the download |
| Copernicus Land Monitoring Service, CORINE Land Cover (raster) | context zones (§17) | manual download from the CLMS portal | Copernicus data policy; confirm the attribution text with the download |
| GHSL settlement model (SMOD) and population (GHS-POP) | context zones; population in Stage D4 | GHSL data package (templates not verified) | to be confirmed from the data package, with citation |

For GHSL, record the release, epoch and the citation of the GHSL Data Package 2023 with
the outputs: `ghsl_sources.json` lists the epoch, tiles and URLs used.

## 17. Stage D (in progress): industrial sites, dense urban areas, large systems

The approved plan, including the measured runtime impact, is in
[docs/stage_d_plan.md](docs/stage_d_plan.md). Step **D1** is implemented; steps D2–D5
follow.

### 17.1 New inputs (step D1)

- **DSM (`dsm`).** `prepare_dtm.py` inspects and benchmarks every DSM with the code it
  uses for the DTM (`Raster.benchmark_reads`).
  - **Outputs:** `dsm_inspection.csv`, `dsm_benchmark.csv`, `dsm_footprint.gpkg` and
    `work/<c>/dsm/dsm.vrt`.
  - **CRS:** the DSM must share the DTM CRS; the script stops otherwise.
- **OSM (`download_osm.py`).**
  - **Outputs:**
    - `osm_buildings.gpkg`, only without `building_footprints`;
    - `osm_water.gpkg`, only without `water_polygons` (tags `osm_tags_water`);
    - `osm_crossings.gpkg` (`barriers.osm_tags_crossings`);
    - `osm_walls.gpkg`, only when `barriers.osm_tags_walls` is set;
    - `osm_industrial.gpkg` (`context.osm_industrial`).
  - **Areas.** Closed ways tagged `man_made=works` or `waterway=riverbank` count as areas
    (`config/osmconf.ini`).
  - **Choosing the wall tags.** taginfo was not reachable, so the script counts the tags of
    `barriers.osm_wall_candidates` in the country extract and logs the counts; fill
    `osm_tags_walls` from them.
- **GHSL (`download_ghsl.py`).** With templates for `smod` and `pop`, it also mosaics the
  settlement model and population of the same epoch (`ghsl_smod.vrt`, `ghsl_pop.vrt`).
  - **SMOD class codes:** 10 water, 11 very low density rural, 12 low density rural,
    13 rural cluster, 21 suburban or peri-urban, 22 semi-dense urban cluster, 23 dense
    urban cluster, 30 urban centre.
  - **Their status:** confirmed in the Earth Engine catalogue entry of GHS-SMOD P2023A, to be
    confirmed in the GHSL data package.
- **CLMS (`download_clms.py`).** It downloads nothing: the CLMS portal needs a login, and
  no download source could be verified here. It prepares local files:
  - **CORINE:**
    - class index → three-digit CLC code through the legend shipped with the raster
      (columns configurable);
    - verified on known pixels (`corine_check`, required); a mismatch stops the script;
    - clipped to the corridor: `corine_codes.tif`.
  - **Urban Atlas:** the configured layer of multi-layer files, merged and reprojected,
    keeping the polygons that touch the corridor: `urban_atlas.gpkg` with `ua_code`.

### 17.2 Context zones (`context_zones.py`)

1. **Sources per context:**
   - `industrial`: CORINE industrial codes, Urban Atlas industrial codes, OSM industrial
     areas;
   - `dense_urban`: CORINE dense codes, Urban Atlas dense codes, SMOD dense classes;
   - `urban`: SMOD urban classes;
   - `rural`: the rest of the corridor; it must be last in `context.order`.
2. **Order.**
   - **Construction.** Each context is its sources inside the corridor, minus the zones of
     earlier contexts.
   - **Small parts.** Parts below `min_zone_area_m2` are dropped, and their area falls
     through to later contexts. A small industrial patch inside a town becomes
     `dense_urban`, not `rural`.
3. **Output.** `zones.gpkg` (`zone_id`, `context`, `sources`, `area_m2`). The zones do not
   overlap and tile the corridor exactly (tested).
4. **Use in later steps.** `gis.Vector.zones_to_grid` rasterizes the zones onto a unit grid
   for the steps that follow.
5. **Measured.** 7.3 s on a synthetic corridor of 2829 km² with 50 000 Urban Atlas
   polygons, 5000 OSM areas, 2 M CORINE pixels and SMOD at 1 km.

## 18. Stage E (experimental): implied protection in dense built-up reaches

Plan and measured estimates: [docs/stage_e_plan.md](docs/stage_e_plan.md). Outputs stay separate and never enter
`defenses_clip.gpkg` unless `implied.include_implied = true`. All values of `[implied]` are placeholders.

### 18.1 Thalweg snapping (`snap_rivers.py`, step E1)

Snaps the river lines to the thalweg of the DTM. HAND and the bank subsegments of the next steps use the result.

1. **Rivers.**
   - **Source:** `rivers_sources.gpkg` features with `is_source = true` (Stage B), else the rivers input.
   - **Merging:** features joined end to end are merged into one line; confluences and crossings stay breaks.
   - **Pilot units:** with `implied.units`, the lines are clipped to the unit windows, and chainage starts at
     the clipped ends.
2. **Direction.** Mean DTM elevation over the first and the last `direction_check_m`; a line whose start is
   lower than its end is reversed (`flipped`).
3. **Chunks.** `chunk_length_m` with `chunk_overlap_m` overlap. Per chunk:
   - **Read:** `gdal.Warp` average to `resolution_m` over the bounding box of the corridor (chunk buffered by
     `max_offset_m`); cells outside the corridor are impassable.
   - **Narrow river:** cost = z − minimum filter of z over `min_window_m` + `cost_epsilon`.
     - **Nodata** in the corridor is taken as a water surface removed from the DTM and gets the lowest cost.
     - **Corridor mostly nodata** (more than `max_nodata_share`, e.g. a border river at the DTM edge): the
       original chunk is kept. This key is an addition to the specification.
   - **Wide river** (water polygon area in the corridor / chunk length > `wide_river_m`): cost =
     1 / (distance to the water polygon edge + `cost_epsilon`); outside the polygons impassable. The distance
     comes from the polygons rasterized on a window padded by 2 × `max_offset_m` with a dry border, so the
     window edge never acts as a bank.
   - **Path:** `gis.Hydro.least_cost_path` (Dijkstra, numba, 8-connected, diagonal steps × √2) from the cell
     of the chunk start to the cell of its end, or the nearest passable cells.
   - **Quality check:** simplified path length / chunk length at most `max_length_ratio`; otherwise the
     original chunk is kept (`kept_original`) and the reason is logged. Mean and max offset are recorded.
4. **Stitching and smoothing.** Each chunk keeps the part between the middles of its overlaps. The stitched line
   is simplified with tolerance `resolution_m` / 2 and smoothed with two Chaikin iterations; a vertex that would
   move more than `resolution_m` keeps its place.
5. **Chainage** from `chainage_origin` (default: the downstream end) per line, without M values.
6. **Outputs.**
   - `rivers_snapped.gpkg`:
     - specified fields: `river_id`, `flipped`, `n_chunks`, `n_kept_original`, `mean_offset_m`, `max_offset_m`,
       `wide_river_chunks`;
     - additions: `source_river_ids`, `name`, `length_m`, `length_original_m`, `upstream_elev_m`,
       `downstream_elev_m`.
   - `river_points.gpkg`: points every `step_m` with `river_id`, `chainage_m`, and the unit tangent
     `tangent_dx`, `tangent_dy` pointing downstream.
   - `snap_rivers_performance.json`: timings per step and peak RAM.
7. **Plots** (`plots.py`), in `plots/implied/<RUN_ID>/` with `index.html`:
   - **`snap_profile`** per river: DTM elevation along the original and the snapped line against chainage,
     chunk boundaries, kept-original chunks shaded.
   - **`snap_map`** per chunk: overview with the chunk and a detail crop (8 × `max_offset_m`) at the largest
     offset, both on a hillshade with the original line, the least-cost path and the water polygons.
   - **Order of plotting:** flagged items first (kept original, wide, any nodata in the corridor, path at the
     corridor edge, flipped river), then a random sample (`sample_seed`) up to `max_per_type`. Units in
     `plot_all_in_units` get all their items.
   - **`SHOW_PLOTS`** (run variable instead of `--show-plots`): up to `show_max` figures per type are shown
     without blocking; without a display a warning is logged.
   - **Colours** are fixed per entity: water blue, snapped orange, original aqua and dashed (its contrast is
     below 3:1, so the dash and the legend carry it).
8. **Measured** on the benchmark world of Stage A (24 × 24 km, 1 m, 3 rivers, 47 km, 12 chunks, 2 workers):

   | DTM layout | Reads, storage cold (first run) | Reads, OS cache cleared | Cost + path + stitch | Plots (15) | Peak RAM main / worker |
   |---|---|---|---|---|---|
   | tiled, overviews (option A) | 13.9 s (0.30 s/km) | 0.47 s (0.010 s/km) | 0.90 s (0.019 s/km) | 5.2 s | 572 / 438 MB |
   | strip, no overviews | 75.4 s (1.6 s/km) | 5.4 s (0.11 s/km) | 0.79 s | 7.9 s | 1025 / 764 MB |

   - **The two read columns** show how much storage matters. The first run read files idle since the day
     before; the repeat dropped only the page cache of the machine (`posix_fadvise`), while the storage below
     stayed warm.
   - **Per 1000 km of river:** about 20 s of computation plus 10–300 s of reads on a tiled DTM, or 110–1600 s
     on a strip DTM, divided by the number of workers. Plots add 0.35–0.5 s per figure; the main process
     reads the plot windows with the configured GDAL cache (`gdal_cachemax_mb`).

### 18.2 Bank subsegments (`implied_segments.py`, step E2)

Per unit in the process pool, with a checkpoint per unit (`subsegments_done.json`, hash of the parameters and
input times; `FORCE = True` recomputes).

1. **HAND.** `hand_undef.tif` of stage A. With `recompute_hand = true`, HAND is recomputed with the snapped
   rivers as streams on the undefended DTM of stage A (`dtm_undefended.tif`, kept with `keep_intermediate =
   true`). It uses `gis.Hydro.hand_above_streams`, the code path that `compartments.py` now calls as well, and
   writes `hand_undef_snapped.tif`.
2. **Envelope.** Cells with HAND ≤ `envelope_hand_m`, 4-connected to the cells of the snapped rivers.
3. **Allocation.** River points are rasterized onto the unit grid. An EDT with indices over the whole window
   (core and halo) gives every envelope cell its nearest point: `river_id`, chainage and tangent. Side = sign of
   the cross product of the downstream tangent and the vector to the cell (left or right, looking downstream).
4. **Subsegments.** `k = floor(chainage / subsegment_length_m)`; labels encode (`river_id`, side, `k`) with a
   lookup table. Only core cells get labels, so every cell is counted once; the halo only serves the
   nearest-point search.
5. **Outputs.**
   - **Per unit:** `subsegments.tif` (labels, 0 = none) and `subsegments.csv` (`label`, `river_id`, `side`, `k`,
     `n_cells`, `chainage_from_m`, `chainage_to_m`, `envelope_area_m2`).
   - **Merged:** `subsegments.gpkg`, dissolved across units by (`river_id`, `side`, `k`), with
     `subsegment_id`.
6. **Longer subsegments** of `subsegment_length_grid` are aggregated from the base length (`k // factor`), never
   recomputed. The test checks that the aggregation equals a direct computation.
7. **Measured** on the two pilot units of the benchmark world (30.25 M cells each, 2 m):
   - **Per unit:** about 10 s: HAND read 6 s, envelope 0.5 s, allocation (EDT with indices) 3 s, writing 0.9 s.
   - **Memory:** peak RAM 1.1 GB per worker.
   - **Output:** 80 subsegments of 500 m, envelope 34.4 km².

### 18.3 Built-up curves and levels (`implied_protection.py`, parts A–B, step E3)

Results go to `implied/<RUN_ID>/` (with `run_config.toml`, a copy of the configuration used), plots to
`plots/implied/<RUN_ID>/`.

1. **Built-up fraction per unit cell** (pass 1, per unit in the process pool), saved as `implied_built.tif`:
   - **footprints:**
     - sources: the national footprints (`building_footprints`) and the OSM buildings (`osm_buildings.gpkg`,
       now always extracted by `download_osm.py`, with `building` and `man_made`);
     - features matching `exclude_building_tags` are dropped, with counts logged per tag;
     - rasterized at a quarter of the working cell without all-touched, in tiles of 4096² fine cells, and
       block-averaged.
   - **ghsl:** the GHSL built-up surface (m² per cell) averaged onto the grid as a share of its cell. The file is
     `ghsl_built_s_10m.vrt` if present, else the latest of the epochs.
   - **max:** the larger of the two per cell. Industrial land use is never treated as built-up; `zones.gpkg`
     only labels the context.
2. **Curves.** Bins from 0 to `envelope_hand_m` by `hand_bin_m`. Envelope and built area per (subsegment, bin)
   with `np.bincount`, per unit, summed over units in the main process. Built area is also counted per GHSL
   epoch of `ghsl_epochs`, for the temporal signal of E4. All of it goes to `curves.npz`.
3. **Levels per subsegment** (`gis.Stats`):
   - **quantile levels:** area-weighted HAND quantiles of the built area, interpolated within the bin;
   - **onset:**
     - share = built / envelope per bin, NaN for bins below `min_bin_area_m2`;
     - smoothed with a centred moving average over `smooth_bins`, ignoring NaN;
     - onset = the lower edge of the lowest bin where the smoothed share is at least `onset_fraction` × its
       maximum and stays so for `onset_persist_bins` bins;
   - **sharpness:** (mean smoothed share in [onset, onset + `sharp_window_m`] − mean below) / maximum.
4. **Eligibility.** `built_share_total` ≥ `built_share_min` and built area ≥ `built_area_min_m2`, and an onset
   exists; the reason is stored when not eligible.
5. **Context.** Majority zone by area (`zones.gpkg`), else `none`.
6. **Output.** `implied_subsegments.gpkg`: the subsegments of E2 with `built_area_m2`, `built_share_total`,
   `eligible`, `reason`, `sharpness`, `level_onset_m`, `level_q0_02_m`, … , `level_m` (primary method),
   `context`.
7. **Plots.**
   - **`curve`** per subsegment. Upper panel: envelope area per bin (bars), built share and smoothed share on the
     secondary axis (as specified), lines at every level. Lower panel: cumulative built fraction against HAND.
   - **`summary`:** histograms of levels, sharpness and built share, stacked by context; subsegments per context.
   - **Flagged first:** eligible with low sharpness, eligible with onset and lowest quantile more than
     `sharp_window_m` apart, onset found but not eligible.
   - **Colours:** each level method keeps its colour and dash pattern in every figure.
8. **Measured** on the two pilot units of the benchmark world (30.25 M cells each) with 100 000 synthetic
   buildings per core: about 5 s per unit (reads 1.1 s, rasterization of the footprints at 0.5 m 3.7 s, curves
   0.3 s). Peak RAM 1.05 GB per worker; figures 0.28 s each.

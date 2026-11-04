# flood_defenses — Defenses layer for the fluvial flood hazard pipeline

The package derives **compartments**: areas protected by flood levees. Each compartment
gets a return period (RP) of protection. The layer clips fluvial hazard maps in
post-processing, so it is built **only from terrain, levee lines, the river network and
free open data** (OpenStreetMap, GHSL; GFM and GloFAS in Stage C). No hazard map, hydraulic
model output, local hazard dataset or FABDEM is used anywhere.

**Error priority: recall.** When the evidence is ambiguous, a compartment is kept rather
than dropped. Every design choice below that trades precision for recall says so
explicitly.

| Stage | Content | Status |
|---|---|---|
| A | `gis.py` additions, `prepare_dtm.py`, `build_units.py`, `compartments.py` | done, tested, benchmarked |
| A+ | faster kernels; fixes B1–B4; gap filling for water surfaces | done, tested, benchmarked |
| B | `download_osm.py`, `osm_hydraulics.py`, `download_ghsl.py`, `rp_index.py` | done, tested on synthetic data; download URLs not verified |
| C | `download_gfm.py`, `gfm_evidence.py`, `download_glofas.py`, `rp_bounds.py`, `validate.py` | not started |

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
10. [Measured performance](#10-measured-performance)
11. [Tests](#11-tests)
12. [Bugs, fixes, limitations, open decisions](#12-bugs-fixes-limitations-open-decisions)
13. [Data sources and licences](#13-data-sources-and-licences)

---

## 1. Running the scripts

**No command-line parameters.** Each script starts with a block of *run variables*:

```python
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"
FORCE = False   # recompute finished units / footprint / conversion / extraction
UNITS = []      # compartments.py, rp_index.py: selected unit ids, empty = all
```

Everything else comes from the TOML configuration: data paths, thresholds, resolutions,
URLs and tag lists. Code holds no paths, CRS, URLs or thresholds.

Batch runs and the tests override the run variables through environment variables named
`DEFENSES_<VARIABLE>`, for example `DEFENSES_COUNTRY=CZ` or
`DEFENSES_UNITS=012_034,012_035`. Without them the values in the script apply.

**Scripts run as modules from the repository root.** `python -m package.script` puts the
working directory on `sys.path`, so `import gis` works without `sys.path` manipulation or
packaging:

```
pixi run prepare-dtm        pixi run download-osm       pixi run test-kernels
pixi run build-units        pixi run osm-hydraulics     pixi run test-pipeline
pixi run compartments       pixi run download-ghsl      pixi run test-stage-b
                            pixi run rp-index
```

`pixi run python -m processing.compartments` works the same way. Running a script by its
file path does **not** work, because Python then puts the script folder on the path
instead of the root.

**Dependencies** are listed in `pixi.toml`:
- `python ≥ 3.11` (for `tomllib` and `ProcessPoolExecutor(max_tasks_per_child)`),
- `gdal ≥ 3.11` (`MEM` vector datasets; the OSM driver; SpatiaLite SQL for the OSM
  extraction),
- `numpy`, `scipy`, `numba`, `geopandas ≥ 1.0`, `shapely ≥ 2.0`, `pyproj`,
- `pyflwdir`, only because the existing `gis.py` imports it (it depends on affine,
  numba, numpy and scipy, not on rasterio).

Downloads use the standard library (`urllib`). No rasterio-based library is used.

## 2. Repository layout

```
gis.py                     shared GIS library (existing; new functions appended, nothing changed)
config/defenses.toml       configuration template
config/osmconf.ini         OSM driver configuration exposing the tag keys used
download/
  download_osm.py          Geofabrik extract -> hydraulic points, drains (single-layer GeoPackages)
  download_ghsl.py         GHSL building volume tiles -> VRT mosaics
processing/
  prepare_dtm.py           DTM inspection, footprint, benchmark, estimates, optional conversion
  build_units.py           river corridor, processing units, memory estimate
  osm_hydraulics.py        water sources (closed streams), pumping stations
  compartments.py          per-unit compartments in parallel, checkpoints, merging
  rp_index.py              exposure, score, RP per compartment
tests/
  test_kernels.py          numba kernels on small grids (16 checks)
  test_pipeline.py         end-to-end synthetic countries for Stage A (40 checks)
  test_stage_b.py          end-to-end synthetic country for Stage B (20 checks)
  make_benchmark_world.py  24 x 24 km synthetic DTM for performance runs
validate_defenses.py       existing validation script (moves to processing/validate.py in Stage C)
download_nld.py            US National Levee Database download (optional)
```

### New `gis.py` functions

Functions are appended to the existing classes or added as new classes; no existing
function was modified. The only change to existing lines is the import block.

| Class | Function | Purpose |
|---|---|---|
| `Raster` | `inspect`, `build_vrt` | Storage properties of a raster; VRT over several DTMs (stops on a CRS mismatch) |
| `Raster` | `read_window(src, bounds, res, resample, warp_mb, overview, multithread, dst_wkt)` | `gdal.Warp` of a map window to MEM, optionally reprojected |
| `Raster` | `read_native(src, bounds)` | Direct read of native pixels |
| `Raster` | `array_to_mem`, `write_array`, `best_compression` | Array ↔ GDAL; tiled compressed GeoTIFF; ZSTD → DEFLATE fallback |
| `Raster` | `rasterize_geometries`, `polygonize_labels`, `fill_nodata_cells`, `footprint` | Burn geometries; label raster → polygons; `gdal.FillNodata` on selected cells; valid-data footprint |
| `Vector` | `check_single_layer`, `read_merged`, `write_gpkg`, `length_near_boundary` | Single-layer inputs; merged reads with `source`; single-layer GeoPackage output; levee length near compartment boundaries |
| `Hydro` | `breach_to_streams(dem, stream, burn, burn_buffer_cells)` | River cell reached downstream (numba) |
| `Hydro` | `priority_flood(surface, seeds)` | Spill level, seed index and reach order per cell (numba) |
| `Hydro` | `distance_cells(mask, limit)` | Chamfer distance in cells (numba) |
| `Config` | `load`, `get`, `as_list`, `warn_placeholders` | TOML reading and validation |
| `Run` | `env_overrides`, `logger`, `peak_rss_mb`, `configure_gdal`, `run_units`, `format_seconds` | Run variables, run logs, memory measurement (Linux/macOS/Windows), process pool |
| `Http` | `download(url, path, retries, pause_s, headers)` | urllib download with size check, skip and backoff; `file://` supported |
| `Stats` | `percentile_rank`, `quantile_classes` | Ranks and quantile mapping for the RP score |

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
| `compartments.stream_burn_m`, `stream_burn_buffer_m` | B1 fix (§12.1) |
| `compartments.preserve_river_cells` | B2 fix |
| `compartments.cap_river_holes`, `max_fill_distance_m` | Gap filling for water surfaces (§5.5) |
| `rp.value_total_tolerance` | Tolerance of the value-warp totals check (spec: 1 %) |
| `sources.request_pause_s`, `max_retries` | Polite downloads with backoff |
| `sources.osm_tags_closing`, `osm_tags_pumping` | Which hydraulic structures close a river or pump water |
| `sources.ghsl_url_templates` (table `total`, `non_residential`), `ghsl_tile_schema`, `ghsl_tile_field`, `ghsl_nodata_as_zero` | GHSL tiles; a table instead of the list in the original schema, so that products cannot be swapped |

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
| 8 | `rp_index.py` | per unit `rp_index.csv`, `rp_index.json`; merged `compartments_rp.gpkg` |

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
     run along the channel and reach the river only at its lowest point (§12.1).
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
per 30 M-cell unit (§10): terrain read **8–9 s** in strip layout, **1.7–2.3 s** after
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
     specification says ("from the cap only"). See the open decisions (§12.4).
7. **Interval:**
   - `rp_low` = the previous value in `classes` (the same value if `rp` is the lowest);
   - `rp_high` = the next value in `classes`, bounded by `rp_cap` and `max_rp`, never
     below `rp`.
8. **Overrides.** A compartment containing a point of the overrides file gets its `rp`,
   with `rp_low = rp_high = rp` and `rp_source = "override"`. If points disagree, the
   first is used and a warning is logged.
9. **Stage C columns.** `obs_rp_low`, `obs_rp_high`, `n_held`, `n_flooded`,
   `obs_conflict` and `model_obs_conflict` are present and empty, to be filled by
   `rp_bounds.py`. `value_exposed` is an added column.

## 10. Measured performance

### 10.1 Where the numbers come from

All numbers are from this cloud container: 4 vCPUs, 15 GB RAM, virtual disk. The OS page
cache was dropped before timed runs. The synthetic benchmark DTM comes from
`tests/make_benchmark_world.py`:
- two GeoTIFFs of 12 000 × 24 000 px at 1 m;
- strip layout, no overviews.

Working resolution 2 m, cores 7 km, halo 2 km. That gives windows of
5500² = 30.25 M cells, close to a real unit of 20 + 2 × 3 km at 5 m (27 M).

**Disk and decompression speeds on your machines will differ.**

### 10.2 DTM reads (20 windows per file)

| Layout | Native 4096² window | `gdal.Warp` to 2 m (same ground) | Sequential | Overview used |
|---|---|---|---|---|
| Strip, no overviews | 0.31–0.45 s (148–214 MB/s) | 0.74–0.81 s (82–90 MB/s native eq.) | 530–640 MB/s | no |
| Option A (tiled, ZSTD, overviews) | 0.14–0.20 s (342–468 MB/s) | 0.34–0.38 s (178–197 MB/s) | 368–501 MB/s | yes (factor 2) |

Conversion results:
- **Size.** ZSTD with predictor 3 compressed the synthetic DTM to 0.47 of its raw size.
- **Time.** Option A took 8–11 s per 1.15 GB file plus 5–6 s for overviews.

### 10.3 Kernels on one unit (30.25 M cells)

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

### 10.4 `compartments.py` per unit (30.25 M window cells, 24.4 M domain cells, cold cache, option A source)

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

### 10.5 `rp_index.py` per unit

With synthetic GHSL over the benchmark world it takes **1.6–4.6 s** with 0.6–0.9 GB peak.
Warped and source volume agree within 0.05 %.

## 11. Tests

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

`DEFENSES_TEST_DIR` keeps a test folder for inspection.

## 12. Bugs, fixes, limitations, open decisions

### 12.1 Bugs in the original method: reported in Stage A, fixed now (switchable)

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

### 12.2 Problems found and fixed in the new code

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

### 12.3 Limitations

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

### 12.4 Open decisions (assumptions taken, please confirm)

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

### 12.5 Issues in the existing `gis.py` (reported, not changed)

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

## 13. Data sources and licences

| Data | Use | Source (to be verified) | Licence |
|---|---|---|---|
| OpenStreetMap | Gates, weirs, pumping stations, drains | Geofabrik country extracts (`.osm.pbf`) | ODbL 1.0; © OpenStreetMap contributors |
| GHSL built-up volume, total and non-residential, R2023A, 100 m | Exposure value | JRC GHSL data package 2023 | to be confirmed from the data package |
| JRC global depth-damage functions | Maximum damage and damage curves | provided by the user as CSV | per the JRC publication |
| CEMS GFM, GloFAS | Stage C | — | — |

For GHSL, record the release, epoch and the citation of the GHSL Data Package 2023 with
the outputs: `ghsl_sources.json` lists the epoch, tiles and URLs used.

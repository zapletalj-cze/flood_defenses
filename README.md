# flood_defenses — Defenses layer for the fluvial flood hazard pipeline

The package derives **compartments**: areas protected by flood levees. Each compartment
will later get a return period (RP) of protection. The layer clips fluvial hazard maps in
post-processing, so it is built **only from terrain, levee lines, the river network and
free open data**. No hazard map, hydraulic model output or FABDEM is used anywhere.

**Error priority: recall.** When the evidence is ambiguous, a compartment is kept rather
than dropped. Every design choice below that trades precision for recall says so
explicitly.

Status: **Stage A** is implemented and tested:
- `gis.py` additions,
- `prepare_dtm.py`,
- `build_units.py`,
- the `compartments.py` refactor.

Stages B (OSM, GHSL, RP index) and C (GFM, GloFAS, RP bounds, validation) follow after
review; see [Implementation stages](#11-implementation-stages).

---

## Contents

1. [Running the scripts](#1-running-the-scripts)
2. [Repository layout](#2-repository-layout)
3. [Configuration](#3-configuration)
4. [Run order and outputs](#4-run-order-and-outputs)
5. [Method: from levee lines to protected compartments](#5-method-from-levee-lines-to-protected-compartments)
6. [Large DTMs: corridors, storage layout, units](#6-large-dtms-corridors-storage-layout-units)
7. [Parallel processing, checkpoints, merging](#7-parallel-processing-checkpoints-merging)
8. [Measured performance](#8-measured-performance)
9. [Tests](#9-tests)
10. [Bugs found, limitations, open decisions](#10-bugs-found-limitations-open-decisions)
11. [Implementation stages](#11-implementation-stages)
12. [Data sources and licences](#12-data-sources-and-licences)

---

## 1. Running the scripts

**No command-line parameters.** Each script starts with a block of *run variables*:

```python
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"
FORCE = False   # recompute finished units / footprint / conversion
UNITS = []      # compartments.py: selected unit ids, empty = all
```

Everything else comes from the TOML configuration: data paths, thresholds, resolutions,
URLs and CRS-independent settings. Code holds no paths, CRS, URLs or thresholds.

Batch runs and the tests override the run variables through environment variables named
`DEFENSES_<VARIABLE>`, for example `DEFENSES_COUNTRY=CZ` or
`DEFENSES_UNITS=012_034,012_035`. Without them the values in the script apply.

**Scripts run as modules from the repository root.** `python -m package.script` puts the
working directory on `sys.path`, so `import gis` works without `sys.path` manipulation or
packaging:

```
pixi run prepare-dtm      # = python -m processing.prepare_dtm
pixi run build-units      # = python -m processing.build_units
pixi run compartments     # = python -m processing.compartments
pixi run test-kernels
pixi run test-pipeline
```

`pixi run python -m processing.compartments` works the same way. Running a script by its
file path (`python processing/compartments.py`) does **not** work, because Python then
puts `processing/` on the path instead of the root.

**Dependencies** are listed in `pixi.toml`:
- `python ≥ 3.11` (for `tomllib` and `ProcessPoolExecutor(max_tasks_per_child)`),
- `gdal ≥ 3.11` (in-memory vector datasets through the `MEM` driver),
- `numpy`, `scipy`, `numba`, `geopandas ≥ 1.0`, `shapely ≥ 2.0`, `pyproj`,
- `pyflwdir`, only because the existing `gis.py` imports it. It depends on affine, numba,
  numpy and scipy and does not pull in rasterio.

No rasterio-based library is used.

## 2. Repository layout

```
gis.py                     shared GIS library (existing; new functions appended, nothing changed)
config/defenses.toml       configuration template
processing/
  prepare_dtm.py           DTM inspection, footprint, benchmark, estimates, optional conversion
  build_units.py           river corridor, processing units, memory estimate
  compartments.py          per-unit compartments in parallel, checkpoints, merging
download/                  Stage B/C download scripts (not yet written)
tests/
  test_kernels.py          numba kernels on small grids (16 checks)
  test_pipeline.py         end-to-end synthetic country (28 checks)
  make_benchmark_world.py  24 x 24 km synthetic DTM for performance runs
validate_defenses.py       existing validation script (moves to processing/validate.py in Stage C)
download_nld.py            US National Levee Database download (optional)
```

### New `gis.py` functions

Functions are appended to the existing classes; no existing function was modified. The
only change to existing lines is the import block.

| Class | Function | Purpose |
|---|---|---|
| `Raster` | `inspect(path)` | Size, pixel size, type, nodata, compression, block size, overviews, file and raw size |
| `Raster` | `build_vrt(paths, out)` | VRT over several DTMs; stops on a CRS mismatch |
| `Raster` | `read_window(src, bounds, res, resample, warp_mb, overview, multithread)` | `gdal.Warp` of a map window to MEM on a grid aligned with the bounds |
| `Raster` | `read_native(src, bounds)` | Direct read of native pixels |
| `Raster` | `array_to_mem`, `write_array`, `best_compression` | Array ↔ GDAL; tiled compressed GeoTIFF with the right predictor; ZSTD → DEFLATE fallback |
| `Raster` | `rasterize_geometries(geoms, values, shape, gt, wkt, all_touched, dtype)` | Burns shapely geometries onto a grid |
| `Raster` | `polygonize_labels(labels, gt, wkt)` | One dissolved polygon per label (`gdal.Polygonize`, 4-connected) |
| `Raster` | `fill_nodata_cells(array, to_fill, gt, wkt, max_cells)` | `gdal.FillNodata` restricted to the bounding box of the cells to fill |
| `Raster` | `footprint(path, res)` | Valid-data polygon from a decimated read |
| `Vector` | `check_single_layer`, `read_merged`, `write_gpkg`, `length_near_boundary` | Single-layer inputs, merged reads with `source`, single-layer GeoPackage output, levee length near compartment boundaries |
| `Hydro` | `breach_to_streams(dem, stream, burn, burn_buffer_cells)` | Stream cell reached downstream (numba) |
| `Hydro` | `priority_flood(surface, seeds)` | Spill level, seed index and reach order per cell (numba) |
| `Config` | `load`, `get`, `as_list`, `warn_placeholders` | TOML reading and validation |
| `Run` | `env_overrides`, `logger`, `peak_rss_mb`, `configure_gdal`, `run_units`, `format_seconds` | Run variables, run logs, memory measurement (Linux/macOS/Windows), process pool |

## 3. Configuration

`config/defenses.toml` is the template:
- `[general]`: directories, workers, memory limits, GDAL cache and threads.
- `[benchmark]`: benchmark windows.
- `[countries.<code>]`: inputs, resolutions, units, conversion, with sub-tables
  `compartments` and `rp`.
- `[sources]` and `[evidence]`: Stages B/C.

All path keys accept a string or a list. Every script validates exactly the keys it uses,
before doing any work:
- **presence and type;**
- **file existence;**
- **allowed values** (`dtm_conversion` A/B, `dtm_compress`, `log_level`);
- **consistency between keys**: share sums, list lengths, ascending classes, and
  `unit_halo_m ≥ corridor_buffer_m` (reason in §6.4).

All problems are listed in one message and the script stops.

**Placeholders.** The top-level array `placeholders` lists keys whose values are not
confirmed yet. Every script logs a warning for each placeholder it uses. Remove a key
from the list once its value is confirmed.

New keys beyond the original schema, with the reason for each:

| Key | Reason |
|---|---|
| `general.gdal_threads_per_worker` | Avoid oversubscription: n workers × GDAL threads ≤ cores |
| `benchmark.*` | Window count, size and seed of the benchmark (reproducible positions) |
| `country.dtm_conversion` | `"A"` whole file tiled, `"B"` corridor copy (§6.3) |
| `country.dtm_compress`, `footprint_resolution_m` | Conversion compression; footprint cell size |
| `country.bytes_per_cell` | Memory model for unit sizing; **replace with the value logged by `compartments.py`** |
| `country.crest_box_size_m` | Size of the boxes read at native resolution around levees |
| `country.merge_overlap_ratio`, `drop_contradicted_truncated`, `contradiction_ratio` | Merging of units (§7.3) |
| `country.id_rounding_m` | Rounding of the representative point in `compartment_id` |
| `country.keep_intermediate` | Write DTM, spill, margin, barrier and crest rasters per unit |
| `compartments.touch_tolerance_m`, `self_close_min_length_factor` | Former code constants of the bridging |
| `compartments.gated_crossings`, `fill_nodata_holes` | Former command-line switches |
| `compartments.stream_burn_m`, `stream_burn_buffer_m` | Fix for bug B1 (§10.1), **0 = previous behaviour** |

## 4. Run order and outputs

| Step | Script | Writes (`output_dir/<country>/` unless noted) |
|---|---|---|
| 1 | `prepare_dtm.py` | `dtm_inspection.csv`, `dtm_footprint.gpkg`, `dtm_benchmark.csv`, `dtm_estimates.json`, `work/<c>/dtm/dtm.vrt` |
| 2 | `build_units.py` | `corridor.gpkg`, `units_core.gpkg`, `units_window.gpkg`, `units.csv`, `work/<c>/units_domain.gpkg` |
| 3 (optional) | `prepare_dtm.py` with `convert_dtm = true` | tiled copies in `work/<c>/dtm/`, `dtm.vrt` repointed to them |
| 4 | `compartments.py`, first `UNITS = [two units]` | per unit in `work/<c>/units/<id>/`: `hand_undef.tif`, `labels.tif`, `compartments.gpkg`, `barrier_crest.csv`, `done.json`, `unit.log`; merged: `compartments.gpkg`, `barriers.gpkg`, `enclosed_streams.gpkg`, `performance.json` |
| 5 | `compartments.py`, `UNITS = []` | remaining units; finished units are skipped |

Run logs go to `output_dir/<country>/logs/<script>_<timestamp>.log`.

**Notes on the steps:**
- **Order of steps 1–3.** Step 1 runs before `build_units.py` because the corridor needs
  the footprint. Option B needs the corridor, so conversion (step 3) comes after step 2.
  Rerunning `prepare_dtm.py` after step 2 also adds the per-country read estimate and the
  option B disk estimate.
- **`compartments.gpkg` columns:** `compartment_id`, `area_m2`, `protection_level_m`,
  `margin_median_m`, `margin_max_m`, `river_ref_median_m`, `entry_x`, `entry_y`,
  `source_x`, `source_y`, `truncated`, `touches_corridor_edge`, `has_pumping`,
  `no_levee`, `unit_id`, `unit_label`, geometry.
- **`unit_label`** (an addition) points to the compartment's cells in the unit's
  `labels.tif`; Stage B/C scripts need it to aggregate rasters per compartment.

## 5. Method: from levee lines to protected compartments

### 5.1 Principle

A levee protects an area if water from the river has to rise **higher** to reach that
area with the levee than without it. Both cases are computed on the same
height-above-river surface:

* `hand_undef` (Height Above Nearest Drainage on the *undefended* DTM): elevation of a cell
  minus the elevation of the river cell its drainage path reaches (`river_ref`).
* **Spill level** of a cell: the lowest level above the river at which water starting at
  any river cell reaches the cell, moving between 4-neighbours. It is computed by a
  priority flood (Dijkstra with the max-operator):
  `spill(n) = min over paths of max(surface along the path)`, starting from 0 at river
  cells. The heap processes cells in increasing spill order, so each cell is final when
  it is first taken from the heap.
* `spill_undef`: the flood over `hand_undef`. `spill_def`: the same flood over the
  *defended surface*, which equals `hand_undef` except on barrier cells, where it is
  `max(crest − river_ref, hand_undef)`.
* `margin = spill_def − spill_undef ≥ 0`. A cell is **protected** if
  `margin ≥ min_margin_m` and `spill_undef ≤ max_level_m`. The second condition excludes
  terrain so high that the river does not reach it in any realistic event.
* 4-connected protected regions form compartments; regions smaller than `min_area_m2` are
  dropped.
* **`protection_level_m`** is the minimum of `spill_def` in the compartment: the level
  above the river at which water first enters.

Consequences that look odd at first but are correct:
- **Upslope parts are not protected.** Parts of a ring that lie higher than its lowest
  crest are not protected, because water would reach them at their own height anyway.
- **Nested rings.** The inner compartment of nested rings has the higher of the two
  crests as its protection level. The test `nested rings` checks both orderings.

### 5.2 Barriers

1. **Merging.** All levee files are read, reprojected to the DTM CRS and split into
   single lines. The column `source` holds the file stem. Bridging and line ids are
   computed **once per country** (`work/<c>/barrier_lines.gpkg`), not per unit, so all
   units see identical barriers and seams cannot disagree.
2. **Bridging** (vectorised with a spatial index):
   - **Rule.** Each line endpoint with another line within `gap_tolerance_m` is joined
     to the nearest such line: to that line's endpoint if one is within the tolerance,
     otherwise to the nearest point on the line (T-junctions).
   - **Already connected.** An endpoint closer than `touch_tolerance_m` to another line
     counts as connected and gets no bridge.
   - **Self-closing.** A single line whose own ends are within the tolerance closes on
     itself (an almost-closed ring digitised as one line) if it is longer than
     `self_close_min_length_factor × gap_tolerance_m`. The length condition stops short
     lines from closing.
   - Bridges have `bridged = True` and `source = "a|b"`.
3. **Rasterization** on the working grid uses `ALL_TOUCHED`. Diagonal-only contacts are
   then closed by adding the higher of the two corner cells, so that barrier cells form
   4-connected lines. A 4-connected flood cannot pass between two barrier cells that
   touch only at a corner, so this is extra safety. It also guarantees that a river
   crossing a staircase levee cuts a barrier cell.
4. **Crest elevation** is the maximum of the *native* DTM near the line, because line
   geometry can be offset from the real crest:
   - **Where native pixels are read.** Only in boxes of `crest_box_size_m` around the
     levees (buffer = search radius + one cell). Boxes in one row that are separated by
     at most one empty box are read together.
   - **How each box is read.** `gdal.Warp` with `resampleAlg = "max"` gives the maximum
     native elevation per working cell. **Overviews are disabled for this read**
     (`overviewLevel = "NONE"`): overviews built with `AVERAGE` would lower the crest.
   - **Search radius.** A circular maximum filter with `crest_search_radius_m` follows.
     The effective search radius is up to one working cell larger than configured, which
     raises crests slightly; this favours recall.
5. **Bridged cells** get crests interpolated linearly along the bridge between the crests
   at its two ends. If one end is unknown, the other end's crest is used.
6. **Rivers crossing barriers.** By default the river wins: crossing cells are not
   barrier. With `gated_crossings = true` crossings stay closed and river cells in areas
   enclosed by barriers are not water sources.
   - Such areas are listed in `enclosed_streams.gpkg` in both modes.
   - In Stage B, `rivers_sources.gpkg` (from OSM sluices, gates and pumping stations)
     replaces this blanket switch with per-feature decisions.

### 5.3 Undefended DTM

All barrier cells, including 4-connectivity additions and bridges, are dilated by
`barrier_removal_width_m`. The result is set to nodata and interpolated with
`gdal.FillNodata`. The dilation uses a Euclidean distance transform restricted to the
bounding box of the barriers, to save memory. This is the terrain as it would be without
levees.

### 5.4 HAND on the undefended DTM

1. **Stream cells** are the river lines rasterized along their centre line (not
   all-touched), so that rivers cut as few barrier cells as possible.
2. **Flow routing** runs on a *conditioned copy* of the undefended DTM:
   - **Tapered stream burning** (§10.1, bug B1): river cells are lowered by
     `stream_burn_m`, decreasing linearly to 0 at `stream_burn_buffer_m`. This gives
     channel cells a gradient towards the rasterized river.
   - **Complete breaching:** a priority flood from the river cells. When a lower cell (a
     pit) is reached from a higher one, the path back towards the river is carved below
     the pit. Flats get an ε gradient (1e-5 m) towards the outlet. Breaching is preferred
     to filling because it keeps the conditioned surface close to the real one.
   - **D8 steepest descent** on the conditioned surface. Every non-river cell has a
     strictly lower neighbour, so processing cells by increasing elevation resolves every
     path in a single pass.
3. **`river_ref`** is the elevation of the river cell reached, taken from the
   *unconditioned* undefended DTM. `hand_undef = dtm_undefended − river_ref`.
4. **Fallback.** Cells not connected to any river cell take the nearest river cell; their
   number is logged.
5. **Only modelled rivers are outlets.** Side catchments are therefore routed to a
   modelled river. This is intended: every cell gets a reference of a river that the
   hazard maps contain.

### 5.5 Seeds, entry and source

The priority flood also returns, for every cell:
- **the seed** (river cell) where its minimal spill path starts,
- **the order** in which cells were reached.

For each compartment:
- **`entry_x`, `entry_y`:** the compartment cell reached first by the defended flood, at
  its cell centre. This is where water first enters.
- **`source_x`, `source_y`:** the seed of the entry cell, i.e. the river cell whose
  water enters first. Stage C matches it with GloFAS upstream areas.

Taking the cell with the minimal `spill_def` would be ambiguous: inside a compartment,
every cell below the crest has the same spill level. The reach order resolves this.

### 5.6 Flags

| Flag | Meaning |
|---|---|
| `truncated` | The compartment touches the edge of the unit window, where it may continue (§7.3) |
| `touches_corridor_edge` | It is adjacent to cells outside the processing domain (corridor + halo) |
| `has_pumping` | A pumping station lies inside or within `levee_adjacent_distance_m` of it (needs `osm_pumping.gpkg`, Stage B; otherwise False) |
| `no_levee` | No barrier line within `levee_adjacent_distance_m` of its boundary |

## 6. Large DTMs: corridors, storage layout, units

National DTMs of 150–600 GB are never loaded whole. All reads go through `gis.Raster`
window functions on processing units.

### 6.1 Footprint

`Raster.footprint` reads the DTM decimated to `footprint_resolution_m` in a single
`ReadAsArray(buf_xsize, buf_ysize)` call:
- **with overviews**, GDAL serves the read from the closest overview;
- **without overviews**, it decodes only the strips or tiles that contain the sampled rows
  and columns.

The valid mask is dilated by one coarse cell, so that narrow data strips lost by
decimation are kept (recall), and then polygonized. The footprint is cached;
`FORCE = True` recomputes it. On the 24 × 24 km benchmark DTM in strip layout it took
0.1 s per file.

### 6.2 Corridor

Corridor = rivers buffered by `corridor_buffer_m`, intersected with the footprint.
`build_units.py` reports its share of the valid DTM area; it was 30 % on the benchmark
world. Compartments far from the modelled rivers cannot be flooded by them, so they are
irrelevant for clipping fluvial hazard maps.

### 6.3 Storage layout and conversion

`prepare_dtm.py` reports for every DTM:
- size, pixel size, data type, nodata, compression, block size, overviews, file size;
- a **benchmark**:
  - 20 random 4096 × 4096 native windows inside the footprint, read with `ReadAsArray`;
  - the same ground windows read with `gdal.Warp` to the working resolution;
  - one sequential full-width read;
  - the compression ratio and speed of the target layout, measured on one window.

**Overview use is verified, not assumed.** One extra, untimed, single-threaded warp runs
with `CPL_DEBUG` enabled. The warp kernel reports the source window it reads, and a
window smaller than the native one proves that an overview was used. A multithreaded
warp cannot be inspected this way, because its debug messages come from worker threads,
whose error handlers are separate per thread.

**Why layout matters.** In a strip layout (block height below 64, or block width equal to
the raster width) without overviews:
- **window reads** decode whole strips;
- **every warp to the working resolution** decodes every native pixel.

Measured on the synthetic DTM (§8):
- **warp time per window:** 0.74–0.81 s in strip layout, 0.34–0.38 s after conversion;
- **terrain read per unit:** 8–9 s before conversion, 1.7–2.3 s after.

Conversion only runs with `convert_dtm = true`. Both options read the source once,
sequentially, and build `AVERAGE` overviews at factors 2, 4, 8, 16, 32 with the same
compression and predictor:
- **Option A:** the whole file as a tiled GeoTIFF (512 × 512, ZSTD or DEFLATE, predictor
  3 for float or 2 for integer, BigTIFF, all CPUs).
- **Option B:** only the corridor **plus the halo** (the processing domain), written in
  one sequential pass of row blocks sized by `warp_memory_limit_mb`; cells outside are
  nodata.

Before any conversion the script prints, from the measurements:
- **disk space:** raw size × measured compression ratio × 4/3 for overviews; for B, also
  × the domain share;
- **time:** sequential read time plus compression time.

It also checks the free disk space and stops if it is insufficient.

**Which option to choose.** Choose B when the domain covers clearly less than the whole
DTM (55 % on the benchmark world) and disk space is tight. Choose A when the converted
DTM will be reused for other purposes.

### 6.4 Units

- **Cores.** Square, `unit_core_size_m`, aligned to the working grid anchored at the DTM
  origin. Only cores that intersect the corridor are kept.
- **Windows.** Each core plus `unit_halo_m`, clipped to the DTM extent.
- **Domain.** Inside the window, cells outside the corridor plus halo are set to nodata
  and skipped by every algorithm.
- **Unit ids.** `<column>_<row>` in the grid of cores, with rows counted from the north.

**Halo rule: `unit_halo_m ≥ corridor_buffer_m`** (checked in the configuration).
1. A core cell can lie up to `corridor_buffer_m` from its river.
2. The window extends `unit_halo_m` beyond the core in every direction.
3. So the river point nearest to every core cell is inside the window only if
   halo ≥ buffer.
4. Without it, HAND has no valid reference.

The test valley showed the failure mode: windows without the river produced no
compartments at all. `build_units.py` additionally lists windows that contain no river.

**Memory model.** RAM per unit ≈ window cells × `bytes_per_cell`.
- **Check:** `build_units.py` stops if the largest window would exceed
  `max_ram_gb_per_worker` and prints the largest core size that fits.
- **Measured value:** `compartments.py` measures peak RSS per unit and logs bytes per
  window cell: **106–108 B/cell** on the benchmark, after the optimisations.
- **Where it comes from:** the per-unit arrays and the numba kernels' internal arrays.
- **Worker count:** total RAM ≈ `n_workers × max unit RAM`. Also keep
  `n_workers × gdal_threads_per_worker` at or below the number of cores.

## 7. Parallel processing, checkpoints, merging

### 7.1 Pool

`gis.Run.run_units` uses a `ProcessPoolExecutor` with the `spawn` start method and
`max_tasks_per_child = 1`:

- **A fresh process per unit.** Peak RSS is measured per unit rather than accumulated,
  memory is returned to the system after each unit, and Linux and Windows behave the
  same.
- **Per-process GDAL handles.** Each worker opens its own GDAL datasets, because GDAL
  handles cannot be pickled, and sets `GDAL_CACHEMAX` and `GDAL_NUM_THREADS` from the
  configuration in the initializer.
- **Compilation once.** The numba kernels are compiled in the main process before the
  pool starts, so workers load them from the cache.
- **Submission order.** Units are submitted largest first (estimated domain cells).
- **ETA.** It is based on the domain cells finished so far.
- **Failures.** A crashed worker, e.g. one killed for lack of memory, is reported, and
  the remaining units go on. A rerun retries failed units only.

### 7.2 Checkpoints

Each unit writes its outputs to `work_dir/<country>/units/<unit_id>/` and finishes with
`done.json`, written atomically through a temporary file and a rename. `done.json` holds:
- timings per step,
- peak RSS and bytes per cell,
- counts.

Reruns skip units that have `done.json` unless `FORCE = True`. `UNITS` restricts a run to
the listed units. Merging always uses **all** finished units, so partial runs can be
combined step by step.

### 7.3 Merging

The same compartment can be found by several units, because windows overlap by twice the
halo.

1. **Grouping.** Candidates from different units that overlap by at least
   `merge_overlap_ratio` of the smaller candidate's area form one group (connected
   components of the overlap graph).
2. **Groups with a complete candidate.** If any candidate in a group is not truncated, the
   group keeps the candidate whose representative point (`point_on_surface`) lies in its
   own unit's core. If several qualify, it keeps the one farthest from its window edge.
   This implements and refines the "representative point in the core" rule: a
   compartment straddling a core boundary is seen complete by both units and is kept
   exactly once (test `A_seam`).
3. **Groups that are truncated in every window** (compartments longer than a window):
   - the pieces are **united** (recall);
   - the protection level, entry and source come from the piece with the **lowest**
     protection level. A truncated piece can only overestimate protection, because its
     window edge acts as a wall.
   - The result keeps `truncated = True`. Test `G_long`: 2.8 km long, pieces from several units,
     merged into one.
4. **Contradiction rule** (`drop_contradicted_truncated`):
   - **What it catches:** window edges can create artificial protection. Example: an
     area behind a levee that is longer than the window looks enclosed because the
     window edges close it, while in reality water flows around the levee ends.
   - **Rule:** a compartment that is truncated in every window is dropped if finished
     units whose windows see the same area found **no** compartment there, and their
     windows cover at least `contradiction_ratio` of its area.
   - **Why this is safe:** those units see more of the levee, and a window edge can only
     add walls, never remove protection.
   - **In the test:** this removed exactly one artifact; true long compartments are kept,
     because every window finds a piece of them.
   - The number dropped is logged.
5. **Identifiers.** `compartment_id = <country>_<x>_<y>`, from the representative point
   rounded to `id_rounding_m`. It is stable between runs as long as the geometry is
   stable. Duplicates get a suffix `_2`, `_3`.
6. **Report.** The log lists the number of truncated compartments and their share. Many
   truncated compartments mean the halo is too small.

Barrier crest statistics (`crest_min_m`, `crest_hand_min_m` in `barriers.gpkg`) are
sampled inside each unit's **core** only and combined as the minimum over units. Every
line section is therefore evaluated in exactly one unit.

## 8. Measured performance

### 8.1 Where the numbers come from

All numbers are from this cloud container: 4 vCPUs, 15 GB RAM, virtual disk. The OS page
cache was dropped before every timed run (`echo 3 > /proc/sys/vm/drop_caches`).

The data is the synthetic DTM from `tests/make_benchmark_world.py`:
- two GeoTIFFs of 12 000 × 24 000 px at 1 m (1.15 GB raw each, 0.54 GB DEFLATE);
- strip layout with blocks of 12 000 × 1, no overviews;
- a meandering river, two tributaries, 27 levee lines (seven rings with embankments).

Working resolution 2 m, cores 7 km, halo 2 km. That gives windows of
5500 × 5500 = 30.25 M cells, close to a real unit of 20 km + 3 km halo at 5 m
(5200 × 5200 = 27 M).

**Disk and decompression speeds on your machines will differ.** These numbers show the
*structure* of the costs, not your runtimes.

### 8.2 DTM reads (20 windows per file)

| Layout | Native 4096² window | `gdal.Warp` to 2 m (same ground) | Sequential read | Overview used |
|---|---|---|---|---|
| Strip, no overviews | 0.31–0.45 s (148–214 MB/s) | 0.74–0.81 s (82–90 MB/s native equivalent) | 530–640 MB/s | no (source 4096 px) |
| Option A (tiled 512, ZSTD, overviews) | 0.14–0.20 s (342–468 MB/s) | 0.34–0.38 s (178–197 MB/s) | 368–501 MB/s | yes (source 2049 px, factor 2) |

Conversion results:
- **Size.** ZSTD with predictor 3 compressed the synthetic DTM to 0.47 of its raw size.
  This includes 3 cm of noise; real lidar may compress differently, and the benchmark
  measures it per file.
- **Time.** Option A took 8–11 s per 1.15 GB file plus 5–6 s for overviews; each output
  is 0.71 GB.

### 8.3 One unit (30.25 M window cells, 24.4 M domain cells, cold cache)

| Step | Strip source | Option A source |
|---|---|---|
| Terrain read (`gdal.Warp` average) | 8.0–9.4 s | 1.7–2.3 s |
| Crest boxes (`gdal.Warp` max, 6–13 boxes) | 0.7–1.7 s | 0.6–1.1 s |
| HAND (taper burning, breaching, D8) | 29–31 s | 31–33 s |
| Two priority floods | 26–31 s | 28–33 s |
| Compartment statistics | 2–3 s (after optimisation; 18 s with `scipy.ndimage` label statistics) | |
| Write (two GeoTIFFs, GeoPackages) | 1.6–3.2 s | |
| **Total** | **95–99 s** (before optimisation) | **76–81 s** |
| Peak RSS | 4.2 GB | **3.0 GB (106–108 B per window cell)** |

**Runtime per country:**
- **Formula:** worker time ≈ domain cells × **3.2–4.0 µs** per cell; wall time ≈ worker
  time ÷ `n_workers`. Domain cells per country are in `units.csv` (`est_domain_cells`).
- **Example:** at 5 m, 1000 km² of domain is 40 M cells, about 2.2–2.7 worker-minutes.
- **Measured split:** HAND and the two floods take about 80 % of the compute once the
  DTM is tiled. Terrain reading takes about 10 % from a strip-layout source and 2–3 %
  after conversion.
- **Your machine:** `compartments.py` logs this extrapolation from your own test units
  and writes it to `performance.json`.

## 9. Tests

`pixi run test-kernels` runs 16 checks of the numba kernels:
- closed ring, gap, diagonal (staircase) barrier;
- nested rings, both crest orders;
- seed and order propagation, entry cell;
- breaching, nodata;
- bug B1 reproduced without burning and fixed with burning.

`pixi run test-pipeline` runs 28 checks on a synthetic 3 × 2 km country. The DTM is two
strip-layout files at 1 m; the working resolution is 2 m with 1 km cores and a 500 m halo:
- `prepare_dtm` (inspection, footprint, benchmark) and conversions A and B, giving tiled
  files with 5 overviews;
- `build_units`;
- unit selection, checkpointing (resume, skip);
- a ring straddling a core boundary appears exactly once;
- a bridged 10 m gap and a 40 m gap (no protection);
- a staircase levee (barrier rasters 4-connected, no leak);
- nested rings, inner level = higher crest (to 0.01 m);
- a compartment longer than any window, merged from truncated pieces, plus the dropped
  window-edge artifact;
- unique ids; source cells on the river, entry cells in the compartment;
- an interpolated nodata hole over the river;
- empty levee input: margin exactly 0 everywhere, no compartments.

`DEFENSES_TEST_DIR` keeps the test folder for inspection.

`python -m tests.make_benchmark_world` writes the performance world and its config
(`DEFENSES_BENCH_DIR`).

## 10. Bugs found, limitations, open decisions

### 10.1 Bugs in the existing method (reported; decision needed)

**B1. HAND drainage paths run parallel to the rasterized river.** (Fix implemented as an
option; default = previous behaviour.)
- **Cause.** The river is rasterized along its centre line. Channel cells beside it are
  just as low, and on lidar the centre line often misses the thalweg. D8 then prefers the
  down-valley slope along the channel and reaches the river only at the lowest river cell
  in the window.
- **Effect.** The whole window takes the river elevation at its downstream edge as
  `river_ref`. HAND is too high by up to the river slope × window length; for a 26 km
  window at 0.5 m/km that is up to 13 m.
- **How it was found.** In the test valley this produced a 0.7 m error in nested-ring
  protection levels. `tests/test_kernels.py` reproduces it.
- **Fix.** Tapered stream burning in the conditioned copy only
  (`stream_burn_m`, `stream_burn_buffer_m`; AGREE-style). Elevations are unaffected.
  The template keeps `stream_burn_m = 0.0` (previous behaviour, with a warning).
- **Proposed setting:** `stream_burn_m = 5.0`, with `stream_burn_buffer_m` larger than
  half the widest modelled channel (e.g. 100 m).

**B2. Levee removal can overwrite river cells.** (Reported, not changed.)
- **Cause.** The removal zone (barrier dilated by 10 m) can reach the river channel next
  to a levee. `FillNodata` then interpolates the channel bed from the banks, which raises
  `river_ref` locally.
- **Proposed fix.** Exclude river cells, plus one cell, from the removal zone.

**B3. DTM nodata at the coverage edge acts as a wall.** (Reported, not changed.)
- **Effect.** An area enclosed by levees and missing DTM coverage can appear protected.
- **Proposed fix.** Flag such compartments (`touches_dtm_edge`).

**B4. Gated mode ignores rivers leaving the DTM coverage.** (Reported, not changed.) A
river that leaves the DTM coverage inside a window is not recognised as entering from
outside; only the window border counts.

### 10.2 Problems found and fixed in the new code during Stage A

- Option B first clipped the copy to the corridor without the halo, so units read nodata
  in the halo. It now clips to corridor + halo.
- The search radius for filling nodata holes was derived from all nodata areas, including
  the huge area outside the coverage. It now uses only the enclosed holes.
- The overview check reported an overview *scan* as overview *use*. It now parses the
  source window size.
- A configuration with `unit_halo_m < corridor_buffer_m` silently produced units without
  rivers. It is now a validation error (§6.4).
- `scipy.ndimage` label statistics scanned the whole window: 18 s → 3 s per unit.

### 10.3 Limitations

- **L1. Water-surface nodata holes.**
  - **Effect.** In lidar, water is often nodata across the full river width. Holes are
    filled by inverse-distance interpolation from the banks, which raises the river-cell
    elevations. In the synthetic test this was about +1 m, and the resulting jump in
    `river_ref` produced a spurious compartment with a 0.67 m margin.
  - **Test setup.** The test hole is therefore narrower than the channel bed.
  - **Proposal.** Cap filled river cells by the local minimum of the surrounding valid
    cells (water level ≈ lowest bank).
- **L2. Mixed river references.** HAND levels of cells draining to different rivers are
  compared in one flood. Near confluences, and where a tributary crosses a main-river
  levee through a gate, the margin can be underestimated. This is inherent to HAND; a
  per-origin reference was tested in an earlier session and created false protection
  downstream of floodplain obstacles.
- **L3. Truncated compartments.** Where no neighbouring window contradicts them, window
  edges can still close areas artificially. They are flagged `truncated`; their share is
  reported.
- **L4. Crests from the DTM.** Floodwalls narrower than the DTM resolution are
  underestimated. The crest search can take a crest up to one working cell beyond the
  configured radius.
- **L5. Platform coverage.** Everything was tested on Linux. On Windows, the spawn pool
  and the `peak_rss_mb` implementation through the Windows API are written but untested.
- **L6. Processing order.** The levee merge reads all levee lines within the bounding box
  of all unit windows into memory at once (vector data, small compared with the DTM).

### 10.4 Issues in the existing `gis.py` (reported, not changed)

| Function | Issue |
|---|---|
| (module import) | Imports `pyflwdir`, added to `pixi.toml` |
| (module import) | Sets `gdal.SetCacheMax(4 GB)`; workers override it from the config |
| `Raster.update_nodata_vrt` | Calls `FlushCache()` on `None` |
| `Vector.load_ogr_ds_to_geodataframe` | Returns a tuple because of a trailing comma |
| `Raster.raster_to_vect`, `raster_to_mem_vector` | Fail with `filter_values=None` and with a path argument (`drivername` undefined) |
| `Files.get_root_path` | Loops forever without an `IF-HydroSim` parent folder |
| `Vector.extract_intersection` | Creates an exception without `raise` |
| Several functions | Assume EPSG codes; the new functions use WKT |

## 11. Implementation stages

- **Stage A (done):** `gis.py` additions, `prepare_dtm.py`, `build_units.py`,
  `compartments.py` refactor; tests; synthetic benchmark.
- **Stage B:** `download_osm.py`, `osm_hydraulics.py` (water sources, `osm_pumping.gpkg`),
  `download_ghsl.py`, `rp_index.py`.
- **Stage C:** `download_gfm.py`, `gfm_evidence.py`, `download_glofas.py`,
  `rp_bounds.py`, `validate.py` (moved from `validate_defenses.py`, with RP confusion
  matrix).

## 12. Data sources and licences

The download URLs are **not verified yet**. This session's network blocked the official
hosts, so the config entries are empty and marked `NOT VERIFIED`. They will be filled
from official documentation in Stages B and C, with release, citation and licence:
- OSM via Geofabrik (ODbL),
- GHSL built-up volume R2023A (JRC; licence to be confirmed),
- CEMS GFM via its STAC API,
- GloFAS historical discharge via EWDS.

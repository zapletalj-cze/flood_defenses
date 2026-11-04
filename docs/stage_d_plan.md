<!-- Author: Jakub Zapletal -->
<!-- Date: 2026-11-04 -->
<!-- Stage D, step D.1: repository inventory, optional inputs, plan and measured runtime impact (for approval). -->

# Stage D, step D.1: inventory, inputs and plan

Status: **approved.** Step D1 is implemented (README §17); D2–D5 follow. The prototype kernels of section 4 were
written only in a scratch folder, to measure the runtime figures.

---

## 1. Repository state after stages A–C: what differs from the main prompt

### 1.1 Form and layout
- **Run variables instead of command-line flags.** The scripts take run variables (`FORCE`, `UNITS`, `HINTS`,
  `DATE_FROM`, `DATE_TO`) that `DEFENSES_<NAME>` environment variables override, as you asked. `--hints` of
  `gfm_evidence.py` became `HINTS`.
- **Execution.** Scripts run as modules (`python -m processing.x`), so `import gis` works without `sys.path`
  changes.
- **Moved script.** `validate_defenses.py` became `processing/validate.py`.
- **Additions outside the main prompt:**
  - `download_nld.py` (earlier request);
  - `processing/HAND_YZ.py` and `download/download_worldcover.py` (land cover exposure supplementing GHSL);
  - `docs/navod_cs.md` (short Czech guide to the scripts).
- **Tests.** 109 checks: kernels 16, pipeline 40, Stage B 20, Stage C 26, HAND_YZ 7. All pass from scratch.

### 1.2 Configuration
80 keys beyond the 57 of the original schema (137 in total), each listed with its reason in README §3.
- **Stage A.** DTM conversion, footprint resolution, memory model, merging and contradiction rules, gap filling,
  stream burning (fix B1), river cell preservation (fix B2).
- **Stage B.** OSM region and tag lists, including the closing and pumping subsets, and GHSL templates as a table.
- **Stage C.**
  - **GFM:** collection, asset keys, page size, token variable, screening resolution, download threads, pixel
    codes.
  - **GloFAS:** request table, time and area keys, years, local files, variables, upstream area, day convention.
  - **Hints:** keys including `hint_exclusion_m`.
  - **Validation:** `validation_extent`.
- **HAND_YZ.** The `[countries.XX.lulc]` table and the WorldCover sources.

### 1.3 Method changes and interpretations
- **compartments.py**
  - **Bug fixes B1–B4:**
    - tapered stream burning;
    - river cells kept in the levee removal zone;
    - `touches_dtm_edge` flag;
    - gated rivers entering from missing coverage.
  - **Gap filling for water surfaces:** coverage closing, river profile interpolation and a cap at the waterline.
  - **OSM water sources.** Closed streams are excluded from both the seeds and the HAND network.
  - **Entry and source** come from the reach order of the defended flood.
  - **New columns:** `unit_label`, `touches_dtm_edge`.
  - **Merging.**
    - A complete candidate inside its own core wins.
    - Pieces truncated everywhere are united and take the lowest level.
    - A contradiction rule drops truncated compartments that better-placed units do not see.
  - **Barrier lines** are bridged **once per country** in `compartments.py`; Stage D moves this into
    `barriers_merge.py`.
- **rp_index.py**
  - totals check with exact areas;
  - valuation only in the bounding box of each unit's compartments;
  - Stage C columns created empty;
  - `value_ghsl + value_lulc` (HAND_YZ).
- **Stage C.**
  - **`download_gfm.py`** screens each scene once for all units, in threads, with a checkpoint (the prompt implied
    a screening per unit).
  - **Hints:** connectivity through permanent water and river lines; new exclusion distance around levees.
  - **`rp_bounds.py`:** the combination rules are my interpretation (README §11.5).
  - **GloFAS:** read from local files; `cdsapi` is not installed.
- **Not verified.** These entries are empty and flagged:
  - Geofabrik and GHSL URL templates;
  - GFM endpoint, collection, assets and pixel codes;
  - GloFAS upstream area URL;
  - Urban Atlas class field and codes.

### 1.4 Where Stage D touches existing code

| Stage D item | Current state |
|---|---|
| `barriers_merge.py` | Merging and bridging live in `compartments.py` (`work/<c>/barrier_lines.gpkg`, fields `source`, `bridged`) |
| `rasterize_lines_max` | Barriers are rasterized with `gdal.RasterizeLayer` (`ALL_TOUCHED`), diagonal contacts are closed, crests are read as the native maximum in boxes per unit |
| Controller propagation | `gis.Hydro.priority_flood` returns spill, seed and order (numba); compartments are 4-connected protected regions |
| Per-context units | One working resolution and core size per country |
| Population, context distributions, `defenses_clip.gpkg` | One country distribution; no population; no clip layer |
| GFM context, EMS | Not present |
| OSM extensions | Points and drains only |
| Urban Atlas | Read by `HAND_YZ.py` (`Vector.single_layer_copy` for multi-layer files), reusable for `context_zones.py` |

---

## 2. Optional inputs per country (D.1 item 2)

- **No country data is available in this session.** The repository holds only the template country `XX`, with
  every optional input empty: `dsm`, `building_footprints`, `manual_barriers`, `industrial_registry`,
  `ems_flood_extents` and `water_polygons`.
- **No access to your machine either.** No computer folder is shared with this session (Remote Control), so I
  cannot inspect or benchmark your DSMs.
- **Proposal.** In D1, `prepare_dtm.py` treats `dsm` exactly like `dtm`: same report (size, pixel size, type,
  nodata, compression, block size, overviews, file size) and the same benchmark (20 native 4096² windows, the same
  windows through `gdal.Warp`, a sequential read). You can then:
  - run it and send me `dtm_inspection.csv` and `dtm_benchmark.csv`, or
  - share the data folder through Remote Control (Claude desktop app: Settings, Claude Code, Remote Control), so
    I can run it.
- **Reference until then.** A DSM with the same storage layout as the benchmark DTM reads like it (README §13.2).
  - **Strip layout:** 0.74–0.81 s per 4096² warp window.
  - **Tiled with overviews:** 0.34–0.38 s.

---

## 3. Plan

### 3.1 Files per step

| Step | File | Change |
|---|---|---|
| D1 | `download/download_osm.py` | adds `osm_buildings.gpkg`, `osm_water.gpkg`, `osm_crossings.gpkg`, `osm_walls.gpkg` (one SQL per output); logs counts per tag value |
| D1 | `download/download_ghsl.py` | adds SMOD and GHS-POP of the configured epoch (`ghsl_smod.vrt`, `ghsl_pop.vrt`) |
| D1 | `download/download_clms.py` (new) | CORINE raster and Urban Atlas from **verified** sources only; otherwise it validates and prepares local files (CORINE class index → code through the legend, checked on a known pixel; single-layer Urban Atlas copies) |
| D1 | `processing/context_zones.py` (new) | `zones.gpkg` (field `context`, order `industrial` > `dense_urban` > `urban` > `rural`, small zones dropped) |
| D1 | `processing/prepare_dtm.py` | inspection and benchmark of `dsm` |
| D2 | `processing/crest_lines.py` (new) | `crest_lines.gpkg` (transects, numba detection, building and roughness exclusions, linking) |
| D2 | `processing/barriers_merge.py` (new) | `barriers_all.gpkg` (dedupe by priority, crest Z, closure gaps at OSM crossings, regular gaps) |
| D3 | `processing/compartments.py` | reads `barriers_all.gpkg`; `rasterize_lines_max`; compartments by controller cell; `parent_id`, `nesting_level`; context fields; per-context resolution and core size; reconciliation with the coarse pass |
| D3 | `processing/build_units.py` | units of `dense_urban` and `industrial` use their context's resolution and core size |
| D3 | `processing/coarse_connectivity.py` (new) | corridor at `coarse_resolution_m`, split by connected corridor components if RAM is short; `compartments_coarse.gpkg` |
| D4 | `processing/rp_index.py` | population (GHS-POP, totals preserved), per-context weights and distributions (country fallback below 30 compartments), `has_regulated_site`, nesting monotonicity, `defenses_clip.gpkg` |
| D4 | `processing/review_queue.py` (new) | `review_queue.gpkg`, `review_queue.csv` |
| D5 | `processing/gfm_evidence.py` | `context` per test, `no_test_urban`, EMS scenes (`source`) |
| D5 | `processing/raised_sites.py` (new) | `raised_sites.gpkg` (diagnostic only) |
| D1–D5 | `tests/test_stage_d.py` (new) | the checks of D.6, added step by step |
| D5 | `README.md`, `docs/navod_cs.md` | method, reasoning, review loop, sources, licences (the Czech guide only gets the new steps) |

### 3.2 gis.py additions

All additions are extensions. Existing calls keep their behaviour.

| Function | Purpose |
|---|---|
| `Hydro.rasterize_lines_max(shape, geotransform, x0, y0, z0, x1, y1, z1)` | numba grid traversal: cells form a 4-connected path, Z is interpolated along segments, the maximum is kept where lines overlap |
| `Hydro.priority_flood(surface, seeds, controller=False)` | with `controller=True` also returns `ctrl` (`ctrl(n) = n` where the surface lifts the level, else `ctrl(c)`); default output unchanged |
| `Hydro.detect_crests(profiles, step, max_bank_cells, min_drop, landward_cells, width_drop)` | numba: first and second crest candidate per transect, height, riverward drop, width |
| `Raster.sample_points(raster, x, y, chunk)` | nearest sampling of many points with one `ReadAsArray` per window of the raster's own block layout |
| `Raster.max_along_lines(raster, lines, radius)` | native maximum within a radius, sampled along lines at pixel spacing (crest Z of ML and OSM lines) |
| `Raster.map_legend(array, lookup)` | class index → code (CORINE) |
| `Vector.zones_to_grid(zones, order, shape, geotransform, wkt)` | context raster for a unit (earlier context wins) |

### 3.3 Configuration
- **New tables.** The tables of D.2, verbatim. The context parameter tables (`[countries.XX.context.params.<context>]`)
  exist for `industrial`, `dense_urban`, `urban` and `rural`; `rural` defaults to today's country values.
- **Additional keys**, with placeholder values:
  - `rp.min_context_compartments = 30`;
  - `evidence.urban_min_valid_share = 0.5`;
  - `[raised_sites] height_m`, `ring_width_m`.
- **Unverified codes stay as placeholders.** CORINE, Urban Atlas and SMOD codes are logged as placeholders until
  verified:
  - **SMOD:** class codes 10, 11, 12, 13, 21, 22, 23 and 30, as listed in the Earth Engine catalogue entry of GHS-SMOD P2023A. Your
    placeholders `[30]` and `[21, 22, 23]` agree with it; please confirm them in the GHSL data package.

---

## 4. Runtime impact, from measurements

**Environment and data.** 4 vCPU and 15 GB RAM, on the synthetic benchmark world of Stage A (24 × 24 km at 1 m,
strip layout and option A copies).
- **Cold reads.** For each run, the files were evicted from the page cache (`posix_fadvise`), because
  `drop_caches` is not writable in this container.
- **Prototypes.** The kernels are prototypes; the final code is re-measured in each step.

| Item | Measurement | Impact |
|---|---|---|
| Controller propagation | flood on 30.25 M cells: 11.35 s without, 11.64 s with (+2.6 %); spill identical | +0.6 s per unit (two floods) |
| `rasterize_lines_max` | 1652 segments (57 km): 0.4 ms against 39.5 ms for GDAL `ALL_TOUCHED`; same 15 834 cells; no diagonal-only contacts | negligible; replaces the GDAL burn and the diagonal closing |
| Transect sampling and detection (300 m, 1 m samples, both banks, one window per 500 m of river) | 9430 transects: strip layout 1520–1770 /s, tiled (option A) 12 000–12 250 /s; detection 0.02 s in total | 1000 km of river, both banks every 10 m = 200 000 transects: about 2 min (strip) or 17 s (tiled) per raster; the DSM doubles it |
| Coarse terrain read at 25 m | 576 km²: 21.9 s strip without overviews (38 s per 1000 km²), 0.61 s tiled with overviews (1.1 s per 1000 km²) | 20 000 km² of corridor: 13 min (strip) or 21 s (tiled) |
| Coarse kernels at 25 m | measured per cell: HAND 0.46 µs, flood 0.38 µs (+3 % with controller) | 20 000 km² = 32 M cells: about 40 s; RAM at 95 B/cell: 3 GB, so corridors above about 60 000 km² are split by components |
| Fine units at 2 m in `dense_urban` and `industrial` | `compartments.py` 2.4 µs per domain cell, measured at 2 m | 0.6 s per km² of domain, against 0.1 s at 5 m; 1000 km² of such domain = 10 worker-minutes |
| Population warp in `rp_index.py` | GHSL warps measured at 1.6–4.6 s per unit for two products | about +1–2 s per unit |
| `context_zones`, `barriers_merge`, `review_queue` | vector operations; not measured yet | measured in D1, D2 and D4 on the synthetic world before any country estimate |

---

## 5. Points to decide before D1

1. **Data access.** DSM, building footprints, manual barriers and an industrial registry are not in this session.
   Please share a folder or run the D1 inspection yourself.
2. **CLMS (CORINE, Urban Atlas).** The CLMS and EEA hosts are blocked here, and their downloads need a login.
   - **Proposal:** `download_clms.py` prepares and validates local files.
   - **API download:** only after you approve the method and it can be verified.
3. **GHSL SMOD and GHS-POP.** Same download path as `BUILT_V`; the URL templates are still unverified, so local
   files are used.
4. **Industrial registry** (EEA industrial reporting). Its licence needs checking; it is used from a local file.
5. **OSM flood wall tags.** taginfo is blocked here. Proposal: `download_osm.py` counts the candidate tags in the
   country extract itself, e.g. `barrier=wall` with `wall=flood_wall`, or `man_made=dyke`, and
   `osm_tags_walls` is filled from those counts.
6. **Compartment definition changes in D3.** Compartments sharing a controller replace connected protected
   regions. Nested systems split into parent and child, so some Stage A expectations change, and the affected
   tests will be updated, not loosened.
7. **HAND_YZ and context zones.** Land cover value stays in `HAND_YZ.py` (WorldCover, Urban Atlas). Zoning comes
   from CORINE, Urban Atlas, SMOD and OSM as specified. Both share the Urban Atlas reader.
8. **Stops.** After each of D1–D5, with results and measured timings.

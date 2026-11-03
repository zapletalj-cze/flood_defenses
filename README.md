# flood_defenses

Defenses layer for the fluvial flood hazard pipeline: areas protected by flood levees, derived from terrain,
levee lines and the modelled river network only (no hazard maps or hydraulic model output).

## Scripts

`defenses.py` derives the protected compartments.

```
pixi run python defenses.py --dtm tile1.tif tile2.tif --levees levees_a.gpkg levees_b.shp \
    --rivers rivers.gpkg --extent pilot.gpkg --output-dir out [--resolution 5] [--keep-intermediate]
```

Parameters (defaults): `--gap-tolerance-m 20`, `--crest-search-radius-m 5`, `--barrier-removal-width-m 10`,
`--min-margin-m 0.3`, `--max-level-m 20`, `--min-area-m2 10000`, `--resolution` (native),
`--buffer-m 1000` (context around `--extent`), `--gated-crossings` (off: a river crossing a levee wins).

Outputs in `--output-dir`:
- `protected_areas.gpkg`: `id`, `area_m2`, `protection_level_m`, `margin_median_m`, `margin_max_m`,
  `river_ref_median_m`, `edge` (touches the processing boundary), `rp` (empty).
- `enclosed_streams.gpkg`: areas enclosed by levees that contain modelled river cells, with `area_m2`,
  `stream_cells`, `gated`.
- `barriers.gpkg`: merged and bridged levee lines with `source`, `bridged`, `crest_min_m`, `crest_hand_min_m`.
- with `--keep-intermediate`: `dtm_undefended.tif`, `hand_undef.tif`, `river_ref.tif`, `spill_undef.tif`,
  `spill_def.tif`, `margin.tif`, `barrier.tif`, `crest.tif`.

`validate_defenses.py` compares the result with reference protected-area polygons.

```
pixi run python validate_defenses.py --result out/protected_areas.gpkg --reference reference.gpkg \
    --output out/validation.gpkg [--extent pilot.gpkg]
```

`download_nld.py` downloads the USACE National Levee Database from its public OGC WFS service with GDAL, one
GeoPackage per layer (levee centrelines and floodwalls for `--levees`, leveed areas for `--reference`).

```
pixi run python download_nld.py --output-dir nld --list
pixi run python download_nld.py --output-dir nld [--layers embankment floodwall leveed] [--extent pilot.gpkg]
```

`tests/synthetic_checks.py` runs both scripts on synthetic terrain (`pixi run python tests/synthetic_checks.py`).

## Method

1. The DTM is read for the extent bounding box plus `--buffer-m`. Interior nodata holes (e.g. water surfaces in
   lidar) are filled with `gdal.FillNodata`; nodata touching the grid edge is treated as missing coverage.
2. Levee lines are merged and split into single parts. Each endpoint closer than `gap_tolerance_m` to another
   line is connected to it: to that line's endpoint when it is within the tolerance, otherwise to its nearest
   point. A single line whose own ends are within the tolerance (and which is longer than 3 x tolerance) is
   closed. All such segments have `bridged=True`.
3. Barriers are rasterized with all-touched; diagonal-only contacts are closed by adding the higher corner cell.
   Crest elevation is the maximum of the native DTM within `crest_search_radius_m`; bridged cells interpolate
   linearly between the endpoint crests.
4. Zones separated by barriers are labelled. River cells in zones whose streams do not reach the grid edge are
   enclosed by levees and written to `enclosed_streams.gpkg`. By default a river crossing a levee removes the
   barrier there (river wins), so such zones are unprotected; with `--gated-crossings` the crossing stays closed
   and the enclosed streams are not water sources.
5. The barrier mask dilated by `barrier_removal_width_m` is filled with `gdal.FillNodata` (undefended DTM).
6. Flow routing on the undefended DTM: complete breaching (priority flood from the stream cells, depressions are
   carved along the flood path), then D8 steepest descent. `river_ref` is the unconditioned elevation of the stream
   cell reached; `hand_undef = dtm_undefended - river_ref`.
7. 4-connected priority flood from the water-source stream cells on `hand_undef` and on the defended surface
   (barrier cells set to `crest - river_ref`) gives `spill_undef` and `spill_def`.
8. Protected: `spill_def - spill_undef >= min_margin_m` and `spill_undef <= max_level_m`, barrier cells excluded.
   4-connected regions smaller than `min_area_m2` are dropped; with `--extent`, compartments intersecting it are
   kept whole.

## Limitations

- HAND levels are relative to each cell's own drainage stream. Where two rivers meet at a levee (e.g. a tributary
  crossing a main-river levee through a gate), the crest and the area behind it can be measured against different
  rivers, which underestimates the margin there. Measuring every flood path against its source river instead
  removes this but creates false protected "shadows" downstream of floodplain obstacles, so it is not used.
- Only modelled rivers are outlets of the flow routing; side catchments are routed (carved) to a modelled river.
  This is intentional (every cell gets a modelled-river reference); use `--buffer-m` so that their real drainage
  path lies inside the processed grid.
- Crest elevations come from the DTM. Floodwalls narrower than the DTM resolution are underestimated.
- Processing is in memory (no tiling); `rp` is not assigned yet; transport embankments are not an input.
- `download_nld.py` relies on the NLD WFS service; its layer names and paging behaviour should be checked with
  `--list` on first use.

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
`--min-margin-m 0.3`, `--max-level-m 20`, `--min-area-m2 10000`, `--resolution` (native).

Outputs in `--output-dir`:
- `protected_areas.gpkg`: `id`, `area_m2`, `protection_level_m`, `margin_median_m`, `margin_max_m`,
  `river_ref_median_m`, `rp` (empty).
- `barriers.gpkg`: merged and bridged levee lines with `source`, `bridged`, `crest_min_m`, `crest_hand_min_m`.
- with `--keep-intermediate`: `dtm_undefended.tif`, `hand_undef.tif`, `river_ref.tif`, `spill_undef.tif`,
  `spill_def.tif`, `margin.tif`.

`validate_defenses.py` compares the result with reference protected-area polygons.

```
pixi run python validate_defenses.py --result out/protected_areas.gpkg --reference reference.gpkg \
    --output out/validation.gpkg [--extent pilot.gpkg]
```

`tests/synthetic_checks.py` runs both scripts on synthetic terrain (`pixi run python tests/synthetic_checks.py`).

## Method

1. Levee lines are merged and split into single parts. Each endpoint closer than `gap_tolerance_m` to another
   line is connected to it: to that line's endpoint when it is within the tolerance, otherwise to its nearest
   point (`bridged=True`).
2. Barriers are rasterized with all-touched; diagonal-only contacts are closed by adding the higher corner cell.
   Crest elevation is the maximum of the native DTM within `crest_search_radius_m`; bridged cells interpolate
   linearly between the endpoint crests. River cells are removed from the barrier.
3. The barrier mask dilated by `barrier_removal_width_m` is filled with `gdal.FillNodata` (undefended DTM).
4. Flow routing on the undefended DTM: complete breaching (priority flood from the stream cells, depressions are
   carved along the flood path), then D8 steepest descent. `river_ref` is the unconditioned elevation of the stream
   cell reached; `hand_undef = dtm_undefended - river_ref`.
5. 4-connected priority flood from the stream cells on `hand_undef` and on the defended surface (barrier cells set
   to `crest - river_ref`) gives `spill_undef` and `spill_def`.
6. Protected: `spill_def - spill_undef >= min_margin_m` and `spill_undef <= max_level_m`, barrier cells excluded.
   4-connected regions smaller than `min_area_m2` are dropped.

## Limitations

- Stream cells inside a levee-enclosed area (a modelled tributary behind a levee, or a river crossing a levee)
  make the whole area unprotected; such areas are reported as warnings in the log.
- Only modelled rivers are outlets of the flow routing, so flow is routed (carved) to a modelled river even across
  divides; cells not connected to any stream cell take the nearest stream cell.
- HAND values relative to different rivers are compared in one priority flood, which is approximate near
  confluences of rivers with different water levels.
- Processing is in memory on the bounding box of the extent; compartments touching its edge may be incomplete.
- `rp` is not assigned yet; transport embankments are not an input.

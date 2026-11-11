# Review of the reports (Stage G, step 1)

Author: Jakub Zapletal
Date: 2026-11-11

## Scope
- **Read:** the only report in the repository, `reports/POL/report.md` and
  `reports/POL/overview.png` (run of 2026-11-06). Both are unchanged.
- **Also checked:** the code paths that explain the numbers (`processing/compartments.py`:
  barrier bridging, water sources, river crossings).

## 1. What the report shows

### Inputs
| Item | Value |
|---|---|
| DTM | `PL_DEM_10m_borders_20220913.tif`, 83,187 × 73,848 px, 10 m, EPSG:2180, Float32 |
| DTM storage | strip layout (block 83,187 × 1), no overviews, LZW |
| DTM read speed | 3.5 s per 4096² window (19 MB/s); option A conversion estimated at 18 GB and about 2 min |
| Footprint | 366,546 km² (Poland plus a border buffer); 38 % nodata in random windows (outside the coverage) |
| Levees | the OSM Europe file `man_made=embankment or dyke`, **not a national levee dataset** |
| Levees in the footprint | 77,408 lines, 22,920 km |
| Levee tags | `man_made=embankment` 73,906 lines (95 %), `dyke` 3,502 (4.5 %) |
| Levee length | median 125 m per line |
| Rivers | MERIT, 8,268 lines in the footprint, 69,160 km, `uparea` 25 km² to 194,369 km² |
| `river_upstream_area_field` | still `'50'`; the column is `uparea` |
| OSM hydraulics | `rivers_sources.gpkg` and `osm_pumping.gpkg` were not produced (the `osm_hydraulics` log has 0 lines) |

The levee file also contains lines that are clearly not flood levees: named `Transzeja` (trenches, 99), `Niemieckie okopy z II WŚ` (WWII trenches), lines with `highway=track/path/motorway`, `barrier=fence/wall` and `waterway=dam`.

### Outputs
| Item | Value |
|---|---|
| Units | 971 units of 20 km core and 3 km halo at a 10 m working cell; 22 min with 16 workers, 114 B per cell |
| Compartments | **3,255, 824.9 km² in total**, in 440 of the 971 units |
| Size per compartment | median 4 ha, p75 13 ha, p95 116 ha, max 11.6 km² |
| `protection_level_m` | median 3.2 m, p95 16.4 m, max 27.3 m |
| `margin_median_m` | median 0.70 m; p5 0.34 m, just above `min_margin_m` 0.3 |
| Flags | truncated 21, corridor edge 2, DTM edge 2, pumping 0 |
| Barriers | 151,273 segments, 26,050 km |
| Barrier gaps | 62,141 segments (41 %) are bridged gaps added by the gap bridging |
| `crest_hand_min_m` | median 13.9 m, p75 34 m, p95 108 m; null for 19,619 segments |
| `enclosed_streams.gpkg` | 894 levee-enclosed zones with river cells, 55.7 km² (largest 28.4 km²), all treated as unprotected ("river wins") |
| Run time | flood 25 %, crest read 24 %, HAND 23 %; crest reading is slow because of the strip layout and the 77,000 short lines |

**Compared with reality.** The commonly cited figures for Poland are about 8,500 km of flood levees protecting about 1.09–1.25 million ha (10,900–12,500 km²). The result, 825 km², is about **7 %** of that. On the map, protected areas appear as small spots along the levees of the Odra, Vistula, Warta and Bug. Long continuous protected floodplains, such as Żuławy or the Odra valley below Wrocław, are absent.

## 2. Visible problems
1. **Levee input is OSM, mostly not flood levees.**
   - 95 % of the lines are generic `man_made=embankment`: road, rail and field embankments, trenches and earth walls.
   - Flood levees should be tagged `man_made=dyke` in OSM (3,502 lines). Their length is not in the report.
   - The input is very fragmented: median line 125 m, and 41 % of the barrier segments are bridging connectors.
2. **Most barrier length lies far above the rivers.** Half of the segments have their crest 14 m or more above the nearest river, and a quarter 34 m or more. They can never hold back a river flood. They only cost crest reading time (24 % of the run).
3. **Implausible protection levels.** `protection_level_m` reaches 16–27 m. These are embankments across valleys (roads, dams, `waterway=dam`) that act like dams, not levees.
4. **Many tiny compartments.** The median of 4 ha sits close to `min_area_m2` (1 ha). They are typical for short embankments that enclose small depressions.
5. **River network.** MERIT lines are derived from a roughly 90 m grid and are offset from the real channels. Every MERIT reach with 25 km² or more of upstream area is a water source, because `rivers_sources.gpkg` is missing.
6. **Configuration.**
   - `river_upstream_area_field = "50"` is not a column name, so `rp_index.py` cannot run.
   - `bytes_per_cell` is still at 180 (114 measured).
   - The DTM has not been converted (strip layout).
7. **Missing steps.** No `download_osm` or `osm_hydraulics` outputs, so there are no closures or pumping stations.

## 3. Why protected areas are few and narrow (assessment, most likely first)
1. **Rivers open the levees and flood the polders from inside.**
   - With `gated_crossings = false`, every cell where a river line crosses a barrier stops being a barrier ("river wins", §5.2.6). Every river cell is also a water source at level 0.
   - **Effect of MERIT geometry:** the coarse lines cross levee lines at places where the real river does not, and they run on the landward side for stretches. Each such crossing is a hole through which the protected side fills, and every source cell inside a polder floods it from within.
   - **Effect of the 25 km² threshold:** tributaries and parallel drains behind the levees count as open water sources too, while in reality they pass the levee through sluices or flap gates or are pumped. A protected area then ends at the first internal stream or drain. This matches the narrow strips along the levees on the map and the 894 enclosed zones.
2. **Leaks through gaps.**
   - Gaps longer than `gap_tolerance_m` (20 m) stay open. OSM levee lines are often interrupted at road and rail crossings, sluices and bridges.
   - One open gap removes the protection of everything behind it, because the spill level behind becomes the gap level. The heavy fragmentation of the input (41 % bridging connectors) makes such gaps likely along many levees.
3. **Crests underestimated on a 10 m DTM.**
   - A levee crest is 3–6 m wide. A 10 m DTM (resampled, probably by averaging) lowers it, often by decimetres to over a metre.
   - The margin median is 0.7 m with a 0.3 m threshold, so a lowered crest removes area directly, and makes gaps at crossings lower still.
4. **Missing barriers.** Road and railway embankments in the floodplain that close a polder laterally are absent. Where the levee ties into such an embankment instead of high ground, the polder stays open at its end. This is what Stage G targets.
5. **Not the cause.** The corridor buffer and the unit size: only 2 compartments touch the corridor edge and 21 are truncated (0.6 %).

The ranking of 1 and 2 against 3 and 4 is an assessment from the statistics and the map; the pilot units of Stage G will measure it (see §4).

## 4. What this means for Stage G
- **Transport embankments alone will not fix the area.**
  - While rivers open every barrier they cross and every MERIT reach is a source, a new barrier class is cut by the same crossings. The ~95 % gap to reality would shrink only partly.
  - The opening logic of Stage G (bridges, culverts, switchable) must therefore apply to river crossings of all barriers, and the source and crossing rule for small tributaries must be decided first (question 3).
  - Otherwise the comparison "with and without transport embankments" measures the leaks rather than the embankments.
- **The current run already treats transport embankments as levees, without any check.**
  - Most of the 73,906 `man_made=embankment` lines are road, rail or field embankments.
  - Stage G should separate the classes:
    - levees from the levee input (`dyke`, or a national dataset);
    - transport embankments only where the DTM confirms them, with their openings.
  - Using OSM tags to classify would contradict the specification, so the levee input itself has to be cleaner.
- **DTM resolution.**
  - **Embankment widths:** road crests are 7–15 m wide, railway crests 5–8 m, footprints 15–40 m.
  - **What 10 m resolves:** crests are averaged out, side slopes are not resolved, and cuttings or one-sided slopes are hard to tell from embankments. Only large embankments (about 1.5 m high or more, footprint 30 m or more) are classifiable, and openings narrower than about 20 m are invisible.
  - The analysis needs about 1–2 m.
  - **Test data:** two open 1 m lidar terrain models are reachable from this session: USGS 3DEP (public domain) and LINZ New Zealand (CC BY 4.0). Polish, Czech, Dutch and German services are blocked by the network policy.
- **DTM storage.** Cross-sections read many small native windows. The strip-layout file reads at 19 MB/s per window, so conversion (option A, about 18 GB, a few minutes) should precede Stage G on Poland.
- **Speed.** Barriers whose crest lies above `max_level_m` over the nearest river can be skipped before crest reading, which saves most of the 24 % crest read time. This applies equally to transport embankments.

## 5. Questions
1. **Levees.** You mentioned a local dataset, but the run used the OSM Europe file (95 % `embankment`). Is there a national levee dataset (e.g. from Wody Polskie or ISOK) for Poland? Otherwise: only `man_made=dyke` as levees, and everything else left to the DTM analysis of Stage G?
2. **Rivers.** MERIT lines are offset from the channels. Options:
   - snap them to the thalweg (`snap_rivers.py`, Stage E1, with `use_snapped_rivers = true`);
   - use a better network (MPHP);
   - keep MERIT.
3. **Water sources and crossings.** I propose:
   - only rivers above an upstream area threshold (a new key, e.g. 500 km², your choice) are open water sources;
   - smaller streams that cross a barrier count as closed (sluice or flap gate) unless the DTM shows an opening; Stage G applies the same opening detection to them.

   Do you agree, and what threshold?
4. **Reference for "reality".** What do you compare against? For example ISOK maps of areas at risk when a levee fails, or your own model. With such a layer, `validate.py` measures recall per variant.
5. **DTM.** How was the 10 m DTM made (from NMT 1 m by averaging)? Can you use NMT 1 m for the pilot units? For my internal tests: is a USGS 3DEP or LINZ 1 m area acceptable? I would pick a reach with levees and road and rail embankments in the floodplain.
6. **Pilot units.** Proposal:
   - Odra around Wrocław (units 020_028, 020_029 with the most compartments);
   - Vistula below Sandomierz;
   - one unit in Żuławy (missing today).

   Or your choice.
7. **Hints from roads and railways.** OSM geometry only as a search hint. Two ways to get it:
   - your `osm_pbf` through the existing GDAL OSM driver: no new dependency, works on your machine;
   - for my tests, Overture Maps transport (reachable here, GeoParquet): this needs the GDAL Parquet plugin `libgdal-arrow-parquet`, a new dependency that needs your approval.

   Alternatively, detection from the DTM alone, without hints.
8. **Renaming** "compartment" → "protected area":
   - files: `compartments.gpkg` → `protected_areas.gpkg`, script `compartments.py` → `protected_areas.py`, pixi task `compartments` → `protected-areas`;
   - fields: `compartment_id` → `protected_area_id`;
   - config table: `[…compartments]` → `[…protected_areas]`.

   Old configuration keys would be read with a warning, as in Stage F. Agreed?
9. **`osm_hydraulics`.** Its log is empty. Did `download_osm` run? If closures and pumping stations from OSM are wanted, they need the OSM extract first.

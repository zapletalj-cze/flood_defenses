# Data report POL

Created 2026-11-06 19:51; configuration config/defenses.toml; output folder D:\01_Projects\997_SmallProjects\99_DEFENSES\01_MD\outputs\POL

## Configuration

- general.n_workers = 16
- general.max_ram_gb_per_worker = 32.0
- country.dtm = ['D:/01_Projects/997_SmallProjects/99_DEFENSES/02_DATA/Poland/PL_DEM_10m_borders_20220913.tif']
- country.levees = ['D:/01_Projects/997_SmallProjects/99_DEFENSES/02_DATA/OSM/levees/europe-latest_man_made_embankment_or_dyke_lines.gpkg']
- country.rivers = ['R:/01_Projects/900_FL_Europe/01_MD/01_HAZARD/02_RN/MERIT/MERIT_3035/MERIT_RN_with_IFLakeType_260728_addFrom_To_point_IFAddArea.gpkg']
- country.river_upstream_area_field = '50'
- country.upstream_area_raster = ''
- country.osm_pbf = 'D:/01_Projects/900_FL_Europe/01_MD/01_HAZARD/01_DTM/GEDTM/europe-latest.osm.pbf'
- country.geofabrik_region = ''
- country.working_resolution_m = 5.0
- country.corridor_buffer_m = 3000.0
- country.unit_core_size_m = 20000.0
- country.unit_halo_m = 3000.0
- country.bytes_per_cell = 180.0
- country.compartments: gap_tolerance_m = 20.0, touch_tolerance_m = 0.01, self_close_min_length_factor = 3.0, crest_search_radius_m = 5.0, barrier_removal_width_m = 10.0, min_margin_m = 0.3, max_level_m = 20.0, min_area_m2 = 10000.0, gated_crossings = False, fill_nodata_holes = True, stream_burn_m = 5.0, stream_burn_buffer_m = 100.0, preserve_river_cells = True, cap_river_holes = True, max_fill_distance_m = 500.0, use_snapped_rivers = False

## Inputs

### DTM

- PL_DEM_10m_borders_20220913.tif: 83,187 x 73,848 px, 10 m, EPSG:2180, Float32, nodata -32768.0, block 83187 x 1, overviews 0, compression LZW
  - 24 random windows: nodata share 38.2%; elevation p0 23.18, p5 25.78, p25 99.57, p50 122.9, p75 168, p95 829.8, p100 1036
- dtm_footprint.gpkg: area 366,546 km2

### levees

- europe-latest_man_made_embankment_or_dyke_lines.gpkg: 300,407 features in the file, LineString, EPSG:4326, columns ['osm_id', 'name', 'highway', 'waterway', 'aerialway', 'barrier', 'man_made', 'railway', 'z_order', 'other_tags']
- inside the DTM footprint:
- features: 77,408; geometry types: {'LineString': 77408}; empty or invalid: 0 / 0
- bounds: 152,939, 127,417, 873,444, 774,677
- line length km: total 2.292e+04; per feature p0 0.0004487, p5 0.01565, p25 0.05218, p50 0.1247, p75 0.2859, p95 1.035, p100 33.59
- `osm_id` (text): 77,408 distinct; '686576284' 1, '686576285' 1, '686576286' 1, '686576287' 1, '686576289' 1, '686576288' 1, '616610173' 1, '616610174' 1
- `name` (text): 173 distinct; 'Transzeja' 99, 'Niemieckie okopy z II WŚ' 21, 'Grobla Swojczycko-Kowalska' 17, 'Transzeje' 14, 'Grobla Janowicko-Swojczycka' 12, 'Grobla Szczytnicko-Bartoszowic' 10, 'Grobla Kozanowska' 7, 'Wał Kotowice-Siedlce' 7
- `highway` (text): 15 distinct; 'track' 1131, 'path' 622, 'service' 68, 'cycleway' 66, 'motorway' 47, 'footway' 29, 'residential' 23, 'construction' 18
- `waterway` (text): 3 distinct; 'dam' 14, 'ditch' 3, 'dyke' 1
- `aerialway` (text): 0 distinct; 
- `barrier` (text): 12 distinct; 'fence' 48, 'guard_rail' 36, 'wall' 34, 'retaining_wall' 34, 'yes' 9, 'berm' 7, 'hedge' 4, 'handrail' 4
- `man_made` (text): 2 distinct; 'embankment' 73906, 'dyke' 3502
- `railway` (text): 2 distinct; 'abandoned' 13, 'razed' 1
- `z_order` (int32): null 0; p0 -20, p5 0, p25 0, p50 0, p75 0, p95 0, p100 20
- `other_tags` (text): 1,057 distinct; '"embankment"=>"yes"' 988, '"source:geometry"=>"Geopoz: Ma' 441, '"source:geometry"=>"geoportal.' 400, '"source:geometry"=>"geoportal.' 366, '"source:geometry"=>"geoportal.' 325, '"source:geometry"=>"geoportal.' 299, '"source:geometry"=>"geoportal.' 271, '"area"=>"no"' 256

### rivers

- MERIT_RN_with_IFLakeType_260728_addFrom_To_point_IFAddArea.gpkg: 164,845 features in the file, MultiLineString, EPSG:3035, columns ['COMID', 'lengthkm', 'lengthdir', 'sinuosity', 'slope', 'uparea', 'order', 'strmDrop_t', 'slope_taud', 'NextDownID', 'maxup', 'up1', 'up2', 'up3', 'up4', 'IFLakeType', 'IFLakeLenRatio', 'Hylak_from', 'Hylak_list', 'Hylak_to', 'IFAddArea', 'pfaf', 'DF_Region', 'uparea_og', 'IFAddArea_where', 'from_point', 'to_point', 'from_area', 'to_area', 'from_area_og']
- inside the DTM footprint:
- features: 8,268; geometry types: {'MultiLineString': 8268}; empty or invalid: 0 / 0
- bounds: 150,748, 105,812, 892,145, 774,823
- line length km: total 6.916e+04; per feature p0 0.05455, p5 0.8896, p25 3.23, p50 6.504, p75 11.33, p95 22.38, p100 74.18
- `COMID` (text): 8,268 distinct; '22015356' 1, '22016203' 1, '22015227' 1, '22015628' 1, '22017302' 1, '22017088' 1, '22016013' 1, '22016439' 1
- `lengthkm` (float64): null 0; p0 0.05454, p5 0.8905, p25 3.231, p50 6.504, p75 11.33, p95 22.38, p100 74.17
- `lengthdir` (float64): null 0; p0 0.0544, p5 0.7162, p25 2.436, p50 4.743, p75 7.936, p95 14.71, p100 45.7
- `sinuosity` (float64): null 0; p0 1, p5 1.097, p25 1.22, p50 1.334, p75 1.491, p95 1.887, p100 4.441
- `slope` (float64): null 0; p0 0, p5 5.027e-05, p25 0.0004887, p50 0.001114, p75 0.002394, p95 0.008088, p100 0.06196
- `uparea` (float64): null 0; p0 25.07, p5 27.75, p25 45.07, p50 106.6, p75 451, p95 1.177e+04, p100 194,369
- `order` (int64): null 0; p0 1, p5 1, p25 1, p50 1, p75 3, p95 5, p100 7
- `strmDrop_t` (float64): null 0; p0 0, p5 0, p25 1.7, p50 6.4, p75 16.5, p95 55.73, p100 480.1
- `slope_taud` (float64): null 0; p0 0, p5 0, p25 0.0003111, p50 0.0009873, p75 0.002233, p95 0.007795, p100 0.03444
- `NextDownID` (int64): null 0; p0 0, p5 24,027,796, p25 24,030,963, p50 24,033,737, p75 24,035,950, p95 24,037,810, p100 24,038,071
- `maxup` (int64): null 0; p0 0, p5 0, p25 0, p50 0, p75 2, p95 2, p100 3
- `up1` (int64): null 0; p0 0, p5 0, p25 0, p50 0, p75 24,033,758, p95 24,037,233, p100 24,038,422
- `up2` (float64): null 0; p0 0, p5 0, p25 0, p50 0, p75 24,034,541, p95 24,037,564, p100 24,038,427
- `up3` (int64): null 0; p0 0, p5 0, p25 0, p50 0, p75 0, p95 0, p100 24,034,910
- `up4` (int64): null 0; p0 0, p5 0, p25 0, p50 0, p75 0, p95 0, p100 0
- `IFLakeType` (int64): null 0; p0 0, p5 0, p25 0, p50 0, p75 0, p95 2, p100 2
- `IFLakeLenRatio` (float64): null 0; p0 0, p5 0, p25 0, p50 0, p75 0, p95 46.27, p100 100
- `Hylak_from` (text): 560 distinct; '' 7552, '13580' 11, '1186' 8, '13403' 7, '13579' 6, '164251' 5, '13456' 5, '13422' 4
- `Hylak_list` (text): 818 distinct; '' 7443, '1335398' 2, '166229' 2, '1334559' 2, '1320701' 2, '1338864' 2, '1330774' 2, '1333706' 2
- `Hylak_to` (text): 310 distinct; '' 7468, '13580' 22, '1186' 14, '164251' 10, '13579' 10, '13472' 8, '13403' 8, '13839' 6
- `IFAddArea` (float64): null 8,259; p0 -1450, p5 -1365, p25 -442, p50 58, p75 442, p95 1365, p100 1450
- `pfaf` (int64): null 0; p0 22, p5 24, p25 24, p50 24, p75 24, p95 24, p100 24
- `DF_Region` (text): 1 distinct; 'Baltic' 371
- `uparea_og` (float64): null 8,268; no finite values
- `IFAddArea_where` (text): 1 distinct; 't' 4
- `from_point` (text): 8,268 distinct; '22_15404' 1, '22_16266' 1, '22_15272' 1, '22_15679' 1, '22_17386' 1, '22_17172' 1, '22_16072' 1, '22_16502' 1
- `to_point` (text): 4,187 distinct; '24_18205' 3, '24_12291' 3, '24_15320' 3, '24_14911' 3, '24_14924' 3, '22_15404' 2, '22_15679' 2, '22_15272' 2
- `from_area` (float64): null 0; p0 0, p5 0, p25 0, p50 0, p75 413.9, p95 1.175e+04, p100 193,921
- `to_area` (float64): null 0; p0 13.74, p5 27.74, p25 45.07, p50 106.6, p75 450.3, p95 1.177e+04, p100 193,927
- `from_area_og` (text): 203 distinct; '0.0' 4, '76.50450911879119' 1, '29111.44398699139' 1, '28979.229134280125' 1, '28869.94596940129' 1, '28765.227717099682' 1, '28536.24905501525' 1, '28334.09241398078' 1
- `to_area_og` (text): 206 distinct; '222.28775081558064' 1, '82.18219551734546' 1, '29160.794257025987' 1, '28980.79788762067' 1, '28892.664531716142' 1, '28782.865599917222' 1, '28723.04052423964' 1, '28343.34773007572' 1
- `from_area_diff` (text): 14 distinct; '0' 8062, '-442.0' 111, '0.0' 68, '442.0' 7, '-152.0' 3, '-94.0' 3, '-1450.1480353929' 3, '1237.9999999999998' 3
- `to_area_diff` (text): 18 distinct; '0' 8062, '-442.0' 107, '0.0' 64, '-442.00000000001455' 5, '-1450.1480353929' 4, '1238.0000000000002' 4, '442.0' 4, '-152.0' 3

## Outputs

- 12 files, 0.08 GB (without logs and plots); top level: barriers.gpkg 48.6 MB, compartments.gpkg 21.5 MB, corridor.gpkg 3.4 MB, dtm_benchmark.csv 0.0 MB, dtm_estimates.json 0.0 MB, dtm_footprint.gpkg 0.4 MB, dtm_inspection.csv 0.0 MB, enclosed_streams.gpkg 1.2 MB, performance.json 0.0 MB, units.csv 0.1 MB, units_core.gpkg 0.3 MB, units_window.gpkg 0.3 MB

### compartments.gpkg

- features: 3,255; geometry types: {'Polygon': 3255}; empty or invalid: 0 / 0
- bounds: 164,280, 175,770, 855,510, 749,595
- polygon area km2: total 824.9; per feature p0 0.01, p5 0.0112, p25 0.01816, p50 0.0395, p75 0.1337, p95 1.159, p100 11.57
- `compartment_id` (text): 3,255 distinct; 'POL_725180_480580' 1, 'POL_725510_480210' 1, 'POL_729660_478160' 1, 'POL_705560_317580' 1, 'POL_706970_316480' 1, 'POL_705680_316060' 1, 'POL_708140_315790' 1, 'POL_706390_315880' 1
- `area_m2` (float64): null 0; p0 1e+04, p5 1.12e+04, p25 1.816e+04, p50 3.95e+04, p75 133,712, p95 1,159,420, p100 11,571,625
- `protection_level_m` (float64): null 0; p0 0.3026, p5 0.6119, p25 1.64, p50 3.173, p75 5.877, p95 16.44, p100 27.27
- `margin_median_m` (float64): null 0; p0 0.3, p5 0.3399, p25 0.484, p50 0.7, p75 1.105, p95 2.347, p100 9.202
- `margin_max_m` (float64): null 0; p0 0.3, p5 0.3693, p25 0.6686, p50 1.098, p75 1.836, p95 3.64, p100 10.7
- `river_ref_median_m` (float64): null 0; p0 -2.323, p5 17.06, p25 73.2, p50 120.9, p75 170.2, p95 250.9, p100 584.8
- `entry_x` (float64): null 0; p0 165,712, p5 214,222, p25 340,292, p50 482,122, p75 621,615, p95 707,729, p100 855,442
- `entry_y` (float64): null 0; p0 177,122, p5 238,255, p25 295,870, p50 397,908, p75 501,368, p95 684,924, p100 749,592
- `source_x` (float64): null 0; p0 164,908, p5 212,813, p25 340,900, p50 482,442, p75 621,442, p95 707,438, p100 855,452
- `source_y` (float64): null 0; p0 177,498, p5 238,398, p25 296,440, p50 398,428, p75 501,228, p95 685,570, p100 750,698
- `truncated` (bool): true 21
- `touches_corridor_edge` (bool): true 2
- `touches_dtm_edge` (bool): true 2
- `has_pumping` (bool): true 0
- `no_levee` (bool): true 0
- `unit_id` (text): 440 distinct; '020_029' 95, '013_023' 89, '027_027' 65, '020_028' 62, '022_029' 49, '021_029' 49, '012_023' 45, '023_029' 45
- `unit_label` (int64): null 0; p0 1, p5 1, p25 5, p50 12, p75 26, p95 56.3, p100 110
- largest 10:

| compartment_id | area_m2 | protection_level_m | margin_median_m | truncated | unit_id |
|---|---|---|---|---|---|
| POL_700200_536240 | 11,571,625 | 2.758 | 0.7237 | False | 029_014 |
| POL_253770_544160 | 10,761,875 | 3.691 | 1.018 | False | 007_014 |
| POL_613100_506020 | 10,656,475 | 5.056 | 1.566 | False | 025_016 |
| POL_547460_516780 | 9,629,950 | 3.415 | 1.858 | False | 022_015 |
| POL_695220_538800 | 9,208,775 | 1.729 | 1.142 | False | 029_014 |
| POL_209380_658760 | 9,201,050 | 0.3557 | 0.3525 | False | 005_008 |
| POL_204260_528860 | 9,103,200 | 1.915 | 1.365 | False | 005_015 |
| POL_270200_466380 | 8,938,150 | 1.674 | 0.8008 | False | 008_018 |
| POL_204870_623670 | 7,814,450 | 0.9395 | 0.8703 | False | 005_010 |
| POL_645240_516660 | 7,789,350 | 0.3668 | 0.3668 | False | 027_015 |


### barriers.gpkg

- features: 151,273; geometry types: {'LineString': 151273}; empty or invalid: 0 / 0
- bounds: 8.9e+04, 8.874e+04, 963,581, 790,439
- line length km: total 2.605e+04; per feature p0 1.111e-05, p5 0.003295, p25 0.01081, p50 0.03129, p75 0.1521, p95 0.6808, p100 33.59
- `line_id` (int64): null 0; p0 0, p5 7564, p25 3.782e+04, p50 7.564e+04, p75 113,454, p95 143,708, p100 151,272
- `source` (text): 1 distinct; 'europe-latest_man_made_embankm' 151273
- `bridged` (bool): true 62,141
- `crest_min_m` (float64): null 13,909; p0 -1.988, p5 29.63, p25 92.88, p50 153.4, p75 254.2, p95 543.8, p100 2180
- `crest_hand_min_m` (float64): null 19,619; p0 -17.23, p5 0.644, p25 4.749, p50 13.94, p75 34.35, p95 108.2, p100 1335

### enclosed_streams.gpkg

- features: 894; geometry types: {'Polygon': 894}; empty or invalid: 0 / 0
- bounds: 153,210, 147,610, 825,615, 754,475
- polygon area km2: total 55.67; per feature p0 2.5e-05, p5 2.5e-05, p25 2.5e-05, p50 0.000125, p75 0.001813, p95 0.02108, p100 28.44
- `zone` (int64): null 0; p0 1, p5 1, p25 2, p50 5, p75 17, p95 58, p100 80
- `area_m2` (float64): null 0; p0 25, p5 25, p25 25, p50 125, p75 1812, p95 2.108e+04, p100 28,444,625
- `stream_cells` (int64): null 0; p0 1, p5 1, p25 1, p50 1, p75 3, p95 20.35, p100 240
- `gated` (bool): true 0
- `unit_id` (text): 230 distinct; '020_028' 72, '012_024' 51, '032_022' 38, '029_026' 35, '021_007' 32, '005_009' 30, '011_024' 25, '020_029' 24

### units_core.gpkg

- features: 971; geometry types: {'Polygon': 971}; empty or invalid: 0 / 0
- bounds: 141,950, 109,700, 901,950, 789,700
- polygon area km2: total 388,400; per feature p0 400, p5 400, p25 400, p50 400, p75 400, p95 400, p100 400
- `unit_id` (text): 971 distinct; '036_023' 1, '031_017' 1, '030_019' 1, '030_025' 1, '030_026' 1, '030_027' 1, '030_029' 1, '005_016' 1
- `est_window_cells` (int64): null 0; p0 6,760,000, p5 6,760,000, p25 6,760,000, p50 6,760,000, p75 6,760,000, p95 6,760,000, p100 6,760,000
- `est_domain_cells` (int64): null 0; p0 563,028, p5 2,652,206, p25 6,599,304, p50 6,731,467, p75 6,759,826, p95 6,760,000, p100 6,760,000
- `est_ram_gb` (float64): null 0; p0 1.217, p5 1.217, p25 1.217, p50 1.217, p75 1.217, p95 1.217, p100 1.217

### units.csv

- rows 971
- `unit_id`: '036_023' 1, '031_017' 1, '030_019' 1, '030_025' 1, '030_026' 1, '030_027' 1, '030_029' 1, '005_016' 1
- `core_xmin`: p0 141,950, p5 201,950, p25 361,950, p50 521,950, p75 681,950, p95 801,950, p100 881,950
- `core_xmax`: p0 161,950, p5 221,950, p25 381,950, p50 541,950, p75 701,950, p95 821,950, p100 901,950
- `core_ymin`: p0 109,700, p5 169,700, p25 309,700, p50 449,700, p75 589,700, p95 709,700, p100 769,700
- `core_ymax`: p0 129,700, p5 189,700, p25 329,700, p50 469,700, p75 609,700, p95 729,700, p100 789,700
- `win_xmin`: p0 138,950, p5 198,950, p25 358,950, p50 518,950, p75 678,950, p95 798,950, p100 878,950
- `win_xmax`: p0 164,950, p5 224,950, p25 384,950, p50 544,950, p75 704,950, p95 824,950, p100 904,950
- `win_ymin`: p0 106,700, p5 166,700, p25 306,700, p50 446,700, p75 586,700, p95 706,700, p100 766,700
- `win_ymax`: p0 132,700, p5 192,700, p25 332,700, p50 472,700, p75 612,700, p95 732,700, p100 792,700
- `win_col0`: p0 3700, p5 9700, p25 2.57e+04, p50 4.17e+04, p75 5.77e+04, p95 6.97e+04, p100 7.77e+04
- `win_row0`: p0 3700, p5 9700, p25 2.17e+04, p50 3.57e+04, p75 4.97e+04, p95 6.37e+04, p100 6.97e+04
- `win_ncols`: p0 2600, p5 2600, p25 2600, p50 2600, p75 2600, p95 2600, p100 2600
- `win_nrows`: p0 2600, p5 2600, p25 2600, p50 2600, p75 2600, p95 2600, p100 2600
- `native_col0`: p0 3700, p5 9700, p25 2.57e+04, p50 4.17e+04, p75 5.77e+04, p95 6.97e+04, p100 7.77e+04
- `native_row0`: p0 3700, p5 9700, p25 2.17e+04, p50 3.57e+04, p75 4.97e+04, p95 6.37e+04, p100 6.97e+04
- `native_ncols`: p0 2600, p5 2600, p25 2600, p50 2600, p75 2600, p95 2600, p100 2600
- `native_nrows`: p0 2600, p5 2600, p25 2600, p50 2600, p75 2600, p95 2600, p100 2600
- `est_window_cells`: p0 6,760,000, p5 6,760,000, p25 6,760,000, p50 6,760,000, p75 6,760,000, p95 6,760,000, p100 6,760,000
- `est_domain_cells`: p0 563,028, p5 2,652,206, p25 6,599,304, p50 6,731,467, p75 6,759,826, p95 6,760,000, p100 6,760,000
- `est_ram_gb`: p0 1.217, p5 1.217, p25 1.217, p50 1.217, p75 1.217, p95 1.217, p100 1.217

### performance.json

```
{"units_measured": 971, "seconds_per_domain_cell": 3.508324787485751e-06, "bytes_per_window_cell": 113.8094201183432, "country_seconds": 1335.6303230349904, "remaining_seconds": 0.0, "time_shares": {"read_terrain": 0.07683605852280367, "domain": 0.038251703824282064, "crest_read": 0.23891773752961712, "barriers": 0.02601866298525148, "rivers": 0.02491610846650123, "undefended": 0.05766008145695036, "hand": 0.2268474718396896, "flood": 0.2510903240409032, "compartments": 0.016976167140285534, "barrier_crest": 0.0013900050210074724, "write": 0.04107306565588229}}
```

### dtm_estimates.json

```
{"country": "POL", "dtm": [{"path": "D:/01_Projects/997_SmallProjects/99_DEFENSES/02_DATA/Poland/PL_DEM_10m_borders_20220913.tif", "windows": 20, "native_s_median": 3.5473704500000167, "native_s_p90": 3.9629959100006094, "native_mb_s_median": 18.91933661693833, "warp_s_median": 3.9982572000008076, "warp_s_p90": 4.402466130003449, "warp_native_equivalent_mb_s_median": 16.78500379255366, "window_ground_km2": 1677.7216, "sequential_mb_s": 543.0957870084188, "compression_ratio": 0.5508424043655396, "compression_mb_s": 309.41085994746, "raw_gb": 24.572774304, "window_px": 4096, "option_a_disk_gb": 18.047634772729474, "option_a_read_s": 45.24574649963006, "option_a_compress_s": 79.4179438568272}], "dsm": []}
```

## Logs (newest per script)

### build_units (build_units_20261106_135245.log, 9 lines)

- 1 x WARNING Placeholder value in use: country.bytes_per_cell = #
- last lines:
```
2026-11-06 14:10:57,845 INFO build_units: Units: 971 of 1554 cores intersect the corridor (core 20000 m, halo 3000 m, working cell 10.0 m).
2026-11-06 14:11:57,916 INFO build_units: Largest window 2600 x 2600 cells needs about 1.2 GB (bytes_per_cell 180); 16 workers need about 19.5 GB of RAM in total. Domain cells to process: 6091.3 million.
2026-11-06 14:11:58,148 INFO pyogrio._io: Created 971 records
2026-11-06 14:11:58,377 INFO pyogrio._io: Created 971 records
2026-11-06 14:12:01,322 INFO pyogrio._io: Created 971 records
2026-11-06 14:12:01,428 INFO build_units: Written corridor.gpkg, units_core.gpkg, units_window.gpkg and units.csv to D:\01_Projects\997_SmallProjects\99_DEFENSES\01_MD\outputs\POL.
```

### compartments (compartments_20261106_142457.log, 985 lines)

- 1 x WARNING Placeholder value in use: country.bytes_per_cell = #
- 1 x WARNING # levee-enclosed areas contain river cells (see enclosed_streams.gpkg).
- last lines:
```
2026-11-06 16:24:14,557 INFO compartments: Compartments: 3255 from 5566 unit candidates in 971 finished units; truncated 21 (0.6 %; many truncated compartments mean the halo is too small), dropped as contradicted by neighbouring windows 105, touching the corridor edge 2, no levee 0, with pumping 0. 
2026-11-06 16:24:17,466 INFO pyogrio._io: Created 151,273 records
2026-11-06 16:24:25,041 INFO pyogrio._io: Created 894 records
2026-11-06 16:24:25,099 WARNING compartments: 894 levee-enclosed areas contain river cells (see enclosed_streams.gpkg).
2026-11-06 16:24:25,385 INFO compartments: Measured on 971 units: 3.51 microseconds per domain cell, up to 114 bytes per window cell (configured 180), peak RSS up to 2.87 GB. Time shares: flood 25 %, crest_read 24 %, hand 23 %, read_terrain 8 %, undefended 6 %, write 4 %, domain 4 %, barriers 3 %, r
2026-11-06 16:24:25,503 INFO compartments: Extrapolation (measured s/cell x estimated domain cells / 16 workers): whole country 0:22:16, remaining 0 units 0:00:00; largest unit window needs about 0.8 GB.
```

### osm_hydraulics (osm_hydraulics_20261106_135212.log, 0 lines)

- last lines:
```
```

### prepare_dtm (prepare_dtm_20261106_134745.log, 12 lines)

- 1 x WARNING PL_DEM_#m_borders_#.tif: strip layout; window reads at working resolution will be slow (see the benchmark and dtm_conversion).
- last lines:
```
2026-11-06 13:50:27,226 INFO prepare_dtm: Benchmark DTM PL_DEM_10m_borders_20220913.tif (20 windows of 4096 x 4096 px = 1677.72 km2): native read 3.55 s/window (p90 3.96), 19 MB/s; warp to 10.0 m 4.00 s/window (p90 4.40), 17 MB/s native equivalent; sequential read 543 MB/s; ZSTD compression ratio 0.
2026-11-06 13:50:27,229 INFO prepare_dtm: Note: the operating system file cache can make repeated benchmarks faster than a cold read.
2026-11-06 13:50:27,229 INFO prepare_dtm: Terrain read estimate for one unit window of 26.0 x 26.0 km: 2 s (median warp s/window x area ratio). Run build_units.py for the per-country total.
2026-11-06 13:50:27,230 INFO prepare_dtm: Option A for PL_DEM_10m_borders_20220913.tif: about 18.0 GB on disk (raw 24.6 GB x measured ratio 0.55 + 1/3 overviews); one sequential pass 0:00:45 (read at 543 MB/s) to 0:02:05 (plus compression at 309 MB/s); overview building not included.
2026-11-06 13:50:27,230 INFO prepare_dtm: Option B estimate needs corridor.gpkg; run build_units.py first.
2026-11-06 13:50:27,246 INFO prepare_dtm: Processing DTM: D:\01_Projects\997_SmallProjects\99_DEFENSES\01_MD\work\POL\dtm\dtm.vrt -> original files
```

Map: overview.png (364 KB)

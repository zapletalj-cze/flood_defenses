# Data report POL

Created 2026-11-13 10:07; configuration D:\90_PersonalFoldlers\JZa\00_GIT\flood_defenses\config\defenses.toml; output folder D:\01_Projects\997_SmallProjects\99_DEFENSES\01_MD\outputs\POL

## Configuration

- general.n_workers = 16
- general.max_ram_gb_per_worker = 32.0
- country.dtm = ['D:/01_Projects/997_SmallProjects/99_DEFENSES/02_DATA/Poland/PL_DEM_10m_borders_20220913.tif']
- country.levees = ['D:/01_Projects/997_SmallProjects/99_DEFENSES/02_DATA/OSM/levees/europe-latest_man_made_embankment_or_dyke_lines.gpkg']
- country.rivers = ['R:/01_Projects/900_FL_Europe/01_MD/01_HAZARD/02_RN/MERIT/MERIT_3035/MERIT_RN_with_IFLakeType_260728_addFrom_To_point_IFAddArea.gpkg']
- country.river_upstream_area_field = 'uparea'
- country.upstream_area_raster = ''
- country.river_guides = ['Q:/01_Projects/165_Global_Flood_maps/01_MD/01_Europe/01_GIS/02_OSM/Europe_3035_80_waterway.gpkg']
- country.working_resolution_m = 10.0
- country.corridor_buffer_m = 3000.0
- country.unit_core_size_m = 20000.0
- country.unit_halo_m = 3000.0
- country.bytes_per_cell = 180.0
- country.protected_areas: gap_tolerance_m = 20.0, touch_tolerance_m = 0.01, self_close_min_length_factor = 3.0, crest_search_radius_m = 5.0, barrier_removal_width_m = 10.0, min_margin_m = 0.3, max_level_m = 20.0, min_area_m2 = 10000.0, gated_crossings = False, fill_nodata_holes = True, stream_burn_m = 5.0, stream_burn_buffer_m = 100.0, preserve_river_cells = True, cap_river_holes = True, max_fill_distance_m = 500.0, use_snapped_rivers = False

## Inputs

### DTM

- PL_DEM_10m_borders_20220913.tif: 83,187 x 73,848 px, 10 m, EPSG:2180, Float32, nodata -32768.0, block 83187 x 1, overviews 0, compression LZW
  - 24 random windows: nodata share 38.2%; elevation p0 23.18, p5 25.78, p25 99.57, p50 122.9, p75 168, p95 829.8, p100 1036

### levees

- europe-latest_man_made_embankment_or_dyke_lines.gpkg: 300,407 features in the file, LineString, EPSG:4326, columns ['osm_id', 'name', 'highway', 'waterway', 'aerialway', 'barrier', 'man_made', 'railway', 'z_order', 'other_tags']
- inside the DTM footprint (no footprint yet: whole file):
- features: 300,407; geometry types: {'LineString': 300407}; empty or invalid: 0 / 0
- bounds: -3,851,016, -1,406,909, 2,743,259, 3,433,607
- line length km: total 7.562e+04; per feature p0 1.125e-05, p5 0.01357, p25 0.04151, p50 0.1005, p75 0.2342, p95 0.8729, p100 54.76
- `osm_id` (text): 300,407 distinct; '3499086' 1, '3500100' 1, '3500101' 1, '3501443' 1, '3514543' 1, '3549107' 1, '3549139' 1, '3567716' 1
- `name` (text): 1,775 distinct; 'Transzeja' 99, "Wat's Dyke" 94, "Offa's Dyke" 70, 'Rheinhauptdeich' 69, 'Val obrov' 57, 'Ligne de la Lauter' 53, 'Levée de la Loire' 51, 'Elbdeich' 50
- `highway` (text): 24 distinct; 'track' 3888, 'path' 2365, 'footway' 643, 'unclassified' 593, 'primary' 423, 'service' 407, 'secondary' 403, 'tertiary' 363
- `waterway` (text): 8 distinct; 'dam' 153, 'stream' 24, 'ditch' 19, 'drain' 13, 'weir' 12, 'dyke' 3, 'artificial' 1, 'river' 1
- `aerialway` (text): 0 distinct; 
- `barrier` (text): 32 distinct; 'retaining_wall' 601, 'wall' 595, 'fence' 516, 'hedge' 167, 'berm' 128, 'guard_rail' 90, 'ditch' 67, 'city_wall' 35
- `man_made` (text): 2 distinct; 'embankment' 280416, 'dyke' 19991
- `railway` (text): 8 distinct; 'abandoned' 112, 'rail' 93, 'razed' 11, 'disused' 2, 'miniature' 2, 'construction' 2, 'light_rail' 1, 'loading_ramp' 1
- `z_order` (int32): null 0; p0 -20, p5 0, p25 0, p50 0, p75 0, p95 0, p100 30
- `other_tags` (text): 10,422 distinct; '"embankment"=>"yes"' 4568, '"embankment"=>"yes","mml:class' 3224, '"natural"=>"earth_bank"' 2403, '"area"=>"no"' 1654, '"source:geometry"=>"Geopoz: Ma' 441, '"mml:class"=>"30100"' 432, '"material"=>"soil"' 408, '"source:geometry"=>"geoportal.' 400

### rivers

- MERIT_RN_with_IFLakeType_260728_addFrom_To_point_IFAddArea.gpkg: 164,845 features in the file, MultiLineString, EPSG:3035, columns ['COMID', 'lengthkm', 'lengthdir', 'sinuosity', 'slope', 'uparea', 'order', 'strmDrop_t', 'slope_taud', 'NextDownID', 'maxup', 'up1', 'up2', 'up3', 'up4', 'IFLakeType', 'IFLakeLenRatio', 'Hylak_from', 'Hylak_list', 'Hylak_to', 'IFAddArea', 'pfaf', 'DF_Region', 'uparea_og', 'IFAddArea_where', 'from_point', 'to_point', 'from_area', 'to_area', 'from_area_og']
- inside the DTM footprint (no footprint yet: whole file):
- features: 164,845; geometry types: {'MultiLineString': 164845}; empty or invalid: 0 / 0
- bounds: -1,999,586, -1,996,061, 2,543,261, 3,636,643
- line length km: total 1,424,506; per feature p0 0.001113, p5 0.8589, p25 3.25, p50 6.587, p75 11.76, p95 23.57, p100 112.3
- `COMID` (text): 164,845 distinct; '21000001' 1, '21000002' 1, '21000003' 1, '21000004' 1, '21000005' 1, '21000006' 1, '21000007' 1, '21000008' 1
- `lengthkm` (float64): null 0; p0 0.001094, p5 0.8486, p25 3.209, p50 6.51, p75 11.61, p95 23.28, p100 111.9
- `lengthdir` (float64): null 0; p0 0.001111, p5 0.7015, p25 2.453, p50 4.744, p75 8.118, p95 15.42, p100 80.21
- `sinuosity` (float64): null 0; p0 1, p5 1.081, p25 1.196, p50 1.315, p75 1.492, p95 1.922, p100 7.658
- `slope` (float64): null 0; p0 0, p5 4.25e-05, p25 0.0006676, p50 0.002123, p75 0.006243, p95 0.02802, p100 0.9134
- `uparea` (float64): null 0; p0 25.03, p5 27.65, p25 44.36, p50 101.4, p75 433.9, p95 8649, p100 795,594
- `order` (int64): null 0; p0 1, p5 1, p25 1, p50 1, p75 2, p95 4, p100 8
- `strmDrop_t` (float64): null 0; p0 0, p5 0, p25 2.5, p50 11.9, p75 40.4, p95 215.1, p100 2206
- `slope_taud` (float64): null 0; p0 0, p5 0, p25 0.000467, p50 0.001901, p75 0.005933, p95 0.02743, p100 0.4357
- `NextDownID` (int64): null 0; p0 0, p5 21,001,651, p25 22,009,704, p50 22,051,537, p75 24,006,871, p95 25,002,059, p100 25,011,633
- `maxup` (int64): null 0; p0 0, p5 0, p25 0, p50 0, p75 2, p95 2, p100 4
- `up1` (int64): null 0; p0 0, p5 0, p25 0, p50 0, p75 22,051,652, p95 24,031,783, p100 25,011,956
- `up2` (float64): null 7; p0 0, p5 0, p25 0, p50 0, p75 22,052,454, p95 24,032,865, p100 25,011,972
- `up3` (int64): null 0; p0 0, p5 0, p25 0, p50 0, p75 0, p95 0, p100 25,010,834
- `up4` (int64): null 0; p0 0, p5 0, p25 0, p50 0, p75 0, p95 0, p100 24,015,899
- `IFLakeType` (int64): null 0; p0 0, p5 0, p25 0, p50 0, p75 0, p95 2, p100 3
- `IFLakeLenRatio` (float64): null 0; p0 0, p5 0, p25 0, p50 0, p75 0, p95 95.91, p100 100
- `Hylak_from` (text): 10,770 distinct; '' 145789, '100' 272, '10' 265, '99' 247, '105' 150, '97' 94, '103' 92, '128' 82
- `Hylak_list` (text): 13,919 distinct; '' 150874, '1335260' 3, '1351099' 3, '124146' 3, '134274' 3, '1362109' 2, '1367638' 2, '14801' 2
- `Hylak_to` (text): 5,771 distinct; '' 141927, '100' 380, '10' 315, '99' 301, '105' 191, '128' 134, '103' 133, '97' 128
- `IFAddArea` (float64): null 164,785; p0 -1.325e+04, p5 -5114, p25 -449.2, p50 104.3, p75 690.7, p95 5910, p100 3.663e+04
- `pfaf` (int64): null 0; p0 21, p5 21, p25 22, p50 23, p75 24, p95 25, p100 25
- `DF_Region` (text): 3 distinct; 'Nordics' 22563, 'Baltic' 7556, 'Denmark' 762
- `uparea_og` (float64): null 164,832; p0 30.19, p5 61.48, p25 90.76, p50 244.8, p75 380.9, p95 525.1, p100 593.4
- `IFAddArea_where` (text): 2 distinct; 't' 39, 'f' 1
- `from_point` (text): 164,845 distinct; '21_1' 1, '21_3' 1, '21_5' 1, '21_7' 1, '21_8' 1, '21_9' 1, '21_10' 1, '21_11' 1
- `to_point` (text): 85,375 distinct; '24_1622' 4, '21_1404' 3, '21_3525' 3, '21_3710' 3, '21_5677' 3, '21_8302' 3, '21_8518' 3, '21_9312' 3
- `from_area` (float64): null 0; p0 0, p5 0, p25 0, p50 0, p75 397.3, p95 8622, p100 795,535
- `to_area` (float64): null 0; p0 13.74, p5 27.65, p25 44.36, p50 101.5, p75 434, p95 8650, p100 795,594
- `from_area_og` (text): 1,414 distinct; '0.0' 18, '523.458013375317' 1, '483.1019822677213' 1, '343.37265468202355' 1, '303.79410821986926' 1, '267.06724694502736' 1, '209.64533180338208' 1, '158.36319186350102' 1
- `to_area_og` (text): 1,430 distinct; '70.03' 2, '535.7969170529597' 1, '486.7530268684776' 1, '359.1395054408034' 1, '310.20026362051107' 1, '271.5144655038161' 1, '218.37192689735215' 1, '166.07389825767902' 1
- `from_area_diff` (text): 107 distinct; '0' 163414, '0.0' 680, '-442.0' 111, '470.65453532796164' 54, '442.0' 35, '-5674.787433206671' 28, '-166.92143680917798' 25, '114.6273461488945' 24
- `to_area_diff` (text): 130 distinct; '0' 163414, '0.0' 643, '-442.0' 107, '470.65453532796164' 53, '442.0' 32, '-5674.787433206671' 30, '-166.92143680917889' 25, '-166.92143680917798' 25

## Outputs

- 18 files, 0.45 GB (without logs and plots); top level: barriers.gpkg 90.7 MB, corridor.gpkg 3.4 MB, dtm_benchmark.csv 0.0 MB, dtm_estimates.json 0.0 MB, dtm_inspection.csv 0.0 MB, enclosed_streams.gpkg 0.6 MB, leaks.gpkg 17.4 MB, levees_crest.gpkg 64.5 MB, performance.json 0.0 MB, protected_areas.gpkg 20.9 MB, protected_areas_by_group.gpkg 188.7 MB, river_groups.gpkg 31.5 MB, rivers_channel.gpkg 34.7 MB, units.csv 0.1 MB, units_core.gpkg 0.3 MB, units_window.gpkg 0.3 MB, worldcover.vrt 0.0 MB, worldcover_sources.json 0.0 MB

### protected_areas.gpkg

- features: 3,760; geometry types: {'Polygon': 3760}; empty or invalid: 0 / 0
- bounds: 166,630, 175,890, 855,500, 749,670
- polygon area km2: total 887.5; per feature p0 0.01, p5 0.0112, p25 0.0182, p50 0.0414, p75 0.1371, p95 1.047, p100 15.47
- `protected_area_id` (text): 3,760 distinct; 'POL_725510_480200' 1, 'POL_705620_317560' 1, 'POL_705220_316740' 1, 'POL_706970_316480' 1, 'POL_705280_316500' 1, 'POL_707860_315900' 1, 'POL_711440_315680' 1, 'POL_707000_315540' 1
- `area_m2` (float64): null 0; p0 1e+04, p5 1.12e+04, p25 1.82e+04, p50 4.14e+04, p75 137,100, p95 1,046,510, p100 15,467,900
- `protection_level_m` (float64): null 0; p0 0.3001, p5 0.7928, p25 2.017, p50 3.598, p75 6.428, p95 17.12, p100 27.22
- `margin_median_m` (float64): null 0; p0 0.3001, p5 0.344, p25 0.4852, p50 0.6813, p75 1.083, p95 2.534, p100 9.897
- `margin_max_m` (float64): null 0; p0 0.3001, p5 0.3791, p25 0.6844, p50 1.098, p75 1.792, p95 3.898, p100 11.67
- `river_ref_median_m` (float64): null 0; p0 -1.786, p5 14.02, p25 71.85, p50 118.9, p75 168.3, p95 249.6, p100 583
- `entry_x` (float64): null 0; p0 166,915, p5 212,998, p25 342,528, p50 483,635, p75 622,532, p95 705,698, p100 855,445
- `entry_y` (float64): null 0; p0 176,135, p5 240,174, p25 296,810, p50 403,370, p75 507,770, p95 688,976, p100 749,605
- `source_x` (float64): null 0; p0 164,915, p5 211,884, p25 342,095, p50 483,135, p75 622,895, p95 705,548, p100 855,445
- `source_y` (float64): null 0; p0 175,715, p5 239,685, p25 297,765, p50 403,595, p75 507,808, p95 688,916, p100 750,685
- `truncated` (bool): true 26
- `touches_corridor_edge` (bool): true 3
- `touches_dtm_edge` (bool): true 2
- `no_levee` (bool): true 0
- `unit_id` (text): 461 distinct; '020_028' 97, '013_023' 78, '020_029' 75, '027_027' 68, '022_029' 57, '029_025' 56, '029_026' 53, '023_029' 51
- `unit_label` (int64): null 0; p0 1, p5 1, p25 5, p50 13, p75 28, p95 61, p100 131
- largest 10:

| protected_area_id | area_m2 | protection_level_m | margin_median_m | truncated | unit_id |
|---|---|---|---|---|---|
| POL_172770_575520 | 15,467,900 | 2.36 | 1.437 | False | 003_012 |
| POL_253560_544200 | 10,456,500 | 3.777 | 0.9382 | False | 007_014 |
| POL_613090_505960 | 10,330,100 | 4.526 | 1.569 | False | 025_016 |
| POL_693360_409100 | 10,036,500 | 4.009 | 0.5694 | False | 029_021 |
| POL_547460_516780 | 9,632,500 | 3.417 | 1.872 | False | 022_015 |
| POL_662420_441820 | 7,912,200 | 4.724 | 0.8333 | False | 027_019 |
| POL_204870_623660 | 7,726,800 | 0.8402 | 0.8356 | False | 005_010 |
| POL_645240_516660 | 7,521,200 | 0.3433 | 0.3433 | False | 027_015 |
| POL_696100_368680 | 7,393,900 | 3.943 | 0.9337 | False | 029_023 |
| POL_687380_300900 | 7,008,000 | 5.972 | 2.33 | False | 029_026 |


### protected_areas_by_group.gpkg

- features: 21,769; geometry types: {'Polygon': 21769}; empty or invalid: 0 / 0
- bounds: 156,650, 158,000, 855,500, 754,500
- polygon area km2: total 1.569e+04; per feature p0 0.01, p5 0.0114, p25 0.0196, p50 0.0497, p75 0.2018, p95 2.285, p100 421.1
- `protected_area_id` (text): 21,769 distinct; 'POL_725510_480200_g0' 1, 'POL_705650_317540_g0' 1, 'POL_706830_316580_g0' 1, 'POL_707860_315900_g0' 1, 'POL_711440_315680_g0' 1, 'POL_707000_315540_g0' 1, 'POL_711880_314940_g0' 1, 'POL_709670_314980_g0' 1
- `area_m2` (float64): null 0; p0 1e+04, p5 1.14e+04, p25 1.96e+04, p50 4.97e+04, p75 201,800, p95 2,285,500, p100 421,060,000
- `protection_level_m` (float64): null 0; p0 0.3001, p5 1.093, p25 3.186, p50 5.904, p75 10.16, p95 18.44, p100 31.99
- `margin_median_m` (float64): null 0; p0 0.3, p5 0.3347, p25 0.4604, p50 0.6452, p75 1.031, p95 2.418, p100 12.95
- `margin_max_m` (float64): null 0; p0 0.3001, p5 0.3598, p25 0.6125, p50 0.9895, p75 1.692, p95 3.602, p100 12.95
- `river_ref_median_m` (float64): null 0; p0 -2.663, p5 11.05, p25 70.25, p50 117.2, p75 158.2, p95 239.4, p100 583
- `entry_x` (float64): null 0; p0 156,965, p5 222,505, p25 349,805, p50 485,905, p75 633,435, p95 710,655, p100 855,445
- `entry_y` (float64): null 0; p0 158,065, p5 241,325, p25 298,355, p50 384,135, p75 508,105, p95 676,735, p100 754,405
- `source_x` (float64): null 0; p0 160,625, p5 219,975, p25 347,845, p50 485,065, p75 635,925, p95 706,175, p100 861,895
- `source_y` (float64): null 0; p0 158,815, p5 242,005, p25 300,325, p50 388,225, p75 510,165, p95 675,585, p100 752,695
- `truncated` (bool): true 115
- `touches_corridor_edge` (bool): true 25
- `touches_dtm_edge` (bool): true 9
- `no_levee` (bool): true 0
- `unit_id` (text): 528 distinct; '013_023' 699, '027_027' 681, '020_028' 677, '020_029' 613, '029_025' 542, '030_025' 454, '025_028' 454, '007_014' 366
- `unit_label` (int64): null 0; p0 1, p5 1, p25 4, p50 10, p75 22, p95 50, p100 113
- `group_id` (int64): null 0; p0 0, p5 0, p25 119, p50 818, p75 1859, p95 3064, p100 3454
- largest 10:

| protected_area_id | group_id | area_m2 | protection_level_m | unit_id |
|---|---|---|---|---|
| POL_484160_719700_g2563 | 2563 | 421,060,000 | 5.791 | 019_005 |
| POL_484160_719700_g1061 | 1061 | 421,060,000 | 6.506 | 019_005 |
| POL_484160_719700_g1917 | 1917 | 421,059,900 | 6.45 | 019_005 |
| POL_173380_559700_g27 | 27 | 222,170,600 | 0.6135 | 003_013 |
| POL_173380_559700_g26 | 26 | 222,170,600 | 3.881 | 003_013 |
| POL_173380_559700_g2976 | 2976 | 222,170,600 | 0.7797 | 003_013 |
| POL_479170_461620_g227 | 227 | 222,150,200 | 2.207 | 018_018 |
| POL_478130_617520_g3157 | 3157 | 196,321,400 | 2.162 | 018_010 |
| POL_193840_528700_g51 | 51 | 158,761,500 | 2.48 | 004_015 |
| POL_456400_579320_g2910 | 2910 | 118,401,000 | 1.234 | 017_012 |


### leaks.gpkg

- features: 89,733; geometry types: {'Point': 89733}; empty or invalid: 0 / 0
- bounds: 152,981, 127,707, 873,442, 774,673
- `leak_id` (text): 89,733 distinct; 'POL_L1' 1, 'POL_L2' 1, 'POL_L3' 1, 'POL_L4' 1, 'POL_L5' 1, 'POL_L6' 1, 'POL_L7' 1, 'POL_L8' 1
- `cause` (text): 5 distinct; 'terrain_above_crest' 56251, 'low_crest' 27938, 'line_end' 4596, 'river_crossing' 945, 'water_source_behind' 3
- `x` (float64): null 0; p0 152,981, p5 222,451, p25 331,135, p50 451,296, p75 566,885, p95 744,947, p100 873,442
- `y` (float64): null 0; p0 127,707, p5 240,335, p25 305,768, p50 371,926, p75 547,444, p95 716,876, p100 774,673
- `levee_length_m` (float64): null 0; p0 0.007137, p5 8.156, p25 40.38, p50 52.93, p75 111.8, p95 407.9, p100 3.88e+04
- `n_samples` (int64): null 0; p0 1, p5 1, p25 1, p50 1, p75 3, p95 9, p100 780
- `line_ids` (text): 74,410 distinct; '59017' 75, '59030' 51, '66152' 44, '28526' 42, '85739' 40, '20905' 36, '4373' 33, '85820' 27
- `crest_hand_median_m` (float64): null 0; p0 -15.73, p5 1.109, p25 5.792, p50 17.29, p75 39.62, p95 116.2, p100 1339
- `level_behind_median_m` (float64): null 0; p0 -6.35, p5 0.8746, p25 5.805, p50 17.58, p75 39.65, p95 117.3, p100 1345
- `lost_margin_median_m` (float64): null 0; p0 -370.9, p5 -4.846, p25 -0.987, p50 -0.03264, p75 0.7181, p95 4.372, p100 113.1
- `unit_id` (text): 824 distinct; '020_028' 2371, '020_027' 2000, '012_025' 1436, '012_026' 1426, '026_017' 1333, '021_028' 1254, '013_023' 1245, '011_024' 1187
- largest 10:

| leak_id | cause | x | y | levee_length_m | lost_margin_median_m | line_ids |
|---|---|---|---|---|---|---|
| POL_L1 | river_crossing | 687,065 | 311,195 | 3.88e+04 | 3.09 | 24199,24200,24202,25191,25194,25202,25203,25204,25206,25207,25210,25211,25212,25214,25215,25218,25219,25220,25223,25224,25225,26017,26020,26021,26022,26097,26468,106341,106356,106365,106662,140204 |
| POL_L2 | low_crest | 303,715 | 423,235 | 2.483e+04 | 2.581 | 61449,66282,66283,88745,88746,88747,88748,88749,88751,88752,88755,88756,88757,88758,88759,88760,88761,88762,88763,88764,88765,88767,88768,88769,88770,88771,88772,88776,88962,88963,88964,88965,88968,88972,88973,88974,88977,88978,88979,89018,89023,89024,89028,89029,89030,89034,89035,89036,89037,89038,89039,89040,89041,89042,89043,89044,89045,89046,89104,89106,89116,89118,89119,89120,135561,135569,135570,135691,135714,135733,135739,135746,135781,152831,152834,152835,152837,152838,152839,152899,152904,152907,171413,186889,186933,186935,186936 |
| POL_L3 | low_crest | 271,745 | 468,485 | 2.36e+04 | 2.569 | 59221,64794,64795,64796,64797,64798,64799,64800,64801,64803,64804,64805,64814,64815,64817,64818,64825,64826,64827,64828,64829,64830,64831,64832,64894,64900,64901,64904,64905,64919,64925,64926,64931,64932,124276,124293,124294,124300,124337,147494,147502,147511,147515,147516,147541,166302 |
| POL_L4 | river_crossing | 404,085 | 331,255 | 2.177e+04 | 2.217 | 11764,53832,53833,53834,53835,53837,53838,53839,53840,53841,53842,53843,54389,54390,54391,54392,54393,54394,54402,54422,54423,82948,145606 |
| POL_L5 | river_crossing | 210,145 | 646,285 | 2.112e+04 | 0.786 | 4899,4903,4908,4910,4911,4912,4923,60097,60098,60099,60100,60168,60169,60170,60171,60185,60186,60187,92410,122181 |
| POL_L6 | low_crest | 324,725 | 402,355 | 1.984e+04 | 2.238 | 70542,70581,70582,70663,70664,70665,70668,70669,70670,70671,70672,70674,70679,70680,70699,70700,70701,70702,70703,126869,126887,148803,148804,148809,148810,148813 |
| POL_L7 | river_crossing | 699,915 | 320,455 | 1.98e+04 | 2.908 | 24431,24434,24435,24436,24439,24441,24442,24535,24536,24539,24540,24544,24545,24546,24547,24548,25230,25231,25232,25233,25234,25237,25238,25241,25247,25248,25249,25250,25317,25321,25322,25323,105899,106374,106376,106379,139796,139839,140079 |
| POL_L8 | line_end | 194,002 | 602,408 | 1.869e+04 | 0.7487 | 4166,4167,4170,4171,4173,4174,55680,55681,55682,55685,55687,55691,55701,55714,97371,97373,119998,120001 |
| POL_L9 | river_crossing | 599,845 | 245,885 | 1.711e+04 | 3.967 | 25836,28816,28831,28832,28833,28834,28835,28836,28837,28838,28839,28840,28841,28842,28843,28844,28845,28846,28847,28848,28849,28850,28868,28873,28874,28880,108130,108132,108133,108140,108148,108159,140754 |
| POL_L10 | low_crest | 698,825 | 409,985 | 1.59e+04 | 2.976 | 26503,26504,26505,41995,42001,42002 |


### barriers.gpkg

- features: 186,937; geometry types: {'LineString': 186937}; empty or invalid: 0 / 0
- bounds: 8.9e+04, 8.874e+04, 963,581, 790,439
- line length km: total 2.886e+04; per feature p0 1.025e-05, p5 0.002602, p25 0.01226, p50 0.03925, p75 0.1309, p95 0.5934, p100 34.85
- `line_id` (int64): null 0; p0 0, p5 9347, p25 4.673e+04, p50 9.347e+04, p75 140,202, p95 177,589, p100 186,936
- `source` (text): 1 distinct; 'europe-latest_man_made_embankm' 186937
- `ridge_share` (float64): null 97,805; p0 0, p5 0, p25 0, p50 0.2364, p75 0.7895, p95 1, p100 1
- `crest_height_m` (float64): null 129,169; p0 0.5, p5 0.5721, p25 0.825, p50 1.322, p75 2.415, p95 4.92, p100 14.17
- `shift_median_m` (float64): null 97,805; p0 0, p5 0, p25 0, p50 0, p75 10, p95 10, p100 10
- `no_ridge` (bool): true 55,030
- `bridged` (bool): true 97,805
- `bridge_type` (text): 3 distinct; '' 89132, 'short' 63808, 'dtm' 33997
- `crest_min_m` (float64): null 13,986; p0 -2.81, p5 29.85, p25 94.22, p50 157.8, p75 256, p95 516, p100 2180
- `crest_hand_min_m` (float64): null 20,053; p0 -16.75, p5 0.8339, p25 5.283, p50 15.04, p75 35.7, p95 110.4, p100 1333

### levees_crest.gpkg

- features: 89,132; geometry types: {'LineString': 89132}; empty or invalid: 0 / 0
- bounds: 8.9e+04, 8.874e+04, 963,581, 790,439
- line length km: total 2.616e+04; per feature p0 0.0001833, p5 0.01635, p25 0.05371, p50 0.1246, p75 0.2836, p95 1.02, p100 34.85
- `source` (text): 1 distinct; 'europe-latest_man_made_embankm' 89132
- `ridge_share` (float64): null 0; p0 0, p5 0, p25 0, p50 0.2364, p75 0.7895, p95 1, p100 1
- `crest_height_m` (float64): null 31,364; p0 0.5, p5 0.5721, p25 0.825, p50 1.322, p75 2.415, p95 4.92, p100 14.17
- `shift_median_m` (float64): null 0; p0 0, p5 0, p25 0, p50 0, p75 10, p95 10, p100 10
- `no_ridge` (bool): true 55,030

### enclosed_streams.gpkg

- features: 590; geometry types: {'Polygon': 590}; empty or invalid: 0 / 0
- bounds: 153,320, 136,160, 816,100, 744,500
- polygon area km2: total 62.14; per feature p0 0.0001, p5 0.0001, p25 0.0001, p50 0.0004, p75 0.0032, p95 0.112, p100 28.33
- `zone` (int64): null 0; p0 1, p5 1, p25 2, p50 6, p75 21, p95 54, p100 72
- `area_m2` (float64): null 0; p0 100, p5 100, p25 100, p50 400, p75 3200, p95 112,000, p100 28,329,800
- `stream_cells` (int64): null 0; p0 1, p5 1, p25 1, p50 2, p75 5, p95 28, p100 154
- `gated` (bool): true 0
- `unit_id` (text): 128 distinct; '020_028' 55, '012_024' 50, '032_022' 49, '027_028' 27, '020_029' 26, '027_016' 22, '013_025' 15, '019_028' 15

### rivers_channel.gpkg

- features: 13,799; geometry types: {'MultiLineString': 13799}; empty or invalid: 0 / 0
- bounds: 4.939e+04, 104.4, 997,958, 888,675
- line length km: total 109,325; per feature p0 0.0001341, p5 0.8552, p25 3.029, p50 6.065, p75 10.7, p95 21.26, p100 76.54
- `COMID` (text): 13,799 distinct; '22022599' 1, '22021263' 1, '22020231' 1, '22020249' 1, '22022464' 1, '22023000' 1, '22022869' 1, '22023060' 1
- `lengthkm` (float64): null 0; p0 0.05321, p5 0.891, p25 3.211, p50 6.467, p75 11.48, p95 22.74, p100 82.74
- `lengthdir` (float64): null 0; p0 0.05312, p5 0.7376, p25 2.427, p50 4.727, p75 8.009, p95 14.88, p100 47.32
- `sinuosity` (float64): null 0; p0 1, p5 1.095, p25 1.215, p50 1.33, p75 1.491, p95 1.887, p100 4.76
- `slope` (float64): null 0; p0 0, p5 4.698e-05, p25 0.0004972, p50 0.001228, p75 0.00297, p95 0.01101, p100 0.06196
- `uparea` (float64): null 0; p0 25.06, p5 27.92, p25 45, p50 106.1, p75 469, p95 9306, p100 194,369
- `order` (int64): null 0; p0 1, p5 1, p25 1, p50 1, p75 3, p95 5, p100 7
- `strmDrop_t` (float64): null 0; p0 0, p5 0, p25 1.8, p50 7, p75 19.7, p95 82.12, p100 599.8
- `slope_taud` (float64): null 0; p0 0, p5 0, p25 0.0003271, p50 0.00109, p75 0.002803, p95 0.01081, p100 0.05213
- `NextDownID` (int64): null 0; p0 0, p5 22,014,000, p25 23,010,612, p50 24,030,961, p75 24,034,282, p95 24,036,745, p100 24,038,071
- `maxup` (int64): null 0; p0 0, p5 0, p25 0, p50 0, p75 2, p95 2, p100 3
- `up1` (int64): null 0; p0 0, p5 0, p25 0, p50 0, p75 24,030,972, p95 24,036,584, p100 24,038,422
- `up2` (float64): null 0; p0 0, p5 0, p25 0, p50 0, p75 24,032,116, p95 24,037,276, p100 24,038,427
- `up3` (int64): null 0; p0 0, p5 0, p25 0, p50 0, p75 0, p95 0, p100 24,034,910
- `up4` (int64): null 0; p0 0, p5 0, p25 0, p50 0, p75 0, p95 0, p100 0
- `IFLakeType` (int64): null 0; p0 0, p5 0, p25 0, p50 0, p75 0, p95 2, p100 2
- `IFLakeLenRatio` (float64): null 0; p0 0, p5 0, p25 0, p50 0, p75 0, p95 40.65, p100 100
- `Hylak_from` (text): 824 distinct; '' 12762, '13580' 11, '1186' 8, '13403' 7, '13579' 6, '13902' 5, '164251' 5, '13456' 5
- `Hylak_list` (text): 1,288 distinct; '' 12502, '1335398' 2, '1338995' 2, '166229' 2, '1334559' 2, '164980' 2, '1320701' 2, '1338864' 2
- `Hylak_to` (text): 447 distinct; '' 12651, '13580' 22, '1186' 14, '164251' 10, '13579' 10, '13636' 10, '13799' 8, '13472' 8
- `IFAddArea` (float64): null 13,790; p0 -1450, p5 -1365, p25 -442, p50 58, p75 442, p95 1365, p100 1450
- `pfaf` (int64): null 0; p0 22, p5 22, p25 23, p50 24, p75 24, p95 24, p100 24
- `DF_Region` (text): 3 distinct; 'Baltic' 1542, 'Nordics' 36, 'Denmark' 7
- `uparea_og` (float64): null 13,799; no finite values
- `IFAddArea_where` (text): 1 distinct; 't' 4
- `from_point` (text): 13,799 distinct; '22_22640' 1, '22_21285' 1, '22_20297' 1, '22_20315' 1, '22_22504' 1, '22_23045' 1, '22_22914' 1, '22_23106' 1
- `to_point` (text): 6,979 distinct; '24_18205' 3, '24_12291' 3, '24_15320' 3, '24_8873' 3, '24_8795' 3, '24_14911' 3, '23_6368' 3, '24_11799' 3
- `from_area` (float64): null 0; p0 0, p5 0, p25 0, p50 0, p75 433.4, p95 9280, p100 193,921
- `to_area` (float64): null 0; p0 13.74, p5 27.92, p25 45, p50 106.2, p75 468.9, p95 9306, p100 193,927
- `from_area_og` (text): 228 distinct; '0.0' 4, '76.50450911879119' 1, '29111.44398699139' 1, '28979.229134280125' 1, '28869.94596940129' 1, '28765.227717099682' 1, '28536.24905501525' 1, '28334.09241398078' 1
- `to_area_og` (text): 231 distinct; '222.28775081558064' 1, '82.18219551734546' 1, '29160.794257025987' 1, '28980.79788762067' 1, '28892.664531716142' 1, '28782.865599917222' 1, '28723.04052423964' 1, '28343.34773007572' 1
- `from_area_diff` (text): 14 distinct; '0' 13568, '-442.0' 111, '0.0' 68, '442.0' 32, '-152.0' 3, '-94.0' 3, '-1450.1480353929' 3, '1237.9999999999998' 3
- `to_area_diff` (text): 18 distinct; '0' 13568, '-442.0' 107, '0.0' 64, '442.0' 29, '-442.00000000001455' 5, '-1450.1480353929' 4, '1238.0000000000002' 4, '-152.0' 3
- `source` (text): 1 distinct; 'MERIT_RN_with_IFLakeType_26072' 13799

### river_groups.gpkg

- features: 13,799; geometry types: {'MultiLineString': 13799}; empty or invalid: 0 / 0
- bounds: 4.939e+04, 104.4, 997,958, 888,675
- line length km: total 109,325; per feature p0 0.0001341, p5 0.8552, p25 3.029, p50 6.065, p75 10.7, p95 21.26, p100 76.54
- `group_id` (int64): null 0; p0 0, p5 0, p25 0, p50 60, p75 1165, p95 2774, p100 3459
- `uparea` (float64): null 0; p0 25.06, p5 27.92, p25 45, p50 106.1, p75 469, p95 9306, p100 194,369
- `length_m` (float64): null 0; p0 0.1341, p5 855.2, p25 3029, p50 6065, p75 1.07e+04, p95 2.126e+04, p100 7.654e+04

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
{"units_measured": 971, "seconds_per_domain_cell": 1.0670808918468384e-05, "bytes_per_window_cell": 199.71938461538463, "country_seconds": 4062.410645005446, "remaining_seconds": 0.0, "time_shares": {"read_terrain": 0.028336382932275903, "domain": 0.005818948844709968, "crest_read": 0.026827167201999715, "barriers": 0.003160584784000367, "rivers": 0.0033809755063290486, "undefended": 0.009063847270720769, "hand": 0.04565999822147744, "flood": 0.06868402177422996, "protected_areas": 0.004017658106068129, "leaks": 0.0032574508675142708, "groups": 0.7933221482951972, "barrier_crest": 0.00046297896069842474, "write": 0.00798821379071025}}
```

### dtm_estimates.json

```
{"country": "POL", "dtm": [{"path": "D:/01_Projects/997_SmallProjects/99_DEFENSES/02_DATA/Poland/PL_DEM_10m_borders_20220913.tif", "windows": 20, "native_s_median": 3.554456299985759, "native_s_p90": 3.842615399992792, "native_mb_s_median": 18.88155032550034, "warp_s_median": 3.927105050039245, "warp_s_p90": 4.187321840017103, "warp_native_equivalent_mb_s_median": 17.093877080945916, "window_ground_km2": 1677.7216, "sequential_mb_s": 609.6333440130032, "compression_ratio": 0.5508424043655396, "compression_mb_s": 289.4448697479267, "raw_gb": 24.572774304, "window_px": 4096, "option_a_disk_gb": 18.047634772729474, "option_a_read_s": 40.30746438875212, "option_a_compress_s": 84.89621642076457, "corridor_share": 0.9698088273036742, "option_b_disk_gb": 17.502755514545782}], "dsm": [], "terrain_read_s_total": 1536.4504137191539}
```

## Logs (newest per script)

### build_units (build_units_20261111_232918.log, 9 lines)

- 1 x WARNING Placeholder value in use: country.bytes_per_cell = #
- last lines:
```
2026-11-11 23:42:43,412 INFO build_units: Units: 971 of 1554 cores intersect the corridor (core 20000 m, halo 3000 m, working cell 10.0 m).
2026-11-11 23:43:22,284 INFO build_units: Largest window 2600 x 2600 cells needs about 1.2 GB (bytes_per_cell 180); 16 workers need about 19.5 GB of RAM in total. Domain cells to process: 6091.3 million.
2026-11-11 23:43:22,337 INFO pyogrio._io: Created 971 records
2026-11-11 23:43:22,382 INFO pyogrio._io: Created 971 records
2026-11-11 23:43:24,046 INFO pyogrio._io: Created 971 records
2026-11-11 23:43:24,084 INFO build_units: Written corridor.gpkg, units_core.gpkg, units_window.gpkg and units.csv to D:\01_Projects\997_SmallProjects\99_DEFENSES\01_MD\outputs\POL.
```

### download_worldcover (download_worldcover_20261112_223235.log, 18 lines)

- 1 x WARNING Download of https://esa-worldcover.s#.eu-central-#.amazonaws.com/v#/#/map/ESA_WorldCover_#m_#_v#_N#E#_Map.tif failed (<urlopen error [WinError #] An existing connection was forcibly closed by the remote host>); retry in 
- last lines:
```
2026-11-12 22:34:36,027 INFO gis: Downloaded https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_N51E024_Map.tif (77.9 MB).
2026-11-12 22:34:38,351 INFO gis: Downloaded https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_N54E012_Map.tif (41.9 MB).
2026-11-12 22:34:39,852 INFO gis: Downloaded https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_N54E015_Map.tif (24.2 MB).
2026-11-12 22:34:41,233 INFO gis: Downloaded https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_N54E018_Map.tif (21.5 MB).
2026-11-12 22:34:47,431 INFO gis: Downloaded https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_N54E021_Map.tif (80.2 MB).
2026-11-12 22:34:47,520 INFO download_worldcover: worldcover.vrt over 14 tiles; licence: CC-BY 4.0 - https://creativecommons.org/licenses/by/4.0/.
```

### prepare_dtm (prepare_dtm_20261111_223947.log, 11 lines)

- 1 x WARNING PL_DEM_#m_borders_#.tif: strip layout; window reads at working resolution will be slow (see the benchmark and dtm_conversion).
- last lines:
```
2026-11-11 22:42:18,368 INFO prepare_dtm: Benchmark DTM PL_DEM_10m_borders_20220913.tif (20 windows of 4096 x 4096 px = 1677.72 km2): native read 3.55 s/window (p90 3.84), 19 MB/s; warp to 10.0 m 3.93 s/window (p90 4.19), 17 MB/s native equivalent; sequential read 610 MB/s; ZSTD compression ratio 0.
2026-11-11 22:42:18,371 INFO prepare_dtm: Note: the operating system file cache can make repeated benchmarks faster than a cold read.
2026-11-11 22:42:18,375 INFO prepare_dtm: Terrain read estimate: 971 units, 2 s per unit on average, 0:01:36 in total on 16 workers (source: median warp s/window x unit window area / benchmark window area; assumes reads scale with workers, which the disk may not allow).
2026-11-11 22:42:26,138 INFO prepare_dtm: Option A for PL_DEM_10m_borders_20220913.tif: about 18.0 GB on disk (raw 24.6 GB x measured ratio 0.55 + 1/3 overviews); one sequential pass 0:00:40 (read at 610 MB/s) to 0:02:05 (plus compression at 289 MB/s); overview building not included.
2026-11-11 22:42:26,192 INFO prepare_dtm: Option B for PL_DEM_10m_borders_20220913.tif: corridor plus halo covers 97.0 % of the valid data; about 17.5 GB on disk; the source is still read once in full (0:00:40), compression of the corridor part only (0:01:22).
2026-11-11 22:42:26,208 INFO prepare_dtm: Processing DTM: D:\01_Projects\997_SmallProjects\99_DEFENSES\01_MD\work\POL\dtm\dtm.vrt -> original files
```

### protected_areas (protected_areas_20261113_095554.log, 21 lines)

- 1 x WARNING Placeholder value in use: country.bytes_per_cell = #
- 1 x WARNING # levee-enclosed areas contain river cells (see enclosed_streams.gpkg).
- last lines:
```
2026-11-13 10:07:27,262 INFO pyogrio._io: Created 89,733 records
2026-11-13 10:07:27,509 INFO protected_areas: Levees that protect nothing (leaks.gpkg), by cause: terrain_above_crest 56251 places, 5629.1 km of levee, low_crest 27938 places, 4715.3 km of levee, line_end 4596 places, 686.8 km of levee, river_crossing 945 places, 516.4 km of levee, water_source_behi
2026-11-13 10:07:32,337 INFO pyogrio._io: Created 590 records
2026-11-13 10:07:32,394 WARNING protected_areas: 590 levee-enclosed areas contain river cells (see enclosed_streams.gpkg).
2026-11-13 10:07:32,600 INFO protected_areas: Measured on 971 units: 10.67 microseconds per domain cell, up to 200 bytes per window cell (configured 180), peak RSS up to 1.26 GB. Time shares: groups 79 %, flood 7 %, hand 5 %, read_terrain 3 %, crest_read 3 %, undefended 1 %, write 1 %, domain 1 %, p
2026-11-13 10:07:32,601 INFO protected_areas: Extrapolation (measured s/cell x estimated domain cells / 16 workers): whole country 1:07:42, remaining 0 units 0:00:00; largest unit window needs about 1.4 GB.
```

Map: overview.png (241 KB)

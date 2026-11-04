# Author: Jakub Zapletal
# Date: 2026-11-04
# Writes a large synthetic DTM (two strip-layout GeoTIFFs without overviews, like delivered national mosaics),
# rivers, levees and a config for performance runs of prepare_dtm.py, build_units.py and compartments.py.
# Run: pixi run python -m tests.make_benchmark_world  (DEFENSES_BENCH_DIR sets the folder)

import json
import os
from pathlib import Path

import geopandas as gpd
import numpy as np
from osgeo import gdal, osr
from shapely.geometry import LineString

gdal.UseExceptions()

EPSG = 32633  # synthetic test data only
SIZE_M, NATIVE, ROWS_PER_PASS = 24000, 1.0, 1000
X0, Y_TOP = 400000.0, 5524000.0

root = Path(os.environ.get("DEFENSES_BENCH_DIR", "benchmark_world"))
data = root / "data"
data.mkdir(parents=True, exist_ok=True)
n = int(SIZE_M / NATIVE)
rng = np.random.default_rng(7)

# Meandering main river along x and two tributaries; valley floor rising away from the river.
xs = np.linspace(0, SIZE_M, 400)
main = np.column_stack([xs, SIZE_M / 2 + 1500 * np.sin(xs / 3000)])
trib_1 = np.column_stack([np.full(50, 7000.0), np.linspace(SIZE_M, main[117, 1], 50)])
trib_2 = np.column_stack([np.full(50, 16000.0), np.linspace(0, main[267, 1], 50)])
rivers = [LineString(main), LineString(trib_1), LineString(trib_2)]

# Levees: both banks of the main river 120 m from the centre line with gaps, plus closed rings on the floodplain.
levees = []
for side in (-1, 1):
    bank = main + np.array([0, side * 120.0])
    for start in range(0, len(bank) - 30, 40):
        levees.append(LineString(bank[start : start + 36]))
for cx in range(2000, SIZE_M - 2000, 3000):
    cy = SIZE_M / 2 + 1500 * np.sin(cx / 3000) + 600
    levees.append(
        LineString(
            [(cx, cy), (cx + 500, cy), (cx + 500, cy + 400), (cx, cy + 400), (cx, cy)]
        )
    )

srs = osr.SpatialReference()
srs.ImportFromEPSG(EPSG)
files = []
centre = np.interp(np.arange(n) * NATIVE, main[:, 0], main[:, 1])
for k, (c0, c1) in enumerate([(0, n // 2), (n // 2, n)]):
    path = data / f"dtm_block_{k}.tif"
    ds = gdal.GetDriverByName("GTiff").Create(
        str(path), c1 - c0, n, 1, gdal.GDT_Float32, ["COMPRESS=DEFLATE", "PREDICTOR=3"]
    )  # strip layout, no overviews
    ds.SetGeoTransform([X0 + c0 * NATIVE, NATIVE, 0, Y_TOP, 0, -NATIVE])
    ds.SetProjection(srs.ExportToWkt())
    ds.GetRasterBand(1).SetNoDataValue(-9999.0)
    x = (np.arange(c0, c1) + 0.5) * NATIVE
    for r in range(0, n, ROWS_PER_PASS):
        y = SIZE_M - (np.arange(r, min(r + ROWS_PER_PASS, n)) + 0.5) * NATIVE
        distance = np.abs(y[:, None] - centre[None, c0:c1])
        z = 150 + 0.0008 * x[None, :] + 0.004 * distance - 3.0 * (distance < 20)
        z += 0.5 * np.sin(x[None, :] / 230.0) * np.cos(y[:, None] / 310.0)
        z += rng.normal(0, 0.03, z.shape)  # lidar-like noise
        ds.GetRasterBand(1).WriteArray(z.astype(np.float32), 0, r)
    ds = None
    files.append(str(path))
    print(f"Written {path}")

# Embankments (2.5 m) on the rings only; the bank lines stay without embankment, like false detections.
for k, (c0, c1) in enumerate([(0, n // 2), (n // 2, n)]):
    ds = gdal.Open(files[k], gdal.GA_Update)
    band = ds.GetRasterBand(1)
    for ring in levees[-7:]:
        xmin, ymin, xmax, ymax = ring.bounds
        w0, w1 = max(int(xmin - 5) - c0, 0), min(int(xmax + 5) - c0, c1 - c0)
        r0, r1 = max(int(SIZE_M - ymax - 5), 0), min(int(SIZE_M - ymin + 5), n)
        if w1 <= w0:
            continue
        block = band.ReadAsArray(w0, r0, w1 - w0, r1 - r0)
        gx, gy = np.meshgrid(
            (np.arange(w0, w1) + c0 + 0.5) * NATIVE,
            SIZE_M - (np.arange(r0, r1) + 0.5) * NATIVE,
        )
        near = __import__("shapely").dwithin(
            ring, __import__("shapely").points(gx.ravel(), gy.ravel()), 2.0
        )
        block[near.reshape(block.shape)] += 2.5
        band.WriteArray(block, w0, r0)
    ds = None

world = lambda g: __import__("shapely").affinity.translate(g, X0, Y_TOP - SIZE_M)
gpd.GeoDataFrame(geometry=[world(g) for g in rivers], crs=EPSG).to_file(
    data / "rivers.gpkg", driver="GPKG"
)
gpd.GeoDataFrame(geometry=[world(g) for g in levees], crs=EPSG).to_file(
    data / "levees.gpkg", driver="GPKG"
)
(root / "benchmark.toml").write_text(f"""placeholders = ["country.bytes_per_cell"]
[general]
work_dir = "{(root / 'work').resolve()}"
output_dir = "{(root / 'output').resolve()}"
n_workers = 2
max_ram_gb_per_worker = 6.0
gdal_cachemax_mb = 512
warp_memory_limit_mb = 512
gdal_threads_per_worker = 2
log_level = "INFO"

[benchmark]
n_windows = 20
window_px = 4096
sequential_rows = 512
seed = 42

[countries.BM]
dtm = {json.dumps([str(Path(f).resolve()) for f in files])}
levees = "{(data / 'levees.gpkg').resolve()}"
rivers = "{(data / 'rivers.gpkg').resolve()}"
working_resolution_m = 2.0
corridor_buffer_m = 2000.0
unit_core_size_m = 7000.0
unit_halo_m = 2000.0
convert_dtm = false
dtm_conversion = "A"
dtm_compress = "ZSTD"
footprint_resolution_m = 100.0
bytes_per_cell = 180.0
crest_box_size_m = 500.0
merge_overlap_ratio = 0.5
drop_contradicted_truncated = true
contradiction_ratio = 0.5
id_rounding_m = 10.0
keep_intermediate = false

[countries.BM.compartments]
gap_tolerance_m = 20.0
touch_tolerance_m = 0.01
self_close_min_length_factor = 3.0
crest_search_radius_m = 5.0
barrier_removal_width_m = 10.0
min_margin_m = 0.3
max_level_m = 20.0
min_area_m2 = 10000.0
gated_crossings = false
fill_nodata_holes = true
stream_burn_m = 5.0
stream_burn_buffer_m = 100.0
preserve_river_cells = true
cap_river_holes = true
max_fill_distance_m = 500.0

[countries.BM.rp]
levee_adjacent_distance_m = 10.0
""")
print(f"Config: {root / 'benchmark.toml'}")
